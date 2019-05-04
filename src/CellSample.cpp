/*
    A spartial constraint tumour growth model
    Copyright (C) 2018 Timon Heide (timon.heide@icr.ac.uk)
                     & Kate Chkhaidze (Kate.Chkhaidze@icr.ac.uk)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "CellSample.hpp"
#include "extern_global_variables.hpp"
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/beta_distribution.hpp>
#include <math.h>       /* pow */
#include <stdlib.h>     /* abs */
#include <iomanip>      // std::setw
#include <fstream>
#include <iostream>

#define DEPTH_OVERDISPERSION 0.08


namespace ngs_simulation {
  int simulate_depth(double depth, int depth_model) {
    // Sample depth depending upon selected depth model:
    int seq_depth;
    switch(depth_model) {

      case 0: {// Equal depth model.
        seq_depth = depth;
        break;
      }

      case 1: {// Poisson distributed depth model.
        boost::random::poisson_distribution<int> dist_depth1(depth);
        seq_depth = dist_depth1(rng);
        break;
      }

      case 2: {// overdispersed beta binomial model
        double mu = 0.6;
        double rho = DEPTH_OVERDISPERSION;
        double sh1 = mu * (1.0 / rho - 1.0); // scale param 1 (beta)
        double sh2 = (mu - 1.0) * (rho - 1.0) / rho; // scale param 2 (beta)

        // beta component
        boost::random::beta_distribution <double> beta_component_depth2(sh1, sh2);
        double cbd = beta_component_depth2(rng);

        // binomial component
        boost::random::binomial_distribution<int> dist_depth2(depth / mu, cbd);
        seq_depth = dist_depth2(rng);
        break;
      }

      case 3: {// fixed depth model
        seq_depth = depth;
        break;
      }

      default: {
        std::cerr << "Error in function AddNoiseToVaf:" << std::endl;
        std::cerr << "    Unknow depth_noise_model." << std::endl;
        exit(EXIT_FAILURE);
        break;
      }
    }
    return seq_depth;
  }


  int simulate_sequencing(int seq_depth, double seq_vaf) {
    // Do bionomial sampling at vaf:
    boost::random::binomial_distribution<int> binom(seq_depth, seq_vaf);
    int seq_reads = binom(rng);
    return seq_reads;
  }
}

// Base class TumourSample /////////////////////////////////////////////////////

CellSample::CellSample()
  : mTotalNumberCells(0)
  {}

void CellSample::AddMutation(std::string mutation_id,
                             unsigned int clone_id,
                             unsigned int number_mutated)
{
  mvMutationId.push_back(mutation_id);
  mvCloneId.push_back(clone_id);
  mvNumberMutatedCells.push_back(number_mutated);
}

void CellSample::SetSampledCellNumber(unsigned int n_cells) {
  mTotalNumberCells = n_cells;
}

SequencingResult CellSample::Sequence(const SimulationParameterSet& params) const {
  // Get options from parameter object:
  double depth = params.SequencingDepth();
  int depth_model = params.SequencingDepthModel();
  int min_reads = params.SequencingMinReads();
  double min_vaf = params.SequencingMinVaf();

  return(this->Sequence(depth, depth_model, min_reads, min_vaf));
}


SequencingResult CellSample::Sequence(double depth, int depth_model,
                                      int min_reads, double min_vaf) const
{
  SequencingResult results;

  for (std::vector<unsigned int>::size_type i = 0;
       i < mvNumberMutatedCells.size();
       i++)
  {

    // Simulate ngs:
    double exp_vaf = 0.5 * mvNumberMutatedCells[i] / mTotalNumberCells;
    int sim_depth = ngs_simulation::simulate_depth(depth, depth_model);
    int sim_reads = ngs_simulation::simulate_sequencing(sim_depth, exp_vaf);

    // Add to result object if more equal min reads:
    if (sim_reads >= min_reads && (sim_reads * 1.0 / sim_depth) > min_vaf){
      results.AddMutation(mvMutationId[i], mvCloneId[i], sim_reads, sim_depth);
    }
  }

  return results;
}



// Derived class SequencingResult //////////////////////////////////////////////


SequencingResult::SequencingResult()
  : mNumberMutations(0)
  {}

void SequencingResult::AddMutation(std::string mutation_id,
                                   unsigned int clone_id,
                                   unsigned int alt_count,
                                   unsigned int depth)
{
  mNumberMutations++;
  mvMutationId.push_back(mutation_id);
  mvCloneId.push_back(clone_id);
  mvAltCounts.push_back(alt_count);
  mvSequencedDepth.push_back(depth);
}


void SequencingResult::PrintVAFs() const {

}

unsigned int SequencingResult::MutationBurden() const {
  return mNumberMutations;
}


std::array <int, 100> SequencingResult::BinnedVaf() const {
  std::array <int, 100> binnedCounts;
  binnedCounts.fill(0);

  for (int i = 0; i < mNumberMutations; i++) {
    int bin_mutation = 100 * mvAltCounts[i] / mvSequencedDepth[i];
    binnedCounts[bin_mutation]++;
  }

  return binnedCounts;
}


double SequencingResult::Distance(const SequencingResult& x) const{
  std::array <int, 100> counts_this = this->BinnedVaf();
  std::array <int, 100> counts_x = x.BinnedVaf();
  double distance = 0;
  for (int i = 6; i < 100; i++) {
    distance += abs(counts_this[i] - counts_x[i]);
  }
  return distance;
}


void SequencingResult::Hist() const {
  // Bin the current seq results:
  std::array <int, 100> counts_this = this->BinnedVaf();

  // Determine maximum "height" of binned counts:
  int max_height = 0;
  for (int i = 6; i < 60; i ++)
    if (max_height < counts_this[i])
      max_height = counts_this[i];

  // Use this to normalize now ...
  int target_chars = 15;
  for (int i = 1; i < 60; i ++)
      counts_this[i] = counts_this[i] * (target_chars - 5) / max_height;

  // y-axis & hist body.
  std::cout << max_height << std::setw(40) << "VAF-Histogram" << std::endl;
  for (int i = target_chars; i >= 0; i--) {
    std::cout << "  | ";
    for (int j = 1; j < 60; j++)
      counts_this[j] > i ? std::cout << "#" : std::cout << " ";
    std::cout << std::endl;
  }

  // x-axis
  std::cout << " -|-";
  for (int j = 1; j < 70; j++)
    std::cout << "-";

  // x-axis label
  std::cout << std::endl;
  std::cout << "    " << "       ";
  for (int j = 10; j <= 60; j+=10)
    std::cout << std::left << std::setw(10) << j / 100.0;
  std::cout << std::endl << std::endl;
}

double* SequencingResult::Vafs() const {
  double* pVafArray = new double [mNumberMutations];

  for (int i = 0; i < mNumberMutations; i++) {
    pVafArray[i] = mvAltCounts[i] * 1.0 / mvSequencedDepth[i];
  }

  return pVafArray;
}

std::vector <double> SequencingResult::VafsVector() const {
  std::vector <double> vector;

  for (int i = 0; i < mNumberMutations; i++) {
    vector.push_back(mvAltCounts[i] * 1.0 / mvSequencedDepth[i]);
  }
  return vector;
}

std::vector <unsigned int> SequencingResult::AltVector() const {
  return mvAltCounts;
}

std::vector <unsigned int> SequencingResult::CloneIdVector() const {
  return mvCloneId;
}

std::vector <unsigned int> SequencingResult::DepthsVector() const {
  return mvSequencedDepth;
}

std::vector <std::string> SequencingResult::MutationIdVector() const {
  return mvMutationId;
}

bool SequencingResult::SaveToFile(std::string output_file) const {
  std::string delim = ",";

  // Try to open output file ...
  std::ofstream output_vaf;
  output_vaf.open(output_file, std::ios::out | std::ios::trunc);
  // ... and detect if this fails:
  if(output_vaf.is_open()){

    // Write a header line:
    output_vaf << "mutation_id" << delim
               << "clone_id" << delim
               << "alt_count" << delim
               << "depth" << delim
               << "vaf" << std::endl;

    // Write all mutation to the output file:
    for(int i = 0; i < mNumberMutations; i++) {
      output_vaf << mvMutationId[i] << delim
                 << mvCloneId[i] << delim
                 << mvAltCounts[i] << delim
                 << mvSequencedDepth[i] << delim
                 << mvAltCounts[i] * 1.0 / mvSequencedDepth[i] << std::endl;
    }
    output_vaf.close();
  } else { // Well ... opening of the output file failed.
    // Print a error message and return false. It's probably fine?
    std::cerr << "Error in function SequencingResult::SaveToFile:" << std::endl;
    std::cerr << "   Couldn't open output file." << std::endl;
    std::cerr << "   " << output_file << std::endl;
    return false;
  }
  return true;
}


// Class SingleCellSample //////////////////////////////////////////////////////

SingleCellSample::SingleCellSample() : mMaxLengthGenotype(0) {};

void SingleCellSample::appendCell(std::vector<int> geno, std::string id) {
  // Append the data to the member vectors:
  mvCellGenotypeVectors.push_back(geno);
  mvCellIds.push_back(id);

  // Set maximum genotype vector length:
  if (geno.size() > mMaxLengthGenotype) {
    mMaxLengthGenotype = geno.size();
  }
}

void SingleCellSample::print() const {

  // Print header with cell ids:
  std::cout << std::endl;
  std::cout << "## Sampled cell ids:" << std::endl;
  for (int i = 0; i < mvCellIds.size(); i++)
    std::cout << i + 1 << ": " << mvCellIds[i] << std::endl;
  std::cout << std::endl;


  // Print the result header:
  std::cout << "## Genotypes:" << std::endl;
  for (int i = 1; i <= mvCellIds.size(); i++)
    std::cout << std::setw(10) << i;
  std::cout << std::endl << std::endl;


  // Print all the genotypes:
  for (int i = 0; i < mMaxLengthGenotype; i++) { // for each cell
    for (int j = 0; j < mvCellGenotypeVectors.size(); j++) { // and node level
      if (mvCellGenotypeVectors[j].size() < i + 1) { // if passed the length of this cell genotype
        std::cout << std::setw(10) << 0; // fill with zeros
      } else {
        std::cout << std::setw(10) << mvCellGenotypeVectors[j][i];
      }
    }
    std::cout << std::endl;
  }
}
