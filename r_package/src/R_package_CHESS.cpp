#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "boost/multi_array.hpp"
#include <iostream>
#include <cassert>
#include <time.h>
#include "Phylogeny.hpp"
#include "Universe.hpp"
#include "Cell.hpp"
#include "CellType.hpp"
#include "Shapes.hpp"
#include "CellSample.hpp"
#include "SimulationParameterSet.hpp"
#include "R_package_CHESS.hpp"
#include "extern_global_variables.hpp"
#include "Rcpp.h"
using namespace Rcpp;

// Default params:
#define SIZE_X 100 // size of the space [int]
#define SIZE_Y 100 // size of the space [int]
#define SIZE_Z 1 // size of the space [int]
#define MU 10          // mutation rate [double]
#define WTBR 1         // wildtype birth rate [double]
#define MUTBR 1.5        // mutant birth rate   [double]
#define CLST 1.0       // clone start time [double]
#define KLRT -1.0       // kill-regrow time [double]
#define DISPFRQ 0.05   // display frequency [double]
#define ALPHA 0.0      // alpha [double]
#define BETA 1.0       // beta [double]
#define AGGR 1.0       // cell aggression [double]
#define CENTER(X) ((X + 1) / 2) - 1
#define OUTPUT_DIR "."
#define SPACE_OUPUT_FILE_PREFIX "space"
#define TYPES_OUPUT_FILE_PREFIX "types"
#define PHYLO_OUPUT_FILE_PREFIX "phylo"
#define HISTORY_OUPUT_FILE_PREFIX "cell_number.csv"

boost::random::mt19937_64 rng;

Universe_class_wrap::Universe_class_wrap(int size_x, int size_y, int size_z,
                                         Rcpp::NumericMatrix clone_params,
                                         int seed,
                                         bool verbose,
                                         int clonal_mutations) {

  std::vector<double> vBirthrates, vDeathrates, vAggressions,
                      vMutation_rates, vClone_start_times, vKill_regrow_times;

  std::vector <unsigned int> vFather;

  std::vector <unsigned int> vPushpowers;

  for (int i = 0; i < clone_params.ncol(); i++) {
    vBirthrates.push_back(clone_params(0,i));
    vDeathrates.push_back(clone_params(1,i));
    vAggressions.push_back(clone_params(2,i));
    vPushpowers.push_back(clone_params(3,i));
    vMutation_rates.push_back(clone_params(4,i));
    vClone_start_times.push_back(clone_params(5,i));
    vKill_regrow_times.push_back(clone_params(6,i));
    vFather.push_back(clone_params(7,i));
  }

  SimulationParameterSet params(size_x, size_y, size_z,
                                 vMutation_rates,
                                 vBirthrates,
                                 vDeathrates,
                                 vAggressions,
                                 vPushpowers,
                                 vClone_start_times,
                                 vKill_regrow_times,
                                 vFather,
                                 seed,
                                 clonal_mutations);

  // Print parameters if called in verbose mode:
  if (verbose)
    params.Print();

  // Create the universe and run simulation:
  mpUniverse = new Universe(params);
  mpUniverse->RunSimulation();
}


Rcpp::List Universe_class_wrap::TakeSample(
  int minx, int miny, int minz, int maxx, int maxy, int maxz,
  double dp, int dp_model, int min_reads, double min_vaf,
  int seed)
{

  // Sample
  Shape* sample_shape = new Box(minx, miny, minz, maxx, maxy, maxz);
  SequencingResult result =
    mpUniverse->TakeSample(0,sample_shape).Sequence(dp, dp_model, min_reads, min_vaf);
  delete sample_shape;

  // Get std::vectors containing result elements:
  std::vector <unsigned int> alt_cpp_type = result.AltVector();
  std::vector <unsigned int> clone_id_cpp_type = result.CloneIdVector();
  std::vector <unsigned int> depth_cpp_type = result.DepthsVector();
  std::vector <std::string> mutation_id_cpp_type = result.MutationIdVector();

  // Allocate R vector types:
  Rcpp::NumericVector alt(alt_cpp_type.size());
  Rcpp::NumericVector clone_id(clone_id_cpp_type.size());
  Rcpp::NumericVector depth_res(depth_cpp_type.size());
  Rcpp::CharacterVector mutation_id(mutation_id_cpp_type.size());

  // Fill R vectors with the data:
  for (size_t i = 0; i < alt_cpp_type.size(); i++) { // vectors should all have the same length ...
    alt[i] = alt_cpp_type[i];
    depth_res[i] = depth_cpp_type[i];
    clone_id[i] = clone_id_cpp_type[i];
    mutation_id[i] = mutation_id_cpp_type[i];
  }

  std::vector <unsigned long> cellcounts;
  cellcounts = mpUniverse->CellCountsPerType();

  Rcpp::NumericVector cell_counts(cellcounts.size());
  for(int i = 0; i < cellcounts.size(); i++){
    cell_counts[i] = cellcounts[i];
  }

  // Return results as list containing a data.frame:
  return Rcpp::List::create(
    Rcpp::Named("number_mutations") =
      alt_cpp_type.size(),
    Rcpp::Named("mutation_data") =
      Rcpp::DataFrame::create(
        Rcpp::Named("clone") = clone_id,
        Rcpp::Named("alt") = alt,
        Rcpp::Named("depth") = depth_res,
        Rcpp::Named("id") = mutation_id,
        _["stringsAsFactors"] = false),
    Rcpp::Named("cell_counts") =
      Rcpp::DataFrame::create(
        Rcpp::Named("cellcounts") = cell_counts
      )
    );
}


Rcpp::List Universe_class_wrap::SingleCellTree(
  std::vector <int> x, std::vector <int> y, std::vector <int> z)
{

  // Sample
  Shape* sample_shape = new PositionVector(x, y, z);
  std::vector <std::string> result = mpUniverse->SingleCellTrees(0,sample_shape);
  delete sample_shape;

  // Allocate R vector types:
  Rcpp::CharacterVector tree_strings(result.size());

  // Fill R vectors with the data:
  for (size_t i = 0; i < result.size(); i++) {
    tree_strings[i] = result[i];
  }

  std::vector <unsigned long> cellcounts;
  cellcounts = mpUniverse->CellCountsPerType();

  Rcpp::NumericVector cell_counts(cellcounts.size());
  for(int i = 0; i < cellcounts.size(); i++){
    cell_counts[i] = cellcounts[i];
  }

  // Return results as list containing a data.frame:
  //return tree_strings;
  return Rcpp::List::create(
    Rcpp::Named("tree_string") = tree_strings,
    Rcpp::Named("cell_counts") =
      Rcpp::DataFrame::create(
        Rcpp::Named("cellcounts") = cell_counts
      )
  );
}


Rcpp::List Universe_class_wrap::SingleCellTreeN(
  int nsamples, int grid_x, int grid_y, int grid_z)
{
  std::vector <int> x;
  std::vector <int> y;
  std::vector <int> z;

  int cnt = 0;
  while (cnt < nsamples) {

    boost::random::uniform_int_distribution<int> unif_coord_x(0, grid_x);
    int coord_x = unif_coord_x(rng);
    boost::random::uniform_int_distribution<int> unif_coord_y(0, grid_y);
    int coord_y = unif_coord_y(rng);
    boost::random::uniform_int_distribution<int> unif_coord_z(0, grid_z);
    int coord_z = unif_coord_z(rng);

    if (mpUniverse->ContainsCoordinate(0, coord_x, coord_y, coord_z) &&
        mpUniverse->ContainsCell(0, coord_x, coord_y, coord_z)) {

      x.push_back(coord_x);
      y.push_back(coord_y);
      z.push_back(coord_z);
      cnt++;
    }
  }

  // Sample
  Shape* sample_shape = new PositionVector(x, y, z);
  std::vector <std::string> result = mpUniverse->SingleCellTrees(0,sample_shape);
  delete sample_shape;

  // Allocate R vector types:
  Rcpp::CharacterVector tree_strings(result.size());

  // Fill R vectors with the data:
  for (size_t i = 0; i < result.size(); i++) {
    tree_strings[i] = result[i];
  }

  std::vector <unsigned long> cellcounts;
  cellcounts = mpUniverse->CellCountsPerType();

  Rcpp::NumericVector cell_counts(cellcounts.size());
  for(int i = 0; i < cellcounts.size(); i++){
    cell_counts[i] = cellcounts[i];
  }

  // Return results as list containing a data.frame:
  //return tree_strings;
  return Rcpp::List::create(
    Rcpp::Named("tree_string") = tree_strings,
    Rcpp::Named("cell_counts") =
      Rcpp::DataFrame::create(
        Rcpp::Named("cellcounts") = cell_counts
      )
  );
}

Universe_class_wrap::~Universe_class_wrap() {
  delete mpUniverse;
}

RCPP_MODULE(CHESS_Universe_object) {
  using namespace Rcpp;

  class_<Universe_class_wrap>("CHESS_Universe_object")
  .constructor<int, int, int, Rcpp::NumericMatrix, int, bool, int>()
  .method("TakeSample", &Universe_class_wrap::TakeSample)
  .method("SingleCellTree", &Universe_class_wrap::SingleCellTree)
  .method("SingleCellTreeN", &Universe_class_wrap::SingleCellTreeN)
  ;
}
