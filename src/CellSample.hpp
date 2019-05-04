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

#ifndef CELL_SAMPLE_H
#define CELL_SAMPLE_H

// Forward declerations: ///////////////////////////////////////////////////////
class SequencingResult;

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>
#include <string>
#include <array>
#include "SimulationParameterSet.hpp"


// Base class CellSample /////////////////////////////////////////////////////
class CellSample {
    unsigned int mTotalNumberCells;
    std::vector <std::string> mvMutationId;
    std::vector <unsigned int> mvCloneId;
    std::vector <unsigned int> mvNumberMutatedCells;

   public:
     // Constructor
     CellSample();

     // Setter functions:
     void AddMutation(std::string, unsigned int, unsigned int);
     void SetSampledCellNumber(unsigned int);

     // Getter functions:
     class SequencingResult Sequence(const SimulationParameterSet&) const;
     class SequencingResult Sequence(double, int, int, double) const;
};

 // Derived class SequencingResult //////////////////////////////////////////////
 class SequencingResult {
    unsigned int mNumberMutations;
    std::vector <std::string> mvMutationId;
    std::vector <unsigned int> mvCloneId;
    std::vector <unsigned int> mvAltCounts;
    std::vector <unsigned int> mvSequencedDepth;

    public:
      SequencingResult();
      void AddMutation(std::string, unsigned int, unsigned int, unsigned int);
      void PrintVAFs() const;
      unsigned int MutationBurden() const;
      std::array <int, 100> BinnedVaf() const;
      double Distance(const SequencingResult&) const;
      void Hist() const;
      double* Vafs() const;
      bool SaveToFile(std::string) const;

      std::vector <double> VafsVector() const;
      std::vector <unsigned int> AltVector() const;
      std::vector <unsigned int> CloneIdVector() const;
      std::vector <unsigned int> DepthsVector() const;
      std::vector <std::string> MutationIdVector() const;

};

class SingleCellSample {
  unsigned long mMaxLengthGenotype;
  std::vector <std::string> mvCellIds;
  std::vector <std::vector<int> > mvCellGenotypeVectors;

  public:
    SingleCellSample();
    void appendCell(std::vector<int>, std::string);
    void print() const;
};

#endif // CELL_SAMPLE_H
