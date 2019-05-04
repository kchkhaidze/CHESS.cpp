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

#ifndef UNIVERSE_H
#define UNIVERSE_H

// Forward declerations: ///////////////////////////////////////////////////////
//class Phylogeny_Node;
class Cell;
class CellType;
class PhylogenyRoot;
class Shape;

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>
#include "boost/multi_array.hpp"
#include "SimulationParameterSet.hpp"
#include "CellSample.hpp"

#ifdef _DISPLAY_
  // Glut:
  #include <GL/glut.h>
#endif // _DISPLAY_

// Universe ////////////////////////////////////////////////////////////////////

class Universe {
    // Store parameters:
    SimulationParameterSet mParameters;

    // Other parameters:
    long double mTime;
    std::vector< boost::multi_array<Cell*, 3> > mvSpace; // the space
    int mMaxSpaceDim;
    std::vector<CellType*> mpTypes;
    std::vector<PhylogenyRoot*> mpPhylogenies;
    bool mLimitReached;

    // Universe history:
    std::vector<long double> mvHistoryTime;
    std::vector< std::vector <unsigned long> > mvCellCounts;

    // Display related parameters:
    #ifdef _DISPLAY_
      double mLastPrintTime;
    #endif // _DISPLAY_

  public:
    // Constructor:
    Universe(SimulationParameterSet);

    // Destructor:
    //:ToDo: Add destructor
    ~Universe();

    // Getter functions:
    long double Time() const;
    bool LimitIsReached() const;
    bool ContainsCell(int, int, int, int) const;
    bool ContainsCoordinate(int, int, int, int) const;
    void Size(int, int&, int&, int&) const;
    class Cell* GetCell (int i, int x, int y, int z);
    class Cell* RemoveCell (Cell*);
    class Cell* RemoveCell (int, int, int, int);
    class CellType* NextReactionType(long double*, int*);
    std::vector< std::array<int, 3> > FreeNeighbours(int, int, int, int) const;
    std::vector <unsigned long> CellCountsPerType() const;

    // Sampling functions (getters):
    class CellSample TakeSample(int, Shape*);
    class SingleCellSample TakeSingleCellSample(int, Shape*);
    std::vector <CellSample> TakeMultipleSamples(int);
    std::vector <std::string> SingleCellTrees(int, Shape*);
    // void TakeSquares(double, int, int, int, int, int, int,
    //                  bool, std::string, std::string);

    // Setter functions:
    void IncrementTimeBy(long double);
    void MarkLimitIsReached();
    bool InsertCell(int, int, int, int, Cell*);
    bool InsertCell(int, int, int, int, Cell*, bool);
    void RegisterType(CellType *);
    void RecordToHistory();

    // Output Functions:
    void Print() const;
    void SpacesToCsvFiles(std::string) const;
    void SpaceToCsvFile(std::string) const;
    void SpaceToCsvFile(int, std::string) const;
    void SpaceGenerationToFile(int, std::string) const;
    void SampGridToCsvFile(std::string, boost::multi_array<int, 3>, int) const;
    void TypesToCsvFile(std::string) const;
    void PhylogeniesToFile(std::string) const;
    bool HistoryToFile(std::string) const;

    // Other functions:
    bool PushCell(int, double, double, double, double, double, double, unsigned int);
    bool RunSimulation(bool);
    bool RunSimulation();
    void SaveResults() const;

    // Optional features:
    #ifdef _DISPLAY_
        void Display() const;
        void Display(int) const;
        void WaitTillDisplayClosed() const;
        void SaveDisplayToFile(std::string) const;
    #endif // _DISPLAY_
};

#endif // UNIVERSE_H
