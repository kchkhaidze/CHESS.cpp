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

#ifndef SIMULATIONARAMETERSET_H
#define SIMULATIONARAMETERSET_H

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>
#include <string>

class SimulationParameterSet {
  public:
    //
    std::vector <int> mvSizeX;
    std::vector <int> mvSizeY;
    std::vector <int> mvSizeZ;
    int mNumberOfClones;
    int mNumberOfSpaces;

    int mNumberClonalMutations;
    double mSequencingDepth;
    double mMinSequencingVaf;
    int mMinSequencingReads;
    int mDepthModel;

    std::vector <double> mMutationrates;
    std::vector <double> mBirthrates;
    std::vector <double> mDeathrates;
    std::vector <double> mAggressions;
    std::vector <unsigned int> mPushPower;
    std::vector <double> mCloneStartTimes;
    std::vector <double> mKillRegrowTime;
    std::vector <unsigned int> mFathers;
    std::vector <unsigned int> mUniverses;

    int mSeed;
    int mSimId;

    double mDisplayFrequency;

    std::string mOutputPrefix;

  public:

  // Constructor:
  SimulationParameterSet(int argc, char* argv[]);

  SimulationParameterSet(int, int, int,             // size_x, size_y, size_z
                          std::vector<double>,       // mutation_rates
                          std::vector<double>,       // birthrates
                          std::vector<double>,       // deathrates
                          std::vector<double>,       // aggressions
                          std::vector<unsigned int>, // push_power
                          std::vector<double>,       // clone_start_times
                          std::vector<double>,       // kill_regrow_times
                          std::vector<unsigned int>, // fathers
                          int,                       // seed
                          int);                      // clonal mutations


    // Getters:
    std::vector <int> SizeX() const;
    std::vector <int> SizeY() const;
    std::vector <int> SizeZ() const;
    int NumberOfClones() const;
    int NumberOfSpaces() const;
    int Seed() const;

    int NumberClonalMutations() const;
    double SequencingDepth() const;
    double SequencingMinVaf() const;
    int SequencingMinReads() const;
    int SequencingDepthModel() const;

    std::vector <double> Mutationrates() const;
    std::vector <double> Birthrates() const;
    std::vector <double> Deathrates() const;
    std::vector <double> Aggressions() const;
    std::vector <unsigned int> PushPower() const;
    std::vector <double> CloneStartTimes() const;
    std::vector <double> KillRegrowTime() const;
    std::vector <unsigned int> Fathers() const;
    std::vector <unsigned int> Spaces() const;
    double DisplayFrequency() const;
    std::string OutputPrefix() const;
    int SimulationId() const;


  // Output functions:
  void Print() const;
  bool SaveToFile(std::string) const;
};

#endif // SIMULATIONARAMETERSET_H
