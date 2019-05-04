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

#include "SimulationParameterSet.hpp"
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <algorithm>

#ifdef _DISPLAY_
#define DI(x) x
#else
#define DI(x)
#endif

// Params of old function simulateTumorGrowth2D
#define SIM_ID 1  // simulation run id [int]
#define SIZE_X 100 // size of the space [int]
#define SIZE_Y 100 // size of the space [int]
#define SIZE_Z 1 // size of the space [int]
#define MU 10.0          // mutation rate [double]
#define BR 1.0         // wildtype birth rate [double]
#define DR 0.0        // wildtype death rate [double]
#define CLST 1.0      // clone start time [double]
#define KLRT 0.0       // kill regrow time [double]
#define AGGR 1.0    // cell aggression [double] for wt
#define POWER SIZE_X
#define DISPFRQ 0.05   // display frequency [double]
#define CENTER(X) ((X + 1) / 2) - 1
#define OUTPUT_DIR "./"
#define SEED time(NULL)
#define MUTS_CLONAL 0    // clonal mutations [int]
#define DEPTH 100.0
#define DEPTHMODEL 1 // (1: poisson, 2: overdispersed beta-binomial)
#define MIN_READS 2
#define MIN_VAF 0.0


namespace {
  void printHelp(char* argv[]) {
    std::cout <<  "\n";
    std::cout << "A spartial constraint tumour growth model\n";
    std::cout << "Copyright (C) 2018, Timon Heide (timon.heide@icr.ac.uk)\n";
    std::cout << "                  & Kate Chkhaidze (Kate.Chkhaidze@icr.ac.uk)\n";
    std::cout << "\n";
    std::cout << "This program is free software; you can redistribute it and/or modify\n";
    std::cout << "it under the terms of the GNU General Public License as published by\n";
    std::cout << "the Free Software Foundation; either version 3 of the License, or\n";
    std::cout << "(at your option) any later version.";
    std::cout << "\n";
    std::cout << "This program is distributed in the hope that it will be useful,\n";
    std::cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
    std::cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
    std::cout << "GNU General Public License for more details.\n";
    std::cout << "\n";
    std::cout << "You should have received a copy of the GNU General Public License\n";
    std::cout << "along with this program. If not, see http://www.gnu.org/licenses/.\n";
    std::cout << "\n";
    std::cout << "\n";
    std::cout <<  "Usage: " << argv[0] << " [options]" << std::endl;
    std::cout <<  "  Options: " << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -i ID, --sim_id=ID\n";
    std::cout <<  "        Simulation Id\n";
    std::cout <<  "          default: " << SIM_ID << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -x N, --size_x=N\n";
    std::cout <<  "        Size of dimension x.\n";
    std::cout <<  "          default: " << SIZE_X << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -y N, --size_y=N\n";
    std::cout <<  "        Size of dimension y.\n";
    std::cout <<  "          default: " << SIZE_Y << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -z N, --size_z=N\n";
    std::cout <<  "        Size of dimension z.\n";
    std::cout <<  "          default: " << SIZE_Y << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -m MU, --mutation_rates=MU\n";
    std::cout <<  "        Mutation rates of cell types.\n";
    std::cout <<  "        The number of new mutations\n";
    std::cout <<  "        will drawn randomly from a Poisson \n";
    std::cout <<  "        distribution with lambda = MU during\n";
    std::cout <<  "        each division.\n";
    std::cout <<  "          default: " << MU << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -b F, --birthrates=F\n";
    std::cout <<  "        Birth rate of cell types.\n";
    std::cout <<  "          default: " << BR << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -d F, --deathrates=F\n";
    std::cout <<  "        Death rate of cell types.\n";
    std::cout <<  "          default: " << DR << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -r F, --aggressions=F\n";
    std::cout <<  "        Probabilty of cell types to push.\n";
    std::cout <<  "          default: " << AGGR << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -P N, --push_powers=N\n";
    std::cout <<  "        Push power of cell types.\n";
    std::cout <<  "          default: " << POWER << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -t F, --clone_start_times=F\n";
    std::cout <<  "        Time cell types are introduced.\n";
    std::cout <<  "          default: " << CLST << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -w F, --kill_regrow_times=F\n";
    std::cout <<  "        Time to kill 99% of cells and regrow the tumour.\n";
    std::cout <<  "          default: " << KLRT << " (never kill)" << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -o DIR, --output_prefix=DIR\n";
    std::cout <<  "        Output prefix.\n";
    std::cout <<  "          default: " << OUTPUT_DIR << std::endl;
    std::cout <<  "\n";
    DI(std::cout <<  "    -f FREQ, --display_frequency=FREQ\n";)
    DI(std::cout <<  "        Update frequency of the display.\n";)
    DI(std::cout <<  "          default: " << DISPFRQ << std::endl;)
    DI(std::cout <<  "\n";)
    std::cout <<  "    -s SEED, --seed=SEED\n";
    std::cout <<  "        Random seed.\n";
    std::cout <<  "          default: time(NULL)" << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -M N, --number_clonal=N\n";
    std::cout <<  "        Number of clonal mutations.\n";
    std::cout <<  "          default: " << MUTS_CLONAL << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -D F, --depth=F\n";
    std::cout <<  "        Average sequencing depth.\n";
    std::cout <<  "            default: " << DEPTH << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -X N, --depth_model=N\n";
    std::cout <<  "        Model for the sequencing depth.\n";
    std::cout <<  "          1) Poisson distributed depth.\n";
    std::cout <<  "          2) Beta-binomial distributed depth.\n";
    std::cout <<  "          3) Fixed depth.\n";
    std::cout <<  "            default: " << DEPTHMODEL << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -R N, --min_reads=N\n";
    std::cout <<  "        Minimum number of alternate reads.\n";
    std::cout <<  "            default: " << MIN_READS << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -V F, --min_vaf=F\n";
    std::cout <<  "        Minimum variant allele frequency.\n";
    std::cout <<  "            default: " << MIN_VAF << std::endl;
    std::cout <<  "\n";
    std::cout <<  "    -h, --help\n";
    std::cout <<  "        Print this help\n";
    std::cout << std::endl;
  }
}


// Constructor:
SimulationParameterSet::SimulationParameterSet(int argc, char* argv[])
  : mNumberClonalMutations(MUTS_CLONAL),
    mSequencingDepth(DEPTH),
    mMinSequencingVaf(MIN_VAF),
    mMinSequencingReads(MIN_READS),
    mDepthModel(DEPTHMODEL),
    mSeed(SEED),
    mSimId(0),
    mDisplayFrequency(DISPFRQ),
    mOutputPrefix(OUTPUT_DIR)
{

    mvSizeX.push_back(SIZE_X);
    mvSizeY.push_back(SIZE_Y);
    mvSizeZ.push_back(SIZE_Z);


    // Tell user about the help page:
    if (argc <= 1) {
      std::cout << "\n";
      std::cout << "A spartial constraint tumour growth model\n";
      std::cout << "Copyright (C) 2018, Timon Heide (timon.heide@icr.ac.uk)\n";
      std::cout << "Copyright (C) 2018, Kate Chkhaidze (Kate.Chkhaidze@icr.ac.uk)\n";
      std::cout << "\n";
      std::cout << "Programm will run with default parameters.\n";
      std::cout << "Please type '" << argv[0] << " -h' for help.\n";
      std::cout << std::endl;
    }

    // The possible arguments for getopts
     static struct option long_options[] = {
       {"size_x",  required_argument, nullptr, 'x'},
       {"size_y",  required_argument, nullptr, 'y'},
       {"size_z",  required_argument, nullptr, 'z'},

       {"mutation_rates", required_argument, nullptr, 'm'},
       {"birthrates",  required_argument, nullptr, 'b'},
       {"deathrates",  required_argument, nullptr, 'd'},
       {"aggressions",    required_argument, nullptr,  'r'},
       {"push_powers",    required_argument, nullptr,  'P'},
       {"clone_start_times",    required_argument, nullptr, 't'},
       {"kill_regrow_times", required_argument, nullptr, 'w'},
       {"father_types",    required_argument, nullptr, 'p'},
       {"universes",    required_argument, nullptr, 'u'},
       {"display_frequency",    required_argument, nullptr,  'f'},
       {"sim_id",  required_argument, nullptr, 'i'},
       {"output_dir",    required_argument, nullptr,  'o'},
       {"seed",    required_argument, nullptr,  's'},

       {"number_clonal",    required_argument, nullptr,  'M'},
       {"depth",    required_argument, nullptr,  'D'},
       {"depth_model",    required_argument, nullptr,  'X'},
       {"min_reads",    required_argument, nullptr,  'R'},
       {"min_vaf",    required_argument, nullptr,  'V'},

       {"help",    no_argument, nullptr,  'h'},
       {NULL, 0, nullptr, 0}
     };

     // Parse commandline arguments (from getop man-page):
     while (true) {
      int c = getopt_long(argc, argv, "x:y:z:m:b:d:r:P:t:w:p:u:f:i:o:s:M:D:X:R:V:h",
                          long_options, NULL);

      if (c == -1)
        break;

      std::stringstream ss;
      int elm_i;
      double elm_d;

      switch (c) {
        case 'x':
        ss.str(optarg);
        mvSizeX.clear();
        while (ss >> elm_i){
          if (elm_i < 1) {
            std::cout << "Error: " << std::endl;
            std::cout << "  Universe size (x) should be > 0." << std::endl;
            exit(EXIT_FAILURE);
          }
          mvSizeX.push_back(elm_i);

          if (ss.peek() == ',') {
            ss.ignore();
          }
        }
        break;

        case 'y':
        ss.str(optarg);
        mvSizeY.clear();
        while (ss >> elm_i){
          if (elm_i < 0) {
            std::cout << "Error: " << std::endl;
            std::cout << "  Universe size (y) should be >= 0." << std::endl;
            exit(EXIT_FAILURE);
          }
          mvSizeY.push_back(elm_i);

          if (ss.peek() == ',') {
            ss.ignore();
          }
        }
        break;

        case 'z':
        ss.str(optarg);
        mvSizeZ.clear();
        while (ss >> elm_i){
          if (elm_i < 1) {
            std::cout << "Error: " << std::endl;
            std::cout << "  Universe size (z) should be >= 0." << std::endl;
            exit(EXIT_FAILURE);
          }
          mvSizeZ.push_back(elm_i);

          if (ss.peek() == ',') {
            ss.ignore();
          }
        }
        break;

         case 'm':
           ss.str(optarg);
           mMutationrates.clear();
           while (ss >> elm_d){
             if (elm_d < 0.0) {
               std::cout << "Error: " << std::endl;
               std::cout << "  Mutation rate (m) should be >= 0." << std::endl;
               exit(EXIT_FAILURE);
             }
             mMutationrates.push_back(elm_d);

             if (ss.peek() == ',') {
               ss.ignore();
             }
           }
           break;
         case 'b':
           ss.str(optarg);
           mBirthrates.clear();
           while (ss >> elm_d){
             if (elm_d < 0.0) {
               std::cout << "Error: " << std::endl;
               std::cout << "  Birth rates should be >= 0." << std::endl;
               exit(EXIT_FAILURE);
             }
             mBirthrates.push_back(elm_d);

             if (ss.peek() == ',') {
               ss.ignore();
             }
           }
           break;
         case 'd':
           ss.str(optarg);
           mDeathrates.clear();
           while (ss >> elm_d){
             if (elm_d < 0.0) {
               std::cout << "Error: " << std::endl;
               std::cout << "  Death rates should be >= 0." << std::endl;
               exit(EXIT_FAILURE);
             }
             mDeathrates.push_back(elm_d);

             if (ss.peek() == ',') {
               ss.ignore();
             }
           }
           break;
         case 'r':
           ss.str(optarg);
           mAggressions.clear();
           while (ss >> elm_d){
             if (elm_d < 0.0 || elm_d > 1.0 ) {
               std::cout << "Error: " << std::endl;
               std::cout << "  Parameter aggression (r) should be 0 < R < 1.";
               std::cout << std::endl;
               exit(EXIT_FAILURE);
             }
             mAggressions.push_back(elm_d);

             if (ss.peek() == ',') {
               ss.ignore();
             }
           }
           break;
         case 'P':
         ss.str(optarg);
         mPushPower.clear();
         while (ss >> elm_i){
           if (elm_i <= 0) {
             std::cout << "Error: " << std::endl;
             std::cout << "  Push power should be >= 0." << std::endl;
             exit(EXIT_FAILURE);
           }
           mPushPower.push_back(elm_i);
           if (ss.peek() == ',') {
             ss.ignore();
           }
         }
         break;
         case 't':
           ss.str(optarg);
           mCloneStartTimes.clear();
           while (ss >> elm_d){
             if (elm_d < 0.0) {
               std::cout << "Error: " << std::endl;
               std::cout << "  Clone start time should be >= 0." << std::endl;
               exit(EXIT_FAILURE);
             }
             mCloneStartTimes.push_back(elm_d);

             if (ss.peek() == ',') {
               ss.ignore();
             }
           }
           break;
         case 'w':
           ss.str(optarg);
           mKillRegrowTime.clear();
           while (ss >> elm_d){
             if (elm_d < 0.0) {
               std::cout << "Error: " << std::endl;
               std::cout << "  Kill and regrow time should be >= 0." << std::endl;
               exit(EXIT_FAILURE);
             }
             mKillRegrowTime.push_back(elm_d);

             if (ss.peek() == ',') {
               ss.ignore();
             }
           }
           break;
         //
         case 'p':
           ss.str(optarg);
           mFathers.clear();
           while (ss >> elm_i){
             mFathers.push_back(elm_i);
             if (ss.peek() == ',') {
               ss.ignore();
             }
           }
           break;
          case 'u':
           ss.str(optarg);
           mUniverses.clear();
           while (ss >> elm_i){
             mUniverses.push_back(elm_i);
             if (ss.peek() == ',') {
               ss.ignore();
             }
           }
           break;
         case 'f':
           DI(mDisplayFrequency = atof(optarg);)
           DI(if (mDisplayFrequency < 0.0) {)
           DI(  std::cout << "Error: " << std::endl;)
           DI(  std::cout << "  Display frequency should be >= 0." << std::endl;)
           DI(  exit(EXIT_FAILURE);)
           DI(})
           break;
         case 'i':
           mSimId = strtol(optarg, NULL, 10);
           break;
         case 'o':
           mOutputPrefix = optarg;
           break;
         case 's':
           mSeed = atoi(optarg);
           break;

        case 'M':
          ss.str(optarg);
          ss >> elm_i;
          mNumberClonalMutations = elm_i;
          if (elm_i < 0) {
            std::cout << "Error: " << std::endl;
            std::cout << "  Number of clonal muts. (M) must be >=0.\n";
            exit(EXIT_FAILURE);
          }
        break;
        //
        case 'D':
          ss.str(optarg);
          ss >> elm_d;
          mSequencingDepth = elm_d;
          if (elm_d < 0) {
            std::cout << "Error: " << std::endl;
            std::cout << "  Depth (D) must be >0.\n";
            exit(EXIT_FAILURE);
          }
        break;
        //
        case 'X':
          ss.str(optarg);
          ss >> elm_i;
          mDepthModel = elm_i;
          if (elm_i < 0 || elm_i > 3) {
            std::cout << "Error: " << std::endl;
            std::cout << "  Invalid sequencing depth model." << std::endl;
            exit(EXIT_FAILURE);
          }
        break;
        //
        case 'R':
          ss.str(optarg);
          ss >> elm_i;
          mMinSequencingReads = elm_i;
          if (elm_i < 0) {
            std::cout << "Error: " << std::endl;
            std::cout << "  Min. number of sequencing reads should be >= 0.\n";
            exit(EXIT_FAILURE);
          }
        break;
        //
        case 'V':
          ss.str(optarg);
          ss >> elm_d;
          mMinSequencingVaf = elm_d;
          if (elm_d < 0.0 || elm_d > 1.0) {
            std::cout << "Error: " << std::endl;
            std::cout << "  Min. VAF should be between 0 and 1.\n";
            exit(EXIT_FAILURE);
          }
        break;
         case ':':
           break;
         case 'h':
         case '?':
           printHelp(argv);
           exit(EXIT_FAILURE);
   	      break;
         default:
           std::cout << " returned character code 0" << c << std::endl;
       }
     }

    if (mvSizeX.size() != mvSizeY.size() || mvSizeX.size() != mvSizeZ.size()) {
      std::cerr << "Length of space dimensions have to be of same length!\n";
      exit(EXIT_FAILURE);
    }
    mNumberOfSpaces = mvSizeX.size();

    // Check if all universes for cell types are defined:
    for (int i = 0; i < mUniverses.size(); i++) {
      if (mUniverses[i] >= mNumberOfSpaces) {
        std::cout << "Error: " << std::endl;
        std::cout << "  Universe vector defines a non existent universe\n";
        exit(EXIT_FAILURE);
      }
    }


    // Find the maximum clone param vector length:
    std::vector <int> param_v_sizes;
    param_v_sizes.push_back(static_cast<int>(mMutationrates.size()));
    param_v_sizes.push_back(static_cast<int>(mBirthrates.size()));
    param_v_sizes.push_back(static_cast<int>(mDeathrates.size()));
    param_v_sizes.push_back(static_cast<int>(mAggressions.size()));
    param_v_sizes.push_back(static_cast<int>(mPushPower.size()));
    param_v_sizes.push_back(static_cast<int>(mCloneStartTimes.size()));
    param_v_sizes.push_back(static_cast<int>(mKillRegrowTime.size()));
    param_v_sizes.push_back(static_cast<int>(mFathers.size()));
    param_v_sizes.push_back(static_cast<int>(mUniverses.size()));
    mNumberOfClones = *max_element(param_v_sizes.begin(), param_v_sizes.end());

    // Test if the length of all clone parameters is equal:
    for (std::vector<int>::size_type i = 0; i < param_v_sizes.size(); i++) {
      if (param_v_sizes[i] != 0 && param_v_sizes[i] != mNumberOfClones) {
        std::cerr << "Error during creation of SimulationParameterSet:\n";
        std::cerr << "  Clone parameter vectors not of same length.\n";
        exit(EXIT_FAILURE);
      }
    }

    // If there are no defined types then simulate two neutral clones:
    if (mNumberOfClones == 0) {
      mNumberOfClones = 2;
    }

    // Fill non defined param vectors with default values:
    if (mMutationrates.size() == 0)
      for (int i = 0; i < mNumberOfClones; i++)
        mMutationrates.push_back(MU);

    if (mBirthrates.size() == 0)
      for (int i = 0; i < mNumberOfClones; i++)
        mBirthrates.push_back(BR);

    if (mDeathrates.size() == 0)
      for (int i = 0; i < mNumberOfClones; i++)
        mDeathrates.push_back(DR);

    if (mAggressions.size() == 0)
      for (int i = 0; i < mNumberOfClones; i++)
        mAggressions.push_back(AGGR);

    if (mPushPower.size() == 0)
      for (int i = 0; i < mNumberOfClones; i++)
        mPushPower.push_back(POWER);

    if (mCloneStartTimes.size() == 0)
      for (int i = 0; i < mNumberOfClones; i++)
        mCloneStartTimes.push_back(CLST * i);

    if (mKillRegrowTime.size() == 0)
      for (int i = 0; i < mNumberOfClones; i++)
        mKillRegrowTime.push_back(KLRT * i);

    if (mFathers.size() == 0)
      for (int i = 0; i < mNumberOfClones; i++)
        mFathers.push_back(0);

    if (mUniverses.size() == 0)
      for (int i = 0; i < mNumberOfClones; i++)
        mUniverses.push_back(0);
}




//
SimulationParameterSet::SimulationParameterSet(
    int size_x, int size_y, int size_z,
    std::vector<double> mutation_rates,
    std::vector<double>  birthrates,
    std::vector<double>  deathrates,
    std::vector<double>  aggressions,
    std::vector<unsigned int> push_power,
    std::vector<double> clone_start_times,
    std::vector<double> kill_regrow_times,
    std::vector<unsigned int> fathers,
    int seed,
    int clonal_mutations
  )
  :
    mNumberOfSpaces(1),

    mNumberClonalMutations(clonal_mutations),
    mSequencingDepth(-1.0),
    mMinSequencingVaf(-1.0),
    mMinSequencingReads(-1),
    mDepthModel(-1),

    mMutationrates(mutation_rates),
    mBirthrates(birthrates),
    mDeathrates(deathrates),
    mAggressions(aggressions),
    mPushPower(push_power),
    mCloneStartTimes(clone_start_times),
    mKillRegrowTime(kill_regrow_times),
    mFathers(fathers),

    mSeed(seed),
    mSimId(0),

    mDisplayFrequency(0),
    mOutputPrefix("")


  {
    // Init space size vectors:
    mvSizeX.push_back(size_x);
    mvSizeY.push_back(size_y);
    mvSizeZ.push_back(size_z);

    // Find the maximum clone param vector length:
    std::vector <int> param_v_sizes;
    param_v_sizes.push_back(static_cast<int>(mMutationrates.size()));
    param_v_sizes.push_back(static_cast<int>(mBirthrates.size()));
    param_v_sizes.push_back(static_cast<int>(mDeathrates.size()));
    param_v_sizes.push_back(static_cast<int>(mAggressions.size()));
    param_v_sizes.push_back(static_cast<int>(mPushPower.size()));
    param_v_sizes.push_back(static_cast<int>(mCloneStartTimes.size()));
    param_v_sizes.push_back(static_cast<int>(mKillRegrowTime.size()));
    param_v_sizes.push_back(static_cast<int>(mFathers.size()));
    mNumberOfClones = *max_element(param_v_sizes.begin(), param_v_sizes.end());

    // Test if the length of all clone parameters is equal:
    for (std::vector<int>::size_type i = 0; i < param_v_sizes.size(); i++) {
      if (param_v_sizes[i] != 0 && param_v_sizes[i] != mNumberOfClones) {
        std::cerr << "Error during creation of SimulationParameterSet:\n";
        std::cerr << "  Clone parameter vectors not of same length.\n";
        exit(EXIT_FAILURE);
      }
    }

    // Init the mUniverses params with default 0:
    for (int i = 0; i < mNumberOfClones; i++) {
      mUniverses.push_back(0);
    }
  }

// Getters:
std::vector <int> SimulationParameterSet::SizeX() const {return mvSizeX;};
std::vector <int> SimulationParameterSet::SizeY() const {return mvSizeY;};
std::vector <int> SimulationParameterSet::SizeZ() const {return mvSizeZ;};
int SimulationParameterSet::NumberOfClones() const {return mNumberOfClones;};
int SimulationParameterSet::NumberOfSpaces() const {return mNumberOfSpaces;};
int SimulationParameterSet::Seed() const {return mSeed;};


int SimulationParameterSet::NumberClonalMutations() const {
  return mNumberClonalMutations;
}

double SimulationParameterSet::SequencingDepth() const {
  return mSequencingDepth;
}

double SimulationParameterSet::SequencingMinVaf() const {
  return mMinSequencingVaf;
}

int SimulationParameterSet::SequencingMinReads() const {
  return mMinSequencingReads;
}

int SimulationParameterSet::SequencingDepthModel() const {
  return mDepthModel;
}


std::vector <double> SimulationParameterSet::Mutationrates() const {
  return mMutationrates;
}
std::vector <double> SimulationParameterSet::Birthrates() const {
  return mBirthrates;
}
std::vector <double> SimulationParameterSet::Deathrates() const {
  return mDeathrates;
}
std::vector <double> SimulationParameterSet::Aggressions() const {
  return mAggressions;
}
std::vector <unsigned int> SimulationParameterSet::PushPower() const {
  return mPushPower;
}
std::vector <double> SimulationParameterSet::CloneStartTimes() const {
  return mCloneStartTimes;
}
std::vector <double> SimulationParameterSet::KillRegrowTime() const {
  return mKillRegrowTime;
}
std::vector <unsigned int> SimulationParameterSet::Fathers() const {
  return mFathers;
}
std::vector <unsigned int> SimulationParameterSet::Spaces() const {
  return mUniverses;
}
double SimulationParameterSet::DisplayFrequency() const {
  return mDisplayFrequency;
}
std::string SimulationParameterSet::OutputPrefix() const {return mOutputPrefix; };
int SimulationParameterSet::SimulationId() const {return mSimId;};


// Output functions:
void SimulationParameterSet::Print() const {
  // Print arguments:
  std::cout << std::endl;
  std::cout << "########## Options #############" << std::endl;
  for (int i = 0; i < mNumberOfSpaces; i++) {
    std::cout << "  Size " << i << ": "
              << mvSizeX[i] << "x" << mvSizeY[i] << "x" << mvSizeZ[i] << "\n";
  }
  // clone level data
  for (int i = 0; i < mNumberOfClones; i++) { // for each clone seperatly:
    std::cout << std::endl;
    std::cout << "  Mutation rate: " << mMutationrates[i] << std::endl;
    std::cout << "  Birth rate: " << mBirthrates[i] << std::endl;
    std::cout << "  Death rate: " << mDeathrates[i] << std::endl;
    std::cout << "  Aggression: " << mAggressions[i] << std::endl;
    std::cout << "  Push power: " << mPushPower[i] << std::endl;
    std::cout << "  Clone start time: " << mCloneStartTimes[i] << std::endl;
    std::cout << "  Kill & regrow time: " << mKillRegrowTime[i] << std::endl;
    std::cout << "  Father: " << mFathers[i] << std::endl;
    std::cout << "  Universe: " << mUniverses[i] << std::endl;

  } // end printing of clone level data
  std::cout << std::endl;
  std::cout << "  Clonal mutations: " << mNumberClonalMutations << "\n";
  if (mSequencingDepth > 0.0)
    std::cout << "  Sequencing depth: " << mSequencingDepth << "\n";
  if (mMinSequencingReads >= 0)
    std::cout << "  Min variant reads: " << mMinSequencingReads << "\n";
  if (mMinSequencingVaf >= 0)
    std::cout << "  Min VAF: " << mMinSequencingVaf << "\n";
  if (mDepthModel >= 0)
    std::cout << "  Sequencing depth model: " << mDepthModel << "\n";
  std::cout << std::endl;
  std::cout << "  Random seed: " << mSeed << std::endl;
  std::cout << "  Output prefix: " << mOutputPrefix << std::endl;
  std::cout << "################################" << std::endl;
  std::cout << std::endl;
}

bool SimulationParameterSet::SaveToFile(std::string out_file) const {
  std::ofstream output_stream;
  output_stream.open(out_file, std::ios::out | std::ios::trunc);

  if (output_stream.is_open()) {
    // First put the grid size:
    for (int i = 0; i < mNumberOfSpaces; i++) {
      output_stream << "# Universe size" << i << ": "
                    << mvSizeX[i] << "x"
                    << mvSizeY[i] << "x"
                    << mvSizeZ[i] << std::endl;
    }
    // Then other singlular params:
    output_stream << "# Seed: " << mSeed << std::endl;
    output_stream << "# N clones: " << mNumberOfClones << std::endl;
    output_stream << "# N clonal muts: " << mNumberClonalMutations << std::endl;
    if (mSequencingDepth > 0.0)
      output_stream << "# Sequencing depth: " << mSequencingDepth << std::endl;
    if (mMinSequencingReads >= 0)
      output_stream << "# Min variant reads: " << mMinSequencingReads << std::endl;
    if (mMinSequencingVaf >= 0)
      output_stream << "# Min VAF: " << mMinSequencingVaf << std::endl;
    if (mDepthModel >= 0)
      output_stream << "# Sequencing depth model: " << mDepthModel << std::endl;
    output_stream << "# Output prefix: " << mOutputPrefix << std::endl;

    // Now the header:
    char delim = '\t'; // TSV seperated.
    output_stream << "type_number" << delim
                  <<  "mutation_rate" << delim
                  <<  "birth_rate" << delim
                  <<  "death_rates" << delim
                  <<  "aggression" << delim
                  <<  "start_time" << delim
                  <<  "father" << delim
                  <<  "universe" << std::endl;

    // And now clone types line by line:
    for (int i = 0; i < mNumberOfClones; i++) {
      output_stream << i << delim // Column by column ...
                    << mMutationrates[i] << delim
                    << mBirthrates[i] << delim
                    << mDeathrates[i] << delim
                    << mAggressions[i] << delim
                    << mCloneStartTimes[i] << delim
                    << mFathers[i] << delim
                    << mUniverses[i] << std::endl;
    }
  } else { // failed to open out_file
    return false;
  }
  output_stream.close();
  return true;
}
