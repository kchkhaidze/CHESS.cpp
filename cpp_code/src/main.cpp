//#define DEBUG
#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

#define DISPLAY
#ifdef DISPLAY
#define DI(x) x
#else
#define DI(x)
#endif

#include <getopt.h>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "boost/multi_array.hpp"
#include <iostream>
#include <cassert>
#include "CImg.h"
#include <time.h>
#include "Phylogeny.hpp"
#include "Universe.hpp"
#include "Cell.hpp"
#include "CellType.hpp"



// Params of old function simulateTumorGrowth2D
#define SIZE_X 100 // size of the space [int]
#define SIZE_Y 100 // size of the space [int]
#define SIZE_Z 1 // size of the space [int]
#define MU 10          // mutation rate [float]
#define WTBR 1         // wildtype birth rate [float]
#define MUTBR 1.5        // mutant birth rate   [float]
#define CLST 5.0       // clone start time [float]
#define DISPFRQ 0.05   // display frequency [float]
#define ALPHA 0.0      // alpha [float]
#define BETA 1.0       // beta [float]
#define AGGR 1.0       // cell aggression [float]
#define CENTER(X) ((X + 1) / 2) - 1
#define OUTPUT_DIR "output"
#define SPACE_OUPUT_FILE_PREFIX "space"
#define TYPES_OUPUT_FILE_PREFIX "types"
#define PHYLO_OUPUT_FILE_PREFIX "phylo"
#define HISTORY_OUPUT_FILE_PREFIX "cell_number.csv"

int main(int argc, char* argv[]) {

  // Define possible arguments and init with their corresponding defaults:
  int size_x = SIZE_X; // size of the space [int]
  int size_y = SIZE_Y;  // size of the space [int]
  int size_z = SIZE_Z;  // size of the space [int]
  float mutation_rate = MU;       // mutation rate [int]
  float wildtype_birthrate = WTBR;
  float mutant_birthrate = MUTBR;
  float clone_start_time = CLST; //
  float display_frequency = DISPFRQ;
  float alpha = ALPHA;
  float beta = BETA;
  float aggression = AGGR;
  std::string output_dir = OUTPUT_DIR;
  time_t seed = time(NULL);

 // The possible arguments for getopts
  static struct option long_options[] = {
    {"size_x",  optional_argument, 0, 'x'},
    {"size_y",  optional_argument, 0, 'y'},
    {"size_z",  optional_argument, 0, 'z'},
    {"mutation_rate", optional_argument, 0, 'M'},
    {"wildtype_birthrate",  optional_argument, 0, 'w'},
    {"mutant_birthrate",  optional_argument, 0, 'm'},
    {"clone_start_time",    optional_argument, 0, 't'},
    {"display_frequency",    optional_argument, 0,  'f'},
    {"alpha",    optional_argument, 0,  'A'},
    {"beta",    optional_argument, 0,  'B'},
    {"aggression",    optional_argument, 0,  'R'},
    {"output_dir",    optional_argument, 0,  'o'},
    {"seed",    optional_argument, 0,  's'},
    {NULL, 0, NULL, 0}
  };

  // Parse commandline arguments (copied from getop man-page):
  while (true) {
    int c = getopt_long(argc, argv, "x:y:z:M:w:m:t:f:A:B:R:o:s:", long_options, NULL);

    if (c == -1)
      break;
    switch (c) {
      case 'x':
      size_x = atoi(optarg);
      break;
      case 'y':
      size_y = atoi(optarg);
      break;
      case 'z':
      size_z = atoi(optarg);
      break;
      case 'M':
      mutation_rate = atof(optarg);
      break;
      case 'w':
      wildtype_birthrate = atof(optarg);
      break;
      case 'm':
      mutant_birthrate = atof(optarg);
      break;
      case 't':
      clone_start_time = atof(optarg);
      break;
      case 'f':
      display_frequency = atof(optarg);
      break;
      case 'A':
      alpha = atof(optarg);
      break;
      case 'B':
      beta = atof(optarg);
      break;
      case 'R':
      aggression = atof(optarg);
      break;
      case 'o':
      output_dir = atof(optarg);
      break;
      case 's':
      seed = atoi(optarg);
      break;
      case ':':
        break;
      case '?':
        break;
        std::cerr <<  "Usage: " << argv[0] << " [please fill]" << std::endl;
        exit(EXIT_FAILURE);
      default:
        std::cout << "getopt returned character code 0" << c << std::endl;
    }
  }

  // Print arguments:
  std::cout << std::endl;
  std::cout << "########## Options #############" << std::endl;
  std::cout << "  Size: " << size_x << "x" << size_y << "x" << size_z << "\n";
  std::cout << "  Mutation rate: " << mutation_rate << std::endl;
  std::cout << "  Wildtype birth rate: " << wildtype_birthrate << std::endl;
  std::cout << "  Mutant birth rate: " << mutant_birthrate << std::endl;
  std::cout << "  Clone start time: " << clone_start_time << std::endl;
  std::cout << "  Alpha: " << alpha << std::endl;
  std::cout << "  Beta: " << beta << std::endl;
  std::cout << "  Aggression: " << aggression << std::endl;
  std::cout << std::endl;
  std::cout << "  Display frequency: " << display_frequency << std::endl;
  std::cout << "  Output dir: " << output_dir << std::endl;
  std::cout << "  Random seed: " << seed << std::endl;
  std::cout << "################################" << std::endl;

  // randomize seed
  srand(seed);

  // Display:
  DI(cimg_library::CImg<unsigned char> buffer(size_x, size_y, size_z, 3, 0);)
  DI(cimg_library::CImgDisplay disp(buffer, "Universe");)
  DI(double last_print_time = 0.0;)

  // Create the universe and a reusable pointer to handled cells:
  Universe universe(size_x, size_y, size_z);
  double dt = 0.0;
  Cell* p_cell;

  // Create two cell types (red and blue):
  CellType cell_type_blue(mutant_birthrate, alpha, beta, aggression,
                          mutation_rate, 0, 255, 0);
  CellType cell_type_red(mutant_birthrate, alpha, beta, aggression,
                         mutation_rate, 255, 0, 0);

  // Create and insert a new blue cell:
  p_cell = new Cell(&cell_type_blue);
  universe.InsertCell(CENTER(size_x), CENTER(size_y), CENTER(size_z), p_cell);

  // Use these vectors to track number of alive cells during each cyle:
  std::vector<float> track_time_gen;
  std::vector<int> track_cells_blue, track_cells_red;

  // Run till insertion of new cell type after CLST should occure:
  while(universe.Time() <= clone_start_time && !universe.LimitIsReached()) {

    // Print new universe to screen:
    DI(if (universe.Time() - last_print_time > display_frequency) {)
    DI(last_print_time = universe.Time();)
    DI(universe.Display(&buffer, &disp);)
    DI(})

    // Track current state of the universe:
    track_time_gen.push_back(universe.Time());
    track_cells_blue.push_back(cell_type_blue.NumMembers());
    track_cells_red.push_back(cell_type_red.NumMembers());

    // Select reaction type and members:
    universe.NextReactionType(dt)->RandomMember()->Divide();
    universe.IncrementTimeBy(dt);
  } // stop running after reaching CLST



  // Make conversion of a random blue cell:
  D(std::cout << "########## Comment #############" << std::endl;)
  D(std::cout << " Conversion of one blue -> red" << std::endl;)
  D(std::cout << "################################" << std::endl;)
  cell_type_blue.RandomMember()->Type(&cell_type_red);



  // Run
  while(!universe.LimitIsReached()) {

    // Print new universe to screen:
    DI(if (universe.Time() - last_print_time > display_frequency) {)
    DI(last_print_time = universe.Time();)
    DI(universe.Display(&buffer, &disp);)
    DI(})

    // Track current state of the universe:
    track_time_gen.push_back(universe.Time());
    track_cells_blue.push_back(cell_type_blue.NumMembers());
    track_cells_red.push_back(cell_type_red.NumMembers());

    // Select reaction type and members:
    universe.NextReactionType(dt)->RandomMember()->Divide();
    universe.IncrementTimeBy(dt);
  } // stop running after reaching limit of universe


  cell_type_blue.RandomMember()->AssociatedNode()->PrintAncestry();


  // Write all results to output files:
  std::cout << "Dumping output" << std::endl;
  universe.SpaceToCsvFile(output_dir + SPACE_OUPUT_FILE_PREFIX);
  universe.TypesToCsvFile(output_dir + TYPES_OUPUT_FILE_PREFIX);
  universe.PhylogeniesToFile(output_dir + PHYLO_OUPUT_FILE_PREFIX);

  // Write history to another output file:
  std::ofstream output_stream;
  output_stream.open(output_dir + HISTORY_OUPUT_FILE_PREFIX,
                     std::ios::out | std::ios::trunc);

  if (output_stream.is_open()) {
    output_stream << "time,blue_cells,red_cells" << std::endl;
    for(std::vector<float>::size_type i = 0; i < track_time_gen.size(); i++) {
      output_stream << track_time_gen[i] << ","
        << track_cells_blue[i] << "," << track_cells_red[i] << std::endl;
    }
    output_stream.close();
  } else {
    std::cerr << "Couldn't write histroy to csv file." << std::endl;
  }
  return EXIT_SUCCESS;
}
