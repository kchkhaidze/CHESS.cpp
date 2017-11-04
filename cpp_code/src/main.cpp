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

#include <random>
#include <cmath>
#include <iostream>
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
#define SIZE_X 500 // size of the space [int]
#define SIZE_Y 500 // size of the space [int]
#define SIZE_Z 1 // size of the space [int]
#define MU 10          // mutation rate [int]
#define SD 1           // [int]
#define WTBR 1         // wildtype birth rate [int]
#define MUTBR 1.5        // mutant birth rate   [int]
#define CLST 5.0       // clone start time [float]
#define DISPFRQ 0.05   // display frequency [float]
#define ALPHA 0.0      // alpha [float]
#define BETA 1.0       // beta [float]
#define AGGR 1.0       // cell aggression [float]
#define CENTER(X) ((X + 1) / 2) - 1
#define SPACE_OUPUT_FILE_PREFIX "output/space"
#define TYPES_OUPUT_FILE_PREFIX "output/types"

int main() {

  // randomize seed
  time_t seed = time(NULL);
  std::cerr << "Random seed: " << seed << std::endl;
  srand(seed);

  // Display:
  DI(cimg_library::CImg<unsigned char> buffer(SIZE_X,SIZE_Y,SIZE_Z,3,0);)
  DI(cimg_library::CImgDisplay disp(buffer, "Universe");)
  DI(double last_print_time = 0.0;)

  // Create the universe and a reusable pointer to handled cells:
  Universe universe(SIZE_X, SIZE_Y, SIZE_Z);
  double dt = 0.0;
  Cell* p_cell;

  // Create two cell types (red and blue):
  CellType cell_type_blue(WTBR, ALPHA, BETA, AGGR, 0, 255, 0);
  CellType cell_type_red(MUTBR, ALPHA, BETA, AGGR, 255, 0, 0);

  // Create and insert a new blue cell:
  p_cell = new Cell(&cell_type_blue);
  universe.InsertCell(CENTER(SIZE_X), CENTER(SIZE_Y), CENTER(SIZE_Z), p_cell);

  // Use these vectors to track number of alive cells during each cyle:
  std::vector<int> track_time_gen, track_cells_blue, track_cells_red;


  // Run till insertion of new cell type after CLST should occure:
  while(universe.Time() <= CLST && !universe.LimitIsReached()) {

    // Print new universe to screen:
    DI(if (universe.Time() - last_print_time > DISPFRQ) {)
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
    DI(if (universe.Time() - last_print_time > DISPFRQ) {)
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


  // Write all results to output files:
  std::cerr << "Dumping output" << std::endl;
  universe.SpaceToCsvFile(SPACE_OUPUT_FILE_PREFIX);
  universe.TypesToCsvFile(TYPES_OUPUT_FILE_PREFIX);

  //
  track_time_gen.push_back(universe.Time());
  track_cells_blue.push_back(cell_type_blue.NumMembers());
  track_cells_red.push_back(cell_type_red.NumMembers());


  return EXIT_SUCCESS;
}
