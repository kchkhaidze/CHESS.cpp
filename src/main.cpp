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

#ifdef _DEBUG_
#define D(x) x
#else
#define D(x)
#endif

#ifdef _DISPLAY_
#define DI(x) x
#else
#define DI(x)
#endif


#include "Universe.hpp"
#include "Shapes.hpp"
#include "extern_global_variables.hpp"
#include "SimulationParameterSet.hpp"

// Params of old function simulateTumorGrowth2D
#define WVAF_OUTPUT_FILE_PREFIX "vaf_wholetumour"
#define BVAF_OUTPUT_FILE_PREFIX "vaf_bulk"
#define SAMPGRID_OUTPUT_FILE_PREFIX "samp_grid"


boost::random::mt19937_64 rng;


int main(int argc, char* argv[]) {

  // Parse command line sarguments:
  SimulationParameterSet sim_params(argc, argv);


  // General simulation:
  Universe universe(sim_params);
  universe.RunSimulation(true);
  universe.SaveResults();


  // Sampling model:
  for (int i = 0; i < sim_params.NumberOfSpaces(); i++) {
    Shape* square = new Box(0, 0, 0,
                            sim_params.SizeX()[i]-1,
                            sim_params.SizeY()[i]-1,
                            sim_params.SizeZ()[i]-1);

    // This is the sample:
    SequencingResult sample = universe.TakeSample(i, square).Sequence(sim_params);
    sample.Hist();

    // Save sampled whole tumour vaf to file:
    sample.SaveToFile(sim_params.OutputPrefix() + WVAF_OUTPUT_FILE_PREFIX
                       + "_" + std::to_string(sim_params.SimulationId())
                       +".csv");
  }




  // Wait for user to close output screen.
  DI(universe.WaitTillDisplayClosed();)

  return EXIT_SUCCESS;
}
