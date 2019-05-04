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

#define CENTER(X) ((X + 1) / 2) - 1
#define KILL_REGROW_PERCENT 0.99

#define SPACE_OUPUT_FILE_PREFIX "space"
#define GENOT_OUPUT_FILE_PREFIX "space_generation"
#define TYPES_OUPUT_FILE_PREFIX "types"
#define PHYLO_OUPUT_FILE_PREFIX "phylo"
#define HISTORY_OUPUT_FILE_PREFIX "timestep_cellcount"
#define DIPLAY_FILE "space_image"
#define PARAMS_OUTPUT_FILE_PREFIX "sim_params"

#include "extern_global_variables.hpp"
#include "CellType.hpp"
#include "Shapes.hpp"
#include "Cell.hpp"
#include "Universe.hpp"
#include "Phylogeny.hpp"
#include <stdlib.h>
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <cmath>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include "HSL2RGB.hpp"

#ifdef _DISPLAY_

void init_display(double dist, bool three_d, int max_dim, int n_space) {

  if (three_d) {
    /* Enable one positional light. */
    GLfloat light_diffuse1[] = {1.0, 1.0, 1.0, 1.0};  /* White light. */
    GLfloat light_specular1[] = {1.0, 1.0, 1.0, 1.0};  /* White specular light. */
    GLfloat light_position1[] = {1000.0, 600.0, 600.0, 0.0};  /* Positional light. */
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);
    glLightfv(GL_LIGHT1, GL_POSITION, light_position1);

    /* Enable one directional light. */
    GLfloat light_diffuse2[] = {0.2, 0.2, 0.2, 1.0};  /* Dim white light. */
    GLfloat light_specular2[] = {0, 0, 0, 0};  /* No specular light. */
    GLfloat light_position2[] = {0.0, 0.0, 1.0, 1.0};  /* Infinite light along z.*/
    glLightfv(GL_LIGHT2, GL_DIFFUSE, light_diffuse2);
    glLightfv(GL_LIGHT2, GL_SPECULAR, light_specular2);
    glLightfv(GL_LIGHT2, GL_POSITION, light_position2);

    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);
  } else {
    /* Enable one directional light. */
    GLfloat light_diffuse3[] = {1.0, 1.0, 1.0, 1.0};  /* Dim white light. */
    GLfloat light_specular3[] = {0.0, 0.0, 0.0, 0.0};  /* No specular light. */
    GLfloat light_position3[] = {0.0, 0.0, 1.0, 1.0};  /* Infinite light along z.*/
    glLightfv(GL_LIGHT3, GL_DIFFUSE, light_diffuse3);
    glLightfv(GL_LIGHT3, GL_SPECULAR, light_specular3);
    glLightfv(GL_LIGHT3, GL_POSITION, light_position3);

    glEnable(GL_LIGHT3);
  }
  glEnable(GL_LIGHTING);


  /* Use depth buffering for hidden surface elimination. */
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_SMOOTH);

  /* Setup the view of the cube. */
  glMatrixMode(GL_PROJECTION);
  gluPerspective( /* field of view in degree */ 40.0,
    /* aspect ratio */ n_space, /* Z near */ 2.0, /* Z far */ 5000.0);
  glMatrixMode(GL_MODELVIEW);

  gluLookAt(max_dim * n_space / 2, 0.0, dist,  /* eye at (0,0,dist) */
    max_dim * n_space / 2, max_dim / 2, max_dim / 2,
    0.0, 1.0, 0.0);          /* up is in positive Y direction */

  //if (three_d) {
    /* Adjust position to be asthetic angle. */
    //glTranslatef(0.0, 0.0, -1.0);
    //glRotatef(60, 1.0, 0.0, 0.0);
    //glRotatef(-20s 0.0, 0.0, 1.0);
  //}
}

#endif // _DISPLAY_

// Universe ////////////////////////////////////////////////////////////////////

// Constructor:
Universe::Universe(SimulationParameterSet params)
  : mParameters(params),
    mTime(0.0L),
    mLimitReached(false)
    #ifdef _DISPLAY_
      /* continue with optional params if display mode is on */,
      mLastPrintTime(0.0)
    #endif // _DISPLAY_
{

  // Construct the space vector:
  std::vector <int> vx = params.SizeX();
  std::vector <int> vy = params.SizeY();
  std::vector <int> vz = params.SizeZ();
  for (int n = 0; n < vx.size(); n++) {
    boost::multi_array<Cell*, 3> space(boost::extents[vx[n]][vy[n]][vz[n]]);
    // Make the space empty:
    for(size_t i = 0; i < space.shape()[0]; i++) {
      for(size_t j = 0; j < space.shape()[1]; j++) {
        for(size_t k = 0; k < space.shape()[2]; k++) {
          space[i][j][k] = 0;
        }
      }
    }
    mvSpace.push_back(space);
  }

  #ifdef _DISPLAY_



    // Only create a display if display freq is positive:
    int m = 0;
    glutInitWindowSize(300 * params.NumberOfSpaces(), 300);
    glutInit(&m, 0);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutCreateWindow("Universe");

    // Determine a good spectator distance.
    std::vector <int> max_dims;
    max_dims.push_back(*max_element(vx.begin(), vx.end()));
    max_dims.push_back(*max_element(vy.begin(), vy.end()));
    max_dims.push_back(*max_element(vz.begin(), vz.end()));
    mMaxSpaceDim = *max_element(max_dims.begin(), max_dims.end());

    max_dims.push_back(150);
    int view_dist = *max_element(max_dims.begin(), max_dims.end());
    //view_dist *= params.NumberOfSpaces();

    // Init the display:
    if (mvSpace[0].shape()[2] > 1) {
      init_display(view_dist * 1.5, true, mMaxSpaceDim, params.NumberOfSpaces());
    } else {
      init_display(view_dist * 1.5, false, mMaxSpaceDim, params.NumberOfSpaces());
    }

    // Display:
    glFlush();
    this->Display();
    glutPostRedisplay();

  #endif // _DISPLAY_
}

// Destructor:
Universe::~Universe(){
  D(std::cout << "Universe::~Universe() " << this << std::endl;)

  // Remove all Cells:
  for (int n = 0; n < mvSpace.size(); n++) {
    for(int x = 0; x < mvSpace[n].shape()[0]; x++) {
      for(int y = 0; y < mvSpace[n].shape()[1]; y++) {
        for(int z = 0; z < mvSpace[n].shape()[2]; z++) {
          if (this->ContainsCell(n, x, y, z)) {
            delete mvSpace[n][x][y][z];
          }
        }
      }
    }
  }


  // Remove all CellTypes:
  for(std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++) {
    delete mpTypes[i];
  }

  // Remove the complete phylogeny:
  for(std::vector<PhylogenyRoot*>::size_type i = 0;
     i < mpPhylogenies.size(); i++) {
    delete mpPhylogenies[i];
  }
}


// Getter functions:
long double Universe::Time() const {return mTime;};

bool  Universe::LimitIsReached() const {return mLimitReached;};

bool  Universe::ContainsCell(int i, int x, int y, int z) const {
  // Returns true if the position i, x, y, z contains a cell.
  return mvSpace[i][x][y][z] != 0;
}

bool Universe::ContainsCoordinate(int i, int x, int y, int z) const {
  // Returns true if the position x, y, z is within the limits of the universe.

  return x < mvSpace[i].shape()[0] &&
         y < mvSpace[i].shape()[1] &&
         z < mvSpace[i].shape()[2] &&
         x >= 0 &&
         y >= 0 &&
         z >= 0;
}

void  Universe::Size(int n, int& x, int& y, int& z) const {
  x = mvSpace[n].shape()[0];
  y = mvSpace[n].shape()[1];
  z = mvSpace[n].shape()[2];
}

Cell *Universe::GetCell(int i, int x, int y, int z) {
  return mvSpace[i][x][y][z];
}

Cell* Universe::RemoveCell (Cell* pCell) {
  return this->RemoveCell(pCell->N(), pCell->X(), pCell->Y(), pCell->Z());
}

Cell* Universe::RemoveCell (int i, int x, int y, int z) {
  // Returns a cell after removal from the space of the universe.

  Cell* p_removed_cell = mvSpace[i][x][y][z];
  mvSpace[i][x][y][z] = 0;
  p_removed_cell->Location(0, 0, 0, 0);
  return p_removed_cell;
}

CellType* Universe::NextReactionType(long double *r_delta_time, int *action){
  // Samples and returns the next reaction type that occures.

  // Return null if the universe contains no types:
  if (mpTypes.size() == 0) {
    std::cerr << "Error in Universe::NextReactionType(double *, int *):\n";
    std::cerr << "   > No types to choose from!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Debug messages:
  D(std::cout << std::endl;)
  D(std::cout << "########## Sampling ###########" << std::endl;)
  D(std::cout << "  Next reaction type:" << std::endl;)
  D(std::cout << std::endl;)

  // Store reaction specs in these vectors:
  std::vector<long double> deltas;
  std::vector<int> reaction_types;
  std::vector<std::vector<CellType*>::size_type> cell_types;
  D(int n_react = 0;)

  boost::uniform_real<double> runi_dist(0.0, 1.0);

  // Calculate dt for all types in the universe:
  for(std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++) {
    for(int reaction_type = 1; reaction_type <= 2; reaction_type++) {

      // :ToDo: Check this is correct.
      double runi = runi_dist(rng);
      double s_number = mpTypes[i]->NumMembers();

      double reaction_rate;
      switch(reaction_type){
        case 1:
          reaction_rate = mpTypes[i]->Deathrate();
          break;
        case 2:
          reaction_rate = mpTypes[i]->Birthrate();
          break;
        default: //undefined action requested
          std::cerr << "Error in Universe::NextReactionType:\n";
          std::cerr << "   > Action " <<  reaction_type << " undefined.\n";
          std::cerr << std::endl;
          exit(EXIT_FAILURE);
      }

      long double delta_t = (log(1) - log(runi)) / (s_number * reaction_rate);

      // Debug messages:
      D(std::cout << "    Candidate reaction: " << ++n_react << std::endl;)
      D(std::cout << "      Species: " << i << std::endl;)
      D(std::cout << "      Species ID: " << mpTypes[i]->Id() << std::endl;)
      D(std::cout << "      Number of members : " << s_number << std::endl;)
      D(std::cout << "      Reaction type: " << reaction_type << std::endl;)
      D(std::cout << "      Reaction rate: " << reaction_rate << std::endl;)
      D(std::cout << "      Random uint: " << runi << std::endl;)
      D(std::cout << "      Delta t: " << delta_t << std::endl;)
      D(std::cout << std::endl;)

      // Collect candidate reactions:
      deltas.push_back(delta_t);
      reaction_types.push_back(reaction_type);
      cell_types.push_back(i);
    }
  }

  // Find minimum delta t of all reactions:
  int i_min;
  long double delta_time_minimum;
  for(std::vector<double>::size_type i = 0; i < deltas.size(); i++){
    if(i == 0 || delta_time_minimum > deltas[i]){
      i_min = i;
      delta_time_minimum = deltas[i];
    }
  }

  // Debug messages:
  D(std::cout << "  Next reaction:" << std::endl;)
  D(std::cout << "    Reaction id: " << reaction_types[i_min] << std::endl;)
  D(std::cout << "    Species id: " << cell_types[i_min] << std::endl;)
  D(std::cout << "    dt: " << deltas[i_min] << std::endl;)
  D(std::cout << "###############################" << std::endl;)

  // Return results:
  *r_delta_time = deltas[i_min];
  *action = reaction_types[i_min];
  return(mpTypes[cell_types[i_min]]);
}

std::vector< std::array<int, 3> > Universe::FreeNeighbours(int i, int x, int y, int z) const {
  std::vector< std::array<int, 3> > result;
  std::array<int, 3> cords;

  int xstart, xend;
  if (mvSpace[i].shape()[0] > 1) {
    xstart = x - 1;
    xend = x + 1;
  } else {
    xstart = x;
    xend = x;
  }

  int ystart, yend;
  if (mvSpace[i].shape()[1] > 1) {
    ystart = y - 1;
    yend = y + 1;
  } else {
    ystart = y;
    yend = y;
  }

  int zstart, zend;
  if (mvSpace[i].shape()[2] > 1) {
    zstart = z - 1;
    zend = z + 1;
  } else {
    zstart = z;
    zend = z;
  }

  // Get all free neighbours:
  for (int xn = xstart; xn <= xend; xn++) {
    for (int yn = ystart; yn <= yend; yn++) {
      for (int zn = zstart; zn <= zend; zn++) {
        if (this->ContainsCoordinate(i, xn, yn, zn) &&
              !this->ContainsCell(i, xn, yn, zn) &&
                (xn != x || yn != y || zn != z))
        {
          cords[0] = xn; cords[1] = yn; cords[2] = zn;
          result.push_back(cords);
        }
      }
    }
  }
  return result;
}

// get total cell counts per cell type
std::vector <unsigned long> Universe::CellCountsPerType() const {
  std::vector <unsigned long> result;
  for(int i = 0; i < mvCellCounts.size(); i++){
    result.push_back(mvCellCounts[i].back());
  }
  return result;
}

// Sampling functions (getters):
SingleCellSample Universe::TakeSingleCellSample(int n, Shape* pShape){

  // Ints to store the next cell to pick:
  int x, y, z;

  // Stage all samples for sequencing:
  while (pShape->next_coordinate(x, y, z)) { // Get next location from shape.
    if (mvSpace[n][x][y][z] != 0) { //selected location is not empty
        mvSpace[n][x][y][z]->AssociatedNode()->StageNodeForSequencing();
    }
  }

  //Now sequence the staged cells:
  SingleCellSample sample;
  std::vector<int> genotype_stack;

  for (std::vector<PhylogenyRoot*>::size_type i=0; i<mpPhylogenies.size(); i++){
  // for all stored phylo indices do:
    mpPhylogenies[i]->Root()->SampleStagedSingleCells(sample, genotype_stack);
    mpPhylogenies[i]->Root()->UnstageNodes();
  }

  return sample;
}

// Sampling functions (getters):
CellSample Universe::TakeSample(int n, Shape* pShape){

  // Counter for the number of sampled nodes:
  unsigned long int sampled_cells_cnt = 0;

  // Ints to store the next cell to pick:
  int x, y, z;

  // Stage all samples for sequencing:
  while (pShape->next_coordinate(x, y, z)) { // Get next location from shape.
    if (mvSpace[n][x][y][z] != 0) { //selected location is not empty
      // Increment counters of cells:
      sampled_cells_cnt++;
      mvSpace[n][x][y][z]->AssociatedNode()->StageNodeForSequencing();
    }
  }

  // The resulting sample will be stored here:
  CellSample sample;
  sample.SetSampledCellNumber(sampled_cells_cnt);

  // Now sequence the staged cells:
  for (std::vector<PhylogenyRoot*>::size_type i = 0;
       i < mpPhylogenies.size(); i++)
  {
    mpPhylogenies[i]->Root()->SequenceStagedNodes(sample);
    mpPhylogenies[i]->Root()->UnstageNodes();
  }

  return sample;
}

std::vector <CellSample> Universe::TakeMultipleSamples(int sample_strategy) {

  switch(sample_strategy){
    case 1: // Take squares low depth (20X)
      // ...
      break;
    case 2: // Take squares high depth (100X)
      // ...
      //TakeSquares(radius, sim_id, size_x, size_y, size_z, depth, min_reads,
      //  pois_depth, sampgrid_output_file, bulk_output_file);
      break;
    case 20:
      // to do: TakeStripes
    break;
  }

  std::vector <CellSample> pEmptyVector;
  return pEmptyVector;
}

std::vector <std::string> Universe::SingleCellTrees(int n, Shape *pShape) {

  std::vector <std::string> result;
  std::stringstream outstream;

  // Ints to store the next cell to pick:
  int x, y, z;

  // Stage all samples for sequencing:
  while (pShape->next_coordinate(x, y, z)) { // Get next location from shape.

    mvSpace[n][x][y][z]->AssociatedNode()->StageNodeForSequencing();

  }
  // Now sequence the staged cells:
  for (std::vector<PhylogenyRoot*>::size_type i = 0;
       i < mpPhylogenies.size(); i++)
  {
    outstream.str("");
    outstream << "tree" << i;
    mpPhylogenies[i]->Root()->StagedNodesToStream(outstream);
    outstream << ";";
    mpPhylogenies[i]->Root()->UnstageNodes();
    result.push_back(outstream.str());
  }

  return result;
}


// Setter functions:
void Universe::IncrementTimeBy(long double Delta){ mTime += Delta; }

void Universe::MarkLimitIsReached() { mLimitReached = true; }

bool Universe::InsertCell(int n, int x, int y, int z, Cell* pCell) {
  return InsertCell(n, x, y, z, pCell, true);
}

bool Universe::InsertCell(int n, int x, int y, int z, Cell* pCell,
                          bool is_new_lineage){

  // Debug messages:
  D(std::cout << std::endl;)
  D(std::cout << "########## Insertion ##########" << std::endl;)
  D(std::cout << "   ID: " << pCell->Id() << std::endl;)
  D(std::cout << "   Location:" << std::endl;)
  D(std::cout << "       Universe: " << this << std::endl;)
  D(std::cout << "       x: " << x << std::endl;)
  D(std::cout << "       y: " << y << std::endl;)
  D(std::cout << "       z: " << z << std::endl;)
  D(std::cout << "###############################" << std::endl;)

  // Dont allow inserts outside of the array extends:
  if (this->ContainsCoordinate(n,x,y,z) == false) {
    std::cout << "Error: Insertion of cell failed:" << std::endl;
    std::cout << "Error: Coordinates (i+1) were: " << n << x + 1 << ",";
    std::cout << y + 1 << "," << z + 1 << std::endl;
    std::cout << "Error:  Universe size was: "
              << mvSpace[n].shape()[0] << ","
              << mvSpace[n].shape()[1] << ","
              << mvSpace[n].shape()[2] << std::endl;
    //exit EXIT_FAILURE;
    return false;
  }

  // Dont allow inserts on top of cell:
  if (this->ContainsCell(n,x,y,z)) {
    std::cout << "Error: Insertion of cell failed:" << std::endl;
    std::cout << "Error: Coordinates (i) were: " << x << ",";
    std::cout << y << "," << z  << std::endl;
    std::cout << "Error:  Position was taken by" << mvSpace[n][x][y][z] << std::endl;
    exit(EXIT_FAILURE);
    //return false;
  }

  // Insert type into universe
  if (is_new_lineage) {
    this->RegisterType(pCell->Type());
  }

  // Update position and universe of the cell:
  pCell->Location(this, n, x, y, z);

  // Update universe:
  mvSpace[n][x][y][z] = pCell;

  // Register new phylogeny:
  if (pCell->AssociatedNode() == 0) {
    PhylogenyRoot* pNewPhylo = new PhylogenyRoot(pCell);
    mpPhylogenies.push_back(pNewPhylo);
  }

  return true;
}

void Universe::RegisterType(CellType* p_new_type){
  std::vector<CellType *>::size_type i = 0;

  while(i < mpTypes.size()) {
    D(std::cout << "Comparing new cell type " <<  p_new_type;)
    D(std::cout << " with " << mpTypes[i] << std::endl;)
    if (mpTypes[i] == p_new_type)
      break;
    i++;
  }

  if (i == mpTypes.size()) { // reached end.
    D(std::cout << "Registering new cell type " <<  p_new_type << std::endl;)

    // Insert type itself:
    mpTypes.push_back(p_new_type);

    // Recreate the history record up to the present date:
    std::vector<unsigned long> v_new_history;
    v_new_history.resize(mvHistoryTime.size());
    std::fill(v_new_history.begin(), v_new_history.end(), 0UL);
    mvCellCounts.push_back(v_new_history);
  }
}

void Universe::RecordToHistory() {
  for (int i = 0; i < mpTypes.size(); i++) {
    mvCellCounts[i].push_back(mpTypes[i]->NumMembers());
  }
  mvHistoryTime.push_back(mTime);
}

// Output Functions:
void Universe::Print() const {

  for (int n = 0; n < mvSpace.size(); n++) {
    std::cout << std::endl;
    std::cout << "############ Space " << n << " ############" << std::endl;
    std::cout << std::endl;

    for(int k = mvSpace[n].shape()[2] - 1; k >= 0; k--) {
      std::cout << "  z " << k << ":" << std::endl;
      for(int j = mvSpace[n].shape()[1] - 1; j >= 0; j--) {
        for(int i = 0; i < mvSpace[n].shape()[0]; i++) {
          std::cout << " " <<  mvSpace[n][i][j][k]->Type()->Id() << std::endl;
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
  }
  std::cout << "###############################" << std::endl;
  std::cout << std::endl;
}


void Universe::SpacesToCsvFiles(std::string outfile_prefix) const {
  for (int n = 0; n < mvSpace.size(); n++) {
    std::string outfile = outfile_prefix + std::to_string(n) + ".csv";
    this->SpaceToCsvFile(n, outfile);
  }
}

void Universe::SpaceToCsvFile(std::string outfile_name) const {
  this->SpaceToCsvFile(0, outfile_name);
}

void Universe::SpaceToCsvFile(int n, std::string outfile_name) const {

  for(int z = 0; z < mvSpace[n].shape()[2]; z++) {
    std::ofstream output_stream;
    output_stream.open(outfile_name, std::ios::out | std::ios::trunc);

    if (output_stream.is_open()) {
      for(int y = mvSpace[n].shape()[1] - 1; y >= 0; y--) {
        for(int x = 0; x < mvSpace[n].shape()[0]; x++) {
          if (this->ContainsCell(n,x,y,z)) {
            output_stream << std::to_string(mvSpace[n][x][y][z]->Type()->Id());
          } else {
            output_stream <<  "NA";

          }
          if (x < mvSpace[n].shape()[0] - 1) {
            output_stream << ",";
          }
        }
        output_stream << std::endl;
      }
    } else {
      std::cerr << "Couldn't write space to csv file." << std::endl;
    }
    output_stream.close();
  }
}

void Universe::SpaceGenerationToFile(int n, std::string outfile_name) const {

  for(int z = 0; z < mvSpace[n].shape()[2]; z++) {
    std::ofstream output_stream;
    output_stream.open(outfile_name, std::ios::out | std::ios::trunc);

    if (output_stream.is_open()) {
      for(int y = mvSpace[n].shape()[1] - 1; y >= 0; y--) {
        for(int x = 0; x < mvSpace[n].shape()[0]; x++) {
          if (this->ContainsCell(n,x,y,z)) {
            output_stream << std::to_string(
              mvSpace[n][x][y][z]->AssociatedNode()->Generation());
          } else {
            output_stream <<  "NA";

          }
          if (x < mvSpace[n].shape()[0] - 1) {
            output_stream << ",";
          }
        }
        output_stream << std::endl;
      }
    } else {
      std::cerr << "Couldn't write space to csv file." << std::endl;
    }
    output_stream.close();
  }
}

void Universe::SampGridToCsvFile(std::string out_prefix,
                                  boost::multi_array<int, 3> samp_grid,
                                  int sim_id) const {

  std::string suffix = ".csv";
  std::ofstream output_stream;
  std::string outfile_name;

  outfile_name = out_prefix + "_" + std::to_string(sim_id) + suffix;
  output_stream.open(outfile_name, std::ios::out | std::ios::trunc);

   if (output_stream.is_open()) {
     for(int y = samp_grid.shape()[1] - 1; y >= 0; y--) {
       for(int x = 0; x < samp_grid.shape()[0]; x++) {
         output_stream << std::to_string(samp_grid[x][y][0]);
         if (x < samp_grid.shape()[0] - 1) {
           output_stream << ",";
         }
       }
       output_stream << std::endl;
     }
   } else {
     std::cerr << "Couldn't write space to csv file." << std::endl;
   }
  output_stream.close();
}


void Universe::TypesToCsvFile(std::string outfile_name) const {

  std::ofstream output_stream;
  output_stream.open(outfile_name, std::ios::out | std::ios::trunc);

  if (output_stream.is_open()) {
    output_stream << "ID,birthrate,aggression,mu,pushpower,N\n";

    for(std::vector<CellType *>::size_type i = 0; i < mpTypes.size(); i++) {
      output_stream << std::to_string(mpTypes[i]->Id()) << ","
        << std::setprecision(6) << std::fixed
        << mpTypes[i]->Birthrate() << ","
        << mpTypes[i]->Aggression() << ","
        << mpTypes[i]->Mu() << ","
        << mpTypes[i]->PushPower() << ","
        << mpTypes[i]->NumMembers() << std::endl;
    }
    output_stream.close();
  } else {
    std::cerr << "Couldn't write types to csv file." << std::endl;
  }
}

void Universe::PhylogeniesToFile(std::string out_prefix) const {
  std::string suffix = ".tree";
  std::ofstream output_stream;
  std::string outfile_name;

  for(std::vector<PhylogenyRoot*>::size_type i=0; i<mpPhylogenies.size(); i++) {

      outfile_name = out_prefix + suffix;
      output_stream.open(outfile_name, std::ios::out | std::ios::trunc);

      if (output_stream.is_open()) {
        output_stream << (*mpPhylogenies[i]->Root()) << ";" << std::endl;
      } else {
        std::cerr << "Couldn't write space to csv file." << std::endl;
      }
    output_stream.close();
  }
}

bool Universe::HistoryToFile(std::string output_file) const {

  // Open file:
  std::ofstream output_stream;
  output_stream.open(output_file, std::ios::out | std::ios::trunc);

  if (output_stream.is_open()) { // Guess what ... now we write to the file!

    // Start with a header:
    output_stream << "time";
    for (int i = 0; i < mvCellCounts.size(); i++) {
      output_stream << "," << "type" << i + 1;
    }
    output_stream << std::endl;

    // then write the cell counts:
    for (int i = 0; i < mvHistoryTime.size(); i++) {
      output_stream << mvHistoryTime[i];
      for (int j = 0; j < mvCellCounts.size(); j++) {
        output_stream << "," << mvCellCounts[j][i];
      }
      output_stream << std::endl;
    }

    // Close the stream. Done!
    output_stream.close();

  } else { // Be kind and let people know if you couldn't open the file ...
    std::cerr << "Error in function Universe::HistoryToFile:" << std::endl;
    std::cerr << "   Couldn't open output file." << std::endl;
    std::cerr << "   " << output_file << std::endl;
    return false;
  }
  return true;
}


// Other functions:
bool Universe::PushCell(int n, double newX, double newY, double newZ,
                        double dx, double dy, double dz,
                        unsigned int push) {


  // Determine integer coordinates of current cell:
  int intCurX = round(newX);
  int intCurY = round(newY);
  int intCurZ = round(newZ);
  int intNewX;
  int intNewY;
  int intNewZ;

  // Increase till new integer coords jumps:
  do {
    newX += dx;
    newY += dy;
    newZ += dz;
    intNewX = round(newX);
    intNewY = round(newY);
    intNewZ = round(newZ);
  } while (intNewX == intCurX && intNewY == intCurY && intNewZ == intCurZ);


  // If you hit a border stop pushes.
  if (this->ContainsCoordinate(n, intNewX, intNewY, intNewZ) == false) {
    D(std::cout << "  Hit border during push." << std::endl;)
    //this->MarkLimitIsReached();
    return false;
  }

  // Test if next space is free.
  if (this->ContainsCell(n, intNewX, intNewY, intNewZ) == false ||
      (push > 0 && PushCell(n, newX, newY, newZ, dx, dy, dz, push - 1)))
  {
    // Debug messages:
    D(std::cout << "  Push sucessfull:" << std::endl;)
    D(std::cout << "    id: " << mvSpace[n][intCurX][intCurY][intCurZ] << std::endl;)
    D(std::cout << "    x: " << intCurX << "->" << intNewX << std::endl;)
    D(std::cout << "    y: " << intCurY << "->" << intNewY << std::endl;)
    D(std::cout << "    z: " << intCurZ << "->" << intNewZ << std::endl;)
    D(std::cout << "    push: " << push << "->" << intNewZ << std::endl;)
    D(std::cout << std::endl;)

    // Move cell:
    mvSpace[n][intNewX][intNewY][intNewZ] = mvSpace[n][intCurX][intCurY][intCurZ];
    mvSpace[n][intCurX][intCurY][intCurZ]->Location(intNewX, intNewY, intNewZ);
    mvSpace[n][intCurX][intCurY][intCurZ] = 0;
    return true;
  }

  return false;
}

// Simulate tumours:
bool Universe::RunSimulation(bool verbose)
{

  if (verbose) { // Print arguments:
    mParameters.Print();
  } // end of argument printing:

  // Set random seed
  rng.seed(mParameters.Seed());

  // Basic checks of the input:
  double last_start_time = 0.0;
  for (int i = 0; i < mParameters.NumberOfClones(); i++) {
    // Ensure that the clone start times are ordered:
    if (mParameters.CloneStartTimes()[i] < last_start_time) {
      std::cerr << "Defined clone start times have to be sorted!" << std::endl;
      return 0;
    }
    last_start_time = mParameters.CloneStartTimes()[i];

    // Ensure that fathers are defined:
    if ( (i == 0 && mParameters.Fathers()[i] != 0) ||
         (i > 0 && mParameters.Fathers()[i] >= i)
       )
    {
      std::cerr << "Fathers of first clone has to be 0 and" << std::endl;
      std::cerr << "Fathers of clones have to be introduced first!";
      std::cerr << std::endl;
      return 0;
    }
  }

  // Create the universe and reusable pointers to handled cells and types:
  Cell* pCell;
  CellType* pType;
  long double dt = 0.0L;
  int action;

  // Create all cell types and put them into a vector:
  std::vector <CellType*> pvCelltypes; // vector of all types
  for (int i = 0; i < mParameters.NumberOfClones(); i++) {

    // Generate rainbow colors:
    unsigned int hue_deg = 360 * i / mParameters.NumberOfClones();
    unsigned int color = HSLtoRGB(hue_deg, 100, 50);

    // New cell type:
    pType = new CellType(i+1,
                         mParameters.Birthrates()[i],
                         mParameters.Deathrates()[i],
                         mParameters.Aggressions()[i],
                         mParameters.PushPower()[i],
                         mParameters.Mutationrates()[i],
                         GetRValue(color),
                         GetGValue(color),
                         GetBValue(color));

    // Insert cell types:
    pvCelltypes.push_back(pType);
  }

  // Create and insert a cell of the first type:
  pCell = new Cell(pvCelltypes[0]);
  this->InsertCell(mParameters.Spaces()[0],
                   CENTER(mParameters.SizeX()[mParameters.Spaces()[0]]),
                   CENTER(mParameters.SizeY()[mParameters.Spaces()[0]]),
                   CENTER(mParameters.SizeZ()[mParameters.Spaces()[0]]),
                   pCell);

  pCell->MutateCellFixedNumber(mParameters.NumberClonalMutations());

  // universe time to start time of fist clone:
  this->IncrementTimeBy(mParameters.CloneStartTimes()[0]);
  int next_clone = 1; // count number of next type to introduce:

  // Debug messages:
  if (verbose) {
    std::cout << std::endl;
    std::cout << "########## Comment #############" << std::endl;
    std::cout << " New clone: " << next_clone - 1 << std::endl;
    std::cout << " Time: " << this->Time() << std::endl;
    std::cout << "################################" << std::endl;
    std::cout << std::endl;
  }

  bool kill_regrow_done=false;
  // Run till limit of universe has been reached:
  while(!this->LimitIsReached()) {

    // Print new universe to screen:
    DI(if (this->Time() - mLastPrintTime > mParameters.DisplayFrequency()) {)
    DI(mLastPrintTime = this->Time();)
    DI(this->Display();)
    DI(})

    // Check if new clones need to be introduced in each cycle:
    for (; next_clone < mParameters.NumberOfClones() &&
           mParameters.CloneStartTimes()[next_clone] <= this->Time();
           next_clone++)
    {

      // Debug messages:
      if (verbose) {
        std::cout << std::endl;
        std::cout << "########## Comment #############" << std::endl;
        std::cout << " New clone: " << next_clone << std::endl;
        std::cout << " Time: " << this->Time() << std::endl;
        std::cout << "################################" << std::endl;
        std::cout << std::endl;
      }

      // Introduce a new clone:
      int c_father = mParameters.Fathers()[next_clone];
      int n_space = mParameters.Spaces()[next_clone];
      pCell = pvCelltypes[c_father]->RandomMember();

      if (pCell->N() != n_space) {
        pCell = this->RemoveCell(pCell);
        pCell->Type(pvCelltypes[next_clone]);

        this->InsertCell(n_space,
                         CENTER(mParameters.SizeX()[n_space]),
                         CENTER(mParameters.SizeY()[n_space]),
                         CENTER(mParameters.SizeZ()[n_space]),
                         pCell);
      } else {
        pCell->Type(pvCelltypes[next_clone]);
      }
    } // end of for loop for introduction of new clones

    if (mParameters.KillRegrowTime()[0] > 0 && !kill_regrow_done &&
      mParameters.KillRegrowTime()[0] <= this->Time()) {

      kill_regrow_done = true;

      int curr_popsize = 0;
      for (int i = 0; i < mpTypes.size(); i++) {
        curr_popsize += mpTypes[i]->NumMembers();
      }

      int popsize_tokill = KILL_REGROW_PERCENT*curr_popsize;

      int popsize_killed = 0;
      while (popsize_killed < popsize_tokill) {

        int rand_celltype_i = rand() % static_cast<int>(mpTypes.size());

        if (mpTypes[rand_celltype_i]->NumMembers() > 1) {

          popsize_killed += 1;
          mpTypes[rand_celltype_i]->RandomMember()->Kill();

        }
      }

    } else {

      this->RecordToHistory();
      this->NextReactionType(&dt, &action)->RandomMember()->DoAction(&action);
      this->IncrementTimeBy(dt);

    }

  } // stop running after reaching limit
  next_clone = 0;
  // Return the memory location of the Universe array:
  return true;
}

bool Universe::RunSimulation() {
  return RunSimulation(false);
}


void Universe::SaveResults() const {
  std::string output_prefix = mParameters.OutputPrefix();
  int sim_id = mParameters.SimulationId();
  int n = 0;

  // Save simultation latice to file:
  this->SpaceToCsvFile(n, output_prefix + SPACE_OUPUT_FILE_PREFIX
                        + "_" + std::to_string(sim_id) + ".csv");

  // save total number of mutations per cell (from the latice) to file:
  this->SpaceGenerationToFile(n, output_prefix + GENOT_OUPUT_FILE_PREFIX
                                 + "_" + std::to_string(sim_id) + ".csv");

  // Save cell type properties to file:
  this->TypesToCsvFile(output_prefix + TYPES_OUPUT_FILE_PREFIX
                         + "_" + std::to_string(sim_id) + ".csv");

  // Save the whole phylogeny as nexus file:
  this->PhylogeniesToFile(output_prefix + PHYLO_OUPUT_FILE_PREFIX
                            + "_" + std::to_string(sim_id));

  // Last but not least, save the universe history to a file:
  this->HistoryToFile(output_prefix + HISTORY_OUPUT_FILE_PREFIX
                       + "_" + std::to_string(sim_id) + std::string(".csv"));

  // Save simulation parameters to file:
  mParameters.SaveToFile(output_prefix + PARAMS_OUTPUT_FILE_PREFIX
                          + "_" + std::to_string(sim_id) + std::string(".csv"));


  // Save display to image:
  DI(this->SaveDisplayToFile(output_prefix + DIPLAY_FILE
                              + "_" + std::to_string(sim_id) + ".bmp");)


}


// Optional features:
#ifdef _DISPLAY_

void Universe::Display() const {

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glClearColor(0.4f, 0.4f, 0.4f, 0.0f); // gray Background
  glEnable(GL_COLOR_MATERIAL);

  // Determine maximum number of generations:
  double MaxGen = 0.0;
  for (int n = 0; n < mvSpace.size(); n++) {
    for(int j = 0; j < mvSpace[n].shape()[0]; j++) {
      for(int k = 0; k < mvSpace[n].shape()[1]; k++) {
        for(int l = 0; l < mvSpace[n].shape()[2]; l++) {
          if (this->ContainsCell(n,j,k,l) &&
              MaxGen < mvSpace[n][j][k][l]->AssociatedNode()->Generation())
          {
            MaxGen = mvSpace[n][j][k][l]->AssociatedNode()->Generation();
          }
        }
      }
    }
  }

  for (int n = 0; n < mvSpace.size(); n++) {
    for(int j = 0; j < mvSpace[n].shape()[0]; j++) {
      for(int k = 0; k < mvSpace[n].shape()[1]; k++) {
        for(int l = 0; l < mvSpace[n].shape()[2]; l++) {
          if (this->ContainsCell(n,j,k,l)) {

            glMatrixMode(GL_MODELVIEW);
            glPushMatrix(); //remember current matrix
            glTranslatef(j + n * mMaxSpaceDim, k, l);

            unsigned char* cid = mvSpace[n][j][k][l]->Type()->Color();
            double cov = mvSpace[n][j][k][l]->AssociatedNode()->Generation()/MaxGen;
            glColor3d(cid[0]/255.0*cov, cid[1]/255.0*cov, cid[2]/255.0*cov);

            GLfloat mat_specular[4] = {0.2,0.2,0.2,1.0};
            glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
            GLfloat low_shininess[] = { 100.0 };
            glMaterialfv(GL_FRONT, GL_SHININESS, low_shininess);

            if (mvSpace[n].shape()[2] == 1) {
              glutSolidCube(1.0);
            } else {
              glutSolidSphere(2.0,10,10);
            }
            glPopMatrix(); //restore matrix
            // if (max_delta < j - static_cast <int>(mSpace.shape()[0]) / 2)
            //   max_delta = j - static_cast <int>(mSpace.shape()[0]) / 2;
          }
        }
      }
    }
  }

  glutSwapBuffers();
  glutPostRedisplay();
}

void Universe::WaitTillDisplayClosed() const {

  // Update display once.
  this->Display();

  // ToDo: Add detection of display state:
  // Ask for user input to continue.
  std::cout << "Press any key to continue program." << std::endl;
  std::cin.get();
}

void Universe::SaveDisplayToFile(std::string file) const {
  //mpBuffer->save(file.c_str());
}

#endif // _DISPLAY_
