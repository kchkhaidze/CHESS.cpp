//#define DEBUG
#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

#include "CellType.hpp"
#include "Shapes.hpp"
#include "Cell.hpp"
#include "Universe.hpp"
#include "Phylogeny.hpp"
#include <random>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <fstream>

// Universe ////////////////////////////////////////////////////////////////////


// Constructor:
Universe::Universe(int X, int Y, int Z)
  : mTime(0.0), mSpace(boost::extents[X][Y][Z]), mLimitReached(false)
  {
    // Make the space empty:
    for(size_t i = 0; i < mSpace.shape()[0]; i++) {
      for(size_t j = 0; j < mSpace.shape()[1]; j++) {
        for(size_t k = 0; k < mSpace.shape()[2]; k++) {
          mSpace[i][j][k] = 0;
        }
      }
    }
  }

// Destructor:
//:ToDo: Add destructor


// Getter functions:
float Universe::Time() { return mTime; }
bool  Universe::LimitIsReached() {return mLimitReached;};

bool  Universe::ContainsCell(int x, int y, int z) {
  return mSpace[x][y][z] != 0;
}

bool Universe::ContainsCoordinate(int x, int y, int z) {
  return x < mSpace.shape()[0] &&
         y < mSpace.shape()[1] &&
         z < mSpace.shape()[2] &&
         x >= 0 &&
         y >= 0 &&
         z >= 0;
}

Cell *Universe::GetCell(int x, int y, int z) { return mSpace[x][y][z]; }

Cell* Universe::RemoveCell (int x, int y, int z) {
  Cell* p_removed_cell = mSpace[x][y][z];
  mSpace[x][y][z] = 0;
  p_removed_cell->Location(0, 0, 0, 0);
  return p_removed_cell;
}

CellType* Universe::NextReactionType(double& r_delta_time){

  // Return null if the universe contains no types:
  if (mpTypes.size() == 0) {
    r_delta_time = 0;
    return 0;
  }

  // Debug messages:
  D(std::cout << std::endl;)
  D(std::cout << "########## Sampling ###########" << std::endl;)
  D(std::cout << "  Next reaction type:" << std::endl;)
  D(std::cout << std::endl;)


  // Keep track of minum index and dt:
  double delta_time_minimum;
  int index_minimum;

  // Calculate dt for all types in the universe:
  for(std::vector<CellType*>::size_type i = 0; i < mpTypes.size(); i++) {
    // Debug messages:
    D(std::cout << "    Candidate " << i <<  ":" << std::endl;)
    D(std::cout << "      No: " << i << std::endl;)
    D(std::cout << "      Type: " << mpTypes[i]->Id() << std::endl;)
    D(std::cout << "      BR: " << mpTypes[i]->Birthrate() << std::endl;)
    D(std::cout << "      N: " << mpTypes[i]->NumMembers() << std::endl;)

    // :ToDo: Add better sampling method.
    // :ToDo: Check this is correct.
    double runi = ((double) rand() / (RAND_MAX));
    double c_birthrate = mpTypes[i]->Birthrate();
    double c_number = mpTypes[i]->NumMembers();
    double delta_time_current = (log(1) - log(runi)) / (c_number * c_birthrate);

    // Debug messages:
    D(std::cout << "      X: " << runi << std::endl;)
    D(std::cout << "      dt: " <<  delta_time_current <<  std::endl;)
    D(std::cout << std::endl;)

    // Keep track of minum index and dt:
    if (i == 0 || delta_time_minimum > delta_time_current) {
      index_minimum= i;
      delta_time_minimum = delta_time_current;
    }

  }

  // Debug messages:
  D(std::cout << "  Type No: #" << index_minimum << std::endl;)
  D(std::cout << "  Type: " << mpTypes[index_minimum]->Id() << std::endl;)
  D(std::cout << "  dt: " << delta_time_minimum << std::endl;)
  D(std::cout << "###############################" << std::endl;)

  // Return types pointer:
  r_delta_time = delta_time_minimum;
  return(mpTypes[index_minimum]);
}


// :ToDo: Replace this with new, better version ...
std::vector< std::array<int, 3> > Universe::FreeNeighbours(int x, int y, int z){
  std::vector< std::array<int, 3> > result;
  std::array<int, 3> cords;

  int xstart, xend;
  if (mSpace.shape()[0] > 1) {
    xstart = x - 1;
    xend = x + 1;
  } else {
    xstart = x;
    xend = x;
  }

  int ystart, yend;
  if (mSpace.shape()[1] > 1) {
    ystart = y - 1;
    yend = y + 1;
  } else {
    ystart = y;
    yend = y;
  }

  int zstart, zend;
  if (mSpace.shape()[2] > 1) {
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
        if (this->ContainsCoordinate(xn, yn, zn) &&
              !this->ContainsCell(xn, yn, zn) &&
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


bool Universe::TakeSample(Shape* shape,
                          std::vector <std::string>& rvOutputId,
                          std::vector <float>& rvOutputVaf){

  // Part 1 - Collection of samples. ////////////////////////////

  // Debug messages:
  D(std::cout << std::endl;)
  D(std::cout << "###### Cell sampling ##########" << std::endl;)
  D(std::cout << "  Collecting cells:" << std::endl;)

  // Vector to store all current ancestries of sampled cells:
  std::vector< std::vector <PhylogenyNode*> > vpAncestries;

  // Counter for the number of sampled nodes:
  int normal_samples = 0;
  int tumour_samples = 0;

  // Ints to store the next cell to pick:
  int x, y, z;

  while (shape->next_coordinate(x, y, z)) { // Get next location from shape.
    if (mSpace[x][y][z] == 0) { // Check if location is empty
      normal_samples++;

      // Debug messages:
      D(std::cout << "  Normal cell " << normal_samples << ": " << std::endl;)
      D(std::cout << "    x: " << x << " y: " << y << " z:" << z << std::endl;)
    } else { // selected location is not empty
      // Push current ancestry into the vector
      tumour_samples++;
      vpAncestries.push_back(mSpace[x][y][z]->AssociatedNode()->NodeAncestry());

      // Debug messages:
      D(std::cout << "  Tumour cell " << tumour_samples << ": " << std::endl;)
      D(std::cout << "    x: " << x << " y: " << y << " z:" << z << std::endl;)
    }
  }
  // Debug messages:
  D(std::cout << "###############################" << std::endl;)

  // Part 2 - Interpretation of samples. ////////////////////////

  // Debug messages:
  D(std::cout << std::endl;)
  D(std::cout << "#### Sample interpretation ####" << std::endl;)
  D(std::cout << " Number of cells: " << std::endl;)
  D(std::cout << "   Tumour: " << tumour_samples << std::endl;)
  D(std::cout << "   Normal: " << normal_samples << std::endl;)
  D(std::cout << std::endl;)

  // By tracking back along the phylogenies from the end to the fron we will be
  // able to calculate the VAF of all mutations at each sample generation.

  // For each level store the mutation related data in the following arrays:
  unsigned int node_ids[tumour_samples + 1];
  unsigned int mutation_number[tumour_samples + 1];
  unsigned int affected_cells[tumour_samples + 1];


  bool any_new = true; // During each iteration use this mark new muts present
  while (any_new) {
    any_new = false;

    // Reset all arrays for the current level of mutation related data:
    for (int i = 0; i < tumour_samples + 1; i++) {
      mutation_number[i] = 0;
      affected_cells[i] = 0;
      node_ids[i] = 0;
    }

    // Iterate over all ancestry vectors and collect the data:
    for (std::vector< std::vector<PhylogenyNode*> >::size_type i = 0;
         i < vpAncestries.size();
         i++)
    {

      // Check that there are data in the current cells anyestry:
      if (vpAncestries[i].size() > 0) {
        // Mark that new data were present at the last level:
        any_new = true;

        // Index for the end of the currents cells ancestry:
        std::vector<PhylogenyNode*>::size_type j = vpAncestries[i].size() - 1;

        // Extract id of the last node of current ancestry:
        unsigned int cid = vpAncestries[i][j]->Id();

        // Insert data into the appropiate array elements:
        for (int k = 0; k <= tumour_samples; k++) { // Search through the array
          if (node_ids[k] == cid) { // Found at position k of arrays.
            affected_cells[k]++;
            break;
          } else if (node_ids[k] == 0){ // Not present in array yet. Insert at end.
            affected_cells[k]++;
            node_ids[k] = cid;
            mutation_number[k] = vpAncestries[i][j]->NumMutations();
            break;
          }
        }

        // Remove last element of vector:
        vpAncestries[i].pop_back();
      }
    }

    // Interpret all node sets at the current generation level:
    for (int k = 0; k <= tumour_samples && node_ids[k] != 0; k++) {
      //std::cout << "k: " << k << ", node_ids: " << node_ids[k]
      //  << ", affected_cells: " <<  affected_cells[k]
      //  << ", mutation_number: "<< mutation_number[k] << std::endl;

      // Convert node id to a hex
      std::stringstream out;
      out << std::hex << node_ids[k];
      std::string hex_node_id = out.str();

      // Binomial sampling related:
      // std::default_random_engine generator;
      // int depth = 100;
      float expected_vaf = 0.5 * affected_cells[k] / tumour_samples;
      // std::binomial_distribution<int> distribution(depth, expect_freq);

      // Sample each mutation in current node once:
      for (int l = 0; l < mutation_number[k]; l++) {
        // int alt = distribution(generator);
        // int ref = depth - alt;
        // float sampled_freq = alt * 1.0 / depth;

        // Convert mutation id to a hex string:
        out.str("");
        out << std::hex << l;
        std::string hex_mutation_id = out.str();

        // Append mutation information to output vectors:
        //std::cout << hex_node_id << "X" << hex_mutation_id << std::endl;
        rvOutputId.push_back(hex_node_id + "X" + hex_mutation_id);
        rvOutputVaf.push_back(expected_vaf);
      }
    }
  }
  return true;
}

// Setter functions:
void Universe::IncrementTimeBy(float Delta){ mTime += Delta; }
void Universe::MarkLimitIsReached() { mLimitReached = true; }

bool Universe::InsertCell(int x, int y, int z, Cell* pCell) {
  return InsertCell(x, y, z, pCell, true);
}

bool Universe::InsertCell(int x, int y, int z, Cell* pCell,
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
  if (this->ContainsCoordinate(x,y,z) == false) {
    std::cout << "Error: Insertion of cell failed:" << std::endl;
    std::cout << "Error: Coordinates (i+1) were: " << x + 1 << ",";
    std::cout << y + 1 << "," << z + 1 << std::endl;
    std::cout << "Error:  Universe size was: " << mSpace.shape()[0] << ",";
    std::cout << mSpace.shape()[1] << "," << mSpace.shape()[2] << std::endl;
    //exit EXIT_FAILURE;
    return false;
  }

  // Dont allow inserts on top of cell:
  if (this->ContainsCell(x,y,z)) {
    std::cout << "Error: Insertion of cell failed:" << std::endl;
    std::cout << "Error: Coordinates (i) were: " << x << ",";
    std::cout << y << "," << z  << std::endl;
    std::cout << "Error:  Position was taken by" << mSpace[x][y][z] << std::endl;
    exit(EXIT_FAILURE);
    //return false;
  }

  // Insert type into universe
  if (is_new_lineage) {
    this->RegisterType(pCell->Type());
  }

  // Update position and universe of the cell:
  pCell->Location(this, x, y, z);

  // Update universe:
  mSpace[x][y][z] = pCell;

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
    mpTypes.push_back(p_new_type);
  }
}

// Output Functions:
void Universe::Print() {
  std::cout << std::endl;
  std::cout << "############ Space ############" << std::endl;
  std::cout << std::endl;

  for(int k = mSpace.shape()[2] - 1; k >= 0; k--) {
    std::cout << "  z " << k << ":" << std::endl;
    for(int j = mSpace.shape()[1] - 1; j >= 0; j--) {
      for(int i = 0; i < mSpace.shape()[0]; i++) {
        std::cout << " " <<  mSpace[i][j][k]->Type()->Id() << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << "###############################" << std::endl;
}

void Universe::Display(cimg_library::CImg<unsigned char>* pBuffer,
                      cimg_library::CImgDisplay* pDisp) {

  if (!pDisp->is_closed()) {

    // Make buffer empty
    pBuffer->fill(0);

    for(int j = 0; j < mSpace.shape()[0]; j++) {
      for(int k = 0; k < mSpace.shape()[1]; k++) {
        for(int l = 0; l < mSpace.shape()[2]; l++) {

          if (this->ContainsCell(j,k,l)) {
            unsigned char* cid = mSpace[j][k][l]->Type()->Color();
            *pBuffer->data(j,k,l,0) = cid[0];
            *pBuffer->data(j,k,l,1) = cid[1];
            *pBuffer->data(j,k,l,2) = cid[2];
          }
        }
      }
    }
    pBuffer->display(*pDisp);
  }
}

void Universe::SpaceToCsvFile(std::string out_prefix) {
  std::string suffix = ".csv";
  std::ofstream output_stream;
  std::string outfile_name;

  for(int z = 0; z < mSpace.shape()[2]; z++) {
    outfile_name = out_prefix + std::string("_z_") + std::to_string(z) + suffix;
    output_stream.open(outfile_name, std::ios::out | std::ios::trunc);

    if (output_stream.is_open()) {
      for(int y = mSpace.shape()[1] - 1; y >= 0; y--) {
        for(int x = 0; x < mSpace.shape()[0]; x++) {
          if (this->ContainsCell(x,y,z)) {
            output_stream << std::to_string(mSpace[x][y][z]->Type()->Id());
          } else {
            output_stream <<  "NA";

          }
          if (x < mSpace.shape()[0] - 1) {
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

void Universe::TypesToCsvFile(std::string out_prefix) {
  std::string suffix = ".csv";
  std::string outfile_name = out_prefix + suffix;
  std::ofstream output_stream;
  output_stream.open(outfile_name, std::ios::out | std::ios::trunc);

  if (output_stream.is_open()) {
    output_stream << "ID,birthrate,alpha,beta,aggression,mu,pushpower,N"
      << std::endl;

    for(std::vector<CellType *>::size_type i = 0; i < mpTypes.size(); i++) {
      output_stream << std::to_string(mpTypes[i]->Id()) << ","
        << std::to_string(mpTypes[i]->Birthrate()) << ","
        << std::to_string(mpTypes[i]->Alpha()) << ","
        << std::to_string(mpTypes[i]->Beta()) << ","
        << std::to_string(mpTypes[i]->Aggression()) << ","
        << std::to_string(mpTypes[i]->Mu()) << ","
        << std::to_string(mpTypes[i]->PushPower()) << ","
        << std::to_string(mpTypes[i]->NumMembers()) << std::endl;
    }
    output_stream.close();
  } else {
    std::cerr << "Couldn't write types to csv file." << std::endl;
  }
}

void Universe::PutNodesToStream(PhylogenyNode* current_node,
                                std::ofstream& outstream) {
  // Print this node:
  //PhylogenyNode* up = current_node->UpNode();
  PhylogenyNode* left = current_node->LeftNode();
  PhylogenyNode* right = current_node->RightNode();
  //Cell* cell = current_node->AssociatedCell();
  //unsigned int up_id = (up == 0) ? 0 : up->Id();
  //unsigned int left_id = (left == 0) ? 0 : left->Id();
  //unsigned int right_id = (right == 0) ? 0 : right->Id();
  //unsigned int cell_id = (cell == 0) ? 0 : cell->Id();

  //  outstream << current_node->Id() <<  up_id << "," << left_id << ","
  //  << right_id <<  "," << current_node->Generation() << ","
  //  << current_node->NumMutations() << "," << cell_id << "," << std::endl;

  std::string spacer =  std::string(current_node->Generation(), ' ');
  if (left == 0 && right == 0) { // "Print as Leaf"
    outstream << spacer << current_node->Id() << ":" << current_node->NumMutations();
  }  else { // Print as node?
    outstream << spacer << "(" << std::endl;
    if (left != 0) {
      PutNodesToStream(left, outstream);
      outstream << "," << std::endl;
    }
    if (right != 0){
      PutNodesToStream(right, outstream);
      outstream << std::endl;
    }
    outstream << spacer << "):" << current_node->NumMutations();
  }
}

void Universe::PhylogeniesToFile(std::string out_prefix) {
  std::string suffix = ".tree";
  std::ofstream output_stream;
  std::string outfile_name;

  for(std::vector<PhylogenyRoot*>::size_type i=0; i<mpPhylogenies.size(); i++) {


      outfile_name = out_prefix + "_phylogeny_" + std::to_string(i) + suffix;
      output_stream.open(outfile_name, std::ios::out | std::ios::trunc);

      if (output_stream.is_open()) {
        //output_stream << "id,up,left,right,gen,num_muts,cell" << std::endl;
        PhylogenyNode* current_node = mpPhylogenies[i]->Root();
        PutNodesToStream(current_node, output_stream);
        output_stream << std::endl;
      } else {
        std::cerr << "Couldn't write space to csv file." << std::endl;
      }
    output_stream.close();
  }
}


// Other functions:
bool Universe::PushCell(int x, int y, int z, int dx, int dy, int dz,
                        unsigned int push) {
  int new_x = x + dx;
  int new_y = y + dy;
  int new_z = z + dz;

  // If you hit a border stop pushes.
  if (this->ContainsCoordinate(new_x, new_y, new_z) == false) {
    D(std::cout << "  Hit border during push." << std::endl;)
    this->MarkLimitIsReached();
    return false;
  }

  // Test if next space is free.
  if (this->ContainsCell(new_x, new_y, new_z) == false ||
      (push > 0 && PushCell(new_x, new_y, new_z, dx, dy, dz, push - 1)))
  {
    // Debug messages:
    D(std::cout << "  Push sucessfull:" << std::endl;)
    D(std::cout << "    id: " << mSpace[x][y][z] << std::endl;)
    D(std::cout << "    x: " << x << "->" << new_x << std::endl;)
    D(std::cout << "    y: " << y << "->" << new_y << std::endl;)
    D(std::cout << "    z: " << z << "->" << new_z << std::endl;)

    // Move cell:
    mSpace[new_x][new_y][new_z] = mSpace[x][y][z];
    mSpace[x][y][z]->Location(new_x, new_y, new_z);
    mSpace[x][y][z] = 0;
    return true;
  }

  return false;
}
