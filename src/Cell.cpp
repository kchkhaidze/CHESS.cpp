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

#ifdef _ALLOWDIEOUT_
#define DIEOUT true
#else
#define DIEOUT false
#endif

#ifdef _DEBUG_
#define D(x) x
#else
#define D(x)
#endif

#include "extern_global_variables.hpp"
#include "Cell.hpp"
#include "CellType.hpp"
#include "Phylogeny.hpp"
#include "Universe.hpp"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <iostream>
#include <array>
#include <vector>


// Statics:
unsigned int Cell::msNextId = 1;


// Constructors:
Cell::Cell () // Default definition without given location and types.
  : mId(msNextId++),
    mpUniverse(0),
    mX(0),
    mY(0),
    mZ(0),
    mTypeIndex(0),
    mpType(0),
    mpNode(0),
    mTotalMutNum(0)
  {}

Cell::Cell (CellType* pType) // Mostly default but a given cell type.
  : mId(msNextId++),
    mpUniverse(0),
    mX(0),
    mY(0),
    mZ(0),
    mTypeIndex(0),
    mpType(pType),
    mpNode(0),
    mTotalMutNum(0)
  {
    mpType->RegisterMember(this);
  }


// Destructor:
Cell::~Cell(){ // Proper deletion of a cell.
  D(std::cout << "Cell::~Cell()" << std::endl;)
  D(std::cout << "Cell " << mId << " dies!" << std::endl;) // Debug message
  mpType->DeregisterMember(this);      // Let the associated type forget.
  mpUniverse->RemoveCell(mN,mX,mY,mZ);  // Remove the cell from the universe.
  mpNode->AssociatedCell(0);           // Unlink from the associated node.
  // Removal of cell complete.
}


// Getter functions:
unsigned int Cell::Id() const {return(mId);}
int Cell::TypeIndex() const {return mTypeIndex;};
int Cell::N() const {return mN;};
int Cell::X() const {return mX;};
int Cell::Y() const {return mY;};
int Cell::Z() const {return mZ;};
int Cell::TotalMutNum() const {return mTotalMutNum;};

bool Cell::TriesPush() const { // Random realization to push with prob aggression.
  double aggression = mpType->Aggression();
  boost::random::bernoulli_distribution<> bern_aggression(aggression);
  bool res = bern_aggression(rng);
  return res;
}

CellType* Cell::Type() {return mpType;};
PhylogenyNode* Cell::AssociatedNode() {return mpNode;}; // Assoc. node



// Setter functions:
void Cell::Location(Universe *pUniverse, int N, int X, int Y, int Z) {
  mpUniverse = pUniverse; mN = N, mX = X; mY = Y; mZ = Z;
}

void Cell::Location(int N, int X, int Y, int Z) { mN=N; mX=X; mY=Y; mZ=Z; }

void Cell::Location(int X, int Y, int Z) { mX=X; mY=Y; mZ=Z; }

void Cell::Type(CellType* pNewType) {
  // Change memberships:
  if (mpUniverse != 0) { // If the cell is member of a universe,
    mpUniverse->RegisterType(pNewType); // insert new cell type into universe.
  }

  // Change cell type memberships:
  mpType->DeregisterMember(this);
  pNewType->RegisterMember(this);

  // Update the TypeId of the associated node:
  if (mpNode != 0) { // If the cell has a associated node, member of a universe,
    mpNode->TypeId(pNewType->Id()); // update its TypeId.
  }
  // Update the TypeId of all ancestral nodes to 0 (undertermined/mix type):
  PhylogenyNode *pCurrentNode = mpNode;
  while ((pCurrentNode = pCurrentNode->UpNode()) != 0) {
    pCurrentNode->TypeId(0);
  }


  // Last but not least, update the type of the cell:
  mpType = pNewType;
};

void Cell::AssociatedNode(PhylogenyNode* new_node) { mpNode = new_node; };

void Cell::TypeIndex(int newIndex) {mTypeIndex = newIndex;};

void Cell::TotalMutNum(int newTotalMutNum) {mTotalMutNum = newTotalMutNum;};

// Other functions:
void Cell::MutateCell() {
  Cell::MutateCell(mpType->Mu());
}

void Cell::MutateCell(double mu){
  // Sample number of new muts from poisson:
  boost::random::poisson_distribution<int> dist_mutations(mu);
  int new_muts = dist_mutations(rng);

  // Append muts to node:
  mpNode->AddNewMutations(new_muts);
  // increment total mutation count for the cell
  this->TotalMutNum((this->TotalMutNum()) + new_muts);
}

void Cell::MutateCellFixedNumber(int new_muts){
  // Append muts to node:
  mpNode->AddNewMutations(new_muts);
}

void Cell::Print() const {
  std::cout << std::endl;
  std::cout << "############ Cell ##############" << std::endl;
  std::cout << "   ID: " << mId << std::endl;
  std::cout << "   Location:" << std::endl;
  std::cout << "       Universe: " << mpUniverse << std::endl;
  std::cout << "       x: " << mX << std::endl;
  std::cout << "       y: " << mY << std::endl;
  std::cout << "       z: " << mZ << std::endl;
  std::cout << "   Type: " << mpType->Id() << std::endl;
  std::cout << "       Birth rate: " << mpType->Birthrate() << std::endl;
  std::cout << "###############################" << std::endl;
}

void Cell::Divide() {

  // Store current associated node, then branch(turn in) to the left and mutate:
  PhylogenyNode* old_node = this->mpNode;
  PhylogenyNode* new_left_node = new PhylogenyNode(this, old_node);
  old_node->LeftNode(new_left_node);
  this->MutateCell();

  // Get all free spaces around current pos:
  int new_x, new_y, new_z;
  std::vector< std::array<int, 3> > free_neighbours;
  free_neighbours = mpUniverse->FreeNeighbours(mN, mX, mY, mZ);

  // Determine the number of free neighbours:
  std::vector< std::vector<int> >::size_type size = free_neighbours.size();

  if (size > 0){ // If there is free space, choose a random free position
    boost::random::uniform_int_distribution<int> unif_neighbour(0, size - 1);
    int sel_neighbour = unif_neighbour(rng);
    std::array<int, 3> sel_free_neighbour = free_neighbours[sel_neighbour];

    new_x = sel_free_neighbour[0];
    new_y = sel_free_neighbour[1];
    new_z = sel_free_neighbour[2];

    // Test if the new cell touches border:
    if (!mpUniverse->ContainsCoordinate(mN,
                                        new_x + (new_x - mX),
                                        new_y + (new_y - mY),
                                        new_z + (new_z - mZ)))
    {
      mpUniverse->MarkLimitIsReached();
    }

  } else if (this->TriesPush() && this->RandomPush(new_x, new_y, new_z)){
    // or if there is no free space, try a push
    // Note: New random location are stored in params passed to "RandomPush"
    ;
  } else { // If the push failed as well:
    return; // No devision takes place;
  }
  //Create second daughter, insert on branch to the right of old node and mutate
  Cell* pDaughter = new Cell(mpType);
  PhylogenyNode* new_right_node = new PhylogenyNode(pDaughter, old_node);
  old_node->RightNode(new_right_node);
  mpUniverse->InsertCell(mN, new_x, new_y, new_z, pDaughter, false);
  pDaughter->MutateCell();

}

void Cell::Kill() {
  
  delete this;
}

void Cell::DoAction(int *action){

  if (mpUniverse == 0) {return;}

  switch(*action){
    case 1:  //cell dies
      if (DIEOUT || mpType->NumMembers() > 1) {
        delete this;
      }
      break;
    case 2:  //cell divides
      this->Divide();
      break;
    default: //undefined action requested
      std::cerr << "Error in void Cell::DoAction(int *):\n";
      std::cerr << "   > Action " <<  *action << " undefined." << std::endl;
      exit(EXIT_FAILURE);
  }
}

// Other functions:
bool Cell::RandomPush(int& newx, int& newy, int& newz) {

   // Cells that are not in a universe can't push:
   if (mpUniverse == 0) {
     newx = mX;
     newy = mY;
     newz= mZ;
     return false;
   }

  // Init. direction vectors:
  // Create a random vector:
  boost::random::uniform_real_distribution<double> rand_phi(-M_PI, M_PI);
  boost::random::uniform_real_distribution<double> rand_u(-1.0, 1.0);

  double phi = rand_phi(rng);
  double u = rand_u(rng);

  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;

  if (mpUniverse->ContainsCoordinate(mN,1,0,0))
    dx = std::sqrt(1.0 - std::pow(u, 2.0)) * std::cos(phi);
  if (mpUniverse->ContainsCoordinate(mN,0,1,0))
    dy = std::sqrt(1.0 - std::pow(u, 2.0)) * std::sin(phi);
  if (mpUniverse->ContainsCoordinate(mN,0,0,1))
    dz = u;


  // Normalize so that largest delta value is 1.0:
  double max_component = fmax(std::abs(dx), fmax(std::abs(dy), std::abs(dz)));
  dx /= max_component;
  dy /= max_component;
  dz /= max_component;

  // Debug messages:
  D(std::cout << " Trying random push: " << std::endl;)
  D(std::cout << "    id: " << mId << std::endl;)
  D(std::cout << "    x: " << mX << " delta: " << dx << std::endl;)
  D(std::cout << "    y: " << mY << " delta: " << dy << std::endl;)
  D(std::cout << "    z: " << mZ << " delta: " << dz << std::endl;)


  // Then try to push the cell towards this position.
  // If this position is taken than the next cell(s) will move as well,
  // if none hits a border.
  double newX = mX;
  double newY = mY;
  double newZ = mZ;
  while (round(newX) == mX && round(newy) == mY && round(newZ) == mZ) {
    newX += dx;
    newY += dx;
    newZ += dx;
  }

  newx = round(mX);
  newy = round(mY);
  newz = round(mZ);
  return mpUniverse->PushCell(mN, mX, mY, mZ,
                              dx, dy, dz,
                              mpType->PushPower());
}
