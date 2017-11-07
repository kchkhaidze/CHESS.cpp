//#define DEBUG
#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

#include "Cell.hpp"
#include "CellType.hpp"
#include "Phylogeny.hpp"
#include "Universe.hpp"
#include <random>
#include <iostream>
#include <array>
#include <vector>


// Statics:
unsigned int Cell::msNextId = 1;


// Constructors:
Cell::Cell () // Default definition without location and types.
  : mId(msNextId++),
    mpUniverse(0),
    mX(0),
    mY(0),
    mZ(0),
    mpType(0),
    mpNode(0)
  {}

Cell::Cell (CellType* pType) // Defaults except of cell type.
  : mId(msNextId++),
    mpUniverse(0),
    mX(0),
    mY(0),
    mZ(0),
    mpType(pType),
    mpNode(0)
  {
    mpType->RegisterMember(this);
  }


// Destructor:
Cell::~Cell(){ // Proper deletion of a cell.
  D(std::cout << "Cell " << mId << " dies!" << std::cout;)
  mpType->DeregisterMember(this);
  mpUniverse->RemoveCell(mX, mY, mZ);
  mpNode->AssociatedCell(0);
}


// Getter functions:
unsigned int Cell::Id() {return(mId);}
CellType* Cell::Type() {return mpType;};

bool Cell::AsProgenitorDies() {
  return ((double) rand() / (RAND_MAX)) <= mpType->Alpha();
}

bool Cell::AsOffspringDies() {
  return ((double) rand() / (RAND_MAX)) <= (1.0 - mpType->Beta());
}

bool Cell::TriesPush() {
  return ((double) rand() / (RAND_MAX)) <= mpType->Aggression();
}

PhylogenyNode* Cell::AssociatedNode() {return mpNode;};


// Setter functions:
void Cell::Location(Universe *pUniverse, int X, int Y, int Z) {
  mpUniverse = pUniverse; mX = X; mY = Y; mZ = Z;
}

void Cell::Location(int X, int Y, int Z) { mX = X; mY = Y; mZ = Z; }

void Cell::Type(CellType* pNewType) {
  mpUniverse->RegisterType(pNewType);
  mpType->DeregisterMember(this);
  pNewType->RegisterMember(this);
  mpType = pNewType;
};

void Cell::AssociatedNode(PhylogenyNode* new_node) { mpNode = new_node; };


// Other functions:
void Cell::MutateCell() {
  Cell::MutateCell(mpType->Mu());
}

void Cell::MutateCell(float mu){
  //:ToDo: Reimplement mutations of cells!
  // Sample number of new muts from poisson:
  std::random_device rd;
  std::mt19937 gen(rd());
  std::poisson_distribution<> d(mu);
  int new_muts = d(gen);

  // Append muts to node:
  mpNode->AddNewMutations(new_muts);
}


void Cell::Print() {
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

  // Cell without a universe can't divide!
  if (mpUniverse == 0) {
    return;
  }

  // Debug messages:
  D(std::cout << std::endl;)
  D(std::cout << "########## Division ###########" << std::endl;)
  D(std::cout << "   ID: " << mId << std::endl;)
  D(std::cout << "   Location:" << std::endl;)
  D(std::cout << "       Universe: " << mpUniverse << std::endl;)
  D(std::cout << "       x: " << mX << std::endl;)
  D(std::cout << "       y: " << mY << std::endl;)
  D(std::cout << "       z: " << mZ << std::endl;)
  D(std::cout << std::endl;)

  if (this->AsProgenitorDies()) {
    D(std::cout << " Progenitor died" << std::endl;)
    D(std::cout << "###############################" << std::endl;)
    delete this;
  } else {
    // Store current associated node, then branch to the left and mutate:
    PhylogenyNode* old_node = this->mpNode;
    PhylogenyNode* new_left_node = new PhylogenyNode(this, old_node);
    old_node->LeftNode(new_left_node);
    this->MutateCell();

    // Get all free spaces around current pos:
    int new_x, new_y, new_z;
    std::vector< std::array<int, 3> > free_neighbours;
    free_neighbours = mpUniverse->FreeNeighbours(mX, mY, mZ);

    std::vector< std::vector<int> >::size_type size = free_neighbours.size();
    if (size > 0){
      int r_unif = rand() % size;
      std::array<int, 3> sel_free_neighbour = free_neighbours[r_unif];

      new_x = sel_free_neighbour[0];
      new_y = sel_free_neighbour[1];
      new_z = sel_free_neighbour[2];
    } else if (this->TriesPush() && this->RandomPush(new_x, new_y, new_z)){
      ;
    } else {
      D(std::cout << "  No free locations!" << std::endl;)
      D(std::cout << "  Cell could not push." << std::endl;)
      D(std::cout << "###############################" << std::endl;)
      // No devision takes place; RETURN! //////////////////////////////////////
      return;
    }

    // Debug messages:
    D(std::cout << std::endl;)
    D(std::cout << "  Choosen location:" << std::endl;)
    D(std::cout << "    x: " << new_x;)
    D(std::cout << "    y: " << new_y;)
    D(std::cout << "    z: " << new_z << std::endl;)

    Cell* pDaughter = new Cell(mpType);
    PhylogenyNode* new_right_node = new PhylogenyNode(pDaughter, old_node);
    old_node->RightNode(new_right_node);
    mpUniverse->InsertCell(new_x, new_y, new_z, pDaughter, false);
    pDaughter->MutateCell();


    if (pDaughter->AsOffspringDies()) {
      D(std::cout << " Daughter died" << std::endl;)
      D(std::cout << "###############################" << std::endl;)
      delete pDaughter;
    } else {
      D(std::cout << "###############################" << std::endl;)
    }
  }

}

// Other functions:
bool Cell::RandomPush(int& newx, int& newy, int& newz) {
   // Cell that are not in a universe can't push:
   if (mpUniverse == 0) {
     newx = mX;
     newy = mY;
     newz= mZ;
     return false;
   }

  int dx = 0;
  int dy = 0;
  int dz = 0;

  // :ToDo: replace this with a better solution
  while ((dx == 0 && dy == 0 && dz == 0) ||
          !mpUniverse->ContainsCoordinate(mX + dx, mY + dy, mZ + dz)) {
    dx = rand() % 3 - 1;
    dy = rand() % 3 - 1;
    dz = rand() % 3 - 1;
  }


  D(std::cout << " Trying random push: " << std::endl;)
  D(std::cout << "    id: " << mId << std::endl;)
  D(std::cout << "    x: " << mX << "->" << mX+dx << std::endl;)
  D(std::cout << "    y: " << mY << "->" << mY+dy << std::endl;)
  D(std::cout << "    z: " << mZ << "->" << mZ+dz << std::endl;)


  newx = mX + dx;
  newy = mY + dy;
  newz= mZ + dz;
  return mpUniverse->PushCell(mX + dx, mY + dy, mZ + dz,
                              dx, dy, dz,
                              mpType->PushPower());

}
