//#define DEBUG
#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

#include "Phylogeny.hpp"
#include "Cell.hpp"
#include <iostream>

// PhylogenyRoot ///////////////////////////////////////////////////////////////

// Statics:
unsigned int PhylogenyRoot::msNextId = 0;

// Constructor:
PhylogenyRoot::PhylogenyRoot(Cell* pCell) : mId(msNextId++) {
  mpRoot = new PhylogenyNode(pCell);
}

PhylogenyNode* PhylogenyRoot::Root(){ return mpRoot; }

// PhylogenyNode ///////////////////////////////////////////////////////////////

// Statics:
unsigned int PhylogenyNode::msNextId = 0;

// Constructors
PhylogenyNode::PhylogenyNode(Cell* pCell) //
  : mId(msNextId++),
    mpCell(pCell),
    mpUp(0),
    mpLeft(0),
    mpRight(0),
    mGeneration(1),
    mNumMutsGeneration(0)
  {
    if (pCell->AssociatedNode() != 0) {
      pCell->AssociatedNode()->AssociatedCell(0);
    }
    pCell->AssociatedNode(this);
  }

PhylogenyNode::PhylogenyNode(Cell* pCell, PhylogenyNode* pUp)
  : mId(msNextId++),
  mpCell(pCell),
  mpUp(pUp),
  mpLeft(0),
  mpRight(0),
  mNumMutsGeneration(0)
  {
    mGeneration = pUp->Generation() + 1;
    if (pCell->AssociatedNode() != 0) {
      pCell->AssociatedNode()->AssociatedCell(0);
    }
    pCell->AssociatedNode(this);
  }

// Getters:
unsigned int PhylogenyNode::Id() { return mId; };
unsigned int PhylogenyNode::Generation(){ return mGeneration; }
unsigned int PhylogenyNode::NumMutations(){ return mNumMutsGeneration; }
Cell* PhylogenyNode::AssociatedCell() { return mpCell; };
PhylogenyNode* PhylogenyNode::UpNode() { return mpUp; };
PhylogenyNode* PhylogenyNode::LeftNode() { return mpLeft; };
PhylogenyNode* PhylogenyNode::RightNode() { return mpRight; };

// Setters:
void PhylogenyNode::AddNewMutations(int num_new_mutations) {
  mNumMutsGeneration += num_new_mutations;
}

void PhylogenyNode::AssociatedCell(Cell* pCell) { mpCell = pCell; }

void PhylogenyNode::UpNode(PhylogenyNode* pNode) {mpUp = pNode;};
void PhylogenyNode::LeftNode(PhylogenyNode* pNode) {mpLeft = pNode;};
void PhylogenyNode::RightNode(PhylogenyNode* pNode) {mpRight = pNode;};



// Output:
void PhylogenyNode::Print() {
  unsigned int up = (mpUp == 0) ? 0 : mpUp->Id();
  unsigned int left = (mpLeft == 0) ? 0 : mpLeft->Id();
  unsigned int right = (mpRight == 0) ? 0 : mpRight->Id();

  std::cout << "###### Phylogenetic Node ######" << std::endl;
  std::cout << "   ID: " << mId << std::endl;
  std::cout << "   Up: " << up << std::endl;
  std::cout << "   Left: " << left << std::endl;
  std::cout << "   Right: " << right << std::endl;
  std::cout << "   Generation: " << mGeneration << std::endl;
  std::cout << "   Mutations: " << mNumMutsGeneration << std::endl;
  if (mpCell != 0) { std::cout << "   Cell: " << mpCell << std::endl; }
  std::cout << "###############################" << std::endl;
}

void PhylogenyNode::PrintAncestry() {

  if (mpLeft == 0 && mpRight == 0) { // If this is a leaf
    std::string cell_id = (mpCell == 0) ? "NA" : std::to_string(mpCell->Id());
    std::cout << std::endl;
    std::cout << "###### Ancestry of cell #######" << std::endl;
    std::cout << "  Cell ID: " << cell_id << std::endl;
    std::cout << std::endl;
  }

  this->Print();

  if (mpUp != 0) { // If next is not root
    std::cout << "              | " << std::endl;
    std::cout << "              | " << std::endl;
    mpUp->PrintAncestry();
  } else {
    std::cout << std::endl;
    std::cout << "###############################" << std::endl;
  }
}
