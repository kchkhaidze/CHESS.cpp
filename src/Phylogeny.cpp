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

#include "Phylogeny.hpp"
#include "Cell.hpp"
#include "CellType.hpp"
#include <iostream>

// PhylogenyRoot ///////////////////////////////////////////////////////////////

// Statics:
unsigned int PhylogenyRoot::msNextId = 0;

// Constructor:
PhylogenyRoot::PhylogenyRoot(Cell* pCell) : mId(msNextId++) {
  mpRoot = new PhylogenyNode(pCell);
}

PhylogenyNode* PhylogenyRoot::Root(){ return mpRoot; }

//Destructor
PhylogenyRoot::~PhylogenyRoot() {
  delete mpRoot;
}

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
    mNumMutsGeneration(0),
    mNumStagedCells(0)
  {
    if (pCell->AssociatedNode() != 0) {
      pCell->AssociatedNode()->AssociatedCell(0);
    }
    mTypeId = pCell->Type()->Id();
    pCell->AssociatedNode(this);
  }

PhylogenyNode::PhylogenyNode(Cell* pCell, PhylogenyNode* pUp)
  : mId(msNextId++),
  mpCell(pCell),
  mpUp(pUp),
  mpLeft(0),
  mpRight(0),
  mNumMutsGeneration(0),
  mNumStagedCells(0)
  {
    mGeneration = pUp->Generation() + 1;
    if (pCell->AssociatedNode() != 0) {
      pCell->AssociatedNode()->AssociatedCell(0);
    }
    mTypeId = pCell->Type()->Id();
    pCell->AssociatedNode(this);
  }

//Destructor
PhylogenyNode::~PhylogenyNode() {
  D(std::cout << "PhylogenyNode::~PhylogenyNode()" << std::endl;)

  // Delete left node:
  if (mpLeft != 0) {
    delete mpLeft;
    mpLeft = 0;
  }

  // Delete right node:
  if (mpRight != 0) {
    delete mpRight;
    mpRight = 0;
  }

  // Delete associated cell:
  if (mpCell != 0) {
    delete mpCell;
    mpCell = 0;
  }

  // Clean pointers in the up node:
  if (mpUp != 0) {
    if (mpUp->LeftNode() == this) {
      mpUp->LeftNode(0);
    } else if (mpUp->RightNode() == this) {
      mpUp->RightNode(0);
    }
  }
}


// Getters:
unsigned int PhylogenyNode::Id() const {return mId;};
unsigned int PhylogenyNode::TypeId() const {return mTypeId;};
unsigned int PhylogenyNode::Generation() const {return mGeneration;};
unsigned int PhylogenyNode::NumMutations() const {return mNumMutsGeneration;};
unsigned int PhylogenyNode::NumStagedCells() const {return mNumStagedCells;};
Cell* PhylogenyNode::AssociatedCell() const {return mpCell;};
PhylogenyNode* PhylogenyNode::UpNode() const {return mpUp;};
PhylogenyNode* PhylogenyNode::LeftNode() const {return mpLeft;};
PhylogenyNode* PhylogenyNode::RightNode() const {return mpRight;};

std::vector <PhylogenyNode*> PhylogenyNode::NodeAncestry() {
  std::vector <PhylogenyNode*> result;
  PhylogenyNode *pCurrent = this;
  do {
    result.push_back(pCurrent);
  } while ((pCurrent = pCurrent->UpNode()) != 0);
  return result;
};


void PhylogenyNode::SequenceStagedNodes(CellSample &sample) const {
  // Skip if there are no staged nodes at this level anymore:
  if (mNumStagedCells == 0)
    return;

  // Convert node id to a hex
  std::stringstream hex_node_id;
  hex_node_id << std::hex << mId;

  // Sample each mutation in current node once:
  for (int i = 0; i < mNumMutsGeneration; i++) {
    // Convert mutation id to a hex string:
    std::stringstream hex_mutation_id;
    hex_mutation_id.str("");
    hex_mutation_id << std::hex << i;
    std::string full_mut_id = hex_node_id.str() + "X" + hex_mutation_id.str();

    // Append mutation information to output vectors:
    sample.AddMutation(full_mut_id, mTypeId, mNumStagedCells);
  }

  // Now parse the nodes at the next level if these nodes are not empty:
  if (mpLeft != 0)
    mpLeft->SequenceStagedNodes(sample);

  if (mpRight != 0)
    mpRight->SequenceStagedNodes(sample);

  return;
};


void PhylogenyNode::SampleStagedSingleCells(SingleCellSample &sample,
                                            std::vector<int> &stack) const
  {

  // Skip if there are no staged nodes at this level anymore:
  if (mNumStagedCells == 0)
    return;

  // Append to stack:
  int stack_pos = 0;
  if (mNumMutsGeneration > 0) {
    stack.push_back(mNumMutsGeneration);
    stack_pos = stack.size();
  }

  // First the right subtree.
  if (mpLeft != 0) {
    mpLeft->SampleStagedSingleCells(sample, stack);
  }

  // Then the left subtree.
  if (mpRight != 0) {
    mpRight->SampleStagedSingleCells(sample, stack);
  }

  // Check if a single cell has to be added.
  if (mpRight == 0 && mpLeft == 0 && mpCell != 0) { // At a true tip!

    // Create a cell id:
    std::stringstream cell_id;
    cell_id << std::hex << mpCell->Id() << std::dec << " ("
            << mpCell->X() << "," << mpCell->Y() << "," << mpCell->Z() << ")";

    // Insert a copy of the current stack as cell.
    sample.appendCell(stack, cell_id.str());
  }

  // Clean up and return.
  if (stack_pos != 0) {
    stack[stack_pos - 1] = 0;
  }

  return;
}


void PhylogenyNode::StagedNodesToStream(std::ostream& os) const {
  // Skip if there are no staged nodes at this level anymore:
  if (mNumStagedCells == 0)
    return;
  if (mpLeft == 0 && mpRight == 0) { // as leaf
    if (mpCell != 0 && mpCell->Type() != 0) { // proper cell
      os << "cell" << mpCell->Id() << "_"
         << "type" << mpCell->Type()->Id() << "_"
         << "x" << mpCell->X() << "_"
         << "y" << mpCell->Y() << "_"
         << "z" << mpCell->Z();
    } else {
      os << this;
    }
    os << ":" << mNumMutsGeneration;
  } else { // as node
    bool putLeft = mpLeft != 0 && mpLeft->NumStagedCells() > 0;
    bool putRight = mpRight != 0 && mpRight->NumStagedCells() > 0;
    os << "(";
    if (putLeft)
      mpLeft->StagedNodesToStream(os);
    if (putLeft && putRight)
      os << ",";
    if (putRight)
      mpRight->StagedNodesToStream(os);
    os << "):" << mNumMutsGeneration;
  }
}


std::ostream& operator<<(std::ostream& os, const PhylogenyNode& node) {
  if (node.mpLeft == 0 && node.mpRight == 0) { // "Print as Leaf"
    Cell* cell = node.mpCell;
    unsigned int cell_id = (cell==0)? 0 : cell->Id();
    unsigned int type_id = (cell==0 || cell->Type()==0)? 0 : cell->Type()->Id();
    os << cell_id << "_" << type_id << ":" << node.mNumMutsGeneration;
  } else { // Print as node?
    os << "(";
    if (node.mpLeft != 0)
      os << (*node.mpLeft);
    if (node.mpLeft != 0 && node.mpRight != 0)
      os << ",";
    if (node.mpRight != 0)
      os << (*node.mpRight);
    os << "):" << node.mNumMutsGeneration;
  }
  return os;
}


// Setters:
void PhylogenyNode::AddNewMutations(int num_new_mutations) {
  mNumMutsGeneration += num_new_mutations;
}


void PhylogenyNode::StageNodeForSequencing(){
  mNumStagedCells++;
  if (mpUp != 0) // Not root node
    mpUp->StageNodeForSequencing();
};


void PhylogenyNode::UnstageNodes() {
  if (mNumStagedCells > 0) {
    mNumStagedCells = 0;
    if (mpLeft != 0)
      mpLeft->UnstageNodes();
    if (mpRight != 0)
      mpRight->UnstageNodes();
  }
}


void PhylogenyNode::AssociatedCell(Cell* pCell) { mpCell = pCell; }
void PhylogenyNode::UpNode(PhylogenyNode* pNode) {mpUp = pNode;};
void PhylogenyNode::LeftNode(PhylogenyNode* pNode) {mpLeft = pNode;};
void PhylogenyNode::RightNode(PhylogenyNode* pNode) {mpRight = pNode;};
void PhylogenyNode::TypeId(unsigned int NewType) {mTypeId = NewType;};


// Output:
void PhylogenyNode::Print() const {
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

void PhylogenyNode::PrintAncestry() const {

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
