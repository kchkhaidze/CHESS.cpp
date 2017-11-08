#ifndef PHYLOGENY_H
#define PHYLOGENY_H


// Forward declerations: ///////////////////////////////////////////////////////
class PhylogenyNode;
class Cell;
class Universe;

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>
#include "boost/multi_array.hpp"
#include "CImg.h"

// Phylo ///////////////////////////////////////////////////////////////////////

class PhylogenyRoot {
    const unsigned int mId;
    static unsigned int msNextId; // Increments from one.
    PhylogenyNode* mpRoot;

  public:
    //Constructor
    PhylogenyRoot(Cell*); // Create root node with cell.
    PhylogenyNode* Root();
};


class PhylogenyNode {
  const unsigned int mId;
  static unsigned int msNextId; // Increments from one.

  Cell* mpCell;
  PhylogenyNode* mpUp;
  PhylogenyNode* mpLeft;
  PhylogenyNode* mpRight;

  unsigned int mGeneration;
  unsigned int mNumMutsGeneration;

  public:
    // Constructors:
    PhylogenyNode(Cell*);
    PhylogenyNode(Cell*, PhylogenyNode*);

    // Getters:
    unsigned int Id();
    unsigned int Generation();
    unsigned int NumMutations();
    Cell* AssociatedCell();
    PhylogenyNode* UpNode();
    PhylogenyNode* LeftNode();
    PhylogenyNode* RightNode();
    std::vector <PhylogenyNode*> NodeAncestry();

    // Setters:
    void AddNewMutations(int);
    void AssociatedCell(Cell*);
    void UpNode(PhylogenyNode*);
    void LeftNode(PhylogenyNode*);
    void RightNode(PhylogenyNode*);

    // Output:
    void Print();
    void PrintAncestry();

};

#endif // PHYLOGENY_H
