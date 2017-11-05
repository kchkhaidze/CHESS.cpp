#ifndef CELL_H
#define CELL_H

// Forward declerations: ///////////////////////////////////////////////////////
class CellType;
class Universe;
class PhylogenyNode;

//class Phylogeny_Node;

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>
#include <array>
#include "boost/multi_array.hpp"

// Cell ////////////////////////////////////////////////////////////////////////

class Cell {
    // Identifiers:
    const unsigned int mId;
    static unsigned int msNextId;

    // Location related:
    Universe* mpUniverse;
    int mX, mY, mZ;

    // Properties:
    CellType* mpType;
    PhylogenyNode* mpNode;

  public:
    // Constructors:
    Cell ();
    Cell (CellType*);

    // Destructor:
    ~Cell ();

    // Getter functions:
    unsigned int Id();
    CellType* Type();
    bool AsProgenitorDies();
    bool AsOffspringDies();
    bool TriesPush();
    PhylogenyNode* AssociatedNode();

    // Setter functions:
    void Location(Universe*, int, int, int);
    void Location(int, int, int);
    void Type(CellType*);
    void AssociatedNode(PhylogenyNode*);

    // Other functions:
    void MutateCell();
    void MutateCell(float);
    void Print();
    void Divide();
    bool RandomPush(int&, int&, int&);

};

#endif // CELL_H
