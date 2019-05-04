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
    int mN, mX, mY, mZ;
    int mTypeIndex;

    // Properties:
    CellType* mpType;
    PhylogenyNode* mpNode;
    int mTotalMutNum;

  public:
    // Constructors:
    Cell ();
    Cell (CellType*);

    // Destructor:
    ~Cell ();

    // Getter functions:
    unsigned int Id() const;
    int TypeIndex() const;
    int N() const;
    int X() const;
    int Y() const;
    int Z() const;
    int TotalMutNum() const;
    bool AsProgenitorDies() const;
    bool TriesPush() const;

    CellType* Type();
    PhylogenyNode* AssociatedNode();

    // Setter functions:
    void Location(Universe*, int, int, int, int);
    void Location(int, int, int, int);
    void Location(int, int, int);
    void Type(CellType*);
    void AssociatedNode(PhylogenyNode*);
    void TypeIndex(int);
    void TotalMutNum(int);

    // Other functions:
    void MutateCell();
    void MutateCell(double);
    void MutateCellFixedNumber(int);
    void Print() const;
    void Divide();
    void Kill();
    void DoAction(int*);
    bool RandomPush(int&, int&, int&);
};

#endif // CELL_H
