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

#ifndef PHYLOGENY_H
#define PHYLOGENY_H

// Forward declerations: ///////////////////////////////////////////////////////
class PhylogenyNode;
class Cell;
class Universe;

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>
#include "boost/multi_array.hpp"
#include "CellSample.hpp"

// Phylo ///////////////////////////////////////////////////////////////////////

class PhylogenyRoot {
    const unsigned int mId;
    static unsigned int msNextId; // Increments from one.
    PhylogenyNode* mpRoot;

  public:
    //Constructor
    PhylogenyRoot(Cell*); // Create root node with cell.
    PhylogenyNode* Root();
    //Destructor
    ~PhylogenyRoot(); // Create root node with cell.
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
  unsigned int mTypeId;
  int mNumStagedCells;

  public:
    // Constructors:
    PhylogenyNode(Cell*);
    PhylogenyNode(Cell*, PhylogenyNode*);


    //Destructor
    ~PhylogenyNode();

    // Getters:
    unsigned int Id() const;
    unsigned int TypeId() const;
    unsigned int Generation() const;
    unsigned int NumMutations() const;
    unsigned int NumStagedCells() const;
    Cell* AssociatedCell() const;
    PhylogenyNode* UpNode() const;
    PhylogenyNode* LeftNode() const;
    PhylogenyNode* RightNode() const;
    std::vector <PhylogenyNode*> NodeAncestry();
    friend std::ostream& operator<<(std::ostream&, const PhylogenyNode&);

    void SequenceStagedNodes(CellSample&) const;
    void SampleStagedSingleCells(SingleCellSample&, std::vector<int>&) const;
    void StagedNodesToStream(std::ostream& os) const;


    // Setters:
    void AddNewMutations(int);
    void AssociatedCell(Cell*);
    void UpNode(PhylogenyNode*);
    void LeftNode(PhylogenyNode*);
    void RightNode(PhylogenyNode*);
    void TypeId(unsigned int);
    void StageNodeForSequencing();
    void UnstageNodes();


    // Output:
    void Print() const;
    void PrintAncestry() const;
};

#endif // PHYLOGENY_H
