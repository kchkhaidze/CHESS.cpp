#ifndef PHYLOGENY_H
#define PHYLOGENY_H

/*
// Forward declerations: ///////////////////////////////////////////////////////
class Phylogeny_Node;
class Cell;
class Universe;

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>
#include "boost/multi_array.hpp"
#include "CImg.h"

// Phylo ///////////////////////////////////////////////////////////////////////

class Phylogeny_Node {
  static int N, next_id;
  Cell *cell;
  Phylogeny_Node *up, *left, *right;
  int gen;
  int n_muts_gen, cum_sum_muts;

  public:
    // Constructors
    Phylogeny_Node(Cell *cell);
    Phylogeny_Node(Cell *cell, Phylogeny_Node *up);
    // Functions
    void record_new_muts(int n);
};

class Phylogeny {
  int id;
  Phylogeny_Node *root;

  public:
    Phylogeny(Cell *cell); // Create root node with cell.
};

*/

#endif // PHYLOGENY_H
