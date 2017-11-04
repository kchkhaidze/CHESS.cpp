/*
//#define DEBUG
#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

#include "Cell.hpp"
#include <random>
#include <iostream>

// Phylogeny_Node //////////////////////////////////////////////////////////////

// Statics:
int Phylogeny_Node::N = 0;
int Phylogeny_Node::next_id = 1;


// Constructors
Phylogeny_Node::Phylogeny_Node(Cell *cell) // Create root node with a cell.
: cell(cell), up(0), left(0), right(0), gen(0), n_muts_gen(0)  {
}

Phylogeny_Node::Phylogeny_Node(Cell *cell, Phylogeny_Node *up)
: cell(cell), up(up), left(0), right(0) {
}

// Functions:
void Phylogeny_Node::record_new_muts(int n) {
  n_muts_gen+=n;
  cum_sum_muts+=n;
}



// Constructors:
Phylogeny::Phylogeny(Cell *cell) { // Create root node with cell.
  root = new Phylogeny_Node(cell);
}
*/
