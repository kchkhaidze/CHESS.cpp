#ifndef CANCERSIM_H
#define CANCERSIM_H

#include "Rcpp.h"

class Universe; // forward declare

class Universe_class_wrap {
  Universe* mpUniverse;

  public:

  // Creation:
  Universe_class_wrap(int, int, int, Rcpp::NumericMatrix, int, bool, int);

  // Destruction:
  ~Universe_class_wrap();

  // Sampling:
  Rcpp::List TakeSample(int,int,int, int,int,int, double,int,int,double, int);

  Rcpp::List SingleCellTree(std::vector<int>,
                            std::vector<int>,
                            std::vector<int>);

  Rcpp::List SingleCellTreeN(int, int, int, int);
};

#endif // CANCERSIM_H
