#' @title Uniform perturbation kernel
#'
#' @description Perturbs a given parameter value with a corresponding
#'    uniform kernel.
#'
#' @param value Parameter value to perturb.
#' @param maxval Upper limit for the uniform distribution.
#' @param minval Lower limit for the uniform distribution.
#' @param frac Scale factor for sigma.
#'
#' @return Perturbed parameter value.
#'
#' @author Kate Chkhaidze, \email{kate.chkhaidze@icr.ac.uk}
#'
#' @export
PerturbUniKernel <- function(value, maxval, minval, frac=2){

  sigma <- (maxval - minval)/frac
  try.result <- runif(1, value - sigma, value + sigma)

  while (try.result < 0) {
    try.result <- runif(1, value - sigma, value + sigma)
  }

  return(try.result)
}

#' @title Make unrooted non-binary phylogenetic tree rooted and binary
#'
#' @description Checks if a given phylo object is rooted and binary. If not, converts and returns an updated tree.
#'
#' @param tree ape phylo object
#'
#' @return Rooted binary version of a given phylogenetic tree.
#'
#' @author Kate Chkhaidze, \email{kate.chkhaidze@icr.ac.uk}
#'
#' @export
CheckRootBinary <- function(tree){
  
  if (!(is.rooted(tree))) {
    tree <- root(tree, Ntip(tree) + 1, resolve.root = T)
  }
  
  if (!(is.binary.tree(tree))) {
    tree <- multi2di(tree)
  }
  
  return(tree)
}


#' @title Find mode of a numeric or character vector
#'
#' @description Given a vector of numeric or character values, returns its mode (most frequent values). If all elements have a same frequency, returns NA.
#'
#' @param x Numeric or character vector.
#'
#' @return Mode of x.
#'
#' @author Kate Chkhaidze, \email{kate.chkhaidze@icr.ac.uk}
#'
#' @export
GetMode <- function(x) { 
  
  x.tbl <- table(x)
  x.max <- max(x.tbl)
  
  if (all(x.tbl == x.max)) {
    
    mode <- NA
    
  } else {
    
    if (is.numeric(x)) {
      
      mode <- as.numeric(names(x.tbl)[x.tbl == x.max])
      
    } else {
      
      mode <- names(x.tbl)[x.tbl == x.max]
    }
  }
  return(mode)
}
