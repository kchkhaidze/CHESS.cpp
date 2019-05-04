#' @title Get needle biopsy coordinates
#'
#' @description Prepares coordinates for cells inside a given needle shaped
#'      biopsies to be then sampled from a simulated tumour.
#'
#' @param radius Needle width radius.
#' @param grid Tumour space grid size.
#' @param nsamples Number of biopsies to sample.
#' @param dimension Tumour space grid dimension.
#'
#' @return List of cell coordinates for sampling.
#'
#' @author Kate Chkhaidze, \email{kate.chkhaidze@icr.ac.uk}
#'
#' @export
PrepareBiopsyNeedles <- function(radius, grid, nsamples, dimension) {

  
  x <- y <- grid
  xc <- c(x/4 + 10, x/2, 3*x/4 - 10)
  
  xc.lower <- xc - radius
  xc.upper <- xc + radius
  
  yc.lower <- rep(1, nsamples)
  yc.upper <- rep(grid, nsamples)
  
  if (dimension == '2D') {
    
    zc.lower <- zc.upper <- rep(1, length(xc)) 
    
  } else if (dimension == '3D') {
    
    zc <- rep(grid/2, length(xc))
    zc.lower <- zc - radius
    zc.upper <- zc + radius
    
  } else {
    print('Unrecognized input for dimension')
  }
  
  sample.strategy <- list()
  
  for (i in 1:nsamples) {
    sample.strategy[[length(sample.strategy) + 1]] <-
      c(xc.lower[i], yc.lower[i], zc.lower[i],
        xc.upper[i], yc.upper[i], zc.upper[i])
  }

  return(sample.strategy)
}


#' @title Get punch biopsy coordinates
#'
#' @description Prepares coordinates for cells inside a given square shaped
#'      biopsies to be then sampled from a simulated tumour.
#'
#' @param radius Square side times 1/2.
#' @param grid Tumour space grid size.
#' @param dimension Tumour space grid dimension.
#' 
#' @return List of cell coordinates for sampling.
#'
#' @author Kate Chkhaidze, \email{kate.chkhaidze@icr.ac.uk}
#'
#' @export
PrepareBiopsySquares <- function(radius, grid, dimension) {

  x <- y <- grid
  xc <- c(x/4, x/2, 3*x/4, x/2, x/2)
  yc <- c(y/2, 3*y/4, y/2, y/4, y/2)
  
  xc.lower <- xc - radius
  xc.upper <- xc + radius
  
  yc.lower <- yc - radius
  yc.upper <- yc + radius
  
  if (dimension == '2D') {
    
    zc.lower <- zc.upper <- rep(1, length(xc)) 
    
  } else if (dimension == '3D') {
    
    zc <- rep(grid/2, length(xc))
    zc.lower <- zc - radius
    zc.upper <- zc + radius
  
  } else {
    print('Unrecognized input for dimension')
  }

  sample.strategy <- list()
  
  for (i in 1:length(xc)) {
    sample.strategy[[length(sample.strategy) + 1]] <-
      c(xc.lower[i], yc.lower[i], zc.lower[i],
        xc.upper[i], yc.upper[i], zc.upper[i])
  }

  return(sample.strategy)
}


#' @title Get single cell sample coordinates
#'
#' @description Prepares coordinates for single cells to be sampled
#'    from a simulated tumour.
#'
#' @param x Tumour space grid dimension 1.
#' @param y Tumour space grid dimension 2.
#' @param z Tumour space grid dimension 3.
#' @param nsamples Number of single cells to sample.
#' @param dimension Tumour space grid dimension.
#'
#' @return List of single cell coordinates to sample.
#'
#' @author Kate Chkhaidze, \email{kate.chkhaidze@icr.ac.uk}
#'
#' @export
PrepareSingleCellSamples <- function(x, y, z, nsamples, dimension) {

  radius <- x/2
  sample.strategy <- data.frame()
  
  if (dimension == '2D') {
    while (nrow(sample.strategy) < nsamples) {
      
      x.sample <- sample(1:x, 1)
      y.sample <- sample(1:y, 1)
      
      if ((x.sample - (x - radius))^2 + 
          (y.sample - (x - radius))^2 < radius^2) {
        
        sample.strategy <- rbind(
          sample.strategy, data.frame(x.sample, y.sample))
      }
    }
    sample.strategy$z.sample <- z
    
  } else if (dimension == '3D') {
    while (nrow(sample.strategy) < nsamples) {
      
      x.sample <- sample(1:x, 1)
      y.sample <- sample(1:y, 1)
      z.sample <- sample(1:z, 1)
      
      if ((x.sample - (x - radius))^2 + 
          (y.sample - (x - radius))^2 + 
          (z.sample - (z - radius))^2 < radius^2) {
        
        sample.strategy <- rbind(
          sample.strategy, data.frame(x.sample, y.sample, z.sample))
      }
    }
  } else {
    print('Unrecognized input for dimension')
  }

  return(sample.strategy)
}
