#' @title Multivariate distribution of multiple bulk samples
#'
#' @description Calculates multivariate joint distribution for a set of sampled
#'    bulk mutations.
#'
#' @param bulk.data.list List of sampled bulk mutations.
#' @param nbins Number of bins to build a multivariate histogram with.
#' @param nbulks Number of bulk samples.
#' @param add.clonal.muts Boolean to add random clonal mutations or not.
#'
#' @return Multivariate distribution of VAFs obtained for a given list of bulk samples.
#'
#' @author Kate Chkhaidze, \email{kate.chkhaidze@icr.ac.uk}
#'
#' @export
GetMultivrtDistr <- function(bulk.data.list,
                             nbins,
                             nbulks,
                             add.clonal.muts = F) {

  if (length(bulk.data.list) == 0) {
    return(NULL)
  }

  vafs.merged <- data.frame(matrix(
    ncol = 1, nrow = 0, dimnames = list(NULL, 'mut_id')))

  which.bulks <- sample(1:length(bulk.data.list), nbulks, replace = F)

  for (b in which.bulks) {

    mut.data <- bulk.data.list[[b]]$mutation_data

    if (nrow(mut.data) == 0) {

      vafs.merged[, paste('vaf_', b, sep = '')] <- rep(NA, nrow(vafs.merged))

    } else {

      mut.data$vaf <- mut.data$alt/mut.data$depth

      if (add.clonal.muts) {

        clonal.depth <- rpois(100, dph)
        clonal.vaf <- rbinom(100, clonal.depth, .5)/clonal.depth
        clonal.mutid <- paste0('clonal_', 1:length(clonal.vaf))
        mutid.vaf <- data.frame(c(mut.data$id, clonal.mutid),
                                c(mut.data$vaf, clonal.vaf))
      }

      mutid.vaf <- data.frame(mut.data$id, mut.data$vaf)
      colnames(mutid.vaf) <- c('mut_id', paste0('vaf_', b))
      vafs.merged <- merge(vafs.merged,
                           mutid.vaf,
                           by = 'mut_id',
                           all = TRUE)
    }
  }

  vafs.merged <- data.frame(apply(
    vafs.merged[, -1], 2, function(x){ x[is.na(x)] <- 0; x}))

  if (ncol(vafs.merged) > 1) {

    vafs.binned <-  data.frame(apply(
      vafs.merged, 2, function(x) cut(x,
                                      breaks = seq(0, 1, 1/nbins),
                                      right = 0)))

    vafs.binned.merged <- data.frame(table(apply(
      vafs.binned, 1, paste0, collapse = ',')))

    return(vafs.binned.merged)

  } else {

    return(NULL)
  }
}


#' @title VAF cumulative sums
#'
#' @description Calculates cumulative sums for a given VAF array.
#'
#' @param vaf Variant Allele Frequency of a sample.
#' @param fmin Minimum frequency to include.
#' @param fmax Maximum frequency to include.
#' @param step.size Step size for a cumulative walk.
#'
#' @return Vector of cumulative sums.
#' 
#' @author Kate Chkhaidze, \email{kate.chkhaidze@icr.ac.uk}
#'
#'
#' @export
GetCumSums <- function(vaf, fmin = .05, fmax = .7, step.size = .01) {

  steps <- seq(fmin, fmax, step.size)
  cumsum <- sapply(steps, function(x) sum(vaf <= x))
  cumsum <- cumsum - cumsum[1]

  return(cumsum)
}


#' @title VAF density counts
#'
#' @description Extracts bin counts from a given VAF histogram.
#'
#' @param vaf Variant Allele Frequency of a sample.
#' @param nbreaks Number of bins.
#'
#' @return Vector of histogram counts.
#'
#' @author Kate Chkhaidze, \email{kate.chkhaidze@icr.ac.uk}
#'
#' @export
GetDensityCounts <- function(vaf,
                             nbreaks = seq(0, 1, .01)) {

  return(hist(vaf, breaks = nbreaks, plot = F)$counts)
}
