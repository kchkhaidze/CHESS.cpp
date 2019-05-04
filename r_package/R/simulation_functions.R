#' @title Simulate tumour bulk samples
#' 
#' @description Simulates tumour growth followed by bulk sampling and exporting sequenced list of mutations.
#'
#' @param x Space size in dimension one.
#' @param y Space size in dimension two. If missing model is 1D.
#' @param z Space size in dimension three. If missing model is 2D.
#' @param birthrates Birth rates.
#' @param deathrates Death rates.
#' @param aggressions Aggression probability allowing cells to push or not.
#' @param push_power Cells' pushing probability scaled by their location on the grid.
#' @param mutation_rates Mutation rates.
#' @param clone_start_times Times to new subclone introduction.
#' @param kill_regrow_times Times to kill 99 percent of tumour cell population and regrow the tumour.
#' @param father Father cell ID.
#' @param clonal_mutations Number of clonal mutations to add post simulation and sampling.
#' @param seed Random state.
#' @param verbose Prints simulation progress.
#' @param depth Sequencing depth.
#' @param min_reads Minimum number of detectable reads.
#' @param min_vaf Filters mutations by given minimum VAF.
#' @param depth_model Depth model - fixed vs sampling sequencing depth from a Poisson distribuition with mean given by the 'depth' argument.
#' @param sample_strategy Bulk sample coordinates.
#'
#' @return List of matrices (corresponding to each bulk) of sampled and sequenced mutations.
#' 
#' @author Timon Heide, \email{timon.heide@icr.ac.uk}, Kate Chkhaidze, \email{kate.chkhaidze@icr.ac.uk}
#'
#' @export
SimulateTumourSample <-
  function(x = 100,
           y = 100,
           z = 1,
           birthrates = 1,
           deathrates = 0,
           aggressions = 0,
           push_power = max(x, y, z),
           mutation_rates = 10,
           clone_start_times = 0,
           kill_regrow_times = 1,
           father = 0,
           clonal_mutations = 0,
           seed = 42,
           verbose = FALSE,
           depth = 100, 
           min_reads = 2,
           min_vaf = 0.0,
           depth_model = 1,
           sample_strategy = list(c(1, 1, 1, x, y, z)))
  {
    
    param.matrix <- rbind(birthrates, 
                          deathrates, 
                          aggressions, 
                          push_power,
                          mutation_rates,
                          clone_start_times,
                          kill_regrow_times, 
                          father)
    
    univ <- new(CHESS_Universe_object,
                x, y, z, param.matrix, seed, verbose, clonal_mutations)
    
    res <- lapply(sample_strategy, function(coord) {
      a <- coord - 1
      univ$TakeSample(a[1], a[2], a[3], a[4], a[5], a[6],
                      depth, depth_model, min_reads, min_vaf, seed)
    })
    
    rm(univ)
    return(res)
  }

#' @title Simulate tumour single cell sampling phylogenetic trees
#' 
#' @description Simulates tumour growth followed by single cell sampling and constructing corresponding phylogenetyic trees.
#'
#' @param x Space size in dimension one.
#' @param y Space size in dimension two. If missing model is 1D.
#' @param z Space size in dimension three. If missing model is 2D.
#' @param birthrates Birth rates.
#' @param deathrates Death rates.
#' @param aggressions Aggression probability allowing cells to push or not.
#' @param push_power Cells' pushing probability scaled by their location on the grid.
#' @param mutation_rates Mutation rates.
#' @param clone_start_times Times to new subclone introduction.
#' @param kill_regrow_times Times to kill 99 percent of tumour cell population and regrow the tumour.
#' @param father Father cell ID.
#' @param clonal_mutations Number of clonal mutations to add post simulation and sampling.
#' @param seed Random state.
#' @param verbose Prints simulation progress.
#' @param strategy Single-cell sample coordinates.
#' @param nsamples Number of single-cell samples.
#'
#' @return List of matrices (corresponding to each bulk) of sampled and sequenced mutations.
#' 
#' @author Timon Heide, \email{timon.heide@icr.ac.uk}, Kate Chkhaidze, \email{kate.chkhaidze@icr.ac.uk}
#'
#' @export
SimulateTumourTree <-
  function(x = 100,
           y = 100,
           z = 1,
           birthrates = 1,
           deathrates = 0,
           aggressions = 1,
           push_power = max(x, y, z),
           mutation_rates = 10,
           clone_start_times = 0, 
           kill_regrow_times = 1, 
           father = 0,
           clonal_mutations = 0,
           seed = 42,
           verbose = FALSE,
           strategy = NULL,
           nsamples = 10)
  {
    
    param.matrix <- rbind(birthrates,
                          deathrates, 
                          aggressions,
                          push_power,
                          mutation_rates,
                          clone_start_times, 
                          kill_regrow_times,
                          father)
    
    univ <- new(CHESS_Universe_object,
                x, y, z, param.matrix, seed, verbose, clonal_mutations)
    
    if (is.null(strategy)) {

      res <- univ$SingleCellTreeN(nsamples, x, y, z)

    } else {

      res <- univ$SingleCellTree(strategy$x - 1,
                                 strategy$y - 1,
                                 strategy$z - 1)
    }

    trees <- lapply(res$tree_string, function(x) ape::read.tree(text = x))
    
    rm(univ)
    return(list(tree = trees, cellcount = res$cell_counts))
  }
