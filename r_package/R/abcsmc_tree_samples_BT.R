#' @title ABC SMC with single cell sample trees using tree Branching Times as summary statistics
#'
#' @description Performs ABC SMC simulations to recover a given parameter set
#'    for a tumour simulation using its single cell sampling phylogenetic trees.
#'    
#' @param sim.id Simulated tumour ID.
#' @param grid Tumour space grid size.
#' @param mu.rl Mutation rate.
#' @param s.rl Selective advantage.
#' @param t.rl Time to a new subclone introduction.
#' @param d.rl Death rate.
#' @param a.rl Aggression (cell pushing probability scaled by its location on the grid).
#' @param random.state Random state.
#' @param nsamples Number of single cells to sample.
#' @param dimension Tumour space grid dimension.
#' @param output.dir Path to a directory where the output files will be exported.
#' @param sampling Sampling scheme.
#' @param params.to.rec String of parameter names (will be set as a directory name where the inferred particle files will be saved per SMC round).
#' @param ncores Number of cores to run in parallel with.
#' @param sigma Permutation kernel scale factor.
#' @param N Number of simulations to run per SMC round.
#' @param rounds Number of SMC rounds.
#' @param params Data.table containing parameter names, uniform prior distribution limits and flags to fix each parameter at given (previously inferred or already known) values or not. For exact format of the table, please see the examples section. 
#' @param target.tree Phylogenetic tree (phylo object) to recover the parameters for. 
#' @param target.popsize Estimated size (as number of cells) of the tumour to infer the parameters for.
#' 
#' @return Exports inferred set of parameters by each SMC round (as .csv files) into a given directory.
#'
#' @author Kate Chkhaidze, \email{kate.chkhaidze@icr.ac.uk}
#'
#' @export
ABCSMCwithTreeSamplesBT <- 
  function(sim.id = 1,
           grid = 50,
           mu.rl = 10,
           s.rl = 2,
           t.rl = 3,
           d.rl = 0,
           a.rl = 1,
           random.state = 123,
           nsamples = 10,
           dimension = '2D',
           output.dir = 'data/inferred_params/',
           sampling = 'tree',
           params.to.rec = 'mu_t_s_d_a',
           sigma = .5,
           N = 100,
           rounds = 10,
           ncores = 2,
           params = data.table(
             name = c('mu', 't', 's', 'd', 'a'),
             min = c(1, 1, 1, 0, 0),
             max = c(100, 20, 5, 1, 1),
             fix = c(F, F, F, F, F)),
           target.tree = NULL,
           target.popsize = NULL) {
    
  # Set path for the output data
  output.subdir <- paste0(
    output.dir,sampling,'/rec_',params.to.rec,'/T',sim.id,'/')
    
  dir.create(output.subdir, showWarnings = F, recursive = T)
  
  # Set value for z as the third dimension
  z <- ifelse(dimension == '3D', grid, 1)
  
  # Get a sampling strategy
  # sampling.strategy <- PrepareSingleCellSamples(
  #   x = grid, y = grid, z = z, nsamples = nsamples, dimension = dimension)

  if (is.null(target.tree)) {

    # SimulateTumourTree first simulates a tumour with a given parameter set,
    # then samples single cells and constructs a phylogenetic tree
    target.data <- SimulateTumourTree(
      x = grid,
      y = grid,
      z = z,
      birthrates = c(1, s.rl),
      deathrates = c(d.rl, d.rl),
      aggressions = c(1, 1),
      push_power = c(a.rl, a.rl)*grid,
      mutation_rates = c(mu.rl, mu.rl),
      clone_start_times = c(0, t.rl),
      father = c(0, 0),
      seed = random.state,
      nsamples = nsamples
    )

    # Collapse single node branches in order to have a binary tree
    target.tree <- ape::collapse.singles(target.data$tree[[1]])
    
    # Get population size
    target.popsize <- sum(target.data$cellcount)

  } else {

    # Collapse single node branches in order to have a binary tree
    target.tree <- ape::collapse.singles(target.tree)
  }
  
  # Replace 0 length branches by a relatively small number
  target.tree$edge.length[which(target.tree$edge.length <= 0)] <- 1e-010
  
  # Get tree branching times
  target.branchingtimes <-
    ape::branching.times(ape::chronopl(target.tree, lambda = .1))
  
  # Get cumulative sums of barnching times
  target.branchingtimes.cums <- GetCumSums(target.branchingtimes,
                                           fmin = 0, fmax = 1)

  keep.quantile <- c(1, rep(.5, rounds - 1)) # percentage of params to keep
  epsilion <- rep(NA, rounds) # ABC acceptance threshold

  particles.perround <- list()
  for (r in 1:rounds) {

    if (r == 1) {

      # Store proposed particles set per round
      proposed.particles <- do.call(
        rbind, parallel::mclapply(1:N, function(i) {
          cat('i ------------------------------> ', i, '\n')
          
          # Make sure only successful simulations are loged
          while (TRUE) {
  
            # Make sure simulated tumour population size is
            # at least the size of target population
            repeat {
              
              mu <- ifelse(params[name == 'mu']$fix,
                           mu.rl, 
                           runif(1,
                                 params[name == 'mu']$min,
                                 params[name == 'mu']$max))
              s <- ifelse(params[name == 's']$fix,
                          s.rl, 
                          runif(1,
                                params[name == 's']$min,
                                params[name == 's']$max))
              t <- ifelse(params[name == 't']$fix,
                          t.rl, 
                          runif(1,
                                params[name == 't']$min,
                                params[name == 't']$max))
              d <- ifelse(params[name == 'd']$fix,
                          d.rl, 
                          runif(1, 
                                params[name == 'd']$min,
                                params[name == 'd']$max))
              a <- ifelse(params[name == 'a']$fix,
                          a.rl, 
                          runif(1, 
                                params[name == 'a']$min,
                                params[name == 'a']$max))
  
              seed <- round(runif(1)*10^9)
  
              sim.data <- SimulateTumourTree(
                x = grid,
                y = grid,
                z = z,
                birthrates = c(1, s),
                deathrates = c(d, d),
                aggressions = c(1, 1),
                push_power = round(c(a, a))*grid,
                mutation_rates = c(mu, mu),
                clone_start_times = c(0, t),
                father = c(0, 0),
                seed = seed,
                nsamples = nsamples
              )
  
              sim.popsize <- sum(sim.data$cellcount)
              sim.tree.size <- length(sim.data$tree[[1]]$tip.label)
              cat('sim.popsize =', sim.popsize, 
                  '  sim.tree.size =', sim.tree.size, '\n')
              if (sim.popsize >= target.popsize) {
                break
              }
            }
  
            # Tidy simulated sample tree data same way as the target tree sample
            sim.tree <- ape::collapse.singles(sim.data$tree[[1]])
            sim.tree$edge.length[which(sim.tree$edge.length <= 0)] <- 1e-010
            sim.branchingtimes <-
              ape::branching.times(ape::chronopl(sim.tree, lambda = .1))
            
            # Calculate the Euclidean distance between 
            # summary statistics of the simulated and target tree samples
            dist <- as.numeric(dist(rbind(
              target.branchingtimes.cums, 
              GetCumSums(sim.branchingtimes, fmin = 0, fmax = 1))))
  
            # Store number of WT and mutant cells along with the simulation params
            wt.cells <- sim.data$cellcount$cellcounts[1]
            mt.cells <- sim.data$cellcount$cellcounts[2]
  
            curr.core.params <- data.frame(
              seed = seed,
              dist = dist,
              mu = mu,
              s = s,
              t = t,
              d = d,
              a = a,
              wt.cells = ifelse(is.na(wt.cells), 0, wt.cells),
              mt.cells = ifelse(is.na(mt.cells), 0, mt.cells)
            )
            return(curr.core.params)
          }
      }, mc.cores = ncores))

      # Simulation parameter sets are referred as particles
      previous.particles <-
        proposed.particles[complete.cases(proposed.particles),]

      # Set weights per particle (all weights are set to 1 in the first round)
      previous.particles$weight <- rep(1, nrow(previous.particles))

      # Order particles by their distance values
      previous.particles <- previous.particles[order(previous.particles$dist,
                                                     na.last = NA),]

      # Set epsilion to the distance up to the keep.quantile -th value
      epsilion[r] <- previous.particles$dist[N*keep.quantile[r]]

      # Update particles' list per round
      particles.perround[[r]] <- previous.particles

      # Update rumber of particles per round
      nprevparticles <- nrow(previous.particles)

      # Save particles per round
      fwrite(previous.particles,
             paste0(output.subdir, 'smc_particles_r', r, '.csv'))

    } else {
      # Is not the first ARC round

      # Get the number of particles from previous round
      N <- nprevparticles

      # Get the corresponding magnitude of permutation per pamater
      perturb.mu <- sigma*(max(
        particles.perround[[r - 1]]$mu) - min(particles.perround[[r - 1]]$mu))
      perturb.t <- sigma*(max(
        particles.perround[[r - 1]]$t) - min(particles.perround[[r - 1]]$t))
      perturb.s <- sigma*(max(
        particles.perround[[r - 1]]$s) - min(particles.perround[[r - 1]]$s))
      perturb.d <- sigma*(max(
        particles.perround[[r - 1]]$d) - min(particles.perround[[r - 1]]$d))
      perturb.a <- sigma*(max(
        particles.perround[[r - 1]]$a) - min(particles.perround[[r - 1]]$a))

      perturb.kern <- data.frame(
        mu = perturb.mu,
        t = perturb.t,
        s = perturb.s,
        d = perturb.d,
        a = perturb.a
      )

      # Repeat until the number of new proposed particles hits N
      proposed.particles <- data.frame()
      while (nrow(proposed.particles) < N) {

        nrunsleft <- N - nrow(proposed.particles)

        proposed.particles.part <- do.call(
          rbind, parallel::mclapply(1:nrunsleft, function(i) {
            cat('i ------------------------------> ', i, '\n')
           
            # Sample particles from previous round of particles' set,
            # perturb each particle parameter and run a new simulation
            while (TRUE) { 
              # Make sure simulated tumour population size is
              # at least the size of target population
              repeat {

                sample.particle.position <- sample(
                  N, 1, prob = particles.perround[[r - 1]]$weight)

                theta <-
                  particles.perround[[r - 1]][sample.particle.position,
                                              params$name]

                repeat {
                  proposed.mu <- 
                    runif(1, 
                          theta$mu - perturb.kern$mu,
                          theta$mu + perturb.kern$mu)
                  if (dunif(proposed.mu,
                            params[name == 'mu']$min,
                            params[name == 'mu']$max) > 0) {
                    break
                  }
                }
                repeat {
                  proposed.t <-
                    runif(1,
                          theta$t - perturb.kern$t,
                          theta$t + perturb.kern$t)
                  if (dunif(proposed.t,
                            params[name == 't']$min, 
                            params[name == 't']$max) > 0) {
                    break
                  }
                }
                repeat {
                  proposed.s <-
                    runif(1, 
                          theta$s - perturb.kern$s,
                          theta$s + perturb.kern$s)
                  if (dunif(proposed.s, 
                            params[name == 's']$min, 
                            params[name == 's']$max) > 0) {
                    break
                  }
                }
                repeat {
                  proposed.d <-
                    runif(1, 
                          theta$d - perturb.kern$d,
                          theta$d + perturb.kern$d)
                  if (dunif(proposed.d, 
                            params[name == 'd']$min, 
                            params[name == 'd']$max) > 0) {
                    break
                  }
                }
                repeat {
                  proposed.a <- 
                    runif(1,
                          theta$a - perturb.kern$a,
                          theta$a + perturb.kern$a)
                  if (dunif(proposed.a,
                            params[name == 'a']$min, 
                            params[name == 'a']$max) > 0) {
                    break
                  }
                }
                
                seed <- round(runif(1)*10^9)

                sim.data <- SimulateTumourTree(
                  x = grid,
                  y = grid,
                  z = z,
                  birthrates = c(1, proposed.s),
                  deathrates = c(proposed.d, proposed.d),
                  aggressions = c(1, 1),
                  push_power = round(c(proposed.a, proposed.a))*grid,
                  mutation_rates = c(proposed.mu, proposed.mu),
                  clone_start_times = c(0, proposed.t),
                  father = c(0, 0),
                  seed = seed,
                  nsamples = nsamples
                )

                sim.popsize <- sum(sim.data$cellcount)
                sim.tree.size <- length(sim.data$tree[[1]]$tip.label)
                cat('sim.popsize =', sim.popsize, 
                    '  sim.tree.size =', sim.tree.size, '\n')
                if (sim.popsize >= target.popsize) {
                  break
                }
              }

              # Process the sampled tree
              sim.tree <- ape::collapse.singles(sim.data$tree[[1]])
              sim.tree$edge.length[which(sim.tree$edge.length <= 0)] <- 1e-10
              sim.branchingtimes <-
                ape::branching.times(ape::chronopl(sim.tree, lambda = .1))
              
              # Calculate the Euclidean distance between 
              # summary statistics of the simulated and target tree samples
              potential.dist <- as.numeric(dist(rbind(
                target.branchingtimes.cums, 
                GetCumSums(sim.branchingtimes, fmin = 0, fmax = 1))))

              # If the potential distance is less than epsilion,
              # store the current particle
              if (potential.dist < epsilion[r - 1]) {

                wt.cells <- sim.data$cellcount$cellcounts[1]
                mt.cells <- sim.data$cellcount$cellcounts[2]

                curr.core.data <- data.frame(
                  seed = seed,
                  dist = potential.dist,
                  mu = proposed.mu,
                  s = proposed.s,
                  t = proposed.t,
                  d = proposed.d,
                  a = proposed.a,
                  wt.cells = ifelse(is.na(wt.cells), 0, wt.cells),
                  mt.cells = ifelse(is.na(mt.cells), 0, mt.cells))

                return(curr.core.data)
              } else {
                return(NULL)
              }
            }
          }, mc.cores = ncores))

        rownames(proposed.particles.part) <- NULL
        proposed.particles <- rbind(proposed.particles,
                                    proposed.particles.part)
      }

      # Archive proposed particles into previous particles
      proposed.particles[, 1:ncol(proposed.particles)] %<>%
        lapply(function(x) as.numeric(as.character(x)))
      previous.particles <- na.omit(proposed.particles)

      previous.particles <-
        previous.particles[order(previous.particles$dist, na.last = NA), ]

      # Update epsilion
      nprevparticles <- nrow(previous.particles)
      rowstokeep <- nprevparticles*keep.quantile[r]
      epsilion[r] <- previous.particles$dist[rowstokeep]
      particles.perround[[r]] <- previous.particles

      # Store variable param info separately
      params.var <- params[fix == FALSE]
      perturb.kern.var <- perturb.kern[
        colnames(perturb.kern) %in% params.var$name]
      
      # Calculate weights per particle
      weights.r <- rep(0, nprevparticles)
      for (n in 1:nprevparticles) {

        curr.particle <- as.numeric(
          particles.perround[[r]][n, params.var$name])

        w.numerator <- prod(dunif(
          curr.particle,
          params.var$min,
          params.var$max))

        w.denominator <- 0
        for (j in 1:nprevparticles) {

          prev.particle <- as.numeric(
            particles.perround[[r - 1]][j, params.var$name])
          
          w.denominator <- w.denominator +
            prod(dunif(curr.particle,
                       prev.particle - as.numeric(perturb.kern.var),
                       prev.particle + as.numeric(perturb.kern.var)))
        }

        weights.fraction <-
          w.numerator/(w.denominator*particles.perround[[r - 1]]$weight[n])

        if (!is.na(weights.fraction) & is.finite(weights.fraction)) {
          weights.r[n] <- weights.fraction
        }
      }
      
      # Set new weights for the new particle set
      particles.perround[[r]]$weight <- weights.r/sum(weights.r)

      # Save particles per round
      fwrite(particles.perround[[r]],
             paste0(output.subdir, 'smc_particles_r', r, '.csv'))
      
      # Update saved epsilions
      fwrite(data.frame(epsilion = epsilion),
             paste0(output.subdir, 'smc_epsilion.csv'))
    }
  }
}


