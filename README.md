Cancer Heterogeneity with Spatial Simulations (CHESS)

Copyright (C) 2018, Timon Heide (timon.heide@icr.ac.uk)
                  & Kate Chkhaidze (Kate.Chkhaidze@icr.ac.uk)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.


Compilation: Install dependencies (make, gcc, boost, GLUT) and run "make".

Note: In the makefile and r_package/src/Makevars (or Makevars.win for windows users), please indicate the path where these dependencies will be installed. 

Usage: ./cancer_gillespie_simulation [options]
  Options: 

    -i ID, --sim_id=ID
        Simulation Id
          default: 1

    -x N, --size_x=N
        Size of dimension x.
          default: 100

    -y N, --size_y=N
        Size of dimension y.
          default: 100

    -z N, --size_z=N
        Size of dimension z.
          default: 100

    -m MU, --mutation_rates=MU
        Mutation rates of cell types.
        The number of new mutations
        will drawn randomly from a Poisson 
        distribution with lambda = MU during
        each division.
          default: 10

    -b F, --birthrates=F
        Birth rate of cell types.
          default: 1

    -d F, --deathrates=F
        Death rate of cell types.
          default: 0

    -r F, --aggressions=F
        Probabilty of cell types to push.
          default: 1

    -P N, --push_powers=N
        Push power of cell types.
          default: 4294967295

    -t F, --clone_start_times=F
        Time cell types are introduced.
          default: 1
          
    -w F, --kill_regrow_times=F
        Time to kill 99% of cells and regrow the tumour.
          default: 0 (never kill)

    -o DIR, --output_prefix=DIR
        Output prefix.
          default: ./

    -f FREQ, --display_frequency=FREQ
        Update frequency of the display.
          default: 0.05

    -s SEED, --seed=SEED
        Random seed.
          default: time(NULL)

    -M N, --number_clonal=N
        Number of clonal mutations.
          default: 0

    -D F, --depth=F
        Average sequencing depth.
            default: 100

    -X N, --depth_model=N
        Model for the sequencing depth.
          1) Poisson distributed depth.
          2) Beta-binomial distributed depth.
          3) Fixed depth.
            default: 1

    -R N, --min_reads=N
        Minimum number of alternate reads.
            default: 2

    -V F, --min_vaf=F
        Minimum variant allele frequency.
            default: 0

    -h, --help
        Print this help

Simulations can also be run from the CHESS R package.

To install and load the package, after compiling the C++ code:

```
R CMD INSTALL CHESS_0.0.0.9000.tar.gz
R
> library(CHESS)
```

The package depends on (and imports) the following R packages: Rcpp (>= 0.12.14), ape (>= 5.0), magrittr (>= 1.5), data.table (>= 1.11.4), transport (>= 0.9-4)

To run a tumour simulation that simulates a spatial tumour growth, followed by bulk and/or single cell sampling data generation, use the functions:

```SimulateTumourSample()``` and/or ```SimulateTumourTree()```

The functions output sample VAF and/or phylogenetic tree files, respectively.

The package implements three strategies of a tumour growth parameter inference (using the ABC SMC algorithm) by approximating a target tumour, using its:

```
1. Bulk samples (punch or needle biopsies) - ABCSMCwithBulkSamples()
2. Single cell sample phylogenetic trees - ABCSMCwithTreeSampleBL() and ABCSMCwithTreeSampleBT() (using Branch Lengths or Branching Times as summary statistics)
3. Whole tumour bulk sample - ABCSMCwithWholeTumour()
```

Depending on the strategy, a user would need to provide either target tumour bulk sample VAFs (list of R data.frames where each row should correspond to a unique mutation with the following columns: clone (Clone type label), alt (Number of reads), depth (Sequencing depth), id (Unique ID)), an array of whole tumour sample VAFs or single cell sampling phylogenetic trees. Alternatively, a user can provide a set of parameters (please refer to the package documentation for the details of each input parameter format) to simulate a synthetic target tumour to then recover these input parameters.  

The functions output sequence of files containing sets of inferred parameters corresponding to each SMC round (that can then be used to construct the posterior distributions of each parameter).
