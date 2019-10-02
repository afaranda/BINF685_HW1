## Peptide Lattice Model Overview

There are two major modules in this pipeline, pep_lat.py and darwin.py
The module "pep_lat.py" contains classes that implement the data structures
required to build the lattice model and generate boltzman distrubutions

The module "darwin.py" implements the "darwin" class which provides
the genetic algorithm used to find optimal peptide configurations. 

The module "main.py" invokes several monte-carlo simulations (at different
temperatures) to generate distributions.  It also invokes the genetic algorithim

The helper script "myrand.py" accounts for the lack of "random.choices" in python
versions < 3.6

The helper script "plot_histograms.R" generates histograms from simulated data

The wrapper script bigsim.sh is a shell script that will run the simulations, and
then generate histograms by calling "plot_histograms.R".  The standard out from
main.py is redirected to the file "sim_result.txt" which captures the output of
the genetic algorithm

## Instructions

1. Download and copy all files to a folder or Clone the git hub repoistory git@github.com:afaranda/BINF685_HW1


2. Update executable permissions on bigsim.sh (chmod a+x bigsim.sh) at the command prompt


3. Invoke bigsim.sh at the command prompt 

## Note
Occasionaly the python interpreter complains about double scalar overflow when
running the genetic algorithm.  This does not seem to interfere with the algorithm