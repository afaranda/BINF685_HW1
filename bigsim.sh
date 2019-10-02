#!/bin/sh
#SBATCH --job-name=Lattice_Sim
#SBATCH --mem=64000
#SBATCH --exclude=biomix38
echo $(date)
python3 main.py > sim_result.txt
echo $(date)

R CMD BATCH plot_histograms.R
