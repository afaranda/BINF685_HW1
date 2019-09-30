#!/bin/sh
#SBATCH --job-name=Lattice_Sim
#SBATCH --mem=64000

echo $(date)
python3 main.py
echo $(date)
