#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=def-wan
#SBATCH --mem=250G
#SBATCH --cpus-per-task=1
# python experiments.py < scaletest.txt
bash compile_and_run.sh