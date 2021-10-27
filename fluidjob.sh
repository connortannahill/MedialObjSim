#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --account=def-wan
#SBATCH --mem=250G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=connor.tannahill@uwaterloo.ca
#SBATCH --mail-type=ALL
python experiments.py < scaletest3d.txt
# python experiments.py < objscale.txt
# bash compile_and_run.sh