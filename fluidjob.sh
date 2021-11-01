#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --account=def-wan
#SBATCH --mem=250G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=connor.tannahill@uwaterloo.ca
#SBATCH --mail-type=ALL

module load python/3.8.10
module load scipy-stack


args="Q3ManyObject 3 1"
# args="objscale.txt"
echo "Args = $args"

scaletestcmd="python experiments.py < $args"
singletestcmd="bash compile_and_run.sh $args"


# cmd=$scaletestcmd
cmd=$singletestcmd
echo "Running cmds $cmd"
$cmd
