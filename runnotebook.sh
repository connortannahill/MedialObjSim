module load python/3.8.10

salloc --time=1:0:0 --ntasks=1 --cpus-per-task=2 --mem-per-cpu=1024M --account=def-wan srun $VIRTUAL_ENV/bin/notebook.sh