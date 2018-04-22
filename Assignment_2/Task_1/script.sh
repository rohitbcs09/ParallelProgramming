#!/bin/bash

# Example submission for SLURM system
# Loop through multiple job submissions
# Pass environment variables to job script
NUMBERS=$(seq 1 68)

for NUM in ${NUMBERS}
do
	#sbatch batch_script_kij_i ${NUM}
	#sbatch batch_script_kij_ij ${NUM}
	#./a_radix_sort ${NUM} 1048576
	./a_rand_quicksort ${NUM} 1048576 1024
        #echo ${NUM}
done
