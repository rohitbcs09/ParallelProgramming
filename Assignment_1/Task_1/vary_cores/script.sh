#!/bin/bash

# Example submission for SLURM system
# Loop through multiple job submissions
# Pass environment variables to job script
NUMBERS=$(seq 41 68)

for NUM in ${NUMBERS}
do
	#sbatch batch_script_kij_i ${NUM}
	#sbatch batch_script_kij_ij ${NUM}
	sbatch batch_script_kij_j ${NUM}
done
