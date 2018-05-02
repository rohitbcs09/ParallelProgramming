#!/bin/bash
#SBATCH -J rot_A_rot_B           # job name
#SBATCH -o rot_A_rot_B_job_%j.o  # output and error file name (%j expands to jobID)
#SBATCH -N 2              
#SBATCH -n 2                    # total number of mpi tasks requested
#SBATCH -p skx-dev               # queue (partition) -- normal, development, etc.
#SBATCH -t 00:30:00              # run time (hh:mm:ss) - 0.5 hours
#SBATCH --mail-user=rokkumar@cs.stonybrook.edu
#SBATCH --mail-type=begin        # email me when the job starts
#SBATCH --mail-type=end          # email me when the job finishes
ibrun -n 2 ./rot_A_rot_B 8 
