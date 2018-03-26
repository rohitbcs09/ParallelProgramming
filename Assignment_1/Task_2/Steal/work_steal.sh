#!/bin/bash
#SBATCH -J work_steal # job name
#SBATCH -o work_steal%j.o       # output and error file name (%j expands to jobID)
#SBATCH -n 1              # total number of mpi tasks requested
#SBATCH -p normal     # queue (partition) -- normal, development, etc.
#SBATCH -t 00:30:00        # run time (hh:mm:ss) - 0.5 hours
#SBATCH --mail-user=rokkumar@cs.stonybrook.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
ibrun ./steal 16 2 > work_steal_16_nodes_2.out 
