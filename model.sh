#!/bin/bash
# Bash script used to run the Lionfish source model on HPC computer (5 nodes 100 CPUs)

#To Run use command: sbatch Cope_ERGM_Bash.shn
#SBATCH -J Number_LF          # Name of the job
#SBATCH -o Lion_number.out       # Name of file that will have program output
#SBATCH -e Lion_number.err       # Name of the file that will have job errors, if any
#SBATCH -N 1           # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 20              # Number of cores ( my test allocated 2 per node )
#SBATCH -p normal           # Partition    
                  # (see available partitions and their number of nodes with sinfo --long -Node command )
#SBATCH -t 4-00:00:00	#Set time limit

module load openmpi
module load R/openmpi/intel/3.2.2

mpirun -np 1 R CMD BATCH model.R 
