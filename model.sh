#!/bin/bash

#SBATCH -J SeaInv                 # Name of the job
#SBATCH -o SeaInvaders_%j.out     # Name of file that will have program output
#SBATCH -e SeaInvaders_%j.err     # Name of the file that will have job errors, if any
#SBATCH -N 1                      # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 40                     # Number of cores ( my test allocated 2 per node )
#SBATCH -p cbirdq                 # Partition    
                                  # (see available partitions and their number of nodes with sinfo --long -Node command )
#SBATCH -t 4-00:00:00	#Set time limit

module load openmpi
module load R/openmpi/intel/3.2.2

#To Run use command: sbatch model.sh NameOfParameters.R

TIME=$(date +%s%N | cut -b1-13)
OUTDIR=output_${1%.*}_$TIME
mkdir $OUTDIR
cp $1 $OUTDIR/$1

Rscript --verbose model.R $1 $OUTDIR > $OUTDIR/model.Rout