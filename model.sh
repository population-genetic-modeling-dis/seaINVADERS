#!/bin/bash

#SBATCH -J seaNVDRS                 # Name of the job
#SBATCH -o seaNVDRS_%j.out     # Name of file that will have program output
#SBATCH -e seaNVDRS_%j.err     # Name of the file that will have job errors, if any
#SBATCH -N 1                      # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 20                     # Number of cores ( my test allocated 2 per node )
#SBATCH -p normal                 # Job queue    
                                  # (see available partitions and their number of nodes with sinfo --long -Node command )
#SBATCH -t 4-00:00:00	#Set time limit

module load openmpi
module load R/openmpi/intel/3.2.2

#To Run use command: sbatch model.sh NameOfParameters.R

PARAMS=$1
BS=$(grep '^NUM_BOOTSTRAPS' $PARAMS | sed 's/NUM_BOOTSTRAPS.*<-//g' | sed 's/#Number.*$//g' | sed -e "s/[[:space:]]//g")
YRS=$(grep '^MODEL_DURATION' $PARAMS | sed 's/MODEL_DURATION.*<-//g' | sed 's/#Number.*$//g' | sed -e "s/[[:space:]]//g")
SPECIES=$(grep '^species' $PARAMS | sed 's/species.*<-//g' | sed 's/\"//g' | sed -e "s/[[:space:]]//g")
SRC=$(grep '^source\.name' $PARAMS | sed 's/source\.name.*<-//g' | sed 's/\"//g' | sed -e "s/[[:space:]]//g")
DEST=$(grep '^destination\.name' $PARAMS | sed 's/destination\.name.*<-//g' | sed 's/\"//g' | sed -e "s/[[:space:]]//g")
THETA=$(grep '^source\.thetas' $PARAMS)
HAPS=$(grep '^source\.hap' $PARAMS)
if [[ ! -z $THETA && ! -z $HAPS ]]; then
	HAPsTHETA=ALL
elif [[ ! -z $THETA && -z $HAPS ]]; then
	HAPsTHETA=THETA
elif [[ -z $THETA && ! -z $HAPS ]]; then
	HAPsTHETA=HAPS
else
	echo "no source population provided, exiting"
	exit 1
fi
OUTDIR=$(paste -d _ <(echo output) <(echo $SLURM_JOB_ID) <(echo $SPECIES) <(echo $SRC) <(echo $DEST) <(echo $HAPsTHETA) <(echo ${YRS}YRS) <(echo ${BS}BS) )

mkdir $OUTDIR
cp $PARAMS $OUTDIR/$PARAMS

Rscript --verbose model.R $PARAMS $OUTDIR > $OUTDIR/seaINVADERS.Rout

mv seaNVDRS_${SLURM_JOB_ID}.out $OUTDIR
mv seaNVDRS_${SLURM_JOB_ID}.err $OUTDIR