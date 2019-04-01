#### Parameters for space_invaders.R model #####
NP <- 20  			#Number of processors to use
MONTHS <- 12   		#Number of "months" or portions to divide a year into
MODEL_DURATION <- 2	#Number of years to run the model


#### Demographic Parameters for Model ####

np <- NP-1 #Number of processors to use
NUM_BOOTSTRAPS <- 10 #Number of simulations to run with each source

##Native Range Demographics
#from arlequin

source.theta<-7.63797
source.theta.sd<-2.65535

## Haplotype information
source.hap<-c(3,1,0,0,0,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,1,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,0,0,0,2,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,1,1,0,1,0,0,0,1,0,0,1,0,0,1,1,0,1,1,0,1)
destination.hap<-c(10,3,1,2,2,3,2,9,3,3,0,0,3,1,2,0,0,1,0,2,4,2,1,4,3,0,4,3,1,1,2,2,0,0,0,6,4,2,1,2,0,2,1,1,1,0,2,2,3,0,1,1,1,1,0,1,1,3,1,4,1,1,1,3,1,1,3,0,5,1,1,2,1,1,0,1,1,1,1,4,1,1,2,1,1,1,2,1,4,1,1,1,1,2,2,1,1,1,1,3,1,1,1,1,1,1,1,1,3,1,1,1,1,1,1,1,1,1,1,3,1,1,1,1,1,2,0,4,1,2,1,1,2,0,0,3,0,4,18,0,0,1,0,1,2,0,5,1,2,1,0,2,3,1,2,2,4,0,0,3,0,3,4,1,0,2,3,0,2,1,0,0,1,4,1,1,0)


## Demographic Parameters
BIN             <- 12     # Number of different age-classes
JUVI.MORT       <- 1-0.165  # Juvenile mortality 
ADULT.MORT      <- 1-0.052  # Adult mortality 
ADULT.FRAC      <- 1-0.96   # Fraction of starting population that are adults* 
# * Based on the empirical estimates of Belize 2014 lionfish sample
#   on the forereef.

## Calculating Recruit per Individual
# This section takes egg and larval demographic parameters and calculates
# the monthly number of recruits per individual adult lionfish.
ADULT.FEM.FRAC  <- 0.49    # Proportion of adults that are females
ADULT.FEM.MAT   <- 0.79    # Proportion of females that are mature
FE              <- 194577  # Fecundity - number of eggs per female per month
ME              <- 0.31    # Egg mortality (days)
DE              <- 3       # Egg duration (days)
ML              <- 0.35    # Larval mortality (days)
DL              <- 27      # Larval duration (days)


## Vectors and List for above ***** DO NOT CHANGE ****
#Mortalities and time in month for recovery
Demo.param   <- c(ADULT.MORT,JUVI.MORT,ADULT.FRAC,BIN)

#Demographic parameters used to calculate recruit per individual in monthly
#time steps
RPR          <- c(ADULT.FEM.FRAC,ADULT.FEM.MAT,FE,ME,DE,ML,DL)
