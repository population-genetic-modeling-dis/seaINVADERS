#### Parameters for space_invaders.R model #####
NP <- 20  			#Number of processors to use
MONTHS <- 12   		#Number of "months" or portions to divide a year into
MODEL_DURATION <- 56	#Number of years to run the model


#### Demographic Parameters for Model ####

np <- NP-1 #Number of processors to use
NUM_BOOTSTRAPS <- 100 #Number of simulations to run with each source
thin <- FALSE # Save only first and last month of simulated pops?
source.name <- "Society"
destination.name <- "Hawaii"

#proportion.successful.recruits variables
min_prop <- .25 #This is the minimum proportion of successful recruits.
max_prop <- 1 #This is the maximum proportion of successful recruits.
prop_increment <- 2 #ex. If prop_increment=4, then the model will count (.25,.5,.75,1)

#initial.females variables
min_f_number <- 25 #This is the minimum # of females.
max_f_number <- 1525 #This is the maximum # of females.
f_increment <- 125 #ex. If f_increment=5, then model will count (5,10,15,20...etc)

##Native Range Demographics
#from arlequin

source.theta<-8.52248
source.theta.sd<-2.69315

## Haplotype information
source.hap<-c(4,1,1,1,1,6,2,0,0,0,0,0,0,0,4,3,0,0,1,3,1,1,0,1,0,2,0,0,1,0,1,0,0,1,1,1,1,1,1)
destination.hap<-c(8,0,2,0,2,13,0,2,2,1,5,5,2,1,8,0,1,1,0,16,0,6,1,0,1,0,2,1,4,3,0,1,1,0,0,1,0,0,0)

## Demographic Parameters
BIN             <- 12     # Number of different age-classes
JUVI.MORT       <- 1-0.395  # Juvenile mortality 
ADULT.MORT      <- 1-0.0083  # Adult mortality 
ADULT.FRAC      <- 1-0   # Fraction of starting population that are adults* 
# * Based on the empirical estimates of Belize 2014 lionfish sample
#   on the forereef.

## Calculating Recruit per Individual
# This section takes egg and larval demographic parameters and calculates
# the monthly number of recruits per individual adult lionfish.
ADULT.FEM.FRAC  <- 0.54128    # Proportion of adults that are females
ADULT.FEM.MAT   <- 0.6004    # Proportion of females that are mature
FE              <- 11986.11  # Fecundity - number of eggs per female per month
ME              <- 15.83   # Egg mortality (days)
DE              <- 0.778       # Egg duration (days)
ML              <- 0.34884    # Larval mortality (days)
DL              <- 20.4      # Larval duration (days)

## Standard deviations of each parameter (not in Jason's model)
# Most recently updated by Gray on 4_8_2019
#SR.sd              <- ?          #Spawning rate standard deviation
#FE.sd              <- 6279.28    #Fecundinty standard deviation
#ME.sd              <- 2.5        #Egg mortality standard deviation
#DE.sd              <- 0.339      #Egg duration standard deviation
#ML.sd              <- 0.03453    #Larval mortality standard deviation
#DL.sd              <- 0.5        #Larval duration standard deviation
#JUVI.MORT.sd       <- ?          #Juvenile mortality standard deviation
#JD.sd              <- ?          #Juvenile duration standard deviation
#ARM.sd             <- ?          #Age of reproductive maturity standard deviation
#ADULT.MORT.sd      <- ?          #Adult mortality standard deviation
#ADULT.FRAC.sd      <- 0          #Proportion of starting population that are adults standard deviation
#ADULT.FEM.MAT.sd   <- ?          #Proportion of females that are mature standard deviation
#ADULT.FEM.FRAC.sd  <- ?          #Fraction of population that are females

## Vectors and List for above ***** DO NOT CHANGE ****
#Mortalities and time in month for recovery
Demo.param   <- c(ADULT.MORT,JUVI.MORT,ADULT.FRAC,BIN)

#Demographic parameters used to calculate recruit per individual in monthly
#time steps
RPR          <- c(ADULT.FEM.FRAC,ADULT.FEM.MAT,FE,ME,DE,ML,DL)
