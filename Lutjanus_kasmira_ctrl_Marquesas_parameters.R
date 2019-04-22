#### Parameters for space_invaders.R model #####
NP <- 20  			#Number of processors to use
MONTHS <- 12   		#Number of "months" or portions to divide a year into
MODEL_DURATION <- 56 #4	#Number of years to run the model (haplotype data from 2014, introduction event in 1958)


#### Demographic Parameters for Model ####

np <- NP-1 #Number of processors to use
NUM_BOOTSTRAPS <- 100 #Number of simulations to run with each source
thin <- FALSE # Save only first and last month of simulated pops?
source.name <- "Marq"
destination.name <- "Hawaii"

#proportion.successful.recruits variables
min_prop <- .25 #This is the minimum proportion of successful recruits.
max_prop <- 1 #This is the maximum proportion of successful recruits.
prop_increment <- 2 #ex. If prop_increment=4, then the model will count (.25,.5,.75,1)

#initial.females variables
min_f_number <- 25 #This is the minimum # of females.
max_f_number <- 3025 #This is the maximum # of females.
f_increment <- 250 #ex. If f_increment=5, then model will count (5,10,15,20...etc)

##Native Range Demographics
#from arlequin

source.theta<-14.95801
source.theta.sd<-4.45541

## Haplotype information
source.hap<-c(3,1,0,0,0,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,1,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,0,0,0,2,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,1,1,0,1,0,0,0,1,0,0,1,0,0,1,1,0,1,1,0,1)
destination.hap<-c(10,3,1,2,2,3,2,9,3,3,0,0,3,1,2,0,0,1,0,2,4,2,1,4,3,0,4,3,1,1,2,2,0,0,0,6,4,2,1,2,0,2,1,1,1,0,2,2,3,0,1,1,1,1,0,1,1,3,1,4,1,1,1,3,1,1,3,0,5,1,1,2,1,1,0,1,1,1,1,4,1,1,2,1,1,1,2,1,4,1,1,1,1,2,2,1,1,1,1,3,1,1,1,1,1,1,1,1,3,1,1,1,1,1,1,1,1,1,1,3,1,1,1,1,1,2,0,4,1,2,1,1,2,0,0,3,0,4,18,0,0,1,0,1,2,0,5,1,2,1,0,2,3,1,2,2,4,0,0,3,0,3,4,1,0,2,3,0,2,1,0,0,1,4,1,1,0)

## Demographic Parameters 
# Most recently updated by Gray on 4_8_2019
BIN             <- 12     # Number of different age-classes
JUVI.MORT       <- 1-0.395  # Juvenile mortality 
ADULT.MORT      <- 1-0.0083  # Adult mortality 
ADULT.FRAC      <- 1-0   # Fraction of starting population that are adults* 
# * originally: Based on the empirical estimates of Belize 2014 lionfish sample
#   on the forereef. (1-0.96)
#update: assumption that fish released were all adults


## Calculating Recruit per Individual
# This section takes egg and larval demographic parameters and calculates
# the monthly number of recruits per individual adult lionfish.
# Most recently updated by Gray on 4_8_2019
ADULT.FEM.FRAC  <- 0.54128    # Proportion of adults that are females
ADULT.FEM.MAT   <- 0.6004    # Proportion of females that are mature
FE              <- 11986.11  # Fecundity - number of eggs per female per month
ME              <- 15.83    # Egg mortality (days)
DE              <- 0.778       # Egg duration (days)
ML              <- 0.34884    # Larval mortality (days)
DL              <- 20.4      # Larval duration (days)

## Additional parameters for genomic modeling DIS 2019 (not originally in this file)
#Most recently updated by Gray on 4_8_2019
SR              <- 0.0833    #Spawning rate (spawning events per month)
JD              <- 26.93     #Juvenile duration (months)
ARM             <- 2.3       #Age of reproductive maturity (years)

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
