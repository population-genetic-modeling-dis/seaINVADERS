#### Parameters for space_invaders.R model #####
NP               <- 40 	#Number of processors to use
MODEL_DURATION   <- 15	#Number of years to run the model
NUM_BOOTSTRAPS   <- 2 #Number of simulations to run with each source

#### Population Genetic Parameters ####
species          <- "PtrVol"
source.name      <- "Indonesia"
#source.thetas    <- c(1.789662 + (seq(-2,2,2)*0.608607))  #Theta(S), n=36
#source.thetas    <- c(0.728371,20.218042,2.773071)  #Theta(k), n=36
source.thetas    <- c(7.716828,5.493651,20.218042)  #Theta(S,Pi,k), n=36
#source.hap      <-c()

destination.name <- "Atlantic"
destination.hap <-c(194,172,24,19,14,19,15,1,1)
#dest.gendiv      <- c(0.6742,0.0137)   #(mean,stdev)

#### Demographic Parameters for Model ####

MONTHS           <- 12	#Number of "months" or portions to divide a year into
thin             <- FALSE 	#Save only first and last month of simulated pops?

min_prop        <- 1 #This is the minimum proportion of successful recruits
max_prop        <- 1 #This is the maximum proportion of successful
prop_bins       <- 1 #ex. If prop_increment=4, then 4 different props will be run, ex: (.25,.5,.75,1)

min_f_number    <- 2 #This is the minimum # of females.
max_f_number    <- 1000 #This is the maximum # of females.
f_bins          <- 39 #ex. If f_increment=3, then three different numbers of female colonists will be run (min, min+(max-min)/2,max)

BIN             <- 12     # Number of different age-classes
JUVI.MORT       <- 0.165  # Juvenile mortality 
ADULT.MORT      <- 0.052  # Adult mortality 
ADULT.FRAC      <- 1   # Fraction of starting population that are adults* 
# * Based on the empirical estimates of Belize 2014 lionfish sample
#   on the forereef.

## Calculating Recruit per Individual
# This section takes egg and larval demographic parameters and calculates
# the monthly number of recruits per individual adult lionfish.
ADULT.FEM.FRAC  <- 0.49    # Proportion of adults that are females
ADULT.FEM.MAT   <- 0.79    # Proportion of females that are mature
FE              <- 194577  # Fecundity: number of eggs per female per month
ME              <- 0.31    # Egg mortality (days)
DE              <- 3       # Egg duration (days)
ML              <- 0.35    # Larval mortality (days)
DL              <- 27      # Larval duration (days)

K               <- 1000000	#Carrying capacity, used to modulate "birth rate" for logistic pop growth

#FE.sd           <- 1       # Standard dev in fecundity, if enabled and >0 then fecundity will be stochastic within a bootstrap (same fecundity per year for a rep)