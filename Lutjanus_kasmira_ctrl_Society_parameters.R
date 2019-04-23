#### Parameters for space_invaders.R model #####
NP               <- 7		#Number of processors to use
MODEL_DURATION   <- 56		#Number of years to run the model
NUM_BOOTSTRAPS   <- 10		#Number of simulations to run with each source

#### Population Genetic Parameters ####
source.name      <- "Society"
source.theta     <- 8.52248
source.theta.sd  <- 2.69315
source.hap       <-c(4,1,1,1,1,6,2,0,0,0,0,0,0,0,4,3,0,0,1,3,1,1,0,1,0,2,0,0,1,0,1,0,0,1,1,1,1,1,1)   #haplotype freqs

destination.name <- "Hawaii"
destination.hap  <-c(8,0,2,0,2,13,0,2,2,1,5,5,2,1,8,0,1,1,0,16,0,6,1,0,1,0,2,1,4,3,0,1,1,0,0,1,0,0,0)


#### Demographic Parameters for Model ####
MONTHS           <- 12   	#Number of "months" or portions to divide a year into (age repro mat)
thin             <- FALSE	# Save only first and last month of simulated pops?

min_prop         <- 8		#This is the minimum proportion of successful recruits relative to that of the demographic parameters in the literature
max_prop         <- 10		#This is the maximum proportion of successful recruits relative to that of the demographic parameters in the literature
prop_bins        <- 2 		#ex. If prop_increment=4, then 4 different props will be run, ex: (.25,.5,.75,1)

min_f_number     <- 500 	#This is the minimum # of colonizing females.
max_f_number     <- 3025 	#This is the maximum # of colonizing females.
f_bins           <- 3 		#ex. If f_increment=3, then three different numbers of female colonists will be run (min, min+(max-min)/2,max)

BIN             <- 12     	# Number of different age-classes
JUVI.MORT       <- 0.395  	# Juvenile mortality ***too high?
ADULT.MORT      <- 0.0083  	# Adult mortality 
ADULT.FRAC      <- 1   	# Fraction of starting population that are adults* 

ADULT.FEM.FRAC  <- 0.54128	# Proportion of adults that are females
ADULT.FEM.MAT   <- 0.6004	# Proportion of females that are mature
FE              <- 11986.11	# Fecundity - number of eggs per female per month ***too low?
ME              <- 0.1583	# Egg mortality (per day) 
DE              <- 0.778	# Egg duration (days)
ML              <- 0.34884	# Larval mortality (per day)  
DL              <- 20.4		# Larval duration (days) 

K               <- 1000000	#Carrying capacity, used to modulate "birth rate" for logistic pop growth

## Standard deviations of each parameter (not in Jason's model)
# Most recently updated by Gray on 4_8_2019
#SR.sd             <- ?          #Spawning rate standard deviation
#FE.sd              <- 6279.28    #Fecundinty standard deviation
#ME.sd             <- 2.5        #Egg mortality standard deviation
#DE.sd             <- 0.339      #Egg duration standard deviation
#ML.sd             <- 0.03453    #Larval mortality standard deviation
#DL.sd             <- 0.5        #Larval duration standard deviation
#JUVI.MORT.sd      <- ?          #Juvenile mortality standard deviation
#JD.sd             <- ?          #Juvenile duration standard deviation
#ARM.sd            <- ?          #Age of reproductive maturity standard deviation
#ADULT.MORT.sd     <- ?          #Adult mortality standard deviation
#ADULT.FRAC.sd     <- 0          #Proportion of starting population that are adults standard deviation
#ADULT.FEM.MAT.sd  <- ?          #Proportion of females that are mature standard deviation
#ADULT.FEM.FRAC.sd <- ?          #Fraction of population that are females


