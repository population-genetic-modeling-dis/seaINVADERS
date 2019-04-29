#### Parameters for seaINVADERS.R model #####
NP               <- 20		#Number of processors to use
MODEL_DURATION   <- 2		#Number of years to run the model
NUM_BOOTSTRAPS   <- 2		#Number of simulations to run with each source

#### Population Genetic Parameters ####
species          <- "LutKas"
source.name      <- "Marquesas"
#source.hap       <- c(3,1,0,0,0,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,1,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,0,0,0,2,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,1,1,0,1,0,0,0,1,0,0,1,0,0,1,1,0,1,1,0,1)
source.hap       <- c(1,2,0,0,0,0,3,1,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,1,0,1,0,0,0,0,1,1,1,0,1,0,0,0,1,0,1,0,0,1,1,0,1,1,0,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,1,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#source.thetas     <- c(14.958010 + (seq(-2,2,2)* 4.455406)   #theta(S)
#source.thetas    <- c(149,376,1062)  #list of thetas to run, must comment out source.theta and source.theta.sd or this will be ignored, use Theta(k) from Arlequin

destination.name <- "Hawaii"
#destination.hap  <- c(10,3,1,2,2,3,2,9,3,3,0,0,3,1,2,0,0,1,0,2,4,2,1,4,3,0,4,3,1,1,2,2,0,0,0,6,4,2,1,2,0,2,1,1,1,0,2,2,3,0,1,1,1,1,0,1,1,3,1,4,1,1,1,3,1,1,3,0,5,1,1,2,1,1,0,1,1,1,1,4,1,1,2,1,1,1,2,1,4,1,1,1,1,2,2,1,1,1,1,3,1,1,1,1,1,1,1,1,3,1,1,1,1,1,1,1,1,1,1,3,1,1,1,1,1,2,0,4,1,2,1,1,2,0,0,3,0,4,18,0,0,1,0,1,2,0,5,1,2,1,0,2,3,1,2,2,4,0,0,3,0,3,4,1,0,2,3,0,2,1,0,0,1,4,1,1,0)
destination.hap  <- c(1,2,0,0,0,0,3,1,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,1,0,1,0,0,0,0,1,1,1,0,1,0,0,0,1,0,1,0,0,1,1,0,1,1,0,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,1,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
dest.gendiv      <- c(0.9926,0.0011)   #(mean,stdev)
#dest.theta       <- c(119.017997)  #c() list of thetas to run, use Theta(k) from Arlequin, comment out if haplotype vector adequately describes diversity of entire population

#### Demographic Parameters for Model ####
MONTHS           <- 12   	#Number of "months" or portions to divide a year into (age repro mat)
thin             <- FALSE	# Save only first and last month of simulated pops?

min_prop         <- 4		#This is the minimum proportion of successful recruits relative to that of the demographic parameters in the literature
max_prop         <- 10		#This is the maximum proportion of successful recruits relative to that of the demographic parameters in the literature
prop_bins        <- 2 		#ex. If prop_increment=4, then 4 different props will be run, ex: (.25,.5,.75,1)

min_f_number     <- 20 	#This is the minimum # of colonizing females.
max_f_number     <- 100 	#This is the maximum # of colonizing females.
f_bins           <- 19 		#ex. If f_increment=3, then three different numbers of female colonists will be run (min, min+(max-min)/2,max)

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
