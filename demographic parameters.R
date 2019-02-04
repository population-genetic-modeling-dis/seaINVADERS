#### Demographic Parameters for Model ####

##Native Range Demographics
#from arlequin
japan.theta<-5.59781
japan.theta.sd<-2.73165

indonesia.theta<-7.63797
indonesia.theta.sd<-2.65535

atlantic.theta<-1.78966
atlantic.theta.sd<-0.60442

## Haplotype information
#All haplotype info below from JJohnson Spreadsheet
atlantic.hap<-c(194,172,24,19,14,19,15,1,1)
carib.hap<-c(159,405,3,34,0,0,0,0,0)
gomx.hap<-c(52,113,0,23,0,0,0,0,0)

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
