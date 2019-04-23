rm(list=ls())
setwd("C:/Users/cbird/Documents/GCL/scripts/space_invaders")

parameterfile <- "Pterois_volitans_ctrl_parameters.R"
outputdir <- "outdir"

#tanslate input into variable used by model
arguments <- c(parameterfile,outputdir)

#____________________________________________________________________________
#goto model.R and run from line 9
#__________________________________________________________________________

parameters <- arguments[1]
outDir <- arguments[2]

library(parallel)
library(abind)
library(snow)
library(ggplot2)
#### Source Functions ####
source("model_functions.R")
source(parameters)

#### User Inputs from Parameters File ####
BS <- NUM_BOOTSTRAPS
RUN.MONTH <- MONTHS*MODEL_DURATION     # Number of 12 months * number of years to run model
initial.females <- seq(min_f_number,max_f_number,by=f_increment) # vector of number of starting female lionfish
proportion.successful.recruits<-seq(min_prop,max_prop,length.out=prop_increment)

if(exists("source.hap")){haplotype.sources<-list(source.name=source.hap)}

all.thetas<-source.theta+(2:-2)*source.theta.sd

for(i in all.thetas){
  if(i < 0.5){
    all.thetas[which(all.thetas == i)] <- 0.5
  }
}

source.dist<-as.list(all.thetas)
names(source.dist)<-paste('Sim',source.name,'theta',all.thetas,sep='.')

if(exists("source.hap")){
  if(exists("source.theta")){
    sources<-c(source.dist,haplotype.sources)
  } else {
    sources<-c(haplotype.sources)
  }
} else {
  if(exists("source.theta")){
    sources<-c(source.dist)
  } else {
    stop("Space Invaders Haulted, no source data provided")
  }
}

haplotype.destinations<-list(destination.name=destination.hap)

destinations<-c(haplotype.destinations)

for.workers<-ls()

#### Build MPI Cluster ####
cluster <- makeCluster(np,type="SOCK",outfile = paste(outDir,"/debug.txt",sep=""))
clusterExport(cl=cluster, as.list(for.workers),
              envir=environment())
clusterEvalQ(cluster, {library(abind)})

#### Run the model ####
#CEB bypassing for loop, setting prop succ rec to 1
v=0.1
clusterExport(cl=cluster, list('v'),envir=environment())
#CEB bypassing for loop, setting source pop to 1
S=2
clusterExport(cl=cluster, list('S'),envir=environment())

#_____________________________________________________________
#run.Model function is called, will only run for first number of starting females and 1 BS, running line by line
#_____________________________________________________________

FEMALE.START <- initial.females[1]
hap.num.start.freq <- sources[[S]]
#RUN.MONTH <- 
#Demo.param <-
#RPR <-
verbose <- FALSE
variable.RPR <- v
THIN <- thin

#   This constant represents the fraction of surviving larvae (recruits) per adult lionfish, based on the egg and 
#   larval mortalities, the duration of both stages, and the expected fecundity of an individual female lionfish.  
REC.PER.IND <- variable.RPR*RPR[1]*RPR[2]*RPR[3]*exp(-(RPR[4]*RPR[5]+RPR[6]*RPR[7])) # Recruits per individual

if(length(hap.num.start.freq)>1){
  hap.num.init<-length(hap.num.start.freq)
  if(sum(hap.num.start.freq!=1)){hap.num.start.freq<-hap.num.start.freq/sum(hap.num.start.freq)}
  s.f.sfreq<-c(rmultinom(1,FEMALE.START,hap.num.start.freq)) #Here is where drawing from initial population can be vectorized - change from 1 to >1
} else if (length(hap.num.start.freq)==1){
  s.f.sfreq<-rinfall(hap.num.start.freq,FEMALE.START)
  s.f.sfreq<-c(s.f.sfreq,rep(0,300-length(s.f.sfreq)))
  hap.num.init<-300
}

# Intialize the 4-d array where model output will be stored
model.output <- array(dim=c(hap.num.init,Demo.param[4],RUN.MONTH))
model.summary <- array(dim=c(hap.num.init,RUN.MONTH))
demo.summary <- array(dim=c(Demo.param[4],RUN.MONTH))



#### Model with monthly removals ####

for( k in 1:RUN.MONTH) {
  #full model with individuals output
  model.output[,,k] <- model.func(model.output,s.f.sfreq,k,REC.PER.IND,Demo.param[1],Demo.param[2],RPR[1],hap.num.init,verbose)
  #summary of all age classes by month
  model.summary[,k] <- apply(model.output[,,k],1,function(x) sum(x))
  demo.summary[,k] <- apply(model.output[,,k],2,function(x) sum(x))
  #adjusting individual summary to a proportion
  #model.summary.adj[,k] <- model.summary[,k]/sum(model.summary[,k])
  if(verbose){
    print(paste(".....Analysis",round(100*k/RUN.MONTH,digits=1),"% complete"))
  }
}

sapply(seq(1:RUN.MONTH), function(x) sum(model.output[,1:11,x]))
sapply(seq(1:RUN.MONTH), function(x) sum(model.output[,12,x]))
plot(sapply(seq(1:RUN.MONTH), function(x) sum(model.output[,1:12,x])))
model.output[,,100]
model.summary[,13]



#### Model with monthly removals ####

  #_____________________________________________________
  #model.func
  #_____________________________________________________

  M <- 1   # month of simulation
  REC.IND <- REC.PER.IND
  AM <- Demo.param[1]
  JM <- Demo.param[2]
  fem.perc <- RPR[1]
  hap.num <- hap.num.init
  
  temp.m <- model.output[,,M]
  # Initializes the matrix at month one with user defined starting female lionfish distributed
  # randomly across haplotype according to defined hap freqs.
  age_bins=dim(model.output)[2]  
  
  if(M==1){
    #temp.m[,1:11] <- 0  
    temp.m[,1:(age_bins-1)] <- 0
    #temp.m[,12]   <- s.f.sfreq , CEB: 12 is the number fed to BIN in parameters, dim(model.output)[2] is the same thing
    temp.m[,age_bins]   <- s.f.sfreq
    if(verbose){
      print("Model Initialization...")
    }
  }
  
  model.output[,,M] <- temp.m
  
  M <- 2
  
  temp.e <- model.output[,,M-1]
  
  for( j in 1:age_bins) { #For each age bin
    if(j == 1){
      var.temp<-list(fem.perc,REC.IND)
      numb.produced<-n.size(temp.e[,age_bins],var.temp)
      temp.m[,j]<-c(rmultinom(1,numb.produced,temp.e[,age_bins]))
    } else{
      if(j == age_bins){
        temp.m[,j]<-0
        
        if(sum(temp.e[,j-1]) > 0){
          # numb.juv.maturing <- n.size(temp.e[,j-1],JM)
          # temp.m[,j]<-temp.m[,j]+c(rmultinom(1,numb.juv.maturing,temp.e[,j-1]))
          temp.m[,j] <- n.size(temp.e[,j-1],JM)
        }
        
        if(sum(temp.e[,j]) > 0){
          #numb.adults.surviving<-n.size(temp.e[,j],AM)
          #if(numb.adults.surviving > 1000000){numb.adults.surviving <- 1000000}
          #temp.m[,j]<-temp.m[,j]+c(rmultinom(1,numb.adults.surviving,temp.e[,j]))
          temp.m[,j] <- temp.m[,j] + n.size(temp.e[,j],AM)
        }
        
      } else {
        #Mortality calculation with no stochasticity
        #temp.m[,,j] <- temp.e[,,j-1] * JM
        #Stochastically adds mortality across haplotypes based on previous months hap freq
        #hap.freq.e<-freq.convert(temp.e[,j-1])
        if(sum(temp.e[,j-1]) > 0){
          #t<-sample(1:hap.num,size=n.size(temp.e[l,,j-1],JM),replace=TRUE,prob = hap.freq.e)
          #temp.m[l,,j]<-table(factor(t,levels=1:hap.num)) 
          
          # numb.surviving<-n.size(temp.e[,j-1],JM)
          # temp.m[,j]<-c(rmultinom(1,numb.surviving,temp.e[,j-1])) 
          temp.m[,j] <- n.size(temp.e[,j-1],JM)
        } else {
          temp.m[,j]<-0
        }
      }
    }
  }

  model.output[,,M] <- temp.m
        
  
  #_____________________________________________________
  #end model.func
  #_____________________________________________________
  
  #summary of all age classes by month
  model.summary[,k] <- apply(model.output[,,k],1,function(x) sum(x))
  demo.summary[,k] <- apply(model.output[,,k],2,function(x) sum(x))
  #adjusting individual summary to a proportion
  #model.summary.adj[,k] <- model.summary[,k]/sum(model.summary[,k])
  if(verbose){
    print(paste(".....Analysis",round(100*k/RUN.MONTH,digits=1),"% complete"))
  }

