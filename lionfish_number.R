#### Recolonization model run script ####
## Jason Selwyn
## Last Mod: 2016-Sept-21

rm(list=ls())
library(parallel)
library(abind)
#library(Rmpi)
library(snow)
library(ggplot2)
#### Source Functions ####
source("model_functions.R")
source("demographic parameters.R")


#### User Inputs for lionfish_number.R ####
np <- 40-1 #Number of processors to use
BS <- 10 #Number of simulations to run with each source
RUN.MONTH <- 12*2     # Number of 12 months * number of years to run model
initial.females <- seq(1,8,by=1) # vector of number of starting female lionfish
var.rpr<-seq(0.25,1,length.out=4)


#### Setup Data ####
region.sources<-list(Atlantic=atlantic.hap,Caribbean=carib.hap,GoMx=gomx.hap)

all.thetas<-indonesia.theta+(2:-2)*indonesia.theta.sd
native.dist<-as.list(all.thetas)
names(native.dist)<-paste('Sim.Native','theta',all.thetas,sep='.')

sources<-c(native.dist,region.sources[-3])
destinations<-c(region.sources)

for.workers<-ls()


#### Build MPI Cluster ####
cluster <- makeCluster(np,type="SOCK",outfile = "debug.txt")
clusterExport(cl=cluster, as.list(for.workers),
              envir=environment())
clusterEvalQ(cluster, {library(abind)})

#### Run the model ####
for(v in var.rpr){
  clusterExport(cl=cluster, list('v'),envir=environment())
  for(S in 1:length(sources)){
    clusterExport(cl=cluster, list('S'),envir=environment())
    
    s<-parallel::parSapply(cl=cluster, initial.females, function(x) replicate(BS, run.Model(FEMALE.START=x,hap.num.start.freq=sources[[S]],RUN.MONTH,Demo.param,RPR,F,variable.RPR=v), 
                                                                              simplify = "array"), simplify = 'array')
    
    print('Finished Cluster Simulations')
    print(paste("s",format(object.size(s),units="Mb"),sep=' = '))
    
    s<-s[,!apply(s,2,remove.0.haps),,]
    print('Summary part 1 finished')
    print(paste("s",format(object.size(s),units="Mb"),sep=' = '))
    
    s<-aperm(s, c(1,2,4,3))
    print('Summary part 2 finished')
    print(paste("s",format(object.size(s),units="Mb"),sep=' = '))
    
    names(dim(s))<-c('Months','Haplotypes','Initial Females','Bootstrap')
    print('Summary part 3 finished')
    print(paste("s",format(object.size(s),units="Mb"),sep=' = '))
    
    f<-array(apply(s,c(3,4),freq.summary),dim=dim(s))
    print('Summary finished')
    print(paste("f",format(object.size(f),units="Mb"),sep=' = '))
    
    save.image(paste('./',names(sources)[S],'-rpr-',v,'_source.RData',sep=''))
    print('Model Image Saved')
    #### Make Model Overview Plots ####
    #Problems with this resulting in this error: 
    #/home/apps/R/gcc/3.3.2/lib64/R/bin/BATCH: line 60: 94950 Killed                  ${R_HOME}/bin/R -f ${in} ${opts} ${R_BATCH_OPTIONS} > ${out} 2>&1
    
    # pdf(paste('./',names(sources)[S],'-plots.pdf',sep=''),width=100,height=100,onefile = T)
    # 
    # total.plot<-plotting.model(s,initial.females,type='total',mem.redux=T)
    # print(total.plot)
    # frequency.plot<-plotting.model(f,initial.females,type='freq',mem.redux=T)
    # print(frequency.plot)
    # 
    # dev.off()
    # 
    # rm(s)
    # rm(total.plot)
    # rm(frequency.plot)
    # gc()
    # 
    # print('Finished simulation plots')
    
    #### Statistical Analysis ####
    
    pdf(paste('./',names(sources)[S],'-rpr-',v,'-statistical_plots.pdf',sep=''),onefile = T)
    
    for(D in 1:length(destinations)){
      stats.output<-model.statistics(destination.haplotypes=destinations[[D]],model.freq=f)
      
      statistical.plots<-plotting.statistics(stats.output,start.females=initial.females,title=paste(names(sources)[S],'to',names(destinations)[D]))
      
      print(statistical.plots[[1]])
      print(statistical.plots[[2]])
      print(statistical.plots[[3]])
      print(statistical.plots[[4]])
      print(statistical.plots[[5]])
      
      write.csv(statistical.plots[[6]],paste('./',names(sources)[S],'-to-',names(destinations)[D],'-rpr-',v,'-plot.data.csv',sep=''))
    }
    
    dev.off()
    print('Finished statistical plots')
    rm(s)
    rm(f)
    rm(statistical.plots)
    gc()
    
    print(paste('Finished',names(sources)[S]))
  }
}

#### Stop Cluster ####
stopCluster(cluster)
