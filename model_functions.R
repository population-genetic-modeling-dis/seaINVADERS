#### Model script ####
## Alan Downey-Wall
## Last Mod: 2015-12-08

## Modified by Jason Selwyn
## 2016-Sep-19
#optimize to use rmultinom and remove self recruitment (assume full self-recruitment at the ocean basic scale)

# Description: Script for modeling the impact of self-recruitment and removal 
#              variability and there impact on recolonization on an arbitrary 
#              local scale over time.  Demographic parameters are default set 
#              to those of the red lionfish (Pterois volitans). 

run.Model <- function (FEMALE.START,hap.num.start.freq,RUN.MONTH,Demo.param,RPR,verbose,variable.RPR=1,THIN=T) {
  
  ## New Calculation
  #   This constant represents the fraction of surviving larvae (recruits) per adult, based on the egg and 
  #   larval mortalities, the duration of both stages, and the expected fecundity of an individual female lionfish.  
  if(RPR[8] > 0){
	FECUNDITY <- rnorm(1,RPR[3],RPR[8])
  } else { 
	FECUNDITY <- RPR[3]
  }
  REC.PER.IND <- variable.RPR*RPR[1]*RPR[2]*FECUNDITY*exp(-(RPR[4]*RPR[5]+RPR[6]*RPR[7])) # Recruits per individual
  
  #either read in hap freqs from parameters or calculate hap freqs from theta
  if(length(hap.num.start.freq)>1){
    hap.num.init<-length(hap.num.start.freq)
    if(sum(hap.num.start.freq!=1)){hap.num.start.freq<-hap.num.start.freq/sum(hap.num.start.freq)}
    
    #### Model Initializations ####
    
    ## Randomize initial female haplotype
    
    #s.f<-sample(1:hap.info.vec[2],size = hap.info.vec[1],replace = TRUE, prob= hap.num.start.freq)
    #s.f.sfreq<-table(factor(s.f, levels=1:hap.info.vec[2]))
    
    s.f.sfreq<-c(rmultinom(1,FEMALE.START,hap.num.start.freq)) #Here is where drawing from initial population can be vectorized - change from 1 to >1
    #print(s.f.sfreq)
    #Also need to change to have matrix of outputs not vector
    #s.f.sfreq<-s.f/hap.info.vec[1]
  } else if (length(hap.num.start.freq)==1){
    s.f.sfreq<-rinfall(hap.num.start.freq,FEMALE.START)
    s.f.sfreq<-c(s.f.sfreq,rep(0,300-length(s.f.sfreq)))
    hap.num.init<-300
  }

  
  # Intialize the 4-d array where model output will be stored
  model.output <- array(dim=c(hap.num.init,Demo.param[4],RUN.MONTH))
  
  
  #### Model with monthly removals ####

  for( k in 1:RUN.MONTH) {
    #full model with individuals output
    model.output[,,k] <- model.func(model.output,s.f.sfreq,k,REC.PER.IND,Demo.param[1],Demo.param[2],RPR[1],hap.num.init,verbose, RPR[9])
    #summary of all age classes by month
    #model.summary[,k] <- apply(model.output[,,k],1,function(x) sum(x))
    #adjusting individual summary to a proportion
    #model.summary.adj[,k] <- model.summary[,k]/sum(model.summary[,k])
    if(verbose){
      print(paste(".....Analysis",round(100*k/RUN.MONTH,digits=1),"% complete"))
    }
  }

    #names(hap.num.start.freq)
    # dimnames(model.summary)[[1]]<-SR.range
    # dimnames(model.summary)[[2]]<-names(hap.num.start.freq)
    # dimnames(model.summary)[[3]]<-seq(from=1,to=Demo.param[3],by=1)
    # dimnames(model.summary.adj)[[1]]<-SR.range
    # dimnames(model.summary.adj)[[2]]<-names(hap.num.start.freq)
    # dimnames(model.summary.adj)[[3]]<-seq(from=1,to=Demo.param[3],by=1)
    #list.models<-list(model.output,model.summary,model.summary.adj)
    #names <- c("model.output","model.summary","model.summary.adj")
    #names(list.models) <- names
  
    if(verbose){
      print("Analysis complete: Model run successful.") 
    }
  #return(list.models)
  #model.output
  output<-apply(model.output,c(3,1),sum)
  
  if(THIN){
    output<-month.thinning(output,RUN.MONTH)
  }
  output
}

#### Function that populates model output array for month M ####

# Description: Takes list of demographically relate/local distrubance parameters
# and models population recovery for month M given variable distance (removal) and 
# self-recruitment rates.

model.func <- function(model.output,s.f.sfreq,M,REC.IND,AM,JM,fem.perc,hap.num,verbose,K){
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
  } else {
    temp.e <- model.output[,,M-1]
    for( j in 1:age_bins) { #For each age bin
      if(j == 1){
         #Stochastically adds mortality across haplotypes based on previous months hap freq
        
		logistic_REC.IND <- REC.IND - (((REC.IND - 1)/K) * sum(temp.e[,age_bins]))
		var.temp<-list(fem.perc,logistic_REC.IND)
		
		#CEB: To have non-reproductive seasons (i.e. Months 5,6,7) Skip the next 2 lines when M==5:7
        #numb.produced<-n.size(temp.e[,age_bins],REC.IND)
		numb.produced <- n.size(temp.e[,age_bins],var.temp)
		if(sum(temp.e[,age_bins]) > 0 ){
			temp.m[,j]<-c(rmultinom(1,numb.produced,temp.e[,age_bins]))
		} else {
			temp.m[,j] <- 0
		}
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
  }
  return(temp.m)
}

n.size<-function(x,m){
  if(is.list(m)){
  #Creation of new fish - m[[1]] = percent females, m[[2]] = recruits per individual
    y <- sum(x)*m[[1]]*m[[2]]
	if(y < 1){
		y <- rbinom(1,1,y)  #use binom dist to determine if any female offspring are created
	} else {
		y <- round(y)
	}
	#stochastic alternative to the above code, whether offspring are male or female determined by bionomial dist
	#if(sum(x)*m[[2]] >= 1){
		#y <- rbinom(1,round(sum(x)*m[[2]]),m[[1]])
	#} else {
		#y <- rbinom(1,1,sum(x)*m[[1]]*m[[2]])
	#}
  }
  else{
  #Mortality rate of a given mortality rate = m
    # y <- sum(x)*m
	# if(y < 1){
		# y <- rbinom(1,1,y)  #use binom dist to determine if any females survive
	# } else {
		# y <- round(y)
	# }
	#stochastic alternative to the above code, survival determined by binomial dist
	#y <- rbinom(1, sum(x), prob = m)
	y <- sapply(x, function(x) rbinom(1, x, m))
  }
  return(y)
}

#### Plot outputs ####
#model.array=full.model.sum;start.females=initial.females;type='overall'
plotting.model<-function(model.array,start.females,type,bootstrap_culling=0.95,mem.redux=T,months=RUN.MONTH){
  #type can either be sum or freq or overall
  #start.females is a vector of starting number of females
  #assumed to be bootstrapped and with multiple female starts
  library(ggplot2,quietly=T)
  
  if(mem.redux){
    culling<-round(seq(1,length(start.females),length.out=10))
    
    model.array<-model.array[,,culling,bootstrap_culling<runif(dim(model.array)[4],0,1)]
    start.females<-start.females[culling]
  }
  
  if(type!='overall'){
    plot.data<-data.frame(month=rep(c(1,seq(2,months,by=1)),dim(model.array)[2]*dim(model.array)[3]*dim(model.array)[4]),
                          haplotype.frequency=matrix(model.array),
                          haplotype=as.factor(rep(rep(rep(1:dim(model.array)[2],each=dim(model.array)[1]),dim(model.array)[3]),dim(model.array)[4])),
                          females=rep(as.factor(rep(start.females,each=dim(model.array)[2]*dim(model.array)[1])),dim(model.array)[4]),
                          bootstrap=rep(1:dim(model.array)[4],each=prod(dim(model.array)[1:3])))
    
    #Remove haplotypes with only 0 in all simulations
    for(i in 1:nlevels(plot.data$haplotype)){
      if(all(plot.data$haplotype.frequency[plot.data$haplotype==levels(plot.data$haplotype)[i]]==0)){
        plot.data<-plot.data[plot.data$haplotype!=levels(plot.data$haplotype)[i],]
      }
    }
    plot.data$haplotype<-as.factor(as.numeric(plot.data$haplotype))
    
    p1<-ggplot(data=plot.data,aes(x=month,y=haplotype.frequency,col=haplotype,group=paste(bootstrap,haplotype)))+
      geom_line()+
      facet_wrap(~haplotype+females,labeller = "label_both",nrow=dim(model.array)[2],ncol=dim(model.array)[3])+
      theme_bw()+theme(legend.position="none")+xlab('Month')+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    if(type=='freq'){p1<-p1+ylab('Haplotype Frequency')}
    if(type=='total'){p1<-p1+ylab('Total Lionfish')+scale_y_log10()}
  }
  
  p1
}


#### Statistics ####
#destination.haplotypes=atlantic.hap;source.haplotypes = hap.num.start.freq;model.freq=full.model.freq
model.statistics<-function(destination.haplotypes, model.freq, dest.params){

  #Draw sample from bootstrap distributions with number = destination haplotype count
  sample.sum<-0*model.freq[dim(model.freq)[1],,,]
  for(invaders in 1:dim(model.freq)[3]){
    for(boot in 1:dim(model.freq)[4]){
		#handle NaNs
		if(model.freq[dim(model.freq)[1],,invaders,boot][1] == "NaN"){
			model.freq[dim(model.freq)[1],,invaders,boot]<-0
			model.freq[dim(model.freq)[1],,invaders,boot][1]<-1
		}
		#this line simulates the genetic sample from the colonizing population
      sample.sum[,invaders,boot]<-rmultinom(1,sum(destination.haplotypes),model.freq[dim(model.freq)[1],,invaders,boot])
    }
  }
  
  #comment out?  CEB
  if(length(destination.haplotypes)!=dim(model.freq)[2]){
    destination.haplotypes<-c(destination.haplotypes,rep(0,abs(dim(model.freq)[2]-length(destination.haplotypes))))
  }
  
  #Bootstrap to simulate drawing genetic sample from colonizing population for observed gen div   
  boot.dest.hap<-rmultinom(1000,sum(destination.haplotypes),destination.haplotypes)
  
  #Richness
  HN.calc<-function(x){
    sum(x>0)
  }
  
  obs.Hn<-HN.calc(destination.haplotypes)
  
  sim.Hn<-apply(sample.sum,c(2,3),HN.calc)
  prob.Hn<-rowSums(sim.Hn==obs.Hn)/ncol(sim.Hn)
  #if(all(prob.Hn!=0)){prob.Hn<-prob.Hn/sum(prob.Hn)}
  hap<-list(prob.Hn,sim.Hn,obs.Hn)
  
  #Diversity
  # HS.calc<-function(x){
    # 1-sum((x/sum(x))^2)
  # }
  HS.calc<-function(x){
	p <- x/sum(x)
	sum(x)*(1-sum(p^2))/(sum(x)-1)
  } 
  
  obs.Hs<-HS.calc(destination.haplotypes)
  
  #decide how to define confidence window around observed destination gen div
  if(dest.params == NULL){
	  #assume that dest haps represent whole pop, not just a sample
	  boot.Hs <- apply(boot.dest.hap,2,HS.calc)
	  Hs.CI <- quantile(boot.Hs,c(0.025,0.975))
  } else if(length(dest.params == 2)) {
	  #use the est gendiv and stdev of the est
	  #dest.params holds gen div and stdev of gen div
	  Hs.CI <- c(obs.Hs - 1*dest.params[2], obs.Hs + 1*dest.params[2])
  } else {
	  #use est of theta(k) to simulate drawing samples from dest population
	  #dummy code, for now
  	  boot.Hs <- apply(boot.dest.hap,2,HS.calc)
	  Hs.CI <- quantile(boot.Hs,c(0.025,0.975))
  }

  #simulate samples from destination population for calculating gen div
  sim.Hs<-apply(sample.sum,c(2,3),HS.calc)
  
  prob.Hs<-rowSums(sim.Hs>=Hs.CI[1] & sim.Hs<=Hs.CI[2])/ncol(sim.Hs)
  #if(all(prob.Hs!=0)){prob.Hs<-prob.Hs/sum(prob.Hs)}
  
  diversity<-list(prob.Hs,sim.Hs,unname(c(obs.Hs,Hs.CI)))
  
  #Joint probability Richness + Diversity
  joint<-rowSums(sim.Hn==obs.Hn & sim.Hs>=Hs.CI[1] & sim.Hs<=Hs.CI[2])/ncol(sim.Hn)
  #if(all(joint!=0)){joint<-joint/sum(joint)}
  
  output<-list(hap,diversity,joint)
  output
  
}


#### Plot Statistics ####
#start.females<-initial.females
plotting.statistics<-function(stats.output,start.females,title){
  library(ggplot2,quietly=T)
  

  #Output[[1]][[1]] = logLikelihood (Haplotype Richness)
  
  #Output[[1]][[2]] = bootstrap Haplotype Richness

  #Output[[1]][[3]] = observed Haplotype Richness

   
  bounds.Hn<-t(apply(t(stats.output[[1]][[2]]),2,quantile,c(0.025,0.975)))
  bounds.Hs<-t(apply(t(stats.output[[2]][[2]]),2,quantile,c(0.025,0.975)))
   
   plot.data<-data.frame(starting.females=start.females,
                        
                         mean.Hn=apply(stats.output[[1]][[2]],1,mean),
                         median.Hn=apply(stats.output[[1]][[2]],1,median),
                         lower.Hn=bounds.Hn[,1],upper.Hn=bounds.Hn[,2],
                         logLik.Hn=stats.output[[1]][[1]],
                         
                         mean.Hs=apply(stats.output[[2]][[2]],1,mean),
                         median.Hs=apply(stats.output[[2]][[2]],1,median),
                         lower.Hs=bounds.Hs[,1],upper.Hs=bounds.Hs[,2],
                         logLik.Hs=stats.output[[2]][[1]],
                         
                         joint=stats.output[[3]])
   
   #Parameter estimate
   max.prob.Hn<-plot.data$starting.females[which(plot.data$logLik.Hn==max(plot.data$logLik.Hn))]
   max.prob.Hs<-plot.data$starting.females[which(plot.data$logLik.Hs==max(plot.data$logLik.Hs))]
   max.prob.joint<-plot.data$starting.females[which(plot.data$joint==max(plot.data$joint))]
   
   Hn.plot<-ggplot(data=plot.data,aes(x=starting.females,y=mean.Hn))+
     geom_point()+geom_errorbar(aes(ymin=lower.Hn, ymax=upper.Hn))+
     geom_point(aes(y=median.Hn),col='red')+
     geom_hline(yintercept=stats.output[[1]][[3]],col='red',lty=2)+
     theme_bw()+theme(legend.position="none")+xlab('Number of Female Colonists')+ylab('Haplotype Richness')+ggtitle(title)+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
   #Hn.plot
   
   
   Hn.probability.plot<-ggplot(data=plot.data,aes(x=starting.females,y=logLik.Hn))+
     geom_point()+
     theme_bw()+theme(legend.position="none")+
     xlab(paste('Number of Female Colonists','\n','Most probable number of female colonists = ',max.prob.Hn,sep=''))+
     ylab('Probability of Having Observed # of Haplotypes')+ggtitle(title)+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
   #Hn.probability.plot
   
   Hs.plot<-ggplot(data=plot.data,aes(x=starting.females,y=mean.Hs))+
     geom_point()+geom_errorbar(aes(ymin=lower.Hs, ymax=upper.Hs))+
     geom_point(aes(y=median.Hs),col='red')+
     geom_hline(yintercept=stats.output[[2]][[3]][1],col='red',lty=2)+
     geom_ribbon(ymax=stats.output[[2]][[3]][3],ymin=stats.output[[2]][[3]][2],fill='red',alpha=0.2)+
     theme_bw()+theme(legend.position="none")+xlab('Number of Female Colonists')+ylab('Haplotype Diversity')+ggtitle(title)+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
   #Hs.plot
   
   
   Hs.probability.plot<-ggplot(data=plot.data,aes(x=starting.females,y=logLik.Hs))+
     geom_point()+
     theme_bw()+theme(legend.position="none")+
     xlab(paste('Number of Female Colonists','\n','Most probable number of female colonists = ',max.prob.Hs,sep=''))+
     ylab('Probability of Having Observed Haplotype Diversity')+ggtitle(title)+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
   #Hs.probability.plot
   
   # Fst.plot<-ggplot(data=plot.data,aes(x=starting.females,y=mean.Fst))+
   #   geom_point()+geom_errorbar(aes(ymin=lower.Fst, ymax=upper.Fst))+
   #   geom_point(aes(y=median.Fst),col='red')+
   #   #geom_hline(yintercept=stats.output[[1]][[3]],col='red',lty=2)+
   #   theme_bw()+theme(legend.position="none")+xlab('Starting Number of Females')+ylab('Fst')+ggtitle(title)+
   #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
   # #Fst.plot
   
   joint.probability.plot<-ggplot(data=plot.data,aes(x=starting.females,y=joint))+
     geom_point()+
     theme_bw()+theme(legend.position="none")+
     xlab(paste('Number of Female Colonists','\n','Most probable number of female colonists = ',max.prob.joint,sep=''))+
     ylab('Probability of Having Observed Haplotype Diversity and Richness')+ggtitle(title)+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
   
   list(Hn.plot,Hn.probability.plot,Hs.plot,Hs.probability.plot,joint.probability.plot,plot.data)
}

#### Simulate Haplotype Frequency ####
rinfall<-function(theta,N){
  #Theta is the estimate of theta from somewhere else where theta = 4Nu
  #N is the population size (not sampled but overall)
  CRP<-theta/(theta+(1:N)-1)
  new.haps<-runif(N) < CRP
  
  sim.haps<-rep(NA,N)
  sim.haps[new.haps]<-1:sum(new.haps)
  
  for(n in 1:N){
    if(is.na(sim.haps[n])){
      sim.haps[n]<-sim.haps[sample(n-1,1)]
    }
  }
  
  unname(table(sim.haps)[order(table(sim.haps),decreasing=T)])
}

#### Add in -9 to make arrays all same size ####
balloon<-function(sim.res,dimension){
  supplement.dimensions<-dimension-dim(sim.res)
  if(!all(supplement.dimensions==0)){
    supplement.dimensions<-ifelse(supplement.dimensions==0,dimension,supplement.dimensions)
    
    supplement<-array(-9,dim=supplement.dimensions)
    out<-abind(sim.res,supplement,along=1)
  }
  if(all(supplement.dimensions==0)){
    out<-sim.res
  }
  out
}

#Summarize model sum to frequency
freq.summary<-function(x){
  x/rowSums(x)
}

month.thinning<-function(x,RUN.MONTH){
  #x<-x[c(1,seq(6,RUN.MONTH,by=6)),]
  x<-x[c(1,RUN.MONTH),]
  x
}

remove.0.haps<-function(x){
  all(x==0)
}
