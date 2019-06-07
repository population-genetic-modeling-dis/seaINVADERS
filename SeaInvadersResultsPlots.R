rm(list=ls())
library(ggplot2)
library(Rmisc)

setwd(dir = "Masters/DIS/Results/")
SeaInvaders <- read.table("jobs_to_run_ForR2.csv", sep=",", header=TRUE)


 


###########################################################################################################################


Pvol_Atl_Carib <- SeaInvaders[SeaInvaders$Model==1,] 
ggplot(data=Pvol_Atl_Carib,aes(x=MostProbColonists,y=RPI,colour=as.factor(Theta))) +
  geom_errorbar(aes(ymin=RPI-Lower50CI, ymax=RPI+Upper50CI,width=0)) +
  ggtitle(label = "Pvol_Atl_Carib_RPIvMP") +
  theme(legend.text=element_text(colour="#000000",size=18)) + 
  theme(legend.title=element_text(colour="#000000",size=18)) +
  theme(axis.title.x=element_text(face="bold",colour="#000000",size=24)) + 
  theme(axis.title.y=element_text(face="bold",colour="#000000",size=24)) +
  theme(axis.text.y=element_text(face="bold",colour="#000000",size=12)) + 
  theme(axis.text.x=element_text(face="bold",colour="#000000",size=12)) + 
  geom_point() 


###########################################################################################################################


Pvol_Indo_Atl <- SeaInvaders[SeaInvaders$Model==2,] 
ggplot(data=Pvol_Indo_Atl,aes(x=MostProbColonists,y=RPI,colour=as.factor(Theta))) +
  geom_errorbar(aes(ymin=RPI-Lower50CI, ymax=RPI+Upper50CI,width=0)) + 
  ggtitle(label = "Pvol_Indo_Carib_RPIvMP") +
  theme(legend.text=element_text(colour="#000000",size=18)) + 
  theme(legend.title=element_text(colour="#000000",size=18)) +
  theme(axis.title.x=element_text(face="bold",colour="#000000",size=24)) + 
  theme(axis.title.y=element_text(face="bold",colour="#000000",size=24)) +
  theme(axis.text.y=element_text(face="bold",colour="#000000",size=12)) + 
  theme(axis.text.x=element_text(face="bold",colour="#000000",size=12)) + 
  geom_point() 


###########################################################################################################################


Lkas_Marq_Haw <- SeaInvaders[SeaInvaders$Model==3,] 
ggplot(data=Lkas_Marq_Haw,aes(x=RPI,y=MostProbColonists,colour=as.factor(Theta))) +
  geom_errorbar(aes(ymin=MostProbColonists-Lower50CI, ymax=MostProbColonists+Upper50CI),width=0) +  
  ggtitle(label = "Lkas_Marq_Haw_RPIvMP") +
  theme(legend.text=element_text(colour="#000000",size=18)) + 
  theme(legend.title=element_text(colour="#000000",size=18)) +
  theme(axis.title.x=element_text(face="bold",colour="#000000",size=24)) + 
  theme(axis.title.y=element_text(face="bold",colour="#000000",size=24)) +
  theme(axis.text.y=element_text(face="bold",colour="#000000",size=12)) + 
  theme(axis.text.x=element_text(face="bold",colour="#000000",size=12)) + 
  geom_point() +
  stat_smooth(method = 'nls', formula = 'y~a*x^b' , method.args = list(start= c(a = 200,b=-0.3)),se=FALSE)


###########################################################################################################################


Lkas_Soc_Haw <- SeaInvaders[SeaInvaders$Model==4,] 
ggplot(data=Lkas_Soc_Haw,aes(x=RPI,y=MostProbColonists,colour=as.factor(Theta))) +
  geom_errorbar(aes(ymin=MostProbColonists-Lower50CI, ymax=MostProbColonists+Upper50CI),width=0) +  
  ggtitle(label = "Lkas_Soc_Haw_RPIvMP") +
  theme(legend.text=element_text(colour="#000000",size=18)) + 
  theme(legend.title=element_text(colour="#000000",size=18)) +
  theme(axis.title.x=element_text(face="bold",colour="#000000",size=24)) + 
  theme(axis.title.y=element_text(face="bold",colour="#000000",size=24)) +
  theme(axis.text.y=element_text(face="bold",colour="#000000",size=12)) + 
  theme(axis.text.x=element_text(face="bold",colour="#000000",size=12)) + 
  geom_point() +
  stat_smooth(method = 'nls', formula = 'y~a*x^b',  method.args = list(start= c(a = 200,b=-0.3)),se=FALSE)
