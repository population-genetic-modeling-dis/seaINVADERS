rm(list=ls())
setwd("C:/Users/cbird/OneDrive - Texas A&M University - Corpus Christi/PopGenSim2019/model_output_cbird/output_Lutjanus_kasmira_ctrl_Society_parameters_1556050188672")
lsdatas<-load("Sim.Society.theta.13.90878-rpr-4_source.RData")

# rm(list=ls())
# setwd("C:/Users/cbird/OneDrive - Texas A&M University - Corpus Christi/PopGenSim2019/model_output_cbird/output_Pterois_volitans_ctrl_parameters_1555982999780")
# lsdatas<-load("Sim.IndoP.theta.10.29332-rpr-0.25_source.RData")


dim(s)
tmp_s<-s[,,,1]
dim(tmp_s)
plot(apply(tmp_s[,,1], 1, sum))

tmp_s <- s[,,1,]
dim(tmp_s)
plot(apply(tmp_s[,,], 1, sum)/dim(tmp_s)[3])

pdf(file="testdemoplot.pdf", onefile=TRUE)
for(i in 1:30){
  tmp_s <- s[,,i,]
  dim(tmp_s)
  plot(apply(tmp_s[,,], 1, sum)/dim(tmp_s)[3])
}
dev.off()


library(ggplot2)

tmp_s<-s[,,,40]
tmp_s[,1:5,7]

tmp_s<-s[,,,25]
tmp_s[,1:5,7]

rbinom(1, 0, prob = -.1)
rmultinom(1,10,c(0,0,0))
