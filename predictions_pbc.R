library(joineR)
library(nlme)
library(lcmm)
library(pec)
library(timeROC)
library(dplyr)
library(JM)


data("pbc2")
pbc2$id=as.numeric(pbc2$id)
pbc2$status3[pbc2$status=="alive"]<-0
pbc2$status3[pbc2$status=="transplanted"]<-1
pbc2$status3[pbc2$status=="dead"]<-2
pbc2$drug2[pbc2$drug=="D-penicil"]<-1
pbc2$drug2[pbc2$drug=="placebo"]<-0

landmark.times=c(2,4,6,8)
window.t=5




jlcmFit5 <- Jointlcmm(
  fixed = serBilir ~ drug * year,
  random = ~ year,
  mixture = ~ year+ drug:year,
  subject = "id",
  survival = Surv(years, status3) ~ cause(drug),
  hazard = "Weibull",
  hazardtype = "PH",
  ng = 5,
  #partialH = FALSE,
  logscale = TRUE,  
  data = pbc2)

###create a dataset wth unique id and times
pbcn=distinct(pbc2,id,.keep_all=TRUE)
pbcn=pbcn[,c(1,2,21,22)]
AUCst.dose <- rep(NA,length(landmark.times))
sdAUCst.dose <- rep(NA,length(landmark.times))
matiidAUCst.dose <- matrix(NA,nrow=length(landmark.times),ncol=nrow(pbcn))
rownames(matiidAUCst.dose) <- (paste0("s=",landmark.times))

matrixstat <- matrix(NA,nrow=length(landmark.times),ncol=4)
colnames(matrixstat) <- c("Cases (s,s+t)","survivors (s,s+t)","OtherEvent (s,s+t)","Censored  (s,s+t)")
rownames(matrixstat) <- (paste0("s=",landmark.times))
ComTimeAUCst <- rep(NA,length(landmark.times))
CumIncst <- rep(NA,length(landmark.times))
Survst <- rep(NA,length(landmark.times))


BrierS.s.dose <- rep(NA,length(landmark.times))
se.BS.s.dose <- rep(NA,length(landmark.times))
Mat.iid.BS.dose<-  matrix(NA,nrow=length(landmark.times),ncol=nrow(pbcn))
for (s in landmark.times){
  start.s <- Sys.time()
  predictions=dynpred(jlcmFit5,landmark = s,horizon = window.t,
                      var.time = "year",newdata = pbc2,
                      event = 1)
  
  res=data.frame(predictions$pred)
  
  # Create landmark data set
  d.s <- pbcn[,c("years","status3")]
  d.s$Pred.s.dose <- res[,4]
  d.s <- d.s[d.s$years>s,]
  d.s$time.s <- d.s$years-s
  d.s=na.omit(d.s)
  # AUC and BS for prediction based on IST
  # estimate ROC curve and AUC
  ROC.s.dose <- timeROC(T=d.s$time.s,
                        delta=d.s$status3,
                        marker=d.s$Pred.s.dose,
                        cause=1,weighting="marginal",
                        times=c(window.t),
                        iid=TRUE)
  # pec(as.matrix(d.s$Pred.s.dose),Surv(with.time,status)~treat,data=epil,cause=1,exact = TRUE)
  BS.s.dose <- BS(timepoints=c(window.t),
                  times=d.s$time.s,
                  status=d.s$status3,
                  pred=as.matrix(d.s$Pred.s.dose),
                  cause=1)
  
  BrierS.s.dose[which(s==landmark.times)] <- BS.s.dose$BS # BS estimate
  se.BS.s.dose[which(s==landmark.times)] <- BS.s.dose$sd  # BS s.e. estimate
  AUCst.dose[which(s==landmark.times)] <- ROC.s.dose$AUC_2[2] # AUC estimate
  sdAUCst.dose[which(s==landmark.times)] <- ROC.s.dose$inference$vect_sd_2[2]  # AUC s.e. estimate
  matrixstat[which(s==landmark.times),] <- ROC.s.dose$Stats[2,] # proportions of cases, controls, censored subjects within the prediction window and survivors at s+t
  CumIncst[which(s==landmark.times)] <- ROC.s.dose$CumulativeIncidence[2]  # marginal cumulative incidence  = P(T<s+t,eta=1|T>s)
  Survst[which(s==landmark.times)] <- ROC.s.dose$survProb[2] # marginal survival probability = P(T>s+t|T>s)
  
}
#AUC plots
plot(landmark.times,AUCst.dose,type="l",xlab="Landmark time s(years)",ylab="estimates of AUC(s,t)",
     main="Time dependent AUC curve(Death)",col="black",lwd=2,ylim = c(0,1))
points(landmark.times,AUCst.dose,pch=20,lwd=2)

lines(landmark.times,AUCst.dose-1.96*sdAUCst.dose,col="grey60",lty=2,lwd=2)
lines(landmark.times,AUCst.dose+1.96*sdAUCst.dose,col="grey60",lty=2,lwd=2)

#Brier score plot
plot(landmark.times,BrierS.s.dose,type="l",xlab="Landmark Time s(years)",ylab="Estimates of BS(s,t)",
     main="Dynamic Brier Score",col="black",ylim=c(0.1,0.4),lwd=2)
points(landmark.times,BrierS.s.dose,pch=20,lwd=2)
# axis(side=1,at=c(0.5,1,1.5,2,2.5,3,3.5),labels=c("0.5","1","1.5","2","2.5","3","3.5"))
# axis(side=2,at=c(0.05),labels=c("0.05"))
lines(landmark.times,BrierS.s.dose-1.96*se.BS.s.dose,col="black",lty=2,lwd=2)
lines(landmark.times,BrierS.s.dose+1.96*se.BS.s.dose,col="black",lty=2,lwd=2)

