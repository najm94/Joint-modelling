library(JM)
library(lcmm)
data("pbc2")
pbc2$id=as.numeric(pbc2$id)
pbc2$status3[pbc2$status=="alive"]<-0
pbc2$status3[pbc2$status=="transplanted"]<-1
pbc2$status3[pbc2$status=="dead"]<-2
jlcmFit1 <- Jointlcmm(
  fixed = serBilir ~ drug * year,
  random = ~ year,
  subject = "id",
  survival = Surv(years, status3) ~ cause(drug),
  hazard = "Weibull",
  hazardtype = "PH",
  ng = 1,
  logscale = TRUE,
  data = pbc2)
summary(jlcmFit1)
jlcmFit2 <- Jointlcmm(
  fixed = serBilir ~ drug * year,
  random = ~ year,
  mixture = ~ year+ drug:year,
  subject = "id",
  survival = Surv(years, status3) ~ cause(drug),
  hazard = "Weibull",
  hazardtype = "PH",
  ng = 2,
  #partialH = FALSE,
  logscale = TRUE,  
  data = pbc2)
summary(jlcmFit2)
jlcmFit3 <- Jointlcmm(
  fixed = serBilir ~ drug * year,
  random = ~ year,
  mixture = ~ year+ drug:year,
  subject = "id",
  survival = Surv(years, status3) ~ cause(drug),
  hazard = "Weibull",
  hazardtype = "PH",
  ng = 3,
  #partialH = FALSE,
  logscale = TRUE,  
  data = pbc2)
summary(jlcmFit3)
confint(jlcmFit3)
?Jointlcmm
jlcmFit4 <- Jointlcmm(
  fixed = serBilir ~ drug * year,
  random = ~ year,
  mixture = ~ year+ drug:year,
  subject = "id",
  survival = Surv(years, status3) ~ cause(drug),
  hazard = "Weibull",
  hazardtype = "PH",
  ng = 4,
  #partialH = FALSE,
  logscale = TRUE,  
  data = pbc2)
pbc2$serb1=log(pbc2$serBilir)
fit1=Jointlcmm(
  fixed = serb1 ~ drug * year,
  random = ~ year,
  mixture = ~ year+ drug:year,
  subject = "id",
  survival = Surv(years, status3) ~ cause(drug),
  hazard = "Weibull",
  hazardtype = "PH",
  ng = 4,
  #partialH = FALSE,
  logscale = TRUE,  
  data = pbc2)
plot(fit1)
summary(fit1)
summary(jlcmFit4)
plot(jlcmFit3)
confint(jlcmFit4)
jlcmFit5 <- Jointlcmm(
  fixed = serb1 ~ drug * year,
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
plot(jlcmFit5)
predictions=dynpred(jlcmFit5,landmark = c(2,4,6,8),horizon = 5,var.time = "year",newdata =pbc2[3:11,],
                    event = 1)
predictions2=dynpred(jlcmFit5,landmark = c(2,4,6,8),horizon = 5,var.time = "year",newdata =pbc2[3:11,],
                    event = 2)
plot(predictions)
plot(predictions2)
summary(jlcmFit5)
options(scipen=999)
?confint

confint(jlcmFit5)
jlcmFit6 <- Jointlcmm(
  fixed = serBilir ~ drug * year,
  random = ~ year,
  mixture = ~ year+ drug:year,
  subject = "id",
  survival = Surv(years, status3) ~ cause(drug),
  hazard = "Weibull",
  hazardtype = "PH",
  ng = 6,
  #partialH = FALSE,
  logscale = TRUE,  
  data = pbc2)
summary(jlcmFit6)
jlcmFit4$scoretest
jlcmFit5.gs <- gridsearch(
  m = Jointlcmm(
    fixed = serBilir ~ drug * year,
    random = ~ year,
    mixture = ~ year + drug:year,
    subject = "id",
    survival = Surv(years, status3) ~ cause(drug),
    hazard = "Weibull",
    hazardtype = "PH",
    ng = 5,
    logscale = TRUE,    
    data = pbc2),
  rep = 40,
  maxiter = 35,
  minit = jlcmFit1
)
summary(jlcmFit5.gs)
jlcmFit5.gs$scoretest

summarytable(jlcmFit1,jlcmFit2,jlcmFit3,jlcmFit4,jlcmFit5)
plot(jlcmFit5)

dev.off()
summary(pbc2$year)
plot(cuminc(jlcmFit5, time = seq(0, 14, 2), drug = c(0, 1)), profil = 1, event = 1,
     main = "", bty = "n", ylim = c(0, 1), xlim = c(0, 14), lwd = 2,
     legend.loc = "topright",
     xlab = "Time from randomization (years)",
     ylab = "Cumulative incidence of transplantation")
plot(cuminc(jlcmFit5, time = seq(0, 14, 2), drug = c(0, 1)), profil = 1, event = 2,
     main = "", bty = "n", ylim = c(0, 1), xlim = c(0, 14), lwd = 2,
     legend.loc = "topright",
     xlab = "Time from randomization (years)",
     ylab = "Cumulative incidence of death")