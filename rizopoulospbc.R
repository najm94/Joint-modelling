setwd("C:/Users/USER/Desktop/Semester4/thesis/MyCodes")
library(JM)
library(lme4)
library(nlme)
library(ggplot2)
library(dplyr)
data(pbc2)
summary(pbc2$year)
#pbc2$drug=as.numeric(pbc2$drug)
pbc2$drug=as.factor(pbc2$drug)
pbc2.idCR <- crLong(pbc2.id, statusVar = "status",
                    censLevel = "alive", nameStrata = "CR")
pbcn=distinct(pbc2,id,.keep_all=TRUE)
summary(pbcn)
table(pbcn$drug,pbcn$status)
prop.table(table(pbcn$drug,pbcn$status),1)
#profile plots
ggplot(aes(x = year, y = serBilir), data = pbc2) +
  geom_line(aes(group = id), colour = "grey", size = 0.8) +
  facet_grid(status ~ drug) +
  #theme_bw() +
  labs(
    x = "Time (years)", 
    y = "Serum Bilirubin(mg/dL)"
  ) +
  geom_smooth(aes(group = 1), colour = "red", se = FALSE, size = 1.5,
              method = "loess", span = 0.7) +
  scale_x_continuous(breaks=seq(0,14,2))+
  theme(text = element_text(size = 16))
ggplot(aes(x = year, y = serBilir), data = pbc2) +
  geom_line(aes(group = id), colour = "grey", size = 0.8) +
  facet_wrap(~status) +
  #theme_bw() +
  labs(
    x = "Time(years)", 
    y = "Serum Bilirubin(mg/dL)"
  ) +
  geom_smooth(aes(group = 1), colour = "red", se = FALSE, size = 1.5,
              method = "loess", span = 0.7) +scale_x_continuous(breaks=seq(0,14,2))+
  theme(text = element_text(size = 16))
#linear mixed model with time,treatment and interaction
lmeFit.pbc <- lme(log(serBilir) ~ drug * year,
                  random = ~ year | id, data = pbc2)
summary(lmeFit.pbc)
intervals(lmeFit.pbc)
#survival submodel
coxFit4.pbc <-
  coxph(Surv(years, status2) ~ drug  * CR + strata(CR),
        data = pbc2.idCR, x = TRUE)
summary(coxFit4.pbc)
confint(coxFit4.pbc)
# Current value parameterization
ptm <- proc.time()
jointFit1 <- jointModel(lmeFit.pbc, coxFit4.pbc,
                        timeVar = "year",
                        method = "spline-PH-aGH",
                        CompRisk = TRUE,
                        interFact = list(value = ~ CR, data = pbc2.idCR),
                        #GHk = 9,
                        #GKk = 15,
                        #iter.qN = 500,
                        #numeriDeriv = "cd",
                        verbose = FALSE)

summary(jointFit1)
confint(jointFit1)
plot(jointFit1)
#time dependent slopes
# Time-dependent slopes parameterization
ptm <- proc.time()
dform <- list(
  fixed = ~ 1 + drug, indFixed = 3:4,
  random = ~ 1, indRandom = 2)

jointFit2 <- jointModel(lmeFit.pbc, coxFit4.pbc,
                        timeVar = "year",
                        method = "spline-PH-aGH",
                        CompRisk = TRUE,
                        interFact = list(value = ~ CR, slope = ~ CR, 
                                         data =pbc2.idCR),
                        #GHk = 9,
                        #GKk = 15,
                        #iter.qN = 500,
                        #numeriDeriv = "cd",
                        parameterization = "both",
                        derivForm = dform,
                        verbose = FALSE)
summary(jointFit2)
confint(jointFit2)
anova(jointFit1,jointFit2,test = FALSE)
?anova
##lagged effects
jointFit3 <- update(jointFit1, lag = 0.25)
summary(jointFit3)
confint(jointFit3)
## Cummulative effects parameterization
ptm <- proc.time()
iform <- list(
  fixed = ~ -1 + year + I(drug * year) + 
    I(year^2 / 2) + I(drug * year^2 / 2), 
  indFixed = 1:4,
  random = ~ -1 + year + I(year^2 / 2), 
  indRandom = 1:2)

jointFit4 <- jointModel(lmeFit.pbc, coxFit4.pbc,
                        timeVar = "year",
                        method = "spline-PH-aGH",
                        CompRisk = TRUE,
                        interFact = list(slope = ~ CR, 
                                         data = pbc2.idCR),
                        #GHk = 9,
                        #GKk = 15,
                        #iter.qN = 500,
                        #numeriDeriv = "cd",
                        parameterization = "slope",
                        derivForm = iform,
                        verbose = FALSE)
summary(jointFit4)
confint(jointFit4)
str(pbc2)
?dnorm
##weighted cumulative
g <- function(u, pow = 0) {
  f <- function(t) integrate(function(s) s^pow * dnorm(t - s), 0, t)$value
  sapply(u, f)
}

iformW <- list(
  fixed = ~ -1 + I(g(year)) + I(drug * g(year)) + 
    I(g(year, 1)) + I(drug * g(year, 1)), 
  indFixed = 1:4,
  random = ~ -1 + I(g(year)) + I(g(year, 1)), 
  indRandom = 1:2)

jointFit5 <- jointModel(lmeFit.pbc, coxFit4.pbc,
                        timeVar = "year",
                        method = "spline-PH-aGH",
                        CompRisk = TRUE,
                        interFact = list(slope = ~ CR, 
                                         data = pbc2.idCR),
                        #GHk = 9,
                        #GKk = 15,
                        #iter.qN = 500,
                        #numeriDeriv = "cd",
                        parameterization = "slope",
                        derivForm = iformW,
                        verbose = FALSE)
summary(jointFit5)
confint(jointFit5)
proc.time() - ptm # CPU time

