library(joineR)
data(epileptic)
epileptic$interaction <- with(epileptic, time * (treat == "LTG"))
longitudinal <- epileptic[, c(1:3, 13)]
survival <- UniqueVariables(epileptic, c(4, 6), "id")
baseline <- UniqueVariables(epileptic, "treat", "id")
data <- jointdata(longitudinal = longitudinal,
                  survival = survival,
                  baseline = baseline,
                  id.col = "id",
                  time.col = "time")
fit2 <- joint(data = data,
              long.formula = dose ~ time + treat + interaction,
              surv.formula = Surv(with.time, with.status2) ~ treat,
              longsep = FALSE, survsep = FALSE,
              gpt = 3)

?joint
summary(fit2)

data("pbc2")
pbc2$drugF[pbc2$drug=="D-penicil"]<-1

pbc2$drugF[pbc2$drug=="placebo"]<-0
pbc2$status3[pbc2$status=="alive"]<-0
pbc2$status3[pbc2$status=="transplanted"]<-1
pbc2$status3[pbc2$status=="dead"]<-2
pbc2$interaction=with(pbc2,year*(drug=="D-penicil"))
longitpbc=pbc2[,c(1,7,22,12)]
survpbc=UniqueVariables(pbc2,c(2,23),"id")
baselinepbc=UniqueVariables(pbc2,"drug","id")
data <- jointdata(longitudinal = longitpbc,
                  survival = survpbc,
                  baseline = baselinepbc,
                  id.col = "id",
                  time.col = "year")
fit2 <- joint(data = data,
              long.formula = serBilir ~ year + drug+ interaction,
              surv.formula = Surv(years, status3) ~ drug,
              longsep = FALSE, survsep = FALSE,
              gpt = 4,max.it = 1000)
summary(fit2)
