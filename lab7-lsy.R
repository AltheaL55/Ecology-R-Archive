library(lavaan)
library(AICcmodavg)
setwd("D:/ELsy/2020fall/EDA/7/Labnotes")

#final model in the class handout
#read data
streamdata <- read.csv("Streamproductivity.csv", header = T)
#preprocessing data
streamdata$Patch <- streamdata$PatchDiversity/10
streamdata$Stream <- streamdata$StreamDiversity/100
streamdata$O2 <- streamdata$O2Production*100
#fit original model
model1 <- 'Patch ~ Stream + logNutrient + logNutrient2
            Biomass ~ Patch + logNutrient
            logNutrient ~~ logNutrient2
            O2 ~ Biomass'
fit1 <- sem(model1, data = streamdata)
summary(fit1)
subset(modindices(fit1), mi > 3.8)

#delete lognutrient for patch(p>0.05), add none(no variance's mi is greater than 3.8)
model2 <- 'Patch ~ Stream + logNutrient2
            Biomass ~ Patch + logNutrient
            logNutrient ~~ logNutrient2
            O2 ~ Biomass'
fit2 <- sem(model2, data = streamdata)
summary(fit2)
subset(modindices(fit2), mi > 3.8)
anova(fit2,fit1)
#fit2 is slightly better than fit1(lower AIC) and could be regard as a satisfying model.
#According to anova test, the two shows no significant difference.


#Fit model for forest fire and biodiversity
#read data
dat <- read.csv("LCMdata.csv", header = T)
#preprocessing data
#ppt94-98 = annual precipitation; covtot94-98 = herb cover; sr94-98 = species richness; 
#fireint1 = fire severity; abio = abiotic conditions
data2 <- with(dat, data.frame(age))
data2$r1 <- dat$sr94adj
data2$r2 <- dat$sr95adj
data2$r3 <- dat$sr96adj
data2$r4 <- dat$sr97adj
data2$r5 <- dat$sr98adj
data2$c1 <- dat$covtot94
data2$c2 <- dat$covtot95
data2$c3 <- dat$covtot96
data2$c4 <- dat$covtot97
data2$c5 <- dat$covtot98
data2$abio <- dat$abio
data2$fire <- dat$fireint1

#acoording to Fig 5.10 & 5.11 in the original article, r1¨Cr5 with time-varying covariates and c1-c5, as well as time-invariant covariate, abio (abiotic conditions) are included in the model.
fmodel1 <- '
#intercept and slope with fixed coefficients
i =~ 1*r1 +1*r2 +1*r3 +1*r4 +1*r5
s =~ 0*r1 +1*r2 +2*r3 +3*r4 +4*r5
# time-varying regressions
r1 ~ c1
r2 ~ c2
r3 ~ c3
r4 ~ c4
r5 ~ c5
# autoregressive effects
r2 ~ r1
# time-invariant effects
i ~ abio
s ~ abio
# fire severity effects
r1 ~ fire +abio
c1 ~ fire +abio
fire ~ age +abio'

fit.fm1 <- sem(fmodel1, data=data2)
print(fit.fm1)
summary(fit.fm1, rsq=TRUE)
subset(modindices(fit.fm1), mi > 3.8)


fmodel2 <- '
#intercept and slope with fixed coefficients
i =~ 1*r1 +1*r2 +1*r3 +1*r4 +1*r5
s =~ 0*r1 +1*r2 +2*r3 +3*r4 +4*r5
r1 ~ c1
r2 ~ c2
r3 ~ c3
r2 ~ r1
i ~ abio
r1 ~ fire
c1 ~ fire
fire ~ age
r2 ~~  c1
'

fit.fm2 <- sem(fmodel2, data=data2)
print(fit.fm2)
summary(fit.fm2, rsq=TRUE)
subset(modindices(fit.fm2), mi > 3.8)

#add i  ~  c2; c1  ~  c2;
fmodel3 <- '
i =~ 1*r1 +1*r2 +1*r3 +1*r4 +1*r5
s =~ 0*r1 +1*r2 +2*r3 +3*r4 +4*r5
r1 ~ c1
r2 ~ c2
r3 ~ c3
r2 ~ r1
i ~ abio
r1 ~ fire
c1 ~ fire
fire ~ age
r2 ~~  c1
c1  ~  c2
i  ~  c2'
fit.fm3 <- growth(fmodel3, data=data2)
print(fit.fm3)
summary(fit.fm3, rsq=TRUE)
subset(modindices(fit.fm3),mi > 3.8)

#add i ~ r1; delete r2 ~ c2, c1 ~ fire.
fmodel4 <- '
i =~ 1*r1 +1*r2 +1*r3 +1*r4 +1*r5
s =~ 0*r1 +1*r2 +2*r3 +3*r4 +4*r5
r1 ~ c1
r3 ~ c3
r2 ~ r1
i ~ abio
r1 ~ fire
fire ~ age
r2 ~~  c1
c1  ~  c2
i  ~  c2
i ~ r1'
fit.fm4 <- growth(fmodel4, data=data2)
print(fit.fm4)
summary(fit.fm4, rsq=TRUE)
subset(modindices(fit.fm4), mi > 3.8)
anova(fit.fm3,fit.fm4)

#add c1 ~ r1; c1 ~ r1; c1 ~ abio;c2 ~ c1
fmodel5 <- '
i =~ 1*r1 +1*r2 +1*r3 +1*r4 +1*r5
s =~ 0*r1 +1*r2 +2*r3 +3*r4 +4*r5
r1 ~ c1
r3 ~ c3
r2 ~ r1
i ~ abio
r1 ~ fire
fire ~ age
r2 ~~  c1
c1  ~  c2
i  ~  c2
i ~ r1
c1 ~ r1
c1 ~ abio
c2 ~ c1'
fit.fm5 <- growth(fmodel5, data=data2)
print(fit.fm5)
summary(fit.fm5, rsq=TRUE)
subset(modindices(fit.fm5), mi > 3.8)
anova(fit.fm5,fit.fm4)
#delete c1 ~ r1

#include p(mean of precipitation) into model
fmodel6 <- '
i =~ 1*r1 +1*r2 +1*r3 +1*r4 +1*r5
p =~ 0*r1 +1*r2 +0*r4 +2*r5
s =~ 0*r1 +1*r2 +2*r3 +3*r4 +4*r5
p ~ r3
r1 ~ c1
r2 ~ c2
r3 ~ c3
r4 ~ c4
r5 ~ c5
r3 ~ r2
r1 ~ fire
c1 ~ fire
i ~ abio
fire ~ age
r2 ~~  c1
i  ~  c2
i ~ r1
c1 ~ abio
c2 ~ c1'
fit.fm6 <- growth(fmodel6, data=data2)
print(fit.fm6)
summary(fit.fm6, rsq=TRUE)
subset(modindices(fit.fm6), mi > 3.8)
anova(fit.fm5,fit.fm6)

#Then the P-value (Chi-square) keep showing 0.000 in these following models.
#Although the p-value of latent variables are lower then 0.05, there seems to be too manny potential regression equations.
