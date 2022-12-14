library(corrplot)
library(dplyr)
library(ggplot2)

#read in table
plot_stats <- read.csv("Q:\\courseProject\\R\\stats.csv", sep = " ",header = T)
#format digits
plot_stats<-lapply(plot_stats,round,4)
#transpose
stats_tr<- data.frame(t(plot_stats))
#correlation 
png(filename="Q:\\courseProject\\R\\corr.png")
corrplot(cor(stats_tr))
dev.off()

#categorize aspect
stats_tr$asp_alph <- 1:length(plot_stats[1,])
stats_tr <- stats_tr %>%
  mutate(
    asp_alph = case_when(
      between(aspect_mean, 0, 22.5) ~ "N",
      between(aspect_mean, 22.5, 67.5) ~ "NE",
      between(aspect_mean, 67.5, 112.5) ~ "E",
      between(aspect_mean, 112.5, 157.5) ~ "SE",
      between(aspect_mean, 157.5, 202.5) ~ "S",
      between(aspect_mean, 202.5, 247.5) ~ "SW",
      between(aspect_mean, 247.5, 292.5) ~ "W",
      between(aspect_mean, 292.5, 337.5) ~ "NW",
      between(aspect_mean, 337.5, 360) ~ "N",
    )
  )
source("jackLM.R")
rich.jack <- jackLM(k.dat[,8],k.dat[,1:7])
print(rich.jack$jtable, digits=3)
#build uni-variate models; dR2without is r^2 that lost if delete that variable
# a parsimonious model?
rich.step <- step(rich.lm)
summary(rich.step)

# end normal regression.

# prep for SEM ...
# fix the scale of the data ...
# lavaan wants the variables to scale to the same order of magnitude:
# z score is in a same unit (so comparable) and unique to certain data, not migratable.
summary(k.dat)
k.dat$distance <- k.dat$distance/100
k.dat$elev <- k.dat$elev/1000
k.dat$abiotic <- k.dat$abiotic/100
k.dat$age <- k.dat$age/100
k.dat$firesev <- k.dat$firesev/10
k.dat$rich <- k.dat$rich/100

# load library:
library(lavaan)

model2 <- '
  hetero ~ a1*distance
  abiotic ~ b1*distance
  age ~ c1*distance
  firesev ~ c2*age
  cover ~ c3*firesev
  rich ~ a2*hetero
  rich ~ b2*abiotic
  rich ~ c4*cover
  rich ~ d1*distance
  # the effects (will collect after model fitted):
  direct:= d1
  ind.het:= a1*a2 
  ind.abio:= b1*b2 
  ind.fire:= c1*c2*c3*c4
  total:= direct + ind.het + ind.abio + ind.fire
'
# refit
model2.sem <- sem(model2, data=k.dat)
summary(model2.sem,rsq=T)
(mi3 <- modindices(model3.sem))
anova(model2.sem,model3,sem,test="Chi")
# look at the standardized solutions/parameters
model2.std <- standardizedSolution(model2.sem)
model2.std
library(semPlot)
# play with the plot parameters ...
semPaths(model2.sem, what='std',style='lisrel',layout='tree2',rotation=1,sizeMan=8,residuals=F)

