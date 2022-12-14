library(ggplot2)
library(sjPlot)
library(dplyr)
library(ggiraphExtra)
library(corrplot)
library(GGally)
library(ggpubr)
library(hrbrthemes)
library(huxtable)
library(jtools)
df <- read.csv("Q:\\courseProject\\R\\plot_stats_50_NNA.csv",header = T)
df <- read.csv("/Users/Lsy/Library/CloudStorage/Box-Box/Duke/fall 2022/RS_558/final/R/plot_stats_50_NNA.csv",header = T)

#add "Northness"
df <- df %>%
  mutate(
    asp_N = case_when(
      between(aspect_mean, 0, 22.5) ~ 1,
      between(aspect_mean, 22.5, 67.5) ~ 0.5,
      between(aspect_mean, 67.5, 112.5) ~ 0,
      between(aspect_mean, 112.5, 157.5) ~ -0.5,
      between(aspect_mean, 157.5, 202.5) ~ -1,
      between(aspect_mean, 202.5, 247.5) ~ -0.5,
      between(aspect_mean, 247.5, 292.5) ~ 0,
      between(aspect_mean, 292.5, 337.5) ~ 0.5,
      between(aspect_mean, 337.5, 360) ~ 1,
    )
  )
#save file
write.csv(df,"Q:\\courseProject\\R\\plot_stats_50_aspN.csv")
corrplot(cor(df[,-1]))

#ggpairs(data.frame(t(df)))

#observe distribution
ggplot(df, aes(x = asp_N) ) +
  geom_point(aes(y = TreeCount),color='red',alpha=.15) +
  #geom_smooth(y = df$TreeCount,method = "glm",alpha = .15)+
  geom_point(aes(x = asp_N, y = TreeH_mean),color='blue',alpha=.15) +
  scale_y_continuous(
    # Features of the first axis
    name = "Tree Count",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*1, name="Crown Area SD")
  )+ 
  
  theme_ipsum() +
  theme(
    axis.title.y = element_text(color = 'red', size=13),
    axis.title.y.right = element_text(color = 'blue', size=13)
  )  
  #scale_color_gradient2(midpoint=mean(df$slope_mean), low="red", mid="white",
  #                      high="blue", space ="Lab" )
ggplot(df, aes(x = elev_mean, y = CrownArea_SD, color = asp_N) ) +
  geom_point() +
  geom_smooth(method = "glm",alpha = .15, aes(fill = asp_N))+
  scale_color_gradient2(midpoint=mean(df$asp_N), low="red", mid="white",
                      high="blue", space ="Lab" )
ggplot(df, aes(x = elev_mean, y = TreeH_mean, color = asp_N) ) +
  geom_point() +
  scale_color_gradient2(midpoint=mean(df$asp_N), low="red", mid="white",
                      high="blue", space ="Lab" )
#build glms
####for treecount
par(mfrow=c(2,2))
plot(glm.tc.full)
ggqqplot(df$TreeCount)
ggqqplot(df$CrownArea_mean)
ggqqplot(df$TreeH_mean)
glm.tc.full <- glm(df, formula = TreeCount ~ elev_mean + slope_mean + asp_N + elev_mean:slope_mean + slope_mean:asp_N + elev_mean:asp_N, family = poisson())
summary(glm.tc.full)
step(glm.tc.full)
#plot
par(mfrow=c(2,2))
plot(glm.tc.full)
autoplot(glm.tc.full,label.size = 3)
####split slope at 40 degree
df_steep <- df[df$slope_mean>40,]
df_gentle <- df[df$slope_mean<=40,]
#fit 
glm.tc.steep <- glm(df_steep, formula = TreeCount ~ elev_mean + slope_mean + asp_N, family = poisson()) #no elev_mean:slope_mean
summary(glm.tc.steep)
#step(glm.tc.steep)
#plot
par(mfrow=c(2,2))
plot(glm.tc.steep)
#gentle
glm.tc.gentle <- glm(df_gentle, formula = TreeCount ~ elev_mean + slope_mean +asp_N , family = poisson())
summary(glm.tc.gentle)
#step(glm.tc.steep)
#plot
par(mfrow=c(2,2))
plot(glm.tc.gentle)
#plot splitted curve TC to slope
cfs_stp <- coef(glm.tc.steep)
cfs_gtl <- coef(glm.tc.gentle)
with(df, plot(x = slope_mean, y = TreeCount, pch = 4, xlab = "slope_mean", ylab = "TreeCount"))
x <- with(df, seq(min(slope_mean), max(slope_mean), length = 100)) 
curve(exp(cfs_gtl[1] + cfs_gtl[2]*mean(df_gentle$elev_mean) + cfs_gtl[3]*x + cfs_gtl[4]*mean(df_gentle$asp_N)), add = T, col = "blue", lwd = 2, xlim = c(0,40))
curve(exp(cfs_stp[1] + cfs_stp[2]*mean(df_steep$elev_mean) + cfs_stp[3]*x + cfs_stp[4]*mean(df_steep$asp_N)), add = T, col = "blue", lwd = 2, xlim = c(40,max(df$slope_mean)))



##for crwon area
glm.cr.full <- glm(df, formula = CrownArea_mean ~ elev_mean + slope_mean + asp_N + elev_mean:slope_mean + slope_mean:asp_N + elev_mean:asp_N, family = gaussian())
summary(glm.cr.full)
step(glm.cr.full)
#plot
par(mfrow=c(2,2))
plot(glm.cr.full)
autoplot(glm.cr.full,label.size = 3)

#for treeheight
glm.th.full <- glm(df, formula = TreeH_mean ~ I(elev_mean^2) + slope_mean + elev_mean:slope_mean + slope_mean:asp_N, family = gaussian())
summary(glm.th.full)
step(glm.th.full)
#plot
par(mfrow=c(2,2))
plot(glm.th.full)
autoplot(glm.th.full,label.size = 3)
#check collinearity  (less than 5, multicollinearity is low)
car::vif(glm.tc.full)
# check model fit
pchisq(glm.th.full$deviance, glm.th.full$df.residual, lower.tail=F)
#do not reject null, residual deviance is not significant different from chi-square distribution, model
pscl::pR2(glm.th.full)[4]
##ENV710 Lab9

###export model summaries

summ(glm.tc.full,digits = 4)
summ(glm.cr.full,digits = 4)
summ(glm.th.full,digits = 4)
