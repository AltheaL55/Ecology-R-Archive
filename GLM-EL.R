setwd("/Users/Lsy/Box Sync/Duke/spring2022/724-LandAnalys/labs/birds")
#get data
birds.data <- read.csv("birds.csv")
habitat.data <- read.csv("habitat.csv")

#remove the bad samples #198 and 207
habitat.data <- habitat.data[c(-198, -207),]
##!!If deleting rows several seperate times, remember that the column number would change!!
#make the bird data match:
birds.data <- birds.data[c(-198, -207),]

#check correlations among habitat variables:
habitat.cor <- cor(habitat.data[,-1])
habitat.cor
#and test one correlation:
cor.test(habitat.data$grndcov,habitat.data$canpcov)
# cull out a focal species:
revi <- birds.data$revi
cor2m(habitat.data[,-1],as.matrix(revi))

#coef with bird; p-value sorted
source("spphab_cor.R")
shc <- spphab.cor(revi,habitat.data[,-1])
var.list <- as.numeric(rownames(shc))
tmp <- habitat.data[,-1]
revi.habitat <- tmp[,c(var.list)]
names(revi.habitat) # to see them in order

# check for redundant variables, with default correlation |0.7|
source("screen_cor.R")
screen.cor(revi.habitat)
# remove the redundant ones, favoring the ones with higher correlations
# with revi:
names(revi.habitat)
##delete the ones with |r|>0.7
revi.habitat <- revi.habitat[,-c(2,3,4,5,6,7,13)]
names(revi.habitat)
screen.cor(revi.habitat)
dim(revi.habitat)
# end of EDA
rhcn<- spphab.cor(revi,revi.habitat)
rhcn

###GLM
# cull out a focal species: American Robin

# now we merge these for convenience:
revi.data <- cbind(revi, revi.habitat)

# WORKFLOW:
# 1. prep data
# 2. fit model
# 3. evaluate model:
# 	a. significant? (P-value)
#	b. explanatory power?
#	c. variable importance
#	d. classification success
#		tuning/calibration?
# 4. (optional) validate on independent data
# 5. apply:  predict/forecast other scenarios (management, restoration, climate change ...)

# fit the GLM (this is a generic/default logistic GLM; the link defaults to logistic):
# revi.glm <- glm(as.factor(revi)~., data=revi.data, family=binomial)
# BUT
# since the sample sizes (0/1) are quite unbalanced, 
# refit with weighted obs, to make the effective sample sizes equal ...
# weight is (# presences/# absences) for absence points
# (putting the command in parens prints the result as it runs):
(p <- length(revi[revi==1]))  	# presences 
(a <- length(revi[revi==0]))  	# absences
(w <- p/a)			# assign relative weight of each absence
revi.wts <- rep(1,length(revi))
revi.wts[revi==0] <- w		# downweight the absences
# fit:
revi.glm <- glm(as.factor(revi)~., data=revi.data, family=binomial,weights=revi.wts)

# print the fitted model (i.e., the coefficients):
print(revi.glm)
# a concise summary, with tests of the linear predictors (t-tests)
summary(revi.glm)
# this is where you check the ecological plausibility of the model:
# do the significant predictors make sense? do they have the expected sign?

# ANOVA test of the linked model (which variables are significant?):
anova(revi.glm, test="Chi")

# an intuitive plot, with the fitted model and actual data:
plot(revi.glm$linear.predictor, revi.glm$fitted.values, ylim=c(0,1))
points(revi.glm$linear.predictor, revi)
# a prettier plot ...
plot(revi.glm$linear.predictor, revi.glm$fitted.values, ylim=c(0,1),pch=19,col="blue", xlab="Linear Predictor", ylab="P(habitat)")
revi.col <- rep("red",length(revi)) # a vector of 270 "red" values
revi.col[revi==1] <- "green" # turn the presences green
points(revi.glm$linear.predictor, revi, pch=19, col=revi.col)
legend("topleft",legend=c("absences","presences"),text.col=c("red","green"), bty="n",horiz=F)
box(lwd=2)
# a nice trick for better graphics, if you need to reuse them in PPT or elsewhere. 
# the jpeg() function allows you to size the plot (in inches) and specify the resolution;
# the resolution is important cuz most pubs won't accept low-res plots (300 is good; 600, better).
jpeg(filename="revi_glm.jpg", units="in", width=6, height=5, res=600)
# then re-run the plot/points/legend/box commands above, then 
# dev.off() finishes the plot and writes the JPG to your workspace. 
dev.off() 
# there are similar functions for PNGs and TIFFs, if you prefer.

# is it significant?
# a P value, based on null (initial) and residual deviance, Chi-square,
# and the degrees of freedom (# of estimated coefficients)...
# DF is # of predictors + 1 (for the intercept term):
nc.revi <- length(revi.glm$coefficients) # count them
revi.glm.p <- 1 - pchisq(revi.glm$null - revi.glm$deviance, nc.revi)   
#"null" is the full deviance, "deviance" is the residual deviance
revi.glm.p # print the P-value (for the whole model)

# approximate explanatory value (as proportion of deviance):
revi.glm.d <- 1 - (revi.glm$deviance/revi.glm$null)
revi.glm.d # print it (how much explained)
# this begs a lot of questions about how to properly estimate a
# pseudo-R2 for a GLM or other non-least-squares regression. 

# build and compare 2 nested models:
revi.glm1 <- glm(as.factor(revi)~rd.2000, data=revi.data, family=binomial)
revi.glm2 <- glm(as.factor(revi)~rd.2000+psize, data=revi.data, family=binomial) 
anova(revi.glm1, revi.glm2, test="Chi")
# we would like to see a significant change with the addition of the second predictor. 

# the (marginal) effect of a single predictor ... 
# library(sjPlot) # this package loads a million dependencies!
library(ggplot2)
# plot_model(revi.glm, type="pred",terms="rd.500")
# there are other ways to get partial response plots but they are a pain.
# a dirt-ball, approximate solution:
# plot the predictor(s) against the fitted values from the GLM...
# if it looks nonlinear, try a quadratic in the GLM, or a GAM (or other nonlinear model).
# e.g.,
plot(revi.glm$fitted.values, revi.data$rd.2000, pch=19)
# e.g., a GLM with a quadratic on evergreenness:
# revi.glmq <- glm(as.factor(revi)~poly(rd.500,2)+rd1.500+d.rivers, data=revi.data, family=binomial, weights=revi.wts)
#summary(revi.glmq)
##Noticed nonlinear, use polynomial concave? create independent polynomial term...Not intuitive to interpret the whole transformation
# the coefficients are reported for "poly(rd.500,2)1" (linear) and "poly(rd.500,2)2" (quadratic);
# here, the quadratic term is not significant (as we might have guessed from the previous plot)

# if you don't really care about the variables, and just want a parsimonious stepwise model, 
# post-process stepwise:
revi.glm.step <- step(revi.glm)
# summary of the final model:
summary(revi.glm.step)
anova(revi.glm.step,test="Chi")
#explanatory power:
1 - (revi.glm.step$deviance/revi.glm.step$null)

plot(revi.glm.step$linear.predictor, revi.glm.step$fitted.values, ylim=c(0,1),pch=19,col="blue", xlab="Linear Predictor", ylab="P(habitat)")
revi.col <- rep("red",length(revi)) # a vector of 270 "red" values
revi.col[revi==1] <- "green" # turn the presences green
points(revi.glm.step$linear.predictor, revi, pch=19, col=revi.col)
legend("topleft",legend=c("absences","presences"),text.col=c("red","green"), bty="n",horiz=F)
box(lwd=2)


# for more information on variables:  a jackknifing estimate of variable importance ...
# this is borrowed from the maxent software.
source("jack_wGLM.R")
# this function expects weights for the observations.
revi.jtable <- jackGLM(revi,revi.habitat,revi.wts)
revi.jtable
# easier to scan:
round(revi.jtable$jtable,3) # 3 significant digits

# to dump this from R to use in a report:
write.csv(revi.jtable$jtable, "revi_jtable.csv")

# model evaluation in terms of classification success ...
# to do this, we need to convert the prediction to binary, by thresholding. 
# to motivate this, consider thresholding to 0/1 in various ways ...

# threshold to a binary prediction, here at p(hab)=0.5 (which is arbitrary!):
revi.glm.p50 <- revi.glm$fitted.value
revi.glm.p50[revi.glm.p50<0.50] <- 0
revi.glm.p50[revi.glm.p50>=0.50] <- 1

# a simple confusion matrix, cut-off at p=0.5:
revi.cm.5 <- table(revi.glm.p50,revi)
revi.cm.5
# percent correct is (sum of diagonals)/(sum over full table),
# but note asymmetry of success rates of presences and absences.
sum(diag(revi.cm.5))/sum(revi.cm.5)

# look at the frequency distributions of the data w/re: the
# fitted values of the GLM:
par(mfrow=c(2,1))
revi.glm.pred <- revi.glm$fitted.values
hist(revi.glm.pred[revi==0],ylim=c(0,40),col="red",main="Absences",xlab="P(habitat)")
hist(revi.glm.pred[revi==1],ylim=c(0,40),col="green",main="Presences",xlab="P(habitat)")
# the task is to choose a cut-off value (on the X axis) to give the
# best binary predictions
dev.off()

# an alternative version of same, and an easier plot to view:
plot(density(revi.glm.pred[revi==0]),lwd=2,col="red",xlim=c(0,1),ylim=c(0,2),xlab="P(habitat)",main="")
lines(density(revi.glm.pred[revi==1]),lwd=2,col="dark green")
legend("topleft",legend=c("absences","presences"),text.col=c("red","dark green"), bty="n",horiz=T)
box(lwd=2)

# tuning the model using ROC curves ...
# attach the library:
library(ROCR)
# create a "prediction" object (table of predicted & actual values):
revi.pred <- prediction(revi.glm$fitted.values, revi)
revi.pred.step <- prediction(revi.glm.step$fitted.values, revi)
# create a "performance" object with true positives, false positives:
revi.perf <- performance (revi.pred, "tpr", "fpr")
revi.perf.step <- performance (revi.pred.step, "tpr", "fpr")

# plot a ROC curve:
plot(revi.perf.step, colorize=TRUE,lwd=3)
abline(0,1, lty=2) # adds the 1:1 line
box(lwd=2)
##tangent point: get most cases right

# compute the area under the ROC curve (AUC):
revi.auc.step <- performance(revi.pred.step, "auc")
revi.auc.step@y.values[[1]] # I don't know why it's so hard to get the answer back!
##0.80 is great (great is subjective)

# plot sensitivity and specificity together:
plot(performance(revi.pred, "sens"),col="blue", lwd=2,ylab="Performance")
plot(performance(revi.pred, "spec"), add=TRUE, col="red",lwd=2)
legend("top", horiz=T,legend=c("Sensitivity","Specificity"), text.col=c("blue","red"),bty="n")
box(lwd=2) # make the bounding box thicker
##get intersect: balance sensitivity and specificity

# ROCR doesn't tell you what cut-off value to use ...
# to retrieve actual cut-off values, use our function 
# in script cutoff.ROCR.R
# either File->Open Script [browse to this file], or
source("cutoff.ROCR.R")
# the "help" is cutoff.ROCR.txt; open it in notepad to read it


# see how you did in tuning:
#max, maximizes the SUM of TPR and TNR
cutoff.max <- cutoff.ROCR(revi.pred, method="max")
revi.glm.max <- revi.glm$fitted.value
revi.glm.max[revi.glm.px<cutoff.max] <- 0
revi.glm.max[revi.glm.px>=cutoff.max] <- 1
revi.cm.max <- table(revi.glm.max,revi)
revi.cm.max
# percent correct:
sum(diag(revi.cm.max))/sum(revi.cm.max)

#(intersect of sensitivity and specificity) balances TPR and TNR
cutoff.int <- cutoff.ROCR(revi.pred, method="x")
revi.glm.px <- revi.glm$fitted.value
revi.glm.px[revi.glm.px<cutoff.int] <- 0
revi.glm.px[revi.glm.px>=cutoff.int] <- 1
revi.cm.int <- table(revi.glm.px,revi)
revi.cm.int
# percent correct:
sum(diag(revi.cm.int))/sum(revi.cm.int)

#stepwise intersect
cutoff.x <- cutoff.ROCR(revi.pred.step, method="x")
revi.glm.x <- revi.glm.step$fitted.value
revi.glm.x[revi.glm.x<cutoff.x] <- 0
revi.glm.x[revi.glm.x>=cutoff.x] <- 1
revi.cm.x <- table(revi.glm.x,revi)
revi.cm.x
# percent correct:
sum(diag(revi.cm.x))/sum(revi.cm.x)

# hit a target TPR (default is 95%):
cutoff.tpr <- cutoff.ROCR(revi.pred, method="tpr")
revi.glm.tpr <- revi.glm$fitted.value
revi.glm.tpr[revi.glm.tpr<cutoff.tpr] <- 0
revi.glm.tpr[revi.glm.tpr>=cutoff.tpr] <- 1
revi.cm.tpr <- table(revi.glm.tpr,revi)
revi.cm.tpr
# percent correct:
sum(diag(revi.cm.tpr))/sum(revi.cm.tpr)


cutoff.ROCR(revi.pred, "tpr", target=0.90)	# change the TPR target
#small size: maximize absent, true positive or balance "intersect" for tuning would better for ecological meaning 

#????note we didn't improve the classification of TPs at the cost of TNs.
# more typically, we want the % correct for presences and absences;
# these are tallied by column in the table (e.g., TPR = 51/(51+11)).
##turn into percentage and flip to "1 0" order when writing lab

# try a rough approximation for tuning without ROC ... the mean prediction:
(cut.glm <- mean(revi.glm$fitted.values))
revi.glm.px2 <- revi.glm$fitted.value
revi.glm.px2[revi.glm.px2<cut.glm] <- 0
revi.glm.px2[revi.glm.px2>=cut.glm] <- 1
table(revi.glm.px2,revi)
# the tuned version does the same for TPs, worse for TNs.
(70+125)/270

##pay more attention to TP (true absences in real data are often pseudo-absences (generated))
##"prevalence"
##"generalist"(average species, low niche position) "specialist"(distinct species, far from center, high niche position)

# end of GLM

citation() # returns citation for R itself
citation("ROCR") # returns citation for a package
citation("ecodist")

--------
##final!! removed weak!
# delete insignificant variables:
revi.glm.fin <- update(revi.glm.fin, .~.-h2o.x5-d.road)
# summary of the final model:
summary(revi.glm.fin)
anova(revi.glm.fin,test="Chi")
#explanatory power:
1 - (revi.glm.fin$deviance/revi.glm.fin$null)
#plot fitted value
plot(revi.glm.fin$linear.predictor, revi.glm.fin$fitted.values, ylim=c(0,1),pch=19,col="blue", xlab="Linear Predictor", ylab="P(habitat)")
revi.col <- rep("red",length(revi)) # a vector of 270 "red" values
revi.col[revi==1] <- "green" # turn the presences green
points(revi.glm.fin$linear.predictor, revi, pch=19, col=revi.col)
legend("topleft",legend=c("absences","presences"),text.col=c("red","green"), bty="n",horiz=F)
box(lwd=2)
# create a "prediction" object (table of predicted & actual values):
revi.pred.fin <- prediction(revi.glm.fin$fitted.values, revi)
# create a "performance" object with true positives, false positives:
revi.perf.fin <- performance (revi.pred.fin, "tpr", "fpr")

# plot a ROC curve:
plot(revi.perf.fin, colorize=TRUE,lwd=3)
abline(0,1, lty=2) # adds the 1:1 line
box(lwd=2)
# compute the area under the ROC curve (AUC):
revi.auc.fin <- performance(revi.pred.fin, "auc")
revi.auc.fin@y.values[[1]] # I don't know why it's so hard to get the answer back!

# plot sensitivity and specificity together:
plot(performance(revi.pred, "sens"),col="blue", lwd=2,ylab="Performance")
plot(performance(revi.pred, "spec"), add=TRUE, col="red",lwd=2)
legend("top", horiz=T,legend=c("Sensitivity","Specificity"), text.col=c("blue","red"),bty="n")
box(lwd=2) # make the bounding box thicker

#confusion matrix
#(intersect of sensitivity and specificity) balances TPR and TNR
cutoff.int <- cutoff.ROCR(revi.pred, method="x")
revi.glm.px <- revi.glm$fitted.value
revi.glm.px[revi.glm.px<cutoff.int] <- 0
revi.glm.px[revi.glm.px>=cutoff.int] <- 1
revi.cm.int <- table(revi.glm.px,revi)
revi.cm.int
# percent correct:
sum(diag(revi.cm.int))/sum(revi.cm.int)

