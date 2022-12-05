{\rtf1\ansi\ansicpg1252\cocoartf2639
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fmodern\fcharset0 CourierNewPSMT;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\ri0\partightenfactor0

\f0\fs20 \cf0 setwd("/Users/Lsy/Box Sync/Duke/spring2022/724-LandAnalys/labs/sequoia")\
#import data\
spp17.data <- read.csv("spp17.csv")\
env21.data <- read.csv("env21.csv")\
xy.data <- read.csv("plotxy.csv")\
spp.summary <- read.csv("spp_summary.csv")\
#packages\
library(vegan) \
library(ecodist)\
library(corrplot)\
\
#remove rare species with Frequency < 5 plots\
spp.summary$Species[spp.summary$Frequency < 5]\
spp11.data <- spp17.data[, c(-4, -5, -7, -14, -17, -18)]\
names(spp11.data)\
\
#data relativizations, "Wisconsin double", by column maximum and the row sums\
source("reldata.R")\
spp11.rel <- reldata(spp11.data[,-1], byrow=T, bycol=T, rowfirst=F)\
\
#remove redundant ENV variables ...\
names(env21.data)\
env21.cor <- cor(env21.data[,-1])\
print(env21.cor,digits=3)\
source("screen_cor.R")\
screen.cor(env21.data[,-1],threshold=0.8)\
# a possible configuration (other solutions are possible): \
# cull N, Ca, Mg, BS, & Sand as redundant:\
env.cull1 <- env21.data[,-c(11,14,15,19,22)]\
names(env.cull1)\
screen.cor(env.cull1[,-1]) # defaults to a threshold of 0.7\
env.cull2 <- env.cull1[,-c(13,14,17)]\
dim(env.cull2)\
names(env.cull2)\
env13.data <- env.cull2\
rm(env.cull1)\
rm(env.cull2)\
\
#print correlations\
env.cor2m <- cor2m(env13.data[,-1], env13.data[,-1])\
print(as.matrix(env.cor2m), digits=3, na.print="")\
write.csv(env.cor2m,"envcorm.csv")\
#convert the ENV data to z-scores (for later):\
env13.z <- scale(as.matrix(env13.data[,-1]))\
#convert into comparable unit\
# end EDA\
\
#visualize the correlation matrix:\
env13.cor <- cor(env13.data[,-1])\
corrplot(env13.cor, method="color")\
corrplot(env13.cor, method="color", type="upper", order="hclust")\
\
#check normal distribution? PCA is robust to non-normally distributed ones...\
\
#PCA \
#correlation matrix (covariance matrix.) \
env13.pca <- princomp(env13.data[,-1], cor=T, scores=T)\
\
#print eigenvalues\
print(env13.pca)\
#with raw and cumulative variance\
summary(env13.pca) \
# the proportion of variance for PC1 is its SDev^2/13,(its variance/total variance of the dataset(number of variables))\
\
# the actual eigenvalues:\
env13.pca$sdev^2\
#eigenvalue > 1.0 is contributing to data (first 4 PCs)\
\
#print loadings, to see which variables go with each PC (the first 4, 3 significant digits)\
print(round(env13.pca$loadings[,1:4],digits=3))\
#eigenvectors (i.e. regression coefficients)\
\
#biplot\
biplot(env13.pca, xlabs=rep("*",99))\
# the same, but for axes 1 & 3 instead of the default 1 & 2:\
biplot(env13.pca, choices=c(1,3), xlabs=rep("*",99))\
\
# *** the ordination (samples in PC space) ...\
\
#plot samples in PC space\
plot(env13.pca$scores) \
env13.pca.scores <- as.data.frame(env13.pca$scores)\
\
#check SPP correlations with the ENV axes\
print(as.matrix(cor2m(spp17.data[,-1],env13.pca.scores[,1:4])),digits=3,na.print="")\
\
#put the species into PC space (weighted averaging)\
env13.pca.spp11wa <- wascores(env13.pca.scores[,1:4],spp11.rel)\
env13.pca.spp11wa\
write.csv(env13.pca.spp11wa,"spp_pca_WA.csv")\
#plot samples in PCA (factor) space (ignore weird outlier) (1&2)\
plot(env13.pca.scores[,2], env13.pca.scores[,1], xlim=c(-3.5,3.5), ylim=c(-3.5,3.5), xlab="PC 2", ylab="PC1 (Elevation)", pch=19)\
text(env13.pca.spp11wa[,2:1],rownames(env13.pca.spp11wa),cex=0.8,col="blue")\
# correlation vectors (axes in order 2:1):\
pca21.vf <- vf(env13.pca$scores[,c(2,1)],env13.data[,-1],nperm=1000)\
#nperm: random permutation test\
plot.vf(pca21.vf, pval=0.01, ascale=1.0, lwd=2, col="red", length=0.125)\
\
# PC 3:\
plot(env13.pca.scores[,3], env13.pca.scores[,1], xlim=c(-3.5,3.5), ylim=c(-3.5,3.5), xlab="PC 3", ylab="PC1 (Elevation)", pch=19)\
text(env13.pca.spp11wa[,c(3,1)],rownames(env13.pca.spp11wa),cex=0.8,col="blue")\
# correlation vectors (axes in order 3:1):\
pca21.vf <- vf(env13.pca$scores[,c(3,1)],env13.data[,-1],nperm=1000)\
#nperm: random permutation test\
plot.vf(pca21.vf, pval=0.01, ascale=1.0, lwd=2, col="red", length=0.125)\
\
# PC 4:\
plot(env13.pca.scores[,4], env13.pca.scores[,1], xlim=c(-3.5,3.5), ylim=c(-3.5,3.5), xlab="PC 4", ylab="PC1 (Elevation)", pch=19)\
# correlation vectors (axes in order 4:1):\
pca21.vf <- vf(env13.pca$scores[,c(4,1)],env13.data[,-1],nperm=1000)\
#nperm: random permutation test\
plot.vf(pca21.vf, pval=0.01, ascale=1.0, lwd=2, col="red", length=0.125)\
text(env13.pca.spp11wa[,c(4,1)],rownames(env13.pca.spp11wa),cex=0.8,col="blue")\
\
#color bubble plot (1&2)\
plot(env13.pca.scores[,2], env13.pca.scores[,1], xlab="PC2", xlim=c(-3.5,3.5), ylim=c(-3.5,3.5), ylab="PC1 (Elevation)", pch=19, col="gray")\
#highlight a few species, dot size proportional to basal area:\
points(env13.pca.scores[,2:1],pch=19,col="green", cex=3.0*spp11.rel$ABco) # "3" scales, by trial and error\
points(env13.pca.scores[,2:1],pch=19,col="red", cex=3.0*spp11.rel$ABma)\
points(env13.pca.scores[,2:1],pch=19,col="orange", cex=3.0*spp11.rel$PIpo)\
points(env13.pca.scores[,2:1],pch=19,col="blue", cex=3.0*spp11.rel$PIco)\
legend("topleft",horiz=F,legend=c("ABco","ABma","PIpo","PIco"),text.col=c("green","red","orange","blue"),bty="n")\
text(env13.pca.spp11wa[,2:1],rownames(env13.pca.spp11wa),cex=0.8)\
box(lwd=2)\
}
