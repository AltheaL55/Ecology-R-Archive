# NMS ordinations in R ...
# this focuses on demos with ECODIST, plus some VEGAN bits.

# workflow for an NMS ordination analysis ...
# (1) do data transformations/relativizations as appropriate
# (2) create distance matrix
# (3) do step-down run to choose number of axes K
# (4) redo ordinations at dimension K, and choose the best one
# (5) rotate the solution for ease of interpretation
# (6) compute R2 for ordination axes
# (7) overlay species scores by weighted averaging
# (8) overlay environmental factors as correlations and/or vectors


# Alternatives ...
# You can do the ordination in ECODIST or VEGAN,
# but we're only using VEGAN helper functions, for now ... 
library(vegan)
# Build a distance matrix in VEGAN:
# spp11.bcd<-vegdist(spp11.rel, method="bray")
# is exactly the same as the ECODIST call:
library(ecodist)
spp11.bcd<-distance(as.matrix(spp11.rel), method="bray")

# does the distance measure saturate on this gradient?
length(spp11.bcd) # how big is the distance matrix? (the distance matrix is stored as a vector)
length(spp11.bcd[spp11.bcd==1]) # how many saturated?
length(spp11.bcd[spp11.bcd==1])/length(spp11.bcd) # proportion?
# Extended distances ...
# (this is in VEGAN and only in VEGAN):
# for every sample pair with D=1, replace D with stepping-stone distance:
spp11.xbcd<-stepacross(spp11.bcd) # defaults to shortest paths

# a step-down procedure in ECODIST ...
# Generate 60 ordinations, stepdown, from 1 to 6 dimensions;
# this does min:max=1:6 dimensions, 10 iterations each:
spp11.nms.step <- nmds(spp11.xbcd, nits=10, mindim=1, maxdim=6)
attributes(spp11.nms.step)

# Pull out the stress values and R2 values for each
# dimension (# of axes), and plot them
# (new function plot.nmds in ecodist does all this):
plot.nmds(spp11.nms.step)
# this is a scree plot of stress and a plot of R^2 vs # axes K.
# choose K based on either/both.
# (export plot for your Supporting Materials for explaining why choose N(e.g.2) axis)
dev.off()

# Do a final configuration ...
# does 20 reps (iterations), for 2 dimensions only.
spp11.nmds <- nmds(spp11.xbcd,mindim=2,maxdim=2,nits=20)
# choose the best of these 20 replicate ordinations:
(s.min <- which.min(spp11.nmds$stress)) # returns might be different for you
spp11.nmds$stress[s.min] # returns the stress for the best one;
# note this is scaled on [0,1]:
spp11.nms <- nmds.min(spp11.nmds,dims=2) # grabs the best of 20 reps
#STRESS?

# Rotate the NMS ordination ...
# this forces axis 1 to have the most variance, etc.
# (make sure you understand why this works):
nms.pca<-princomp(spp11.nms)
print(nms.pca)
summary(nms.pca)
spp11.nms<-nms.pca$scores
colnames(spp11.nms) <- c("NMS1", "NMS2")
#PCA scores does not matter but needed in species composition

# Plot the Shepard diagram ...
# compute distances apart for sample pairs in ordination space:
nms2.xod<-dist(spp11.nms)
plot(nms2.xod,spp11.xbcd,pch="*",xlab="Ordination Distance", ylab="Extended B-C Distance")
abline(0,1,col="red") # put in the 1:1 line (intercept=0, slope=1)
title("Shepard Diagram")
box(lwd=2)
# in the figure, xBC distances > 1.0 have been extended (stepped across).
#(if not extended, all ecological distance would end up with 1)
# the correlation between ordination distance and ecological distance ...
# a Mantel correlation is a correlation between distance matrices:
ecodist::mantel(nms2.xod~spp11.xbcd,nperm=10000,nboot=0)

# Get R2 for NMS ordinations ...
# You have to do this by differencing for NMS, cuz
# the axes are computed simultaneously (not as sequential residuals):
# Get OD for axes ...
nms.od1 <- dist(spp11.nms[,1])
nms.od2 <- dist(spp11.nms[,1:2])
# axis 1 is OK as is:
r1<-cor(spp11.xbcd,nms.od1)
r2.1<-r1^2; r2.1
# axis 2 is 2-D minus 1-D solution:
r2<-cor(spp11.xbcd,nms.od2)
r2.2<-r2^2; r2.2-r2.1; r2.2
#r2 individually and accumulative with r1
# capture these numbers for your report: % variance on each axis,
# and cumulative for both axes. 


# Post-processing and interpretation ...

# Correlate ordination axes with ENV:
cor2m(env13.data[,-1],spp11.nms)
# Flip axis 1, if you need to for cosmetic reasons
# (the correlation of NMS1 with elevation should be positive;
# use the next command ONLY IF this correlation is negative):
spp11.nms[,1] <- spp11.nms[,1]* -1.0
# save the table:
spp11.nms.cor2m<-as.matrix(cor2m(env13.data[,-1],spp11.nms))
print(spp11.nms.cor2m,digits=3,na.print="")
write.csv(spp11.nms.cor2m,"spp11_nms_cor2m.csv")

# back to VEGAN for wgt'd avg scores for the species
# (note, VEGAN functions are still attached):
# library(vegan)
# use the raw data to compute the SPP scores 
spp11.nms.wa <- wascores(spp11.nms,spp11.data[,-1])
spp11.nms.wa

# plot axes 1 & 2, with NMS1 being the vertical axis (= elevation):
plot(spp11.nms[,2:1],pch=19, xlab="NMS 2",ylab="NMS 1")
#points(spp11.nms.wa[,2:1],pch="+",lwd=2, cex=1.25,col="blue")
text(spp11.nms.wa[,2:1],rownames(spp11.nms.wa),cex=0.8, col="blue")

# Correlation vectors with ENV:
spp11.nms.vf<-vf(spp11.nms[,2:1],env13.data[,-1],nperm=1000)
plot.vf(spp11.nms.vf, pval=0.01, col="red", lwd=2, length=0.067)
# (length specifies the size of the arrowheads, in inches, and it annoys R)
box(lwd=2)
# (export this graph for your report)

# if you wanna go crazy, a color bubble plot of species abundances:
# re-issue the PLOT command above, then:
points(spp11.nms[,2:1],pch=19,col="green", cex=2.0*spp11.rel$ABco) # "2" scales, by trial and error
points(spp11.nms[,2:1],pch=19,col="red", cex=2.0*spp11.rel$ABma)
points(spp11.nms[,2:1],pch=19,col="orange", cex=2.0*spp11.rel$PIpo)
points(spp11.nms[,2:1],pch=19,col="blue", cex=2.0*spp11.rel$PIco)
points(spp11.nms[,2:1],pch=19,col="cyan", cex=2.0*spp11.rel$PImo)
points(spp11.nms[,2:1],pch=19,col="brown", cex=2.0*spp11.rel$PIje)
points(spp11.nms[,2:1],pch=19,col="dark green", cex=2.0*spp11.rel$PIla)
points(spp11.nms[,2:1],pch=19,col="gray", cex=2.0*spp11.rel$CAde)
points(spp11.nms[,2:1],pch=19,col="sea green2", cex=2.0*spp11.rel$SEgi)

# edit the legend names and colors if you use different species:
legend.txt <- c("ABco","ABma","PIpo","PIco","PImo","PIje","PIla","CAde","SEgi")
legend.col <- c("green","red","orange","blue","cyan","brown","dark green","gray","sea green2")
legend("topright",horiz=F,legend=legend.txt,text.col=legend.col,bty="n")
box(lwd=2)
# export this as a JPG or EMF to use it in a report,
#the location of sample sites in the ordination space could help find control-treatment group
#also could change composition of the sample and see where the dot move in the graph

# explore this elevation gradient ...
# plot species along axis 1:
summary(spp11.nms) # get min and max for axis 1
# see arg xlim and how it effectively removes all the 0's in plotting:
plot(spp11.nms[,1],spp11.rel$ABco,xlim=c(-1.0,1.5),xlab="NMS 1",ylim=c(0.1,1.2),ylab="Abundance",pch=19,col="green")
points(spp11.nms[,1],spp11.rel$ABma,pch=19, col="red")
points(spp11.nms[,1],spp11.rel$PIpo,pch=19, col="orange")
points(spp11.nms[,1],spp11.rel$PIco,pch=19, col="blue")
points(spp11.nms[,1],spp11.rel$PImo,pch=19, col="cyan")
points(spp11.nms[,1],spp11.rel$PIje,pch=19, col="brown")
points(spp11.nms[,1],spp11.rel$PIla,pch=19, col="dark green")
points(spp11.nms[,1],spp11.rel$CAde,pch=19, col="gray")
points(spp11.nms[,1],spp11.rel$SEgi,pch=19, col="sea green2")

# if you use the same species and colors, you can reuse the legend from above
legend("top",horiz=T,legend=legend.txt,text.col=legend.col,bty="n", cex=0.9)
box(lwd=2)
# (export and save, if you like)

# export NMS scores for use elsewhere ...
#write.csv(spp11.nms, "spp11_nms.csv")

# end of NMS


# back-track to previous ordinations and check Shepard diagrams ...
#used to be used for evaluating the restoration progress
# library(ecodist)
# PCA ...
# ordination distance:
env13.pc2.od <- dist(env13.pca.scores[,1:2]) # just on 2 PCs
# environmental dissimilarities as Mahalanobis distances:
env13.md <- distance(env13.data[,-1],"mahal")
# Shepard diagram:  ordination distance vs environmental distance:
plot(env13.pc2.od,env13.md, pch="*",xlab="Distance in PC space",ylab="Environmental Distance (MD)")
# the correlation:
ecodist::mantel(env13.pc2.od~env13.md,nperm=10000,nboot=0)
env13.pc4.od <- dist(env13.pca.scores[,1:4]) # on 4 PCs
plot(env13.pc4.od,env13.md, pch="*",xlab="Distance in PC space",ylab="Environmental Distance (MD)")
ecodist::mantel(env13.pc4.od~env13.md,nperm=10000,nboot=0)

# and to illustrate what we've learned about PCA and species data ...
spp11.pca <- princomp(spp11.rel,cor=T,scores=T)
plot(spp11.pca)
summary(spp11.pca)
spp11.pc2.od <- dist(spp11.pca$scores[,1:2])
plot(spp11.pc2.od,spp11.bcd,pch="*",xlab="Distance on 2 PCs",ylab="Compositional Distance (BC)")
mantel(spp11.pc2.od~spp11.bcd,nperm=10000,nboot=0)
# try extended BC distances ...
plot(spp11.pc2.od,spp11.xbcd,pch="*",xlab="Distance on 2 PCs",ylab="Compositional Distance (BC)")
mantel(spp11.pc2.od~spp11.xbcd,nperm=10000,nboot=0)

# end of all this


