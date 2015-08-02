source("tcor.R")  # thresholded correlation prototype function

library(biclust)  # for the EisenYeast and BicatYeast example data
# A tiny example dataset
data(BicatYeast)  # in gene by sample orientation (a tiny example for testing)
A = t(BicatYeast) # to correlate genes to genes

t = 0.95

C = cor(A)                     # actual correlation matrix for reference
C[lower.tri(C,diag=TRUE)] = 0  # cut out symmetric part and ones on diagonal

# Compare solutions in several dimensions
x2  = tcor(A,t,p=2)
x4  = tcor(A,t,p=4)
x8  = tcor(A,t,p=8)
x20 = tcor(A,t,p=20)
x24 = tcor(A,t,p=24)
x32 = tcor(A,t,p=32)

d = function(i) diag(s$v[,1:i] %*% diag(s$d[1:i]^2) %*% t(s$v[,1:i]))



par(mfrow=c(2,2))
plot(x2$dnum,ylab="",xlab="index difference",main=paste("2-d projection, num. candidates =",nrow(x2$candidates)))
plot(x4$dnum,ylab="",xlab="index difference",main=paste("4-d projection, num. candidates =",nrow(x4$candidates)))
plot(x8$dnum,ylab="",xlab="index difference",main=paste("8-d projection, num. candidates =",nrow(x8$candidates)))
plot(x24$dnum,ylab="",xlab="index difference",main=paste("24-d projection, num. candidates =",nrow(x24$candidates)))
