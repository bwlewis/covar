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
x24 = tcor(A,t,p=24)

d = function(i) diag(s$v[,1:i] %*% diag(s$d[1:i]^2) %*% t(s$v[,1:i]))



p = par(mfrow=c(2,2))
plot(x2$dnum,ylab="",xlab="index difference",main=paste("2-d projection, num. candidates =",nrow(x2$candidates)))
plot(x4$dnum,ylab="",xlab="index difference",main=paste("4-d projection, num. candidates =",nrow(x4$candidates)))
plot(x8$dnum,ylab="",xlab="index difference",main=paste("8-d projection, num. candidates =",nrow(x8$candidates)))
plot(x24$dnum,ylab="",xlab="index difference",main=paste("24-d projection, num. candidates =",nrow(x24$candidates)))
par(p)


# Shows how conservative the sqrt(2*(1-t))/d[1] bound is
#A = matrix(rnorm(100*500),nrow=100)
#C = cor(A)
#C[lower.tri(C,diag=TRUE)] = 0  # cut out symmetric part and ones on diagonal
#t = 0.80*max(C)
#i = which(C>=t,arr.ind=TRUE)
#A = scale(A,center=TRUE,scale=sqrt(apply(A,2,crossprod)))
#s = svd(A)
#lim = sqrt(2*(1-t))/s$d[1]
#x = abs(s$v[i[,1],1] - s$v[i[,2],1])
#plot(x, ylim=c(0,lim*1.2))
#abline(h=lim,col=2)
