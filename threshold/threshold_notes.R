library(s4vd)    # for the lung data
library(biclust) # for the EisenYeast and BicatYeast data
library(zoo)
#data(lung200)
#A = t(lung200)
#data(EisenYeast)
#A = t(EisenYeast)
data(BicatYeast)
A = t(BicatYeast)

C = cor(A) # for reference

# Explicitly center and scale the data matrix
A = sweep(A,2,apply(A,2,mean))
A = sweep(A,2,sqrt(apply(A,2,crossprod)),FUN=`/`)
s = svd(A)

t = 0.95 # correlation threshold

# permutation ordering of the first basis vector
i = order(s$v[,1])

# Naive and inefficient way to find the longest run of 1-d values within
# the correlation tolerance (IMPROVE ?).
width = 1000
limit = sqrt(2*(1-t))/s$d[1]   # 1-d tolerance
group_sizes = rollapply(s$v[i,1],width=width, FUN=function(x) sum(x - x[1] < limit), align="left", partial=TRUE)
max_group_size = max(group_sizes)
# if max_group_size=width, then we mis-estimated. Increase width and try again until
# max_group_size < width to be sure that we've found the maximum.

cat("Maximum group size",max_group_size,"\n")


# This plot shows each vector projected down to a point and ordered in
# increasing order...
plot(s$v[i,1], main="1-d ordered", ylab="v")
# ...and the interval beyond which we can say vectors corresponding to those
# points won't be correlated above the threshold value, here shown at zero.
abline(h=c(0,sqrt(2*(1-t))/s$d[1]),col=2)




# This plot shows which sequentially-adjacent points (lag=1) fall above/below
# the threshold and are therefore correlated or not. Each point in this plot
# represents the difference between a *pair* of ordered points (a pair of
# projected vectors).
k = 10    # subspace dimension
d = diff(s$v[i,1:k,drop=FALSE],lag=1)^2 %*% s$d[1:k]^2
ylim = c(0,max(max(d), 1.01*2*(1-t)))
plot(d,ylim=ylim)
abline(h=2*(1 - t),col=2,lwd=1)
j = which(d < 2*(1-t))
points(j,d[j],pch=19,col="#ff000044")
print(length(j))
