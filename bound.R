# Set up an example matrix x with sampled iid rows and known covariance structure
v = qr.Q(qr(matrix(rnorm(100),10)))
d = exp(-(1:10))
sigma = v %*% diag(d) %*% t(v)
library("mvtnorm")
x = rmvnorm(n=500,mean=1:10,sigma=sigma)

# Generate some row samples
set.seed(1)
i = lapply(rep(seq(from=2,to=500,by=1),1), function(j) list(j=j,idx=sample(500,j)))

# Plot the ratio of row sample covariance norms to the full matrix covariance norm
ncx = norm(cov(x),"2")
a = t(sapply(i, function(j) c(j$j,norm(cov(x[j$idx,]),"2"))))
plot(a[,1], a[,2]/norm(cov(x),"2"), main="Covariance norm ratio", xlab="sample size", ylab="ratio")

# Our upper bound
s = svd(x)
alpha = t(sapply(i, function(j) c(j$j,(499/(j$j-1))*norm(crossprod(s$u[j$idx,]),"2"))))
lines(alpha[,1],alpha[,2], col=2,lwd=2)
