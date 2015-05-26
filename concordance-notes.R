# Set up an example matrix x with sampled iid rows and known covariance structure
v = qr.Q(qr(matrix(rnorm(100),10)))
d = exp(-(1:10))
sigma = v %*% diag(d) %*% t(v)
library("mvtnorm")
x = rmvnorm(n=500,mean=1:10,sigma=sigma)

# Generate some row samples
set.seed(1)
i = lapply(rep(seq(from=2,to=500,by=1),1), function(j) list(j=j,idx=sample(500,j)))

# I think that 'concordance' is trying to quantify this
# kind of relationship:
nrm = "2"
ncx = norm(cov(x),nrm)
a = t(sapply(i, function(j) c(j$j,norm(cov(x[j$idx,]),nrm))))
#plot(a[,1], abs(a[,2] - ncx)/ncx, main="rel difference")

# or maybe this ratio relationship:
plot(a[,1], abs(a[,2]/norm(cov(x),nrm)), main="ratio")


# A somewhat crude bound for this is:
s = svd(x)
alpha = t(sapply(i, function(j) c(j$j,(499/(j$j-1))*norm(crossprod(s$u[j$idx,]),nrm))))

lines(alpha[,1],alpha[,2], col=2,lwd=2)

# Note that althought ratio looks like 1/x, it's not. You can see this by plotting:
# plot(alpha[,1], 1/alpha[,2])
# That funky shape is from the || U_k^T U_k|| term (see the LaTeX notes).
