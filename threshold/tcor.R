# EXPERIMENTAL PROTOTYPE FUNCTION
#
# As of 20-July, 2015 this needs the development version of irlba available
# from GitHub: devtools::install_github("bwlewis/IRL")
require(irlba)
require(zoo)   # for rollapply

#' Compute the thresholded correlations between columns of a matrix.
#'
#' Increase p to cut down the total number of candate pairs.
#'
#' @param A an m by n real-valued dense or sparse matrix
#' @param t a threshold value for correlation
#' @param p projected subspace dimension
#'
#' @return
#' A two-column matrix. Each row of the matrix are the
#' indices of candidate vectors that meet the correlation threshold.
tcor = function(A, t=0.99, p=5)
{
  mu = colMeans(A)
  s  = sqrt(apply(A,2,crossprod))
  L  = irlba(A, p, center=mu, scale=s)
  P  = order(L$v[,1])  # order the entries of v1 (the permutation in the paper)
  limit = sqrt(2*(1-t))/L$d[1]
# XXX this is stupid, find a better way ('l' in the paper):
  max_group_size  = max(rollapply(L$v[P,1], width=nrow(A), FUN=function(x) sum(x-x[1] <= limit), align="left", partial=TRUE))
cat("longest run of adjacent points within threshold l =",max_group_size,"\n")

# this is the big union in step 4 of algorithm 2.1, combined with step 6 to
# convert back to original indices:
  candidates = lapply(1:max_group_size, function(i)
  {
    d = diff(L$v[P,1:p,drop=FALSE],lag=i)^2 %*% L$d[1:p]^2
    j = which(d <= 2*(1-t))  # These ordered indices meet the threshold
    # return original un-permuted column indices
    cbind(P[j], P[j+i])
  })
  Reduce(rbind,candidates) # Return a two-column array of index pairs
}

# EXAMPLE USE


library(biclust)  # for the EisenYeast and BicatYeast data
data(BicatYeast)  # in gene by sample orientation (a tiny example for testing)
A = t(BicatYeast) # to correlate genes to genes

C = cor(A)                     # actual correlation matrix for reference
C[lower.tri(C,diag=TRUE)] = 0  # cut out symmetric part and ones on diagonal

# Compute candidate pairs with the new algorithm
x = tcor(A, t=0.95, p=10)

# Compare
# First note, there are 115 pairs with better than 0.95 correlation:
sum(C>0.95)

# and our algorithm with p=10 returns a list of 431 possible candidates:
dim(x)

# We need to verify that x at least contains the 115 good pairs! Here is
# one cheesy way to tell. First note that the entries of x are in somewhat
# arbitrary order defined by the first singular vector. For example, the
# pair (15,3) might be listed but not the equivalent redundant pair (3,15).
# We account for that and build the following list of candidate
# pairs in order so that the first coordinate is less than the second:
i = which(x[,2] < x[,1])
x[i,] = x[i,2:1]
# Express the coordinates as strings
candidates = apply(x,1,paste,collapse=",")

# Now the actual list of correlated pairs as strings (already ordered):
actual = apply(which(C>0.95, arr.ind=TRUE),1,paste,collapse=",")

# Now check that all the actual correlated pairs are accounted for in the
# candidates:
length(intersect(actual, candidates))
# It's 115, so all the actual correlated pairs were found (plus some extra)!

# Note that even many of the extras tend to be highly correlated (just not
# meeting the threshold)!
summary(C[x])
