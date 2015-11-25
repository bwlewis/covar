require(irlba)
require(parallel)
#' Compute the thresholded correlations between columns of a matrix.
#'
#' Increase p to cut down the total number of candidate pairs evaluated,
#' at the expense of costlier truncated SVDs.
#'
#' @param A an m by n real-valued dense or sparse matrix
#' @param t a threshold value for correlation
#' @param p projected subspace dimension
#'
#' @return A three-column matrix. The first two columns contain
#' indices of vectors meeting the correlation threshold \code{t},
#' and the third column contains the corresponding correlation value.
tcor = function(A, t=0.99, p=5)
{
  mu = colMeans(A)
  s  = sqrt(apply(A,2,crossprod) - nrow(A)*mu^2) # column norms of centered matrix
  L  = irlba(A, p, center=mu, scale=s)
  P  = order(L$v[,1])  # order the entries of v1 (the permutation in the paper)
  limit = sqrt(2*(1-t))/L$d[1]

  v = L$v[P,1]
# linear time longest run search (A. Poliakov):
  lower = 1
  ell = 1
  for(upper in 2:length(v))
  {
    if(v[upper] - v[lower] <= limit)
    {
      ell = max(ell, upper - lower + 1)
    } else
    {
      while(lower < upper && v[upper] - v[lower] > limit) lower = lower + 1
    }
  }


# This is the big union in step 4 of algorithm 2.1, combined with step 6 to
# convert back to original indices, and step 7 to evaluate the candiadtes.
# Each step from 1 to \ell is independent of the others; the steps can run
# in parallel:
  indices = mclapply(1:ell, function(i)
  {
    d = diff(L$v[P,1:p,drop=FALSE],lag=i)^2 %*% L$d[1:p]^2
    j = which(d <= 2*(1-t))  # These ordered indices meet the threshold
    # return original un-permuted column indices that meet threshold (step 7)
    if(length(j) == 0) return()
    v = vapply(j, function(k) cor(A[, P[k]], A[, P[k+i]]), 1)
    h = v >= t
    j = j[h]
    v = v[h]
    cbind(i=P[j], j=P[j+i], cor=v)
  })
  Reduce(rbind,indices)
}
