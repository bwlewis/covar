# EXPERIMENTAL PROTOTYPE FUNCTION
# For an example use of this, see the file example.R.
#
# As of 20-July, 2015 this needs the development version of irlba available
# from GitHub: devtools::install_github("bwlewis/IRL")
require(irlba)

#' Compute the thresholded correlations between columns of a matrix.
#'
#' Increase p to cut down the total number of candate pairs.
#'
#' @param A an m by n real-valued dense or sparse matrix
#' @param t a threshold value for correlation
#' @param p projected subspace dimension
#'
#' @return
#' A list with entries:
#' candidates: A two-column matrix. Each row of the matrix are the
#'             indices of candidate vectors that meet the correlation threshold.
#' dnum: a vector of length equal to the longest run of adjacent points
#'       in v_1 meeting the threshold that contains the number of candidate
#'       entries making the threshold cut at each index difference.
tcor = function(A, t=0.99, p=5, reduce=1)
{
  mu = colMeans(A)
  s  = sqrt(apply(A,2,crossprod) - nrow(A)*mu^2) # column norms of centered matrix
  L  = irlba(A, p, center=mu, scale=s)
  P  = order(L$v[,1])  # order the entries of v1 (the permutation in the paper)
  limit = sqrt(2*(1-t))/L$d[1]
# XXX EXPERIMENTAL reduce limit by fraction of variance in 1-d
#if(reduce) limit = limit*L$d[1]^2/ncol(A) 
limit = limit*reduce
# sum of squared singular values of centered/scaled A = ncol(A)

# XXX this is stupid, find a better way ('l' or LaTeX '\ell' in the paper):
  v = L$v[P,1]
  ell = max(vapply(1:length(v), function(i) {x=v[-(1:i)];sum(x-x[1] <= limit)}, 1))

cat("longest run of adjacent points within threshold l =",ell,"\n")

  dnum = rep(0,ell)
# this is the big union in step 4 of algorithm 2.1, combined with step 6 to
# convert back to original indices:
  candidates = lapply(1:ell, function(i)
  {
    d = diff(L$v[P,1:p,drop=FALSE],lag=i)^2 %*% L$d[1:p]^2
    j = which(d <= 2*(1-t))  # These ordered indices meet the threshold
    dnum[i] <<- length(j)    # Record the number of them
    # return original un-permuted column indices
    cbind(P[j], P[j+i])
  })
# Return candidates and dnum
  list(candidates=Reduce(rbind,candidates), dnum=dnum)
}
