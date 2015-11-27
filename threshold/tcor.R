require(irlba)
require(doParallel)

#' Compute the thresholded correlations between columns of a matrix.
#'
#' Increase p to cut down the total number of candidate pairs evaluated,
#' at the expense of costlier truncated SVDs.
#'
#' @param A an m by n real-valued dense or sparse matrix
#' @param t a threshold value for correlation
#' @param p projected subspace dimension
#' @param cores number of CPU cores to use in explicitly parallel code sections
#'  (note that available BLAS libraries may implicitly use multiple threads)
#' @param additional arguments passed to \code{\link{irlba}}
#'
#' @return A list with three elements:
#' \enumerate{
#'   \item \code{indices} A three-column matrix such that the  first two columns contain
#'         indices of vectors meeting the correlation threshold \code{t},
#'         and the third column contains the corresponding correlation value.
#'   \item \code{longest_run} The largest number of successive entries in the
#'     ordered first singular vector within a projected distance defined by the
#'     correlation threshold.
#'   \item \code{n} The total number of _possible_ vectors that meet
#'     the correlation threshold identified by the algorithm.
#'   \item \code{total_time} Total run time.
#' }
tcor = function(A, t=0.99, p=10, cores=detectCores(), ...)
{
  if(ncol(A) < p) p = max(1, floor(ncol(A) / 2 - 1))
  t0 = proc.time()
  mu = colMeans(A)
  s  = sqrt(apply(A, 2, crossprod) - nrow(A) * mu ^ 2) # col norms of centered matrix
  if(any(s < 10*.Machine$double.eps)) stop("the standard deviation is zero for some columns")
  L  = irlba(A, p, center=mu, scale=s, ...)
  t1 = (proc.time() - t0)[[3]]
  P  = order(L$v[, 1])  # order the entries of v1 (the permutation in the paper)
  limit = sqrt(2 * (1 - t)) / L$d[1]
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
# in parallel, for which we use foreach with an R socket cluster.
  combine = function(x, y)
  {
    list(idx=rbind(x$idx, y$idx), n=x$n + y$n)
  }
  cl = makePSOCKcluster(cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))

  indices = foreach(i=1:ell, .combine=combine, .inorder=FALSE) %dopar%
  {
    d = diff(L$v[P, 1:p, drop=FALSE], lag=i) ^ 2 %*% L$d[1:p] ^ 2
    # These ordered indices meet the projected threshold:
    j = which(d <= 2 * (1 - t))
    n = length(j)
    # return original un-permuted column indices that meet true threshold
    # (step 7), including the number of possible candidates for info.:
    if(n == 0)
    {
      ans = vector("list", 2)
      names(ans) = c("idx", "n")
      ans$n = n
      return(ans)
    }
    v = vapply(j, function(k) cor(A[, P[k]], A[, P[k+i]]), 1)
    h = v >= t
    j = j[h]
    v = v[h]
    return(list(idx=cbind(i=P[j], j=P[j + i], cor=v), n=n))
  }
  c(indices, longest_run=ell, irlb_time=t1, total_time=(proc.time() - t0)[[3]])
}
