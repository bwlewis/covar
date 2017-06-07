# Example efficient Pearson's sample correlation matrix computation
#
# You should have R linked against a high-performance BLAS library for best
# results. See http://illposed.net/r-on-linux.html

#' @param x input m by n data matrix
#' @result n x n symmetric Pearson's sample correlation matrix
#' @notes compare with cor(x, method="pearson")
pearson = function(x)
{
  m = nrow(x)
  z = .colMeans(x, nrow(x), ncol(x))
  w = apply(x, 2, crossprod)
  W = 1 / sqrt(w - m * z * z)
  e = rep(1.0, m)
  r = (2 * z) %*% crossprod(e, x)
  r = m * tcrossprod(z) - r
  ans = crossprod(x) + r
  sweep((W * ans), 2, W, '*')
}
