# Example parallel computation of correlation matrix entries
#
# The test here does not do anything with the results. See the note
# below within the foreach loop.
#
library('foreach')
library('doSNOW')
nodelist = c("n1","n2","n3","n4")

m = 200     # Number of rows
n = 100000  # Number of columns

A = matrix(0,nrow=m, ncol=n)
t1 = proc.time()
set.seed(1)
A[,] = runif(m*n)
dt = (proc.time() - t1)[[3]]
cat("Time to create random entries",dt,"\n")

t1 = proc.time()
x = apply(A,2,mean)
s = apply(A,2,sd)
for(j in 1:length(x)) A[,j] = A[,j] - x[j]

# In practice, we'll want a copy of the transpose and the original matrix for
# speed. A 1M x 1K matrix consumes only about 8GB memory.  So, with this
# extra copy we need at least about 16GB for that sized problem.
TA = t(A)
# We also scale the columns of A (but not TA) in advance for efficiency's sake
A = t(t(A)/s)
dt = (proc.time() - t1)[[3]]
cat("Time to center matrix and compute sd",dt,"\n")

# Save A and distribute to other nodes
t1 = proc.time()
save(A, file="/dev/shm/A.mat")
save(TA, file="/dev/shm/TA.mat")
for(j in nodelist[-1]) {
  cmd = paste("scp -r /dev/shm/A.mat ",j,":/dev/shm/",sep="")
  system(cmd)
  cmd = paste("scp -r /dev/shm/TA.mat ",j,":/dev/shm/",sep="")
  system(cmd)
}
dt = (proc.time() - t1)[[3]]
cat("Time to distribute data",dt,"\n")


matinit = function()
{
  if(!exists("A",envir=globalenv())) {
    load("/dev/shm/A.mat", envir=globalenv())
    load("/dev/shm/TA.mat", envir=globalenv())
  }
}

# The subcor function computes a subset of the correlation matrix of A, where
# A is a matrix and assumes columns are centered and scaled by sd,
# TA is the transpose of unscaled A (either matrix or big.matrix),
# s is a vector of columnwise standard deviations of A,
# rows is the desired correlation matrix row subset,
# cols is the desired correlation matrix column subset.
# A matrix of size length(rows) x length(cols) is returned.
subcor = function(A, TA, s, rows, cols)
{
  n = nrow(A) - 1
  X = TA[rows,] %*% A[,cols]
  X = diag(1/(n*s[rows])) %*% X
  dimnames(X)[[1]] = rows
  dimnames(X)[[2]] = cols
  X
}
# The subcor function may be used in at least two parallel code sections:
# 1. Compute blocks of rows of the big correlation matrix for subsequent
#    thesholding and evaluation.
# 2. Computation of many mini-adjacency matrices corresponding to a given
#    row of the thresholded large correlation matrix.

# Let's just consider #1 for now, computed in blocks of rows:
blk = 500
# First, in parallel w/all nodes
cl = makeCluster(nodelist)
registerDoSNOW(cl)
t1 = proc.time()
w=foreach(j=seq(1,n,by=blk),.noexport=c("A","TA")) %dopar%
{
t1=proc.time()
  matinit()
  cls = 1:n
  k = min(n,j+blk-1)
  S = subcor(A, TA, s, j:k, cls)
# Here is where we threshold the subset of the correlation matrix and
# then form adjacency matrices, etc., which is another task easily
# computed in parallel. However, SNOW is not the best framework for
# nested parallelism.
  NULL
}
dt = (proc.time() - t1)[[3]]
cat("Wall time, 4 nodes:",dt,"\n")
stopCluster(cl)

# Now with 2 nodes:
cl = makeCluster(nodelist[c(3,4)])
registerDoSNOW(cl)
t1 = proc.time()
w=foreach(j=seq(1,n,by=blk),.noexport=c("A","TA")) %dopar%
{
  matinit()
  cls = 1:n
  k = min(n,j+blk-1)
  S = subcor(A, TA, s, j:k, cls)
  NULL
}
dt = (proc.time() - t1)[[3]]
cat("Wall time, 2 nodes:",dt,"\n")
stopCluster(cl)

# And finally, sequentially:
registerDoSEQ()
t1 = proc.time()
w=foreach(j=seq(1,n,by=blk),.noexport=c("A","TA")) %dopar%
{
  matinit()
  cls = 1:n
  k = min(n,j+blk-1)
  S = subcor(A, TA, s, j:k, cls)
  NULL
}
dt = (proc.time() - t1)[[3]]
cat("Wall time, 1 node :",dt,"\n")
