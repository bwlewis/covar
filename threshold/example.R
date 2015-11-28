source("tcor.R")

library(biclust)  # for the EisenYeast and BicatYeast data
data(BicatYeast)  # in gene by sample orientation (a tiny example for testing)
A = t(BicatYeast) # to correlate genes to genes
threshold = 0.95


# Compute using tcor():
m1  = sum(gc()[,2])  # Memory in megabytes currently in use
t1  = proc.time()
tx  = tcor(A, t=threshold, p=10)
t1  = (proc.time() - t1)[3] # Time in seconds
m1  = sum(gc()[,6]) - m1    # Peak excess memory use during the test

cat("tcor time (s)                      ", t1, "\n")
cat("tcor peak excess memory (MB)       ", m1, "\n")

# Compute using R's native cor() function:
m2  = sum(gc()[,2])
t2  = proc.time()
cx  = cor(A)                # full correlation matrix
diag(cx) = -Inf
pairs = which(cx >= threshold, arr.ind=TRUE)
t2  = (proc.time() - t2)[3]
m2  = sum(gc()[,6]) - m2

cat("brute force time (s)               ", t2, "\n")
cat("brute force peak excess memory (MB)", m2, "\n")

# Compare results!
cat("Number of vector pairs meeting threshold found by brute force", nrow(pairs)/2, "\n")
cat("Number of vector pairs meeting threshold found by tcor       ", nrow(tx[[1]]), "\n")
