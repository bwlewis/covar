# EXAMPLE USE of tcor
source("tcor.R")


library(biclust)  # for the EisenYeast and BicatYeast data
data(BicatYeast)  # in gene by sample orientation (a tiny example for testing)
A = t(BicatYeast) # to correlate genes to genes

C = cor(A)                     # actual correlation matrix for reference
C[lower.tri(C,diag=TRUE)] = 0  # cut out symmetric part and ones on diagonal

# Compute candidate pairs with the new algorithm
t = 0.95
x = tcor(A, t=t, p=10)

# Compare
# First note, there are 115 pairs with better than 0.95 correlation:
sum(C>t)

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
ans = apply(x,1,paste,collapse=",")

# Now the actual list of correlated pairs as strings (already ordered):
actual = apply(which(C>t, arr.ind=TRUE),1,paste,collapse=",")

# Now check that all the actual correlated pairs are accounted for in the
# candidates:
length(intersect(actual, ans$candidates))
# It's 115, so all the actual correlated pairs were found (plus some extra)!
