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
sum(C > t)
dim(x)
