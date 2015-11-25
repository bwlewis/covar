# EXAMPLE USE of tcor
source("tcor.R")

library(biclust)  # for the EisenYeast and BicatYeast data
#data(BicatYeast)  # in gene by sample orientation (a tiny example for testing)
#A = t(BicatYeast) # to correlate genes to genes
data(EisenYeast)  # slightly larger example
A = t(EisenYeast) # to correlate genes to genes

t = 0.95

t1 = proc.time()
C = cor(A)                     # actual correlation matrix for reference
C[lower.tri(C,diag=TRUE)] = 0  # cut out symmetric part and ones on diagonal
C = C[C>t]
print(proc.time() - t1)

t1 = proc.time()
# Compute candidate pairs with the new algorithm
x = tcor(A, t=t, p=10)
print(proc.time() - t1)

# Compare
print(sum(C > t))
print(nrow(x$cor))
