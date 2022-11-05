# Load functions
source("R/dcmase.R")
source("R/comdet-dcmase.R")
source("R/SpectralMethods.R")
source("R/make_dcsbm_plots.R")
library(Matrix)


load("Data/US_airport_data.RData")

# Symmetrice networks (total number of flights between locations)
Adj_list <- lapply(Adj_list, function(A) ((A+t(A))))


# Check that matrix is symmetric
A = (Adj_list[[1]])
isSymmetric(A)
plot(svd(A)$d, main = "Singular values of first graph")

# Fit a 4-community multilayer DCSBM
comdcmase.res <- comdet_dcmase(Adj_list, K = 4, clustering_method = "kmeans", scaled = T)

# Plot map with airports colored by community
plot_usmap(comdcmase.res$community_memberships)

# Plot time series of degree corrections
plot_degree_corrections(comdcmase.res)

# Plot degree corrections of only one community
plot_degree_corrections_one_community(comdcmase.res, 1, select_vertices = c("LGA", "EWR", "DFW", "ORD", "ATL"))

# Plot block connectivity matrices
plot_pairs_Bmats(comdcmase.res)




