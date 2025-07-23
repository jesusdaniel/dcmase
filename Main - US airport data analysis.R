#######################################
# Main data analysis  results from 
# "Joint Spectral Clustering in
# Multilayer Degree Corrected Blockmodles"
#######################################

# Load functions
source("R/dcmase.R")
source("R/comdet-dcmase.R")
source("R/comdetmethods.R")
source("R/SpectralMethods.R")
source("R/make_dcsbm_plots.R")
source("R/getElbows.R")
library(Matrix)
library(mclust)

load("Data/US_airport_data.RData")

# Symmetric networks (total number of flights between locations)
Adj_list <- lapply(Adj_list, function(A) ((A+t(A))))

# Fit a 4-community multilayer DCSBM
embedding = dcmase(Adj_list = Adj_list)
dim(embedding) # The method suggests a 4-dimensional embedding

comdcmase.res <- comdet_dcmase(Adj_list, K = 4, clustering_method = "kmeans")

# Plot map with airports colored by community
plot_usmap(comdcmase.res$community_memberships)

# Plot time series of degree corrections
plot_degree_corrections(comdcmase.res)

# Plot degree corrections of only one community
plot_degree_corrections_one_community(comdcmase.res, 1, select_vertices = c("LGA", "EWR", "DFW", "ORD", "ATL"))

# Plot block connectivity matrices
plot_pairs_Bmats(comdcmase.res)


#######################################
# Supplement: Compare agreement of communities month vs overall
#######################################
## Communities month by month
comdcmase.res_bymonth = list()
comdcmase.res_bymonth  <- lapply(1:69, function(i) comdet_dcmase(Adj_list[i], K = 4, clustering_method = "kmeans"))
cluster_diff_month = matrix(0, 69, 69)

# compare month to overall
cluster_diff_month_ov = rep(NA, 69)
for(i in 1:69) {
    cluster_diff_month_ov[i] = mclust::classError(comdcmase.res_bymonth[[i]]$community_memberships,
                                                  comdcmase.res$community_memberships)$errorRate
}


plot(cluster_diff_month_ov, xlab = "Month", ylab = "Proportion of nodes", )
abline(v = 50, col = "red")
