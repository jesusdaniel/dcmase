#######################################
# Supplementary analysis  results from "Joint Spectral Clustering in
# Multilayer Degree Corrected Blockmodels"
#######################################

# This code compares the performance of different clustering methods in real data
# using out-of-sample prediction error

measure_out_of_sample_error <- function(Atest, community_memberships) {
  n <- nrow(Atest)
  K <- max(community_memberships)
  Z <- matrix(0, n, K)
  Z[cbind(1:n, community_memberships)] <- 1
  
  # Estimate model parameters given the clusters -------------------------------

  
  # For degree correct. use d_i/sum(d_i per community) * n_k
  n_k <- colSums(Z)
  degrees <- rowSums(Atest)
  normalizations <- crossprod(Z, degrees)/n_k
  degree_corrections <- degrees/(Z %*% normalizations)
  
  # For connectivity, calculate average edges per cell
  nk_nk <- tcrossprod(n_k)
  B <- as.matrix(crossprod(Z, crossprod(Atest, Z)) / nk_nk)
  
  P <- (diag(degree_corrections[,1]) %*% Z) %*% B %*% t(diag(degree_corrections[,1]) %*% Z)
  
  diff <- Atest - P
  return(sum(diff^2)/n)
}

comparisons_out_of_sample <- function(K) {
  print(K)
  m.oos <- length(Adj_list)
  # Methods to run
  methods = c("dcmase", "ave_spherical", "sq-bias-adjusted",
              "mase-spherical")
  
  out_of_sample <- matrix(NA, length(methods), m.oos)
  
  for(i in 1:m.oos) {
    out_of_sample[, i] <- sapply(methods, function(method) {
      comdet1 <- comdetmethods(Adj_list[-i], K, method = method)
      measure_out_of_sample_error(Adj_list[[i]], comdet1)  
    })
    
    print(rowMeans(out_of_sample, na.rm= T))
  }
  return(out_of_sample)
}



#######################################

# Load functions
source("R/dcmase.R")
source("R/comdet-dcmase.R")
source("R/comdetmethods.R")
source("R/SpectralMethods.R")
source("R/make_dcsbm_plots.R")
source("R/getElbows.R")
source("R/dcmase.R")
source("R/make_ggplot.R")
library(ggthemes)
library(Matrix)
library(mclust)
library(ggplot2)

load("Data/US_airport_data.RData")
methods = c("dcmase", "ave_spherical", "sq-bias-adjusted",
            "mase-spherical")
# Symmetric networks (total number of flights between locations)
Adj_list <- lapply(Adj_list, function(A) ((A+t(A))))


#######################################
# Run comparisons
#######################################
results_oos <- lapply(2:10, comparisons_out_of_sample)  

#save(results_oos, file ="oos-results")

#######################################
# Generate figures
# To only make figures, uncomment next line and skip previous chunk
#load("Experiments/oos-results")
table_oos <- sapply(results_oos, function(u) rowMeans(u))   

colnames(table_oos) = 2:10
rownames(table_oos) = methods
library(reshape2)
oos_melted <-melt(table_oos, varnames = c("Method", "K"))
oos_melted$K <- factor(oos_melted$K)
oos_melted$Method <- factor(oos_melted$Method)

methodnames = c("DC-MASE", "Sum of adj. matrices", "Bias-adjusted SoS",  "MASE")

oos_melted_all <- Reduce(rbind, lapply(1:9, function(i) {
  df <- results_oos[[i]]
  colnames(df) = 1:length(Adj_list)
  rownames(df) = 1:4
  df_melted <- melt(df, varnames = c("Method", "Time"))
  df_melted$K <- i+1
  df_melted
}))


oos_melted_all$Method = factor(oos_melted_all$Method)
oos_melted_all$K = factor(oos_melted_all$K)



df2 <- data_summary(data = oos_melted_all, varname="value", 
                    groupnames=c("Method", "K"))



ggplot(df2, aes(x = K, y = value, fill = Method)) + 
  geom_bar(stat = "identity",  position = position_dodge(),  color="black") +
  scale_fill_manual(labels = methodnames, 
                    values = colorblind_pal()(8)[c(7, 2, 4, 6, 3, 8)]) +
  scale_shape_manual(labels = methodnames, 
                     values = c(19, 17,15, 7, 3, 8)) +
  scale_linetype_manual(labels = methodnames, 
                        values = c(1:6)) +
  theme(legend.position="top", legend.text.align = 0) +
  geom_errorbar(aes(ymin=value-2*se, ymax=value+2*se), width=.2,
                position=position_dodge(.9)) 





### comparisons vs dcmase


n <- ncol(Adj_list[[1]])

oos_melted_all_paired <- Reduce(rbind, lapply(1:9, function(i) {
  df <- results_oos[[i]]
  dfnew <- df[2:4,] - rbind(df[1,],df[1,],df[1,])
  colnames(dfnew) = 1:length(Adj_list)
  rownames(dfnew) = 2:4
  df_melted <- melt(dfnew, varnames = c("Method", "Time"))
  df_melted$K <- i+1
  df_melted
}))

oos_melted_all_paired$Method = factor(oos_melted_all_paired$Method)
oos_melted_all_paired$K = factor(oos_melted_all_paired$K)

gg_paired_dif <- ggplot(oos_melted_all_paired, aes(x = K, y = value/n^2, fill = Method)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  geom_boxplot() +
  scale_fill_manual(labels = c("Sum of adj. vs DC-MASE", "Bias-adjust. SoS vs DC-MASE", "MASE vs DC-MASE"), 
                    values = colorblind_pal()(8)[c(2, 4, 6, 3, 8)]) +
  theme_bw() +
  theme(legend.position="top", legend.text.align = 0) +
  xlab("Number of communities") +
  ylab("Paired MSE difference")

#png("../Figures/USairports-paired-diff.png", width = 1200, height = 800, res = 180)
gg_paired_dif
#dev.off()
