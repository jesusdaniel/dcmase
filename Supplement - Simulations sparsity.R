#######################################
# Supplementary simulation results from "Joint Spectral Clustering in
# Multilayer Degree Corrected Blockmodles"
#######################################

# These simulations compare the performance
# of different multilayer clustering methods
# in terms of edge density

source("Experiments/run_all_methods.R")
source("Experiments/simulations sparsity.R")
library(Matrix)
library(igraph)
library(mclust)



num_replications <- 100
ave_degs <- seq(2, 24, 2)

# Same B same theta
parameters_list <- ave_degs
param_iter = ave_degs/150
results_simulation1 <- iterate_parameters(sim_setting = simulation1s, parameters_list, param_iter, num_replications)
make_ggplot_res_error(results_simulation1,  "Edge density", ylim = c(0,1), xbreaks = c(0.01, 0.04, 0.07, 0.1, 0.13, 0.16))

# Different B same theta
parameters_list <- ave_degs
param_iter = ave_degs/150
results_simulation2 <- iterate_parameters(sim_setting = simulation2s, parameters_list, param_iter, num_replications)
make_ggplot_res_error(results_simulation2,  "Edge density", ylim = c(0,1), xbreaks = c(0.01, 0.03, 0.05, 0.07, 0.09, 0.11))

# Diff B diff theta
parameters_list <- ave_degs
param_iter = ave_degs/150
results_simulation3 <- iterate_parameters(simulation3s, parameters_list, param_iter, num_replications)
make_ggplot_res_error(results_simulation3,  "Edge density", ylim = c(0,1), xbreaks = c(0.01, 0.03, 0.05, 0.07, 0.09, 0.11))

# Same B different theta
parameters_list <- ave_degs
param_iter = ave_degs/150
results_simulation4 <- iterate_parameters(simulation4s, parameters_list, param_iter, num_replications)
make_ggplot_res_error(results_simulation4,  "Edge density", ylim = c(0,1), xbreaks = c(0.01, 0.03, 0.05, 0.07, 0.09, 0.11))



#############################################

# Alternating degrees
param_iter = ave_degs/150
parameters_list <- lapply(ave_degs, 
                          function(x) c(x, 150, 3))
results_simulation6 <- iterate_parameters(simulation6s, parameters_list, param_iter, num_replications)
make_ggplot_res_error(results_simulation6,  "Edge density", ylim = c(0,1), xbreaks = c(0.01, 0.03, 0.05, 0.07, 0.09, 0.11))

results_simulation7 <- iterate_parameters(simulation7s, parameters_list, param_iter, num_replications)
make_ggplot_res_error(results_simulation7,  "Edge density", ylim = c(0,1),xbreaks = c(0.01, 0.03, 0.05, 0.07, 0.09, 0.11))




#######################################
results_simulation1$scenarioB <- "Same B"
results_simulation2$scenarioB <- "Different B"
results_simulation3$scenarioB <- "Different B"
results_simulation4$scenarioB <- "Same B"
results_simulation6$scenarioB <- "Same B"
results_simulation7$scenarioB <- "Different B"

results_simulation1$scenarioT <- "Same \u0398"
results_simulation2$scenarioT <- "Same \u0398"
results_simulation3$scenarioT <- "Different \u0398"
results_simulation4$scenarioT <- "Different \u0398"
results_simulation6$scenarioT <- "Alternating \u0398"
results_simulation7$scenarioT <- "Alternating \u0398"

different_scenarios <- rbind(results_simulation1, results_simulation2, results_simulation4, results_simulation3, 
                             results_simulation6, results_simulation7)

different_scenarios$scenarioB <- factor(different_scenarios$scenarioB,
                                        levels = c("Same B", "Different B"))
different_scenarios$scenarioT <- factor(different_scenarios$scenarioT,
                                        levels = c("Same \u0398", "Different \u0398", "Alternating \u0398"))
#save(different_scenarios, file = "Results-allmethods-rep100-sparsity.RData")
# load("Experiments/Results-allmethods-rep100-sparsity.RData")

#png("Figures/Simulation-rep100-6scenarios-sparsity.png", width = 1200, height = 750, res = 150)
make_ggplot_multipleBT2(different_scenarios, "Edge density", xbreaks = c(0.01, 0.04, 0.07, 0.1, 0.13, 0.16),ylim = c(0, 0.65),
                        methodnames = c("DC-MASE", "Sum of adj. matrices", "Bias-adjusted SoS",  "MASE", "OLMF", "graph-tool")) 
#dev.off()


