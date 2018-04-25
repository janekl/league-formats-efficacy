library("doMC")
library("stringi")
source("league_formats.R")
source("extra_functions.R")
source("simulation_functions.R")
source("config.R")
source("rating_systems/rating_systems.R")
source("schedule_functions.R")

# CLI:
simulation_params <- parse_cl_arguments(commandArgs(trailingOnly = TRUE))
simulation_params[["results_save_folder"]] <- results_save_folder
# Hack [error prone]: Change specific_result_folder dump folder based on the drift_option parameter:
#simulation_params[["specific_result_folder"]] <- ifelse(simulation_params[["drift_option"]] == "fixed", "results_104_fixed_xeon01", "results_104_float_xeon01")
#print(simulation_params)
do.call(run_simulations, simulation_params)

# ------------------------------------------------------------------------
# DEBUG (interactive in RStudio)
# n_cores <- 1
# n_sim <- 10^2 # How to choose this number?
# simulation_params <- list("n" = 16, "model" = "poisson_correlated",
#                          "n_sim" = n_sim, "shape" = 2,
#                          "sigma" = 0.5, "n_cores" = 1, "log2file" = F,
#                          "results_save_folder" = results_save_folder)
# simulation_params[["results_save_folder"]] <- results_save_folder
# simulation_params[["specific_result_folder"]] <- specific_result_folder
# 
# do.call(run_simulations, simulation_params)
