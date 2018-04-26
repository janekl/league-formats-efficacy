## MIT License
## 
## Copyright (c) 2018 Jan Lasek
## 
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
##   
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
## 
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
## -----------------------------------------------------------------------------

# Functions to run simulation for different league formats
rinvgamma <- function(n = 1, shape = 1, scale = 1) {
  1/rgamma(n, shape = shape, scale = 1/scale)
}

get_model_setup <- function(model, n, sigma, shape, seed) {
  set.seed(seed)
  # These constants are defined from the analysis of empirical data.
  # In this analysis they were set as averages of these parameters
  # for three leagues - German, Polish, Scottish for 2016/17 season.
  HTA_OLR <- 0.345
  INTERCEPT_OLR <- 0.582
  HTA_POI <- 0.312
  INTERCEPT_POI <- 0.091
  
  if(model == "olr") {
    ratings <- structure(list("ratings" = rnorm(n, 0, sigma), 
                              "intercept" = INTERCEPT_OLR, "hta" = HTA_OLR), 
                         class = "olr_model")
  }
  if(model == "poisson") {
    ratings <- structure(list("ratings" = matrix(rnorm(2 * n, 0, sigma), ncol = 2, dimnames = list(NULL, c("Att", "Def"))), 
                              "intercept" = INTERCEPT_POI, "hta" = HTA_POI), 
                         class = "poisson_model")
  }
  if(model == "poisson_correlated") {
    corr <- 0.45
    Sigma <- matrix(c(sigma**2, corr*sigma**2)[c(1, 2, 2, 1)], nrow = 2, ncol = 2)
    initial <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
    ratings <- structure(list("ratings" = matrix(initial, ncol = 2, dimnames = list(NULL, c("Att", "Def"))), 
                              "intercept" = INTERCEPT_POI, "hta" = HTA_POI, "corr" = corr), 
                         class = "poisson_model_correlated")
  }
  if(shape == -1) {
    drift <- rep(0, n)
  } else {
    drift <- rinvgamma(n, shape, 1)
  }
  return(list("ratings" = ratings, "drift" = drift))
}

run_simulations <- function(n = 16, model = "olr", sigma_priors = 1, shape_params = 5, drift_option = "fixed", team_strength_fun = function(x) apply(x, 2, mean),
                            n_sim = 10, n_cores = 1, log2file = F, results_save_folder = ".", specific_result_folder = "results", verbose = F) {
  ##############
  # Basic function to run simulations
  # Note:
  # * team properties are new for each simulation
  # INPUT:
  # * drift_option - whether total drift in the league should be constant (fixed). If a numeric value is given then
  #                  the drift is divided over this number
  ##############
  if(log2file) {
    sink(file.path(results_save_folder, "logs", paste0("run_", stri_replace_all_regex(as.character(Sys.time()), ":|\\s", "-"), ".txt")))
    on.exit(sink())
  }
  registerDoMC(n_cores)
  total <- prod(sapply(list(sigma_priors, shape_params), length)) # Total number of iterations
  cat("Parameter grid size / n_teams / n_sim / n_cores:", total, "/", n, "/", n_sim, "/", n_cores, "\n")
  cat("* model:       ", model, "\n")
  cat("* sigma priors:", sigma_priors, "\n")
  cat("* shape params:", shape_params, "\n")

  start_time <- Sys.time()
  i <- 0
  cat("START!\n")
  for(sigma in sigma_priors) {
    for(shape in shape_params) {
      # Different leagues
      ## k-round-robin leagues
      for(k in 1:3) {
        if(verbose) print(c("TYPE K, n:", k, n))
        schedules <- prepare_schedule_kRRn(k, n)
        simulation <- foreach(j = 1:n_sim) %dopar% {
          model_params <- get_model_setup(model, n, sigma, shape, j)
          drift_per_round <- model_params[["drift"]]
          if(drift_option == "fixed") {
            drift_per_round <- drift_per_round / sqrt(length(schedules[["rounds"]]))
          } else {
            drift_per_round <- drift_per_round / drift_option
          }
          results <- one_stage_league(model_params[["ratings"]], rep(0, n),
                                      schedules[["schedule"]], schedules[["rounds"]],
                                      drift = drift_per_round)
          results
        }
        write_results(simulation, results_save_folder, specific_result_folder, n, model, paste0(k, "RR"), team_strength_fun, sigma, shape)
        gc()
      }

      # Two stage systems in the new way
      for(schedule_base in c("poland", "kazakhstan", "scotland")) {
        #league_name <- tail(stri_split_fixed(schedule_generator, "_")[[1]], 1)
        schedules <- do.call(paste0("prepare_schedule_", schedule_base), list("n" = n))
        determine_schedule <- schedule_base == "scotland"
        for(division in c(TRUE, FALSE)) {
          if(verbose) print(c("TYPE", schedule_base, "Divisiom:", division))
          simulation <- foreach(j = 1:n_sim) %dopar% {
            model_params <- get_model_setup(model, n, sigma, shape, j)
            drift_per_round <- model_params[["drift"]]
            if(drift_option == "fixed") {
              drift_per_round <- drift_per_round / sqrt(length(schedules[["rounds_1"]]) + length(schedules[["rounds_2"]]))
            } else {
              drift_per_round <- drift_per_round / drift_option
            }
            results <- two_stage_league(model_params[["ratings"]], rep(0, n),
                                        schedules[["schedule_1"]], schedules[["rounds_1"]],
                                        schedules[["schedule_2"]], schedules[["rounds_2"]],
                                        determine_schedule = determine_schedule,
                                        drift = drift_per_round, divide_points = division)
            results
          }
          write_results(simulation, results_save_folder, specific_result_folder, n, model, paste0(schedule_base, "_", ifelse(division, "half", "full")), team_strength_fun, sigma, shape)
          gc()
        }
      }
      cat(sprintf("\r%.f%%", 100 * (i <- i + 1)/total))
    }
  }
  cat("\nDONE.\n")
  print(Sys.time() - start_time)
  return(invisible())
}
