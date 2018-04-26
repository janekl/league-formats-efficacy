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

# Grid search for model parameters

library("stringi")

# -----------------------------------------------------------------------------
# Parameters - specify settings in this code block
league <- "sco"
year_start <- 15 # That is 2015
experiment_name <- "grid_search"
n_cores <- 3
eval_funcs <- c("logloss", "accuracy")
train_frac <- 0.4 # Fraction of games to start predictions from
lambda_min <- 0   # Minimal reg. param. on the grid
lambda_max <- 50  # Maximal reg. param. on the grid
models <- c(
  "estimate_olr_regularized",
  "estimate_poisson_regularized",
  "estimate_poisson_correlated_regularized"
)
# -----------------------------------------------------------------------------

season_train <- sprintf("%02d%02d", year_start %% 100, (year_start + 1) %% 100)
season_test <- sprintf("%02d%02d", (year_start + 1) %% 100, (year_start + 2) %% 100)

leagues <- c(
  "England",
  "France", 
  "Germany",
  "Italy", 
  "Poland",
  "Scotland", 
  "Spain"
)

available_leagues <- stri_sub(stri_trans_tolower(leagues), 1, 3)
if(!(league %in% available_leagues))
  stop(paste(c("Run the script with a single parameter for short league name from:", available_leagues), collapse = " "))

dir_results <- file.path("results", experiment_name)
if(!dir.exists(dir_results))
  dir.create(dir_results, recursive = TRUE)
log_file <- file.path("results", experiment_name, paste0(league, '_', season_train, ".log"))
results_file <- file.path("results", experiment_name, paste0(league, '_', season_train, ".csv"))

source("rating_systems.R")
source("prediction_functions.R")

get_filename <- function(league, season, leagues, available_leagues) {
  paste0(leagues[available_leagues == league], season, ".csv")
}

drop_columns <- function(df_list, columns) {
  library("dplyr")
  lapply(df_list, function(df) {
    select(df, -one_of(columns))
  })
}

sink(log_file)
# Evaluation for a given league
best_params_train <- full_evaluation(models = models,
                                     # TEST:
                                     #param_grids = lapply(get_parameter_grid_for_league(league, models), function(x) {i <- sample(1:nrow(x), size = 3); return(x[i,,drop = F])}),
                                     param_grids = get_parameter_grid_for_league(league, models, lambda_min = lambda_min, lambda_max = lambda_max),
                                     games_csv = get_filename(league, season_train, leagues, available_leagues),
                                     train_frac = train_frac,
                                     eval_funcs = eval_funcs,
                                     n_cores = n_cores)

# Test set performance
best_params_test <- full_evaluation(models = models,
                                    param_grids = drop_columns(best_params_train, eval_funcs),
                                    games_csv = get_filename(league, season_test, leagues, available_leagues),
                                    train_frac = train_frac,
                                    eval_funcs = eval_funcs)

# Betting odds performance
#odds_benchmark(get_filename(league, season_test, leagues, available_leagues), train_frac, as.list(eval_funcs), F)

save_report <- function(league, season, rows, results_file) {
  results <- data.frame(
    league = character(), 
    season = character(), 
    model = character(), 
    lambda = numeric(), 
    corr = numeric(), 
    logloss = numeric(), 
    accuracy = numeric()
  )
  for(model in names(rows)) {
    x <- rows[[model]]
    x["model"] <- model
    x["league"] <- league
    x["season"] <- season_test
    for(col in colnames(results))
      if(!(col %in% names(x)))
        x[col] <- NA
    results <- rbind(results, x[colnames(results)])
  }
  write.csv(results, results_file, row.names = F)
}

save_report(league, season_test, best_params_test, results_file)
sink()
