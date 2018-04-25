# The script is intended to be used from command line with a single parameter
# indicating short country name as shown below
cmd_args <- commandArgs(TRUE)
experiment_name <- cmd_args[1]
models <- cmd_args[2]
league <- cmd_args[3]
year_start <- as.numeric(cmd_args[4]) # "93"
n_cores <- cmd_args[5] # Whatever is available
library('stringi')

# Specified by cmd arguments - here just for testing:
#---------------------------
#league <- "ita"            #|
#year_start <- 1            #|
#experiment_name <- 'test'  #|
#models <- 'estimate_poisson_unconstrained_regularized'
#n_cores <- 1               #|
#---------------------------
# Settings
#season_test <- as.character(as.numeric(season_train) + 101)
season_train <- sprintf("%02d%02d", year_start %% 100, (year_start + 1) %% 100)
season_test <- sprintf("%02d%02d", (year_start + 1) %% 100, (year_start + 2) %% 100)
eval_funcs <- c("logloss", "accuracy")
train_frac <- 0.4 # Fraction of games start predictions from
lambda_min <- 0  # Minimal reg. param. on the grid
lambda_max <- 50 # Maximal reg. param. on the grid
# models <- c(
#   #"estimate_olr_model_regularized",
#   #"estimate_olr_model_unconstrained_regularized"
#   #"estimate_poisson_regularized",
#   #"estimate_poisson_unconstrained_regularized"
#   #"estimate_poisson_regularized_prime",
#   "estimate_poisson_unconstrained_regularized_prime"
#   #"estimate_poisson_DC_regularized",
#   #"estimate_poisson_DC_regularized_prime"
# ) 
leagues <- c(
  "Austria",
  "Denmark",
  "Belgium",
  "England",
  "Finland",
  "France", 
  "Germany",
  #"Greece",
  "Italy", 
  "Ireland",
  "Netherlands",
  "Norway",
  "Poland",
  "Portugal", 
  "Russia",
  "Scotland", 
  "Spain",
  "Romania",
  "Sweden",
  "Switzerland",
  "Turkey"
)

available_leagues <- stri_sub(stri_trans_tolower(leagues), 1, 3)
if(!(league %in% available_leagues))
  stop(paste(c("Run the script with a single parameter for short league name from:", available_leagues), collapse=" "))

dir_results <- file.path('results', experiment_name)
if(!dir.exists(dir_results))
  dir.create(dir_results)
log_file <- file.path('results', experiment_name, paste0(league, '_', season_train, '.log'))

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
                                     #param_grids = lapply(get_parameter_grid_for_league(league, models), function(x) {i <- sample(1:nrow(x), size = 1); return(x[i,,drop = F])}),
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

odds_benchmark(get_filename(league, season_test, leagues, available_leagues), train_frac, as.list(eval_funcs), F)

append_to_report <- function(league, season, rows, experiment_name) {
  results_file <- file.path('results', experiment_name, 'results.csv')
  if(!file.exists(results_file))
    results <- read.csv(file.path('results', 'results_template.csv'))
  else
    results <- read.csv(results_file)
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

append_to_report(league, season_test, best_params_test, experiment_name)
sink()
