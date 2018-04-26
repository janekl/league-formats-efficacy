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

# Various routines for model optimization

library("dplyr")
library("doMC")

read_data <- function(dataset, directory = file.path('..', '..', 'data')) {
  #cat(filepath, "\n")
  games <- read.csv(file.path(directory, dataset), header = TRUE, stringsAsFactors = FALSE)
  games <- rename(games, home_team = HomeTeam, away_team = AwayTeam)
  games[["result"]] <- ifelse(games$FTHG > games$FTAG, 1, ifelse(games$FTHG < games$FTAG, -1, 0))
  if("Info" %in% colnames(games))
    games <- games[games$Info != "wo",]
  teams <- sort(unique(c(games$home_team, games$away_team)))
  games$home_team <- as.numeric(factor(games$home_team, levels = teams))
  games$away_team <- as.numeric(factor(games$away_team, levels = teams))
  return(games)
}

read_data2 <- function(dataset, train_frac = 0.4, directory = file.path('..', '..', 'data')) {
  # This function extends `read_data()` for the number of teams and round start for predictions
  games <- read_data(dataset, directory = directory)
  n_teams <- length(unique(games$home_team))
  round_start <- 2
  while(mean(games$Round < round_start) < train_frac) {
    round_start <- round_start + 1
  }
  list(games = games, n_teams = n_teams, round_start = round_start)
}

get_parameter_grid_for_league <- function(league, models, lambda_min = 0, lambda_max = 50) {
  param_grids <- list(
    "estimate_olr_regularized" = data.frame("lambda" = seq(lambda_min, lambda_max, by = 0.5)),
    "estimate_poisson_regularized" = data.frame("lambda" = seq(lambda_min, lambda_max, by = 0.5)),
    "estimate_poisson_correlated_regularized" = expand.grid("lambda" = seq(lambda_min, lambda_max, by = 0.5), "corr" = 0.45) # "corr" = c(0.0, 0.01, seq(0.05, 0.95, by = 0.05), 0.99)),
  )
  for(model in names(param_grids)) {
    if(!(model %in% models))
      param_grids[[model]] <- NULL # Removes from list
  }
  stopifnot(all(models %in% names(param_grids)))
  return(param_grids)
}

accuracy <- function(predictions, results) {
  return(mean(2 - apply(predictions, 1, which.max) == results))
}

logloss <- function(predictions, results) {
  prob <- predictions[matrix(c(1:nrow(predictions), 2 - results), ncol = 2)]
  return(-mean(log(prob)))
}

evaluation <- function(model, param_grid, games_csv, train_frac, eval_funcs, n_cores, sort_metric = "logloss") {
  cat(model, "\n")
  results_all <- list()
  for(f in games_csv) {
    data_all <- read_data2(f, train_frac = train_frac)
    results_model <- evaluate_model_parallel(model, param_grid, data_all[['games']], data_all[['round_start']], data_all[['n_teams']], eval_funcs, n_cores)
    results_all[[f]] <- results_model
  }
  results_all <- do.call("rbind", results_all)
  results_all <- arrange_(results_all, .dots = sort_metric)
  return(results_all)
}

full_evaluation <- function(models, param_grids, games_csv, train_frac, eval_funcs, n_cores = 1, print_top = 50, return_top = 1) {
  # Function to compte predictions in sliding window manner.
  # INPUT:
  # * return_top - either integer value or string 'all' - how many results to return (sorted according to a given metric)
  print(games_csv)
  best_params <- list()
  for(model in models) {
    results <- evaluation(model = model, param_grid = param_grids[[model]], games_csv = games_csv, train_frac = train_frac, eval_funcs = eval_funcs, n_cores = n_cores)
    print(head(results, print_top))
    if(return_top == 'all') {
      best_params[[model]] <- results
    } else {
      best_params[[model]] <- head(results, return_top)
    }
  }
  return(invisible(best_params))
}

evaluate_model_sequential <- function(rating_model, parameter_grid, games, round_start, n_teams, eval_funcs, n_cores = NULL) {
  results <- matrix(NA, nrow(parameter_grid), length(eval_funcs), dimnames = list(NULL, eval_funcs))
  for(i in 1:nrow(parameter_grid)) {
    # Get parameters
    params <- list("alpha" = 0, "teams" = 1:n_teams)
    for(p in colnames(parameter_grid)) 
      params[p] <- parameter_grid[i, p]
    # Estimate model
    #print(params)
    results[i, ] <- league_prediction(games, rating_model, params, eval_funcs, round_start = round_start) #[[1]]
  }
  results <- cbind(parameter_grid, results)
  return(results)
}

evaluate_model_parallel <- function(rating_model, parameter_grid, games, round_start, n_teams, eval_funcs, n_cores) {
  registerDoMC(n_cores)
  n <- nrow(parameter_grid)
  results <- foreach(i = 1:n) %dopar% {
    params <- list("alpha" = 0, "teams" = 1:n_teams) # Only L2 regularization: alpha == 0
    for(p in colnames(parameter_grid)) params[p] <- parameter_grid[i, p]
    output <- league_prediction(games, rating_model, params, eval_funcs, round_start = round_start)
    for(p in colnames(parameter_grid)) output[p] <- parameter_grid[i, p]
    output
  }
  results <- do.call(rbind, results)
  data.frame(results, row.names = NULL)
}

league_prediction <- function(games, model, params, eval_funcs = c("logloss", "accuracy"), round_start = 17, odds_benchmark = F) {
  # Allocate matrix for predictions
  predictions <- matrix(NA, nrow = sum(games$Round >= round_start), ncol = 3, dimnames = list(rownames(games[games$Round >= round_start,]), c("h", "d", "a")))
  n_rounds <- max(games$Round)
  for(r in round_start:n_rounds) {
    params[["games"]] <- games[games$Round < r, c("home_team", "away_team", "FTHG", "FTAG", "result")]
    test_set <- games[games$Round == r, c("home_team", "away_team")]
    ratings <- do.call(model, params)
    predictions[rownames(test_set),] <- t(apply(test_set, 1, function(x) hda_probs(ratings, x[1], x[2])))
  }
  err_model <- numeric(length(eval_funcs))
  for(i in seq_along(eval_funcs))
    err_model[i] <- do.call(eval_funcs[i], list(predictions, games[games$Round >= round_start, "result"]))
  names(err_model) <- eval_funcs
  # Przeniesc
  #err_distr <- do.call(eval_fun, list(matrix(table(games$result)[3:1]/nrow(games), nrow = sum(games$Round >= round_start), ncol = 3, byrow = TRUE), 
  #                                    games[games$Round >= round_start, "result"]))
  #cat("Model logloss:\n")
  #cat(str(err_model))
  #cat("Distr logloss:", err_distr, "\n")
  return(err_model)
}

odds_benchmark <- function(games_csv, train_frac, eval_funcs, normalize_odds = TRUE) {
  print(games_csv)
  data_all <- read_data2(games_csv, train_frac = train_frac)
  games <- data_all[['games']]
  round_start <- data_all[['round_start']]
  if(normalize_odds) # For Ekstraklasa we have raw odds which should be normalized
    games[, c("H", "D", "A")] <- 1 / games[,c("H", "D", "A")] / apply(1 / games[,c("H", "D", "A")], 1, sum)
  games[["result"]] <- ifelse(games$FTHG > games$FTAG, 1, ifelse(games$FTHG < games$FTAG, -1, 0))
  in_test <- games$Round >= round_start
  frac_test <- mean(in_test)
  total_test <- sum(in_test)
  games <- filter(games, Round >= round_start)
  odds_results <- c()
  for(eval_fun in eval_funcs) {
    x <- do.call(eval_fun, list(games[,c("H", "D", "A")], games[,"result"]))
    odds_results[eval_fun] <- x
    #cat(eval_fun, ":", x, "\n")
  }
  odds_results["frac_test_games"] <- frac_test
  odds_results["total_test_games"] <- total_test
  return(odds_results)
}

estimate_model <- function(games_csv, rating_system, parameters) {
  games <- read_data(games_csv)
  parameters[["games"]] <- games
  parameters[["teams"]] <- 1:max(games$home_team) # A bit hacky
  parameters[["alpha"]] <- 0
  model <- do.call(rating_system, parameters)
  model
}
