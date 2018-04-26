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

# Results evaluation
library("stringi")
library("dplyr")
library("xtable")
source("config.R")
source("evaluation_functions.R")

evaluation <- function(result_folder, n_teams, model, team_strength_fun, infinity_number, verbose = T) {
  # Compute csv with results
  results <- list()
  league_formats <- list.dirs(file.path(results_save_folder, result_folder, n_teams, model), recursive = F, full.names = F)
  result_files_pattern <- "0\\.\\d+_-?\\d+(\\.\\d)?\\.rds"
  for(league_format in league_formats) {
    if(verbose) cat("League:", league_format, "\n")
    almost_full_path <- file.path(results_save_folder, result_folder, n_teams, model, league_format)
    for(results_file in dir(file.path(almost_full_path, "points"), pattern = result_files_pattern)) {
      if(verbose) cat("Processing:", results_file, "\n")
      ratings <- readRDS(file.path(almost_full_path, "ratings", results_file))
      # Here comes the definition of the strongest team. If it was already computed when saving results 
      # in write_results() function during simulations, this function should be set to identity
      team_strength <- lapply(ratings, team_strength_fun)
      points_team <- readRDS(file.path(almost_full_path, "points", results_file))
      params <- suppressWarnings(as.numeric(stri_split_fixed(stri_replace_all_regex(results_file, "\\.rds$", ""), "_")[[1]]))
      results[[length(results) + 1]] <- list("system" = league_format, 
                                             "sigma" = params[1], 
                                             "shape" = params[2], 
                                             "the_best_win" = the_best_win(points_team, team_strength),
                                             "kendall_tau" = kendall_tau(points_team, team_strength),
                                             "spearman_footrule" = spearman_footrule(points_team, team_strength))
    }
  }
  results <- data.frame(do.call(rbind, results))
  # Something strange - values in columns are... lists: converting to vectors manually
  for(col in colnames(results)) results[[col]] <- as.vector(results[[col]], mode = class(results[[col]][[1]]))
  results$shape[results$shape == -1] <- infinity_number
  return(results)
}

evaluate_all <- function(specific_result_folders, models, teams_number) {
  for(model in models) {
    for(teams in teams_number) {
      for(results in specific_result_folders) {
        filename_out <- paste0(results, "_", model, "_", teams, ".csv")
        cat(filename_out, "...\n")
        results <- evaluation(results, teams, model, team_strength_fun = function(x) x, infinity_number = 1000, verbose = F)
        write.csv(results, file.path(results_save_folder, filename_out), row.names = F)
      }
    }
  }
}

evaluate_all("results_all", "poisson_correlated", "12")
