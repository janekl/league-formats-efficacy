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

# Number of points for win/draw
POINTS_WIN <- 3
POINTS_DRAW <- 1

# Functions to run league simulations
# Standard k-RR_n tournament
one_stage_league <- function(ratings_start, points_start, schedule, rounds, resolve_ties = T, drift = 0.0, update_last = F) {
  n <- length(points_start)
  # Array to store match and league results - see description in  function compute_league_table2():
  results <- array(0, dim = c(n, n, 2)) #matrix(NA, nrow = n, ncol = n)
  matches_at_home <- numeric(n)
  league_points <- points_start
  ranks_tied <- numeric(n - 1)
  teams <- 1:length(ratings_start)
  ratings_rounds <- matrix(NA, nrow = length(rounds), ncol = n)
  ratings <- ratings_start
  #------------------
  # Compute results
  for(r in 1:length(rounds)) {
    for(m in rounds[[r]]) { # For every game in round
      # Sample result of a match
      team_i <- schedule[m, "Team1"]
      team_j <- schedule[m, "Team2"]
      # Store the result and add points for the teams
      updated_results <- sample_result(ratings, team_i, team_j, results, league_points)
      results <- updated_results[["results"]]
      league_points <- updated_results[["league_points"]]
      # For checking balance in Scotland
      matches_at_home[team_i] <- matches_at_home[team_i] + 1
    }
    # Save and update ratings
    ratings_rounds[r,] <- get_ratings(ratings)
    if (r < length(rounds) || update_last) { # Do not change after the last round
      ratings <- update_ratings(ratings, drift)  
      #drift <- exp(rnorm(1, log(drift), tau)) # Stochastic variance model as in Glickman "Dynamic paired comparison models with stochastic variances"
    }
  }
  # Resolve ties
  # Currently done on head-to-head result basis. 
  # For Possion family we have information on goals scored that could be used.
  if(resolve_ties) {
    ties <- ties_solver(league_points, n, results)
    league_points <- ties[["league_points"]]
    ranks_tied <- ranks_tied + ties[["ranks_tied"]]
  }
  # Home, Draw, Away counts - The representation of results array doesn't allow for this: 
  # Information who plays at home is lost - in fact no longer useful.
  #count_1x2 <- as.table(results)
  return(list("ratings_start" = ratings_start,
              "league_points" = league_points, 
              "ratings_final" = ratings,
              "ratings_rounds" = ratings_rounds,
              "results" = results,
              "ranks_tied" = ranks_tied,
              "matches_at_home" = matches_at_home))
}

# Two stage formats (for example: Ekstraklasa since 2013/14)
two_stage_league <- function(ratings_start, points_start, schedule_1, rounds_1, schedule_2, rounds_2, determine_schedule = F, resolve_ties = T, divide_points = T, drift = 0.0) {
  n <- length(points_start)
  #stopifnot(n %% 2 == 0)
  league_points <- points_start #rep(0, n)
  ratings <- ratings_start
  ratings_rounds <- matrix(NA, nrow = length(rounds_1) + length(rounds_2), ncol = n)
  # First stage
  stage_first <- one_stage_league(ratings, points_start, schedule_1, rounds_1, resolve_ties = T, drift = drift, update_last = T)
  points_all <- stage_first[["league_points"]]
  ratings <- stage_first[["ratings_final"]]
  ratings_rounds[1:length(rounds_1),] <- stage_first[["ratings_rounds"]]
  team_order <- order(points_all, decreasing = T)
  points_all <- floor(points_all) # After ordering we ignore the fraction added to points to resolve ties in the first round
  if(divide_points)
    points_all <- round(points_all / 2)
  n2 <- n / 2
  upper_teams <- team_order[1:n2]
  lower_teams <- team_order[(n2 + 1):n]
  points_champ <- points_all[upper_teams]
  ratings_champ <- get_model_subset(ratings, upper_teams)
  drift_champ <- drift[upper_teams]
  points_releg <- points_all[lower_teams]
  ratings_releg <- get_model_subset(ratings, lower_teams)
  drift_releg <- drift[lower_teams]
  # Final round
  if(determine_schedule) {
    # This is temporal solution: for Scotland, schedule for final round 
    # is determined (in particular, variables 'schedule_2' and 'rounds_2')
    # only later. But the number of rounds is known in advance: this is 
    # a hack, rounds_2 is overwritten later anyway.
    second_stage <- determine_schedule_scotland(schedule_1, upper_teams, lower_teams)
    schedule_2_champ <- second_stage[["schedule_2_champ"]]
    schedule_2_releg <- second_stage[["schedule_2_releg"]]
    rounds_2 <- second_stage[["rounds_2"]]
    stage_champ <- one_stage_league(ratings_champ, points_champ, schedule_2_champ, rounds_2, resolve_ties = F, drift = drift_champ, update_last = F)
    stage_releg <- one_stage_league(ratings_releg, points_releg, schedule_2_releg, rounds_2, resolve_ties = F, drift = drift_releg, update_last = F)
  } else {
    stage_champ <- one_stage_league(ratings_champ, points_champ, schedule_2, rounds_2, resolve_ties = F, drift = drift_champ, update_last = F)
    stage_releg <- one_stage_league(ratings_releg, points_releg, schedule_2, rounds_2, resolve_ties = F, drift = drift_releg, update_last = F)
  }
  # Points gained 
  # NOTE: max points to gain in second stage in case of 2RR tournament is 3 * 2 * (n / 2 - 1) =  3 * (n - 2).
  # This is added + 1 for upper teams' tally. It assures that both halves are not mixed.
  points_all[upper_teams] <- stage_champ[["league_points"]] + 3 * (n - 2) + 1 
  points_all[lower_teams] <- stage_releg[["league_points"]]
  ratings <- replace_ratings(ratings, upper_teams, stage_champ[["ratings_final"]])
  ratings <- replace_ratings(ratings, lower_teams, stage_releg[["ratings_final"]])
  ratings_rounds[length(rounds_1) + 1:length(rounds_2), upper_teams] <- stage_champ[["ratings_rounds"]] # rounds_2 start from 34 for Scotland!!!!
  ratings_rounds[length(rounds_1) + 1:length(rounds_2), lower_teams] <- stage_releg[["ratings_rounds"]]
  # Merge results from both rounds
  results <- stage_first[["results"]]
  results[upper_teams, upper_teams, ] <- results[upper_teams, upper_teams, ] + stage_champ[["results"]]
  results[lower_teams, lower_teams, ] <- results[lower_teams, lower_teams, ] + stage_releg[["results"]]
  matches_at_home <- stage_first[["matches_at_home"]]
  matches_at_home[upper_teams] <- matches_at_home[upper_teams] + stage_champ[["matches_at_home"]]
  matches_at_home[lower_teams] <- matches_at_home[lower_teams] + stage_releg[["matches_at_home"]]
  # Ties resolution
  if(resolve_ties) {
    # Ties are resolved based on head-to-head results in both stages
    ties <- ties_solver(points_all, n, results)
    points_all <- ties[["league_points"]]
    ranks_tied <- ties[["ranks_tied"]]
  }
  # Home, Draw, Away counts - not used anymore.
  #count_1x2 <- add_tables(stage_first[["count_1x2"]], stage_champ[["count_1x2"]], stage_releg[["count_1x2"]])
  return(list("ratings_start" = ratings_start,
              "league_points_first" = stage_first[["league_points"]],
              "league_points" = points_all, 
              "ratings_final" = ratings,
              "results_stage01" = stage_first[["results"]],
              "results_stage02" = list("results_champ" = stage_champ[["results"]], "results_releg" = stage_releg[["results"]]),
              "results" = results,
              "ratings_rounds" = ratings_rounds,
              "upper_teams" = upper_teams,
              "lower_teams" = lower_teams,
              "ranks_tied" = ranks_tied,
              "matches_at_home" = matches_at_home))
}

compute_league_table2 <- function(results) {
  ###################################
  # INPUT:
  # results - results of matches given in array of dim n_teams x n_teams x 2. Entry results[t1, t2, 1] stores
  #           the number of matches won by t1 over t2 and results[t1, t2, 2] stores the number of draws.
  #           The table should be initialized to 0.
  # OUTPUT:
  # league_points - a numeric vector
  ###################################
  league_points <- sapply(1:nrow(results), function(i){
    # Wins + draws at entries results[t1, t2, 2] and results[t2, t1, 2]
    return(3 * sum(results[i,,1], na.rm = TRUE) + sum(results[i,,2], na.rm = TRUE) + sum(results[,i,2], na.rm = TRUE))
  })
  return(league_points)
}

flip_within_group <- function(schedule, team_mapping, n_games_home) {
  flip <- logical(nrow(schedule))
  for(k in 1:nrow(schedule)) {
    i_org <- team_mapping[schedule$Team1[k]]
    j_org <- team_mapping[schedule$Team2[k]]
    flip[k] <- n_games_home[i_org, j_org] == 2
  }
  if(any(flip))
    schedule[flip, c("Team1", "Team2")] <- schedule[flip, c("Team2", "Team1")]
  return(schedule)
}

determine_schedule_scotland <- function(schedule_1, upper_teams, lower_teams) {
  # Assumption: The schedule is balanced so that 4th game is played in such a way
  # that teams i and j play equal number of games home and away (2). This is different 
  # somewhat than the official rules.
  stopifnot(length(upper_teams) == length(lower_teams))
  n2 <- length(upper_teams)
  n <- 2 * n2
  # Determine number of games played at home so far
  n_games_home <- matrix(0, nrow = n, ncol = n)
  for(k in 1:nrow(schedule_1)) {
    i <- schedule_1$Team1[k]
    j <- schedule_1$Team2[k]
    n_games_home[i, j] <- n_games_home[i, j] + 1
  }
  schedule_final <- prepare_schedule_1RRn(n2)
  schedule_2_champ <- flip_within_group(schedule_final[["schedule"]], upper_teams, n_games_home)
  schedule_2_releg <- flip_within_group(schedule_final[["schedule"]], lower_teams, n_games_home)
  return(list("schedule_2_champ" = schedule_2_champ,
              "schedule_2_releg" = schedule_2_releg,
              "rounds_2" = schedule_final[["rounds"]]))
}

ties_solver <- function(league_points, n, results, eps = 0.05) {
  #################### 
  # INPUT:
  # league_points
  # results - as specified in compute_league_table2() function.
  # OUTPUT 
  # league_points - updated points with added eps points to break ties. Note: eps should be small enough
  #                 so that teams do not leapfrog in the table. For n <= 20 teams eps of 0.05 is safe since there will
  #                 be maximum of 20 teams tied so at most 19 * 0.5 = 0.95 changes in points.
  # ranks_tied    - vector with updated counts with rank tied (i-th entry stores counts for a series of i + 1 tied teams).
  # NOTE:
  # Ties are resolved on points head-to-head basis, goal difference is ignored (if not resolved on points, then
  # random guess is applied). It is implemented by adding a fraction (eps) of point for teams ranked higher. This might
  # be optimized depending on later use. It is also worth reporting the number of tied teams in each simulation.
  #################### 
  ranks_tied <- numeric(n - 1)
  if(anyDuplicated(league_points)) {
    team_order <- order(league_points)
    league_points_sorted <- league_points[team_order]
    i <- 1
    j <- i
    while(j < n) {
      while(j < n & league_points_sorted[i] == league_points_sorted[j + 1]) {
        j <- j + 1
      }
      if(i != j) {
        ranks_tied[j-i] <- ranks_tied[j-i] + 1
        tied_teams <- team_order[i:j]
        head2head <- compute_league_table2(results[tied_teams, tied_teams, ])
        # We use 'rank' function since it can resolve teams with equal number of points randomy (as opposed to 'order')
        head2head <- order(rank(head2head, ties.method = "random"), decreasing = TRUE)
        # Add epsilon of points to distinguish between teams
        league_points[tied_teams[head2head]] <- league_points[tied_teams[head2head]] + eps * ((length(tied_teams) : 1) - 1)
      }
      i <- j + 1
      j <- i
    }
  }
  return(list("league_points" = league_points, "ranks_tied" = ranks_tied))
}

sample_result <- function(rating_model, team_i, team_j, results, league_points) {
  # Match result
  game_result <- sample(c(1, 0, -1), size = 1, prob = hda_probs(rating_model, team_i, team_j))
  # Store the result and add points for the teams
  if(game_result == 0) {
    results[team_i, team_j, 2] <- results[team_i, team_j, 2] + POINTS_DRAW
    league_points[c(team_i, team_j)] <- league_points[c(team_i, team_j)] + POINTS_DRAW
  } else {
    winner <- c(team_i, team_j)[game_result]
    loser <- c(team_i, team_j)[-game_result]
    results[winner, loser, 1] <- results[winner, loser, 1] + 1
    league_points[winner] <- league_points[winner] + POINTS_WIN
  }
  return(list("results" = results, "league_points" = league_points))
}
