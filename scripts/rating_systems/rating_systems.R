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

# Functions for estimating different rating systems and generating predictions

library("skellam")

logistic <- function(x){
  return (1/(1 + exp(-x)))
}

estimate_olr_regularized <- function(games, teams, lambda, alpha, newcomers = c()) {
  which_home <- which(games[,"result"] == 1)
  which_draw <- which(games[,"result"] == 0)
  which_away <- which(games[,"result"] == -1)
  n_teams <- length(teams)
  # Log-likelihood function  
  loglik <- function(parameters) {
    params_temp <- parameters[1:n_teams] # No hard constraint on sum of ratings!
    zeta1 <- parameters[n_teams + 1]
    zeta2 <- parameters[n_teams + 2]
    diff_rating <- params_temp[games[,"home_team"]] - params_temp[games[,"away_team"]]
    home <- 1/(1 + exp(zeta1 - diff_rating))
    away <- 1 - 1/(1 + exp(zeta2 - diff_rating))
    draw <- pmax(1 - home - away, 0.0001)
    return(-((sum(log(home[which_home])) + sum(log(draw[which_draw])) + sum(log(away[which_away]))) #/ nrow(games)
             - lambda*((1 - alpha)/2 * sum(params_temp^2)))) # + alpha*sum(abs(params_temp))))) - Only L2 is employed here
  }
  ll_max <- optim(par = c(rep(0, n_teams + 1), -1), fn = loglik, control = list(maxit = 100, trace = F), method = "BFGS")
  # Get parameters and scale such that strengths sum up to 0
  ratings <- ll_max$par[1:n_teams]
  ratings <- ratings - mean(ratings) # Scaling parameters so that mean(rating) == 0
  intercept <- 0.5*(ll_max$par[n_teams + 1] - ll_max$par[n_teams + 2])
  hta <- -0.5*(ll_max$par[n_teams + 2] + ll_max$par[n_teams + 1])
  return(structure(list("ratings" = ratings, "intercept" = intercept, "hta" = hta), class = "olr_model"))
}

estimate_poisson_regularized <- function(games, teams, lambda, alpha) {
  n_teams <- length(teams)
  # Log-likeligood function  
  loglik <- function(parameters) {
    params_temp <- parameters[1:(2 * n_teams)] # Defensive skills of last team IS NOT set to 0
    hta <- parameters[length(parameters)]
    diff_home <- hta + params_temp[games[,"home_team"]] - params_temp[n_teams + games[,"away_team"]]
    diff_away <- params_temp[games[,"away_team"]] - params_temp[n_teams + games[,"home_team"]]
    rate1 <- exp(diff_home)
    rate2 <- exp(diff_away)
    ll <- sum(games[,"FTHG"] * diff_home + games[,"FTAG"] * diff_away - rate1 - rate2) - lambda*((1 - alpha)/ 2 * sum(params_temp^2)) # + alpha * sum(abs(params_temp)))
    return(-ll)
  }
  ll_max <- optim(par=rep(0, 2 * n_teams + 1), fn = loglik, control = list(maxit = 100, trace = F), method = "BFGS")
  # Get parameters and scale such that strengths sum up to 0
  att <- ll_max$par[1:n_teams]
  def <- ll_max$par[n_teams + 1:n_teams]
  mean_att <- mean(att)
  mean_def <- mean(def)
  params <- matrix(c(att - mean_att, def - mean_def), ncol = 2, dimnames = list(teams, c("Att", "Def")))
  intercept <- mean_att - mean_def
  hta <- ll_max$par[length(ll_max$par)]
  return(structure(list("ratings" = params, "intercept" = intercept, "hta" = hta), class = "poisson_model"))
}

estimate_poisson_correlated_regularized <- function(games, teams, lambda, alpha, corr) {
  lambda <- lambda / (1 - corr^2)
  n_teams <- length(teams)
  # Log-likeligood function  
  loglik <- function(parameters) {
    params_temp <- parameters[1:(2 * n_teams)]
    hta <- parameters[length(parameters)]
    diff_home <- hta + params_temp[games[,"home_team"]] - params_temp[n_teams + games[,"away_team"]]
    diff_away <- params_temp[games[,"away_team"]] - params_temp[n_teams + games[,"home_team"]]
    rate1 <- exp(diff_home)
    rate2 <- exp(diff_away)
    ll <- (sum(games[,"FTHG"] * diff_home + games[,"FTAG"] * diff_away - rate1 - rate2)
           - lambda*((1 - alpha)/2 * sum(params_temp^2) - corr *  sum(params_temp[1:n_teams] * params_temp[n_teams + 1:n_teams])))
    return(-ll)
  }
  ll_max <- optim(par=rep(0, 2 * n_teams + 1), fn = loglik, control = list(maxit = 100, trace = F), method = "BFGS")
  # Get parameters and scale such that strengths sum up to 0
  att <- ll_max$par[1:n_teams]
  def <- ll_max$par[n_teams + 1:n_teams]
  mean_att <- mean(att)
  mean_def <- mean(def)
  params <- matrix(c(att - mean_att, def - mean_def), ncol = 2, dimnames = list(teams, c("Att", "Def")))
  intercept <- mean_att - mean_def
  hta <- ll_max$par[length(ll_max$par)]
  return(structure(list("ratings" = params, "intercept" = intercept, "hta" = hta, "corr" = corr), class = "poisson_model"))
}

hda_probs <- function(model, ...) {
  UseMethod("hda_probs", model)
}

# Specific prediction functions
hda_probs.random_model <- function(i, j) {
  return(rep(1, 3) / 3)
}

hda_probs.olr_model <- function(rating_system, i, j) {
  rating1 <- rating_system$ratings[i]
  rating2 <- rating_system$ratings[j]
  intercept <- rating_system$intercept # 0.6
  hta <- rating_system$hta # 0.4 # Home team advantage
  rating_diff <- hta + rating1 - rating2
  prob1 <- logistic(-intercept + rating_diff)
  prob3 <- 1 - logistic(intercept + rating_diff)
  return(c("team1" = prob1, "draw" = 1 - prob1 - prob3, "team2" = prob3))
}

hda_probs.poisson_model <- function(rating_system, i, j) {
  rate1 <- exp(rating_system$intercept + rating_system$hta + rating_system$ratings[i, "Att"] - rating_system$ratings[j, "Def"])
  rate2 <- exp(rating_system$intercept + rating_system$ratings[j, "Att"] - rating_system$ratings[i, "Def"])
  # Implementation using 'skellam' library
  win1 <- 1 - pskellam(0, lambda1 = rate1, lambda2 = rate2)
  draw <- dskellam(0, lambda1 = rate1, lambda2 = rate2)
  win2 <- max(c(1 - win1 - draw, 0)) # Use max to avoid rounding errors near 0
  return(c("team1" = win1, "draw" = draw, "team2" = win2))
}

hda_probs.poisson_model_correlated <- hda_probs.poisson_model

get_ratings <- function(rating_model, ...) {
  UseMethod("get_ratings", rating_model)
}

get_ratings.olr_model <- function(rating_model) {
  return(rating_model[["ratings"]])
}

get_ratings.poisson_model <- function(rating_model) {
  # This is sum of (a, d) strengths
  return(apply(rating_model[["ratings"]], 1, sum))
}

# Same as get_ratings.poisson_model:
get_ratings.poisson_model_correlated <- get_ratings.poisson_model

update_ratings <- function(rating_model, ...) {
  UseMethod("update_ratings", rating_model)
}

update_ratings.olr_model <- function(rating_model, drift) {
  rating_model[["ratings"]] <- rating_model[["ratings"]] + rnorm(length(rating_model[["ratings"]]), sd = drift)
  return(rating_model)
}

update_ratings.poisson_model <- function(rating_model, drift) {
  rating_model[["ratings"]] <- rating_model[["ratings"]] + matrix(rnorm(length(rating_model[["ratings"]]), sd = drift), ncol = 2)
  return(rating_model)
}

update_ratings.poisson_model_correlated <- function(rating_model, drift) {
  corr <- rating_model[["corr"]]
  n <- 2 * length(drift)
  # Create covariance matrix for sampling
  Sigma <- diag(rep(drift**2, each = 2))
  skip <- 2*n+2 # Need to fill every 2*n+2 correlation starting from 2 and n+1
  last <- n**2  # Number of entries in the correlation matrix
  corr_times_var <- corr*drift**2
  Sigma[seq(n+1, last, by = skip)] <- corr_times_var
  Sigma[seq(2, last, by = skip)] <- corr_times_var
  # Generate updates
  updates <- mvrnorm(n = 1, mu = rep(0, n), Sigma = Sigma)
  updates <- matrix(updates, ncol = 2, byrow = TRUE)
  # This is done now w/o looping
  #for(i in seq_along(drift)) {
  #  Sigma <- matrix(c(drift[i]**2, corr*drift[i]**2)[c(1, 2, 2, 1)], nrow = 2, ncol = 2)
  #  updates <- mvrnorm(n = 1, mu = c(0, 0), Sigma = Sigma)
  #  ratings_updated[i,] <- ratings_updated[i,] + updates
  #}
  rating_model[["ratings"]] <- rating_model[["ratings"]] + updates
  return(rating_model)
}

get_model_subset <- function(rating_model, ...) {
  UseMethod("get_model_subset", rating_model)
}

get_model_subset.olr_model <- function(rating_model, subset_teams) {
  rating_model[["ratings"]] <- rating_model[["ratings"]][subset_teams]
  return(rating_model)
}

get_model_subset.poisson_model <- function(rating_model, subset_teams) {
  rating_model[["ratings"]] <- rating_model[["ratings"]][subset_teams,]
  return(rating_model)
}

get_model_subset.poisson_model_correlated <- get_model_subset.poisson_model

replace_ratings <- function(rating_model, ...) {
  UseMethod("replace_ratings", rating_model)
}

replace_ratings.olr_model <- function(rating_model, subset_teams, rating_update) {
  stopifnot(length(subset_teams) == length(rating_update$ratings))
  rating_model$ratings[subset_teams] <- rating_update$ratings
  return(rating_model)
}

replace_ratings.poisson_model <- function(rating_model, subset_teams, rating_update) {
  stopifnot(2 * length(subset_teams) == length(rating_update$ratings))
  rating_model$ratings[subset_teams,] <- rating_update$ratings
  return(rating_model)
}

replace_ratings.poisson_model_correlated <- replace_ratings.poisson_model
