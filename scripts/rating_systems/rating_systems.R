library("MASS")
library("skellam")

#-------------------------------------------------
# Auxiliary functions
# This is borrowed (copy+paste) from "FootballLeagueTablePrediction" project
#-------------------------------------------------

logistic <- function(x){
  return (1/(1 + exp(-x)))
}

# Ordinal logistic regression ratings
estimate_olr_model <- function(games, teams, newcomers) {
  ###################################
  # INPUT:
  # games   - data.frame with columns: HomeTeam, AwayTeam, FTHG, FTAG, Result
  # teams   - teams encoding vector
  # home_team, away_team, result - column names in games for data
  # OUTPUT:
  # Object of class 'olr_model'
  ###################################
  X <- as.data.frame(matrix(0, nrow = nrow(games), ncol = length(teams)))
  rating_team <- paste0("rating_", teams)
  newcomer_id <- paste0("rating_", newcomers)
  colnames(X) <- rating_team
  # Construct a design matrix:
  # Home team indicator
  X[matrix(c(1:nrow(games), games[,"home_team"]), ncol = 2)] <- 1
  # Away team indicator
  X[matrix(c(1:nrow(games), games[,"away_team"]), ncol = 2)] <- -1
  keep_teams <- rep(TRUE, ncol(X))
  for(i in newcomers) {
    # If a newcomer is not present in training set then we just drop a column as 
    # it is full of 0. If it is present (that happens when it gets relegated and 
    # comes back another season), we just keep it.
    keep_teams[i] <- any(X[,i] != 0)
  }
  X <- X[,keep_teams]
  # Reference level: the def strengh of the last team
  X <- data.frame('result' = factor(games[,"result"], ordered = TRUE), X[,-ncol(X)])
  model <- polr(result ~ ., data=X)
  ratings <- numeric(length(teams))
  names(ratings) <- rating_team
  ratings[names(model$coef)] <- model$coef 
  #ratings[length(ratings)]# Fill last team rating = 0
  #names(ratings)[length(teams)] <- paste0("rating_", length(teams)) # What for?
  ratings <- ratings - mean(c(model$coef, 0)) # Scaling parameters so that mean(rating) == 0
  # Assign ratings to newcomers
  ratings[newcomer_id] <- min(ratings[names(model$coef)], 0)
  intercept <- 0.5 * (model$zeta[2] - model$zeta[1])
  hta <- -0.5 * (model$zeta[1] + model$zeta[2])
  return(structure(list("ratings" = ratings, "intercept" = intercept, "hta" = hta), class = "olr_model"))
}

estimate_olr_model_regularized <- function(games, teams, newcomers, lambda, alpha) {
  ###################################
  # INPUT:
  # games   - data.frame with columns: HomeTeam, AwayTeam, FTHG, FTAG, Result
  # teams   - teams encoding vector
  # lambda, alpha - regularization parameters for Elastic Net
  # OUTPUT:
  # Object of class 'olr_model'
  ###################################
  # Precompute indicators for particular events
  which_home <- which(games[,"result"] == 1)
  which_draw <- which(games[,"result"] == 0)
  which_away <- which(games[,"result"] == -1)
  # Log-likelihood function  
  n_teams <- length(teams)
  rating_team <- paste0("rating_", teams)
  #newcomer_id <- paste0("rating_", newcomers)
  loglik <- function(parameters) {
    params_temp <- c(parameters[1:(n_teams-1)], 0) # rating of last team set to 0
    zeta1 <- parameters[n_teams]
    zeta2 <- parameters[n_teams + 1]
    diff_rating <- params_temp[games[,"home_team"]] - params_temp[games[,"away_team"]]
    home <- 1/(1 + exp(zeta1 - diff_rating))
    away <- 1 - 1/(1 + exp(zeta2 - diff_rating))
    draw <- pmax(1 - home - away, 0.0001)
    #print(cbind(home, draw, away))
    return(-((sum(log(home[which_home])) + sum(log(draw[which_draw])) + sum(log(away[which_away]))) #/ nrow(games) 
             - lambda*((1 - alpha) / 2 * sum(params_temp^2) + alpha * sum(abs(params_temp)))))
  }
  # Finding parameters by max likelihood: for n - 1 ratings and two intercepts: length == n + 1
  ll_max <- optim(par = c(rep(0, n_teams), -1), fn = loglik, control = list(maxit = 100, trace = F), method = "BFGS")
  # Get parameters and scale such that strengths sum up to 0 (excluding newcomers)
  ratings <- c(ll_max$par[1:(n_teams-1)], 0)
  #if(length(newcomers)) {
  #  ratings[-newcomers] <- ratings[-newcomers] - mean(ratings[-newcomers]) # Scaling parameters so that mean(rating) == 0
  #  # Assigning ratings to newcomers
  #  ratings[newcomers] <- min(ratings[-newcomers])
  #}
  intercept <- 0.5 * (ll_max$par[n_teams] - ll_max$par[n_teams + 1])
  hta <- -0.5 * (ll_max$par[n_teams + 1] + ll_max$par[n_teams])
  return(structure(list("ratings" = ratings, "intercept" = intercept, "hta" = hta), class = "olr_model"))
}

estimate_olr_model_unconstrained_regularized <- function(games, teams, lambda, alpha, newcomers = c()) {
  ###################################
  # INPUT & OUTPUT: see function 'estimate_olr_model_regularized'
  ###################################
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
             - lambda*((1 - alpha)/2 * sum(params_temp^2)))) # + alpha*sum(abs(params_temp))))) - This is not used in this project
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
  # Function to compute Dixon & Coles model with regularization  
  ###################################
  # INPUT & OUTPUT: see function 'estimate_poisson_DC' + regularization params
  ###################################
  n_teams <- length(teams)
  # Log-likeligood function  
  loglik <- function(parameters) {
    params_temp <- c(parameters[1:(2 * n_teams - 1)], 0) # Defensive skills of last team set to 0
    hta <- parameters[length(parameters)]
    diff_home <- hta + params_temp[games[,"home_team"]] - params_temp[n_teams + games[,"away_team"]]
    diff_away <- params_temp[games[,"away_team"]] - params_temp[n_teams + games[,"home_team"]]
    rate1 <- exp(diff_home)
    rate2 <- exp(diff_away)
    ll <- sum(games[,"FTHG"] * diff_home + games[,"FTAG"] * diff_away - rate1 - rate2) - lambda*((1 - alpha)/ 2 * sum(params_temp^2) + alpha * sum(abs(params_temp)))
    return(-ll)
  }
  ll_max <- optim(par=rep(0, 2 * n_teams), fn = loglik, control = list(maxit = 100, trace = F), method = "BFGS")
  # Get parameters and scale such that strengths sum up to 0
  att <- ll_max$par[1:n_teams]
  def <- c(ll_max$par[(1 + n_teams):(2 * n_teams - 1)], 0)
  mean_att <- mean(att)
  mean_def <- mean(def)
  params <- matrix(c(att - mean_att, def - mean_def), ncol = 2, dimnames = list(teams, c("Att", "Def"))) #mozna odjac?
  intercept <- mean_att - mean_def
  hta <- ll_max$par[length(ll_max$par)]
  return(structure(list("ratings" = params, "intercept" = intercept, "hta" = hta), class = "poisson_model"))
}

estimate_poisson_unconstrained_regularized <- function(games, teams, lambda, alpha) {
  # Function to compute Dixon & Coles model with regularization  
  ###################################
  # INPUT & OUTPUT: see function 'estimate_poisson_DC' + regularization params
  ###################################
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
  params <- matrix(c(att - mean_att, def - mean_def), ncol = 2, dimnames = list(teams, c("Att", "Def"))) #mozna odjac?
  intercept <- mean_att - mean_def
  hta <- ll_max$par[length(ll_max$par)]
  return(structure(list("ratings" = params, "intercept" = intercept, "hta" = hta), class = "poisson_model"))
}


estimate_poisson_unconstrained_regularized_prime <- function(games, teams, lambda, alpha, corr) {
  # Function to compute Dixon & Coles model with regularization  
  ###################################
  # INPUT & OUTPUT: see function 'estimate_poisson_DC' + regularization params
  ###################################
  lambda <- lambda / (1 - corr^2)
  n_teams <- length(teams)
  # Log-likeligood function  
  loglik <- function(parameters) {
    params_temp <- parameters[1:(2 * n_teams)] # Defensive skills of last team IS NOT set to 0
    #print(data.frame(att = params_temp[1:n_teams], def = params_temp[n_teams + 1:n_teams]))
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
  params <- matrix(c(att - mean_att, def - mean_def), ncol = 2, dimnames = list(teams, c("Att", "Def"))) #mozna odjac?
  intercept <- mean_att - mean_def
  hta <- ll_max$par[length(ll_max$par)]
  return(structure(list("ratings" = params, "intercept" = intercept, "hta" = hta, "corr" = corr), class = "poisson_model"))
}

estimate_poisson_regularized_prime <- function(games, teams, lambda, alpha, corr) {
  # Function to compute Dixon & Coles model with regularization  
  ###################################
  # INPUT & OUTPUT: see function 'estimate_poisson_DC' + regularization params
  ###################################
  lambda <- lambda / (1 - corr^2)
  n_teams <- length(teams)
  # Log-likeligood function  
  loglik <- function(parameters) {
    params_temp <- c(parameters[1:(2 * n_teams - 1)], 0) # Defensive skills of last team set to 0
    hta <- parameters[length(parameters)]
    diff_home <- hta + params_temp[games[,"home_team"]] - params_temp[n_teams + games[,"away_team"]]
    diff_away <- params_temp[games[,"away_team"]] - params_temp[n_teams + games[,"home_team"]]
    rate1 <- exp(diff_home)
    rate2 <- exp(diff_away)
    ll <- (sum(games[,"FTHG"] * diff_home + games[,"FTAG"] * diff_away - rate1 - rate2)
           - lambda*((1 - alpha)/2 * sum(params_temp^2) - corr *  sum(params_temp[1:(n_teams - 1)] * params_temp[n_teams + 1:(n_teams - 1)])))
    return(-ll)
  }
  ll_max <- optim(par=rep(0, 2 * n_teams), fn = loglik, control = list(maxit = 100, trace = F), method = "BFGS")
  # Get parameters and scale such that strengths sum up to 0
  att <- ll_max$par[1:n_teams]
  def <- c(ll_max$par[(1 + n_teams):(2 * n_teams - 1)], 0)
  mean_att <- mean(att)
  mean_def <- mean(def)
  params <- matrix(c(att - mean_att, def - mean_def), ncol = 2, dimnames = list(teams, c("Att", "Def"))) #mozna odjac?
  intercept <- mean_att - mean_def
  hta <- ll_max$par[length(ll_max$par)]
  return(structure(list("ratings" = params, "intercept" = intercept, "hta" = hta, "corr" = corr), class = "poisson_model"))
}

estimate_poisson_DC_regularized <- function(games, teams, lambda, alpha) {
  # Function to compute Dixon & Coles model with regularization  
  ###################################
  # INPUT & OUTPUT: see function 'estimate_poisson_DC' + regularization params
  ###################################
  n_teams <- length(teams)
  parameters <- rep(0, 2 * n_teams + 1) # for def/off strength, HTA and rho
  # Log-tau for modelling the correlation of low scores
  log_tau <- function(x, y, lambda, mu, rho, eps = 10^-8) { # Look at DC paper
    val <- ifelse(x > 1 | y > 1, 1, 
                  ifelse(x == 1 & y == 0, 1 + mu * rho, 
                         ifelse(x == 0 & y == 1, 1 + lambda * rho,  
                                ifelse(x == 1 & y == 1, 1 - rho, 1 - lambda * mu * rho))))
    return(log(pmax(val, eps))) # ! log(0) !
  }
  # Log-likeligood function  
  loglik <- function(parameters) {
    params_temp <- c(parameters[1:(2 * n_teams - 1)], 0) # Defensive skills of last team set to 0
    hta <- parameters[length(parameters) - 1]
    rho <- parameters[length(parameters)]
    diff_home <- hta + params_temp[games[,"home_team"]] - params_temp[n_teams + games[,"away_team"]]
    diff_away <- params_temp[games[,"away_team"]] - params_temp[n_teams + games[,"home_team"]]
    rate1 <- exp(diff_home)
    rate2 <- exp(diff_away)
    ll <- (sum(log_tau(games[,"FTHG"], games[,"FTAG"], rate1, rate2, rho) + games[,"FTHG"] * diff_home + games[,"FTAG"] * diff_away 
              - rate1 - rate2) - lambda*((1 - alpha)/ 2 * sum(params_temp^2) + alpha * sum(abs(params_temp))))
    return(-ll)
  }
  ll_max <- optim(par=rep(0, 2 * n_teams + 1), fn = loglik, control = list(maxit = 100, trace = F), method = "BFGS")
  # Get parameters and scale such that strengths sum up to 0
  att <- ll_max$par[1:n_teams]
  def <- c(ll_max$par[(1 + n_teams):(2 * n_teams - 1)], 0)
  mean_att <- mean(att)
  mean_def <- mean(def)
  params <- matrix(c(att - mean_att, def - mean_def), ncol = 2, dimnames = list(teams, c("Att", "Def"))) #mozna odjac?
  intercept <- mean_att - mean_def
  hta <- ll_max$par[length(ll_max$par)-1]
  rho <- ll_max$par[length(ll_max$par)]
  return(structure(list("ratings" = params, "intercept" = intercept, "hta" = hta, "rho" = rho), class = "poisson_DC_model"))
}

estimate_poisson_DC_regularized_prime <- function(games, teams, lambda, alpha, corr) {
  # Function to compute Dixon & Coles model with regularization  
  ###################################
  # INPUT & OUTPUT: see function 'estimate_poisson_DC' + regularization params
  ###################################
  lambda <- lambda / (1 - corr^2)
  n_teams <- length(teams)
  parameters <- rep(0, 2 * n_teams + 1) # for def/off strength, HTA and rho
  # Log-tau for modelling the correlation of low scores
  log_tau <- function(x, y, lambda, mu, rho, eps = 10^-8) { # Look at DC paper
    val <- ifelse(x > 1 | y > 1, 1, 
                  ifelse(x == 1 & y == 0, 1 + mu * rho, 
                         ifelse(x == 0 & y == 1, 1 + lambda * rho,  
                                ifelse(x == 1 & y == 1, 1 - rho, 1 - lambda * mu * rho))))
    return(log(pmax(val, eps))) # ! log(0) !
  }
  # Log-likeligood function  
  loglik <- function(parameters) {
    params_temp <- c(parameters[1:(2 * n_teams - 1)], 0) # Defensive skills of last team set to 0
    hta <- parameters[length(parameters) - 1]
    rho <- parameters[length(parameters)]
    diff_home <- hta + params_temp[games[,"home_team"]] - params_temp[n_teams + games[,"away_team"]]
    diff_away <- params_temp[games[,"away_team"]] - params_temp[n_teams + games[,"home_team"]]
    rate1 <- exp(diff_home)
    rate2 <- exp(diff_away)
    ll <- (sum(log_tau(games[,"FTHG"], games[,"FTAG"], rate1, rate2, rho) + games[,"FTHG"] * diff_home + games[,"FTAG"] * diff_away 
              - rate1 - rate2) - lambda*((1 - alpha)/ 2 * sum(params_temp^2) - corr *  sum(params_temp[1:(n_teams - 1)] * params_temp[n_teams + 1:(n_teams - 1)])))
    return(-ll)
  }
  ll_max <- optim(par=rep(0, 2 * n_teams + 1), fn = loglik, control = list(maxit = 100, trace = F), method = "BFGS")
  # Get parameters and scale such that strengths sum up to 0
  att <- ll_max$par[1:n_teams]
  def <- c(ll_max$par[(1 + n_teams):(2 * n_teams - 1)], 0)
  mean_att <- mean(att)
  mean_def <- mean(def)
  params <- matrix(c(att - mean_att, def - mean_def), ncol = 2, dimnames = list(teams, c("Att", "Def"))) #mozna odjac?
  intercept <- mean_att - mean_def
  hta <- ll_max$par[length(ll_max$par)-1]
  rho <- ll_max$par[length(ll_max$par)]
  return(structure(list("ratings" = params, "intercept" = intercept, "hta" = hta, "rho" = rho, "corr" = corr), class = "poisson_DC_model"))
}

# ---------------------------------------------------

hda_probs <- function(model, ...) {
  #print(list(...))
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
  # Old implementation:
  #prob1 <- dpois(0:20, rate1)
  #prob2 <- dpois(0:20, rate2)
  #if(sum(prob1) < 0.99 || sum(prob2) < 0.99) { # check for too big distrepancy
  #  # Debug info: usually happens when one team has overwhelmingly better rate
  #  print(list("rate1" = rate1, "rate2" = rate2))
  #  print(list("prob1" = sum(prob1), "prob2" = sum(prob2)))
  #  #print(rating_system)
  #  stop("sum(prob1) << 1!")
  #}
  #goals <- prob1 %*% t(prob2)
  # New implementation using 'skellam' library
  win1 <- 1 - pskellam(0, lambda1 = rate1, lambda2 = rate2)
  draw <- dskellam(0, lambda1 = rate1, lambda2 = rate2)
  win2 <- max(c(1 - win1 - draw, 0)) # We take max to avoid rounding errors near 0
  return(c("team1" = win1, "draw" = draw, "team2" = win2))
}

hda_probs.poisson_model_correlated <- hda_probs.poisson_model

hda_probs.poisson_DC_model <- function(rating_system, i, j) {
  stop("Use Skellam distribution for HTA computation.")
  rate1 <- exp(rating_system$intercept + rating_system$hta + rating_system$ratings[i, "Att"] - rating_system$ratings[j, "Def"])
  rate2 <- exp(rating_system$intercept + rating_system$ratings[j, "Att"] - rating_system$ratings[i, "Def"])
  prob1 <- dpois(0:25, rate1)
  prob2 <- dpois(0:25, rate2)
  goals <- prob1 %*% t(prob2)
  rho <- rating_system$rho
  tau <- matrix(c(1 - rate1 * rate2 * rho, 1 + rate1 * rho, 1 + rate2 * rho, 1 - rho), ncol = 2, nrow = 2, byrow = TRUE)
  goals[1:2, 1:2] <- tau * goals[1:2, 1:2]
  if(sum(goals) < 0.99) {# check for too big distrepancy
    print(rate1)
    print(rate2)
    stop("sum(goals) << 1!")
  }
  win1 <- sum(goals[lower.tri(goals)])
  draw <- sum(diag(goals))
  return(c("team1" = win1, "draw" = draw, "team2" = 1 - win1 - draw))
}

### Functions for sampling results (and other related to rating systems)
# REIMPLEMENTED - can be deleted
#sample_result <- function(rating_model, ...) {
#  UseMethod("sample_result", rating_model)
#}
#
#sample_result.olr_model <- function(rating_model, team_i, team_j, results, league_points) {
#  # Match result
#  game_result <- sample(c(1, 0, -1), size = 1, prob = hda_probs.olr_model(rating_model, team_i, team_j))
#  # Store the result and add points for the teams
#  updated <- update_results_basic(results, league_points, game_result)
#  return(list("results" = updated[["results"]], "league_points" = updated[["league_points"]]))
#}
#
#sample_result.poisson_model <- function(rating_model, team_i, team_j, results, league_points) {
#  # Match result
#  game_result <- sample(c(1, 0, -1), size = 1, prob = hda_probs.poisson_model(rating_model, team_i, team_j))
#  # Store the result and add points for the teams
#  updated <- update_results_basic(results, league_points, game_result)
#  return(list("results" = updated[["results"]], "league_points" = updated[["league_points"]]))
#}

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

# Same as get_ratings.poisson_model():
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
