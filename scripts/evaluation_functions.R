# Evaluation

the_best_win <- function(points_team, team_strength, compute_mean = TRUE) {
  max_points_team <- sapply(points_team, which.max)
  max_rating_team <- sapply(team_strength, which.max)
  if(compute_mean)
    return(mean(max_points_team == max_rating_team))
  else
    return(max_points_team == max_rating_team)
}

kendall_tau0 <- function(x, y) cor(x, y, method = "kendall")
spearman_footrule0 <- function(x, y) mean(abs(rank(x, ties.method = "random") - rank(y, ties.method = "random")))

kendall_tau <- function(points_team, team_strength, compute_mean = TRUE) {
  corr_all <- numeric(length(points_team))
  for(i in 1:length(points_team)) {
    corr_all[i] <- cor(points_team[[i]], team_strength[[i]], method = "kendall")
  }
  if(compute_mean)
    return(mean(corr_all))
  else
    return(corr_all)
}

spearman_footrule <- function(points_team, team_strength, compute_mean = TRUE) {
  rank_points <- sapply(points_team, rank)
  rank_ratings <- sapply(team_strength, rank)
  rank_diff <- abs(rank_points - rank_ratings)
  if(compute_mean)
    return(sum(rank_diff) / length(points_team) / length(points_team[[1]])) # mean?
  else
    return(apply(rank_diff, 2, sum) / length(points_team[[1]]))
}
