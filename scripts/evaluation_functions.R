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

# Evaluation metrics

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
