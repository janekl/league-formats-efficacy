prepare_schedule_1RRn <- function(n) {
  # Single round robin schedule generation.
  # Algorithm: D. de Werra: "Scheduling in Sports" (1981) - almost https://en.wikipedia.org/wiki/Round-robin_tournament
  # NOTE: this algorithm produces 2 "conflicts" from round to round, that is team that repeteadly plays home or away game
  # (1) Initialize first round
  # (2) Rotate teams clockwise
  # (3) Balance home and away rounds by flipping 1st match (from top) in even rounds
  # (4) Introduce flips in even matches in each round to finally balance the schedule
  stopifnot((n %% 2) == 0, n >= 4) # This algorithm works for even n and for sure for n >= 6. TODO:check n = 2 and n = 4
  n2 <- n / 2 # Number of matches per round
  # (1) Initialize first round
  single_rr <- data.frame("Round" = rep(1:(n - 1), each = n2), "Team1" = NA, "Team2" = NA)
  single_rr[1:n2, c("Team1", "Team2")] <- matrix(c(1:n2, n:(n2+1)), ncol = 2, byrow = FALSE)
  n_rounds <- n - 1
  rounds <- vector("list", n_rounds)
  for(r in 1:n_rounds) # For fast access to matches on a given round
    rounds[[r]] <- 1:n2 + n2 * (r - 1) # == which(single_rr$Round == r)
  # (2) Rotate teams clockwise
  for(r in 2:(n - 1)) {
    home_team_prev <- single_rr[rounds[[r - 1]], "Team1"]
    away_team_prev <- single_rr[rounds[[r - 1]], "Team2"]
    home_team_next <- c(home_team_prev[2:n2], away_team_prev[n2])
    away_team_next <- c(away_team_prev[1], home_team_prev[1], away_team_prev[2:(n2 - 1)])
    single_rr[rounds[[r]], c("Team1", "Team2")] <- c(home_team_next, away_team_next)
  }
  # (3) Balance home and away rounds by flipping 1st match (from top) in even rounds
  for(r in seq(2, n - 1, by = 2)) {
    index <- rounds[[r]]
    single_rr[index[1], c("Team1", "Team2")] <- single_rr[index[1], c("Team2", "Team1")]
  }
  # (4) Introduce flips in even matches in each round to finally balance the schedule
  if(n2 >= 3) { # Does it work for n== 4? DOUBLE CHECK!
    for(r in 1:(n - 1)) {
      index <- rounds[[r]][seq(2, n2 - 1, by = 2)]
      single_rr[index, c("Team1", "Team2")] <- single_rr[index, c("Team2", "Team1")]
    }
  }
  schedule <- list()
  schedule[["schedule"]] <- single_rr
  schedule[["rounds"]] <- rounds 
  return(schedule)
}

prepare_schedule_1RRn_OLD <- function(n) {
  stop("DEPRECATED: use prepare_schedule_1RRn()")
  # 1RR schedule
  # Algorithm: https://en.wikipedia.org/wiki/Round-robin_tournament [original citation needed]
  # Remark: this algorithm produces 2 "conflicts" from round to round, that is team that repeteadly plays home or away game
  # (1) Initialize first round
  # (2) Rotate teams counterclockwise
  # (3) Balance home and away rounds by flipping teams in even rounds
  # (4) Introduce flips for some teams to balance schedule
  # stopifnot((n %% 4) == 0) # This algorithm works only for special case n = 4m
  # This is connected to the flips introduced - they are applicable only in this case
  # Otherwise do not expect balanced schedule
  n2 <- n / 2
  # (1) Initialize first round
  single_rr <- data.frame("Round" = rep(1:(n - 1), each = n2), "Team1" = NA, "Team2" = NA)
  single_rr[1:n2, c("Team1", "Team2")] <- matrix(1:n, ncol = 2, byrow = FALSE)
  n_rounds <- n - 1
  rounds <- vector("list", n_rounds)
  for(r in 1:n_rounds) # For fast access to matches on a given round
    rounds[[r]] <- which(single_rr$Round == r)
  
  # (2) Rotate teams "counterclockwise"
  for(r in 2:(n - 1)) {
    home_team_prev <- single_rr[rounds[[r - 1]], "Team1"]
    away_team_prev <- single_rr[rounds[[r - 1]], "Team2"]
    home_team_next <- c(home_team_prev[1], away_team_prev[1], home_team_prev[2:(n2 - 1)])
    away_team_next <- c(away_team_prev[2:n2], home_team_prev[n2])
    single_rr[rounds[[r]], c("Team1", "Team2")] <- c(home_team_next, away_team_next)
  }
  # (3)
  for(r in seq(2, n - 1, by = 2)) {
    index <- rounds[[r]]
    single_rr[index, c("Team1", "Team2")] <- single_rr[index, c("Team2", "Team1")]
  }
  # (4) Introduce flips for some teams to balance schedule
  for(r in seq(n2 + 1, n - 1, by = 2)) {
    index <- rounds[[r]][n2]
    single_rr[index, c("Team1", "Team2")] <- single_rr[index, c("Team2", "Team1")]
  }
  schedule <- list()
  schedule[["schedule"]] <- single_rr
  schedule[["rounds"]] <- rounds 
  return(schedule)
}

prepare_schedule_kRRn <- function(k, n) {
  # kRR schedule
  single_rr <- prepare_schedule_1RRn(n)[["schedule"]]
  if(k > 1) {
    full_rr <- single_rr
    for(i in 2:k) {
      next_rr <- single_rr
      if(i %% 2 == 0) {
        # Flip home and away (2nd leg)
        next_rr[, c("Team1", "Team2")] <- next_rr[, c("Team2", "Team1")]
      }
      full_rr <- rbind(full_rr, next_rr)
    }
  } else {
    full_rr <- single_rr
  }
  full_rr[,"Round"] <- rep(1:(k * (n - 1)), each = n / 2)
  schedule <- list()
  schedule[["schedule"]] <- full_rr
  n_rounds <- k * (n - 1)
  rounds <- vector("list", n_rounds)
  for(r in 1:n_rounds) # For fast access to matches on a given round
    rounds[[r]] <- which(full_rr$Round == r)
  schedule[["rounds"]] <- rounds 
  return(schedule)
}

prepare_schedule_poland <- function(n) {
  stopifnot(n %in% c(12, 16))
  # 1st stage
  first_stage <- prepare_schedule_kRRn(2, n)
  # 2nd stage
  n_rounds_final <- n / 2 - 1
  if(n == 12) {
    schedule_final <- matrix(c(6,1, 2,5, 3,4,  
                               5,4, 6,2, 1,3,
                               1,5, 3,6, 2,4,
                               4,1, 2,3, 5,6,
                               3,5, 4,6, 1,2), ncol = 2, byrow = TRUE)
  }
  if(n == 16) {
    schedule_final <- matrix(c(1,6,  2,5,  3,8,  4,7,
                               8,1,  5,4,  6,2,  7,3,
                               1,5,  4,8,  2,7,  3,6,
                               7,1,  3,5,  8,2,  4,6,
                               1,3,  5,7,  2,4,  6,8,
                               4,1,  2,3,  6,7,  8,5,
                               1,2,  3,4,  5,6,  7,8), ncol = 2, byrow = TRUE)
  }
  
  schedule_final <- data.frame(rep(1:n_rounds_final, each = n / 4), schedule_final)
  colnames(schedule_final) <- c("Round", "Team1", "Team2")
  rounds_final <- vector("list", n_rounds_final)
  for(r in 1:n_rounds_final) # For fast access to matches on a given round
    rounds_final[[r]] <- which(schedule_final$Round == r)
  return(list("schedule_1" = first_stage[["schedule"]],
              "rounds_1" = first_stage[["rounds"]],
              "schedule_2" = schedule_final,
              "rounds_2" = rounds_final))
}

prepare_schedule_kazakhstan <- function(n) {
  stopifnot(n %% 2 == 0)
  # 1st stage
  first_stage <- prepare_schedule_kRRn(2, n)
  # 2nd stage
  #n_rounds_final <- n / 2 - 1
  second_stage <- prepare_schedule_kRRn(2, n / 2)
  return(list("schedule_1" = first_stage[["schedule"]],
              "rounds_1" = first_stage[["rounds"]],
              "schedule_2" = second_stage[["schedule"]],
              "rounds_2" = second_stage[["rounds"]]))
}

prepare_schedule_scotland <- function(n) {
  stopifnot(n %% 2 == 0)
  # 1st stage
  first_stage <- prepare_schedule_kRRn(3, n)
  #if(runif(1) < 0.5) {
  #  # Rename teams: shift by 1. I wanted to do this to make totally independent
  #  # the number of games at gome and the team "name" (index). Apparently, this 
  #  # doesn't help - well, it is clear since the name doesn't influece distribution.  
  #  schedule <- first_stage$schedule
  #  schedule$Team1 <- schedule$Team1 %% n + 1
  #  schedule$Team2 <- schedule$Team2 %% n + 1
  #  first_stage$schedule <- schedule
  #}
  # 2nd stage
  #n_rounds_final <- n / 2 - 1
  #second_stage <- prepare_schedule_kRR(n / 2, 1)
  # Look at 'rounds_2' - temporal solution for Scotland, schedule for final round 
  # is determined later (in particular, variables 'schedule_2' and 'rounds_2')
  # based on results of 1st round. But the number of rounds is known in advance: this is 
  # a hack, rounds_2 is overwritten by simulation later anyway.
  return(list("schedule_1" = first_stage[["schedule"]],
              "rounds_1" = first_stage[["rounds"]],
              "schedule_2" = NULL, # second_stage[["schedule"]],
              "rounds_2" = 1:(n/2 - 1))) # second_stage[["rounds"]]))
}

prepare_schedule <- function(games){
  warning("Function prepare_schedule() is deprecated - use prepare_schedule_ekstraklasa() instead.")
  # Preprocessing data for simulation
  # All games
  teams <- sort(unique(c(games$Team1, games$Team2)))
  schedule <- data.frame("Round" = games$Round, 
                         #"Day" = games$Date_index, 
                         "Team1" = as.numeric(factor(games$Team1, levels = teams)), 
                         "Team2" = as.numeric(factor(games$Team2, levels = teams)))
  # Runda zasadnicza
  rounds <- vector("list", 30)
  for(r in 1:30) # For fast access to matches on a given round
    rounds[[r]] <- which(schedule$Round == r)
  # Runda finalowa
  schedule_final <- matrix(c(1,6,  2,5,  3,8,  4,7,
                             8,1,  5,4,  6,2,  7,3,
                             1,5,  4,8,  2,7,  3,6,
                             7,1,  3,5,  8,2,  4,6,
                             1,3,  5,7,  2,4,  6,8,
                             4,1,  2,3,  6,7,  8,5,
                             1,2,  3,4,  5,6,  7,8), ncol=2, byrow=TRUE)
  schedule_final <- data.frame(rep(1:7, each=4), schedule_final)
  colnames(schedule_final) <- c("Round", "Team1", "Team2")
  rounds_final <- vector("list", 7)
  for(r in 1:7) # For fast access to matches on a given round
    rounds_final[[r]] <- which(schedule_final$Round == r)
  return(list("schedule_1" = schedule,
              "rounds_1" = rounds,
              "schedule_2" = schedule_final,
              "rounds_2" = rounds_final))
}
