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

# Functions to preprocess raw data downloaded from http://www.football-data.co.uk/.
# There are two main functions for preprocessing:
# * preprocess() - for preprocessing main leagues given on the website saved - single file per season 'G_1314.csv'
# * preprocess_new_leagues() - for minor leagues: single file for multiple seasons like 'POL1217.csv'
# There are some differences between these two files including: header names (passed as parameter to a function)

preprocess_odds <- function(bookies) {
  n_bookmakers <- ncol(bookies) / 3
  cat("No of bookies:", n_bookmakers, "\n")
  bookies <- 1/bookies
  hda <- matrix(0, nrow = nrow(bookies), ncol = 3)
  for(i in 1:n_bookmakers) {
    bookie <- bookies[,(i-1) * 3 + 1:3]
    normalize_rows <- apply(bookie, 1, sum)
    bookie <- bookie / normalize_rows
    bookie[is.na(bookie)] <- 0
    hda <- hda + bookie
  }
  # Average (might be missing data)
  hda <- hda / apply(bookies, 1, function(x) sum(!is.na(x))/3)
  colnames(hda) <- c("H", "D", "A")
  hda
}

determine_rounds <- function(dates, max_games_per_round, same_round_within_days = 1, minimal_number_of_games_in_round = 3) {
  rounds <- numeric(length(dates))
  rounds[1] <- 1
  dates_diff <- diff(dates)
  current_round <- 1
  n_games_round <- 1
  # Assign rounds assuming that the gap can be within `same_round_within_days` days apart
  for(i in 2:length(rounds)) {
    if((n_games_round <= max_games_per_round) & (dates_diff[i-1] <= same_round_within_days)) {
      rounds[i] <- current_round
      n_games_round <- n_games_round + 1
    } else {
      current_round <- current_round + 1
      rounds[i] <- current_round
      n_games_round <- 1
    }
  }
  # Map "small" rounds (1 - 2 postponed matches) to previous round
  rounds_rle <- rle(rounds)
  if(any(rounds_rle$lengths < minimal_number_of_games_in_round)) {
    for(i in seq_along(rounds_rle$lengths)) {
      if(rounds_rle$lengths[i] < minimal_number_of_games_in_round) {
        # Merge round with the previous one (for the first round - with the second)
        if(i > 1) {
          rounds_rle$values[i] <- rounds_rle$values[i-1]
          if(i < length(rounds_rle$lengths)) {
            # Decrease round id by 1
            next_round_indexes <- (i+1):length(rounds_rle$lengths)
            rounds_rle$values[next_round_indexes] <- rounds_rle$values[next_round_indexes] - 1
          }
        } else {
          rounds_rle$values <- rounds_rle$values - 1
          rounds_rle$values[i] <- rounds_rle$values[i+1]
        }
      }
    }
  }
  inverse.rle(rounds_rle)
}

filter_empty_lines <- function(dataframe, possible_non_missing_count = 3) {
  # Function removes rows which contains almost either empty strings or NA values.
  # There is a possibility that `possible_non_missing_count` values will be present
  # but the row is ignored anyway.
  rows_to_remove <- apply(dataframe, 1, function(x) (sum(is.na(x)) + sum(x == "", na.rm = T)) >= (ncol(dataframe) - possible_non_missing_count))
  cols_to_remove <- apply(dataframe, 2, function(x) (sum(is.na(x)) + sum(x == "", na.rm = T)) == nrow(dataframe))
  dataframe[!rows_to_remove,!cols_to_remove]
}

determine_bookies_cols <- function(games) {
  browser()
}

preprocess <- function(games, home_team = "HomeTeam", away_team = "AwayTeam", home_goals = "FTHG", away_goals = "FTAG", 
                       bookies_cols = 'B365H:WHA', date_format = '%d/%m/%Y') {
  require("dplyr")
  games <- filter_empty_lines(games)
  games$Date <- as.Date(games$Date, format = date_format)
  games <- games[order(games$Date),]
  n_teams <- length(unique(games[,home_team]))
  cat("No of teams:", n_teams, "\n")
  Round <- determine_rounds(games$Date, floor(n_teams / 2))
  games <- cbind(Round, games)
  # Get oods, odds to probability and average over bookies
  # TODO:
  #bookies_cols <- determine_bookies_cols(games)
  #bookies <- select_(games, bookies_cols)
  #hda <- preprocess_odds(bookies)
  # Final data selection
  #games <- cbind(games[,c("Date", "Round", home_team, away_team, home_goals, away_goals)], hda)
  games <- cbind(games[,c("Date", "Round", home_team, away_team, home_goals, away_goals)])
  games
}

filter_playoff_teams <- function(games, home_team, away_team, limit = 4) {
  n_games_per_team <- table(unlist(games[,c(home_team, away_team)]))
  n_games_per_team <- n_games_per_team[n_games_per_team <= limit]
  games[!(games[,home_team] %in% names(n_games_per_team) | games[,away_team] %in% names(n_games_per_team)),]
}

preprocess_new_leagues <- function(football_data_co_uk_csv, seasons, filename2country,
              home_team = "Home", away_team = "Away", home_goals = "HG", away_goals = "AG", 
              rename_columns = c("Home" = "HomeTeam", "Away" = "AwayTeam", "HG" = "FTHG", "AG" = "FTAG")) {
  # Input:
  # * tournament_types        - either a vector of length 1 or a mapping like c(season = tournament) type in a vector/list
  # * home_team, away_team    - column names for home and away team, respectively
  # * home_goals, away_goals  - column names for home and away team goals, respectively
  library("dplyr")
  games_all <- read.csv(football_data_co_uk_csv, header = TRUE, stringsAsFactors = FALSE)
  for(season in seasons) {
    games <- games_all %>% filter(Season == season | Season == substr(season, 1, 4))
    games <- filter_playoff_teams(games, home_team, away_team)
    games <- preprocess(games, home_team = home_team, away_team = away_team, home_goals = home_goals, away_goals = away_goals, bookies_cols = 'AvgH:AvgA', date_format = '%d/%m/%Y')
    fileout <- paste0(filename2country[football_data_co_uk_csv], paste0(unlist(strsplit(season, ""))[c(3,4,8,9)], collapse = ""), '.csv')
    for(col in names(rename_columns))
      colnames(games)[colnames(games) == col] <- rename_columns[col]
    write.csv(games, fileout, row.names = F)
  }
}

# Preprocessing first group of leagues
league_long_name <- c(
  D  = "Germany",
  E  = "England",
  F  = "France",
  I  = "Italy",
  SC = "Scotland",
  SP = "Spain"
)
 
for(start_year in seq(93, 116)) {
  for(league in names(league_long_name)) {
    print(c(start_year, league))
    n <- sprintf('%02d', start_year %% 100)
    m <- sprintf('%02d', (start_year + 1) %% 100)
    season <- paste0(n, m)
    fileout <- paste0(league_long_name[league], season, ".csv")
    file_raw <- paste0(league, "_", season, ".csv")
    if(!file.exists(file_raw)) 
      next
    games_raw <- read.csv(file_raw, header = TRUE, stringsAsFactors = FALSE)
    games <- preprocess(games_raw)
    write.csv(games, fileout, row.names = F)
  }
}

# Preprocessing other leagues
filename2country <- c(
  'POL.csv' = 'Poland'
)
seasons <- c('2012/2013', '2013/2014', '2014/2015', '2015/2016', '2016/2017')

for(f in names(filename2country)) {
  print(f)
  preprocess_new_leagues(f, seasons, filename2country)
}
