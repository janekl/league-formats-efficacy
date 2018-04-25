# Some auxiliary functions 

parse_cl_arguments <- function(args) {
  # Basic parser for command line arguments like "--parameter=10".
  require("stringi")
  #print(args)
  args_list <- list()
  parse_arg <- function(arg_string, args_list) {
    arg_name <- stri_extract_first_regex(arg_string, "(?<=(--)).+(?=(=))")
    arg_value <- stri_extract_first_regex(arg_string, "(?<=(=)).+$")
    # Convert to numeric if possible
    if(!is.na(suppressWarnings(arg_value_num <- as.numeric(arg_value)))) 
      arg_value <- arg_value_num
    args_list[[arg_name]] <- arg_value
    return(args_list)
  }
  for(arg in args)
    args_list <- parse_arg(arg, args_list)
  #print(args_list)
  return(args_list)
}

write_results <- function(simulation, results_save_folder, specific_result_folder, n, model, league_format, team_strength_fun, ...) {
  ### Write results to .rds files for later processing ###  
  rts <- lapply(simulation, function(x) x[["ratings_rounds"]])
  # Aggregate team strength
  rts <- lapply(rts, team_strength_fun)
  pts <- lapply(simulation, function(x) x[["league_points"]])
  params <- paste0(paste(..., sep = "_"), ".rds")
  saveRDS(rts, file.path(results_save_folder, specific_result_folder, n, model, league_format, "ratings", params))
  saveRDS(pts, file.path(results_save_folder, specific_result_folder, n, model, league_format, "points", params))
  rm("simulation", "rts", "pts")
  gc()
  return(invisible())
}

add_tables <- function(...) {
  ### Auxiliary R function to add tables with 1x2 results
  count_1x2 <- as.table(c("-1" = 0, "0" = 0, "1" = 0))
  for(cnt in list(...)) {
    count_1x2[names(cnt)] <- count_1x2[names(cnt)] + cnt
  }
  return(count_1x2)
}
