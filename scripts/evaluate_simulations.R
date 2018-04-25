# Evaluation
library("stringi")
library("dplyr")
library("xtable")
source("config.R")
source("evaluation_functions.R")

# ----------------------------- F U N C T I O N S -----------------------------
evaluation <- function(result_folder, n_teams, model, team_strength_fun, infinity_number, verbose = T) {
  # Compute csv with results
  results <- list()
  league_formats <- list.dirs(file.path(results_save_folder, result_folder, n_teams, model), recursive = F, full.names = F)
  result_files_pattern <- "0\\.\\d+_-?\\d+(\\.\\d)?\\.rds"
  for(league_format in league_formats) {
    if(verbose) cat("League:", league_format, "\n")
    #dir(file.path("../data/results_temp", league_format, "points"))
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
  # Something strange - values in columns are... lists: Converting to vectors manually
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

# Ranking of tournaments
ranking_tournaments <- function(results, metric, league_formats, aggr_fun) {
  metric_aggr <- numeric(length(league_formats))
  names(metric_aggr) <- league_formats
  total <- 0
  for(sigma_val in unique(results$sigma)) {
    for(shape_val in unique(results$shape)) {
      results_temp <- results %>% filter(sigma == sigma_val, shape == shape_val) %>% arrange(system)
      metric_aggr_temp <- aggr_fun(results_temp[,metric], metric) # rank(ifelse(metric == "spearman_footrule", 1, -1) * results_temp[,metric])
      names(metric_aggr_temp) <- results_temp[,"system"]
      metric_aggr <- metric_aggr + metric_aggr_temp[names(metric_aggr)]
      total <- total + 1
    }
  }
  metric_aggr / total
}

summary_over_settings <- function(results, metrics, league_formats, aggr_fun) {
  ranking <- data.frame(matrix(NA, nrow = length(league_formats), ncol = length(metrics), dimnames = list(league_formats, metrics)))
  for(metric in metrics) {
    ranking[, metric] <- ranking_tournaments(results, metric, league_formats, aggr_fun)
  }
  ranking
}

report_table_rankings <- function(ranking, digits = 3) {
  print(xtable(ranking, align = rep("c", ncol(ranking) + 1), digits = digits, caption = paste0("Aaaa")), label = label, 
        include.rownames = T,
        sanitize.text.function = function(x) stri_replace_all_fixed(x, "_", "-"),
        sanitize.colnames.function = function(x) paste0("\\textbf{", stri_replace_all_fixed(x, "_", "-"), "}"), 
        sanitize.rownames.function = function(x)  stri_replace_all_fixed(x, "_", "-"),
        hline.after = c(0, nrow(ranking)))
}

# Plots of metrics against parameters
plot_metric <- function(results, tournaments, metric, param_x, param_fixed, param_value, ylim, legend_pos, legend_lab, log_x = F, cex.axis = 1.25,
                        xlab = "", ylab = "", xylab_cex = 1, mar = c(3.0, 3.5, 0.5, 1.5), param_x_filter = NA, save_pdf = NULL, save_width = 0, save_height = 0) {
  if(!is.null(save_pdf)) {
    pdf(save_pdf, width = save_width, height = save_height)
    on.exit(dev.off())
  }
  par(mar = mar)
  results2 <- results[results[,param_fixed] == param_value, c("system", param_x, metric)]
  results2 <- results2 %>% filter(system %in% tournaments)
  if(!any(is.na(param_x_filter)))
    results2 <- results2[results2[,param_x] %in% param_x_filter,]
  # Base plot
  tournament <- tournaments[1]
  param_x_values <- sort(unique(results2[,param_x]))
  results2_subset <- results2 %>% filter(system == tournament) %>% arrange_(param_x)
  if(log_x)
    x <- log(param_x_values)
  else
    x <- param_x_values
  plot(x = x, y = results2_subset[,metric], type = "l", ylim = ylim, axes = F, xlab = "", ylab = "")
  i <- 2
  for(tournament in tournaments[-1]) {
    results2_subset <- results2 %>% filter(system == tournament) %>% arrange_(param_x)
    #if(length(x) != length(results2_subset[,metric])) browser()
    lines(x = x, y = results2_subset[,metric], lty = i)
    i <- i + 1
  }
  axis(1, at = x, labels = param_x_values, tck = -0.01, cex.axis = cex.axis, mgp=c(2, 0.3, 0))
  #ticks <- log10(as.numeric(rownames(r)))
  axis(2, at = seq(ylim[1], ylim[2], by = 0.1), labels = seq(ylim[1], ylim[2], by = 0.1), tck = -0.01, cex.axis = cex.axis, mgp = c(2, 0.3, 0), las = 2)
  mtext(side = 1, text = xlab, cex = xylab_cex, line = 2)
  mtext(side = 2, text = ylab, cex = xylab_cex, line = 2.5)
  n <- length(tournaments)
  legend(legend_pos, legend = legend_lab, lty = 1:n, lwd = rep(2, n), bty = "n")
}

# Plots of metrics against number of rounds
regression_line_plot <- function(x, y, legend_pos, axis_y_at, fmt_labels = "%.2f", xlab = "", ylab = "", ylab_pos = 3.0, xylab_cex = 1.25, mar = c(1.5, 2.5, 0.5, 1.5), save_pdf = NULL, save_width = 0, save_height = 0) {
  if(!is.null(save_pdf)) {
    pdf(save_pdf, width = save_width, height = save_height)
    on.exit(dev.off())
  }
  par(mar = mar)
  print(summary(m <- lm(y~log(x))))
  plot(x, y, log = "x", pch = 1, cex = 1.5, axes = F, ylim = range(axis_y_at), xlab = "", ylab = "")
  lines(x, m$coefficients[1] + m$coefficients[2]*log(x))
  legend(legend_pos, lty=1, col=1, sprintf("y = %.2f% + .2f log(x)", m$coefficients[1], m$coefficients[2]), bty = "n", cex=1.25)
  axis(1, at = c(10.04, 20, 50, 110), labels = c(10, 20, 50, 100), tck = -0.01, cex.axis = 1.25, mgp=c(2, 0.3, 0))
  axis(2, at = axis_y_at, labels = sprintf(fmt_labels, axis_y_at), tck = -0.01, cex.axis = 1.25, mgp=c(2, 0.3, 0), las = 2)
  mtext(side = 1, text = xlab, cex = xylab_cex, line = 1.5)
  mtext(side = 2, text = ylab, cex = xylab_cex, line = ylab_pos)
}

# Points allocation
plot_xy <- function(x, y, at_y, xlab = "", ylab = "", ylab_pos = 3.0, xylab_cex = 1.25, labels_y = TRUE, mar = c(3.0, 4.5, 0.5, 1.5), save_pdf = NULL, save_width = 0, save_height = 0) {
  if(!is.null(save_pdf)) {
    pdf(save_pdf, width = save_width, height = save_height)
    on.exit(dev.off())
  }
  par(mar = mar)
  plot(x, y, pch = 1, cex = 1.5, axes = F, ylim = range(at_y), xlab = "", ylab = "")
  #legend(legend_pos, lty=1, col=1, sprintf("y = %.2f% + .2f log(x)", m$coefficients[1], m$coefficients[2]), bty = "n", cex=1.25)
  axis(1, tck = -0.01, cex.axis = 1.25, mgp=c(2, 0.3, 0))
  axis(2, at = at_y, labels = labels_y, tck = -0.01, cex.axis = 1.25, mgp = c(2, 0.3, 0), las = 2)
  mtext(side = 1, text = xlab, cex = xylab_cex, line = 2.0)
  mtext(side = 2, text = ylab, cex = xylab_cex, line = ylab_pos)
}

report_table <- function(results, metrics, sigma_values, shape_values, tournament_order, digits, add_column) {
  results2 <- matrix(nrow = length(tournament_order) * length(shape_values), ncol = length(metrics) * length(sigma_values))
  for(i in seq_along(shape_values)) {
    for(j in seq_along(sigma_values)) {
      results_temp <- results %>% filter(sigma == sigma_values[j], shape == shape_values[i])
      for(k in seq_along(metrics)) {
        values <- results_temp[, metrics[k]]
        names(values) <- results_temp[,"system"]
        results2[1:length(tournament_order) + (i - 1) * length(tournament_order), j + (k-1) * length(sigma_values)] <- values[tournament_order]
      }
    }  
  }
  results2 <- round(results2, 3)
  colnames(results2) <- rep(sigma_values, length(metrics))
  results2 <- data.frame(results2)
  if(!missing(add_column))
    results2 <- cbind(results2, "extra" = add_column)
  
  print(xtable(results2, align = c("l", rep("c", ncol(results2))), digits = digits, caption = paste0("Aaaa")), label = label, 
        include.rownames = T,
        sanitize.text.function = function(x) x, #stri_replace_all_fixed(x, "_", "-"),
        sanitize.colnames.function = function(x) paste0("\\textbf{", x, "}"), 
        sanitize.rownames.function = function(x) "", #stri_replace_all_fixed(x, "_", "-"),
        hline.after = c(0, seq(length(tournament_order), nrow(results2), length.out = length(shape_values))))
}

# -------------------------------------- R E S U L T S -------------------------------------
# kRR for k = 1:10
evaluate_all("results_fixed_104_krr", "poisson_correlated", "12")

# Points allocation
evaluate_all("results_fixed_104_allocation", "poisson_correlated", "12")
evaluate_all("results_fixed_105_allocation", "poisson_correlated", "12")

# Everything and everything plus
evaluate_all("results_fixed_105_all", "poisson_correlated", "12")
evaluate_all("results_fixed_104_all_plus", "poisson_correlated", c("12", "16"))

# Load results

results <- read.csv(file.path(results_save_folder, "results", "results_fixed_105_all_poisson_correlated_12.csv"))
league_formats <- as.character(unique(results$system))

# -------------------------------- Special case and overall analysis --------------------------------
metrics <- c("kendall_tau", "spearman_footrule", "the_best_win")
ranking0 <- summary_over_settings(results[results$sigma == 0.3 & results$shape == 20, ], metrics, league_formats, function(x, metric) x)
ranking0 <- ranking0[order(-ranking0$kendall_tau),]
ranking <- ranking0

ranking <- summary_over_settings(results, metrics, league_formats, function(x, metric) (x - mean(x)) / sd(x))
ranking2 <- summary_over_settings(results, metrics, league_formats, function(x, metric) rank(ifelse(metric == "spearman_footrule", 1, -1) * x))
ranking2 <- round(ranking2[order(ranking2$kendall_tau),], 3)

ranking[order(-ranking$kendall_tau),]

ranking <- ranking[order(-ranking$kendall_tau),]

report_table_rankings(ranking)

# --------------------------------  Joint report on parameter grid --------------------------------
results <- read.csv(file.path(results_save_folder, "results", "results_fixed_105_all_poisson_correlated_12.csv"))

metrics <- c("kendall_tau", "spearman_footrule", "the_best_win")
sigma_values <- sort(unique(results$sigma))[c(1, 4, 7)]
shape_values <- sort(unique(results$shape))[c(1, 2, 4, 7)]
tournament_order <- c("scotland_full", "scotland_half", "3RR", "kazakhstan_full", "kazakhstan_half", "ekstraklasa_full", "ekstraklasa_half", "2RR", "1RR")

report_table(results, metrics, sigma_values, shape_values, tournament_order, digits = c(1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1),
             add_column = rep(paste0("$", c("a_1", "a_2", "b", "c_1", "c_2", "d_1", "d_2", "e", "f"), "$"), 4))


# -------------------------------- Plots of metrics against parameters --------------------------------
tournaments <- c("scotland_full", "3RR", "kazakhstan_full", "ekstraklasa_full", "2RR", "1RR")
legend_lab <- c(expression("a"[1]), "b", expression("c"[1]), expression("d"[1]), "e", "f")

plot_metric(results, tournaments, "kendall_tau", "sigma", "shape", 20, ylim = c(0.2, 0.8), legend_pos = "topleft", legend_lab = legend_lab,
            log_x = F, xlab = expression("Prior ratings" ~ sigma), ylab = expression("Kendall's" ~ tau), mar = c(3, 4.0, 0.5, 1.0), xylab_cex = 1.4, save_pdf = "kendallsigma.pdf", save_width = 6, save_height = 5) 

plot_metric(results, tournaments, "kendall_tau", "shape", "sigma", 0.1, ylim = c(0.2, 0.5), legend_pos = "topright", legend_lab = legend_lab, log_x = T,
            param_x_filter = c(10, 20, 50, 100, 200), xlab = expression("Temporal drift" ~ alpha), xylab_cex = 1.4, mar = c(3, 3.5, 0.5, 1.0), save_pdf = "kendallshape.pdf", save_width = 6, save_height = 5)

# -------------------------------- Contour plots --------------------------------
coutour_plot <- function(results, league_format, metric, shapes_to_filter) {
  results_temp <- results %>% 
    filter(system == league_format, !(shape %in% shapes_to_filter)) %>% 
    dplyr::select(sigma, shape, one_of(metric))
  sigma_values <- sort(unique(results_temp$sigma))
  shape_values <- sort(unique(results_temp$shape))
  r <- matrix(NA, nrow = length(shape_values), ncol = length(sigma_values), dimnames = list(shape_values, sigma_values)) 
  for(sigma in sigma_values) {
    for(shape in shape_values) {
      r[as.character(shape), as.character(sigma)] <- results_temp[results_temp[,"sigma"] == sigma & results_temp[,"shape"] == shape, metric]
    }
  }
  # Plot
  par(mar = c(1.5, 1.5, 0.5, 1.5))
  x <- as.numeric(colnames(r))
  y <- log10(as.numeric(rownames(r)))
  contour(x, y, t(r), axes = F, ylim = range(y, finite = TRUE))
  axis(1, at = x, labels = x, tck = -0.01, cex.axis = 1, mgp=c(2, 0.3, 0))
  ticks <- log10(as.numeric(rownames(r)))
  axis(2, at = ticks, labels = round(10**ticks),
       tck = -0.01, cex.axis = 1, mgp = c(2, 0.3, 0), las = 2)
}
results <- read.csv(file.path(results_save_folder, "results_fixed_104_all_plus_poisson_correlated_12.csv"))
coutour_plot(results, "scotland_full", "kendall_tau", c(50, 100, 200, 500, 1000))
coutour_plot(results, "scotland_full", "spearman_footrule", c(200, 500, 1000))
coutour_plot(results, "scotland_full", "the_best_win", c(200, 500, 1000))


# -------------------------------- RR and metrics --------------------------------
results2 <- read.csv(file.path(results_save_folder, "results", "results_fixed_104_krr_poisson_correlated_12.csv"))
#results2 <- results[(results$shape ==  20) & (results$sigma == 0.3),]
n_games_rr <- 1:10 * 11 # * 6 
names(n_games_rr) <- paste0(1:10, "RR")
results2$n_games <- n_games_rr[as.character(results2$system)]

regression_line_plot(results2$n_games, results2[,"kendall_tau"], "topleft", seq(0.55, 0.85, 0.05), ylab = expression("Kendall's" ~ tau), mar = c(2.5, 4.0, 0.5, 1.5))#, save_pdf = "rrkendall.pdf", save_width = 5, save_height = 5)
regression_line_plot(results2$n_games, results2[,"spearman_footrule"], "topright", seq(0.8, 2., 0.2), xlab = "Rounds", ylab = expression("Spearman's" ~ rho), ylab_pos = 2.5, mar = c(2.5, 4.0, 0.5, 1.5), fmt_labels = "%.1f")#, save_pdf = "rrspearman.pdf", save_width = 5, save_height = 5)
regression_line_plot(results2$n_games, results2[,"the_best_win"], "topleft", seq(0.45, 0.8, 0.05), ylab = "Best team win", mar = c(2.5, 4.0, 0.5, 1.5))#, save_pdf = "rrwin.pdf", save_width = 5, save_height = 5)

# -------------------------------- Points allocation --------------------------------
results2 <- read.csv(file.path(results_save_folder, "results", "results_fixed_105_allocation_poisson_correlated_12.csv"))
plot_xy(seq(1.5, 5, 0.5), results2$kendall_tau, seq(0.722, 0.734, 0.002), ylab = expression("Kendall's" ~ tau), ylab_pos = 3.5, save_pdf = "allocationk.pdf", save_width = 5, save_height = 5)
plot_xy(seq(1.5, 5, 0.5), results2$spearman_footrule, seq(1.24, 1.29, 0.01), xlab = "Points", ylab = expression("Spearman's" ~ rho), ylab_pos = 3.0, save_pdf = "allocations.pdf", save_width = 5, save_height = 5)
plot_xy(seq(1.5, 5, 0.5), results2$the_best_win, seq(0.635, 0.650, 0.005), ylab = "Best team win", ylab_pos = 3.5, save_pdf = "allocationw.pdf", save_width = 5, save_height = 5)


# -------------------------------- Differences between 12 and 16 teams --------------------------------
results12 <- read.csv(file.path(results_save_folder, "results_fixed_104_all_plus_poisson_correlated_12.csv"))
results16 <- read.csv(file.path(results_save_folder, "results_fixed_104_all_plus_poisson_correlated_16.csv"))

summary_over_settings(results12[results12$sigma == 0.3 & results12$shape == 20, ], metrics, league_formats, function(x, metric) x)
summary_over_settings(results16[results16$sigma == 0.3 & results16$shape == 20, ], metrics, league_formats, function(x, metric) x)
