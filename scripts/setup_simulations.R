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

# Setup necessary packages, directories to write results, clean up old resutls (if desired)
if(!require("doMC"))
  install.packages("doMC")

if(!require("stringi"))
  install.packages("stringi")

if(!require("skellam"))
  install.packages("skellam")

################## EDIT THIS PART ##################
# And the config.R file for setting general path to store results (`results_save_folder` variable)
source("config.R")
specific_result_folder <- "results_all"
models <- "poisson_correlated" # c("olr", "poisson", "poisson_correlated")
n_teams <- c(12, 16)

league_formats <- c(
  c("1RR", "2RR", "3RR", "4RR"), # Round-robins
  paste0(rep(c("poland", "kazakhstan", "scotland"), each = 2), c("_half", "_full"))
)

# Some more tests and results
#league_formats <- paste0(1:10, "RR") # Influence of additional round-robin rounds
#league_formats <- paste0("scotland_full_", seq(1.5, 5, 0.5), "_", 1) # Special analysis for different point allocation
#league_formats <- c("4RR") # Additinal test
####################################################

dir.create(results_save_folder, recursive = T, showWarnings = F)
dir.create(file.path(results_save_folder, "logs"), recursive = F, showWarnings = F)

setups <- expand.grid(n = n_teams, model = models, league_format = league_formats, data_type = c("points", "ratings"))
for(i in seq_along(rownames(setups))) {
  dir_path <- file.path(results_save_folder, specific_result_folder, setups$n[i], setups$model[i], setups$league_format[i], setups$data_type[i])
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  if(length(list.files(dir_path))) {
    stop(paste0("Result directory is not empty: ", dir_path))
  }
}
