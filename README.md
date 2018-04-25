# The Efficacy of League Formats in Ranking Teams

This repository contains the code used for the analysis of football data
in the paper *The efficacy of league formats in ranking teams*
(submitted to the special issue on sports analytics in
*Statistical Modelling Journal*) by [Jan Lasek](http://lasek.rexamine.com)
and [Marek Gagolewski](http://gagolewski.com/).
The efficacy of a given league format is understood as the accuracy to
reproduce the true ranking of competing teams. The content of this repository,
steps to reproduce the main results and supplementary materials are outlined below.

## Contents

```
MatchPrediction
├── data
│    ├── Poland{1415,1516}.csv               # Data from http://www.90minut.pl/
│    ├── download_data.sh                    # Script for downloading data from http://www.football-data.co.uk/
│    └── preprocessing_football_data_co_uk.R # For processing raw CSV files from http://www.football-data.co.uk/
└── scripts
     ├── config.R
     ├── setup_simulations.R
     ├── run_simulations.R
     ├── run_all_final_104.sh                # 10^4 simulations for all settings for Poisson correlated model
     ├── run_all_final_105.sh                # 10^5 simulations for chosen (most important for analysis) settings
     ├── run_all_olr.sh                      # Simulations for OLR model
     ├── run_all_poisson.sh                  # Base Poisson model
     ├── run_simulations_krr.sh              # Round impact in the kRR tournament
     ├── run_simulations_scotland_points_allocation.sh # Points allocations
     ├── schedule_functions.R
     ├── simulation_functions.R
     ├── evaluation_functions.R
     ├── extra_functions.R
     ├── evaluate_simulations.R
     └── rating_systems
           ├── rating_systems.R
           ├── prediction_functions.R
           └── optimise_models.R
```

### Requirements

The project was developed in [R](https://www.r-project.org/) using
the following packages:

* doMC_1.3.5
* dplyr_0.5.0
* iterators_1.0.8
* foreach_1.4.3
* skellam_0.2.0
* stringi_1.1.7
* xtable_1.8-2

Compatibility of the scripts for other versions of the aforementioned packages
is not guaranteed.



## Reproducing the Main Results

### Step 1

First, create file **scripts/config.R** and set appropriate paths
for keeping data and saving the results. This depends on local platform
settings. File contents should look like:

```R
results_save_folder <- "where/you/want/to/store/data/and/results"
```

All the results will be saved in **results_save_folder**.
This is also the default path for storing data.



### Step 2

Second, to install necessary libraries and setup the
output data folders, run:

```
$ Rscript setup_simulations.R
```

Additionally, you can edit the script and define the following:

* the team rating model to use (`model` variable),
* the league formats to analyse (`x` variable),
* directory for storing the results of a specific experiment
(`specific_result_folder`).


An error will be raised (`directory not empty`), if there are some old
results stored. Clean them first (or move to some other location).


### Step 3

Third, to produce the simulation results with given parameters,
edit the **run_all.sh** script and execute:

```
$ ./run_all.sh
```

Alternatively, to produce the results for a given parameter setup, run:
```
$ Rscript run_simulations.R --n=16 --model=poisson_correlated --n_sim=100 --shape=20 --sigma=0.3 --n_cores=3 --log2file=1
```

with appropriate parameters (please consult the script). These operations
are performed for a parameters grid in the **run_all.sh** script. The results
will be saved in the **data** folder with appropriate paths as
specified in the **config.R** and **setup_simulations.R** scripts.

Results are saved to **/mnt/ml-team/experiments/Leagues**
and their backup copy is at **ml100:~/Leagues/results_backup** (29-01-2018).

## Reproducing Intermediate Results

### Parameter Setting for the Model

This is specified in the **scripts/rating_systems/parameter_setting_for_paper.sh**
script. Make sure to set ρ=0.45 (or as desired) for the
correlated Poisson model in the **rating_systems/prediction_functions.R**
script. To reproduce results, run

```
$ ./parameter_setting_for_paper.sh
```

in the appropriate folder. The results will be saved in the **results/** folder.
For a complete set of tables, parameters, results and a comparison
with the odds model, consult **prediction_results_report.R**
(interactively in, e.g., RStudio). You will need to transfer the estimates
for the λ there.

## Appendices

### League Formats in UEFA

The first appendix consists of a listing of formats that are in
operation in the UEFA countries in the 2017/2018 (or 2018) season.

Click to download:
[PDF](Appendix_League_Formats_in_UEFA.pdf) |
[MARKDOWN](Appendix_League_Formats_in_UEFA.md)


### Tournament Metrics for Several Parameter Combinations

The second appendix gives the detailed estimates of different tournament
metrics considered in our study.

Click to download:
[PDF](Appendix_Tournament_Metrics.pdf) |
[MARKDOWN](Appendix_Tournament_Metrics.md)


### Schedules in Two Stage Systems

Two example schedules of the final round in the *2RR + (1RR/1RR)* league format
employed in the championship and the relegation group in the case of 12 and 16
teams are listed in the third appendix.

Click to download:
[PDF](Appendix_Schedules_in_Two_Stage_Systems.pdf) |
[MARKDOWN](Appendix_Schedules_in_Two_Stage_Systems.md)
