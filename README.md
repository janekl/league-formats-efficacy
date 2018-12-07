# The Efficacy of League Formats in Ranking Teams

This repository contains the code used for the analysis of football data
in the paper *The efficacy of league formats in ranking teams*
that appeared in the [special issue on sports analytics in
*Statistical Modelling Journal](https://journals.sagepub.com/doi/abs/10.1177/1471082X18798426)*) by [Jan Lasek](http://lasek.rexamine.com)
and [Marek Gagolewski](http://gagolewski.com/).
The efficacy of a given league format is defined as the accuracy to
reproduce the true ranking of competing teams. The content of this repository,
steps to reproduce the main results and supplementary materials are outlined below.

## Contents

```
league-formats-efficacy
├── data
│    ├── download_data.sh # Script for downloading data from http://www.football-data.co.uk/
│    └── preprocessing_football_data_co_uk.R
├── scripts
│    ├── config.R # It should be defined locally
│    ├── setup_simulations.R
│    ├── run_simulations.R
│    ├── run_all.sh # Simulations for all settings
│    ├── schedule_functions.R
│    ├── simulation_functions.R
│    ├── evaluation_functions.R
│    ├── extra_functions.R
│    ├── evaluate_simulations.R
│    └── rating_systems
│         ├── rating_systems.R
│         ├── prediction_functions.R
│         └── optimise_models.R
└── appendix
     ├── League_Formats_in_UEFA(.md|.pdf)
     ├── Schedules_in_Two_Stage_Systems(.md|.pdf)
     └── Tournament_Metrics(.md|.pdf)
```

### Requirements

The project was developed in [R](https://www.r-project.org/) using
the following packages:

* doMC_1.3.5
* iterators_1.0.8
* foreach_1.4.3
* dplyr_0.5.0
* stringi_1.1.7
* skellam_0.2.0
* MASS_7.3-45

Compatibility of the scripts for other versions of the aforementioned packages
is not guaranteed.



## Reproducing the Main Results

### Step 1

First, create file **scripts/config.R** and set appropriate paths
for keeping data and saving the results. This depends on local platform
settings. File contents should look like:

```R
results_save_folder <- "where/you/want/to/store/simulation/results"
```

All the results will be saved in **results_save_folder**.



### Step 2

Second, to install necessary libraries and setup 
the output data folders, run:

```
$ Rscript setup_simulations.R
```

Additionally, you can edit the script and define the following:

* the team rating model to use (`model` variable),
* the league formats to analyse (`league_formats`),
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
$ Rscript run_simulations.R --n=12 --model=poisson_correlated --n_sim=100 --shape=20 --sigma=0.3 --n_cores=3 --specific_result_folder=results_all --log2file=1
```

with appropriate parameters (please consult the script). These operations
are performed for a parameters grid in the **run_all.sh** script. The results
will be saved in the folder specified in the **config.R** script.

### Step 4

Finally, for measuring the agreement between the final league standings 
and latent team strength parameters, execute:

```
$ Rscript evaluate_simulations.R
```

This produces a *csv* file with results in the respective folder.


## Reproducing Intermediate Results

### Rating Systems Performance

In order to generate and evaluate the predictions of different rating systems 
first download data from http://www.football-data.co.uk/.
There is script **data/download_data.sh** to assist you with it.
The data need to be first preprocessed by running

```
$ Rscript preprocessing_football_data_co_uk.R
```

in the **data** folder.

The correlation parameter is set to ρ=0.45 (or as desired) for 
the correlated Poisson model in the **scripts/rating_systems/prediction_functions.R** script. 
To reproduce results, go to **scripts/rating_systems** and run

```
$ Rscript optimize_models.R
```
The script runs grid search for regularization parameter λ for different models 
and a league season of choice (this can be specified directly in the script) and 
saves some logs and results into **results** folder in the same directory.

## Appendices

### League Formats in UEFA

The first appendix consists of a listing of formats that are in
operation in the UEFA countries in the 2017/2018 (or 2018) season.

Click to download:
[PDF](appendix/League_Formats_in_UEFA.pdf) |
[MARKDOWN](appendix/League_Formats_in_UEFA.md)


### Tournament Metrics for Several Parameter Combinations

The second appendix gives the detailed estimates of different tournament
metrics considered in our study.

Click to download:
[PDF](appendix/Tournament_Metrics.pdf) |
[MARKDOWN](appendix/Tournament_Metrics.md)


### Schedules in Two Stage Systems

Two example schedules of the final round in the *2RR + (1RR/1RR)* league format
employed in the championship and the relegation group in the case of 12 and 16
teams are listed in the third appendix.

Click to download:
[PDF](appendix/Schedules_in_Two_Stage_Systems.pdf) |
[MARKDOWN](appendix/Schedules_in_Two_Stage_Systems.md)
