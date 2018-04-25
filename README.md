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
(PDF)[Appendix_League_formats_in_UEFA.pdf] |
(MARKDOWN)[Appendix_League_formats_in_UEFA.md]


### Tournament Metrics for Several Parameter Combinations

Table below presents the detailed estimates of different tournament metrics considered in this study.
Each entry in the table corresponds to a different pair of parameters $`(\alpha, \sigma)`$ given in the rows and columns, respectively.
The table is organised in blocks.
Each block corresponds to nine tournament designs studied. **They are presented in the same order as given in Tab. \ref{LeaguesUnderStudy} --
the last column gives the short name for the particular league design scheme.**

We suggest that the significance of differences between different formats is compared based
on the confidence intervals resulting from the normal approximation. That is,
the $1-\bar{\alpha}$ confidence interval for a sample of observations $`\mathbf{x} = (x_1, x_2, \dots, x_n)`$ is $`\bar{\mathbf{x}} \pm \frac{z_{1-\bar{\alpha}/2} \cdot sd(\mathbf{x})}{\sqrt{n}}`$,

$`
\left(\bar{\mathbf{x}} - \frac{z_{1-\bar{\alpha}/2} \cdot sd(\mathbf{x})}{\sqrt{n}},
\bar{\mathbf{x}} + \frac{z_{1-\bar{\alpha}/2} \cdot sd(\mathbf{x})}{\sqrt{n}}\right),
`$

where $`z`$ denotes a given quantile of the standard normal distribution and
$`\bar{\mathbf{x}}`$ and $`sd(\mathbf{x})`$ are the sample mean and the sample standard deviation, respectively.
The width of this interval is $`\frac{2}{\sqrt{n}} z_{1-\bar{\alpha}/2} \cdot sd(\mathbf{x})`$.
Assuming the significance level of $`\bar{\alpha} = 0.05`$, we suggest the three given metrics
considered -- Kendall's $\tau$, Spearman's Footrule distance and the fraction of the best team wins -- should be considered with error margins of ca. $`\pm 0.001`$, $`\pm 0.004`$ and $`\pm 0.003`$, respectively.
Accordingly, the metric values are rounded to the third decimal place.
These margins should be taken into account when considering the significance of differences between the reported numbers.

|           |    <-     |  Kendall  |     ->    |     <-    |  Spearman |     ->    |     <-    | Best wins |    ->     |  Format   |
|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|
| $`\alpha \setminus \sigma`$  |  **0.1**  | **0.25**  |  **0.4**  |  **0.1**  | **0.25**  |  **0.4**  |  **0.1**  | **0.25**  |  **0.4**  |        |
|**10**     |0.463      |0.692      |0.789      |2.319      |1.409      |1.008      |0.383      |0.603      |0.704      |$`a_1`$    |
|**10**     |0.458      |0.685      |0.781      |2.338      |1.442      |1.039      |0.376      |0.589      |0.691      |$`a_2`$    |
|**10**     |0.445      |0.676      |0.775      |2.386      |1.477      |1.067      |0.360      |0.578      |0.682      |$`b`$      |
|**10**     |0.421      |0.662      |0.768      |2.484      |1.535      |1.095      |0.363      |0.585      |0.691      |$`c_1`$    |
|**10**     |0.418      |0.654      |0.760      |2.495      |1.565      |1.127      |0.360      |0.568      |0.677      |$`c_2`$    |
|**10**     |0.405      |0.644      |0.753      |2.540      |1.606      |1.159      |0.337      |0.553      |0.666      |$`d_1`$    |
|**10**     |0.403      |0.638      |0.747      |2.548      |1.630      |1.184      |0.331      |0.539      |0.655      |$`d_2`$    |
|**10**     |0.382      |0.618      |0.731      |2.625      |1.712      |1.253      |0.307      |0.517      |0.631      |$`e`$      |
|**10**     |0.285      |0.510      |0.639      |2.983      |2.135      |1.627      |0.234      |0.412      |0.531      |$`f`$      |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
|**20**     |0.423      |0.686      |0.787      |2.472      |1.435      |1.015      |0.342      |0.598      |0.703      |$`a_1`$    |
|**20**     |0.416      |0.678      |0.780      |2.498      |1.471      |1.046      |0.332      |0.583      |0.690      |$`a_2`$    |
|**20**     |0.405      |0.669      |0.773      |2.539      |1.503      |1.074      |0.323      |0.573      |0.683      |$`b`$      |
|**20**     |0.385      |0.657      |0.766      |2.617      |1.555      |1.101      |0.325      |0.579      |0.692      |$`c_1`$    |
|**20**     |0.380      |0.649      |0.759      |2.635      |1.588      |1.134      |0.316      |0.562      |0.678      |$`c_2`$    |
|**20**     |0.369      |0.638      |0.751      |2.675      |1.630      |1.167      |0.302      |0.547      |0.667      |$`d_1`$    |
|**20**     |0.365      |0.631      |0.745      |2.688      |1.658      |1.193      |0.298      |0.533      |0.654      |$`d_2`$    |
|**20**     |0.345      |0.612      |0.729      |2.761      |1.737      |1.261      |0.276      |0.513      |0.631      |$`e`$      |
|**20**     |0.255      |0.503      |0.637      |3.090      |2.162      |1.634      |0.213      |0.405      |0.531      |$`f`$      |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
|**100**    |0.410      |0.685      |0.787      |2.519      |1.441      |1.016      |0.334      |0.599      |0.704      |$`a_1`$    |
|**100**    |0.403      |0.676      |0.779      |2.547      |1.477      |1.048      |0.321      |0.582      |0.691      |$`a_2`$    |
|**100**    |0.392      |0.668      |0.773      |2.585      |1.511      |1.076      |0.317      |0.574      |0.683      |$`b`$      |
|**100**    |0.375      |0.656      |0.766      |2.656      |1.560      |1.103      |0.315      |0.579      |0.692      |$`c_1`$    |
|**100**    |0.368      |0.647      |0.758      |2.680      |1.594      |1.136      |0.304      |0.561      |0.676      |$`c_2`$    |
|**100**    |0.357      |0.636      |0.751      |2.719      |1.639      |1.169      |0.294      |0.547      |0.664      |$`d_1`$    |
|**100**    |0.352      |0.630      |0.744      |2.736      |1.665      |1.196      |0.285      |0.534      |0.651      |$`d_2`$    |
|**100**    |0.334      |0.610      |0.728      |2.804      |1.745      |1.263      |0.270      |0.510      |0.632      |$`e`$      |
|**100**    |0.247      |0.502      |0.637      |3.122      |2.169      |1.637      |0.208      |0.406      |0.532      |$`f`$      |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
|$`\infty`$ |0.409      |0.685      |0.787      |2.522      |1.441      |1.015      |0.334      |0.598      |0.704      |$`a_1`$    |
|$`\infty`$ |0.402      |0.676      |0.779      |2.549      |1.477      |1.047      |0.322      |0.582      |0.693      |$`a_2`$    |
|$`\infty`$ |0.392      |0.668      |0.773      |2.587      |1.510      |1.075      |0.315      |0.574      |0.681      |$`b`$      |
|$`\infty`$ |0.375      |0.655      |0.767      |2.655      |1.562      |1.099      |0.317      |0.577      |0.692      |$`c_1`$    |
|$`\infty`$ |0.369      |0.647      |0.759      |2.678      |1.595      |1.132      |0.305      |0.560      |0.676      |$`c_2`$    |
|$`\infty`$ |0.357      |0.636      |0.751      |2.719      |1.638      |1.166      |0.296      |0.548      |0.665      |$`d_1`$    |
|$`\infty`$ |0.352      |0.630      |0.745      |2.736      |1.665      |1.193      |0.286      |0.534      |0.652      |$`d_2`$    |
|$`\infty`$ |0.334      |0.610      |0.729      |2.805      |1.742      |1.260      |0.272      |0.512      |0.632      |$`e`$      |
|$`\infty`$ |0.247      |0.502      |0.637      |3.120      |2.166      |1.636      |0.208      |0.407      |0.532      |$`f`$      |


### Schedules in Two Stage Systems

Tables below give example schedules of the final round in the $`2RR + (1RR/1RR)`$ league format
employed in the championship and the relegation group in the case of 12 and 16 teams, respectively.
The integer codes represent a team's rank after the first stage of the tournament.
The schedule given in the second table was originally applied in the Polish league in the 2013/14 season.
The schedule given in the first table is its modification for 12 teams.
Notably, according to these schedules, the top half teams after the initial stage of the
competition play one more match at their home field than the bottom half. Moreover, the match
between the first and the second team after the first phase is played at the first team's home field.


| Round     |  Matches  |           |           |
|:---------:|:---------:|:---------:|:---------:|
|23         |6 - 1      |2 - 5      |3 - 4      |
|24         |1 - 3      |6 - 2      |5 - 4      |
|25         |1 - 5      |2 - 4      |3 - 6      |
|26         |4 - 1      |2 - 3      |5 - 6      |
|27         |1 - 2      |3 - 5      |4 - 6      |


| Round     |  Matches  |           |           |           |
|:---------:|:---------:|:---------:|:---------:|:---------:|
|31         |1 - 6      |2 - 5      |3 - 8      |4 - 7      |
|32         |8 - 1      |6 - 2      |7 - 3      |5 - 4      |
|33         |1 - 5      |2 - 7      |3 - 6      |4 - 8      |
|34         |7 - 1      |8 - 2      |3 - 5      |4 - 6      |
|35         |1 - 3      |2 - 4      |5 - 7      |6 - 8      |
|36         |4 - 1      |2 - 3      |6 - 7      |8 - 5      |
|37         |1 - 2      |3 - 4      |5 - 6      |7 - 8      |

