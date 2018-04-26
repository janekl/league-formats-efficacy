i=1
n_sim=100000
n_cores=24
folder=results_all
model=poisson_correlated
drift=fixed
teams=12

for shape in -1 5 7.5 10 20 50 100 200 500
do
  for sigma in 0.1 0.15 0.2 0.25 0.30 0.35 0.4
  do
    echo $i: "("$drift, $teams, $shape, $sigma")"
    nice -20 Rscript run_simulations.R --n=$teams --model=$model --n_sim=$n_sim --shape=$shape --sigma=$sigma --drift_option=$drift --n_cores=$n_cores --specific_result_folder=$folder --log2file=1
    i=$((i+1))
  done
done
