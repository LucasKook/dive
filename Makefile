401k:
	Rscript --vanilla inst/code/401k.R
	Rscript --vanilla inst/code/401k-comparison.R
	Rscript --vanilla inst/code/ivqr-401k.R

schooling:
	Rscript --vanilla inst/code/schooling.R
	Rscript --vanilla inst/code/schooling-comparison.R
	Rscript --vanilla inst/code/ivqr-schooling.R

run-simulations:
	Rscript --vanilla inst/code/simulations.R

vis-simulations:
	Rscript --vanilla inst/code/vis-sim-results.R

loss-landscape:
	Rscript --vanilla inst/code/loss-landscape.R

figures:
	Rscript --vanilla inst/code/plot-examples.R 1
	Rscript --vanilla inst/code/plot-examples.R 2
	Rscript --vanilla inst/code/plot-examples.R 3
	Rscript --vanilla inst/code/plot-examples.R 4

all: 401k schooling run-simulations vis-simulations figures loss-landscape
