
SETTING = cond

foster-all:
	Rscript --vanilla inst/scripts/break-foster.R cond
	Rscript --vanilla inst/scripts/break-foster.R cond-noX
	Rscript --vanilla inst/scripts/break-foster.R marg
	Rscript --vanilla inst/scripts/break-foster.R marg-noX

ex1:
	Rscript --vanilla inst/scripts/example1-continuous.R 0 0 0
	Rscript --vanilla inst/scripts/example1-continuous.R 1 0 0
	Rscript --vanilla inst/scripts/example1-continuous.R 0 1 0
	Rscript --vanilla inst/scripts/example1-continuous.R 0 0 1
	Rscript --vanilla inst/scripts/example1-continuous.R 1 1 0
	Rscript --vanilla inst/scripts/example1-continuous.R 1 0 1
	Rscript --vanilla inst/scripts/example1-continuous.R 0 1 1
	Rscript --vanilla inst/scripts/example1-continuous.R 1 1 1
