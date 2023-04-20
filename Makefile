
SETTING = cond

foster-all:
	Rscript --vanilla inst/scripts/break-foster.R cond
	Rscript --vanilla inst/scripts/break-foster.R cond-noX
	Rscript --vanilla inst/scripts/break-foster.R marg
	Rscript --vanilla inst/scripts/break-foster.R marg-noX
