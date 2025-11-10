.PHONY: pull push

MSG ?= Auto-commit from makefile

pull:
	git pull
	R CMD INSTALL .

push:
	Rscript -e "roxygen2::roxygenize()"
	R CMD INSTALL .
	git add .
	git commit -m "$(MSG)"
	git push
