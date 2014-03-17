EPStoPDF = epstopdf

ALL: CompHydroTutorial.pdf

DIRS := intro \
        finite-volume \
        advection \
        burgers \
        Euler \
        multigrid \
        diffusion \
        multiphysics \
        incompressible

TEXS := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.tex))
EPSS := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.eps))


CompHydroTutorial.pdf: CompHydroTutorial.tex $(TEXS) $(EPSS) refs.bib
	pdflatex CompHydroTutorial  < /dev/null
	bibtex CompHydroTutorial.aux
	pdflatex CompHydroTutorial  < /dev/null
	pdflatex CompHydroTutorial  < /dev/null


clean:
	$(RM) *.aux *.log *.dvi *.bbl *.blg 
	$(RM) *~

.PHONY: clean
