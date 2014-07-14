EPStoPDF = epstopdf

ALL: CompHydroTutorial.pdf

DIRS := preface \
        intro \
        finite-volume \
        advection \
        burgers \
        Euler \
        multigrid \
        diffusion \
        multiphysics \
        incompressible \
        low_mach \
        pyro \
        hydro_examples \
        software-engineering

TEXS := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.tex))
EPSS := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.eps))

#TEXS += git_info.tex

#git_info.tex:


CompHydroTutorial.pdf: CompHydroTutorial.tex $(TEXS) $(EPSS) refs.bib
	git rev-parse HEAD > git_info.tex
	pdflatex CompHydroTutorial  < /dev/null
	bibtex CompHydroTutorial.aux
	pdflatex CompHydroTutorial  < /dev/null
	pdflatex CompHydroTutorial  < /dev/null
	$(RM) git_info.tex


clean:
	$(RM) *.aux *.log *.dvi *.bbl *.blg *.lof *.toc *.exc *.out
	$(RM) *~

realclean: clean
	$(RM) CompHydroTutorial.pdf

.PHONY: clean
