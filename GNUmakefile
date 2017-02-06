# by default, typing make will build the `production' version of the 
# Computational Hydrodynamics notes.
#
# if you instead for 'make DEBUG=t', then it will make the draft version.
# this will enable margin comments and build in any chapters that are
# not yet complete.

ALL: CompHydroTutorial.pdf

DEBUG := 

DIRS := preface \
        intro \
	pde-classes \
        finite-volume \
        advection \
        burgers \
        Euler \
        hydro-test-problems \
        multigrid \
        diffusion \
        multiphysics \
        incompressible \
        low_mach \
        pyro \
        hydro_examples \
        software-engineering \
        symbols \
        radiation

# dependencies
TEXS := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.tex))
EPSS := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.eps))

# for PDFs, pdflatex will automagically convert file.eps to
# file-converted-to.pdf at build time.  These lines create
# a variable, PDFS, that just contains the original PDF files
# as dependencies
tPDFS := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.pdf))
cPDFS := $(foreach dir, $(DIRS), $(wildcard $(dir)/*converted-to.pdf))
PDFS := $(filter-out $(cPDFS), $(tPDFS)) 


conditional.tex:
ifdef DEBUG
	echo "\def\debugmode{}" > conditional.tex
else
	echo "" > conditional.tex
endif

CompHydroTutorial.pdf: CompHydroTutorial.tex conditional.tex $(TEXS) $(EPSS) $(PDFS) refs.bib
	git rev-parse --short=12 HEAD > git_info.tex
	pdflatex CompHydroTutorial  < /dev/null
	bibtex CompHydroTutorial.aux
	pdflatex CompHydroTutorial  < /dev/null
	pdflatex CompHydroTutorial  < /dev/null
	pdflatex CompHydroTutorial  < /dev/null
	$(RM) git_info.tex
	$(RM) conditional.tex


clean:
	$(RM) *.aux *.log *.dvi *.bbl *.blg *.lof *.toc *.exc *.out *.brf
	$(RM) *~ conditional.tex

realclean: clean
	$(RM) CompHydroTutorial.pdf
	find . -name "*-converted-to.pdf" -exec $(RM) {} \;

.PHONY: clean


print-%: ; @echo $* is $($*)
