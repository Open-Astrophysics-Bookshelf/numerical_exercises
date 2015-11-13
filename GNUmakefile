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
        software-engineering \
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


CompHydroTutorial.pdf: CompHydroTutorial.tex $(TEXS) $(EPSS) $(PDFS) refs.bib
	git rev-parse --short=12 HEAD > git_info.tex
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
	find . -name "*-converted-to.pdf" -exec $(RM) {} \;

.PHONY: clean


print-%: ; @echo $* is $($*)
