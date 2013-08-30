EPStoPDF = epstopdf

ALL: CompHydroTutorial.pdf

TEXS := intro/intro.tex \
        finite-volume/finite-volume.tex \
        advection/advection.tex \
        Euler/Euler.tex \
        multigrid/multigrid.tex \
        diffusion/diffusion.tex \
        multiphysics/multiphysics.tex \
        incompressible/incompressible.tex 

CompHydroTutorial.pdf: CompHydroTutorial.tex $(TEXS)
	pdflatex CompHydroTutorial  < /dev/null
	bibtex CompHydroTutorial.aux
	pdflatex CompHydroTutorial  < /dev/null
	pdflatex CompHydroTutorial  < /dev/null


clean:
	$(RM) *.aux *.log *.dvi *.bbl *.blg 
	$(RM) *~

.PHONY: clean
