EPStoPDF = epstopdf

ALL: incompressible.pdf


incompressible.pdf: main.tex incompressible.tex
	pdflatex -jobname=incompressible main.tex < /dev/null
	bibtex incompressible.aux
	pdflatex -jobname=incompressible main.tex < /dev/null
	pdflatex -jobname=incompressible main.tex < /dev/null


clean:
	$(RM) *.aux *.log *.dvi *.bbl *.blg 
	$(RM) *~

.PHONY: clean
