EPStoPDF = epstopdf

ALL: advection.pdf


advection.pdf: main.tex advection.tex
	pdflatex -jobname=advection main.tex < /dev/null
	bibtex advection.aux
	pdflatex -jobname=advection main.tex < /dev/null
	pdflatex -jobname=advection main.tex < /dev/null


clean:
	$(RM) *.aux *.log *.dvi *.bbl *.blg 
	$(RM) *~

.PHONY: clean
