EPStoPDF = epstopdf

ALL: Euler.pdf


Euler.pdf: main.tex Euler.tex
	pdflatex -jobname=Euler main.tex < /dev/null
	bibtex Euler.aux
	pdflatex -jobname=Euler main.tex < /dev/null
	pdflatex -jobname=Euler main.tex < /dev/null


clean:
	$(RM) *.aux *.log *.dvi *.bbl *.blg 
	$(RM) *~

.PHONY: clean
