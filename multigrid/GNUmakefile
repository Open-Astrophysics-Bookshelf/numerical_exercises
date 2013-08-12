EPStoPDF = epstopdf

ALL: multigrid.pdf


multigrid.pdf: main.tex multigrid.tex
	pdflatex -jobname=multigrid main.tex < /dev/null
	bibtex multigrid.aux
	pdflatex -jobname=multigrid main.tex < /dev/null
	pdflatex -jobname=multigrid main.tex < /dev/null


clean:
	$(RM) *.aux *.log *.dvi *.bbl *.blg 
	$(RM) *~

.PHONY: clean
