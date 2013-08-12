EPStoPDF = epstopdf

ALL: finite-volume.pdf


finite-volume.pdf: main.tex finite-volume.tex
	pdflatex -jobname=finite-volume main.tex < /dev/null
	bibtex finite-volume.aux
	pdflatex -jobname=finite-volume main.tex < /dev/null
	pdflatex -jobname=finite-volume main.tex < /dev/null


clean:
	$(RM) *.aux *.log *.dvi *.bbl *.blg 
	$(RM) *~

.PHONY: clean
