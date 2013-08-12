EPStoPDF = epstopdf

ALL: diffusion.pdf


diffusion.pdf: main.tex diffusion.tex
	pdflatex -jobname=diffusion main.tex < /dev/null
	bibtex diffusion.aux
	pdflatex -jobname=diffusion main.tex < /dev/null
	pdflatex -jobname=diffusion main.tex < /dev/null


clean:
	$(RM) *.aux *.log *.dvi *.bbl *.blg 
	$(RM) *~

.PHONY: clean
