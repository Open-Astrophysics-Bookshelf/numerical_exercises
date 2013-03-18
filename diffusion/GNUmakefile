EPStoPDF = epstopdf

ALL: diffusion.ps diffusion.pdf


%.dvi: %.tex refs.bib
	latex $< < /dev/null
	bibtex diffusion
	latex $< < /dev/null
	latex $< < /dev/null

%.ps: %.dvi
	dvips -t letter -Ppdf -o $@ $<

%.pdf: %.ps
	ps2pdf $<

clean:
	$(RM) *.aux *.log *.dvi
	$(RM) *~

.PHONY: clean
