EPStoPDF = epstopdf

ALL: multigrid.ps multigrid.pdf


%.dvi: %.tex refs.bib
	latex $< < /dev/null
	bibtex multigrid
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
