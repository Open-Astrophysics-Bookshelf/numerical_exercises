EPStoPDF = epstopdf

ALL: incompressible.ps incompressible.pdf


%.dvi: %.tex refs.bib
	latex $< < /dev/null
	bibtex incompressible
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
