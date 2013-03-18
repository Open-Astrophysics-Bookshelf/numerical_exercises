EPStoPDF = epstopdf

ALL: Euler.ps Euler.pdf


%.dvi: %.tex refs.bib
	latex $< < /dev/null
	bibtex Euler
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
