EPStoPDF = epstopdf

ALL: finite-volume.ps finite-volume.pdf


%.dvi: %.tex refs.bib
	latex $< < /dev/null
	bibtex finite-volume
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
