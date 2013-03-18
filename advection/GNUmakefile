EPStoPDF = epstopdf

ALL: advection.ps advection.pdf


%.dvi: %.tex
	latex $< < /dev/null
	bibtex advection       
	latex $< < /dev/null
	latex $< < /dev/null

%.ps: %.dvi
	dvips -t letter -Ppdf -o $@ $<

%.pdf: %.ps
	ps2pdf $<

clean:
	$(RM) *.aux *.log *.dvi *.bbl *.blg 
	$(RM) *~

.PHONY: clean
