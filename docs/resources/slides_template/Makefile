PDFS = slides.pdf
BIBFILES = pds.bib

PDFLATEX=pdflatex --halt-on-error
BIBTEX=bibtex

all : $(PDFS)

# Main compile with latexmk, if available, otherwise with brute force
ifndef LATEXMK
LATEXMK=
endif
ifeq (${LATEXMK},)
%.pdf : %.tex $(BIBFILES)
	$(PDFLATEX) $<
	$(BIBTEX) $(<:.tex=)
	$(PDFLATEX) $<
	$(PDFLATEX) $<
else
%.pdf : %.tex $(BIBFILES)
	$(LATEXMK) $<
endif

# Remove intermediate files
clean :
	rm -f *.aux *.log *.bbl *.blg *-blx.bib *.nav *.snm *.toc *.vrb *.run.xml *.out *.spl *.fls *.fdb_latexmk
	rm -f $(PDFS)

.PHONY : all clean
