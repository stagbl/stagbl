PDFS = logo.pdf

PDFLATEX=pdflatex --halt-on-error

all : $(PDFS)

%.pdf : %.tex
	$(PDFLATEX) $<

# Remove intermediate files
clean :
	rm -f *.aux *.log *.bbl *.blg *-blx.bib *.nav *.snm *.toc *.vrb *.run.xml *.out *.spl *.fls *.fdb_latexmk

.PHONY : all clean
