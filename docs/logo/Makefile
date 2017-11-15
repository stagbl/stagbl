PDFS = logo.pdf

PDFLATEX=pdflatex --halt-on-error

all : $(PDFS) logo.svg

%.pdf : %.tex
	$(PDFLATEX) $<

# requires ImageMagick
logo.svg : logo.pdf
	convert logo.pdf logo.svg 

# Remove intermediate files
clean :
	rm -f *.aux *.log *.bbl *.blg *-blx.bib *.nav *.snm *.toc *.vrb *.run.xml *.out *.spl *.fls *.fdb_latexmk

.PHONY : all clean
