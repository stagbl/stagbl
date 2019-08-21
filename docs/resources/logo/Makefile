PDFS = logo.pdf logo_slides.pdf
SVGS = logo.svg

PDFLATEX=pdflatex --halt-on-error

all : $(PDFS) $(SVGS)

%.pdf : %.tex
	$(PDFLATEX) $<

# requires ImageMagick
%.svg : %.pdf
	convert $^ $@

# Remove intermediate files
clean :
	rm -f *.aux *.log *.bbl *.blg *-blx.bib *.nav *.snm *.toc *.vrb *.run.xml *.out *.spl *.fls *.fdb_latexmk

.PHONY : all clean
