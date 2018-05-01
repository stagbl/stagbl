PDFS = slides.pdf

PDFLATEX=pdflatex --halt-on-error
BIBTEX=bibtex
PYTHON3=python3

PLOT_PDFS= \
			 images/A11_breakdown_plot1.pdf \
			 images/A11_breakdown_plot2.pdf \

all : $(PDFS) $(PLOT_PDFS)

$(PDFS) : $(PLOT_PDFS)

%.pdf : %.tex
	$(PDFLATEX) $<
	$(BIBTEX) $(<:%.tex=%)
	$(PDFLATEX) $<
	$(PDFLATEX) $<

images/A11_breakdown_plot%.pdf : images/A11_breakdown_plot.py
	$(PYTHON3) $<
	mv A11_breakdown_plot*.pdf images

# Keep the intermediate images
.SECONDARY : $(PLOT_PDFS)

# Remove intermediate files
clean :
	rm -f *.aux *.log *.bbl *.blg *-blx.bib *.nav *.snm *.toc *.vrb *.run.xml *.out *.spl *.fls *.fdb_latexmk
	rm -f $(PDFS)

.PHONY : all clean
