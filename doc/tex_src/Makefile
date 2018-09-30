MAIN = snowglobes_1.3-draft

TEXFILES = $(MAIN).tex $(wildcard sec_*.tex) $(wildcard *.bib)
FIGURES = $(shell ls figures/*)

PDFFILE = $(MAIN).pdf

LATEX = pdflatex
BIBTEX = bibtex

all : $(PDFFILE)

$(PDFFILE) : $(TEXFILES)
	$(LATEX) $(MAIN).tex
	$(LATEX) $(MAIN).tex
	$(BIBTEX) $(MAIN)
	$(LATEX) $(MAIN).tex
	$(LATEX) $(MAIN).tex
	$(LATEX) $(MAIN).tex

clean :
	$(RM) *~ $(PDFFILE) *.aux *.log *.toc */*.aux *.lot *.lof *.blg *.bbl *.out
clena :
	$(RM) *~ $(PDFFILE) *.aux *.log *.toc */*.aux *.lot *.lof *.blg *.bbl *.out
