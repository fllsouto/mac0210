#!/bin/bash
# LaTeX Makefile

FILE=main

all: $(FILE).pdf

.PHONY: clean
clean:
	rm -fv *.aux *.out *.log *.bbl *.blg *.brf *.idx *.toc
	rm -fv chapters/*.aux
	rm -fv bib/*.bbl bib/*.blg bib/*.brf


clean_all:
	rm -fv *.aux *.out *.log *.bbl *.blg *.brf *.idx *.toc *.pdf *.lof *.lot
	rm -fv partes/*.aux
	rm -fv bib/*.bbl bib/*.blg bib/*.brf bib/*.aux
	rm -fv front/*.bbl front/*.blg front/*.brf front/*.aux

$(FILE).pdf: $(FILE).tex
	pdflatex $(FILE)
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)