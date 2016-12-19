
DOCS = sd1d-manual.pdf

.PHONY: all

all: $(DOCS)

%.pdf: %.tex 
	pdflatex $< 
	pdflatex $<

