# Build PDF versions of the manuscript and supplement.
#
#   make            build pdf/index.pdf and pdf/supplemental.pdf
#   make index      main text only
#   make supplement supplement only
#   make clean      remove built PDFs
#
# Requires pandoc and a XeLaTeX install (TinyTeX works).
# Layout tweaks: pdf/pdf-preamble.tex (ragged right) and
# pdf/pagebreaks.lua (which sections start on a new page).

PANDOC   ?= pandoc
OUTDIR   := pdf
FILTER   := $(OUTDIR)/pagebreaks.lua
PREAMBLE := $(OUTDIR)/pdf-preamble.tex

PDF_OPTS := --pdf-engine=xelatex \
            -V title= \
            -V geometry:margin=1in \
            -V mainfont="Palatino" \
            -V mathfont="STIX Two Math" \
            -V linestretch=1.15 \
            -V colorlinks=true \
            --lua-filter=$(FILTER) \
            -H $(PREAMBLE)

MEDIA := $(wildcard media/*)

.PHONY: all index supplement clean

all: index supplement

index: $(OUTDIR)/index.pdf
supplement: $(OUTDIR)/supplemental.pdf

$(OUTDIR)/%.pdf: %.md %.meta.yml $(FILTER) $(PREAMBLE) $(MEDIA)
	$(PANDOC) $< --metadata-file=$*.meta.yml $(PDF_OPTS) -o $@

clean:
	rm -f $(OUTDIR)/index.pdf $(OUTDIR)/supplemental.pdf
