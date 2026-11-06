########################################################
#
# 'make doc-html' (synonym for 'make' and 'make all') builds the html documentation.
# 'make doc-pdf' builds the PDF documentation.
#
# The "sage" script must be in PATH, and the Sage library must be installed for these to work.
# These commands are intended to be called by build/make/Makefile via
# build/pkgs/sagemath_doc_{html,pdf}/spkg-install.
#
# 'make clean' removes build artifacts; the "sage" script is not needed for this to work.

# Documentation build parallelism:
# PDF builds are memory-intensive, so limit them to 2 jobs by default.
SAGE_DOC_JOBS ?=
SAGE_DOC_HTML_JOBS ?= $(SAGE_DOC_JOBS)
SAGE_DOC_PDF_JOBS ?= $(if $(SAGE_DOC_JOBS),$(SAGE_DOC_JOBS),2)

all: doc-html

clean:
	rm -rf "$(SAGE_ROOT)/build/sage-distro/src/doc"
	# and a few more things that were previously generated
	rm -rf en/reference/*/sage
	rm -rf en/reference/documentation/sage_docbuild
	rm -rf en/reference/sage
	rm -f common/*.pyc

doc-html:
	$(if $(SAGE_DOC_NUM_THREADS),SAGE_NUM_THREADS=$(SAGE_DOC_NUM_THREADS) )meson compile -C ../../build/sage-distro $(if $(SAGE_DOC_HTML_JOBS),-j $(SAGE_DOC_HTML_JOBS) )doc-html

doc-pdf:
	SAGE_NUM_THREADS=$(SAGE_DOC_PDF_THREADS) meson compile -C ../../build/sage-distro $(if $(SAGE_DOC_PDF_JOBS),-j $(SAGE_DOC_PDF_JOBS) )doc-pdf

.PHONY: all clean \
	doc-src \
	doc-html doc-pdf \
	doc-inventory-reference doc-html-reference doc-pdf-reference \
	doc-html-other doc-pdf-other
