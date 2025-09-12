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

all: doc-html

clean:
	rm -rf en/reference/*/sage
	rm -rf en/reference/documentation/sage_docbuild
	rm -rf en/reference/sage
	rm -f common/*.pyc
	# clean the meson build directory
	rm -rf ../../build/sage-distro/src/doc

doc-html:
	meson compile -C ../../build/sage-distro doc-html
	mkdir -p $(SAGE_DOC)
	cp -rf ../../build/sage-distro/src/doc/html $(SAGE_DOC)
	cp -f ../../build/sage-distro/src/doc/index.html $(SAGE_DOC)

doc-pdf:
	meson compile -C ../../build/sage-distro doc-pdf
	mkdir -p $(SAGE_DOC)
	cp -rf ../../build/sage-distro/src/doc/pdf $(SAGE_DOC)

.PHONY: all clean \
	doc-src \
	doc-html doc-pdf \
	doc-inventory-reference doc-html-reference doc-pdf-reference \
	doc-html-other doc-pdf-other
