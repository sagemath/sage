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
	rm -rf "$(SAGE_ROOT)/build/sage-distro/src/doc"
	# and a few more things that were previously generated
	rm -rf en/reference/*/sage
	rm -rf en/reference/documentation/sage_docbuild
	rm -rf en/reference/sage
	rm -f common/*.pyc

doc-html:
	meson compile -C ../../build/sage-distro doc-html

doc-pdf:
	meson compile -C ../../build/sage-distro doc-pdf

.PHONY: all clean \
	doc-src \
	doc-html doc-pdf \
	doc-inventory-reference doc-html-reference doc-pdf-reference \
	doc-html-other doc-pdf-other
