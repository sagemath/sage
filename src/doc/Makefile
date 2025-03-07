########################################################
#
# 'make doc-html' (synonym for 'make' and 'make all') builds the html documentation.
# 'make doc-pdf' builds the PDF documentation.
#
# The "sage" script must be in PATH, and the Sage library and the package sage_docbuild
# must be installed for these to work.
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

# Sources generated at build time. (For sources generated at bootstrap time, see bootstrap.)
doc-src:
	mkdir -p en/reference/repl
	sage -advanced > en/reference/repl/options.txt

# Matches doc-inventory--reference-manifolds etc.
doc-inventory--%:
	sage --docbuild --no-pdf-links $(subst -,/,$(subst doc-inventory--,,$@)) inventory $(SAGE_DOCBUILD_OPTS)

# Matches doc-html--developer, doc-html--reference-manifolds etc.
doc-html--%:
	sage --docbuild --no-pdf-links $(subst -,/,$(subst doc-html--,,$@)) html $(SAGE_DOCBUILD_OPTS)

# reference manual, inventory
ifndef SAGE_ROOT
doc-inventory-reference:
	$(error SAGE_ROOT undefined. This Makefile needs to be invoked by build/make/install)
else
doc-inventory-reference: doc-src
	$(eval DOCS = $(shell sage --docbuild --all-documents reference))
	@if [ -z "$(DOCS)" ]; then echo "Error: 'sage --docbuild --all-documents' failed"; exit 1; fi
	$(eval BIBLIO = $(firstword $(DOCS)))
	$(eval OTHER_DOCS = $(wordlist 2, 100, $(DOCS)))
	$(MAKE) doc-inventory--$(subst /,-,$(BIBLIO))
	$(MAKE) $(foreach doc, $(OTHER_DOCS), doc-inventory--$(subst /,-,$(doc)))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" doc-inventory--reference_top
endif

# sub docs of reference manual, html
doc-html-reference-sub: doc-inventory-reference
	$(eval DOCS = $(shell sage --docbuild --all-documents reference))
	@if [ -z "$(DOCS)" ]; then echo "Error: 'sage --docbuild --all-documents' failed"; exit 1; fi
	$(eval BIBLIO = $(firstword $(DOCS)))
	$(eval OTHER_DOCS = $(wordlist 2, 100, $(DOCS)))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" doc-html--$(subst /,-,$(BIBLIO))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" $(foreach doc, $(OTHER_DOCS), doc-html--$(subst /,-,$(doc)))

# reference manual, html; reference_top is built after sub docs
doc-html-reference: doc-html-reference-sub
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" doc-html--reference_top

# other documentation, html
doc-html-other: doc-html-reference
	$(eval DOCS = $(shell sage --docbuild --all-documents all))
	@if [ -z "$(DOCS)" ]; then echo "Error: 'sage --docbuild --all-documents' failed"; exit 1; fi
	$(MAKE) $(foreach doc, $(wordlist 2, 100, $(DOCS)), doc-html--$(subst /,-,$(doc)))

doc-html: doc-html-reference doc-html-other
	SAGE_DOC=$$(sage --python -c "from sage.env import SAGE_DOC; print(SAGE_DOC)")
	find $${SAGE_DOC}/html -type d -path "*/jupyter_execute" -exec rm -rf {} +

# Matches doc-pdf--developer, doc-pdf--reference-manifolds etc.
doc-pdf--%:
	LATEXOPTS="--file-line-error --interaction=batchmode" sage --docbuild $(subst -,/,$(subst doc-pdf--,,$@)) pdf $(SAGE_DOCBUILD_OPTS)

# reference manual, pdf
doc-pdf-reference: doc-inventory-reference
	$(eval DOCS = $(shell sage --docbuild --all-documents reference))
	@if [ -z "$(DOCS)" ]; then echo "Error: 'sage --docbuild --all-documents' failed"; exit 1; fi
	$(eval BIBLIO = $(firstword $(DOCS)))
	$(eval OTHER_DOCS = $(wordlist 2, 100, $(DOCS)))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" doc-pdf--$(subst /,-,$(BIBLIO))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" $(foreach doc, $(OTHER_DOCS), doc-pdf--$(subst /,-,$(doc)))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" doc-pdf--reference_top

# other documentation, pdf
doc-pdf-other: doc-pdf-reference
	$(eval DOCS = $(shell sage --docbuild --all-documents all))
	@if [ -z "$(DOCS)" ]; then echo "Error: 'sage --docbuild --all-documents' failed"; exit 1; fi
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" $(foreach doc, $(wordlist 2, 100, $(DOCS)), doc-pdf--$(subst /,-,$(doc)))

doc-pdf: doc-pdf-reference doc-pdf-other
	SAGE_DOC=$$(sage --python -c "from sage.env import SAGE_DOC; print(SAGE_DOC)")
	find $${SAGE_DOC}/latex -type d -path "*/jupyter_execute" -exec rm -rf {} +

.PHONY: all clean \
	doc-src \
	doc-html doc-pdf \
	doc-inventory-reference doc-html-reference doc-pdf-reference \
	doc-html-other doc-pdf-other
