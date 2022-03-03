########################################################
#
# 'make doc-html' (synonym for 'make' and 'make all') builds the html documentation.
# 'make doc-pdf' builds the PDF documentation.
#
# SAGE_ROOT must be defined for these to work, so these commands
# should be called by build/make/Makefile.
#
# 'make clean' removes build artifacts; SAGE_ROOT need not be defined for this to work.

all: doc-html

clean:
	rm -rf en/reference/*/sage
	rm -rf en/reference/sage
	rm -f common/*.pyc

# Matches doc-inventory--reference-manifolds etc.
doc-inventory--%:
	cd $(SAGE_ROOT) && ./sage --docbuild --no-pdf-links $(subst -,/,$(subst doc-inventory--,,$@)) inventory $(SAGE_DOCBUILD_OPTS)

# Matches doc-html--developer, doc-html--reference-manifolds etc.
doc-html--%:
	cd $(SAGE_ROOT) && ./sage --docbuild --no-pdf-links $(subst -,/,$(subst doc-html--,,$@)) html $(SAGE_DOCBUILD_OPTS)

# reference manual, inventory
ifndef SAGE_ROOT
doc-inventory-reference:
	$(error SAGE_ROOT undefined. This Makefile needs to be invoked by build/make/install)
else
doc-inventory-reference:
	$(eval DOCS = $(shell cd $(SAGE_ROOT) && ./sage --docbuild --all-documents reference))
	@if [ -z "$(DOCS)" ]; then echo "Error: './sage --docbuild --all-documents' failed"; exit 1; fi
	$(eval BIBLIO = $(firstword $(DOCS)))
	$(eval OTHER_DOCS = $(wordlist 2, 100, $(DOCS)))
	$(MAKE) doc-inventory--$(subst /,-,$(BIBLIO))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" $(foreach doc, $(OTHER_DOCS), doc-inventory--$(subst /,-,$(doc)))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" doc-inventory--reference_top
endif

# reference manual, html
doc-html-reference: doc-inventory-reference
	$(eval DOCS = $(shell cd $(SAGE_ROOT) && ./sage --docbuild --all-documents reference))
	@if [ -z "$(DOCS)" ]; then echo "Error: './sage --docbuild --all-documents' failed"; exit 1; fi
	$(eval BIBLIO = $(firstword $(DOCS)))
	$(eval OTHER_DOCS = $(wordlist 2, 100, $(DOCS)))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" doc-html--$(subst /,-,$(BIBLIO))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" $(foreach doc, $(OTHER_DOCS), doc-html--$(subst /,-,$(doc)))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" doc-html--reference_top

# other documentation, html
doc-html-other: doc-html-reference
	$(eval DOCS = $(shell cd $(SAGE_ROOT) && ./sage --docbuild --all-documents all))
	@if [ -z "$(DOCS)" ]; then echo "Error: './sage --docbuild --all-documents' failed"; exit 1; fi
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" $(foreach doc, $(wordlist 2, 100, $(DOCS)), doc-html--$(subst /,-,$(doc)))

doc-html: doc-html-reference doc-html-other

# Matches doc-pdf--developer, doc-pdf--reference-manifolds etc.
doc-pdf--%:
	cd $(SAGE_ROOT) && ./sage --docbuild $(subst -,/,$(subst doc-pdf--,,$@)) pdf $(SAGE_DOCBUILD_OPTS)

# reference manual, pdf
doc-pdf-reference: doc-inventory-reference
	$(eval DOCS = $(shell cd $(SAGE_ROOT) && ./sage --docbuild --all-documents reference))
	@if [ -z "$(DOCS)" ]; then echo "Error: './sage --docbuild --all-documents' failed"; exit 1; fi
	$(eval BIBLIO = $(firstword $(DOCS)))
	$(eval OTHER_DOCS = $(wordlist 2, 100, $(DOCS)))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" doc-pdf--$(subst /,-,$(BIBLIO))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" $(foreach doc, $(OTHER_DOCS), doc-pdf--$(subst /,-,$(doc)))
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" doc-pdf--reference_top

# other documentation, pdf
doc-pdf-other: doc-html-reference
	$(eval DOCS = $(shell cd $(SAGE_ROOT) && ./sage --docbuild --all-documents all))
	@if [ -z "$(DOCS)" ]; then echo "Error: './sage --docbuild --all-documents' failed"; exit 1; fi
	$(MAKE) SAGE_DOCBUILD_OPTS="$(SAGE_DOCBUILD_OPTS) --no-prune-empty-dirs" $(foreach doc, $(wordlist 2, 100, $(DOCS)), doc-pdf--$(subst /,-,$(doc)))

# website with pdf links
doc-pdf-website:
	cd $(SAGE_ROOT) && ./sage --docbuild website html $(SAGE_DOCBUILD_OPTS)

doc-pdf: doc-pdf-reference doc-pdf-other doc-pdf-website

.PHONY: all clean \
	doc-html doc-pdf \
	doc-inventory-reference doc-html-reference doc-pdf-reference \
	doc-html-other doc-pdf-other
