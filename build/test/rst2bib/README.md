# Comparing documentation before and after the BibTeX migration

These commands compare existing documentation outputs. They accept file or
directory paths and do not need an installed Sage library. Use the same Sphinx
version, theme, and extensions for the before and after builds; the recorded
migration comparisons used Sphinx 9.1.0.

In a Python environment for these tools, install the dependencies:

```sh
python -m pip install -r build/test/rst2bib/requirements.txt
```

Run the examples below from the Sage checkout. Replace the `before` and `after`
paths with your saved outputs. The commands return zero when the comparison
passes and nonzero otherwise. Without `--output`, the full report is JSON on
standard output. With `--output`, the full report is saved to that file and its
summary is printed. All three commands also provide `--help`.

## Reference catalogue

`compare-html-dom.py` compares each citation's anchor, normalized text, links and
their captions, mathematics, emphasis, bold text, and inline code. It excludes
backlink labels and normalizes whitespace. Missing citations always fail;
additional citations require `--allow-added`.

```sh
python build/test/rst2bib/compare-html-dom.py \
  before/html/en/reference/references/index.html \
  after/html/en/reference/references/index.html \
  --allow-added \
  --expected-differences build/test/rst2bib/catalogue-differences.json \
  --output catalogue-comparison.json
```

`catalogue-differences.json` records the exact field changes for `GrLe1996` and
`VT2001` described in the [migration report](../../../BIBTEX_MIGRATION_REPORT.md).
The entire differences object must match this file; it does not allow arbitrary
changes to those two entries. Omit `--expected-differences` when no changed
reference is expected.

For changed citation keys, `--key-map` accepts a JSON object mapping old keys to
new keys. `--ignore-ids` excludes anchors when the before page intentionally uses
different IDs. Both options are explicit so the catalogue's stable anchors are
checked by default.

The separate 177-entry docstring fixture used in the migration report has labels
`LocalRef000` through `LocalRef176`. For that fixture:

```sh
python build/test/rst2bib/compare-html-dom.py \
  before-local/index.html after/html/en/reference/references/index.html \
  --allow-added --ignore-ids \
  --key-map build/test/rst2bib/local-key-map.json \
  --expected-differences build/test/rst2bib/local-differences.json \
  --output local-comparison.json
```

`local-key-map.json` belongs to that aggregate fixture, not ordinary module pages.
The fixture HTML and other generated build outputs are not included here. For
other inputs, provide your own key map and expected differences as needed.

## Module prose and examples

`compare-module-pages.py` accepts either two HTML files or two documentation
directories. Directory mode compares all `.html` files unless `--page` selects
specific relative paths. Added or missing selected pages fail. Both prose and
example differences cause a nonzero exit status and are included in the report.

To compare one module, allowing the bibliography relocation and the source-path
heading change:

```sh
python build/test/rst2bib/compare-module-pages.py \
  before/module.html after/module.html \
  --ignore-title --ignore-bibliography \
  --key-map build/test/rst2bib/citation-key-aliases.json \
  --output module-comparison.json
```

The exclusions are optional. `--ignore-title` removes the first `h1`;
`--ignore-bibliography` removes citation blocks, citation-only list items, and
the bibliography rubric. Citation-key renames are applied once to the old page's
prose. Prose whitespace is normalized; highlighted example blocks are compared
exactly. This check does not compare module-page hyperlinks or visual layout;
the citation comparer handles reference links separately.

To reproduce the three-module comparison on the saved single-file build layout:

```sh
python build/test/rst2bib/compare-module-pages.py before/modules after/modules \
  --page difference_family/output/html/index.html \
         dyck_word/output/html/index.html \
         diagram_algebras/output/html/index.html \
  --ignore-title --ignore-bibliography \
  --key-map build/test/rst2bib/citation-key-aliases.json \
  --output module-pages.json
```

The expected result is three unchanged pages with 546 unchanged example blocks.

## Splitting the database into several files

`compare-multibib.py` compares parsed BibTeX records and preambles, citation HTML
in document order, and the catalogue's main HTML. `--before-bib` and `--after-bib`
accept one or more files, or directories whose `.bib` files are combined.
Duplicate keys, including case-insensitive duplicates across files, are rejected.

```sh
python build/test/rst2bib/compare-multibib.py \
  --before-bib before/reference-database.bib \
  --after-bib src/doc/references \
  --before-html before/index.html \
  --after-html after/html/en/reference/references/index.html \
  --before-latex before/references.tex \
  --after-latex after/latex/en/reference/references/references.tex \
  --before-pdf before/references.pdf \
  --after-pdf after/pdf/en/reference/references/references.pdf \
  --output multibib-comparison.json
```

The LaTeX and PDF arguments are optional; each before/after pair must be supplied
together. LaTeX source is compared byte-for-byte. PDF text is extracted with
Ghostscript, which must be available as `gs` on `PATH`, or specified with `--gs`.
The report records whole-HTML and PDF binary equality, but those flags do not
affect success: HTML outside the main content and PDF build metadata may differ.
The required checks still enforce identical main/citation HTML and, when
requested, identical LaTeX and extracted PDF text.

For the recorded two-file split, all 1,979 records, citation blocks, complete HTML,
and LaTeX match; PDF text matches, while PDF bytes differ between builds.
