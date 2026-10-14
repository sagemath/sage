# Sage BibTeX migration: changes and validation

The bibliography migration moves module bibliographies into shared BibTeX files and provides
local bibliographies for independently built manuals. The original two-file split
produced the same catalogue HTML and LaTeX as the preceding single-file version;
the subsequent semantic audit and its corrections are recorded below.
The [14 October duplicate and URL audit](build/test/rst2bib/URL_AUDIT.md)
consolidates 24 repeated entries, updates their citations, and repairs or upgrades
433 URLs in #42966. It leaves #42968 unchanged. Other citation keys remain stable.

The initial validation covered code through commit
[`2b440e3f44c`](https://github.com/dimpase/sage/commit/2b440e3f44c), on
[`docs/bibtex-references`](https://github.com/dimpase/sage/tree/docs/bibtex-references).
The original documentation baseline is Sage 11.0.beta2, commit
[`0835358cba6`](https://github.com/sagemath/sage/commit/0835358cba6).
Validation was performed on 12–13 October 2026.

## Separate review of cleanup and conversion

The work is now split into two pull requests:

1. [#42968](https://github.com/sagemath/sage/pull/42968) corrects inconsistent
   references and moves module bibliographies into the existing RST database.
   Separate commits distinguish corrections from relocation.
2. [#42966](https://github.com/sagemath/sage/pull/42966) converts that consolidated
   database to BibTeX and supplies the Sphinx integration, dependency/build
   metadata, developer documentation, and comparison tools.

The original-baseline comparisons below cover both parts together. The cleanup
was additionally tested independently with the unchanged RST generator and
Sphinx 9.1.0. Its catalogue and three standalone module pages build successfully.
All 1,807 original shared entries remain. After the semantic audit described
below, 1,792 render identically in the RST cleanup; 15 have documented corrections.
Comparison of the 177 relocated definitions finds six documented changes.
All 546 example blocks and the module prose remain unchanged after
excluding the relocated bibliography and normalizing renamed keys. Executable
AST/token checks also pass for all 62 edited source files.

Before the duplicate and URL audit, comparing the cleanup's 1,979-entry RST
catalogue with the BibTeX catalogue found 1,978 identical entries and only the
documented `VT2001` mathematical-symbol difference. The current comparison maps
all 1,979 old keys to 1,955 canonical entries: 1,835 descriptions are identical
and 144 have documented changes. The exact differences are recorded in
[`conversion-differences.json`](build/test/rst2bib/conversion-differences.json).
The semantic audit adds structured metadata corrections without changing the
conversion implementation.

## Bibliography contents and preservation

| Item | Result |
| --- | --- |
| Original shared citation keys | All 1,807 accounted for, including documented duplicate-key mappings |
| Bibliographic definitions migrated from docstrings | 177, from 60 Python and two Cython files |
| New shared entries | 172; matching existing references reused |
| Final shared entries | 1,955 (24 duplicates consolidated after the initial 1,979-entry migration) |
| Final DOI fields | 594, all unique; the earlier semantic audit had 597 fields before consolidation |
| Original DOI fields | 364; corrected identifiers are recorded below and in the semantic audit |

Migrated entries have an `annotation` field containing the complete original
description in TeX, including mathematical expressions, formatting, notes, and
hyperlinks. Sage renders that description. The semantic audit supplies titles for every record and corrects the ordinary
author/editor and publication fields; new entries can use those fields without an
`annotation`.

DOI metadata was retrieved during migration to supply publication details.
Bibliography rendering uses the checked-in files and makes no metadata requests.
Descriptions and hyperlinks were checked both before conversion and after parsing
the serialized BibTeX. The subsequent HTML comparisons are described below.

The executable ASTs of all 60 edited Python files and the executable token streams
of both edited Cython files are unchanged when docstrings are excluded. Five
numbered bibliographic footnotes became named citations. Two malformed
`.. [GriRei18]_` definitions became citations. The two explanatory footnotes in
`src/sage/combinat/tutorial.py` remain local footnotes.

## Conflicting citation keys

Three local keys referred to different works from the same keys in the shared
database. The shared keys retain their original meanings; only citations to the
conflicting local works were changed.

| Local key | Meaning of the existing shared key | Local work and resulting key |
| --- | --- | --- |
| `GL1996` | Golub–van Loan, *Matrix Computations* | Graham–Lehrer, *Cellular algebras*: reused existing `GrLe1996` |
| `BK2001` | Bruns–Koch, *Computing the integral closure of an affine semigroup* | Bandlow–Killpatrick, *An area-to_inv bijection between Dyck paths and 312-avoiding permutations*: `BandlowKillpatrick2001` |
| `CK2008` | Corneil–Krueger, *A Unified View of Graph Searching* | Claesson–Kitaev, *Classification of bijections between 321- and 132-avoiding permutations*: `ClaessonKitaev2008` |

These changes are in the
[docstring-bibliography migration commit](https://github.com/dimpase/sage/commit/5f9bc16e632).
Installed-wheel HTML checks verify that all six works resolve to their respective
descriptions.

## Corrected identifiers and added information

| Entry | Original problem or missing information | Treatment |
| --- | --- | --- |
| `OSShimo03` | `10.1.1.192.2095` was a CiteSeer identifier labelled as a DOI; arXiv `0810.5067` referred to a different paper | Replaced with DOI `10.1090/S1088-4165-03-00155-9` and arXiv `math/0105017` |
| `GrLe1996` | Shared description abbreviated Graham's initials as `J.`, whereas the matching local reference used `J.J.`; publisher DOI was absent | Displayed initials changed to `J.J.`; DOI `10.1007/BF01232365` added; existing MathSciNet link retained |
| `BandlowKillpatrick2001` | Publisher DOI was absent from the local description | Added DOI `10.37236/1584` |
| `CS2022` | A sentence-ending full stop had been included in the DOI identifier | Changed `10.1016/j.jcta.2021.105572.` to `10.1016/j.jcta.2021.105572`, placing punctuation outside the hyperlink |
| `VT2001` | A Unicode mathematical italic *n* caused a PDF font problem | Encoded the symbol as inline mathematics |
| `Rees00`, `Rees93` | The author followed the paper title | Moved Rolf S. Rees before the title in both the RST cleanup and BibTeX descriptions |

For `OSShimo03`, the intended work is the Okado–Schilling–Shimozono paper in
*Representation Theory* 7 (2003), 101–163, identified by
[arXiv:math/0105017](https://arxiv.org/abs/math/0105017) and its
[journal DOI](https://doi.org/10.1090/S1088-4165-03-00155-9).
The former [arXiv:0810.5067](https://arxiv.org/abs/0810.5067) instead identifies
Fourier–Okado–Schilling's *Kirillov–Reshetikhin crystals for nonexceptional types*.
The bibliography contains an
[explicit correction comment](https://github.com/dimpase/sage/blob/2b440e3f44c/src/doc/references/reference-database.bib#L14410)
recording both replaced identifiers.

The additional identifiers resolve to the publisher records for
[Graham–Lehrer](https://doi.org/10.1007/BF01232365),
[Bandlow–Killpatrick](https://doi.org/10.37236/1584), and
[Crew–Spirkl](https://doi.org/10.1016/j.jcta.2021.105572).
For `GrLe1996`, both the displayed description and structured author now use
`J. J. Graham`.

## Semantic audit after conversion

The [structured-metadata audit](build/test/rst2bib/SEMANTIC_AUDIT.md) records the
method, limitations, and source evidence for 16 additional corrections. It found
that 453 converted entries had neither an author nor a title; the preserved
`annotation` had hidden these omissions. All 1,979 entries now have titles and
named people in the appropriate author/editor fields. Twenty authorless web
resources are listed explicitly.

The audit corrects five DOIs that identified different works or versions, restores
truncated mathematical titles, separates chapter authors from editors, and fixes
misclassified publication types. It adds 207 DOI fields. The source corrections
are also applied to the preceding RST cleanup PR, so the conversion-only HTML
comparison still has just the existing `VT2001` difference.

The offline `check-bibtex.py` check passes for all entries. With `annotation`
removed, 1,972 entries format using ordinary `unsrt`; seven chapter links omit
a publication date and are explicit exceptions. It also checks serialized
structure, malformed names, HTML in TeX fields, and duplicate keys. Deliberately
malformed inputs confirmed that missing fields hidden by an annotation, duplicate
keys, and initials parsed as family names are rejected. This does not claim
independent verification of every bibliographic fact.

The [per-key change list](build/test/rst2bib/metadata-changes.tsv) identifies every
changed field. Comparing rendered HTML with the pre-audit catalogue finds exactly
the 16 documented source corrections and no other changes, additions, or removals.

## Several BibTeX files

All `.bib` files in [`src/doc/references/`](src/doc/references/) form one shared
citation database. Keys must be unique across all files, including
case-insensitive duplicates. Entry order within a file and the order in which
files are loaded do not determine the rendered order; Sage's bibliography style
sorts the combined references case-insensitively by citation key.

The final database is split as follows:

| File | Entries |
| --- | ---: |
| [`reference-database.bib`](src/doc/references/reference-database.bib) | 1,855 |
| [`reference-database-extra.bib`](src/doc/references/reference-database-extra.bib) | 100 |

The second file contains the last 100 entries of the former single file, from
`VB1996` through `ZZ2005`. At the time of the split, their entry blocks were moved without reserialization,
and all 1,979 parsed records were identical to the single-file version. Subsequent
metadata corrections keep the same file division.

Meson copies the files listed in
[`src/doc/references/meson.build`](src/doc/references/meson.build). The
`sage-docbuild` wheel includes `references/*.bib` through the shared-directory
symlink. Sphinx configuration and Sage's citation scanner read all these files.
The developer guide documents how to add another file.

Both `[Key]_` and `:cite:p:` citation syntax work. Each independently built manual
includes the references it cites, so its HTML and PDF do not require another
manual's citation cache or a separately built reference catalogue.

## Build and dependency metadata

The change updates the Meson documentation build, the `sage-docbuild` wheel's
dependencies and package data, the root documentation dependency group and
`uv.lock`, and all ten existing conda environment files for Python 3.13 and 3.14.
Sage package metadata provides standard wheel packages for `latexcodec`,
`pybtex`, `pybtex_docutils`, and `sphinxcontrib_bibtex`, with their dependency order
and downloaded-wheel checksums checked.

The former RST reference generator is removed. The catalogue uses a static RST
page containing a bibliography directive. Documentation builds no longer generate
that page from `reference-database.rst`.

## Validation environment and commands

Validation used macOS on Apple Silicon, Python 3.14.8, Sphinx 9.1.0,
`sphinxcontrib-bibtex` 2.7.0, and Pybtex 0.26.1. The branch's `sage-docbuild` wheel
was installed into a test environment using the already installed Homebrew native
Sage 11.0.beta1 library. Documentation sources came from the branch based on
11.0.beta2. Sage's native library was not rebuilt for this documentation change.

A genuine top-level Sage Meson build was configured with `-Dbuild-docs=true` and
the host's external dependency paths. The following targets passed, including
their reference-manual inventory prerequisites; `BUILD` denotes that configured
build directory:

```sh
export SAGE_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
python -m mesonbuild.mesonmain compile -C BUILD -j4 \
    doc-html-reference-references doc-pdf-reference-references
```

The PDF target produced a 98-page catalogue. Its LaTeX log contained no missing
characters or undefined references.

All **48 documentation regression tests passed**. The selected modules were
`main_test.py`, `bootstrap_test.py`, and `bibtex_test.py` under `src/sage_docbuild/`.
In the test environment described above, native Sage was imported before adding
the branch source directory to Python's import path:

```sh
export SAGE_SRC="$PWD/src"
export SAGE_DOC_SRC="$SAGE_SRC/doc"
python -c 'import sage.all, sys; sys.path.insert(0, "src"); import pytest; raise SystemExit(pytest.main(["--confcutdir=src/sage_docbuild", "src/sage_docbuild/main_test.py", "src/sage_docbuild/bootstrap_test.py", "src/sage_docbuild/bibtex_test.py", "-q"]))'
```

The bibliography integration tests run with one and two input files and cover
reverse file order, duplicate keys spanning files, parallel reading, missing-key
errors, local RST definitions and their removal, stale citations from another
manual, incremental builds, typography, heading-anchor collisions, and the shared
inventory/HTML doctree cache.

Separate installed-wheel checks verified both bundled files against the source
files by SHA256 and built standalone HTML, LaTeX, and PDF outside the checkout.
They exercised references from both files, the six formerly ambiguous works, and
a real docstring from `q_bernoulli.pyx`. The changed extension, citation scanner,
and tests passed Ruff; `git diff --check` passed. Conda environment files parsed
and their added package pins were checked against conda-forge metadata.

## Comparison with the original HTML

Before/after HTML was built using the original RST bibliography and generator at
the baseline commit, and the migrated bibliography on this branch, with the same
Sphinx and installed native Sage library.

The citation DOM comparison checked anchor IDs, normalized text, links and their
captions, mathematical expressions, emphasis, bold text, and inline code. It
normalized whitespace and excluded backlink labels.

| Comparison | Result |
| --- | --- |
| Original catalogue entries | 1,791 of 1,807 match all checked fields exactly |
| Original catalogue anchors | All 1,807 identical |
| Migrated docstring bibliography definitions | 171 of 177 match all checked fields exactly |
| Complete module-page prose | Identical after excluding relocated bibliography blocks and normalizing the three renamed citation keys |
| Module example blocks | All 546 identical: 112 in `difference_family`, 196 in `dyck_word`, and 238 in `diagram_algebras` |

The 16 differing catalogue entries consist of `GrLe1996`, `VT2001`, and the
14 additional shared-reference corrections listed in the semantic audit. The six
differing docstring definitions are `GrLe1996`, `BandlowKillpatrick2001`,
`OSShimo03`, `Rees00`, `Rees93`, and `MAC`. Exact text, link, formatting, and
mathematics differences are recorded in the comparison fixtures.

This comparison also caught and led to fixes for smart-quote changes, lost escaped
quotes, changed DOI/arXiv link captions, the malformed `CS2022` DOI URL, a GSL
heading/bibliography-anchor collision, and stale HTML after an inventory builder
had consumed a bibliography-only edit.

The comparison evaluates content and document structure. The catalogue's former
alphabetical jump list was removed, and module bibliographies move to the end of
their pages. Those layout changes are excluded from the content comparisons.

## Comparison before and after splitting the BibTeX file

The single-file outputs were rebuilt and saved immediately before the split.
The real Sage Meson HTML and PDF targets were then run with the two-file database.

| Check | Result |
| --- | --- |
| Parsed entry records | All 1,979 identical |
| Rendered citation blocks, including order, anchors, text, and links | All 1,979 identical |
| Entire catalogue HTML | Byte-for-byte identical |
| Generated `references.tex` | Byte-for-byte identical |
| Extracted PDF text | Identical |
| PDF binary | Differs between builds; binary equality is not claimed |

The SHA256 of both complete catalogue HTML files is:

```text
b18adb82bfe4bf9b8722ffdc776549659abb0204c2bb97b9ab4c54116e9a6a14
```

A trailing blank separator in the primary `.bib` file was removed after rendering
to satisfy `git diff --check`. Parsed records and installed-file checksums were
verified again after this whitespace-only edit.

## Scope and retained evidence

The checks cover the complete bibliography, representative independent manuals,
and all required reference inventories. They do not constitute a rebuild of every
Sage HTML/PDF manual or a new native Sage build. The HTML comparisons are DOM and
content comparisons, rather than screenshot comparisons.

The host lacked `xindy`, so a module PDF test used `makeindex`. One earlier
single-file-module PDF run compiled successfully but encountered an existing
`str.glob` error during the wrapper's final export; its generated PDF and LaTeX
log were checked directly. The catalogue PDF built through Meson successfully.
The ten pre-existing Ruff findings in `conf.py` also occur in the baseline version.

Build logs, generated HTML/PDF files, metadata-audit results, and before/after
comparison JSON files are retained in the local validation workspace. They are
not attached to this report. This document publishes their aggregate results and
the correction record without depending on machine-specific file paths.

The reusable comparison commands, exact expected differences, citation-key maps,
and usage instructions are included in
[`build/test/rst2bib/`](build/test/rst2bib/README.md). They accept explicit input
paths and produce JSON reports. Generated documentation outputs remain separate
inputs to these commands.

The implementation history is available in the commits for the
[shared bibliography](https://github.com/dimpase/sage/commit/3d680951d13),
[docstring bibliography migration](https://github.com/dimpase/sage/commit/5f9bc16e632),
[rendering and cache fixes](https://github.com/dimpase/sage/commit/4bb2ecbd043),
[removal of the source-order restriction](https://github.com/dimpase/sage/commit/1835786f11e),
and [multiple-file support](https://github.com/dimpase/sage/commit/2b440e3f44c).
