# Structured bibliography audit

This records the audit completed on 13 October 2026. The subsequent
[duplicate and URL audit](URL_AUDIT.md) updates the counts and comparison
fixtures, and applies the review remarks from #42968 only to #42966.

This audit follows the initial conversion at `28cf1e290b0`. That conversion
preserved rendered descriptions in `annotation`, but those descriptions hid
incomplete or incorrectly assigned ordinary BibTeX fields. In particular, 453
records had neither `author` nor `title`, and 465 lacked an author. Successful
Sphinx rendering alone did not establish that the structured fields were useful.

## Scope and method

All 1,979 records were parsed and checked for title/author separation, malformed
names, truncated mathematical titles, HTML/MathML in TeX fields, entry types,
publication fields, and ordinary formatting without `annotation`. Ambiguous
boundaries and suspicious matches were inspected against the existing reference
descriptions. Publisher metadata was used to check and supplement the records;
search matches were accepted individually, not simply by similarity score.

Book reviews were rejected as metadata for the books they review. Chapters were
not substituted for whole books, and later editions and journal versions were
not substituted for cited editions, preprints, or conference papers. Online-first
dates were distinguished from the dates of the cited print volumes. Dates were
not inferred from citation keys.

This is a field-structure audit with targeted bibliographic verification, not a
claim that every fact in all 1,979 publications has been independently verified.
Initials and incomplete author lists from the source are retained unless there
is evidence to expand them. The preserved descriptions also retain information
that does not fit the ordinary fields, such as secondary references and errata.

## Results

- Every record has a structured `title`.
- Named people are in `author` or `editor`; 20 authorless web resources are
  explicitly listed in `metadata-exceptions.json`.
- There are 597 DOI fields: 207 added in this audit, and five existing identifiers
  corrected to identify the cited work or version. DOI retrieval is not required
  at documentation build time.
- Ordinary `unsrt` formatting is exercised after removing `annotation`.
  It succeeds for 1,972 records.
  Seven chapter citations omit an edition/publication date: the three
  Abramowitz–Stegun chapter links and four continuously maintained DLMF chapter
  links. Their missing years are explicit exceptions. Undated miscellaneous web
  resources and preprints remain undated where the source supplies no date.
- Citation keys and the 1,879/100 division between the two files are unchanged.
- Only the 16 source corrections below change `annotation`. All other changes
  concern structured metadata and leave Sage's displayed descriptions intact.

[`metadata-changes.tsv`](metadata-changes.tsv) lists every affected key, its
changed fields, and any added/corrected DOI. The Git diff contains the actual
before/after values. Run `check-bibtex.py --tsv fields.tsv` to export all current
author, editor, title, and publication fields for further inspection.

## Corrections to the source bibliography

These belong to the RST cleanup PR #42968 as well as the BibTeX follow-up #42966.
Keys are preserved even when a corrected date differs from the key.

| Key | Correction and evidence |
| --- | --- |
| `AH2002` | The chapter author is Bernhard von Stengel; Aumann and Hart are the volume editors, and Elsevier is the publisher. [Author's chapter](https://ranger.uta.edu/~weems/NOTES6319/PAPERSONE/vonStengel.pdf). |
| `BKLPPRSV2007` | The DOI identified a 2011 encyclopedia entry rather than the 2007 PRESENT paper. Use `10.1007/978-3-540-74735-2_31`; LNCS volume 4727, not 7427. [Original paper, first page](https://crypto.orange-labs.fr/papers/ches2007-450.pdf). |
| `CHPSS18` | This is a EUROCRYPT 2018 paper, not an article in ToSC 2017(4). [Author's university record](https://pure.royalholloway.ac.uk/en/publications/boomerang-connectivity-table-a-new-cryptanalysis-tool/), [DOI](https://doi.org/10.1007/978-3-319-78375-8_22). |
| `Ive2012` | Remove the unrelated Lapointe–Morse paper title accidentally prepended to Iveson's thesis title. [Thesis](https://escholarship.org/uc/item/7pd1v1b5). |
| `Jack1970` | The title says “symmetric polynomials”, not “symmetric functions”. [Publisher DOI](https://doi.org/10.1017/S0080454100008517). |
| `JS2010` | The old DOI identified a rat-gait paper in *Behavioural Brain Research*. Correct it to Jones–Schilling, [Journal of Algebra](https://doi.org/10.1016/j.jalgebra.2010.07.041). |
| `KS2012` | Separate Kelleher–O'Sullivan's preprint from Merca's unrelated journal metadata. Retain the intended [2009 preprint](https://arxiv.org/abs/0909.2331), with DOI `10.48550/arXiv.0909.2331`. |
| `Lim2001` | The cited CRYPTON paper appeared at FSE 1999, not FSE 2001. [Publisher DOI](https://doi.org/10.1007/3-540-48519-8_3). |
| `MAC` | The 1995 second edition, ISBN 0-19-853489-2, is *Symmetric Functions and Hall Polynomials*, not *Symmetric Functions and Orthogonal Polynomials*. [Publisher's review and ISBN](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/6A7C2A63691B1FD7F9EF7E61C25153B5/S0013091500023592a.pdf/i_g_macdonald_symmetric_functions_and_hall_polynomials_2nd_edition_clarendon_press_oxford_1995_x_475_pp_0_19_853489_2_55.pdf). |
| `Motsak2006` | Correct Oleksandr's spelling and identify the work as a diploma thesis, not a Ph.D. thesis. [Title page](https://kluedo.ub.rptu.de/frontdoor/deliver/index/docId/2308/file/Thesis.pdf). |
| `NT2007` | The existing DOI is correct, but the journal, volume, pages, and year were not: *J. Combinatorial Theory B* 98(5), 1038–1055 (2008). [Publisher](https://www.sciencedirect.com/science/article/pii/S0095895608000038). |
| `Rees00` | Restore the missing final variable `n` in the title. Populate author/title and publication fields for both Rees entries. [2000 paper](https://doi.org/10.1006/jcta.1999.3031), [1993 paper](https://doi.org/10.1002/jcd.3180010104). |
| `Reg09` | Replace the STOC 2005 DOI with the DOI of the cited JACM 2009 version, [10.1145/1568318.1568324](https://doi.org/10.1145/1568318.1568324). |
| `Sch2013` | The old DOI identified a Robbins paper from 1939. Use Schmidt's [2013 connectivity paper](https://doi.org/10.1016/j.ipl.2013.01.016). |
| `Sil2007` | Correct “Dynamics Systems” to “Dynamical Systems”. [Author's book page](https://www.math.brown.edu/johsilve/ADSHome.html). |
| `Ste1996` | The first page is 353, not 355. [Publisher DOI](https://doi.org/10.1023/A:1022452717148). |

Institution/date additions to incomplete thesis metadata use the cited title
pages or institutional records, including [Pons](https://perso.lisn.upsaclay.fr/pons/pages/publications.html),
[Molnar](https://people.math.sfu.ca/~numthry/),
[de Cannière](https://cosicdatabase.esat.kuleuven.be/backend/publications/files/these/139),
[Leroux](https://defeo.lu/students), and
[Neurohr](https://oops.uni-oldenburg.de/3607/1/neueff18.pdf).

Some publisher titles expand or phrase the title differently from the retained
description. For `Fil2017`, the structured author, Ivana Protrka, follows the
[published PDF's title page](https://doi.org/10.5592/CO/CCD.2016.05); the migrated
description retains the name Ivana Filipan. Such differences are not treated as
proof of a mismatched DOI without checking the actual work.

## Reproducing the checks

```sh
python build/test/rst2bib/check-bibtex.py src/doc/references/*.bib \
    --output metadata-check.json --tsv fields.tsv
```

This checks the serialized files, including case-insensitive duplicate keys,
not just the in-memory records used during conversion. An exact field/person/type
round trip was also checked when generating the corrections. The checker rejects
unexpected formatting failures and stale exception entries.

The HTML comparison fixtures in this directory distinguish intended source
corrections from conversion changes. At this stage, the RST cleanup catalogue
and the BibTeX catalogue differed only in the previously documented `VT2001`
mathematical-symbol encoding. The fixtures now also include the subsequent
duplicate and URL corrections. See the [comparison instructions](README.md) and
[migration report](../../../BIBTEX_MIGRATION_REPORT.md) for build commands and
rendering validation.

After these corrections, the real Sage Meson HTML and PDF catalogue targets
passed, including the reference-manual inventories. The PDF contains 98 pages;
its LaTeX log has no missing characters or undefined references. All 48
documentation regression tests passed. HTML comparisons checked all 1,979
entries against the corrected RST catalogue and the pre-audit BibTeX output;
the observed differences are exactly those described above.
