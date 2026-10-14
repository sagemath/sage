# Duplicate references and URL review, 14 October 2026

This follow-up applies the remarks on
[#42968](https://github.com/sagemath/sage/pull/42968) to **#42966 only**.
The cleanup branch and PR #42968 are unchanged.

## Results

- Consolidated 24 duplicate entries, leaving 1,955 records (1,855 plus 100).
  [The key map](deduplicated-keys.json) records every retired key and its
  replacement. All affected citations in 34 Python/Cython files were updated.
- Checked 1,387 distinct URLs extracted from all fields, including annotations,
  of the original 1,979 records. The initial pass used four concurrent requests,
  followed redirects, verified TLS certificates, and tried HTTPS for HTTP links.
  HEAD requests rejected with 400/405/501 were retried with GET. Failures received
  targeted retries, searches, and publication/institutional checks.
- Changed 433 distinct URLs; [url-changes.tsv](url-changes.tsv) records the old
  and new addresses, canonical citation keys, and basis for each change.
  Both ordinary fields and rendered annotations were updated.
  [url-relocations.tsv](url-relocations.tsv) contains only the 36 substantive
  replacements, excluding 81 HTTP-to-HTTPS-only upgrades and 316
  `dx.doi.org` → `doi.org` normalizations with unchanged DOI identifiers.
  The `ST94` and `ST1994` annotations retain stable DOI links; publisher
  redirect targets inadvertently recorded in the first audit were removed.
- Kept 594 unique DOI fields. Duplicate records contribute their complementary
  metadata to the retained entry. DOI values are not removed just because a
  publisher denies automated access.
- The offline metadata check now rejects repeated DOI values as well as repeated
  keys. Exact-title, normalized-title, author, DOI, and arXiv comparisons were
  used to identify candidates; titles alone were not grounds for merging.

## Reviewer remarks

| Reference | Applied correction |
| --- | --- |
| `Deli2012` | Use the working [Bourbaki PDF](https://www.bourbaki.fr/TEXTES/1048.pdf); the former ENS hostname times out over HTTPS here. |
| `Delsarte48` | Use [JSTOR over HTTPS](https://www.jstor.org/stable/1969047). Crossref/JSTOR identify volume **49**, not 45; add DOI `10.2307/1969047`. |
| `HKP2015` | Use the suggested [DOI](https://doi.org/10.37236/5026), verified against its author/title metadata, and populate the journal/volume/issue/year fields. |
| `HT04` | Use HTTPS for the [MuPAD-Combinat site](https://mupad-combinat.sf.net/). |
| `Kauf1999` | Use the author's [HTTPS PDF](https://homepages.math.uic.edu/~kauffman/VKT.pdf). |
| `Las` | Replace the dead Marne-la-Vallée host with the [Lascoux archive PDF](https://alain-lascoux.github.io/files/lascoux__ProcCrac.pdf). Its title page gives September 2000 and “Youngs's Representations of the Symmetric Group”. |
| `SchaThi1994` | Consolidate into `ST94`, retaining full author names and the verified [DOI](https://doi.org/10.1006/aima.1994.1019); remove the obsolete FTP address. |
| `StaCat98` | Use [math.mit.edu over HTTPS](https://math.mit.edu/~rstan/ec/catalan.pdf); the old `www-math.mit.edu` name does not match the HTTPS certificate. |
| Stanley, *Enumerative Combinatorics*, volume 2 | Consolidate `STA` and `Sta-EC2` into `EnumComb2`. Preserve the 1999 publication and note the 2001 paperback printing, both confirmed by [Stanley's page](https://math.mit.edu/~rstan/ec/). Volume 1 and the separate Catalan exercise excerpt remain distinct. |

## Other duplicates and distinctions

The key map also consolidates repeated Cohen, Haglund, Haiman, Hanani,
Gianni–Trager, Lam, Lam–Schilling–Shimozono, Leinster, Lothaire, Macdonald,
Malvenuto–Reutenauer, Morse, Okado–Schilling–Shimozono,
Okado–Sakamoto–Schilling, Hivert–Schilling–Thiéry, Krob–Thibon,
Haglund–Luoto–Mason–van Willigenburg, Rodríguez Villegas, and TIDES references.
An undated, incomplete `LLMSSZ` reference is consolidated into `LLMSSZ2013`.

Similar titles with different authors, volumes, editions, or explicitly cited
versions remain separate. Examples include `GR2018v5sol`/`GriRei18` (different
arXiv versions and a solutions supplement), `Coh1993`/`Cohen2` (different printings),
`Lot1983`/`Lot1997` (different editions/publishers), `Rea2009`/`Read2011` and
`Reg2005`/`Reg09` (conference/journal versions), and preprint/publication pairs
such as `Ar2006`/`Arm06`, `BMFPR`/`BMFPR2011`, `Cha18`/`Cha2020`,
`NCSF`/`NCSF1`, `LT1998`/`LT2000`, and `BZ05`/`BerZab05`.

`BJL99` is used for both volumes of *Design Theory*, including chapter IX of
volume II; it is not interchangeable with the volume-I-only `BJL99-1`.
The mistakenly attached volume-I DOI was removed from the general reference
and remains on `BJL99-1`.

`RW2008` and `RW2012` are different papers. The latter mistakenly repeated the
former title's “smooth points”; its [arXiv record](https://arxiv.org/abs/1009.5715)
and [journal PDF](https://combinatorialpress.com/article/ojac/vol6/Raichev_2011.pdf)
say **“multiple points”**. Both the title field and displayed annotation were
corrected.

## Moved and malformed links

Replacements include six malformed arXiv identifiers (missing legacy `math/`
prefixes, or an extraneous prefix on a modern identifier), the retired UC Davis
arXiv frontend, Princeton technical reports, the La Jolla Covering Repository,
DesignTheory.org's institutional mirror, Beezer's textbook, a Leiden thesis,
Morrison's manuscript, the HAL record of INRIA report 1779, CWI reports,
Malvenuto–Reutenauer's article, and the Rose-Hulman journal repository.
The redundant `sudokumin.php` suffix in an existing Wayback URL was removed.
Other archive timestamps and embedded HTTP targets were preserved.

The [NIST PDF](https://nvlpubs.nist.gov/nistpubs/legacy/sp/nistspecialpublication800-17.pdf)
replaces the failing MOVS DOI hyperlink; the publisher-confirmed DOI remains
in the `doi` field. PDF replacements were checked for the intended publication;
an HTTP 200 response alone does not establish bibliographic identity.

## Limits and unresolved checks

[unverified-urls.tsv](unverified-urls.tsv) records 168 remaining inconclusive or
error responses. Of these, 128 are HTTP 403 responses, often publisher bot
protection. They are **not** classified as broken publications. Other cases
include timeouts, certificate errors, and unavailable archived snapshots.
Some sites still work over HTTP but fail over HTTPS; those were left unchanged.
In particular, 26 HTTP addresses remain rather than replacing working resources
with unverified HTTPS addresses.

Eleven addresses still returned 404/410 without a verified replacement being
applied: `Beck`, `Ha2005`, `Oha2011`, `CHK2001`, `BHKP2008`, `AHK1994`,
`Connell1999`, `Fogg`, `Wam1999b`, `KG2016`, and `Roberts2015`. The two ACM `10.5555` links
need further investigation. A later Rice technical report and an unrelated
HAL record linked from Karpman's page were not substituted for the cited works.
The La Jolla replacement was readable via the web browser tool but timed out
in local curl; that limitation is included in the TSV.

These are observations on 14 October 2026, not a guarantee of permanent link
availability or an assertion that every successful response was manually read.
No network checks are required when building documentation.

## Validation

The offline metadata/formatting audit passes for all 1,955 entries: 1,948 format
without annotations; the same seven explicitly undated chapters retain their
existing exceptions. All 594 DOI fields are unique. Every retired citation key
has an existing canonical target, and no old citation labels remain in tracked
Python/Cython/RST sources. Executable Python ASTs and Cython non-string tokens
are unchanged in all 34 edited source files.

The catalogue comparison retains anchor checks and uses the explicit key map.
Unchanged annotations must produce identical text, links, and formatting;
changes in rendered annotations are accounted for by the reviewed URL/metadata
corrections and duplicate consolidation. The comparison fixtures document each
expected field difference. Real Meson HTML/PDF catalogue builds and all 48
documentation regression tests passed. The final PDF has 97 pages, with no
missing-character or undefined-reference warnings; its corrected titles and
replacement hyperlinks were also checked directly.

The metadata checker accepts new entries without annotations and rejects both
case-insensitive duplicate DOIs and annotations copied verbatim into `note`.
