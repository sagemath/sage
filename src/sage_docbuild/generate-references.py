r"""
Generate the "Bibliographic References" ReST document

This script generates the Sage bibliography from a list of references
(its first argument) and a template (its second argument). The usage
is simply::

  $ python src/sage_docbuild/generate-references.py <references> <infile> <outfile>

In most circumstances,

  <references> = src/doc/reference-database.rst
  <infile>     = src/doc/en/reference/references/index.rst.in
  <outfile>    = src/doc/en/reference/references/index.rst

because this script is executed with those parameters by the build system.

The goal here is to ensure that our reference database remains sorted
and syntactically valid. The "Bibliographic References" page is
free-form ReST, so it is difficult to parse and validate as a
whole. Instead, we keep (only) the references in a dedicated file;
this allows us to easily check whether or not that file is sorted and
valid. If it is, we populate the ``<infile>`` template and write it to
``<outfile>``. If not, an error is raised.

Sorted in this case means grouped by label (because that's how the
template displays them) and case-insensitively within each group.
"""

from collections import defaultdict
from pathlib import Path
from string import ascii_uppercase
from sys import argv

from sage.misc.sagedoc import process_extlinks

reference_db = argv[1]
infile = argv[2]
outfile = argv[3]

def parse_rst(rst):
    r"""
    Parse RST text into a docutils document.
    """
    from docutils.frontend import get_default_settings
    from docutils.parsers.rst import Parser
    from docutils.utils import new_document
    settings = get_default_settings(Parser)
    document = new_document('<rst-doc>', settings=settings)
    Parser().parse(rst, document)
    return document

def indent(s):
    r"""
    Indent all lines of the given string except the first by three
    spaces.

    In RST, the indentation level is used to indicate the end of
    a citation. The three spaces ensure alignment with a leading
    ``.. [label]``.
    """
    return "\n".join("   " + l for l in s.splitlines()).lstrip()

def citation_label(c):
    r"""
    Return the contents (sans brackets) of a citation label.
    """
    return c.children[0].children[0]

def citation_label_upper(c):
    r"""
    Return the contents (sans brackets) of a citation label, in
    uppercase.
    """
    return citation_label(c).upper()

def citation_to_rst(c):
    r"""
    Convert a citation element back to RST.
    """
    return f".. [{citation_label(c)}] {indent(c.rawsource)}"

# Read in the database, and preprocess any nonstandard extlinks
# like ":arxiv:"
doc = parse_rst(process_extlinks(Path(reference_db).read_text()))
sorted_citations = sorted(doc.citations, key=citation_label_upper)

# Citations grouped by the first letter of their label, ignoring case.
groups = defaultdict(list)

# We're also going to build (grouped) lists of citation labels as we
# go, so that we can check them for duplicates.
seen_labels = defaultdict(list)

for (c, d) in zip(sorted_citations, doc.citations):
    # Build the groups, check for mis-ordering, and check for
    # duplicate labels all at once.
    if c != d:
        raise ValueError(
            "reference database is out of order "
            f"({citation_label_upper(c)})"
        )

    # Duplicate labels must be checked case-INsensitively, because the
    # generated HTML uses the lowercase label as the "id" of a <div>.
    # Those ids are used to generate links, which must be unambiguous.
    uclabel = citation_label_upper(c)
    x = uclabel[0]

    if uclabel in seen_labels[x]:
        # Checking for duplicate *citations* would be far too lenient.
        raise ValueError("reference database contains duplicate labels")

    seen_labels[x].append(uclabel)
    groups[x].append(c)


template = Path(infile).read_text()

for x in ascii_uppercase:
    these_refs = "\n\n".join(map(citation_to_rst, groups[x]))
    template = template.replace(f"@REFS_{x}@", these_refs)

Path(outfile).write_text(template)
