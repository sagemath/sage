r"""Shared BibTeX references, including Sage's existing ``[key]_`` citations.

Each independently built manual includes its own cited references.  In
particular, LaTeX never has to refer to a bibliography in another PDF.
The reference catalogue includes all entries and retains their citation keys.
"""

import codecs
import json
import re
from pathlib import Path

import latexcodec  # noqa: F401 -- register the ulatex codec
from docutils import nodes
from docutils.utils import new_document
from pybtex.database.input.bibtex import Parser
from pybtex.plugin import register_plugin
from pybtex.richtext import HRef, Tag, Text
from pybtex.style import FormattedEntry
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.style.labels import BaseLabelStyle
from pybtex.style.sorting import BaseSortingStyle
from sphinx.addnodes import pending_xref
from sphinx.errors import ExtensionError
from sphinx.transforms import SphinxSmartQuotes
from sphinxcontrib.bibtex.citation_target import CitationTarget
from sphinxcontrib.bibtex.directives import (
    BibliographyDirective as BibtexBibliographyDirective,
)
from sphinxcontrib.bibtex.domain import BibtexDomain
from sphinxcontrib.bibtex.roles import CitationRef


class KeyLabels(BaseLabelStyle):
    """Display existing citation keys, rather than generated numeric labels."""

    def format_labels(self, entries):
        return (entry.key for entry in entries)


class KeySorting(BaseSortingStyle):
    """Keep the catalogue sorted case-insensitively by key."""

    def sort(self, entries):
        return sorted(entries, key=lambda entry: entry.key.upper())


class DatabaseParser(Parser):
    """Validate keys before Pybtex's case-insensitive cache drops duplicates."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.seen = set()
        self.data.want_entry = lambda key: True

    def process_entry(self, entry_type, key, fields):
        if key.upper() in self.seen:
            raise ExtensionError('bibliography contains case-insensitive duplicate keys')
        self.seen.add(key.upper())
        super().process_entry(entry_type, key, fields)

    def handle_error(self, error):
        raise ExtensionError(f'invalid bibliography: {error}')


class BibliographyDirective(BibtexBibliographyDirective):
    """Keep citation anchors from colliding with headings or local citations.

    sphinxcontrib-bibtex reserves IDs for every database entry before it knows
    which entries are cited. Keep the traditional IDs when possible, notably
    in the reference catalogue, and prefix them when this document already
    uses one of those IDs.
    """

    def run(self):
        config = self.env.config
        original = config.bibtex_cite_id
        template = original
        count = self.env.temp_data.get('bibtex_bibliography_count', 0) + 1
        keyprefix = self.options.get('keyprefix', '')
        entries = self.env.get_domain('cite').bibdata.data.entries
        while any(nodes.make_id(template.format(bibliography_count=count, key=keyprefix + key))
                  in self.state.document.ids for key in entries):
            template = f'bibtex-{count}-' + template
        config.bibtex_cite_id = template
        try:
            return super().run()
        finally:
            config.bibtex_cite_id = original


def _group(text, start):
    """Read a balanced TeX group, including escaped literal braces."""
    depth = 1
    end = start
    while end < len(text):
        if text[end] == '\\':
            end += 2
            continue
        if text[end] == '{':
            depth += 1
        elif text[end] == '}':
            depth -= 1
            if not depth:
                return text[start:end], end + 1
        end += 1
    raise ValueError('unbalanced braces in bibliography annotation')


def annotation_text(text):
    r"""Render the TeX in a migrated citation without losing maths or links.

    Pybtex's general-purpose LaTeX decoder handles accents, but not TeX
    formatting commands, hyperlinks, or mathematical expressions.  Migrated
    annotations use only these standard commands; new entries can instead
    use ordinary BibTeX fields and the usual Pybtex formatting.
    """
    parts = []
    position = 0
    pattern = re.compile(r'\\(emph|textbf|texttt|href|url)\{|\\\(')
    while match := pattern.search(text, position):
        parts.append(_decode(text[position:match.start()]))
        command = match[1]
        if command is None:
            end = text.find(r'\)', match.end())
            if end == -1:
                raise ValueError('unclosed maths in bibliography annotation')
            parts.append(Tag('math', text[match.end():end]))
            position = end + 2
        else:
            content, position = _group(text, match.end())
            if command == 'href':
                label, position = _group(text, position + 1)
                parts.append(HRef(content, annotation_text(label)))
            elif command == 'url':
                parts.append(HRef(content, content))
            else:
                tag = {'emph': 'em', 'textbf': 'strong', 'texttt': 'tt'}[command]
                parts.append(Tag(tag, annotation_text(content)))
    parts.append(_decode(text[position:]))
    return Text(*parts)


def _decode(text):
    # Keep whitespace at the boundaries of formatting commands.  A TeX
    # decoder otherwise treats the leading spaces of each fragment as the
    # leading spaces of a whole document and discards them.
    literals = []

    def protect(match):
        command = match[0]
        literals.append({
            r'\textasciitilde{}': '~',
            r'\textasciicircum{}': '^',
            r'\textbackslash{}': '\\',
            r'\textquotesingle{}': "\x00'",
            r'\textquotedbl{}': '\x00"',
        }.get(command, command[1:]))
        return chr(0xF0000 + len(literals) - 1)

    text = re.sub(r'\\(?:textasciitilde|textasciicircum|textbackslash|textquotesingle|textquotedbl)\{\}|\\[{}%&#_$]',
                  protect, text)
    text = codecs.decode('x' + text + 'x', 'ulatex')[1:-1]
    for index, literal in enumerate(literals):
        text = text.replace(chr(0xF0000 + index), literal)
    return text


class SageStyle(UnsrtStyle):
    """Preserve migrated descriptions, and format new structured entries."""

    default_label_style = 'sage'
    default_sorting_style = 'sage'

    def format_entry(self, label, entry, bib_data=None):
        if 'annotation' in entry.fields:
            return FormattedEntry(entry.key, annotation_text(entry.fields['annotation']), label)
        return super().format_entry(label, entry, bib_data)


def validate_database(app):
    """Reject keys which would collide in the legacy citation index."""
    parser = DatabaseParser(encoding=app.config.bibtex_encoding)
    parser.parse_files(app.env.get_domain('cite').bibdata.bibfiles)


def add_bibliography(app, docname, source):
    """Give standalone manuals a bibliography of their cited entries."""
    if (not app.config.sage_bibtex_auto_bibliography
            or getattr(app.config, 'multidocs_subdoc_list', [])
            or docname != app.config.root_doc
            or re.search(r'^\s*\.\. bibliography::', source[0], re.M)):
        return
    source[0] += '\n\n.. rubric:: Bibliography\n\n.. bibliography::\n'


def bibliography_signature(app):
    """Identify the bibliography version used by the current builder."""
    return {str(path): bibfile.mtime
            for path, bibfile in app.env.get_domain('cite').bibdata.bibfiles.items()}


def bibliography_stamp(app):
    """Keep output-format stamps beside Sphinx's other cached build data."""
    return Path(app.doctreedir) / f'sage-bibtex-{app.builder.name}.json'


def outdated_bibliography_documents(app, env, added, changed, removed):
    """Rebuild citations after another builder refreshes the shared cache.

    Sage builds inventories before HTML with the same doctrees. An inventory
    pass can consume a bibliography-only change without updating HTML, so
    compare with the bibliography used by this output format's last build.
    """
    try:
        previous = json.loads(bibliography_stamp(app).read_text())
    except (FileNotFoundError, ValueError):
        previous = None
    if previous == bibliography_signature(app):
        return []
    domain = env.get_domain('cite')
    documents = {key.docname for key in domain.bibliographies}
    documents.update(reference.docname for reference in domain.citation_refs)
    documents.update(domain.data.get('sage_legacy_refs', {}))
    return sorted(documents)


def record_bibliography_build(app, exception):
    """Record the bibliography only after a successful build of this format."""
    if exception is None:
        bibliography_stamp(app).write_text(json.dumps(bibliography_signature(app)))


def smartquote_bibliography(app, doctree, docname):
    """Apply Sphinx's typography to entries generated after its transforms.

    Native citations have already passed through SphinxSmartQuotes. BibTeX
    entries are inserted by a post-transform, so process only those entries,
    respecting the language and the builder's smartquotes configuration.
    """
    domain = app.env.get_domain('cite')
    generated = {(key.docname, node['ids'][0])
                 for key, bibliography in domain.bibliographies.items()
                 for node in bibliography.citation_nodes.values()}
    citations = [node for node in doctree.findall(nodes.citation)
                 if (node.get('docname'), node['ids'][0]) in generated]
    document = new_document(docname, settings=doctree.settings)
    document.extend(node.deepcopy() for node in citations)
    SphinxSmartQuotes(document, None).apply()
    for original, formatted in zip(citations, document.children):
        original[:] = formatted.children


def collect_legacy_references(app, doctree):
    """Include native reST citations in sphinxcontrib-bibtex's cited set."""
    domain = app.env.get_domain('cite')
    entries = domain.bibdata.data.entries
    references = []
    for node in doctree.findall(pending_xref):
        if node.get('refdomain') != 'citation':
            continue
        key = node['reftarget']
        if key not in entries:
            continue
        # A distinct id identifies these records when excluding locally
        # defined reST citations below.  It also provides a real backlink.
        citation_id = f'sage-citation-{len(references)}'
        while citation_id in doctree.ids:
            citation_id += '-1'
        node['ids'].append(citation_id)
        doctree.ids[citation_id] = node
        references.append(CitationRef(
            citation_id, app.env.docname, node.line or 0,
            [CitationTarget(key=key, pre='', post='')],
        ))
    domain.data.setdefault('sage_legacy_refs', {})[app.env.docname] = references


def purge_legacy_references(app, env, docname):
    env.get_domain('cite').data.get('sage_legacy_refs', {}).pop(docname, None)


def merge_legacy_references(app, env, docnames, other):
    references = env.get_domain('cite').data.setdefault('sage_legacy_refs', {})
    incoming = other.domaindata['cite'].get('sage_legacy_refs', {})
    for docname in docnames:
        if docname in incoming:
            references[docname] = incoming[docname]


def prefer_local_citations(app, env):
    """Keep references defined in docstrings local, even for shared keys."""
    domain = env.get_domain('cite')
    native = env.get_domain('citation')
    # Remove only entries we published ourselves on the previous build.
    # Otherwise an incremental build would mistake them for docstring
    # definitions, and stop including their references in the bibliography.
    for key, value in native.data.get('sage_bibtex_citations', {}).items():
        if native.citations.get(key) == value:
            del native.citations[key]
    # multidocs also loads citations from *other* manuals into this domain.
    # They must not prevent a local bibliography (especially in a PDF).
    local = {key for key, value in native.citations.items()
             if value[0] in env.found_docs}
    domain.citation_refs[:] = [reference for reference in domain.citation_refs
                              if not reference.citation_ref_id.startswith('sage-citation-')]
    domain.citation_refs.extend(
        reference for references in domain.data.get('sage_legacy_refs', {}).values()
        for reference in references if reference.targets[0].key not in local
    )


def publish_legacy_citations(app, env):
    """Publish BibTeX entries to the citation index used by Sage's manuals.

    sphinxcontrib-bibtex creates citation nodes during post-transforms, after
    Sphinx normally collects the native citation domain.  Register them once
    its labels are known, before multidocs merges or saves that domain.
    """
    domain = env.get_domain('citation')
    published = {}
    for citation in env.get_domain('cite').citations:
        key = citation.key
        value = (citation.bibliography_key.docname, citation.citation_id, 0)
        if key not in domain.citations or domain.citations[key][0] not in env.found_docs:
            domain.citations[key] = value
            published[key] = value
        # A catalogue deliberately includes uncited entries.  This does not
        # suppress warnings about misspelled or otherwise missing citations.
        domain.citation_refs.setdefault(key, set())
    domain.data['sage_bibtex_citations'] = published


def setup(app):
    register_plugin('pybtex.style.labels', 'sage', KeyLabels)
    register_plugin('pybtex.style.sorting', 'sage', KeySorting)
    register_plugin('pybtex.style.formatting', 'sage', SageStyle)
    BibtexDomain.backend.tags['math'] = nodes.math
    app.setup_extension('sphinxcontrib.bibtex')
    app.add_directive('bibliography', BibliographyDirective, override=True)
    app.add_config_value('sage_bibtex_auto_bibliography', True, 'env')
    app.connect('builder-inited', validate_database, priority=600)
    app.connect('source-read', add_bibliography)
    app.connect('env-get-outdated', outdated_bibliography_documents)
    app.connect('build-finished', record_bibliography_build)
    app.connect('doctree-resolved', smartquote_bibliography)
    app.connect('doctree-read', collect_legacy_references)
    app.connect('env-purge-doc', purge_legacy_references)
    app.connect('env-merge-info', merge_legacy_references)
    app.connect('env-updated', prefer_local_citations, priority=400)
    app.connect('env-updated', publish_legacy_citations)
    return {'env_version': 2, 'parallel_read_safe': True, 'parallel_write_safe': True}
