"""Integration tests for separately built HTML and LaTeX bibliographies."""

import re
import subprocess
import sys
from html import unescape
from pathlib import Path

import pytest

pytest.importorskip('sphinxcontrib.bibtex')


def build(source, output, builder='html', fresh=False, doctrees=None):
    command = [sys.executable, '-m', 'sphinx', '-b', builder, '-W', '--keep-going',
               '-j', '2', str(source), str(output)]
    if fresh:
        command.insert(3, '-E')
    if doctrees is not None:
        command[3:3] = ['-d', str(doctrees)]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.fixture(params=[False, True], ids=['one-bibfile', 'two-bibfiles'])
def manual(tmp_path, request):
    source = tmp_path / 'source'
    source.mkdir()
    (source / 'refs.bib').write_text(r'''
@article{ExistingKey,
  author = {Gauss, Carl Friedrich},
  title = {{A title with preserved capitals}},
  journal = {Journal of tests},
  year = {1801},
  doi = {10.1234/example_a}
}
@misc{OtherKey,
  annotation = {A. Author. \emph{A title in \(\mathbb{Z}\)}.
                A note, \url{https://example.org/a_b}.}
}
''')
    # pytest can import the source package without installing sage-docbuild.
    # Give Sphinx's separate interpreter the same package location, so these
    # tests exercise this checkout rather than another installed version.
    (source / 'conf.py').write_text(
        f'import sys\nsys.path.insert(0, {str(Path(__file__).resolve().parent.parent)!r})\n' + '''
extensions = ['sage_docbuild.ext.bibtex']
bibtex_bibfiles = ['refs.bib']
bibtex_default_style = 'sage'
bibtex_cite_id = '{key}'
project = 'Bibliography test'
latex_documents = [('index', 'test.tex', 'Bibliography test', 'Sage', 'manual')]
''')
    if request.param:
        bibfile = source / 'refs.bib'
        text = bibfile.read_text()
        cut = text.index('@misc{OtherKey,')
        (source / 'extra.bib').write_text(text[cut:])
        bibfile.write_text(text[:cut])
        # Read OtherKey before ExistingKey: file order must not dictate the
        # rendered order, and duplicate validation must span both files.
        with (source / 'conf.py').open('a') as file:
            file.write("\nbibtex_bibfiles = ['extra.bib', 'refs.bib']\n")
    (source / 'index.rst').write_text('''
A manual
========

An existing citation [ExistingKey]_ and a new citation :cite:p:`OtherKey`.

.. toctree::

   chapter
   parallel1
   parallel2
   parallel3
''')
    (source / 'chapter.rst').write_text('''
A chapter
=========

Repeated citation [OtherKey]_ and :cite:p:`ExistingKey`.
''')
    for number in range(1, 4):
        (source / f'parallel{number}.rst').write_text(
            f'Chapter {number}\n=========\n\n[ExistingKey]_.\n')
    return source


def test_html_citations_and_incremental_build(manual, tmp_path):
    output = tmp_path / 'html'
    build(manual, output)
    page = (output / 'index.html').read_text()
    chapter = (output / 'chapter.html').read_text()
    assert 'href="#existingkey"' in page
    assert 'href="index.html#otherkey"' in chapter
    assert 'https://doi.org/10.1234/example_a' in page
    assert 'https://example.org/a_b' in page
    assert r'\mathbb{Z}' in page
    assert '[ExistingKey]' in page
    assert '[OtherKey]' in unescape(re.sub(r'<[^>]*>', '', page))
    assert page.index('id="existingkey"') < page.index('id="otherkey"')
    # Reuse the environment, then change a citing document. Cached BibTeX
    # entries must not be mistaken for local reST citation definitions.
    build(manual, output)
    with (manual / 'chapter.rst').open('a') as file:
        file.write('\nAn edit, still citing [OtherKey]_.\n')
    build(manual, output)
    assert 'href="index.html#otherkey"' in (output / 'chapter.html').read_text()


def test_separate_manuals_and_latex(manual, tmp_path):
    # No objects.inv or citation cache from the reference catalogue is
    # required, and both old and new citation syntax have local PDF targets.
    build(manual, tmp_path / 'first-html')
    (manual / 'index.rst').write_text('''
Another independent manual
==========================

See [ExistingKey]_ and :cite:p:`OtherKey`.
''')
    for page in manual.glob('*.rst'):
        if page.name != 'index.rst':
            page.unlink()
    build(manual, tmp_path / 'second-html')
    build(manual, tmp_path / 'latex', 'latex')
    latex = (tmp_path / 'latex' / 'test.tex').read_text()
    assert 'ExistingKey' in latex
    assert 'OtherKey' in latex
    assert r'\mathbb{Z}' in latex
    assert r'\bibitem' in latex


def test_missing_citation_is_an_error(manual, tmp_path):
    (manual / 'index.rst').write_text('Title\n=====\n\n[MisspelledKey]_.\n')
    result = subprocess.run(
        [sys.executable, '-m', 'sphinx', '-b', 'html', '-W', str(manual), str(tmp_path / 'html')],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0
    assert 'citation not found: MisspelledKey' in result.stderr


def test_case_collision_is_an_error(manual, tmp_path):
    with (manual / 'refs.bib').open('a') as file:
        file.write('\n@misc{otherkey, note={A different entry}}\n')
    result = subprocess.run(
        [sys.executable, '-m', 'sphinx', '-b', 'html', str(manual), str(tmp_path / 'html')],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0
    assert 'case-insensitive duplicate keys' in result.stderr


def test_local_definition_and_its_removal(manual, tmp_path):
    (manual / 'chapter.rst').write_text('''
A local definition
==================

.. [ExistingKey] The local source takes precedence.
''')
    output = tmp_path / 'html'
    build(manual, output)
    assert 'href="chapter.html#existingkey"' in (output / 'index.html').read_text()
    (manual / 'chapter.rst').write_text('A chapter\n=========\n\n[ExistingKey]_.\n')
    build(manual, output)
    assert 'href="#existingkey"' in (output / 'index.html').read_text()


def test_cached_citations_from_another_manual(manual, tmp_path):
    # Sage's multidocs cache puts other manuals' citation definitions into
    # the native domain. They must not be treated as docstring definitions
    # when a later edit starts using one of those keys locally.
    with (manual / 'conf.py').open('a') as file:
        file.write('''
def cached_citations(app, env):
    domain = env.get_domain('citation')
    domain.citations.setdefault('OtherKey', ('../references/index', 'otherkey', 0))
    domain.citation_refs.setdefault('OtherKey', set())

def setup(app):
    app.connect('env-updated', cached_citations, priority=700)
''')
    (manual / 'index.rst').write_text('Title\n=====\n\n[ExistingKey]_.\n')
    for page in manual.glob('*.rst'):
        if page.name != 'index.rst':
            page.unlink()
    output = tmp_path / 'html'
    build(manual, output)
    with (manual / 'index.rst').open('a') as file:
        file.write('\nNow also [OtherKey]_.\n')
    build(manual, output)
    page = (output / 'index.html').read_text()
    assert 'href="#otherkey"' in page
    assert 'https://example.org/a_b' in page


@pytest.mark.parametrize('smartquotes', [True, False])
def test_bibtex_typography_matches_native_citations(manual, tmp_path, smartquotes):
    # BibTeX inserts its entries after Sphinx's normal typography transform.
    # Compare actual rendered paragraphs with native reST citations, including
    # quotes crossing inline markup and text which must stay literal.
    with (manual / 'refs.bib').open('a') as file:
        file.write(r'''
@misc{QuoteKey,
  annotation = {A. Author. "\emph{An author's title}".
                \texttt{'literal'} and \(x'\).
                \href{https://example.org/a_b}{An author's website}.
                Escaped \textquotesingle{}quote\textquotesingle{} and
                \textquotedbl{}double\textquotedbl{}.}
}
''')
    with (manual / 'conf.py').open('a') as file:
        file.write(f'\nsmartquotes = {smartquotes!r}\n')
    index = manual / 'index.rst'
    index.write_text(r'''
A manual
========

[NativeQuoteKey]_.

.. [NativeQuoteKey] A. Author. "*An author's title*".
   ``'literal'`` and :math:`x'`.
   `An author's website <https://example.org/a_b>`_.
   Escaped \'quote\' and \"double\".
''')
    for page in manual.glob('*.rst'):
        if page != index:
            page.unlink()
    output = tmp_path / 'html'

    def paragraph(key):
        page = (output / 'index.html').read_text()
        citation = re.search(rf'<div class="citation" id="{key}".*?</div>', page, re.S)[0]
        body = re.search(r'<p>(.*?)</p>', citation, re.S)[1]
        return re.sub(r'\s+', ' ', unescape(body)).strip()

    build(manual, output)
    native = paragraph('nativequotekey')
    index.write_text('A manual\n========\n\n[QuoteKey]_.\n')
    build(manual, output)
    assert paragraph('quotekey') == native


def test_heading_colliding_with_uncited_entry(manual, tmp_path):
    # Even uncited entries reserve IDs when parsing the bibliography. This
    # used to break Sage's GSL heading while building the libraries manual.
    (manual / 'index.rst').write_text('ExistingKey\n===========\n\n[OtherKey]_.\n')
    for page in manual.glob('*.rst'):
        if page.name != 'index.rst':
            page.unlink()
    output = tmp_path / 'html'
    build(manual, output)
    page = (output / 'index.html').read_text()
    assert 'id="existingkey"' in page
    assert 'href="#bibtex-1-otherkey"' in page


def test_root_local_definition_and_removal(manual, tmp_path):
    index = manual / 'index.rst'
    index.write_text('''Title
=====

[ExistingKey]_ and [OtherKey]_.

.. [ExistingKey] A local definition in the bibliography's document.
''')
    for page in manual.glob('*.rst'):
        if page != index:
            page.unlink()
    output = tmp_path / 'html'
    build(manual, output)
    page = (output / 'index.html').read_text()
    assert 'href="#existingkey"' in page
    assert 'href="#bibtex-1-otherkey"' in page
    index.write_text('Title\n=====\n\n[ExistingKey]_ and [OtherKey]_.\n')
    build(manual, output)
    page = (output / 'index.html').read_text()
    assert 'href="#existingkey"' in page
    assert 'href="#otherkey"' in page
    assert 'A title with preserved capitals' in page


def test_bibliography_edit_after_inventory_build(manual, tmp_path):
    # Meson's inventory and HTML targets reuse the same doctrees. Reading the
    # changed .bib during inventory must not leave the HTML output stale.
    with (manual / 'conf.py').open('a') as file:
        file.write("\nextensions.append('sage_docbuild.ext.inventory_builder')\n")
    doctrees = tmp_path / 'doctrees'
    output = tmp_path / 'html'
    build(manual, output, doctrees=doctrees)
    bibfile = manual / 'refs.bib'
    bibfile.write_text(bibfile.read_text().replace('A title with preserved capitals', 'An updated title'))
    build(manual, tmp_path / 'inventory', 'inventory', doctrees=doctrees)
    build(manual, output, doctrees=doctrees)
    assert 'An updated title' in (output / 'index.html').read_text()
    assert 'href="index.html#existingkey"' in (output / 'chapter.html').read_text()
    previous = (output / 'index.html').stat().st_mtime_ns
    build(manual, output, doctrees=doctrees)
    assert (output / 'index.html').stat().st_mtime_ns == previous


def test_citation_scanner_loads_all_bibliography_files(tmp_path):
    from sage_docbuild.sphinxbuild import _documentation_citation_labels

    directory = tmp_path / 'references'
    directory.mkdir()
    (directory / 'first.bib').write_text('@misc{FirstKey, note={First reference}}')
    (directory / 'second.bib').write_text('@misc{SecondKey, note={Second reference}}')
    (directory / 'unrelated.txt').write_text('Not a bibliography.')
    labels, have_bibliography = _documentation_citation_labels([tmp_path])
    assert have_bibliography
    assert labels == {'FirstKey', 'SecondKey'}
