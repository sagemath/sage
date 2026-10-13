"""Integration tests for separately built HTML and LaTeX bibliographies."""

import re
import subprocess
import sys
from html import unescape

import pytest

pytest.importorskip('sphinxcontrib.bibtex')


def build(source, output, builder='html', fresh=False):
    command = [sys.executable, '-m', 'sphinx', '-b', builder, '-W', '--keep-going',
               '-j', '2', str(source), str(output)]
    if fresh:
        command.insert(3, '-E')
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.fixture
def manual(tmp_path):
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
    (source / 'conf.py').write_text('''
extensions = ['sage_docbuild.ext.bibtex']
bibtex_bibfiles = ['refs.bib']
bibtex_default_style = 'sage'
bibtex_cite_id = '{key}'
project = 'Bibliography test'
latex_documents = [('index', 'test.tex', 'Bibliography test', 'Sage', 'manual')]
''')
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
        file.write('\n@misc{existingkey, note={A different entry}}\n')
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
