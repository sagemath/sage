# sage.doctest: needs sphinx
"""
Documentation builders

.. NOTE::

   If you are a developer and want to build the SageMath documentation from source,
   refer to `developer's guide <../../../developer/sage_manuals.html>`_.

This module is the starting point for building documentation, and is
responsible to figure out what to build and with which options. The actual
documentation build for each individual document is then done in a subprocess
call to Sphinx, see :func:`builder_helper`. Note that

* The builders are configured with ``build_options.py``;
* The Sphinx subprocesses are configured in ``conf.py``.

:class:`DocBuilder` is the base class of all Builders. It has builder helpers
:meth:`html()`, :meth:`latex`, :meth:`pdf`, :meth:`inventory`, etc, which are
invoked depending on the output type. Each type corresponds with the Sphinx
builder format, except that ``pdf`` is Sphinx latex builder plus compiling
latex to pdf. Note that Sphinx inventory builder is not native to Sphinx
but provided by Sage. See ``sage_docbuild/ext/inventory_builder.py``. The
Sphinx inventory builder is a dummy builder with no actual output but produces
doctree files in ``local/share/doctree`` and ``inventory.inv`` inventory files
in ``local/share/inventory``.

The reference manual is built in two passes, first by :class:`ReferenceBuilder`
with ``inventory`` output type and secondly with ``html`` output type. The
:class:`ReferenceBuilder` itself uses :class:`ReferenceTopBuilder` and
:class:`ReferenceSubBuilder` to build subcomponents of the reference manual.
The :class:`ReferenceSubBuilder` examines the modules included in the
subcomponent by comparing the modification times of the module files with the
times saved in ``local/share/doctree/reference.pickle`` from the previous
build. Then new rst files are generated for new and updated modules. See
:meth:`get_new_and_updated_modules()`.

After :issue:`31948`, when Sage is built, :class:`ReferenceBuilder` is not used
and its responsibility is now taken by the ``Makefile`` in ``$SAGE_ROOT/src/doc``.
"""

# ****************************************************************************
#       Copyright (C) 2008-2009 Mike Hansen <mhansen@gmail.com>
#                     2009-2010 Mitesh Patel <qed777@gmail.com>
#                     2009-2015 J. H. Palmieri <palmieri@math.washington.edu>
#                     2009 Carl Witty <cwitty@newtonlabs.com>
#                     2010-2017 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#                     2012 William Stein <wstein@gmail.com>
#                     2012-2014 Nicolas M. Thiery <nthiery@users.sf.net>
#                     2012-2015 André Apitzsch <andre.apitzsch@etit.tu-chemnitz.de>
#                     2012 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#                     2013-2014 Volker Braun <vbraun.name@gmail.com>
#                     2013 R. Andrew Ohana <andrew.ohana@gmail.com>
#                     2015 Thierry Monteil <sage@lma.metelu.net>
#                     2015 Marc Mezzarobba <marc@mezzarobba.net>
#                     2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#                     2016-2017 Frédéric Chapoton <chapoton@math.univ-lyon1.fr>
#                     2016 Erik M. Bray <erik.bray@lri.fr>
#                     2017 Kwankyu Lee <ekwankyu@gmail.com>
#                     2017 François Bissey <frp.bissey@gmail.com>
#                     2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import logging
import os
import pickle
import re
import shlex
import shutil
import subprocess
import sys
import time
import warnings
from collections.abc import Generator
from pathlib import Path
from typing import Literal

from . import build_options
from .build_options import BuildOptions
from .utils import build_many as _build_many

logger = logging.getLogger(__name__)


##########################################
#      Parallel Building Ref Manual      #
##########################################

def build_ref_doc(args):
    doc = args[0]
    lang = args[1]
    format = args[2]
    kwds = args[3]
    args = args[4:]
    if format == 'inventory':  # you must not use the inventory to build the inventory
        kwds['use_multidoc_inventory'] = False
    getattr(ReferenceSubBuilder(doc, lang), format)(*args, **kwds)


##########################################
#             Builders                   #
##########################################

def builder_helper(type):
    """
    Return a function which builds the documentation for
    output type ``type``.
    """
    def f(self, *args, **kwds):
        output_dir = self._output_dir(type)

        options = build_options.ALLSPHINXOPTS

        if self.name == 'website':
            # WEBSITESPHINXOPTS is either empty or " -A hide_pdf_links=1 "
            options += build_options.WEBSITESPHINXOPTS

        if kwds.get('use_multidoc_inventory', True) and type != 'inventory':
            options += ' -D multidoc_first_pass=0'
        else:
            options += ' -D multidoc_first_pass=1'

        build_command = '-b %s -d %s %s %s %s' % (type, self._doctrees_dir(),
                                                  options, self.dir,
                                                  output_dir)

        # Provide "pdf" tag to be used with "only" directive as an alias of "latex"
        if type == 'latex':
            build_command = '-t pdf ' + build_command

        logger.debug(build_command)

        # Run Sphinx with Sage's special logger
        sys.argv = ["sphinx-build"] + build_command.split()
        from .sphinxbuild import runsphinx
        try:
            runsphinx()
        except Exception:
            if build_options.ABORT_ON_ERROR:
                raise
        except BaseException as e:
            # We need to wrap a BaseException that is not an Exception in a
            # regular Exception. Otherwise multiprocessing.Pool.get hangs, see
            # #25161
            if build_options.ABORT_ON_ERROR:
                raise Exception("Non-exception during docbuild: %s" % (e,), e)

        if type == 'latex':
            logger.warning(f"LaTeX files can be found in {output_dir}.")
        elif type != 'inventory':
            logger.warning(f"Build finished. The built documents can be found in {output_dir}.")

    f.is_output_format = True
    return f


class DocBuilder():
    def __init__(self, name: str, options: BuildOptions):
        """
        INPUT:

        - ``name`` -- the name of a subdirectory in ``doc/<lang>``, such as
          'tutorial' or 'installation'
        """
        self.name = name
        self.dir = options.source_dir / self.name
        self._options = options

    def _output_dir(self, type):
        """
        Return the directory where the output of type ``type`` is stored.

        If the directory does not exist, then it will automatically be
        created.

        EXAMPLES::

            sage: from sage_docbuild.builders import DocBuilder
            sage: from sage_docbuild.build_options import BuildOptions
            sage: import tempfile
            sage: with tempfile.TemporaryDirectory() as directory:
            ....:   options = BuildOptions(output_dir=Path(directory), source_dir=Path('src/doc'))
            ....:   builder = DocBuilder('en/tutorial', options)
            ....:   builder._output_dir('html')
            ...Path('.../html/en/tutorial')
        """
        dir = self._options.output_dir / type / self.name
        dir.mkdir(parents=True, exist_ok=True)
        return dir

    def _doctrees_dir(self) -> Path:
        """
        Return the directory where the doctrees are stored.

        If the directory does not exist, then it will automatically be
        created.

        EXAMPLES::

            sage: from sage_docbuild.builders import DocBuilder
            sage: from sage_docbuild.build_options import BuildOptions
            sage: import tempfile
            sage: with tempfile.TemporaryDirectory() as directory:
            ....:   options = BuildOptions(output_dir=Path(directory), source_dir=Path('src/doc'))
            ....:   builder = DocBuilder('en/tutorial', options)
            ....:   builder._doctrees_dir()
            ...Path('.../doctrees/en/tutorial')
        """
        dir = self._options.output_dir / 'doctrees' / self.name
        dir.mkdir(parents=True, exist_ok=True)
        return dir

    def _output_formats(self):
        """
        Return a list of the possible output formats.

        EXAMPLES::

            sage: from sage_docbuild.builders import DocBuilder
            sage: from sage_docbuild.build_options import BuildOptions
            sage: options = BuildOptions(source_dir=Path('src/doc'))
            sage: builder = DocBuilder('tutorial', options)
            sage: builder._output_formats()
            ['changes', 'html', 'htmlhelp', 'inventory', 'json', 'latex', 'linkcheck', 'pickle', 'web']
        """
        # Go through all the attributes of self and check to
        # see which ones have an 'is_output_format' attribute.  These
        # are the ones created with builder_helper.
        output_formats = []
        for attr in dir(self):
            if hasattr(getattr(self, attr), 'is_output_format'):
                output_formats.append(attr)
        output_formats.sort()
        return output_formats

    def pdf(self):
        """
        Build the PDF files for this document.

        This is done by first (re)-building the LaTeX output, going
        into that LaTeX directory, and running 'make all-pdf' there.

        EXAMPLES::

            sage: from sage_docbuild.builders import DocBuilder
            sage: from sage_docbuild.build_options import BuildOptions
            sage: options = BuildOptions(source_dir = Path('src/doc'))
            sage: builder = DocBuilder('tutorial', options)
            sage: builder.pdf() #not tested
        """
        self.latex()
        tex_dir = self._output_dir('latex')
        pdf_dir = self._output_dir('pdf')

        if self.name == 'reference':
            # recover maths in tex, undoing what Sphinx did (trac #29993)
            tex_file = tex_dir / 'reference.tex'
            with open(tex_file) as f:
                ref = f.read()
                ref = re.sub(r'\\textbackslash{}', r'\\', ref)
                ref = re.sub(r'\\textbackslash{}', r'\\', ref)
                ref = re.sub(r'\\{', r'{', ref)
                ref = re.sub(r'\\}', r'}', ref)
                ref = re.sub(r'\\_', r'_', ref)
                ref = re.sub(r'\\textasciicircum{}', r'^', ref)
            with open(tex_file, 'w') as f:
                f.write(ref)

        make_cmd = os.environ.get('MAKE', 'make')
        command = shlex.split(make_cmd) + ['all-pdf']
        logger.debug(f"Running {' '.join(command)} in {tex_dir}")

        proc = subprocess.run(
            command,
            check=False, cwd=tex_dir,
            capture_output=True,
            text=True,
        )

        if proc.returncode != 0:
            logger.error(f"stdout from {make_cmd}:\n{proc.stdout}")
            logger.error(f"stderr from {make_cmd}:\n{proc.stderr}")
            raise RuntimeError(f"failed to run {' '.join(command)} in {tex_dir}")

        if proc.stdout:
            logger.debug(f"make stdout:\n{proc.stdout}")
        if proc.stderr:
            # Still surface stderr even on success, but at debug level
            logger.debug(f"make stderr:\n{proc.stderr}")

        # Move generated PDFs
        for pdf in tex_dir.glob("*.pdf"):
            try:
                dst_pdf = os.path.join(pdf_dir, os.path.basename(pdf))
                shutil.move(str(pdf), dst_pdf)
            except Exception as e:
                logger.error(f"Failed moving {pdf} to {dst_pdf}: {e}")
                raise

        logger.info(f"Build finished. The built documents can be found in {pdf_dir}.")

    def clean(self, *args):
        shutil.rmtree(self._doctrees_dir())
        output_formats = list(args) if args else self._output_formats()
        for format in output_formats:
            shutil.rmtree(self._output_dir(format), ignore_errors=True)

    html = builder_helper('html')
    pickle = builder_helper('pickle')
    web = pickle
    json = builder_helper('json')
    htmlhelp = builder_helper('htmlhelp')
    latex = builder_helper('latex')
    changes = builder_helper('changes')
    linkcheck = builder_helper('linkcheck')
    # import the customized builder for object.inv files
    inventory = builder_helper('inventory')


def build_many(target, args, processes=None):
    """
    Thin wrapper around :func:`sage_docbuild.utils.build_many` which uses the
    docbuild settings ``NUM_THREADS`` and ``ABORT_ON_ERROR``.
    """
    if processes is None:
        processes = build_options.NUM_THREADS
    try:
        _build_many(target, args, processes=processes)
    except BaseException:
        if build_options.ABORT_ON_ERROR:
            raise


##########################################
#      Parallel Building Ref Manual      #
##########################################
class WebsiteBuilder(DocBuilder):
    def html(self):
        """
        After we have finished building the website index page, we copy
        everything one directory up, that is, to the base diectory ``html/en``.

        In addition, an index file is installed into the root doc directory.

        Thus we have three index.html files:

            html/en/website/index.html  (not used)
            html/en/index.html  (base directory)
            index.html  (root doc directory)
        """
        super().html()
        html_output_dir = self._output_dir('html')

        # This file is used by src/doc/common/static/jupyter-sphinx-furo.js
        # for doc version selector
        shutil.copy2(os.path.join(self.dir, 'versions.txt'), html_output_dir)

        for f in os.listdir(html_output_dir):
            src = os.path.join(html_output_dir, f)
            dst = os.path.join(html_output_dir, '..', f)
            if os.path.isdir(src):
                shutil.rmtree(dst, ignore_errors=True)
                shutil.copytree(src, dst)
            else:
                shutil.copy2(src, dst)

        shutil.copy2(os.path.join(self.dir, 'root_index.html'),
                     os.path.join(html_output_dir, '../../../index.html'))

    def pdf(self):
        """
        Build the website hosting pdf docs.
        """
        super().pdf()

        # If the website exists, update it.

        from sage.env import SAGE_DOC
        website_dir = os.path.join(SAGE_DOC, 'html', 'en', 'website')

        if os.path.exists(os.path.join(website_dir, 'index.html')):
            # Rebuild WITHOUT --no-pdf-links, which is translated to
            # "-A hide_pdf_links=1" Sphinx argument. Thus effectively
            # the index page SHOWS links to pdf docs.
            self.html()

    def clean(self):
        """
        When we clean the output for the website index, we need to
        remove all of the HTML that were placed in the parent
        directory.

        In addition, remove the index file installed into the root doc directory.
        """
        html_output_dir = self._output_dir('html')
        parent_dir = os.path.realpath(os.path.join(html_output_dir, '..'))
        for filename in os.listdir(html_output_dir):
            parent_filename = os.path.join(parent_dir, filename)
            if not os.path.exists(parent_filename):
                continue
            if os.path.isdir(parent_filename):
                shutil.rmtree(parent_filename, ignore_errors=True)
            else:
                os.unlink(parent_filename)

        root_index_file = os.path.join(html_output_dir, '../../../index.html')
        if os.path.exists(root_index_file):
            os.remove(root_index_file)

        DocBuilder.clean(self)


class ReferenceBuilder():
    """
    This class builds the reference manual. It uses DocBuilder to
    build the top-level page and ReferenceSubBuilder for each
    sub-component.
    """
    def __init__(self, name:str, options: BuildOptions):
        """
        Record the reference manual's name, in case it's not
        identical to 'reference'.
        """
        self.name = name
        self.options = options

    def _output_dir(self, type: Literal['html', 'latex', 'pdf']) -> Path:
        """
        Return the directory where the output of type ``type`` is stored.

        If the directory does not exist, then it will automatically be
        created.

        EXAMPLES::

            sage: from sage_docbuild.builders import ReferenceBuilder
            sage: from sage_docbuild.build_options import BuildOptions
            sage: import tempfile
            sage: with tempfile.TemporaryDirectory() as directory:
            ....:   options = BuildOptions(output_dir = Path(directory))
            ....:   builder = ReferenceBuilder('reference', options)
            ....:   builder._output_dir('html')
            ...Path('.../html/reference')
        """
        dir = self.options.output_dir / type / self.name
        dir.mkdir(parents=True, exist_ok=True)
        return dir

    def _source_dir(self) -> Path:
        return self.options.source_dir / self.name

    def _build_bibliography(self, format, *args, **kwds):
        """
        Build the bibliography only

        The bibliography references.aux is referenced by the other
        manuals and needs to be built first.
        """
        references = [
            (doc, 'en', format, kwds) + args for doc in get_all_documents(self._source_dir())
            if doc == 'reference/references'
        ]
        build_many(build_ref_doc, references)

    def _build_everything_except_bibliography(self, format, *args, **kwds):
        """
        Build the entire reference manual except the bibliography
        """
        non_references = [
            (doc, 'en', format, kwds) + args for doc in get_all_documents(self._source_dir())
            if doc != Path('reference/references')
        ]
        build_many(build_ref_doc, non_references)

    def _build_top_level(self, format, *args, **kwds):
        """
        Build top-level document.
        """
        getattr(ReferenceTopBuilder('reference', self.options), format)(*args, **kwds)

    def _wrapper(self, format, *args, **kwds):
        """
        Build reference manuals: build the top-level document and its components.
        """
        logger.info('Building bibliography')
        self._build_bibliography(format, *args, **kwds)
        logger.info('Bibliography finished, building dependent manuals')
        self._build_everything_except_bibliography(format, *args, **kwds)
        # The html refman must be built at the end to ensure correct
        # merging of indexes and inventories.
        # Sphinx is run here in the current process (not in a
        # subprocess) and the IntersphinxCache gets populated to be
        # used for the second pass of the reference manual and for
        # the other documents.
        self._build_top_level(format, *args, **kwds)

class ReferenceTopBuilder(DocBuilder):
    """
    This class builds the top-level page of the reference manual.
    """
    def __init__(self, name: str, options: BuildOptions):
        DocBuilder.__init__(self, 'en/reference', options)

    def html(self):
        """
        Build the top-level document.
        """
        super().html()

        # We want to build master index file which lists all of the PDF file.
        # We modify the file index.html from the "reference_top" target, if it
        # exists. Otherwise, we are done.
        output_dir = self._output_dir('html')

        # Install in output_dir a symlink to the directory containing static files.
        # Prefer relative path for symlinks.
        relpath = output_dir.relative_to(self._options.output_dir)
        try:
            (output_dir / '_static').symlink_to(relpath / '_static')
        except FileExistsError:
            pass

        # Now modify top reference index.html page and write it to output_dir.
        with open(output_dir / 'index.html') as f:
            html = f.read()
        # Fix links in navigation bar
        html = re.sub(r'<a href="(.*)">Sage(.*)Documentation</a>',
                      r'<a href="../../../html/en/index.html">Sage\2Documentation</a>',
                      html)
        html = re.sub(r'<li class="right"(.*)>', r'<li class="right" style="display: none" \1>',
                      html)
        html = re.sub(r'<div class="sphinxsidebar"(.*)>', r'<div class="sphinxsidebar" style="display: none" \1>',
                      html)

        # From index.html, we want the preamble and the tail.
        html_end_preamble = html.find(r'<section')
        html_bottom = html.rfind(r'</section>') + len(r'</section>')

        # For the content, we modify doc/en/reference/index.rst, which
        # has two parts: the body and the table of contents.
        with open(self.dir / 'index.rst') as f:
            rst = f.read()
        # Get rid of todolist and miscellaneous rst markup.
        rst = rst.replace('.. _reference-manual:\n\n', '')
        rst = re.sub(r'\\\\', r'\\', rst)
        # Replace rst links with html links. There are three forms:
        #
        #   `blah`__    followed by __ LINK
        #
        #   `blah <LINK>`_
        #
        #   :doc:`blah <module/index>`
        #
        # Change the first and the second forms to
        #
        #   <a href="LINK">blah</a>
        #
        # Change the third form to
        #
        #   <a href="module/module.pdf"><img src="_static/pdf.png">blah</a>
        #
        rst = re.sub(r'`([^`\n]*)`__.*\n\n__ (.*)',
                     r'<a href="\2">\1</a>.', rst)
        rst = re.sub(r'`([^<\n]*)\s+<(.*)>`_',
                     r'<a href="\2">\1</a>', rst)
        rst = re.sub(r':doc:`([^<]*?)\s+<(.*)/index>`',
                     r'<a title="PDF" class="pdf" href="../../../pdf/en/reference/\2/\2.pdf"><img src="_static/pdf.png"></a><a href="\2/index.html">\1</a> ', rst)
        # Body: add paragraph <p> markup.
        start = rst.rfind('*\n') + 1
        end = rst.find('\nUser Interfaces')
        rst_body = rst[start:end]
        rst_body = rst_body.replace('\n\n', '</p>\n<p>')
        # TOC: don't include the indices
        start = rst.find('\nUser Interfaces')
        end = rst.find('Indices and Tables')
        rst_toc = rst[start:end]
        # change * to <li>; change rst headers to html headers
        rst_toc = re.sub(r'\*(.*)\n',
                         r'<li>\1</li>\n', rst_toc)
        rst_toc = re.sub(r'\n([A-Z][a-zA-Z, ]*)\n[=]*\n',
                         r'</ul>\n\n\n<h2>\1</h2>\n\n<ul>\n', rst_toc)
        rst_toc = re.sub(r'\n([A-Z][a-zA-Z, ]*)\n[-]*\n',
                         r'</ul>\n\n\n<h3>\1</h3>\n\n<ul>\n', rst_toc)
        # now write the file.
        with open(output_dir / 'index-pdf.html', 'w') as new_index:
            new_index.write(html[:html_end_preamble])
            new_index.write('<h1>Sage Reference Manual</h1>')
            new_index.write(rst_body)
            new_index.write('<ul>')
            new_index.write(rst_toc)
            new_index.write('</ul>\n\n')
            new_index.write(html[html_bottom:])


class ReferenceSubBuilder(DocBuilder):
    """
    This class builds sub-components of the reference manual. It is
    responsible for making sure that the auto generated reST files for the
    Sage library are up to date.

    When building any output, we must first go through and check
    to see if we need to update any of the autogenerated reST
    files. There are two cases where this would happen:

    1. A new module gets added to one of the toctrees.
    2. The actual module gets updated and possibly contains a new title.
    """
    _cache = None

    def __init__(self, name: str, options: BuildOptions):
        DocBuilder.__init__(self, "en/" + name, options)
        self._wrap_builder_helpers()

    def _wrap_builder_helpers(self):
        from functools import partial, update_wrapper
        for attr in dir(self):
            if hasattr(getattr(self, attr), 'is_output_format'):
                f = partial(self._wrapper, attr)
                f.is_output_format = True
                update_wrapper(f, getattr(self, attr))
                setattr(self, attr, f)

    def _wrapper(self, build_type, *args, **kwds):
        """
        This is the wrapper around the builder_helper methods that
        goes through and makes sure things are up to date.
        """
        # Force regeneration of all modules if the inherited
        # and/or underscored members options have changed.
        cache = self.get_cache()
        force = False
        try:
            if (cache['option_inherited'] != self._options.inherited or
                    cache['option_underscore'] != self._options.underscore):
                logger.info("Detected change(s) in inherited and/or underscored members option(s).")
                force = True
        except KeyError:
            force = True
        cache['option_inherited'] = self._options.inherited
        cache['option_underscore'] = self._options.underscore
        self.save_cache()

        # Refresh the reST file mtimes in environment.pickle
        if self._options.update_mtimes:
            logger.info("Checking for reST file mtimes to update...")
            self.update_mtimes()

        if force:
            # Write reST files for all modules from scratch.
            self.clean_auto()
            for module_name in self.get_all_included_modules():
                self.write_auto_rest_file(module_name)
        else:
            # Write reST files for new and updated modules.
            for module_name in self.get_new_and_updated_modules():
                self.write_auto_rest_file(module_name)

        # Copy over the custom reST files from _sage
        _sage = self.dir / '_sage'
        if _sage.exists():
            logger.info(f"Copying over custom reST files from {_sage} ...")
            shutil.copytree(_sage, self.dir / 'sage')

        getattr(DocBuilder, build_type)(self, *args, **kwds)

    def cache_file(self) -> Path:
        """
        Return the filename where the pickle of the reference cache
        is stored.
        """
        return self._doctrees_dir() / 'reference.pickle'

    def get_cache(self):
        """
        Retrieve the reference cache which contains the options previously used
        by the reference builder.

        If it doesn't exist, then we just return an empty dictionary. If it
        is corrupted, return an empty dictionary.
        """
        if self._cache is not None:
            return self._cache

        cache_file = self.cache_file()
        if not cache_file.exists():
            return {}
        try:
            with cache_file.open('rb') as file:
                cache = pickle.load(file)
        except Exception:
            logger.debug(f"Cache file '{cache_file}' is corrupted; ignoring it...")
            cache = {}
        else:
            logger.debug(f"Loaded the reference cache: {cache_file}")
        self._cache = cache
        return cache

    def save_cache(self):
        """
        Pickle the current reference cache for later retrieval.
        """
        cache = self.get_cache()
        try:
            with open(self.cache_file(), 'wb') as file:
                pickle.dump(cache, file)
            logger.debug("Saved the reference cache: %s", self.cache_file())
        except PermissionError:
            logger.debug("Permission denied for the reference cache: %s", self.cache_file())

    def get_sphinx_environment(self):
        """
        Return the Sphinx environment for this project.
        """
        env_pickle = os.path.join(self._doctrees_dir(), 'environment.pickle')
        try:
            with open(env_pickle, 'rb') as f:
                env = pickle.load(f)
                logger.debug("Opened Sphinx environment: %s", env_pickle)
                return env
        except (OSError, EOFError):
            logger.debug(
                f"Failed to open Sphinx environment '{env_pickle}'", exc_info=True)

    def update_mtimes(self):
        """
        Update the modification times for reST files in the Sphinx
        environment for this project.
        """
        env = self.get_sphinx_environment()
        if env is not None:
            for doc in env.all_docs:
                env.all_docs[doc] = time.time()
            logger.info("Updated %d reST file mtimes", len(env.all_docs))

            # This is the only place we need to save (as opposed to
            # load) Sphinx's pickle, so we do it right here.
            env_pickle = os.path.join(self._doctrees_dir(), 'environment.pickle')

            # remove unpicklable attributes
            env.set_warnfunc(None)
            with open(env_pickle, 'wb') as picklefile:
                pickle.dump(env, picklefile, pickle.HIGHEST_PROTOCOL)

            logger.debug("Saved Sphinx environment: %s", env_pickle)

    def get_modified_modules(self):
        """
        Return an iterator for all the modules that have been modified
        since the documentation was last built.
        """
        env = self.get_sphinx_environment()
        if env is None:
            logger.debug("Stopped check for modified modules.")
            return
        try:
            added, changed, removed = env.get_outdated_files(False)
            logger.info("Sphinx found %d modified modules", len(changed))
        except OSError as err:
            logger.debug("Sphinx failed to determine modified modules: %s", err)
            return
        for name in changed:
            # Only pay attention to files in a directory sage/... In
            # particular, don't treat a file like 'sagetex.rst' in
            # doc/en/reference/misc as an autogenerated file: see
            # #14199.
            if name.startswith('sage' + os.sep):
                yield name

    def print_modified_modules(self):
        """
        Print a list of all the modules that have been modified since
        the documentation was last built.
        """
        for module_name in self.get_modified_modules():
            print(module_name)

    def get_all_rst_files(self) -> Generator[Path, None, None]:
        """
        Return an iterator for all rst files which are not autogenerated.
        """
        for file in self.dir.rglob('*.rst'):
            if 'sage' in file.relative_to(self.dir).parts:
                continue
            yield file

    def get_all_included_modules(self):
        """
        Return an iterator for all modules which are included in the
        reference manual.
        """
        for file in self.get_all_rst_files():
            for module in self.get_modules(file):
                yield module

    def get_new_and_updated_modules(self):
        """
        Return an iterator for all new and updated modules that appear in
        the toctrees, and remove obsolete old modules.
        """
        env = self.get_sphinx_environment()
        if env is None:
            all_docs = {}
        else:
            all_docs = env.all_docs

        new_modules = []
        updated_modules = []
        old_modules = []
        for module_name in self.get_all_included_modules():
            docname = module_name.replace('.', os.path.sep)

            if docname not in all_docs:
                new_modules.append(module_name)
                yield module_name
                continue

            # get the modification timestamp of the reST doc for the module
            mtime = all_docs[docname]
            try:
                with warnings.catch_warnings():
                    # primarily intended to ignore deprecation warnings
                    warnings.simplefilter("ignore")
                    __import__(module_name)
            except ImportError as err:
                logger.error("Warning: Could not import %s %s", module_name, err)
                raise

            module_filename = sys.modules[module_name].__file__
            if module_filename is None:
                # Namespace package
                old_modules.append(module_name)
                continue
            if (module_filename.endswith('.pyc') or module_filename.endswith('.pyo')):
                source_filename = module_filename[:-1]
                if (os.path.exists(source_filename)):
                    module_filename = source_filename
            newtime = os.path.getmtime(module_filename)

            if newtime > mtime:
                updated_modules.append(module_name)
                yield module_name
            else:  # keep good old module
                old_modules.append(module_name)

        removed_modules = []
        for docname in all_docs.keys():
            if docname.startswith('sage' + os.path.sep):
                module_name = docname.replace(os.path.sep, '.')
                if not (module_name in old_modules or module_name in updated_modules):
                    try:
                        os.remove(os.path.join(self.dir, docname) + '.rst')
                    except OSError:  # already removed
                        pass
                    logger.debug("Deleted auto-generated reST file {}".format(docname))
                    removed_modules.append(module_name)

        logger.info("Found %d new modules", len(new_modules))
        logger.info("Found %d updated modules", len(updated_modules))
        logger.info("Removed %d obsolete modules", len(removed_modules))

    def print_new_and_updated_modules(self):
        """
        Print all the modules that appear in the toctrees that
        are newly included or updated.
        """
        for module_name in self.get_new_and_updated_modules():
            print(module_name)

    def get_modules(self, file: Path) -> Generator[str, None, None]:
        """
        Given a reST file, return an iterator for
        all of the autogenerated reST files that it includes.
        """
        # Create the regular expression used to detect an autogenerated file
        auto_re = re.compile(r'^\s*(..\/)*(sage(_docbuild)?\/[\w\/]*)\s*$')

        # Read the lines
        with file.open(encoding='utf-8') as f:
            lines = f.readlines()

        for line in lines:
            match = auto_re.match(line)
            if match:
                yield match.group(2).replace('/', '.')

    def get_module_docstring_title(self, module_name):
        """
        Return the title of the module from its docstring.
        """
        # Try to import the module
        try:
            __import__(module_name)
        except ImportError as err:
            logger.error("Warning: Could not import %s %s", module_name, err)
            return "UNABLE TO IMPORT MODULE"
        module = sys.modules[module_name]

        # Get the docstring
        doc = module.__doc__
        if doc is None:
            doc = module.doc if hasattr(module, 'doc') else ""

        # Extract the title
        i = doc.find('\n')
        if i != -1:
            return doc[i + 1:].lstrip().splitlines()[0]
        else:
            return doc

    def auto_rest_filename(self, module_name: str) -> Path:
        """
        Return the name of the file associated to a given module

        EXAMPLES::

            sage: from sage_docbuild.builders import ReferenceSubBuilder
            sage: from sage_docbuild.build_options import BuildOptions
            sage: options = BuildOptions(source_dir = Path('src/doc'))
            sage: ReferenceSubBuilder("reference", options).auto_rest_filename("sage.combinat.partition")
            ...Path('src/doc/en/reference/sage/combinat/partition.rst')
        """
        return self.dir / (module_name.replace('.', os.path.sep) + '.rst')

    def write_auto_rest_file(self, module_name: str):
        """
        Write the autogenerated reST file for module_name.
        """
        if not module_name.startswith('sage'):
            return

        title = self.get_module_docstring_title(module_name)
        if title == '':
            logger.error("Warning: Missing title for %s", module_name)
            title = "MISSING TITLE"

        rst_file = self.auto_rest_filename(module_name)
        rst_file.parent.mkdir(parents=True, exist_ok=True)
        with rst_file.open('w') as outfile:
            # Don't doctest the autogenerated file.
            outfile.write(".. nodoctest\n\n")
            # Now write the actual content.
            outfile.write(".. _%s:\n\n" % (module_name.replace(".__init__", "")))
            outfile.write(title + '\n')
            outfile.write('=' * len(title) + "\n\n")
            outfile.write('.. This file has been autogenerated.\n\n')

            inherited = ':inherited-members:' if self._options.inherited else ''

            automodule = '''
.. automodule:: %s
   :members:
   :undoc-members:
   :show-inheritance:
   %s

'''
            outfile.write(automodule % (module_name, inherited))

    def clean_auto(self):
        """
        Remove all autogenerated reST files.
        """
        try:
            shutil.rmtree(os.path.join(self.dir, 'sage'))
            logger.debug("Deleted auto-generated reST files in: %s",
                         os.path.join(self.dir, 'sage'))
        except OSError:
            pass

    def get_unincluded_modules(self):
        """
        Return an iterator for all the modules in the Sage library
        which are not included in the reference manual.
        """
        # Make a dictionary of the included modules
        included_modules = {}
        for module_name in self.get_all_included_modules():
            included_modules[module_name] = True

        base_path = os.path.join(SAGE_SRC, 'sage')
        for directory, subdirs, files in os.walk(base_path):
            for filename in files:
                if not (filename.endswith('.py') or
                        filename.endswith('.pyx')):
                    continue

                path = os.path.join(directory, filename)

                # Create the module name
                module_name = path[len(base_path):].replace(os.path.sep, '.')
                module_name = 'sage' + module_name
                module_name = module_name[:-4] if module_name.endswith('pyx') else module_name[:-3]

                # Exclude some ones  -- we don't want init the manual
                if module_name.endswith('__init__') or module_name.endswith('all'):
                    continue

                if module_name not in included_modules:
                    yield module_name

    def print_unincluded_modules(self):
        """
        Print all of the modules which are not included in the Sage
        reference manual.
        """
        for module_name in self.get_unincluded_modules():
            print(module_name)

    def print_included_modules(self):
        """
        Print all of the modules that are included in the Sage reference
        manual.
        """
        for module_name in self.get_all_included_modules():
            print(module_name)


class SingleFileBuilder(DocBuilder):
    """
    This is the class used to build the documentation for a single
    user-specified file. If the file is called 'foo.py', then the
    documentation is built in ``DIR/foo/`` if the user passes the
    command line option "-o DIR", or in ``DOT_SAGE/docbuild/foo/``
    otherwise.
    """
    def __init__(self, path):
        """
        INPUT:

        - ``path`` -- the path to the file for which documentation
          should be built
        """
        self.lang = 'en'
        self.name = 'single_file'
        path = os.path.abspath(path)

        # Create docbuild and relevant subdirectories, e.g.,
        # the static and templates directories in the output directory.
        # By default, this is DOT_SAGE/docbuild/MODULE_NAME, but can
        # also be specified at the command line.
        module_name = os.path.splitext(os.path.basename(path))[0]
        latex_name = module_name.replace('_', r'\\_')

        if self._options.output_dir:
            base_dir = os.path.join(self._options.output_dir, module_name)
            if os.path.exists(base_dir):
                logger.warning('Warning: Directory %s exists. It is safer to build in a new directory.' % base_dir)
        else:
            base_dir = os.path.join(DOT_SAGE, 'docbuild', module_name)
            try:
                shutil.rmtree(base_dir)
            except OSError:
                pass
        self.dir = os.path.join(base_dir, 'source')

        os.makedirs(os.path.join(self.dir, "static"), exist_ok=True)
        os.makedirs(os.path.join(self.dir, "templates"), exist_ok=True)
        # Write self.dir/conf.py
        conf = r"""# This file is automatically generated by {}, do not edit!

import sys, os, contextlib
sys.path.append({!r})

from sage.docs.conf import *
html_static_path = [] + html_common_static_path

project = 'Documentation for {}'
release = 'unknown'
name = {!r}
html_title = project
html_short_title = project
htmlhelp_basename = name

with contextlib.suppress(ValueError):
    extensions.remove('multidocs') # see #29651
    extensions.remove('inventory_builder')

latex_domain_indices = False
latex_documents = [
  ('index', name + '.tex', 'Documentation for {}',
   'unknown', 'manual'),
]
""".format(__file__, self.dir, module_name, module_name, latex_name)

        if 'SAGE_DOC_UNDERSCORE' in os.environ:
            conf += r"""
def setup(app):
    app.connect('autodoc-skip-member', skip_member)
"""

        with open(os.path.join(self.dir, 'conf.py'), 'w') as conffile:
            conffile.write(conf)

        # Write self.dir/index.rst
        title = 'Docs for file %s' % path
        heading = title + "\n" + ("=" * len(title))
        index = r"""{}

.. This file is automatically generated by {}, do not edit!

.. automodule:: {}
   :members:
   :undoc-members:
   :show-inheritance:
""".format(heading, __file__, module_name)
        with open(os.path.join(self.dir, 'index.rst'), 'w') as indexfile:
            indexfile.write(index)

        # Create link from original file to self.dir. Note that we
        # append self.dir to sys.path in conf.py. This is reasonably
        # safe (but not perfect), since we just created self.dir.
        try:
            os.symlink(path, os.path.join(self.dir, os.path.basename(path)))
        except OSError:
            pass

    def _output_dir(self, type):
        """
        Return the directory where the output of type ``type`` is stored.

        If the directory does not exist, then it will automatically be
        created.
        """
        base_dir = os.path.split(self.dir)[0]
        d = os.path.join(base_dir, "output", type)
        os.makedirs(d, exist_ok=True)
        return d

    def _doctrees_dir(self):
        """
        Return the directory where the doctrees are stored.

        If the directory does not exist, then it will automatically be
        created.
        """
        return self._output_dir('doctrees')


def get_builder(name: str, options: BuildOptions) -> DocBuilder | ReferenceBuilder:
    """
    Return an appropriate *Builder* object for the document ``name``.

    DocBuilder and its subclasses do all the real work in building the
    documentation.
    """
    if name == 'reference_top':
        return ReferenceTopBuilder('reference', options)
    elif name.endswith('reference'):
        return ReferenceBuilder(name, options)
    elif 'reference' in name and (options.source_dir / 'en' / name).exists():
        return ReferenceSubBuilder(name, options)
    elif name.endswith('website'):
        return WebsiteBuilder(name, options)
    elif name.startswith('file='):
        path = name[5:]
        if path.endswith('.sage') or path.endswith('.pyx'):
            raise NotImplementedError('Building documentation for a single file only works for Python files.')
        return SingleFileBuilder(path)
    elif Path(name) in get_all_documents(options.source_dir):
        return DocBuilder(name, options)
    else:
        print("'%s' is not a recognized document. Type 'sage --docbuild -D' for a list" % name)
        print("of documents, or 'sage --docbuild --help' for more help.")
        sys.exit(1)


def get_all_documents(source: Path) -> list[Path]:
    """
    Return a list of all of the documents, relative to the source
    directory.

    A document is a directory within one of the language
    subdirectories of ``doc``.

    EXAMPLES::

        sage: from sage_docbuild.builders import get_all_documents
        sage: from sage.env import SAGE_DOC_SRC
        sage: documents = get_all_documents(Path(SAGE_DOC_SRC))
        sage: Path('en/tutorial') in documents
        True
    """
    documents = []
    for lang in [path for path in source.iterdir() if path.is_dir()]:
        if not re.match('^[a-z][a-z]$', lang.name):
            # Skip non-language directories
            continue
        for document in lang.iterdir():
            if (document.name not in build_options.OMIT
                    and document.is_dir()):
                documents.append(document.relative_to(source))

    # Top-level reference document is build seperately
    if Path('en/reference') in documents:
        documents.remove(Path('en/reference'))

    return documents

def get_all_reference_documents(source: Path) -> list[Path]:
    """
    Return a list of all reference manual documents to build, relative to the
    specified source directory.

    We add a document if it's a subdirectory of the manual's
    directory and contains a file named 'index.rst'.

    The order corresponds to the order in which the documents should be built.

    EXAMPLES::

        sage: from sage_docbuild.builders import get_all_reference_documents
        sage: from sage.env import SAGE_DOC_SRC
        sage: documents = get_all_reference_documents(Path(SAGE_DOC_SRC) / 'en')
        sage: Path('reference/algebras') in documents
        True
    """
    documents: list[tuple[int, Path]] = []

    for directory in (source / 'reference').iterdir():
        if (directory / 'index.rst').exists():
            n = len(list(directory.iterdir()))
            documents.append((-n, directory.relative_to(source)))

    # Sort largest component (most subdirectory entries) first since 
    # they will take the longest to build
    docs = [doc[1] for doc in sorted(documents)]
    # Put the bibliography first, because it needs to be built first:
    docs.remove(Path('reference/references'))
    docs.insert(0, Path('reference/references'))

    # Add the top-level reference document
    docs.append(Path('reference_top'))

    return docs
