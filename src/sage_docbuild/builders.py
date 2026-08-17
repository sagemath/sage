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
:meth:`~DocBuilder.html`, :meth:`~DocBuilder.latex`,
:meth:`~DocBuilder.pdf`, :meth:`~DocBuilder.inventory`, etc, which are
invoked depending on the output type. Each type corresponds with the Sphinx
builder format, except that :meth:`~DocBuilder.pdf` is Sphinx latex builder
plus compiling latex to pdf. Note that Sphinx inventory builder is not native
to Sphinx but provided by Sage. See
:mod:`sage_docbuild.ext.inventory_builder`. The
Sphinx inventory builder is a dummy builder with no actual output but produces
doctree files in ``$SAGE_DOC/doctrees`` and ``objects.inv`` inventory files
in ``$SAGE_DOC/inventory``.

The reference manual is built in two passes, first by :class:`ReferenceBuilder`
with ``inventory`` output type and secondly with ``html`` output type. The
:class:`ReferenceBuilder` itself uses :class:`ReferenceTopBuilder` and
:class:`ReferenceSubBuilder` to build subcomponents of the reference manual.
The :class:`ReferenceSubBuilder` examines the modules included in the
subcomponent by comparing the modification times of the module files with the
times saved in ``local/share/doctree/reference.pickle`` from the previous
build. Then new rst files are generated for new and updated modules. See
:meth:`~ReferenceSubBuilder.get_new_and_updated_modules`.

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

import functools
import importlib.machinery
import logging
import os
import pickle
import re
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
import threading
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
    format = args[1]
    options = args[2]
    kwds = args[3]
    args = args[4:]
    if format == 'inventory':  # you must not use the inventory to build the inventory
        kwds['use_multidoc_inventory'] = False
    getattr(ReferenceSubBuilder(str(doc), options), format)(*args, **kwds)


##########################################
#             Builders                   #
##########################################

def builder_helper(type):
    """
    Return a function which builds the documentation for
    output type ``type``.
    """
    def f(self, *args, **kwds):
        single_file = self.documents_single_file
        load_lock_acquired = False
        succeeded = False
        try:
            if single_file:
                # The alternate finder remains live after ``load_single_file``
                # returns, while Sphinx imports and inspects the module.  Keep
                # another same-process single-file build from clearing that
                # transaction until this Sphinx lifecycle is complete.  The
                # lock is reentrant because conf.py calls ``load_single_file``
                # synchronously in this thread.
                _SINGLE_FILE_LOAD_LOCK.acquire()
                load_lock_acquired = True
                _clear_single_file_imports()
            output_dir = self._output_dir(type)

            options = list(build_options.ALLSPHINXOPTS)

            # The name carries the language, as in 'en/website'.
            if self.name.split('/')[-1] == 'website':
                options += build_options.WEBSITESPHINXOPTS

            # An inventory build, and the first pass of the reference manual,
            # run before the inventories of the other documents exist.
            first_pass = (type == 'inventory'
                          or not kwds.get('use_multidoc_inventory', True))
            options += ['-D', f'multidoc_first_pass={int(first_pass)}']

            # Cross-references are legitimately unresolvable in a first pass,
            # so nitpicky mode is only meaningful afterwards.
            if build_options.WARN_LINKS and type != 'inventory':
                options.append('-n')

            # Provide the ``pdf`` tag as an alias of ``latex``.
            tags = ['-t', 'pdf'] if type == 'latex' else []

            argv = [*tags, '-b', type, '-d', str(self._doctrees_dir()),
                    *options, str(self.dir), str(output_dir)]

            logger.debug('sphinx-build %s', ' '.join(argv))

            # Run Sphinx with Sage's special logger.  The build type decides
            # how its diagnostics are treated; see
            # :func:`~sage_docbuild.sphinxbuild.runsphinx`.
            from .sphinxbuild import runsphinx
            try:
                runsphinx(argv,
                          prefix=os.path.basename(output_dir),
                          warnings_are_errors=type != 'latex',
                          first_pass=first_pass,
                          is_inventory=type == 'inventory',
                          single_file=self.documents_single_file,
                          single_file_path=getattr(self, 'single_file_path', None),
                          single_file_source_root=getattr(
                              self, 'single_file_source_root', None))
                succeeded = True
            except Exception:
                if build_options.ABORT_ON_ERROR:
                    raise
            except BaseException as e:
                # We need to wrap a BaseException that is not an Exception in a
                # regular Exception. Otherwise multiprocessing.Pool.get hangs, see
                # #25161
                if build_options.ABORT_ON_ERROR:
                    raise Exception("Non-exception during docbuild: %s" % (e,), e)

            if succeeded:
                mark_success = getattr(self, '_mark_build_success', None)
                if mark_success is not None:
                    mark_success()

            if not succeeded:
                return
            if type == 'latex':
                logger.warning(f"LaTeX files can be found in {output_dir}.")
            elif type != 'inventory':
                logger.warning(
                    f"Build finished. The built documents can be found in {output_dir}.")
        finally:
            try:
                if single_file and load_lock_acquired:
                    _clear_single_file_imports()
            finally:
                if load_lock_acquired:
                    _SINGLE_FILE_LOAD_LOCK.release()

    f.is_output_format = True
    f.__name__ = type
    f.__qualname__ = f"DocBuilder.{type}"
    f.__doc__ = f"Build the documentation for output type ``{type}``."
    return f


class DocBuilder():
    #: Whether this builder documents one file of its own rather than a manual;
    #: see :class:`SingleFileBuilder` and :func:`builder_helper`.
    documents_single_file = False

    def __init__(self, name: str, options: BuildOptions):
        """
        INPUT:

        - ``name`` -- the name of a document directory below ``doc``, such as
          'en/tutorial' or 'fr/tutorial'
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


@functools.cache
def _output_formats() -> frozenset:
    """
    Return the names of the output formats that a builder can produce.

    EXAMPLES::

        sage: from sage_docbuild.builders import _output_formats
        sage: sorted(_output_formats())                                 # needs sphinx
        ['changes', 'html', 'htmlhelp', 'inventory', 'json', 'latex',
         'linkcheck', 'pdf', 'pickle', 'web']
    """
    formats = {name for name in dir(DocBuilder)
               if getattr(getattr(DocBuilder, name, None), 'is_output_format', False)}
    # pdf is a method of its own: it builds the latex output and compiles it.
    return frozenset(formats | {'pdf'})


def _library_modules():
    """
    Return an iterator over the modules of the Sage library.

    The modules that no manual documents are left out: a package's
    ``__init__`` and the ``all`` modules gathering a namespace.

    EXAMPLES::

        sage: from sage_docbuild.builders import _library_modules
        sage: modules = set(_library_modules())
        sage: 'sage.graphs.graph' in modules
        True
        sage: any(name.endswith(('__init__', '.all')) for name in modules)
        False
    """
    from sage.env import SAGE_SRC

    base_path = os.path.join(SAGE_SRC, 'sage')
    for directory, subdirs, files in os.walk(base_path):
        for filename in files:
            if not (filename.endswith('.py') or filename.endswith('.pyx')):
                continue

            path = os.path.join(directory, filename)

            # Create the module name
            module_name = path[len(base_path):].replace(os.path.sep, '.')
            module_name = 'sage' + module_name
            module_name = module_name[:-4] if module_name.endswith('pyx') else module_name[:-3]

            # Exclude some ones  -- we don't want init the manual
            if module_name.endswith('__init__') or module_name.endswith('all'):
                continue

            yield module_name


def _reference_commands() -> frozenset:
    """
    Return the commands that can be run over the reference manual.

    A command reports on the modules of a sub-document instead of building it;
    the reference manual runs it on each of its sub-documents in turn.

    EXAMPLES::

        sage: from sage_docbuild.builders import _reference_commands
        sage: sorted(_reference_commands())
        ['print_included_modules',
         'print_modified_modules',
         'print_new_and_updated_modules',
         'print_unincluded_modules']
    """
    return frozenset(name for name in dir(ReferenceSubBuilder)
                     if name.startswith('print_'))


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

    #: The bibliography, which the other manuals reference.
    _bibliography = Path('reference/references')

    def _sub_documents(self) -> list[Path]:
        """
        Return the sub-documents of the reference manual.

        The top-level document is not one of them: it is built last, by
        :meth:`_build_top_level`.
        """
        return [doc
                for doc in get_all_reference_documents(self.options.source_dir / 'en')
                if doc != Path('reference_top')]

    def __getattr__(self, attr):
        """
        Return a function running ``attr`` over the whole reference manual.

        ``attr`` is either an output format of :class:`DocBuilder`, which
        builds each sub-document (see :func:`builder_helper`), or a command of
        :class:`ReferenceSubBuilder` reporting on each of them (see
        :func:`_reference_commands`).

        EXAMPLES::

            sage: from sage_docbuild.builders import ReferenceBuilder
            sage: from sage_docbuild.build_options import BuildOptions
            sage: builder = ReferenceBuilder('reference', BuildOptions())
            sage: builder.html                       # needs sphinx
            functools.partial(<bound method ReferenceBuilder._wrapper of ...>, 'html')
            sage: builder.print_included_modules
            functools.partial(<bound method ReferenceBuilder._wrapper of ...>, 'print_included_modules')

        A command the manual answers for as a whole is a method of its own::

            sage: builder.print_unincluded_modules
            <bound method ReferenceBuilder.print_unincluded_modules of ...>

        Anything else is not an attribute of the builder, so that a mere
        lookup does not answer with a function that fails when called::

            sage: builder.no_such_format
            Traceback (most recent call last):
            ...
            AttributeError: 'ReferenceBuilder' object has no attribute 'no_such_format'
            sage: hasattr(builder, '_ipython_canary_method_should_not_exist_')
            False
        """
        if attr not in _output_formats() and attr not in _reference_commands():
            raise AttributeError(
                f"{type(self).__name__!r} object has no attribute {attr!r}")
        from functools import partial
        return partial(self._wrapper, attr)

    def _build_bibliography(self, format, *args, **kwds):
        """
        Build the bibliography only

        The bibliography references.aux is referenced by the other
        manuals and needs to be built first.
        """
        references = [(self._bibliography, format, self.options, kwds) + args]
        build_many(build_ref_doc, references)

    def _build_everything_except_bibliography(self, format, *args, **kwds):
        """
        Build the entire reference manual except the bibliography
        """
        non_references = [
            (doc, format, self.options, kwds) + args for doc in self._sub_documents()
            if doc != self._bibliography
        ]
        build_many(build_ref_doc, non_references)

    def _build_top_level(self, format, *args, **kwds):
        """
        Build top-level document.
        """
        getattr(ReferenceTopBuilder('reference', self.options), format)(*args, **kwds)

    def print_unincluded_modules(self):
        """
        Print the modules of the Sage library that the reference manual does
        not include.

        A sub-document knows the modules of its own toctrees only, so running
        the command of :class:`ReferenceSubBuilder` over each of them in turn
        would report every module that any *other* sub-document documents. The
        answer is the complement of their union.
        """
        included = set()
        for document in self._sub_documents():
            builder = ReferenceSubBuilder(document.as_posix(), self.options)
            included.update(builder.get_all_included_modules())
        for module_name in _library_modules():
            if module_name not in included:
                print(module_name)

    def _wrapper(self, format, *args, **kwds):
        """
        Run ``format`` over the reference manual: over the top-level document
        and each of its components if it is an output format, over the
        components alone if it is one of the :func:`_reference_commands`, which
        the top-level document does not answer to.
        """
        is_format = format in _output_formats()
        if is_format and format != 'inventory':
            # The sub-documents cross-reference each other through the
            # inventories, so the manual is built twice: the first pass writes
            # the inventories, the second resolves the references against them.
            # Under meson the two passes are separate targets, but nothing else
            # schedules them for a build started here.
            logger.warning('Building reference manual, first pass.\n')
            self._wrapper('inventory', *args, **kwds)
            logger.warning('Building reference manual, second pass.\n')
        logger.info('Building bibliography')
        self._build_bibliography(format, *args, **kwds)
        logger.info('Bibliography finished, building dependent manuals')
        self._build_everything_except_bibliography(format, *args, **kwds)
        if is_format:
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
        from functools import partial
        for attr in dir(self):
            if hasattr(getattr(self, attr), 'is_output_format'):
                f = partial(self._wrapper, attr)
                f.is_output_format = True
                wrapped = getattr(self, attr)
                f.__doc__ = wrapped.__doc__
                f.__name__ = attr
                f.__qualname__ = f"{type(self).__qualname__}.{attr}"
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
        which are not included in this document.

        Beware that a module the reference manual documents elsewhere is not
        included *here*; :meth:`ReferenceBuilder.print_unincluded_modules`
        answers for the manual as a whole.
        """
        included_modules = set(self.get_all_included_modules())
        for module_name in _library_modules():
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


def _module_path(name):
    """
    Return the file that importing ``name`` would read, or ``None``.

    Nothing is imported: a package on the way may well be, but the module
    itself is only looked up.

    EXAMPLES::

        sage: from sage_docbuild.builders import _module_path
        sage: _module_path('sage.graphs.graph').endswith('graphs/graph.py')
        True
        sage: _module_path('no_such_module_2718') is None
        True
    """
    import importlib.util

    try:
        spec = importlib.util.find_spec(name)
    except (AttributeError, ImportError, ValueError):
        return None
    if spec is None or spec.origin is None:
        return None
    return os.path.realpath(spec.origin)


class _SiblingFinder:
    """
    Find the descendants of ``package`` below ``directory``.

    A file documented on its own may be one of a copy of a package that the
    import system does not lead to - a second checkout of Sage, most often -
    and the modules beside it are then the ones it means to import.  Answering
    for them takes a finder of one's own: an entry in ``__path__`` is only
    consulted by the path finder, which the finder of an editable install
    comes before.

    Delegating to :class:`importlib.machinery.PathFinder` makes all of Python's
    ordinary file loaders available, rather than only the source loader.  In
    particular, packages below a sibling package and compiled extension
    modules are found as well.
    """
    def __init__(self, package: str, directory: str):
        self.package = package
        self.directory = os.path.realpath(directory)
        self.seen = set()

    def find_spec(self, fullname, path=None, target=None):
        prefix = self.package + '.'
        if not fullname.startswith(prefix):
            return None
        relative = fullname[len(prefix):].split('.')
        parent = os.path.join(self.directory, *relative[:-1])
        spec = importlib.machinery.PathFinder.find_spec(fullname, [parent])
        if spec is not None:
            self.seen.add(fullname)
        return spec

    def contains(self, module) -> bool:
        r"""
        Return whether ``module`` was loaded from this finder's tree.

        Inspecting an unrelated module does not invoke its PEP 562 hooks::

            sage: from types import ModuleType
            sage: from sage_docbuild.builders import _SiblingFinder
            sage: module = ModuleType('_sage_docbuild_contains_test')
            sage: calls = []
            sage: module.__getattr__ = lambda name: calls.append(name)
            sage: finder = _SiblingFinder('example', '/no/such/tree')
            sage: finder.contains(module), calls
            (False, [])
        """
        if not isinstance(module, type(sys)):
            return False
        # A module can define PEP 562 attribute hooks.  Cleanup must not invoke
        # arbitrary code merely to find where an entry of sys.modules came
        # from, especially while unwinding a failed import.
        namespace = type(sys).__getattribute__(module, '__dict__')
        paths = []
        filename = namespace.get('__file__')
        if filename:
            paths.append(filename)
        search = namespace.get('__path__', ())
        if isinstance(search, (str, bytes, os.PathLike)):
            paths.append(search)
        else:
            try:
                paths.extend(search)
            except Exception:
                pass
        for path in paths:
            try:
                if os.path.commonpath((self.directory, os.path.realpath(path))) == self.directory:
                    return True
            except Exception:
                pass
        return False


def _finder_packages_overlap(left: str, right: str) -> bool:
    """Return whether two alternate-package finders can answer one name."""
    return (left == right or left.startswith(right + '.')
            or right.startswith(left + '.'))


_SINGLE_FILE_IMPORT_LOCK = threading.RLock()
_SINGLE_FILE_LOAD_LOCK = threading.RLock()
_SINGLE_FILE_IMPORT_TRANSACTIONS = []
_MISSING = object()


def _in_import_scope(fullname: str, scope: str) -> bool:
    """Return whether ``fullname`` is ``scope``, its ancestor, or descendant."""
    return (fullname == scope or fullname.startswith(scope + '.')
            or scope.startswith(fullname + '.'))


def _scoped_module_snapshot(scope: str):
    """Snapshot modules and namespaces at and around ``scope``."""
    modules = {}
    namespaces = {}
    for module_name, module in tuple(sys.modules.items()):
        if (not isinstance(module_name, str) or module is None
                or not _in_import_scope(module_name, scope)):
            continue
        modules[module_name] = module
        if isinstance(module, type(sys)):
            namespace = type(sys).__getattribute__(module, '__dict__')
            namespaces[module_name] = namespace.copy()
    return modules, namespaces


class _SingleFileImportTransaction:
    """State displaced while one explicitly loaded file is being documented."""

    def __init__(self, scope: str):
        self.scope = scope
        self.modules_before, self.namespaces_before = (
            _scoped_module_snapshot(scope))
        self.finder = None
        self.removed_finders = []
        self.touched = set()


def _record_changed_scoped_modules(transaction) -> None:
    """Record module changes made while importing the containing package."""
    current, _ = _scoped_module_snapshot(transaction.scope)
    for fullname in transaction.modules_before.keys() | current.keys():
        if (transaction.modules_before.get(fullname, _MISSING)
                is not current.get(fullname, _MISSING)):
            transaction.touched.add(fullname)


def _restore_import_transaction(transaction, *, restore_finders=True) -> None:
    """Restore only the module bindings displaced by ``transaction``."""
    finder = transaction.finder
    if finder is not None:
        transaction.touched.update(finder.seen)
        sys.meta_path[:] = [item for item in sys.meta_path if item is not finder]

    if restore_finders:
        for index, old_finder in sorted(transaction.removed_finders):
            if not any(item is old_finder for item in sys.meta_path):
                sys.meta_path.insert(min(index, len(sys.meta_path)), old_finder)

    # Restore just the names that the parent import, the alternate finder, or
    # the explicit target displaced.  Imports outside this scope, including
    # modules imported by user code before it failed, remain intact.
    for fullname in sorted(
            transaction.touched, key=lambda item: item.count('.'), reverse=True):
        previous = transaction.modules_before.get(fullname, _MISSING)
        if previous is _MISSING:
            sys.modules.pop(fullname, None)
        else:
            sys.modules[fullname] = previous

    # Repair only canonical child bindings.  Clearing whole package
    # dictionaries would discard unrelated attributes created while the target
    # ran; ignoring these bindings would leave ``pkg.child`` inconsistent with
    # ``sys.modules['pkg.child']``.
    for fullname in sorted(transaction.touched, key=lambda item: item.count('.')):
        parent_name, _, child = fullname.rpartition('.')
        if not parent_name:
            continue
        parent = sys.modules.get(parent_name)
        if not isinstance(parent, type(sys)):
            continue
        parent_namespace = type(sys).__getattribute__(parent, '__dict__')
        previous_parent = transaction.modules_before.get(parent_name, _MISSING)
        previous_namespace = transaction.namespaces_before.get(parent_name)
        if parent is previous_parent and previous_namespace is not None:
            if child in previous_namespace:
                parent_namespace[child] = previous_namespace[child]
            else:
                parent_namespace.pop(child, None)
            continue
        desired = sys.modules.get(fullname, _MISSING)
        if desired is not _MISSING:
            parent_namespace[child] = desired
        else:
            current = parent_namespace.get(child, _MISSING)
            if (isinstance(current, type(sys))
                    and type(sys).__getattribute__(current, '__dict__').get(
                        '__name__') == fullname):
                parent_namespace.pop(child, None)


def _clear_single_file_imports(scope=None) -> None:
    r"""
    End active alternate-tree transactions and remove stray finders.

    A stale finder only owns names in its package, not every module whose file
    happens to sit below the same filesystem directory::

        sage: from types import ModuleType
        sage: from sage_docbuild.builders import (_SiblingFinder,
        ....:     _clear_single_file_imports)
        sage: unrelated = ModuleType('_sage_docbuild_unrelated_test')
        sage: unrelated.__file__ = '/tmp/shared-tree/tool.py'
        sage: sys.modules[unrelated.__name__] = unrelated
        sage: finder = _SiblingFinder('_sage_docbuild_package_test',
        ....:                         '/tmp/shared-tree')
        sage: sys.meta_path.insert(0, finder)
        sage: _clear_single_file_imports()
        sage: sys.modules[unrelated.__name__] is unrelated
        True
        sage: del sys.modules[unrelated.__name__]
    """
    # Serializing explicit loads does not block ordinary imports.  It only
    # prevents a second single-file transaction from taking a baseline while
    # the first is temporarily executing against alternate module bindings.
    with _SINGLE_FILE_LOAD_LOCK:
        with _SINGLE_FILE_IMPORT_LOCK:
            for transaction in reversed(
                    _SINGLE_FILE_IMPORT_TRANSACTIONS.copy()):
                if (scope is not None
                        and not _finder_packages_overlap(
                            transaction.scope, scope)):
                    continue
                _restore_import_transaction(transaction)
                _SINGLE_FILE_IMPORT_TRANSACTIONS.remove(transaction)

            # A finder can predate this transaction machinery (or survive an
            # interrupted caller).  It has no trustworthy snapshot to restore,
            # but both it and modules visibly loaded from its tree must stop
            # influencing a later active-tree build.
            stale_finders = [
                item for item in sys.meta_path
                if (isinstance(item, _SiblingFinder)
                    and (scope is None
                         or _finder_packages_overlap(item.package, scope)))
            ]
            sys.meta_path[:] = [item for item in sys.meta_path
                                if item not in stale_finders]
            stale_modules = {
                fullname for fullname, module in tuple(sys.modules.items())
                if isinstance(fullname, str) and module is not None
                and any(
                    (fullname == finder.package
                     or fullname.startswith(finder.package + '.'))
                    and finder.contains(module)
                    for finder in stale_finders)
            }
            _drop_modules(stale_modules)


def _serialized_single_file_load(function):
    """Serialize explicit loads without taking Python's global import lock."""
    @functools.wraps(function)
    def wrapped(*args, **kwds):
        with _SINGLE_FILE_LOAD_LOCK:
            return function(*args, **kwds)
    return wrapped


def _shadowed_module_names(finder, old_finders, target, package_is_target):
    """Return loaded module names that an alternate tree must replace."""
    package = finder.package
    names = {target}
    modules = tuple(sys.modules.items())
    for fullname, module in modules:
        if not isinstance(fullname, str):
            continue
        if fullname == package:
            if package_is_target:
                names.add(fullname)
            continue
        if not fullname.startswith(package + '.'):
            continue
        if package_is_target:
            names.add(fullname)
            continue
        if finder.find_spec(fullname) is not None:
            names.add(fullname)
            continue
        if module is not None and any(old.contains(module) for old in old_finders):
            names.add(fullname)

    # A previous import can leave a child on its parent package even after its
    # sys.modules entry was removed.  Account for those stale attributes too.
    for _, parent in modules:
        if not isinstance(parent, type(sys)):
            continue
        namespace = type(sys).__getattribute__(parent, '__dict__')
        for value in tuple(namespace.values()):
            if not isinstance(value, type(sys)):
                continue
            fullname = type(sys).__getattribute__(value, '__dict__').get(
                '__name__', '')
            if (not isinstance(fullname, str)
                    or not (fullname == package
                            or fullname.startswith(package + '.'))):
                continue
            if (package_is_target or fullname == target
                    or finder.find_spec(fullname) is not None
                    or any(old.contains(value) for old in old_finders)):
                names.add(fullname)
    return names


def _drop_modules(names) -> None:
    """Drop modules in ``names`` and their attributes on parent packages."""
    for fullname in sorted(names, key=lambda item: item.count('.'), reverse=True):
        parent_name, _, child = fullname.rpartition('.')
        parent = sys.modules.get(parent_name)
        namespace = (type(sys).__getattribute__(parent, '__dict__')
                     if isinstance(parent, type(sys)) else None)
        value = namespace.get(child) if namespace is not None else None
        if (isinstance(value, type(sys))
                and type(sys).__getattribute__(value, '__dict__').get(
                    '__name__') == fullname):
            namespace.pop(child, None)
        sys.modules.pop(fullname, None)


@_serialized_single_file_load
def load_single_file(name: str, path: str):
    r"""
    Import the file ``path`` under the module name ``name``, and return it.

    This is what the configuration of a single-file build runs when importing
    ``name`` would read another file - the same module of a second checkout,
    say - or no file at all.  The modules beside ``path`` come before the ones
    of that other copy, so that what the file imports is what sits next to it;
    a module of that copy which is loaded already is dropped, since an import
    would answer with it and never read this directory.

    Modules that ``path`` reaches by any other route are the ones already
    loaded, which is what the caller warns about.

    For an ordinary module, the already-active containing package is retained;
    its initializer is not re-executed from the alternate tree.  When ``path``
    is itself ``__init__.py``, that package initializer is the explicit target
    and is replaced together with its children.

    EXAMPLES::

        sage: import os, sys, tempfile
        sage: from sage_docbuild.builders import (load_single_file,
        ....:     _clear_single_file_imports)
        sage: with tempfile.TemporaryDirectory() as directory:
        ....:     package = os.path.join(directory, 'pkg_2718')
        ....:     os.mkdir(package)
        ....:     _ = open(os.path.join(package, '__init__.py'), 'w').close()
        ....:     with open(os.path.join(package, 'helper.py'), 'w') as f:
        ....:         _ = f.write('answer = 42\n')
        ....:     with open(os.path.join(package, 'main.py'), 'w') as f:
        ....:         _ = f.write('import _imp\nfrom .helper import answer\n'
        ....:                     'lock_held = _imp.lock_held()\n')
        ....:     sys.path.insert(0, directory)
        ....:     module = load_single_file('pkg_2718.main',
        ....:                               os.path.join(package, 'main.py'))
        ....:     _clear_single_file_imports()
        ....:     sys.path.remove(directory)
        sage: module.answer
        42
        sage: module.__name__
        'pkg_2718.main'
        sage: module.lock_held
        False

    Already-loaded children and descendants of another copy are replaced, and
    a failed replacement restores them::

        sage: import importlib
        sage: from sage_docbuild.builders import (_SiblingFinder,
        ....:     _module_path)
        sage: with tempfile.TemporaryDirectory() as active, \
        ....:      tempfile.TemporaryDirectory() as alternate:
        ....:     active_pkg = Path(active) / 'pkg_shadow_2718'
        ....:     alternate_pkg = Path(alternate) / 'pkg_shadow_2718'
        ....:     (active_pkg / 'nested').mkdir(parents=True)
        ....:     (alternate_pkg / 'nested').mkdir(parents=True)
        ....:     for package in (active_pkg, alternate_pkg,
        ....:                     active_pkg / 'nested', alternate_pkg / 'nested'):
        ....:         _ = (package / '__init__.py').write_text('')
        ....:     _ = (active_pkg / 'helper.py').write_text("VALUE = 'active'\n")
        ....:     _ = (active_pkg / 'nested' / 'helper.py').write_text(
        ....:         "VALUE = 'active nested'\n")
        ....:     _ = (active_pkg / 'broken.py').write_text("VALUE = 'old'\n")
        ....:     _ = (active_pkg / 'fallback.py').write_text(
        ....:         "VALUE = 'active fallback'\n")
        ....:     _ = (Path(active) / 'unrelated_tx_2718.py').write_text(
        ....:         "VALUE = 'preserved'\n")
        ....:     _ = (alternate_pkg / 'helper.py').write_text(
        ....:         "VALUE = 'alternate'\n")
        ....:     _ = (alternate_pkg / 'nested' / 'helper.py').write_text(
        ....:         "VALUE = 'alternate nested'\n")
        ....:     _ = (alternate_pkg / 'main.py').write_text(
        ....:         'from . import helper\n'
        ....:         'from .nested.helper import VALUE as NESTED\n'
        ....:         'VALUE = helper.VALUE\n')
        ....:     sys.path.insert(0, active)
        ....:     package = importlib.import_module('pkg_shadow_2718')
        ....:     old_helper = importlib.import_module('pkg_shadow_2718.helper')
        ....:     old_nested = importlib.import_module(
        ....:         'pkg_shadow_2718.nested.helper')
        ....:     main = load_single_file(
        ....:         'pkg_shadow_2718.main', str(alternate_pkg / 'main.py'))
        ....:     values = (main.VALUE, main.NESTED,
        ....:               package.helper is not old_helper,
        ....:               sys.modules['pkg_shadow_2718.nested.helper'] is not old_nested)
        ....:     _clear_single_file_imports()
        ....:     active_restored = (
        ....:         package.helper is old_helper
        ....:         and sys.modules['pkg_shadow_2718.nested.helper'] is old_nested
        ....:         and _module_path('pkg_shadow_2718.helper')
        ....:             == os.path.realpath(active_pkg / 'helper.py')
        ....:         and not any(isinstance(item, _SiblingFinder)
        ....:                     and item.package == 'pkg_shadow_2718'
        ....:                     for item in sys.meta_path))
        ....:     old_broken = importlib.import_module('pkg_shadow_2718.broken')
        ....:     _ = sys.modules.pop('unrelated_tx_2718', None)
        ....:     _ = (alternate_pkg / 'broken.py').write_text(
        ....:         "import unrelated_tx_2718\n"
        ....:         "from . import fallback\n"
        ....:         "VALUE = 'partial'\nraise RuntimeError('broken')\n")
        ....:     meta_path = list(sys.meta_path)
        ....:     try:
        ....:         load_single_file('pkg_shadow_2718.broken',
        ....:                          str(alternate_pkg / 'broken.py'))
        ....:     except RuntimeError:
        ....:         restored = (sys.modules['pkg_shadow_2718.broken'] is old_broken
        ....:                     and package.broken is old_broken
        ....:                     and sys.meta_path == meta_path)
        ....:     unrelated_preserved = 'unrelated_tx_2718' in sys.modules
        ....:     fallback_preserved = (
        ....:         package.fallback
        ....:         is sys.modules['pkg_shadow_2718.fallback'])
        ....:     _ = sys.modules.pop('unrelated_tx_2718', None)
        ....:     sys.path.remove(active)
        ....:     for loaded in list(sys.modules):
        ....:         if loaded == 'pkg_shadow_2718' or loaded.startswith(
        ....:                 'pkg_shadow_2718.'):
        ....:             del sys.modules[loaded]
        ....:     sys.meta_path[:] = [finder for finder in sys.meta_path
        ....:                         if not (isinstance(finder, _SiblingFinder)
        ....:                                 and finder.package == 'pkg_shadow_2718')]
        sage: (values, active_restored, restored, unrelated_preserved,
        ....:  fallback_preserved)
        (('alternate', 'alternate nested', True, True), True, True, True, True)

    A package initializer replaces its loaded children as well::

        sage: with tempfile.TemporaryDirectory() as active, \
        ....:      tempfile.TemporaryDirectory() as alternate:
        ....:     for root, value in ((active, 'active'),
        ....:                         (alternate, 'alternate')):
        ....:         package = Path(root) / 'pkg_init_2718'
        ....:         package.mkdir()
        ....:         _ = (package / '__init__.py').write_text(
        ....:             'from .helper import VALUE\n')
        ....:         _ = (package / 'helper.py').write_text(
        ....:             f"VALUE = {value!r}\n")
        ....:     sys.path.insert(0, active)
        ....:     old_package = importlib.import_module('pkg_init_2718')
        ....:     old_helper = importlib.import_module('pkg_init_2718.helper')
        ....:     package = load_single_file(
        ....:         'pkg_init_2718',
        ....:         str(Path(alternate) / 'pkg_init_2718' / '__init__.py'))
        ....:     answer = (package.VALUE,
        ....:               sys.modules['pkg_init_2718.helper'] is not old_helper)
        ....:     _clear_single_file_imports()
        ....:     sys.path.remove(active)
        ....:     for loaded in list(sys.modules):
        ....:         if loaded == 'pkg_init_2718' or loaded.startswith(
        ....:                 'pkg_init_2718.'):
        ....:             del sys.modules[loaded]
        ....:     sys.meta_path[:] = [finder for finder in sys.meta_path
        ....:                         if not (isinstance(finder, _SiblingFinder)
        ....:                                 and finder.package == 'pkg_init_2718')]
        sage: answer
        ('alternate', True)

    If importing a previously absent containing package is part of a failed
    load, that package and its canonical child binding are rolled back too::

        sage: with tempfile.TemporaryDirectory() as directory:
        ....:     package = Path(directory) / 'pkg_parent_2718'
        ....:     package.mkdir()
        ....:     _ = (package / '__init__.py').write_text(
        ....:         'from . import side\n')
        ....:     _ = (package / 'side.py').write_text('VALUE = 42\n')
        ....:     broken = package / 'broken.py'
        ....:     _ = broken.write_text("raise RuntimeError('broken')\n")
        ....:     sys.path.insert(0, directory)
        ....:     try:
        ....:         load_single_file('pkg_parent_2718.broken', str(broken))
        ....:     except RuntimeError:
        ....:         rolled_back = not any(
        ....:             loaded == 'pkg_parent_2718'
        ....:             or loaded.startswith('pkg_parent_2718.')
        ....:             for loaded in sys.modules)
        ....:     sys.path.remove(directory)
        sage: rolled_back
        True
    """
    import importlib
    import importlib.util

    path = os.path.abspath(path)
    directory = os.path.dirname(path)
    package_is_target = os.path.basename(path) == '__init__.py'
    containing_package = name.rpartition('.')[0]
    package = name if package_is_target else containing_package
    scope = package or name

    # A preceding explicit load of this package belongs to an earlier build.
    # Restore its active-tree baseline before taking the next snapshot.
    _clear_single_file_imports(scope)
    with _SINGLE_FILE_IMPORT_LOCK:
        # This snapshot deliberately precedes the parent import.  Importing a
        # previously absent package is itself part of the transaction.
        transaction = _SingleFileImportTransaction(scope)

    parent_import_recorded = False
    try:
        # Import the containing package before replacing anything below it.  A
        # package initializer is itself replaced, so only its parent is kept.
        parent_to_import = containing_package if package_is_target else package
        if parent_to_import:
            importlib.import_module(parent_to_import)

        # Only shared import-state surgery is serialized.  In particular, the
        # target's top-level code below runs without Python's global import lock
        # and without this dedicated state lock.
        with _SINGLE_FILE_IMPORT_LOCK:
            _record_changed_scoped_modules(transaction)
            parent_import_recorded = True
            importlib.invalidate_caches()
            if package:
                finder = _SiblingFinder(package, directory)
                transaction.finder = finder
                transaction.removed_finders = [
                    (index, item) for index, item in enumerate(sys.meta_path)
                    if (isinstance(item, _SiblingFinder)
                        and _finder_packages_overlap(item.package, package))
                ]
                old_finders = [item for _, item
                               in transaction.removed_finders]
                sys.meta_path[:] = [
                    item for item in sys.meta_path
                    if not any(item is old for old in old_finders)
                ]
                sys.meta_path.insert(0, finder)
                shadowed = _shadowed_module_names(
                    finder, old_finders, name, package_is_target)
            else:
                shadowed = {name}
            transaction.touched.update(shadowed)
            _drop_modules(shadowed)

            spec = importlib.util.spec_from_file_location(name, path)
            if spec is None or spec.loader is None:
                raise ImportError(f'cannot load {name!r} from {path!r}')
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            parent_name, _, child = name.rpartition('.')
            parent = sys.modules.get(parent_name)
            if isinstance(parent, type(sys)):
                namespace = type(sys).__getattribute__(parent, '__dict__')
                namespace[child] = module

        spec.loader.exec_module(module)
    except BaseException:
        with _SINGLE_FILE_IMPORT_LOCK:
            if not parent_import_recorded:
                _record_changed_scoped_modules(transaction)
            _restore_import_transaction(transaction)
        raise

    with _SINGLE_FILE_IMPORT_LOCK:
        transaction.touched.add(name)
        if transaction.finder is not None:
            transaction.touched.update(transaction.finder.seen)
        _SINGLE_FILE_IMPORT_TRANSACTIONS.append(transaction)
    return module


def _extend_over_namespace_packages(directory, parts):
    r"""
    Carry the walk of :func:`_single_file_module` over namespace packages.

    A namespace package has no initializer, so the walk up the directories
    stops below it and the module ends up named as if it stood alone, which a
    relative import in it does not survive.  Such a package is only a directory
    that a directory on :data:`sys.path` leads to, so the walk goes on for as
    long as one does, and no further: a directory that nothing imports from
    names nothing.

    This intentionally does not guess a namespace root for an arbitrary path:
    any number of its identifier-named parent directories could be namespace
    packages.  A single-file build of such a module must therefore have its
    import root on :data:`sys.path` already (for example through
    :envvar:`PYTHONPATH`).

    EXAMPLES::

        sage: import os, sys, tempfile
        sage: from sage_docbuild.builders import _extend_over_namespace_packages
        sage: with tempfile.TemporaryDirectory() as root:
        ....:     inner = os.path.join(root, 'ns_2718', 'inner')
        ....:     os.makedirs(inner)
        ....:     sys.path.insert(0, root)
        ....:     answer = _extend_over_namespace_packages(inner, ['mod'])
        ....:     sys.path.remove(root)
        sage: answer[0] == root, answer[1]
        (True, ['mod', 'inner', 'ns_2718'])

    A directory that no entry of ``sys.path`` leads to is left alone::

        sage: _extend_over_namespace_packages('/no/such/tree', ['mod'])
        ('/no/such/tree', ['mod'])
    """
    if not parts:
        return directory, parts
    roots = {os.path.realpath(entry or os.curdir) for entry in sys.path}
    walked = []
    current = os.path.realpath(directory)
    while current not in roots:
        parent, package = os.path.split(current)
        if not package or parent == current or not package.isidentifier():
            return directory, parts   # no entry of sys.path leads here
        walked.append(package)
        current = parent
    return current, parts + walked


def _single_file_module(path):
    r"""
    Return how the file ``path`` has to be imported to be documented.

    OUTPUT:

    a triple ``(name, directory, from_path)``, where ``name`` is the name to
    hand to the ``automodule`` directive, ``directory`` is the directory that
    makes the import work when added to :data:`sys.path`, and ``from_path``
    tells whether the module has to be loaded from ``path`` explicitly instead
    of by name.

    A file that belongs to a package is imported under its qualified name. Its
    base name alone would be ambiguous, and answering it is not even enough:
    ``sage/graphs/generators/random.py`` would be documented as the
    :mod:`random` module of the standard library, which any documentation
    build has imported long before it gets here.

    EXAMPLES::

        sage: import os
        sage: from sage.env import SAGE_SRC
        sage: from sage_docbuild.builders import _single_file_module
        sage: _single_file_module(os.path.join(SAGE_SRC, 'sage', 'graphs',
        ....:                                  'generators', 'random.py'))
        ('sage.graphs.generators.random', '...', False)

    A file outside of any package is loaded from its path, under a name that
    nothing else answers to::

        sage: import tempfile
        sage: with tempfile.TemporaryDirectory() as directory:
        ....:     path = os.path.join(directory, 'random.py')
        ....:     open(path, 'w').close()
        ....:     _single_file_module(path)[0::2]
        ('single_file_random', True)
        sage: with tempfile.TemporaryDirectory() as directory:
        ....:     path = os.path.join(directory, 'no_such_module_2718.py')
        ....:     open(path, 'w').close()
        ....:     _single_file_module(path)[0::2]
        ('no_such_module_2718', True)

    A filename need not itself be a Python identifier::

        sage: with tempfile.TemporaryDirectory() as directory:
        ....:     path = os.path.join(directory, 'module.part-2718.py')
        ....:     open(path, 'w').close()
        ....:     _single_file_module(path)[0::2]
        ('single_file_module_part_2718', True)

    The initializer of a package is the package itself, not a ``__init__``
    submodule of it, which would run the initializer a second time under a
    name of its own::

        sage: _single_file_module(os.path.join(SAGE_SRC, 'sage', 'graphs',
        ....:                                  '__init__.py'))[0]
        'sage.graphs'
    """
    import importlib.util

    directory, filename = os.path.split(path)
    name = os.path.splitext(filename)[0]
    # The initializer of a package carries the name of the package.
    initializer = name == '__init__'
    if not initializer and not name.isidentifier():
        name = 'single_file_' + re.sub(r'\W', '_', name)
    parts = [] if initializer else [name]
    while os.path.exists(os.path.join(directory, '__init__.py')):
        directory, package = os.path.split(directory)
        if not package:
            break
        if not package.isidentifier():
            raise ValueError(
                f'invalid package directory for single-file documentation: '
                f'{package!r}')
        parts.append(package)
    directory, parts = _extend_over_namespace_packages(directory, parts)
    if len(parts) > 1:
        return '.'.join(reversed(parts)), directory, False
    if parts and initializer:
        # A package that is not itself inside one.
        return parts[0], directory, False

    if name not in sys.modules:
        try:
            taken = importlib.util.find_spec(name) is not None
        except (ImportError, ValueError):
            taken = False
    else:
        taken = True
    if taken:
        name = 'single_file_' + name
    return name, directory, True


_SINGLE_FILE_OWNER = '.sage-docbuild-owned'
_SINGLE_FILE_OWNER_CONTENT = 'sage-docbuild single-file output version 1\n'
_SINGLE_FILE_SUCCESS = 'documented-file'


def _validated_single_file(path) -> Path:
    r"""
    Return ``path`` as an absolute path to a regular Python source file.

    Validation happens before a :class:`SingleFileBuilder` creates, removes,
    or stamps an output directory.  In particular, a directory or a path with
    an empty base name cannot turn the output root itself into the build
    directory.

    EXAMPLES::

        sage: import tempfile
        sage: from sage_docbuild.builders import _validated_single_file
        sage: with tempfile.TemporaryDirectory() as directory:
        ....:     path = Path(directory) / 'module.py'
        ....:     _ = path.write_text('answer = 42\n')
        ....:     _validated_single_file(path) == path.absolute()
        True
        sage: _validated_single_file('/')
        Traceback (most recent call last):
        ...
        ValueError: single-file documentation requires an existing regular .py file: /

    A directory containing ``__init__.py`` must have a usable package name;
    only the standalone filename is assigned a synthetic name::

        sage: with tempfile.TemporaryDirectory() as directory:
        ....:     package = Path(directory) / 'bad.package'
        ....:     package.mkdir()
        ....:     _ = (package / '__init__.py').write_text('')
        ....:     module = package / 'module.py'
        ....:     _ = module.write_text('answer = 42\n')
        ....:     _validated_single_file(module)
        Traceback (most recent call last):
        ...
        ValueError: invalid package directory for single-file documentation: 'bad.package'
    """
    path = Path(os.path.abspath(os.fspath(path)))
    valid = path.suffix == '.py'
    if valid:
        try:
            valid = stat.S_ISREG(path.stat().st_mode)
        except OSError:
            valid = False
    if not valid:
        raise ValueError(
            'single-file documentation requires an existing regular .py file: '
            f'{path}')
    directory = path.parent
    while (directory / '__init__.py').exists():
        if not directory.name.isidentifier():
            raise ValueError(
                f'invalid package directory for single-file documentation: '
                f'{directory.name!r}')
        directory = directory.parent
    return path


def _read_regular_text(path: Path):
    """Read a small regular file without following a symbolic link."""
    try:
        before = path.lstat()
    except OSError:
        return None
    if not stat.S_ISREG(before.st_mode):
        return None
    flags = os.O_RDONLY | getattr(os, 'O_NOFOLLOW', 0)
    try:
        descriptor = os.open(path, flags)
    except OSError:
        return None
    try:
        after = os.fstat(descriptor)
        if ((before.st_dev, before.st_ino) != (after.st_dev, after.st_ino)
                or not stat.S_ISREG(after.st_mode)):
            return None
        with os.fdopen(descriptor, encoding='utf-8') as file:
            descriptor = None
            text = file.read(4097)
        if len(text) > 4096:
            return None
        return text
    except (OSError, UnicodeError):
        return None
    finally:
        if descriptor is not None:
            os.close(descriptor)


def _single_file_output_owned(directory: Path) -> bool:
    r"""
    Return whether ``directory`` carries a regular docbuild ownership marker.

    A symbolic link is never accepted as a marker.  The ``documented-file``
    success stamp is deliberately insufficient: an unowned directory can
    contain an attacker-chosen file of that name.

    EXAMPLES::

        sage: import tempfile
        sage: from sage_docbuild.builders import (_single_file_output_owned,
        ....:     _SINGLE_FILE_OWNER, _SINGLE_FILE_OWNER_CONTENT,
        ....:     _SINGLE_FILE_SUCCESS)
        sage: with tempfile.TemporaryDirectory() as root:
        ....:     directory = Path(root) / 'output'
        ....:     directory.mkdir()
        ....:     _ = (directory / _SINGLE_FILE_SUCCESS).write_text(
        ....:         '/tmp/module.py\n')
        ....:     legacy_is_owned = _single_file_output_owned(directory)
        ....:     marker = directory / _SINGLE_FILE_OWNER
        ....:     _ = marker.write_text(_SINGLE_FILE_OWNER_CONTENT)
        ....:     answer = _single_file_output_owned(directory)
        sage: legacy_is_owned, answer
        (False, True)
    """
    return (_read_regular_text(directory / _SINGLE_FILE_OWNER)
            == _SINGLE_FILE_OWNER_CONTENT)


def _prepare_single_file_output(directory: Path, *, explicit: bool) -> None:
    r"""
    Create a fresh owned output directory without following a symlink at the
    output leaf. Parent components use ordinary filesystem path resolution.

    An explicit directory that has no ownership marker is left untouched::

        sage: import stat, tempfile
        sage: from sage_docbuild.builders import (_prepare_single_file_output,
        ....:     _single_file_output_owned, _SINGLE_FILE_SUCCESS)
        sage: with tempfile.TemporaryDirectory() as root:
        ....:     directory = Path(root) / 'module'
        ....:     directory.mkdir()
        ....:     sentinel = directory / 'valuable.txt'
        ....:     _ = sentinel.write_text('keep me')
        ....:     _ = (directory / _SINGLE_FILE_SUCCESS).write_text(
        ....:         '/tmp/module.py\n')
        ....:     try:
        ....:         _prepare_single_file_output(directory, explicit=True)
        ....:     except FileExistsError:
        ....:         preserved = sentinel.read_text()
        sage: preserved
        'keep me'

    A new directory is marked as owned before it becomes visible, and carries
    the permissions of an ordinary directory rather than ``mkdtemp``'s
    private mode::

        sage: with tempfile.TemporaryDirectory() as root:
        ....:     directory = Path(root) / 'module'
        ....:     ordinary = Path(root) / 'ordinary'
        ....:     ordinary.mkdir()
        ....:     _prepare_single_file_output(directory, explicit=True)
        ....:     owned = _single_file_output_owned(directory)
        ....:     same_mode = (stat.S_IMODE(directory.stat().st_mode)
        ....:                  == stat.S_IMODE(ordinary.stat().st_mode))
        sage: owned, same_mode
        (True, True)

    The output root is created when it does not exist yet::

        sage: with tempfile.TemporaryDirectory() as root:
        ....:     directory = Path(root) / 'new-root' / 'module'
        ....:     _prepare_single_file_output(directory, explicit=False)
        ....:     answer = (directory.is_dir(),
        ....:               _single_file_output_owned(directory))
        sage: answer
        (True, True)

    The same holds for an explicitly requested output root::

        sage: with tempfile.TemporaryDirectory() as root:
        ....:     directory = Path(root) / 'explicit-root' / 'module'
        ....:     _prepare_single_file_output(directory, explicit=True)
        ....:     answer = (directory.is_dir(),
        ....:               _single_file_output_owned(directory))
        sage: answer
        (True, True)
    """
    directory.parent.mkdir(parents=True, exist_ok=True)
    try:
        status = directory.lstat()
    except FileNotFoundError:
        status = None
    if status is not None:
        if not stat.S_ISDIR(status.st_mode):
            raise FileExistsError(
                f'refusing to replace non-directory single-file output {directory}')
        if explicit and not _single_file_output_owned(directory):
            raise FileExistsError(
                f'refusing to remove unowned single-file output directory '
                f'{directory}; remove it or choose another -o directory')
        shutil.rmtree(directory)

    temporary = Path(tempfile.mkdtemp(
        prefix=f'.{directory.name}.preparing-', dir=directory.parent))
    try:
        with (temporary / _SINGLE_FILE_OWNER).open('x', encoding='utf-8') as file:
            file.write(_SINGLE_FILE_OWNER_CONTENT)
        # ``mkdtemp`` always uses 0700.  Probe the mode that a normal mkdir
        # receives under the process umask without temporarily changing that
        # process-global umask.
        mode_probe = temporary / '.mkdir-mode'
        mode_probe.mkdir()
        normal_mode = stat.S_IMODE(mode_probe.stat().st_mode)
        mode_probe.rmdir()
        temporary.chmod(normal_mode)
        temporary.rename(directory)
    except BaseException:
        shutil.rmtree(temporary)
        raise


class SingleFileBuilder(DocBuilder):
    """
    This is the class used to build the documentation for a single
    user-specified file. If the file is called 'foo.py', then the
    documentation is built in ``DIR/foo/`` if the user passes the
    command line option "-o DIR", or in ``DOT_SAGE/docbuild/foo/``
    otherwise.
    """
    documents_single_file = True

    def __init__(self, path: str, options: BuildOptions):
        r"""
        INPUT:

        - ``path`` -- the path to the file for which documentation
          should be built

        - ``options`` -- the build options
        """
        from sage.env import DOT_SAGE

        self.lang = 'en'
        self.name = 'single_file'
        self._options = options
        path_object = _validated_single_file(path)
        path = str(path_object)
        self.single_file_path = path_object
        base_name = os.path.splitext(os.path.basename(path))[0]
        explicit_output = getattr(options, 'output_dir_given', False)
        if explicit_output:
            base_dir = Path(options.output_dir) / base_name
        else:
            dot_sage = Path(DOT_SAGE)
            # Match sage.misc.misc's protection for history and configuration;
            # as there, do not change the permissions of an existing DOT_SAGE.
            dot_sage.mkdir(mode=0o700, parents=True, exist_ok=True)
            base_dir = dot_sage / 'docbuild' / base_name
        self._single_file_base_dir = base_dir

        # Create docbuild and relevant subdirectories, e.g.,
        # the static and templates directories in the output directory.
        # By default, this is DOT_SAGE/docbuild/MODULE_NAME, but can
        # also be specified at the command line.

        # ``_module_path`` must see the active import tree, not a finder
        # retained by an earlier single-file build in this process.  Resetting
        # that shared state is itself serialized; see
        # :func:`_clear_single_file_imports`.
        _clear_single_file_imports()
        module_name, module_dir, from_path = _single_file_module(path)
        self.single_file_source_root = Path(module_dir).absolute()
        latex_name = module_name.replace('_', r'\\_')

        if not from_path:
            # A module of a package is imported under its qualified name,
            # and a relative import in it only works that way. Should
            # importing that name read another file - a second checkout,
            # the installed copy, or nothing at all when the module is new
            # - the file that was asked for is loaded under the name
            # explicitly instead.
            found = _module_path(module_name)
            if found != os.path.realpath(path):
                logger.warning(
                    'Warning: importing %s reads %s, not the file given. '
                    'The file given is the one documented, and the modules '
                    'beside it are read from there as well; anything else '
                    'it imports is what this process has loaded already.',
                    module_name, found or 'no file')
                from_path = True

        _prepare_single_file_output(base_dir, explicit=explicit_output)
        self.dir = os.fspath(base_dir / 'source')

        os.makedirs(os.path.join(self.dir, "static"), exist_ok=True)
        os.makedirs(os.path.join(self.dir, "templates"), exist_ok=True)

        if from_path:
            # Importing the name would not read this file - either nothing
            # answers to it, or another copy of the package does - so hand the
            # file itself to the import system, under the name that the members
            # it defines carry.  A package initializer keeps its search path:
            # spec_from_file_location() takes it from the loader.
            load = """
from sage_docbuild.builders import load_single_file
load_single_file({name!r}, {path!r})
""".format(name=module_name, path=path)
        else:
            load = ''

        # Write self.dir/conf.py.  The directory that makes the import work
        # goes in front of sys.path, so that a package reached only from there
        # is read from there and not from a copy that shares its name.
        conf = r"""# This file is automatically generated by {}, do not edit!

import sys, os, contextlib
sys.path.append({!r})
sys.path.insert(0, {!r})
{}
from sage_docbuild.conf import *
html_static_path = [] + html_common_static_path

project = 'Documentation for {}'
release = 'unknown'
name = {!r}
html_title = project
html_short_title = project
htmlhelp_basename = name

with contextlib.suppress(ValueError):
    extensions.remove('sage_docbuild.ext.multidocs') # see #29651
    extensions.remove('sage_docbuild.ext.inventory_builder')

latex_domain_indices = False
latex_documents = [
  ('index', name + '.tex', 'Documentation for {}',
   'unknown', 'manual'),
]
""".format(__file__, self.dir, module_dir, load, module_name, module_name, latex_name)

        # Note that the members to document are selected by the extension
        # sage_docbuild.ext.members, which reads SAGE_DOC_UNDERSCORE itself:
        # defining setup() here would override the one of sage_docbuild.conf.

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

        # Create link from original file to self.dir, so that the modules that
        # sit next to it can be imported.  Note that we append self.dir to
        # sys.path in conf.py.  An output directory named on the command line
        # is reused as it is, so a link left there by an earlier build has to
        # be replaced: it may well point at another file of the same name.
        link = os.path.join(self.dir, os.path.basename(path))
        try:
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(path, link)
        except OSError as error:
            logger.warning('Warning: could not link %s to %s: %s', path, link, error)

    def pdf(self):
        """Build a PDF, recording success only once the whole build ends."""
        self._single_file_defer_completion = True
        succeeded = False
        try:
            super().pdf()
            succeeded = True
        finally:
            # The LaTeX phase uses ``builder_helper``.  Its ordinary completion
            # hook is deferred so that a failed ``make all-pdf`` does not record
            # the file as documented on the strength of the LaTeX phase alone.
            self._single_file_defer_completion = False
            if succeeded:
                self._mark_build_success()

    def _mark_build_success(self):
        """Atomically record the documented file after Sphinx succeeds."""
        if getattr(self, '_single_file_defer_completion', False):
            return
        descriptor, temporary = tempfile.mkstemp(
            prefix='.documented-file-', dir=self._single_file_base_dir)
        try:
            with os.fdopen(descriptor, 'w', encoding='utf-8') as file:
                descriptor = None
                file.write(os.fspath(self.single_file_path) + '\n')
            os.replace(temporary,
                       self._single_file_base_dir / _SINGLE_FILE_SUCCESS)
        finally:
            if descriptor is not None:
                os.close(descriptor)
            try:
                os.unlink(temporary)
            except FileNotFoundError:
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
    r"""
    Return an appropriate *Builder* object for the document ``name``.

    DocBuilder and its subclasses do all the real work in building the
    documentation.

    The name is one of those that :func:`get_documents` lists, so an English
    document may be named without its ``en/`` prefix.

    EXAMPLES::

        sage: from sage.env import SAGE_DOC_SRC
        sage: from sage_docbuild.build_options import BuildOptions
        sage: from sage_docbuild.builders import get_builder
        sage: options = BuildOptions(source_dir=Path(SAGE_DOC_SRC))
        sage: get_builder('developer', options).name
        'en/developer'
        sage: get_builder('en/developer', options).name
        'en/developer'
        sage: get_builder('fr/tutorial', options).name
        'fr/tutorial'

    The website has a builder of its own, and reaches it under either name::

        sage: type(get_builder('website', options)).__name__
        'WebsiteBuilder'
        sage: get_builder('website', options).name
        'en/website'

    A document merely ending in ``reference`` is not mistaken for the Sage
    reference manual::

        sage: import contextlib, io, tempfile
        sage: with tempfile.TemporaryDirectory() as directory:
        ....:     source = Path(directory)
        ....:     document = source / 'fr' / 'myreference'
        ....:     document.mkdir(parents=True)
        ....:     _ = (document / 'index.rst').write_text('Title\n=====\n')
        ....:     english = source / 'en'
        ....:     english.mkdir()
        ....:     _ = (english / 'index.rst').write_text('English\n=======\n')
        ....:     local = BuildOptions(source_dir=source)
        ....:     answer = get_builder('fr/myreference', local)
        ....:     with contextlib.redirect_stdout(io.StringIO()):
        ....:         try:
        ....:             get_builder('reference/..', local)
        ....:         except SystemExit:
        ....:             traversal_rejected = True
        sage: type(answer).__name__, answer.name, traversal_rejected
        ('DocBuilder', 'fr/myreference', True)
    """
    reference = options.source_dir / 'en' / 'reference'
    reference_exists = (reference / 'index.rst').is_file()
    if name == 'reference_top' and reference_exists:
        return ReferenceTopBuilder('reference', options)
    if name in ('reference', 'en/reference') and reference_exists:
        return ReferenceBuilder(name, options)
    document = Path(name)
    if (len(document.parts) == 2 and document.parts[0] == 'reference'
            and document.parts[1].isidentifier()
            and (options.source_dir / 'en' / document / 'index.rst').is_file()):
        return ReferenceSubBuilder(name, options)
    if name.startswith('file='):
        path = name[5:]
        if path.endswith('.sage') or path.endswith('.pyx'):
            raise NotImplementedError('Building documentation for a single file only works for Python files.')
        return SingleFileBuilder(path, options)
    documents = get_all_documents(options.source_dir)
    for document in (Path(name), Path('en') / name):
        if document in documents:
            if document.name == 'website':
                return WebsiteBuilder(document.as_posix(), options)
            return DocBuilder(document.as_posix(), options)
    print("'%s' is not a recognized document. Type 'sage --docbuild -D' for a list" % name)
    print("of documents, or 'sage --docbuild --help' for more help.")
    sys.exit(1)


def get_documents(source: Path) -> list[str]:
    """
    Return the documents that the documentation builder accepts as
    command-line arguments, in the order in which they are built.

    These are the documents of :func:`get_all_documents`, named the way the
    command line names them: the English ones without their language prefix,
    and the reference manual - which :func:`get_all_documents` leaves out, its
    top level having a builder of its own - in front.

    EXAMPLES::

        sage: from sage.env import SAGE_DOC_SRC
        sage: from sage_docbuild.builders import get_documents
        sage: documents = get_documents(Path(SAGE_DOC_SRC))
        sage: documents[0]
        'reference'
        sage: 'tutorial' in documents
        True
        sage: 'en/tutorial' in documents
        False
    """
    reference_index = source / 'en' / 'reference' / 'index.rst'
    documents = ['reference'] if reference_index.is_file() else []
    for document in get_all_documents(source):
        if document.parts[0] == 'en':
            documents.append(document.relative_to('en').as_posix())
        else:
            documents.append(document.as_posix())
    return documents


def get_all_documents(source: Path) -> list[Path]:
    r"""
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

    A directory without the root source that Sphinx would build is not a
    document::

        sage: import tempfile
        sage: with tempfile.TemporaryDirectory() as directory:
        ....:     source = Path(directory)
        ....:     empty = source / 'en' / 'empty'
        ....:     complete = source / 'en' / 'complete'
        ....:     empty.mkdir(parents=True)
        ....:     complete.mkdir()
        ....:     _ = (complete / 'index.rst').write_text('Complete\n========\n')
        ....:     get_all_documents(source) == [Path('en/complete')]
        True
    """
    documents = []
    for lang in [path for path in source.iterdir() if path.is_dir()]:
        if not re.match('^[a-z][a-z]$', lang.name):
            # Skip non-language directories
            continue
        for document in lang.iterdir():
            if (document.name not in build_options.OMIT
                    and document.is_dir()
                    and (document / 'index.rst').is_file()):
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
