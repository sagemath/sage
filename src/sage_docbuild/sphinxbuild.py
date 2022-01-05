# -*- coding: utf-8 -*-
r"""
This is Sage's version of the sphinx-build script

We redirect stdout and stderr to our own logger, and remove some unwanted chatter.
"""
# ****************************************************************************
#       Copyright (C) 2013-2014 Volker Braun <vbraun.name@gmail.com>
#                     2013-2017 J. H. Palmieri <<palmieri@math.washington.edu>
#                     2013-2017 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#                     2014      Christopher Schwan <cschwan@students.uni-mainz.de>
#                     2014      Nicolas M. Thiéry <nthiery@users.sf.net>
#                     2015      Marc Mezzarobba <marc@mezzarobba.net>
#                     2015      André Apitzsch <andre.apitzsch@etit.tu-chemnitz.de>
#                     2018      Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

import os, sys, re, sphinx
import sphinx.cmd.build

# override the fancy multi-line formatting
def term_width_line(text):
    return text + '\n'

sphinx.util.console.term_width_line = term_width_line


class SageSphinxLogger(object):
    r"""
    This implements the file object interface to serve as
    ``sys.stdout``/``sys.stderr`` replacement.
    """
    ansi_color = re.compile(r'\x1b\[[0-9;]*m')
    ansi_reset = re.compile(r'\x1b\[39;49;00m')
    prefix_len = 9

    def __init__(self, stream, prefix):
        self._init_chatter()
        self._stream = stream
        self._color = stream.isatty()
        prefix = prefix[0:self.prefix_len]
        prefix = ('[{0:'+str(self.prefix_len)+'}]').format(prefix)
        self._is_stdout = (stream.fileno() == 1)
        self._is_stderr = (stream.fileno() == 2)
        if self._is_stdout:
            color = 'darkgreen'
        elif self._is_stderr:
            color = 'red'
        else:
            color = 'lightgray'
        self._prefix = sphinx.util.console.colorize(color, prefix)
        # When we see an error in the log, we store it here and raise it at the
        # end of the file (sometimes the lines following the error still
        # contain valuable information.)
        self._error = None

    def _init_chatter(self):
        # We drop any messages from the output that match these regular
        # expressions. These just bloat the output and do not contain any
        # information that we care about.
        self._useless_chatter = (
            re.compile('^$'),
            re.compile('^Running Sphinx v'),
            re.compile('^loading intersphinx inventory from '),
            re.compile('^loading pickled environment... done'),
            re.compile(r'^loading cross citations... done \([0-9]* citations\).'),
            re.compile('^Compiling a sub-document'),
            re.compile('^updating environment: 0 added, 0 changed, 0 removed'),
            re.compile('^looking for now-outdated files... none found'),
            re.compile(r'^building \[.*\]: targets for 0 source files that are out of date'),
            re.compile(r'^building \[.*\]: targets for 0 po files that are out of date'),
            re.compile(r'^building \[.*\]: targets for 0 mo files that are out of date'),
            re.compile('^pickling environment... done'),
            re.compile('^dumping object inventory... done'),
            # We still have "Build finished."
            re.compile('^build succeeded.'),
            re.compile('^checking consistency... done'),
            re.compile('^preparing documents... done'),
            re.compile('^copying extra files... done'),
            re.compile('^writing additional pages... search'),
            re.compile('^Writing js search indexes...writing additional pages... .*'),
            re.compile('^generating indices... .*'),
            re.compile('^dumping search index in .* ... done'),
            re.compile('^linking _static directory'),
            re.compile('^copying static files... done'),
            re.compile('^copying extra files... done'),
            re.compile(r'^loading translations \[.*\]... done'),
            re.compile('^Compiling the master document'),
            re.compile('^Saved pickle file: citations.pickle'),
            re.compile(r'^writing output... \[.*\] '),
            re.compile(r'^copying images... \[.*\] '),
            re.compile(r'^reading sources... \[.*\] '),
            re.compile('language "hu" not supported'),
            re.compile('^$'),
            re.compile('^WARNING:$'),
            )

        # We fail whenever a line starts with "WARNING:", however, we ignore
        # these warnings, as they are not relevant.
        self._ignored_warnings = (
            re.compile("WARNING: favicon file 'favicon.ico' does not exist"),
            re.compile('WARNING: html_static_path entry .* does not exist'),
            re.compile('WARNING: while setting up extension'),
            re.compile('WARNING: Any IDs not assiend for figure node'),
            re.compile('WARNING: .* is not referenced'),
            re.compile('WARNING: Build finished'),
            re.compile('WARNING: rST localisation for language .* not found')
        )
        # The warning "unknown config value 'multidoc_first_pass'..."
        # should only appear when building the documentation for a
        # single file (SingleFileBuilder from __init__.py), and it
        # needs to be ignored in that case. See #29651.
        self._ignored_warnings += (re.compile('WARNING: unknown config value \'multidoc_first_pass\''),)
        self._useless_chatter += self._ignored_warnings

        # replacements: pairs of regular expressions and their replacements,
        # to be applied to Sphinx output.
        self.replacements = [(re.compile('build succeeded, [0-9]+ warning[s]?.'),
                              'build succeeded.')]

        if 'inventory' in sys.argv:
            # When building the inventory, ignore warnings about missing
            # citations and the search index.
            ignored = (
                re.compile('WARNING: citation not found:'),
                re.compile("WARNING: search index couldn't be loaded, but not all documents will be built: the index will be incomplete.")
                )
            self._ignored_warnings += ignored
            self._useless_chatter += ignored

        # Regular expressions indicating a problem with docbuilding. Raise an
        # exception if any of these occur.
        self._error_patterns = (re.compile('Segmentation fault'),
                    re.compile('SEVERE'),
                    re.compile('ERROR'),
                    re.compile('^make.*Error'),
                    re.compile('Exception occurred'),
                    re.compile('Sphinx error'))

        # We want all warnings to actually be errors.
        # Exceptions:
        # - warnings upon building the LaTeX documentation
        # - undefined labels upon the first pass of the compilation: some
        #   cross links may legitimately not yet be resolvable at this point.
        if 'latex' not in sys.argv:
            self._error_patterns += (re.compile('WARNING:'),)
            if 'multidoc_first_pass=1' in sys.argv:
                ignore = (re.compile('WARNING: undefined label'),)
                self._ignored_warnings += ignore
                self._useless_chatter += ignore

    def _filter_out(self, line):
        if self._error is not None and self._is_stdout:
            # swallow non-errors after an error occurred
            return True
        line = re.sub(self.ansi_color, '', line)
        line = line.strip()
        for regex in self._useless_chatter:
            if regex.search(line) is not None:
                return True
        return False

    def _check_errors(self, line):
        r"""
        Search for errors in line.

        EXAMPLES::

            sage: from sys import stdout
            sage: from sage_docbuild.sphinxbuild import SageSphinxLogger
            sage: logger = SageSphinxLogger(stdout, "doctesting")
            sage: logger._log_line("Segmentation fault!\n") # indirect doctest
            [doctestin] Segmentation fault!
            sage: logger.raise_errors()
            Traceback (most recent call last):
            ...
            OSError: Segmentation fault!

        """
        if self._error is not None:
            return  # we already have found an error
        for error in self._error_patterns:
            if error.search(line) is not None:
                for ignored in self._ignored_warnings:
                    if ignored.search(line) is not None:
                        break
                else:
                    self._error = line
                    return

    def _log_line(self, line):
        r"""
        Write ``line`` to the output stream with some mangling.

        EXAMPLES::

            sage: from sys import stdout
            sage: from sage_docbuild.sphinxbuild import SageSphinxLogger
            sage: logger = SageSphinxLogger(stdout, "doctesting")
            sage: logger._log_line("building documentation…\n")
            [doctestin] building documentation…

        TESTS:

        Verify that :trac:`25160` has been resolved::

            sage: logger = SageSphinxLogger(stdout, "#25160")
            sage: import traceback
            sage: try:
            ....:     raise Exception("artificial exception")
            ....: except Exception:
            ....:     for line in traceback.format_exc().split('\n'):
            ....:         logger._log_line(line)
            [#25160   ] Traceback (most recent call last):
            [#25160   ]   File ...
            [#25160   ]     raise Exception("artificial exception")
            [#25160   ] Exception: artificial exception
        """
        skip_this_line = self._filter_out(line)
        self._check_errors(line)
        for (old, new) in self.replacements:
            line = old.sub(new, line)
        line = self._prefix + ' ' + line.rstrip() + '\n'
        if not self._color:
            line = self.ansi_color.sub('', line)
        if not skip_this_line:
            # sphinx does produce messages in the current locals which
            # could be non-ascii
            # see https://trac.sagemath.org/ticket/27706
            self._stream.write(line if isinstance(line, str) else line.encode('utf8'))
            self._stream.flush()

    def raise_errors(self):
        r"""
        Raise an exceptions if any errors have been found while parsing the
        Sphinx output.

        EXAMPLES::

            sage: from sys import stdout
            sage: from sage_docbuild.sphinxbuild import SageSphinxLogger
            sage: logger = SageSphinxLogger(stdout, "doctesting")
            sage: logger._log_line("This is a SEVERE error\n")
            [doctestin] This is a SEVERE error
            sage: logger.raise_errors()
            Traceback (most recent call last):
            ...
            OSError: This is a SEVERE error

        """
        if self._error is not None:
            raise OSError(self._error)

    _line_buffer = ''

    def _write(self, string):
        self._line_buffer += string
        lines = self._line_buffer.splitlines()
        for i, line in enumerate(lines):
            last = (i == len(lines)-1)
            if last and not self._line_buffer.endswith('\n'):
                self._line_buffer = line
                return
            self._log_line(line)
        self._line_buffer = ''


    # file object interface follows

    closed = False
    encoding = None
    mode = 'w'
    name = '<log>'
    newlines = None
    softspace = 0

    def isatty(self):
        return True

    def close(self):
        if self._line_buffer != '':
            self._log_line(self._line_buffer)
            self._line_buffer = ''

    def flush(self):
        self._stream.flush()

    def write(self, str):
        try:
            self._write(str)
        except OSError:
            raise
        except Exception:
            import traceback
            traceback.print_exc(file=self._stream)

    def writelines(self, sequence):
        for line in sequence:
            self.write(line)


def runsphinx():
    output_dir = sys.argv[-1]

    saved_stdout = sys.stdout
    saved_stderr = sys.stderr

    try:
        sys.stdout = SageSphinxLogger(sys.stdout, os.path.basename(output_dir))
        sys.stderr = SageSphinxLogger(sys.stderr, os.path.basename(output_dir))
        # Note that this call as of early 2018 leaks memory. So make sure that
        # you don't call runsphinx() several times in a row. (i.e., you want to
        # fork() somewhere before this call.)
        # We don't use subprocess here, as we don't want to re-initialize Sage
        # for every docbuild as this takes a while.
        sphinx.cmd.build.main(sys.argv[1:])
        sys.stderr.raise_errors()
        sys.stdout.raise_errors()
    finally:
        sys.stdout = saved_stdout
        sys.stderr = saved_stderr
        sys.stdout.flush()
        sys.stderr.flush()
