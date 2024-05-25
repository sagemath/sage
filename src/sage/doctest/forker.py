# sage_setup: distribution = sagemath-repl
"""
Processes for running doctests

This module controls the processes started by Sage that actually run
the doctests.

EXAMPLES:

The following examples are used in doctesting this file::

    sage: doctest_var = 42; doctest_var^2
    1764
    sage: R.<a> = ZZ[]
    sage: a + doctest_var
    a + 42

AUTHORS:

- David Roe (2012-03-27) -- initial version, based on Robert Bradshaw's code.

- Jeroen Demeyer (2013 and 2015) -- major improvements to forking and logging
"""

# ****************************************************************************
#       Copyright (C) 2012-2013 David Roe <roed.math@gmail.com>
#                     2012      Robert Bradshaw <robertwb@gmail.com>
#                     2012      William Stein <wstein@gmail.com>
#                     2013      R. Andrew Ohana
#                     2013-2018 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#                     2013-2020 John H. Palmieri
#                     2013-2017 Volker Braun
#                     2014      André Apitzsch
#                     2014      Darij Grinberg
#                     2016-2021 Frédéric Chapoton
#                     2017-2019 Erik M. Bray
#                     2018      Julian Rüth
#                     2020      Jonathan Kliem
#                     2020-2023 Matthias Koeppe
#                     2022      Markus Wageringel
#                     2022      Michael Orlitzky
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


import os
import sys
import time
import signal
import linecache
import hashlib
import warnings
import re
import errno
import doctest
import traceback
import tempfile
from collections import defaultdict
from dis import findlinestarts
from queue import Empty
import gc
import IPython.lib.pretty

import sage.misc.randstate as randstate
from sage.misc.timing import walltime
from .util import Timer, RecordingDict, count_noun
from .sources import DictAsObject
from .parsing import OriginalSource, reduce_hex
from sage.structure.sage_object import SageObject
from .parsing import SageOutputChecker, pre_hash, get_source, unparse_optional_tags
from sage.repl.user_globals import set_globals
from sage.cpython.atexit import restore_atexit
from sage.cpython.string import bytes_to_str, str_to_bytes

# With OS X, Python 3.8 defaults to use 'spawn' instead of 'fork' in
# multiprocessing, and Sage doctesting doesn't work with 'spawn'. See
# trac #27754.
if os.uname().sysname == 'Darwin':
    import multiprocessing
    multiprocessing.set_start_method('fork', force=True)

try:
    import multiprocessing
    from multiprocessing import Process
except ImportError:
    multiprocessing = None
    class Process:
        pass


def _sorted_dict_pprinter_factory(start, end):
    """
    Modified version of :func:`IPython.lib.pretty._dict_pprinter_factory`
    that sorts the keys of dictionaries for printing.

    EXAMPLES::

        sage: {2: 0, 1: 0} # indirect doctest
        {1: 0, 2: 0}
    """
    def inner(obj, p, cycle):
        if cycle:
            return p.text('{...}')
        step = len(start)
        p.begin_group(step, start)
        keys = obj.keys()
        keys = IPython.lib.pretty._sorted_for_pprint(keys)
        for idx, key in p._enumerate(keys):
            if idx:
                p.text(',')
                p.breakable()
            p.pretty(key)
            p.text(': ')
            p.pretty(obj[key])
        p.end_group(step, end)
    return inner


def init_sage(controller=None):
    """
    Import the Sage library.

    This function is called once at the beginning of a doctest run
    (rather than once for each file).  It imports the Sage library,
    sets DOCTEST_MODE to True, and invalidates any interfaces.

    EXAMPLES::

        sage: from sage.doctest.forker import init_sage
        sage: sage.doctest.DOCTEST_MODE = False
        sage: init_sage()
        sage: sage.doctest.DOCTEST_MODE
        True

    Check that pexpect interfaces are invalidated, but still work::

        sage: gap.eval("my_test_var := 42;")
        '42'
        sage: gap.eval("my_test_var;")
        '42'
        sage: init_sage()
        sage: gap('Group((1,2,3)(4,5), (3,4))')
        Group( [ (1,2,3)(4,5), (3,4) ] )
        sage: gap.eval("my_test_var;")
        Traceback (most recent call last):
        ...
        RuntimeError: Gap produced error output...

    Check that SymPy equation pretty printer is limited in doctest
    mode to default width (80 chars)::

        sage: # needs sympy
        sage: from sympy import sympify
        sage: from sympy.printing.pretty.pretty import PrettyPrinter
        sage: s = sympify('+x^'.join(str(i) for i in range(30)))
        sage: print(PrettyPrinter(settings={'wrap_line': True}).doprint(s))
         29    28    27    26    25    24    23    22    21    20    19    18    17...
        x   + x   + x   + x   + x   + x   + x   + x   + x   + x   + x   + x   + x...
        <BLANKLINE>
        ... 16    15    14    13    12    11    10    9    8    7    6    5    4    3...
        ...x   + x   + x   + x   + x   + x   + x   + x  + x  + x  + x  + x  + x  + x...
        <BLANKLINE>
        ...

    The displayhook sorts dictionary keys to simplify doctesting of
    dictionary output::

        sage: {'a':23, 'b':34, 'au':56, 'bbf':234, 'aaa':234}
        {'a': 23, 'aaa': 234, 'au': 56, 'b': 34, 'bbf': 234}
    """
    try:
        # We need to ensure that the Matplotlib font cache is built to
        # avoid spurious warnings (see Issue #20222).
        import matplotlib.font_manager
    except ImportError:
        # Do not require matplotlib for running doctests (Issue #25106).
        pass
    else:
        # Make sure that the agg backend is selected during doctesting.
        # This needs to be done before any other matplotlib calls.
        matplotlib.use('agg')

    # Do this once before forking off child processes running the tests.
    # This is more efficient because we only need to wait once for the
    # Sage imports.
    import sage.doctest
    sage.doctest.DOCTEST_MODE = True

    # Set the Python PRNG class to the Python 2 implementation for consistency
    # of 'random' test results that use it; see
    # https://github.com/sagemath/sage/issues/24508
    # We use the baked in copy of the random module for both Python 2 and 3
    # since, although the upstream copy is unlikely to change, this further
    # ensures consistency of test results
    import sage.misc.randstate
    from sage.cpython._py2_random import Random
    sage.misc.randstate.DEFAULT_PYTHON_RANDOM = Random

    # IPython's pretty printer sorts the repr of dicts by their keys by default
    # (or their keys' str() if they are not otherwise orderable).  However, it
    # disables this for CPython 3.6+ opting to instead display dicts' "natural"
    # insertion order, which is preserved in those versions).
    # However, this order is random in some instances.
    # Also modifications of code may affect the order.
    # So here we fore sorted dict printing.
    IPython.lib.pretty.for_type(dict, _sorted_dict_pprinter_factory('{', '}'))

    if controller is None:
        import sage.repl.ipython_kernel.all_jupyter
    else:
        controller.load_environment()

    try:
        from sage.interfaces.quit import invalidate_all
        invalidate_all()
    except ModuleNotFoundError:
        pass

    # Disable cysignals debug messages in doctests: this is needed to
    # make doctests pass when cysignals was built with debugging enabled
    from cysignals.signals import set_debug_level
    set_debug_level(0)

    # Use the rich output backend for doctest
    from sage.repl.rich_output import get_display_manager
    dm = get_display_manager()
    from sage.repl.rich_output.backend_doctest import BackendDoctest
    dm.switch_backend(BackendDoctest())

    # Switch on extra debugging
    from sage.structure.debug_options import debug
    debug.refine_category_hash_check = True

    # We import readline before forking, otherwise Pdb doesn't work
    # on OS X: https://github.com/sagemath/sage/issues/14289
    try:
        import readline
    except ModuleNotFoundError:
        # Do not require readline for running doctests (Issue #31160).
        pass

    try:
        import sympy
    except ImportError:
        # Do not require sympy for running doctests (Issue #25106).
        pass
    else:
        # Disable SymPy terminal width detection
        from sympy.printing.pretty.stringpict import stringPict
        stringPict.terminal_width = lambda self: 0


def showwarning_with_traceback(message, category, filename, lineno, file=None, line=None):
    r"""
    Displays a warning message with a traceback.

    INPUT: see :func:`warnings.showwarning`.

    OUTPUT: none

    EXAMPLES::

        sage: from sage.doctest.forker import showwarning_with_traceback
        sage: showwarning_with_traceback("bad stuff", UserWarning, "myfile.py", 0)
        doctest:warning...
          File "<doctest sage.doctest.forker.showwarning_with_traceback[1]>", line 1, in <module>
            showwarning_with_traceback("bad stuff", UserWarning, "myfile.py", Integer(0))
        :
        UserWarning: bad stuff
    """
    # Flush stdout to get predictable ordering of output and warnings
    sys.stdout.flush()

    # Get traceback to display in warning
    tb = traceback.extract_stack()
    tb = tb[:-1]  # Drop this stack frame for showwarning_with_traceback()
    for i, frame_summary in enumerate(tb):
        if frame_summary.filename.endswith('sage/doctest/forker.py') and frame_summary.name == 'compile_and_execute':
            tb = tb[i + 1:]
            break

    # Format warning
    lines = ["doctest:warning\n"]  # Match historical warning messages in doctests
    lines.extend(traceback.format_list(tb))
    lines.append(":\n")            # Match historical warning messages in doctests
    lines.extend(traceback.format_exception_only(category, category(message)))

    if file is None:
        file = sys.stderr
    try:
        file.writelines(lines)
        file.flush()
    except OSError:
        pass  # the file is invalid


class SageSpoofInOut(SageObject):
    r"""
    We replace the standard :class:`doctest._SpoofOut` for three reasons:

    - we need to divert the output of C programs that don't print
      through sys.stdout,
    - we want the ability to recover partial output from doctest
      processes that segfault.
    - we also redirect stdin (usually from /dev/null) during doctests.

    This class defines streams ``self.real_stdin``, ``self.real_stdout``
    and ``self.real_stderr`` which refer to the original streams.

    INPUT:

    - ``outfile`` -- (default: ``tempfile.TemporaryFile()``) a seekable open file
      object to which stdout and stderr should be redirected

    - ``infile`` -- (default: ``open(os.devnull)``) an open file object
      from which stdin should be redirected

    EXAMPLES::

        sage: import subprocess, tempfile
        sage: from sage.doctest.forker import SageSpoofInOut
        sage: O = tempfile.TemporaryFile()
        sage: S = SageSpoofInOut(O)
        sage: try:
        ....:     S.start_spoofing()
        ....:     print("hello world")
        ....: finally:
        ....:     S.stop_spoofing()
        ....:
        sage: S.getvalue()
        'hello world\n'
        sage: _ = O.seek(0)
        sage: S = SageSpoofInOut(outfile=sys.stdout, infile=O)
        sage: try:
        ....:     S.start_spoofing()
        ....:     _ = subprocess.check_call("cat")
        ....: finally:
        ....:     S.stop_spoofing()
        ....:
        hello world
        sage: O.close()
    """
    def __init__(self, outfile=None, infile=None):
        """
        Initialization.

        TESTS::

            sage: from tempfile import TemporaryFile
            sage: from sage.doctest.forker import SageSpoofInOut
            sage: with TemporaryFile() as outfile:
            ....:     with TemporaryFile() as infile:
            ....:         SageSpoofInOut(outfile, infile)
            <sage.doctest.forker.SageSpoofInOut object at ...>
        """
        if infile is None:
            self.infile = open(os.devnull)
            self._close_infile = True
        else:
            self.infile = infile
            self._close_infile = False
        if outfile is None:
            self.outfile = tempfile.TemporaryFile()
            self._close_outfile = True
        else:
            self.outfile = outfile
            self._close_outfile = False
        self.spoofing = False
        self.real_stdin = os.fdopen(os.dup(sys.stdin.fileno()), "r")
        self.real_stdout = os.fdopen(os.dup(sys.stdout.fileno()), "w")
        self.real_stderr = os.fdopen(os.dup(sys.stderr.fileno()), "w")
        self.position = 0

    def __del__(self):
        """
        Stop spoofing.

        TESTS::

            sage: from sage.doctest.forker import SageSpoofInOut
            sage: spoof = SageSpoofInOut()
            sage: spoof.start_spoofing()
            sage: print("Spoofed!")  # No output
            sage: del spoof
            sage: print("Not spoofed!")
            Not spoofed!
        """
        self.stop_spoofing()
        if self._close_infile:
            self.infile.close()
        if self._close_outfile:
            self.outfile.close()
        for stream in ('stdin', 'stdout', 'stderr'):
            getattr(self, 'real_' + stream).close()

    def start_spoofing(self):
        r"""
        Set stdin to read from ``self.infile`` and stdout to print to
        ``self.outfile``.

        EXAMPLES::

            sage: import os, tempfile
            sage: from sage.doctest.forker import SageSpoofInOut
            sage: O = tempfile.TemporaryFile()
            sage: S = SageSpoofInOut(O)
            sage: try:
            ....:     S.start_spoofing()
            ....:     print("this is not printed")
            ....: finally:
            ....:     S.stop_spoofing()
            ....:
            sage: S.getvalue()
            'this is not printed\n'
            sage: _ = O.seek(0)
            sage: S = SageSpoofInOut(infile=O)
            sage: try:
            ....:     S.start_spoofing()
            ....:     v = sys.stdin.read()
            ....: finally:
            ....:     S.stop_spoofing()
            ....:
            sage: v
            'this is not printed\n'

        We also catch non-Python output::

            sage: try:
            ....:     S.start_spoofing()
            ....:     retval = os.system('''echo "Hello there"\nif [ $? -eq 0 ]; then\necho "good"\nfi''')
            ....: finally:
            ....:     S.stop_spoofing()
            ....:
            sage: S.getvalue()
            'Hello there\ngood\n'
            sage: O.close()
        """
        if not self.spoofing:
            sys.stdout.flush()
            sys.stderr.flush()
            self.outfile.flush()
            os.dup2(self.infile.fileno(), sys.stdin.fileno())
            os.dup2(self.outfile.fileno(), sys.stdout.fileno())
            os.dup2(self.outfile.fileno(), sys.stderr.fileno())
            self.spoofing = True

    def stop_spoofing(self):
        """
        Reset stdin and stdout to their original values.

        EXAMPLES::

            sage: from sage.doctest.forker import SageSpoofInOut
            sage: S = SageSpoofInOut()
            sage: try:
            ....:     S.start_spoofing()
            ....:     print("this is not printed")
            ....: finally:
            ....:     S.stop_spoofing()
            ....:
            sage: print("this is now printed")
            this is now printed
        """
        if self.spoofing:
            sys.stdout.flush()
            sys.stderr.flush()
            self.real_stdout.flush()
            self.real_stderr.flush()
            os.dup2(self.real_stdin.fileno(), sys.stdin.fileno())
            os.dup2(self.real_stdout.fileno(), sys.stdout.fileno())
            os.dup2(self.real_stderr.fileno(), sys.stderr.fileno())
            self.spoofing = False

    def getvalue(self):
        r"""
        Get the value that has been printed to ``outfile`` since the
        last time this function was called.

        EXAMPLES::

            sage: from sage.doctest.forker import SageSpoofInOut
            sage: S = SageSpoofInOut()
            sage: try:
            ....:     S.start_spoofing()
            ....:     print("step 1")
            ....: finally:
            ....:     S.stop_spoofing()
            ....:
            sage: S.getvalue()
            'step 1\n'
            sage: try:
            ....:     S.start_spoofing()
            ....:     print("step 2")
            ....: finally:
            ....:     S.stop_spoofing()
            ....:
            sage: S.getvalue()
            'step 2\n'
        """
        sys.stdout.flush()
        self.outfile.seek(self.position)
        result = self.outfile.read()
        self.position = self.outfile.tell()
        if not result.endswith(b"\n"):
            result += b"\n"
        return bytes_to_str(result)


from collections import namedtuple
TestResults = namedtuple('TestResults', 'failed attempted')


class SageDocTestRunner(doctest.DocTestRunner):
    def __init__(self, *args, **kwds):
        """
        A customized version of DocTestRunner that tracks dependencies
        of doctests.

        INPUT:

        - ``stdout`` -- an open file to restore for debugging

        - ``checker`` -- ``None``, or an instance of
          :class:`doctest.OutputChecker`

        - ``verbose`` -- boolean, determines whether verbose printing
          is enabled

        - ``optionflags`` -- controls the comparison with the expected
          output.  See :mod:`testmod` for more information

        - ``baseline`` -- dictionary, the ``baseline_stats`` value

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: DTR
            <sage.doctest.forker.SageDocTestRunner object at ...>
        """
        O = kwds.pop('outtmpfile', None)
        self.msgfile = kwds.pop('msgfile', None)
        self.options = kwds.pop('sage_options')
        self.baseline = kwds.pop('baseline', {})
        doctest.DocTestRunner.__init__(self, *args, **kwds)
        self._fakeout = SageSpoofInOut(O)
        if self.msgfile is None:
            self.msgfile = self._fakeout.real_stdout
        self.history = []
        self.references = []
        self.setters = defaultdict(dict)
        self.running_global_digest = hashlib.md5()
        self.total_walltime_skips = 0
        self.total_performed_tests = 0
        self.total_walltime = 0

    def _run(self, test, compileflags, out):
        """
        This function replaces :meth:`doctest.DocTestRunner.__run`.

        It changes the following behavior:

        - We call :meth:`SageDocTestRunner.execute` rather than just
          exec

        - We don't truncate _fakeout after each example since we want
          the output file to be readable by the calling
          :class:`SageWorker`.

        Since it needs to be able to read stdout, it should be called
        while spoofing using :class:`SageSpoofInOut`.

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = sage.doctest.forker.__file__
            sage: FDS = FileDocTestSource(filename, DD)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: DTR.run(doctests[0], clear_globs=False) # indirect doctest
            TestResults(failed=0, attempted=4)

        TESTS:

        Check that :issue:`26038` is fixed::

            sage: a = 1
            ....: b = 2
            Traceback (most recent call last):
            ...
            SyntaxError: multiple statements found while compiling a single statement
            sage: a = 1
            ....: @syntax error
            Traceback (most recent call last):
            ...
            SyntaxError: multiple statements found while compiling a single statement
        """
        # Ensure that injecting globals works as expected in doctests
        set_globals(test.globs)

        # Keep track of the number of failures and tries.
        failures = tries = walltime_skips = 0
        quiet = False

        # Save the option flags (since option directives can be used
        # to modify them).
        original_optionflags = self.optionflags

        SUCCESS, FAILURE, BOOM = range(3)  # `outcome` state

        check = self._checker.check_output

        # Process each example.
        for examplenum, example in enumerate(test.examples):
            if failures:
                # If exitfirst is set, abort immediately after a
                # failure.
                if self.options.exitfirst:
                    break

                # If REPORT_ONLY_FIRST_FAILURE is set, then suppress
                # reporting after the first failure (but continue
                # running the tests).
                quiet |= (self.optionflags & doctest.REPORT_ONLY_FIRST_FAILURE)

            # Merge in the example's options.
            self.optionflags = original_optionflags
            if example.options:
                for (optionflag, val) in example.options.items():
                    if val:
                        self.optionflags |= optionflag
                    else:
                        self.optionflags &= ~optionflag

            # Skip this test if we exceeded our --short budget of walltime for
            # this doctest
            if self.options.target_walltime != -1 and self.total_walltime >= self.options.target_walltime:
                walltime_skips += 1
                self.optionflags |= doctest.SKIP

            # If 'SKIP' is set, then skip this example.
            if self.optionflags & doctest.SKIP:
                continue

            # Record that we started this example.
            tries += 1

            # We print the example we're running for easier debugging
            # if this file times out or crashes.
            with OriginalSource(example):
                print("sage: " + example.source[:-1] + " ## line %s ##" % (test.lineno + example.lineno + 1))
            # Update the position so that result comparison works
            self._fakeout.getvalue()
            if not quiet:
                self.report_start(out, test, example)

            # Flush files before running the example, so we know for
            # sure that everything is reported properly if the test
            # crashes.
            sys.stdout.flush()
            sys.stderr.flush()
            self.msgfile.flush()

            # Use a special filename for compile(), so we can retrieve
            # the source code during interactive debugging (see
            # __patched_linecache_getlines).
            filename = '<doctest %s[%d]>' % (test.name, examplenum)

            # Run the example in the given context (globs), and record
            # any exception that gets raised. But for SystemExit, we
            # simply propagate the exception.
            exception = None

            def compiler(example):
                # Compile mode "single" is meant for running a single
                # statement like on the Python command line. It implies
                # in particular that the resulting value will be printed.
                code = compile(example.source, filename, "single",
                               compileflags, 1)

                # Python 2 ignores everything after the first complete
                # statement in the source code. To verify that we really
                # have just a single statement and nothing more, we also
                # compile in "exec" mode and verify that the line
                # numbers are the same.
                execcode = compile(example.source, filename, "exec",
                                   compileflags, 1)

                # findlinestarts() returns pairs (index, lineno) where
                # "index" is the index in the bytecode where the line
                # number changes to "lineno".
                linenumbers1 = {lineno for (index, lineno)
                                in findlinestarts(code)}
                linenumbers2 = {lineno for (index, lineno)
                                in findlinestarts(execcode)}
                if linenumbers1 != linenumbers2:
                    raise SyntaxError("doctest is not a single statement")

                return code

            if not self.options.gc:
                pass
            elif self.options.gc > 0:
                if gc.isenabled():
                    gc.collect()
            elif self.options.gc < 0:
                gc.disable()

            from cysignals.signals import SignalError
            try:
                # Don't blink!  This is where the user's code gets run.
                self.compile_and_execute(example, compiler, test.globs)
            except (SignalError, SystemExit):
                # Tests can be killed by signals in unexpected places.
                raise
            except BaseException:
                exception = sys.exc_info()
            finally:
                if self.debugger is not None:
                    self.debugger.set_continue()  # ==== Example Finished ====
            check_timer = Timer().start()
            got = self._fakeout.getvalue()

            outcome = FAILURE   # guilty until proved innocent or insane

            probed_tags = getattr(example, 'probed_tags', False)

            # If the example executed without raising any exceptions,
            # verify its output.
            if exception is None:
                if check(example.want, got, self.optionflags):
                    if probed_tags and probed_tags is not True:
                        example.warnings.append(
                            f"The tag '{unparse_optional_tags(probed_tags)}' "
                            f"may no longer be needed; these features are not present, "
                            f"but we ran the doctest anyway as requested by --probe, "
                            f"and it succeeded.")
                    outcome = SUCCESS

            # The example raised an exception: check if it was expected.
            else:
                exc_msg = traceback.format_exception_only(*exception[:2])[-1]

                if example.exc_msg is not None:
                    # On Python 3 the exception repr often includes the
                    # exception's full module name (for non-builtin
                    # exceptions), whereas on Python 2 does not, so we
                    # normalize Python 3 exceptions to match tests written to
                    # Python 2
                    # See https://github.com/sagemath/sage/issues/24271
                    exc_cls = exception[0]
                    exc_name = exc_cls.__name__
                    if exc_cls.__module__:
                        exc_fullname = (exc_cls.__module__ + '.' +
                                        exc_cls.__qualname__)
                    else:
                        exc_fullname = exc_cls.__qualname__

                    if (example.exc_msg.startswith(exc_name) and
                            exc_msg.startswith(exc_fullname)):
                        exc_msg = exc_msg.replace(exc_fullname, exc_name, 1)

                if not quiet:
                    got += doctest._exception_traceback(exception)

                # If `example.exc_msg` is None, then we weren't expecting
                # an exception.
                if example.exc_msg is None:
                    outcome = BOOM

                # We expected an exception: see whether it matches.
                elif check(example.exc_msg, exc_msg, self.optionflags):
                    if probed_tags and probed_tags is not True:
                        example.warnings.append(
                            f"The tag '{unparse_optional_tags(example.probed_tags)}' "
                            f"may no longer be needed; these features are not present, "
                            f"but we ran the doctest anyway as requested by --probe, "
                            f"and it succeeded (raised the expected exception).")
                    outcome = SUCCESS

                # Another chance if they didn't care about the detail.
                elif self.optionflags & doctest.IGNORE_EXCEPTION_DETAIL:
                    m1 = re.match(r'(?:[^:]*\.)?([^:]*:)', example.exc_msg)
                    m2 = re.match(r'(?:[^:]*\.)?([^:]*:)', exc_msg)
                    if m1 and m2 and check(m1.group(1), m2.group(1),
                                           self.optionflags):
                        if probed_tags and probed_tags is not True:
                            example.warnings.append(
                                f"The tag '{unparse_optional_tags(example.probed_tags)}' "
                                f"may no longer be needed; these features are not present, "
                                f"but we ran the doctest anyway as requested by --probe, "
                                f"and it succeeded (raised an exception as expected).")
                        outcome = SUCCESS

            check_timer.stop()
            self.total_walltime += example.walltime + check_timer.walltime

            # Report the outcome.
            if example.warnings:
                for warning in example.warnings:
                    out(self._failure_header(test, example, f'Warning: {warning}'))
            if outcome is SUCCESS:
                if self.options.warn_long > 0 and example.cputime + check_timer.cputime > self.options.warn_long:
                    self.report_overtime(out, test, example, got,
                                         check_timer=check_timer)
                elif example.warnings:
                    pass
                elif not quiet:
                    self.report_success(out, test, example, got,
                                        check_timer=check_timer)
            elif probed_tags:
                pass
            elif outcome is FAILURE:
                if not quiet:
                    self.report_failure(out, test, example, got, test.globs)
                failures += 1
            elif outcome is BOOM:
                if not quiet:
                    self.report_unexpected_exception(out, test, example,
                                                     exception)
                failures += 1
            else:
                assert False, ("unknown outcome", outcome)

        # Restore the option flags (in case they were modified)
        self.optionflags = original_optionflags

        # Record and return the number of failures and tries.
        self._DocTestRunner__record_outcome(test, failures, tries)
        self.total_walltime_skips += walltime_skips
        self.total_performed_tests += tries
        return TestResults(failures, tries)

    def run(self, test, compileflags=0, out=None, clear_globs=True):
        """
        Run the examples in a given doctest.

        This function replaces :class:`doctest.DocTestRunner.run`
        since it needs to handle spoofing. It also leaves the display
        hook in place.

        INPUT:

        - ``test`` -- an instance of :class:`doctest.DocTest`

        - ``compileflags`` -- integer (default: 0) the set of compiler flags
          used to execute examples (passed in to the :func:`compile`)

        - ``out`` -- a function for writing the output (defaults to
          :func:`sys.stdout.write`)

        - ``clear_globs`` -- boolean (default: ``True``); whether to clear
          the namespace after running this doctest

        OUTPUT:

        - ``f`` -- integer, the number of examples that failed

        - ``t`` -- the number of examples tried

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), verbose=False, sage_options=DD,
            ....:                         optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = sage.doctest.forker.__file__
            sage: FDS = FileDocTestSource(filename, DD)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: DTR.run(doctests[0], clear_globs=False)
            TestResults(failed=0, attempted=4)
        """
        self.setters = defaultdict(dict)
        randstate.set_random_seed(self.options.random_seed)
        warnings.showwarning = showwarning_with_traceback
        self.running_doctest_digest = hashlib.md5()
        self.test = test
        # We use this slightly modified version of Pdb because it
        # interacts better with the doctesting framework (like allowing
        # doctests for sys.settrace()). Since we already have output
        # spoofing in place, there is no need for redirection.
        if self.options.debug:
            self.debugger = doctest._OutputRedirectingPdb(sys.stdout)
            self.debugger.reset()
        else:
            self.debugger = None
        self.save_linecache_getlines = linecache.getlines
        linecache.getlines = self._DocTestRunner__patched_linecache_getlines
        if out is None:
            def out(s):
                self.msgfile.write(s)
                self.msgfile.flush()

        self._fakeout.start_spoofing()
        # If self.options.initial is set, we show only the first failure in each doctest block.
        self.no_failure_yet = True
        try:
            return self._run(test, compileflags, out)
        finally:
            self._fakeout.stop_spoofing()
            linecache.getlines = self.save_linecache_getlines
            if clear_globs:
                test.globs.clear()

    def summarize(self, verbose=None):
        """
        Print results of testing to ``self.msgfile`` and return number
        of failures and tests run.

        INPUT:

        - ``verbose`` -- whether to print lots of stuff

        OUTPUT:

        - returns ``(f, t)``, a :class:`doctest.TestResults` instance
          giving the number of failures and the total number of tests
          run.

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: DTR._name2ft['sage.doctest.forker'] = (1,120)
            sage: results = DTR.summarize()
            **********************************************************************
            1 item had failures:
                1 of 120 in sage.doctest.forker
            sage: results
            TestResults(failed=1, attempted=120)
        """
        if verbose is None:
            verbose = self._verbose
        m = self.msgfile
        notests = []
        passed = []
        failed = []
        totalt = totalf = 0
        for x in self._name2ft.items():
            name, (f, t) = x
            assert f <= t
            totalt += t
            totalf += f
            if not t:
                notests.append(name)
            elif not f:
                passed.append((name, t))
            else:
                failed.append(x)
        if verbose:
            if notests:
                print(count_noun(len(notests), "item"), "had no tests:", file=m)
                notests.sort()
                for thing in notests:
                    print("    %s" % thing, file=m)
            if passed:
                print(count_noun(len(passed), "item"), "passed all tests:", file=m)
                passed.sort()
                for thing, count in passed:
                    print(" %s in %s" % (count_noun(count, "test", pad_number=3, pad_noun=True), thing), file=m)
        if failed:
            print(self.DIVIDER, file=m)
            print(count_noun(len(failed), "item"), "had failures:", file=m)
            failed.sort()
            for thing, (f, t) in failed:
                print(" %3d of %3d in %s" % (f, t, thing), file=m)
        if verbose:
            print(count_noun(totalt, "test") + " in " + count_noun(len(self._name2ft), "item") + ".", file=m)
            print("%s passed and %s failed." % (totalt - totalf, totalf), file=m)
            if totalf:
                print("***Test Failed***", file=m)
            else:
                print("Test passed.", file=m)
        m.flush()
        return doctest.TestResults(totalf, totalt)

    def update_digests(self, example):
        """
        Update global and doctest digests.

        Sage's doctest runner tracks the state of doctests so that
        their dependencies are known.  For example, in the following
        two lines ::

            sage: R.<x> = ZZ[]
            sage: f = x^2 + 1

        it records that the second line depends on the first since the
        first INSERTS ``x`` into the global namespace and the second
        line RETRIEVES ``x`` from the global namespace.

        This function updates the hashes that record these
        dependencies.

        INPUT:

        - ``example`` -- a :class:`doctest.Example` instance

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os, hashlib
            sage: DTR = SageDocTestRunner(SageOutputChecker(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = sage.doctest.forker.__file__
            sage: FDS = FileDocTestSource(filename, DD)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: DTR.running_global_digest.hexdigest()
            'd41d8cd98f00b204e9800998ecf8427e'
            sage: DTR.running_doctest_digest = hashlib.md5()
            sage: ex = doctests[0].examples[0]; ex.predecessors = None
            sage: DTR.update_digests(ex)
            sage: DTR.running_global_digest.hexdigest()
            '3cb44104292c3a3ab4da3112ce5dc35c'
        """
        s = str_to_bytes(pre_hash(get_source(example)), 'utf-8')
        self.running_global_digest.update(s)
        self.running_doctest_digest.update(s)
        if example.predecessors is not None:
            digest = hashlib.md5(s)
            gen = (e.running_state for e in example.predecessors)
            digest.update(str_to_bytes(reduce_hex(gen), 'ascii'))
            example.running_state = digest.hexdigest()

    def compile_and_execute(self, example, compiler, globs):
        """
        Run the given example, recording dependencies.

        Rather than using a basic dictionary, Sage's doctest runner
        uses a :class:`sage.doctest.util.RecordingDict`, which records
        every time a value is set or retrieved.  Executing the given
        code with this recording dictionary as the namespace allows
        Sage to track dependencies between doctest lines.  For
        example, in the following two lines ::

            sage: R.<x> = ZZ[]
            sage: f = x^2 + 1

        the recording dictionary records that the second line depends
        on the first since the first INSERTS ``x`` into the global
        namespace and the second line RETRIEVES ``x`` from the global
        namespace.

        INPUT:

        - ``example`` -- a :class:`doctest.Example` instance

        - ``compiler`` -- a callable that, applied to example,
          produces a code object

        - ``globs`` -- dictionary in which to execute the code

        OUTPUT: the output of the compiled code snippet

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.util import RecordingDict
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os, hashlib
            sage: DTR = SageDocTestRunner(SageOutputChecker(), verbose=False, sage_options=DD,
            ....:           optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: DTR.running_doctest_digest = hashlib.md5()
            sage: filename = sage.doctest.forker.__file__
            sage: FDS = FileDocTestSource(filename, DD)
            sage: globs = RecordingDict(globals())
            sage: 'doctest_var' in globs
            False
            sage: doctests, extras = FDS.create_doctests(globs)
            sage: ex0 = doctests[0].examples[0]
            sage: flags = 32768 if sys.version_info.minor < 8 else 524288
            sage: def compiler(ex):
            ....:     return compile(ex.source, '<doctest sage.doctest.forker[0]>',
            ....:                    'single', flags, 1)
            sage: DTR.compile_and_execute(ex0, compiler, globs)
            1764
            sage: globs['doctest_var']
            42
            sage: globs.set
            {'doctest_var'}
            sage: globs.got
            {'Integer'}

        Now we can execute some more doctests to see the dependencies. ::

            sage: ex1 = doctests[0].examples[1]
            sage: def compiler(ex):
            ....:     return compile(ex.source, '<doctest sage.doctest.forker[1]>',
            ....:                    'single', flags, 1)
            sage: DTR.compile_and_execute(ex1, compiler, globs)
            sage: sorted(list(globs.set))
            ['R', 'a']
            sage: globs.got
            {'ZZ'}
            sage: ex1.predecessors
            []

        ::

            sage: ex2 = doctests[0].examples[2]
            sage: def compiler(ex):
            ....:     return compile(ex.source, '<doctest sage.doctest.forker[2]>',
            ....:                    'single', flags, 1)
            sage: DTR.compile_and_execute(ex2, compiler, globs)
            a + 42
            sage: list(globs.set)
            []
            sage: sorted(list(globs.got))
            ['a', 'doctest_var']
            sage: set(ex2.predecessors) == set([ex0,ex1])
            True
        """
        if isinstance(globs, RecordingDict):
            globs.start()
        example.sequence_number = len(self.history)
        if not hasattr(example, 'warnings'):
            example.warnings = []
        self.history.append(example)
        timer = Timer().start()
        try:
            compiled = compiler(example)
            timer.start()    # reset timer
            exec(compiled, globs)
        finally:
            timer.stop().annotate(example)
            if isinstance(globs, RecordingDict):
                example.predecessors = []
                for name in globs.got:
                    setters_dict = self.setters.get(name)  # setter_optional_tags -> setter
                    if setters_dict:
                        was_set = False
                        for setter_optional_tags, setter in setters_dict.items():
                            if setter_optional_tags.issubset(example.optional_tags):  # was set in a less constrained doctest
                                was_set = True
                                example.predecessors.append(setter)
                        if not was_set:
                            if example.probed_tags:
                                # Probing confusion.
                                # Do not issue the "was set only in doctest marked" warning;
                                # and also do not issue the "may no longer be needed" notice
                                example.probed_tags = True
                            else:
                                f_setter_optional_tags = "; ".join("'"
                                                                   + unparse_optional_tags(setter_optional_tags)
                                                                   + "'"
                                                                   for setter_optional_tags in setters_dict)
                                example.warnings.append(f"Variable '{name}' referenced here "
                                                        f"was set only in doctest marked {f_setter_optional_tags}")
                for name in globs.set:
                    self.setters[name][example.optional_tags] = example
            else:
                example.predecessors = None
            self.update_digests(example)
            example.total_state = self.running_global_digest.hexdigest()
            example.doctest_state = self.running_doctest_digest.hexdigest()

    def _failure_header(self, test, example, message='Failed example:'):
        """
        We strip out ``sage:`` prompts, so we override
        :meth:`doctest.DocTestRunner._failure_header` for better
        reporting.

        INPUT:

        - ``test`` -- a :class:`doctest.DocTest` instance

        - ``example`` -- a :class:`doctest.Example` instance in ``test``

        OUTPUT: string used for reporting that the given example failed

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = sage.doctest.forker.__file__
            sage: FDS = FileDocTestSource(filename, DD)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: ex = doctests[0].examples[0]
            sage: print(DTR._failure_header(doctests[0], ex))
            **********************************************************************
            File ".../sage/doctest/forker.py", line 12, in sage.doctest.forker
            Failed example:
                doctest_var = 42; doctest_var^2
            <BLANKLINE>

        Without the source swapping::

            sage: import doctest
            sage: print(doctest.DocTestRunner._failure_header(DTR, doctests[0], ex))
            **********************************************************************
            File ".../sage/doctest/forker.py", line 12, in sage.doctest.forker
            Failed example:
                doctest_var = Integer(42); doctest_var**Integer(2)
            <BLANKLINE>

        The ``'Failed example:'`` message can be customized::

            sage: print(DTR._failure_header(doctests[0], ex, message='Hello there!'))
            **********************************************************************
            File ".../sage/doctest/forker.py", line 12, in sage.doctest.forker
            Hello there!
                doctest_var = 42; doctest_var^2
            <BLANKLINE>
        """
        out = [self.DIVIDER]
        with OriginalSource(example):
            if self.options.format == 'sage':
                if test.filename:
                    if test.lineno is not None and example.lineno is not None:
                        lineno = test.lineno + example.lineno + 1
                    else:
                        lineno = '?'
                    out.append('File "%s", line %s, in %s' %
                               (test.filename, lineno, test.name))
                else:
                    out.append('Line %s, in %s' % (example.lineno + 1, test.name))
                out.append(message)
            elif self.options.format == 'github':
                # https://docs.github.com/en/actions/using-workflows/workflow-commands-for-github-actions#using-workflow-commands-to-access-toolkit-functions
                if message.startswith('Warning: '):
                    command = f'::warning title={message}'
                    message = message[len('Warning: '):]
                elif self.baseline.get('failed', False):
                    command = f'::notice title={message}'
                    message += ' [failed in baseline]'
                else:
                    command = f'::error title={message}'
                if extra := getattr(example, 'extra', None):
                    message += f': {extra}'
                if test.filename:
                    command += f',file={test.filename}'
                    if test.lineno is not None and example.lineno is not None:
                        lineno = test.lineno + example.lineno + 1
                        command += f',line={lineno}'
                    lineno = None
                else:
                    command += f',line={example.lineno + 1}'
                #
                # Urlencoding trick for multi-line annotations
                # https://github.com/actions/starter-workflows/issues/68#issuecomment-581479448
                #
                # This only affects the display in the workflow Summary, after clicking "Show more";
                # the message needs to be long enough so that "Show more" becomes available.
                # https://github.com/actions/toolkit/issues/193#issuecomment-1867084340
                #
                # Unfortunately, this trick does not make the annotations in the diff view multi-line.
                #
                if '\n' in message:
                    message = message.replace('\n', '%0A')
                    # The actual threshold for "Show more" to appear depends on the window size.
                    show_more_threshold = 500
                    if (pad := show_more_threshold - len(message)) > 0:
                        message += ' ' * pad
                command += f'::{message}'
                out.append(command)
            else:
                raise ValueError(f'unknown format option: {self.options.format}')
            source = example.source
            out.append(doctest._indent(source))
            return '\n'.join(out)

    def report_start(self, out, test, example):
        """
        Called when an example starts.

        INPUT:

        - ``out`` -- a function for printing

        - ``test`` -- a :class:`doctest.DocTest` instance

        - ``example`` -- a :class:`doctest.Example` instance in ``test``

        OUTPUT: prints a report to ``out``

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), verbose=True, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = sage.doctest.forker.__file__
            sage: FDS = FileDocTestSource(filename, DD)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: ex = doctests[0].examples[0]
            sage: DTR.report_start(sys.stdout.write, doctests[0], ex)
            Trying (line 12):    doctest_var = 42; doctest_var^2
            Expecting:
                1764
        """
        # We completely replace doctest.DocTestRunner.report_start so that we can include line numbers
        with OriginalSource(example):
            if self._verbose:
                start_txt = ('Trying (line %s):' % (test.lineno + example.lineno + 1)
                             + doctest._indent(example.source))
                if example.want:
                    start_txt += 'Expecting:\n' + doctest._indent(example.want)
                else:
                    start_txt += 'Expecting nothing\n'
                out(start_txt)

    def report_success(self, out, test, example, got, *, check_timer=None):
        """
        Called when an example succeeds.

        INPUT:

        - ``out`` -- a function for printing

        - ``test`` -- a :class:`doctest.DocTest` instance

        - ``example`` -- a :class:`doctest.Example` instance in ``test``

        - ``got`` -- string; the result of running ``example``

        - ``check_timer`` -- a :class:`sage.doctest.util.Timer` (default:
          ``None``) that measures the time spent checking whether or not
          the output was correct

        OUTPUT: prints a report to ``out``; if in debugging mode, starts an
        IPython prompt at the point of the failure

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: from sage.doctest.util import Timer
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), verbose=True, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = sage.doctest.forker.__file__
            sage: FDS = FileDocTestSource(filename, DD)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: ex = doctests[0].examples[0]
            sage: ex.cputime = 1.01
            sage: ex.walltime = 1.12
            sage: check = Timer()
            sage: check.cputime = 2.14
            sage: check.walltime = 2.71
            sage: DTR.report_success(sys.stdout.write, doctests[0], ex, '1764',
            ....:                    check_timer=check)
            ok [3.83s wall]
        """
        # We completely replace doctest.DocTestRunner.report_success
        # so that we can include time taken for the test
        if self._verbose:
            out("ok [%.2fs wall]\n" %
                (example.walltime + check_timer.walltime))

    def report_failure(self, out, test, example, got, globs):
        r"""
        Called when a doctest fails.

        INPUT:

        - ``out`` -- a function for printing

        - ``test`` -- a :class:`doctest.DocTest` instance

        - ``example`` -- a :class:`doctest.Example` instance in ``test``

        - ``got`` -- string, the result of running ``example``

        - ``globs`` -- dictionary of globals, used if in debugging mode

        OUTPUT: prints a report to ``out``

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), verbose=True, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = sage.doctest.forker.__file__
            sage: FDS = FileDocTestSource(filename, DD)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: ex = doctests[0].examples[0]
            sage: DTR.no_failure_yet = True
            sage: DTR.report_failure(sys.stdout.write, doctests[0], ex, 'BAD ANSWER\n', {})
            **********************************************************************
            File ".../sage/doctest/forker.py", line 12, in sage.doctest.forker
            Failed example:
                doctest_var = 42; doctest_var^2
            Expected:
                1764
            Got:
                BAD ANSWER

        If debugging is turned on this function starts an IPython
        prompt when a test returns an incorrect answer::

            sage: sage0.quit()
            sage: _ = sage0.eval("import doctest, sys, os, multiprocessing, subprocess")
            sage: _ = sage0.eval("from sage.doctest.parsing import SageOutputChecker")
            sage: _ = sage0.eval("import sage.doctest.forker as sdf")
            sage: _ = sage0.eval("from sage.doctest.control import DocTestDefaults")
            sage: _ = sage0.eval("DD = DocTestDefaults(debug=True)")
            sage: _ = sage0.eval("ex1 = doctest.Example('a = 17', '')")
            sage: _ = sage0.eval("ex2 = doctest.Example('2*a', '1')")
            sage: _ = sage0.eval("DT = doctest.DocTest([ex1,ex2], globals(), 'doubling', None, 0, None)")
            sage: _ = sage0.eval("DTR = sdf.SageDocTestRunner(SageOutputChecker(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)")
            sage: print(sage0.eval("sdf.init_sage(); DTR.run(DT, clear_globs=False)")) # indirect doctest
            **********************************************************************
            Line 1, in doubling
            Failed example:
                2*a
            Expected:
                1
            Got:
                34
            **********************************************************************
            Previously executed commands:
            sage: sage0._expect.expect('sage: ')   # sage0 just mis-identified the output as prompt, synchronize
            0
            sage: sage0.eval("a")
            '...17'
            sage: sage0.eval("quit")
            'Returning to doctests...TestResults(failed=1, attempted=2)'
        """
        if not self.options.initial or self.no_failure_yet:
            self.no_failure_yet = False
            example.extra = f'Got: {got}'
            returnval = doctest.DocTestRunner.report_failure(self, out, test, example, got)
            if self.options.debug:
                self._fakeout.stop_spoofing()
                restore_tcpgrp = None
                try:
                    if os.isatty(0):
                        # In order to read from the terminal, we need
                        # to make the current process group the
                        # foreground group.
                        restore_tcpgrp = os.tcgetpgrp(0)
                        signal.signal(signal.SIGTTIN, signal.SIG_IGN)
                        signal.signal(signal.SIGTTOU, signal.SIG_IGN)
                        os.tcsetpgrp(0, os.getpgrp())
                    print("*" * 70)
                    print("Previously executed commands:")
                    for ex in test.examples:
                        if ex is example:
                            break
                        if hasattr(ex, 'sage_source'):
                            src = '    sage: ' + ex.sage_source
                        else:
                            src = '    sage: ' + ex.source
                        if src[-1] == '\n':
                            src = src[:-1]
                        src = src.replace('\n', '\n    ....: ')
                        print(src)
                        if ex.want:
                            print(doctest._indent(ex.want[:-1]))
                    from sage.repl.configuration import sage_ipython_config
                    from IPython.terminal.embed import InteractiveShellEmbed
                    cfg = sage_ipython_config.default()
                    # Currently this doesn't work: prompts only work in pty
                    # We keep simple_prompt=True, prompts will be "In [0]:"
                    # cfg.InteractiveShell.prompts_class = DebugPrompts
                    # cfg.InteractiveShell.simple_prompt = False
                    shell = InteractiveShellEmbed(config=cfg, banner1='', user_ns=dict(globs))
                    shell(header='', stack_depth=2)
                except KeyboardInterrupt:
                    # Assume this is a *real* interrupt. We need to
                    # escalate this to the master doctesting process.
                    if not self.options.serial and multiprocessing:
                        os.kill(os.getppid(), signal.SIGINT)
                    raise
                finally:
                    # Restore the foreground process group.
                    if restore_tcpgrp is not None:
                        os.tcsetpgrp(0, restore_tcpgrp)
                        signal.signal(signal.SIGTTIN, signal.SIG_DFL)
                        signal.signal(signal.SIGTTOU, signal.SIG_DFL)
                    print("Returning to doctests...")
                    self._fakeout.start_spoofing()
            return returnval

    def report_overtime(self, out, test, example, got, *, check_timer=None):
        r"""
        Called when the ``warn_long`` option flag is set and a doctest
        runs longer than the specified time.

        INPUT:

        - ``out`` -- a function for printing

        - ``test`` -- a :class:`doctest.DocTest` instance

        - ``example`` -- a :class:`doctest.Example` instance in ``test``

        - ``got`` -- string; the result of running ``example``

        - ``check_timer`` -- a :class:`sage.doctest.util.Timer` (default:
          ``None``) that measures the time spent checking whether or not
          the output was correct

        OUTPUT: prints a report to ``out``

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: from sage.doctest.util import Timer
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), verbose=True, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = sage.doctest.forker.__file__
            sage: FDS = FileDocTestSource(filename, DD)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: ex = doctests[0].examples[0]
            sage: ex.cputime = 1.23
            sage: ex.walltime = 2.50
            sage: check = Timer()
            sage: check.cputime = 2.34
            sage: check.walltime = 3.12
            sage: DTR.report_overtime(sys.stdout.write, doctests[0], ex, 'BAD ANSWER\n', check_timer=check)
            **********************************************************************
            File ".../sage/doctest/forker.py", line 12, in sage.doctest.forker
            Warning: slow doctest:
                doctest_var = 42; doctest_var^2
            Test ran for 1.23s cpu, 2.50s wall
            Check ran for 2.34s cpu, 3.12s wall
        """
        out(self._failure_header(test, example, 'Warning: slow doctest:') +
            ('Test ran for %.2fs cpu, %.2fs wall\nCheck ran for %.2fs cpu, %.2fs wall\n'
             % (example.cputime,
                example.walltime,
                check_timer.cputime,
                check_timer.walltime)))

    def report_unexpected_exception(self, out, test, example, exc_info):
        r"""
        Called when a doctest raises an exception that's not matched by the expected output.

        If debugging has been turned on, starts an interactive debugger.

        INPUT:

        - ``out`` -- a function for printing

        - ``test`` -- a :class:`doctest.DocTest` instance

        - ``example`` -- a :class:`doctest.Example` instance in ``test``

        - ``exc_info`` -- the result of ``sys.exc_info()``

        OUTPUT: prints a report to ``out``

        - if in debugging mode, starts PDB with the given traceback

        EXAMPLES::

            sage: from sage.interfaces.sage0 import sage0
            sage: sage0.quit()
            sage: _ = sage0.eval("import doctest, sys, os, multiprocessing, subprocess")
            sage: _ = sage0.eval("from sage.doctest.parsing import SageOutputChecker")
            sage: _ = sage0.eval("import sage.doctest.forker as sdf")
            sage: _ = sage0.eval("from sage.doctest.control import DocTestDefaults")
            sage: _ = sage0.eval("DD = DocTestDefaults(debug=True)")
            sage: _ = sage0.eval("ex = doctest.Example('E = EllipticCurve([0,0]); E', 'A singular Elliptic Curve')")
            sage: _ = sage0.eval("DT = doctest.DocTest([ex], globals(), 'singular_curve', None, 0, None)")
            sage: _ = sage0.eval("DTR = sdf.SageDocTestRunner(SageOutputChecker(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)")
            sage: old_prompt = sage0._prompt
            sage: sage0._prompt = r"\(Pdb\) "
            sage: sage0.eval("DTR.run(DT, clear_globs=False)") # indirect doctest
            '... ArithmeticError(self._equation_string() + " defines a singular curve")'
            sage: sage0.eval("l")
            '...if self.discriminant() == 0:...raise ArithmeticError...'
            sage: sage0.eval("u")
            '...-> super().__init__(R, data, category=category)'
            sage: sage0.eval("u")
            '...EllipticCurve_field.__init__(self, K, ainvs)'
            sage: sage0.eval("p ainvs")
            '(0, 0, 0, 0, 0)'
            sage: sage0._prompt = old_prompt
            sage: sage0.eval("quit")
            'TestResults(failed=1, attempted=1)'
        """
        if not self.options.initial or self.no_failure_yet:
            self.no_failure_yet = False

            example.extra = "Exception raised:\n" + "".join(traceback.format_exception(*exc_info))
            returnval = doctest.DocTestRunner.report_unexpected_exception(self, out, test, example, exc_info)
            if self.options.debug:
                self._fakeout.stop_spoofing()
                restore_tcpgrp = None
                try:
                    if os.isatty(0):
                        # In order to read from the terminal, we need
                        # to make the current process group the
                        # foreground group.
                        restore_tcpgrp = os.tcgetpgrp(0)
                        signal.signal(signal.SIGTTIN, signal.SIG_IGN)
                        signal.signal(signal.SIGTTOU, signal.SIG_IGN)
                        os.tcsetpgrp(0, os.getpgrp())

                    exc_type, exc_val, exc_tb = exc_info
                    if exc_tb is None:
                        raise RuntimeError(
                            "could not start the debugger for an unexpected "
                            "exception, probably due to an unhandled error "
                            "in a C extension module")
                    self.debugger.reset()
                    self.debugger.interaction(None, exc_tb)
                except KeyboardInterrupt:
                    # Assume this is a *real* interrupt. We need to
                    # escalate this to the master doctesting process.
                    if not self.options.serial and multiprocessing:
                        os.kill(os.getppid(), signal.SIGINT)
                    raise
                finally:
                    # Restore the foreground process group.
                    if restore_tcpgrp is not None:
                        os.tcsetpgrp(0, restore_tcpgrp)
                        signal.signal(signal.SIGTTIN, signal.SIG_DFL)
                        signal.signal(signal.SIGTTOU, signal.SIG_DFL)
                    self._fakeout.start_spoofing()
            return returnval

    def update_results(self, D):
        """
        When returning results we pick out the results of interest
        since many attributes are not pickleable.

        INPUT:

        - ``D`` -- dictionary to update with cputime and walltime

        OUTPUT: the number of failures (or ``False`` if there is no failure
        attribute)

        EXAMPLES::

            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.sources import FileDocTestSource, DictAsObject
            sage: from sage.doctest.control import DocTestDefaults; DD = DocTestDefaults()
            sage: import doctest, sys, os
            sage: DTR = SageDocTestRunner(SageOutputChecker(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: filename = sage.doctest.forker.__file__
            sage: FDS = FileDocTestSource(filename, DD)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: from sage.doctest.util import Timer
            sage: T = Timer().start()
            sage: DTR.run(doctests[0])
            TestResults(failed=0, attempted=4)
            sage: T.stop().annotate(DTR)
            sage: D = DictAsObject({'cputime': [], 'walltime': [], 'err': None})
            sage: DTR.update_results(D)
            0
            sage: sorted(list(D.items()))
            [('cputime', [...]), ('err', None), ('failures', 0), ('tests', 4),
             ('walltime', [...]), ('walltime_skips', 0)]
        """
        for key in ["cputime", "walltime"]:
            if key not in D:
                D[key] = []
            if hasattr(self, key):
                D[key].append(self.__dict__[key])
        D['tests'] = self.total_performed_tests
        D['walltime_skips'] = self.total_walltime_skips
        if hasattr(self, 'failures'):
            D['failures'] = self.failures
            return self.failures
        else:
            return False


def dummy_handler(sig, frame):
    """
    Dummy signal handler for SIGCHLD (just to ensure the signal
    isn't ignored).

    TESTS::

        sage: import signal
        sage: from sage.doctest.forker import dummy_handler
        sage: _ = signal.signal(signal.SIGUSR1, dummy_handler)
        sage: os.kill(os.getpid(), signal.SIGUSR1)
        sage: signal.signal(signal.SIGUSR1, signal.SIG_DFL)
        <function dummy_handler at ...>
    """
    pass


class DocTestDispatcher(SageObject):
    """
    Create parallel :class:`DocTestWorker` processes and dispatches
    doctesting tasks.
    """
    def __init__(self, controller):
        """
        INPUT:

        - ``controller`` -- a :class:`sage.doctest.control.DocTestController` instance

        EXAMPLES::

            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.doctest.forker import DocTestDispatcher
            sage: DocTestDispatcher(DocTestController(DocTestDefaults(), []))
            <sage.doctest.forker.DocTestDispatcher object at ...>
        """
        self.controller = controller
        init_sage(controller)

    def serial_dispatch(self):
        """
        Run the doctests from the controller's specified sources in series.

        There is no graceful handling for signals, no possibility of
        interrupting tests and no timeout.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.doctest.forker import DocTestDispatcher
            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.util import Timer
            sage: import os
            sage: homset = os.path.join(SAGE_SRC, 'sage', 'rings', 'homset.py')
            sage: ideal = os.path.join(SAGE_SRC, 'sage', 'rings', 'ideal.py')
            sage: DC = DocTestController(DocTestDefaults(), [homset, ideal])
            sage: DC.expand_files_into_sources()
            sage: DD = DocTestDispatcher(DC)
            sage: DR = DocTestReporter(DC)
            sage: DC.reporter = DR
            sage: DC.dispatcher = DD
            sage: DC.timer = Timer().start()
            sage: DD.serial_dispatch()
            sage -t .../rings/homset.py
                [... tests, ...s wall]
            sage -t .../rings/ideal.py
                [... tests, ...s wall]
        """
        for source in self.controller.sources:
            heading = self.controller.reporter.report_head(source)
            baseline = self.controller.source_baseline(source)
            if not self.controller.options.only_errors:
                self.controller.log(heading)

            with tempfile.TemporaryFile() as outtmpfile:
                result = DocTestTask(source)(self.controller.options,
                                             outtmpfile, self.controller.logger,
                                             baseline=baseline)
                outtmpfile.seek(0)
                output = bytes_to_str(outtmpfile.read())

            self.controller.reporter.report(source, False, 0, result, output)
            if self.controller.options.exitfirst and result[1].failures:
                break

    def parallel_dispatch(self):
        r"""
        Run the doctests from the controller's specified sources in parallel.

        This creates :class:`DocTestWorker` subprocesses, while the master
        process checks for timeouts and collects and displays the results.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.doctest.forker import DocTestDispatcher
            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.util import Timer
            sage: import os
            sage: crem = os.path.join(SAGE_SRC, 'sage', 'databases', 'cremona.py')
            sage: bigo = os.path.join(SAGE_SRC, 'sage', 'rings', 'big_oh.py')
            sage: DC = DocTestController(DocTestDefaults(), [crem, bigo])
            sage: DC.expand_files_into_sources()
            sage: DD = DocTestDispatcher(DC)
            sage: DR = DocTestReporter(DC)
            sage: DC.reporter = DR
            sage: DC.dispatcher = DD
            sage: DC.timer = Timer().start()
            sage: DD.parallel_dispatch()
            sage -t .../databases/cremona.py
                [... tests, ...s wall]
            sage -t .../rings/big_oh.py
                [... tests, ...s wall]

        If the ``exitfirst=True`` option is given, the results for a failing
        module will be immediately printed and any other ongoing tests
        canceled::

            sage: from tempfile import NamedTemporaryFile as NTF
            sage: with NTF(suffix='.py', mode='w+t') as f1, \
            ....:      NTF(suffix='.py', mode='w+t') as f2:
            ....:     _ = f1.write("'''\nsage: import time; time.sleep(60)\n'''")
            ....:     f1.flush()
            ....:     _ = f2.write("'''\nsage: True\nFalse\n'''")
            ....:     f2.flush()
            ....:     DC = DocTestController(DocTestDefaults(exitfirst=True,
            ....:                                            nthreads=2),
            ....:                            [f1.name, f2.name])
            ....:     DC.expand_files_into_sources()
            ....:     DD = DocTestDispatcher(DC)
            ....:     DR = DocTestReporter(DC)
            ....:     DC.reporter = DR
            ....:     DC.dispatcher = DD
            ....:     DC.timer = Timer().start()
            ....:     DD.parallel_dispatch()
            sage -t ...
            **********************************************************************
            File "...", line 2, in ...
            Failed example:
                True
            Expected:
                False
            Got:
                True
            **********************************************************************
            1 item had failures:
               1 of   1 in ...
                [1 test, 1 failure, ...s wall]
            Killing test ...
        """
        opt = self.controller.options

        job_client = None
        try:
            from gnumake_tokenpool import JobClient, NoJobServer
        except ImportError:
            pass
        else:
            try:
                job_client = JobClient(use_cysignals=True)
            except NoJobServer:
                pass

        source_iter = iter(self.controller.sources)

        # If timeout was 0, simply set a very long time
        if opt.timeout <= 0:
            opt.timeout = 2**60
        # Timeout we give a process to die (after it received a SIGQUIT
        # signal). If it doesn't exit by itself in this many seconds, we
        # SIGKILL it. This is 5% of doctest timeout, with a maximum of
        # 10 minutes and a minimum of 60 seconds.
        die_timeout = opt.timeout * 0.05
        if die_timeout > 600:
            die_timeout = 600
        elif die_timeout < 60:
            die_timeout = 60
        # allow override via cmdline option
        if opt.die_timeout >= 0:
            die_timeout = opt.die_timeout

        # If we think that we can not finish running all tests until
        # target_endtime, we skip individual tests. (Only enabled with
        # --short.)
        if opt.target_walltime == -1:
            target_endtime = None
        else:
            target_endtime = time.time() + opt.target_walltime
        pending_tests = len(self.controller.sources)

        # List of alive DocTestWorkers (child processes). Workers which
        # are done but whose messages have not been read are also
        # considered alive.
        workers = []

        # List of DocTestWorkers which have finished running but
        # whose results have not been reported yet.
        finished = []

        # If exitfirst is set and we got a failure.
        abort_now = False

        # One particular worker that we are "following": we report the
        # messages while it's running. For other workers, we report the
        # messages if there is no followed worker.
        follow = None

        # Install signal handler for SIGCHLD
        signal.signal(signal.SIGCHLD, dummy_handler)

        # Logger
        log = self.controller.log

        from cysignals.pselect import PSelecter
        try:
            # Block SIGCHLD and SIGINT except during the pselect() call
            with PSelecter([signal.SIGCHLD, signal.SIGINT]) as sel:
                # Function to execute in the child process which exits
                # this "with" statement (which restores the signal mask)
                # and resets to SIGCHLD handler to default.
                # Since multiprocessing.Process is implemented using
                # fork(), signals would otherwise remain blocked in the
                # child process.
                def sel_exit():
                    signal.signal(signal.SIGCHLD, signal.SIG_DFL)
                    sel.__exit__(None, None, None)

                while True:
                    # To avoid calling time.time() all the time while
                    # checking for timeouts, we call it here, once per
                    # loop. It's not a problem if this isn't very
                    # precise, doctest timeouts don't need millisecond
                    # precision.
                    now = time.time()

                    # If there were any substantial changes in the state
                    # (new worker started or finished worker reported),
                    # restart this while loop instead of calling pselect().
                    # This ensures internal consistency and a reasonably
                    # accurate value for "now".
                    restart = False

                    # Process all workers. Check for timeouts on active
                    # workers and move finished/crashed workers to the
                    # "finished" list.
                    # Create a new list "new_workers" containing the active
                    # workers (to avoid updating "workers" in place).
                    new_workers = []
                    for w in workers:
                        if w.rmessages is not None or w.is_alive():
                            if now >= w.deadline:
                                # Timeout => (try to) kill the process
                                # group (which normally includes
                                # grandchildren) and close the message
                                # pipe.
                                # We don't report the timeout yet, we wait
                                # until the process has actually died.
                                w.kill()
                                w.deadline = now + die_timeout
                            if not w.is_alive():
                                # Worker is done but we haven't read all
                                # messages (possibly a grandchild still
                                # has the messages pipe open).
                                # Adjust deadline to read all messages:
                                newdeadline = now + die_timeout
                                w.deadline = min(w.deadline, newdeadline)
                            new_workers.append(w)
                        else:
                            # Save the result and output of the worker
                            # and close the associated file descriptors.
                            # It is important to do this now. If we
                            # would leave them open until we call
                            # report(), parallel testing can easily fail
                            # with a "Too many open files" error.
                            w.save_result_output()
                            # In python3 multiprocessing.Process also
                            # opens a pipe internally, which has to be
                            # closed here, as well.
                            # But afterwards, exitcode and pid are
                            # no longer available.
                            w.copied_exitcode = w.exitcode
                            w.copied_pid = w.pid
                            w.close()
                            finished.append(w)
                            if job_client:
                                job_client.release()

                    workers = new_workers

                    # Similarly, process finished workers.
                    new_finished = []
                    for w in finished:
                        if opt.exitfirst and w.result[1].failures:
                            abort_now = True
                        elif follow is not None and follow is not w:
                            # We are following a different worker, so
                            # we cannot report now.
                            new_finished.append(w)
                            continue

                        # Report the completion of this worker
                        log(w.messages, end="")
                        self.controller.reporter.report(
                            w.source,
                            w.killed,
                            w.copied_exitcode,
                            w.result,
                            w.output,
                            pid=w.copied_pid)

                        pending_tests -= 1

                        restart = True
                        follow = None

                    finished = new_finished

                    if abort_now:
                        break

                    # Start new workers if possible
                    while (source_iter is not None and len(workers) < opt.nthreads
                           and (not job_client or job_client.acquire())):
                        try:
                            source = next(source_iter)
                        except StopIteration:
                            source_iter = None
                            if job_client:
                                job_client.release()
                        else:
                            # Start a new worker.
                            import copy
                            worker_options = copy.copy(opt)
                            baseline = self.controller.source_baseline(source)
                            if target_endtime is not None:
                                worker_options.target_walltime = (target_endtime - now) / (max(1, pending_tests / opt.nthreads))
                            w = DocTestWorker(source, options=worker_options, funclist=[sel_exit], baseline=baseline)
                            heading = self.controller.reporter.report_head(w.source)
                            if not self.controller.options.only_errors:
                                w.messages = heading + "\n"
                            # Store length of heading to detect if the
                            # worker has something interesting to report.
                            w.heading_len = len(w.messages)
                            w.start()  # This might take some time
                            w.deadline = time.time() + opt.timeout
                            workers.append(w)
                            restart = True

                    # Recompute state if needed
                    if restart:
                        continue

                    # We are finished if there are no DocTestWorkers left
                    if len(workers) == 0:
                        # If there are no active workers, we should have
                        # reported all finished workers.
                        assert len(finished) == 0
                        break

                    # The master pselect() call
                    rlist = [w.rmessages for w in workers if w.rmessages is not None]
                    tmout = min(w.deadline for w in workers) - now
                    tmout = min(tmout, 5)
                    rlist, _, _, _ = sel.pselect(rlist, timeout=tmout)

                    # Read messages
                    for w in workers:
                        if w.rmessages is not None and w.rmessages in rlist:
                            w.read_messages()

                    # Find a worker to follow: if there is only one worker,
                    # always follow it. Otherwise, take the worker with
                    # the earliest deadline of all workers whose
                    # messages are more than just the heading.
                    if follow is None:
                        if len(workers) == 1:
                            follow = workers[0]
                        else:
                            for w in workers:
                                if len(w.messages) > w.heading_len:
                                    if follow is None or w.deadline < follow.deadline:
                                        follow = w

                    # Write messages of followed worker
                    if follow is not None:
                        log(follow.messages, end="")
                        follow.messages = ""
        finally:
            # Restore SIGCHLD handler (which is to ignore the signal)
            signal.signal(signal.SIGCHLD, signal.SIG_DFL)

            # Kill all remaining workers (in case we got interrupted)
            for w in workers:
                if w.kill():
                    log("Killing test %s" % w.source.printpath)
            # Fork a child process with the specific purpose of
            # killing the remaining workers.
            if len(workers) > 0 and os.fork() == 0:
                # Block these signals
                with PSelecter([signal.SIGQUIT, signal.SIGINT]):
                    try:
                        from time import sleep
                        sleep(die_timeout)
                        for w in workers:
                            w.kill()
                            if job_client:
                                job_client.release()
                    finally:
                        os._exit(0)

            # Hack to ensure multiprocessing leaves these processes
            # alone (in particular, it doesn't wait for them when we
            # exit).
            p = multiprocessing.process
            assert hasattr(p, '_children')
            p._children = set()

    def dispatch(self):
        """
        Run the doctests for the controller's specified sources,
        by calling :meth:`parallel_dispatch` or :meth:`serial_dispatch`
        according to the ``--serial`` option.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.doctest.forker import DocTestDispatcher
            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.util import Timer
            sage: import os
            sage: freehom = os.path.join(SAGE_SRC, 'sage', 'modules', 'free_module_homspace.py')
            sage: bigo = os.path.join(SAGE_SRC, 'sage', 'rings', 'big_oh.py')
            sage: DC = DocTestController(DocTestDefaults(), [freehom, bigo])
            sage: DC.expand_files_into_sources()
            sage: DD = DocTestDispatcher(DC)
            sage: DR = DocTestReporter(DC)
            sage: DC.reporter = DR
            sage: DC.dispatcher = DD
            sage: DC.timer = Timer().start()
            sage: DD.dispatch()
            sage -t .../sage/modules/free_module_homspace.py
                [... tests, ...s wall]
            sage -t .../sage/rings/big_oh.py
                [... tests, ...s wall]
        """
        if self.controller.options.serial or not multiprocessing:
            self.serial_dispatch()
        else:
            self.parallel_dispatch()


class DocTestWorker(Process):
    """
    The DocTestWorker process runs one :class:`DocTestTask` for a given
    source. It returns messages about doctest failures (or all tests if
    verbose doctesting) through a pipe and returns results through a
    ``multiprocessing.Queue`` instance (both these are created in the
    :meth:`start` method).

    It runs the task in its own process-group, such that killing the
    process group kills this process together with its child processes.

    The class has additional methods and attributes for bookkeeping
    by the master process. Except in :meth:`run`, nothing from this
    class should be accessed by the child process.

    INPUT:

    - ``source`` -- a :class:`DocTestSource` instance

    - ``options`` -- an object representing doctest options

    - ``funclist`` -- list of callables to be called at the start of
      the child process

    - ``baseline`` -- dictionary, the ``baseline_stats`` value

    EXAMPLES::

        sage: from sage.doctest.forker import DocTestWorker, DocTestTask
        sage: from sage.doctest.sources import FileDocTestSource
        sage: from sage.doctest.reporting import DocTestReporter
        sage: from sage.doctest.control import DocTestController, DocTestDefaults
        sage: filename = sage.doctest.util.__file__
        sage: DD = DocTestDefaults()
        sage: FDS = FileDocTestSource(filename, DD)
        sage: W = DocTestWorker(FDS, DD)
        sage: W.start()
        sage: DC = DocTestController(DD, filename)
        sage: reporter = DocTestReporter(DC)
        sage: W.join()  # Wait for worker to finish
        sage: result = W.result_queue.get()
        sage: reporter.report(FDS, False, W.exitcode, result, "")
            [... tests, ...s wall]
    """
    def __init__(self, source, options, funclist=[], baseline=None):
        """
        Initialization.

        TESTS::

            sage: run_doctests(sage.rings.big_oh) # indirect doctest
            Running doctests with ID ...
            Doctesting 1 file.
            sage -t .../sage/rings/big_oh.py
                [... tests, ...s wall]
            ----------------------------------------------------------------------
            All tests passed!
            ----------------------------------------------------------------------
            Total time for all tests: ... seconds
                cpu time: ... seconds
                cumulative wall time: ... seconds
            Features detected...
        """
        Process.__init__(self)

        self.source = source
        self.options = options
        self.funclist = funclist
        self.baseline = baseline

        # Open pipe for messages. These are raw file descriptors,
        # not Python file objects!
        self.rmessages, self.wmessages = os.pipe()

        # Create Queue for the result. Since we're running only one
        # doctest, this "queue" will contain only 1 element.
        self.result_queue = multiprocessing.Queue(1)

        # Temporary file for stdout/stderr of the child process.
        # Normally, this isn't used in the master process except to
        # debug timeouts/crashes.
        self.outtmpfile = tempfile.NamedTemporaryFile(delete=False)

        # Create string for the master process to store the messages
        # (usually these are the doctest failures) of the child.
        # These messages are read through the pipe created above.
        self.messages = ""

        # Has this worker been killed (because of a time out)?
        self.killed = False

    def run(self):
        """
        Run the :class:`DocTestTask` under its own PGID.

        TESTS::

            sage: run_doctests(sage.symbolic.units)  # indirect doctest                 # needs sage.symbolic
            Running doctests with ID ...
            Doctesting 1 file.
            sage -t .../sage/symbolic/units.py
                [... tests, ...s wall]
            ----------------------------------------------------------------------
            All tests passed!
            ----------------------------------------------------------------------
            Total time for all tests: ... seconds
                cpu time: ... seconds
                cumulative wall time: ... seconds
            Features detected...
        """
        os.setpgid(os.getpid(), os.getpid())

        # Run functions
        for f in self.funclist:
            f()

        # Write one byte to the pipe to signal to the master process
        # that we have started properly.
        os.write(self.wmessages, b"X")

        task = DocTestTask(self.source)

        # Ensure the Python stdin is the actual stdin
        # (multiprocessing redirects this).
        # We will do a more proper redirect of stdin in SageSpoofInOut.
        try:
            sys.stdin = os.fdopen(0, "r")
        except OSError:
            # We failed to open stdin for reading, this might happen
            # for example when running under "nohup" (Issue #14307).
            # Simply redirect stdin from /dev/null and try again.
            with open(os.devnull) as f:
                os.dup2(f.fileno(), 0)
            sys.stdin = os.fdopen(0, "r")

        # Close the reading end of the pipe (only the master should
        # read from the pipe) and open the writing end.
        os.close(self.rmessages)
        msgpipe = os.fdopen(self.wmessages, "w")
        try:
            task(self.options, self.outtmpfile, msgpipe, self.result_queue,
                 baseline=self.baseline)
        finally:
            msgpipe.close()
            self.outtmpfile.close()

    def start(self):
        """
        Start the worker and close the writing end of the message pipe.

        TESTS::

            sage: from sage.doctest.forker import DocTestWorker, DocTestTask
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: filename = sage.doctest.util.__file__
            sage: DD = DocTestDefaults()
            sage: FDS = FileDocTestSource(filename, DD)
            sage: W = DocTestWorker(FDS, DD)
            sage: W.start()
            sage: try:
            ....:     os.fstat(W.wmessages)
            ....: except OSError:
            ....:     print("Write end of pipe successfully closed")
            Write end of pipe successfully closed
            sage: W.join()  # Wait for worker to finish
        """
        super().start()

        # Close the writing end of the pipe (only the child should
        # write to the pipe).
        os.close(self.wmessages)

        # Read one byte from the pipe as a sign that the child process
        # has properly started (to avoid race conditions). In particular,
        # it will have its process group changed.
        os.read(self.rmessages, 1)

    def read_messages(self):
        """
        In the master process, read from the pipe and store the data
        read in the ``messages`` attribute.

        .. NOTE::

            This function may need to be called multiple times in
            order to read all of the messages.

        EXAMPLES::

            sage: from sage.doctest.forker import DocTestWorker, DocTestTask
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: filename = sage.doctest.util.__file__
            sage: DD = DocTestDefaults(verbose=True,nthreads=2)
            sage: FDS = FileDocTestSource(filename, DD)
            sage: W = DocTestWorker(FDS, DD)
            sage: W.start()
            sage: while W.rmessages is not None:
            ....:     W.read_messages()
            sage: W.join()
            sage: len(W.messages) > 0
            True
        """
        # It's absolutely important to execute only one read() system
        # call, more might block. Assuming that we used pselect()
        # correctly, one read() will not block.
        if self.rmessages is not None:
            s = os.read(self.rmessages, 4096)
            self.messages += bytes_to_str(s)
            if len(s) == 0:  # EOF
                os.close(self.rmessages)
                self.rmessages = None

    def save_result_output(self):
        """
        Annotate ``self`` with ``self.result`` (the result read through
        the ``result_queue`` and with ``self.output``, the complete
        contents of ``self.outtmpfile``. Then close the Queue and
        ``self.outtmpfile``.

        EXAMPLES::

            sage: from sage.doctest.forker import DocTestWorker, DocTestTask
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: filename = sage.doctest.util.__file__
            sage: DD = DocTestDefaults()
            sage: FDS = FileDocTestSource(filename, DD)
            sage: W = DocTestWorker(FDS, DD)
            sage: W.start()
            sage: W.join()
            sage: W.save_result_output()
            sage: sorted(W.result[1].keys())
            ['cputime', 'err', 'failures', 'optionals', 'tests', 'walltime', 'walltime_skips']
            sage: len(W.output) > 0
            True

        .. NOTE::

            This method is called from the parent process, not from the
            subprocess.
        """
        try:
            self.result = self.result_queue.get(block=False)
        except Empty:
            self.result = (0, DictAsObject({'err': 'noresult'}))
        del self.result_queue

        self.outtmpfile.seek(0)
        self.output = bytes_to_str(self.outtmpfile.read())
        self.outtmpfile.close()
        try:
            # Now it is safe to delete the outtmpfile; we manage this manually
            # so that the file does not get deleted via TemporaryFile.__del__
            # in the worker process
            os.unlink(self.outtmpfile.name)
        except OSError as exc:
            if exc.errno != errno.ENOENT:
                raise

        del self.outtmpfile

    def kill(self):
        """
        Kill this worker.  Return ``True`` if the signal(s) are sent
        successfully or ``False`` if the worker process no longer exists.

        This method is only called if there is something wrong with the
        worker. Under normal circumstances, the worker is supposed to
        exit by itself after finishing.

        The first time this is called, use ``SIGQUIT``. This will trigger
        the cysignals ``SIGQUIT`` handler and try to print an enhanced
        traceback.

        Subsequent times, use ``SIGKILL``.  Also close the message pipe
        if it was still open.

        EXAMPLES::

            sage: import time
            sage: from sage.doctest.forker import DocTestWorker, DocTestTask
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','tests','99seconds.rst')
            sage: DD = DocTestDefaults()
            sage: FDS = FileDocTestSource(filename, DD)

        We set up the worker to start by blocking ``SIGQUIT``, such that
        killing will fail initially::

            sage: from cysignals.pselect import PSelecter
            sage: import signal
            sage: def block_hup():
            ....:     # We never __exit__()
            ....:     PSelecter([signal.SIGQUIT]).__enter__()
            sage: W = DocTestWorker(FDS, DD, [block_hup])
            sage: W.start()
            sage: W.killed
            False
            sage: W.kill()
            True
            sage: W.killed
            True
            sage: time.sleep(float(0.2))  # Worker doesn't die
            sage: W.kill()         # Worker dies now
            True
            sage: time.sleep(1)
            sage: W.is_alive()
            False
        """

        if self.rmessages is not None:
            os.close(self.rmessages)
            self.rmessages = None

        try:
            if not self.killed:
                self.killed = True
                os.killpg(self.pid, signal.SIGQUIT)
            else:
                os.killpg(self.pid, signal.SIGKILL)
        except OSError as exc:
            # Handle a race condition where the process has exited on
            # its own by the time we get here, and ESRCH is returned
            # indicating no processes in the specified process group
            if exc.errno != errno.ESRCH:
                raise

            return False

        return True


class DocTestTask:
    """
    This class encapsulates the tests from a single source.

    This class does not insulate from problems in the source
    (e.g. entering an infinite loop or causing a segfault), that has to
    be dealt with at a higher level.

    INPUT:

    - ``source`` -- a :class:`sage.doctest.sources.DocTestSource` instance

    - ``verbose`` -- boolean, controls reporting of progress by :class:`doctest.DocTestRunner`

    EXAMPLES::

        sage: from sage.doctest.forker import DocTestTask
        sage: from sage.doctest.sources import FileDocTestSource
        sage: from sage.doctest.control import DocTestDefaults, DocTestController
        sage: import os
        sage: filename = sage.doctest.sources.__file__
        sage: DD = DocTestDefaults()
        sage: FDS = FileDocTestSource(filename, DD)
        sage: DTT = DocTestTask(FDS)
        sage: DC = DocTestController(DD,[filename])
        sage: ntests, results = DTT(options=DD)
        sage: ntests >= 300 or ntests
        True
        sage: sorted(results.keys())
        ['cputime', 'err', 'failures', 'optionals', 'tests', 'walltime', 'walltime_skips']
    """

    def __init__(self, source):
        """
        Initialization.

        TESTS::

            sage: from sage.doctest.forker import DocTestTask
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults
            sage: import os
            sage: filename = sage.doctest.sources.__file__
            sage: FDS = FileDocTestSource(filename, DocTestDefaults())
            sage: DocTestTask(FDS)
            <sage.doctest.forker.DocTestTask object at ...>
        """
        self.source = source

    def __call__(self, options, outtmpfile=None, msgfile=None, result_queue=None, *,
                 baseline=None):
        """
        Calling the task does the actual work of running the doctests.

        INPUT:

        - ``options`` -- an object representing doctest options

        - ``outtmpfile`` -- a seekable file that's used by the doctest
          runner to redirect stdout and stderr of the doctests

        - ``msgfile`` -- a file or pipe to send doctest messages about
          doctest failures (or all tests in verbose mode)

        - ``result_queue`` -- an instance of :class:`multiprocessing.Queue`
          to store the doctest result. For testing, this can also be ``None``

        - ``baseline`` -- dictionary, the ``baseline_stats`` value

        OUTPUT:

        - ``(doctests, result_dict)`` where ``doctests`` is the number of
          doctests and ``result_dict`` is a dictionary annotated with
          timings and error information.

        - Also put ``(doctests, result_dict)`` onto the ``result_queue``
          if the latter isn't None.

        EXAMPLES::

            sage: from sage.doctest.forker import DocTestTask
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: import os
            sage: filename = sage.doctest.parsing.__file__
            sage: DD = DocTestDefaults()
            sage: FDS = FileDocTestSource(filename, DD)
            sage: DTT = DocTestTask(FDS)
            sage: DC = DocTestController(DD, [filename])
            sage: ntests, runner = DTT(options=DD)
            sage: runner.failures
            0
            sage: ntests >= 200 or ntests
            True
        """
        result = None
        try:
            runner = SageDocTestRunner(
                SageOutputChecker(),
                verbose=options.verbose,
                outtmpfile=outtmpfile,
                msgfile=msgfile,
                sage_options=options,
                optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS,
                baseline=baseline)
            runner.basename = self.source.basename
            runner.filename = self.source.path
            N = options.file_iterations
            results = DictAsObject({'walltime': [], 'cputime': [],
                                    'err': None, 'walltime_skips': 0})

            # multiprocessing.Process instances don't run exit
            # functions, so we run the functions added by doctests
            # when exiting this context.
            with restore_atexit(run=True):
                for it in range(N):
                    doctests, extras = self._run(runner, options, results)
                    runner.summarize(options.verbose)
                    if runner.update_results(results):
                        break

            if extras['tab']:
                results.err = 'tab'
                results.tab_linenos = extras['tab']
            if extras['line_number']:
                results.err = 'line_number'
            results.optionals = extras['optionals']
            # We subtract 1 to remove the sig_on_count() tests
            result = (sum(max(0, len(test.examples) - 1) for test in doctests),
                      results)

        except BaseException:
            exc_info = sys.exc_info()
            tb = "".join(traceback.format_exception(*exc_info))
            result = (0, DictAsObject({'err': exc_info[0], 'tb': tb}))

        if result_queue is not None:
            result_queue.put(result, False)

        return result

    def _run(self, runner, options, results):
        """
        Actually run the doctests with the right set of globals
        """
        # Import Jupyter globals to doctest the Jupyter
        # implementation of widgets and interacts
        from importlib import import_module
        sage_all = import_module(options.environment)
        dict_all = sage_all.__dict__
        # When using global environments other than sage.all,
        # make sure startup is finished so we don't get "Resolving lazy import"
        # warnings.
        from sage.misc.lazy_import import ensure_startup_finished
        ensure_startup_finished()
        # Remove '__package__' item from the globals since it is not
        # always in the globals in an actual Sage session.
        dict_all.pop('__package__', None)
        sage_namespace = RecordingDict(dict_all)
        sage_namespace['__name__'] = '__main__'
        doctests, extras = self.source.create_doctests(sage_namespace)
        timer = Timer().start()

        for test in doctests:
            result = runner.run(test)
            if options.exitfirst and result.failed:
                break

        timer.stop().annotate(runner)
        return doctests, extras
