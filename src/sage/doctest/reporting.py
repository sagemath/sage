# -*- coding: utf-8 -*-
r"""
Reporting doctest results

This module determines how doctest results are reported to the user.

It also computes the exit status in the ``error_status`` attribute of
:class:`DocTestReporter`. This is a bitwise OR of the following bits:

- 1: Doctest failure
- 2: Bad command line syntax or invalid options
- 4: Test timed out
- 8: Test exited with non-zero status
- 16: Test crashed with a signal (e.g. segmentation fault)
- 32: TAB character found
- 64: Internal error in the doctesting framework
- 128: Testing interrupted, not all tests run
- 256: Doctest contains explicit source line number

AUTHORS:

- David Roe (2012-03-27) -- initial version, based on Robert Bradshaw's code.
"""

# ****************************************************************************
#       Copyright (C) 2012 David Roe <roed.math@gmail.com>
#                          Robert Bradshaw <robertwb@gmail.com>
#                          William Stein <wstein@gmail.com>
#       Copyright (C) 2013 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sys import stdout
from signal import (SIGABRT, SIGALRM, SIGBUS, SIGFPE, SIGHUP, SIGILL,
                    SIGINT, SIGKILL, SIGPIPE, SIGQUIT, SIGSEGV, SIGTERM)
from sage.structure.sage_object import SageObject
from sage.doctest.util import count_noun
from sage.doctest.sources import DictAsObject
from .external import available_software

def signal_name(sig):
    """
    Return a string describing a signal number.

    EXAMPLES::

        sage: from signal import SIGSEGV
        sage: from sage.doctest.reporting import signal_name
        sage: signal_name(SIGSEGV)
        'segmentation fault'
        sage: signal_name(9)
        'kill signal'
        sage: signal_name(12345)
        'signal 12345'
    """
    if sig == SIGHUP:
        return "hangup"
    if sig == SIGINT:
        return "interrupt"
    if sig == SIGQUIT:
        return "quit"
    if sig == SIGILL:
        return "illegal instruction"
    if sig == SIGABRT:
        return "abort"
    if sig == SIGFPE:
        return "floating point exception"
    if sig == SIGKILL:
        return "kill signal"
    if sig == SIGSEGV:
        return "segmentation fault"
    if sig == SIGPIPE:
        return "broken pipe"
    if sig == SIGALRM:
        return "alarm"
    if sig == SIGTERM:
        return "terminate"
    if sig == SIGBUS:
        return "bus error"
    return "signal %s"%sig

class DocTestReporter(SageObject):
    """
    This class reports to the users on the results of doctests.
    """
    def __init__(self, controller):
        """
        Initialize the reporter.

        INPUT:

        - ``controller`` -- a
          :class:`sage.doctest.control.DocTestController` instance.
          Note that some methods assume that appropriate tests have
          been run by the controller.

        EXAMPLES::

            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','reporting.py')
            sage: DC = DocTestController(DocTestDefaults(),[filename])
            sage: DTR = DocTestReporter(DC)
        """
        self.controller = controller
        self.postscript = dict(lines=[], cputime=0, walltime=0)
        self.sources_completed = 0
        self.stats = {}
        self.error_status = 0

    def were_doctests_with_optional_tag_run(self, tag):
        r"""
        Return whether doctests marked with this tag were run.

        INPUT:

        - ``tag`` -- string

        EXAMPLES::

            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','reporting.py')
            sage: DC = DocTestController(DocTestDefaults(),[filename])
            sage: DTR = DocTestReporter(DC)

        ::

            sage: DTR.were_doctests_with_optional_tag_run('sage')
            True
            sage: DTR.were_doctests_with_optional_tag_run('nice_unavailable_package')
            False

        When latex is available, doctests marked with optional tag
        ``latex`` are run by default since :trac:`32174`::

            sage: filename = os.path.join(SAGE_SRC,'sage','misc','latex.py')
            sage: DC = DocTestController(DocTestDefaults(),[filename])
            sage: DTR = DocTestReporter(DC)
            sage: DTR.were_doctests_with_optional_tag_run('latex')   # optional - latex
            True

        """
        if self.controller.options.optional is True or tag in self.controller.options.optional:
            return True
        if tag in available_software.seen():
            return True
        return False

    def report_head(self, source):
        """
        Return the "sage -t [options] file.py" line as string.

        INPUT:

        - ``source`` -- a source from :mod:`sage.doctest.sources`

        EXAMPLES::

            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.env import SAGE_SRC
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','reporting.py')
            sage: DD = DocTestDefaults()
            sage: FDS = FileDocTestSource(filename,DD)
            sage: DC = DocTestController(DD, [filename])
            sage: DTR = DocTestReporter(DC)
            sage: print(DTR.report_head(FDS))
            sage -t .../sage/doctest/reporting.py

        The same with various options::

            sage: DD.long = True
            sage: print(DTR.report_head(FDS))
            sage -t --long .../sage/doctest/reporting.py
        """
        cmd = "sage -t"
        if self.controller.options.long:
            cmd += " --long"

        warnlong = self.controller.options.warn_long
        if warnlong >= 0:
            cmd += " --warn-long"
            if warnlong != 1.0:
                cmd += " %.1f" % (warnlong)
        seed = self.controller.options.random_seed
        cmd += " --random-seed={}".format(seed)
        cmd += " " + source.printpath
        return cmd

    def report(self, source, timeout, return_code, results, output, pid=None):
        """
        Report on the result of running doctests on a given source.

        This doesn't print the :meth:`report_head`, which is assumed
        to be printed already.

        INPUT:

        - ``source`` -- a source from :mod:`sage.doctest.sources`

        - ``timeout`` -- a boolean, whether doctests timed out

        - ``return_code`` -- an int, the return code of the process
          running doctests on that file.

        - ``results`` -- (irrelevant if ``timeout`` or
          ``return_code``), a tuple

          - ``ntests`` -- the number of doctests

          - ``timings`` -- a
            :class:`sage.doctest.sources.DictAsObject` instance
            storing timing data.

        - ``output`` -- a string, printed if there was some kind of
          failure

        - ``pid`` -- optional integer (default: ``None``). The pid of
          the worker process.

        EXAMPLES::

            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource, DictAsObject
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.util import Timer
            sage: from sage.env import SAGE_SRC
            sage: import os, sys, doctest
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','reporting.py')
            sage: DD = DocTestDefaults()
            sage: FDS = FileDocTestSource(filename,DD)
            sage: DC = DocTestController(DD,[filename])
            sage: DTR = DocTestReporter(DC)

        You can report a timeout::

            sage: DTR.report(FDS, True, 0, None, "Output so far...", pid=1234)
                Timed out
            **********************************************************************
            Tests run before process (pid=1234) timed out:
            Output so far...
            **********************************************************************
            sage: DTR.stats
            {'sage.doctest.reporting': {'failed': True, 'walltime': 1000000.0}}

        Or a process that returned a bad exit code::

            sage: DTR.report(FDS, False, 3, None, "Output before trouble")
                Bad exit: 3
            **********************************************************************
            Tests run before process failed:
            Output before trouble
            **********************************************************************
            sage: DTR.stats
            {'sage.doctest.reporting': {'failed': True, 'walltime': 1000000.0}}

        Or a process that segfaulted::

            sage: from signal import SIGSEGV
            sage: DTR.report(FDS, False, -SIGSEGV, None, "Output before trouble")
                Killed due to segmentation fault
            **********************************************************************
            Tests run before process failed:
            Output before trouble
            **********************************************************************
            sage: DTR.stats
            {'sage.doctest.reporting': {'failed': True, 'walltime': 1000000.0}}

        Report a timeout with results and a ``SIGKILL``::

            sage: from signal import SIGKILL
            sage: DTR.report(FDS, True, -SIGKILL, (1,None), "Output before trouble")
                Timed out after testing finished (and interrupt failed)
            **********************************************************************
            Tests run before process timed out:
            Output before trouble
            **********************************************************************
            sage: DTR.stats
            {'sage.doctest.reporting': {'failed': True, 'walltime': 1000000.0}}

        This is an internal error since results is None::

            sage: DTR.report(FDS, False, 0, None, "All output")
                Error in doctesting framework (bad result returned)
            **********************************************************************
            Tests run before error:
            All output
            **********************************************************************
            sage: DTR.stats
            {'sage.doctest.reporting': {'failed': True, 'walltime': 1000000.0}}

        Or tell the user that everything succeeded::

            sage: doctests, extras = FDS.create_doctests(globals())
            sage: runner = SageDocTestRunner(SageOutputChecker(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: Timer().start().stop().annotate(runner)
            sage: D = DictAsObject({'err':None})
            sage: runner.update_results(D)
            0
            sage: DTR.report(FDS, False, 0, (sum([len(t.examples) for t in doctests]), D), "Good tests")
                [... tests, ... s]
            sage: DTR.stats
            {'sage.doctest.reporting': {'walltime': ...}}

        Or inform the user that some doctests failed::

            sage: runner.failures = 1
            sage: runner.update_results(D)
            1
            sage: DTR.report(FDS, False, 0, (sum([len(t.examples) for t in doctests]), D), "Doctest output including the failure...")
                [... tests, 1 failure, ... s]

        If the user has requested that we report on skipped doctests,
        we do so::

            sage: DC.options = DocTestDefaults(show_skipped=True)
            sage: from collections import defaultdict
            sage: optionals = defaultdict(int)
            sage: optionals['magma'] = 5; optionals['long time'] = 4; optionals[''] = 1; optionals['not tested'] = 2
            sage: D = DictAsObject(dict(err=None,optionals=optionals))
            sage: runner.failures = 0
            sage: runner.update_results(D)
            0
            sage: DTR.report(FDS, False, 0, (sum([len(t.examples) for t in doctests]), D), "Good tests")
                1 unlabeled test not run
                4 long tests not run
                5 magma tests not run
                2 not tested tests not run
                0 tests not run because we ran out of time
                [... tests, ... s]

        Test an internal error in the reporter::

            sage: DTR.report(None, None, None, None, None)
            Traceback (most recent call last):
            ...
            AttributeError: 'NoneType' object has no attribute 'basename'

        The only-errors mode does not output anything on success::

            sage: DD = DocTestDefaults(only_errors=True)
            sage: FDS = FileDocTestSource(filename, DD)
            sage: DC = DocTestController(DD, [filename])
            sage: DTR = DocTestReporter(DC)
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: runner = SageDocTestRunner(SageOutputChecker(), verbose=False, sage_options=DD, optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: Timer().start().stop().annotate(runner)
            sage: D = DictAsObject({'err':None})
            sage: runner.update_results(D)
            0
            sage: DTR.report(FDS, False, 0, (sum([len(t.examples) for t in doctests]), D), "Good tests")

        However, failures are still output in the errors-only mode::

            sage: runner.failures = 1
            sage: runner.update_results(D)
            1
            sage: DTR.report(FDS, False, 0, (sum([len(t.examples) for t in doctests]), D), "Failed test")
                [... tests, 1 failure, ... s]
        """
        log = self.controller.log
        process_name = 'process (pid={0})'.format(pid) if pid else 'process'
        try:
            postscript = self.postscript
            stats = self.stats
            basename = source.basename
            cmd = self.report_head(source)
            try:
                ntests, result_dict = results
            except (TypeError, ValueError):
                ntests = 0
                result_dict = DictAsObject(dict(err='badresult'))
            if timeout:
                fail_msg = "Timed out"
                if ntests > 0:
                    fail_msg += " after testing finished"
                if return_code > 0:
                    fail_msg += " (with error after interrupt)"
                elif return_code < 0:
                    sig = -return_code
                    if sig == SIGQUIT:
                        pass  # and interrupt succeeded
                    elif sig == SIGKILL:
                        fail_msg += " (and interrupt failed)"
                    else:
                        fail_msg += " (with %s after interrupt)"%signal_name(sig)
                log("    %s\n%s\nTests run before %s timed out:"%(fail_msg, "*"*70, process_name))
                log(output)
                log("*"*70)
                postscript['lines'].append(cmd + "  # %s"%fail_msg)
                stats[basename] = dict(failed=True, walltime=1e6)
                self.error_status |= 4
            elif return_code:
                if return_code > 0:
                    fail_msg = "Bad exit: %s"%return_code
                else:
                    fail_msg = "Killed due to %s"%signal_name(-return_code)
                if ntests > 0:
                    fail_msg += " after testing finished"
                log("    %s\n%s\nTests run before %s failed:"%(fail_msg,"*"*70, process_name))
                log(output)
                log("*"*70)
                postscript['lines'].append(cmd + "  # %s" % fail_msg)
                stats[basename] = dict(failed=True, walltime=1e6)
                self.error_status |= (8 if return_code > 0 else 16)
            else:
                if hasattr(result_dict, 'walltime') and hasattr(result_dict.walltime, '__len__') and len(result_dict.walltime) > 0:
                    wall = sum(result_dict.walltime) / len(result_dict.walltime)
                else:
                    wall = 1e6
                if hasattr(result_dict, 'cputime') and hasattr(result_dict.cputime, '__len__') and len(result_dict.cputime) > 0:
                    cpu = sum(result_dict.cputime) / len(result_dict.cputime)
                else:
                    cpu = 1e6
                if result_dict.err == 'badresult':
                    log("    Error in doctesting framework (bad result returned)\n%s\nTests run before error:"%("*"*70))
                    log(output)
                    log("*"*70)
                    postscript['lines'].append(cmd + "  # Testing error: bad result")
                    self.error_status |= 64
                elif result_dict.err == 'noresult':
                    log("    Error in doctesting framework (no result returned)\n%s\nTests run before error:"%("*"*70))
                    log(output)
                    log("*"*70)
                    postscript['lines'].append(cmd + "  # Testing error: no result")
                    self.error_status |= 64
                elif result_dict.err == 'tab':
                    if len(result_dict.tab_linenos) > 5:
                        result_dict.tab_linenos[3:-1] = "..."
                    tabs = " " + ",".join(result_dict.tab_linenos)
                    if len(result_dict.tab_linenos) > 1:
                        tabs = "s" + tabs
                    log("    Error: TAB character found at line%s"%(tabs))
                    postscript['lines'].append(cmd + "  # Tab character found")
                    self.error_status |= 32
                elif result_dict.err == 'line_number':
                    log("    Error: Source line number found")
                    postscript['lines'].append(cmd + "  # Source line number found")
                    self.error_status |= 256
                elif result_dict.err is not None:
                    # This case should not occur
                    if result_dict.err is True:
                        fail_msg = "Error in doctesting framework"
                    else:
                        if hasattr(result_dict.err, '__name__'):
                            err = result_dict.err.__name__
                        else:
                            err = repr(result_dict.err)
                        fail_msg = "%s in doctesting framework"%err

                    log("    %s\n%s"%(fail_msg, "*"*70))
                    if output:
                        log("Tests run before doctest exception:\n" + output)
                        log("*"*70)
                    postscript['lines'].append(cmd + "  # %s"%fail_msg)
                    if hasattr(result_dict, 'tb'):
                        log(result_dict.tb)
                    if hasattr(result_dict, 'walltime'):
                        stats[basename] = dict(failed=True, walltime=wall)
                    else:
                        stats[basename] = dict(failed=True, walltime=1e6)
                    self.error_status |= 64
                if result_dict.err is None or result_dict.err == 'tab':
                    f = result_dict.failures
                    if f:
                        postscript['lines'].append(cmd + "  # %s failed" % (count_noun(f, "doctest")))
                        self.error_status |= 1
                    if f or result_dict.err == 'tab':
                        stats[basename] = dict(failed=True, walltime=wall)
                    else:
                        stats[basename] = dict(walltime=wall)
                    postscript['cputime'] += cpu
                    postscript['walltime'] += wall

                    try:
                        optionals = result_dict.optionals
                    except AttributeError:
                        optionals = dict()
                    for tag in sorted(optionals):
                        nskipped = optionals[tag]
                        if tag == "long time":
                            if not self.controller.options.long:
                                if self.controller.options.show_skipped:
                                    log("    %s not run"%(count_noun(nskipped, "long test")))
                        elif tag == "not tested":
                            if self.controller.options.show_skipped:
                                log("    %s not run"%(count_noun(nskipped, "not tested test")))
                        elif tag == "not implemented":
                            if self.controller.options.show_skipped:
                                log("    %s for not implemented functionality not run"%(count_noun(nskipped, "test")))
                        else:
                            if not self.were_doctests_with_optional_tag_run(tag):
                                if tag == "bug":
                                    if self.controller.options.show_skipped:
                                        log("    %s not run due to known bugs"%(count_noun(nskipped, "test")))
                                elif tag == "":
                                    if self.controller.options.show_skipped:
                                        log("    %s not run"%(count_noun(nskipped, "unlabeled test")))
                                else:
                                    if self.controller.options.show_skipped:
                                        log("    %s not run"%(count_noun(nskipped, tag + " test")))

                    nskipped = result_dict.walltime_skips
                    if self.controller.options.show_skipped:
                        log("    %s not run because we ran out of time"%(count_noun(nskipped, "test")))

                    if nskipped != 0:
                        # It would be nice to report "a/b tests run" instead of
                        # the percentage that is printed here.  However, it is
                        # not clear how to pull out the actual part of "ntests"
                        # that has been run for a variety of reasons, such as
                        # the sig_on_count() tests, the possibility to run
                        # tests multiple times, and some other unclear mangling
                        # of these numbers that was not clear to the author.
                        ntests_run = result_dict.tests
                        total = "%d%% of tests run"%(round(100*ntests_run/float(ntests_run + nskipped)))
                    else:
                        total = count_noun(ntests, "test")
                    if not (self.controller.options.only_errors and not f):
                        log("    [%s, %s%.2f s]" % (total, "%s, "%(count_noun(f, "failure")) if f else "", wall))

            self.sources_completed += 1

        except Exception:
            import traceback
            log(traceback.format_exc(), end="")


    def finalize(self):
        """
        Print out the postscript that summarizes the doctests that were run.

        EXAMPLES:

        First we have to set up a bunch of stuff::

            sage: from sage.doctest.reporting import DocTestReporter
            sage: from sage.doctest.control import DocTestController, DocTestDefaults
            sage: from sage.doctest.sources import FileDocTestSource, DictAsObject
            sage: from sage.doctest.forker import SageDocTestRunner
            sage: from sage.doctest.parsing import SageOutputChecker
            sage: from sage.doctest.util import Timer
            sage: from sage.env import SAGE_SRC
            sage: import os, sys, doctest
            sage: filename = os.path.join(SAGE_SRC,'sage','doctest','reporting.py')
            sage: DD = DocTestDefaults()
            sage: FDS = FileDocTestSource(filename,DD)
            sage: DC = DocTestController(DD,[filename])
            sage: DTR = DocTestReporter(DC)

        Now we pretend to run some doctests::

            sage: DTR.report(FDS, True, 0, None, "Output so far...", pid=1234)
                Timed out
            **********************************************************************
            Tests run before process (pid=1234) timed out:
            Output so far...
            **********************************************************************
            sage: DTR.report(FDS, False, 3, None, "Output before bad exit")
                Bad exit: 3
            **********************************************************************
            Tests run before process failed:
            Output before bad exit
            **********************************************************************
            sage: doctests, extras = FDS.create_doctests(globals())
            sage: runner = SageDocTestRunner(SageOutputChecker(), verbose=False, sage_options=DD,optionflags=doctest.NORMALIZE_WHITESPACE|doctest.ELLIPSIS)
            sage: t = Timer().start().stop()
            sage: t.annotate(runner)
            sage: DC.timer = t
            sage: D = DictAsObject({'err':None})
            sage: runner.update_results(D)
            0
            sage: DTR.report(FDS, False, 0, (sum([len(t.examples) for t in doctests]), D), "Good tests")
                [... tests, ... s]
            sage: runner.failures = 1
            sage: runner.update_results(D)
            1
            sage: DTR.report(FDS, False, 0, (sum([len(t.examples) for t in doctests]), D), "Doctest output including the failure...")
                [... tests, 1 failure, ... s]

        Now we can show the output of finalize::

            sage: DC.sources = [None] * 4 # to fool the finalize method
            sage: DTR.finalize()
            ----------------------------------------------------------------------
            sage -t .../sage/doctest/reporting.py  # Timed out
            sage -t .../sage/doctest/reporting.py  # Bad exit: 3
            sage -t .../sage/doctest/reporting.py  # 1 doctest failed
            ----------------------------------------------------------------------
            Total time for all tests: 0.0 seconds
                cpu time: 0.0 seconds
                cumulative wall time: 0.0 seconds

        If we interrupted doctests, then the number of files tested
        will not match the number of sources on the controller::

            sage: DC.sources = [None] * 6
            sage: DTR.finalize()
            <BLANKLINE>
            ----------------------------------------------------------------------
            sage -t .../sage/doctest/reporting.py  # Timed out
            sage -t .../sage/doctest/reporting.py  # Bad exit: 3
            sage -t .../sage/doctest/reporting.py  # 1 doctest failed
            Doctests interrupted: 4/6 files tested
            ----------------------------------------------------------------------
            Total time for all tests: 0.0 seconds
                cpu time: 0.0 seconds
                cumulative wall time: 0.0 seconds
        """
        log = self.controller.log
        postscript = self.postscript
        if self.sources_completed < len(self.controller.sources) * self.controller.options.global_iterations:
            postscript['lines'].append("Doctests interrupted: %s/%s files tested"%(self.sources_completed, len(self.controller.sources)))
            self.error_status |= 128
        elif not postscript['lines']:
            postscript['lines'].append("All tests passed!")
        log('-' * 70)
        log("\n".join(postscript['lines']))
        log('-' * 70)
        log("Total time for all tests: %.1f seconds" % self.controller.timer.walltime)
        log("    cpu time: %.1f seconds" % postscript['cputime'])
        log("    cumulative wall time: %.1f seconds" % postscript['walltime'])
        stdout.flush()
