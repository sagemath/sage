# sage_setup: distribution = sagemath-repl
"""
Classes involved in doctesting

This module controls the various classes involved in doctesting.

AUTHORS:

- David Roe (2012-03-27) -- initial version, based on Robert Bradshaw's code.
"""
# ****************************************************************************
#       Copyright (C) 2012-2013 David Roe <roed.math@gmail.com>
#                     2012-2013 Robert Bradshaw <robertwb@gmail.com>
#                     2012      William Stein <wstein@gmail.com>
#                     2013      R. Andrew Ohana
#                     2013-2014 Volker Braun
#                     2013-2018 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#                     2013-2021 John H. Palmieri
#                     2017      Erik M. Bray
#                     2017-2021 Frédéric Chapoton
#                     2018      Sébastien Labbé
#                     2019      François Bissey
#                     2020-2023 Matthias Koeppe
#                     2022      Michael Orlitzky
#                     2022      Sebastian Oehms
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import importlib
import random
import os
import sys
import time
import json
import shlex
import types
import sage.misc.flatten
import sage.misc.randstate as randstate
from sage.structure.sage_object import SageObject
from sage.env import DOT_SAGE, SAGE_LIB, SAGE_SRC, SAGE_VENV, SAGE_EXTCODE
from sage.misc.temporary_file import tmp_dir
from cysignals.signals import AlarmInterrupt, init_cysignals

from .sources import FileDocTestSource, DictAsObject, get_basename
from .forker import DocTestDispatcher
from .reporting import DocTestReporter
from .util import Timer, count_noun, dict_difference
from .external import available_software
from .parsing import parse_optional_tags, parse_file_optional_tags, unparse_optional_tags, \
     nodoctest_regex, optionaltag_regex, optionalfiledirective_regex


# Optional tags which are always automatically added

auto_optional_tags = set()

class DocTestDefaults(SageObject):
    """
    This class is used for doctesting the Sage doctest module.

    INPUT:

    - ``runtest_default`` -- (boolean, default ``False``); if ``True``,
      fills in attribute to be the same as the defaults defined in
      ``sage-runtests``. If ``False``, change defaults in a few places
      for use in doctests of the doctester, which is mostly to make
      doctesting more predictable.

    - ``**kwds`` -- attributes to override defaults

    EXAMPLES::

        sage: from sage.doctest.control import DocTestDefaults
        sage: D = DocTestDefaults(); D
        DocTestDefaults()
        sage: D.timeout
        -1

    Keyword arguments become attributes::

        sage: D = DocTestDefaults(timeout=100); D
        DocTestDefaults(timeout=100)
        sage: D.timeout
        100

    The defaults for ``sage-runtests``::

        sage: D = DocTestDefaults(runtest_default=True); D
        DocTestDefaults(abspath=False, file_iterations=0, global_iterations=0,
                        optional='sage,optional', random_seed=None,
                        stats_path='.../timings2.json')
    """
    def __init__(self, runtest_default=False, **kwds):
        """
        Edit these parameters after creating an instance.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: D = DocTestDefaults()
            sage: 'sage' in D.optional
            True
        """
        # NOTE that these are NOT the defaults used by the sage-runtests
        # script (which is what gets invoked when running `sage -t`).
        # These are only basic defaults when invoking the doctest runner
        # from Python, which is not the typical use case.
        self.nthreads = 1
        self.serial = False
        self.timeout = -1
        self.die_timeout = -1
        self.all = False
        self.installed = False
        self.logfile = None
        self.long = False
        self.warn_long = -1.0
        self.randorder = None
        self.random_seed = None if runtest_default else 0
        self.global_iterations = 0 if runtest_default else 1
        self.file_iterations = 0 if runtest_default else 1
        self.environment = "sage.repl.ipython_kernel.all_jupyter"
        self.initial = False
        self.exitfirst = False
        self.force_lib = False
        self.if_installed = False
        self.abspath = not runtest_default
        self.verbose = False
        self.debug = False
        self.only_errors = False
        self.gdb = False
        self.lldb = False
        self.valgrind = False
        self.massif = False
        self.cachegrind = False
        self.omega = False
        self.failed = False
        self.new = False
        self.show_skipped = False
        self.target_walltime = -1
        self.baseline_stats_path = None
        self.format = "sage"

        # sage-runtests contains more optional tags. Technically, adding
        # auto_optional_tags here is redundant, since that is added
        # automatically anyway. However, this default is still used for
        # displaying user-defined optional tags and we don't want to see
        # the auto_optional_tags there.
        if runtest_default:
            self.optional = ','.join(['sage', 'optional'])
        else:
            self.optional = {'sage'} | auto_optional_tags
        self.hide = ''
        self.probe = ''

        # > 0: always run GC before every test
        # < 0: disable GC
        self.gc = 0

        # We don't want to use the real stats file by default so that
        # we don't overwrite timings for the actual running doctests.
        self.stats_path = os.path.join(
            DOT_SAGE, "timings2.json" if runtest_default else "timings_dt_test.json")
        self.__dict__.update(kwds)

    def _repr_(self):
        """
        Return the print representation.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: DocTestDefaults(timeout=100, foobar='hello')
            DocTestDefaults(foobar='hello', timeout=100)
        """
        s = "DocTestDefaults("
        for k in sorted(dict_difference(self.__dict__, DocTestDefaults().__dict__).keys()):
            if s[-1] != "(":
                s += ", "
            s += str(k) + "=" + repr(getattr(self, k))
        s += ")"
        return s

    def __eq__(self, other):
        """
        Comparison by __dict__.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: DD1 = DocTestDefaults(long=True)
            sage: DD2 = DocTestDefaults(long=True)
            sage: DD1 == DD2
            True
        """
        if not isinstance(other, DocTestDefaults):
            return False
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        """
        Test for non-equality.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults
            sage: DD1 = DocTestDefaults(long=True)
            sage: DD2 = DocTestDefaults(long=True)
            sage: DD1 != DD2
            False
        """
        return not (self == other)


def skipdir(dirname):
    """
    Return ``True`` if and only if the directory ``dirname`` should not be
    doctested.

    EXAMPLES::

        sage: from sage.doctest.control import skipdir
        sage: skipdir(sage.env.SAGE_SRC)
        False
        sage: skipdir(os.path.join(sage.env.SAGE_SRC, "sage", "doctest", "tests"))
        True
    """
    if os.path.exists(os.path.join(dirname, "nodoctest.py")) or os.path.exists(os.path.join(dirname, "nodoctest")):
        return True
    return False


def skipfile(filename, tested_optional_tags=False, *,
             if_installed=False, log=None):
    """
    Return ``True`` if and only if the file ``filename`` should not be doctested.

    INPUT:

    - ``filename`` -- name of a file

    - ``tested_optional_tags`` -- list or tuple or set of optional tags to test,
      or ``False`` (no optional test) or ``True`` (all optional tests)

    - ``if_installed`` -- boolean (default: ``False``); whether to skip Python/Cython files
      that are not installed as modules

    - ``log`` -- function to call with log messages, or ``None``

    If ``filename`` contains a line of the form ``"# sage.doctest:
    optional - xyz")``, then this will return ``False`` if "xyz" is in
    ``tested_optional_tags``. Otherwise, it returns the matching tag
    ("optional - xyz").

    EXAMPLES::

        sage: from sage.doctest.control import skipfile
        sage: skipfile("skipme.c")
        True
        sage: filename = tmp_filename(ext='.pyx')
        sage: skipfile(filename)
        False
        sage: with open(filename, "w") as f:
        ....:     _ = f.write("# nodoctest")
        sage: skipfile(filename)
        True
        sage: with open(filename, "w") as f:
        ....:     _ = f.write("# sage.doctest: "    # broken in two source lines to avoid the pattern
        ....:                 "optional - xyz")     # of relint (multiline_doctest_comment)
        sage: skipfile(filename, False)
        'optional - xyz'
        sage: bool(skipfile(filename, False))
        True
        sage: skipfile(filename, ['abc'])
        'optional - xyz'
        sage: skipfile(filename, ['abc', 'xyz'])
        False
        sage: skipfile(filename, True)
        False
    """
    if filename.endswith('.rst.txt'):
        ext = '.rst.txt'
    if filename.endswith('__main__.py'):
        if log:
            log(f"Skipping '{filename}' because it is a __main__.py file")
        return True
    _ , ext = os.path.splitext(filename)
    # .rst.txt appear in the installed documentation in subdirectories named "_sources"
    if ext not in ('.py', '.pyx', '.pxd', '.pxi', '.sage', '.spyx', '.rst', '.tex', '.rst.txt'):
        if log:
            log(f"Skipping '{filename}' because it does not have one of the recognized file name extensions")
        return True
    if if_installed and ext not in ('.py', '.pyx'):
        if log:
            log(f"Skipping '{filename}' because it is not the source file of a Python module")
        return True
    if "jupyter_execute" in filename:
        if log:
            log(f"Skipping '{filename}' because it is created by the jupyter-sphinx extension for internal use and should not be tested")
        return True
    if if_installed:
        module_name = get_basename(filename)
        try:
            if not importlib.util.find_spec(module_name):  # tries to import the containing package
                if log:
                    log(f"Skipping '{filename}' because module {module_name} is not present in the venv")
                return True
        except ModuleNotFoundError as e:
            if log:
                log(f"Skipping '{filename}' because module {e.name} cannot be imported")
            return True

    with open(filename) as F:
        file_optional_tags = parse_file_optional_tags(enumerate(F))

    if 'not tested' in file_optional_tags:
        if log:
            log(f"Skipping '{filename}' because it is marked 'nodoctest'")
        return True

    if tested_optional_tags is False:
        if file_optional_tags:
            file_tag_string = unparse_optional_tags(file_optional_tags, prefix='')
            if log:
                log(f"Skipping '{filename}' because it is marked '# {file_tag_string}'")
            return file_tag_string

    elif tested_optional_tags is not True:
        extra = {tag for tag in file_optional_tags
                 if tag not in tested_optional_tags}
        if extra:
            file_tag_string = unparse_optional_tags(file_optional_tags, prefix='')
            if log:
                log(f"Skipping '{filename}' because it is marked '{file_tag_string}'")
            return file_tag_string

    return False


class Logger:
    r"""
    File-like object which implements writing to multiple files at
    once.

    EXAMPLES::

        sage: from sage.doctest.control import Logger
        sage: with open(tmp_filename(), "w+") as t:
        ....:     L = Logger(sys.stdout, t)
        ....:     _ = L.write("hello world\n")
        ....:     _ = t.seek(0)
        ....:     t.read()
        hello world
        'hello world\n'
    """
    def __init__(self, *files):
        r"""
        Initialize the logger for writing to all files in ``files``.

        TESTS::

            sage: from sage.doctest.control import Logger
            sage: Logger().write("hello world\n")  # no-op
        """
        self.files = list(files)

    def write(self, x):
        r"""
        Write ``x`` to all files.

        TESTS::

            sage: from sage.doctest.control import Logger
            sage: Logger(sys.stdout).write("hello world\n")
            hello world
        """
        for f in self.files:
            f.write(x)

    def flush(self):
        """
        Flush all files.

        TESTS::

            sage: from sage.doctest.control import Logger
            sage: Logger(sys.stdout).flush()
        """
        for f in self.files:
            f.flush()


class DocTestController(SageObject):
    """
    This class controls doctesting of files.

    After creating it with appropriate options, call the :meth:`run` method to run the doctests.
    """
    def __init__(self, options, args):
        """
        Initialization.

        INPUT:

        - ``options`` -- either options generated from the command line by sage-runtests
          or a DocTestDefaults object (possibly with some entries modified)
        - ``args`` -- list of filenames to doctest

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: DC = DocTestController(DocTestDefaults(), [])
            sage: DC
            DocTest Controller
        """
        # First we modify options to take environment variables into
        # account and check compatibility of the user's specified
        # options.
        if options.timeout < 0:
            if options.gdb or options.lldb or options.debug:
                # Interactive debuggers: "infinite" timeout
                options.timeout = 0
            elif options.valgrind or options.massif or options.cachegrind or options.omega:
                # Non-interactive debuggers: 48 hours
                options.timeout = int(os.getenv('SAGE_TIMEOUT_VALGRIND', 48 * 60 * 60))
            elif options.long:
                options.timeout = int(os.getenv('SAGE_TIMEOUT_LONG', 30 * 60))
            else:
                options.timeout = int(os.getenv('SAGE_TIMEOUT', 10 * 60))
            # For non-default GC options, double the timeout
            if options.gc:
                options.timeout *= 2
        if options.nthreads == 0:
            options.nthreads = int(os.getenv('SAGE_NUM_THREADS_PARALLEL', 1))
        if options.failed and not (args or options.new):
            # If the user doesn't specify any files then we rerun all failed files.
            options.all = True
        if options.global_iterations == 0:
            options.global_iterations = int(os.environ.get('SAGE_TEST_GLOBAL_ITER', 1))
        if options.file_iterations == 0:
            options.file_iterations = int(os.environ.get('SAGE_TEST_ITER', 1))
        if options.debug:
            if options.nthreads > 1:
                print("Debugging requires single-threaded operation, setting number of threads to 1.")
            if options.logfile:
                print("Debugging is not compatible with logging, disabling logfile.")
            options.serial = True
            options.logfile = None
        if options.serial:
            options.nthreads = 1
        if options.verbose:
            options.show_skipped = True

        options.hidden_features = set()
        if isinstance(options.hide, str):
            if not len(options.hide):
                options.hide = set()
            else:
                s = options.hide.lower()
                options.hide = set(s.split(','))
                for h in options.hide:
                    if not optionaltag_regex.search(h):
                        raise ValueError('invalid optional tag {!r}'.format(h))
            if 'all' in options.hide:
                options.hide.discard('all')
                from sage.features.all import all_features
                feature_names = {f.name for f in all_features() if not f.is_standard()}
                options.hide = options.hide.union(feature_names)
            if 'optional' in options.hide:
                options.hide.discard('optional')
                from sage.features.all import all_features
                feature_names = {f.name for f in all_features() if f.is_optional()}
                options.hide = options.hide.union(feature_names)

        options.disabled_optional = set()
        if isinstance(options.optional, str):
            s = options.optional.lower()
            options.optional = set(s.split(','))
            if "all" in options.optional:
                # Special case to run all optional tests
                options.optional = True
            else:
                # We replace the 'optional' tag by all optional
                # packages for which the installed version matches the
                # latest available version (this implies in particular
                # that the package is actually installed).
                if 'optional' in options.optional:
                    options.optional.discard('optional')
                    from sage.misc.package import list_packages
                    for pkg in list_packages('optional', local=True).values():
                        if pkg.name in options.hide:
                            continue
                        # Skip features for which we have a more specific runtime feature test.
                        if pkg.name in ['bliss', 'coxeter3', 'mcqd', 'meataxe', 'sirocco', 'tdlib']:
                            continue
                        if pkg.is_installed() and pkg.installed_version == pkg.remote_version:
                            options.optional.add(pkg.name)

                    from sage.features import package_systems
                    options.optional.update(system.name
                                            for system in package_systems())
                # Check that all tags are valid
                for o in options.optional:
                    if o.startswith('!'):
                        if not optionaltag_regex.search(o[1:]):
                            raise ValueError('invalid optional tag {!r}'.format(o))
                        options.disabled_optional.add(o[1:])
                    elif not optionaltag_regex.search(o):
                        raise ValueError('invalid optional tag {!r}'.format(o))

                options.optional |= auto_optional_tags
                options.optional -= options.disabled_optional

        if isinstance(options.probe, str):
            if options.probe == 'none':
                options.probe = ''
            s = options.probe.lower()
            if not s:
                options.probe = set()
            else:
                options.probe = set(s.split(','))
                if "all" in options.probe:
                    # Special case to probe all features that are not present
                    options.probe = True
                else:
                    # Check that all tags are valid
                    for o in options.probe:
                        if not optionaltag_regex.search(o):
                            raise ValueError('invalid optional tag {!r}'.format(o))

        self.options = options

        self.files = args
        if options.logfile:
            if not isinstance(options.logfile, str):
                # file from sage-runtests
                self.logfile = options.logfile
            else:
                # string from DocTestDefaults
                try:
                    self.logfile = open(options.logfile, 'a')
                except OSError:
                    print("Unable to open logfile {!r}\nProceeding without logging.".format(options.logfile))
                    self.logfile = None
        else:
            self.logfile = None

        # Flush any diagnostic messages we just printed
        sys.stdout.flush()
        sys.stderr.flush()

        # In serial mode, we run just one process. Then the doctests
        # will interfere with the output logging (both use stdout).
        # To solve this, we create real_stdout which will always
        # write to the actual standard output, regardless of
        # redirections.
        if options.serial:
            self._real_stdout = os.fdopen(os.dup(sys.stdout.fileno()), "w")
            self._close_stdout = True
        else:
            # Parallel mode: no special tricks needed
            self._real_stdout = sys.stdout
            self._close_stdout = False

        if self.logfile is None:
            self.logger = self._real_stdout
        else:
            self.logger = Logger(self._real_stdout, self.logfile)

        self.stats = {}
        self.load_stats(options.stats_path)
        self.baseline_stats = {}
        if options.baseline_stats_path:
            self.load_baseline_stats(options.baseline_stats_path)
        self._init_warn_long()

        if self.options.random_seed is None:
            randstate.set_random_seed()
            self.options.random_seed = randstate.initial_seed()

    def __del__(self):
        if getattr(self, 'logfile', None) is not None:
            self.logfile.close()

        if getattr(self, '_close_stdout', False):
            self._real_stdout.close()

    def _init_warn_long(self):
        """
        Pick a suitable default for the ``--warn-long`` option if not
        specified.

        It is desirable to have all tests (even ``# long`` ones)
        finish in less than about 5 seconds. Longer tests typically
        don't add coverage, they just make testing slow.

        The default used here is 5 seconds, unless `--long` was used,
        in which case it is 30 seconds.

        TESTS:

        Ensure that the user's command-line options are not changed::

            sage: from sage.doctest.control import (DocTestDefaults,
            ....:                                   DocTestController)
            sage: DC = DocTestController(DocTestDefaults(), [])
            sage: DC.options.warn_long = 5.0
            sage: DC._init_warn_long()
            sage: DC.options.warn_long
            5.00000000000000
        """
        # default is -1.0
        if self.options.warn_long >= 0:     # Specified on the command line
            return

        # The developer's guide says that even a "long time" test
        # should ideally complete in under five seconds, so we're
        # being rather generous here.
        self.options.warn_long = 5.0
        if self.options.long:
            self.options.warn_long = 30.0

    def second_on_modern_computer(self):
        """
        Return the wall time equivalent of a second on a modern computer.

        OUTPUT:

        Float. The wall time on your computer that would be equivalent
        to one second on a modern computer. Unless you have kick-ass
        hardware this should always be >= 1.0. This raises a
        :exc:`RuntimeError` if there are no stored timings to use as
        benchmark.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: DC = DocTestController(DocTestDefaults(), [])
            sage: DC.second_on_modern_computer()   # not tested
        """
        from sage.misc.superseded import deprecation
        deprecation(32981, "this method is no longer used by the sage library and will eventually be removed")

        if len(self.stats) == 0:
            raise RuntimeError('no stored timings available')
        success = []
        failed = []
        for mod in self.stats.values():
            if mod.get('failed', False):
                failed.append(mod['walltime'])
            else:
                success.append(mod['walltime'])
        if len(success) < 2500:
            raise RuntimeError('too few successful tests, not using stored timings')
        if len(failed) > 20:
            raise RuntimeError('too many failed tests, not using stored timings')
        expected = 12800.0       # Core i7 Quad-Core 2014
        return sum(success) / expected

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: DC = DocTestController(DocTestDefaults(), [])
            sage: repr(DC) # indirect doctest
            'DocTest Controller'
        """
        return "DocTest Controller"

    def load_environment(self):
        """
        Return the module that provides the global environment.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: DC = DocTestController(DocTestDefaults(), [])
            sage: 'BipartiteGraph' in DC.load_environment().__dict__
            True
            sage: DC = DocTestController(DocTestDefaults(environment='sage.doctest.all'), [])
            sage: 'BipartiteGraph' in  DC.load_environment().__dict__
            False
            sage: 'run_doctests' in DC.load_environment().__dict__
            True
        """
        from importlib import import_module
        return import_module(self.options.environment)

    def load_baseline_stats(self, filename):
        """
        Load baseline stats.

        This must be a JSON file in the same format that :meth:`load_stats`
        expects.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: DC = DocTestController(DocTestDefaults(), [])
            sage: import json
            sage: filename = tmp_filename()
            sage: with open(filename, 'w') as stats_file:
            ....:     json.dump({'sage.doctest.control':{'failed':True}}, stats_file)
            sage: DC.load_baseline_stats(filename)
            sage: DC.baseline_stats['sage.doctest.control']
            {'failed': True}

        If the file doesn't exist, nothing happens. If there is an
        error, print a message. In any case, leave the stats alone::

            sage: d = tmp_dir()
            sage: DC.load_baseline_stats(os.path.join(d))  # Cannot read a directory
            Error loading baseline stats from ...
            sage: DC.load_baseline_stats(os.path.join(d, "no_such_file"))
            sage: DC.baseline_stats['sage.doctest.control']
            {'failed': True}
        """
        # Simply ignore non-existing files
        if not os.path.exists(filename):
            return

        try:
            with open(filename) as stats_file:
                self.baseline_stats.update(json.load(stats_file))
        except Exception as e:
            self.log("Error loading baseline stats from %s: %s" % (filename, e))

    def load_stats(self, filename):
        """
        Load stats from the most recent run(s).

        Stats are stored as a JSON file, and include information on
        which files failed tests and the walltime used for execution
        of the doctests.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: DC = DocTestController(DocTestDefaults(), [])
            sage: import json
            sage: filename = tmp_filename()
            sage: with open(filename, 'w') as stats_file:
            ....:     json.dump({'sage.doctest.control': {'walltime': 1.0r}}, stats_file)
            sage: DC.load_stats(filename)
            sage: DC.stats['sage.doctest.control']
            {'walltime': 1.0}

        If the file doesn't exist, nothing happens. If there is an
        error, print a message. In any case, leave the stats alone::

            sage: d = tmp_dir()
            sage: DC.load_stats(os.path.join(d))  # Cannot read a directory
            Error loading stats from ...
            sage: DC.load_stats(os.path.join(d, "no_such_file"))
            sage: DC.stats['sage.doctest.control']
            {'walltime': 1.0}
        """
        # Simply ignore non-existing files
        if not os.path.exists(filename):
            return

        try:
            with open(filename) as stats_file:
                self.stats.update(json.load(stats_file))
        except Exception:
            self.log("Error loading stats from %s" % filename)

    def save_stats(self, filename):
        """
        Save stats from the most recent run as a JSON file.

        WARNING: This function overwrites the file.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: DC = DocTestController(DocTestDefaults(), [])
            sage: DC.stats['sage.doctest.control'] = {'walltime': 1.0r}
            sage: filename = tmp_filename()
            sage: DC.save_stats(filename)
            sage: import json
            sage: with open(filename) as f:
            ....:     D = json.load(f)
            sage: D['sage.doctest.control']
            {'walltime': 1.0}
        """
        from sage.misc.temporary_file import atomic_write
        with atomic_write(filename) as stats_file:
            json.dump(self.stats, stats_file, sort_keys=True, indent=4)

    def log(self, s, end='\n'):
        """
        Log the string ``s + end`` (where ``end`` is a newline by default)
        to the logfile and print it to the standard output.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: DD = DocTestDefaults(logfile=tmp_filename())
            sage: DC = DocTestController(DD, [])
            sage: DC.log("hello world")
            hello world
            sage: DC.logfile.close()
            sage: with open(DD.logfile) as f:
            ....:     print(f.read())
            hello world

        In serial mode, check that logging works even if ``stdout`` is
        redirected::

            sage: DD = DocTestDefaults(logfile=tmp_filename(), serial=True)
            sage: DC = DocTestController(DD, [])
            sage: from sage.doctest.forker import SageSpoofInOut
            sage: with open(os.devnull, 'w') as devnull:
            ....:     S = SageSpoofInOut(devnull)
            ....:     S.start_spoofing()
            ....:     DC.log("hello world")
            ....:     S.stop_spoofing()
            hello world
            sage: DC.logfile.close()
            sage: with open(DD.logfile) as f:
            ....:     print(f.read())
            hello world

        Check that no duplicate logs appear, even when forking (:issue:`15244`)::

            sage: DD = DocTestDefaults(logfile=tmp_filename())
            sage: DC = DocTestController(DD, [])
            sage: DC.log("hello world")
            hello world
            sage: if os.fork() == 0:
            ....:     DC.logfile.close()
            ....:     os._exit(0)
            sage: DC.logfile.close()
            sage: with open(DD.logfile) as f:
            ....:     print(f.read())
            hello world
        """
        self.logger.write(s + end)
        self.logger.flush()

    def create_run_id(self):
        """
        Create the run id.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: DC = DocTestController(DocTestDefaults(), [])
            sage: DC.create_run_id()
            Running doctests with ID ...
        """
        self.run_id = time.strftime('%Y-%m-%d-%H-%M-%S-') + "%08x" % random.getrandbits(32)
        self.log("Running doctests with ID %s." % self.run_id)

    def add_files(self):
        r"""
        Check for the flags '--all' and '--new'.

        For each one present, this function adds the appropriate directories and files to the todo list.

        EXAMPLES::

            sage: from sage.doctest.control import (DocTestDefaults,
            ....:                                   DocTestController)
            sage: from sage.env import SAGE_SRC
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile() as f:
            ....:     DD = DocTestDefaults(all=True, logfile=f.name)
            ....:     DC = DocTestController(DD, [])
            ....:     DC.add_files()
            Doctesting ...
            sage: os.path.join(SAGE_SRC, 'sage') in DC.files
            True

        ::

            sage: DD = DocTestDefaults(new = True)
            sage: DC = DocTestController(DD, [])
            sage: DC.add_files()
            Doctesting ...
        """
        opj = os.path.join
        from sage.env import SAGE_SRC, SAGE_DOC_SRC, SAGE_ROOT, SAGE_ROOT_GIT, SAGE_DOC
        # SAGE_ROOT_GIT can be None on distributions which typically
        # only have the SAGE_LOCAL install tree but not SAGE_ROOT
        if SAGE_ROOT_GIT is not None:
            have_git = os.path.exists(SAGE_ROOT_GIT)
        else:
            have_git = False

        def all_installed_modules():
            self.log("Doctesting all installed modules of the Sage library.")
            import sage
            self.files.extend(sage.__path__)
            try:
                import sage_setup
                self.files.extend(sage_setup.__path__)
            except ImportError:
                pass
            try:
                import sage_docbuild
                self.files.extend(sage_docbuild.__path__)
            except ImportError:
                pass

        def all_installed_doc():
            if SAGE_DOC and os.path.isdir(SAGE_DOC):
                self.log("Doctesting all installed documentation sources.")
                self.files.append(SAGE_DOC)

        def all_files():
            if not SAGE_SRC:
                return all_installed_modules()
            self.log("Doctesting entire Sage library.")
            self.files.append(opj(SAGE_SRC, 'sage'))
            # Only test sage_setup and sage_docbuild if the relevant
            # imports work. They may not work if not in a build
            # environment or if the documentation build has been
            # disabled.
            try:
                import sage_setup
                self.files.append(opj(SAGE_SRC, 'sage_setup'))
            except ImportError:
                pass
            try:
                import sage_docbuild
                self.files.append(opj(SAGE_SRC, 'sage_docbuild'))
            except ImportError:
                pass

        def all_doc_sources():
            if SAGE_DOC_SRC and os.path.isdir(SAGE_DOC_SRC):
                self.log("Doctesting all documentation sources.")
                self.files.append(SAGE_DOC_SRC)
            else:
                all_installed_doc()

        if self.options.installed:
            all_installed_modules()
            all_installed_doc()

        elif self.options.all or (self.options.new and not have_git):
            all_files()
            all_doc_sources()

        elif self.options.new and have_git:
            # Get all files changed in the working repo.
            self.log("Doctesting files changed since last git commit")
            import subprocess
            change = subprocess.check_output(["git",
                                              "--git-dir=" + SAGE_ROOT_GIT,
                                              "--work-tree=" + SAGE_ROOT,
                                              "status",
                                              "--porcelain"])
            change = change.decode('utf-8')
            for line in change.split("\n"):
                if not line:
                    continue
                data = line.strip().split(' ')
                status, filename = data[0], data[-1]
                if (set(status).issubset("MARCU")
                        and filename.startswith("src/sage")
                        and (filename.endswith(".py") or
                             filename.endswith(".pyx") or
                             filename.endswith(".rst"))
                        and not skipfile(opj(SAGE_ROOT, filename),
                                         bool(self.options.optional),
                                         if_installed=self.options.if_installed)):
                    self.files.append(os.path.relpath(opj(SAGE_ROOT, filename)))

    def expand_files_into_sources(self):
        r"""
        Expand ``self.files``, which may include directories, into a
        list of :class:`sage.doctest.FileDocTestSource`

        This function also handles the optional command line option.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: dirname = os.path.join(SAGE_SRC, 'sage', 'doctest')
            sage: DD = DocTestDefaults(optional='all')
            sage: DC = DocTestController(DD, [dirname])
            sage: DC.expand_files_into_sources()
            sage: len(DC.sources)
            15
            sage: DC.sources[0].options.optional
            True

        ::

            sage: DD = DocTestDefaults(optional='magma,guava')
            sage: DC = DocTestController(DD, [dirname])
            sage: DC.expand_files_into_sources()
            sage: all(t in DC.sources[0].options.optional for t in ['magma','guava'])
            True

        We check that files are skipped appropriately::

            sage: dirname = tmp_dir()
            sage: filename = os.path.join(dirname, 'not_tested.py')
            sage: with open(filename, 'w') as f:
            ....:     _ = f.write("#"*80 + "\n\n\n\n## nodoctest\n    sage: 1+1\n    4")
            sage: DC = DocTestController(DD, [dirname])
            sage: DC.expand_files_into_sources()
            sage: DC.sources
            []

        The directory ``sage/doctest/tests`` contains ``nodoctest.py``
        but the files should still be tested when that directory is
        explicitly given (as opposed to being recursed into)::

            sage: DC = DocTestController(DD, [os.path.join(SAGE_SRC, 'sage', 'doctest', 'tests')])
            sage: DC.expand_files_into_sources()
            sage: len(DC.sources) >= 10
            True
        """
        def expand():
            for path in self.files:
                if os.path.isdir(path):
                    for root, dirs, files in os.walk(path):
                        for dir in list(dirs):
                            if dir[0] == "." or skipdir(os.path.join(root, dir)):
                                dirs.remove(dir)
                        for file in files:
                            if not skipfile(os.path.join(root, file),
                                            bool(self.options.optional),
                                            if_installed=self.options.if_installed):
                                yield os.path.join(root, file)
                else:
                    if not skipfile(path, bool(self.options.optional),
                                    if_installed=self.options.if_installed, log=self.log):  # log when directly specified filenames are skipped
                        yield path
        self.sources = [FileDocTestSource(path, self.options) for path in expand()]

    def filter_sources(self):
        """

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: dirname = os.path.join(SAGE_SRC, 'sage', 'doctest')
            sage: DD = DocTestDefaults(failed=True)
            sage: DC = DocTestController(DD, [dirname])
            sage: DC.expand_files_into_sources()
            sage: for i, source in enumerate(DC.sources):
            ....:     DC.stats[source.basename] = {'walltime': 0.1r * (i+1)}
            sage: DC.stats['sage.doctest.control'] = {'failed': True, 'walltime': 1.0r}
            sage: DC.filter_sources()
            Only doctesting files that failed last test.
            sage: len(DC.sources)
            1
        """
        # Filter the sources to only include those with failing doctests if the --failed option is passed
        if self.options.failed:
            self.log("Only doctesting files that failed last test.")

            def is_failure(source):
                basename = source.basename
                return basename not in self.stats or self.stats[basename].get('failed')
            self.sources = [x for x in self.sources if is_failure(x)]

    def sort_sources(self):
        r"""
        This function sorts the sources so that slower doctests are run first.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: dirname = os.path.join(SAGE_SRC, 'sage', 'doctest')
            sage: DD = DocTestDefaults(nthreads=2)
            sage: DC = DocTestController(DD, [dirname])
            sage: DC.expand_files_into_sources()
            sage: DC.sources.sort(key=lambda s:s.basename)
            sage: for i, source in enumerate(DC.sources):
            ....:     DC.stats[source.basename] = {'walltime': 0.1r * (i+1)}
            sage: DC.sort_sources()
            Sorting sources by runtime so that slower doctests are run first....
            sage: print("\n".join(source.basename for source in DC.sources))
            sage.doctest.util
            sage.doctest.test
            sage.doctest.sources
            sage.doctest.rif_tol
            sage.doctest.reporting
            sage.doctest.parsing_test
            sage.doctest.parsing
            sage.doctest.marked_output
            sage.doctest.forker
            sage.doctest.fixtures
            sage.doctest.external
            sage.doctest.control
            sage.doctest.check_tolerance
            sage.doctest.all
            sage.doctest
        """
        if self.options.nthreads > 1 and len(self.sources) > self.options.nthreads:
            self.log("Sorting sources by runtime so that slower doctests are run first....")
            default = {'walltime': 0}

            def sort_key(source):
                basename = source.basename
                return -self.stats.get(basename, default).get('walltime', 0), basename
            self.sources = sorted(self.sources, key=sort_key)

    def source_baseline(self, source):
        r"""
        Return the ``baseline_stats`` value of ``source``.

        INPUT:

        - ``source`` -- a :class:`DocTestSource` instance

        OUTPUT: a dictionary

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: filename = sage.doctest.util.__file__
            sage: DD = DocTestDefaults()
            sage: DC = DocTestController(DD, [filename])
            sage: DC.expand_files_into_sources()
            sage: DC.source_baseline(DC.sources[0])
            {}
        """
        if self.baseline_stats:
            basename = source.basename
            return self.baseline_stats.get(basename, {})
        return {}

    def run_doctests(self):
        """
        Actually run the doctests.

        This function is called by :meth:`run`.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: dirname = os.path.join(SAGE_SRC, 'sage', 'rings', 'homset.py')
            sage: DD = DocTestDefaults()
            sage: DC = DocTestController(DD, [dirname])
            sage: DC.expand_files_into_sources()
            sage: DC.run_doctests()
            Doctesting 1 file.
            sage -t .../sage/rings/homset.py
                [... tests, ...s wall]
            ----------------------------------------------------------------------
            All tests passed!
            ----------------------------------------------------------------------
            Total time for all tests: ... seconds
                cpu time: ... seconds
                cumulative wall time: ... seconds...
        """
        nfiles = 0
        nother = 0
        for F in self.sources:
            if isinstance(F, FileDocTestSource):
                nfiles += 1
            else:
                nother += 1
        if self.sources:
            filestr = ", ".join(([count_noun(nfiles, "file")] if nfiles else []) +
                                ([count_noun(nother, "other source")] if nother else []))
            threads = " using %s threads" % (self.options.nthreads) if self.options.nthreads > 1 else ""
            iterations = []
            if self.options.global_iterations > 1:
                iterations.append("%s global iterations" % (self.options.global_iterations))
            if self.options.file_iterations > 1:
                iterations.append("%s file iterations" % (self.options.file_iterations))
            iterations = ", ".join(iterations)
            if iterations:
                iterations = " (%s)" % (iterations)
            if self.baseline_stats:
                self.log(f"Using --baseline-stats-path={self.options.baseline_stats_path}")
            self.log("Doctesting %s%s%s." % (filestr, threads, iterations))
            self.reporter = DocTestReporter(self)
            self.dispatcher = DocTestDispatcher(self)
            N = self.options.global_iterations
            for _ in range(N):
                try:
                    self.timer = Timer().start()
                    self.dispatcher.dispatch()
                except KeyboardInterrupt:
                    break
                finally:
                    self.timer.stop()
                    self.reporter.finalize()
                    self.cleanup(False)
        else:
            self.log("No files to doctest")
            self.reporter = DictAsObject({'error_status': 0, 'stats': {}})

    def cleanup(self, final=True):
        """
        Run cleanup activities after actually running doctests.

        In particular, saves the stats to disk and closes the logfile.

        INPUT:

        - ``final`` -- whether to close the logfile

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: dirname = os.path.join(SAGE_SRC, 'sage', 'rings', 'all.py')
            sage: DD = DocTestDefaults()

            sage: DC = DocTestController(DD, [dirname])
            sage: DC.expand_files_into_sources()
            sage: DC.sources.sort(key=lambda s:s.basename)

            sage: for i, source in enumerate(DC.sources):
            ....:     DC.stats[source.basename] = {'walltime': 0.1r * (i+1)}
            ....:

            sage: DC.run()
            Running doctests with ID ...
            Doctesting 1 file.
            sage -t .../rings/all.py
                [... tests, ...s wall]
            ----------------------------------------------------------------------
            All tests passed!
            ----------------------------------------------------------------------
            Total time for all tests: ... seconds
                cpu time: ... seconds
                cumulative wall time: ... seconds
            Features detected...
            0
            sage: DC.cleanup()
        """
        self.stats.update(self.reporter.stats)
        self.save_stats(self.options.stats_path)
        # Close the logfile
        if final and self.logfile is not None:
            self.logfile.close()
            self.logfile = None

    def _optional_tags_string(self):
        """
        Return a string describing the optional tags used.

        OUTPUT: string with comma-separated tags (without spaces, so
        it can be used to build a command-line)

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: DC = DocTestController(DocTestDefaults(), [])
            sage: DC._optional_tags_string()
            'sage'
            sage: DC = DocTestController(DocTestDefaults(optional='all,and,some,more'), [])
            sage: DC._optional_tags_string()
            'all'
            sage: DC = DocTestController(DocTestDefaults(optional='sage,openssl'), [])
            sage: DC._optional_tags_string()
            'openssl,sage'
        """
        tags = self.options.optional
        if tags is True:
            return "all"
        else:
            return ",".join(sorted(tags - auto_optional_tags))

    def _assemble_cmd(self):
        """
        Assemble a shell command used in running tests under gdb, lldb, or valgrind.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: DC = DocTestController(DocTestDefaults(timeout=123), ["hello_world.py"])
            sage: print(DC._assemble_cmd())
            ...python... -m sage.doctest --serial... --timeout=123... hello_world.py
        """
        cmd = f"{shlex.quote(sys.executable)} -m sage.doctest --serial "
        opt = dict_difference(self.options.__dict__, DocTestDefaults(runtest_default=True).__dict__)
        # Options with no argument
        for o in ("all", "installed", "long", "initial", "exitfirst",
                  "force_lib", "if_installed", "abspath", "verbose",
                  "debug", "only_errors", "failed", "new",
                  "show_skipped"):
            if o in opt:
                cmd += "--%s " % o.replace('_', '-')
        # Options with one argument
        for o in ("timeout", "die_timeout", "logfile", "warn_long", "randorder",
                  "random_seed", "global_iterations", "file_iterations",
                  "environment", "baseline_stats_path", "stats_path"):
            if o in opt:
                cmd += "--%s=%s " % (o.replace('_', '-'), opt[o])
        # One with a different dest
        if "target_walltime" in opt:
            cmd += "--%s=%s " % ("short", opt[o])
        if "optional" in opt:
            cmd += "--optional={} ".format(self._optional_tags_string())
        return cmd + " ".join(self.files)

    def run_val_gdb(self, testing=False):
        """
        Spawns a subprocess to run tests under the control of gdb, lldb, or valgrind.

        INPUT:

        - ``testing`` -- boolean (default: ``False``); if ``True`` then the
          command to be run will be printed rather than a subprocess started

        EXAMPLES:

        Note that the command lines include unexpanded environment
        variables. It is safer to let the shell expand them than to
        expand them here and risk insufficient quoting. ::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: DD = DocTestDefaults(gdb=True)
            sage: DC = DocTestController(DD, ["hello_world.py"])
            sage: DC.run_val_gdb(testing=True)
            exec gdb --eval-command="run" --args ...python... -m sage.doctest --serial... --timeout=0... hello_world.py

        ::

            sage: DD = DocTestDefaults(valgrind=True, optional='all', timeout=172800)
            sage: DC = DocTestController(DD, ["hello_world.py"])
            sage: DC.run_val_gdb(testing=True)
            exec valgrind --tool=memcheck --leak-resolution=high --leak-check=full --num-callers=25 --suppressions=.../valgrind/pyalloc.supp --suppressions=.../valgrind/sage.supp --suppressions=.../valgrind/sage-additional.supp --suppressions=.../valgrind/valgrind-python.supp  --log-file=.../valgrind/sage-memcheck.%p ...python... -m sage.doctest --serial... --timeout=172800... --optional=all hello_world.py
        """
        try:
            sage_cmd = self._assemble_cmd()
        except ValueError:
            self.log(sys.exc_info()[1])
            return 2
        opt = self.options

        if opt.gdb:
            cmd = f'''exec gdb --eval-command="run" --args '''
            flags = ""
            if opt.logfile:
                sage_cmd += f" --logfile {shlex.quote(opt.logfile)}"
        elif opt.lldb:
            cmd = f'''exec lldb --one-line "process launch" --one-line "cont" -- '''
            flags = ""
        else:
            if opt.logfile is None:
                default_log = os.path.join(DOT_SAGE, "valgrind")
                if os.path.exists(default_log):
                    if not os.path.isdir(default_log):
                        self.log(f"{default_log} must be a directory")
                        return 2
                else:
                    os.makedirs(default_log)
                logfile = os.path.join(default_log, "sage-%s")
            else:
                logfile = opt.logfile
            if opt.valgrind:
                toolname = "memcheck"
                flags = os.getenv("SAGE_MEMCHECK_FLAGS")
                if flags is None:
                    flags = "--leak-resolution=high --leak-check=full --num-callers=25 "
                    for supp in ["pyalloc.supp", "sage.supp", "sage-additional.supp", "valgrind-python.supp"]:
                        fname = os.path.join(SAGE_EXTCODE, "valgrind", supp)
                        flags += f"--suppressions={shlex.quote(fname)} "
            elif opt.massif:
                toolname = "massif"
                flags = os.getenv("SAGE_MASSIF_FLAGS", "--depth=6 ")
            elif opt.cachegrind:
                toolname = "cachegrind"
                flags = os.getenv("SAGE_CACHEGRIND_FLAGS", "")
            elif opt.omega:
                toolname = "exp-omega"
                flags = os.getenv("SAGE_OMEGA_FLAGS", "")
            cmd = "exec valgrind --tool=%s " % (toolname)
            flags += f''' --log-file={shlex.quote(logfile)} '''
            if opt.omega:
                toolname = "omega"
            if "%s" in flags:
                flags %= toolname + ".%p"  # replace %s with toolname
        cmd += flags + sage_cmd

        sys.stdout.flush()
        sys.stderr.flush()
        self.log(cmd)

        if testing:
            return

        # Setup signal handlers.
        # Save crash logs in temporary directory.
        os.putenv('CYSIGNALS_CRASH_LOGS', tmp_dir("crash_logs_"))
        init_cysignals()

        import signal
        import subprocess
        p = subprocess.Popen(cmd, shell=True)

        if opt.timeout > 0:
            signal.alarm(opt.timeout)
        try:
            return p.wait()
        except AlarmInterrupt:
            self.log("    Timed out")
            return 4
        except KeyboardInterrupt:
            self.log("    Interrupted")
            return 128
        finally:
            signal.alarm(0)
            if p.returncode is None:
                p.terminate()

    def run(self):
        """
        This function is called after initialization to set up and run all doctests.

        EXAMPLES::

            sage: from sage.doctest.control import DocTestDefaults, DocTestController
            sage: from sage.env import SAGE_SRC
            sage: import os
            sage: DD = DocTestDefaults()
            sage: filename = os.path.join(SAGE_SRC, "sage", "sets", "non_negative_integers.py")
            sage: DC = DocTestController(DD, [filename])
            sage: DC.run()
            Running doctests with ID ...
            Doctesting 1 file.
            sage -t .../sage/sets/non_negative_integers.py
                [... tests, ...s wall]
            ----------------------------------------------------------------------
            All tests passed!
            ----------------------------------------------------------------------
            Total time for all tests: ... seconds
                cpu time: ... seconds
                cumulative wall time: ... seconds
            Features detected...
            0

        We check that :issue:`25378` is fixed (testing external packages
        while providing a logfile does not raise a ValueError: I/O
        operation on closed file)::

            sage: logfile = tmp_filename(ext='.log')
            sage: DD = DocTestDefaults(optional=set(['sage', 'external']), logfile=logfile)
            sage: filename = tmp_filename(ext='.py')
            sage: DC = DocTestController(DD, [filename])
            sage: DC.run()
            Running doctests with ID ...
            Using --optional=external,sage
            Features to be detected: ...
            Doctesting 1 file.
            sage -t ....py
                [0 tests, ...s wall]
            ----------------------------------------------------------------------
            All tests passed!
            ----------------------------------------------------------------------
            Total time for all tests: ... seconds
                cpu time: ... seconds
                cumulative wall time: ... seconds
            Features detected...
            0

        We test the ``--hide`` option (:issue:`34185`)::

            sage: from sage.doctest.control import test_hide
            sage: filename = tmp_filename(ext='.py')
            sage: with open(filename, 'w') as f:
            ....:     f.write(test_hide)
            ....:     f.close()
            714
            sage: DF = DocTestDefaults(hide='buckygen,all')
            sage: DC = DocTestController(DF, [filename])
            sage: DC.run()
            Running doctests with ID ...
            Using --optional=sage...
            Features to be detected: ...
            Doctesting 1 file.
            sage -t ....py
                [4 tests, ...s wall]
            ----------------------------------------------------------------------
            All tests passed!
            ----------------------------------------------------------------------
            Total time for all tests: ... seconds
                cpu time: ... seconds
                cumulative wall time: ... seconds
            Features detected...
            0

            sage: DF = DocTestDefaults(hide='benzene,optional')
            sage: DC = DocTestController(DF, [filename])
            sage: DC.run()
            Running doctests with ID ...
            Using --optional=sage
            Features to be detected: ...
            Doctesting 1 file.
            sage -t ....py
                [4 tests, ...s wall]
            ----------------------------------------------------------------------
            All tests passed!
            ----------------------------------------------------------------------
            Total time for all tests: ... seconds
                cpu time: ... seconds
                cumulative wall time: ... seconds
            Features detected...
            0

        Test *Features that have been hidden* message::

            sage: DC.run()                              # optional - meataxe
            Running doctests with ID ...
            Using --optional=sage
            Features to be detected: ...
            Doctesting 1 file.
            sage -t ....py
                [4 tests, ...s wall]
            ----------------------------------------------------------------------
            All tests passed!
            ----------------------------------------------------------------------
            Total time for all tests: ... seconds
                cpu time: ... seconds
                cumulative wall time: ... seconds
            Features detected...
            Features that have been hidden: ...meataxe...
            0
        """
        opt = self.options
        L = (opt.gdb, opt.lldb, opt.valgrind, opt.massif, opt.cachegrind, opt.omega)
        if any(L):
            if L.count(True) > 1:
                self.log("You may only specify one of gdb, valgrind/memcheck, massif, cachegrind, omega")
                return 2
            return self.run_val_gdb()
        else:
            self.create_run_id()
            from sage.env import SAGE_ROOT_GIT, SAGE_LOCAL, SAGE_VENV
            # SAGE_ROOT_GIT can be None on distributions which typically
            # only have the SAGE_LOCAL install tree but not SAGE_ROOT
            if (SAGE_ROOT_GIT is not None) and os.path.isdir(SAGE_ROOT_GIT):
                import subprocess
                try:
                    branch = subprocess.check_output(["git",
                                                      "--git-dir=" + SAGE_ROOT_GIT,
                                                      "rev-parse",
                                                      "--abbrev-ref",
                                                      "HEAD"])
                    branch = branch.decode('utf-8')
                    self.log("Git branch: " + branch, end="")
                except subprocess.CalledProcessError:
                    pass
                try:
                    ref = subprocess.check_output(["git",
                                                   "--git-dir=" + SAGE_ROOT_GIT,
                                                   "describe",
                                                   "--always",
                                                   "--dirty"])
                    ref = ref.decode('utf-8')
                    self.log("Git ref: " + ref, end="")
                except subprocess.CalledProcessError:
                    pass

            self.log(f"Running with {SAGE_LOCAL=} and {SAGE_VENV=}")

            self.log("Using --optional=" + self._optional_tags_string())
            available_software._allow_external = self.options.optional is True or 'external' in self.options.optional

            for h in self.options.hide:
                try:
                    i = available_software._indices[h]
                except KeyError:
                    pass
                else:
                    f = available_software._features[i]
                    f.hide()
                    self.options.hidden_features.add(f)
                    for g in f.joined_features():
                        if g.name in self.options.optional:
                            self.options.optional.discard(g.name)

            for o in self.options.disabled_optional:
                try:
                    i = available_software._indices[o]
                except KeyError:
                    pass
                else:
                    available_software._seen[i] = -1

            self.log("Features to be detected: " + ','.join(available_software.detectable()))
            if self.options.probe:
                self.log("Features to be probed: " + ('all' if self.options.probe is True
                                                      else ','.join(self.options.probe)))
            self.add_files()
            self.expand_files_into_sources()
            self.filter_sources()
            self.sort_sources()
            self.run_doctests()

            self.log("Features detected for doctesting: "
                     + ','.join(available_software.seen()))
            if self.options.hidden_features:
                for f in self.options.hidden_features:
                    f.unhide()
                self.log("Features that have been hidden: " + ','.join(available_software.hidden()))
            self.cleanup()
            return self.reporter.error_status


def run_doctests(module, options=None):
    """
    Run the doctests in a given file.

    INPUT:

    - ``module`` -- a Sage module, a string, or a list of such

    - ``options`` -- a DocTestDefaults object or ``None``

    EXAMPLES::

        sage: run_doctests(sage.rings.all)
        Running doctests with ID ...
        Doctesting 1 file.
        sage -t .../sage/rings/all.py
            [... tests, ...s wall]
        ----------------------------------------------------------------------
        All tests passed!
        ----------------------------------------------------------------------
        Total time for all tests: ... seconds
            cpu time: ... seconds
            cumulative wall time: ... seconds
        Features detected...
    """
    import sys
    sys.stdout.flush()

    def stringify(x):
        if isinstance(x, (list, tuple)):
            F = [stringify(a) for a in x]
            return sage.misc.flatten.flatten(F)
        elif isinstance(x, types.ModuleType):
            F = x.__file__.replace(SAGE_LIB, SAGE_SRC)
            base, pyfile = os.path.split(F)
            file, ext = os.path.splitext(pyfile)
            if ext == ".pyc":
                ext = ".py"
            elif ext == ".so":
                ext = ".pyx"
            if file == "__init__":
                return [base]
            else:
                return [os.path.join(base, file) + ext]
        elif isinstance(x, str):
            return [os.path.abspath(x)]
    F = stringify(module)
    if options is None:
        options = DocTestDefaults()
    DC = DocTestController(options, F)

    # Determine whether we're in doctest mode
    save_dtmode = sage.doctest.DOCTEST_MODE

    # We need the following if we're not in DOCTEST_MODE
    # Tell IPython to avoid colors: it screws up the output checking.
    if not save_dtmode:
        if options.debug:
            raise ValueError("You should not try to run doctests with a debugger from within Sage: IPython objects to embedded shells")
        from IPython.core.getipython import get_ipython
        IP = get_ipython()
        if IP is not None:
            old_color = IP.colors
            IP.run_line_magic('colors', 'NoColor')
            old_config_color = IP.config.TerminalInteractiveShell.colors
            IP.config.TerminalInteractiveShell.colors = 'NoColor'

    try:
        DC.run()
    finally:
        sage.doctest.DOCTEST_MODE = save_dtmode
        if not save_dtmode and IP is not None:
            IP.run_line_magic('colors', old_color)
            IP.config.TerminalInteractiveShell.colors = old_config_color


###############################################################################
# Declaration of doctest strings
###############################################################################

test_hide = r"""r{quotmark}
{prompt}: next(graphs.fullerenes(20))
Traceback (most recent call last):
 ...
FeatureNotPresentError: buckygen is not available.
...
{prompt}: next(graphs.fullerenes(20))   # optional - buckygen
Graph on 20 vertices

{prompt}: len(list(graphs.fusenes(2)))
Traceback (most recent call last):
 ...
FeatureNotPresentError: benzene is not available.
...
{prompt}: len(list(graphs.fusenes(2)))  # optional - benzene
1
{prompt}: from sage.matrix.matrix_space import get_matrix_class
{prompt}: get_matrix_class(GF(25,'x'), 4, 4, False, 'meataxe')
Failed lazy import:
meataxe is not available.
...
{prompt}: get_matrix_class(GF(25,'x'), 4, 4, False, 'meataxe')  # optional - meataxe
<class 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>
{quotmark}
""".format(quotmark='"""', prompt='sage')  # using prompt to hide these lines from _test_enough_doctests
