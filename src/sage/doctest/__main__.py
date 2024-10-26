import argparse
import os
import sys

# Note: the DOT_SAGE and SAGE_STARTUP_FILE environment variables have already been set by sage-env
DOT_SAGE = os.environ.get('DOT_SAGE', os.path.join(os.environ.get('HOME'),
                                                   '.sage'))

# Override to not pick up user configuration, see Issue #20270
os.environ['SAGE_STARTUP_FILE'] = os.path.join(DOT_SAGE, 'init-doctests.sage')


def _get_optional_defaults():
    """Return the default value for the --optional flag."""
    optional = ['sage', 'optional']

    return ','.join(optional)


def _make_parser():
    r"""
    Return the :class:`argparse.ArgumentParser`.

    TESTS:

    Test that the defaults are the consistent::

        sage: from sage.doctest.control import DocTestDefaults
        sage: from sage.doctest.__main__ import _make_parser
        sage: os.environ.pop('SAGE_DOCTEST_RANDOM_SEED', None)
        ...
        sage: parser = _make_parser()
        sage: args = parser.parse_args([])
        sage: DD = DocTestDefaults(runtest_default=True); DD
        DocTestDefaults(abspath=False, file_iterations=0, global_iterations=0,
                        optional='sage,optional', random_seed=None,
                        stats_path='.../timings2.json')
        sage: D = copy(args.__dict__)
        sage: del D['filenames']
        sage: DA = DocTestDefaults(runtest_default=True, **D); DA
        DocTestDefaults(abspath=False, file_iterations=0, global_iterations=0,
                        optional='sage,optional', random_seed=None,
                        stats_path='.../timings2.json')
    """
    parser = argparse.ArgumentParser(usage="sage -t [options] filenames",
                                     description="Run all tests in a file or a list of files whose extensions "
                                                 "are one of the following: "
                                                 ".py, .pyx, .pxd, .pxi, .sage, .spyx, .tex, .rst.")
    parser.add_argument("-p", "--nthreads", dest="nthreads",
                        type=int, nargs='?', const=0, default=1, metavar="N",
                        help="test in parallel using N threads, with 0 interpreted as max(2, min(8, cpu_count())); "
                        "when run under the control of the GNU make jobserver (make -j), request as most N job slots")
    parser.add_argument("-T", "--timeout", type=int, default=-1, help="timeout (in seconds) for doctesting one file, 0 for no timeout")
    what = parser.add_mutually_exclusive_group()
    what.add_argument("-a", "--all", action="store_true", default=False, help="test all files in the Sage library")
    what.add_argument("--installed", action="store_true", default=False, help="test all installed modules of the Sage library")
    parser.add_argument("--logfile", type=argparse.FileType('a'), metavar="FILE", help="log all output to FILE")

    parser.add_argument("--format", choices=["sage", "github"], default="sage",
                        help="set format of error messages and warnings")
    parser.add_argument("-l", "--long", action="store_true", default=False, help="include lines with the phrase 'long time'")
    parser.add_argument("-s", "--short", dest="target_walltime", nargs='?',
                        type=int, default=-1, const=300, metavar="SECONDS",
                        help="run as many doctests as possible in about 300 seconds (or the number of seconds given as an optional argument)")
    parser.add_argument("--warn-long", dest="warn_long", nargs='?',
                        type=float, default=-1.0, const=1.0, metavar="SECONDS",
                        help="warn if tests take more CPU time than SECONDS")
    # By default, include all tests marked 'sagemath_doc_html' -- see
    # https://github.com/sagemath/sage/issues/25345 and
    # https://github.com/sagemath/sage/issues/26110:
    parser.add_argument("--optional", metavar="FEATURES", default=_get_optional_defaults(),
        help='only run tests including one of the "# optional" tags listed in FEATURES (separated by commas); '
             'if "sage" is listed, will also run the standard doctests; '
             'if "sagemath_doc_html" is listed, will also run the tests relying on the HTML documentation; '
             'if "optional" is listed, will also run tests for installed optional packages or detected features; '
             'if "external" is listed, will also run tests for available external software; '
             'if set to "all", then all tests will be run; '
             'use "!FEATURE" to disable tests marked "# optional - FEATURE". '
             'Note that "!" needs to be quoted or escaped in the shell.')
    parser.add_argument("--hide", metavar="FEATURES", default="",
        help='run tests pretending that the software listed in FEATURES (separated by commas) is not installed; '
             'if "all" is listed, will also hide features corresponding to all optional or experimental packages; '
             'if "optional" is listed, will also hide features corresponding to optional packages.')
    parser.add_argument("--probe", metavar="FEATURES", default="",
        help='run tests that would not be run because one of the given FEATURES (separated by commas) is not installed; '
             'report the tests that pass nevertheless')
    parser.add_argument("--randorder", type=int, metavar="SEED", help="randomize order of tests")
    parser.add_argument("--random-seed", dest="random_seed", type=int, metavar="SEED", help="random seed (integer) for fuzzing doctests",
                        default=os.environ.get("SAGE_DOCTEST_RANDOM_SEED"))
    parser.add_argument("--global-iterations", "--global_iterations", type=int, default=0, help="repeat the whole testing process this many times")
    parser.add_argument("--file-iterations", "--file_iterations", type=int, default=0, help="repeat each file this many times, stopping on the first failure")
    parser.add_argument("--environment", type=str, default="sage.repl.ipython_kernel.all_jupyter", help="name of a module that provides the global environment for tests")

    parser.add_argument("-i", "--initial", action="store_true", default=False, help="only show the first failure in each file")
    parser.add_argument("--exitfirst", action="store_true", default=False, help="end the test run immediately after the first failure or unexpected exception")
    parser.add_argument("--force_lib", "--force-lib", action="store_true", default=False, help="do not import anything from the tested file(s)")
    parser.add_argument("--if-installed", action="store_true", default=False, help="skip Python/Cython files that are not installed as modules")
    parser.add_argument("--abspath", action="store_true", default=False, help="print absolute paths rather than relative paths")
    parser.add_argument("--verbose", action="store_true", default=False, help="print debugging output during the test")
    parser.add_argument("-d", "--debug", action="store_true", default=False, help="drop into a python debugger when an unexpected error is raised")
    parser.add_argument("--only-errors", action="store_true", default=False, help="only output failures, not test successes")

    parser.add_argument("--gdb", action="store_true", default=False, help="run doctests under the control of gdb")
    parser.add_argument("--lldb", action="store_true", default=False, help="run doctests under the control of lldb")
    parser.add_argument("--valgrind", "--memcheck", action="store_true", default=False,
                        help="run doctests using Valgrind's memcheck tool.  The log "
                        "files are named sage-memcheck.PID and can be found in " +
                        os.path.join(DOT_SAGE, "valgrind"))
    parser.add_argument("--massif", action="store_true", default=False,
                        help="run doctests using Valgrind's massif tool.  The log "
                        "files are named sage-massif.PID and can be found in " +
                        os.path.join(DOT_SAGE, "valgrind"))
    parser.add_argument("--cachegrind", action="store_true", default=False,
                        help="run doctests using Valgrind's cachegrind tool.  The log "
                        "files are named sage-cachegrind.PID and can be found in " +
                        os.path.join(DOT_SAGE, "valgrind"))
    parser.add_argument("--omega", action="store_true", default=False,
                        help="run doctests using Valgrind's omega tool.  The log "
                        "files are named sage-omega.PID and can be found in " +
                        os.path.join(DOT_SAGE, "valgrind"))

    parser.add_argument("-f", "--failed", action="store_true", default=False,
        help="doctest only those files that failed in the previous run")
    what.add_argument("-n", "--new", action="store_true", default=False,
        help="doctest only those files that have been changed in the repository and not yet been committed")
    parser.add_argument("--show-skipped", "--show_skipped", action="store_true", default=False,
        help="print a summary at the end of each file of optional tests that were skipped")

    parser.add_argument("--stats_path", "--stats-path", default=os.path.join(DOT_SAGE, "timings2.json"),
                        help="path to a json dictionary for timings and failure status for each file from previous runs; it will be updated in this run")
    parser.add_argument("--baseline_stats_path", "--baseline-stats-path", default=None,
                        help="path to a json dictionary for timings and failure status for each file, to be used as a baseline; it will not be updated")

    class GCAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            gcopts = dict(DEFAULT=0, ALWAYS=1, NEVER=-1)
            new_value = gcopts[values]
            setattr(namespace, self.dest, new_value)

    parser.add_argument("--gc",
                        choices=["DEFAULT", "ALWAYS", "NEVER"],
                        default=0,
                        action=GCAction,
                        help="control garbarge collection "
                        "(ALWAYS: collect garbage before every test; NEVER: disable gc; DEFAULT: Python default)")

    # The --serial option is only really for internal use, better not
    # document it.
    parser.add_argument("--serial", action="store_true", default=False, help=argparse.SUPPRESS)
    # Same for --die_timeout
    parser.add_argument("--die_timeout", type=int, default=-1, help=argparse.SUPPRESS)

    parser.add_argument("filenames", help="file names", nargs='*')
    return parser


def main():
    parser = _make_parser()
    # custom treatment to separate properly
    # one or several file names at the end
    new_arguments = []
    need_filenames = True
    in_filenames = False
    afterlog = False
    for arg in sys.argv[1:]:
        if arg in ('-n', '--new', '-a', '--all', '--installed'):
            need_filenames = False
        elif need_filenames and not (afterlog or in_filenames) and os.path.exists(arg):
            in_filenames = True
            new_arguments.append('--')
        new_arguments.append(arg)
        afterlog = arg in ['--logfile', '--stats_path', '--stats-path',
                           '--baseline_stats_path', '--baseline-stats-path']

    args = parser.parse_args(new_arguments)

    if not args.filenames and not (args.all or args.new or args.installed):
        print('either use --new, --all, --installed, or some filenames')
        return 2

    # Limit the number of threads to 2 to save system resources.
    # See Issue #23713, #23892, #30351
    if sys.platform == 'darwin':
        os.environ["OMP_NUM_THREADS"] = "1"
    else:
        os.environ["OMP_NUM_THREADS"] = "2"

    os.environ["SAGE_NUM_THREADS"] = "2"

    from sage.doctest.control import DocTestController
    DC = DocTestController(args, args.filenames)
    err = DC.run()

    # Issue #33521: Do not run pytest if the pytest configuration is not available.
    # This happens when the source tree is not available and SAGE_SRC falls back
    # to SAGE_LIB.
    from sage.env import SAGE_SRC
    if not all(os.path.isfile(os.path.join(SAGE_SRC, f))
               for f in ["conftest.py", "tox.ini"]):
        return err

    try:
        exit_code_pytest = 0
        import pytest
        pytest_options = []
        if args.verbose:
            pytest_options.append("-v")

        # #35999: no filename in arguments defaults to "src"
        if not args.filenames:
            filenames = [SAGE_SRC]
        else:
            # #31924: Do not run pytest on individual Python files unless
            # they match the pytest file pattern.  However, pass names
            # of directories. We use 'not os.path.isfile(f)' for this so that
            # we do not silently hide typos.
            filenames = [f for f in args.filenames
                         if f.endswith("_test.py") or not os.path.isfile(f)]
        if filenames:
            print(f"Running pytest on {filenames} with options {pytest_options}")
            exit_code_pytest = pytest.main(filenames + pytest_options)
            if exit_code_pytest == 5:
                # Exit code 5 means there were no test files, pass in this case
                exit_code_pytest = 0

    except ModuleNotFoundError:
        print("pytest is not installed in the venv, skip checking tests that rely on it")

    if err == 0:
        return exit_code_pytest
    return err


if __name__ == "__main__":
    sys.exit(main())
