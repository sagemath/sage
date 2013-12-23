Invoking Sage
=============

To run Sage, you basically just need to type ``sage`` from the
command-line prompt to start the Sage interpreter.  See the Sage
Installation Guide for information about making sure your
:envvar:`$PATH` is set correctly, etc.

Command-line options for Sage
-----------------------------

.. rubric:: Running Sage, the most common options

- ``file.[sage|py|spyx]`` -- run the given .sage, .py or .spyx
  files (as in ``sage my_file.sage``)
- ``-h``, ``-?``, ``--help`` -- print a short help message
- ``-v``, ``--version`` -- print the Sage version
- ``--advanced`` -- print (essentially this) list of Sage options
- ``-c cmd`` -- evaluate ``cmd`` as sage code.  For example, ``sage
  -c 'print factor(35)'`` will print "5 * 7".

.. rubric:: Running Sage, other options

- ``--preparse file.sage`` -- preparse ``file.sage``, a file of
  Sage code, and produce the corresponding Python file
  ``file.sage.py``.  See the Sage tutorial for more about preparsing
  and the differences between Sage and Python.
- ``-q`` -- quiet; start with no banner
- ``--grep [options] <string>`` -- grep through all the Sage library
  code for ``string``. Any options will get passed to the "grep"
  command; for example, ``sage --grep -i epstein`` will search for
  ``epstein``, and the ``-i`` flag tells grep to ignore case when
  searching. Note that while running Sage, you can also use the
  function :func:`search_src <sage.misc.sagedoc.search_src>` to
  accomplish the same thing.
- ``--grepdoc [options] <string>`` -- grep through all the Sage
  documentation for ``string``.  Note that while running Sage, you can
  also use the function :func:`search_doc
  <sage.misc.sagedoc.search_doc>` to accomplish the same thing.
- ``--min [...]`` -- do not populate global namespace (must be first
  option)
- ``-gthread``, ``-qthread``, ``-q4thread``, ``-wthread``,
  ``-pylab`` -- pass the option through to IPython
- ``--nodotsage`` -- run Sage without using the user's
  :file:`.sage` directory: create and use a temporary :file:`.sage`
  directory instead.  Warning: notebooks are stored in the
  :file:`.sage` directory, so any notebooks created while running with
  ``--nodotsage`` will be temporary also.

.. rubric:: Running the notebook

- ``-n``, ``--notebook`` -- start the Sage notebook, passing all
  remaining arguments to the 'notebook' command in Sage
- ``-bn [...]``, ``--build-and-notebook [...]`` -- build the Sage
  library (as by running ``sage -b``) then start the Sage notebook
- ``--inotebook [...]`` -- start the *insecure* Sage notebook

.. rubric:: Running external programs and utilities

- ``--cython [...]`` -- run Cython with the given arguments
- ``--ecl [...]``, ``--lisp [...]`` -- run Sage's copy of ECL
  (Embeddable Common Lisp) with the given arguments
- ``--gap [...]`` -- run Sage's Gap with the given arguments
- ``--gp [...]`` -- run Sage's PARI/GP calculator with the given arguments
- ``--hg [...]`` -- run Sage's Mercurial with the given arguments
- ``--ipython [...]`` -- run Sage's IPython using the default
  environment (not Sage), passing additional options to IPython
- ``--kash [...]`` -- run Sage's Kash with the given arguments
- ``--M2 [...]`` -- run Sage's Macaulay2 with the given arguments
- ``--maxima [...]`` -- run Sage's Maxima with the given arguments
- ``--mwrank [...]`` -- run Sage's mwrank with the given arguments
- ``--python [...]`` -- run the Python interpreter
- ``-R [...]`` -- run Sage's R with the given arguments
- ``--scons [...]`` -- run Sage's scons
- ``--singular [...]`` -- run Sage's singular with the given arguments
- ``--twistd [...]`` -- run Twisted server
- ``--sh [...]`` -- run a shell with Sage environment variables set
- ``--gdb`` -- run Sage under the control of gdb
- ``--gdb-ipython`` -- run Sage's IPython under the control of gdb
- ``--cleaner`` -- run the Sage cleaner.  This cleans up after Sage,
  removing temporary directories and spawned processes.  (This gets
  run by Sage automatically, so it is usually not necessary to run
  it separately.)

.. rubric:: Installing packages and upgrading

- ``-i [options] [packages]`` -- install the given Sage packages (unless
  they are already installed); if no packages are given, print
  a list of all installed packages. Options:

  - ``-c`` -- run the packages' test suites, overriding the settings of
    :envvar:`SAGE_CHECK` and :envvar:`SAGE_CHECK_PACKAGES`.
  - ``-f`` -- force build: install the packages even if they are
    already installed.
  - ``-s`` -- do not delete the ``spkg/build`` directories after a
    successful build -- useful for debugging.

- ``-f [options] [packages]`` -- shortcut for ``-i -f``: force build of
  the given Sage packages.
- ``--info [packages]`` -- display the ``SPKG.txt`` file of the given
  Sage packages.
- ``--standard`` -- list all standard packages that can be installed
- ``--optional`` -- list all optional packages that can be installed
- ``--experimental`` -- list all experimental packages that can be installed
- ``--upgrade [url]`` -- download, build and install standard
  packages from given url.  If url not given, automatically selects
  a suitable mirror.  If url='ask', it lets you select the mirror.

.. rubric:: The Sage-combinat package manager

Sage-combinat is a collection of experimental patches
(i.e. extensions) on top of Sage, developed by a community of
researchers, with a focus, at least to some extent, in
combinatorics. Many of those patches get eventually integrated into
Sage as soon as they are mature enough, but you can install the
still-experimental ones by running ``sage -combinat install``.  This
creates a new branch, called ``sage-combinat`` by default, containing
the new patches. More information on sage-combinat is available at the
`Sage wiki`__.  More details on the ``--combinat`` command-line option
for Sage:

__ http://wiki.sagemath.org/combinat

- ``--combinat [options] command`` -- run the ``sage-combinat``
  patch management script.  Commands:

  - ``config`` -- show current configuration (Sage command, path, version, ...)
  - ``install`` -- install the sage-combinat patches
  - ``update`` -- update to the latest sage-combinat patches
  - ``upgrade`` -- upgrade Sage and update to the latest sage-combinat patches
  - ``status`` -- show changed files in the working directory and in
    the patch queue
  - ``qselect`` -- choose appropriate guards for the current version of Sage

  Options:

  - ``-h``, ``--help`` -- print a help message
  - ``-b BRANCH``, ``--branch=BRANCH`` -- use ``sage-BRANCH``
    instead of ``sage-combinat``
  - ``--sage=/opt/bin/sage`` -- specify the path to Sage
  - ``-f``, ``--force`` -- force proceeding, skipping any relevant queries
  - ``-v`` -- Be verbose; print status messages
  - ``-q``, ``--quiet`` -- don't print status messages
  - ``-s URL``, ``--server=URL`` -- set the URL for the
    sage-combinat server; the default is
    ``http://combinat.sagemath.org/patches``
  - ``-n`` -- after qselect: disable all previous non version guards

.. rubric:: Building and testing the Sage library

- ``--root`` -- print the Sage root directory
- ``--branch`` -- print the current Sage branch
- ``--clone [new branch]`` -- clone a new branch of the Sage library from the
  current branch
- ``-b [branch]`` -- build Sage library -- do this if you have modified
  any source code files in :file:`$SAGE_ROOT/devel/sage/`.  If
  ``branch`` is given, switch to the branch in
  :file:`$SAGE_ROOT/devel/sage-branch` and build it.
- ``-ba [branch]`` -- same as ``-b``, but rebuild *all* Cython
  code.  This could take a while, so you will be asked if you want
  to proceed.
- ``-ba-force [branch]`` -- same as ``-ba``, but don't query before
  rebuilding
- ``--br [branch]`` -- switch to, build, and run Sage with the given
  branch
- ``-t [options] <files|dir>`` -- test examples in .py, .pyx, .sage
  or .tex files.  Options:

  - ``--long``  -- include lines with the phrase 'long time'
  - ``--verbose`` -- print debugging output during the test
  - ``--optional`` -- also test all examples labeled ``# optional``
  - ``--only-optional[=tags]`` -- if no ``tags`` are specified, only
    run blocks of tests containing a line labeled ``# optional``. If
    a comma separated list of tags is specified, only run blocks containing
    a line labeled ``# optional tag`` for any of the tags given and in these blocks only
    run the lines which are unlabeled or labeled ``#optional`` or labeled
    ``#optional tag`` for any of the tags given.
  - ``--randorder[=seed]`` -- randomize order of tests

- ``-tnew [...]`` -- like ``-t`` above, but only tests files
  modified since last commit
- ``-tp <N> [...]`` -- like ``-t`` above, but tests in parallel
  using ``N`` threads with 0 interpreted as ``minimum(8, cpu_count())``
- ``--testall [options]`` -- test all source files, docs, and
  examples; options are the same as for ``-t``.
- ``-bt [...]`` -- build and test, options like ``-t`` above
- ``-btp <N> [...]`` -- build and test in parallel, options like
  ``-tp`` above
- ``-btnew [...]`` -- build and test modified files, options like ``-tnew``
- ``--fixdoctests file.py [output_file] [--long]`` -- writes a new
  version of ``file.py`` to ``output_file`` (default: ``file.py.out``)
  that will pass the doctests. With the optional ``--long`` argument
  the long time tests are also checked. A patch for the new file is
  printed to stdout.
- ``--startuptime [module]`` -- display how long each component of Sage takes
  to start up. Optionally specify a module (e.g., "sage.rings.qqbar") to get
  more details about that particular module.
- ``--coverage <files>`` -- give information about doctest coverage
  of files
- ``--coverageall`` -- give summary info about doctest coverage of
  all files in the Sage library
- ``--sync-build`` -- delete any files in :file:`$SAGE_ROOT/devel/sage/build/`
  which don't have a corresponding source file in
  :file:`$SAGE_ROOT/devel/sage/sage/`

.. rubric:: Documentation

- ``--docbuild [options] document (format | command)`` -- build or
  return information about the Sage documentation.

  - ``document`` -- name of the document to build
  - ``format`` -- document output format
  - ``command`` -- document-specific command

  A ``document`` and either a ``format`` or a ``command`` are required, unless a
  list of one or more of these is requested.

  Options:

  - ``help``, ``-h``, ``--help`` -- print a help message
  - ``-H``, ``--help-all`` -- print an extended help message,
    including the output from the options ``-h``, ``-D``, ``-F``,
    ``-C all``, and a short list of examples.
  - ``-D``, ``--documents`` -- list all available documents
  - ``-F``, ``--formats`` -- list all output formats
  - ``-C DOC``, ``--commands=DOC`` -- list all commands for document
    ``DOC``; use ``-C all`` to list all
  - ``-i``, ``--inherited`` -- include inherited members in
    reference manual; may be slow, may fail for PDF output
  - ``-u``, ``--underscore`` -- include variables prefixed with
    ``_`` in reference manual; may be slow, may fail for PDF output
  - ``-j``, ``--jsmath`` -- render math using jsMath; formats:
    ``html``, ``json``, ``pickle``, ``web``
  - ``--no-pdf-links`` -- do not include PDF links in document
    ``website``; formats: ``html``, ``json``, ``pickle``, ``web``
  - ``--check-nested`` -- check picklability of nested classes in
    document ``reference``
  - ``-N``, ``--no-colors`` -- do not color output; does not affect
    children
  - ``-q``, ``--quiet`` -- work quietly; same as ``--verbose=0``
  - ``-v LEVEL``, ``--verbose=LEVEL`` -- report progress at level 0
    (quiet), 1 (normal), 2 (info), or 3 (debug); does not affect
    children

  Advanced -- use these options with care:

  - ``-S OPTS``, ``--sphinx-opts=OPTS`` -- pass comma-separated ``OPTS``
    to sphinx-build
  - ``-U``, ``--update-mtimes`` -- before building reference manual,
    update modification times for auto-generated ReST files

.. rubric:: Making Sage packages or distributions

- ``--pkg dir`` -- create the Sage package ``dir.spkg`` from the
  directory ``dir``
- ``--pkg_nc dir`` -- as ``--pkg``, but do not compress the package
- ``--merge`` -- run Sage's automatic merge and test script
- ``--bdist VER`` -- build a binary distribution of Sage, with
  version ``VER``
- ``--sdist`` -- build a source distribution of Sage
- ``--crap sage-ver.tar`` -- detect suspicious garbage in the Sage
  source tarball

.. rubric:: Valgrind memory debugging

- ``--cachegrind`` -- run Sage using Valgrind's cachegrind tool
- ``--callgrind`` -- run Sage using Valgrind's callgrind tool
- ``--massif`` -- run Sage using Valgrind's massif tool
- ``--memcheck`` -- run Sage using Valgrind's memcheck tool
- ``--omega`` -- run Sage using Valgrind's omega tool
- ``--valgrind`` -- this is an alias for ``--memcheck``
