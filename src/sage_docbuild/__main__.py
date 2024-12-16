# sage.doctest: needs sphinx
r"""
Sage docbuild main

This module defines the Sage documentation build command::

    sage --docbuild [OPTIONS] DOCUMENT (FORMAT | COMMAND)

If ``FORMAT`` is given, it builds ``DOCUMENT`` in ``FORMAT``. If ``COMMAND`` is
given, it returns information about ``DOCUMENT``.

Run ``sage --docbuild`` to get detailed explanations about
arguments and options.

Positional arguments::

  DOCUMENT              name of the document to build. It can be either one of
                        the documents listed by -D or 'file=/path/to/FILE' to
                        build documentation for this specific file.
  FORMAT or COMMAND     document output format (or command)

Standard options::

  -h, --help            show a help message and exit
  -H, --help-all        show an extended help message and exit
  -D, --documents       list all available DOCUMENTs
  -F, --formats         list all output FORMATs
  -C DOC, --commands DOC
                        list all COMMANDs for DOCUMENT DOC; use 'all' to list all
  -i, --inherited       include inherited members in reference manual; may be
                        slow, may fail for PDF output
  -u, --underscore      include variables prefixed with '_' in reference
                        manual; may be slow, may fail for PDF output
  -j, --mathjax, --jsmath
                        ignored for backwards compatibility
  --no-plot             do not include graphics auto-generated using the '.. plot' markup
  --no-preparsed-examples
                        do not show preparsed versions of EXAMPLES blocks
  --include-tests-blocks
                        include TESTS blocks in the reference manual
  --no-pdf-links        do not include PDF links in DOCUMENT 'website';
                        FORMATs: html, json, pickle, web
  --live-doc            make Sage code blocks live for html FORMAT
  --warn-links          issue a warning whenever a link is not properly
                        resolved; equivalent to '--sphinx-opts -n' (sphinx
                        option: nitpicky)
  --check-nested        check picklability of nested classes in DOCUMENT 'reference'
  --no-prune-empty-dirs
                        do not prune empty directories in the documentation source
  --use-cdns            assume internet connection and use CDNs; in particular,
                        use MathJax CDN
  -N, --no-colors       do not color output; does not affect children
  -q, --quiet           work quietly; same as --verbose=0
  -v LEVEL, --verbose LEVEL
                        report progress at LEVEL=0 (quiet), 1 (normal), 2
                        (info), or 3 (debug); does not affect children
  -o DIR, --output DIR  if DOCUMENT is a single file ('file=...'), write output
                        to this directory

Advanced options::

  Use these options with care.

  -S OPTS, --sphinx-opts OPTS
                        pass comma-separated OPTS to sphinx-build; must precede
                        OPTS with '=', as in '-S=-q,-aE' or '-S="-q,-aE"'
  -U, --update-mtimes   before building reference manual, update modification
                        times for auto-generated reST files
  -k, --keep-going      Do not abort on errors but continue as much as possible
                        after an error
  --all-documents ARG   if ARG is 'reference', list all subdocuments of
                        en/reference. If ARG is 'all', list all main documents
"""

import logging
import argparse
import os
import shlex
import sys
import sphinx.ext.intersphinx
from sage.env import SAGE_DOC_SRC
from .builders import DocBuilder, ReferenceBuilder, get_builder, get_documents
from . import build_options

logger = logging.getLogger(__name__)


def format_columns(lst, align='<', cols=None, indent=4, pad=3, width=80):
    """
    Utility function that formats a list as a simple table and returns
    a Unicode string representation.

    The number of columns is
    computed from the other options, unless it's passed as a keyword
    argument.  For help on Python's string formatter, see

    https://docs.python.org/library/string.html#format-string-syntax
    """
    # Can we generalize this (efficiently) to other / multiple inputs
    # and generators?
    size = max(map(len, lst)) + pad
    if cols is None:
        import math
        cols = math.trunc((width - indent) / size)
    s = " " * indent
    for i in range(len(lst)):
        if i != 0 and i % cols == 0:
            s += "\n" + " " * indent
        s += "{0:{1}{2}}".format(lst[i], align, size)
    s += "\n"
    return s


def help_usage(s="", compact=False):
    """
    Append and return a brief usage message for the Sage documentation builder.

    If 'compact' is False, the function adds a final newline character.
    """
    s += "sage --docbuild [OPTIONS] DOCUMENT (FORMAT | COMMAND)"
    if not compact:
        s += "\n"
    return s


def help_description(s="", compact=False):
    """
    Append and return a brief description of the Sage documentation builder.

    If 'compact' is ``False``, the function adds a final newline character.
    """
    s += "Build or return information about Sage documentation. "
    s += "A DOCUMENT and either a FORMAT or a COMMAND are required."
    if not compact:
        s += "\n"
    return s


def help_examples(s=""):
    """
    Append and return some usage examples for the Sage documentation builder.
    """
    s += "Examples:\n"
    s += "    sage --docbuild -C all\n"
    s += "    sage --docbuild constructions pdf\n"
    s += "    sage --docbuild reference html -jv3\n"
    s += "    sage --docbuild reference print_unincluded_modules\n"
    s += "    sage --docbuild developer html --sphinx-opts='-q,-aE' --verbose 2"
    return s


def help_documents():
    """
    Append and return a tabular list of documents, including a
    shortcut 'all' for all documents, available to the Sage
    documentation builder.
    """
    docs = get_documents()
    s = "DOCUMENTs:\n"
    s += format_columns(docs)
    s += "\n"
    if 'reference' in docs:
        s += "Other valid document names take the form 'reference/DIR', where\n"
        s += "DIR is a subdirectory of SAGE_DOC_SRC/en/reference/.\n"
        s += "This builds just the specified part of the reference manual.\n"
    s += "DOCUMENT may also have the form 'file=/path/to/FILE', which builds\n"
    s += "the documentation for the specified file.\n"
    return s


def get_formats():
    """
    Return a list of output formats the Sage documentation builder
    will accept on the command-line.
    """
    tut_b = DocBuilder('en/tutorial')
    formats = tut_b._output_formats()
    formats.remove('html')
    return ['html', 'pdf'] + formats


def help_formats():
    """
    Append and return a tabular list of output formats available to
    the Sage documentation builder.
    """
    return "FORMATs:\n" + format_columns(get_formats())


def help_commands(name='all'):
    """
    Append and return a tabular list of commands, if any, the Sage
    documentation builder can run on the indicated document.  The
    default is to list all commands for all documents.
    """
    # To do: Generate the lists dynamically, using class attributes,
    # as with the Builders above.
    s = ""
    command_dict = {'reference': [
        'print_included_modules', 'print_modified_modules        (*)',
        'print_unincluded_modules', 'print_new_and_updated_modules (*)']}
    for doc in command_dict:
        if name == 'all' or doc == name:
            s += "COMMANDs for the DOCUMENT '" + doc + "':\n"
            s += format_columns(command_dict[doc])
            s += "(*) Since the last build.\n"
    return s


class help_message_long(argparse.Action):
    """
    Print an extended help message for the Sage documentation builder
    and exits.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        help_funcs = [help_usage, help_description, help_documents,
                      help_formats, help_commands]
        for f in help_funcs:
            print(f())
        parser.print_help()
        print(help_examples())
        sys.exit(0)


class help_message_short(argparse.Action):
    """
    Print a help message for the Sage documentation builder.

    The message includes command-line usage and a list of options.
    The message is printed only on the first call.  If error is True
    during this call, the message is printed only if the user hasn't
    requested a list (e.g., documents, formats, commands).
    """
    def __call__(self, parser, namespace, values, option_string=None):
        if not hasattr(namespace, 'printed_help'):
            parser.print_help()
            setattr(namespace, 'printed_help', 1)
        sys.exit(0)


class help_wrapper(argparse.Action):
    """
    A helper wrapper for command-line options to the Sage
    documentation builder that print lists, such as document names,
    formats, and document-specific commands.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        if option_string in ['-D', '--documents']:
            print(help_documents(), end="")
        if option_string in ['-F', '--formats']:
            print(help_formats(), end="")
        if self.dest == 'commands':
            print(help_commands(values), end="")
        if self.dest == 'all_documents':
            if values == 'reference':
                b = ReferenceBuilder('reference')
                refdir = os.path.join(os.environ['SAGE_DOC_SRC'], 'en', b.name)
                s = b.get_all_documents(refdir)
                # Put the bibliography first, because it needs to be built first:
                s.remove('reference/references')
                s.insert(0, 'reference/references')
            elif values == 'all':
                s = get_documents()
                # Put the reference manual first, because it needs to be built first:
                s.remove('reference')
                s.insert(0, 'reference')
            for d in s:
                print(d)
        setattr(namespace, 'printed_list', 1)
        sys.exit(0)


def setup_parser():
    """
    Set up and return a command-line ArgumentParser instance for the
    Sage documentation builder.
    """
    # Documentation: https://docs.python.org/library/argparse.html
    parser = argparse.ArgumentParser(usage=help_usage(compact=True),
                                     description=help_description(compact=True),
                                     add_help=False)
    # Standard options. Note: We use explicit option.dest names
    # to avoid ambiguity.
    standard = parser.add_argument_group("Standard")
    standard.add_argument("-h", "--help", nargs=0, action=help_message_short,
                          help="show a help message and exit")
    standard.add_argument("-H", "--help-all", nargs=0, action=help_message_long,
                          help="show an extended help message and exit")
    standard.add_argument("-D", "--documents", nargs=0, action=help_wrapper,
                          help="list all available DOCUMENTs")
    standard.add_argument("-F", "--formats", nargs=0, action=help_wrapper,
                          help="list all output FORMATs")
    standard.add_argument("-C", "--commands", dest="commands",
                          type=str, metavar="DOC", action=help_wrapper,
                          help="list all COMMANDs for DOCUMENT DOC; use 'all' to list all")
    standard.add_argument("-i", "--inherited", dest="inherited",
                          action="store_true",
                          help="include inherited members in reference manual; may be slow, may fail for PDF output")
    standard.add_argument("-u", "--underscore", dest="underscore",
                          action="store_true",
                          help="include variables prefixed with '_' in reference manual; may be slow, may fail for PDF output")
    standard.add_argument("-j", "--mathjax", "--jsmath", dest="mathjax",
                          action="store_true",
                          help="ignored for backwards compatibility")
    standard.add_argument("--no-plot", dest="no_plot",
                          action="store_true",
                          help="do not include graphics auto-generated using the '.. plot' markup")
    standard.add_argument("--no-preparsed-examples", dest="no_preparsed_examples",
                          action="store_true",
                          help="do not show preparsed versions of EXAMPLES blocks")
    standard.add_argument("--include-tests-blocks", dest="skip_tests", default=True,
                          action="store_false",
                          help="include TESTS blocks in the reference manual")
    standard.add_argument("--no-pdf-links", dest="no_pdf_links",
                          action="store_true",
                          help="do not include PDF links in DOCUMENT 'website'; FORMATs: html, json, pickle, web")
    standard.add_argument("--live-doc", dest="live_doc",
                          action="store_true",
                          help="make Sage code blocks live for html FORMAT")
    standard.add_argument("--warn-links", dest="warn_links",
                          action="store_true",
                          help="issue a warning whenever a link is not properly resolved; equivalent to '--sphinx-opts -n' (sphinx option: nitpicky)")
    standard.add_argument("--check-nested", dest="check_nested",
                          action="store_true",
                          help="check picklability of nested classes in DOCUMENT 'reference'")
    standard.add_argument("--no-prune-empty-dirs", dest="no_prune_empty_dirs",
                          action="store_true",
                          help="do not prune empty directories in the documentation source")
    standard.add_argument("--use-cdns", dest="use_cdns", default=False,
                          action="store_true",
                          help="assume internet connection and use CDNs; in particular, use MathJax CDN")
    standard.add_argument("-N", "--no-colors", dest="color",
                          action="store_false",
                          help="do not color output; does not affect children")
    standard.add_argument("-q", "--quiet", dest="verbose",
                          action="store_const", const=0,
                          help="work quietly; same as --verbose=0")
    standard.add_argument("-v", "--verbose", dest="verbose",
                          type=int, default=1, metavar="LEVEL",
                          action="store",
                          help="report progress at LEVEL=0 (quiet), 1 (normal), 2 (info), or 3 (debug); does not affect children")
    standard.add_argument("-o", "--output", dest="output_dir", default=None,
                          metavar="DIR", action="store",
                          help="if DOCUMENT is a single file ('file=...'), write output to this directory")

    # Advanced options.
    advanced = parser.add_argument_group("Advanced",
                                         "Use these options with care.")
    advanced.add_argument("-S", "--sphinx-opts", dest="sphinx_opts",
                          type=str, metavar="OPTS",
                          action="store",
                          help="pass comma-separated OPTS to sphinx-build; must precede OPTS with '=', as in '-S=-q,-aE' or '-S=\"-q,-aE\"'")
    advanced.add_argument("-U", "--update-mtimes", dest="update_mtimes",
                          action="store_true",
                          help="before building reference manual, update modification times for auto-generated reST files")
    advanced.add_argument("-k", "--keep-going", dest="keep_going",
                          action="store_true",
                          help="Do not abort on errors but continue as much as possible after an error")
    advanced.add_argument("--all-documents", dest="all_documents",
                          type=str, metavar="ARG",
                          choices=['all', 'reference'],
                          action=help_wrapper,
                          help="if ARG is 'reference', list all subdocuments"
                          " of en/reference. If ARG is 'all', list all main"
                          " documents")
    parser.add_argument("document", nargs='?', type=str, metavar="DOCUMENT",
                        help="name of the document to build. It can be either one of the documents listed by -D or 'file=/path/to/FILE' to build documentation for this specific file.")
    parser.add_argument("format", nargs='?', type=str,
                        metavar="FORMAT or COMMAND", help='document output format (or command)')
    return parser


def setup_logger(verbose=1, color=True):
    r"""
    Set up a Python Logger instance for the Sage documentation builder.

    The optional argument sets logger's level and message format.

    EXAMPLES::

        sage: from sage_docbuild.__main__ import setup_logger, logger
        sage: setup_logger()
        sage: type(logger)
        <class 'logging.Logger'>
    """
    # Set up colors. Adapted from sphinx.cmdline.
    import sphinx.util.console as c
    if not color or not sys.stdout.isatty() or not c.color_terminal():
        c.nocolor()

    # Available colors: black, darkgray, (dark)red, dark(green),
    # brown, yellow, (dark)blue, purple, fuchsia, turquoise, teal,
    # lightgray, white.  Available styles: reset, bold, faint,
    # standout, underline, blink.

    # Set up log record formats.
    format_std = "%(message)s"
    formatter = logging.Formatter(format_std)

    # format_debug = "%(module)s #%(lineno)s %(funcName)s() %(message)s"
    fields = ['%(module)s', '#%(lineno)s', '%(funcName)s()', '%(message)s']
    colors = ['darkblue', 'darkred', 'brown', 'reset']
    styles = ['reset', 'reset', 'reset', 'reset']
    format_debug = ""
    for i in range(len(fields)):
        format_debug += c.colorize(styles[i], c.colorize(colors[i], fields[i]))
        if i != len(fields):
            format_debug += " "

    # Note: There's also Handler.setLevel().  The argument is the
    # lowest severity message that the respective logger or handler
    # will pass on.  The default levels are DEBUG, INFO, WARNING,
    # ERROR, and CRITICAL.  We use "WARNING" for normal verbosity and
    # "ERROR" for quiet operation.  It's possible to define custom
    # levels.  See the documentation for details.
    if verbose == 0:
        logger.setLevel(logging.ERROR)
    if verbose == 1:
        logger.setLevel(logging.WARNING)
    if verbose == 2:
        logger.setLevel(logging.INFO)
    if verbose == 3:
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(format_debug)

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)


class IntersphinxCache:
    """
    Replace sphinx.ext.intersphinx.fetch_inventory by an in-memory
    cached version.
    """
    def __init__(self):
        self.inventories = {}
        self.real_fetch_inventory = sphinx.ext.intersphinx.fetch_inventory
        sphinx.ext.intersphinx.fetch_inventory = self.fetch_inventory

    def fetch_inventory(self, app, uri, inv):
        """
        Return the result of ``sphinx.ext.intersphinx.fetch_inventory()``
        from a cache if possible. Otherwise, call
        ``sphinx.ext.intersphinx.fetch_inventory()`` and cache the result.
        """
        t = (uri, inv)
        try:
            return self.inventories[t]
        except KeyError:
            i = self.real_fetch_inventory(app, uri, inv)
            self.inventories[t] = i
            return i


def main():
    # Parse the command-line.
    parser = setup_parser()
    args = parser.parse_args()
    DocBuilder._options = args

    # Get the name and type (target format) of the document we are
    # trying to build.
    name, typ = args.document, args.format
    if not name or not typ:
        parser.print_help()
        sys.exit(1)
    elif name == 'all':
        sys.exit(os.system(f'cd {shlex.quote(SAGE_DOC_SRC)} '
                           f'&& ${{MAKE:-make}} -j${{SAGE_NUM_THREADS_PARALLEL:-1}} doc-{typ}'))

    # Set up module-wide logging.
    setup_logger(args.verbose, args.color)

    def excepthook(*exc_info):
        logger.error('Error building the documentation.', exc_info=exc_info)
        if build_options.INCREMENTAL_BUILD:
            logger.error('''
    Note: incremental documentation builds sometimes cause spurious
    error messages. To be certain that these are real errors, run
    "make doc-clean doc-uninstall" first and try again.''')

    sys.excepthook = excepthook

    # Process selected options.
    if args.check_nested:
        os.environ['SAGE_CHECK_NESTED'] = 'True'

    if args.underscore:
        os.environ['SAGE_DOC_UNDERSCORE'] = "True"

    if args.sphinx_opts:
        build_options.ALLSPHINXOPTS += args.sphinx_opts.replace(',', ' ') + " "
    if args.no_pdf_links:
        build_options.WEBSITESPHINXOPTS = " -A hide_pdf_links=1 "
    if args.warn_links:
        build_options.ALLSPHINXOPTS += "-n "
    if args.no_plot:
        os.environ['SAGE_SKIP_PLOT_DIRECTIVE'] = 'yes'
    if args.no_preparsed_examples:
        os.environ['SAGE_PREPARSED_DOC'] = 'no'
    if args.live_doc:
        os.environ['SAGE_LIVE_DOC'] = 'yes'
    if args.skip_tests:
        os.environ['SAGE_SKIP_TESTS_BLOCKS'] = 'True'
    if args.use_cdns:
        os.environ['SAGE_USE_CDNS'] = 'yes'

    build_options.ABORT_ON_ERROR = not args.keep_going

    # Set up Intersphinx cache
    _ = IntersphinxCache()

    builder = get_builder(name)

    if not args.no_prune_empty_dirs:
        # Delete empty directories. This is needed in particular for empty
        # directories due to "git checkout" which never deletes empty
        # directories it leaves behind. See Issue #20010.
        # Issue #31948: This is not parallelization-safe; use the option
        # --no-prune-empty-dirs to turn it off
        for dirpath, dirnames, filenames in os.walk(builder.dir, topdown=False):
            if not dirnames + filenames:
                logger.warning('Deleting empty directory {0}'.format(dirpath))
                os.rmdir(dirpath)

    build = getattr(builder, typ)
    build()


if __name__ == '__main__':
    sys.exit(main())
