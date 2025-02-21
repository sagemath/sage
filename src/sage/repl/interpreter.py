# sage_setup: distribution = sagemath-repl
r"""
Sage's IPython Modifications

This module contains all of Sage's customizations to the IPython
interpreter.  These changes consist of the following major components:

  - :class:`SageTerminalApp`
  - :class:`SageInteractiveShell`
  - :class:`SageTerminalInteractiveShell`
  - :func:`interface_shell_embed`

SageTerminalApp
---------------

This is the main application object.  It is used by the
``$SAGE_LOCAL/bin/sage-ipython`` script to start the Sage
command-line.  It's primary purpose is to

  - Initialize the :class:`SageTerminalInteractiveShell`.

  - Provide default configuration options for the shell, and its
    subcomponents.  These work with (and can be overridden by)
    IPython's configuration system.

  - Load the Sage ipython extension (which does things like preparsing,
    add magics, etc.).

  - Provide a custom :class:`SageCrashHandler` to give the user
    instructions on how to report the crash to the Sage support
    mailing list.

SageInteractiveShell
--------------------

The :class:`SageInteractiveShell` object is the object responsible for
accepting input from the user and evaluating it.  From the command-line,
this object can be retrieved by running::

    sage: shell = get_ipython()   # not tested

Any input is preprocessed and evaluated inside the ``shell.run_cell``
method. If the command line processing does not do what you want it to
do, you can step through it in the debugger::

    sage: %debug shell.run_cell('?')        # not tested

The :class:`SageInteractiveShell` provides the following
customizations:

  - Modify the libraries before calling system commands. See
    :meth:`~SageInteractiveShell.system_raw`.

SageTerminalInteractiveShell
----------------------------

The :class:`SageTerminalInteractiveShell` is a close relative of
:class:`SageInteractiveShell` that is specialized for running in a
terminal. In particular, running commands like ``!ls`` will directly
write to stdout. Technically, the ``system`` attribute will point to
``system_raw`` instead of ``system_piped``.

Interface Shell
---------------

The function :func:`interface_shell_embed` takes a
:class:`~sage.interfaces.interface.Interface` object and returns an
embeddable IPython shell which can be used to directly interact with
that shell.  The bulk of this functionality is provided through
:class:`InterfaceShellTransformer`.

TESTS:

Check that Cython source code appears in tracebacks::

    sage: from sage.repl.interpreter import get_test_shell
    sage: shell = get_test_shell()
    sage: shell.run_cell('1/0') # known bug (meson doesn't include the Cython source code)
    ...
    ZeroDivisionError...Traceback (most recent call last)
    ...
    ----> 1 Integer(1)/Integer(0)
    ...
    ...integer.pyx... in sage.rings.integer.Integer...div...
    ...
    -> ...                  raise ZeroDivisionError("rational division by zero")
       ....:            x = <Rational> Rational.__new__(Rational)
       ....:            mpq_div_zz(x.value, ....value, (<Integer>right).value)
    <BLANKLINE>
    ZeroDivisionError: rational division by zero
    sage: shell.quit()

Test prompt transformer::

    sage: from sage.repl.interpreter import SagePromptTransformer
    sage: spt = SagePromptTransformer
    sage: spt(["sage: 2 + 2"])
    ['2 + 2']
    sage: spt([''])
    ['']
    sage: spt(["....: 2+2"])
    ['2+2']

This should strip multiple prompts: see :issue:`16297`::

    sage: spt(["sage:   sage: 2+2"])
    ['2+2']
    sage: spt(["   sage: ....: 2+2"])
    ['2+2']

The prompt contains a trailing space. Extra spaces between the
last prompt and the remainder should not be stripped::

    sage: spt(["   sage: ....:    2+2"])
    ['   2+2']

We test that the input transformer is enabled on the Sage command
line::

    sage: from sage.repl.interpreter import get_test_shell
    sage: shell = get_test_shell()
    sage: shell.run_cell('sage: a = 123')              # single line
    sage: shell.run_cell('sage: a = [\n... 123]')      # old-style multi-line
    sage: shell.run_cell('sage: a = [\n....: 123]')    # new-style multi-line

We test that :issue:`16196` is resolved::

    sage: shell.run_cell('    sage: 1+1')
    2
    sage: shell.quit()
"""

# ****************************************************************************
#       Copyright (C) 2004-2012 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import re
from traitlets import Bool, Type

from sage.repl.preparse import preparse, containing_block
from sage.repl.prompts import InterfacePrompts
from sage.repl.configuration import sage_ipython_config, SAGE_EXTENSION

from IPython.core.interactiveshell import InteractiveShell
from IPython.terminal.interactiveshell import TerminalInteractiveShell
from IPython.core.inputtransformer2 import PromptStripper
from IPython.core.prefilter import PrefilterTransformer
from IPython.terminal.embed import InteractiveShellEmbed
from IPython.terminal.ipapp import TerminalIPythonApp, IPAppCrashHandler
from IPython.core.crashhandler import CrashHandler

from ctypes import pythonapi, c_int, c_void_p
# The following functions are part of the stable ABI since python 3.2
# See: https://docs.python.org/3/c-api/sys.html#c.PyOS_getsig

# PyOS_sighandler_t PyOS_getsig(int i)
pythonapi.PyOS_getsig.restype = c_void_p
pythonapi.PyOS_getsig.argtypes = c_int,

# PyOS_sighandler_t PyOS_setsig(int i, PyOS_sighandler_t h)
pythonapi.PyOS_setsig.restype = c_void_p
pythonapi.PyOS_setsig.argtypes = c_int, c_void_p,


# TODO: This global variable _do_preparse should be associated with an
# IPython InteractiveShell as opposed to a global variable in this
# module.
_do_preparse = True


def preparser(on=True):
    """
    Turn on or off the Sage preparser.

    - ``on`` -- boolean; whether to turn on preparsing

    EXAMPLES::

        sage: 2/3
        2/3
        sage: preparser(False)
        sage: 2/3  # not tested since doctests are always preparsed
        0
        sage: preparser(True)
        sage: 2^3
        8
    """
    global _do_preparse
    _do_preparse = on is True


##############################
# Sage[Terminal]InteractiveShell
##############################
class SageShellOverride:
    """
    Mixin to override methods in IPython's [Terminal]InteractiveShell
    classes.
    """

    def show_usage(self):
        """
        Print the basic Sage usage.

        This method ends up being called when you enter ``?`` and
        nothing else on the command line.

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('?')
            Welcome to Sage ...
            sage: shell.quit()
        """
        from sage.misc.sagedoc import help
        help()

    def system_raw(self, cmd):
        """
        Run a system command.

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.system_raw('false')
            sage: shell.user_ns['_exit_code'] > 0
            True
            sage: shell.system_raw('true')
            sage: shell.user_ns['_exit_code']
            0
            sage: shell.system_raw('R --version')   # optional - r
            R version ...
            sage: shell.user_ns['_exit_code']       # optional - r
            0
            sage: shell.quit()
        """
        return super().system_raw(cmd)


class SageNotebookInteractiveShell(SageShellOverride, InteractiveShell):
    """
    IPython Shell for the Sage IPython Notebook.

    The doctests are not tested since they would change the current
    rich output backend away from the doctest rich output backend.

    EXAMPLES::

        sage: from sage.repl.interpreter import SageNotebookInteractiveShell
        sage: SageNotebookInteractiveShell()   # not tested
        <sage.repl.interpreter.SageNotebookInteractiveShell object at 0x...>
    """

    def init_display_formatter(self):
        """
        Switch to the Sage IPython notebook rich output backend.

        EXAMPLES::

            sage: from sage.repl.interpreter import SageNotebookInteractiveShell
            sage: SageNotebookInteractiveShell().init_display_formatter()   # not tested
        """
        from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
        backend = BackendIPythonNotebook()
        backend.get_display_manager().switch_backend(backend, shell=self)


class SageTerminalInteractiveShell(SageShellOverride, TerminalInteractiveShell):
    """
    IPython Shell for the Sage IPython Commandline Interface.

    The doctests are not tested since they would change the current
    rich output backend away from the doctest rich output backend.

    EXAMPLES::

        sage: from sage.repl.interpreter import SageTerminalInteractiveShell
        sage: SageTerminalInteractiveShell()   # not tested
        <sage.repl.interpreter.SageNotebookInteractiveShell object at 0x...>
    """

    def init_display_formatter(self):
        """
        Switch to the Sage IPython commandline rich output backend.

        EXAMPLES::

            sage: from sage.repl.interpreter import SageTerminalInteractiveShell
            sage: SageTerminalInteractiveShell().init_display_formatter()   # not tested
        """
        from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
        backend = BackendIPythonCommandline()
        backend.get_display_manager().switch_backend(backend, shell=self)

    def prompt_for_code(self):
        # save sigint handlers (python and os level)
        # https://github.com/prompt-toolkit/python-prompt-toolkit/issues/1576
        # https://github.com/sagemath/sage/issues/33428
        # https://github.com/sagemath/sage/pull/35251
        import signal
        sigint = signal.getsignal(signal.SIGINT)
        sigint_os = pythonapi.PyOS_getsig(signal.SIGINT)
        text = TerminalInteractiveShell.prompt_for_code(self)
        signal.signal(signal.SIGINT, sigint)
        pythonapi.PyOS_setsig(signal.SIGINT, sigint_os)
        return text


class SageTestShell(SageShellOverride, TerminalInteractiveShell):
    """
    Test Shell.

    Care must be taken in these doctests to quit the test shell in
    order to switch back the rich output display backend to the
    doctest backend.

    EXAMPLES::

        sage: from sage.repl.interpreter import get_test_shell
        sage: shell = get_test_shell();  shell
        <sage.repl.interpreter.SageTestShell object at 0x...>
        sage: shell.quit()
    """

    def init_display_formatter(self):
        """
        Switch to the Sage IPython commandline rich output backend.

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell();  shell
            <sage.repl.interpreter.SageTestShell object at 0x...>
            sage: shell.quit()
            sage: shell.init_display_formatter()
            sage: shell.quit()
        """
        from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
        self._ipython_backend = backend = BackendIPythonCommandline()
        self._display_manager = backend.get_display_manager()
        self._doctest_backend = self._display_manager.switch_backend(backend, shell=self)

    def quit(self):
        """
        Quit the test shell.

        To make the test shell as realistic as possible, we switch to
        the
        :class:`~sage.repl.rich_output.backend_ipython.BackendIPythonCommandline`
        display backend. This method restores the previous display
        backend, which is the
        :class:`~sage.repl.rich_output.backend_doctest.BackendDoctest`
        during doctests.

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: from sage.repl.rich_output import get_display_manager
            sage: get_display_manager()
            The Sage display manager using the doctest backend

            sage: shell = get_test_shell()
            sage: get_display_manager()
            The Sage display manager using the IPython command line backend

            sage: shell.quit()
            sage: get_display_manager()
            The Sage display manager using the doctest backend
        """
        self._display_manager.switch_backend(self._doctest_backend)

    def _restart(self):
        """
        Restart the test shell (after :meth:`quit`).

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.quit()
            sage: shell._restart()
            sage: shell.quit()
            sage: from sage.repl.rich_output import get_display_manager
            sage: get_display_manager()
            The Sage display manager using the doctest backend
        """
        self._display_manager.switch_backend(self._ipython_backend, shell=self)

    def run_cell(self, *args, **kwds):
        """
        Run IPython cell.

        Starting with IPython-3.0, this returns an success/failure
        information. Since it is more convenient for doctests, we
        ignore it.

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: rc = shell.run_cell('2^50')
            1125899906842624
            sage: rc is None
            True
            sage: shell.quit()
        """
        super().run_cell(*args, **kwds)


###################################################################
# Transformers used in the SageInputSplitter
###################################################################


def SagePreparseTransformer(lines):
    r"""
    EXAMPLES::

        sage: from sage.repl.interpreter import SagePreparseTransformer
        sage: spt = SagePreparseTransformer
        sage: spt(['1+1r+2.3^2.3r\n'])
        ["Integer(1)+1+RealNumber('2.3')**2.3\n"]
        sage: preparser(False)
        sage: spt(['2.3^2\n'])
        ['2.3^2\n']

    TESTS:

    Check that syntax errors in the preparser do not crash IPython,
    see :issue:`14961`. ::

        sage: preparser(True)
        sage: bad_syntax = "R.<t> = QQ{]"
        sage: preparse(bad_syntax)
        Traceback (most recent call last):
        ...
        SyntaxError: mismatched ']'
        sage: from sage.repl.interpreter import get_test_shell
        sage: shell = get_test_shell()
        sage: shell.run_cell(bad_syntax)
          File...<string>...
        SyntaxError: mismatched ']'
        <BLANKLINE>
        sage: shell.quit()

    Make sure the quote state is carried over across subsequent lines in order
    to avoid interfering with multi-line strings, see :issue:`30417`. ::

        sage: SagePreparseTransformer(["'''\n", 'abc-1-2\n', "'''\n"])
        ["'''\n", 'abc-1-2\n', "'''\n"]
        sage: # instead of ["'''\n", 'abc-Integer(1)-Integer(2)\n', "'''\n"]

    .. NOTE::

        IPython may call this function more than once for the same input lines.
        So when debugging the preparser, print outs may be duplicated. If using
        IPython >= 7.17, try:
        ``sage.repl.interpreter.SagePreparseTransformer.has_side_effects = True``
    """
    if _do_preparse:
        # IPython ensures the input lines end with a newline, and it expects
        # the same of the output lines.
        lines = preparse(''.join(lines)).splitlines(keepends=True)
    return lines


SagePromptTransformer = PromptStripper(prompt_re=re.compile(r'^(\s*(:?sage: |\.\.\.\.: ))+'))


###################
# Interface shell #
###################

class logstr(str):
    """
    For use by :meth`~InterfaceShellTransformer.transform`.
    This provides a ``_latex_`` method which is just the string
    wrapped in a ``\\verb`` environment.
    """
    def __repr__(self):
        """
        EXAMPLES::
            sage: from sage.repl.interpreter import logstr
            sage: logstr("ABC")
            ABC
        """
        return self

    def _latex_(self):
        r"""
        EXAMPLES::
            sage: from sage.repl.interpreter import logstr
            sage: logstr("A")
            A
            sage: latex(logstr("A"))
            \verb#A#
            sage: latex(logstr("A#B"))
            \verb@A#B@
        """
        # return "\\begin{verbatim}%s\\end{verbatim}"%self
        if '#' not in self:
            delim = '#'
        elif '@' not in self:
            delim = '@'
        elif '~' not in self:
            delim = '~'
        return r"""\verb%s%s%s""" % (delim, self.replace('\n\n', '\n').replace('\n', '; '), delim)


class InterfaceShellTransformer(PrefilterTransformer):
    priority = 50

    def __init__(self, *args, **kwds):
        """
        Initialize this class.  All of the arguments get passed to
        :meth:`PrefilterTransformer.__init__`.

        .. attribute:: temporary_objects

           a list of hold onto interface objects and keep them from being
           garbage collected

        .. SEEALSO:: :func:`interface_shell_embed`

        TESTS::

            sage: # needs sage.symbolic
            sage: from sage.repl.interpreter import interface_shell_embed
            sage: shell = interface_shell_embed(maxima)
            sage: ift = shell.prefilter_manager.transformers[0]
            sage: ift.temporary_objects
            set()
            sage: ift._sage_import_re.findall('sage(a) + maxima(b)')
            ['sage(', 'maxima(']
        """
        super().__init__(*args, **kwds)
        self.temporary_objects = set()
        self._sage_import_re = re.compile(r'(?:sage|%s)\('
                                          % self.shell.interface.name())

    def preparse_imports_from_sage(self, line):
        """
        Finds occurrences of strings such as ``sage(object)`` in
        *line*, converts ``object`` to :attr:`shell.interface`,
        and replaces those strings with their identifier in the new
        system.  This also works with strings such as
        ``maxima(object)`` if :attr:`shell.interface` is
        ``maxima``.

        - ``line`` -- string; the line to transform

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.repl.interpreter import interface_shell_embed, InterfaceShellTransformer
            sage: shell = interface_shell_embed(maxima)
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config,
            ....:     prefilter_manager=shell.prefilter_manager)
            sage: ift.shell.ex('a = 3')
            sage: ift.preparse_imports_from_sage('2 + sage(a)')
            '2 + sage0 '
            sage: maxima.eval('sage0')
            '3'
            sage: ift.preparse_imports_from_sage('2 + maxima(a)') # maxima calls set_seed on startup which is why 'sage0' will becomes 'sage4' and not just 'sage1'
            '2 +  sage4 '
            sage: ift.preparse_imports_from_sage('2 + gap(a)')
            '2 + gap(a)'

        Since :issue:`28439`, this also works with more complicated expressions
        containing nested parentheses::

            sage: # needs sage.libs.gap sage.symbolic
            sage: shell = interface_shell_embed(gap)
            sage: shell.user_ns = locals()
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config,
            ....:     prefilter_manager=shell.prefilter_manager)
            sage: line = '2 + sage((1+2)*gap(-(5-3)^2).sage()) - gap(1+(2-1))'
            sage: line = ift.preparse_imports_from_sage(line)
            sage: gap.eval(line)
            '-12'
        """
        new_line = []
        pos = 0
        while True:
            m = self._sage_import_re.search(line, pos)
            if not m:
                new_line.append(line[pos:])
                break
            expr_start, expr_end = containing_block(line, m.end() - 1,
                                                    delimiters=['()'])
            expr = preparse(line[expr_start + 1:expr_end - 1])
            result = self.shell.interface(eval(expr, self.shell.user_ns))
            self.temporary_objects.add(result)
            new_line += [line[pos:m.start()], result.name()]
            pos = expr_end
        return ' '.join(new_line)

    def transform(self, line, continue_prompt):
        r'''
        Evaluates *line* in :attr:`shell.interface` and returns a
        string representing the result of that evaluation.

        - ``line`` -- string; the line to be transformed *and evaluated*
        - ``continue_prompt`` -- boolean; whether this line is a continuation in a
          sequence of multiline input

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.repl.interpreter import interface_shell_embed, InterfaceShellTransformer
            sage: shell = interface_shell_embed(maxima)
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config,
            ....:     prefilter_manager=shell.prefilter_manager)
            sage: ift.transform('2+2', False)   # note: output contains triple quotation marks
            'sage.repl.interpreter.logstr(r"""4""")'
            sage: ift.shell.ex('a = 4')
            sage: ift.transform(r'sage(a)+4', False)
            'sage.repl.interpreter.logstr(r"""8""")'
            sage: ift.temporary_objects
            set()
            sage: shell = interface_shell_embed(gap)                                    # needs sage.libs.gap
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config,     # needs sage.libs.gap
            ....:     prefilter_manager=shell.prefilter_manager)
            sage: ift.transform('2+2', False)
            'sage.repl.interpreter.logstr(r"""4""")'

        TESTS:

        Check that whitespace is not stripped and that special characters are
        escaped (:issue:`28439`)::

            sage: shell = interface_shell_embed(gap)                                    # needs sage.libs.gap sage.symbolic
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config,     # needs sage.libs.gap sage.symbolic
            ....:     prefilter_manager=shell.prefilter_manager)
            sage: ift.transform(r'Print("  -\n\\\\-  ");', False)                       # needs sage.symbolic
            'sage.repl.interpreter.logstr(r"""  -\n\\\\-""")'

            sage: # optional - macaulay2
            sage: shell = interface_shell_embed(macaulay2)
            sage: ift = InterfaceShellTransformer(shell=shell, config=shell.config,
            ....:     prefilter_manager=shell.prefilter_manager)
            sage: ift.transform('net(ZZ^2)', False)
            'sage.repl.interpreter.logstr(r"""  2\nZZ""")'
        '''
        line = self.preparse_imports_from_sage(line)

        try:
            t = self.shell.interface.eval(line)
        finally:
            # Once we've evaluated the lines, we can clear the
            # temporary objects
            self.temporary_objects = set()
        # We do not strip whitespace from t here as the individual interface is
        # responsible for that
        return 'sage.repl.interpreter.logstr(r"""%s""")' % t


def interface_shell_embed(interface):
    """
    Return an IPython shell which uses a Sage interface on the
    backend to perform the evaluations.  It uses
    :class:`InterfaceShellTransformer` to transform the input into the
    appropriate ``interface.eval(...)`` input.

    INPUT:

    - ``interface`` -- a Sage ``PExpect`` interface instance

    EXAMPLES::

        sage: from sage.repl.interpreter import interface_shell_embed
        sage: shell = interface_shell_embed(gap)                                        # needs sage.libs.gap
        sage: shell.run_cell('List( [1..10], IsPrime )')                                # needs sage.libs.gap
        [ false, true, true, false, true, false, true, false, false, false ]
        <ExecutionResult object at ..., execution_count=None error_before_exec=None error_in_exec=None ...result=[ false, true, true, false, true, false, true, false, false, false ]>
    """
    cfg = sage_ipython_config.copy()
    ipshell = InteractiveShellEmbed(config=cfg,
                                    banner1='\n  --> Switching to %s <--\n\n' % interface,
                                    exit_msg='\n  --> Exiting back to Sage <--\n')
    ipshell.interface = interface
    ipshell.prompts = InterfacePrompts(interface.name())

    while ipshell.prefilter_manager.transformers:
        ipshell.prefilter_manager.transformers.pop()
    while ipshell.prefilter_manager.checkers:
        ipshell.prefilter_manager.checkers.pop()
    ipshell.ex('import sage.misc.all')

    InterfaceShellTransformer(shell=ipshell,
                              prefilter_manager=ipshell.prefilter_manager,
                              config=cfg)
    return ipshell


def get_test_shell():
    """
    Return a IPython shell that can be used in testing the functions
    in this module.

    OUTPUT: an IPython shell

    EXAMPLES::

        sage: from sage.repl.interpreter import get_test_shell
        sage: shell = get_test_shell(); shell
        <sage.repl.interpreter.SageTestShell object at 0x...>
        sage: shell.parent.shell_class
        <class 'sage.repl.interpreter.SageTestShell'>
        sage: shell.parent.test_shell
        True
        sage: shell.quit()

    TESTS:

    Check that :issue:`14070` has been resolved::

        sage: from sage.tests.cmdline import test_executable
        sage: cmd = 'from sage.repl.interpreter import get_test_shell; shell = get_test_shell()'
        sage: (out, err, ret) = test_executable(["sage", "-c", cmd])
        sage: out + err
        ''
    """
    config = sage_ipython_config.default()
    config.TerminalIPythonApp.test_shell = True
    config.TerminalIPythonApp.shell_class = SageTestShell
    app = SageTerminalApp.instance(config=config)
    if app.shell is None:
        app.initialize(argv=[])
    else:
        try:
            app.shell._restart()
        except AttributeError:
            pass
    # No quit noise
    app.shell.verbose_quit = False
    return app.shell


#######################
# IPython TerminalApp #
#######################

class SageCrashHandler(IPAppCrashHandler):
    def __init__(self, app):
        """
        A custom :class:`CrashHandler` which gives the user
        instructions on how to post the problem to sage-support.

        EXAMPLES::

            sage: from sage.repl.interpreter import SageTerminalApp, SageCrashHandler
            sage: app = SageTerminalApp.instance()
            sage: sch = SageCrashHandler(app); sch
            <sage.repl.interpreter.SageCrashHandler object at 0x...>
            sage: sorted(sch.info.items())
            [('app_name', 'Sage'),
             ('bug_tracker', 'https://github.com/sagemath/sage/issues'),
             ('contact_email', 'sage-support@googlegroups.com'),
             ('contact_name', 'sage-support'),
             ('crash_report_fname', 'Crash_report_Sage.txt')]
        """
        contact_name = 'sage-support'
        contact_email = 'sage-support@googlegroups.com'
        bug_tracker = 'https://github.com/sagemath/sage/issues'
        CrashHandler.__init__(self,
                              app, contact_name, contact_email,
                              bug_tracker, show_crash_traceback=True)
        self.crash_report_fname = 'Sage_crash_report.txt'


class SageTerminalApp(TerminalIPythonApp):
    name = 'Sage'
    crash_handler_class = SageCrashHandler

    test_shell = Bool(False, help='Whether the shell is a test shell')
    test_shell.tag(config=True)
    shell_class = Type(InteractiveShell, help='Type of the shell')
    shell_class.tag(config=True)

    def load_config_file(self, *args, **kwds):
        r"""
        Merges a config file with the default sage config.

        .. NOTE::

            This code is based on :meth:`Application.update_config`.

        TESTS:

        Test that :issue:`15972` has been fixed::

            sage: import tempfile
            sage: from sage.repl.interpreter import SageTerminalApp
            sage: from IPython.paths import get_ipython_dir
            sage: with tempfile.TemporaryDirectory() as d:
            ....:     IPYTHONDIR = get_ipython_dir()
            ....:     os.environ['IPYTHONDIR'] = d
            ....:     SageTerminalApp().load_config_file()
            ....:     os.environ['IPYTHONDIR'] = IPYTHONDIR
        """
        super().load_config_file(*args, **kwds)
        newconfig = sage_ipython_config.default()
        # merge in the config loaded from file
        newconfig.merge(self.config)
        self.config = newconfig

    def init_shell(self):
        r"""
        Initialize the :class:`SageInteractiveShell` instance.

        .. NOTE::

            This code is based on
            :meth:`TerminalIPythonApp.init_shell`.

        EXAMPLES::

            sage: from sage.repl.interpreter import SageTerminalApp
            sage: app = SageTerminalApp.instance()
            sage: app.shell
            <sage.repl.interpreter.SageTestShell object at 0x...>
        """
        # Shell initialization
        self.shell = self.shell_class.instance(
            parent=self,
            config=self.config,
            profile_dir=self.profile_dir,
            ipython_dir=self.ipython_dir)
        self.shell.configurables.append(self)
        self.shell.has_sage_extensions = SAGE_EXTENSION in self.extensions

        # Load the %lprun extension if available
        try:
            import line_profiler
            assert line_profiler  # silence pyflakes
        except ImportError:
            pass
        else:
            self.extensions.append('line_profiler')

        if self.shell.has_sage_extensions:
            self.extensions.remove(SAGE_EXTENSION)

            # load sage extension here to get a crash if
            # something is wrong with the sage library
            self.shell.extension_manager.load_extension(SAGE_EXTENSION)
