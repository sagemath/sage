# sage_setup: distribution = sagemath-repl
r"""
Sage's IPython Extension

A Sage extension which adds sage-specific features:

* magics

  - ``%crun``

  - ``%runfile``

  - ``%attach``

  - ``%display``

  - ``%mode`` (like ``%maxima``, etc.)

  - ``%%cython``

  - ``%%fortran``

* preparsing of input

* loading Sage library

* running init.sage

* changing prompt to Sage prompt

* Display hook

TESTS:

We test that preparsing is off for ``%runfile``, on for ``%time``::

    sage: import os, re
    sage: from sage.repl.interpreter import get_test_shell
    sage: from sage.misc.temporary_file import tmp_dir
    sage: shell = get_test_shell()
    sage: TMP = tmp_dir()
    sage: TMP = os.path.join(TMP, "12345", "temp")
    sage: os.makedirs(TMP)

The temporary directory should have a name of the form
``.../12345/...``, to demonstrate that file names are not
preparsed when calling ``%runfile`` ::

    sage: bool(re.search('/12345/', TMP))
    True
    sage: tmp = os.path.join(TMP, 'run_cell.py')
    sage: with open(tmp, 'w') as f:
    ....:     _ = f.write('a = 2\n')
    sage: shell.run_cell('%runfile '+tmp)
    sage: shell.run_cell('a')
    2

In contrast, input to the ``%time`` magic command is preparsed::

    sage: shell.run_cell('%time 594.factor()')
    CPU times: user ...
    Wall time: ...
    2 * 3^3 * 11
    sage: shell.quit()
"""

from IPython.core.magic import Magics, magics_class, line_magic, cell_magic

from sage.repl.load import load_wrap
from sage.env import SAGE_IMPORTALL, SAGE_STARTUP_FILE
from sage.misc.lazy_import import LazyImport
from sage.misc.misc import run_once


@magics_class
class SageMagics(Magics):

    @line_magic
    def crun(self, s):
        r"""
        Profile C function calls.

        INPUT:

        - ``s`` -- string; Sage command to profile

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('%crun sum(1/(1+n^2) for n in range(100))')   # optional - gperftools
            PROFILE: interrupts/evictions/bytes = ...
            Using local file ...
            Using local file ...
            sage: shell.quit()
        """
        import sage.misc.gperftools
        sage.misc.gperftools.crun(s, evaluator=self.shell.ex)

    @line_magic
    def runfile(self, s):
        r"""
        Execute the code contained in the file ``s``.

        This is designed to be used from the command line as
        ``%runfile /path/to/file``.

        - ``s`` -- string; the file to be loaded

        EXAMPLES::

            sage: import os
            sage: from sage.repl.interpreter import get_test_shell
            sage: from sage.misc.temporary_file import tmp_dir
            sage: shell = get_test_shell()
            sage: tmp = os.path.join(tmp_dir(), 'run_cell.py')
            sage: with open(tmp, 'w') as f:
            ....:     _ = f.write('a = 2\n')
            sage: shell.run_cell('%runfile '+tmp)
            sage: shell.run_cell('a')
            2
            sage: shell.quit()
        """
        return self.shell.ex(load_wrap(s, attach=False))

    @line_magic
    def attach(self, s):
        r"""
        Attach the code contained in the file ``s``.

        This is designed to be used from the command line as ``%attach
        /path/to/file``.

        - ``s`` -- string. The file to be attached

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: from tempfile import NamedTemporaryFile as NTF
            sage: with NTF(mode='w+t', suffix='.py', delete=False) as f:
            ....:     _ = f.write('a = 2\n')
            sage: shell.run_cell('%attach ' + f.name)
            sage: shell.run_cell('a')
            2
            sage: sleep(1)  # filesystem timestamp granularity
            sage: with open(f.name, 'w') as f: _ = f.write('a = 3\n')

        Note that the doctests are never really at the command prompt, so
        we call the input hook manually::

            sage: shell.run_cell('from sage.repl.attach import reload_attached_files_if_modified')
            sage: shell.run_cell('reload_attached_files_if_modified()')
            ### reloading attached file ... modified at ... ###

            sage: shell.run_cell('a')
            3
            sage: shell.run_cell('detach(%r)' % f.name)
            sage: shell.run_cell('attached_files()')
            []
            sage: os.remove(f.name)
            sage: shell.quit()
        """
        return self.shell.ex(load_wrap(s, attach=True))

    @line_magic
    def iload(self, args):
        """
        A magic command to interactively load a file as in MAGMA.

        - ``args`` -- string. The file to be interactively loaded

        .. NOTE::

            Currently, this cannot be completely doctested as it
            relies on :func:`raw_input`.

        EXAMPLES::

            sage: ip = get_ipython()           # not tested: works only in interactive shell
            sage: ip.magic_iload('/dev/null')  # not tested: works only in interactive shell
            Interactively loading "/dev/null"  # not tested: works only in interactive shell
        """
        content = self.shell.find_user_code(args).splitlines()

        # we create a stack so e.g. having an iload inside of an iload
        # will process the inner iload and then resume the outer iload
        orig_readline = self.shell.pre_readline

        def pre_readline():
            if self.shell.rl_next_input is None:
                self.shell.rl_next_input = content.pop(0)
                self.shell.rl_do_indent = False
            orig_readline()
            if not content:
                # restore original hook
                self.shell.readline_startup_hook(orig_readline)
                self.shell.pre_readline = orig_readline

        self.shell.readline_startup_hook(pre_readline)
        self.shell.pre_readline = pre_readline

        print('Interactively loading "%s"' % args)

    _magic_display_status = 'simple'

    @line_magic
    def display(self, args):
        r"""
        A magic command to switch between simple display and ASCII art display.

        - ``args`` -- string.  See
          :mod:`sage.repl.rich_output.preferences`
          for allowed values. If the mode is ``ascii_art``, it can
          optionally be followed by a width.

        How to use: if you want to activate the ASCII art mode::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('%display ascii_art')

        That means you do not have to use :func:`ascii_art` to get an ASCII art
        output::

            sage: shell.run_cell("i = var('i')")                                        # needs sage.symbolic
            sage: shell.run_cell('sum(i^2*x^i, i, 0, 10)')                              # needs sage.symbolic
                 10       9       8       7       6       5       4      3      2
            100*x   + 81*x  + 64*x  + 49*x  + 36*x  + 25*x  + 16*x  + 9*x  + 4*x  + x

        Then when you want to return to 'textual mode'::

            sage: shell.run_cell('%display text plain')
            sage: shell.run_cell('%display plain')        # shortcut for "text plain"
            sage: shell.run_cell('sum(i^2*x^i, i, 0, 10)')                              # needs sage.symbolic
            100*x^10 + 81*x^9 + 64*x^8 + 49*x^7 + 36*x^6 + 25*x^5 + 16*x^4 + 9*x^3 + 4*x^2 + x

        Sometime you could have to use a special output width and you
        could specify it::

            sage: shell.run_cell('%display ascii_art')
            sage: shell.run_cell('StandardTableaux(4).list()')                          # needs sage.combinat
            [
            [                                                                  1  4    1  3
            [                 1  3  4    1  2  4    1  2  3    1  3    1  2    2       2
            [   1  2  3  4,   2      ,   3      ,   4      ,   2  4,   3  4,   3   ,   4   ,
            <BLANKLINE>
                       1 ]
               1  2    2 ]
               3       3 ]
               4   ,   4 ]
            sage: shell.run_cell('%display ascii_art 50')
            sage: shell.run_cell('StandardTableaux(4).list()')                          # needs sage.combinat
            [
            [
            [                 1  3  4    1  2  4    1  2  3
            [   1  2  3  4,   2      ,   3      ,   4      ,
            <BLANKLINE>
                                                      1 ]
                              1  4    1  3    1  2    2 ]
              1  3    1  2    2       2       3       3 ]
              2  4,   3  4,   3   ,   4   ,   4   ,   4 ]

        As yet another option, typeset mode. This is used in the emacs
        interface::

            sage: shell.run_cell('%display text latex')
            sage: shell.run_cell('1/2')
            1/2

        Switch back::

            sage: shell.run_cell('%display default')

        Switch graphics to default to vector or raster graphics file
        formats::

            sage: shell.run_cell('%display graphics vector')

        TESTS::

            sage: shell.run_cell('%display invalid_mode')
            value must be unset (None) or one of ('plain', 'ascii_art', 'unicode_art', 'latex'), got invalid_mode
            sage: shell.quit()
        """
        from sage.repl.rich_output import get_display_manager
        dm = get_display_manager()
        args = args.strip().split()
        if not args:
            print(dm.preferences)
            return
        arg0 = args[0]
        # deprecated values
        if arg0 == 'simple':
            dm.preferences.text = 'plain'
        elif arg0 == 'typeset':
            dm.preferences.text = 'latex'
        elif arg0 in ['ascii_art', 'unicode_art'] and len(args) > 1:
            try:
                max_width = int(args[1])
            except ValueError:
                max_width = 0
            if max_width <= 0:
                raise ValueError(
                        "max width must be a positive integer")
            import sage.typeset.character_art as character_art
            character_art.MAX_WIDTH = max_width
            dm.preferences.text = arg0
        # Unset all
        elif arg0 in ['default', 'None']:  # un-stringify "%display None"
            for option in map(str, dm.preferences.available_options()):
                delattr(dm.preferences, option)
        # Normal argument handling
        elif arg0 in map(str, dm.preferences.available_options()) and len(args) <= 2:
            if len(args) == 1:
                # "%display text" => get current value
                print(getattr(dm.preferences, arg0))
            else:
                # "%display text latex" => set new value
                assert len(args) == 2
                if args[1] in ['default', 'None']:
                    delattr(dm.preferences, arg0)
                else:
                    try:
                        setattr(dm.preferences, arg0, args[1])
                    except ValueError as err:
                        print(err)  # do not show traceback
        # If all else fails: assume text
        else:
            try:
                dm.preferences.text = arg0
            except ValueError as err:
                print(err)  # do not show traceback

    @cell_magic
    def cython(self, line, cell):
        """
        Cython cell magic.

        This is syntactic sugar on the
        :func:`~sage.misc.cython.cython_compile` function.

        INPUT:

        - ``line`` -- parsed as keyword arguments. The allowed arguments are:

          - ``--verbose N`` / ``-v N``
          - ``--compile-message`` / ``-m``
          - ``--use-cache`` / ``-c``
          - ``--create-local-c-file`` / ``-l``
          - ``--annotate`` / ``-a``
          - ``--sage-namespace`` / ``-s``
          - ``--create-local-so-file`` / ``-o``
          - ``--no-compile-message``, ``--no-use-cache``, etc. (there is no short form for the ``--no-*`` flags)

          See :func:`~sage.misc.cython.cython` for details.

        - ``cell`` -- string; the Cython source code to process

        OUTPUT: none; the Cython code is compiled and loaded

        EXAMPLES::

            sage: # needs sage.misc.cython
            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell(
            ....: '''
            ....: %%cython -v1 --annotate --no-sage-namespace
            ....: def f():
            ....:     print('test')
            ....: ''')
            Compiling ....pyx because it changed.
            [1/1] Cythonizing ....pyx
            sage: f()
            test

        TESTS:

        See :mod:`sage.repl.interpreter` for explanation of the dummy line.

        Test unrecognized arguments::

            sage: # needs sage.misc.cython
            sage: shell.run_cell('''
            ....: %%cython --some-unrecognized-argument
            ....: print(1)
            ....: ''')
            UsageError: unrecognized arguments: --some-unrecognized-argument

        Test ``--help`` is disabled::

            sage: # needs sage.misc.cython
            sage: shell.run_cell('''
            ....: %%cython --help
            ....: print(1)
            ....: ''')
            UsageError: unrecognized arguments: --help

        Test invalid quotes::

            sage: # needs sage.misc.cython
            sage: print("dummy line"); shell.run_cell('''
            ....: %%cython --a='
            ....: print(1)
            ....: ''')
            dummy line
            ...
            ValueError...Traceback (most recent call last)
            ...
            ValueError: No closing quotation
        """
        from sage.misc.cython import cython_compile
        import shlex
        import argparse

        class ExitCatchingArgumentParser(argparse.ArgumentParser):
            def error(self, message):
                # exit_on_error=False does not work completely in some Python versions
                # see https://stackoverflow.com/q/67890157
                # we raise UsageError to make the interface similar to what happens when e.g.
                # IPython's ``%run`` gets unrecognized arguments
                from IPython.core.error import UsageError
                raise UsageError(message)

        parser = ExitCatchingArgumentParser(prog="%%cython", add_help=False)
        parser.add_argument("--verbose", "-v", type=int)
        parser.add_argument("-m", "--compile-message", action=argparse.BooleanOptionalAction)
        parser.add_argument("-c", "--use-cache", action=argparse.BooleanOptionalAction)
        parser.add_argument("-l", "--create-local-c-file", action=argparse.BooleanOptionalAction)
        parser.add_argument("-a", "--annotate", action=argparse.BooleanOptionalAction)
        parser.add_argument("-s", "--sage-namespace", action=argparse.BooleanOptionalAction)
        parser.add_argument("-o", "--create-local-so-file", action=argparse.BooleanOptionalAction)
        args = parser.parse_args(shlex.split(line))
        return cython_compile(cell, **{k: v for k, v in args.__dict__.items() if v is not None})

    @cell_magic
    def fortran(self, line, cell):
        """
        Fortran cell magic.

        This is syntactic sugar on the
        :func:`~sage.misc.inline_fortran.fortran` function.

        INPUT:

        - ``line`` -- ignored

        - ``cell`` -- string; the Cython source code to process

        OUTPUT: none; the Fortran code is compiled and loaded

        EXAMPLES::

            sage: # needs numpy
            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('''
            ....: %%fortran
            ....: C FILE: FIB1.F
            ....:       SUBROUTINE FIB(A,N)
            ....: C
            ....: C     CALCULATE FIRST N FIBONACCI NUMBERS
            ....: C
            ....:       INTEGER N
            ....:       REAL*8 A(N)
            ....:       DO I=1,N
            ....:          IF (I.EQ.1) THEN
            ....:             A(I) = 0.0D0
            ....:          ELSEIF (I.EQ.2) THEN
            ....:             A(I) = 1.0D0
            ....:          ELSE
            ....:             A(I) = A(I-1) + A(I-2)
            ....:          ENDIF
            ....:       ENDDO
            ....:       END
            ....: C END FILE FIB1.F
            ....: ''')
            sage: fib
            <fortran ...>
            sage: from numpy import array
            sage: a = array(range(10), dtype=float)
            sage: fib(a, 10)
            sage: a
            array([  0.,   1.,   1.,   2.,   3.,   5.,   8.,  13.,  21.,  34.])
        """
        from sage.misc.inline_fortran import fortran
        return fortran(cell)


class SageCustomizations:

    def __init__(self, shell=None):
        """
        Initialize the Sage plugin.
        """
        self.shell = shell

        self.auto_magics = SageMagics(shell)
        self.shell.register_magics(self.auto_magics)

        self.shell.set_hook('editor', LazyImport("sage.misc.edit_module", "edit_devel"))

        self.init_inspector()
        self.init_line_transforms()

        try:
            import sage.all # until sage's import hell is fixed
        except ImportError:
            import sage.all__sagemath_repl

        self.shell.verbose_quit = True

        self.register_interface_magics()

        if SAGE_IMPORTALL == 'yes':
            self.init_environment()

    def register_interface_magics(self):
        """
        Register magics for each of the Sage interfaces
        """
        from sage.repl.interface_magic import InterfaceMagic
        InterfaceMagic.register_all(self.shell)

    @staticmethod
    def all_globals():
        """
        Return a Python module containing all globals which should be
        made available to the user.

        EXAMPLES::

            sage: from sage.repl.ipython_extension import SageCustomizations
            sage: SageCustomizations.all_globals()
            <module 'sage.all_cmdline' ...>
        """
        try:
            from sage import all_cmdline
        except ImportError:
            from sage import all__sagemath_repl as all_cmdline
        return all_cmdline

    def init_environment(self):
        """
        Set up Sage command-line environment
        """
        # import outside of cell so we don't get a traceback
        from sage.repl.user_globals import initialize_globals
        initialize_globals(self.all_globals(), self.shell.user_ns)
        self.run_init()

    def run_init(self):
        """
        Run Sage's initial startup file.
        """
        try:
            with open(SAGE_STARTUP_FILE) as f:
                self.shell.run_cell(f.read(), store_history=False)
        except OSError:
            pass

    def init_inspector(self):
        # Ideally, these would just be methods of the Inspector class
        # that we could override; however, IPython looks them up in
        # the global :class:`IPython.core.oinspect` module namespace.
        # Thus, we have to monkey-patch.
        import IPython.core.oinspect
        IPython.core.oinspect.getdoc = LazyImport("sage.misc.sageinspect", "sage_getdoc")
        IPython.core.oinspect.getsource = LazyImport("sage.misc.sagedoc", "my_getsource")
        IPython.core.oinspect.find_file = LazyImport("sage.misc.sageinspect", "sage_getfile")
        IPython.core.oinspect.getargspec = LazyImport("sage.misc.sageinspect", "sage_getargspec")

    def init_line_transforms(self):
        """
        Set up transforms (like the preparser).

        TESTS:

        Check that :issue:`31951` is fixed::

             sage: from IPython import get_ipython
             sage: ip = get_ipython()
             sage: ip.input_transformer_manager.check_complete('''  # indirect doctest
             ....: for i in [1 .. 2]:
             ....:     a = 2''')
             ('incomplete', 2)
             sage: ip.input_transformer_manager.check_complete('''
             ....: def foo(L)
             ....:     K.<a> = L''')
             ('invalid', None)
             sage: ip.input_transformer_manager.check_complete('''
             ....: def foo(L):
             ....:     K.<a> = L''')
             ('incomplete', 4)
             sage: ip.input_transformer_manager.check_complete('''
             ....: def foo(L):
             ....:     K.<a> = L''')
             ('incomplete', 4)
             sage: ip.input_transformer_manager.check_complete('''
             ....: def foo(R):
             ....:     a = R.0''')
             ('incomplete', 4)
             sage: ip.input_transformer_manager.check_complete('''
             ....: def foo(a):
             ....:     b = 2a''')
             ('invalid', None)
             sage: implicit_multiplication(True)
             sage: ip.input_transformer_manager.check_complete('''
             ....: def foo(a):
             ....:     b = 2a''')
             ('incomplete', 4)
             sage: ip.input_transformer_manager.check_complete('''
             ....: def foo():
             ....:     f(x) = x^2''')
             ('incomplete', 4)
             sage: ip.input_transformer_manager.check_complete('''
             ....: def foo():
             ....:     2.factor()''')
             ('incomplete', 4)
        """
        from IPython.core.inputtransformer2 import TransformerManager
        from .interpreter import SagePromptTransformer, SagePreparseTransformer

        self.shell.input_transformer_manager.cleanup_transforms.insert(1, SagePromptTransformer)
        self.shell.input_transformers_post.append(SagePreparseTransformer)

        # Create an input transformer that does Sage's special syntax in the first step.
        # We append Sage's preparse to the cleanup step, so that ``check_complete`` recognizes
        # Sage's special syntax.
        # Behaviour is somewhat inconsistent, but the syntax is recognized as desired.
        M = TransformerManager()
        M.token_transformers = self.shell.input_transformer_manager.token_transformers
        M.cleanup_transforms.insert(1, SagePromptTransformer)
        M.cleanup_transforms.append(SagePreparseTransformer)
        self.shell._check_complete_transformer = M
        self.shell.input_transformer_manager.check_complete = M.check_complete


class SageJupyterCustomizations(SageCustomizations):
    @staticmethod
    def all_globals():
        """
        Return a Python module containing all globals which should be
        made available to the user when running the Jupyter notebook.

        EXAMPLES::

            sage: from sage.repl.ipython_extension import SageJupyterCustomizations
            sage: SageJupyterCustomizations.all_globals()
            <module 'sage.repl.ipython_kernel.all_jupyter' ...>
        """
        from .ipython_kernel import all_jupyter
        return all_jupyter


@run_once
def load_ipython_extension(ip):
    """
    Load the extension in IPython.
    """
    # this modifies ip
    SageCustomizations(shell=ip)
