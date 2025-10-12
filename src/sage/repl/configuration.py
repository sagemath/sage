# sage_setup: distribution = sagemath-repl
r"""
Sage's IPython Configuration

TESTS:

We check that Sage stdin can be piped in even if stdout is a tty; In that case
the IPython simple prompt is being used::

    sage: cmd = 'print([sys.stdin.isatty(), sys.stdout.isatty()])'
    sage: import pexpect                                                                # needs pexpect
    sage: output = pexpect.run(                                                         # needs pexpect
    ....:     'bash -c \'echo "{0}" | sage\''.format(cmd),
    ....: ).decode('utf-8', 'surrogateescape')
    sage: 'sage: [False, True]' in output                                               # needs pexpect
    True
"""
# ****************************************************************************
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


import copy
import sys

from traitlets.config.loader import Config

from sage.repl.prompts import SagePrompts

# Name of the Sage IPython extension
SAGE_EXTENSION = 'sage'


class SageIpythonConfiguration:

    def _doctest_mode(self):
        """
        Whether we are in doctest mode.

        This returns ``True`` during doctests.

        EXAMPLES::

            sage: from sage.repl.configuration import sage_ipython_config
            sage: sage_ipython_config._doctest_mode()
            True
        """
        from sage.doctest import DOCTEST_MODE
        return DOCTEST_MODE

    def _allow_ansi(self):
        """
        Whether to allow ANSI escape sequences.

        This returns ``False`` during doctests to avoid ANSI escape
        sequences.

        EXAMPLES::

            sage: from sage.repl.configuration import sage_ipython_config
            sage: sage_ipython_config._allow_ansi()
            False
        """
        return (not self._doctest_mode()) and sys.stdin.isatty() and sys.stdout.isatty()

    def colors(self):
        """
        Return the IPython color palette.

        This returns ``'nocolor'`` during doctests to avoid ANSI escape
        sequences.

        EXAMPLES::

            sage: from sage.repl.configuration import sage_ipython_config
            sage: sage_ipython_config.simple_prompt()
            True
        """
        if not self._allow_ansi():
            return 'nocolor'
        from sage.repl.interpreter import SageTerminalInteractiveShell
        return SageTerminalInteractiveShell.colors.default()

    def simple_prompt(self):
        """
        Return whether to use the simple prompt.

        This returns ``True`` during doctests to avoid ANSI escape sequences.

        EXAMPLES::

            sage: from sage.repl.configuration import sage_ipython_config
            sage: sage_ipython_config.simple_prompt()
            True
        """
        return not self._allow_ansi()

    def term_title(self):
        """
        Return whether to set the terminal title.

        This returns false during doctests to avoid ANSI escape sequences.

        EXAMPLES::

            sage: from sage.repl.configuration import sage_ipython_config
            sage: sage_ipython_config.term_title()
            False
        """
        return self._allow_ansi()

    def default(self):
        """
        Return a new default configuration object.

        EXAMPLES::

            sage: from sage.repl.configuration import sage_ipython_config
            sage: conf = sage_ipython_config.default()
            sage: type(conf)
            <class 'traitlets.config.loader.Config'>
            sage: 'InteractiveShell' in conf
            True
        """
        from sage.repl.interpreter import SageTerminalInteractiveShell

        # Use the same config for both InteractiveShell, and its subclass
        # TerminalInteractiveShell (note: in fact some configs like term_title
        # only apply to the latter, but we can still use the same config for
        # both for simplicity's sake; see Issue #28289)
        InteractiveShell = Config(
            prompts_class=SagePrompts,
            ast_node_interactivity='all',
            colors=self.colors(),
            simple_prompt=self.simple_prompt(),
            term_title=self.term_title(),
            confirm_exit=False,
            separate_in=''
        )

        cfg = Config(
            TerminalIPythonApp=Config(
                display_banner=False,
                verbose_crash=True,
                test_shell=False,
                shell_class=SageTerminalInteractiveShell,
            ),
            InteractiveShell=InteractiveShell,
            TerminalInteractiveShell=InteractiveShell,
            InteractiveShellApp=Config(extensions=[SAGE_EXTENSION]),
            # TODO: jedi is disabled by default because it causes too many troubles
            # disabling issue: https://github.com/sagemath/sage/issues/31648
            # reenabling issue: https://github.com/sagemath/sage/issues/31649
            IPCompleter=Config(use_jedi=False),
        )
        if self._doctest_mode():
            # Using the file-backed history causes problems in parallel tests
            cfg.HistoryManager = Config(hist_file=':memory:')
        return cfg

    def copy(self):
        """
        Return a copy of the current configuration.

        EXAMPLES::

            sage: from sage.repl.configuration import sage_ipython_config
            sage: conf = sage_ipython_config.copy()
            sage: type(conf)
            <class 'traitlets.config.loader.Config'>
            sage: 'InteractiveShell' in conf
            True
        """
        try:
            return copy.deepcopy(get_ipython().config)
        except NameError:
            return self.default()


sage_ipython_config = SageIpythonConfiguration()
