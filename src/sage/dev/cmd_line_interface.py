r"""
Command line user interface

This module provides :class:`CmdLineInterface`, an implementation of a
:class:`sage.dev.user_interface.UserInterface` on the command line.

AUTHORS:

- David Roe, Julian Rueth: initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 David Roe <roed.math@gmail.com>
#                          Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function
from subprocess import check_call, CalledProcessError
from getpass import getpass

import os
import textwrap
import itertools

from user_interface import UserInterface
from user_interface import ERROR, WARNING, NORMAL, INFO, DEBUG

try:
    from sage.doctest import DOCTEST_MODE
except ImportError:
    DOCTEST_MODE = False


class CmdLineInterface(UserInterface):
    r"""
    An implementation of a :class:`sage.dev.user_interface.UserInterface` on
    the command line.

    EXAMPLES::

        sage: from sage.dev.test.config import DoctestConfig
        sage: from sage.dev.cmd_line_interface import CmdLineInterface
        sage: CmdLineInterface(DoctestConfig()["UI"])
        CmdLineInterface()
    """
    def __repr__(self):
        r"""
        Return a printable representation of this object.

        EXAMPLES::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.cmd_line_interface import CmdLineInterface
            sage: CmdLineInterface(DoctestConfig()["UI"])
            CmdLineInterface()
        """
        return "CmdLineInterface()"

    def _std_values(self, prompt, options, default):
        r"""
        Helper method to generate prompts.

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.cmd_line_interface import CmdLineInterface
            sage: UI = CmdLineInterface(DoctestConfig()["UI"])
            sage: UI._std_values("Should I delete your home directory?",
            ....:                ("yes", "no", "maybe"), default=1)
            ('Should I delete your home directory? [yes/No/maybe] ', ('yes', 'no', 'maybe'), 'no')
        """
        if options is not None:
            options = tuple(opt.lower() for opt in options)
            if options:
                prompt += " ["
                prompt += "/".join(opt if i != default else opt.capitalize()
                                    for i, opt in enumerate(options))
                prompt += "]"
                if default is not None:
                    default = options[default]
            else:
                options = None
        prompt = self._color_code('prompt') + prompt + self._color_code() + " "
        return prompt, options, default

    def _get_input(self, prompt, options=None, default=None, input_func=raw_input):
        r"""
        Helper method for :meth:`switch`, :meth:`get_input`, and :meth:`get_password`.

        INPUT:

        - ``prompt`` -- a string

        - ``options`` -- a list of strings or ``None`` (default: ``None``)

        - ``default`` -- an integer or ``None`` (deault: ``None``), the default
          option as an index into ``options``

        - ``input_func`` -- a function to get input from user (default:
          ``raw_input``)

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.cmd_line_interface import CmdLineInterface
            sage: UI = CmdLineInterface(DoctestConfig()["UI"])
            sage: def input_func(prompt):
            ....:     print(prompt)
            ....:     return "no"
            sage: UI._get_input("Should I delete your home directory?",
            ....:     ("yes","no","maybe"), default=0, input_func = input_func)
            Should I delete your home directory? [Yes/no/maybe]
            'no'
        """
        try:
            prompt, options, default = self._std_values(prompt, options, default)
            while True:
                s = input_func(prompt)
                if options is None:
                    return s
                if len(s.strip()) == 0:
                    if default is None:
                        self.show("Please enter an option.")
                        continue
                    else:
                        return default

                itr = (opt for opt in options if opt.lower().startswith(s.lower()))
                try:
                    ret = next(itr)
                except StopIteration:
                    self.show("Please specify an allowable option.")
                    continue

                try:
                    ret = next(itr)
                    self.show("Please disambiguate between options.")
                except StopIteration:
                    return ret
        except KeyboardInterrupt:
            from user_interface_error import OperationCancelledError
            raise OperationCancelledError("cancelled by keyboard interrupt")

    def select(self, prompt, options, default=None):
        r"""
        Ask user to select from ``options`` and return selected.

        INPUT:

        - ``prompt`` -- a string, the prompt to display

        - ``options`` -- an iterable that contains the options

        - ``default`` -- an integer or ``None`` (default: ``None``), the
          default option as an index into ``options``

        .. SEEALSO::

            :meth:`sage.dev.user_interface.UserInterface.select`

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.cmd_line_interface import CmdLineInterface
            sage: UI = CmdLineInterface(DoctestConfig()["UI"])
            sage: UI.select("Should I delete your home directory?",    # not tested
            ....:           ("yes", "no", "maybe"), default=2)
            Should I delete your home directory? [yes/no/Maybe] m
            'maybe'

            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: UI = DoctestUserInterface(DoctestConfig()["UI"])
            sage: UI.append("n")
            sage: UI.select("Should I delete your home directory?",    # indirect doctest
            ....:           ("yes","no","maybe"), default=2)
            Should I delete your home directory? [yes/no/Maybe] n
            'no'
        """
        return self._get_input(prompt, options, default)

    def get_input(self, prompt):
        r"""
        Prompt for input.

        INPUT:

        - ``prompt`` -- a string, the prompt to display

        .. SEEALSO::

            :meth:`sage.dev.user_interface.UserInterface.get_input`

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.cmd_line_interface import CmdLineInterface
            sage: UI = CmdLineInterface(DoctestConfig()["UI"])
            sage: UI.get_input("What do you want for dinner?") # not tested
            What do you want for dinner? filet mignon
            'filet mignon'

            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: UI = DoctestUserInterface(DoctestConfig()["UI"])
            sage: UI.append("filet mignon")
            sage: UI.get_input("What do you want for dinner?") # indirect doctest
            What do you want for dinner? filet mignon
            'filet mignon'
        """
        return self._get_input(prompt)

    def get_password(self, prompt):
        r"""
        Prompt for password.

        INPUT:

        - ``prompt`` -- a string, the prompt to display

        .. SEEALSO::

            :meth:`sage.dev.user_interface.UserInterface.get_password`

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.cmd_line_interface import CmdLineInterface
            sage: UI = CmdLineInterface(DoctestConfig()["UI"])
            sage: UI.get_password("What is the key combo for your safe?") # not tested
            What is the key combo for your safe?
            '9247'

            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: UI = DoctestUserInterface(DoctestConfig()["UI"])
            sage: UI.append('9247')
            sage: UI.get_password("What is the key combo for your safe?") # indirect doctest
            What is the key combo for your safe?
            '9247'
        """
        return self._get_input(prompt, input_func=getpass)

    def _get_dimensions(self):
        r"""
        Return the dimensions of the terminal.

        OUTPUT:

        A pair of ints, the number of rows and columns of the terminal.

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.cmd_line_interface import CmdLineInterface
            sage: UI = CmdLineInterface(DoctestConfig()["UI"])
            sage: UI._get_dimensions()
            (25, 80)
        """
        dim = self._ioctl_GWINSZ(0) or self._ioctl_GWINSZ(1) or self._ioctl_GWINSZ(2)
        if dim is None:
             fd = os.open(os.ctermid(), os.O_RDONLY)
             try:
                dim = self._ioctl_GWINSZ(fd)
             finally:
                os.close(fd)
        if dim is None:
            ret = (0,)
        else:
            ret = tuple(int(x) for x in dim)
        if all(ret):
            return ret
        else:
            # fallback values
            return (25, 80)

    def _ioctl_GWINSZ(self, fd):
        r"""
        Return the window size of the terminal at the file descriptor ``fd``.

        OUTPUT:

        A pair of ints, the terminal's rows and columns. If the size could not
        be determined, returns ``None``.

        TESTS::

            sage: import os
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.cmd_line_interface import CmdLineInterface
            sage: CmdLineInterface(DoctestConfig()["UI"])._ioctl_GWINSZ(0)
            (25, 80)
        """
        if DOCTEST_MODE:
            return (25, 80)
        try:
            import struct
            import fcntl
            import termios
        except ImportError:
            return None

        try:
            return struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
        except IOError:
            return None

    def _color_code(self, color=None):
        """
        Return an ansi color code.

        INPUT:

        - ``color`` -- ``None``, ``'prompt'``, or one of the constants
          ``ERROR``, ``WARNING``, ``NORMAL``, ``INFO``, or ``DEBUG``.

        OUTPUT:

        String, possibly containing a color code.
        """
        if DOCTEST_MODE:
            return ''
        if color is None:
            return '\033[0m'
        elif color == ERROR:
            return '\033[0;31m'
        elif color == WARNING:
            return '\033[0;33m'
        elif color == INFO or color == DEBUG:
            return '\033[0;36m'
        elif color == 'prompt':
            return '\033[0;34m'
        else:
            return ''

    def _show(self, message, log_level, *args):
        r"""
        Display ``message``.

        .. SEEALSO::

            :meth:`sage.dev.user_interface.UserInterface.show`

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.cmd_line_interface import CmdLineInterface
            sage: UI = CmdLineInterface(DoctestConfig()["UI"])
            sage: UI.info("I ate {0} for dinner, a whole {1} for dessert, and then took a nice walk around the lake.", 'filet mignon', 'apple pie') # indirect doctest
            #  I ate filet mignon for dinner, a whole apple pie for dessert, and then took a
            #  nice walk around the lake.
        """
        if len(args) > 0:
            message = message.format(*args)
        height, width = self._get_dimensions()
        kwds = {'width': width}
        if log_level == INFO:
            kwds['initial_indent'] = kwds['subsequent_indent'] = '#  '

        wrap = textwrap.TextWrapper(**kwds).wrap
        message = list(itertools.chain.from_iterable(
            wrap(line) for line in message.splitlines()))

        if len(message) <= height:
            print(self._color_code(log_level) +
                  '\n'.join(message) +
                  self._color_code())
        else:
            message = '\n'.join(message)+'\n'
            try:
                self._pager(message)
            except AttributeError:
                import pydoc
                self._pager = pydoc.getpager()
                self._pager(message)

    def edit(self, filename):
        r"""
        Drop user into editor with filename open.

        .. SEEALSO::

            :meth:`sage.dev.user_interface.UserInterface.edit`

        TESTS::

            sage: tmp = tmp_filename()
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.cmd_line_interface import CmdLineInterface
            sage: UI = CmdLineInterface(DoctestConfig()["UI"])
            sage: UI.edit(tmp) # not tested
            sage: print open(tmp,'r').read() # not tested
            Some
            lines
            <BLANKLINE>
            sage: os.unlink(tmp)
        """
        try:
            editor = os.environ.get('EDITOR', 'nano')
            check_call(['sage-native-execute', editor, filename])
        except CalledProcessError:
            from user_interface_error import OperationCancelledError
            raise OperationCancelledError("Editor returned non-zero exit value")
