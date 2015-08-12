r"""
User Interface

This module provides an abstract base class that is used for displaying messages
and prompting for user input.

.. SEEALSO::

    :class:`CmdLineInterface` for an implementation of this class

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

import six

# log levels
ERROR = -2
WARNING = -1
NORMAL = 0
INFO = 1 # display informational messages, such as git commands executed
DEBUG = 2 # display additional debug information

class UserInterface(object):
    r"""
    An abstract base class for displaying messages and prompting the user for
    input.

    TESTS::

        sage: from sage.dev.user_interface import UserInterface
        sage: from sage.dev.test.config import DoctestConfig
        sage: UserInterface(DoctestConfig())
        <sage.dev.user_interface.UserInterface object at 0x...>
    """

    def __init__(self, config):
        r"""
        Initialization.

        TESTS:

            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: type(UserInterface(DoctestConfig()))
            <class 'sage.dev.user_interface.UserInterface'>
        """
        self._config = config

    def select(self, prompt, options, default=None):
        r"""
        Ask the user to select from list of options and return selected option.

        INPUT:

        - ``prompt`` -- a string, prompt to display

        - ``options`` -- iterable of strings, the options

        - ``default`` -- an integer or ``None`` (default: ``None``),
          the index of the default option

        TESTS::

            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: UI = UserInterface(DoctestConfig())
            sage: UI.select("Should I delete your home directory?", ("yes","no","maybe"), default=1)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def confirm(self, question, default=None):
        r"""
        Ask a yes/no question and return the response as a boolean.

        INPUT:

        - ``question`` -- a string

        - ``default`` -- a boolean or ``None`` (default: ``None``), the
          default value

        TESTS::

            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: UI = UserInterface(DoctestConfig())
            sage: UI.confirm("Should I delete your home directory?")
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if default is None:
            default_option = None
        elif default is True:
            default_option = 0
        elif default is False:
            default_option = 1
        else:
            raise ValueError

        return self.select(question, ("yes", "no"), default_option) == "yes"

    def get_input(self, prompt):
        r"""
        Read input after displaying ``prompt``.

        TESTS::

            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: UI = UserInterface(DoctestConfig())
            sage: UI.get_input("What do you want for dinner?")
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def get_password(self, prompt):
        r"""
        Ask for a password after displaying ``prompt``.

        TESTS::

            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: UI = UserInterface(DoctestConfig())
            sage: UI.get_password("What is the passphrase for your safe?")
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def show(self, message, *args, **kwds):
        r"""
        Display ``message``.

        INPUT:

        - ``message`` -- a string or list/tuple/iterable of
          strings. Each individual message will be wrapped to the
          terminal width.

        - ``*args`` -- optional positional arguments that will be
          passed to ``message.format()``.

        - ``log_level=INFO`` -- one of ``ERROR``, ``WARNING``,
          ``NORMAL``, ``INFO``, or ``DEBUG`` (default:
          ``INFO``). Optional keyword argument.

        TESTS::

            sage: from sage.dev.user_interface import UserInterface, DEBUG
            sage: from sage.dev.test.config import DoctestConfig
            sage: UI = UserInterface(DoctestConfig())
            sage: UI.show("I ate filet mignon for dinner.")
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: UI.show("I ate filet mignon for dinner.", log_level=DEBUG)
        """
        log_level = kwds.pop('log_level', NORMAL)
        if log_level not in (ERROR, WARNING, NORMAL, INFO, DEBUG):
            raise ValueError('log_level is invalid')
        if len(kwds) != 0:
            raise ValueError('the only allowed keyword argument is "log_level"')
        config_log_level = int(self._config.get("log_level", str(INFO)))
        if config_log_level >= log_level:
            if isinstance(message, six.string_types):
                self._show(message, log_level, *args)
            else:
                for msg in message:
                    self._show(msg, log_level, *args)

    def debug(self, message, *args):
        r"""
        Display ``message``.

        INPUT:

        - ``message`` -- a string or list/tuple/iterable of strings.

        - ``*args`` -- optional positional arguments that will be
          passed to ``message.format()``.

        TESTS:

        Debug messages are not displayed in doctests::

            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: UI = UserInterface(DoctestConfig())
            sage: UI.debug("I ate filet mignon for dinner.")
        """
        self.show(message, *args, log_level=DEBUG)

    def info(self, message, *args):
        r"""
        Display ``message``.

        These are informational messages to be shown to the user,
        typically indicating how to proceed.

        INPUT:

        - ``message`` -- a string or list/tuple/iterable of strings.

        - ``*args`` -- optional positional arguments that will be
          passed to ``message.format()``.

        TESTS:

        Info messages are not displayed in doctests::

            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: UI = UserInterface(DoctestConfig())
            sage: UI.info("I ate filet mignon for dinner.")
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        self.show(message, *args, log_level=INFO)

    def warning(self, message, *args):
        r"""
        Display ``message``.

        INPUT:

        - ``message`` -- a string or list/tuple/iterable of strings.

        - ``*args`` -- optional positional arguments that will be
          passed to ``message.format()``.

        TESTS:

            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: UI = UserInterface(DoctestConfig())
            sage: UI.warning("I ate filet mignon for dinner.")
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        self.show(message, *args, log_level=WARNING)

    def error(self, message, *args):
        r"""
        Display ``message``.

        INPUT:

        - ``message`` -- a string or list/tuple/iterable of strings.

        - ``*args`` -- optional positional arguments that will be
          passed to ``message.format()``.

        TESTS:

            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: UI = UserInterface(DoctestConfig())
            sage: UI.error("I ate filet mignon for dinner.")
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        self.show(message, *args, log_level=ERROR)

    def _show(self, message, log_level, *args):
        r"""
        Display ``message``.

        INPUT:

        - ``message`` -- a string

        - ``log_level`` -- one of ``ERROR``, ``WARNING``, ``NORMAL``,
          ``INFO``, or ``DEBUG`` (default: ``NORMAL``)

        - ``*args`` -- optional positional arguments that will be
          passed to ``message.format()``.

        TESTS::

            sage: from sage.dev.user_interface import UserInterface, NORMAL
            sage: from sage.dev.test.config import DoctestConfig
            sage: UI = UserInterface(DoctestConfig())
            sage: UI._show("I ate filet mignon for dinner.", NORMAL)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def edit(self, filename):
        r"""
        Drop user into an editor with ``filename`` open.

        OUTPUT:

        Raises a :class:`sage.dev.user_interface_error.OperationCancelledError`
        if the editor exits with non-zero exit code.

        TESTS::

            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: UI = UserInterface(DoctestConfig())
            sage: UI.edit("filename")
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
