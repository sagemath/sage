r"""
Command line user interface for doctesting

This module provides :class:`DoctestUserInterface`, an implementation of a
:class:`sage.dev.user_interface.UserInterface` on the command line which can be
used for (non-interactive) doctesting.

AUTHORS:

- Julian Rueth: initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import six
from getpass import getpass

from sage.dev.cmd_line_interface import CmdLineInterface, INFO

class DoctestUserInterface(CmdLineInterface, list):
    r"""
    A :class:`sage.dev.user_interface.UserInterface` which can be used for
    doctesting. Whenever the user would normally be prompted for input, the
    answer to the question is taken to be `self.pop()`, i.e., answers can be
    provided by appending them via :meth:`append`.

    EXAMPLES::

        sage: from sage.dev.test.config import DoctestConfig
        sage: from sage.dev.test.user_interface import DoctestUserInterface
        sage: UI = DoctestUserInterface(DoctestConfig()["UI"])
        sage: UI.append("Answer")
        sage: UI._get_input("Question?")
        Question? Answer
        'Answer'
        sage: UI._get_input("Question?")
        Traceback (most recent call last):
        ...
        RuntimeError: no answers left in DoctestUserInterface for question "Question? ".
    """
    def _get_input(self, prompt, options=None, default=None, input_func=six.moves.input):
        r"""
        Overwrites
        :meth:`sage.dev.cmd_line_interface.CmdLineInterface._get_input` for
        doctesting.

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: UI = DoctestUserInterface(DoctestConfig()["UI"])
            sage: UI.append('')
            sage: UI._get_input("Should I delete your home directory?", ("yes","no","maybe"), default=0)
            Should I delete your home directory? [Yes/no/maybe]
            'yes'
            sage: UI._get_input("Should I delete your home directory?", ("yes","no","maybe"), default=0)
            Traceback (most recent call last):
            ...
            RuntimeError: no answers left in DoctestUserInterface for question "Should
            I delete your home directory? [Yes/no/maybe] ".
        """
        old_input_func = input_func
        def input_func(prompt):
            if not self:
                raise RuntimeError('no answers left in DoctestUserInterface for'
                                   ' question "{0}".'.format(prompt))
            ret = self.pop()
            if old_input_func is getpass:
                self.show(prompt)
            else:
                self.show(prompt+ret)
            return ret
        return CmdLineInterface._get_input(self, prompt, options, default, input_func)

    def _show(self, message, log_level, *args):
        r"""
        Overwrites :meth:`sage.dev.cmd_line_interface.CmdLineInterface._show`
        for doctesting. In particular, no line wrapping is performed.

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: UI = DoctestUserInterface(DoctestConfig()["UI"])
            sage: UI.info("I ate {0} for dinner, a whole {1} for dessert, and then took a nice walk around the lake.", 'filet mignon', 'apple pie') # indirect doctest
            #  I ate filet mignon for dinner, a whole apple pie for dessert, and then took a nice walk around the lake.
        """
        if len(args) > 0:
            message = message.format(*args)

        if log_level == INFO:
            prefix = '# '
        else:
            prefix = ''

        message = (prefix + line for line in message.splitlines())

        print(self._color_code(log_level) +
              '\n'.join(message) +
              self._color_code())

    def edit(self, filename):
        r"""
        Overwrites :meth:`sage.dev.cmd_line_interface.CmdLineInterface.edit`
        for doctesting.

        TESTS::

            sage: tmp = tmp_filename()
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: UI = DoctestUserInterface(DoctestConfig()["UI"])
            sage: UI.append("Some\nlines\n")
            sage: UI.edit(tmp)
            sage: print open(tmp,'r').read()
            Some
            lines
            <BLANKLINE>
            sage: UI.edit(tmp)
            Traceback (most recent call last):
            ...
            RuntimeError: no answers left in DoctestUserInterface
            sage: os.unlink(tmp)
        """
        if not self:
            raise RuntimeError("no answers left in DoctestUserInterface")
        open(filename, 'w').write(self.pop())
