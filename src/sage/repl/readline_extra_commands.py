r"""
Extra Readline Commands

.. WARNING::

    The feature described here is no longer available in Sage, as IPython upon which
    Sage's command line interface is based adopted prompt_toolkit as a replacement
    of readline as of IPython version 5.0

The following extra readline commands are available in Sage:

    - ``operate-and-get-next``
    - ``history-search-backward-and-save``
    - ``history-search-forward-and-save``

The ``operate-and-get-next`` command accepts the input line and fetches the next line
from the history. This is the same command with the same name in the Bash shell.

The ``history-search-backward-and-save`` command searches backward in the history
for the string of characters from the start of the input line to the current cursor
position, and fetches the first line found. If the cursor is at the start of the line, the previous line
is fetched. The position of the fetched line is saved internally, and the next search begins at the
saved position.

The ``history-search-forward-and-save`` command behaves similarly but forward.

The previous two commands is best used in tandem to fetch a block of lines from the history,
by searching backward the first line of the block and then issuing the forward command as many times as needed.
They are intended to replace the ``history-search-backward`` command and the ``history-search-forward`` command
provided by the GNU readline library used in Sage.

To bind these commands with keys, insert the relevant lines into the IPython configuration file
``$DOT_SAGE/ipython-*/profile_default/ipython_config.py``.  Note that ``$DOT_SAGE`` is ``$HOME/.sage``
by default. For example,
::

    c = get_config()

    c.InteractiveShell.readline_parse_and_bind = [
        '"\C-o": operate-and-get-next',
        '"\e[A": history-search-backward-and-save',
        '"\e[B": history-search-forward-and-save'
        ]

binds the three commands with the control-o key, the up arrow key, and the down arrow key,
respectively. *Warning:* Sometimes, these keys may be bound to do other actions by the terminal and does not
reach to the readline properly (check this by running ``stty -a`` and reading the ``cchars`` section). Then
you may need to turn off these bindings before the new readline commands work fine . A prominent case is when
control-o is bound to ``discard`` by the terminal. You can turn this off by running ``stty discard undef``.

AUTHORS:

    - Kwankyu Lee (2010-11-23): initial version
    - Kwankyu Lee (2013-06-05): updated for the new IPython configuration format.
"""

#*****************************************************************************
#       Copyright (C) 2010 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.superseded import deprecation
deprecation(21342, "This module and the feature it provides is not available anymore in Sage")
