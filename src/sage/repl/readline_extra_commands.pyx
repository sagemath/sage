r"""
Extra readline commands

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
``$DOT_SAGE/ipython-*/profile_sage/ipython_config.py``.  Note that ``$DOT_SAGE`` is ``$HOME/.sage``
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

# functions and variables declared in "readline.h"
# of the GNU readline library
ctypedef int rl_hook_func_t()
ctypedef int rl_command_func_t(int,int)

cdef extern int rl_newline(int,int)
cdef extern int rl_get_next_history(int,int)
cdef extern int rl_get_previous_history(int,int)
cdef extern int rl_replace_line(char *,int)
cdef extern int rl_add_defun(char *,rl_command_func_t *,int)
cdef extern void rl_ding()
cdef extern int rl_readline_state
cdef extern rl_command_func_t *rl_last_func
cdef extern char *rl_line_buffer
cdef extern int rl_point, rl_end, rl_mark
cdef extern rl_hook_func_t *rl_startup_hook

DEF RL_STATE_SEARCH = 0x000200

# functions, datatypes and variables declared in "history.h"
# of the GNU readline library
ctypedef void *histdata_t
ctypedef struct HIST_ENTRY:
    char *line
    char *timestamp
    histdata_t data

cdef extern HIST_ENTRY *current_history()
cdef extern int history_is_stifled()
cdef extern int where_history()
cdef extern int history_set_pos(int)
cdef extern int history_search_prefix(char *,int)
cdef extern int history_search(char *,int)
cdef extern int history_length
cdef extern int history_max_entries

#
# The operate-and-get-next command
#

cdef rl_hook_func_t *old_rl_startup_hook = <rl_hook_func_t *> NULL
cdef int saved_history_line_to_use = -1

cdef int set_saved_history():
    """
    Custom startup hook for the operate-and-get-next command
    """
    global rl_startup_hook
    global history_length
    global old_rl_startup_hook
    global saved_history_line_to_use

    if saved_history_line_to_use >= 0:
        rl_get_previous_history(history_length - saved_history_line_to_use, 0)
    saved_history_line_to_use = -1
    rl_startup_hook = old_rl_startup_hook
    return 0

cdef int operate_and_get_next(int count, int key):
    """
    Accepts the current input line and fetches the next line from the history.
    """
    global rl_startup_hook
    global old_rl_startup_hook
    global history_length
    global history_max_entries
    global saved_history_line_to_use

    cdef int where

    # Accept the current line.
    rl_newline(1, key)

    # Find the current line, and find the next line to use.
    where = where_history()

    if history_is_stifled() and history_length >= history_max_entries or where >= history_length - 1:
        saved_history_line_to_use = where
    else:
        saved_history_line_to_use = where + 1

    old_rl_startup_hook = rl_startup_hook
    rl_startup_hook = &set_saved_history

    return 0

#
# The history-search-backward(forward)-and-save commands
#

cdef int history_search_pos = -1
_history_search_string = ""

# DIR < 0 means to search backwards through the history list, DIR >= 0 means to search forward.
cdef int history_search_internal(int count, int dir):
    global _history_search_string
    global history_search_pos

    global rl_readline_state

    cdef int pos, oldpos, ret
    cdef char *temp

    # Search COUNT times through the history for a line whose prefix
    # matches _history_search_string.  When this loop finishes, TEMP,
    # if non-null, is the history line to copy into the line buffer.
    temp = NULL
    oldpos = where_history()
    while count > 0:
        pos = history_search_pos + dir

        if history_set_pos(pos) == 0:
            break

        rl_readline_state |= RL_STATE_SEARCH # readline is in search state
        ret = history_search_prefix(_history_search_string,dir)
        rl_readline_state &= RL_STATE_SEARCH # readline is not in search state

        if ret == 0:
            history_search_pos = where_history()
            temp = current_history().line
        else:
            history_set_pos(oldpos)
            break
        count = count - 1

    #  If we didn't find anything at all, return.
    if temp == NULL:
        rl_ding()
        return 1
    else: # Copy the line we found into the current line buffer.
        rl_replace_line(temp,1)
        return 0

cdef int history_search_forward_and_save(int count, int key):
    """
    Search forward through the history for the string of characters
    from the start of the input line to the cursor position.
    """
    global rl_last_func
    global rl_point
    global rl_line_buffer
    global history_search_pos
    global _history_search_string

    cdef int ret

    if count == 0: return 0

    # initialization for the first time
    if rl_last_func != history_search_forward_and_save and rl_last_func != history_search_backward_and_save:
        if rl_point > 0:
            history_search_pos = where_history()
        _history_search_string = rl_line_buffer[:rl_point]

    # step or search forward through history
    history_set_pos(history_search_pos)
    if len(_history_search_string) == 0:
        ret = rl_get_next_history(count, key)
        history_search_pos = where_history()
    else:
        ret = history_search_internal(abs(count), 1 if count > 0 else -1)

    return ret

cdef int history_search_backward_and_save(int count,int key):
    """
    Search backward through the history for the string of characters
    from the start of the input line to the cursor position.
    """
    global rl_last_func
    global rl_point
    global rl_line_buffer
    global history_search_pos
    global _history_search_string

    cdef int ret

    if count == 0: return 0

    # initialization for the first time
    if rl_last_func != history_search_forward_and_save and rl_last_func != history_search_backward_and_save:
        history_search_pos = where_history()
        _history_search_string = rl_line_buffer[:rl_point]

    # step or search backward through history
    history_set_pos(history_search_pos)
    if len(_history_search_string) == 0:
        ret = rl_get_previous_history(count, key)
        history_search_pos = where_history()
    else:
        ret = history_search_internal(abs(count), -1 if count > 0 else 1)

    return ret

# Install the commands, unbound.
#
# Recommended key bindings are to bind the operate-and-get-next with the control-o key,
# and the history-search-backward-and-save with the up-arrow key, and
# the history-search-forward-and-save with the down-arrow key.
#
rl_add_defun("operate-and-get-next", &operate_and_get_next, -1)
rl_add_defun("history-search-backward-and-save", &history_search_backward_and_save, -1)
rl_add_defun("history-search-forward-and-save", &history_search_forward_and_save, -1)

