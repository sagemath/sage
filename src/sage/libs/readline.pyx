"""
Readline

This is the library behind the command line input, it takes keypresses
until you hit Enter and then returns it as a string to Python. We hook
into it so we can make it redraw the input area.

EXAMPLES::

    sage: from sage.libs.readline import *
    sage: replace_line('foobar', 0)
    sage: set_point(3)
    sage: print 'current line:', repr(copy_text(0, get_end()))
    current line: 'foobar'
    sage: print 'cursor position:', get_point()
    cursor position: 3

When printing with :class:`interleaved_output` the prompt and current
line is removed::

    sage: with interleaved_output():
    ....:     print 'output'
    ....:     print 'current line:', repr(copy_text(0, get_end()))
    ....:     print 'cursor position:', get_point()
    output
    current line: ''
    cursor position: 0

After the interleaved output, the line and cursor is restored to the
old value::

    sage: print 'current line:', repr(copy_text(0, get_end()))
    current line: 'foobar'
    sage: print 'cursor position:', get_point()
    cursor position: 3

Finally, clear the current line for the remaining doctests::

    sage: replace_line('', 1)
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun  <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************



cdef extern from 'readline/readline.h':
    int rl_forced_update_display()
    int rl_redisplay()
    int rl_message(char* msg)
    int rl_clear_message()
    void rl_replace_line(char* line, int pos)
    char* rl_copy_text(int begin, int end)
    void rl_save_prompt()
    void rl_restore_prompt()
    int rl_read_key()
    int rl_set_signals()
    int rl_clear_signals()
    int rl_crlf()
    int rl_initialize()
    int rl_catch_signals
    int rl_catch_sigwinch
    int rl_point
    int rl_end


def print_status():
    """
    Print readline status for debug purposes

    EXAMPLES::

        sage: from sage.libs.readline import print_status
        sage: print_status()
        catch_signals: 1
        catch_sigwinch: 1
    """
    print 'catch_signals:', rl_catch_signals
    print 'catch_sigwinch:', rl_catch_sigwinch

def set_signals():
    """
    Install the readline signal handlers

    Install Readline's signal handler for SIGINT, SIGQUIT, SIGTERM,
    SIGALRM, SIGTSTP, SIGTTIN, SIGTTOU, and SIGWINCH, depending on the
    values of rl_catch_signals and rl_catch_sigwinch.

    EXAMPLES::

        sage: from sage.libs.readline import set_signals
        sage: set_signals()
        0
    """
    return rl_set_signals()


def clear_signals():
    """
    Remove the readline signal handlers

    Remove all of the Readline signal handlers installed by
    :func:`set_signals`

    EXAMPLES::

        sage: from sage.libs.readline import clear_signals
        sage: clear_signals()
        0
    """
    return rl_clear_signals()

def get_point():
    """
    Get the cursor position

    OUTPUT:

    Integer

    EXAMPLES::

        sage: from sage.libs.readline import get_point, set_point
        sage: get_point()
        0
        sage: set_point(5)
        sage: get_point()
        5
        sage: set_point(0)
    """
    return rl_point

def get_end():
    """
    Get the end position of the current input

    OUTPUT:

    Integer

    EXAMPLES::

        sage: from sage.libs.readline import get_end
        sage: get_end()
        0
    """
    return rl_end

def set_point(point):
    """
    Set the cursor position

    INPUT:

    - ``point`` -- integer. The new cursor position.

    EXAMPLES::

        sage: from sage.libs.readline import get_point, set_point
        sage: get_point()
        0
        sage: set_point(5)
        sage: get_point()
        5
        sage: set_point(0)
    """
    global rl_point
    rl_point = point

def forced_update_display():
    """
    Force the line to be updated and redisplayed, whether or not
    Readline thinks the screen display is correct.

    EXAMPLES::

        sage: from sage.libs.readline import forced_update_display
        sage: forced_update_display()
        0
    """
    return rl_forced_update_display()

def copy_text(pos_start, pos_end):
    """
    Return a copy of the text between start and end in the current line.

    INPUT:

    - ``pos_start``, ``pos_end`` -- integer. Start and end position.

    OUTPUT:

    String.

    EXAMPLES::

        sage: from sage.libs.readline import copy_text, replace_line
        sage: replace_line('foobar', 0)
        sage: copy_text(1, 5)
        'ooba'
    """
    return rl_copy_text(pos_start, pos_end)

def replace_line(text, clear_undo):
    """
    Replace the contents of rl_line_buffer with text.

    The point and mark are preserved, if possible.

    INPUT:

    - ``text`` -- the new content of the line.

    - ``clear_undo`` -- integer. If non-zero, the undo list associated
      with the current line is cleared.

    EXAMPLES::

        sage: from sage.libs.readline import copy_text, replace_line
        sage: replace_line('foobar', 0)
        sage: copy_text(1, 5)
        'ooba'
    """
    rl_replace_line(text, clear_undo)

def initialize():
    """
    Initialize or re-initialize Readline’s internal state. It’s not
    strictly necessary to call this; readline() calls it before
    reading any input.

    EXAMPLES::

        sage: from sage.libs.readline import initialize
        sage: initialize()
        0
    """
    return rl_initialize()



class interleaved_output:

    def __init__(self):
        """
        Context manager for asynchronous output

        This allows you to show output while at the readline
        prompt. When the block is left, the prompt is restored even if
        it was clobbered by the output.

        EXAMPLES::

            sage: from sage.libs.readline import interleaved_output
            sage: with interleaved_output():
            ....:     print 'output'
            output
        """
        pass

    def __enter__(self):
        """
        Called when entering the with block

        EXAMPLES::

            sage: from sage.libs.readline import interleaved_output
            sage: with interleaved_output():
            ....:     print 'output'
            output
        """
        self._saved_point = rl_point;
        self._saved_line = rl_copy_text(0, rl_end)
        rl_save_prompt()
        rl_replace_line('', 0)
        rl_redisplay()
        rl_clear_signals()

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Called when entering the with block

        EXAMPLES::

            sage: from sage.libs.readline import interleaved_output
            sage: with interleaved_output():
            ....:     print 'output'
            output
        """
        rl_set_signals()
        rl_replace_line(self._saved_line, 0)
        global rl_point
        rl_point = self._saved_point
        rl_restore_prompt()
        rl_forced_update_display()
        return False





