#*****************************************************************************
#       Copyright (C) 2009 Carl Witty <Carl.Witty@gmail.com>
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

"""Miscellaneous utility routines used for the interpreter generator."""

from __future__ import print_function, absolute_import

import os
import textwrap

from jinja2 import Environment
from jinja2.runtime import StrictUndefined


# We share a single jinja2 environment among all templating in this
# file.  We use trim_blocks=True (which means that we ignore white
# space after "%}" jinja2 command endings), and set undefined to
# complain if we use an undefined variable.
JINJA_ENV = Environment(trim_blocks=True, undefined=StrictUndefined)

# Allow 'i' as a shorter alias for the built-in 'indent' filter.
JINJA_ENV.filters['i'] = JINJA_ENV.filters['indent']


def je(template, **kwargs):
    r"""
    A convenience method for creating strings with Jinja templates.

    The name je stands for "Jinja evaluate".

    The first argument is the template string; remaining keyword
    arguments define Jinja variables.

    If the first character in the template string is a newline, it is
    removed (this feature is useful when using multi-line templates defined
    with triple-quoted strings -- the first line doesn't have to be on
    the same line as the quotes, which would screw up the indentation).

    (This is very inefficient, because it recompiles the Jinja
    template on each call; don't use it in situations where
    performance is important.)

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters.internal import je
        sage: je("{{ a }} > {{ b }} * {{ c }}", a='"a suffusion of yellow"', b=3, c=7)
        '"a suffusion of yellow" > 3 * 7'
    """
    if template and template[0] == '\n':
        template = template[1:]

    # It looks like Jinja2 automatically removes one trailing newline?
    if template and template[-1] == '\n':
        template = template + '\n'

    tmpl = JINJA_ENV.from_string(template)
    return tmpl.render(kwargs)


def indent_lines(n, text):
    r"""
    Indent each line in text by ``n`` spaces.

    INPUT:

    - ``n`` -- indentation amount
    - ``text`` -- text to indent

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters.internal import indent_lines
        sage: indent_lines(3, "foo")
        '   foo'
        sage: indent_lines(3, "foo\nbar")
        '   foo\n   bar'
        sage: indent_lines(3, "foo\nbar\n")
        '   foo\n   bar\n'
    """
    lines = text.splitlines(True)
    spaces = ' ' * n
    return ''.join((spaces if line.strip() else '') + line
                   for line in lines)


def reindent_lines(n, text):
    r"""
    Strip any existing indentation on the given text (while keeping
    relative indentation) then re-indents the text by ``n`` spaces.

    INPUT:

    - ``n`` -- indentation amount
    - ``text`` -- text to indent

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters.internal import reindent_lines
        sage: print(reindent_lines(3, " foo\n  bar"))
           foo
             bar
    """

    return indent_lines(n, textwrap.dedent(text))


def write_if_changed(fn, value):
    r"""
    Write value to the file named fn, if value is different than
    the current contents.

    EXAMPLES::

        sage: from sage_setup.autogen.interpreters.internal import *
        sage: def last_modification(fn): return os.stat(fn).st_mtime
        sage: fn = tmp_filename('gen_interp')
        sage: write_if_changed(fn, 'Hello, world')
        sage: t1 = last_modification(fn)
        sage: open(fn).read()
        'Hello, world'
        sage: sleep(2)            # long time
        sage: write_if_changed(fn, 'Goodbye, world')
        sage: t2 = last_modification(fn)
        sage: open(fn).read()
        'Goodbye, world'
        sage: sleep(2)            # long time
        sage: write_if_changed(fn, 'Goodbye, world')
        sage: t3 = last_modification(fn)
        sage: open(fn).read()
        'Goodbye, world'
        sage: t1 == t2            # long time
        False
        sage: t2 == t3
        True
    """

    old_value = None
    try:
        with open(fn) as file:
            old_value = file.read()
    except OSError:
        pass

    if value != old_value:
        # We try to remove the file, in case it exists.  This is to
        # automatically break hardlinks... see #5350 for motivation.
        try:
            os.remove(fn)
        except OSError:
            pass

        with open(fn, 'w') as file:
            file.write(value)
