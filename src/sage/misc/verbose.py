# sage_setup: distribution = sagemath-objects
r"""
Verbosity System and Logging in SageMath

Howto: Logging
==============

Using Python's Logging Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Import it::

    sage: import logging
    sage: logging.basicConfig()  # only needed once

Setting the level::

    sage: logging.getLogger().setLevel(logging.INFO)

Log something::

    sage: logger = logging.getLogger(__name__)
    sage: logger.info('Hello. I am talking to you.')
    INFO:__main__:Hello. I am talking to you.

If we haven't set the logging level to ``logging.INFO``, then the previous
wouldn't have been shown.
::

    sage: logger.debug('Hello. I am really talking a lot.')

The latter is not shown as the current logging level is only
``logging.INFO`` and not ``logging.DEBUG``.

Reset the level::

    sage: logging.getLogger().setLevel(logging.WARNING)

Warnings are still shown at this default level (``logging.WARNING``)::

    sage: logger.warning('Hello. I am warning you.')
    WARNING:__main__:Hello. I am warning you.

And that's all.

There are a lot more features, see
:python:`Logging facility for Python<library/logging.html>`.


Using SageMath's Verbosity System
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, this module provides
:func:`verbose`, :func:`set_verbose`, :func:`get_verbose` which can
be used as follows::

    sage: from sage.misc.verbose import verbose, set_verbose, get_verbose
    sage: set_verbose(1)
    sage: t = verbose("This is SageMath.", level=0)
    verbose 0 (<module>) This is SageMath.
    sage: t = verbose("This is SageMath.", level=1)
    verbose 1 (<module>) This is SageMath.
    sage: t = verbose("This is SageMath.", level=2)


Logging Levels of SageMath and Python
=====================================

.. csv-table::
    :class: contentstable
    :widths: 20, 20
    :delim: |

    SageMath | Python
    `-2` | ``logging.CRITICAL``
    `-1` | ``logging.ERROR``
    `0` | ``logging.WARNING``
    `1` | ``logging.INFO``
    `2` | ``logging.DEBUG``


Various
=======

AUTHORS:

- Daniel Krenn (2016)


Functions
=========
"""
# *****************************************************************************
#       Copyright (C) 2006, 2007 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 Gonzalo Tornaria
#       Copyright (C) 2008 John H. Palmieri
#       Copyright (C) 2009 Mike Hansen
#       Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#       Copyright (C) 2018 Frédéric Chapoton
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import sys
import os

LEVEL = 0  # default

verbose_files = []


def verbose(mesg='', t=0, level=1, caller_name=None):
    """
    Print a message if the current verbosity is at least level.

    INPUT:

    - ``mesg`` -- string; a message to print

    - ``t`` -- integer (optional); if included, will also print ``cputime(t)``,
      which is the time since time ``t``. Thus ``t`` should have been obtained
      with ``t=cputime()``

    - ``level`` -- integer (default: 1); the verbosity level of
      what we are printing

    - ``caller_name`` -- string (default: ``None``); the name
      of the calling function. In most cases Python can deduce this, so
      it need not be provided.

    OUTPUT: possibly prints a message to stdout; also returns ``cputime()``

    EXAMPLES::

        sage: set_verbose(1)
        sage: t = cputime()
        sage: t = verbose("This is Sage.", t, level=1, caller_name="william")       # not tested
        VERBOSE1 (william): This is Sage. (time = 0.0)
        sage: set_verbose(0)
    """
    from sage.misc.timing import cputime
    if level > LEVEL:
        return cputime()

    frame = sys._getframe(1).f_code
    file_name = frame.co_filename
    lineno = frame.co_firstlineno
    if 'all' in verbose_files or level <= 0:
        show = True
    else:
        show = False
        for X in verbose_files:
            if file_name.find(X) != -1:
                show = True
                break

    if not show:
        return cputime()

    if t != 0 and mesg == "":
        mesg = "Finished."

    # see recipe 14.7 in Python Cookbook
    if caller_name is None:
        caller_name = frame.co_name
        if caller_name == "?: ":
            caller_name = ""
    short_file_name = os.path.split(frame.co_filename)[1]
    if '<' in short_file_name and '>' in short_file_name:
        s = "verbose %s (%s) %s" % (level, caller_name, mesg)
    else:
        s = "verbose %s (%s: %s, %s) %s" % (level, lineno,
                                            short_file_name, caller_name, mesg)
    if t != 0:
        s = s + " (time = %s)" % cputime(t)
    print(s)
    sys.stdout.flush()
    return cputime()


def set_verbose(level, files='all'):
    """
    Set the global Sage verbosity level.

    INPUT:

    - ``level`` -- integer between 0 and 2, inclusive

    - ``files`` -- (default: ``'all'``) list of files to make verbose, or
      'all' to make ALL files verbose (the default)

    OUTPUT: changes the state of the verbosity flag and possibly
    appends to the list of files that are verbose

    EXAMPLES::

        sage: set_verbose(2)
        sage: verbose("This is Sage.", level=1)  # not tested
        VERBOSE1 (?): This is Sage.
        sage: verbose("This is Sage.", level=2)  # not tested
        VERBOSE2 (?): This is Sage.
        sage: verbose("This is Sage.", level=3)  # not tested
        [no output]
        sage: set_verbose(0)
    """
    if level is None:
        level = -1
    if isinstance(level, str):
        set_verbose_files([level])
    global LEVEL
    LEVEL = level
    if isinstance(files, str):
        files = [files]
    set_verbose_files(files)


def set_verbose_files(file_name):
    if not isinstance(file_name, list):
        file_name = [file_name]
    global verbose_files
    verbose_files = file_name


def get_verbose_files():
    return verbose_files


def unset_verbose_files(file_name):
    if not isinstance(file_name, list):
        file_name = [file_name]
    for X in file_name:
        verbose_files.remove(X)


def get_verbose():
    """
    Return the global Sage verbosity level.

    EXAMPLES::

        sage: get_verbose()
        0
        sage: set_verbose(2)
        sage: get_verbose()
        2
        sage: set_verbose(0)
    """
    global LEVEL
    return LEVEL
