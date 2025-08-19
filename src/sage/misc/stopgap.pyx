"""
Stopgaps
"""

# ****************************************************************************
#       Copyright (C) 2012 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import warnings


from sage.doctest import DOCTEST_MODE
cdef bint ENABLED = not DOCTEST_MODE


def set_state(bint mode):
    """
    Enable or disable the stopgap warnings.

    INPUT:

    - ``mode`` -- boolean; if ``True``, enable stopgaps. Otherwise, disable.

    EXAMPLES::

        sage: import sage.misc.stopgap
        sage: sage.misc.stopgap.set_state(False)
        sage: sage.misc.stopgap.stopgap("Displays nothing.", 12345)
        sage: sage.misc.stopgap.set_state(True)
        sage: sage.misc.stopgap.stopgap("Displays something.", 123456)
        doctest:...:
        ...StopgapWarning: Displays something.
        This issue is being tracked at https://github.com/sagemath/sage/issues/123456.
        sage: sage.misc.stopgap.set_state(False)
    """
    global ENABLED
    ENABLED = mode


class StopgapWarning(Warning):
    """
    This class is used to warn users of a known issue that may produce
    mathematically incorrect results.
    """
    pass


warnings.filterwarnings('always', category=StopgapWarning)


cdef set _stopgap_cache = set()


def stopgap(message, int issue_no):
    r"""
    Issue a stopgap warning.

    INPUT:

    - ``message`` -- an explanation of how an incorrect answer might be produced

    - ``issue_no`` -- integer; giving the number of the Github issue tracking
      the underlying issue

    EXAMPLES::

        sage: import sage.misc.stopgap
        sage: sage.misc.stopgap.set_state(True)
        sage: sage.misc.stopgap.stopgap("Computation of heights on elliptic curves over number fields can return very imprecise results.", 12509)
        doctest:...
        ...StopgapWarning: Computation of heights on elliptic curves over number fields can return very imprecise results.
        This issue is being tracked at https://github.com/sagemath/sage/issues/12509.
        sage: sage.misc.stopgap.stopgap("This is not printed", 12509)
        sage: sage.misc.stopgap.set_state(False)  # so rest of doctesting fine
    """
    if not ENABLED or issue_no in _stopgap_cache:
        return
    # We reset show_warning so that the message is not printed twice.
    old_format = warnings.formatwarning

    def my_format(message, category, filename, lineno, line=None):
        return "%s:%s:\n%s\n%s\n%s\n" % (filename, lineno,
                                         "*" * 80, message, "*" * 80)
    warnings.formatwarning = my_format
    message = message + "\nThis issue is being tracked at https://github.com/sagemath/sage/issues/%s." % issue_no
    warnings.warn(StopgapWarning(message), stacklevel=2)
    warnings.formatwarning = old_format
    _stopgap_cache.add(issue_no)
