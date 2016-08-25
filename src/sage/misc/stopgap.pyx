"""
Stopgaps
"""

########################################################################
#       Copyright (C) 2012 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

import warnings


from sage.doctest import DOCTEST_MODE
cdef bint ENABLED = not DOCTEST_MODE

def set_state(bint mode):
    """
    Enable or disable the stopgap warnings.

    INPUT:

    - ``mode`` -- (bool); if True, enable stopgaps; otherwise, disable.

    EXAMPLES::

        sage: import sage.misc.stopgap
        sage: sage.misc.stopgap.set_state(False)
        sage: sage.misc.stopgap.stopgap("Displays nothing.", 12345)
        sage: sage.misc.stopgap.set_state(True)
        sage: sage.misc.stopgap.stopgap("Displays something.", 123456)
        doctest:...:
        ********************************************************************************
        Displays something.
        This issue is being tracked at http://trac.sagemath.org/sage_trac/ticket/123456.
        ********************************************************************************
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

cdef set _stopgap_cache = set([])

def stopgap(message, int ticket_no):
    r"""
    Issue a stopgap warning.

    INPUT:

     - ``message`` - an explanation of how an incorrect answer might be produced.

     - ``ticket_no`` - an integer, giving the number of the Trac ticket tracking the underlying issue.

    EXAMPLES::

        sage: import sage.misc.stopgap
        sage: sage.misc.stopgap.set_state(True)
        sage: sage.misc.stopgap.stopgap("Computation of heights on elliptic curves over number fields can return very imprecise results.", 12509)
        doctest:...:
        ********************************************************************************
        Computation of heights on elliptic curves over number fields can return very imprecise results.
        This issue is being tracked at http://trac.sagemath.org/sage_trac/ticket/12509.
        ********************************************************************************
        sage: sage.misc.stopgap.stopgap("This is not printed", 12509)
        sage: sage.misc.stopgap.set_state(False)  # so rest of doctesting fine
    """
    if not ENABLED or ticket_no in _stopgap_cache:
        return
    # We reset show_warning so that the message is not printed twice.
    old_format = warnings.formatwarning
    def my_format(message, category, filename, lineno, line=None):
        return "%s:%s:\n%s\n%s\n%s\n" % (filename, lineno, "*"*80, message, "*"*80)
    warnings.formatwarning = my_format
    message = message + "\nThis issue is being tracked at http://trac.sagemath.org/sage_trac/ticket/%s."%ticket_no
    warnings.warn(StopgapWarning(message), stacklevel=2)
    warnings.formatwarning = old_format
    _stopgap_cache.add(ticket_no)
