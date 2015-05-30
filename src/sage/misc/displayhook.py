# -*- encoding: utf-8 -*-
"""
Entrypoint for the SageNB Display System

Just for compatibility with the notebook, you should not use this any
more. Look into ``sage/repl/`` instead.
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.repl.rich_output import get_display_manager
from sage.repl.rich_output.backend_sagenb import BackendSageNB


def DisplayHook():
    """
    This function is called by SageNB to set up its displayhook.

    OUTPUT:

    The new displayhook that will be used by SageNB

    EXAMPLES::

        sage: from sage.misc.displayhook import DisplayHook
        sage: d = DisplayHook()
        sage: d
        <bound method DisplayManager.displayhook of The 
        Sage display manager using the SageNB backend>

        sage: d(set([1, 2, 3]))       # Sage commandline output
        {1, 2, 3}

        sage: from sage.repl.rich_output import get_display_manager
        sage: get_display_manager()
        The Sage display manager using the SageNB backend
    """
    display_manager = get_display_manager()
    backend = BackendSageNB()
    display_manager.switch_backend(backend)
    return display_manager.displayhook

