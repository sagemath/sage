# -*- encoding: utf-8 -*-
"""
The Interact Decorator
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import inspect
import textwrap
from sage.misc.html import html


def interact(func, **kwds):
    """
    The interact decorator

    This is the global @interact-decorator. It should always return an
    interact implementation appropriate to the UI. Right now only
    SageNB is supported.

    INPUT:

    - ``func`` -- function. The interact function.

    - ``**kwds`` -- optional keyword arguments. To be passed to the
      UI interact decorator.

    EXAMPLES::

        sage: @interact
        ....: def f(x=[1,2,3]):
        ....:     print(x)
        <html>...</html>
    """
    from sagenb.notebook.interact import interact as sagenb_interact
    return sagenb_interact(func, **kwds)





