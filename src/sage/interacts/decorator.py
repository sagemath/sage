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


from sage.misc.superseded import deprecation
deprecation(22636, "the module sage.interacts.decorator is deprecated")


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
        sage: from sage.interacts.decorator import interact
        doctest:...: DeprecationWarning: the module sage.interacts.decorator is deprecated
        See http://trac.sagemath.org/22636 for details.
    """
    from sagenb.notebook.interact import interact as sagenb_interact
    return sagenb_interact(func, **kwds)
