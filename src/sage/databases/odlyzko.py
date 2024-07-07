"""
Database of the zeros of the Riemann zeta function

The main access function to the database of the zeros of the Riemann zeta
function is :func:`zeta_zeros`. In order to use ``zeta_zeros()``, you need to
install the optional :ref:`database_odlyzko_zeta <spkg_database_odlyzko_zeta>`
package::

    sage -i database_odlyzko_zeta

AUTHORS:

- William Stein: initial version
- Jeroen Demeyer (2015-01-20): converted ``database_odlyzko_zeta`` to new-style
  package
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os

from sage.misc.persist import load
from sage.env import SAGE_SHARE

def zeta_zeros():
    r"""
    List of the imaginary parts of the first 2,001,052 zeros of the
    Riemann zeta function, accurate to within 4e-9.

    REFERENCES:

    - http://www.dtc.umn.edu/~odlyzko/zeta_tables/index.html

    EXAMPLES:

    The following example shows the imaginary part of the 13th
    nontrivial zero of the Riemann zeta function::

        sage: # optional - database_odlyzko_zeta
        sage: zz = zeta_zeros()
        sage: zz[12]
        59.347044003
        sage: len(zz)
        2001052
    """
    from sage.misc.verbose import verbose
    sobj = os.path.join(SAGE_SHARE, 'odlyzko', 'zeros.sobj')
    verbose("Loading Odlyzko database from " + sobj)
    return load(sobj)
