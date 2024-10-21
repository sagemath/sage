# sage_setup: distribution = sagemath-symbolics
"""
Tests for the Sage/Pynac interaction
"""

# ****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


def rational_powers_memleak():
    """
    Check that there is no memory leak in rational powers.

    OUTPUT: boolean; whether the memory leak was detected

    See :issue:`9129`.

    EXAMPLES::

        sage: from sage.symbolic.tests import rational_powers_memleak
        sage: rational_powers_memleak()
        False
    """
    from sage.rings.integer_ring import ZZ
    import gc
    gc.collect()
    c0 = sum(1 for obj in gc.get_objects())
    for i in range(1000):
        _ = ZZ(2).sqrt()
    gc.collect()
    c1 = sum(1 for obj in gc.get_objects())
    # Test that we do not leak an object at each iteration
    return (c1 - c0) >= 1000
