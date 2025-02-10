# sage.doctest: needs sage.libs.pari
r"""
Integer factorization using PARI

AUTHORS:

- Jeroen Demeyer (2015)
"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.pari.all import pari
from sage.rings.integer cimport Integer


def factor_using_pari(n, int_=False, debug_level=0, proof=None):
    r"""
    Factor this integer using PARI.

    This function returns a list of pairs, not a :class:`Factorization`
    object. The first element of each pair is the factor, of type
    ``Integer`` if ``int_`` is ``False`` or ``int`` otherwise,
    the second element is the positive exponent, of type ``int``.

    INPUT:

    - ``int_`` -- (default: ``False``), whether the factors are
      of type ``int`` instead of ``Integer``

    - ``debug_level`` -- (default: 0), debug level of the call
      to PARI

    - ``proof`` -- (default: ``None``), whether the factors are
      required to be proven prime;  if ``None``, the global default
      is used

    OUTPUT: list of pairs

    EXAMPLES::

        sage: factor(-2**72 + 3, algorithm='pari')  # indirect doctest
        -1 * 83 * 131 * 294971519 * 1472414939

    Check that PARI's debug level is properly reset (:issue:`18792`)::

        sage: from sage.doctest.util import ensure_interruptible_after
        sage: with ensure_interruptible_after(0.5): factor(2^1000 - 1, verbose=5)
        ...
        doctest:warning...
        RuntimeWarning: cypari2 leaked ... bytes on the PARI stack
        sage: pari.get_debug_level()
        0
    """
    if proof is None:
        from sage.structure.proof.proof import get_flag
        proof = get_flag(proof, "arithmetic")

    prev = pari.get_debug_level()

    cdef Py_ssize_t i
    try:
        if prev != debug_level:
            pari.set_debug_level(debug_level)

        p, e = n.__pari__().factor(proof=proof)
        if int_:
            return [(int(p[i]), int(e[i])) for i in range(len(p))]
        else:
            return [(Integer(p[i]), int(e[i])) for i in range(len(p))]
    finally:
        if prev != debug_level:
            pari.set_debug_level(prev)
