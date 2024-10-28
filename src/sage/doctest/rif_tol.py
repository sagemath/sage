# sage_setup: distribution = sagemath-repl
"""
Helpers for tolerance checking in doctests
"""

# ****************************************************************************
#       Copyright (C) 2012-2018 David Roe <roed.math@gmail.com>
#                     2012      Robert Bradshaw <robertwb@gmail.com>
#                     2012      William Stein <wstein@gmail.com>
#                     2013      R. Andrew Ohana
#                     2013      Volker Braun
#                     2013-2018 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#                     2016-2021 Frédéric Chapoton
#                     2017-2018 Erik M. Bray
#                     2020      Marc Mezzarobba
#                     2020-2023 Matthias Koeppe
#                     2022      John H. Palmieri
#                     2022      Sébastien Labbé
#                     2023      Kwankyu Lee
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.doctest.marked_output import MarkedOutput


_RIFtol = None


def RIFtol(*args):
    """
    Create an element of the real interval field used for doctest tolerances.

    It allows large numbers like 1e1000, it parses strings with spaces
    like ``RIF(" - 1 ")`` out of the box and it carries a lot of
    precision. The latter is useful for testing libraries using
    arbitrary precision but not guaranteed rounding such as PARI. We use
    1044 bits of precision, which should be good to deal with tolerances
    on numbers computed with 1024 bits of precision.

    The interval approach also means that we do not need to worry about
    rounding errors and it is also very natural to see a number with
    tolerance as an interval.

    EXAMPLES::

        sage: from sage.doctest.parsing import RIFtol
        sage: RIFtol(-1, 1)
        0.?
        sage: RIFtol(" - 1 ")
        -1
        sage: RIFtol("1e1000")
        1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000?e1000
    """
    global _RIFtol
    if _RIFtol is None:
        try:
            # We need to import from sage.all to avoid circular imports.
            from sage.rings.real_mpfi import RealIntervalField
        except ImportError:
            from warnings import warn
            warn("RealIntervalField not available, ignoring all tolerance specifications in doctests")

            def fake_RIFtol(*args):
                return 0
            _RIFtol = fake_RIFtol
        else:
            _RIFtol = RealIntervalField(1044)
    return _RIFtol(*args)


def add_tolerance(wantval, want: MarkedOutput):
    """
    Enlarge the real interval element ``wantval`` according to
    the tolerance options in ``want``.

    INPUT:

    - ``wantval`` -- a real interval element
    - ``want`` -- a :class:`MarkedOutput` describing the tolerance

    OUTPUT: an interval element containing ``wantval``

    EXAMPLES::

        sage: from sage.doctest.parsing import MarkedOutput, SageOutputChecker
        sage: from sage.doctest.rif_tol import add_tolerance
        sage: want_tol = MarkedOutput().update(tol=0.0001)
        sage: want_abs = MarkedOutput().update(abs_tol=0.0001)
        sage: want_rel = MarkedOutput().update(rel_tol=0.0001)
        sage: add_tolerance(RIF(pi.n(64)), want_tol).endpoints()                 # needs sage.symbolic
        (3.14127849432443, 3.14190681285516)
        sage: add_tolerance(RIF(pi.n(64)), want_abs).endpoints()                 # needs sage.symbolic
        (3.14149265358979, 3.14169265358980)
        sage: add_tolerance(RIF(pi.n(64)), want_rel).endpoints()                 # needs sage.symbolic
        (3.14127849432443, 3.14190681285516)
        sage: add_tolerance(RIF(1e1000), want_tol)
        1.000?e1000
        sage: add_tolerance(RIF(1e1000), want_abs)
        1.000000000000000?e1000
        sage: add_tolerance(RIF(1e1000), want_rel)
        1.000?e1000
        sage: add_tolerance(0, want_tol)
        0.000?
        sage: add_tolerance(0, want_abs)
        0.000?
        sage: add_tolerance(0, want_rel)
        0
    """
    if want.tol:
        if wantval == 0:
            return RIFtol(want.tol) * RIFtol(-1, 1)
        else:
            return wantval * (1 + RIFtol(want.tol) * RIFtol(-1, 1))
    elif want.abs_tol:
        return wantval + RIFtol(want.abs_tol) * RIFtol(-1, 1)
    elif want.rel_tol:
        return wantval * (1 + RIFtol(want.rel_tol) * RIFtol(-1, 1))
    else:
        return wantval
