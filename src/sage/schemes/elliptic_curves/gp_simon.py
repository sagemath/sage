# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.pari
"""
Denis Simon's PARI scripts
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from pathlib import Path

from cypari2.handle_error import PariError

from sage.env import SAGE_EXTCODE
from sage.libs.pari import pari
from sage.misc.randstate import current_randstate
from sage.misc.superseded import deprecation
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.parent_gens import localvars


simon_dir = Path(SAGE_EXTCODE) / 'pari' / 'simon'


def simon_two_descent(E, verbose=0, lim1=None, lim3=None, limtriv=None,
                      maxprob=20, limbigprime=30, known_points=[]):
    """
    Interface to Simon's gp script for two-descent.

    .. NOTE::

        Users should instead run E.simon_two_descent()

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.gp_simon
        sage: E = EllipticCurve('389a1')
        sage: sage.schemes.elliptic_curves.gp_simon.simon_two_descent(E)
        doctest:warning...:
        DeprecationWarning: please use the 2-descent algorithm over QQ inside pari
        See https://github.com/sagemath/sage/issues/38461 for details.
        (2, 2, [(5/4 : 5/8 : 1), (-3/4 : 7/8 : 1)])

    TESTS::

        sage: # needs sage.rings.number_field
        sage: E = EllipticCurve('37a1').change_ring(QuadraticField(-11,'x'))
        sage: E.simon_two_descent()
        (1, 1, [(0 : 0 : 1)])

    An example with an elliptic curve defined over a relative number field::

        sage: # needs sage.rings.number_field
        sage: F.<a> = QuadraticField(29)
        sage: x = QQ['x'].gen()
        sage: K.<b> = F.extension(x^2-1/2*a+1/2)
        sage: E = EllipticCurve(K,[1, 0, 5/2*a + 27/2, 0, 0])   # long time (about 3 s)
        sage: E.simon_two_descent(lim1=2, limtriv=3)
        (1, 1, ...)

    Check that :issue:`16022` is fixed::

        sage: # needs sage.rings.number_field
        sage: K.<y> = NumberField(x^4 + x^2 - 7)
        sage: E = EllipticCurve(K, [1, 0, 5*y^2 + 16, 0, 0])
        sage: E.simon_two_descent(lim1=2, limtriv=3)            # long time (about 3 s)
        (1, 1, ...)

    An example that checks that :issue:`9322` is fixed (it should take less than a second to run)::

        sage: # needs sage.rings.number_field
        sage: K.<w> = NumberField(x^2 - x - 232)
        sage: E = EllipticCurve([2 - w, 18 + 3*w, 209 + 9*w, 2581 + 175*w, 852 - 55*w])
        sage: E.simon_two_descent()                             # long time
        (0, 2, [])
    """
    pari.read(simon_dir / "ellQ.gp")
    pari.read(simon_dir / "ell.gp")
    pari.read(simon_dir / "qfsolve.gp")
    pari.read(simon_dir / "resultant3.gp")

    current_randstate().set_seed_pari()

    K = E.base_ring()
    K_orig = K
    E_orig = E

    # The following is to correct the bug at #5204: the gp script
    # fails when K is a number field whose generator is called 'x'.
    # It also deals with relative number fields.

    if K is not QQ:
        K = K_orig.absolute_field('a')
        y = K.gen()
        from_K, to_K = K.structure()
        E = E_orig.change_ring(to_K)
        with localvars(K.polynomial().parent(), 'y'):
            # Simon's program requires that this name be y.
            K_pari = pari.bnfinit(K.polynomial())
        known_points = [P.change_ring(to_K) for P in known_points]
    else:
        deprecation(38461, "please use the 2-descent algorithm over QQ inside pari")
        from_K = lambda x: x

    # The block below mimics the defaults in Simon's scripts.
    # They need to be changed when these are updated.
    if K is QQ:
        over_QQ = True
        if lim1 is None:
            lim1 = 5
        if lim3 is None:
            lim3 = 50
        if limtriv is None:
            limtriv = 3
    else:
        over_QQ = False
        if lim1 is None:
            lim1 = 2
        if lim3 is None:
            lim3 = 4
        if limtriv is None:
            limtriv = 2

    pari('DEBUGLEVEL_ell=%s; LIM1=%s; LIM3=%s; LIMTRIV=%s; MAXPROB=%s; LIMBIGPRIME=%s;' % (
        verbose, lim1, lim3, limtriv, maxprob, limbigprime))

    try:
        if over_QQ:
            ans = pari("ellQ_ellrank")(E, known_points)
        else:
            ans = pari("bnfellrank")(K_pari, E, known_points)
    except PariError as err:
        raise RuntimeError("an error occurred while running Simon's 2-descent program") from err

    loc = {} if over_QQ else {'y': y}
    lower, upper, pts = ans.sage(locals=loc)
    lower = ZZ(lower)
    upper = ZZ(upper)
    points = [E_orig([from_K(c) for c in P]) for P in pts]
    points = [P for P in points if P.has_infinite_order()]
    return lower, upper, points
