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

from sage.structure.parent_gens import localvars

from sage.interfaces.gp import Gp
from sage.misc.sage_eval import sage_eval
from sage.misc.randstate import current_randstate
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


gp = None
def init():
    """
    Function to initialize the gp process
    """
    global gp
    if gp is None:
        import os
        from sage.env import DOT_SAGE
        logfile = os.path.join(DOT_SAGE, 'gp-simon.log')
        gp = Gp(script_subdirectory='simon', logfile=logfile)
        gp.read("ellQ.gp")
        gp.read("ell.gp")
        gp.read("qfsolve.gp")
        gp.read("resultant3.gp")


def simon_two_descent(E, verbose=0, lim1=None, lim3=None, limtriv=None,
                      maxprob=20, limbigprime=30, known_points=[]):
    """
    Interface to Simon's gp script for two-descent.

    .. NOTE::

       Users should instead run E.simon_two_descent()

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.gp_simon
        sage: E=EllipticCurve('389a1')
        sage: sage.schemes.elliptic_curves.gp_simon.simon_two_descent(E)
        (2, 2, [(1 : 0 : 1), (-11/9 : 28/27 : 1)])

    TESTS::

        sage: E = EllipticCurve('37a1').change_ring(QuadraticField(-11,'x'))
        sage: E.simon_two_descent()
        (1, 1, [(0 : 0 : 1)])

    An example with an elliptic curve defined over a relative number field::

        sage: F.<a> = QuadraticField(29)
        sage: x = QQ['x'].gen()
        sage: K.<b> = F.extension(x^2-1/2*a+1/2)
        sage: E = EllipticCurve(K,[1, 0, 5/2*a + 27/2, 0, 0]) # long time (about 3 s)
        sage: E.simon_two_descent(lim1=2, limtriv=3)
        (1, 1, ...)

    Check that :trac:`16022` is fixed::

        sage: K.<y> = NumberField(x^4 + x^2 - 7)
        sage: E = EllipticCurve(K, [1, 0, 5*y^2 + 16, 0, 0])
        sage: E.simon_two_descent(lim1=2, limtriv=3)  # long time (about 3 s)
        (1, 1, ...)

    An example that checks that :trac:`9322` is fixed (it should take less than a second to run)::

        sage: K.<w> = NumberField(x^2-x-232)
        sage: E = EllipticCurve([2-w,18+3*w,209+9*w,2581+175*w,852-55*w])
        sage: E.simon_two_descent()
        (0, 2, [])
    """
    init()

    current_randstate().set_seed_gp(gp)

    K = E.base_ring()
    K_orig = K
    # The following is to correct the bug at #5204: the gp script
    # fails when K is a number field whose generator is called 'x'.
    # It also deals with relative number fields.
    E_orig = E
    if K is not QQ:
        K = K_orig.absolute_field('a')
        from_K, to_K = K.structure()
        E = E_orig.change_ring(to_K)
        known_points = [P.change_ring(to_K) for P in known_points]
        # Simon's program requires that this name be y.
        with localvars(K.polynomial().parent(), 'y'):
            gp.eval("K = bnfinit(%s);" % K.polynomial())
            if verbose >= 2:
                print("K = bnfinit(%s);" % K.polynomial())
        gp.eval("%s = Mod(y,K.pol);" % K.gen())
        if verbose >= 2:
            print("%s = Mod(y,K.pol);" % K.gen())
    else:
        from_K = lambda x: x
        to_K = lambda x: x

    # The block below mimics the defaults in Simon's scripts, and needs to be changed
    # when these are updated.
    if K is QQ:
        cmd = 'ellrank(%s, %s);' % (list(E.ainvs()), [P.__pari__() for P in known_points])
        if lim1 is None:
            lim1 = 5
        if lim3 is None:
            lim3 = 50
        if limtriv is None:
            limtriv = 3
    else:
        cmd = 'bnfellrank(K, %s, %s);' % (list(E.ainvs()), [P.__pari__() for P in known_points])
        if lim1 is None:
            lim1 = 2
        if lim3 is None:
            lim3 = 4
        if limtriv is None:
            limtriv = 2

    gp('DEBUGLEVEL_ell=%s; LIM1=%s; LIM3=%s; LIMTRIV=%s; MAXPROB=%s; LIMBIGPRIME=%s;'%(
       verbose, lim1, lim3, limtriv, maxprob, limbigprime))

    if verbose >= 2:
        print(cmd)
    s = gp.eval('ans=%s;'%cmd)
    if s.find(" *** ") != -1:
        raise RuntimeError("\n%s\nAn error occurred while running Simon's 2-descent program"%s)
    if verbose > 0:
        print(s)
    v = gp.eval('ans')
    if v=='ans': # then the call to ellrank() or bnfellrank() failed
        raise RuntimeError("An error occurred while running Simon's 2-descent program")
    if verbose >= 2:
        print("v = %s" % v)

    # pari represents field elements as Mod(poly, defining-poly)
    # so this function will return the respective elements of K
    def _gp_mod(*args):
        return args[0]
    ans = sage_eval(v, {'Mod': _gp_mod, 'y': K.gen(0)})
    lower = ZZ(ans[0])
    upper = ZZ(ans[1])
    points = [E_orig([from_K(c) for c in list(P)]) for P in ans[2]]
    points = [P for P in points if P.has_infinite_order()]
    return lower, upper, points
