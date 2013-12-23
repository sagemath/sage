"""
qfsolve: Programme de resolution des equations quadratiques.

Interface to the GP quadratic forms code of Denis Simon.

AUTHORS:

 * Denis Simon (GP code)

 * Nick Alexander (Sage interface)
"""

#*****************************************************************************
#       Copyright (C) 2008 Nick Alexander
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.interfaces.gp import Gp
from sage.rings.all import ZZ, QQ

_gp_for_simon_interpreter = None    # Global GP interpreter for Denis Simon's code
def _gp_for_simon():
    r"""
    Start a GP interpreter for the use of Denis Simon's Qfsolve and Qfparam
    if it is not started already.

    EXAMPLE ::

        sage: from sage.quadratic_forms.qfsolve import _gp_for_simon
        sage: _gp_for_simon()
        PARI/GP interpreter
    """
    global _gp_for_simon_interpreter
    if _gp_for_simon_interpreter is None:
        _gp_for_simon_interpreter = Gp(script_subdirectory='simon')
        _gp_for_simon_interpreter.read("qfsolve.gp")
    return _gp_for_simon_interpreter

# \\ - Qfsolve(G,factD): pour resoudre l'equation quadratique X^t*G*X = 0
# \\ G doit etre une matrice symetrique n*n, a coefficients dans Z.
# \\ S'il n'existe pas de solution, la reponse est un entier
# \\ indiquant un corps local dans lequel aucune solution n'existe
# \\ (-1 pour les reels, p pour Q_p).
# \\ Si on connait la factorisation de -abs(2*matdet(G)),
# \\ on peut la passer par le parametre factD pour gagner du temps.
# \\
# \\ - Qfparam(G,sol,fl): pour parametrer les solutions de la forme
# \\ quadratique ternaire G, en utilisant la solution particuliere sol.
# \\ si fl>0, la 'fl'eme forme quadratique est reduite.

def qfsolve(G, factD=None):
    r"""
    Find a solution `x = (x_0,...,x_n)` to `x G x^t = 0` for an
    `n \times n`-matrix ``G`` over `\QQ`.

    If a solution exists, returns a tuple of rational numbers `x`.
    Otherwise, returns `-1` if no solutions exists over the reals or a
    prime `p` if no solution exists over the `p`-adic field `\QQ_p`.

    EXAMPLES::

        sage: from sage.quadratic_forms.qfsolve import qfsolve
        sage: M = Matrix(QQ, [[0, 0, -12], [0, -12, 0], [-12, 0, -1]]); M
        [  0   0 -12]
        [  0 -12   0]
        [-12   0  -1]
        sage: sol = qfsolve(M); sol
        (1, 0, 0)
        sage: sol[0].parent() is QQ
        True

        sage: M = Matrix(QQ, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        sage: ret = qfsolve(M); ret
        -1
        sage: ret.parent() is ZZ
        True

        sage: M = Matrix(QQ, [[1, 0, 0], [0, 1, 0], [0, 0, -7]])
        sage: qfsolve(M)
        7

        sage: M = Matrix(QQ, [[3, 0, 0, 0], [0, 5, 0, 0], [0, 0, -7, 0], [0, 0, 0, -11]])
        sage: qfsolve(M)
        (-3, 4, 3, 2)
    """
    gp = _gp_for_simon()
    if factD is not None:
        raise NotImplementedError, "qfsolve not implemented with parameter factD"
    ret = gp('Qfsolve(%s)' % G._pari_init_())
    if str(ret.type()) == 't_COL':  # Need explicit str(), see #15522
        return tuple(QQ(r) for r in ret)
    return ZZ(ret)

def qfparam(G, sol):
    r"""
    Parametrizes the conic defined by the matrix ``G``.

    INPUT:

     - ``G`` -- a `3 \times 3`-matrix over `\QQ`.

     - ``sol`` -- a triple of rational numbers providing a solution
       to sol*G*sol^t = 0.

    OUTPUT:

    A triple of polynomials that parametrizes all solutions of
    x*G*x^t = 0 up to scaling.

    ALGORITHM:

    Uses Denis Simon's pari script Qfparam.

    EXAMPLES::

        sage: from sage.quadratic_forms.qfsolve import qfsolve, qfparam
        sage: M = Matrix(QQ, [[0, 0, -12], [0, -12, 0], [-12, 0, -1]]); M
        [  0   0 -12]
        [  0 -12   0]
        [-12   0  -1]
        sage: sol = qfsolve(M);
        sage: ret = qfparam(M, sol); ret
        (-t^2 - 12, 24*t, 24*t^2)
        sage: ret[0].parent() is QQ['t']
        True
    """
    gp = _gp_for_simon()
    R = QQ['t']
    t = R.gen()
    s = 'Qfparam((%s), (%s)~)*[t^2,t,1]~' % (G._pari_init_(), list(sol))
    ret = gp(s)
    return tuple(R(str(r)) for r in ret)
