# sage.doctest: needs sage.libs.pari
"""
Conjectural slopes of Hecke polynomials

Interface to Kevin Buzzard's PARI program for computing conjectural
slopes of characteristic polynomials of Hecke operators.

AUTHORS:

- William Stein (2006-03-05): Sage interface

- Kevin Buzzard: PARI program that implements underlying functionality
"""
#############################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
#############################################################################
from pathlib import Path

from sage.env import SAGE_EXTCODE
from sage.libs.pari import pari

buzzard_dir = Path(SAGE_EXTCODE) / "pari" / "buzzard"

# def buzzard_dimension_cusp_forms(eps, k):
#     r"""
#     eps is [N, i x 3 matrix], where eps[2][,1] is the primes dividing
#     N, eps[2][,2] is the powers of these primes that divide N, and eps[2][,3]
#     is the following: for p odd, p^n||N, it's t such that znprimroot(p^n)
#     gets sent to exp(2*pi*i/phi(p^n))^t. And for p=2, it's
#     0 for 2^1, it's 0 (trivial) or -1 (non-trivial) for 2^2, and for p^n>=8
#     it's either t>=0 for the even char sending 5 to exp(2*pi*i/p^(n-2))^t,
#     or t<=-1 for the odd char sending 5 to exp(2*pi*i/p^(n-2))^(-1-t).
#     (so either 0<=t<2^(n-2) or -1>=t>-1-2^(n-2) )

#     EXAMPLES::

#         sage: buzzard_dimension_cusp_forms('TrivialCharacter(100)', 4)

#     Next we compute a dimension for the character of level 45 which is
#     the product of the character of level 9 sending znprimroot(9)=2 to
#     $e^{2 \pi i/6}^1$ and the character of level 5 sending
#     \code{znprimroot(5)=2} to $e^{2 \pi i/4}^2=-1$.

#         sage: buzzard_dimension_cusp_forms('DirichletCharacter(45,[1,2])', 4)
#         <boom!>  which is why this is commented out!
#     """
#     s = gp().eval('DimensionCuspForms(%s, %s)'%(eps,k))
#     print(s)
#     return Integer(s)


def buzzard_tpslopes(p, N, kmax):
    r"""
    Return a vector of length kmax, whose `k`-th entry
    (`0 \leq k \leq k_{max}`) is the conjectural sequence
    of valuations of eigenvalues of `T_p` on forms of level
    `N`, weight `k`, and trivial character.

    This conjecture is due to Kevin Buzzard, and is only made assuming
    that `p` does not divide `N` and if `p` is
    `\Gamma_0(N)`-regular.

    EXAMPLES::

        sage: from sage.modular.buzzard import buzzard_tpslopes
        sage: c = buzzard_tpslopes(2,1,50)
        ...
        sage: c[50]
        [4, 8, 13]

    Hence Buzzard would conjecture that the `2`-adic valuations
    of the eigenvalues of `T_2` on cusp forms of level 1 and
    weight `50` are `[4,8,13]`, which indeed they are,
    as one can verify by an explicit computation using, e.g., modular
    symbols::

        sage: M = ModularSymbols(1,50, sign=1).cuspidal_submodule()
        sage: T = M.hecke_operator(2)
        sage: f = T.charpoly('x')
        sage: f.newton_slopes(2)
        [13, 8, 4]

    AUTHORS:

    - Kevin Buzzard: several PARI/GP scripts

    - William Stein (2006-03-17): small Sage wrapper of Buzzard's scripts
    """
    pari.read(buzzard_dir / "DimensionSk.g")
    pari.read(buzzard_dir / "genusn.g")
    pari.read(buzzard_dir / "Tpprog.g")
    # v = pari.tpslopes(p, N, kmax).sage()
    v = pari('tpslopes(%s, %s, %s)' % (p, N, kmax)).sage()
    v.insert(0, [])   # so that v[k] = info about weight k
    return v
