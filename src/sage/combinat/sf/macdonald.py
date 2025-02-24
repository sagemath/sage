# sage.doctest: needs sage.combinat sage.modules
r"""
Macdonald Polynomials

Notation used in the definitions follows mainly [Mac1995]_.

The integral forms of the bases `H` and `Ht` do not appear in
Macdonald's book.  They correspond to the two bases
`H_\mu[X;q,t] = \sum_{\nu} K_{\nu\mu}(q,t) s_\mu[X]` and
`{\tilde H}_\mu[X;q,t] = t^{n(\mu)} \sum_{\nu} K_{\nu\mu}(q,1/t) s_\nu[X]`
where `K_{\mu\nu}(q,t)` are the Macdonald `q,t`-Koskta coefficients.

The `Ht` in this case is short for `{\tilde H}` and is the basis which is
the graded Frobenius image of the Garsia-Haiman modules [GH1993]_.

REFERENCES:

- [Mac1995]_

.. [GH1993] \A. Garsia, M. Haiman, A graded representation module for Macdonald's
   polynomials, Proc. Nat. Acad. U.S.A. no. 90, 3607--3610.

.. [BGHT1999] \F. Bergeron, A. M. Garsia, M. Haiman, and G. Tesler, Identities and
   positivity conjectures for some remarkable operators in the theory of symmetric
   functions, Methods Appl. Anal. 6 (1999), no. 3, 363--420.

.. [LLM1998] \L. Lapointe, A. Lascoux, J. Morse, Determinantal Expressions for
   Macdonald Polynomials, IRMN no. 18 (1998).
   :arxiv:`math/9808050`.

.. [BH2013] \F. Bergeron, M. Haiman, Tableaux Formulas for Macdonald Polynomials,
   Special edition in honor of Christophe Reutenauer 60 birthday, International
   Journal of Algebra and Computation, Volume 23, Issue 4, (2013), pp. 833-852.
"""

#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.categories.modules_with_basis import ModulesWithBasis
from . import sfa
from sage.combinat.partition import Partitions_n, _Partitions
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.rational_field import QQ
from sage.misc.misc_c import prod
from sage.misc.cachefunc import cached_function
import functools

# cache in q,t globally and subs locally with q and t values
# these caches are stored in self._self_to_s_cache and self._s_to_self_cache
#J basis cache
_j_to_s_cache = {}
_s_to_j_cache = {}

#Ht basis cache
_ht_to_m_cache = {}

#S basis cache
_S_to_s_cache = {}
_s_to_S_cache = {}

_qt_kostka_cache = {}


class Macdonald(UniqueRepresentation):

    def __repr__(self):
        r"""
        The family of Macdonald symmetric function bases.

        INPUT:

        - ``self`` -- a family of Macdonald symmetric function bases

        OUTPUT: string representing the Macdonald symmetric function family

        EXAMPLES::

            sage: t = QQ['t'].gen(); SymmetricFunctions(QQ['t'].fraction_field()).macdonald(q=t,t=1)
            Macdonald polynomials with q=t and t=1 over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        return self._name

    @staticmethod
    def __classcall__(cls, Sym, q='q', t='t'):
        """
        Normalize the arguments.

        TESTS::

            sage: R.<q, t> = QQ[]
            sage: B1 = SymmetricFunctions(R).macdonald().P()
            sage: B2 = SymmetricFunctions(R).macdonald(q, t).P()
            sage: B3 = SymmetricFunctions(R).macdonald(t, q).P()
            sage: B1 is B2
            True
            sage: B1 == B3
            False
        """
        return super().__classcall__(cls, Sym, Sym.base_ring()(q), Sym.base_ring()(t))

    def __init__(self, Sym, q, t):
        r"""
        Macdonald Symmetric functions including `P`, `Q`, `J`, `H`, `Ht` bases
        also including the S basis which is the plethystic transformation
        of the Schur basis (that which is dual to the Schur basis
        with respect to the Macdonald `q,t`-scalar product)

        INPUT:

        - ``self`` -- a family of Macdonald symmetric function bases

        EXAMPLES::

            sage: t = QQ['t'].gen(); SymmetricFunctions(QQ['t'].fraction_field()).macdonald(q=t,t=1)
            Macdonald polynomials with q=t and t=1 over Fraction Field of Univariate Polynomial Ring in t over Rational Field
            sage: Sym = SymmetricFunctions(FractionField(QQ['t'])).macdonald()
            Traceback (most recent call last):
            ...
            TypeError: unable to evaluate 'q' in Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        self._sym = Sym
        self._s = Sym.s()
        self.q = q
        self.t = t
        self._name_suffix = ""
        if str(q) != 'q':
            self._name_suffix += " with q=%s" % q
            if str(t) != 't':
                self._name_suffix += " and "
        if str(t) != 't':
            if str(q) == 'q':
                self._name_suffix += " with "
            self._name_suffix += "t=%s" % t
        self._name = "Macdonald polynomials"+self._name_suffix+" over "+repr(Sym.base_ring())

    def base_ring( self ):
        r"""
        Return the base ring of the symmetric functions where the
        Macdonald symmetric functions live.

        INPUT:

        - ``self`` -- a family of Macdonald symmetric function bases

        OUTPUT: the base ring associated to the corresponding symmetric function ring

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ['q'].fraction_field())
            sage: Mac = Sym.macdonald(t=0)
            sage: Mac.base_ring()
            Fraction Field of Univariate Polynomial Ring in q over Rational Field
        """
        return self._sym.base_ring()

    def symmetric_function_ring( self ):
        r"""
        Return the base ring of the symmetric functions where the
        Macdonald symmetric functions live.

        INPUT:

        - ``self`` -- a family of Macdonald symmetric function bases

        OUTPUT: the symmetric function ring associated to the Macdonald bases

        EXAMPLES::

            sage: Mac = SymmetricFunctions(QQ['q'].fraction_field()).macdonald(t=0)
            sage: Mac.symmetric_function_ring()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in q over Rational Field
        """
        return self._sym

    def P(self):
        r"""
        Return Macdonald polynomials in `P` basis.
        The `P` basis is defined here as a normalized form of the `J` basis.

        INPUT:

        - ``self`` -- a family of Macdonald symmetric function bases

        OUTPUT: the `P` Macdonald basis of symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: P = Sym.macdonald().P(); P
            Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald P basis
            sage: P[2]
            McdP[2]

        The `P` Macdonald basis is upper triangularly related to the monomial symmetric functions and are
        orthogonal with respect to the `qt`-Hall scalar product::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: P = Sym.macdonald().P(); P
            Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald P basis
            sage: m = Sym.monomial()
            sage: P.transition_matrix(m,2)
            [                          1 (q*t - q + t - 1)/(q*t - 1)]
            [                          0                           1]
            sage: P([1,1]).scalar_qt(P([2]))
            0
            sage: P([2]).scalar_qt(P([2]))
            (-q^3 + q^2 + q - 1)/(-q*t^2 + q*t + t - 1)
            sage: P([1,1]).scalar_qt(P([1,1]))
            (-q^2*t + q*t + q - 1)/(-t^3 + t^2 + t - 1)

        When `q = 0`, the Macdonald polynomials on the `P` basis are the same
        as the Hall-Littlewood polynomials on the `P` basis.

        ::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: P = Sym.macdonald(q=0).P(); P
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Macdonald P with q=0 basis
            sage: P([2])^2
            (t+1)*McdP[2, 2] + (-t+1)*McdP[3, 1] + McdP[4]
            sage: HLP = Sym.hall_littlewood().P()
            sage: HLP([2])^2
            (t+1)*HLP[2, 2] + (-t+1)*HLP[3, 1] + HLP[4]

        Coercions from the `Q` and `J` basis (proportional) are
        implemented::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: P = Sym.macdonald().P()
            sage: Q = Sym.macdonald().Q()
            sage: J = Sym.macdonald().J()
            sage: s = Sym.schur()

        ::

            sage: P(Q([2]))
            ((q*t^2-q*t-t+1)/(q^3-q^2-q+1))*McdP[2]
            sage: P(Q([2,1]))
            ((-q*t^4+2*q*t^3-q*t^2+t^2-2*t+1)/(-q^4*t+2*q^3*t-q^2*t+q^2-2*q+1))*McdP[2, 1]

        ::

            sage: P(J([2]))
            (q*t^2-q*t-t+1)*McdP[2]
            sage: P(J([2,1]))
            (-q*t^4+2*q*t^3-q*t^2+t^2-2*t+1)*McdP[2, 1]

        By transitivity, one get coercions from the classical bases::

            sage: P(s([2]))
            ((q-t)/(q*t-1))*McdP[1, 1] + McdP[2]
            sage: P(s([2,1]))
            ((q*t-t^2+q-t)/(q*t^2-1))*McdP[1, 1, 1] + McdP[2, 1]

        ::

            sage: Sym = SymmetricFunctions(QQ['x','y','z'].fraction_field())
            sage: (x,y,z) = Sym.base_ring().gens()
            sage: Macxy = Sym.macdonald(q=x,t=y)
            sage: Macyz = Sym.macdonald(q=y,t=z)
            sage: Maczx = Sym.macdonald(q=z,t=x)
            sage: P1 = Macxy.P()
            sage: P2 = Macyz.P()
            sage: P3 = Maczx.P()
            sage: m(P1[2,1])
            ((-2*x*y^2+x*y-y^2+x-y+2)/(-x*y^2+1))*m[1, 1, 1] + m[2, 1]
            sage: m(P2[2,1])
            ((-2*y*z^2+y*z-z^2+y-z+2)/(-y*z^2+1))*m[1, 1, 1] + m[2, 1]
            sage: m(P1(P2(P3[2,1])))
            ((-2*x^2*z-x^2+x*z-x+z+2)/(-x^2*z+1))*m[1, 1, 1] + m[2, 1]
            sage: P1(P2[2])
            ((-x*y^2+2*x*y*z-y^2*z-x+2*y-z)/(x*y^2*z-x*y-y*z+1))*McdP[1, 1] + McdP[2]
            sage: m(z*P1[2]+x*P2[2])
            ((x^2*y^2*z+x*y^2*z^2-x^2*y^2+x^2*y*z-x*y*z^2+y^2*z^2-x^2*y-2*x*y*z-y*z^2+x*y-y*z+x+z)/(x*y^2*z-x*y-y*z+1))*m[1, 1] + (x+z)*m[2]
        """
        return MacdonaldPolynomials_p(self)

    def Q(self):
        r"""
        Return the Macdonald polynomials on the `Q` basis. These are dual to
        the Macdonald polynomials on the P basis with respect to the
        `qt`-Hall scalar product.
        The `Q` basis is defined to be a normalized form of the `J` basis.

        INPUT:

        - ``self`` -- a family of Macdonald symmetric function bases

        OUTPUT: the `Q` Macdonald basis of symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: Q = Sym.macdonald().Q(); Q
            Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald Q basis
            sage: P = Sym.macdonald().P()
            sage: Q([2]).scalar_qt(P([2]))
            1
            sage: Q([2]).scalar_qt(P([1,1]))
            0
            sage: Q([1,1]).scalar_qt(P([2]))
            0
            sage: Q([1,1]).scalar_qt(P([1,1]))
            1
            sage: Q(P([2]))
            ((q^3-q^2-q+1)/(q*t^2-q*t-t+1))*McdQ[2]
            sage: Q(P([1,1]))
            ((q^2*t-q*t-q+1)/(t^3-t^2-t+1))*McdQ[1, 1]


        Coercions from the `P` and `J` basis (proportional) are implemented::

            sage: P = Sym.macdonald().P()
            sage: Q = Sym.macdonald().Q()
            sage: J = Sym.macdonald().J()
            sage: s = Sym.schur()

        ::

            sage: Q(J([2]))
            (q^3-q^2-q+1)*McdQ[2]

        ::

            sage: Q(P([2]))
            ((q^3-q^2-q+1)/(q*t^2-q*t-t+1))*McdQ[2]
            sage: P(Q(P([2])))
            McdP[2]
            sage: Q(P(Q([2])))
            McdQ[2]

        By transitivity, one gets coercions from the classical bases::

            sage: Q(s([2]))
            ((q^2-q*t-q+t)/(t^3-t^2-t+1))*McdQ[1, 1] + ((q^3-q^2-q+1)/(q*t^2-q*t-t+1))*McdQ[2]
        """
        return MacdonaldPolynomials_q(self)

    def J(self):
        r"""
        Return the Macdonald polynomials on the `J` basis also known as the
        integral form of the Macdonald polynomials. These are scalar
        multiples of both the `P` and `Q` bases. When expressed in the `P` or `Q`
        basis, the scaling coefficients are polynomials in `q` and `t` rather
        than rational functions.

        The `J` basis is calculated using determinantal formulas of
        Lapointe-Lascoux-Morse giving the action on the S-basis [LLM1998]_.

        INPUT:

        - ``self`` -- a family of Macdonald symmetric function bases

        OUTPUT: the `J` Macdonald basis of symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: J = Sym.macdonald().J(); J
            Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald J basis
            sage: P = Sym.macdonald().P()
            sage: Q = Sym.macdonald().Q()
            sage: P(J([2]))
            (q*t^2-q*t-t+1)*McdP[2]
            sage: P(J([1,1]))
            (t^3-t^2-t+1)*McdP[1, 1]
            sage: Q(J([2]))
            (q^3-q^2-q+1)*McdQ[2]
            sage: Q(J([1,1]))
            (q^2*t-q*t-q+1)*McdQ[1, 1]

        Coercions from the `Q` and `J` basis (proportional) and to/from
        the Schur basis are implemented::

            sage: P = Sym.macdonald().P()
            sage: Q = Sym.macdonald().Q()
            sage: J = Sym.macdonald().J()
            sage: s = Sym.schur()

        ::

            sage: J(P([2]))
            (1/(q*t^2-q*t-t+1))*McdJ[2]

        ::

            sage: J(Q([2]))
            (1/(q^3-q^2-q+1))*McdJ[2]

        ::

            sage: s(J([2]))
            (-q*t+t^2+q-t)*s[1, 1] + (q*t^2-q*t-t+1)*s[2]
            sage: J(s([2]))
            ((q-t)/(q*t^4-q*t^3-q*t^2-t^3+q*t+t^2+t-1))*McdJ[1, 1] + (1/(q*t^2-q*t-t+1))*McdJ[2]
        """
        return MacdonaldPolynomials_j(self)

    def H(self):
        r"""
        Return the Macdonald polynomials on the H basis. When the `H` basis
        is expanded on the Schur basis, the coefficients are the `qt`-Kostka
        numbers.

        INPUT:

        - ``self`` -- a family of Macdonald symmetric function bases

        OUTPUT: the `H` Macdonald basis of symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: H = Sym.macdonald().H(); H
            Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald H basis
            sage: s = Sym.schur()
            sage: s(H([2]))
            q*s[1, 1] + s[2]
            sage: s(H([1,1]))
            s[1, 1] + t*s[2]

        Coercions to/from the Schur basis are implemented::

            sage: H = Sym.macdonald().H()
            sage: s = Sym.schur()
            sage: H(s([2]))
            (q/(q*t-1))*McdH[1, 1] - (1/(q*t-1))*McdH[2]
        """
        return MacdonaldPolynomials_h(self)

    def Ht(self):
        r"""
        Return the Macdonald polynomials on the `Ht` basis. The elements of
        the `Ht` basis are eigenvectors of the `nabla` operator. When expanded
        on the Schur basis, the coefficients are the modified `qt`-Kostka
        numbers.

        INPUT:

        - ``self`` -- a family of Macdonald symmetric function bases

        OUTPUT: the `Ht` Macdonald basis of symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: Ht = Sym.macdonald().Ht(); Ht
            Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald Ht basis
            sage: [Ht(p).nabla() for p in Partitions(3)]
            [q^3*McdHt[3], q*t*McdHt[2, 1], t^3*McdHt[1, 1, 1]]

        ::

            sage: s = Sym.schur()
            sage: from sage.combinat.sf.macdonald import qt_kostka
            sage: q,t = Ht.base_ring().gens()
            sage: s(Ht([2,1]))
            q*t*s[1, 1, 1] + (q+t)*s[2, 1] + s[3]
            sage: qt_kostka([1,1,1],[2,1]).subs(t=1/t)*t^Partition([2,1]).weighted_size()
            q*t
            sage: qt_kostka([2,1],[2,1]).subs(t=1/t)*t^Partition([2,1]).weighted_size()
            q + t
            sage: qt_kostka([3],[2,1]).subs(t=1/t)*t^Partition([2,1]).weighted_size()
            1

        Coercions to/from the Schur basis are implemented::

            sage: Ht = Sym.macdonald().Ht()
            sage: s = Sym.schur()
            sage: Ht(s([2,1]))
            (q/(q*t^2-t^3-q^2+q*t))*McdHt[1, 1, 1] + ((-q^2-q*t-t^2)/(q^2*t^2-q^3-t^3+q*t))*McdHt[2, 1] + (t/(-q^3+q^2*t+q*t-t^2))*McdHt[3]
            sage: Ht(s([2]))
            ((-q)/(-q+t))*McdHt[1, 1] + (t/(-q+t))*McdHt[2]
        """
        return MacdonaldPolynomials_ht(self)

    def S(self):
        r"""
        Return the modified Schur functions defined by the plethystic
        substitution `S_{\mu} = s_{\mu}[X(1-t)/(1-q)]`. When the
        Macdonald polynomials in the J basis are expressed in terms of the
        modified Schur functions at `q=0`, the coefficients are `qt`-Kostka numbers.

        INPUT:

        - ``self`` -- a family of Macdonald symmetric function bases

        OUTPUT: the `S` Macdonald basis of symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: S = Sym.macdonald().S(); S
            Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the Macdonald S basis
            sage: p = Sym.power()
            sage: p(S[2,1])
            ((1/3*t^3-t^2+t-1/3)/(q^3-3*q^2+3*q-1))*p[1, 1, 1] + ((-1/3*t^3+1/3)/(q^3-1))*p[3]
            sage: J = Sym.macdonald().J()
            sage: S(J([2]))
            (q^3-q^2-q+1)*McdS[2]
            sage: S(J([1,1]))
            (q^2*t-q*t-q+1)*McdS[1, 1] + (q^2-q*t-q+t)*McdS[2]
            sage: S = Sym.macdonald(q=0).S()
            sage: S(J[1,1])
            McdS[1, 1] + t*McdS[2]
            sage: S(J[2])
            q*McdS[1, 1] + McdS[2]
            sage: p(S[2,1])
            (-1/3*t^3+t^2-t+1/3)*p[1, 1, 1] + (1/3*t^3-1/3)*p[3]

            sage: from sage.combinat.sf.macdonald import qt_kostka
            sage: qt_kostka([2],[1,1])
            t
            sage: qt_kostka([1,1],[2])
            q

        Coercions to/from the Schur basis are implemented::

            sage: S = Sym.macdonald().S()
            sage: s = Sym.schur()
            sage: S(s([2]))
            ((q^2-q*t-q+t)/(t^3-t^2-t+1))*McdS[1, 1] + ((-q^2*t+q*t+q-1)/(-t^3+t^2+t-1))*McdS[2]
            sage: s(S([1,1]))
            ((-q*t^2+q*t+t-1)/(-q^3+q^2+q-1))*s[1, 1] + ((q*t-t^2-q+t)/(-q^3+q^2+q-1))*s[2]
        """
        return MacdonaldPolynomials_s(self)


QQqt = QQ['q,t'].fraction_field()

##############################################


def c1(part, q, t):
    r"""
    This function returns the qt-Hall scalar product between ``J(part)``
    and ``P(part)``.

    This coefficient is `c_\lambda` in equation (8.1') p. 352 of
    Macdonald's book [Mac1995]_.

    INPUT:

    - ``part`` -- a partition
    - ``q``, ``t`` -- parameters

    OUTPUT: a polynomial of the scalar product between the `J` and `P` bases

    EXAMPLES::

        sage: from sage.combinat.sf.macdonald import c1
        sage: R.<q,t> = QQ[]
        sage: c1(Partition([2,1]),q,t)
        -q^4*t + 2*q^3*t - q^2*t + q^2 - 2*q + 1
        sage: c1(Partition([1,1]),q,t)
        q^2*t - q*t - q + 1
    """
    R = q.parent()
    arms = part.arm_lengths(flat=True)
    legs = part.leg_lengths(flat=True)
    return R.prod(1 - q**(a + 1) * t**l for a, l in zip(arms, legs))


def c2(part, q, t):
    r"""
    This function returns the qt-Hall scalar product between J(part)
    and Q(part).

    This coefficient is `c_\lambda` in equation (8.1) p. 352 of
    Macdonald's book [Mac1995]_.

    INPUT:

    - ``part`` -- a partition
    - ``q``, ``t`` -- parameters

    OUTPUT: a polynomial of the scalar product between the `J` and `P` bases

    EXAMPLES::

        sage: from sage.combinat.sf.macdonald import c2
        sage: R.<q,t> = QQ[]
        sage: c2(Partition([1,1]),q,t)
        t^3 - t^2 - t + 1
        sage: c2(Partition([2,1]),q,t)
        -q*t^4 + 2*q*t^3 - q*t^2 + t^2 - 2*t + 1
    """
    R = q.parent()
    arms = part.arm_lengths(flat=True)
    legs = part.leg_lengths(flat=True)
    return R.prod(1 - q**a * t**(l + 1) for a, l in zip(arms, legs))


@cached_function
def cmunu1(mu, nu):
    r"""
    Return the coefficient of `{\tilde H}_\nu` in `h_1^\perp {\tilde H}_\mu`.

    INPUT:

    - ``mu``, ``nu`` -- partitions with ``nu`` precedes ``mu``

    OUTPUT: an element of the fraction field of polynomials in `q` and `t`

    EXAMPLES::

        sage: from sage.combinat.sf.macdonald import cmunu1
        sage: cmunu1(Partition([2,1]),Partition([2]))
        (-t^2 + q)/(q - t)
        sage: cmunu1(Partition([2,1]),Partition([1,1]))
        (-q^2 + t)/(-q + t)
        sage: Sym = SymmetricFunctions(QQ['q','t'].fraction_field())
        sage: h = Sym.h()
        sage: Ht = Sym.macdonald().Ht()
        sage: all(Ht[3,2,1].skew_by(h[1]).coefficient(nu)
        ....:     == cmunu1(Partition([3,2,1]),nu)
        ....:     for nu in Partition([3,2,1]).down_list())
        True
    """
    q,t = QQqt.gens()
    # The following for loop is equivalent to getting the cell:
    #    SkewPartition([mu,nu]).cells()[0]
    for i, val in enumerate(nu._list):
        if val < mu._list[i]:
            A = prod((t**mu.leg_length(i, s) - q**(mu.arm_length(i, s)+1))
                     / (t**nu.leg_length(i, s) - q**(nu.arm_length(i, s)+1))
                     for s in range(val))
            B = prod((q**mu.arm_length(*s) - t**(mu.leg_length(*s)+1))
                     / (q**nu.arm_length(*s) - t**(nu.leg_length(*s)+1))
                     for s in nu.cells() if s[1] == val)
            return QQqt(A * B)

    return QQqt(prod( (q**mu.arm_length(s, 0) - t**(mu.leg_length(s, 0)+1))
                      / (q**nu.arm_length(s, 0) - t**(nu.leg_length(s, 0)+1))
                      for s in range(len(nu._list)) ))


@cached_function
def cmunu(mu, nu):
    r"""
    Return the coefficient of `{\tilde H}_\nu` in `h_r^\perp {\tilde H}_\mu`.

    Proposition 5 of F. Bergeron and M. Haiman [BH2013]_ states

    .. MATH::

        c_{\mu\nu} = \sum_{\alpha \leftarrow \nu} c_{\mu\alpha}
        c_{\alpha\nu} B_{\alpha/\nu}/B_{\mu/\nu}

    where `c_{\mu\nu}` is the coefficient of `{\tilde H}_\nu` in
    `h_r^\perp {\tilde H}_\mu` and `B_{\mu/\nu}` is the bi-exponent generator
    implemented in the function :func:`sage.combinat.sf.macdonald.Bmu`.

    INPUT:

    - ``mu``, ``nu`` -- partitions with ``nu`` contained in ``mu``

    OUTPUT: an element of the fraction field of polynomials in `q` and `t`

    EXAMPLES::

        sage: from sage.combinat.sf.macdonald import cmunu
        sage: cmunu(Partition([2,1]),Partition([1]))
        q + t + 1
        sage: cmunu(Partition([2,2]),Partition([1,1]))
        (-q^3 - q^2 + q*t + t)/(-q + t)
        sage: Sym = SymmetricFunctions(QQ['q','t'].fraction_field())
        sage: h = Sym.h()
        sage: Ht = Sym.macdonald().Ht()
        sage: all(Ht[2,2].skew_by(h[r]).coefficient(nu)
        ....:     == cmunu(Partition([2,2]),nu)
        ....:     for r in range(1,5) for nu in Partitions(4-r))
        True
    """
    if not nu:
        return QQqt.one()
    if not mu.contains(nu):
        return QQqt.zero()
    if mu.size() == nu.size() + 1:
        return cmunu1(mu, nu)

    # This is equivalent to:
    #   Bmu(SkewPartition([outer, inner]))
    def Bmu_skew(outer, inner):
        inner = list(inner) # This makes a (shallow) copy of inner
        inner += [0]*(len(outer)-len(inner))
        q,t = QQqt.gens()
        res = QQqt.zero()
        for i, val in enumerate(outer):
            for j in range(inner[i], val):
                res += t**i * q**j
        return res

    nulist = nu._list
    return (sum(cmunu(mu, al) * cmunu1(al, nu) * Bmu_skew(al, nulist)
                for al in nu.up()) / Bmu_skew(mu, nulist))

#Generic MacdonaldPolynomials


class MacdonaldPolynomials_generic(sfa.SymmetricFunctionAlgebra_generic):

    def __init__(self, macdonald):
        r"""
        A class for methods for one of the Macdonald bases of the symmetric functions.

        INPUT:

        - ``self`` -- a Macdonald basis
        - ``macdonald`` -- a family of Macdonald symmetric function bases

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q,t'])); Sym.rename('Sym'); Sym
            Sym
            sage: Sym.macdonald().P()
            Sym in the Macdonald P basis
            sage: Sym.macdonald(t=2).P()
            Sym in the Macdonald P with t=2 basis
            sage: Sym.rename()

        TESTS::

            sage: Sym.macdonald().P()._prefix
            'McdP'
            sage: Sym.macdonald().Ht()._prefix
            'McdHt'
        """
        s = self.__class__.__name__[21:].capitalize()
        sfa.SymmetricFunctionAlgebra_generic.__init__(
            self, macdonald._sym,
            basis_name="Macdonald " + s + macdonald._name_suffix,
            prefix="Mcd" + s)
        self.q = macdonald.q
        self.t = macdonald.t
        self._macdonald = macdonald
        self._s = self._macdonald._s

        # Bases defined by orthotriangularity should inherit from some
        # common category BasesByOrthotriangularity (shared with Jack, HL, orthotriang, Mcdo)
        if hasattr(self, "_s_cache"):
            # temporary until Hom(GradedHopfAlgebrasWithBasis work better)
            category = ModulesWithBasis(self.base_ring())
            self.register_coercion(SetMorphism(Hom(self._s, self, category), self._s_to_self))
            self._s.register_coercion(SetMorphism(Hom(self, self._s, category), self._self_to_s))

    def construction(self):
        """
        Return a pair ``(F, R)``, where ``F`` is a
        :class:`SymmetricFunctionsFunctor` and `R` is a ring, such
        that ``F(R)`` returns ``self``.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q']))
            sage: J = Sym.macdonald(t=2).J()
            sage: J.construction()
            (SymmetricFunctionsFunctor[Macdonald J with t=2],
             Fraction Field of Univariate Polynomial Ring in q over Rational Field)
        """
        return (sfa.SymmetricFunctionsFamilyFunctor(self, Macdonald,
                                                    self.basis_name(),
                                                    self.q, self.t),
                self.base_ring())

    def _s_to_self(self, x):
        r"""
        Isomorphism from the Schur basis into ``self``.

        INPUT:

        - ``self`` -- a Macdonald basis
        - ``x`` -- an element of the Schur basis

        OUTPUT: the basis element ``x`` in the basis ``self``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q']))
            sage: J = Sym.macdonald(t=2).J()
            sage: s = Sym.schur()
            sage: J._s_to_self(s[2,1])
            ((-1/28*q+1/14)/(q-1/4))*McdJ[1, 1, 1] - (1/4/(q-1/4))*McdJ[2, 1]

        This is for internal use only. Please use instead::

            sage: J(s[2,1])
            ((-1/28*q+1/14)/(q-1/4))*McdJ[1, 1, 1] - (1/4/(q-1/4))*McdJ[2, 1]
        """
        return self._from_cache(x, self._s_cache, self._s_to_self_cache,
                                q=self.q, t=self.t)

    def _self_to_s(self, x):
        r"""
        Isomorphism from ``self`` to the Schur basis.

        INPUT:

        - ``self`` -- a Macdonald basis
        - ``x`` -- an element of a Macdonald basis

        OUTPUT: the basis element ``x`` in the Schur functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q']))
            sage: J = Sym.macdonald(t=2).J()
            sage: s = Sym.schur()
            sage: J._self_to_s(J[2,1])
            (3*q-6)*s[1, 1, 1] + (-4*q+1)*s[2, 1]

        This is for internal use only. Please use instead::

            sage: s(J[2,1])
            (3*q-6)*s[1, 1, 1] + (-4*q+1)*s[2, 1]
        """
        return self._s._from_cache(x, self._s_cache, self._self_to_s_cache,
                                   q=self.q, t=self.t)

    def c1(self, part):
        r"""
        Return the `qt`-Hall scalar product between ``J(part)`` and ``P(part)``.

        INPUT:

        - ``self`` -- a Macdonald basis
        - ``part`` -- a partition

        OUTPUT:

        - returns the `qt`-Hall scalar product between ``J(part)`` and ``P(part)``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: P = Sym.macdonald().P()
            sage: P.c1(Partition([2,1]))
            -q^4*t + 2*q^3*t - q^2*t + q^2 - 2*q + 1
        """
        return c1(part, self.q, self.t)

    def c2(self, part):
        r"""
        Return the `qt`-Hall scalar product between ``J(part)`` and ``Q(part)``.

        INPUT:

        - ``self`` -- a Macdonald basis
        - ``part`` -- a partition

        OUTPUT:

        - returns the `qt`-Hall scalar product between ``J(part)`` and ``Q(part)``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: P = Sym.macdonald().P()
            sage: P.c2(Partition([2,1]))
            -q*t^4 + 2*q*t^3 - q*t^2 + t^2 - 2*t + 1
        """
        return c2(part, self.q, self.t)

    def product(self, left, right):
        r"""
        Multiply an element of the Macdonald symmetric function
        basis ``self`` and another symmetric function

        Convert to the Schur basis, do the multiplication there, and
        convert back to ``self`` basis.

        INPUT:

        - ``self`` -- a Macdonald symmetric function basis
        - ``left`` -- an element of the basis ``self``
        - ``right`` -- another symmetric function

        OUTPUT: the product of ``left`` and ``right`` expanded in the basis ``self``

        EXAMPLES::

            sage: Mac = SymmetricFunctions(FractionField(QQ['q','t'])).macdonald()
            sage: H = Mac.H()
            sage: J = Mac.J()
            sage: P = Mac.P()
            sage: Q = Mac.Q()
            sage: Ht = Mac.Ht()
            sage: J([1])^2 #indirect doctest
            ((q-1)/(q*t-1))*McdJ[1, 1] + ((t-1)/(q*t-1))*McdJ[2]
            sage: J.product( J[1], J[2] )
            ((-q^2+1)/(-q^2*t+1))*McdJ[2, 1] + ((-t+1)/(-q^2*t+1))*McdJ[3]
            sage: H.product( H[1], H[2] )
            ((q^2-1)/(q^2*t-1))*McdH[2, 1] + ((-t+1)/(-q^2*t+1))*McdH[3]
            sage: P.product( P[1], P[2] )
            ((-q^3*t^2+q*t^2+q^2-1)/(-q^3*t^2+q^2*t+q*t-1))*McdP[2, 1] + McdP[3]
            sage: Q.product(Q[1],Q[2])
            McdQ[2, 1] + ((q^2*t-q^2+q*t-q+t-1)/(q^2*t-1))*McdQ[3]
            sage: Ht.product(Ht[1],Ht[2])
            ((q^2-1)/(q^2-t))*McdHt[2, 1] + ((t-1)/(-q^2+t))*McdHt[3]
        """
        return self(self._s(left) * self._s(right))

    def macdonald_family(self):
        r"""
        Return the family of Macdonald bases associated to the basis ``self``.

        INPUT:

        - ``self`` -- a Macdonald basis

        OUTPUT: the family of Macdonald symmetric functions associated to ``self``

        EXAMPLES::

            sage: MacP = SymmetricFunctions(QQ['q'].fraction_field()).macdonald(t=0).P()
            sage: MacP.macdonald_family()
            Macdonald polynomials with t=0 over Fraction Field of Univariate Polynomial Ring in q over Rational Field
        """
        return self._macdonald

    class Element(sfa.SymmetricFunctionAlgebra_generic.Element):

        def nabla(self, q=None, t=None, power=1):
            r"""
            Return the value of the nabla operator applied to ``self``.

            The eigenvectors of the nabla operator are the Macdonald
            polynomials in the `Ht` basis.  For more information
            see: [BGHT1999]_.

            The operator nabla acts on symmetric functions and has the
            Macdonald `Ht` basis as eigenfunctions and the eigenvalues
            are `q^{n(\mu')} t^{n(\mu)}` where
            `n(\mu) = \sum_{i} (i-1) \mu_i` and `\mu'` is the conjugate
            shape of `\mu`.

            If the parameter ``power`` is an integer then it calculates
            nabla to that integer.  The default value of ``power`` is 1.

            INPUT:

            - ``self`` -- an element of a Macdonald basis
            - ``q``, ``t`` -- (optional) parameters to specialize
            - ``power`` -- integer (default: 1)

            OUTPUT: the symmetric function of `\nabla` acting on ``self``

            EXAMPLES::

                sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
                sage: P = Sym.macdonald().P()
                sage: P([1,1]).nabla()
                ((q^2*t+q*t^2-2*t)/(q*t-1))*McdP[1, 1] + McdP[2]
                sage: P([1,1]).nabla(t=1)
                ((q^2*t+q*t-t-1)/(q*t-1))*McdP[1, 1] + McdP[2]
                sage: H = Sym.macdonald().H()
                sage: H([1,1]).nabla()
                t*McdH[1, 1] + (-t^2+1)*McdH[2]
                sage: H([1,1]).nabla(q=1)
                ((t^2+q-t-1)/(q*t-1))*McdH[1, 1] + ((-t^3+t^2+t-1)/(q*t-1))*McdH[2]
                sage: H(0).nabla()
                0
                sage: H([2,2,1]).nabla(t=1/H.t)
                ((-q^2)/(-t^4))*McdH[2, 2, 1]
                sage: H([2,2,1]).nabla(t=1/H.t,power=-1)
                ((-t^4)/(-q^2))*McdH[2, 2, 1]
            """
            parent = self.parent()
            if (q is None and t is None):
                Ht = parent._macdonald.Ht()
            else:
                if q is None:
                    q = parent.q
                if t is None:
                    t = parent.t
                Ht = parent.realization_of().macdonald(q=q,t=t).Ht()
            return parent(Ht(self).nabla(power=power))

#P basis


class MacdonaldPolynomials_p(MacdonaldPolynomials_generic):
    def __init__(self, macdonald):
        r"""
        The `P` basis is defined here as the `J` basis times a
        normalizing coefficient `c2`.

        INPUT:

        - ``self`` -- a Macdonald `P` basis
        - ``macdonald`` -- a family of Macdonald bases

        TESTS::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: P = Sym.macdonald().P()
            sage: TestSuite(P).run(skip=["_test_associativity","_test_distributivity","_test_prod"])  # long time (20s on sage.math, 2012)
            sage: TestSuite(P).run(elements = [P.t*P[1,1]+P.q*P[2], P[1]+(P.q+P.t)*P[1,1]])  # long time (depends on previous)
        """
        MacdonaldPolynomials_generic.__init__(self, macdonald)

        self._J = macdonald.J()
        # temporary until Hom(GradedHopfAlgebrasWithBasis work better)
        category = ModulesWithBasis(self.base_ring())
        phi = self._J.module_morphism(diagonal=self.c2,
                                      codomain=self, category=category)
        self.register_coercion(phi)
        self._J.register_coercion(~phi)

    def scalar_qt_basis(self, part1, part2=None):
        r"""
        Return the scalar product of `P(part1)` and `P(part2)`
        This scalar product formula is given in equation (4.11) p.323
        and (6.19) p.339 of Macdonald's book [Mac1995]_.

        INPUT:

        - ``self`` -- a Macdonald `P` basis
        - ``part1``, ``part2`` -- partitions

        OUTPUT:

        - returns the scalar product of ``P(part1)`` and ``P(part2)``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: P = Sym.macdonald().P()
            sage: P.scalar_qt_basis(Partition([2,1]), Partition([1,1,1]))
            0
            sage: f = P.scalar_qt_basis(Partition([3,2,1]), Partition([3,2,1]))
            sage: factor(f.numerator())
            (q - 1)^3 * (q^2*t - 1)^2 * (q^3*t^2 - 1)
            sage: factor(f.denominator())
            (t - 1)^3 * (q*t^2 - 1)^2 * (q^2*t^3 - 1)

        With a single argument, takes `part2 = part1`::

            sage: P.scalar_qt_basis(Partition([2,1]), Partition([2,1]))
            (-q^4*t + 2*q^3*t - q^2*t + q^2 - 2*q + 1)/(-q*t^4 + 2*q*t^3 - q*t^2 + t^2 - 2*t + 1)
        """
        if part2 is not None and part1 != part2:
            return self.base_ring().zero()
        return self.c1(part1) / self.c2(part1)

    class Element(MacdonaldPolynomials_generic.Element):
        pass


#Q basis
class MacdonaldPolynomials_q(MacdonaldPolynomials_generic):
    def __init__(self, macdonald):
        r"""
        The `Q` basis is defined here as the `J` basis times a
        normalizing coefficient.

        INPUT:

        - ``self`` -- a Macdonald `Q` basis
        - ``macdonald`` -- a family of Macdonald bases

        TESTS::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: Q = Sym.macdonald().Q()
            sage: TestSuite(Q).run(skip=["_test_associativity","_test_distributivity","_test_prod"])  # long time (29s on sage.math, 2012)
            sage: TestSuite(Q).run(elements = [Q.t*Q[1,1]+Q.q*Q[2], Q[1]+(Q.q+Q.t)*Q[1,1]])  # long time (depends on previous)
        """
        MacdonaldPolynomials_generic.__init__(self, macdonald)

        self._J = macdonald.J()
        self._P = macdonald.P()

        # temporary until Hom(GradedHopfAlgebrasWithBasis) works better
        category = ModulesWithBasis(self.base_ring())
        phi = self._P.module_morphism(diagonal=self._P.scalar_qt_basis,
                                      codomain=self, category=category)
        self.register_coercion(phi)
        self._P.register_coercion(~phi)

    class Element(MacdonaldPolynomials_generic.Element):
        pass


class MacdonaldPolynomials_j(MacdonaldPolynomials_generic):
    def __init__(self, macdonald):
        r"""
        The `J` basis is calculated using determinantal formulas of
        Lapointe-Lascoux-Morse giving the action on the `S`-basis.

        INPUT:

        - ``self`` -- a Macdonald `J` basis
        - ``macdonald`` -- a family of Macdonald bases

        TESTS::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: J = Sym.macdonald().J()
            sage: TestSuite(J).run(skip=["_test_associativity","_test_distributivity","_test_prod"])  # long time (19s on sage.math, 2012)
            sage: TestSuite(J).run(elements = [J.t*J[1,1]+J.q*J[2], J[1]+(J.q+J.t)*J[1,1]])  # long time (depends on previous)
        """
        self._self_to_s_cache = _j_to_s_cache
        self._s_to_self_cache = _s_to_j_cache
        MacdonaldPolynomials_generic.__init__(self, macdonald)

    def _s_cache(self, n):
        r"""
        Compute the change of basis and its inverse between the Macdonald
        polynomials on the `J` basis and the Schur functions
        these computations are completed with coefficients in fraction
        field of polynomials in `q` and `t`

        INPUT:

        - ``self`` -- a Macdonald `J` basis
        - ``n`` -- nonnegative integer

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: J = Sym.macdonald().J()
            sage: J._s_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l( J._s_to_self_cache[2] )
            [([1, 1], [([1, 1], 1/(t^3 - t^2 - t + 1))]),
             ([2],
              [([1, 1], (q - t)/(q*t^4 - q*t^3 - q*t^2 - t^3 + q*t + t^2 + t - 1)),
               ([2], 1/(q*t^2 - q*t - t + 1))])]
            sage: l( J._self_to_s_cache[2] )
            [([1, 1], [([1, 1], t^3 - t^2 - t + 1)]),
             ([2], [([1, 1], -q*t + t^2 + q - t), ([2], q*t^2 - q*t - t + 1)])]
        """
        self._invert_morphism(n, QQqt, self._self_to_s_cache,
                              self._s_to_self_cache,
                              to_other_function=self._to_s,
                              upper_triangular=False)

    def _to_s(self, part):
        r"""
        Return a function which gives the coefficient of a partition in
        the Schur expansion of self(part).

        These computations are completed with coefficients in fraction
        field of polynomials in `q` and `t`

        INPUT:

        - ``self`` -- a Macdonald `J` basis
        - ``part`` -- a partition

        OUTPUT:

        - returns a function which returns the coefficients of the expansion of `J`
          in the Schur basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: J = Sym.macdonald().J()
            sage: f21 = J._to_s(Partition([2,1]))
            sage: [f21(part) for part in Partitions(3)]
            [0,
             -q*t^4 + 2*q*t^3 - q*t^2 + t^2 - 2*t + 1,
             q*t^3 - t^4 - q*t^2 + t^3 - q*t + t^2 + q - t]
             sage: Sym.schur()( J[2,1] )
             (q*t^3-t^4-q*t^2+t^3-q*t+t^2+q-t)*s[1, 1, 1] + (-q*t^4+2*q*t^3-q*t^2+t^2-2*t+1)*s[2, 1]
        """
        q, t = QQqt.gens()
        S = self._macdonald.S()
        res = S(1)
        for k in reversed(part):
            res = res.creation(k)
        res = res._omega_qt_in_schurs()
        res = res.map_coefficients(lambda c: c(t, q))
        f = res.coefficient
        return f

    class Element(MacdonaldPolynomials_generic.Element):
        pass


class MacdonaldPolynomials_h(MacdonaldPolynomials_generic):
    def __init__(self, macdonald):
        r"""
        The `H` basis is defined as `H_\mu = \sum_{\lambda} K_{\lambda\mu}(q,t) s_\lambda`
        where `K_{\lambda\mu}(q,t)` are the Macdonald Kostka coefficients.

        In this implementation, it is calculated by using the Macdonald `Ht` basis and
        substituting `t \rightarrow 1/t` and multiplying by `t^{n(\mu)}`.

        INPUT:

        - ``self`` -- a Macdonald `H` basis
        - ``macdonald`` -- a family of Macdonald bases

        TESTS::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: H = Sym.macdonald().H()
            sage: TestSuite(H).run(skip=["_test_associativity","_test_distributivity","_test_prod"])
            sage: TestSuite(H).run(elements = [H.t*H[1,1]+H.q*H[2], H[1]+(H.q+H.t)*H[1,1]])  # long time (26s on sage.math, 2012)
        """
        MacdonaldPolynomials_generic.__init__(self, macdonald)
        self._m = self._sym.m()
        self._Lmunu = macdonald.Ht()._Lmunu
        if not self.t:
            self._Qp = self._sym.hall_littlewood(t=self.q).Qp()
        category = ModulesWithBasis(self.base_ring())
        self._s.register_coercion(SetMorphism(Hom(self, self._s, category), self._self_to_s))
        self.register_coercion(SetMorphism(Hom(self._s, self, category), self._s_to_self))
        self._m.register_coercion(SetMorphism(Hom(self, self._m, category), self._self_to_m))
        self.register_coercion(SetMorphism(Hom(self._m, self, category), self._m_to_self))

    def _self_to_s(self, x):
        r"""
        Return an element ``x`` of ``self`` expanded in the Schur basis.

        If `t=0` then the expression for `H_\mu[X;q,0]` is equal to
        `\omega Q'_{\mu'}[X;q]` where `Q'` is the Hall-Littlewood basis.

        If `t \neq 0` then the conversion uses the monomial basis.

        INPUT:

        - ``x`` -- an element of ``H`` basis

        OUTPUT: an element of the Schur basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: H = Sym.macdonald().H()
            sage: s = Sym.s()
            sage: s(H[2,1])
            q*s[1, 1, 1] + (q*t+1)*s[2, 1] + t*s[3]
            sage: H2 = Sym.macdonald(t=0).H()
            sage: s(H2[2,1])
            q*s[1, 1, 1] + s[2, 1]

            sage: Sym = SymmetricFunctions(FractionField(QQ['x']))
            sage: x = Sym.base_ring().gen()
            sage: H = Sym.macdonald(q=x,t=1/x).H()
            sage: s = Sym.s()
            sage: s(H[2,1])
            x*s[1, 1, 1] + 2*s[2, 1] + 1/x*s[3]
        """
        if self.t:
            return self._s(self._self_to_m(x))
        else:
            return sum(cmu*self._s(self._Qp(mu.conjugate())) for mu,cmu in x).omega()

    def _s_to_self(self, x):
        r"""
        Return an element of the Schur basis ``x`` expanded in the ``H`` basis.

        If `t=0` then the expression for `H_\mu[X;q,0]` is equal to
        `\omega Q'_{\mu'}[X;q]` where `Q'` is the Hall-Littlewood basis.

        If `t \neq 0` then the conversion uses the monomial basis.

        INPUT:

        - ``x`` -- an element of the Schur basis

        OUTPUT: an element of the ``H`` basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: H = Sym.macdonald().H()
            sage: s = Sym.s()
            sage: H(s[1,1])
            -(1/(q*t-1))*McdH[1, 1] + (t/(q*t-1))*McdH[2]
            sage: (q,t) = Sym.base_ring().gens()
            sage: H(q*s[1, 1, 1] + (q*t+1)*s[2, 1] + t*s[3])
            McdH[2, 1]
            sage: H2 = Sym.macdonald(t=0).H()
            sage: H2(q*s[1, 1, 1] + (q*t+1)*s[2, 1] + t*s[3])
            (-q^2*t+1)*McdH[2, 1] + t*McdH[3]

            sage: Sym = SymmetricFunctions(FractionField(QQ['x']))
            sage: x = Sym.base_ring().gen()
            sage: H = Sym.macdonald(q=x,t=1/x).H()
            sage: s = Sym.s()
            sage: H(x*s[1, 1, 1] + 2*s[2, 1] + 1/x*s[3])
            McdH[2, 1]
        """
        if self.t:
            return self._m_to_self(self._m(x))
        else:
            return self._from_dict({mu.conjugate() : cmu for mu,cmu in self._Qp(x.omega())})

    def _self_to_m(self, x):
        r"""
        Return an element ``x`` of ``self`` expanded in the monomial basis.

        INPUT:

        - ``x`` -- an element of ``H`` basis

        OUTPUT: an element of the monomial basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: H = Sym.macdonald().H()
            sage: m = Sym.m()
            sage: m(H[2,1])
            (2*q*t+q+t+2)*m[1, 1, 1] + (q*t+t+1)*m[2, 1] + t*m[3]

            sage: Sym = SymmetricFunctions(FractionField(QQ['x']))
            sage: x = Sym.base_ring().gen()
            sage: H = Sym.macdonald(q=x,t=1/x).H()
            sage: m = Sym.m()
            sage: m(H[2,1])
            ((x^2+4*x+1)/x)*m[1, 1, 1] + ((2*x+1)/x)*m[2, 1] + 1/x*m[3]
            sage: H2 = Sym.macdonald(q=x,t=1).H()
            sage: H2((3*x+3)*m[1, 1, 1] + (x+2)*m[2, 1] + m[3])
            McdH[2, 1]
            sage: H3 = Sym.macdonald(q=x,t=0).H()
            sage: H3((3*x+3)*m[1, 1, 1] + (x+2)*m[2, 1] + m[3])
            (-x^2+1)*McdH[2, 1] + McdH[3]
        """
        if self.t:
            tinv = ~self.t
            part_coeff = lambda x, d: sorted((mu,c) for mu,c in x if sum(mu) == d)
            return self._m._from_dict({ part2:
                self._base( sum(c * self.t**mu.weighted_size()
                                * self._Lmunu(part2, mu).subs(q=self.q, t=tinv)
                                for mu,c in part_coeff(x, d)) )
                for d in range(x.degree()+1) for part2 in Partitions_n(d) })
        else:
            return self._m(self._self_to_s(x))

    def _m_to_self( self, f ):
        r"""
        Convert an element ``f`` from the monomial basis to the ``H`` basis.

        This calculation is performed by using the fact that `H_\mu[X(1-t)]`
        is `c_\mu m_\mu` plus terms which are smaller in dominance order.
        The leading coefficient of the expansion of ``f`` in the `H_\mu`
        basis is equal to the leading coefficient of `c_\mu^{-1} f[X(1-t)]`.

        If `t=1`, we must appeal to another triangularity since ``Ht``
        is (usually) still a basis, however `H_\mu[X(1-t)]=0`.  In this
        case we assume that it is a basis and `H_\mu[X(q-1)]` is
        `c_\mu m_{\mu'}` plus terms which are lower in dominance order.

        INPUT:

        - ``f`` -- an element of the monomial basis

        OUTPUT: an element of the ``H`` basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: H = Sym.macdonald().H()
            sage: m = Sym.m()
            sage: H(m[1,1])
            -(1/(q*t-1))*McdH[1, 1] + (t/(q*t-1))*McdH[2]
            sage: (q,t) = Sym.base_ring().gens()
            sage: H((2*q*t+q+t+2)*m[1, 1, 1] + (q*t+t+1)*m[2, 1] + t*m[3])
            McdH[2, 1]

            sage: Sym = SymmetricFunctions(FractionField(QQ['x']))
            sage: x = Sym.base_ring().gen()
            sage: H = Sym.macdonald(q=x,t=1).H()
            sage: m = Sym.m()
            sage: H((3*x+3)*m[1, 1, 1] + (x+2)*m[2, 1] + m[3])
            McdH[2, 1]
        """
        if self.t == 1:
            g = f.omega_qt(q=self.q, t=0)
            fl = lambda x: x.conjugate()
            mu_to_H = lambda mu: self._self_to_m(self(mu.conjugate())).omega_qt(q=self.q, t=0)
        else:
            g = f.theta_qt(q=self.t, t=0)
            fl = lambda x: x
            mu_to_H = lambda mu: self._self_to_m(self(mu)).theta_qt(q=self.t, t=0)
        out = {}
        while not g.is_zero():
            sprt = sorted(g.support())
            Hmu = mu_to_H(sprt[-1])
            fl_sprt = fl(sprt[-1])
            out[fl_sprt] = self._base(g.coefficient(sprt[-1]) / Hmu.coefficient(sprt[-1]))
            g -= out[fl_sprt] * Hmu
        return self._from_dict(out)

    class Element(MacdonaldPolynomials_generic.Element):
        pass


class MacdonaldPolynomials_ht(MacdonaldPolynomials_generic):
    def __init__(self, macdonald):
        r"""
        The `Ht` basis is defined as `{\tilde H}_\mu = t^{n(\mu)} \sum_{\lambda}
        K_{\lambda\mu}(q,t^{-1}) s_\lambda` where `K_{\lambda\mu}(q,t)` are the
        Macdonald `(q,t)`-Kostka coefficients and `n(\mu) = \sum_{i} (i-1) \mu_i`.

        It is implemented here by using a Pieri formula due to F. Bergeron
        and M. Haiman [BH2013]_.

        INPUT:

        - ``self`` -- a Macdonald `Ht` basis
        - ``macdonald`` -- a family of Macdonald bases

        TESTS::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: Ht = Sym.macdonald().Ht()
            sage: TestSuite(Ht).run(skip=["_test_associativity","_test_distributivity","_test_prod"])  # long time (26s on sage.math, 2012)
            sage: TestSuite(Ht).run(elements = [Ht.t*Ht[1,1]+Ht.q*Ht[2], Ht[1]+(Ht.q+Ht.t)*Ht[1,1]])  # long time (depends on previous)
        """
        MacdonaldPolynomials_generic.__init__(self, macdonald)
        self._self_to_m_cache = _ht_to_m_cache
        self._m = self._sym.m()
        category = ModulesWithBasis(self.base_ring())
        self._s.register_coercion(SetMorphism(Hom(self, self._s, category), self._self_to_s))
        self.register_coercion(SetMorphism(Hom(self._s, self, category), self._s_to_self))
        self._m.register_coercion(SetMorphism(Hom(self, self._m, category), self._self_to_m))
        self.register_coercion(SetMorphism(Hom(self._m, self, category), self._m_to_self))

    def _self_to_s(self, x):
        r"""
        Convert an element of the ``Ht`` basis to the Schur basis.

        This function is here to force the coercion path to the Schur basis
        because these bases are computed using their monomial expansion.

        INPUT:

        - ``x`` -- an element of ``self``

        OUTPUT: an element of the Schur basis

        EXAMPLES::

            sage: Ht = SymmetricFunctions(FractionField(QQ['q','t'])).macdonald().Ht()
            sage: s = Ht.symmetric_function_ring().s()
            sage: Ht._self_to_s(Ht[2,1])
            q*t*s[1, 1, 1] + (q+t)*s[2, 1] + s[3]
        """
        return self._s(self._self_to_m(x))

    def _s_to_self( self, x ):
        r"""
        Convert an element of either the Schur basis to the ``Ht`` basis.

        This function is here to force the coercion path from the Schur basis
        because these bases are computed using the monomial expansion.

        INPUT:

        - ``x`` -- an element of ``s`` basis

        OUTPUT: an element of the basis ``self``

        EXAMPLES::

            sage: s = SymmetricFunctions(FractionField(QQ['q','t'])).s()
            sage: Ht = s.symmetric_function_ring().macdonald().Ht()
            sage: Ht._s_to_self(s[2])
            ((-q)/(-q+t))*McdHt[1, 1] + (t/(-q+t))*McdHt[2]
        """
        return self._m_to_self(self._m(x))

    def _Lmunu(self, nu, mu):
        r"""
        Return the coefficient of `m_\nu` in `{\tilde H}_\mu`.

        The coefficient is a `(q,t)`-analogue of `n` choose `\nu`.  Let
        `c_{\mu\gamma}` be the coefficient of `{\tilde H}_\gamma` in
        `h_r^\perp {\tilde H}_\mu`.  The coefficient of `m_\nu` in
        `{\tilde H}_\mu` is

        .. MATH::

            L_{\mu\nu} = \sum_{\gamma} c_{\mu\gamma} L_{\gamma\overline{\nu}},

        where the sum is over partitions `\gamma \subseteq \mu` with
        `\gamma \vdash |(\nu_1, \ldots, \nu_{\ell(\nu)-1})|`.  There is a
        recursive formula for `c_{\mu\gamma}` which is due to F. Bergeron
        and M. Haiman [BH2013]_ and this is implemented in
        :func:`sage.combinat.sf.macdonald.cmunu`.

        INPUT:

        - ``nu``, ``mu`` -- partitions of the same size

        OUTPUT: a polynomial in `q` and `t`

        EXAMPLES::

            sage: Lmunu = SymmetricFunctions(FractionField(QQ['q','t'])).macdonald().Ht()._Lmunu
            sage: Lmunu(Partition([3,1]),Partition([3,1]))
            q^2 + q + t + 1
            sage: Lmunu(Partition([3,1]),Partition([2,2]))
            q*t + q + t + 1
            sage: Lmunu(Partition([3,1]),Partition([2,1,1]))
            t^2 + q + t + 1
            sage: Lmunu(Partition([2,2]),Partition([2,1,1]))
            q*t + 2*t^2 + q + t + 1
        """
        if not mu:
            if not nu:
                return QQqt.one()
            else:
                return QQqt.zero()
        if (mu,nu) in self._self_to_m_cache:
            return self._self_to_m_cache[(mu,nu)]
        if len(nu) == 1:
            return QQqt.one()
        short_nu = _Partitions(nu[:-1])
        if nu[-1] == 1:
            self._self_to_m_cache[(mu,nu)] = QQqt( sum(cmunu1(mu,ga) * self._Lmunu(short_nu, ga)
                                                       for ga in mu.down()) )
        else:
            self._self_to_m_cache[(mu,nu)] = QQqt( sum(cmunu(mu,ga) * self._Lmunu(short_nu, ga)
                        for ga in Partitions_n(short_nu.size()) if mu.contains(ga) ) )
        return self._self_to_m_cache[(mu,nu)]

    def _self_to_m(self, x):
        r"""
        Take an element of the ``Ht`` basis and return the expansion in the
        monomial basis.

        INPUT:

        - ``x`` -- an element of ``Ht`` basis

        OUTPUT: an element of the monomial basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: Ht = Sym.macdonald().Ht()
            sage: m = Sym.m()
            sage: m(Ht[2,1])
            (q*t+2*q+2*t+1)*m[1, 1, 1] + (q+t+1)*m[2, 1] + m[3]

            sage: Sym = SymmetricFunctions(FractionField(QQ['x']))
            sage: x = Sym.base_ring().gen()
            sage: Ht = Sym.macdonald(q=x,t=1/x).Ht()
            sage: m = Sym.m()
            sage: m(Ht[2,1])
            ((2*x^2+2*x+2)/x)*m[1, 1, 1] + ((x^2+x+1)/x)*m[2, 1] + m[3]
        """
        part_coeff = lambda x, d: sorted((mu,c) for mu,c in x if sum(mu) == d)
        return self._m._from_dict({ part2:
            self._base( sum(c * self._Lmunu(part2, mu).subs(q=self.q, t=self.t)
                            for mu,c in part_coeff(x, d)) )
                    for d in range(x.degree()+1) for part2 in Partitions_n(d) })

    def _m_to_self( self, f ):
        r"""
        Convert an element ``f`` from the monomial basis to the ``Ht`` basis.

        This calculation is performed by using the fact that
        `{\tilde H}_\mu[X(t-1)]` is `c_\mu m_\mu` plus terms which are
        smaller in dominance order.  The leading coefficient of the
        expansion of ``f`` in the `{\tilde H}_\mu` basis is equal to
        the leading coefficient of `c_\mu^{-1} f[X(t-1)]`.

        If `t=1`, we must appeal to another triangularity since ``Ht``
        is (usually) still a basis, however `{\tilde H}_\mu[X(t-1)]=0`.
        In this case we assume that it is a basis and `{\tilde H}_\mu[X(q-1)]`
        is `c_\mu m_{\mu'}` plus terms which are lower in dominance order.

        INPUT:

        - ``f`` -- an element of the monomial basis

        OUTPUT: an element of the ``Ht`` basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: Ht = Sym.macdonald().Ht()
            sage: m = Sym.m()
            sage: Ht(m[1,1])
            (1/(-q+t))*McdHt[1, 1] - (1/(-q+t))*McdHt[2]
            sage: (q,t) = Sym.base_ring().gens()
            sage: Ht((q*t+2*q+2*t+1)*m[1, 1, 1] + (q+t+1)*m[2, 1] + m[3])
            McdHt[2, 1]

            sage: Sym = SymmetricFunctions(FractionField(QQ['x']))
            sage: x = Sym.base_ring().gen()
            sage: Ht = Sym.macdonald(q=x,t=1).Ht()
            sage: m = Sym.m()
            sage: Ht((3*x+3)*m[1, 1, 1] + (x+2)*m[2, 1] + m[3])
            McdHt[2, 1]
        """
        if self.t == 1:
            subsval = self.q
            fl = lambda x: x.conjugate()
        else:
            subsval = self.t
            fl = lambda x: x
        g = f.omega_qt(q=subsval, t=0)
        out = {}
        while not g.is_zero():
            sprt = sorted(g.support())
            Htmu = self._self_to_m(self(fl(sprt[-1]))).omega_qt(q=subsval, t=0)
            out[fl(sprt[-1])] = self._base(g.coefficient(sprt[-1]) / Htmu.coefficient(sprt[-1]))
            g -= out[fl(sprt[-1])] * Htmu
        return self._from_dict(out)

    class Element(MacdonaldPolynomials_generic.Element):
        def nabla(self, q=None, t=None, power=1):
            r"""
            Return the value of the nabla operator applied to ``self``. The
            eigenvectors of the `nabla` operator are the Macdonald polynomials in
            the `Ht` basis.  For more information see: [BGHT1999]_.

            The operator `nabla` acts on symmetric functions and has the
            Macdonald `Ht` basis as eigenfunctions and the eigenvalues
            are `q^{n(\mu')} t^{n(\mu)}` where `n(\mu) = \sum_{i} (i-1) \mu_i`.

            If the parameter ``power`` is an integer then it calculates
            nabla to that integer.  The default value of ``power`` is 1.

            INPUT:

            - ``self`` -- an element of the Macdonald `Ht` basis
            - ``q``, ``t`` -- (optional) parameters to specialize
            - ``power`` -- integer (default: 1)

            OUTPUT: the symmetric function of `\nabla` acting on ``self``

            EXAMPLES::

                sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
                sage: Ht = Sym.macdonald().Ht()
                sage: t = Ht.t; q = Ht.q
                sage: s = Sym.schur()
                sage: a = sum(Ht(p) for p in Partitions(3))
                sage: Ht(0).nabla()
                0
                sage: a.nabla() == t^3*Ht([1,1,1])+q*t*Ht([2,1]) + q^3*Ht([3])
                True
                sage: a.nabla(t=3) == 27*Ht([1,1,1])+3*q*Ht([2,1]) + q^3*Ht([3])
                True
                sage: a.nabla(q=3) == t^3*Ht([1,1,1])+3*t*Ht([2,1]) + 27*Ht([3])
                True
                sage: Ht[2,1].nabla(power=-1)
                1/(q*t)*McdHt[2, 1]
                sage: Ht[2,1].nabla(power=4)
                q^4*t^4*McdHt[2, 1]
                sage: s(a.nabla(q=3))
                (t^6+27*q^3+3*q*t^2)*s[1, 1, 1] + (t^5+t^4+27*q^2+3*q*t+3*t^2+27*q)*s[2, 1] + (t^3+3*t+27)*s[3]
                sage: Ht = Sym.macdonald(q=3).Ht()
                sage: a = sum(Ht(p) for p in Partitions(3))
                sage: s(a.nabla())
                (t^6+9*t^2+729)*s[1, 1, 1] + (t^5+t^4+3*t^2+9*t+324)*s[2, 1] + (t^3+3*t+27)*s[3]
            """
            P = self.parent()
            Ht = P._macdonald.Ht()
            selfHt = Ht(self)
            if self == Ht.zero():
                return Ht.zero()
            if q is None:
                q = Ht.q
            if t is None:
                t = Ht.t
            f = lambda part: t**(part.weighted_size()*power)*q**(part.conjugate().weighted_size()*power)*Ht(part)
            return P(Ht._apply_module_morphism(selfHt, f))


class MacdonaldPolynomials_s(MacdonaldPolynomials_generic):
    def __init__(self, macdonald):
        r"""
        An implementation of the basis `s_\lambda[(1-t)X/(1-q)]`.

        This is perhaps misnamed as a 'Macdonald' basis for
        the symmetric functions but is used in the calculation
        of the Macdonald `J` basis (see method 'creation' below)
        but does use both of the two parameters and can be
        specialized to `s_\lambda[(1-t)X]` and `s_\lambda[X/(1-t)]`.

        INPUT:

        - ``self`` -- a Macdonald `S` basis
        - ``macdonald`` -- a family of Macdonald bases

        TESTS::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: S = Sym.macdonald().S()
            sage: TestSuite(S).run(skip=["_test_associativity","_test_distributivity","_test_prod"])
            sage: TestSuite(S).run(elements = [S.t*S[1,1]+S.q*S[2], S[1]+(S.q+S.t)*S[1,1]])
        """
        MacdonaldPolynomials_generic.__init__(self, macdonald)
        self._s = macdonald._s
        self._self_to_s_cache = _S_to_s_cache
        self._s_to_self_cache = _s_to_S_cache

    def product(self, left, right):
        r"""
        The multiplication of the modified Schur functions behaves the same
        as the multiplication of the Schur functions.

        INPUT:

        - ``self`` -- a Macdonald `S` basis
        - ``left``, ``right`` -- a symmetric functions

        OUTPUT: the product of ``left`` and ``right``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: S = Sym.macdonald().S()
            sage: S([2])^2 #indirect doctest
            McdS[2, 2] + McdS[3, 1] + McdS[4]
        """
        s_left = self._s._from_element(left)
        s_right = self._s._from_element(right)
        product = s_left * s_right
        return self._from_element(product)

    def _to_s(self, part):
        r"""
        Return a function which gives the coefficient of a partition in
        the Schur expansion of ``self(part)``.
        these computations are completed with coefficients in fraction
        field of polynomials in `q` and `t`

        INPUT:

        - ``self`` -- a Macdonald `S` basis
        - ``part`` -- a partition

        OUTPUT:

        - returns a function which accepts a partition ``part2`` and
          this function returns the coefficient of the Schur function
          indexed by ``part2`` in ``S(part)``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: S = Sym.macdonald().S()
            sage: S2 = S._to_s(Partition([2]))
            sage: S2(Partition([2]))
            (-q*t^2 + q*t + t - 1)/(-q^3 + q^2 + q - 1)
            sage: S2(Partition([1,1]))
            (q*t - t^2 - q + t)/(-q^3 + q^2 + q - 1)
        """
        # Convert to the power sum
        (q, t) = QQqt.gens()
        p = self._sym.p()
        s = self._s
        p_x = p(s(part))
        f = lambda m, c: (m, c * prod([(1 - t**k) / (1 - q**k) for k in m]))
        res = s(p_x.map_item(f))
        f = res.coefficient
        return f

    def _s_cache(self, n):
        r"""
        Compute the change of basis and its inverse between the Macdonald
        polynomials on the S basis and the Schur functions.

        These computations are completed with coefficients in fraction
        field of polynomials in `q` and `t`.

        INPUT:

        - ``self`` -- a Macdonald `S` basis
        - ``n`` -- positive integer

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: S = Sym.macdonald().S()
            sage: S._s_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l( S._s_to_self_cache[2] )
            [([1, 1], [([1, 1], (q^2*t - q*t - q + 1)/(t^3 - t^2 - t + 1)), ([2], (q^2 - q*t - q + t)/(t^3 - t^2 - t + 1))]), ([2], [([1, 1], (q^2 - q*t - q + t)/(t^3 - t^2 - t + 1)), ([2], (-q^2*t + q*t + q - 1)/(-t^3 + t^2 + t - 1))])]
            sage: l( S._self_to_s_cache[2] )
            [([1, 1], [([1, 1], (-q*t^2 + q*t + t - 1)/(-q^3 + q^2 + q - 1)), ([2], (q*t - t^2 - q + t)/(-q^3 + q^2 + q - 1))]), ([2], [([1, 1], (q*t - t^2 - q + t)/(-q^3 + q^2 + q - 1)), ([2], (-q*t^2 + q*t + t - 1)/(-q^3 + q^2 + q - 1))])]
        """
        self._invert_morphism(n, QQqt, self._self_to_s_cache,
                              self._s_to_self_cache,
                              to_other_function=self._to_s)

    class Element(MacdonaldPolynomials_generic.Element):

        def _creation_by_determinant_helper(self, k, part):
            r"""
            Formula from [LLM1998]_ Corollary 4.3 p. 970.

            This is part of a formula for a column adding creation operator
            for the `J` basis and its action on the `S` basis.

            INPUT:

            - ``self`` -- an element of the Macdonald `S` basis
            - ``k`` -- positive integer at least as big as the length of ``part``
            - ``part`` -- a partition

            OUTPUT:

            - returns a symmetric function which is the action of a
              column adding operator on the `J` basis

            EXAMPLES::

                sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
                sage: S = Sym.macdonald().S()
                sage: a = S([2,1])
                sage: a._creation_by_determinant_helper(2,[1])
                (q^3*t-q^2*t-q+1)*McdS[2, 1] + (q^3-q^2*t-q+t)*McdS[3]
            """
            (q,t) = QQqt.gens()
            from sage.combinat.sf.sf import SymmetricFunctions
            S = SymmetricFunctions(QQqt).macdonald().S()

            part += [0]*(k-len(part))

            if len(part) > k:
                raise ValueError("the column to add is too small")

            #Create the matrix over the homogeneous symmetric
            #functions and take its determinant
            h = S._sym.homogeneous()
            MS = MatrixSpace(h, k, k)
            m = []
            for i in range(k):
                row = [0]*max(0, (i+1)-2-part[i])
                for j in range(max(0, (i+1)-2-part[i]),k):
                    value = part[i]+j-i+1
                    p = [value] if value > 0 else []
                    row.append( (1-q**(part[i]+j-i+1)*t**(k-(j+1)))*h(p) )
                m.append(row)
            M = MS(m)
            res = M.det()

            #Convert to the Schurs
            res = S._s( res )
            return S._from_element(res)

        def _creation_by_determinant(self, k):
            r"""
            This function is a creation operator for the `J`-basis
            for which the action is known on the Macdonald `S`-basis
            by formula from [LLM1998]_.

            INPUT:

            - ``self`` -- an element of the Macdonald `S` basis
            - ``k`` -- positive integer

            OUTPUT: the column adding operator on the `J` basis on ``self``

            EXAMPLES::

                sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
                sage: S = Sym.macdonald().S()
                sage: a = S(1)
                sage: a = a._creation_by_determinant(1); a
                (-q+1)*McdS[1]
                sage: a = a._creation_by_determinant(3)
                sage: Sym.macdonald().J()(a)
                McdJ[2, 1, 1]
            """
            S = self.parent()
            f = functools.partial(self._creation_by_determinant_helper,k)
            return S._apply_module_morphism(self, f)

        def creation(self, k):
            r"""
            This function is a creation operator for the J-basis
            for which the action is known on the 'Macdonald' S-basis
            by formula from [LLM1998]_.

            INPUT:

            - ``self`` -- an element of the Macdonald `S` basis
            - ``k`` -- positive integer

            OUTPUT: the column adding operator on the `J` basis on ``self``

            EXAMPLES::

                sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
                sage: S = Sym.macdonald().S()
                sage: a = S(1)
                sage: a.creation(1)
                (-q+1)*McdS[1]
                sage: a.creation(2)
                (q^2*t-q*t-q+1)*McdS[1, 1] + (q^2-q*t-q+t)*McdS[2]
            """
            return self._creation_by_determinant(k)

        def _omega_qt_in_schurs(self):
            r"""
            Return the image of ``self`` under the omega_qt automorphism in the
            Schur basis.

            INPUT:

            - ``self`` -- an element of the Macdonald `S` basis

            OUTPUT:

            - the action of the `\omega_{qt}` operator changes the
              Macdonald `S` basis indexed by the partition `\mu` to
              the Schur indexed by the conjugate of `\mu`

            EXAMPLES::

                sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
                sage: S = Sym.macdonald().S()
                sage: a = S([2,1]) + S([1,1,1])
                sage: a._omega_qt_in_schurs()
                s[2, 1] + s[3]
            """
            S = self.parent()
            f = lambda part: S._s(part.conjugate())
            return S._s._apply_module_morphism(self, f)


def qt_kostka(lam, mu):
    r"""
    Return the `K_{\lambda\mu}(q,t)` by computing the change
    of basis from the Macdonald H basis to the Schurs.

    INPUT:

    - ``lam``, ``mu`` -- partitions of the same size

    OUTPUT:

    - returns the `q,t`-Kostka polynomial indexed by the
      partitions ``lam`` and ``mu``

    EXAMPLES::

        sage: from sage.combinat.sf.macdonald import qt_kostka
        sage: qt_kostka([2,1,1],[1,1,1,1])
        t^3 + t^2 + t
        sage: qt_kostka([1,1,1,1],[2,1,1])
        q
        sage: qt_kostka([1,1,1,1],[3,1])
        q^3
        sage: qt_kostka([1,1,1,1],[1,1,1,1])
        1
        sage: qt_kostka([2,1,1],[2,2])
        q^2*t + q*t + q
        sage: qt_kostka([2,2],[2,2])
        q^2*t^2 + 1
        sage: qt_kostka([4],[3,1])
        t
        sage: qt_kostka([2,2],[3,1])
        q^2*t + q
        sage: qt_kostka([3,1],[2,1,1])
        q*t^3 + t^2 + t
        sage: qt_kostka([2,1,1],[2,1,1])
        q*t^2 + q*t + 1
        sage: qt_kostka([2,1],[1,1,1,1])
        0
    """
    lam = _Partitions(lam)
    mu = _Partitions(mu)

    if lam.size() != mu.size():
        return QQqt.zero()

    if (lam,mu) in _qt_kostka_cache:
        return _qt_kostka_cache[(lam,mu)]

    from sage.combinat.sf.sf import SymmetricFunctions
    Sym = SymmetricFunctions(QQqt)
    H = Sym.macdonald().H()
    s = Sym.schur()

    parts = Partitions_n(mu.size())

    for p2 in parts:
        res = s(H(p2))
        for p1 in parts:
            _qt_kostka_cache[(p1, p2)] = QQqt(res.coefficient(p1).numerator())

    return _qt_kostka_cache[(lam, mu)]


# Backward compatibility for unpickling
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.sf.macdonald', 'MacdonaldPolynomial_h', MacdonaldPolynomials_h.Element)
register_unpickle_override('sage.combinat.sf.macdonald', 'MacdonaldPolynomial_ht', MacdonaldPolynomials_ht.Element)
register_unpickle_override('sage.combinat.sf.macdonald', 'MacdonaldPolynomial_j', MacdonaldPolynomials_j.Element)
register_unpickle_override('sage.combinat.sf.macdonald', 'MacdonaldPolynomial_p', MacdonaldPolynomials_p.Element)
register_unpickle_override('sage.combinat.sf.macdonald', 'MacdonaldPolynomial_q', MacdonaldPolynomials_q.Element)
register_unpickle_override('sage.combinat.sf.macdonald', 'MacdonaldPolynomial_s', MacdonaldPolynomials_s.Element)
