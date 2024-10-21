# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Hall-Littlewood Polynomials

Notation used in the definitions follows mainly [Mac1995]_.
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
from . import sfa
import sage.combinat.partition
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.misc.lazy_import import lazy_import
from sage.rings.rational_field import QQ

lazy_import('sage.libs.symmetrica.all', 'hall_littlewood')
lazy_import('sage.matrix.constructor', 'matrix')


# P basis cache
p_to_s_cache = {}
s_to_p_cache = {}
# Qp basis cache
qp_to_s_cache = {}
s_to_qp_cache = {}

QQt = QQ['t'].fraction_field()

# TODO: optimize! which is the fastest way of computing HL's and kostka-polynomials?
# Qp basis is computed using symmetrica, while P basis is computed using rigged
# configurations


class HallLittlewood(UniqueRepresentation):
    r"""
    The family of Hall-Littlewood symmetric function bases.

    The Hall-Littlewood symmetric functions are a family of symmetric
    functions that depend on a parameter `t`.

    INPUT:

    By default the parameter for these functions is `t`, and
    whatever the parameter is, it must be in the base ring.

    EXAMPLES::

        sage: SymmetricFunctions(QQ).hall_littlewood(1)
        Hall-Littlewood polynomials with t=1 over Rational Field
        sage: SymmetricFunctions(QQ['t'].fraction_field()).hall_littlewood()
        Hall-Littlewood polynomials over Fraction Field of Univariate Polynomial Ring in t over Rational Field
    """

    def __repr__(self):
        r"""
        A string representing the family of Hall-Littlewood symmetric function bases.

        OUTPUT: string representing the class

        EXAMPLES::

            sage: SymmetricFunctions(QQ).hall_littlewood(1)
            Hall-Littlewood polynomials with t=1 over Rational Field
        """
        return self._name + " over %s" % self._sym.base_ring()

    @staticmethod
    def __classcall__(cls, Sym, t='t'):
        """
        Normalize the arguments.

        TESTS::

            sage: R.<q, t> = QQ[]
            sage: B1 = SymmetricFunctions(R).hall_littlewood()
            sage: B2 = SymmetricFunctions(R).hall_littlewood(t)
            sage: B3 = SymmetricFunctions(R).hall_littlewood(q)
            sage: B1 is B2
            True
            sage: B1 == B3
            False
        """
        return super().__classcall__(cls, Sym, Sym.base_ring()(t))

    def __init__(self, Sym, t):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: HL = SymmetricFunctions(FractionField(QQ['t'])).hall_littlewood()
            sage: TestSuite(HL).run()
        """
        self._sym = Sym
        self.t = Sym.base_ring()(t)
        self._name_suffix = ""
        if str(t) != 't':
            self._name_suffix += " with t=%s" % t
        self._name = "Hall-Littlewood polynomials"+self._name_suffix

    def symmetric_function_ring( self ):
        r"""
        Return the ring of symmetric functions associated to the class of
        Hall-Littlewood symmetric functions.

        INPUT:

        - ``self`` -- a class of Hall-Littlewood symmetric function bases

        EXAMPLES::

            sage: HL = SymmetricFunctions(FractionField(QQ['t'])).hall_littlewood()
            sage: HL.symmetric_function_ring()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        return self._sym

    def base_ring( self ):
        r"""
        Return the base ring of the symmetric functions where the
        Hall-Littlewood symmetric functions live.

        INPUT:

        - ``self`` -- a class of Hall-Littlewood symmetric function bases

        EXAMPLES::

            sage: HL = SymmetricFunctions(QQ['t'].fraction_field()).hall_littlewood(t=1)
            sage: HL.base_ring()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        return self._sym.base_ring()

    def P(self):
        r"""
        Return the algebra of symmetric functions in the Hall-Littlewood
        `P` basis. This is the same as the `HL` basis in John Stembridge's
        SF examples file.

        INPUT:

        - ``self`` -- a class of Hall-Littlewood symmetric function bases

        OUTPUT: the class of the Hall-Littlewood `P` basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLP = Sym.hall_littlewood().P(); HLP
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Hall-Littlewood P basis
            sage: SP = Sym.hall_littlewood(t=-1).P(); SP
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Hall-Littlewood P with t=-1 basis
            sage: s = Sym.schur()
            sage: s(HLP([2,1]))
            (-t^2-t)*s[1, 1, 1] + s[2, 1]

        The Hall-Littlewood polynomials in the `P` basis at `t = 0` are the
        Schur functions::

            sage: Sym = SymmetricFunctions(QQ)
            sage: HLP = Sym.hall_littlewood(t=0).P()
            sage: s = Sym.schur()
            sage: s(HLP([2,1])) == s([2,1])
            True

        The Hall-Littlewood polynomials in the `P` basis at `t = 1` are the
        monomial symmetric functions::

            sage: Sym = SymmetricFunctions(QQ)
            sage: HLP = Sym.hall_littlewood(t=1).P()
            sage: m = Sym.monomial()
            sage: m(HLP([2,2,1])) == m([2,2,1])
            True

        We end with some examples of coercions between:

            1. Hall-Littlewood `P` basis.

            2. Hall-Littlewood polynomials in the `Q` basis

            3. Hall-Littlewood polynomials in the `Q^\prime` basis (via the Schurs)

            4. Classical symmetric functions

        ::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLP  = Sym.hall_littlewood().P()
            sage: HLQ  = Sym.hall_littlewood().Q()
            sage: HLQp = Sym.hall_littlewood().Qp()
            sage: s = Sym.schur()
            sage: p = Sym.power()
            sage: HLP(HLQ([2])) # indirect doctest
            (-t+1)*HLP[2]
            sage: HLP(HLQp([2]))
            t*HLP[1, 1] + HLP[2]
            sage: HLP(s([2]))
            t*HLP[1, 1] + HLP[2]
            sage: HLP(p([2]))
            (t-1)*HLP[1, 1] + HLP[2]
            sage: s = HLQp.symmetric_function_ring().s()
            sage: HLQp.transition_matrix(s,3)
            [      1       0       0]
            [      t       1       0]
            [    t^3 t^2 + t       1]
            sage: s.transition_matrix(HLP,3)
            [      1       t     t^3]
            [      0       1 t^2 + t]
            [      0       0       1]

        The method :meth:`sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.hl_creation_operator`
        is a creation operator for the `Q` basis::

            sage: HLQp[1].hl_creation_operator([3]).hl_creation_operator([3])
            HLQp[3, 3, 1]

        Transitions between bases with the parameter `t` specialized::

            sage: Sym = SymmetricFunctions(FractionField(QQ['y','z']))
            sage: (y,z) = Sym.base_ring().gens()
            sage: HLy = Sym.hall_littlewood(t=y)
            sage: HLz = Sym.hall_littlewood(t=z)
            sage: Qpy = HLy.Qp()
            sage: Qpz = HLz.Qp()
            sage: s = Sym.schur()
            sage: s( Qpy[3,1] + z*Qpy[2,2] )
            z*s[2, 2] + (y*z+1)*s[3, 1] + (y^2*z+y)*s[4]
            sage: s( Qpy[3,1] + y*Qpz[2,2] )
            y*s[2, 2] + (y*z+1)*s[3, 1] + (y*z^2+y)*s[4]
            sage: s( Qpy[3,1] + y*Qpy[2,2] )
            y*s[2, 2] + (y^2+1)*s[3, 1] + (y^3+y)*s[4]

            sage: Qy = HLy.Q()
            sage: Qz = HLz.Q()
            sage: Py = HLy.P()
            sage: Pz = HLz.P()
            sage: Pz(Qpy[2,1])
            (y*z^3+z^2+z)*HLP[1, 1, 1] + (y*z+1)*HLP[2, 1] + y*HLP[3]
            sage: Pz(Qz[2,1])
            (z^2-2*z+1)*HLP[2, 1]
            sage: Qz(Py[2])
            ((-y+z)/(z^3-z^2-z+1))*HLQ[1, 1] + (1/(-z+1))*HLQ[2]
            sage: Qy(Pz[2])
            ((y-z)/(y^3-y^2-y+1))*HLQ[1, 1] + (1/(-y+1))*HLQ[2]
            sage: Qy.hall_littlewood_family() == HLy
            True
            sage: Qy.hall_littlewood_family() == HLz
            False
            sage: Qz.symmetric_function_ring() == Qy.symmetric_function_ring()
            True

            sage: Sym = SymmetricFunctions(FractionField(QQ['q']))
            sage: q = Sym.base_ring().gen()
            sage: HL = Sym.hall_littlewood(t=q)
            sage: HLQp = HL.Qp()
            sage: HLQ = HL.Q()
            sage: HLP = HL.P()
            sage: s = Sym.schur()
            sage: s(HLQp[3,2].plethysm((1-q)*s[1]))/(1-q)^2
            (-q^5-q^4)*s[1, 1, 1, 1, 1] + (q^3+q^2)*s[2, 1, 1, 1] - q*s[2, 2, 1] - q*s[3, 1, 1] + s[3, 2]
            sage: s(HLP[3,2])
            (-q^5-q^4)*s[1, 1, 1, 1, 1] + (q^3+q^2)*s[2, 1, 1, 1] - q*s[2, 2, 1] - q*s[3, 1, 1] + s[3, 2]

        The `P` and `Q`-Schur at `t=-1` indexed by strict partitions are a basis for
        the space algebraically generated by the odd power sum symmetric functions::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q']))
            sage: SP = Sym.hall_littlewood(t=-1).P()
            sage: SQ = Sym.hall_littlewood(t=-1).Q()
            sage: p = Sym.power()
            sage: SP(SQ[3,2,1])
            8*HLP[3, 2, 1]
            sage: SP(SQ[2,2,1])
            0
            sage: p(SP[3,2,1])
            1/45*p[1, 1, 1, 1, 1, 1] - 1/9*p[3, 1, 1, 1] - 1/9*p[3, 3] + 1/5*p[5, 1]
            sage: SP(p[3,3])
            -4*HLP[3, 2, 1] + 2*HLP[4, 2] - 2*HLP[5, 1] + HLP[6]
            sage: SQ( SQ[1]*SQ[3] -2*(1-q)*SQ[4] )
            HLQ[3, 1] + 2*q*HLQ[4]

        TESTS::

            sage: HLP(s[[]])
            HLP[]
            sage: HLQ(s[[]])
            HLQ[]
            sage: HLQp(s[[]])
            HLQp[]
        """
        return HallLittlewood_p(self)

    def Q(self):
        r"""
        Return the algebra of symmetric functions in Hall-Littlewood `Q`
        basis. This is the same as the `Q` basis in John Stembridge's SF
        examples file.

        More extensive examples can be found in the documentation for the
        Hall-Littlewood `P` basis.

        INPUT:

        - ``self`` -- a class of Hall-Littlewood symmetric function bases

        OUTPUT: the class of the Hall-Littlewood `Q` basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLQ = Sym.hall_littlewood().Q(); HLQ
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Hall-Littlewood Q basis
            sage: SQ = SymmetricFunctions(QQ).hall_littlewood(t=-1).Q(); SQ
            Symmetric Functions over Rational Field in the Hall-Littlewood Q with t=-1 basis
        """
        return HallLittlewood_q(self)

    def Qp(self):
        r"""
        Return the algebra of symmetric functions in Hall-Littlewood `Q^\prime` (Qp)
        basis. This is dual to the Hall-Littlewood `P` basis with respect to
        the standard scalar product.

        More extensive examples can be found in the documentation for the
        Hall-Littlewood P basis.

        INPUT:

        - ``self`` -- a class of Hall-Littlewood symmetric function bases

        OUTPUT: the class of the Hall-Littlewood `Qp`-basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLQp = Sym.hall_littlewood().Qp(); HLQp
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Hall-Littlewood Qp basis
        """
        return HallLittlewood_qp(self)


class HallLittlewood_generic(sfa.SymmetricFunctionAlgebra_generic):
    def __init__(self, hall_littlewood):
        r"""
        A class with methods for working with Hall-Littlewood symmetric functions which
        are common to all bases.

        INPUT:

        - ``self`` -- a Hall-Littlewood symmetric function basis
        - ``hall_littlewood`` -- a class of Hall-Littlewood bases

        TESTS::

            sage: SymmetricFunctions(QQ['t'].fraction_field()).hall_littlewood().P()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Hall-Littlewood P basis
            sage: SymmetricFunctions(QQ).hall_littlewood(t=2).P()
            Symmetric Functions over Rational Field in the Hall-Littlewood P with t=2 basis
        """
        s = self.__class__.__name__[15:].capitalize()
        sfa.SymmetricFunctionAlgebra_generic.__init__(
            self, hall_littlewood._sym,
            basis_name="Hall-Littlewood " + s + hall_littlewood._name_suffix,
            prefix="HL" + s)
        self.t = hall_littlewood.t
        self._sym = hall_littlewood._sym
        self._hall_littlewood = hall_littlewood
        self._s = self._sym.schur()

        # This coercion is broken: HLP = HallLittlewoodP(QQ); HLP(HLP._s[1])

        # Bases defined by orthotriangularity should inherit from some
        # common category BasesByOrthotriangularity (shared with Jack, HL, orthotriang, Mcdo)
        if hasattr(self, "_s_cache"):
            # temporary until Hom(GradedHopfAlgebrasWithBasis work better)
            category = sage.categories.all.ModulesWithBasis(self._sym.base_ring())
            self   .register_coercion(SetMorphism(Hom(self._s, self, category), self._s_to_self))
            self._s.register_coercion(SetMorphism(Hom(self, self._s, category), self._self_to_s))

    def construction(self):
        """
        Return a pair ``(F, R)``, where ``F`` is a
        :class:`SymmetricFunctionsFunctor` and `R` is a ring, such
        that ``F(R)`` returns ``self``.

        EXAMPLES::

            sage: P = SymmetricFunctions(QQ).hall_littlewood(t=2).P()
            sage: P.construction()
            (SymmetricFunctionsFunctor[Hall-Littlewood P with t=2], Rational Field)
        """

        return (sfa.SymmetricFunctionsFamilyFunctor(self,
                                                    HallLittlewood,
                                                    self.basis_name(),
                                                    self.t),
                self.base_ring())

    def _s_to_self(self, x):
        r"""
        Isomorphism from the Schur basis into ``self``.

        INPUT:

        - ``self`` -- a Hall-Littlewood symmetric function basis
        - ``x`` -- an element of the Schur basis

        OUTPUT: an element of ``self`` equivalent to ``x``

        EXAMPLES::

            sage: P = SymmetricFunctions(QQ).hall_littlewood(t=2).P()
            sage: s = SymmetricFunctions(QQ).schur()
            sage: P._s_to_self(s[2,1])
            6*HLP[1, 1, 1] + HLP[2, 1]

        This is for internal use only. Please use instead::

            sage: P(s[2,1])
            6*HLP[1, 1, 1] + HLP[2, 1]
        """
        return self._from_cache(x, self._s_cache, self._s_to_self_cache,
                                t=self.t)

    def _self_to_s(self, x):
        r"""
        Isomorphism from ``self`` to the Schur basis.

        INPUT:

        - ``self`` -- a Hall-Littlewood symmetric function basis
        - ``x`` -- an element of the basis ``self``

        OUTPUT: an element of the Schur basis equivalent to ``x``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: P = Sym.hall_littlewood(t=2).P()
            sage: s = Sym.schur()
            sage: P._self_to_s(P[2,1])
            -6*s[1, 1, 1] + s[2, 1]

        This is for internal use only. Please use instead::

            sage: s(P[2,1])
            -6*s[1, 1, 1] + s[2, 1]
        """
        return self._s._from_cache(x, self._s_cache, self._self_to_s_cache,
                                   t=self.t)

    def transition_matrix(self, basis, n):
        r"""
        Return the transitions matrix between ``self`` and ``basis`` for the
        homogeneous component of degree ``n``.

        INPUT:

        - ``self`` -- a Hall-Littlewood symmetric function basis
        - ``basis`` -- another symmetric function basis
        - ``n`` -- nonnegative integer representing the degree

        OUTPUT:

        - Returns a `r \times r` matrix of elements of the base ring of ``self``
          where `r` is the number of partitions of ``n``.
          The entry corresponding to row `\mu`, column `\nu` is the
          coefficient of ``basis`` `(\nu)` in ``self`` `(\mu)`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLP = Sym.hall_littlewood().P()
            sage: s   = Sym.schur()
            sage: HLP.transition_matrix(s, 4)
            [             1             -t              0            t^2           -t^3]
            [             0              1             -t             -t      t^3 + t^2]
            [             0              0              1             -t            t^3]
            [             0              0              0              1 -t^3 - t^2 - t]
            [             0              0              0              0              1]
            sage: HLQ = Sym.hall_littlewood().Q()
            sage: HLQ.transition_matrix(s,3)
            [                        -t + 1                        t^2 - t                     -t^3 + t^2]
            [                             0                  t^2 - 2*t + 1           -t^4 + t^3 + t^2 - t]
            [                             0                              0 -t^6 + t^5 + t^4 - t^2 - t + 1]
            sage: HLQp = Sym.hall_littlewood().Qp()
            sage: HLQp.transition_matrix(s,3)
            [      1       0       0]
            [      t       1       0]
            [    t^3 t^2 + t       1]
        """
        P = sage.combinat.partition.Partitions_n(n)
        Plist = P.list()
        m = []
        for row_part in Plist:
            z = basis(self(row_part))
            m.append( [z.coefficient(col_part) for col_part in Plist] )
        return matrix(m)

    def product(self, left, right):
        r"""
        Multiply an element of the Hall-Littlewood symmetric function
        basis ``self`` and another symmetric function

        Convert to the Schur basis, do the multiplication there, and
        convert back to ``self`` basis.

        INPUT:

        - ``self`` -- a Hall-Littlewood symmetric function basis
        - ``left`` -- an element of the basis ``self``
        - ``right`` -- another symmetric function

        OUTPUT: the product of ``left`` and ``right`` expanded in the basis ``self``

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLP = Sym.hall_littlewood().P()
            sage: HLP([2])^2 # indirect doctest
            (t+1)*HLP[2, 2] + (-t+1)*HLP[3, 1] + HLP[4]

            sage: HLQ = Sym.hall_littlewood().Q()
            sage: HLQ([2])^2 # indirect doctest
            HLQ[2, 2] + (-t+1)*HLQ[3, 1] + (-t+1)*HLQ[4]

            sage: HLQp = Sym.hall_littlewood().Qp()
            sage: HLQp([2])^2 # indirect doctest
            HLQp[2, 2] + (-t+1)*HLQp[3, 1] + (-t+1)*HLQp[4]
        """
        return self(self._s(left) * self._s(right))

    def hall_littlewood_family(self):
        r"""
        The family of Hall-Littlewood bases associated to ``self``.

        INPUT:

        - ``self`` -- a Hall-Littlewood symmetric function basis

        OUTPUT: the class of Hall-Littlewood bases

        EXAMPLES::

            sage: HLP = SymmetricFunctions(FractionField(QQ['t'])).hall_littlewood(1).P()
            sage: HLP.hall_littlewood_family()
            Hall-Littlewood polynomials with t=1 over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        return self._hall_littlewood

    class Element(sfa.SymmetricFunctionAlgebra_generic.Element):
        r"""
        Methods for elements of a Hall-Littlewood basis that are common to all bases.
        """

        def expand(self, n, alphabet='x'):
            r"""
            Expand the symmetric function as a symmetric polynomial in ``n`` variables.

            INPUT:

            - ``self`` -- an element of a Hall-Littlewood basis
            - ``n`` -- positive integer
            - ``alphabet`` -- string representing a variable name (default: ``'x'``)

            OUTPUT: a symmetric polynomial of ``self`` in ``n`` variables

            EXAMPLES::

                sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
                sage: HLP = Sym.hall_littlewood().P()
                sage: HLQ = Sym.hall_littlewood().Q()
                sage: HLQp = Sym.hall_littlewood().Qp()
                sage: HLP([2]).expand(2)
                x0^2 + (-t + 1)*x0*x1 + x1^2
                sage: HLQ([2]).expand(2)
                (-t + 1)*x0^2 + (t^2 - 2*t + 1)*x0*x1 + (-t + 1)*x1^2
                sage: HLQp([2]).expand(2)
                x0^2 + x0*x1 + x1^2
                sage: HLQp([2]).expand(2, 'y')
                y0^2 + y0*y1 + y1^2
                sage: HLQp([2]).expand(1)
                x^2
            """
            s = self.parent().realization_of().schur()
            return s(self).expand(n, alphabet=alphabet)

        def scalar(self, x, zee=None):
            r"""
            Return standard scalar product between ``self`` and ``x``.

            This is the default implementation that converts both ``self`` and ``x``
            into Schur functions and performs the scalar product that basis.

            The Hall-Littlewood `P` basis is dual to the `Qp` basis with respect to
            this scalar product.

            INPUT:

            - ``self`` -- an element of a Hall-Littlewood basis
            - ``x`` -- another symmetric element of the symmetric functions

            OUTPUT: the scalar product between ``self`` and ``x``

            EXAMPLES::

                sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
                sage: HLP = Sym.hall_littlewood().P()
                sage: HLQ = Sym.hall_littlewood().Q()
                sage: HLQp = Sym.hall_littlewood().Qp()
                sage: HLP([2]).scalar(HLQp([2]))
                1
                sage: HLP([2]).scalar(HLQp([1,1]))
                0
                sage: HLP([2]).scalar(HLQ([2]), lambda mu: mu.centralizer_size(t = HLP.t))
                1
                sage: HLP([2]).scalar(HLQ([1,1]), lambda mu: mu.centralizer_size(t = HLP.t))
                0
            """
            s = self.parent().realization_of().schur()
            s_self = s(self)
            s_x = s(x)
            return s_self.scalar(s_x, zee)

        def scalar_hl(self, x, t=None):
            r"""
            Return the Hall-Littlewood (with parameter ``t``) scalar product
            of ``self`` and ``x``.

            The Hall-Littlewood scalar product is defined in Macdonald's
            book [Mac1995]_.  The power sum basis is orthogonal and
            `\langle p_\mu, p_\mu \rangle = z_\mu \prod_{i} 1/(1-t^{\mu_i})`

            The Hall-Littlewood `P` basis is dual to the `Q` basis with respect to
            this scalar product.

            INPUT:

            - ``self`` -- an element of a Hall-Littlewood basis
            - ``x`` -- another symmetric element of the symmetric functions
            - ``t`` -- an optional parameter, if this parameter is not specified then
              the value of the ``t`` from the basis is used in the calculation

            OUTPUT: the Hall-Littlewood scalar product between ``self`` and ``x``

            EXAMPLES::

                sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
                sage: HLP = Sym.hall_littlewood().P()
                sage: HLQ = Sym.hall_littlewood().Q()
                sage: HLP([2]).scalar_hl(HLQ([2]))
                1
                sage: HLP([2]).scalar_hl(HLQ([1,1]))
                0
                sage: HLQ([2]).scalar_hl(HLQ([2]))
                -t + 1
                sage: HLQ([2]).scalar_hl(HLQ([1,1]))
                0
                sage: HLP([2]).scalar_hl(HLP([2]))
                -1/(t - 1)
            """
            parent = self.parent()
            if t is None:
                t = parent.t
            p = parent.realization_of().power()
            f = lambda part1, part2: part1.centralizer_size(t=t)
            return parent._apply_multi_module_morphism(p(self), p(x), f,
                                                       orthogonal=True)


###########
# P basis #
###########

class HallLittlewood_p(HallLittlewood_generic):
    r"""
    A class representing the Hall-Littlewood `P` basis of symmetric functions
    """

    class Element(HallLittlewood_generic.Element):
        pass

    def __init__(self, hall_littlewood):
        r"""
        A class with methods for working with the Hall-Littlewood `P` basis.

        The `P` basis is calculated from the Schur basis using the functions
        in :meth:`sage.combinat.sf.kfpoly`.  These functions calculate Kostka-Foulkes polynomials
        using rigged configuration formulas.

        This change of basis is inverted to convert to the Schur basis.

        INPUT:

        - ``self`` -- an instance of the Hall-Littlewood `P` basis
        - ``hall_littlewood`` -- a class for the family of Hall-Littlewood bases

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: P = Sym.hall_littlewood().P()
            sage: TestSuite(P).run(skip=['_test_associativity', '_test_distributivity', '_test_prod']) # products are too expensive
            sage: TestSuite(P).run(elements = [P.t*P[1,1]+P[2], P[1]+(1+P.t)*P[1,1]])
        """
        HallLittlewood_generic.__init__(self, hall_littlewood)
        self._self_to_s_cache = p_to_s_cache
        self._s_to_self_cache = s_to_p_cache

    def _q_to_p_normalization(self, m):
        r"""
        The coefficient relating the `Q` and the `P` bases.

        Returns the scalar coefficient that is used when converting from the
        `Q` basis to the `P` basis. Note that this assumes that ``m`` is a
        Partition object.

        INPUT:

        - ``self`` -- an instance of the Hall-Littlewood `P` basis
        - ``m`` -- a partition

        OUTPUT:

        - returns the coefficient equal to `Q(m)/P(m)`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLP = Sym.hall_littlewood().P()
            sage: HLP._q_to_p_normalization(Partition([2,1]))
            t^2 - 2*t + 1
        """
        t = self.t
        coeff = (1-t)**len(m)
        for i in m.to_exp():
            for j in range(1,i+1):
                coeff *= (1-t**j)/(1-t)
        return coeff

    def _s_to_self_base(self, part):
        r"""
        Return a function which gives the coefficient of a partition
        in the expansion of the Schur functions ``s(part)`` in the Hall-Littlewood
        `P` basis.

        INPUT:

        - ``self`` -- an instance of the Hall-Littlewood `P` basis
        - ``part`` -- a partition

        OUTPUT:

        - returns a function which accepts a partition ``part2`` and returns
          the coefficient of ``P(part2)`` in ``s(part)``
          This coefficient is the t-Kostka-Foulkes polynomial  `K_{part,part2}(t)`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLP = Sym.hall_littlewood().P()
            sage: f21 = HLP._s_to_self_base(Partition([2,1]))
            sage: [f21(p) for p in Partitions(3)]
            [0, 1, t^2 + t]
        """
        from sage.combinat.sf.kfpoly import schur_to_hl
        t = QQt.gen()
        zero = self.base_ring().zero()
        res_dict = schur_to_hl(part, t)
        f = lambda part2: res_dict.get(part2,zero)
        return f

    def _s_cache(self, n):
        r"""
        Compute the change of basis between the `P` polynomials and the
        Schur functions for partitions of size ``n``.

        Uses the fact that the transformation matrix is upper-triangular in
        order to obtain the inverse transformation.

        INPUT:

        - ``self`` -- an instance of the Hall-Littlewood `P` basis
        - ``n`` -- positive integer

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLP = Sym.hall_littlewood().P()
            sage: HLP._s_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(HLP._s_to_self_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], t), ([2], 1)])]
            sage: l(HLP._self_to_s_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], -t), ([2], 1)])]
            sage: HLP = Sym.hall_littlewood(10).P()
            sage: HLP._s_cache(2)
            sage: l(HLP._s_to_self_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], t), ([2], 1)])]
        """
        self._invert_morphism(n, QQt, self._self_to_s_cache,
                              self._s_to_self_cache, to_self_function=self._s_to_self_base,
                              upper_triangular=True, ones_on_diagonal=True)


###########
# Q basis #
###########

class HallLittlewood_q(HallLittlewood_generic):
    class Element(HallLittlewood_generic.Element):
        pass

    def __init__(self, hall_littlewood):
        r"""
        The `Q` basis is defined as a normalization of the `P` basis.

        INPUT:

        - ``self`` -- an instance of the Hall-Littlewood `P` basis
        - ``hall_littlewood`` -- a class for the family of Hall-Littlewood bases

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: Q = Sym.hall_littlewood().Q()
            sage: TestSuite(Q).run(skip=['_test_associativity', '_test_distributivity', '_test_prod']) # products are too expensive, long time (3s on sage.math, 2012)
            sage: TestSuite(Q).run(elements = [Q.t*Q[1,1]+Q[2], Q[1]+(1+Q.t)*Q[1,1]])  # long time (depends on previous)

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLP = Sym.hall_littlewood().P()
            sage: HLQ = Sym.hall_littlewood().Q()
            sage: HLQp = Sym.hall_littlewood().Qp()
            sage: s = Sym.schur(); p = Sym.power()
            sage: HLQ( HLP([2,1]) + HLP([3]) )
            (1/(t^2-2*t+1))*HLQ[2, 1] - (1/(t-1))*HLQ[3]
            sage: HLQ(HLQp([2])) # indirect doctest
            (t/(t^3-t^2-t+1))*HLQ[1, 1] - (1/(t-1))*HLQ[2]
            sage: HLQ(s([2]))
            (t/(t^3-t^2-t+1))*HLQ[1, 1] - (1/(t-1))*HLQ[2]
            sage: HLQ(p([2]))
            (1/(t^2-1))*HLQ[1, 1] - (1/(t-1))*HLQ[2]
        """
        HallLittlewood_generic.__init__(self, hall_littlewood)

        self._P = self._hall_littlewood.P()
        # temporary until Hom(GradedHopfAlgebrasWithBasis work better)
        category = sage.categories.all.ModulesWithBasis(self.base_ring())

        phi = self.module_morphism(diagonal=self._P._q_to_p_normalization,
                                   codomain=self._P, category=category)
        self._P.register_coercion(phi)
        self.register_coercion(~phi)

    def _p_to_q_normalization(self, m):
        r"""
        Return the scalar coefficient on self(m) when converting from the
        `Q` basis to the `P` basis. Note that this assumes that ``m`` is a
        Partition object.

        Note: this is not used anymore!

        Returns the scalar coefficient that is used when converting from the
        `P` basis to the `Q` basis. Note that this assumes that ``m`` is a
        Partition object.

        INPUT:

        - ``self`` -- an instance of the Hall-Littlewood `P` basis
        - ``m`` -- a partition

        OUTPUT:

        - returns the coefficient equal to `P(m)/Q(m)`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLQ = Sym.hall_littlewood().Q()
            sage: HLQ._p_to_q_normalization(Partition([2,1]))
            1/(t^2 - 2*t + 1)
        """
        t = self.t
        coeff = 1/(1-t)**len(m)
        for i in m.to_exp():
            for j in range(1,i+1):
                coeff *= (1-t)/(1-t**j)
        return coeff


############
# Qp basis #
############

class HallLittlewood_qp(HallLittlewood_generic):

    class Element(HallLittlewood_generic.Element):
        pass

    def __init__(self, hall_littlewood):
        r"""
        The Hall-Littlewood `Qp` basis is calculated through the symmetrica
        library (see the function :meth:`HallLittlewood_qp._to_s`).

        INPUT:

        - ``self`` -- an instance of the Hall-Littlewood `P` basis
        - ``hall_littlewood`` -- a class for the family of Hall-Littlewood bases

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: Qp = Sym.hall_littlewood().Q()
            sage: TestSuite(Qp).run(skip=['_test_passociativity', '_test_distributivity', '_test_prod']) # products are too expensive, long time (3s on sage.math, 2012)
            sage: TestSuite(Qp).run(elements = [Qp.t*Qp[1,1]+Qp[2], Qp[1]+(1+Qp.t)*Qp[1,1]])  # long time (depends on previous)

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLP = Sym.hall_littlewood().P()
            sage: HLQ = Sym.hall_littlewood().Q()
            sage: HLQp = Sym.hall_littlewood().Qp()
            sage: s = Sym.schur(); p = Sym.power()
            sage: HLQp(HLP([2])) # indirect doctest
            -t*HLQp[1, 1] + (t^2+1)*HLQp[2]
            sage: HLQp(s(HLQ([2]))) # work around bug reported in issue #12969
            (t^2-t)*HLQp[1, 1] + (-t^3+t^2-t+1)*HLQp[2]
            sage: HLQp(s([2]))
            HLQp[2]
            sage: HLQp(p([2]))
            -HLQp[1, 1] + (t+1)*HLQp[2]
            sage: s = HLQp.symmetric_function_ring().s()
            sage: HLQp.transition_matrix(s,3)
            [      1       0       0]
            [      t       1       0]
            [    t^3 t^2 + t       1]
            sage: s.transition_matrix(HLP,3)
            [      1       t     t^3]
            [      0       1 t^2 + t]
            [      0       0       1]
        """
        HallLittlewood_generic.__init__(self, hall_littlewood)
        self._self_to_s_cache = qp_to_s_cache
        self._s_to_self_cache = s_to_qp_cache

    def _to_s(self, part):
        r"""
        Return a function which gives the coefficient of a partition
        in the Schur expansion of ``self(part)``.

        INPUT:

        - ``self`` -- an instance of the Hall-Littlewood `P` basis
        - ``part`` -- a partition

        OUTPUT:

        - returns a function which accepts a second partition ``part2``
          and returns the coefficient of the expansion of the `Qp`
          in the Schur basis.  This is the `t`-Kostka-Foulkes polynomial
          `K_{part2,part}(t)`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLQp = Sym.hall_littlewood().Qp()
            sage: f21 = HLQp._to_s(Partition([2,1]))
            sage: [f21(p) for p in Partitions(3)]
            [t, 1, 0]
        """
        t = QQt.gen()

        if not part:
            return lambda part2: QQt.one()

        res = hall_littlewood(part) # call to symmetrica (returns in variable x)
        f = lambda part2: res.coefficient(part2).subs(x=t)
        return f

    def _s_cache(self, n):
        r"""
        Compute the change of basis between the `Q^\prime` polynomials and the
        Schur functions for partitions of size ``n``.

        Uses the fact that the transformation matrix is lower-triangular in
        order to obtain the inverse transformation.

        INPUT:

        - ``self`` -- an instance of the Hall-Littlewood `P` basis
        - ``n`` -- positive integer

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLQp = Sym.hall_littlewood().Qp()
            sage: HLQp._s_cache(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(HLQp._s_to_self_cache[2])
            [([1, 1], [([1, 1], 1), ([2], -t)]), ([2], [([2], 1)])]
            sage: l(HLQp._self_to_s_cache[2])
            [([1, 1], [([1, 1], 1), ([2], t)]), ([2], [([2], 1)])]
        """
        self._invert_morphism(n, QQt, self._self_to_s_cache,
                              self._s_to_self_cache,
                              to_other_function=self._to_s,
                              lower_triangular=True, ones_on_diagonal=True)


# Unpickling backward compatibility
sage.misc.persist.register_unpickle_override('sage.combinat.sf.hall_littlewood', 'HallLittlewoodElement_p', HallLittlewood_p.Element)
sage.misc.persist.register_unpickle_override('sage.combinat.sf.hall_littlewood', 'HallLittlewoodElement_q', HallLittlewood_q.Element)
sage.misc.persist.register_unpickle_override('sage.combinat.sf.hall_littlewood', 'HallLittlewoodElement_qp', HallLittlewood_qp.Element)
