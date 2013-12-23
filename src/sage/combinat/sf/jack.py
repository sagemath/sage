r"""
Jack Symmetric Functions

Jack's symmetric functions appear in [Ma1995]_ Chapter VI, section 10.
Zonal polynomials are the subject of [Ma1995]_ Chapter VII.
The parameter `\alpha` in that reference is the parameter `t` in this
implementation in sage.

REFERENCES:

.. [Jack1970] H. Jack,
   *A class of symmetric functions with a parameter*,
   Proc. R. Soc. Edinburgh (A), 69, 1-18.

.. [Ma1995] I. G. Macdonald,
   *Symmetric functions and Hall polynomials*,
   second ed.,
   The Clarendon Press, Oxford University Press, New York, 1995, With contributions
   by A. Zelevinsky, Oxford Science Publications.
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2012 Mike Zabrocki <mike.zabrocki@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.unique_representation import UniqueRepresentation
from sage.calculus.var import var
import sage.categories.all
from sage.rings.all import Integer, gcd, lcm, QQ
from sage.rings.fraction_field import is_FractionField
from sage.misc.misc import prod
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom, End
from sage.rings.fraction_field import FractionField
import sfa

QQt = FractionField(QQ['t'])

p_to_m_cache = {}
m_to_p_cache = {}
class Jack(UniqueRepresentation):

    def __init__(self, Sym, t='t'):
        r"""
        The family of Jack symmetric functions including the `P`, `Q`, `J`, `Qp`
        bases.  The default parameter is ``t``.

        INPUT:

        - ``self`` -- the family of Jack symmetric function bases
        - ``Sym`` -- a ring of symmetric functions
        - ``t`` -- an optional parameter (default : 't')

        EXAMPLES::

            sage: SymmetricFunctions(FractionField(QQ['t'])).jack()
            Jack polynomials over Fraction Field of Univariate Polynomial Ring in t over Rational Field
            sage: SymmetricFunctions(QQ).jack(1)
            Jack polynomials with t=1 over Rational Field
        """
        self._sym = Sym
        if not (t in Sym.base_ring() or var(t) in Sym.base_ring()):
            raise ValueError, "parameter t must be in the base ring"
        self.t = Sym.base_ring()(t)
        self._name_suffix = ""
        if str(t) !='t':
            self._name_suffix += " with t=%s"%t
        self._name = "Jack polynomials"+self._name_suffix+" over "+Sym.base_ring().__repr__()

    def __repr__(self):
        r"""
        The string representation for the family of Jack symmetric function bases

        INPUT:

        - ``self`` -- the family of Jack symmetric function bases

        OUTPUT:

        - returns the name of the family of bases

        EXAMPLES::

            sage: SymmetricFunctions(QQ).jack(1)
            Jack polynomials with t=1 over Rational Field
        """
        return self._name

    def base_ring( self ):
        r"""
        Returns the base ring of the symmetric functions in which the
        Jack symmetric functions live

        INPUT:

        - ``self`` -- the family of Jack symmetric function bases

        OUTPUT:

        - the base ring of the symmetric functions ring of ``self``

        EXAMPLES::

            sage: J2 = SymmetricFunctions(QQ).jack(t=2)
            sage: J2.base_ring()
            Rational Field
        """
        return self._sym.base_ring()

    def symmetric_function_ring( self ):
        r"""
        Returns the base ring of the symmetric functions of the Jack symmetric
        function bases

        INPUT:

        - ``self`` -- the family of Jack symmetric function bases

        OUTPUT:

        - the symmetric functions ring of ``self``

        EXAMPLES::

            sage: Jacks = SymmetricFunctions(FractionField(QQ['t'])).jack()
            sage: Jacks.symmetric_function_ring()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        return self._sym

    def P(self):
        r"""
        Returns the algebra of Jack polynomials in the `P` basis.

        INPUT:

        - ``self`` -- the family of Jack symmetric function bases

        OUTPUT:

        - the `P` basis of the Jack symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: JP = Sym.jack().P(); JP
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack P basis
            sage: Sym.jack(t=-1).P()
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack P with t=-1 basis

        At `t = 1`, the Jack polynomials in the `P` basis are the Schur
        symmetric functions.

        ::

            sage: Sym = SymmetricFunctions(QQ)
            sage: JP = Sym.jack(t=1).P()
            sage: s = Sym.schur()
            sage: s(JP([2,2,1]))
            s[2, 2, 1]
            sage: JP(s([2,2,1]))
            JackP[2, 2, 1]
            sage: JP([2,1])^2
            JackP[2, 2, 1, 1] + JackP[2, 2, 2] + JackP[3, 1, 1, 1] + 2*JackP[3, 2, 1] + JackP[3, 3] + JackP[4, 1, 1] + JackP[4, 2]

        At `t = 2`, the Jack polynomials in the `P` basis are the zonal
        polynomials.

        ::

            sage: Sym = SymmetricFunctions(QQ)
            sage: JP = Sym.jack(t=2).P()
            sage: Z = Sym.zonal()
            sage: Z(JP([2,2,1]))
            Z[2, 2, 1]
            sage: JP(Z[2, 2, 1])
            JackP[2, 2, 1]
            sage: JP([2])^2
            64/45*JackP[2, 2] + 16/21*JackP[3, 1] + JackP[4]
            sage: Z([2])^2
            64/45*Z[2, 2] + 16/21*Z[3, 1] + Z[4]

        :::

            sage: Sym = SymmetricFunctions(QQ['a','b'].fraction_field())
            sage: (a,b) = Sym.base_ring().gens()
            sage: Jacka = Sym.jack(t=a)
            sage: Jackb = Sym.jack(t=b)
            sage: m = Sym.monomial()
            sage: JPa = Jacka.P()
            sage: JPb = Jackb.P()
            sage: m(JPa[2,1])
            (6/(a+2))*m[1, 1, 1] + m[2, 1]
            sage: m(JPb[2,1])
            (6/(b+2))*m[1, 1, 1] + m[2, 1]
            sage: m(a*JPb([2,1]) + b*JPa([2,1]))
            ((6*a^2+6*b^2+12*a+12*b)/(a*b+2*a+2*b+4))*m[1, 1, 1] + (a+b)*m[2, 1]
            sage: JPa(JPb([2,1]))
            ((6*a-6*b)/(a*b+2*a+2*b+4))*JackP[1, 1, 1] + JackP[2, 1]

        ::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: JQ = Sym.jack().Q()
            sage: JP = Sym.jack().P()
            sage: JJ = Sym.jack().J()

        ::

            sage: JP(JQ([2,1]))
            ((t+2)/(2*t^3+t^2))*JackP[2, 1]
            sage: JP(JQ([3]))
            ((2*t^2+3*t+1)/(6*t^3))*JackP[3]
            sage: JP(JQ([1,1,1]))
            (6/(t^3+3*t^2+2*t))*JackP[1, 1, 1]

        ::

            sage: JP(JJ([3]))
            (2*t^2+3*t+1)*JackP[3]
            sage: JP(JJ([2,1]))
            (t+2)*JackP[2, 1]
            sage: JP(JJ([1,1,1]))
            6*JackP[1, 1, 1]

        ::

            sage: s = Sym.schur()
            sage: JP(s([2,1]))
            ((2*t-2)/(t+2))*JackP[1, 1, 1] + JackP[2, 1]
            sage: s(_)
            s[2, 1]
        """
        return JackPolynomials_p(self)

    def Q(self):
        r"""
        Returns the algebra of Jack polynomials in the `Q` basis.

        INPUT:

        - ``self`` -- the family of Jack symmetric function bases

        OUTPUT:

        - the `Q` basis of the Jack symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: JQ = Sym.jack().Q(); JQ
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack Q basis
            sage: Sym = SymmetricFunctions(QQ)
            sage: Sym.jack(t=-1).Q()
            Symmetric Functions over Rational Field in the Jack Q with t=-1 basis

        ::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: JQ = Sym.jack().Q()
            sage: JP = Sym.jack().P()
            sage: JQ(sum(JP(p) for p in Partitions(3)))
            (1/6*t^3+1/2*t^2+1/3*t)*JackQ[1, 1, 1] + ((2*t^3+t^2)/(t+2))*JackQ[2, 1] + (6*t^3/(2*t^2+3*t+1))*JackQ[3]

        ::

            sage: s = Sym.schur()
            sage: JQ(s([3])) # indirect doctest
            (1/6*t^3-1/2*t^2+1/3*t)*JackQ[1, 1, 1] + ((2*t^3-2*t^2)/(t+2))*JackQ[2, 1] + (6*t^3/(2*t^2+3*t+1))*JackQ[3]
            sage: JQ(s([2,1]))
            (1/3*t^3-1/3*t)*JackQ[1, 1, 1] + ((2*t^3+t^2)/(t+2))*JackQ[2, 1]
            sage: JQ(s([1,1,1]))
            (1/6*t^3+1/2*t^2+1/3*t)*JackQ[1, 1, 1]
        """
        return JackPolynomials_q(self)

    def J(self):
        r"""
        Returns the algebra of Jack polynomials in the `J` basis.

        INPUT:

        - ``self`` -- the family of Jack symmetric function bases

        OUTPUT: the `J` basis of the Jack symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: JJ = Sym.jack().J(); JJ
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack J basis
            sage: Sym = SymmetricFunctions(QQ)
            sage: Sym.jack(t=-1).J()
            Symmetric Functions over Rational Field in the Jack J with t=-1 basis

        At `t = 1`, the Jack polynomials in the `J` basis are scalar multiples
        of the Schur functions with the scalar given by a Partition's
        :meth:`~sage.combinat.partition.Partition.hook_product` method at 1::

            sage: Sym = SymmetricFunctions(QQ)
            sage: JJ = Sym.jack(t=1).J()
            sage: s = Sym.schur()
            sage: p = Partition([3,2,1,1])
            sage: s(JJ(p)) == p.hook_product(1)*s(p)  # long time (4s on sage.math, 2012)
            True

        At `t = 2`, the Jack polynomials in the `J` basis are scalar multiples
        of the zonal polynomials with the scalar given by a Partition's
        :meth:`~sage.combinat.partition.Partition.hook_product` method at 2.

        ::

            sage: Sym = SymmetricFunctions(QQ)
            sage: JJ = Sym.jack(t=2).J()
            sage: Z = Sym.zonal()
            sage: p = Partition([2,2,1])
            sage: Z(JJ(p)) == p.hook_product(2)*Z(p)
            True

        ::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: JJ = Sym.jack().J()
            sage: JP = Sym.jack().P()
            sage: JJ(sum(JP(p) for p in Partitions(3)))
            1/6*JackJ[1, 1, 1] + (1/(t+2))*JackJ[2, 1] + (1/(2*t^2+3*t+1))*JackJ[3]

        ::

            sage: s = Sym.schur()
            sage: JJ(s([3])) # indirect doctest
            ((t^2-3*t+2)/(6*t^2+18*t+12))*JackJ[1, 1, 1] + ((2*t-2)/(2*t^2+5*t+2))*JackJ[2, 1] + (1/(2*t^2+3*t+1))*JackJ[3]
            sage: JJ(s([2,1]))
            ((t-1)/(3*t+6))*JackJ[1, 1, 1] + (1/(t+2))*JackJ[2, 1]
            sage: JJ(s([1,1,1]))
            1/6*JackJ[1, 1, 1]
        """
        return JackPolynomials_j(self)

    def Qp(self):
        r"""
        Returns the algebra of Jack polynomials in the `Qp`, which is dual to
        the `P` basis with respect to the standard scalar product.

        INPUT:

        - ``self`` -- the family of Jack symmetric function bases

        OUTPUT:

        - the `Q'` basis of the Jack symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: JP = Sym.jack().P()
            sage: JQp = Sym.jack().Qp(); JQp
            Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack Qp basis
            sage: a = JQp([2])
            sage: a.scalar(JP([2]))
            1
            sage: a.scalar(JP([1,1]))
            0
            sage: JP(JQp([2]))                        # todo: missing auto normalization
            ((t-1)/(t+1))*JackP[1, 1] + JackP[2]
            sage: JP._normalize(JP(JQp([2])))
            ((t-1)/(t+1))*JackP[1, 1] + JackP[2]
        """
        return JackPolynomials_qp(self)

###################################################################
def c1(part, t):
    r"""
    Returns the `t`-Jack scalar product between ``J(part)`` and ``P(part)``.

    INPUT:

    - ``part`` -- a partition
    - ``t`` -- an optional parameter (default: uses the parameter `t` from the
      Jack basis)

    OUTPUT:

    - a polynomial in the parameter ``t`` which is equal to the scalar
      product of ``J(part)`` and ``P(part)``

    EXAMPLES::

        sage: from sage.combinat.sf.jack import c1
        sage: t = QQ['t'].gen()
        sage: [c1(p,t) for p in Partitions(3)]
        [2*t^2 + 3*t + 1, t + 2, 6]
    """
    return prod([1+t*part.arm_lengths(flat=True)[i]+part.leg_lengths(flat=True)[i] for i in range(sum(part))],
                t.parent().one())

def c2(part, t):
    r"""
    Returns the t-Jack scalar product between ``J(part)`` and ``Q(part)``.

    INPUT:

    - ``self`` -- a Jack basis of the symmetric functions
    - ``part`` -- a partition
    - ``t`` -- an optional parameter (default: uses the parameter `t` from the
      Jack basis)

    OUTPUT:

    - a polynomial in the parameter ``t`` which is equal to the scalar
      product of ``J(part)`` and ``Q(part)``

    EXAMPLES::

        sage: from sage.combinat.sf.jack import c2
        sage: t = QQ['t'].gen()
        sage: [c2(p,t) for p in Partitions(3)]
        [6*t^3, 2*t^3 + t^2, t^3 + 3*t^2 + 2*t]
    """
    return prod([t+t*part.arm_lengths(flat=True)[i]+part.leg_lengths(flat=True)[i] for i in range(sum(part))],
                t.parent().one())

def normalize_coefficients(self, c):
    r"""
    If our coefficient ring is the field of fractions over a univariate
    polynomial ring over the rationals, then we should clear both the
    numerator and denominator of the denominators of their
    coefficients.

    INPUT:

    - ``self`` -- a Jack basis of the symmetric functions
    - ``c`` -- a coefficient in the base ring of ``self``

    OUTPUT:

    - divide numerator and denominator by the greatest common divisor

    EXAMPLES::

        sage: JP = SymmetricFunctions(FractionField(QQ['t'])).jack().P()
        sage: t = JP.base_ring().gen()
        sage: a = 2/(1/2*t+1/2)
        sage: JP._normalize_coefficients(a)
        4/(t + 1)
        sage: a = 1/(1/3+1/6*t)
        sage: JP._normalize_coefficients(a)
        6/(t + 2)
        sage: a = 24/(4*t^2 + 12*t + 8)
        sage: JP._normalize_coefficients(a)
        6/(t^2 + 3*t + 2)
    """
    BR = self.base_ring()
    if is_FractionField(BR) and BR.base_ring() == QQ:
        denom = c.denominator()
        numer = c.numerator()

        #Clear the denominators
        a = lcm([i.denominator() for i in denom.coeffs()])
        b = lcm([i.denominator() for i in numer.coeffs()])
        l = Integer(a).lcm(Integer(b))
        denom *= l
        numer *= l

        #Divide through by the gcd of the numerators
        a = gcd([i.numerator() for i in denom.coeffs()])
        b = gcd([i.numerator() for i in numer.coeffs()])
        l = Integer(a).gcd(Integer(b))

        denom = denom / l
        numer = numer / l

        return c.parent()(numer, denom)
    else:
        return c

####################################################################

class JackPolynomials_generic(sfa.SymmetricFunctionAlgebra_generic):
    def __init__(self, jack):
        r"""
        A class of methods which are common to all Jack bases of the symmetric functions

        INPUT:

        - ``self`` -- a Jack basis of the symmetric functions
        - ``jack`` -- a family of Jack symmetric function bases

        EXAMPLES ::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: JP = Sym.jack().P(); JP.base_ring()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
            sage: Sym = SymmetricFunctions(QQ)
            sage: JP = Sym.jack(t=2).P(); JP.base_ring()
            Rational Field
        """
        s = self.__class__.__name__[16:].capitalize()
        sfa.SymmetricFunctionAlgebra_generic.__init__(
            self, jack._sym,
            basis_name = "Jack " + s + jack._name_suffix,
            prefix = "Jack"+s)
        self.t = jack.t
        self._sym = jack._sym
        self._jack = jack

        # Bases defined by orthotriangularity should inherit from some
        # common category BasesByOrthotriangularity (shared with Jack, HL, orthotriang, Mcdo)
        if hasattr(self, "_m_cache"):
            # temporary until Hom(GradedHopfAlgebrasWithBasis work better)
            category = sage.categories.all.ModulesWithBasis(self._sym.base_ring())
            self._m = self._sym.monomial()
            self   .register_coercion(SetMorphism(Hom(self._m, self, category), self._m_to_self))
            self._m.register_coercion(SetMorphism(Hom(self, self._m, category), self._self_to_m))
        if hasattr(self, "_h_cache"):
            # temporary until Hom(GradedHopfAlgebrasWithBasis work better)
            category = sage.categories.all.ModulesWithBasis(self._sym.base_ring())
            self._h = self._sym.homogeneous()
            self   .register_coercion(SetMorphism(Hom(self._h, self, category), self._h_to_self))
            self._h.register_coercion(SetMorphism(Hom(self, self._h, category), self._self_to_h))

    def _m_to_self(self, x):
        r"""
        Isomorphism from the monomial basis into ``self``

        INPUT:

        - ``self`` -- a Jack basis of the symmetric functions
        - ``x`` -- element of the monomial basis

        OUTPUT:

        - an element of ``self`` equivalent to ``x``

        EXAMPLES ::

            sage: Sym = SymmetricFunctions(QQ)
            sage: JP = Sym.jack(t=2).P()
            sage: m = Sym.monomial()
            sage: JP._m_to_self(m[2,1])
            -3/2*JackP[1, 1, 1] + JackP[2, 1]

        This is for internal use only. Please use instead::

            sage: JP(m[2,1])
            -3/2*JackP[1, 1, 1] + JackP[2, 1]
        """
        return self._from_cache(x, self._m_cache, self._m_to_self_cache, t = self.t)

    def _self_to_m(self, x):
        r"""
        Isomorphism from self to the monomial basis

        INPUT:

        - ``self`` -- a Jack basis of the symmetric functions
        - ``x`` -- an element of ``self``

        OUTPUT:

        - an element of the monomial basis equivalent to ``x``

        EXAMPLES ::

            sage: Sym = SymmetricFunctions(QQ)
            sage: JP = Sym.jack(t=2).P()
            sage: m = Sym.monomial()
            sage: JP._self_to_m(JP[2,1])
            3/2*m[1, 1, 1] + m[2, 1]

        This is for internal use only. Please use instead::

            sage: m(JP[2,1])
            3/2*m[1, 1, 1] + m[2, 1]
        """
        return self._m._from_cache(x, self._m_cache, self._self_to_m_cache, t = self.t)

    def c1(self, part):
        r"""
        Returns the `t`-Jack scalar product between ``J(part)`` and ``P(part)``.

        INPUT:

        - ``self`` -- a Jack basis of the symmetric functions
        - ``part`` -- a partition
        - ``t`` -- an optional parameter (default: uses the parameter `t` from the
          Jack basis)

        OUTPUT:

        - a polynomial in the parameter ``t`` which is equal to the scalar
          product of ``J(part)`` and ``P(part)``

        EXAMPLES ::

            sage: JP = SymmetricFunctions(FractionField(QQ['t'])).jack().P()
            sage: JP.c1(Partition([2,1]))
            t + 2
        """
        return c1(part, self.t)

    def c2(self, part):
        r"""
        Returns the `t`-Jack scalar product between ``J(part)`` and ``Q(part)``.

        INPUT:

        - ``self`` -- a Jack basis of the symmetric functions
        - ``part`` -- a partition
        - ``t`` -- an optional parameter (default: uses the parameter `t` from the
            Jack basis)

        OUTPUT:

        - a polynomial in the parameter ``t`` which is equal to the scalar
          product of ``J(part)`` and ``Q(part)``

        EXAMPLES::

            sage: JP = SymmetricFunctions(FractionField(QQ['t'])).jack().P()
            sage: JP.c2(Partition([2,1]))
            2*t^3 + t^2
        """
        return c2(part, self.t)

    _normalize_coefficients = normalize_coefficients

    def _normalize(self, x):
        r"""
        Normalize the coefficients of ``x``

        INPUT:

        - ``self`` -- a Jack basis of the symmetric functions
        - ``x`` -- an element of ``self``

        OUTPUT:

        - returns ``x`` with _normalize_coefficient applied to each of the coefficients

        EXAMPLES::

            sage: JP = SymmetricFunctions(FractionField(QQ['t'])).jack().P()
            sage: t = JP.base_ring().gen()
            sage: a = 2/(1/2*t+1/2)
            sage: b = 1/(1/3+1/6*t)
            sage: c = 24/(4*t^2 + 12*t + 8)
            sage: JP._normalize( a*JP[1] + b*JP[2] + c*JP[2,1] )
            (4/(t+1))*JackP[1] + (6/(t+2))*JackP[2] + (6/(t^2+3*t+2))*JackP[2, 1]

        .. todo:: this should be a method on the elements (what's the standard name for such methods?)
        """
        return x.map_coefficients(self._normalize_coefficients)

    def _normalize_morphism(self, category):
        r"""
        Returns the normalize morphism

        INPUT:

        - ``self`` -- a Jack basis of the symmetric functions
        - ``category`` -- a category

        OUTPUT:

        - the normalized morphism

        EXAMPLES::

            sage: JP = SymmetricFunctions(FractionField(QQ['t'])).jack().P()
            sage: normal = JP._normalize_morphism(AlgebrasWithBasis(JP.base_ring()))
            sage: normal.parent()
            Set of Homomorphisms from Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack P basis to Symmetric Functions over Fraction Field of Univariate Polynomial Ring in t over Rational Field in the Jack P basis
            sage: normal.category_for()
            Category of algebras with basis over Fraction Field of Univariate Polynomial Ring in t over Rational Field

            sage: t = JP.t
            sage: a = 2/(1/2*t+1/2)
            sage: b = 1/(1/3+1/6*t)
            sage: c = 24/(4*t^2 + 12*t + 8)
            sage: normal( a*JP[1] + b*JP[2] + c*JP[2,1] )
            (4/(t+1))*JackP[1] + (6/(t+2))*JackP[2] + (6/(t^2+3*t+2))*JackP[2, 1]

        .. TODO::

            This method should not be needed once short idioms to
            construct morphisms are available
        """
        return SetMorphism(End(self, category), self._normalize)

    def _multiply(self, left, right):
        r"""
        The product of two Jack symmetric functions is done by multiplying the
        elements in the `P` basis and then expressing the elements
        the basis ``self``.

        INPUT:

        - ``self`` -- a Jack basis of the symmetric functions
        - ``left``, ``right`` -- symmetric function elements

        OUTPUT:

        - returns the product of ``left`` and ``right`` expanded in the basis ``self``

        EXAMPLES::

            sage: JJ = SymmetricFunctions(FractionField(QQ['t'])).jack().J()
            sage: JJ([1])^2              # indirect doctest
            (t/(t+1))*JackJ[1, 1] + (1/(t+1))*JackJ[2]
            sage: JJ([2])^2
            (2*t^2/(2*t^2+3*t+1))*JackJ[2, 2] + (4*t/(3*t^2+4*t+1))*JackJ[3, 1] + ((t+1)/(6*t^2+5*t+1))*JackJ[4]
            sage: JQ = SymmetricFunctions(FractionField(QQ['t'])).jack().Q()
            sage: JQ([1])^2              # indirect doctest
            JackQ[1, 1] + (2/(t+1))*JackQ[2]
            sage: JQ([2])^2
            JackQ[2, 2] + (2/(t+1))*JackQ[3, 1] + ((6*t+6)/(6*t^2+5*t+1))*JackQ[4]
        """
        return self( self._P(left)*self._P(right) )

    def jack_family( self ):
        r"""
        Returns the family of Jack bases associated to the basis ``self``

        INPUT:

        - ``self`` -- a Jack basis of the symmetric functions

        OUTPUT:

        - the family of Jack symmetric functions associated to ``self``

        EXAMPLES::

            sage: JackP = SymmetricFunctions(QQ).jack(t=2).P()
            sage: JackP.jack_family()
            Jack polynomials with t=2 over Rational Field
        """
        return self._jack

    def coproduct_by_coercion(self, elt):
        r"""
        Returns the coproduct of the element ``elt`` by coercion to the Schur basis.

        INPUT:

        - ``self`` -- a Jack symmetric function basis
        - ``elt`` -- an instance of this basis

        OUTPUT:

        - The coproduct acting on ``elt``, the result is an element of the
          tensor squared of the Jack symmetric function basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ['t'].fraction_field())
            sage: Sym.jack().P()[2,2].coproduct() #indirect doctest
            JackP[] # JackP[2, 2] + (2/(t+1))*JackP[1] # JackP[2, 1] + ((8*t+4)/(t^3+4*t^2+5*t+2))*JackP[1, 1] # JackP[1, 1] + JackP[2] # JackP[2] + (2/(t+1))*JackP[2, 1] # JackP[1] + JackP[2, 2] # JackP[]
        """
        from sage.categories.tensor import tensor
        s = self.realization_of().schur()
        g = self.tensor_square().sum(coeff*tensor([self(s[x]), self(s[y])])
                                        for ((x,y), coeff) in s(elt).coproduct())
        normalize = self._normalize_coefficients
        return self.tensor_square().sum(normalize(coeff)*tensor([self(x), self(y)])
                    for ((x,y), coeff) in g)


    class Element(sfa.SymmetricFunctionAlgebra_generic.Element):
        def scalar_jack(self, x, t=None):
            r"""
            A scalar product where the power sums are orthogonal and
            `\langle p_\mu, p_\mu \rangle = z_\mu t^{length(\mu)}`

            INPUT:

            - ``self`` -- an element of a Jack basis of the symmetric functions
            - ``x`` -- an element of the symmetric functions
            - ``t`` -- an optional parameter (default : None uses the parameter from
                the basis)

            OUTPUT:

            - returns the Jack scalar product between ``x`` and ``self``

            EXAMPLES::

                sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
                sage: JP = Sym.jack().P()
                sage: JQ = Sym.jack().Q()
                sage: p = Partitions(3).list()
                sage: matrix([[JP(a).scalar_jack(JQ(b)) for a in p] for b in p])
                [1 0 0]
                [0 1 0]
                [0 0 1]
            """
            parent = self.parent()
            p = parent.realization_of().power()
            res = p(self).scalar_jack(p(x), t)

            return parent._normalize_coefficients(res)


def part_scalar_jack(part1, part2, t):
    r"""
    Returns the Jack scalar product between ``p(part1)`` and ``p(part2)`` where
    `p` is the power-sum basis.

    INPUT:

    - ``part1``, ``part2`` -- two partitions
    - ``t`` -- a parameter

    OUTPUT:

    - returns the scalar product between the power sum indexed by ``part1`` and ``part2``

    EXAMPLES::

        sage: Q.<t> = QQ[]
        sage: from sage.combinat.sf.jack import part_scalar_jack
        sage: matrix([[part_scalar_jack(p1,p2,t) for p1 in Partitions(4)] for p2 in Partitions(4)])
        [   4*t      0      0      0      0]
        [     0  3*t^2      0      0      0]
        [     0      0  8*t^2      0      0]
        [     0      0      0  4*t^3      0]
        [     0      0      0      0 24*t^4]
    """
    if part1 != part2:
        return 0
    else:
        return part1.centralizer_size()*t**len(part1)

#P basis

class JackPolynomials_p(JackPolynomials_generic):

    def __init__(self, jack):
        r"""
        The `P` basis is uni-triangularly related to the monomial basis and
        orthogonal with respect to the Jack scalar product.

        INPUT:

        - ``self`` -- an instance of the Jack `P` basis of the symmetric functions
        - ``jack`` -- a family of Jack symmetric function bases

        EXAMPLES::

            sage: P = SymmetricFunctions(FractionField(QQ['t'])).jack().P()
            sage: TestSuite(P).run(skip=['_test_associativity', '_test_distributivity', '_test_prod']) # products are too expensive
            sage: TestSuite(P).run(elements = [P.t*P[1,1]+P[2], P[1]+(1+P.t)*P[1,1]])
        """
        self._name = "Jack polynomials in the P basis"
        self._prefix = "JackP"

        self._m_to_self_cache = m_to_p_cache
        self._self_to_m_cache = p_to_m_cache
        JackPolynomials_generic.__init__(self, jack)

    def _m_cache(self, n):
        r"""
        Computes the change of basis between the Jack polynomials in the `P`
        basis and the monomial symmetric functions. This uses Gram-Schmidt
        to go to the monomials, and then that matrix is simply inverted.

        INPUT:

        - ``self`` -- an instance of the Jack `P` basis of the symmetric functions
        - ``n`` -- a positive integer indicating the degree

        EXAMPLES::

            sage: JP = SymmetricFunctions(FractionField(QQ['t'])).jack().P()
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: JP._m_cache(2)
            sage: l(JP._self_to_m_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], 2/(t + 1)), ([2], 1)])]
            sage: l(JP._m_to_self_cache[2])
            [([1, 1], [([1, 1], 1)]), ([2], [([1, 1], -2/(t + 1)), ([2], 1)])]
            sage: JP._m_cache(3)
            sage: l(JP._m_to_self_cache[3])
            [([1, 1, 1], [([1, 1, 1], 1)]),
             ([2, 1], [([1, 1, 1], -6/(t + 2)), ([2, 1], 1)]),
             ([3], [([1, 1, 1], 6/(t^2 + 3*t + 2)), ([2, 1], -3/(2*t + 1)), ([3], 1)])]
            sage: l(JP._self_to_m_cache[3])
            [([1, 1, 1], [([1, 1, 1], 1)]),
             ([2, 1], [([1, 1, 1], 6/(t + 2)), ([2, 1], 1)]),
             ([3], [([1, 1, 1], 6/(2*t^2 + 3*t + 1)), ([2, 1], 3/(2*t + 1)), ([3], 1)])]
        """
        if n in self._self_to_m_cache:
            return
        else:
            self._self_to_m_cache[n] = {}
        t = QQt.gen()
        monomial = sage.combinat.sf.sf.SymmetricFunctions(QQt).monomial()
        JP = sage.combinat.sf.sf.SymmetricFunctions(QQt).jack().P()
        JP._gram_schmidt(n, monomial, lambda p: part_scalar_jack(p,p,t), \
                           self._self_to_m_cache[n], upper_triangular=True)
        JP._invert_morphism(n, QQt, self._self_to_m_cache, \
                              self._m_to_self_cache, to_other_function = self._to_m)

    def _to_m(self, part):
        r"""
        Return a function that takes in a partition lambda that returns the
        coefficient of lambda in the expansion of self(part) in the
        monomial basis.

        This assumes that the cache from the Jack polynomials in the `P`
        basis to the monomial symmetric functions has already been
        computed.

        INPUT:

        - ``self`` -- an instance of the Jack `P` basis of the symmetric functions
        - ``part`` -- a partition

        OUTPUT:

        - returns a function that accepts a partition and returns the coefficients
          of the expansion of the element of ``P(part)`` in the monomial basis

        EXAMPLES::

            sage: JP = SymmetricFunctions(FractionField(QQ['t'])).jack().P()
            sage: JP._m_cache(3)
            sage: f = JP._to_m(Partition([2,1]))
            sage: [f(part) for part in Partitions(3)]
            [0, 1, 6/(t + 2)]
            sage: JP.symmetric_function_ring().m()(JP[2,1])
            (6/(t+2))*m[1, 1, 1] + m[2, 1]
        """
        f = lambda part2: self._self_to_m_cache[sum(part)][part].get(part2, 0)
        return f

    def _multiply(self, left, right):
        r"""
        The product of two Jack symmetric functions is done by multiplying the
        elements in the monomial basis and then expressing the elements
        the basis ``self``.

        INPUT:

        - ``self`` -- a Jack basis of the symmetric functions
        - ``left``, ``right`` -- symmetric function elements

        OUTPUT:

        - returns the product of ``left`` and ``right`` expanded in the basis ``self``

        EXAMPLES::

            sage: JP = SymmetricFunctions(FractionField(QQ['t'])).jack().P()
            sage: m = JP.symmetric_function_ring().m()
            sage: JP([1])^2 # indirect doctest
            (2*t/(t+1))*JackP[1, 1] + JackP[2]
            sage: m(_)
            2*m[1, 1] + m[2]
            sage: JP = SymmetricFunctions(QQ).jack(t=2).P()
            sage: JP([2,1])^2
            125/63*JackP[2, 2, 1, 1] + 25/12*JackP[2, 2, 2] + 25/18*JackP[3, 1, 1, 1] + 12/5*JackP[3, 2, 1] + 4/3*JackP[3, 3] + 4/3*JackP[4, 1, 1] + JackP[4, 2]
            sage: m(_)
            45*m[1, 1, 1, 1, 1, 1] + 51/2*m[2, 1, 1, 1, 1] + 29/2*m[2, 2, 1, 1] + 33/4*m[2, 2, 2] + 9*m[3, 1, 1, 1] + 5*m[3, 2, 1] + 2*m[3, 3] + 2*m[4, 1, 1] + m[4, 2]
        """
        return self( self._m(left)*self._m(right) )


    def scalar_jack_basis(self, part1, part2 = None):
        r"""
        Returns the scalar product of `P(part1)` and `P(part2)`.

        This is equation (10.16) of [Mc1995]_ on page 380.

        INPUT:

        - ``self`` -- an instance of the Jack `P` basis of the symmetric functions
        - ``part1`` -- a partition
        - ``part2`` -- an optional partition (default : None)

        OUTPUT:

        - the scalar product between `P(part1)` and `P(part2)` (or itself if `part2` is None)

        REFRENCES:

            .. [Mc1995] I. G. Macdonald, Symmetric functions and Hall polynomials, second ed.,
               The Clarendon Press, Oxford University Press, New York, 1995, With contributions
               by A. Zelevinsky, Oxford Science Publications.

        EXAMPLES::

            sage: JP = SymmetricFunctions(FractionField(QQ['t'])).jack().P()
            sage: JJ = SymmetricFunctions(FractionField(QQ['t'])).jack().J()
            sage: JP.scalar_jack_basis(Partition([2,1]), Partition([1,1,1]))
            0
            sage: JP._normalize_coefficients(JP.scalar_jack_basis(Partition([3,2,1]), Partition([3,2,1])))
            (12*t^6 + 20*t^5 + 11*t^4 + 2*t^3)/(2*t^3 + 11*t^2 + 20*t + 12)
            sage: JJ(JP[3,2,1]).scalar_jack(JP[3,2,1])
            (12*t^6 + 20*t^5 + 11*t^4 + 2*t^3)/(2*t^3 + 11*t^2 + 20*t + 12)

        With a single argument, takes `part2 = part1`::

            sage: JP.scalar_jack_basis(Partition([2,1]), Partition([2,1]))
            (2*t^3 + t^2)/(t + 2)
            sage: JJ(JP[2,1]).scalar_jack(JP[2,1])
            (2*t^3 + t^2)/(t + 2)
        """
        if part2 is not None and part1 != part2:
            return self.base_ring().zero()
        return self.c2(part1) / self.c1(part1)


    class Element(JackPolynomials_generic.Element):
        def scalar_jack(self, x, t=None):
            r"""
            The scalar product on the symmetric functions where the power sums
            are orthogonal and `\langle p_\mu, p_\mu \rangle = z_\mu t^{length(mu)}`
            where the t parameter from the Jack symmetric function family.

            INPUT:

            - ``self`` -- an element of the Jack `P` basis
            - ``x`` -- an element of the `P` basis

            EXAMPLES ::

                sage: JP = SymmetricFunctions(FractionField(QQ['t'])).jack().P()
                sage: l = [JP(p) for p in Partitions(3)]
                sage: matrix([[a.scalar_jack(b) for a in l] for b in l])
                [  6*t^3/(2*t^2 + 3*t + 1)                         0                         0]
                [                        0     (2*t^3 + t^2)/(t + 2)                         0]
                [                        0                         0 1/6*t^3 + 1/2*t^2 + 1/3*t]
            """
            if isinstance(x, JackPolynomials_p) and t is None:
                P = self.parent()
                return P._apply_multi_module_morphism(self, x, P.scalar_jack_basis, orthogonal=True)
            else:
                return JackPolynomials_generic.Element.scalar_jack(self, x, t)

#J basis

class JackPolynomials_j(JackPolynomials_generic):

    def __init__(self, jack):
        r"""
        The `J` basis is a defined as a normalized form of the `P` basis

        INPUT:

        - ``self`` -- an instance of the Jack `P` basis of the symmetric functions
        - ``jack`` -- a family of Jack symmetric function bases

        EXAMPLES::

            sage: J = SymmetricFunctions(FractionField(QQ['t'])).jack().J()
            sage: TestSuite(J).run(skip=['_test_associativity', '_test_distributivity', '_test_prod']) # products are too expensive
            sage: TestSuite(J).run(elements = [J.t*J[1,1]+J[2], J[1]+(1+J.t)*J[1,1]])  # long time (3s on sage.math, 2012)
        """
        self._name = "Jack polynomials in the J basis"
        self._prefix = "JackJ"
        JackPolynomials_generic.__init__(self, jack)

        # Should be shared with _q (and possibly other bases in Macdo/HL) as BasesByRenormalization
        self._P = self._jack.P()
        # temporary until Hom(GradedHopfAlgebrasWithBasis) works better
        category = sage.categories.all.ModulesWithBasis(self.base_ring())
        phi = self.module_morphism(diagonal = self.c1, codomain = self._P, category = category)
        # should use module_morphism(on_coeffs = ...) once it exists
        self._P.register_coercion(self._P._normalize_morphism(category) * phi)
        self   .register_coercion(self   ._normalize_morphism(category) *~phi)

    class Element(JackPolynomials_generic.Element):
        pass



#Q basis
class JackPolynomials_q(JackPolynomials_generic):

    def __init__(self, jack):
        r"""
        The `Q` basis is defined as a normalized form of the `P` basis

        INPUT:

        - ``self`` -- an instance of the Jack `Q` basis of the symmetric functions
        - ``jack`` -- a family of Jack symmetric function bases

        EXAMPLES::

            sage: Q = SymmetricFunctions(FractionField(QQ['t'])).jack().Q()
            sage: TestSuite(Q).run(skip=['_test_associativity', '_test_distributivity', '_test_prod']) # products are too expensive
            sage: TestSuite(Q).run(elements = [Q.t*Q[1,1]+Q[2], Q[1]+(1+Q.t)*Q[1,1]])  # long time (3s on sage.math, 2012)
        """
        self._name = "Jack polynomials in the Q basis"
        self._prefix = "JackQ"
        JackPolynomials_generic.__init__(self, jack)

        # Should be shared with _j (and possibly other bases in Macdo/HL) as BasesByRenormalization
        self._P = self._jack.P()
        # temporary until Hom(GradedHopfAlgebrasWithBasis) works better
        category = sage.categories.all.ModulesWithBasis(self.base_ring())
        phi = self._P.module_morphism(diagonal = self._P.scalar_jack_basis, codomain = self, category = category)
        self   .register_coercion(self   ._normalize_morphism(category) *  phi)
        self._P.register_coercion(self._P._normalize_morphism(category) * ~phi)

    class Element(JackPolynomials_generic.Element):
        pass

qp_to_h_cache = {}
h_to_qp_cache = {}
class JackPolynomials_qp(JackPolynomials_generic):
    def __init__(self, jack):
        r"""
        The `Qp` basis is the dual basis to the `P` basis with respect to the
        standard scalar product

        INPUT:

        - ``self`` -- an instance of the Jack `Qp` basis of the symmetric functions
        - ``jack`` -- a family of Jack symmetric function bases

        EXAMPLES::

            sage: Qp = SymmetricFunctions(FractionField(QQ['t'])).jack().Qp()
            sage: TestSuite(Qp).run(skip=['_test_associativity', '_test_distributivity', '_test_prod']) # products are too expensive
            sage: TestSuite(Qp).run(elements = [Qp.t*Qp[1,1]+Qp[2], Qp[1]+(1+Qp.t)*Qp[1,1]])  # long time (3s on sage.math, 2012)
        """
        self._name = "Jack polynomials in the Qp basis"
        self._prefix = "JackQp"
        JackPolynomials_generic.__init__(self, jack)
        self._P = self._jack.P()
        self._self_to_h_cache = qp_to_h_cache
        self._h_to_self_cache = h_to_qp_cache

    def _multiply(self, left, right):
        r"""
        The product of two Jack symmetric functions is done by multiplying the
        elements in the monomial basis and then expressing the elements
        the basis ``self``.

        INPUT:

        - ``self`` -- an instance of the Jack `Qp` basis of the symmetric functions
        - ``left``, ``right`` -- symmetric function elements

        OUTPUT:

        - returns the product of ``left`` and ``right`` expanded in the basis ``self``

        EXAMPLES::

            sage: JQp = SymmetricFunctions(FractionField(QQ['t'])).jack().Qp()
            sage: h = JQp.symmetric_function_ring().h()
            sage: JQp([1])^2 # indirect doctest
            JackQp[1, 1] + (2/(t+1))*JackQp[2]
            sage: h(_)
            h[1, 1]
            sage: JQp = SymmetricFunctions(QQ).jack(t=2).Qp()
            sage: h = SymmetricFunctions(QQ).h()
            sage: JQp([2,1])^2
            JackQp[2, 2, 1, 1] + 2/3*JackQp[2, 2, 2] + 2/3*JackQp[3, 1, 1, 1] + 48/35*JackQp[3, 2, 1] + 28/75*JackQp[3, 3] + 128/225*JackQp[4, 1, 1] + 28/75*JackQp[4, 2]
            sage: h(_)
            h[2, 2, 1, 1] - 6/5*h[3, 2, 1] + 9/25*h[3, 3]
        """
        return self( self._h(left)*self._h(right) )

    def _h_cache(self, n):
        r"""
        Computes the change of basis between the Jack polynomials in the `Qp`
        basis and the homogeneous symmetric functions. This uses the coefficients
        in the change of basis between the Jack `P` basis and the monomial basis.

        INPUT:

        - ``self`` -- an instance of the Jack `Qp` basis of the symmetric functions
        - ``n`` -- a positive integer indicating the degree

        EXAMPLES::

            sage: JQp = SymmetricFunctions(FractionField(QQ['t'])).jack().Qp()
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: JQp._h_cache(2)
            sage: l(JQp._self_to_h_cache[2])
            [([1, 1], [([1, 1], 1), ([2], -2/(t + 1))]), ([2], [([2], 1)])]
            sage: l(JQp._h_to_self_cache[2])
            [([1, 1], [([1, 1], 1), ([2], 2/(t + 1))]), ([2], [([2], 1)])]
            sage: JQp._h_cache(3)
            sage: l(JQp._h_to_self_cache[3])
            [([1, 1, 1], [([1, 1, 1], 1), ([2, 1], 6/(t + 2)), ([3], 6/(2*t^2 + 3*t + 1))]), ([2, 1], [([2, 1], 1), ([3], 3/(2*t + 1))]), ([3], [([3], 1)])]
            sage: l(JQp._self_to_h_cache[3])
            [([1, 1, 1], [([1, 1, 1], 1), ([2, 1], -6/(t + 2)), ([3], 6/(t^2 + 3*t + 2))]), ([2, 1], [([2, 1], 1), ([3], -3/(2*t + 1))]), ([3], [([3], 1)])]
        """
        if n in self._self_to_h_cache:
            return
        else:
            self._self_to_h_cache[n] = {}
            self._h_to_self_cache[n] = {}
        self._P._m_cache(n)
        from_cache_1 = self._P._self_to_m_cache[n]
        to_cache_1 = self._self_to_h_cache[n]
        from_cache_2 = self._P._m_to_self_cache[n]
        to_cache_2 = self._h_to_self_cache[n]
        for mu in from_cache_1.keys():
            for la in from_cache_1[mu].keys():
                if not la in to_cache_1:
                    to_cache_1[la] = {}
                    to_cache_2[la] = {}
                to_cache_2[la][mu] = from_cache_1[mu][la]
                to_cache_1[la][mu] = from_cache_2[mu][la]

    def _self_to_h( self, x ):
        r"""
        Isomorphism from self to the homogeneous basis

        INPUT:

        - ``self`` -- a Jack `Qp` basis of the symmetric functions
        - ``x`` -- an element of the Jack `Qp` basis

        OUTPUT:

        - an element of the homogeneous basis equivalent to ``x``

        EXAMPLES ::

            sage: Sym = SymmetricFunctions(QQ)
            sage: JQp = Sym.jack(t=2).Qp()
            sage: h = Sym.homogeneous()
            sage: JQp._self_to_h(JQp[2,1])
            h[2, 1] - 3/5*h[3]

        This is for internal use only. Please use instead::

            sage: h(JQp[2,1])
            h[2, 1] - 3/5*h[3]
        """
        return self._h._from_cache(x, self._h_cache, self._self_to_h_cache, t = self.t)

    def _h_to_self( self, x ):
        r"""
        Isomorphism from the homogeneous basis into ``self``

        INPUT:

        - ``self`` -- a Jack `Qp` basis of the symmetric functions
        - ``x`` -- element of the homogeneous basis

        OUTPUT:

        - an element of the Jack `Qp` basis equivalent to ``x``

        EXAMPLES ::

            sage: Sym = SymmetricFunctions(QQ)
            sage: JQp = Sym.jack(t=2).Qp()
            sage: h = Sym.homogeneous()
            sage: JQp._h_to_self(h[2,1])
            JackQp[2, 1] + 3/5*JackQp[3]

        This is for internal use only. Please use instead::

            sage: JQp(h[2,1])
            JackQp[2, 1] + 3/5*JackQp[3]
        """
        return self._from_cache(x, self._h_cache, self._h_to_self_cache, t = self.t)

    def coproduct_by_coercion( self, elt ):
        r"""
        Returns the coproduct of the element ``elt`` by coercion to the Schur basis.

        INPUT:

        - ``elt`` -- an instance of the ``Qp`` basis

        OUTPUT:

        - The coproduct acting on ``elt``, the result is an element of the
          tensor squared of the ``Qp`` symmetric function basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ['t'].fraction_field())
            sage: JQp = Sym.jack().Qp()
            sage: JQp[2,2].coproduct()   #indirect doctest
            JackQp[] # JackQp[2, 2] + (2*t/(t+1))*JackQp[1] # JackQp[2, 1] + JackQp[1, 1] # JackQp[1, 1] + ((4*t^3+8*t^2)/(2*t^3+5*t^2+4*t+1))*JackQp[2] # JackQp[2] + (2*t/(t+1))*JackQp[2, 1] # JackQp[1] + JackQp[2, 2] # JackQp[]
        """
        h = elt.parent().realization_of().h()
        parent = elt.parent()
        from sage.categories.tensor import tensor
        cfunc = lambda x, y: tensor([parent(x), parent(y)])
        cprod = h(elt).coproduct().apply_multilinear_morphism( cfunc )
        normalize = lambda c: normalize_coefficients( parent, c )
        return cprod.parent().sum(normalize(coeff)*tensor([parent(x), parent(y)])
                        for ((x,y), coeff) in cprod)

    class Element(JackPolynomials_generic.Element):
        pass

#Zonal polynomials ( =P(at t=2) )
class SymmetricFunctionAlgebra_zonal(sfa.SymmetricFunctionAlgebra_generic):
    def __init__(self, Sym):
        r"""
        Returns the algebra of zonal polynomials.

        INPUT:

        - ``self`` -- a zonal basis of the symmetric functions
        - ``Sym`` -- a ring of the symmetric functions

        EXAMPLES ::

            sage: Z = SymmetricFunctions(QQ).zonal()
            sage: Z([2])^2
            64/45*Z[2, 2] + 16/21*Z[3, 1] + Z[4]
            sage: Z = SymmetricFunctions(QQ).zonal()
            sage: TestSuite(Z).run(skip=['_test_associativity', '_test_distributivity', '_test_prod']) # products are too expensive
            sage: TestSuite(Z).run(elements = [Z[1,1]+Z[2], Z[1]+2*Z[1,1]])
        """
        self._sym = Sym
        self._jack = self._sym.jack(t=2)
        self._P = self._jack.P()
        #self._m_to_self_cache = {} Now that we compute Jacks once, there is a global cache
        #self._self_to_m_cache = {} and we don't need to compute it separately for zonals
        sfa.SymmetricFunctionAlgebra_generic.__init__(self, self._sym,
                                                      prefix="Z", basis_name="zonal")
        category = sage.categories.all.ModulesWithBasis(self._sym.base_ring())
        self   .register_coercion(SetMorphism(Hom(self._P, self, category), self.sum_of_terms))
        self._P.register_coercion(SetMorphism(Hom(self, self._P, category), self._P.sum_of_terms))

    def _multiply(self, left, right):
        r"""
        The product of two zonal symmetric functions is done by multiplying the
        elements in the monomial basis and then expressing the elements
        the basis ``self``.

        INPUT:

        - ``self`` -- a zonal basis of the symmetric functions
        - ``left``, ``right`` -- symmetric function elements

        OUTPUT:

        - returns the product of ``left`` and ``right`` expanded in the basis ``self``

        EXAMPLES ::

            sage: Sym = SymmetricFunctions(QQ)
            sage: Z = Sym.zonal()
            sage: JP = Sym.jack(t=1).P()
            sage: Z([2])*Z([3])                    # indirect doctest
            192/175*Z[3, 2] + 32/45*Z[4, 1] + Z[5]
            sage: Z([2])*JP([2])
            10/27*Z[2, 1, 1] + 64/45*Z[2, 2] + 23/21*Z[3, 1] + Z[4]
            sage: JP = Sym.jack(t=2).P()
            sage: Z([2])*JP([2])
            64/45*Z[2, 2] + 16/21*Z[3, 1] + Z[4]
        """
        return self( self._P(left)*self._P(right) )


    class Element(sfa.SymmetricFunctionAlgebra_generic.Element):
        def scalar_zonal(self, x):
            r"""
            The zonal scalar product has the power sum basis and the zonal
            symmetric functions are orthogonal. In particular,
            `\langle p_\mu, p_\mu \rangle = z_\mu 2^{length(\mu)}`.

            INPUT:

            - ``self`` -- an element of the zonal basis
            - ``x`` -- an element of the symmetric function

            OUTPUT:

            - the scalar product between ``self`` and ``x``

            EXAMPLES ::

                sage: Sym = SymmetricFunctions(QQ)
                sage: Z = Sym.zonal()
                sage: parts = Partitions(3).list()
                sage: matrix([[Z(a).scalar_zonal(Z(b)) for a in parts] for b in parts])
                [16/5    0    0]
                [   0    5    0]
                [   0    0    4]
                sage: p = Z.symmetric_function_ring().power()
                sage: matrix([[Z(p(a)).scalar_zonal(p(b)) for a in parts] for b in parts])
                [ 6  0  0]
                [ 0  8  0]
                [ 0  0 48]
            """
            P = self.parent()._P
            return P(self).scalar_jack(P(x),2)

#############
#   Cache   #
#############
#from sage.misc.cache import Cache
#cache_p = Cache(JackPolynomials_p)
#cache_j = Cache(JackPolynomials_j)
#cache_q = Cache(JackPolynomials_q)
#cache_z = Cache(SymmetricFunctionAlgebra_zonal)

# Backward compatibility for unpickling
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.sf.jack', 'JackPolynomial_qp', JackPolynomials_qp.Element)
register_unpickle_override('sage.combinat.sf.jack', 'JackPolynomial_j', JackPolynomials_j.Element)
register_unpickle_override('sage.combinat.sf.jack', 'JackPolynomial_p', JackPolynomials_p.Element)
register_unpickle_override('sage.combinat.sf.jack', 'JackPolynomial_q', JackPolynomials_q.Element)
#register_unpickle_override('sage.combinat.sf.jack', 'SymmetricFunctionAlgebra_zonal',  SymmetricFunctionAlgebra_zonal.Element)
