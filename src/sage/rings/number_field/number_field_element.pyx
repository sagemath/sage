"""
Number Field Elements

AUTHORS:

- William Stein: version before it got Cython'd

- Joel B. Mohler (2007-03-09): First reimplementation in Cython

- William Stein (2007-09-04): add doctests

- Robert Bradshaw (2007-09-15): specialized classes for relative and
  absolute elements

- John Cremona (2009-05-15): added support for local and global
  logarithmic heights.

- Robert Harron (2012-08): conjugate() now works for all fields contained in
  CM fields

"""
#*****************************************************************************
#       Copyright (C) 2004, 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator

include 'sage/ext/interrupt.pxi'
from cpython.int cimport *
include "sage/ext/stdsage.pxi"

import sage.rings.field_element
import sage.rings.infinity
import sage.rings.polynomial.polynomial_element
import sage.rings.rational_field
import sage.rings.rational
import sage.rings.integer_ring
import sage.rings.integer
import sage.rings.arith

import number_field

from sage.rings.integer_ring cimport IntegerRing_class
from sage.rings.rational cimport Rational

from sage.modules.free_module_element import vector

from sage.libs.pari.all import pari_gen
from sage.structure.element cimport Element, generic_power_c
from sage.structure.element import canonical_coercion, parent

QQ = sage.rings.rational_field.QQ
ZZ = sage.rings.integer_ring.ZZ
Integer_sage = sage.rings.integer.Integer

from sage.rings.real_mpfi import RealInterval

from sage.rings.complex_field import ComplexField
CC = ComplexField(53)

# this is a threshold for the charpoly() methods in this file
# for degrees <= this threshold, pari is used
# for degrees > this threshold, sage matrices are used
# the value was decided by running a tuning script on a number of
# architectures; you can find this script attached to trac
# ticket 5213
TUNE_CHARPOLY_NF = 25

def is_NumberFieldElement(x):
    """
    Return True if x is of type NumberFieldElement, i.e., an element of
    a number field.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field_element import is_NumberFieldElement
        sage: is_NumberFieldElement(2)
        False
        sage: k.<a> = NumberField(x^7 + 17*x + 1)
        sage: is_NumberFieldElement(a+1)
        True
    """
    return PY_TYPE_CHECK(x, NumberFieldElement)

def __create__NumberFieldElement_version0(parent, poly):
    """
    Used in unpickling elements of number fields pickled under very old Sage versions.

    EXAMPLE::

        sage: k.<a> = NumberField(x^3 - 2)
        sage: R.<z> = QQ[]
        sage: sage.rings.number_field.number_field_element.__create__NumberFieldElement_version0(k, z^2 + z + 1)
        a^2 + a + 1
    """
    return NumberFieldElement(parent, poly)

def __create__NumberFieldElement_version1(parent, cls, poly):
    """
    Used in unpickling elements of number fields.

    EXAMPLES:

    Since this is just used in unpickling, we unpickle.

    ::

        sage: k.<a> = NumberField(x^3 - 2)
        sage: loads(dumps(a+1)) == a + 1 # indirect doctest
        True

    This also gets called for unpickling order elements; we check that #6462 is
    fixed::

        sage: L = NumberField(x^3 - x - 1,'a'); OL = L.maximal_order(); w = OL.0
        sage: loads(dumps(w)) == w # indirect doctest
        True
    """
    return cls(parent, poly)

def _inverse_mod_generic(elt, I):
    r"""
    Return an inverse of elt modulo the given ideal. This is a separate
    function called from each of the OrderElement_xxx classes, since
    otherwise we'd have to have the same code three times over (there
    is no OrderElement_generic class - no multiple inheritance). See
    trac 4190.

    EXAMPLES::

        sage: OE = NumberField(x^3 - x + 2, 'w').ring_of_integers()
        sage: w = OE.ring_generators()[0]
        sage: from sage.rings.number_field.number_field_element import _inverse_mod_generic
        sage: _inverse_mod_generic(w, 13*OE)
        6*w^2 - 6
    """
    from sage.matrix.constructor import matrix
    R = elt.parent()
    try:
        I = R.ideal(I)
    except ValueError:
        raise ValueError, "inverse is only defined modulo integral ideals"
    if I == 0:
        raise ValueError, "inverse is not defined modulo the zero ideal"
    n = R.absolute_degree()
    m = matrix(ZZ, map(R.coordinates, I.integral_basis() + [elt*s for s in R.gens()]))
    a, b = m.echelon_form(transformation=True)
    if a[0:n] != 1:
        raise ZeroDivisionError, "%s is not invertible modulo %s" % (elt, I)
    v = R.coordinates(1)
    y = R(0)
    for j in xrange(n):
        if v[j] != 0:
            y += v[j] * sum([b[j,i+n] * R.gen(i) for i in xrange(n)])
    return I.small_residue(y)

__pynac_pow = False

cdef class NumberFieldElement(FieldElement):
    """
    An element of a number field.

    EXAMPLES::

        sage: k.<a> = NumberField(x^3 + x + 1)
        sage: a^3
        -a - 1
    """
    cdef _new(self):
        """
        Quickly creates a new initialized NumberFieldElement with the same
        parent as self.
        """
        cdef NumberFieldElement x
        x = <NumberFieldElement>PY_NEW_SAME_TYPE(self)
        x._parent = self._parent
        x.__fld_numerator = self.__fld_numerator
        x.__fld_denominator = self.__fld_denominator
        return x

    cdef number_field(self):
        r"""

        Return the number field of self. Only accessible from Cython.
        EXAMPLE::

            sage: K.<a> = NumberField(x^3 + 3)
            sage: a._number_field() # indirect doctest
            Number Field in a with defining polynomial x^3 + 3
        """
        return self._parent

    def _number_field(self):
        r"""
        EXAMPLE::

            sage: K.<a> = NumberField(x^3 + 3)
            sage: a._number_field()
            Number Field in a with defining polynomial x^3 + 3
        """
        return self.number_field()

    def __init__(self, parent, f):
        """
        INPUT:


        -  ``parent`` - a number field

        -  ``f`` - defines an element of a number field.


        EXAMPLES:

        The following examples illustrate creation of elements of
        number fields, and some basic arithmetic.

        First we define a polynomial over Q::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^2 + 1

        Next we use f to define the number field::

            sage: K.<a> = NumberField(f); K
            Number Field in a with defining polynomial x^2 + 1
            sage: a = K.gen()
            sage: a^2
            -1
            sage: (a+1)^2
            2*a
            sage: a^2
            -1
            sage: z = K(5); 1/z
            1/5

        We create a cube root of 2::

            sage: K.<b> = NumberField(x^3 - 2)
            sage: b = K.gen()
            sage: b^3
            2
            sage: (b^2 + b + 1)^3
            12*b^2 + 15*b + 19

        We can create number field elements from PARI::

            sage: K.<a> = NumberField(x^3 - 17)
            sage: K(pari(42))
            42
            sage: K(pari("5/3"))
            5/3
            sage: K(pari("[3/2, -5, 0]~"))    # Uses Z-basis
            -5/3*a^2 + 5/3*a - 1/6

        From a PARI polynomial or ``POLMOD``, note that the variable
        name does not matter::

            sage: K(pari("-5/3*q^2 + 5/3*q - 1/6"))
            -5/3*a^2 + 5/3*a - 1/6
            sage: K(pari("Mod(-5/3*q^2 + 5/3*q - 1/6, q^3 - 17)"))
            -5/3*a^2 + 5/3*a - 1/6
            sage: K(pari("x^5/17"))
            a^2
            sage: K(pari("Mod(-5/3*q^2 + 5/3*q - 1/6, q^3 - 999)"))  # Wrong modulus
            Traceback (most recent call last):
            ...
            TypeError: Coercion of PARI polmod with modulus q^3 - 999 into number field with defining polynomial x^3 - 17 failed

        This example illustrates save and load::

            sage: K.<a> = NumberField(x^17 - 2)
            sage: s = a^15 - 19*a + 3
            sage: loads(s.dumps()) == s
            True

        TESTS:

        Test round-trip conversion to PARI and back::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^3 - 1/2*x + 1/3)
            sage: b = K.random_element()
            sage: K(pari(b)) == b
            True
        """
        sage.rings.field_element.FieldElement.__init__(self, parent)
        self.__fld_numerator, self.__fld_denominator = parent.absolute_polynomial_ntl()

        cdef ZZ_c coeff
        if isinstance(f, (int, long, Integer_sage)):
            # set it up and exit immediately
            # fast pathway
            (<Integer>ZZ(f))._to_ZZ(&coeff)
            ZZX_SetCoeff( self.__numerator, 0, coeff )
            ZZ_conv_from_int( self.__denominator, 1 )
            return

        elif isinstance(f, NumberFieldElement):
            if PY_TYPE(self) is PY_TYPE(f):
                self.__numerator = (<NumberFieldElement>f).__numerator
                self.__denominator = (<NumberFieldElement>f).__denominator
                return
            else:
                f = f.polynomial()

        from sage.rings.number_field import number_field_rel
        if isinstance(parent, number_field_rel.NumberField_relative):
            ppr = parent.base_field().polynomial_ring()
        else:
            ppr = parent.polynomial_ring()

        cdef long i
        if isinstance(f, pari_gen):
            if f.type() in ["t_INT", "t_FRAC", "t_POL"]:
                pass
            elif f.type() == "t_POLMOD":
                # Check whether we are dealing with a *relative*
                # number field element
                if parent.is_relative():
                    # If the modulus is a polynomial with polynomial
                    # coefficients, then the element is relative.
                    fmod = f.mod()
                    for i from 0 <= i <= fmod.poldegree():
                        if fmod.polcoeff(i).type() in ["t_POL", "t_POLMOD"]:
                            # Convert relative element to absolute
                            # This returns a polynomial, not a polmod
                            f = parent.pari_rnf().rnfeltreltoabs(f)
                            break
                # Check that the modulus is actually the defining polynomial
                # of the number field.
                # Unfortunately, this check only works for absolute elements
                # since the rnfeltreltoabs() destroys all information about
                # the number field.
                if f.type() == "t_POLMOD":
                    fmod = f.mod()
                    if fmod != parent.pari_polynomial(fmod.variable()):
                        raise TypeError("Coercion of PARI polmod with modulus %s into number field with defining polynomial %s failed"%(fmod, parent.pari_polynomial()))
                    f = f.lift()
            else:
                f = parent.pari_nf().nfbasistoalg_lift(f)
        f = ppr(f)
        if f.degree() >= parent.absolute_degree():
            from sage.rings.number_field import number_field_rel
            if isinstance(parent, number_field_rel.NumberField_relative):
                f %= ppr(parent.absolute_polynomial())
            else:
                f %= parent.polynomial()

        # Set Denominator
        den = f.denominator()
        (<Integer>ZZ(den))._to_ZZ(&self.__denominator)

        num = f * den
        for i from 0 <= i <= num.degree():
            (<Integer>ZZ(num[i]))._to_ZZ(&coeff)
            ZZX_SetCoeff( self.__numerator, i, coeff )

    def __cinit__(self):
        r"""
        Initialise C variables.

        EXAMPLE::

            sage: K.<b> = NumberField(x^3 - 2)
            sage: b = K.gen(); b # indirect doctest
            b
        """
        ZZX_construct(&self.__numerator)
        ZZ_construct(&self.__denominator)

    def __dealloc__(self):
        ZZX_destruct(&self.__numerator)
        ZZ_destruct(&self.__denominator)

    def _lift_cyclotomic_element(self, new_parent, bint check=True, int rel=0):
        """
        Creates an element of the passed field from this field. This is
        specific to creating elements in a cyclotomic field from elements
        in another cyclotomic field, in the case that
        self.number_field()._n() divides new_parent()._n(). This
        function aims to make this common coercion extremely fast!

        More general coercion (i.e. of zeta6 into CyclotomicField(3)) is
        implemented in the _coerce_from_other_cyclotomic_field method
        of a CyclotomicField.

        EXAMPLES::

            sage: C.<zeta5>=CyclotomicField(5)
            sage: CyclotomicField(10)(zeta5+1)  # The function _lift_cyclotomic_element does the heavy lifting in the background
            zeta10^2 + 1
            sage: (zeta5+1)._lift_cyclotomic_element(CyclotomicField(10))  # There is rarely a purpose to call this function directly
            zeta10^2 + 1
            sage: cf4 = CyclotomicField(4)
            sage: cf1 = CyclotomicField(1) ; one = cf1.0
            sage: cf4(one)
            1
            sage: type(cf4(1))
            <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
            sage: cf33 = CyclotomicField(33) ; z33 = cf33.0
            sage: cf66 = CyclotomicField(66) ; z66 = cf66.0
            sage: z33._lift_cyclotomic_element(cf66)
            zeta66^2
            sage: z66._lift_cyclotomic_element(cf33)
            Traceback (most recent call last):
            ...
            TypeError: The zeta_order of the new field must be a multiple of the zeta_order of the original.
            sage: cf33(z66)
            -zeta33^17

        AUTHORS:

        - Joel B. Mohler

        - Craig Citro (fixed behavior for different representation of
          quadratic field elements)
        """
        if check:
            if not isinstance(self.number_field(), number_field.NumberField_cyclotomic) \
                   or not isinstance(new_parent, number_field.NumberField_cyclotomic):
                raise TypeError, "The field and the new parent field must both be cyclotomic fields."

        if rel == 0:
            small_order = self.number_field()._n()
            large_order = new_parent._n()

            try:
                rel = ZZ(large_order / small_order)
            except TypeError:
                raise TypeError, "The zeta_order of the new field must be a multiple of the zeta_order of the original."

        ## degree 2 is handled differently, because elements are
        ## represented differently
        if new_parent.degree() == 2:
            if rel == 1:
                return new_parent._element_class(new_parent, self)
            else:
                return self.polynomial()(new_parent.gen()**rel)

        cdef NumberFieldElement x = <NumberFieldElement>PY_NEW_SAME_TYPE(self)
        x._parent = <ParentWithBase>new_parent
        x.__fld_numerator, x.__fld_denominator = new_parent.polynomial_ntl()
        x.__denominator = self.__denominator
        cdef ZZX_c result
        cdef ZZ_c tmp
        cdef int i
        cdef ntl_ZZX _num
        cdef ntl_ZZ _den
        for i from 0 <= i <= ZZX_deg(self.__numerator):
            tmp = ZZX_coeff(self.__numerator, i)
            ZZX_SetCoeff(result, i*rel, tmp)
        ZZX_rem(x.__numerator, result, x.__fld_numerator.x)
        return x

    def __reduce__(self):
        """
        Used in pickling number field elements.

        Note for developers: If this is changed, please also change the doctests of __create__NumberFieldElement_version1.

        EXAMPLES::

            sage: k.<a> = NumberField(x^3 - 17*x^2 + 1)
            sage: t = a.__reduce__(); t
            (<built-in function __create__NumberFieldElement_version1>, (Number Field in a with defining polynomial x^3 - 17*x^2 + 1, <type 'sage.rings.number_field.number_field_element.NumberFieldElement_absolute'>, x))
            sage: t[0](*t[1]) == a
            True
        """
        return __create__NumberFieldElement_version1, \
               (self.parent(), type(self), self.polynomial())

    def _repr_(self):
        """
        String representation of this number field element, which is just a
        polynomial in the generator.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 2)
            sage: b = (2/3)*a + 3/5
            sage: b._repr_()
            '2/3*a + 3/5'
        """
        x = self.polynomial()
        K = self.number_field()
        return str(x).replace(x.parent().variable_name(), K.variable_name())

    def _im_gens_(self, codomain, im_gens):
        """
        This is used in computing homomorphisms between number fields.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 - 2)
            sage: m.<b> = NumberField(x^4 - 2)
            sage: phi = k.hom([b^2])
            sage: phi(a+1)
            b^2 + 1
            sage: (a+1)._im_gens_(m, [b^2])
            b^2 + 1
        """
        # NOTE -- if you ever want to change this so relative number
        # fields are in terms of a root of a poly.  The issue is that
        # elements of a relative number field are represented in terms
        # of a generator for the absolute field.  However the morphism
        # gives the image of gen, which need not be a generator for
        # the absolute field.  The morphism has to be *over* the
        # relative element.
        return codomain(self.polynomial()(im_gens[0]))

    def _latex_(self):
        """
        Returns the latex representation for this element.

        EXAMPLES::

            sage: C,zeta12=CyclotomicField(12).objgen()
            sage: latex(zeta12^4-zeta12) # indirect doctest
            \zeta_{12}^{2} - \zeta_{12} - 1
        """
        return self.polynomial()._latex_(name=self.number_field().latex_variable_name())

    def _gap_init_(self):
        """
        Return gap string representation of self.

        EXAMPLES::

            sage: F=CyclotomicField(8)
            sage: F.gen()
            zeta8
            sage: F._gap_init_()
            'CyclotomicField(8)'
            sage: f = gap(F)
            sage: f.GeneratorsOfDivisionRing()
            [ E(8) ]
            sage: p=F.gen()^2+2*F.gen()-3
            sage: p
            zeta8^2 + 2*zeta8 - 3
            sage: p._gap_init_() # The variable name $sage2 belongs to the gap(F) and is somehow random
            'GeneratorsOfField($sage2)[1]^2 + 2*GeneratorsOfField($sage2)[1] - 3'
            sage: gap(p._gap_init_())
            -3+2*E(8)+E(8)^2
        """
        s = self._repr_()
        return s.replace(str(self.parent().gen()), 'GeneratorsOfField(%s)[1]'%sage.interfaces.gap.gap(self.parent()).name())

    def _libgap_(self):
        """
        Return a LibGAP representation of ``self``.

        EXAMPLES::

            sage: F = CyclotomicField(8)
            sage: F.gen()._libgap_()
            E(8)
            sage: libgap(F.gen())   # syntactic sugar
            E(8)
            sage: E8 = F.gen()
            sage: libgap(E8 + 3/2*E8^2 + 100*E8^7)
            E(8)+3/2*E(8)^2-100*E(8)^3
            sage: type(_)
            <type 'sage.libs.gap.element.GapElement_Cyclotomic'>
        """
        n = self.parent()._n()
        from sage.libs.gap.libgap import libgap
        En = libgap(self.parent()).GeneratorsOfField()[0]
        return self.polynomial()(En)

    def _pari_(self, name='y'):
        r"""
        Return PARI representation of self.

        The returned element is a PARI ``POLMOD`` in the variable
        ``name``, which is by default 'y' - not the name of the generator
        of the number field.

        INPUT:

        -  ``name`` -- (default: 'y') the PARI variable name used.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + 2)
            sage: K(1)._pari_()
            Mod(1, y^3 + 2)
            sage: (a + 2)._pari_()
            Mod(y + 2, y^3 + 2)
            sage: L.<b> = K.extension(x^2 + 2)
            sage: (b + a)._pari_()
            Mod(24/101*y^5 - 9/101*y^4 + 160/101*y^3 - 156/101*y^2 + 397/101*y + 364/101, y^6 + 6*y^4 - 4*y^3 + 12*y^2 + 24*y + 12)

        ::

            sage: k.<j> = QuadraticField(-1)
            sage: j._pari_('j')
            Mod(j, j^2 + 1)
            sage: pari(j)
            Mod(y, y^2 + 1)

        By default the variable name is 'y'. This allows 'x' to be used
        as polynomial variable::

            sage: P.<a> = PolynomialRing(QQ)
            sage: K.<b> = NumberField(a^2 + 1)
            sage: R.<x> = PolynomialRing(K)
            sage: pari(b*x)
            Mod(y, y^2 + 1)*x

        In PARI many variable names are reserved, for example ``theta``
        and ``I``::

            sage: R.<theta> = PolynomialRing(QQ)
            sage: K.<theta> = NumberField(theta^2 + 1)
            sage: theta._pari_('theta')
            Traceback (most recent call last):
            ...
            PariError: theta already exists with incompatible valence
            sage: theta._pari_()
            Mod(y, y^2 + 1)
            sage: k.<I> = QuadraticField(-1)
            sage: I._pari_('I')
            Traceback (most recent call last):
            ...
            PariError: I already exists with incompatible valence

        Instead, request the variable be named different for the coercion::

            sage: pari(I)
            Mod(y, y^2 + 1)
            sage: I._pari_('i')
            Mod(i, i^2 + 1)
            sage: I._pari_('II')
            Mod(II, II^2 + 1)

        Examples with relative number fields, which always yield an
        *absolute* representation of the element::

            sage: y = QQ['y'].gen()
            sage: k.<j> = NumberField([y^2 - 7, y^3 - 2])
            sage: pari(j)
            Mod(42/5515*y^5 - 9/11030*y^4 - 196/1103*y^3 + 273/5515*y^2 + 10281/5515*y + 4459/11030, y^6 - 21*y^4 + 4*y^3 + 147*y^2 + 84*y - 339)
            sage: j^2
            7
            sage: pari(j)^2
            Mod(7, y^6 - 21*y^4 + 4*y^3 + 147*y^2 + 84*y - 339)
            sage: (j^2)._pari_('x')
            Mod(7, x^6 - 21*x^4 + 4*x^3 + 147*x^2 + 84*x - 339)

        A tower of three number fields::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2 + 2)
            sage: L.<b> = NumberField(polygen(K)^2 + a)
            sage: M.<c> = NumberField(polygen(L)^3 + b)
            sage: L(b)._pari_()
            Mod(y, y^4 + 2)
            sage: M(b)._pari_('c')
            Mod(-c^3, c^12 + 2)
            sage: c._pari_('c')
            Mod(c, c^12 + 2)
        """
        try:
            return self.__pari[name]
        except KeyError:
            pass
        except TypeError:
            self.__pari = {}
        f = self.polynomial()._pari_or_constant(name)
        g = self.number_field().pari_polynomial(name)
        h = f.Mod(g)
        self.__pari[name] = h
        return h

    def _pari_init_(self, name='y'):
        """
        Return PARI/GP string representation of self.

        The returned string defines a PARI ``POLMOD`` in the variable
        ``name``, which is by default 'y' - not the name of the generator
        of the number field.

        INPUT:

        -  ``name`` -- (default: 'y') the PARI variable name used.

        EXAMPLES::

            sage: K.<a> = NumberField(x^5 - x - 1)
            sage: ((1 + 1/3*a)^4)._pari_init_()
            'Mod(1/81*y^4 + 4/27*y^3 + 2/3*y^2 + 4/3*y + 1, y^5 - y - 1)'
            sage: ((1 + 1/3*a)^4)._pari_init_('a')
            'Mod(1/81*a^4 + 4/27*a^3 + 2/3*a^2 + 4/3*a + 1, a^5 - a - 1)'

        Note that _pari_init_ can fail because of reserved words in
        PARI, and since it actually works by obtaining the PARI
        representation of something::

            sage: K.<theta> = NumberField(x^5 - x - 1)
            sage: b = (1/2 - 2/3*theta)^3; b
            -8/27*theta^3 + 2/3*theta^2 - 1/2*theta + 1/8
            sage: b._pari_init_('theta')
            Traceback (most recent call last):
            ...
            PariError: theta already exists with incompatible valence

        Fortunately pari_init returns everything in terms of y by
        default::

            sage: pari(b)
            Mod(-8/27*y^3 + 2/3*y^2 - 1/2*y + 1/8, y^5 - y - 1)
        """
        return repr(self._pari_(name=name))

    def __getitem__(self, n):
        """
        Return the n-th coefficient of this number field element, written
        as a polynomial in the generator.

        Note that `n` must be between 0 and `d-1`, where
        `d` is the degree of the number field.

        EXAMPLES::

            sage: m.<b> = NumberField(x^4 - 1789)
            sage: c = (2/3-4/5*b)^3; c
            -64/125*b^3 + 32/25*b^2 - 16/15*b + 8/27
            sage: c[0]
            8/27
            sage: c[2]
            32/25
            sage: c[3]
            -64/125

        We illustrate bounds checking::

            sage: c[-1]
            Traceback (most recent call last):
            ...
            IndexError: index must be between 0 and degree minus 1.
            sage: c[4]
            Traceback (most recent call last):
            ...
            IndexError: index must be between 0 and degree minus 1.

        The list method implicitly calls ``__getitem__``::

            sage: list(c)
            [8/27, -16/15, 32/25, -64/125]
            sage: m(list(c)) == c
            True
        """
        if n < 0 or n >= self.number_field().degree():     # make this faster.
            raise IndexError, "index must be between 0 and degree minus 1."
        return self.polynomial()[n]

    def __richcmp__(left, right, int op):
        r"""
        EXAMPLE::

            sage: K.<a> = NumberField(x^3 - 3*x + 8)
            sage: a  + 1 > a # indirect doctest
            True
            sage: a + 1 < a # indirect doctest
            False
        """
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, sage.structure.element.Element right) except -2:
        cdef NumberFieldElement _right = right
        return not (ZZX_equal(left.__numerator, _right.__numerator) and ZZ_equal(left.__denominator, _right.__denominator))

    def _random_element(self, num_bound=None, den_bound=None, distribution=None):
        """
        Return a new random element with the same parent as self.

        INPUT:

        - ``num_bound`` - Bound for the numerator of coefficients of result

        - ``den_bound`` - Bound for the denominator of coefficients of result

        - ``distribution`` - Distribution to use for coefficients of result

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-2)
            sage: a._random_element()
            -1/2*a^2 - 4
            sage: K.<a> = NumberField(x^2-5)
            sage: a._random_element()
            -2*a - 1
        """
        cdef NumberFieldElement elt = self._new()
        elt._randomize(num_bound, den_bound, distribution)
        return elt

    cdef int _randomize(self, num_bound, den_bound, distribution) except -1:
        cdef int i
        cdef Integer denom_temp = PY_NEW(Integer)
        cdef Integer tmp_integer = PY_NEW(Integer)
        cdef ZZ_c ntl_temp
        cdef list coeff_list
        cdef Rational tmp_rational

        # It seems like a simpler approach would be to simply generate
        # random integers for each coefficient of self.__numerator
        # and an integer for self.__denominator. However, this would
        # generate things with a fairly fixed shape: in particular,
        # we'd be very unlikely to get elements like 1/3*a^3 + 1/7,
        # or anything where the denominators are actually unrelated
        # to one another. The extra code below is to make exactly
        # these kinds of results possible.

        if den_bound == 1:
            # in this case, we can skip all the business with LCMs,
            # storing a list of rationals, etc. this gives a factor of
            # two or so speedup ...

            # set the denominator
            mpz_set_si(denom_temp.value, 1)
            denom_temp._to_ZZ(&self.__denominator)
            for i from 0 <= i < ZZX_deg(self.__fld_numerator.x):
                tmp_integer = <Integer>(ZZ.random_element(x=num_bound,
                                                   distribution=distribution))
                tmp_integer._to_ZZ(&ntl_temp)
                ZZX_SetCoeff(self.__numerator, i, ntl_temp)

        else:
            coeff_list = []
            mpz_set_si(denom_temp.value, 1)
            tmp_integer = PY_NEW(Integer)

            for i from 0 <= i < ZZX_deg(self.__fld_numerator.x):
                tmp_rational = <Rational>(QQ.random_element(num_bound=num_bound,
                                                            den_bound=den_bound,
                                                            distribution=distribution))
                coeff_list.append(tmp_rational)
                mpz_lcm(denom_temp.value, denom_temp.value,
                        mpq_denref(tmp_rational.value))

            # now denom_temp has the denominator, and we just need to
            # scale the numerators and set everything appropriately

            # first, the denominator (easy)
            denom_temp._to_ZZ(&self.__denominator)

            # now the coefficients themselves.
            for i from 0 <= i < ZZX_deg(self.__fld_numerator.x):
                # calculate the new numerator. if our old entry is
                # p/q, and the lcm is k, it's just pk/q, which we
                # also know is integral -- so we can use mpz_divexact
                # below
                tmp_rational = <Rational>(coeff_list[i])
                mpz_mul(tmp_integer.value, mpq_numref(tmp_rational.value),
                        denom_temp.value)
                mpz_divexact(tmp_integer.value, tmp_integer.value,
                             mpq_denref(tmp_rational.value))

                # now set the coefficient of self
                tmp_integer._to_ZZ(&ntl_temp)
                ZZX_SetCoeff(self.__numerator, i, ntl_temp)

        return 0  # No error


    def __abs__(self):
        r"""
        Return the numerical absolute value of this number field element
        with respect to the first archimedean embedding, to double
        precision.

        This is the ``abs( )`` Python function. If you want a
        different embedding or precision, use
        ``self.abs(...)``.

        EXAMPLES::

            sage: k.<a> = NumberField(x^3 - 2)
            sage: abs(a)
            1.25992104989487
            sage: abs(a)^3
            2.00000000000000
            sage: a.abs(prec=128)
            1.2599210498948731647672106072782283506
        """
        return self.abs(prec=53, i=0)

    def abs(self, prec=53, i=0):
        r"""
        Return the absolute value of this element with respect to the
        `i`-th complex embedding of parent, to the given precision.

        If prec is 53 (the default), then the complex double field is
        used; otherwise the arbitrary precision (but slow) complex
        field is used.

        INPUT:


        -  ``prec`` - (default: 53) integer bits of precision

        -  ``i`` - (default: ) integer, which embedding to
           use


        EXAMPLES::

            sage: z = CyclotomicField(7).gen()
            sage: abs(z)
            1.00000000000000
            sage: abs(z^2 + 17*z - 3)
            16.0604426799931
            sage: K.<a> = NumberField(x^3+17)
            sage: abs(a)
            2.57128159065824
            sage: a.abs(prec=100)
            2.5712815906582353554531872087
            sage: a.abs(prec=100,i=1)
            2.5712815906582353554531872087
            sage: a.abs(100, 2)
            2.5712815906582353554531872087

        Here's one where the absolute value depends on the embedding.

        ::

            sage: K.<b> = NumberField(x^2-2)
            sage: a = 1 + b
            sage: a.abs(i=0)
            0.414213562373095
            sage: a.abs(i=1)
            2.41421356237309
        """
        P = self.number_field().complex_embeddings(prec)[i]
        return abs(P(self))

    def abs_non_arch(self, P, prec=None):
        r"""
        Return the non-archimedean absolute value of this element with
        respect to the prime `P`, to the given precision.

        INPUT:

        -  ``P`` - a prime ideal of the parent of self

        - ``prec`` (int) -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        (real) the non-archimedean absolute value of this element with
        respect to the prime `P`, to the given precision.  This is the
        normalised absolute value, so that the underlying prime number
        `p` has absolute value `1/p`.


        EXAMPLES::

            sage: K.<a> = NumberField(x^2+5)
            sage: [1/K(2).abs_non_arch(P) for P in K.primes_above(2)]
            [2.00000000000000]
            sage: [1/K(3).abs_non_arch(P) for P in K.primes_above(3)]
            [3.00000000000000, 3.00000000000000]
            sage: [1/K(5).abs_non_arch(P) for P in K.primes_above(5)]
            [5.00000000000000]

        A relative example::

            sage: L.<b> = K.extension(x^2-5)
            sage: [b.abs_non_arch(P) for P in L.primes_above(b)]
            [0.447213595499958, 0.447213595499958]
        """
        from sage.rings.real_mpfr import RealField
        if prec is None:
            R = RealField()
        else:
            R = RealField(prec)

        if self.is_zero():
            return R.zero_element()
        val = self.valuation(P)
        nP = P.residue_class_degree()*P.absolute_ramification_index()
        return R(P.absolute_norm()) ** (-R(val) / R(nP))

    def coordinates_in_terms_of_powers(self):
        r"""
        Let `\alpha` be self. Return a callable object (of type
        :class:`~CoordinateFunction`) that takes any element of the
        parent of self in `\QQ(\alpha)` and writes it in terms of the
        powers of `\alpha`: `1, \alpha, \alpha^2, ...`.

        (NOT CACHED).

        EXAMPLES:

        This function allows us to write elements of a number
        field in terms of a different generator without having to construct
        a whole separate number field.

        ::

            sage: y = polygen(QQ,'y'); K.<beta> = NumberField(y^3 - 2); K
            Number Field in beta with defining polynomial y^3 - 2
            sage: alpha = beta^2 + beta + 1
            sage: c = alpha.coordinates_in_terms_of_powers(); c
            Coordinate function that writes elements in terms of the powers of beta^2 + beta + 1
            sage: c(beta)
            [-2, -3, 1]
            sage: c(alpha)
            [0, 1, 0]
            sage: c((1+beta)^5)
            [3, 3, 3]
            sage: c((1+beta)^10)
            [54, 162, 189]

        This function works even if self only generates a subfield of this
        number field.

        ::

            sage: k.<a> = NumberField(x^6 - 5)
            sage: alpha = a^3
            sage: c = alpha.coordinates_in_terms_of_powers()
            sage: c((2/3)*a^3 - 5/3)
            [-5/3, 2/3]
            sage: c
            Coordinate function that writes elements in terms of the powers of a^3
            sage: c(a)
            Traceback (most recent call last):
            ...
            ArithmeticError: vector is not in free module
        """
        K = self.number_field()
        V, from_V, to_V = K.absolute_vector_space()
        h = K(1)
        B = [to_V(h)]
        f = self.absolute_minpoly()
        for i in range(f.degree()-1):
            h *= self
            B.append(to_V(h))
        W = V.span_of_basis(B)
        return CoordinateFunction(self, W, to_V)

    def complex_embeddings(self, prec=53):
        """
        Return the images of this element in the floating point complex
        numbers, to the given bits of precision.

        INPUT:


        -  ``prec`` - integer (default: 53) bits of precision


        EXAMPLES::

            sage: k.<a> = NumberField(x^3 - 2)
            sage: a.complex_embeddings()
            [-0.629960524947437 - 1.09112363597172*I, -0.629960524947437 + 1.09112363597172*I, 1.25992104989487]
            sage: a.complex_embeddings(10)
            [-0.63 - 1.1*I, -0.63 + 1.1*I, 1.3]
            sage: a.complex_embeddings(100)
            [-0.62996052494743658238360530364 - 1.0911236359717214035600726142*I, -0.62996052494743658238360530364 + 1.0911236359717214035600726142*I, 1.2599210498948731647672106073]
        """
        phi = self.number_field().complex_embeddings(prec)
        return [f(self) for f in phi]

    def complex_embedding(self, prec=53, i=0):
        """
        Return the i-th embedding of self in the complex numbers, to the
        given precision.

        EXAMPLES::

            sage: k.<a> = NumberField(x^3 - 2)
            sage: a.complex_embedding()
            -0.629960524947437 - 1.09112363597172*I
            sage: a.complex_embedding(10)
            -0.63 - 1.1*I
            sage: a.complex_embedding(100)
            -0.62996052494743658238360530364 - 1.0911236359717214035600726142*I
            sage: a.complex_embedding(20, 1)
            -0.62996 + 1.0911*I
            sage: a.complex_embedding(20, 2)
            1.2599
        """
        return self.number_field().complex_embeddings(prec)[i](self)

    def is_unit(self):
        """
        Return ``True`` if ``self`` is a unit in the ring where it is defined.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - x - 1)
            sage: OK = K.ring_of_integers()
            sage: OK(a).is_unit()
            True
            sage: OK(13).is_unit()
            False
            sage: K(13).is_unit()
            True

        It also works for relative fields and orders::

            sage: K.<a,b> = NumberField([x^2 - 3, x^4 + x^3 + x^2 + x + 1])
            sage: OK = K.ring_of_integers()
            sage: OK(b).is_unit()
            True
            sage: OK(a).is_unit()
            False
            sage: a.is_unit()
            True
        """
        if self.parent().is_field():
            return bool(self)
        return self.norm().is_unit()

    def is_norm(self, L, element=False, proof=True):
        r"""
        Determine whether self is the relative norm of an element
        of L/K, where K is self.parent().

        INPUT:

         - L -- a number field containing K=self.parent()
         - element -- True or False, whether to also output an element
           of which self is a norm
         - proof -- If True, then the output is correct unconditionally.
           If False, then the output is correct under GRH.

        OUTPUT:

        If element is False, then the output is a boolean B, which is
        True if and only if self is the relative norm of an element of L
        to K.
        If element is False, then the output is a pair (B, x), where
        B is as above. If B is True, then x is an element of L such that
        self == x.norm(K). Otherwise, x is None.

        ALGORITHM:

        Uses PARI's rnfisnorm. See self._rnfisnorm().

        EXAMPLES::

            sage: K.<beta> = NumberField(x^3+5)
            sage: Q.<X> = K[]
            sage: L = K.extension(X^2+X+beta, 'gamma')
            sage: (beta/2).is_norm(L)
            False
            sage: beta.is_norm(L)
            True

        With a relative base field::

            sage: K.<a, b> = NumberField([x^2 - 2, x^2 - 3])
            sage: L.<c> = K.extension(x^2 - 5)
            sage: (2*a*b).is_norm(L)
            True
            sage: _, v = (2*b*a).is_norm(L, element=True)
            sage: v.norm(K) == 2*a*b
            True

        Non-Galois number fields::

            sage: K.<a> = NumberField(x^2 + x + 1)
            sage: Q.<X> = K[]
            sage: L.<b> = NumberField(X^4 + a + 2)
            sage: (a/4).is_norm(L)
            True
            sage: (a/2).is_norm(L)
            Traceback (most recent call last):
            ...
            NotImplementedError: is_norm is not implemented unconditionally for norms from non-Galois number fields
            sage: (a/2).is_norm(L, proof=False)
            False

            sage: K.<a> = NumberField(x^3 + x + 1)
            sage: Q.<X> = K[]
            sage: L.<b> = NumberField(X^4 + a)
            sage: t = (-a).is_norm(L, element=True); t
            (True, -b^3 + 1)
            sage: t[1].norm(K)
            -a

        AUTHORS:

        - Craig Citro (2008-04-05)

        - Marco Streng (2010-12-03)
        """
        if not element:
            return self.is_norm(L, element=True, proof=proof)[0]

        K = self.parent()
        from sage.rings.number_field.all import is_AbsoluteNumberField, \
                                                is_NumberField
        if not is_NumberField(L):
            raise ValueError, "L (=%s) must be a NumberField in is_norm" % L

        if is_AbsoluteNumberField(L):
            Lrel = L.relativize(K.hom(L), ('a', 'b'))
            b, x = self.is_norm(Lrel, element=True, proof=proof)
            h = Lrel.structure()[0]
            return b, h(x)

        if L.relative_degree() == 1 or self.is_zero():
            return True, L(self)

        a, b = self._rnfisnorm(L, proof=proof)
        if b == 1:
            assert a.norm(K) == self
            return True, a

        if L.is_galois_relative():
            return False, None

        # The following gives the galois closure of K/QQ, but the galois
        # closure of K/self.parent() would suffice.
        M = L.galois_closure('a')
        from sage.functions.log import log
        from sage.functions.other import floor
        extra_primes = floor(12*log(abs(M.discriminant()))**2)
        a, b = self._rnfisnorm(L, proof=proof, extra_primes=extra_primes)
        if b == 1:
            assert a.norm(K) == self
            return True, a

        if proof:
            raise NotImplementedError, "is_norm is not implemented unconditionally for norms from non-Galois number fields"
        return False, None

    def _rnfisnorm(self, L, proof=True, extra_primes=0):
        r"""
        Gives the output of the PARI function rnfisnorm.

        This tries to decide whether the number field element self is
        the norm of some x in the extension L/K (with K = self.parent()).

        The output is a pair (x, q), where self = Norm(x)*q. The
        algorithm looks for a solution x that is an S-integer, with S
        a list of places of L containing at least the ramified primes,
        the generators of the class group of L, as well as those primes
        dividing self.

        If L/K is Galois, then this is enough; otherwise,
        extra_primes is used to add more primes to S: all the places
        above the primes p <= extra_primes (resp. p|extra_primes) if
        extra_primes > 0 (resp. extra_primes < 0).

        The answer is guaranteed (i.e., self is a norm iff q = 1) if the
        field is Galois, or, under GRH, if S contains all primes less
        than 12log^2|\disc(M)|, where M is the normal closure of L/K.

        INPUT:

         - L -- a relative number field with base field self.parent()
         - proof -- whether to certify outputs of PARI init functions.
           If false, truth of the output depends on GRH.
         - extra_primes -- an integer as explained above.

        OUTPUT:

        A pair (x, q) with x in L and q in K as explained above
        such that self == x.norm(K)*q.

        ALGORITHM:

        Uses PARI's rnfisnorm.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + x^2 - 2*x - 1, 'a')
            sage: P.<X> = K[]
            sage: L = NumberField(X^2 + a^2 + 2*a + 1, 'b')
            sage: K(17)._rnfisnorm(L)
            ((a^2 - 2)*b - 4, 1)

            sage: K.<a> = NumberField(x^3 + x + 1)
            sage: Q.<X> = K[]
            sage: L.<b> = NumberField(X^4 + a)
            sage: t = (-a)._rnfisnorm(L); t
            (-b^3 + 1, 1)
            sage: t[0].norm(K)
            -a
            sage: t = K(3)._rnfisnorm(L); t
            ((a^2 + 1)*b^3 + b^2 - a*b + a^2 + 1, -3*a)
            sage: t[0].norm(K)*t[1]
            3

        An example where the base field is a relative field::

            sage: K.<a, b> = NumberField([x^2 - 2, x^2 - 3])
            sage: L.<c> = K.extension(x^3 + 2)
            sage: s = 2*a + b
            sage: t = s._rnfisnorm(L)
            sage: t[1] == 1 # True iff s is a norm
            False
            sage: s == t[0].norm(K)*t[1]
            True

        AUTHORS:

        - Craig Citro (2008-04-05)

        - Marco Streng (2010-12-03)

        - Francis Clarke (2010-12-26)
        """
        K = self.parent()
        from sage.rings.number_field.all import is_RelativeNumberField
        if (not is_RelativeNumberField(L)) or L.base_field() != K:
            raise ValueError, "L (=%s) must be a relative number field with base field K (=%s) in rnfisnorm" % (L, K)

        rnf_data = K.pari_rnfnorm_data(L, proof=proof)
        x, q = self._pari_().rnfisnorm(rnf_data)
        return L(x), K(q)

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 1)
            sage: RR(a^2)
            -1.00000000000000
            sage: RR(a)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce a to a rational
            sage: (a^2)._mpfr_(RR)
            -1.00000000000000

        Verify that :trac:`#13005` has been fixed::

            sage: K.<a> = NumberField(x^2-5)
            sage: RR(K(1))
            1.00000000000000
            sage: RR(a)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce a to a rational
            sage: K.<a> = NumberField(x^3+2, embedding=-1.25)
            sage: RR(a)
            -1.25992104989487
            sage: RealField(prec=100)(a)
            -1.2599210498948731647672106073
        """
        if self.parent().coerce_embedding() is None:
            return R(self.base_ring()(self))
        else:
            return R(R.complex_field()(self))

    def __float__(self):
        """
        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 1)
            sage: float(a^2)
            -1.0
            sage: float(a)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce a to a rational
            sage: (a^2).__float__()
            -1.0
            sage: k.<a> = NumberField(x^2 + 1,embedding=I)
            sage: float(a)
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce to a real number
        """
        if self.parent().coerce_embedding() is None:
            return float(self.base_ring()(self))
        else:
            c = complex(self)
            if c.imag == 0:
                return c.real
            raise TypeError('unable to coerce to a real number')

    def _complex_double_(self, CDF):
        """
        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 1)
            sage: abs(CDF(a))
            1.0
        """
        return CDF(CC(self))

    def __complex__(self):
        """
        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 1)
            sage: complex(a)
            1j
            sage: a.__complex__()
            1j
        """
        return complex(CC(self))

    def factor(self):
        """
        Return factorization of this element into prime elements and a unit.

        OUTPUT:

        (Factorization) If all the prime ideals in the support are
        principal, the output is a Factorization as a product of prime
        elements raised to appropriate powers, with an appropriate
        unit factor.

        Raise ValueError if the factorization of the
        ideal (self) contains a non-principal prime ideal.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2+1)
            sage: (6*i + 6).factor()
            (-i) * (i + 1)^3 * 3

        In the following example, the class number is 2.  If a factorization
        in prime elements exists, we will find it::

            sage: K.<a> = NumberField(x^2-10)
            sage: factor(169*a + 531)
            (-6*a - 19) * (-3*a - 1) * (-2*a + 9)
            sage: factor(K(3))
            Traceback (most recent call last):
            ...
            ValueError: Non-principal ideal in factorization

        Factorization of 0 is not allowed::

            sage: K.<i> = QuadraticField(-1)
            sage: K(0).factor()
            Traceback (most recent call last):
            ...
            ArithmeticError: Factorization of 0 not defined.

        """
        if self.is_zero():
            raise ArithmeticError, "Factorization of 0 not defined."

        K = self.parent()
        fac = K.ideal(self).factor()
        # Check whether all prime ideals in `fac` are principal
        for P,e in fac:
            if not P.is_principal():
                raise ValueError, "Non-principal ideal in factorization"
        element_fac = [(P.gens_reduced()[0],e) for P,e in fac]
        # Compute the product of the p^e to figure out the unit
        from sage.misc.all import prod
        element_product = prod([p**e for p,e in element_fac], K(1))
        from sage.structure.all import Factorization
        return Factorization(element_fac, unit=self/element_product)

    def is_totally_positive(self):
        """
        Returns True if self is positive for all real embeddings of its
        parent number field. We do nothing at complex places, so e.g. any
        element of a totally complex number field will return True.

        EXAMPLES::

            sage: F.<b> = NumberField(x^3-3*x-1)
            sage: b.is_totally_positive()
            False
            sage: (b^2).is_totally_positive()
            True

        TESTS:

        Check that the output is correct even for numbers that are
        very close to zero (ticket #9596)::

            sage: K.<sqrt2> = QuadraticField(2)
            sage: a = 30122754096401; b = 21300003689580
            sage: (a/b)^2 > 2
            True
            sage: (a/b+sqrt2).is_totally_positive()
            True
            sage: r = RealField(3020)(2).sqrt()*2^3000
            sage: a = floor(r)/2^3000
            sage: b = ceil(r)/2^3000
            sage: (a+sqrt2).is_totally_positive()
            False
            sage: (b+sqrt2).is_totally_positive()
            True

        Check that 0 is handled correctly::

            sage: K.<a> = NumberField(x^5+4*x+1)
            sage: K(0).is_totally_positive()
            False
        """
        for v in self.number_field().embeddings(sage.rings.qqbar.AA):
            if v(self) <= 0:
                return False
        return True

    def is_square(self, root=False):
        """
        Return True if self is a square in its parent number field and
        otherwise return False.

        INPUT:


        -  ``root`` - if True, also return a square root (or
           None if self is not a perfect square)


        EXAMPLES::

            sage: m.<b> = NumberField(x^4 - 1789)
            sage: b.is_square()
            False
            sage: c = (2/3*b + 5)^2; c
            4/9*b^2 + 20/3*b + 25
            sage: c.is_square()
            True
            sage: c.is_square(True)
            (True, 2/3*b + 5)

        We also test the functional notation.

        ::

            sage: is_square(c, True)
            (True, 2/3*b + 5)
            sage: is_square(c)
            True
            sage: is_square(c+1)
            False
        """
        v = self.sqrt(all=True)
        t = len(v) > 0
        if root:
            if t:
                return t, v[0]
            else:
                return False, None
        else:
            return t

    def sqrt(self, all=False):
        """
        Returns the square root of this number in the given number field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 - 3)
            sage: K(3).sqrt()
            a
            sage: K(3).sqrt(all=True)
            [a, -a]
            sage: K(a^10).sqrt()
            9*a
            sage: K(49).sqrt()
            7
            sage: K(1+a).sqrt()
            Traceback (most recent call last):
            ...
            ValueError: a + 1 not a square in Number Field in a with defining polynomial x^2 - 3
            sage: K(0).sqrt()
            0
            sage: K((7+a)^2).sqrt(all=True)
            [a + 7, -a - 7]

        ::

            sage: K.<a> = CyclotomicField(7)
            sage: a.sqrt()
            a^4

        ::

            sage: K.<a> = NumberField(x^5 - x + 1)
            sage: (a^4 + a^2 - 3*a + 2).sqrt()
            a^3 - a^2

        ALGORITHM: Use PARI to factor `x^2` - ``self`` in `K`.
        """
        # For now, use pari's factoring abilities
        R = self.number_field()['t']
        f = R([-self, 0, 1])
        roots = f.roots()
        if all:
            return [r[0] for r in roots]
        elif len(roots) > 0:
            return roots[0][0]
        else:
            try:
                # This is what integers, rationals do...
                from sage.all import SR, sqrt
                return sqrt(SR(self))
            except TypeError:
                raise ValueError, "%s not a square in %s"%(self, self._parent)

    def nth_root(self, n, all=False):
        r"""
        Return an nth root of self in the given number field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4-7)
            sage: K(7).nth_root(2)
            a^2
            sage: K((a-3)^5).nth_root(5)
            a - 3

        ALGORITHM: Use PARI to factor `x^n` - ``self`` in `K`.
        """
        R = self.number_field()['t']
        if not self:
            return [self] if all else self
        f = (R.gen(0) << (n-1)) - self
        roots = f.roots()
        if all:
            return [r[0] for r in roots]
        elif len(roots) > 0:
            return roots[0][0]
        else:
            raise ValueError, "%s not a %s-th root in %s"%(self, n, self._parent)

    def __pow__(base, exp, dummy):
        """
        EXAMPLES::

            sage: K.<sqrt2> = QuadraticField(2)
            sage: sqrt2^2
            2
            sage: sqrt2^5
            4*sqrt2
            sage: (1+sqrt2)^100
            66992092050551637663438906713182313772*sqrt2 + 94741125149636933417873079920900017937
            sage: (1+sqrt2)^-1
            sqrt2 - 1

        If the exponent is not integral, perform this operation in
        the symbolic ring::

            sage: sqrt2^(1/5)
            2^(1/10)
            sage: sqrt2^sqrt2
            2^(1/2*sqrt(2))

        Sage follows Python's convention 0^0 = 1::

            sage: a = K(0)^0; a
            1
            sage: a.parent()
            Number Field in sqrt2 with defining polynomial x^2 - 2

        TESTS::

            sage: 2^I
            2^I
        """
        if (PY_TYPE_CHECK(base, NumberFieldElement) and
            (PY_TYPE_CHECK(exp, Integer) or PY_TYPE_CHECK_EXACT(exp, int) or exp in ZZ)):
            return generic_power_c(base, exp, None)
        else:
            cbase, cexp = canonical_coercion(base, exp)
            if not isinstance(cbase, NumberFieldElement):
                return cbase ** cexp
            # Return a symbolic expression.
            # We use the hold=True keyword argument to prevent the
            # symbolics library from trying to simplify this expression
            # again. This would lead to infinite loops otherwise.
            from sage.symbolic.ring import SR
            try:
                res = QQ(base)**exp
            except TypeError:
                pass
            else:
                if res.parent() is not SR:
                    return parent(cbase)(res)
                return res
            sbase = SR(base)
            if sbase.operator() is operator.pow:
                nbase, pexp = sbase.operands()
                return nbase.power(pexp * exp, hold=True)
            else:
                return sbase.power(exp, hold=True)

    cdef void _reduce_c_(self):
        """
        Pull out common factors from the numerator and denominator!
        """
        cdef ZZ_c gcd
        cdef ZZ_c t1
        cdef ZZX_c t2
        ZZX_content(t1, self.__numerator)
        ZZ_GCD(gcd, t1, self.__denominator)
        if ZZ_sign(gcd) != ZZ_sign(self.__denominator):
            ZZ_negate(t1, gcd)
            gcd = t1
        ZZX_div_ZZ(t2, self.__numerator, gcd)
        ZZ_div(t1, self.__denominator, gcd)
        self.__numerator = t2
        self.__denominator = t1

    cpdef ModuleElement _add_(self, ModuleElement right):
        r"""
        EXAMPLE::

            sage: K.<s> = QuadraticField(2)
            sage: s + s # indirect doctest
            2*s
            sage: s + ZZ(3) # indirect doctest
            s + 3
        """
        cdef NumberFieldElement x
        cdef NumberFieldElement _right = right
        x = self._new()
        ZZ_mul(x.__denominator, self.__denominator, _right.__denominator)
        cdef ZZX_c t1, t2
        ZZX_mul_ZZ(t1, self.__numerator, _right.__denominator)
        ZZX_mul_ZZ(t2, _right.__numerator, self.__denominator)
        ZZX_add(x.__numerator, t1, t2)
        x._reduce_c_()
        return x

    cpdef ModuleElement _sub_(self, ModuleElement right):
        r"""
        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + 2)
            sage: (a/2) - (a + 3) # indirect doctest
            -1/2*a - 3
        """
        cdef NumberFieldElement x
        cdef NumberFieldElement _right = right
        x = self._new()
        ZZ_mul(x.__denominator, self.__denominator, _right.__denominator)
        cdef ZZX_c t1, t2
        ZZX_mul_ZZ(t1, self.__numerator, _right.__denominator)
        ZZX_mul_ZZ(t2, _right.__numerator, self.__denominator)
        ZZX_sub(x.__numerator, t1, t2)
        x._reduce_c_()
        return x

    cpdef RingElement _mul_(self, RingElement right):
        """
        Returns the product of self and other as elements of a number
        field.

        EXAMPLES::

            sage: C.<zeta12>=CyclotomicField(12)
            sage: zeta12*zeta12^11
            1
            sage: G.<a> = NumberField(x^3 + 2/3*x + 1)
            sage: a^3 # indirect doctest
            -2/3*a - 1
            sage: a^3+a # indirect doctest
            1/3*a - 1
        """
        cdef NumberFieldElement x
        cdef NumberFieldElement _right = right
        cdef ZZX_c temp
        cdef ZZ_c temp1
        x = self._new()
        sig_on()
        # MulMod doesn't handle non-monic polynomials.
        # Therefore, we handle the non-monic case entirely separately.

        if ZZ_IsOne(ZZX_LeadCoeff(self.__fld_numerator.x)):
            ZZ_mul(x.__denominator, self.__denominator, _right.__denominator)
            ZZX_MulMod(x.__numerator, self.__numerator, _right.__numerator, self.__fld_numerator.x)
        else:
            ZZ_mul(x.__denominator, self.__denominator, _right.__denominator)
            ZZX_mul(x.__numerator, self.__numerator, _right.__numerator)
            if ZZX_deg(x.__numerator) >= ZZX_deg(self.__fld_numerator.x):
                ZZX_mul_ZZ( x.__numerator, x.__numerator, self.__fld_denominator.x )
                ZZX_mul_ZZ( temp, self.__fld_numerator.x, x.__denominator )
                ZZ_power(temp1,ZZX_LeadCoeff(temp),ZZX_deg(x.__numerator)-ZZX_deg(self.__fld_numerator.x)+1)
                ZZX_PseudoRem(x.__numerator, x.__numerator, temp)
                ZZ_mul(x.__denominator, x.__denominator, self.__fld_denominator.x)
                ZZ_mul(x.__denominator, x.__denominator, temp1)
        sig_off()
        x._reduce_c_()
        return x

        #NOTES: In LiDIA, they build a multiplication table for the
        #number field, so it's not necessary to reduce modulo the
        #defining polynomial every time:
        #     src/number_fields/algebraic_num/order.cc: compute_table
        # but asymptotically fast poly multiplication means it's
        # actually faster to *not* build a table!?!

    cpdef RingElement _div_(self, RingElement right):
        """
        Returns the quotient of self and other as elements of a number
        field.

        EXAMPLES::

            sage: C.<I>=CyclotomicField(4)
            sage: 1/I # indirect doctest
            -I
            sage: I/0 # indirect doctest
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        ::

            sage: G.<a> = NumberField(x^3 + 2/3*x + 1)
            sage: a/a # indirect doctest
            1
            sage: 1/a # indirect doctest
            -a^2 - 2/3
            sage: a/0 # indirect doctest
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Number field element division by zero
        """
        cdef NumberFieldElement x
        cdef NumberFieldElement _right = right
        cdef ZZX_c inv_num
        cdef ZZ_c inv_den
        cdef ZZX_c temp
        cdef ZZ_c temp1
        if not _right:
            raise ZeroDivisionError, "Number field element division by zero"
        x = self._new()
        sig_on()
        _right._invert_c_(&inv_num, &inv_den)
        if ZZ_IsOne(ZZX_LeadCoeff(self.__fld_numerator.x)):
            ZZ_mul(x.__denominator, self.__denominator, inv_den)
            ZZX_MulMod(x.__numerator, self.__numerator, inv_num, self.__fld_numerator.x)
        else:
            ZZ_mul(x.__denominator, self.__denominator, inv_den)
            ZZX_mul(x.__numerator, self.__numerator, inv_num)
            if ZZX_deg(x.__numerator) >= ZZX_deg(self.__fld_numerator.x):
                ZZX_mul_ZZ( x.__numerator, x.__numerator, self.__fld_denominator.x )
                ZZX_mul_ZZ( temp, self.__fld_numerator.x, x.__denominator )
                ZZ_power(temp1,ZZX_LeadCoeff(temp),ZZX_deg(x.__numerator)-ZZX_deg(self.__fld_numerator.x)+1)
                ZZX_PseudoRem(x.__numerator, x.__numerator, temp)
                ZZ_mul(x.__denominator, x.__denominator, self.__fld_denominator.x)
                ZZ_mul(x.__denominator, x.__denominator, temp1)
        x._reduce_c_()
        sig_off()
        return x

    def __floordiv__(self, other):
        """
        Return the quotient of self and other. Since these are field
        elements the floor division is exactly the same as usual division.

        EXAMPLES::

            sage: m.<b> = NumberField(x^4 + x^2 + 2/3)
            sage: c = (1+b) // (1-b); c
            3/4*b^3 + 3/4*b^2 + 3/2*b + 1/2
            sage: (1+b) / (1-b) == c
            True
            sage: c * (1-b)
            b + 1
        """
        return self / other

    def __nonzero__(self):
        """
        Return True if this number field element is nonzero.

        EXAMPLES::

            sage: m.<b> = CyclotomicField(17)
            sage: m(0).__nonzero__()
            False
            sage: b.__nonzero__()
            True

        Nonzero is used by the bool command::

            sage: bool(b + 1)
            True
            sage: bool(m(0))
            False
        """
        return not IsZero_ZZX(self.__numerator)

    cpdef ModuleElement _neg_(self):
        r"""
        EXAMPLE::

            sage: K.<a> = NumberField(x^3 + 2)
            sage: -a # indirect doctest
            -a
        """
        cdef NumberFieldElement x
        x = self._new()
        ZZX_mul_long(x.__numerator, self.__numerator, -1)
        x.__denominator = self.__denominator
        return x

    def __copy__(self):
        r"""
        EXAMPLE::

            sage: K.<a> = NumberField(x^3 + 2)
            sage: b = copy(a)
            sage: b == a
            True
            sage: b is a
            False
        """
        cdef NumberFieldElement x
        x = self._new()
        x.__numerator = self.__numerator
        x.__denominator = self.__denominator
        return x

    def __int__(self):
        """
        Attempt to convert this number field element to a Python integer,
        if possible.

        EXAMPLES::

            sage: C.<I>=CyclotomicField(4)
            sage: int(1/I)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to int
            sage: int(I*I)
            -1

        ::

            sage: K.<a> = NumberField(x^10 - x - 1)
            sage: int(a)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to int
            sage: int(K(9390283))
            9390283

        The semantics are like in Python, so the value does not have to
        preserved.

        ::

            sage: int(K(393/29))
            13
        """
        return int(self.polynomial())

    def __long__(self):
        """
        Attempt to convert this number field element to a Python long, if
        possible.

        EXAMPLES::

            sage: K.<a> = NumberField(x^10 - x - 1)
            sage: long(a)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to long
            sage: long(K(1234))
            1234L

        The value does not have to be preserved, in the case of fractions.

        ::

            sage: long(K(393/29))
            13L
        """
        return long(self.polynomial())

    cdef void _invert_c_(self, ZZX_c *num, ZZ_c *den):
        """
        Computes the numerator and denominator of the multiplicative
        inverse of this element.

        Suppose that this element is x/d and the parent mod'ding polynomial
        is M/D. The NTL function XGCD( r, s, t, a, b ) computes r,s,t such
        that `r=s*a+t*b`. We compute XGCD( r, s, t, x\*D, M\*d )
        and set num=s\*D\*d den=r

        EXAMPLES:

        I'd love to, but since we are dealing with c-types, I
        can't at this level. Check __invert__ for doc-tests that rely
        on this functionality.
        """
        cdef ZZX_c t # unneeded except to be there
        cdef ZZX_c a, b
        ZZX_mul_ZZ( a, self.__numerator, self.__fld_denominator.x )
        ZZX_mul_ZZ( b, self.__fld_numerator.x, self.__denominator )
        ZZX_XGCD( den[0], num[0],  t, a, b, 1 )
        ZZX_mul_ZZ( num[0], num[0], self.__fld_denominator.x )
        ZZX_mul_ZZ( num[0], num[0], self.__denominator )

    def __invert__(self):
        """
        Returns the multiplicative inverse of self in the number field.

        EXAMPLES::

            sage: C.<I>=CyclotomicField(4)
            sage: ~I
            -I
            sage: (2*I).__invert__()
            -1/2*I
        """
        if IsZero_ZZX(self.__numerator):
            raise ZeroDivisionError
        cdef NumberFieldElement x
        x = self._new()
        self._invert_c_(&x.__numerator, &x.__denominator)
        x._reduce_c_()
        return x

    def _integer_(self, Z=None):
        """
        Returns an integer if this element is actually an integer.

        EXAMPLES::

            sage: C.<I>=CyclotomicField(4)
            sage: (~I)._integer_()
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce -I to an integer
            sage: (2*I*I)._integer_()
            -2
        """
        if ZZX_deg(self.__numerator) >= 1:
            raise TypeError, "Unable to coerce %s to an integer"%self
        return ZZ(self._rational_())

    def _rational_(self):
        """
        Returns a rational number if this element is actually a rational
        number.

        EXAMPLES::

            sage: C.<I>=CyclotomicField(4)
            sage: (~I)._rational_()
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce -I to a rational
            sage: (I*I/2)._rational_()
            -1/2
        """
        if ZZX_deg(self.__numerator) >= 1:
            raise TypeError, "Unable to coerce %s to a rational"%self
        cdef Integer num
        num = PY_NEW(Integer)
        ZZX_getitem_as_mpz(&num.value, &self.__numerator, 0)
        return num / (<IntegerRing_class>ZZ)._coerce_ZZ(&self.__denominator)

    def _symbolic_(self, SR):
        """
        If an embedding into CC is specified, then a representation of this
        element can be made in the symbolic ring (assuming roots of the
        minimal polynomial can be found symbolically).

        EXAMPLES::

            sage: K.<a> = QuadraticField(2)
            sage: SR(a) # indirect doctest
            sqrt(2)
            sage: SR(3*a-5) # indirect doctest
            3*sqrt(2) - 5
            sage: K.<a> = QuadraticField(2, embedding=-1.4)
            sage: SR(a) # indirect doctest
            -sqrt(2)
            sage: K.<a> = NumberField(x^2 - 2)
            sage: SR(a) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: An embedding into RR or CC must be specified.

        Now a more complicated example::

            sage: K.<a> = NumberField(x^3 + x - 1, embedding=0.68)
            sage: b = SR(a); b # indirect doctest
            1/3*(3*(1/18*sqrt(31)*sqrt(3) + 1/2)^(2/3) - 1)/(1/18*sqrt(31)*sqrt(3) + 1/2)^(1/3)

            sage: (b^3 + b - 1).simplify_radical()
            0

        Make sure we got the right one::

            sage: CC(a)
            0.682327803828019
            sage: CC(b)
            0.682327803828019

        Special case for cyclotomic fields::

            sage: K.<zeta> = CyclotomicField(19)
            sage: SR(zeta) # indirect doctest
            e^(2/19*I*pi)
            sage: CC(zeta)
            0.945817241700635 + 0.324699469204683*I
            sage: CC(SR(zeta))
            0.945817241700635 + 0.324699469204683*I

            sage: SR(zeta^5 + 2)
            e^(10/19*I*pi) + 2

        For degree greater than 5, sometimes Galois theory prevents a
        closed-form solution.  In this case, a numerical approximation
        is used::

            sage: K.<a> = NumberField(x^5-x+1, embedding=-1)
            sage: SR(a)
            -1.1673040153

        ::

            sage: K.<a> = NumberField(x^6-x^3-1, embedding=1)
            sage: SR(a)
            (1/2*sqrt(5) + 1/2)^(1/3)
        """
        if self.__symbolic is None:

            K = self._parent.fraction_field()

            gen = K.gen()
            if not self is gen:
                try:
                    # share the hard work...
                    gen_image = gen._symbolic_(SR)
                    self.__symbolic = self.polynomial()(gen_image)
                    return self.__symbolic
                except TypeError:
                    pass # we may still be able to do this particular element...

            embedding = K.specified_complex_embedding()
            if embedding is None:
                raise TypeError, "An embedding into RR or CC must be specified."

            if isinstance(K, number_field.NumberField_cyclotomic):
                # solution by radicals may be difficult, but we have a closed form
                from sage.all import exp, I, pi, ComplexField, RR
                CC = ComplexField(53)
                two_pi_i = 2 * pi * I
                k = ( K._n()*CC(K.gen()).log() / CC(two_pi_i) ).real().round() # n ln z / (2 pi i)
                gen_image = exp(k*two_pi_i/K._n())
                if self is gen:
                    self.__symbolic = gen_image
                else:
                    self.__symbolic = self.polynomial()(gen_image)
            else:
                # try to solve the minpoly and choose the closest root
                poly = self.minpoly()
                roots = []
                var = SR(poly.variable_name())
                for soln in SR(poly).solve(var, to_poly_solve=True):
                    if soln.lhs() == var:
                        roots.append(soln.rhs())
                if len(roots) != poly.degree():
                    raise TypeError, "Unable to solve by radicals."
                from number_field_morphisms import matching_root
                from sage.rings.complex_field import ComplexField
                gen_image = matching_root(roots, self, ambient_field=ComplexField(53), margin=2)
                if gen_image is not None:
                    self.__symbolic = gen_image
                else:
                    # should be rare, e.g. if there is insufficient precision
                    raise TypeError, "Unable to determine which root in SR is this element."

        return self.__symbolic

    def galois_conjugates(self, K):
        r"""
        Return all Gal(Qbar/Q)-conjugates of this number field element in
        the field K.

        EXAMPLES:

        In the first example the conjugates are obvious::

            sage: K.<a> = NumberField(x^2 - 2)
            sage: a.galois_conjugates(K)
            [a, -a]
            sage: K(3).galois_conjugates(K)
            [3]

        In this example the field is not Galois, so we have to pass to an
        extension to obtain the Galois conjugates.

        ::

            sage: K.<a> = NumberField(x^3 - 2)
            sage: c = a.galois_conjugates(K); c
            [a]
            sage: K.<a> = NumberField(x^3 - 2)
            sage: c = a.galois_conjugates(K.galois_closure('a1')); c
            [1/84*a1^4 + 13/42*a1, -1/252*a1^4 - 55/126*a1, -1/126*a1^4 + 8/63*a1]
            sage: c[0]^3
            2
            sage: parent(c[0])
            Number Field in a1 with defining polynomial x^6 + 40*x^3 + 1372
            sage: parent(c[0]).is_galois()
            True

        There is only one Galois conjugate of `\sqrt[3]{2}` in
        `\QQ(\sqrt[3]{2})`.

        ::

            sage: a.galois_conjugates(K)
            [a]

        Galois conjugates of `\sqrt[3]{2}` in the field
        `\QQ(\zeta_3,\sqrt[3]{2})`::

            sage: L.<a> = CyclotomicField(3).extension(x^3 - 2)
            sage: a.galois_conjugates(L)
            [a, (-zeta3 - 1)*a, zeta3*a]
        """
        f = self.absolute_minpoly()
        g = K['x'](f)
        return [a for a,_ in g.roots()]

    def conjugate(self):
        """
        Return the complex conjugate of the number field element.

        This is only well-defined for fields contained in CM fields
        (i.e. for totally real fields and CM fields). Recall that a CM
        field is a totally imaginary quadratic extension of a totally
        real field. For other fields, a ValueError is raised.

        EXAMPLES::

            sage: k.<I> = QuadraticField(-1)
            sage: I.conjugate()
            -I
            sage: (I/(1+I)).conjugate()
            -1/2*I + 1/2
            sage: z6 = CyclotomicField(6).gen(0)
            sage: (2*z6).conjugate()
            -2*zeta6 + 2

        The following example now works.

        ::

            sage: F.<b> = NumberField(x^2 - 2)
            sage: K.<j> = F.extension(x^2 + 1)
            sage: j.conjugate()
            -j

        Raise a ValueError if the field is not contained in a CM field.

        ::

            sage: K.<b> = NumberField(x^3 - 2)
            sage: b.conjugate()
            Traceback (most recent call last):
            ...
            ValueError: Complex conjugation is only well-defined for fields contained in CM fields.

        An example of a non-quadratic totally real field.

        ::

            sage: F.<a> = NumberField(x^4 + x^3 - 3*x^2 - x + 1)
            sage: a.conjugate()
            a

        An example of a non-cyclotomic CM field.

        ::

            sage: K.<a> = NumberField(x^4 - x^3 + 2*x^2 + x + 1)
            sage: a.conjugate()
            -1/2*a^3 - a - 1/2
            sage: (2*a^2 - 1).conjugate()
            a^3 - 2*a^2 - 2

        """

        nf = self.number_field()
        return nf.complex_conjugation()(self)

    def polynomial(self, var='x'):
        """
        Return the underlying polynomial corresponding to this number field
        element.

        The resulting polynomial is currently *not* cached.

        EXAMPLES::

            sage: K.<a> = NumberField(x^5 - x - 1)
            sage: f = (-2/3 + 1/3*a)^4; f
            1/81*a^4 - 8/81*a^3 + 8/27*a^2 - 32/81*a + 16/81
            sage: g = f.polynomial(); g
            1/81*x^4 - 8/81*x^3 + 8/27*x^2 - 32/81*x + 16/81
            sage: parent(g)
            Univariate Polynomial Ring in x over Rational Field

        Note that the result of this function is not cached (should this be
        changed?)::

            sage: g is f.polynomial()
            False
        """
        return QQ[var](self._coefficients())

    def __hash__(self):
        """
        Return hash of this number field element, which is just the
        hash of the underlying polynomial.

        EXAMPLE::

            sage: K.<b> = NumberField(x^3 - 2)
            sage: hash(b^2 + 1) == hash((b^2 + 1).polynomial()) # indirect doctest
            True
        """
        return hash(self.polynomial())

    def _coefficients(self):
        """
        Return the coefficients of the underlying polynomial corresponding
        to this number field element.

        OUTPUT:

        - a list whose length corresponding to the degree of this
          element written in terms of a generator.

        EXAMPLES:

            sage: K.<b> = NumberField(x^3 - 2)
            sage: (b^2 + 1)._coefficients()
            [1, 0, 1]
        """
        coeffs = []
        cdef Integer den = (<IntegerRing_class>ZZ)._coerce_ZZ(&self.__denominator)
        cdef Integer numCoeff
        cdef int i
        for i from 0 <= i <= ZZX_deg(self.__numerator):
            numCoeff = PY_NEW(Integer)
            ZZX_getitem_as_mpz(&numCoeff.value, &self.__numerator, i)
            coeffs.append( numCoeff / den )
        return coeffs

    cdef void _ntl_coeff_as_mpz(self, mpz_t* z, long i):
        if i > ZZX_deg(self.__numerator):
            mpz_set_ui(z[0], 0)
        else:
            ZZX_getitem_as_mpz(z, &self.__numerator, i)

    cdef void _ntl_denom_as_mpz(self, mpz_t* z):
        cdef Integer denom = (<IntegerRing_class>ZZ)._coerce_ZZ(&self.__denominator)
        mpz_set(z[0], denom.value)

    def denominator(self):
        """
        Return the denominator of this element, which is by definition the
        denominator of the corresponding polynomial representation. I.e.,
        elements of number fields are represented as a polynomial (in
        reduced form) modulo the modulus of the number field, and the
        denominator is the denominator of this polynomial.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(3)
            sage: a = 1/3 + (1/5)*z
            sage: print a.denominator()
            15
        """
        return (<IntegerRing_class>ZZ)._coerce_ZZ(&self.__denominator)

    def _set_multiplicative_order(self, n):
        """
        Set the multiplicative order of this number field element.

        .. warning::

           Use with caution - only for internal use! End users should
           never call this unless they have a very good reason to do
           so.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + x + 1)
            sage: a._set_multiplicative_order(3)
            sage: a.multiplicative_order()
            3

        You can be evil with this so be careful. That's why the function
        name begins with an underscore.

        ::

            sage: a._set_multiplicative_order(389)
            sage: a.multiplicative_order()
            389
        """
        self.__multiplicative_order = n

    def multiplicative_order(self):
        """
        Return the multiplicative order of this number field element.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(5)
            sage: z.multiplicative_order()
            5
            sage: (-z).multiplicative_order()
            10
            sage: (1+z).multiplicative_order()
            +Infinity

            sage: x = polygen(QQ)
            sage: K.<a>=NumberField(x^40 - x^20 + 4)
            sage: u = 1/4*a^30 + 1/4*a^10 + 1/2
            sage: u.multiplicative_order()
            6
            sage: a.multiplicative_order()
            +Infinity

        An example in a relative extension::

            sage: K.<a, b> = NumberField([x^2 + x + 1, x^2 - 3])
            sage: z = (a - 1)*b/3
            sage: z.multiplicative_order()
            12
            sage: z^12==1 and z^6!=1 and z^4!=1
            True

        """
        if self.__multiplicative_order is not None:
            return self.__multiplicative_order

        one = self.number_field().one_element()
        infinity = sage.rings.infinity.infinity

        if self == one:
            self.__multiplicative_order = ZZ(1)
            return self.__multiplicative_order
        if self == -one:
            self.__multiplicative_order = ZZ(2)
            return self.__multiplicative_order

        if isinstance(self.number_field(), number_field.NumberField_cyclotomic):
            t = self.number_field()._multiplicative_order_table()
            f = self.polynomial()
            if t.has_key(f):
                self.__multiplicative_order = t[f]
                return self.__multiplicative_order
            else:
                self.__multiplicative_order = sage.rings.infinity.infinity
                return self.__multiplicative_order

        if self.is_rational_c() or not self.is_integral() or not self.norm() ==1:
            self.__multiplicative_order = infinity
            return self.__multiplicative_order

        # Now we have a unit of norm 1, and check if it is a root of unity

        n = self.number_field().zeta_order()
        if not self**n ==1:
            self.__multiplicative_order = infinity
            return self.__multiplicative_order
        from sage.groups.generic import order_from_multiple
        self.__multiplicative_order = order_from_multiple(self,n,operation='*')
        return self.__multiplicative_order

    def additive_order(self):
        r"""
        Return the additive order of this element (i.e. infinity if
        self != 0, 1 if self == 0)

        EXAMPLES::

            sage: K.<u> = NumberField(x^4 - 3*x^2 + 3)
            sage: u.additive_order()
            +Infinity
            sage: K(0).additive_order()
            1
            sage: K.ring_of_integers().characteristic() # implicit doctest
            0
        """
        if self == 0: return 1
        else: return sage.rings.infinity.infinity

    cdef bint is_rational_c(self):
        return ZZX_deg(self.__numerator) == 0

    def trace(self, K=None):
        """
        Return the absolute or relative trace of this number field
        element.

        If K is given then K must be a subfield of the parent L of self, in
        which case the trace is the relative trace from L to K. In all
        other cases, the trace is the absolute trace down to QQ.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 -132/7*x^2 + x + 1); K
            Number Field in a with defining polynomial x^3 - 132/7*x^2 + x + 1
            sage: a.trace()
            132/7
            sage: (a+1).trace() == a.trace() + 3
            True

        If we are in an order, the trace is an integer::

            sage: K.<zeta> = CyclotomicField(17)
            sage: R = K.ring_of_integers()
            sage: R(zeta).trace().parent()
            Integer Ring

        TESTS::

            sage: F.<z> = CyclotomicField(5) ; t = 3*z**3 + 4*z**2 + 2
            sage: t.trace(F)
            3*z^3 + 4*z^2 + 2
        """
        if K is None:
            trace = self._pari_('x').trace()
            return QQ(trace) if self._parent.is_field() else ZZ(trace)
        return self.matrix(K).trace()

    def norm(self, K=None):
        """
        Return the absolute or relative norm of this number field element.

        If K is given then K must be a subfield of the parent L of self, in
        which case the norm is the relative norm from L to K. In all other
        cases, the norm is the absolute norm down to QQ.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + x^2 + x - 132/7); K
            Number Field in a with defining polynomial x^3 + x^2 + x - 132/7
            sage: a.norm()
            132/7
            sage: factor(a.norm())
            2^2 * 3 * 7^-1 * 11
            sage: K(0).norm()
            0

        Some complicated relatives norms in a tower of number fields.

        ::

            sage: K.<a,b,c> = NumberField([x^2 + 1, x^2 + 3, x^2 + 5])
            sage: L = K.base_field(); M = L.base_field()
            sage: a.norm()
            1
            sage: a.norm(L)
            1
            sage: a.norm(M)
            1
            sage: a
            a
            sage: (a+b+c).norm()
            121
            sage: (a+b+c).norm(L)
            2*c*b - 7
            sage: (a+b+c).norm(M)
            -11

        We illustrate that norm is compatible with towers::

            sage: z = (a+b+c).norm(L); z.norm(M)
            -11

        If we are in an order, the norm is an integer::

            sage: K.<a> = NumberField(x^3-2)
            sage: a.norm().parent()
            Rational Field
            sage: R = K.ring_of_integers()
            sage: R(a).norm().parent()
            Integer Ring

        TESTS::

            sage: F.<z> = CyclotomicField(5)
            sage: t = 3*z**3 + 4*z**2 + 2
            sage: t.norm(F)
            3*z^3 + 4*z^2 + 2
        """
        if K is None:
            norm = self._pari_('x').norm()
            return QQ(norm) if self._parent.is_field() else ZZ(norm)
        return self.matrix(K).determinant()

    def absolute_norm(self):
        """
        Return the absolute norm of this number field element.

        EXAMPLES::

            sage: K1.<a1> = CyclotomicField(11)
            sage: K2.<a2> = K1.extension(x^2 - 3)
            sage: K3.<a3> = K2.extension(x^2 + 1)
            sage: (a1 + a2 + a3).absolute_norm()
            1353244757701

            sage: QQ(7/5).absolute_norm()
            7/5
        """
        return self.norm()

    def relative_norm(self):
        """
        Return the relative norm of this number field element over the next field
        down in some tower of number fields.

        EXAMPLES::

            sage: K1.<a1> = CyclotomicField(11)
            sage: K2.<a2> = K1.extension(x^2 - 3)
            sage: (a1 + a2).relative_norm()
            a1^2 - 3
            sage: (a1 + a2).relative_norm().relative_norm() == (a1 + a2).absolute_norm()
            True

            sage: K.<x,y,z> = NumberField([x^2 + 1, x^3 - 3, x^2 - 5])
            sage: (x + y + z).relative_norm()
            y^2 + 2*z*y + 6
        """
        return self.norm(self.parent().base_field())

    def vector(self):
        """
        Return vector representation of self in terms of the basis for the
        ambient number field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: (2/3*a - 5/6).vector()
            (-5/6, 2/3)
            sage: (-5/6, 2/3)
            (-5/6, 2/3)
            sage: O = K.order(2*a)
            sage: (O.1).vector()
            (0, 2)
            sage: K.<a,b> = NumberField([x^2 + 1, x^2 - 3])
            sage: (a + b).vector()
            (b, 1)
            sage: O = K.order([a,b])
            sage: (O.1).vector()
            (-b, 1)
            sage: (O.2).vector()
            (1, -b)
        """
        return self.number_field().relative_vector_space()[2](self)

    def charpoly(self, var='x'):
        r"""
        Return the characteristic polynomial of this number field element.

        EXAMPLE::

            sage: K.<a> = NumberField(x^3 + 7)
            sage: a.charpoly()
            x^3 + 7
            sage: K(1).charpoly()
            x^3 - 3*x^2 + 3*x - 1
        """
        raise NotImplementedError, "Subclasses of NumberFieldElement must override charpoly()"

    def minpoly(self, var='x'):
        """
        Return the minimal polynomial of this number field element.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+3)
            sage: a.minpoly('x')
            x^2 + 3
            sage: R.<X> = K['X']
            sage: L.<b> = K.extension(X^2-(22 + a))
            sage: b.minpoly('t')
            t^2 - a - 22
            sage: b.absolute_minpoly('t')
            t^4 - 44*t^2 + 487
            sage: b^2 - (22+a)
            0
        """
        return self.charpoly(var).radical() # square free part of charpoly

    def is_integral(self):
        r"""
        Determine if a number is in the ring of integers of this number
        field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23)
            sage: a.is_integral()
            True
            sage: t = (1+a)/2
            sage: t.is_integral()
            True
            sage: t.minpoly()
            x^2 - x + 6
            sage: t = a/2
            sage: t.is_integral()
            False
            sage: t.minpoly()
            x^2 + 23/4

        An example in a relative extension::

            sage: K.<a,b> = NumberField([x^2+1, x^2+3])
            sage: (a+b).is_integral()
            True
            sage: ((a-b)/2).is_integral()
            False
        """
        return all([a in ZZ for a in self.absolute_minpoly()])

    def matrix(self, base=None):
        r"""
        If base is None, return the matrix of right multiplication by the
        element on the power basis `1, x, x^2, \ldots, x^{d-1}` for
        the number field. Thus the *rows* of this matrix give the images of
        each of the `x^i`.

        If base is not None, then base must be either a field that embeds
        in the parent of self or a morphism to the parent of self, in which
        case this function returns the matrix of multiplication by self on
        the power basis, where we view the parent field as a field over
        base.

        Specifying base as the base field over which the parent of self
        is a relative extension is equivalent to base being None

        INPUT:


        -  ``base`` - field or morphism


        EXAMPLES:

        Regular number field::

            sage: K.<a> = NumberField(QQ['x'].0^3 - 5)
            sage: M = a.matrix(); M
            [0 1 0]
            [0 0 1]
            [5 0 0]
            sage: M.base_ring() is QQ
            True

        Relative number field::

            sage: L.<b> = K.extension(K['x'].0^2 - 2)
            sage: M = b.matrix(); M
            [0 1]
            [2 0]
            sage: M.base_ring() is K
            True

        Absolute number field::

            sage: M = L.absolute_field('c').gen().matrix(); M
            [  0   1   0   0   0   0]
            [  0   0   1   0   0   0]
            [  0   0   0   1   0   0]
            [  0   0   0   0   1   0]
            [  0   0   0   0   0   1]
            [-17 -60 -12 -10   6   0]
            sage: M.base_ring() is QQ
            True

        More complicated relative number field::

            sage: L.<b> = K.extension(K['x'].0^2 - a); L
            Number Field in b with defining polynomial x^2 - a over its base field
            sage: M = b.matrix(); M
            [0 1]
            [a 0]
            sage: M.base_ring() is K
            True

        An example where we explicitly give the subfield or the embedding::

            sage: K.<a> = NumberField(x^4 + 1); L.<a2> = NumberField(x^2 + 1)
            sage: a.matrix(L)
            [ 0  1]
            [a2  0]

        Notice that if we compute all embeddings and choose a different
        one, then the matrix is changed as it should be::

            sage: v = L.embeddings(K)
            sage: a.matrix(v[1])
            [  0   1]
            [-a2   0]

        The norm is also changed::

            sage: a.norm(v[1])
            a2
            sage: a.norm(v[0])
            -a2

        TESTS::

            sage: F.<z> = CyclotomicField(5) ; t = 3*z**3 + 4*z**2 + 2
            sage: t.matrix(F)
            [3*z^3 + 4*z^2 + 2]
            sage: x=QQ['x'].gen()
            sage: K.<v>=NumberField(x^4 + 514*x^2 + 64321)
            sage: R.<r>=NumberField(x^2 + 4*v*x + 5*v^2 + 514)
            sage: r.matrix()
            [           0            1]
            [-5*v^2 - 514         -4*v]
            sage: r.matrix(K)
            [           0            1]
            [-5*v^2 - 514         -4*v]
            sage: r.matrix(R)
            [r]
            sage: foo=R.random_element()
            sage: foo.matrix(R) == matrix(1,1,[foo])
            True
        """
        from sage.matrix.constructor import matrix
        if base is self.parent():
            return matrix(1,1,[self])
        if base is not None and base is not self.base_ring():
            if number_field.is_NumberField(base):
                return self._matrix_over_base(base)
            else:
                return self._matrix_over_base_morphism(base)
        # Multiply each power of field generator on
        # the left by this element; make matrix
        # whose rows are the coefficients of the result,
        # and transpose.
        if self.__matrix is None:
            K = self.number_field()
            v = []
            x = K.gen()
            a = K(1)
            d = K.relative_degree()
            for n in range(d):
                v += (a*self).list()
                a *= x
            k = K.base_ring()
            import sage.matrix.matrix_space
            M = sage.matrix.matrix_space.MatrixSpace(k, d)
            self.__matrix = M(v)
        return self.__matrix

    def valuation(self, P):
        """
        Returns the valuation of self at a given prime ideal P.

        INPUT:


        -  ``P`` - a prime ideal of the parent of self


        .. note::

           The function ``ord()`` is an alias for ``valuation()``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: P = K.ideal(61).factor()[0][0]
            sage: b = a^2 + 30
            sage: b.valuation(P)
            1
            sage: b.ord(P)
            1
            sage: type(b.valuation(P))
            <type 'sage.rings.integer.Integer'>

        The function can be applied to elements in relative number fields::

            sage: L.<b> = K.extension(x^2 - 3)
            sage: [L(6).valuation(P) for P in L.primes_above(2)]
            [4]
            sage: [L(6).valuation(P) for P in L.primes_above(3)]
            [2, 2]
        """
        from number_field_ideal import is_NumberFieldIdeal
        from sage.rings.infinity import infinity
        if not is_NumberFieldIdeal(P):
            if is_NumberFieldElement(P):
                P = self.number_field().fractional_ideal(P)
            else:
                raise TypeError, "P must be an ideal"
        if not P.is_prime():
            raise ValueError, "P must be prime"
        if self == 0:
            return infinity
        return Integer_sage(self.number_field().pari_nf().elementval(self._pari_(), P.pari_prime()))

    ord = valuation

    def local_height(self, P, prec=None, weighted=False):
        r"""
        Returns the local height of self at a given prime ideal `P`.

        INPUT:


        -  ``P`` - a prime ideal of the parent of self

        - ``prec`` (int) -- desired floating point precision (defult:
          default RealField precision).

        - ``weighted`` (bool, default False) -- if True, apply local
          degree weighting.

        OUTPUT:

        (real) The local height of this number field element at the
        place `P`.  If ``weighted`` is True, this is multiplied by the
        local degree (as required for global heights).

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: P = K.ideal(61).factor()[0][0]
            sage: b = 1/(a^2 + 30)
            sage: b.local_height(P)
            4.11087386417331
            sage: b.local_height(P, weighted=True)
            8.22174772834662
            sage: b.local_height(P, 200)
            4.1108738641733112487513891034256147463156817430812610629374
            sage: (b^2).local_height(P)
            8.22174772834662
            sage: (b^-1).local_height(P)
            0.000000000000000

        A relative example::

            sage: PK.<y> = K[]
            sage: L.<c> = NumberField(y^2 + a)
            sage: L(1/4).local_height(L.ideal(2, c-a+1))
            1.38629436111989
        """
        if self.valuation(P) >= 0: ## includes the case self=0
            from sage.rings.real_mpfr import RealField
            if prec is None:
                return RealField().zero_element()
            else:
                return RealField(prec).zero_element()
        ht = self.abs_non_arch(P,prec).log()
        if not weighted:
            return ht
        nP = P.residue_class_degree()*P.absolute_ramification_index()
        return nP*ht

    def local_height_arch(self, i, prec=None, weighted=False):
        r"""
        Returns the local height of self at the `i`'th infinite place.

        INPUT:


        - ``i`` (int) - an integer in ``range(r+s)`` where `(r,s)` is the
           signature of the parent field (so `n=r+2s` is the degree).

        - ``prec`` (int) -- desired floating point precision (default:
          default RealField precision).

        - ``weighted`` (bool, default False) -- if True, apply local
          degree weighting, i.e. double the value for complex places.

        OUTPUT:

        (real) The archimedean local height of this number field
        element at the `i`'th infinite place.  If ``weighted`` is
        True, this is multiplied by the local degree (as required for
        global heights), i.e. 1 for real places and 2 for complex
        places.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: [p.codomain() for p in K.places()]
            [Real Field with 106 bits of precision,
            Real Field with 106 bits of precision,
            Complex Field with 53 bits of precision]
            sage: [a.local_height_arch(i) for i in range(3)]
            [0.5301924545717755083366563897519,
            0.5301924545717755083366563897519,
            0.886414217456333]
            sage: [a.local_height_arch(i, weighted=True) for i in range(3)]
            [0.5301924545717755083366563897519,
            0.5301924545717755083366563897519,
            1.77282843491267]

        A relative example::

            sage: L.<b, c> = NumberFieldTower([x^2 - 5, x^3 + x + 3])
            sage: [(b + c).local_height_arch(i) for i in range(4)]
            [1.238223390757884911842206617439,
            0.02240347229957875780769746914391,
            0.780028961749618,
            1.16048938497298]
        """
        K = self.number_field()
        emb = K.places(prec=prec)[i]
        a = emb(self).abs()
        Kv = emb.codomain()
        if a <= Kv.one_element():
            return Kv.zero_element()
        ht = a.log()
        from sage.rings.real_mpfr import is_RealField
        if weighted and not is_RealField(Kv):
            ht*=2
        return ht

    def global_height_non_arch(self, prec=None):
        """
        Returns the total non-archimedean component of the height of self.

        INPUT:

        - ``prec`` (int) -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        (real) The total non-archimedean component of the height of
        this number field element; that is, the sum of the local
        heights at all finite places, weighted by the local degrees.

        ALGORITHM:

        An alternative formula is `\log(d)` where `d` is the norm of
        the denominator ideal; this is used to avoid factorization.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: b = a/6
            sage: b.global_height_non_arch()
            7.16703787691222

        Check that this is equal to the sum of the non-archimedean
        local heights::

            sage: [b.local_height(P) for P in b.support()]
            [0.000000000000000, 0.693147180559945, 1.09861228866811, 1.09861228866811]
            sage: [b.local_height(P, weighted=True) for P in b.support()]
            [0.000000000000000, 2.77258872223978, 2.19722457733622, 2.19722457733622]
            sage: sum([b.local_height(P,weighted=True) for P in b.support()])
            7.16703787691222

        A relative example::

            sage: PK.<y> = K[]
            sage: L.<c> = NumberField(y^2 + a)
            sage: (c/10).global_height_non_arch()
            18.4206807439524
        """
        from sage.rings.real_mpfr import RealField
        if prec is None:
            R = RealField()
        else:
            R = RealField(prec)
        if self.is_zero():
            return R.zero_element()
        return R(self.denominator_ideal().absolute_norm()).log()

    def global_height_arch(self, prec=None):
        """
        Returns the total archimedean component of the height of self.

        INPUT:

        - ``prec`` (int) -- desired floating point precision (defult:
          default RealField precision).

        OUTPUT:

        (real) The total archimedean component of the height of
        this number field element; that is, the sum of the local
        heights at all infinite places.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: b = a/2
            sage: b.global_height_arch()
            0.38653407379277...
        """
        r,s = self.number_field().signature()
        hts = [self.local_height_arch(i, prec, weighted=True) for i in range(r+s)]
        return sum(hts, hts[0].parent().zero_element())

    def global_height(self, prec=None):
        """
        Returns the absolute logarithmic height of this number field element.

        INPUT:

        - ``prec`` (int) -- desired floating point precision (defult:
          default RealField precision).

        OUTPUT:

        (real) The absolute logarithmic height of this number field
        element; that is, the sum of the local heights at all finite
        and infinite places, scaled by the degree to make the result independent of
        the parent field.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: b = a/2
            sage: b.global_height()
            0.789780699008...
            sage: b.global_height(prec=200)
            0.78978069900813892060267152032141577237037181070060784564457

        The global height of an algebraic number is absolute, i.e. it
        does not depend on the parent field::

            sage: QQ(6).global_height()
            1.79175946922805
            sage: K(6).global_height()
            1.79175946922805

            sage: L.<b> = NumberField((a^2).minpoly())
            sage: L.degree()
            2
            sage: b.global_height() # element of L (degree 2 field)
            1.41660667202811
            sage: (a^2).global_height() # element of K (degree 4 field)
            1.41660667202811

        And of course every element has the same height as it's inverse::

            sage: K.<s> = QuadraticField(2)
            sage: s.global_height()
            0.346573590279973
            sage: (1/s).global_height()   #make sure that 11758 is fixed
            0.346573590279973

        """
        return (self.global_height_non_arch(prec)+self.global_height_arch(prec))/self.number_field().absolute_degree()

    def numerator_ideal(self):
        """
        Return the numerator ideal of this number field element.

        .. note::

           A ValueError will be raised if this function is called on
           0.

        .. seealso::

           :meth:`~denominator_ideal`

        OUTPUT:

        (integral ideal) The numerator ideal `N` of this element,
        where for a non-zero number field element `a`, the principal
        ideal generated by `a` has the form `N/D` where `N` and `D`
        are coprime integral ideals.  An error is raised if the
        element is zero.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+5)
            sage: b = (1+a)/2
            sage: b.norm()
            3/2
            sage: N = b.numerator_ideal(); N
            Fractional ideal (3, a + 1)
            sage: N.norm()
            3
            sage: (1/b).numerator_ideal()
            Fractional ideal (2, a + 1)

        TESTS:

        Undefined for 0::

            sage: K(0).numerator_ideal()
            Traceback (most recent call last):
            ...
            ValueError: numerator ideal of 0 is not defined.
        """
        if self.is_zero():
            raise ValueError, "numerator ideal of 0 is not defined."
        return self.number_field().ideal(self).numerator()

    def denominator_ideal(self):
        """
        Return the denominator ideal of this number field element.

        .. note::

           A ValueError will be raised if this function is called on
           0.

        .. seealso::

           :meth:`~numerator_ideal`

        OUTPUT:

        (integral ideal) The denominator ideal `D` of this element,
        where for a non-zero number field element `a`, the principal
        ideal generated by `a` has the form `N/D` where `N` and `D`
        are coprime integral ideals.  An error is raised if the
        element is zero.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+5)
            sage: b = (1+a)/2
            sage: b.norm()
            3/2
            sage: D = b.denominator_ideal(); D
            Fractional ideal (2, a + 1)
            sage: D.norm()
            2
            sage: (1/b).denominator_ideal()
            Fractional ideal (3, a + 1)

        TESTS:

        Undefined for 0::

            sage: K(0).denominator_ideal()
            Traceback (most recent call last):
            ...
            ValueError: denominator ideal of 0 is not defined.
        """
        if self.is_zero():
            raise ValueError, "denominator ideal of 0 is not defined."
        return self.number_field().ideal(self).denominator()

    def support(self):
        """
        Return the support of this number field element.

        OUTPUT: A sorted list of the primes ideals at which this number
        field element has nonzero valuation. An error is raised if the
        element is zero.

        EXAMPLES::

            sage: x = ZZ['x'].gen()
            sage: F.<t> = NumberField(x^3 - 2)

        ::

            sage: P5s = F(5).support()
            sage: P5s
            [Fractional ideal (-t^2 - 1), Fractional ideal (t^2 - 2*t - 1)]
            sage: all(5 in P5 for P5 in P5s)
            True
            sage: all(P5.is_prime() for P5 in P5s)
            True
            sage: [ P5.norm() for P5 in P5s ]
            [5, 25]

        TESTS:

        It doesn't make sense to factor the ideal (0)::

            sage: F(0).support()
            Traceback (most recent call last):
            ...
            ArithmeticError: Support of 0 is not defined.
        """
        if self.is_zero():
            raise ArithmeticError, "Support of 0 is not defined."
        return self.number_field().primes_above(self)

    def _matrix_over_base(self, L):
        """
        Return the matrix of self over the base field L.

        EXAMPLES::

            sage: K.<a> = NumberField(ZZ['x'].0^3-2, 'a')
            sage: L.<b> = K.extension(ZZ['x'].0^2+3, 'b')
            sage: L(a)._matrix_over_base(K) == L(a).matrix()
            True
        """
        K = self.number_field()
        E = L.embeddings(K)
        if len(E) == 0:
            raise ValueError, "no way to embed L into parent's base ring K"
        phi = E[0]
        return self._matrix_over_base_morphism(phi)

    def _matrix_over_base_morphism(self, phi):
        """
        Return the matrix of self over a specified base, where phi gives a
        map from the specified base to self.parent().

        EXAMPLES::

            sage: F.<alpha> = NumberField(ZZ['x'].0^5-2)
            sage: h = Hom(QQ,F)([1])
            sage: alpha._matrix_over_base_morphism(h) == alpha.matrix()
            True
            sage: alpha._matrix_over_base_morphism(h) == alpha.matrix(QQ)
            True
        """
        L = phi.domain()

        ## the code below doesn't work if the morphism is
        ## over QQ, since QQ.primitive_element() doesn't
        ## make sense
        if L is QQ:
            K = phi.codomain()
            if K != self.number_field():
                raise ValueError, "codomain of phi must be parent of self"
            ## the variable name is irrelevant below, because the
            ## matrix is over QQ
            F = K.absolute_field('alpha')
            from_f, to_F = F.structure()
            return to_F(self).matrix()

        alpha = L.primitive_element()
        beta = phi(alpha)
        K = phi.codomain()
        if K != self.number_field():
            raise ValueError, "codomain of phi must be parent of self"

        # Construct a relative extension over L (= QQ(beta))
        M = K.relativize(beta, ('a','b'))
                     # variable name a is OK, since this is temporary

        # Carry self over to M.
        from_M, to_M = M.structure()
        try:
            z = to_M(self)
        except Exception:
            return to_M, self, K, beta

        # Compute the relative matrix of self, but in M
        R = z.matrix()

        # Map back to L.
        psi = M.base_field().hom([alpha])
        return R.apply_morphism(psi)


    def list(self):
        """
        Return the list of coefficients of self written in terms of a power
        basis.

        EXAMPLE::

            sage: K.<a> = NumberField(x^3 - x + 2); ((a + 1)/(a + 2)).list()
            [1/4, 1/2, -1/4]
            sage: K.<a, b> = NumberField([x^3 - x + 2, x^2 + 23]); ((a + b)/(a + 2)).list()
            [3/4*b - 1/2, -1/2*b + 1, 1/4*b - 1/2]
        """
        raise NotImplementedError

    def inverse_mod(self, I):
        """
        Returns the inverse of self mod the integral ideal I.

        INPUT:

        -  ``I`` - may be an ideal of self.parent(), or an element or list
           of elements of self.parent() generating a nonzero ideal. A ValueError
           is raised if I is non-integral or zero. A ZeroDivisionError is
           raised if I + (x) != (1).

        NOTE: It's not implemented yet for non-integral elements.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 23)
            sage: N = k.ideal(3)
            sage: d = 3*a + 1
            sage: d.inverse_mod(N)
            1

        ::

            sage: k.<a> = NumberField(x^3 + 11)
            sage: d = a + 13
            sage: d.inverse_mod(a^2)*d - 1 in k.ideal(a^2)
            True
            sage: d.inverse_mod((5, a + 1))*d - 1 in k.ideal(5, a + 1)
            True
            sage: K.<b> = k.extension(x^2 + 3)
            sage: b.inverse_mod([37, a - b])
            7
            sage: 7*b - 1 in K.ideal(37, a - b)
            True
            sage: b.inverse_mod([37, a - b]).parent() == K
            True
        """
        R = self.number_field().ring_of_integers()
        try:
            return _inverse_mod_generic(R(self), I)
        except TypeError: # raised by failure of R(self)
            raise NotImplementedError, "inverse_mod is not implemented for non-integral elements"


    def residue_symbol(self, P, m, check=True):
        r"""
        The m-th power residue symbol for an element self and proper ideal P.

        .. math:: \left(\frac{\alpha}{\mathbf{P}}\right) \equiv \alpha^{\frac{N(\mathbf{P})-1}{m}} \operatorname{mod} \mathbf{P}

        .. note:: accepts m=1, in which case returns 1

        .. note:: can also be called for an ideal from sage.rings.number_field_ideal.residue_symbol

        .. note:: self is coerced into the number field of the ideal P

        .. note:: if m=2, self is an integer, and P is an ideal of a number field of absolute degree 1 (i.e. it is a copy of the rationals), then this calls kronecker_symbol, which is implemented using GMP.

        INPUT:

        - ``P`` - proper ideal of the number field (or an extension)

        - ``m`` - positive integer

        OUTPUT:

        - an m-th root of unity in the number field

        EXAMPLES:

        Quadratic Residue (11 is not a square modulo 17)::

            sage: K.<a> = NumberField(x - 1)
            sage: K(11).residue_symbol(K.ideal(17),2)
            -1
            sage: kronecker_symbol(11,17)
            -1

        The result depends on the number field of the ideal::

            sage: K.<a> = NumberField(x - 1)
            sage: L.<b> = K.extension(x^2 + 1)
            sage: K(7).residue_symbol(K.ideal(11),2)
            -1
            sage: K(7).residue_symbol(L.ideal(11),2)
            1

        Cubic Residue::

            sage: K.<w> = NumberField(x^2 - x + 1)
            sage: (w^2 + 3).residue_symbol(K.ideal(17),3)
            -w

        The field must contain the m-th roots of unity::

            sage: K.<w> = NumberField(x^2 - x + 1)
            sage: (w^2 + 3).residue_symbol(K.ideal(17),5)
            Traceback (most recent call last):
            ...
            ValueError: The residue symbol to that power is not defined for the number field

        """
        return P.residue_symbol(self,m,check)



cdef class NumberFieldElement_absolute(NumberFieldElement):

    def _magma_init_(self, magma):
        """
        Return Magma version of this number field element.

        INPUT:


        -  ``magma`` - a Magma interpreter


        OUTPUT: MagmaElement that has parent the Magma object corresponding
        to the parent number field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + 2)
            sage: a._magma_init_(magma)            # optional - magma
            '(_sage_[...]![0, 1, 0])'
            sage: m = magma((2/3)*a^2 - 17/3); m   # optional - magma
            1/3*(2*a^2 - 17)
            sage: m.sage()                         # optional - magma
            2/3*a^2 - 17/3

        An element of a cyclotomic field.

        ::

            sage: K = CyclotomicField(9)
            sage: K.gen()
            zeta9
            sage: K.gen()._magma_init_(magma)     # optional - magma
            '(_sage_[...]![0, 1, 0, 0, 0, 0])'
            sage: magma(K.gen())                  # optional - magma
            zeta9
            sage: _.sage()                        # optional - magma
            zeta9
        """
        K = magma(self.parent())
        return '(%s!%s)'%(K.name(), self.list())

    def absolute_charpoly(self, var='x', algorithm=None):
        r"""
        Return the characteristic polynomial of this element over `\QQ`.

        For the meaning of the optional argument ``algorithm``, see :meth:`charpoly`.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: K.<a> = NumberField(x^4 + 2, 'a')
            sage: a.absolute_charpoly()
            x^4 + 2
            sage: a.absolute_charpoly('y')
            y^4 + 2
            sage: (-a^2).absolute_charpoly()
            x^4 + 4*x^2 + 4
            sage: (-a^2).absolute_minpoly()
            x^2 + 2

            sage: a.absolute_charpoly(algorithm='pari') == a.absolute_charpoly(algorithm='sage')
            True
        """
        # this hack is necessary because quadratic fields override
        # charpoly(), and they don't take the argument 'algorithm'
        if algorithm is None:
            return self.charpoly(var)
        return self.charpoly(var, algorithm)

    def absolute_minpoly(self, var='x', algorithm=None):
        r"""
        Return the minimal polynomial of this element over
        `\QQ`.

        For the meaning of the optional argument algorithm, see :meth:`charpoly`.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: f = x^10 - 5*x^9 + 15*x^8 - 68*x^7 + 81*x^6 - 221*x^5 + 141*x^4 - 242*x^3 - 13*x^2 - 33*x - 135
            sage: K.<a> = NumberField(f, 'a')
            sage: a.absolute_charpoly()
            x^10 - 5*x^9 + 15*x^8 - 68*x^7 + 81*x^6 - 221*x^5 + 141*x^4 - 242*x^3 - 13*x^2 - 33*x - 135
            sage: a.absolute_charpoly('y')
            y^10 - 5*y^9 + 15*y^8 - 68*y^7 + 81*y^6 - 221*y^5 + 141*y^4 - 242*y^3 - 13*y^2 - 33*y - 135
            sage: b = -79/9995*a^9 + 52/9995*a^8 + 271/9995*a^7 + 1663/9995*a^6 + 13204/9995*a^5 + 5573/9995*a^4 + 8435/1999*a^3 - 3116/9995*a^2 + 7734/1999*a + 1620/1999
            sage: b.absolute_charpoly()
            x^10 + 10*x^9 + 25*x^8 - 80*x^7 - 438*x^6 + 80*x^5 + 2950*x^4 + 1520*x^3 - 10439*x^2 - 5130*x + 18225
            sage: b.absolute_minpoly()
            x^5 + 5*x^4 - 40*x^2 - 19*x + 135

            sage: b.absolute_minpoly(algorithm='pari') == b.absolute_minpoly(algorithm='sage')
            True
        """
        # this hack is necessary because quadratic fields override
        # minpoly(), and they don't take the argument 'algorithm'
        if algorithm is None:
            return self.minpoly(var)
        return self.minpoly(var, algorithm)

    def charpoly(self, var='x', algorithm=None):
        r"""
        The characteristic polynomial of this element, over
        `\QQ` if self is an element of a field, and over
        `\ZZ` is self is an element of an order.

        This is the same as ``self.absolute_charpoly`` since
        this is an element of an absolute extension.

        The optional argument algorithm controls how the
        characteristic polynomial is computed: 'pari' uses PARI,
        'sage' uses charpoly for Sage matrices.  The default value
        None means that 'pari' is used for small degrees (up to the
        value of the constant TUNE_CHARPOLY_NF, currently at 25),
        otherwise 'sage' is used.  The constant TUNE_CHARPOLY_NF
        should give reasonable performance on all architectures;
        however, if you feel the need to customize it to your own
        machine, see trac ticket 5213 for a tuning script.

        EXAMPLES:

        We compute the characteristic polynomial of the cube root of `2`.

        ::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-2)
            sage: a.charpoly('x')
            x^3 - 2
            sage: a.charpoly('y').parent()
            Univariate Polynomial Ring in y over Rational Field

        TESTS::

            sage: R = K.ring_of_integers()
            sage: R(a).charpoly()
            x^3 - 2
            sage: R(a).charpoly().parent()
            Univariate Polynomial Ring in x over Integer Ring

            sage: R(a).charpoly(algorithm='pari') == R(a).charpoly(algorithm='sage')
            True
        """
        if algorithm is None:
            if self._parent.degree() <= TUNE_CHARPOLY_NF:
                algorithm = 'pari'
            else:
                algorithm = 'sage'
        R = self._parent.base_ring()[var]
        if algorithm == 'pari':
            return R(self._pari_('x').charpoly())
        if algorithm == 'sage':
            return R(self.matrix().charpoly())

    def minpoly(self, var='x', algorithm=None):
        """
        Return the minimal polynomial of this number field element.

        For the meaning of the optional argument algorithm, see charpoly().

        EXAMPLES:

        We compute the characteristic polynomial of cube root of `2`.

        ::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-2)
            sage: a.minpoly('x')
            x^3 - 2
            sage: a.minpoly('y').parent()
            Univariate Polynomial Ring in y over Rational Field

        TESTS::

            sage: R = K.ring_of_integers()
            sage: R(a).minpoly()
            x^3 - 2
            sage: R(a).minpoly().parent()
            Univariate Polynomial Ring in x over Integer Ring

            sage: R(a).minpoly(algorithm='pari') == R(a).minpoly(algorithm='sage')
            True

        """
        return self.charpoly(var, algorithm).radical() # square free part of charpoly

    def list(self):
        """
        Return the list of coefficients of self written in terms of a power
        basis.

        EXAMPLE::

            sage: K.<z> = CyclotomicField(3)
            sage: (2+3/5*z).list()
            [2, 3/5]
            sage: (5*z).list()
            [0, 5]
            sage: K(3).list()
            [3, 0]
        """
        n = self.number_field().degree()
        v = self._coefficients()
        z = sage.rings.rational.Rational(0)
        return v + [z]*(n - len(v))

    def lift(self, var='x'):
        """
        Return an element of QQ[x], where this number field element
        lives in QQ[x]/(f(x)).

        EXAMPLES::

            sage: K.<a> = QuadraticField(-3)
            sage: a.lift()
            x

        """
        R = self.number_field().base_field()[var]
        return R(self.list())

    def is_real_positive(self, min_prec=53):
        r"""
        Using the ``n`` method of approximation, return ``True`` if
        ``self`` is a real positive number and ``False`` otherwise.
        This method is completely dependent of the embedding used by
        the ``n`` method.

        The algorithm first checks that ``self`` is not a strictly
        complex number. Then if ``self`` is not zero, by approximation
        more and more precise, the method answers True if the
        number is positive. Using `RealInterval`, the result is
        guaranteed to be correct.

        For CyclotomicField, the embedding is the natural one
        sending `zetan` on `cos(2*pi/n)`.

        EXAMPLES::

            sage: K.<a> = CyclotomicField(3)
            sage: (a+a^2).is_real_positive()
            False
            sage: (-a-a^2).is_real_positive()
            True
            sage: K.<a> = CyclotomicField(1000)
            sage: (a+a^(-1)).is_real_positive()
            True
            sage: K.<a> = CyclotomicField(1009)
            sage: d = a^252
            sage: (d+d.conjugate()).is_real_positive()
            True
            sage: d = a^253
            sage: (d+d.conjugate()).is_real_positive()
            False
            sage: K.<a> = QuadraticField(3)
            sage: a.is_real_positive()
            True
            sage: K.<a> = QuadraticField(-3)
            sage: a.is_real_positive()
            False
            sage: (a-a).is_real_positive()
            False
        """
        if self != self.conjugate() or self.is_zero():
            return False
        else:
            approx = RealInterval(self.n(min_prec).real())
            if approx.lower() > 0:
                return True
            else:
                if approx.upper() < 0:
                    return False
                else:
                    return self.is_real_positive(min_prec+20)

cdef class NumberFieldElement_relative(NumberFieldElement):
    r"""
    The current relative number field element implementation
    does everything in terms of absolute polynomials.

    All conversions from relative polynomials, lists, vectors, etc
    should happen in the parent.
    """
    def __init__(self, parent, f):
        r"""
        EXAMPLE::

            sage: L.<a, b> = NumberField([x^2 + 1, x^2 + 2])
            sage: type(a) # indirect doctest
            <type 'sage.rings.number_field.number_field_element.NumberFieldElement_relative'>

        We can create relative number field elements from PARI::

            sage: y = polygen(QQ)
            sage: K.<a> = NumberField(y^2 + y + 1)
            sage: x = polygen(K)
            sage: L.<b> = NumberField(x^4 + a*x + 2)
            sage: e = pari(a*b); e
            Mod(-y^4 - 2, y^8 - y^5 + 4*y^4 + y^2 - 2*y + 4)
            sage: L(e)  # Conversion from PARI absolute number field element
            a*b
            sage: e = L.pari_rnf().rnfeltabstorel(e); e
            Mod(Mod(y, y^2 + y + 1)*x, x^4 + y*x + 2)
            sage: L(e)  # Conversion from PARI relative number field element
            a*b
            sage: e = pari('Mod(0, x^8 + 1)'); L(e)  # Wrong modulus
            Traceback (most recent call last):
            ...
            TypeError: Coercion of PARI polmod with modulus x^8 + 1 into number field with defining polynomial x^8 - x^5 + 4*x^4 + x^2 - 2*x + 4 failed

        We test a relative number field element created "by hand"::

            sage: e = pari("Mod(Mod(y, y^2 + y + 1)*x^2 + Mod(1, y^2 + y + 1), x^4 + y*x + 2)")
            sage: L(e)
            a*b^2 + 1

        Currently, conversions of PARI relative number fields are not checked::

            sage: e = pari('Mod(y*x, x^4 + y^2*x + 2)'); L(e)  # Wrong modulus, but succeeds anyway
            a*b
        """
        NumberFieldElement.__init__(self, parent, f)

    def __getitem__(self, n):
        """
        Return the n-th coefficient of this relative number field element, written
        as a polynomial in the generator.

        Note that `n` must be between 0 and `d-1`, where
        `d` is the relative degree of the number field.

        EXAMPLES::

            sage: K.<a, b> = NumberField([x^3 - 5, x^2 + 3])
            sage: c = (a + b)^3; c
            3*b*a^2 - 9*a - 3*b + 5
            sage: c[0]
            -3*b + 5

        We illustrate bounds checking::

            sage: c[-1]
            Traceback (most recent call last):
            ...
            IndexError: index must be between 0 and the relative degree minus 1.
            sage: c[4]
            Traceback (most recent call last):
            ...
            IndexError: index must be between 0 and the relative degree minus 1.

        The list method implicitly calls ``__getitem__``::

            sage: list(c)
            [-3*b + 5, -9, 3*b]
            sage: K(list(c)) == c
            True
        """
        if n < 0 or n >= self.parent().relative_degree():
            raise IndexError, "index must be between 0 and the relative degree minus 1."
        return self.vector()[n]

    def _magma_init_(self, magma):
        """
        EXAMPLES::

            sage: K.<a, b> = NumberField([x^3 - 5, x^2 + 3])
            sage: a._magma_init_(magma)
            Traceback (most recent call last):
            ...
            TypeError: coercion of relative number field elements to Magma is not implemented
        """
        raise TypeError, "coercion of relative number field elements to Magma is not implemented"

    def list(self):
        """
        Return the list of coefficients of self written in terms of a power
        basis.

        EXAMPLES::

            sage: K.<a,b> = NumberField([x^3+2, x^2+1])
            sage: a.list()
            [0, 1, 0]
            sage: v = (K.base_field().0 + a)^2 ; v
            a^2 + 2*b*a - 1
            sage: v.list()
            [-1, 2*b, 1]
        """
        return self.vector().list()

    def lift(self, var='x'):
        """
        Return an element of K[x], where this number field element
        lives in the relative number field K[x]/(f(x)).

        EXAMPLES::

            sage: K.<a> = QuadraticField(-3)
            sage: x = polygen(K)
            sage: L.<b> = K.extension(x^7 + 5)
            sage: u = L(1/2*a + 1/2 + b + (a-9)*b^5)
            sage: u.lift()
            (a - 9)*x^5 + x + 1/2*a + 1/2

        """
        K = self.number_field()
        # Compute representation of self in terms of relative vector space.
        R = K.base_field()[var]
        return R(self.list())

    def _repr_(self):
        r"""
        EXAMPLE::

            sage: L.<a, b> = NumberField([x^3 - x + 1, x^2 + 23])
            sage: repr(a^4*b) # indirect doctest
            'b*a^2 - b*a'
        """
        K = self.number_field()
        # Compute representation of self in terms of relative vector space.
        R = K.base_field()[K.variable_name()]
        return repr(R(self.list()))

    def _latex_(self):
        r"""
        Returns the latex representation for this element.

        EXAMPLES::

            sage: C.<zeta> = CyclotomicField(12)
            sage: PC.<x> = PolynomialRing(C)
            sage: K.<alpha> = NumberField(x^2 - 7)
            sage: latex((alpha + zeta)^4) # indirect doctest
            \left(4 \zeta_{12}^{3} + 28 \zeta_{12}\right) \alpha + 43 \zeta_{12}^{2} + 48
            sage: PK.<y> = PolynomialRing(K)
            sage: L.<beta> = NumberField(y^3 + y + alpha)
            sage: latex((beta + zeta)^3) # indirect doctest
            3 \zeta_{12} \beta^{2} + \left(3 \zeta_{12}^{2} - 1\right) \beta - \alpha + \zeta_{12}^{3}
        """
        K = self.number_field()
        R = K.base_field()[K.variable_name()]
        return R(self.list())._latex_()

    def charpoly(self, var='x'):
        r"""
        The characteristic polynomial of this element over its base field.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: K.<a, b> = QQ.extension([x^2 + 2, x^5 + 400*x^4 + 11*x^2 + 2])
            sage: a.charpoly()
            x^2 + 2
            sage: b.charpoly()
            x^2 - 2*b*x + b^2
            sage: b.minpoly()
            x - b

            sage: K.<a, b> = NumberField([x^2 + 2, x^2 + 1000*x + 1])
            sage: y = K['y'].0
            sage: L.<c> = K.extension(y^2 + a*y + b)
            sage: c.charpoly()
            x^2 + a*x + b
            sage: c.minpoly()
            x^2 + a*x + b
            sage: L(a).charpoly()
            x^2 - 2*a*x - 2
            sage: L(a).minpoly()
            x - a
            sage: L(b).charpoly()
            x^2 - 2*b*x - 1000*b - 1
            sage: L(b).minpoly()
            x - b
        """
        R = self._parent.base_ring()[var]
        return R(self.matrix().charpoly())

    def absolute_charpoly(self, var='x', algorithm=None):
        r"""
        The characteristic polynomial of this element over
        `\QQ`.

        We construct a relative extension and find the characteristic
        polynomial over `\QQ`.

        The optional argument algorithm controls how the
        characteristic polynomial is computed: 'pari' uses PARI,
        'sage' uses charpoly for Sage matrices.  The default value
        None means that 'pari' is used for small degrees (up to the
        value of the constant TUNE_CHARPOLY_NF, currently at 25),
        otherwise 'sage' is used.  The constant TUNE_CHARPOLY_NF
        should give reasonable performance on all architectures;
        however, if you feel the need to customize it to your own
        machine, see trac ticket 5213 for a tuning script.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-2)
            sage: S.<X> = K[]
            sage: L.<b> = NumberField(X^3 + 17); L
            Number Field in b with defining polynomial X^3 + 17 over its base field
            sage: b.absolute_charpoly()
            x^9 + 51*x^6 + 867*x^3 + 4913
            sage: b.charpoly()(b)
            0
            sage: a = L.0; a
            b
            sage: a.absolute_charpoly('x')
            x^9 + 51*x^6 + 867*x^3 + 4913
            sage: a.absolute_charpoly('y')
            y^9 + 51*y^6 + 867*y^3 + 4913

            sage: a.absolute_charpoly(algorithm='pari') == a.absolute_charpoly(algorithm='sage')
            True
        """
        if algorithm is None:
            # this might not be the optimal condition; maybe it should
            # be .degree() instead of .absolute_degree()
            # there are too many bugs in relative number fields to
            # figure this out now
            if self._parent.absolute_degree() <= TUNE_CHARPOLY_NF:
                algorithm = 'pari'
            else:
                algorithm = 'sage'
        R = QQ[var]
        if algorithm == 'pari':
            return R(self._pari_().charpoly())
        if algorithm == 'sage':
            return R(self.matrix(QQ).charpoly())

    def absolute_minpoly(self, var='x', algorithm=None):
        r"""
        Return the minimal polynomial over `\QQ` of this element.

        For the meaning of the optional argument ``algorithm``, see :meth:`absolute_charpoly`.

        EXAMPLES::

            sage: K.<a, b> = NumberField([x^2 + 2, x^2 + 1000*x + 1])
            sage: y = K['y'].0
            sage: L.<c> = K.extension(y^2 + a*y + b)
            sage: c.absolute_charpoly()
            x^8 - 1996*x^6 + 996006*x^4 + 1997996*x^2 + 1
            sage: c.absolute_minpoly()
            x^8 - 1996*x^6 + 996006*x^4 + 1997996*x^2 + 1
            sage: L(a).absolute_charpoly()
            x^8 + 8*x^6 + 24*x^4 + 32*x^2 + 16
            sage: L(a).absolute_minpoly()
            x^2 + 2
            sage: L(b).absolute_charpoly()
            x^8 + 4000*x^7 + 6000004*x^6 + 4000012000*x^5 + 1000012000006*x^4 + 4000012000*x^3 + 6000004*x^2 + 4000*x + 1
            sage: L(b).absolute_minpoly()
            x^2 + 1000*x + 1
        """
        return self.absolute_charpoly(var, algorithm).radical()

    def valuation(self, P):
        """
        Returns the valuation of self at a given prime ideal P.

        INPUT:


        -  ``P`` - a prime ideal of relative number field which is the parent of self


        EXAMPLES::

            sage: K.<a, b, c> = NumberField([x^2 - 2, x^2 - 3, x^2 - 5])
            sage: P = K.prime_factors(5)[0]
            sage: (2*a + b - c).valuation(P)
            1
        """
        P_abs = P.absolute_ideal()
        abs = P_abs.number_field()
        to_abs = abs.structure()[1]
        return to_abs(self).valuation(P_abs)


cdef class OrderElement_absolute(NumberFieldElement_absolute):
    """
    Element of an order in an absolute number field.

    EXAMPLES::

        sage: K.<a> = NumberField(x^2 + 1)
        sage: O2 = K.order(2*a)
        sage: w = O2.1; w
        2*a
        sage: parent(w)
        Order in Number Field in a with defining polynomial x^2 + 1

        sage: w.absolute_charpoly()
        x^2 + 4
        sage: w.absolute_charpoly().parent()
        Univariate Polynomial Ring in x over Integer Ring
        sage: w.absolute_minpoly()
        x^2 + 4
        sage: w.absolute_minpoly().parent()
        Univariate Polynomial Ring in x over Integer Ring
    """
    def __init__(self, order, f):
        r"""
        EXAMPLE::

            sage: K.<a> = NumberField(x^3 + 2)
            sage: O2 = K.order(2*a)
            sage: type(O2.1) # indirect doctest
            <type 'sage.rings.number_field.number_field_element.OrderElement_absolute'>
        """
        K = order.number_field()
        NumberFieldElement_absolute.__init__(self, K, f)
        self._number_field = K
        (<Element>self)._parent = order

    cdef _new(self):
        """
        Quickly creates a new initialized NumberFieldElement with the same
        parent as self.

        EXAMPLES:

        This is called implicitly in multiplication::

            sage: O = EquationOrder(x^3 + 18, 'a')
            sage: O.1 * O.1 * O.1
            -18
        """
        cdef OrderElement_absolute x
        x = <OrderElement_absolute>PY_NEW_SAME_TYPE(self)
        x._parent = self._parent
        x._number_field = self._parent.number_field()
        x.__fld_numerator = self.__fld_numerator
        x.__fld_denominator = self.__fld_denominator
        return x

    cdef number_field(self):
        r"""
        Return the number field of self. Only accessible from Cython.

        EXAMPLE::

            sage: K = NumberField(x^3 - 17, 'a')
            sage: OK = K.ring_of_integers()
            sage: a = OK(K.gen())
            sage: a._number_field() is K # indirect doctest
            True
        """
        return self._number_field

    cpdef RingElement _div_(self, RingElement other):
        r"""
        Implement division, checking that the result has the right parent.
        It's not so crucial what the parent actually is, but it is crucial
        that the returned value really is an element of its supposed
        parent! This fixes trac #4190.

        EXAMPLES::

            sage: K = NumberField(x^3 - 17, 'a')
            sage: OK = K.ring_of_integers()
            sage: a = OK(K.gen())
            sage: (17/a) in OK # indirect doctest
            True
            sage: (17/a).parent() is K # indirect doctest
            True
            sage: (17/(2*a)).parent() is K # indirect doctest
            True
            sage: (17/(2*a)) in OK # indirect doctest
            False
        """
        cdef NumberFieldElement_absolute x
        x = NumberFieldElement_absolute._div_(self, other)
        return self._parent.number_field()(x)

    def inverse_mod(self, I):
        r"""
        Return an inverse of self modulo the given ideal.

        INPUT:


        -  ``I`` - may be an ideal of self.parent(), or an
           element or list of elements of self.parent() generating a nonzero
           ideal. A ValueError is raised if I is non-integral or is zero.
           A ZeroDivisionError is raised if I + (x) != (1).


        EXAMPLES::

            sage: OE = NumberField(x^3 - x + 2, 'w').ring_of_integers()
            sage: w = OE.ring_generators()[0]
            sage: w.inverse_mod(13*OE)
            6*w^2 - 6
            sage: w * (w.inverse_mod(13)) - 1 in 13*OE
            True
            sage: w.inverse_mod(13).parent() == OE
            True
            sage: w.inverse_mod(2*OE)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: w is not invertible modulo Fractional ideal (2)
        """
        R = self.parent()
        return R(_inverse_mod_generic(self, I))

    def __invert__(self):
        r"""
        Implement inversion, checking that the return value has the right
        parent. See trac #4190.

        EXAMPLE::

            sage: K = NumberField(x^3 -x + 2, 'a')
            sage: OK = K.ring_of_integers()
            sage: a = OK(K.gen())
            sage: (~a).parent() is K
            True
            sage: (~a) in OK
            False
            sage: a**(-1) in OK
            False
        """
        return self._parent.number_field()(NumberFieldElement_absolute.__invert__(self))

cdef class OrderElement_relative(NumberFieldElement_relative):
    """
    Element of an order in a relative number field.

    EXAMPLES::

        sage: O = EquationOrder([x^2 + x + 1, x^3 - 2],'a,b')
        sage: c = O.1; c
        (-2*b^2 - 2)*a - 2*b^2 - b
        sage: type(c)
        <type 'sage.rings.number_field.number_field_element.OrderElement_relative'>
    """
    def __init__(self, order, f):
        r"""
        EXAMPLE::

            sage: O = EquationOrder([x^2 + x + 1, x^3 - 2],'a,b')
            sage: type(O.1) # indirect doctest
            <type 'sage.rings.number_field.number_field_element.OrderElement_relative'>
        """
        K = order.number_field()
        NumberFieldElement_relative.__init__(self, K, f)
        (<Element>self)._parent = order
        self._number_field = K

    cdef number_field(self):
        return self._number_field

    cdef _new(self):
        """
        Quickly creates a new initialized NumberFieldElement with the same
        parent as self.

        EXAMPLES:

        This is called implicitly in multiplication::

            sage: O = EquationOrder([x^2 + 18, x^3 + 2], 'a,b')
            sage: c = O.1 * O.2; c
            (-23321*b^2 - 9504*b + 10830)*a + 10152*b^2 - 104562*b - 110158
            sage: parent(c) == O
            True
        """
        cdef OrderElement_relative x
        x = <OrderElement_relative>PY_NEW_SAME_TYPE(self)
        x._parent = self._parent
        x._number_field = self._parent.number_field()
        x.__fld_numerator = self.__fld_numerator
        x.__fld_denominator = self.__fld_denominator
        return x

    cpdef RingElement _div_(self, RingElement other):
        r"""
        Implement division, checking that the result has the right parent.
        It's not so crucial what the parent actually is, but it is crucial
        that the returned value really is an element of its supposed
        parent. This fixes trac #4190.

        EXAMPLES::

            sage: K1.<a> = NumberField(x^3 - 17)
            sage: R.<y> = K1[]
            sage: K2 = K1.extension(y^2 - a, 'b')
            sage: OK2 = K2.order(K2.gen()) # (not maximal)
            sage: b = OK2.gens()[1]; b
            b
            sage: (17/b).parent() is K2 # indirect doctest
            True
            sage: (17/b) in OK2 # indirect doctest
            True
            sage: (17/b^7) in OK2 # indirect doctest
            False
        """
        cdef NumberFieldElement_relative x
        x = NumberFieldElement_relative._div_(self, other)
        return self._parent.number_field()(x)

    def __invert__(self):
        r"""
        Implement division, checking that the result has the right parent.
        See trac #4190.

        EXAMPLES::

            sage: K1.<a> = NumberField(x^3 - 17)
            sage: R.<y> = K1[]
            sage: K2 = K1.extension(y^2 - a, 'b')
            sage: OK2 = K2.order(K2.gen()) # (not maximal)
            sage: b = OK2.gens()[1]; b
            b
            sage: b.parent() is OK2
            True
            sage: (~b).parent() is K2
            True
            sage: (~b) in OK2 # indirect doctest
            False
            sage: b**(-1) in OK2 # indirect doctest
            False
        """
        return self._parent.number_field()(NumberFieldElement_relative.__invert__(self))

    def inverse_mod(self, I):
        r"""
        Return an inverse of self modulo the given ideal.

        INPUT:


        -  ``I`` - may be an ideal of self.parent(), or an
           element or list of elements of self.parent() generating a nonzero
           ideal. A ValueError is raised if I is non-integral or is zero.
           A ZeroDivisionError is raised if I + (x) != (1).


        EXAMPLES::

            sage: E.<a,b> = NumberField([x^2 - x + 2, x^2 + 1])
            sage: OE = E.ring_of_integers()
            sage: t = OE(b - a).inverse_mod(17*b)
            sage: t*(b - a) - 1 in E.ideal(17*b)
            True
            sage: t.parent() == OE
            True
        """
        R = self.parent()
        return R(_inverse_mod_generic(self, I))

    def charpoly(self, var='x'):
        r"""
        The characteristic polynomial of this order element over its base ring.

        This special implementation works around bug \#4738.  At this
        time the base ring of relative order elements is ZZ; it should
        be the ring of integers of the base field.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: K.<a,b> = NumberField([x^2 + 1, x^2 - 3])
            sage: OK = K.maximal_order(); OK.basis()
            [1, 1/2*a - 1/2*b, -1/2*b*a + 1/2, a]
            sage: charpoly(OK.1)
            x^2 + b*x + 1
            sage: charpoly(OK.1).parent()
            Univariate Polynomial Ring in x over Maximal Order in Number Field in b with defining polynomial x^2 - 3
            sage: [ charpoly(t) for t in OK.basis() ]
            [x^2 - 2*x + 1, x^2 + b*x + 1, x^2 - x + 1, x^2 + 1]
        """
        R = self.parent().number_field().base_field().ring_of_integers()[var]
        return R(self.matrix().charpoly(var))

    def minpoly(self, var='x'):
        r"""
        The minimal polynomial of this order element over its base ring.

        This special implementation works around bug \#4738.  At this
        time the base ring of relative order elements is ZZ; it should
        be the ring of integers of the base field.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: K.<a,b> = NumberField([x^2 + 1, x^2 - 3])
            sage: OK = K.maximal_order(); OK.basis()
            [1, 1/2*a - 1/2*b, -1/2*b*a + 1/2, a]
            sage: minpoly(OK.1)
            x^2 + b*x + 1
            sage: charpoly(OK.1).parent()
            Univariate Polynomial Ring in x over Maximal Order in Number Field in b with defining polynomial x^2 - 3
            sage: _, u, _, v = OK.basis()
            sage: t = 2*u - v; t
            -b
            sage: t.charpoly()
            x^2 + 2*b*x + 3
            sage: t.minpoly()
            x + b

            sage: t.absolute_charpoly()
            x^4 - 6*x^2 + 9
            sage: t.absolute_minpoly()
            x^2 - 3
        """
        K = self.parent().number_field()
        R = K.base_field().ring_of_integers()[var]
        return R(K(self).minpoly(var))

    def absolute_charpoly(self, var='x'):
        r"""
        The absolute characteristic polynomial of this order element over ZZ.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: K.<a,b> = NumberField([x^2 + 1, x^2 - 3])
            sage: OK = K.maximal_order()
            sage: _, u, _, v = OK.basis()
            sage: t = 2*u - v; t
            -b
            sage: t.absolute_charpoly()
            x^4 - 6*x^2 + 9
            sage: t.absolute_minpoly()
            x^2 - 3
            sage: t.absolute_charpoly().parent()
            Univariate Polynomial Ring in x over Integer Ring
        """
        K = self.parent().number_field()
        R = ZZ[var]
        return R(K(self).absolute_charpoly(var))

    def absolute_minpoly(self, var='x'):
        r"""
        The absolute minimal polynomial of this order element over ZZ.

        EXAMPLES::

            sage: x = ZZ['x'].0
            sage: K.<a,b> = NumberField([x^2 + 1, x^2 - 3])
            sage: OK = K.maximal_order()
            sage: _, u, _, v = OK.basis()
            sage: t = 2*u - v; t
            -b
            sage: t.absolute_charpoly()
            x^4 - 6*x^2 + 9
            sage: t.absolute_minpoly()
            x^2 - 3
            sage: t.absolute_minpoly().parent()
            Univariate Polynomial Ring in x over Integer Ring
        """
        K = self.parent().number_field()
        R = ZZ[var]
        return R(K(self).absolute_minpoly(var))



class CoordinateFunction:
    r"""
    This class provides a callable object which expresses
    elements in terms of powers of a fixed field generator `\alpha`.

    EXAMPLE::

        sage: K.<a> = NumberField(x^2 + x + 3)
        sage: f = (a + 1).coordinates_in_terms_of_powers(); f
        Coordinate function that writes elements in terms of the powers of a + 1
        sage: f.__class__
        <class sage.rings.number_field.number_field_element.CoordinateFunction at ...>
        sage: f(a)
        [-1, 1]
        sage: f == loads(dumps(f))
        True
    """
    def __init__(self, NumberFieldElement alpha, W, to_V):
        r"""
        EXAMPLE::

            sage: K.<a> = NumberField(x^2 + x + 3)
            sage: f = (a + 1).coordinates_in_terms_of_powers(); f # indirect doctest
            Coordinate function that writes elements in terms of the powers of a + 1
        """
        self.__alpha = alpha
        self.__W = W
        self.__to_V = to_V
        self.__K = alpha.number_field()

    def __repr__(self):
        r"""
        EXAMPLE::

            sage: K.<a> = NumberField(x^2 + x + 3)
            sage: f = (a + 1).coordinates_in_terms_of_powers(); repr(f) # indirect doctest
            'Coordinate function that writes elements in terms of the powers of a + 1'
        """
        return "Coordinate function that writes elements in terms of the powers of %s"%self.__alpha

    def alpha(self):
        r"""
        EXAMPLE::

            sage: K.<a> = NumberField(x^3 + 2)
            sage: (a + 2).coordinates_in_terms_of_powers().alpha()
            a + 2
        """
        return self.__alpha

    def __call__(self, x):
        r"""
        EXAMPLE::

            sage: K.<a> = NumberField(x^3 + 2)
            sage: f = (a + 2).coordinates_in_terms_of_powers()
            sage: f(1/a)
            [-2, 2, -1/2]
            sage: f(ZZ(2))
            [2, 0, 0]
            sage: L.<b> = K.extension(x^2 + 7)
            sage: g = (a + b).coordinates_in_terms_of_powers()
            sage: g(a/b)
            [-3379/5461, -371/10922, -4125/38227, -15/5461, -14/5461, -9/76454]
            sage: g(a)
            [4459/10922, -4838/5461, -273/5461, -980/5461, -9/10922, -42/5461]
            sage: f(b)
            Traceback (most recent call last):
            ...
            TypeError: Cannot coerce element into this number field
        """
        from sage.all import parent
        if not self.__K.has_coerce_map_from(parent(x)):
            raise TypeError, "Cannot coerce element into this number field"
        return self.__W.coordinates(self.__to_V(self.__K(x)))

    def __cmp__(self, other):
        r"""
        EXAMPLE::

            sage: K.<a> = NumberField(x^4 + 1)
            sage: f = (a + 1).coordinates_in_terms_of_powers()
            sage: f == loads(dumps(f))
            True
            sage: f == (a + 2).coordinates_in_terms_of_powers()
            False
            sage: f == NumberField(x^2 + 3,'b').gen().coordinates_in_terms_of_powers()
            False
        """
        return cmp(self.__class__, other.__class__) or cmp(self.__K, other.__K) or cmp(self.__alpha, other.__alpha)



#################

cdef void _ntl_poly(f, ZZX_c *num, ZZ_c *den):
    cdef long i
    cdef ZZ_c coeff
    cdef ntl_ZZX _num
    cdef ntl_ZZ _den

    __den = f.denominator()
    (<Integer>ZZ(__den))._to_ZZ(den)

    __num = f * __den
    for i from 0 <= i <= __num.degree():
        (<Integer>ZZ(__num[i]))._to_ZZ(&coeff)
        ZZX_SetCoeff( num[0], i, coeff )


