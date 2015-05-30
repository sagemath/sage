"""
Finite field elements implemented via PARI's POLMOD type (deprecated)

EXAMPLES::

    sage: K = FiniteField(2)
    sage: V = VectorSpace(K,3)
    sage: w = V([0,1,2])
    sage: K(1)*w
    (0, 1, 0)

We do some arithmetic involving a bigger field and a Conway
polynomial, i.e., we verify compatibility condition.

::

    sage: f = conway_polynomial(2,63)
    sage: K.<a> = GF(2**63, name='a', modulus=f, impl='pari_mod')
    doctest:...: DeprecationWarning: The "pari_mod" finite field implementation is deprecated
    See http://trac.sagemath.org/17297 for details.
    sage: n = f.degree()
    sage: m = 3;
    sage: e = (2^n - 1) / (2^m - 1)
    sage: c = a^e
    sage: conway = conway_polynomial(2,m)
    sage: conway(c) == 0
    True
"""


import operator

import sage.structure.element as element
import sage.rings.integer_ring as integer_ring
from sage.rings.integer import Integer
from sage.libs.pari.all import pari, pari_gen
from sage.rings.finite_rings.element_base import FinitePolyExtElement
import sage.rings.field_element as field_element
import sage.rings.finite_rings.integer_mod as integer_mod
from element_base import is_FiniteFieldElement
from sage.modules.free_module_element import FreeModuleElement
from sage.structure.dynamic_class import dynamic_class
from sage.categories.finite_fields import FiniteFields

class FiniteField_ext_pariElement(FinitePolyExtElement):
    """
    An element of a finite field.

    Create elements by first defining the finite field ``F``, then use the
    notation ``F(n)``, for ``n`` an integer, or let ``a = F.gen()`` and
    write the element in terms of ``a``.

    EXAMPLES::

        sage: K = FiniteField(10007^10, 'a', impl='pari_mod')
        sage: a = K.gen(); a
        a
        sage: loads(a.dumps()) == a
        True
        sage: K = GF(10007)
        sage: a = K(938); a
        938
        sage: loads(a.dumps()) == a
        True

    TESTS::

        sage: K.<a> = GF(2^16, impl='pari_mod')
        sage: K(0).is_zero()
        True
        sage: (a - a).is_zero()
        True
        sage: a - a
        0
        sage: a == a
        True
        sage: a - a == 0
        True
        sage: a - a == K(0)
        True
        sage: TestSuite(a).run()
    """
    def __init__(self, parent, value, value_from_pari=False):
        """
        Create element of a finite field.

        INPUT:

        - ``parent`` -- A finite field, the parent of this element.

        - ``value`` -- Anything that can be converted into the PARI
          interface and can be interpreted as an element of ``parent``.

        - ``value_from_pari`` -- optional bool, default ``False``. If it
          evaluates as ``True``, then ``value`` *must* be an element of the
          PARI interface, and there will be no conversion.

        .. NOTE::

            If the given value is a list or an element of the vector space
            associated with the given parent, then it is interpreted as
            the list of coefficients of a polynomial over the prime subfield,
            and that polynomial is interpreted as an element of the given
            parent. The empty list results in zero.

            If ``value_from_pari`` is ``True`` then it is assumed that the
            given value is a suitable representation of the element in PARI,
            and there is no conversion. Hence, it is very fast, but must be
            used with care.

        EXAMPLES::

            sage: k = FiniteField(9, 'a', impl='pari_mod')
            sage: a = k(11); a
            2
            sage: a.parent()
            Finite Field in a of size 3^2
            sage: V = k.vector_space(); v = V((1,2))
            sage: k(v)
            2*a + 1

        We create elements using a list and verify that :trac:`10486` has
        been fixed::

            sage: k = FiniteField(3^11, 't', impl='pari_mod')
            sage: x = k([1,0,2,1]); x
            t^3 + 2*t^2 + 1
            sage: x + x + x
            0
            sage: pari(x)
            Mod(Mod(1, 3)*a^3 + Mod(2, 3)*a^2 + Mod(1, 3), Mod(1, 3)*a^11 + Mod(2, 3)*a^2 + Mod(1, 3))

        If the list is longer than the degree, we just get the result
        modulo the modulus::

            sage: R.<a> = PolynomialRing(GF(5))
            sage: k = FiniteField(5^2, 't', impl='pari_mod', modulus=a^2-2)
            sage: x = k([0,0,0,1]); x
            2*t
            sage: pari(x)
            Mod(Mod(2, 5)*a, Mod(1, 5)*a^2 + Mod(3, 5))

        When initializing from a list, the elements are first coerced
        to the prime field (:trac:`11685`)::

            sage: k = FiniteField(3^11, 't', impl='pari_mod')
            sage: k([ 0, 1/2 ])
            2*t
            sage: k([ k(0), k(1) ])
            t
            sage: k([ GF(3)(2), GF(3^5,'u')(1) ])
            t + 2
            sage: R.<x> = PolynomialRing(k)
            sage: k([ R(-1), x/x ])
            t + 2

        We demonstrate the creation of an element via polynomials::

            sage: k.polynomial()
            t^11 + 2*t^2 + 1
            sage: P = k.polynomial_ring()
            sage: k(P.0^11)
            t^2 + 2

        We demonstrate the creation of an element via a vector::

            sage: V = k.vector_space()
            sage: V
            Vector space of dimension 11 over Finite Field of size 3
            sage: v = V([0,1,2,0,1,2,0,1,2,0,1])
            sage: k(v)
            t^10 + 2*t^8 + t^7 + 2*t^5 + t^4 + 2*t^2 + t

        TESTS:

        Check that zeros are created correctly (:trac:`11685`)::

            sage: K = FiniteField(3^11, 't', impl='pari_mod'); a = K.0
            sage: v = 0; pari(K(v)).lift()
            Mod(0, 3)
            sage: v = Mod(0,3); pari(K(v)).lift()
            Mod(0, 3)
            sage: v = pari(0); pari(K(v)).lift()
            Mod(0, 3)
            sage: v = pari("Mod(0,3)"); pari(K(v)).lift()
            Mod(0, 3)
            sage: v = []; pari(K(v)).lift()
            Mod(0, 3)
            sage: v = [0]; pari(K(v)).lift()
            Mod(0, 3)
            sage: v = [0,0]; pari(K(v)).lift()
            Mod(0, 3)
            sage: v = pari("Pol(0)"); pari(K(v)).lift()
            Mod(0, 3)
            sage: v = pari("Mod(0, %s)"%K._pari_modulus()); pari(K(v)).lift()
            Mod(0, 3)
            sage: v = pari("Mod(Pol(0), %s)"%K._pari_modulus()); pari(K(v)).lift()
            Mod(0, 3)
            sage: v = K(1) - K(1); pari(K(v)).lift()
            Mod(0, 3)
            sage: v = K([1]) - K([1]); pari(K(v)).lift()
            Mod(0, 3)
            sage: v = a - a; pari(K(v)).lift()
            Mod(0, 3)
            sage: v = K(1)*0; pari(K(v)).lift()
            Mod(0, 3)
            sage: v = K([1])*K([0]); pari(K(v)).lift()
            Mod(0, 3)
            sage: v = a*0; pari(K(v)).lift()
            Mod(0, 3)

        The following test documents the optional argument ``value_from_pari``.
        It is for internal use only and greatly improves the speed in
        arithmetic operations. However, the example shows why it must only be
        used carefully::

            sage: from sage.rings.finite_rings.element_ext_pari import FiniteField_ext_pariElement
            sage: a = FiniteField_ext_pariElement(K,pari(0),value_from_pari=True)
            sage: a._pari_().type()
            't_INT'
            sage: K(0)._pari_().type()
            't_POLMOD'
        """
        field_element.FieldElement.__init__(self, parent)
        self.__class__ = dynamic_FiniteField_ext_pariElement

        # If value_from_pari is True, directly set self.__value to value.
        # This assumes that value is a POLMOD with the correct modulus
        # whose lift is either an INTMOD or a POL with INTMOD
        # coefficients. In practice, this means that value comes from a
        # PARI calculation with other elements of this finite field.
        # This assumption is not checked.
        if value_from_pari:
            self.__value = value
            return

        try:
            if isinstance(value, pari_gen):
                try:
                    # In some cases we get two different versions of
                    # the 0 element of a finite field.  This has to do
                    # with how PARI's Mod() function works -- it treats
                    # polynomials different from integers.
                    # Also, we need to fix things like ``1/Mod(3,5)`` when
                    # we want ``Mod(2,5)``.
                    # These issues are solved by first simplifying the
                    # given value.
                    value = value.simplify()
                    if value.type() != "t_POLMOD":
                        value = value.Mod(parent._pari_modulus())
                    self.__value = value * parent._pari_one()
                except RuntimeError:
                    raise TypeError("no possible coercion implemented")
            elif isinstance(value, FiniteField_ext_pariElement):
                if parent != value.parent():
                    raise TypeError("no coercion implemented")
                else:
                    self.__value = value.__value
            elif isinstance(value, FreeModuleElement):
                if parent.vector_space() != value.parent():
                    raise TypeError("e.parent must match self.vector_space")
                self.__value = pari(0).Mod(parent._pari_modulus())*parent._pari_one()
                for i in range(len(value)):
                    self.__value = self.__value + pari(int(value[i])).Mod(parent._pari_modulus())*pari("a^%s"%i)
            elif isinstance(value, list):
                # AUTHOR: Jeroen Demeyer
                # See Trac #10486 and #11685 for comments about this
                # implementation

                # We need a special case for the empty list because
                # PARI's ``Pol([])`` returns an exact 0 while we want
                # ``Mod(0,3)``.
                if not value:
                    value = [0]
                try:
                    # Convert the list to PARI, then mod out the
                    # characteristic (PARI can do this directly for lists).
                    # If we get only INTMODs, we can do the conversion
                    # directly.
                    parilist = pari(value).Mod(parent.characteristic())
                    if not all(c.type() == "t_INTMOD" for c in parilist):
                        raise RuntimeError
                except RuntimeError:
                    # That didn't work, do it in a more general but also
                    # slower way: first convert all list elements to the
                    # prime field.
                    GFp = parent.prime_subfield()
                    parilist = pari([GFp(c) for c in value])
                self.__value = parilist.Polrev("a").Mod(parent._pari_modulus())
            elif isinstance(value, str):
                raise TypeError("value must not be a string")
            else:
                try:
                    self.__value = pari(value).Mod(parent._pari_modulus())*parent._pari_one()
                except RuntimeError:
                    raise TypeError("no coercion implemented")

        except (AttributeError, TypeError):
            raise TypeError("unable to coerce")

    def __hash__(self):
        """
        The hash of this element is the hash of the underlying polynomial.

        EXAMPLES::

            sage: k.<a> = GF(3^15, impl='pari_mod')
            sage: R = GF(3)['a']; aa = R.gen()
            sage: hash(a^2 + 1) == hash(aa^2 + 1)
            True
        """
        return hash(self.polynomial())

    def polynomial(self):
        """
        Elements of a finite field are represented as a polynomial modulo a
        modulus. This function returns the representing polynomial as an
        element of the polynomial ring over the prime finite field, with
        the same variable as the finite field.

        EXAMPLES:

            sage: k = FiniteField(3**2, 't', impl='pari_mod')
            sage: k.gen().polynomial()
            t

        The variable can be any string::

            sage: k = FiniteField(3**4, "alpha", impl='pari_mod')
            sage: a = k.gen()
            sage: a.polynomial()
            alpha
            sage: (a**2 + 1).polynomial()
            alpha^2 + 1
            sage: (a**2 + 1).polynomial().parent()
            Univariate Polynomial Ring in alpha over Finite Field of size 3
        """
        return self.parent().polynomial_ring()(self.__value.lift())

    def __copy__(self):
        """
        Return a copy of this element.

        EXAMPLES::

            sage: k = FiniteField(3**3, 'a', impl='pari_mod')
            sage: a = k(5)
            sage: a
            2
            sage: copy(a)
            2
            sage: b = copy(a)
            sage: a == b
            True
            sage: a is b
            False
            sage: a is a
            True
        """
        return FiniteField_ext_pariElement(self.parent(), self.__value, value_from_pari=True)

    def _pari_(self, var=None):
        """
        Return PARI object corresponding to this finite field element.

        EXAMPLES::

            sage: k.<a> = FiniteField(3**3, 'a', impl='pari_mod')
            sage: b = a**2 + 2*a + 1
            sage: b._pari_()
            Mod(Mod(1, 3)*a^2 + Mod(2, 3)*a + Mod(1, 3), Mod(1, 3)*a^3 + Mod(2, 3)*a + Mod(1, 3))

        Looking at the PARI representation of a finite field element, it's
        no wonder people find PARI difficult to work with directly. Compare
        our representation::

            sage: b
            a^2 + 2*a + 1
            sage: b.parent()
            Finite Field in a of size 3^3
        """
        if var is None:
            var = self.parent().variable_name()
        if var == 'a':
            return self.__value
        else:
            return self.__value.subst('a', var)

    def _pari_init_(self):
        """
        The string producing this finite field element in PARI.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_mod')
            sage: a._pari_init_()
            'Mod(Mod(1, 3)*a, Mod(1, 3)*a^17 + Mod(2, 3)*a + Mod(1, 3))'
        """
        return str(self.__value)

    def _magma_init_(self, magma):
        """
        Return a string representation of self that Magma can understand.

        EXAMPLES::

            sage: GF(7)(3)._magma_init_(magma)            # optional - magma
            'GF(7)!3'
        """
        km = magma(self.parent())
        vn = km.gen(1).name()
        return ("%s"%(self.__value.lift().lift())).replace('a',vn)

    def _gap_init_(self):
        r"""
        Supports returning corresponding GAP object. This can be slow since
        non-prime GAP finite field elements are represented as powers of a
        generator for the multiplicative group, so the discrete log problem
        must be solved.

        .. NOTE::

           The order of the parent field must be `\leq 65536`.

        EXAMPLES::

            sage: F.<a> = FiniteField(8, 'a', impl='pari_mod', modulus='primitive')
            sage: gap(a) # indirect doctest
            Z(2^3)
            sage: b = F.multiplicative_generator()
            sage: a = b^3
            sage: gap(a)
            Z(2^3)^3
            sage: gap(a^3)
            Z(2^3)^2

        You can specify the instance of the Gap interpreter that is used::

            sage: F.<a> = FiniteField(next_prime(200)^2, 'a', impl='pari_mod', modulus='primitive')
            sage: a._gap_ (gap)
            Z(211^2)
            sage: (a^20)._gap_(gap)
            Z(211^2)^20

        Gap only supports relatively small finite fields::

            sage: F.<a> = FiniteField(next_prime(1000)^2, 'a', impl='pari_mod', modulus='primitive')
            sage: gap._coerce_(a)
            Traceback (most recent call last):
            ...
            TypeError: order must be at most 65536
        """
        F = self.parent()
        if F.order() > 65536:
            raise TypeError("order must be at most 65536")

        if self == 0:
            return '0*Z(%s)'%F.order()
        assert F.degree() > 1
        g = F.multiplicative_generator()
        n = self.log(g)
        return 'Z(%s)^%s'%(F.order(), n)

    def _repr_(self):
        """
        String representation of this element.

        EXAMPLES::

            sage: k.<c> = GF(3^17, impl='pari_mod')
            sage: c^20 # indirect doctest
            c^4 + 2*c^3
        """
        return ("%s"%(self.__value.lift().lift())).replace('a',self.parent().variable_name())

    def _add_(self, right):
        """
        Addition.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_mod')
            sage: a + a^2 # indirect doctest
            a^2 + a
        """
        return FiniteField_ext_pariElement(self.parent(), self.__value + right.__value, value_from_pari=True)

    def _sub_(self, right):
        """
        Subtraction.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_mod')
            sage: a - a # indirect doctest
            0
        """
        return FiniteField_ext_pariElement(self.parent(), self.__value - right.__value, value_from_pari=True)

    def _mul_(self, right):
        """
        Multiplication.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_mod')
            sage: (a^12 + 1)*(a^15 - 1) # indirect doctest
            a^15 + 2*a^12 + a^11 + 2*a^10 + 2
        """
        return FiniteField_ext_pariElement(self.parent(), self.__value * right.__value, value_from_pari=True)

    def _div_(self, right):
        """
        Division.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_mod')
            sage: (a - 1) / (a + 1) # indirect doctest
            2*a^16 + a^15 + 2*a^14 + a^13 + 2*a^12 + a^11 + 2*a^10 + a^9 + 2*a^8 + a^7 + 2*a^6 + a^5 + 2*a^4 + a^3 + 2*a^2 + a + 1
        """
        if right.__value == 0:
            raise ZeroDivisionError
        return FiniteField_ext_pariElement(self.parent(), self.__value / right.__value, value_from_pari=True)

    def __int__(self):
        """
        Lifting to a python int, if possible.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_mod'); b = k(2)
            sage: int(b)
            2
            sage: int(a)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce PARI a to an Integer
        """
        try:
            return int(self.__value.lift().lift())
        except ValueError:
            raise TypeError("cannot coerce to int")

    def _integer_(self, ZZ=None):
        """
        Lifting to a sage integer if possible.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_mod'); b = k(2)
            sage: b._integer_()
            2
            sage: a._integer_()
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce PARI a to an Integer
        """
        return self.lift()

    def __long__(self):
        """
        Lifting to a python long, if possible.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_mod'); b = k(2)
            sage: long(b)
            2L
        """
        try:
            return long(self.__value.lift().lift())
        except ValueError:
            raise TypeError("cannot coerce to long")

    def __float__(self):
        """
        Lifting to a python float, if possible.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_mod'); b = k(2)
            sage: float(b)
            2.0
        """
        try:
            return float(self.__value.lift().lift())
        except ValueError:
            raise TypeError("cannot coerce to float")

    def __pow__(self, _right):
        """
        TESTS::

            sage: K.<a> = GF(5^10, impl='pari_mod')
            sage: n = (2*a)/a

        Naively compute `n^{-15}` in PARI, note that the result is `1/3`.
        This is mathematically correct (modulo 5), but not what we want.
        In particular, comparisons will fail::

            sage: pari(n)^-15
            Mod(1/Mod(3, 5), Mod(1, 5)*a^10 + Mod(3, 5)*a^5 + Mod(3, 5)*a^4 + Mod(2, 5)*a^3 + Mod(4, 5)*a^2 + Mod(1, 5)*a + Mod(2, 5))

        We need to :meth:`simplify()` the result (which is done in the
        :class:`FiniteField_ext_pariElement` constructor::

            sage: n^-15
            2

        Large exponents are not a problem::

            sage: e = 3^10000
            sage: a^e
            2*a^9 + a^5 + 4*a^4 + 4*a^3 + a^2 + 3*a
            sage: a^(e % (5^10 - 1))
            2*a^9 + a^5 + 4*a^4 + 4*a^3 + a^2 + 3*a
        """
        right = int(_right)
        if right != _right:
             raise ValueError
        # It is important to set value_from_pari=False, see doctest above!
        return FiniteField_ext_pariElement(self.parent(), self.__value**right, value_from_pari=False)

    def __neg__(self):
        """
        Negation.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_mod')
            sage: -a
            2*a
        """
        return FiniteField_ext_pariElement(self.parent(), -self.__value, value_from_pari=True)

    def __pos__(self):
        """
        Unitary positive operator...

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_mod')
            sage: +a
            a
        """
        return self

    def __abs__(self):
        """
        Absolute value, which is not defined.

        EXAMPLES::

            sage: k.<a> = GF(3^17, impl='pari_mod')
            sage: abs(a)
            Traceback (most recent call last):
            ...
            ArithmeticError: absolute value not defined
        """
        raise ArithmeticError("absolute value not defined")

    def __invert__(self):
        """
        EXAMPLES::

            sage: k.<a> = FiniteField(9, impl='pari_mod')
            sage: ~a
            a + 2
            sage: (a+1)*a
            2*a + 1
            sage: ~((2*a)/a)
            2
        """

        if self.__value == 0:
            raise ZeroDivisionError("Cannot invert 0")
        return FiniteField_ext_pariElement(self.parent(), ~self.__value, value_from_pari=True)

    def lift(self):
        """
        If this element lies in a prime finite field, return a lift of this
        element to an integer.

        EXAMPLES::

            sage: k.<t> = GF(next_prime(10^10)^2, impl='pari_mod')
            sage: a = k(17)/k(19)
            sage: b = a.lift(); b
            7894736858
            sage: b.parent()
            Integer Ring
        """
        return integer_ring.IntegerRing()(self.__value.lift().lift())


    def __cmp__(self, other):
        """
        Compare an element of a finite field with other. If other is not an
        element of a finite field, an attempt is made to coerce it so it is
        one.

        EXAMPLES::

            sage: k.<a> = FiniteField(3**3, impl='pari_mod')
            sage: a == 1
            False
            sage: a**0 == 1
            True
            sage: a == a
            True
            sage: a < a**2
            True
            sage: a > a**2
            False
        """
        return cmp(self.__value.lift().lift().Pol(), other.__value.lift().lift().Pol())

    def log(self, base):
        """
        Return `x` such that `b^x = a`, where `x`
        is `a` and `b` is the base.

        INPUT:

        -  ``b`` -- finite field element that generates the
           multiplicative group.

        OUTPUT:

        Integer `x` such that `a^x = b`, if it exists. Raises a
        ``ValueError`` exception if no such `x` exists.

        EXAMPLES::

            sage: F = GF(17)
            sage: F(3^11).log(F(3))
            11
            sage: F = GF(113)
            sage: F(3^19).log(F(3))
            19
            sage: F = GF(next_prime(10000))
            sage: F(23^997).log(F(23))
            997

        ::

            sage: F = FiniteField(2^10, 'a', impl='pari_mod')
            sage: g = F.gen()
            sage: b = g; a = g^37
            sage: a.log(b)
            37
            sage: b^37; a
            a^8 + a^7 + a^4 + a + 1
            a^8 + a^7 + a^4 + a + 1

        AUTHORS:

        - David Joyner and William Stein (2005-11)
        """
        from  sage.groups.generic import discrete_log

        b = self.parent()(base)
        # TODO: This function is TERRIBLE!
        return discrete_log(self, b)

dynamic_FiniteField_ext_pariElement = None
def _late_import():
    """
    Used to reset the class of PARI finite field elements in their initialization.

    EXAMPLES::

        sage: from sage.rings.finite_rings.element_ext_pari import FiniteField_ext_pariElement
        sage: k.<a> = GF(3^17, impl='pari_mod')
        sage: a.__class__ is FiniteField_ext_pariElement # indirect doctest
        False
    """
    global dynamic_FiniteField_ext_pariElement
    dynamic_FiniteField_ext_pariElement = dynamic_class("%s_with_category"%FiniteField_ext_pariElement.__name__, (FiniteField_ext_pariElement, FiniteFields().element_class), doccls=FiniteField_ext_pariElement)

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.rings.finite_field_element', 'FiniteField_ext_pariElement', FiniteField_ext_pariElement)
