r"""
Elements lying in extension of rings

AUTHOR:

- Xavier Caruso (2019)
"""

#############################################################################
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#                  2022 Julian Rüth <julian.rueth@fsfe.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************


from sage.ext.stdsage cimport PY_NEW
from sage.misc.cachefunc import cached_method
from sage.cpython.getattr cimport AttributeErrorMessage
from sage.cpython.getattr import dir_with_other_class
from sage.misc.latex import latex

from sage.structure.category_object import normalize_names
from sage.structure.element cimport CommutativeAlgebraElement
from sage.structure.element cimport Element, parent
from sage.rings.integer_ring import ZZ
from sage.categories.fields import Fields
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.rings.ring_extension_conversion cimport backend_parent
from sage.rings.ring_extension_conversion cimport to_backend, from_backend


cdef class RingExtensionElement(CommutativeAlgebraElement):
    r"""
    Generic class for elements lying in ring extensions.

    TESTS::

        sage: K = GF(5^4).over()
        sage: x = K.random_element()
        sage: TestSuite(x).run()

    """
    def __init__(self, extension, x):
        r"""
        Initialize this element.

        INPUT:

        - ``extension`` -- the parent of this element

        - ``x`` -- an element of the backend

        TESTS::

            sage: Q = QQ.over(ZZ)
            sage: x = Q(1/2)
            sage: x
            1/2
        """
        from sage.rings.ring_extension import RingExtension_generic
        if not isinstance(extension, RingExtension_generic):
            raise TypeError(f"{extension} is not a ring extension")

        if x.parent() is not extension._backend:
            # We are strict here since the base could be the same as the backend so we can't rely on coercion
            raise ValueError("x must be an element of the backend")
        CommutativeAlgebraElement.__init__(self, extension)
        self._backend = x

    def __reduce__(self):
        """
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: K = GF(5^3).over()
            sage: x = K.random_element()
            sage: type(x)
            <class 'sage.rings.ring_extension_element.RingExtensionWithBasisElement'>
            sage: loads(dumps(x)) == x
            True
        """
        return self._parent._from_backend_morphism, (self._backend,)

    def __getattr__(self, name):
        """
        If the parent of this element was created with ``import_methods = True``,
        return a wrapper to the corresponding method of the backend element
        (if it exists).

        EXAMPLES::

            sage: A.<a> = QQ.extension(x^2 - 2)
            sage: K.<a> = A.over()  # over QQ

            sage: hasattr(a, 'continued_fraction')
            True
            sage: a.continued_fraction()
            [1; (2)*]
        """
        try:
            return self.getattr_from_category(name)
        except AttributeError:
            pass
        method = None
        if self._parent._import_methods and hasattr(self._backend, name):
            method = getattr(self._backend, name)
        if not callable(method):
            raise AttributeError(AttributeErrorMessage(self, name))

        def wrapper(*args, **kwargs):
            output = method(*to_backend(args), **to_backend(kwargs))
            return from_backend(output, self._parent)
        wrapper.__doc__ = method.__doc__
        return wrapper

    def __dir__(self):
        """
        Return the list of all the attributes of this element;
        if the parent of this element was created with ``import_methods = True``,
        concatenate this list with the list of all the methods of the backend
        element.

        EXAMPLES::

            sage: A.<a> = QQ.extension(x^2 - 2)
            sage: K.<a> = A.over()

            sage: dir(a)
            ['__abs__',
             '__add__',
             ...
             'complex_embeddings',
             'conjugate',
             'continued_fraction',
             'continued_fraction_list',
             ...
             'trace',
             'valuation',
             'vector',
             'xgcd']
        """
        d = dir_with_other_class(self, self._parent.category().element_class)
        if not self._parent._import_methods:
            return d
        for name in dir(self._backend):
            try:
                attribute = getattr(self._backend, name)
                if callable(attribute):
                    d.append(name)
            except AttributeError:
                pass
        return sorted(set(d))

    def __hash__(self):
        """
        Return a hash of this element.

        EXAMPLES:

            sage: E.<a> = GF(5^3).over()
            sage: hash(a)
            5
        """
        return hash(self._backend)

    def _repr_(self, **options):
        r"""
        Return a string representation of this element.

        Do not override this method in subclasses;
        instead override the method :meth:`_repr_extension`.

        TESTS::

            sage: K.<a> = GF(5^2).over()
            sage: L.<b> = GF(5^4).over(K)
            sage: b._repr_()
            'b'
        """
        parent = self._parent
        if 'print_elements_as' in options:
            print_as = options.pop('print_elements_as')
        else:
            print_as = parent._print_options.get('print_elements_as')
        if print_as is not None:
            return print_as(self._backend)._repr_(**options)
        print_options = parent._print_options.copy()
        for (name, value) in options.items():
            method = None
            if hasattr(parent, '_print_option_' + name):
                method = getattr(parent, '_print_option_' + name)
            if not callable(method):
                raise ValueError("option '%s' does not exist" % name)
            print_options[name] = method(value)
        return self._repr_extension(**print_options)

    def _repr_extension(self, **options):
        r"""
        Return a string representation of this element.

        TESTS::

            sage: K = QQ.over(ZZ)
            sage: x = K(1/2)
            sage: x._repr_extension()
            '1/2'
        """
        return str(self._backend)

    def _latex_(self, **options):
        r"""
        Return a LaTeX representation of this element.

        Do not override this method in subclasses;
        instead override the method :meth:`_latex_extension`.

        TESTS::

            sage: K.<a> = GF(5^2).over()
            sage: L.<b> = GF(5^4).over(K)
            sage: b._latex_()
            'b'
        """
        parent = self._parent
        if 'print_elements_as' in options:
            print_as = options.pop('print_elements_as')
        else:
            print_as = parent._print_options.get('print_elements_as')
        if print_as is not None:
            return print_as(self._backend)._latex_(**options)
        print_options = parent._print_options.copy()
        for (name, value) in options.items():
            method = None
            if hasattr(parent, '_print_option_' + name):
                method = getattr(parent, '_print_option_' + name)
            if not callable(method):
                raise ValueError("option '%s' does not exist" % name)
            print_options[name] = method(value)
        return self._latex_extension(**print_options)

    def _latex_extension(self, **options):
        r"""
        Return a LaTeX representation of this element.

        TESTS::

            sage: K = QQ.over(ZZ)
            sage: x = K(1/2)
            sage: x._latex_extension()
            \frac{1}{2}
        """
        return latex(self._backend)

    cpdef _richcmp_(left, right, int op):
        r"""
        Compare this element with ``right`` according to
        the rich comparison operator ``op``.

        The comparison is performed by comparing the backend
        elements.

        INPUT:

        - ``right`` -- an element in the same parent

        - ``op`` -- the comparison operator

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()
            sage: x = K.random_element()
            sage: x == x
            True
            sage: x == x + 1
            False
            sage: x == x^25
            True
        """
        return left._backend._richcmp_(right._backend, op)

    cpdef _add_(self,other):
        r"""
        Return the sum of this element and ``other``.

        TESTS::

            sage: K = GF(5^4).over(GF(5^2))
            sage: x = K.random_element()
            sage: y = K.random_element()

            sage: (x+y).parent() is K
            True
            sage: x + y == y + x
            True
        """
        cdef RingExtensionElement ans = PY_NEW(type(self))
        ans._parent = self._parent
        ans._backend = self._backend + (<RingExtensionElement>other)._backend
        return ans

    cpdef _neg_(self):
        r"""
        Return the opposite of this element.

        TESTS::

            sage: K = GF(5^4).over(GF(5^2))
            sage: x = K.random_element()

            sage: y = -x
            sage: y.parent() is K
            True
            sage: x + y == 0
            True
        """
        cdef RingExtensionElement ans = PY_NEW(type(self))
        ans._parent = self._parent
        ans._backend = -self._backend
        return ans

    cpdef _sub_(self,other):
        r"""
        Return the difference of this element and ``other``.

        TESTS::

            sage: K = GF(5^4).over(GF(5^2))
            sage: x = K.random_element()
            sage: y = K.random_element()

            sage: (x-y).parent() is K
            True
            sage: x - y == x + (-y)
            True
        """
        cdef RingExtensionElement ans = PY_NEW(type(self))
        ans._parent = self._parent
        ans._backend = self._backend - (<RingExtensionElement>other)._backend
        return ans

    cpdef _mul_(self,other):
        r"""
        Return the product of this element and ``other``.

        TESTS::

            sage: K = GF(5^4).over(GF(5^2))
            sage: x = K.random_element()
            sage: y = K.random_element()

            sage: (x*y).parent() is K
            True
            sage: x * y == y * x
            True
        """
        cdef RingExtensionElement ans = PY_NEW(type(self))
        ans._parent = self._parent
        ans._backend = self._backend * (<RingExtensionElement>other)._backend
        return ans

    cpdef _div_(self,other):
        r"""
        Return the quotient of this element by ``other``,
        considered as an element of the fraction field.

        TESTS::

            sage: A.<a> = ZZ.extension(x^2 - 2)
            sage: OK = A.over()
            sage: a = OK(a)

            sage: b = 1/a; b
            a/2
            sage: b.parent()
            Fraction Field of Order in Number Field in a with defining polynomial x^2 - 2 over its base
            sage: a*b
            1
        """
        cdef RingExtensionElement ans
        parent = self._parent
        K = self.parent().fraction_field()
        return K.element_class(K, self._backend / (<RingExtensionElement>other)._backend)

    def additive_order(self):
        r"""
        Return the additive order of this element.

        EXAMPLES::

            sage: K.<a> = GF(5^4).over(GF(5^2))
            sage: a.additive_order()
            5
        """
        return self._backend.additive_order()

    def multiplicative_order(self):
        r"""
        Return the multiplicite order of this element.

        EXAMPLES::

            sage: K.<a> = GF(5^4).over(GF(5^2))
            sage: a.multiplicative_order()
            624
        """
        return self._backend.multiplicative_order()

    def is_unit(self):
        r"""
        Return whether if this element is a unit in this ring.

        EXAMPLES::

            sage: A.<x> = PolynomialRing(QQ)
            sage: E = A.over(QQ)
            sage: E(4).is_unit()
            True
            sage: E(x).is_unit()
            False
        """
        return self._backend.is_unit()

    def is_nilpotent(self):
        r"""
        Return whether if this element is nilpotent in this ring.

        EXAMPLES::

            sage: A.<x> = PolynomialRing(QQ)
            sage: E = A.over(QQ)
            sage: E(0).is_nilpotent()
            True
            sage: E(x).is_nilpotent()
            False
        """
        return self._backend.is_nilpotent()

    def is_prime(self):
        r"""
        Return whether this element is a prime element in this ring.

        EXAMPLES::

            sage: A.<x> = PolynomialRing(QQ)
            sage: E = A.over(QQ)
            sage: E(x^2+1).is_prime()
            True
            sage: E(x^2-1).is_prime()
            False
        """
        return self._backend.is_prime()

    def is_square(self, root=False):
        r"""
        Return whether this element is a square in this ring.

        INPUT:

        - ``root`` -- a boolean (default: ``False``); if ``True``,
          return also a square root

        EXAMPLES::

            sage: K.<a> = GF(5^3).over()
            sage: a.is_square()
            False
            sage: a.is_square(root=True)
            (False, None)

            sage: b = a + 1
            sage: b.is_square()
            True
            sage: b.is_square(root=True)
            (True, 2 + 3*a + a^2)
        """
        is_sq = self._backend.is_square()
        sq = None
        if root and is_sq:
            sq = self.sqrt(extend=False, all=False)
        if root:
            return is_sq, sq
        else:
            return is_sq

    def sqrt(self, extend=True, all=False, name=None):
        r"""
        Return a square root or all square roots of this element.

        INPUT:

        - ``extend`` -- a boolean (default: ``True``); if "True",
          return a square root in an extension ring, if necessary.
          Otherwise, raise a ``ValueError`` if the root is not in
          the ring

        - ``all`` -- a boolean (default: ``False``); if ``True``,
          return all square roots of this element, instead of just one.

        - ``name`` -- Required when ``extend=True`` and ``self`` is not a
          square. This will be the name of the generator extension.

        .. NOTE::

            The option `extend=True` is often not implemented.

        EXAMPLES::

            sage: K.<a> = GF(5^3).over()
            sage: b = a + 1
            sage: b.sqrt()
            2 + 3*a + a^2
            sage: b.sqrt(all=True)
            [2 + 3*a + a^2, 3 + 2*a - a^2]
        """
        sq = self._backend.sqrt(extend=extend, all=all)
        if all:
            gen = sq[0]
        else:
            gen = sq
        parent = self._parent
        gen_parent = gen.parent()

        parent_backend, from_parent_backend, _ = backend_parent(parent, map=True)

        if gen_parent is not parent_backend:
            from sage.rings.ring_extension import RingExtension
            if name is None:
                raise ValueError("you must specify a variable name")
            names = normalize_names(1, name)
            from sage.rings.ring_extension import RingExtensionWithGen
            constructor = (RingExtensionWithGen,
                           {'gen': gen, 'name': names[0]})
            parent = RingExtension(gen_parent, parent, (gen,), names, canonical_backend=False, constructors=[constructor])
        if all:
            return [ from_parent_backend(s) for s in sq ]
        else:
            return from_parent_backend(sq)


# Fraction fields
#################

cdef class RingExtensionFractionFieldElement(RingExtensionElement):
    r"""
    A class for elements lying in fraction fields of ring extensions.

    TESTS::

        sage: Z = ZZ.over()
        sage: Q = Z.fraction_field()
        sage: x = Q.random_element()
        sage: type(x)
        <class 'sage.rings.ring_extension_element.RingExtensionFractionFieldElement'>
        sage: TestSuite(x).run()
    """
    def __hash__(self):
        """
        Return a hash of this element.

        EXAMPLES:

            sage: E.<a> = GF(5^3).over()
            sage: hash(a)
            5
        """
        return hash(self._backend)

    def _repr_extension(self, **options):
        r"""
        Return a string representation of this element.

        TESTS::

            sage: Z = ZZ.over()
            sage: Q = Z.fraction_field()
            sage: x = Q(1/2)
            sage: x._repr_extension()
            '1/2'
            sage: R = QQ['x'].over()
            sage: K = R.fraction_field()
            sage: x = R.gen()
            sage: (x^2 + 1) / (x^2 - 1)
            (x^2 + 1)/(x^2 - 1)
            sage: x / (x + 1)
            x/(x + 1)
            sage: (x + 1)/(-x)
            (-x - 1)/x
        """
        num = self.numerator()
        denom = self.denominator()
        if denom == 1:
            sd = ""
        elif denom == -1:
            num = -num
            sd = ""
        elif denom._is_atomic():
            sd = "/%s" % denom
        elif (-denom)._is_atomic():
            sd = "/%s" % (-denom)
            num = -num
        else:
            sd = "/(%s)" % denom
        if num._is_atomic():
            return "%s%s" % (num, sd)
        else:
            return "(%s)%s" % (num, sd)

    def _latex_extension(self, **options):
        r"""
        Return a LaTeX representation of this element.

        TESTS::

            sage: Z = ZZ.over()
            sage: Q = Z.fraction_field()
            sage: x = Q(1/2)
            sage: x._latex_extension()
            '\\frac{1}{2}'
        """
        num = self.numerator()
        denom = self.denominator()
        if denom == -1:
            denom = 1
            num = -num
        from sage.rings.ring_extension import RingExtension_generic
        if isinstance(self._parent._ring, RingExtension_generic):
            snum = num._latex_(**options)
            sdenom = denom._latex_(**options)
        else:
            snum = latex(num)
            sdenom = latex(denom)
        if denom == 1:
            return snum
        else:
            return r"\frac{%s}{%s}" % (snum, sdenom)

    def numerator(self):
        r"""
        Return the numerator of this element.

        EXAMPLES::

            sage: A.<a> = ZZ.extension(x^2 - 2)
            sage: OK = A.over()  # over ZZ
            sage: K = OK.fraction_field()
            sage: K
            Fraction Field of Order in Number Field in a with defining polynomial x^2 - 2 over its base

            sage: x = K(1/a); x
            a/2
            sage: num = x.numerator(); num
            a

        The numerator is an element of the ring which was used
        to construct the fraction field::

            sage: num.parent()
            Order in Number Field in a with defining polynomial x^2 - 2 over its base
            sage: num.parent() is OK
            True

        TESTS::

            sage: x = K.random_element()
            sage: x == x.numerator() / x.denominator()
            True
        """
        R = self.parent().ring()
        if R is not self.base_ring():
            raise NotImplementedError("cannot determine numerator in this fraction field yet")

        return R._from_backend_morphism(self._backend.numerator())

    def denominator(self):
        r"""
        Return the denominator of this element.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: A.<a> = ZZ.extension(x^2 - 2)
            sage: OK = A.over()  # over ZZ
            sage: K = OK.fraction_field()
            sage: K
            Fraction Field of Order in Number Field in a with defining polynomial x^2 - 2 over its base

            sage: x = K(1/a); x
            a/2
            sage: denom = x.denominator(); denom
            2

        The denominator is an element of the ring which was used
        to construct the fraction field::

            sage: denom.parent()
            Order in Number Field in a with defining polynomial x^2 - 2 over its base
            sage: denom.parent() is OK
            True

        TESTS::

            sage: x = K.random_element()
            sage: x == x.numerator() / x.denominator()
            True
        """
        R = self.parent().ring()
        if R is not self.base_ring():
            raise NotImplementedError("cannot determine denominator in this fraction field yet")

        return R._from_backend_morphism(self._backend.denominator())


# Finite free extensions
########################

cdef class RingExtensionWithBasisElement(RingExtensionElement):
    r"""
    A class for elements lying in finite free extensions.

    TESTS::

        sage: K.<a> = GF(5^3).over()
        sage: L.<b> = GF(5^9).over(K)
        sage: type(b)
        <class 'sage.rings.ring_extension_element.RingExtensionWithBasisElement'>
        sage: TestSuite(b).run()
    """
    def __hash__(self):
        """
        Return a hash of this element.

        EXAMPLES:

            sage: E.<a> = GF(5^3).over()
            sage: hash(a)
            5
        """
        return hash(self._backend)

    def _repr_extension(self, base=None, **options):
        r"""
        Return a string representation of this element written as
        a linear combination over ``base`` in the basis provided by
        the method :meth:`basis_over`.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None``

        EXAMPLES::

            sage: K.<a> = GF(5^3).over()
            sage: L.<b> = GF(5^9).over(K)
            sage: u = 1/(a+b)

            sage: u._repr_extension(base=K)
            '(2 + 2*a) + (-1 + a - a^2)*b + (2 + 3*a + 3*a^2)*b^2'
            sage: u._repr_extension(base=GF(5))
            '2 + 2*a - b + a*b - a^2*b + 2*b^2 + 3*a*b^2 + 3*a^2*b^2'
        """
        parent = self._parent
        base = parent._check_base(base)
        coeffs = self._vector(base)
        names = parent._basis_names
        b = parent._base
        while b is not base:
            new_names = [ ]
            for y in names:
                for x in b._basis_names:
                    if x == "":
                        new_names.append(y)
                    elif y == "":
                        new_names.append(x)
                    else:
                        new_names.append(x + "*" + y)
            names = new_names
            b = b._base
        s = ""
        for i in range(len(names)):
            c = coeffs[i]
            if c.is_zero():
                continue
            sign = 1
            ss = ""
            if c == 1:
                pass # here for characteristic 2
            elif c == -1:
                sign = -1
            else:
                atomic = c._is_atomic()
                if not atomic and (-c)._is_atomic():
                    c = -c
                    sign = -sign
                    atomic = True
                sc = str(c)
                if atomic:
                    ss += sc
                else:
                    ss += "(" + sc + ")"
                if names[i] != "":
                    ss += "*"
            if ss and ss[0] == "-":
                ss = ss[1:]
                sign *= -1
            if s == "":
                if sign == -1:
                    s = "-"
            else:
                s += " + " if sign == 1 else " - "
            ss += names[i]
            if ss == "":
                ss += "1"
            s += ss
        if s == "":
            return "0"
        if s[0] == "(" and s[-1] == ")":
            s = s[1:-1]
        return s

    def _latex_extension(self, base, **options):
        r"""
        Return a LaTeX representation of this element written as
        a linear combination over ``base`` in the basis provided by
        the method :meth:`basis_over`.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None``

        EXAMPLES::

            sage: K.<a> = GF(5^3).over()
            sage: L.<b> = GF(5^9).over(K)
            sage: u = 1/(a+b)

            sage: u._latex_extension(base=K)
            \left( 2 + 2 a \right) + \left( -1 + a - a^2 \right) b + \left( 2 + 3 a + 3 a^2 \right) b^2
            sage: u._latex_extension(base=GF(5))
            2 + 2 a - b + ab - a^2b + 2 b^2 + 3 ab^2 + 3 a^2b^2

        """
        parent = self._parent
        coeffs = self._vector(base)
        names = parent._basis_latex_names
        b = parent._base
        while b is not base:
            names = [ x + y for y in names for x in b._basis_latex_names ]
            b = b._base
        s = ""
        for i in range(len(names)):
            c = coeffs[i]
            if c.is_zero():
                continue
            sign = 1
            ss = ""
            if c == -1:
                sign = -1
            elif c != 1:
                atomic = c._is_atomic()
                if not atomic and (-c)._is_atomic():
                    c = -c
                    sign = -sign
                    atomic = True
                sc = latex(c)
                if atomic:
                    ss += sc
                else:
                    ss += r"\left(" + sc + r"\right)"
            if ss != "" and ss[0] == "-":
                ss = ss[1:]
                sign *= -1
            if s == "":
                if sign == -1:
                    s = "-"
            else:
                s += " + " if sign == 1 else " - "
            ss += names[i]
            if ss == "":
                ss += "1"
            s += ss
        if s == "":
            return "0"
        if s[:6] == r"\left(" and s[-7] == r"\right)":
            s = s[6:-7]
        return s

    def vector(self, base=None):
        r"""
        Return the vector of coordinates of this element over ``base``
        (in the basis output by the method :meth:`basis_over`).

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None``

        EXAMPLES::

            sage: F = GF(5)
            sage: K.<a> = GF(5^2).over()  # over F
            sage: L.<b> = GF(5^6).over(K)
            sage: x = (a+b)^4; x
            (-1 + a) + (3 + a)*b + (1 - a)*b^2

            sage: x.vector(K)  # basis is (1, b, b^2)
            (-1 + a, 3 + a, 1 - a)

            sage: x.vector(F)  # basis is (1, a, b, a*b, b^2, a*b^2)
            (4, 1, 3, 1, 1, 4)

        If ``base`` is omitted, it is set to its default which is the
        base of the extension::

            sage: x.vector()
            (-1 + a, 3 + a, 1 - a)

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: x.vector(GF(5^3))
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field in z3 of size 5^3
        """
        base = self._parent._check_base(base)
        return self._vector(base)

    def _vector_(self, reverse=False, base=None):
        """
        Return the vector of coordinates of this element over a base ring.

        The arguments are chosen for compatibility with the analogous methods of finite fields.

        INPUT:

        - ``reverse`` -- a boolean, whether to reverse the coordinates of the output
        - ``base`` -- a commutative ring; if specified it must be an explicit base over which
            this extension has been defined.

        TESTS::

            sage: K = GF(11^10).over(GF(11^2))
            sage: x = K.random_element()
            sage: coeffs = x._vector_()
            sage: rcoeffs = x._vector_(reverse=True)
            sage: list(coeffs) == list(reversed(rcoeffs))
            True
        """
        v = self.vector(base)
        if reverse:
            v = v.parent()(list(reversed(v)))
        return v

    cdef _vector(self, CommutativeRing base):
        r"""
        Return the vector of coordinates of this element over ``base``
        (in the basis output by the method :meth:`basis_over`).

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None``

        TESTS::

            sage: K = GF(11^10).over(GF(11^2))
            sage: x = K.random_element()
            sage: coeffs = x.vector()
            sage: basis = K.basis_over()
            sage: x == sum(coeffs[i]*basis[i] for i in range(5))
            True
        """
        _, _, j = self._parent._free_module(base, map=True)
        return j(self)

    def polynomial(self, var='x', base=None):
        r"""
        Return a polynomial (in one or more variables) over ``base``
        whose evaluation at the generators of the parent equals this
        element.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None``

        EXAMPLES::

            sage: F.<a> = GF(5^2).over()  # over GF(5)
            sage: K.<b> = GF(5^4).over(F)
            sage: L.<c> = GF(5^12).over(K)
            sage: u = 1/(a + b + c); u
            (2 + (-1 - a)*b) + ((2 + 3*a) + (1 - a)*b)*c + ((-1 - a) - a*b)*c^2

            sage: P = u.polynomial(base=K); P
            ((-1 - a) - a*b)*x^2 + ((2 + 3*a) + (1 - a)*b)*x + 2 + (-1 - a)*b
            sage: P.base_ring() is K
            True
            sage: P(c) == u
            True

        When the base is `F`, we obtain a bivariate polynomial::

            sage: P = u.polynomial(base=F); P
            (-a)*x0^2*x1 + (-1 - a)*x0^2 + (1 - a)*x0*x1 + (2 + 3*a)*x0 + (-1 - a)*x1 + 2

        We check that its value at the generators is the element we started with::

            sage: L.gens(F)
            (c, b)
            sage: P(c, b) == u
            True

        Similarly, when the base is ``GF(5)``, we get a trivariate polynomial:

            sage: P = u.polynomial(base=GF(5)); P
            -x0^2*x1*x2 - x0^2*x2 - x0*x1*x2 - x0^2 + x0*x1 - 2*x0*x2 - x1*x2 + 2*x0 - x1 + 2
            sage: P(c, b, a) == u
            True

        Different variable names can be specified::

            sage: u.polynomial(base=GF(5), var='y')
            -y0^2*y1*y2 - y0^2*y2 - y0*y1*y2 - y0^2 + y0*y1 - 2*y0*y2 - y1*y2 + 2*y0 - y1 + 2
            sage: u.polynomial(base=GF(5), var=['x','y','z'])
            -x^2*y*z - x^2*z - x*y*z - x^2 + x*y - 2*x*z - y*z + 2*x - y + 2

        If ``base`` is omitted, it is set to its default which is the
        base of the extension::

            sage: u.polynomial()
            ((-1 - a) - a*b)*x^2 + ((2 + 3*a) + (1 - a)*b)*x + 2 + (-1 - a)*b

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: u.polynomial(base=GF(5^3))
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field in z3 of size 5^3
        """
        base = self._parent._check_base(base)
        degrees = [ ]
        b = self._parent
        degree = 1
        from sage.rings.ring_extension import RingExtensionWithGen
        while b is not base:
            if not isinstance(b, RingExtensionWithGen):
                raise NotImplementedError
            reldeg = b.relative_degree()
            degree *= reldeg
            degrees.append(reldeg)
            b = b.base_ring()
        degrees.reverse()
        coeffs = { }
        v = self._vector(base)
        S = PolynomialRing(base, len(degrees), names=var)
        for i in range(degree):
            ii = ZZ(i)
            exponents = [ ]
            for d in degrees:
                ii, exponent = ii.quo_rem(d)
                exponents.append(exponent)
            coeffs[tuple(reversed(exponents))] = v[i]
        return S(coeffs)

    def matrix(self, base=None):
        r"""
        Return the matrix of the multiplication by this element (in
        the basis output by :meth:`basis_over`).

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None``

        EXAMPLES::

            sage: K.<a> = GF(5^3).over()  # over GF(5)
            sage: L.<b> = GF(5^6).over(K)
            sage: u = a/(1+b)

            sage: u
            (2 + a + 3*a^2) + (3 + 3*a + a^2)*b
            sage: b*u
            (3 + 2*a^2) + (2 + 2*a - a^2)*b
            sage: u.matrix(K)
            [2 + a + 3*a^2 3 + 3*a + a^2]
            [    3 + 2*a^2 2 + 2*a - a^2]

            sage: u.matrix(GF(5))
            [2 1 3 3 3 1]
            [1 3 1 2 0 3]
            [2 3 3 1 3 0]
            [3 0 2 2 2 4]
            [4 2 0 3 0 2]
            [0 4 2 4 2 0]

        If ``base`` is omitted, it is set to its default which is the
        base of the extension::

            sage: u.matrix()
            [2 + a + 3*a^2 3 + 3*a + a^2]
            [    3 + 2*a^2 2 + 2*a - a^2]

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: u.matrix(GF(5^2))
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field in z2 of size 5^2
        """
        parent = self._parent
        base = parent._check_base(base)
        if not (parent._is_finite_over(base) and parent._is_free_over(base)):
            raise ValueError("the extension is not finite free")
        return self._matrix(base)

    cdef _matrix(self, CommutativeRing base):
        r"""
        Return the matrix of the multiplication by this element (in
        the basis output by :meth:`basis_over`).

        This method does not check its input.
        Do not call it directly; use :meth:`matrix` instead.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        TESTS::

            sage: F = GF(11^2)
            sage: K = GF(11^6).over(F)
            sage: L = GF(11^18).over(K)

            sage: for base in L.bases():
            ....:     x = L.random_element()
            ....:     y = L.random_element()
            ....:     assert((x+y).matrix(base) == x.matrix(base) + y.matrix(base))
            ....:     assert((x*y).matrix(base) == x.matrix(base) * y.matrix(base))
        """
        from sage.matrix.matrix_space import MatrixSpace
        parent = self._parent
        _, _, j = parent._free_module(base, map=True)
        x = self._backend
        M = [ j(parent._from_backend_morphism(x * (<RingExtensionElement>b)._backend)) for b in parent._basis_over(base) ]
        return MatrixSpace(base, len(M))(M)

    def trace(self, base=None):
        r"""
        Return the trace of this element over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None``

        EXAMPLES::

            sage: F = GF(5)
            sage: K.<a> = GF(5^3).over(F)
            sage: L.<b> = GF(5^6).over(K)
            sage: u = a/(1+b)

            sage: tr = u.trace(K); tr
            -1 + 3*a + 2*a^2

        We check that the trace lives in the base ring::

            sage: tr.parent()
            Field in a with defining polynomial x^3 + 3*x + 3 over its base
            sage: tr.parent() is K
            True

        Similarly, one can compute the trace over F::

            sage: u.trace(F)
            0

        We check the transitivity of the trace::

            sage: u.trace(F) == tr.trace(F)
            True

        If ``base`` is omitted, it is set to its default which is the
        base of the extension::

            sage: u.trace()
            -1 + 3*a + 2*a^2

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: u.trace(GF(5^2))
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field in z2 of size 5^2
        """
        parent = self._parent
        base = parent._check_base(base)
        if not (parent._is_finite_over(base) and parent._is_free_over(base)):
            raise ValueError("the extension is not finite free")
        return self._trace(base)

    cdef _trace(self, CommutativeRing base):
        r"""
        Return the trace of this element over ``base``.

        This method does not check its input.
        Do not call it directly; use :meth:`trace` instead.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        TESTS::

            sage: F = GF(11^2)
            sage: K = GF(11^6).over(F)
            sage: L = GF(11^18).over(K)

            sage: x = L.random_element()
            sage: x.trace(F) == x.trace().trace()
            True

            sage: for base in L.bases():
            ....:     x = L.random_element()
            ....:     y = L.random_element()
            ....:     assert(x.trace(base) == x.matrix(base).trace())
            ....:     assert((x+y).trace(base) == x.trace(base) + y.trace(base))
        """
        parent = self._parent
        cdef CommutativeRing b
        if base is parent:
            return self
        b = parent._base
        t = self._matrix(b).trace()
        if base is b:
            return t
        return (<RingExtensionWithBasisElement>t)._trace(base)

    def norm(self, base=None):
        r"""
        Return the norm of this element over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None``

        EXAMPLES::

            sage: F = GF(5)
            sage: K.<a> = GF(5^3).over(F)
            sage: L.<b> = GF(5^6).over(K)
            sage: u = a/(1+b)

            sage: nr = u.norm(K); nr
            3 + 2*a^2

        We check that the norm lives in the base ring::

            sage: nr.parent()
            Field in a with defining polynomial x^3 + 3*x + 3 over its base
            sage: nr.parent() is K
            True

        Similarly, one can compute the norm over F::

            sage: u.norm(F)
            4

        We check the transitivity of the norm::

            sage: u.norm(F) == nr.norm(F)
            True

        If ``base`` is omitted, it is set to its default which is the
        base of the extension::

            sage: u.norm()
            3 + 2*a^2

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: u.norm(GF(5^2))
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field in z2 of size 5^2
        """
        parent = self._parent
        base = parent._check_base(base)
        if not (parent._is_finite_over(base) and parent._is_free_over(base)):
            raise ValueError("the extension is not finite free")
        return self._norm(base)

    cdef _norm(self, CommutativeRing base):
        r"""
        Return the norm of this element over ``base``.

        This method does not check its input.
        Do not call it directly; use :meth:`norm` instead.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension)

        TESTS::

            sage: F = GF(11^2)
            sage: K = GF(11^6).over(F)
            sage: L = GF(11^18).over(K)

            sage: x = L.random_element()
            sage: x.norm(F) == x.norm().norm()
            True

            sage: for base in L.bases():
            ....:     x = L.random_element()
            ....:     y = L.random_element()
            ....:     assert(x.norm(base) == x.matrix(base).determinant())
            ....:     assert((x*y).norm(base) == x.norm(base) * y.norm(base))
        """
        parent = self._parent
        cdef CommutativeRing b
        if base is parent:
            return self
        b = parent._base
        n = self._matrix(b).determinant()
        if base is b:
            return n
        return (<RingExtensionWithBasisElement>n)._norm(base)

    def charpoly(self, var='x', base=None):
        r"""
        Return the characteristic polynomial of this element over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None``

        EXAMPLES::

            sage: F = GF(5)
            sage: K.<a> = GF(5^3).over(F)
            sage: L.<b> = GF(5^6).over(K)
            sage: u = a/(1+b)

            sage: chi = u.charpoly(base=K); chi
            x^2 + (1 + 2*a + 3*a^2)*x + 3 + 2*a^2

        We check that the charpoly has coefficients in the base ring::

            sage: chi.base_ring()
            Field in a with defining polynomial x^3 + 3*x + 3 over its base
            sage: chi.base_ring() is K
            True

        and that it annihilates u::

            sage: chi(u)
            0

        Similarly, one can compute the characteristic polynomial over F::

            sage: u.charpoly(base=F)
            x^6 + x^4 + 2*x^3 + 3*x + 4

        A different variable name can be specified::

            sage: u.charpoly(base=F, var='t')
            t^6 + t^4 + 2*t^3 + 3*t + 4

        If ``base`` is omitted, it is set to its default which is the
        base of the extension::

            sage: u.charpoly()
            x^2 + (1 + 2*a + 3*a^2)*x + 3 + 2*a^2

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: u.charpoly(base=GF(5^2))
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field in z2 of size 5^2

        TESTS:

        We check that the characteristic polynomial of an element in the base
        ring is a power of a polynomial of degree 1::

            sage: S.<x> = K[]
            sage: u = K.random_element()
            sage: L(u).charpoly() == (x - u)^2
            True
        """
        return self.matrix(base).charpoly(var)

    cpdef minpoly(self, var='x', base=None):
        r"""
        Return the minimal polynomial of this element over ``base``.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None``

        EXAMPLES::

            sage: F = GF(5)
            sage: K.<a> = GF(5^3).over(F)
            sage: L.<b> = GF(5^6).over(K)
            sage: u = 1 / (a+b)

            sage: chi = u.minpoly(base=K); chi
            x^2 + (2*a + a^2)*x - 1 + a

        We check that the minimal polynomial has coefficients in the base ring::

            sage: chi.base_ring()
            Field in a with defining polynomial x^3 + 3*x + 3 over its base
            sage: chi.base_ring() is K
            True

        and that it annihilates u::

            sage: chi(u)
            0

        Similarly, one can compute the minimal polynomial over F::

            sage: u.minpoly(base=F)
            x^6 + 4*x^5 + x^4 + 2*x^2 + 3

        A different variable name can be specified::

            sage: u.minpoly(base=F, var='t')
            t^6 + 4*t^5 + t^4 + 2*t^2 + 3

        If ``base`` is omitted, it is set to its default which is the
        base of the extension::

            sage: u.minpoly()
            x^2 + (2*a + a^2)*x - 1 + a

        Note that ``base`` must be an explicit base over which the
        extension has been defined (as listed by the method :meth:`bases`)::

            sage: u.minpoly(base=GF(5^2))
            Traceback (most recent call last):
            ...
            ValueError: not (explicitly) defined over Finite Field in z2 of size 5^2

        TESTS:

        We check that the minimal polynomial of an element in the base
        ring has degree 1::

            sage: S.<x> = K[]
            sage: u = K.random_element()
            sage: L(u).minpoly() == x - u
            True

        In a similar fashion, the minimal polynomial over `F` of an element
        of `K` should have degree 1 or 3::

            sage: L(u).minpoly(base=F).degree() in [ 1, 3 ]
            True
        """
        from sage.modules.free_module import FreeModule
        if not isinstance(var, str):
            raise ValueError("Variable name must be a string, but it is %s" % (type(var)))
        from sage.modules.free_module import FreeModule
        parent = self._parent
        if base is None:
            mod = parent.modulus()
            S = mod.parent().quotient(mod)
            try:
                return S(list(self.vector())).minpoly()
            except NotImplementedError:
                pass  # fall back to generic code below

        base = parent._check_base(base)
        if not (parent._is_finite_over(base) and parent._is_free_over(base)):
            raise ValueError("the extension is not finite free")
        if base not in Fields():
            raise NotImplementedError("minpoly is only implemented when the base is a field")
        base_backend, from_base_backend, _ = backend_parent(base, map=True)
        degree = parent._degree_over(base)
        _, _, j = parent._free_module(base, map=True)
        V = FreeModule(base_backend, degree)
        vector = [base_backend(1)] + (degree-1)*[base_backend(0)]
        vectors = [vector]
        W = V.span(vectors)
        elt = self
        while True:
            vector = V(j.backend_coefficients(elt))
            if vector in W: break
            vectors.append(vector)
            W += V.span([vector])
            elt *= self
        W = V.span_of_basis(vectors)
        coeffs = [ -c for c in W.coordinate_vector(vector) ] + [base_backend(1)]
        coeffs = [ from_base_backend(c) for c in coeffs ]
        return PolynomialRing(base, name=var)(coeffs)

    def minimal_polynomial(self, var='x', base=None):
        """
        Return the minimal polynomial of this element over a specified base ring.

        INPUT:

        - ``base`` -- a commutative ring (which might be itself an
          extension) or ``None`` (the base of this ring extension)

        .. SEEALSO::

            :meth:`minpoly`

        EXAMPLES::

            sage: k = GF(3^2); R.<x> = k[]
            sage: l.<b> = GF(3^6).over(k)
            sage: (b+1).minimal_polynomial()
            x^3 + (2*z2 + 2)*x^2 + 2*x + 1
        """
        return self.minpoly(var, base)

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        The image of this element under the morphism defined by
        ``im_gens`` in ``codomain``, where elements of the
        base ring are mapped by ``base_map``.

        INPUT:

        - ``codomain`` -- another ring

        - ``im_gens`` -- a list of elements in the codomain, giving the images of the generators of the parent

        - ``base_map`` -- a homomorphism on the base ring

        EXAMPLES::

            sage: A.<a> = QQ.extension(x^2 - 2)
            sage: B.<b> = QQ.extension(x^6 - 2)
            sage: f = A.hom([b^3])
            sage: E = B.over(f)
            sage: c = b^2 + b + 1
            sage: R.<x> = QQ[]
            sage: K.<z> = NumberField(x^6 - 2)
            sage: c._im_gens_(K, [-z], base_map=A.hom([-z^3]))
            z^2 - z + 1
        """
        if base_map is None:
            R = self._parent.base_ring()
            S = codomain.base_ring()
            if S.has_coerce_map_from(R):
                base_map = S.coerce_map_from(R)
            elif codomain.has_coerce_map_from(R):
                base_map = codomain.coerce_map_from(R)
            else:
                raise ValueError(f"There is no coercion from {R} to {codomain}, so you must explicitly give a base map")
        from sage.rings.ring_extension import RingExtensionWithGen
        if isinstance(self._parent, RingExtensionWithGen) and len(im_gens) == 1:
            f = self.polynomial().change_ring(base_map)
            return f(im_gens[0])
        elif len(im_gens) != self._parent.ngens():
            raise ValueError(f"{len(im_gens)} images given but parent has {self._parent.ngens()} generators")
        else:
            v = self.vector()
            ans = codomain(0)
            for c, g in zip(v, im_gens):
                ans += base_map(c) * g
            return ans
