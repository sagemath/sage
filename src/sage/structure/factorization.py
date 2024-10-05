# sage_setup: distribution = sagemath-objects
r"""
Factorizations

The :class:`Factorization` class provides a structure for holding quite
general lists of objects with integer multiplicities.  These may hold
the results of an arithmetic or algebraic factorization, where the
objects may be primes or irreducible polynomials and the
multiplicities are the (nonzero) exponents in the factorization.  For
other types of examples, see below.

:class:`Factorization` class objects contain a ``list``, so can be
printed nicely and be manipulated like a list of prime-exponent pairs,
or easily turned into a plain list.  For example, we factor the
integer `-45`::

    sage: F = factor(-45)

This returns an object of type :class:`Factorization`::

    sage: type(F)
    <class 'sage.structure.factorization_integer.IntegerFactorization'>

It prints in a nice factored form::

    sage: F
    -1 * 3^2 * 5

There is an underlying list representation, which ignores the unit part::

    sage: list(F)
    [(3, 2), (5, 1)]

A :class:`Factorization` is not actually a list::

    sage: isinstance(F, list)
    False

However, we can access the :class:`Factorization` F itself as if it were a list::

    sage: F[0]
    (3, 2)
    sage: F[1]
    (5, 1)

To get at the unit part, use the :meth:`Factorization.unit` function::

    sage: F.unit()
    -1

All factorizations are immutable, up to ordering with ``sort()`` and
simplifying with ``simplify()``.  Thus if you write a function that
returns a cached version of a factorization, you do not have to return
a copy.

::

    sage: F = factor(-12); F
    -1 * 2^2 * 3
    sage: F[0] = (5,4)
    Traceback (most recent call last):
    ...
    TypeError: 'Factorization' object does not support item assignment

EXAMPLES:

This more complicated example involving polynomials also illustrates
that the unit part is not discarded from factorizations::

    sage: # needs sage.libs.pari
    sage: x = QQ['x'].0
    sage: f = -5*(x-2)*(x-3)
    sage: f
    -5*x^2 + 25*x - 30
    sage: F = f.factor(); F
    (-5) * (x - 3) * (x - 2)
    sage: F.unit()
    -5
    sage: F.value()
    -5*x^2 + 25*x - 30

The underlying list is the list of pairs `(p_i, e_i)`, where each
`p_i` is a 'prime' and each `e_i` is an integer. The unit part
is discarded by the list::

    sage: # needs sage.libs.pari
    sage: list(F)
    [(x - 3, 1), (x - 2, 1)]
    sage: len(F)
    2
    sage: F[1]
    (x - 2, 1)

In the ring `\ZZ[x]`, the integer `-5` is not a unit, so the
factorization has three factors::

    sage: # needs sage.libs.pari
    sage: x = ZZ['x'].0
    sage: f = -5*(x-2)*(x-3)
    sage: f
    -5*x^2 + 25*x - 30
    sage: F = f.factor(); F
    (-1) * 5 * (x - 3) * (x - 2)
    sage: F.universe()
    Univariate Polynomial Ring in x over Integer Ring
    sage: F.unit()
    -1
    sage: list(F)
    [(5, 1), (x - 3, 1), (x - 2, 1)]
    sage: F.value()
    -5*x^2 + 25*x - 30
    sage: len(F)
    3

On the other hand, -1 is a unit in `\ZZ`, so it is included in the unit::

    sage: # needs sage.libs.pari
    sage: x = ZZ['x'].0
    sage: f = -1 * (x-2) * (x-3)
    sage: F = f.factor(); F
    (-1) * (x - 3) * (x - 2)
    sage: F.unit()
    -1
    sage: list(F)
    [(x - 3, 1), (x - 2, 1)]

Factorizations can involve fairly abstract mathematical objects::

    sage: # needs sage.modular
    sage: F = ModularSymbols(11,4).factorization(); F
    (Modular Symbols subspace of dimension 2 of Modular Symbols space
      of dimension 6 for Gamma_0(11) of weight 4 with sign 0 over Rational Field) *
    (Modular Symbols subspace of dimension 2 of Modular Symbols space
      of dimension 6 for Gamma_0(11) of weight 4 with sign 0 over Rational Field) *
    (Modular Symbols subspace of dimension 2 of Modular Symbols space
      of dimension 6 for Gamma_0(11) of weight 4 with sign 0 over Rational Field)
    sage: type(F)
    <class 'sage.structure.factorization.Factorization'>


    sage: # needs sage.rings.number_field
    sage: x = ZZ['x'].0
    sage: K.<a> = NumberField(x^2 + 3); K
    Number Field in a with defining polynomial x^2 + 3
    sage: f = K.factor(15); f
    (Fractional ideal (-a))^2 * (Fractional ideal (5))
    sage: f.universe()
    Monoid of ideals of Number Field in a with defining polynomial x^2 + 3
    sage: f.unit()
    Fractional ideal (1)
    sage: g = K.factor(9); g
    (Fractional ideal (-a))^4
    sage: f.lcm(g)
    (Fractional ideal (-a))^4 * (Fractional ideal (5))
    sage: f.gcd(g)
    (Fractional ideal (-a))^2
    sage: f.is_integral()
    True

TESTS::

    sage: F = factor(-20); F
    -1 * 2^2 * 5
    sage: G = loads(dumps(F)); G
    -1 * 2^2 * 5
    sage: G == F
    True
    sage: G is F
    False

AUTHORS:

- William Stein (2006-01-22): added unit part as suggested by David Kohel.

- William Stein (2008-01-17): wrote much of the documentation and
  fixed a couple of bugs.

- Nick Alexander (2008-01-19): added support for non-commuting factors.

- John Cremona (2008-08-22): added division, lcm, gcd, is_integral and
  universe functions
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.structure.element import Element
from sage.structure.sequence import Sequence
from sage.structure.richcmp import richcmp_method, richcmp, richcmp_not_equal
from sage.misc.cachefunc import cached_method


@richcmp_method
class Factorization(SageObject):
    """
    A formal factorization of an object.

    EXAMPLES::

        sage: N = 2006
        sage: F = N.factor(); F
        2 * 17 * 59
        sage: F.unit()
        1
        sage: F = factor(-2006); F
        -1 * 2 * 17 * 59
        sage: F.unit()
        -1
        sage: loads(F.dumps()) == F
        True
        sage: F = Factorization([(x, 1/3)])                                             # needs sage.symbolic
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer
    """
    def __init__(self, x, unit=None, cr=False, sort=True, simplify=True):
        """
        Create a :class:`Factorization` object.

        INPUT:

        - ``x`` -- list of pairs (p, e) with e an integer
          otherwise a :exc:`TypeError` is raised

        - ``unit`` -- (default: 1) the unit part of the factorization

        - ``cr`` -- (default: ``False``) if ``True``, print the factorization
          with carriage returns between factors

        - ``sort`` -- (default: ``True``) if ``True``, sort the factors by
          calling the sort function ``self.sort()`` after creating
          the factorization

        - ``simplify`` -- (default: ``True``) if ``True``, remove duplicate
          factors from the factorization.  See the documentation for
          self.simplify.

        OUTPUT: a Factorization object

        EXAMPLES:

        We create a factorization with all the default options::

            sage: Factorization([(2,3), (5, 1)])
            2^3 * 5

        We create a factorization with a specified unit part::

            sage: Factorization([(2,3), (5, 1)], unit=-1)
            -1 * 2^3 * 5

        We try to create a factorization but with a string an exponent, which
        results in a TypeError::

            sage: Factorization([(2,3), (5, 'x')])
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'x' to an integer

        We create a factorization that puts newlines after each multiply sign
        when printing.  This is mainly useful when the primes are large::

            sage: Factorization([(2,3), (5, 2)], cr=True)
            2^3 *
            5^2

        Another factorization with newlines and nontrivial unit part, which
        appears on a line by itself::

            sage: Factorization([(2,3), (5, 2)], cr=True, unit=-2)
            -2 *
            2^3 *
            5^2

        A factorization, but where we do not sort the factors::

            sage: Factorization([(5,3), (2, 3)], sort=False)
            5^3 * 2^3

        By default, in the commutative case, factorizations are sorted by the
        prime base::

            sage: Factorization([(2, 7), (5,2), (2, 5)])
            2^12 * 5^2
            sage: R.<a,b> = FreeAlgebra(QQ, 2)                                          # needs sage.combinat sage.modules
            sage: Factorization([(a,1), (b,1), (a,2)])                                  # needs sage.combinat sage.modules
            a * b * a^2

        Autosorting (the default) swaps around the factors below::

            sage: F = Factorization([(ZZ^3, 2), (ZZ^2, 5)], cr=True); F                 # needs sage.modules
            (Ambient free module of rank 2 over the principal ideal domain Integer Ring)^5 *
            (Ambient free module of rank 3 over the principal ideal domain Integer Ring)^2
        """
        from sage.rings.integer import Integer
        x = [(p, Integer(e)) for (p, e) in x]

        try:
            self.__universe = Sequence(t[0] for t in x).universe()
        except TypeError:
            self.__universe = None

        self.__x = [(t[0], int(t[1])) for t in x]
        if unit is None:
            if x:
                try:
                    unit = self.__universe(1)
                except (AttributeError, TypeError):
                    unit = Integer(1)
            else:
                unit = Integer(1)
        self.__unit = unit
        self.__cr = cr
        if sort and self.is_commutative():
            self.sort()
        if simplify:
            self.simplify()

    def __getitem__(self, i):
        """
        Return `i`-th factor of ``self``.

        EXAMPLES::

            sage: a = factor(-75); a
            -1 * 3 * 5^2
            sage: a[0]
            (3, 1)
            sage: a[1]
            (5, 2)
            sage: a[-1]
            (5, 2)
            sage: a[5]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self.__x[i]

    def __setitem__(self, i, v):
        """
        Set the `i`-th factor of ``self``.

        .. warning::

           NOT ALLOWED -- Factorizations are immutable.

        EXAMPLES::

            sage: a = factor(-75); a
            -1 * 3 * 5^2
            sage: a[0] = (2,3)
            Traceback (most recent call last):
            ...
            TypeError: 'Factorization' object does not support item assignment
        """
        raise TypeError("'Factorization' object does not support item assignment")

    def __len__(self):
        """
        Return the number of prime factors of ``self``, not counting
        the unit part.

        EXAMPLES::

            sage: len(factor(15))
            2

        Note that the unit part is not included in the count::

            sage: a = factor(-75); a
            -1 * 3 * 5^2
            sage: len(a)
            2
            sage: list(a)
            [(3, 1), (5, 2)]
            sage: len(list(a))
            2
        """
        return len(self.__x)

    def __richcmp__(self, other, op):
        """
        Compare ``self`` and ``other``.

        This first compares the values.

        If values are equal, this compares the units.

        If units are equal, this compares the underlying lists of
        ``self`` and ``other``.

        EXAMPLES:

        We compare two contrived formal factorizations::

            sage: a = Factorization([(2, 7), (5,2), (2, 5)])
            sage: b = Factorization([(2, 7), (5,10), (7, 3)])
            sage: a
            2^12 * 5^2
            sage: b
            2^7 * 5^10 * 7^3
            sage: a < b
            True
            sage: b < a
            False
            sage: a.value()
            102400
            sage: b.value()
            428750000000

        We compare factorizations of some polynomials::

            sage: x = polygen(QQ)
            sage: x^2 - 1 > x^2 - 4
            True
            sage: factor(x^2 - 1) > factor(x^2 - 4)                                     # needs sage.libs.pari
            True
        """
        if not isinstance(other, Factorization):
            return NotImplemented

        lx = self.value()
        rx = other.value()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx = self.__unit
        rx = other.__unit
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        return richcmp(self.__x, other.__x, op)

    def __copy__(self):
        r"""
        Return a copy of ``self``.

        This is *not* a deepcopy -- only references to the factors are
        returned, not copies of them.  Use ``deepcopy(self)`` if you need
        a deep copy of ``self``.

        EXAMPLES:

        We create a factorization that has mutable primes::

            sage: F = Factorization([([1,2], 5), ([5,6], 10)]); F
            ([1, 2])^5 * ([5, 6])^10

        We make a copy of it::

            sage: G = copy(F); G
            ([1, 2])^5 * ([5, 6])^10
            sage: G is F
            False

        Note that if we change one of the mutable "primes" of F, this does
        change G::

            sage: F[1][0][0] = 'hello'
            sage: G
            ([1, 2])^5 * (['hello', 6])^10
        """
        # No need to sort, since the factorization is already sorted
        # in whatever order is desired.
        return Factorization(self.__x, unit=self.__unit, cr=self.__cr,
                             sort=False, simplify=False)

    def __deepcopy__(self, memo):
        r"""
        Return a deep copy of ``self``.

        EXAMPLES:

        We make a factorization that has mutable entries::

            sage: F = Factorization([([1,2], 5), ([5,6], 10)]); F
            ([1, 2])^5 * ([5, 6])^10

        Now we make a copy of it and a deep copy::

            sage: K = copy(F)
            sage: G = deepcopy(F); G
            ([1, 2])^5 * ([5, 6])^10

        We change one of the mutable entries of F::

            sage: F[0][0][0] = 10

        This of course changes F::

            sage: F
            ([10, 2])^5 * ([5, 6])^10

        It also changes the copy K of F::

            sage: K
            ([10, 2])^5 * ([5, 6])^10

        It does *not* change the deep copy G::

            sage: G
            ([1, 2])^5 * ([5, 6])^10
        """
        import copy
        return Factorization(copy.deepcopy(list(self), memo),
                             cr=self.__cr, sort=False, simplify=False)

    def universe(self):
        r"""
        Return the parent structure of my factors.

        .. NOTE::

           This used to be called ``base_ring``, but the universe
           of a factorization need not be a ring.

        EXAMPLES::

            sage: F = factor(2006)
            sage: F.universe()
            Integer Ring

            sage: R.<x,y,z> = FreeAlgebra(QQ, 3)                                        # needs sage.combinat sage.modules
            sage: F = Factorization([(z, 2)], 3)                                        # needs sage.combinat sage.modules
            sage: (F*F^-1).universe()                                                   # needs sage.combinat sage.modules
            Free Algebra on 3 generators (x, y, z) over Rational Field

            sage: F = ModularSymbols(11,4).factorization()                              # needs sage.modular
            sage: F.universe()                                                          # needs sage.modular
        """
        try:
            return self.__universe
        except AttributeError:
            return None

    def base_change(self, U):
        """
        Return the factorization ``self``, with its factors (including the
        unit part) coerced into the universe `U`.

        EXAMPLES::

            sage: F = factor(2006)
            sage: F.universe()
            Integer Ring
            sage: P.<x> = ZZ[]
            sage: F.base_change(P).universe()
            Univariate Polynomial Ring in x over Integer Ring

        This method will return a :exc:`TypeError` if the coercion is not
        possible::

            sage: g = x^2 - 1
            sage: F = factor(g); F                                                      # needs sage.libs.pari
            (x - 1) * (x + 1)
            sage: F.universe()                                                          # needs sage.libs.pari
            Univariate Polynomial Ring in x over Integer Ring
            sage: F.base_change(ZZ)                                                     # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            TypeError: Impossible to coerce the factors of (x - 1) * (x + 1) into Integer Ring
        """
        if len(self) == 0:
            return self
        try:
            return Factorization([(U(f[0]), f[1]) for f in list(self)], unit=U(self.unit()))
        except TypeError:
            raise TypeError("Impossible to coerce the factors of %s into %s" % (self, U))

    def is_commutative(self) -> bool:
        """
        Return whether the factors commute.

        EXAMPLES::

            sage: F = factor(2006)
            sage: F.is_commutative()
            True

            sage: # needs sage.rings.number_field
            sage: K = QuadraticField(23, 'a')
            sage: F = K.factor(13)
            sage: F.is_commutative()
            True

            sage: # needs sage.combinat sage.modules
            sage: R.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: F = Factorization([(z, 2)], 3)
            sage: F.is_commutative()
            False
            sage: (F*F^-1).is_commutative()
            False
        """
        try:
            return self.universe().is_commutative()
        except Exception:
            # This is not the mathematically correct default, but agrees with
            # history -- we've always assumed factored things commute
            return True

    def _set_cr(self, cr):
        """
        Change whether or not the factorization is printed with
        carriage returns after each factor.

        EXAMPLES::

            sage: x = polygen(QQ,'x')
            sage: F = factor(x^6 - 1); F                                                # needs sage.libs.pari
            (x - 1) * (x + 1) * (x^2 - x + 1) * (x^2 + x + 1)
            sage: F._set_cr(True); F                                                    # needs sage.libs.pari
            (x - 1) *
            (x + 1) *
            (x^2 - x + 1) *
            (x^2 + x + 1)
            sage: F._set_cr(False); F                                                   # needs sage.libs.pari
            (x - 1) * (x + 1) * (x^2 - x + 1) * (x^2 + x + 1)
        """
        self.__cr = bool(cr)

    def simplify(self):
        """
        Combine adjacent products as much as possible.

        TESTS::

            sage: # needs sage.combinat sage.modules
            sage: R.<x,y> = FreeAlgebra(ZZ, 2)
            sage: F = Factorization([(x,3), (y, 2), (y,2)], simplify=False); F
            x^3 * y^2 * y^2
            sage: F.simplify(); F
            x^3 * y^4
            sage: F * Factorization([(y, -2)], 2)
            (2) * x^3 * y^2
        """
        repeat = False
        simp = []
        import itertools
        for obj, agroup in itertools.groupby(list(self), lambda x: x[0]):
            xs = list(agroup)
            if len(xs) > 1:
                repeat = True
            n = sum([x[1] for x in xs])
            if n != 0:
                simp.append((obj, n))
        self.__x[0:] = simp
        if repeat:
            self.simplify()

    def sort(self, key=None):
        r"""
        Sort the factors in this factorization.

        INPUT:

        - ``key`` -- (default: ``None``) comparison key

        OUTPUT: changes this factorization to be sorted (inplace)

        If ``key`` is ``None``, we determine the comparison key as
        follows:

        If the prime in the first factor has a dimension
        method, then we sort based first on *dimension* then on
        the exponent.

        If there is no dimension method, we next
        attempt to sort based on a degree method, in which case, we
        sort based first on *degree*, then exponent to break ties
        when two factors have the same degree, and if those match
        break ties based on the actual prime itself.

        Otherwise, we sort according to the prime itself.

        EXAMPLES:

        We create a factored polynomial::

            sage: x = polygen(QQ, 'x')
            sage: F = factor(x^3 + 1); F                                                # needs sage.libs.pari
            (x + 1) * (x^2 - x + 1)

        We sort it by decreasing degree::

            sage: F.sort(key=lambda x: (-x[0].degree(), x))                             # needs sage.libs.pari
            sage: F                                                                     # needs sage.libs.pari
            (x^2 - x + 1) * (x + 1)
        """
        if len(self) == 0:
            return

        if key is not None:
            self.__x.sort(key=key)
            return

        a = self.__x[0][0]
        sort_key = None
        if hasattr(a, 'dimension'):
            try:
                a.dimension()

                def sort_key(f):
                    return (f[0].dimension(), f[1], f[0])
            except (AttributeError, NotImplementedError, TypeError):
                pass
        elif hasattr(a, 'degree'):
            try:
                a.degree()

                def sort_key(f):
                    return (f[0].degree(), f[1], f[0])
            except (AttributeError, NotImplementedError, TypeError):
                pass

        if sort_key is None:

            def sort_key(f):
                return f[0]

        self.__x.sort(key=sort_key)

    def unit(self):
        r"""
        Return the unit part of this factorization.

        EXAMPLES:

        We create a polynomial over the real double field and factor it::

            sage: x = polygen(RDF, 'x')
            sage: F = factor(-2*x^2 - 1); F                                             # needs numpy
            (-2.0) * (x^2 + 0.5000000000000001)

        Note that the unit part of the factorization is `-2.0`::

            sage: F.unit()                                                              # needs numpy
            -2.0

            sage: F = factor(-2006); F
            -1 * 2 * 17 * 59
            sage: F.unit()
            -1
        """
        return self.__unit

    def _cr(self):
        """
        Return whether or not factorizations are printed with carriage
        returns between factors.

        EXAMPLES:

        Our first example involves factoring an integer::

            sage: F = factor(-93930); F
            -1 * 2 * 3 * 5 * 31 * 101
            sage: F._cr()
            False
            sage: F._set_cr(True)
            sage: F._cr()
            True

        This of course looks funny::

            sage: F
            -1 *
            2 *
            3 *
            5 *
            31 *
            101

        Next we factor a modular symbols space::

            sage: F = ModularSymbols(11).factor(); F                                    # needs sage.modular
            (Modular Symbols subspace of dimension 1 of ...) *
            (Modular Symbols subspace of dimension 1 of ...) *
            (Modular Symbols subspace of dimension 1 of ...)
        """
        try:
            return self.__cr
        except AttributeError:
            self.__cr = False
            return False

    def _repr_(self):
        """
        Return the string representation of this factorization.

        EXAMPLES::

            sage: f = factor(-100); f
            -1 * 2^2 * 5^2
            sage: f._repr_()
            '-1 * 2^2 * 5^2'

        Note that the default printing of a factorization can be overloaded
        using the rename method::

            sage: f.rename('factorization of -100')
            sage: f
            factorization of -100

        However ``_repr_`` always prints normally::

            sage: f._repr_()
            '-1 * 2^2 * 5^2'

        EXAMPLES::

           sage: x = polygen(QQ)
           sage: Factorization([(x-1,1), (x-2,2)])
           (x - 1) * (x - 2)^2
           sage: Factorization([(x + 1, -3)])
           (x + 1)^-3
        """
        cr = self._cr()
        if len(self) == 0:
            return repr(self.__unit)
        s = ''
        mul = ' * '
        if cr:
            mul += '\n'
        x = self.__x[0][0]
        try:
            atomic = (isinstance(x, int) or
                      self.universe()._repr_option('element_is_atomic'))
        except AttributeError:
            atomic = False

        if isinstance(x, Element):
            one = x.parent()(1)
        else:
            one = 1

        for i in range(len(self)):
            t = repr(self.__x[i][0])
            n = self.__x[i][1]
            if not atomic and (n != 1 or len(self) > 1 or self.__unit != one):
                if '+' in t or '-' in t or ' ' in t:
                    t = '(%s)' % t
            if n != 1:
                t += '^%s' % n
            s += t
            if i < len(self) - 1:
                s += mul
        if self.__unit != one:
            if atomic:
                u = repr(self.__unit)
            else:
                u = '(%s)' % self.__unit
            s = u + mul + s
        return s

    def _latex_(self):
        r"""
        Return the LaTeX representation of this factorization.

        EXAMPLES::

            sage: f = factor(-100); f
            -1 * 2^2 * 5^2
            sage: latex(f)
            -1 \cdot 2^{2} \cdot 5^{2}
            sage: f._latex_()
            '-1 \\cdot 2^{2} \\cdot 5^{2}'
            sage: x = AA['x'].0; factor(x^2 + x + 1)._latex_()  # Issue #12178          # needs sage.rings.number_field
            '(x^{2} + x + 1.000000000000000?)'
        """
        if len(self) == 0:
            return self.__unit._latex_()
        try:
            atomic = (isinstance(self.__x[0][0], int) or
                      self.universe()._repr_option('element_is_atomic'))
        except AttributeError:
            atomic = False
        s = ''
        for i in range(len(self)):
            t = self.__x[i][0]._latex_()
            if not atomic and ('+' in t or '-' in t or ' ' in t):
                t = '(%s)' % t
            n = self.__x[i][1]
            if n != 1:
                t += '^{%s}' % n
            s += t
            if i < len(self) - 1:
                s += ' \\cdot '
        if self.__unit != 1:
            if atomic:
                u = self.__unit._latex_()
            else:
                u = '\\left(%s\\right)' % self.__unit._latex_()
            s = u + ' \\cdot ' + s
        return s

    @cached_method
    def __pari__(self):
        """
        Return the PARI factorization matrix corresponding to ``self``.

        EXAMPLES::

            sage: f = factor(-24)
            sage: pari(f)                                                               # needs sage.libs.pari
            [-1, 1; 2, 3; 3, 1]

            sage: R.<x> = QQ[]
            sage: g = factor(x^10 - 1)                                                  # needs sage.libs.pari
            sage: pari(g)                                                               # needs sage.libs.pari
            [x - 1, 1; x + 1, 1; x^4 - x^3 + x^2 - x + 1, 1; x^4 + x^3 + x^2 + x + 1, 1]
        """
        from sage.libs.pari.all import pari
        from itertools import chain

        n = len(self)
        if self.__unit == 1:
            init = ()
        else:
            init = (self.__unit, 1)
            n += 1
        # concatenate (p, e) tuples
        entries = init + tuple(chain.from_iterable(self))
        return pari.matrix(n, 2, entries)

    def __add__(self, other):
        """
        Return the (unfactored) sum of ``self`` and ``other``.

        EXAMPLES::

            sage: factor(-10) + 16
            6
            sage: factor(10) - 16
            -6
            sage: factor(100) + factor(19)
            119
        """
        if isinstance(other, Factorization):
            other = other.value()
        return self.value() + other

    def __sub__(self, other):
        """
        Return the (unfactored) difference of ``self`` and ``other``.

        EXAMPLES::

            sage: factor(-10) + 16
            6
            sage: factor(10) - 16
            -6
        """
        if isinstance(other, Factorization):
            other = other.value()
        return self.value() - other

    def __radd__(self, left):
        """
        Return the (unfactored) sum of ``self`` and ``left``.

        EXAMPLES::

            sage: 16 + factor(-10)
            6
        """
        return self.value() + left

    def __rsub__(self, left):
        """
        Return the (unfactored) difference of ``left`` and ``self``.

        EXAMPLES::

            sage: 16 - factor(10)
            6
        """
        return left - self.value()

    def __neg__(self):
        """
        Return negative of this factorization.

        EXAMPLES::

            sage: a = factor(-75); a
            -1 * 3 * 5^2
            sage: -a
            3 * 5^2
            sage: (-a).unit()
            1
        """
        unit = -self.__unit
        return Factorization(list(self), unit, self.__cr,
                             sort=False, simplify=False)

    def __rmul__(self, left):
        """
        Return the product ``left * self``, where ``left`` is not a Factorization.

        EXAMPLES::

            sage: a = factor(15); a
            3 * 5
            sage: -2 * a
            -2 * 3 * 5
            sage: a * -2
            -2 * 3 * 5
            sage: R.<x,y> = FreeAlgebra(QQ, 2)                                          # needs sage.combinat sage.modules
            sage: f = Factorization([(x,2), (y,3)]); f                                  # needs sage.combinat sage.modules
            x^2 * y^3
            sage: x * f                                                                 # needs sage.combinat sage.modules
            x^3 * y^3
            sage: f * x                                                                 # needs sage.combinat sage.modules
            x^2 * y^3 * x

        Note that this does not automatically factor ``left``::

            sage: F = Factorization([(5,3), (2,3)])
            sage: 46 * F
            2^3 * 5^3 * 46
        """
        return Factorization([(left, 1)]) * self

    def __mul__(self, other):
        r"""
        Return the product of two factorizations, which is obtained by
        combining together like factors.

        If the two factorizations have different universes, this
        method will attempt to find a common universe for the
        product.  A :exc:`TypeError` is raised if this is impossible.

        EXAMPLES::

            sage: factor(-10) * factor(-16)
            2^5 * 5
            sage: factor(-10) * factor(16)
            -1 * 2^5 * 5

            sage: # needs sage.combinat sage.modules
            sage: R.<x,y> = FreeAlgebra(ZZ, 2)
            sage: F = Factorization([(x,3), (y, 2), (x,1)]); F
            x^3 * y^2 * x
            sage: F*F
            x^3 * y^2 * x^4 * y^2 * x
            sage: -1 * F
            (-1) * x^3 * y^2 * x

            sage: P.<x> = ZZ[]
            sage: f = 2*x + 2
            sage: c = f.content(); g = f//c
            sage: Fc = factor(c); Fc.universe()
            Integer Ring
            sage: Fg = factor(g); Fg.universe()
            Univariate Polynomial Ring in x over Integer Ring
            sage: F = Fc * Fg; F.universe()
            Univariate Polynomial Ring in x over Integer Ring
            sage: [type(a[0]) for a in F]
            [<... 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>,
             <... 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>]
        """
        if not isinstance(other, Factorization):
            return self * Factorization([(other, 1)])

        if len(self) and len(other):
            try:
                # since self is a factorization, all its factors
                # are in the same universe.
                # the same is true for the factorization other.
                # so if we want to put the factorizations together we just
                # need to find a common universe for the first factor of
                # self and the first factor of other
                U = Sequence([self[0][0], other[0][0]]).universe()
                self = self.base_change(U)
                other = other.base_change(U)
            except TypeError:
                raise TypeError("Cannot multiply %s and %s because they cannot be coerced into a common universe" % (self, other))

        if self.is_commutative() and other.is_commutative():
            d1 = dict(self)
            d2 = dict(other)
            s = {}
            for a in set(d1).union(set(d2)):
                s[a] = d1.get(a, 0) + d2.get(a, 0)
            return Factorization(list(s.items()), unit=self.unit() * other.unit())
        else:
            return Factorization(list(self) + list(other), unit=self.unit() * other.unit())

    def __pow__(self, n):
        """
        Return the `n`-th power of a factorization, which is got by
        combining together like factors.

        EXAMPLES::

            sage: f = factor(-100); f
            -1 * 2^2 * 5^2
            sage: f^3
            -1 * 2^6 * 5^6
            sage: f^4
            2^8 * 5^8

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 39*x - 91)                                  # needs sage.rings.number_field
            sage: F = K.factor(7); F                                                    # needs sage.rings.number_field
            (Fractional ideal (7, a)) * (Fractional ideal (7, a + 2)) * (Fractional ideal (7, a - 2))
            sage: F^9                                                                   # needs sage.rings.number_field
            (Fractional ideal (7, a))^9 * (Fractional ideal (7, a + 2))^9 * (Fractional ideal (7, a - 2))^9

            sage: R.<x,y> = FreeAlgebra(ZZ, 2)                                          # needs sage.combinat sage.modules
            sage: F = Factorization([(x,3), (y, 2), (x,1)]); F                          # needs sage.combinat sage.modules
            x^3 * y^2 * x
            sage: F**2                                                                  # needs sage.combinat sage.modules
            x^3 * y^2 * x^4 * y^2 * x
        """
        from sage.rings.integer import Integer
        if not isinstance(n, Integer):
            try:
                n = Integer(n)
            except TypeError:
                raise TypeError("Exponent n (= %s) must be an integer." % n)
        if n == 1:
            return self
        if n == 0:
            return Factorization([])
        if self.is_commutative():
            return Factorization([(p, n * e) for p, e in self], unit=self.unit()**n,
                                 cr=self.__cr, sort=False, simplify=False)
        if n < 0:
            self = ~self
            n = -n
        from sage.arith.power import generic_power
        return generic_power(self, n)

    def __invert__(self):
        r"""
        Return the formal inverse of the factors in the factorization.

        EXAMPLES::

            sage: F = factor(2006); F
            2 * 17 * 59
            sage: F^-1
            2^-1 * 17^-1 * 59^-1

            sage: R.<x,y> = FreeAlgebra(QQ, 2)                                          # needs sage.combinat sage.modules
            sage: F = Factorization([(x,3), (y, 2), (x,1)], 2); F                       # needs sage.combinat sage.modules
            (2) * x^3 * y^2 * x
            sage: F^-1                                                                  # needs sage.combinat sage.modules
            (1/2) * x^-1 * y^-2 * x^-3
        """
        return Factorization([(p, -e) for p, e in reversed(self)],
                             cr=self._cr(), unit=self.unit()**(-1))

    def __truediv__(self, other):
        r"""
        Return the quotient of two factorizations, which is obtained by
        multiplying the first by the inverse of the second.

        EXAMPLES::

            sage: factor(-10) / factor(-16)
            2^-3 * 5
            sage: factor(-10) / factor(16)
            -1 * 2^-3 * 5

            sage: # needs sage.combinat sage.modules
            sage: R.<x,y> = FreeAlgebra(QQ, 2)
            sage: F = Factorization([(x,3), (y, 2), (x,1)]); F
            x^3 * y^2 * x
            sage: G = Factorization([(y, 1), (x,1)],1); G
            y * x
            sage: F / G
            x^3 * y
        """
        if not isinstance(other, Factorization):
            return self / Factorization([(other, 1)])
        return self * other**-1

    def __call__(self, *args, **kwds):
        """
        Implement the substitution.

        This is assuming that each term can be substituted.

        There is another mechanism for substitution
        in symbolic products.

        EXAMPLES::

            sage: # needs sage.combinat sage.modules
            sage: R.<x,y> = FreeAlgebra(QQ, 2)
            sage: F = Factorization([(x,3), (y, 2), (x,1)])
            sage: F(x=4)
            4^3 * y^2 * 4
            sage: F.subs({y:2})
            x^3 * 2^2 * x

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: F = Factorization([(x,3), (y, 2), (x,1)])
            sage: F(x=4)
            4 * 4^3 * y^2
            sage: F.subs({y:x})
            x * x^2 * x^3
            sage: F(x=y+x)
            (x + y) * y^2 * (x + y)^3

        TESTS::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: F = Factorization([(x-2,3), (y+3, 2)])
            sage: F(x=2)
            0

            sage: QQt = QQ['t'].fraction_field()
            sage: t = QQt.gen()
            sage: R.<x> = PolynomialRing(QQt, 1)
            sage: F = Factorization([(x,3), (x+t, 2)], unit=QQt.gen())
            sage: F(t=0)
            0

            sage: # needs sage.libs.pari sage.modules
            sage: R.<x> = LaurentPolynomialRing(QQ, 1)
            sage: F = ((x+2)/x**3).factor()
            sage: F(x=4)
            1/64 * 6
        """
        unit = self.__unit.subs(*args, **kwds)
        if unit == 0:
            return self.universe().zero()
        data = [(p.subs(*args, **kwds), e) for p, e in self.__x]
        if any(p == 0 for p, _ in data):
            return self.universe().zero()
        return Factorization(data, unit=unit, simplify=False)

    subs = __call__

    def value(self):
        """
        Return the product of the factors in the factorization, multiplied out.

        EXAMPLES::

            sage: F = factor(-2006); F
            -1 * 2 * 17 * 59
            sage: F.value()
            -2006

            sage: R.<x,y> = FreeAlgebra(ZZ, 2)                                          # needs sage.combinat sage.modules
            sage: F = Factorization([(x,3), (y, 2), (x,1)]); F                          # needs sage.combinat sage.modules
            x^3 * y^2 * x
            sage: F.value()                                                             # needs sage.combinat sage.modules
            x^3*y^2*x
        """
        from sage.misc.misc_c import prod
        return prod([p**e for p, e in self.__x], self.__unit)

    # Two aliases for ``value(self)``.
    expand = value
    prod = value

    def gcd(self, other):
        r"""
        Return the gcd of two factorizations.

        If the two factorizations have different universes, this
        method will attempt to find a common universe for the
        gcd.  A :exc:`TypeError` is raised if this is impossible.

        EXAMPLES::

            sage: factor(-30).gcd(factor(-160))
            2 * 5
            sage: factor(gcd(-30,160))
            2 * 5

            sage: R.<x> = ZZ[]
            sage: (factor(-20).gcd(factor(5*x+10))).universe()                          # needs sage.libs.pari
            Univariate Polynomial Ring in x over Integer Ring
        """
        if not isinstance(other, Factorization):
            raise NotImplementedError("can't take gcd of factorization and non-factorization")

        if len(self) and len(other):
            try:
                # first get the two factorizations to have the same
                # universe
                U = Sequence([self[0][0], other[0][0]]).universe()
                self = self.base_change(U)
                other = other.base_change(U)
            except TypeError:
                raise TypeError("Cannot take the gcd of %s and %s because they cannot be coerced into a common universe" % (self, other))

        if self.is_commutative() and other.is_commutative():
            d1 = dict(self)
            d2 = dict(other)
            s = {}
            for a in set(d1).intersection(set(d2)):
                s[a] = min(d1[a], d2[a])
            return Factorization(list(s.items()))
        else:
            raise NotImplementedError("gcd is not implemented for non-commutative factorizations")

    def lcm(self, other):
        r"""
        Return the lcm of two factorizations.

        If the two factorizations have different universes, this
        method will attempt to find a common universe for the
        lcm.  A :exc:`TypeError` is raised if this is impossible.

        EXAMPLES::

            sage: factor(-10).lcm(factor(-16))
            2^4 * 5
            sage: factor(lcm(-10,16))
            2^4 * 5

            sage: R.<x> = ZZ[]
            sage: (factor(-20).lcm(factor(5*x + 10))).universe()                        # needs sage.libs.pari
            Univariate Polynomial Ring in x over Integer Ring
        """
        if not isinstance(other, Factorization):
            raise NotImplementedError("can't take lcm of factorization and non-factorization")

        if len(self) and len(other):
            try:
                # first get the two factorizations to have the same
                # universe
                U = Sequence([self[0][0], other[0][0]]).universe()
                self = self.base_change(U)
                other = other.base_change(U)
            except TypeError:
                raise TypeError("Cannot take the lcm of %s and %s because they cannot be coerced into a common universe" % (self, other))

        if self.is_commutative() and other.is_commutative():
            d1 = dict(self)
            d2 = dict(other)
            s = {}
            for a in set(d1).union(set(d2)):
                s[a] = max(d1.get(a, 0), d2.get(a, 0))
            return Factorization(list(s.items()))
        else:
            raise NotImplementedError("lcm is not implemented for non-commutative factorizations")

    def is_integral(self) -> bool:
        r"""
        Return whether all exponents of this Factorization are nonnegative.

        EXAMPLES::

            sage: F = factor(-10); F
            -1 * 2 * 5
            sage: F.is_integral()
            True

            sage: F = factor(-10) / factor(16); F
            -1 * 2^-3 * 5
            sage: F.is_integral()
            False
        """
        return all(e >= 0 for p, e in self.__x)

    def radical(self):
        """
        Return the factorization of the radical of the value of ``self``.

        First, check that all exponents in the factorization are
        positive, raise :exc:`ValueError` otherwise.  If all exponents are
        positive, return ``self`` with all exponents set to 1 and with the
        unit set to 1.

        EXAMPLES::

            sage: F = factor(-100); F
            -1 * 2^2 * 5^2
            sage: F.radical()
            2 * 5
            sage: factor(1/2).radical()
            Traceback (most recent call last):
            ...
            ValueError: all exponents in the factorization must be positive
        """
        if not all(e > 0 for _, e in self.__x):
            raise ValueError("all exponents in the factorization must be positive")
        return Factorization([(p, 1) for p, _ in self.__x], unit=self.unit().parent()(1),
                             cr=self.__cr, sort=False, simplify=False)

    def radical_value(self):
        """
        Return the product of the prime factors in ``self``.

        First, check that all exponents in the factorization are
        positive, raise :exc:`ValueError` otherwise.  If all exponents are
        positive, return the product of the prime factors in ``self``.
        This should be functionally equivalent to
        ``self.radical().value()``.

        EXAMPLES::

            sage: F = factor(-100); F
            -1 * 2^2 * 5^2
            sage: F.radical_value()
            10
            sage: factor(1/2).radical_value()
            Traceback (most recent call last):
            ...
            ValueError: all exponents in the factorization must be positive
        """
        if not all(e > 0 for _, e in self.__x):
            raise ValueError("all exponents in the factorization must be positive")
        from sage.misc.misc_c import prod
        return prod([p for p, _ in self.__x])
