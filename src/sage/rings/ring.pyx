"""
Rings

This module provides the abstract base class :class:`Ring` from which
all rings in Sage (used to) derive, as well as a selection of more
specific base classes.

.. WARNING::

    Those classes, except maybe for the lowest ones like :class:`Ring`,
    :class:`CommutativeRing`, :class:`Algebra` and :class:`CommutativeAlgebra`,
    are being progressively deprecated in favor of the corresponding
    categories. which are more flexible, in particular with respect to multiple
    inheritance.

The class inheritance hierarchy is:

- :class:`Ring`

  - :class:`Algebra`
  - :class:`CommutativeRing`

    - :class:`NoetherianRing`
    - :class:`CommutativeAlgebra`
    - :class:`IntegralDomain`

      - :class:`DedekindDomain`
      - :class:`PrincipalIdealDomain`

Subclasses of :class:`PrincipalIdealDomain` are

- :class:`EuclideanDomain`
- :class:`Field`

  - :class:`~sage.rings.finite_rings.finite_field_base.FiniteField`

Some aspects of this structure may seem strange, but this is an unfortunate
consequence of the fact that Cython classes do not support multiple
inheritance. Hence, for instance, :class:`Field` cannot be a subclass of both
:class:`NoetherianRing` and :class:`PrincipalIdealDomain`, although all fields
are Noetherian PIDs.

(A distinct but equally awkward issue is that sometimes we may not know *in
advance* whether or not a ring belongs in one of these classes; e.g. some
orders in number fields are Dedekind domains, but others are not, and we still
want to offer a unified interface, so orders are never instances of the
:class:`DedekindDomain` class.)

AUTHORS:

- David Harvey (2006-10-16): changed :class:`CommutativeAlgebra` to derive from
  :class:`CommutativeRing` instead of from :class:`Algebra`.
- David Loeffler (2009-07-09): documentation fixes, added to reference manual.
- Simon King (2011-03-29): Proper use of the category framework for rings.
- Simon King (2011-05-20): Modify multiplication and _ideal_class_ to support
  ideals of non-commutative rings.

"""

#*****************************************************************************
#       Copyright (C) 2005, 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method

from sage.structure.coerce cimport coercion_model
from sage.structure.parent cimport Parent
from sage.structure.category_object import check_default_category
from sage.structure.sequence import Sequence
from sage.misc.prandom import randint
from sage.categories.rings import Rings
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.integral_domains import IntegralDomains
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.categories.euclidean_domains import EuclideanDomains

_Rings = Rings()
_CommutativeRings = CommutativeRings()

cdef class Ring(ParentWithGens):
    """
    Generic ring class.

    TESTS:

    This is to test against the bug fixed in :trac:`9138`::

        sage: R.<x> = QQ[]
        sage: R.sum([x,x])
        2*x
        sage: R.<x,y> = ZZ[]
        sage: R.sum([x,y])
        x + y
        sage: TestSuite(QQ['x']).run(verbose=True)
        running ._test_additive_associativity() . . . pass
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_characteristic() . . . pass
        running ._test_construction() . . . pass
        running ._test_distributivity() . . . pass
        running ._test_divides() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_nonzero_equal() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_euclidean_degree() . . . pass
        running ._test_fraction_field() . . . pass
        running ._test_gcd_vs_xgcd() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_quo_rem() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_zero() . . . pass
        running ._test_zero_divisors() . . . pass
        sage: TestSuite(QQ['x','y']).run()
        sage: TestSuite(ZZ['x','y']).run()
        sage: TestSuite(ZZ['x','y']['t']).run()

    Test against another bug fixed in :trac:`9944`::

        sage: QQ['x'].category()
        Join of Category of euclidean domains and Category of commutative algebras over
         (number fields and quotient fields and metric spaces) and Category of infinite sets
        sage: QQ['x','y'].category()
        Join of Category of unique factorization domains and Category of commutative algebras over
         (number fields and quotient fields and metric spaces) and Category of infinite sets
        sage: PolynomialRing(MatrixSpace(QQ,2),'x').category()
        Category of infinite algebras over (finite dimensional algebras with basis over
         (number fields and quotient fields and metric spaces) and infinite sets)
        sage: PolynomialRing(SteenrodAlgebra(2),'x').category()
        Category of infinite algebras over (super hopf algebras with basis
         over Finite Field of size 2 and supercocommutative super coalgebras
         over Finite Field of size 2)

    TESTS::

        sage: Zp(7)._repr_option('element_is_atomic')
        False
        sage: QQ._repr_option('element_is_atomic')
        True
        sage: CDF._repr_option('element_is_atomic')
        False

    Check that categories correctly implement `is_finite` and `cardinality`::

        sage: QQ.is_finite()
        False
        sage: GF(2^10,'a').is_finite()
        True
        sage: R.<x> = GF(7)[]
        sage: R.is_finite()
        False
        sage: S.<y> = R.quo(x^2+1)
        sage: S.is_finite()
        True

        sage: Integers(7).cardinality()
        7
        sage: QQ.cardinality()
        +Infinity
     """
    def __init__(self, base, names=None, normalize=True, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: ZZ
            Integer Ring
            sage: R.<x,y> = QQ[]
            sage: R
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        # Unfortunately, ParentWithGens inherits from sage.structure.parent_old.Parent.
        # Its __init__ method does *not* call Parent.__init__, since this would somehow
        # yield an infinite recursion. But when we call it from here, it works.
        # This is done in order to ensure that __init_extra__ is called.
        #
        # ParentWithGens.__init__(self, base, names=names, normalize=normalize)
        #
        # This is a low-level class. For performance, we trust that the category
        # is fine, if it is provided. If it isn't, we use the category of rings.
        if category is None:
            category=_Rings
        Parent.__init__(self, base=base, names=names, normalize=normalize,
                        category=category)

    def __iter__(self):
        r"""
        Return an iterator through the elements of ``self``.
        Not implemented in general.

        EXAMPLES::

            sage: sage.rings.ring.Ring.__iter__(ZZ)
            Traceback (most recent call last):
            ...
            NotImplementedError: object does not support iteration
        """
        raise NotImplementedError("object does not support iteration")

    def __len__(self):
        r"""
        Return the cardinality of this ring if it is finite, else raise
        a ``NotImplementedError``.

        EXAMPLES::

            sage: len(Integers(24))
            24
            sage: len(RR)
            Traceback (most recent call last):
            ...
            NotImplementedError: len() of an infinite set
        """
        if self.is_finite():
            return self.cardinality()
        raise NotImplementedError('len() of an infinite set')

    def __xor__(self, n):
        r"""
        Trap the operation ``^``.

        EXAMPLES::

            sage: eval('RR^3')
            Traceback (most recent call last):
            ...
            RuntimeError: use ** for exponentiation, not '^', which means xor in Python, and has the wrong precedence
        """
        raise RuntimeError("use ** for exponentiation, not '^', which means xor "
              "in Python, and has the wrong precedence")

    def base_extend(self, R):
        """
        EXAMPLES::

            sage: QQ.base_extend(GF(7))
            Traceback (most recent call last):
            ...
            TypeError: no base extension defined
            sage: ZZ.base_extend(GF(7))
            Finite Field of size 7
        """
        if R.has_coerce_map_from(self):
            return R
        raise TypeError('no base extension defined')

    def category(self):
        """
        Return the category to which this ring belongs.

        .. NOTE::

            This method exists because sometimes a ring is its own base ring.
            During initialisation of a ring `R`, it may be checked whether the
            base ring (hence, the ring itself) is a ring. Hence, it is
            necessary that ``R.category()`` tells that ``R`` is a ring, even
            *before* its category is properly initialised.

        EXAMPLES::

            sage: FreeAlgebra(QQ, 3, 'x').category() # todo: use a ring which is not an algebra!
            Category of algebras with basis over Rational Field

        Since a quotient of the integers is its own base ring, and during
        initialisation of a ring it is tested whether the base ring belongs
        to the category of rings, the following is an indirect test that the
        ``category()`` method of rings returns the category of rings
        even before the initialisation was successful::

            sage: I = Integers(15)
            sage: I.base_ring() is I
            True
            sage: I.category()
            Join of Category of finite commutative rings
                and Category of subquotients of monoids
                and Category of quotients of semigroups
                and Category of finite enumerated sets
        """
        # Defining a category method is deprecated for parents.
        # For rings, however, it is strictly needed that self.category()
        # returns (a sub-category of) the category of rings before
        # initialisation has finished.
        return self._category or _Rings

    def ideal_monoid(self):
        """
        Return the monoid of ideals of this ring.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(ZZ, 3)
            sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F
            sage: Q = F.quotient(I)
            sage: Q.ideal_monoid()
            Monoid of ideals of Quotient of Free Algebra on 3 generators (x, y, z) over Integer Ring by the ideal (x*y + y*z, x^2 + x*y - y*x - y^2)
            sage: F.<x,y,z> = FreeAlgebra(ZZ, implementation='letterplace')
            sage: I = F*[x*y+y*z,x^2+x*y-y*x-y^2]*F
            sage: Q = F.quo(I)
            sage: Q.ideal_monoid()
            Monoid of ideals of Quotient of Free Associative Unital Algebra on 3 generators (x, y, z) over Integer Ring by the ideal (x*y + y*z, x*x + x*y - y*x - y*y)

        """
        if self.__ideal_monoid is not None:
            return self.__ideal_monoid
        else:
            from sage.rings.noncommutative_ideals import IdealMonoid_nc
            M = IdealMonoid_nc(self)
            self.__ideal_monoid = M
            return M

    def ideal(self, *args, **kwds):
        """
        Return the ideal defined by ``x``, i.e., generated by ``x``.

        INPUT:

        - ``*x`` -- list or tuple of generators (or several input arguments)

        - ``coerce`` -- bool (default: ``True``); this must be a keyword
          argument. Only set it to ``False`` if you are certain that each
          generator is already in the ring.

        - ``ideal_class`` -- callable (default: ``self._ideal_class_()``);
          this must be a keyword argument. A constructor for ideals, taking
          the ring as the first argument and then the generators.
          Usually a subclass of :class:`~sage.rings.ideal.Ideal_generic` or
          :class:`~sage.rings.noncommutative_ideals.Ideal_nc`.

        - Further named arguments (such as ``side`` in the case of
          non-commutative rings) are forwarded to the ideal class.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: R.ideal(x,y)
            Ideal (x, y) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: R.ideal(x+y^2)
            Ideal (y^2 + x) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: R.ideal( [x^3,y^3+x^3] )
            Ideal (x^3, x^3 + y^3) of Multivariate Polynomial Ring in x, y over Rational Field

        Here is an example over a non-commutative ring::

            sage: A = SteenrodAlgebra(2)
            sage: A.ideal(A.1,A.2^2)
            Twosided Ideal (Sq(2), Sq(2,2)) of mod 2 Steenrod algebra, milnor basis
            sage: A.ideal(A.1,A.2^2,side='left')
            Left Ideal (Sq(2), Sq(2,2)) of mod 2 Steenrod algebra, milnor basis

        TESTS:

        Make sure that :trac:`11139` is fixed::

            sage: R.<x> = QQ[]
            sage: R.ideal([])
            Principal ideal (0) of Univariate Polynomial Ring in x over Rational Field
            sage: R.ideal(())
            Principal ideal (0) of Univariate Polynomial Ring in x over Rational Field
            sage: R.ideal()
            Principal ideal (0) of Univariate Polynomial Ring in x over Rational Field
        """
        if 'coerce' in kwds:
            coerce = kwds['coerce']
            del kwds['coerce']
        else:
            coerce = True

        from sage.rings.ideal import Ideal_generic
        from sage.structure.parent import is_Parent
        gens = args
        while isinstance(gens, (list, tuple)) and len(gens) == 1:
            first = gens[0]
            if isinstance(first, Ideal_generic):
                R = first.ring()
                m = self.convert_map_from(R)
                if m is not None:
                    gens = [m(g) for g in first.gens()]
                    coerce = False
                else:
                    m = R.convert_map_from(self)
                    if m is not None:
                        raise NotImplementedError
                    else:
                        raise TypeError
                break
            elif isinstance(first, (list, tuple)):
                gens = first
            elif is_Parent(first) and self.has_coerce_map_from(first):
                gens = first.gens() # we have a ring as argument
            else:
                break

        if len(gens) == 0:
            gens = [self.zero()]

        if coerce:
            gens = [self(g) for g in gens]
        if isinstance(self, PrincipalIdealDomain):
            # Use GCD algorithm to obtain a principal ideal
            g = gens[0]
            if len(gens) == 1:
                try:
                    g = g.gcd(g) # note: we set g = gcd(g, g) to "canonicalize" the generator: make polynomials monic, etc.
                except (AttributeError, NotImplementedError):
                    pass
            else:
                for h in gens[1:]:
                    g = g.gcd(h)
            gens = [g]
        if 'ideal_class' in kwds:
            C = kwds['ideal_class']
            del kwds['ideal_class']
        else:
            C = self._ideal_class_(len(gens))
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        return C(self, gens, **kwds)

    def __mul__(self, x):
        """
        Return the ideal ``x*R`` generated by ``x``, where ``x`` is either an
        element or tuple or list of elements.

        EXAMPLES::

            sage: R.<x,y,z> = GF(7)[]
            sage: (x+y)*R
            Ideal (x + y) of Multivariate Polynomial Ring in x, y, z over Finite Field of size 7
            sage: (x+y,z+y^3)*R
            Ideal (x + y, y^3 + z) of Multivariate Polynomial Ring in x, y, z over Finite Field of size 7

        The following was implemented in :trac:`7797`::

            sage: A = SteenrodAlgebra(2)
            sage: A*[A.1+A.2,A.1^2]
            Left Ideal (Sq(2) + Sq(4), Sq(1,1)) of mod 2 Steenrod algebra, milnor basis
            sage: [A.1+A.2,A.1^2]*A
            Right Ideal (Sq(2) + Sq(4), Sq(1,1)) of mod 2 Steenrod algebra, milnor basis
            sage: A*[A.1+A.2,A.1^2]*A
            Twosided Ideal (Sq(2) + Sq(4), Sq(1,1)) of mod 2 Steenrod algebra, milnor basis

        """
        if isinstance(self, Ring):
            if self.is_commutative():
                return self.ideal(x)
            try:
                side = x.side()
            except AttributeError:
                return self.ideal(x, side='left')
            # presumably x is an ideal...
            try:
                x = x.gens()
            except (AttributeError, NotImplementedError):
                pass # ... not an ideal
            if side in ['left','twosided']:
                return self.ideal(x,side=side)
            elif side=='right':
                return self.ideal(x,side='twosided')
            else: # duck typing failed
                raise TypeError("Don't know how to transform %s into an ideal of %s" % (x, self))
        else: # the sides are switched because this is a Cython / extension class
            if x.is_commutative():
                return x.ideal(self)
            try:
                side = self.side()
            except AttributeError:
                return x.ideal(self, side='right')
            # presumably self is an ideal...
            try:
                self = self.gens()
            except (AttributeError, NotImplementedError):
                pass # ... not an ideal
            if side in ['right','twosided']:
                return x.ideal(self,side='twosided')
            elif side=='left':
                return x.ideal(self,side='twosided')
            else:
                raise TypeError("Don't know how to transform %s into an ideal of %s" % (self, x))

    def _ideal_class_(self, n=0):
        r"""
        Return a callable object that can be used to create ideals in this
        ring. For generic rings, this returns the factory function
        :func:`sage.rings.ideal.Ideal`, which does its best to be clever about
        what is required.

        This class can depend on `n`, the number of generators of the ideal.
        The default input of `n=0` indicates an unspecified number of generators,
        in which case a class that works for any number of generators is returned.

        EXAMPLES::

            sage: ZZ._ideal_class_()
            <class 'sage.rings.ideal.Ideal_pid'>
            sage: RR._ideal_class_()
            <class 'sage.rings.ideal.Ideal_pid'>
            sage: R.<x,y> = GF(5)[]
            sage: R._ideal_class_(1)
            <class 'sage.rings.polynomial.multi_polynomial_ideal.MPolynomialIdeal'>
            sage: S = R.quo(x^3-y^2)
            sage: S._ideal_class_(1)
            <class 'sage.rings.quotient_ring.QuotientRingIdeal_principal'>
            sage: S._ideal_class_(2)
            <class 'sage.rings.quotient_ring.QuotientRingIdeal_generic'>
            sage: T.<z> = S[]
            sage: T._ideal_class_(5)
            <class 'sage.rings.ideal.Ideal_generic'>
            sage: T._ideal_class_(1)
            <class 'sage.rings.ideal.Ideal_principal'>

        Since :trac:`7797`, non-commutative rings have ideals as well::

            sage: A = SteenrodAlgebra(2)
            sage: A._ideal_class_()
            <class 'sage.rings.noncommutative_ideals.Ideal_nc'>

        """
        # One might need more than just n, but I can't think of an example.
        from sage.rings.noncommutative_ideals import Ideal_nc
        try:
            if not self.is_commutative():
                return Ideal_nc
        except (NotImplementedError, AttributeError):
            return Ideal_nc
        from sage.rings.ideal import Ideal_generic, Ideal_principal
        if n == 1:
            return Ideal_principal
        else:
            return Ideal_generic

    def principal_ideal(self, gen, coerce=True):
        """
        Return the principal ideal generated by gen.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: R.principal_ideal(x+2*y)
            Ideal (x + 2*y) of Multivariate Polynomial Ring in x, y over Integer Ring
        """
        C = self._ideal_class_(1)
        if coerce:
            gen = self(gen)
        return C(self, [gen])

    def unit_ideal(self):
        """
        Return the unit ideal of this ring.

        EXAMPLES::

            sage: Zp(7).unit_ideal()
            Principal ideal (1 + O(7^20)) of 7-adic Ring with capped relative precision 20
        """
        if self._unit_ideal is None:
            I = Ring.ideal(self, [self(1)], coerce=False)
            self._unit_ideal = I
            return I
        return self._unit_ideal

    def zero_ideal(self):
        """
        Return the zero ideal of this ring (cached).

        EXAMPLES::

            sage: ZZ.zero_ideal()
            Principal ideal (0) of Integer Ring
            sage: QQ.zero_ideal()
            Principal ideal (0) of Rational Field
            sage: QQ['x'].zero_ideal()
            Principal ideal (0) of Univariate Polynomial Ring in x over Rational Field

        The result is cached::

            sage: ZZ.zero_ideal() is ZZ.zero_ideal()
            True

        TESTS:

        Make sure that :trac:`13644` is fixed::

            sage: K = Qp(3)
            sage: R.<a> = K[]
            sage: L.<a> = K.extension(a^2-3)
            sage: L.ideal(a)
            Principal ideal (1 + O(a^40)) of 3-adic Eisenstein Extension Field in a defined by a^2 - 3

        """
        if self._zero_ideal is None:
            I = Ring.ideal(self, [self.zero()], coerce=False)
            self._zero_ideal = I
            return I
        return self._zero_ideal

    def zero(self):
        """
        Return the zero element of this ring (cached).

        EXAMPLES::

            sage: ZZ.zero()
            0
            sage: QQ.zero()
            0
            sage: QQ['x'].zero()
            0

        The result is cached::

            sage: ZZ.zero() is ZZ.zero()
            True
        """
        if self._zero_element is None:
            x = self(0)
            self._zero_element = x
            return x
        return self._zero_element

    def one(self):
        """
        Return the one element of this ring (cached), if it exists.

        EXAMPLES::

            sage: ZZ.one()
            1
            sage: QQ.one()
            1
            sage: QQ['x'].one()
            1

        The result is cached::

            sage: ZZ.one() is ZZ.one()
            True
        """
        if self._one_element is None:
            x = self(1)
            self._one_element = x
            return x
        return self._one_element

    def is_commutative(self):
        """
        Return ``True`` if this ring is commutative.

        EXAMPLES::

            sage: QQ.is_commutative()
            True
            sage: QQ['x,y,z'].is_commutative()
            True
            sage: Q.<i,j,k> = QuaternionAlgebra(QQ, -1,-1)
            sage: Q.is_commutative()
            False
        """
        if self.is_zero():
            return True
        raise NotImplementedError

    def is_field(self, proof = True):
        """
        Return ``True`` if this ring is a field.

        INPUT:

        - ``proof`` -- (default: ``True``) Determines what to do in unknown
          cases

        ALGORITHM:

        If the parameter ``proof`` is set to ``True``, the returned value is
        correct but the method might throw an error.  Otherwise, if it is set
        to ``False``, the method returns True if it can establish that self is
        a field and False otherwise.

        EXAMPLES::

            sage: QQ.is_field()
            True
            sage: GF(9,'a').is_field()
            True
            sage: ZZ.is_field()
            False
            sage: QQ['x'].is_field()
            False
            sage: Frac(QQ['x']).is_field()
            True

        This illustrates the use of the ``proof`` parameter::

            sage: R.<a,b> = QQ[]
            sage: S.<x,y> = R.quo((b^3))
            sage: S.is_field(proof = True)
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: S.is_field(proof = False)
            False
        """
        if self.is_zero():
            return False

        if proof:
            raise NotImplementedError("No way to prove that %s is an integral domain!" % self)
        else:
            return False

    cpdef bint is_exact(self) except -2:
        """
        Return ``True`` if elements of this ring are represented exactly, i.e.,
        there is no precision loss when doing arithmetic.

        .. NOTE::

            This defaults to ``True``, so even if it does return ``True`` you
            have no guarantee (unless the ring has properly overloaded this).

        EXAMPLES::

            sage: QQ.is_exact()    # indirect doctest
            True
            sage: ZZ.is_exact()
            True
            sage: Qp(7).is_exact()
            False
            sage: Zp(7, type='capped-abs').is_exact()
            False
        """
        return True

    def is_subring(self, other):
        """
        Return ``True`` if the canonical map from ``self`` to ``other`` is
        injective.

        Raises a ``NotImplementedError`` if not known.

        EXAMPLES::

            sage: ZZ.is_subring(QQ)
            True
            sage: ZZ.is_subring(GF(19))
            False

        TESTS::

            sage: QQ.is_subring(QQ['x'])
            True
            sage: QQ.is_subring(GF(7))
            False
            sage: QQ.is_subring(CyclotomicField(7))
            True
            sage: QQ.is_subring(ZZ)
            False

        Every ring is a subring of itself, :trac:`17287`::

            sage: QQbar.is_subring(QQbar)
            True
            sage: RR.is_subring(RR)
            True
            sage: CC.is_subring(CC)
            True
            sage: K.<a> = NumberField(x^3-x+1/10)
            sage: K.is_subring(K)
            True
            sage: R.<x> = RR[]
            sage: R.is_subring(R)
            True
        """
        if self is other:
            return True
        try:
            return self.Hom(other).natural_map().is_injective()
        except (TypeError, AttributeError):
            return False

    def is_prime_field(self):
        r"""
        Return ``True`` if this ring is one of the prime fields `\QQ` or
        `\GF{p}`.

        EXAMPLES::

            sage: QQ.is_prime_field()
            True
            sage: GF(3).is_prime_field()
            True
            sage: GF(9,'a').is_prime_field()
            False
            sage: ZZ.is_prime_field()
            False
            sage: QQ['x'].is_prime_field()
            False
            sage: Qp(19).is_prime_field()
            False
        """
        return False

    def is_integral_domain(self, proof = True):
        """
        Return ``True`` if this ring is an integral domain.

        INPUT:

        - ``proof`` -- (default: ``True``) Determines what to do in unknown
          cases

        ALGORITHM:

        If the parameter ``proof`` is set to ``True``, the returned value is
        correct but the method might throw an error.  Otherwise, if it is set
        to ``False``, the method returns ``True`` if it can establish that self
        is an integral domain and ``False`` otherwise.

        EXAMPLES::

            sage: QQ.is_integral_domain()
            True
            sage: ZZ.is_integral_domain()
            True
            sage: ZZ['x,y,z'].is_integral_domain()
            True
            sage: Integers(8).is_integral_domain()
            False
            sage: Zp(7).is_integral_domain()
            True
            sage: Qp(7).is_integral_domain()
            True
            sage: R.<a,b> = QQ[]
            sage: S.<x,y> = R.quo((b^3))
            sage: S.is_integral_domain()
            False

        This illustrates the use of the ``proof`` parameter::

            sage: R.<a,b> = ZZ[]
            sage: S.<x,y> = R.quo((b^3))
            sage: S.is_integral_domain(proof = True)
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: S.is_integral_domain(proof = False)
            False

        TESTS:

        Make sure :trac:`10481` is fixed::

            sage: var('x')
            x
            sage: R.<a> = ZZ['x'].quo(x^2)
            sage: R.fraction_field()
            Traceback (most recent call last):
            ...
            TypeError: self must be an integral domain.
            sage: R.is_integral_domain()
            False

        Forward the proof flag to ``is_field``, see :trac:`22910`::

            sage: R1.<x> = GF(5)[]
            sage: F1 = R1.quotient_ring(x^2+x+1)
            sage: R2.<x> = F1[]
            sage: F2 = R2.quotient_ring(x^2+x+1)
            sage: F2.is_integral_domain(False)
            False
        """
        if self.is_field(proof):
            return True

        if self.is_zero():
            return False

        if proof:
            raise NotImplementedError
        else:
            return False

    def is_noetherian(self):
        """
        Return ``True`` if this ring is Noetherian.

        EXAMPLES::

            sage: QQ.is_noetherian()
            True
            sage: ZZ.is_noetherian()
            True
        """
        raise NotImplementedError

    def order(self):
        """
        The number of elements of ``self``.

        EXAMPLES::

            sage: GF(19).order()
            19
            sage: QQ.order()
            +Infinity
        """
        if self.is_zero():
            return 1
        raise NotImplementedError

    def zeta(self, n=2, all=False):
        """
        Return a primitive ``n``-th root of unity in ``self`` if there
        is one, or raise a ``ValueError`` otherwise.

        INPUT:

        - ``n`` -- positive integer

        - ``all`` -- bool (default: False) - whether to return
          a list of all primitive `n`-th roots of unity. If True, raise a ``ValueError``
          if ``self`` is not an integral domain.

        OUTPUT:

        Element of ``self`` of finite order

        EXAMPLES::

            sage: QQ.zeta()
            -1
            sage: QQ.zeta(1)
            1
            sage: CyclotomicField(6).zeta(6)
            zeta6
            sage: CyclotomicField(3).zeta(3)
            zeta3
            sage: CyclotomicField(3).zeta(3).multiplicative_order()
            3
            sage: a = GF(7).zeta(); a
            3
            sage: a.multiplicative_order()
            6
            sage: a = GF(49,'z').zeta(); a
            z
            sage: a.multiplicative_order()
            48
            sage: a = GF(49,'z').zeta(2); a
            6
            sage: a.multiplicative_order()
            2
            sage: QQ.zeta(3)
            Traceback (most recent call last):
            ...
            ValueError: no n-th root of unity in rational field
            sage: Zp(7, prec=8).zeta()
            3 + 4*7 + 6*7^2 + 3*7^3 + 2*7^5 + 6*7^6 + 2*7^7 + O(7^8)

        TESTS::

            sage: from sage.rings.ring import Ring
            sage: Ring.zeta(QQ, 1)
            1
            sage: Ring.zeta(QQ, 2)
            -1
            sage: Ring.zeta(QQ, 3)
            Traceback (most recent call last):
            ...
            ValueError: no 3rd root of unity in Rational Field
            sage: IntegerModRing(8).zeta(2, all = True)
            Traceback (most recent call last):
            ...
            ValueError: ring is not an integral domain
        """
        if all and not self.is_integral_domain():
            raise ValueError("ring is not an integral domain")
        if n == 2:
            if all:
                return [self(-1)]
            else:
                return self(-1)
        elif n == 1:
            if all:
                return [self(1)]
            else:
                return self(1)
        else:
            f = self['x'].cyclotomic_polynomial(n)
            if all:
                return [-P[0] for P, e in f.factor() if P.degree() == 1]
            for P, e in f.factor():
                if P.degree() == 1:
                    return -P[0]
            from sage.rings.integer_ring import ZZ
            raise ValueError("no %s root of unity in %r" % (ZZ(n).ordinal_str(), self))

    def zeta_order(self):
        """
        Return the order of the distinguished root of unity in ``self``.

        EXAMPLES::

            sage: CyclotomicField(19).zeta_order()
            38
            sage: GF(19).zeta_order()
            18
            sage: GF(5^3,'a').zeta_order()
            124
            sage: Zp(7, prec=8).zeta_order()
            6
        """
        return self.zeta().multiplicative_order()

    def random_element(self, bound=2):
        """
        Return a random integer coerced into this ring, where the
        integer is chosen uniformly from the interval ``[-bound,bound]``.

        INPUT:

        - ``bound`` -- integer (default: 2)

        ALGORITHM:

        Uses Python's randint.

        TESTS:

        The following example returns a ``NotImplementedError`` since the
        generic ring class ``__call__`` function returns a
        ``NotImplementedError``. Note that
        ``sage.rings.ring.Ring.random_element`` performs a call in the generic
        ring class by a random integer::

            sage: R = sage.rings.ring.Ring(ZZ); R
            <sage.rings.ring.Ring object at ...>
            sage: R.random_element()
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot construct elements of <sage.rings.ring.Ring object at ...>
        """
        return self(randint(-bound,bound))

    def _random_nonzero_element(self, *args, **kwds):
        """
        Returns a random non-zero element in this ring.

        The default behaviour of this method is to repeatedly call the
        ``random_element`` method until a non-zero element is obtained.
        In this implementation, all parameters are simply pushed forward
        to the ``random_element`` method.

        INPUT:

        -  ``*args``, ``**kwds`` - Parameters that can be forwarded to
           the ``random_element`` method

        OUTPUT:

        - Random non-zero element

        EXAMPLES::

            sage: ZZ._random_nonzero_element() != 0
            True
        """
        while True:
            x = self.random_element(*args, **kwds)
            if not x.is_zero():
                return x

    def ideal_monoid(self):
        """
        Return the monoid of ideals of this ring.

        EXAMPLES::

            sage: ZZ.ideal_monoid()
            Monoid of ideals of Integer Ring
            sage: R.<x>=QQ[]; R.ideal_monoid()
            Monoid of ideals of Univariate Polynomial Ring in x over Rational Field
        """
        if self.__ideal_monoid is not None:
            return self.__ideal_monoid
        else:
            from sage.rings.ideal_monoid import IdealMonoid
            M = IdealMonoid(self)
            self.__ideal_monoid = M
            return M

    @cached_method
    def epsilon(self):
        """
        Return the precision error of elements in this ring.

        EXAMPLES::

            sage: RDF.epsilon()
            2.220446049250313e-16
            sage: ComplexField(53).epsilon()
            2.22044604925031e-16
            sage: RealField(10).epsilon()
            0.0020

        For exact rings, zero is returned::

            sage: ZZ.epsilon()
            0

        This also works over derived rings::

            sage: RR['x'].epsilon()
            2.22044604925031e-16
            sage: QQ['x'].epsilon()
            0

        For the symbolic ring, there is no reasonable answer::

            sage: SR.epsilon()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        one = self.one()
        try:
            return one.ulp()
        except AttributeError:
            pass

        try:
            eps = one.real().ulp()
        except AttributeError:
            pass
        else:
            return self(eps)

        B = self._base
        if B is not None and B is not self:
            eps = self.base_ring().epsilon()
            return self(eps)
        if self.is_exact():
            return self.zero()
        raise NotImplementedError

cdef class CommutativeRing(Ring):
    """
    Generic commutative ring.
    """
    def __init__(self, base_ring, names=None, normalize=True, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Integers(389)['x,y']
            Multivariate Polynomial Ring in x, y over Ring of integers modulo 389
        """
        try:
            if not base_ring.is_commutative():
                raise TypeError("base ring %s is no commutative ring" % base_ring)
        except AttributeError:
            raise TypeError("base ring %s is no commutative ring" % base_ring)
        # This is a low-level class. For performance, we trust that
        # the category is fine, if it is provided. If it isn't, we use
        # the category of commutative rings.
        if category is None:
            category=_CommutativeRings
        Ring.__init__(self, base_ring, names=names, normalize=normalize,
                      category=category)

    def localization(self, additional_units, names=None, normalize=True, category=None):
        """
        Return the localization of ``self`` at the given additional units.

        EXAMPLES::

            sage: R.<x, y> = GF(3)[]
            sage: R.localization((x*y, x**2+y**2))
            Multivariate Polynomial Ring in x, y over Finite Field of size 3 localized at (y, x, x^2 + y^2)
            sage: ~y in _
            True
        """
        if not self.is_integral_domain():
            raise TypeError("self must be an integral domain.")

        from sage.rings.localization import Localization
        return Localization(self, additional_units, names=names, normalize=normalize, category=category)

    def fraction_field(self):
        """
        Return the fraction field of ``self``.

        EXAMPLES::

            sage: R = Integers(389)['x,y']
            sage: Frac(R)
            Fraction Field of Multivariate Polynomial Ring in x, y over Ring of integers modulo 389
            sage: R.fraction_field()
            Fraction Field of Multivariate Polynomial Ring in x, y over Ring of integers modulo 389
        """
        try:
            if self.is_field():
                return self
        except NotImplementedError:
            pass

        if not self.is_integral_domain():
            raise TypeError("self must be an integral domain.")

        if self.__fraction_field is not None:
            return self.__fraction_field
        else:
            import sage.rings.fraction_field
            K = sage.rings.fraction_field.FractionField_generic(self)
            self.__fraction_field = K
        return self.__fraction_field

    def _pseudo_fraction_field(self):
        r"""
        This method is used by the coercion model to determine if `a / b`
        should be treated as `a * (1/b)`, for example when dividing an element
        of `\ZZ[x]` by an element of `\ZZ`.

        The default is to return the same value as ``self.fraction_field()``,
        but it may return some other domain in which division is usually
        defined (for example, ``\ZZ/n\ZZ`` for possibly composite `n`).

        EXAMPLES::

            sage: ZZ._pseudo_fraction_field()
            Rational Field
            sage: ZZ['x']._pseudo_fraction_field()
            Fraction Field of Univariate Polynomial Ring in x over Integer Ring
            sage: Integers(15)._pseudo_fraction_field()
            Ring of integers modulo 15
            sage: Integers(15).fraction_field()
            Traceback (most recent call last):
            ...
            TypeError: self must be an integral domain.
        """
        try:
            return self.fraction_field()
        except (NotImplementedError,TypeError):
            return coercion_model.division_parent(self)

    def is_commutative(self):
        """
        Return ``True``, since this ring is commutative.

        EXAMPLES::

            sage: QQ.is_commutative()
            True
            sage: ZpCA(7).is_commutative()
            True
            sage: A = QuaternionAlgebra(QQ, -1, -3, names=('i','j','k')); A
            Quaternion Algebra (-1, -3) with base ring Rational Field
            sage: A.is_commutative()
            False
        """
        return True

    def krull_dimension(self):
        """
        Return the Krull dimension of this commutative ring.

        The Krull dimension is the length of the longest ascending chain
        of prime ideals.

        TESTS:

        ``krull_dimension`` is not implemented for generic commutative
        rings. Fields and PIDs, with Krull dimension equal to 0 and 1,
        respectively, have naive implementations of ``krull_dimension``.
        Orders in number fields also have Krull dimension 1::

            sage: R = CommutativeRing(ZZ)
            sage: R.krull_dimension()
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: QQ.krull_dimension()
            0
            sage: ZZ.krull_dimension()
            1
            sage: type(R); type(QQ); type(ZZ)
            <class 'sage.rings.ring.CommutativeRing'>
            <class 'sage.rings.rational_field.RationalField_with_category'>
            <class 'sage.rings.integer_ring.IntegerRing_class'>

        All orders in number fields have Krull dimension 1, including
        non-maximal orders::

            sage: K.<i> = QuadraticField(-1)
            sage: R = K.maximal_order(); R
            Gaussian Integers in Number Field in i with defining polynomial x^2 + 1 with i = 1*I
            sage: R.krull_dimension()
            1
            sage: R = K.order(2*i); R
            Order in Number Field in i with defining polynomial x^2 + 1 with i = 1*I
            sage: R.is_maximal()
            False
            sage: R.krull_dimension()
            1
        """
        raise NotImplementedError

    def ideal_monoid(self):
        """
        Return the monoid of ideals of this ring.

        EXAMPLES::

            sage: ZZ.ideal_monoid()
            Monoid of ideals of Integer Ring
            sage: R.<x>=QQ[]; R.ideal_monoid()
            Monoid of ideals of Univariate Polynomial Ring in x over Rational Field
        """
        if self.__ideal_monoid is not None:
            return self.__ideal_monoid
        else:
            from sage.rings.ideal_monoid import IdealMonoid
            M = IdealMonoid(self)
            self.__ideal_monoid = M
            return M

    def extension(self, poly, name=None, names=None, **kwds):
        """
        Algebraically extends self by taking the quotient ``self[x] / (f(x))``.

        INPUT:

        - ``poly`` -- A polynomial whose coefficients are coercible into
          ``self``

        - ``name`` -- (optional) name for the root of `f`

        .. NOTE::

            Using this method on an algebraically complete field does *not*
            return this field; the construction ``self[x] / (f(x))`` is done
            anyway.

        EXAMPLES::

            sage: R = QQ['x']
            sage: y = polygen(R)
            sage: R.extension(y^2 - 5, 'a')
            Univariate Quotient Polynomial Ring in a over Univariate Polynomial Ring in x over Rational Field with modulus a^2 - 5

        ::

            sage: P.<x> = PolynomialRing(GF(5))
            sage: F.<a> = GF(5).extension(x^2 - 2)
            sage: P.<t> = F[]
            sage: R.<b> = F.extension(t^2 - a); R
            Univariate Quotient Polynomial Ring in b over Finite Field in a of size 5^2 with modulus b^2 + 4*a
        """
        from sage.rings.polynomial.polynomial_element import Polynomial
        if not isinstance(poly, Polynomial):
            try:
                poly = poly.polynomial(self)
            except (AttributeError, TypeError):
                raise TypeError("polynomial (=%s) must be a polynomial." % repr(poly))
        if not names is None:
            name = names
        if isinstance(name, tuple):
            name = name[0]
        if name is None:
            name = str(poly.parent().gen(0))
        for key, val in kwds.items():
            if key not in ['structure', 'implementation', 'prec', 'embedding', 'latex_name', 'latex_names']:
                raise TypeError("extension() got an unexpected keyword argument '%s'" % key)
            if not (val is None or isinstance(val, list) and all(c is None for c in val)):
                raise NotImplementedError("ring extension with prescribed %s is not implemented" % key)
        R = self[name]
        I = R.ideal(R(poly.list()))
        return R.quotient(I, name)

    def frobenius_endomorphism(self, n=1):
        """
        INPUT:

        - ``n`` -- a nonnegative integer (default: 1)

        OUTPUT:

        The `n`-th power of the absolute arithmetic Frobenius
        endomorphism on this finite field.

        EXAMPLES::

            sage: K.<u> = PowerSeriesRing(GF(5))
            sage: Frob = K.frobenius_endomorphism(); Frob
            Frobenius endomorphism x |--> x^5 of Power Series Ring in u over Finite Field of size 5
            sage: Frob(u)
            u^5

        We can specify a power::

            sage: f = K.frobenius_endomorphism(2); f
            Frobenius endomorphism x |--> x^(5^2) of Power Series Ring in u over Finite Field of size 5
            sage: f(1+u)
            1 + u^25
        """
        from .morphism import FrobeniusEndomorphism_generic
        return FrobeniusEndomorphism_generic(self, n)

    def derivation_module(self, codomain=None, twist=None):
        r"""
        Returns the module of derivations over this ring.

        INPUT:

        - ``codomain`` -- an algebra over this ring or a ring homomorphism
          whose domain is this ring or ``None`` (default: ``None``); if it
          is a morphism, the codomain of derivations will be the codomain
          of the morphism viewed as an algebra over ``self`` through the
          given morphism; if ``None``, the codomain will be this ring

        - ``twist`` -- a morphism from this ring to ``codomain``
          or ``None`` (default: ``None``); if ``None``, the coercion
          map from this ring to ``codomain`` will be used

        .. NOTE::

            A twisted derivation with respect to `\theta` (or a
            `\theta`-derivation for short) is an additive map `d`
            satisfying the following axiom for all `x, y` in the domain:

            .. MATH::

                d(xy) = \theta(x) d(y) + d(x) y.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: M = R.derivation_module(); M
            Module of derivations over Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: M.gens()
            (d/dx, d/dy, d/dz)

        We can specify a different codomain::

            sage: K = R.fraction_field()
            sage: M = R.derivation_module(K); M
            Module of derivations from Multivariate Polynomial Ring in x, y, z over
             Rational Field to Fraction Field of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: M.gen() / x
            1/x*d/dx

        Here is an example with a non-canonical defining morphism::

            sage: ev = R.hom([QQ(0), QQ(1), QQ(2)])
            sage: ev
            Ring morphism:
              From: Multivariate Polynomial Ring in x, y, z over Rational Field
              To:   Rational Field
              Defn: x |--> 0
                    y |--> 1
                    z |--> 2
            sage: M = R.derivation_module(ev)
            sage: M
            Module of derivations from Multivariate Polynomial Ring in x, y, z over Rational Field to Rational Field

        Elements in `M` acts as derivations at `(0,1,2)`::

            sage: Dx = M.gen(0); Dx
            d/dx
            sage: Dy = M.gen(1); Dy
            d/dy
            sage: Dz = M.gen(2); Dz
            d/dz
            sage: f = x^2 + y^2 + z^2
            sage: Dx(f)  # = 2*x evaluated at (0,1,2)
            0
            sage: Dy(f)  # = 2*y evaluated at (0,1,2)
            2
            sage: Dz(f)  # = 2*z evaluated at (0,1,2)
            4

        An example with a twisting homomorphism::

            sage: theta = R.hom([x^2, y^2, z^2])
            sage: M = R.derivation_module(twist=theta); M
            Module of twisted derivations over Multivariate Polynomial Ring in x, y, z
             over Rational Field (twisting morphism: x |--> x^2, y |--> y^2, z |--> z^2)

        .. SEEALSO::

            :meth:`derivation`

        """
        from sage.rings.derivation import RingDerivationModule
        if codomain is None:
            codomain = self
        return RingDerivationModule(self, codomain, twist)

    def derivation(self, arg=None, twist=None):
        r"""
        Return the twisted or untwisted derivation over this ring
        specified by ``arg``.

        .. NOTE::

            A twisted derivation with respect to `\theta` (or a
            `\theta`-derivation for short) is an additive map `d`
            satisfying the following axiom for all `x, y` in the domain:

            .. MATH::

                d(xy) = \theta(x) d(y) + d(x) y.

        INPUT:

        - ``arg`` -- (optional) a generator or a list of coefficients
          that defines the derivation

        - ``twist`` -- (optional) the twisting homomorphism

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: R.derivation()
            d/dx

        In that case, ``arg`` could be a generator::

            sage: R.derivation(y)
            d/dy

        or a list of coefficients::

            sage: R.derivation([1,2,3])
            d/dx + 2*d/dy + 3*d/dz

        It is not possible to define derivations with respect to a
        polynomial which is not a variable::

            sage: R.derivation(x^2)
            Traceback (most recent call last):
            ...
            ValueError: unable to create the derivation

        Here is an example with twisted derivations::

            sage: R.<x,y,z> = QQ[]
            sage: theta = R.hom([x^2, y^2, z^2])
            sage: f = R.derivation(twist=theta); f
            0
            sage: f.parent()
            Module of twisted derivations over Multivariate Polynomial Ring in x, y, z
             over Rational Field (twisting morphism: x |--> x^2, y |--> y^2, z |--> z^2)

        Specifying a scalar, the returned twisted derivation is the
        corresponding multiple of `\theta - id`::

            sage: R.derivation(1, twist=theta)
            [x |--> x^2, y |--> y^2, z |--> z^2] - id
            sage: R.derivation(x, twist=theta)
            x*([x |--> x^2, y |--> y^2, z |--> z^2] - id)

        """
        if isinstance(arg, (list, tuple)):
            codomain = Sequence([self(0)] + list(arg)).universe()
        else:
            codomain = self
        return self.derivation_module(codomain, twist=twist)(arg)


cdef class IntegralDomain(CommutativeRing):
    """
    Generic integral domain class.

    This class is deprecated. Please use the
    :class:`sage.categories.integral_domains.IntegralDomains`
    category instead.
    """
    _default_category = IntegralDomains()

    def __init__(self, base_ring, names=None, normalize=True, category=None):
        """
        Initialize ``self``.

        INPUT:

         - ``category`` (default: ``None``) -- a category, or ``None``

        This method is used by all the abstract subclasses of
        :class:`IntegralDomain`, like :class:`NoetherianRing`,
        :class:`PrincipalIdealDomain`, :class:`DedekindDomain`,
        :class:`EuclideanDomain`, :class:`Field`, ... in order to
        avoid cascade calls Field.__init__ ->
        PrincipalIdealDomain.__init__ -> IntegralDomain.__init__ ->
        ...

        EXAMPLES::

            sage: F = IntegralDomain(QQ)
            sage: F.category()
            Category of integral domains

            sage: F = PrincipalIdealDomain(QQ)
            sage: F.category()
            Category of principal ideal domains

            sage: F = EuclideanDomain(QQ)
            sage: F.category()
            Category of euclidean domains

            sage: F = Field(QQ)
            sage: F.category()
            Category of fields

        If a category is specified, then the category is set to the
        join of that category with the default category::

            sage: F = PrincipalIdealDomain(QQ, category=EnumeratedSets())

        The default value for the category is specified by the class
        attribute ``default_category``::

            sage: IntegralDomain._default_category
            Category of integral domains

            sage: PrincipalIdealDomain._default_category
            Category of principal ideal domains

            sage: EuclideanDomain._default_category
            Category of euclidean domains

            sage: Field._default_category
            Category of fields

        """
        category = check_default_category(self._default_category, category)
        CommutativeRing.__init__(self, base_ring, names=names, normalize=normalize,
                                 category=category)

    def is_integral_domain(self, proof = True):
        """
        Return ``True``, since this ring is an integral domain.

        (This is a naive implementation for objects with type
        ``IntegralDomain``)

        EXAMPLES::

            sage: ZZ.is_integral_domain()
            True
            sage: QQ.is_integral_domain()
            True
            sage: ZZ['x'].is_integral_domain()
            True
            sage: R = ZZ.quotient(ZZ.ideal(10)); R.is_integral_domain()
            False
        """
        return True

    def is_integrally_closed(self):
        r"""
        Return ``True`` if this ring is integrally closed in its field of
        fractions; otherwise return ``False``.

        When no algorithm is implemented for this, then this
        function raises a ``NotImplementedError``.

        Note that ``is_integrally_closed`` has a naive implementation
        in fields. For every field `F`, `F` is its own field of fractions,
        hence every element of `F` is integral over `F`.

        EXAMPLES::

            sage: ZZ.is_integrally_closed()
            True
            sage: QQ.is_integrally_closed()
            True
            sage: QQbar.is_integrally_closed()
            True
            sage: GF(5).is_integrally_closed()
            True
            sage: Z5 = Integers(5); Z5
            Ring of integers modulo 5
            sage: Z5.is_integrally_closed()
            Traceback (most recent call last):
            ...
            AttributeError: 'IntegerModRing_generic_with_category' object has no attribute 'is_integrally_closed'
        """
        raise NotImplementedError

    def is_field(self, proof = True):
        r"""
        Return ``True`` if this ring is a field.

        EXAMPLES::

            sage: GF(7).is_field()
            True

        The following examples have their own ``is_field`` implementations::

            sage: ZZ.is_field(); QQ.is_field()
            False
            True
            sage: R.<x> = PolynomialRing(QQ); R.is_field()
            False
        """
        if self.is_finite():
            return True
        if proof:
            raise NotImplementedError("unable to determine whether or not is a field.")
        else:
            return False

cdef class NoetherianRing(CommutativeRing):
    """
    Generic Noetherian ring class.

    A Noetherian ring is a commutative ring in which every ideal is
    finitely generated.

    This class is deprecated, and not actually used anywhere in the
    Sage code base.  If you think you need it, please create a
    category :class:`NoetherianRings`, move the code of this class
    there, and use it instead.
    """
    def is_noetherian(self):
        """
        Return ``True`` since this ring is Noetherian.

        EXAMPLES::

            sage: ZZ.is_noetherian()
            True
            sage: QQ.is_noetherian()
            True
            sage: R.<x> = PolynomialRing(QQ)
            sage: R.is_noetherian()
            True
        """
        return True

cdef class DedekindDomain(IntegralDomain):
    """
    Generic Dedekind domain class.

    A Dedekind domain is a Noetherian integral domain of Krull
    dimension one that is integrally closed in its field of fractions.

    This class is deprecated, and not actually used anywhere in the
    Sage code base.  If you think you need it, please create a
    category :class:`DedekindDomains`, move the code of this class
    there, and use it instead.
    """
    def krull_dimension(self):
        """
        Return 1 since Dedekind domains have Krull dimension 1.

        EXAMPLES:

        The following are examples of Dedekind domains (Noetherian integral
        domains of Krull dimension one that are integrally closed over its
        field of fractions)::

            sage: ZZ.krull_dimension()
            1
            sage: K = NumberField(x^2 + 1, 's')
            sage: OK = K.ring_of_integers()
            sage: OK.krull_dimension()
            1

        The following are not Dedekind domains but have
        a ``krull_dimension`` function::

            sage: QQ.krull_dimension()
            0
            sage: T.<x,y> = PolynomialRing(QQ,2); T
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: T.krull_dimension()
            2
            sage: U.<x,y,z> = PolynomialRing(ZZ,3); U
            Multivariate Polynomial Ring in x, y, z over Integer Ring
            sage: U.krull_dimension()
            4

            sage: K.<i> = QuadraticField(-1)
            sage: R = K.order(2*i); R
            Order in Number Field in i with defining polynomial x^2 + 1 with i = 1*I
            sage: R.is_maximal()
            False
            sage: R.krull_dimension()
            1
        """
        return 1

    def is_integrally_closed(self):
        """
        Return ``True`` since Dedekind domains are integrally closed.

        EXAMPLES:

        The following are examples of Dedekind domains (Noetherian integral
        domains of Krull dimension one that are integrally closed over its
        field of fractions).

        ::

            sage: ZZ.is_integrally_closed()
            True
            sage: K = NumberField(x^2 + 1, 's')
            sage: OK = K.ring_of_integers()
            sage: OK.is_integrally_closed()
            True

        These, however, are not Dedekind domains::

            sage: QQ.is_integrally_closed()
            True
            sage: S = ZZ[sqrt(5)]; S.is_integrally_closed()
            False
            sage: T.<x,y> = PolynomialRing(QQ,2); T
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: T.is_integral_domain()
            True
        """
        return True

    def integral_closure(self):
        r"""
        Return ``self`` since Dedekind domains are integrally closed.

        EXAMPLES::

            sage: K = NumberField(x^2 + 1, 's')
            sage: OK = K.ring_of_integers()
            sage: OK.integral_closure()
            Gaussian Integers in Number Field in s with defining polynomial x^2 + 1
            sage: OK.integral_closure() == OK
            True

            sage: QQ.integral_closure() == QQ
            True
        """
        return self

    def is_noetherian(self):
        r"""
        Return ``True`` since Dedekind domains are Noetherian.

        EXAMPLES:

        The integers, `\ZZ`, and rings of integers of number
        fields are Dedekind domains::

            sage: ZZ.is_noetherian()
            True
            sage: K = NumberField(x^2 + 1, 's')
            sage: OK = K.ring_of_integers()
            sage: OK.is_noetherian()
            True
            sage: QQ.is_noetherian()
            True
        """
        return True


cdef class PrincipalIdealDomain(IntegralDomain):
    """
    Generic principal ideal domain.

    This class is deprecated. Please use the
    :class:`~sage.categories.principal_ideal_domains.PrincipalIdealDomains`
    category instead.
    """
    _default_category = PrincipalIdealDomains()

    def is_noetherian(self):
        """
        Every principal ideal domain is noetherian, so we return ``True``.

        EXAMPLES::

            sage: Zp(5).is_noetherian()
            True
        """
        return True

    def class_group(self):
        """
        Return the trivial group, since the class group of a PID is trivial.

        EXAMPLES::

            sage: QQ.class_group()
            Trivial Abelian group
        """
        from sage.groups.abelian_gps.abelian_group import AbelianGroup
        return AbelianGroup([])

    def gcd(self, x, y, coerce=True):
        r"""
        Return the greatest common divisor of ``x`` and ``y``, as elements
        of ``self``.

        EXAMPLES:

        The integers are a principal ideal domain and hence a GCD domain::

            sage: ZZ.gcd(42, 48)
            6
            sage: 42.factor(); 48.factor()
            2 * 3 * 7
            2^4 * 3
            sage: ZZ.gcd(2^4*7^2*11, 2^3*11*13)
            88
            sage: 88.factor()
            2^3 * 11

        In a field, any nonzero element is a GCD of any nonempty set
        of nonzero elements. In previous versions, Sage used to return
        1 in the case of the rational field. However, since :trac:`10771`,
        the rational field is considered as the
        *fraction field* of the integer ring. For the fraction field
        of an integral domain that provides both GCD and LCM, it is
        possible to pick a GCD that is compatible with the GCD of the
        base ring::

            sage: QQ.gcd(ZZ(42), ZZ(48)); type(QQ.gcd(ZZ(42), ZZ(48)))
            6
            <class 'sage.rings.rational.Rational'>
            sage: QQ.gcd(1/2, 1/3)
            1/6

        Polynomial rings over fields are GCD domains as well. Here is a simple
        example over the ring of polynomials over the rationals as well as
        over an extension ring. Note that ``gcd`` requires x and y to be
        coercible::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = NumberField(x^2 - 2, 'a')
            sage: f = (x - a)*(x + a); g = (x - a)*(x^2 - 2)
            sage: print(f); print(g)
            x^2 - 2
            x^3 - a*x^2 - 2*x + 2*a
            sage: f in R
            True
            sage: g in R
            False
            sage: R.gcd(f,g)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce 2*a to a rational
            sage: R.base_extend(S).gcd(f,g)
            x^2 - 2
            sage: R.base_extend(S).gcd(f, (x - a)*(x^2 - 3))
            x - a
        """
        if coerce:
            x = self(x)
            y = self(y)
        return x.gcd(y)

    def content(self, x, y, coerce=True):
        r"""
        Return the content of `x` and `y`, i.e. the unique element `c` of
        ``self`` such that `x/c` and `y/c` are coprime and integral.

        EXAMPLES::

            sage: QQ.content(ZZ(42), ZZ(48)); type(QQ.content(ZZ(42), ZZ(48)))
            6
            <class 'sage.rings.rational.Rational'>
            sage: QQ.content(1/2, 1/3)
            1/6
            sage: factor(1/2); factor(1/3); factor(1/6)
            2^-1
            3^-1
            2^-1 * 3^-1
            sage: a = (2*3)/(7*11); b = (13*17)/(19*23)
            sage: factor(a); factor(b); factor(QQ.content(a,b))
            2 * 3 * 7^-1 * 11^-1
            13 * 17 * 19^-1 * 23^-1
            7^-1 * 11^-1 * 19^-1 * 23^-1

        Note the changes to the second entry::

            sage: c = (2*3)/(7*11); d = (13*17)/(7*19*23)
            sage: factor(c); factor(d); factor(QQ.content(c,d))
            2 * 3 * 7^-1 * 11^-1
            7^-1 * 13 * 17 * 19^-1 * 23^-1
            7^-1 * 11^-1 * 19^-1 * 23^-1
            sage: e = (2*3)/(7*11); f = (13*17)/(7^3*19*23)
            sage: factor(e); factor(f); factor(QQ.content(e,f))
            2 * 3 * 7^-1 * 11^-1
            7^-3 * 13 * 17 * 19^-1 * 23^-1
            7^-3 * 11^-1 * 19^-1 * 23^-1
        """
        if coerce:
            x = self(x)
            y = self(y)
        return x.content(y)

    def _ideal_class_(self, n=0):
        """
        Ideals in PIDs have their own special class.

        EXAMPLES::

            sage: ZZ._ideal_class_()
            <class 'sage.rings.ideal.Ideal_pid'>
        """
        from sage.rings.ideal import Ideal_pid
        return Ideal_pid

cdef class EuclideanDomain(PrincipalIdealDomain):
    """
    Generic Euclidean domain class.

    This class is deprecated. Please use the
    :class:`~sage.categories.euclidean_domains.EuclideanDomains`
    category instead.
    """
    _default_category = EuclideanDomains()

    def parameter(self):
        """
        Return an element of degree 1.

        EXAMPLES::

            sage: R.<x>=QQ[]
            sage: R.parameter()
            x
       """
        raise NotImplementedError

cpdef bint _is_Field(x) except -2:
    """
    Return ``True`` if ``x`` is a field.

    EXAMPLES::

        sage: from sage.rings.ring import _is_Field
        sage: _is_Field(QQ)
        True
        sage: _is_Field(ZZ)
        False
        sage: _is_Field(pAdicField(2))
        True
        sage: _is_Field(5)
        False

    NOTE:

    ``_is_Field(R)`` is of internal use. It is better (and faster) to
    use ``R in Fields()`` instead.
    """
    # The result is not immediately returned, since we want to refine
    # x's category, so that calling x in Fields() will be faster next time.
    try:
        result = isinstance(x, Field) or x.is_field()
    except AttributeError:
        result = False
    if result:
        x._refine_category_(_Fields)
    return result

from sage.categories.algebras import Algebras
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.fields import Fields
_Fields = Fields()

cdef class Field(PrincipalIdealDomain):
    """
    Generic field
    """
    _default_category = _Fields

    def fraction_field(self):
        """
        Return the fraction field of ``self``.

        EXAMPLES:

        Since fields are their own field of fractions, we simply get the
        original field in return::

            sage: QQ.fraction_field()
            Rational Field
            sage: RR.fraction_field()
            Real Field with 53 bits of precision
            sage: CC.fraction_field()
            Complex Field with 53 bits of precision

            sage: F = NumberField(x^2 + 1, 'i')
            sage: F.fraction_field()
            Number Field in i with defining polynomial x^2 + 1
        """
        return self

    def _pseudo_fraction_field(self):
        """
        The fraction field of ``self`` is always available as ``self``.

        EXAMPLES::

            sage: QQ._pseudo_fraction_field()
            Rational Field
            sage: K = GF(5)
            sage: K._pseudo_fraction_field()
            Finite Field of size 5
            sage: K._pseudo_fraction_field() is K
            True
        """
        return self

    def divides(self, x, y, coerce=True):
        """
        Return ``True`` if ``x`` divides ``y`` in this field (usually ``True``
        in a field!).  If ``coerce`` is ``True`` (the default), first coerce
        ``x`` and ``y`` into ``self``.

        EXAMPLES::

            sage: QQ.divides(2, 3/4)
            True
            sage: QQ.divides(0, 5)
            False
        """
        if coerce:
            x = self(x)
            y = self(y)
        if x.is_zero():
            return y.is_zero()
        return True

    def ideal(self, *gens, **kwds):
        """
        Return the ideal generated by gens.

        EXAMPLES::

            sage: QQ.ideal(2)
            Principal ideal (1) of Rational Field
            sage: QQ.ideal(0)
            Principal ideal (0) of Rational Field
        """
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        if not isinstance(gens, (list, tuple)):
            gens = [gens]
        for x in gens:
            if not self(x).is_zero():
                return self.unit_ideal()
        return self.zero_ideal()

    def integral_closure(self):
        """
        Return this field, since fields are integrally closed in their
        fraction field.

        EXAMPLES::

            sage: QQ.integral_closure()
            Rational Field
            sage: Frac(ZZ['x,y']).integral_closure()
            Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring
        """
        return self

    def is_field(self, proof = True):
        """
        Return ``True`` since this is a field.

        EXAMPLES::

            sage: Frac(ZZ['x,y']).is_field()
            True
        """
        return True

    def is_integrally_closed(self):
        """
        Return ``True`` since fields are trivially integrally closed in
        their fraction field (since they are their own fraction field).

        EXAMPLES::

            sage: Frac(ZZ['x,y']).is_integrally_closed()
            True
        """
        return True

    def is_noetherian(self):
        """
        Return ``True`` since fields are Noetherian rings.

        EXAMPLES::

            sage: QQ.is_noetherian()
            True
        """
        return True

    def krull_dimension(self):
        """
        Return the Krull dimension of this field, which is 0.

        EXAMPLES::

            sage: QQ.krull_dimension()
            0
            sage: Frac(QQ['x,y']).krull_dimension()
            0
        """
        return 0

    def prime_subfield(self):
        """
        Return the prime subfield of ``self``.

        EXAMPLES::

            sage: k = GF(9, 'a')
            sage: k.prime_subfield()
            Finite Field of size 3
        """
        if self.characteristic() == 0:
            import sage.rings.rational_field
            return sage.rings.rational_field.RationalField()
        else:
            from sage.rings.finite_rings.finite_field_constructor import GF
            return GF(self.characteristic())

    def algebraic_closure(self):
        """
        Return the algebraic closure of ``self``.

        .. NOTE::

           This is only implemented for certain classes of field.

        EXAMPLES::

            sage: K = PolynomialRing(QQ,'x').fraction_field(); K
            Fraction Field of Univariate Polynomial Ring in x over Rational Field
            sage: K.algebraic_closure()
            Traceback (most recent call last):
            ...
            NotImplementedError: Algebraic closures of general fields not implemented.
        """
        raise NotImplementedError("Algebraic closures of general fields not implemented.")


cdef class Algebra(Ring):
    """
    Generic algebra
    """
    def __init__(self, base_ring, names=None, normalize=True, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: A = Algebra(ZZ); A
            <sage.rings.ring.Algebra object at ...>
        """
        # This is a low-level class. For performance, we trust that the category
        # is fine, if it is provided. If it isn't, we use the category of Algebras(base_ring).
        if category is None:
            category = Algebras(base_ring)
        Ring.__init__(self,base_ring, names=names, normalize=normalize,
                      category=category)

    def characteristic(self):
        r"""
        Return the characteristic of this algebra, which is the same
        as the characteristic of its base ring.

        See objects with the ``base_ring`` attribute for additional examples.
        Here are some examples that explicitly use the :class:`Algebra` class.

        EXAMPLES::

            sage: A = Algebra(ZZ); A
            <sage.rings.ring.Algebra object at ...>
            sage: A.characteristic()
            0
            sage: A = Algebra(GF(7^3, 'a'))
            sage: A.characteristic()
            7
        """
        return self.base_ring().characteristic()

    def has_standard_involution(self):
        r"""
        Return ``True`` if the algebra has a standard involution and ``False`` otherwise.
        This algorithm follows Algorithm 2.10 from John Voight's *Identifying the Matrix Ring*.
        Currently the only type of algebra this will work for is a quaternion algebra.
        Though this function seems redundant, once algebras have more functionality, in particular
        have a method to construct a basis, this algorithm will have more general purpose.

        EXAMPLES::

            sage: B = QuaternionAlgebra(2)
            sage: B.has_standard_involution()
            True
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<u> = NumberField(x**2 - 2)
            sage: A = QuaternionAlgebra(K,-2,5)
            sage: A.has_standard_involution()
            True
            sage: L.<a,b> = FreeAlgebra(QQ,2)
            sage: L.has_standard_involution()
            Traceback (most recent call last):
            ...
            NotImplementedError: has_standard_involution is not implemented for this algebra
            """
        field = self.base_ring()
        try:
            basis = self.basis()
        except AttributeError:
            raise AttributeError("Basis is not yet implemented for this algebra.")
        try:
            # TODO: The following code is specific to the quaternion algebra
            #   and should belong there
            #step 1
            for i in range(1,4):
                ei = basis[i]
                a = ei**2
                coef = a.coefficient_tuple()
                ti = coef[i]
                ni = a - ti*ei
                if ni not in field:
                    return False
            #step 2
            for i in range(1,4):
                for j in range(2,4):
                    ei = basis[i]
                    ej = basis[j]
                    a = ei**2
                    coef = a.coefficient_tuple()
                    ti = coef[i]
                    b = ej**2
                    coef = b.coefficient_tuple()
                    tj = coef[j]
                    nij = (ei + ej)**2 - (ti + tj)*(ei + ej)
                    if nij not in field:
                        return False
        except AttributeError:
            raise NotImplementedError("has_standard_involution is not implemented for this algebra")
        return True

cdef class CommutativeAlgebra(CommutativeRing):
    """
    Generic commutative algebra
    """
    def __init__(self, base_ring, names=None, normalize=True, category=None):
        r"""
        Standard init function. This just checks that the base is a commutative
        ring and then passes the buck.

        EXAMPLES::

            sage: sage.rings.ring.CommutativeAlgebra(QQ)
            <sage.rings.ring.CommutativeAlgebra object at ...>

            sage: sage.rings.ring.CommutativeAlgebra(QuaternionAlgebra(QQ,-1,-1))
            Traceback (most recent call last):
            ...
            TypeError: base ring must be a commutative ring
        """
        # TODO: use the idiom base_ring in CommutativeRings()
        try:
            if not base_ring.is_commutative():
                raise TypeError("base ring must be a commutative ring")
        except (AttributeError, NotImplementedError):
            raise TypeError("base ring must be a commutative ring")
        # This is a low-level class. For performance, we trust that
        # the category is fine, if it is provided. If it isn't, we use
        # the category of commutative algebras.
        if category is None:
            category = CommutativeAlgebras(base_ring)
        CommutativeRing.__init__(self, base_ring, names=names, normalize=normalize, category=category)

    def is_commutative(self):
        """
        Return ``True`` since this algebra is commutative.

        EXAMPLES:

        Any commutative ring is a commutative algebra over itself::

            sage: A = sage.rings.ring.CommutativeAlgebra
            sage: A(ZZ).is_commutative()
            True
            sage: A(QQ).is_commutative()
            True

        Trying to create a commutative algebra over a non-commutative ring
        will result in a ``TypeError``.
        """
        return True


def is_Ring(x):
    """
    Return ``True`` if ``x`` is a ring.

    EXAMPLES::

        sage: from sage.rings.ring import is_Ring
        sage: is_Ring(ZZ)
        True
        sage: MS = MatrixSpace(QQ,2)
        sage: is_Ring(MS)
        True
    """
    return x in _Rings
