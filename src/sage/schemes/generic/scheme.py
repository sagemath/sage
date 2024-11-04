"""
Schemes

AUTHORS:

- William Stein, David Kohel, Kiran Kedlaya (2008): added zeta_series

- Volker Braun (2011-08-11): documenting, improving, refactoring.
"""
# ****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2008 Kiran Kedlaya <kedlaya@mit.edu>
#       Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2005 William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.commutative_rings import CommutativeRings
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.rings.ideal import Ideal_generic
from sage.rings.integer_ring import ZZ
from sage.schemes.generic.point import SchemeTopologicalPoint_prime_ideal
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

lazy_import('sage.schemes.generic.morphism', 'SchemeMorphism')
lazy_import('sage.schemes.elliptic_curves.ell_generic', 'EllipticCurve_generic', as_='EllipticCurve')


def is_Scheme(x):
    """
    Test whether ``x`` is a scheme.

    INPUT:

    - ``x`` -- anything

    OUTPUT: boolean; whether ``x`` derives from :class:`Scheme`

    EXAMPLES::

        sage: from sage.schemes.generic.scheme import is_Scheme
        sage: is_Scheme(5)
        doctest:warning...
        DeprecationWarning: The function is_Scheme is deprecated; use 'isinstance(..., Scheme)' or categories instead.
        See https://github.com/sagemath/sage/issues/38022 for details.
        False
        sage: X = Spec(QQ)
        sage: is_Scheme(X)
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38022, "The function is_Scheme is deprecated; use 'isinstance(..., Scheme)' or categories instead.")
    return isinstance(x, Scheme)


class Scheme(Parent):
    r"""
    The base class for all schemes.

    INPUT:

    - ``X`` -- a scheme, scheme morphism, commutative ring,
      commutative ring morphism, or ``None`` (optional). Determines
      the base scheme. If a commutative ring is passed, the spectrum
      of the ring will be used as base.

    - ``category`` -- the category (optional); will be automatically
      constructed by default

    EXAMPLES::

        sage: from sage.schemes.generic.scheme import Scheme
        sage: Scheme(ZZ)
        <sage.schemes.generic.scheme.Scheme_with_category object at ...>

    A scheme is in the category of all schemes over its base::

        sage: ProjectiveSpace(4, QQ).category()
        Category of schemes over Rational Field

    There is a special and unique `\mathrm{Spec}(\ZZ)` that is the default base
    scheme::

        sage: Spec(ZZ).base_scheme() is Spec(QQ).base_scheme()
        True
    """

    def __init__(self, X=None, category=None):
        """
        Construct a scheme.

        TESTS:

        The full test suite works since :issue:`7946`::

            sage: R.<x, y> = QQ[]
            sage: I = (x^2 - y^2)*R
            sage: RmodI = R.quotient(I)
            sage: X = Spec(RmodI)
            sage: TestSuite(X).run()                                                    # needs sage.libs.singular
        """
        from sage.categories.map import Map
        from sage.categories.rings import Rings

        if X is None:
            self._base_ring = ZZ
        elif isinstance(X, Scheme):
            self._base_scheme = X
        elif isinstance(X, SchemeMorphism):
            self._base_morphism = X
        elif X in CommutativeRings():
            self._base_ring = X
        elif isinstance(X, Map) and X.category_for().is_subcategory(Rings()):
            # X is a morphism of Rings
            self._base_ring = X.codomain()
        else:
            raise ValueError('The base must be define by a scheme, '
                             'scheme morphism, or commutative ring.')

        from sage.categories.schemes import Schemes
        if X is None:
            default_category = Schemes()
        else:
            default_category = Schemes(self.base_scheme())
        if category is None:
            category = default_category
        else:
            assert category.is_subcategory(default_category), \
                "%s is not a subcategory of %s" % (category, default_category)

        Parent.__init__(self, self.base_ring(), category=category)

    def union(self, X):
        """
        Return the disjoint union of the schemes ``self`` and ``X``.

        EXAMPLES::

            sage: S = Spec(QQ)
            sage: X = AffineSpace(1, QQ)
            sage: S.union(X)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    __add__ = union

    def _morphism(self, *args, **kwds):
        """
        Construct a morphism determined by action on points of ``self``.

        EXAMPLES::

            sage: X = Spec(QQ)
            sage: X._morphism()
            Traceback (most recent call last):
            ...
            NotImplementedError

        TESTS:

        This shows that issue at :issue:`7389` is solved::

            sage: S = Spec(ZZ)
            sage: f = S.identity_morphism()
            sage: from sage.schemes.generic.glue import GluedScheme
            sage: T = GluedScheme(f, f)
            sage: S.hom([1],T)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def base_extend(self, Y):
        """
        Extend the base of the scheme.

        Derived classes must override this method.

        EXAMPLES::

            sage: from sage.schemes.generic.scheme import Scheme
            sage: X = Scheme(ZZ)
            sage: X.base_scheme()
            Spectrum of Integer Ring
            sage: X.base_extend(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __call__(self, *args):
        """
        Call syntax for schemes.

        INPUT/OUTPUT: the arguments must be one of the following:

        - a ring or a scheme `S`. Output will be the set `X(S)` of
          `S`-valued points on `X`.

        - If `S` is a list or tuple or just the coordinates, return a
          point in `X(T)`, where `T` is the base scheme of ``self``.

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)

        We create some point sets::

            sage: A(QQ)
            Set of rational points of Affine Space of dimension 2 over Rational Field
            sage: A(RR)                                                                 # needs sage.rings.real_mpfr
            Set of rational points of Affine Space of dimension 2
             over Real Field with 53 bits of precision

        Space of dimension 2 over Rational Field::

            sage: R.<x> = PolynomialRing(QQ)
            sage: A(NumberField(x^2 + 1, 'a'))                                          # needs sage.rings.number_field
            Set of rational points of Affine Space of dimension 2
             over Number Field in a with defining polynomial x^2 + 1
            sage: A(GF(7))
            Traceback (most recent call last):
            ...
            ValueError: There must be a natural map S --> R, but
            S = Rational Field and R = Finite Field of size 7

        We create some points::

            sage: A(QQ)([1, 0])
            (1, 0)

        We create the same point by giving the coordinates of the point
        directly::

            sage: A(1, 0)
            (1, 0)

        Check that :issue:`16832` is fixed::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: X(P([4, 4, 1]))
            (4 : 4 : 1)
        """
        if len(args) == 1:
            from sage.schemes.generic.morphism import SchemeMorphism_point
            S = args[0]
            if S in CommutativeRings():
                return self.point_homset(S)
            elif isinstance(S, Scheme):
                return S.Hom(self)
            elif isinstance(S, (list, tuple)):
                args = S
            elif isinstance(S, SchemeMorphism_point):
                if S.codomain() is self:
                    return S
                args = S
        return self.point(args)

    @cached_method
    def point_homset(self, S=None):
        r"""
        Return the set of S-valued points of this scheme.

        INPUT:

        - ``S`` -- a commutative ring

        OUTPUT: the set of morphisms `\mathrm{Spec}(S) \to X`

        EXAMPLES::

            sage: P = ProjectiveSpace(ZZ, 3)
            sage: P.point_homset(ZZ)
            Set of rational points of Projective Space of dimension 3 over Integer Ring
            sage: P.point_homset(QQ)
            Set of rational points of Projective Space of dimension 3 over Rational Field
            sage: P.point_homset(GF(11))
            Set of rational points of Projective Space of dimension 3 over
             Finite Field of size 11

        TESTS::

            sage: P = ProjectiveSpace(QQ, 3)
            sage: P.point_homset(GF(11))
            Traceback (most recent call last):
            ...
            ValueError: There must be a natural map S --> R, but
            S = Rational Field and R = Finite Field of size 11
        """
        if S is None:
            S = self.base_ring()
        SpecS = AffineScheme(S, self.base_ring())
        from sage.schemes.generic.homset import SchemeHomset
        return SchemeHomset(SpecS, self, as_point_homset=True)

    def point(self, v, check=True):
        """
        Create a point.

        INPUT:

        - ``v`` -- anything that defines a point

        - ``check`` -- boolean (default: ``True``); whether
          to check the defining data for consistency

        OUTPUT: a point of the scheme

        EXAMPLES::

            sage: A2 = AffineSpace(QQ, 2)
            sage: A2.point([4, 5])
            (4, 5)

            sage: R.<t> = PolynomialRing(QQ)
            sage: E = EllipticCurve([t + 1, t, t, 0, 0])                                # needs sage.schemes
            sage: E.point([0, 0])                                                       # needs sage.schemes
            (0 : 0 : 1)
        """
        # todo: update elliptic curve stuff to take point_homset as argument
        if isinstance(self, EllipticCurve):
            try:
                return self._point(self.point_homset(), v, check=check)
            except AttributeError:  # legacy code without point_homset
                return self._point(self, v, check=check)

        return self.point_homset()(v, check=check)

    def _point(self):
        """
        Return the Hom-set from some affine scheme to ``self``.

        OUTPUT:

        A scheme Hom-set, see :mod:`~sage.schemes.generic.homset`.

        EXAMPLES::

            sage: X = Spec(QQ)
            sage: X._point()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _point_homset(self, *args, **kwds):
        """
        Return the Hom-set from ``self`` to another scheme.

        EXAMPLES::

            sage: from sage.schemes.generic.scheme import Scheme
            sage: X = Scheme(QQ)
            sage: X._point_homset()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __truediv__(self, Y):
        """
        Return the base extension of ``self`` to Y.

        See :meth:`base_extend` for details.

        EXAMPLES::

            sage: A = AffineSpace(3, ZZ)
            sage: A
            Affine Space of dimension 3 over Integer Ring
            sage: A/QQ
            Affine Space of dimension 3 over Rational Field
            sage: A/GF(7)
            Affine Space of dimension 3 over Finite Field of size 7
        """
        return self.base_extend(Y)

    def base_ring(self):
        """
        Return the base ring of the scheme ``self``.

        OUTPUT: a commutative ring

        EXAMPLES::

            sage: A = AffineSpace(4, QQ)
            sage: A.base_ring()
            Rational Field

            sage: X = Spec(QQ)
            sage: X.base_ring()
            Integer Ring
        """
        try:
            return self._base_ring
        except AttributeError:
            if hasattr(self, '_base_morphism'):
                self._base_ring = self._base_morphism.codomain().coordinate_ring()
            elif hasattr(self, '_base_scheme'):
                self._base_ring = self._base_scheme.coordinate_ring()
            else:
                self._base_ring = ZZ
            return self._base_ring

    def base_scheme(self):
        """
        Return the base scheme.

        OUTPUT: a scheme

        EXAMPLES::

            sage: A = AffineSpace(4, QQ)
            sage: A.base_scheme()
            Spectrum of Rational Field

            sage: X = Spec(QQ)
            sage: X.base_scheme()
            Spectrum of Integer Ring
        """
        try:
            return self._base_scheme
        except AttributeError:
            if hasattr(self, '_base_morphism'):
                self._base_scheme = self._base_morphism.codomain()
            elif hasattr(self, '_base_ring'):
                self._base_scheme = AffineScheme(self._base_ring)
            else:
                from sage.schemes.generic.spec import SpecZ
                self._base_scheme = SpecZ
            return self._base_scheme

    def base_morphism(self):
        """
        Return the structure morphism from ``self`` to its base
        scheme.

        OUTPUT: a scheme morphism

        EXAMPLES::

            sage: A = AffineSpace(4, QQ)
            sage: A.base_morphism()
            Scheme morphism:
              From: Affine Space of dimension 4 over Rational Field
              To:   Spectrum of Rational Field
              Defn: Structure map

            sage: X = Spec(QQ)
            sage: X.base_morphism()
            Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Structure map
        """
        try:
            return self._base_morphism
        except AttributeError:
            from sage.categories.schemes import Schemes
            from sage.schemes.generic.spec import SpecZ
            SCH = Schemes()
            if hasattr(self, '_base_scheme'):
                self._base_morphism = self.Hom(self._base_scheme, category=SCH).natural_map()
            elif hasattr(self, '_base_ring'):
                self._base_morphism = self.Hom(AffineScheme(self._base_ring), category=SCH).natural_map()
            else:
                self._base_morphism = self.Hom(SpecZ, category=SCH).natural_map()
            return self._base_morphism

    structure_morphism = base_morphism

    def coordinate_ring(self):
        """
        Return the coordinate ring.

        OUTPUT:

        The global coordinate ring of this scheme, if
        defined. Otherwise this raises a :exc:`ValueError`.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: I = (x^2 - y^2)*R
            sage: X = Spec(R.quotient(I))
            sage: X.coordinate_ring()
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field
             by the ideal (x^2 - y^2)
        """
        try:
            return self._coordinate_ring
        except AttributeError:
            raise ValueError("This scheme has no associated coordinated ring (defined).")

    def dimension_absolute(self):
        """
        Return the absolute dimension of this scheme.

        OUTPUT: integer

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: I = (x^2 - y^2)*R
            sage: X = Spec(R.quotient(I))
            sage: X.dimension_absolute()
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: X.dimension()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError # override in derived class

    dimension = dimension_absolute

    def dimension_relative(self):
        """
        Return the relative dimension of this scheme over its base.

        OUTPUT: integer

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: I = (x^2 - y^2)*R
            sage: X = Spec(R.quotient(I))
            sage: X.dimension_relative()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError # override in derived class

    def identity_morphism(self):
        """
        Return the identity morphism.

        OUTPUT: the identity morphism of the scheme ``self``

        EXAMPLES::

            sage: X = Spec(QQ)
            sage: X.identity_morphism()
            Scheme endomorphism of Spectrum of Rational Field
              Defn: Identity map
        """
        from sage.schemes.generic.morphism import SchemeMorphism_id
        return SchemeMorphism_id(self)

    def hom(self, x, Y=None, check=True):
        """
        Return the scheme morphism from ``self`` to ``Y`` defined by ``x``.

        INPUT:

        - ``x`` -- anything that determines a scheme morphism; if
          ``x`` is a scheme, try to determine a natural map to ``x``

        - ``Y`` -- the codomain scheme (optional); if ``Y`` is not
          given, try to determine ``Y`` from context

        - ``check`` -- boolean (default: ``True``); whether
          to check the defining data for consistency

        OUTPUT: the scheme morphism from ``self`` to ``Y`` defined by ``x``

        EXAMPLES::

            sage: P = ProjectiveSpace(ZZ, 3)
            sage: P.hom(Spec(ZZ))
            Scheme morphism:
              From: Projective Space of dimension 3 over Integer Ring
              To:   Spectrum of Integer Ring
              Defn: Structure map
        """
        if Y is None:
            if isinstance(x, Scheme):
                return self.Hom(x).natural_map()
            else:
                raise TypeError("unable to determine codomain")
        return self.Hom(Y)(x, check=check)

    def _Hom_(self, Y, category=None, check=True):
        """
        Return the set of scheme morphisms from ``self`` to ``Y``.

        INPUT:

        - ``Y`` -- a scheme; the codomain of the Hom-set

        - ``category`` -- a category (optional); the category of the
          Hom-set

        - ``check`` -- boolean (default: ``True``); whether
          to check the defining data for consistency

        OUTPUT: the set of morphisms from ``self`` to ``Y``

        EXAMPLES::

            sage: P = ProjectiveSpace(ZZ, 3)
            sage: S = Spec(ZZ)
            sage: S._Hom_(P)
            Set of morphisms
              From: Spectrum of Integer Ring
              To:   Projective Space of dimension 3 over Integer Ring

        TESTS::

            sage: S._Hom_(P).__class__
            <class 'sage.schemes.generic.homset.SchemeHomset_generic_with_category'>

            sage: Hom(Spec(ZZ), Spec(ZZ)).__class__
            <class 'sage.schemes.generic.homset.SchemeHomset_generic_with_category_with_equality_by_id'>
        """
        from sage.schemes.generic.homset import SchemeHomset
        return SchemeHomset(self, Y, category=category, check=check)

    point_set = point_homset

    def count_points(self, n):
        r"""
        Count points over finite fields.

        INPUT:

        - ``n`` -- integer

        OUTPUT:

        An integer. The number of points over `\GF{q}, \ldots,
        \GF{q^n}` on a scheme over a finite field `\GF{q}`.

        EXAMPLES::

            sage: # needs sage.schemes
            sage: P.<x> = PolynomialRing(GF(3))
            sage: C = HyperellipticCurve(x^3 + x^2 + 1)
            sage: C.count_points(4)
            [6, 12, 18, 96]
            sage: C.base_extend(GF(9,'a')).count_points(2)                              # needs sage.rings.finite_rings
            [12, 96]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(4, 't'), 2)                            # needs sage.rings.finite_rings
            sage: X = P.subscheme([y^2*z - x^3 - z^3])                                  # needs sage.rings.finite_rings
            sage: X.count_points(2)                                                     # needs sage.libs.singular sage.rings.finite_rings
            [5, 17]
        """
        F = self.base_ring()
        if not F.is_finite():
            raise TypeError("Point counting only defined for schemes over finite fields")
        a = [len(self.rational_points())]
        for i in range(2, n+1):
            F1, psi = F.extension(i, map=True)
            S1 = self.change_ring(psi)
            a.append(len(S1.rational_points()))
        return a

    def zeta_function(self):
        r"""
        Compute the zeta function of a generic scheme.

        Derived classes should override this method.

        OUTPUT: rational function in one variable

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(GF(4, 't'), 2)                            # needs sage.rings.finite_rings
            sage: X = P.subscheme([y^2*z - x^3 - z^3])                                  # needs sage.rings.finite_rings
            sage: X.zeta_function()                                                     # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def zeta_series(self, n, t):
        """
        Return the zeta series.

        Compute a power series approximation to the zeta function of a
        scheme over a finite field.

        INPUT:

        - ``n`` -- the number of terms of the power series to compute

        - ``t`` -- the variable which the series should be returned

        OUTPUT: a power series approximating the zeta function of ``self``

        EXAMPLES::

            sage: P.<x> = PolynomialRing(GF(3))
            sage: C = HyperellipticCurve(x^3 + x^2 + 1)                                 # needs sage.schemes
            sage: R.<t> = PowerSeriesRing(Integers())
            sage: C.zeta_series(4, t)                                                   # needs sage.schemes
            1 + 6*t + 24*t^2 + 78*t^3 + 240*t^4 + O(t^5)
            sage: (1+2*t+3*t^2)/(1-t)/(1-3*t) + O(t^5)
            1 + 6*t + 24*t^2 + 78*t^3 + 240*t^4 + O(t^5)

        If the scheme has a method ``zeta_function``, this is used to
        provide the required approximation.
        Otherwise this function depends on ``count_points``, which is only
        defined for prime order fields for general schemes.
        Nonetheless, since :issue:`15108` and :issue:`15148`, it supports
        hyperelliptic curves over non-prime fields::

            sage: C.base_extend(GF(9, 'a')).zeta_series(4, t)                           # needs sage.rings.finite_rings sage.schemes
            1 + 12*t + 120*t^2 + 1092*t^3 + 9840*t^4 + O(t^5)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(4, 't'), 2)                            # needs sage.rings.finite_rings
            sage: X = P.subscheme([y^2*z - x^3 - z^3])                                  # needs sage.rings.finite_rings

            sage: R.<t> = PowerSeriesRing(Integers())
            sage: X.zeta_series(2, t)                                                   # needs sage.libs.singular sage.rings.finite_rings
            1 + 5*t + 21*t^2 + O(t^3)

        TESTS::

            sage: P.<x> = PolynomialRing(ZZ)
            sage: C = HyperellipticCurve(x^3 + x + 1)                                   # needs sage.schemes
            sage: R.<t> = PowerSeriesRing(Integers())
            sage: C.zeta_series(4, t)                                                   # needs sage.schemes
            Traceback (most recent call last):
            ...
            TypeError: zeta functions only defined for schemes
            over finite fields
        """
        F = self.base_ring()
        if not F.is_finite():
            raise TypeError('zeta functions only defined for schemes over finite fields')
        R = t.parent()
        u = t.O(n + 1)
        try:
            return self.zeta_function()(u)
        except (AttributeError, NotImplementedError):
            pass
        try:
            a = self.count_points(n)
        except AttributeError:
            raise NotImplementedError('count_points() required but not implemented')
        temp = R.sum(a[i - 1] * u**i / i for i in range(1, n + 1))
        return temp.exp()


def is_AffineScheme(x):
    """
    Return ``True`` if `x` is an affine scheme.

    EXAMPLES::

        sage: from sage.schemes.generic.scheme import is_AffineScheme
        sage: is_AffineScheme(5)
        doctest:warning...
        DeprecationWarning: The function is_AffineScheme is deprecated; use 'isinstance(..., AffineScheme)' instead.
        See https://github.com/sagemath/sage/issues/38022 for details.
        False
        sage: E = Spec(QQ)
        sage: is_AffineScheme(E)
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38022, "The function is_AffineScheme is deprecated; use 'isinstance(..., AffineScheme)' instead.")
    return isinstance(x, AffineScheme)


class AffineScheme(UniqueRepresentation, Scheme):
    """
    Class for general affine schemes.

    TESTS::

        sage: from sage.schemes.generic.scheme import AffineScheme
        sage: A = QQ['t']
        sage: X_abs = AffineScheme(A); X_abs
        Spectrum of Univariate Polynomial Ring in t over Rational Field
        sage: X_rel = AffineScheme(A, QQ); X_rel
        Spectrum of Univariate Polynomial Ring in t over Rational Field

        sage: X_abs == X_rel
        False
        sage: X_abs.base_ring()
        Integer Ring
        sage: X_rel.base_ring()
        Rational Field

    .. SEEALSO::

        For affine spaces over a base ring and subschemes thereof, see
        :class:`sage.schemes.generic.algebraic_scheme.AffineSpace`.
    """
    def __init__(self, R, S=None, category=None):
        """
        Construct the affine scheme with coordinate ring `R`.

        INPUT:

        - ``R`` -- commutative ring

        - ``S`` -- (optional) commutative ring admitting a natural map
          to ``R``

        OUTPUT:

        The spectrum of `R`, i.e. the unique affine scheme with
        coordinate ring `R` as a scheme over the base ring `S`.

        EXAMPLES::

            sage: from sage.schemes.generic.scheme import AffineScheme
            sage: A.<x, y> = PolynomialRing(QQ)
            sage: X = AffineScheme(A, QQ)
            sage: X
            Spectrum of Multivariate Polynomial Ring in x, y over Rational Field
            sage: X.category()
            Category of schemes over Rational Field

        The standard way to construct an affine scheme is to use the
        :func:`~sage.schemes.generic.spec.Spec` functor::

            sage: S = Spec(ZZ)
            sage: S
            Spectrum of Integer Ring
            sage: S.category()
            Category of schemes
            sage: type(S)
            <class 'sage.schemes.generic.scheme.AffineScheme_with_category'>
        """
        if R not in CommutativeRings():
            raise TypeError("R (={}) must be a commutative ring".format(R))
        self.__R = R
        if S is not None:
            if S not in CommutativeRings():
                raise TypeError("S (={}) must be a commutative ring".format(S))
            if not R.has_coerce_map_from(S):
                raise ValueError("There must be a natural map S --> R, but S = {} and R = {}".format(S, R))
        Scheme.__init__(self, S, category=category)

    def __setstate__(self, state):
        """
        Needed to unpickle old Spec objects.

        The name-mangled attribute ``__R`` used to be in a class
        called ``Spec``; we have to translate this mangled name.

        TESTS::

            sage: S = Spec(QQ)
            sage: loads(dumps(S))
            Spectrum of Rational Field
        """
        if '_Spec__R' in state:
            state['_AffineScheme__R'] = state.pop('_Spec__R')
        super().__setstate__(state)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: Spec(PolynomialRing(QQ, 3, 'x'))
            Spectrum of Multivariate Polynomial Ring in x0, x1, x2 over Rational Field

        TESTS::

            sage: Spec(PolynomialRing(QQ, 3, 'x'))._repr_()
            'Spectrum of Multivariate Polynomial Ring in x0, x1, x2 over Rational Field'
        """
        return "Spectrum of {}".format(self.__R)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: S = Spec(PolynomialRing(ZZ, 2, 'x'))
            sage: S
            Spectrum of Multivariate Polynomial Ring in x0, x1 over Integer Ring
            sage: S._latex_()
            '\\mathrm{Spec}(\\Bold{Z}[x_{0}, x_{1}])'
        """
        return "\\mathrm{{Spec}}({})".format(self.__R._latex_())

    def __call__(self, *args):
        """
        Construct a scheme-valued or topological point of ``self``.

        INPUT/OUTPUT: the argument ``x`` must be one of the following:

        - a prime ideal of the coordinate ring; the output will
          be the corresponding point of `X`

        - a ring or a scheme `S`; the output will be the set `X(S)` of
          `S`-valued points on `X`

        EXAMPLES::

            sage: S = Spec(ZZ)
            sage: P = S(ZZ.ideal(3)); P
            Point on Spectrum of Integer Ring defined by the Principal ideal (3) of Integer Ring
            sage: type(P)
            <class 'sage.schemes.generic.scheme.AffineScheme_with_category.element_class'>
            sage: S(ZZ.ideal(next_prime(1000000)))                                      # needs sage.libs.pari
            Point on Spectrum of Integer Ring
             defined by the Principal ideal (1000003) of Integer Ring

            sage: R.<x, y, z> = QQ[]
            sage: S = Spec(R)
            sage: P = S(R.ideal(x, y, z)); P
            Point on Spectrum of Multivariate Polynomial Ring
             in x, y, z over Rational Field defined by the Ideal (x, y, z)
              of Multivariate Polynomial Ring in x, y, z over Rational Field

        This indicates the fix of :issue:`12734`::

            sage: S = Spec(ZZ)
            sage: S(ZZ)
            Set of rational points of Spectrum of Integer Ring

        Note the difference between the previous example and the
        following one::

            sage: S(S)
            Set of morphisms
              From: Spectrum of Integer Ring
              To:   Spectrum of Integer Ring

        For affine or projective varieties, passing the correct number
        of elements of the base ring constructs the rational point
        with these elements as coordinates::

            sage: S = AffineSpace(ZZ, 1)
            sage: S(0)
            (0)

        To prevent confusion with this usage, topological points must
        be constructed by explicitly specifying a prime ideal, not
        just generators::

            sage: R = S.coordinate_ring()
            sage: S(R.ideal(0))
            Point on Affine Space of dimension 1 over Integer Ring
             defined by the Ideal (0) of Multivariate Polynomial Ring in x over Integer Ring

        This explains why the following example raises an error rather
        than constructing the topological point defined by the prime
        ideal `(0)` as one might expect::

            sage: S = Spec(ZZ)
            sage: S(0)
            Traceback (most recent call last):
            ...
            TypeError: cannot call Spectrum of Integer Ring with arguments (0,)
        """
        if len(args) == 1:
            x = args[0]
            if ((isinstance(x, self.element_class) and (x.parent() is self or x.parent() == self))
                or (isinstance(x, Ideal_generic) and x.ring() is self.coordinate_ring())):
                # Construct a topological point from x.
                return self._element_constructor_(x)
        try:
            # Construct a scheme homset or a scheme-valued point from
            # args using the generic Scheme.__call__() method.
            return super().__call__(*args)
        except NotImplementedError:
            # This arises from self._morphism() not being implemented.
            # We must convert it into a TypeError to keep the coercion
            # system working.
            raise TypeError('cannot call %s with arguments %s' % (self, args))

    Element = SchemeTopologicalPoint_prime_ideal

    def _element_constructor_(self, x):
        """
        Construct a topological point from `x`.

        TESTS::

            sage: S = Spec(ZZ)
            sage: S(ZZ.ideal(0))
            Point on Spectrum of Integer Ring defined by the Principal ideal (0) of Integer Ring
        """
        if isinstance(x, self.element_class):
            if x.parent() is self:
                return x
            elif x.parent() == self:
                return self.element_class(self, x.prime_ideal())
        elif isinstance(x, Ideal_generic) and x.ring() is self.coordinate_ring():
            return self.element_class(self, x)
        raise TypeError('cannot convert %s to a topological point of %s' % (x, self))

    def _an_element_(self):
        r"""
        Return an element of the spectrum of the ring.

        OUTPUT: a point of the affine scheme ``self``

        EXAMPLES::

            sage: Spec(QQ).an_element()
            Point on Spectrum of Rational Field defined by the Principal ideal (0) of Rational Field
            sage: Spec(ZZ).an_element()    # random output
            Point on Spectrum of Integer Ring defined by the Principal ideal (811) of Integer Ring
        """
        if self.coordinate_ring() is ZZ:
            from sage.arith.misc import random_prime
            return self(ZZ.ideal(random_prime(1000)))
        return self(self.coordinate_ring().zero_ideal())

    def coordinate_ring(self):
        """
        Return the underlying ring of this scheme.

        OUTPUT: a commutative ring

        EXAMPLES::

            sage: Spec(QQ).coordinate_ring()
            Rational Field
            sage: Spec(PolynomialRing(QQ, 3, 'x')).coordinate_ring()
            Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        """
        return self.__R

    def is_noetherian(self):
        """
        Return ``True`` if ``self`` is Noetherian, ``False`` otherwise.

        EXAMPLES::

            sage: Spec(ZZ).is_noetherian()
            True
        """
        return self.__R.is_noetherian()

    def dimension_absolute(self):
        """
        Return the absolute dimension of this scheme.

        OUTPUT: integer

        EXAMPLES::

            sage: S = Spec(ZZ)
            sage: S.dimension_absolute()
            1
            sage: S.dimension()
            1
        """
        return self.__R.krull_dimension()

    dimension = dimension_absolute

    def dimension_relative(self):
        """
        Return the relative dimension of this scheme over its base.

        OUTPUT: integer

        EXAMPLES::

            sage: S = Spec(ZZ)
            sage: S.dimension_relative()
            0
        """
        return self.__R.krull_dimension() - self.base_ring().krull_dimension()

    def base_extend(self, R):
        """
        Extend the base ring/scheme.

        INPUT:

        - ``R`` -- an affine scheme or a commutative ring

        EXAMPLES::

            sage: Spec_ZZ = Spec(ZZ);  Spec_ZZ
            Spectrum of Integer Ring
            sage: Spec_ZZ.base_extend(QQ)
            Spectrum of Rational Field

            sage: Spec(ZZ['x']).base_extend(Spec(QQ))
            Spectrum of Univariate Polynomial Ring in x over Rational Field
        """
        if R in CommutativeRings():
            return AffineScheme(self.coordinate_ring().base_extend(R), self.base_ring())
        if not self.base_scheme() == R.base_scheme():
            raise ValueError('the new base scheme must be a scheme over the old base scheme')
        return AffineScheme(self.coordinate_ring().base_extend(R.coordinate_ring()),
                            self.base_ring())

    def _point_homset(self, *args, **kwds):
        """
        Construct a point Hom-set.

        For internal use only. See :mod:`morphism` for more details.

        EXAMPLES::

            sage: Spec(QQ)._point_homset(Spec(QQ), Spec(ZZ))
            Set of rational points of Spectrum of Integer Ring
        """
        from sage.schemes.affine.affine_homset import SchemeHomset_points_spec
        return SchemeHomset_points_spec(*args, **kwds)

    def hom(self, x, Y=None):
        r"""
        Return the scheme morphism from ``self`` to ``Y`` defined by ``x``.

        INPUT:

        - ``x`` -- anything that determines a scheme morphism; if
          ``x`` is a scheme, try to determine a natural map to ``x``

        - ``Y`` -- the codomain scheme (optional); if ``Y`` is not
          given, try to determine ``Y`` from context

        - ``check`` -- boolean (default: ``True``); whether
          to check the defining data for consistency

        OUTPUT: the scheme morphism from ``self`` to ``Y`` defined by ``x``

        EXAMPLES:

        We construct the inclusion from `\mathrm{Spec}(\QQ)` into
        `\mathrm{Spec}(\ZZ)` induced by the inclusion from `\ZZ` into
        `\QQ`::

            sage: X = Spec(QQ)
            sage: X.hom(ZZ.hom(QQ))
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Natural morphism:
                      From: Integer Ring
                      To:   Rational Field

        TESTS:

        We can construct a morphism to an affine curve (:issue:`7956`)::

            sage: S.<p,q> = QQ[]
            sage: A1.<r> = AffineSpace(QQ, 1)
            sage: A1_emb = Curve(p - 2)                                                 # needs sage.schemes
            sage: A1.hom([2, r], A1_emb)                                                # needs sage.schemes
            Scheme morphism:
              From: Affine Space of dimension 1 over Rational Field
              To:   Affine Plane Curve over Rational Field defined by p - 2
              Defn: Defined on coordinates by sending (r) to (2, r)
        """
        from sage.categories.map import Map
        from sage.categories.rings import Rings

        if isinstance(x, Scheme):
            return self.Hom(x).natural_map()
        if Y is None and isinstance(x, Map) and x.category_for().is_subcategory(Rings()):
            # x is a morphism of Rings
            Y = AffineScheme(x.domain())
        return Scheme.hom(self, x, Y)
