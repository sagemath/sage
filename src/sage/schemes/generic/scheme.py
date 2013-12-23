"""
Schemes

AUTHORS:

- William Stein, David Kohel, Kiran Kedlaya (2008): added zeta_series

- Volker Braun (2011-08-11): documenting, improving, refactoring.
"""


#*****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2008 Kiran Kedlaya <kedlaya@mit.edu>
#       Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2005 William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.structure.parent import Parent
from sage.misc.all import cached_method
from sage.rings.all import (IntegerRing,
                            ZZ, GF, PowerSeriesRing,
                            Rationals)

from sage.rings.commutative_ring import is_CommutativeRing
from sage.rings.morphism import is_RingHomomorphism

def is_Scheme(x):
    """
    Test whether ``x`` is a scheme.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean. Whether ``x`` derives from :class:`Scheme`.

    EXAMPLES::

        sage: from sage.schemes.generic.scheme import is_Scheme
        sage: is_Scheme(5)
        False
        sage: X = Spec(QQ)
        sage: is_Scheme(X)
        True
    """
    return isinstance(x, Scheme)



class Scheme(Parent):
    """
    The base class for all schemes.

    INPUT:

    - ``X`` -- a scheme, scheme morphism, commutative ring,
      commutative ring morphism, or ``None`` (optional). Determines
      the base scheme. If a commutative ring is passed, the spectrum
      of the ring will be used as base.

    - ``category`` -- the category (optional). Will be automatically
      construted by default.

    EXAMPLES::

        sage: from sage.schemes.generic.scheme import Scheme
        sage: Scheme(ZZ)
        <class 'sage.schemes.generic.scheme.Scheme_with_category'>

    A scheme is in the category of all schemes over its base::

        sage: ProjectiveSpace(4, QQ).category()
        Category of schemes over Rational Field

    There is a special and unique `Spec(\ZZ)` that is the default base
    scheme::

        sage: Spec(ZZ).base_scheme() is Spec(QQ).base_scheme()
        True
    """

    def __init__(self, X=None, category=None):
        """
        Construct a scheme.

        TESTS::

            sage: R.<x, y> = QQ[]
            sage: I = (x^2 - y^2)*R
            sage: RmodI = R.quotient(I)
            sage: X = Spec(RmodI)
            sage: TestSuite(X).run(skip = ["_test_an_element", "_test_elements",
            ...                            "_test_some_elements", "_test_category"]) # See #7946
        """
        from sage.schemes.generic.spec import is_Spec
        from sage.schemes.generic.morphism import is_SchemeMorphism

        if X is None:
            try:
                from sage.schemes.generic.spec import SpecZ
                self._base_scheme = SpecZ
            except ImportError:  # we are currently constructing SpecZ
                self._base_ring = ZZ
        elif is_Scheme(X):
            self._base_scheme = X
        elif is_SchemeMorphism(X):
            self._base_morphism = X
        elif is_CommutativeRing(X):
            self._base_ring = X
        elif is_RingHomomorphism(X):
            self._base_ring = X.codomain()
        else:
            raise ValueError('The base must be define by a scheme, '
                             'scheme morphism, or commutative ring.')

        from sage.categories.schemes import Schemes
        if not X:
            default_category = Schemes()
        else:
            default_category = Schemes(self.base_scheme())
        if category is None:
            category = default_category
        else:
            assert category.is_subcategory(default_category), \
                "%s is not a subcategory of %s"%(category, default_category)

        Parent.__init__(self, self.base_ring(), category = category)

    def __cmp__(left, right):
        """
        Compare two schemes.

        INPUT:

        - ``right`` -- anything. To compare against the scheme
          ``left``.

        OUTPUT:

        ``+1``, ``0``, or ``-1``.

        EXAMPLES::

            sage: X = Spec(QQ);  Y = Spec(QQ)
            sage: X == Y
            True
            sage: X is Y
            False
        """
        if not is_Scheme(right):
            return -1
        return left._cmp_(right)

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

        This shows that issue at trac ticket 7389 is solved::

            sage: S = Spec(ZZ)
            sage: f = S.identity_morphism()
            sage: from sage.schemes.generic.glue import GluedScheme
            sage: T = GluedScheme(f,f)
            sage: S.hom([1],T)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def base_extend(self, Y):
        """
        Extend the base of the scheme.

        Derived clases must override this method.

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

        INPUT/OUTPUT:

        The arguments must be one of the following:

        - a ring or a scheme `S`. Output will be the set `X(S)` of
          `S`-valued points on `X`.

        - If `S` is a list or tuple or just the coordinates, return a
          point in `X(T)`, where `T` is the base scheme of self.

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)

        We create some point sets::

            sage: A(QQ)
            Set of rational points of Affine Space of dimension 2 over Rational Field
            sage: A(RR)
            Set of rational points of Affine Space of dimension 2 over Real Field
            with 53 bits of precision

        Space of dimension 2 over Rational Field::

            sage: R.<x> = PolynomialRing(QQ)
            sage: A(NumberField(x^2+1, 'a'))
            Set of rational points of Affine Space of dimension 2 over Number Field
            in a with defining polynomial x^2 + 1
            sage: A(GF(7))
            Traceback (most recent call last):
            ...
            ValueError: There must be a natural map S --> R, but
            S = Rational Field and R = Finite Field of size 7

        We create some points::

            sage: A(QQ)([1,0])
            (1, 0)

        We create the same point by giving the coordinates of the point
        directly::

            sage: A( 1,0 )
            (1, 0)
        """
        if len(args) == 0:
            raise TypeError('You need to specify at least one argument.')

        S = args[0]
        if is_CommutativeRing(S):
            return self.point_homset(S)
        if is_Scheme(S):
            return S.Hom(self)
        from sage.schemes.generic.morphism import SchemeMorphism_point
        if isinstance(S, (list, tuple)):
            args = S
        elif isinstance(S, SchemeMorphism_point):
            if S.codomain() == self:
                return S
        else:
            # TODO: fix circular import resulting from non-multiple inheritance
            from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
            if isinstance(S, EllipticCurvePoint_field):
                if S.codomain() == self:
                    return S
                else:
                    return self.point(S)
        return self.point(args)

    @cached_method
    def point_homset(self, S=None):
        """
        Return the set of S-valued points of this scheme.

        INPUT:

        - ``S`` -- a commutative ring.

        OUTPUT:

        The set of morphisms `Spec(S)\to X`.

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

            sage: P = ProjectiveSpace(QQ,3)
            sage: P.point_homset(GF(11))
            Traceback (most recent call last):
            ...
            ValueError: There must be a natural map S --> R, but
            S = Rational Field and R = Finite Field of size 11
        """
        if S is None:
            S = self.base_ring()
        from sage.schemes.generic.spec import Spec
        SpecS = Spec(S, self.base_ring())
        from sage.schemes.generic.homset import SchemeHomset
        return SchemeHomset(SpecS, self)

    def point(self, v, check=True):
        """
        Create a point.

        INPUT:

        - ``v`` -- anything that defines a point.

        - ``check`` -- boolean (optional, default=``True``). Whether
          to check the defining data for consistency.

        OUTPUT:

        A point of the scheme.

        EXAMPLES::

            sage: A2 = AffineSpace(QQ,2)
            sage: A2.point([4,5])
            (4, 5)

            sage: R.<t> = PolynomialRing(QQ)
            sage: E = EllipticCurve([t + 1, t, t, 0, 0])
            sage: E.point([0, 0])
            (0 : 0 : 1)
        """
        # todo: update elliptic curve stuff to take point_homset as argument
        from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
        if is_EllipticCurve(self):
            try:
                return self._point(self.point_homset(), v, check=check)
            except AttributeError:  # legacy code without point_homset
                return self._point(self, v, check=check)

        return self.point_homset() (v, check=check)

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

    def __div__(self, Y):
        """
        Return the base extension of self to Y.

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
        Return the base ring of the scheme self.

        OUTPUT:

        A commutative ring.

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

        OUTPUT:

        A scheme.

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
                from sage.schemes.generic.spec import Spec
                self._base_scheme = Spec(self._base_ring)
            else:
                from sage.schemes.generic.spec import SpecZ
                self._base_scheme = SpecZ
            return self._base_scheme

    def base_morphism(self):
        """
        Return the structure morphism from ``self`` to its base
        scheme.

        OUTPUT:

        A scheme morphism.

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
            from sage.schemes.generic.spec import Spec, SpecZ
            SCH = Schemes()
            if hasattr(self, '_base_scheme'):
                self._base_morphism = self.Hom(self._base_scheme, category=SCH).natural_map()
            elif hasattr(self, '_base_ring'):
                self._base_morphism = self.Hom(Spec(self._base_ring), category=SCH).natural_map()
            else:
                self._base_morphism = self.Hom(SpecZ, category=SCH).natural_map()
            return self._base_morphism

    structure_morphism = base_morphism

    def coordinate_ring(self):
        """
        Return the coordinate ring.

        OUTPUT:

        The global coordinate ring of this scheme, if
        defined. Otherwise raise a ``ValueError``.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: I = (x^2 - y^2)*R
            sage: X = Spec(R.quotient(I))
            sage: X.coordinate_ring()
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 - y^2)
        """
        try:
            return self._coordinate_ring
        except AttributeError:
            raise ValueError, "This scheme has no associated coordinated ring (defined)."

    def dimension_absolute(self):
        """
        Return the absolute dimension of this scheme.

        OUTPUT:

        Integer.

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

        OUTPUT:

        Integer.

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

        OUTPUT:

        The identity morphism of the scheme ``self``.

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

        - ``x`` -- anything hat determines a scheme morphism. If ``x``
          is a scheme, try to determine a natural map to ``x``.

        - ``Y`` -- the codomain scheme (optional). If ``Y`` is not
          given, try to determine ``Y`` from context.

        - ``check`` -- boolean (optional, default=``True``). Whether
          to check the defining data for consistency.

        OUTPUT:

        The scheme morphism from ``self`` to ``Y`` defined by ``x``.

        EXAMPLES::

            sage: P = ProjectiveSpace(ZZ, 3)
            sage: P.hom(Spec(ZZ))
            Scheme morphism:
              From: Projective Space of dimension 3 over Integer Ring
              To:   Spectrum of Integer Ring
              Defn: Structure map
        """
        if Y is None:
            if is_Scheme(x):
                return self.Hom(x).natural_map()
            else:
                raise TypeError, "unable to determine codomain"
        return self.Hom(Y)(x, check)

    def _Hom_(self, Y, category=None, check=True):
        """
        Return the set of scheme morphisms from ``self`` to ``Y``.

        INPUT:

        - ``Y`` -- a scheme. The codomain of the Hom-set.

        - ``category`` -- a category (optional). The category of the
          Hom-set.

        - ``check`` -- boolean (optional, default=``True``). Whether
          to check the defining data for consistency.

        OUTPUT:

        The set of morphisms from ``self`` to ``Y``.

        EXAMPLES::

            sage: P = ProjectiveSpace(ZZ, 3)
            sage: S = Spec(ZZ)
            sage: S._Hom_(P)
            Set of rational points of Projective Space of dimension 3 over Integer Ring

        TESTS::

            sage: S._Hom_(P).__class__
            <class 'sage.schemes.projective.projective_homset.SchemeHomset_points_projective_ring_with_category'>

            sage: E = EllipticCurve('37a1')
            sage: Hom(E, E).__class__
            <class 'sage.schemes.generic.homset.SchemeHomset_generic_with_category'>

            sage: Hom(Spec(ZZ), Spec(ZZ)).__class__
            <class 'sage.schemes.affine.affine_homset.SchemeHomset_points_spec_with_category'>
        """
        from sage.schemes.generic.homset import SchemeHomset
        return SchemeHomset(self, Y, category=category, check=check)

    point_set = point_homset

    def count_points(self, n):
        r"""
        Count points over finite fields.

        INPUT:

        - ``n`` -- integer.

        OUTPUT:

        An integer. The number of points over `\GF{q}, \ldots,
        \GF{q^n}` on a scheme over a finite field `\GF{q}`.

        .. note::

           This is currently only implemented for schemes over prime
           order finite fields.

        EXAMPLES::

            sage: P.<x> = PolynomialRing(GF(3))
            sage: C = HyperellipticCurve(x^3+x^2+1)
            sage: C.count_points(4)
            [6, 12, 18, 96]
            sage: C.base_extend(GF(9,'a')).count_points(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Point counting only implemented for schemes over prime fields
        """
        F = self.base_ring()
        if not F.is_finite():
            raise TypeError, "Point counting only defined for schemes over finite fields"
        q = F.cardinality()
        if not q.is_prime():
            raise NotImplementedError, "Point counting only implemented for schemes over prime fields"
        a = []
        for i in range(1, n+1):
            F1 = GF(q**i, name='z')
            S1 = self.base_extend(F1)
            a.append(len(S1.rational_points()))
        return(a)

    def zeta_series(self, n, t):
        """
        Return the zeta series.

        Compute a power series approximation to the zeta function of a
        scheme over a finite field.

        INPUT:

        -  ``n`` -- the number of terms of the power series to
           compute

        -  ``t`` -- the variable which the series should be
           returned


        OUTPUT:

        A power series approximating the zeta function of self

        EXAMPLES::

            sage: P.<x> = PolynomialRing(GF(3))
            sage: C = HyperellipticCurve(x^3+x^2+1)
            sage: R.<t> = PowerSeriesRing(Integers())
            sage: C.zeta_series(4,t)
            1 + 6*t + 24*t^2 + 78*t^3 + 240*t^4 + O(t^5)
            sage: (1+2*t+3*t^2)/(1-t)/(1-3*t) + O(t^5)
            1 + 6*t + 24*t^2 + 78*t^3 + 240*t^4 + O(t^5)

        Note that this function depends on count_points, which is only
        defined for prime order fields::

            sage: C.base_extend(GF(9,'a')).zeta_series(4,t)
            Traceback (most recent call last):
            ...
            NotImplementedError: Point counting only implemented for schemes over prime fields
        """

        F = self.base_ring()
        if not F.is_finite():
            raise TypeError('zeta functions only defined for schemes over finite fields')
        try:
            a = self.count_points(n)
        except AttributeError:
            raise NotImplementedError('count_points() required but not implemented')
        R = PowerSeriesRing(Rationals(), 'u')
        u = R.gen()
        temp = sum(a[i-1]*(u.O(n+1))**i/i for i in range(1,n+1))
        temp2 = temp.exp()
        return(temp2(t).O(n+1))

def is_AffineScheme(x):
    """
    Return True if `x` is an affine scheme.

    EXAMPLES::

        sage: from sage.schemes.generic.scheme import is_AffineScheme
        sage: is_AffineScheme(5)
        False
        sage: E = Spec(QQ)
        sage: is_AffineScheme(E)
        True
    """
    return isinstance(x, AffineScheme)

class AffineScheme(Scheme):
    """
    An abstract affine scheme.
    """
    def hom(self, x, Y=None):
        r"""
        Return the scheme morphism from ``self`` to ``Y`` defined by ``x``.

        INPUT:

        - ``x`` -- anything hat determines a scheme morphism. If ``x``
          is a scheme, try to determine a natural map to ``x``.

        - ``Y`` -- the codomain scheme (optional). If ``Y`` is not
          given, try to determine ``Y`` from context.

        - ``check`` -- boolean (optional, default=``True``). Whether
          to check the defining data for consistency.

        OUTPUT:

        The scheme morphism from ``self`` to ``Y`` defined by ``x``.

        EXAMPLES:

        We construct the inclusion from `\mathrm{Spec}(\QQ)` into
        `\mathrm{Spec}(\ZZ)` induced by the inclusion from `\ZZ` into
        `\QQ`::

            sage: X = Spec(QQ)
            sage: X.hom(ZZ.hom(QQ))
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Ring Coercion morphism:
                      From: Integer Ring
                      To:   Rational Field

        TESTS:

        We can construct a morphism to an affine curve (trac #7956)::

            sage: S.<p,q> = QQ[]
            sage: A1.<r> = AffineSpace(QQ,1)
            sage: A1_emb = Curve(p-2)
            sage: A1.hom([2,r],A1_emb)
            Scheme morphism:
              From: Affine Space of dimension 1 over Rational Field
              To:   Affine Curve over Rational Field defined by p - 2
              Defn: Defined on coordinates by sending (r) to
                    (2, r)

        """
        if is_Scheme(x):
            return self.Hom(x).natural_map()
        if Y is None:
            if is_RingHomomorphism(x):
                import spec
                Y = spec.Spec(x.domain())
        return Scheme.hom(self, x, Y)

