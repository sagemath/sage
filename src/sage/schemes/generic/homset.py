r"""
Set of homomorphisms between two schemes

For schemes `X` and `Y`, this module implements the set of morphisms
`Hom(X,Y)`. This is done by :class:`SchemeHomset_generic`.

As a special case, the Hom-sets can also represent the points of a
scheme. Recall that the `K`-rational points of a scheme `X` over `k`
can be identified with the set of morphisms `Spec(K) \to X`. In Sage
the rational points are implemented by such scheme morphisms. This is
done by :class:`SchemeHomset_points` and its subclasses.

.. note::

    You should not create the Hom-sets manually. Instead, use the
    :meth:`~sage.structure.parent.Hom` method that is inherited by all
    schemes.

AUTHORS:

- William Stein (2006): initial version.

- Volker Braun (2011-08-11): significant improvement and refactoring.

- Ben Hutz (June 2012): added support for projective ring
"""


#*****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.categories.homset import HomsetWithBase
from sage.structure.factory import UniqueFactory

from sage.rings.all import ( gcd, ZZ, QQ )

from sage.rings.morphism import is_RingHomomorphism
from sage.rings.rational_field import is_RationalField
from sage.rings.finite_rings.constructor import is_FiniteField
from sage.rings.commutative_ring import is_CommutativeRing


from sage.schemes.generic.scheme import is_Scheme
from sage.schemes.generic.spec import Spec, is_Spec
from sage.schemes.generic.morphism import (
    SchemeMorphism,
    SchemeMorphism_structure_map,
    SchemeMorphism_spec )


def is_SchemeHomset(H):
    r"""
    Test whether ``H`` is a scheme Hom-set.

    EXAMPLES::

        sage: f = Spec(QQ).identity_morphism();  f
        Scheme endomorphism of Spectrum of Rational Field
          Defn: Identity map
        sage: from sage.schemes.generic.homset import is_SchemeHomset
        sage: is_SchemeHomset(f)
        False
        sage: is_SchemeHomset(f.parent())
        True
        sage: is_SchemeHomset('a string')
        False
    """
    return isinstance(H, SchemeHomset_generic)


#*******************************************************************
# Factory for Hom sets of schemes
#*******************************************************************
class SchemeHomsetFactory(UniqueFactory):
    """
    Factory for Hom-sets of schemes.

    EXAMPLES::

        sage: A2 = AffineSpace(QQ,2)
        sage: A3 = AffineSpace(QQ,3)
        sage: Hom = A3.Hom(A2)

    The Hom-sets are unique::

        sage: Hom is copy(Hom)
        True
        sage: Hom is A3.Hom(A2)
        True
        sage: loads(Hom.dumps()) is Hom
        True

    Here is a tricky point. The Hom-sets are not identical if
    domains/codomains are isomorphic but not identiacal::

        sage: A3_iso = AffineSpace(QQ,3)
        sage: [ A3_iso is A3, A3_iso == A3 ]
        [False, True]
        sage: Hom_iso = A3_iso.Hom(A2)
        sage: Hom_iso is Hom
        False
        sage: Hom_iso == Hom
        True

    TESTS::

        sage: Hom.base()
        Integer Ring
        sage: Hom.base_ring()
        Integer Ring
    """

    def create_key_and_extra_args(self, X, Y, category=None, base=ZZ,
                                  check=True):
        """
        Create a key that uniquely determines the Hom-set.

        INPUT:

        - ``X`` -- a scheme. The domain of the morphisms.

        - ``Y`` -- a scheme. The codomain of the morphisms.

        - ``category`` -- a category for the Hom-sets (default: schemes over
          given base).

        - ``base`` -- a scheme or a ring. The base scheme of domain
          and codomain schemes. If a ring is specified, the spectrum
          of that ring will be used as base scheme.

        - ``check`` -- boolean (default: ``True``).

        EXAMPLES::

            sage: A2 = AffineSpace(QQ,2)
            sage: A3 = AffineSpace(QQ,3)
            sage: A3.Hom(A2)    # indirect doctest
            Set of morphisms
              From: Affine Space of dimension 3 over Rational Field
              To:   Affine Space of dimension 2 over Rational Field
            sage: from sage.schemes.generic.homset import SchemeHomsetFactory
            sage: SHOMfactory = SchemeHomsetFactory('test')
            sage: key, extra = SHOMfactory.create_key_and_extra_args(A3,A2,check=False)
            sage: key
            (..., ..., Category of schemes over Integer Ring)
            sage: extra
            {'Y': Affine Space of dimension 2 over Rational Field,
             'X': Affine Space of dimension 3 over Rational Field,
             'base_ring': Integer Ring, 'check': False}
        """
        if not is_Scheme(X) and is_CommutativeRing(X):
            X = Spec(X)
        if not is_Scheme(Y) and is_CommutativeRing(Y):
            Y = Spec(Y)
        if is_Spec(base):
            base_spec = base
            base_ring = base.coordinate_ring()
        elif is_CommutativeRing(base):
            base_spec = Spec(base)
            base_ring = base
        else:
            raise ValueError(
                        'The base must be a commutative ring or its spectrum.')
        if not category:
            from sage.categories.schemes import Schemes
            category = Schemes(base_spec)
        key = tuple([id(X), id(Y), category])
        extra = {'X':X, 'Y':Y, 'base_ring':base_ring, 'check':check}
        return key, extra

    def create_object(self, version, key, **extra_args):
        """
        Create a :class:`SchemeHomset_generic`.

        INPUT:

        - ``version`` -- object version. Currently not used.

        - ``key`` -- a key created by :meth:`create_key_and_extra_args`.

        - ``extra_args`` -- a dictionary of extra keyword arguments.

        EXAMPLES::

            sage: A2 = AffineSpace(QQ,2)
            sage: A3 = AffineSpace(QQ,3)
            sage: A3.Hom(A2) is A3.Hom(A2)   # indirect doctest
            True
            sage: from sage.schemes.generic.homset import SchemeHomsetFactory
            sage: SHOMfactory = SchemeHomsetFactory('test')
            sage: SHOMfactory.create_object(0, [id(A3),id(A2),A3.category()], check=True,
            ...                             X=A3, Y=A2, base_ring=QQ)
            Set of morphisms
              From: Affine Space of dimension 3 over Rational Field
              To:   Affine Space of dimension 2 over Rational Field
        """
        category = key[2]
        X = extra_args.pop('X')
        Y = extra_args.pop('Y')
        base_ring = extra_args.pop('base_ring')
        if is_Spec(X):
            return Y._point_homset(X, Y, category=category, base=base_ring, **extra_args)
        try:
            return X._homset(X, Y, category=category, base=base_ring, **extra_args)
        except AttributeError:
            return SchemeHomset_generic(X, Y, category=category, base=base_ring, **extra_args)


SchemeHomset = SchemeHomsetFactory('sage.schemes.generic.homset.SchemeHomset')



#*******************************************************************
# Base class
#*******************************************************************
class SchemeHomset_generic(HomsetWithBase):
    r"""
    The base class for Hom-sets of schemes.

    INPUT:

    - ``X`` -- a scheme. The domain of the Hom-set.

    - ``Y`` -- a scheme. The codomain of the Hom-set.

    - ``category`` -- a category (optional). The category of the
      Hom-set.

    - ``check`` -- boolean (optional, default=``True``). Whether to
      check the defining data for consistency.

    EXAMPLES::

        sage: from sage.schemes.generic.homset import SchemeHomset_generic
        sage: A2 = AffineSpace(QQ,2)
        sage: Hom = SchemeHomset_generic(A2, A2); Hom
        Set of morphisms
          From: Affine Space of dimension 2 over Rational Field
          To:   Affine Space of dimension 2 over Rational Field
        sage: Hom.category()
        Category of hom sets in Category of Schemes
    """
    Element = SchemeMorphism

    def __call__(self, *args, **kwds):
        r"""
        Make Hom-sets callable.

        See the ``_call_()`` method of the derived class. All
        arguments are handed through.

        EXAMPLES::

            sage: A2 = AffineSpace(QQ,2)
            sage: A2(4,5)
            (4, 5)
        """
        # Homset (base of HomsetWithBase) overrides __call__ @#$
        from sage.structure.parent import Set_generic
        return Set_generic.__call__(self, *args, **kwds)

    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: A = AffineSpace(4, QQ)
            sage: print A.structure_morphism()._repr_()
            Scheme morphism:
              From: Affine Space of dimension 4 over Rational Field
              To:   Spectrum of Rational Field
              Defn: Structure map
        """
        s = 'Set of morphisms'
        s += '\n  From: %s' % self.domain()
        s += '\n  To:   %s' % self.codomain()
        return s

    def natural_map(self):
        r"""
        Return a natural map in the Hom space.

        OUTPUT:

        A :class:`SchemeMorphism` if there is a natural map from
        domain to codomain. Otherwise, a ``NotImplementedError`` is
        raised.

        EXAMPLES::

            sage: A = AffineSpace(4, QQ)
            sage: A.structure_morphism()   # indirect doctest
            Scheme morphism:
              From: Affine Space of dimension 4 over Rational Field
              To:   Spectrum of Rational Field
              Defn: Structure map
        """
        X = self.domain()
        Y = self.codomain()
        if is_Spec(Y) and Y.coordinate_ring() == X.base_ring():
            return SchemeMorphism_structure_map(self)
        raise NotImplementedError

    def _element_constructor_(self, x, check=True):
        """
        Construct a scheme morphism.

        INPUT:

        - `x` -- a ring morphism, or a list or a tuple that define a
          ring morphism.

        - ``check`` -- boolean (default: ``True``) passed onto
          functions called by this one to be more careful about input
          argument type checking.

        EXAMPLES::

            sage: f = ZZ.hom(QQ); f
            Ring Coercion morphism:
              From: Integer Ring
              To:   Rational Field

            sage: H = Hom(Spec(QQ, ZZ), Spec(ZZ)); H
            Set of rational points of Spectrum of Integer Ring

            sage: phi = H(f); phi
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Ring Coercion morphism:
                      From: Integer Ring
                      To:   Rational Field

        TESTS:

            sage: H._element_constructor_(f)
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Ring Coercion morphism:
                      From: Integer Ring
                      To:   Rational Field

        We illustrate input type checking::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: C = A.subscheme(x*y-1)
            sage: H = C.Hom(C); H
            Set of morphisms
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x*y - 1
              To:   Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x*y - 1
            sage: H(1)
            Traceback (most recent call last):
            ...
            TypeError: x must be a ring homomorphism, list or tuple
        """
        if isinstance(x, (list, tuple)):
            return self.domain()._morphism(self, x, check=check)

        if is_RingHomomorphism(x):
            return SchemeMorphism_spec(self, x, check=check)

        raise TypeError, "x must be a ring homomorphism, list or tuple"


#*******************************************************************
# Base class for points
#*******************************************************************
class SchemeHomset_points(SchemeHomset_generic):
    """
    Set of rational points of the scheme.

    Recall that the `K`-rational points of a scheme `X` over `k` can
    be identified with the set of morphisms `Spec(K) \to X`. In Sage,
    the rational points are implemented by such scheme morphisms.

    If a scheme has a finite number of points, then the homset is
    supposed to implement the Python iterator interface. See
    :class:`~sage.schemes.toric.homset.SchemeHomset_points_toric_field`
    for example.

    INPUT:

    See :class:`SchemeHomset_generic`.

    EXAMPLES::

        sage: from sage.schemes.generic.homset import SchemeHomset_points
        sage: SchemeHomset_points(Spec(QQ), AffineSpace(ZZ,2))
        Set of rational points of Affine Space of dimension 2 over Rational Field
    """

    def __init__(self, X, Y, category=None, check=True, base=ZZ):
        """
        Python constructor.

        INPUT:

        See :class:`SchemeHomset_generic`.

        EXAMPLES::

            sage: from sage.schemes.generic.homset import SchemeHomset_points
            sage: SchemeHomset_points(Spec(QQ), AffineSpace(ZZ,2))
            Set of rational points of Affine Space of dimension 2 over Rational Field
        """
        if check and not is_Spec(X):
            raise ValueError('The domain must be an affine scheme.')
        SchemeHomset_generic.__init__(self, X, Y, category=category, check=check, base=base)

    def _element_constructor_(self, *v, **kwds):
        """
        The element contstructor.

        INPUT:

        - ``v`` -- anything that determines a scheme morphism in the
          Hom-set.

        OUTPUT:

        The scheme morphism determined by ``v``.

        EXAMPLES::

            sage: A2 = AffineSpace(ZZ,2)
            sage: F = GF(3)
            sage: F_points = A2(F);  type(F_points)
            <class 'sage.schemes.affine.affine_homset.SchemeHomset_points_affine_with_category'>
            sage: F_points([2,5])
            (2, 2)

            sage: P2 = ProjectiveSpace(GF(3),2)
            sage: F.<a> = GF(9,'a')
            sage: F_points = P2(F)
            sage: type(F_points)
            <class 'sage.schemes.projective.projective_homset.SchemeHomset_points_projective_field_with_category'>
            sage: F_points([4,2*a])
            (1 : 2*a : 1)

        TESTS::

            sage: F_points._element_constructor_([4,2*a])
            (1 : 2*a : 1)
        """
        if len(v) == 1:
            v = v[0]
        return self.codomain()._point(self, v, **kwds)

    def extended_codomain(self):
        """
        Return the codomain with extended base, if necessary.

        OUTPUT:

        The codomain scheme, with its base ring extended to the
        codomain. That is, the codomain is of the form `Spec(R)` and
        the base ring of the domain is extended to `R`.

        EXAMPLES::

            sage: P2 = ProjectiveSpace(QQ,2)
            sage: K.<a> = NumberField(x^2 + x - (3^3-3))
            sage: K_points = P2(K);  K_points
            Set of rational points of Projective Space of dimension 2
            over Number Field in a with defining polynomial x^2 + x - 24

            sage: K_points.codomain()
            Projective Space of dimension 2 over Rational Field

            sage: K_points.extended_codomain()
            Projective Space of dimension 2 over Number Field in a with
            defining polynomial x^2 + x - 24
        """
        if '_extended_codomain' in self.__dict__:
            return self._extended_codomain
        R = self.domain().coordinate_ring()
        if R is not self.codomain().base_ring():
            X = self.codomain().base_extend(R)
        else:
            X = self.codomain()
        self._extended_codomain = X
        return X

    def _repr_(self):
        """
        Return a string representation of ``self``.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: P2 = ProjectiveSpace(ZZ,2)
            sage: P2(QQ)._repr_()
            'Set of rational points of Projective Space of dimension 2 over Rational Field'
        """
        return 'Set of rational points of '+str(self.extended_codomain())

    def value_ring(self):
        """
        Return `R` for a point Hom-set `X(Spec(R))`.

        OUTPUT:

        A commutative ring.

        EXAMPLES::

            sage: P2 = ProjectiveSpace(ZZ,2)
            sage: P2(QQ).value_ring()
            Rational Field
        """
        dom = self.domain()
        if not is_Spec(dom):
            raise ValueError("value rings are defined for Spec domains only!")
        return dom.coordinate_ring()

    def cardinality(self):
        """
        Return the number of points.

        OUTPUT:

        An integer or infinity.

        EXAMPLES::

            sage: toric_varieties.P2().point_set().cardinality()
            +Infinity

            sage: P2 = toric_varieties.P2(base_ring=GF(3))
            sage: P2.point_set().cardinality()
            13
        """
        if hasattr(self, 'is_finite') and not self.is_finite():
            from sage.rings.infinity import Infinity
            return Infinity
        return sum(ZZ.one() for point in self)

    __len__ = cardinality

    def list(self):
        """
        Return a tuple containing all points.

        OUTPUT:

        A tuple containing all points of the toric variety.

        EXAMPLE::

            sage: P1 = toric_varieties.P1(base_ring=GF(3))
            sage: P1.point_set().list()
            ([0 : 1], [1 : 0], [1 : 1], [1 : 2])
        """
        return tuple(self)
