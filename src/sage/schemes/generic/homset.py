# sage_setup: distribution = sagemath-categories
r"""
Set of homomorphisms between two schemes

For schemes `X` and `Y`, this module implements the set of morphisms
`\mathrm{Hom}(X,Y)`. This is done by :class:`SchemeHomset_generic`.

As a special case, the Hom-sets can also represent the points of a scheme.
Recall that the `K`-rational points of a scheme `X` over `k` can be identified
with the set of morphisms `\mathrm{Spec}(K) \to X`. In Sage the rational points
are implemented by such scheme morphisms. This is done by
:class:`SchemeHomset_points` and its subclasses.

.. NOTE::

    You should not create the Hom-sets manually. Instead, use the
    :meth:`~sage.structure.parent.Hom` method that is inherited by all
    schemes.

AUTHORS:

- William Stein (2006): initial version.

- Volker Braun (2011-08-11): significant improvement and refactoring.

- Ben Hutz (June 2012): added support for projective ring
"""

# *****************************************************************************
#        Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#        Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.homset import HomsetWithBase
from sage.misc.lazy_import import lazy_import
from sage.structure.factory import UniqueFactory
from sage.structure.parent import Set_generic

from sage.rings.integer_ring import ZZ
from sage.rings.ring import CommutativeRing
from sage.categories.commutative_rings import CommutativeRings

from sage.schemes.generic.scheme import AffineScheme, is_AffineScheme
from sage.schemes.generic.morphism import (
    SchemeMorphism,
    SchemeMorphism_structure_map,
    SchemeMorphism_spec )

lazy_import('sage.schemes.affine.affine_space', 'AffineSpace_generic', as_='AffineSpace')
lazy_import('sage.schemes.generic.algebraic_scheme', 'AlgebraicScheme_subscheme')
lazy_import('sage.schemes.product_projective.space', 'ProductProjectiveSpaces_ring', as_='ProductProjectiveSpaces')
lazy_import('sage.schemes.projective.projective_space', 'ProjectiveSpace_ring', as_='ProjectiveSpace')


def is_SchemeHomset(H):
    r"""
    Test whether ``H`` is a scheme Hom-set.

    EXAMPLES::

        sage: f = Spec(QQ).identity_morphism();  f
        Scheme endomorphism of Spectrum of Rational Field
          Defn: Identity map
        sage: from sage.schemes.generic.homset import is_SchemeHomset
        sage: is_SchemeHomset(f)
        doctest:warning...
        DeprecationWarning: The function is_SchemeHomset is deprecated; use 'isinstance(..., SchemeHomset_generic)' instead.
        See https://github.com/sagemath/sage/issues/38022 for details.
        False
        sage: is_SchemeHomset(f.parent())
        True
        sage: is_SchemeHomset('a string')
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(38022, "The function is_SchemeHomset is deprecated; use 'isinstance(..., SchemeHomset_generic)' instead.")
    return isinstance(H, SchemeHomset_generic)


# *******************************************************************
#  Factory for Hom sets of schemes
# *******************************************************************

class SchemeHomsetFactory(UniqueFactory):
    """
    Factory for Hom-sets of schemes.

    EXAMPLES::

        sage: A2 = AffineSpace(QQ, 2)
        sage: A3 = AffineSpace(QQ, 3)
        sage: Hom = A3.Hom(A2)

    The Hom-sets are uniquely determined by domain and codomain::

        sage: Hom is copy(Hom)
        True
        sage: Hom is A3.Hom(A2)
        True

    The Hom-sets are identical if the domains and codomains are
    identical::

        sage: loads(Hom.dumps()) is Hom
        True
        sage: A3_iso = AffineSpace(QQ, 3)
        sage: A3_iso is A3
        True
        sage: Hom_iso = A3_iso.Hom(A2)
        sage: Hom_iso is Hom
        True

    TESTS::

        sage: Hom.base()
        Rational Field
        sage: Hom.base_ring()
        Rational Field
    """

    def create_key_and_extra_args(self, X, Y, category=None, base=None,
                                  check=True, as_point_homset=False):
        """
        Create a key that uniquely determines the Hom-set.

        INPUT:

        - ``X`` -- a scheme; the domain of the morphisms

        - ``Y`` -- a scheme; the codomain of the morphisms

        - ``category`` -- a category for the Hom-sets (default: schemes over
          given base)

        - ``base`` -- a scheme or a ring; the base scheme of domain
          and codomain schemes. If a ring is specified, the spectrum
          of that ring will be used as base scheme.

        - ``check`` -- boolean (default: ``True``)

        EXAMPLES::

            sage: A2 = AffineSpace(QQ, 2)
            sage: A3 = AffineSpace(QQ, 3)
            sage: A3.Hom(A2)    # indirect doctest
            Set of morphisms
              From: Affine Space of dimension 3 over Rational Field
              To:   Affine Space of dimension 2 over Rational Field
            sage: from sage.schemes.generic.homset import SchemeHomsetFactory
            sage: SHOMfactory = SchemeHomsetFactory('test')
            sage: key, extra = SHOMfactory.create_key_and_extra_args(A3, A2, check=False)
            sage: key
            (..., ..., Category of schemes over Rational Field, False)
            sage: extra
            {'X': Affine Space of dimension 3 over Rational Field,
             'Y': Affine Space of dimension 2 over Rational Field,
             'base_ring': Rational Field,
             'check': False}
        """
        _CommRings = CommutativeRings()
        if X in _CommRings:
            X = AffineScheme(X)
        if Y in _CommRings:
            Y = AffineScheme(Y)
        if base is None:
            from sage.structure.element import coercion_model
            base = coercion_model.common_parent(X.base_ring(), Y.base_ring())
        if isinstance(base, AffineScheme):
            base_spec = base
            base_ring = base.coordinate_ring()
        elif base in _CommRings:
            base_spec = AffineScheme(base)
            base_ring = base
        else:
            raise ValueError('base must be a commutative ring or its spectrum')
        if not category:
            from sage.categories.schemes import Schemes
            category = Schemes(base_spec)
        key = (id(X), id(Y), category, as_point_homset)
        extra = {'X':X, 'Y':Y, 'base_ring':base_ring, 'check':check}
        return key, extra

    def create_object(self, version, key, **extra_args):
        """
        Create a :class:`SchemeHomset_generic`.

        INPUT:

        - ``version`` -- object version; currently not used

        - ``key`` -- a key created by :meth:`create_key_and_extra_args`

        - ``extra_args`` -- dictionary of extra keyword arguments

        EXAMPLES::

            sage: A2 = AffineSpace(QQ, 2)
            sage: A3 = AffineSpace(QQ, 3)
            sage: A3.Hom(A2) is A3.Hom(A2)   # indirect doctest
            True
            sage: from sage.schemes.generic.homset import SchemeHomsetFactory
            sage: SHOMfactory = SchemeHomsetFactory('test')
            sage: SHOMfactory.create_object(0, [id(A3), id(A2), A3.category(), False],
            ....:                           check=True, X=A3, Y=A2, base_ring=QQ)
            Set of morphisms
              From: Affine Space of dimension 3 over Rational Field
              To:   Affine Space of dimension 2 over Rational Field
        """
        category = key[2]
        X = extra_args.pop('X')
        Y = extra_args.pop('Y')
        base_ring = extra_args.pop('base_ring')
        if len(key) >= 4 and key[3]:  # as_point_homset=True
            return Y._point_homset(X, Y, category=category, base=base_ring, **extra_args)
        try:
            return X._homset(X, Y, category=category, base=base_ring, **extra_args)
        except AttributeError:
            return SchemeHomset_generic(X, Y, category=category, base=base_ring, **extra_args)


SchemeHomset = SchemeHomsetFactory('sage.schemes.generic.homset.SchemeHomset')


# *******************************************************************
#  Base class
# *******************************************************************

class SchemeHomset_generic(HomsetWithBase):
    r"""
    The base class for Hom-sets of schemes.

    INPUT:

    - ``X`` -- a scheme; the domain of the Hom-set

    - ``Y`` -- a scheme; the codomain of the Hom-set

    - ``category`` -- a category (optional); the category of the
      Hom-set

    - ``check`` -- boolean (default: ``True``); whether to
      check the defining data for consistency

    EXAMPLES::

        sage: from sage.schemes.generic.homset import SchemeHomset_generic
        sage: A2 = AffineSpace(QQ,2)
        sage: Hom = SchemeHomset_generic(A2, A2); Hom
        Set of morphisms
          From: Affine Space of dimension 2 over Rational Field
          To:   Affine Space of dimension 2 over Rational Field
        sage: Hom.category()
        Category of endsets of schemes over Rational Field
    """
    Element = SchemeMorphism

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: A2 = AffineSpace(QQ, 2)
            sage: A3 = AffineSpace(QQ, 3)
            sage: Hom = A3.Hom(A2)
            sage: loads(Hom.dumps()) == Hom
            True
        """
        return SchemeHomset, (self.domain(), self.codomain(), self.homset_category(),
                              self.base_ring(), False, False)

    def __call__(self, *args, **kwds):
        r"""
        Make Hom-sets callable.

        See the ``_call_()`` method of the derived class. All
        arguments are handed through.

        EXAMPLES::

            sage: A2 = AffineSpace(QQ, 2)
            sage: A2(4,5)
            (4, 5)
        """
        # Homset (base of HomsetWithBase) overrides __call__ @#$
        return Set_generic.__call__(self, *args, **kwds)

    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: A = AffineSpace(4, QQ)
            sage: print(A.structure_morphism()._repr_())
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
        domain to codomain. Otherwise, a :exc:`NotImplementedError` is raised.

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
        if isinstance(Y, AffineScheme) and Y.coordinate_ring() == X.base_ring():
            return SchemeMorphism_structure_map(self)
        raise NotImplementedError

    def _element_constructor_(self, x, check=True):
        """
        Construct a scheme morphism.

        INPUT:

        - ``x`` -- a ring morphism, or a list or a tuple that define a
          ring morphism

        - ``check`` -- boolean (default: ``True``); passed onto
          functions called by this one to be more careful about input
          argument type checking

        EXAMPLES::

            sage: f = ZZ.hom(QQ); f
            Natural morphism:
              From: Integer Ring
              To:   Rational Field

            sage: H = Hom(Spec(QQ, ZZ), Spec(ZZ)); H
            Set of morphisms
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring

            sage: phi = H(f); phi
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Natural morphism:
                      From: Integer Ring
                      To:   Rational Field

        TESTS::

            sage: H._element_constructor_(f)
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Natural morphism:
                      From: Integer Ring
                      To:   Rational Field

        We illustrate input type checking::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: C = A.subscheme(x*y - 1)
            sage: H = C.Hom(C); H
            Set of morphisms
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field
                    defined by: x*y - 1
              To:   Closed subscheme of Affine Space of dimension 2 over Rational Field
                    defined by: x*y - 1
            sage: H(1)
            Traceback (most recent call last):
            ...
            TypeError: x must be a ring homomorphism, list or tuple
        """
        if isinstance(x, (list, tuple)):
            return self.domain()._morphism(self, x, check=check)

        from sage.categories.map import Map
        from sage.categories.rings import Rings
        if isinstance(x, Map) and x.category_for().is_subcategory(Rings()):
            # x is a morphism of Rings
            return SchemeMorphism_spec(self, x, check=check)

        raise TypeError("x must be a ring homomorphism, list or tuple")


# *******************************************************************
#  Base class for points
# *******************************************************************

class SchemeHomset_points(SchemeHomset_generic):
    r"""
    Set of rational points of the scheme.

    Recall that the `K`-rational points of a scheme `X` over `k` can be
    identified with the set of morphisms `\mathrm{Spec}(K) \to X`. In Sage, the
    rational points are implemented by such scheme morphisms.

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
        if check and not isinstance(X, AffineScheme):
            raise ValueError('The domain must be an affine scheme.')
        SchemeHomset_generic.__init__(self, X, Y, category=category, check=check, base=base)

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: A2 = AffineSpace(QQ, 2)
            sage: Hom = A2(QQ)
            sage: loads(Hom.dumps()) == Hom
            True
        """
        return SchemeHomset, (self.domain(), self.codomain(), self.homset_category(),
                              self.base_ring(), False, True)

    def _coerce_map_from_(self, other):
        r"""
        Return true if ``other`` canonically coerces to ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: P = ProjectiveSpace(QQ, 1, 'x')
            sage: P2 = ProjectiveSpace(R, 1, 'x')
            sage: P2(R)._coerce_map_from_(P(QQ))
            True
            sage: P(QQ)._coerce_map_from_(P2(R))
            False

        ::

            sage: P = ProjectiveSpace(QQ, 1, 'x')
            sage: P2 = ProjectiveSpace(CC, 1, 'y')                                      # needs sage.rings.real_mpfr
            sage: P2(CC)._coerce_map_from_(P(QQ))                                       # needs sage.rings.real_mpfr
            False

        ::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: H = A.subscheme(z)
            sage: L = A.subscheme([z, y + z])
            sage: A(QQ)._coerce_map_from_(H(QQ))
            True
            sage: H(QQ)._coerce_map_from_(L(QQ))                                        # needs sage.libs.singular
            True
            sage: L(QQ).has_coerce_map_from(H(QQ))                                      # needs sage.libs.singular
            False
            sage: A(CC)._coerce_map_from_(H(QQ))                                        # needs sage.rings.real_mpfr
            True
            sage: H(CC)._coerce_map_from_(L(RR))                                        # needs sage.libs.singular sage.rings.real_mpfr
            True

        ::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: A2.<u,v> = AffineSpace(QQ, 2)
            sage: A(QQ).has_coerce_map_from(A2(QQ))
            False

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: P.<u,v,w> = ProjectiveSpace(QQ, 2)
            sage: A(QQ).has_coerce_map_from(P(QQ))
            False

        ::

            sage: A = AffineSpace(QQ, 1)
            sage: A(QQ)._coerce_map_from_(ZZ)
            True

        ::

            sage: PS = ProjectiveSpace(ZZ, 1, 'x')
            sage: PS2 = ProjectiveSpace(Zp(7), 1, 'x')                                  # needs sage.rings.padics
            sage: PS(ZZ).has_coerce_map_from(PS2(Zp(7)))                                # needs sage.rings.padics
            False
            sage: PS2(Zp(7)).has_coerce_map_from(PS(ZZ))                                # needs sage.rings.padics
            True

        ::

            sage: PP1 = ProductProjectiveSpaces(ZZ, [2,1], 'x')
            sage: PP1(QQ)._coerce_map_from_(PP1(ZZ))
            True
            sage: PP2 = ProductProjectiveSpaces(QQ, [1,2], 'x')
            sage: PP2(QQ)._coerce_map_from_(PP1(ZZ))
            False
            sage: PP3 = ProductProjectiveSpaces(QQ, [2,1], 'y')
            sage: PP3(QQ)._coerce_map_from_(PP1(ZZ))
            False

        ::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(2)
            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: H = A.subscheme(z)
            sage: A(K).has_coerce_map_from(H(QQ))
            True

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: X = P.subscheme([x - y])
            sage: P(1,1) == X(1,1)
            True

        ::

            sage: A = AffineSpace(QQ, 1, 'x')
            sage: AC = AffineSpace(CC, 1, 'x')                                          # needs sage.rings.real_mpfr
            sage: A(3/2) == AC(3/2)                                                     # needs sage.rings.real_mpfr
            True

        ::

            sage: A = AffineSpace(QQ, 1)
            sage: A(0) == 0
            True
        """
        target = self.codomain()
        # ring elements can be coerced to a space if we're affine dimension 1
        # and the base rings are coercible
        if isinstance(other, CommutativeRing):
            try:
                if (isinstance(target.ambient_space(), AffineSpace)
                        and target.ambient_space().dimension_relative() == 1):
                    return target.base_ring().has_coerce_map_from(other)
                else:
                    return False
            except AttributeError:  # no .ambient_space
                return False
        elif isinstance(other, SchemeHomset_points):
        #we are converting between scheme points
            source = other.codomain()
            if isinstance(target, AlgebraicScheme_subscheme):
                #subscheme coerce when there is containment
                if not isinstance(source, AlgebraicScheme_subscheme):
                    return False
                if target.ambient_space() == source.ambient_space():
                    if all(g in source.defining_ideal()
                           for g in target.defining_polynomials()):
                        return self.domain().coordinate_ring().has_coerce_map_from(other.domain().coordinate_ring())
            else:
                #if the target is an ambient space, we can coerce if the base rings coerce
                #and they are the same type: affine, projective, etc and have the same
                #variable names
                try:
                    ta = target.ambient_space()
                    sa = source.ambient_space()
                except AttributeError: #no .ambient_space
                    return False
                #for projective and affine varieties, we check dimension
                #and matching variable names
                if ((isinstance(ta, ProjectiveSpace) and isinstance(sa, ProjectiveSpace))
                        or (isinstance(ta, AffineSpace) and isinstance(sa, AffineSpace))):
                    if (ta.variable_names() == sa.variable_names()):
                        return self.domain().coordinate_ring().has_coerce_map_from(other.domain().coordinate_ring())
                    else:
                        return False
                #for products of projective spaces, we check dimension of
                #components and matching variable names
                elif isinstance(ta, ProductProjectiveSpaces) and isinstance(sa, ProductProjectiveSpaces):
                    if (ta.dimension_relative_components() == sa.dimension_relative_components()) \
                      and (ta.variable_names() == sa.variable_names()):
                        return self.domain().coordinate_ring().has_coerce_map_from(other.domain().coordinate_ring())
                    else:
                        return False

    def _element_constructor_(self, *v, **kwds):
        """
        The element constructor.

        INPUT:

        - ``v`` -- anything that determines a scheme morphism in the
          Hom-set

        OUTPUT: the scheme morphism determined by ``v``

        EXAMPLES::

            sage: A2 = AffineSpace(ZZ, 2)
            sage: F = GF(3)
            sage: F_points = A2(F);  type(F_points)
            <class 'sage.schemes.affine.affine_homset.SchemeHomset_points_affine_with_category'>
            sage: F_points([2,5])
            (2, 2)

            sage: # needs sage.rings.finite_rings
            sage: P2 = ProjectiveSpace(GF(3), 2)
            sage: F.<a> = GF(9, 'a')
            sage: F_points = P2(F)
            sage: type(F_points)
            <class 'sage.schemes.projective.projective_homset.SchemeHomset_points_projective_field_with_category'>
            sage: F_points([4,2*a])
            (1 : 2*a : 1)

        TESTS::

            sage: F_points._element_constructor_([4,2*a])                               # needs sage.rings.finite_rings
            (1 : 2*a : 1)
        """
        if len(v) == 1:
            v = v[0]
        return self.extended_codomain()._point(self, v, **kwds)

    def __iter__(self):
        r"""
        Return an iterator for the set of rational points on this scheme.

        By default, this calls the :meth:`points` method, which is implemented
        when the base ring is a field

        - for affine homsets at :meth:`sage.schemes.affine.affine_homset.SchemeHomset_points_affine.points`;
        - for projective homsets at :meth:`sage.schemes.projective.projective_homset.SchemeHomset_points_projective_field.points`;
        - and toric homsets at :meth:`sage.schemes.toric.homset.SchemeHomset_points_toric_field._enumerator`.

        OUTPUT: iterator over points

        TESTS::

            sage: E = EllipticCurve(GF(19), [1, 0])
            sage: list(E.point_homset())
            [(0 : 1 : 0), (0 : 0 : 1), (3 : 7 : 1), (3 : 12 : 1), (4 : 7 : 1),
             (4 : 12 : 1), (5 : 4 : 1), (5 : 15 : 1), (8 : 8 : 1), (8 : 11 : 1),
             (9 : 4 : 1), (9 : 15 : 1), (12 : 7 : 1), (12 : 12 : 1), (13 : 5 : 1),
             (13 : 14 : 1), (17 : 3 : 1), (17 : 16 : 1), (18 : 6 : 1), (18 : 13 : 1)]
            sage: _ == list(E)
            True
            sage: E.point_homset().cardinality()
            20

        ::

            sage: A.<x, y> = AffineSpace(2, GF(5))
            sage: list(A.point_homset())
            [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4),
             (1, 0), (1, 1), (1, 2), (1, 3), (1, 4),
             (2, 0), (2, 1), (2, 2), (2, 3), (2, 4),
             (3, 0), (3, 1), (3, 2), (3, 3), (3, 4),
             (4, 0), (4, 1), (4, 2), (4, 3), (4, 4)]
            sage: _ == list(A)
            True
            sage: A.point_homset().cardinality()
            25

        ::

            sage: P1 = toric_varieties.P1(base_ring=GF(3))
            sage: list(P1.point_homset())
            [[0 : 1], [1 : 0], [1 : 1], [1 : 2]]
            sage: P1.point_homset().cardinality()
            4
        """
        yield from self.points()

    def extended_codomain(self):
        r"""
        Return the codomain with extended base, if necessary.

        OUTPUT:

        The codomain scheme, with its base ring extended to the
        codomain. That is, the codomain is of the form `\mathrm{Spec}(R)` and
        the base ring of the domain is extended to `R`.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: P2 = ProjectiveSpace(QQ, 2)
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + x - (3^3-3))
            sage: K_points = P2(K);  K_points
            Set of rational points of Projective Space of dimension 2
             over Number Field in a with defining polynomial x^2 + x - 24
            sage: K_points.codomain()
            Projective Space of dimension 2 over Rational Field
            sage: K_points.extended_codomain()
            Projective Space of dimension 2
             over Number Field in a with defining polynomial x^2 + x - 24
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

    def zero(self):
        """
        Return the identity of the codomain with extended base, if necessary.
        """
        return self.extended_codomain().zero()

    def _repr_(self):
        """
        Return a string representation of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: P2 = ProjectiveSpace(ZZ, 2)
            sage: P2(QQ)._repr_()
            'Set of rational points of Projective Space of dimension 2 over Rational Field'
        """
        return 'Set of rational points of '+str(self.extended_codomain())

    def value_ring(self):
        r"""
        Return `R` for a point Hom-set `X(\mathrm{Spec}(R))`.

        OUTPUT: a commutative ring

        EXAMPLES::

            sage: P2 = ProjectiveSpace(ZZ, 2)
            sage: P2(QQ).value_ring()
            Rational Field
        """
        dom = self.domain()
        if not isinstance(dom, AffineScheme):
            raise ValueError("value rings are defined for affine domains only")
        return dom.coordinate_ring()

    def cardinality(self):
        """
        Return the number of points.

        OUTPUT: integer or infinity

        EXAMPLES::

            sage: toric_varieties.P2().point_set().cardinality()                        # needs sage.geometry.polyhedron sage.graphs
            +Infinity

            sage: P2 = toric_varieties.P2(base_ring=GF(3))                              # needs sage.geometry.polyhedron sage.graphs
            sage: P2.point_set().cardinality()                                          # needs sage.geometry.polyhedron sage.graphs
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

        OUTPUT: a tuple containing all points of the toric variety

        EXAMPLES::

            sage: P1 = toric_varieties.P1(base_ring=GF(3))                              # needs sage.geometry.polyhedron sage.graphs
            sage: P1.point_set().list()                                                 # needs sage.geometry.polyhedron sage.graphs
            ([0 : 1], [1 : 0], [1 : 1], [1 : 2])
        """
        return tuple(self)
