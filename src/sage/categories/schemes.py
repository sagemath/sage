# sage_setup: distribution = sagemath-categories
r"""
Schemes
"""

# ****************************************************************************
#       Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2012 Nicolas M. Thiery <nthiery@users.sf.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.category import Category
from sage.categories.category_types import Category_over_base
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.sets_cat import Sets
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.categories.rings import Rings
from sage.categories.fields import Fields
from sage.categories.homsets import HomsetsCategory
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_import import lazy_import

lazy_import('sage.categories.map', 'Map')
lazy_import('sage.schemes.generic.morphism', 'SchemeMorphism')
lazy_import('sage.schemes.generic.scheme', 'Scheme')


class Schemes(Category):
    """
    The category of all schemes.

    EXAMPLES::

        sage: Schemes()
        Category of schemes

    ``Schemes`` can also be used to construct the category of schemes
    over a given base::

        sage: Schemes(Spec(ZZ))
        Category of schemes over Integer Ring

        sage: Schemes(ZZ)
        Category of schemes over Integer Ring

    .. TODO::

        Make ``Schemes()`` a singleton category (and remove
        :class:`Schemes` from the workaround in
        :meth:`.category_types.Category_over_base._test_category_over_bases`).

        This is currently incompatible with the dispatching below.

    TESTS::

        sage: TestSuite(Schemes()).run()

    Check that Hom sets of schemes are in the correct category::

        sage: Schemes().Homsets().super_categories()
        [Category of homsets]
    """
    @staticmethod
    def __classcall_private__(cls, X=None):
        """
        Implement the dispatching ``Schemes(ZZ)`` -> ``Schemes_over_base``.

        EXAMPLES::

            sage: Schemes()
            Category of schemes

            sage: Schemes(Spec(ZZ))
            Category of schemes over Integer Ring

            sage: Schemes(ZZ)
            Category of schemes over Integer Ring
        """
        if X is not None:
            if X not in Schemes():
                X = Schemes()(X)
            return Schemes_over_base(X)
        return super().__classcall__(cls)

    def super_categories(self):
        """
        EXAMPLES::

            sage: Schemes().super_categories()
            [Category of sets]
        """
        return [Sets()]

    def _call_(self, x):
        """
        Construct a scheme from the data in ``x``.

        EXAMPLES:

        Let us first construct the category of schemes::

            sage: S = Schemes(); S
            Category of schemes

        We create a scheme from a ring::

            sage: X = S(ZZ); X                  # indirect doctest
            Spectrum of Integer Ring

        We create a scheme from a scheme (do nothing)::

            sage: S(X)
            Spectrum of Integer Ring

        We create a scheme morphism from a ring homomorphism.x::

            sage: phi = ZZ.hom(QQ); phi
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
            sage: f = S(phi); f                 # indirect doctest
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Natural morphism:
                      From: Integer Ring
                      To:   Rational Field

            sage: f.domain()
            Spectrum of Rational Field
            sage: f.codomain()
            Spectrum of Integer Ring
            sage: S(f)                          # indirect doctest
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Natural morphism:
                      From: Integer Ring
                      To:   Rational Field
        """
        if isinstance(x, (SchemeMorphism, Scheme)):
            return x
        from sage.categories.commutative_rings import CommutativeRings
        from sage.schemes.generic.spec import Spec
        if x in CommutativeRings():
            return Spec(x)
        elif isinstance(x, Map) and x.category_for().is_subcategory(Rings()):
            # x is a morphism of Rings
            A = Spec(x.codomain())
            return A.hom(x)
        else:
            raise TypeError("No way to create an object or morphism in %s from %s" % (self, x))


class Schemes_over_base(Category_over_base):
    """
    The category of schemes over a given base scheme.

    EXAMPLES::

        sage: Schemes(Spec(ZZ))
        Category of schemes over Integer Ring

    TESTS::

        sage: C = Schemes(ZZ)
        sage: TestSuite(C).run()
    """
    def base_scheme(self):
        """
        EXAMPLES::

            sage: Schemes(Spec(ZZ)).base_scheme()
            Spectrum of Integer Ring
        """
        return self.base()

    def super_categories(self):
        """
        EXAMPLES::

            sage: Schemes(Spec(ZZ)).super_categories()
            [Category of schemes]
        """
        return [Schemes()]

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Schemes(Spec(ZZ)) # indirect doctest
            Category of schemes over Integer Ring
        """
        from sage.schemes.generic.scheme import AffineScheme
        base = self.base()
        if isinstance(base, AffineScheme):
            base = base.coordinate_ring()
        return f"schemes over {base}"


class AbelianVarieties(Schemes_over_base):
    r"""
    The category of abelian varieties over a given field scheme.

    EXAMPLES::

        sage: AbelianVarieties(QQ)
        Category of abelian varieties over Rational Field
        sage: AbelianVarieties(ZZ)
        Traceback (most recent call last):
        ...
        ValueError: category of abelian varieties is only defined over fields
    """
    @staticmethod
    def __classcall__(cls, base):
        r"""
        Normalise arguments for the constructor of the
        :class:`AbelianVarieties` category.

        EXAMPLES::

            sage: AbelianVarieties(QQ) is AbelianVarieties(Spec(QQ))  # indirect doctest
            True
        """
        from sage.schemes.generic.scheme import Scheme, AffineScheme
        if not isinstance(base, Schme):
            base = Schemes()(base)
        if not (isinstance(base, AffineScheme) and base.coordinate_ring() in Fields()):
            raise ValueError('category of abelian varieties is only defined over fields')
        return UniqueRepresentation.__classcall__(cls, base)

    def __init__(self, base):
        r"""
        Cosntructor for the :class:`AbelianVarieties` category.

        EXAMPLES::

            sage: A = AbelianVarieties(QQ); A
            Category of abelian varieties over Rational Field
            sage: B = AbelianVarieties(Spec(QQ)); B
            Category of abelian varieties over Rational Field
            sage: A is B
            True

        The category of abelian varieties is only defined over fields::

            sage: AbelianVarieties(ZZ)
            Traceback (most recent call last):
            ...
            ValueError: category of abelian varieties is only defined over fields
        """
        super().__init__(base)

    def base_scheme(self):
        """
        EXAMPLES::

            sage: Schemes(Spec(ZZ)).base_scheme()
            Spectrum of Integer Ring
        """
        base = self.base()
        if base not in Schemes():
            base = Schemes()(base)
        return base

    def super_categories(self):
        """
        EXAMPLES::

            sage: AbelianVarieties(QQ).super_categories()
            [Category of schemes over Rational Field]
        """
        return [Schemes(self.base_scheme())]

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: AbelianVarieties(Spec(QQ))  # indirect doctest
            Category of abelian varieties over Rational Field
        """
        from sage.schemes.generic.scheme import is_AffineScheme
        if is_AffineScheme(self.base_scheme()):
            return "abelian varieties over %s" % self.base_scheme().coordinate_ring()
        else:
            return "abelian varieties over %s" % self.base_scheme()

    class ParentMethods:
        def zero(self):
            """
            Return the additive identity of this abelian variety.

            EXAMPLES::

                sage: E = EllipticCurve([1, 3])
                sage: E.zero()
                (0 : 1 : 0)
                sage: E.zero() == E(0)
                True
            """
            try:
                return self(0)
            except Exception:
                try:
                    return self.point_homset()(0)
                except Exception:
                    raise NotImplementedError

    class Homsets(HomsetsCategory):
        r"""
        Overloaded ``Homsets`` class to register the homset
        as an additive abelian group.

        EXAMPLES::

            sage: AbelianVarieties(QQ).Homsets().is_subcategory(CommutativeAdditiveGroups())
            True
        """
        def extra_super_categories(self):
            r"""
            Register the homset as an additive abelian group.

            EXAMPLES::

                sage: Hom(EllipticCurve(j=1), EllipticCurve(j=2)) in CommutativeAdditiveGroups()
                True
            """
            return [CommutativeAdditiveGroups()]

        class Endset(CategoryWithAxiom):
            r"""
            Overloaded ``Endset`` class to register the endset
            as a ring.

            sage: AbelianVarieties(QQ).Endsets().is_subcategory(Rings())
            True
            """
            def extra_super_categories(self):
                r"""
                Register the endset as a ring.

                EXAMPLES::

                    sage: End(EllipticCurve(j=1)) in Rings()
                    True
                """
                return [Rings()]


class Jacobians(Schemes_over_base):
    """
    The category of Jacobians attached to curves or function fields.

    EXAMPLES::

        sage: Jacobians(QQ)
        Category of Jacobians over Rational Field

    TESTS::

        sage: TestSuite(Jacobians(QQ)).run()
    """
    def __init__(self, base):
        r"""
        Constructor of this category.

        EXAMPLES::

            sage: Jacobians(QQ)
            Category of Jacobians over Rational Field
            sage: Jacobians(Spec(QQ))
            Category of Jacobians over Rational Field
        """
        from sage.schemes.generic.scheme import AffineScheme
        if isinstance(base, AffineScheme):
            base = base.coordinate_ring()
        if base not in Fields():
            raise ValueError('category of Jacobians is only defined over fields')
        super().__init__(base)

    def base_scheme(self):
        """
        Return the base scheme of this Jacobians category.

        EXAMPLES::

            sage: Jacobians(QQ).base_scheme()
            Spectrum of Rational Field
        """
        base = self.base()
        if base not in Schemes():
            base = Schemes()(base)
        return base

    def super_categories(self):
        """
        Return the super categories of this Jacobians category.

        EXAMPLES::

            sage: Jacobians(QQ).super_categories()
            [Category of abelian varieties over Rational Field]
        """
        return [AbelianVarieties(self.base_scheme())]

    def _repr_object_names(self):
        """
        Return the string representation of this category.

        EXAMPLES::

            sage: Jacobians(Spec(QQ))  # indirect doctest
            Category of Jacobians over Rational Field
        """
        return "Jacobians over %s" % self.base()

    class ParentMethods:

        @abstract_method
        def base_curve(self):
            """
            Return the curve to which this Jacobian is attached.

            EXAMPLES::

                sage: # needs sage.rings.function_field
                sage: K.<x> = FunctionField(GF(2))
                sage: J = K.jacobian()
                sage: J.base_curve()
                Rational function field in x over Finite Field of size 2
            """
