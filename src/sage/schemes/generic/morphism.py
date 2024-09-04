r"""
Scheme morphism

.. NOTE::

    You should never create the morphisms directly. Instead, use the
    :meth:`~sage.schemes.generic.scheme.hom` and
    :meth:`~sage.structure.parent.Hom` methods that are inherited by
    all schemes.

If you want to extend the Sage library with some new kind of scheme,
your new class (say, ``MyScheme``) should provide a method

* ``MyScheme._morphism(*args, **kwds)`` returning a morphism
  between two schemes in your category, usually defined via
  polynomials. Your morphism class should derive from
  :class:`SchemeMorphism_polynomial`. These morphisms will usually be
  elements of the Hom-set
  :class:`~sage.schemes.generic.homset.SchemeHomset_generic`.

Optionally, you can also provide a special Hom-set class for your
subcategory of schemes. If you want to do this, you should also
provide a method

* ``MyScheme._homset(*args, **kwds)`` returning a
  Hom-set, which must be an element of a derived class of
  :class:`~sage.schemes.generic.homset.SchemeHomset_generic`. If your
  new Hom-set class does not use ``MyScheme._morphism`` then you
  do not have to provide it.

Note that points on schemes are morphisms `\mathrm{Spec}(K)\to X`, too. But we
typically use a different notation, so they are implemented in a different
derived class. For this, you should implement a method

* ``MyScheme._point(*args, **kwds)`` returning a point, that is, a morphism
  `\mathrm{Spec}(K)\to X`. Your point class should derive from
  :class:`SchemeMorphism_point`.

Optionally, you can also provide a special Hom-set for the points, for
example the point Hom-set can provide a method to enumerate all
points. If you want to do this, you should also provide a method

* ``MyScheme._point_homset(*args, **kwds)`` returning
  the :mod:`~sage.schemes.generic.homset` of points. The Hom-sets of
  points are implemented in classes named ``SchemeHomset_points_...``.
  If your new Hom-set class does not use ``MyScheme._point`` then
  you do not have to provide it.

AUTHORS:

- David Kohel, William Stein

- William Stein (2006-02-11): fixed bug where P(0,0,0) was allowed as
  a projective point.

- Volker Braun (2011-08-08): Renamed classes, more documentation, misc
  cleanups.

- Ben Hutz (June 2012): added support for projective ring

- Simon King (2013-10): copy the changes of :class:`~sage.categories.morphism.Morphism`
  that have been introduced in :issue:`14711`.
"""

# ****************************************************************************
#       Copyright (C) 2013 Simon King <simon.king@uni-jena.de>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import operator

from sage.arith.power import generic_power
from sage.categories.homset import Homset, Hom, End
from sage.categories.map import FormalCompositeMap, Map
from sage.categories.morphism import SetMorphism
from sage.misc.constant_function import ConstantFunction
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.fraction_field import FractionField_generic
from sage.rings.fraction_field_element import FractionFieldElement
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme
from sage.structure.element import Element, parent, coercion_model
from sage.structure.richcmp import richcmp
from sage.structure.sequence import Sequence


def is_SchemeMorphism(f):
    """
    Test whether ``f`` is a scheme morphism.

    INPUT:

    - ``f`` -- anything

    OUTPUT:

    boolean; return ``True`` if ``f`` is a scheme morphism or a point
    on an elliptic curve.

    EXAMPLES::

        sage: A.<x,y> = AffineSpace(QQ, 2); H = A.Hom(A)
        sage: f = H([y, x^2 + y]); f
        Scheme endomorphism of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to (y, x^2 + y)
        sage: from sage.schemes.generic.morphism import is_SchemeMorphism
        sage: is_SchemeMorphism(f)
        doctest:warning...
        DeprecationWarning: The function is_SchemeMorphism is deprecated;
        use 'isinstance(..., SchemeMorphism)' instead.
        See https://github.com/sagemath/sage/issues/38296 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38296,
                "The function is_SchemeMorphism is deprecated; "
                "use 'isinstance(..., SchemeMorphism)' instead.")
    return isinstance(f, SchemeMorphism)


class SchemeMorphism(Element):
    """
    Base class for scheme morphisms.

    INPUT:

    - ``parent`` -- the parent of the morphism

    .. TODO::

        For historical reasons, :class:`SchemeMorphism` copies code from
        :class:`~sage.categories.map.Map` rather than inheriting from it.
        Proper inheritance should be used instead. See :issue:`14711`.

    EXAMPLES::

        sage: X = Spec(ZZ)
        sage: Hom = X.Hom(X)
        sage: from sage.schemes.generic.morphism import SchemeMorphism
        sage: f = SchemeMorphism(Hom)
        sage: type(f)
        <class 'sage.schemes.generic.morphism.SchemeMorphism'>

    TESTS::

        sage: A2 = AffineSpace(QQ, 2)
        sage: A2.structure_morphism().domain()
        Affine Space of dimension 2 over Rational Field
        sage: A2.structure_morphism().category()
        Category of homsets of schemes
    """

    def __init__(self, parent, codomain=None):
        """
        The Python constructor.

        EXAMPLES::

            sage: X = Spec(ZZ)
            sage: Hom = X.Hom(X)
            sage: from sage.schemes.generic.morphism import SchemeMorphism
            sage: f = SchemeMorphism(Hom)
            sage: type(f)
            <class 'sage.schemes.generic.morphism.SchemeMorphism'>
        """
        if codomain is not None:
            parent = Hom(parent, codomain)
        if not isinstance(parent, Homset):
            raise TypeError("parent (=%s) must be a Homspace" % parent)
        Element.__init__(self, parent)
        self._codomain = parent.codomain()

    @lazy_attribute
    def domain(self):
        r"""
        The constant function from the domain.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ['x,y'])
            sage: H = A.Hom(A)
            sage: f = H([y, x^2 + y])
            sage: f.domain() is A
            True
        """
        return ConstantFunction(self.parent().domain())

    @lazy_attribute
    def codomain(self):
        r"""
        The constant function from the codomain.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ['x,y'])
            sage: H = A.Hom(A)
            sage: f = H([y, x^2 + y])
            sage: f.codomain() is A
            True
        """
        return ConstantFunction(self._codomain)

    # We copy methods of sage.categories.map.Map, to make
    # a future transition of SchemeMorphism to a sub-class of Morphism
    # easier.
    def __call__(self, x, *args, **kwds):
        """
        Do not override this method!

        For implementing application of maps, implement a method
        ``_call_(self, x)`` and/or a method ``_call_with_args(x, args, kwds)`.
        In these methods, you can assume that ``x`` belongs to the domain of
        this morphism, ``args`` is a tuple and ``kwds`` is a dict.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: H = A.Hom(A)
            sage: f = H([y, x^2 + y])
            sage: f([2,3])    # indirect doctest
            (3, 7)

        An example with optional arguments::

            sage: PS.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(PS, PS)
            sage: f = H([x^3, x*y^2])
            sage: P = PS(0, 1)
            sage: f(P, check=False)     # indirect doctest
            (0 : 0)
        """
        P = parent(x)
        D = self.domain()
        if P is D: # we certainly want to call _call_/with_args
            if not args and not kwds:
                return self._call_(x)
            return self._call_with_args(x, args, kwds)
        # Is there coercion?
        converter = D._internal_coerce_map_from(P)
        if converter is None:
            try:
                return self.pushforward(x,*args,**kwds)
            except (AttributeError, TypeError, NotImplementedError):
                pass # raise TypeError, "%s must be coercible into %s"%(x, self.domain())
            # Here, we would like to do
            ##try:
            ##    x = D(x).
            ##except (TypeError, NotImplementedError):
            ##    raise TypeError, "%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(x, self.domain())
            # However, this would involve a test whether x.codomain() ==
            # self. This would trigger a Groebner basis computation, that
            # (1) could be slow and (2) could involve an even slower toy
            # implementation, resulting in a warning.
            #
            # Contract: If x is a scheme morphism point, then _call_ knows
            # what to do with it (e.g., use the _coords attribute). Otherwise,
            # we can try a conversion into the domain (e.g., if x is a list),
            # WITHOUT to trigger a Groebner basis computation.
            if kwds.get('check', True):
                if not isinstance(x, SchemeMorphism_point):
                    try:
                        x = D(x)
                    except (TypeError, NotImplementedError):
                        raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented" % (x, self.domain()))
                elif self.domain() != x.codomain():
                    raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented" % (x, self.domain()))
        else:
            x = converter(x)
        if not args and not kwds:
            return self._call_(x)
        return self._call_with_args(x, args, kwds)

    def _repr_defn(self):
        r"""
        Return a string representation of the definition of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: X = Spec(ZZ)
            sage: Hom = X.Hom(X)
            sage: from sage.schemes.generic.morphism import SchemeMorphism
            sage: f = SchemeMorphism(Hom)
            sage: f._repr_defn()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _repr_type(self):
        r"""
        Return a string representation of the type of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: A2 = AffineSpace(QQ, 2)
            sage: A2.structure_morphism()  # indirect doctest
            Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Spectrum of Rational Field
              Defn: Structure map
        """
        return "Scheme"

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: X = Spec(ZZ)
            sage: Hom = X.Hom(X)
            sage: from sage.schemes.generic.morphism import SchemeMorphism
            sage: f = SchemeMorphism(Hom)
            sage: f._repr_()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if self.is_endomorphism():
            s = "%s endomorphism of %s" % (self._repr_type(), self.domain())
        else:
            s = "%s morphism:" % self._repr_type()
            s += "\n  From: %s" % self.domain()
            s += "\n  To:   %s" % self._codomain
        d = self._repr_defn()
        if d != '':
            s += "\n  Defn: %s" % ('\n        '.join(self._repr_defn().split('\n')))
        return s

    def __mul__(self, right):
        """
        We can currently only multiply scheme morphisms.

        If one factor is an identity morphism, the other is returned.
        Otherwise, a formal composition of maps obtained from the scheme
        morphisms is returned.

        EXAMPLES:

        Identity maps do not contribute to the product::

            sage: X = AffineSpace(QQ, 2)
            sage: id = X.identity_morphism()
            sage: id^0    # indirect doctest
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Identity map
            sage: id^2
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Identity map

        Here, we see a formal composition::

            sage: X = AffineSpace(QQ, 2)
            sage: f = X.structure_morphism()
            sage: Y = Spec(QQ)
            sage: g = Y.structure_morphism()
            sage: g * f    # indirect doctest
            Composite map:
              From: Affine Space of dimension 2 over Rational Field
              To:   Spectrum of Integer Ring
              Defn:   Generic morphism:
                      From: Affine Space of dimension 2 over Rational Field
                      To:   Spectrum of Rational Field
                    then
                      Generic morphism:
                      From: Spectrum of Rational Field
                      To:   Spectrum of Integer Ring

        Of course, the codomain of the first factor must coincide with the
        domain of the second factor::

            sage: f * g
            Traceback (most recent call last):
            ...
            TypeError: self (=Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Spectrum of Rational Field
              Defn: Structure map) domain must equal right (=Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Structure map) codomain
        """
        if not isinstance(right, SchemeMorphism):
            return coercion_model.bin_op(self, right, operator.mul)
        if right.codomain() != self.domain():
            raise TypeError("self (=%s) domain must equal right (=%s) codomain" % (self, right))
        if isinstance(self, SchemeMorphism_id):
            return right
        if isinstance(right, SchemeMorphism_id):
            return self
        return self._composition(right)

    def __pow__(self, n, dummy=None):
        """
        Exponentiate an endomorphism.

        INPUT:

        - ``n`` -- integer; the exponent

        OUTPUT: a composite map that belongs to the same endomorphism set as ``self``

        EXAMPLES::

            sage: X = AffineSpace(QQ, 2)
            sage: id = X.identity_morphism()
            sage: id^0
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Identity map
            sage: id^2
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Identity map
        """
        if not self.is_endomorphism():
            raise TypeError("self must be an endomorphism.")
        if n == 0:
            return self.domain().identity_morphism()
        return generic_power(self, n)

    def category(self):
        """
        Return the category of the Hom-set.

        OUTPUT: a category

        EXAMPLES::

            sage: A2 = AffineSpace(QQ, 2)
            sage: A2.structure_morphism().category()
            Category of homsets of schemes
        """
        return self.parent().category()

    def category_for(self):
        """
        Return the category which this morphism belongs to.

        EXAMPLES::

            sage: A2 = AffineSpace(QQ, 2)
            sage: A2.structure_morphism().category_for()
            Category of schemes
        """
        return self.parent().homset_category()

    def is_endomorphism(self) -> bool:
        """
        Return whether the morphism is an endomorphism.

        OUTPUT: boolean; whether the domain and codomain are identical

        EXAMPLES::

            sage: X = AffineSpace(QQ, 2)
            sage: X.structure_morphism().is_endomorphism()
            False
            sage: X.identity_morphism().is_endomorphism()
            True
        """
        return self.parent().is_endomorphism_set()

    def base_ring(self):
        r"""
        Return the base ring of ``self``, that is, the ring over which
        the defining polynomials of ``self`` are defined.

        OUTPUT: ring

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([3/5*x^2, 6*y^2])
            sage: f.base_ring()
            Rational Field

        ::

            sage: R.<t> = PolynomialRing(ZZ, 1)
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: H = Hom(P, P)
            sage: f = H([3*x^2, y^2])
            sage: f.base_ring()
            Multivariate Polynomial Ring in t over Integer Ring

        Points have correct base rings too (:issue:`34336`)::

            sage: x = P(t, 5); x
            (t : 5)
            sage: x.base_ring()
            Multivariate Polynomial Ring in t over Integer Ring

        ::

            sage: # needs sage.rings.finite_rings sage.schemes
            sage: E = EllipticCurve(GF((17,2)), [1,2,3,4,5])
            sage: P = E.random_point()
            sage: P.base_ring()
            Finite Field in z2 of size 17^2
        """
        return self.domain().base_ring()

    def _composition(self, right):
        """
        A helper for multiplying maps by composition.

        .. WARNING::

            Do not override this method! Override :meth:`_composition_`
            instead.

        EXAMPLES::

            sage: X = AffineSpace(QQ, 2)
            sage: f = X.structure_morphism()
            sage: Y = Spec(QQ)
            sage: g = Y.structure_morphism()
            sage: g * f    # indirect doctest
             Composite map:
              From: Affine Space of dimension 2 over Rational Field
              To:   Spectrum of Integer Ring
              Defn:   Generic morphism:
                      From: Affine Space of dimension 2 over Rational Field
                      To:   Spectrum of Rational Field
                    then
                      Generic morphism:
                      From: Spectrum of Rational Field
                      To:   Spectrum of Integer Ring

            sage: f * g
            Traceback (most recent call last):
            ...
            TypeError: self (=Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Spectrum of Rational Field
              Defn: Structure map) domain
            must equal right (=Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Structure map) codomain
        """
        category = self.category_for()._meet_(right.category_for())
        H = Hom(right.domain(), self._codomain, category)
        return self._composition_(right, H)

    def _composition_(self, right, homset):
        """
        Helper to construct the composition of two morphisms.

        Override this if you want to have a different behaviour of composition

        INPUT:

        - ``right`` -- a map or callable
        - ``homset`` -- a homset containing the composed map

        OUTPUT:

        An element of ``homset``. The output is obtained by converting the
        arguments to :class:`~sage.categories.morphism.SetMorphism` if
        necessary, and then forming a :class:`~sage.categories.map.FormalCompositeMap`

        EXAMPLES::

            sage: X = AffineSpace(QQ, 2)
            sage: f = X.structure_morphism()
            sage: Y = Spec(QQ)
            sage: g = Y.structure_morphism()
            sage: g * f    # indirect doctest
            Composite map:
              From: Affine Space of dimension 2 over Rational Field
              To:   Spectrum of Integer Ring
              Defn:   Generic morphism:
                      From: Affine Space of dimension 2 over Rational Field
                      To:   Spectrum of Rational Field
                    then
                      Generic morphism:
                      From: Spectrum of Rational Field
                      To:   Spectrum of Integer Ring
        """
        if not isinstance(right, Map):
            right = SetMorphism(right.parent(), right)
        return FormalCompositeMap(homset, right, SetMorphism(self.parent(),self))

    def glue_along_domains(self, other):
        r"""
        Glue two morphisms.

        INPUT:

        - ``other`` -- a scheme morphism with the same domain

        OUTPUT:

        Assuming that ``self`` and ``other`` are open immersions with the same
        domain, return scheme obtained by gluing along the images.

        EXAMPLES:

        We construct a scheme isomorphic to the projective line over
        `\mathrm{Spec}(\QQ)` by gluing two copies of `\mathbb{A}^1`
        minus a point::

            sage: # needs sage.libs.singular
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S.<xbar, ybar> = R.quotient(x*y - 1)
            sage: Rx = PolynomialRing(QQ, 'x')
            sage: i1 = Rx.hom([xbar])
            sage: Ry = PolynomialRing(QQ, 'y')
            sage: i2 = Ry.hom([ybar])
            sage: Sch = Schemes()
            sage: f1 = Sch(i1)
            sage: f2 = Sch(i2)

        Now f1 and f2 have the same domain, which is a
        `\mathbb{A}^1` minus a point. We glue along the domain::

            sage: # needs sage.libs.singular
            sage: P1 = f1.glue_along_domains(f2); P1
            Scheme obtained by gluing X and Y along U, where
              X: Spectrum of Univariate Polynomial Ring in x over Rational Field
              Y: Spectrum of Univariate Polynomial Ring in y over Rational Field
              U: Spectrum of Quotient of Multivariate Polynomial Ring in x, y
                 over Rational Field by the ideal (x*y - 1)
            sage: a, b = P1.gluing_maps()
            sage: a
            Affine Scheme morphism:
              From: Spectrum of Quotient of Multivariate Polynomial Ring in x, y
                    over Rational Field by the ideal (x*y - 1)
              To:   Spectrum of Univariate Polynomial Ring in x over Rational Field
              Defn: Ring morphism:
                      From: Univariate Polynomial Ring in x over Rational Field
                      To:   Quotient of Multivariate Polynomial Ring in x, y over
                            Rational Field by the ideal (x*y - 1)
                      Defn: x |--> xbar
            sage: b
            Affine Scheme morphism:
              From: Spectrum of Quotient of Multivariate Polynomial Ring in x, y
                    over Rational Field by the ideal (x*y - 1)
              To:   Spectrum of Univariate Polynomial Ring in y over Rational Field
              Defn: Ring morphism:
                      From: Univariate Polynomial Ring in y over Rational Field
                      To:   Quotient of Multivariate Polynomial Ring in x, y over
                            Rational Field by the ideal (x*y - 1)
                      Defn: y |--> ybar
        """
        from . import glue
        return glue.GluedScheme(self, other)


class SchemeMorphism_id(SchemeMorphism):
    """
    Return the identity morphism from `X` to itself.

    INPUT:

    - ``X`` -- the scheme

    EXAMPLES::

        sage: X = Spec(ZZ)
        sage: X.identity_morphism()  # indirect doctest
        Scheme endomorphism of Spectrum of Integer Ring
          Defn: Identity map
    """
    def __init__(self, X):
        """
        The Python constructor.

        See :class:`SchemeMorphism_id` for details.

        TESTS::

            sage: Spec(ZZ).identity_morphism()
            Scheme endomorphism of Spectrum of Integer Ring
              Defn: Identity map
        """
        SchemeMorphism.__init__(self, X.Hom(X))

    def _repr_defn(self):
        r"""
        Return a string representation of the definition of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: Spec(ZZ).identity_morphism()._repr_defn()
            'Identity map'
        """
        return 'Identity map'


class SchemeMorphism_structure_map(SchemeMorphism):
    r"""
    The structure morphism.

    INPUT:

    - ``parent`` -- Hom-set with codomain equal to the base scheme of
      the domain

    EXAMPLES::

        sage: Spec(ZZ).structure_morphism()    # indirect doctest
        Scheme endomorphism of Spectrum of Integer Ring
          Defn: Structure map
    """
    def __init__(self, parent, codomain=None):
        """
        The Python constructor.

        See :class:`SchemeMorphism_structure_map` for details.

        TESTS::

            sage: from sage.schemes.generic.morphism import SchemeMorphism_structure_map
            sage: SchemeMorphism_structure_map( Spec(QQ).Hom(Spec(ZZ)) )
            Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Structure map
        """
        SchemeMorphism.__init__(self, parent, codomain=None)
        if self.domain().base_scheme() != self._codomain:
            raise ValueError("parent must have codomain equal the base scheme of domain.")

    def _repr_defn(self):
        r"""
        Return a string representation of the definition of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: Spec(ZZ).structure_morphism()._repr_defn()
            'Structure map'
        """
        return 'Structure map'


class SchemeMorphism_spec(SchemeMorphism):
    """
    Morphism of spectra of rings.

    INPUT:

    - ``parent`` -- Hom-set whose domain and codomain are affine schemes

    - ``phi`` -- a ring morphism with matching domain and codomain

    - ``check`` -- boolean (default: ``True``); whether to
      check the input for consistency

    EXAMPLES::

        sage: R.<x> = PolynomialRing(QQ)
        sage: phi = R.hom([QQ(7)]); phi
        Ring morphism:
          From: Univariate Polynomial Ring in x over Rational Field
          To:   Rational Field
          Defn: x |--> 7

        sage: X = Spec(QQ); Y = Spec(R)
        sage: f = X.hom(phi); f
        Affine Scheme morphism:
          From: Spectrum of Rational Field
          To:   Spectrum of Univariate Polynomial Ring in x over Rational Field
          Defn: Ring morphism:
                  From: Univariate Polynomial Ring in x over Rational Field
                  To:   Rational Field
                  Defn: x |--> 7

        sage: f.ring_homomorphism()
        Ring morphism:
          From: Univariate Polynomial Ring in x over Rational Field
          To:   Rational Field
          Defn: x |--> 7
    """
    def __init__(self, parent, phi, check=True):
        """
        The Python constructor.

        See :class:`SchemeMorphism_structure_map` for details.

        TESTS::

            sage: from sage.schemes.generic.morphism import SchemeMorphism_spec
            sage: SchemeMorphism_spec(Spec(QQ).Hom(Spec(ZZ)), ZZ.hom(QQ))
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Natural morphism:
                      From: Integer Ring
                      To:   Rational Field
        """
        SchemeMorphism.__init__(self, parent)
        if check:
            from sage.categories.rings import Rings
            if not (isinstance(phi, Map) and phi.category_for().is_subcategory(Rings())):
                raise TypeError("phi (=%s) must be a ring homomorphism" % phi)
            if phi.domain() != parent.codomain().coordinate_ring():
                raise TypeError("phi (=%s) must have domain %s"
                                % (phi, parent.codomain().coordinate_ring()))
            if phi.codomain() != parent.domain().coordinate_ring():
                raise TypeError("phi (=%s) must have codomain %s"
                                % (phi, parent.domain().coordinate_ring()))
        self.__ring_homomorphism = phi

    def _call_(self, x):
        r"""
        Make morphisms callable.

        INPUT:

        - ``x`` -- a scheme point

        OUTPUT: the image scheme point

        EXAMPLES:

        The following fails because inverse images of prime ideals
        under ring homomorphisms are not yet implemented::

            sage: R.<x> = PolynomialRing(QQ)
            sage: phi = R.hom([QQ(7)])
            sage: X = Spec(QQ); Y = Spec(R)
            sage: f = X.hom(phi)
            sage: f(X.an_element())    # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: inverse not implemented for morphisms of
            Rational Field
        """
        # By virtue of argument preprocessing in __call__, we can assume that
        # x is a topological scheme point of self
        S = self.ring_homomorphism().inverse_image(x.prime_ideal())
        return self._codomain(S)

    def _repr_type(self):
        r"""
        Return a string representation of the type of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: phi = R.hom([QQ(7)])
            sage: X = Spec(QQ); Y = Spec(R)
            sage: f = X.hom(phi)
            sage: f._repr_type()
            'Affine Scheme'
        """
        return "Affine Scheme"

    def _repr_defn(self):
        r"""
        Return a string representation of the definition of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: phi = R.hom([QQ(7)])
            sage: X = Spec(QQ); Y = Spec(R)
            sage: f = X.hom(phi)
            sage: print(f._repr_defn())
            Ring morphism:
              From: Univariate Polynomial Ring in x over Rational Field
              To:   Rational Field
              Defn: x |--> 7
        """
        return repr(self.ring_homomorphism())

    def ring_homomorphism(self):
        """
        Return the underlying ring homomorphism.

        OUTPUT: a ring homomorphism

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: phi = R.hom([QQ(7)])
            sage: X = Spec(QQ); Y = Spec(R)
            sage: f = X.hom(phi)
            sage: f.ring_homomorphism()
            Ring morphism:
              From: Univariate Polynomial Ring in x over Rational Field
              To:   Rational Field
              Defn: x |--> 7
        """
        return self.__ring_homomorphism


############################################################################
# Morphisms between schemes given on points. The _affine and _projective below
# refer to the CODOMAIN. The domain can be either affine or projective
# regardless of the class
############################################################################

class SchemeMorphism_polynomial(SchemeMorphism):
    r"""
    A morphism of schemes determined by polynomials that define what
    the morphism does on points in the ambient space.

    INPUT:

    - ``parent`` -- Hom-set whose domain and codomain are affine or
      projective schemes

    - ``polys`` -- list/tuple/iterable of polynomials defining the
      scheme morphism

    - ``check`` -- boolean (default: ``True``); whether to
      check the input for consistency

    EXAMPLES:

    An example involving the affine plane::

        sage: R.<x,y> = QQ[]
        sage: A2 = AffineSpace(R)
        sage: H = A2.Hom(A2)
        sage: f = H([x - y, x*y])
        sage: f([0, 1])
        (-1, 0)

    An example involving the projective line::

        sage: R.<x,y> = QQ[]
        sage: P1 = ProjectiveSpace(R)
        sage: H = P1.Hom(P1)
        sage: f = H([x^2 + y^2, x*y])
        sage: f([0, 1])
        (1 : 0)

    Some checks are performed to make sure the given polynomials
    define a morphism::

        sage: f = H([exp(x),exp(y)])                                                    # needs sage.symbolic
        Traceback (most recent call last):
        ...
        TypeError: polys (=[e^x, e^y]) must be elements of Multivariate
        Polynomial Ring in x, y over Rational Field
    """
    def __init__(self, parent, polys, check=True):
        """
        The Python constructor.

        See :class:`SchemeMorphism_polynomial` for details.

        EXAMPLES::

            sage: A2.<x,y> = AffineSpace(QQ, 2)
            sage: H = A2.Hom(A2)
            sage: H([x - y, x*y])
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to (x - y, x*y)
        """
        if check:
            if not isinstance(polys, (list, tuple)):
                raise TypeError("polys (=%s) must be a list or tuple" % polys)
            source_ring = parent.domain().ambient_space().coordinate_ring()
            target = parent._codomain.ambient_space()
            if len(polys) != target.ngens():
                raise ValueError("there must be %s polynomials" % target.ngens())
            F = []
            for poly in polys:
                try:
                    p = source_ring(poly)
                except TypeError:
                    try:
                        p = source_ring(poly.lift())
                    except (TypeError, AttributeError):
                        try:
                            p = source_ring(poly.numerator()) / source_ring(poly.denominator())
                        except (TypeError, AttributeError):
                            raise TypeError("polys (=%s) must be elements of %s" % (polys, source_ring))
                F.append(p)
            polys = Sequence(F)

        self._polys = tuple(polys)

        SchemeMorphism.__init__(self, parent)

    def __eq__(self, other):
        """
        Check equality of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a morphism

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(2, QQ)
            sage: I = A.identity_morphism()
            sage: I.parent().identity() == I
            True
        """
        if isinstance(other, SchemeMorphism_polynomial):
            if self.parent() == other.parent() and self._polys == other._polys:
                return True
        raise TypeError('cannot determine equality')

    def defining_polynomials(self):
        """
        Return the defining polynomials.

        OUTPUT:

        An immutable sequence of polynomials that defines this scheme
        morphism.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: H = A.Hom(A)
            sage: H([x^3 + y, 1 - x - y]).defining_polynomials()
            (x^3 + y, -x - y + 1)
        """
        return self._polys

    def _call_(self, x):
        """
        Apply this morphism to a point in the domain.

        INPUT:

        - ``x`` -- a point in the domain or a list or tuple that
          defines a point in the domain

        OUTPUT: a point in the codomain

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: H = A.Hom(A)
            sage: f = H([y, x^2 + y])
            sage: f([2, 3])    # indirect doctest
            (3, 7)

        An example with algebraic schemes::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: X = A.subscheme(x)
            sage: Y = A.subscheme(y)
            sage: Hom_XY = X.Hom(Y)
            sage: f = Hom_XY([y,0])   # (0,y) |-> (y,0)
            sage: f
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field
                    defined by: x
              To:   Closed subscheme of Affine Space of dimension 2 over Rational Field
                    defined by: y
              Defn: Defined on coordinates by sending (x, y) to (y, 0)
            sage: f([0,3])
            (3, 0)

        The input must be convertible into the map's domain::

            sage: f(0)
            Traceback (most recent call last):
            ...
            TypeError: 0 fails to convert into the map's domain Closed
            subscheme of Affine Space of dimension 2 over Rational Field
            defined by:
              x, but a `pushforward` method is not properly implemented

        It is possible to avoid the checks on the resulting point which can be
        useful for indeterminacies, but be careful!! ::

            sage: PS.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(PS, PS)
            sage: f = H([x^3, x*y^2])
            sage: P = PS(0,1)
            sage: f(P, check=False)
            (0 : 0)

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: H = Hom(X, X)
            sage: f = H([x^2, y^2, z^2])
            sage: f([4,4,1])
            (16 : 16 : 1)

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: H = Hom(X, X)
            sage: f = H([x^2, y^2, z^2])
            sage: f(P([4,4,1]))
            (16 : 16 : 1)
        """
        # Checks were done in __call__
        P = [f(x._coords) for f in self.defining_polynomials()]
        return self._codomain.point(P,check=True)

    def _call_with_args(self, x, args, kwds):
        """
        Apply this morphism to a point in the domain, with additional arguments.

        INPUT:

        - ``x`` -- a point in the domain or a list or tuple that defines a
          point in the domain
        - ``check`` -- boolean; either provided by position or name

        OUTPUT: a point in the codomain

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: H = A.Hom(A)
            sage: f = H([y, x^2 + y])
            sage: f([2,3])
            (3, 7)

        An example with algebraic schemes::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: X = A.subscheme(x)
            sage: Y = A.subscheme(y)
            sage: Hom_XY = X.Hom(Y)
            sage: f = Hom_XY([y,0])   # (0,y) |-> (y,0)
            sage: f
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field
                    defined by: x
              To:   Closed subscheme of Affine Space of dimension 2 over Rational Field
                    defined by: y
              Defn: Defined on coordinates by sending (x, y) to (y, 0)
            sage: f([0,3])
            (3, 0)

        As usual, if the input does not belong to a map's domain, it is first
        attempted to convert it::

            sage: f(0)
            Traceback (most recent call last):
            ...
            TypeError: 0 fails to convert into the map's domain Closed subscheme of
            Affine Space of dimension 2 over Rational Field defined by:
              x, but a `pushforward` method is not properly implemented

        It is possible to avoid the checks on the resulting point which can be
        useful for indeterminacies, but be careful!! ::

            sage: PS.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(PS,PS)
            sage: f = H([x^3, x*y^2])
            sage: P = PS(0,1)
            sage: f(P, check=False)     # indirect doctest
            (0 : 0)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: H = Hom(X, X)
            sage: f = H([x^2, y^2, z^2])
            sage: f([4,4,1])
            (16 : 16 : 1)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: P2.<u,v,w,t> = ProjectiveSpace(ZZ, 3)
            sage: X = P.subscheme(x^2 - y^2)
            sage: H = Hom(X, X)
            sage: f = H([x^2, y^2, z^2])
            sage: f(P2([4,4,1,1]))
            Traceback (most recent call last):
            ...
            TypeError: (4 : 4 : 1 : 1) fails to convert into the map's domain Closed subscheme of
            Projective Space of dimension 2 over Integer Ring defined by:
              x^2 - y^2, but a `pushforward` method is not properly implemented
        """
        if args:
            check = args[0]
        else:
            check = kwds.get("check", False)
        # containment of x in the domain has already been checked, in __call__
        P = [f(x._coords) for f in self.defining_polynomials()]
        return self._codomain.point(P,check)

    def _repr_defn(self):
        """
        Return a string representation of the definition of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: H = A.Hom(A)
            sage: f = H([y, x^2 + y])
            sage: print(f._repr_defn())
            Defined on coordinates by sending (x, y) to (y, x^2 + y)
        """
        i = self.domain().ambient_space()._repr_generic_point()
        o = self._codomain.ambient_space()._repr_generic_point(self.defining_polynomials())
        return "Defined on coordinates by sending %s to\n%s" % (i,o)

    def __getitem__(self, i):
        """
        Return the i-th poly with ``self[i]``.

        INPUT:

        - ``i`` -- integer

        OUTPUT: element of the coordinate ring of the domain

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([3/5*x^2, 6*y^2])
            sage: f[1]
            6*y^2
        """
        return self._polys[i]

    def __copy__(self):
        r"""
        Return a copy of ``self``.

        OUTPUT: :class:`SchemeMorphism_polynomial`

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([3/5*x^2, 6*y^2])
            sage: g = copy(f)
            sage: f == g
            True
            sage: f is g
            False

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: Q = X(23, 23, 46)
            sage: P = X(1, 1, 1)
            sage: P != Q
            True
        """
        return self.parent()(self._polys)

    def coordinate_ring(self):
        r"""
        Return the coordinate ring of the ambient projective space.

        OUTPUT: a multivariable polynomial ring over the base ring

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([3/5*x^2, 6*y^2])
            sage: f.coordinate_ring()
            Multivariate Polynomial Ring in x, y over Rational Field

        ::

            sage: R.<t> = PolynomialRing(ZZ, 1)
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: H = Hom(P, P)
            sage: f = H([3*x^2, y^2])
            sage: f.coordinate_ring()
            Multivariate Polynomial Ring in x, y over Multivariate Polynomial Ring
            in t over Integer Ring
        """
        return self._polys[0].parent()

    def change_ring(self, R, check=True):
        r"""
        Return a new :class:`SchemeMorphism_polynomial` which is this map
        coerced to ``R``.

        If ``check`` is ``True``, then the initialization checks are performed.

        INPUT:

        - ``R`` -- ring or morphism

        - ``check`` -- boolean

        OUTPUT: a new :class:`SchemeMorphism_polynomial` which is this map
        coerced to ``R``

        TESTS::

            sage: # needs sage.rings.number_field
            sage: R.<t> = QQ[]
            sage: K.<v> = QuadraticField(2)
            sage: K2.<w> = NumberField(t**4 - 2)
            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: phi = K.embeddings(K2)[0]
            sage: f = DynamicalSystem_projective([x**2 + 3*y**2, y**2])                 # needs sage.schemes
            sage: f.change_ring(phi)                                                    # needs sage.schemes
            Dynamical System of Projective Space of dimension 1 over
             Number Field in w with defining polynomial t^4 - 2
              Defn: Defined on coordinates by sending (x : y) to (x^2 + 3*y^2 : y^2)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<t> = QQ[]
            sage: K.<u> = QuadraticField(2)
            sage: K1.<v> = NumberField(t^4 - 2)
            sage: K2.<w> = NumberField(t^8 - 2)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: phi = K1.embeddings(K2)[0]
            sage: f = DynamicalSystem_projective([x^2 + 3*y^2, y^2])                    # needs sage.schemes
            sage: f.change_ring(phi)                                                    # needs sage.schemes
            Traceback (most recent call last):
            ...
            ValueError: no canonical coercion of base ring of morphism to domain of embedding

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: H = Hom(P, P)
            sage: f = H([3*x^2, y^2])
            sage: f.change_ring(GF(3))
            Scheme endomorphism of Projective Space of dimension 1 over Finite Field of size 3
              Defn: Defined on coordinates by sending (x : y) to (0 : y^2)

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P, P)
            sage: f = H([5/2*x^3 + 3*x*y^2 - y^3, 3*z^3 + y*x^2, x^3 - z^3])
            sage: f.change_ring(GF(3))
            Scheme endomorphism of Projective Space of dimension 2 over Finite Field of size 3
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^3 - y^3 : x^2*y : x^3 - z^3)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: X = P.subscheme([5*x^2 - y^2])
            sage: H = Hom(X, X)
            sage: f = H([x, y])
            sage: f.change_ring(GF(3))
            Scheme endomorphism of Closed subscheme of Projective Space of dimension 1
            over Finite Field of size 3 defined by: -x^2 - y^2
              Defn: Defined on coordinates by sending (x : y) to (x : y)


        Check that :issue:`16834` is fixed::

            sage: # needs sage.rings.real_mpfr
            sage: A.<x,y,z> = AffineSpace(RR, 3)
            sage: h = Hom(A, A)
            sage: f = h([x^2 + 1.5, y^3, z^5 - 2.0])
            sage: f.change_ring(CC)
            Scheme endomorphism of Affine Space of dimension 3 over
             Complex Field with 53 bits of precision
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x^2 + 1.50000000000000, y^3, z^5 - 2.00000000000000)

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: B.<u,v> = ProjectiveSpace(QQ, 1)
            sage: h = Hom(A,B)
            sage: f = h([x^2, y^2])
            sage: f.change_ring(QQ)
            Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to (x^2 : y^2)

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: H = Hom(A, A)
            sage: f = H([3*x^2/y, y^2/x])
            sage: f.change_ring(RR)                                                     # needs sage.rings.real_mpfr
            Scheme endomorphism of Affine Space of dimension 2 over Real Field with
            53 bits of precision
              Defn: Defined on coordinates by sending (x, y) to
                    (3.00000000000000*x^2/y, y^2/x)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^3 - x + 1)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = End(P)
            sage: f = H([x^2 + a*x*y + a^2*y^2, y^2])
            sage: emb = K.embeddings(QQbar)
            sage: f.change_ring(emb[0])
            Scheme endomorphism of Projective Space of dimension 1
             over Algebraic Field
               Defn: Defined on coordinates by sending (x : y) to
                     (x^2 + (-1.324717957244746?)*x*y + 1.754877666246693?*y^2 : y^2)
            sage: f.change_ring(emb[1])
            Scheme endomorphism of Projective Space of dimension 1
             over Algebraic Field
               Defn: Defined on coordinates by sending (x : y) to
                     (x^2 + (0.6623589786223730? - 0.5622795120623013?*I)*x*y
                      + (0.1225611668766537? - 0.744861766619745?*I)*y^2 : y^2)

        ::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: K.<v> = QuadraticField(2, embedding=QQbar(sqrt(2)))
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = End(P)
            sage: f = H([x^2 + v*y^2, y^2])
            sage: f.change_ring(QQbar)
            Scheme endomorphism of Projective Space of dimension 1
             over Algebraic Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + 1.414213562373095?*y^2 : y^2)

        ::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: from sage.misc.verbose import set_verbose
            sage: set_verbose(-1)
            sage: K.<w> = QuadraticField(2, embedding=QQbar(-sqrt(2)))
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: X = P.subscheme(x - y)
            sage: H = End(X)
            sage: f = H([6*x^2 + 2*x*y + 16*y^2, -w*x^2 - 4*x*y - 4*y^2])
            sage: f.change_ring(QQbar)
            Scheme endomorphism of Closed subscheme of Projective Space of dimension 1
             over Algebraic Field defined by: x - y
              Defn: Defined on coordinates by sending (x : y) to
                    (6*x^2 + 2*x*y + 16*y^2 : 1.414213562373095?*x^2 + (-4)*x*y + (-4)*y^2)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: f = x^6 - 2
            sage: L.<b> = NumberField(f, embedding=f.roots(QQbar)[1][0])
            sage: A.<x,y> = AffineSpace(L, 2)
            sage: H = Hom(A, A)
            sage: F = H([b*x/y, 1 + y])
            sage: F.change_ring(QQbar)
            Scheme endomorphism of Affine Space of dimension 2 over Algebraic Field
              Defn: Defined on coordinates by sending (x, y) to
                    (1.122462048309373?*x/y, y + 1)

        ::

            sage: # needs sage.rings.number_field
            sage: K.<a> = QuadraticField(-1)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: H = End(A)
            sage: phi = H([x/y, y])
            sage: emb = K.embeddings(QQbar)[0]
            sage: phi.change_ring(emb)
            Scheme endomorphism of Affine Space of dimension 2 over Algebraic Field
              Defn: Defined on coordinates by sending (x, y) to (x/y, y)
        """
        T = self.domain().change_ring(R)
        if self.is_endomorphism():
            H = End(T)
        else:
            S = self.codomain().change_ring(R)
            H = Hom(T, S)

        if isinstance(R, Map):
            from sage.structure.coerce_maps import CallableConvertMap
            if R.domain() == self.base_ring():
                S = self.domain().ambient_space().coordinate_ring()
                T = T.ambient_space().coordinate_ring()
                phi = CallableConvertMap(S, T, lambda self, g:T(g.map_coefficients(R)))
                G = []
                for f in self:
                    if isinstance(f, FractionFieldElement):
                        G.append(phi(f.numerator())/phi(f.denominator()))
                    else:
                        G.append(phi(f))
            elif R.domain().coerce_map_from(self.base_ring()) is not None:
                R = R*R.domain().coerce_map_from(self.base_ring())
                S = self.domain().ambient_space().coordinate_ring()
                T = T.ambient_space().coordinate_ring()
                phi = CallableConvertMap(S, T, lambda self, g:T(g.map_coefficients(R)))
                G = []
                for f in self:
                    if isinstance(f, FractionFieldElement):
                        G.append(phi(f.numerator())/phi(f.denominator()))
                    else:
                        G.append(phi(f))
            else:
                raise ValueError("no canonical coercion of base ring of morphism to domain of embedding")

        else:
            G = []
            for f in self:
                if isinstance(f, FractionFieldElement):
                    G.append(f.numerator().change_ring(R) / f.denominator().change_ring(R))
                else:
                    G.append(f.change_ring(R))
        return H(G, check)

    def specialization(self, D=None, phi=None, homset=None):
        r"""
        Specialization of this map.

        Given a family of maps defined over a polynomial ring. A specialization
        is a particular member of that family. The specialization can be specified either
        by a dictionary or a :class:`SpecializationMorphism`.

        INPUT:

        - ``D`` -- dictionary (optional)

        - ``phi`` -- SpecializationMorphism (optional)

        - ``homset`` -- homset of specialized map (optional)

        OUTPUT: :class:`SchemeMorphism_polynomial`

        EXAMPLES::

            sage: R.<c> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: H = End(P)
            sage: f = H([x^2 + c*y^2, y^2])
            sage: f.specialization({c: 1})
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to (x^2 + y^2 : y^2)

        ::

            sage: R.<a,b> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: H = End(P)
            sage: f = H([x^3 + a*x*y^2 + b*y^3, y^3])
            sage: from sage.rings.polynomial.flatten import SpecializationMorphism
            sage: phi = SpecializationMorphism(P.coordinate_ring(), {a: 2, b: -1})
            sage: F = f.specialization(phi=phi); F
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^3 + 2*x*y^2 - y^3 : y^3)
            sage: g = H([x^2 + a*y^2, y^2])
            sage: G = g.specialization(phi=phi)
            sage: G.parent() is F.parent()
            True
            sage: G = g.specialization(phi=phi, homset=F.parent())
            sage: G.parent() is F.parent()
            True

        ::

            sage: R.<c> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: X = P.subscheme([x - c*y])
            sage: H = End(X)
            sage: f = H([x^2, c*y^2])
            sage: f.specialization({c: 2})
            Scheme endomorphism of Closed subscheme of Projective Space of dimension 1
             over Rational Field defined by: x - 2*y
                  Defn: Defined on coordinates by sending (x : y) to (x^2 : 2*y^2)

        ::

            sage: R.<c> = QQ[]
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: f = DynamicalSystem_projective([x^2 + c*y^2, y^2], domain=P)
            sage: F = f.dynatomic_polynomial(3)                                         # needs sage.libs.pari
            sage: g = F.specialization({c: 1}); g
            x^6 + x^5*y + 4*x^4*y^2 + 3*x^3*y^3 + 7*x^2*y^4 + 4*x*y^5 + 5*y^6
            sage: g == f.specialization({c:1}).dynatomic_polynomial(3)                  # needs sage.libs.pari
            True

        ::

            sage: R1.<alpha, beta> = QQ[]
            sage: A.<x> = AffineSpace(Frac(R1), 1)
            sage: f = DynamicalSystem_affine([alpha/(x^2 + 1/alpha)/(x - 1/beta^2)])
            sage: f.specialization({alpha: 5, beta: 10})
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                      (5/(x^3 - 1/100*x^2 + 1/5*x - 1/500))
            sage: f_5_10 = f.specialization({alpha: 5}).specialization({beta: 10})
            sage: f_5_10 == f.specialization({alpha: 5, beta: 10})
            True
        """
        if D is None:
            if phi is None:
                raise ValueError("either the dictionary or the specialization must be provided")
        else:
            if isinstance(self[0].parent(), FractionField_generic):
                from sage.rings.polynomial.flatten import FractionSpecializationMorphism
                phi = FractionSpecializationMorphism(self[0].parent(), D)
            else:
                from sage.rings.polynomial.flatten import SpecializationMorphism
                phi = SpecializationMorphism(self[0].parent(), D)
        if homset is None:
            domain = self.domain()
            if isinstance(domain, AlgebraicScheme_subscheme):
                domain = domain.specialization(phi=phi)
            else:
                domain = domain.change_ring(phi.codomain().base_ring())
            if self.is_endomorphism():
                homset = End(domain)
            else:
                codomain = self.codomain()
                if isinstance(codomain, AlgebraicScheme_subscheme):
                    codomain = codomain.specialization(phi=phi)
                else:
                    codomain = codomain.change_ring(phi.codomain().base_ring())
                homset = Hom(domain, codomain)
        return homset([phi(g) for g in self])

    def _composition_(self, other, homset):
        r"""
        Straightforward implementation of composition for scheme morphisms
        defined by polynomials.

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(P, P)
            sage: f = H([x^2 - 29/16*y^2, y^2])
            sage: g = H([y, x + y])
            sage: h = f*g
            sage: h
            Scheme endomorphism of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (-29/16*x^2 - 29/8*x*y - 13/16*y^2 : x^2 + 2*x*y + y^2)
            sage: p = P((1,3))
            sage: h(p) == f(g(p))
            True

            sage: Q = ProjectiveSpace(QQ, 2)
            sage: H2 = Hom(P, Q)
            sage: h2 = H2([x^2 + y^2, x^2, y^2 + 2*x^2])
            sage: h2 * f
            Scheme morphism:
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^4 - 29/8*x^2*y^2 + 1097/256*y^4
                     : x^4 - 29/8*x^2*y^2 + 841/256*y^4
                     : 2*x^4 - 29/4*x^2*y^2 + 969/128*y^4)

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: A1.<z> = AffineSpace(QQ, 1)
            sage: H = End(A)
            sage: f = H([x^2 + y^2, y^2/x])
            sage: H1 = Hom(A, A1)
            sage: g = H1([x + y^2])
            sage: g*f
            Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    ((x^4 + x^2*y^2 + y^4)/x^2)
            sage: f*g
            Traceback (most recent call last):
            ...
            TypeError: self (=Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    (x^2 + y^2, y^2/x)) domain
            must equal right (=Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    (y^2 + x)) codomain

        Not both defined by polynomials::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2 - 2)
            sage: p1, p2 = K.Hom(K)
            sage: R.<x,y> = K[]
            sage: q1 = R.Hom(R)(p1)
            sage: A = AffineSpace(R)
            sage: f1 = A.Hom(A)(q1)
            sage: g = A.Hom(A)([x^2 - y, y + 1])
            sage: g*f1
            Composite map:
              From: Affine Space of dimension 2 over Number Field in a
                    with defining polynomial x^2 - 2
              To:   Affine Space of dimension 2 over Number Field in a
                    with defining polynomial x^2 - 2
              Defn:   Generic endomorphism of Affine Space of dimension 2
                       over Number Field in a with defining polynomial x^2 - 2
                    then
                      Generic endomorphism of Affine Space of dimension 2
                       over Number Field in a with defining polynomial x^2 - 2
        """
        try:
            opolys = tuple(other._polys)
        except AttributeError:
            return super()._composition_(other, homset)
        return homset([p(*opolys) for p in self._polys])


class SchemeMorphism_polynomial_id(SchemeMorphism_id, SchemeMorphism_polynomial):
    """
    Return the identity morphism from `X` to itself.

    INPUT:

    - ``X`` -- an affine or projective scheme

    EXAMPLES::

        sage: X = Spec(ZZ)
        sage: X.identity_morphism()  # indirect doctest
        Scheme endomorphism of Spectrum of Integer Ring
          Defn: Identity map
    """
    def __init__(self, X):
        """
        Initialize.

        TESTS::

            sage: A = AffineSpace(2, GF(3))
            sage: A.identity_morphism().defining_polynomials()
            (x0, x1)
        """
        super().__init__(X)
        variables = X.ambient_space().coordinate_ring().gens()
        SchemeMorphism_polynomial.__init__(self, X.Hom(X), variables, check=False)


############################################################################
# Rational points on schemes, which we view as morphisms determined
# by coordinates.
############################################################################

class SchemeMorphism_point(SchemeMorphism):
    r"""
    Base class for rational points on schemes.

    Recall that the `K`-rational points of a scheme `X` over `k` can
    be identified with the set of morphisms `\mathrm{Spec}(K) \to X`. In Sage,
    the rational points are implemented by such scheme morphisms.

    EXAMPLES::

        sage: from sage.schemes.generic.morphism import SchemeMorphism
        sage: f = SchemeMorphism(Spec(ZZ).Hom(Spec(ZZ)))
        sage: type(f)
        <class 'sage.schemes.generic.morphism.SchemeMorphism'>
    """
    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: a._repr_()
            '(1, 2)'
        """
        return self._codomain.ambient_space()._repr_generic_point(self._coords)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: latex(a) == a._latex_()
            True
            sage: a._latex_()
            '\\left(1, 2\\right)'
        """
        return self._codomain.ambient_space()._latex_generic_point(self._coords)

    def __getitem__(self, n):
        """
        Return the ``n``-th coordinate.

        OUTPUT: the coordinate values as an element of the base ring

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: a[0]
            1
            sage: a[1]
            2
        """
        return self._coords[n]

    def __iter__(self):
        """
        Iterate over the coordinates of the point.

        OUTPUT: an iterator

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: iter = a.__iter__()
            sage: next(iter)
            1
            sage: next(iter)
            2
            sage: list(a)
            [1, 2]
        """
        return iter(self._coords)

    def __tuple__(self):
        """
        Return the coordinates as a tuple.

        OUTPUT: a tuple

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: tuple(a)
            (1, 2)
        """
        return self._coords

    def __len__(self):
        """
        Return the number of coordinates.

        OUTPUT: integer. The number of coordinates used to describe the point

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: len(a)
            2
        """
        return len(self._coords)

    def _richcmp_(self, other, op):
        """
        Compare two scheme morphisms.

        INPUT:

        - ``other`` -- anything; to compare against the scheme
          morphism ``self``

        OUTPUT: boolean

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: b = A(3,4)
            sage: a < b
            True
            sage: a != b
            True
        """
        if not isinstance(other, SchemeMorphism_point):
            try:
                other = self._codomain.ambient_space()(other)
            except TypeError:
                return NotImplemented
        return richcmp(self._coords, other._coords, op)

    def scheme(self):
        """
        Return the scheme whose point is represented.

        OUTPUT: a scheme

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: a.scheme()
            Affine Space of dimension 2 over Rational Field
        """
        return self._codomain

    def change_ring(self, R, check=True):
        r"""
        Return a new :class:`SchemeMorphism_point` which is this point coerced
        to ``R``.

        If ``check`` is true, then the initialization checks are performed.

        INPUT:

        - ``R`` -- ring or morphism

        - ``check`` -- boolean

        OUTPUT: :class:`SchemeMorphism_point`

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: X = P.subscheme(x^2 - y^2)
            sage: X(23,23,1).change_ring(GF(13))
            (10 : 10 : 1)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P(-2/3,1).change_ring(CC)                                             # needs sage.rings.real_mpfr
            (-0.666666666666667 : 1.00000000000000)

        ::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: P(152,113).change_ring(Zp(5))                                         # needs sage.rings.padics
            (2 + 5^2 + 5^3 + O(5^20) : 3 + 2*5 + 4*5^2 + O(5^20))

        ::

            sage: # needs sage.rings.number_field
            sage: K.<v> = QuadraticField(-7)
            sage: O = K.maximal_order()
            sage: P.<x,y> = ProjectiveSpace(O, 1)
            sage: H = End(P)
            sage: F = H([x^2 + O(v)*y^2, y^2])
            sage: F.change_ring(K).change_ring(K.embeddings(QQbar)[0])
            Scheme endomorphism of Projective Space of dimension 1
             over Algebraic Field
              Defn: Defined on coordinates by sending (x : y) to
                    (x^2 + (-2.645751311064591?*I)*y^2 : y^2)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^2 - x + 1)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: Q = P([a + 1, 1])
            sage: emb = K.embeddings(QQbar)
            sage: Q.change_ring(emb[0])
            (1.5000000000000000? - 0.866025403784439?*I : 1)
            sage: Q.change_ring(emb[1])
            (1.5000000000000000? + 0.866025403784439?*I : 1)

        ::

            sage: # needs sage.rings.number_field
            sage: K.<v> = QuadraticField(2)
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: Q = P([v,1])
            sage: Q.change_ring(K.embeddings(QQbar)[0])
            (-1.414213562373095? : 1)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: f = x^6 - 2
            sage: L.<b> = NumberField(f, embedding=f.roots(QQbar)[1][0])
            sage: A.<x,y> = AffineSpace(L, 2)
            sage: P = A([b,1])
            sage: P.change_ring(QQbar)
            (1.122462048309373?, 1)
        """
        S = self.codomain().change_ring(R)
        Q = [R(t) for t in self]
        return S.point(Q, check=check)

    def __copy__(self):
        r"""
        Return a copy of the :class:`SchemeMorphism_point` ``self`` coerced to
        `R`.

        OUTPUT: :class:`SchemeMorphism_point`

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: Q = P(152, 113)
            sage: Q2 = copy(Q)
            sage: Q2 is Q
            False
            sage: Q2 == Q
            True
        """
        return self._codomain.point(self._coords, check=False)

    def specialization(self, D=None, phi=None, ambient=None):
        r"""
        Specialization of this point.

        Given a family of points defined over a polynomial ring. A specialization
        is a particular member of that family. The specialization can be specified either
        by a dictionary or a :class:`SpecializationMorphism`.

        INPUT:

        - ``D`` -- dictionary (optional)

        - ``phi`` -- SpecializationMorphism (optional)

        - ``ambient`` -- ambient space of specialized point (optional)

        OUTPUT: :class:`SchemeMorphism_polynomial`

        EXAMPLES::

            sage: R.<c> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: Q = P([c,1])
            sage: Q.specialization({c: 1})
            (1 : 1)

        ::

            sage: R.<a,b> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: Q = P([a^2 + 2*a*b + 34, 1])
            sage: from sage.rings.polynomial.flatten import SpecializationMorphism
            sage: phi = SpecializationMorphism(P.coordinate_ring(), {a: 2, b: -1})
            sage: T = Q.specialization(phi=phi); T
            (34 : 1)
            sage: Q2 = P([a,1])
            sage: T2 = Q2.specialization(phi=phi)
            sage: T2.codomain() is T.codomain()
            True
            sage: T3 = Q2.specialization(phi=phi, ambient=T.codomain())
            sage: T3.codomain() is T.codomain()
            True

        ::

            sage: R.<c> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: X = P.subscheme([x - c*y])
            sage: Q = X([c, 1])
            sage: Q2 = Q.specialization({c:2}); Q2
            (2 : 1)
            sage: Q2.codomain()
            Closed subscheme of Projective Space of dimension 1 over Rational Field
             defined by: x - 2*y

        ::

            sage: R.<l> = PolynomialRing(QQ)
            sage: S.<k,j> = PolynomialRing(R)
            sage: K.<a,b,c,d> = S[]
            sage: P.<x,y> = ProjectiveSpace(K, 1)
            sage: H = End(P)
            sage: Q = P([a^2, b^2])
            sage: Q.specialization({a: 2})
            (4 : b^2)
        """
        if D is None:
            if phi is None:
                raise ValueError("either the dictionary or the specialization must be provided")
        else:
            from sage.rings.polynomial.flatten import SpecializationMorphism
            phi = SpecializationMorphism(self.codomain().ambient_space().coordinate_ring(), D)
        if ambient is None:
            ambient = self.codomain()
            if isinstance(ambient, AlgebraicScheme_subscheme):
                ambient = ambient.specialization(phi=phi)
            else:
                ambient = ambient.change_ring(phi.codomain().base_ring())
        psi = ambient.ambient_space().coordinate_ring().hom([0 for i in range(ambient.ambient_space().ngens())], ambient.base_ring())
        return ambient([psi(phi(t)) for t in self])
