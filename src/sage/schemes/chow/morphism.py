# -*- coding: utf-8 -*-
r"""
Morphisms between ChowSchemes

A morphism between ChowSchemes is given by a ring map between its Chow rings.
For example, the Veronese embedding is given as follows::

    sage: X.<x> = Proj(2, 'x', name='X')
    sage: Y.<y> = Proj(5, 'y', name='Y')
    sage: f = X.hom([2*x], Y); f
    ChowScheme morphism:
      From: X
      To:   Y
      Defn: y |--> 2*x

The domain and codomain are retrieved as usual::

    sage: f.domain() == X
    True
    sage: f.codomain() == Y
    True

Compute the images of elements `y` in the Chow ring of the codomain, either
by ``f.upperstar(y)`` or simply ``f(y)``::

    sage: f.upperstar(y)
    2*x
    sage: f(y)
    2*x
    sage: f(y^2)
    4*x^2
    sage: f(2*y^3)
    0

The composition of two morphisms is given by applying the * operator::

    sage: Z.<z> = Proj(56, 'z', name='Z')
    sage: g = Y.hom([3*y], Z)
    sage: h = g * f
    sage: k = X.hom([6*x], Z)  # Direct definition of the composition
    sage: h == k
    True

If the point_classes of the domain and codomain are defined, e.g. domain and
codomain are representing projective non singular varieties, as they  are for
:class:`library.proj.Proj`, one can compute ``f.lowerstar(x)`` for elements `x`
in the Chow ring of the domain::

    sage: f.lowerstar(x)
    2*y^4
    sage: f.lowerstar(x^2)
    y^5
    sage: f.lowerstar(1 + x^2)
    y^5 + 4*y^3

If E and F are sheaves on X and Y respectively, compute the pullback and
pushforward as follows::

    sage: TY = Y.tangent_bundle()
    sage: f.upperstar(TY)
    Bundle(X, 5, [1, 12*x, 60*x^2])
    sage: TX = X.tangent_bundle()
    sage: f.lowerstar(TX)
    Sheaf(Y, 0, [1, 0, 0, 16*y^3, 72*y^4, 240*y^5])

TESTS::

    sage: TestSuite(f).run(skip=["_test_category"])

AUTHORS:

- Manfred Lehn (2013)
- Christoph Sorger (2013)
"""

# ****************************************************************************
#       Copyright (C) 2013 Manfred Lehn <lehn@mathematik.uni-mainz.de>
#       Copyright (C) 2013 Christoph Sorger <christoph.sorger@univ-nantes.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.all import Rings
_Rings = Rings() # type: ignore
from sage.categories.map import Map
from sage.structure.element import Element
from sage.arith.power import generic_power


def is_chowSchemeMorphism(f):
    """
    Test whether ``f`` is a chowscheme morphism.

    INPUT:

    - ``f`` -- anything.

    OUTPUT:

    Boolean. Return ``True`` if ``f`` is a chowscheme morphism.

    EXAMPLES::

        sage: P1.<h> = ChowScheme(1, 'h', 1, 'h^2', 'h', name='P1')
        sage: P3.<k> = ChowScheme(3, 'k', 1, 'k^4', 'k^3', name='P3')
        sage: H = P1.Hom(P3)
        sage: f = H([3*h])
        sage: is_chowSchemeMorphism(f)
        True
    """
    # from chow import morphism
    # return isinstance(f, morphism.ChowSchemeMorphism)
    return isinstance(f, ChowSchemeMorphism)


class ChowSchemeMorphism(Element):
    r"""
    Base class for morphisms between ChowSchemes.
    """
    def __init__(self, parent, phi):
        """
        Construct a :class:`ChowSchemeMorphism`

        INPUT:

        - ``parent`` -- the parent of the morphism.
        - ``phi`` -- a ring morphism from f.codomain().chowring() to
          f.domain().chowring() or a list of elements representing such
          a ring morphism.

        OUTPUT:

        - :class:`ChowRing <ChowRing_generic>`.

        EXAMPLES::

            sage: from sage.schemes.chow.morphism import ChowSchemeMorphism
            sage: X = ChowScheme(1, 'h', 1, 'h^2')
            sage: Hom = X.Hom(X)
            sage: f = ChowSchemeMorphism(Hom, ['h'])
        """
        from sage.categories.homset import Homset
        if not isinstance(parent, Homset):
            raise TypeError("parent (=%s) must be a hom-space" % parent)
        Element.__init__(self, parent)
        X, Y = parent.domain(), parent.codomain()
        AX, AY = X.chowring(), Y.chowring()
        self._domain, self._codomain = X, Y
        self._category = parent.category()
        self._is_endomorphism_set = parent.is_endomorphism_set()
        if isinstance(phi, Map) and phi.category_for().is_subcategory(_Rings):
            msg = "Expect a ring morphism from %s to %s" % (str(AY), str(AX))
            if not (phi.domain() == AY and phi.codomain() == AX):
                raise TypeError(msg)
        else:
            phi = phi if type(phi) is list or tuple else [phi]
            phi = [AX(x) for x in phi]
            phi = AY.hom(phi, AX)
        self._chowring_morphism = phi

    def __call__(self, y):
        r"""
        Return the image of the element ``y`` of the Chow ring of the codomain
        in the Chow ring of the domain under this ChowSchemeMorphism

        EXAMPLES::

            sage: P1.<h> = ChowScheme(1, 'h', 1, 'h^2', 'h', name='P1')
            sage: P3.<k> = ChowScheme(3, 'k', 1, 'k^4', 'k^3', name='P3')
            sage: f = P1.hom([3*h], P3)
            sage: f(3*k^2+4*k)
            12*h
        """
        return self._chowring_morphism(y)

    def __eq__(self, other):
        """
        Return True if this ChowSchemeMorphism and other have isomorphic
        domain, codomain and Chow ring morphism.

        OUTPUT:

        - Boolean

        EXAMPLES::

            sage: P1.<h> = ChowScheme(1, 'h', 1, 'h^2', 'h', name='P1')
            sage: P3.<k> = ChowScheme(3, 'k', 1, 'k^4', 'k^3', name='P3')
            sage: f = P1.hom([3*h], P3)
            sage: g = P1.hom([5*h^2 + 3*h], P3)
            sage: f == g
            True

        TESTS::

            sage: P1.<h> = ChowScheme(1, 'h', 1, 'h^2', 'h', name='P1')
            sage: P3.<k> = ChowScheme(3, 'k', 1, 'k^4', 'k^3', name='P3')
            sage: f = P1.hom([3*h], P3)
            sage: g = loads(f.dumps())
            sage: f == g
            True
        """
        if other is None:
            return False
        if self.domain() == other.domain() and \
                self.codomain() == other.codomain():
            if self.chowring_morphism() != other.chowring_morphism():
                return False
            return True
        return False

    def __ne__(self, other):
        r"""
        Return True if this ChowSchemesMorphism and other are not equal.

        OUTPUT:

        - Boolean

        EXAMPLES::

            sage: P1.<h> = ChowScheme(1, 'h', 1, 'h^2', 'h', name='P1')
            sage: P3.<k> = ChowScheme(3, 'k', 1, 'k^4', 'k^3', name='P3')
            sage: f = P1.hom([3*h], P3)
            sage: g = P1.hom([2*h], P3)
            sage: f != g
            True
        """
        return not self.__eq__(other)

    def __mul__(self, other):
        """
        Return the composition of this ChowSchemeMorphism with another
        ChowSchemeMorphism.

        OUTPUT:

        - A ChowSchemeMorphism

        EXAMPLES::

            sage: X.<h> = ChowScheme(2, 'h', 1, 'h^3', 'h^2')
            sage: Y.<k> = ChowScheme(5, 'k', 1, 'k^6', 'k^5')
            sage: Z.<l> = ChowScheme(56, 'l', 1, 'l^57', 'l^56')
            sage: f = X.hom([2*h], Y)
            sage: g = Y.hom([3*k], Z)
            sage: h = X.hom([6*h], Z)  # Direct definition of the composition
            sage: h == (g * f)
            True

        TESTS::

            sage: f * g
            Traceback (most recent call last):
            ...
            TypeError: The codomain of first map and the domain of second differ.
        """
        g, f = self, other
        X, Z = f.domain(), g.codomain()
        AZ = Z.chowring()
        if not f.codomain() is g.domain():
            msg = "The codomain of first map and the domain of second differ."
            raise TypeError(msg)
        return X.hom([other(self(x)) for x in AZ.gens()], Z)

    def __pow__(self, n):
        r"""
        If this ChowSchemeMorphism is an endomorphism return its `n`-th power.

        OUTPUT:

        - A ChowSchemeMorphism

        EXAMPLES::

            sage: P2.<h> = ChowScheme(2, 'h', 1, 'h^3', name='P2')
            sage: f = P2.hom([h], P2)  # Identity
            sage: f^2 == f
            True
            sage: g = f^0; g
            ChowScheme endomorphism of P2
              Defn: Identity map
            sage: g(h)
            h
            sage: f = P2.hom([2*h], P2)  # Identity
            sage: f^2 == f
            False
            sage: (f^3)(h)
            8*h

        TESTS::

            sage: X.<h> = ChowScheme(2, 'h', 1, 'h^3', 'h^2')
            sage: Y.<k> = ChowScheme(5, 'k', 1, 'k^6', 'k^5')
            sage: f = X.hom([2*h], Y)
            sage: f^3
            Traceback (most recent call last):
            ...
            TypeError: Can take powers only of endomorphisms.
        """
        if not self.is_endomorphism():
            raise TypeError("Can take powers only of endomorphisms.")
        if n == 0:
            return self.domain().identity_morphism()
        return generic_power(self, n)

    def _repr_defn(self):
        r"""
        Return a string representation of the definition of this
        ChowSchemeMorphism.

        OUTPUT:

        - a string.

        TESTS::

            sage: X.<h> = ChowScheme(2, 'h', 1, 'h^3', 'h^2', name='X')
            sage: Y.<k> = ChowScheme(5, 'k', 1, 'k^6', 'k^5', name='Y')
            sage: f = X.hom([2*h], Y)
            sage: f._repr_defn()
            'k |--> 2*h'
        """
        D = self._chowring_morphism.domain().gens()
        Im = self._chowring_morphism.morphism_from_cover().im_gens()
        return '\n'.join(['%s |--> %s' % (D[i], Im[i]) for i in range(len(D))])

    def _repr_type(self):
        r"""
        Return a string representation of the type of this ChowSchemeMorphism.

        OUTPUT:

        - a string.

        TESTS::

            sage: X.<h> = ChowScheme(2, 'h', 1, 'h^3', 'h^2', name='X')
            sage: Y.<k> = ChowScheme(5, 'k', 1, 'k^6', 'k^5', name='Y')
            sage: f = X.hom([2*h], Y)
            sage: f._repr_type()
            'ChowScheme'
        """
        return 'ChowScheme'

    def _repr_(self):
        r"""
        Return a string representation of this ChowSchemeMorphism.

        OUTPUT:

        - a string.

        TESTS::

            sage: X.<h> = ChowScheme(2, 'h', 1, 'h^3', 'h^2', name='X')
            sage: Y.<k> = ChowScheme(5, 'k', 1, 'k^6', 'k^5', name='Y')
            sage: f = X.hom([2*h], Y)
            sage: f._repr_()
            'ChowScheme morphism:\n  From: X\n  To:   Y\n  Defn: k |--> 2*h'
        """
        if self.is_endomorphism():
            s = "%s endomorphism of %s" % (self._repr_type(), self.domain())
        else:
            s = "%s morphism:" % self._repr_type()
            s += "\n  From: %s" % self.domain()
            s += "\n  To:   %s" % self.codomain()
        d = self._repr_defn()
        if d != '':
            s += "\n  Defn: %s" % ('\n        '.join(d.split('\n')))
        return s

    def chowring_morphism(self):
        r"""
        Return the morphism on the level of Chow rings for this
        ChowSchemeMorphism.

        EXAMPLES::

            sage: P1.<h> = ChowScheme(1, 'h', 1, 'h^2', 'h', name='P1')
            sage: P3.<k> = ChowScheme(3, 'k', 1, 'k^4', 'k^3', name='P3')
            sage: f = P1.hom([3*h], P3)
            sage: f.chowring_morphism()
            Ring morphism:
              From: A(P3)
              To:   A(P1)
              Defn: k |--> 3*h
        """
        return self._chowring_morphism

    def domain(self):
        r"""
        Return the domain of this ChowSchemeMorphism.

        OUTPUT:

        A Chowscheme. The domain of this ChowSchemeMorphism

        EXAMPLES::

            sage: X.<w> = Proj(1, 'w')
            sage: Y.<h> = Proj(3, 'h')
            sage: i = X.hom(['3 * w'], Y)
            sage: i.domain()
            ChowScheme(1, 'w', 1, 'w^2', 'w')
        """
        return self._domain

    def codomain(self):
        """
        Return the codomain of this ChowSchemeMorphism.

        OUTPUT:

        A ChowsScheme. The codomain of this ChowSchemeMorphism.

        EXAMPLES::

            sage: X.<w> = Proj(1, 'w')
            sage: Y.<h> = Proj(3, 'h')
            sage: i = X.hom(['3 * w'], Y)
            sage: i.codomain()
            ChowScheme(3, 'h', 1, 'h^4', 'h^3')
        """
        return self._codomain

    def category(self):
        """
        Return the category of the SHom-set.

        OUTPUT:

        A category.

        EXAMPLES::

            sage: X.<w> = Proj(1, 'w')
            sage: Y.<h> = Proj(3, 'h')
            sage: i = X.hom(['3 * w'], Y)
            sage: i.category()
            Category of homsets of chowschemes over PointChowScheme
        """
        return self._category

    def is_endomorphism(self):
        """
        Return whether the morphism is an endomorphism.

        OUTPUT:

        Boolean. Return true if the domain and codomain are identical.

        EXAMPLES::

            sage: X.<w> = Proj(1, 'w')
            sage: Y.<h> = Proj(3, 'h')
            sage: i = X.hom(['3 * w'], Y)
            sage: i.is_endomorphism()
            False
            sage: j = X.hom(['2 * w'], X)
            sage: j.is_endomorphism()
            True
        """
        return self._is_endomorphism_set

    def upperstar(self, v):
        r"""
        Return `f^{*}(v)` where ``v`` can be a class in `A^{*}(Y)` or a
        sheaf on `Y`.

        INPUT:

        - ``v`` -- a class in the Chow ring of the codomain or a sheaf on the
          codomain.

        OUTPUT:

        A class in the Chow ring of the domain or a sheaf on the domain
        depending whether ``v`` is a class or a sheaf.

        EXAMPLES::

            sage: X.<w> = Proj(1, 'w', name='X')
            sage: Y.<h> = Proj(3, 'h', name='Y')
            sage: i = X.hom(['3 * w'], Y)
            sage: i.upperstar(h)
            3*w
            sage: i.upperstar(Y.o(1))
            Bundle(X, 1, [1, 3*w])
            sage: i.upperstar(Y.tangent_bundle())
            Bundle(X, 3, [1, 12*w])
        """
        from .sheaf import Sheaf, is_sheaf
        if is_sheaf(v):
            # Return a sheaf if the argument is a sheaf
            if not v.chowscheme() is self.codomain():
                raise ValueError("Expect a sheaf on codomain.")
            cc = v.chern_classes()
            upperstar_cc = self.upperstar(cc)
            if len(upperstar_cc) > self.domain().dimension():
                upperstar_cc = upperstar_cc[0:self.domain().dimension() + 1]
            result = Sheaf(self.domain(), v.rank(), upperstar_cc)
            # If called from a derived class like bundle return a bundle
            result.__class__ = v.__class__
            return result
        if type(v) is list:
            return [self.upperstar(y) for y in v]
        return self(v)

    def lowerstar(self, v, normal_bundle=None, verbose=False):
        r"""
        Return `f_{*}(v)` where `v` can be a class in `A^{*}(X)` or a
        sheaf on the domain `X`. Both, the domain `X` and the codomain `Y` are
        supposed to be projective, e.g. their point_classes respectively
        have to be known.

        In case, `v` is a sheaf, the tangent bundles on `X` and `Y` have to
        be known as Grothendieck-Hirzebruch-Riemann-Roch is
        used for the computation. Note, as indicated already above, the result
        should be understood as a virtual sheaf
        (e.g. the alternating sum of `R^{i}f_{*}v`).

        If the normal bundle is specified, `f` is supposed to be an embedding
        and Riemann-Roch without denominators is used for the computation. In
        this case no need for the tangent bundles on `X` and `Y` and the result
        is not virtual.

        INPUT:

        - ``v`` -- a class in the Chow ring of the domain or a sheaf on the
          domain.

        - ``normal_bundle`` -- the (optional) normal bundle of this
          ChowSchemeMorphism (supposed to be an embedding in this case).

        OUTPUT:

        A class in the Chow ring of the codomain or a sheaf on the codomain
        depending whether ``v`` is a class or a sheaf.

        EXAMPLES::

            sage: X.<w> = Proj(1, 'w', name='X')
            sage: Y.<h> = Proj(3, 'h', name='Y')
            sage: i = X.hom([3*w], Y)
            sage: i.lowerstar(w)
            h^3
            sage: i.lowerstar(X.o(1))
            Sheaf(Y, 0, [1, 0, -3*h^2, -8*h^3])
            sage: i.lowerstar(X.tangent_bundle())
            Sheaf(Y, 0, [1, 0, -3*h^2, -6*h^3])

            sage: P2.<h> = ChowScheme(2, 'h', 1, 'h^3', 'h^2', name='P2')
            sage: P5.<k> = ChowScheme(5, 'k', 1, 'k^6', 'k^5', name='P5')
            sage: f = P2.hom([2*h], P5)
            sage: E = Bundle(P2, 2, [1, 3*h, 3*h^2])
            sage: N = Bundle(P2, 3, [1, 9*h, 30*h^2])  # Normal bundle
            sage: f.lowerstar(E, normal_bundle=N)
            Sheaf(P5, 0, [1, 0, 0, 16*k^3, 72*k^4, 240*k^5])

            sage: P2 = Proj(2)
            sage: f = P2.base_morphism()  # P2 -> Pt
            sage: f.lowerstar(P2.o(1)).rank()  # = dim H^0(P2, O(1)) par GHRR
            3
        """
        X, Y = self.domain(), self.codomain()
        from .sheaf import Sheaf, is_sheaf
        if is_sheaf(v):
            # Return a sheaf if the argument is a sheaf
            if not v.chowscheme() is self.domain():
                raise ValueError("Expect a sheaf on domain.")
            if normal_bundle is None:
                # Use Grothendieck-Hirzebruch-Riemann-Roch
                TX, TY = X.tangent_bundle(), Y.tangent_bundle()
                Tf = TX - self.upperstar(TY)
                ch = self.lowerstar(v.chern_character() * Tf.todd_class(),
                                    verbose=verbose)
                return Sheaf(Y, ch=ch)
            else:
                # Normal bundle given, use Riemann-Roch without denominators
                cd = Y.dimension() - X.dimension()
                ls = self.lowerstar(normal_bundle.sans_denominateurs(v),
                                    verbose=verbose)
                lowerstar_cc = [1] + [0] * (cd - 1) + ls
                if len(lowerstar_cc) > Y.dimension():
                    lowerstar_cc = lowerstar_cc[0:Y.dimension() + 1]
                result = Sheaf(Y, 0, lowerstar_cc)
                return result
        if type(v) is list:
            return [self.lowerstar(x, verbose=verbose) for x in v]
        AX, AY = X.chowring(), Y.chowring()
        v = AX(v) if isinstance(v, str) else v
        rb, dbd = AY.basis(), AY.dual_basis_dict(verbose=verbose)
        return sum((self(dbd[e]) * v).integral() * e for e in rb)


class ChowSchemeMorphism_id(ChowSchemeMorphism):
    """
    Return the identity morphism from `X` to itself.

    INPUT:

    - ``X`` -- the ChowScheme.

    EXAMPLES::

            sage: from sage.schemes.chow.morphism import ChowSchemeMorphism_id
            sage: X = ChowScheme(1, 'h', 1, 'h^2', name='P1')
            sage: ChowSchemeMorphism_id(X)
            ChowScheme endomorphism of P1
              Defn: Identity map


    """
    def __init__(self, X):
        r"""
        Construct a :class:`ChowSchemeMorphism_id`.

        TESTS:
            sage: X = ChowScheme(1, 'h', 1, 'h^2', name='P1')
            sage: X.identity_morphism()
            ChowScheme endomorphism of P1
              Defn: Identity map
        """
        AX = X.chowring()
        f = AX.identity_morphism()
        ChowSchemeMorphism.__init__(self, X.Hom(X), f)

    def _repr_defn(self):
        r"""
        Return a string representation of the definition of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.schemes.chow.morphism import ChowSchemeMorphism_id
            sage: X = ChowScheme(1, 'h', 1, 'h^2', name='P1')
            sage: ChowSchemeMorphism_id(X)._repr_defn()
            'Identity map'
        """
        return 'Identity map'


class ChowSchemeMorphism_structure_map(ChowSchemeMorphism):
    r"""
    The structure morphism

    INPUT:

    - ``parent`` -- SHom-set with codomain equal to the base ChowScheme of
      the domain.
    """
    def __init__(self, parent, f):
        """
        Construct a :class:`ChowSchemeMorphism_structure_map`.

        EXAMPLES::

            sage: X = ChowScheme(1, 'h', 1, 'h^2', name='P1')
            sage: X.base_morphism()
            ChowScheme morphism:
              From: P1
              To:   PointChowScheme
              Defn: Structure map
        """
        ChowSchemeMorphism.__init__(self, parent, f)
        if self.domain().base_chowscheme() != self.codomain():
            msg = "parent must have codomain equal the base scheme of domain."
            raise ValueError(msg)

    def _repr_defn(self):
        r"""
        Return a string representation of the definition of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: X = ChowScheme(1, 'h', 1, 'h^2', name='P1')
            sage: X.base_morphism()._repr_defn()
            'Structure map'
        """
        return 'Structure map'
