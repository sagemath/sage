# -*- coding: utf-8 -*-
r"""
ChowSchemes

A ChowScheme is thought about as a smooth algebraic variety having a given
Chow ring. A morphism between ChowSchemes is thought about as a morphism
between these varieties represented by the induced morphism on the level of
Chow rings.

As an example, consider `\mathbb{P}^2` and `\mathbb{P}^5` and the morphism
`f:\mathbb{P}^2\longrightarrow\mathbb{P}^5` given by the Veronese embedding
`[x:y:z]\mapsto [x^2:y^2:z^2:yz:xz:xy]` which we encode as follows::

    sage: P2 = ChowScheme(2, 'h', 1, 'h^3')
    sage: P5 = ChowScheme(5, 'k', 1, 'k^6')
    sage: f = P2.hom(['2*h'], P5)

In the first two lines, we define `\mathbb{P}^2` and `\mathbb{P}^5` by their
(rational) Chow rings respectively `A^*(\mathbb{P}^2)=\QQ[h]/(h^3)` and
`A^*(\mathbb{P}^5)=\QQ[k]/(k^6)`. The Veronese embedding `f` induced on
the level of Chow rings

.. math::

    f^*:\ A^*(\mathbb{P}^5)\ \longrightarrow\ A^*(\mathbb{P}^2)

is given by sending the generator `k` to `2*h` which explains the third line.

Actually the projective space is already defined in the libraries and one
might code as follows::

    sage: P2 = Proj(2)
    sage: P5 = Proj(5, hyperplane_class='k')
    sage: f = P2.hom(['2*h'], P5)

TESTS::

    sage: X = ChowScheme(2, ['H', 'E'], [1, 1], ['E*H', 'E^2+H^2'], 'H^2')
    sage: TestSuite(X).run(skip=["_test_an_element", "_test_some_elements"])

AUTHORS:

- Manfred Lehn (2013)
- Christoph Sorger (2013)
"""
# ****************************************************************************
#       Copyright (C) 2013 Manfred Lehn <lehn@mathematik.uni-mainz.de>
#       Copyright (C) 2023 Christoph Sorger <christoph.sorger@univ-nantes.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.all import Rings
_Rings = Rings() # type: ignore
from sage.categories.map import Map
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.schemes.chow.ring import ChowRing, PointChowRing, is_chowRing
from sage.schemes.chow.schemes import ChowSchemes
from sage.schemes.chow.sheaf import Sheaf
from sage.structure.parent import Parent


def is_chowScheme(x):
    """
    Test whether ``x`` is a chow scheme.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean. Whether ``x`` derives from :class:`ChowScheme_generic`.

    TESTS::

        sage: X = ChowScheme(1)
        sage: is_chowScheme(X)
        True
        sage: is_chowScheme(1)
        False
        sage: is_chowScheme(PointChowScheme)
        True
    """
    return isinstance(x, ChowScheme_generic)


def ChowScheme(dimension,
               generators=None, degrees=None, relations=None,
               point_class=None,
               names=None, name=None, latex_name=None):
    r"""
    Returns a ChowScheme, given its dimension and (optional) generators,
    degrees, relations, point_class, names, name, latex_name.

    INPUT:

    - ``dimension`` -- The dimension of the ChowScheme
    - ``generators`` -- An optional (list of) strings, the generators
    - ``degrees`` -- An optional (list of) integers, the degrees of the generators
    - ``relations`` -- An optional (list of) strings, the relations between the generators
    - ``point_class`` -- An optional point_class of the ChowScheme.
    - ``names`` -- An optional (list of) strings, names of the generators in the ChowScheme
    - ``name`` -- An optional string, the name of the ChowScheme
    - ``latex_name``-- An optional string, the latex representation of the ChowScheme

    OUTPUT:

    - A ChowScheme

    EXAMPLES::

        sage: P2.<h> = ChowScheme(2, 'h', 1, 'h^3')
        sage: h^3
        0

    Note that the optional name is not the unique string representation of the
    ChowScheme (as given by python's `:meth:repr()`) but is rather used in the
    string representations of other objects as morphisms or bundles on the
    ChowScheme in order to simplify the output.

    EXAMPLES::

        sage: P2.<h> = ChowScheme(2, 'h', 1, 'h^3', 'h^2', name='P2'); P2
        ChowScheme(2, 'h', 1, 'h^3', 'h^2')
        sage: str(P2)
        'P2'
        sage: P2.<h> = Proj(2, name='P2'); P2
        ChowScheme(2, 'h', 1, 'h^3', 'h^2')
        sage: P2.o(-1)  # The Hopf bundle on P2
        Bundle(P2, 1, [1, -h])


    TESTS::
        sage: P2.objgen()
        (ChowScheme(2, 'h', 1, 'h^3', 'h^2'), h)
        sage: P2.objgens()
        (ChowScheme(2, 'h', 1, 'h^3', 'h^2'), (h,))
        sage: P2.inject_variables()
        Defining h
        sage: P2.latex_variable_names()
        ['h']
        sage: P2.variable_names()
        ('h',)
    """
    # Build the ChowRing
    a_name = 'A(%s)' % name if name else None
    a_latex_name = 'A^{*}(%s)' % latex_name if latex_name else None
    A = ChowRing(generators, degrees, relations, names, a_name, a_latex_name)
    # Set dimension and point_class
    A.set_dimension(dimension)
    if point_class is not None:
        A.set_point_class(point_class)
    return ChowScheme_generic(A, name, latex_name)


class ChowScheme_generic(Parent):
    """
    The base class for all ChowSchemes.
    """

    def __init__(self, R=None, name=None, latex_name=None):
        """
        Construct a :class:`ChowScheme <ChowScheme_generic>`.

        INPUT:

        - ``R`` -- a ChowRing or a ring morphism between ChowRings

        - ``name`` -- an optional name for the ChowScheme

        - ``latex_name`` -- an optional latex name for the ChowScheme

        OUTPUT:

        - :class:`ChowScheme <ChowScheme_generic>`.

        TESTS::

            sage: from sage.schemes.chow.scheme import ChowScheme_generic
            sage: X = ChowScheme_generic(); X
            ChowScheme(0, '1')
        """

        # Options
        self._name = name
        self._latex_name = latex_name
        # Private variables
        self._sheaves = dict()
        # Ring or ring morphism check
        from sage.schemes.chow.morphism import is_chowSchemeMorphism
        if R is None:
            self._chowring = PointChowRing
            self._base_chowring = PointChowRing
            self._base_chowring_morphism = PointChowRing.identity_morphism()
        elif is_chowRing(R):
            self._chowring = R
            self._base_chowring = PointChowRing
            self._base_chowring_morphism = self._base_chowring.hom([], R)
        elif isinstance(R, Map) and R.category_for().is_subcategory(_Rings):
            self._chowring = R.codomain()
            self._base_chowring = R.domain()
            self._base_chowring_morphism = R
        elif is_chowSchemeMorphism(R):
            self._chowring = R.domain().chowring()
            self._base_chowring = R.codomain().chowring()
            self._base_chowring_morphism = R.chowring_morphism()
            self._base_chowscheme = R.codomain()
            # Keep sheaves
            for k in R.domain()._sheaves:
                E = Sheaf(self, ch=R.domain().sheaves[k].chern_character())
                self._sheaves[k] = E
        else:
            raise TypeError("Expect ChowRing or a morphism between ChowRings")
        category = ChowSchemes(self.base_chowscheme())
        Parent.__init__(self, self.base_chowscheme(), category=category)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - a string.

        TESTS::

            sage: P1 = ChowScheme(1, 'h', 1, 'h^2', 'h'); P1._repr_()
            "ChowScheme(1, 'h', 1, 'h^2', 'h')"
        """
        if self is PointChowScheme:
            return "PointChowScheme"

        result = "ChowScheme(%d" % self.dimension()

        if self.ngens() == 0:
            pass
        elif self.ngens() == 1:
            result += ", '%s'" % self.variable_name()
            result += ", %s" % self.deg()
        else:
            result += ", %s" % list(self.variable_names())
            result += ", %s" % list(self.degs())
        if self.nrels() == 0:
            pass
        elif self.nrels() == 1:
            result += ", '%s'" % self.rel()
        else:
            result += ", %s" % [str(result) for result in self.rels()]

        pc = self.point_class()
        result += ", '%s')" % pc if pc is not None else ")"

        return result

    def __str__(self):
        r"""
        Return a text representation of this ChowScheme. This is generally
        the name specified during definition or the string representation
        if a name has not been specified.

        EXAMPLES::

            sage: X = ChowScheme(1, 'h', 1, 'h^2', name='X')
            sage: X
            ChowScheme(1, 'h', 1, 'h^2')
            sage: str(X)
            'X'


        """
        return self._name if self._name else repr(self)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        OUTPUT:

        - a string.

        TESTS::

            sage: P1 = ChowScheme(1, 'h', 1, 'h^2', latex_name='\\mathbb{P}^1')
            sage: P1._latex_()
            '\\mathbb{P}^1'

        """
        if self._latex_name:
            return self._latex_name
        return self._repr_()

    def _cache_key(self):
        return(self.parent(),str(self))

    def __hash__(self):
        return hash(type(self))

    def __eq__(self, other):
        """
        Return True if the two ChowScheme are "equal", e.g. if their underlying
        ChowRings coincide, else False.

        TESTS::

            sage: P1 = ChowScheme(1, 'h', 1, 'h^2')
            sage: S1 = ChowScheme(1, 'k', 1, 'k^2')
            sage: P1 == S1
            True
            sage: S1 = ChowScheme(1, 'k', 1, 'k^2', 'k')
            sage: P1 == S1
            False
            sage: S1 = ChowScheme(2, 'k', 2, 'k^2')
            sage: P1 == S1
            False

            sage: P4 = ChowScheme(4, 'w', 1, 'w^5', name='P4')  # P4
            sage: X = ChowScheme(6, ['e1', 'e2', 'w'], [1, 2, 1], ['w^5'])
            sage: Y = X.base_change(X.hom('w', P4))
            sage: X.chowring() == Y.chowring()
            True
            sage: X == Y
            False
            sage: Z = X.base_change(X.hom([0], P4))
            sage: Z.chowring() == Y.chowring()
            True
            sage: Z.base_chowring() == Y.base_chowring()
            True
            sage: Z == Y
            False
        """
        # Check instances
        if not isinstance(self, ChowScheme_generic):
            return False
        if not isinstance(other, ChowScheme_generic):
            return False
        # Check Chow rings
        if self.chowring() != other.chowring():
            return False
        # Check base Chow rings
        if self.base_chowring() != other.base_chowring():
            return False
        # Check base Chow ring morphism
        if self.base_chowring_morphism() != other.base_chowring_morphism():
            return False
        # All checks done, return True
        return True

    def __ne__(self, other):
        """
        Return True if the two ChowScheme are "not equal", else False.

        TESTS::

            sage: P1 = ChowScheme(1, 'h', 1, 'h^2')
            sage: S1 = ChowScheme(1, 'k', 1, 'k^2')
            sage: P1 != S1
            False
        """
        return not self.__eq__(other)

    def __call__(self, *args):
        r"""
        EXAMPLES::

            sage: P1 = ChowScheme(1, 'h', 1, 'h^2', name='P1')
            sage: P1(PointChowScheme)
            Set of morphisms
              From: PointChowScheme
              To:   P1
        """
        if len(args) == 0:
            raise TypeError('You need to specify at least one argument.')
        S = args[0]
        if is_chowScheme(S):
            return S.Hom(self)

    def __mul__(self, other):
        r"""
        Return the product of two ChowSchemes with affine stratification.

        .. Warning::

            Recall that there is no Künneth formula for the Chow ring of a
            product of varieties. Already for the product of two smooth curves
            with genus `geq 1`, there is no algorithm for calculating `A^1`,
            and we do not know what A^2 is, except that it can't be
            finite-dimensional. There is however, a Künneth formula for
            varieties with affine stratification as projective spaces
            or more generally Grassmannians. This what this method returns.

        EXAMPLES:

            sage: P2xP2 = Proj(2, 'h') * Proj(2, 'k')
            sage: P8 = Proj(8, 'l')
            sage: f = P2xP2.hom(['h+k'], P8)  # Segre map P2xP2 -> P8
            sage: g = Blowup(f)
            sage: B = g.codomain()
            sage: B.betti_numbers()
            [1, 2, 4, 7, 8, 7, 4, 2, 1]
        """
        from .bundle import Bundle
        dimension = self.dimension() + other.dimension()
        generators = self.variable_names() + other.variable_names()
        degrees = self.degs() + other.degs()
        relations = [str(r) for r in self.rels()] + [str(r) for r in other.rels()]
        point_class = '%s*%s' % (self.point_class(), other.point_class())
        # TODO names, name, latex_name
        prod = ChowScheme(dimension, generators, degrees, relations,
                          point_class)
        TS = self.tangent_bundle()
        TO = other.tangent_bundle()
        T1 = Bundle(prod, TS.rank(), [str(c) for c in TS.chern_classes()])
        T2 = Bundle(prod, TO.rank(), [str(c) for c in TO.chern_classes()])
        prod.sheaves['tangent'] = T1 + T2
        return prod

    @property
    def sheaves(self):
        r"""
        Return a dictionary with sheaves on this ChowScheme. The naming
        convention is to suppress the word "bundle" in the keys.
        For instance, the *tangent bundle* has the key 'tangent', the
        *universal sub bundle* the key 'universal_sub' and the
        *universal quotient bundle* the key 'universal_quotient'.


        EXAMPLES::

            sage: P4 = Proj(4, name='P4')  # P4 in the sens of Grothendieck
            sage: P4.sheaves["tangent"]
            Bundle(P4, 4, [1, 5*h, 10*h^2, 10*h^3, 5*h^4])
            sage: P4.sheaves["universal_sub"]
            Bundle(P4, 4, [1, -h, h^2, -h^3, h^4])
            sage: P4.sheaves["universal_quotient"]
            Bundle(P4, 1, [1, h])

        """
        return self._sheaves

    @sheaves.setter
    def sheaves(self, value):
        r"""
        Dictionary. Allows to add sheaves to this ChowScheme.

        EXAMPLES::

            sage: P1.<h> = Proj(1, 'h', name='P1')
            sage: P3.<k> = Proj(3, 'k', name='P3')
            sage: f = P1.hom([3*h], P3)  # Twisted cubic
            sage: NIX = f.upperstar(P3.tangent_bundle()) - P1.tangent_bundle()
            sage: P1.sheaves["normal"] = NIX
            sage: P1.sheaves["normal"]
            Sheaf(P1, 2, [1, 10*h])
        """
        self._sheaves = value

    def _morphism(self, *args, **kwds):
        r"""
        Internal. Construct a morphism determined the arguments.

        TESTS::

            sage: P1.<h> = ChowScheme(1, 'h', 1, 'h^2', 'h', name='P1')
            sage: P3.<k> = ChowScheme(3, 'k', 1, 'k^4', 'k^3', name='P3')
            sage: h = P1._morphism(P1.Hom(P3), [3*h]); h
            ChowScheme morphism:
              From: P1
              To:   P3
              Defn: k |--> 3*h
        """
        from sage.schemes.chow.morphism import ChowSchemeMorphism
        return ChowSchemeMorphism(*args, **kwds)

    def identity_morphism(self):
        """
        Return the identity morphism.

        OUTPUT:

        The identity morphism of the ChowScheme ``self``.

        EXAMPLES::

            sage: X = ChowScheme(1, 'h', 1, 'h^2', name='P1')
            sage: X.identity_morphism()
            ChowScheme endomorphism of P1
              Defn: Identity map
        """
        from sage.schemes.chow.morphism import ChowSchemeMorphism_id
        return ChowSchemeMorphism_id(self)

    def hom(self, x, Y=None, check=True):
        """
        Return the ChowScheme morphism from ``self`` to ``Y`` defined by ``x``.

        INPUT:

        - ``x`` -- anything hat determines a ChowScheme morphism, typically the
          images of the generators of the ChowRing of Y in the ChowRing of X.
          If ``x`` is a scheme, try to determine a natural map to ``x``.

        - ``Y`` -- the codomain scheme (optional). If ``Y`` is not
          given, try to determine ``Y`` from context.

        - ``check`` -- boolean (optional, default=``True``). Whether
          to check the defining data for consistency.

        OUTPUT:

        The ChowScheme morphism from ``self`` to ``Y`` defined by ``x``.

        EXAMPLES::

            sage: P1.<h> = ChowScheme(1, 'h', 1, 'h^2', 'h')
            sage: P3.<k> = ChowScheme(3, 'k', 1, 'k^4', 'k^3')
            sage: f = P1.hom([3*h], P3)  # Twisted cubic
            sage: f(k)
            3*h

        TESTS::

            sage: P = PointChowScheme
            sage: P2.<h> = Proj(2)
            sage: f = P.hom([0], P2)  # PointChowScheme in P2
            sage: f(h)
            0
        """
        if Y is None:
            if is_chowScheme(x):
                return self.Hom(x).natural_map()
            else:
                raise TypeError("unable to determine codomain")
        return self.Hom(Y)(x, check)

    def _Hom_(self, Y, category=None, check=True):
        r"""
        Return the set of ChowScheme morphisms from ``self`` to ``Y``.

        INPUT:

        - ``Y`` -- a ChowScheme. The codomain of the SHom-set.

        - ``category`` -- a category (optional). The category of the
          SHom-set.

        - ``check`` -- boolean (optional, default=``True``). Whether
          to check the defining data for consistency.

        OUTPUT:

        The set of morphisms from ``self`` to ``Y``.

        EXAMPLES::

            sage: P1.<h> = ChowScheme(1, 'h', 1, 'h^2', 'h', name='P1')
            sage: P3.<k> = ChowScheme(3, 'k', 1, 'k^4', 'k^3', name='P3')
            sage: P1._Hom_(P3)
            Set of morphisms
              From: P1
              To:   P3
        """
        from sage.schemes.chow.scheme_homset import ChowSchemeHomset
        return ChowSchemeHomset(self, Y, category=category, check=check)

    def chowring(self):
        """
        Return the Chow Ring of this ChowScheme.

        EXAMPLES::

            sage: P1 = ChowScheme(1, 'h', 1, 'h^2', 'h')
            sage: P1.chowring()
            Quotient of Multivariate Polynomial Ring in h over Rational Field by the ideal (h^2)

        """
        return self._chowring

    def base_chowring(self):
        """
        Return the Chow Ring of the Base Chow Scheme of this ChowScheme.

        EXAMPLES::

            sage: P1 = ChowScheme(1, 'h', 1, 'h^2', 'h')
            sage: P1.base_chowring()
            Quotient of Multivariate Polynomial Ring in no variables over Rational Field by the ideal (0)
        """
        return self._base_chowring

    def base_chowring_morphism(self):
        """
        Return the Chow Ring morphism from the Base Chow Ring to the Chow Ring

        EXAMPLES::

            sage: P1 = ChowScheme(1, 'h', 1, 'h^2', 'h')
            sage: f = P1.base_chowring_morphism()
        """
        return self._base_chowring_morphism

    def base_ring(self):
        r"""
        Return the Chow ring of the base of this ChowScheme.
        Synonym of `:meth:base_chowring()`.

        EXAMPLES::

            sage: P1 = ChowScheme(1, 'h', 1, 'h^2', 'h')
            sage: P1.base_ring()
            Quotient of Multivariate Polynomial Ring in no variables over Rational Field by the ideal (0)
        """
        return self.base_chowring()

    def dimension(self):
        """
        Return the dimension of this ChowScheme.

        OUTPUT:

        An integer, the dimension of this ChowScheme

        EXAMPLES::

            sage: X = ChowScheme(2, ['c1', 'c2'], [1, 2])
            sage: X.dimension()
            2

        TESTS::

            sage: X = PointChowScheme
            sage: X.dimension()
            0
        """
        return self.chowring().dimension()

    def gen(self, i=0):
        r"""
        Return the i-th generator for the Chow ring of this ChowScheme.

        INPUT:

        - ``i`` -- an (optional) integer

        OUTPUT:

        The element of the Chow ring of this ChowScheme corresponding
        to the i-th generator.

        EXAMPLES:

            sage: X = ChowScheme(2, ['c1', 'c2'], [1, 2])
            sage: X.gen(1)
            c2
        """
        return self.chowring().gen(i)

    def gens(self):
        r"""
        Return the generators for the Chow ring of this ChowScheme.

        EXAMPLES:

            sage: X = ChowScheme(2, ['c1', 'c2'], [1, 2])
            sage: X.gens()
            (c1, c2)
        """
        return self.chowring().gens()

    def ngens(self):
        r"""
        Return the number of generators for the Chow ring of this ChowScheme.

        EXAMPLES:

            sage: X = ChowScheme(2, ['c1', 'c2'], [1, 2])
            sage: X.ngens()
            2
        """
        return self.chowring().ngens()

    def gens_dict(self):
        r"""
        Return the dictionary variable_name -- generator of the Chow ring
        of this ChowScheme.

        EXAMPLES::

            sage: X = ChowScheme(2, ['c1', 'c2'], [1, 2])
            sage: gens = X.gens_dict()
            sage: [(key, gens[key]) for key in sorted(gens)]
            [('c1', c1), ('c2', c2)]

        """
        return self.chowring().gens_dict()

    def gens_dict_recursive(self):
        r"""
        Return recursively the dictionary variable_name -- generator of the Chow
        ring of this ChowScheme.

        EXAMPLES::

            sage: X = ChowScheme(2, ['c1', 'c2'], [1, 2])
            sage: gens = X.gens_dict_recursive()
            sage: [(key, gens[key]) for key in sorted(gens)]
            [('c1', c1), ('c2', c2)]
        """
        return self.chowring().gens_dict_recursive()

    def deg(self, i=0):
        """
        Return the degree of the i-th generator of the Chow ring of this
        ChowScheme.  If i is unspecified the degree of the first generator
        (eg i=0) is returned if exists.

        INPUT:

        - ``i`` -- an (optional) integer.

        OUTPUT:

        The degree of the i-th generator.

        EXAMPLES::

            sage: X = ChowScheme(1, 'h', 1, 'h^2', 'h')  # ChowScheme of P1
            sage: X.deg()
            1
            sage: X = ChowScheme(2, ['c1', 'c2'], [1, 2])
            sage: X.deg()
            1
            sage: X.deg(1)
            2

        CORNER CASE::

            sage: X = PointChowScheme
            sage: X.deg()
            Traceback (most recent call last):
            ...
            IndexError: tuple index out of range
        """
        return self.chowring().deg(i)

    def degs(self):
        """
        Return the degrees of the generators of the Chow ring of this
        ChowScheme.

        OUTPUT:

        A tuple of integers representing the degrees of the generators.

        EXAMPLES::

            sage: X = ChowScheme(1, 'h', 1, 'h^2', 'h')  # ChowScheme of P1
            sage: X.degs()
            (1,)
            sage: X = ChowScheme(2, ['c1', 'c2'], [1, 2])
            sage: X.degs()
            (1, 2)

        CORNER CASE::

            sage: X = PointChowScheme
            sage: X.degs()
            ()
        """
        return self.chowring().degs()

    def rel(self, i=0):
        """
        Return the i-th relation of the Chow ring of this ChowScheme.
        If i is unspecified the first relation (eg i=0) is returned if exists.

        INPUT:

        - ``i`` -- an (optional) integer.

        OUTPUT:

        The i-th relation.

        EXAMPLES::

            sage: X = ChowScheme(1, 'h', 1, 'h^2', 'h')  # P1
            sage: X.rel()
            h^2

        TESTS::

            sage: X = ChowScheme(2, ['c1', 'c2'], [1, 2])
            sage: X.rel()
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self.chowring().rel(i)

    def nrels(self):
        """
        Return the number of relations of the Chow ring of this ChowScheme.

        OUTPUT:

        An integer, the number of relations.

        EXAMPLES::

            sage: X = ChowScheme(1, 'h', 1, 'h^2', 'h')
            sage: X.nrels()
            1
            sage: X = ChowScheme(2, ['c1', 'c2'], [1, 2])
            sage: X.nrels()
            0
        """
        return self.chowring().nrels()

    def rels(self):
        """
        Return the relations of Chow ring of this ChowScheme.

        OUTPUT:

        A list of ring elements representing the relations.

        EXAMPLES::

            sage: X = ChowScheme(1, 'h', 1, 'h^2', 'h')
            sage: X.rels()
            [h^2]
            sage: X = ChowScheme(2, ['c1', 'c2'], [1, 2])
            sage: X.rels()
            []

        TESTS::

            sage: X = PointChowScheme
            sage: X.rels()
            []
        """
        return self.chowring().rels()

    def variable_name(self):
        r"""
        Return the (first) variable name of the Chow ring of this ChowScheme.

        EXAMPLES::

            sage: X = ChowScheme(2, ['c1', 'c2'], [1, 2])
            sage: X.variable_name()
            'c1'

        """
        return self.chowring().variable_name()

    def variable_names(self):
        r"""
        Return the variable names of the Chow ring of this ChowScheme.

        EXAMPLES::

            sage: X = ChowScheme(2, ['c1', 'c2'], [1, 2])
            sage: X.variable_names()
            ('c1', 'c2')
        """
        return self.chowring().variable_names()

    def point_class(self):
        """
        Return the point_class of this ChowScheme

        OUTPUT:

        The point_class of this ChowScheme.

        EXAMPLES::

            sage: X = ChowScheme(3, 'h', 1, 'h^4', 'h^3')  # ChowScheme of P3
            sage: X.point_class()
            h^3

            sage: X = ChowScheme(1, 'h', 1, 'h^2')  # P1, no point_class
            sage: X.point_class()

        TESTS::

            sage: X = PointChowScheme
            sage: X.point_class()
            1
        """
        return self.chowring().point_class()

    def set_point_class(self, value):
        """
        Set the point_class of this ChowScheme.

        INPUT:

        - ``v`` -- an element (or its string representation) of the point_class

        EXAMPLES::

            sage: X = ChowScheme(3, 'h', 1, 'h^4')  # ChowScheme of P3
            sage: h = X.gen()
            sage: X.set_point_class(h^3)

            sage: X = ChowScheme(3, 'h', 1, 'h^4')  # ChowScheme of P3
            sage: X.set_point_class('h^3')

        TESTS::

            sage: X = ChowScheme(3, 'h', 1, 'h^4')  # ChowRing of P3
            sage: X.set_point_class('k^3')
            Traceback (most recent call last):
            ...
            TypeError: Can't coerce pointclass to ring element.
        """
        self.chowring().set_point_class(value)

    def is_point(self):
        r"""
        Return True if this ChowScheme corresponds to a point.

        EXAMPLES::

            sage: P2 = Proj(2)
            sage: P2.base().is_point()
            True

        TESTS::

            sage: X = ChowScheme(0, point_class='1')
            sage: X.is_point()
            True
        """
        return self.ngens() == 0 and self.dimension() == 0

    def is_point_chowscheme(self):
        r"""
        Return True if this ChowScheme is the instance PointChowScheme

        EXAMPLES::

            sage: P2 = Proj(2)
            sage: P2.is_point_chowscheme()
            False
            sage: P2.base().is_point_chowscheme()
            True

        TESTS::

            sage: X = ChowScheme(0, point_class='1')
            sage: X.is_point_chowscheme()
            False
        """
        return self is PointChowScheme

    def base_change(self, f):
        r"""
        Return this ChowScheme with base changes to f.codomain().

        INPUT:

        - ``f`` -- a morphism from this ChowScheme to another ChowScheme

        OUTPUT:

        This ChowScheme based changed to f.codomain()

        EXAMPLES::

            sage: P4 = ChowScheme(4, 'w', 1, 'w^5')  # P4
            sage: X = ChowScheme(6, ['e1', 'e2', 'w'], [1, 2, 1], ['w^5'])
            sage: X.base()
            PointChowScheme
            sage: X = X.base_change(X.hom('w', P4))
            sage: X.base()
            ChowScheme(4, 'w', 1, 'w^5')
        """

        from sage.schemes.chow.morphism import is_chowSchemeMorphism
        if is_chowSchemeMorphism(f):
            if self is f.domain():
                return ChowScheme_generic(f, name=self._name,
                                          latex_name=self._latex_name)
            else:
                return TypeError("Domain of ChowScheme morphism wrong.")
        raise TypeError("Expect morphism of Chow schemes, got", f)

    def base_chowscheme(self):
        """
        Return the base ChowScheme of this ChowScheme.

        OUTPUT:

        A chowscheme.

        EXAMPLES::

            sage: P2 = Proj(2)
            sage: P2.base_chowscheme()
            PointChowScheme

        TESTS::

            sage: P = PointChowScheme
            sage: P.base_chowscheme()
            PointChowScheme
        """

        try:
            return self._base_chowscheme
        except AttributeError:

            # During initialization of the PointChowScheme we return
            # None in order to avoid infinite recursion. This implies
            # that the category of the PointChowScheme is simply the
            # Category of ChowSchemes
            # Then, while defining the PointChowScheme outside the class we
            # explicitly set the _base_chowscheme attribute to itself.
            # There is probably a more elegant solution to this.
            # Feel free to patch.

            if self.chowring().is_point_chowring():
                return None

            if self.base_chowring().is_point_chowring():
                return PointChowScheme

            self._base_chowscheme = ChowScheme_generic(self.base_chowring())
            return self._base_chowscheme

    def base(self):
        r"""
        Return the base ChowScheme of this ChowScheme. Synonym for
         `meth:base_chowscheme()`

        OUTPUT :

        A ChowScheme.

        Examples::

            sage: G = Grass(4, 2)
            sage: G.base()
            PointChowScheme
            sage: Q = G.sheaves["universal_quotient"]
            sage: ProjQ = ProjBundle(Q.symm(2).dual())  # Proj(Q)->G
            sage: ProjQ.base()
            ChowScheme(4, ['c1', 'c2'], [1, 2], ['c2^3', 'c1*c2^2', 'c1^2*c2 - c2^2', 'c1^3 - 2*c1*c2'], 'c2^2')
        """

        return self.base_chowscheme()

    def base_morphism(self):
        """
        Return the structure morphism from ``self`` to its base chow scheme.

        OUTPUT:

        A ChowScheme morphism.

        EXAMPLES::

            sage: G = Grass(4, 3)
            sage: Q = G.sheaves["universal_quotient"]
            sage: PG = ProjBundle(Q.symm(2).dual())
            sage: PG.base_morphism()
            ChowScheme morphism:
              From: Proj(Bundle(Grass(4, 3), 6, [1, -4*c, 10*c^2, -20*c^3]))
              To:   Grass(4, 3)
              Defn: Structure map
        """
        try:
            return self._base_morphism
        except AttributeError:
            from sage.schemes.chow.morphism import \
                ChowSchemeMorphism_structure_map
            bcs = self.base_chowscheme()
            bcm = self.base_chowring_morphism()
            sm = ChowSchemeMorphism_structure_map(self.Hom(bcs), bcm)
            self._base_morphism = sm
            return self._base_morphism

    def relative_dimension(self):
        r"""
        Return the relative dimension of this ChowScheme.

        EXAMPLES::

            sage: G = Grass(4, 3)
            sage: Q = G.sheaves["universal_quotient"]
            sage: Q.symm(2).dual().rank()
            6
            sage: ProjG = ProjBundle(Q.symm(2).dual())
            sage: ProjG.relative_dimension()
            5
        """
        self_dim, base_dim = self.chowring().dimension(), self.base_chowring().dimension()
        if not (isinstance(self_dim, int) and isinstance(base_dim, int)):
            raise ValueError("The dimensions of the Chow scheme and its base have to be defined before computing its relative dimension.")
        return self_dim - base_dim

    @cached_method
    def betti_numbers(self):
        """
        Return the Betti number of this ChowScheme.

        EXAMPLES:

            sage: P2 = Proj(2)
            sage: P = PointChowScheme
            sage: f = P.hom([0], P2)
            sage: B = Blowup(f).codomain()  # P2 blown up in a point
            sage: B.betti_numbers()
            [1, 2, 1]

        TESTS::

            sage: X = PointChowScheme
            sage: X.betti_numbers()
            [1]
        """
        R = self.chowring()
        dim, rbbd = R.dimension(), R.basis_by_degree()
        if not isinstance(dim, int):
            raise ValueError("The dimension of the Chow scheme has to be defined before computing its Betti numbers.")
        return [len(rbbd[i]) for i in range(dim + 1)]

    @cached_method
    def euler_number(self):
        """
        Return the Euler number of this ChowScheme.

        EXAMPLES:

            sage: G = Proj(2)
            sage: G.euler_number()
            3

        TESTS::

            sage: X = PointChowScheme
            sage: X.euler_number()
            1
        """
        return len(self.chowring().basis())

    def o(self, n=0):
        r"""
        Return `\mathcal{O}(n)` for this ChowScheme.

        INPUT:

        - ``n`` -- an integer

        EXAMPLES::

            sage: G = Grass(4,3)
            sage: G.o(1)
            Bundle(Grass(4, 3), 1, [1, c])
            sage: G.sheaves["universal_quotient"].determinant()
            Bundle(Grass(4, 3), 1, [1, c])

        TESTS::

            sage: X = ChowScheme(1, 'h', 1, 'h^2', 'h')  # P1
            sage: X.o(-2)
            Traceback (most recent call last):
            ...
            RuntimeError: o(1) is undefined yet.

        """
        n = int(n) if isinstance(n, Integer) else n
        if not isinstance(n, int):
            raise ValueError('Only integers allowed here.')
        if n == 0:
            from .bundle import Bundle
            return Bundle(self, 1, [1])

        if 'o1' not in self.sheaves:
            raise RuntimeError('o(1) is undefined yet.')
        return self.sheaves['o1'] ** n

    def tangent_bundle(self):
        r"""
        Return the tangent bundle of this ChowScheme if defined.

        EXAMPLES::

            sage: G = Grass(4,2)
            sage: G.tangent_bundle()
            Bundle(Grass(4, 2), 4, [1, 4*c1, 7*c1^2, 12*c1*c2, 6*c2^2])

        TESTS::

            sage: P = PointChowScheme
            sage: P.tangent_bundle()
            Bundle(PointChowScheme, 0, [1])
        """
        if 'tangent' not in self.sheaves:
            raise RuntimeError('Tangent bundle is undefined yet.')
        return self.sheaves['tangent']

    @cached_method
    def canonical_class(self):
        r"""
        Return the canonical class of this ChowScheme.

        EXAMPLES::

            sage: P2 = Proj(2)
            sage: P2.canonical_class()
            -3*h
        """
        return self.tangent_bundle().dual().determinant().chern_classes()[1]

    @cached_method
    def todd_class(self):
        r"""
        Return the todd class of this ChowScheme.

        EXAMPLES::

            sage: Proj(3).todd_class()
            h^3 + 11/6*h^2 + 2*h + 1

        TESTS::

            sage: P2 = ChowScheme(2, 'h', 1, 'h^3', 'h^2')
            sage: P2.todd_class()
            Traceback (most recent call last):
            ...
            RuntimeError: Tangent bundle is undefined yet.

        """
        return self.tangent_bundle().todd_class()


PointChowScheme_generic = ChowScheme_generic(name='PointChowScheme', latex_name='PointChowScheme')
PointChowScheme = ChowScheme_generic(name='PointChowScheme', latex_name='PointChowScheme')
PointChowScheme._base_chowscheme = PointChowScheme
