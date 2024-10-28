# sage.doctest: needs sage.combinat sage.graphs sage.modules
r"""
Grossman-Larson Hopf Algebras

AUTHORS:

- Frédéric Chapoton (2017)
"""

# ****************************************************************************
#       Copyright (C) 2017 Frédéric Chapoton
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.hopf_algebras import HopfAlgebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.alphabet import Alphabet
from sage.combinat.rooted_tree import (RootedTrees, RootedTree,
                                       LabelledRootedTrees,
                                       LabelledRootedTree)
from sage.misc.cachefunc import cached_method
from sage.categories.rings import Rings
from sage.sets.family import Family
from sage.rings.integer_ring import ZZ
from itertools import combinations, product

# we use a fixed special symbol for the fake root
ROOT = '#'


class GrossmanLarsonAlgebra(CombinatorialFreeModule):
    r"""
    The Grossman-Larson Hopf Algebra.

    The Grossman-Larson Hopf Algebras are Hopf algebras with a basis
    indexed by forests of decorated rooted trees. They are the
    universal enveloping algebras of free pre-Lie algebras, seen
    as Lie algebras.

    The Grossman-Larson Hopf algebra on a given set `E` has an
    explicit description using rooted forests. The underlying vector
    space has a basis indexed by finite rooted forests endowed with a
    map from their vertices to `E` (called the "labeling").
    In this basis, the product of two
    (decorated) rooted forests `S * T` is a sum over all maps from
    the set of roots of `T` to the union of a singleton `\{\#\}` and
    the set of vertices of `S`. Given such a map, one defines a new
    forest as follows. Starting from the disjoint union of all rooted trees
    of `S` and `T`, one adds an edge from every root of `T` to its
    image when this image is not the fake vertex labelled ``#``.
    The coproduct sends a rooted forest `T` to the sum of all tensors
    `T_1 \otimes T_2` obtained by splitting the connected components
    of `T` into two subsets and letting `T_1` be the forest formed
    by the first subset and `T_2` the forest formed by the second.
    This yields a connected graded Hopf algebra (the degree of a
    forest is its number of vertices).

    See [Pana2002]_ (Section 2) and [GroLar1]_.
    (Note that both references use rooted trees rather than rooted
    forests, so think of each rooted forest grafted onto a new root.
    Also, the product is reversed, so they are defining the opposite
    algebra structure.)

    .. WARNING::

        For technical reasons, instead of using forests as labels for
        the basis, we use rooted trees. Their root vertex should be
        considered as a fake vertex. This fake root vertex is labelled
        ``'#'`` when labels are present.

    EXAMPLES::

        sage: G = algebras.GrossmanLarson(QQ, 'xy')
        sage: x, y = G.single_vertex_all()
        sage: ascii_art(x*y)
        B  + B
         #      #_
         |     / /
         x    x y
         |
         y

        sage: ascii_art(x*x*x)
        B  + B     + 3*B     + B
         #      #         #_      _#__
         |      |        / /     / / /
         x      x_      x x     x x x
         |     / /        |
         x    x x         x
         |
         x

    The Grossman-Larson algebra is associative::

        sage: z = x * y
        sage: x * (y * z) == (x * y) * z
        True

    It is not commutative::

        sage: x * y == y * x
        False

    When ``None`` is given as input, unlabelled forests are used instead;
    this corresponds to a `1`-element set `E`::

        sage: G = algebras.GrossmanLarson(QQ, None)
        sage: x = G.single_vertex_all()[0]
        sage: ascii_art(x*x)
        B  + B
         o      o_
         |     / /
         o    o o
         |
         o

    .. NOTE::

        Variables names can be ``None``, a list of strings, a string
        or an integer. When ``None`` is given, unlabelled rooted
        forests are used. When a single string is given, each letter is taken
        as a variable. See
        :func:`sage.combinat.words.alphabet.build_alphabet`.

    .. WARNING::

        Beware that the underlying combinatorial free module is based
        either on ``RootedTrees`` or on ``LabelledRootedTrees``, with no
        restriction on the labellings. This means that all code calling
        the :meth:`basis` method would not give meaningful results, since
        :meth:`basis` returns many "chaff" elements that do not belong to
        the algebra.

    REFERENCES:

    - [Pana2002]_

    - [GroLar1]_
    """
    @staticmethod
    def __classcall_private__(cls, R, names=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: F1 = algebras.GrossmanLarson(QQ, 'xyz')
            sage: F2 = algebras.GrossmanLarson(QQ, ['x','y','z'])
            sage: F3 = algebras.GrossmanLarson(QQ, Alphabet('xyz'))
            sage: F1 is F2 and F1 is F3
            True
        """
        if names is not None:
            if names not in ZZ and ',' in names:
                names = [u for u in names if u != ',']
            names = Alphabet(names)

        if R not in Rings():
            raise TypeError("argument R must be a ring")
        return super().__classcall__(cls, R, names)

    def __init__(self, R, names=None):
        """
        Initialize ``self``.

        TESTS::

            sage: A = algebras.GrossmanLarson(QQ, '@'); A
            Grossman-Larson Hopf algebra on one generator ['@']
            over Rational Field
            sage: TestSuite(A).run() # long time

            sage: F = algebras.GrossmanLarson(QQ, 'xy')
            sage: TestSuite(F).run() # long time

            sage: A = algebras.GrossmanLarson(QQ, None); A
            Grossman-Larson Hopf algebra on one generator ['o'] over
            Rational Field

            sage: F = algebras.GrossmanLarson(QQ, ['x','y']); F
            Grossman-Larson Hopf algebra on 2 generators ['x', 'y']
            over Rational Field

            sage: A = algebras.GrossmanLarson(QQ, []); A
            Grossman-Larson Hopf algebra on 0 generators [] over
            Rational Field
        """
        if names is None:
            Trees = RootedTrees()
            key = RootedTree.sort_key
            self._alphabet = Alphabet(['o'])
        else:
            Trees = LabelledRootedTrees()
            key = LabelledRootedTree.sort_key
            self._alphabet = names
        # Here one would need LabelledRootedTrees(names)
        # so that one can restrict the labels to some fixed set

        cat = HopfAlgebras(R).WithBasis().Graded()
        CombinatorialFreeModule.__init__(self, R, Trees,
                                         latex_prefix='',
                                         sorting_key=key,
                                         category=cat)

    def variable_names(self):
        r"""
        Return the names of the variables.

        This returns the set `E` (as a family).

        EXAMPLES::

            sage: R = algebras.GrossmanLarson(QQ, 'xy')
            sage: R.variable_names()
            {'x', 'y'}

            sage: R = algebras.GrossmanLarson(QQ, ['a','b'])
            sage: R.variable_names()
            {'a', 'b'}

            sage: R = algebras.GrossmanLarson(QQ, 2)
            sage: R.variable_names()
            {0, 1}

            sage: R = algebras.GrossmanLarson(QQ, None)
            sage: R.variable_names()
            {'o'}
        """
        return self._alphabet

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: algebras.GrossmanLarson(QQ, '@')  # indirect doctest
            Grossman-Larson Hopf algebra on one generator ['@'] over Rational Field
            sage: algebras.GrossmanLarson(QQ, None)  # indirect doctest
            Grossman-Larson Hopf algebra on one generator ['o'] over Rational Field
            sage: algebras.GrossmanLarson(QQ, ['a','b'])
            Grossman-Larson Hopf algebra on 2 generators ['a', 'b'] over Rational Field
        """
        n = len(self.single_vertex_all())
        if n == 1:
            gen = "one generator"
        else:
            gen = "{} generators".format(n)
        s = "Grossman-Larson Hopf algebra on {} {} over {}"
        try:
            return s.format(gen, self._alphabet.list(), self.base_ring())
        except NotImplementedError:
            return s.format(gen, self._alphabet, self.base_ring())

    def single_vertex(self, i):
        r"""
        Return the ``i``-th rooted forest with one vertex.

        This is the rooted forest with just one vertex, labelled by the
        ``i``-th element of the label list.

        .. SEEALSO:: :meth:`single_vertex_all`.

        INPUT:

        - ``i`` -- nonnegative integer

        EXAMPLES::

            sage: F = algebras.GrossmanLarson(ZZ, 'xyz')
            sage: F.single_vertex(0)
            B[#[x[]]]

            sage: F.single_vertex(4)
            Traceback (most recent call last):
            ...
            IndexError: argument i (= 4) must be between 0 and 2
        """
        G = self.single_vertex_all()
        n = len(G)
        if i < 0 or not i < n:
            m = "argument i (= {}) must be between 0 and {}".format(i, n - 1)
            raise IndexError(m)
        return G[i]

    def single_vertex_all(self):
        """
        Return the rooted forests with one vertex in ``self``.

        They freely generate the Lie algebra of primitive elements
        as a pre-Lie algebra.

        .. SEEALSO:: :meth:`single_vertex`.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(ZZ, 'fgh')
            sage: A.single_vertex_all()
            (B[#[f[]]], B[#[g[]]], B[#[h[]]])

            sage: A = algebras.GrossmanLarson(QQ, ['x1','x2'])
            sage: A.single_vertex_all()
            (B[#[x1[]]], B[#[x2[]]])

            sage: A = algebras.GrossmanLarson(ZZ, None)
            sage: A.single_vertex_all()
            (B[[[]]],)
        """
        Trees = self.basis().keys()
        return tuple(Family(self._alphabet,
                            lambda a: self.monomial(Trees([Trees([], a)], ROOT))))

    def _first_ngens(self, n):
        """
        Return the first generators.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(QQ, ['x1','x2'])
            sage: A._first_ngens(2)
            (B[#[x1[]]], B[#[x2[]]])

            sage: A = algebras.GrossmanLarson(ZZ, None)
            sage: A._first_ngens(1)
            (B[[[]]],)
        """
        return self.single_vertex_all()[:n]

    def change_ring(self, R):
        """
        Return the Grossman-Larson algebra in the same variables over `R`.

        INPUT:

        - ``R`` -- a ring

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(ZZ, 'fgh')
            sage: A.change_ring(QQ)
            Grossman-Larson Hopf algebra on 3 generators ['f', 'g', 'h']
            over Rational Field
        """
        return GrossmanLarsonAlgebra(R, names=self.variable_names())

    def degree_on_basis(self, t):
        """
        Return the degree of a rooted forest in the Grossman-Larson algebra.

        This is the total number of vertices of the forest.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(QQ, '@')
            sage: RT = A.basis().keys()
            sage: A.degree_on_basis(RT([RT([])]))
            1
        """
        return t.node_number() - 1

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(QQ, 'xy')
            sage: A.an_element()
            B[#[x[]]] + 2*B[#[x[x[]]]] + 2*B[#[x[], x[]]]
        """
        o = self.single_vertex(0)
        return o + 2 * o * o

    def some_elements(self):
        """
        Return some elements of the Grossman-Larson Hopf algebra.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(QQ, None)
            sage: A.some_elements()
            [B[[[]]], B[[]] + B[[[[]]]] + B[[[], []]],
            4*B[[[[]]]] + 4*B[[[], []]]]

        With several generators::

            sage: A = algebras.GrossmanLarson(QQ, 'xy')
            sage: A.some_elements()
            [B[#[x[]]],
             B[#[]] + B[#[x[x[]]]] + B[#[x[], x[]]],
             B[#[x[x[]]]] + 3*B[#[x[y[]]]] + B[#[x[], x[]]] + 3*B[#[x[], y[]]]]
        """
        o = self.single_vertex(0)
        o1 = self.single_vertex_all()[-1]
        x = o * o
        y = o * o1
        return [o, 1 + x, x + 3 * y]

    def product_on_basis(self, x, y):
        """
        Return the product of two forests `x` and `y`.

        This is the sum over all possible ways for the components
        of the forest `y` to either fall side-by-side with components
        of `x` or be grafted on a vertex of `x`.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(QQ, None)
            sage: RT = A.basis().keys()
            sage: x = RT([RT([])])
            sage: A.product_on_basis(x, x)
            B[[[[]]]] + B[[[], []]]

        Check that the product is the correct one::

            sage: A = algebras.GrossmanLarson(QQ, 'uv')
            sage: RT = A.basis().keys()
            sage: Tu = RT([RT([],'u')],'#')
            sage: Tv = RT([RT([],'v')],'#')
            sage: A.product_on_basis(Tu, Tv)
            B[#[u[v[]]]] + B[#[u[], v[]]]
        """
        return self.sum(self.basis()[x.single_graft(y, graftingFunction)]
                        for graftingFunction in
                        product(list(x.paths()), repeat=len(y)))

    def one_basis(self):
        """
        Return the empty rooted forest.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(QQ, 'ab')
            sage: A.one_basis()
            #[]

            sage: A = algebras.GrossmanLarson(QQ, None)
            sage: A.one_basis()
            []
        """
        Trees = self.basis().keys()
        return Trees([], ROOT)

    def coproduct_on_basis(self, x):
        """
        Return the coproduct of a forest.

        EXAMPLES::

            sage: G = algebras.GrossmanLarson(QQ,2)
            sage: x, y = G.single_vertex_all()
            sage: ascii_art(G.coproduct(x))  # indirect doctest
            1 # B  + B  # 1
                 #    #
                 |    |
                 0    0

            sage: Delta_xy = G.coproduct(y*x)
            sage: ascii_art(Delta_xy)  # random indirect doctest
            1 # B     + 1 # B  + B  # B  + B     # 1 + B  # B  + B  # 1
                   #_        #    #    #      #_        #    #    #
                  / /        |    |    |     / /        |    |    |
                 0 1         1    0    1    0 1         1    0    1
                             |                                    |
                             0                                    0

        TESTS::

            sage: Delta_xy.coefficients()
            [1, 1, 1, 1, 1, 1]
            sage: sortkey = G.print_options()['sorting_key']
            sage: doublekey = lambda tt: (sortkey(tt[0]), sortkey(tt[1]))
            sage: sorted(Delta_xy.monomial_coefficients(), key=doublekey)
            [(#[], #[1[0[]]]),
             (#[], #[0[], 1[]]),
             (#[0[]], #[1[]]),
             (#[1[]], #[0[]]),
             (#[1[0[]]], #[]),
             (#[0[], 1[]], #[])]
        """
        B = self.basis()
        Trees = B.keys()
        subtrees = list(x)
        num_subtrees = len(subtrees)
        indx = list(range(num_subtrees))
        return sum(B[Trees([subtrees[i] for i in S], ROOT)].tensor(
                   B[Trees([subtrees[i] for i in indx if i not in S], ROOT)])
                   for k in range(num_subtrees + 1)
                   for S in combinations(indx, k))

    def counit_on_basis(self, x):
        """
        Return the counit on a basis element.

        This is zero unless the forest `x` is empty.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(QQ, 'xy')
            sage: RT = A.basis().keys()
            sage: x = RT([RT([],'x')],'#')
            sage: A.counit_on_basis(x)
            0
            sage: A.counit_on_basis(RT([],'#'))
            1
        """
        if x.node_number() == 1:
            return self.base_ring().one()
        return self.base_ring().zero()

    def antipode_on_basis(self, x):
        """
        Return the antipode of a forest.

        EXAMPLES::

            sage: G = algebras.GrossmanLarson(QQ,2)
            sage: x, y = G.single_vertex_all()
            sage: G.antipode(x)  # indirect doctest
            -B[#[0[]]]

            sage: G.antipode(y*x)  # indirect doctest
            B[#[0[1[]]]] + B[#[0[], 1[]]]
        """
        B = self.basis()
        Trees = B.keys()
        subtrees = list(x)
        if not subtrees:
            return self.one()
        num_subtrees = len(subtrees)
        indx = list(range(num_subtrees))
        return sum(- self.antipode_on_basis(Trees([subtrees[i] for i in S], ROOT))
                   * B[Trees([subtrees[i] for i in indx if i not in S], ROOT)]
                   for k in range(num_subtrees)
                   for S in combinations(indx, k))

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: R = algebras.GrossmanLarson(QQ, 'xy')
            sage: x, y = R.single_vertex_all()
            sage: R(x)
            B[#[x[]]]
            sage: R(x+4*y)
            B[#[x[]]] + 4*B[#[y[]]]

            sage: Trees = R.basis().keys()
            sage: R(Trees([],'#'))
            B[#[]]

            sage: D = algebras.GrossmanLarson(ZZ, 'xy')
            sage: X, Y = D.single_vertex_all()
            sage: R(X-Y).parent()
            Grossman-Larson Hopf algebra on 2 generators ['x', 'y'] over Rational Field

        TESTS::

            sage: Trees = R.basis().keys()
            sage: R(Trees([],'x'))
            Traceback (most recent call last):
            ...
            ValueError: incorrect root label

            sage: R.<x,y> = algebras.GrossmanLarson(QQ)
            sage: R(x) is x
            True
            sage: S.<z> = algebras.GrossmanLarson(GF(3))
            sage: R(z)
            Traceback (most recent call last):
            ...
            TypeError: not able to convert this to this algebra
        """
        if (isinstance(x, (RootedTree, LabelledRootedTree))
                and x in self.basis().keys()):
            if hasattr(x, 'label') and x.label() != ROOT:
                raise ValueError('incorrect root label')
            return self.monomial(x)
        try:
            P = x.parent()
            if isinstance(P, GrossmanLarsonAlgebra):
                if P is self:
                    return x
                if self._coerce_map_from_(P):
                    return self.element_class(self, x.monomial_coefficients())
        except AttributeError:
            raise TypeError('not able to convert this to this algebra')
        else:
            raise TypeError('not able to convert this to this algebra')
        # Ok, not an element (or should not be viewed as one).

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - Grossman-Larson Hopf algebras whose set `E` of labels is
          a subset of the corresponding set of ``self``, and whose base
          ring has a coercion map into ``self.base_ring()``

        EXAMPLES::

            sage: F = algebras.GrossmanLarson(GF(7), 'xyz'); F
            Grossman-Larson Hopf algebra on 3 generators ['x', 'y', 'z']
             over Finite Field of size 7

        Elements of the Grossman-Larson Hopf algebra canonically coerce in::

            sage: x, y, z = F.single_vertex_all()
            sage: F.coerce(x+y) == x+y
            True

        The Grossman-Larson Hopf algebra over `\ZZ` on `x, y, z`
        coerces in, since `\ZZ` coerces to `\GF{7}`::

            sage: G = algebras.GrossmanLarson(ZZ, 'xyz')
            sage: Gx,Gy,Gz = G.single_vertex_all()
            sage: z = F.coerce(Gx+Gy); z
            B[#[x[]]] + B[#[y[]]]
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the Grossman-Larson
        algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

            sage: G.coerce(y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Grossman-Larson Hopf algebra
             on 3 generators ['x', 'y', 'z'] over Finite Field of size
             7 to Grossman-Larson Hopf algebra on 3 generators ['x', 'y', 'z']
             over Integer Ring

        TESTS::

            sage: F = algebras.GrossmanLarson(ZZ, 'xyz')
            sage: G = algebras.GrossmanLarson(QQ, 'xyz')
            sage: H = algebras.GrossmanLarson(ZZ, 'y')
            sage: F._coerce_map_from_(G)
            False
            sage: G._coerce_map_from_(F)
            True
            sage: F._coerce_map_from_(H)
            True
            sage: F._coerce_map_from_(QQ)
            False
            sage: G._coerce_map_from_(QQ)
            False
            sage: F.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
        """
        # Grossman-Larson algebras containing the same variables
        # over any base that coerces in:
        if isinstance(R, GrossmanLarsonAlgebra):
            if all(x in self.variable_names() for x in R.variable_names()):
                if self.base_ring().has_coerce_map_from(R.base_ring()):
                    return True
        return False
