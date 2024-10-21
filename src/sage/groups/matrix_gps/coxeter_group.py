# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.graphs
"""
Coxeter Groups As Matrix Groups

This implements a general Coxeter group as a matrix group by using the
reflection representation.

AUTHORS:

- Travis Scrimshaw (2013-08-28): Initial version
"""

##############################################################################
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

import sage.rings.abc

from sage.categories.coxeter_groups import CoxeterGroups
from sage.combinat.root_system.coxeter_matrix import CoxeterMatrix
from sage.groups.matrix_gps.finitely_generated import FinitelyGeneratedMatrixGroup_generic
from sage.groups.matrix_gps.group_element import MatrixGroupElement_generic
from sage.matrix.args import SparseEntry
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.sets.family import Family
from sage.structure.unique_representation import UniqueRepresentation


class CoxeterMatrixGroup(UniqueRepresentation, FinitelyGeneratedMatrixGroup_generic):
    r"""
    A Coxeter group represented as a matrix group.

    Let `(W, S)` be a Coxeter system. We construct a vector space `V`
    over `\RR` with a basis of `\{ \alpha_s \}_{s \in S}` and inner product

    .. MATH::

        B(\alpha_s, \alpha_t) = -\cos\left( \frac{\pi}{m_{st}} \right)

    where we have `B(\alpha_s, \alpha_t) = -1` if `m_{st} = \infty`. Next we
    define a representation `\sigma_s : V \to V` by

    .. MATH::

        \sigma_s \lambda = \lambda - 2 B(\alpha_s, \lambda) \alpha_s.

    This representation is faithful so we can represent the Coxeter group `W`
    by the set of matrices `\sigma_s` acting on `V`.

    INPUT:

    - ``data`` -- a Coxeter matrix or graph or a Cartan type
    - ``base_ring`` -- (default: the universal cyclotomic field or
      a number field) the base ring which contains all values
      `\cos(\pi/m_{ij})` where `(m_{ij})_{ij}` is the Coxeter matrix
    - ``index_set`` -- (optional) an indexing set for the generators

    For finite Coxeter groups, the default base ring is taken to be `\QQ` or
    a quadratic number field when possible.

    For more on creating Coxeter groups, see
    :meth:`~sage.combinat.root_system.coxeter_group.CoxeterGroup`.

    .. TODO::

        Currently the label `\infty` is implemented as `-1` in the Coxeter
        matrix.

    EXAMPLES:

    We can create Coxeter groups from Coxeter matrices::

        sage: # needs sage.libs.gap sage.rings.number_field
        sage: W = CoxeterGroup([[1, 6, 3], [6, 1, 10], [3, 10, 1]]); W
        Coxeter group over Universal Cyclotomic Field with Coxeter matrix:
        [ 1  6  3]
        [ 6  1 10]
        [ 3 10  1]
        sage: W.gens()
        (
        [                 -1 -E(12)^7 + E(12)^11                   1]
        [                  0                   1                   0]
        [                  0                   0                   1],
        <BLANKLINE>
        [                  1                   0                   0]
        [-E(12)^7 + E(12)^11                  -1     E(20) - E(20)^9]
        [                  0                   0                   1],
        <BLANKLINE>
        [              1               0               0]
        [              0               1               0]
        [              1 E(20) - E(20)^9              -1]
        )

        sage: m = matrix([[1,3,3,3], [3,1,3,2], [3,3,1,2], [3,2,2,1]])
        sage: W = CoxeterGroup(m)
        sage: W.gens()
        (
        [-1  1  1  1]  [ 1  0  0  0]  [ 1  0  0  0]  [ 1  0  0  0]
        [ 0  1  0  0]  [ 1 -1  1  0]  [ 0  1  0  0]  [ 0  1  0  0]
        [ 0  0  1  0]  [ 0  0  1  0]  [ 1  1 -1  0]  [ 0  0  1  0]
        [ 0  0  0  1], [ 0  0  0  1], [ 0  0  0  1], [ 1  0  0 -1]
        )
        sage: a,b,c,d = W.gens()
        sage: (a*b*c)^3
        [ 5  1 -5  7]
        [ 5  0 -4  5]
        [ 4  1 -4  4]
        [ 0  0  0  1]
        sage: (a*b)^3
        [1 0 0 0]
        [0 1 0 0]
        [0 0 1 0]
        [0 0 0 1]
        sage: b*d == d*b
        True
        sage: a*c*a == c*a*c
        True

    We can create the matrix representation over different base rings and with
    different index sets. Note that the base ring must contain all
    `2*\cos(\pi/m_{ij})` where `(m_{ij})_{ij}` is the Coxeter matrix::

        sage: W = CoxeterGroup(m, base_ring=RR, index_set=['a','b','c','d'])
        sage: W.base_ring()
        Real Field with 53 bits of precision
        sage: W.index_set()
        ('a', 'b', 'c', 'd')

        sage: CoxeterGroup(m, base_ring=ZZ)
        Coxeter group over Integer Ring with Coxeter matrix:
        [1 3 3 3]
        [3 1 3 2]
        [3 3 1 2]
        [3 2 2 1]
        sage: CoxeterGroup([[1,4],[4,1]], base_ring=QQ)                                 # needs sage.symbolic
        Traceback (most recent call last):
        ...
        TypeError: unable to convert sqrt(2) to a rational

    Using the well-known conversion between Coxeter matrices and Coxeter
    graphs, we can input a Coxeter graph. Following the standard convention,
    edges with no label (i.e. labelled by ``None``) are treated as 3::

        sage: # needs sage.libs.gap sage.rings.number_field
        sage: G = Graph([(0,3,None), (1,3,15), (2,3,7), (0,1,3)])
        sage: W = CoxeterGroup(G); W
        Coxeter group over Universal Cyclotomic Field with Coxeter matrix:
        [ 1  3  2  3]
        [ 3  1  2 15]
        [ 2  2  1  7]
        [ 3 15  7  1]
        sage: G2 = W.coxeter_diagram()
        sage: CoxeterGroup(G2) is W
        True

    Because there currently is no class for `\ZZ \cup \{ \infty \}`, labels
    of `\infty` are given by `-1` in the Coxeter matrix::

        sage: # needs sage.libs.gap sage.rings.number_field
        sage: G = Graph([(0,1,None), (1,2,4), (0,2,oo)])
        sage: W = CoxeterGroup(G)
        sage: W.coxeter_matrix()
        [ 1  3 -1]
        [ 3  1  4]
        [-1  4  1]

    We can also create Coxeter groups from Cartan types using the
    ``implementation`` keyword::

        sage: W = CoxeterGroup(['D',5], implementation='reflection'); W
        Finite Coxeter group over Integer Ring with Coxeter matrix:
        [1 3 2 2 2]
        [3 1 3 2 2]
        [2 3 1 3 3]
        [2 2 3 1 2]
        [2 2 3 2 1]
        sage: W = CoxeterGroup(['H',3], implementation='reflection'); W                 # needs sage.libs.gap sage.rings.number_field
        Finite Coxeter group over
         Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
         with Coxeter matrix:
         [1 3 2]
         [3 1 5]
         [2 5 1]
    """
    @staticmethod
    def __classcall_private__(cls, data, base_ring=None, index_set=None):
        """
        Normalize arguments to ensure a unique representation.

        EXAMPLES::

            sage: W1 = CoxeterGroup(['A',2], implementation='reflection', base_ring=ZZ)
            sage: W2 = CoxeterGroup([[1,3],[3,1]], index_set=(1,2))
            sage: W1 is W2
            True
            sage: G1 = Graph([(1,2)])
            sage: W3 = CoxeterGroup(G1)
            sage: W1 is W3
            True
            sage: G2 = Graph([(1,2,3)])
            sage: W4 = CoxeterGroup(G2)
            sage: W1 is W4
            True
        """
        data = CoxeterMatrix(data, index_set=index_set)

        if base_ring is None:
            if data.is_simply_laced():
                base_ring = ZZ
            elif data.is_finite():
                from sage.rings.number_field.number_field import QuadraticField
                letter = data.coxeter_type().cartan_type().type()
                if letter in ['B', 'C', 'F']:
                    base_ring = QuadraticField(2)
                elif letter == 'G':
                    base_ring = QuadraticField(3)
                elif letter == 'H':
                    base_ring = QuadraticField(5)
                else:
                    from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField
                    base_ring = UniversalCyclotomicField()
            else:
                from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField
                base_ring = UniversalCyclotomicField()
        return super().__classcall__(cls, data, base_ring, data.index_set())

    def __init__(self, coxeter_matrix, base_ring, index_set):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup([[1,3,2],[3,1,3],[2,3,1]])
            sage: TestSuite(W).run()  # long time

            sage: # long time, needs sage.rings.number_field sage.symbolic
            sage: W = CoxeterGroup([[1,3,2],[3,1,4],[2,4,1]], base_ring=QQbar)
            sage: TestSuite(W).run()
            sage: W = CoxeterGroup([[1,3,2],[3,1,6],[2,6,1]])
            sage: TestSuite(W).run(max_runs=30)
            sage: W = CoxeterGroup([[1,3,2],[3,1,-1],[2,-1,1]])
            sage: TestSuite(W).run(max_runs=30)

        We check that :issue:`16630` is fixed::

            sage: CoxeterGroup(['D',4], base_ring=QQ).category()
            Category of finite irreducible Coxeter groups

            sage: # needs sage.rings.number_field
            sage: CoxeterGroup(['H',4], base_ring=QQbar).category()
            Category of finite irreducible Coxeter groups
            sage: F = CoxeterGroups().Finite()
            sage: all(CoxeterGroup([letter,i]) in F
            ....:     for i in range(2,5) for letter in ['A','B','D'])
            True
            sage: all(CoxeterGroup(['E',i]) in F for i in range(6,9))
            True
            sage: CoxeterGroup(['F',4]).category()
            Category of finite irreducible Coxeter groups
            sage: CoxeterGroup(['G',2]).category()
            Category of finite irreducible Coxeter groups
            sage: all(CoxeterGroup(['H',i]) in F for i in range(3,5))
            True
            sage: all(CoxeterGroup(['I',i]) in F for i in range(2,5))
            True
        """
        self._matrix = coxeter_matrix
        n = coxeter_matrix.rank()
        # Compute the matrix with entries `2 \cos( \pi / m_{ij} )`.
        MS = MatrixSpace(base_ring, n, sparse=True)
        one = MS.one()
        # FIXME: Hack because there is no ZZ \cup \{ \infty \}: -1 represents \infty
        if isinstance(base_ring, sage.rings.abc.UniversalCyclotomicField):
            E = base_ring.gen

            def val(x):
                if x == -1:
                    return 2
                else:
                    return E(2 * x) + ~E(2 * x)
        elif isinstance(base_ring, sage.rings.abc.NumberField_quadratic):
            from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField

            E = UniversalCyclotomicField().gen

            def val(x):
                if x == -1:
                    return 2
                else:
                    return base_ring((E(2 * x) + ~E(2 * x)).to_cyclotomic_field())
        else:
            def val(x):
                if x == -1:
                    return 2
                elif x == 1:
                    return -2
                elif x == 2:
                    return 0
                elif x == 3:
                    return 1
                else:
                    from sage.functions.trig import cos
                    from sage.symbolic.constants import pi
                    return base_ring(2 * cos(pi / x))
        gens = [one + MS([SparseEntry(i, j, val(coxeter_matrix[index_set[i], index_set[j]]))
                          for j in range(n)])
                for i in range(n)]
        # Make the generators dense matrices for consistency and speed
        gens = [g.dense_matrix() for g in gens]
        category = CoxeterGroups()
        # Now we shall see if the group is finite, and, if so, refine
        # the category to ``category.Finite()``. Otherwise the group is
        # infinite and we refine the category to ``category.Infinite()``.
        if self._matrix.is_finite():
            category = category.Finite()
        else:
            category = category.Infinite()
        if all(self._matrix._matrix[i, j] == 2
               for i in range(n) for j in range(i)):
            category = category.Commutative()
        if self._matrix.is_irreducible():
            category = category.Irreducible()
        self._index_set_inverse = {i: ii
                                   for ii, i in enumerate(self._matrix.index_set())}
        FinitelyGeneratedMatrixGroup_generic.__init__(self, ZZ(n), base_ring,
                                                      gens, category=category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: CoxeterGroup([[1,3,2],[3,1,4],[2,4,1]])                               # needs sage.libs.gap sage.rings.number_field
            Finite Coxeter group over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? with Coxeter matrix:
            [1 3 2]
            [3 1 4]
            [2 4 1]
        """
        rep = "Finite " if self.is_finite() else ""
        rep += "Coxeter group over {} with Coxeter matrix:\n{}".format(self.base_ring(), self._matrix)
        return rep

    def _coerce_map_from_(self, P):
        """
        Return ``True`` if ``P`` is a Coxeter group of the same
        Coxeter type and ``False`` otherwise.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: W = CoxeterGroup(["A",4])
            sage: W2 = WeylGroup(["A",4])
            sage: W._coerce_map_from_(W2)
            True
            sage: W3 = WeylGroup(["A",4], implementation='permutation')
            sage: W._coerce_map_from_(W3)
            True
            sage: W4 = WeylGroup(["A",3])
            sage: W.has_coerce_map_from(W4)
            False
        """
        if P in CoxeterGroups() and P.coxeter_type() is self.coxeter_type():
            return True
        return super()._coerce_map_from_(P)

    def coxeter_matrix(self):
        """
        Return the Coxeter matrix of ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup([[1,3],[3,1]])
            sage: W.coxeter_matrix()
            [1 3]
            [3 1]
            sage: W = CoxeterGroup(['H',3])                                             # needs sage.libs.gap sage.rings.number_field
            sage: W.coxeter_matrix()                                                    # needs sage.libs.gap sage.rings.number_field
            [1 3 2]
            [3 1 5]
            [2 5 1]
        """
        return self._matrix

    def bilinear_form(self):
        r"""
        Return the bilinear form associated to ``self``.

        Given a Coxeter group `G` with Coxeter matrix `M = (m_{ij})_{ij}`,
        the associated bilinear form `A = (a_{ij})_{ij}` is given by

        .. MATH::

            a_{ij} = -\cos\left( \frac{\pi}{m_{ij}} \right).

        If `A` is positive definite, then `G` is of finite type (and so
        the associated Coxeter group is a finite group). If `A` is
        positive semidefinite, then `G` is affine type.

        EXAMPLES::

            sage: W = CoxeterGroup(['D',4])
            sage: W.bilinear_form()                                                     # needs sage.symbolic
            [   1 -1/2    0    0]
            [-1/2    1 -1/2 -1/2]
            [   0 -1/2    1    0]
            [   0 -1/2    0    1]
        """
        return self._matrix.bilinear_form(self.base_ring().fraction_field())

    def is_finite(self):
        """
        Return ``True`` if this group is finite.

        EXAMPLES::

            sage: # needs sage.libs.gap sage.rings.number_field
            sage: [l for l in range(2, 9) if
            ....:  CoxeterGroup([[1,3,2],[3,1,l],[2,l,1]]).is_finite()]
            [2, 3, 4, 5]
            sage: [l for l in range(2, 9) if
            ....:  CoxeterGroup([[1,3,2,2],[3,1,l,2],[2,l,1,3],[2,2,3,1]]).is_finite()]
            [2, 3, 4]
            sage: [l for l in range(2, 9) if
            ....:  CoxeterGroup([[1,3,2,2,2], [3,1,3,3,2], [2,3,1,2,2],
            ....:                [2,3,2,1,l], [2,2,2,l,1]]).is_finite()]
            [2, 3]
            sage: [l for l in range(2, 9) if
            ....:  CoxeterGroup([[1,3,2,2,2], [3,1,2,3,3], [2,2,1,l,2],
            ....:                [2,3,l,1,2], [2,3,2,2,1]]).is_finite()]
            [2, 3]
            sage: [l for l in range(2, 9) if
            ....:  CoxeterGroup([[1,3,2,2,2,2], [3,1,l,2,2,2], [2,l,1,3,l,2],
            ....:                [2,2,3,1,2,2], [2,2,l,2,1,3], [2,2,2,2,3,1]]).is_finite()]
            [2, 3]
        """
        # Finite Coxeter groups are marked as finite in
        # their ``__init__`` method, so we can just check
        # the category of ``self``.
        return "Finite" in self.category().axioms()

    def is_commutative(self):
        """
        Return whether ``self`` is commutative.

        EXAMPLES::

            sage: CoxeterGroup(['A', 2]).is_commutative()
            False
            sage: W = CoxeterGroup(['I',2])
            sage: W.is_commutative()
            True

        TESTS::

            sage: CoxeterGroup([['A', 2], ['A', 1]]).is_commutative()
            False
            sage: CoxeterGroup([['A', 1]] * 3).is_commutative()
            True
        """
        return "Commutative" in self.category().axioms()

    @cached_method
    def order(self):
        """
        Return the order of ``self``.

        If the Coxeter group is finite, this uses an iterator.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: W = CoxeterGroup([[1,3],[3,1]])
            sage: W.order()
            6
            sage: W = CoxeterGroup([[1,-1],[-1,1]])                                     # needs sage.libs.gap
            sage: W.order()                                                             # needs sage.libs.gap
            +Infinity
        """
        if self.is_finite():
            return len(self)
        return infinity

    def canonical_representation(self):
        r"""
        Return the canonical faithful representation of ``self``, which
        is ``self``.

        EXAMPLES::

            sage: W = CoxeterGroup([[1,3],[3,1]])
            sage: W.canonical_representation() is W
            True
        """
        return self

    def simple_reflection(self, i):
        """
        Return the simple reflection `s_i`.

        INPUT:

        - ``i`` -- an element from the index set

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], implementation='reflection')
            sage: W.simple_reflection(1)
            [-1  1  0]
            [ 0  1  0]
            [ 0  0  1]
            sage: W.simple_reflection(2)
            [ 1  0  0]
            [ 1 -1  1]
            [ 0  0  1]
            sage: W.simple_reflection(3)
            [ 1  0  0]
            [ 0  1  0]
            [ 0  1 -1]
        """
        return self.gen(self._index_set_inverse[i])

    @cached_method
    def _positive_roots_reflections(self):
        """
        Return a family whose keys are the positive roots
        and values are the reflections.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 2])
            sage: F = W._positive_roots_reflections()
            sage: F.keys()
            [(1, 0), (1, 1), (0, 1)]
            sage: list(F)
            [
            [-1  1]  [ 0 -1]  [ 1  0]
            [ 0  1], [-1  0], [ 1 -1]
            ]
        """
        if not self.is_finite():
            raise NotImplementedError('not available for infinite groups')

        word = self.long_element(as_word=True)
        N = len(word)

        from sage.modules.free_module import FreeModule
        simple_roots = FreeModule(self.base_ring(), self.ngens()).gens()

        refls = self.simple_reflections()
        resu = []
        d = {}
        for i in range(1, N + 1):
            segment = word[:i]
            last = segment.pop()
            ref = refls[last]
            rt = simple_roots[last - 1]
            while segment:
                last = segment.pop()
                cr = refls[last]
                ref = cr * ref * cr
                rt = refls[last] * rt
            rt.set_immutable()
            resu += [rt]
            d[rt] = ref
        return Family(resu, lambda rt: d[rt])

    def positive_roots(self):
        """
        Return the positive roots.

        These are roots in the Coxeter sense, that all have the
        same norm. They are given by their coefficients in the
        base of simple roots, also taken to have all the same
        norm.

        .. SEEALSO::

            :meth:`reflections`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], implementation='reflection')
            sage: W.positive_roots()
            ((1, 0, 0), (1, 1, 0), (0, 1, 0), (1, 1, 1), (0, 1, 1), (0, 0, 1))

            sage: # needs sage.libs.gap sage.rings.number_field
            sage: W = CoxeterGroup(['I',5], implementation='reflection')
            sage: W.positive_roots()
            ((1, 0),
             (-E(5)^2 - E(5)^3, 1),
             (-E(5)^2 - E(5)^3, -E(5)^2 - E(5)^3),
             (1, -E(5)^2 - E(5)^3),
             (0, 1))
        """
        return tuple(self._positive_roots_reflections().keys())

    def reflections(self):
        """
        Return the set of reflections.

        The order is the one given by :meth:`positive_roots`.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',2], implementation='reflection')
            sage: list(W.reflections())
            [
            [-1  1]  [ 0 -1]  [ 1  0]
            [ 0  1], [-1  0], [ 1 -1]
            ]
        """
        return self._positive_roots_reflections()

    @cached_method
    def roots(self):
        """
        Return the roots.

        These are roots in the Coxeter sense, that all have the
        same norm. They are given by their coefficients in the
        base of simple roots, also taken to have all the same
        norm.

        The positive roots are listed first, then the negative roots
        in the same order. The order is the one given by :meth:`roots`.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], implementation='reflection')
            sage: W.roots()
            ((1, 0, 0),
             (1, 1, 0),
             (0, 1, 0),
             (1, 1, 1),
             (0, 1, 1),
             (0, 0, 1),
             (-1, 0, 0),
             (-1, -1, 0),
             (0, -1, 0),
             (-1, -1, -1),
             (0, -1, -1),
             (0, 0, -1))

            sage: # needs sage.libs.gap sage.rings.number_field
            sage: W = CoxeterGroup(['I',5], implementation='reflection')
            sage: len(W.roots())
            10
        """
        if not self.is_finite():
            raise NotImplementedError('not available for infinite groups')
        positive = self.positive_roots()
        return positive + tuple([-v for v in positive])

    def simple_root_index(self, i):
        r"""
        Return the index of the simple root `\alpha_i`.

        This is the position of `\alpha_i` in the list of all roots
        as given be :meth:`roots`.

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], implementation='reflection')
            sage: [W.simple_root_index(i) for i in W.index_set()]
            [0, 2, 5]
        """
        roots = self.roots()
        rt = roots[0].parent().gen(self._index_set_inverse[i])
        return roots.index(rt)

    @cached_method
    def fundamental_weights(self):
        """
        Return the fundamental weights for ``self``.

        This is the dual basis to the basis of simple roots.

        The base ring must be a field.

        .. SEEALSO:: :meth:`fundamental_weight`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], implementation='reflection')
            sage: W.fundamental_weights()                                               # needs sage.symbolic
            Finite family {1: (3/2, 1, 1/2), 2: (1, 2, 1), 3: (1/2, 1, 3/2)}
        """
        simple_weights = self.bilinear_form().inverse()
        I = self.index_set()
        D = {i: simple_weights[k] for k, i in enumerate(I)}
        return Family(I, D.__getitem__)

    def fundamental_weight(self, i):
        r"""
        Return the fundamental weight with index ``i``.

        .. SEEALSO:: :meth:`fundamental_weights`

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], implementation='reflection')
            sage: W.fundamental_weight(1)                                               # needs sage.symbolic
            (3/2, 1, 1/2)
        """
        return self.fundamental_weights()[i]

    class Element(MatrixGroupElement_generic):
        """
        A Coxeter group element.
        """
        def first_descent(self, side='right', index_set=None, positive=False):
            """
            Return the first left (resp. right) descent of ``self``, as
            ane element of ``index_set``, or ``None`` if there is none.

            See :meth:`descents` for a description of the options.

            EXAMPLES::

                sage: W = CoxeterGroup(['A',3], implementation='reflection')
                sage: a,b,c = W.gens()
                sage: elt = b*a*c
                sage: elt.first_descent()
                1
                sage: elt.first_descent(side='left')
                2
            """
            M = self.matrix()
            if side != 'right':
                M = ~M
            I = self.parent().index_set()
            n = len(I)
            zero = M.base_ring().zero()
            if index_set is None:
                index_set = range(n)
            else:
                I_inv = self.parent()._index_set_inverse
                index_set = [I_inv[i] for i in index_set]
            if positive:
                for i in index_set:
                    if not _matrix_test_right_descent(M, i, n, zero):
                        return I[i]
            else:
                for i in index_set:
                    if _matrix_test_right_descent(M, i, n, zero):
                        return I[i]
            return None

        def descents(self, side='right', index_set=None, positive=False):
            """
            Return the descents of ``self``, as a list of elements of the
            ``index_set``.

            INPUT:

            - ``index_set`` -- (default: all of them) a subset (as a list
              or iterable) of the nodes of the Dynkin diagram
            - ``side`` -- (default: ``'right'``) ``'left'`` or ``'right'``
            - ``positive`` -- boolean (default: ``False``)

            EXAMPLES::

                sage: W = CoxeterGroup(['A',3], implementation='reflection')
                sage: a,b,c = W.gens()
                sage: elt = b*a*c
                sage: elt.descents()
                [1, 3]
                sage: elt.descents(positive=True)
                [2]
                sage: elt.descents(index_set=[1,2])
                [1]
                sage: elt.descents(side='left')
                [2]
            """
            M = self.matrix()
            if side != 'right':
                M = ~M
            I = self.parent().index_set()
            n = len(I)
            zero = M.base_ring().zero()
            if index_set is None:
                index_set = range(n)
            else:
                I_inv = self.parent()._index_set_inverse
                index_set = [I_inv[i] for i in index_set]
            if positive:
                return [I[i] for i in index_set if not _matrix_test_right_descent(M, i, n, zero)]
            return [I[i] for i in index_set if _matrix_test_right_descent(M, i, n, zero)]

        def has_right_descent(self, i):
            r"""
            Return whether ``i`` is a right descent of ``self``.

            A Coxeter system `(W, S)` has a root system defined as
            `\{ w(\alpha_s) \}_{w \in W}` and we define the positive
            (resp. negative) roots `\alpha = \sum_{s \in S} c_s \alpha_s`
            by all `c_s \geq 0` (resp. `c_s \leq 0`). In particular, we note
            that if `\ell(w s) > \ell(w)` then `w(\alpha_s) > 0` and if
            `\ell(ws) < \ell(w)` then `w(\alpha_s) < 0`.
            Thus `i \in I` is a right descent if `w(\alpha_{s_i}) < 0`
            or equivalently if the matrix representing `w` has all entries
            of the `i`-th column being nonpositive.

            INPUT:

            - ``i`` -- an element in the index set

            EXAMPLES::

                sage: W = CoxeterGroup(['A',3], implementation='reflection')
                sage: a,b,c = W.gens()
                sage: elt = b*a*c
                sage: [elt.has_right_descent(i) for i in [1, 2, 3]]
                [True, False, True]
            """
            i = self.parent()._index_set_inverse[i]
            n = len(self.parent().index_set())
            M = self.matrix()
            zero = M.base_ring().zero()
            return _matrix_test_right_descent(M, i, n, zero)

        def canonical_matrix(self):
            r"""
            Return the matrix of ``self`` in the canonical faithful
            representation, which is ``self`` as a matrix.

            EXAMPLES::

                sage: W = CoxeterGroup(['A',3], implementation='reflection')
                sage: a,b,c = W.gens()
                sage: elt = a*b*c
                sage: elt.canonical_matrix()
                [ 0  0 -1]
                [ 1  0 -1]
                [ 0  1 -1]
            """
            return self.matrix()

        @cached_method
        def action_on_root_indices(self, i, side='left'):
            """
            Return the action on the set of roots.

            The roots are ordered as in the output of the method :meth:`roots`.

            EXAMPLES::

                sage: W = CoxeterGroup(['A',3], implementation='reflection')
                sage: w = W.w0
                sage: w.action_on_root_indices(0)
                11
            """
            if side == "left":
                w = self
            elif side == "right":
                w = ~self
            else:
                raise ValueError('side must be "left" or "right"')
            roots = self.parent().roots()
            rt = w * roots[i]
            return roots.index(rt)


def _matrix_test_right_descent(M, i, n, zero):
    """
    Test if the matrix ``M`` has a right ``i``-descent.

    INPUT:

    - ``M`` -- the matrix
    - ``i`` -- the index
    - ``n`` -- the size of the matrix
    - ``zero`` -- the zero element in the base ring of ``M``

    .. NOTE::

        This is a helper function for :class:`CoxeterMatrixGroup.Element`
        and optimized for speed. Specifically, it is called often and
        there is no need to recompute ``n`` (and ``zero``) each time this
        function is called.

    .. TODO::

        Cythonize this function.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.coxeter_group import _matrix_test_right_descent
        sage: W = CoxeterGroup(['A',3], implementation='reflection')
        sage: a,b,c = W.gens()
        sage: elt = b*a*c
        sage: zero = W.base_ring().zero()
        sage: [_matrix_test_right_descent(elt.matrix(), i, 3, zero)
        ....:  for i in range(3)]
        [True, False, True]
    """
    for j in range(n):
        c = M[j, i]
        if c < zero:
            return True
        elif c > zero:
            return False
    raise AssertionError('a zero column, so there must be a bug')
