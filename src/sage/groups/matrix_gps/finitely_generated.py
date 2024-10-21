# sage_setup: distribution = sagemath-modules
"""
Finitely Generated Matrix Groups

This class is designed for computing with matrix groups defined by a
finite set of generating matrices.

EXAMPLES::

    sage: F = GF(3)
    sage: gens = [matrix(F, 2, [1,0, -1,1]), matrix(F, 2, [1,1,0,1])]
    sage: G = MatrixGroup(gens)
    sage: G.conjugacy_classes_representatives()                                         # needs sage.libs.gap
    (
    [1 0]  [0 2]  [0 1]  [2 0]  [0 2]  [0 1]  [0 2]
    [0 1], [1 1], [2 1], [0 2], [1 2], [2 2], [1 0]
    )

The finitely generated matrix groups can also be constructed as
subgroups of matrix groups::

    sage: SL2Z = SL(2, ZZ)
    sage: S, T = SL2Z.gens()                                                            # needs sage.libs.gap
    sage: SL2Z.subgroup([T^2])                                                          # needs sage.libs.gap
    Subgroup with 1 generators (
    [1 2]
    [0 1]
    ) of Special Linear Group of degree 2 over Integer Ring

AUTHORS:

- William Stein: initial version

- David Joyner (2006-03-15): degree, base_ring, _contains_, list,
  random, order methods; examples

- William Stein (2006-12): rewrite

- David Joyner (2007-12): Added invariant_generators (with Martin
  Albrecht and Simon King)

- David Joyner (2008-08): Added module_composition_factors (interface
  to GAP's MeatAxe implementation) and as_permutation_group (returns
  isomorphic PermutationGroup).

- Simon King (2010-05): Improve invariant_generators by using GAP
  for the construction of the Reynolds operator in Singular.

- Volker Braun (2013-1) port to new Parent, libGAP.

- Sebastian Oehms (2018-07): Added _permutation_group_element_ (Issue #25706)
- Sebastian Oehms (2019-01): Revision of :issue:`25706` (:issue:`26903` and :issue:`27143`).
"""

# #############################################################################
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#                     2009 Mike Hansen
#                     2012 Rob Beezer
#                     2013 Volker Braun <vbraun.name@gmail.com>
#                     2013 Nathann Cohen
#                     2013 Travis Scrimshaw
#                     2017 Peter Bruin
#                     2021 Frédéric Chapoton
#                     2023 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# #############################################################################

from sage.groups.matrix_gps.group_element import is_MatrixGroupElement
from sage.groups.matrix_gps.matrix_group import MatrixGroup_generic
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.structure.element import Matrix
from sage.structure.sequence import Sequence


def normalize_square_matrices(matrices):
    """
    Find a common space for all matrices.

    OUTPUT: list of matrices, all elements of the same matrix space

    EXAMPLES::

        sage: from sage.groups.matrix_gps.finitely_generated import normalize_square_matrices
        sage: m1 = [[1,2], [3,4]]
        sage: m2 = [2, 3, 4, 5]
        sage: m3 = matrix(QQ, [[1/2,1/3], [1/4,1/5]])
        sage: m4 = MatrixGroup(m3).gen(0)
        sage: normalize_square_matrices([m1, m2, m3, m4])
        [
        [1 2]  [2 3]  [1/2 1/3]  [1/2 1/3]
        [3 4], [4 5], [1/4 1/5], [1/4 1/5]
        ]
    """
    deg = []
    gens = []
    for m in matrices:
        if is_MatrixGroupElement(m):
            deg.append(m.parent().degree())
            gens.append(m.matrix())
            continue
        if isinstance(m, Matrix):
            if not m.is_square():
                raise TypeError('matrix must be square')
            deg.append(m.ncols())
            gens.append(m)
            continue
        try:
            m = list(m)
        except TypeError:
            gens.append(m)
            continue
        if isinstance(m[0], (list, tuple)):
            m = [list(_) for _ in m]
            degree = ZZ(len(m))
        else:
            degree, rem = ZZ(len(m)).sqrtrem()
            if rem != 0:
                raise ValueError('list of plain numbers must have square integer length')
        deg.append(degree)
        gens.append(matrix(degree, degree, m))
    deg = set(deg)
    if len(set(deg)) != 1:
        raise ValueError('not all matrices have the same size')
    gens = Sequence(gens, immutable=True)
    MS = gens.universe()
    if not isinstance(MS, MatrixSpace):
        raise TypeError('all generators must be matrices')
    if MS.nrows() != MS.ncols():
        raise ValueError('matrices must be square')
    return gens


def QuaternionMatrixGroupGF3():
    r"""
    The quaternion group as a set of `2\times 2` matrices over `\GF{3}`.

    OUTPUT:

    A matrix group consisting of `2\times 2` matrices with
    elements from the finite field of order 3.  The group is
    the quaternion group, the nonabelian group of order 8 that
    is not isomorphic to the group of symmetries of a square
    (the dihedral group `D_4`).

    .. NOTE::

        This group is most easily available via ``groups.matrix.QuaternionGF3()``.

    EXAMPLES:

    The generators are the matrix representations of the
    elements commonly called `I` and `J`, while `K`
    is the product of `I` and `J`. ::

        sage: # needs sage.libs.gap
        sage: from sage.groups.matrix_gps.finitely_generated import QuaternionMatrixGroupGF3

        sage: # needs sage.libs.gap
        sage: Q = QuaternionMatrixGroupGF3()
        sage: Q.order()
        8
        sage: aye = Q.gens()[0]; aye
        [1 1]
        [1 2]
        sage: jay = Q.gens()[1]; jay
        [2 1]
        [1 1]
        sage: kay = aye*jay; kay
        [0 2]
        [1 0]

    TESTS::

        sage: groups.matrix.QuaternionGF3()
        Matrix group over Finite Field of size 3 with 2 generators (
        [1 1]  [2 1]
        [1 2], [1 1]
        )

        sage: # needs sage.groups
        sage: Q = QuaternionMatrixGroupGF3()
        sage: QP = Q.as_permutation_group()
        sage: QP.is_isomorphic(QuaternionGroup())
        True
        sage: H = DihedralGroup(4)
        sage: H.order()
        8
        sage: QP.is_abelian(), H.is_abelian()
        (False, False)
        sage: QP.is_isomorphic(H)
        False
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    from sage.matrix.matrix_space import MatrixSpace
    MS = MatrixSpace(FiniteField(3), 2)
    aye = MS([1,1,1,2])
    jay = MS([2,1,1,1])
    return MatrixGroup([aye, jay])


def MatrixGroup(*gens, **kwds):
    r"""
    Return the matrix group with given generators.

    INPUT:

    - ``*gens`` -- matrices, or a single list/tuple/iterable of
      matrices, or a matrix group

    - ``check`` -- boolean keyword argument (default: ``True``);
      whether to check that each matrix is invertible

    EXAMPLES::

        sage: F = GF(5)
        sage: gens = [matrix(F, 2, [1,2, -1,1]), matrix(F,2, [1,1, 0,1])]
        sage: G = MatrixGroup(gens); G
        Matrix group over Finite Field of size 5 with 2 generators (
        [1 2]  [1 1]
        [4 1], [0 1]
        )

    In the second example, the generators are a matrix over
    `\ZZ`, a matrix over a finite field, and the integer
    `2`. Sage determines that they both canonically map to
    matrices over the finite field, so creates that matrix group
    there::

        sage: gens = [matrix(2, [1,2, -1,1]), matrix(GF(7), 2, [1,1, 0,1]), 2]
        sage: G = MatrixGroup(gens); G
        Matrix group over Finite Field of size 7 with 3 generators (
        [1 2]  [1 1]  [2 0]
        [6 1], [0 1], [0 2]
        )

    Each generator must be invertible::

        sage: G = MatrixGroup([matrix(ZZ, 2, [1,2,3,4])])
        Traceback (most recent call last):
        ...
        ValueError: each generator must be an invertible matrix

        sage: F = GF(5); MS = MatrixSpace(F, 2, 2)
        sage: MatrixGroup([MS.0])
        Traceback (most recent call last):
        ...
        ValueError: each generator must be an invertible matrix
        sage: MatrixGroup([MS.0], check=False)  # works formally but is mathematical nonsense
        Matrix group over Finite Field of size 5 with 1 generators (
        [1 0]
        [0 0]
        )

    Some groups are not supported, or do not have much functionality
    implemented::

        sage: G = SL(0, QQ)
        Traceback (most recent call last):
        ...
        ValueError: the degree must be at least 1

        sage: SL2C = SL(2, CC);  SL2C
        Special Linear Group of degree 2 over Complex Field with 53 bits of precision
        sage: SL2C.gens()
        Traceback (most recent call last):
        ...
        AttributeError: 'LinearMatrixGroup_generic_with_category' object has no attribute 'gens'...
    """
    if isinstance(gens[-1], dict):   # hack for unpickling
        kwds.update(gens[-1])
        gens = gens[:-1]
    check = kwds.get('check', True)
    if len(gens) == 1:
        if isinstance(gens[0], (list, tuple)):
            gens = list(gens[0])
        else:
            try:
                gens = [g.matrix() for g in gens[0]]
            except AttributeError:
                pass
    if len(gens) == 0:
        raise ValueError('need at least one generator')
    gens = normalize_square_matrices(gens)
    if check and any(not g.is_invertible() for g in gens):
        raise ValueError('each generator must be an invertible matrix')
    MS = gens.universe()
    base_ring = MS.base_ring()
    degree = ZZ(MS.ncols())   # == MS.nrows()
    category = kwds.get('category', None)
    try:
        from sage.libs.gap.libgap import libgap
        from .finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
    except ImportError:
        pass
    else:
        try:
            gap_gens = [libgap(matrix_gen) for matrix_gen in gens]
            gap_group = libgap.Group(gap_gens)
            return FinitelyGeneratedMatrixGroup_gap(degree, base_ring, gap_group,
                                                    category=category)
        except (TypeError, ValueError):
            pass

    return FinitelyGeneratedMatrixGroup_generic(degree, base_ring, gens,
                                                category=category)

###################################################################
#
# Matrix group over a generic ring
#
###################################################################


class FinitelyGeneratedMatrixGroup_generic(MatrixGroup_generic):
    """
    TESTS::

        sage: # needs sage.symbolic
        sage: m1 = matrix(SR, [[1,2], [3,4]])
        sage: m2 = matrix(SR, [[1,3], [-1,0]])
        sage: MatrixGroup(m1) == MatrixGroup(m1)
        True
        sage: MatrixGroup(m1) == MatrixGroup(m1.change_ring(QQ))
        False
        sage: MatrixGroup(m1) == MatrixGroup(m2)
        False
        sage: MatrixGroup(m1, m2) == MatrixGroup(m2, m1)
        False

        sage: m1 = matrix(QQ, [[1,2], [3,4]])
        sage: m2 = matrix(QQ, [[1,3], [-1,0]])
        sage: MatrixGroup(m1) == MatrixGroup(m1)
        True
        sage: MatrixGroup(m1) == MatrixGroup(m2)
        False
        sage: MatrixGroup(m1, m2) == MatrixGroup(m2, m1)
        False

        sage: # needs sage.libs.gap
        sage: G = GL(2, GF(3))
        sage: H = G.as_matrix_group()
        sage: H == G, G == H
        (True, True)
    """

    def __init__(self, degree, base_ring, generator_matrices, category=None):
        """
        Matrix group generated by a finite number of matrices.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: m1 = matrix(SR, [[1,2], [3,4]])
            sage: m2 = matrix(SR, [[1,3], [-1,0]])
            sage: G = MatrixGroup(m1, m2)
            sage: TestSuite(G).run()
            sage: type(G)
            <class 'sage.groups.matrix_gps.finitely_generated.FinitelyGeneratedMatrixGroup_generic_with_category'>

            sage: from sage.groups.matrix_gps.finitely_generated import \
            ....:     FinitelyGeneratedMatrixGroup_generic
            sage: G = FinitelyGeneratedMatrixGroup_generic(2, QQ, [matrix(QQ, [[1,2], [3,4]])])
            sage: G.gens()
            (
            [1 2]
            [3 4]
            )
        """
        self._gens_matrix = generator_matrices
        MatrixGroup_generic.__init__(self, degree, base_ring, category=category)

    @cached_method
    def gens(self) -> tuple:
        """
        Return the generators of the matrix group.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F, 2, 2)
            sage: gens = [MS([[1,0], [0,1]]), MS([[1,1], [0,1]])]
            sage: G = MatrixGroup(gens)
            sage: gens[0] in G
            True
            sage: gens = G.gens()
            sage: gens[0] in G
            True
            sage: gens = [MS([[1,0], [0,1]]), MS([[1,1], [0,1]])]

            sage: F = GF(5); MS = MatrixSpace(F, 2, 2)
            sage: G = MatrixGroup([MS(1), MS([1,2, 3,4])])
            sage: G
            Matrix group over Finite Field of size 5 with 2 generators (
            [1 0]  [1 2]
            [0 1], [3 4]
            )
            sage: G.gens()
            (
            [1 0]  [1 2]
            [0 1], [3 4]
            )
        """
        return tuple(self.element_class(self, x, check=False, convert=False)
                     for x in self._gens_matrix)

    def gen(self, i):
        """
        Return the `i`-th generator.

        OUTPUT: the `i`-th generator of the group

        EXAMPLES::

            sage: # needs sage.libs.gap
            sage: H = GL(2, GF(3))
            sage: h1, h2 = H([[1,0], [2,1]]), H([[1,1], [0,1]])
            sage: G = H.subgroup([h1, h2])
            sage: G.gen(0)
            [1 0]
            [2 1]
            sage: G.gen(0).matrix() == h1.matrix()
            True
        """
        return self.gens()[i]

    def ngens(self):
        """
        Return the number of generators.

        OUTPUT: integer; the number of generators

        EXAMPLES::

            sage: # needs sage.libs.gap
            sage: H = GL(2, GF(3))
            sage: h1, h2 = H([[1,0], [2,1]]), H([[1,1], [0,1]])
            sage: G = H.subgroup([h1, h2])
            sage: G.ngens()
            2
        """
        return len(self._gens_matrix)

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::

            sage: G = MatrixGroup([matrix(CC, [[1,2], [3,4]]),
            ....:                  matrix(CC, [[1,3], [-1,0]])])
            sage: loads(dumps(G)) == G
            True

        Check that :issue:`22128` is fixed::

            sage: # needs sage.symbolic
            sage: R = MatrixSpace(SR, 2)
            sage: G = MatrixGroup([R([[1, 1], [0, 1]])])
            sage: G.register_embedding(R)
            sage: loads(dumps(G))
            Matrix group over Symbolic Ring with 1 generators (
            [1 1]
            [0 1]
            )
        """
        return MatrixGroup, (self._gens_matrix, {'check': False})

    def _test_matrix_generators(self, **options):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: m1 = matrix(SR, [[1,2], [3,4]])
            sage: m2 = matrix(SR, [[1,3], [-1,0]])
            sage: G = MatrixGroup(m1, m2)
            sage: G._test_matrix_generators()
        """
        tester = self._tester(**options)
        for g,h in zip(self.gens(), MatrixGroup(self.gens()).gens()):
            tester.assertEqual(g.matrix(), h.matrix())
