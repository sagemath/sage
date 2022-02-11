"""
Finitely Generated Matrix Groups

This class is designed for computing with matrix groups defined by a
finite set of generating matrices.

EXAMPLES::

    sage: F = GF(3)
    sage: gens = [matrix(F,2, [1,0, -1,1]), matrix(F,2, [1,1,0,1])]
    sage: G = MatrixGroup(gens)
    sage: G.conjugacy_classes_representatives()
    (
    [1 0]  [0 2]  [0 1]  [2 0]  [0 2]  [0 1]  [0 2]
    [0 1], [1 1], [2 1], [0 2], [1 2], [2 2], [1 0]
    )

The finitely generated matrix groups can also be constructed as
subgroups of matrix groups::

    sage: SL2Z = SL(2,ZZ)
    sage: S, T = SL2Z.gens()
    sage: SL2Z.subgroup([T^2])
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

- Sebastian Oehms (2018-07): Added _permutation_group_element_ (Trac #25706)
- Sebastian Oehms (2019-01): Revision of :trac:`25706` (:trac:`26903` and :trac:`27143`).
"""

##############################################################################
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
##############################################################################

from sage.rings.integer_ring import ZZ
from sage.rings.all import QQbar
from sage.structure.element import is_Matrix
from sage.matrix.matrix_space import MatrixSpace, is_MatrixSpace
from sage.matrix.constructor import matrix
from sage.structure.sequence import Sequence
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.fraction_field import FractionField
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.number_field.number_field import CyclotomicField
from sage.combinat.integer_vector import IntegerVectors

from sage.groups.matrix_gps.matrix_group import (MatrixGroup_generic,
                                                 MatrixGroup_gap)
from sage.groups.matrix_gps.group_element import is_MatrixGroupElement


def normalize_square_matrices(matrices):
    """
    Find a common space for all matrices.

    OUTPUT:

    A list of matrices, all elements of the same matrix space.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.finitely_generated import normalize_square_matrices
        sage: m1 = [[1,2],[3,4]]
        sage: m2 = [2, 3, 4, 5]
        sage: m3 = matrix(QQ, [[1/2,1/3],[1/4,1/5]])
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
        if is_Matrix(m):
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
    if not is_MatrixSpace(MS):
        raise TypeError('all generators must be matrices')
    if MS.nrows() != MS.ncols():
        raise ValueError('matrices must be square')
    return gens


def QuaternionMatrixGroupGF3():
    r"""
    The quaternion group as a set of `2\times 2` matrices over `GF(3)`.

    OUTPUT:

    A matrix group consisting of `2\times 2` matrices with
    elements from the finite field of order 3.  The group is
    the quaternion group, the nonabelian group of order 8 that
    is not isomorphic to the group of symmetries of a square
    (the dihedral group `D_4`).

    .. note::
        This group is most easily available via ``groups.matrix.QuaternionGF3()``.

    EXAMPLES:

    The generators are the matrix representations of the
    elements commonly called `I` and `J`, while `K`
    is the product of `I` and `J`. ::

        sage: from sage.groups.matrix_gps.finitely_generated import QuaternionMatrixGroupGF3
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
      matrices, or a matrix group.

    - ``check`` -- boolean keyword argument (optional, default:
      ``True``). Whether to check that each matrix is invertible.

    EXAMPLES::

        sage: F = GF(5)
        sage: gens = [matrix(F,2,[1,2, -1, 1]), matrix(F,2, [1,1, 0,1])]
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

        sage: gens = [matrix(2,[1,2, -1, 1]), matrix(GF(7), 2, [1,1, 0,1]), 2]
        sage: G = MatrixGroup(gens); G
        Matrix group over Finite Field of size 7 with 3 generators (
        [1 2]  [1 1]  [2 0]
        [6 1], [0 1], [0 2]
        )

    Each generator must be invertible::

        sage: G = MatrixGroup([matrix(ZZ,2,[1,2,3,4])])
        Traceback (most recent call last):
        ...
        ValueError: each generator must be an invertible matrix

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
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
        AttributeError: 'LinearMatrixGroup_generic_with_category' object has no attribute 'gens'
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
    from sage.libs.gap.libgap import libgap
    category = kwds.get('category', None)
    try:
        gap_gens = [libgap(matrix_gen) for matrix_gen in gens]
        gap_group = libgap.Group(gap_gens)
        return FinitelyGeneratedMatrixGroup_gap(degree, base_ring, gap_group,
                                                category=category)
    except (TypeError, ValueError):
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

        sage: m1 = matrix(SR, [[1,2],[3,4]])
        sage: m2 = matrix(SR, [[1,3],[-1,0]])
        sage: MatrixGroup(m1) == MatrixGroup(m1)
        True
        sage: MatrixGroup(m1) == MatrixGroup(m1.change_ring(QQ))
        False
        sage: MatrixGroup(m1) == MatrixGroup(m2)
        False
        sage: MatrixGroup(m1, m2) == MatrixGroup(m2, m1)
        False

        sage: m1 = matrix(QQ, [[1,2],[3,4]])
        sage: m2 = matrix(QQ, [[1,3],[-1,0]])
        sage: MatrixGroup(m1) == MatrixGroup(m1)
        True
        sage: MatrixGroup(m1) == MatrixGroup(m2)
        False
        sage: MatrixGroup(m1, m2) == MatrixGroup(m2, m1)
        False

        sage: G = GL(2, GF(3))
        sage: H = G.as_matrix_group()
        sage: H == G, G == H
        (True, True)
    """

    def __init__(self, degree, base_ring, generator_matrices, category=None):
        """
        Matrix group generated by a finite number of matrices.

        EXAMPLES::

            sage: m1 = matrix(SR, [[1,2],[3,4]])
            sage: m2 = matrix(SR, [[1,3],[-1,0]])
            sage: G = MatrixGroup(m1, m2)
            sage: TestSuite(G).run()
            sage: type(G)
            <class 'sage.groups.matrix_gps.finitely_generated.FinitelyGeneratedMatrixGroup_generic_with_category'>

            sage: from sage.groups.matrix_gps.finitely_generated import \
            ....:     FinitelyGeneratedMatrixGroup_generic
            sage: G = FinitelyGeneratedMatrixGroup_generic(2, QQ, [matrix(QQ,[[1,2],[3,4]])])
            sage: G.gens()
            (
            [1 2]
            [3 4]
            )
        """
        self._gens_matrix = generator_matrices
        MatrixGroup_generic.__init__(self, degree, base_ring, category=category)

    @cached_method
    def gens(self):
        """
        Return the generators of the matrix group.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]), MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: gens[0] in G
            True
            sage: gens = G.gens()
            sage: gens[0] in G
            True
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
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
        Return the `i`-th generator

        OUTPUT:

        The `i`-th generator of the group.

        EXAMPLES::

            sage: H = GL(2, GF(3))
            sage: h1, h2 = H([[1,0],[2,1]]), H([[1,1],[0,1]])
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
        Return the number of generators

        OUTPUT:

        An integer. The number of generators.

        EXAMPLES::

            sage: H = GL(2, GF(3))
            sage: h1, h2 = H([[1,0],[2,1]]), H([[1,1],[0,1]])
            sage: G = H.subgroup([h1, h2])
            sage: G.ngens()
            2
        """
        return len(self._gens_matrix)

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::

            sage: G = MatrixGroup([matrix(CC, [[1,2],[3,4]]),
            ....:                  matrix(CC, [[1,3],[-1,0]])])
            sage: loads(dumps(G)) == G
            True

        Check that :trac:`22128` is fixed::

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

            sage: m1 = matrix(SR, [[1,2],[3,4]])
            sage: m2 = matrix(SR, [[1,3],[-1,0]])
            sage: G = MatrixGroup(m1, m2)
            sage: G._test_matrix_generators()
        """
        tester = self._tester(**options)
        for g,h in zip(self.gens(), MatrixGroup(self.gens()).gens()):
            tester.assertEqual(g.matrix(), h.matrix())

###################################################################
#
# Matrix group over a ring that GAP understands
#
###################################################################


class FinitelyGeneratedMatrixGroup_gap(MatrixGroup_gap):
    """
    Matrix group generated by a finite number of matrices.

    EXAMPLES::

        sage: m1 = matrix(GF(11), [[1,2],[3,4]])
        sage: m2 = matrix(GF(11), [[1,3],[10,0]])
        sage: G = MatrixGroup(m1, m2);  G
        Matrix group over Finite Field of size 11 with 2 generators (
        [1 2]  [ 1  3]
        [3 4], [10  0]
        )
        sage: type(G)
        <class 'sage.groups.matrix_gps.finitely_generated.FinitelyGeneratedMatrixGroup_gap_with_category'>
        sage: TestSuite(G).run()
    """

    def __reduce__(self):
        """
        Implement pickling.

        EXAMPLES::

            sage: m1 = matrix(QQ, [[1,2],[3,4]])
            sage: m2 = matrix(QQ, [[1,3],[-1,0]])
            sage: loads(MatrixGroup(m1, m2).dumps())
            Matrix group over Rational Field with 2 generators (
            [1 2]  [ 1  3]
            [3 4], [-1  0]
            )
        """
        return (MatrixGroup,
                tuple(g.matrix() for g in self.gens()) + ({'check':False},))

    def as_permutation_group(self, algorithm=None, seed=None):
        r"""
        Return a permutation group representation for the group.

        In most cases occurring in practice, this is a permutation
        group of minimal degree (the degree being determined from
        orbits under the group action). When these orbits are hard to
        compute, the procedure can be time-consuming and the degree
        may not be minimal.

        INPUT:

        - ``algorithm`` -- ``None`` or ``'smaller'``. In the latter
          case, try harder to find a permutation representation of
          small degree.
        - ``seed`` -- ``None`` or an integer specifying the seed
          to fix results depending on pseudo-random-numbers. Here
          it makes sense to be used with respect to the ``'smaller'``
          option, since gap produces random output in that context.

        OUTPUT:

        A permutation group isomorphic to ``self``. The
        ``algorithm='smaller'`` option tries to return an isomorphic
        group of low degree, but is not guaranteed to find the
        smallest one and must not even differ from the one obtained
        without the option. In that case repeating the invocation
        may help (see the example below).

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2), 5, 5)
            sage: A = MS([[0,0,0,0,1],[0,0,0,1,0],[0,0,1,0,0],[0,1,0,0,0],[1,0,0,0,0]])
            sage: G = MatrixGroup([A])
            sage: G.as_permutation_group().order()
            2

        A finite subgroup of  GL(12,Z) as a permutation group::

            sage: imf = libgap.function_factory('ImfMatrixGroup')
            sage: GG = imf( 12, 3 )
            sage: G = MatrixGroup(GG.GeneratorsOfGroup())
            sage: G.cardinality()
            21499084800
            sage: P = G.as_permutation_group()
            sage: Psmaller = G.as_permutation_group(algorithm="smaller", seed=6)
            sage: P == Psmaller  # see the note below
            True
            sage: Psmaller = G.as_permutation_group(algorithm="smaller")
            sage: P == Psmaller
            False
            sage: P.cardinality()
            21499084800
            sage: P.degree()
            144
            sage: Psmaller.cardinality()
            21499084800
            sage: Psmaller.degree()
            80

        .. NOTE::

            In this case, the "smaller" option returned an isomorphic
            group of lower degree. The above example used GAP's library
            of irreducible maximal finite ("imf") integer matrix groups
            to construct the MatrixGroup G over GF(7). The section
            "Irreducible Maximal Finite Integral Matrix Groups" in the
            GAP reference manual has more details.

        .. NOTE::

            Concerning the option ``algorithm='smaller'`` you should note
            the following from GAP documentation: "The methods used might
            involve the use of random elements and the permutation
            representation (or even the degree of the representation) is
            not guaranteed to be the same for different calls of
            SmallerDegreePermutationRepresentation."

            To obtain a reproducible result the optional argument ``seed``
            may be used as in the example above.

        TESTS::

            sage: A= matrix(QQ, 2, [0, 1, 1, 0])
            sage: B= matrix(QQ, 2, [1, 0, 0, 1])
            sage: a, b= MatrixGroup([A, B]).as_permutation_group().gens()
            sage: a.order(), b.order()
            (2, 1)

        The above example in GL(12,Z), reduced modulo 7::

            sage: MS = MatrixSpace(GF(7), 12, 12)
            sage: G = MatrixGroup([MS(g) for g in GG.GeneratorsOfGroup()])
            sage: G.cardinality()
            21499084800
            sage: P = G.as_permutation_group()
            sage: P.cardinality()
            21499084800

        Check that large degree is still working::

            sage: Sp(6,3).as_permutation_group().cardinality()
            9170703360

        Check that :trac:`25706` still works after :trac:`26903`::

            sage: MG = GU(3,2).as_matrix_group()
            sage: PG = MG.as_permutation_group()
            sage: mg = MG.an_element()
            sage: PG(mg).order() # particular element depends on the set of GAP packages installed
            6
        """
        # Note that the output of IsomorphismPermGroup() depends on
        # memory locations and will change if you change the order of
        # doctests and/or architecture
        from sage.groups.perm_gps.permgroup import PermutationGroup
        if not self.is_finite():
            raise NotImplementedError("Group must be finite.")
        if seed is not None:
            from sage.libs.gap.libgap import libgap
            libgap.set_seed(ZZ(seed))
        iso = self._libgap_().IsomorphismPermGroup()
        if algorithm == "smaller":
            iso = iso.Image().SmallerDegreePermutationRepresentation()
        return PermutationGroup(iso.Image().GeneratorsOfGroup().sage(),
                                canonicalize=False)

    def module_composition_factors(self, algorithm=None):
        r"""
        Return a list of triples consisting of [base field, dimension,
        irreducibility], for each of the Meataxe composition factors
        modules. The ``algorithm="verbose"`` option returns more information,
        but in Meataxe notation.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,4,4)
            sage: M = MS(0)
            sage: M[0,1]=1;M[1,2]=1;M[2,3]=1;M[3,0]=1
            sage: G = MatrixGroup([M])
            sage: G.module_composition_factors()
            [(Finite Field of size 3, 1, True),
             (Finite Field of size 3, 1, True),
             (Finite Field of size 3, 2, True)]
            sage: F = GF(7); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[0,1],[-1,0]]),MS([[1,1],[2,3]])]
            sage: G = MatrixGroup(gens)
            sage: G.module_composition_factors()
            [(Finite Field of size 7, 2, True)]

        Type ``G.module_composition_factors(algorithm='verbose')`` to get a
        more verbose version.

        For more on MeatAxe notation, see
        https://www.gap-system.org/Manuals/doc/ref/chap69.html
        """
        from sage.libs.gap.libgap import libgap
        F = self.base_ring()
        if not F.is_finite():
            raise NotImplementedError("Base ring must be finite.")
        n = self.degree()
        MS = MatrixSpace(F, n, n)
        mats = [MS(g.matrix()) for g in self.gens()]
        # initializing list of mats by which the gens act on self
        mats_gap = libgap(mats)
        M = mats_gap.GModuleByMats(F)
        compo = libgap.function_factory('MTX.CompositionFactors')
        MCFs = compo(M)
        if algorithm == "verbose":
            print(str(MCFs) + "\n")
        return sorted((MCF['field'].sage(),
                       MCF['dimension'].sage(),
                       MCF['IsIrreducible'].sage()) for MCF in MCFs)

    def invariant_generators(self):
        r"""
        Return invariant ring generators.

        Computes generators for the polynomial ring
        `F[x_1,\ldots,x_n]^G`, where `G` in `GL(n,F)` is a finite matrix
        group.

        In the "good characteristic" case the polynomials returned
        form a minimal generating set for the algebra of `G`-invariant
        polynomials.  In the "bad" case, the polynomials returned
        are primary and secondary invariants, forming a not
        necessarily minimal generating set for the algebra of
        `G`-invariant polynomials.

        ALGORITHM:

        Wraps Singular's ``invariant_algebra_reynolds`` and ``invariant_ring``
        in ``finvar.lib``.

        EXAMPLES::

            sage: F = GF(7); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[0,1],[-1,0]]),MS([[1,1],[2,3]])]
            sage: G = MatrixGroup(gens)
            sage: G.invariant_generators()
            [x1^7*x2 - x1*x2^7,
             x1^12 - 2*x1^9*x2^3 - x1^6*x2^6 + 2*x1^3*x2^9 + x2^12,
             x1^18 + 2*x1^15*x2^3 + 3*x1^12*x2^6 + 3*x1^6*x2^12 - 2*x1^3*x2^15 + x2^18]

            sage: q = 4; a = 2
            sage: MS = MatrixSpace(QQ, 2, 2)
            sage: gen1 = [[1/a,(q-1)/a],[1/a, -1/a]]; gen2 = [[1,0],[0,-1]]; gen3 = [[-1,0],[0,1]]
            sage: G = MatrixGroup([MS(gen1),MS(gen2),MS(gen3)])
            sage: G.cardinality()
            12
            sage: G.invariant_generators()
            [x1^2 + 3*x2^2, x1^6 + 15*x1^4*x2^2 + 15*x1^2*x2^4 + 33*x2^6]

            sage: F = CyclotomicField(8)
            sage: z = F.gen()
            sage: a = z+1/z
            sage: b = z^2
            sage: MS = MatrixSpace(F,2,2)
            sage: g1 = MS([[1/a, 1/a], [1/a, -1/a]])
            sage: g2 = MS([[-b, 0], [0, b]])
            sage: G=MatrixGroup([g1,g2])
            sage: G.invariant_generators()
            [x1^4 + 2*x1^2*x2^2 + x2^4,
             x1^5*x2 - x1*x2^5,
             x1^8 + 28/9*x1^6*x2^2 + 70/9*x1^4*x2^4 + 28/9*x1^2*x2^6 + x2^8]

        AUTHORS:

        - David Joyner, Simon King and Martin Albrecht.

        REFERENCES:

        - Singular reference manual

        - [Stu1993]_

        - S. King, "Minimal Generating Sets of non-modular invariant
          rings of finite groups", :arxiv:`math/0703035`.
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.interfaces.singular import singular
        gens = self.gens()
        singular.LIB("finvar.lib")
        n = self.degree()  # len((gens[0].matrix()).rows())
        F = self.base_ring()
        q = F.characteristic()
        # test if the field is admissible
        if F.gen() == 1:  # we got the rationals or GF(prime)
            FieldStr = str(F.characteristic())
        elif hasattr(F,'polynomial'):  # we got an algebraic extension
            if len(F.gens()) > 1:
                raise NotImplementedError("can only deal with finite fields and (simple algebraic extensions of) the rationals")
            FieldStr = '(%d,%s)' % (F.characteristic(), str(F.gen()))
        else:  # we have a transcendental extension
            FieldStr = '(%d,%s)' % (F.characteristic(),
                                    ','.join(str(p) for p in F.gens()))

        # Setting Singular's variable names
        # We need to make sure that field generator and variables get different names.
        if str(F.gen())[0] == 'x':
            VarStr = 'y'
        else:
            VarStr = 'x'
        VarNames = '(' + ','.join((VarStr+str(i) for i in range(1, n+1)))+')'
        # The function call and affectation below have side-effects. Do not remove!
        # (even if pyflakes say so)
        R = singular.ring(FieldStr, VarNames, 'dp')
        if hasattr(F, 'polynomial') and F.gen() != 1:
            # we have to define minpoly
            singular.eval('minpoly = '+str(F.polynomial()).replace('x',str(F.gen())))
        A = [singular.matrix(n,n,str((x.matrix()).list())) for x in gens]
        Lgens = ','.join((x.name() for x in A))
        PR = PolynomialRing(F, n, [VarStr+str(i) for i in range(1,n+1)])

        if q == 0 or (q > 0 and self.cardinality() % q):
            from sage.all import Matrix
            try:
                elements = [g.matrix() for g in self.list()]
            except (TypeError, ValueError):
                elements
            if elements is not None:
                ReyName = 't'+singular._next_var_name()
                singular.eval('matrix %s[%d][%d]' % (ReyName,
                                                     self.cardinality(), n))
                for i in range(1,self.cardinality()+1):
                    M = Matrix(F, elements[i-1])
                    D = [{} for foobar in range(self.degree())]
                    for x,y in M.dict().items():
                        D[x[0]][x[1]] = y
                    for row in range(self.degree()):
                        for t in D[row].items():
                            singular.eval('%s[%d,%d]=%s[%d,%d]+(%s)*var(%d)'
                                          % (ReyName,i,row+1,ReyName,i,row+1, repr(t[1]),t[0]+1))
                IRName = 't'+singular._next_var_name()
                singular.eval('matrix %s = invariant_algebra_reynolds(%s)' % (IRName,ReyName))
            else:
                ReyName = 't'+singular._next_var_name()
                singular.eval('list %s=group_reynolds((%s))' % (ReyName, Lgens))
                IRName = 't'+singular._next_var_name()
                singular.eval('matrix %s = invariant_algebra_reynolds(%s[1])' % (IRName, ReyName))

            OUT = [singular.eval(IRName+'[1,%d]' % (j))
                   for j in range(1, 1+int(singular('ncols('+IRName+')')))]
            return [PR(gen) for gen in OUT]
        if self.cardinality() % q == 0:
            PName = 't' + singular._next_var_name()
            SName = 't' + singular._next_var_name()
            singular.eval('matrix %s,%s=invariant_ring(%s)' % (PName, SName, Lgens))
            OUT = [singular.eval(PName+'[1,%d]' % (j))
                   for j in range(1,1+singular('ncols('+PName+')'))]
            OUT += [singular.eval(SName+'[1,%d]' % (j))
                    for j in range(2,1+singular('ncols('+SName+')'))]
            return [PR(gen) for gen in OUT]

    def molien_series(self, chi=None, return_series=True, prec=20, variable='t'):
        r"""
        Compute the Molien series of this finite group with respect to the
        character ``chi``. It can be returned either as a rational function
        in one variable or a power series in one variable. The base field
        must be a finite field, the rationals, or a cyclotomic field.

        Note that the base field characteristic cannot divide the group
        order (i.e., the non-modular case).

        ALGORITHM:

        For a finite group `G` in characteristic zero we construct the Molien series as

        .. MATH::

            \frac{1}{|G|}\sum_{g \in G} \frac{\chi(g)}{\text{det}(I-tg)},

        where `I` is the identity matrix and `t` an indeterminate.

        For characteristic `p` not dividing the order of `G`, let `k` be the base field
        and `N` the order of `G`. Define `\lambda` as a primitive `N`-th root of unity over `k`
        and `\omega` as a primitive `N`-th root of unity over `\QQ`. For each `g \in G`
        define `k_i(g)` to be the positive integer such that
        `e_i = \lambda^{k_i(g)}` for each eigenvalue `e_i` of `g`. Then the Molien series
        is computed as

        .. MATH::

            \frac{1}{|G|}\sum_{g \in G} \frac{\chi(g)}{\prod_{i=1}^n(1 - t\omega^{k_i(g)})},

        where `t` is an indeterminant. [Dec1998]_

        INPUT:

        - ``chi`` -- (default: trivial character) a linear group character of this group

        - ``return_series`` -- boolean (default: ``True``) if ``True``, then returns
          the Molien series as a power series, ``False`` as a rational function

        - ``prec`` -- integer (default: 20); power series default precision

        - ``variable`` -- string (default: ``'t'``); Variable name for the Molien series

        OUTPUT: single variable rational function or power series with integer coefficients

        EXAMPLES::

            sage: MatrixGroup(matrix(QQ,2,2,[1,1,0,1])).molien_series()
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for finite groups
            sage: MatrixGroup(matrix(GF(3),2,2,[1,1,0,1])).molien_series()
            Traceback (most recent call last):
            ...
            NotImplementedError: characteristic cannot divide group order

        Tetrahedral Group::

            sage: K.<i> = CyclotomicField(4)
            sage: Tetra =  MatrixGroup([(-1+i)/2,(-1+i)/2, (1+i)/2,(-1-i)/2], [0,i, -i,0])
            sage: Tetra.molien_series(prec=30)
            1 + t^8 + 2*t^12 + t^16 + 2*t^20 + 3*t^24 + 2*t^28 + O(t^30)
            sage: mol = Tetra.molien_series(return_series=False); mol
            (t^8 - t^4 + 1)/(t^16 - t^12 - t^4 + 1)
            sage: mol.parent()
            Fraction Field of Univariate Polynomial Ring in t over Integer Ring
            sage: chi = Tetra.character(Tetra.character_table()[1])
            sage: Tetra.molien_series(chi, prec=30, variable='u')
            u^6 + u^14 + 2*u^18 + u^22 + 2*u^26 + 3*u^30 + 2*u^34 + O(u^36)
            sage: chi = Tetra.character(Tetra.character_table()[2])
            sage: Tetra.molien_series(chi)
            t^10 + t^14 + t^18 + 2*t^22 + 2*t^26 + O(t^30)

        ::

            sage: S3 = MatrixGroup(SymmetricGroup(3))
            sage: mol = S3.molien_series(prec=10); mol
            1 + t + 2*t^2 + 3*t^3 + 4*t^4 + 5*t^5 + 7*t^6 + 8*t^7 + 10*t^8 + 12*t^9 + O(t^10)
            sage: mol.parent()
            Power Series Ring in t over Integer Ring

        Octahedral Group::

            sage: K.<v> = CyclotomicField(8)
            sage: a = v-v^3 #sqrt(2)
            sage: i = v^2
            sage: Octa = MatrixGroup([(-1+i)/2,(-1+i)/2, (1+i)/2,(-1-i)/2], [(1+i)/a,0, 0,(1-i)/a])
            sage: Octa.molien_series(prec=30)
            1 + t^8 + t^12 + t^16 + t^18 + t^20 + 2*t^24 + t^26 + t^28 + O(t^30)

        Icosahedral Group::

            sage: K.<v> = CyclotomicField(10)
            sage: z5 = v^2
            sage: i = z5^5
            sage: a = 2*z5^3 + 2*z5^2 + 1 #sqrt(5)
            sage: Ico = MatrixGroup([[z5^3,0, 0,z5^2], [0,1, -1,0], [(z5^4-z5)/a, (z5^2-z5^3)/a, (z5^2-z5^3)/a, -(z5^4-z5)/a]])
            sage: Ico.molien_series(prec=40)
            1 + t^12 + t^20 + t^24 + t^30 + t^32 + t^36 + O(t^40)

        ::

            sage: G = MatrixGroup(CyclicPermutationGroup(3))
            sage: chi = G.character(G.character_table()[1])
            sage: G.molien_series(chi, prec=10)
            t + 2*t^2 + 3*t^3 + 5*t^4 + 7*t^5 + 9*t^6 + 12*t^7 + 15*t^8 + 18*t^9 + 22*t^10 + O(t^11)

        ::

            sage: K = GF(5)
            sage: S = MatrixGroup(SymmetricGroup(4))
            sage: G = MatrixGroup([matrix(K,4,4,[K(y) for u in m.list() for y in u])for m in S.gens()])
            sage: G.molien_series(return_series=False)
            1/(t^10 - t^9 - t^8 + 2*t^5 - t^2 - t + 1)

        ::

            sage: i = GF(7)(3)
            sage: G = MatrixGroup([[i^3,0,0,-i^3],[i^2,0,0,-i^2]])
            sage: chi = G.character(G.character_table()[4])
            sage: G.molien_series(chi)
            3*t^5 + 6*t^11 + 9*t^17 + 12*t^23 + O(t^25)
        """
        if not self.is_finite():
            raise NotImplementedError("only implemented for finite groups")
        if chi is None:
            chi = self.trivial_character()
        M = self.matrix_space()
        R = FractionField(self.base_ring())
        N = self.order()
        if R.characteristic() == 0:
            P = PolynomialRing(R, variable)
            t = P.gen()
            # it is possible the character is over a larger cyclotomic field
            K = chi.values()[0].parent()
            if K.degree() != 1:
                if R.degree() != 1:
                    L = K.composite_fields(R)[0]
                else:
                    L = K
            else:
                L = R
            mol = P(0)
            for g in self:
                mol += L(chi(g)) / (M.identity_matrix()-t*g.matrix()).det().change_ring(L)
        elif R.characteristic().divides(N):
            raise NotImplementedError("characteristic cannot divide group order")
        else:  # char p>0
            # find primitive Nth roots of unity over base ring and QQ
            F = cyclotomic_polynomial(N).change_ring(R)
            w = F.roots(ring=R.algebraic_closure(), multiplicities=False)[0]
            # don't need to extend further in this case since the order of
            # the roots of unity in the character divide the order of the group
            L = CyclotomicField(N, 'v')
            v = L.gen()
            # construct Molien series
            P = PolynomialRing(L, variable)
            t = P.gen()
            mol = P(0)
            for g in self:
                # construct Phi
                phi = L(chi(g))
                for e in g.matrix().eigenvalues():
                    # find power such that w**n  = e
                    n = 1
                    while w**n != e and n < N+1:
                        n += 1
                    # raise v to that power
                    phi *= (1-t*v**n)
                mol += P(1)/phi
        # We know the coefficients will be integers
        mol = mol.numerator().change_ring(ZZ) / mol.denominator().change_ring(ZZ)
        # divide by group order
        mol /= N
        if return_series:
            PS = PowerSeriesRing(ZZ, variable, default_prec=prec)
            return PS(mol)
        return mol

    def reynolds_operator(self, poly, chi=None):
        r"""
        Compute the Reynolds operator of this finite group `G`.

        This is the projection from a polynomial ring to the ring of
        relative invariants [Stu1993]_. If possible, the invariant is
        returned defined over the base field of the given polynomial
        ``poly``, otherwise, it is returned over the compositum of the
        fields involved in the computation.
        Only implemented for absolute fields.

        ALGORITHM:

        Let `K[x]` be a polynomial ring and `\chi` a linear character for `G`. Let

        .. MATH:

            K[x]^G_{\chi} = \{f \in K[x] | \pi f = \chi(\pi) f \forall \pi\in G\}

        be the ring of invariants of `G` relative to `\chi`. Then the Reynold's operator
        is a map `R` from `K[x]` into `K[x]^G_{\chi}` defined by

        .. MATH:

            f \mapsto \frac{1}{|G|} \sum_{ \pi \in G} \chi(\pi) f.

        INPUT:

        - ``poly`` -- a polynomial

        - ``chi`` -- (default: trivial character) a linear group character of this group

        OUTPUT: an invariant polynomial relative to `\chi`

        AUTHORS:

        Rebecca Lauren Miller and Ben Hutz

        EXAMPLES::

            sage: S3 = MatrixGroup(SymmetricGroup(3))
            sage: R.<x,y,z> = QQ[]
            sage: f = x*y*z^3
            sage: S3.reynolds_operator(f)
            1/3*x^3*y*z + 1/3*x*y^3*z + 1/3*x*y*z^3

        ::

            sage: G = MatrixGroup(CyclicPermutationGroup(4))
            sage: chi = G.character(G.character_table()[3])
            sage: K.<v> = CyclotomicField(4)
            sage: R.<x,y,z,w> = K[]
            sage: G.reynolds_operator(x, chi)
            1/4*x + (1/4*v)*y - 1/4*z + (-1/4*v)*w
            sage: chi = G.character(G.character_table()[2])
            sage: R.<x,y,z,w> = QQ[]
            sage: G.reynolds_operator(x*y, chi)
            1/4*x*y + (-1/4*zeta4)*y*z + (1/4*zeta4)*x*w - 1/4*z*w

        ::

            sage: K.<i> = CyclotomicField(4)
            sage: G =  MatrixGroup(CyclicPermutationGroup(3))
            sage: chi = G.character(G.character_table()[1])
            sage: R.<x,y,z> = K[]
            sage: G.reynolds_operator(x*y^5, chi)
            1/3*x*y^5 + (-2/3*izeta3^3 - izeta3^2 - 8/3*izeta3 - 4/3)*x^5*z + (2/3*izeta3^3 + izeta3^2 + 8/3*izeta3 + 1)*y*z^5
            sage: R.<x,y,z> = QQbar[]
            sage: G.reynolds_operator(x*y^5, chi)
             1/3*x*y^5 + (-0.1666666666666667? + 0.2886751345948129?*I)*x^5*z + (-0.1666666666666667? - 0.2886751345948129?*I)*y*z^5

        ::

            sage: K.<i> = CyclotomicField(4)
            sage: Tetra =  MatrixGroup([(-1+i)/2,(-1+i)/2, (1+i)/2,(-1-i)/2], [0,i, -i,0])
            sage: chi = Tetra.character(Tetra.character_table()[4])
            sage: L.<v> = QuadraticField(-3)
            sage: R.<x,y> = L[]
            sage: Tetra.reynolds_operator(x^4)
            0
            sage: Tetra.reynolds_operator(x^4, chi)
            1/4*x^4 + (1/2*v)*x^2*y^2 + 1/4*y^4
            sage: R.<x>=L[]
            sage: LL.<w> = L.extension(x^2+v)
            sage: R.<x,y> = LL[]
            sage: Tetra.reynolds_operator(x^4, chi)
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for absolute fields

        ::

            sage: G =  MatrixGroup(DihedralGroup(4))
            sage: chi = G.character(G.character_table()[1])
            sage: R.<x,y> = QQ[]
            sage: f = x^4
            sage: G.reynolds_operator(f, chi)
            Traceback (most recent call last):
            ...
            TypeError: number of variables in polynomial must match size of matrices
            sage: R.<x,y,z,w> = QQ[]
            sage: f = x^3*y
            sage: G.reynolds_operator(f, chi)
            1/8*x^3*y - 1/8*x*y^3 + 1/8*y^3*z - 1/8*y*z^3 - 1/8*x^3*w + 1/8*z^3*w +
            1/8*x*w^3 - 1/8*z*w^3

        Characteristic p>0 examples::

            sage: G = MatrixGroup([[0,1,1,0]])
            sage: R.<w,x> = GF(2)[]
            sage: G.reynolds_operator(x)
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented when characteristic divides group order

        ::

            sage: i = GF(7)(3)
            sage: G = MatrixGroup([[i^3,0,0,-i^3],[i^2,0,0,-i^2]])
            sage: chi = G.character(G.character_table()[4])
            sage: R.<w,x> = GF(7)[]
            sage: f = w^5*x + x^6
            sage: G.reynolds_operator(f, chi)
            Traceback (most recent call last):
            ...
            NotImplementedError: nontrivial characters not implemented for characteristic > 0
            sage: G.reynolds_operator(f)
            x^6

        ::

            sage: K = GF(3^2,'t')
            sage: G = MatrixGroup([matrix(K,2,2, [0,K.gen(),1,0])])
            sage: R.<x,y> = GF(3)[]
            sage: G.reynolds_operator(x^8)
            -x^8 - y^8

        ::

            sage: K = GF(3^2,'t')
            sage: G = MatrixGroup([matrix(GF(3),2,2, [0,1,1,0])])
            sage: R.<x,y> = K[]
            sage: f = -K.gen()*x
            sage: G.reynolds_operator(f)
            t*x + t*y
        """
        if poly.parent().ngens() != self.degree():
            raise TypeError("number of variables in polynomial must match size of matrices")
        R = FractionField(poly.base_ring())
        C = FractionField(self.base_ring())
        if chi is None:  # then this is the trivial character
            if R.characteristic() == 0:
                # non-modular case
                if C == QQbar or R == QQbar:
                    L = QQbar
                elif not C.is_absolute() or not R.is_absolute():
                    raise NotImplementedError("only implemented for absolute fields")
                else:  # create the compositum
                    if C.absolute_degree() == 1:
                        L = R
                    elif R.absolute_degree() == 1:
                        L = C
                    else:
                        L = C.composite_fields(R)[0]
            elif not R.characteristic().divides(self.order()):
                if R.characteristic() != C.characteristic():
                    raise ValueError("base fields must have same characteristic")
                else:
                    if R.degree() >= C.degree():
                        L = R
                    else:
                        L = C
            else:
                raise NotImplementedError("not implemented when characteristic divides group order")
            poly = poly.change_ring(L)
            poly_gens = vector(poly.parent().gens())
            F = L.zero()
            for g in self:
                F += poly(*g.matrix()*vector(poly.parent().gens()))
            F /= self.order()
            return F
        # non-trivial character case
        K = chi.values()[0].parent()
        if R.characteristic() == 0:
            # extend base_ring to compositum
            if C == QQbar or K == QQbar or R == QQbar:
                L = QQbar
            elif not C.is_absolute() or not K.is_absolute() or not R.is_absolute():
                raise NotImplementedError("only implemented for absolute fields")
            else:
                fields = []
                for M in [R,K,C]:
                    if M.absolute_degree() != 1:
                        fields.append(M)
                l = len(fields)
                if l == 0:
                    # all are QQ
                    L = R
                elif l == 1:
                    # only one is an extension
                    L = fields[0]
                elif l == 2:
                    # only two are extensions
                    L = fields[0].composite_fields(fields[1])[0]
                else:
                    # all three are extensions
                    L1 = fields[0].composite_fields(fields[1])[0]
                    L = L1.composite_fields(fields[2])[0]
        else:
            raise NotImplementedError("nontrivial characters not implemented for characteristic > 0")
        poly = poly.change_ring(L)
        poly_gens = vector(poly.parent().gens())
        F = L.zero()
        for g in self:
            F += L(chi(g)) * poly(*g.matrix().change_ring(L)*poly_gens)
        F /= self.order()
        try:  # attempt to move F to base_ring of polynomial
            F = F.change_ring(R)
        except (TypeError, ValueError):
            pass
        return F

    def invariants_of_degree(self, deg, chi=None, R=None):
        r"""
        Return the (relative) invariants of given degree for this group.

        For this group, compute the invariants of degree ``deg``
        with respect to the group character ``chi``. The method
        is to project each possible monomial of degree ``deg`` via
        the Reynolds operator. Note that if the polynomial ring ``R``
        is specified it's base ring may be extended if the resulting
        invariant is defined over a bigger field.

        INPUT:

        - ``degree`` -- a positive integer

        - ``chi`` -- (default: trivial character) a linear group character of this group

        - ``R`` -- (optional) a polynomial ring

        OUTPUT: list of polynomials

        EXAMPLES::

            sage: Gr = MatrixGroup(SymmetricGroup(2))
            sage: sorted(Gr.invariants_of_degree(3))
            [x0^2*x1 + x0*x1^2, x0^3 + x1^3]
            sage: R.<x,y> = QQ[]
            sage: sorted(Gr.invariants_of_degree(4, R=R))
            [x^2*y^2, x^3*y + x*y^3, x^4 + y^4]

        ::

            sage: R.<x,y,z> = QQ[]
            sage: Gr = MatrixGroup(DihedralGroup(3))
            sage: ct = Gr.character_table()
            sage: chi = Gr.character(ct[0])
            sage: all(f(*(g.matrix()*vector(R.gens()))) == chi(g)*f
            ....: for f in Gr.invariants_of_degree(3, R=R, chi=chi) for g in Gr)
            True

        ::

            sage: i = GF(7)(3)
            sage: G = MatrixGroup([[i^3,0,0,-i^3],[i^2,0,0,-i^2]])
            sage: G.invariants_of_degree(25)
            []

        ::

            sage: G = MatrixGroup(SymmetricGroup(5))
            sage: R = QQ['x,y']
            sage: G.invariants_of_degree(3, R=R)
            Traceback (most recent call last):
            ...
            TypeError: number of variables in polynomial ring must match size of matrices

        ::

            sage: K.<i> = CyclotomicField(4)
            sage: G =  MatrixGroup(CyclicPermutationGroup(3))
            sage: chi = G.character(G.character_table()[1])
            sage: R.<x,y,z> = K[]
            sage: sorted(G.invariants_of_degree(2, R=R, chi=chi))
            [x*y + (-2*izeta3^3 - 3*izeta3^2 - 8*izeta3 - 4)*x*z + (2*izeta3^3 + 3*izeta3^2 + 8*izeta3 + 3)*y*z,
             x^2 + (2*izeta3^3 + 3*izeta3^2 + 8*izeta3 + 3)*y^2 + (-2*izeta3^3 - 3*izeta3^2 - 8*izeta3 - 4)*z^2]

        ::

            sage: S3 = MatrixGroup(SymmetricGroup(3))
            sage: chi = S3.character(S3.character_table()[0])
            sage: sorted(S3.invariants_of_degree(5, chi=chi))
            [x0^3*x1^2 - x0^2*x1^3 - x0^3*x2^2 + x1^3*x2^2 + x0^2*x2^3 - x1^2*x2^3,
            x0^4*x1 - x0*x1^4 - x0^4*x2 + x1^4*x2 + x0*x2^4 - x1*x2^4]
        """
        D = self.degree()
        deg = int(deg)
        if deg <= 0:
            raise ValueError("degree must be a positive integer")
        if R is None:
            R = PolynomialRing(self.base_ring(), 'x', D)
        elif R.ngens() != D:
            raise TypeError("number of variables in polynomial ring must match size of matrices")

        ms = self.molien_series(prec=deg+1,chi=chi)
        if ms[deg].is_zero():
            return []
        inv = set()
        for e in IntegerVectors(deg, D):
            F = self.reynolds_operator(R.monomial(*e), chi=chi)
            if not F.is_zero():
                F = F / F.lc()
                inv.add(F)
                if len(inv) == ms[deg]:
                    break
        return list(inv)
