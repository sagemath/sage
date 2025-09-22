r"""
Bubble and Shuffle lattices

Shuffle lattices were defined by Greene in [Gre1988]_.

Bubble lattices were introduced by McConville and Mühle in [MacCM2022]_.

The Bubble lattice `B_{m,n}` and the Shuffle lattice `S_{m,n}` share
the same underlying set, namely all shuffles of two subwords of the
two words `X = (x_1,x_2,\ldots,x_{m})` and `Y = (y_1,y_2,\ldots,y_n)`.

.. NOTE::

    In the implementation here, the underlying set is the set of all shuffles
    of subsets of `\{-m,\ldots,-1\}` with subsets of `\{1,\ldots,n\}`.
"""
from collections.abc import Iterator

from sage.categories.finite_lattice_posets import FiniteLatticePosets
from sage.combinat.posets.lattices import LatticePoset
from sage.combinat.subset import subsets
from sage.combinat.shuffle import ShuffleProduct
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ


def bubble_cardinality(m, n) -> Integer:
    r"""
    Return the cardinality of the Bubble lattice `B_{m,n}`.

    We have

    .. MATH::

        |B_{m,n}| = \sum_{i=0}^m \sum_{j=0}^n \binom{i+j}{j} \binom{m}{i} \binom{n}{j}.

    This is also the cardinality of the Shuffle lattice `S_{m,n}`.

    INPUT:

    - ``m`` -- integer
    - ``n`` -- integer

    EXAMPLES::

        sage: from sage.combinat.posets.bubble_shuffle import bubble_cardinality
        sage: bubble_cardinality(2,1)
        12
    """
    return ZZ.sum(ZZ(i + j).binomial(j) * ZZ(m).binomial(i) * ZZ(n).binomial(j)
                  for i in range(m + 1) for j in range(n + 1))


def bubble_set(m, n) -> Iterator[tuple[int, ...]]:
    r"""
    Return the underlying set of the Bubble lattice `B_{m,n}`.

    This is the set of all shuffles of subsets of `\{-m,\ldots,-1\}`
    with subsets of `\{1,\ldots,n\}`.

    This is also the underlying set of the Shuffle lattice `S_{m,n}`.

    INPUT:

    - ``m`` -- integer
    - ``n`` -- integer

    EXAMPLES::

        sage: from sage.combinat.posets.bubble_shuffle import bubble_set
        sage: list(bubble_set(2,1))
        [(),
         (1,),
         (-1,),
         (-1, 1),
         (1, -1),
         (-2,),
         (-2, 1),
         (1, -2),
         (-1, -2),
         (-1, -2, 1),
         (1, -1, -2),
         (-1, 1, -2)]
    """
    X = [-i for i in range(1, m + 1)]
    Y = list(range(1, n + 1))
    for subX in subsets(X):
        wX = tuple(subX)
        for subY in subsets(Y):
            for shu in ShuffleProduct(wX, tuple(subY)):
                yield tuple(shu)


def bubble_coverings(m, n, mot, transpose=True) -> Iterator[tuple[int, ...]]:
    """
    Return generating relations of the Bubble lattice `B_{m,n}`.

    Note that these relations include the cover relations, but not only them.

    This can also produce covers in the Shuffle lattice `S_{m,n}`.

    INPUT:

    - ``m`` -- integer
    - ``n`` -- integer
    - ``mot`` -- element of `B_{m,n}` as a tuple
    - ``transpose`` -- boolean (default: ``True``) whether to return covers
      in the Bubble lattice or in the Shuffle lattice

    EXAMPLES::

        sage: from sage.combinat.posets.bubble_shuffle import bubble_coverings
        sage: list(bubble_coverings(2, 1, (-2, 1)))
        [(1,), (1, -2)]
        sage: list(bubble_coverings(2, 1, (-2, 1), False))
        [(1,)]
    """
    # removal of one x
    for j, letter in enumerate(mot):
        if letter < 0:
            yield tuple(mot[:j] + mot[j + 1:])

    # insertion of one y
    for j in range(len(mot) + 1):
        debut, fin = tuple(mot[:j]), tuple(mot[j:])
        yavant = max((letter for letter in debut if letter > 0), default=0)
        yapres = min((letter for letter in fin if letter > 0), default=n + 1)
        for y in range(yavant + 1, yapres):
            yield debut + (y,) + fin

    if not transpose:
        return

    # exchange from xy to yx
    for j in range(len(mot) - 1):
        if mot[j] < 0 and mot[j + 1] > 0:
            mot2 = list(mot)
            mot2[j] = mot[j + 1]
            mot2[j + 1] = mot[j]
            yield tuple(mot2)


def BubblePoset(m, n) -> LatticePoset:
    r"""
    Return the Bubble lattice `B_{m,n}`.

    Bubble lattices were introduced by McConville and Mühle in [MacCM2022]_.

    The Bubble lattice `B_{m,n}` and the Shuffle lattice `S_{m,n}` share
    the same underlying set, namely all shuffles of two subwords of the
    two words `X = (x_1,x_2,\ldots,x_{m})` and `Y = (y_1,y_2,\ldots,y_n)`.

    The Bubble poset is an extension of the Shuffle poset, by adding the
    exchange of adjacent letters from `X` and `Y`, from `xy` to `yx`.

    .. SEEALSO::

        :func:`ShufflePoset`, :func:`noncrossing_bipartite_complex`

    EXAMPLES::

        sage: P = posets.BubblePoset(2,1); P
        Finite lattice containing 12 elements
        sage: P.zeta_polynomial()
        1/40*q^5 + 7/24*q^4 + 23/24*q^3 - 7/24*q^2 + 1/60*q
    """
    bubbles = list(bubble_set(m, n))

    dg = DiGraph([(x, y) for x in bubbles for y in bubble_coverings(m, n, x)])
    # here we have more than just the cover relations
    cat = FiniteLatticePosets().CongruenceUniform().Extremal()
    return LatticePoset(dg, category=cat)


def ShufflePoset(m, n) -> LatticePoset:
    r"""
    Return the Shuffle lattice `S_{m,n}`.

    Shuffle lattices were defined by Greene in [Gre1988]_.

    The Bubble lattice `B_{m,n}` and the Shuffle lattice `S_{m,n}` share
    the same underlying set, namely all shuffles of two subwords of the
    two words `X = (x_1,x_2,\ldots,x_{m})` and `Y = (y_1,y_2,\ldots,y_n)`.

    The partial order in the Shuffle poset is defined by either inserting a
    letter from `Y` or deleting a letter from `X`.

    .. SEEALSO:: :func:`BubblePoset`

    EXAMPLES::

        sage: P = posets.ShufflePoset(2,1); P
        Finite lattice containing 12 elements
        sage: P.zeta_polynomial()
        2*q^3 - q^2
    """
    bubbles = list(bubble_set(m, n))

    dg = DiGraph([(x, y) for x in bubbles
                  for y in bubble_coverings(m, n, x, transpose=False)])
    # here we just have the cover relations
    return LatticePoset(dg, cover_relations=True)


def noncrossing_bipartite_complex(m, n):
    """
    Return a simplicial complex related to the Bubble lattice `B_{m,n}`.

    This is a pure spherical simplicial complex, whose flip graph
    is isomorphic to the Hasse diagram of `B_{m,n}`.

    .. SEEALSO:: :func:`BubblePoset`

    EXAMPLES::

        sage: C = simplicial_complexes.NoncrossingBipartiteComplex(2,1)
        sage: H = C.flip_graph()
        sage: P = posets.BubblePoset(2,1)
        sage: H.is_isomorphic(P.hasse_diagram().to_undirected())
        True
    """
    vertices: list[tuple] = [("x", i) for i in range(1, m + 1)]
    vertices.extend(("y", i) for i in range(1, n + 1))
    vertices.extend(("xy", i, j) for i in range(m + 1) for j in range(n + 1)
                    if i or j)

    def compatible(v: tuple, w: tuple) -> bool:
        if v == w:
            return False
        if v[0] != "xy" and w[0] != "xy":
            return True
        if v[0] == "xy" and w[0] == "xy":
            return not ((w[1] < v[1] and w[2] > v[2])
                        or (v[1] < w[1] and v[2] > w[2]))
        if v[0] == "xy":
            if w[0] == "x":
                return v[1] != w[1]
            return v[2] != w[1]
        if v[0] == "x":
            return w[1] != v[1]
        return w[2] != v[1]

    return Graph([vertices, compatible]).clique_complex()
