r"""
Hochschild lattices

Hochschild lattices were defined in [Cha2020]_ as posets and
shown to be congruence uniform lattices in [Com2021]_.

The name comes from the Hochschild polytopes introduced by Sanebdlidze
in [San2009]_ for the study of free loop spaces in algebraic topology.
The Hasse diagram of the Hochschild lattice is an orientation of the
1-skeleton of the Hochschild polytope.

For `n \geq 1`, the cardinality of the Hochschild lattice `H_n` is
`2^{n - 2} \times (n + 3)`, starting with
`2, 5, 12, 28, 64, 144, 320, 704, 1536, 3328, 7168, 15360, \ldots`.

The underlying set of `H_n` consists of some words in the alphabet
`(0,1,2)`, whose precise description can be found in [Com2021]_.
"""
from collections.abc import Iterator

from sage.categories.finite_lattice_posets import FiniteLatticePosets
from sage.combinat.posets.lattices import LatticePoset
from sage.graphs.digraph import DiGraph
from sage.topology.simplicial_complex import SimplicialComplex


def hochschild_lattice(n) -> LatticePoset:
    r"""
    Return the Hochschild lattice `H_n`.

    INPUT:

    - `n \geq 1` -- an integer

    The cardinality of `H_n` is `2^{n - 2} \times (n + 3)`.

    .. SEEALSO:: :func:`hochschild_simplicial_complex`, :func:`hochschild_fan`

    EXAMPLES::

        sage: P = posets.HochschildLattice(5); P
        Finite lattice containing 64 elements
        sage: P.degree_polynomial()
        x^5 + 9*x^4*y + 22*x^3*y^2 + 22*x^2*y^3 + 9*x*y^4 + y^5

    TESTS::

        sage: posets.HochschildLattice(0)
        Traceback (most recent call last):
        ...
        ValueError: this requires n >= 1
    """
    if n <= 0:
        raise ValueError("this requires n >= 1")

    def iterator_True(n) -> Iterator[tuple[int, ...]]:
        """
        Iterator over part of the underlying set.
        """
        if n == 0:
            yield ()
            return
        for w in iterator_True(n - 1):
            yield (0,) + w
            yield (2,) + w

    def iterator_False(n) -> Iterator[tuple[int, ...]]:
        """
        Iterator over rest of the underlying set.
        """
        if n == 0:
            yield ()
            return
        for w in iterator_True(n - 1):
            yield (0,) + w
        for w in iterator_False(n - 1):
            yield (1,) + w
            yield (2,) + w

    def sommets(n) -> Iterator[tuple[int, ...]]:
        """
        Iterator over the full underlying set.
        """
        for w in iterator_True(n - 1):
            yield (0,) + w
        for w in iterator_False(n - 1):
            yield (1,) + w

    verts = frozenset(sommets(n))

    def cover_relations(a):
        for i, ai in enumerate(a):
            if not ai:
                continue
            b = list(a)
            for k in range(ai - 1, -1, -1):
                b[i] = k
                tb = tuple(b)
                if tb in verts:
                    yield tb
                    break

    dg = DiGraph({a: list(cover_relations(a)) for a in verts},
                 format="dict_of_lists")
    return LatticePoset(dg.reverse(), cover_relations=True, check=False,
                        category=FiniteLatticePosets().CongruenceUniform())


def hochschild_fan(n):
    """
    Return Saneblidze's fan for the Hochschild polytope.

    The dual polytope is obtained from a standard simplex
    by a sequence of truncations.

    .. SEEALSO::

        :func:`hochschild_simplicial_complex`, :func:`hochschild_lattice`

    EXAMPLES::

        sage: from sage.combinat.posets.hochschild_lattice import hochschild_fan
        sage: F = hochschild_fan(4); F
        Rational polyhedral fan in 4-d lattice N
        sage: F.f_vector()
        (1, 11, 39, 56, 28)
    """
    from sage.geometry.cone import Cone
    from sage.geometry.fan import Fan
    from sage.modules.free_module_element import vector

    rays = [vector([1] * n)]
    rays.extend(vector([0 if j != i else -1 for j in range(n)])
                for i in range(n))

    cones = [Cone([rays[j] for j in range(n + 1) if j != i])
             for i in range(n + 1)]

    # standard fan of projective space Pn
    F = Fan(cones, check=False)

    # double sequence of blowups
    subdiv = [sum(r for r in rays[k + 1:]) for k in range(n - 1)]
    subdiv.extend(sum(r for r in rays[:n - k]) for k in range(n - 1))

    return F.subdivide(subdiv)


def hochschild_simplicial_complex(n) -> SimplicialComplex:
    """
    Return a simplicial complex related to the Hochschild lattice `H_n`.

    This is a pure spherical simplicial complex, whose flip graph
    is isomorphic to the Hasse diagram of `H_n`.

    .. SEEALSO:: :func:`hochschild_fan`, :func:`hochschild_lattice`

    EXAMPLES::

        sage: C = simplicial_complexes.HochschildSphere(3); C
        Simplicial complex with 8 vertices and 12 facets
        sage: H = C.flip_graph()
        sage: P = posets.HochschildLattice(3)
        sage: H.is_isomorphic(P.hasse_diagram().to_undirected())
        True
    """
    F = hochschild_fan(n)
    return SimplicialComplex([C.rays() for C in F.cones(n)])
