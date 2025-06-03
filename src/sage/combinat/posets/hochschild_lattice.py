r"""
Hochschild lattices

Hochschild lattices were defined in [Cha2020]_ as posets and
shown to be congruence uniform lattices in [Com2021]_.

The name comes from the Hochschild polytopes introduced by Sanebdlidze
in [San2009]_ for the study of free loop spaces in algebraic topology.
The Hasse diagram of the Hochschild lattice is an orientation of the
1-skeleton of the Hochschild polytope.

For `n \geq 1`, the cardinality of the Hochschild lattice `H(n)` is
`2^{n - 2} \times (n + 3)`, starting with
`2, 5, 12, 28, 64, 144, 320, 704, 1536, 3328, 7168, 15360, \ldots`.

The underlying set of `H(n)` consists of some words in the alphabet `(0,1,2)`, whose precise description can be found in [Com2021]_.
"""
from typing import Iterator

from sage.combinat.posets.lattices import LatticePoset
from sage.misc.lazy_import import lazy_import
from sage.modules.free_module_element import vector
from sage.topology.simplicial_complex import SimplicialComplex

lazy_import("sage.geometry.cone", "Cone")
lazy_import("sage.geometry.fan", "Fan")


def hochschild_lattice(n) -> LatticePoset:
    r"""
    Return the Hochschild lattice `H(n)`.

    INPUT:

    - `n \geq 1` -- an integer

    The cardinality of `H(n)` is `2^{n - 2} \times (n + 3)`.

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

    verts = list(sommets(n))

    def compare(a, b) -> bool:
        return all(ai <= bi for ai, bi in zip(a, b))

    return LatticePoset([verts, compare])


def hochschild_fan(n):
    """
    Return Saneblidze's fan for the Hochschild polytope.

    The dual polytope is obtained from a standard simplex
    by a sequence of truncations.

    EXAMPLES::

        sage: from sage.combinat.posets.hochschild_lattice import hochschild_fan
        sage: F = hochschild_fan(4); F
        Rational polyhedral fan in 4-d lattice N
        sage: F.f_vector()
        (1, 11, 39, 56, 28)
    """
    rays = [vector([1] * n)]
    rays.extend(vector([0 if j != i else -1 for j in range(n)])
                for i in range(n))

    cones = [Cone([rays[j] for j in range(n + 1) if j != i])
             for i in range(n + 1)]

    # standard fan of projective space Pn
    F = Fan(cones)

    # double sequence of blowups
    subdiv = [sum(r for r in rays[k + 1:]) for k in range(n - 1)]
    subdiv.extend(sum(r for r in rays[:n - k]) for k in range(n - 1))

    return F.subdivide(subdiv)


def hochschild_simplicial_complex(n) -> SimplicialComplex:
    """
    Return a simplicial complex related to the Hochschild lattice `H_{n}`.

    This is a pure spherical simplicial complex, whose flip graph
    is isomorphic to the Hasse diagram of `H_{n}`.

    EXAMPLES::

        sage: from sage.combinat.posets.hochschild_lattice import hochschild_simplicial_complex
        sage: C = hochschild_simplicial_complex(3); C
        Simplicial complex with 8 vertices and 12 facets
        sage: H = C.flip_graph()
        sage: P = posets.HochschildLattice(3)
        sage: H.is_isomorphic(P.hasse_diagram().to_undirected())
        True
    """
    F = hochschild_fan(n)
    return SimplicialComplex([C.rays() for C in F.cones(n)])
