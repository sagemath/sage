r"""
Lattices of sashes

These lattices were introduced by S. Law in [Law2014]_. They are
lattice quotients of the weak order on the symmetric groups. This
implies that they are congruence-uniform.

There is a lattice of sashes `\Sigma_n` for every integer `n \geq 1`.
The underlying set of `\Sigma_n` is the set of words of length `n`
in the letters ``□``, ``▨`` and ``■■`` of respective lengths `1,1` and `2`.

The cardinalities of the lattices `\Sigma_n` are therefore given by
the Pell numbers (:oeis:`A000129`).

The implementation describes sashes as tuples of strings.

This file also implements the polytopes known as pellytopes, and their fans.
The skeletons of pellytopes are isomorphic to the undirected Hasse
diagrams of the lattices of sashes.

REFERENCES:

- [Law2014]_
"""

from itertools import pairwise
from collections.abc import Iterator

from sage.categories.finite_lattice_posets import FiniteLatticePosets
from sage.combinat.posets.lattices import LatticePoset
from sage.geometry.cone import Cone
from sage.geometry.fan import Fan
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.misc.cachefunc import cached_function
from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ

B, N, BB = "□", "▨", "■■"


def sashes(n: int) -> Iterator[tuple[str, ...]]:
    """
    Iterate over the sashes of length ``n``.

    INPUT:

    - ``n`` -- nonnegative integer

    EXAMPLES::

        sage: from sage.combinat.posets.sashes import sashes
        sage: [''.join(s) for s in sashes(2)]
        ['□□', '□▨', '▨□', '▨▨', '■■']

    TESTS::

        sage: [''.join(s) for s in sashes(0)]
        ['']
        sage: [''.join(s) for s in sashes(1)]
        ['□', '▨']
        sage: list(sashes(-1))
        Traceback (most recent call last):
        ...
        ValueError: n must be nonnegative
    """
    if n < 0:
        raise ValueError("n must be nonnegative")
    if n == 0:
        yield ()
        return
    if n == 1:
        yield (B,)
        yield (N,)
        return
    for s in sashes(n - 1):
        yield s + (B,)
        yield s + (N,)
    for s in sashes(n - 2):
        yield s + (BB,)


def cover_relations(s: tuple[str, ...]) -> Iterator[tuple[str, ...]]:
    """
    Iterate over the cover relations of the given sash.

    INPUT:

    - ``s`` -- sash, as a tuple of strings

    EXAMPLES::

        sage: from sage.combinat.posets.sashes import cover_relations
        sage: s = ('□', '▨', '▨')
        sage: [''.join(t) for t in cover_relations(s)]
        ['□□▨', '□▨□']
        sage: s = ('□', '■■', '▨')
        sage: [''.join(t) for t in cover_relations(s)]
        ['□□□▨', '□■■□']
    """
    for i, letter in enumerate(s):
        if letter == BB:
            yield s[:i] + (B, B) + s[i + 1:]
    for i, (l1, l2) in enumerate(pairwise(s)):
        if l1 == N:
            if l2 == B:
                yield s[:i] + (BB,) + s[i + 2:]
            else:
                yield s[:i] + (B,) + s[i + 1:]
    if s[-1] == N:
        yield s[:-1] + (B,)


def lattice_of_sashes(n: int) -> LatticePoset:
    """
    Return the lattice of sashes of length ``n``.

    INPUT:

    - ``n`` -- positive integer

    EXAMPLES::

        sage: L = posets.Sashes(4); L
        Finite lattice containing 29 elements
        sage: L.hasse_diagram().to_undirected().is_regular(4)
        True

    TESTS::

        sage: posets.Sashes(0)
        Traceback (most recent call last):
        ...
        ValueError: n must be positive
    """
    if n <= 0:
        raise ValueError("n must be positive")
    cat = FiniteLatticePosets().CongruenceUniform()
    return LatticePoset({s: list(cover_relations(s)) for s in sashes(n)},
                        cover_relations=True, check=False,
                        category=cat)


@cached_function
def pellytope_fan(n: int) -> Fan:
    """
    Return the fan of the pellytope of dimension ``n``.

    This is defined by induction.

    INPUT:

    - ``n`` -- integer

    EXAMPLES::

        sage: from sage.combinat.posets.sashes import pellytope_fan
        sage: pellytope_fan(3).f_vector()
        (1, 8, 18, 12)

    TESTS::

        sage: pellytope_fan(1)
        Rational polyhedral fan in 1-d lattice N
        sage: pellytope_fan(0)
        Traceback (most recent call last):
        ...
        ValueError: n must be positive

    REFERENCES:

    - [BTTM2024]_
    """
    if n <= 0:
        raise ValueError("n must be positive")
    dim_one = Fan([Cone([[-1]]), Cone([[1]])])
    if n == 1:
        return dim_one
    G = pellytope_fan(n - 1).cartesian_product(dim_one)
    v = vector([0] * (n - 2) + [-1, 1])
    return G.subdivide(new_rays=[v])


def pellytope(n: int) -> Polyhedron:
    r"""
    Return the pellytope of dimension ``n``.

    This is defined as a Minkowski sum of the unit volume `n`
    dimensional hypercube in the positive orthant and the
    triangles with vertices `(0, e_i, e_i + e_{i+1})` for
    all `1 \leq i < n`.

    INPUT:

    - ``n`` -- integer

    EXAMPLES::

        sage: P3 = polytopes.pellytope(3); P3
        A 3-dimensional polyhedron in ZZ^3 defined as
        the convex hull of 12 vertices
        sage: P3.f_vector()
        (1, 12, 18, 8, 1)

    TESTS::

        sage: polytopes.pellytope(0)
        Traceback (most recent call last):
        ...
        ValueError: n must be positive

        sage: G = posets.Sashes(3).hasse_diagram().to_undirected()
        sage: polytopes.pellytope(3).graph().is_isomorphic(G)
        True

    REFERENCES:

    - [BTTM2024]_
    """
    from sage.geometry.polyhedron.library import Polytopes
    if n <= 0:
        raise ValueError("n must be positive")
    M = FreeModule(ZZ, n)
    v = M.basis()
    zero = M.zero()

    resu = Polytopes().hypercube(n, intervals='zero_one')

    return resu + sum(Polyhedron(vertices=[zero, v[i], v[i] + v[i + 1]])
                      for i in range(n - 1))
