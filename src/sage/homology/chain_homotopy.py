r"""
Chain homotopies and chain contractions

Chain homotopies are standard constructions in homological algebra:
given chain complexes `C` and `D` and chain maps `f, g: C \to D`, say
with differential of degree `-1`, a *chain homotopy* `H` between `f` and
`g` is a collection of maps `H_n: C_n \to D_{n+1}` satisfying

.. MATH::

   \partial_D H + H \partial_C = f - g.

The presence of a chain homotopy defines an equivalence relation
(*chain homotopic*) on chain maps. If `f` and `g` are chain homotopic,
then one can show that `f` and `g` induce the same map on homology.

Chain contractions are not as well known. The papers [MAR2009]_, [RMA2009]_,
and [PR2015]_ provide some references. Given two chain complexes `C` and
`D`, a *chain contraction* is a chain homotopy `H: C \to C` for which
there are chain maps `\pi: C \to D` ("projection") and `\iota: D \to
C` ("inclusion") such that

- `H` is a chain homotopy between `1_C` and `\iota \pi`,
- `\pi \iota = 1_D`,
- `\pi H = 0`,
- `H \iota = 0`,
- `H H = 0`.

Such a chain homotopy provides a strong relation between the chain
complexes `C` and `D`; for example, their homology groups are
isomorphic.
"""

########################################################################
#       Copyright (C) 2015 John H. Palmieri <palmieri@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
########################################################################

from sage.categories.morphism import Morphism
from sage.categories.homset import Hom
from sage.homology.chain_complex_morphism import ChainComplexMorphism


# In a perfect world, this would inherit from something like
# "TwoMorphism" rather than "Morphism"...
class ChainHomotopy(Morphism):
    r"""
    A chain homotopy.

    A chain homotopy `H` between chain maps `f, g: C \to D` is a sequence
    of maps `H_n: C_n \to D_{n+1}` (if the chain complexes are graded
    homologically) satisfying

    .. MATH::

       \partial_D H + H \partial_C = f - g.

    INPUT:

    - ``matrices`` -- dictionary of matrices, keyed by dimension
    - ``f`` -- chain map `C \to D`
    - ``g`` -- (optional) chain map `C \to D`

    The dictionary ``matrices`` defines ``H`` by specifying the matrix
    defining it in each degree: the entry `m` corresponding to key `i`
    gives the linear transformation `C_i \to D_{i+1}`.

    If `f` is specified but not `g`, then `g` can be recovered from
    the defining formula. That is, if `g` is not specified, then it
    is defined to be `f - \partial_D H - H \partial_C`.

    Note that the degree of the differential on the chain complex `C`
    must agree with that for `D`, and those degrees determine the
    "degree" of the chain homotopy map: if the degree of the
    differential is `d`, then the chain homotopy consists of a
    sequence of maps `C_n \to C_{n-d}`. The keys in the dictionary
    ``matrices`` specify the starting degrees.

    EXAMPLES::

        sage: from sage.homology.chain_homotopy import ChainHomotopy
        sage: C = ChainComplex({0: identity_matrix(ZZ, 1)})
        sage: D = ChainComplex({0: zero_matrix(ZZ, 1)})
        sage: f = Hom(C,D)({0: identity_matrix(ZZ, 1), 1: zero_matrix(ZZ, 1)})
        sage: g = Hom(C,D)({0: zero_matrix(ZZ, 1), 1: zero_matrix(ZZ, 1)})
        sage: H = ChainHomotopy({0: zero_matrix(ZZ, 0, 1), 1: identity_matrix(ZZ, 1)}, f, g)

    Note that the maps `f` and `g` are stored in the attributes ``H._f``
    and ``H._g``::

        sage: H._f
        Chain complex morphism:
          From: Chain complex with at most 2 nonzero terms over Integer Ring
          To: Chain complex with at most 2 nonzero terms over Integer Ring
        sage: H._f.in_degree(0)
        [1]
        sage: H._g.in_degree(0)
        [0]

    A non-example::

        sage: H = ChainHomotopy({0: zero_matrix(ZZ, 0, 1), 1: zero_matrix(ZZ, 1)}, f, g)
        Traceback (most recent call last):
        ...
        ValueError: the data do not define a valid chain homotopy
    """
    def __init__(self, matrices, f, g=None):
        r"""
        Create a chain homotopy between the given chain maps
        from a dictionary of matrices.

        EXAMPLES:

        If ``g`` is not specified, it is set equal to
        `f - (H \partial + \partial H)`. ::

            sage: from sage.homology.chain_homotopy import ChainHomotopy
            sage: C = ChainComplex({1: matrix(ZZ, 1, 2, (1,0)), 2: matrix(ZZ, 2, 1, (0, 2))},
            ....:                  degree_of_differential=-1)
            sage: D = ChainComplex({2: matrix(ZZ, 1, 1, (6,))}, degree_of_differential=-1)
            sage: f_d = {1: matrix(ZZ, 1, 2, (0,3)), 2: identity_matrix(ZZ, 1)}
            sage: f = Hom(C,D)(f_d)
            sage: H_d = {0: identity_matrix(ZZ, 1), 1: matrix(ZZ, 1, 2, (2,2))}
            sage: H = ChainHomotopy(H_d, f)
            sage: H._g.in_degree(0)
            []
            sage: H._g.in_degree(1)
            [-13  -9]
            sage: H._g.in_degree(2)
            [-3]

        TESTS:

        Try to construct a chain homotopy in which the maps do not
        have matching domains and codomains::

            sage: g = Hom(C,C)({}) # the zero chain map
            sage: H = ChainHomotopy(H_d, f, g)
            Traceback (most recent call last):
            ...
            ValueError: the chain maps are not compatible
        """
        domain = f.domain()
        codomain = f.codomain()
        deg = domain.degree_of_differential()
        # Check that the chain complexes are compatible. This should
        # never arise, because first there should be errors in
        # constructing the chain maps. But just in case...
        if domain.degree_of_differential() != codomain.degree_of_differential():
            raise ValueError('the chain complexes are not compatible')
        if g is not None:
            # Check that the chain maps are compatible.
            if not (domain == g.domain() and codomain ==
                    g.codomain()):
                raise ValueError('the chain maps are not compatible')
            # Check that the data define a chain homotopy.
            for i in domain.differential():
                if i in matrices and i+deg in matrices:
                    if not (codomain.differential(i-deg) * matrices[i] + matrices[i+deg] * domain.differential(i) == f.in_degree(i) - g.in_degree(i)):
                        raise ValueError('the data do not define a valid chain homotopy')
                elif i in matrices:
                    if not (codomain.differential(i-deg) * matrices[i] == f.in_degree(i) - g.in_degree(i)):
                        raise ValueError('the data do not define a valid chain homotopy')
                elif i+deg in matrices:
                    if not (matrices[i+deg] * domain.differential(i) == f.in_degree(i) - g.in_degree(i)):
                        raise ValueError('the data do not define a valid chain homotopy')
        else:
            # Define g.
            g_data = {}
            for i in domain.differential():
                if i in matrices and i+deg in matrices:
                    g_data[i] = f.in_degree(i) - matrices[i+deg] * domain.differential(i) - codomain.differential(i-deg) * matrices[i]
                elif i in matrices:
                    g_data[i] = f.in_degree(i) - codomain.differential(i-deg) * matrices[i]
                elif i+deg in matrices:
                    g_data[i] = f.in_degree(i) - matrices[i+deg] * domain.differential(i)
            g = ChainComplexMorphism(g_data, domain, codomain)
        self._matrix_dictionary = {}
        for i in matrices:
            m = matrices[i]
            # Use immutable matrices because they're hashable.
            m.set_immutable()
            self._matrix_dictionary[i] = m
        self._f = f
        self._g = g
        Morphism.__init__(self, Hom(domain, codomain))

    def is_algebraic_gradient_vector_field(self):
        r"""
        An algebraic gradient vector field is a linear map
        `H: C \to C` such that `H H = 0`.

        (Some authors also require that `H \partial H = H`, whereas
        some make this part of the definition of "homology gradient
        vector field. We have made the second choice.) See
        Molina-Abril and Réal [MAR2009]_ and Réal and Molina-Abril
        [RMA2009]_ for this and related terminology.

        See also :meth:`is_homology_gradient_vector_field`.

        EXAMPLES::

            sage: from sage.homology.chain_homotopy import ChainHomotopy
            sage: C = ChainComplex({0: zero_matrix(ZZ, 1), 1: identity_matrix(ZZ, 1)})

        The chain complex `C` is chain homotopy equivalent to a copy of
        `\ZZ` in degree 0. Two chain maps `C \to C` will be chain
        homotopic as long as they agree in degree 0. ::

            sage: f = Hom(C,C)({0: identity_matrix(ZZ, 1),
            ....:               1: matrix(ZZ, 1, 1, [3]),
            ....:               2: matrix(ZZ, 1, 1, [3])})
            sage: g = Hom(C,C)({0: identity_matrix(ZZ, 1),
            ....:               1: matrix(ZZ, 1, 1, [2]),
            ....:               2: matrix(ZZ, 1, 1, [2])})
            sage: H = ChainHomotopy({0: zero_matrix(ZZ, 0, 1),
            ....:                    1: zero_matrix(ZZ, 1),
            ....:                    2: identity_matrix(ZZ, 1)}, f, g)
            sage: H.is_algebraic_gradient_vector_field()
            True

        A chain homotopy which is not an algebraic gradient vector field::

            sage: H = ChainHomotopy({0: zero_matrix(ZZ, 0, 1),
            ....:                    1: identity_matrix(ZZ, 1),
            ....:                    2: identity_matrix(ZZ, 1)}, f, g)
            sage: H.is_algebraic_gradient_vector_field()
            False
        """
        if self.domain() != self.codomain():
            return False
        deg = self.domain().degree_of_differential()
        matrices = self._matrix_dictionary
        for i in matrices:
            if i-deg in matrices:
                if matrices[i-deg] * matrices[i] != 0:
                    return False
        return True

    def is_homology_gradient_vector_field(self):
        r"""
        A homology gradient vector field is an algebraic gradient vector
        field `H: C \to C` (i.e., a chain homotopy satisfying `H
        H = 0`) such that `\partial H \partial = \partial` and `H
        \partial H = H`.

        See Molina-Abril and Réal [MAR2009]_ and Réal and Molina-Abril
        [RMA2009]_ for this and related terminology.

        See also :meth:`is_algebraic_gradient_vector_field`.

        EXAMPLES::

            sage: from sage.homology.chain_homotopy import ChainHomotopy
            sage: C = ChainComplex({0: zero_matrix(ZZ, 1), 1: identity_matrix(ZZ, 1)})

            sage: f = Hom(C,C)({0: identity_matrix(ZZ, 1),
            ....:               1: matrix(ZZ, 1, 1, [3]),
            ....:               2: matrix(ZZ, 1, 1, [3])})
            sage: g = Hom(C,C)({0: identity_matrix(ZZ, 1),
            ....:               1: matrix(ZZ, 1, 1, [2]),
            ....:               2: matrix(ZZ, 1, 1, [2])})
            sage: H = ChainHomotopy({0: zero_matrix(ZZ, 0, 1),
            ....:                    1: zero_matrix(ZZ, 1),
            ....:                    2: identity_matrix(ZZ, 1)}, f, g)
            sage: H.is_homology_gradient_vector_field()
            True
        """
        if not self.is_algebraic_gradient_vector_field():
            return False
        deg = self.domain().degree_of_differential()
        matrices = self._matrix_dictionary
        for i in matrices:
            if i+deg in matrices:
                diff_i = self.domain().differential(i)
                if diff_i * matrices[i+deg] * diff_i != diff_i:
                    return False
            if matrices[i] * self.domain().differential(i-deg) * matrices[i] != matrices[i]:
                return False
        return True

    def in_degree(self, n):
        """
        The matrix representing this chain homotopy in degree ``n``.

        INPUT:

        - ``n`` -- degree

        EXAMPLES::

            sage: from sage.homology.chain_homotopy import ChainHomotopy
            sage: C = ChainComplex({1: matrix(ZZ, 0, 2)})  # one nonzero term in degree 1
            sage: D = ChainComplex({0: matrix(ZZ, 0, 1)})  # one nonzero term in degree 0
            sage: f = Hom(C, D)({})
            sage: H = ChainHomotopy({1: matrix(ZZ, 1, 2, (3,1))}, f, f)
            sage: H.in_degree(1)
            [3 1]

        This returns an appropriately sized zero matrix if the chain
        homotopy is not defined in degree n::

            sage: H.in_degree(-3)
            []
        """
        try:
            return self._matrix_dictionary[n]
        except KeyError:
            from sage.matrix.constructor import zero_matrix
            deg = self.domain().degree_of_differential()
            rows = self.codomain().free_module_rank(n-deg)
            cols = self.domain().free_module_rank(n)
            return zero_matrix(self.domain().base_ring(), rows, cols)

    def dual(self):
        r"""
        Dual chain homotopy to this one.

        That is, if this one is a chain homotopy between chain maps
        `f, g: C \to D`, then its dual is a chain homotopy between the
        dual of `f` and the dual of `g`, from `D^*` to `C^*`. It is
        represented in each degree by the transpose of the
        corresponding matrix.

        EXAMPLES::

            sage: from sage.homology.chain_homotopy import ChainHomotopy
            sage: C = ChainComplex({1: matrix(ZZ, 0, 2)})  # one nonzero term in degree 1
            sage: D = ChainComplex({0: matrix(ZZ, 0, 1)})  # one nonzero term in degree 0
            sage: f = Hom(C, D)({})
            sage: H = ChainHomotopy({1: matrix(ZZ, 1, 2, (3,1))}, f, f)
            sage: H.in_degree(1)
            [3 1]
            sage: H.dual().in_degree(0)
            [3]
            [1]
        """
        matrix_dict = self._matrix_dictionary
        deg = self.domain().degree_of_differential()
        matrices = {i-deg: matrix_dict[i].transpose() for i in matrix_dict}
        return ChainHomotopy(matrices, self._f.dual(), self._g.dual())

    def __hash__(self):
        """
        TESTS::

            sage: from sage.homology.chain_homotopy import ChainHomotopy
            sage: C = ChainComplex({1: matrix(ZZ, 0, 2)})  # one nonzero term in degree 1
            sage: D = ChainComplex({0: matrix(ZZ, 0, 1)})  # one nonzero term in degree 0
            sage: f = Hom(C, D)({})
            sage: H = ChainHomotopy({1: matrix(ZZ, 1, 2, (3,1))}, f, f)
            sage: hash(H)  # random
            314159265358979
        """
        return hash(self._f) ^ hash(self._g) ^ hash(tuple(self._matrix_dictionary.items()))

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: from sage.homology.chain_homotopy import ChainHomotopy
            sage: C = ChainComplex({0: identity_matrix(ZZ, 1)})
            sage: D = ChainComplex({0: zero_matrix(ZZ, 1)})
            sage: f = Hom(C,D)({0: identity_matrix(ZZ, 1), 1: zero_matrix(ZZ, 1)})
            sage: g = Hom(C,D)({0: zero_matrix(ZZ, 1), 1: zero_matrix(ZZ, 1)})
            sage: ChainHomotopy({0: zero_matrix(ZZ, 0, 1), 1: identity_matrix(ZZ, 1)}, f, g)
            Chain homotopy between:
              Chain complex morphism:
                From: Chain complex with at most 2 nonzero terms over Integer Ring
                To: Chain complex with at most 2 nonzero terms over Integer Ring
              and Chain complex morphism:
                From: Chain complex with at most 2 nonzero terms over Integer Ring
                To: Chain complex with at most 2 nonzero terms over Integer Ring
        """
        s = 'Chain homotopy between:'
        s += '\n  {}'.format('\n  '.join(self._f._repr_().split('\n')))
        s += '\n  and {}'.format('\n  '.join(self._g._repr_().split('\n')))
        return s


class ChainContraction(ChainHomotopy):
    r"""
    A chain contraction.

    An algebraic gradient vector field `H: C \to C` (that is a chain
    homotopy satisfying `H H = 0`) for which there are chain
    maps `\pi: C \to D` ("projection") and `\iota: D \to C`
    ("inclusion") such that

    - `H` is a chain homotopy between `1_C` and `\iota \pi`,
    - `\pi \iota = 1_D`,
    - `\pi H = 0`,
    - `H \iota = 0`.

    ``H`` is defined by a dictionary ``matrices`` of matrices.

    INPUT:

    - ``matrices`` -- dictionary of matrices, keyed by dimension
    - ``pi`` -- a chain map `C \to D`
    - ``iota`` -- a chain map `D \to C`

    EXAMPLES::

        sage: from sage.homology.chain_homotopy import ChainContraction
        sage: C = ChainComplex({0: zero_matrix(ZZ, 1), 1: identity_matrix(ZZ, 1)})
        sage: D = ChainComplex({0: matrix(ZZ, 0, 1)})

    The chain complex `C` is chain homotopy equivalent to `D`, which is just
    a copy of `\ZZ` in degree 0, and we construct a chain contraction::

        sage: pi = Hom(C,D)({0: identity_matrix(ZZ, 1)})
        sage: iota = Hom(D,C)({0: identity_matrix(ZZ, 1)})
        sage: H = ChainContraction({0: zero_matrix(ZZ, 0, 1),
        ....:                       1: zero_matrix(ZZ, 1),
        ....:                       2: identity_matrix(ZZ, 1)}, pi, iota)
    """
    def __init__(self, matrices, pi, iota):
        r"""
        Create a chain contraction from the given data.

        EXAMPLES::

            sage: from sage.homology.chain_homotopy import ChainContraction
            sage: C = ChainComplex({0: zero_matrix(ZZ, 1), 1: identity_matrix(ZZ, 1)})
            sage: D = ChainComplex({0: matrix(ZZ, 0, 1)})

        The chain complex `C` is chain homotopy equivalent to `D`,
        which is just a copy of `\ZZ` in degree 0, and we try
        construct a chain contraction, but get the map `\iota` wrong::

            sage: pi = Hom(C,D)({0: identity_matrix(ZZ, 1)})
            sage: iota = Hom(D,C)({0: zero_matrix(ZZ, 1)})
            sage: H = ChainContraction({0: zero_matrix(ZZ, 0, 1), 1: zero_matrix(ZZ, 1), 2: identity_matrix(ZZ, 1)}, pi, iota)
            Traceback (most recent call last):
            ...
            ValueError: the composite 'pi iota' is not the identity

        Another bad `\iota`::

            sage: iota = pi  # wrong domain, codomain
            sage: H = ChainContraction({0: zero_matrix(ZZ, 0, 1), 1: zero_matrix(ZZ, 1), 2: identity_matrix(ZZ, 1)}, pi, iota)
            Traceback (most recent call last):
            ...
            ValueError: the chain maps are not composable

        `\iota` is okay, but wrong data defining `H`::

            sage: iota = Hom(D,C)({0: identity_matrix(ZZ, 1)})
            sage: H = ChainContraction({0: zero_matrix(ZZ, 0, 1), 1: identity_matrix(ZZ, 1), 2: identity_matrix(ZZ, 1)}, pi, iota)
            Traceback (most recent call last):
            ...
            ValueError: not an algebraic gradient vector field
        """
        from sage.matrix.constructor import identity_matrix
        from .chain_complex_morphism import ChainComplexMorphism

        if not (pi.domain() == iota.codomain()
                and pi.codomain() == iota.domain()):
            raise ValueError('the chain maps are not composable')
        C = pi.domain()
        D = pi.codomain()
        base_ring = C.base_ring()

        # Check that the composite 'pi iota' is 1.
        for i in D.nonzero_degrees():
            if pi.in_degree(i) * iota.in_degree(i) != identity_matrix(base_ring, D.free_module_rank(i)):
                raise ValueError("the composite 'pi iota' is not the identity")

        # Construct the chain map 'id_C'.
        id_C_dict = {}
        for i in C.nonzero_degrees():
            id_C_dict[i] = identity_matrix(base_ring, C.free_module_rank(i))
        id_C = ChainComplexMorphism(id_C_dict, C, C)

        # Now check that
        # - `H` is a chain homotopy between `id_C` and `\iota \pi`
        # - `HH = 0`
        ChainHomotopy.__init__(self, matrices, id_C, iota * pi)
        if not self.is_algebraic_gradient_vector_field():
            raise ValueError('not an algebraic gradient vector field')
        # Check that `\pi H = 0`:
        deg = C.degree_of_differential()
        for i in matrices:
            if pi.in_degree(i-deg) * matrices[i] != 0:
                raise ValueError('the data do not define a valid chain contraction: pi H != 0')
        # Check that `H \iota = 0`:
        for i in iota._matrix_dictionary:
            if i in matrices:
                if matrices[i] * iota.in_degree(i) != 0:
                    raise ValueError('the data do not define a valid chain contraction: H iota != 0')
        self._pi = pi
        self._iota = iota

    def pi(self):
        r"""
        The chain map `\pi` associated to this chain contraction.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: S2 = simplicial_complexes.Sphere(2)
            sage: phi, M = S2.algebraic_topological_model(QQ)
            sage: phi.pi()
            Chain complex morphism:
              From: Chain complex with at most 3 nonzero terms over Rational Field
              To: Chain complex with at most 3 nonzero terms over Rational Field
            sage: phi.pi().in_degree(0)  # Every vertex represents a homology class.
            [1 1 1 1]
            sage: phi.pi().in_degree(1)  # No homology in degree 1.
            []

        The degree 2 homology generator is detected on a single simplex::

            sage: phi.pi().in_degree(2)                                                 # needs sage.graphs
            [0 0 0 1]
        """
        return self._pi

    def iota(self):
        r"""
        The chain map `\iota` associated to this chain contraction.

        EXAMPLES::

            sage: S2 = simplicial_complexes.Sphere(2)                                   # needs sage.graphs
            sage: phi, M = S2.algebraic_topological_model(QQ)                           # needs sage.graphs
            sage: phi.iota()                                                            # needs sage.graphs
            Chain complex morphism:
              From: Chain complex with at most 3 nonzero terms over Rational Field
              To: Chain complex with at most 3 nonzero terms over Rational Field

        Lifting the degree zero homology class gives a single vertex::

            sage: phi.iota().in_degree(0)                                               # needs sage.graphs
            [0]
            [0]
            [0]
            [1]

        Lifting the degree two homology class gives the signed sum of
        all of the 2-simplices::

            sage: phi.iota().in_degree(2)                                               # needs sage.graphs
            [-1]
            [ 1]
            [-1]
            [ 1]
        """
        return self._iota

    def dual(self):
        """
        The chain contraction dual to this one.

        This is useful when switching from homology to cohomology.

        EXAMPLES::

            sage: S2 = simplicial_complexes.Sphere(2)                                   # needs sage.graphs
            sage: phi, M = S2.algebraic_topological_model(QQ)                           # needs sage.graphs
            sage: phi.iota()                                                            # needs sage.graphs
            Chain complex morphism:
              From: Chain complex with at most 3 nonzero terms over Rational Field
              To: Chain complex with at most 3 nonzero terms over Rational Field

        Lifting the degree zero homology class gives a single vertex,
        but the degree zero cohomology class needs to be detected on
        every vertex, and vice versa for degree 2::

            sage: # needs sage.graphs
            sage: phi.iota().in_degree(0)
            [0]
            [0]
            [0]
            [1]
            sage: phi.dual().iota().in_degree(0)
            [1]
            [1]
            [1]
            [1]
            sage: phi.iota().in_degree(2)
            [-1]
            [ 1]
            [-1]
            [ 1]
            sage: phi.dual().iota().in_degree(2)
            [0]
            [0]
            [0]
            [1]
        """
        matrix_dict = self._matrix_dictionary
        deg = self.domain().degree_of_differential()
        matrices = {i-deg: matrix_dict[i].transpose() for i in matrix_dict}
        return ChainContraction(matrices, self.iota().dual(), self.pi().dual())
