"""
Moment-angle complexes

AUTHORS:

- Ognjen Petrov (2023-06-25): initial version
"""

# ****************************************************************************
#       Copyright (C) 2023 Ognjen Petrov <ognjenpetrov@yahoo.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.algebras import Algebras
from sage.categories.modules import Modules
from sage.categories.algebras import Algebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.homology.homology_group import HomologyGroup
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.sets.family import Family
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from .cubical_complex import CubicalComplex, cubical_complexes
from .simplicial_complex import SimplicialComplex, copy
from sage.topology import simplicial_complex_catalog as simplicial_complexes
from itertools import combinations


def _cubical_complex_union(c1, c2):
    """
    Return the union of cubical complexes.

    This method returns a cubical complex whose set of maximal faces
    is the union of sets of maximal faces of ``c1`` and ``c2``.

    INPUT:

    - ``c1``, ``c2`` -- a cubical complex

    OUTPUT: the union of cubical complexes ``c1`` and ``c2``

    .. WARNING::

        This is regular union, not disjoint union. One should be careful
        with the nomenclature of the vertices.

    EXAMPLES::

        sage: from sage.topology.moment_angle_complex import (
        ....:      _cubical_complex_union as union
        ....: )
        sage: C1 = CubicalComplex([([0,0], [2,3]), ([0,1], [3,3]),
        ....:                      ([0,1], [2,2]), ([1,1], [2,3])]); C1
        Cubical complex with 4 vertices and 8 cubes
        sage: C2 = CubicalComplex([([0,0], [2,3]), ([0,1], [3,3]),
        ....:                      ([0,1], [2,2]), ([2,2], [2,3])]); C2
        Cubical complex with 6 vertices and 10 cubes
        sage: union(C1, C2)
        Cubical complex with 6 vertices and 11 cubes
        sage: union(C1, C1) == C1
        True
    """
    facets = list(c1.maximal_cells())
    facets.extend(c2.maximal_cells())
    return CubicalComplex(facets)


class MomentAngleComplex(UniqueRepresentation, SageObject):
    r"""
    A moment-angle complex.

    Given a simplicial complex `K`, with a set of vertices
    `V = \{v_1, v_2, \ldots, v_n\}`, a moment-angle complex over `K` is a
    topological space `Z`, which is a union of `X_{\sigma}`, where
    `\sigma \in K`, and `X_{\sigma} = Y_{v_1} \times Y_{v_2} \times \cdots
    \times Y_{v_n}` and `Y_{v_i}` is a 2-disk (a 2-simplex) if
    `v_i \in \sigma` , or a 1-sphere otherwise.

    .. MATH::

        Y_{v_i} =
        \begin{cases}
            D^2, &v_i \in \sigma,\\
            S^1, &v_i \notin \sigma.
        \end{cases}

    .. NOTE::

        The mentioned union is not a disjoint union of topological spaces.
        The unit disks and the unit spheres are considered subsets of `\CC`,
        so the union is just a normal union of subsets of `\CC^n`.

    Here we view moment-angle complexes as cubical complexes and
    try to compute mostly things which would not require computing
    the moment-angle complex itself, but rather work with the
    corresponding simplicial complex.

    .. NOTE::

        One of the more useful properties will be the
        :meth:`bigraded Betti numbers
        <sage.topology.simplicial_complex.bigraded_betti_numbers>`,
        and the underlying theorem which makes this possible is Hochter's
        formula, which can be found on page 104 of [BP2014]_.

    INPUT:

    - ``simplicial_complex`` -- an instance of ``SimplicialComplex``,
      or an object from which an instance of ``SimplicialComplex`` can be
      created (e.g., list of facets), which represents the associated
      simplicial complex over which this moment-angle complex is created

    EXAMPLES::

        sage: MomentAngleComplex([[1,2,3], [2,4], [3,4]])
        Moment-angle complex of Simplicial complex with vertex set
        (1, 2, 3, 4) and facets {(2, 4), (3, 4), (1, 2, 3)}
        sage: X = SimplicialComplex([[0,1], [1,2], [1,3], [2,3]])
        sage: Z = MomentAngleComplex(X); Z
        Moment-angle complex of Simplicial complex with vertex set
        (0, 1, 2, 3) and facets {(0, 1), (1, 2), (1, 3), (2, 3)}
        sage: M = MomentAngleComplex([[1], [2]]); M
        Moment-angle complex of Simplicial complex with vertex set
        (1, 2) and facets {(1,), (2,)}

    We can perform a number of operations, such as find the dimension or
    compute the homology::

        sage: M.homology()
        {0: 0, 1: 0, 2: 0, 3: Z}
        sage: Z.dimension()
        6
        sage: Z.homology()
        {0: 0, 1: 0, 2: 0, 3: Z x Z, 4: Z, 5: Z, 6: Z}

    If the associated simplicial complex is an `n`-simplex, then the
    corresponding moment-angle complex is a polydisc (a complex ball) of
    complex dimension `n+1`::

        sage: Z = MomentAngleComplex([[0, 1, 2]]); Z
        Moment-angle complex of Simplicial complex with vertex set (0, 1, 2)
        and facets {(0, 1, 2)}

    This can be seen by viewing the components used in the construction
    of this moment-angle complex by calling :meth:`components()`::

        sage: Z.components()
        {(0, 1, 2): [The 2-simplex, The 2-simplex, The 2-simplex]}

    If the associated simplicial complex is a disjoint union of 2 points,
    then the corresponding moment-angle complex is homeomorphic to a boundary
    of a 3-sphere::

        sage: Z = MomentAngleComplex([[0], [1]]); Z
        Moment-angle complex of Simplicial complex with vertex set
        (0, 1) and facets {(0,), (1,)}
        sage: dict(sorted(Z.components().items()))
        {(0,): [The 2-simplex, Minimal triangulation of the 1-sphere],
         (1,): [Minimal triangulation of the 1-sphere, The 2-simplex]}

    The moment-angle complex passes all the tests of the test suite relative
    to its category::

        sage: TestSuite(Z).run()
    """
    @staticmethod
    def __classcall_private__(cls, simplicial_complex):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: Z = MomentAngleComplex([[0,2], [1,2,3]])
            sage: W = MomentAngleComplex([[0,2], [1,2,3]])
            sage: Z is W
            True
            sage: Z is MomentAngleComplex(Z)
            True
        """
        if simplicial_complex:
            if isinstance(simplicial_complex, MomentAngleComplex):
                # Allows for copy constructor
                immutable_complex = SimplicialComplex(simplicial_complex._simplicial_complex, is_mutable=False)
            elif not isinstance(simplicial_complex, SimplicialComplex):
                # Try to create a SimplicialComplex out of simplicial_complex
                # in case that simplicial_complex is a list of facets, or
                # something that can generate a SimplicialComplex
                immutable_complex = SimplicialComplex(simplicial_complex, is_mutable=False)
            elif simplicial_complex.is_mutable():
                immutable_complex = SimplicialComplex(simplicial_complex, is_mutable=False)
            else:
                immutable_complex = simplicial_complex
        else:
            immutable_complex = SimplicialComplex(is_mutable=False)
        return super().__classcall__(cls, immutable_complex)

    def __init__(self, simplicial_complex):
        """
        Initialize ``self``.

        TESTS::

            sage: Z = MomentAngleComplex([[0,1,2], [1,2,3], [0, 3]])
            sage: TestSuite(Z).run()
        """
        # The underlying simplicial complex
        self._simplicial_complex = copy(simplicial_complex)
        # A dictionary of components indexed by facets
        self._components = {}

        vertices = self._simplicial_complex.vertices()
        # it suffices to perform union only over facets
        for facet in self._simplicial_complex.maximal_faces():
            Y = []
            for j in vertices:
                if j in facet:
                    Y.append(simplicial_complexes.Simplex(2))
                else:
                    Y.append(simplicial_complexes.Sphere(1))

            self._components[facet] = Y

    @lazy_attribute
    def _moment_angle_complex(self):
        """
        Create the moment-angle complex as a cubical complex.

        If this lazy attribute is accessed, we explicitly compute
        the moment-angle complex, viewed as a cubical complex.

        .. WARNING::

            The construction can be very slow, it is not reccomended unless
            the corresponding simplicial complex has 5 or less vertices.

        TESTS::

            sage: Z = MomentAngleComplex([[0], [1], [2]]); Z
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2) and facets {(0,), (1,), (2,)}
            sage: Z._moment_angle_complex
            Cubical complex with 64 vertices and 705 cubes

        This is called by :meth:`cubical_complex()`::

            sage: Z.cubical_complex()
            Cubical complex with 64 vertices and 705 cubes
            sage: Z.cubical_complex() == Z._moment_angle_complex
            True
        """
        n = len(self._simplicial_complex.vertices())
        D = [cubical_complexes.Cube(2)] * n
        S = [cubical_complexes.Sphere(1)] * n

        moment_angle_complex = CubicalComplex()
        for component in self._components.values():
            x = D[0] if component[0] == simplicial_complexes.Simplex(2) else S[0]
            for j in range(1, len(component)):
                y = D[j] if component[j] == simplicial_complexes.Simplex(2) else S[j]
                x = x.product(y)
            moment_angle_complex = _cubical_complex_union(moment_angle_complex, x)

        return moment_angle_complex

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: Z = MomentAngleComplex([[0,1], [1,2], [2,0]])
            sage: Z._repr_()
            'Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2) and facets {(0, 1), (0, 2), (1, 2)}'
            sage: repr(Z)
            'Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2) and facets {(0, 1), (0, 2), (1, 2)}'
            sage: Z
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2) and facets {(0, 1), (0, 2), (1, 2)}
            sage: Z = MomentAngleComplex([[i for i in range(20)]])
            sage: Z._repr_()
            'Moment-angle complex of Simplicial complex with
            20 vertices and 1 facets'
        """
        return "Moment-angle complex of " + repr(self._simplicial_complex)

    def cubical_complex(self):
        """
        Return the cubical complex that represents ``self``.

        This method returns returns a cubical complex which is
        derived by explicitly computing products and unions in the
        definition of a moment-angle complex.

        .. WARNING::

            The construction can be very slow, it is not reccomended unless
            the corresponding simplicial complex has 5 or less vertices.

        EXAMPLES::

            sage: Z = MomentAngleComplex([[0,1,2], [1,3]]); Z
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2, 3) and facets {(1, 3), (0, 1, 2)}
            sage: Z.cubical_complex()
            Cubical complex with 256 vertices and 6409 cubes
            sage: dim(Z.cubical_complex()) == dim(Z)
            True
            sage: Z = MomentAngleComplex([[0,1], [1,2], [2,0], [1,3]]); Z
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2, 3) and facets {(0, 1), (0, 2), (1, 2), (1, 3)}
            sage: Z.betti() == Z.cubical_complex().betti()  # long time
            True

        We can now work with moment-angle complexes as concrete cubical
        complexes. Though, it can be very slow, due to the size of the
        complex. However, for some smaller moment-angle complexes, this
        may be possible::

            sage: Z = MomentAngleComplex([[0], [1]]); Z
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1) and facets {(0,), (1,)}
            sage: Z.cubical_complex().f_vector()
            [1, 16, 32, 24, 8]
        """
        return self._moment_angle_complex

    def simplicial_complex(self):
        """
        Return the simplicial complex that defines ``self``.

        EXAMPLES::

            sage: RP2 = simplicial_complexes.RealProjectivePlane()
            sage: Z = MomentAngleComplex(RP2)
            sage: Z.simplicial_complex()
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and 10 facets
            sage: Z = MomentAngleComplex([[0], [1], [2]])
            sage: Z.simplicial_complex()
            Simplicial complex with vertex set (0, 1, 2)
            and facets {(0,), (1,), (2,)}
        """
        return self._simplicial_complex

    def components(self):
        r"""
        Return the dictionary of components of ``self``, indexed by facets
        of the associated simplicial complex.

        OUTPUT:

        A dictonary, whose values are lists, representing spheres
        and disks described in the construction of the moment-angle
        complex. ``The 2-simplex`` represents a 2-disk, and
        ``Minimal triangulation of the 1-sphere`` represents a 1-sphere.

        EXAMPLES::

            sage: M = MomentAngleComplex([[0, 1, 2]]); M
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2) and facets {(0, 1, 2)}
            sage: M.components()
            {(0, 1, 2): [The 2-simplex, The 2-simplex, The 2-simplex]}
            sage: Z = MomentAngleComplex([[0], [1]]); Z
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1) and facets {(0,), (1,)}
            sage: sorted(Z.components().items())
            [((0,), [The 2-simplex, Minimal triangulation of the 1-sphere]),
             ((1,), [Minimal triangulation of the 1-sphere, The 2-simplex])]

        We interpret the output of this method by taking the product
        of all the elements in each list, and then taking the union
        of all products. From the previous example, we have
        `\mathcal{Z} = S^1 \times D^2 \cup D^2 \times S^1
        = \partial (D^2 \times D^2) = \partial D^4 = S^3`::

            sage: Z = MomentAngleComplex([[0,1], [1,2], [2,3], [3,0]])
            sage: sorted(Z.components().items())
            [((0, 1),
              [The 2-simplex,
               The 2-simplex,
               Minimal triangulation of the 1-sphere,
               Minimal triangulation of the 1-sphere]),
             ((0, 3),
              [The 2-simplex,
               Minimal triangulation of the 1-sphere,
               Minimal triangulation of the 1-sphere,
               The 2-simplex]),
             ((1, 2),
              [Minimal triangulation of the 1-sphere,
               The 2-simplex,
               The 2-simplex,
               Minimal triangulation of the 1-sphere]),
             ((2, 3),
              [Minimal triangulation of the 1-sphere,
               Minimal triangulation of the 1-sphere,
               The 2-simplex,
               The 2-simplex])]

        It is not that difficult to prove that the previous
        moment-angle complex is homeomorphic to a product of
        two 3-spheres. We can look at the cohomologies to try
        and validate whether this makes sense::

            sage: S3 = simplicial_complexes.Sphere(3)
            sage: product_of_spheres = S3.product(S3)
            sage: Z.cohomology()
            {0: 0, 1: 0, 2: 0, 3: Z x Z, 4: 0, 5: 0, 6: Z}
            sage: Z.cohomology() == product_of_spheres.cohomology()  # long time
            True
        """
        return self._components

    def dimension(self):
        r"""
        The dimension of ``self``.

        The dimension of a moment-angle complex is the dimension
        of the constructed (cubical) complex. It is not difficult
        to see that this turns out to be `m+n+1`, where `m` is the
        number of vertices and `n` is the dimension of the associated
        simplicial complex.

        EXAMPLES::

            sage: Z = MomentAngleComplex([[0,1], [1,2,3]])
            sage: Z.dimension()
            7
            sage: Z = MomentAngleComplex([[0, 1, 2]])
            sage: Z.dimension()
            6
            sage: dim(Z)
            6

        We can construct the cubical complex and compare whether
        the dimensions coincide::

            sage: dim(Z) == dim(Z.cubical_complex())
            True
        """
        number_of_vertices = len(self._simplicial_complex.vertices())
        dim = self._simplicial_complex.dimension()
        return number_of_vertices + dim + 1

    @cached_method  # maybe ignore the algorithm?
    def _homology_group(self, i, base_ring, cohomology, algorithm, verbose, reduced):
        """
        The `i`-th (reduced) homology group of ``self``.

        .. SEEALSO::

            :meth:`homology`,
            :meth:`.cell_complex.GenericCellComplex.homology`.

        TESTS::

            sage: Z = MomentAngleComplex([[0,1,2], [1,2,3]]); Z
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2, 3) and facets {(0, 1, 2), (1, 2, 3)}
            sage: Z._homology_group(3, base_ring=ZZ,
            ....:                   reduced=True, verbose=False,
            ....:                   cohomology=False, algorithm='pari')
            Z
            sage: Z._homology_group(4, base_ring=ZZ,
            ....:                   reduced=True, verbose=False,
            ....:                   cohomology=False, algorithm='pari')
            0
            sage: Z.homology()
            {0: 0, 1: 0, 2: 0, 3: Z, 4: 0, 5: 0, 6: 0, 7: 0}
            sage: RP2 = simplicial_complexes.RealProjectivePlane()
            sage: Z = MomentAngleComplex(RP2)
            sage: Z._homology_group(8, base_ring=ZZ,
            ....:                   reduced=True, verbose=False,
            ....:                   cohomology=False, algorithm='pari')
            C2
        """
        if i == 0:
            # This is a special case when computing (co)homology
            if reduced:
                return HomologyGroup(0, base_ring)
            return HomologyGroup(1, base_ring)

        vertices = self._simplicial_complex.vertices()
        n = len(vertices)
        invfac = []

        for j in range(n+1):
            for x in combinations(vertices, j):
                S = self._simplicial_complex.generated_subcomplex(x)
                if base_ring.is_field():
                    invfac.append(S.homology(i-j-1, base_ring=base_ring,
                                             cohomology=cohomology, algorithm=algorithm,
                                             verbose=verbose, reduced=True).dimension())
                else:
                    invfac.extend(S.homology(i-j-1, base_ring=base_ring,
                                             cohomology=cohomology, algorithm=algorithm,
                                             verbose=verbose, reduced=True)._original_invts)

        if base_ring.is_field():
            return HomologyGroup(sum(invfac), base_ring)

        m = len(invfac)
        return HomologyGroup(m, base_ring, invfac)

    def homology(self, dim=None, base_ring=ZZ, cohomology=False,
                 algorithm='pari', verbose=False, reduced=True):
        r"""
        The (reduced) homology of ``self``.

        INPUT:

        - ``dim`` -- integer, or a list of integers; represents the
          homology (or homologies) we want to compute
        - ``base_ring`` -- commutative ring (default: ``ZZ``); must be ``ZZ``
          or a field
        - ``cohomology`` -- boolean (default: ``False``);
          if ``True``, compute cohomology rather than homology
        - ``algorithm`` -- string (default: ``'pari'``); the options are
          ``'auto'``, ``'dhsw'``, or ``'pari'``; see
          :meth:`.cell_complex.GenericCellComplex.homology` documentation
          for a description of what they mean
        - ``verbose`` -- boolean (default: ``False``); if ``True``,
          print some messages as the homology is computed
        - ``reduced`` -- boolean (default: ``True``); if ``True``,
          return the reduced homology

        ALGORITHM:

        This algorithm is adopted from Theorem 4.5.8 of [BP2014]_.

        The (co)homology of the moment-angle complex is closely related
        to the (co)homologies of certain full subcomplexes of the
        associated simplicial complex. More specifically, we know that:

        .. MATH::

            H_l(\mathcal{Z}_\mathcal{K}) \cong
            \bigoplus_{J \subseteq [m]} \widetilde{H}_{l-|J|-1}(\mathcal{K}_J),

        where `\mathcal{Z}_\mathcal{K}` denotes the moment-angle complex
        associated to a simplicial complex `\mathcal{K}`, on the set of
        vertices `\{1, 2, 3, \ldots, m\} =: [m]`. `\mathcal{K}_J` denotes the
        full subcomplex of `\mathcal{K}`, generated by a set of vertices `J`.
        The same formula holds true for cohomology groups as well.

        .. SEEALSO::

            :meth:`.cell_complex.GenericCellComplex.homology`

        EXAMPLES::

            sage: Z = MomentAngleComplex([[0,1,2], [1,2,3], [3,0]]); Z
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2, 3) and facets {(0, 3), (0, 1, 2), (1, 2, 3)}
            sage: Z = MomentAngleComplex([[0,1,2], [1,2,3], [3,0]])
            sage: Z.homology()
            {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: Z x Z, 6: Z, 7: 0}
            sage: Z.homology(base_ring=GF(2))
            {0: Vector space of dimension 0 over Finite Field of size 2,
             1: Vector space of dimension 0 over Finite Field of size 2,
             2: Vector space of dimension 0 over Finite Field of size 2,
             3: Vector space of dimension 0 over Finite Field of size 2,
             4: Vector space of dimension 0 over Finite Field of size 2,
             5: Vector space of dimension 2 over Finite Field of size 2,
             6: Vector space of dimension 1 over Finite Field of size 2,
             7: Vector space of dimension 0 over Finite Field of size 2}
           sage: RP = simplicial_complexes.RealProjectivePlane()
           sage: Z = MomentAngleComplex(RP)
           sage: Z.homology()
           {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: Z^10, 6: Z^15, 7: Z^6, 8: C2, 9: 0}

        This yields the same result as creating a cubical complex
        from this moment-angle complex, and then computing its (co)homology,
        but that is incomparably slower and is really only possible when
        the associated simplicial complex is very small::

            sage: Z = MomentAngleComplex([[0,1], [1,2], [2,0]]); Z
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2) and facets {(0, 1), (0, 2), (1, 2)}
            sage: Z.cubical_complex()
            Cubical complex with 64 vertices and 729 cubes
            sage: Z.cubical_complex().homology() == Z.homology()
            True

        Meanwhile, the homology computation used here is quite efficient
        and works well even with significantly larger underlying simplicial
        complexes::

            sage: Z = MomentAngleComplex([[0,1,2,3,4,5], [0,1,2,3,4,6],
            ....:                         [0,1,2,3,5,7], [0,1,2,3,6,8,9]])
            sage: Z.homology()  # long time
            {0: 0,
             1: 0,
             2: 0,
             3: Z^9,
             4: Z^17,
             5: Z^12,
             6: Z x Z x Z,
             7: 0,
             8: 0,
             9: 0,
             10: 0,
             11: 0,
             12: 0,
             13: 0,
             14: 0,
             15: 0,
             16: 0,
             17: 0}
            sage: Z = MomentAngleComplex([[0,1,2,3], [0,1,2,4], [0,1,3,5],
            ....:                         [0,1,4,5], [0,2,3,6], [0,2,4,6]])
            sage: Z.homology(dim=range(0,5), reduced=True)
            {0: 0, 1: 0, 2: 0, 3: Z x Z x Z x Z, 4: Z x Z}
            sage: Z.homology(dim=range(0,5), reduced=False)
            {0: Z, 1: 0, 2: 0, 3: Z x Z x Z x Z, 4: Z x Z}
            sage: all(Z.homology(i,reduced=True) == Z.homology(i,reduced=False)
            ....:     for i in range(1, dim(Z)))
            True
            sage: all(Z.homology(i,reduced=True) == Z.homology(i,reduced=False)
            ....:     for i in range(0, dim(Z)))
            False
        """
        if dim is not None:
            if isinstance(dim, (list, tuple, range)):
                low = min(dim)
                high = max(dim)
            else:
                low = dim
                high = dim
            dims = range(low, high + 1)
        else:
            dims = range(self.dimension()+1)

        answer = {i: self._homology_group(i, base_ring=base_ring, cohomology=cohomology,
                                          algorithm=algorithm, verbose=verbose, reduced=reduced) for i in dims}
        return answer

    def cohomology(self, dim=None, base_ring=ZZ, algorithm='pari',
                   verbose=False, reduced=True):
        r"""
        The reduced cohomology of ``self``.

        This is equivalent to calling the ``homology()`` method,
        with ``cohomology=True`` as an argument.

        .. SEEALSO::

            :meth:`homology`.

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1],[1,2],[2,3],[3,0]])
            sage: Z = MomentAngleComplex(X)

        It is known that the previous moment-angle complex is homeomorphic
        to a product of two 3-spheres (which can be seen by looking at the
        output of ``components()``)::

            sage: S3 = simplicial_complexes.Sphere(3)
            sage: product_of_spheres = S3.product(S3)
            sage: Z.cohomology()
            {0: 0, 1: 0, 2: 0, 3: Z x Z, 4: 0, 5: 0, 6: Z}
            sage: Z.cohomology() == product_of_spheres.cohomology()  # long time
            True
        """
        return self.homology(dim=dim, cohomology=True, base_ring=base_ring,
                             algorithm=algorithm, verbose=verbose, reduced=reduced)

    def betti(self, dim=None):
        r"""
        Return the Betti number (or numbers) of ``self``.

        The the `i`-th Betti number is the rank of the `i`-th homology group.

        INPUT:

        - ``dim`` -- (optional) an integer or a list of integers

        OUTPUT:

        If ``dim`` is an integer or a list of integers, then return
        a dictionary of Betti numbers for each given dimension, indexed
        by dimension. Otherwise, return all Betti numbers.

        EXAMPLES::

            sage: Z = MomentAngleComplex([[0,1], [1,2], [2,0], [1,2,3]])
            sage: Z.betti()
            {0: 1, 1: 0, 2: 0, 3: 1, 4: 0, 5: 1, 6: 1, 7: 0}
            sage: Z = MomentAngleComplex([[0,1], [1,2], [2,0], [1,2,3], [3,0]])
            sage: Z.betti(dim=6)
            {6: 2}
        """
        dict = {}
        H = self.homology(dim=dim, base_ring=QQ)
        try:
            for n in H.keys():
                dict[n] = H[n].dimension()
                if n == 0:
                    dict[n] += 1
            return dict
        except AttributeError:
            return H.dimension()

    def euler_characteristic(self):
        """
        Return the Euler characteristic of ``self``.

        The Euler characteristic is defined as the alternating sum
        of the Betti numbers of ``self``.

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1,2,3,4,5], [0,1,2,3,4,6],
            ....:                        [0,1,2,3,5,7], [0,1,2,3,6,8,9]])
            sage: M = MomentAngleComplex(X)
            sage: M.euler_characteristic()  # long time
            0
            sage: Z = MomentAngleComplex([[0,1,2,3,4]])
            sage: Z.euler_characteristic()
            1
        """
        betti_numbers = self.betti()
        return ZZ.sum((-1)**n * betti_numbers[n] for n in range(self.dimension() + 1))

    def product(self, other):
        """
        Return the product of ``self`` with ``other``.

        It is known that the product of two moment-angle complexes
        is a moment-angle complex over the join of the two corresponding
        simplicial complexes. This result can be found on page 138 of
        [BP2014]_.

        OUTPUT: a moment-angle complex which is the product of the
        parsed moment-angle complexes

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1,2,3], [1,4], [3,2,4]])
            sage: Y = SimplicialComplex([[1,2,3],[1,2,4],[3,5],[4,5]])
            sage: Z = MomentAngleComplex(X)
            sage: M = MomentAngleComplex(Y)
            sage: Z.product(M)
            Moment-angle complex of Simplicial complex with
            10 vertices and 12 facets
            sage: Z.product(M) == MomentAngleComplex(X*Y)
            True
        """
        simplicial_complex = self._simplicial_complex.join(other._simplicial_complex, rename_vertices=True)
        return MomentAngleComplex(simplicial_complex)

    def has_trivial_lowest_deg_massey_product(self):
        """
        Return whether ``self`` has a non-trivial lowest degree
        triple Massey product.

        This is the Massey product in the cohomology of this
        moment-angle complex. This relies on the theorem which was
        proven in [GL2019]_.

        ALGORITHM:

        We obtain the one-skeleton from the associated simplicial complex,
        which we consider to be a graph. We then perform ``subgraph_search``,
        searching for any subgraph isomorphic to one of the 8 obstruction
        graphs listed in the mentioned paper.

        EXAMPLES:

        A simplex will not have a trivial triple lowest-degree
        Massey product, because its one-skeleton certainly does
        contain a subcomplex isomorphic to one of the 8 mentioned
        in the paper::

            sage: Z = MomentAngleComplex([[1,2,3,4,5,6]])
            sage: Z.has_trivial_lowest_deg_massey_product()
            False

        The following is one of the 8 obstruction graphs::

            sage: Z = MomentAngleComplex([[1, 2], [1, 4], [2, 3], [3, 5],
            ....:                         [5, 6], [4, 5], [1, 6]])
            sage: Z.has_trivial_lowest_deg_massey_product()
            False

        A hexagon is not isomorphic to any of the 8 obstruction graphs::

            sage: Z = MomentAngleComplex([[0,1], [1,2], [2,3],
            ....:                         [3,4], [4,5], [5,0]])
            sage: Z.has_trivial_lowest_deg_massey_product()
            True
        """
        from sage.graphs.graph import Graph

        one_skeleton = self._simplicial_complex.graph()

        obstruction_graphs = [
            Graph([(1, 2), (1, 4), (2, 3), (3, 5), (5, 6), (4, 5), (1, 6)]),
            Graph([(1, 2), (1, 4), (2, 3), (3, 5), (5, 6), (4, 5), (1, 6), (2, 6)]),
            Graph([(1, 2), (1, 4), (2, 3), (3, 5), (5, 6), (4, 5), (1, 6), (4, 6)]),
            Graph([(1, 2), (1, 4), (2, 3), (3, 5), (5, 6), (4, 5), (1, 6), (2, 6), (4, 6)]),
            Graph([(1, 2), (1, 4), (2, 3), (3, 5), (5, 6), (3, 4), (2, 6), (1, 6), (4, 5)]),
            Graph([(1, 2), (1, 4), (2, 3), (3, 5), (5, 6), (3, 4), (2, 6), (1, 6), (4, 5), (4, 6)]),
            Graph([(1, 2), (1, 4), (2, 3), (3, 5), (5, 6), (3, 4), (2, 6), (4, 5), (4, 6)]),
            Graph([(1, 2), (1, 4), (2, 3), (3, 5), (5, 6), (3, 4), (2, 6), (4, 6)]),
        ]

        return not any(one_skeleton.subgraph_search(g) is not None for g in obstruction_graphs)

    def cohomology_ring(self, base_ring=QQ):
        r"""
        Return the unreduced cohomology of ``self`` with coefficients in
        ``base_ring``.

        This offers additional information about the cohomology,
        and is intended to be used only when in need of specific
        cohomology operations, such as the cup product, which
        this does offer.

        INPUT:

        - ``base_ring`` -- commutative ring (default: ``QQ``); must be
          ``ZZ`` or a field

        The basis elements in dimension ``dim`` are named ``'h^{dim,i}'``
        where `i` ranges between 0 and `r-1`, where `r` is the rank of
        the cohomology group.

        .. SEEALSO::

            For more information on what this offers, see
            :class:`.moment_angle_complex.CohomologyRing`.

        EXAMPLES::

            sage: Z = MomentAngleComplex([[1,2], [2,3], [3,4], [4,5], [5,1]])
            sage: H = Z.cohomology_ring()
            sage: sorted(H.basis())
            [h^{0,0},
             h^{3,0},
             h^{3,1},
             h^{3,2},
             h^{3,3},
             h^{3,4},
             h^{4,0},
             h^{4,1},
             h^{4,2},
             h^{4,3},
             h^{4,4},
             h^{7,0}]

        Notice (by looking at the dimension) that this does
        indeed coincide with the reduced cohomology::

            sage: Z.cohomology(reduced=False, base_ring=QQ)
            {0: Vector space of dimension 1 over Rational Field,
             1: Vector space of dimension 0 over Rational Field,
             2: Vector space of dimension 0 over Rational Field,
             3: Vector space of dimension 5 over Rational Field,
             4: Vector space of dimension 5 over Rational Field,
             5: Vector space of dimension 0 over Rational Field,
             6: Vector space of dimension 0 over Rational Field,
             7: Vector space of dimension 1 over Rational Field}
            sage: a = H.basis()[3, 0]; a
            h^{3,0}
            sage: b = H.basis()[4, 4]; b
            h^{4,4}
            sage: a.cup_product(b)
            -h^{7,0}
            sage: a * b  # alternative notation
            -h^{7,0}
            sage: a * a
            0

        We can lift cohomology classes to their cocycle
        represetnatives and also acquire their originating
        subcomplexes::

            sage: a.to_cycle()
            \chi_(1,)
            sage: a.get_simplicial_complex()
            Simplicial complex with vertex set (1, 3) and facets {(1,), (3,)}
        """

        return CohomologyRing(base_ring=base_ring, moment_angle_complex=self)


class CohomologyRing(CombinatorialFreeModule):
    r"""
    Cohomology ring of a moment-angle complex.

    Here we don't explicitly compute the cohomology ring
    of a moment-angle complex as a topological space, but
    we use the following result to compute its cohomology
    ring (Theorem 4.5.8 of [BP2014]_):

    .. MATH::

        H^*(\mathcal{Z}_\mathcal{K}) \cong \bigoplus_{J \subseteq [m]}
        \widetilde{H}^*(\mathcal{K}_J),

    where `\widetilde{H}^{-1}(\mathcal{K}_\emptyset)` is the
    base ring over which we compute all of the cohomologies.

    .. NOTE::

        The given isomorphism does not rely on the default
        ring structure given by the direct sum of cohomologies,
        but rather it is a vector space isomorphism, where the ring
        structure (multiplication on the right-hand side) is defined
        in a certain way. See
        :meth:`.moment_angle_complex.CohomologyRing.Element.cup_product`
        for more information.

    .. NOTE::

        This is not intended to be created directly by the user, but instead
        via the :meth:`.moment_angle_complex.MomentAngleComplex.cohomology_ring`
        method.

    INPUT:

    - ``base_ring`` -- the base_ring over which we compute the cohomology
    - ``moment_angle_complex`` -- the moment-angle complex whose homology
      we are computing

    EXAMPLES::

        sage: Z = MomentAngleComplex([[1,2,3], [1,2,4], [3,5], [4,5]])
        sage: H = Z.cohomology_ring(); H
        Cohomology module of Moment-angle complex of Simplicial complex with
        vertex set (1, 2, 3, 4, 5) and facets {(3, 5), (4, 5), (1, 2, 3),
        (1, 2, 4)} over Rational Field
        sage: a = H.an_element(); a
        2*h^{0,0} + 2*h^{3,0} + 3*h^{3,1}
        sage: b = H.basis()[3,2]

    We can compute the cup product::

        sage: a.cup_product(b)
        2*h^{3,2} - 5*h^{6,0}
        sage: a * a
        4*h^{0,0} + 8*h^{3,0} + 12*h^{3,1}

        sage: RP2 = simplicial_complexes.RealProjectivePlane()
        sage: Z = MomentAngleComplex(RP2)
        sage: H = Z.cohomology_ring(GF(2))
        sage: x = H.basis(5)[5, 0]
        sage: y = H.basis(5)[5, 1]
        sage: x * y
        0
    """
    def __init__(self, base_ring, moment_angle_complex):
        """
        Initialize ``self``.

        TESTS::

            sage: S2 = simplicial_complexes.Sphere(2)
            sage: Z = MomentAngleComplex(S2)
            sage: H = Z.cohomology_ring(GF(5))
            sage: TestSuite(H).run()
            sage: H = Z.cohomology_ring(ZZ)
            sage: TestSuite(H).run()
        """
        self._complex = moment_angle_complex
        self._base_ring = base_ring

        vertices = moment_angle_complex._simplicial_complex.vertices()
        n = len(vertices)
        self._graded_indices = {}

        # Will be used for storing information about the subcomplexes
        # from which we compute the cohomology
        self._gens = {}
        indices = []
        for deg in range(moment_angle_complex.dimension() + 1):
            num_of_gens = 0
            self._gens[deg] = []
            for i in range(n+1):
                for x in combinations(vertices, i):
                    S = moment_angle_complex._simplicial_complex.generated_subcomplex(x, is_mutable=False)
                    # Because of the empty combination
                    if len(S.vertices()) > 0 and isinstance(S.cohomology(deg-i-1, base_ring, generators=True), list):
                        chmlgy = S.cohomology(deg-i-1, base_ring, generators=True)
                        for y in chmlgy:
                            self._gens[deg].append((set(x), deg-i-1, y))
                        num_of_gens += len(chmlgy)
                    elif len(S.vertices()) == 0 and deg == 0:
                        num_of_gens = 1

            indices.extend([(deg, k) for k in range(num_of_gens)])
            self._graded_indices[deg] = range(num_of_gens)

        cat = Algebras(base_ring).WithBasis().Graded().FiniteDimensional()
        CombinatorialFreeModule.__init__(self, base_ring, indices, category=cat)

    def basis(self, d=None):
        """
        Return (the degree ``d`` homogeneous component of) the basis
        of ``self``.

        INPUT:

        - ``d`` -- (optional) the degree

        EXAMPLES::

            sage: S2 = simplicial_complexes.Sphere(2)
            sage: Z = MomentAngleComplex(S2)
            sage: H = Z.cohomology_ring()
            sage: H.basis()
            Finite family {(0, 0): h^{0,0}, (7, 0): h^{7,0}}
            sage: H.basis(5)
            Finite family {}
            sage: H.basis(7)
            Finite family {(7, 0): h^{7,0}}
            sage: H.basis()[7, 0]
            h^{7,0}
            sage: H.basis(8)
            Finite family {}
        """
        if d is None:
            return Family(self._indices, self.monomial)
        else:
            indices = [(d, i) for i in self._graded_indices.get(d, [])]
            return Family(indices, self.monomial)

    def degree_on_basis(self, i):
        """
        Return the degree of the basis element indexed by ``i``.

        EXAMPLES::

            sage: Z = MomentAngleComplex([[1,2], [3,4,5], [1,5], [2,3]])
            sage: H = Z.cohomology_ring(GF(5))
            sage: H.degree_on_basis((4, 0))
            4
        """
        return i[0]

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: MomentAngleComplex([[1,2,3], [3,4]]).cohomology_ring()
            Cohomology module of Moment-angle complex of Simplicial complex
            with vertex set (1, 2, 3, 4) and facets {(3, 4), (1, 2, 3)}
            over Rational Field
        """
        return "Cohomology module of {} over {}".format(self._complex, self.base_ring())

    def _repr_term(self, i):
        """
        Return ``'h^{i[0],i[1]}'``, for the basis element indexed by ``i``.

        EXAMPLES::

            sage: Z = MomentAngleComplex([[1,2], [3,4,5], [1,5], [2,3]])
            sage: Z = MomentAngleComplex([[1,2], [3,4,5], [1,4], [2,3,5]])
            sage: H = Z.cohomology_ring(GF(7))
            sage: H.basis()[3, 2]
            h^{3,2}
        """
        return 'h^{{{},{}}}'.format(i[0], i[1])

    _latex_term = _repr_term

    def one(self):
        """
        Return the multiplicative identity element.

        EXAMPLES::

            sage: Z = MomentAngleComplex([[1,2], [3,4], [2,3,5], [1,5], [4,5]])
            sage: H = Z.cohomology_ring()
            sage: H.one()
            h^{0,0}
            sage: all(H.one() * x == x == x * H.one() for x in H.basis())
            True
        """
        one = self._base_ring.one()
        d = {(0, i): one for i in self._graded_indices[0]}
        return self._from_dict(d, remove_zeros=False)

    def complex(self):
        """
        Return the moment-angle complex associated with ``self``.

        EXAMPLES::

            sage: S2 = simplicial_complexes.Sphere(2)
            sage: H = MomentAngleComplex(S2).cohomology_ring()
            sage: H.complex()
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2, 3) and facets {(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)}
        """
        return self._complex

    @cached_method
    def _to_cycle_on_basis(self, i):
        r"""
        Return the cocycle representative of the basis element
        indexed by ``i``.

        .. SEEALSO::

            :meth:`Element.to_cycle`, :meth:`Element.get_simplicial_complex`.

        EXAMPLES::

            sage: Z = MomentAngleComplex([[1,2], [3,4,5], [1,5], [2,3]])
            sage: H = Z.cohomology_ring()
            sage: H._to_cycle_on_basis((3, 2))
            \chi_(2,)
            sage: H._to_cycle_on_basis((4, 2))
            \chi_(2,) + \chi_(4,)
        """
        # The multiplicative identity must be special-cased here
        if i == (0, 0):
            subcomplex = self._complex.simplicial_complex().generated_subcomplex(set(), is_mutable=False)
            cochains = subcomplex.n_chains(-1, base_ring=self._base_ring, cochains=True)
            return cochains.basis().first()

        subcomplex = self._complex.simplicial_complex().generated_subcomplex(self._gens[i[0]][i[1]][0], is_mutable=False)
        cochains = subcomplex.n_chains(self._gens[i[0]][i[1]][1], base_ring=self._base_ring, cochains=True)
        cochain = self._gens[i[0]][i[1]][2][1]
        return cochains.from_vector(cochain.to_vector())

    @cached_method
    def product_on_basis(self, li, ri):
        """
        Return the cup product of the basis elements indexed by
        ``li`` and ``ri`` in ``self``.

        INPUT:

        - ``li``, ``ri`` -- index of a cohomology class

        .. SEEALSO::

            See :meth:`CohomologyRing.Element.cup_product`
            documentation, which describes the algorithm.

        EXAMPLES::

            sage: Z = MomentAngleComplex([[0,1], [1,2], [2,3], [3,0]])
            sage: H = Z.cohomology_ring()
            sage: a, b = H.basis(3)
            sage: b.cup_product(a)
            h^{6,0}
            sage: a.cup_product(b)
            -h^{6,0}

            sage: Z = MomentAngleComplex([[1,2], [1,4], [2,3],
            ....:                         [3,5], [5,6], [4,5], [1,6]])
            sage: H = Z.cohomology_ring(GF(2))
            sage: one = H.one()
            sage: a = H.basis()[4, 3]
            sage: b = H.basis()[4, 6]
            sage: a.cup_product(b)
            h^{8,1}
            sage: b.cup_product(a)
            h^{8,1}
            sage: one.cup_product(a) == a.cup_product(one)
            True
        """
        # First we check whether either one of the elements
        # represents a unit of this cohomology ring;
        # This is necessary because the multiplicative identity
        # is the only element whose corresponding cochain has degree -1
        if li == (0, 0):
            return self.basis()[ri]
        elif ri == (0, 0):
            return self.basis()[li]

        from sage.topology.simplicial_complex import Simplex
        from sage.homology.chain_homotopy import ChainContraction

        left = self._gens[li[0]][li[1]]
        right = self._gens[ri[0]][ri[1]]

        # Sets of vertices of subcomplexes from which
        # these cohomology classes originate
        set_left = left[0]
        set_right = right[0]
        # Because of the definition of multiplication
        if not set_left.isdisjoint(set_right):
            return self.zero()

        # We extract the cocycle monomials in
        # order to loop over them
        left_cocycles = left[2][1].monomials()
        right_cocycles = right[2][1].monomials()
        res = self.zero()
        for left_cocycle in left_cocycles:
            for right_cocycle in right_cocycles:
                # We acquire the set of vertices of cocyles,
                # because of the union
                left_cocycle_vertices = left_cocycle.leading_support().set()
                right_cocycle_vertices = right_cocycle.leading_support().set()
                union = set_left.union(set_right)
                res_union = left_cocycle_vertices.union(right_cocycle_vertices)
                subcomplex_union = self._complex._simplicial_complex.generated_subcomplex(union, is_mutable=False)

                # Because we special-cased the multiplicative identity
                # at the beginning of the method, we know that both
                # left[1] and right[1] are at least 0 (they represent
                # the degree of the corresponding cochains). Therefore,
                # deg is at least 1.
                deg = left[1] + right[1] + 1
                cochains_basis = subcomplex_union.n_chains(deg, cochains=True).basis()
                # If the resulting cochain exists, then we
                # add it to the result; otherwise, we add 0
                if cochains_basis.has_key(Simplex(res_union)):
                    res_part = self.zero()
                    res_cochain = cochains_basis[Simplex(res_union)]
                    # The algebraic_topological_model() works with unreduced cohomology,
                    # but we require reduced cohomology here. This only causes problems
                    # in dimension 0, but because of excluding the unit element at
                    # the beginning of the method, we know that deg is at least 1
                    # See comments at the beginning of the method and when defining deg
                    phi, _ = subcomplex_union.algebraic_topological_model(self._base_ring)
                    coeff_vec = phi.dual().pi().in_degree(deg) * res_cochain.to_vector()
                    # We compute the sign of the result
                    for i in range(len(coeff_vec)):
                        res_part += coeff_vec[i] * self.basis()[li[0]+ri[0], i]
                    zeta = 1
                    for k in set_left.difference(left_cocycle_vertices):
                        zeta *= eps({k}, set_right.union({k}).difference(right_cocycle_vertices))
                    epsilon = (eps(left_cocycle_vertices, set_left)
                               * eps(right_cocycle_vertices, set_right)
                               * eps(res_union, union)
                               * zeta)
                    if epsilon == -1:
                        res_part = -res_part
                    res += res_part

        return res

    def cup_length(self):
        """
        Return the cup length of ``self``.

        The cup length of a cohomology ring is defined as the maximal
        number of elements whose product is non-trivial.

        EXAMPLES::

            sage: Z = MomentAngleComplex([[0], [1], [2]])
            sage: H = Z.cohomology_ring(GF(2))
            sage: H.cup_length()
            1
            sage: Z = MomentAngleComplex([[1,2], [2,3], [3,4], [4,5], [5,1]])
            sage: H = Z.cohomology_ring()
            sage: H.cup_length()
            2
        """
        # It suffices to check this for base elements
        elements = set(self.basis())
        elements.remove(self.one())
        # We create a dictionary, in which the keys are
        # base elements, and values are lists, which
        # represent the elements whose product with
        # the key gives us zero
        memo = {}
        for base_element in elements:
            memo[base_element] = set()
            for x in elements:
                if base_element * x == self.zero():
                    memo[base_element].add(x)

        max_length = 1
        # The multiplication is commutative up to a sign,
        # so we will not worry about the order in which
        # we multiply the elements
        for base_element in elements:
            prod = base_element
            length = 1
            # We will avoid the elements which give us zero
            avoid = memo[base_element]
            for x in elements:
                if x not in avoid and prod * x != self.zero():
                    prod = prod * x
                    length = length + 1
                    # We further restrict the possible elements
                    avoid = avoid.union(memo[x])

            if length > max_length:
                max_length = length

        return max_length

    class Element(CombinatorialFreeModule.Element):
        def to_cycle(self):
            r"""
            Return the cocycle representative of ``self``.

            The cohomology class gets lifted to its cococyle
            representative in the cochain complex of the
            corresponding simplicial complex given in the main
            isomorphism.

            .. NOTE::

                The cocycle is practically useless without
                knowing from which subcomplex it origintes.
                We can use the :meth:`get_simplicial_complex`
                method to retrieve that information.

            EXAMPLES::

                sage: Z = MomentAngleComplex([[0,1], [1,2], [2,0],
                ....:                         [1,3], [3,4], [4,1]])
                sage: H = Z.cohomology_ring()
                sage: a = H.basis()[5, 1]
                sage: a.to_cycle()
                \chi_(3, 4)
                sage: b = H.basis()[6, 3]
                sage: b.to_cycle()
                \chi_(3, 4)
                sage: H.one().to_cycle()
                \chi_()

            We can obtain the originating simplicial complex
            by using :meth:`get_simplicial_complex`::

                sage: a.get_simplicial_complex()
                Simplicial complex with vertex set (1, 3, 4) and
                facets {(1, 3), (1, 4), (3, 4)}
                sage: H.one().get_simplicial_complex()
                Simplicial complex with vertex set () and facets {()}
            """
            if not self.is_homogeneous():
                raise ValueError("only defined for homogeneous elements")
            return sum(c * self.parent()._to_cycle_on_basis(i) for i, c in self)

        def get_simplicial_complex(self, is_mutable=True):
            r"""
            Return the simplicial complex from which ``self`` originates.

            EXAMPLES::

                sage: Z = MomentAngleComplex([[0,1], [1,2], [2,0], [1,3,4]])
                sage: H = Z.cohomology_ring()
                sage: a = H.basis()[3, 3]
                sage: a.get_simplicial_complex()
                Simplicial complex with vertex set (2, 4) and facets {(2,), (4,)}
                sage: b = H.basis()[7, 0]
                sage: b.get_simplicial_complex(is_mutable=False).is_mutable()
                False

            Because `\widetilde{H}^{-1}(\mathcal{K}_\emptyset)` is
            isomorphic to the base ring, the multiplicative identity
            in the cohomology of a moment-angle complex originates from
            here, and its simplicial complex is the empty subcomplex
            of `\mathcal{K}`::

                sage: H.one().get_simplicial_complex()
                Simplicial complex with vertex set () and facets {()}
            """
            if not self.is_homogeneous():
                raise ValueError("only defined for homogeneous elements")
            if self.is_one():
                return self.parent()._complex.simplicial_complex().generated_subcomplex(set(), is_mutable=is_mutable)
            vertex_set = self.parent()._gens[self.leading_support()[0]][self.leading_support()[1]][0]
            return self.parent()._complex.simplicial_complex().generated_subcomplex(vertex_set, is_mutable=is_mutable)

        def cup_product(self, other):
            r"""
            Return the cup product of ``self`` and ``other``.

            We define the cup product on cochains in cochain complexes
            of full subcomplexes of the simplicial complex associated
            with the moment-angle complex of ``self`` (this is
            possible because of the main cohomology ring isomorphism).
            Given two cochains, `\chi_L` and `\chi_M`, from cochain
            complexes of full subcomplexes `\mathcal{K}_I` and `\mathcal{K}_J`,
            respectively, we define

            .. MATH::

                \chi_L \otimes \chi_M =
                \begin{cases}
                    c_{L \cup M} \chi_{L \cup M}, &\text{if } I\cap J =
                    \emptyset,\\
                    0, &\text{otherwise}
                \end{cases}

            where

            .. MATH::

                c_{L \cup M} = \epsilon(L, I) \epsilon(M, J) \zeta
                \epsilon(L\cup M, I\cup J)

            and

            .. MATH::

                \zeta = \prod_{k \in I \setminus L} \epsilon(k, k \cup J
                \setminus M).

            ALGORITHM:

            This algorithm is adopted from [Lin2019]_, p. 16.

            First we lift the cohomology classes to their cocycle
            representatives (which are generators of cohomology groups of all
            full subcomplexes of the simplicial complex associated with the
            moment-angle complex; these are stored during initialization),
            as well as acquire their corresponding simplicial complexes.
            We then multiply them, as described above and extract the
            cohomology class using chain contractions of the appropriate
            simplicial complex.

            .. SEEALSO::

                For more information on the `\epsilon` function,
                see :meth:`sage.topology.moment_angle_complex.eps`.

            EXAMPLES::

                sage: Z = MomentAngleComplex([[1,2], [2,3], [3,4], [4,1]])
                sage: H = Z.cohomology_ring()
                sage: a, b = H.basis(3)
                sage: a.get_simplicial_complex()
                Simplicial complex with vertex set (1, 3) and facets {(1,), (3,)}
                sage: b.get_simplicial_complex()
                Simplicial complex with vertex set (2, 4) and facets {(2,), (4,)}

            We see that the vertex sets are indeed disjoint, so
            we can expect that the product of ``a`` and ``b`` is
            indeed non-trivial in this case::

                sage: c = a.cup_product(b); c
                -h^{6,0}
                sage: b * a
                h^{6,0}
                sage: c.get_simplicial_complex()
                Simplicial complex with vertex set (1, 2, 3, 4) and
                facets {(1, 2), (1, 4), (2, 3), (3, 4)}

                sage: Z = MomentAngleComplex([[1, 2], [1, 4], [2, 3], [3, 5],
                ....:                         [5, 6], [3, 4], [2, 6], [4, 6]])
                sage: H = Z.cohomology_ring(GF(2))
                sage: a = H.basis()[4, 0]; a.get_simplicial_complex()
                Simplicial complex with vertex set (1, 2, 5) and facets {(5,), (1, 2)}
                sage: b = H.basis()[4, 1]; b.get_simplicial_complex()
                Simplicial complex with vertex set (1, 3, 5) and facets {(1,), (3, 5)}

            Here the vertex sets are not disjoint, so the
            cup product is trivial. Because of this definition,
            multiplying any cohomology class with itself yields 0::

                sage: a * b
                0
                sage: all(a ** n == H.zero() for n in range(2, dim(Z)))
                True
            """
            return self * other


# Used for computing coeffeicients when multiplying in cohomology
def eps(subcomplex, simplicial_complex):
    r"""
    Return the coefficient `\epsilon`, used when computing
    the cup product of cohomology classes.

    By definition `\epsilon(j, J) = (-1)^{r-1}`, where `j`
    is the `r`-th element of `J`. For a subset `L\subset J`
    we define

    .. MATH::

        \epsilon(L, J) := \prod_{j \in L} \epsilon(j, J).

    INPUT:

    - ``subcomplex`` -- a set or a simplicial complex,
      represents `L` given above
    - ``simplicial_complex`` -- a set or a simplicial
      complex, represents `J` given above

    REFERENCES:

        Fore more information, see [Lin2019]_, p. 13.

    EXAMPLES::

        sage: from sage.topology.moment_angle_complex import eps
        sage: eps({1,2}, {1,2,3,4})
        -1
        sage: eps({1,3}, {1,2,3,4})
        1
        sage: eps({2}, SimplicialComplex([[1,2], [2,3,4], [4,1]]))
        -1
    """
    def _eps(element, simplicial_complex):
        if element not in simplicial_complex._vertex_to_index:
            raise ValueError("{} is not a vertex of this simplicial complex".format(element))
        return (-1) ** simplicial_complex._vertex_to_index[element]

    if not isinstance(subcomplex, SimplicialComplex):
        # This may completely change the structure,
        # but the ordering of the vertices will
        # remain the same, and that's the only things
        # that matters
        subcomplex = SimplicialComplex([subcomplex])
    if not isinstance(simplicial_complex, SimplicialComplex):
        simplicial_complex = SimplicialComplex([simplicial_complex])
    res = 1
    for element in subcomplex.vertices():
        res *= _eps(element, simplicial_complex)
    return res
