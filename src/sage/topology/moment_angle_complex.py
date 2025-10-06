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
from itertools import combinations

from sage.categories.fields import Fields
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import lazy_import
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.topology import simplicial_complex_catalog as simplicial_complexes
from sage.topology.cubical_complex import CubicalComplex, cubical_complexes
from sage.topology.simplicial_complex import SimplicialComplex, copy

lazy_import('sage.homology.homology_group', 'HomologyGroup')


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

        sage: M.homology()                                                              # needs sage.modules
        {0: 0, 1: 0, 2: 0, 3: Z}
        sage: Z.dimension()
        6
        sage: Z.homology()                                                              # needs sage.modules
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

    def __init__(self, simplicial_complex) -> None:
        """
        Initialize ``self``.

        TESTS::

            sage: Z = MomentAngleComplex([[0,1,2], [1,2,3], [0, 3]])
            sage: TestSuite(Z).run()
        """
        # The underlying simplicial complex
        self._simplicial_complex = copy(simplicial_complex)
        vertices = self._simplicial_complex.vertices()

        disk = simplicial_complexes.Simplex(2)
        circle = simplicial_complexes.Sphere(1)

        # A dictionary of components indexed by facets
        self._components = {facet: [disk if j in facet else circle
                                    for j in vertices]
                            for facet in self._simplicial_complex.maximal_faces()}

    @lazy_attribute
    def _moment_angle_complex(self):
        """
        Create the moment-angle complex as a cubical complex.

        If this lazy attribute is accessed, we explicitly compute
        the moment-angle complex, viewed as a cubical complex.

        .. WARNING::

            The construction can be very slow, it is not recommended unless
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
        cube = cubical_complexes.Cube(2)
        sphere = cubical_complexes.Sphere(1)

        moment_angle_complex = CubicalComplex()
        for component in self._components.values():
            x = cube if component[0] == simplicial_complexes.Simplex(2) else sphere
            for j in range(1, len(component)):
                y = cube if component[j] == simplicial_complexes.Simplex(2) else sphere
                x = x.product(y)
            moment_angle_complex = _cubical_complex_union(moment_angle_complex, x)

        return moment_angle_complex

    def _repr_(self) -> str:
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

            The construction can be very slow, it is not recommended unless
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

    def components(self) -> dict:
        r"""
        Return the dictionary of components of ``self``, indexed by facets
        of the associated simplicial complex.

        OUTPUT:

        A dictionary, whose values are lists, representing spheres
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
            sage: Z.cohomology()                                                        # needs sage.modules
            {0: 0, 1: 0, 2: 0, 3: Z x Z, 4: 0, 5: 0, 6: Z}
            sage: Z.cohomology() == product_of_spheres.cohomology()  # long time        # needs sage.modules
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

            sage: # needs sage.modules
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

        in_field = base_ring in Fields()

        for j in range(n + 1):
            for x in combinations(vertices, j):
                S = self._simplicial_complex.generated_subcomplex(x)
                if in_field:
                    invfac.append(S.homology(i - j - 1, base_ring=base_ring,
                                             cohomology=cohomology, algorithm=algorithm,
                                             verbose=verbose, reduced=True).dimension())
                else:
                    invfac.extend(S.homology(i - j - 1, base_ring=base_ring,
                                             cohomology=cohomology, algorithm=algorithm,
                                             verbose=verbose, reduced=True)._original_invts)

        if in_field:
            return HomologyGroup(sum(invfac), base_ring)

        m = len(invfac)
        return HomologyGroup(m, base_ring, invfac)

    def homology(self, dim=None, base_ring=ZZ, cohomology=False,
                 algorithm='pari', verbose=False, reduced=True) -> dict:
        r"""
        The (reduced) homology of ``self``.

        INPUT:

        - ``dim`` -- integer or a list of integers; represents the
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

            sage: # needs sage.modules
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
            sage: Z.cubical_complex().homology() == Z.homology()                        # needs sage.modules
            True

        Meanwhile, the homology computation used here is quite efficient
        and works well even with significantly larger underlying simplicial
        complexes::

            sage: # needs sage.modules
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
            sage: Z.homology(dim=range(5), reduced=True)
            {0: 0, 1: 0, 2: 0, 3: Z x Z x Z x Z, 4: Z x Z}
            sage: Z.homology(dim=range(5), reduced=False)
            {0: Z, 1: 0, 2: 0, 3: Z x Z x Z x Z, 4: Z x Z}
            sage: all(Z.homology(i,reduced=True) == Z.homology(i,reduced=False)
            ....:     for i in range(1, dim(Z)))
            True
            sage: all(Z.homology(i,reduced=True) == Z.homology(i,reduced=False)
            ....:     for i in range(dim(Z)))
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
            dims = range(self.dimension() + 1)

        return {i: self._homology_group(i, base_ring=base_ring, cohomology=cohomology,
                                        algorithm=algorithm, verbose=verbose, reduced=reduced) for i in dims}

    def cohomology(self, dim=None, base_ring=ZZ, algorithm='pari',
                   verbose=False, reduced=True) -> dict:
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

            sage: # needs sage.modules
            sage: S3 = simplicial_complexes.Sphere(3)
            sage: product_of_spheres = S3.product(S3)
            sage: Z.cohomology()
            {0: 0, 1: 0, 2: 0, 3: Z x Z, 4: 0, 5: 0, 6: Z}
            sage: Z.cohomology() == product_of_spheres.cohomology()  # long time
            True
        """
        return self.homology(dim=dim, cohomology=True, base_ring=base_ring,
                             algorithm=algorithm, verbose=verbose, reduced=reduced)

    def betti(self, dim=None) -> dict:
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

            sage: # needs sage.modules
            sage: Z = MomentAngleComplex([[0,1], [1,2], [2,0], [1,2,3]])
            sage: Z.betti()
            {0: 1, 1: 0, 2: 0, 3: 1, 4: 0, 5: 1, 6: 1, 7: 0}
            sage: Z = MomentAngleComplex([[0,1], [1,2], [2,0], [1,2,3], [3,0]])
            sage: Z.betti(dim=6)
            {6: 2}
        """
        dic = {}
        H = self.homology(dim=dim, base_ring=QQ)
        try:
            for n in H:
                dic[n] = H[n].dimension()
                if n == 0:
                    dic[n] += 1
        except AttributeError:
            return H.dimension()
        else:
            return dic

    def euler_characteristic(self):
        """
        Return the Euler characteristic of ``self``.

        The Euler characteristic is defined as the alternating sum
        of the Betti numbers of ``self``.

        The Euler characteristic of a moment-angle complex is 0
        if the associated simplicial complex is not a simplex.

        EXAMPLES::

            sage: # needs sage.modules
            sage: X = SimplicialComplex([[0,1,2,3,4,5], [0,1,2,3,4,6],
            ....:                        [0,1,2,3,5,7], [0,1,2,3,6,8,9]])
            sage: M = MomentAngleComplex(X)
            sage: M.euler_characteristic()  # long time
            0
            sage: Z = MomentAngleComplex([[0,1,2,3,4]])
            sage: Z.euler_characteristic()
            1
        """
        sc = self.simplicial_complex()
        return (ZZ.one() if sc.dimension() + 1 == len(sc.vertices())
                else ZZ.zero())

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

    def has_trivial_lowest_deg_massey_product(self) -> bool:
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

        return all(one_skeleton.subgraph_search(g) is None for g in obstruction_graphs)
