r"""
Moment angle complexes

AUTHORS:

- Ognjen Petrov (2023-06-25)

This module implements the basic structure of moment-angle complexes.
Given a simplicial complex `K`, with a set of vertices
`V = \{v_1, v_2, \dotso, v_n\}`, a moment-angle complex over `K` is a
topological space `Z`, which is a disjoint union of `X_{\sigma}`, where
`\sigma \in K`, and `X_{\sigma} = Y_{v_1} \times Y_{v_2} \times \dotso \times Y_{v_n}`
and `Y_{v_i}` is a 2-disk (a 2-simplex) if `v_i \in \sigma`, or a 1-sphere otherwise.

.. MATH::
   :nowrap:

    Y_{v_i} =
    \begin{cases}
        D^2, &v_i \in \sigma\\
        S^1, &v_i \notin \sigma
    \end{cases}

They are one of the main topics of resarch in fields such as algebraic and
toric topology, as well as combinatorial algebraic geometry.

Here we view them as simplicial complexes and try to compute mostly
things which would not require computing the moment-angle complex itself,
but rather work with the corresponding simplicial complex.

.. NOTE::

    One of the most useful properties will be the
    :meth:`bigraded Betti numbers<sage.topology.simplicial_complex.bigraded_betti_numbers>`,
    and the underlying theorem which makes this possible is Hochter's formula, which
    can be found on page 104 of :arxiv:`Toric topoloogy<1210.2368>`.

EXAMPLES::

<Lots and lots of examples>

"""

# ****************************************************************************
#       Copyright (C) 2013 Ognjen Petrov <ognjenpetrov@yahoo.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.homology.homology_group import HomologyGroup
from sage.rings.integer_ring import ZZ
from sage.structure.sage_object import SageObject
from .simplicial_complex import SimplicialComplex
from .simplicial_complex_examples import Sphere, Simplex
from itertools import combinations

#TODO's:
# - Documentation(examples and tests)
# - latex, represenation
# - perhaps make _components a dictionary, indexed by \sigma?
# - compute up to homotopy?
# - and a lot more ...

class MomentAngleComplex(SageObject): # should this inherit SimplicialComplex?
    """
    Define a moment-angle complex.

    :param simplicial_complex: the corresponding simplicial complex
    :type simplicial_complex: a simplicial complex
    :param construct: see below
    :type construct: boolean; optional, default ``False``
    :return: the associated moment-angle complex

    ``simplicial_complex`` must be an intance of
    :class:`~sage.topology.simplicial_complex.SimplicialComplex`, which
    represents the simplicial complex whose moment-angle complex we
    wish to create.

    If ``construct`` is ``True``, we also explicitly compute the
    moment-angle complexes, with the construction described above.

    .. WARNING::

        The construction can be very slow, it is not reccomended unless
        the corresponding simplicial complex has less than 5 vertices.

    EXAMPLES::

    <Lots and lots of examples>
    """

    def __init__(self,
                 simplicial_complex,
                 construct=False):
        """
        Define a moment-angle complex.  See :class:`MomentAngleComplex`
        for full documentation.

        EXAMPLES::

        """
        if not isinstance(simplicial_complex, SimplicialComplex):
            raise ValueError("simplicial_complex must be a simplicial complex")

        self._simplicial_complex = simplicial_complex
        self._moment_angle_complex = None

        vertices = simplicial_complex.vertices()
        k = len(vertices)
        self._components = []
        self._symbolic_components = []

        for face in simplicial_complex.faces().values():
            for subcomplex in face:
                Y = []
                Ys = []
                for j in vertices:
                    if j in subcomplex:
                        Y.append(Simplex(2))
                        Ys.append("D^2")
                    else:
                        Y.append(Sphere(1))
                        Ys.append("S^1")

                self._components.append(Y)
                self._symbolic_components.append(Ys)

        self._constructed = False
        if construct:
            self.construct()

    def __eq__(self, other):
        """
        Return ``True`` iff ``self`` is the same moment-angle complex
        as ``other``.

        We consider two moment-angle complexes to be equal if their
        corresponding simplicial complexes are equal.

        :param other: the other moment-angle complex

        EXAMPLES::

        """
        return isinstance(other, MomentAngleComplex) and self._simplicial_complex.__eq__(other._simplicial_complex)

    def __ne__(self, other):
        """
        Return ``True`` iff ``self`` is not equal to ``other``.

        :param other: the other moment-angle complex

        EXAMPLES::

        """
        return not self.__eq__(other)

    def construct(self):
        """
        Construct the moment-angle complex as a simplicial complex.

        EXAMPLES::

        """
        if self._constructed:
            return

        self._moment_angle_complex = SimplicialComplex()
        for component in self._components:
            x = component[0]
            for j in range(1, len(component)-1):
                x = x.product(component[j])

            self._moment_angle_complex = self._moment_angle_complex.disjoint_union(x)

        self._constructed = True

    def vertices(self):
        """
        Return vertices of the moment-angle complex.

        EXAMPLES::

        """
        if not self._constructed:
            raise ValueError("the moment-angle complex is not constructed")

        return self._moment_angle_complex.vertices()

    def homology(self, dim=None, base_ring=ZZ, subcomplex=None,
                 generators=False, cohomology=False, algorithm='pari',
                 verbose=False, reduced=True, **kwds):
        """
        The (reduced) homology of this moment-angle complex.


        .. SEEALSO::
            :meth:`~sage.topology.cell_complex.homology`

        EXAMPLES::

        """
        if not self._constructed:
            raise ValueError("the moment-angle complex is not constructed")

        return self._moment_angle_complex.homology(dim, base_ring, subcomplex,
                 generators, cohomology, algorithm,
                 verbose, reduced, **kwds)
