r"""
Oriented matroid from real hyperplane arrangements

This implements an oriented matroid from real hyperplane arrangements

AUTHORS:

- Aram Dermenjian (2019-07-12): Initial version
"""

# ****************************************************************************
#      Copyright (C) 2019   Aram Dermenjian <aram.dermenjian.math at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matroids.oriented_matroids.covector_oriented_matroid import CovectorOrientedMatroid


class RealHyperplaneArrangementOrientedMatroid(CovectorOrientedMatroid):
    r"""
    An oriented matroid implemented from a real hyperplane arrangement.

    This implements an oriented matroid using hyperplane arrangements.
    Oriented matroids arise from central hyperplane arrangements.

    INPUT:

    - ``data`` -- a :class:`HyperplaneArrangementElement` element.
    - ``groundset`` -- (default: ``None``) the groundset for the data; if not
      provided, we grab the data from the hyperplane arrangement

    EXAMPLES::

        sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
        sage: A = hyperplane_arrangements.braid(3)
        sage: M = OrientedMatroid(A); M
        Hyperplane arrangement oriented matroid of rank 2

        sage: A = hyperplane_arrangements.Catalan(3)
        sage: M = OrientedMatroid(A)
        sage: M.is_valid(certificate=True)
        (False, 'hyperplane arrangement must be central')

        sage: G = Graph({1: [2,4], 2: [3,4]})
        sage: A = hyperplane_arrangements.graphical(G)
        sage: M = OrientedMatroid(A); M
        Hyperplane arrangement oriented matroid of rank 3

    .. SEEALSO::

        - :class:`~sage.oriented_matroids.oriented_matroid.OrientedMatroid`
        - :class:`~sage.oriented_matroids.covector_oriented_matroid.CovectorOrientedMatroid`
        - :class:`sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangementElement`
    """
    def __init__(self, data, groundset=None, category=None):
        """
        Return a ``RealHyperplaneArrangementOrientedMatroid`` object.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: A = hyperplane_arrangements.braid(3)
            sage: M = OrientedMatroid(A); M
            Hyperplane arrangement oriented matroid of rank 2
            sage: TestSuite(M).run()
        """
        self._arrangement = data

        if data and groundset is None:
            groundset = tuple(data.hyperplanes())

        # Set up our covectors after our groundset is made
        faces = [i[0] for i in self._arrangement.closed_faces()]

        CovectorOrientedMatroid.__init__(self, data=faces, groundset=groundset, category=category)

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: A = hyperplane_arrangements.braid(3)
            sage: M = OrientedMatroid(A); M
            Hyperplane arrangement oriented matroid of rank 2
        """
        try:
            rep = "Hyperplane arrangement oriented matroid of rank {}".format(
                self.arrangement().rank())
        except ValueError:
            rep = "Hyperplane arrangement oriented matroid"
        return rep

    def is_valid(self, certificate=False) -> bool | tuple[bool, str]:
        """
        Return whether or not the arrangement is an oriented matroid.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: A = hyperplane_arrangements.braid(3)
            sage: M = OrientedMatroid(A)
            sage: M.is_valid()
            True
        """
        if not self.arrangement().is_central():
            if certificate:
                return (
                    False,
                    "hyperplane arrangement must be central"
                )
            return False

        if certificate:
            return (True, "")
        return True

    def arrangement(self):
        """
        Return the arrangement.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: G = Graph({1: [2,4], 2: [3,4]})
            sage: A = hyperplane_arrangements.graphical(G)
            sage: M = OrientedMatroid(A); M
            Hyperplane arrangement oriented matroid of rank 3
            sage: M.arrangement()
            Arrangement <t1 - t2 | t1 - t3 | t0 - t1 | t0 - t3>
        """
        return self._arrangement

    def deletion(self, hyperplanes):
        """
        Return the hyperplane arrangement oriented matroid with hyperplanes removed.

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: G = Graph({1: [2,4], 2: [3,4], 3: [4], 4:[]})
            sage: A = hyperplane_arrangements.graphical(G)
            sage: H = [A.hyperplanes()[i] for i in range(2, 5)]
            sage: # long time
            sage: M = OrientedMatroid(A); M
            Hyperplane arrangement oriented matroid of rank 3
            sage: M2 = M.deletion(H); M2
            Hyperplane arrangement oriented matroid of rank 2
        """
        A = self.arrangement()
        if isinstance(hyperplanes, list) or isinstance(hyperplanes, tuple):
            for h in hyperplanes:
                A = A.deletion(h)
        else:
            A = A.deletion(hyperplanes)
        return RealHyperplaneArrangementOrientedMatroid(A)