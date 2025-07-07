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

from __future__ import annotations

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
        (False, {'elt': None, 'msg': 'hyperplane arrangement must be central'})

        sage: G = Graph({1: [2,4], 2: [3,4]})
        sage: A = hyperplane_arrangements.graphical(G)
        sage: M = OrientedMatroid(A); M
        Hyperplane arrangement oriented matroid of rank 3

    .. SEEALSO::

        - :class:`~sage.oriented_matroids.oriented_matroid.OrientedMatroid`
        - :class:`~sage.oriented_matroids.covector_oriented_matroid.CovectorOrientedMatroid`
        - :class:`sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangementElement`
    """
    def __init__(self, data, groundset=None):
        """
        Return a ``RealHyperplaneArrangementOrientedMatroid`` object.

        INPUT:

        - ``data`` -- a :class:`HyperplaneArrangementElement` element.
        - ``groundset`` -- (default: ``None``) the groundset for the data; if not
          provided, we grab the data from the hyperplane arrangement

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

        CovectorOrientedMatroid.__init__(self, data=faces, groundset=groundset)

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

    def is_valid(self, certificate=False) -> bool | tuple[bool, dict]:
        """
        Return whether or not the arrangement is an oriented matroid.

        INPUT:

        - ``certificate`` -- (default: ``False``) additional information
          on why the oriented matroid is/not valid

        EXAMPLES::

            sage: from sage.matroids.oriented_matroids.oriented_matroid import OrientedMatroid
            sage: A = hyperplane_arrangements.braid(3)
            sage: M = OrientedMatroid(A)
            sage: M.is_valid()
            True

            sage: A = hyperplane_arrangements.Shi(2)
            sage: M = OrientedMatroid(A)
            sage: M.is_valid(certificate=True)
            (False, {'elt': None, 'msg': 'hyperplane arrangement must be central'})
        """
        if not self.arrangement().is_central():
            if certificate:
                error_info = {
                    'msg': "hyperplane arrangement must be central",
                    'elt': None
                }
                return (
                    False,
                    error_info
                )
            return False

        if certificate:
            return (True, {})
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

    def deletion(self, change_set):
        """
        Return the hyperplane arrangement oriented matroid with hyperplanes removed.

        INPUT:

        - ``change_set`` - hyperplane or list of hyperplanes to take deletion of

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

            sage: A = hyperplane_arrangements.braid(3)
            sage: M = OrientedMatroid(A)
            sage: M.groundset()
            (Hyperplane 0*t0 + t1 - t2 + 0,
            Hyperplane t0 - t1 + 0*t2 + 0,
            Hyperplane t0 + 0*t1 - t2 + 0)
            sage: D = M.deletion(M.groundset()[0]); D
            Hyperplane arrangement oriented matroid of rank 2
            sage: D.groundset()
            (Hyperplane t0 - t1 + 0*t2 + 0, Hyperplane t0 + 0*t1 - t2 + 0)
        """
        A = self.arrangement()
        if isinstance(change_set, list) or isinstance(change_set, tuple):
            for h in change_set:
                A = A.deletion(h)
        else:
            A = A.deletion(change_set)
        return RealHyperplaneArrangementOrientedMatroid(A)
