r"""
Open subset of Euclidian space with coordinates

An open subset of Euclidian space with a specific set of coordinates.  This
is the background on which differential forms can be defined.

AUTHORS:

- Joris Vankerschaver (2010-07-25)

EXAMPLES::

    sage: x, y, z = var('x, y, z')
    sage: S = CoordinatePatch((x, y, z)); S
    doctest:...: DeprecationWarning: Use Manifold instead.
    See http://trac.sagemath.org/24444 for details.
    Open subset of R^3 with coordinates x, y, z

::

    sage: u, v = var('u, v')
    sage: S = CoordinatePatch((u, v)); S
    Open subset of R^2 with coordinates u, v

.. TODO::

    Add functionality for metric tensors

"""

# ****************************************************************************
#    Copyright (C) 2010 Joris Vankerschaver <joris.vankerschaver@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.structure.parent import Parent


class CoordinatePatch(Parent):
    """
        Construct a coordinate patch, i.e. an open subset of
        Euclidian space with a given set of coordinates.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: S = CoordinatePatch((x, y, z)); S
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z

            sage: u, v = var('u, v')
            sage: T = CoordinatePatch((u, v)); T
            Open subset of R^2 with coordinates u, v
            sage: loads(T.dumps()) == T
            True

        In a future release, it will be possible to specify a
        metric tensor on a coordinate patch.  For now, providing
        any kind of metric raises an exception::

            sage: x, y, z = var('x, y, z')
            sage: m = matrix(SR, 3)
            sage: S = CoordinatePatch((x, y, z), metric=m)
            Traceback (most recent call last):
            ...
            NotImplementedError: Metric geometry not supported yet.

    """
    def __init__(self, coordinates, metric = None):
        """
        An open subset of Euclidian space with a specific set of
        coordinates. See ``CoordinatePatch`` for details.

        INPUT:

        - ``coordinates`` -- a set of symbolic variables that serve
          as coordinates on this space.

        - ``metric`` (default: ``None``) -- a metric tensor on this
          coordinate patch.  Providing anything other than ``None``
          is currently not defined.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: S = CoordinatePatch((x, y, z)); S
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
        """
        from sage.symbolic.ring import is_SymbolicVariable
        from sage.misc.superseded import deprecation
        deprecation(24444, 'Use Manifold instead.')

        if not all(is_SymbolicVariable(c) for c in coordinates):
            raise TypeError("%s is not a valid vector of coordinates." %
                            coordinates)

        self._coordinates = tuple(coordinates)

        if metric is not None:
            raise NotImplementedError("Metric geometry not supported yet.")

    def __eq__(self, other):
        """
        Return equality if and only if other has the same coordinates
        as self, in the same order.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: S = CoordinatePatch((x, y, z)); S
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: u, v = var('u, v')
            sage: T = CoordinatePatch((u, v)); T
            Open subset of R^2 with coordinates u, v
            sage: U = CoordinatePatch((x, y, z)); U
            Open subset of R^3 with coordinates x, y, z
            sage: U == S
            True
            sage: U == U
            True
            sage: U == T
            False

        Note that the order of the coordinates matters::

            sage: x, y, z = var('x, y, z')
            sage: S = CoordinatePatch((x, y, z)); S
            Open subset of R^3 with coordinates x, y, z
            sage: T = CoordinatePatch((x, z, y)); T
            Open subset of R^3 with coordinates x, z, y
            sage: S == T
            False
        """

        return str(self._coordinates) == str(other._coordinates)


    def __ne__(self, other):
        """
        Test whether two coordinate patches are not equal.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: S = CoordinatePatch((x, y, z)); S
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: u, v = var('u, v')
            sage: T = CoordinatePatch((u, v)); T
            Open subset of R^2 with coordinates u, v
            sage: S != T
            True
        """
        return not self == other

    def coordinates(self):
        """
        Return coordinates on this coordinate patch.

        OUTPUT:

        - list -- a list of coordinates on this space.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: S = CoordinatePatch((x, y, z)); S
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: S.coordinates()
            (x, y, z)

        """
        return self._coordinates


    def coordinate(self, i=0):
        """
        Return the `i^{th}` coordinate on ``self``

        INPUT:

        - ``i`` - integer (optional, default 0)


        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: S = CoordinatePatch((x, y, z)); S
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: S.coordinate(0)
            x
            sage: S.coordinate(1)
            y
            sage: S.coordinate(2)
            z

        """
        return self._coordinates[i]


    def dim(self):
        """
        Return the dimension of this coordinate patch, i.e. the dimension
        of the Euclidian space of which this coordinate patch is an open
        subset.

        EXAMPLES::

            sage: a, b, c, d, e = var('a, b, c, d, e')
            sage: U = CoordinatePatch((a, b, c, d, e)); U
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^5 with coordinates a, b, c, d, e
            sage: U.dim()
            5
        """
        return len(self._coordinates)


    def _repr_(self):
        r"""
        Return string representation of this coordinate patch.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: S = CoordinatePatch((x, y, z)); S
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: S._repr_()
            'Open subset of R^3 with coordinates x, y, z'
            sage: S.rename('coordinate patch'); S
            coordinate patch
            sage: S.rename(); S
            Open subset of R^3 with coordinates x, y, z
        """
        return r"Open subset of R^%s with coordinates %s" % \
            (self.dim(), ', '.join(x._latex_() for x in self._coordinates))

    def _latex_(self):
        r"""
        Return latex representation of this coordinate patch.

        EXAMPLES::

            sage: x, y, z = var('x, y, z')
            sage: S = CoordinatePatch((x, y, z)); S
            doctest:...: DeprecationWarning: Use Manifold instead.
            See http://trac.sagemath.org/24444 for details.
            Open subset of R^3 with coordinates x, y, z
            sage: latex(S)
            \mathbb{\RR}^3
            sage: latex(S) == S._latex_()
            True
        """
        return "\\mathbb{\\RR}^%s" % self.dim()
