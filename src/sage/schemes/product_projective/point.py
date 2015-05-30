r"""
Points for products of projective spaces

This class builds on the projective space class and its point and morphism classes.

EXAMPLES:

We construct products projective spaces of various dimensions over the same ring.::

    sage: P1xP1.<x,y, u,v> = ProductProjectiveSpaces(QQ, [1,1])
    sage: P1xP1([2,1, 3,1])
    (2 : 1 , 3 : 1)
"""
#*****************************************************************************
# Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#                    Ben Hutz <bn4941@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
# http://www.gnu.org/licenses/
#*****************************************************************************
from copy import copy

from sage.schemes.generic.morphism import SchemeMorphism
from sage.schemes.generic.morphism import SchemeMorphism_point


class ProductProjectiveSpaces_point_ring(SchemeMorphism_point):
    r"""
    The class of points on products of projective spaces.
    The components are projective space points.

    EXAMPLES::

        sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2,1],QQ)
        sage: T.point([1,2,3,4,5]);
        (1/3 : 2/3 : 1 , 4/5 : 1)
    """
    def __init__(self, parent, polys, check=True):
        r"""
        The Python constructor.

        INPUT:

        - ``parent`` -- Homset

        - ``polys`` -- anything that defines a point in the class

        - ``check`` -- Boolean. Whether or not to perform input checks
          (Default: ``True``)

        EXAMPLES::

            sage: P1.<x0,x1,x2> = ProjectiveSpace(QQ,2)
            sage: P2 = ProjectiveSpace(QQ,3,'y')
            sage: T = ProductProjectiveSpaces([P1,P2])
            sage: Q1 = P1(1,1,1)
            sage: Q2 = P2(1,2,3,4)
            sage: T([Q1,Q2])
            (1 : 1 : 1 , 1/4 : 1/2 : 3/4 : 1)

        ::

            sage: T = ProductProjectiveSpaces([2,2,2],GF(5),'x')
            sage: T.point([1,2,3,4,5,6,7,8,9])
            (2 : 4 : 1 , 4 : 0 : 1 , 3 : 2 : 1)

        ::

            sage: T.<x,y,z,w> = ProductProjectiveSpaces([1,1],GF(5))
            sage: X = T.subscheme([x-y,z-2*w])
            sage: X([1,1,2,1])
            (1 : 1 , 2 : 1)
        """
        SchemeMorphism.__init__(self, parent)
        if all(isinstance(P, SchemeMorphism_point) for P in polys):
            if check:
                Q = []
                self._points=[]
                for i in range(len(polys)):
                    if polys[i].codomain() != parent.codomain().ambient_space()[i]:
                        raise ValueError("Points must be in correct projective spaces")
                    Q += list(polys[i])
                    self._points.append(polys[i])
                parent.codomain()._check_satisfies_equations(Q)
            self._points = polys
        else:
            N=parent.codomain().ambient_space().dimension_relative_components()
            if check:
                parent.codomain()._check_satisfies_equations(polys)
            splitpolys=self.codomain().ambient_space()._factors(polys)
            self._points = [parent.codomain().ambient_space()[i].point(splitpolys[i], check) for i in range(len(N))]

    def __getitem__(self, i):
        r"""
        Return the `i`-th coordinate point.

        INPUT:

        - `i` - integer

        OUTPUT:

        The projective space point that is the `i`-th coordinate.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([2,2,2],GF(5),'x')
            sage: P = T([1,0,1,1,0,0,0,0,1])
            sage: P[1]
            (1 : 0 : 0)
            sage: P[1].codomain()
            Projective Space of dimension 2 over Finite Field of size 5
            sage: P[1][0]
            1
        """
        return(self._points[i])

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([2,2],ZZ,'x')
            sage: P = T([1,2,3,4,5,6])
            sage: P._repr_()
            '(1 : 2 : 3 , 4 : 5 : 6)'
        """
        return('(%s)'%(" , ".join((" : ".join([repr(f) for f in Q])) for Q in self._points)))

    def __eq__(self, other):
        r"""
        Tests the projective equality of two points.

        INPUT:

        - ``right`` - a point on a product of projective spaces

        OUTPUT:

        - Boolean - ``True`` if ``self`` and ``right`` define the same point.
          ``False`` otherwise.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([2,2],ZZ,'x')
            sage: P = T([1,2,3,4,5,6])
            sage: Q = T([2,4,6,4,5,6])
            sage: P == Q
            True
        """
        for i in range(self.codomain().num_components()):
           if self[i] != other[i]:
               return False
        return True

    def __ne__(self, other):
        r"""
        Tests the projective inequality of two points.

        INPUT:

        - ``right`` -- a point on a product of projective spaces

        OUTPUT:

        - Boolean - ``False`` if ``self`` and ``right`` define the same point.
          ``True`` otherwise.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1,1,1],ZZ,'x')
            sage: P = T([1,2,3,4,5,6])
            sage: Q = T([2,4,6,4,1,0])
            sage: P != Q
            True
        """
        for i in range(self.codomain().num_components()):
           if self[i] != other[i]:
               return True
        return False

    def __cmp__(self, right):
        r"""
        Compare two points in products of projective spaces.

        INPUT:

        - ``other`` -- anything. To compare against the point ``self``.

        OUTPUT:

        ``+1``, ``0``, or ``-1``.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1,1,1],GF(5),'x')
            sage: P = T([3,2,3,4,1,0])
            sage: Q = T([1,2,3,4,3,1])
            sage: P.__cmp__(Q)
            1

        ::

            sage: T = ProductProjectiveSpaces([1,1,1],GF(5),'x')
            sage: P = T([1,2,3,4,1,0])
            sage: Q = T([1,2,3,4,3,0])
            sage: P.__cmp__(Q)
            0

        ::

            sage: T = ProductProjectiveSpaces([1,1,1],GF(5),'x')
            sage: P = T([1,2,3,4,1,0])
            sage: Q = T([1,2,3,4,3,1])
            sage: P.__cmp__(Q)
            -1
        """
        #needed for Digraph
        if not isinstance(right, (ProductProjectiveSpaces_point_ring)):
            return -1
        else:
            return(cmp(self._points, right._points))

    def __copy__(self):
        r"""
        Return a copy of the point ``self``.

        OUTPUT:

        - a point in the same space as ``self``.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1,1],QQ,'x')
            sage: P = T([2,1,0,1])
            sage: Q = P.__copy__()
            sage: P is Q
            False
            sage: P == Q
            True
        """
        P = [copy(self[i]) for i in range(self.codomain().ambient_space().num_components())]
        return(self.codomain().point(P, False))

    def __iter__(self):
        r"""
        Iterate over the coordinates of the point.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1,1],QQ,'x')
            sage: P = T([2,1,0,1])
            sage: iter = P.__iter__()
            sage: next(iter)
            2
            sage: next(iter)
            1
            sage: list(P)
            [2, 1, 0, 1]
        """
        L = []
        for P in self._points:
            L += P._coords
        return iter(L)

    def normalize_coordinates(self):
        r"""
        Removes common factors (componentwise) from the coordinates of ``self`` (including `-1`).

        OUTPUT:

        None.

        EXAMPLES::

            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([2,2],ZZ)
            sage: P = T.point([5,10,15,4,2,6]);
            sage: P.normalize_coordinates()
            sage: P
            (1 : 2 : 3 , 2 : 1 : 3)
        """
        for i in range(self.codomain().ambient_space().num_components()):
            self[i].normalize_coordinates()

    def scale_by(self, t):
        r"""
        Scale the coordinates of the point ``self`` by `t`, done componentwise.

        A ``TypeError`` occurs if the point is not in the base ring of the
        codomain after scaling.

        INPUT:

        - ``t`` -- a ring element

        EXAMPLES::

            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([1,1,1],ZZ)
            sage: P = T.point([5,10,15,4,2,6]);
            sage: P.scale_by([2,1,1])
            sage: P
            (10 : 20 , 15 : 4 , 2 : 6)
        """
        if not isinstance(t, (tuple, list)):
            raise TypeError("%s must be a list or tuple"%t)
        if len(t) != self.codomain().ambient_space().num_components():
            raise TypeError("%s must have same number of components as %r"%(t, self))
        for i in range(self.codomain().ambient_space().num_components()):
            self[i].scale_by(t[i])

    def change_ring(self, R, **kwds):
        r"""
        Returns a new :class:`ProductProjectiveSpaces_point` which is ``self`` coerced to ``R``.

        If the keyword ``check`` is ``True``, then the initialization checks are performed.
        The user may specify the embedding into ``R`` with a keyword.

        INPUT:

        - ``R`` -- ring

        kwds:

        - ``check`` -- Boolean

        - ``embedding`` -- field embedding from the base ring of ``self`` to ``R``.

        OUTPUT:

        :class:`ProductProjectiveSpaces_point`

        EXAMPLES::

            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([1,1,1],ZZ)
            sage: P = T.point([5,3,15,4,2,6]);
            sage: P.change_ring(GF(3))
            (1 : 0 , 0 : 1 , 1 : 0)
        """
        check = kwds.get('check', True)
        S = self.codomain().change_ring(R)
        Q = [P.change_ring(R,**kwds) for P in self._points]
        return(S.point(Q, check))
