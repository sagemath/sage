r"""
Points for products of projective spaces

This class builds on the projective space class and its point and morphism classes.

EXAMPLES:

We construct products projective spaces of various dimensions over the same ring.::

    sage: P1xP1.<x,y, u,v> = ProductProjectiveSpaces(QQ, [1, 1])
    sage: P1xP1([2, 1, 3, 1])
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
from sage.categories.integral_domains import IntegralDomains
from sage.rings.fraction_field import FractionField
from sage.rings.integer_ring import ZZ
from sage.schemes.generic.morphism import SchemeMorphism
from sage.schemes.generic.morphism import SchemeMorphism_point
from sage.structure.sequence import Sequence


class ProductProjectiveSpaces_point_ring(SchemeMorphism_point):
    r"""
    The class of points on products of projective spaces.

    The components are projective space points.

    EXAMPLES::

        sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
        sage: T.point([1, 2, 3, 4, 5]);
        (1/3 : 2/3 : 1 , 4/5 : 1)
    """
    def __init__(self, parent, polys, check=True):
        r"""
        The Python constructor.

        INPUT:

        - ``parent`` -- Hom-set.

        - ``polys`` -- anything that defines a point in the class.

        - ``check`` -- Boolean. Whether or not to perform input checks.
          (Default: ``True``)

        EXAMPLES::

            sage: P1.<x0,x1,x2> = ProjectiveSpace(QQ, 2)
            sage: P2 = ProjectiveSpace(QQ, 3, 'y')
            sage: T = ProductProjectiveSpaces([P1, P2])
            sage: Q1 = P1(1, 1, 1)
            sage: Q2 = P2(1, 2, 3, 4)
            sage: T([Q1, Q2])
            (1 : 1 : 1 , 1/4 : 1/2 : 3/4 : 1)

        ::

            sage: T = ProductProjectiveSpaces([2, 2, 2], GF(5), 'x')
            sage: T.point([1, 2, 3, 4, 5, 6, 7, 8, 9])
            (2 : 4 : 1 , 4 : 0 : 1 , 3 : 2 : 1)

        ::

            sage: T.<x,y,z,w> = ProductProjectiveSpaces([1, 1], GF(5))
            sage: X = T.subscheme([x-y, z-2*w])
            sage: X([1, 1, 2, 1])
            (1 : 1 , 2 : 1)
        """
        polys = copy(polys)
        SchemeMorphism.__init__(self, parent)
        if all(isinstance(P, SchemeMorphism_point) for P in polys):
            if check:
                Q = []
                self._points = []
                for i in range(len(polys)):
                    if polys[i].codomain() != parent.codomain().ambient_space()[i]:
                        raise ValueError("points must be in correct projective spaces")
                    Q += list(polys[i])
                    self._points.append(polys[i])
                parent.codomain()._check_satisfies_equations(Q)
            self._points = polys
        else:
            R = parent.codomain().ambient_space().base_ring()
            polys = Sequence(polys, R)
            N = parent.codomain().ambient_space().dimension_relative_components()
            if check:
                parent.codomain()._check_satisfies_equations(polys)
            splitpolys=self.codomain().ambient_space()._factors(polys)
            self._points = [parent.codomain().ambient_space()[i].point(splitpolys[i], check) for i in range(len(N))]

    def __getitem__(self, i):
        r"""
        Return the ``i``-th coordinate point.

        INPUT:

        - ``i`` - integer.

        OUTPUT:

        The projective space point that is the ``i``-th coordinate.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([2, 2, 2], GF(5), 'x')
            sage: P = T([1, 0, 1, 1, 0, 0, 0, 0, 1])
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
        Return a string representation of this point.

        OUTPUT: String.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([2, 2], ZZ, 'x')
            sage: P = T([1, 2, 3, 4, 5, 6])
            sage: P._repr_()
            '(1 : 2 : 3 , 4 : 5 : 6)'
        """
        return('(%s)'%(" , ".join((" : ".join([repr(f) for f in Q])) for Q in self._points)))

    def __eq__(self, other):
        r"""
        Tests the projective equality of two points.

        INPUT:

        - ``right`` - a point on a product of projective spaces.

        OUTPUT:

        - Boolean - ``True`` if the two points define the same point.
          ``False`` otherwise.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([2, 2], ZZ, 'x')
            sage: P = T([1, 2, 3, 4, 5, 6])
            sage: Q = T([2, 4, 6, 4, 5, 6])
            sage: P == Q
            True
        """
        for i in range(self.codomain().ambient_space().num_components()):
           if self[i] != other[i]:
               return False
        return True

    def __ne__(self, other):
        r"""
        Tests the projective inequality of two points.

        INPUT:

        - ``right`` -- a point on a product of projective spaces.

        OUTPUT:

        - Boolean - ``False`` if the two points define the same point.
          ``True`` otherwise.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 1, 1], ZZ, 'x')
            sage: P = T([1, 2, 3, 4, 5, 6])
            sage: Q = T([2, 4, 6, 4, 1, 0])
            sage: P != Q
            True
        """
        for i in range(self.codomain().ambient_space().num_components()):
           if self[i] != other[i]:
               return True
        return False

    def __cmp__(self, right):
        r"""
        Compare two points in products of projective spaces.

        INPUT:

        - ``other`` -- anything. To compare against this point.

        OUTPUT:

        ``+1``, ``0``, or ``-1``.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 1, 1], GF(5), 'x')
            sage: P = T([3, 2, 3, 4, 1, 0])
            sage: Q = T([1, 2, 3, 4, 3, 1])
            sage: P.__cmp__(Q)
            1

        ::

            sage: T = ProductProjectiveSpaces([1, 1, 1], GF(5), 'x')
            sage: P = T([1, 2, 3, 4, 1, 0])
            sage: Q = T([1, 2, 3, 4, 3, 0])
            sage: P.__cmp__(Q)
            0

        ::

            sage: T = ProductProjectiveSpaces([1, 1, 1], GF(5), 'x')
            sage: P = T([1, 2, 3, 4, 1, 0])
            sage: Q = T([1, 2, 3, 4, 3, 1])
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
        Return a copy of this point.

        OUTPUT:

        - a point in the same space as third point.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 1], QQ, 'x')
            sage: P = T([2, 1, 0, 1])
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

        OUTPUT: An iterator.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 1], QQ, 'x')
            sage: P = T([2, 1, 0, 1])
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

    def __hash__(self):
        """
        Computes the hash value of this point.

        OUTPUT: Integer.

        EXAMPLES::

            sage: PP = ProductProjectiveSpaces(Zmod(6), [1, 1])
            sage: hash(PP([5, 1, 2, 4]))
            1266382469                            # 32-bit
            -855399699883264379                   # 64-bit

        ::

            sage: PP = ProductProjectiveSpaces(ZZ, [1, 2])
            sage: hash(PP([1, 1, 2, 2, 2]))
            805439612                            # 32-bit
            7267864846446758012                  # 64-bit
            sage: hash(PP([1, 1, 1, 1, 1]))
            805439612                            # 32-bit
            7267864846446758012                  # 64-bit

        ::

            sage: PP = ProductProjectiveSpaces(QQ, [1, 1])
            sage: hash(PP([1/7, 1, 2, 1]))
            1139616004                          # 32-bit
            -7585172175017137916                # 64-bit

        ::

            sage: PP = ProductProjectiveSpaces(GF(7), [1, 1, 1])
            sage: hash(PP([4, 1, 5, 4, 6, 1]))
            1796924635                          # 32-bit
            -4539377540667874085                # 64-bit
        """
        R = self.codomain().base_ring()
        # if there is a fraction field normalize the point so that
        # equal points have equal hash values
        if R in IntegralDomains():
            P = self.change_ring(FractionField(R))
            P.normalize_coordinates()
            return hash(tuple(P))
        # if there is no good way to normalize return
        # a constant value
        return hash(self.codomain())

    def normalize_coordinates(self):
        r"""
        Removes common factors (componentwise) from the coordinates of this point (including `-1`).

        OUTPUT: None.

        EXAMPLES::

            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([2, 2], ZZ)
            sage: P = T.point([5, 10, 15, 4, 2, 6]);
            sage: P.normalize_coordinates()
            sage: P
            (1 : 2 : 3 , 2 : 1 : 3)
        """
        for i in range(self.codomain().ambient_space().num_components()):
            self[i].normalize_coordinates()

    def scale_by(self, t):
        r"""
        Scale the coordinates of the point by ``t``, done componentwise.

        A ``TypeError`` occurs if the point is not in the base ring of the
        codomain after scaling.

        INPUT:

        - ``t`` -- a ring element

        EXAMPLES::

            sage: T.<x, y, z, u, v, w> = ProductProjectiveSpaces([1, 1, 1], ZZ)
            sage: P = T.point([5, 10, 15, 4, 2, 6]);
            sage: P.scale_by([2, 1, 1])
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
        Returns a new :class:`ProductProjectiveSpaces_point` which is this point coerced to ``R``.

        If the keyword ``check`` is ``True``, then the initialization checks are performed.
        The user may specify the embedding into ``R`` with a keyword.

        INPUT:

        - ``R`` -- ring.

        kwds:

        - ``check`` -- Boolean.

        - ``embedding`` -- field embedding from the base ring of this point to ``R``.

        OUTPUT:

        :class:`ProductProjectiveSpaces_point`.

        EXAMPLES::

            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([1, 1, 1], ZZ)
            sage: P = T.point([5, 3, 15, 4, 2, 6]);
            sage: P.change_ring(GF(3))
            (1 : 0 , 0 : 1 , 1 : 0)
        """
        check = kwds.get('check', True)
        S = self.codomain().change_ring(R)
        Q = [P.change_ring(R,**kwds) for P in self._points]
        return(S.point(Q, check))

    def nth_iterate(self, f, n, normalize=False):
        r"""
        For a map of this point and a point `P` in ``self.domain()``
        this function returns the nth iterate of `P` by  this point.

        If ``normalize == True``,
        then the coordinates are automatically normalized.

        INPUT:

        - ``f`` -- a ProductProjectiveSpaces_morphism_ring with ``self`` in ``f.domain()``.

        - ``n`` -- a positive integer.

        - ``normalize`` -- Boolean (optional Default: ``False``).

        OUTPUT:

        - A point in ``self.codomain()``

        EXAMPLES::

            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1, 2], QQ)
            sage: H = End(Z)
            sage: f = H([a^3, b^3 + a*b^2, x^2, y^2 - z^2, z*y])
            sage: P = Z([1, 1, 1, 1, 1])
            sage: P.nth_iterate(f, 3)
            (1/1872 : 1 , 1 : 1 : 0)

        ::

            sage: Z.<a,b,x,y> = ProductProjectiveSpaces([1, 1], ZZ)
            sage: H = End(Z)
            sage: f = H([a*b, b^2, x^3 - y^3, y^2*x])
            sage: P = Z([2, 6, 2, 4])
            sage: P.nth_iterate(f, 2, normalize = True)
            (1 : 3 , 407 : 112)

        .. TODO:: Is there a more efficient way to do this?
        """
        if self.codomain() != f.domain():
            raise TypeError("point is not defined over domain of function")
        if not f.is_endomorphism():
            raise TypeError("domain and codomain of function not equal")
        try:
            n = ZZ(n)
        except TypeError:
            raise TypeError("iterate number must be an integer")
        if n < 0:
            raise TypeError("must be a forward orbit")
        if n == 0:
            return(self)
        else:
            Q = f(self)
            if normalize == True:
                Q.normalize_coordinates()
            for i in range(2,n+1):
                Q = f(Q)
                if normalize == True:
                    Q.normalize_coordinates()
            return(Q)

    def orbit(self, f, N, **kwds):
        r"""
        Returns the orbit this point by ``f``.

        If ``N`` is an integer it returns `[P, self(P), \ldots,self^N(P)]`.

        If ``N`` is a list or tuple `N=[m,k]` it returns `[self^m(P),\ldots,self^k(P)`].
        Automatically normalize the points if ``normalize == True``.
        Perform the checks on point initialization if
        ``check==True``

        INPUT:

        - ``f`` -- a :class:`ProductProjectiveSpaces_morphism_ring` with the orbit of `P` in ``f.domain()``.

        - ``N`` -- a non-negative integer or list or tuple of two non-negative integers.

        kwds:

        - ``check`` -- Boolean (optional - default: ``True``).

        - ``normalize`` -- Boolean (optional - default: ``False``).


        OUTPUT:

        - a list of points in ``self.codomain()``.

        EXAMPLES::

            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1, 2], QQ)
            sage: H = End(Z)
            sage: f = H([a^3, b^3 + a*b^2, x^2, y^2 - z^2, z*y])
            sage: P = Z([1, 1, 1, 1, 1])
            sage: P.orbit(f, 3)
            [(1 : 1 , 1 : 1 : 1), (1/2 : 1 , 1 : 0 : 1), (1/12 : 1 , -1 : 1 : 0), (1/1872 : 1 , 1 : 1 : 0)]

        ::

            sage: Z.<a,b,x,y> = ProductProjectiveSpaces([1, 1], ZZ)
            sage: H = End(Z)
            sage: f = H([a*b, b^2, x^3 - y^3, y^2*x])
            sage: P = Z([2, 6, 2, 4])
            sage: P.orbit(f, 3, normalize = True)
            [(1 : 3 , 1 : 2), (1 : 3 , -7 : 4), (1 : 3 , 407 : 112), (1 : 3 , 66014215 : 5105408)]
        """
        if self.codomain() != f.domain():
            raise TypeError("point is not defined over domain of function")
        if not f.is_endomorphism():
            raise TypeError("domain and codomain of function not equal")
        if (isinstance(N, (list,tuple)) == False):
            N = [0, N]
        try:
            N[0] = ZZ(N[0])
            N[1] = ZZ(N[1])
        except TypeError:
            raise TypeError("orbit bounds must be integers")
        if N[0] < 0 or N[1] < 0:
            raise TypeError("orbit bounds must be non-negative")
        if N[0] > N[1]:
            return([])

        Q=copy(self)
        check = kwds.pop("check", True)
        normalize = kwds.pop("normalize", False)

        if normalize == True:
            Q.normalize_coordinates()
        for i in range(1, N[0]+1):
            Q = f(Q, check)
            if normalize == True:
                Q.normalize_coordinates()
        Orb = [Q]
        for i in range(N[0]+1, N[1]+1):
            Q = f(Q, check)
            if normalize == True:
                Q.normalize_coordinates()
            Orb.append(Q)
        return(Orb)

class ProductProjectiveSpaces_point_field(ProductProjectiveSpaces_point_ring):

    def intersection_multiplicity(self, X):
        r"""
        Return the intersection multiplicity of the codomain of this point and subscheme ``X`` at this point.

        This uses the subscheme implementation of intersection_multiplicity. This point must be a point
        on a subscheme of a product of projective spaces.

        INPUT:

        - ``X`` -- a subscheme in the same ambient space as the codomain of this point.

        OUTPUT: An integer.

        EXAMPLES::

            sage: PP.<x,y,z,u,v> = ProductProjectiveSpaces(QQ, [2,1])
            sage: X = PP.subscheme([y^2*z^3*u - x^5*v])
            sage: Y = PP.subscheme([u^3 - v^3, x - y])
            sage: Q = X([0,0,1,1,1])
            sage: Q.intersection_multiplicity(Y)
            2
        """
        from sage.schemes.product_projective.space import is_ProductProjectiveSpaces
        if is_ProductProjectiveSpaces(self.codomain()):
            raise TypeError("this point must be a point on a subscheme of a product of projective spaces")
        return self.codomain().intersection_multiplicity(X, self)

    def multiplicity(self):
        r"""
        Return the multiplicity of this point on its codomain.

        This uses the subscheme implementation of multiplicity. This point must be a point
        on a subscheme of a product of projective spaces.

        OUTPUT: an integer.

        EXAMPLES::

            sage: PP.<x,y,z,w,u,v,t> = ProductProjectiveSpaces(QQ, [3,2])
            sage: X = PP.subscheme([x^8*t - y^8*t + z^5*w^3*v])
            sage: Q1 = X([1,1,0,0,-1,-1,1])
            sage: Q1.multiplicity()
            1
            sage: Q2 = X([0,0,0,1,0,1,1])
            sage: Q2.multiplicity()
            5
            sage: Q3 = X([0,0,0,1,1,0,0])
            sage: Q3.multiplicity()
            6
        """
        from sage.schemes.product_projective.space import is_ProductProjectiveSpaces
        if is_ProductProjectiveSpaces(self.codomain()):
            raise TypeError("this point must be a point on a subscheme of a product of projective spaces")
        return self.codomain().multiplicity(self)

class ProductProjectiveSpaces_point_finite_field(ProductProjectiveSpaces_point_field):
    pass
