r"""
Base class for polyhedra: Methods related to lattice points
"""

# ****************************************************************************
#       Copyright (C) 2008-2012 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011-2015 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2012-2018 Frederic Chapoton
#       Copyright (C) 2013      Andrey Novoseltsev
#       Copyright (C) 2014-2017 Moritz Firsching
#       Copyright (C) 2014-2019 Thierry Monteil
#       Copyright (C) 2015      Nathann Cohen
#       Copyright (C) 2015-2017 Jeroen Demeyer
#       Copyright (C) 2015-2017 Vincent Delecroix
#       Copyright (C) 2015-2018 Dima Pasechnik
#       Copyright (C) 2015-2020 Jean-Philippe Labbe <labbe at math.huji.ac.il>
#       Copyright (C) 2015-2021 Matthias Koeppe
#       Copyright (C) 2016-2019 Daniel Krenn
#       Copyright (C) 2017      Marcelo Forets
#       Copyright (C) 2017-2018 Mark Bell
#       Copyright (C) 2019      Julian Ritter
#       Copyright (C) 2019-2020 Laith Rastanawi
#       Copyright (C) 2019-2020 Sophia Elia
#       Copyright (C) 2019-2021 Jonathan Kliem <jonathan.kliem@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools

from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from .base1 import Polyhedron_base1


class Polyhedron_base2(Polyhedron_base1):
    """
    Methods related to lattice points.

    See :class:`sage.geometry.polyhedron.base.Polyhedron_base`.

    TESTS::

        sage: from sage.geometry.polyhedron.base2 import Polyhedron_base2
        sage: P = polytopes.cube()
        sage: Polyhedron_base2.is_lattice_polytope.f(P)
        True
        sage: Polyhedron_base2.lattice_polytope(P)
        3-d reflexive polytope in 3-d lattice M
        sage: Polyhedron_base2.integral_points(P)
        ((-1, -1, -1),
        (-1, -1, 0),
        (-1, -1, 1),
        (-1, 0, -1),
        (-1, 0, 0),
        (-1, 0, 1),
        (-1, 1, -1),
        (-1, 1, 0),
        (-1, 1, 1),
        (0, -1, -1),
        (0, -1, 0),
        (0, -1, 1),
        (0, 0, -1),
        (0, 0, 0),
        (0, 0, 1),
        (0, 1, -1),
        (0, 1, 0),
        (0, 1, 1),
        (1, -1, -1),
        (1, -1, 0),
        (1, -1, 1),
        (1, 0, -1),
        (1, 0, 0),
        (1, 0, 1),
        (1, 1, -1),
        (1, 1, 0),
        (1, 1, 1))
    """

    @cached_method
    def is_lattice_polytope(self):
        r"""
        Return whether the polyhedron is a lattice polytope.

        OUTPUT:

        ``True`` if the polyhedron is compact and has only integral
        vertices, ``False`` otherwise.

        EXAMPLES::

            sage: polytopes.cross_polytope(3).is_lattice_polytope()
            True
            sage: polytopes.regular_polygon(5).is_lattice_polytope()                    # needs sage.rings.number_field
            False
        """
        if not self.is_compact():
            return False
        if self.base_ring() is ZZ:
            return True
        return all(v.is_integral() for v in self.vertex_generator())

    @cached_method
    def lattice_polytope(self, envelope=False):
        r"""
        Return an encompassing lattice polytope.

        INPUT:

        - ``envelope`` -- boolean (default: ``False``); if the
          polyhedron has non-integral vertices, this option decides
          whether to return a strictly larger lattice polytope or
          raise a :exc:`ValueError`. This option has no effect if the
          polyhedron has already integral vertices.

        OUTPUT:

        A :class:`LatticePolytope
        <sage.geometry.lattice_polytope.LatticePolytopeClass>`. If the
        polyhedron is compact and has integral vertices, the lattice
        polytope equals the polyhedron. If the polyhedron is compact
        but has at least one non-integral vertex, a strictly larger
        lattice polytope is returned.

        If the polyhedron is not compact, a :exc:`NotImplementedError` is
        raised.

        If the polyhedron is not integral and ``envelope=False``, a
        :exc:`ValueError` is raised.

        ALGORITHM:

        For each non-integral vertex, a bounding box of integral
        points is added and the convex hull of these integral points
        is returned.

        EXAMPLES:

        First, a polyhedron with integral vertices::

            sage: P = Polyhedron(vertices=[(1, 0), (0, 1), (-1, 0), (0, -1)])
            sage: lp = P.lattice_polytope(); lp
            2-d reflexive polytope... in 2-d lattice M
            sage: lp                                                            # optional - polytopes_db, needs palp
            2-d reflexive polytope #3 in 2-d lattice M
            sage: lp.vertices()
            M(-1,  0),
            M( 0, -1),
            M( 0,  1),
            M( 1,  0)
            in 2-d lattice M

        Here is a polyhedron with non-integral vertices::

            sage: P = Polyhedron( vertices = [(1/2, 1/2), (0, 1), (-1, 0), (0, -1)])
            sage: lp = P.lattice_polytope()
            Traceback (most recent call last):
            ...
            ValueError: Some vertices are not integral. You probably want
            to add the argument "envelope=True" to compute an enveloping
            lattice polytope.
            sage: lp = P.lattice_polytope(True)
            sage: lp                                                            # optional - polytopes_db, needs palp
            2-d reflexive polytope #5 in 2-d lattice M
            sage: lp.vertices()
            M(-1,  0),
            M( 0, -1),
            M( 1,  1),
            M( 0,  1),
            M( 1,  0)
            in 2-d lattice M
        """
        if not self.is_compact():
            raise NotImplementedError('only compact lattice polytopes are allowed')

        try:
            vertices = self.vertices_matrix(ZZ).columns()
        except TypeError:
            if not envelope:
                raise ValueError('Some vertices are not integral. '
                    'You probably want to add the argument '
                    '"envelope=True" to compute an enveloping lattice polytope.')
            from sage.arith.misc import integer_ceil as ceil
            from sage.arith.misc import integer_floor as floor
            vertices = []
            for v in self.vertex_generator():
                vbox = [ set([floor(x), ceil(x)]) for x in v ]
                vertices.extend( itertools.product(*vbox) )

        # construct the (enveloping) lattice polytope
        from sage.geometry.lattice_polytope import LatticePolytope
        return LatticePolytope(vertices)

    def _integral_points_PALP(self):
        r"""
        Return the integral points in the polyhedron using PALP.

        This method is for testing purposes and will eventually be removed.

        OUTPUT:

        The list of integral points in the polyhedron. If the
        polyhedron is not compact, a :exc:`ValueError` is raised.

        EXAMPLES::

            sage: Polyhedron(vertices=[(-1,-1),(1,0),(1,1),(0,1)])._integral_points_PALP()                              # needs palp
            [M(-1, -1), M(0, 1), M(1, 0), M(1, 1), M(0, 0)]
            sage: Polyhedron(vertices=[(-1/2,-1/2),(1,0),(1,1),(0,1)]).lattice_polytope(True).points()                  # needs palp
            M(-1, -1),
            M(-1,  0),
            M( 0, -1),
            M( 1,  1),
            M( 0,  1),
            M( 1,  0),
            M( 0,  0)
            in 2-d lattice M
            sage: Polyhedron(vertices=[(-1/2,-1/2),(1,0),(1,1),(0,1)])._integral_points_PALP()                          # needs palp
            [M(1, 1), M(0, 1), M(1, 0), M(0, 0)]
        """
        if not self.is_compact():
            raise ValueError('can only enumerate points in a compact polyhedron')
        lp = self.lattice_polytope(True)
        # remove cached values to get accurate timings
        try:
            del lp._points
            del lp._npoints
        except AttributeError:
            pass
        if self.is_lattice_polytope():
            return list(lp.points())
        return [p for p in lp.points() if self.contains(p)]

    @cached_method(do_pickle=True)
    def h_star_vector(self):
        r"""
        Return the `h^*`-vector of the lattice polytope.

        The `h^*`-vector records the coefficients of the polynomial in the
        numerator of the Ehrhart series of a lattice polytope.

        INPUT:

        - ``self`` -- a lattice polytope

        OUTPUT:

        A list whose entries give the `h^*`-vector.

        .. NOTE::

            The backend of ``self`` should be ``'normaliz'``.
            This function depends on Normaliz (i.e. the ``'pynormaliz'`` optional
            package). See the Normaliz documentation for further details.

        EXAMPLES:

        The `h^*`-vector of a unimodular simplex S (a simplex with
        volume = `\frac{1}{dim(S)!}`) is always 1. Here we test this on
        simplices up to dimension 3::

            sage: # optional - pynormaliz
            sage: s1 = polytopes.simplex(1,backend='normaliz')
            sage: s2 = polytopes.simplex(2,backend='normaliz')
            sage: s3 = polytopes.simplex(3,backend='normaliz')
            sage: [s1.h_star_vector(), s2.h_star_vector(), s3.h_star_vector()]
            [[1], [1], [1]]

        For a less trivial example, we compute the `h^*`-vector of the
        `0/1`-cube, which has the Eulerian numbers `(3,i)` for `i \in [0,2]`
        as an `h^*`-vector::

            sage: cube = polytopes.cube(intervals='zero_one', backend='normaliz')   # optional - pynormaliz
            sage: cube.h_star_vector()                                              # optional - pynormaliz
            [1, 4, 1]
            sage: from sage.combinat.combinat import eulerian_number
            sage: [eulerian_number(3,i) for i in range(3)]
            [1, 4, 1]

        TESTS::

            sage: s3 = polytopes.simplex(3)
            sage: s3.h_star_vector()
            Traceback (most recent call last):
            ...
            TypeError: The backend of self must be normaliz

            sage: t = Polyhedron(vertices=[[0],[1/2]])
            sage: t.h_star_vector()
            Traceback (most recent call last):
            ...
            TypeError: The h_star vector is only defined for lattice polytopes

            sage: t2 = Polyhedron(vertices=[[AA(sqrt(2))], [1/2]])                      # needs sage.rings.number_field sage.symbolic
            sage: t2.h_star_vector()                                                    # needs sage.rings.number_field sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: The h_star vector is only defined for lattice polytopes

        Check that the `h^*`-vector is pickled::

            sage: new_cube = loads(dumps(cube))         # optional - pynormaliz
            sage: new_cube.h_star_vector.is_in_cache()  # optional - pynormaliz
            True
        """
        if self.is_empty():
            return 0
        if not self.is_lattice_polytope():
            raise TypeError('The h_star vector is only defined for lattice polytopes')
        if not self.backend() == 'normaliz':
            raise TypeError('The backend of self must be normaliz')
        return self._h_star_vector_normaliz()

    def _h_star_vector_normaliz(self):
        r"""
        Return the `h^*`-vector of a lattice polytope with backend = 'normaliz'.

        INPUT:

        - ``self`` -- a lattice polytope

        OUTPUT:

        The `h^*`-vector as a list.

        .. NOTE::

            The backend of ``self`` should be ``'normaliz'``.

        TESTS::

            sage: s3 = polytopes.simplex(3)
            sage: s3._h_star_vector_normaliz()
            Traceback (most recent call last):
            ...
            TypeError: the backend should be normaliz
        """
        raise TypeError("the backend should be normaliz")

    def integral_points_count(self, **kwds):
        r"""
        Return the number of integral points in the polyhedron.

        This generic version of this method simply calls :meth:`integral_points`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.integral_points_count()
            27

        We shrink the polyhedron a little bit::

            sage: Q = P*(8/9)
            sage: Q.integral_points_count()
            1

        Same for a polyhedron whose coordinates are not rationals.  Note that
        the answer is an integer even though there are no guarantees for
        exactness::

            sage: Q = P*RDF(8/9)
            sage: Q.integral_points_count()
            1

        Unbounded polyhedra (with or without lattice points) are not supported::

            sage: P = Polyhedron(vertices=[[1/2, 1/3]], rays=[[1, 1]])
            sage: P.integral_points_count()
            Traceback (most recent call last):
            ...
            NotImplementedError: ...
            sage: P = Polyhedron(vertices=[[1, 1]], rays=[[1, 1]])
            sage: P.integral_points_count()
            Traceback (most recent call last):
            ...
            NotImplementedError: ...
        """
        return len(self.integral_points())

    def integral_points(self, threshold=100000):
        r"""
        Return the integral points in the polyhedron.

        Uses either the naive algorithm (iterate over a rectangular
        bounding box) or triangulation + Smith form.

        INPUT:

        - ``threshold`` -- integer (default: 100000); use the naive
          algorithm as long as the bounding box is smaller than this

        OUTPUT:

        The list of integral points in the polyhedron. If the
        polyhedron is not compact, a :exc:`ValueError` is raised.

        EXAMPLES::

            sage: Polyhedron(vertices=[(-1,-1),(1,0),(1,1),(0,1)]).integral_points()
            ((-1, -1), (0, 0), (0, 1), (1, 0), (1, 1))

            sage: simplex = Polyhedron([(1,2,3), (2,3,7), (-2,-3,-11)])
            sage: simplex.integral_points()
            ((-2, -3, -11), (0, 0, -2), (1, 2, 3), (2, 3, 7))

        The polyhedron need not be full-dimensional::

            sage: simplex = Polyhedron([(1,2,3,5), (2,3,7,5), (-2,-3,-11,5)])
            sage: simplex.integral_points()
            ((-2, -3, -11, 5), (0, 0, -2, 5), (1, 2, 3, 5), (2, 3, 7, 5))

            sage: point = Polyhedron([(2,3,7)])
            sage: point.integral_points()
            ((2, 3, 7),)

            sage: empty = Polyhedron()
            sage: empty.integral_points()
            ()

        Here is a simplex where the naive algorithm of running over
        all points in a rectangular bounding box no longer works fast
        enough::

            sage: v = [(1,0,7,-1), (-2,-2,4,-3), (-1,-1,-1,4), (2,9,0,-5), (-2,-1,5,1)]
            sage: simplex = Polyhedron(v); simplex
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 5 vertices
            sage: len(simplex.integral_points())
            49

        A case where rounding in the right direction goes a long way::

            sage: P = 1/10*polytopes.hypercube(14, backend='field')
            sage: P.integral_points()
            ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),)

        Finally, the 3-d reflexive polytope number 4078::

            sage: v = [(1,0,0), (0,1,0), (0,0,1), (0,0,-1), (0,-2,1),
            ....:      (-1,2,-1), (-1,2,-2), (-1,1,-2), (-1,-1,2), (-1,-3,2)]
            sage: P = Polyhedron(v)
            sage: pts1 = P.integral_points()                     # Sage's own code
            sage: all(P.contains(p) for p in pts1)
            True
            sage: pts2 = LatticePolytope(v).points()                                    # needs palp
            sage: for p in pts1: p.set_immutable()
            sage: set(pts1) == set(pts2)                                                # needs palp
            True

            sage: timeit('Polyhedron(v).integral_points()')   # not tested - random
            625 loops, best of 3: 1.41 ms per loop
            sage: timeit('LatticePolytope(v).points()')       # not tested - random
            25 loops, best of 3: 17.2 ms per loop

        TESTS:

        Test some trivial cases (see :issue:`17937`)::

            sage: P = Polyhedron(ambient_dim=1)  # empty polyhedron in 1 dimension
            sage: P.integral_points()
            ()
            sage: P = Polyhedron(ambient_dim=0)  # empty polyhedron in 0 dimensions
            sage: P.integral_points()
            ()
            sage: P = Polyhedron([[3]])  # single point in 1 dimension
            sage: P.integral_points()
            ((3),)
            sage: P = Polyhedron([[1/2]])  # single non-integral point in 1 dimension
            sage: P.integral_points()
            ()
            sage: P = Polyhedron([[]])  # single point in 0 dimensions
            sage: P.integral_points()
            ((),)

        Test unbounded polyhedron::

            sage: P = Polyhedron(rays=[[1,0,0]])
            sage: P.integral_points()
            Traceback (most recent call last):
            ...
            ValueError: can only enumerate points in a compact polyhedron
        """
        from sage.misc.misc_c import prod
        if not self.is_compact():
            raise ValueError('can only enumerate points in a compact polyhedron')
        # Trivial cases: polyhedron with 0 or 1 vertices
        if self.n_vertices() == 0:
            return ()
        if self.n_vertices() == 1:
            v = self.vertices_list()[0]
            try:
                return (vector(ZZ, v),)
            except TypeError:  # vertex not integral
                return ()

        # for small bounding boxes, it is faster to naively iterate over the points of the box
        box_min, box_max = self.bounding_box(integral_hull=True)
        if box_min is None:
            return ()
        box_points = prod(max_coord-min_coord+1 for min_coord, max_coord in zip(box_min, box_max))
        if not self.is_lattice_polytope() or \
                (self.is_simplex() and box_points < 1000) or \
                box_points < threshold:
            from sage.geometry.integral_points import rectangular_box_points
            return rectangular_box_points(list(box_min), list(box_max), self)

        # for more complicate polytopes, triangulate & use smith normal form
        from sage.geometry.integral_points import simplex_points
        if self.is_simplex():
            return simplex_points(self.Vrepresentation())
        triangulation = self.triangulate()
        points = set()
        for simplex in triangulation:
            triang_vertices = [self.Vrepresentation(i) for i in simplex]
            new_points = simplex_points(triang_vertices)
            for p in new_points:
                p.set_immutable()
            points.update(new_points)
        # assert all(self.contains(p) for p in points)   # slow
        return tuple(points)

    def get_integral_point(self, index, **kwds):
        r"""
        Return the ``index``-th integral point in this polyhedron.

        This is equivalent to ``sorted(self.integral_points())[index]``.
        However, so long as :meth:`integral_points_count` does not need to
        enumerate all integral points, neither does this method. Hence it can
        be significantly faster. If the polyhedron is not compact, a
        :exc:`ValueError` is raised.

        INPUT:

        - ``index`` -- integer; the index of the integral point to be found. If
          this is not in [0, ``self.integral_point_count()``), an :exc:`IndexError`
          is raised.

        - ``**kwds`` -- optional keyword parameters that are passed to
          :meth:`integral_points_count`

        ALGORITHM:

        The function computes each of the components of the requested point in
        turn. To compute x_i, the `i`-th component, it bisects the upper and lower
        bounds on x_i given by the bounding box. At each bisection, it uses
        :meth:`integral_points_count` to determine on which side of the
        bisecting hyperplane the requested point lies.

        .. SEEALSO::

            :meth:`integral_points_count`.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(-1,-1),(1,0),(1,1),(0,1)])
            sage: P.get_integral_point(1)
            (0, 0)
            sage: P.get_integral_point(4)
            (1, 1)
            sage: sorted(P.integral_points())
            [(-1, -1), (0, 0), (0, 1), (1, 0), (1, 1)]
            sage: P.get_integral_point(5)
            Traceback (most recent call last):
            ...
            IndexError: ...

            sage: Q = Polyhedron([(1,3), (2, 7), (9, 77)])
            sage: [Q.get_integral_point(i) for i in range(Q.integral_points_count())] == sorted(Q.integral_points())
            True
            sage: Q.get_integral_point(0, explicit_enumeration_threshold=0, triangulation='cddlib')  # optional - latte_int
            (1, 3)
            sage: Q.get_integral_point(0, explicit_enumeration_threshold=0, triangulation='cddlib', foo=True)  # optional - latte_int
            Traceback (most recent call last):
            ...
            RuntimeError: ...

            sage: R = Polyhedron(vertices=[[1/2, 1/3]], rays=[[1, 1]])
            sage: R.get_integral_point(0)
            Traceback (most recent call last):
            ...
            ValueError: ...
        """
        from sage.arith.misc import integer_ceil as ceil
        from sage.arith.misc import integer_floor as floor
        if not self.is_compact():
            raise ValueError('can only enumerate points in a compact polyhedron')

        if not 0 <= index < self.integral_points_count(**kwds):
            raise IndexError('polytope index out of range')

        D = self.ambient_dim()
        lower_bounds, upper_bounds = self.bounding_box()
        coordinate = []
        P = self
        S = self.parent()
        for i in range(D):  # Now compute x_i, the ith component of coordinate.
            lower, upper = ceil(lower_bounds[i]), floor(upper_bounds[i]) + 1  # So lower <= x_i < upper.
            while lower < upper-1:
                guess = (lower + upper) // 2  # > lower.
                # Build new polyhedron by intersecting P with the halfspace {x_i < guess}.
                P_lt_guess = P.intersection(S(None, ([[guess-1] + [0] * i + [-1] + [0] * (D - i - 1)], [])))
                # Avoid computing P_geq_guess = P.intersection({x_i >= guess}) right now, it might not be needed.
                P_lt_guess_count = P_lt_guess.integral_points_count(**kwds)
                if P_lt_guess_count > index:  # Move upper down to guess.
                    upper = guess
                    index -= 0
                    P = P_lt_guess
                else:  # P_lt_guess_count <= index:  # Move lower up to guess.
                    lower = guess
                    index -= P_lt_guess_count
                    P_geq_guess = P.intersection(S(None, ([[-guess] + [0] * i + [1] + [0] * (D - i - 1)], [])))
                    P = P_geq_guess
            coordinate.append(lower)  # Record the new component that we have found.
        point = vector(ZZ, coordinate)
        point.set_immutable()
        return point

    def random_integral_point(self, **kwds):
        r"""
        Return an integral point in this polyhedron chosen uniformly at random.

        INPUT:

        - ``**kwds`` -- optional keyword parameters that are passed to
          :meth:`get_integral_point`

        OUTPUT:

        The integral point in the polyhedron chosen uniformly at random. If the
        polyhedron is not compact, a :exc:`ValueError` is raised. If the
        polyhedron does not contain any integral points, an
        :class:`~sage.categories.sets_cat.EmptySetError` is raised.

        .. SEEALSO::

            :meth:`get_integral_point`.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(-1,-1),(1,0),(1,1),(0,1)])
            sage: P.random_integral_point()                                     # random
            (0, 0)
            sage: P.random_integral_point() in P.integral_points()
            True
            sage: P.random_integral_point(explicit_enumeration_threshold=0,     # random, optional - latte_int
            ....:                         triangulation='cddlib')
            (1, 1)
            sage: P.random_integral_point(explicit_enumeration_threshold=0,     # optional - latte_int
            ....:                         triangulation='cddlib', foo=7)
            Traceback (most recent call last):
            ...
            RuntimeError: ...

            sage: Q = Polyhedron(vertices=[(2, 1/3)], rays=[(1, 2)])
            sage: Q.random_integral_point()
            Traceback (most recent call last):
            ...
            ValueError: ...

            sage: R = Polyhedron(vertices=[(1/2, 0), (1, 1/2), (0, 1/2)])
            sage: R.random_integral_point()
            Traceback (most recent call last):
            ...
            EmptySetError: ...
        """
        from sage.misc.randstate import current_randstate

        if not self.is_compact():
            raise ValueError('can only sample integral points in a compact polyhedron')

        count = self.integral_points_count()
        if count == 0:
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError('polyhedron does not contain any integral points')

        return self.get_integral_point(current_randstate().python_random().randint(0, count-1), **kwds)

    def generating_function_of_integral_points(self, **kwds):
        r"""
        Return the multivariate generating function of the
        integral points of this polyhedron.

        To be precise, this returns

        .. MATH::

            \sum_{(r_0,\dots,r_{d-1}) \in \mathit{polyhedron}\cap \ZZ^d}
            y_0^{r_0} \dots y_{d-1}^{r_{d-1}}.

        This calls
        :func:`~sage.geometry.polyhedron.generating_function.generating_function_of_integral_points`,
        so have a look at the documentation and examples there.

        INPUT:

        The following keyword arguments are passed to
        :func:`~sage.geometry.polyhedron.generating_function.generating_function_of_integral_points`:

        - ``split`` -- boolean (default: ``False``) or list

          - ``split=False`` computes the generating function directly,
            without any splitting.

          - When ``split`` is a list of disjoint polyhedra, then
            for each of these polyhedra, this polyhedron is intersected with it,
            its generating function computed and all these generating functions
            are summed up.

          - ``split=True`` splits into `d!` disjoint polyhedra.

        - ``result_as_tuple`` -- (default: ``None``) a boolean or ``None``

          This specifies whether the output is a (partial) factorization
          (``result_as_tuple=False``) or a sum of such (partial)
          factorizations (``result_as_tuple=True``). By default
          (``result_as_tuple=None``), this is automatically determined.
          If the output is a sum, it is represented as a tuple whose
          entries are the summands.

        - ``indices`` -- (default: ``None``) a list or tuple

          If this
          is ``None``, this is automatically determined.

        - ``name`` -- (default: ``'y'``) a string

          The variable names of the Laurent polynomial ring of the output
          are this string followed by an integer.

        - ``names`` -- list or tuple of names (strings), or a comma separated string

          ``name`` is extracted from ``names``, therefore ``names`` has to contain
          exactly one variable name, and ``name`` and``names`` cannot be specified
          both at the same time.

        - ``Factorization_sort`` (default: ``False``) and
          ``Factorization_simplify`` (default: ``True``) -- booleans

          These are passed on to
          :class:`sage.structure.factorization.Factorization` when creating
          the result.

        - ``sort_factors`` -- boolean (default: ``False``)

          If set, then the factors of the output are sorted such that the
          numerator is first and only then all factors of the denominator. It
          is ensured that the sorting is always the same; use this for
          doctesting.

        OUTPUT:

        The generating function as a (partial)
        :class:`~sage.structure.factorization.Factorization`
        of the result whose factors are Laurent polynomials

        The result might be a tuple of such factorizations
        (depending on the parameter ``result_as_tuple``) as well.

        .. NOTE::

            At the moment, only polyhedra with nonnegative coordinates
            (i.e. a polyhedron in the nonnegative orthant) are handled.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: P2 = (Polyhedron(ieqs=[(0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, -1)]),
            ....:       Polyhedron(ieqs=[(0, -1, 0, 1), (0, 1, 0, 0), (0, 0, 1, 0)]))
            sage: P2[0].generating_function_of_integral_points(sort_factors=True)
            1 * (-y0 + 1)^-1 * (-y1 + 1)^-1 * (-y0*y2 + 1)^-1
            sage: P2[1].generating_function_of_integral_points(sort_factors=True)
            1 * (-y1 + 1)^-1 * (-y2 + 1)^-1 * (-y0*y2 + 1)^-1
            sage: (P2[0] & P2[1]).Hrepresentation()
            (An equation (1, 0, -1) x + 0 == 0,
             An inequality (1, 0, 0) x + 0 >= 0,
             An inequality (0, 1, 0) x + 0 >= 0)
            sage: (P2[0] & P2[1]).generating_function_of_integral_points(sort_factors=True)
            1 * (-y1 + 1)^-1 * (-y0*y2 + 1)^-1

        The number of integer partitions
        `1 \leq r_0 \leq r_1 \leq r_2 \leq r_3 \leq r_4`::

            sage: # needs sage.combinat
            sage: P = Polyhedron(ieqs=[(-1, 1, 0, 0, 0, 0), (0, -1, 1, 0, 0, 0),
            ....:                      (0, 0, -1, 1, 0, 0), (0, 0, 0, -1, 1, 0),
            ....:                      (0, 0, 0, 0, -1, 1)])
            sage: f = P.generating_function_of_integral_points(sort_factors=True); f
            y0*y1*y2*y3*y4 * (-y4 + 1)^-1 * (-y3*y4 + 1)^-1 * (-y2*y3*y4 + 1)^-1 *
            (-y1*y2*y3*y4 + 1)^-1 * (-y0*y1*y2*y3*y4 + 1)^-1
            sage: f = f.value()
            sage: P.<z> = PowerSeriesRing(ZZ)
            sage: c = f.subs({y: z for y in f.parent().gens()}); c
            z^5 + z^6 + 2*z^7 + 3*z^8 + 5*z^9 + 7*z^10 + 10*z^11 + 13*z^12 + 18*z^13 +
            23*z^14 + 30*z^15 + 37*z^16 + 47*z^17 + 57*z^18 + 70*z^19 + 84*z^20 +
            101*z^21 + 119*z^22 + 141*z^23 + 164*z^24 + O(z^25)
            sage: ([Partitions(k, length=5).cardinality() for k in range(5,20)] ==
            ....:     c.truncate().coefficients(sparse=False)[5:20])
            True

        .. SEEALSO::

            More examples can be found at
            :func:`~sage.geometry.polyhedron.generating_function.generating_function_of_integral_points`.
        """
        from .generating_function import generating_function_of_integral_points
        return generating_function_of_integral_points(self, **kwds)
