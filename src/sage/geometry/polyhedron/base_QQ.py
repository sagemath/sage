"""
Base class for polyhedra over `\QQ`
"""
from __future__ import absolute_import

from sage.rings.all import QQ
from sage.misc.all import prod
from .base import Polyhedron_base


class Polyhedron_QQ(Polyhedron_base):
    """
    Base class for Polyhedra over `\QQ`

    TESTS::

        sage: p = Polyhedron([(0,0)], base_ring=QQ);  p
        A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex
        sage: TestSuite(p).run()
    """
    def _is_zero(self, x):
        """
        Test whether ``x`` is zero.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], base_ring=QQ)
            sage: p._is_zero(0)
            True
            sage: p._is_zero(1/100000)
            False
        """
        return x==0

    def _is_nonneg(self, x):
        """
        Test whether ``x`` is nonnegative.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], base_ring=QQ)
            sage: p._is_nonneg(1)
            True
            sage: p._is_nonneg(-1/100000)
            False
        """
        return x>=0

    def _is_positive(self, x):
        """
        Test whether ``x`` is positive.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], base_ring=QQ)
            sage: p._is_positive(1)
            True
            sage: p._is_positive(0)
            False
        """
        return x>0

    _base_ring = QQ

    def integral_points_count(self, verbose=False, use_Hrepresentation=False,
                              explicit_enumeration_threshold=1000, preprocess=True, **kwds):
        r"""
        Return the number of integral points in the polyhedron.

        This method uses the optional package ``latte_int``
        if an estimate for lattice points based on bounding boxes
        exceeds ``explicit_enumeration_threshold``.

        INPUT:

        - ``verbose`` (boolean; ``False`` by default) -- whether to display
          verbose output.

        - ``use_Hrepresentation`` - (boolean; ``False`` by default) -- whether
          to send the H or V representation to LattE

        - ``preprocess`` - (boolean; ``True`` by default) -- whether, if the integral hull
          is known to lie in a coordinate hyperplane, to tighten bounds to reduce dimension

        .. SEEALSO::

            :mod:`~sage.interfaces.latte` the interface to LattE interfaces

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.integral_points_count()
            27
            sage: P.integral_points_count(explicit_enumeration_threshold=0) # optional - latte_int
            27

        We enlarge the polyhedron to force the use of the generating function methods
        implemented in LattE integrale, rather than explicit enumeration.

            sage: (1000000000*P).integral_points_count(verbose=True) # optional - latte_int
            This is LattE integrale...
            ...
            Total time:...
            8000000012000000006000000001

        We shrink the polyhedron a little bit::

            sage: Q = P*(8/9)
            sage: Q.integral_points_count()
            1
            sage: Q.integral_points_count(explicit_enumeration_threshold=0) # optional - latte_int
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

        "Fibonacci" knapsacks (preprocessing helps a lot)::

            sage: def fibonacci_knapsack(d, b, backend=None):
            ....:     lp = MixedIntegerLinearProgram(base_ring=QQ)
            ....:     x = lp.new_variable(nonnegative=True)
            ....:     lp.add_constraint(lp.sum(fibonacci(i+3)*x[i] for i in range(d)) <= b)
            ....:     return lp.polyhedron(backend=backend)
            sage: fibonacci_knapsack(20, 12).integral_points_count() # does not finish with preprocess=False
            33

        TESTS:

        We check that :trac:`21491` is fixed::

            sage: P = Polyhedron(ieqs=[], eqns=[[-10,0,1],[-10,1,0]])
            sage: P.integral_points_count() # optional - latte_int
            1
            sage: P = Polyhedron(ieqs=[], eqns=[[-11,0,2],[-10,1,0]])
            sage: P.integral_points_count() # optional - latte_int
            0
        """
        if self.is_empty():
            return 0
        if self.dimension() == 0:
            return 1 if self.is_lattice_polytope() else 0
        if not self.is_compact():
            # LattE just prints the warning 'readCddExtFile:: Given polyhedron is unbounded!!!'
            # but then returns 0.
            raise NotImplementedError('Unbounded polyhedra are not supported')

        # for small bounding boxes, it is faster to naively iterate over the points of the box
        box_min, box_max = self.bounding_box(integral_hull=True)
        if box_min is None:
            return 0
        box_points = prod(max_coord-min_coord+1 for min_coord, max_coord in zip(box_min, box_max))

        if explicit_enumeration_threshold is None or box_points <= explicit_enumeration_threshold:
            return len(self.integral_points())

        p = self

        if preprocess:
            # If integral hull is known to lie in a coordinate hyperplane,
            # tighten bounds to reduce dimension.
            rat_box_min, rat_box_max = self.bounding_box(integral=False)
            if any(a == b and (ra < a or b < rb)
                   for ra, a, b, rb in zip(rat_box_min, box_min, box_max, rat_box_max)):
                lp, x = self.to_linear_program(return_variable=True)
                for i, a in enumerate(box_min):
                    lp.set_min(x[i], a)
                for i, b in enumerate(box_max):
                    lp.set_max(x[i], b)
                p = lp.polyhedron() # this recomputes the double description, which is wasteful
                if p.is_empty():
                    return 0
                if p.dimension() == 0:
                    return 1 if p.is_lattice_polytope() else 0

        # LattE does not like V-representation of lower-dimensional polytopes.
        if not p.is_full_dimensional():
            use_Hrepresentation = True

        from sage.interfaces.latte import count
        return count(
                p.cdd_Hrepresentation() if use_Hrepresentation else p.cdd_Vrepresentation(),
                cdd=True,
                verbose=verbose,
                **kwds)
