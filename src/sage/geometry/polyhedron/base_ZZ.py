r"""
Base class for polyhedra over `\ZZ`
"""

#*****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from __future__ import absolute_import

from sage.rings.all import ZZ, QQ
from sage.misc.all import cached_method
from sage.modules.free_module_element import vector
from sage.arith.all import gcd
from .constructor import Polyhedron
from .base import Polyhedron_base


#########################################################################
class Polyhedron_ZZ(Polyhedron_base):
    """
    Base class for Polyhedra over `\ZZ`

    TESTS::

        sage: p = Polyhedron([(0,0)], base_ring=ZZ);  p
        A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex
        sage: TestSuite(p).run(skip='_test_pickling')
    """
    def _is_zero(self, x):
        """
        Test whether ``x`` is zero.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], base_ring=ZZ)
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

            sage: p = Polyhedron([(0,0)], base_ring=ZZ)
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

            sage: p = Polyhedron([(0,0)], base_ring=ZZ)
            sage: p._is_positive(1)
            True
            sage: p._is_positive(0)
            False
        """
        return x>0

    _base_ring = ZZ

    def is_lattice_polytope(self):
        r"""
        Return whether the polyhedron is a lattice polytope.

        OUTPUT:

        ``True`` if the polyhedron is compact and has only integral
        vertices, ``False`` otherwise.

        EXAMPLES::

            sage: polytopes.cross_polytope(3).is_lattice_polytope()
            True
            sage: polytopes.regular_polygon(5).is_lattice_polytope()
            False
        """
        return True

    def ehrhart_polynomial(self, verbose=False, dual=None,
            irrational_primal=None, irrational_all_primal=None, maxdet=None,
            no_decomposition=None, compute_vertex_cones=None, smith_form=None,
            dualization=None, triangulation=None, triangulation_max_height=None,
            **kwds):
        r"""
        Return the Ehrhart polynomial of this polyhedron.

        Let `P` be a lattice polytope in `\RR^d` and define `L(P,t) = \# (tP
        \cap \ZZ^d)`. Then E. Ehrhart proved in 1962 that `L` coincides with a
        rational polynomial of degree `d` for integer `t`. `L` is called the
        *Ehrhart polynomial* of `P`. For more information see the
        :wikipedia:`Ehrhart_polynomial`.

        INPUT:

        - ``verbose`` - (boolean, default to ``False``) if ``True``, print the
          whole output of the LattE command.

        The following options are passed to the LattE command, for details you
        should consult `the LattE documentation
        <https://www.math.ucdavis.edu/~latte/software/packages/latte_current/>`__:

        - ``dual`` - (boolean) triangulate and signed-decompose in the dual
          space

        - ``irrational_primal`` - (boolean) triangulate in the dual space,
          signed-decompose in the primal space using irrationalization.

        - ``irrational_all_primal`` - (boolean) Triangulate and signed-decompose
          in the primal space using irrationalization.

        - ``maxdet`` -- (integer) decompose down to an index (determinant) of
          ``maxdet`` instead of index 1 (unimodular cones).

        - ``no_decomposition`` -- (boolean) do not signed-decompose simplicial cones.

        - ``compute_vertex_cones`` -- (string) either 'cdd' or 'lrs' or '4ti2'

        - ``smith_form`` -- (string) either 'ilio' or 'lidia'

        - ``dualization`` -- (string) either 'cdd' or '4ti2'

        - ``triangulation`` - (string) 'cddlib', '4ti2' or 'topcom'

        - ``triangulation_max_height`` - (integer) use a uniform distribution of
          height from 1 to this number

        .. NOTE::

            Any additional argument is forwarded to LattE's executable
            ``count``. All occurrences of '_' will be replaced with a '-'.

        ALGORITHM:

        This method calls the program ``count`` from LattE integrale, a program
        for lattice point enumeration (see
        https://www.math.ucdavis.edu/~latte/).

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(0,0,0),(3,3,3),(-3,2,1),(1,-1,-2)])
            sage: p = P.ehrhart_polynomial()    # optional - latte_int
            sage: p                             # optional - latte_int
            7/2*t^3 + 2*t^2 - 1/2*t + 1
            sage: p(1)                          # optional - latte_int
            6
            sage: len(P.integral_points())
            6
            sage: p(2)                          # optional - latte_int
            36
            sage: len((2*P).integral_points())
            36

        The unit hypercubes::

            sage: from itertools import product
            sage: def hypercube(d):
            ....:     return Polyhedron(vertices=list(product([0,1],repeat=d)))
            sage: hypercube(3).ehrhart_polynomial()   # optional - latte_int
            t^3 + 3*t^2 + 3*t + 1
            sage: hypercube(4).ehrhart_polynomial()   # optional - latte_int
            t^4 + 4*t^3 + 6*t^2 + 4*t + 1
            sage: hypercube(5).ehrhart_polynomial()   # optional - latte_int
            t^5 + 5*t^4 + 10*t^3 + 10*t^2 + 5*t + 1
            sage: hypercube(6).ehrhart_polynomial()   # optional - latte_int
            t^6 + 6*t^5 + 15*t^4 + 20*t^3 + 15*t^2 + 6*t + 1

        An empty polyhedron::

            sage: P = Polyhedron(ambient_dim=3, vertices=[])
            sage: P.ehrhart_polynomial()    # optional - latte_int
            0
            sage: parent(_)                 # optional - latte_int
            Univariate Polynomial Ring in t over Rational Field

        TESTS:

        Test options::

            sage: P = Polyhedron(ieqs=[[1,-1,1,0], [-1,2,-1,0], [1,1,-2,0]], eqns=[[-1,2,-1,-3]], base_ring=ZZ)

            sage: p = P.ehrhart_polynomial(maxdet=5, verbose=True)  # optional - latte_int
            This is LattE integrale ...
            ...
            Invocation: count --ehrhart-polynomial '--redundancy-check=none' '--maxdet=5' --cdd ...
            ...
            sage: p    # optional - latte_int
            1/2*t^2 + 3/2*t + 1

            sage: p = P.ehrhart_polynomial(dual=True, verbose=True)  # optional - latte_int
            This is LattE integrale ...
            ...
            Invocation: count --ehrhart-polynomial '--redundancy-check=none' --dual --cdd ...
            ...
            sage: p   # optional - latte_int
            1/2*t^2 + 3/2*t + 1

            sage: p = P.ehrhart_polynomial(irrational_primal=True, verbose=True)   # optional - latte_int
            This is LattE integrale ...
            ...
            Invocation: count --ehrhart-polynomial '--redundancy-check=none' --irrational-primal --cdd ...
            ...
            sage: p   # optional - latte_int
            1/2*t^2 + 3/2*t + 1

            sage: p = P.ehrhart_polynomial(irrational_all_primal=True, verbose=True)  # optional - latte_int
            This is LattE integrale ...
            ...
            Invocation: count --ehrhart-polynomial '--redundancy-check=none' --irrational-all-primal --cdd ...
            sage: p   # optional - latte_int
            1/2*t^2 + 3/2*t + 1

        Test bad options::

            sage: P.ehrhart_polynomial(bim_bam_boum=19)   # optional - latte_int
            Traceback (most recent call last):
            ...
            RuntimeError: LattE integrale failed with exit code 1 to execute...
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(QQ, 't')
        if self.is_empty():
            return R.zero()

        from sage.misc.misc import SAGE_TMP
        from subprocess import Popen, PIPE

        ine = self.cdd_Hrepresentation()

        args = ['count', '--ehrhart-polynomial']
        if 'redundancy_check' not in kwds:
            args.append('--redundancy-check=none')

        # note: the options below are explicitely written in the function
        # declaration in order to keep tab completion (see #18211).
        kwds.update({
            'dual'                    : dual,
            'irrational_primal'       : irrational_primal,
            'irrational_all_primal'   : irrational_all_primal,
            'maxdet'                  : maxdet,
            'no_decomposition'        : no_decomposition,
            'compute_vertex_cones'    : compute_vertex_cones,
            'smith_form'              : smith_form,
            'dualization'             : dualization,
            'triangulation'           : triangulation,
            'triangulation_max_height': triangulation_max_height})

        for key,value in kwds.items():
            if value is None or value is False:
                continue

            key = key.replace('_','-')
            if value is True:
                args.append('--{}'.format(key))
            else:
                args.append('--{}={}'.format(key, value))
        args += ['--cdd', '/dev/stdin']

        try:
            # The cwd argument is needed because latte
            # always produces diagnostic output files.
            latte_proc = Popen(args,
                               stdin=PIPE, stdout=PIPE,
                               stderr=(None if verbose else PIPE),
                               cwd=str(SAGE_TMP))
        except OSError:
            from sage.misc.package import PackageNotFoundError
            raise PackageNotFoundError('latte_int')

        ans, err = latte_proc.communicate(ine)
        ret_code = latte_proc.poll()
        if ret_code:
            if err is None:
                err = ", see error message above"
            else:
                err = ":\n" + err
            raise RuntimeError("LattE integrale failed with exit code {} to execute {}".format(ret_code, ' '.join(args)) + err.strip())

        p = ans.splitlines()[-2]

        return R(p)

    @cached_method
    def polar(self):
        """
        Return the polar (dual) polytope.

        The polytope must have the IP-property (see
        :meth:`has_IP_property`), that is, the origin must be an
        interior point. In particular, it must be full-dimensional.

        OUTPUT:

        The polytope whose vertices are the coefficient vectors of the
        inequalities of ``self`` with inhomogeneous term normalized to
        unity.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(1,0,0),(0,1,0),(0,0,1),(-1,-1,-1)], base_ring=ZZ)
            sage: p.polar()
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: type(_)
            <class 'sage.geometry.polyhedron.backend_ppl.Polyhedra_ZZ_ppl_with_category.element_class'>
            sage: p.polar().base_ring()
            Integer Ring
        """
        if not self.has_IP_property():
            raise ValueError('The polytope must have the IP property.')

        vertices = [ ieq.A()/ieq.b() for
                     ieq in self.inequality_generator() ]
        if all( all(v_i in ZZ for v_i in v) for v in vertices):
            return Polyhedron(vertices=vertices, base_ring=ZZ)
        else:
            return Polyhedron(vertices=vertices, base_ring=QQ)

    @cached_method
    def is_reflexive(self):
        """
        EXAMPLES::

            sage: p = Polyhedron(vertices=[(1,0,0),(0,1,0),(0,0,1),(-1,-1,-1)], base_ring=ZZ)
            sage: p.is_reflexive()
            True
        """
        return self.polar().is_lattice_polytope()

    @cached_method
    def has_IP_property(self):
        """
        Test whether the polyhedron has the IP property.

        The IP (interior point) property means that

        * ``self`` is compact (a polytope).

        * ``self`` contains the origin as an interior point.

        This implies that

        * ``self`` is full-dimensional.

        * The dual polyhedron is again a polytope (that is, a compact
          polyhedron), though not necessarily a lattice polytope.

        EXAMPLES::

            sage: Polyhedron([(1,1),(1,0),(0,1)], base_ring=ZZ).has_IP_property()
            False
            sage: Polyhedron([(0,0),(1,0),(0,1)], base_ring=ZZ).has_IP_property()
            False
            sage: Polyhedron([(-1,-1),(1,0),(0,1)], base_ring=ZZ).has_IP_property()
            True

        REFERENCES:

        ..  [PALP]
            Maximilian Kreuzer, Harald Skarke:
            "PALP: A Package for Analyzing Lattice Polytopes
            with Applications to Toric Geometry"
            Comput.Phys.Commun. 157 (2004) 87-106
            :arxiv:`math/0204356`
        """
        return self.is_compact() and self.interior_contains(self.ambient_space().zero())

    def fibration_generator(self, dim):
        """
        Generate the lattice polytope fibrations.

        For the purposes of this function, a lattice polytope fiber is
        a sub-lattice polytope. Projecting the plane spanned by the
        subpolytope to a point yields another lattice polytope, the
        base of the fibration.

        INPUT:

        - ``dim`` -- integer. The dimension of the lattice polytope
          fiber.

        OUTPUT:

        A generator yielding the distinct lattice polytope fibers of
        given dimension.

        EXAMPLES::

            sage: P = Polyhedron(toric_varieties.P4_11169().fan().rays(), base_ring=ZZ)
            sage: list( P.fibration_generator(2) )
            [A 2-dimensional polyhedron in ZZ^4 defined as the convex hull of 3 vertices]
        """
        from sage.combinat.combination import Combinations
        if not self.is_compact():
            raise ValueError('Only polytopes (compact polyhedra) are allowed.')

        nonzero_points = [p for p in self.integral_points() if not p.is_zero()]
        origin = [[0]*self.ambient_dim()]
        fibers = set()
        parent = self.parent()

        for points in Combinations(nonzero_points, dim):
                plane = parent.element_class(parent, [origin,[],points], None)
                if plane.dim() != dim:
                    continue
                fiber = self.intersection(plane)
                if fiber.base_ring() is not ZZ:
                    continue
                fiber_vertices = tuple(sorted(tuple(v) for v in fiber.vertex_generator()))
                if fiber_vertices not in fibers:
                    yield fiber
                    fibers.update([fiber_vertices])
                plane._delete()

    def find_translation(self, translated_polyhedron):
        """
        Return the translation vector to ``translated_polyhedron``.

        INPUT:

        - ``translated_polyhedron`` -- a polyhedron.

        OUTPUT:

        A `\ZZ`-vector that translates ``self`` to
        ``translated_polyhedron``. A ``ValueError`` is raised if
        ``translated_polyhedron`` is not a translation of ``self``,
        this can be used to check that two polyhedra are not
        translates of each other.

        EXAMPLES::

            sage: X = polytopes.cube()
            sage: X.find_translation(X + vector([2,3,5]))
            (2, 3, 5)
            sage: X.find_translation(2*X)
            Traceback (most recent call last):
            ...
            ValueError: polyhedron is not a translation of self
        """
        no_translation_exception = ValueError('polyhedron is not a translation of self')
        if ( set(self.rays()) != set(translated_polyhedron.rays()) or
             set(self.lines()) != set(translated_polyhedron.lines()) or
             self.n_vertices() != translated_polyhedron.n_vertices() ):
            raise no_translation_exception
        sorted_vertices = sorted(map(vector, self.vertices()))
        sorted_translated_vertices = sorted(map(vector, translated_polyhedron.vertices()))
        v = sorted_translated_vertices[0] - sorted_vertices[0]
        if any(vertex+v != translated_vertex
               for vertex, translated_vertex in zip(sorted_vertices, sorted_translated_vertices)):
            raise no_translation_exception
        return v

    def _subpoly_parallel_facets(self):
        """
        Generator for all lattice sub-polyhedra with parallel facets.

        In a sub-polyhedron `Y\subset X` not all edges of `Y` need to
        be parallel to `X`. This method iterates over all
        sub-polyhedra where they are parallel, up to an overall
        translation of the sub-polyhedron. Degenerate sub-polyhedra of
        dimension strictly smaller are included.

        OUTPUT:

        A generator yielding `\ZZ`-polyhedra. By construction, each
        facet of the returned polyhedron is parallel to one of the
        facets of ``self``.

        EXAMPLES::

            sage: X = Polyhedron(vertices=[(0,0), (0,1), (1,0), (1,1)])
            sage: X._subpoly_parallel_facets()
            <generator object _subpoly_parallel_facets at 0x...>
            sage: for p in X._subpoly_parallel_facets():
            ....:     print(p.Vrepresentation())
            (A vertex at (0, 0),)
            (A vertex at (0, -1), A vertex at (0, 0))
            (A vertex at (-1, 0), A vertex at (0, 0))
            (A vertex at (-1, -1), A vertex at (-1, 0), A vertex at (0, -1), A vertex at (0, 0))

        TESTS::

            sage: X = Polyhedron(vertices=[(0,), (3,)])
            sage: [ p.vertices() for p in X._subpoly_parallel_facets() ]
            [(A vertex at (0),),
             (A vertex at (-1), A vertex at (0)),
             (A vertex at (-2), A vertex at (0)),
             (A vertex at (-3), A vertex at (0))]
            sage: list( Polyhedron(vertices=[[0,0]])._subpoly_parallel_facets() )
            [A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex]
            sage: list( Polyhedron()._subpoly_parallel_facets() )
            [The empty polyhedron in ZZ^0]
        """
        if self.dim()>2 or not self.is_compact():
            raise NotImplementedError('only implemented for bounded polygons')
        from sage.geometry.polyhedron.plot import cyclic_sort_vertices_2d
        vertices = cyclic_sort_vertices_2d(self.vertices())
        n = len(vertices)
        if n==1:  # single point
            yield self
            return
        edge_vectors = []
        for i in range(0,n):
            v = vertices[(i+1) % n].vector() - vertices[i].vector()
            d = gcd(list(v))
            v_prim = (v/d).change_ring(ZZ)
            edge_vectors.append([ v_prim*i for i in range(d+1) ])
        origin = self.ambient_space().zero()
        parent = self.parent()
        from itertools import product
        for edges in product(*edge_vectors):
            v = []
            point = origin
            for e in edges:
                point += e
                v.append(point)
            if point!=origin:   # does not close up, not a subpolygon
                continue
            yield parent([v, [], []], None)

    @cached_method
    def Minkowski_decompositions(self):
        """
        Return all Minkowski sums that add up to the polyhedron.

        OUTPUT:

        A tuple consisting of pairs `(X,Y)` of `\ZZ`-polyhedra that
        add up to ``self``. All pairs up to exchange of the summands
        are returned, that is, `(Y,X)` is not included if `(X,Y)`
        already is.

        EXAMPLES::

            sage: square = Polyhedron(vertices=[(0,0),(1,0),(0,1),(1,1)])
            sage: square.Minkowski_decompositions()
            ((A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices),
             (A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices,
              A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices))

        Example from http://cgi.di.uoa.gr/~amantzaf/geo/ ::

            sage: Q = Polyhedron(vertices=[(4,0), (6,0), (0,3), (4,3)])
            sage: R = Polyhedron(vertices=[(0,0), (5,0), (8,4), (3,2)])
            sage: (Q+R).Minkowski_decompositions()
            ((A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 7 vertices),
             (A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices),
             (A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 7 vertices),
             (A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 5 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices),
             (A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 7 vertices),
             (A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 5 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices),
             (A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 7 vertices),
             (A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 6 vertices))

           sage: [ len(square.dilation(i).Minkowski_decompositions())
           ...     for i in range(6) ]
           [1, 2, 5, 8, 13, 18]
           sage: [ ceil((i^2+2*i-1)/2)+1 for i in range(10) ]
           [1, 2, 5, 8, 13, 18, 25, 32, 41, 50]
        """
        if self.dim()>2 or not self.is_compact():
            raise NotImplementedError('only implemented for bounded polygons')
        summands = []
        def is_known_summand(poly):
            for summand in summands:
                try:
                    poly.find_translation(summand)
                    return True
                except ValueError:
                    pass
        decompositions = []
        for X in self._subpoly_parallel_facets():
            if is_known_summand(X):
                continue
            Y = self - X
            if X+Y != self:
                continue
            decompositions.append((X,Y))
            summands += [X, Y]
        return tuple(decompositions)

