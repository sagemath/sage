r"""
Base class for polyhedra over `\ZZ`
"""

# ****************************************************************************
#       Copyright (C) 2011-2013 Volker Braun <vbraun.name@gmail.com>
#                     2015      Nathann Cohen
#                     2015      Vincent Delecroix
#                     2017-2018 Frédéric Chapoton
#                     2019      Sophia Elia
#                     2019-2020 Jonathan Kliem
#                     2023      Luze Xu
#                     2023      Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from .base_QQ import Polyhedron_QQ
from sage.arith.misc import gcd


#########################################################################
class Polyhedron_ZZ(Polyhedron_QQ):
    r"""
    Base class for Polyhedra over `\ZZ`.

    TESTS::

        sage: p = Polyhedron([(0,0)], base_ring=ZZ);  p
        A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex
        sage: TestSuite(p).run()
    """
    _base_ring = ZZ

    def __getattribute__(self, name):
        r"""
        TESTS:

        A lattice polytope does not have a Ehrhart quasipolynomial because it
        is always a polynomial::

            sage: P = polytopes.cube()
            sage: P.__getattribute__(name='ehrhart_quasipolynomial')
            Traceback (most recent call last):
            ...
            AttributeError: ehrhart_quasipolynomial
        """
        if name in ['ehrhart_quasipolynomial']:
            raise AttributeError(name)
        return super().__getattribute__(name)

    def __dir__(self):
        r"""
        TESTS:

        Removes the Ehrhart quasipolynomial from the list of methods for the
        lattice polyhedron::

            sage: P = polytopes.cube()
            sage: 'ehrhart_polynomial' in P.__dir__()
            True
            sage: 'ehrhart_quasipolynomial' in P.__dir__()
            False
        """
        orig_dir = (set(dir(self.__class__)) | set(self.__dict__.keys()))
        return sorted(orig_dir - set(['ehrhart_quasipolynomial']))

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

        TESTS:

        Check :issue:`22622`::

            sage: P1 = Polyhedron(vertices = [[1, 0], [0, 1]], rays = [[1, 1]])
            sage: P1.is_lattice_polytope()
            False
        """
        return self.is_compact()

    def _ehrhart_polynomial_latte(self, verbose=False, dual=None,
            irrational_primal=None, irrational_all_primal=None, maxdet=None,
            no_decomposition=None, compute_vertex_cones=None, smith_form=None,
            dualization=None, triangulation=None, triangulation_max_height=None,
            **kwds):
        r"""
        Return the Ehrhart polynomial of this polyhedron using LattE integrale.

        Let `P` be a lattice polytope in `\RR^d` and define `L(P,t) = \# (tP
        \cap \ZZ^d)`. Then E. Ehrhart proved in 1962 that `L` coincides with a
        rational polynomial of degree `d` for integer `t`. `L` is called the
        *Ehrhart polynomial* of `P`. For more information see the
        :wikipedia:`Ehrhart_polynomial`.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``); if ``True``, print the
          whole output of the LattE command

        The following options are passed to the LattE command, for details you
        should consult `the LattE documentation
        <https://www.math.ucdavis.edu/~latte/software/packages/latte_current/>`__:

        - ``dual`` -- boolean; triangulate and signed-decompose in the dual
          space

        - ``irrational_primal`` -- boolean; triangulate in the dual space,
          signed-decompose in the primal space using irrationalization

        - ``irrational_all_primal`` -- boolean; triangulate and signed-decompose
          in the primal space using irrationalization

        - ``maxdet`` -- integer; decompose down to an index (determinant) of
          ``maxdet`` instead of index 1 (unimodular cones)

        - ``no_decomposition`` -- boolean; do not signed-decompose simplicial cones

        - ``compute_vertex_cones`` -- string; either 'cdd' or 'lrs' or '4ti2'

        - ``smith_form`` -- string; either 'ilio' or 'lidia'

        - ``dualization`` -- string; either 'cdd' or '4ti2'

        - ``triangulation`` -- string; 'cddlib', '4ti2' or 'topcom'

        - ``triangulation_max_height`` -- integer; use a uniform distribution of
          height from 1 to this number

        .. NOTE::

            Any additional argument is forwarded to LattE's executable
            ``count``. All occurrences of '_' will be replaced with a '-'.

        ALGORITHM:

        This method calls the program ``count`` from LattE integrale, a program
        for lattice point enumeration (see
        https://www.math.ucdavis.edu/~latte/).

        .. SEEALSO::

            :mod:`~sage.interfaces.latte` the interface to LattE integrale

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(0,0,0),(3,3,3),(-3,2,1),(1,-1,-2)])
            sage: p = P._ehrhart_polynomial_latte()    # optional - latte_int
            sage: p                                    # optional - latte_int
            7/2*t^3 + 2*t^2 - 1/2*t + 1
            sage: p(1)                                 # optional - latte_int
            6
            sage: len(P.integral_points())
            6
            sage: p(2)                                 # optional - latte_int
            36
            sage: len((2*P).integral_points())
            36

        The unit hypercubes::

            sage: # optional - latte_int
            sage: from itertools import product
            sage: def hypercube(d):
            ....:     return Polyhedron(vertices=list(product([0,1],repeat=d)))
            sage: hypercube(3)._ehrhart_polynomial_latte()
            t^3 + 3*t^2 + 3*t + 1
            sage: hypercube(4)._ehrhart_polynomial_latte()
            t^4 + 4*t^3 + 6*t^2 + 4*t + 1
            sage: hypercube(5)._ehrhart_polynomial_latte()
            t^5 + 5*t^4 + 10*t^3 + 10*t^2 + 5*t + 1
            sage: hypercube(6)._ehrhart_polynomial_latte()
            t^6 + 6*t^5 + 15*t^4 + 20*t^3 + 15*t^2 + 6*t + 1

        TESTS:

        Test options::

            sage: P = Polyhedron(ieqs=[[1,-1,1,0], [-1,2,-1,0], [1,1,-2,0]], eqns=[[-1,2,-1,-3]], base_ring=ZZ)

            sage: p = P._ehrhart_polynomial_latte(maxdet=5, verbose=True)  # optional - latte_int
            This is LattE integrale ...
            ...
            Invocation: ...count --ehrhart-polynomial '--redundancy-check=none' --cdd '--maxdet=5' /dev/stdin
            ...
            sage: p    # optional - latte_int
            1/2*t^2 + 3/2*t + 1

            sage: p = P._ehrhart_polynomial_latte(dual=True, verbose=True)  # optional - latte_int
            This is LattE integrale ...
            ...
            Invocation: ...count --ehrhart-polynomial '--redundancy-check=none' --cdd --dual /dev/stdin
            ...
            sage: p   # optional - latte_int
            1/2*t^2 + 3/2*t + 1

            sage: p = P._ehrhart_polynomial_latte(irrational_primal=True, verbose=True)   # optional - latte_int
            This is LattE integrale ...
            ...
            Invocation: ...count --ehrhart-polynomial '--redundancy-check=none' --cdd --irrational-primal /dev/stdin
            ...
            sage: p   # optional - latte_int
            1/2*t^2 + 3/2*t + 1

            sage: p = P._ehrhart_polynomial_latte(irrational_all_primal=True, verbose=True)  # optional - latte_int
            This is LattE integrale ...
            ...
            Invocation: ...count --ehrhart-polynomial '--redundancy-check=none' --cdd --irrational-all-primal /dev/stdin
            ...
            sage: p   # optional - latte_int
            1/2*t^2 + 3/2*t + 1

        Test bad options::

            sage: P._ehrhart_polynomial_latte(bim_bam_boum=19)   # optional - latte_int
            Traceback (most recent call last):
            ...
            RuntimeError: LattE integrale program failed (exit code 1):
            ...
            Invocation: ...count --ehrhart-polynomial '--redundancy-check=none' --cdd '--bim-bam-boum=19' /dev/stdin
            Unknown command/option --bim-bam-boum=19
        """
        # note: the options below are explicitly written in the function
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

        from sage.interfaces.latte import count
        ine = self.cdd_Hrepresentation()
        return count(ine, cdd=True, ehrhart_polynomial=True, verbose=verbose, **kwds)

    def _ehrhart_polynomial_normaliz(self, variable='t'):
        r"""
        Compute the Ehrhart polynomial of a lattice polytope using Normaliz.

        The backend of ``self`` must be 'normaliz'.

        INPUT:

        - ``variable`` -- string (default: ``'t'``); the variable in which the
          Ehrhart polynomial is expressed

        OUTPUT: a univariate polynomial over a rational field

        EXAMPLES::

            sage: c = Polyhedron(vertices = [[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],backend='normaliz') # optional - pynormaliz
            sage: c._ehrhart_polynomial_normaliz()             # optional - pynormaliz
            t^3 + 3*t^2 + 3*t + 1

        Changing the variable works::

            sage: c._ehrhart_polynomial_normaliz(variable='k') # optional - pynormaliz
            k^3 + 3*k^2 + 3*k + 1

        TESTS:

        Receive a type error if the backend is not normaliz::

            sage: c = Polyhedron(vertices = [[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]])
            sage: c._ehrhart_polynomial_normaliz()
            Traceback (most recent call last):
            ...
            TypeError: The polyhedron's backend should be 'normaliz'
        """
        raise TypeError("The polyhedron's backend should be 'normaliz'")

    @cached_method(do_pickle=True)
    def ehrhart_polynomial(self, engine=None, variable='t', verbose=False,
                           dual=None, irrational_primal=None,
                           irrational_all_primal=None, maxdet=None,
                           no_decomposition=None, compute_vertex_cones=None,
                           smith_form=None, dualization=None,
                           triangulation=None, triangulation_max_height=None,
                           **kwds):
        r"""
        Return the Ehrhart polynomial of this polyhedron.

        Let `P` be a lattice polytope in `\RR^d` and define `L(P,t) = \# (tP
        \cap \ZZ^d)`. Then E. Ehrhart proved in 1962 that `L` coincides with a
        rational polynomial of degree `d` for integer `t`. `L` is called the
        *Ehrhart polynomial* of `P`. For more information see the
        :wikipedia:`Ehrhart_polynomial`.

        The Ehrhart polynomial may be computed using either  LattE Integrale
        or Normaliz by setting ``engine``  to 'latte' or 'normaliz' respectively.

        INPUT:

        - ``engine`` -- string; the backend to use. Allowed values are:

          * ``None`` (default); When no input is given the Ehrhart polynomial
            is computed using LattE Integrale (optional)
          * ``'latte'``; use LattE integrale program (optional)
          * ``'normaliz'``; use Normaliz program (optional). The backend of
            ``self`` must be set to 'normaliz'.

        - ``variable`` -- string (default: ``'t'``); the variable in which the
          Ehrhart polynomial should be expressed

        - When the ``engine`` is 'latte' or None, the additional input values are:

          * ``verbose`` -- boolean (default: ``False``); if ``True``, print the
            whole output of the LattE command.

          The following options are passed to the LattE command, for details
          consult `the LattE documentation
          <https://www.math.ucdavis.edu/~latte/software/packages/latte_current/>`__:

          * ``dual`` -- boolean; triangulate and signed-decompose in the dual
            space
          * ``irrational_primal`` -- boolean; triangulate in the dual space,
            signed-decompose in the primal space using irrationalization.
          * ``irrational_all_primal`` -- boolean; triangulate and signed-decompose
            in the primal space using irrationalization.
          * ``maxdet`` -- integer; decompose down to an index (determinant) of
            ``maxdet`` instead of index 1 (unimodular cones).
          * ``no_decomposition`` -- boolean; do not signed-decompose
            simplicial cones.
          * ``compute_vertex_cones`` -- string; either 'cdd' or 'lrs' or '4ti2'
          * ``smith_form`` -- string; either 'ilio' or 'lidia'
          * ``dualization`` -- string; either 'cdd' or '4ti2'
          * ``triangulation`` -- string; 'cddlib', '4ti2' or 'topcom'
          * ``triangulation_max_height`` -- integer; use a uniform distribution of
            height from 1 to this number

        OUTPUT:

        The Ehrhart polynomial as a univariate polynomial in ``variable``
        over a rational field.

        .. SEEALSO::

            :mod:`~sage.interfaces.latte` the interface to LattE Integrale
            `PyNormaliz <https://pypi.org/project/PyNormaliz>`_

        EXAMPLES:

        To start, we find the Ehrhart polynomial of a three-dimensional
        ``simplex``, first using ``engine='latte'``. Leaving the engine
        unspecified sets the ``engine`` to 'latte' by default::

            sage: simplex = Polyhedron(vertices=[(0,0,0),(3,3,3),(-3,2,1),(1,-1,-2)])
            sage: poly = simplex.ehrhart_polynomial(engine = 'latte')  # optional - latte_int
            sage: poly                                                 # optional - latte_int
            7/2*t^3 + 2*t^2 - 1/2*t + 1
            sage: poly(1)                                              # optional - latte_int
            6
            sage: len(simplex.integral_points())
            6
            sage: poly(2)                                              # optional - latte_int
            36
            sage: len((2*simplex).integral_points())
            36

        Now we find the same Ehrhart polynomial, this time using
        ``engine='normaliz'``. To use the Normaliz engine, the ``simplex`` must
        be defined with ``backend='normaliz'``::

            sage: simplex = Polyhedron(vertices=[(0,0,0),(3,3,3),(-3,2,1),(1,-1,-2)], backend='normaliz') # optional - pynormaliz
            sage: poly = simplex.ehrhart_polynomial(engine='normaliz') # optional - pynormaliz
            sage: poly                                                 # optional - pynormaliz
            7/2*t^3 + 2*t^2 - 1/2*t + 1

        If the ``engine='normaliz'``, the backend should be ``'normaliz'``, otherwise
        it returns an error::

            sage: simplex = Polyhedron(vertices=[(0,0,0),(3,3,3),(-3,2,1),(1,-1,-2)])
            sage: simplex.ehrhart_polynomial(engine='normaliz')
            Traceback (most recent call last):
            ...
            TypeError: The polyhedron's backend should be 'normaliz'

        Now we find the Ehrhart polynomials of the unit hypercubes of
        dimensions three through six. They are computed first with
        ``engine='latte'`` and then with ``engine='normaliz'``.
        The degree of the Ehrhart polynomial matches the dimension of the
        hypercube, and the coefficient of the leading monomial equals the
        volume of the unit hypercube::

            sage: # optional - latte_int
            sage: from itertools import product
            sage: def hypercube(d):
            ....:     return Polyhedron(vertices=list(product([0,1],repeat=d)))
            sage: hypercube(3).ehrhart_polynomial()
            t^3 + 3*t^2 + 3*t + 1
            sage: hypercube(4).ehrhart_polynomial()
            t^4 + 4*t^3 + 6*t^2 + 4*t + 1
            sage: hypercube(5).ehrhart_polynomial()
            t^5 + 5*t^4 + 10*t^3 + 10*t^2 + 5*t + 1
            sage: hypercube(6).ehrhart_polynomial()
            t^6 + 6*t^5 + 15*t^4 + 20*t^3 + 15*t^2 + 6*t + 1

            sage: # optional - pynormaliz
            sage: from itertools import product
            sage: def hypercube(d):
            ....:     return Polyhedron(vertices=list(product([0,1],repeat=d)),backend='normaliz')
            sage: hypercube(3).ehrhart_polynomial(engine='normaliz')
            t^3 + 3*t^2 + 3*t + 1
            sage: hypercube(4).ehrhart_polynomial(engine='normaliz')
            t^4 + 4*t^3 + 6*t^2 + 4*t + 1
            sage: hypercube(5).ehrhart_polynomial(engine='normaliz')
            t^5 + 5*t^4 + 10*t^3 + 10*t^2 + 5*t + 1
            sage: hypercube(6).ehrhart_polynomial(engine='normaliz')
            t^6 + 6*t^5 + 15*t^4 + 20*t^3 + 15*t^2 + 6*t + 1

        An empty polyhedron::

            sage: p = Polyhedron(ambient_dim=3, vertices=[])
            sage: p.ehrhart_polynomial()
            0
            sage: parent(_)
            Univariate Polynomial Ring in t over Rational Field

        The polyhedron should be compact::

            sage: C = Polyhedron(rays=[[1,2],[2,1]])
            sage: C.ehrhart_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: Ehrhart polynomial only defined for compact polyhedra

        TESTS:

        The cache of the Ehrhart polynomial is being pickled::

            sage: P = polytopes.cube()
            sage: P.ehrhart_polynomial()  # optional - latte_int
            8*t^3 + 12*t^2 + 6*t + 1
            sage: Q = loads(dumps(P))
            sage: Q.ehrhart_polynomial.is_in_cache()  # optional - latte_int
            True
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.rational_field import QQ
        R = PolynomialRing(QQ, variable)

        if self.is_empty():
            return R.zero()

        if self.dimension() == 0:
            return R.one()

        if not self.is_compact():
            raise ValueError("Ehrhart polynomial only defined for compact polyhedra")

        if engine is None:
            # setting the default to 'latte'
            engine = 'latte'
        if engine == 'latte':
            poly = self._ehrhart_polynomial_latte(verbose, dual,
            irrational_primal, irrational_all_primal, maxdet,
            no_decomposition, compute_vertex_cones, smith_form,
            dualization, triangulation, triangulation_max_height,
            **kwds)
            return poly.change_variable_name(variable)
            # TO DO: replace this change of variable by creating the appropriate
            #        polynomial ring in the latte interface.

        elif engine == 'normaliz':
            return self._ehrhart_polynomial_normaliz(variable)
        else:
            raise ValueError("engine must be 'latte' or 'normaliz'")

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
            <class 'sage.geometry.polyhedron.parent.Polyhedra_ZZ_ppl_with_category.element_class'>
            sage: p.polar().base_ring()
            Integer Ring

        TESTS:

        Test that :issue:`28551` is fixed::

            sage: polytopes.cube(backend='normaliz').polar().backend()  # optional - pynormaliz
            'normaliz'
        """
        if not self.has_IP_property():
            raise ValueError('The polytope must have the IP property.')

        vertices = tuple(ieq.A() / ieq.b() for
                         ieq in self.inequality_generator())

        ieqs = ((1,) + tuple(v[:]) for v in self.vertices())

        pref_rep = 'Hrep' if self.n_vertices() <= self.n_inequalities() else 'Vrep'

        if all(v_i in ZZ for v in vertices for v_i in v):
            parent = self.parent()
            vertices = (v.change_ring(ZZ) for v in vertices)
        else:
            parent = self.parent().change_ring(QQ)

        return parent.element_class(parent, [vertices, [], []], [ieqs, []],
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    @cached_method
    def is_reflexive(self):
        r"""
        A lattice polytope is reflexive if it contains the origin in its interior
        and its polar with respect to the origin is a lattice polytope.

        Equivalently, it is reflexive if it is of the form `\{x \in \mathbb{R}^d: Ax \leq 1\}`
        for some integer matrix `A` and `d` the ambient dimension.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(1,0,0),(0,1,0),(0,0,1),(-1,-1,-1)], base_ring=ZZ)
            sage: p.is_reflexive()
            True
            sage: polytopes.hypercube(4).is_reflexive()
            True
            sage: p = Polyhedron(vertices=[(1,0), (0,2), (-1,0), (0,-1)], base_ring=ZZ)
            sage: p.is_reflexive()
            False
            sage: p = Polyhedron(vertices=[(1,0), (0,2), (-1,0)], base_ring=ZZ)
            sage: p.is_reflexive()
            False

        An error is raised, if the polyhedron is not compact::

            sage: p = Polyhedron(rays=[(1,)], base_ring=ZZ)
            sage: p.is_reflexive()
            Traceback (most recent call last):
            ...
            ValueError: the polyhedron is not compact
        """
        if not self.is_compact():
            raise ValueError("the polyhedron is not compact")

        for H in self.Hrepresentation():
            if H.is_equation():
                return False
            b = H.b()
            if b < 1:
                return False
            if not all(v_i/b in ZZ for v_i in H.A()):
                return False

        return True

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

        - [PALP]_
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

        - ``dim`` -- integer; the dimension of the lattice polytope
          fiber

        OUTPUT:

        A generator yielding the distinct lattice polytope fibers of
        given dimension.

        EXAMPLES::

            sage: P = Polyhedron(toric_varieties.P4_11169().fan().rays(), base_ring=ZZ)             # needs palp sage.graphs
            sage: list(P.fibration_generator(2))                                        # needs palp sage.graphs
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
        r"""
        Return the translation vector to ``translated_polyhedron``.

        INPUT:

        - ``translated_polyhedron`` -- a polyhedron

        OUTPUT:

        A `\ZZ`-vector that translates ``self`` to
        ``translated_polyhedron``. A :exc:`ValueError` is raised if
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
        r"""
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
            <generator object ..._subpoly_parallel_facets at 0x...>
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
        if self.dim() > 2 or not self.is_compact():
            raise NotImplementedError('only implemented for bounded polygons')
        from sage.geometry.polyhedron.plot import cyclic_sort_vertices_2d
        vertices = cyclic_sort_vertices_2d(self.vertices())
        n = len(vertices)
        if n == 1:  # single point
            yield self
            return
        edge_vectors = []
        for i in range(n):
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
            if point != origin:   # does not close up, not a subpolygon
                continue
            yield parent([v, [], []], None)

    @cached_method
    def minkowski_decompositions(self):
        r"""
        Return all Minkowski sums that add up to the polyhedron.

        OUTPUT:

        A tuple consisting of pairs `(X,Y)` of `\ZZ`-polyhedra that
        add up to ``self``. All pairs up to exchange of the summands
        are returned, that is, `(Y,X)` is not included if `(X,Y)`
        already is.

        EXAMPLES::

            sage: square = Polyhedron(vertices=[(0,0),(1,0),(0,1),(1,1)])
            sage: square.minkowski_decompositions()
            ((A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices),
             (A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices,
              A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices))

        Example from http://cgi.di.uoa.gr/~amantzaf/geo/ ::

            sage: Q = Polyhedron(vertices=[(4,0), (6,0), (0,3), (4,3)])
            sage: R = Polyhedron(vertices=[(0,0), (5,0), (8,4), (3,2)])
            sage: (Q+R).minkowski_decompositions()
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

           sage: [ len(square.dilation(i).minkowski_decompositions())
           ....:   for i in range(6) ]
           [1, 2, 5, 8, 13, 18]
           sage: [ integer_ceil((i^2 + 2*i - 1) / 2) + 1 for i in range(10) ]
           [1, 2, 5, 8, 13, 18, 25, 32, 41, 50]
        """
        if self.dim() > 2 or not self.is_compact():
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
            Y = Y.change_ring(ZZ)  # Minkowski difference returns QQ-polyhedron
            if X + Y != self:
                continue
            decompositions.append((X, Y))
            summands += [X, Y]
        return tuple(decompositions)

    def normal_form(self, algorithm='palp_native', permutation=False):
        r"""
        Return the normal form of vertices of the lattice polytope ``self``.

        INPUT:

        - ``algorithm`` -- must be ``'palp_native'``, the default

        - ``permutation`` -- boolean (default: ``False``); if ``True``, the permutation
          applied to vertices to obtain the normal form is returned as well

        For more more detail,
        see :meth:`~sage.geometry.lattice_polytope.LatticePolytopeClass.normal_form`.

        EXAMPLES:

        We compute the normal form of the "diamond"::

            sage: d = Polyhedron([(1,0), (0,1), (-1,0), (0,-1)])
            sage: d.vertices()
            (A vertex at (-1, 0),
             A vertex at (0, -1),
             A vertex at (0, 1),
             A vertex at (1, 0))
            sage: d.normal_form()                                                       # needs sage.groups
            [(1, 0), (0, 1), (0, -1), (-1, 0)]
            sage: d.lattice_polytope().normal_form("palp_native")                       # needs sage.groups
            M( 1,  0),
            M( 0,  1),
            M( 0, -1),
            M(-1,  0)
            in 2-d lattice M

        Using ``permutation=True``::

            sage: d.normal_form(permutation=True)                                       # needs sage.groups
            ([(1, 0), (0, 1), (0, -1), (-1, 0)], ())

        It is not possible to compute normal forms for polytopes which do not
        span the space::

            sage: p = Polyhedron([(1,0,0), (0,1,0), (-1,0,0), (0,-1,0)])
            sage: p.normal_form()
            Traceback (most recent call last):
            ...
            ValueError: normal form is not defined for lower-dimensional polyhedra, got
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices

        The normal form is also not defined for unbounded polyhedra::

            sage: p = Polyhedron(vertices=[[1, 1]], rays=[[1, 0], [0, 1]], base_ring=ZZ)
            sage: p.normal_form()
            Traceback (most recent call last):
            ...
            ValueError: normal form is not defined for unbounded polyhedra, got
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 2 rays

        See :issue:`15280` for proposed extensions to these cases.

        TESTS::

            sage: d.normal_form(algorithm='palp_fiction')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be 'palp_native'
        """
        from sage.geometry.palp_normal_form import _palp_PM_max, _palp_canonical_order

        if algorithm != "palp_native":
            raise ValueError("algorithm must be 'palp_native'")

        if self.dim() < self.ambient_dim():
            raise ValueError("normal form is not defined for lower-dimensional polyhedra, got %s" % self)

        if not self.is_compact():
            raise ValueError("normal form is not defined for unbounded polyhedra, got %s" % self)

        PM = self.slack_matrix().transpose()
        PM_max, permutations = _palp_PM_max(PM, check=True)
        out = _palp_canonical_order(self.vertices(), PM_max, permutations)

        if permutation:
            return out
        else:
            return out[0]
