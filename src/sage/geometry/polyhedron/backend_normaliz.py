# sage.doctest: optional - pynormaliz
"""
The Normaliz backend for polyhedral computations

.. NOTE::

    This backend requires `PyNormaliz <https://pypi.org/project/PyNormaliz>`_.
    To install PyNormaliz, type :code:`sage -i pynormaliz` in the terminal.

AUTHORS:

- Matthias Köppe (2016-12): initial version
- Jean-Philippe Labbé (2019-04): Expose normaliz features and added functionalities
"""
# ****************************************************************************
#  Copyright (C) 2016-2022 Matthias Köppe <mkoeppe at math.ucdavis.edu>
#                2016-2018 Travis Scrimshaw
#                2017      Jeroen Demeyer
#                2018-2020 Jean-Philippe Labbé
#                2019      Vincent Delecroix
#                2019-2021 Jonathan Kliem
#                2019-2021 Sophia Elia
#                2020      Frédéric Chapoton
#                2022      Yuan Zhou
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.element import Element
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.misc.lazy_import import lazy_import
import sage.features.normaliz
lazy_import('PyNormaliz', ['NmzResult', 'NmzCompute', 'NmzCone', 'NmzConeCopy'],
            feature=sage.features.normaliz.PyNormaliz())

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.arith.functions import LCM_list
from sage.misc.functional import denominator
from sage.matrix.constructor import vector

from .base_QQ import Polyhedron_QQ
from .base_ZZ import Polyhedron_ZZ
from .base_number_field import Polyhedron_base_number_field


def _format_function_call(fn_name, *v, **k):
    """
    Return a Python function call as a string.

    Keywords are sorted.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.backend_normaliz import _format_function_call
        sage: _format_function_call('foo', 17, hellooooo='goodbyeeee')
        "foo(17, hellooooo='goodbyeeee')"
    """
    args = [repr(a) for a in v] + ["%s=%r" % (arg, val) for arg, val in sorted(k.items())]
    return "{}({})".format(fn_name, ", ".join(args))


#########################################################################
class Polyhedron_normaliz(Polyhedron_base_number_field):
    """
    Polyhedra with normaliz.

    INPUT:

    - ``parent`` -- :class:`~sage.geometry.polyhedron.parent.Polyhedra`
      the parent

    - ``Vrep`` -- list ``[vertices, rays, lines]`` or ``None``; the
      V-representation of the polyhedron; if ``None``, the polyhedron
      is determined by the H-representation

    - ``Hrep`` -- list ``[ieqs, eqns]`` or ``None``; the
      H-representation of the polyhedron; if ``None``, the polyhedron
      is determined by the V-representation

    - ``normaliz_cone`` -- a PyNormaliz wrapper of a normaliz cone

    Only one of ``Vrep``, ``Hrep``, or ``normaliz_cone`` can be different
    from ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0), (1,0), (0,1)],
        ....:                rays=[(1,1)], lines=[],
        ....:                backend='normaliz')
        sage: TestSuite(p).run()

    Two ways to get the full space::

        sage: Polyhedron(eqns=[[0, 0, 0]], backend='normaliz')
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines
        sage: Polyhedron(ieqs=[[0, 0, 0]], backend='normaliz')
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines

    A lower-dimensional affine cone; we test that there are no mysterious
    inequalities coming in from the homogenization::

        sage: P = Polyhedron(vertices=[(1, 1)], rays=[(0, 1)],
        ....:                backend='normaliz')
        sage: P.n_inequalities()
        1
        sage: P.equations()
        (An equation (1, 0) x - 1 == 0,)

    The empty polyhedron::

        sage: P = Polyhedron(ieqs=[[-2, 1, 1], [-3, -1, -1], [-4, 1, -2]],
        ....:                backend='normaliz')
        sage: P
        The empty polyhedron in QQ^2
        sage: P.Vrepresentation()
        ()
        sage: P.Hrepresentation()
        (An equation -1 == 0,)

    TESTS:

    Tests copied from various methods in :mod:`sage.geometry.polyhedron.base`::

        sage: p = Polyhedron(vertices=[[1,0,0], [0,1,0], [0,0,1]],
        ....:                backend='normaliz')
        sage: p.n_equations()
        1
        sage: p.n_inequalities()
        3

        sage: p = Polyhedron(vertices=[[t,t^2,t^3] for t in range(6)],
        ....:                backend='normaliz')
        sage: p.n_facets()
        8

        sage: p = Polyhedron(vertices=[[1,0], [0,1], [1,1]], rays=[[1,1]],
        ....:                backend='normaliz')
        sage: p.n_vertices()
        2

        sage: p = Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,1]],
        ....:                backend='normaliz')
        sage: p.n_rays()
        1

        sage: p = Polyhedron(vertices=[[0,0]], rays=[[0,1], [0,-1]],
        ....:                backend='normaliz')
        sage: p.n_lines()
        1

    Algebraic polyhedra::

        sage: P = Polyhedron(vertices=[[1], [sqrt(2)]],                                 # needs sage.rings.number_field sage.symbolic
        ....:                backend='normaliz', verbose=True)
        # ----8<---- Equivalent Normaliz input file ----8<----
        amb_space 1
        number_field min_poly (a^2 - 2) embedding [1.414213562373095 +/- 2.99e-16]
        cone 0
        subspace 0
        vertices 2
         1 1
         (a) 1
        # ----8<-------------------8<-------------------8<----
        # Calling PyNormaliz.NmzCone(cone=[], number_field=['a^2 - 2', 'a', '[1.414213562373095 +/- 2.99e-16]'], subspace=[], vertices=[[1, 1], [[[0, 1], [1, 1]], 1]])
        sage: P                                                                         # needs sage.rings.number_field sage.symbolic
        A 1-dimensional polyhedron in (Symbolic Ring)^1 defined as
         the convex hull of 2 vertices
        sage: P.vertices()                                                              # needs sage.rings.number_field sage.symbolic
        (A vertex at (1), A vertex at (sqrt(2)))

        sage: P = polytopes.icosahedron(exact=True,                                     # needs sage.rings.number_field
        ....:                           backend='normaliz'); P
        A 3-dimensional polyhedron in
         (Number Field in sqrt5 with defining polynomial x^2 - 5
          with sqrt5 = 2.236067977499790?)^3
         defined as the convex hull of 12 vertices

        sage: x = polygen(ZZ)
        sage: P = Polyhedron(vertices=[[sqrt(2)],                                       # needs sage.rings.number_field sage.symbolic
        ....:                          [AA.polynomial_root(x^3 - 2, RIF(0,3))]],
        ....:                backend='normaliz', verbose=True)
        # ----8<---- Equivalent Normaliz input file ----8<----
        amb_space 1
        number_field min_poly (a^6 - 2) embedding [1.122462048309373 +/- 5.38e-16]
        cone 0
        subspace 0
        vertices 2
         (a^3) 1
         (a^2) 1
        # ----8<-------------------8<-------------------8<----
        # Calling PyNormaliz.NmzCone(cone=[], number_field=['a^6 - 2', 'a', '[1.122462048309373 +/- 5.38e-16]'], subspace=[], vertices=[[[[0, 1], [0, 1], [0, 1], [1, 1], [0, 1], [0, 1]], 1], [[[0, 1], [0, 1], [1, 1], [0, 1], [0, 1], [0, 1]], 1]])
        sage: P                                                                         # needs sage.rings.number_field sage.symbolic
        A 1-dimensional polyhedron in (Symbolic Ring)^1 defined as
         the convex hull of 2 vertices
        sage: P.vertices()                                                              # needs sage.rings.number_field sage.symbolic
        (A vertex at (2^(1/3)), A vertex at (sqrt(2)))
    """
    def __init__(self, parent, Vrep, Hrep, normaliz_cone=None, normaliz_data=None, internal_base_ring=None, **kwds):
        """
        Initialize the polyhedron.

        See :class:`Polyhedron_normaliz` for a description of the input
        data.

        TESTS::

            sage: p = Polyhedron(backend='normaliz')
            sage: TestSuite(p).run()
            sage: p = Polyhedron(vertices=[(1, 1)], rays=[(0, 1)],
            ....:                backend='normaliz')
            sage: TestSuite(p).run()
            sage: p = Polyhedron(vertices=[(-1,-1), (1,0), (1,1), (0,1)],
            ....:                backend='normaliz')
            sage: TestSuite(p).run()
        """
        if normaliz_cone:
            if Hrep is not None or Vrep is not None or normaliz_data is not None:
                raise ValueError("only one of Vrep, Hrep, normaliz_cone, or normaliz_data can be different from None")
            Element.__init__(self, parent=parent)
            self._init_from_normaliz_cone(normaliz_cone, internal_base_ring)
        elif normaliz_data:
            if Hrep is not None or Vrep is not None:
                raise ValueError("only one of Vrep, Hrep, normaliz_cone, or normaliz_data can be different from None")
            Element.__init__(self, parent=parent)
            self._init_from_normaliz_data(normaliz_data, internal_base_ring)
        else:
            if internal_base_ring:
                raise ValueError("if Vrep or Hrep are given, cannot provide internal_base_ring")
            Polyhedron_base_number_field.__init__(self, parent, Vrep, Hrep, **kwds)

    def _nmz_result(self, normaliz_cone, property):
        """
        Call PyNormaliz's NmzResult function with appropriate conversion between number format.

        TESTS::

            sage: p = Polyhedron(vertices=[(0, 0), (1, 0), (0, 1)], rays=[(1,1)],
            ....:                lines=[], backend='normaliz')
            sage: p._nmz_result(p._normaliz_cone, 'EquivariantXyzzyModuleSeries')
            Traceback (most recent call last):
            ...
            NormalizError: Some error in the normaliz input data detected: Unknown ConeProperty...

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 - 3, embedding=AA(3)**(1/3))
            sage: p = Polyhedron(vertices=[(0, 0), (1, 1), (a, 3), (-1, a**2)],
            ....:                rays=[(-1,-a)], backend='normaliz')
            sage: sorted(p._nmz_result(p._normaliz_cone, 'VerticesOfPolyhedron'))
            [[-1, a^2, 1], [1, 1, 1], [a, 3, 1]]
            sage: triangulation_generators = p._nmz_result(p._normaliz_cone,
            ....:                                          'Triangulation')[1]
            sage: sorted(triangulation_generators)
            [[-a^2, -3, 0], [-1, a^2, 1], [0, 0, 1], [1, 1, 1], [a, 3, 1]]
            sage: p._nmz_result(p._normaliz_cone, 'AffineDim') == 2
            True
            sage: p._nmz_result(p._normaliz_cone, 'EmbeddingDim') == 3
            True
            sage: p._nmz_result(p._normaliz_cone, 'ExtremeRays')
            [[-1/3*a^2, -1, 0]]
            sage: p._nmz_result(p._normaliz_cone, 'MaximalSubspace')
            []
        """
        def rational_handler(list):
            return QQ(tuple(list))

        def nfelem_handler(coords):
            # coords might be too short which is not accepted by Sage number field
            v = list(coords) + [0] * (self._internal_base_ring.degree() - len(coords))
            return self._internal_base_ring(v)
        return NmzResult(normaliz_cone, property,
                         RationalHandler=rational_handler,
                         NumberfieldElementHandler=nfelem_handler)

    def _init_from_normaliz_cone(self, normaliz_cone, internal_base_ring):
        """
        Construct polyhedron from a PyNormaliz wrapper of a normaliz cone.

        TESTS::

            sage: p = Polyhedron(backend='normaliz')
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz
            sage: Polyhedron_normaliz._init_from_Hrepresentation(p, [], [])  # indirect doctest
        """
        if internal_base_ring is None:
            internal_base_ring = QQ
        self._internal_base_ring = internal_base_ring

        if normaliz_cone and self._nmz_result(normaliz_cone, "AffineDim") < 0:
            # Empty polyhedron. Special case because Normaliz defines the
            # recession cone of an empty polyhedron given by an
            # H-representation as the cone defined by the homogenized system.
            self._init_empty_polyhedron()
        else:
            self._normaliz_cone = normaliz_cone
            self._init_Vrepresentation_from_normaliz()
            self._init_Hrepresentation_from_normaliz()

    @staticmethod
    def _convert_to_pynormaliz(x):
        """
        Convert a number or nested lists and tuples of numbers to pynormaliz input format.

        TESTS::

            sage: K.<sqrt2> = QuadraticField(2)                                         # needs sage.rings.number_field
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz as Pn
            sage: Pn._convert_to_pynormaliz(17)
            17
            sage: Pn._convert_to_pynormaliz(901824309821093821093812093810928309183091832091)
            901824309821093821093812093810928309183091832091
            sage: Pn._convert_to_pynormaliz(QQ(17))
            17
            sage: Pn._convert_to_pynormaliz(28/5)
            [[28, 5]]
            sage: Pn._convert_to_pynormaliz(28901824309821093821093812093810928309183091832091/5234573685674784567853456543456456786543456765)
            [[28901824309821093821093812093810928309183091832091, 5234573685674784567853456543456456786543456765]]
            sage: Pn._convert_to_pynormaliz(7 + sqrt2)                                  # needs sage.rings.number_field
            [[7, 1], [1, 1]]
            sage: Pn._convert_to_pynormaliz(7/2 + sqrt2)                                # needs sage.rings.number_field
            [[7, 2], [1, 1]]
            sage: Pn._convert_to_pynormaliz([[1, 2], (3, 4)])
            [[1, 2], [3, 4]]

        Check that :issue:`29836` is fixed::

            sage: P = polytopes.simplex(backend='normaliz')
            sage: K.<sqrt2> = QuadraticField(2)                                         # needs sage.rings.number_field
            sage: P.dilation(sqrt2)                                                     # needs sage.rings.number_field
            A 3-dimensional polyhedron in
             (Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.41...)^4
             defined as the convex hull of 4 vertices
        """
        def _QQ_pair(x):
            x = QQ(x)
            return [int(x.numerator()), int(x.denominator())]
        from sage.rings.rational import Rational
        from types import GeneratorType
        if isinstance(x, (list, tuple, GeneratorType)):
            return [Polyhedron_normaliz._convert_to_pynormaliz(y) for y in x]
        try:
            return int(ZZ(x))
        except TypeError:
            pass

        if isinstance(x, Rational):
            return [_QQ_pair(x)]  # need extra brackets to distinguish from quadratic numberfield element
        # number field
        return [_QQ_pair(c) for c in x.list()]

    def _init_from_normaliz_data(self, data, internal_base_ring=None, verbose=False):
        """
        Construct polyhedron from normaliz ``data`` (a dictionary).

        TESTS::

            sage: p = Polyhedron(backend='normaliz', ambient_dim=2)
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_QQ_normaliz
            sage: data = {'inhom_inequalities': [[-1, 2, 0], [0, 0, 1], [2, -1, 0]]}
            sage: Polyhedron_QQ_normaliz._init_from_normaliz_data(p, data)
            sage: p.inequalities_list()
            [[0, -1, 2], [0, 2, -1]]

            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz
            sage: from sage.rings.qqbar import AA                                       # needs sage.rings.number_field
            sage: from sage.rings.number_field.number_field import QuadraticField       # needs sage.rings.number_field
            sage: data = {'number_field': ['a^2 - 2', 'a', '[1.4 +/- 0.1]'],
            ....: 'inhom_inequalities': [[-1, 2, 0], [0, 0, 1], [2, -1, 0]]}
            sage: from sage.geometry.polyhedron.parent import Polyhedra_normaliz
            sage: parent = Polyhedra_normaliz(AA, 2, 'normaliz')                        # needs sage.rings.number_field
            sage: Polyhedron_normaliz(parent, None, None,                               # needs sage.rings.number_field
            ....:                     normaliz_data=data,
            ....:                     internal_base_ring=QuadraticField(2))
            A 2-dimensional polyhedron in AA^2 defined as the convex hull of 1 vertex and 2 rays
            sage: _.inequalities_list()                                                 # needs sage.rings.number_field
            [[0, -1/2, 1], [0, 2, -1]]
        """
        if internal_base_ring is None:
            internal_base_ring = QQ
        cone = self._cone_from_normaliz_data(data, verbose)
        self._init_from_normaliz_cone(cone, internal_base_ring)

    def _cone_from_normaliz_data(self, data, verbose=False):
        """
        Construct a normaliz cone from ``data`` (a dictionary).

        EXAMPLES::

            sage: p = Polyhedron(backend='normaliz', ambient_dim=2)
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_QQ_normaliz
            sage: data = {'inhom_inequalities': [[-1, 2, 0], [0, 0, 1], [2, -1, 0]]}
            sage: cone = Polyhedron_QQ_normaliz._cone_from_normaliz_data(p, data)
            sage: p._nmz_result(cone,'SupportHyperplanes')
            [[-1, 2, 0], [0, 0, 1], [2, -1, 0]]
        """
        if verbose:
            if isinstance(verbose, str):
                print("# Wrote equivalent Normaliz input file to {}".format(verbose))
                self._normaliz_format(data, file_output=verbose)
            else:
                print("# ----8<---- Equivalent Normaliz input file ----8<----")
                print(self._normaliz_format(data), end='')
                print("# ----8<-------------------8<-------------------8<----")

        for key, value in data.items():
            if key != 'number_field':
                data[key] = self._convert_to_pynormaliz(value)

        if verbose:
            print("# Calling {}".format(_format_function_call('PyNormaliz.NmzCone', **data)))

        cone = NmzCone(**data)
        assert cone, "{} did not return a cone".format(_format_function_call('PyNormaliz.NmzCone', **data))
        return cone

    def _is_zero(self, x) -> bool:
        """
        Test whether ``x`` is zero.

        INPUT:

        - ``x`` -- a number in the base ring

        OUTPUT: boolean

        EXAMPLES::

            sage: p = Polyhedron([(sqrt(3),sqrt(2))], base_ring=AA)                     # needs sage.rings.number_field sage.symbolic
            sage: p._is_zero(0)                                                         # needs sage.rings.number_field sage.symbolic
            True
            sage: p._is_zero(1/100000)                                                  # needs sage.rings.number_field sage.symbolic
            False
        """
        return x == 0

    def _is_nonneg(self, x) -> bool:
        """
        Test whether ``x`` is nonnegative.

        INPUT:

        - ``x`` -- a number in the base ring

        OUTPUT: boolean

        EXAMPLES::

            sage: p = Polyhedron([(sqrt(3),sqrt(2))], base_ring=AA)                     # needs sage.rings.number_field sage.symbolic
            sage: p._is_nonneg(1)                                                       # needs sage.rings.number_field sage.symbolic
            True
            sage: p._is_nonneg(-1/100000)                                               # needs sage.rings.number_field sage.symbolic
            False
        """
        return x >= 0

    def _is_positive(self, x) -> bool:
        """
        Test whether ``x`` is positive.

        INPUT:

        - ``x`` -- a number in the base ring

        OUTPUT: boolean

        EXAMPLES::

            sage: p = Polyhedron([(sqrt(3),sqrt(2))], base_ring=AA)                     # needs sage.rings.number_field sage.symbolic
            sage: p._is_positive(1)                                                     # needs sage.rings.number_field sage.symbolic
            True
            sage: p._is_positive(0)                                                     # needs sage.rings.number_field sage.symbolic
            False
        """
        return x > 0

    def _init_from_Vrepresentation(self, vertices, rays, lines, minimize=True, verbose=False):
        r"""
        Construct polyhedron from V-representation data.

        INPUT:

        - ``vertices`` -- list of point; each point can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``rays`` -- list of rays; each ray can be specified as any
          iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``lines`` -- list of lines; each line can be specified as
          any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``verbose`` -- boolean (default: ``False``); whether to print
          verbose output for debugging purposes

        EXAMPLES::

            sage: p = Polyhedron(backend='normaliz')
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz
            sage: Polyhedron_normaliz._init_from_Vrepresentation(p, [], [], [])
        """

        def vert_ray_line_QQ(vertices, rays, lines):
            nmz_vertices = []
            for v in vertices:
                d = LCM_list([denominator(v_i) for v_i in v])
                dv = [d * v_i for v_i in v]
                nmz_vertices.append(dv + [d])
            nmz_rays = []
            for r in rays:
                d = LCM_list([denominator(r_i) for r_i in r])
                dr = [d * r_i for r_i in r]
                nmz_rays.append(dr)
            nmz_lines = []
            for l in lines:
                d = LCM_list([denominator(l_i) for l_i in l])
                dl = [d * l_i for l_i in l]
                nmz_lines.append(dl)
            return nmz_vertices, nmz_rays, nmz_lines

        def vert_ray_line_NF(vertices, rays, lines):
            h_vertices = [list(v) + [1] for v in vertices]
            return h_vertices, rays, lines

        if vertices is None:
            vertices = []
        if rays is None:
            rays = []
        if lines is None:
            lines = []

        (nmz_vertices, nmz_rays, nmz_lines), internal_base_ring \
            = self._compute_data_lists_and_internal_base_ring(
                (vertices, rays, lines), vert_ray_line_QQ, vert_ray_line_NF)

        if not nmz_vertices and not nmz_rays and not nmz_lines:
            # Special case to avoid:
            #   error: Some error in the normaliz input data detected:
            #   All input matrices empty!
            self._init_empty_polyhedron()
        else:
            data = {"vertices": nmz_vertices,
                    "cone": nmz_rays,
                    "subspace": nmz_lines}
            number_field_data = self._number_field_triple(internal_base_ring)
            if number_field_data:
                data["number_field"] = number_field_data
            self._init_from_normaliz_data(data, internal_base_ring=internal_base_ring, verbose=verbose)

    def _init_from_Hrepresentation(self, ieqs, eqns, minimize=True, verbose=False):
        r"""
        Construct polyhedron from H-representation data.

        INPUT:

        - ``ieqs`` -- list of inequalities; each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``eqns`` -- list of equalities; each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``minimize`` -- boolean (default: ``True``); ignored

        - ``verbose`` -- boolean (default: ``False``); whether to print
          verbose output for debugging purposes

        EXAMPLES::

            sage: p = Polyhedron(backend='normaliz')
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz
            sage: Polyhedron_normaliz._init_from_Hrepresentation(p, [], [])

        TESTS::

            sage: K.<a> = QuadraticField(2)                                             # needs sage.rings.number_field
            sage: p = Polyhedron(ieqs=[(1, a, 0)], backend='normaliz')                  # needs sage.rings.number_field
            sage: p & p == p                                                            # needs sage.rings.number_field
            True

        Check that :issue:`30248` is fixed, that maps as input works::

            sage: # needs sage.rings.number_field
            sage: q = Polyhedron(backend='normaliz', base_ring=AA,
            ....:                rays=[(0, 0, 1), (0, 1, -1), (1, 0, -1)])
            sage: def make_new_Hrep(h):
            ....:     return tuple(x if i == 0 else -1*x
            ....:                  for i, x in enumerate(h._vector))
            sage: new_inequalities = map(make_new_Hrep, q.inequality_generator())
            sage: new_equations = map(make_new_Hrep, q.equation_generator())
            sage: parent = q.parent()
            sage: new_q = parent.element_class(parent, None,
            ....:                              [new_inequalities, new_equations])
            sage: new_q
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 1 vertex and 3 rays
        """

        def nmz_ieqs_eqns_NF(ieqs, eqns):
            nmz_ieqs = [list(ieq[1:]) + [ieq[0]] for ieq in ieqs]
            nmz_eqns = [list(eqn[1:]) + [eqn[0]] for eqn in eqns]
            return nmz_ieqs, nmz_eqns

        def nmz_ieqs_eqns_QQ(ieqs, eqns):
            nmz_ieqs = []
            for ieq in ieqs:
                d = LCM_list([denominator(ieq_i) for ieq_i in ieq])
                dieq = [ZZ(d * ieq_i) for ieq_i in ieq]
                b = dieq[0]
                A = dieq[1:]
                nmz_ieqs.append(A + [b])
            nmz_eqns = []
            for eqn in eqns:
                d = LCM_list([denominator(eqn_i) for eqn_i in eqn])
                deqn = [ZZ(d * eqn_i) for eqn_i in eqn]
                b = deqn[0]
                A = deqn[1:]
                nmz_eqns.append(A + [b])
            return nmz_ieqs, nmz_eqns

        if ieqs is None:
            ieqs = []
        if eqns is None:
            eqns = []

        (nmz_ieqs, nmz_eqns), internal_base_ring \
            = self._compute_data_lists_and_internal_base_ring(
                (ieqs, eqns), nmz_ieqs_eqns_QQ, nmz_ieqs_eqns_NF)
        if not nmz_ieqs:
            # If normaliz gets an empty list of inequalities, it adds
            # nonnegativities. So let's add a tautological inequality to work
            # around this.
            nmz_ieqs.append([0] * self.ambient_dim() + [0])
        data = {"inhom_equations": nmz_eqns,
                "inhom_inequalities": nmz_ieqs}
        number_field_data = self._number_field_triple(internal_base_ring)
        if number_field_data:
            data["number_field"] = number_field_data
        self._init_from_normaliz_data(data, internal_base_ring=internal_base_ring, verbose=verbose)

    def _cone_from_Vrepresentation_and_Hrepresentation(self, vertices, rays, lines, ieqs, eqns=None, verbose=False, homogeneous=False):
        r"""
        Construct cone from V-representation data and H-representation data.

        INPUT:

        - ``vertices`` -- list of point; each point can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``rays`` -- list of rays; each ray can be specified as any
          iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``ieqs`` -- list of inequalities; each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``eqns`` -- list of equalities; each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``verbose`` -- boolean (default: ``False``); whether to print
          verbose output for debugging purposes

        - ``homogeneous`` -- boolean (default: ``False``); if ``True`` set
          up the cone without explicit inhomogenization

        EXAMPLES::

            sage: P = (polytopes.hypercube(4, backend='normaliz')
            ....:       * Polyhedron(rays=[[0,1]])
            ....:       * Polyhedron(lines=[[1,0]])); P
            A 6-dimensional polyhedron in ZZ^8 defined as the convex hull of 16 vertices, 1 ray, 1 line

            sage: cone = P._cone_from_Vrepresentation_and_Hrepresentation(
            ....:     P.vertices(), P.rays(), P.lines(),
            ....:     P.inequalities(), P.equations())
            sage: import PyNormaliz
            sage: PyNormaliz.NmzIsComputed(cone, "VerticesOfPolyhedron")
            True
            sage: PyNormaliz.NmzIsComputed(cone, "ExtremeRays")
            True
            sage: PyNormaliz.NmzIsComputed(cone, "MaximalSubspace")
            True
            sage: PyNormaliz.NmzIsComputed(cone, "SupportHyperplanes")
            True
            sage: PyNormaliz.NmzIsComputed(cone, "Equations")
            False

        All values must be specified::

            sage: cone = P._cone_from_Vrepresentation_and_Hrepresentation(
            ....:     P.vertices(), None, P.lines(),
            ....:     P.inequalities(), P.equations())
            Traceback (most recent call last):
            ...
            ValueError: please specify vertices, rays, lines, inequalities and equations completely

        This method cannot be used for the empty cone::

            sage: P = Polyhedron(backend='normaliz')
            sage: cone = P._cone_from_Vrepresentation_and_Hrepresentation(
            ....:     P.vertices(), P.rays(), P.lines(),
            ....:     P.inequalities(), P.equations())
            Traceback (most recent call last):
            ...
            ValueError: this method cannot be used to initialize the empty cone

        TESTS::

            sage: def test_poly(P):
            ....:     cone = P._cone_from_Vrepresentation_and_Hrepresentation(P.vertices(),P.rays(),P.lines(),P.inequalities(),P.equations())
            ....:     cone2 = P._normaliz_cone
            ....:     args = ['Equations','VerticesOfPolyhedron','ExtremeRays','SupportHyperplanes','MaximalSubspace']
            ....:     return all(P._nmz_result(cone,arg) == P._nmz_result(cone2, arg)
            ....:                for arg in args)
            sage: test_poly(polytopes.simplex(backend='normaliz'))
            True
            sage: test_poly(polytopes.dodecahedron(backend='normaliz'))                 # needs sage.rings.number_field
            True
            sage: test_poly(Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,1]],
            ....:                      backend='normaliz'))
            True
            sage: test_poly(Polyhedron(vertices=[[-1,0], [1,0]], lines=[[0,1]],
            ....:                      backend='normaliz'))
            True
            sage: test_poly(Polyhedron(rays=[[1,0,0],[0,1,0]],
            ....:                      backend='normaliz'))
            True
            sage: test_poly(Polyhedron(vertices=[[1,0,0], [0,1,0]], rays=[[1,0,0], [0,1,0]],
            ....:                      backend='normaliz'))
            True
            sage: test_poly(Polyhedron(vertices=[[0,0,0], [0,1,1], [1,0,1], [-1,-1,1]],
            ....:                      rays=[[0,0,1]],
            ....:                      backend='normaliz'))
            True

        Old input format will give a meaningful error message::

            sage: cone = P._cone_from_Vrepresentation_and_Hrepresentation(
            ....:     P.vertices(), P.rays(),
            ....:     P.inequalities(), P.equations())
            Traceback (most recent call last):
            ...
            ValueError: the specification of this method has changed; please specify the lines as well

            sage: cone = P._cone_from_Vrepresentation_and_Hrepresentation(
            ....:     P.vertices(), P.rays(),
            ....:     P.inequalities(), P.equations(), True)
            Traceback (most recent call last):
            ...
            ValueError: the specification of this method has changed; please specify the lines as well

        Check that :issue:`30891` is fixed::

            sage: p = Polyhedron(vertices=[(-3,-3), (3,0), (3,3), (0,3)],
            ....:                backend='normaliz')
            sage: q = loads(p.dumps())
            sage: q.volume()
            18
            sage: q.ehrhart_series()
            (13*t^2 + 22*t + 1)/(-t^3 + 3*t^2 - 3*t + 1)
        """
        if eqns in (True, False, None):
            # Previously, the method had input ``vertices, rays, ieqs, eqns`` (optionally ``verbose``).
            # Now it requires ``vertices, rays, lines, ieqs, eqns``.
            # Actually, ``eqns`` wouldn't be required, but we keep it to catch deprecated calls.
            # (And it's more stable against changes of normaliz now.)
            raise ValueError("the specification of this method has changed; please specify the lines as well")
        if None in (vertices, rays, lines, ieqs, eqns):
            raise ValueError("please specify vertices, rays, lines, inequalities and equations completely")
        if not vertices:
            raise ValueError("this method cannot be used to initialize the empty cone")

        def rays_subspace_lattice_ieqs_QQ(vertices, rays, lines, ieqs):
            nmz_vertices = []
            for v in vertices:
                d = LCM_list([denominator(v_i) for v_i in v])
                dv = [d * v_i for v_i in v]
                nmz_vertices.append(dv + [d])
            nmz_rays = []
            for r in rays:
                d = LCM_list([denominator(r_i) for r_i in r])
                dr = [d * r_i for r_i in r]
                nmz_rays.append(dr + [0])
            nmz_lines = []
            for l in lines:
                d = LCM_list([denominator(l_i) for l_i in l])
                dl = [d * l_i for l_i in l]
                nmz_lines.append(dl + [0])

            nmz_ieqs = []
            for ieq in ieqs:
                d = LCM_list([denominator(ieq_i) for ieq_i in ieq])
                dieq = [ZZ(d * ieq_i) for ieq_i in ieq]
                b = dieq[0]
                A = dieq[1:]
                nmz_ieqs.append(A + [b])

            from sage.matrix.constructor import Matrix
            lattice = Matrix(ZZ, nmz_vertices + nmz_rays + nmz_lines).saturation()
            nmz_lattice = [list(y) for y in lattice]

            if Matrix(ZZ, nmz_vertices + nmz_rays).rank() == Matrix(ZZ, nmz_rays).rank() + 1:
                # The recession cone is full-dimensional.
                # In this case the homogenized inequalities
                # do not ensure nonnegativy in the last coordinate.
                # In the homogeneous cone the far face is a facet.
                pos_ieq = [ZZ.zero()] * len(nmz_vertices[0])
                pos_ieq[-1] = ZZ.one()
                nmz_ieqs.append(pos_ieq)

            return nmz_vertices + nmz_rays, nmz_lines, nmz_lattice, nmz_ieqs

        def rays_subspace_lattice_ieqs_NF(vertices, rays, lines, ieqs):
            nmz_vertices = [list(v) + [1] for v in vertices]
            nmz_rays = [list(r) + [0] for r in rays]
            nmz_lines = [list(l) + [1] for l in lines]

            nmz_ieqs = []
            for ieq in ieqs:
                b = ieq[0]
                A = ieq[1:]
                nmz_ieqs.append(list(A) + [b])

            from sage.matrix.constructor import Matrix
            lattice = Matrix(nmz_vertices + nmz_rays + nmz_lines).row_space().basis()
            nmz_lattice = [list(y) for y in lattice]

            if Matrix(nmz_vertices + nmz_rays).rank() == Matrix(nmz_rays).rank() + 1:
                # The recession cone is full-dimensional.
                # In this case the homogenized inequalities
                # do not ensure nonnegativy in the last coordinate.
                # In the homogeneous cone the far face is a facet.
                pos_ieq = [0] * len(nmz_vertices[0])
                pos_ieq[-1] = 1
                nmz_ieqs.append(pos_ieq)

            return nmz_vertices + nmz_rays, nmz_lines, nmz_lattice, nmz_ieqs

        (nmz_extreme_rays, nmz_subspace, nmz_lattice, nmz_ieqs), internal_base_ring \
            = self._compute_data_lists_and_internal_base_ring(
                (vertices, rays, lines, ieqs), rays_subspace_lattice_ieqs_QQ,
                rays_subspace_lattice_ieqs_NF)

        data = {"extreme_rays": nmz_extreme_rays,
                "maximal_subspace": nmz_subspace,
                "generated_lattice": nmz_lattice,
                "support_hyperplanes": nmz_ieqs}

        ambient_dim = len(data["extreme_rays"][0])
        if not homogeneous:
            data["dehomogenization"] = [[0] * (ambient_dim - 1) + [1]]

        number_field_data = self._number_field_triple(internal_base_ring)
        if number_field_data:
            data["number_field"] = number_field_data
        return self._cone_from_normaliz_data(data, verbose=verbose)

    def _test_far_facet_condition(self, tester=None, **options):
        """
        Test that we add an extra inequality in the correct cases.

        TESTS::

            sage: P = Polyhedron(rays=[[1,1]], backend='normaliz')
            sage: P._test_far_facet_condition()

            sage: P = Polyhedron(vertices=[[1,0], [0,1]],
            ....:                rays=[[1,1]], backend='normaliz')
            sage: P._test_far_facet_condition()

            sage: P = Polyhedron(rays=[[1,1,0]],
            ....:                lines=[[0,0,1]], backend='normaliz')
            sage: P._test_far_facet_condition()

            sage: P = Polyhedron(vertices=[[1,0,0], [0,1,0]],
            ....:                rays=[[1,1,0]],
            ....:                lines=[[0,0,1]], backend='normaliz')
            sage: P._test_far_facet_condition()
        """
        if tester is None:
            tester = self._tester(**options)

        if self.is_empty():
            return

        nmz_vertices = self._nmz_result(self._normaliz_cone, "VerticesOfPolyhedron")
        nmz_rays = self._nmz_result(self._normaliz_cone, "ExtremeRays")
        nmz_ieqs = self._nmz_result(self._normaliz_cone, "SupportHyperplanes")

        from sage.matrix.constructor import Matrix
        far_facet_condition = Matrix(nmz_vertices + nmz_rays).rank() == Matrix(nmz_rays).rank() + 1

        tester.assertEqual(far_facet_condition, self.n_inequalities() != len(nmz_ieqs))

        if far_facet_condition:
            tester.assertEqual(self.n_inequalities() + 1, len(nmz_ieqs))
            tester.assertTrue(any(ieq == [0] * self.ambient_dim() + [1] for ieq in nmz_ieqs))

    def _init_Vrepresentation_from_normaliz(self):
        r"""
        Create the Vrepresentation objects from the normaliz polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2), (2,0), (4,5/6)],  # indirect doctest
            ....:                backend='normaliz')
            sage: p.Hrepresentation()
            (An inequality (-5, 12) x + 10 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (1, 4) x - 2 >= 0)
            sage: p.Vrepresentation()
            (A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6))
        """
        self._Vrepresentation = []
        parent = self.parent()
        base_ring = self.base_ring()
        cone = self._normaliz_cone
        for g in self._nmz_result(cone, "VerticesOfPolyhedron"):
            d = g[-1]
            if d == 1:
                parent._make_Vertex(self, g[:-1])
            else:
                parent._make_Vertex(self, [base_ring(x) / d for x in g[:-1]])
        for g in self._nmz_result(cone, "ExtremeRays"):
            parent._make_Ray(self, g[:-1])
        for g in self._nmz_result(cone, "MaximalSubspace"):
            parent._make_Line(self, g[:-1])
        self._Vrepresentation = tuple(self._Vrepresentation)

    def _init_Hrepresentation_from_normaliz(self):
        r"""
        Create the Hrepresentation objects from the normaliz polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2), (2,0), (4,5/6)],  # indirect doctest
            ....:                backend='normaliz')
            sage: p.Hrepresentation()
            (An inequality (-5, 12) x + 10 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (1, 4) x - 2 >= 0)
            sage: p.Vrepresentation()
            (A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6))
        """
        self._Hrepresentation = []
        cone = self._normaliz_cone
        parent = self.parent()
        for g in self._nmz_result(cone, "SupportHyperplanes"):
            if all(x == 0 for x in g[:-1]):
                # Ignore vertical inequality
                pass
            else:
                parent._make_Inequality(self, (g[-1],) + tuple(g[:-1]))
        for g in self._nmz_result(cone, "Equations"):
            parent._make_Equation(self, (g[-1],) + tuple(g[:-1]))
        self._Hrepresentation = tuple(self._Hrepresentation)

    def _init_empty_polyhedron(self):
        r"""
        Initialize an empty polyhedron.

        TESTS::

            sage: empty = Polyhedron(backend='normaliz'); empty
            The empty polyhedron in ZZ^0
            sage: empty.Vrepresentation()
            ()
            sage: empty.Hrepresentation()
            (An equation -1 == 0,)
            sage: Polyhedron(vertices=[], backend='normaliz')
            The empty polyhedron in ZZ^0
            sage: Polyhedron(backend='normaliz')._init_empty_polyhedron()
        """
        super()._init_empty_polyhedron()
        # Can't seem to set up an empty _normaliz_cone.
        # For example, PyNormaliz.NmzCone(vertices=[]) gives
        # error: Some error in the normaliz input data detected: All input matrices empty!
        self._normaliz_cone = None

    @classmethod
    def _from_normaliz_cone(cls, parent, normaliz_cone, internal_base_ring=None):
        r"""
        Initialize a polyhedron from a PyNormaliz wrapper of a normaliz cone.

        TESTS::

            sage: P=Polyhedron(ieqs=[[1, 0, 2], [3, 0, -2], [3, 2, -2]],
            ....:              backend='normaliz')
            sage: PI = P.integral_hull()                   # indirect doctest
        """
        return cls(parent, None, None, normaliz_cone=normaliz_cone, internal_base_ring=internal_base_ring)

    @staticmethod
    def _number_field_triple(internal_base_ring) -> list:
        r"""
        Construct the PyNormaliz triple that describes ``internal_base_ring``.

        TESTS::

            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz as Pn
            sage: Pn._number_field_triple(QQ) is None
            True
            sage: Pn._number_field_triple(QuadraticField(5))                            # needs sage.rings.number_field
            ['a^2 - 5', 'a', '[2.236067977499789 +/- 8.06e-16]']
        """
        R = internal_base_ring
        if R is QQ:
            return None
        from sage.rings.real_arb import RealBallField
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        emb = RealBallField(53)(R.gen(0))
        gen = 'a'
        R_a = PolynomialRing(QQ, gen)
        min_poly = R_a(R.polynomial())
        return [str(min_poly), gen, str(emb)]

    @staticmethod
    def _make_normaliz_cone(data, verbose=False):
        r"""
        Return a normaliz cone from ``data``.

        INPUT:

        - ``data`` -- dictionary

        - ``verbose`` -- boolean (default: ``False``)

        TESTS::

            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz
            sage: data = {'inhom_inequalities': [[-1, 2, 0], [0, 0, 1], [2, -1, 0]]}
            sage: nmz_cone = Polyhedron_normaliz._make_normaliz_cone(data,verbose=False)
            sage: from PyNormaliz import NmzResult
            sage: NmzResult(nmz_cone, "ExtremeRays")
            [[1, 2, 0], [2, 1, 0]]
        """
        if verbose:
            print("# Calling PyNormaliz.NmzCone(**{})".format(data))
        cone = NmzCone(**data)
        assert cone, "NmzCone(**{}) did not return a cone".format(data)
        return cone

    def _get_nmzcone_data(self) -> dict:
        r"""
        Get the data necessary to reproduce the normaliz cone.

        OUTPUT: ``data`` -- dictionary

        TESTS:

        The empty polyhedron::

            sage: P = Polyhedron(backend='normaliz')
            sage: P._get_nmzcone_data()
            {}

        Another simple example::

            sage: C = Polyhedron(backend='normaliz', rays=[[1, 2], [2, 1]])
            sage: C._get_nmzcone_data()
            {'cone': [[1, 2], [2, 1]],
             'inhom_equations': [],
             'inhom_inequalities': [[-1, 2, 0], [0, 0, 1], [2, -1, 0]],
             'subspace': [],
             'vertices': [[0, 0, 1]]}
        """
        if self.is_empty():
            return {}

        vertices = self._nmz_result(self._normaliz_cone, "VerticesOfPolyhedron")
        # get rid of the last 0 in rays:
        rays = [r[:-1] for r in self._nmz_result(self._normaliz_cone, "ExtremeRays")]
        lines = self._nmz_result(self._normaliz_cone, "MaximalSubspace")
        ineqs = self._nmz_result(self._normaliz_cone, "SupportHyperplanes")
        eqs = self._nmz_result(self._normaliz_cone, "Equations")

        return {'vertices': vertices,
                'cone': rays,
                'subspace': lines,
                'inhom_equations': eqs,
                'inhom_inequalities': ineqs}

    def _normaliz_format(self, data, file_output=None):
        r"""
        Return a string containing normaliz format.

        INPUT:

        - ``data`` -- dictionary of PyNormaliz cone input properties

        - ``file_output`` -- string (optional); a filename to which the
          representation should be written. If set to ``None`` (default),
          representation is returned as a string.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[0, 0], [0, 1], [1, 0]],  # indirect doctest
            ....:                backend='normaliz', verbose=True)
            # ----8<---- Equivalent Normaliz input file ----8<----
            amb_space 2
            cone 0
            subspace 0
            vertices 3
             0 0 1
             0 1 1
             1 0 1
            # ----8<-------------------8<-------------------8<----
            # Calling ...
        """
        def format_number(x):
            try:
                return '{}'.format(QQ(x))
            except (ValueError, TypeError):
                return '({})'.format(x.polynomial('a'))

        def format_field(key, value):
            if isinstance(value, (list, tuple)):
                s = '{} {}\n'.format(key, len(value))
                for e in value:
                    for x in e:
                        s += ' ' + format_number(x)
                    s += '\n'
                return s
            return '{} {}\n'.format(key, value)

        def format_number_field_data(nf_triple):
            min_poly, gen, emb = nf_triple
            return 'min_poly ({}) embedding {}'.format(min_poly, emb)

        s = format_field('amb_space', self.ambient_dim())
        if 'number_field' in data:
            from copy import copy
            data = copy(data)
            s += 'number_field {}\n'.format(format_number_field_data(data['number_field']))
            del data['number_field']
        for key, value in sorted(data.items()):
            s += format_field(key, value)
        if file_output is not None:
            with open(file_output, 'w') as in_file:
                in_file.write(s)
        else:
            return s

    def __copy__(self):
        r"""
        Return a copy of ``self``.

        TESTS::

            sage: P = polytopes.cube(backend='normaliz')
            sage: Q = copy(P)
            sage: P._normaliz_cone is Q._normaliz_cone
            False
        """
        other = super().__copy__()

        # Make a copy of the cone.
        cone = self._normaliz_cone
        conecopy = NmzConeCopy(cone)
        other._normaliz_cone = conecopy
        return other

    def __getstate__(self):
        r"""
        Remove the normaliz cone for pickling.

        TESTS::

            sage: P = polytopes.simplex(backend='normaliz')
            sage: P.__getstate__()
            (Polyhedra in ZZ^4,
             {'_Hrepresentation': (An inequality (0, 0, 0, 1) x + 0 >= 0,
               An inequality (0, 0, 1, 0) x + 0 >= 0,
               An inequality (0, 1, 0, 0) x + 0 >= 0,
               An inequality (1, 0, 0, 0) x + 0 >= 0,
               An equation (1, 1, 1, 1) x - 1 == 0),
              '_Vrepresentation': (A vertex at (0, 0, 0, 1),
               A vertex at (0, 0, 1, 0),
               A vertex at (0, 1, 0, 0),
               A vertex at (1, 0, 0, 0)),
              '_internal_base_ring': Rational Field,
              '_pickle_equations': [(-1, 1, 1, 1, 1)],
              '_pickle_inequalities': [(0, 0, 0, 0, 1),
               (0, 0, 0, 1, 0),
               (0, 0, 1, 0, 0),
               (0, 1, 0, 0, 0)],
              '_pickle_lines': [],
              '_pickle_rays': [],
              '_pickle_vertices': [(0, 0, 0, 1),
               (0, 0, 1, 0),
               (0, 1, 0, 0),
               (1, 0, 0, 0)]})
        """
        state = super().__getstate__()
        state = (state[0], state[1].copy())
        # Remove the unpicklable entries.
        del state[1]['_normaliz_cone']
        state[1]["_pickle_vertices"] = [v._vector for v in self.vertices()]
        state[1]["_pickle_rays"] = [v._vector for v in self.rays()]
        state[1]["_pickle_lines"] = [v._vector for v in self.lines()]
        state[1]["_pickle_inequalities"] = [v._vector for v in self.inequalities()]
        state[1]["_pickle_equations"] = [v._vector for v in self.equations()]
        return state

    def __setstate__(self, state):
        r"""
        Initialize the normaliz cone after pickling.

        TESTS::

            sage: P = polytopes.permutahedron(4, backend='normaliz')
            sage: P.volume(measure='induced_lattice', engine='normaliz')
            96
            sage: P.volume.clear_cache()
            sage: P1 = loads(dumps(P))                 # indirect doctest
            sage: P1.volume(measure='induced_lattice', engine='normaliz')
            96

        Test that the obtained cone is valid::

            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz
            sage: P = polytopes.permutahedron(4, backend='normaliz')
            sage: P1 = loads(dumps(P))
            sage: P2 = Polyhedron_normaliz(P1.parent(), None, None, P1._normaliz_cone)
            sage: P2 == P
            True

            sage: P = Polyhedron(lines=[[1,0], [0,1]], backend='normaliz')
            sage: P1 = loads(dumps(P))
            sage: P2 = Polyhedron_normaliz(P1.parent(), None, None, P1._normaliz_cone)
            sage: P2 == P
            True

            sage: P = Polyhedron(backend='normaliz')
            sage: P1 = loads(dumps(P))
            sage: P2 = Polyhedron_normaliz(P1.parent(), None, None, P1._normaliz_cone)
            sage: P2 == P
            True

            sage: P = polytopes.permutahedron(4, backend='normaliz') * Polyhedron(lines=[[1]], backend='normaliz')
            sage: P1 = loads(dumps(P))
            sage: P2 = Polyhedron_normaliz(P1.parent(), None, None, P1._normaliz_cone)
            sage: P2 == P
            True

            sage: # needs sage.rings.number_field
            sage: P = polytopes.dodecahedron(backend='normaliz')
            sage: P1 = loads(dumps(P))
            sage: P2 = Polyhedron_normaliz(P1.parent(), None, None, P1._normaliz_cone,
            ....:                          internal_base_ring=P1._internal_base_ring)
            sage: P == P2
            True

        Test that :issue:`31820` is fixed::

            sage: P = polytopes.cube(backend='normaliz')
            sage: v = P.Vrepresentation()[0]
            sage: v1 = loads(v.dumps())
        """
        if "_pickle_vertices" in state[1]:
            vertices = state[1].pop("_pickle_vertices")
            rays = state[1].pop("_pickle_rays")
            lines = state[1].pop("_pickle_lines")
            inequalities = state[1].pop("_pickle_inequalities")
            equations = state[1].pop("_pickle_equations")
        else:
            vertices = None

        super().__setstate__(state)

        if self.is_empty():
            # Special case to avoid.
            self._normaliz_cone = None
            return

        if vertices is None:
            vertices = self.vertices()
            rays = self.rays()
            lines = self.lines()
            inequalities = self.inequalities()
            equations = self.equations()

        self._normaliz_cone = \
            self._cone_from_Vrepresentation_and_Hrepresentation(
                vertices, rays, lines, inequalities, equations)

    def integral_hull(self):
        r"""
        Return the integral hull in the polyhedron.

        This is a new polyhedron that is the convex hull of all integral
        points.

        EXAMPLES:

        Unbounded example from Normaliz manual, "a dull polyhedron"::

            sage: P = Polyhedron(ieqs=[[1, 0, 2], [3, 0, -2], [3, 2, -2]],
            ....:                backend='normaliz')
            sage: PI = P.integral_hull()
            sage: P.plot(color='yellow') + PI.plot(color='green')                       # needs sage.plot
            Graphics object consisting of 10 graphics primitives
            sage: PI.Vrepresentation()
            (A vertex at (-1, 0),
             A vertex at (0, 1),
             A ray in the direction (1, 0))

        Nonpointed case::

            sage: P = Polyhedron(vertices=[[1/2, 1/3]], rays=[[1, 1]],
            ....:                lines=[[-1, 1]], backend='normaliz')
            sage: PI = P.integral_hull()
            sage: PI.Vrepresentation()
            (A vertex at (1, 0),
             A ray in the direction (1, 0),
             A line in the direction (1, -1))

        Empty polyhedron::

            sage: P = Polyhedron(backend='normaliz')
            sage: PI = P.integral_hull()
            sage: PI.Vrepresentation()
            ()
        """
        if self.is_empty():
            return self
        cone = self._nmz_result(self._normaliz_cone, "IntegerHull")
        return self.parent().element_class._from_normaliz_cone(parent=self.parent(),
                                                               normaliz_cone=cone)

    def _h_star_vector_normaliz(self) -> list:
        r"""
        Return the `h^*`-vector of the lattice polytope.

        INPUT:

        - ``self`` -- a lattice polytope with backend ``'normaliz'``

        OUTPUT:

        The `h^*`-vector as a list.

        EXAMPLES:

        The `h^*`-vector of a unimodular simplex is 1::

            sage: s3 = polytopes.simplex(3, backend='normaliz')
            sage: s3._h_star_vector_normaliz()
            [1]

        The `h^*`-vector of the `0/1`-cube is [1,4,1]::

            sage: cube = polytopes.cube(intervals='zero_one', backend='normaliz')
            sage: cube.h_star_vector()
            [1, 4, 1]

        TESTS:

        Check that :issue:`33847` is fixed::

            sage: L = [[1, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0],
            ....:      [1, 0, 0, 1, 0, 0], [1, 0, 0, 0, 1, 0], [1, 0, 0, 1, 2, 3]]
            sage: P = Polyhedron(vertices=L, backend='normaliz')
            sage: P.h_star_vector()
            [1, 0, 2]
        """
        return self.ehrhart_series().numerator().list()

    def _volume_normaliz(self, measure='euclidean'):
        r"""
        Compute the volume of a polytope using normaliz.

        INPUT:

        - ``measure`` -- string. The measure to use. Allowed values are:

          * ``'euclidean'`` (default): corresponds to ``'EuclideanVolume`` in normaliz
          * ``'induced_lattice'``: corresponds to ``'Volume'`` in normaliz
          * ``'ambient'``: Lebesgue measure of ambient space (volume)

        OUTPUT:

        A float value (when ``measure`` is 'euclidean'),
        a rational number (when ``measure`` is 'induced_lattice'),
        a rational number or symbolic number otherwise (dependent on base ring).

        .. NOTE::

            This function depends on Normaliz (i.e., the ``pynormaliz`` optional
            package).

        REFERENCES:

        See section 6.1.1 of [NormalizMan]_.

        EXAMPLES:

        For normaliz, the default is the euclidean volume in the ambient
        space and the result is a float::

            sage: s = polytopes.simplex(3, backend='normaliz')
            sage: s._volume_normaliz()
            0.3333333333333333

        One other possibility is to compute the scaled volume where a unimodular
        simplex has volume 1::

            sage: s._volume_normaliz(measure='induced_lattice')
            1
            sage: v = [[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
            sage: cube = Polyhedron(vertices=v, backend='normaliz')
            sage: cube._volume_normaliz()
            1.0
            sage: cube._volume_normaliz(measure='induced_lattice')
            6

        Or one can calculate the ambient volume, which is the above multiplied by the
        volume of the unimodular simplex (or zero if not full-dimensional)::

            sage: cube._volume_normaliz(measure='ambient')
            1
            sage: s._volume_normaliz(measure='ambient')
            0

        TESTS:

        Check that :issue:`28872` is fixed::

            sage: P = polytopes.dodecahedron(backend='normaliz')                        # needs sage.rings.number_field
            sage: P.volume(measure='induced_lattice')                                   # needs sage.rings.number_field
            -1056*sqrt5 + 2400

        Some sanity checks that the ambient volume works correctly::

            sage: (2*cube)._volume_normaliz(measure='ambient')
            8
            sage: (1/2*cube)._volume_normaliz(measure='ambient')
            1/8
            sage: s._volume_normaliz(measure='ambient')
            0

            sage: P = polytopes.regular_polygon(3, backend='normaliz')                  # needs sage.rings.number_field
            sage: P._volume_normaliz('ambient') == P.volume(engine='internal')          # needs sage.rings.number_field
            True

            sage: P = polytopes.dodecahedron(backend='normaliz')                        # needs sage.rings.number_field
            sage: P._volume_normaliz('ambient') == P.volume(engine='internal')          # needs sage.rings.number_field
            True

            sage: P = Polyhedron(rays=[[1]], backend='normaliz')
            sage: P.volume()
            +Infinity
        """
        cone = self._normaliz_cone
        assert cone
        if measure == 'euclidean':
            return self._nmz_result(cone, 'EuclideanVolume')

        if measure == 'induced_lattice':
            if self._internal_base_ring in (ZZ, QQ):
                return self._nmz_result(cone, 'Volume')
            return self._nmz_result(cone, 'RenfVolume')

        if measure == 'ambient':
            if self.dim() < self.ambient_dim():
                return self.base_ring().zero()
            if not self.is_compact():
                from sage.rings.infinity import infinity
                return infinity

            from sage.arith.misc import factorial
            return self._volume_normaliz('induced_lattice') / factorial(self.dim())

        raise TypeError("the measure should be `ambient`, `euclidean`, or `induced_lattice`")

    def _triangulate_normaliz(self):
        r"""
        Give a triangulation of the polyhedron using normaliz.

        OUTPUT:

        For compact polyhedra a list of simplices
        each represented by indices of their vertices.

        For cones a list of simplicial cones
        each represented by indices of their rays.

        .. NOTE::

            This function depends on Normaliz (i.e. the ``pynormaliz`` optional
            package). See the Normaliz documentation for further details.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[0,0,1], [1,0,1], [0,1,1], [1,1,1]], backend='normaliz')
            sage: P._triangulate_normaliz()
            [(0, 1, 2), (1, 2, 3)]
            sage: C1 = Polyhedron(rays=[[0,0,1], [1,0,1], [0,1,1], [1,1,1]], backend='normaliz')
            sage: C1._triangulate_normaliz()
            [(0, 1, 2), (1, 2, 3)]
            sage: C2 = Polyhedron(rays=[[1,0,1], [0,0,1], [0,1,1], [1,1,10/9]], backend='normaliz')
            sage: C2._triangulate_normaliz()
            [(0, 1, 2), (1, 2, 3)]

        Works only for cones and compact polyhedra::

            sage: P = polytopes.cube(backend='normaliz')
            sage: Q = Polyhedron(rays=[[0,1]], backend='normaliz')
            sage: R = Polyhedron(lines=[[0,1]], backend='normaliz')
            sage: (P*Q)._triangulate_normaliz()
            Traceback (most recent call last):
            ...
            NotImplementedError: triangulation of non-compact polyhedra that are not cones is not supported
            sage: (P*R)._triangulate_normaliz()
            Traceback (most recent call last):
            ...
            NotImplementedError: triangulation of non-compact not pointed polyhedron is not supported

        TESTS:

        Check that :issue:`30531` is fixed::

            sage: P = polytopes.cube(backend='normaliz')*AA(2).sqrt()
            sage: P._triangulate_normaliz()
            [(0, 1, 2, 4),
            (1, 2, 4, 3),
            (1, 3, 4, 5),
            (3, 5, 6, 7),
            (6, 2, 4, 3),
            (6, 3, 4, 5)]

        ::

            sage: C1 = Polyhedron(rays=[[0,0,1], [1,0,AA(2).sqrt()], [0,1,1], [1,1,1]],
            ....:                 backend='normaliz')
            sage: C1._triangulate_normaliz()
            [(0, 1, 3), (0, 3, 2)]
        """
        if self.lines():
            raise NotImplementedError("triangulation of non-compact not pointed polyhedron is not supported")
        if len(self.vertices_list()) >= 2 and self.rays_list():  # A mix of polytope and cone
            raise NotImplementedError("triangulation of non-compact polyhedra that are not cones is not supported")

        if self.is_compact():
            cone = self._normaliz_cone
        else:
            # Make a inhomogeneous copy of the cone.
            cone = self._cone_from_Vrepresentation_and_Hrepresentation(
                self.vertices(), self.rays(), self.lines(),
                self.inequalities(), self.equations(), homogeneous=True)

        # Compute the triangulation.
        assert cone

        # Normaliz does not guarantee that the order of generators is kept during
        # computation of the triangulation.
        # Those are the generators that the indices of the triangulation correspond to:
        nmz_triangulation, nmz_triangulation_generators = self._nmz_result(cone, "Triangulation")

        base_ring = self.base_ring()
        v_list = self.vertices_list()
        r_list = self.rays_list()

        new_to_old = {}
        for i, g in enumerate(nmz_triangulation_generators):
            if self.is_compact():
                d = base_ring(g[-1])
                vertex = [base_ring(x) / d for x in g[:-1]]
                new_to_old[i] = v_list.index(vertex)
            else:
                if g[-1] > 0:
                    new_to_old[i] = None
                else:
                    try:
                        new_to_old[i] = r_list.index([base_ring(x) for x in g[:-1]])
                    except ValueError:
                        # Rays are only unique up to scaling.
                        new_ray = vector(base_ring, g[:-1])

                        for j, r in enumerate(self.rays()):
                            ray = r.vector()
                            try:
                                # Check for colinearity.
                                _ = new_ray / ray
                                new_to_old[i] = j
                                break
                            except (TypeError, ArithmeticError):
                                pass
                        else:
                            raise ValueError("could not match rays after computing triangulation with original rays")

        def new_indices(old_indices):
            for i in old_indices:
                if new_to_old[i] is not None:
                    yield new_to_old[i]

        return [tuple(new_indices(x[0])) for x in nmz_triangulation]


#########################################################################
class Polyhedron_QQ_normaliz(Polyhedron_normaliz, Polyhedron_QQ):
    r"""
    Polyhedra over `\QQ` with normaliz.

    INPUT:

    - ``Vrep`` -- list ``[vertices, rays, lines]`` or ``None``
    - ``Hrep`` -- list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0), (1,0), (0,1)],
        ....:                rays=[(1,1)], lines=[],
        ....:                backend='normaliz', base_ring=QQ)
        sage: TestSuite(p).run()
    """

    @cached_method(do_pickle=True)
    def ehrhart_series(self, variable='t'):
        r"""
        Return the Ehrhart series of a compact rational polyhedron.

        The Ehrhart series is the generating function where the coefficient of
        `t^k` is number of integer lattice points inside the `k`-th dilation of
        the polytope.

        INPUT:

        - ``variable`` -- string (default: ``'t'``)

        OUTPUT: a rational function

        EXAMPLES::

            sage: S = Polyhedron(vertices=[[0,1], [1,0]], backend='normaliz')
            sage: ES = S.ehrhart_series()
            sage: ES.numerator()
            1
            sage: ES.denominator().factor()
            (t - 1)^2

            sage: C = Polyhedron(vertices=[[0,0,0], [0,0,1], [0,1,0], [0,1,1],
            ....:                          [1,0,0], [1,0,1], [1,1,0], [1,1,1]],
            ....:                backend='normaliz')
            sage: ES = C.ehrhart_series()
            sage: ES.numerator()
            t^2 + 4*t + 1
            sage: ES.denominator().factor()
            (t - 1)^4

        The following example is from the Normaliz manual contained in the file
        ``rational.in``::

            sage: rat_poly = Polyhedron(vertices=[[1/2,1/2], [-1/3,-1/3], [1/4,-1/2]],
            ....:                       backend='normaliz')
            sage: ES = rat_poly.ehrhart_series()
            sage: ES.numerator()
            2*t^6 + 3*t^5 + 4*t^4 + 3*t^3 + t^2 + t + 1
            sage: ES.denominator().factor()
            (-1) * (t + 1)^2 * (t - 1)^3 * (t^2 + 1) * (t^2 + t + 1)

        The polyhedron should be compact::

            sage: C = Polyhedron(rays=[[1,2], [2,1]], backend='normaliz')
            sage: C.ehrhart_series()
            Traceback (most recent call last):
            ...
            NotImplementedError: Ehrhart series can only be computed for compact polyhedron

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.backend_normaliz.hilbert_series`

        TESTS:

        Check that the Ehrhart series is pickled::

            sage: new_poly = loads(dumps(rat_poly))
            sage: new_poly.ehrhart_series.is_in_cache()
            True
        """
        if self.is_empty():
            return 0

        if not self.is_compact():
            raise NotImplementedError("Ehrhart series can only be computed for compact polyhedron")

        cone = self._normaliz_cone
        e = self._nmz_result(cone, "EhrhartSeries")
        # The output format of PyNormaliz is a list with 3 things:
        # 1) the coefficients of the h^*-polynomial
        # 2) a list of the exponents e such that (1-t^e) appears as a factor in
        # the denominator
        # 3) a shifting of the generating function.

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        poly_ring = PolynomialRing(ZZ, variable).fraction_field()
        t = poly_ring.gens()[0]
        es = sum([e[0][i] * t**i for i in range(len(e[0]))])
        for expo in range(len(e[1])):
            es = es / (1 - t**e[1][expo])

        # The shift:
        return es * t**e[2]

    def _ehrhart_quasipolynomial_normaliz(self, variable='t'):
        r"""
        Return the Ehrhart quasipolynomial of a compact rational polyhedron
        using Normaliz.

        If it is a polynomial, returns the polynomial. Otherwise, returns a
        tuple of rational polynomials whose length is the quasi-period of the
        quasipolynomial and each rational polynomial describes a residue class.

        INPUT:

        - ``variable`` -- string (default: ``'t'``)

        OUTPUT: a polynomial or tuple of polynomials

        EXAMPLES::

            sage: C = Polyhedron(vertices=[[0,0,0], [0,0,1], [0,1,0], [0,1,1],
            ....:                          [1,0,0], [1,0,1], [1,1,0], [1,1,1]],
            ....:                backend='normaliz')
            sage: C._ehrhart_quasipolynomial_normaliz()
            t^3 + 3*t^2 + 3*t + 1

            sage: P = Polyhedron(vertices=[[0,0], [3/2,0], [0,3/2], [1,1]], backend='normaliz')
            sage: P._ehrhart_quasipolynomial_normaliz()
            (3/2*t^2 + 2*t + 1, 3/2*t^2 + 2*t + 1/2)
            sage: P._ehrhart_quasipolynomial_normaliz('x')
            (3/2*x^2 + 2*x + 1, 3/2*x^2 + 2*x + 1/2)

        The quasipolynomial evaluated at ``i`` counts the integral points
        in the ``i``-th dilate::

            sage: Q = Polyhedron(vertices=[[-1/3], [2/3]], backend='normaliz')
            sage: p0,p1,p2 = Q._ehrhart_quasipolynomial_normaliz()
            sage: r0 = [p0(i) for i in range(15)]
            sage: r1 = [p1(i) for i in range(15)]
            sage: r2 = [p2(i) for i in range(15)]
            sage: result = [None]*15
            sage: result[::3] = r0[::3]
            sage: result[1::3] = r1[1::3]
            sage: result[2::3] = r2[2::3]
            sage: result == [(i*Q).integral_points_count() for i in range(15)]
            True


        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.backend_normaliz.hilbert_series`,
            :meth:`~sage.geometry.polyhedron.backend_normaliz.ehrhart_series`
        """
        cone = self._normaliz_cone
        # Normaliz needs to compute the EhrhartSeries first
        assert NmzCompute(cone, ["EhrhartSeries"])
        e = self._nmz_result(cone, "EhrhartQuasiPolynomial")

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        poly_ring = PolynomialRing(QQ, variable)
        t = poly_ring.gens()[0]
        if len(e) == 2:
            # It is a polynomial
            es = sum([e[0][i] * t**i for i in range(len(e[0]))])
            return es / ZZ(e[1])

        # It is a quasipolynomial
        polynomials = []
        for p in e[:-1]:
            es = sum([p[i] * t**i for i in range(len(p))]) / ZZ(e[-1])
            polynomials += [es]

        return tuple(polynomials)

    _ehrhart_polynomial_normaliz = _ehrhart_quasipolynomial_normaliz

    @cached_method(do_pickle=True, key=lambda self, g, v: (tuple(g), v))
    def hilbert_series(self, grading, variable='t'):
        r"""
        Return the Hilbert series of the polyhedron with respect to ``grading``.

        INPUT:

        - ``grading`` -- vector. The grading to use to form the Hilbert series

        - ``variable`` -- string (default: ``'t'``)

        OUTPUT: a rational function

        EXAMPLES::

            sage: C = Polyhedron(backend='normaliz',
            ....:                rays=[[0,0,1], [0,1,1], [1,0,1], [1,1,1]])
            sage: HS = C.hilbert_series([1,1,1])
            sage: HS.numerator()
            t^2 + 1
            sage: HS.denominator().factor()
            (-1) * (t + 1) * (t - 1)^3 * (t^2 + t + 1)

        By changing the grading, you can get the Ehrhart series of the square
        lifted at height 1::

            sage: C.hilbert_series([0,0,1])
            (t + 1)/(-t^3 + 3*t^2 - 3*t + 1)

        Here is an example ``2cone.in`` from the Normaliz manual::

            sage: C = Polyhedron(backend='normaliz', rays=[[1,3], [2,1]])
            sage: HS = C.hilbert_series([1,1])
            sage: HS.numerator()
            t^5 + t^4 + t^3 + t^2 + 1
            sage: HS.denominator().factor()
            (t + 1) * (t - 1)^2 * (t^2 + 1) * (t^2 + t + 1)

            sage: HS = C.hilbert_series([1,2])
            sage: HS.numerator()
            t^8 + t^6 + t^5 + t^3 + 1
            sage: HS.denominator().factor()
            (t + 1) * (t - 1)^2 * (t^2 + 1) * (t^6 + t^5 + t^4 + t^3 + t^2 + t + 1)

        Here is the magic square example form the Normaliz manual::

            sage: eq = [[0,1,1,1,-1,-1,-1, 0, 0, 0],
            ....:       [0,1,1,1, 0, 0, 0,-1,-1,-1],
            ....:       [0,0,1,1,-1, 0, 0,-1, 0, 0],
            ....:       [0,1,0,1, 0,-1, 0, 0,-1, 0],
            ....:       [0,1,1,0, 0, 0,-1, 0, 0,-1],
            ....:       [0,0,1,1, 0,-1, 0, 0, 0,-1],
            ....:       [0,1,1,0, 0,-1, 0,-1, 0, 0]]
            sage: magic_square = (Polyhedron(eqns=eq, backend='normaliz')
            ....:                 & Polyhedron(rays=identity_matrix(9).rows()))
            sage: grading = [1,1,1,0,0,0,0,0,0]
            sage: magic_square.hilbert_series(grading)
            (t^6 + 2*t^3 + 1)/(-t^9 + 3*t^6 - 3*t^3 + 1)

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.backend_normaliz.ehrhart_series`

        TESTS:

        Check that the Hilbert series is pickled::

            sage: new_magic = loads(dumps(magic_square))
            sage: new_magic.hilbert_series.is_in_cache(grading)
            True
        """
        if self.is_empty():
            return 0

        data = self._get_nmzcone_data()
        data['grading'] = [grading]
        new_cone = self._make_normaliz_cone(data)
        h = self._nmz_result(new_cone, "HilbertSeries")

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        poly_ring = PolynomialRing(ZZ, variable).fraction_field()
        t = poly_ring.gens()[0]
        hs = sum([h[0][i] * t**i for i in range(len(h[0]))])
        for expo in range(len(h[1])):
            hs = hs / (1 - t**h[1][expo])

        # The shift:
        return hs * t**h[2]

    def integral_points(self, threshold=10000) -> tuple:
        r"""
        Return the integral points in the polyhedron.

        Uses either the naive algorithm (iterate over a rectangular
        bounding box) or triangulation + Smith form.

        INPUT:

        - ``threshold`` -- integer (default: 10000); use the naïve
          algorithm as long as the bounding box is smaller than this

        OUTPUT:

        The list of integral points in the polyhedron. If the
        polyhedron is not compact, a :exc:`ValueError` is raised.

        EXAMPLES::

            sage: Polyhedron(vertices=[(-1,-1), (1,0), (1,1), (0,1)],
            ....:            backend='normaliz').integral_points()
            ((-1, -1), (0, 0), (0, 1), (1, 0), (1, 1))

            sage: simplex = Polyhedron([(1,2,3), (2,3,7), (-2,-3,-11)],
            ....:                      backend='normaliz')
            sage: simplex.integral_points()
            ((-2, -3, -11), (0, 0, -2), (1, 2, 3), (2, 3, 7))

        The polyhedron need not be full-dimensional::

            sage: simplex = Polyhedron([(1,2,3,5), (2,3,7,5), (-2,-3,-11,5)],
            ....:                      backend='normaliz')
            sage: simplex.integral_points()
            ((-2, -3, -11, 5), (0, 0, -2, 5), (1, 2, 3, 5), (2, 3, 7, 5))

            sage: point = Polyhedron([(2,3,7)],
            ....:                    backend='normaliz')
            sage: point.integral_points()
            ((2, 3, 7),)

            sage: empty = Polyhedron(backend='normaliz')
            sage: empty.integral_points()
            ()

        Here is a simplex where the naive algorithm of running over
        all points in a rectangular bounding box no longer works fast
        enough::

            sage: v = [(1,0,7,-1), (-2,-2,4,-3), (-1,-1,-1,4), (2,9,0,-5), (-2,-1,5,1)]
            sage: simplex = Polyhedron(v, backend='normaliz'); simplex
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 5 vertices
            sage: len(simplex.integral_points())
            49

        A rather thin polytope for which the bounding box method would
        be a very bad idea (note this is a rational (non-lattice)
        polytope, so the other backends use the bounding box method)::

            sage: P = Polyhedron(vertices=((0, 0), (178933,37121))) + 1/1000*polytopes.hypercube(2)
            sage: P = Polyhedron(vertices=P.vertices_list(),
            ....:                backend='normaliz')
            sage: len(P.integral_points())
            434

        Finally, the 3-d reflexive polytope number 4078::

            sage: v = [(1,0,0), (0,1,0), (0,0,1), (0,0,-1), (0,-2,1),
            ....:      (-1,2,-1), (-1,2,-2), (-1,1,-2), (-1,-1,2), (-1,-3,2)]
            sage: P = Polyhedron(v, backend='normaliz')
            sage: pts1 = P.integral_points()
            sage: all(P.contains(p) for p in pts1)
            True
            sage: pts2 = LatticePolytope(v).points()                                    # needs palp
            sage: for p in pts1: p.set_immutable()
            sage: set(pts1) == set(pts2)                                                # needs palp
            True

            sage: timeit('Polyhedron(v, backend='normaliz').integral_points()')  # not tested - random
            625 loops, best of 3: 1.41 ms per loop
            sage: timeit('LatticePolytope(v).points()')                          # not tested - random
            25 loops, best of 3: 17.2 ms per loop

        TESTS:

        Test some trivial cases (see :issue:`17937`):

        Empty polyhedron in 1 dimension::

            sage: P = Polyhedron(ambient_dim=1, backend='normaliz')
            sage: P.integral_points()
            ()

        Empty polyhedron in 0 dimensions::

            sage: P = Polyhedron(ambient_dim=0, backend='normaliz')
            sage: P.integral_points()
            ()

        Single point in 1 dimension::

            sage: P = Polyhedron([[3]], backend='normaliz')
            sage: P.integral_points()
            ((3),)

        Single non-integral point in 1 dimension::

            sage: P = Polyhedron([[1/2]], backend='normaliz')
            sage: P.integral_points()
            ()

        Single point in 0 dimensions::

            sage: P = Polyhedron([[]], backend='normaliz')
            sage: P.integral_points()
            ((),)

        A polytope with no integral points (:issue:`22938`)::

            sage: ieqs = [[1, 2, -1, 0], [0, -1, 2, -1], [0, 0, -1, 2],
            ....:         [0, -1, 0, 0], [0, 0, -1, 0],  [0, 0, 0, -1],
            ....:         [-1, -1, -1, -1], [1, 1, 0, 0], [1, 0, 1, 0],
            ....:         [1, 0, 0, 1]]
            sage: P = Polyhedron(ieqs=ieqs, backend='normaliz')
            sage: P.bounding_box()
            ((-3/4, -1/2, -1/4), (-1/2, -1/4, 0))
            sage: P.bounding_box(integral_hull=True)
            (None, None)
            sage: P.integral_points()
            ()

        Check the polytopes from :issue:`22984`::

            sage: base = [[0, 2, 0, -1, 0, 0, 0, 0, 0],
            ....:         [0, 0, 2, 0, -1, 0, 0, 0, 0],
            ....:         [1, -1, 0, 2, -1, 0, 0, 0, 0],
            ....:         [0, 0, -1, -1, 2, -1, 0, 0, 0],
            ....:         [0, 0, 0, 0, -1, 2, -1, 0, 0],
            ....:         [0, 0, 0, 0, 0, -1, 2, -1, 0],
            ....:         [1, 0, 0, 0, 0, 0, -1, 2, -1],
            ....:         [0, 0, 0, 0, 0, 0, 0, -1, 2],
            ....:         [0, -1, 0, 0, 0, 0, 0, 0, 0],
            ....:         [0, 0, -1, 0, 0, 0, 0, 0, 0],
            ....:         [0, 0, 0, -1, 0, 0, 0, 0, 0],
            ....:         [0, 0, 0, 0, -1, 0, 0, 0, 0],
            ....:         [0, 0, 0, 0, 0, -1, 0, 0, 0],
            ....:         [0, 0, 0, 0, 0, 0, -1, 0, 0],
            ....:         [0, 0, 0, 0, 0, 0, 0, -1, 0],
            ....:         [0, 0, 0, 0, 0, 0, 0, 0, -1],
            ....:         [-1, -1, -1, -1, -1, -1, -1, -1, -1]]

            sage: ieqs = base + [
            ....:         [2, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:         [4, 0, 1, 0, 0, 0, 0, 0, 0],
            ....:         [4, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:         [7, 0, 0, 0, 1, 0, 0, 0, 0],
            ....:         [6, 0, 0, 0, 0, 1, 0, 0, 0],
            ....:         [4, 0, 0, 0, 0, 0, 1, 0, 0],
            ....:         [2, 0, 0, 0, 0, 0, 0, 1, 0],
            ....:         [1, 0, 0, 0, 0, 0, 0, 0, 1]]
            sage: P = Polyhedron(ieqs=ieqs, backend='normaliz')
            sage: P.integral_points()
            ((-2, -2, -4, -5, -4, -3, -2, -1),
             (-2, -2, -4, -5, -4, -3, -2, 0),
             (-1, -2, -3, -4, -3, -2, -2, -1),
             (-1, -2, -3, -4, -3, -2, -1, 0),
             (-1, -1, -2, -2, -2, -2, -2, -1),
             (-1, -1, -2, -2, -1, -1, -1, 0),
             (-1, -1, -2, -2, -1, 0, 0, 0),
             (-1, 0, -2, -2, -2, -2, -2, -1),
             (0, -1, -1, -2, -2, -2, -2, -1),
             (0, 0, -1, -1, -1, -1, -1, 0))

            sage: ieqs = base + [
            ....:         [3, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:         [4, 0, 1, 0, 0, 0, 0, 0, 0],
            ....:         [6, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:         [8, 0, 0, 0, 1, 0, 0, 0, 0],
            ....:         [6, 0, 0, 0, 0, 1, 0, 0, 0],
            ....:         [4, 0, 0, 0, 0, 0, 1, 0, 0],
            ....:         [2, 0, 0, 0, 0, 0, 0, 1, 0],
            ....:         [1, 0, 0, 0, 0, 0, 0, 0, 1]]
            sage: P = Polyhedron(ieqs=ieqs, backend='normaliz')
            sage: P.integral_points()
            ((-3, -4, -6, -8, -6, -4, -2, -1),
             (-3, -4, -6, -8, -6, -4, -2, 0),
             (-2, -2, -4, -5, -4, -3, -2, -1),
             (-2, -2, -4, -5, -4, -3, -2, 0),
             (-1, -2, -3, -4, -3, -2, -2, -1),
             (-1, -2, -3, -4, -3, -2, -1, 0),
             (-1, -1, -2, -2, -2, -2, -2, -1),
             (-1, -1, -2, -2, -1, -1, -1, 0),
             (-1, -1, -2, -2, -1, 0, 0, 0),
             (-1, 0, -2, -2, -2, -2, -2, -1),
             (0, -1, -1, -2, -2, -2, -2, -1),
             (0, 0, -1, -1, -1, -1, -1, 0))
        """
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
        if threshold > 1:
            box_min, box_max = self.bounding_box(integral_hull=True)
            if box_min is None:
                return ()
            box_points = prod(max_coord - min_coord + 1
                              for min_coord, max_coord in zip(box_min, box_max))
            if box_points < threshold:
                from sage.geometry.integral_points import rectangular_box_points
                return rectangular_box_points(list(box_min), list(box_max), self)
        # Compute with normaliz
        points = []
        cone = self._normaliz_cone
        assert cone
        for g in self._nmz_result(cone, "ModuleGenerators"):
            assert g[-1] == 1
            points.append(vector(ZZ, g[:-1]))
        return tuple(points)

    def integral_points_generators(self):
        r"""
        Return the integral points generators of the polyhedron.

        Every integral point in the polyhedron can be written as a (unique)
        nonnegative linear combination of integral points contained in the three
        defining parts of the polyhedron: the integral points (the compact
        part), the recession cone, and the lineality space.

        OUTPUT:

        A tuple consisting of the integral points, the Hilbert basis of the
        recession cone, and an integral basis for the lineality space.

        EXAMPLES:

        Normaliz gives a nonnegative integer basis of the lineality space::

            sage: P = Polyhedron(backend='normaliz', lines=[[2,2]])
            sage: P.integral_points_generators()
            (((0, 0),), (), ((1, 1),))

        A recession cone generated by two rays::

            sage: C = Polyhedron(backend='normaliz', rays=[[1,2], [2,1]])
            sage: C.integral_points_generators()
            (((0, 0),), ((1, 1), (1, 2), (2, 1)), ())

        Empty polyhedron::

            sage: P = Polyhedron(backend='normaliz')
            sage: P.integral_points_generators()
            ((), (), ())
        """
        # Trivial cases: polyhedron with 0 vertices
        if self.n_vertices() == 0:
            return ((), (), ())
        # Compute with normaliz
        cone = self._normaliz_cone
        compact_part = []
        recession_cone_part = []
        lineality_part = []
        assert cone
        for g in self._nmz_result(cone, "ModuleGenerators"):
            assert g[-1] == 1
            compact_part.append(vector(ZZ, g[:-1]))

        for g in self._nmz_result(cone, "HilbertBasis"):
            assert g[-1] == 0
            recession_cone_part.append(vector(ZZ, g[:-1]))

        for g in self._nmz_result(cone, "MaximalSubspace"):
            assert g[-1] == 0
            lineality_part.append(vector(ZZ, g[:-1]))

        return tuple(compact_part), tuple(recession_cone_part), tuple(lineality_part)

    def _Hstar_function_normaliz(self, acting_group=None, output=None):
        r"""
        Return `H^*` as a rational function in `t` with coefficients in
        the ring of class functions of the ``acting_group`` of the polytope.

        As in [Stap2011]_, when ``self`` is the polytope `P`,
        `H^*(t) = (\sum_{m \geq 0} \chi_{mP} t^m)(\det(I-\rho(t)))`.
        The irreducible characters of ``acting_group`` form an orthonormal basis
        for the ring of class functions with values in `\CC`.
        The coefficients of `H^*(t)` are expressed in this basis.

        INPUT:

        - ``acting_group`` -- (default=None) a permgroup object. A subgroup of
          ``self``'s ``restricted_automorphism_group`` output as a permutation.
          If ``None``, it is set to the full ``restricted_automorphism_group``
          of ``self``. The acting group should always use output='permutation'.

        - ``output`` -- string. an output option. The allowed values are:

            * ``None`` (default): returns the rational function `H^*(t)`. `H^*` is
              a rational function in `t` with coefficients in the ring of
              class functions.
            * ``'e_series_list'``: string. Returns a list of the ehrhart series
              for the fixed subpolytopes of each conjugacy class representative.
            * ``'determinant_vec'``: string. Returns a list of the determinants
              of `Id-\rho*t` for each conjugacy class representative.
            * ``'Hstar_as_lin_comb'``: string. Returns a vector of the coefficients
              of the irreducible representations in the expression of `H^*`.
            * ``'prod_det_es'``: string. Returns a vector of the product of
              determinants and the Ehrhart series.
            * ``'complete'``: string. Returns a dictionary with Hstar,
              Hstar_as_lin_comb, the conjugacy class representatives,
              the character table of the acting group, and
              whether Hstar is effective.

        OUTPUT:

        The default output is the rational function `H^*`. `H^*` is a rational
        function in `t` with coefficients in the ring of class functions.
        There are several output options to see the intermediary outputs of the
        function.

        EXAMPLES:

        The `H^*`-polynomial of the standard `d-1` dimensional simplex
        `S = conv(e_1, \dots, e_d)` under its ``restricted_automorphism_group``
        is equal to 1 = `\chi_{trivial}` (Prop 6.1 [Stap2011]_).
        Here is the computation for the 3-dimensional standard simplex::

            sage: # needs sage.groups
            sage: S = polytopes.simplex(3, backend='normaliz'); S
            A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 4 vertices
            sage: G = S.restricted_automorphism_group(output='permutation')
            sage: G.is_isomorphic(SymmetricGroup(4))
            True
            sage: len(G)
            24
            sage: Hstar = S._Hstar_function_normaliz(G); Hstar
            chi_4
            sage: G.character_table()
            [ 1 -1  1  1 -1]
            [ 3 -1  0 -1  1]
            [ 2  0 -1  2  0]
            [ 3  1  0 -1 -1]
            [ 1  1  1  1  1]

        The next example is Example 7.6 in [Stap2011]_, and shows that `H^*`
        is not always a polynomial. Let P be the polytope with vertices
        `\pm(0,0,1),\pm(1,0,1), \pm(0,1,1), \pm(1,1,1)` and let
        G = `\Zmod{2}` act on P as follows::

            sage: # needs sage.groups
            sage: P = Polyhedron(vertices=[[0,0,1], [0,0,-1], [1,0,1], [-1,0,-1],
            ....:                          [0,1,1], [0,-1,-1], [1,1,1], [-1,-1,-1]],
            ....:                backend='normaliz')
            sage: K = P.restricted_automorphism_group(output='permutation')
            sage: G = K.subgroup(gens=[K([(0,2),(1,3),(4,6),(5,7)])])
            sage: conj_reps = G.conjugacy_classes_representatives()
            sage: Dict = P.permutations_to_matrices(conj_reps, acting_group=G)
            sage: list(Dict.keys())[0]
            (0,2)(1,3)(4,6)(5,7)
            sage: list(Dict.values())[0]
            [-1  0  1  0]
            [ 0  1  0  0]
            [ 0  0  1  0]
            [ 0  0  0  1]
            sage: len(G)
            2
            sage: G.character_table()
            [ 1  1]
            [ 1 -1]

        Then we calculate the rational function `H^*(t)`::

            sage: Hst = P._Hstar_function_normaliz(G); Hst                              # needs sage.groups
            (chi_0*t^4 + (3*chi_0 + 3*chi_1)*t^3
              + (8*chi_0 + 2*chi_1)*t^2 + (3*chi_0 + 3*chi_1)*t + chi_0)/(t + 1)

        To see the exact as written in [Stap2011]_, we can format it as
        ``'Hstar_as_lin_comb'``. The first coordinate is the coefficient of the
        trivial character; the second is the coefficient of the sign character::

            sage: lin = P._Hstar_function_normaliz(G, output='Hstar_as_lin_comb'); lin  # needs sage.groups
            ((t^4 + 3*t^3 + 8*t^2 + 3*t + 1)/(t + 1), (3*t^3 + 2*t^2 + 3*t)/(t + 1))
        """
        from sage.groups.conjugacy_classes import ConjugacyClassGAP
        from sage.rings.number_field.number_field import CyclotomicField
        from sage.rings.qqbar import QQbar
        from sage.matrix.matrix_space import MatrixSpace
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.matrix.special import identity_matrix
        # Setting the group
        G_perm = self.restricted_automorphism_group(output='permutation')

        if acting_group is not None:
            if not acting_group.is_subgroup(G_perm):
                raise TypeError("the 'acting_group' should be a subgroup of the 'restricted_automorphism_group'.")
            G_perm = acting_group
        # Create the Gap group one time only (each creation has different conj reps)
        G_perm_gap = G_perm._libgap_()

        # Fixing the conjugacy classes representatives once and for all
        cls = G_perm_gap.ConjugacyClasses()
        L = [cl.Representative() for cl in cls]
        conj_classes = [ConjugacyClassGAP(G_perm, G_perm.element_class(rep, G_perm, check=False)) for rep in L]
        conj_reps = [cl[0] for cl in conj_classes]

        # Creating the Character Table
        n_classes = len(conj_reps)
        irrG_perm_gap = G_perm_gap.Irr()
        ct = [[irrG_perm_gap[i, j] for j in range(n_classes)] for i in range(n_classes)]
        e = irrG_perm_gap.Flat().Conductor()
        K = CyclotomicField(e)
        ct = [[K(x) for x in v] for v in ct]
        # Finally return the result as a matrix.
        MS = MatrixSpace(K, n_classes)
        char_initial = MS(ct)

        # A check on whether the character table has permuted columns
        tbl = G_perm_gap.CharacterTable()
        perm = tbl.IdentificationOfConjugacyClasses()
        ident_perm = list(range(1, 1 + n_classes))
        assert perm == ident_perm, "The conjugacy classes don't match with the character table"

        # Create fixed subpolytopes and their Ehrhart series
        group_dict = self.permutations_to_matrices(conj_reps, acting_group)
        fix_polys = self.fixed_subpolytopes(conj_reps)
        list_es = [fix_polys[g].ehrhart_series() for g in conj_reps]

        # get the list of the denominators det([Id - rho (t)])
        ring = PolynomialRing(QQbar, 't')
        det_vector = list()
        dim = group_dict[G_perm.gens()[0]].dimensions()[0]
        t = ring.gens()[0]
        ts_matrix = t * identity_matrix(ring, dim)
        identity = identity_matrix(ring, dim)

        # ix the determinant if polytope isn't full dimensional
        codim = self.ambient_dim() - self.dim()
        for perm in conj_reps:
            mat = group_dict[perm]
            mat = mat.change_ring(ring)
            new_matrix = identity - mat * ts_matrix
            det = (1 - t)**-codim * (new_matrix.determinant())
            det_vector.append(det)

        FF = ring.fraction_field()
        initial_result = vector(FF, [a * b for a, b in zip(det_vector, list_es)])
        Char = char_initial.change_ring(FF)
        new_result = Char.solve_left(initial_result)

        new_new_result = self._Hstar_as_rat_fct(new_result)
        if output is None:
            return new_new_result
        if output == 'e_series_list':
            return list_es
        if output == 'determinant_vec':
            return det_vector
        if output == 'Hstar_as_lin_comb':
            return new_result
        if output == 'prod_det_es':
            return initial_result
        if output == 'complete':
            results_dictionary = {}
            results_dictionary['Hstar'] = new_new_result
            results_dictionary['Hstar_as_lin_comb'] = new_result
            results_dictionary['conjugacy_class_reps'] = conj_reps
            results_dictionary['character_table'] = char_initial
            results_dictionary['is_effective'] = self._is_effective_normaliz(new_new_result, new_result)
            return results_dictionary

    def _Hstar_as_rat_fct(self, initial_Hstar):
        r"""
        Rewrite the vector representing `H^*(t)` given as a linear combination
        of the irreducible representations of the acting group as a rational
        function in `t`.

        INPUT:

        - ``initial_Hstar`` -- a vector of rational functions in `t`

        OUTPUT:

        A rational function in `t` with coefficients in the ring of class functions
        of ``self.restricted_automorphism_group()``.

        EXAMPLES:

        The expression of `H^*` as a polynomial in `t` for a 3-dimensional simplex
        is computed as follows::

            sage: simplex = Polyhedron(vertices=[[0,0,0], [1,0,0],
            ....:                                [0,1,0], [0,0,1]], backend='normaliz')
            sage: Hstar = simplex.Hstar_function(); Hstar  # indirect doctest           # needs sage.rings.number_field
            chi_4

        The polynomial is `\chi_4 \cdot t^0`. We can see which irreducible
        representation `\chi_4` corresponds to by looking at the character table::

            sage: G = simplex.restricted_automorphism_group(output='permutation')       # needs sage.groups
            sage: char = G.character_table(); char                                      # needs sage.groups
            [ 1 -1  1  1 -1]
            [ 3 -1  0 -1  1]
            [ 2  0 -1  2  0]
            [ 3  1  0 -1 -1]
            [ 1  1  1  1  1]

        Thus `\chi_4` corresponds to the trivial representation of the group, and
        for every element in the group, it evaluates to 1.

        As another example, we can look at `H^*(t)` for the `\pm 1` square::

            sage: square = Polyhedron(vertices=[[1,1], [-1,1], [-1,-1], [1,-1]],
            ....:                     backend='normaliz')
            sage: Hstar = square.Hstar_function(); Hstar                                # needs sage.rings.number_field
            chi_0*t^2 + (2*chi_0 + chi_2 + chi_3 + chi_4)*t + chi_0

        Plugging in the values from the first column of the character table below
        yields the `h^*`-polynomial of the square, `t^2+6t+1`::

            sage: G = square.restricted_automorphism_group(output='permutation')        # needs sage.groups
            sage: G.character_table()                                                   # needs sage.groups
            [ 1  1  1  1  1]
            [ 1 -1 -1  1  1]
            [ 1 -1  1 -1  1]
            [ 1  1 -1 -1  1]
            [ 2  0  0  0 -2]
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.qqbar import QQbar
        chi_vars = ','.join(f'chi_{i}' for i in range(len(initial_Hstar)))
        Chi_ring = PolynomialRing(QQbar, chi_vars)
        virtual_ring = PolynomialRing(Chi_ring, initial_Hstar.base_ring().gens())
        fraction_virtual_ring = virtual_ring.fraction_field()
        return initial_Hstar.change_ring(fraction_virtual_ring) * vector(fraction_virtual_ring, Chi_ring.gens())

    def _is_effective_normaliz(self, Hstar, Hstar_as_lin_comb):
        r"""
        Test for the effectiveness of the ``Hstar`` series of this polytope.

        The ``Hstar`` series of the polytope is determined by the action of a
        subgroup of the polytope's ``restricted_automorphism_group``. The
        ``Hstar`` series is effective if it is a polynomial in `t` and the
        coefficient of each `t^i` is an effective character in the ring of
        class functions of the acting group. A character `\rho` is effective if
        the coefficients of the irreducible representations in the expression
        of `\rho` are nonnegative integers.

        INPUT:

        - ``Hstar`` -- a rational function in `t` with coefficients in the ring
          of class functions

        - ``Hstar_as_lin_comb`` -- vector. The coefficients of the irreducible
          representations of the acting group in the expression of ``Hstar`` as
          a linear combination of irreducible representations with coefficients
          in the field of rational functions in `t`.

        OUTPUT: boolean; whether the ``Hstar`` series is effective

        EXAMPLES:

        The `H^*` series of the two-dimensional permutahedron under the action
        of the symmetric group is effective::

            sage: # needs sage.groups
            sage: p3 = polytopes.permutahedron(3, backend='normaliz')
            sage: G = p3.restricted_automorphism_group(output='permutation')
            sage: reflection12 = G([(0,2),(1,4),(3,5)])
            sage: reflection23 = G([(0,1),(2,3),(4,5)])
            sage: S3 = G.subgroup(gens=[reflection12, reflection23])
            sage: S3.is_isomorphic(SymmetricGroup(3))
            True
            sage: Hstar = p3.Hstar_function(S3)                                         # needs sage.rings.number_field
            sage: Hlin  = p3.Hstar_function(S3, output='Hstar_as_lin_comb')             # needs sage.rings.number_field
            sage: p3._is_effective_normaliz(Hstar, Hlin)                                # needs sage.rings.number_field
            True

        If the `H^*`-series is not polynomial, then it is not effective::

            sage: # needs sage.groups
            sage: P = Polyhedron(vertices=[[0,0,1], [0,0,-1], [1,0,1], [-1,0,-1],
            ....:                          [0,1,1], [0,-1,-1], [1,1,1], [-1,-1,-1]],
            ....:                backend='normaliz')
            sage: G = P.restricted_automorphism_group(output='permutation')
            sage: H = G.subgroup(gens = [G([(0,2),(1,3),(4,6),(5,7)])])
            sage: Hstar = P.Hstar_function(H); Hstar                                    # needs sage.rings.number_field
            (chi_0*t^4 + (3*chi_0 + 3*chi_1)*t^3
              + (8*chi_0 + 2*chi_1)*t^2 + (3*chi_0 + 3*chi_1)*t + chi_0)/(t + 1)
            sage: Hstar_lin = P.Hstar_function(H, output='Hstar_as_lin_comb')           # needs sage.rings.number_field
            sage: P._is_effective_normaliz(Hstar, Hstar_lin)                            # needs sage.rings.number_field
            False
        """
        if not Hstar.denominator().is_unit():
            return False
        for irrep in range(len(Hstar_as_lin_comb)):
            coeffs = Hstar_as_lin_comb[irrep].numerator().coefficients()
            for i in coeffs:
                if not i.is_integer() or i < 0:
                    return False
        return True


#########################################################################
class Polyhedron_ZZ_normaliz(Polyhedron_QQ_normaliz, Polyhedron_ZZ):
    r"""
    Polyhedra over `\ZZ` with normaliz.

    INPUT:

    - ``Vrep`` -- list ``[vertices, rays, lines]`` or ``None``
    - ``Hrep`` -- list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0), (1,0), (0,1)],
        ....:                rays=[(1,1)], lines=[],
        ....:                backend='normaliz', base_ring=ZZ)
        sage: TestSuite(p).run()
    """
    pass
