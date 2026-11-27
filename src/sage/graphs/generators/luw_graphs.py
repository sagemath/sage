r"""
Graphs defined by systems of equations

This module implements a class of bipartite graphs defined by triangular
systems of equations, popularized by Lazebnik, Ustimenko, and Woldar. More
precisely, let `R` be a finite commutative ring. The graph has point
part `P = R^n` and line part `L = R^n`. For `2 \leq i \leq n`, let
`f_i` be a polynomial function in `2i - 2` variables. A point
`(p_1, p_2, \ldots, p_n)` is adjacent to a line
`(l_1, l_2, \ldots, l_n)` if

.. MATH::

    p_i + l_i = f_i(p_1, l_1, \ldots, p_{i-1}, l_{i-1}),
    \qquad 2 \leq i \leq n.

The class :class:`LUWGraphDescriptor` stores a validated set of equations
which define the graph. This lets users experiment with large examples
without paying the cost of building the full graph. Currently, the
descriptor-returning constructors are
:func:`define_luw_graph`, :func:`define_Akq`, :func:`define_Dkq`, and
:func:`define_WengerGraph`. The graph-returning constructors :func:`Akq`,
:func:`Dkq`, :func:`LUWGraph`, and :func:`WengerGraph` are exposed through
:mod:`sage.graphs.graph_generators`, and so are available as
``graphs.LUWGraph(...)``, ``graphs.WengerGraph(...)``, and so on.

The general construction and the graph families implemented here are surveyed
in [LW2026]_.

EXAMPLES:

The graph constructor is available directly from ``graphs``::

    sage: F = GF(3)
    sage: R = PolynomialRing(F, names=("p1", "l1", "p2", "l2"))
    sage: p1, l1, p2, l2 = R.gens()
    sage: G = graphs.LUWGraph(F, [p1*l1, p1*l2])
    sage: G.order(), G.size(), G.girth()
    (54, 81, 8)

Use :func:`define_luw_graph` when the algebraic description should be kept as
a :class:`LUWGraphDescriptor` without immediately constructing the full graph.

AUTHORS:

- Vladislav Taranchuk (2026-04): initial version
"""

# ****************************************************************************
#       Copyright (C) 2026 Vladislav Taranchuk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from itertools import product
from math import prod
from warnings import warn

from sage.graphs.graph import Graph
from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class LUWGraphDescriptor:
    r"""
    A validated algebraic description of a bipartite graph.

    The graph itself is built from the stored equations and coordinate
    sets. See [LW2026]_ for a survey of these graphs.

    EXAMPLES::

        sage: from sage.graphs.generators.luw_graphs import define_luw_graph
        sage: F = GF(5)
        sage: R = PolynomialRing(F, names=("p1", "l1"))
        sage: p1, l1 = R.gens()
        sage: D = define_luw_graph(F, [p1*l1])
        sage: D.order()
        50
        sage: D.size()
        125

    REFERENCES:

    - [LW2026]_
    """

    def __init__(self, ring, equations, name,
                 point_coordinate_sets, line_coordinate_sets):
        r"""
        Initialize ``self``.

        INPUT:

        - ``ring`` -- finite commutative ring

        - ``equations`` -- tuple of validated defining polynomials

        - ``name`` -- string; name of the algebraic descriptor

        - ``point_coordinate_sets`` -- tuple of allowed coordinate sets on the
          point side

        - ``line_coordinate_sets`` -- tuple of allowed coordinate sets on the
          line side

        TESTS::

            sage: import inspect
            sage: from sage.graphs.generators.luw_graphs import LUWGraphDescriptor
            sage: "equations" in str(inspect.signature(LUWGraphDescriptor))
            True
            sage: "_equations" in str(inspect.signature(LUWGraphDescriptor))
            False
        """
        self.ring = ring
        self._equations = tuple(equations)
        self.name = name
        self.point_coordinate_sets = point_coordinate_sets
        self.line_coordinate_sets = line_coordinate_sets

    def dimension(self):
        r"""
        Return the number of point or line coordinates.

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_luw_graph
            sage: F = GF(3)
            sage: R = PolynomialRing(F, names=("p1", "l1"))
            sage: p1, l1 = R.gens()
            sage: define_luw_graph(F, [p1*l1]).dimension()
            2
        """
        return len(self._equations) + 1

    def base_ring(self):
        r"""
        Return the underlying finite ring.

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_luw_graph
            sage: Z4 = Integers(4)
            sage: R = PolynomialRing(Z4, names=("p1", "l1"))
            sage: p1, l1 = R.gens()
            sage: define_luw_graph(Z4, [p1*l1]).base_ring()
            Ring of integers modulo 4
        """
        return self.ring

    def point_first_coordinates(self):
        r"""
        Return the allowed first coordinates on the point side.

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_luw_graph
            sage: F = GF(3)
            sage: R = PolynomialRing(F, names=("p1", "l1"))
            sage: p1, l1 = R.gens()
            sage: define_luw_graph(F, [p1*l1], A=[0, 1]).point_first_coordinates()
            (0, 1)
        """
        return self.point_coordinate_sets[0]

    def line_first_coordinates(self):
        r"""
        Return the allowed first coordinates on the line side.

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_luw_graph
            sage: F = GF(3)
            sage: R = PolynomialRing(F, names=("p1", "l1"))
            sage: p1, l1 = R.gens()
            sage: define_luw_graph(F, [p1*l1], B=[1, 2]).line_first_coordinates()
            (1, 2)
        """
        return self.line_coordinate_sets[0]

    def point_count(self):
        r"""
        Return the number of point vertices.

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_luw_graph
            sage: F = GF(3)
            sage: R = PolynomialRing(F, names=("p1", "l1"))
            sage: p1, l1 = R.gens()
            sage: define_luw_graph(F, [p1*l1], A=[0, 1]).point_count()
            6
        """
        return prod(len(values) for values in self.point_coordinate_sets)

    def line_count(self):
        r"""
        Return the number of line vertices.

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_luw_graph
            sage: F = GF(3)
            sage: R = PolynomialRing(F, names=("p1", "l1"))
            sage: p1, l1 = R.gens()
            sage: define_luw_graph(F, [p1*l1], B=[1, 2]).line_count()
            6
        """
        return prod(len(values) for values in self.line_coordinate_sets)

    def order(self):
        r"""
        Return the number of vertices of the graph.

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_luw_graph
            sage: F = GF(3)
            sage: R = PolynomialRing(F, names=("p1", "l1"))
            sage: p1, l1 = R.gens()
            sage: define_luw_graph(F, [p1*l1]).order()
            18
        """
        return self.point_count() + self.line_count()

    def size(self):
        r"""
        Return the number of edges of the graph.

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_luw_graph
            sage: F = GF(3)
            sage: R = PolynomialRing(F, names=("p1", "l1"))
            sage: p1, l1 = R.gens()
            sage: define_luw_graph(F, [p1*l1]).size()
            27
        """
        return self.point_count() * len(self.line_first_coordinates())

    def adjacency_matrix(self, vertices=None):
        r"""
        Return the adjacency matrix without first constructing the graph.

        INPUT:

        - ``vertices`` -- list, tuple, or ``None`` (default: ``None``); vertex
          order

        When ``vertices`` is ``None``, the vertex order agrees with Sage's
        default graph adjacency-matrix order, namely sorted vertex labels.

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_luw_graph
            sage: F = GF(3)
            sage: R = PolynomialRing(F, names=("p1", "l1", "p2", "l2"))
            sage: p1, l1, p2, l2 = R.gens()
            sage: D = define_luw_graph(F, [p1*l1, p1*l2 + p2])
            sage: D.adjacency_matrix().dimensions()
            (54, 54)
            sage: D.adjacency_matrix() == D.graph().adjacency_matrix()
            True

        A custom vertex order can also be supplied::

            sage: vertices = list(D.graph().vertices(sort=False))
            sage: D.adjacency_matrix(vertices=vertices) == \
            ....:     D.graph().adjacency_matrix(vertices=vertices)
            True
        """
        point_vertices = tuple(("P", point) for point in _point_tuples(self))
        line_vertices = tuple(("L", line) for line in _line_tuples(self))
        all_vertices = point_vertices + line_vertices

        if vertices is None:
            vertices = tuple(sorted(all_vertices))
        else:
            vertices = tuple(_normalize_vertex(self, vertex)
                             for vertex in vertices)
            if (len(vertices) != len(all_vertices)
                    or set(vertices) != set(all_vertices)):
                raise ValueError(
                    "vertices must contain each graph vertex exactly once"
                )

        # Record where each vertex appears in the chosen matrix order.
        # For example, if ``vertices`` begins ``[v0, v1, v2]``, then this
        # dictionary stores ``{v0: 0, v1: 1, v2: 2}``.
        vertex_indices = {vertex: index for index, vertex in enumerate(vertices)}
        adjacency = matrix(ZZ, len(vertices), len(vertices), sparse=True)

        for point in _point_tuples(self):
            point_vertex = ("P", point)
            point_index = vertex_indices[point_vertex]
            for first_line_value in self.line_first_coordinates():
                line_vertex = _complete_line_from_point(self, point, first_line_value)
                line_index = vertex_indices.get(line_vertex)
                if line_index is None:
                    continue
                adjacency[point_index, line_index] = 1
                adjacency[line_index, point_index] = 1
        return adjacency

    def equations(self):
        r"""
        Return the validated defining polynomials.

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_Akq
            sage: define_Akq(3, 3).equations()
            (p1*l1, p1*l2)
        """
        return self._equations

    def graph(self, immutable=False):
        r"""
        Build the full graph corresponding to ``self``.

        INPUT:

        - ``immutable`` -- boolean (default: ``False``); whether to return an
          immutable or a mutable graph

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_WengerGraph
            sage: G = define_WengerGraph(2, 3).graph()
            sage: G.order(), G.size()
            (54, 81)
            sage: define_WengerGraph(2, 3).graph(immutable=True).is_immutable()
            True
        """
        def edges():
            for point in _point_tuples(self):
                point_vertex = ("P", point)
                for first_line_value in self.line_first_coordinates():
                    yield (point_vertex,
                           _complete_line_from_point(self, point,
                                                     first_line_value))

        G = Graph(edges(), format="list_of_edges",
                  name=self.name, immutable=immutable)
        G._luw_graph_metadata = {
            "definition": self,
            "ring": self.ring,
            "dimension": self.dimension(),
            "equations": self._equations,
            "point_coordinate_sets": self.point_coordinate_sets,
            "line_coordinate_sets": self.line_coordinate_sets,
        }
        return G

    def ball(self, start_vertex, depth, immutable=False):
        r"""
        Build the ball of radius ``depth`` around ``start_vertex``.

        The construction stops early if a new layer adds no vertices.

        INPUT:

        - ``start_vertex`` -- pair ``(side, coordinates)`` where ``side`` is
          either ``"P"`` or ``"L"``

        - ``depth`` -- nonnegative integer; radius of the ball

        - ``immutable`` -- boolean (default: ``False``); whether to return an
          immutable or a mutable graph

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_luw_graph
            sage: F = GF(3)
            sage: R = PolynomialRing(F, names=("p1", "l1"))
            sage: p1, l1 = R.gens()
            sage: D = define_luw_graph(F, [p1*l1])
            sage: B = D.ball(("P", (0, 0)), 1)
            sage: B.order(), B.size()
            (4, 3)
            sage: D.ball(("P", (0, 0)), 1, immutable=True).is_immutable()
            True
            sage: D = define_luw_graph(F, [p1*l1], A=[0], B=[0])
            sage: B = D.ball(("P", (0, 0)), 5)
            sage: B.order(), B._luw_graph_ball_metadata["component_recovered"]
            (2, True)
        """
        if depth < 0:
            raise ValueError("the depth must be nonnegative")

        start_vertex = _normalize_vertex(self, start_vertex)

        edges = []
        distances = {start_vertex: 0}
        frontier = {start_vertex}
        reached_depth = 0
        component_recovered = False

        for current_depth in range(depth):
            next_frontier = set()
            for vertex in frontier:
                for neighbor in self.neighbors(vertex):
                    edges.append((vertex, neighbor))
                    if neighbor not in distances:
                        distances[neighbor] = current_depth + 1
                        next_frontier.add(neighbor)

            reached_depth = current_depth + 1
            if not next_frontier:
                component_recovered = True
                break
            frontier = next_frontier

        G = Graph([distances, edges], format="vertices_and_edges",
                  name=f"Ball of radius {depth} in {self.name}",
                  immutable=immutable)
        G._luw_graph_ball_metadata = {
            "definition": self,
            "ring": self.ring,
            "dimension": self.dimension(),
            "equations": self._equations,
            "point_coordinate_sets": self.point_coordinate_sets,
            "line_coordinate_sets": self.line_coordinate_sets,
            "center": start_vertex,
            "requested_depth": depth,
            "reached_depth": reached_depth,
            "component_recovered": component_recovered,
        }
        G._luw_graph_distances = distances
        return G

    def neighbors(self, vertex):
        r"""
        Return the neighbors of ``vertex`` without building the full graph.

        INPUT:

        - ``vertex`` -- pair ``(side, coordinates)`` where ``side`` is either
          ``"P"`` or ``"L"``

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_luw_graph
            sage: F = GF(3)
            sage: R = PolynomialRing(F, names=("p1", "l1"))
            sage: p1, l1 = R.gens()
            sage: D = define_luw_graph(F, [p1*l1])
            sage: v = ("P", (F(0), F(0)))
            sage: D.neighbors(v)
            (('L', (0, 0)), ('L', (1, 0)), ('L', (2, 0)))
            sage: G = D.graph()
            sage: set(D.neighbors(v)) == set(G.neighbors(v))
            True
        """
        side, coords = _normalize_vertex(self, vertex)

        if side == "P":
            return tuple(
                _complete_line_from_point(self, coords, first_line_value)
                for first_line_value in self.line_first_coordinates()
            )
        return tuple(
            _complete_point_from_line(self, coords, first_point_value)
            for first_point_value in self.point_first_coordinates()
        )

    def lift(self, new_functions, name=None):
        r"""
        Return the lifted algebraic definition obtained by appending equations.

        INPUT:

        - ``new_functions`` -- iterable of Sage polynomials to append

        - ``name`` -- string (default: ``None``); name of the new definition

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_luw_graph
            sage: F = GF(3)
            sage: R = PolynomialRing(F, names=("p1", "l1", "p2", "l2"))
            sage: p1, l1, p2, l2 = R.gens()
            sage: D = define_luw_graph(F, [p1*l1])
            sage: D2 = D.lift([p1*l2])
            sage: D2.equations()
            (p1*l1, p1*l2)
        """
        new_functions = tuple(new_functions)
        if not new_functions:
            raise ValueError("give at least one equation to append")

        return define_luw_graph(
            self.ring,
            self._equations + new_functions,
            name=name or f"Lift of {self.name}",
            point_coordinate_sets=self.point_coordinate_sets
            + (self.base_ring(),) * len(new_functions),
            line_coordinate_sets=self.line_coordinate_sets
            + (self.base_ring(),) * len(new_functions),
        )

    def project_down(self, x, name=None):
        r"""
        Return the projection down obtained by removing the last ``x`` equations.

        INPUT:

        - ``x`` -- positive integer; number of final equations to remove

        - ``name`` -- string (default: ``None``); name of the new definition

        EXAMPLES::

            sage: from sage.graphs.generators.luw_graphs import define_luw_graph
            sage: F = GF(3)
            sage: R = PolynomialRing(F, names=("p1", "l1", "p2", "l2"))
            sage: p1, l1, p2, l2 = R.gens()
            sage: D = define_luw_graph(F, [p1*l1, p1*l2])
            sage: D.project_down(1).equations()
            (p1*l1,)
        """
        if x not in ZZ:
            raise ValueError("x must be a positive integer")

        x = ZZ(x)
        if x <= 0:
            raise ValueError("x must be a positive integer")
        if x >= len(self._equations):
            raise ValueError("cannot remove all defining equations")

        return define_luw_graph(
            self.ring,
            self._equations[:-x],
            name=name or f"Projection of {self.name}",
            point_coordinate_sets=self.point_coordinate_sets[:-x],
            line_coordinate_sets=self.line_coordinate_sets[:-x],
        )


def _canonical_variable_names(dimension):
    r"""
    Return the canonical variable order ``p1, l1, p2, l2, ...``.

    EXAMPLES::

        sage: from sage.graphs.generators import luw_graphs as luw
        sage: luw._canonical_variable_names(3)
        ('p1', 'l1', 'p2', 'l2')
    """
    return tuple(f"{kind}{i}" for i in range(1, dimension) for kind in ("p", "l"))


def _canonical_polynomial_ring(base_ring, dimension):
    r"""
    Return the canonical polynomial ring for the given dimension.

    EXAMPLES::

        sage: from sage.graphs.generators import luw_graphs as luw
        sage: luw._canonical_polynomial_ring(GF(3), 3)
        Multivariate Polynomial Ring in p1, l1, p2, l2 over Finite Field of size 3
    """
    return PolynomialRing(base_ring, names=_canonical_variable_names(dimension))


def _canonicalize_defining_polynomial(base_ring, polynomial, canonical_ring):
    r"""
    Rewrite ``polynomial`` into ``canonical_ring`` by variable name.

    The input polynomial may come from a polynomial ring whose generators are
    ordered differently from the LUW order ``p1, l1, p2, l2, ...``. The
    polynomial is stored in the canonical ring by matching generator names, not
    by matching generator positions. If the used variables appear in a different
    order in the input parent, a warning is issued.

    EXAMPLES::

        sage: from sage.graphs.generators import luw_graphs as luw
        sage: F = GF(3)
        sage: R = PolynomialRing(F, names=("p2", "l1", "p1", "l2"))
        sage: p2, l1, p1, l2 = R.gens()
        sage: C = luw._canonical_polynomial_ring(F, 3)
        sage: luw._canonicalize_defining_polynomial(F, p1*l1, C)
        doctest:...: UserWarning: defining polynomial variables were reordered by name into LUW order p1, l1, p2, l2, ...
        p1*l1
    """
    if not hasattr(polynomial, "parent") or not hasattr(polynomial.parent(), "ngens"):
        raise TypeError("each defining function must be a Sage polynomial")

    parent_ring = polynomial.parent()
    if parent_ring.base_ring() != base_ring:
        raise ValueError("each polynomial must be defined over the input ring")

    target_indices = {
        name: index for index, name in enumerate(canonical_ring.variable_names())
    }
    source_names = parent_ring.variable_names()
    used_source_indices = set()
    terms = {}
    for exponents, coefficient in polynomial.dict().items():
        target_exponents = [0] * canonical_ring.ngens()
        for source_index, exponent in enumerate(exponents):
            if not exponent:
                continue

            name = source_names[source_index]
            if name not in target_indices:
                expected = ", ".join(canonical_ring.variable_names())
                raise ValueError(f"unknown variable {name}; expected variables named {expected}")
            used_source_indices.add(source_index)
            target_exponents[target_indices[name]] += exponent

        target_exponents = tuple(target_exponents)
        terms[target_exponents] = terms.get(target_exponents, base_ring(0))
        terms[target_exponents] += base_ring(coefficient)

    used_names = tuple(source_names[i] for i in sorted(used_source_indices))
    canonical_used_names = tuple(sorted(used_names, key=target_indices.__getitem__))
    if used_names != canonical_used_names:
        warn(
            "defining polynomial variables were reordered by name into "
            "LUW order p1, l1, p2, l2, ...",
            stacklevel=3,
        )

    return canonical_ring(terms)


def _validate_defining_polynomial(polynomial, available_values):
    r"""
    Check that ``polynomial`` uses only coordinates already available.

    TESTS::

        sage: from sage.graphs.generators import luw_graphs as luw
        sage: R = PolynomialRing(GF(3), names=("p1", "l1", "p2", "l2"))
        sage: p1, l1, p2, l2 = R.gens()
        sage: luw._validate_defining_polynomial(p1*l1, 2)
        sage: luw._validate_defining_polynomial(p2, 2)
        Traceback (most recent call last):
        ...
        ValueError: a defining polynomial uses a variable that is not available yet
    """
    ring = polynomial.parent()
    used_indices = [ring.gens().index(variable) for variable in polynomial.variables()]
    if used_indices and max(used_indices) >= available_values:
        raise ValueError("a defining polynomial uses a variable that is not available yet")


def _normalize_coordinate_subset(base_ring, values):
    r"""
    Normalize a user-given coordinate subset to a tuple of ring elements.

    EXAMPLES::

        sage: from sage.graphs.generators import luw_graphs as luw
        sage: luw._normalize_coordinate_subset(GF(3), [0, 1, 0])
        (0, 1)
        sage: luw._normalize_coordinate_subset(GF(3), None)
        Finite Field of size 3
    """
    if values is None or values == base_ring:
        return base_ring

    try:
        values = tuple(values)
    except TypeError:
        values = (values,)

    normalized = []
    for raw_value in values:
        coerced_value = base_ring(raw_value)
        if coerced_value not in normalized:
            normalized.append(coerced_value)
    return tuple(normalized)


def _normalize_coordinate_sets(base_ring, dimension, first_values=None,
                               coordinate_sets=None):
    r"""
    Normalize the allowed first coordinates on one side of the graph.

    EXAMPLES::

        sage: from sage.graphs.generators import luw_graphs as luw
        sage: luw._normalize_coordinate_sets(GF(3), 2, [0, 1])
        ((0, 1), Finite Field of size 3)
        sage: luw._normalize_coordinate_sets(GF(3), 2)
        (Finite Field of size 3, Finite Field of size 3)
    """
    if coordinate_sets is not None:
        first_values = tuple(coordinate_sets)[0]
    first_set = _normalize_coordinate_subset(base_ring, first_values)
    return (first_set,) + (base_ring,) * (dimension - 1)


def _normalize_vertex(luw_graph, vertex):
    r"""
    Normalize a user-given vertex and check that it lies in the domain.

    EXAMPLES::

        sage: from sage.graphs.generators import luw_graphs as luw
        sage: F = GF(3)
        sage: R = PolynomialRing(F, names=("p1", "l1"))
        sage: p1, l1 = R.gens()
        sage: D = luw.define_luw_graph(F, [p1*l1], A=[0, 1])
        sage: luw._normalize_vertex(D, ("P", (1, 2)))
        ('P', (1, 2))
    """
    side, coords = vertex
    coords = tuple(luw_graph.ring(value) for value in coords)
    if side not in ("P", "L") or len(coords) != luw_graph.dimension():
        raise ValueError('the vertex must have side "P" or "L" and the right number of '
                         "coordinates")

    first_coordinates = (luw_graph.point_first_coordinates() if side == "P"
                         else luw_graph.line_first_coordinates())
    if coords[0] not in first_coordinates:
        raise ValueError("the first coordinate is not in the allowed set for this side")
    return side, coords


def define_luw_graph(ring, equations, name=None, A=None, B=None,
                     point_coordinate_sets=None, line_coordinate_sets=None):
    r"""
    Return a validated algebraic descriptor for an LUW graph.

    The general construction is surveyed in [LW2026]_.

    INPUT:

    - ``ring`` -- finite commutative ring

    - ``equations`` -- nonempty iterable of Sage polynomials

      The polynomial variables are interpreted by their Sage generator names.
      Users should label Python variables in the same order as the polynomial
      ring generators, or bind them by name, for example with
      ``R.gens_dict()``. If ``p1`` is accidentally bound to the generator named
      ``p2``, the polynomial will be interpreted as using ``p2``.

    - ``name`` -- string (default: ``None``); name for the resulting
      definition

    - ``A`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the point side

    - ``B`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the line side

    - ``point_coordinate_sets`` -- iterable of coordinate sets or ``None``
      (default: ``None``)

    - ``line_coordinate_sets`` -- iterable of coordinate sets or ``None``
      (default: ``None``)

    OUTPUT:

    A :class:`LUWGraphDescriptor`.

    EXAMPLES::

        sage: from sage.graphs.generators.luw_graphs import define_luw_graph
        sage: F = GF(3)
        sage: R = PolynomialRing(F, names=("p1", "l1", "p2", "l2"))
        sage: p1, l1, p2, l2 = R.gens()
        sage: D = define_luw_graph(F, [p1*l1, p1*l2], A=[0, 1], B=[1, 2])
        sage: D.point_first_coordinates()
        (0, 1)
        sage: D.line_first_coordinates()
        (1, 2)
        sage: Z4 = Integers(4)
        sage: R = PolynomialRing(Z4, names=("p1", "l1"))
        sage: p1, l1 = R.gens()
        sage: define_luw_graph(Z4, [p1*l1]).order()
        32

    TESTS::

        sage: define_luw_graph(GF(3), [])
        Traceback (most recent call last):
        ...
        ValueError: give the list [f_2, ..., f_n]
        sage: R = PolynomialRing(GF(3), names=("p1", "l1", "p2", "l2"))
        sage: p1, l1, p2, l2 = R.gens()
        sage: define_luw_graph(GF(3), [p2])
        Traceback (most recent call last):
        ...
        ValueError: unknown variable p2; expected variables named p1, l1

    REFERENCES:

    - [LW2026]_
    """
    if not hasattr(ring, "is_finite") or not ring.is_finite():
        raise ValueError("the ring must be finite")
    if not hasattr(ring, "is_commutative") or not ring.is_commutative():
        raise ValueError("the ring must be commutative")

    equations = tuple(equations)
    if not equations:
        raise ValueError("give the list [f_2, ..., f_n]")

    dimension = len(equations) + 1
    canonical_ring = _canonical_polynomial_ring(ring, dimension)
    equations = tuple(_canonicalize_defining_polynomial(ring, polynomial, canonical_ring)
                      for polynomial in equations)

    for index, polynomial in enumerate(equations, start=2):
        _validate_defining_polynomial(polynomial, available_values=2 * (index - 1))

    point_coordinate_sets = _normalize_coordinate_sets(
        ring, dimension, A, point_coordinate_sets
    )
    line_coordinate_sets = _normalize_coordinate_sets(
        ring, dimension, B, line_coordinate_sets
    )

    return LUWGraphDescriptor(
        ring,
        equations,
        name or f"LUWGraph_{dimension}({ring})",
        point_coordinate_sets,
        line_coordinate_sets,
    )


def LUWGraph(ring, equations, name=None, A=None, B=None,
             point_coordinate_sets=None, line_coordinate_sets=None,
             immutable=False):
    r"""
    Build a graph directly from a given set of equations.

    This is the graph-returning constructor used by ``graphs.LUWGraph``.
    The construction is surveyed in [LW2026]_.

    INPUT:

    - ``ring`` -- finite commutative ring

    - ``equations`` -- nonempty iterable of Sage polynomials

    - ``name`` -- string (default: ``None``); graph name

    - ``A`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the point side

    - ``B`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the line side

    - ``point_coordinate_sets`` -- iterable of coordinate sets or ``None``
      (default: ``None``)

    - ``line_coordinate_sets`` -- iterable of coordinate sets or ``None``
      (default: ``None``)

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    OUTPUT:

    A Sage :class:`~sage.graphs.graph.Graph`.

    EXAMPLES::

        sage: F = GF(3)
        sage: R = PolynomialRing(F, names=("p1", "l1", "p2", "l2", "p3", "l3"))
        sage: p1, l1, p2, l2, p3, l3 = R.gens()
        sage: G = graphs.LUWGraph(F, [p1*l1, p1*l2, p3*l1])
        sage: G.order(), G.size()
        (162, 243)

    REFERENCES:

    - [LW2026]_
    """
    return define_luw_graph(
        ring,
        equations,
        name=name,
        A=A,
        B=B,
        point_coordinate_sets=point_coordinate_sets,
        line_coordinate_sets=line_coordinate_sets,
    ).graph(immutable=immutable)


def _evaluate_defining_polynomial(polynomial, previous_values):
    r"""
    Evaluate ``polynomial`` from an initial segment of LUW coordinates.

    The values in ``previous_values`` are interpreted in the canonical order
    ``p1, l1, p2, l2, ...``. If fewer values are supplied than the parent
    polynomial ring has generators, only the remaining tail is padded by zero.
    This helper is used after the defining equations have been validated as
    triangular.

    EXAMPLES::

        sage: from sage.graphs.generators import luw_graphs as luw
        sage: R = PolynomialRing(GF(3), names=("p1", "l1", "p2", "l2"))
        sage: p1, l1, p2, l2 = R.gens()
        sage: luw._evaluate_defining_polynomial(p1*l1, (2, 1))
        2
        sage: luw._evaluate_defining_polynomial(p1*l1 + p2, (2, 1, 1, 0))
        0
    """
    parent_ring = polynomial.parent()
    base_ring = parent_ring.base_ring()
    values = list(previous_values[: parent_ring.ngens()])
    values.extend(base_ring(0) for _ in range(parent_ring.ngens() - len(values)))
    return base_ring(polynomial(*values))


def _complete_line_from_point(luw_graph, point, first_line_value):
    r"""
    Reconstruct a line from a point and a first line coordinate.

    EXAMPLES::

        sage: from sage.graphs.generators import luw_graphs as luw
        sage: F = GF(3)
        sage: R = PolynomialRing(F, names=("p1", "l1"))
        sage: p1, l1 = R.gens()
        sage: D = luw.define_luw_graph(F, [p1*l1])
        sage: luw._complete_line_from_point(D, (2, 1), 2)
        ('L', (2, 0))
    """
    ring = luw_graph.ring
    line = [first_line_value]
    previous_values = [point[0], first_line_value]
    for index, fn in enumerate(luw_graph._equations, start=1):
        right_side = _evaluate_defining_polynomial(fn, previous_values)
        next_value = ring(right_side - point[index])
        line.append(next_value)
        previous_values.extend([point[index], next_value])
    return "L", tuple(line)


def _complete_point_from_line(luw_graph, line, first_point_value):
    r"""
    Reconstruct a point from a line and a first point coordinate.

    EXAMPLES::

        sage: from sage.graphs.generators import luw_graphs as luw
        sage: F = GF(3)
        sage: R = PolynomialRing(F, names=("p1", "l1"))
        sage: p1, l1 = R.gens()
        sage: D = luw.define_luw_graph(F, [p1*l1])
        sage: luw._complete_point_from_line(D, (2, 1), 2)
        ('P', (2, 0))
    """
    ring = luw_graph.ring
    point = [first_point_value]
    previous_values = [first_point_value, line[0]]
    for index, fn in enumerate(luw_graph._equations, start=1):
        right_side = _evaluate_defining_polynomial(fn, previous_values)
        next_value = ring(right_side - line[index])
        point.append(next_value)
        previous_values.extend([next_value, line[index]])
    return "P", tuple(point)


def _point_tuples(luw_graph):
    r"""
    Iterate over the point-coordinate tuples.

    EXAMPLES::

        sage: from itertools import islice
        sage: from sage.graphs.generators import luw_graphs as luw
        sage: F = GF(3)
        sage: R = PolynomialRing(F, names=("p1", "l1"))
        sage: p1, l1 = R.gens()
        sage: D = luw.define_luw_graph(F, [p1*l1], A=[0, 1])
        sage: tuple(islice(luw._point_tuples(D), 3))
        ((0, 0), (0, 1), (0, 2))
    """
    yield from product(*luw_graph.point_coordinate_sets)


def _line_tuples(luw_graph):
    r"""
    Iterate over the line-coordinate tuples.

    EXAMPLES::

        sage: from itertools import islice
        sage: from sage.graphs.generators import luw_graphs as luw
        sage: F = GF(3)
        sage: R = PolynomialRing(F, names=("p1", "l1"))
        sage: p1, l1 = R.gens()
        sage: D = luw.define_luw_graph(F, [p1*l1], B=[1, 2])
        sage: tuple(islice(luw._line_tuples(D), 3))
        ((1, 0), (1, 1), (1, 2))
    """
    yield from product(*luw_graph.line_coordinate_sets)


def _coordinate_polynomial_ring(ring, dimension):
    r"""
    Return indexed point and line generators in the canonical ring.

    EXAMPLES::

        sage: from sage.graphs.generators import luw_graphs as luw
        sage: point_vars, line_vars = luw._coordinate_polynomial_ring(GF(3), 3)
        sage: point_vars[1], line_vars[2]
        (p1, l2)
    """
    polynomial_ring = _canonical_polynomial_ring(ring, dimension)
    point_vars = [None] + list(polynomial_ring.gens()[0::2])
    line_vars = [None] + list(polynomial_ring.gens()[1::2])
    return point_vars, line_vars


def Akq(k, q, A=None, B=None, immutable=False):
    r"""
    Return the graph `A(k, q)` described in [LW2026]_.

    The point set and line set are `\GF{q}^k`. A point `(p_1, \ldots, p_k)`
    is adjacent to a line `(l_1, \ldots, l_k)` if

    .. MATH::

        p_i + l_i =
        \begin{cases}
        p_{i-1} l_1, & i \text{ even},\\
        p_1 l_{i-1}, & i \text{ odd},
        \end{cases}
        \qquad 2 \leq i \leq k.

    INPUT:

    - ``k`` -- integer at least `2`

    - ``q`` -- prime power

    - ``A`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the point side

    - ``B`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the line side

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    OUTPUT:

    A Sage :class:`~sage.graphs.graph.Graph`.

    EXAMPLES::

        sage: G = graphs.Akq(5, 3)
        sage: G.order(), G.size(), G.girth()
        (486, 729, 12)

    TESTS::

        sage: graphs.Akq(1, 3)
        Traceback (most recent call last):
        ...
        ValueError: k must be at least 2

    REFERENCES:

    - [LW2026]_
    """
    return define_Akq(k, q, A=A, B=B).graph(immutable=immutable)


def Dkq(k, q, A=None, B=None, immutable=False):
    r"""
    Return the graph `D(k, q)` described in [LW2026]_.

    The point set and line set are `\GF{q}^k`. A point `(p_1, \ldots, p_k)`
    is adjacent to a line `(l_1, \ldots, l_k)` if

    .. MATH::

        p_2 + l_2 = p_1 l_1,

    and, when `k \geq 3`,

    .. MATH::

        p_3 + l_3 = p_1 l_2.

    For `4 \leq i \leq k`, the remaining equations are

    .. MATH::

        p_i + l_i =
        \begin{cases}
        p_{i-2} l_1, & i \equiv 0, 1 \pmod 4,\\
        p_1 l_{i-2}, & i \equiv 2, 3 \pmod 4.
        \end{cases}

    INPUT:

    - ``k`` -- integer at least `2`

    - ``q`` -- prime power

    - ``A`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the point side

    - ``B`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the line side

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    OUTPUT:

    A Sage :class:`~sage.graphs.graph.Graph`.

    EXAMPLES::

        sage: G = graphs.Dkq(5, 3)
        sage: G.order(), G.size(), G.girth()
        (486, 729, 12)

    TESTS::

        sage: graphs.Dkq(1, 3)
        Traceback (most recent call last):
        ...
        ValueError: k must be at least 2

    REFERENCES:

    - [LW2026]_
    """
    return define_Dkq(k, q, A=A, B=B).graph(immutable=immutable)


def WengerGraph(m, q, A=None, B=None, immutable=False):
    r"""
    Return the Wenger graph `W_m(q)` described in [LW2026]_.

    The point set and line set are `\GF{q}^{m+1}`. A point
    `(p_1, \ldots, p_{m+1})` is adjacent to a line
    `(l_1, \ldots, l_{m+1})` if

    .. MATH::

        p_{i+1} + l_{i+1} = p_1 l_i,
        \qquad 1 \leq i \leq m.

    INPUT:

    - ``m`` -- positive integer

    - ``q`` -- prime power

    - ``A`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the point side

    - ``B`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the line side

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    OUTPUT:

    A Sage :class:`~sage.graphs.graph.Graph`.

    EXAMPLES::

        sage: G = graphs.WengerGraph(2, 3)
        sage: G.order(), G.size()
        (54, 81)

    TESTS::

        sage: graphs.WengerGraph(0, 3)
        Traceback (most recent call last):
        ...
        ValueError: m must be positive
        sage: graphs.WengerGraph(3, 3).order()
        doctest:...: UserWarning: WengerGraph(m, q) is disconnected when m > q - 1.
        162

    REFERENCES:

    - [LW2026]_
    """
    return define_WengerGraph(m, q, A=A, B=B).graph(immutable=immutable)


def define_Akq(k, q, A=None, B=None):
    r"""
    Return a descriptor for the graph `A(k, q)` described in [LW2026]_.

    The point set and line set are `\GF{q}^k`. A point `(p_1, \ldots, p_k)`
    is adjacent to a line `(l_1, \ldots, l_k)` if

    .. MATH::

        p_i + l_i =
        \begin{cases}
        p_{i-1} l_1, & i \text{ even},\\
        p_1 l_{i-1}, & i \text{ odd},
        \end{cases}
        \qquad 2 \leq i \leq k.

    INPUT:

    - ``k`` -- integer at least `2`

    - ``q`` -- prime power

    - ``A`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the point side

    - ``B`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the line side

    OUTPUT:

    A :class:`LUWGraphDescriptor`.

    EXAMPLES::

        sage: from sage.graphs.generators.luw_graphs import define_Akq
        sage: D = define_Akq(5, 3)
        sage: D.dimension(), D.order(), D.size()
        (5, 486, 729)
        sage: D.graph().order()
        486

    REFERENCES:

    - [LW2026]_
    """
    k = ZZ(k)
    if k < 2:
        raise ValueError("k must be at least 2")

    field = GF(q)
    point_vars, line_vars = _coordinate_polynomial_ring(field, k)
    equations = [
        point_vars[i - 1] * line_vars[1]
        if i % 2 == 0
        else point_vars[1] * line_vars[i - 1]
        for i in range(2, k + 1)
    ]
    return define_luw_graph(field, equations, name=f"A{k}{field.order()}",
                            A=A, B=B)


def define_Dkq(k, q, A=None, B=None):
    r"""
    Return a descriptor for the graph `D(k, q)` described in [LW2026]_.

    The point set and line set are `\GF{q}^k`. A point `(p_1, \ldots, p_k)`
    is adjacent to a line `(l_1, \ldots, l_k)` if

    .. MATH::

        p_2 + l_2 = p_1 l_1,

    and, when `k \geq 3`,

    .. MATH::

        p_3 + l_3 = p_1 l_2.

    For `4 \leq i \leq k`, the remaining equations are

    .. MATH::

        p_i + l_i =
        \begin{cases}
        p_{i-2} l_1, & i \equiv 0, 1 \pmod 4,\\
        p_1 l_{i-2}, & i \equiv 2, 3 \pmod 4.
        \end{cases}

    INPUT:

    - ``k`` -- integer at least `2`

    - ``q`` -- prime power

    - ``A`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the point side

    - ``B`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the line side

    OUTPUT:

    A :class:`LUWGraphDescriptor`.

    EXAMPLES::

        sage: from sage.graphs.generators.luw_graphs import define_Dkq
        sage: D = define_Dkq(5, 3)
        sage: D.dimension(), D.order(), D.size()
        (5, 486, 729)
        sage: D.graph().size()
        729

    REFERENCES:

    - [LW2026]_
    """
    k = ZZ(k)
    if k < 2:
        raise ValueError("k must be at least 2")

    field = GF(q)
    point_vars, line_vars = _coordinate_polynomial_ring(field, k)
    equations = [point_vars[1] * line_vars[1]]
    if k >= 3:
        equations.append(point_vars[1] * line_vars[2])
    for i in range(4, k + 1):
        if i % 4 in (0, 1):
            equations.append(point_vars[i - 2] * line_vars[1])
        else:
            equations.append(point_vars[1] * line_vars[i - 2])
    return define_luw_graph(field, equations, name=f"D{k}{field.order()}",
                            A=A, B=B)


def define_WengerGraph(m, q, A=None, B=None):
    r"""
    Return a descriptor for the Wenger graph `W_m(q)` described in [LW2026]_.

    The point set and line set are `\GF{q}^{m+1}`. A point
    `(p_1, \ldots, p_{m+1})` is adjacent to a line
    `(l_1, \ldots, l_{m+1})` if

    .. MATH::

        p_{i+1} + l_{i+1} = p_1 l_i,
        \qquad 1 \leq i \leq m.

    INPUT:

    - ``m`` -- positive integer

    - ``q`` -- prime power

    - ``A`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the point side

    - ``B`` -- iterable, single value, or ``None`` (default: ``None``);
      allowed first coordinates on the line side

    OUTPUT:

    A :class:`LUWGraphDescriptor`.

    EXAMPLES::

        sage: from sage.graphs.generators.luw_graphs import define_WengerGraph
        sage: D = define_WengerGraph(2, 3)
        sage: D.dimension(), D.order(), D.size()
        (3, 54, 81)
        sage: D.graph().order()
        54

    REFERENCES:

    - [LW2026]_
    """
    m = ZZ(m)
    if m < 1:
        raise ValueError("m must be positive")

    field = GF(q)
    if m > field.order() - 1:
        warn(
            "WengerGraph(m, q) is disconnected when m > q - 1.",
            stacklevel=2,
        )
    point_vars, line_vars = _coordinate_polynomial_ring(field, m + 1)
    equations = [point_vars[1] * line_vars[i] for i in range(1, m + 1)]
    return define_luw_graph(field, equations, name=f"W_{m}({q})", A=A, B=B)
