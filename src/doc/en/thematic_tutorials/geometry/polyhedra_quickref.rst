.. -*- coding: utf-8 -*-
.. linkall

.. _polyhedra_quickref:

=====================================
Quick reference for polyhedra in Sage
=====================================

.. MODULEAUTHOR:: Jean-Philippe Labbé <labbe@math.fu-berlin.de>
                  Vincent Delecroix <vincent.delecroix@u-bordeaux.fr>

List of Polyhedron methods
==========================

**H and V-representation**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.base_ring` | ring on which the polyhedron is defined
    :meth:`~sage.geometry.polyhedron.base1.Polyhedron_base1.ambient_space` | ambient vector space or free module
    :meth:`~sage.geometry.polyhedron.base1.Polyhedron_base1.Hrepresentation_space` | vector space or free module used for the vectors of the H-representation
    :meth:`~sage.geometry.polyhedron.base1.Polyhedron_base1.Vrepresentation_space` | vector space or free module used for the vectors of the V-representation
    :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.n_Hrepresentation` | number of elements in the H-representation (sum of the number of equations and inequalities)
    :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.n_Vrepresentation` | number of elements in the V-representation (sum of vertices, rays and lines)
    :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.n_equations` | number of equations
    :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.n_inequalities` | number of inequalities
    :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.n_vertices` | number of vertices
    :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.n_rays` | number of rays
    :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.n_lines` | number of lines
    :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.n_facets` | number of facets

**Polyhedron boolean properties:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base1.Polyhedron_base1.is_empty` | tests emptyness
    :meth:`~sage.geometry.polyhedron.base1.Polyhedron_base1.is_universe` | tests whether a polyhedra is the whole ambient space
    :meth:`~sage.geometry.convex_set.ConvexSet_base.is_full_dimensional` | tests if the polyhedron has the same dimension as the ambient space
    :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.is_combinatorially_isomorphic` | tests whether two polyhedra are combinatorially isomorphic
    :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.is_compact` | tests compactness, or boundedness of a polyhedron
    :meth:`~sage.geometry.polyhedron.base2.Polyhedron_base2.is_lattice_polytope` | tests whether a polyhedron is a lattice polytope
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_inscribed` | tests whether the polyhedron is inscribed in a sphere
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_minkowski_summand` | tests if the polyhedron can be used to produce another given polyhedron using a Minkowski sum.
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.is_neighborly` | tests whether the polyhedron has full skeleton until half of the dimension (or up to a certain dimension)
    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.is_reflexive` | tests if the polar of a lattice polytope is also a lattice polytope (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.is_simple` |  checks whether the degree of all vertices is equal to the dimension of the polytope
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.is_simplex` | test whether a polytope is a simplex
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.is_simplicial` |  checks whether all faces of the polyhedron are simplices
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.is_lawrence_polytope` |  tests whether self is a Lawrence polytope
    :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.is_self_dual` |  tests whether the polytope is self-dual
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.is_pyramid` | test whether the polytope is a pyramid over one of its facets
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.is_bipyramid` | test whether the polytope is combinatorially equivalent to a bipyramid over some polytope
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.is_prism` | test whether the polytope is combinatorially equivalent to a prism of some polytope

**Enumerative properties**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base1.Polyhedron_base1.ambient_dim` |  the dimension of the ambient vector space
    :meth:`~sage.geometry.polyhedron.base1.Polyhedron_base1.dim` |  the dimension of the polytope
    :meth:`~sage.geometry.polyhedron.base1.Polyhedron_base1.dimension` |  alias of dim
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.f_vector` |  the `f`-vector (number of faces of each dimension)
    :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.flag_f_vector` |  the flag-`f`-vector (number of chains of faces)
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.neighborliness` | highest cardinality for which all `k`-subsets of the vertices are faces of the polyhedron
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.simpliciality` | highest cardinality for which all `k`-faces are simplices
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.simplicity` | highest cardinality for which the polar is `k`-simplicial

**Implementation properties**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.backend` | gives the backend used
    :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.base_ring` | gives the base ring used
    :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.change_ring` | changes the base ring

**Transforming polyhedra**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.minkowski_sum` | Minkowski sum of two polyhedra
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.minkowski_difference` | Minkowski difference of two polyhedra
    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.minkowski_decompositions` | Minkowski decomposition (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.product` | cartesian product of two polyhedra
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.intersection` | intersection of two polyhedra
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.join` | join of two polyhedra
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.convex_hull` | convex hull of the union of two polyhedra
    :meth:`~sage.geometry.polyhedron.base6.Polyhedron_base6.affine_hull_projection` | constructs an affinely equivalent full-dimensional polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.barycentric_subdivision` | constructs a geometric realization of the barycentric subdivision
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.dilation` |  scalar dilation
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.face_truncation` | truncates a specific face
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.face_split` | returns the face splitting of a face of self
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.one_point_suspension` | the one-point suspension over a vertex of self (face splitting of a vertex)
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.stack` | stack a face of the polyhedron
    :meth:`~sage.geometry.polyhedron.base2.Polyhedron_base2.lattice_polytope` | returns an encompassing lattice polytope.
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.polar` | returns the polar of a polytope (needs to be compact)
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.prism` | prism over a polyhedron (increases both the dimension of the polyhedron and the dimension of the ambient space)
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.pyramid` | pyramid over a polyhedron (increases both the dimension of the polyhedron and the dimension of the ambient space)
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.bipyramid` | bipyramid over a polyhedron (increases both the dimension of the polyhedron and the dimension of the ambient)
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.translation` | translates by a given vector
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.truncation` | truncates all vertices simultaneously
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.lawrence_extension` | returns the Lawrence extension of self on a given point
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.lawrence_polytope` | returns the Lawrence polytope of self
    :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.wedge` | returns the wedge over a face of self

**Combinatorics**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.combinatorial_polyhedron` | the combinatorial polyhedron
    :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.face_lattice` | the face lattice
    :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.hasse_diagram` | the hasse diagram
    :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.combinatorial_automorphism_group` | the automorphism group of the underlying combinatorial polytope
    :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.graph`, :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.vertex_graph` | underlying graph
    :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.vertex_digraph` | digraph (orientation of edges determined by a linear form)
    :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.vertex_facet_graph` | bipartite digraph given vertex-facet adjacency
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.adjacency_matrix` | adjacency matrix
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.incidence_matrix` | incidence matrix
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.slack_matrix` | slack matrix
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.facet_adjacency_matrix` | adjacency matrix of the facets
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.vertex_adjacency_matrix` | adjacency matrix of the vertices

**Integral points**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.ehrhart_polynomial` | the Ehrhart polynomial for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`
    :meth:`~sage.geometry.polyhedron.base_QQ.Polyhedron_QQ.ehrhart_polynomial` | the Ehrhart polynomial for :class:`Polyhedron over QQ <sage.geometry.polyhedron.base_QQ.Polyhedron_QQ>`
    :meth:`~sage.geometry.polyhedron.base_QQ.Polyhedron_QQ.ehrhart_quasipolynomial` | the Ehrhart quasipolynomial for :class:`Polyhedron over QQ <sage.geometry.polyhedron.base_QQ.Polyhedron_QQ>`
    :meth:`~sage.geometry.polyhedron.base2.Polyhedron_base2.h_star_vector` | the `h^*`-vector for polytopes with integral vertices
    :meth:`~sage.geometry.polyhedron.base2.Polyhedron_base2.integral_points` | list of integral points
    :meth:`~sage.geometry.polyhedron.base2.Polyhedron_base2.integral_points_count` | number of integral points
    :meth:`~sage.geometry.polyhedron.base2.Polyhedron_base2.get_integral_point` | get the i-th integral point without computing all interior lattice points
    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.has_IP_property` | checks whether the origin is an interior lattice point and compactness (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)
    :meth:`~sage.geometry.polyhedron.base2.Polyhedron_base2.random_integral_point` | get a random integral point


**Getting related geometric objects**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base6.Polyhedron_base6.affine_hull` | returns the smallest affine subspace containing the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.boundary_complex` | returns the boundary complex of simplicial compact polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.center` | returns the average of the vertices of the polyhedron
    :meth:`~sage.geometry.polyhedron.base7.Polyhedron_base7.centroid` | returns the center of the mass
    :meth:`~sage.geometry.polyhedron.base1.Polyhedron_base1.representative_point` | returns the sum of the center and the rays
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.a_maximal_chain` | returns a maximal chain of faces
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.face_fan` | returns the fan spanned by the faces of the polyhedron
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.face_generator` | a generator over the faces
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.faces` | the list of faces
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.facets` | the list of facets
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.join_of_Vrep`, :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.least_common_superface_of_Vrep` | smallest face containing specified Vrepresentatives
    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.meet_of_Hrep`, :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.greatest_common_subface_of_Hrep` | largest face contained in specified Hrepresentatives
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.normal_fan` | returns the fan spanned by the normals of the supporting hyperplanes of the polyhedron
    :meth:`~sage.geometry.polyhedron.base6.Polyhedron_base6.gale_transform` | returns the (affine) Gale transform of the vertices of the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.hyperplane_arrangement` | returns the hyperplane arrangement given by the defining facets of the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.to_linear_program` | transform the polyhedra into a Linear Program
    :meth:`~sage.geometry.polyhedron.base7.Polyhedron_base7.triangulate` | returns a triangulation of the polyhedron
    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.fibration_generator` | returns an iterator of the fibrations of the lattice polytope (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)

**Other**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |


    :meth:`~sage.geometry.polyhedron.base3.Polyhedron_base3.bounded_edges` | generator for bounded edges
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.bounding_box` | returns the vertices of an encompassing cube
    :meth:`~sage.geometry.polyhedron.base1.Polyhedron_base1.contains` | tests whether the polyhedron contains a vector
    :meth:`~sage.geometry.polyhedron.base1.Polyhedron_base1.interior_contains` | tests whether the polyhedron contains a vector in its interior using the ambient topology
    :meth:`~sage.geometry.polyhedron.base1.Polyhedron_base1.relative_interior_contains` | tests whether the polyhedron contains a vector in its relative interior
    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.find_translation` | returns the translation vector between two translation of two polyhedron (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)
    :meth:`~sage.geometry.polyhedron.base7.Polyhedron_base7.integrate` | computes the integral of a polynomial over the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.radius` | returns the radius of the smallest sphere containing the polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.radius_square` | returns the square of the radius of the smallest sphere containing the polyhedron
    :meth:`~sage.geometry.polyhedron.base7.Polyhedron_base7.volume` | computes different volumes of the polyhedron
    :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.restricted_automorphism_group` | returns the restricted automorphism group
    :meth:`~sage.geometry.polyhedron.ppl_lattice_polytope.LatticePolytope_PPL_class.lattice_automorphism_group` | returns the lattice automorphism group. Only for :class:`PPL Lattice Polytope <sage.geometry.polyhedron.ppl_lattice_polytope.LatticePolytope_PPL_class>`
