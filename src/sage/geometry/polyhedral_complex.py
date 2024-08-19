# sage.doctest: needs sage.graphs
r"""
Finite polyhedral complexes

This module implements the basic structure of finite polyhedral complexes.
For more information, see :class:`PolyhedralComplex`.

AUTHORS:

- Yuan Zhou (2021-05): initial implementation

List of PolyhedralComplex methods
---------------------------------

**Maximal cells and cells**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~PolyhedralComplex.maximal_cells` | Return the dictionary of the maximal cells in this polyhedral complex.
    :meth:`~PolyhedralComplex.maximal_cell_iterator` | Return an iterator over maximal cells in this polyhedral complex.
    :meth:`~PolyhedralComplex.maximal_cells_sorted` | Return the sorted list of all maximal cells in this polyhedral complex.
    :meth:`~PolyhedralComplex.n_maximal_cells` | List the maximal cells of dimension `n` in this polyhedral complex.
    :meth:`~PolyhedralComplex._n_maximal_cells_sorted` | Return the sorted list of maximal cells of dim `n` in this complex.
    :meth:`~PolyhedralComplex.is_maximal_cell` | Return ``True`` if the given cell is a maximal cell in this complex.
    :meth:`~PolyhedralComplex.cells` | Return the dictionary of the cells in this polyhedral complex.
    :meth:`~PolyhedralComplex.cell_iterator` | Return an iterator over cells in this polyhedral complex.
    :meth:`~PolyhedralComplex.cells_sorted` | Return the sorted list of all cells in this polyhedral complex.
    :meth:`~sage.topology.cell_complex.GenericCellComplex.n_cells` | List the cells of dimension `n` in this polyhedral complex.
    :meth:`~PolyhedralComplex._n_cells_sorted` | Return the sorted list of `n`-cells in this polyhedral complex.
    :meth:`~PolyhedralComplex.is_cell` | Return ``True`` if the given cell is in this polyhedral complex.
    :meth:`~PolyhedralComplex.face_poset` | Return the poset of nonempty cells in the polyhedral complex.
    :meth:`~PolyhedralComplex.relative_boundary_cells` | List the maximal cells on the boundary of the polyhedral complex.

**Properties of the polyhedral complex**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~PolyhedralComplex.dimension` | Return the dimension of the polyhedral complex.
    :meth:`~PolyhedralComplex.ambient_dimension` | Return the ambient dimension of the polyhedral complex.
    :meth:`~PolyhedralComplex.is_pure` | Return ``True`` if the polyhedral complex is pure.
    :meth:`~PolyhedralComplex.is_full_dimensional` | Return ``True`` if the polyhedral complex is full dimensional.
    :meth:`~PolyhedralComplex.is_compact` | Return ``True`` if the polyhedral complex is bounded.
    :meth:`~PolyhedralComplex.is_connected` | Return ``True`` if the polyhedral complex is connected.
    :meth:`~PolyhedralComplex.is_subcomplex` | Return ``True`` if this complex is a subcomplex of the other.
    :meth:`~PolyhedralComplex.is_convex` | Return ``True`` if the polyhedral complex is convex.
    :meth:`~PolyhedralComplex.is_mutable` | Return ``True`` if the polyhedral complex is mutable.
    :meth:`~PolyhedralComplex.is_immutable` | Return ``True`` if the polyhedral complex is not mutable.
    :meth:`~PolyhedralComplex.is_simplicial_complex` | Return ``True`` if the polyhedral complex is a simplicial complex.
    :meth:`~PolyhedralComplex.is_polyhedral_fan` | Return ``True`` if the polyhedral complex is a fan.
    :meth:`~PolyhedralComplex.is_simplicial_fan` | Return ``True`` if the polyhedral complex is a simplicial fan.

**New polyhedral complexes from old ones**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~PolyhedralComplex.connected_component` | Return the connected component containing a cell as a subcomplex.
    :meth:`~PolyhedralComplex.connected_components` | Return the connected components of this polyhedral complex.
    :meth:`~PolyhedralComplex.n_skeleton` | Return the `n`-skeleton of this polyhedral complex.
    :meth:`~PolyhedralComplex.stratify` | Return the (pure) subcomplex formed by the maximal cells of dim `n` in this complex.
    :meth:`~PolyhedralComplex.boundary_subcomplex` | Return the boundary subcomplex of this polyhedral complex.
    :meth:`~PolyhedralComplex.product` | Return the (Cartesian) product of this polyhedral complex with another one.
    :meth:`~PolyhedralComplex.disjoint_union` | Return the disjoint union of this polyhedral complex with another one.
    :meth:`~PolyhedralComplex.union` | Return the union of this polyhedral complex with another one.
    :meth:`~PolyhedralComplex.join` | Return the join of this polyhedral complex with another one.
    :meth:`~PolyhedralComplex.subdivide` | Return a new polyhedral complex (with option ``make_simplicial``) subdividing this one.

**Update polyhedral complex**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~PolyhedralComplex.set_immutable` | Make this polyhedral complex immutable.
    :meth:`~PolyhedralComplex.add_cell` | Add a cell to this polyhedral complex.
    :meth:`~PolyhedralComplex.remove_cell` | Remove a cell from this polyhedral complex.

**Miscellaneous**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~PolyhedralComplex.plot` | Return a Graphic object showing the plot of polyhedral complex.
    :meth:`~PolyhedralComplex.graph` | Return a directed graph corresponding to the 1-skeleton of this polyhedral complex, given that it is bounded.
    :meth:`~PolyhedralComplex.union_as_polyhedron` | Return a ``Polyhedron`` which is the union of cells in this polyhedral complex, given that it is convex.

Classes and functions
---------------------
"""

# ****************************************************************************
#       Copyright (C) 2021 Yuan Zhou <yuan.zhou@uky.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy

import sage.geometry.abc

from sage.topology.cell_complex import GenericCellComplex
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.graphs.graph import Graph
from sage.combinat.posets.posets import Poset
from sage.combinat.subset import powerset


class PolyhedralComplex(GenericCellComplex):
    r"""
    A polyhedral complex.

    A **polyhedral complex** `PC` is a collection of polyhedra in a certain
    ambient space `\RR^n` such that the following hold.

    - If a polyhedron `P` is in `PC`, then all the faces of `P` are in `PC`.

    - If polyhedra `P` and `Q` are in `PC`, then `P \cap Q` is either empty
      or a face of both `P` and `Q`.

    In this context, a "polyhedron" means the geometric realization
    of a polyhedron. This is in contrast to :mod:`simplicial complex
    <sage.topology.simplicial_complex>`, whose cells are abstract simplices.
    The concept of a polyhedral complex generalizes that of a **geometric**
    simplicial complex.

    .. NOTE::

       This class derives from
       :class:`~sage.topology.cell_complex.GenericCellComplex`, and so
       inherits its methods.  Some of those methods are not listed here;
       see the :mod:`Generic Cell Complex <sage.topology.cell_complex>`
       page instead.

    INPUT:

    - ``maximal_cells`` -- list, tuple, or dictionary (indexed by
      dimension) of cells of the Complex. Each cell is of class
      :class:`Polyhedron` of the same ambient dimension. To set up a
      :class:PolyhedralComplex, it is sufficient to provide the maximal
      faces. Use keyword argument ``partial=True`` to set up a partial
      polyhedral complex, which is a subset of the faces (viewed as
      relatively open) of a polyhedral complex that is not necessarily
      closed under taking intersection.

    - ``maximality_check`` -- boolean (default: ``True``);
      if ``True``, then the constructor checks that each given
      maximal cell is indeed maximal, and ignores those that are not

    - ``face_to_face_check`` -- boolean (default: ``False``);
      if ``True``, then the constructor checks whether the cells
      are face-to-face, and it raises a :exc:`ValueError` if they are not

    - ``is_mutable`` and ``is_immutable`` -- boolean (default: ``True`` and
      ``False`` respectively); set ``is_mutable=False`` or ``is_immutable=True``
      to make this polyhedral complex immutable

    - ``backend`` -- string (optional); the name of the backend used for
      computations on Sage polyhedra; if it is not given, then each cell has
      its own backend; otherwise it must be one of the following:

      * ``'ppl'`` -- the Parma Polyhedra Library

      * ``'cdd'`` -- CDD

      * ``'normaliz'`` -- normaliz

      * ``'polymake'`` -- polymake

      * ``'field'`` -- a generic Sage implementation

    - ``ambient_dim`` -- integer (optional); used to set up an empty
      complex in the intended ambient space

    EXAMPLES::

        sage: pc = PolyhedralComplex([
        ....:         Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1/7, 2/7)]),
        ....:         Polyhedron(vertices=[(1/7, 2/7), (0, 0), (0, 1/4)])])
        sage: [p.Vrepresentation() for p in pc.cells_sorted()]
        [(A vertex at (0, 0), A vertex at (0, 1/4), A vertex at (1/7, 2/7)),
         (A vertex at (0, 0), A vertex at (1/3, 1/3), A vertex at (1/7, 2/7)),
         (A vertex at (0, 0), A vertex at (0, 1/4)),
         (A vertex at (0, 0), A vertex at (1/7, 2/7)),
         (A vertex at (0, 0), A vertex at (1/3, 1/3)),
         (A vertex at (0, 1/4), A vertex at (1/7, 2/7)),
         (A vertex at (1/3, 1/3), A vertex at (1/7, 2/7)),
         (A vertex at (0, 0),),
         (A vertex at (0, 1/4),),
         (A vertex at (1/7, 2/7),),
         (A vertex at (1/3, 1/3),)]
        sage: pc.plot()                                                                 # needs sage.plot
        Graphics object consisting of 10 graphics primitives
        sage: pc.is_pure()
        True
        sage: pc.is_full_dimensional()
        True
        sage: pc.is_compact()
        True
        sage: pc.boundary_subcomplex()
        Polyhedral complex with 4 maximal cells
        sage: pc.is_convex()
        True
        sage: pc.union_as_polyhedron().Hrepresentation()
        (An inequality (1, -4) x + 1 >= 0,
         An inequality (-1, 1) x + 0 >= 0,
         An inequality (1, 0) x + 0 >= 0)
        sage: pc.face_poset()
        Finite poset containing 11 elements
        sage: pc.is_connected()
        True
        sage: pc.connected_component() == pc
        True

    TESTS:

    Check that non-maximal cells are ignored if ``maximality_check=True``::

        sage: pc = PolyhedralComplex([
        ....:         Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
        ....:         Polyhedron(vertices=[(1, 2), (0, 0)]) ])
        sage: pc.maximal_cells()
        {2: {A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices}}

    Check that non face-to-face can be detected::

        sage: PolyhedralComplex([
        ....:         Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
        ....:         Polyhedron(vertices=[(2, 2), (0, 0)]) ],
        ....:         face_to_face_check=True)
        Traceback (most recent call last):
        ...
        ValueError: the given cells are not face-to-face

    Check that all the cells must have the same ambient dimension::

        sage: PolyhedralComplex([
        ....:         Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
        ....:         Polyhedron(vertices=[[2], [0]]) ])
        Traceback (most recent call last):
        ...
        ValueError: the given cells are not polyhedra in the same ambient space

    Check that backend is passed to all the cells::

        sage: P = Polyhedron(vertices=[(0, 0), (1, 1)])
        sage: P.backend()
        'ppl'
        sage: pc = PolyhedralComplex([P], backend='cdd')
        sage: Q = pc.maximal_cells_sorted()[0]
        sage: Q.backend()
        'cdd'
    """
    def __init__(self, maximal_cells=None, backend=None, maximality_check=True,
                 face_to_face_check=False, is_mutable=True, is_immutable=False,
                 ambient_dim=None):
        r"""
        Define a PolyhedralComplex.

        See ``PolyhedralComplex`` for more information.

        EXAMPLES::

            sage: pc = PolyhedralComplex([Polyhedron(vertices=[(1, 1), (0, 0)])])
            sage: pc
            Polyhedral complex with 1 maximal cell
            sage: TestSuite(pc).run()
        """
        self._backend = backend
        if maximal_cells is None:
            cells_dict = {}
        elif isinstance(maximal_cells, (list, tuple)):
            if backend:
                maximal_cells = [p.base_extend(p.base_ring(), backend)
                                 for p in maximal_cells]
            cells_dict = cells_list_to_cells_dict(maximal_cells)
        elif isinstance(maximal_cells, dict):
            cells_dict = {}
            for (k, l) in maximal_cells.items():
                if backend:
                    cells_dict[k] = set([p.base_extend(p.base_ring(), backend)
                                        for p in l])
                else:
                    cells_dict[k] = set(l)
        else:
            raise ValueError("the maximal cells are not given in correct form")
        if not cells_dict:
            self._dim = -1
            if ambient_dim is None:
                ambient_dim = -1
        else:
            self._dim = max(cells_dict.keys())
            if ambient_dim is None:
                ambient_dim = next(iter(cells_dict[self._dim])).ambient_dim()
        self._ambient_dim = ambient_dim
        self._maximal_cells = cells_dict
        if not all((isinstance(cell, sage.geometry.abc.Polyhedron) and
                   cell.ambient_dim() == self._ambient_dim)
                   for cell in self.maximal_cell_iterator()):
            raise ValueError("the given cells are not polyhedra " +
                             "in the same ambient space")
        # initialize the attributes
        self._is_convex = None
        self._polyhedron = None
        self._maximal_cells_sorted = None    # needed for hash
        self._cells = None
        self._face_poset = None

        if maximality_check:
            self.cells()    # compute self._cells and self._face_poset
            self._maximal_cells = cells_list_to_cells_dict(
                                      self._face_poset.maximal_elements())
        if face_to_face_check:
            poset = self.face_poset()
            maximal_cells = poset.maximal_elements()    # a list
            for i in range(len(maximal_cells)):
                p = maximal_cells[i]
                for j in range(i, len(maximal_cells)):
                    q = maximal_cells[j]
                    r = p.intersection(q)
                    if not (r.is_empty() or (r in poset) and
                            poset.is_gequal(p, r) and poset.is_gequal(q, r)):
                        raise ValueError("the given cells are not face-to-face")
        self._is_immutable = False
        if not is_mutable or is_immutable:
            self.set_immutable()

    def cells(self, subcomplex=None):
        """
        The cells of this polyhedral complex, in the form of a dictionary:
        the keys are integers, representing dimension, and the value
        associated to an integer `d` is the set of `d`-cells.

        INPUT:

        - ``subcomplex`` -- (optional) if a subcomplex is given then
          return the cells which are **not** in this subcomplex

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: list(pc.cells().keys())
            [2, 1, 0]
        """
        if subcomplex is not None:
            raise NotImplementedError("providing subcomplex is not implemented")
        if self._cells is not None:
            return self._cells
        maximal_cells = self.maximal_cells()
        cells = {}
        covers = {}
        for k in range(self._dim, -1, -1):
            if k in maximal_cells:
                if k not in cells:
                    cells[k] = set([])
                cells[k].update(maximal_cells[k])
            if k in cells:
                for cell in cells[k]:
                    if cell not in covers:
                        covers[cell] = []
                    for facet in cell.facets():
                        p = facet.as_polyhedron()
                        if p not in covers:
                            covers[p] = []
                        covers[p].append(cell)
                        if (k-1) not in cells:
                            cells[k-1] = set([])
                        cells[k-1].add(p)
        self._face_poset = Poset(covers)
        self._cells = cells
        return self._cells

    def cell_iterator(self, increasing=True):
        """
        An iterator for the cells in this polyhedral complex.

        INPUT:

        - ``increasing`` -- boolean (default: ``True``); if ``True``, return
          cells in increasing order of dimension, thus starting with the
          zero-dimensional cells; otherwise it returns cells in decreasing
          order of dimension

        .. NOTE::

            Among the cells of a fixed dimension, there is no sorting.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: len(list(pc.cell_iterator()))
            11
        """
        cells = self.cells()
        dim_index = range(0, self.dimension() + 1)
        if not increasing:
            dim_index = reversed(dim_index)
        for d in dim_index:
            if d in cells:
                yield from cells[d]

    def _n_cells_sorted(self, n, subcomplex=None):
        """
        Sorted list of cells of dimension ``n`` of this polyhedral complex.

        INPUT:

        - ``n`` -- nonnegative integer; the dimension
        - ``subcomplex`` -- (optional) if a subcomplex is given then
          return the cells which are **not** in this subcomplex

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: [p.Vrepresentation() for p in pc._n_cells_sorted(1)]
            [(A vertex at (0, 0), A vertex at (0, 2)),
             (A vertex at (0, 0), A vertex at (1, 1)),
             (A vertex at (0, 0), A vertex at (1, 2)),
             (A vertex at (0, 2), A vertex at (1, 2)),
             (A vertex at (1, 1), A vertex at (1, 2))]
            sage: pc._n_cells_sorted(3)
            []
        """
        n_cells = self.n_cells(n, subcomplex)
        return sorted(n_cells,
                      key=lambda p: (p.vertices(), p.rays(), p.lines()))

    def cells_sorted(self, subcomplex=None):
        """
        The sorted list of the cells of this polyhedral complex
        in non-increasing dimensions.

        INPUT:

        - ``subcomplex`` -- (optional) if a subcomplex is given then
          return the cells which are **not** in this subcomplex

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: len(pc.cells_sorted())
            11
            sage: pc.cells_sorted()[0].Vrepresentation()
            (A vertex at (0, 0), A vertex at (0, 2), A vertex at (1, 2))
        """
        cells = []
        for n in range(self._dim, -1, -1):
            cells += self._n_cells_sorted(n, subcomplex)
        return cells

    def maximal_cells(self):
        """
        The maximal cells of this polyhedral complex, in the form of a
        dictionary: the keys are integers, representing dimension, and the
        value associated to an integer `d` is the set of `d`-maximal cells.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: len(pc.maximal_cells()[2])
            2
            sage: 1 in pc.maximal_cells()
            False

        Wrong answer due to ``maximality_check=False``::

            sage: pc_invalid = PolyhedralComplex([p1, p2, p3],
            ....:              maximality_check=False)
            sage: len(pc_invalid.maximal_cells()[1])
            1
        """
        return self._maximal_cells

    def maximal_cell_iterator(self, increasing=False):
        r"""
        An iterator for the maximal cells in this polyhedral complex.

        INPUT:

        - ``increasing`` -- boolean (default: ``False``); if ``True``, return
          maximal cells in increasing order of dimension.
          Otherwise it returns cells in decreasing order of dimension.

        .. NOTE::

            Among the cells of a fixed dimension, there is no sorting.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: len(list(pc.maximal_cell_iterator()))
            2

        Wrong answer due to ``maximality_check=False``::

            sage: pc_invalid = PolyhedralComplex([p1, p2, p3],
            ....:              maximality_check=False)
            sage: len(list(pc_invalid.maximal_cell_iterator()))
            3
        """
        maximal_cells = self.maximal_cells()
        dim_index = range(-1, self.dimension() + 1)
        if not increasing:
            dim_index = reversed(dim_index)
        for d in dim_index:
            if d in maximal_cells:
                yield from maximal_cells[d]

    def n_maximal_cells(self, n):
        r"""
        List of maximal cells of dimension ``n`` of this polyhedral complex.

        INPUT:

        - ``n`` -- nonnegative integer; the dimension

        .. NOTE::

            The resulting list need not be sorted. If you want a sorted
            list of `n`-cells, use :meth:`_n_maximal_cells_sorted`.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: len(pc.n_maximal_cells(2))
            2
            sage: len(pc.n_maximal_cells(1))
            0

        Wrong answer due to ``maximality_check=False``::

            sage: pc_invalid = PolyhedralComplex([p1, p2, p3],
            ....:              maximality_check=False)
            sage: len(pc_invalid.n_maximal_cells(1))
            1
        """
        if n in self.maximal_cells():
            return list(self.maximal_cells()[n])
        else:
            return []

    def _n_maximal_cells_sorted(self, n):
        """
        Sorted list of maximal cells of dimension ``n`` of this polyhedral
        complex.

        INPUT:

        - ``n`` -- nonnegative integer; the dimension

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: pc._n_maximal_cells_sorted(2)[0].vertices_list()
            [[0, 0], [0, 2], [1, 2]]
        """
        n_maximal_cells = self.n_maximal_cells(n)
        return sorted(n_maximal_cells,
                      key=lambda p: (p.vertices(), p.rays(), p.lines()))

    def maximal_cells_sorted(self):
        """
        Return the sorted list of the maximal cells of this polyhedral complex
        by non-increasing dimensions.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: [p.vertices_list() for p in pc.maximal_cells_sorted()]
            [[[0, 0], [0, 2], [1, 2]], [[0, 0], [1, 1], [1, 2]]]
        """
        if self._maximal_cells_sorted is None:
            maximal_cells = []
            for n in range(self._dim, -1, -1):
                maximal_cells += self._n_maximal_cells_sorted(n)
            self._maximal_cells_sorted = maximal_cells
        return self._maximal_cells_sorted

    def is_maximal_cell(self, c):
        """
        Return whether the given cell ``c`` is a maximal cell of ``self``.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: pc.is_maximal_cell(p1)
            True
            sage: pc.is_maximal_cell(p3)
            False

        Wrong answer due to ``maximality_check=False``::

            sage: pc_invalid = PolyhedralComplex([p1, p2, p3],
            ....:              maximality_check=False)
            sage: pc_invalid.is_maximal_cell(p3)
            True
        """
        d = c.dimension()
        # return (c in self.n_maximal_cells(d)) # use set instead of list
        return (d in self.maximal_cells()) and (c in self.maximal_cells()[d])

    def is_cell(self, c):
        """
        Return whether the given cell ``c`` is a cell of ``self``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2])
            sage: pc.is_cell(p3)
            True
            sage: pc.is_cell(Polyhedron(vertices=[(0, 0)]))
            True
        """
        d = c.dimension()
        return (d in self.cells()) and (c in self.cells()[d])

    def dimension(self):
        """
        The dimension of this cell complex: the maximum
        dimension of its cells.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 2)]) ])
            sage: pc.dimension()
            2
            sage: empty_pc = PolyhedralComplex([])
            sage: empty_pc.dimension()
            -1
        """
        return self._dim

    def ambient_dimension(self):
        """
        The ambient dimension of this cell complex: the ambient
        dimension of each of its cells.

        EXAMPLES::

            sage: pc = PolyhedralComplex([Polyhedron(vertices=[(1, 2, 3)])])
            sage: pc.ambient_dimension()
            3
            sage: empty_pc = PolyhedralComplex([])
            sage: empty_pc.ambient_dimension()
            -1
            sage: pc0 = PolyhedralComplex(ambient_dim=2)
            sage: pc0.ambient_dimension()
            2
        """
        return self._ambient_dim

    def plot(self, **kwds):
        """
        Return a plot of the polyhedral complex, if it is of dim at most 3.

        INPUT:

        - ``explosion_factor`` -- (default: 0) if positive, separate the cells of
          the complex by extra space. In this case, the following keyword arguments
          can be passed to :func:`exploded_plot`:

          - ``center`` -- (default: ``None``, denoting the origin) the center of explosion
          - ``sticky_vertices`` -- (default: ``False``) boolean or dict;
            whether to draw line segments between shared vertices of the given polyhedra.
            A dict gives options for :func:`sage.plot.line`.
          - ``sticky_center`` -- (default: ``True``) boolean or dict. When ``center`` is
            a vertex of some of the polyhedra, whether to draw line segments connecting the
            ``center`` to the shifted copies of these vertices.
            A dict gives options for :func:`sage.plot.line`.

        - ``color`` -- (default: ``None``) if ``'rainbow'``, assign a different color
          to every maximal cell; otherwise, passed on to
          :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.plot`.

        - other keyword arguments are passed on to
          :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.plot`.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(0, 0), (0, 2), (-1, 1)])
            sage: pc1 = PolyhedralComplex([p1, p2, p3, -p1, -p2, -p3])
            sage: bb = dict(xmin=-2, xmax=2, ymin=-3, ymax=3, axes=False)
            sage: g0 = pc1.plot(color='rainbow', **bb)                                  # needs sage.plot
            sage: g1 = pc1.plot(explosion_factor=0.5, **bb)                             # needs sage.plot
            sage: g2 = pc1.plot(explosion_factor=1, color='rainbow', alpha=0.5, **bb)   # needs sage.plot
            sage: graphics_array([g0, g1, g2]).show(axes=False)                        # not tested

            sage: pc2 = PolyhedralComplex([polytopes.hypercube(3)])
            sage: pc3 = pc2.subdivide(new_vertices=[(0, 0, 0)])
            sage: g3 = pc3.plot(explosion_factor=1, color='rainbow',                    # needs sage.plot
            ....:               alpha=0.5, axes=False, online=True)
            sage: pc4 = pc2.subdivide(make_simplicial=True)
            sage: g4 = pc4.plot(explosion_factor=1, center=(1, -1, 1), fill='blue',     # needs sage.plot
            ....:              wireframe='white', point={'color':'red', 'size':10},
            ....:              alpha=0.6, online=True)
            sage: pc5 = PolyhedralComplex([
            ....:         Polyhedron(rays=[[1,0,0], [0,1,0], [0,0,-1]]),
            ....:         Polyhedron(rays=[[1,0,0], [0,-1,0], [0,0,-1]]),
            ....:         Polyhedron(rays=[[1,0,0], [0,-1,0], [0,0,1]]),
            ....:         Polyhedron(rays=[[-1,0,0], [0,-1,0], [0,0,-1]]),
            ....:         Polyhedron(rays=[[-1,0,0], [0,-1,0], [0,0,1]]),
            ....:         Polyhedron(rays=[[-1,0,0], [0,1,0], [0,0,-1]]),
            ....:         Polyhedron(rays=[[-1,0,0], [0,1,0], [0,0,1]])])
            sage: g5 = pc5.plot(explosion_factor=0.3, color='rainbow', alpha=0.8,       # needs sage.plot
            ....:               point={'size': 20}, axes=False, online=True)
        """
        if self.dimension() > 3:
            raise ValueError("cannot plot in high dimension")
        if kwds.get('explosion_factor', 0):
            return exploded_plot(self.maximal_cell_iterator(), **kwds)

        from sage.plot.colors import rainbow
        from sage.plot.graphics import Graphics

        color = kwds.get('color')
        polyhedra = self.maximal_cell_iterator()
        if color == 'rainbow':
            polyhedra = list(polyhedra)
            cell_colors_dict = dict(zip(polyhedra,
                                        rainbow(len(polyhedra))))
        g = Graphics()
        for cell in polyhedra:
            options = copy(kwds)
            if color == 'rainbow':
                options['color'] = cell_colors_dict[cell]
            g += cell.plot(**options)
        return g

    def is_pure(self):
        """
        Test if this polyhedral complex is pure.

        A polyhedral complex is pure if and only if all of its maximal cells
        have the same dimension.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: pc.is_pure()
            True

        Wrong answer due to ``maximality_check=False``::

            sage: pc_invalid = PolyhedralComplex([p1, p2, p3],
            ....:              maximality_check=False)
            sage: pc_invalid.is_pure()
            False
        """
        return len(self._maximal_cells) == 1

    def is_full_dimensional(self):
        """
        Return whether this polyhedral complex is full-dimensional:
        its dimension is equal to its ambient dimension.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: pc.is_full_dimensional()
            True
            sage: PolyhedralComplex([p3]).is_full_dimensional()
            False
        """
        return self._dim == self._ambient_dim

    def __hash__(self):
        """
        Compute the hash value of ``self`` using its ``maximal_cells_sorted``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])
            sage: pc1 = PolyhedralComplex([p1, p2], is_mutable=False)
            sage: hash(pc1) == hash(pc1)
            True
            sage: pc2 = PolyhedralComplex([p2, p1], is_mutable=False)
            sage: hash(pc1) == hash(pc2)
            True
            sage: pc3 = PolyhedralComplex([p1, p2])
            sage: hash(pc3)
            Traceback (most recent call last):
            ...
            ValueError: this polyhedral complex must be immutable; call set_immutable()
        """
        if not self._is_immutable:
            raise ValueError("this polyhedral complex must be immutable; " +
                             "call set_immutable()")
        return hash(tuple(self.maximal_cells_sorted()))

    def __eq__(self, right):
        """
        Two polyhedral complexes are equal iff their maximal cells are equal.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])
            sage: pc1 = PolyhedralComplex([p1, p2])
            sage: pc1 == pc1
            True
            sage: pc2 = PolyhedralComplex([p2, p1])
            sage: pc1 == pc2
            True
        """
        return isinstance(right, PolyhedralComplex) and (
               self.maximal_cells_sorted() == right.maximal_cells_sorted())

    def __ne__(self, right):
        """
        Return ``True`` if ``self`` and ``right`` are not equal.

        EXAMPLES::

            sage: pc1 = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)])])
            sage: pc2 = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])])
            sage: pc1 != pc2
            True
        """
        return not self.__eq__(right)

    def __copy__(self):
        """
        Return a mutable copy of ``self``.

        EXAMPLES::

            sage: pc1 = PolyhedralComplex([Polyhedron(vertices=[(0, 0)])])
            sage: pc2 = copy(pc1)
            sage: pc1 == pc2
            True
        """
        return PolyhedralComplex(self._maximal_cells, maximality_check=False,
                                 backend=self._backend)

    def _an_element_(self):
        """
        Return a (maximal) cell of this complex.

        EXAMPLES::

            sage: PolyhedralComplex()._an_element_()
            Traceback (most recent call last):
            ...
            EmptySetError: the complex is empty
            sage: pc = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])])
            sage: element = pc._an_element_().vertices_list()
            sage: element   # random output (one of the two maximal cells)
            [[0, 0], [0, 1/2], [1, 2]]
            sage: element in ([[0, 0], [0, 1/2], [1, 2]], [[0, 0], [1/3, 1/3], [1, 2]])
            True
        """
        try:
            return next(self.maximal_cell_iterator(increasing=False))
        except StopIteration:
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError("the complex is empty")

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is a polyhedron which is contained in this complex.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])
            sage: pc = PolyhedralComplex([p1, p2])
            sage: (p1 in pc) and (p2 in pc)
            True
            sage: Polyhedron(vertices=[(1, 2), (0, 0)]) in pc
            True
            sage: Polyhedron(vertices=[(1, 1), (0, 0)]) in pc
            False
            sage: Polyhedron(vertices=[(0, 0)]) in pc
            True
            sage: (0, 0) in pc  # not a polyhedron
            False
        """
        if not isinstance(x, sage.geometry.abc.Polyhedron):
            return False
        dim = x.dimension()
        return dim in self.cells() and x in self.cells()[dim]

    def __call__(self, x):
        """
        If ``x`` is a polyhedron in this complex, return it.
        Otherwise, raise a :exc:`ValueError`.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])])
            sage: pc(Polyhedron(vertices=[(1, 2), (0, 0)]))
            A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: pc(Polyhedron(vertices=[(1, 1)]))
            Traceback (most recent call last):
            ...
            ValueError: the polyhedron is not in this complex
        """
        if x not in self:
            raise ValueError('the polyhedron is not in this complex')
        return x

    def face_poset(self):
        r"""
        The face poset of this polyhedral complex, the poset of
        nonempty cells, ordered by inclusion.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])])
            sage: poset = pc.face_poset()
            sage: poset
            Finite poset containing 11 elements
            sage: d = {i: i.vertices_matrix() for i in poset}
            sage: poset.plot(element_labels=d)                                          # needs sage.plot
            Graphics object consisting of 28 graphics primitives

        For a nonbounded polyhedral complex::

            sage: pc = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)]),
            ....:         Polyhedron(vertices=[(-1/2, -1/2)], lines=[(1, -1)]),
            ....:         Polyhedron(rays=[(1, 0)])])
            sage: poset = pc.face_poset()
            sage: poset
            Finite poset containing 13 elements
            sage: d = {i:''.join([str(v)+'\n'
            ....:      for v in i.Vrepresentation()]) for i in poset}
            sage: poset.show(element_labels=d, figsize=15)        # not tested
            sage: pc = PolyhedralComplex([
            ....: Polyhedron(rays=[(1,0),(0,1)]),
            ....: Polyhedron(rays=[(-1,0),(0,1)]),
            ....: Polyhedron(rays=[(-1,0),(0,-1)]),
            ....: Polyhedron(rays=[(1,0),(0,-1)])])
            sage: pc.face_poset()
            Finite poset containing 9 elements
        """
        if self._face_poset is None:
            self.cells()    # poset is obtained and cached in cells()
        return self._face_poset

    def is_subcomplex(self, other):
        r"""
        Return whether ``self`` is a subcomplex of ``other``.

        INPUT:

        - ``other`` -- a polyhedral complex

        Each maximal cell of ``self`` must be a cell of ``other``
        for this to be ``True``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])
            sage: p3 = Polyhedron(vertices=[(0, 0), (1, 0)])
            sage: pc = PolyhedralComplex([p1, Polyhedron(vertices=[(1, 0)])])
            sage: pc.is_subcomplex(PolyhedralComplex([p1, p2, p3]))
            True
            sage: pc.is_subcomplex(PolyhedralComplex([p1, p2]))
            False
        """
        other_cells = other.cells()
        for (d, stratum) in self.maximal_cells().items():
            if not stratum.issubset(other_cells.get(d, set([]))):
                return False
        return True

    def is_compact(self):
        """
        Test for boundedness of the polyhedral complex.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])
            sage: p2 = Polyhedron(rays=[(1, 0)])
            sage: PolyhedralComplex([p1]).is_compact()
            True
            sage: PolyhedralComplex([p1, p2]).is_compact()
            False
        """
        return all(p.is_compact() for p in self.maximal_cell_iterator())

    def graph(self):
        """
        Return the 1-skeleton of this polyhedral complex, as a graph.

        The vertices of the graph are of type ``vector``. This raises
        a :exc:`NotImplementedError` if the polyhedral complex is unbounded.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: g = pc.graph(); g
            Graph on 4 vertices
            sage: g.vertices(sort=True)
            [(0, 0), (0, 2), (1, 1), (1, 2)]
            sage: g.edges(sort=True, labels=False)
            [((0, 0), (0, 2)), ((0, 0), (1, 1)), ((0, 0), (1, 2)), ((0, 2), (1, 2)), ((1, 1), (1, 2))]
            sage: PolyhedralComplex([Polyhedron(rays=[(1,1)])]).graph()
            Traceback (most recent call last):
            ...
            NotImplementedError: the polyhedral complex is unbounded

        Wrong answer due to ``maximality_check=False``::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: PolyhedralComplex([p1, p2]).is_pure()
            True
            sage: PolyhedralComplex([p2, p3], maximality_check=True).is_pure()
            True
            sage: PolyhedralComplex([p2, p3], maximality_check=False).is_pure()
            False
        """
        if not self.is_compact():
            raise NotImplementedError("the polyhedral complex is unbounded")
        edges = self.n_cells(1)
        d = {}
        for e in edges:
            v, max_e = sorted(e.vertices_matrix().columns())
            if v in d:
                d[v].append(max_e)
            else:
                d[v] = [max_e]
        for v in self.n_maximal_cells(0):
            d[v] = []
        return Graph(d)

    def is_connected(self):
        """
        Return whether ``self`` is connected.

        EXAMPLES::

            sage: pc1 = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: pc1.is_connected()
            True
            sage: pc2 = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(0, 2)])])
            sage: pc2.is_connected()
            False
            sage: pc3 = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)]),
            ....:         Polyhedron(vertices=[(-1/2, -1/2)], lines=[(1, -1)]),
            ....:         Polyhedron(rays=[(1, 0)])])
            sage: pc3.is_connected()
            False
            sage: pc4 = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)]),
            ....:         Polyhedron(rays=[(1, 0)])])
            sage: pc4.is_connected()
            True
        """
        if self.is_compact():
            return self.graph().is_connected()    # faster than using poset?
        else:
            return self.face_poset().is_connected()

    def connected_component(self, cell=None):
        """
        Return the connected component of this polyhedral complex
        containing a given cell.

        INPUT:

        - ``cell`` -- (default: ``self.an_element()``) a cell of ``self``

        OUTPUT:

        The connected component containing ``cell``. If the polyhedral complex
        is empty or if it does not contain the given cell, raise an error.

        EXAMPLES::

            sage: t1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: t2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: v1 = Polyhedron(vertices=[(1, 1)])
            sage: v2 = Polyhedron(vertices=[(0, 2)])
            sage: v3 = Polyhedron(vertices=[(-1, 0)])
            sage: o =  Polyhedron(vertices=[(0, 0)])
            sage: r = Polyhedron(rays=[(1, 0)])
            sage: l = Polyhedron(vertices=[(-1, 0)], lines=[(1, -1)])
            sage: pc1 = PolyhedralComplex([t1, t2])
            sage: pc1.connected_component() == pc1
            True
            sage: pc1.connected_component(v1) == pc1
            True
            sage: pc2 = PolyhedralComplex([t1, v2])
            sage: pc2.connected_component(t1) == PolyhedralComplex([t1])
            True
            sage: pc2.connected_component(o) == PolyhedralComplex([t1])
            True
            sage: pc2.connected_component(v3)
            Traceback (most recent call last):
            ...
            ValueError: the polyhedral complex does not contain the given cell
            sage: pc2.connected_component(r)
            Traceback (most recent call last):
            ...
            ValueError: the polyhedral complex does not contain the given cell
            sage: pc3 = PolyhedralComplex([t1, t2, r])
            sage: pc3.connected_component(v2) == pc3
            True
            sage: pc4 = PolyhedralComplex([t1, t2, r, l])
            sage: pc4.connected_component(o) == pc3
            True
            sage: pc4.connected_component(v3)
            Traceback (most recent call last):
            ...
            ValueError: the polyhedral complex does not contain the given cell
            sage: pc5 = PolyhedralComplex([t1, t2, r, l, v3])
            sage: pc5.connected_component(v3) == PolyhedralComplex([v3])
            True
            sage: PolyhedralComplex([]).connected_component()
            Traceback (most recent call last):
            ...
            ValueError: the empty polyhedral complex has no connected components
        """
        if self.dimension() == -1:
            raise ValueError(
                "the empty polyhedral complex has no connected components")
        if cell is None:
            cell = self._an_element_()
        if self.is_compact():    # use graph (faster than poset?)
            if not cell.is_compact():
                raise ValueError(
                    "the polyhedral complex does not contain the given cell")
            v = cell.vertices_matrix().columns()[0]
            g = self.graph()
            if v not in g:
                raise ValueError(
                    "the polyhedral complex does not contain the given cell")
            vertices = g.connected_component_containing_vertex(v, sort=False)
            facets = [f for f in self.maximal_cell_iterator()
                      if any(vf in f.vertices_matrix().columns()
                             for vf in vertices)]
        else:    # use face_poset
            g = self.face_poset().hasse_diagram()
            if cell not in g:
                raise ValueError(
                    "the polyhedral complex does not contain the given cell")
            faces = g.connected_component_containing_vertex(cell, sort=False)
            facets = [f for f in self.maximal_cell_iterator()
                      if f in faces]
        return PolyhedralComplex(facets, maximality_check=False,
                                 is_immutable=self._is_immutable,
                                 backend=self._backend)

    def connected_components(self):
        """
        Return the connected components of this polyhedral complex,
        as list of (sub-)PolyhedralComplexes.

        EXAMPLES::

            sage: t1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: t2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: v1 = Polyhedron(vertices=[(1, 1)])
            sage: v2 = Polyhedron(vertices=[(0, 2)])
            sage: v3 = Polyhedron(vertices=[(-1, 0)])
            sage: o =  Polyhedron(vertices=[(0, 0)])
            sage: r = Polyhedron(rays=[(1, 0)])
            sage: l = Polyhedron(vertices=[(-1, 0)], lines=[(1, -1)])
            sage: pc1 = PolyhedralComplex([t1, t2])
            sage: len(pc1.connected_components())
            1
            sage: pc2 = PolyhedralComplex([t1, v2])
            sage: len(pc2.connected_components())
            2
            sage: pc3 = PolyhedralComplex([t1, t2, r])
            sage: len(pc3.connected_components())
            1
            sage: pc4 = PolyhedralComplex([t1, t2, r, l])
            sage: len(pc4.connected_components())
            2
            sage: pc5 = PolyhedralComplex([t1, t2, r, l, v3])
            sage: len(pc5.connected_components())
            3
            sage: PolyhedralComplex([]).connected_components()
            Traceback (most recent call last):
            ...
            ValueError: the empty polyhedral complex has no connected components
        """
        if self.dimension() == -1:
            raise ValueError(
                "the empty polyhedral complex has no connected components")
        if self.is_compact():    # use graph (faster than poset)?
            g = self.graph()
            lists_of_vertices = g.connected_components(sort=False)
            lists_of_facets = [[f for f in self.maximal_cell_iterator()
                                if any(vf in f.vertices_matrix().columns()
                                       for vf in vertices)]
                               for vertices in lists_of_vertices]
        else:    # use face_poset
            g = self.face_poset().hasse_diagram()
            lists_of_faces = g.connected_components(sort=False)
            lists_of_facets = [
                [f for f in self.maximal_cell_iterator() if f in faces]
                for faces in lists_of_faces]
        results = [PolyhedralComplex(facets, maximality_check=False,
                                     is_immutable=self._is_immutable,
                                     backend=self._backend)
                   for facets in lists_of_facets]
        return results

    def n_skeleton(self, n):
        r"""
        The `n`-skeleton of this polyhedral complex.

        The `n`-skeleton of a polyhedral complex is obtained by discarding
        all of the cells in dimensions larger than `n`.

        INPUT:

        - ``n`` -- nonnegative integer; the dimension

        .. SEEALSO::

            :meth:`stratify`

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....:         Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: pc.n_skeleton(2)
            Polyhedral complex with 2 maximal cells
            sage: pc.n_skeleton(1)
            Polyhedral complex with 5 maximal cells
            sage: pc.n_skeleton(0)
            Polyhedral complex with 4 maximal cells
        """
        if n >= self.dimension():
            return copy(self)
        facets = [f for f in self.maximal_cell_iterator() if f.dimension() < n]
        facets.extend(self.n_cells(n))
        return PolyhedralComplex(facets, maximality_check=False,
                                 is_immutable=self._is_immutable,
                                 backend=self._backend)

    def stratify(self, n):
        r"""
        Return the pure sub-polyhedral complex which is constructed from the
        `n`-dimensional maximal cells of this polyhedral complex.

        .. SEEALSO::

            :meth:`n_skeleton`

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: pc.stratify(2) == pc
            True
            sage: pc.stratify(1)
            Polyhedral complex with 0 maximal cells

        Wrong answer due to ``maximality_check=False``::

            sage: pc_invalid = PolyhedralComplex([p1, p2, p3],
            ....:                                maximality_check=False)
            sage: pc_invalid.stratify(1)
            Polyhedral complex with 1 maximal cell
        """
        n_faces = self.n_maximal_cells(n)
        return PolyhedralComplex(n_faces, maximality_check=False,
                                 is_immutable=self._is_immutable,
                                 backend=self._backend)

    def boundary_subcomplex(self):
        """
        Return the sub-polyhedral complex that is the boundary of ``self``.

        A point `P` is on the boundary of a set `S` if `P` is in the
        closure of `S` but not in the interior of `S`.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: bd = PolyhedralComplex([p1, p2]).boundary_subcomplex()
            sage: len(bd.n_maximal_cells(2))
            0
            sage: len(bd.n_maximal_cells(1))
            4
            sage: pt = PolyhedralComplex([p3])
            sage: pt.boundary_subcomplex() == pt
            True

        Test on polyhedral complex which is not pure::

            sage: pc_non_pure = PolyhedralComplex([p1, p3])
            sage: pc_non_pure.boundary_subcomplex() == pc_non_pure.n_skeleton(1)
            True

        Test with ``maximality_check == False``::

            sage: pc_invalid = PolyhedralComplex([p2, p3],
            ....:                                maximality_check=False)
            sage: pc_invalid.boundary_subcomplex() == pc_invalid.n_skeleton(1)
            True

        Test unbounded cases::

            sage: pc1 = PolyhedralComplex([
            ....:         Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,0], [0,1]])])
            sage: pc1.boundary_subcomplex() == pc1.n_skeleton(1)
            True
            sage: pc1b = PolyhedralComplex([Polyhedron(
            ....:         vertices=[[1,0,0], [0,1,0]], rays=[[1,0,0],[0,1,0]])])
            sage: pc1b.boundary_subcomplex() == pc1b
            True
            sage: pc2 = PolyhedralComplex([
            ....:         Polyhedron(vertices=[[-1,0], [1,0]], lines=[[0,1]])])
            sage: pc2.boundary_subcomplex() == pc2.n_skeleton(1)
            True
            sage: pc3 = PolyhedralComplex([
            ....:         Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,0], [0,1]]),
            ....:         Polyhedron(vertices=[[1,0], [0,-1]], rays=[[1,0], [0,-1]])])
            sage: pc3.boundary_subcomplex() == pc3.n_skeleton(1)
            False
        """
        if self.is_full_dimensional():
            return PolyhedralComplex(self.relative_boundary_cells(),
                                     is_immutable=self._is_immutable,
                                     backend=self._backend)
        else:
            ans = copy(self)
            if self._is_immutable:
                ans.set_immutable()
            return ans

    def relative_boundary_cells(self):
        r"""
        Return the maximal cells of the relative-boundary sub-complex.

        A point `P` is in the relative boundary of a set `S` if `P` is in the
        closure of `S` but not in the relative interior of `S`.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: p4 = Polyhedron(vertices=[(2, 2)])
            sage: pc = PolyhedralComplex([p1, p2])
            sage: rbd_cells = pc.relative_boundary_cells()
            sage: len(rbd_cells)
            4
            sage: all(p.dimension() == 1 for p in rbd_cells)
            True
            sage: pc_lower_dim = PolyhedralComplex([p3])
            sage: sorted([p.vertices() for p in pc_lower_dim.relative_boundary_cells()])
            [(A vertex at (0, 2),), (A vertex at (1, 2),)]

        Test on polyhedral complex which is not pure::

            sage: pc_non_pure = PolyhedralComplex([p1, p3, p4])
            sage: (set(pc_non_pure.relative_boundary_cells())
            ....:  == set([f.as_polyhedron() for f in p1.faces(1)] + [p3, p4]))
            True

        Test with ``maximality_check == False``::

            sage: pc_invalid = PolyhedralComplex([p2, p3],
            ....:                                maximality_check=False)
            sage: (set(pc_invalid.relative_boundary_cells())
            ....:  == set([f.as_polyhedron() for f in p2.faces(1)]))
            True

        Test unbounded case::

            sage: pc3 = PolyhedralComplex([
            ....:         Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,0], [0,1]]),
            ....:         Polyhedron(vertices=[[1,0], [0,-1]], rays=[[1,0], [0,-1]])])
            sage: len(pc3.relative_boundary_cells())
            4
        """
        d = self.dimension()
        poset = self.face_poset()
        faces = self.n_cells(d - 1)
        ans = [face for face in faces if len(poset.upper_covers(face)) == 1]
        if not self.is_pure():
            ans += [p for p in poset.maximal_elements() if p.dimension() < d]
        return ans

    def is_convex(self):
        r"""
        Return whether the set of points in ``self`` is a convex set.

        When ``self`` is convex, the union of its cells is a Polyhedron.

        .. SEEALSO::

            :meth:`union_as_polyhedron`

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(0, 0), (1, 1), (2, 0)])
            sage: p4 = Polyhedron(vertices=[(2, 2)])
            sage: PolyhedralComplex([p1, p2]).is_convex()
            True
            sage: PolyhedralComplex([p1, p3]).is_convex()
            False
            sage: PolyhedralComplex([p1, p4]).is_convex()
            False

        Test unbounded cases::

            sage: pc1 = PolyhedralComplex([
            ....:         Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,0], [0,1]])])
            sage: pc1.is_convex()
            True
            sage: pc2 = PolyhedralComplex([
            ....:         Polyhedron(vertices=[[-1,0], [1,0]], lines=[[0,1]])])
            sage: pc2.is_convex()
            True
            sage: pc3 = PolyhedralComplex([
            ....:         Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,0], [0,1]]),
            ....:         Polyhedron(vertices=[[1,0], [0,-1]], rays=[[1,0], [0,-1]])])
            sage: pc3.is_convex()
            False
            sage: pc4 = PolyhedralComplex([Polyhedron(rays=[[1,0], [-1,1]]),
            ....:                          Polyhedron(rays=[[1,0], [-1,-1]])])
            sage: pc4.is_convex()
            False

        The whole 3d space minus the first orthant is not convex::

            sage: pc5 = PolyhedralComplex([
            ....:         Polyhedron(rays=[[1,0,0], [0,1,0], [0,0,-1]]),
            ....:         Polyhedron(rays=[[1,0,0], [0,-1,0], [0,0,-1]]),
            ....:         Polyhedron(rays=[[1,0,0], [0,-1,0], [0,0,1]]),
            ....:         Polyhedron(rays=[[-1,0,0], [0,-1,0], [0,0,-1]]),
            ....:         Polyhedron(rays=[[-1,0,0], [0,-1,0], [0,0,1]]),
            ....:         Polyhedron(rays=[[-1,0,0], [0,1,0], [0,0,-1]]),
            ....:         Polyhedron(rays=[[-1,0,0], [0,1,0], [0,0,1]])])
            sage: pc5.is_convex()
            False

        Test some non-full-dimensional examples::

            sage: l = PolyhedralComplex([Polyhedron(vertices=[(1, 2), (0, 2)])])
            sage: l.is_convex()
            True
            sage: pc1b = PolyhedralComplex([Polyhedron(
            ....:         vertices=[[1,0,0], [0,1,0]], rays=[[1,0,0],[0,1,0]])])
            sage: pc1b.is_convex()
            True
            sage: pc4b = PolyhedralComplex([
            ....:         Polyhedron(rays=[[1,0,0], [-1,1,0]]),
            ....:         Polyhedron(rays=[[1,0,0], [-1,-1,0]])])
            sage: pc4b.is_convex()
            False
        """
        if self._is_convex is not None:
            return self._is_convex
        if not self.is_pure():
            self._is_convex = False
            return False
        d = self.dimension()
        if not self.is_full_dimensional():
            # if max cells must lie in different subspaces, can't be convex.
            from sage.modules.free_module import span
            f = self.n_maximal_cells(d)[0]
            affine_space = span(f.equations_list(), f.base_ring())
            for f in self.n_maximal_cells(d)[1::]:
                if span(f.equations_list(), f.base_ring()) != affine_space:
                    self._is_convex = False
                    return False
        # orient the (relative) boundary halfspaces toward a strict convex
        # combination of the vertices. Then check if all vertices are contained
        # After making sure that the affine hulls of the cells are the same,
        # it does not matter that is not full dimensional.
        boundaries = self.relative_boundary_cells()
        vertices = set([])
        rays = set([])
        lines = set([])
        for cell in boundaries:
            # it suffices to consider only vertices on the boundaries
            # Note that a line (as polyhedron) has vertex too
            for v in cell.vertices_list():
                vv = vector(v)
                vv.set_immutable()
                vertices.add(vv)
        for cell in self.n_maximal_cells(d):
            for r in cell.rays_list():
                rr = vector(r)
                rr.set_immutable()
                rays.add(rr)
            for li in cell.lines_list():
                ll = vector(li)
                ll.set_immutable()
                lines.add(ll)
        center = sum(vertices) / len(vertices)
        for cell in boundaries:
            for equation in cell.equations_list():
                coeff = vector(equation[1::])
                const = equation[0]
                if const + coeff * center == 0:
                    sign = 0
                elif const + coeff * center > 0:
                    sign = 1
                    for v in vertices:
                        if const + coeff * v < 0:
                            self._is_convex = False
                            return False
                elif const + coeff * center < 0:
                    sign = -1
                    for v in vertices:
                        if const + coeff * v > 0:
                            self._is_convex = False
                            return False
                for r in rays:
                    if sign == 0:
                        sign = coeff * r
                    else:
                        if sign * (coeff * r) < 0:
                            self._is_convex = False
                            return False
                # lines are in the affine space of each boundary cell already
        self._is_convex = True
        self._polyhedron = Polyhedron(vertices=vertices, rays=rays, lines=lines,
                                      backend=self._backend)
        return True

    def union_as_polyhedron(self):
        """
        Return ``self`` as a :class:`Polyhedron` if ``self`` is convex.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(0, 0), (1, 1), (2, 0)])
            sage: P = PolyhedralComplex([p1, p2]).union_as_polyhedron()
            sage: P.vertices_list()
            [[0, 0], [0, 2], [1, 1], [1, 2]]
            sage: PolyhedralComplex([p1, p3]).union_as_polyhedron()
            Traceback (most recent call last):
            ...
            ValueError: the polyhedral complex is not convex
        """
        if not self.is_convex():
            raise ValueError("the polyhedral complex is not convex")
        return self._polyhedron

    def product(self, right):
        """
        The (Cartesian) product of this polyhedral complex with another one.

        INPUT:

        - ``right`` -- the other polyhedral complex (the right-hand factor)

        OUTPUT: the product ``self x right``

        EXAMPLES::

            sage: pc = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])])
            sage: pc_square = pc.product(pc)
            sage: pc_square
            Polyhedral complex with 1 maximal cell
            sage: next(pc_square.maximal_cell_iterator()).vertices()
            (A vertex at (0, 0),
             A vertex at (0, 1),
             A vertex at (1, 0),
             A vertex at (1, 1))
        """
        maximal_cells = [f.product(g) for f in self.maximal_cell_iterator()
                         for g in right.maximal_cell_iterator()]
        return PolyhedralComplex(maximal_cells, maximality_check=False,
                                 is_immutable=(self._is_immutable and
                                               right._is_immutable),
                                 backend=self._backend)

    def disjoint_union(self, right):
        """
        The disjoint union of this polyhedral complex with another one.

        INPUT:

        - ``right`` -- the other polyhedral complex (the right-hand factor)

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(-1, 0), (0, 0), (0, 1)])
            sage: p2 = Polyhedron(vertices=[(0, -1), (0, 0), (1, 0)])
            sage: p3 = Polyhedron(vertices=[(0, -1), (1, -1), (1, 0)])
            sage: pc = PolyhedralComplex([p1]).disjoint_union(PolyhedralComplex([p3]))
            sage: set(pc.maximal_cell_iterator()) == set([p1, p3])
            True
            sage: pc.disjoint_union(PolyhedralComplex([p2]))
            Traceback (most recent call last):
            ...
            ValueError: the two complexes are not disjoint
        """
        maximal_cells_self = list(self.maximal_cell_iterator())
        maximal_cells_right = list(right.maximal_cell_iterator())
        for cell in maximal_cells_self:
            for cell_right in maximal_cells_right:
                if not cell.intersection(cell_right).is_empty():
                    raise ValueError("the two complexes are not disjoint")
        return PolyhedralComplex(maximal_cells_self + maximal_cells_right,
                                 maximality_check=False,
                                 face_to_face_check=False,
                                 is_immutable=(self._is_immutable and
                                               right._is_immutable),
                                 backend=self._backend)

    def union(self, right):
        """
        The union of this polyhedral complex with another one.

        INPUT:

        - ``right`` -- the other polyhedral complex (the right-hand factor)

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(-1, 0), (0, 0), (0, 1)])
            sage: p2 = Polyhedron(vertices=[(0, -1), (0, 0), (1, 0)])
            sage: p3 = Polyhedron(vertices=[(0, -1), (1, -1), (1, 0)])
            sage: pc = PolyhedralComplex([p1]).union(PolyhedralComplex([p3]))
            sage: set(pc.maximal_cell_iterator()) == set([p1, p3])
            True
            sage: pc.union(PolyhedralComplex([p2]))
            Polyhedral complex with 3 maximal cells
            sage: p4 = Polyhedron(vertices=[(0, -1), (0, 0), (1, 0), (1, -1)])
            sage: pc.union(PolyhedralComplex([p4]))
            Traceback (most recent call last):
            ...
            ValueError: the given cells are not face-to-face
        """
        maximal_cells = list(self.maximal_cell_iterator()) + list(
                        right.maximal_cell_iterator())
        return PolyhedralComplex(maximal_cells, maximality_check=True,
                                 face_to_face_check=True,
                                 is_immutable=(self._is_immutable and
                                               right._is_immutable),
                                 backend=self._backend)

    def join(self, right):
        """
        The join of this polyhedral complex with another one.

        INPUT:

        - ``right`` -- the other polyhedral complex (the right-hand factor)

        EXAMPLES::

            sage: pc = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])])
            sage: pc_join = pc.join(pc)
            sage: pc_join
            Polyhedral complex with 1 maximal cell
            sage: next(pc_join.maximal_cell_iterator()).vertices()
            (A vertex at (0, 0, 0),
             A vertex at (0, 0, 1),
             A vertex at (0, 1, 1),
             A vertex at (1, 0, 0))
        """
        maximal_cells = [f.join(g) for f in self.maximal_cell_iterator()
                         for g in right.maximal_cell_iterator()]
        return PolyhedralComplex(maximal_cells, maximality_check=False,
                                 is_immutable=(self._is_immutable and
                                               right._is_immutable),
                                 backend=self._backend)

    ############################################################
    # abstract methods not implemented in generic cell complex
    ############################################################

    def wedge(self, right):
        """
        The wedge (one-point union) of ``self`` with ``right``.

        .. TODO::

            Implement the wedge product of two polyhedral complexes.

        EXAMPLES::

            sage: pc = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])])
            sage: pc.wedge(pc)
            Traceback (most recent call last):
            ...
            NotImplementedError: wedge is not implemented for polyhedral complex
        """
        raise NotImplementedError("wedge is not implemented for "
                                  + "polyhedral complex")

    ############################################################
    # chain complexes, homology
    ############################################################
    def chain_complex(self, subcomplex=None, augmented=False,
                      verbose=False, check=True, dimensions=None,
                      base_ring=ZZ, cochain=False):
        """
        The chain complex associated to this polyhedral complex.

        .. TODO::

            Implement chain complexes of a polyhedral complex.

        EXAMPLES::

            sage: pc = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])])
            sage: pc.chain_complex()
            Traceback (most recent call last):
            ...
            NotImplementedError: chain_complex is not implemented for polyhedral complex
        """
        raise NotImplementedError("chain_complex is not implemented for "
                                  + "polyhedral complex")

    def alexander_whitney(self, cell, dim_left):
        """
        The decomposition of ``cell`` in this complex into left and right
        factors, suitable for computing cup products.

        .. TODO::

            Implement :meth:`alexander_whitney` of a polyhedral complex.

        EXAMPLES::

            sage: pc = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])])
            sage: pc.alexander_whitney(None, 1)
            Traceback (most recent call last):
            ...
            NotImplementedError: alexander_whitney is not implemented for polyhedral complex
        """
        raise NotImplementedError("alexander_whitney is not implemented for "
                                  + "polyhedral complex")

    ############################################################
    # end of chain complexes, homology
    ############################################################

    # this function overrides the standard one for GenericCellComplex,
    # this one counts the number of maximal cells, not all cells, to
    # avoid calling and computing self.cells()
    def _repr_(self):
        """
        Print representation.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: PolyhedralComplex([p1, p2, p3])
            Polyhedral complex with 2 maximal cells

        Wrong answer due to ``maximality_check=False``::

            sage: PolyhedralComplex([p1, p2, p3], maximality_check=False)
            Polyhedral complex with 3 maximal cells
        """
        num = len(list(self.maximal_cell_iterator()))
        if num == 1:
            return "Polyhedral complex with %s maximal cell" % num
        else:
            return "Polyhedral complex with %s maximal cells" % num

    def set_immutable(self):
        """
        Make this polyhedral complex immutable.

        EXAMPLES::

            sage: pc = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])])
            sage: pc.is_mutable()
            True
            sage: pc.set_immutable()
            sage: pc.is_mutable()
            False
        """
        self._is_immutable = True

    def is_mutable(self):
        """
        Return whether ``self`` is mutable.

        EXAMPLES::

            sage: pc1 = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])])
            sage: pc1.is_mutable()
            True
            sage: pc2 = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])],
            ....:                        is_mutable=False)
            sage: pc2.is_mutable()
            False
            sage: pc1 == pc2
            True
            sage: pc3 = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])],
            ....:                        is_immutable=True)
            sage: pc3.is_mutable()
            False
            sage: pc2 == pc3
            True
        """
        return not self._is_immutable

    def is_immutable(self):
        """
        Return whether ``self`` is immutable.

        EXAMPLES::

            sage: pc1 = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])])
            sage: pc1.is_immutable()
            False
            sage: pc2 = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])],
            ....:                        is_mutable=False)
            sage: pc2.is_immutable()
            True
            sage: pc3 = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])],
            ....:                        is_immutable=True)
            sage: pc3.is_immutable()
            True
        """
        return self._is_immutable

    def add_cell(self, cell):
        """
        Add a cell to this polyhedral complex.

        INPUT:

        - ``cell`` -- a polyhedron

        This **changes** the polyhedral complex, by adding a new cell and all
        of its subfaces.

        EXAMPLES:

        Set up an empty complex in the intended ambient space, then add a cell::

            sage: pc = PolyhedralComplex(ambient_dim=2)
            sage: pc.add_cell(Polyhedron(vertices=[(1, 2), (0, 2)]))
            sage: pc
            Polyhedral complex with 1 maximal cell

        If you add a cell which is already present, there is no effect::

            sage: pc.add_cell(Polyhedron(vertices=[(1, 2)]))
            sage: pc
            Polyhedral complex with 1 maximal cell
            sage: pc.dimension()
            1

        Add a cell and check that dimension is correctly updated::

            sage: pc.add_cell(Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)]))
            sage: pc.dimension()
            2
            sage: pc.maximal_cells()
            {2: {A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices}}
            sage: pc.is_convex()
            True

        Add another cell and check that the properties are correctly updated::

            sage: pc.add_cell(Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]))
            sage: pc
            Polyhedral complex with 2 maximal cells
            sage: len(pc._cells[1])
            5
            sage: pc._face_poset
            Finite poset containing 11 elements
            sage: pc._is_convex
            True
            sage: pc._polyhedron.vertices_list()
            [[0, 0], [0, 2], [1, 1], [1, 2]]

        Add a ray which makes the complex non convex::

            sage: pc.add_cell(Polyhedron(rays=[(1, 0)]))
            sage: pc
            Polyhedral complex with 3 maximal cells
            sage: len(pc._cells[1])
            6
            sage: (pc._is_convex is False) and (pc._polyhedron is None)
            True

        TESTS::

            sage: pc.add_cell(Polyhedron(vertices=[[0]]))
            Traceback (most recent call last):
            ...
            ValueError: the given cell is not a polyhedron in the same ambient space
            sage: pc.add_cell(Polyhedron(vertices=[(1, 1), (0, 0), (2, 0)]))
            Traceback (most recent call last):
            ...
            ValueError: the cell is not face-to-face with complex
            sage: pc.set_immutable()
            sage: pc.add_cell(Polyhedron(vertices=[(-1, -1)]))
            Traceback (most recent call last):
            ...
            ValueError: this polyhedral complex is not mutable
        """
        if self._is_immutable:
            raise ValueError("this polyhedral complex is not mutable")
        if not isinstance(cell, sage.geometry.abc.Polyhedron) or cell.ambient_dim() != self._ambient_dim:
            raise ValueError("the given cell is not a polyhedron " +
                             "in the same ambient space")
        # if cell is already in self, do nothing.
        if self.is_cell(cell):
            return
        if self._backend:
            cell = cell.base_extend(cell.base_ring(), self._backend)
        # update cells and face poset
        cells = self.cells()
        covers = {p: self.face_poset().upper_covers(p)
                  for p in self.cell_iterator()}
        d = cell.dimension()
        d_cells = [cell]
        if d not in cells:
            cells[d] = set(d_cells)
        else:
            cells[d].add(cell)
        covers[cell] = []
        while d > 0:
            d = d - 1
            new_facets = []
            for c in d_cells:
                for facet in c.facets():
                    p = facet.as_polyhedron()
                    if d not in cells:
                        cells[d] = set([])
                    if p not in cells[d]:
                        cells[d].add(p)
                        covers[p] = [c]
                        new_facets.append(p)
                    else:
                        covers[p].append(c)
            d_cells = new_facets
        self._face_poset = poset = Poset(covers)
        self._cells = cells
        # check face-to-face between cell and previous maximal cells
        for p in self.maximal_cell_iterator():
            r = p.intersection(cell)
            if not (r.is_empty() or (r in poset) and
                    poset.is_gequal(p, r) and poset.is_gequal(cell, r)):
                raise ValueError("the cell is not face-to-face with complex")
        # update dim and maximal cells
        d = cell.dimension()
        if d > self._dim:
            self._dim = d
        maximal_cells = poset.maximal_elements()    # a list
        self._maximal_cells = cells_list_to_cells_dict(maximal_cells)
        # update convexity if self was known to be convex, reset otherwise.
        if self._is_convex:
            try:
                new_complex = PolyhedralComplex([self._polyhedron, cell],
                                                face_to_face_check=True)
            except ValueError:
                self._is_convex = False
                self._polyhedron = None
            else:
                self._is_convex = new_complex.is_convex()
                self._polyhedron = new_complex._polyhedron
        else:
            self._is_convex = None
            self._polyhedron = None
        # reset cached attribute
        self._maximal_cells_sorted = None    # needed for hash

    def remove_cell(self, cell, check=False):
        r"""
        Remove ``cell`` from ``self`` and all the cells that contain ``cell``
        as a subface.

        INPUT:

        - ``cell`` -- a cell of the polyhedral complex

        - ``check`` -- boolean (default: ``False``); if ``True``,
          raise an error if ``cell`` is not a cell of this complex

        This does not return anything; instead, it **changes** the
        polyhedral complex.

        EXAMPLES:

        If you add a cell which is already present, there is no effect::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: r = Polyhedron(rays=[(1, 0)])
            sage: pc = PolyhedralComplex([p1, p2, r])
            sage: pc.dimension()
            2
            sage: pc.remove_cell(Polyhedron(vertices=[(0, 0), (1, 2)]))
            sage: pc.dimension()
            1
            sage: pc
            Polyhedral complex with 5 maximal cells
            sage: pc.remove_cell(Polyhedron(vertices=[(1, 2)]))
            sage: pc.dimension()
            1
            sage: pc
            Polyhedral complex with 3 maximal cells
            sage: pc.remove_cell(Polyhedron(vertices=[(0, 0)]))
            sage: pc.dimension()
            0

        TESTS:

        Check that :exc:`ValueError` and empty complex are treated properly::

            sage: p = Polyhedron(vertices=[[1]])
            sage: pc = PolyhedralComplex([p])
            sage: pc.remove_cell(Polyhedron(vertices=[[0]]), check=True)
            Traceback (most recent call last):
            ...
            ValueError: trying to remove a cell which is not in the polyhedral complex
            sage: pc.remove_cell(Polyhedron(vertices=[(1, 1)]))
            Traceback (most recent call last):
            ...
            ValueError: the given cell is not a polyhedron in the same ambient space
            sage: pc.remove_cell(p)
            sage: pc.dimension()
            -1
            sage: pc = PolyhedralComplex([Polyhedron(vertices=[[0]])], is_mutable=False)
            sage: pc.remove_cell(Polyhedron(vertices=[[0]]))
            Traceback (most recent call last):
            ...
            ValueError: this polyhedral complex is not mutable

        Check that this function is coherent with
        :meth:`~sage.topology.simplicial_complex.SimplicialComplex.remove_face`::

            sage: v1 = (1, 0, 0, 0); v2 = (0, 1, 0, 0); v3 = (0, 0, 1, 0); v4 = (0, 0, 0, 1)
            sage: Z = PolyhedralComplex([Polyhedron(vertices=[v1, v2, v3, v4])]); Z
            Polyhedral complex with 1 maximal cell
            sage: Z.remove_cell(Polyhedron(vertices=[v1, v2]))
            sage: Z
            Polyhedral complex with 2 maximal cells
            sage: [c.vertices_list() for c in Z.maximal_cells_sorted()]
            [[[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]],
             [[0, 0, 0, 1], [0, 0, 1, 0], [1, 0, 0, 0]]]

            sage: v0 = (0, 0, 0, 0)
            sage: S = PolyhedralComplex([Polyhedron(vertices=[v0, v1, v2]), Polyhedron(vertices=[v2, v3])])
            sage: S.maximal_cells()
            {1: {A 1-dimensional polyhedron in ZZ^4 defined as the convex hull of 2 vertices},
             2: {A 2-dimensional polyhedron in ZZ^4 defined as the convex hull of 3 vertices}}
            sage: S.remove_cell(Polyhedron(vertices=[v0, v1, v2]))
            sage: S
            Polyhedral complex with 4 maximal cells
            sage: [c.vertices_list() for c in S.maximal_cells_sorted()]
            [[[0, 0, 0, 0], [0, 1, 0, 0]],
             [[0, 0, 0, 0], [1, 0, 0, 0]],
             [[0, 0, 1, 0], [0, 1, 0, 0]],
             [[0, 1, 0, 0], [1, 0, 0, 0]]]

            sage: T = PolyhedralComplex([Polyhedron(vertices=[[1], [2]]), Polyhedron(vertices=[[1], [-3]])])
            sage: T.remove_cell(Polyhedron(vertices=[[-3], [1]]))
            sage: [c.vertices_list() for c in T.maximal_cells_sorted()]
            [[[1], [2]], [[-3]]]
            sage: [c.vertices_list() for c in T.cells_sorted()]
            [[[1], [2]], [[-3]], [[1]], [[2]]]
        """
        if self._is_immutable:
            raise ValueError("this polyhedral complex is not mutable")
        if not isinstance(cell, sage.geometry.abc.Polyhedron) or cell.ambient_dim() != self._ambient_dim:
            raise ValueError("the given cell is not a polyhedron " +
                             "in the same ambient space")
        # if cell is not in self, delete nothing.
        if not self.is_cell(cell):   # self.cells() is called
            if check:
                raise ValueError("trying to remove a cell which is not " +
                                 "in the polyhedral complex")
            return
        # update cells and face poset
        poset = self._face_poset
        deleting = poset.order_filter([cell])
        for c in deleting:
            d = c.dimension()
            self._cells[d].remove(c)
            if not self._cells[d]:
                del self._cells[d]
        covers = {p: [q for q in poset.upper_covers(p) if q not in deleting]
                  for p in self.cell_iterator()}
        self._face_poset = Poset(covers)
        # update dim and maximal cells
        maximal_cells = self._face_poset.maximal_elements()    # a list
        self._maximal_cells = cells_list_to_cells_dict(maximal_cells)
        if not maximal_cells:
            self._dim = -1
        else:
            self._dim = max(self._maximal_cells.keys())
        # reset cached attributes
        self._maximal_cells_sorted = None    # needed for hash
        self._is_convex = None
        self._polyhedron = None

    def is_simplicial_complex(self):
        """
        Test if this polyhedral complex is a simplicial complex.

        A polyhedral complex is **simplicial** if all of its (maximal) cells
        are simplices, i.e., every cell is a bounded polytope with `d+1`
        vertices, where `d` is the dimension of the polytope.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(0, 0), (1, 1), (1, 2)])
            sage: p2 = Polyhedron(rays=[(1, 0)])
            sage: PolyhedralComplex([p1]).is_simplicial_complex()
            True
            sage: PolyhedralComplex([p2]).is_simplicial_complex()
            False
        """
        return all(p.is_simplex() for p in self.maximal_cell_iterator())

    def is_polyhedral_fan(self):
        """
        Test if this polyhedral complex is a polyhedral fan.

        A polyhedral complex is a **fan** if all of its (maximal) cells
        are cones.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(0, 0), (1, 1), (1, 2)])
            sage: p2 = Polyhedron(rays=[(1, 0)])
            sage: PolyhedralComplex([p1]).is_polyhedral_fan()
            False
            sage: PolyhedralComplex([p2]).is_polyhedral_fan()
            True
            sage: halfplane = Polyhedron(rays=[(1, 0), (-1, 0), (0, 1)])
            sage: PolyhedralComplex([halfplane]).is_polyhedral_fan()
            True
        """
        return all((p.n_vertices() == 1) and (
                   vector(p.vertices_list()[0]) == p.ambient_space().zero())
                   for p in self.maximal_cell_iterator())

    def is_simplicial_fan(self):
        """
        Test if this polyhedral complex is a simplicial fan.

        A polyhedral complex is a **simplicial fan** if all of its (maximal)
        cells are simplical cones, i.e., every cell is a pointed cone (with
        vertex being the origin) generated by `d` linearly independent rays,
        where `d` is the dimension of the cone.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(0, 0), (1, 1), (1, 2)])
            sage: p2 = Polyhedron(rays=[(1, 0)])
            sage: PolyhedralComplex([p1]).is_simplicial_fan()
            False
            sage: PolyhedralComplex([p2]).is_simplicial_fan()
            True
            sage: halfplane = Polyhedron(rays=[(1, 0), (-1, 0), (0, 1)])
            sage: PolyhedralComplex([halfplane]).is_simplicial_fan()
            False
        """
        return self.is_polyhedral_fan() and all(
               (p.n_lines() == 0 and p.n_rays() == p.dimension())
               for p in self.maximal_cell_iterator())

    def subdivide(self, make_simplicial=False,
                  new_vertices=None, new_rays=None):
        """
        Construct a new polyhedral complex by iterative stellar subdivision of
        ``self`` for each new vertex/ray given.

        Currently, subdivision is only supported for bounded polyhedral complex
        or polyhedral fan.

        INPUT:

        - ``make_simplicial`` -- boolean (default: ``False``); if ``True``,
          the returned polyhedral complex is simplicial

        - ``new_vertices``, ``new_rays`` -- list (optional); new generators
          to be added during subdivision

        EXAMPLES::

            sage: square_vertices = [(1, 1, 1), (-1, 1, 1), (-1, -1, 1), (1, -1, 1)]
            sage: pc = PolyhedralComplex([
            ....:         Polyhedron(vertices=[(0, 0, 0)] + square_vertices),
            ....:         Polyhedron(vertices=[(0, 0, 2)] + square_vertices)])
            sage: pc.is_compact() and not pc.is_simplicial_complex()
            True
            sage: subdivided_pc = pc.subdivide(new_vertices=[(0, 0, 1)])
            sage: subdivided_pc
            Polyhedral complex with 8 maximal cells
            sage: subdivided_pc.is_simplicial_complex()
            True
            sage: simplicial_pc = pc.subdivide(make_simplicial=True)
            sage: simplicial_pc
            Polyhedral complex with 4 maximal cells
            sage: simplicial_pc.is_simplicial_complex()
            True

            sage: fan = PolyhedralComplex([Polyhedron(rays=square_vertices)])
            sage: fan.is_polyhedral_fan() and not fan.is_simplicial_fan()
            True
            sage: fan.subdivide(new_vertices=[(0, 0, 1)])
            Traceback (most recent call last):
            ...
            ValueError: new vertices cannot be used for subdivision
            sage: subdivided_fan = fan.subdivide(new_rays=[(0, 0, 1)])
            sage: subdivided_fan
            Polyhedral complex with 4 maximal cells
            sage: subdivided_fan.is_simplicial_fan()
            True
            sage: simplicial_fan = fan.subdivide(make_simplicial=True)
            sage: simplicial_fan
            Polyhedral complex with 2 maximal cells
            sage: simplicial_fan.is_simplicial_fan()
            True

            sage: halfspace = PolyhedralComplex([Polyhedron(rays=[(0, 0, 1)],
            ....:             lines=[(1, 0, 0), (0, 1, 0)])])
            sage: halfspace.is_simplicial_fan()
            False
            sage: subdiv_halfspace = halfspace.subdivide(make_simplicial=True)
            sage: subdiv_halfspace
            Polyhedral complex with 4 maximal cells
            sage: subdiv_halfspace.is_simplicial_fan()
            True
        """
        if self.is_compact():
            if new_rays:
                raise ValueError("rays/lines cannot be used for subdivision")
            # bounded version of `fan.subdivide`; not require rational.
            vertices = set([])
            if make_simplicial and not self.is_simplicial_complex():
                for p in self.maximal_cell_iterator():
                    for v in p.vertices_list():
                        vertices.add(tuple(v))
            if new_vertices:
                for v in new_vertices:
                    vertices.add(tuple(v))
            if not vertices:
                return self    # Nothing has to be done
            # bounded version of `fan._subdivide_stellar`; not require rational.
            cells = list(self.maximal_cell_iterator())
            for v in vertices:
                new = []
                for cell in cells:
                    if v in cell:
                        for cell_facet in cell.facets():
                            facet = cell_facet.as_polyhedron()
                            if v in facet:
                                continue
                            p = facet.convex_hull(Polyhedron(vertices=[v]))
                            new.append(p)
                    else:
                        new.append(cell)
                cells = new
            return PolyhedralComplex(cells, maximality_check=False,
                                     backend=self._backend)
        elif self.is_polyhedral_fan():
            if new_vertices and any(vi != 0 for v in new_vertices for vi in v):
                raise ValueError("new vertices cannot be used for subdivision")
            # mimic :meth:`~sage.geometry.fan <RationalPolyhedralFan>.subdivide`
            # but here we allow for non-pointed cones, and we subdivide them.
            rays_normalized = set([])
            self_rays = []
            cones = []
            for p in self.maximal_cell_iterator():
                prays = p.rays_list()
                for r in prays:
                    r_n = vector(r).normalized()
                    r_n.set_immutable()
                    if r_n not in rays_normalized:
                        rays_normalized.add(r_n)
                        self_rays.append(vector(r))
                plines = p.lines_list()
                if not plines:
                    cones.append(p)
                    continue
                # consider a line as two rays
                for pl in plines:
                    l_plus = vector(pl).normalized()
                    l_plus.set_immutable()
                    if l_plus not in rays_normalized:
                        rays_normalized.add(l_plus)
                        self_rays.append(vector(pl))
                    l_minus = (-vector(pl)).normalized()
                    l_minus.set_immutable()
                    if l_minus not in rays_normalized:
                        rays_normalized.add(l_minus)
                        self_rays.append(-vector(pl))
                # subdivide the non-pointed p into pointed cones
                # we rely on the canonical V-repr of Sage polyhedra.
                num_lines = len(plines)
                for neg_rays in powerset(range(num_lines)):
                    lines = [vector(plines[i]) if i not in neg_rays
                             else -vector(plines[i]) for i in range(num_lines)]
                    cones.append(Polyhedron(rays=(prays + lines),
                                            backend=self._backend))
            rays = []
            if new_rays:
                for r in new_rays:
                    if vector(r).is_zero():
                        raise ValueError("zero cannot be used for subdivision")
                    r_n = vector(r).normalized()
                    r_n.set_immutable()
                    if r_n not in rays_normalized:
                        rays_normalized.add(r_n)
                        rays.append(vector(r))
            if make_simplicial and not self.is_simplicial_fan():
                rays = self_rays + rays
            if not rays:
                return self    # Nothing has to be done
            # mimic :class:`RationalPolyhedralFan`._subdivide_stellar(rays)
            # start with self maximal cells (subdivided into pointed cones)
            for ray in rays:
                new = []
                for cone in cones:
                    if ray in cone:
                        for cone_facet in cone.facets():
                            facet = cone_facet.as_polyhedron()
                            if ray in facet:
                                continue
                            new_cone = facet.convex_hull(Polyhedron(rays=[ray]))
                            new.append(new_cone)
                    else:
                        new.append(cone)
                cones = new
            return PolyhedralComplex(cones, maximality_check=False,
                                     backend=self._backend)
        else:
            # TODO: ``self`` is unbounded, make it projectively simplicial.
            # (1) homogenize self of dim d to fan in space of dim d+1;
            # (2) call fan.subdivide(make_simplicial=True);
            # (3) take section back to the space of dim d.
            raise NotImplementedError('subdivision of a non-compact polyhedral ' +
                                      'complex that is not a fan is not supported')

############################################################
# Helper functions
############################################################


def cells_list_to_cells_dict(cells_list):
    r"""
    Helper function that returns the dictionary whose keys are the dimensions,
    and the value associated to an integer `d` is the set of `d`-dimensional
    polyhedra in the given list.

    EXAMPLES::

        sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
        sage: p2 = Polyhedron(vertices=[(1, 1), (0, 0)])
        sage: p3 = Polyhedron(vertices=[(0, 0)])
        sage: p4 = Polyhedron(vertices=[(1, 1)])
        sage: sage.geometry.polyhedral_complex.cells_list_to_cells_dict([p1, p2, p3, p4])
        {0: {A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex,
          A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex},
         1: {A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices},
         2: {A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices}}
    """
    cells_dict = {}
    for cell in cells_list:
        d = cell.dimension()
        if d in cells_dict:
            cells_dict[d].add(cell)
        else:
            cells_dict[d] = set([cell])
    return cells_dict


def exploded_plot(polyhedra, *,
                  center=None, explosion_factor=1, sticky_vertices=False,
                  sticky_center=True, point=None, **kwds):
    r"""
    Return a plot of several ``polyhedra`` in one figure with extra space
    between them.

    INPUT:

    - ``polyhedra`` -- an iterable of
      :class:`~sage.geometry.polyhedron.base.Polyhedron_base` objects

    - ``center`` -- (default: ``None``, denoting the origin) the center of
      explosion

    - ``explosion_factor`` -- (default: 1) a nonnegative number; translate
      polyhedra by this factor of the distance from ``center`` to their center

    - ``sticky_vertices`` -- (default: ``False``) boolean or dict; whether to
      draw line segments between shared vertices of the given polyhedra. A dict
      gives options for :func:`sage.plot.line`.

    - ``sticky_center`` -- (default: ``True``) boolean or dict. When ``center``
      is a vertex of some of the polyhedra, whether to draw line segments
      connecting the ``center`` to the shifted copies of these vertices. A dict
      gives options for :func:`sage.plot.line`.

    - ``color`` -- (default: ``None``) if ``'rainbow'``, assign a different
      color to every maximal cell and every vertex; otherwise, passed on to
      :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.plot`

    - other keyword arguments are passed on to
      :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.plot`

    EXAMPLES::

        sage: from sage.geometry.polyhedral_complex import exploded_plot
        sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
        sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
        sage: p3 = Polyhedron(vertices=[(0, 0), (1, 1), (2, 0)])
        sage: exploded_plot([p1, p2, p3])                                               # needs sage.plot
        Graphics object consisting of 20 graphics primitives
        sage: exploded_plot([p1, p2, p3], center=(1, 1))                                # needs sage.plot
        Graphics object consisting of 19 graphics primitives
        sage: exploded_plot([p1, p2, p3], center=(1, 1), sticky_vertices=True)          # needs sage.plot
        Graphics object consisting of 23 graphics primitives
    """
    from sage.plot.colors import rainbow
    from sage.plot.graphics import Graphics
    from sage.plot.line import line
    from sage.plot.point import point as plot_point
    import itertools

    polyhedra = list(polyhedra)
    g = Graphics()
    if not polyhedra:
        return g
    dim = polyhedra[0].ambient_dimension()
    if center is None:
        from sage.rings.rational_field import QQ
        center = vector(QQ, dim)
    else:
        center = vector(center)
    translations = [explosion_factor * ((p.center()
                                         + sum(r.vector() for r in p.rays()))
                                        - center)
                    for p in polyhedra]
    vertex_translations_dict = {}
    for P, t in zip(polyhedra, translations):
        for v in P.vertices():
            v = v.vector()
            v.set_immutable()
            vertex_translations_dict[v] = vertex_translations_dict.get(v, [])
            vertex_translations_dict[v].append(v + t)

    color = kwds.get('color')
    if color == 'rainbow':
        cell_colors_dict = dict(zip(polyhedra,
                                    rainbow(len(polyhedra))))
    for p, t in zip(polyhedra, translations):
        options = copy(kwds)
        if color == 'rainbow':
            options['color'] = cell_colors_dict[p]
        g += (p + t).plot(point=False, **options)

    if sticky_vertices or sticky_center:
        if sticky_vertices is True:
            sticky_vertices = dict(color='gray')
        if sticky_center is True:
            sticky_center = dict(color='gray')
        for vertex, vertex_translations in vertex_translations_dict.items():
            if vertex == center:
                if sticky_center:
                    for vt in vertex_translations:
                        g += line((center, vt), **sticky_center)
            else:
                if sticky_vertices:
                    for vt1, vt2 in itertools.combinations(vertex_translations, 2):
                        g += line((vt1, vt2), **sticky_vertices)
    if point is None:
        # default from sage.geometry.polyhedron.plot
        point = dict(size=10)
    if point is not False:
        if color == 'rainbow':
            vertex_colors_dict = dict(zip(vertex_translations_dict.keys(),
                                          rainbow(len(vertex_translations_dict.keys()))))
        for vertex, vertex_translations in vertex_translations_dict.items():
            options = copy(point)
            if color == 'rainbow':
                options['color'] = vertex_colors_dict[vertex]
            g += plot_point(vertex_translations, **options)
    return g
