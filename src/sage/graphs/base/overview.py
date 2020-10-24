r"""
Overview of (di)graph data structures

This module contains no code, and describes Sage's data structures for graphs
and digraphs. They can be used directly at Cython/C level, or through the
:class:`Graph` and :class:`DiGraph` classes (except one)

Data structures
---------------

Four data structures are natively available for (di)graphs in Sage:

- :mod:`~sage.graphs.base.sparse_graph` (default) -- for sparse (di)graphs, with
  a `\log(n)` edge test, and easy enumeration of neighbors. It is the most
  general-purpose data structure, though it can have a high memory cost in
  practice.

  - Supports: Addition/removal of edges/vertices, multiple edges, edge labels
    and loops.

- :mod:`~sage.graphs.base.dense_graph` -- for dense (di)graphs, with a `O(1)`
  edge test, and slow enumeration of neighbors.

  - Supports: addition/removal of edges/vertices, and loops.
  - Does not support: multiple edges and edge labels.

- :mod:`~sage.graphs.base.static_sparse_graph` -- for sparse (di)graphs and very
  intensive computations (at C-level). It is faster than
  :mod:`~sage.graphs.base.sparse_graph` in practice and *much* lighter in
  memory.

  - Supports: multiple edges, edge labels and loops
  - Does not support: addition/removal of edges/vertices.

- :mod:`~sage.graphs.base.static_dense_graph` -- for dense (di)graphs and very
  intensive computations (at C-level). It is mostly a wrapper over bitsets.

  - Supports: addition/removal of edges/vertices, and loops.
  - Does not support: multiple edges and edge labels.

For more information, see the data structures' respective pages.

The backends
------------

The :class:`Graph` and :class:`DiGraph` objects delegate the storage of vertices
and edges to other objects: the :mod:`graph backends
<sage.graphs.base.graph_backends>`::

    sage: Graph()._backend
    <sage.graphs.base.sparse_graph.SparseGraphBackend object at ...>

A (di)graph backend is a simpler (di)graph class having only the most elementary
methods (e.g.: add/remove vertices/edges). Its vertices can be arbitrary
hashable objects.

The only backend available in Sage is
:class:`~sage.graphs.base.c_graph.CGraphBackend`.

CGraph and CGraphBackend
------------------------

:class:`~sage.graphs.base.c_graph.CGraphBackend` is the backend of all native
data structures that can be used by :class:`Graph` and :class:`DiGraph`. It is
extended by:

- :class:`~sage.graphs.base.dense_graph.DenseGraphBackend`
- :class:`~sage.graphs.base.sparse_graph.SparseGraphBackend`
- :class:`~sage.graphs.base.static_sparse_backend.StaticSparseBackend`

While a :class:`~sage.graphs.base.c_graph.CGraphBackend` deals with arbitrary
(hashable) vertices, it contains a ``._cg`` attribute of type
:class:`~sage.graphs.base.c_graph.CGraph` which only deals with integer
vertices.

The :class:`~sage.graphs.base.c_graph.CGraph` data structures available in Sage
are:

- :class:`~sage.graphs.base.dense_graph.DenseGraph`
- :class:`~sage.graphs.base.sparse_graph.SparseGraph`
- :class:`~sage.graphs.base.static_sparse_backend.StaticSparseCGraph`

See the :mod:`~sage.graphs.base.c_graph` module for more information.
"""
