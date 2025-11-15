# -*- coding: utf-8 -*-
r"""
Plumbing Graphs

This file implements objects called *plumbing graphs*. These are graphs 
together with decorations on the vertices and the edges, and they encode
*plumbing manifolds*, three dimensional manifolds whose JSJ-pieces are
Seifert manifolds.

EXAMPLES::

    sage: P = PlumbingGraph()
    sage: P.add_seifert(-2,0,[2,2,2])
    0
    sage: P.intersection_matrix()
    [-2  1  1  1]
    [ 1 -2  0  0]
    [ 1  0 -2  0]
    [ 1  0  0 -2]
    sage: P.delete_vertex(0)
    sage: P.add_bamboo(1,j=1)
    0
    sage: P.make_normal()
    sage: P
    A plumbing graph with 2 vertices

AUTHORS:

- Baldur Sigurðsson (2018)
"""

#*****************************************************************************
#       Copyright (C) 2016 Baldur Sigur\dh sson  <baldursigurds@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************



from sage.structure.sage_object import SageObject
from sage.rings.rational import Rational
from sage.rings.integer_ring import ZZ
from sage.matrix.matrix_space import MatrixSpace
from copy import copy

def genus_hash(a,b):
    r"""
    Calculate the genus of a connected sum of two surfaces.

    The surfaces are not assumed orientable. A surface of genus g>=0 is
    the connected sum of g tori. A surface of genus g<0 is the connected
    sum of g real projective planes. The relation to note is that the
    connected sum of a torus and a real projective plane is the same as
    the connected sum of three copies of the real projective plane.

    INPUT:

    - ``a`` -- an integer
    - ``b`` -- an integer

    OUTPUT:

    An integer, the genus of the connected sum of two surfaces having
    genus a and b. See p. 304 of [Neu1981]_.

    EXAMPLES::

        sage: from sage.geometry.plumbing_graph import genus_hash
        sage: genus_hash(2,3)
        5
        sage: genus_hash(-2,3)
        -8
        sage: genus_hash(2,-3)
        -7
        sage: genus_hash(-2,-3)
        -5
    """
    if a >= 0 and b >= 0:
        return a + b
    if a < 0 and b >= 0:
        return a - 2*b
    if a >= 0 and b < 0:
        return -2*a + b
    if a < 0 and b < 0:
        return a + b

def big_genus_hash(S):
    r"""
    Calculate the genus of a connected sum of a finite number of surfaces.

    The surfaces are not assumed orientable. A surface of genus g>=0 is
    the connected sum of g tori. A surface of genus g<0 is the connected
    sum of g real projective planes. The relation to note is that the
    connected sum of a torus and a real projective plane is the same as
    the connected sum of three copies of the real projective plane.

    INPUT:
    
    - ``S`` -- a list of integers

    OUTPUT:
    
    An integer, the genus of the connected sum of surfaces with
    genera the entries of S.

    EXAMPLES::

        sage: from sage.geometry.plumbing_graph import big_genus_hash
        sage: big_genus_hash([1,2,3])
        6
        sage: big_genus_hash([1,2,-3])
        -9
    """
    r = 0
    for s in S:
        r = genus_hash(r,s)
    return r



class PlumbingGraph(SageObject):
    r"""
    A plumbing graph is a decorated graph which encodes graph manifolds.

    We implement the calculus described in [Neu1981]_. The plumbing
    construction creates a graph manifold (three dimensional) specified
    by a plumbing graph. Each vertex of the graph specifies an
    `S^1`-bundle over a surface, and for each edge, the corresponding
    bundles are glued together in a specific way. Each vertex is
    decorated by the Euler number of the corresponding `S^1`-bundle,
    this is ``mb``, the genus of the base of the fibration, this
    is ``g``, and the number of boundary components of the surface,
    ``r``.  Each edge is decorated by +1 or -1, which specify different
    gluings.
    
    The calculus of Neumann [Neu1981]_ consists of several operations, or
    moves, which can be applied to a plumbing graph to get the same
    manifold. The main result is that two graphs inducing homeomorphic
    (or, equivalently, diffeomorphic) plumbed manifolds are related by a
    finite sequence of these operations, or their inverses. Furthermore,
    each equivalence class of graphs by these operations contains a
    unique element, the minimal representative of this manifold, which
    can be obtained from any representative by applying the moves while
    possible.
    
    As a special case, the resolution graph of an isolated surface
    singularity is a plumbing graph which describes the link of the
    singularity. A plumbing graph is so obtained if and only if all
    genera are positive, all edges are positive, and the associated
    intersection matrix is negative definite.
    
    Another case related to singularity theory: As proved by Némethi and
    Szilárd [NSz2012]_, the boundary of the Milnor fiber of a reduced
    hypersurface singularity in complex 3-space is plumbed manifold.
    These do not necessarily have negative definite plumbing graphs, and
    may have negative edges as well.
    """ 

    def __init__(self, mb=[], g=[], r=[], p_edges=[], n_edges=[]):
        r"""
        Initialize ``self``.
    
        INPUT:
        
        - ``mb`` -- a list of integers, the Euler numbers
        - ``g`` -- a list of integers, the genera
        - ``r`` -- a list of nonnegative integers, the number of
          boundary comonents on the correpsonding base surface
        - ``p_edges`` -- a list of one or two- element sets of integers
        - ``n_edges`` -- a list of one or two- element sets of integers
        
        The lists mb, g and r should have the same lenght, say N. Each
        element in the lists p and n should be a set of one or two integers
        between 0 and N-1.
    
        EXAMPLES::

            sage: mb = [-2,-3]
            sage: g = [1,-1]
            sage: r = [0,1]
            sage: p = [{0,1}]
            sage: n = [{0,1}]
            sage: P = PlumbingGraph(mb,g,r,p,n)
            sage: P
            A plumbing graph with 2 vertices
        """
        if mb:
            if g:
                assert len(mb) == len(g),"string of genera is not of correct length"
            if r:
                assert len(mb) == len(r),"string of boundary components is not of correct length"
        self.vertices = set(range(len(mb)))
        self.mb = { i:mb[i] for i in range(len(mb)) }
        if g:
            self.g = { i:g[i] for i in range(len(mb)) }
        else:
            self.g = { i:0 for i in range(len(mb)) }
            
        if r:
            self.r = { i:r[i] for i in range(len(mb)) }
        else:
            self.r = { i:0 for i in range(len(mb)) }
            
        self.edges = set()
        self.adj = {}
        self.epsilon = {}
        
        for e in p_edges:
            self.add_edge(e, 1)
        for e in n_edges:
            self.add_edge(e, -1)
    
    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: P = PlumbingGraph
            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,0,[2,3/2,5/4])
            0
            sage: P
            A plumbing graph with 8 vertices
        """
        if not self.vertices:
            return "An empty plumbing graph"
        if len(self.vertices) == 1:
            return "A plumbing graph with one vertex"
            
        return "A plumbing graph with %s vertices"%(len(self.vertices))

########################################################################
# basic construction
    
    def _new_vertex_name(self):
        r"""
        Find a nonnegative integer which is not the name of any vertex.
        
        OUTPUT:
        
        A nonnegative integer.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P._new_vertex_name()
            0
            sage: P.add_vertex(-2,0,0)
            0
            sage: P._new_vertex_name()
            1
        """
        i = 0
        while i in self.vertices:
            i = i+1
        return i

    def add_vertex(self, mb, g, r):
        r"""
        Create a new vertex in the graph.
        
        INPUT:
    
        - ``mb`` -- the Euler number associated to the new vertex
        - ``g`` -- the genus associated to the new vertex
        - ``r`` -- the number of boundary components of the base surface
          associated to the new vertex

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P
            An empty plumbing graph
            sage: P.add_vertex(-2,0,0)
            0
            sage: P
            A plumbing graph with one vertex

        .. SEEALSO::
    
            :meth:`add_edge`, :meth:`add_bamboo`, :meth:`add_cycle`
            :meth:`add_seifert`
        """
        n = self._new_vertex_name()
        self.vertices = self.vertices.union({n})
        self.mb.update({n:mb})
        self.g.update({n:g})
        self.r.update({n:r})
        return n
    
    def delete_vertex(self, i):
        r"""
        Delete the vertex ``i``, and all adjacent edges.
        
        INPUT:
        
        - ``i`` -- name of the vertex to be deleted
        
        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_vertex(-2,0,1)
            0
            sage: P.delete_vertex(0)
        """
        if not i in self.vertices:
            print("Cannot delete a vertex which is not a vertex")
            return
        for e in self.edges:
            if i in self.adj[e]:
                self.delete_edge(e)
        del self.mb[i]
        del self.g[i]
        del self.r[i]
        self.vertices -= {i}

    def _new_edge_name(self):
        r"""
        Find a nonnegative integer which is not the name of any edge.

        OUTPUT:

        A nonnegative integer.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P._new_edge_name()
            0
            sage: P.add_vertex(-4,3,0)
            0
            sage: P.add_edge({0,0})
            0
            sage: P._new_edge_name()
            1
        """
        e = 0
        while e in self.edges:
            e = e+1
        return e
    
    def add_edge(self, A, epsilon=1):
        r"""
        Create a new edge in the graph.
        
        INPUT:
    
        - ``A`` -- a set of one or two vertices
        - ``epsilon`` -- either 1 or -1, the decoration of the new edge

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(-2,3,0)
            0
            sage: P.add_vertex(-2,0,0)
            1
            sage: P.add_edge({0})
            0
            sage: P.add_edge({1},-1)
            1
            sage: P.add_edge({0,1},1)
            2

        .. SEEALSO::

            :meth:`add_vertex`, :meth:`add_bamboo`, :meth:`add_cycle`
            :meth:`add_seifert`
        """
        if not A.issubset(self.vertices):
            print("Edge must join two vertices")
            return
        if not len(A) in {1,2}:
            print("Edge must join two vertices")
            return
        e = self._new_edge_name()
        self.edges = self.edges.union({e})
        self.adj.update({ e:A })
        self.epsilon.update({ e:epsilon })
        return e

    def delete_edge(self, e):
        r"""
        Delete an edge from the graph.
    
        INPUT:
        
        - ``e`` -- the name of the edge to be deleted

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(-2,3,0)
            0
            sage: P.add_vertex(-2,0,0)
            1
            sage: P.add_edge({0})
            0
            sage: P.add_edge({1},-1)
            1
            sage: P.add_edge({0,1},1)
            2
            sage: P.delete_edge(2)
        """
        if not e in self.edges:
            print("Cannot delete an edge which is not an edge")
        del self.adj[e]
        del self.epsilon[e]
        self.edges = self.edges.difference({e})

    def add_bamboo(self, x, j=None, r=0):
        r"""
        Create a bamboo in the graph.

        A bamboo is a string of rational (``g=0``) vertices. If the vertex
        ``j`` is specified, then this vertex is connected to the first
        vertex in the bamboo, otherwise the bamboo becomes a connected
        component. The function returns the last vertex in the bamboo.
        
        INPUT:
        
        - ``x`` -- a rational number
        - ``j`` -- a vertex of self
        - ``r`` -- a nonnegative integer

        OUTPUT:

        The last vertex of the bamboo.
        
        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(5/2)
            1
            sage: P
            A plumbing graph with 2 vertices

        .. SEEALSO::

            :meth:`add_vertex`, :meth:`add_edge`, :meth:`add_cycle`
            :meth:`add_seifert`
        """
        # Sage will at some point complain if x is an int. Silly. How do
        # rationals even work in this language?
        x = Rational(x)
        if not x.is_integer():
            i = self.add_vertex(-(x.ceil()),0,0)
        else:
            i = self.add_vertex(-(x.ceil()),0,r)

        if j is not None and j in self.vertices:
            self.add_edge({i,j})
        if not x.is_integer():
            return self.add_bamboo(1/(x.ceil() - x), i)
        else:
            return i

    def add_cycle(self, x, epsilon=1):
        r"""
        Create a cycle as a connected component.

        A cycle is a component of the graph, whose underlying graph is
        a cycle, genera and ``r`` vanish, and the Euler numbers are given
        by the negative continued fraction expansion of the rational
        number ``x``. All edges are positive if ``epsilon`` is 1,
        otherwise one edge is negative, and the other positive.

        INPUT:

        - ``x`` -- a rational number
        - ``epsilon`` -- an integer (1 or -1)

        OUTPUT:

        A vertex on the cycle.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_cycle(9/4)
            3
            sage: P
            A plumbing graph with 4 vertices

        .. SEEALSO::

            :meth:`add_vertex`, :meth:`add_edge`, :meth:`add_bamboo`
            :meth:`add_seifert`
        """
        i = self.add_bamboo(x)
        c = self.component(i)
        if len(c) == 1:
            A = {i}
        else:
            A = {j for j in c if self.degree(j) == 1}
        self.add_edge(A,epsilon)
        return i

    def add_seifert(self, mb, g, X, r=0):
        r"""
        Create a new star-shaped connected component in the graph.

        The central vertex of the star shaped component has Euler
        number ``mb``, genus ``g`` and number of boundary components ``r``.
        The bamboos of the star shaped graph have Euler numbers given
        by the negative continued fractin expansion of the elements
        of the list ``X``.
        
        
        INPUT:
        
        - ``mb`` -- an integer, the Euler number of the central vertex
        - ``g`` -- an integer, the genus of the central vertex
        - ``X`` -- a list of rational numbers, whose negative continued
          fraction expansion specifies the bamboos
        - ``r`` -- a nonnegative integer, the number of connected
          components of the surface associated with the central
          vertex

        OUTPUT:

        The central vertex of the star-shaped graph.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,4,[5/2, 9/4])
            0
            sage: P
            A plumbing graph with 7 vertices

        .. SEEALSO::

            :meth:`add_vertex`, :meth:`add_edge`, :meth:`add_bamboo`
            :meth:`add_cycle`
        """
        central = self.add_vertex(mb, g, r)
        for x in X:
            self.add_bamboo(x, central)
        return central

    def induced_subgraph(self, S):
        r"""
        Return the induced subgraph on the set of vertices ``S``.

        INPUT:

        - ``S`` -- a set of vertices

        OUTPUT:
        
        Return a PlumbingGraph object.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,4,[2,2])
            0
            sage: P.add_bamboo(4/3)
            5
            sage: P.induced_subgraph({3,4,5})
            A plumbing graph with 3 vertices
        """
        R = PlumbingGraph()
        L = list(S)
        D = { L[j]:j for j in range(len(L)) }
        for i in S:
            if i in self.vertices:
                R.add_vertex(self.mb[i], self.g[i], self.r[i])
        for e in self.edges:
            if self.adj[e] <= S:
                R.add_edge({D[i] for i in self.adj[e]}, self.epsilon[e])
        return R

    def export_to_graph(self):
        r"""
        Return a ``Graph()`` with the same vertices and edges as the
        plumbing graph.

        OUTPUT:

        A ``Graph()``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,0,[2,2,2])
            0
            sage: P.export_to_graph()
            Looped multi-graph on 4 vertices

        """
        from sage.graphs.graph import Graph
        G = Graph(loops=True, multiedges=True, weighted=True)
        for v in self.vertices:
            G.add_vertex(v)
        for e in self.edges:
            l = list(self.adj[e])
            if len(l) == 1:
                t = (l[0], l[0], self.epsilon[e])
            else:
                t = (l[0], l[1], self.epsilon[e])
            G.add_edge(t)
        return G

########################################################################
# some useful queries
    
    
    def is_loop(self, e):
        r"""
        Check whether an edge is a loop or not.
        
        INPUT:
        
        - ``e`` -- an edge
        
        OUTPUT:
        
        True or False, depending on whether ``e`` is a loop or not.

        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_vertex(-2,0,0)
            0
            sage: P.add_vertex(-2,0,0)
            1
            sage: P.add_edge({0})
            0
            sage: P.add_edge({0,1})
            1
            sage: P.is_loop(0)
            True
            sage: P.is_loop(1)
            False
        """
        assert e in self.edges,"{} is not an edge".format(e)
        return len(self.adj[e]) == 1

    def has_loop(self, i):
        r"""
        Check whether the vertex ``i`` has adjacent loops.

        INPUT:

        - ``i`` -- a vertex

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(-2,0,0)
            0
            sage: P.add_vertex(-2,0,0)
            1
            sage: P.add_edge({0})
            0
            sage: P.add_edge({0,1})
            1
            sage: P.has_loop(0)
            True
            sage: P.has_loop(1)
            False
        """
        for e in self.adjacent_edges(i):
            if self.is_loop(e):
                return True
        return False

    def degree(self, i):
        r"""
        Return the degree of the vertex ``i``.

        INPUT:
        
        - ``i`` -- a vertex

        OUTPUT:
        
        An integer, the degree of ``i``.

        EXAMPLES:
        
        We construct here the `E_8` graph, and check some degrees::
        
            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,0,[2,3/2,5/3])
            0
            sage: P.degree(0)
            3
            sage: P.degree(1)
            1
            sage: P.degree(2)
            2
        """
        d = 0
        for e in self.edges:
            A = self.adj[e]
            if len(A) == 1 and i in A:
                d += 2
            elif i in A:
                d += 1
        return d

    def neighbors(self, i):
        r"""
        Return the neighbours of a vertex ``i`` in a set.

        INPUT:

        - ``i`` -- a vertex
        
        OUTPUT:
        
        The set containing all vertices of self adjacent to ``i``. This
        does not include ``i``, even if ``i`` has a loop.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(-2,0,0)
            0
            sage: P.add_vertex(-2,0,0)
            1
            sage: P.add_vertex(-2,0,0)
            2
            sage: P.add_edge({0,1})
            0
            sage: P.add_edge({1,2})
            1
            sage: P.neighbors(1)
            {0, 2}
        """
        nbrs = set([])
        for e in self.edges:
            if i in self.adj[e]:
                nbrs |= self.adj[e]
        nbrs -= {i}
        return nbrs

    def neighbor(self, i):
        r"""
        Return some neighbor of the vertex ``i`` if one exists, otherwise
        ``-1``.

        INPUT:

        - ``i`` -- a vertex

        OUTPUT:
        
        A vertex (a nonnegative integer) or -1.

        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_vertex(-2,0,0)
            0
            sage: P.add_vertex(-2,0,0)
            1
            sage: P.add_vertex(-2,0,0)
            2
            sage: P.add_edge({1,2})
            0
            sage: P.neighbor(0)
            -1
            sage: P.neighbor(1)
            2
        """
        L = list(self.neighbors(i))
        if not L:
            return -1
        else:
            return L[0]

    def adjacent_edges(self, i):
        r"""
        Return a set of all edges adjacent to the vertex ``i``.

        INPUT:
        
        - ``i`` -- a vertex

        OUTPUT:
        
        A set containing those edges in the graph adjacent to the
        vertex ``i``.

        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,0,[5/4, 9/2, 5/3, 19/3])
            0
            sage: P.adjacent_edges(0)
            {0, 4, 6, 8}
        """
        return { e for e in self.edges if i in self.adj[e] }

    def intersection_matrix(self):
        r"""
        Return the associated intersection matrix.
        
        If ``i`` and ``j`` are distinct vertices, then the ``i,j`` entry
        of the matrix has value the sum of ``epsilon`` decorations of
        edges joining ``i`` and ``j``. The ``i`` th diagonal entry has value
        the Euler number of the ``i`` th vertex, plus twice the sum of
        ``epsilon`` decorations of loops adjacent to ``i``.

        OUTPUT:

        An `n\times n` integral matrix, where `n` is the number of vertices.

        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,4,[5/4,5/3,5/2])
            0
            sage: P.intersection_matrix()
            [-2  1  0  0  0  1  0  1  0]
            [ 1 -2  1  0  0  0  0  0  0]
            [ 0  1 -2  1  0  0  0  0  0]
            [ 0  0  1 -2  1  0  0  0  0]
            [ 0  0  0  1 -2  0  0  0  0]
            [ 1  0  0  0  0 -2  1  0  0]
            [ 0  0  0  0  0  1 -3  0  0]
            [ 1  0  0  0  0  0  0 -3  1]
            [ 0  0  0  0  0  0  0  1 -2]
        """
        n = len(self.vertices)
        V = list(self.vertices)
        W = { V[i]:i for i in range(0,n) }
        
        M = MatrixSpace(ZZ, n,n) 
        I  = copy(M.zero_matrix())
        for i in range(n):
            I[i,i] += self.mb[V[i]]
        for e in self.edges:
            L = list(self.adj[e])
            if len(L) == 1:
                L.append(L[0])
            I[W[L[0]],W[L[1]]] += self.epsilon[e]
            I[W[L[1]],W[L[0]]] += self.epsilon[e]
        return I

    def det(self):
        r"""
        Return the determinant of the graph.
        
        If `I` is the intersection matrix of the graph ``P``, then
        the determinant of ``P`` is the determinant of `-I`.
        
        OUTPUT:
        
        Integer.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,0,[2,3/2,3/2])
            0
            sage: P.det()
            3
            sage: P = PlumbingGraph()
            sage: P.add_seifert(2,0,[2,3/2,3/2])
            0
            sage: P.det()
            -69
        """
        
        return (-self.intersection_matrix()).det()

    def is_good_resolution(self):
        r"""
        Check whether the plumbing graph comes from a good resolution.

        This means that it is the graph associated with a resolution of
        an isolated complex surface singularity, thus describing the
        link of that singularity. This holds if and only if all genera
        are nonnegative, edges are positive, all ``r`` decorations are
        zero, and ``intersection_matrix()`` is negative definite.

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,0,[2,3/2,3/2])
            0
            sage: P.is_good_resolution()
            True
            sage: Q = PlumbingGraph()
            sage: Q.add_seifert(-2,0,[-2,-2,-2,-2])
            0
            sage: Q.is_good_resolution()
            False
        """
        for e in self.edges:
            if self.epsilon[e] == -1:
                return False

        for i in self.vertices:
            if self.g[i] < 0 or self.r[i] != 0:
                return False

        if not (-self.intersection_matrix()).is_positive_definite():
            return False

        return True

    def component(self, i):
        r"""
        Return the set of vertices in the connected component of ``i``.
        
        INPUT:
        
        - ``i`` -- a vertex
        
        OUTPUT:
        
        The set of vertices in the connected component of the vertex
        ``i``.
        
        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_bamboo(5/3)
            1
            sage: P.add_bamboo(9/8)
            9
            sage: P.component(0)
            {0, 1}
            sage: P.component(4)
            {2, 3, 4, 5, 6, 7, 8, 9}
        """
        if not i in self.vertices:
            return
        S = {i}
        keepgoing = True
        while keepgoing:
            T = { s for s in S }
            for j in S:
                T |= self.neighbors(j)
            if S == T:
                keepgoing = False
            S |= T
        return S

    def components(self):
        r"""
        Return a list of the components of the graph.

        OUTPUT:
        
        A list of sets of vertices.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,4,[5/3,10/3,11])
            0
            sage: P.add_bamboo(14/3)
            8
            sage: P.components()
            [{0, 1, 2, 3, 4, 5, 6}, {7, 8}]
        """
        S = {v for v in self.vertices}
        C = []
        while S:
            c = self.component(list(S)[0])
            S -= c
            C.append(c)
        return C
    
    def delete_component(self, i):
        r"""
        Delete the connected component of self containing the vertex ``i``.
        
        INPUT:
        
        - ``i`` -- a vertex

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,4,[5/3,10/3,11])
            0
            sage: P.add_bamboo(14/3)
            8
            sage: P.components()
            [{0, 1, 2, 3, 4, 5, 6}, {7, 8}]
            sage: P.delete_component(0)
            sage: P
            A plumbing graph with 2 vertices
        """
        C = self.component(i)
        for j in C:
            self.delete_vertex(j)

    def _neighbor_yields_leg(self, i, j):
        r"""
        This function checks whether, following the neighbour ``j`` of
        ``i``, we find a "leg", i.e. a bamboo with an end. This is used
        in e.g. ``N1()`` and ``N3()`` to discard exceptions.

        INPUT:

        - ``i`` -- a vertex
        - ``j`` -- a vertex (neighbor of i)

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_seifert(4,4,[6/5,7/5,4/3])
            0
            sage: P._neighbor_yields_leg(0,1)
            True
            sage: P._neighbor_yields_leg(5,4)
            False
        """
        if not {i, j} <= self.vertices:
            return False
        if not j in self.neighbors(i):
            return False
        S = {i}
        run = j
        while True:
            if self.mb[run] > -2 or self.g[run] != 0 or self.r[run] != 0 or self.degree(run) > 2:
                return False
            if self.degree(run) == 1:
                return True
            S |= {run}
            run = list(self.neighbors(run) - S)[0]
        return False

    def _minus_two_leaves(self, i):
        r"""
        Return the set of those neighbours which are (-2)-leaves

        A (-2)-leaf is a vertex of degree 1, rational, i.e.
        ``g == 0``, and ``r == 0`` and ``mb == -2``.
        It's used in ``N1()`` and ``N3()`` and ``make_normal_step_3()``.

        INPUT:

        - ``i`` -- a vertex

        OUTPUT:

        A set of vertices.

        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_seifert(0,-1,[2,2,3,4])
            0
            sage: P._minus_two_leaves(0)
            {1, 2}
        """
        return { j for j in self.neighbors(i)
            if  self.mb[j] == -2
            and self.g[j] == 0
            and self.r[j] == 0
            and self.degree(j) == 1 }

    def nodes(self):
        r"""
        Return the set of nodes.

        A node is a vertex with ``g != 0`` or ``r != 0`` or degree ``> 2``.

        OUTPUT:

        A set of vertices.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,3,[5/3,5/3,5/3])
            0
            sage: P.nodes()
            {0}
        """
        return { i for i in self.vertices if self.g[i] != 0 or self.r[i] != 0 or self.degree(i) > 2 or self.degree(i) == 0 }

    def max_chain(self, i):
        r"""
        Return the set of vertices on a maximal chain containing
        ``i``, if ``i`` is not a node, otherwise the emptyset.
        
        INPUT:

        - ``i`` -- a vertex

        OUTPUT:
        
        A set of vertices.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,3,[5/3,5/3,5/3])
            0
            sage: P.max_chain(0)
            set()
            sage: P.max_chain(1)
            {1, 2}
        """
        Nds = self.nodes()
        if i not in self.vertices - Nds:
            return set()
        R = {i}
        Nbrs = self.neighbors(i)
        while True:
            if Nbrs <= Nds:
                return R
            R |= (Nbrs - Nds)
            Nbrs = set()
            for r in R:
                Nbrs |= (self.neighbors(r) - R)
    
    def max_chain_is_cycle(self, i):
        r"""
        Check whether the maximal chain on which the vertex ``i``
        lies is a cycle.

        INTPUT:

        - ``i`` -- a vertex

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(8/7)
            6
            sage: P.add_cycle(8/7)
            13
            sage: P.max_chain_is_cycle(0)
            False
            sage: P.max_chain_is_cycle(13)
            True
        """
        mc = self.max_chain(i)
        if not mc:
            return False
        for j in mc:
            if self.degree(j) != 2 or not (self.neighbors(j) <= mc):
                return False
        return True


########################################################################
# The moves

    def find_candidate(self, R):
        r"""
        A wrapper for all the candidate finders.

        INPUT:
        
        - ``R`` -- an integer, 1-8

        OUTPUT:

        Either a candidate for move number ``R``, or ``-1`` if none exists.
        
        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(-2,0,0)
            0
            sage: P.add_vertex(-1,0,0)
            1
            sage: P.find_candidate(1)
            1
        """
        assert R in {1,2,3,4,5,6,7,8},"R should be an integer between 1 and 8"

        if R == 1:
            return self.find_R1_candidate()
        if R == 2:
            return self.find_R2_candidate()
        if R == 3:
            return self.find_R3_candidate()
        if R == 4:
            return self.find_R4_candidate()
        if R == 5:
            return self.find_R5_candidate()
        if R == 6:
            return self.find_R6_candidate()
        if R == 7:
            return self.find_R7_candidate()
        if R == 8:
            return self.find_R8_candidate()

    def R_move(self, R, i, j=-1):
        r"""
        A wrapper for all the R-moves.
        
        INPUT:
         
        - ``R`` -- an integer, 1-8
        - ``i`` -- a vertex
        
        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(-1,0,0)
            0
            sage: P.R_move(1,0)
        """
        assert R in {1,2,3,4,5,6,7,8},"R should be an integer between 1 and 8"

        if R == 1:
            self.blow_down(i)
        if R == 2:
            self.projective_absorb(i)
        if R == 3:
            self.zero_chain_absorb(i)
        if R == 4:
            self.unoriented_hanldle_absorb(i)
        if R == 5:
            self.oriented_hanldle_absorb(j,i)
        if R == 6:
            self.split(i)
        if R == 7:
            self.seifert_graph_exchange(i)
        if R == 8:
            self.annulus_absorb(i)


########################################################################
# R0 - reversing signs on edges

    def reverse_edge_signs_pos(self, i):
        r"""
        Reverse the sign of each edge adjacent to the vertex ``i``,
        except for loops. If ``i`` has genus ``<0``, then the function does
        nothing.
        
        INPUT:
        
        - ``i`` -- a vertex

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,0,[3/2,5/2,7/2])
            0
            sage: P.intersection_matrix()
            [-2  1  0  1  0  1  0]
            [ 1 -2  1  0  0  0  0]
            [ 0  1 -2  0  0  0  0]
            [ 1  0  0 -3  1  0  0]
            [ 0  0  0  1 -2  0  0]
            [ 1  0  0  0  0 -4  1]
            [ 0  0  0  0  0  1 -2]
            sage: P.reverse_edge_signs_pos(0)
            sage: P.intersection_matrix()
            [-2 -1  0 -1  0 -1  0]
            [-1 -2  1  0  0  0  0]
            [ 0  1 -2  0  0  0  0]
            [-1  0  0 -3  1  0  0]
            [ 0  0  0  1 -2  0  0]
            [-1  0  0  0  0 -4  1]
            [ 0  0  0  0  0  1 -2]
        """
        assert i in self.vertices,"{} is not a vertex".format(i)
        assert self.g[i] >= 0,"vertex {} has negative genus".format(i)
        for e in self.edges:
            if i in self.adj[e] and not len(self.adj[e]) == 1:
                self.epsilon[e] *= -1

    def reverse_edge_signs_neg(self, e):
        r"""
        Reverse the sign on an edge ``e`` if it is adjacent to a vertex
        with negative genus.
        
        INPUT:
        
        - ``e`` -- an edge

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,-1,[3/2,5/2,7/2])
            0
            sage: P.intersection_matrix()
            [-2  1  0  1  0  1  0]
            [ 1 -2  1  0  0  0  0]
            [ 0  1 -2  0  0  0  0]
            [ 1  0  0 -3  1  0  0]
            [ 0  0  0  1 -2  0  0]
            [ 1  0  0  0  0 -4  1]
            [ 0  0  0  0  0  1 -2]
            sage: P.reverse_edge_signs_neg(0)
            sage: P.intersection_matrix()
            [-2 -1  0  1  0  1  0]
            [-1 -2  1  0  0  0  0]
            [ 0  1 -2  0  0  0  0]
            [ 1  0  0 -3  1  0  0]
            [ 0  0  0  1 -2  0  0]
            [ 1  0  0  0  0 -4  1]
            [ 0  0  0  0  0  1 -2]
        """
        assert e in self.edges,"{} is not an edge".format(e)
        letsgo = 0
        for v in self.adj[e]:
            if self.g[v] < 0:
                letsgo = 1
        assert letsgo,"can only switch sign of an edge adjacent to a negative genus vertex"
        self.epsilon[e] *= -1
            

########################################################################
# R1 - blowing up and down
    
    def blow_up_disjoint(self, epsilon=-1):
        r"""
        The opposite to R1 in [Neu1981]_.

        Creates a new vertex in the graph with ``g == r == 0`` and
        ``mb == epsilon``.

        INPUT:
        
        - ``epsilon`` -- integer, either -1 or 1

        OUTPUT:

        Vertex.

        EXAMPLES:

        Create a plumbing graph having one (-1)-vertex. This corresponds
        to blowing up the origin in `\mathbb{C}^2`::

            sage: P = PlumbingGraph()
            sage: P.blow_up_disjoint(-1)
            0
            sage: P
            A plumbing graph with one vertex
        """
        assert epsilon*epsilon == 1,"epsilon must be -1 or 1"

        return self.add_vertex(epsilon,0,0)

    def blow_up_vertex(self, i, epsilon=-1):
        r"""
        The opposite to R1 in [Neu1981]_.

        Creates a new vertex in the graph
        with ``g = r = 0`` and ``mb = epsilon``, and an edge connecting it
        to ``i``, and modifies the Euler number of ``i``.

        INPUT:
        
        - ``i`` -- a vertex
        - ``epsilon`` -- -1 or 1

        OUTPUT:
        
        The new vertex.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(5/4)
            3
            sage: P.blow_up_vertex(3)
            4
        """
        assert epsilon*epsilon == 1,"epsilon must be -1 or 1"
        assert i in self.vertices,str(i)+" is not a vertex"

        j = self.add_vertex(epsilon,0,0)
        self.add_edge({i,j}, 1)
        self.mb[i] += epsilon
        return j

    def blow_up_edge(self, e, epsilon=-1):
        r"""
        The opposite to R1 in [Neu1981]_.

        Creates a new vertex in the graph
        with ''g = r = 0'' and ``mb = epsilon`` and edges connecting it to the
        vertices adjacent to ``e``, modifies those Euler numbers, and then
        deletes ``e``.

        INPUT:
        
        - ``e`` -- an edge
        - ``epsilon`` -- -1 or 1

        OUTPUT:

        The new vertex.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(3/2)
            1
            sage: P.blow_up_edge(0, epsilon=1)
            2
            sage: P
            A plumbing graph with 3 vertices
        """
        assert epsilon*epsilon == 1,"epsilon must be -1 or 1"
        assert e in self.edges,str(e)+" is not an edge"

        i = self.add_vertex(epsilon, 0, 0)
        epsilon0 = self.epsilon[e]
        epsilon1 = 1
        epsilon2 = (-1) * epsilon * epsilon0
        if self.is_loop(e):
            for j in self.adj[e]:
                self.mb[j] += 2*epsilon
                self.add_edge({i,j}, epsilon1)
                self.add_edge({i,j}, epsilon2)
        else:
            for j in self.adj[e]:
                self.mb[j] += epsilon
                f = self.add_edge({i,j}, epsilon1)
            self.epsilon[f] = epsilon2
        self.delete_edge(e)
        return i

    def blow_down(self, i):
        r"""
        R1 in [Neu1981]_. If possible, blows down the vertex ``i``.

        INPUT:
        
        - ``i`` -- a vertex

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(3/2)
            1
            sage: P.blow_up_edge(0, epsilon=1)
            2
            sage: P.blow_down(2)
        """
        assert i in self.vertices,str(i)+" is not a vertex"
        assert self.R1_candidate(i),"vertex {} cannot be blown down".format(i)

        if self.degree(i) == 0:
            self.delete_vertex(i)

        if self.degree(i) == 1:
            j = self.neighbor(i)
            self.mb[j] -= self.mb[i]
            self.delete_vertex(i)

        if self.degree(i) == 2:
            epsilon0 = -self.mb[i]
            for e in self.adjacent_edges(i):
                epsilon0 *= self.epsilon[e]
            Nbrs = self.neighbors(i)
            for j in Nbrs:
                self.mb[j] -= (3-len(Nbrs))*self.mb[i]
                
            self.add_edge(Nbrs, epsilon0)
            self.delete_vertex(i)

    def R1_candidate(self, i):
        r"""
        Check whether the vertex ``i`` can be blown down or not.
        
        INPUT:
        
        - ``i`` -- a vertex

        OUTPUT:

        ``True`` if ``i`` can be  blown down, otherwise ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(-1,0,0)
            0
            sage: P.R1_candidate(0)
            True
            sage: P.add_vertex(-1,3,0)
            1
            sage: P.R1_candidate(1)
            False
        """
        assert i in self.vertices,str(i)+" is not a vertex"

        if self.g[i] != 0:
            return False
        if self.mb[i]*self.mb[i] != 1:
            return False
        if self.r[i] != 0:
            return False
        if self.degree(i) > 2:
            return False
        if self.degree(i) == 2 and len(self.adjacent_edges(i)) == 1:
            return False
        return True
     
    def find_R1_candidate(self):
        r"""
        Return a vertex which can be blown down, if it exists,
        otherwise, return ``-1``.
        
        OUTPUT:
        
        Some vertex which can be blown down, if it exists, otherwise
        ``-1``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(-2,0,0)
            0
            sage: P.find_R1_candidate()
            -1
            sage: P.add_vertex(-1,0,0)
            1
            sage: P.find_R1_candidate()
            1
        """
        for i in self.vertices:
            if self.R1_candidate(i):
                return i
        return -1

########################################################################
# R2 - RP^2 absorption

    def projective_absorb(self, j, i):
        r"""
        Apply `RP^2`-absorption from [Neu1981]_, p. 305, to the graph.

        INPUT
        
        - ``j`` -- a vertex
        - ``i`` -- a vertex

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-1,0,[2,2,9/7])
            0
            sage: P
            A plumbing graph with 7 vertices
            sage: P.projective_absorb(0,3)
            sage: P
            A plumbing graph with 4 vertices
        """
        assert i in self.vertices,str(i)+" is not a vertex"
        assert j in self.vertices,str(j)+" is not a vertex"
        assert i in self.R2_candidate_neighbors(j),"vertices {} and {} cannot be projective-absorbed".format(j,i)

        for k in self.neighbors(j) - {i}:
            self.delete_vertex(k)
        self.delete_vertex(j)
        self.g[i] = genus_hash(self.g[i], -1)

    def R2_candidate(self, j):
        r"""
        Check whether R2 can be applied to the vertex ``j``.

        INPUT:

        - ``j`` -- a vertex

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-1,0,[2,2,9/7])
            0
            sage: P
            A plumbing graph with 7 vertices
            sage: P.projective_absorb(0,3)
            sage: P
            A plumbing graph with 4 vertices
            sage: P = PlumbingGraph()
            sage: P.add_seifert(-1,0,[2,2,9/7])
            0
            sage: P.R2_candidate(0)
            True
            sage: P.R2_candidate(1)
            False
        """
        assert j in self.vertices,str(j)+" is not a vertex"

        return self.R2_candidate_neighbors(j) != set()

    def R2_candidate_neighbors(self, j):
        r"""
        Let ``j`` be the central vertex on p. 305 in [Neu1981]_ with
        Euler number delta, and two (-2)-neihbors. The function returns
        the set of neighbors which may play the role of ``i`` in that
        picture.

        INPUT:

        - ``j`` -- a vertex

        OUTPUT:

        A set of vertices.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-1,0,[2,2,9/7])
            0
            sage: P.R2_candidate_neighbors(0)
            {3}
            sage: P.R2_candidate_neighbors(3)
            set()
        """
        assert j in self.vertices,str(j)+" is not a vertex"

        if self.g[j] != 0:
            return set()
        if self.r[j] != 0:
            return set()
        if self.degree(j) != 3 or len(self.neighbors(j)) != 3:
            return set()

        two_neighbors = { k for k in self.neighbors(j) if
                self.degree(k) == 1 and
                self.mb[k]*self.mb[k] == 4 and
                self.g[k] == 0 and
                self.r[k] == 0 }
        if len(two_neighbors) < 2:
            return set()
        R = set()
        for i in self.neighbors(j):
            if len(two_neighbors | {i}) == 3 and 4*self.mb[j] == sum([self.mb[tn] for tn in two_neighbors - {i}]):
                R |= {i}
        return R

    def find_R2_candidate(self):
        r"""
        Return a vertex on which R2 can be applied, if it exists.
        Otherwise return ``-1``.
        
        OUTPUT:
        
        Some vertex on which R2 can be applied, if it exists, otherwise,
        ``-1``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-1,0,[2,2,7/4])
            0
            sage: P.find_R2_candidate()
            0
        """
        for j in self.vertices:
            if self.R2_candidate(j):
                return j
        return -1
    
    def projective_extrude(self, i, g, d, d1, d2):
        r"""
        Inverse to R2.

        INPUT:

        - ``i`` -- a vertex
        - ``g`` -- an integer, the genus of i after extrusion
        - ``d`` -- an integer (delta on p. 305 of [Neu1981]_)
        - ``d1`` -- an integer (delta_1 on p. 305 of [Neu1981]_)
        - ``d2`` -- an integer (delta_2 on p. 305 of [Neu1981]_)
        
        EXAMPLES:
        
        We construct a Seifert graph with a central vertex with genus ``-1``,
        and extrude it::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(4, -1, [4/3, 5/2, 9/4])
            0
            sage: P
            A plumbing graph with 10 vertices
            sage: P.projective_extrude(0, 0, -1, -1, -1)
            sage: P
            A plumbing graph with 13 vertices
        """
        assert i in self.vertices,str(i)+" is not a vertex"
        assert self.g[i] == genus_hash(g, -1),"incompatible genus data"
        assert d1*d1 == 1 or d2*d2 == 1 or 2*d == d1+d2,"incompatible Euler numbers"
        
        self.g[i] = g
        a = self.add_vertex(d, 0, 0)
        self.add_edge({i,a})
        b = self.add_vertex(2*d1, 0, 0)
        self.add_edge({a,b})
        c = self.add_vertex(2*d2, 0, 0)
        self.add_edge({a,c})
       

########################################################################
# R3 - 0-chain absorption

    def zero_chain_absorb(self, k):
        r"""
        R3 in [Neu1981]_.

        Apply zero-chain absorption to the vertex ``k``.
        
        INPUT:
        
        - ``k`` -- a vertex

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(5/4)
            3
            sage: P.add_vertex(0,0,0)
            4
            sage: P.add_edge({3,4})
            3
            sage: P.add_bamboo(5/4, j=4)
            8
            sage: P
            A plumbing graph with 9 vertices
            sage: P.zero_chain_absorb(4)
            sage: P
            A plumbing graph with 7 vertices
        """
        assert k in self.vertices,"{} is not a vertex".format(k)
        assert self.R3_candidate(k),"vertex {} cannot be zero chain absorbed".format(k)

        l = self.add_vertex(
            sum([self.mb[i] for i in self.neighbors(k)]),
            big_genus_hash([self.g[i] for i in self.neighbors(k)]),
            sum([self.r[i] for i in self.neighbors(k)]))
        ij = list(self.neighbors(k))
        PR = 1
        for e in self.adjacent_edges(k):
            PR *= self.epsilon[e]
        epsilon_prime = [1, (-1)*PR]
        for a in range(2):
            for e in self.adjacent_edges(ij[a]):
                if self.is_loop(e):
                    self.add_edge({l}, self.epsilon[e])
                elif k not in self.adj[e]:
                    self.add_edge(
                        self.adj[e] - self.neighbors(k) | {l},
                        self.epsilon[e]*epsilon_prime[a])
        self.delete_vertex(k)
        self.delete_vertex(ij[0])
        self.delete_vertex(ij[1])
    
    def R3_candidate(self, k):
        r"""
        Check whether or not R3 can be applied to the vertex ``k``.

        INPUT:
        
        - ``k`` -- a vertex

        OUTPUT:

        ``True`` or ``False``, depending on whether 0-chain absorption
        is applicable to the vertex ``k``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(5/4)
            3
            sage: P.add_vertex(0,0,0)
            4
            sage: P.add_edge({3,4})
            3
            sage: P.add_bamboo(5/4, j=4)
            8
            sage: P.R3_candidate(0)
            False
            sage: P.R3_candidate(4)
            True
        """
        assert k in self.vertices,"{} is not a vertex".format(k)
        if self.mb[k] != 0:
            return False
        if self.g[k] != 0:
            return False
        if self.r[k] != 0:
            return False
        if self.degree(k) != 2:
            return False
        if len(self.neighbors(k)) != 2:
            return False
        return True

    def find_R3_candidate(self):
        r"""
        Return a vertex on which R3 can applied, if it exists,
        otherwise, return ``-1``.
        
        OUTPUT:
        
        Some vertex on which R3 can be applied, if it exists, otherwise,
        ``-1``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(5/4)
            3
            sage: P.add_vertex(0,0,0)
            4
            sage: P.add_edge({3,4})
            3
            sage: P.add_bamboo(5/4, j=4)
            8
            sage: P.find_R3_candidate()
            4
            sage: P.zero_chain_absorb(4)
            sage: P.find_R3_candidate()
            -1
        """
        for k in self.vertices:
            if self.R3_candidate(k):
                return k
        return -1

    def zero_chain_extrude(self, l, MB, G, R, E):
        r"""
        Inverse to R3 from [Neu1981]_, zero-chain extrusion.

        ``l`` is the central vertex on the picture to the right on p. 305.
        ``i`` and ``j`` are labelled on the picture to the left, and
        ``k`` is the
        zero-vertex between them. ``MB``, ``G``, and ``R`` are lists with the
        corresponding labels for ``i`` and ``j``. The function checks
        whether this data makes sense. ``E`` is a list of
        two sets containing the edges which should be connected to
        ``i`` and ``j``. To make things simple, we decorate the new edges
        connecting ``i``,``k`` and ``k``,``j`` with ``+`` and ``-``,
        respectively.
         
        INPUT:
         
        - ``l`` -- a vertex
        - ``MB`` -- list of two integers
        - ``G`` -- list of two integers
        - ``R`` -- list of two nonnegative integers
        - ``E`` -- a list of two sets of edges
         
        OUTPUT:
         
        A list of the vertices ``[i,k,j]``, where ``i`` and ``j``
        are as in the
        picture on p. 305 of [Neu1981]_, and ``k`` is the newly extruded
        zero-vertex.
        
        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,-3,[3/2,3/2,3/2,3/2],r=3)
            0
            sage: P.add_edge({0})
            8
            sage: P.adjacent_edges(0)
            {0, 2, 4, 6, 8}
            sage: P.zero_chain_extrude(0,[-1,-1],[1,-1],[1,2],[{0,2,8},{4,6,8}])
            [9, 11, 10]
        """
        assert l in self.vertices,"{} is not a vertex".format(l)
        assert len(MB) == 2 and MB[0] + MB[1] == self.mb[l],"Euler numbers do not add up"
        assert len(G) == 2 and genus_hash(G[0], G[1]) == self.g[l],"genera do not add up"
        assert len(R) == 2 and R[0] + R[1] == self.r[l],"r's do not add up"
        assert len(E) == 2 and E[0] | E[1] == self.adjacent_edges(l),"E must provide a covering of the adjacent edges"
        for e in E[0] & E[1]:
            assert self.is_loop(e),"edge {} is not a loop".format(e)

        L = [
            [e for e in E[0]-E[1] if self.is_loop(e)],
            [e for e in E[1]-E[0] if self.is_loop(e)]]
        if len(E[0]) + len(L[0]) + len(E[1]) + len(L[1]) != self.degree(l):
            return []
        i = self.add_vertex(MB[0],G[0],R[0])
        j = self.add_vertex(MB[1],G[1],R[1])
        k = self.add_vertex(0,0,0)
        self.add_edge({i,k},1)
        self.add_edge({k,j},-1)
        for e in E[0]:
            if self.is_loop(e):
                if e in E[1]:
                    self.add_edge({i,j}, self.epsilon[e])
                else:
                    self.add_edge({i}, self.epsilon[e])
            else:
                self.add_edge(self.adj[e] - {l} | {i}, self.epsilon[e])
        for e in E[1]:
            if not (self.is_loop(e) and e in E[0]):
                self.add_edge(self.adj[e] - {l} | {j}, self.epsilon[e])
        self.delete_vertex(l)
        return [i,k,j]

########################################################################
# R4 - unoriented handle absorption

    def unoriented_handle_absorb(self, j):
        r"""
        R4 in [Neu1981]_.

        Apply unoriented handle absorption to the vertex ``j``.

        INTPUT:

        - ``j`` -- a vertex

        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_vertex(5,3,2)
            0
            sage: P.add_vertex(0,0,0)
            1
            sage: P.add_edge({0,1})
            0
            sage: P.add_edge({0,1})
            1
            sage: P.unoriented_handle_absorb(1)
            sage: P
            A plumbing graph with one vertex
        """
        assert j in self.vertices,"{} is not a vertex".format(j)
        assert self.mb[j] == 0,"vertex must have Euler number zero"
        assert self.g[j] == 0,"vertex must have genus zero"
        assert self.r[j] == 0,"vertex must have r zero"
        assert not self.has_loop(j),"vertex cannot have loop"
        assert self.degree(j) == 2,"vertex must have degree 2"
        assert len(self.neighbors(j)) == 1,"vertex must have a unique neighbor"
        assert len({self.epsilon[e] for e in self.adjacent_edges(j)}) == 1,"signs must agree on adjacent edges"

        i = self.neighbor(j)
        self.g[i] = genus_hash(self.g[i], -2)
        self.delete_vertex(j)

    def R4_candidate(self, j):
        r"""
        Check whether or not R4 can be applied to the vertex ``j``.

        INPUT:
        
        - ``j`` -- a vertex

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(5,3,2)
            0
            sage: P.add_vertex(0,0,0)
            1
            sage: P.add_edge({0,1})
            0
            sage: P.add_edge({0,1})
            1
            sage: P.R4_candidate(1)
            True
            sage: P.R4_candidate(0)
            False
        """
        assert j in self.vertices,"{} is not a vertex".format(j)

        if self.mb[j] != 0:
            return False
        if self.g[j] != 0:
            return False
        if self.r[j] != 0:
            return False
        if self.degree(j) != 2:
            return False
        if len(self.neighbors(j)) != 1 or j in self.neighbors(j):
            return False
        if len({self.epsilon[e] for e in self.adjacent_edges(j)}) != 1:
            return False
        return True

    def find_R4_candidate(self):
        r"""
        Return a vertex on which R4 can applied, if it exists,
        otherwise, return ``-1``.
        
        OUTPUT:
        
        Some vertex on which R4 can be applied, if it exists, otherwise,
        ``-1``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(5,3,2)
            0
            sage: P.add_vertex(0,0,0)
            1
            sage: P.add_edge({0,1})
            0
            sage: P.add_edge({0,1})
            1
            sage: P.find_R4_candidate()
            1
        """
        for j in self.vertices:
            if self.R4_candidate(j):
                return j
        return -1


########################################################################
# R5 - oriented handle absorption

    def oriented_handle_absorb(self, j):
        r"""
        R5 in [Neu1981]_.

        Apply oriented handle absorption to the vertex ``j``.
        
        INPUT:
        
        - ``j`` -- a vertex

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(5,3,2)
            0
            sage: P.add_vertex(0,0,0)
            1
            sage: P.add_edge({0,1})
            0
            sage: P.add_edge({0,1}, epsilon=-1)
            1
            sage: P.oriented_handle_absorb(1)
            sage: P
            A plumbing graph with one vertex
        """
        assert j in self.vertices,"{} is not a vertex".format(j)
        assert self.mb[j] == 0,"vertex must have Euler number zero"
        assert self.g[j] == 0,"vertex must have genus zero"
        assert self.r[j] == 0,"vertex must have r zero"
        assert not self.has_loop(j),"vertex cannot have loop"
        assert self.degree(j) == 2,"vertex must have degree 2"
        assert len(self.neighbors(j)) == 1,"vertex must have a unique neighbor"
        assert len({self.epsilon[e] for e in self.adjacent_edges(j)}) == 2,"signs must not agree on adjacent edges"

        i = self.neighbor(j)
        self.g[i] = genus_hash(self.g[i], 1)
        self.delete_vertex(j)

    def R5_candidate(self, j):
        r"""
        Check whether or not R5 can be applied to the vertex ``j``.

        INPUT:
        
        - ``j`` -- a vertex
        
        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(5,3,2)
            0
            sage: P.add_vertex(0,0,0)
            1
            sage: P.add_edge({0,1})
            0
            sage: P.add_edge({0,1}, epsilon=-1)
            1
            sage: P.R5_candidate(0)
            False
            sage: P.R5_candidate(1)
            True
        """
        assert j in self.vertices,"{} is not a vertex".format(j)
        if self.mb[j] != 0:
            return False
        if self.g[j] != 0:
            return False
        if self.r[j] != 0:
            return False
        if self.degree(j) != 2:
            return False
        if len(self.neighbors(j)) != 1 or j in self.neighbors(j):
            return False
        if len({self.epsilon[e] for e in self.adjacent_edges(j)}) != 2:
            return False
        return True

    def find_R5_candidate(self):
        r"""
        Return a vertex on which R5 can be applied, if it exists,
        otherwise, return ``-1``.
        
        OUTPUT:
        
        Some vertex on which R5 can be applied, if it exists, otherwise,
        ``-1``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(5,3,2)
            0
            sage: P.add_vertex(0,0,0)
            1
            sage: P.add_edge({0,1})
            0
            sage: P.add_edge({0,1}, epsilon=-1)
            1
            sage: P.find_R5_candidate()
            1
        """
        for j in self.vertices:
            if self.R5_candidate(j):
                return j
        return -1

        
########################################################################
# R6 - splitting

    def split(self, j):
        r"""
        R6 in [Neu1981]_.

        Apply splitting to he graph.
        ``j`` is the vertex with Euler number zero in the picture on
        page 305 of [Neu1981]_.

        INPUT:
        
        - ``j`` -- a vertex

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,0,[2,2,2])
            0
            sage: P.add_seifert(-2,0,[2,2,2])
            4
            sage: P.add_vertex(-2,3,5)
            8
            sage: P.add_edge({8,0})
            6
            sage: P.add_edge({8,4})
            7
            sage: P.add_vertex(0,0,0)
            9
            sage: P.add_edge({9,8})
            8
            sage: P.split(9)
            sage: P
            A plumbing graph with 19 vertices
        """
        assert j in self.vertices,"{} is not a vertex".format(j)
        assert self.R6_candidate(j),"vertex {} cannot be split".format(j)

        i = self.neighbor(j)
        g = self.g[i]
        r = self.r[i]
        d = self.degree(i)
        self.delete_vertex(j)
        self.delete_vertex(i)
        if g >= 0:
            k = 2*g + (d-1) - len(self.components())
        if g < 0:
            k = -g + (d-1) - len(self.components())
        for a in range(k):
            self.add_bamboo(0)
        for a in range(r):
            self.add_bamboo(0,r=1)

    def R6_candidate(self, j):
        r"""
        Check whether or not R6 can be applied to the vertex ``j``.

        INPUT:
        
        - ``j`` -- a vertex

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-4,2,[0,1,3/4])
            0
            sage: P.R6_candidate(0)
            False
            sage: P.R6_candidate(1)
            True
        """
        assert j in self.vertices,"{} is not a vertex".format(j)

        if self.mb[j] != 0:
            return False
        if self.g[j] != 0:
            return False
        if self.r[j] != 0:
            return False
        if self.degree(j) != 1:
            return False
        i = list(self.neighbors(j))[0]
        if self.has_loop(i):
            return False
        return True

    def find_R6_candidate(self):
        r"""
        Return a vertex on which R6 can applied, if it exists,
        otherwise, return ``-1``.
        
        OUTPUT:
        
        Some vertex on which R6 can be applied, if it exists, otherwise,
        ``-1``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-4,2,[0,1,3/4])
            0
            sage: P.find_R6_candidate()
            1
        """
        for j in self.vertices:
            if self.R6_candidate(j):
                return j
        return -1


########################################################################
# R7 - Seifert graph exchange

    def seifert_graph_exchange(self, i):
        r"""
        R7 in [Neu1981]_.
        
        Replace component with a single vertex and loop with the
        appropriate star-shaped graph, see p. 306 of [Neu1981]_.
        
        INPUT:
        
        - ``i`` -- a vertex

        EXAMPLES:

        Example a::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(-1,0,0)
            0
            sage: P.add_edge({0})
            0
            sage: P.seifert_graph_exchange(0)
            sage: P
            A plumbing graph with 4 vertices

        Example b::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(0,0,0)
            0
            sage: P.add_edge({0})
            0
            sage: P.seifert_graph_exchange(0)
            sage: P
            A plumbing graph with 4 vertices
        
        Example c::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(1,0,0)
            0
            sage: P.add_edge({0})
            0
            sage: P.seifert_graph_exchange(0)
            sage: P
            A plumbing graph with 4 vertices
        
        Example d::
        
            sage: P = PlumbingGraph()
            sage: P.add_vertex(-1,0,0)
            0
            sage: P.add_edge({0},epsilon=-1)
            0
            sage: P.seifert_graph_exchange(0)
            sage: P
            A plumbing graph with 7 vertices
        
        Example e::
        
            sage: P = PlumbingGraph()
            sage: P.add_vertex(0,0,0)
            0
            sage: P.add_edge({0},epsilon=-1)
            0
            sage: P.seifert_graph_exchange(0)
            sage: P
            A plumbing graph with 8 vertices
        
        Example f::
        
            sage: P = PlumbingGraph()
            sage: P.add_vertex(1,0,0)
            0
            sage: P.add_edge({0},epsilon=-1)
            0
            sage: P.seifert_graph_exchange(0)
            sage: P
            A plumbing graph with 9 vertices
        """
        
        assert i in self.vertices,"{} is not a vertex".format(i)
        assert self.R7_candidate(i),"this vertex cannot be Seifert exchanged"

        mb = self.mb[i]
        epsilon = self.epsilon[list(self.adjacent_edges(i))[0]]
        
        self.delete_vertex(i)
        if mb == -1 and epsilon == 1:
            self.add_seifert(-1, 0, [2,3,6])
        if mb == 0 and epsilon == 1:
            self.add_seifert(-1, 0, [2,4,4])
        if mb == 1 and epsilon == 1:
            self.add_seifert(-1, 0, [3,3,3])
        if mb == -1 and epsilon == -1:
            self.add_seifert(-2, 0,
                [Rational(3)/Rational(2)
                ,Rational(3)/Rational(2)
                ,Rational(3)/Rational(2)])
        if mb == 0 and epsilon == -1:
            self.add_seifert(-2, 0,
                [2
                ,Rational(4)/Rational(3)
                ,Rational(4)/Rational(3)])
        if mb == 1 and epsilon == -1:
            self.add_seifert(-2, 0,
                [2
                ,Rational(3)/Rational(2)
                ,Rational(6)/Rational(5)])
            
    def R7_candidate(self, i):
        r"""
        Check whether or not R7 can be applied to the vertex ``i``.

        This means that ``i`` should have Euler number ``mb = -1,0,1``,
        genus ``g = 0`` and ``r = 0``, and one loop, and no other
        adjacent edges.

        INPUT:
        
        - ``i`` -- a vertex

        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_vertex(0,0,0)
            0
            sage: P.add_vertex(-2,0,0)
            1
            sage: P.add_edge({0})
            0
            sage: P.R7_candidate(0)
            True
            sage: P.R7_candidate(1)
            False
        """
        assert i in self.vertices,"{} is not a vertex".format(i)

        if self.mb[i]*self.mb[i] > 1:
            return False
        if self.g[i] != 0:
            return False
        if self.r[i] != 0:
            return False
        if self.degree(i) != 2:
            return False
        if not self.has_loop(i):
            return False
        return True

    def find_R7_candidate(self):
        r"""
        Return a vertex on which R7 can applied, if it exists,
        otherwise, return ``-1``.
        
        OUTPUT:
        
        Some vertex on which R7 can be applied, if it exists, otherwise,
        ``-1``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(0,0,0)
            0
            sage: P.add_vertex(-2,0,0)
            1
            sage: P.add_edge({0})
            0
            sage: P.find_R7_candidate()
            0
        """
        for i in self.vertices:
            if self.R7_candidate(i):
                return i
        return -1


########################################################################
# R8 - annulus absorption

    def annulus_absorb(self, j):
        r"""
        R8 in [Neu1981]_.

        Removes an edge and adds a boundary component to its unique
        neighbor, when appropriate.

        INPUT:
        
        - ``j`` -- a vertex

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(5/4)
            3
            sage: P.add_vertex(0,0,1)
            4
            sage: P.add_edge({4,3})
            3
            sage: P.annulus_absorb(4)
        """
        assert j in self.vertices,"{} is not a vertex".format(j)
        assert self.R8_candidate(j),"vertex {} cannot be annulus absorbed".format(j)

        i = list(self.neighbors(j))[0]
        self.r[i] += 1
        self.delete_vertex(j)

    def R8_candidate(self,j):
        r"""
        Check whether or not R8 can be applied to the vertex ``j``.

        This means that ``j`` should have genus zero and ``r=1``,
        and it should have a unique neighbour and no loop.

        INPUT:
        
        - ``j`` -- a vertex

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(5/4)
            3
            sage: P.add_vertex(0,0,1)
            4
            sage: P.add_edge({4,3})
            3
            sage: P.R8_candidate(0)
            False
            sage: P.R8_candidate(4)
            True
        """
        assert j in self.vertices,"{} is not a vertex".format(j)

        if self.g[j] != 0:
            return False
        if self.r[j] != 1:
            return False
        if self.degree(j) != 1:
            return False
        return True

    def find_R8_candidate(self):
        r"""
        Return a vertex on which R8 can be applied, if it exists,
        otherwise, return ``-1``.
        
        OUTPUT:
        
        Some vertex on which R8 can be applied, if it exists, otherwise,
        ``-1``.

        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_bamboo(5/4)
            3
            sage: P.add_vertex(0,0,1)
            4
            sage: P.add_edge({4,3})
            3
            sage: P.R8_candidate(0)
            False
            sage: P.R8_candidate(4)
            True
            sage: P.find_R8_candidate()
            4
        """
        for j in self.vertices:
            if self.R8_candidate(j):
                return j
        return -1


########################################################################
# Normal form check

    def is_normal(self):
        r"""
        Check whether the graph is in normal form or not.

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES:
        
        The Dynkin diagram `D_4` is in normal form::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,0,[2,2,2])
            0
            sage: P.is_normal()
            True
        
        A minus one vertex is not::
        
            sage: P = PlumbingGraph()
            sage: P.blow_up_disjoint()
            0
            sage: P.is_normal()
            False
        """
        return self.N1() and self.N2() and self.N3() and self.N4() and self.N5() and self.N6()

########################################################################
# N1

    def N1(self):
        r"""
        Check if condition N1 of [Neu1981]_ holds.

        OUTPUT:

        ``True`` or ``False``, depending on whether N1 holds or not.

        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_bamboo(3/2)
            1
            sage: P.blow_up_edge(0)
            2
            sage: P.N1()
            False
        """
        for i in self.vertices:
            if self.R1_candidate(i):
                return False
            if self.R2_candidate(i) and not self.N1_exception(i):
                return False
            if self.R3_candidate(i):
                return False
            if self.R4_candidate(i):
                return False
            if self.R5_candidate(i):
                return False
            if self.R6_candidate(i):
                return False
            if self.R7_candidate(i):
                return False
            if self.R8_candidate(i):
                return False
        return True

    def N1_exception(self, i):
        r"""
        Check whether the vertex ``i`` is the (-1)-vertex in the
        exceptional graph described in N1, p. 311 of [Neu1981]_.

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_seifert(-1,0,[2,2,9/4])
            0
            sage: P.N1_exception(0)
            True
            sage: P.N1_exception(1)
            False
        """
        # first check basic info about the central vertex:
        if self.mb[i] != -1:
            return False
        if self.g[i] != 0:
            return False
        if self.r[i] != 0:
            return False
        if self.degree(i) != 3 or len(self.neighbors(i)) != 3:
            return False
        # next see how many neighbors are rational leaves with Euler
        # number -2 and r=0
        mtl = self._minus_two_leaves(i)
        # if there are three such leaves, then it's good
        if len(mtl) == 3:
            return True
        # if there are fewer than two, then no
        if len(mtl) < 2:
            return False
        # finally, assuming that we have exactly two leaves as above,
        # investigate the third neighbor, it should be the end of a bamboo
        j = list(self.neighbors(i) - mtl)[0]
        return self._neighbor_yields_leg(i,j)

########################################################################
# N2

    def N2(self):
        r"""
        Check if condition N2 from p. 311 of [Neu1981]_ holds.
        
        OUTPUT:
        
        ``True`` or ``False``.
        
        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_seifert(4,0,[5/3,5/2,5/4])
            0
            sage: P.N2()
            True
            sage: P.mb[1] = 1
            sage: P.N2()
            False
        """
        for i in self.vertices - self.nodes():
            if self.mb[i] > -2:
                return False
        return True

########################################################################
# N3

    def N3(self):
        r"""
        Check if condition N3 from [Neu1981]_ holds.

        OUTPUT:
        
        ``True`` or ``False``.

        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,0,[2,2,3/2])
            0
            sage: P.N3()
            True
            sage: P.add_edge({3,4})
            4
            sage: P.N3()
            False
        """
        for i in self.vertices:
            if self.N3_obstruction(i):
                return False
        return True

    def N3_obstruction(self, i):
        r"""
        Check if the vertex ``i`` is the vertex with Euler number
        ``e`` in the picture on p. 311 and 312 in [Neu1981]_ which
        obstructs N3.

        INPUT:

        - ``i`` -- a vertex

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,0,[2,2,3/2])
            0
            sage: P.N3_obstruction(0)
            False
            sage: P.add_edge({3,4})
            4
            sage: P.N3_obstruction(0)
            True
        """
        mtn = self._minus_two_leaves(i)
        nbrs = self.neighbors(i)
        if self.degree(i) == len(nbrs) == 3 and len(mtn) >= 2 and self.g[i] == self.r[i] == 0:
            if len(mtn) == 3:
                return False
            j = list(nbrs - mtn)[0]

            for k in self.component(i):
                if k != i and self.mb[k] > -2:
                    return True

            return not self._neighbor_yields_leg(i,j)
        return False

########################################################################
# N4

    def N4(self):
        r"""
        Check if condition N4 from [Neu1981]_ holds.

        Remark: I interpret the phrase "vertex ``i`` is an interior vertex
        of a chain" to mean that ``i`` is not a node, that is, ``g=r=0``
        and
        degree ``<= 2``. In particular, this may mean that the component we
        are working with has only the two vertices visible on the
        picture on p. 312 of [Neu1981]_. This makes sense because of the
        following reasoning. There are two cases to consider, if the
        vertex ``i`` is not a node. In the first case, ``i`` lies on a
        chain connecting the vertex with ``g=-1`` to a node. In this case,
        step 4 would produce something violating N3. We do not want
        that. The other case is when the chain on which vertex i lies
        does not connect the vertex with ``g=-1`` to a node. This case is
        handled by N5 and step 5.

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES:
        
        We construct an example that is the exception, i.e. the neighbor
        of the vertex ``0`` is an interior vertex in a chain. Then we change
        it::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(0,-1,[4/3])
            0
            sage: P.N4()
            True
            sage: P.g[1] = 1
            sage: P.N4()
            False
        """
        for j in self.vertices:
            if self.N4_obstruction(j):
                return False
        return True

    def N4_obstruction(self, j):
        r"""
        Check whether the vertex ``j`` violates N4.

        INPUT:
        
        -``j`` -- a vertex

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(0,-1,[4/3])
            0
            sage: P.N4_obstruction(0)
            False
            sage: P.g[1] = 1
            sage: P.N4_obstruction(0)
            True
        """
        if self.mb[j] == 0 and self.g[j] == -1 and self.r[j] == 0 and self.degree(j) == 1:
            i = self.neighbor(j)
            if self.g[i] != 0 or self.r[i] != 0 or self.degree(i) > 2:
                return True
        return False

########################################################################
# N5
    def N5(self):
        r"""
        Check if condition N5 from [Neu1981]_ holds.

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(4,-1,[8/3])
            0
            sage: P.N5()
            False
            sage: P.vertices
            {0, 1, 2}
            sage: P.add_edge({2})
            2
            sage: P.N5()
            True
        """
        for j in self.vertices:
            if self.N5_obstruction(j):
                return False
        return True
    
    def N5_obstruction(self, j):
        r"""
        Check if ``j`` is the vertex with genus ``-1`` in the drawing in
        N5 on p. 312 of [Neu1981]_.

        INPUT:
        
        -``j`` -- a vertex

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(4,-1,[8/3])
            0
            sage: P.N5_obstruction(0)
            True
            sage: P.add_edge({2})
            2
            sage: P.N5_obstruction(0)
            False
        """
        if self.g[j] == -1 and self.r[j] == 0 and self.degree(j) == 0:
            return True
        
        if self.g[j] == -1 and self.r[j] == 0 and self.degree(j) == 1:
            i = self.neighbor(j)
            if self._neighbor_yields_leg(j,i):
                return True
        return False

########################################################################
# N6
    
    def N6(self):
        r"""
        Check if condition N6 from [Neu1981]_ holds.

        OUTPUT:

        ``True`` or ``False``.
        
        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_cycle(5/4)
            3
            sage: P.N6()
            False
        """
        
        if self.N6_ab_find_candidate() != -1:
            return False
        if self.N6_c_find_candidate() != -1:
            return False
        if self.N6_d_find_candidate() != -1:
            return False
        if self.N6_e_find_candidate() != -1:
            return False
        if self.N6_f_find_candidate() != -1:
            return False
        if self.N6_g_find_candidate() != -1:
            return False
        return True


    def N6_ab_candidate(self, i):
        r"""
        Check if the component containing the vertex ``i`` is a cycle of
        vertices with ``mb = -2``, ``g=0`` and ``r=0``.

        INPUT:

        - ``i`` -- a vertex

        OUTPUT:

        ``True`` or ``False``.
        
        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_cycle(5/4)
            3
            sage: P.N6_ab_candidate(0)
            True
        """
        if self.mb[i] == -2 and self.g[i] == 0 and self.r[i] == 0 and self.degree(i) == 2 and self.has_loop(i):
            return True
        j = i
        S = {i}
        while self.mb[j] == -2 and self.g[j] == 0 and self.r[j] == 0 and self.degree(i) == 2 and not self.has_loop(i):
            if self.neighbors(j) <= S:
                return True
            else:
                S |= {j}
                j = list(self.neighbors(j) - S)[0]
        return False

    def N6_ab_find_candidate(self):
        r"""
        Return some vertex whose component is (equivalent to) one
        of the first two graphs under N6 on p. 312 of [Neu1981]_, otherwise
        return ``-1``.
        
        OUTPUT:
        
        A vertex (nonnegative integer) or ``-1``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_cycle(5/4)
            3
            sage: P.N6_ab_find_candidate()
            0
        """
        
        for i in self.vertices:
            if self.N6_ab_candidate(i):
                return i
        
        return -1


    def N6_c_candidate(self, i):
        r"""
        Check if the vertex ``i`` is one of the endpoints of a string seen as
        the third graph under N6 on p. 312 of [Neu1981]_, i.e. a string whose
        ends have ``mb =-1``, ``g=-1``, ``r=0`` and in between vertices with
        ``mb=-2``, ``g=0``, ``r=0``.
        
        INPUT:
        
        - ``i`` -- a vertex
        
        OUTPUT:
        
        ``True`` or ``False``.

        EXAMPLES::
        
            sage: P = PlumbingGraph()
            sage: P.add_seifert(-1,-1,[4/3])
            0
            sage: P.add_vertex(-1,-1,0)
            4
            sage: P.add_edge({3,4})
            3
            sage: P.N6_c_candidate(0)
            True
            sage: P.N6_c_candidate(1)
            False
        """
        
        if self.mb[i] != -1 or self.g[i] != -1 or self.r[i] != 0 or self.degree(i) != 1:
            return False
        
        j = self.neighbor(i)
        S = {i}
        while True:
            if self.mb[j] == -1 and self.g[j] == -1 and self.r[j] == 0 and self.degree(j) == 1:
                return True
            elif self.mb[j] == -2 and self.g[j] == 0 and self.r[j] == 0 and self.degree(j) == 2:
                S |= {j}
                j = list(self.neighbors(j) - S)[0]
                continue
            else:
                return False

    def N6_c_find_candidate(self):
        r"""
        Return some vertex which is one of the endpoints of the third graph
        under N6 on p. 312 of [Neu1981]_, otherwise return ``-1``.
        
        OUTPUT:
        
        A vertex (nonnegative integer) or ``-1``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-1,-1,[4/3])
            0
            sage: P.add_vertex(-1,-1,0)
            4
            sage: P.add_edge({3,4})
            3
            sage: P.N6_c_find_candidate()
            0
            sage: P.add_edge({4})
            4
            sage: P.N6_c_find_candidate()
            -1
        """
        
        for i in self.vertices:
            if self.N6_c_candidate(i):
                return i
        
        return -1


    def N6_d_candidate(self, i):
        r"""
        Check if the vertex ``i`` is one of the endpoints of a string seen as
        the fourth graph under N6 on p. 312 of [Neu1981]_, i.e. a string which
        consists of two vertices with ``mb=1``, ``g=-1``, ``r=0``.
        
        INPUT:
        
        - ``i`` -- a vertex
        
        OUTPUT:
        
        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(1,-1,0)
            0
            sage: P.add_vertex(1,-1,0)
            1
            sage: P.N6_d_candidate(0)
            False
            sage: P.add_edge({0,1})
            0
            sage: P.N6_d_candidate(0)
            True
        """
        
        if self.mb[i] != 1 or self.g[i] != -1 or self.r[i] != 0 or self.degree(i) != 1:
            return False
        j = self.neighbor(i)
        if self.mb[j] != 1 or self.g[j] != -1 or self.r[j] != 0 or self.degree(j) != 1:
            return False
        return True

    def N6_d_find_candidate(self):
        r"""
        Return some vertex which is one of the endpoints of the fourth graph
        under N6 on p. 312 of [Neu1981]_, otherwise return ``-1``.
        
        OUTPUT:
        
        A vertex (nonnegative integer) or ``-1``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(1,-1,0)
            0
            sage: P.add_vertex(1,-1,0)
            1
            sage: P.N6_d_find_candidate()
            -1
            sage: P.add_edge({0,1})
            0
            sage: P.N6_d_find_candidate()
            0
        """
        
        for i in self.vertices:
            if self.N6_d_candidate(i):
                return i
        
        return -1


    def N6_e_candidate(self, i):
        r"""
        Check if the vertex ``i`` is one of the endpoints of a string seen as
        the fifth graph under N6 on p. 312 of [Neu1981]_, i.e. a string which
        consists of three vertices, the endpoints having ``mb=1``,
        ``g=-1``, ``r=0``,
        and the middle one ``mb=e`` (some integer), ``g=0``, ``r=0``.
        
        INPUT:
        
        - ``i`` -- a vertex
        
        OUTPUT:
        
        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(0,-1,0)
            0
            sage: P.add_vertex(10,0,0)
            1
            sage: P.add_vertex(0,-1,0)
            2
            sage: P.add_edge({0,1})
            0
            sage: P.N6_e_candidate(0)
            False
            sage: P.add_edge({1,2})
            1
            sage: P.N6_e_candidate(0)
            True
        """
        if self.mb[i] == 0 and self.g[i] == -1 and self.r[i] == 0 and self.degree(i) == 1:
            j = list(self.neighbors(i))[0]
            if self.g[j] == 0 and self.r[j] == 0 and self.degree(j) == 2:
                k = list(self.neighbors(j) - {i})[0]
                if self.mb[k] == 0 and self.g[k] == -1 and self.r[k] == 0 and self.degree(k) == 1:
                    return True
        return False

    def N6_e_find_candidate(self):
        r"""
        Return some vertex which is one of the endpoints of the fifth graph
        under N6 on p. 312 of [Neu1981]_, otherwise return ``-1``.
        
        OUTPUT:
        
        A vertex (nonnegative integer) or ``-1``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(0,-1,0)
            0
            sage: P.add_vertex(10,0,0)
            1
            sage: P.add_vertex(0,-1,0)
            2
            sage: P.add_edge({0,1})
            0
            sage: P.N6_e_find_candidate()
            -1
            sage: P.add_edge({1,2})
            1
            sage: P.N6_e_find_candidate()
            0
        """
        
        for i in self.vertices:
            if self.N6_e_candidate(i):
                return i
        
        return -1


    def N6_f_candidate(self, i):
        r"""
        Check if the vertex ``i`` is the second-last graph on p. 312 of
        [Neu1981]_, i.e. a singleton with ``mb=0``, ``g=-2``, ``r=0``.
        
        INPUT:
        
        - ``i`` -- a vertex
        
        OUTPUT:
        
        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(0,-2,0)
            0
            sage: P.add_vertex(1,-2,0)
            1
            sage: P.add_vertex(0,-3,0)
            2
            sage: P.add_vertex(0,-2,3)
            3
            sage: P.N6_f_candidate(0)
            True
            sage: P.N6_f_candidate(1)
            False
            sage: P.N6_f_candidate(2)
            False
            sage: P.N6_f_candidate(3)
            False
        """
        return self.mb[i] == 0 and self.g[i] == -2 and self.r[i] == 0 and self.degree(i) == 0

    def N6_f_find_candidate(self):
        r"""
        Return some vertex which is forms the second-last graph on
        p. 312 of [Neu1981]_, otherwise returns ``-1``.
        
        OUTPUT:
        
        A vertex (nonnegative integer) or ``-1``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(1,-2,0)
            0
            sage: P.add_vertex(0,-1,0)
            1
            sage: P.add_vertex(0,-2,1)
            2
            sage: P.N6_f_find_candidate()
            -1
            sage: P.add_vertex(0,-2,0)
            3
            sage: P.N6_f_find_candidate()
            3
        """
        
        for i in self.vertices:
            if self.N6_f_candidate(i):
                return i
        
        return -1


    def N6_g_candidate(self, i):
        r"""
        Check if the vertex ``i`` is the last graph on p. 312 of [Neu1981]_,
        i.e. a singleton with ``g=-1``, ``r=1``.
        
        INPUT:
        
        - ``i`` -- a vertex
        
        OUTPUT:
        
        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(4,-1,3)
            0
            sage: P.N6_g_candidate(0)
            False
            sage: P.add_vertex(4,-1,1)
            1
            sage: P.N6_g_candidate(1)
            True
        """
        return self.g[i] == -1 and self.r[i] == 1 and self.degree(i) == 0

    def N6_g_find_candidate(self):
        r"""
        Return some vertex which is forms the last graph on
        p. 312 of [Neu1981]_, otherwise return ``-1``.
        
        OUTPUT:
        
        A vertex (nonnegative integer) or ``-1``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(4,-1,2)
            0
            sage: P.N6_g_find_candidate()
            -1
            sage: P.add_vertex(4,-1,1)
            1
            sage: P.N6_g_find_candidate()
            1
        """
        
        for i in self.vertices:
            if self.N6_g_candidate(i):
                return i
        
        return -1

########################################################################
# Finding minimal graph and normal form

    def make_minimal_good_resolution(self):
        r"""
        Reduce a good resolution graph to the minimal good representative.
        
        I.e. change the graph to its minimal representative, assuming
        that the graph is a resolution graph for a good resolution. This
        is obtained simply by repeatedly blowing down any
        (-1)-curves which are not nodes and have ``g = r = 0``.
        The result coincides with the minimal good
        resolution graph of a corresponding singularity.

        This is not the same as applying the algorithm below for the
        normal form, as this can introduce vertices with
        negative genus, which is something which does not occur for
        resolution graphs.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(5/4)
            3
            sage: P.add_vertex(-1,0,0)
            4
            sage: P.add_edge({3,4})
            3
            sage: P.make_minimal_good_resolution()
            sage: P
            An empty plumbing graph
        """
        if not self.is_good_resolution():
            raise TypeError("this function only works on resolution graphs")
        while True:
            v = self.find_R1_candidate()
            if v == -1:
                return
            else:
                self.blow_down(v)

    def make_normal(self):
        r"""
        Perform steps 1-6 described in [Neu1981]_.

        EXAMPLES:

        A pretty basic example::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(9/7)
            3
            sage: P.mb
            {0: -2, 1: -2, 2: -2, 3: -3}
            sage: P.add_vertex(-1,0,0)
            4
            sage: P.add_edge({4,0})
            3
            sage: P.make_normal()
            sage: P
            A plumbing graph with one vertex

        In the following example, we construct a plumbed manifold which
        is the link of a surface singularity, but whose minimal plumbing
        graph is not the resolution graph of a resolution::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-2,0,[2,2,3/2])
            0
            sage: P.add_seifert(-2,0,[2,2,3/2])
            5
            sage: P.add_vertex(-3,0,0)
            10
            sage: P.add_edge({10,4})
            8
            sage: P.add_edge({10,9})
            9
            sage: P.is_good_resolution()
            True
            sage: P.make_normal()
            sage: P.is_good_resolution()
            False
        """
        self.make_normal_step_1()
        self.make_normal_step_2()
        self.make_normal_step_3()
        self.make_normal_step_4()
        self.make_normal_step_5()
        self.make_normal_step_6()

########################################################################
# Step one

    def make_normal_step_1(self):
        r"""
        Execute step 1 in the proof of thm 4.1 in [Neu1981]_.

        This just performs moves R1-8 until this is no longer possible.

        OUTPUT:

        ``True`` if ``N1()`` holds in the end, otherwise ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(9/7)
            3
            sage: P.add_vertex(-1,0,0)
            4
            sage: P.add_edge({4,0})
            3
            sage: P
            A plumbing graph with 5 vertices
            sage: P.make_normal_step_1()
            True
            sage: P
            A plumbing graph with one vertex
        """
        while True:
            i = self.find_R1_candidate()
            if i != -1:
                self.blow_down(i)
                continue

            j = self.find_R2_candidate()
            if j != -1:
                i = list(self.R2_candidate_neighbors(j))[0]
                self.projective_absorb(j,i)
                continue

            k = self.find_R3_candidate()
            if k != -1:
                self.zero_chain_absorb(k)
                continue

            j = self.find_R4_candidate()
            if j != -1:
                self.unoriented_handle_absorb(j)
                continue

            j = self.find_R5_candidate()
            if j != -1:
                self.oriented_handle_absorb(j)
                continue

            j = self.find_R6_candidate()
            if j != -1:
                self.split(j)
                continue

            i = self.find_R7_candidate()
            if i != -1:
                self.seifert_graph_exchange(i)
                continue

            j = self.find_R8_candidate()
            if j != -1:
                self.annulus_absorb(j)
                continue
            break

        return self.N1()
            
########################################################################
# Step two

    def make_normal_step_2(self):
        r"""
        Execute step 2 in the proof of thm 4.1 in [Neu1981]_.

        This modifies all chains to have only Euler number ``<= -2``.

        OUTPUT:

        ``True`` if ``N2()`` holds in the end, otherwise ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_bamboo(9/7)
            3
            sage: for i in P.vertices:
            ....:     P.mb[i] = -P.mb[i]
            ....:     
            sage: P
            A plumbing graph with 4 vertices
            sage: P.make_normal_step_2()
            sage: P
            A plumbing graph with 2 vertices
            sage: P.mb
            {0: -2, 4: -5}
        """
        while True:
            L = [i for i in self.vertices - self.nodes() if self.mb[i] > -2]
            if not L:
                return
            i = L[0]
            mc = self.max_chain(i)
            plus_twos = {j for j in mc if self.mb[j] == 2}
            if self.max_chain_is_cycle(i) and mc == plus_twos:
                s = len(mc)
                cycle_edges = [e for e in self.edges if self.adj[e] <= mc]
                cycle_sign = 1
                for e in cycle_edges:
                    cycle_sign *= self.epsilon[e]
                for j in mc:
                    self.delete_vertex(j)
                if cycle_sign == 1:
                    self.add_vertex(-s, 1, 0)
                elif cycle_sign == -1:
                    self.add_vertex(-s, -2, 0)
            else:
                e = list(self.adjacent_edges(i))[0]
                for a in range(self.mb[i]-1):
                    new_vx = self.blow_up_edge(e, -1)
                    e = [ ed for ed in self.edges
                        if self.adj[ed] == {new_vx,i} ][0]
                self.blow_down(i)
                self.make_normal_step_1()
        return self.N2()

########################################################################
# Step three
    def make_normal_step_3(self):
        r"""
        Execute step 3 in the proof of thm 4.1 in [Neu1981]_.

        OUTPUT:

        ``True`` if ``N3()`` holds in the end, otherwise ``False``.

        EXAMPLES:

        In the following example, we start by creating a Seifert graph
        with two (-2)-neighbors which are leaves, like in N3. For this to
        really violate N3, we add another edge on the longer bamboo.
        This graph satisfies N1 and N2, and after applying
        ``make_normal_step_3()``, it also satisfies N3::
        
            sage: P = PlumbingGraph()
            sage: P.add_seifert(5, 0, [2,2,7/5])
            0
            sage: P.add_edge({5},1)
            5
            sage: P.N1()
            True
            sage: P.N2()
            True
            sage: P.N3()
            False
            sage: P.make_normal_step_3()
            True
            sage: P.N3()
            True
        """
        if (not self.N1()) or (not self.N2()):
            return False
        while True:
            L = [j for j in self.vertices if self.N3_obstruction(j)]
            if not L:
                return True
            j = L[0]
            i = list(self.neighbors(j) - self._minus_two_leaves(j))[0]
            M = self.zero_chain_extrude(
                j,
                [self.mb[j]+1, -1],
                [0,0],
                [0,0],
                [
                    {e for e in self.adjacent_edges(j) if i in self.adj[e]},
                    {e for e in self.adjacent_edges(j) if
                        len(self._minus_two_leaves(j) & self.adj[e]) != 0}
                ])
            self.projective_absorb(M[2], M[1])
            # the above can break N1 or N2, so let's throw this thing in:
            self.make_normal_step_1()
            self.make_normal_step_2()
        
        return self.N3()

########################################################################
# Step four
    def make_normal_step_4(self):
        r"""
        Carry out step 4 in [Neu1981]_.

        OUTPUT:
        
        ``True`` or ``False``.

        EXAMPLES::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(5,-3,1)
            0
            sage: P.add_vertex(0,-1,0)
            1
            sage: P.add_edge({0,1})
            0
            sage: P.N3()
            True
            sage: P.N4()
            False
            sage: P.make_normal_step_4()
            True
            sage: P.N4()
            True
        """
        if not self.N3():
            return False
        
        while True:
            L = [j for j in self.vertices if self.N4_obstruction(j)]
            if not L:
                break
            j = L[0]
            self.projective_extrude(j, 0, -1, -1, -1)
            self.zero_chain_absorb(j)
        
        return self.N4()

########################################################################
# Step five
    def make_normal_step_5(self):
        r"""
        Carry out step 5 in [Neu1981]_.

        OUTPUT:
        
        ``True`` or ``False``.

        EXAMPLES:

        We test all the examples on p. 315 of [Neu1981]_.

        The first example, only one vertex::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(-3,-1,0)
            0
            sage: P.N4()
            True
            sage: P.N5()
            False
            sage: P.make_normal_step_5()
            True
            sage: P.N5()
            True
        
        The first example, several vertices::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-3,-1,[5/4])
            0
            sage: P.N4()
            True
            sage: P.N5()
            False
            sage: P.make_normal_step_5()
            True
            sage: P.N5()
            True

        The second example (b = 4)::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-1,-1,[11/9])
            0
            sage: P.make_normal_step_5()
            True
            sage: P.N5()
            True
            sage: P
            A plumbing graph with 4 vertices

        The third example, case 1 (b=4)::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-1, -1, [5/4])
            0
            sage: P.N4()
            True
            sage: P.N5()
            False
            sage: P.make_normal_step_5()
            True
            sage: P
            A plumbing graph with 5 vertices
        
        The third example, case 2 (b=0)::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(-1,-1,0)
            0
            sage: P.N5()
            False
            sage: P.make_normal_step_5()
            True
            sage: P
            A plumbing graph with one vertex

        The fourth example::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(0,-1,[7/4])
            0
            sage: P.N5()
            False
            sage: P.make_normal_step_5()
            True

        The fourth example again::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(1,-1,[7/4])
            0
            sage: P.N5()
            False
            sage: P.make_normal_step_5()
            True
            sage: P
            A plumbing graph with 5 vertices

        The fifth example, e=2 (the end result is D_4)::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(2,-1,0)
            0
            sage: P.make_normal_step_5()
            True

        The fifth example, e=0::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(0,-1,0)
            0
            sage: P.make_normal_step_5()
            True
            sage: P
            A plumbing graph with 2 vertices
        """
        L = [ j for j in self.vertices if self.N5_obstruction(j) ]
        for j in L:
            self.projective_extrude(j, 0, -1, -1, -1)

        # blow down as much as possible:
        while True:
            i = self.find_R1_candidate()
            if i != -1:
                self.blow_down(i)
                continue
            # in the case of a single vertex with mb = g = -1, we need to
            # zero-chain absorb:
            i = self.find_R3_candidate()
            if i != -1:
                self.zero_chain_absorb(i)
                continue
            # in the last case, we have to split:
            i = self.find_R6_candidate()
            if i != -1:
                self.split(i)
                continue
            break

        self.make_normal_step_2()
        
        return self.N5()

########################################################################
# Step six
    def make_normal_step_6(self):
        r"""
        Carry out step 6 in [Neu1981]_.

        OUTPUT:
        
        ``True`` or ``False``.
        
        EXAMPLES:

        A cycle with positive edges::

            sage: P = PlumbingGraph()
            sage: P.add_cycle(5/4)
            3
            sage: P.make_normal_step_6()
            True
            sage: P
            A plumbing graph with one vertex

        A cycle with a negative edge::

            sage: P = PlumbingGraph()
            sage: P.add_cycle(5/4, epsilon=-1)
            3
            sage: P.make_normal_step_6()
            True
            sage: P
            A plumbing graph with one vertex

        The third graph::

            sage: P = PlumbingGraph()
            sage: P.add_seifert(-1,-1,[6/5])
            0
            sage: P.add_vertex(-1,-1,0)
            6
            sage: P.vertices
            {0, 1, 2, 3, 4, 5, 6}
            sage: P.add_edge({5,6})
            5
            sage: P.make_normal_step_6()
            True
            sage: P
            A plumbing graph with 5 vertices

        The fourth graph::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(1,-1,0)
            0
            sage: P.add_vertex(1,-1,0)
            1
            sage: P.add_edge({0,1})
            0
            sage: P.make_normal_step_6()
            True
            sage: P
            A plumbing graph with 5 vertices

        The fifth graph::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(0,-1,0)
            0
            sage: P.add_vertex(4,0,0)
            1
            sage: P.add_vertex(0,-1,0)
            2
            sage: P.add_edge({0,1})
            0
            sage: P.add_edge({1,2})
            1
            sage: P.make_normal_step_6()
            True
            sage: P
            A plumbing graph with 5 vertices

        The sixth graph::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(0,-2,0)
            0
            sage: P.make_normal_step_6()
            True
            sage: P
            A plumbing graph with 5 vertices

        The seventh graph::

            sage: P = PlumbingGraph()
            sage: P.add_vertex(0,-1,1)
            0
            sage: P.make_normal_step_6()
            True
            sage: P
            A plumbing graph with 3 vertices
        """
        # We use here the naming convention abcdefg for the graphs on
        # p. 316 [Neu1981]_ established above for N6.
        
        # ab
        while True:
            i = self.N6_ab_find_candidate()
            if i == -1:
                break
            C = self.component(i)
            e = len(C)
            s = len({e for e in self.edges if self.adj[e] <= C and self.epsilon[e] == -1})
            self.delete_component(i)
            if not s%2:
                self.add_vertex(e, 1, 0)
            else:
                self.add_vertex(e, -2, 0)
        
        # c
        while True:
            i = self.N6_c_find_candidate()
            if i == -1:
                break
            C = self.component(i)
            e = len(C) - 2
            self.delete_component(i)
            self.add_seifert(e-1, 0, [2,2,2,2])

        # d
        while True:
            i = self.N6_d_find_candidate()
            if i == -1:
                break
            self.delete_component(i)
            self.add_seifert(-3, 0, [2,2,2,2])

        # e
        while True:
            i = self.N6_e_find_candidate()
            if i == -1:
                break
            j = self.neighbor(i)
            e = self.mb[j]
            self.delete_component(i)
            self.add_seifert(e-2, 0, [2,2,2,2])

        # f
        while True:
            i = self.N6_f_find_candidate()
            if i == -1:
                break
            self.delete_component(i)
            self.add_seifert(-2, 0, [2,2,2,2])

        # g
        while True:
            i = self.N6_g_find_candidate()
            if i == -1:
                break
            self.delete_component(i)
            self.add_seifert(0, 0, [2,2], r=1)
        return self.N6()
