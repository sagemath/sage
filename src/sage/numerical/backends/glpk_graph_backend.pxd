#*****************************************************************************
#       Copyright (C) 2012 Christian Kuper <christian.kuper@t-online.de>
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.glpk.types cimport glp_graph

ctypedef struct c_v_data:
         double rhs
         double pi
         double es
         double ls
         long cut

ctypedef struct c_a_data:
         double low
         double cap
         double cost
         double x


cdef class GLPKGraphBackend():
    cdef glp_graph * graph
    cpdef add_vertex(self, name = ?) noexcept
    cpdef list add_vertices(self, vertices) noexcept
    cpdef __add_vertices_sage(self, g) noexcept
    cpdef dict get_vertex(self, vertex) noexcept
    cpdef dict get_vertices(self, verts) noexcept
    cpdef set_vertex_demand(self, vertex, param) noexcept
    cpdef set_vertices_demand(self, list pairs) noexcept
    cpdef list vertices(self) noexcept
    cpdef add_edge(self, u, v, dict params = ?) noexcept
    cpdef __add_edges_sage(self, g) noexcept
    cpdef list add_edges(self, edges) noexcept
    cpdef delete_edge(self, u, v, dict params = ?) noexcept
    cpdef tuple get_edge(self, u, v) noexcept
    cpdef list edges(self) noexcept
    cpdef delete_vertex(self, vert) noexcept
    cpdef delete_vertices(self, list verts) noexcept
    cpdef int _find_vertex(self, vert) noexcept
    cpdef int write_graph(self, fname) noexcept
    cpdef int write_ccdata(self, fname) noexcept
    cpdef int write_mincost(self, fname) noexcept
    cpdef double mincost_okalg(self) except -1
    cdef int s
    cdef int t
    cpdef int write_maxflow(self, fname) except -1
    cpdef double maxflow_ffalg(self, u = ?, v = ?) except -1
    cpdef double cpp(self) noexcept
