# sage_setup: distribution = sagemath-graphs
from libc.stdint cimport uint32_t

cdef unsigned short * c_shortest_path_all_pairs(G, vertex_list=*) except NULL
cdef unsigned short * c_distances_all_pairs(G, vertex_list=*) noexcept
cdef all_pairs_shortest_path_BFS(gg,
                                 unsigned short * predecessors,
                                 unsigned short * distances,
                                 uint32_t       * eccentricity,
                                 vertex_list=*)

cdef uint32_t * c_eccentricity(G, vertex_list=*) except NULL
