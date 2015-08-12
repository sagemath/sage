from libc.stdint cimport uint32_t

cdef unsigned short * c_shortest_path_all_pairs(G) except NULL
cdef unsigned short * c_distances_all_pairs(G)
cdef all_pairs_shortest_path_BFS(gg,
                                 unsigned short * predecessors,
                                 unsigned short * distances,
                                 uint32_t       * eccentricity)

cdef uint32_t * c_eccentricity(G) except NULL
