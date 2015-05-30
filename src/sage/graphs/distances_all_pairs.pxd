cdef unsigned short * c_shortest_path_all_pairs(G) except NULL
cdef unsigned short * c_distances_all_pairs(G)
cdef all_pairs_shortest_path_BFS(gg,
                                 unsigned short * predecessors,
                                 unsigned short * distances,
                                 int            * eccentricity)

cdef int * c_eccentricity(G) except NULL
