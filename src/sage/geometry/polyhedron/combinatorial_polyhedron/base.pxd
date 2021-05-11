cimport cython
from sage.ext.memory_allocator  cimport MemoryAllocator
from sage.structure.sage_object cimport SageObject
from .face_iterator             cimport FaceIterator, CombinatorialFace
from .list_of_faces             cimport ListOfFaces
from .face_data_structure       cimport face_t
from .polyhedron_face_lattice   cimport PolyhedronFaceLattice

@cython.final
cdef class CombinatorialPolyhedron(SageObject):
    cdef public dict __cached_methods

    # Do not assume any of those attributes to be initialized, use the corresponding methods instead.
    cdef tuple _Vrep                       # the names of VRep, if they exist
    cdef tuple _facet_names                # the names of HRep without equalities, if they exist
    cdef tuple _equalities                 # stores equalities, given on input (might belong to Hrep)
    cdef int _dimension                    # stores dimension, -2 on init
    cdef unsigned int _n_Hrepresentation   # Hrep might include equalities
    cdef unsigned int _n_Vrepresentation   # Vrep might include rays/lines
    cdef size_t _n_facets                  # length Hrep without equalities
    cdef bint _bounded                     # ``True`` iff Polyhedron is bounded
    cdef ListOfFaces _bitrep_facets        # facets in bit representation
    cdef ListOfFaces _bitrep_Vrep          # vertices in bit representation
    cdef face_t _far_face                  # a 'face' containing all none-vertices of Vrep
    cdef tuple _far_face_tuple
    cdef tuple _f_vector

    # Edges, ridges and incidences are stored in a pointer of pointers.
    # The first edge has vertices ``edges[0][0]`` and ``edges[0][1]``,
    # the second edge has vertices ``edges[0][2]`` and ``edges[0][3]``, etc.
    # There are ``_length_edges_list`` edges in ``edges[i]``, so the edge
    # ``_length_edges_list + 1`` has vertices ``edges[1][0]`` and ``edges[1][1]``.
    # Likewise for ridges and incidences.
    cdef size_t _length_edges_list


    cdef size_t **_edges                    # stores edges labeled by vertex indices
    cdef size_t _n_edges
    cdef size_t **_ridges                   # stores ridges labeled by facet indices
    cdef size_t _n_ridges
    cdef size_t **_face_lattice_incidences  # stores incidences in Hasse diagram labeled indices of the faces
    cdef size_t _n_face_lattice_incidences
    cdef PolyhedronFaceLattice _all_faces   # class to generate Hasse diagram incidences

    cdef tuple Vrep(self)
    cdef tuple facet_names(self)
    cdef tuple equalities(self)
    cdef unsigned int n_Vrepresentation(self)
    cdef unsigned int n_Hrepresentation(self)
    cdef bint is_bounded(self)
    cdef ListOfFaces bitrep_facets(self)
    cdef ListOfFaces bitrep_Vrep(self)
    cdef tuple far_face_tuple(self)

    # Methods to obtain a different combinatorial polyhedron.
    cpdef CombinatorialPolyhedron dual(self)
    cpdef CombinatorialPolyhedron pyramid(self, new_vertex=*, new_facet=*)

    # Space for edges, ridges, etc. is allocated with ``MemoryAllocators``.
    # Upon success they are copied to ``_mem_tuple``.
    # Thus deallocation (at the correct time) is taken care of.
    cdef tuple _mem_tuple

    cdef FaceIterator _face_iter(self, bint dual, int dimension)
    cdef int _compute_f_vector(self, bint compute_edges=*, given_dual=*) except -1

    cdef inline int _compute_edges(self, dual) except -1:
        return self._compute_edges_or_ridges(dual, True)

    cdef inline int _compute_ridges(self, dual) except -1:
        return self._compute_edges_or_ridges(dual, False)

    cdef int _compute_edges_or_ridges(self, bint dual, bint do_edges) except -1
    cdef size_t _compute_edges_or_ridges_with_iterator(self, FaceIterator face_iter, bint do_atom_rep, size_t ***edges_pt, size_t *counter_pt, size_t *current_length_pt, MemoryAllocator mem) except -1
    cdef int _compute_face_lattice_incidences(self) except -1

    cdef inline int _set_edge(self, size_t a, size_t b, size_t ***edges_pt, size_t *counter_pt, size_t *current_length_pt, MemoryAllocator mem) except -1
    cdef inline size_t _get_edge(self, size_t **edges, size_t edge_number, size_t vertex) except -1
