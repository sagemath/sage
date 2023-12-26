cimport cython
from sage.data_structures.list_of_pairs cimport ListOfPairs
from sage.structure.sage_object         cimport SageObject
from sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator                     cimport FaceIterator, CombinatorialFace
from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces                     cimport ListOfFaces
from sage.geometry.polyhedron.combinatorial_polyhedron.face_data_structure               cimport face_t
from sage.geometry.polyhedron.combinatorial_polyhedron.polyhedron_face_lattice           cimport PolyhedronFaceLattice


@cython.final
cdef class CombinatorialPolyhedron(SageObject):
    cdef public dict _cached_methods

    # Do not assume any of those attributes to be initialized, use the corresponding methods instead.
    cdef tuple _Vrep                       # the names of VRep, if they exist
    cdef tuple _facet_names                # the names of HRep without equations, if they exist
    cdef tuple _equations                  # stores equations, given on input (might belong to Hrep)
    cdef int _dimension                    # stores dimension, -2 on init
    cdef unsigned int _n_Hrepresentation   # Hrep might include equations
    cdef unsigned int _n_Vrepresentation   # Vrep might include rays/lines
    cdef size_t _n_facets                  # length Hrep without equations
    cdef bint _bounded                     # ``True`` iff Polyhedron is bounded
    cdef ListOfFaces _bitrep_facets        # facets in bit representation
    cdef ListOfFaces _bitrep_Vrep          # vertices in bit representation
    cdef face_t _far_face                  # a 'face' containing all none-vertices of Vrep
    cdef tuple _far_face_tuple
    cdef tuple _f_vector

    cdef ListOfPairs _edges                    # stores edges labeled by vertex indices
    cdef ListOfPairs _ridges                   # stores ridges labeled by facet indices
    cdef ListOfPairs _face_lattice_incidences  # stores incidences in Hasse diagram labeled indices of the faces
    cdef PolyhedronFaceLattice _all_faces     # class to generate Hasse diagram incidences

    cdef tuple Vrep(self) noexcept
    cdef tuple facet_names(self) noexcept
    cdef tuple equations(self) noexcept
    cdef tuple equalities(self) noexcept
    cdef unsigned int n_Vrepresentation(self) noexcept
    cdef unsigned int n_Hrepresentation(self) noexcept
    cdef bint is_bounded(self) noexcept
    cdef ListOfFaces bitrep_facets(self) noexcept
    cdef ListOfFaces bitrep_Vrep(self) noexcept
    cdef tuple far_face_tuple(self) noexcept
    cdef int _algorithm_to_dual(self, algorithm) except -2

    # Methods to initialize the combinatorial polyhedron.
    cdef _init_from_polyhedron(self, data) noexcept
    cdef _init_from_lattice_polytope(self, data) noexcept
    cdef _init_from_cone(self, data) noexcept
    cdef _init_facet_names(self, facets) noexcept
    cdef _init_from_incidence_matrix(self, data) noexcept
    cdef _init_from_list_of_facets(self, data) noexcept
    cdef _init_from_ListOfFaces(self, ListOfFaces facets, ListOfFaces Vrep) noexcept
    cdef _initialize_far_face(self) noexcept
    cdef _init_as_trivial_polyhedron(self, int dimension) noexcept

    # Methods to obtain a different combinatorial polyhedron.
    cpdef CombinatorialPolyhedron dual(self) noexcept
    cpdef CombinatorialPolyhedron pyramid(self, new_vertex=*, new_facet=*) noexcept

    cdef FaceIterator _face_iter(self, bint dual, int dimension) noexcept
    cdef int _compute_f_vector(self, size_t num_threads, size_t parallelization_depth, int dual) except -1
    cdef int _persist_f_vector(self, size_t* input_f_vector, bint input_is_reversed) except -1

    cdef inline int _compute_edges(self, dual) except -1:
        return self._compute_edges_or_ridges(dual, True)

    cdef inline int _compute_ridges(self, dual) except -1:
        return self._compute_edges_or_ridges(dual, False)

    cdef int _compute_edges_or_ridges(self, int dual, bint do_edges) except -1
    cdef size_t _compute_edges_or_ridges_with_iterator(
            self, FaceIterator face_iter, const bint do_atom_rep,
            ListOfPairs edges, size_t* f_vector) except -1

    cdef int _compute_face_lattice_incidences(self) except -1
