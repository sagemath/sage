# sage_setup: distribution = sagemath-polyhedra
cimport cython
from sage.structure.sage_object cimport SageObject
from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces             cimport ListOfFaces
from sage.geometry.polyhedron.combinatorial_polyhedron.face_data_structure       cimport face_t
from sage.geometry.polyhedron.combinatorial_polyhedron.face_list_data_structure  cimport face_list_t
from sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face        cimport CombinatorialFace

cdef enum FaceStatus:
    NOT_INITIALIZED
    INITIALIZED
    IGNORE_SUBSETS
    ONLY_VISIT_SUBSETS

cdef struct iter_s:
    bint dual                  # if 1, then iterate over dual Polyhedron
    face_t face                # the current face of the iterator
    FaceStatus face_status
    size_t *atom_rep           # a place where atom-representaion of face will be stored
    size_t *coatom_rep         # a place where coatom-representaion of face will be stored
    int current_dimension      # dimension of current face, dual dimension if ``dual``
    int dimension              # dimension of the polyhedron
    int output_dimension       # only faces of this (dual?) dimension are considered
    int lowest_dimension       # don't consider faces below this (dual?) dimension
    int highest_dimension      # don't consider faces above this (dual?) dimension
    size_t _index              # this counts the number of seen faces, useful for hasing the faces

    # ``visited_all`` points to faces, of which we have visited all faces already.
    # The number of faces in ``visited_all`` might depend on the current dimension:
    #     Consider we visit the facets A,B of some face F.
    #     We will first visit all faces of A and then add A to visited_all.
    #     Then we visit all faces of B and add B to visited_all.
    #     Then we have visited F completely.
    #     Instead of having A and B in ``visited_all`` we will point to F.
    #     In this way, we will append ``visited_all`` in lower dimension, but
    #     will ignore those changes when going up in dimension again.
    #     This is why the number of faces in ``visited_all``depends on dimension.
    face_list_t* visited_all

    # ``new_faces`` is where the new faces are stored.
    # Needs to be long enough to store all possible intersections of a face with all coatoms.
    face_list_t* new_faces

    # After having visited a face completely, we want to add it to ``visited_all``.
    # ``first_time[i]`` will indicate, whether there is one more face in
    # ``newfaces[i]`` then ``n_newfaces[i]`` suggests
    # that has to be added to ``visited_all``.
    # If ``first_time[i] == False``, we still need to
    # add ``newfaces[i][n_newfaces[i]]`` to ``visited_all``.
    bint *first_time

    # The number of elements in newfaces[current_dimension],
    # that have not been visited yet.
    size_t yet_to_visit
    size_t n_coatoms

ctypedef iter_s iter_t[1]


cdef class FaceIterator_base(SageObject):
    cdef iter_t structure
    cdef readonly bint dual         # if 1, then iterate over dual Polyhedron

    # some copies from ``CombinatorialPolyhedron``
    cdef tuple _Vrep, _facet_names, _equations
    cdef size_t _n_equations, _n_facets
    cdef bint _bounded
    cdef face_t _far_face

    # Atoms and coatoms are the vertices/facets of the Polyedron.
    # If ``dual == 0``, then coatoms are facets, atoms vertices and vice versa.
    cdef ListOfFaces atoms, coatoms, coatoms_coatom_rep

    cdef inline CombinatorialFace next_face(self)
    cdef inline int next_dimension(self) except -1
    cdef inline int next_face_loop(self) except -1
    cdef size_t n_atom_rep(self) except -1
    cdef size_t set_coatom_rep(self) except -1
    cdef size_t set_atom_rep(self) except -1
    cdef int ignore_subsets(self) except -1
    cdef int only_subsets(self) except -1
    cdef int find_face(self, face_t face) except -1


@cython.final
cdef class FaceIterator(FaceIterator_base):
    pass


@cython.final
cdef class FaceIterator_geom(FaceIterator_base):
    cdef int _trivial_faces     # Whether to yield the trivial faces.
    cdef object _requested_dim  # Dimension requested on init.
    cdef readonly object P      # The original polyhedron.

cdef int parallel_f_vector(iter_t* structures, size_t num_threads, size_t parallelization_depth, size_t *f_vector) except -1

# Nogil definitions of crucial functions.

cdef int next_dimension(iter_t structure, size_t parallelization_depth=?) except -1 nogil
cdef int next_face_loop(iter_t structure) except -1 nogil
cdef size_t n_atom_rep(iter_t structure) except -1 nogil
