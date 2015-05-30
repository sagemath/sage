from sage.structure.sage_object cimport SageObject

cdef class Matroid(SageObject):
    cdef public __custom_name
    cdef public _custom_name
    cdef public _cached_info
    cdef int _stored_full_rank
    cdef int _stored_size

    # virtual methods
    cpdef groundset(self)
    cpdef _rank(self, X)

    # internal methods, assuming verified input
    cpdef _max_independent(self, X)
    cpdef _circuit(self, X)
    cpdef _fundamental_circuit(self, B, e)
    cpdef _closure(self, X)
    cpdef _corank(self, X)
    cpdef _max_coindependent(self, X)
    cpdef _cocircuit(self, X)
    cpdef _fundamental_cocircuit(self, B, e)
    cpdef _coclosure(self, X)
    cpdef _augment(self, X, Y)

    cpdef _is_independent(self, X)
    cpdef _is_basis(self, X)
    cpdef _is_circuit(self, X)
    cpdef _is_closed(self, X)
    cpdef _is_coindependent(self, X)
    cpdef _is_cobasis(self, X)
    cpdef _is_cocircuit(self, X)
    cpdef _is_coclosed(self, X)

    cpdef _minor(self, contractions, deletions)
    cpdef _has_minor(self, N)
    cpdef _line_length(self, F)
    cpdef _extension(self, element, hyperplanes)

    # ** user-facing methods **

    # cpdef _latex_(self)  # Disabled, because not overridden by current subclasses
    # cpdef show(self)  # Disabled, because not implemented yet
    cpdef size(self)

    # matroid oracle
    cpdef rank(self, X=*)
    cpdef full_rank(self)
    cpdef basis(self)
    cpdef max_independent(self, X)
    cpdef circuit(self, X=*)
    cpdef fundamental_circuit(self, B, e)
    cpdef closure(self, X)
    cpdef k_closure(self, X, k)

    cpdef augment(self, X, Y=*)

    cpdef corank(self, X=*)
    cpdef full_corank(self)
    cpdef cobasis(self)
    cpdef max_coindependent(self, X)
    cpdef cocircuit(self, X=*)
    cpdef fundamental_cocircuit(self, B, e)
    cpdef coclosure(self, X)

    cpdef loops(self)
    cpdef is_independent(self, X)
    cpdef is_dependent(self, X)
    cpdef is_basis(self, X)
    cpdef is_circuit(self, X)
    cpdef is_closed(self, X)
    cpdef is_subset_k_closed(self, X, k)

    cpdef coloops(self)
    cpdef is_coindependent(self, X)
    cpdef is_codependent(self, X)
    cpdef is_cobasis(self, X)
    cpdef is_cocircuit(self, X)
    cpdef is_coclosed(self, X)

    # verification
    cpdef is_valid(self)

    # enumeration
    cpdef circuits(self)
    cpdef nonspanning_circuits(self)
    cpdef cocircuits(self)
    cpdef noncospanning_cocircuits(self)
    cpdef circuit_closures(self)
    cpdef nonspanning_circuit_closures(self)
    cpdef bases(self)
    cpdef independent_r_sets(self, long r)
    cpdef nonbases(self)
    cpdef dependent_r_sets(self, long r)
    cpdef _extend_flags(self, flags)
    cpdef _flags(self, r)
    cpdef flats(self, r)
    cpdef coflats(self, r)
    cpdef hyperplanes(self)
    cpdef f_vector(self)
    cpdef broken_circuits(self, ordering=*)
    cpdef no_broken_circuits_sets(self, ordering=*)

    # isomorphism
    cpdef is_isomorphic(self, other)
    cpdef _is_isomorphic(self, other)
    cpdef isomorphism(self, other)
    cpdef _isomorphism(self, other)
    cpdef equals(self, other)
    cpdef is_isomorphism(self, other, morphism)
    cpdef _is_isomorphism(self, other, morphism)

    # minors, dual, trucation
    cpdef minor(self, contractions=*, deletions=*)
    cpdef contract(self, X)
    cpdef delete(self, X)
    cpdef _backslash_(self, X)
    cpdef dual(self)
    cpdef truncation(self)
    cpdef has_minor(self, N)
    cpdef has_line_minor(self, k, hyperlines=*)
    cpdef _has_line_minor(self, k, hyperlines)

    # extension
    cpdef extension(self, element=*, subsets=*)
    cpdef coextension(self, element=*, subsets=*)
    cpdef modular_cut(self, subsets)
    cpdef linear_subclasses(self, line_length=*, subsets=*)
    cpdef extensions(self, element=*, line_length=*, subsets=*)

    # connectivity
    cpdef simplify(self)
    cpdef cosimplify(self)
    cpdef is_simple(self)
    cpdef is_cosimple(self)
    cpdef components(self)
    cpdef is_connected(self)
    cpdef is_3connected(self)

    # matroid k-closed
    cpdef is_k_closed(self, k)

    # optimization
    cpdef max_weight_independent(self, X=*, weights=*)
    cpdef max_weight_coindependent(self, X=*, weights=*)
    cpdef intersection(self, other, weights=*)
    cpdef _intersection(self, other, weights)
    cpdef _intersection_augmentation(self, other, weights, Y)

    # invariants
    cpdef _internal(self, B)
    cpdef _external(self, B)
    cpdef tutte_polynomial(self, x=*, y=*)
    cpdef flat_cover(self)
    
    # visualization
    cpdef plot(self,B=*,lineorders=*,pos_method=*,pos_dict=*,save_pos=*)
    cpdef show(self,B=*,lineorders=*,pos_method=*,pos_dict=*,save_pos=*,lims=*)
    cpdef _fix_positions(self,pos_dict=*,lineorders=*)
