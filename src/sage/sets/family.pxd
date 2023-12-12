from sage.structure.parent cimport Parent


cdef class AbstractFamily(Parent):
    cdef public __custom_name


cdef class FiniteFamily(AbstractFamily):
    cdef public dict _dictionary
    cdef public object _keys
