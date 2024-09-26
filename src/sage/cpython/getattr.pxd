# sage_setup: distribution = sagemath-objects
# A class that allows for a more efficient creation
# of attribute errors, so that raising them requires
# less time.
cdef class AttributeErrorMessage:
    cdef public cls
    cdef public name

cpdef getattr_from_other_class(self, cls, name)
