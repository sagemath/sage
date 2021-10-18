"""
Abstract base classes for rings
"""

cdef class RealField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.real_mpfr.RealField_class`.
    """

    pass


cdef class RealDoubleField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.real_double.RealDoubleField_class`.
    """

    pass


cdef class ComplexField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.complex_mpfr.ComplexField_class`.
    """

    pass


cdef class ComplexDoubleField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.complex_double.ComplexDoubleField_class`.
    """

    pass
