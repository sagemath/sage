"""
Abstract base classes for rings
"""

cdef class RealField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.real_mpfr.RealField_class`.
    """

    pass


class RealBallField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.real_arb.RealBallField`.
    """

    pass


cdef class RealIntervalField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.real_mpfi.RealIntervalField_class`.
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


class ComplexBallField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.complex_arb.ComplexBallField`.
    """

    pass


class ComplexIntervalField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.complex_interval_field.ComplexIntervalField_class`.
    """

    pass


cdef class ComplexDoubleField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.complex_double.ComplexDoubleField_class`.
    """

    pass


class IntegerModRing:
    r"""
    Abstract base class for :class:`~sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic`.
    """

    pass

class Order:
    r"""
    Abstract base class for :class:`~sage.rings.number_field.order.Order`.
    """

    pass


cdef class SymbolicRing(CommutativeRing):
    r"""
    Abstract base class for :class:`~sage.rings.symbolic.ring.SymbolicRing`.
    """

    pass


class CallableSymbolicExpressionRing(SymbolicRing):
    r"""
    Abstract base class for :class:`~sage.rings.symbolic.callable.CallableSymbolicExpressionRing_class`.
    """

    pass
