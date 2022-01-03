"""
Abstract base classes for rings
"""
from sage.rings.ring import EuclideanDomain


class NumberField_quadratic(Field):
    r"""
    Abstract base class for :class:`~sage.rings.number_field.number_field.NumberField_quadratic`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: K.<sqrt2> = QuadraticField(2)                                     # optional - sage.rings.number_field
        sage: isinstance(K, sage.rings.abc.NumberField_quadratic)               # optional - sage.rings.number_field
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.NumberField_quadratic.__subclasses__()             # optional - sage.rings.number_field
        [<class 'sage.rings.number_field.number_field.NumberField_quadratic'>]

        sage: len(sage.rings.abc.NumberField_quadratic.__subclasses__()) <= 1
        True
    """

    pass


class NumberField_cyclotomic(Field):
    r"""
    Abstract base class for :class:`~sage.rings.number_field.number_field.NumberField_cyclotomic`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: K.<zeta> = CyclotomicField(15)                                    # optional - sage.rings.number_field
        sage: isinstance(K, sage.rings.abc.NumberField_cyclotomic)              # optional - sage.rings.number_field
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.NumberField_cyclotomic.__subclasses__()            # optional - sage.rings.number_field
        [<class 'sage.rings.number_field.number_field.NumberField_cyclotomic'>]

        sage: len(sage.rings.abc.NumberField_cyclotomic.__subclasses__()) <= 1
        True
    """

    pass


class AlgebraicField_common(Field):
    r"""
    Abstract base class for :class:`~sage.rings.qqbar.AlgebraicField_common`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(QQbar, sage.rings.abc.AlgebraicField_common)           # optional - sage.rings.number_field
        True
        sage: isinstance(AA, sage.rings.abc.AlgebraicField_common)              # optional - sage.rings.number_field
        True

    By design, other than the abstract subclasses :class:`~sage.rings.abc.AlgebraicField`
    and :class:`~sage.rings.abc.AlgebraicRealField`, there is only one direct implementation
    subclass::

        sage: sage.rings.abc.AlgebraicField_common.__subclasses__()             # optional - sage.rings.number_field
        [<class 'sage.rings.abc.AlgebraicField'>,
         <class 'sage.rings.abc.AlgebraicRealField'>,
         <class 'sage.rings.qqbar.AlgebraicField_common'>]

        sage: len(sage.rings.abc.AlgebraicField_common.__subclasses__()) <= 3
        True
    """

    pass


class AlgebraicField(AlgebraicField_common):
    r"""
    Abstract base class for :class:`~sage.rings.qqbar.AlgebraicField`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(QQbar, sage.rings.abc.AlgebraicField)                  # optional - sage.rings.number_field
        True
        sage: isinstance(AA, sage.rings.abc.AlgebraicField)                     # optional - sage.rings.number_field
        False

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.AlgebraicField.__subclasses__()                    # optional - sage.rings.number_field
        [<class 'sage.rings.qqbar.AlgebraicField'>]

        sage: len(sage.rings.abc.AlgebraicField.__subclasses__()) <= 1
        True
    """

    pass


class AlgebraicRealField(AlgebraicField_common):
    r"""
    Abstract base class for :class:`~sage.rings.qqbar.AlgebraicRealField`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(QQbar, sage.rings.abc.AlgebraicRealField)              # optional - sage.rings.number_field
        False
        sage: isinstance(AA, sage.rings.abc.AlgebraicRealField)                 # optional - sage.rings.number_field
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.AlgebraicRealField.__subclasses__()                # optional - sage.rings.number_field
        [<class 'sage.rings.qqbar.AlgebraicRealField'>]

        sage: len(sage.rings.abc.AlgebraicRealField.__subclasses__()) <= 1
        True
    """

    pass


cdef class RealField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.real_mpfr.RealField_class`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(RR, sage.rings.abc.RealField)
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.RealField.__subclasses__()
        [<class 'sage.rings.real_mpfr.RealField_class'>]

        sage: len(sage.rings.abc.RealField.__subclasses__()) <= 1
        True
    """

    pass


class RealBallField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.real_arb.RealBallField`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(RBF, sage.rings.abc.RealBallField)
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.RealBallField.__subclasses__()
        [<class 'sage.rings.real_arb.RealBallField'>]

        sage: len(sage.rings.abc.RealBallField.__subclasses__()) <= 1
        True
    """

    pass


cdef class RealIntervalField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.real_mpfi.RealIntervalField_class`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(RIF, sage.rings.abc.RealIntervalField)
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.RealIntervalField.__subclasses__()
        [<class 'sage.rings.real_mpfi.RealIntervalField_class'>]

        sage: len(sage.rings.abc.RealIntervalField.__subclasses__()) <= 1
        True
    """

    pass


cdef class RealDoubleField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.real_double.RealDoubleField_class`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(RDF, sage.rings.abc.RealDoubleField)
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.RealDoubleField.__subclasses__()
        [<class 'sage.rings.real_double.RealDoubleField_class'>]

        sage: len(sage.rings.abc.RealDoubleField.__subclasses__()) <= 1
        True
    """

    pass


cdef class ComplexField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.complex_mpfr.ComplexField_class`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(CC, sage.rings.abc.ComplexField)
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.ComplexField.__subclasses__()
        [<class 'sage.rings.complex_mpfr.ComplexField_class'>]

        sage: len(sage.rings.abc.ComplexField.__subclasses__()) <= 1
        True
    """

    pass


class ComplexBallField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.complex_arb.ComplexBallField`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(CBF, sage.rings.abc.ComplexBallField)
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.ComplexBallField.__subclasses__()
        [<class 'sage.rings.complex_arb.ComplexBallField'>]

        sage: len(sage.rings.abc.ComplexBallField.__subclasses__()) <= 1
        True
    """

    pass


class ComplexIntervalField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.complex_interval_field.ComplexIntervalField_class`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(CIF, sage.rings.abc.ComplexIntervalField)
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.ComplexIntervalField.__subclasses__()
        [<class 'sage.rings.complex_interval_field.ComplexIntervalField_class'>]

        sage: len(sage.rings.abc.ComplexIntervalField.__subclasses__()) <= 1
        True
    """

    pass


cdef class ComplexDoubleField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.complex_double.ComplexDoubleField_class`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(CDF, sage.rings.abc.ComplexDoubleField)
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.ComplexDoubleField.__subclasses__()
        [<class 'sage.rings.complex_double.ComplexDoubleField_class'>]

        sage: len(sage.rings.abc.ComplexDoubleField.__subclasses__()) <= 1
        True
    """

    pass


class IntegerModRing:
    r"""
    Abstract base class for :class:`~sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(Integers(7), sage.rings.abc.IntegerModRing)
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.IntegerModRing.__subclasses__()
        [<class 'sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic'>]

        sage: len(sage.rings.abc.IntegerModRing.__subclasses__()) <= 1
        True
    """

    pass


class Order:
    r"""
    Abstract base class for :class:`~sage.rings.number_field.order.Order`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: K.<a> = NumberField(x^2 + 1); O = K.order(2*a)                    # optional - sage.rings.number_field
        sage: isinstance(O, sage.rings.abc.Order)                               # optional - sage.rings.number_field
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.Order.__subclasses__()                             # optional - sage.rings.number_field
        [<class 'sage.rings.number_field.order.Order'>]

        sage: len(sage.rings.abc.Order.__subclasses__()) <= 1
        True
    """

    pass


class pAdicRing(EuclideanDomain):
    r"""
    Abstract base class for :class:`~sage.rings.padics.generic_nodes.pAdicRingGeneric`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(Zp(5), sage.rings.abc.pAdicRing)
        True
        sage: isinstance(Qp(5), sage.rings.abc.pAdicRing)
        False

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.pAdicRing.__subclasses__()
        [<class 'sage.rings.padics.generic_nodes.pAdicRingGeneric'>]

        sage: len(sage.rings.abc.pAdicRing.__subclasses__()) <= 1
        True
    """

    pass


class pAdicField(Field):
    r"""
    Abstract base class for :class:`~sage.rings.padics.generic_nodes.pAdicFieldGeneric`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(Zp(5), sage.rings.abc.pAdicField)
        False
        sage: isinstance(Qp(5), sage.rings.abc.pAdicField)
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.pAdicField.__subclasses__()
        [<class 'sage.rings.padics.generic_nodes.pAdicFieldGeneric'>]

        sage: len(sage.rings.abc.pAdicField.__subclasses__()) <= 1
        True
    """

    pass


cdef class SymbolicRing(CommutativeRing):
    r"""
    Abstract base class for :class:`~sage.rings.symbolic.ring.SymbolicRing`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: isinstance(SR, sage.rings.abc.SymbolicRing)                       # optional - sage.symbolic
        True

    By design, other than the abstract subclass :class:`~sage.rings.abc.CallableSymbolicExpressionRing`,
    there is only one direct implementation subclass::

        sage: sage.rings.abc.SymbolicRing.__subclasses__()                      # optional - sage.symbolic
        [<class 'sage.rings.abc.CallableSymbolicExpressionRing'>,
         <class 'sage.symbolic.ring.SymbolicRing'>]

        sage: len(sage.rings.abc.SymbolicRing.__subclasses__()) <= 2
        True
    """

    pass


class CallableSymbolicExpressionRing(SymbolicRing):
    r"""
    Abstract base class for :class:`~sage.rings.symbolic.callable.CallableSymbolicExpressionRing_class`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.rings.abc
        sage: f = x.function(x).parent()                                        # optional - sage.symbolic
        sage: isinstance(f, sage.rings.abc.CallableSymbolicExpressionRing)      # optional - sage.symbolic
        True

    By design, there is a unique direct subclass::

        sage: sage.rings.abc.CallableSymbolicExpressionRing.__subclasses__()    # optional - sage.symbolic
        [<class 'sage.symbolic.callable.CallableSymbolicExpressionRing_class'>]

        sage: len(sage.rings.abc.CallableSymbolicExpressionRing.__subclasses__()) <= 1
        True
    """

    pass
