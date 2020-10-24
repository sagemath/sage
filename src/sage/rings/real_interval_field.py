"""
Field of Arbitrary Precision Real Number Intervals
"""

from sage.misc.superseded import deprecation
deprecation(24371, "sage.rings.real_interval_field is deprecated")

from sage.rings.real_mpfi import RealIntervalField_class, RealIntervalFieldElement

def is_RealIntervalField(x):
    """
    Check if ``x`` is a :class:`RealIntervalField_class`.

    EXAMPLES::

        sage: from sage.rings.real_interval_field import is_RealIntervalField as is_RIF
        doctest:...: DeprecationWarning: sage.rings.real_interval_field is deprecated
        See http://trac.sagemath.org/24371 for details.
        sage: is_RIF(RIF)
        True
    """
    return isinstance(x, RealIntervalField_class)

def is_RealIntervalFieldElement(x):
    """
    Check if ``x`` is a :class:`RealIntervalFieldElement`.

    EXAMPLES::

        sage: from sage.rings.real_interval_field import is_RealIntervalFieldElement as is_RIFE
        sage: is_RIFE(RIF(2.5))
        True
    """
    return isinstance(x, RealIntervalFieldElement)

def __reduce__RealIntervalField(prec, sci_not):
    """
    For pickling.

    EXAMPLES::

        sage: from sage.rings.real_interval_field import __reduce__RealIntervalField
        sage: R = RealIntervalField(sci_not=1, prec=200)
        sage: __reduce__RealIntervalField(200, 1) == R
        True
    """
    return RealIntervalField_class(prec, sci_not)

def __reduce__RealIntervalFieldElement(parent, x, base=10):
    """
    For pickling.

    EXAMPLES::

        sage: from sage.rings.real_interval_field import __reduce__RealIntervalFieldElement
        sage: R = RealIntervalField(sci_not=1, prec=200)
        sage: elt = R(2.5)
        sage: __reduce__RealIntervalFieldElement(R, 2.5) == elt
        True
    """
    return RealIntervalFieldElement(parent, x, base=base)
