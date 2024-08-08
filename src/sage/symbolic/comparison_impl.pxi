"""
Comparison of Symbolic Expressions

There are two useful ways to compare symbolic expressions:

* :func:`print_order` is how the terms are ordered. This is always
  defined. If you need a fast comparison, this is it.

* :func:`math_order` is the "mathematical" comparison. This may raise
  an exception if the answer is unknown (to Sage) or cannot, in
  principle, evaluated to a boolean (for example, if it involves
  symbolic variables). Can be very slow as it potentially calls
  Maxima to prove the inequality.

There is also a mixed version:

* :func:`mixed_order` which is print order if variables are present,
  and mathematical/numeric if not. This should enable quick and
  correct results.
"""

# ****************************************************************************
#       Copyright (C) 2014      Volker Braun
#       Copyright (C) 2016      Ralf Stephan
#       Copyright (C) 2017      Frédéric Chapoton
#       Copyright (C) 2021      Michael Orlitzky
#       Copyright (C) 2021      Jonathan Kliem
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython cimport *


cdef int print_order_c(Expression lhs, Expression rhs) noexcept:
    """
    Print comparison.

    See :meth:`print_order` for details.
    """
    return print_order_compare((<Expression>lhs)._gobj, (<Expression>rhs)._gobj)


cpdef int print_order(lhs, rhs) except -2:
    """
    Comparison in the print order

    INPUT:

    - ``lhs``, ``rhs`` -- two symbolic expressions or something that
      can be converted to one

    OUTPUT:

    Either `-1`, `0`, or `+1` indicating the comparison. An exception
    is raised if the arguments cannot be converted into the symbolic
    ring.

    EXAMPLES::

        sage: from sage.symbolic.expression import print_order
        sage: print_order(1, oo)
        1
        sage: print_order(e, oo)
        -1
        sage: print_order(pi, oo)
        1
        sage: print_order(1, sqrt(2))
        1

    Check that :issue:`12967` is fixed::

        sage: print_order(SR(oo), sqrt(2))
        1
    """
    if not isinstance(lhs, Expression):
        from sage.symbolic.ring import SR
        lhs = SR(lhs)
    if not isinstance(rhs, Expression):
        from sage.symbolic.ring import SR
        rhs = SR(rhs)
    return print_order_c(lhs, rhs)


class _print_key():

    def __init__(self, ex):
        """
        Sort key to sort in print order.

        INPUT:

        - ``ex`` -- symbolic expression or something that can be
          converted into one

        EXAMPLES::

            sage: from sage.symbolic.expression import _print_key
            sage: _print_key(1)
            <sage.symbolic.expression._print_key object at 0x...>
        """
        if not isinstance(ex, Expression):
            from sage.symbolic.ring import SR
            ex = SR(ex)
        self.ex = ex

    def __lt__(self, other):
        """
        Implement "less than" to make the key comparable.

        INPUT:

        - ``other`` -- another :class:`_print_key` instance

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.symbolic.expression import print_order, _print_key
            sage: print_order(1, 2)
            -1
            sage: _print_key(1) < _print_key(2)
            True
            sage: print_order(1, sqrt(2))
            1
            sage: _print_key(1) < _print_key(sqrt(2))
            False
        """
        return print_order_c(self.ex, other.ex) < 0


cpdef print_sorted(expressions):
    """
    Sort a list in print order.

    INPUT:

    - ``expressions`` -- list/tuple/iterable of symbolic
      expressions, or something that can be converted to one

    OUTPUT: the list sorted by :meth:`print_order`

    EXAMPLES::

        sage: from sage.symbolic.expression import print_sorted
        sage: print_sorted([SR(1), SR(e), SR(pi), sqrt(2)])
        [e, sqrt(2), pi, 1]
    """
    return sorted(expressions, key=_print_key)


class _math_key():

    def __init__(self, ex):
        """
        Sort key to sort in "Mathematics" order.

        INPUT:

        - ``ex`` -- symbolic expression or something that can be
          converted into one

        EXAMPLES::

            sage: from sage.symbolic.expression import _math_key
            sage: _math_key(1)
            <sage.symbolic.expression._math_key object at 0x...>
        """
        if not isinstance(ex, Expression):
            from sage.symbolic.ring import SR
            ex = SR(ex)
        self.ex = ex

    def __lt__(self, other):
        """
        Implement "less than" to make the key comparable.

        INPUT:

        - ``other`` -- another :class:`_print_key` instance

        OUTPUT: boolean; a :exc:`ValueError` is raised if we do not know how to
        perform the comparison

        EXAMPLES::

            sage: from sage.symbolic.expression import _math_key
            sage: _math_key(1) < _math_key(2)
            True
            sage: _math_key(1) < _math_key(sqrt(2))
            True

        Check that :issue:`12967` is fixed::

            sage: _math_key(1) < _math_key(oo)
            True
        """
        less_than = bool(self.ex < other.ex)
        greater_than = bool(self.ex > other.ex)
        if less_than:
            if not greater_than:
                return True
            else:
                assert False     # unreachable
        else:
            if greater_than:
                return False
            else:
                raise ValueError('cannot compare {0} and {1}'.format(self.ex, other.ex))


cpdef math_sorted(expressions):
    """
    Sort a list of symbolic numbers in the "Mathematics" order.

    INPUT:

    - ``expressions`` -- list/tuple/iterable of symbolic
      expressions, or something that can be converted to one

    OUTPUT:

    The list sorted by ascending (real) value. If an entry does not
    define a real value (or plus/minus infinity), or if the comparison
    is not known, a :exc:`ValueError` is raised.

    EXAMPLES::

        sage: from sage.symbolic.expression import math_sorted
        sage: math_sorted([SR(1), SR(e), SR(pi), sqrt(2)])
        [1, sqrt(2), e, pi]
    """
    return sorted(expressions, key=_math_key)


cpdef int mixed_order(lhs, rhs) except -2:
    """
    Comparison in the mixed order

    INPUT:

    - ``lhs``, ``rhs`` -- two symbolic expressions or something that
      can be converted to one

    OUTPUT:

    Either `-1`, `0`, or `+1` indicating the comparison. An exception
    is raised if the arguments cannot be converted into the symbolic
    ring.

    EXAMPLES::

        sage: from sage.symbolic.expression import mixed_order
        sage: mixed_order(1, oo)
        -1
        sage: mixed_order(e, oo)
        -1
        sage: mixed_order(pi, oo)
        -1
        sage: mixed_order(1, sqrt(2))
        -1
        sage: mixed_order(x + x^2, x*(x+1))
        -1

    Check that :issue:`12967` is fixed::

        sage: mixed_order(SR(oo), sqrt(2))
        1

    Ensure that :issue:`32185` is fixed::

        sage: mixed_order(pi, 0)
        1
        sage: mixed_order(golden_ratio, 0)
        1
        sage: mixed_order(log2, 0)
        1
    """
    if lhs is rhs:
        return 0
    if not isinstance(lhs, Expression):
        from sage.symbolic.ring import SR
        lhs = SR(lhs)
    if not isinstance(rhs, Expression):
        from sage.symbolic.ring import SR
        rhs = SR(rhs)
    less_than = _mixed_key(lhs) < _mixed_key(rhs)
    if less_than:
        return -1
    greater_than = _mixed_key(lhs) > _mixed_key(rhs)
    if greater_than:
        return 1
    else:
        return 0


class _mixed_key():

    def __init__(self, ex):
        """
        Sort key to sort in mixed order.

        Mixed order is print order if variables are present,
        mathematical/numeric if not. This should enable quick
        and correct results.

        INPUT:

        - ``ex`` -- symbolic expression or something that can be
          converted into one

        EXAMPLES::

            sage: from sage.symbolic.expression import _mixed_key
            sage: _mixed_key(1)
            <sage.symbolic.expression._mixed_key object at 0x...>
        """
        if not isinstance(ex, Expression):
            from sage.symbolic.ring import SR
            ex = SR(ex)
        self.ex = ex

    def __lt__(self, other):
        """
        Implement "less than" to make the key comparable.

        INPUT:

        - ``other`` -- another :class:`_mixed_key` instance

        OUTPUT: boolean; a :exc:`ValueError` is raised if we do not know how to
        perform the comparison

        EXAMPLES::

            sage: from sage.symbolic.expression import _mixed_key
            sage: _mixed_key(1) < _mixed_key(2)
            True
            sage: _mixed_key(1) < _mixed_key(sqrt(2))
            True

        Check that :issue:`12967` is fixed::

            sage: _mixed_key(1) < _mixed_key(oo)
            True
        """
        from sage.rings.real_mpfi import RIF
        selfv = len(self.ex.variables())
        otherv = len(other.ex.variables())
        if selfv:
            if otherv:
                return _print_key(self.ex) < _print_key(other.ex)
            else:
                return False
        else:
            if otherv:
                return True

        # no variables involved from here on
        rel = self.ex < other.ex
        if (self.ex.is_infinity() or other.ex.is_infinity()):
            pynac_result = decide_relational((<Expression>rel)._gobj)
            if pynac_result == relational_undecidable:
                raise ValueError('cannot compare {0} and {1}'.format(self.ex, other.ex))
            return pynac_result == relational_true

        det_ex = self.ex - other.ex
        if not has_symbol_or_function((<Expression>rel)._gobj) and not det_ex.is_constant():
            while hasattr(det_ex, 'pyobject') and isinstance(det_ex, Expression):
                try:
                    det_ex = det_ex.pyobject()
                except TypeError:
                    break
            if not isinstance(det_ex, Expression):
                return det_ex < 0
            from sage.rings.qqbar import QQbar
            try:
                from sage.rings.qqbar import QQbar
                num = QQbar(det_ex)
            except (TypeError, AttributeError,ValueError,NotImplementedError):
                try:
                    num = det_ex.expand().n(RIF.prec()+5)
                except (TypeError, AttributeError):
                    raise ValueError('cannot compare {0} and {1}'.format(self.ex, other.ex))
                else:
                    return num < 0
            else:
                return num < 0

        # here we have expressions containing functions
        try:
            num = det_ex.expand().n(RIF.prec()+5)
        except (TypeError, AttributeError):
            raise ValueError('cannot compare {0} and {1}'.format(self.ex, other.ex))
        else:
            return num < 0


cpdef mixed_sorted(expressions):
    """
    Sort a list of symbolic numbers in the "Mixed" order.

    INPUT:

    - ``expressions`` -- list/tuple/iterable of symbolic
      expressions, or something that can be converted to one

    OUTPUT:

    In the list the numeric values are sorted by ascending (real) value,
    and the expressions with variables according to print order.
    If an entry does not
    define a real value (or plus/minus infinity), or if the comparison
    is not known, a :exc:`ValueError` is raised.

    EXAMPLES::

        sage: from sage.symbolic.expression import mixed_sorted
        sage: mixed_sorted([SR(1), SR(e), SR(pi), sqrt(2), x, sqrt(x), sin(1/x)])
        [1, sqrt(2), e, pi, sin(1/x), sqrt(x), x]
    """
    return sorted(expressions, key=_mixed_key)
