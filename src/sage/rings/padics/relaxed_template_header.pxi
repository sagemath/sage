# sage_setup: distribution = sagemath-pari
"""
This file provides the declaration for the RelaxedElement class,
which collects common functionality for the different relaxed `p`-adic
template classes.

It is included in padic_relaxed_element.pxd and should be included
in any pxd file implementing relaxed `p`-adics.

AUTHORS:

- Xavier Caruso (2021-02) -- initial version
"""

#*****************************************************************************
#       Copyright (C) 2021 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer cimport Integer
from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement

cdef enum expansion_mode:
    simple_mode, smallest_mode, teichmuller_mode

cdef class RelaxedElement(pAdicGenericElement):
    cdef long _valuation
    cdef long _precrel
    cdef long _precbound
    cdef PowComputer_class prime_pow

    cdef cdigit_ptr _getdigit_relative(self, long i) noexcept
    cdef cdigit_ptr _getdigit_absolute(self, long i) noexcept
    cdef void _getslice_relative(self, celement slice, long start, long length) noexcept

    cdef int _init_jump(self) except -1
    cdef int _jump_c(self, long prec) noexcept
    cdef int _jump_relative_c(self, long prec, long halt) noexcept
    cdef int _next_c(self) noexcept

    cdef long valuation_c(self, long halt=*) noexcept
    cdef bint _is_equal(self, RelaxedElement right, long prec, bint permissive) except -1

cdef class RelaxedElement_abandon(RelaxedElement):
    pass
cdef relaxedelement_abandon

cdef class RelaxedElementWithDigits(RelaxedElement):
    cdef celement _digits


# Assignment

cdef class RelaxedElement_zero(RelaxedElement):
    pass

cdef class RelaxedElement_one(RelaxedElementWithDigits):
    pass

cdef class RelaxedElement_bound(RelaxedElement):
    cdef RelaxedElement _x

cdef class RelaxedElement_value(RelaxedElementWithDigits):
    cdef long _valuebound
    cdef long _shift
    cdef _value

cdef class RelaxedElement_random(RelaxedElementWithDigits):
    cdef randgen _generator
    # for pickling
    cdef long _initialvaluation
    cdef long _seed


# Operations

cdef class RelaxedElement_slice(RelaxedElement):
    cdef RelaxedElement _x
    cdef long _start
    cdef long _stop
    cdef long _shift

cdef class RelaxedElement_add(RelaxedElementWithDigits):
    cdef RelaxedElement _x
    cdef RelaxedElement _y

cdef class RelaxedElement_sub(RelaxedElementWithDigits):
    cdef RelaxedElement _x
    cdef RelaxedElement _y

cdef class RelaxedElement_mul(RelaxedElementWithDigits):
    cdef RelaxedElement _x
    cdef cdigit _lastdigit_x
    cdef RelaxedElement _y
    cdef cdigit _lastdigit_y
    cdef int _update_last_digit(self) noexcept

cdef class RelaxedElement_muldigit(RelaxedElementWithDigits):
    cdef cdigit_ptr _x
    cdef RelaxedElement _y

cdef class RelaxedElement_div(RelaxedElementWithDigits):
    cdef long _maxprec
    cdef cdigit _inverse
    cdef RelaxedElement _num
    cdef RelaxedElement _denom
    cdef RelaxedElement _definition
    cdef int _bootstrap_c(self) noexcept
    cdef bint _bootstraping

cdef class RelaxedElement_sqrt(RelaxedElementWithDigits):
    cdef RelaxedElement _x
    cdef RelaxedElement _definition
    cdef int _bootstrap_c(self) noexcept

cdef class RelaxedElement_teichmuller(RelaxedElementWithDigits):
    cdef bint _ready
    cdef bint _trivial
    cdef list _xns
    cdef RelaxedElement _xbar
    cdef RelaxedElement _xp

# Self-referent numbers

cdef class RelaxedElement_unknown(RelaxedElementWithDigits):
    cdef RelaxedElement _definition
    cdef long _next
    cpdef set(self, RelaxedElement definition)
    # for pickling
    cdef long _initialvaluation
    cdef long _initialprecrel

# Expansion

cdef class RelaxedElement_zeroone(RelaxedElementWithDigits):
    cdef void _setdigit_to_zero(self) noexcept
    cdef void _setdigit_to_one(self) noexcept

cdef class ExpansionIter():
    cdef RelaxedElement elt
    cdef expansion_mode mode
    cdef long start
    cdef long stop
    cdef long current
    cdef cdigit digit
    # simple mode
    cdef _next_simple(self)
    # smallest mode
    cdef cdigit carry
    cdef _next_smallest(self)
    # teichmuller mode
    cdef RelaxedElement tail
    cdef dict coefficients
    cdef _next_teichmuller(self)
