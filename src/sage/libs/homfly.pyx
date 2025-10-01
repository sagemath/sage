# distutils: libraries = homfly gc
r"""
Cython wrapper for libhomfly library

This is used to call the libhomfly library directly from python. Knots
and links are passed following the convention in libhomfly. It is basically
the oriented Gauss code, represented as a string of integers separated
by spaces as follows:

- how many strings,

  - for each string, how many crossings, then

    - for each crossing, the cross name, then `1` if over, `-1` if under

- for each crossing, the name of the crossing and `1` if right, `-1` if left.

If there are `n` crossings, they must be named `0, 1, ..., n-1`.

AUTHORS:

- Miguel Marco (2015-03-24): initial version.
"""

# ****************************************************************************
#       Copyright (C) 2015 Miguel Marco  <mmarco@unizar.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at youroption) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from sage.cpython.string cimport str_to_bytes, char_to_str

cdef extern from "homfly.h":
    ctypedef int  word
    ctypedef signed long int sb4
    ctypedef unsigned short int ub2
    ctypedef signed short int sb2
    struct Term:
        sb4 coef
        sb2 m
        sb2 l
    struct Poly:
        Term* term
        sb4 len
    Poly* homfly(char *argv)
    char* homfly_str(char *argv)


def homfly_polynomial_string(link):
    r"""
    Return the HOMFLY polynomial of a link.

    INPUT:

    - ``link`` -- string of space-separated integers representing the link

    OUTPUT: string with the HOMFLY polynomial in the variables `M` and `L`

    EXAMPLES::

        sage: from sage.libs.homfly import homfly_polynomial_string
        sage: trefoil = '1 6 0 1  1 -1  2 1  0 -1  1 1  2 -1 0 1 1 1 2 1'
        sage: homfly_polynomial_string(trefoil)
        ' - L^-4 - 2L^-2 + M^2L^-2'
    """
    link = str_to_bytes(link)
    sig_on()
    cdef char* c_output = homfly_str(link)
    sig_off()
    return char_to_str(c_output)


def homfly_polynomial_dict(link):
    """
    Return a dictionary representing the HOMFLY polynomial of a link.

    INPUT:

    - ``link`` -- string of space-separated integers representing the link

    EXAMPLES::

        sage: from sage.libs.homfly import homfly_polynomial_dict
        sage: trefoil = '1 6 0 1  1 -1  2 1  0 -1  1 1  2 -1 0 1 1 1 2 1'
        sage: homfly_polynomial_dict(trefoil)
        {(-4, 0): -1, (-2, 0): -2, (-2, 2): 1}
    """
    link = str_to_bytes(link)
    cdef Term ter
    sig_on()
    cdef Poly* c_output = homfly(link)
    sig_off()
    cdef int l = c_output.len
    d = {}
    for i in range(l):
        ter = c_output.term[i]
        d[(int(ter.l), int(ter.m))] = int(ter.coef)
    return d
