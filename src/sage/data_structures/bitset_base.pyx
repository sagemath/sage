# sage_setup: distribution = sagemath-categories
"""
Few functions from ``bitset_base.pxd`` that are not inlined.
"""

#*****************************************************************************
#     Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef char* bitset_chars(char* s, fused_bitset_t bits, char zero=c'0', char one=c'1') noexcept:
    """
    Return a string representation of the bitset in s, using zero for
    the character representing the items not in the bitset and one for
    the character representing the items in the bitset.

    The string is both stored in s and returned.  If s is NULL, then a
    new string is allocated.
    """
    cdef mp_bitcnt_t i
    if s == NULL:
        s = <char *>sig_malloc(bits.size + 1)
    for i in range(bits.size):
        s[i] = one if bitset_in(bits, i) else zero
    s[bits.size] = 0
    return s

cdef int bitset_from_char(bitset_t bits, char* s, char zero=c'0', char one=c'1') except -1:
    """
    Initialize a bitset with a set derived from the C string s, where one
    represents the character indicating set membership.
    """
    bitset_init(bits, strlen(s))
    cdef mp_bitcnt_t i
    for i in range(bits.size):
        bitset_set_to(bits, i, s[i] == one)
    return 0

cdef int bitset_from_str(bitset_t bits, object s, char zero=c'0', char one=c'1') except -1:
    """
    Initialize a bitset with a set derived from the Python str s, where one
    represents the character indicating set membership.
    """
    cdef bytes b = str_to_bytes(s)
    return bitset_from_char(bits, b, zero, one)

cdef bitset_string(fused_bitset_t bits):
    """
    Return a python string representing the bitset.
    """
    return bytes_to_str(bitset_bytes(bits))

cdef bitset_bytes(fused_bitset_t bits):
    """
    Return a python bytes string representing the bitset.

    On Python 2 this is equivalent to bitset_string.
    """

    cdef char* s = bitset_chars(NULL, bits)
    cdef object py_s
    py_s = s
    sig_free(s)
    return py_s

cdef list bitset_list(fused_bitset_t bits):
    """
    Return a list of elements in the bitset.
    """
    cdef list elts = []
    cdef long elt = bitset_first(bits)
    while elt >= 0:
        elts.append(elt)
        elt = bitset_next(bits, elt + 1)
    return elts

cdef bitset_pickle(bitset_t bs):
    """
    Convert ``bs`` to a reasonably compact Python structure.

    Useful for pickling objects using bitsets as internal data structure.
    To ensure this works on 32-bit and 64-bit machines, the size of a long
    is stored too.
    """
    version = 0
    data = []
    for i in range(bs.limbs):
        data.append(bs.bits[i])
    return (version, bs.size, bs.limbs, sizeof(unsigned long), tuple(data))

cdef bitset_unpickle(bitset_t bs, tuple input):
    """
    Convert the data into a bitset.

    Companion of ``bitset_pickle()``. Assumption: ``bs`` has been initialized.
    """
    version, size, limbs, longsize, data = input
    if version != 0:
        raise TypeError("bitset was saved with newer version of Sage. Please upgrade.")
    if bs.size != size:
        bitset_realloc(bs, size)
    if sizeof(unsigned long) == longsize and bs.limbs == limbs:
        for i in range(bs.limbs):
            bs.bits[i] = data[i]
    else:
        storage = 8 * longsize  # number of elements encoded in one limb
        adder = 0
        bitset_clear(bs)
        for i in range(limbs):
            for j in range(storage):
                if (data[i] >> j) & 1:
                    bitset_add(bs, <mp_bitcnt_t>(j + adder))
            adder += storage
