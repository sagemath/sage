# sage_setup: distribution = sagemath-objects
# ****************************************************************************
#       Copyright (C) 2017 Erik M. Bray <erik.bray@lri.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cdef extern from "string_impl.h":
    str _cstr_to_str(const char* c, encoding, errors)
    bytes _str_to_bytes(s, encoding, errors)


cdef inline str char_to_str(const char* c, encoding=None, errors=None):
    r"""
    Convert a C string to a Python ``str``.
    """
    # Implemented in C to avoid relying on PY_MAJOR_VERSION
    # compile-time variable. We keep the Cython wrapper to deal with
    # the default arguments.
    return _cstr_to_str(c, encoding, errors)


cpdef inline str bytes_to_str(b, encoding=None, errors=None):
    r"""
    Convert ``bytes`` to ``str``.

    This decodes the given ``bytes`` to a Python 3 unicode ``str`` using
    the specified encoding.  It is a no-op on ``str`` input.

    EXAMPLES::

        sage: from sage.cpython.string import bytes_to_str
        sage: s = bytes_to_str(b'\xcf\x80')
        sage: s == u'π'
        True
        sage: bytes_to_str([])
        Traceback (most recent call last):
        ...
        TypeError: expected bytes, list found
    """
    if isinstance(b, str):
        return b
    if type(b) is not bytes:
        raise TypeError(f"expected bytes, {type(b).__name__} found")

    return _cstr_to_str(<bytes>b, encoding, errors)


cpdef inline bytes str_to_bytes(s, encoding=None, errors=None):
    r"""
    Convert ``str`` or ``unicode`` to ``bytes``.

    It encodes the given ``str`` to a Python 3 ``bytes``
    using the specified encoding.  It is a no-op on ``bytes`` input.

    EXAMPLES::

        sage: from sage.cpython.string import str_to_bytes
        sage: bs = [str_to_bytes(u'π')]
        sage: all(b == b'\xcf\x80' for b in bs)
        True
        sage: str_to_bytes([])
        Traceback (most recent call last):
        ...
        TypeError: expected str... list found
    """
    # Implemented in C to avoid relying on PY_MAJOR_VERSION
    # compile-time variable. We keep the Cython wrapper to deal with
    # the default arguments.
    if isinstance(s, bytes):
        return s
    return _str_to_bytes(s, encoding, errors)
