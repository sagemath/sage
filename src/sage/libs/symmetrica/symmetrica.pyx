# distutils: libraries = symmetrica
"""
Symmetrica library
"""

# According to https://gitlab.com/sagemath/symmetrica/-/blob/master/README.md#resolving-minmax-name-conflicts,
# we need to make sure that min and max macros are not defined on Windows.
# Usually, this can be done by setting NOMINMAX before including any Windows headers;
# however, we actually don't include any Windows headers here, but still get the error (even with NOMINMAX passed to the compiler).
# So just undefine min and max here.
cdef extern from *:
    """
    #if defined(_WIN32) || defined(WIN32) || defined(MS_WINDOWS)
    #if defined(min)
    #undef min
    #endif
    #if defined(max)
    #undef max
    #endif
    #endif
    """
    pass

include "symmetrica.pxi"

include "kostka.pxi"

include "sab.pxi"

include "sc.pxi"

include "sb.pxi"

include "part.pxi"

include "schur.pxi"

include "plet.pxi"
