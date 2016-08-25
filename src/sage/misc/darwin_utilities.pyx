# darwin_utilities extension module built only on OS_X >= 10.5 (see end of module_list.py)
from __future__ import print_function

cdef extern from "darwin_memory_usage.h":
    cdef unsigned long long darwin_virtual_size()

def darwin_memory_usage():
    r"""
    On Darwin >= 9 (i.e. OS X >= 10.5), returns the virtual size of the process in bytes.
    This will match what "top" reports as the VSIZE.
    Raises on all other platforms

    EXAMPLES
        sage: uname = os.uname()
        sage: if uname[0] == 'Darwin' and not uname[2].startswith('8.'):
        ...       from sage.misc.darwin_utilities import darwin_memory_usage
        ...       memory_usage = darwin_memory_usage()


        sage: uname = os.uname()
        sage: if uname[0] == 'Darwin' and uname[2].startswith('8.'):
        ...       try:
        ...           from sage.misc.darwin_utilities import darwin_memory_usage
        ...           memory_usage = darwin_memory_usage()
        ...           print("doctest failure!")
        ...           print("darwin_memory_usage() not implemented on platform Darwin 8 (i.e. OS X 10.4)")
        ...       except ImportError:
        ...           pass


        sage: uname = os.uname()
        sage: if uname[0] != 'Darwin':
        ...       try:
        ...           from sage.misc.darwin_utilities import darwin_memory_usage
        ...           memory_usage = darwin_memory_usage()
        ...           print("doctest failure!")
        ...           print("darwin_memory_usage() not implemented on platform %s" % uname[0])
        ...       except ImportError:
        ...           pass


    """
    return darwin_virtual_size()

