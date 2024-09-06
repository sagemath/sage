# distutils: depends = NTL/ZZ.h

from cysignals.memory cimport sig_free
from cysignals.signals cimport sig_off

# Unset the signal handler and create a string from the buffer,
# then free the memory in the buffer.
cdef extern from *:
    void del_charstar "delete[]"(char*)

cdef object string(char* s):
    """
    Take a char* allocated using malloc, and converts it to a Python
    string, then deletes the allocated memory.  Also unsets the signal
    handler, so you *must* call sig_on() right before calling this!
    """
    sig_off()
    # Makes a python string and deletes what is pointed to by s.
    t = str(s)
    sig_free(s)
    return t

cdef object string_delete(char* s):
    """
    Take a char* allocated using C++ new, and converts it to a Python
    string, then deletes the allocated memory.  Also unsets the signal
    handler, so you *must* call sig_on() right before calling this!
    """
    sig_off()
    t = str(s)
    del_charstar(s)
    return t
