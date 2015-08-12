from sage.ext.interrupt.interrupt cimport *

cdef extern from 'sage/ext/interrupt/pxi.h':
    int import_sage__ext__interrupt__interrupt() except -1

# This *must* be done for every module using interrupt functions
# otherwise you will get segmentation faults.
import_sage__ext__interrupt__interrupt()
