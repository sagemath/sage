# sage_setup: distribution = sagemath-modules
# sage.doctest: optional - numpy
r"""
Dense real double vectors using a NumPy backend

EXAMPLES::

    sage: # needs sage.symbolic
    sage: v = vector(RDF, [1, pi, sqrt(2)]); v
    (1.0, 3.141592653589793, 1.414213562373095)
    sage: type(v)
    <class 'sage.modules.vector_real_double_dense.Vector_real_double_dense'>
    sage: parent(v)
    Vector space of dimension 3 over Real Double Field
    sage: v[0] = 5
    sage: v
    (5.0, 3.141592653589793, 1.414213562373095)
    sage: loads(dumps(v)) == v
    True

AUTHORS:
    -- Jason Grout, Oct 2008: switch to numpy backend, factored out
       Vector_double_dense class
"""

#*****************************************************************************
#       Copyright (C) 2008 Jason Grout <jason-sage@creativetrax.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.real_double import RDF

cimport numpy


cdef class Vector_real_double_dense(Vector_double_dense):
    """
    Vectors over the Real Double Field.  These are supposed to be fast
    vector operations using C doubles. Most operations are implemented
    using numpy which will call the underlying BLAS, if needed, on the
    system.

    EXAMPLES::

        sage: v = vector(RDF, [1,2,3,4]); v
        (1.0, 2.0, 3.0, 4.0)
        sage: v*v
        30.0
    """
    def __cinit__(self, parent, entries, coerce=True, copy=True):
        self._numpy_dtype = numpy.dtype('float64')
        self._numpy_dtypeint = numpy.NPY_DOUBLE
        self._python_dtype = float
        # TODO: Make RealDoubleElement instead of RDF for speed
        self._sage_dtype = RDF
        self.__create_vector__()
        return

    def stats_skew(self):
        """
        Compute the skewness of a data set.

        For normally distributed data, the skewness should be about
        0. A skewness value > 0 means that there is more weight in the
        left tail of the distribution. (Paragraph from the scipy.stats
        docstring.)

        EXAMPLES::

            sage: v = vector(RDF, range(9))
            sage: v.stats_skew()                                                        # needs scipy
            0.0
        """
        import scipy.stats
        return self._sage_dtype(scipy.stats.skew(self._vector_numpy))

    def __reduce__(self):
        """
        Pickling.

        EXAMPLES::

            sage: a = vector(RDF, range(9))
            sage: loads(dumps(a)) == a
            True
        """
        return (unpickle_v1, (self._parent, self.list(), self._degree,
                              not self._is_immutable))


# For backwards compatibility, we must keep the function unpickle_v0
def unpickle_v0(parent, entries, degree):
    """
    Create a real double vector containing the entries.

    EXAMPLES::

        sage: v = vector(RDF, [1,2,3])
        sage: w = sage.modules.vector_real_double_dense.unpickle_v0(v.parent(), list(v), v.degree())
        sage: v == w
        True
    """
    return unpickle_v1(parent, entries, degree)


def unpickle_v1(parent, entries, degree, is_mutable=None):
    """
    Create a real double vector with the given parent, entries,
    degree, and mutability.

    EXAMPLES::

        sage: v = vector(RDF, [1,2,3])
        sage: w = sage.modules.vector_real_double_dense.unpickle_v1(v.parent(), list(v), v.degree(), v.is_immutable())
        sage: v == w
        True
    """
    cdef Vector_real_double_dense v = Vector_real_double_dense(parent, entries)
    if is_mutable is not None:
        v._is_immutable = not is_mutable
    return v
