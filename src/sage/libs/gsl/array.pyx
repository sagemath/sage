"""
GSL arrays
"""

from cysignals.memory cimport sig_malloc, sig_free

cdef class GSLDoubleArray:
    """
    EXAMPLES::

        sage: a = WaveletTransform(128,'daubechies',4)
        sage: for i in range(1, 11):
        ....:     a[i] = 1
        sage: a[:6:2]
        [0.0, 1.0, 1.0]
    """
    def __init__(self, size_t n, size_t stride=1, data=None):
        """
        EXAMPLES::

            sage: from sage.libs.gsl.array import GSLDoubleArray
            sage: a = GSLDoubleArray(10)
        """
        cdef int i

        self.n = n
        self.stride = stride
        self.data = <double *> sig_malloc(sizeof(double)*n)
        if data is not None:
            for i from 0 <= i < n:
                self.data[i] = data[i]
        else:
            for i from 0 <= i < n:
                self.data[i] = 0

    def __dealloc__(self):
        """
        EXAMPLES::

            sage: from sage.libs.gsl.array import GSLDoubleArray
            sage: a = GSLDoubleArray(10)
            sage: del a
        """
        sig_free(self.data)

    def __len__(self):
        """
        EXAMPLES::

            sage: from sage.libs.gsl.array import GSLDoubleArray
            sage: a = GSLDoubleArray(10)
            sage: len(a)
            10
        """
        return self.n

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.libs.gsl.array import GSLDoubleArray
            sage: a = GSLDoubleArray(10)
            sage: a
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        """
        return str(list(self))

    def __setitem__(self, size_t i, x):
        """
        EXAMPLES::

            sage: from sage.libs.gsl.array import GSLDoubleArray
            sage: a = GSLDoubleArray(10)
            sage: a[5] = 3
            sage: a
            [0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0]
        """
        # just set real for now
        if i < 0 or i >= self.n:
            raise IndexError
        self.data[i] = x

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: from sage.libs.gsl.array import GSLDoubleArray
            sage: a = GSLDoubleArray(10)
            sage: for i in range(10):
            ....:     a[i] = i
            sage: a[3:7]
            [3.0, 4.0, 5.0, 6.0]
        """
        if isinstance(i, slice):
            start, stop, step = i.indices(len(self))
            # TODO -- make this actually fast.
            return list(self)[start:stop:step]
        else:
            if i < 0 or i >= self.n:
                raise IndexError
            return self.data[i]
