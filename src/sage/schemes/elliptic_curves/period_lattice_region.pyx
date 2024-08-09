# sage.doctest: needs sage.rings.complex_double sage.symbolic
r"""
Regions in fundamental domains of period lattices

This module is used to represent sub-regions of a fundamental parallelogram
of the period lattice of an elliptic curve, used in computing minimum height
bounds.

In particular, these are the approximating sets ``S^{(v)}`` in section 3.2 of
Thotsaphon Thongjunthug's Ph.D. Thesis and paper [Tho2010]_.

AUTHORS:

- Robert Bradshaw (2010): initial version

- John Cremona (2014): added some docstrings and doctests
"""

#*****************************************************************************
#       Copyright (C) 2010 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import numpy as np
cimport numpy as np

from sage.rings.cif import CIF
from cpython.object cimport Py_EQ, Py_NE


cdef class PeriodicRegion:

    # The generators of the lattice, a complex numbers.
    cdef public w1
    cdef public w2

    # A two-dimensional array of (essentially) booleans specifying the subset
    # of parallelogram tiles that cover this region.
    cdef readonly data

    # Flag specifying whether bitmap represents the whole fundamental
    # parallelogram, or only half (in which case it is symmetric about
    # the center).
    cdef readonly bint full

    def __init__(self, w1, w2, data, full=True):
        """
        EXAMPLES::

            sage: import numpy as np
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: S = PeriodicRegion(CDF(2), CDF(2*I), np.zeros((4, 4)))
            sage: S.plot()                                                              # needs sage.plot
            Graphics object consisting of 1 graphics primitive
            sage: data = np.zeros((4, 4))
            sage: data[1,1] = True
            sage: S = PeriodicRegion(CDF(2), CDF(2*I+1), data)
            sage: S.plot()                                                              # needs sage.plot
            Graphics object consisting of 5 graphics primitives
        """
        if data.dtype is not np.int8:
            data = data.astype(np.int8)
        self.w1 = w1
        self.w2 = w2
        self.data = data
        self.full = full

    def is_empty(self):
        """
        Return whether this region is empty.

        EXAMPLES::

            sage: import numpy as np
            sage: if int(np.version.short_version[0]) > 1:
            ....:     np.set_printoptions(legacy="1.25")
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((4, 4))
            sage: PeriodicRegion(CDF(2), CDF(2*I), data).is_empty()
            True
            sage: data[1,1] = True
            sage: PeriodicRegion(CDF(2), CDF(2*I), data).is_empty()
            False
        """
        return self.data.sum() == 0

    def _ensure_full(self):
        """
        Ensure the bitmap in ``self.data`` represents the entire fundamental
        parallelogram, expanding symmetry by duplicating data if necessary.

        Mutates and returns ``self``.

        EXAMPLES::

            sage: import numpy as np
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((4, 4))
            sage: S = PeriodicRegion(CDF(2), CDF(2*I), data, full=False)
            sage: S.data.shape
            (4, 4)
            sage: S.full
            False
            sage: _ = S._ensure_full()
            sage: S.full
            True
            sage: S.data.shape
            (4, 8)
            sage: _ = S._ensure_full()
            sage: S.data.shape
            (4, 8)
        """
        if not self.full:
            rows, cols = self.data.shape
            new_data = np.ndarray((rows, 2*cols), self.data.dtype)
            new_data[:,0:cols] = self.data
            new_data[::-1,-1:cols-1:-1] = self.data
            self.data = new_data
            self.full = True
        return self

    def ds(self):
        """
        Return the sides of each parallelogram tile.

        EXAMPLES::

            sage: import numpy as np
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((4, 4))
            sage: S = PeriodicRegion(CDF(2), CDF(2*I), data, full=False)
            sage: S.ds()
            (0.5, 0.25*I)
            sage: _ = S._ensure_full()
            sage: S.ds()
            (0.5, 0.25*I)

            sage: data = np.zeros((8, 8))
            sage: S = PeriodicRegion(CDF(1), CDF(I + 1/2), data)
            sage: S.ds()
            (0.125, 0.0625 + 0.125*I)
        """
        m, n = self.data.shape
        return self.w1/m, self.w2/(n * (2-self.full))

    def verify(self, condition):
        r"""
        Given a condition that should hold for every line segment on
        the boundary, verify that it actually does so.

        INPUT:

        - ``condition`` -- boolean-valued function on `\CC`

        OUTPUT:

        boolean according to whether the condition holds for all lines on the
        boundary.

        EXAMPLES::

            sage: import numpy as np
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((4, 4))
            sage: data[1, 1] = True
            sage: S = PeriodicRegion(CDF(1), CDF(I), data)
            sage: S.border()
            [(1, 1, 0), (2, 1, 0), (1, 1, 1), (1, 2, 1)]
            sage: condition = lambda z: z.real().abs()<1/2
            sage: S.verify(condition)
            False
            sage: condition = lambda z: z.real().abs()<1
            sage: S.verify(condition)
            True
        """
        for line in self.border(False):
            if not condition(line):
                return False
        return True

    def refine(self, condition=None, int times=1):
        r"""
        Recursive function to refine the current tiling.

        INPUT:

        - ``condition`` -- function (default: ``None``); if not ``None``, only
          keep tiles in the refinement which satisfy the condition

        - ``times`` -- integer (default: 1); the number of times to refine.
          Each refinement step halves the mesh size.

        OUTPUT: the refined PeriodicRegion

        EXAMPLES::

            sage: import numpy as np
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((4, 4))
            sage: S = PeriodicRegion(CDF(2), CDF(2*I), data, full=False)
            sage: S.ds()
            (0.5, 0.25*I)
            sage: S = S.refine()
            sage: S.ds()
            (0.25, 0.125*I)
            sage: S = S.refine(2)
            sage: S.ds()
            (0.125, 0.0625*I)
        """
        if times <= 0:
            return self
        cdef int i, j, m, n
        cdef np.ndarray[np.npy_int8, ndim=2] less, fuzz, new_data
        m, n = self.data.shape
        if condition is None:
            new_data = np.ndarray((2*m, 2*n), self.data.dtype)
            new_data[0::2,0::2] = self.data
            new_data[1::2,0::2] = self.data
            new_data[0::2,1::2] = self.data
            new_data[1::2,1::2] = self.data
        else:
            more = self.expand().data
            less = self.contract().data
            fuzz = less ^ more
            new_data = np.zeros((2*m, 2*n), self.data.dtype)
            dw1, dw2 = self.ds()
            dw1 /= 2
            dw2 /= 2
            for i in range(2*m):
                for j in range(2*n):
                    if less[i//2, j//2]:
                        new_data[i,j] = True
                    elif fuzz[i//2, j//2]:
                        new_data[i,j] = condition(dw1*(i+.5) + dw2*(j+.5))
        return PeriodicRegion(self.w1, self.w2, new_data, self.full).refine(condition, times-1)

    def expand(self, bint corners=True):
        """
        Return a region containing this region by adding all neighbors of
        internal tiles.

        EXAMPLES::

            sage: import numpy as np
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((4, 4))
            sage: data[1,1] = True
            sage: S = PeriodicRegion(CDF(1), CDF(I + 1/2), data)
            sage: S.plot()                                                              # needs sage.plot
            Graphics object consisting of 5 graphics primitives
            sage: S.expand().plot()                                                     # needs sage.plot
            Graphics object consisting of 13 graphics primitives
            sage: S.expand().data
            array([[1, 1, 1, 0],
                   [1, 1, 1, 0],
                   [1, 1, 1, 0],
                   [0, 0, 0, 0]], dtype=int8)
            sage: S.expand(corners=False).plot()                                        # needs sage.plot
            Graphics object consisting of 13 graphics primitives
            sage: S.expand(corners=False).data
            array([[0, 1, 0, 0],
                   [1, 1, 1, 0],
                   [0, 1, 0, 0],
                   [0, 0, 0, 0]], dtype=int8)
        """
        cdef int i, j, m, n
        m, n = self.data.shape
        cdef np.ndarray[np.npy_int8, ndim=2] framed, new_data
        framed = frame_data(self.data, self.full)
        new_data = np.zeros((m+2, n+2), self.data.dtype)
        for i in range(m):
            for j in range(n):
                if framed[i, j]:
                    new_data[i, j] = True
                    new_data[i - 1, j] = True
                    new_data[i + 1, j] = True
                    new_data[i, j - 1] = True
                    new_data[i, j + 1] = True
                    if corners:
                        new_data[i - 1, j - 1] = True
                        new_data[i + 1, j - 1] = True
                        new_data[i + 1, j + 1] = True
                        new_data[i - 1, j + 1] = True
        return PeriodicRegion(self.w1, self.w2, unframe_data(new_data, self.full), self.full)

    def contract(self, corners=True):
        """
        Opposite (but not inverse) of expand; removes neighbors of complement.

        EXAMPLES::

            sage: import numpy as np
            sage: if int(np.version.short_version[0]) > 1:
            ....:     np.set_printoptions(legacy="1.25")
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((10, 10))
            sage: data[1:4,1:4] = True
            sage: S = PeriodicRegion(CDF(1), CDF(I + 1/2), data)
            sage: S.plot()                                                              # needs sage.plot
            Graphics object consisting of 13 graphics primitives
            sage: S.contract().plot()                                                   # needs sage.plot
            Graphics object consisting of 5 graphics primitives
            sage: S.contract().data.sum()
            1
            sage: S.contract().contract().is_empty()
            True
        """
        return ~(~self).expand(corners)

    def __contains__(self, z):
        """
        Return whether this region contains the given point.

        EXAMPLES::

            sage: import numpy as np
            sage: if int(np.version.short_version[0]) > 1:
            ....:     np.set_printoptions(legacy="1.25")
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((4, 4))
            sage: data[1,1] = True
            sage: S = PeriodicRegion(CDF(2), CDF(2*I), data, full=False)
            sage: CDF(0, 0) in S
            False
            sage: CDF(0.6, 0.6) in S
            True
            sage: CDF(1.6, 0.6) in S
            False
            sage: CDF(6.6, 4.6) in S
            True
            sage: CDF(6.6, -1.4) in S
            True

            sage: w1 = CDF(1.4)
            sage: w2 = CDF(1.2 * I - .3)
            sage: S = PeriodicRegion(w1, w2, data)
            sage: z = w1/2 + w2/2; z in S
            False
            sage: z = w1/3 + w2/3; z in S
            True
            sage: z = w1/3 + w2/3 + 5*w1 - 6*w2; z in S
            True
        """
        CC = self.w1.parent()
        RR = CC.construction()[1]
        basis_matrix = (RR**(2,2))((tuple(self.w1), tuple(self.w2))).transpose()
        i, j = basis_matrix.solve_right((RR**2)(tuple(CC(z))))
        # Get the fractional part.
        i -= int(i)
        j -= int(j)
        if i < 0:
            i += 1
        if j < 0:
            j += 1
        if not self.full and j > 0.5:
            j = 1 - j
            i = 1 - i
        m, n = self.data.shape
        return self.data[int(m * i), int(n * j)]

    def __truediv__(self, unsigned int n):
        """
        Return a new region of the same resolution that is the image
        of this region under the map z -> z/n.

        The resolution is the same, so some detail may be lost.  The result is
        always at worse a superset of the true image.

        EXAMPLES::

            sage: import numpy as np
            sage: if int(np.version.short_version[0]) > 1:
            ....:     np.set_printoptions(legacy="1.25")
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion

            sage: data = np.zeros((20, 20))
            sage: data[2, 2:12] = True
            sage: data[2:6, 2] = True
            sage: data[3, 3] = True
            sage: S = PeriodicRegion(CDF(1), CDF(I + 1/2), data)
            sage: S.plot()                                                              # needs sage.plot
            Graphics object consisting of 29 graphics primitives
            sage: (S / 2).plot()                                                        # needs sage.plot
            Graphics object consisting of 57 graphics primitives
            sage: (S / 3).plot()                                                        # needs sage.plot
            Graphics object consisting of 109 graphics primitives
            sage: (S / 2 / 3) == (S / 6) == (S / 3 / 2)
            True

            sage: data = np.zeros((100, 100))
            sage: data[2, 3] = True
            sage: data[7, 9] = True
            sage: S = PeriodicRegion(CDF(1), CDF(I), data)
            sage: inside = [.025 + .035j, .075 + .095j]
            sage: all(z in S for z in inside)
            True
            sage: all(z/2 + j/2 + k/2 in S/2 for z in inside for j in range(2) for k in range(2))
            True
            sage: all(z/3 + j/3 + k/3 in S/3 for z in inside for j in range(3) for k in range(3))
            True
            sage: outside = [.025 + .095j, .075 + .035j]
            sage: any(z in S for z in outside)
            False
            sage: any(z/2 + j/2 + k/2 in S/2 for z in outside for j in range(2) for k in range(2))
            False
            sage: any(z/3 + j/3 + k/3 in S/3 for z in outside for j in range(3) for k in range(3))
            False

            sage: (S / 1) is S
            True
            sage: S / 0
            Traceback (most recent call last):
            ...
            ValueError: divisor must be positive
            sage: S / (-1)
            Traceback (most recent call last):
            ...
            OverflowError: can...t convert negative value to unsigned int
        """
        cdef unsigned int i, j, a, b, rows, cols
        if n <= 1:
            if n == 1:
                return self
            raise ValueError("divisor must be positive")
        self._ensure_full()
        cdef np.ndarray[np.npy_int8, ndim=2] data, new_data
        data = self.data
        rows, cols = self.data.shape

        new_data = np.zeros(self.data.shape, self.data.dtype)
        for i in range(rows):
            for j in range(cols):
                if data[i,j]:
                    for a in range(n):
                        for b in range(n):
                            new_data[(a*rows+i)//n, (b*cols+j)//n] = data[i,j]
        return PeriodicRegion(self.w1, self.w2, new_data)

    def __invert__(self):
        """
        Return the complement of this region.

        EXAMPLES::

            sage: import numpy as np
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((4, 4))
            sage: data[1,1] = True
            sage: S = PeriodicRegion(CDF(1), CDF(I), data)
            sage: .3 + .3j in S
            True
            sage: .3 + .3j in ~S
            False
            sage: 0 in S
            False
            sage: 0 in ~S
            True
        """
        return PeriodicRegion(self.w1, self.w2, 1-self.data, self.full)

    def __and__(left, right):
        """
        Return the intersection of ``left`` and ``right``.

        EXAMPLES::

            sage: import numpy as np
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((4, 4))
            sage: data[1, 1:3] = True
            sage: S1 = PeriodicRegion(CDF(1), CDF(I), data)
            sage: data = np.zeros((4, 4))
            sage: data[1:3, 1] = True
            sage: S2 = PeriodicRegion(CDF(1), CDF(I), data)
            sage: (S1 & S2).data
            array([[0, 0, 0, 0],
                   [0, 1, 0, 0],
                   [0, 0, 0, 0],
                   [0, 0, 0, 0]], dtype=int8)
        """
        assert isinstance(left, PeriodicRegion) and isinstance(right, PeriodicRegion)
        assert left.w1 == right.w1 and left.w2 == right.w2
        if left.full ^ right.full:
            left._ensure_full()
            right._ensure_full()
        return PeriodicRegion(left.w1, left.w2, left.data & right.data, left.full)

    def __xor__(left, right):
        """
        Return the union of ``left`` and ``right`` minus the intersection of
        ``left`` and ``right``.

        EXAMPLES::

            sage: import numpy as np
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((4, 4))
            sage: data[1, 1:3] = True
            sage: S1 = PeriodicRegion(CDF(1), CDF(I), data)
            sage: data = np.zeros((4, 4))
            sage: data[1:3, 1] = True
            sage: S2 = PeriodicRegion(CDF(1), CDF(I), data)
            sage: (S1 ^^ S2).data
            array([[0, 0, 0, 0],
                   [0, 0, 1, 0],
                   [0, 1, 0, 0],
                   [0, 0, 0, 0]], dtype=int8)
        """
        assert isinstance(left, PeriodicRegion) and isinstance(right, PeriodicRegion)
        assert left.w1 == right.w1 and left.w2 == right.w2
        if left.full ^ right.full:
            left._ensure_full()
            right._ensure_full()
        return PeriodicRegion(left.w1, left.w2, left.data ^ right.data, left.full)

    def __richcmp__(left, right, op):
        """
        Compare two regions.

        .. NOTE:: This is good for equality but not an ordering relation.

        TESTS::

            sage: import numpy as np
            sage: if int(np.version.short_version[0]) > 1:
            ....:     np.set_printoptions(legacy="1.25")
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((4, 4))
            sage: data[1, 1] = True
            sage: S1 = PeriodicRegion(CDF(1), CDF(I), data)
            sage: data = np.zeros((4, 4))
            sage: data[1, 1] = True
            sage: S2 = PeriodicRegion(CDF(1), CDF(I), data)
            sage: data = np.zeros((4, 4))
            sage: data[2, 2] = True
            sage: S3 = PeriodicRegion(CDF(1), CDF(I), data)
            sage: S1 == S2
            True
            sage: S2 == S3
            False
        """
        if type(left) is not type(right) or op not in [Py_EQ, Py_NE]:
            return NotImplemented

        if (left.w1, left.w2) != (right.w1, right.w2):
            equal = False
        else:
            if left.full ^ right.full:
                left._ensure_full()
                right._ensure_full()
            equal = (left.data == right.data).all()

        if op is Py_EQ:
            return equal
        return not equal

    def border(self, raw=True):
        """
        Return the boundary of this region as set of tile boundaries.

        If raw is true, returns a list with respect to the internal bitmap,
        otherwise returns complex intervals covering the border.

        EXAMPLES::

            sage: import numpy as np
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((4, 4))
            sage: data[1, 1] = True
            sage: PeriodicRegion(CDF(1), CDF(I), data).border()
            [(1, 1, 0), (2, 1, 0), (1, 1, 1), (1, 2, 1)]
            sage: PeriodicRegion(CDF(2), CDF(I-1/2), data).border()
            [(1, 1, 0), (2, 1, 0), (1, 1, 1), (1, 2, 1)]

            sage: PeriodicRegion(CDF(1), CDF(I), data).border(raw=False)
            [0.25000000000000000? + 1.?*I,
             0.50000000000000000? + 1.?*I,
             1.? + 0.25000000000000000?*I,
             1.? + 0.50000000000000000?*I]
            sage: PeriodicRegion(CDF(2), CDF(I-1/2), data).border(raw=False)
            [0.3? + 1.?*I,
             0.8? + 1.?*I,
             1.? + 0.25000000000000000?*I,
             1.? + 0.50000000000000000?*I]

            sage: data[1:3, 2] = True
            sage: PeriodicRegion(CDF(1), CDF(I), data).border()
            [(1, 1, 0), (2, 1, 0), (1, 1, 1), (1, 2, 0), (1, 3, 1), (3, 2, 0), (2, 2, 1), (2, 3, 1)]
        """
        cdef np.ndarray[np.npy_int8, ndim=2] framed = frame_data(self.data, self.full)
        cdef int m, n
        m, n = self.data.shape
        cdef int i, j
        L = []
        for i in range(m):
            for j in range(n):
                if framed[i,j]:
                    if not framed[i-1, j]:
                        L.append((i, j, 0))
                    if not framed[i+1, j]:
                        L.append((i+1, j, 0))
                    if not framed[i, j-1]:
                        L.append((i, j, 1))
                    if not framed[i, j+1]:
                        L.append((i, j+1, 1))
        if not raw:
            dw1, dw2 = self.ds()
            for ix, (i, j, dir) in enumerate(L):
                L[ix] = CIF(i*dw1 + j*dw2).union(CIF((i+dir)*dw1 + (j+(1-dir))*dw2))
        return L

    def innermost_point(self):
        """
        Return a point well inside the region, specifically the center of
        (one of) the last tile(s) to be removed on contraction.

        EXAMPLES::

            sage: import numpy as np
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((10, 10))
            sage: data[1:4, 1:4] = True
            sage: data[1, 0:8] = True
            sage: S = PeriodicRegion(CDF(1), CDF(I+1/2), data)
            sage: S.innermost_point()
            0.375 + 0.25*I
            sage: S.plot() + point(S.innermost_point())                                 # needs sage.plot
            Graphics object consisting of 24 graphics primitives
        """
        if self.is_empty():
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError("region is empty")
        inside = self
        while not inside.is_empty():
            self, inside = inside, inside.contract(corners=False)
        i, j = tuple(np.transpose(self.data.nonzero())[0])
        dw1, dw2 = self.ds()
        return float(i + 0.5) * dw1 + float(j + 0.5) * dw2

    def plot(self, **kwds):
        """
        Plot this region in the fundamental lattice.  If ``full`` is ``False``, plots
        only the lower half.  Note that the true nature of this region is periodic.

        EXAMPLES::

            sage: import numpy as np
            sage: from sage.schemes.elliptic_curves.period_lattice_region import PeriodicRegion
            sage: data = np.zeros((10, 10))
            sage: data[2, 2:8] = True
            sage: data[2:5, 2] = True
            sage: data[3, 3] = True
            sage: S = PeriodicRegion(CDF(1), CDF(I + 1/2), data)
            sage: plot(S) + plot(S.expand(), rgbcolor=(1, 0, 1), thickness=2)           # needs sage.plot
            Graphics object consisting of 46 graphics primitives
        """
        from sage.plot.line import line
        dw1, dw2 = self.ds()
        L = []
        F = line([(0,0), tuple(self.w1),
                  tuple(self.w1+self.w2), tuple(self.w2), (0,0)])
        if not self.full:
            F += line([tuple(self.w2/2), tuple(self.w1+self.w2/2)])
        if 'rgbcolor' not in kwds:
            kwds['rgbcolor'] = 'red'
        for i, j, dir in self.border():
            ii, jj = i+dir, j+(1-dir)
            L.append(line([tuple(i*dw1 + j*dw2), tuple(ii*dw1 + jj*dw2)], **kwds))
        return sum(L, F)


cdef frame_data(data, bint full=True):
    """
    Helper function for PeriodicRegion.expand() and
    PeriodicRegion.border().  This makes "wrapping around" work
    transparently for symmetric regions (full=False) as well.

    Idea:

    [[a, b, c]        [[a, b, c, a, c],
     [d, e, f]   -->   [d, e, f, d, f],
     [g, h, i]]        [g, h, i, g, i],
                       [a, b, c, a, c],
                       [g, h, i, g, i]]
    """
    m, n = data.shape
    framed = np.empty((m+2, n+2), data.dtype)
    # center
    framed[:-2,:-2] = data
    # top and bottom
    if full:
        framed[:-2,-1] = data[:,-1]
        framed[:-2,-2] = data[:, 0]
    else:
        framed[:-2,-1] = data[::-1, 0]
        framed[:-2,-2] = data[::-1,-1]
    # left and right
    framed[-2,:] = framed[0,:]
    framed[-1,:] = framed[-3,:]
    return framed

cdef unframe_data(framed, bint full=True):
    """
    Helper function for PeriodicRegion.expand().  This glues the
    borders together using the "or" operator.
    """
    framed = framed.copy()
    framed[0,:] |= framed[-2,:]
    framed[-3,:] |= framed[-1,:]
    if full:
        framed[:-2,-3] |= framed[:-2,-1]
        framed[:-2, 0] |= framed[:-2,-2]
    else:
        framed[-3::-1, 0] |= framed[:-2,-1]
        framed[-3::-1,-3] |= framed[:-2,-2]
    return framed[:-2,:-2]
