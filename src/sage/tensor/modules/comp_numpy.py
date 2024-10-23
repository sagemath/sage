r"""
Tensor Component Backend using numpy.ndarray

AUTHORS:

- Aman Moon (2024-07-02): initial version

"""

# ****************************************************************************
#       Copyright (C) 2024 Aman Moon
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.tensor.modules.comp import Components
import numpy as np


class ComponentNumpy(SageObject):

    def __init__(self, ring, frame, nb_indices, shape=None, start_index=0,
                 output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: ComponentNumpy(ZZ, [1,2,3], 2)
            2-indices numpy components w.r.t. (1, 2, 3)

        """
        self._ring = ring
        self._frame = tuple(frame) if all(isinstance(i, (list, tuple)) for i in frame) else (tuple(frame),)
        self._nid = nb_indices
        self._shape = (len(self._frame[0]),) * nb_indices if shape is None else tuple(shape)
        self._sindex = tuple(start_index) if isinstance(start_index, (list, tuple)) else (start_index,) * nb_indices
        self._output_formatter = output_formatter
        self._comp = np.zeros(shape=self._shape, dtype=np.float64)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: c = ComponentNumpy(ZZ, [1,2,3], 2)
            sage: c._repr_()
            '2-indices numpy components w.r.t. (1, 2, 3)'

        """
        description = str()
        if not all(self._shape[0] == dim for dim in self._shape):
            description += "{}-shaped ".format(self._shape)
        description += str(self._nid)
        if self._nid == 1:
            description += "-index "
        else:
            description += "-indices "
        description += "numpy components w.r.t. "
        if len(self._frame) == 1:
            description += str(self._frame[0])
        else:
            description += str(self._frame)

        return description

    def _new_instance(self):
        r"""
        Creates a :class:`ComponentNumpy` instance of the same number of indices
        and w.r.t. the same frame.

        This method must be redefined by derived classes of
        :class:`ComponentsNumpy`.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: c = ComponentNumpy(ZZ, [1,2,3], 2)
            sage: c._new_instance()
            2-indices numpy components w.r.t. (1, 2, 3)

        """
        return self.__class__(self._ring, self._frame, self._nid,
                    self._shape, self._sindex, self._output_formatter)

    def copy(self):
        r"""
        Return an exact copy of ``self``.

        EXAMPLES:

        Copy of a set of components with a single index::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: V = VectorSpace(QQ,3)
            sage: a = ComponentNumpy(QQ, V.basis(), 1)
            sage: a[:] = -2, 1, 5
            sage: b = a.copy() ; b
            1-index numpy components w.r.t. ((1, 0, 0), (0, 1, 0), (0, 0, 1))
            sage: b[:]
            array([-2.,  1.,  5.])
            sage: b == a
            True
            sage: b is a  # b is a distinct object
            False

        """
        result = self._new_instance()
        result._comp = np.copy(self._comp)
        return result

    @property
    def T(self):
        """
        Transpose of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: c = ComponentNumpy(ZZ, [1,2,3], 2)
            sage: c.__setitem__(slice(None), [[0, 1, 2], [3, 4, 5], [6, 7, 8]])
            sage: t = c.T; t
            2-indices numpy components w.r.t. (1, 2, 3)
            sage: t[:]
            array([[0., 3., 6.],
                   [1., 4., 7.],
                   [2., 5., 8.]])
            sage: c == t.T
            True

        """
        result = self.__class__(self._ring, self._frame, self._nid, self._shape[::-1],
                          self._sindex[::-1], self._output_formatter)
        result._comp = np.copy(self._comp.T)
        return result

    def _check_indices(self, indices):
        r"""
        Check the validity of a list of indices and returns a tuple from it

        INPUT:

        - ``indices`` -- list of indices (possibly a single integer if
          self is a 1-index object)

        OUTPUT:

        - a tuple containing valid indices

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: c = ComponentNumpy(ZZ, [1,2,3], 2)
            sage: c._check_indices((0,1))
            (0, 1)
            sage: c._check_indices([0,1])
            (0, 1)
            sage: c._check_indices([2,1])
            (2, 1)
            sage: c._check_indices([2,3])
            Traceback (most recent call last):
            ...
            IndexError: index out of range: 3 not in [0, 2]
            sage: c._check_indices(1)
            Traceback (most recent call last):
            ...
            ValueError: wrong number of indices: 2 expected, while 1 are provided
            sage: c._check_indices([1,2,3])
            Traceback (most recent call last):
            ...
            ValueError: wrong number of indices: 2 expected, while 3 are provided
            sage: c = ComponentNumpy(ZZ, [1,2,3], 2, start_index=[2,3])
            sage: c._check_indices((2,3))
            (0, 0)
            sage: c._check_indices((2,5))
            (0, 2)
            sage: c._check_indices((1,3))
            Traceback (most recent call last):
            ...
            IndexError: index out of range: 1 not in [2, 4]

        """
        if isinstance(indices, (int, Integer)):
            ind = (indices,)
        else:
            ind = tuple(indices)
        shape = self._shape
        if len(ind) != len(shape):
            raise ValueError(("wrong number of indices: {} expected,"
                             " while {} are provided").format(self._nid, len(ind)))
        start_index = self._sindex
        np_ind = tuple()
        for i in range(len(shape)):
            si = start_index[i]
            np_ind += (indices[i] - si,)
            imax = shape[i] + si - 1
            if ind[i] < si or ind[i] > imax:
                raise IndexError("index out of range: " +
                                 "{} not in [{}, {}]".format(ind[i], si, imax))
        return np_ind

    def __getitem__(self, args):
        r"""
        Returns the component corresponding to the given indices.

        INPUT:

        - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object) or the character ``:`` for the full list
          of components

        OUTPUT:

        - the component corresponding to ``args`` or, if ``args`` = ``:``,
          the full list of components, in the form ``T[i][j]...`` for the
          components `T_{ij...}` (for a 2-indices object, a matrix is returned)

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: c = ComponentNumpy(ZZ, [1,2,3], 2, start_index=[1,2])
            sage: c[1,2]    # unset components are zero
            0.0
            sage: c.__getitem__((1,2))
            0.0
            sage: c.__getitem__([1,2])
            0.0
            sage: c[1,2] = -4
            sage: c[1,2]
            -4.0
            sage: c.__getitem__((1,2))
            -4.0
            sage: c[:]
            array([[-4.,  0.,  0.],
                   [ 0.,  0.,  0.],
                   [ 0.,  0.,  0.]])
            sage: c.__getitem__(slice(None))
            array([[-4.,  0.,  0.],
                   [ 0.,  0.,  0.],
                   [ 0.,  0.,  0.]])
            sage: c = ComponentNumpy(ZZ, [1,2,3], 2, start_index=[1,2])
            sage: c[1,2] = 4; c[3,2] = 2
            sage: c[:]
            array([[4., 0., 0.],
                  [0., 0., 0.],
                  [2., 0., 0.]])
            sage: c[1:4:2]
            array([[4., 0., 0.],
                  [2., 0., 0.]])
            sage: c[1:5]
            Traceback (most recent call last):
            ...
            IndexError: [start:stop] not in range [1,4]

        """
        no_format = self._output_formatter is None
        format_type = None # default value, possibly redefined below
        if isinstance(args, list):  # case of [[...]] syntax
            no_format = True
            if isinstance(args[0], slice):
                indices = args[0]
            elif isinstance(args[0], (tuple, list)): # to ensure equivalence between
                indices = args[0]           # [[(i,j,...)]] or [[[i,j,...]]] and [[i,j,...]]
            else:
                indices = tuple(args)
        else:
            # Determining from the input the list of indices and the format
            if isinstance(args, (int, Integer, slice)):
                indices = args
            elif isinstance(args[0], slice):
                indices = args[0]
            elif len(args) == self._nid:
                indices = args
            else:
                format_type = args[-1]
                indices = args[:-1]
        if isinstance(indices, slice):
            start, stop = indices.start, indices.stop
            range_start, range_end = self._sindex[0], self._sindex[0] + self._shape[0]

            if indices.start is not None:
                start = indices.start - range_start
            else:
                start = 0
            if indices.stop is not None:
                stop = indices.stop - range_start
            else:
                stop = range_end - range_start
            if not ((0 <= start <= range_end - range_start - 1) and
                    (0 <= stop <= range_end - range_start)):
                raise IndexError("[start:stop] not in range [{},{}]"
                                 .format(range_start, range_end))

            return self._comp[slice(start, stop, indices.step)]
        else:
            ind = self._check_indices(indices)
            elem = self._comp[ind]
            from sage.rings.real_mpfr import RR
            if no_format:
                return elem
            elif format_type is None:
                return self._output_formatter(RR(elem))
            else:
                return self._output_formatter(RR(elem), format_type)

    def __setitem__(self, args, value):
        r"""
        Sets the component corresponding to the given indices.

        INPUT:

        - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object); if ``[:]`` is provided, all the
          components are set
        - ``value`` -- the value to be set or a list of values if
          ``args = [:]`` (``slice(None)``)

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: c = ComponentNumpy(ZZ, [1,2,3], 2)
            sage: c.__setitem__((0,1), -4)
            sage: c[:]
            array([[ 0., -4.,  0.],
                   [ 0.,  0.,  0.],
                   [ 0.,  0.,  0.]])
            sage: c[0,1] = -4
            sage: c[:]
            array([[ 0., -4.,  0.],
                   [ 0.,  0.,  0.],
                   [ 0.,  0.,  0.]])
            sage: c.__setitem__(slice(None), [[0, 1, 2], [3, 4, 5], [6, 7, 8]])
            sage: c[:]
            array([[0., 1., 2.],
                  [3., 4., 5.],
                  [6., 7., 8.]])
            sage: c = ComponentNumpy(ZZ, [1,2,3], 2, start_index=[1,2])
            sage: c.__setitem__(slice(None), [[0, 1, 2], [3, 4, 5], [6, 7, 8]])
            sage: c[1:4:2] = [[2, 1, 0], [5, 4, 3]]
            sage: c[:]
            array([[2., 1., 0.],
                   [3., 4., 5.],
                   [5., 4., 3.]])
            sage: c[2:5] = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
            Traceback (most recent call last):
            ...
            IndexError: [start:stop] not in range [1,4]

        """
        if isinstance(args, list):  # case of [[...]] syntax
            if isinstance(args[0], slice):
                indices = args[0]
            elif isinstance(args[0], (tuple, list)): # to ensure equivalence between
                indices = args[0]           # [[(i,j,...)]] or [[[i,j,...]]] and [[i,j,...]]
            else:
                indices = tuple(args)
        else:
            # Determining from the input the list of indices and the format
            if isinstance(args, (int, Integer, slice)):
                indices = args
            elif isinstance(args[0], slice):
                indices = args[0]
            elif len(args) == self._nid:
                indices = args

        if isinstance(indices, slice):
            start, stop = indices.start, indices.stop
            range_start, range_end = self._sindex[0], self._sindex[0] + self._shape[0]

            if indices.start is not None:
                start = indices.start - range_start
            else:
                start = 0
            if indices.stop is not None:
                stop = indices.stop - range_start
            else:
                stop = range_end - range_start
            if not ((0 <= start <= range_end - range_start - 1) and
                    (0 <= stop <= range_end - range_start)):
                raise IndexError("[start:stop] not in range [{},{}]"
                    .format(range_start, range_end))

            self._comp[start:stop:indices.step] = value
        else:
            ind = self._check_indices(indices)
            self._comp[ind] = value

    def _broadcast(self, other):
        r"""
        Broadcast self with other. refer ``https://data-apis.org/array-api/2021.12/API_specification/broadcasting.html``

        INPUT:

        - ``other`` -- components, on the same frame as ``self``

        OUTPUT:

        - the broadcasted component of ``self`` by ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: c_1 = ComponentNumpy(ZZ, [1,2,3], 2)
            sage: c_2 = ComponentNumpy(ZZ, [1,2,3], 3, (1,3,3))
            sage: c_1._broadcast(c_2)
            (1, 3, 3)-shaped 3-indices numpy components w.r.t. (1, 2, 3)
            sage: c_2._broadcast(c_1)
            (1, 3, 3)-shaped 3-indices numpy components w.r.t. (1, 2, 3)
            sage: c_2 = ComponentNumpy(ZZ, [1,2,3], 3, (1,2,3))
            sage: c_1._broadcast(c_2)
            Traceback (most recent call last):
            ...
            ValueError: Shapes (3, 3), (1, 2, 3) cannot be broadcasted
        """
        s_shape = self._shape
        o_shape = other._shape
        N1 = len(s_shape)
        N2 = len(o_shape)
        N = max(N1, N2)
        shape = [0] * N

        i = N - 1
        while i >= 0:
            n1 = N1 - N + i
            d1 = s_shape[n1] if n1 >= 0 else 1
            n2 = N2 - N + i
            d2 = o_shape[n2] if n2 >= 0 else 1
            if d1 == 1:
                shape[i] = d2
            elif d2 == 1:
                shape[i] = d1
            elif d1 == d2:
                shape[i] = d1
            else:
                raise ValueError("Shapes {}, {} cannot be broadcasted"
                                 .format(self._shape, other._shape))
            i -= 1

        frame = tuple(set(self._frame + other._frame))
        output_formatter = None
        if self._output_formatter == other._output_formatter: # output formatter is None if self and other had diff output formatter
            output_formatter = self._output_formatter

        return ComponentNumpy(self._ring, frame, len(shape), shape, output_formatter=output_formatter)

## Arithmetic Operators

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT:

        - an exact copy of ``self``

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: c = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: c[:] = 5, 0, -4
            sage: a = c.__pos__() ; a
            1-index numpy components w.r.t. (1, 2, 3)
            sage: a[:]
            array([ 5.,  0., -4.])
            sage: a == +c
            True
            sage: a == c
            True

        """
        return self.copy()

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - the opposite of the components represented by ``self``

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: c = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: c[:] = 5, 0, -4
            sage: a = c.__neg__() ; a
            1-index numpy components w.r.t. (1, 2, 3)
            sage: a[:]
            array([-5., -0.,  4.])
            sage: a == -c
            True

        """
        neg_com = self._new_instance()
        neg_com._comp = (- np.copy(self._comp))
        return neg_com

    def __add__(self, other):
        r"""
        Component addition.

        INPUT:

        - ``other`` -- components of the same number of indices and defined
          on the same frame as ``self``

        OUTPUT:

        - components resulting from the addition of ``self`` and ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__add__(b) ; s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s[:]
            array([5., 5., 3.])
            sage: s == a+b
            True
        """
        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        if not isinstance(other, (Components, ComponentNumpy)):
            raise TypeError("the second argument for the addition must be a " +
                            "an instance of Components")
        # if isinstance(other, CompWithSymNumpy):
        #     return other + self     # to deal properly with symmetries
        if self._shape != other._shape: # use tensor broadcasting
            ret = self._broadcast(other)
        else:
            ret = self.copy()
        ret._comp = self._comp + other._comp
        return ret

    def __radd__(self, other):
        r"""
        Reflected addition (addition on the right to `other``)

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__radd__(b) ; s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s[:]
            array([5., 5., 3.])
            sage: s == a+b
            True
            sage: s = 0 + a ; s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s == a
            True

        """
        return self + other

    def __sub__(self, other):
        r"""
        Component subtraction.

        INPUT:

        - ``other`` -- components, of the same type as ``self``

        OUTPUT:

        - components resulting from the subtraction of ``other`` from ``self``

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__sub__(b) ; s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s[:]
            array([-3., -5., -9.])
            sage: s == a - b
            True

        """
        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        return self + (-other)

    def __rsub__(self, other):
        r"""
        Reflected subtraction (subtraction from ``other``).

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__rsub__(b) ; s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s[:]
            array([3., 5., 9.])
            sage: s == b - a
            True
            sage: s = 0 - a ; s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s[:]
            array([-1., -0.,  3.])
            sage: s == -a
            True

        """
        return (-self) + other

    def __mul__(self, other):
        r"""
        Component tensor product.

        INPUT:

        - ``other`` -- components, on the same frame as ``self``

        OUTPUT:

        - the tensor product of ``self`` by ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: s = a.__mul__(3); s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s[:]
            array([  3.,   0., -9.])
            sage: b = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__mul__(b) ; s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s[:]
            array([  4.,   0., -18.])
            sage: s == a*b
            True
            sage: t = b*a
            sage: t == s
            True

        """
        if isinstance(other, (int, Integer)):
            ret = self.copy()
            ret._comp = self._comp * other
        else:
            if self._shape == other._shape:
                ret = self.copy()
            else:
                ret = self._broadcast(other)
            ret._comp = self._comp * other._comp
        return ret

    def __rmul__(self, other):
        r"""
        Reflected multiplication (multiplication on the left by ``other``).

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: s = a.__rmul__(2) ; s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s[:]
            array([ 2.,  0., -6.])
            sage: s == 2*a
            True
            sage: a.__rmul__(0) == 0
            True

        """
        return self * other

    def __truediv__(self, other):
        r"""
        Division (by a scalar).

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(QQ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: s = a.__truediv__(3) ; s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s[:]
            array([ 0.33333333,  0.        , -1.        ])
            sage: s == a/3
            True
            sage: 3*s == a
            True

        """
        if isinstance(other, (int, Integer)):
            ret = self.copy()
            ret._comp = self._comp / other
        else:
            if self._shape == other._shape:
                ret = self.copy()
            else:
                ret = self._broadcast(other)
            ret._comp = self._comp / other._comp
        return ret

    def __floordiv__(self, other):
        r"""
        Evaluates self_i // other_i for each element of Component
        ``self`` with the respective element of ``other``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(QQ, [1,2,3], 1)
            sage: a[:] = 2, 0, -3
            sage: s = a.__floordiv__(2); s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s[:]
            array([ 1.,  0., -2.])
            sage: b = ComponentNumpy(QQ, [1,2,3], 1)
            sage: b[:] = 1, 4, -2
            sage: s = a.__floordiv__(b); s[:]
            array([2., 0., 1.])

        """
        if isinstance(other, (int, Integer)):
            ret = self.copy()
            ret._comp = self._comp // other
        else:
            if self._shape == other._shape:
                ret = self.copy()
            else:
                ret = self._broadcast(other)
            ret._comp = self._comp // other._comp
        return ret

    def __mod__(self, other):
        r"""
        Evaluates self_i % other_i for each element of Component
        ``self`` with the respective element of ``other``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(QQ, [1,2,3], 1)
            sage: a[:] = 1, 2, -3
            sage: s = a.__mod__(2); s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s[:]
            array([1., 0., 1.])
            sage: b = ComponentNumpy(QQ, [1,2,3], 1)
            sage: b[:] = 1, 4, -2
            sage: s = a.__mod__(b); s[:]
            array([ 0.,  2., -1.])

        """
        if isinstance(other, (int, Integer)):
            ret = self.copy()
            ret._comp = self._comp % other
        else:
            if self._shape == other._shape:
                ret = self.copy()
            else:
                ret = self._broadcast(other)
            ret._comp = self._comp % other._comp
        return ret

    def __pow__(self, other):
        r"""
        Evaluates self_i ^ other_i for each element of Component
        ``self`` with the respective element of ``other``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(QQ, [1,2,3], 1)
            sage: a[:] = 3, 0, -3
            sage: s = a.__pow__(2); s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s[:]
            array([9., 0., 9.])
            sage: b = ComponentNumpy(QQ, [1,2,3], 1)
            sage: b[:] = -1, 4, -2
            sage: s = a.__pow__(b); s[:]
            array([0.33333333, 0.        , 0.11111111])

        """
        if isinstance(other, (int, Integer)):
            ret = self.copy()
            ret._comp = self._comp ** other
        else:
            if self._shape == other._shape:
                ret = self.copy()
            else:
                ret = self._broadcast(other)
            ret._comp = self._comp ** other._comp
        return ret

    def __abs__(self):
        r"""
        Evaluates |self_i| for each element of Component ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(QQ, [1,2,3], 1)
            sage: a[:] = 3, 0, -3
            sage: s = a.__abs__(); s
            1-index numpy components w.r.t. (1, 2, 3)
            sage: s[:]
            array([3., 0., 3.])

        """
        if np.all(self._comp >= 0):
            return self.copy()
        else:
            ret = self.copy()
            ret._comp = np.absolute(ret._comp)
            return ret

# Comparison Operators

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a set of components or 0

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other``,  or ``False`` otherwise

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: c = ComponentNumpy(ZZ, [1,2,3], 2)
            sage: c.__eq__(0)  # uninitialized components are zero
            True
            sage: c[0,1], c[1,2] = 5, -4
            sage: c.__eq__(0)
            False
            sage: c1 = ComponentNumpy(ZZ, [1,2,3], 2)
            sage: c1[0,1] = 5
            sage: c.__eq__(c1)
            False
            sage: c1[1,2] = -4
            sage: c.__eq__(c1)
            True
            sage: v = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: c.__eq__(v)
            False

        """
        if isinstance(other, (int, Integer)): # other is 0
            if other == 0:
                return np.all(self._comp == 0)
            else:
                raise TypeError("cannot compare a set of components to a number")
        else: # other is another Components
            if set(other._frame) != set(self._frame):
                return False
            if other._nid != self._nid:
                return False
            if other._sindex != self._sindex:
                return False
            if other._output_formatter != self._output_formatter:
                return False
            return np.all(other._comp == self._comp)

## In-place Operators

    def __iadd__(self, other):
        r"""
        Inplace addition operator.

        OUTPUT:

        - components resulting from the addition of ``self`` and ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: a+=b; a
            1-index numpy components w.r.t. (1, 2, 3)
            sage: a[:]
            array([5., 5., 3.])
            sage: b[:]
            array([4., 5., 6.])

        """
        return self + other

    def __isub__(self, other):
        r"""
        Inplace subtraction operator.

        OUTPUT:

        - components resulting from the subtraction of ``self`` and ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: a -= b; a
            1-index numpy components w.r.t. (1, 2, 3)
            sage: a[:]
            array([-3., -5., -9.])
            sage: b[:]
            array([4., 5., 6.])

        """
        return self - other

    def __imul__(self, other):
        r"""
        Inplace component tensor product.

        INPUT:

        - ``other`` -- components, on the same frame as ``self``

        OUTPUT:

        - the tensor product of ``self`` by ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: s = a.__mul__(3)
            sage: a *= 3; a
            1-index numpy components w.r.t. (1, 2, 3)
            sage: a[:]
            array([  3.,   0., -9.])
            sage: a == s
            True
            sage: b = ComponentNumpy(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b[:] = 4, 5, 6
            sage: s = a.__mul__(b)
            sage: a *= b; a
            1-index numpy components w.r.t. (1, 2, 3)
            sage: a[:]
            array([  4.,   0., -18.])
            sage: a == s
            True

        """
        return self * other

    def __itruediv__(self, other):
        r"""
        In-place division (by a scalar).

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(QQ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: s = a.__truediv__(3)
            sage: a /= 3; a
            1-index numpy components w.r.t. (1, 2, 3)
            sage: a[:]
            array([ 0.33333333,  0.        , -1.        ])
            sage: a == s
            True

        """
        return self / other

    def __ifloordiv__(self, other):
        r"""
        In-place operator for evaluating self_i // other_i for each element of Component
        ``self`` with the respective element of ``other``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(QQ, [1,2,3], 1)
            sage: a[:] = 2, 0, -3
            sage: s = a.__floordiv__(2)
            sage: a //= 2; a
            1-index numpy components w.r.t. (1, 2, 3)
            sage: a[:]
            array([ 1.,  0., -2.])
            sage: a == s
            True
            sage: b = ComponentNumpy(QQ, [1,2,3], 1)
            sage: a[:] = 2, 0, -3
            sage: b[:] = 1, 4, -2
            sage: s = a.__floordiv__(b)
            sage: a //= b; a
            1-index numpy components w.r.t. (1, 2, 3)
            sage: a[:]
            array([2., 0., 1.])
            sage: a == s
            True

        """
        return self // other

    def __ipow__(self, other):
        r"""
        In-place operator for evaluating self_i ^ other_i for each element of Component
        ``self`` with the respective element of ``other``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(QQ, [1,2,3], 1)
            sage: a[:] = 3, 0, -3
            sage: s = a.__pow__(2)
            sage: a ^= 2; a
            1-index numpy components w.r.t. (1, 2, 3)
            sage: a[:]
            array([9., 0., 9.])
            sage: a == s
            True
            sage: b = ComponentNumpy(QQ, [1,2,3], 1)
            sage: a[:] = 3, 0, -3
            sage: b[:] = -1, 4, -2
            sage: s = a.__pow__(b)
            sage: a ^= b; a[:]
            array([0.33333333, 0.        , 0.11111111])
            sage: a == s
            True

        """
        return self ** other

    def __imod__(self, other):
        r"""
        In-place operator for evaluating self_i % other_i for each element of Component
        ``self`` with the respective element of ``other``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: a = ComponentNumpy(QQ, [1,2,3], 1)
            sage: a[:] = 1, 2, -3
            sage: s = a.__mod__(2)
            sage: a %= 2; a
            1-index numpy components w.r.t. (1, 2, 3)
            sage: a[:]
            array([1., 0., 1.])
            sage: a == s
            True
            sage: b = ComponentNumpy(QQ, [1,2,3], 1)
            sage: a[:] = 1, 2, -3
            sage: b[:] = 1, 4, -2
            sage: s = a.__mod__(b)
            sage: a %= b
            sage: a[:]
            array([ 0.,  2., -1.])
            sage: a == s
            True

        """
        return self % other

    def trace(self, pos1, pos2):
        r"""
        Index contraction.

        INPUT:

        - ``pos1`` -- position of the first index for the contraction (with the
          convention position=0 for the first slot)
        - ``pos2`` -- position of the second index for the contraction

        OUTPUT:

        - set of components resulting from the (pos1, pos2) contraction

        EXAMPLES:

        Self-contraction of a set of components with 2 indices::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: V = VectorSpace(QQ, 3)
            sage: c = ComponentNumpy(QQ, V.basis(), 2)
            sage: c[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: c.trace(0,1)
            15.0
            sage: c[0,0] + c[1,1] + c[2,2]  # check
            15.0

        """
        if self._nid < 2:
            raise ValueError("contraction can be performed only on " +
                             "components with at least 2 indices")
        if pos1 < 0 or pos1 > self._nid - 1:
            raise IndexError("pos1 out of range")
        if pos2 < 0 or pos2 > self._nid - 1:
            raise IndexError("pos2 out of range")
        if pos1 == pos2:
            raise IndexError("the two positions must differ for the " +
                             "contraction to be meaningful")

        comp = self._comp
        if self._nid == 2:
            ret = 0
            for i in range(self._shape[0]):
                ret += comp[i,i]
        else:
            ret = ComponentNumpy(self._ring, self._frame, self._nid - 2,
                                    self._sindex, self._output_formatter)
            ret._comp = np.trace(self._comp, axis1=pos1, axis2=pos2)

        return ret

    def _khatri_rao_product(*args):
        r"""
        Compute the Khatri-Rao product of input numpy-arrays.

        EXAMPLES:

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: import numpy as np
            sage: a = np.array([[1, 2],[3, 4]])
            sage: b = np.array([[5, 6],[7, 8]])
            sage: c = np.array([[9, 10], [11, 12]])
            sage: s = ComponentNumpy._khatri_rao_product(a, b, c); s
            array([[ 45., 120.],
                   [ 55., 144.],
                   [ 63., 160.],
                   [ 77., 192.],
                   [135., 240.],
                   [165., 288.],
                   [189., 320.],
                   [231., 384.]])

        """
        cols = args[0].shape[1]
        rows = np.prod([arg.shape[0] for arg in args])
        ret = np.zeros((rows, cols))
        for i in range(cols):
            temp = args[0][:, i]
            for matrix in args[1:]:
                temp = np.einsum('i,j->ij', temp, matrix[:, i]).ravel()
            ret[:, i] = temp
        return ret

    def _svd(matrix, rank=None):
        r"""
        Compute svd of input numpy-matrix.

        EXAMPLES:

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: import numpy as np
            sage: a = np.array([[1, 2],[3, 4]])
            sage: u, s, v = ComponentNumpy._svd(a)
            sage: u
            array([[-0.40455358, -0.9145143 ],
                   [-0.9145143 ,  0.40455358]])
            sage: s
            array([5.4649857 , 0.36596619])
            sage: v
            array([[-0.57604844, -0.81741556],
                   [ 0.81741556, -0.57604844]])

        """
        import scipy
        dim_1, dim_2 = matrix.shape
        if dim_1 <= dim_2:
            mdim = dim_1
        else:
            mdim = dim_2
        if rank is None or rank >= mdim:
            U, S, V = scipy.linalg.svd(matrix)
            U, S, V = U[:, :rank], S[:rank], V[:rank, :]
            return U, S, V
        else:
            if dim_1 < dim_2:
                S, U = scipy.sparse.linalg.eigsh(np.dot(matrix, matrix.T), k=rank, which='LM')
                S = np.sqrt(S)
                V = np.dot(matrix.T, U * 1 / S[None, :])
            else:
                S, V = scipy.sparse.linalg.eigsh(np.dot(matrix.T, matrix), k=rank, which='LM')
                S = np.sqrt(S)
                U = np.dot(matrix, V) * 1 / S[None, :]

            U, S, V = U[:, ::-1], S[::-1], V[:, ::-1]
            return U, S, V.T

    def CPD(self, rank, iterations=1000, epsilon=10e-5, tolerance=10e-8, factor_matrix=None, algo='svd', return_error=False, return_reconstruction=False):
        r"""
        Canonical polyadic decomposition (CPD) using ALS.

        INPUT:

        - ``rank`` -- Desired rank for the output tensor. This is the
          number of rank-one components the tensor is decomposed into.
        - ``iterations`` -- (default: ``1000``) The number of iterations
          to perform the ALS algorithm.
        - ``epsilon`` -- (default: ``10e-5``) Convergence criterion for
          the ALS algorithm. The algorithm stops if the relative change
          in the factor matrices between iterations is less than this value.
        - ``tolerance`` -- (default: ``10e-8``) Numerical tolerance for
          the ALS optimization process. It sets the threshold for small
          values to be considered as zero, improving numerical stability.
        - ``factor_matrix`` -- (default: ``None``) Initial factor matrices
          for the decomposition. If ``None``, they are initialized randomly
          or via the specified ``algo`` method.
        - ``algo`` -- (default: ``svd``) Initialization algorithm for
          factor matrices. Common options include:
            ``'svd'``: Singular Value Decomposition for a better starting point.
        -   ``'random'``: Random initialization of factor matrices.
        - ``return_error`` -- (default: ``False``) If ``True``, returns the
          reconstruction error after the decomposition.
        - ``return_reconstruction`` -- (default: ``False``) If ``True``,
          returns the reconstructed tensor from the factor matrices.

        OUTPUT:

        - If ``return_error`` is ``True``, returns a tuple with the factor
          matrices and the reconstruction error.
        - If ``return_reconstruction`` is ``True``, returns a tuple with the
          factor matrices and the reconstructed tensor.
        - Otherwise, returns only the factor matrices as a list, where each
          element is a factor matrix corresponding to one mode of the tensor.

        """
        comp = self._comp
        if np.all(comp == 0):
            raise ValueError("all elements in component are zero")
        if factor_matrix is None:
            fmat = [np.array([]) for _ in range(self._nid)]

            if (np.array(self._shape) >= rank).sum() == self._nid and algo == 'svd':
                for mode in range(self._nid):
                    k = np.reshape(np.moveaxis(comp, mode, 0), (comp.shape[mode], -1))

                    fmat[mode], _, _ = ComponentNumpy._svd(k, rank)
            else:
                fmat = [np.random.randn(mode_size, rank) for mode_size in self._shape]
        else:
            if not isinstance(factor_matrix, list):
                raise TypeError("factor_matrix must be of type list")
            if not all(isinstance(m, np.ndarray) for m in factor_matrix):
                raise TypeError("elements of factor matrix must be of type np.ndarray")

            if len(factor_matrix) != self._nid:
                raise ValueError("Invalid length of factor_matrix, expecting {}, got {}".format(self._nid, len(factor_matrix)))
            if not all(m.shape == (mode, rank[0]) for m, mode in zip(factor_matrix, self._shape)):
                raise ValueError("factor matrix with incorrect shape passed")
            fmat = factor_matrix.copy()

        from functools import reduce
        initial_core = np.repeat(np.array([1]), rank)
        norm = np.linalg.norm(comp.data)
        new_cost = 0
        for _ in range(iterations):
            for mode in range(self._nid):
                array = [fmat[i] for i in range(len(fmat)) if i != mode]
                khatri_rao = ComponentNumpy._khatri_rao_product(*array)
                hadamard = reduce(np.multiply, [np.dot(mat.T, mat) for i, mat in enumerate(fmat) if i != mode])

                soln = reduce(np.dot, [np.reshape(np.moveaxis(comp, mode, 0), (comp.shape[mode], -1)).data,
                                        khatri_rao,
                                        np.linalg.pinv(hadamard)])
                fmat[mode] = soln

            core_shape = (fmat[0].shape[1],) * len(fmat)
            order = len(core_shape)
            _rank = core_shape[0]

            t = np.zeros(core_shape)
            t[np.diag_indices(_rank, ndim=order)] = initial_core

            weights = np.array(initial_core)
            for _mode, f_mat in enumerate(fmat):
                core_shape = list(core_shape)
                core_shape[_mode] = f_mat.shape[0]
                shape = list(core_shape)
                mode_dim = shape.pop(_mode)
                shape.insert(0, mode_dim)
                t = np.moveaxis(np.reshape(np.dot(f_mat, np.reshape(np.moveaxis(t, _mode, 0), (t.shape[_mode], -1))), shape), 0, _mode)

            residual = comp - t
            previous_cost = new_cost
            new_cost = abs(np.linalg.norm(residual.data) / norm)

            converged = abs(new_cost - previous_cost) <= tolerance
            if new_cost <= epsilon:
                break
            if converged:
                break

        out = [weights, fmat]
        if return_reconstruction:
            out.append(t)
        if return_error:
            out.append(new_cost)
        return tuple(out)

    def TT(self, rank, return_error=False, return_reconstruction=False):
        r"""
        Tensor Train (TT) Decomposition.

        INPUT:

        - ``rank`` -- Desired rank for the decomposition. This rank defines
          the size of the intermediate tensor cores and controls the trade-off
          between the accuracy of the decomposition and the computational cost.
        - ``return_error`` -- (default: ``False``) If ``True``, returns the
          reconstruction error after performing the decomposition.
        - ``return_reconstruction`` -- (default: ``False``) If ``True``, returns
          the reconstructed tensor from the decomposed TT cores.

        OUTPUT:

        - If ``return_error`` is ``True``, returns a tuple containing the TT
          cores and the reconstruction error.
        - If ``return_reconstruction`` is ``True``, returns a tuple containing
          the TT cores and the reconstructed tensor.
        - Otherwise, returns only the TT cores as a list, where each element is
          a core tensor corresponding to one mode of the original tensor.

        """
        import scipy
        comp = self._comp
        if np.all(comp == 0):
            raise ValueError("all elements in component are zero")
        cores = []
        sizes = self._shape
        C = comp
        for k in range(self._nid - 1):
            rows = rank[k] * sizes[k]
            C = np.reshape(C, [rows, -1], order='F')
            U, S, V = scipy.linalg.svd(C)
            U = U[:, :rank[k + 1]]
            S = S[:rank[k + 1]]
            V = V[:rank[k + 1], :].T

            if k == 0:
                new_core = np.reshape(U, [sizes[k], rank[k+1]], order='F')
            else:
                new_core = np.reshape(U, [rank[k], sizes[k], rank[k+1]], order='F')

            cores.append(new_core)
            C = np.dot(V, np.diag(S)).T
        new_core = C
        cores.append(new_core)

        if return_error or return_reconstruction:
            res = [core.copy() for core in cores]
            rank = tuple(core_values.shape[-1] for core_values in res[:-1]) + (1,)
            core = res[0]
            data = core

            shape = [None] * len(res)
            shape[0] = res[0].shape[0]
            shape[-1] = res[-1].shape[1]
            for i in range(1, len(res) - 1):
                shape[i] = res[i].shape[1]

            for i, core in enumerate(res[1:]):
                shape_2d = [rank[i], rank[i+1] * shape[i + 1]]
                core_flat = np.reshape(core, shape_2d, order='F')
                data = np.reshape(data, [-1, rank[i]], order='F')
                data = np.dot(data, core_flat)
            data = np.reshape(data, shape, order='F')

        out = list([cores])

        if return_reconstruction:
            out.append(data)
        if return_error:
            residual = comp - data
            cost = abs(np.linalg.norm(residual.data) / np.linalg.norm(comp.data))
            out.append(cost)
        return out

    def Tucker(self, rank, process=(), return_error=False, return_reconstruction=False):
        r"""
        Tucker Decomposition.

        The Tucker decomposition, also known as Tucker decomposition or higher-order singular value decomposition (HOSVD), is a method for decomposing a tensor into a core tensor multiplied by a matrix along each mode. This decomposition generalizes the matrix singular value decomposition (SVD) to higher-order tensors.

        INPUT:

        - ``rank`` -- Desired rank for the decomposition. This can be  a tuple/list
          specifying the rank for each mode separately.
        - ``process`` -- (default: ``()``) A tuple or list of mode indices to be
          processed. If empty, all modes are processed. This allows for partial
          decomposition along specific modes if needed.
        - ``return_error`` -- (default: ``False``) If ``True``, returns the
          reconstruction error after performing the decomposition.
        - ``return_reconstruction`` -- (default: ``False``) If ``True``, returns
          the reconstructed tensor from the decomposed core tensor and factor matrices.

        OUTPUT:

        - If ``return_error`` is ``True``, returns a tuple containing the core tensor,
          factor matrices, and the reconstruction error.
        - If ``return_reconstruction`` is ``True``, returns a tuple containing the
          core tensor, factor matrices, and the reconstructed tensor.
        - Otherwise, returns a tuple containing the core tensor and the factor matrices.
          The core tensor is a lower-dimensional representation of the original tensor,
          and the factor matrices represent the mappings along each mode.

    """
        comp = self._comp
        if np.all(comp == 0):
            raise ValueError("all elements in component are zero")
        fmat = [np.array([])] * self._nid
        core = comp.copy()

        if not process:
            process = tuple(range(self._nid))
        for mode in range(self._nid):
            if mode not in process:
                fmat[mode] = np.eye(self._shape[mode])
                continue
            tensor_unfolded = np.reshape(np.moveaxis(comp, mode, 0), (comp.shape[mode], -1))

            U, _, _, = ComponentNumpy._svd(tensor_unfolded, rank[mode])
            fmat[mode] = U

            new_shape = list(core.shape)
            new_shape[mode] = U.T.shape[0]
            full_shape = list(new_shape)
            full_shape.insert(0, full_shape.pop(mode))
            core = np.moveaxis(np.reshape(np.dot(U.T, np.reshape(np.moveaxis(core, mode, 0), (core.shape[mode], -1))), full_shape), 0, mode)

        out = (core, fmat)
        if return_error or return_reconstruction:
            t = core
            res = [mat.copy() for mat in fmat]
            for mode, fmat in enumerate(res):
                orig_shape = list(t.shape)
                new_shape = orig_shape
                new_shape[mode] = fmat.shape[0]
                full_shape = list(new_shape)
                mode_dim = full_shape.pop(mode)
                full_shape.insert(0, mode_dim)
                t = np.moveaxis(np.reshape(np.dot(fmat, np.reshape(np.moveaxis(t, mode, 0), (t.shape[mode], -1))), full_shape), 0, mode)
        if return_reconstruction:
            out += (t,)
        if return_error:
            residual = comp - t
            error = abs(np.linalg.norm(residual.data) / np.linalg.norm(comp.data))
            out += (error,)
        return out

class CompNumpyWithSym(ComponentNumpy):

    def __init__(self, ring, frame, nb_indices, shape=None, start_index=0,
                 output_formatter=None, sym=None, antisym=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp_numpy import CompNumpyWithSym
            sage: c = CompNumpyWithSym(ZZ, [1,2,3], 4, sym=(0,1), antisym=(2,3)); c
            4-indices numpy components w.r.t. (1, 2, 3)

        """
        ComponentNumpy.__init__(self, ring, frame, nb_indices, shape, start_index, output_formatter)
        if any(self._shape[0] != dim for dim in self._shape):
            raise KeyError("symmetry can only be defined for cubic tensors")
        from .comp import CompFullySym
        self._sym, self._antisym = CompFullySym._canonicalize_sym_antisym(
            nb_indices, sym, antisym)

    def _new_instance(self):
        r"""
        Creates a :class:`CompNumpyWithSym` instance of the same number of indices
        and w.r.t. the same frame.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import CompNumpyWithSym
            sage: c = CompNumpyWithSym(ZZ, [1,2,3], 2)
            sage: c._new_instance()
            2-indices numpy components w.r.t. (1, 2, 3)

        """
        return self.__class__(self._ring, self._frame, self._nid, self._shape,
                          self._sindex, self._output_formatter, self._sym, self._antisym)

    def __setitem__(self, args, value):
        r"""
        Sets the component corresponding to the given indices.

        INPUT:

        - ``args`` -- list of indices (possibly a single integer if
          self is a 1-index object); if ``[:]`` is provided, all the
          components are set
        - ``value`` -- the value to be set or a list of values if
          ``args = [:]`` (``slice(None)``)

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import CompNumpyWithSym
            sage: c = CompNumpyWithSym(ZZ, [1,2,3], 2, sym=(0,1))
            sage: c.__setitem__((0,1), -4)
            sage: c[:]
            array([[ 0., -4.,  0.],
                   [-4.,  0.,  0.],
                   [ 0.,  0.,  0.]])
            sage: c[0,1] = 4
            sage: c[:]
            array([[0., 4., 0.],
                   [4., 0., 0.],
                   [0., 0., 0.]])

        """
        if isinstance(args, list):  # case of [[...]] syntax
            if isinstance(args[0], slice):
                indices = args[0]
            elif isinstance(args[0], (tuple, list)): # to ensure equivalence between
                indices = args[0]           # [[(i,j,...)]] or [[[i,j,...]]] and [[i,j,...]]
            else:
                indices = tuple(args)
        else:
            # Determining from the input the list of indices and the format
            if isinstance(args, (int, Integer, slice)):
                indices = args
            elif isinstance(args[0], slice):
                indices = args[0]
            elif len(args) == self._nid:
                indices = args

        if isinstance(indices, slice):
            start, stop = indices.start, indices.stop
            range_start, range_end = self._sindex[0], self._sindex[0] + self._shape[0]

            if indices.start is not None:
                start = indices.start - range_start
            else:
                start = 0
            if indices.stop is not None:
                stop = indices.stop - range_start
            else:
                stop = range_end - range_start
            if not ((0 <= start <= range_end - range_start - 1) and
                    (0 <= stop <= range_end - range_start)):
                raise IndexError("[start:stop] not in range [{},{}]"
                    .format(range_start, range_end))
            #TODO: defining sym tensor for non-sym tensor input
            self._comp[start:stop:indices.step] = value
        else:
            indices = self._check_indices(indices)
            self._comp[indices] = value
            from itertools import permutations
            permuted_indices = [list(indices)]
            for dims in self._sym:
                symmetric_permutations = set(permutations([indices[i] for i in dims]))
                for perm in symmetric_permutations:
                    new_indices = list(indices)
                    for i, dim in enumerate(dims):
                        new_indices[dim] = perm[i]
                    permuted_indices.append(new_indices)
                    self._comp[tuple(new_indices)] = value

            for ind_1, ind_2 in self._antisym:
                ind = list(indices).copy()
                ind[ind_1], ind[ind_2] = ind[ind_2], ind[ind_1]
                self._comp[tuple(ind)] = - value

class CompNumpyFullySym(CompNumpyWithSym):

    def __init__(self, ring, frame, nb_indices, start_index=0,
                 output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp_numpy import CompNumpyFullySym
            sage: c = CompNumpyFullySym(ZZ, [1,2,3], 4, sym=(0,1), antisym=(2,3)); c
            4-indices numpy components w.r.t. (1, 2, 3)

        """
        CompNumpyWithSym.__init__(self, ring, frame, nb_indices, start_index,
                             output_formatter, sym=range(nb_indices))

    def _new_instance(self):
        r"""
        Creates a :class:`CompNumpyFullySym` instance of the same number of indices
        and w.r.t. the same frame.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import CompNumpyFullySym
            sage: c = CompNumpyFullySym(ZZ, [1,2,3], 2)
            sage: c._new_instance()
            2-indices numpy components w.r.t. (1, 2, 3)

        """
        return self.__class__(self._ring, self._frame, self._nid, self._sindex,
                            self._output_formatter)

class CompNumpyFullyAntiSym(CompNumpyWithSym):

    def __init__(self, ring, frame, nb_indices, start_index=0,
                 output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp_numpy import CompNumpyFullyAntiSym
            sage: c = CompNumpyFullyAntiSym(ZZ, [1,2,3], 4, sym=(0,1), antisym=(2,3)); c
            4-indices numpy components w.r.t. (1, 2, 3)

        """
        CompNumpyWithSym.__init__(self, ring, frame, nb_indices, start_index,
                             output_formatter, antisym=range(nb_indices))

    def _new_instance(self):
        r"""
        Creates a :class:`CompNumpyFullyAntiSym` instance of the same number of indices
        and w.r.t. the same frame.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import CompNumpyFullyAntiSym
            sage: c = CompNumpyFullyAntiSym(ZZ, [1,2,3], 2)
            sage: c._new_instance()
            2-indices numpy components w.r.t. (1, 2, 3)

        """
        return self.__class__(self._ring, self._frame, self._nid, self._sindex,
                                self._output_formatter)

class KroneckerDeltaNumpy(CompNumpyFullySym):
    def __init__(self, ring, frame, start_index=0, output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp_numpy import KroneckerDeltaNumpy
            sage: d = KroneckerDeltaNumpy(ZZ, (1,2,3)); d

        """
        CompNumpyFullySym.__init__(self, ring, frame, 2, start_index,
                              output_formatter)
        for i in range(self._sindex, self._dim + self._sindex):
            self._comp[(i,i)] = self._ring(1)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import KroneckerDeltaNumpy
            sage: KroneckerDeltaNumpy(ZZ, (1,2,3))
            Kronecker delta of size 3x3

        """
        n = str(self._dim)
        return "Kronecker delta of size " + n + "x" + n

    def __setitem__(self, args, value):
        r"""
        Should not be used (the components of a Kronecker delta are constant)

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import KroneckerDeltaNumpy
            sage: d = KroneckerDeltaNumpy(ZZ, (1,2,3))
            sage: d.__setitem__((0,0), 1)
            Traceback (most recent call last):
            ...
            TypeError: the components of a Kronecker delta cannot be changed

        """
        raise TypeError("the components of a Kronecker delta cannot be changed")
