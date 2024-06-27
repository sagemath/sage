from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.tensor.modules.comp import Components
import numpy as np


class ComponentNumpy(SageObject):
    def __init__(self, ring, frame, nb_indices, shape=None, start_index=0,
                 output_formatter=None):
        self._ring = ring
        self._frame = tuple(frame) if all(isinstance(i, (list, tuple)) for i in frame)\
              else (frame,) * nb_indices
        self._nid = nb_indices
        self._shape = (len(self._frame[0]),) * nb_indices if shape is None else tuple(shape)
        self._sindex = tuple(start_index) if isinstance(start_index, (list, tuple))\
              else (start_index,) * nb_indices
        self._output_formatter = output_formatter
        self._comp = np.zeros(shape=self._shape, dtype=np.float64)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp_numpy import ComponentNumpy
            sage: c = ComponentNumpy(ZZ, [1,2,3], 2)
            sage: c._repr_()
            '2-indices components w.r.t. [1, 2, 3]'

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: CompWithSym(ZZ, [1,2,3], 4, sym=(0,1))
            4-indices components w.r.t. [1, 2, 3],
             with symmetry on the index positions (0, 1)
            sage: CompWithSym(ZZ, [1,2,3], 4, sym=(0,1), antisym=(2,3))
            4-indices components w.r.t. [1, 2, 3],
             with symmetry on the index positions (0, 1),
             with antisymmetry on the index positions (2, 3)

            sage: from sage.tensor.modules.comp import CompFullySym
            sage: CompFullySym(ZZ, (1,2,3), 4)
            Fully symmetric 4-indices components w.r.t. (1, 2, 3)

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: CompFullyAntiSym(ZZ, (1,2,3), 4)
            Fully antisymmetric 4-indices components w.r.t. (1, 2, 3)
        """
        description = str()
        if all(self._shape[0] == dim for dim in self._shape):
            description += "{}-shaped".format(self._shape)
        description += " numpy"
        description += " components w.r.t. " + str(self._frame)
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
            2-indices components w.r.t. [1, 2, 3]

        """
        return self.__class__(self._ring, self._frame, self._nid, self._shape,
                          self._sindex, self._output_formatter)

    def copy(self):
        r"""
        Return an exact copy of ``self``.

        EXAMPLES:

        Copy of a set of components with a single index::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ,3)
            sage: a = ComponentNumpy(QQ, V.basis(), 1)
            sage: a[:] = -2, 1, 5
            sage: b = a.copy() ; b
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: b[:]
            [-2, 1, 5]
            sage: b == a
            True
            sage: b is a  # b is a distinct object
            False

        """
        result = self._new_instance()
        result._comp = np.copy(self._comp)
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
            sage: c = ComponentNumpy(ZZ, [1,2,3], 2)
            sage: c[1,2]    # unset components are zero
            0
            sage: c.__getitem__((1,2))
            0
            sage: c.__getitem__([1,2])
            0
            sage: c[1,2] = -4
            sage: c[1,2]
            -4
            sage: c.__getitem__((1,2))
            -4
            sage: c[:]
            [ 0  0  0]
            [ 0  0 -4]
            [ 0  0  0]
            sage: c.__getitem__(slice(None))
            [ 0  0  0]
            [ 0  0 -4]
            [ 0  0  0]

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
                if len(args) == 2:
                    format_type = args[1]
            elif len(args) == self._nid:
                indices = args
            else:
                format_type = args[-1]
                indices = args[:-1]
        if isinstance(indices, slice):
            indices = slice(indices.start - self._sindex[0], indices.stop - self._sindex[0])
            return self._comp[indices] # to be implemented ## output formatter
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

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c.__setitem__((0,1), -4)
            sage: c[:]
            [ 0 -4  0]
            [ 0  0  0]
            [ 0  0  0]
            sage: c[0,1] = -4
            sage: c[:]
            [ 0 -4  0]
            [ 0  0  0]
            [ 0  0  0]
            sage: c.__setitem__(slice(None), [[0, 1, 2], [3, 4, 5], [6, 7, 8]])
            sage: c[:]
            [0 1 2]
            [3 4 5]
            [6 7 8]

        """
        format_type = None # default value, possibly redefined below
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
                if len(args) == 2:
                    format_type = args[1]
            elif len(args) == self._nid:
                indices = args
            else:
                format_type = args[-1]
                indices = args[:-1]
        if isinstance(indices, slice):
            self._set_list(indices, format_type, value) # to be implemented ## _setlist
        else:
            ind = self._check_indices(indices)
            self._comp[ind] = value

    def __pos__(self):
        return self.copy()

    def __neg__(self):
        neg_com = self._new_instance()
        neg_com._comp = - np.copy(self._comp)
    
    def __add__(self, other):
        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        if not isinstance(other, (Components, ComponentNumpy)):
            raise TypeError("the second argument for the addition must be " +
                            "an instance of Components")
        if isinstance(other, CompWithSym):
            return other + self     # to deal properly with symmetries
        if other._frame != self._frame:
            raise ValueError("the two sets of components are not defined on " +
                             "the same frame")
        if other._nid != self._nid:
            raise ValueError("the two sets of components do not have the " +
                             "same number of indices")
        if other._sindex != self._sindex:
            raise ValueError("the two sets of components do not have the " +
                             "same starting index")