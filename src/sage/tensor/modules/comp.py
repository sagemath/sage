r"""
Components as indexed sets of ring elements

The class :class:`Components` is a technical class to take in charge the
storage and manipulation of **indexed elements of a commutative ring** that
represent the components of some "mathematical entity" with respect to some
"frame". Examples of *entity/frame* are *vector/vector-space basis* or
*vector field/vector frame on some manifold*. More generally, the components
can be those of a tensor on a free module or those of a tensor field on a
manifold. They can also be non-tensorial quantities, like connection
coefficients or structure coefficients of a vector frame.

The individual components are assumed to belong to a given commutative ring
and are labelled by *indices*, which are *tuples of integers*.
The following operations are implemented on components with respect
to a given frame:

* arithmetics (addition, subtraction, multiplication by a ring element)
* handling of symmetries or antisymmetries on the indices
* symmetrization and antisymmetrization
* tensor product
* contraction

Various subclasses of class :class:`Components` are

* :class:`CompWithSym` for components with symmetries or antisymmetries w.r.t.
  index permutations

  * :class:`CompFullySym` for fully symmetric components w.r.t. index
    permutations

    * :class:`KroneckerDelta` for the Kronecker delta symbol

  * :class:`CompFullyAntiSym` for fully antisymmetric components w.r.t. index
    permutations

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Joris Vankerschaver (2010): for the idea of storing only the non-zero
  components as dictionaries, whose keys are the component indices (implemented
  in the old class ``DifferentialForm``; see :issue:`24444`)
- Marco Mancini (2015) : parallelization of some computations

EXAMPLES:

Set of components with 2 indices on a 3-dimensional vector space, the frame
being some basis of the vector space::

    sage: from sage.tensor.modules.comp import Components
    sage: V = VectorSpace(QQ,3)
    sage: basis = V.basis() ; basis
    [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    sage: c = Components(QQ, basis, 2) ; c
    2-indices components w.r.t. [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1)
    ]

Actually, the frame can be any object that has some length, i.e. on which
the function :func:`len()` can be called::

    sage: basis1 = V.gens() ; basis1
    ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    sage: c1 = Components(QQ, basis1, 2) ; c1
    2-indices components w.r.t. ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    sage: basis2 = ['a', 'b' , 'c']
    sage: c2 = Components(QQ, basis2, 2) ; c2
    2-indices components w.r.t. ['a', 'b', 'c']

A just created set of components is initialized to zero::

    sage: c.is_zero()
    True
    sage: c == 0
    True

This can also be checked on the list of components, which is returned by
the operator ``[:]``::

    sage: c[:]
    [0 0 0]
    [0 0 0]
    [0 0 0]

Individual components are accessed by providing their indices inside
square brackets::

    sage: c[1,2] = -3
    sage: c[:]
    [ 0  0  0]
    [ 0  0 -3]
    [ 0  0  0]
    sage: v = Components(QQ, basis, 1)
    sage: v[:]
    [0, 0, 0]
    sage: v[0]
    0
    sage: v[:] = (-1,3,2)
    sage: v[:]
    [-1, 3, 2]
    sage: v[0]
    -1

Sets of components with 2 indices can be converted into a matrix::

    sage: matrix(c)
    [ 0  0  0]
    [ 0  0 -3]
    [ 0  0  0]
    sage: matrix(c).parent()
    Full MatrixSpace of 3 by 3 dense matrices over Rational Field

By default, the indices range from `0` to `n-1`, where `n` is the length
of the frame. This can be changed via the argument ``start_index`` in
the :class:`Components` constructor::

    sage: v1 = Components(QQ, basis, 1, start_index=1)
    sage: v1[:]
    [0, 0, 0]
    sage: v1[0]
    Traceback (most recent call last):
    ...
    IndexError: index out of range: 0 not in [1, 3]
    sage: v1[1]
    0
    sage: v1[:] = v[:]  # list copy of all components
    sage: v1[:]
    [-1, 3, 2]
    sage: v1[1], v1[2], v1[3]
    (-1, 3, 2)
    sage: v[0], v[1], v[2]
    (-1, 3, 2)

If some formatter function or unbound method is provided via the argument
``output_formatter`` in the :class:`Components` constructor, it is used to
change the output of the access operator ``[...]``::

    sage: a = Components(QQ, basis, 2, output_formatter=Rational.numerical_approx)
    sage: a[1,2] = 1/3
    sage: a[1,2]
    0.333333333333333

The format can be passed to the formatter as the last argument of the
access operator ``[...]``::

    sage: a[1,2,10] # here the format is 10, for 10 bits of precision
    0.33
    sage: a[1,2,100]
    0.33333333333333333333333333333

The raw (unformatted) components are then accessed by the double bracket
operator::

    sage: a[[1,2]]
    1/3

For sets of components declared without any output formatter, there is no
difference between ``[...]`` and ``[[...]]``::

    sage: c[1,2] = 1/3
    sage: c[1,2], c[[1,2]]
    (1/3, 1/3)

The formatter is also used for the complete list of components::

    sage: a[:]
    [0.000000000000000 0.000000000000000 0.000000000000000]
    [0.000000000000000 0.000000000000000 0.333333333333333]
    [0.000000000000000 0.000000000000000 0.000000000000000]
    sage: a[:,10] # with a format different from the default one (53 bits)
    [0.00 0.00 0.00]
    [0.00 0.00 0.33]
    [0.00 0.00 0.00]

The complete list of components in raw form can be recovered by the double
bracket operator, replacing ``:`` by ``slice(None)`` (since ``a[[:]]``
generates a Python syntax error)::

    sage: a[[slice(None)]]
    [  0   0   0]
    [  0   0 1/3]
    [  0   0   0]

Another example of formatter: the Python built-in function :func:`str`
to generate string outputs::

    sage: b = Components(QQ, V.basis(), 1, output_formatter=str)
    sage: b[:] = (1, 0, -4)
    sage: b[:]
    ['1', '0', '-4']

For such a formatter, 2-indices components are no longer displayed as a
matrix::

    sage: b = Components(QQ, basis, 2, output_formatter=str)
    sage: b[0,1] = 1/3
    sage: b[:]
    [['0', '1/3', '0'], ['0', '0', '0'], ['0', '0', '0']]

But unformatted outputs still are::

    sage: b[[slice(None)]]
    [  0 1/3   0]
    [  0   0   0]
    [  0   0   0]

Internally, the components are stored as a dictionary (:attr:`_comp`) whose
keys are the indices; only the nonzero components are stored::

    sage: a[:]
    [0.000000000000000 0.000000000000000 0.000000000000000]
    [0.000000000000000 0.000000000000000 0.333333333333333]
    [0.000000000000000 0.000000000000000 0.000000000000000]
    sage: a._comp
    {(1, 2): 1/3}
    sage: v[:] = (-1, 0, 3)
    sage: v._comp  # random output order of the component dictionary
    {(0,): -1, (2,): 3}

In case of symmetries, only non-redundant components are stored::

    sage: from sage.tensor.modules.comp import CompFullyAntiSym
    sage: c = CompFullyAntiSym(QQ, basis, 2)
    sage: c[0,1] = 3
    sage: c[:]
    [ 0  3  0]
    [-3  0  0]
    [ 0  0  0]
    sage: c._comp
    {(0, 1): 3}
"""

# *****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2015 Marco Mancini <marco.mancini@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.parallel.decorate import parallel
from sage.parallel.parallelism import Parallelism
from operator import itemgetter


class Components(SageObject):
    r"""
    Indexed set of ring elements forming some components with respect
    to a given "frame".

    The "frame" can be a basis of some vector space or a vector frame on some
    manifold (i.e. a field of bases).
    The stored quantities can be tensor components or non-tensorial quantities,
    such as connection coefficients or structure coefficients. The symmetries
    over some indices are dealt by subclasses of the class :class:`Components`.

    INPUT:

    - ``ring`` -- commutative ring in which each component takes its value
    - ``frame`` -- frame with respect to which the components are defined;
      whatever type ``frame`` is, it should have a method ``__len__()``
      implemented, so that ``len(frame)`` returns the dimension, i.e. the size
      of a single index range
    - ``nb_indices`` -- number of integer indices labeling the components
    - ``start_index`` -- (default: 0) first value of a single index;
      accordingly a component index i must obey
      ``start_index <= i <= start_index + dim - 1``, where ``dim = len(frame)``.
    - ``output_formatter`` -- (default: ``None``) function or unbound
      method called to format the output of the component access
      operator ``[...]`` (method __getitem__); ``output_formatter`` must take
      1 or 2 arguments: the 1st argument must be an element of ``ring`` and
      the second one, if any, some format specification.

    EXAMPLES:

    Set of components with 2 indices on a 3-dimensional vector space, the frame
    being some basis of the vector space::

        sage: from sage.tensor.modules.comp import Components
        sage: V = VectorSpace(QQ,3)
        sage: basis = V.basis() ; basis
        [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        sage: c = Components(QQ, basis, 2) ; c
        2-indices components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]

    Actually, the frame can be any object that has some length, i.e. on which
    the function :func:`len()` can be called::

        sage: basis1 = V.gens() ; basis1
        ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        sage: c1 = Components(QQ, basis1, 2) ; c1
        2-indices components w.r.t. ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        sage: basis2 = ['a', 'b' , 'c']
        sage: c2 = Components(QQ, basis2, 2) ; c2
        2-indices components w.r.t. ['a', 'b', 'c']

    By default, the indices range from `0` to `n-1`, where `n` is the length
    of the frame. This can be changed via the argument ``start_index``::

        sage: c1 = Components(QQ, basis, 2, start_index=1)
        sage: c1[0,1]
        Traceback (most recent call last):
        ...
        IndexError: index out of range: 0 not in [1, 3]
        sage: c[0,1]  # for c, the index 0 is OK
        0
        sage: c[0,1] = -3
        sage: c1[:] = c[:] # list copy of all components
        sage: c1[1,2]  # (1,2) = (0,1) shifted by 1
        -3

    If some formatter function or unbound method is provided via the argument
    ``output_formatter``, it is used to change the output of the access
    operator ``[...]``::

        sage: a = Components(QQ, basis, 2, output_formatter=Rational.numerical_approx)
        sage: a[1,2] = 1/3
        sage: a[1,2]
        0.333333333333333

    The format can be passed to the formatter as the last argument of the
    access operator ``[...]``::

        sage: a[1,2,10] # here the format is 10, for 10 bits of precision
        0.33
        sage: a[1,2,100]
        0.33333333333333333333333333333

    The raw (unformatted) components are then accessed by the double bracket
    operator::

        sage: a[[1,2]]
        1/3

    For sets of components declared without any output formatter, there is no
    difference between ``[...]`` and ``[[...]]``::

        sage: c[1,2] = 1/3
        sage: c[1,2], c[[1,2]]
        (1/3, 1/3)

    The formatter is also used for the complete list of components::

        sage: a[:]
        [0.000000000000000 0.000000000000000 0.000000000000000]
        [0.000000000000000 0.000000000000000 0.333333333333333]
        [0.000000000000000 0.000000000000000 0.000000000000000]
        sage: a[:,10] # with a format different from the default one (53 bits)
        [0.00 0.00 0.00]
        [0.00 0.00 0.33]
        [0.00 0.00 0.00]

    The complete list of components in raw form can be recovered by the double
    bracket operator, replacing ``:`` by ``slice(None)`` (since ``a[[:]]``
    generates a Python syntax error)::

        sage: a[[slice(None)]]
        [  0   0   0]
        [  0   0 1/3]
        [  0   0   0]

    Another example of formatter: the Python built-in function :func:`str`
    to generate string outputs::

        sage: b = Components(QQ, V.basis(), 1, output_formatter=str)
        sage: b[:] = (1, 0, -4)
        sage: b[:]
        ['1', '0', '-4']

    For such a formatter, 2-indices components are no longer displayed as a
    matrix::

        sage: b = Components(QQ, basis, 2, output_formatter=str)
        sage: b[0,1] = 1/3
        sage: b[:]
        [['0', '1/3', '0'], ['0', '0', '0'], ['0', '0', '0']]

    But unformatted outputs still are::

        sage: b[[slice(None)]]
        [  0 1/3   0]
        [  0   0   0]
        [  0   0   0]

    Internally, the components are stored as a dictionary (:attr:`_comp`) whose
    keys are the indices; only the nonzero components are stored::

        sage: a[:]
        [0.000000000000000 0.000000000000000 0.000000000000000]
        [0.000000000000000 0.000000000000000 0.333333333333333]
        [0.000000000000000 0.000000000000000 0.000000000000000]
        sage: a._comp
        {(1, 2): 1/3}
        sage: v = Components(QQ, basis, 1)
        sage: v[:] = (-1, 0, 3)
        sage: v._comp  # random output order of the component dictionary
        {(0,): -1, (2,): 3}


    .. RUBRIC:: ARITHMETIC EXAMPLES:

    Unary plus operator::

        sage: a = Components(QQ, basis, 1)
        sage: a[:] = (-1, 0, 3)
        sage: s = +a ; s[:]
        [-1, 0, 3]
        sage: +a == a
        True

    Unary minus operator::

        sage: s = -a ; s[:]
        [1, 0, -3]

    Addition::

        sage: b = Components(QQ, basis, 1)
        sage: b[:] = (2, 1, 4)
        sage: s = a + b ; s[:]
        [1, 1, 7]
        sage: a + b == b + a
        True
        sage: a + (-a) == 0
        True

    Subtraction::

        sage: s = a - b ; s[:]
        [-3, -1, -1]
        sage: s + b == a
        True
        sage: a - b == - (b - a)
        True

    Multiplication by a scalar::

        sage: s = 2*a ; s[:]
        [-2, 0, 6]

    Division by a scalar::

        sage: s = a/2 ; s[:]
        [-1/2, 0, 3/2]
        sage: 2*(a/2) == a
        True

    Tensor product (by means of the operator ``*``)::

        sage: c = a*b ; c
        2-indices components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: a[:], b[:]
        ([-1, 0, 3], [2, 1, 4])
        sage: c[:]
        [-2 -1 -4]
        [ 0  0  0]
        [ 6  3 12]
        sage: d = c*a ; d
        3-indices components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: d[:]
        [[[2, 0, -6], [1, 0, -3], [4, 0, -12]],
         [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
         [[-6, 0, 18], [-3, 0, 9], [-12, 0, 36]]]
        sage: d[0,1,2] == a[0]*b[1]*a[2]
        True
    """
    def __init__(self, ring, frame, nb_indices, start_index=0,
                 output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp import Components
            sage: Components(ZZ, [1,2,3], 2)
            2-indices components w.r.t. [1, 2, 3]
        """
        # For efficiency, no test is performed regarding the type and range of
        # the arguments:
        self._ring = ring
        self._frame = frame
        self._nid = nb_indices
        self._dim = len(frame)
        self._sindex = start_index
        self._output_formatter = output_formatter
        self._comp = {} # the dictionary of components, with the index tuples
                        # as keys

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
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
        prefix, suffix = self._repr_symmetry()
        description = prefix
        description += str(self._nid)
        if self._nid == 1:
            description += "-index"
        else:
            description += "-indices"
        description += " components w.r.t. " + str(self._frame)
        description += suffix
        return description

    def _repr_symmetry(self):
        r"""
        Return a prefix and a suffix string describing the symmetry of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c._repr_symmetry()
            ('', '')
        """
        return "", ""

    def _new_instance(self):
        r"""
        Create a :class:`Components` instance of the same number of indices
        and w.r.t. the same frame.

        This method must be redefined by derived classes of
        :class:`Components`.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c._new_instance()
            2-indices components w.r.t. [1, 2, 3]
        """
        return Components(self._ring, self._frame, self._nid, self._sindex,
                          self._output_formatter)

    def copy(self):
        r"""
        Return an exact copy of ``self``.

        EXAMPLES:

        Copy of a set of components with a single index::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ,3)
            sage: a = Components(QQ, V.basis(), 1)
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
        for ind, val in self._comp.items():
            if isinstance(val, SageObject) and hasattr(val, 'copy'):
                result._comp[ind] = val.copy()
            else:
                result._comp[ind] = val
        return result

    def _del_zeros(self):
        r"""
        Delete all the zeros in the dictionary :attr:`_comp`.

        NB: The use case of this method must be rare because zeros are not
            stored in :attr:`_comp`.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c._comp = {(0,1): 3, (0,2): 0, (1,2): -5, (2,2): 0}  # enforcing zero storage
            sage: c._del_zeros()
            sage: c._comp
            {(0, 1): 3, (1, 2): -5}
        """
        # The zeros are first searched; they are deleted in a second stage, to
        # avoid changing the dictionary while it is read
        zeros = []
        for ind, value in self._comp.items():
            if value == 0:
                zeros.append(ind)
        for ind in zeros:
            del self._comp[ind]

    def _check_indices(self, indices):
        r"""
        Check the validity of a list of indices and returns a tuple from it.

        INPUT:

        - ``indices`` -- list of indices (possibly a single integer if
          ``self`` is a 1-index object)

        OUTPUT: a tuple containing valid indices

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
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
        """
        if isinstance(indices, (int, Integer)):
            ind = (indices,)
        else:
            ind = tuple(indices)
        if len(ind) != self._nid:
            raise ValueError(("wrong number of indices: {} expected,"
                             " while {} are provided").format(self._nid, len(ind)))
        si = self._sindex
        imax = self._dim - 1 + si
        for k in range(self._nid):
            i = ind[k]
            if i < si or i > imax:
                raise IndexError("index out of range: " +
                                 "{} not in [{}, {}]".format(i, si, imax))
        return ind

    def __getitem__(self, args):
        r"""
        Return the component corresponding to the given indices.

        INPUT:

        - ``args`` -- list of indices (possibly a single integer if
          ``self`` is a 1-index object) or the character ``:`` for the full list
          of components

        OUTPUT:

        - the component corresponding to ``args`` or, if ``args`` = ``:``,
          the full list of components, in the form ``T[i][j]...`` for the
          components `T_{ij...}` (for a 2-indices object, a matrix is returned)

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
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
            return self._get_list(indices, no_format, format_type)
        else:
            ind = self._check_indices(indices)
            if ind in self._comp:
                if no_format:
                    return self._comp[ind]
                elif format_type is None:
                    return self._output_formatter(self._comp[ind])
                else:
                    return self._output_formatter(self._comp[ind], format_type)
            else:  # if the value is not stored in self._comp, it is zero:
                if no_format:
                    return self._ring.zero()
                elif format_type is None:
                    return self._output_formatter(self._ring.zero())
                else:
                    return self._output_formatter(self._ring.zero(),
                                                 format_type)

    def _get_list(self, ind_slice, no_format=True, format_type=None):
        r"""
        Return the list of components (as nested list or matrix).

        INPUT:

        - ``ind_slice`` -- a slice object; unless the dimension is 1,
          this must be ``[:]``
        - ``no_format`` -- boolean (default: ``True``); determines whether some
          formatting of the components is to be performed
        - ``format_type`` -- (default: ``None``) argument to be passed
          to the formatting function ``self._output_formatter``, as the
          second (optional) argument

        OUTPUT:

        - general case: the nested list of components in the form
          ``T[i][j]...`` for the components `T_{ij...}`.

        - in the 1-dim case, a slice of that list if
          ``ind_slice = [a:b]``.

        - in the 2-dim case, a matrix (over the base ring of the components or
          of the formatted components if ``no_format`` is ``False``) is
          returned instead, except if the formatted components do not belong
          to any ring (for instance if they are strings).

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c[0,1], c[1,2] = 5, -4
            sage: c._get_list(slice(None))
            [ 0  5  0]
            [ 0  0 -4]
            [ 0  0  0]
            sage: v = Components(ZZ, [1,2,3], 1)
            sage: v[:] = 4, 5, 6
            sage: v._get_list(slice(None))
            [4, 5, 6]
            sage: v._get_list(slice(0,1))
            [4]
            sage: v._get_list(slice(0,2))
            [4, 5]
            sage: v._get_list(slice(2,3))
            [6]
        """
        si = self._sindex
        nsi = si + self._dim
        if self._nid == 1:
            if ind_slice.start is None:
                start = si
            else:
                start = ind_slice.start
            if ind_slice.stop is None:
                stop = nsi
            else:
                stop = ind_slice.stop
            if ind_slice.step is not None:
                raise NotImplementedError("function [start:stop:step] not implemented")
            if no_format:
                return [self[[i]] for i in range(start, stop)]
            else:
                return [self[i, format_type] for i in range(start, stop)]
        if ind_slice.start is not None or ind_slice.stop is not None:
            raise NotImplementedError("function [start:stop] not " +
                      "implemented for components with {} indices".format(self._nid))
        resu = [self._gen_list([i], no_format, format_type)
                for i in range(si, nsi)]
        if self._nid == 2:
            # 2-dim case: convert to matrix for a nicer output
            from sage.matrix.constructor import matrix
            from sage.structure.element import parent
            from sage.categories.rings import Rings
            if parent(resu[0][0]) in Rings():
                return matrix(resu)
        return resu

    def _gen_list(self, ind, no_format=True, format_type=None):
        r"""
        Recursive function to generate the list of values.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c[0,1], c[1,2] = 5, -4
            sage: c._gen_list([])
            [[0, 5, 0], [0, 0, -4], [0, 0, 0]]
            sage: c._gen_list([0])
            [0, 5, 0]
            sage: c._gen_list([1])
            [0, 0, -4]
            sage: c._gen_list([2])
            [0, 0, 0]
            sage: c._gen_list([0,1])
            5
        """
        if len(ind) == self._nid:
            if no_format:
                return self[ind]
            else:
                args = tuple(ind + [format_type])
                return self[args]
        else:
            si = self._sindex
            nsi = si + self._dim
            return [self._gen_list(ind + [i], no_format, format_type)
                    for i in range(si, nsi)]

    def __setitem__(self, args, value):
        r"""
        Set the component corresponding to the given indices.

        INPUT:

        - ``args`` -- list of indices (possibly a single integer if
          ``self`` is a 1-index object); if ``[:]`` is provided, all the
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
            self._set_list(indices, format_type, value)
        else:
            ind = self._check_indices(indices)
            # Check for a zero value
            #   The fast method is_trivial_zero() is employed preferably
            #   to the (possibly expensive) direct comparison to zero:
            if hasattr(value, 'is_trivial_zero'):
                zero_value = value.is_trivial_zero()
            else:
                zero_value = value == 0
            if zero_value:
                # if the component has been set previously, it is deleted,
                # otherwise nothing is done (zero components are not stored):
                if ind in self._comp:
                    del self._comp[ind]
            else:
                if format_type is None:
                    self._comp[ind] = self._ring(value)
                else:
                    self._comp[ind] = self._ring({format_type: value})
                    # NB: the writing
                    #   self._comp[ind] = self._ring(value, format_type)
                    # is not allowed when ring is an algebra and value some
                    # element of the algebra's base ring, cf. the discussion at
                    # https://github.com/sagemath/sage/issues/16054

    def _set_list(self, ind_slice, format_type, values):
        r"""
        Set the components from a list.

        INPUT:

        - ``ind_slice`` -- a slice object
        - ``format_type`` -- format possibly used to construct a ring element
        - ``values`` -- list of values for the components : the full list if
          ``ind_slice = [:]``, in the form ``T[i][j]...`` for the
          component `T_{ij...}`; in the 1-D case, ``ind_slice`` can be
          a slice of the full list, in the form  ``[a:b]``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c._set_list(slice(None), None, [[0, 1, 2], [3, 4, 5], [6, 7, 8]])
            sage: c[:]
            [0 1 2]
            [3 4 5]
            [6 7 8]
        """
        si = self._sindex
        nsi = si + self._dim
        if self._nid == 1:
            if ind_slice.start is None:
                start = si
            else:
                start = ind_slice.start
            if ind_slice.stop is None:
                stop = nsi
            else:
                stop = ind_slice.stop
            if ind_slice.step is not None:
                raise NotImplementedError("function [start:stop:step] not implemented")
            for i in range(start, stop):
                self[i, format_type] = values[i-start]
        else:
            if ind_slice.start is not None or ind_slice.stop is not None:
                raise NotImplementedError("function [start:stop] not " +
                      "implemented for components with {} indices".format(self._nid))
            for i in range(si, nsi):
                self._set_value_list([i], format_type, values[i-si])

    def _set_value_list(self, ind, format_type, val):
        r"""
        Recursive function to set a list of values to ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c._set_value_list([], None, [[1,2,3], [4,5,6], [7,8,9]])
            sage: c[:]
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: c._set_value_list([0], None, [-1,-2,-3])
            sage: c[:]
            [-1 -2 -3]
            [ 4  5  6]
            [ 7  8  9]
            sage: c._set_value_list([2,1], None, -8)
            sage: c[:]
            [-1 -2 -3]
            [ 4  5  6]
            [ 7 -8  9]
        """
        if len(ind) == self._nid:
            if format_type is not None:
                ind = tuple(ind + [format_type])
            self[ind] = val
        else:
            si = self._sindex
            nsi = si + self._dim
            for i in range(si, nsi):
                self._set_value_list(ind + [i], format_type, val[i-si])

    def items(self):
        r"""
        Return an iterable of ``(indices, value)`` elements.

        This may (but is not guaranteed to) suppress zero values.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components, CompWithSym

            sage: c = Components(ZZ, (ZZ^2).basis(), 3)
            sage: c[0,1,0], c[1,0,1], c[1,1,1] = -2, 5, 3
            sage: list(c.items())
            [((0, 1, 0), -2), ((1, 0, 1), 5), ((1, 1, 1), 3)]

            sage: c = CompWithSym(ZZ, (ZZ^2).basis(), 3, sym=((1, 2)))
            sage: c[0,1,0], c[1,0,1], c[1,1,1] = -2, 5, 3
            sage: list(c.items())
            [((0, 0, 1), -2),
            ((0, 1, 0), -2),
            ((1, 0, 1), 5),
            ((1, 1, 0), 5),
            ((1, 1, 1), 3)]
        """
        for ind in self.index_generator():
            val = self[ind]
            if hasattr(val, 'is_trivial_zero'):
                zero_value = val.is_trivial_zero()
            else:
                zero_value = val == 0
            if not zero_value:
                yield ind, val

    def display(self, symbol, latex_symbol=None, index_positions=None,
                index_labels=None, index_latex_labels=None,
                format_spec=None, only_nonzero=True, only_nonredundant=False):
        r"""
        Display all the components, one per line.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``symbol`` -- string (typically a single letter) specifying the
          symbol for the components
        - ``latex_symbol`` -- (default: ``None``) string specifying the LaTeX
          symbol for the components; if ``None``, ``symbol`` is used
        - ``index_positions`` -- (default: ``None``) string of length the
          number of indices of the components and composed of characters 'd'
          (for "down") or 'u' (for "up") to specify the position of each index:
          'd' corresponds to a subscript and 'u' to a superscript. If
          ``index_positions`` is ``None``, all indices are printed as
          subscripts
        - ``index_labels`` -- (default: ``None``) list of strings representing
          the labels of each of the individual indices within the index range
          defined at the construction of the object; if ``None``, integer
          labels are used
        - ``index_latex_labels`` -- (default: ``None``) list of strings
          representing the LaTeX labels of each of the individual indices
          within the index range defined at the construction of the object; if
          ``None``, integers labels are used
        - ``format_spec`` -- (default: ``None``) format specification passed
          to the output formatter declared at the construction of the object
        - ``only_nonzero`` -- boolean (default: ``True``); if ``True``, only
          nonzero components are displayed
        - ``only_nonredundant`` -- boolean (default: ``False``); if ``True``,
          only nonredundant components are displayed in case of symmetries

        EXAMPLES:

        Display of 3-indices components w.r.t. to the canonical basis of the
        free module `\ZZ^2` over the integer ring::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, (ZZ^2).basis(), 3)
            sage: c[0,1,0], c[1,0,1], c[1,1,1] = -2, 5, 3
            sage: c.display('c')
            c_010 = -2
            c_101 = 5
            c_111 = 3

        By default, only nonzero components are shown; to display all the
        components, it suffices to set the parameter ``only_nonzero`` to
        ``False``::

            sage: c.display('c', only_nonzero=False)
            c_000 = 0
            c_001 = 0
            c_010 = -2
            c_011 = 0
            c_100 = 0
            c_101 = 5
            c_110 = 0
            c_111 = 3

        By default, all indices are printed as subscripts, but any index
        position can be specified::

            sage: c.display('c', index_positions='udd')
            c^0_10 = -2
            c^1_01 = 5
            c^1_11 = 3
            sage: c.display('c', index_positions='udu')
            c^0_1^0 = -2
            c^1_0^1 = 5
            c^1_1^1 = 3
            sage: c.display('c', index_positions='ddu')
            c_01^0 = -2
            c_10^1 = 5
            c_11^1 = 3

        The LaTeX output is performed as an array, with the symbol adjustable
        if it differs from the text symbol::

            sage: latex(c.display('c', latex_symbol=r'\Gamma', index_positions='udd'))
            \begin{array}{lcl}
             \Gamma_{\phantom{\, 0}\,1\,0}^{\,0\phantom{\, 1}\phantom{\, 0}} & = & -2 \\
             \Gamma_{\phantom{\, 1}\,0\,1}^{\,1\phantom{\, 0}\phantom{\, 1}} & = & 5 \\
             \Gamma_{\phantom{\, 1}\,1\,1}^{\,1\phantom{\, 1}\phantom{\, 1}} & = & 3
            \end{array}

        The index labels can differ from integers::

            sage: c.display('c', index_labels=['x','y'])
            c_xyx = -2
            c_yxy = 5
            c_yyy = 3

        If the index labels are longer than a single character, they are
        separated by a comma::

            sage: c.display('c', index_labels=['r', 'th'])
            c_r,th,r = -2
            c_th,r,th = 5
            c_th,th,th = 3

        The LaTeX labels for the indices can be specified if they differ
        from the text ones::

            sage: c.display('c', index_labels=['r', 'th'],
            ....:           index_latex_labels=['r', r'\theta'])
            c_r,th,r = -2
            c_th,r,th = 5
            c_th,th,th = 3

        The display of components with symmetries is governed by the parameter
        ``only_nonredundant``::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: c = CompWithSym(ZZ, (ZZ^2).basis(), 3, sym=(1,2)) ; c
            3-indices components w.r.t. [
            (1, 0),
            (0, 1)
            ], with symmetry on the index positions (1, 2)
            sage: c[0,0,1] = 2
            sage: c.display('c')
            c_001 = 2
            c_010 = 2
            sage: c.display('c', only_nonredundant=True)
            c_001 = 2

        If some nontrivial output formatter has been set, the format can be
        specified by means of the argument ``format_spec``::

            sage: c = Components(QQ, (QQ^3).basis(), 2,
            ....:                output_formatter=Rational.numerical_approx)
            sage: c[0,1] = 1/3
            sage: c[2,1] = 2/7
            sage: c.display('C')  # default format (53 bits of precision)
            C_01 = 0.333333333333333
            C_21 = 0.285714285714286
            sage: c.display('C', format_spec=10)  # 10 bits of precision
            C_01 = 0.33
            C_21 = 0.29

        Check that the bug reported in :issue:`22520` is fixed::

            sage: c = Components(SR, [1, 2], 1)                                         # needs sage.symbolic
            sage: c[0] = SR.var('t', domain='real')                                     # needs sage.symbolic
            sage: c.display('c')                                                        # needs sage.symbolic
            c_0 = t
        """
        from sage.misc.latex import latex
        from sage.tensor.modules.format_utilities import FormattedExpansion
        si = self._sindex
        nsi = si + self._dim
        if latex_symbol is None:
            latex_symbol = symbol
        if index_positions is None:
            index_positions = self._nid * 'd'
        elif len(index_positions) != self._nid:
            raise ValueError("the argument 'index_positions' must contain " +
                             "{} characters".format(self._nid))
        if index_labels is None:
            index_labels = [str(i) for i in range(si, nsi)]
        elif len(index_labels) != self._dim:
            raise ValueError("the argument 'index_labels' must contain " +
                             "{} items".format(self._dim))
        # Index separator:
        max_len_symbols = max(len(s) for s in index_labels)
        if max_len_symbols == 1:
            sep = ''
        else:
            sep = ','
        if index_latex_labels is None:
            index_latex_labels = index_labels
        elif len(index_latex_labels) != self._dim:
            raise ValueError("the argument 'index_latex_labels' must " +
                             "contain {} items".format(self._dim))
        if only_nonredundant:
            # To simplify the implementation of the non-redundant
            # index generator, it generates indices in a different
            # order than the redundant index generator. For the
            # display, we sort the indices again.
            generator = list(self.non_redundant_index_generator())
            generator.sort()
        else:
            generator = self.index_generator()
        rtxt = ''
        rlatex = r'\begin{array}{lcl}'
        for ind in generator:
            ind_arg = ind + (format_spec,)
            val = self[ind_arg]
            # Check whether the value is zero, preferably via the
            # fast method is_trivial_zero():
            if hasattr(val, 'is_trivial_zero'):
                zero_value = val.is_trivial_zero()
            else:
                zero_value = val == 0
            if not zero_value or not only_nonzero:
                indices = ''  # text indices
                d_indices = '' # LaTeX down indices
                u_indices = '' # LaTeX up indices
                previous = None  # position of previous index
                for k in range(self._nid):
                    i = ind[k] - si
                    if index_positions[k] == 'd':
                        if previous == 'd':
                            indices += sep + index_labels[i]
                        else:
                            indices += '_' + index_labels[i]
                        d_indices += r'\,' + index_latex_labels[i]
                        u_indices += r'\phantom{{\, {}}}'.format(index_latex_labels[i])
                        previous = 'd'
                    else:
                        if previous == 'u':
                            indices += sep + index_labels[i]
                        else:
                            indices += '^' + index_labels[i]
                        d_indices += r'\phantom{{\, {}}}'.format(index_latex_labels[i])
                        u_indices += r'\,' + index_latex_labels[i]
                        previous = 'u'
                rtxt += symbol + indices + ' = {} \n'.format(val)
                rlatex += (latex_symbol + r'_{' + d_indices + r'}^{'
                           + u_indices + r'} & = & ' + latex(val) + r'\\')
        if rtxt == '':
            # no component has been displayed
            rlatex = ''
        else:
            # closing the display
            rtxt = rtxt[:-1]  # remove the last new line
            rlatex = rlatex[:-2] + r'\end{array}'
        return FormattedExpansion(rtxt, rlatex)

    def swap_adjacent_indices(self, pos1, pos2, pos3):
        r"""
        Swap two adjacent sets of indices.

        This method is essentially required to reorder the covariant and
        contravariant indices in the computation of a tensor product.

        INPUT:

        - ``pos1`` -- position of the first index of set 1 (with the convention
          ``position=0`` for the first slot)
        - ``pos2`` -- position of the first index of set 2 equals 1 plus the
          position of the last index of set 1 (since the two sets are adjacent)
        - ``pos3`` -- 1 plus position of the last index of set 2

        OUTPUT: components with index set 1 permuted with index set 2

        EXAMPLES:

        Swap of the two indices of a 2-indices set of components::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ, 3)
            sage: c = Components(QQ, V.basis(), 2)
            sage: c[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: c1 = c.swap_adjacent_indices(0,1,2)
            sage: c[:], c1[:]
            (
            [1 2 3]  [1 4 7]
            [4 5 6]  [2 5 8]
            [7 8 9], [3 6 9]
            )

        Swap of two pairs of indices on a 4-indices set of components::

            sage: d = c*c1 ; d
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: d1 = d.swap_adjacent_indices(0,2,4)
            sage: d[0,1,1,2]
            16
            sage: d1[1,2,0,1]
            16
            sage: d1[0,1,1,2]
            24
            sage: d[1,2,0,1]
            24
        """
        result = self._new_instance()
        for ind, val in self._comp.items():
            new_ind = ind[:pos1] + ind[pos2:pos3] + ind[pos1:pos2] + ind[pos3:]
            result._comp[new_ind] = val
            # the above writing is more efficient than result[new_ind] = val
            # it does not work for the derived class CompWithSym, but for the
            # latter, the function CompWithSym.swap_adjacent_indices will be
            # called and not the present function.
        return result

    def is_zero(self):
        r"""
        Return ``True`` if all the components are zero and ``False`` otherwise.

        EXAMPLES:

        A just-created set of components is initialized to zero::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ,3)
            sage: c = Components(QQ, V.basis(), 1)
            sage: c.is_zero()
            True
            sage: c[:]
            [0, 0, 0]
            sage: c[0] = 1 ; c[:]
            [1, 0, 0]
            sage: c.is_zero()
            False
            sage: c[0] = 0 ; c[:]
            [0, 0, 0]
            sage: c.is_zero()
            True

        It is equivalent to use the operator == to compare to zero::

            sage: c == 0
            True
            sage: c != 0
            False

        Comparing to a nonzero number is meaningless::

            sage: c == 1
            Traceback (most recent call last):
            ...
            TypeError: cannot compare a set of components to a number
        """
        if not self._comp:
            return True

        #!# What follows could be skipped since _comp should not contain
        # any zero value
        # In other words, the full method should be
        #   return self.comp == {}
        for val in self._comp.values():
            if not (val == 0):
                return False
        return True

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- set of components or 0

        OUTPUT: ``True`` if ``self`` is equal to ``other``,  or ``False`` otherwise

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 2)
            sage: c.__eq__(0)  # uninitialized components are zero
            True
            sage: c[0,1], c[1,2] = 5, -4
            sage: c.__eq__(0)
            False
            sage: c1 = Components(ZZ, [1,2,3], 2)
            sage: c1[0,1] = 5
            sage: c.__eq__(c1)
            False
            sage: c1[1,2] = -4
            sage: c.__eq__(c1)
            True
            sage: v = Components(ZZ, [1,2,3], 1)
            sage: c.__eq__(v)
            False
        """
        if isinstance(other, (int, Integer)): # other is 0
            if other == 0:
                return self.is_zero()
            else:
                raise TypeError("cannot compare a set of components to a number")
        else: # other is another Components
            if not isinstance(other, Components):
                raise TypeError("an instance of Components is expected")
            if other._frame != self._frame:
                return False
            if other._nid != self._nid:
                return False
            if other._sindex != self._sindex:
                return False
            if other._output_formatter != self._output_formatter:
                return False
            return (self - other).is_zero()

    def __ne__(self, other):
        r"""
        Non-equality operator.

        INPUT:

        - ``other`` -- set of components or 0

        OUTPUT: ``True`` if ``self`` is different from ``other``, or ``False``
        otherwise

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 1)
            sage: c.__ne__(0)  # uninitialized components are zero
            False
            sage: c1 = Components(ZZ, [1,2,3], 1)
            sage: c.__ne__(c1)  # c and c1 are both zero
            False
            sage: c[0] = 4
            sage: c.__ne__(0)
            True
            sage: c.__ne__(c1)
            True
        """
        return not self == other

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT: an exact copy of ``self``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 1)
            sage: c[:] = 5, 0, -4
            sage: a = c.__pos__() ; a
            1-index components w.r.t. [1, 2, 3]
            sage: a[:]
            [5, 0, -4]
            sage: a == +c
            True
            sage: a == c
            True
        """
        return self.copy()

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT: the opposite of the components represented by ``self``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: c = Components(ZZ, [1,2,3], 1)
            sage: c[:] = 5, 0, -4
            sage: a = c.__neg__() ; a
            1-index components w.r.t. [1, 2, 3]
            sage: a[:]
            [-5, 0, 4]
            sage: a == -c
            True
        """
        result = self._new_instance()
        for ind, val in self._comp.items():
            result._comp[ind] = - val
        return result

    def __add__(self, other):
        r"""
        Component addition.

        INPUT:

        - ``other`` -- components of the same number of indices and defined
          on the same frame as ``self``

        OUTPUT: components resulting from the addition of ``self`` and ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = Components(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__add__(b) ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [5, 5, 3]
            sage: s == a+b
            True

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: s_par = a.__add__(b) ; s_par
            1-index components w.r.t. [1, 2, 3]
            sage: s_par[:]
            [5, 5, 3]
            sage: s_par == s
            True
            sage: b.__add__(a) == s  # test of commutativity of parallel comput.
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization
        """
        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        if not isinstance(other, Components):
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
        # Initialization of the result to self.copy(), so that there remains
        # only to add other:
        result = self.copy()
        nproc = Parallelism().get('tensor')
        if nproc != 1:
            # Parallel computation
            lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
            ind_list = list(other._comp)
            ind_step = max(1, int(len(ind_list)/nproc/2))
            local_list = lol(ind_list, ind_step)
            # list of input parameters
            listParalInput = [(self, other, ind_part) for ind_part in local_list]

            @parallel(p_iter='multiprocessing', ncpus=nproc)
            def paral_sum(a, b, local_list_ind):
                partial = []
                for ind in local_list_ind:
                    partial.append([ind, a[[ind]]+b[[ind]]])
                return partial

            for ii, val in paral_sum(listParalInput):
                for jj in val:
                    result[[jj[0]]] = jj[1]

        else:
            # Sequential computation
            for ind, val in other._comp.items():
                result[[ind]] += val

        return result

    def __radd__(self, other):
        r"""
        Reflected addition (addition on the right to ``other``).

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = Components(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__radd__(b) ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [5, 5, 3]
            sage: s == a+b
            True
            sage: s = 0 + a ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s == a
            True
        """
        return self + other

    def __sub__(self, other):
        r"""
        Component subtraction.

        INPUT:

        - ``other`` -- components, of the same type as ``self``

        OUTPUT: components resulting from the subtraction of ``other`` from ``self``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = Components(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__sub__(b) ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [-3, -5, -9]
            sage: s == a - b
            True

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: Parallelism().get('tensor')
            2
            sage: s_par = a.__sub__(b) ; s_par
            1-index components w.r.t. [1, 2, 3]
            sage: s_par[:]
            [-3, -5, -9]
            sage: s_par == s
            True
            sage: b.__sub__(a) == -s # test of parallel comput.
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization
        """
        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        return self + (-other)  # ! # correct, deals properly with
                                # symmetries, but is probably not optimal

    def __rsub__(self, other):
        r"""
        Reflected subtraction (subtraction from ``other``).

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = Components(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__rsub__(b) ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [3, 5, 9]
            sage: s == b - a
            True
            sage: s = 0 - a ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [-1, 0, 3]
            sage: s == -a
            True
        """
        return (-self) + other

    def __mul__(self, other):
        r"""
        Component tensor product.

        INPUT:

        - ``other`` -- components, on the same frame as ``self``

        OUTPUT: the tensor product of ``self`` by ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: b = Components(ZZ, [1,2,3], 1)
            sage: b[:] = 4, 5, 6
            sage: s = a.__mul__(b) ; s
            2-indices components w.r.t. [1, 2, 3]
            sage: s[:]
            [  4   5   6]
            [  0   0   0]
            [-12 -15 -18]
            sage: s == a*b
            True
            sage: t = b*a
            sage: aa = a*a ; aa
            Fully symmetric 2-indices components w.r.t. [1, 2, 3]
            sage: aa[:]
            [ 1  0 -3]
            [ 0  0  0]
            [-3  0  9]

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: Parallelism().get('tensor')
            2
            sage: s_par = a.__mul__(b) ; s_par
            2-indices components w.r.t. [1, 2, 3]
            sage: s_par[:]
            [  4   5   6]
            [  0   0   0]
            [-12 -15 -18]
            sage: s_par == s
            True
            sage: b.__mul__(a) == t  # test of parallel comput.
            True
            sage: a.__mul__(a) == aa  # test of parallel comput.
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        Tensor product with a set of symmetric components::

            sage: s = b*aa; s
            3-indices components w.r.t. [1, 2, 3], with symmetry on the index positions (1, 2)
            sage: s[:]
            [[[4, 0, -12], [0, 0, 0], [-12, 0, 36]],
             [[5, 0, -15], [0, 0, 0], [-15, 0, 45]],
             [[6, 0, -18], [0, 0, 0], [-18, 0, 54]]]
            sage: Parallelism().set('tensor', nproc=2)
            sage: s_par = b*aa; s_par
            3-indices components w.r.t. [1, 2, 3], with symmetry on the index positions (1, 2)
            sage: s_par[:] == s[:]
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization
        """
        if not isinstance(other, Components):
            raise TypeError("the second argument for the tensor product " +
                            "must be an instance of Components")
        if other._frame != self._frame:
            raise ValueError("the two sets of components are not defined on " +
                             "the same frame")
        if other._sindex != self._sindex:
            raise ValueError("the two sets of components do not have the " +
                             "same starting index")
        if isinstance(other, CompWithSym):
            sym = []
            if other._sym:
                for s in other._sym:
                    ns = tuple(s[i]+self._nid for i in range(len(s)))
                    sym.append(ns)
            antisym = []
            if other._antisym:
                for s in other._antisym:
                    ns = tuple(s[i]+self._nid for i in range(len(s)))
                    antisym.append(ns)
            result = CompWithSym(self._ring, self._frame, self._nid + other._nid,
                                 self._sindex, self._output_formatter, sym,
                                 antisym)
        elif self._nid == 1 and other._nid == 1:
            if self is other:  # == would be dangerous here
                # The result is symmetric:
                result = CompFullySym(self._ring, self._frame, 2, self._sindex,
                                      self._output_formatter)
                # The loop below on self._comp.items() and
                # other._comp.items() cannot be used in the present case
                # (it would not deal correctly with redundant indices)
                # So we use a loop specific to the current case and return the
                # result:

                nproc = Parallelism().get('tensor')
                if nproc != 1:
                    # Parallel computation
                    lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
                    ind_list = list(result.non_redundant_index_generator())
                    ind_step = max(1, int(len(ind_list)/nproc))
                    local_list = lol(ind_list, ind_step)
                    # list of input parameters:
                    listParalInput = [(self, ind_part) for ind_part in local_list]

                    @parallel(p_iter='multiprocessing',ncpus=nproc)
                    def paral_mul(a, local_list_ind):
                        partial = []
                        for ind in local_list_ind:
                            partial.append([ind, a[[ind[0]]]*a[[ind[1]]]])
                        return partial

                    for ii,val in paral_mul(listParalInput):
                        for jj in val:
                            result[[jj[0]]] = jj[1]
                else:
                    # Sequential computation
                    for ind in result.non_redundant_index_generator():
                        result[[ind]] = self[[ind[0]]] * self[[ind[1]]]
                return result
            else:
                result = Components(self._ring, self._frame, 2, self._sindex,
                                    self._output_formatter)
        else:
            result = Components(self._ring, self._frame, self._nid + other._nid,
                                self._sindex, self._output_formatter)
        nproc = Parallelism().get('tensor')
        if nproc != 1:
            # Parallel computation
            lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
            ind_list = list(self._comp)
            ind_step = max(1, len(ind_list) // nproc)
            local_list = lol(ind_list, ind_step)
            # list of input parameters:
            listParalInput = [(self, other, ind_part) for ind_part in local_list]

            @parallel(p_iter='multiprocessing', ncpus=nproc)
            def paral_mul(a, b, local_list_ind):
                partial = []
                for ind in local_list_ind:
                    for ind_o, val_o in b._comp.items():
                        partial.append([ind + ind_o, a._comp[ind]*val_o])
                return partial

            for ii,val in paral_mul(listParalInput):
                for jj in val:
                    result._comp[jj[0]] = jj[1]
        else:
            # Sequential computation
            for ind_s, val_s in self._comp.items():
                for ind_o, val_o in other._comp.items():
                    result._comp[ind_s + ind_o] = val_s * val_o
        return result

    def __rmul__(self, other):
        r"""
        Reflected multiplication (multiplication on the left by ``other``).

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(ZZ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: s = a.__rmul__(2) ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [2, 0, -6]
            sage: s == 2*a
            True
            sage: a.__rmul__(0) == 0
            True
        """
        if isinstance(other, Components):
            raise NotImplementedError("left tensor product not implemented")
        # Left multiplication by a "scalar":
        result = self._new_instance()
        if other == 0:
            return result   # because a just created Components is zero
        for ind, val in self._comp.items():
            result._comp[ind] = other * val
        return result

    def __truediv__(self, other):
        r"""
        Division (by a scalar).

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: a = Components(QQ, [1,2,3], 1)
            sage: a[:] = 1, 0, -3
            sage: s = a.__truediv__(3) ; s
            1-index components w.r.t. [1, 2, 3]
            sage: s[:]
            [1/3, 0, -1]
            sage: s == a/3
            True
            sage: 3*s == a
            True
        """
        if isinstance(other, Components):
            raise NotImplementedError("division by an object of type " +
                                      "Components not implemented")
        result = self._new_instance()
        for ind, val in self._comp.items():
            result._comp[ind] = val / other
        return result

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

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ, 3)
            sage: c = Components(QQ, V.basis(), 2)
            sage: c[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: c.trace(0,1)
            15
            sage: c[0,0] + c[1,1] + c[2,2]  # check
            15

        Three self-contractions of a set of components with 3 indices::

            sage: v = Components(QQ, V.basis(), 1)
            sage: v[:] = (-1,2,3)
            sage: a = c*v ; a
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s = a.trace(0,1) ; s  # contraction on the first two indices
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s[:]
            [-15, 30, 45]
            sage: [sum(a[j,j,i] for j in range(3)) for i in range(3)]  # check
            [-15, 30, 45]
            sage: s = a.trace(0,2) ; s[:]  # contraction on the first and last indices
            [28, 32, 36]
            sage: [sum(a[j,i,j] for j in range(3)) for i in range(3)]  # check
            [28, 32, 36]
            sage: s = a.trace(1,2) ; s[:] # contraction on the last two indices
            [12, 24, 36]
            sage: [sum(a[i,j,j] for j in range(3)) for i in range(3)]  # check
            [12, 24, 36]
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
        si = self._sindex
        nsi = si + self._dim
        if self._nid == 2:
            res = 0
            for i in range(si, nsi):
                res += self[[i,i]]
            return res
        else:
            # More than 2 indices
            result = Components(self._ring, self._frame, self._nid - 2,
                                self._sindex, self._output_formatter)
            if pos1 > pos2:
                pos1, pos2 = (pos2, pos1)
            for ind, val in self._comp.items():
                if ind[pos1] == ind[pos2]:
                    # there is a contribution to the contraction
                    ind_res = ind[:pos1] + ind[pos1+1:pos2] + ind[pos2+1:]
                    result[[ind_res]] += val
            return result

    def contract(self, *args):
        r"""
        Contraction on one or many indices with another instance of
        :class:`Components`.

        INPUT:

        - ``pos1`` -- positions of the indices in ``self`` involved in the
          contraction; ``pos1`` must be a sequence of integers, with 0 standing
          for the first index position, 1 for the second one, etc. If ``pos1``
          is not provided, a single contraction on the last index position of
          ``self`` is assumed
        - ``other`` -- the set of components to contract with
        - ``pos2`` -- positions of the indices in ``other`` involved in the
          contraction, with the same conventions as for ``pos1``. If ``pos2``
          is not provided, a single contraction on the first index position of
          ``other`` is assumed

        OUTPUT: set of components resulting from the contraction

        EXAMPLES:

        Contraction of a 1-index set of components with a 2-index one::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ, 3)
            sage: a = Components(QQ, V.basis(), 1)
            sage: a[:] = (-1, 2, 3)
            sage: b = Components(QQ, V.basis(), 2)
            sage: b[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: s0 = a.contract(0, b, 0) ; s0
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s0[:]
            [28, 32, 36]
            sage: s0[:] == [sum(a[j]*b[j,i] for j in range(3)) for i in range(3)]  # check
            True
            sage: s1 = a.contract(0, b, 1) ; s1[:]
            [12, 24, 36]
            sage: s1[:] == [sum(a[j]*b[i,j] for j in range(3)) for i in range(3)]  # check
            True

        Parallel computations (see
        :class:`~sage.parallel.parallelism.Parallelism`)::

            sage: Parallelism().set('tensor', nproc=2)
            sage: Parallelism().get('tensor')
            2
            sage: s0_par = a.contract(0, b, 0) ; s0_par
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s0_par[:]
            [28, 32, 36]
            sage: s0_par == s0
            True
            sage: s1_par = a.contract(0, b, 1) ; s1_par[:]
            [12, 24, 36]
            sage: s1_par == s1
            True
            sage: Parallelism().set('tensor', nproc = 1)  # switch off parallelization

        Contraction on 2 indices::

            sage: c = a*b ; c
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s = c.contract(1,2, b, 0,1) ; s
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s[:]
            [-285, 570, 855]
            sage: [sum(sum(c[i,j,k]*b[j,k] for k in range(3)) # check
            ....:      for j in range(3)) for i in range(3)]
            [-285, 570, 855]

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: c_par = a*b ; c_par
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: c_par == c
            True
            sage: s_par = c_par.contract(1,2, b, 0,1) ; s_par
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s_par[:]
            [-285, 570, 855]
            sage: s_par == s
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        Consistency check with :meth:`trace`::

            sage: b = a*a ; b   # the tensor product of a with itself
            Fully symmetric 2-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: b[:]
            [ 1 -2 -3]
            [-2  4  6]
            [-3  6  9]
            sage: b.trace(0,1)
            14
            sage: a.contract(0, a, 0) == b.trace(0,1)
            True

        TESTS:

        Check that :issue:`32355` is fixed::

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: a = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: a[0,1] = 1
            sage: b = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: b[0,1], b[0,2] = 2, 3
            sage: a.contract(0, 1, b, 0, 1)
            4
            sage: a.contract(0, 1, b, 1, 0)
            -4
            sage: a.contract(1, 0, b, 0, 1)
            -4
            sage: a.contract(1, 0, b, 1, 0)
            4
            sage: Parallelism().set('tensor', nproc=2)  # same tests with parallelization
            sage: a.contract(0, 1, b, 0, 1)
            4
            sage: a.contract(0, 1, b, 1, 0)
            -4
            sage: a.contract(1, 0, b, 0, 1)
            -4
            sage: a.contract(1, 0, b, 1, 0)
            4
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization
        """
        #
        # Treatment of the input
        #
        nargs = len(args)
        for i, arg in enumerate(args):
            if isinstance(arg, Components):
                other = arg
                it = i
                break
        else:
            raise ValueError("a set of components must be provided in the " +
                             "argument list")
        if it == 0:
            pos1 = (self._nid - 1,)
        else:
            pos1 = args[:it]
        if it == nargs-1:
            pos2 = (0,)
        else:
            pos2 = args[it+1:]
        ncontr = len(pos1) # number of contractions
        if len(pos2) != ncontr:
            raise TypeError("different number of indices for the contraction")
        if other._frame != self._frame:
            raise TypeError("the two sets of components are not defined on " +
                            "the same frame")
        if other._sindex != self._sindex:
            raise TypeError("the two sets of components do not have the " +
                            "same starting index")
        contractions = [(pos1[i], pos2[i]) for i in range(ncontr)]
        res_nid = self._nid + other._nid - 2*ncontr
        #
        # Special case of a scalar result
        #
        if res_nid == 0:
            # Pairs of indices tuples for the contraction:
            ind_pairs = []
            for ind_s in self.index_generator():
                ind_o = [None] * ncontr
                for pos_s, pos_o in contractions:
                    ind_o[pos_o] = ind_s[pos_s]
                ind_pairs.append((ind_s, ind_o))

            if Parallelism().get('tensor') != 1:
                # parallel computation

                @parallel(p_iter='multiprocessing', ncpus=Parallelism().get('tensor'))
                def compprod(a, b):
                    return a*b

                # parallel list of inputs
                partial = list(compprod([(self[[ind_s]], other[[ind_o]])
                                         for ind_s, ind_o in ind_pairs]))
                res = sum(map(itemgetter(1), partial))
            else:
                # sequential computation
                res = 0
                for ind_s, ind_o in ind_pairs:
                    res += self[[ind_s]] * other[[ind_o]]

            return res

        #
        # Positions of self and other indices in the result
        #  (None = the position is involved in a contraction and therefore
        #   does not appear in the final result)
        #
        pos_s = [None for i in range(self._nid)]  # initialization
        pos_o = [None for i in range(other._nid)] # initialization
        shift = 0
        for pos in range(self._nid):
            for contract_pair in contractions:
                if pos == contract_pair[0]:
                    shift += 1
                    break
            else:
                pos_s[pos] = pos - shift
        for pos in range(other._nid):
            for contract_pair in contractions:
                if pos == contract_pair[1]:
                    shift += 1
                    break
            else:
                pos_o[pos] = self._nid + pos - shift
        rev_s = [pos_s.index(i) for i in range(self._nid-ncontr)]
        rev_o = [pos_o.index(i) for i in range(self._nid-ncontr, res_nid)]
        #
        # Determination of the symmetries of the result
        #
        max_len_sym = 0 # maximum length of symmetries in the result
        max_len_antisym = 0 # maximum length of antisymmetries in the result
        if res_nid > 1:  # no need to search for symmetries if res_nid == 1
            if isinstance(self, CompWithSym):
                s_sym = self._sym
                s_antisym = self._antisym
            else:
                s_sym = []
                s_antisym = []
            if isinstance(other, CompWithSym):
                o_sym = other._sym
                o_antisym = other._antisym
            else:
                o_sym = []
                o_antisym = []
            res_sym = []
            res_antisym = []
            for isym in s_sym:
                r_isym = []
                for pos in isym:
                    if pos_s[pos] is not None:
                        r_isym.append(pos_s[pos])
                if len(r_isym) > 1:
                    res_sym.append(r_isym)
                    max_len_sym = max(max_len_sym, len(r_isym))
            for isym in s_antisym:
                r_isym = []
                for pos in isym:
                    if pos_s[pos] is not None:
                        r_isym.append(pos_s[pos])
                if len(r_isym) > 1:
                    res_antisym.append(r_isym)
                    max_len_antisym = max(max_len_antisym, len(r_isym))
            for isym in o_sym:
                r_isym = []
                for pos in isym:
                    if pos_o[pos] is not None:
                        r_isym.append(pos_o[pos])
                if len(r_isym) > 1:
                    res_sym.append(r_isym)
                    max_len_sym = max(max_len_sym, len(r_isym))
            for isym in o_antisym:
                r_isym = []
                for pos in isym:
                    if pos_o[pos] is not None:
                        r_isym.append(pos_o[pos])
                if len(r_isym) > 1:
                    res_antisym.append(r_isym)
                    max_len_antisym = max(max_len_antisym, len(r_isym))
        #
        # Construction of the result object in view of the remaining symmetries:
        #
        if max_len_sym == 0 and max_len_antisym == 0:
            res = Components(self._ring, self._frame, res_nid,
                             start_index=self._sindex,
                             output_formatter=self._output_formatter)
        elif max_len_sym == res_nid:
            res = CompFullySym(self._ring, self._frame, res_nid,
                               start_index=self._sindex,
                               output_formatter=self._output_formatter)
        elif max_len_antisym == res_nid:
            res = CompFullyAntiSym(self._ring, self._frame, res_nid,
                                   start_index=self._sindex,
                                   output_formatter=self._output_formatter)
        else:
            res = CompWithSym(self._ring, self._frame, res_nid,
                              start_index=self._sindex,
                              output_formatter=self._output_formatter,
                              sym=res_sym, antisym=res_antisym)
        #
        # Performing the contraction
        #
        # To generate the indices tuples (of size ncontr) involved in the
        # the contraction, we create an empty instance of Components with
        # ncontr indices and call the method index_generator() on it:
        comp_for_contr = Components(self._ring, self._frame, ncontr,
                                    start_index=self._sindex)
        shift_o = self._nid - ncontr

        if Parallelism().get('tensor') != 1:
            # parallel computation
            nproc = Parallelism().get('tensor')
            lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
            ind_list = list(res.non_redundant_index_generator())
            ind_step = max(1, int(len(ind_list)/nproc/2))
            local_list = lol(ind_list, ind_step)

            listParalInput = []
            for ind_part in local_list:
                listParalInput.append((self, other, ind_part, rev_s, rev_o,
                                       shift_o, contractions, comp_for_contr))

            # definition of the parallel function
            @parallel(p_iter='multiprocessing', ncpus=nproc)
            def make_Contraction(this, other, local_list, rev_s, rev_o,
                                 shift_o,contractions, comp_for_contr):
                local_res = []
                for ind in local_list:
                    ind_s = [None for i in range(this._nid)]  # initialization
                    ind_o = [None for i in range(other._nid)] # initialization
                    for i, pos in enumerate(rev_s):
                        ind_s[pos] = ind[i]
                    for i, pos in enumerate(rev_o):
                        ind_o[pos] = ind[shift_o+i]
                    sm = 0
                    for ind_c in comp_for_contr.index_generator():
                        ic = 0
                        for pos_s, pos_o in contractions:
                            k = ind_c[ic]
                            ind_s[pos_s] = k
                            ind_o[pos_o] = k
                            ic += 1
                        sm += this[[ind_s]] * other[[ind_o]]
                    local_res.append([ind, sm])
                return local_res

            for ii, val in make_Contraction(listParalInput):
                for jj in val:
                    res[[jj[0]]] = jj[1]
        else:
            # sequential computation
            for ind in res.non_redundant_index_generator():
                ind_s = [None for i in range(self._nid)]  # initialization
                ind_o = [None for i in range(other._nid)] # initialization
                for i, pos in enumerate(rev_s):
                    ind_s[pos] = ind[i]
                for i, pos in enumerate(rev_o):
                    ind_o[pos] = ind[shift_o+i]
                sm = 0
                for ind_c in comp_for_contr.index_generator():
                    ic = 0
                    for pos_s, pos_o in contractions:
                        k = ind_c[ic]
                        ind_s[pos_s] = k
                        ind_o[pos_o] = k
                        ic += 1
                    sm += self[[ind_s]] * other[[ind_o]]
                res[[ind]] = sm

        return res

    def index_generator(self):
        r"""
        Generator of indices.

        OUTPUT: an iterable index

        EXAMPLES:

        Indices on a 3-dimensional vector space::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ,3)
            sage: c = Components(QQ, V.basis(), 1)
            sage: list(c.index_generator())
            [(0,), (1,), (2,)]
            sage: c = Components(QQ, V.basis(), 1, start_index=1)
            sage: list(c.index_generator())
            [(1,), (2,), (3,)]
            sage: c = Components(QQ, V.basis(), 2)
            sage: list(c.index_generator())
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0),
             (2, 1), (2, 2)]
        """
        si = self._sindex
        imax = self._dim - 1 + si
        ind = [si for k in range(self._nid)]
        while True:
            yield tuple(ind)
            for pos in range(self._nid-1,-1,-1):
                if ind[pos] != imax:
                    ind[pos] += 1
                    break
                elif pos != 0:
                    ind[pos] = si
                else:
                    return # end point reached

    def non_redundant_index_generator(self):
        r"""
        Generator of non redundant indices.

        In the absence of declared symmetries, all possible indices are
        generated. So this method is equivalent to :meth:`index_generator`.
        Only versions for derived classes with symmetries or antisymmetries
        are not trivial.

        OUTPUT: an iterable index

        EXAMPLES:

        Indices on a 3-dimensional vector space::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ,3)
            sage: c = Components(QQ, V.basis(), 2)
            sage: list(c.non_redundant_index_generator())
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0),
             (2, 1), (2, 2)]
            sage: c = Components(QQ, V.basis(), 2, start_index=1)
            sage: list(c.non_redundant_index_generator())
            [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1),
             (3, 2), (3, 3)]
        """
        yield from self.index_generator()

    def symmetrize(self, *pos):
        r"""
        Symmetrization over the given index positions.

        INPUT:

        - ``pos`` -- list of index positions involved in the
          symmetrization (with the convention position=0 for the first slot);
          if none, the symmetrization is performed over all the indices

        OUTPUT:

        - an instance of :class:`CompWithSym` describing the symmetrized
          components

        EXAMPLES:

        Symmetrization of 2-indices components::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ, 3)
            sage: c = Components(QQ, V.basis(), 2)
            sage: c[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: s = c.symmetrize() ; s
            Fully symmetric 2-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: c[:], s[:]
            (
            [1 2 3]  [1 3 5]
            [4 5 6]  [3 5 7]
            [7 8 9], [5 7 9]
            )
            sage: c.symmetrize() == c.symmetrize(0,1)
            True

        Full symmetrization of 3-indices components::

            sage: c = Components(QQ, V.basis(), 3)
            sage: c[:] = [[[1,2,3], [4,5,6], [7,8,9]], [[10,11,12], [13,14,15], [16,17,18]], [[19,20,21], [22,23,24], [25,26,27]]]
            sage: s = c.symmetrize() ; s
            Fully symmetric 3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: c[:], s[:]
            ([[[1, 2, 3], [4, 5, 6], [7, 8, 9]],
              [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
              [[19, 20, 21], [22, 23, 24], [25, 26, 27]]],
             [[[1, 16/3, 29/3], [16/3, 29/3, 14], [29/3, 14, 55/3]],
              [[16/3, 29/3, 14], [29/3, 14, 55/3], [14, 55/3, 68/3]],
              [[29/3, 14, 55/3], [14, 55/3, 68/3], [55/3, 68/3, 27]]])
            sage: all(s[i,j,k] == (c[i,j,k]+c[i,k,j]+c[j,k,i]+c[j,i,k]+c[k,i,j]+c[k,j,i])/6  # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: c.symmetrize() == c.symmetrize(0,1,2)
            True

        Partial symmetrization of 3-indices components::

            sage: s = c.symmetrize(0,1) ; s   # symmetrization on the first two indices
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1)
            sage: c[:], s[:]
            ([[[1, 2, 3], [4, 5, 6], [7, 8, 9]],
              [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
              [[19, 20, 21], [22, 23, 24], [25, 26, 27]]],
             [[[1, 2, 3], [7, 8, 9], [13, 14, 15]],
              [[7, 8, 9], [13, 14, 15], [19, 20, 21]],
              [[13, 14, 15], [19, 20, 21], [25, 26, 27]]])
            sage: all(s[i,j,k] == (c[i,j,k]+c[j,i,k])/2   # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s = c.symmetrize(1,2) ; s   # symmetrization on the last two indices
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (1, 2)
            sage: c[:], s[:]
            ([[[1, 2, 3], [4, 5, 6], [7, 8, 9]],
              [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
              [[19, 20, 21], [22, 23, 24], [25, 26, 27]]],
             [[[1, 3, 5], [3, 5, 7], [5, 7, 9]],
              [[10, 12, 14], [12, 14, 16], [14, 16, 18]],
              [[19, 21, 23], [21, 23, 25], [23, 25, 27]]])
            sage: all(s[i,j,k] == (c[i,j,k]+c[i,k,j])/2   # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s = c.symmetrize(0,2) ; s   # symmetrization on the first and last indices
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 2)
            sage: c[:], s[:]
            ([[[1, 2, 3], [4, 5, 6], [7, 8, 9]],
              [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
              [[19, 20, 21], [22, 23, 24], [25, 26, 27]]],
             [[[1, 6, 11], [4, 9, 14], [7, 12, 17]],
              [[6, 11, 16], [9, 14, 19], [12, 17, 22]],
              [[11, 16, 21], [14, 19, 24], [17, 22, 27]]])
            sage: all(s[i,j,k] == (c[i,j,k]+c[k,j,i])/2   # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        if not pos:
            pos = tuple(range(self._nid))
        else:
            if len(pos) < 2:
                raise ValueError("at least two index positions must be given")
            if len(pos) > self._nid:
                raise ValueError("number of index positions larger than the "
                                 "total number of indices")
        n_sym = len(pos) # number of indices involved in the symmetry
        if n_sym == self._nid:
            result = CompFullySym(self._ring, self._frame, self._nid, self._sindex,
                                  self._output_formatter)
        else:
            result = CompWithSym(self._ring, self._frame, self._nid, self._sindex,
                                 self._output_formatter, sym=pos)
        sym_group = SymmetricGroup(n_sym)
        for ind in result.non_redundant_index_generator():
            sum = 0
            for perm in sym_group.list():
                # action of the permutation on [0,1,...,n_sym-1]:
                perm_action = [x - 1 for x in perm.domain()]
                ind_perm = list(ind)
                for k in range(n_sym):
                    ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                sum += self[[ind_perm]]
            result[[ind]] = sum / sym_group.order()
        return result

    def antisymmetrize(self, *pos):
        r"""
        Antisymmetrization over the given index positions.

        INPUT:

        - ``pos`` -- list of index positions involved in the antisymmetrization
          (with the convention position=0 for the first slot); if none, the
          antisymmetrization is performed over all the indices

        OUTPUT: an instance of :class:`CompWithSym` describing the
        antisymmetrized components

        EXAMPLES:

        Antisymmetrization of 2-indices components::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ, 3)
            sage: c = Components(QQ, V.basis(), 2)
            sage: c[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: s = c.antisymmetrize() ; s
            Fully antisymmetric 2-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: c[:], s[:]
            (
            [1 2 3]  [ 0 -1 -2]
            [4 5 6]  [ 1  0 -1]
            [7 8 9], [ 2  1  0]
            )
            sage: c.antisymmetrize() == c.antisymmetrize(0,1)
            True

        Full antisymmetrization of 3-indices components::

            sage: c = Components(QQ, V.basis(), 3)
            sage: c[:] = [[[-1,-2,3], [4,-5,4], [-7,8,9]], [[10,10,12], [13,-14,15], [-16,17,19]], [[-19,20,21], [1,2,3], [-25,26,27]]]
            sage: s = c.antisymmetrize() ; s
            Fully antisymmetric 3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: c[:], s[:]
            ([[[-1, -2, 3], [4, -5, 4], [-7, 8, 9]],
              [[10, 10, 12], [13, -14, 15], [-16, 17, 19]],
              [[-19, 20, 21], [1, 2, 3], [-25, 26, 27]]],
             [[[0, 0, 0], [0, 0, -13/6], [0, 13/6, 0]],
              [[0, 0, 13/6], [0, 0, 0], [-13/6, 0, 0]],
              [[0, -13/6, 0], [13/6, 0, 0], [0, 0, 0]]])
            sage: all(s[i,j,k] == (c[i,j,k]-c[i,k,j]+c[j,k,i]-c[j,i,k]+c[k,i,j]-c[k,j,i])/6  # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: c.symmetrize() == c.symmetrize(0,1,2)
            True

        Partial antisymmetrization of 3-indices components::

            sage: s = c.antisymmetrize(0,1) ; s  # antisymmetrization on the first two indices
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 1)
            sage: c[:], s[:]
            ([[[-1, -2, 3], [4, -5, 4], [-7, 8, 9]],
              [[10, 10, 12], [13, -14, 15], [-16, 17, 19]],
              [[-19, 20, 21], [1, 2, 3], [-25, 26, 27]]],
             [[[0, 0, 0], [-3, -15/2, -4], [6, -6, -6]],
              [[3, 15/2, 4], [0, 0, 0], [-17/2, 15/2, 8]],
              [[-6, 6, 6], [17/2, -15/2, -8], [0, 0, 0]]])
            sage: all(s[i,j,k] == (c[i,j,k]-c[j,i,k])/2  # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s = c.antisymmetrize(1,2) ; s  # antisymmetrization on the last two indices
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (1, 2)
            sage: c[:], s[:]
            ([[[-1, -2, 3], [4, -5, 4], [-7, 8, 9]],
              [[10, 10, 12], [13, -14, 15], [-16, 17, 19]],
              [[-19, 20, 21], [1, 2, 3], [-25, 26, 27]]],
             [[[0, -3, 5], [3, 0, -2], [-5, 2, 0]],
              [[0, -3/2, 14], [3/2, 0, -1], [-14, 1, 0]],
              [[0, 19/2, 23], [-19/2, 0, -23/2], [-23, 23/2, 0]]])
            sage: all(s[i,j,k] == (c[i,j,k]-c[i,k,j])/2  # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True
            sage: s = c.antisymmetrize(0,2) ; s  # antisymmetrization on the first and last indices
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 2)
            sage: c[:], s[:]
            ([[[-1, -2, 3], [4, -5, 4], [-7, 8, 9]],
              [[10, 10, 12], [13, -14, 15], [-16, 17, 19]],
              [[-19, 20, 21], [1, 2, 3], [-25, 26, 27]]],
             [[[0, -6, 11], [0, -9, 3/2], [0, 12, 17]],
              [[6, 0, -4], [9, 0, 13/2], [-12, 0, -7/2]],
              [[-11, 4, 0], [-3/2, -13/2, 0], [-17, 7/2, 0]]])
            sage: all(s[i,j,k] == (c[i,j,k]-c[k,j,i])/2  # Check of the result:
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True

        The order of index positions in the argument does not matter::

            sage: c.antisymmetrize(1,0) == c.antisymmetrize(0,1)
            True
            sage: c.antisymmetrize(2,1) == c.antisymmetrize(1,2)
            True
            sage: c.antisymmetrize(2,0) == c.antisymmetrize(0,2)
            True
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        if not pos:
            pos = tuple(range(self._nid))
        else:
            if len(pos) < 2:
                raise ValueError("at least two index positions must be given")
            if len(pos) > self._nid:
                raise ValueError("number of index positions larger than the "
                                 "total number of indices")
        n_sym = len(pos) # number of indices involved in the antisymmetry
        if n_sym == self._nid:
            result = CompFullyAntiSym(self._ring, self._frame, self._nid,
                                      self._sindex, self._output_formatter)
        else:
            result = CompWithSym(self._ring, self._frame, self._nid, self._sindex,
                                 self._output_formatter, antisym=pos)
        sym_group = SymmetricGroup(n_sym)
        for ind in result.non_redundant_index_generator():
            sum = 0
            for perm in sym_group.list():
                # action of the permutation on [0,1,...,n_sym-1]:
                perm_action = [x - 1 for x in perm.domain()]
                ind_perm = list(ind)
                for k in range(n_sym):
                    ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                if perm.sign() == 1:
                    sum += self[[ind_perm]]
                else:
                    sum -= self[[ind_perm]]
            result[[ind]] = sum / sym_group.order()
        return result

    def _matrix_(self):
        r"""
        Convert a set of ring components with 2 indices into a matrix.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import Components
            sage: V = VectorSpace(QQ, 3)
            sage: c = Components(QQ, V.basis(), 2, start_index=1)
            sage: c[:] = [[-1,2,3], [4,-5,6], [7,8,-9]]
            sage: c._matrix_()
            [-1  2  3]
            [ 4 -5  6]
            [ 7  8 -9]

            sage: matrix(c) == c._matrix_()
            True
        """
        from sage.matrix.constructor import matrix
        if self._nid != 2:
            raise ValueError("the set of components must have 2 indices")
        si = self._sindex
        nsi = self._dim + si
        tab = [[self[[i,j]] for j in range(si, nsi)] for i in range(si, nsi)]
        return matrix(tab)


#******************************************************************************

class CompWithSym(Components):
    r"""
    Indexed set of ring elements forming some components with respect to a
    given "frame", with symmetries or antisymmetries regarding permutations
    of the indices.

    The "frame" can be a basis of some vector space or a vector frame on some
    manifold (i.e. a field of bases).
    The stored quantities can be tensor components or non-tensorial quantities,
    such as connection coefficients or structure coefficients.

    Subclasses of :class:`CompWithSym` are

    * :class:`CompFullySym` for fully symmetric components.
    * :class:`CompFullyAntiSym` for fully antisymmetric components.

    INPUT:

    - ``ring`` -- commutative ring in which each component takes its value
    - ``frame`` -- frame with respect to which the components are defined;
      whatever type ``frame`` is, it should have some method ``__len__()``
      implemented, so that ``len(frame)`` returns the dimension, i.e. the size
      of a single index range
    - ``nb_indices`` -- number of indices labeling the components
    - ``start_index`` -- (default: 0) first value of a single index;
      accordingly a component index i must obey
      ``start_index <= i <= start_index + dim - 1``, where ``dim = len(frame)``.
    - ``output_formatter`` -- (default: ``None``) function or unbound
      method called to format the output of the component access
      operator ``[...]`` (method __getitem__); ``output_formatter`` must take
      1 or 2 arguments: the 1st argument must be an instance of ``ring`` and
      the second one, if any, some format specification.
    - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries among
      the indices: each symmetry is described by a tuple containing the
      positions of the involved indices, with the convention ``position=0``
      for the first slot; for instance:

        * ``sym = (0, 1)`` for a symmetry between the 1st and 2nd indices
        * ``sym = [(0,2), (1,3,4)]`` for a symmetry between the 1st and 3rd
          indices and a symmetry between the 2nd, 4th and 5th indices.

    - ``antisym`` -- (default: ``None``) antisymmetry or list of antisymmetries
      among the indices, with the same convention as for ``sym``

    EXAMPLES:

    Symmetric components with 2 indices::

        sage: from sage.tensor.modules.comp import Components, CompWithSym
        sage: V = VectorSpace(QQ,3)
        sage: c = CompWithSym(QQ, V.basis(), 2, sym=(0,1))  # for demonstration only: it is preferable to use CompFullySym in this case
        sage: c[0,1] = 3
        sage: c[:]  # note that c[1,0] has been set automatically
        [0 3 0]
        [3 0 0]
        [0 0 0]

    Antisymmetric components with 2 indices::

        sage: c = CompWithSym(QQ, V.basis(), 2, antisym=(0,1))  # for demonstration only: it is preferable to use CompFullyAntiSym in this case
        sage: c[0,1] = 3
        sage: c[:]  # note that c[1,0] has been set automatically
        [ 0  3  0]
        [-3  0  0]
        [ 0  0  0]

    Internally, only non-redundant components are stored::

        sage: c._comp
        {(0, 1): 3}

    Components with 6 indices, symmetric among 3 indices (at position
    `(0, 1, 5)`) and antisymmetric among 2 indices (at position `(2, 4)`)::

        sage: c = CompWithSym(QQ, V.basis(), 6, sym=(0,1,5), antisym=(2,4))
        sage: c[0,1,2,0,1,2] = 3
        sage: c[1,0,2,0,1,2]  # symmetry between indices in position 0 and 1
        3
        sage: c[2,1,2,0,1,0]  # symmetry between indices in position 0 and 5
        3
        sage: c[0,2,2,0,1,1]  # symmetry between indices in position 1 and 5
        3
        sage: c[0,1,1,0,2,2]  # antisymmetry between indices in position 2 and 4
        -3

    Components with 4 indices, antisymmetric with respect to the first pair of
    indices as well as with the second pair of indices::

        sage: c = CompWithSym(QQ, V.basis(), 4, antisym=[(0,1),(2,3)])
        sage: c[0,1,0,1] = 3
        sage: c[1,0,0,1]  # antisymmetry on the first pair of indices
        -3
        sage: c[0,1,1,0]  # antisymmetry on the second pair of indices
        -3
        sage: c[1,0,1,0]  # consequence of the above
        3

    .. RUBRIC:: ARITHMETIC EXAMPLES

    Addition of a symmetric set of components with a non-symmetric one: the
    symmetry is lost::

        sage: V = VectorSpace(QQ, 3)
        sage: a = Components(QQ, V.basis(), 2)
        sage: a[:] = [[1,-2,3], [4,5,-6], [-7,8,9]]
        sage: b = CompWithSym(QQ, V.basis(), 2, sym=(0,1))  # for demonstration only: it is preferable to declare b = CompFullySym(QQ, V.basis(), 2)
        sage: b[0,0], b[0,1], b[0,2] = 1, 2, 3
        sage: b[1,1], b[1,2] = 5, 7
        sage: b[2,2] = 11
        sage: s = a + b ; s
        2-indices components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: a[:], b[:], s[:]
        (
        [ 1 -2  3]  [ 1  2  3]  [ 2  0  6]
        [ 4  5 -6]  [ 2  5  7]  [ 6 10  1]
        [-7  8  9], [ 3  7 11], [-4 15 20]
        )
        sage: a + b == b + a
        True

    Addition of two symmetric set of components: the symmetry is preserved::

        sage: c = CompWithSym(QQ, V.basis(), 2, sym=(0,1)) # for demonstration only: it is preferable to declare c = CompFullySym(QQ, V.basis(), 2)
        sage: c[0,0], c[0,1], c[0,2] = -4, 7, -8
        sage: c[1,1], c[1,2] = 2, -4
        sage: c[2,2] = 2
        sage: s = b + c ; s
        2-indices components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ], with symmetry on the index positions (0, 1)
        sage: b[:], c[:], s[:]
        (
        [ 1  2  3]  [-4  7 -8]  [-3  9 -5]
        [ 2  5  7]  [ 7  2 -4]  [ 9  7  3]
        [ 3  7 11], [-8 -4  2], [-5  3 13]
        )
        sage: b + c == c + b
        True

    Check of the addition with counterparts not declared symmetric::

        sage: bn = Components(QQ, V.basis(), 2)
        sage: bn[:] = b[:]
        sage: bn == b
        True
        sage: cn = Components(QQ, V.basis(), 2)
        sage: cn[:] = c[:]
        sage: cn == c
        True
        sage: bn + cn == b + c
        True

    Addition of an antisymmetric set of components with a non-symmetric one:
    the antisymmetry is lost::

        sage: d = CompWithSym(QQ, V.basis(), 2, antisym=(0,1))  # for demonstration only: it is preferable to declare d = CompFullyAntiSym(QQ, V.basis(), 2)
        sage: d[0,1], d[0,2], d[1,2] = 4, -1, 3
        sage: s = a + d ; s
        2-indices components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: a[:], d[:], s[:]
        (
        [ 1 -2  3]  [ 0  4 -1]  [ 1  2  2]
        [ 4  5 -6]  [-4  0  3]  [ 0  5 -3]
        [-7  8  9], [ 1 -3  0], [-6  5  9]
        )
        sage: d + a == a + d
        True

    Addition of two antisymmetric set of components: the antisymmetry is preserved::

        sage: e = CompWithSym(QQ, V.basis(), 2, antisym=(0,1))  # for demonstration only: it is preferable to declare e = CompFullyAntiSym(QQ, V.basis(), 2)
        sage: e[0,1], e[0,2], e[1,2] = 2, 3, -1
        sage: s = d + e ; s
        2-indices components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ], with antisymmetry on the index positions (0, 1)
        sage: d[:], e[:], s[:]
        (
        [ 0  4 -1]  [ 0  2  3]  [ 0  6  2]
        [-4  0  3]  [-2  0 -1]  [-6  0  2]
        [ 1 -3  0], [-3  1  0], [-2 -2  0]
        )
        sage: e + d == d + e
        True
    """
    def __init__(self, ring, frame, nb_indices, start_index=0,
                 output_formatter=None, sym=None, antisym=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: C = CompWithSym(ZZ, [1,2,3], 4, sym=(0,1), antisym=(2,3))
            sage: TestSuite(C).run()
        """
        Components.__init__(self, ring, frame, nb_indices, start_index,
                            output_formatter)
        self._sym, self._antisym = self._canonicalize_sym_antisym(
            nb_indices, sym, antisym)

    @staticmethod
    def _canonicalize_sym_or_antisym(nb_indices, sym_or_antisym):
        r"""
        Bring ``sym`` or ``antisym`` to its canonical form.

        INPUT:

        - ``nb_indices`` -- number of integer indices labeling the components

        - ``sym_or_antisym`` -- (default: ``None``) a symmetry/antisymmetry
          or an iterable of symmetries or an iterable of antisymmetries
          among the tensor arguments: each symmetry/antisymmetry is described
          by a tuple containing the positions of the involved arguments, with
          the convention ``position = 0`` for the first argument.

        TESTS::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: CompWithSym._canonicalize_sym_or_antisym(3, [0, -1])
            Traceback (most recent call last):
            ...
            IndexError: invalid index position: -1 not in [0,2]
            sage: CompWithSym._canonicalize_sym_or_antisym(3, [3, 1])
            Traceback (most recent call last):
            ...
            IndexError: invalid index position: 3 not in [0,2]
        """
        if not sym_or_antisym:
            return ()
        # Handle the case that sym_or_antisym is an iterator
        sym_or_antisym = tuple(sym_or_antisym)
        result_sym_or_antisym = []
        if isinstance(sym_or_antisym[0], (int, Integer)):
            # a single symmetry is provided as a tuple or a range object;
            # it is converted to a 1-item list:
            sym_or_antisym = (tuple(sym_or_antisym),)
        for isym in sym_or_antisym:
            if len(isym) < 2:
                # Drop trivial symmetry
                continue
            isym = tuple(sorted(isym))
            if isym[0] < 0:
                raise IndexError("invalid index position: " + str(isym[0]) +
                                 " not in [0," + str(nb_indices-1) + "]")
            if isym[-1] > nb_indices - 1:
                raise IndexError("invalid index position: " + str(isym[-1]) +
                                 " not in [0," + str(nb_indices-1) + "]")
            result_sym_or_antisym.append(isym)
        # Canonicalize sort order, make tuples
        return tuple(sorted(result_sym_or_antisym))

    @staticmethod
    def _canonicalize_sym_antisym(nb_indices, sym=None, antisym=None):
        r"""
        Bring ``sym`` and ``antisym`` into their canonical form.

        INPUT:

        - ``nb_indices`` -- number of integer indices labeling the components

        - ``sym`` -- (default: ``None``) a symmetry or an iterable of symmetries
          among the tensor arguments: each symmetry is described by a tuple
          containing the positions of the involved arguments, with the
          convention ``position = 0`` for the first argument. For instance:

          * ``sym = (0,1)`` for a symmetry between the 1st and 2nd arguments
          * ``sym = [(0,2), (1,3,4)]`` for a symmetry between the 1st and 3rd
            arguments and a symmetry between the 2nd, 4th and 5th arguments.

        - ``antisym`` -- (default: ``None``) antisymmetry or iterable of
          antisymmetries among the arguments, with the same convention
          as for ``sym``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: CompWithSym._canonicalize_sym_antisym(6, [(2, 1)])
            (((1, 2),), ())
        """
        if not sym and not antisym:
            # fast path
            return (), ()
        result_sym = CompWithSym._canonicalize_sym_or_antisym(nb_indices, sym)
        result_antisym = CompWithSym._canonicalize_sym_or_antisym(nb_indices, antisym)
        # Final consistency check:
        index_list = []
        for isym in result_sym:
            index_list.extend(isym)
        for isym in result_antisym:
            index_list.extend(isym)
        if len(index_list) != len(set(index_list)):
            # There is a repeated index position:
            raise IndexError("incompatible lists of symmetries: the same " +
                             "index position appears more than once")
        return result_sym, result_antisym

    def _repr_symmetry(self):
        r"""
        Return a prefix and a suffix string describing the symmetry of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: cp = CompWithSym(QQ, [1,2,3], 4, sym=(0,1))
            sage: cp._repr_symmetry()
            ('', ', with symmetry on the index positions (0, 1)')
            sage: cp = CompWithSym(QQ, [1,2,3], 4, sym=((0,1),), antisym=((2,3),))
            sage: cp._repr_symmetry()
            ('',
             ', with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)')
        """
        description = ""
        for isym in self._sym:
            description += ", with symmetry on the index positions " + \
                           str(tuple(isym))
        for isym in self._antisym:
            description += ", with antisymmetry on the index positions " + \
                           str(tuple(isym))
        return "", description

    def _new_instance(self):
        r"""
        Create a :class:`CompWithSym` instance w.r.t. the same frame,
        and with the same number of indices and the same symmetries.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: c = CompWithSym(ZZ, [1,2,3], 4, sym=(0,1))
            sage: a = c._new_instance() ; a
            4-indices components w.r.t. [1, 2, 3], with symmetry on the index positions (0, 1)
        """
        return CompWithSym(self._ring, self._frame, self._nid, self._sindex,
                          self._output_formatter, self._sym, self._antisym)

    def _ordered_indices(self, indices):
        r"""
        Given a set of indices, return a set of indices with the indices
        at the positions of symmetries or antisymmetries being ordered,
        as well as some antisymmetry indicator.

        INPUT:

        - ``indices`` -- list of indices (possibly a single integer if
          ``self`` is a 1-index object)

        OUTPUT:

        - a pair ``(s,ind)`` where ``ind`` is a tuple that differs from the
          original list of indices by a reordering at the positions of
          symmetries and antisymmetries and

          * ``s = 0`` if the value corresponding to ``indices`` vanishes by
            antisymmetry (repeated indices); `ind` is then set to ``None``
          * ``s = 1`` if the value corresponding to ``indices`` is the same as
            that corresponding to ``ind``
          * ``s = -1`` if the value corresponding to ``indices`` is the
            opposite of that corresponding to ``ind``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: c = CompWithSym(ZZ, [1,2,3], 4, sym=(0,1), antisym=(2,3))
            sage: c._ordered_indices([0,1,1,2])
            (1, (0, 1, 1, 2))
            sage: c._ordered_indices([1,0,1,2])
            (1, (0, 1, 1, 2))
            sage: c._ordered_indices([0,1,2,1])
            (-1, (0, 1, 1, 2))
            sage: c._ordered_indices([0,1,2,2])
            (0, None)
        """
        from sage.combinat.permutation import Permutation
        ind = list(self._check_indices(indices))
        for isym in self._sym:
            indsym = []
            for pos in isym:
                indsym.append(ind[pos])
            indsym_ordered = sorted(indsym)
            for k, pos in enumerate(isym):
                ind[pos] = indsym_ordered[k]
        sign = 1
        for isym in self._antisym:
            indsym = []
            for pos in isym:
                indsym.append(ind[pos])
            # Returns zero if some index appears twice:
            if len(indsym) != len(set(indsym)):
                return (0, None)
            # From here, all the indices in indsym are distinct and we need
            # to determine whether they form an even permutation of their
            # ordered series
            indsym_ordered = sorted(indsym)
            for k, pos in enumerate(isym):
                ind[pos] = indsym_ordered[k]
            if indsym_ordered != indsym:
                # Permutation linking indsym_ordered to indsym:
                #  (the +1 is required to fulfill the convention of Permutation)
                perm = [indsym.index(i) + 1 for i in indsym_ordered]
                #c#     Permutation(perm).signature()
                sign *= Permutation(perm).signature()
        ind = tuple(ind)
        return (sign, ind)

    def __getitem__(self, args):
        r"""
        Return the component corresponding to the given indices.

        INPUT:

        - ``args`` -- list of indices (possibly a single integer if
          ``self`` is a 1-index object) or the character ``:`` for the full list
          of components

        OUTPUT:

        - the component corresponding to ``args`` or, if ``args`` = ``:``,
          the full list of components, in the form ``T[i][j]...`` for the components
          `T_{ij...}` (for a 2-indices object, a matrix is returned).

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: c = CompWithSym(ZZ, [1,2,3], 4, sym=(0,1), antisym=(2,3))
            sage: c.__getitem__((0,1,1,2)) # uninitialized components are zero
            0
            sage: c[0,1,1,2] = 5
            sage: c.__getitem__((0,1,1,2))
            5
            sage: c.__getitem__((1,0,1,2))
            5
            sage: c.__getitem__((0,1,2,1))
            -5
            sage: c[0,1,2,1]
            -5
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
            return self._get_list(indices, no_format, format_type)
        else:
            sign, ind = self._ordered_indices(indices)
            if (sign == 0) or (ind not in self._comp): # the value is zero:
                if no_format:
                    return self._ring.zero()
                elif format_type is None:
                    return self._output_formatter(self._ring.zero())
                else:
                    return self._output_formatter(self._ring.zero(),
                                                 format_type)
            else: # nonzero value
                if no_format:
                    if sign == 1:
                        return self._comp[ind]
                    else: # sign = -1
                        return -self._comp[ind]
                elif format_type is None:
                    if sign == 1:
                        return self._output_formatter(self._comp[ind])
                    else: # sign = -1
                        return self._output_formatter(-self._comp[ind])
                else:
                    if sign == 1:
                        return self._output_formatter(
                                                 self._comp[ind], format_type)
                    else: # sign = -1
                        return self._output_formatter(
                                                -self._comp[ind], format_type)

    def __setitem__(self, args, value):
        r"""
        Set the component corresponding to the given indices.

        INPUT:

        - ``args`` -- list of indices (possibly a single integer if
          ``self`` is a 1-index object) ; if ``[:]`` is provided, all the
          components are set
        - ``value`` -- the value to be set or a list of values if
          ``args = [:]``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: c = CompWithSym(ZZ, [1,2,3], 2, sym=(0,1))
            sage: c.__setitem__((1,2), 5)
            sage: c[:]
            [0 0 0]
            [0 0 5]
            [0 5 0]
            sage: c = CompWithSym(ZZ, [1,2,3], 2, antisym=(0,1))
            sage: c.__setitem__((1,2), 5)
            sage: c[:]
            [ 0  0  0]
            [ 0  0  5]
            [ 0 -5  0]
            sage: c.__setitem__((2,2), 5)
            Traceback (most recent call last):
            ...
            ValueError: by antisymmetry, the component cannot have a nonzero value for the indices (2, 2)
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
            self._set_list(indices, format_type, value)
        else:
            sign, ind = self._ordered_indices(indices)
            # Check for a zero value
            #   The fast method is_trivial_zero() is employed preferably
            #   to the (possibly expensive) direct comparison to zero:
            if hasattr(value, 'is_trivial_zero'):
                zero_value = value.is_trivial_zero()
            else:
                zero_value = value == 0
            if sign == 0:
                if not zero_value:
                    raise ValueError("by antisymmetry, the component cannot " +
                                     "have a nonzero value for the indices " +
                                     str(indices))
                if ind in self._comp:
                    del self._comp[ind]  # zero values are not stored
            elif zero_value:
                if ind in self._comp:
                    del self._comp[ind]  # zero values are not stored
            else:
                if format_type is None:
                    if sign == 1:
                        self._comp[ind] = self._ring(value)
                    else:   # sign = -1
                        self._comp[ind] = -self._ring(value)
                else:
                    if sign == 1:
                        self._comp[ind] = self._ring({format_type: value})
                    else:   # sign = -1
                        self._comp[ind] = -self._ring({format_type: value})

    def swap_adjacent_indices(self, pos1, pos2, pos3):
        r"""
        Swap two adjacent sets of indices.

        This method is essentially required to reorder the covariant and
        contravariant indices in the computation of a tensor product.

        The symmetries are preserved and the corresponding indices are adjusted
        consequently.

        INPUT:

        - ``pos1`` -- position of the first index of set 1 (with the convention
          position=0 for the first slot)
        - ``pos2`` -- position of the first index of set 2 = 1 + position of
          the last index of set 1 (since the two sets are adjacent)
        - ``pos3`` -- 1 + position of the last index of set 2

        OUTPUT: components with index set 1 permuted with index set 2

        EXAMPLES:

        Swap of the index in position 0 with the pair of indices in position
        (1,2) in a set of components antisymmetric with respect to the indices
        in position (1,2)::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: V = VectorSpace(QQ, 3)
            sage: c = CompWithSym(QQ, V.basis(), 3, antisym=(1,2))
            sage: c[0,0,1], c[0,0,2], c[0,1,2] = (1,2,3)
            sage: c[1,0,1], c[1,0,2], c[1,1,2] = (4,5,6)
            sage: c[2,0,1], c[2,0,2], c[2,1,2] = (7,8,9)
            sage: c[:]
            [[[0, 1, 2], [-1, 0, 3], [-2, -3, 0]],
             [[0, 4, 5], [-4, 0, 6], [-5, -6, 0]],
             [[0, 7, 8], [-7, 0, 9], [-8, -9, 0]]]
            sage: c1 = c.swap_adjacent_indices(0,1,3)
            sage: c._antisym   # c is antisymmetric with respect to the last pair of indices...
            ((1, 2),)
            sage: c1._antisym  #...while c1 is antisymmetric with respect to the first pair of indices
            ((0, 1),)
            sage: c[0,1,2]
            3
            sage: c1[1,2,0]
            3
            sage: c1[2,1,0]
            -3
        """
        result = self._new_instance()
        # The symmetries:
        lpos = list(range(self._nid))
        new_lpos = lpos[:pos1] + lpos[pos2:pos3] + lpos[pos1:pos2] + lpos[pos3:]
        result._sym = []
        for s in self._sym:
            new_s = [new_lpos.index(pos) for pos in s]
            result._sym.append(tuple(sorted(new_s)))
        result._antisym = []
        for s in self._antisym:
            new_s = [new_lpos.index(pos) for pos in s]
            result._antisym.append(tuple(sorted(new_s)))
        result._sym, result._antisym = self._canonicalize_sym_antisym(
            self._nid, result._sym, result._antisym)
        # The values:
        for ind, val in self._comp.items():
            new_ind = ind[:pos1] + ind[pos2:pos3] + ind[pos1:pos2] + ind[pos3:]
            result[new_ind] = val
        return result

    def __add__(self, other):
        r"""
        Component addition.

        INPUT:

        - ``other`` -- components of the same number of indices and defined
          on the same frame as ``self``

        OUTPUT: components resulting from the addition of ``self`` and ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: a = CompWithSym(ZZ, [1,2,3], 2, sym=(0,1))
            sage: a[0,1], a[1,2] = 4, 5
            sage: b = CompWithSym(ZZ, [1,2,3], 2, sym=(0,1))
            sage: b[0,1], b[2,2] = 2, -3
            sage: s = a.__add__(b) ; s  # the symmetry is kept
            2-indices components w.r.t. [1, 2, 3], with symmetry on the index positions (0, 1)
            sage: s[:]
            [ 0  6  0]
            [ 6  0  5]
            [ 0  5 -3]
            sage: s == a + b
            True

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: s_par = a + b ; s_par
            2-indices components w.r.t. [1, 2, 3], with symmetry on the index positions (0, 1)
            sage: s_par[:] == s[:]
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        Addition with different symmetries::

            sage: c = CompWithSym(ZZ, [1,2,3], 2, antisym=(0,1))
            sage: c[0,1], c[0,2] = 3, 7
            sage: s = a.__add__(c) ; s  # the symmetry is lost
            2-indices components w.r.t. [1, 2, 3]
            sage: s[:]
            [ 0  7  7]
            [ 1  0  5]
            [-7  5  0]

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: s_par = a + c ; s_par
            2-indices components w.r.t. [1, 2, 3]
            sage: s_par[:] == s[:]
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization
        """
        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        if not isinstance(other, Components):
            raise TypeError("the second argument for the addition must be a " +
                            "an instance of Components")
        if other._frame != self._frame:
            raise ValueError("the two sets of components are not defined on " +
                             "the same frame")
        if other._nid != self._nid:
            raise ValueError("the two sets of components do not have the " +
                             "same number of indices")
        if other._sindex != self._sindex:
            raise ValueError("the two sets of components do not have the " +
                             "same starting index")
        if isinstance(other, CompWithSym):
            # Are the symmetries of the same type ?
            diff_sym = set(self._sym).symmetric_difference(set(other._sym))
            diff_antisym = \
                set(self._antisym).symmetric_difference(set(other._antisym))
            if diff_sym == set() and diff_antisym == set():
                # The symmetries/antisymmetries are identical:
                result = self.copy()
                nproc = Parallelism().get('tensor')
                if nproc != 1:
                    # Parallel computation
                    lol = lambda lst, sz: [lst[i:i+sz] for i in
                                           range(0, len(lst), sz)]
                    ind_list = list(other._comp)
                    ind_step = max(1, int(len(ind_list)/nproc/2))
                    local_list = lol(ind_list, ind_step)
                    # list of input parameters
                    listParalInput = [(self, other, ind_part) for ind_part in
                                      local_list]

                    @parallel(p_iter='multiprocessing', ncpus=nproc)
                    def paral_sum(a, b, local_list_ind):
                        partial = []
                        for ind in local_list_ind:
                            partial.append([ind, a[[ind]]+b[[ind]]])
                        return partial

                    for ii, val in paral_sum(listParalInput):
                        for jj in val:
                            result[[jj[0]]] = jj[1]

                else:
                    # Sequential computation
                    for ind, val in other._comp.items():
                        result[[ind]] += val
                return result
            else:
                # The symmetries/antisymmetries are different: only the
                # common ones are kept
                common_sym = []
                for isym in self._sym:
                    for osym in other._sym:
                        com = tuple(set(isym).intersection(set(osym)))
                        if len(com) > 1:
                            common_sym.append(com)
                common_antisym = []
                for isym in self._antisym:
                    for osym in other._antisym:
                        com = tuple(set(isym).intersection(set(osym)))
                        if len(com) > 1:
                            common_antisym.append(com)
                if common_sym or common_antisym:
                    result = CompWithSym(self._ring, self._frame, self._nid,
                                         self._sindex, self._output_formatter,
                                         common_sym, common_antisym)
                else:
                    # no common symmetry -> the result is a generic Components:
                    result = Components(self._ring, self._frame, self._nid,
                                        self._sindex, self._output_formatter)
        else:
            # other has no symmetry at all:
            result = Components(self._ring, self._frame, self._nid,
                                self._sindex, self._output_formatter)
        nproc = Parallelism().get('tensor')
        if nproc != 1:
            # Parallel computation
            lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
            ind_list = list(result.non_redundant_index_generator())
            ind_step = max(1, int(len(ind_list)/nproc/2))
            local_list = lol(ind_list, ind_step)
            # definition of the list of input parameters
            listParalInput = [(self, other, ind_part) for ind_part in local_list]

            @parallel(p_iter='multiprocessing', ncpus=nproc)
            def paral_sum(a, b, local_list_ind):
                partial = []
                for ind in local_list_ind:
                    partial.append([ind, a[[ind]]+b[[ind]]])
                return partial

            for ii,val in paral_sum(listParalInput):
                for jj in val:
                    result[[jj[0]]] = jj[1]
        else:
            # Sequential computation
            for ind in result.non_redundant_index_generator():
                result[[ind]] = self[[ind]] + other[[ind]]
        return result

    def __mul__(self, other):
        r"""
        Component tensor product.

        INPUT:

        - ``other`` -- components, on the same frame as ``self``

        OUTPUT: the tensor product of ``self`` by ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompWithSym
            sage: a = CompWithSym(ZZ, [1,2,3], 2, sym=(0,1))
            sage: a[0,1], a[1,2] = 4, 5
            sage: b = CompWithSym(ZZ, [1,2,3], 2, sym=(0,1))
            sage: b[0,1], b[2,2] = 2, -3
            sage: s1 = a.__mul__(b) ; s1
            4-indices components w.r.t. [1, 2, 3], with symmetry on the index positions (0, 1), with symmetry on the index positions (2, 3)
            sage: s1[1,0,0,1]
            8
            sage: s1[1,0,0,1] == a[1,0] * b[0,1]
            True
            sage: s1 == a*b
            True
            sage: c = CompWithSym(ZZ, [1,2,3], 2, antisym=(0,1))
            sage: c[0,1], c[0,2] = 3, 7
            sage: s2 = a.__mul__(c) ; s2
            4-indices components w.r.t. [1, 2, 3], with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)
            sage: s2[1,0,2,0]
            -28
            sage: s2[1,0,2,0] == a[1,0] * c[2,0]
            True
            sage: s2 == a*c
            True

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: Parallelism().get('tensor')
            2
            sage: s1_par = a.__mul__(b) ; s1_par
            4-indices components w.r.t. [1, 2, 3], with symmetry on the index positions (0, 1), with symmetry on the index positions (2, 3)
            sage: s1_par[1,0,0,1]
            8
            sage: s1_par[:] == s1[:]
            True
            sage: s2_par = a.__mul__(c) ; s2_par
            4-indices components w.r.t. [1, 2, 3], with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)
            sage: s2_par[1,0,2,0]
            -28
            sage: s2_par[:] == s2[:]
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization
        """
        if not isinstance(other, Components):
            raise TypeError("the second argument for the tensor product " +
                            "be an instance of Components")
        if other._frame != self._frame:
            raise ValueError("the two sets of components are not defined on " +
                             "the same frame")
        if other._sindex != self._sindex:
            raise ValueError("the two sets of components do not have the " +
                             "same starting index")
        sym = list(self._sym)
        antisym = list(self._antisym)
        if isinstance(other, CompWithSym):
            if other._sym:
                for s in other._sym:
                    ns = tuple(s[i]+self._nid for i in range(len(s)))
                    sym.append(ns)
            if other._antisym:
                for s in other._antisym:
                    ns = tuple(s[i]+self._nid for i in range(len(s)))
                    antisym.append(ns)
        result = CompWithSym(self._ring, self._frame, self._nid + other._nid,
                             self._sindex, self._output_formatter, sym, antisym)
        nproc = Parallelism().get('tensor')
        if nproc != 1:
            # Parallel computation
            lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
            ind_list = list(self._comp)
            ind_step = max(1, int(len(ind_list)/nproc))
            local_list = lol(ind_list, ind_step)
            # list of input parameters:
            listParalInput = [(self, other, ind_part) for ind_part in local_list]

            @parallel(p_iter='multiprocessing', ncpus=nproc)
            def paral_mul(a, b, local_list_ind):
                partial = []
                for ind in local_list_ind:
                    for ind_o, val_o in b._comp.items():
                        partial.append([ind + ind_o, a._comp[ind]*val_o])
                return partial

            for ii,val in paral_mul(listParalInput):
                for jj in val:
                    result._comp[jj[0]] = jj[1]
        else:
            # Sequential computation
            for ind_s, val_s in self._comp.items():
                for ind_o, val_o in other._comp.items():
                    result._comp[ind_s + ind_o] = val_s * val_o
        return result

    def trace(self, pos1, pos2):
        r"""
        Index contraction, taking care of the symmetries.

        INPUT:

        - ``pos1`` -- position of the first index for the contraction (with
          the convention position=0 for the first slot)
        - ``pos2`` -- position of the second index for the contraction

        OUTPUT:

        - set of components resulting from the (pos1, pos2) contraction

        EXAMPLES:

        Self-contraction of symmetric 2-indices components::

            sage: from sage.tensor.modules.comp import Components, CompWithSym, \
            ....:   CompFullySym, CompFullyAntiSym
            sage: V = VectorSpace(QQ, 3)
            sage: a = CompFullySym(QQ, V.basis(), 2)
            sage: a[:] = [[1,2,3],[2,4,5],[3,5,6]]
            sage: a.trace(0,1)
            11
            sage: a[0,0] + a[1,1] + a[2,2]
            11

        Self-contraction of antisymmetric 2-indices components::

            sage: b = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: b[0,1], b[0,2], b[1,2] = (3, -2, 1)
            sage: b.trace(0,1)  # must be zero by antisymmetry
            0

        Self-contraction of 3-indices components with one symmetry::

            sage: v = Components(QQ, V.basis(), 1)
            sage: v[:] = (-2, 4, -8)
            sage: c = v*b ; c
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (1, 2)
            sage: s = c.trace(0,1) ; s
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s[:]
            [-28, 2, 8]
            sage: [sum(v[k]*b[k,i] for k in range(3)) for i in range(3)] # check
            [-28, 2, 8]
            sage: s = c.trace(1,2) ; s
            1-index components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s[:] # is zero by antisymmetry
            [0, 0, 0]
            sage: c = b*v ; c
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 1)
            sage: s = c.trace(0,1)
            sage: s[:]  # is zero by antisymmetry
            [0, 0, 0]
            sage: s = c.trace(1,2) ; s[:]
            [28, -2, -8]
            sage: [sum(b[i,k]*v[k] for k in range(3)) for i in range(3)]  # check
            [28, -2, -8]

        Self-contraction of 4-indices components with two symmetries::

            sage: c = a*b ; c
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)
            sage: s = c.trace(0,1) ; s  # the symmetry on (0,1) is lost:
            Fully antisymmetric 2-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s[:]
            [  0  33 -22]
            [-33   0  11]
            [ 22 -11   0]
            sage: [[sum(c[k,k,i,j] for k in range(3)) for j in range(3)] for i in range(3)]  # check
            [[0, 33, -22], [-33, 0, 11], [22, -11, 0]]
            sage: s = c.trace(1,2) ; s  # both symmetries are lost by this contraction
            2-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s[:]
            [ 0  0  0]
            [-2  1  0]
            [-3  3 -1]
            sage: [[sum(c[i,k,k,j] for k in range(3)) for j in range(3)] for i in range(3)]  # check
            [[0, 0, 0], [-2, 1, 0], [-3, 3, -1]]
        """
        if self._nid < 2:
            raise TypeError("contraction can be performed only on " +
                            "components with at least 2 indices")
        if pos1 < 0 or pos1 > self._nid - 1:
            raise IndexError("pos1 out of range")
        if pos2 < 0 or pos2 > self._nid - 1:
            raise IndexError("pos2 out of range")
        if pos1 == pos2:
            raise IndexError("the two positions must differ for the " +
                             "contraction to take place")
        si = self._sindex
        nsi = si + self._dim
        if self._nid == 2:
            res = 0
            for i in range(si, nsi):
                res += self[[i,i]]
            return res
        else:
            # More than 2 indices
            if pos1 > pos2:
                pos1, pos2 = (pos2, pos1)
            # Determination of the remaining symmetries:
            sym_res = list(self._sym)
            for isym in self._sym:
                isym_res = list(isym)
                if pos1 in isym:
                    isym_res.remove(pos1)
                if pos2 in isym:
                    isym_res.remove(pos2)
                if len(isym_res) < 2:       # the symmetry is lost
                    sym_res.remove(isym)
                else:
                    sym_res[sym_res.index(isym)] = tuple(isym_res)
            antisym_res = list(self._antisym)
            for isym in self._antisym:
                isym_res = list(isym)
                if pos1 in isym:
                    isym_res.remove(pos1)
                if pos2 in isym:
                    isym_res.remove(pos2)
                if len(isym_res) < 2:       # the symmetry is lost
                    antisym_res.remove(isym)
                else:
                    antisym_res[antisym_res.index(isym)] = tuple(isym_res)
            # Shift of the index positions to take into account the
            # suppression of 2 indices:
            max_sym = 0
            for k in range(len(sym_res)):
                isym_res = []
                for pos in sym_res[k]:
                    if pos < pos1:
                        isym_res.append(pos)
                    elif pos < pos2:
                        isym_res.append(pos-1)
                    else:
                        isym_res.append(pos-2)
                max_sym = max(max_sym, len(isym_res))
                sym_res[k] = tuple(isym_res)
            max_antisym = 0
            for k in range(len(antisym_res)):
                isym_res = []
                for pos in antisym_res[k]:
                    if pos < pos1:
                        isym_res.append(pos)
                    elif pos < pos2:
                        isym_res.append(pos-1)
                    else:
                        isym_res.append(pos-2)
                max_antisym = max(max_antisym, len(isym_res))
                antisym_res[k] = tuple(isym_res)
            # Construction of the appropriate object in view of the
            # remaining symmetries:
            nid_res = self._nid - 2
            if max_sym == 0 and max_antisym == 0:
                result = Components(self._ring, self._frame, nid_res, self._sindex,
                                    self._output_formatter)
            elif max_sym == nid_res:
                result = CompFullySym(self._ring, self._frame, nid_res,
                                      self._sindex, self._output_formatter)
            elif max_antisym == nid_res:
                result = CompFullyAntiSym(self._ring, self._frame, nid_res,
                                          self._sindex, self._output_formatter)
            else:
                result = CompWithSym(self._ring, self._frame, nid_res,
                                     self._sindex, self._output_formatter,
                                     sym=sym_res, antisym=antisym_res)
            # The contraction itself:
            for ind_res in result.non_redundant_index_generator():
                ind = list(ind_res)
                ind.insert(pos1, 0)
                ind.insert(pos2, 0)
                res = 0
                for i in range(si, nsi):
                    ind[pos1] = i
                    ind[pos2] = i
                    res += self[[ind]]
                result[[ind_res]] = res
            return result

    def non_redundant_index_generator(self):
        r"""
        Generator of indices, with only ordered indices in case of symmetries,
        so that only non-redundant indices are generated.

        OUTPUT: an iterable index

        EXAMPLES:

        Indices on a 2-dimensional space::

            sage: from sage.tensor.modules.comp import Components, CompWithSym, \
            ....:  CompFullySym, CompFullyAntiSym
            sage: V = VectorSpace(QQ, 2)
            sage: c = CompFullySym(QQ, V.basis(), 2)
            sage: list(c.non_redundant_index_generator())
            [(0, 0), (0, 1), (1, 1)]
            sage: c = CompFullySym(QQ, V.basis(), 2, start_index=1)
            sage: list(c.non_redundant_index_generator())
            [(1, 1), (1, 2), (2, 2)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: list(c.non_redundant_index_generator())
            [(0, 1)]

        Indices on a 3-dimensional space::

            sage: V = VectorSpace(QQ, 3)
            sage: c = CompFullySym(QQ, V.basis(), 2)
            sage: list(c.non_redundant_index_generator())
            [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
            sage: c = CompFullySym(QQ, V.basis(), 2, start_index=1)
            sage: list(c.non_redundant_index_generator())
            [(1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: list(c.non_redundant_index_generator())
            [(0, 1), (0, 2), (1, 2)]
            sage: c = CompWithSym(QQ, V.basis(), 3, sym=(1,2))  # symmetry on the last two indices
            sage: list(c.non_redundant_index_generator())
            [(0, 0, 0), (1, 0, 0), (2, 0, 0), (0, 0, 1), (1, 0, 1),
             (2, 0, 1), (0, 0, 2), (1, 0, 2), (2, 0, 2), (0, 1, 1),
             (1, 1, 1), (2, 1, 1), (0, 1, 2), (1, 1, 2), (2, 1, 2),
             (0, 2, 2), (1, 2, 2), (2, 2, 2)]
            sage: c = CompWithSym(QQ, V.basis(), 3, antisym=(1,2))  # antisymmetry on the last two indices
            sage: list(c.non_redundant_index_generator())
            [(0, 0, 1), (1, 0, 1), (2, 0, 1), (0, 0, 2), (1, 0, 2),
             (2, 0, 2), (0, 1, 2), (1, 1, 2), (2, 1, 2)]
            sage: c = CompFullySym(QQ, V.basis(), 3)
            sage: list(c.non_redundant_index_generator())
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 1), (0, 1, 2), (0, 2, 2),
             (1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 3)
            sage: list(c.non_redundant_index_generator())
            [(0, 1, 2)]

        Indices on a 4-dimensional space::

            sage: V = VectorSpace(QQ, 4)
            sage: c = Components(QQ, V.basis(), 1)
            sage: list(c.non_redundant_index_generator())
            [(0,), (1,), (2,), (3,)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: list(c.non_redundant_index_generator())
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 3)
            sage: list(c.non_redundant_index_generator())
            [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 4)
            sage: list(c.non_redundant_index_generator())
            [(0, 1, 2, 3)]
            sage: c = CompFullyAntiSym(QQ, V.basis(), 5)
            sage: list(c.non_redundant_index_generator())  # nothing since c is identically zero in this case (for 5 > 4)
            []
        """
        if self._nid == 0 or self._dim == 0:
            return
        si = self._sindex
        imax = self._dim - 1 + si
        ind = [si for k in range(self._nid)]
        sym = list(self._sym)  # we may modify this in the following
        antisym = self._antisym
        for pos in range(self._nid):
            for isym in antisym:
                for k in range(1, len(isym)):
                    if pos == isym[k]:
                        if ind[isym[k-1]] == imax:
                            return
                        ind[pos] = ind[isym[k-1]] + 1
            if not any(pos in isym for isym in sym) and not any(pos in isym for isym in antisym):
                sym.append([pos]) # treat non-symmetrized indices as being symmetrized with themselves
        while True:
            yield tuple(ind)
            step_finished = False # each step generates a new index
            for i in range(len(sym)-1,-1,-1):
                # start with symmetrized indices, loop until we find
                # an index which we can increase without going over
                # the maximum
                isym = sym[i]
                if not step_finished:
                    for k in range(len(isym)-1,-1,-1):
                        if ind[isym[k]] != imax:
                            # we have found an index which we can
                            # increase; adjust other indices in the
                            # `isym` symmetrization
                            ind[isym[k]] += 1
                            for l in range(k+1, len(isym)):
                                ind[isym[l]] = ind[isym[l-1]]
                            step_finished = True
                            break
                        else:
                            # this index is at the maximum and we have
                            # to reset it
                            ind[isym[k]] = si
                if not step_finished and i == 0 and len(antisym) == 0:
                    return # we went through all indices and didn't
                           # find one which we can increase, thus we
                           # have generated all indices
            for i in range(len(antisym)-1,-1,-1):
                # the antisymmetrized indices work similar to the
                # symmetrized ones
                isym = antisym[i]
                if not step_finished:
                    for k in range(len(isym)-1,-1,-1):
                        if ind[isym[k]] + len(isym)-1-k != imax:
                            ind[isym[k]] += 1
                            for l in range(k+1, len(isym)):
                                # adjust antisymmetrized index
                                ind[isym[l]] = ind[isym[l-1]] + 1
                            step_finished = True
                            break
                        else:
                            # reset the index
                            ind[isym[k]] = si + k
                if not step_finished and i == 0:
                    return # end point reach

    def symmetrize(self, *pos):
        r"""
        Symmetrization over the given index positions.

        INPUT:

        - ``pos`` -- list of index positions involved in the
          symmetrization (with the convention ``position=0`` for the first
          slot); if none, the symmetrization is performed over all the indices

        OUTPUT:

        - an instance of :class:`CompWithSym` describing the symmetrized
          components

        EXAMPLES:

        Symmetrization of 3-indices components on a 3-dimensional space::

            sage: from sage.tensor.modules.comp import Components, CompWithSym, \
            ....:   CompFullySym, CompFullyAntiSym
            sage: V = VectorSpace(QQ, 3)
            sage: c = Components(QQ, V.basis(), 3)
            sage: c[:] = [[[1,2,3], [4,5,6], [7,8,9]], [[10,11,12], [13,14,15], [16,17,18]], [[19,20,21], [22,23,24], [25,26,27]]]
            sage: cs = c.symmetrize(0,1) ; cs
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1)
            sage: s = cs.symmetrize() ; s
            Fully symmetric 3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: cs[:], s[:]
            ([[[1, 2, 3], [7, 8, 9], [13, 14, 15]],
              [[7, 8, 9], [13, 14, 15], [19, 20, 21]],
              [[13, 14, 15], [19, 20, 21], [25, 26, 27]]],
             [[[1, 16/3, 29/3], [16/3, 29/3, 14], [29/3, 14, 55/3]],
              [[16/3, 29/3, 14], [29/3, 14, 55/3], [14, 55/3, 68/3]],
              [[29/3, 14, 55/3], [14, 55/3, 68/3], [55/3, 68/3, 27]]])
            sage: s == c.symmetrize() # should be true
            True
            sage: s1 = cs.symmetrize(0,1) ; s1   # should return a copy of cs
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1)
            sage: s1 == cs    # check that s1 is a copy of cs
            True

        Let us now start with a symmetry on the last two indices::

            sage: cs1 = c.symmetrize(1,2) ; cs1
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (1, 2)
            sage: s2 = cs1.symmetrize() ; s2
            Fully symmetric 3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s2 == c.symmetrize()
            True

        Symmetrization alters pre-existing symmetries: let us symmetrize w.r.t.
        the index positions `(1, 2)` a set of components that is symmetric
        w.r.t. the index positions `(0, 1)`::

            sage: cs = c.symmetrize(0,1) ; cs
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1)
            sage: css = cs.symmetrize(1,2)
            sage: css # the symmetry (0,1) has been lost:
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (1, 2)
            sage: css[:]
            [[[1, 9/2, 8], [9/2, 8, 23/2], [8, 23/2, 15]],
             [[7, 21/2, 14], [21/2, 14, 35/2], [14, 35/2, 21]],
             [[13, 33/2, 20], [33/2, 20, 47/2], [20, 47/2, 27]]]
            sage: cs[:]
            [[[1, 2, 3], [7, 8, 9], [13, 14, 15]],
             [[7, 8, 9], [13, 14, 15], [19, 20, 21]],
             [[13, 14, 15], [19, 20, 21], [25, 26, 27]]]
            sage: css == c.symmetrize() # css differs from the full symmetrized version
            False
            sage: css.symmetrize() == c.symmetrize() # one has to symmetrize css over all indices to recover it
            True

        Another example of symmetry alteration: symmetrization over `(0, 1)` of
        a 4-indices set of components that is symmetric w.r.t. `(1, 2, 3)`::

            sage: v = Components(QQ, V.basis(), 1)
            sage: v[:] = (-2,1,4)
            sage: a = v*s ; a
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (1, 2, 3)
            sage: a1 = a.symmetrize(0,1) ; a1 # the symmetry (1,2,3) has been reduced to (2,3):
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1), with symmetry on the index positions (2, 3)
            sage: a1._sym  # a1 has two distinct symmetries:
            ((0, 1), (2, 3))
            sage: a[0,1,2,0] == a[0,0,2,1]  # a is symmetric w.r.t. positions 1 and 3
            True
            sage: a1[0,1,2,0] == a1[0,0,2,1] # a1 is not
            False
            sage: a1[0,1,2,0] == a1[1,0,2,0] # but it is symmetric w.r.t. position 0 and 1
            True
            sage: a[0,1,2,0] == a[1,0,2,0] # while a is not
            False

        Partial symmetrization of 4-indices components with an antisymmetry on
        the last two indices::

            sage: a = Components(QQ, V.basis(), 2)
            sage: a[:] = [[-1,2,3], [4,5,-6], [7,8,9]]
            sage: b = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: b[0,1], b[0,2], b[1,2] = (2, 4, 8)
            sage: c = a*b ; c
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (2, 3)
            sage: s = c.symmetrize(0,1) ; s  # symmetrization on the first two indices
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1), with antisymmetry on the index positions (2, 3)
            sage: s[0,1,2,1] == (c[0,1,2,1] + c[1,0,2,1]) / 2 # check of the symmetrization
            True
            sage: s = c.symmetrize() ; s  # symmetrization over all the indices
            Fully symmetric 4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s == 0    # the full symmetrization results in zero due to the antisymmetry on the last two indices
            True
            sage: s = c.symmetrize(2,3) ; s
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (2, 3)
            sage: s == 0    # must be zero since the symmetrization has been performed on the antisymmetric indices
            True
            sage: s = c.symmetrize(0,2) ; s
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 2)
            sage: s != 0  # s is not zero, but the antisymmetry on (2,3) is lost because the position 2 is involved in the new symmetry
            True

        Partial symmetrization of 4-indices components with an antisymmetry on
        the last three indices::

            sage: a = Components(QQ, V.basis(), 1)
            sage: a[:] = (1, -2, 3)
            sage: b = CompFullyAntiSym(QQ, V.basis(), 3)
            sage: b[0,1,2] = 4
            sage: c = a*b ; c
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (1, 2, 3)
            sage: s = c.symmetrize(0,1) ; s
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1),
               with antisymmetry on the index positions (2, 3)

        Note that the antisymmetry on `(1, 2, 3)` has been reduced to
        `(2, 3)` only::

            sage: s = c.symmetrize(1,2) ; s
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (1, 2)
            sage: s == 0 # because (1,2) are involved in the original antisymmetry
            True
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        if not pos:
            pos = tuple(range(self._nid))
        else:
            if len(pos) < 2:
                raise ValueError("at least two index positions must be given")
            if len(pos) > self._nid:
                raise ValueError("number of index positions larger than the "
                                 "total number of indices")
            pos = tuple(pos)
        pos_set = set(pos)
        # If the symmetry is already present, there is nothing to do:
        for isym in self._sym:
            if pos_set.issubset(set(isym)):
                return self.copy()
        #
        # Interference of the new symmetry with existing ones:
        #
        sym_res = [pos]  # starting the list of symmetries of the result
        for isym in self._sym:
            inter = pos_set.intersection(set(isym))
            # if len(inter) == len(isym), isym is included in the new symmetry
            # and therefore has not to be included in sym_res
            if len(inter) != len(isym):
                if len(inter) >= 1:
                    # some part of isym is lost
                    isym_set = set(isym)
                    for k in inter:
                        isym_set.remove(k)
                    if len(isym_set) > 1:
                        # some part of isym remains and must be included in sym_res:
                        isym_res = tuple(isym_set)
                        sym_res.append(isym_res)
                else:
                    # case len(inter)=0: no interference: the existing symmetry is
                    # added to the list of symmetries for the result:
                    sym_res.append(isym)
        #
        # Interference of the new symmetry with existing antisymmetries:
        #
        antisym_res = []  # starting the list of antisymmetries of the result
        zero_result = False
        for iasym in self._antisym:
            inter = pos_set.intersection(set(iasym))
            if len(inter) > 1:
                # If at least two of the symmetry indices are already involved
                # in the antisymmetry, the outcome is zero:
                zero_result = True
            elif len(inter) == 1:
                # some piece of antisymmetry is lost
                k = inter.pop()  # the symmetry index position involved in the
                                 # antisymmetry
                iasym_set = set(iasym)
                iasym_set.remove(k)
                if len(iasym_set) > 1:
                    iasym_res = tuple(iasym_set)
                    antisym_res.append(iasym_res)
                # if len(iasym_set) == 1, the antisymmetry is fully lost, it is
                # therefore not appended to antisym_res
            else:
                # case len(inter)=0: no interference: the antisymmetry is
                # added to the list of antisymmetries for the result:
                antisym_res.append(iasym)
        #
        # Creation of the result object
        #
        max_sym = 0
        for isym in sym_res:
            max_sym = max(max_sym, len(isym))
        if max_sym == self._nid:
            result = CompFullySym(self._ring, self._frame, self._nid, self._sindex,
                                  self._output_formatter)
        else:
            result = CompWithSym(self._ring, self._frame, self._nid, self._sindex,
                                 self._output_formatter, sym=sym_res,
                                 antisym=antisym_res)
        if zero_result:
            return result   # since a just created instance is zero
        #
        # Symmetrization
        #
        n_sym = len(pos) # number of indices involved in the symmetry
        sym_group = SymmetricGroup(n_sym)
        for ind in result.non_redundant_index_generator():
            sum = 0
            for perm in sym_group.list():
                # action of the permutation on [0,1,...,n_sym-1]:
                perm_action = [x - 1 for x in perm.domain()]
                ind_perm = list(ind)
                for k in range(n_sym):
                    ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                sum += self[[ind_perm]]
            result[[ind]] = sum / sym_group.order()
        return result

    def antisymmetrize(self, *pos):
        r"""
        Antisymmetrization over the given index positions.

        INPUT:

        - ``pos`` -- list of index positions involved in the antisymmetrization
          (with the convention ``position=0`` for the first slot); if none, the
          antisymmetrization is performed over all the indices

        OUTPUT:

        - an instance of :class:`CompWithSym` describing the antisymmetrized
          components

        EXAMPLES:

        Antisymmetrization of 3-indices components on a 3-dimensional space::

            sage: from sage.tensor.modules.comp import Components, CompWithSym, \
            ....:  CompFullySym, CompFullyAntiSym
            sage: V = VectorSpace(QQ, 3)
            sage: a = Components(QQ, V.basis(), 1)
            sage: a[:] = (-2,1,3)
            sage: b = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: b[0,1], b[0,2], b[1,2] = (4,1,2)
            sage: c = a*b ; c   # tensor product of a by b
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (1, 2)
            sage: s = c.antisymmetrize() ; s
            Fully antisymmetric 3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: c[:], s[:]
            ([[[0, -8, -2], [8, 0, -4], [2, 4, 0]],
              [[0, 4, 1], [-4, 0, 2], [-1, -2, 0]],
              [[0, 12, 3], [-12, 0, 6], [-3, -6, 0]]],
             [[[0, 0, 0], [0, 0, 7/3], [0, -7/3, 0]],
              [[0, 0, -7/3], [0, 0, 0], [7/3, 0, 0]],
              [[0, 7/3, 0], [-7/3, 0, 0], [0, 0, 0]]])

        Check of the antisymmetrization::

            sage: all(s[i,j,k] == (c[i,j,k]-c[i,k,j]+c[j,k,i]-c[j,i,k]+c[k,i,j]-c[k,j,i])/6
            ....:     for i in range(3) for j in range(3) for k in range(3))
            True

        Antisymmetrization over already antisymmetric indices does not change anything::

            sage: s1 = s.antisymmetrize(1,2) ; s1
            Fully antisymmetric 3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s1 == s
            True
            sage: c1 = c.antisymmetrize(1,2) ; c1
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (1, 2)
            sage: c1 == c
            True

        But in general, antisymmetrization may alter previous antisymmetries::

            sage: c2 = c.antisymmetrize(0,1) ; c2  # the antisymmetry (2,3) is lost:
            3-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 1)
            sage: c2 == c
            False
            sage: c = s*a ; c
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 1, 2)
            sage: s = c.antisymmetrize(1,3) ; s
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 2),
               with antisymmetry on the index positions (1, 3)
            sage: s._antisym  # the antisymmetry (0,1,2) has been reduced to (0,2), since 1 is involved in the new antisymmetry (1,3):
            ((0, 2), (1, 3))

        Partial antisymmetrization of 4-indices components with a symmetry on
        the first two indices::

            sage: a = CompFullySym(QQ, V.basis(), 2)
            sage: a[:] = [[-2,1,3], [1,0,-5], [3,-5,4]]
            sage: b = Components(QQ, V.basis(), 2)
            sage: b[:] = [[1,2,3], [5,7,11], [13,17,19]]
            sage: c = a*b ; c
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1)
            sage: s = c.antisymmetrize(2,3) ; s
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with symmetry on the index positions (0, 1),
               with antisymmetry on the index positions (2, 3)

        Some check of the antisymmetrization::

            sage: all(s[2,2,i,j] == (c[2,2,i,j] - c[2,2,j,i])/2
            ....:     for i in range(3) for j in range(i,3))
            True

        The full antisymmetrization results in zero because of the symmetry on the
        first two indices::

            sage: s = c.antisymmetrize() ; s
            Fully antisymmetric 4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: s == 0
            True

        Similarly, the partial antisymmetrization on the first two indices results in zero::

            sage: s = c.antisymmetrize(0,1) ; s
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 1)
            sage: s == 0
            True

        The partial antisymmetrization on the positions `(0, 2)` destroys
        the symmetry on `(0, 1)`::

            sage: s = c.antisymmetrize(0,2) ; s
            4-indices components w.r.t. [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ], with antisymmetry on the index positions (0, 2)
            sage: s != 0
            True
            sage: s[0,1,2,1]
            27/2
            sage: s[1,0,2,1]  # the symmetry (0,1) is lost
            -2
            sage: s[2,1,0,1]  # the antisymmetry (0,2) holds
            -27/2
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        if not pos:
            pos = tuple(range(self._nid))
        else:
            if len(pos) < 2:
                raise ValueError("at least two index positions must be given")
            if len(pos) > self._nid:
                raise ValueError("number of index positions larger than the "
                                 "total number of indices")
            pos = tuple(pos)
        pos_set = set(pos)
        # If the antisymmetry is already present, there is nothing to do:
        for iasym in self._antisym:
            if pos_set.issubset(set(iasym)):
                return self.copy()
        #
        # Interference of the new antisymmetry with existing ones
        #
        antisym_res = [pos]  # starting the list of symmetries of the result
        for iasym in self._antisym:
            inter = pos_set.intersection(set(iasym))
            # if len(inter) == len(iasym), iasym is included in the new
            # antisymmetry and therefore has not to be included in antisym_res
            if len(inter) != len(iasym):
                if len(inter) >= 1:
                    # some part of iasym is lost
                    iasym_set = set(iasym)
                    for k in inter:
                        iasym_set.remove(k)
                    if len(iasym_set) > 1:
                        # some part of iasym remains and must be included in
                        # antisym_res:
                        iasym_res = tuple(iasym_set)
                        antisym_res.append(iasym_res)
                else:
                    # case len(inter)=0: no interference: the existing
                    # antisymmetry is added to the list of antisymmetries for
                    # the result:
                    antisym_res.append(iasym)
        #
        # Interference of the new antisymmetry with existing symmetries
        #
        sym_res = []  # starting the list of symmetries of the result
        zero_result = False
        for isym in self._sym:
            inter = pos_set.intersection(set(isym))
            if len(inter) > 1:
                # If at least two of the antisymmetry indices are already
                # involved in the symmetry, the outcome is zero:
                zero_result = True
            elif len(inter) == 1:
                # some piece of the symmetry is lost
                k = inter.pop()  # the antisymmetry index position involved in
                                 # the symmetry
                isym_set = set(isym)
                isym_set.remove(k)
                if len(isym_set) > 1:
                    isym_res = tuple(isym_set)
                    sym_res.append(isym_res)
                # if len(isym_set) == 1, the symmetry is fully lost, it is
                # therefore not appended to sym_res
            else:
                # case len(inter)=0: no interference: the symmetry is
                # added to the list of symmetries for the result:
                sym_res.append(isym)
        #
        # Creation of the result object
        #
        max_sym = 0
        for isym in antisym_res:
            max_sym = max(max_sym, len(isym))
        if max_sym == self._nid:
            result = CompFullyAntiSym(self._ring, self._frame, self._nid,
                                      self._sindex, self._output_formatter)
        else:
            result = CompWithSym(self._ring, self._frame, self._nid, self._sindex,
                                 self._output_formatter, sym=sym_res,
                                 antisym=antisym_res)
        if zero_result:
            return result   # since a just created instance is zero
        #
        # Antisymmetrization
        #
        n_sym = len(pos) # number of indices involved in the antisymmetry
        sym_group = SymmetricGroup(n_sym)
        for ind in result.non_redundant_index_generator():
            sum = 0
            for perm in sym_group.list():
                # action of the permutation on [0,1,...,n_sym-1]:
                perm_action = [x - 1 for x in perm.domain()]
                ind_perm = list(ind)
                for k in range(n_sym):
                    ind_perm[pos[perm_action[k]]] = ind[pos[k]]
                if perm.sign() == 1:
                    sum += self[[ind_perm]]
                else:
                    sum -= self[[ind_perm]]
            result[[ind]] = sum / sym_group.order()
        return result


#******************************************************************************

class CompFullySym(CompWithSym):
    r"""
    Indexed set of ring elements forming some components with respect to a
    given "frame" that are fully symmetric with respect to any permutation
    of the indices.

    The "frame" can be a basis of some vector space or a vector frame on some
    manifold (i.e. a field of bases).
    The stored quantities can be tensor components or non-tensorial quantities.

    INPUT:

    - ``ring`` -- commutative ring in which each component takes its value
    - ``frame`` -- frame with respect to which the components are defined;
      whatever type ``frame`` is, it should have some method ``__len__()``
      implemented, so that ``len(frame)`` returns the dimension, i.e. the size
      of a single index range
    - ``nb_indices`` -- number of indices labeling the components
    - ``start_index`` -- (default: 0) first value of a single index;
      accordingly a component index i must obey
      ``start_index <= i <= start_index + dim - 1``, where ``dim = len(frame)``.
    - ``output_formatter`` -- (default: ``None``) function or unbound
      method called to format the output of the component access
      operator ``[...]`` (method __getitem__); ``output_formatter`` must take
      1 or 2 arguments: the 1st argument must be an instance of ``ring`` and
      the second one, if any, some format specification.

    EXAMPLES:

    Symmetric components with 2 indices on a 3-dimensional space::

        sage: from sage.tensor.modules.comp import CompFullySym, CompWithSym
        sage: V = VectorSpace(QQ, 3)
        sage: c = CompFullySym(QQ, V.basis(), 2)
        sage: c[0,0], c[0,1], c[1,2] = 1, -2, 3
        sage: c[:] # note that c[1,0] and c[2,1] have been updated automatically (by symmetry)
        [ 1 -2  0]
        [-2  0  3]
        [ 0  3  0]

    Internally, only non-redundant and nonzero components are stored::

        sage: c._comp  # random output order of the component dictionary
        {(0, 0): 1, (0, 1): -2, (1, 2): 3}

    Same thing, but with the starting index set to 1::

        sage: c1 = CompFullySym(QQ, V.basis(), 2, start_index=1)
        sage: c1[1,1], c1[1,2], c1[2,3] = 1, -2, 3
        sage: c1[:]
        [ 1 -2  0]
        [-2  0  3]
        [ 0  3  0]

    The values stored in ``c`` and ``c1`` are equal::

        sage: c1[:] == c[:]
        True

    but not ``c`` and ``c1``, since their starting indices differ::

        sage: c1 == c
        False

    Fully symmetric components with 3 indices on a 3-dimensional space::

        sage: a = CompFullySym(QQ, V.basis(), 3)
        sage: a[0,1,2] = 3
        sage: a[:]
        [[[0, 0, 0], [0, 0, 3], [0, 3, 0]],
         [[0, 0, 3], [0, 0, 0], [3, 0, 0]],
         [[0, 3, 0], [3, 0, 0], [0, 0, 0]]]
        sage: a[0,1,0] = 4
        sage: a[:]
        [[[0, 4, 0], [4, 0, 3], [0, 3, 0]],
         [[4, 0, 3], [0, 0, 0], [3, 0, 0]],
         [[0, 3, 0], [3, 0, 0], [0, 0, 0]]]

    The full symmetry is preserved by the arithmetics::

        sage: b = CompFullySym(QQ, V.basis(), 3)
        sage: b[0,0,0], b[0,1,0], b[1,0,2], b[1,2,2] = -2, 3, 1, -5
        sage: s = a + 2*b ; s
        Fully symmetric 3-indices components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: a[:], b[:], s[:]
        ([[[0, 4, 0], [4, 0, 3], [0, 3, 0]],
          [[4, 0, 3], [0, 0, 0], [3, 0, 0]],
          [[0, 3, 0], [3, 0, 0], [0, 0, 0]]],
         [[[-2, 3, 0], [3, 0, 1], [0, 1, 0]],
          [[3, 0, 1], [0, 0, 0], [1, 0, -5]],
          [[0, 1, 0], [1, 0, -5], [0, -5, 0]]],
         [[[-4, 10, 0], [10, 0, 5], [0, 5, 0]],
          [[10, 0, 5], [0, 0, 0], [5, 0, -10]],
          [[0, 5, 0], [5, 0, -10], [0, -10, 0]]])

    It is lost if the added object is not fully symmetric::

        sage: b1 = CompWithSym(QQ, V.basis(), 3, sym=(0,1))  # b1 has only symmetry on index positions (0,1)
        sage: b1[0,0,0], b1[0,1,0], b1[1,0,2], b1[1,2,2] = -2, 3, 1, -5
        sage: s = a + 2*b1 ; s  # the result has the same symmetry as b1:
        3-indices components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ], with symmetry on the index positions (0, 1)
        sage: a[:], b1[:], s[:]
        ([[[0, 4, 0], [4, 0, 3], [0, 3, 0]],
          [[4, 0, 3], [0, 0, 0], [3, 0, 0]],
          [[0, 3, 0], [3, 0, 0], [0, 0, 0]]],
         [[[-2, 0, 0], [3, 0, 1], [0, 0, 0]],
          [[3, 0, 1], [0, 0, 0], [0, 0, -5]],
          [[0, 0, 0], [0, 0, -5], [0, 0, 0]]],
         [[[-4, 4, 0], [10, 0, 5], [0, 3, 0]],
          [[10, 0, 5], [0, 0, 0], [3, 0, -10]],
          [[0, 3, 0], [3, 0, -10], [0, 0, 0]]])
        sage: s = 2*b1 + a ; s
        3-indices components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ], with symmetry on the index positions (0, 1)
        sage: 2*b1 + a == a + 2*b1
        True
    """
    def __init__(self, ring, frame, nb_indices, start_index=0,
                 output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp import CompFullySym
            sage: C = CompFullySym(ZZ, (1,2,3), 2)
            sage: TestSuite(C).run()
        """
        CompWithSym.__init__(self, ring, frame, nb_indices, start_index,
                             output_formatter, sym=range(nb_indices))

    def _repr_symmetry(self):
        r"""
        Return a prefix and a suffix string describing the symmetry of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompFullySym
            sage: c = CompFullySym(ZZ, (1,2,3), 4)
            sage: c._repr_symmetry()
            ('Fully symmetric ', '')
        """
        return "Fully symmetric ", ""

    def _new_instance(self):
        r"""
        Create a :class:`CompFullySym` instance w.r.t. the same frame,
        and with the same number of indices.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompFullySym
            sage: c = CompFullySym(ZZ, (1,2,3), 4)
            sage: c._new_instance()
            Fully symmetric 4-indices components w.r.t. (1, 2, 3)
        """
        return CompFullySym(self._ring, self._frame, self._nid, self._sindex,
                            self._output_formatter)

    def __getitem__(self, args):
        r"""
        Return the component corresponding to the given indices of ``self``.

        INPUT:

        - ``args`` -- list of indices (possibly a single integer if
          ``self`` is a 1-index object) or the character ``:`` for the full list
          of components

        OUTPUT:

        - the component corresponding to ``args`` or, if ``args`` = ``:``,
          the full list of components, in the form ``T[i][j]...`` for the
          components `T_{ij...}` (for a 2-indices object, a matrix is returned)

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompFullySym
            sage: c = CompFullySym(ZZ, (1,2,3), 2)
            sage: c[0,1] = 4
            sage: c.__getitem__((0,1))
            4
            sage: c.__getitem__((1,0))
            4
            sage: c.__getitem__(slice(None))
            [0 4 0]
            [4 0 0]
            [0 0 0]
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
            return self._get_list(indices, no_format, format_type)

        ind = self._ordered_indices(indices)[1]  # [0]=sign is not used
        if ind in self._comp: # nonzero value
            if no_format:
                return self._comp[ind]
            elif format_type is None:
                return self._output_formatter(self._comp[ind])
            else:
                return self._output_formatter(self._comp[ind], format_type)

        # the value is zero
        if no_format:
            return self._ring.zero()
        elif format_type is None:
            return self._output_formatter(self._ring.zero())
        else:
            return self._output_formatter(self._ring.zero(),
                                             format_type)

    def __setitem__(self, args, value):
        r"""
        Set the component corresponding to the given indices.

        INPUT:

        - ``indices`` -- list of indices (possibly a single integer if
          ``self`` is a 1-index object); if [:] is provided, all the components
          are set
        - ``value`` -- the value to be set or a list of values if
          ``args = [:]``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompFullySym
            sage: c = CompFullySym(ZZ, (1,2,3), 2)
            sage: c.__setitem__((0,1), 4)
            sage: c[:]
            [0 4 0]
            [4 0 0]
            [0 0 0]
            sage: c.__setitem__((2,1), 5)
            sage: c[:]
            [0 4 0]
            [4 0 5]
            [0 5 0]
            sage: c.__setitem__(slice(None), [[1, 2, 3], [2, 4, 5], [3, 5, 6]])
            sage: c[:]
            [1 2 3]
            [2 4 5]
            [3 5 6]
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
            self._set_list(indices, format_type, value)
        else:
            ind = self._ordered_indices(indices)[1]  # [0]=sign is not used
            # Check for a zero value
            #   The fast method is_trivial_zero() is employed preferably
            #   to the (possibly expensive) direct comparison to zero:
            if hasattr(value, 'is_trivial_zero'):
                zero_value = value.is_trivial_zero()
            else:
                zero_value = value == 0
            if zero_value:
                # if the component has been set previously, it is deleted,
                # otherwise nothing is done (zero components are not stored):
                if ind in self._comp:
                    del self._comp[ind]  # zero values are not stored
            else:
                if format_type is None:
                    self._comp[ind] = self._ring(value)
                else:
                    self._comp[ind] = self._ring({format_type: value})

    def __add__(self, other):
        r"""
        Component addition.

        INPUT:

        - ``other`` -- components of the same number of indices and defined
          on the same frame as ``self``

        OUTPUT: components resulting from the addition of ``self`` and ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompFullySym
            sage: a = CompFullySym(ZZ, (1,2,3), 2)
            sage: a[0,1], a[1,2] = 4, 5
            sage: b = CompFullySym(ZZ, (1,2,3), 2)
            sage: b[0,1], b[2,2] = 2, -3
            sage: s = a.__add__(b) ; s  # the symmetry is kept
            Fully symmetric 2-indices components w.r.t. (1, 2, 3)
            sage: s[:]
            [ 0  6  0]
            [ 6  0  5]
            [ 0  5 -3]
            sage: s == a + b
            True

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: Parallelism().get('tensor')
            2
            sage: s_par = a.__add__(b) ; s_par
            Fully symmetric 2-indices components w.r.t. (1, 2, 3)
            sage: s_par[:]
            [ 0  6  0]
            [ 6  0  5]
            [ 0  5 -3]
            sage: s_par == s
            True
            sage: s_par == b.__add__(a)  # test of commutativity of parallel comput.
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization

        Addition with a set of components having different symmetries::

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: c = CompFullyAntiSym(ZZ, (1,2,3), 2)
            sage: c[0,1], c[0,2] = 3, 7
            sage: s = a.__add__(c) ; s  # the symmetry is lost
            2-indices components w.r.t. (1, 2, 3)
            sage: s[:]
            [ 0  7  7]
            [ 1  0  5]
            [-7  5  0]
            sage: s == a + c
            True

        Parallel computation::

            sage: Parallelism().set('tensor', nproc=2)
            sage: Parallelism().get('tensor')
            2
            sage: s_par = a.__add__(c) ; s_par
            2-indices components w.r.t. (1, 2, 3)
            sage: s_par[:]
            [ 0  7  7]
            [ 1  0  5]
            [-7  5  0]
            sage: s_par[:] == s[:]
            True
            sage: s_par == c.__add__(a)  # test of commutativity of parallel comput.
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization
        """
        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        if not isinstance(other, Components):
            raise TypeError("the second argument for the addition must be a " +
                            "an instance of Components")
        if isinstance(other, CompFullySym):
            if other._frame != self._frame:
                raise ValueError("the two sets of components are not defined " +
                                 "on the same frame")
            if other._nid != self._nid:
                raise ValueError("the two sets of components do not have the " +
                                 "same number of indices")
            if other._sindex != self._sindex:
                raise ValueError("the two sets of components do not have the " +
                                 "same starting index")
            # Initialization of the result to self.copy(), so that there
            # remains only to add other:
            result = self.copy()
            nproc = Parallelism().get('tensor')
            if nproc != 1:
                # parallel sum
                lol = lambda lst, sz: [lst[i:i+sz] for i
                                       in range(0, len(lst), sz)]
                ind_list = list(other._comp)
                ind_step = max(1, int(len(ind_list)/nproc/2))
                local_list = lol(ind_list, ind_step)
                # definition of the list of input parameters
                listParalInput = [(self, other, ind_part) for ind_part
                                  in local_list]

                @parallel(p_iter='multiprocessing', ncpus=nproc)
                def paral_sum(a, b, local_list_ind):
                    partial = []
                    for ind in local_list_ind:
                        partial.append([ind, a[[ind]]+b[[ind]]])
                    return partial

                for ii,val in paral_sum(listParalInput):
                    for jj in val:
                        result[[jj[0]]] = jj[1]

            else:
                # sequential sum
                for ind, val in other._comp.items():
                    result[[ind]] += val
            return result
        else:
            return CompWithSym.__add__(self, other)


#******************************************************************************

class CompFullyAntiSym(CompWithSym):
    r"""
    Indexed set of ring elements forming some components with respect to a
    given "frame" that are fully antisymmetric with respect to any permutation
    of the indices.

    The "frame" can be a basis of some vector space or a vector frame on some
    manifold (i.e. a field of bases).
    The stored quantities can be tensor components or non-tensorial quantities.

    INPUT:

    - ``ring`` -- commutative ring in which each component takes its value
    - ``frame`` -- frame with respect to which the components are defined;
      whatever type ``frame`` is, it should have some method ``__len__()``
      implemented, so that ``len(frame)`` returns the dimension, i.e. the size
      of a single index range
    - ``nb_indices`` -- number of indices labeling the components
    - ``start_index`` -- (default: 0) first value of a single index;
      accordingly a component index i must obey
      ``start_index <= i <= start_index + dim - 1``, where ``dim = len(frame)``.
    - ``output_formatter`` -- (default: ``None``) function or unbound
      method called to format the output of the component access
      operator ``[...]`` (method __getitem__); ``output_formatter`` must take
      1 or 2 arguments: the 1st argument must be an instance of ``ring`` and
      the second one, if any, some format specification.

    EXAMPLES:

    Antisymmetric components with 2 indices on a 3-dimensional space::

        sage: from sage.tensor.modules.comp import CompWithSym, CompFullyAntiSym
        sage: V = VectorSpace(QQ, 3)
        sage: c = CompFullyAntiSym(QQ, V.basis(), 2)
        sage: c[0,1], c[0,2], c[1,2] = 3, 1/2, -1
        sage: c[:]  # note that all components have been set according to the antisymmetry
        [   0    3  1/2]
        [  -3    0   -1]
        [-1/2    1    0]

    Internally, only non-redundant and nonzero components are stored::

        sage: c._comp # random output order of the component dictionary
        {(0, 1): 3, (0, 2): 1/2, (1, 2): -1}

    Same thing, but with the starting index set to 1::

        sage: c1 = CompFullyAntiSym(QQ, V.basis(), 2, start_index=1)
        sage: c1[1,2], c1[1,3], c1[2,3] = 3, 1/2, -1
        sage: c1[:]
        [   0    3  1/2]
        [  -3    0   -1]
        [-1/2    1    0]

    The values stored in ``c`` and ``c1`` are equal::

        sage: c1[:] == c[:]
        True

    but not ``c`` and ``c1``, since their starting indices differ::

        sage: c1 == c
        False

    Fully antisymmetric components with 3 indices on a 3-dimensional space::

        sage: a = CompFullyAntiSym(QQ, V.basis(), 3)
        sage: a[0,1,2] = 3  # the only independent component in dimension 3
        sage: a[:]
        [[[0, 0, 0], [0, 0, 3], [0, -3, 0]],
         [[0, 0, -3], [0, 0, 0], [3, 0, 0]],
         [[0, 3, 0], [-3, 0, 0], [0, 0, 0]]]

    Setting a nonzero value incompatible with the antisymmetry results in an
    error::

        sage: a[0,1,0] = 4
        Traceback (most recent call last):
        ...
        ValueError: by antisymmetry, the component cannot have a nonzero value for the indices (0, 1, 0)
        sage: a[0,1,0] = 0   # OK
        sage: a[2,0,1] = 3   # OK

    The full antisymmetry is preserved by the arithmetics::

        sage: b = CompFullyAntiSym(QQ, V.basis(), 3)
        sage: b[0,1,2] = -4
        sage: s = a + 2*b ; s
        Fully antisymmetric 3-indices components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]
        sage: a[:], b[:], s[:]
        ([[[0, 0, 0], [0, 0, 3], [0, -3, 0]],
          [[0, 0, -3], [0, 0, 0], [3, 0, 0]],
          [[0, 3, 0], [-3, 0, 0], [0, 0, 0]]],
         [[[0, 0, 0], [0, 0, -4], [0, 4, 0]],
          [[0, 0, 4], [0, 0, 0], [-4, 0, 0]],
          [[0, -4, 0], [4, 0, 0], [0, 0, 0]]],
         [[[0, 0, 0], [0, 0, -5], [0, 5, 0]],
          [[0, 0, 5], [0, 0, 0], [-5, 0, 0]],
          [[0, -5, 0], [5, 0, 0], [0, 0, 0]]])

    It is lost if the added object is not fully antisymmetric::

        sage: b1 = CompWithSym(QQ, V.basis(), 3, antisym=(0,1))  # b1 has only antisymmetry on index positions (0,1)
        sage: b1[0,1,2] = -4
        sage: s = a + 2*b1 ; s  # the result has the same symmetry as b1:
        3-indices components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ], with antisymmetry on the index positions (0, 1)
        sage: a[:], b1[:], s[:]
        ([[[0, 0, 0], [0, 0, 3], [0, -3, 0]],
          [[0, 0, -3], [0, 0, 0], [3, 0, 0]],
          [[0, 3, 0], [-3, 0, 0], [0, 0, 0]]],
         [[[0, 0, 0], [0, 0, -4], [0, 0, 0]],
          [[0, 0, 4], [0, 0, 0], [0, 0, 0]],
          [[0, 0, 0], [0, 0, 0], [0, 0, 0]]],
         [[[0, 0, 0], [0, 0, -5], [0, -3, 0]],
          [[0, 0, 5], [0, 0, 0], [3, 0, 0]],
          [[0, 3, 0], [-3, 0, 0], [0, 0, 0]]])
        sage: s = 2*b1 + a ; s
        3-indices components w.r.t. [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ], with antisymmetry on the index positions (0, 1)
        sage: 2*b1 + a == a + 2*b1
        True
    """
    def __init__(self, ring, frame, nb_indices, start_index=0,
                 output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: C = CompFullyAntiSym(ZZ, (1,2,3), 2)
            sage: TestSuite(C).run()
        """
        CompWithSym.__init__(self, ring, frame, nb_indices, start_index,
                             output_formatter, antisym=range(nb_indices))

    def _repr_symmetry(self):
        r"""
        Return a prefix and a suffix string describing the symmetry of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: c = CompFullyAntiSym(ZZ, (1,2,3), 4)
            sage: c._repr_symmetry()
            ('Fully antisymmetric ', '')
        """
        return "Fully antisymmetric ", ""

    def _new_instance(self):
        r"""
        Create a :class:`CompFullyAntiSym` instance w.r.t. the same frame,
        and with the same number of indices.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: c = CompFullyAntiSym(ZZ, (1,2,3), 4)
            sage: c._new_instance()
            Fully antisymmetric 4-indices components w.r.t. (1, 2, 3)
        """
        return CompFullyAntiSym(self._ring, self._frame, self._nid, self._sindex,
                                self._output_formatter)

    def __add__(self, other):
        r"""
        Component addition.

        INPUT:

        - ``other`` -- components of the same number of indices and defined
          on the same frame as ``self``

        OUTPUT: components resulting from the addition of ``self`` and ``other``

        EXAMPLES::

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: a = CompFullyAntiSym(ZZ, (1,2,3), 2)
            sage: a[0,1], a[1,2] = 4, 5
            sage: b = CompFullyAntiSym(ZZ, (1,2,3), 2)
            sage: b[0,1], b[0,2] = 2, -3
            sage: s = a.__add__(b) ; s  # the antisymmetry is kept
            Fully antisymmetric 2-indices components w.r.t. (1, 2, 3)
            sage: s[:]
            [ 0  6 -3]
            [-6  0  5]
            [ 3 -5  0]
            sage: s == a + b
            True
            sage: from sage.tensor.modules.comp import CompFullySym
            sage: c = CompFullySym(ZZ, (1,2,3), 2)
            sage: c[0,1], c[0,2] = 3, 7
            sage: s = a.__add__(c) ; s  # the antisymmetry is lost
            2-indices components w.r.t. (1, 2, 3)
            sage: s[:]
            [ 0  7  7]
            [-1  0  5]
            [ 7 -5  0]
            sage: s == a + c
            True

        Parallel computation::

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: Parallelism().set('tensor', nproc=2)
            sage: a = CompFullyAntiSym(ZZ, (1,2,3), 2)
            sage: a[0,1], a[1,2] = 4, 5
            sage: b = CompFullyAntiSym(ZZ, (1,2,3), 2)
            sage: b[0,1], b[0,2] = 2, -3
            sage: s_par = a.__add__(b) ; s_par  # the antisymmetry is kept
            Fully antisymmetric 2-indices components w.r.t. (1, 2, 3)
            sage: s_par[:]
            [ 0  6 -3]
            [-6  0  5]
            [ 3 -5  0]
            sage: s_par == a + b
            True
            sage: from sage.tensor.modules.comp import CompFullySym
            sage: c = CompFullySym(ZZ, (1,2,3), 2)
            sage: c[0,1], c[0,2] = 3, 7
            sage: s_par = a.__add__(c) ; s_par  # the antisymmetry is lost
            2-indices components w.r.t. (1, 2, 3)
            sage: s_par[:]
            [ 0  7  7]
            [-1  0  5]
            [ 7 -5  0]
            sage: s_par == a + c
            True
            sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization
        """
        if isinstance(other, (int, Integer)) and other == 0:
            return +self
        if not isinstance(other, Components):
            raise TypeError("the second argument for the addition must be a " +
                            "an instance of Components")
        if isinstance(other, CompFullyAntiSym):
            if other._frame != self._frame:
                raise ValueError("the two sets of components are not defined " +
                                 "on the same frame")
            if other._nid != self._nid:
                raise ValueError("the two sets of components do not have the " +
                                 "same number of indices")
            if other._sindex != self._sindex:
                raise ValueError("the two sets of components do not have the " +
                                 "same starting index")
            # Initialization of the result to self.copy(), so that there remains
            # only to add other:
            result = self.copy()
            nproc = Parallelism().get('tensor')
            if nproc != 1:
                # Parallel computation
                lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
                ind_list = list(other._comp)
                ind_step = max(1, int(len(ind_list)/nproc/2))
                local_list = lol(ind_list, ind_step)
                # list of input parameters
                listParalInput = [(self, other, ind_part) for ind_part in local_list]

                @parallel(p_iter='multiprocessing', ncpus=nproc)
                def paral_sum(a, b, local_list_ind):
                    partial = []
                    for ind in local_list_ind:
                        partial.append([ind, a[[ind]]+b[[ind]]])
                    return partial

                for ii, val in paral_sum(listParalInput):
                    for jj in val:
                        result[[jj[0]]] = jj[1]

            else:
                # Sequential computation
                for ind, val in other._comp.items():
                    result[[ind]] += val

            return result
        else:
            return CompWithSym.__add__(self, other)

    def interior_product(self, other):
        r"""
        Interior product with another set of fully antisymmetric components.

        The interior product amounts to a contraction over all the `p` indices
        of ``self`` with the first `p` indices of ``other``, assuming that
        the number `q` of indices of ``other`` obeys `q\geq p`.

        .. NOTE::

            ``self.interior_product(other)`` yields the same result as
            ``self.contract(0,..., p-1, other, 0,..., p-1)``
            (cf. :meth:`~sage.tensor.modules.comp.Components.contract`), but
            ``interior_product`` is more efficient, the antisymmetry of ``self``
            being not used to reduce the computation in
            :meth:`~sage.tensor.modules.comp.Components.contract`.

        INPUT:

        - ``other`` -- fully antisymmetric components defined on the same frame
          as ``self`` and with a number of indices at least equal to that of
          ``self``

        OUTPUT:

        - base ring element (case `p=q`) or set of components (case `p<q`)
          resulting from the contraction over all the `p` indices of ``self``
          with the first `p` indices of ``other``

        EXAMPLES:

        Interior product of a set of components ``a`` with ``p`` indices with a
        set of components ``b`` with ``q`` indices on a 4-dimensional vector
        space.

        Case ``p=2`` and ``q=2``::

            sage: from sage.tensor.modules.comp import CompFullyAntiSym
            sage: V = VectorSpace(QQ, 4)
            sage: a = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: a[0,1], a[0,2], a[0,3] = -2, 4, 3
            sage: a[1,2], a[1,3], a[2,3] = 5, -3, 1
            sage: b = CompFullyAntiSym(QQ, V.basis(), 2)
            sage: b[0,1], b[0,2], b[0,3] = 3, -4, 2
            sage: b[1,2], b[1,3], b[2,3] = 2, 5, 1
            sage: c = a.interior_product(b)
            sage: c
            -40
            sage: c == a.contract(0, 1, b, 0, 1)
            True

        Case ``p=2`` and ``q=3``::

            sage: b = CompFullyAntiSym(QQ, V.basis(), 3)
            sage: b[0,1,2], b[0,1,3], b[0,2,3], b[1,2,3] = 3, -4, 2, 5
            sage: c = a.interior_product(b)
            sage: c[:]
            [58, 10, 6, 82]
            sage: c == a.contract(0, 1, b, 0, 1)
            True

        Case ``p=2`` and ``q=4``::

            sage: b = CompFullyAntiSym(QQ, V.basis(), 4)
            sage: b[0,1,2,3] = 5
            sage: c = a.interior_product(b)
            sage: c[:]
            [  0  10  30  50]
            [-10   0  30 -40]
            [-30 -30   0 -20]
            [-50  40  20   0]
            sage: c == a.contract(0, 1, b, 0, 1)
            True

        Case ``p=3`` and ``q=3``::

            sage: a = CompFullyAntiSym(QQ, V.basis(), 3)
            sage: a[0,1,2], a[0,1,3], a[0,2,3], a[1,2,3] = 2, -1, 3, 5
            sage: b = CompFullyAntiSym(QQ, V.basis(), 3)
            sage: b[0,1,2], b[0,1,3], b[0,2,3], b[1,2,3] = -2, 1, 4, 2
            sage: c = a.interior_product(b)
            sage: c
            102
            sage: c == a.contract(0, 1, 2, b, 0, 1, 2)
            True

        Case ``p=3`` and ``q=4``::

            sage: b = CompFullyAntiSym(QQ, V.basis(), 4)
            sage: b[0,1,2,3] = 5
            sage: c = a.interior_product(b)
            sage: c[:]
            [-150, 90, 30, 60]
            sage: c == a.contract(0, 1, 2, b, 0, 1, 2)
            True

        Case ``p=4`` and ``q=4``::

            sage: a = CompFullyAntiSym(QQ, V.basis(), 4)
            sage: a[0,1,2,3] = 3
            sage: c = a.interior_product(b)
            sage: c
            360
            sage: c == a.contract(0, 1, 2, 3, b, 0, 1, 2, 3)
            True
        """
        from sage.arith.misc import factorial
        # Sanity checks:
        if not isinstance(other, CompFullyAntiSym):
            raise TypeError("{} is not a fully antisymmetric ".format(other) +
                            "set of components")
        if other._frame != self._frame:
            raise ValueError("The {} are not defined on the ".format(other) +
                             "same frame as the {}".format(self))
        if other._nid < self._nid:
            raise ValueError("The {} have less indices than ".format(other) +
                             "the {}".format(self))
        # Number of indices of the result:
        res_nid = other._nid - self._nid
        #
        # Case of a scalar result
        #
        if res_nid == 0:
            res = 0
            for ind in self.non_redundant_index_generator():
                res += self[[ind]] * other[[ind]]
            return factorial(self._nid)*res
        #
        # Case of component result
        #
        if res_nid == 1:
            res = Components(self._ring, self._frame, res_nid,
                             start_index=self._sindex,
                             output_formatter=self._output_formatter)
        else:
            res = CompFullyAntiSym(self._ring, self._frame, res_nid,
                                   start_index=self._sindex,
                                   output_formatter=self._output_formatter)
        factorial_s = factorial(self._nid)
        for ind in res.non_redundant_index_generator():
            sm = 0
            for ind_s, cmp_s in self._comp.items():
                ind_o = ind_s + ind
                sm += cmp_s * other[[ind_o]]
            res[[ind]] = factorial_s*sm
        return res

#******************************************************************************


class KroneckerDelta(CompFullySym):
    r"""
    Kronecker delta `\delta_{ij}`.

    INPUT:

    - ``ring`` -- commutative ring in which each component takes its value
    - ``frame`` -- frame with respect to which the components are defined;
      whatever type ``frame`` is, it should have some method ``__len__()``
      implemented, so that ``len(frame)`` returns the dimension, i.e. the size
      of a single index range
    - ``start_index`` -- (default: 0) first value of a single index;
      accordingly a component index i must obey
      ``start_index <= i <= start_index + dim - 1``, where ``dim = len(frame)``.
    - ``output_formatter`` -- (default: ``None``) function or unbound
      method called to format the output of the component access
      operator ``[...]`` (method ``__getitem__``); ``output_formatter`` must
      take 1 or 2 arguments: the first argument must be an instance of
      ``ring`` and the second one, if any, some format specification

    EXAMPLES:

    The Kronecker delta on a 3-dimensional space::

        sage: from sage.tensor.modules.comp import KroneckerDelta
        sage: V = VectorSpace(QQ,3)
        sage: d = KroneckerDelta(QQ, V.basis()) ; d
        Kronecker delta of size 3x3
        sage: d[:]
        [1 0 0]
        [0 1 0]
        [0 0 1]

    One can read, but not set, the components of a Kronecker delta::

        sage: d[1,1]
        1
        sage: d[1,1] = 2
        Traceback (most recent call last):
        ...
        TypeError: the components of a Kronecker delta cannot be changed

    Examples of use with output formatters::

        sage: d = KroneckerDelta(QQ, V.basis(), output_formatter=Rational.numerical_approx)
        sage: d[:]  # default format (53 bits of precision)
        [ 1.00000000000000 0.000000000000000 0.000000000000000]
        [0.000000000000000  1.00000000000000 0.000000000000000]
        [0.000000000000000 0.000000000000000  1.00000000000000]
        sage: d[:,10] # format = 10 bits of precision
        [ 1.0 0.00 0.00]
        [0.00  1.0 0.00]
        [0.00 0.00  1.0]
        sage: d = KroneckerDelta(QQ, V.basis(), output_formatter=str)
        sage: d[:]
        [['1', '0', '0'], ['0', '1', '0'], ['0', '0', '1']]
    """
    def __init__(self, ring, frame, start_index=0, output_formatter=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.comp import KroneckerDelta
            sage: d = KroneckerDelta(ZZ, (1,2,3))
            sage: TestSuite(d).run()
        """
        CompFullySym.__init__(self, ring, frame, 2, start_index,
                              output_formatter)
        for i in range(self._sindex, self._dim + self._sindex):
            self._comp[(i,i)] = self._ring(1)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.comp import KroneckerDelta
            sage: KroneckerDelta(ZZ, (1,2,3))
            Kronecker delta of size 3x3
        """
        n = str(self._dim)
        return "Kronecker delta of size " + n + "x" + n

    def __setitem__(self, args, value):
        r"""
        Should not be used (the components of a Kronecker delta are constant).

        EXAMPLES::

            sage: from sage.tensor.modules.comp import KroneckerDelta
            sage: d = KroneckerDelta(ZZ, (1,2,3))
            sage: d.__setitem__((0,0), 1)
            Traceback (most recent call last):
            ...
            TypeError: the components of a Kronecker delta cannot be changed
        """
        raise TypeError("the components of a Kronecker delta cannot be changed")
