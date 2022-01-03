"""
Matrix/Vector-Valued Linear Functions: Parents

In Sage, matrices assume that the base is a ring. Hence, we cannot
construct matrices whose entries are linear functions in Sage. Really,
they should be thought of as the tensor product of the `R`-module of
linear functions and the `R`-vector/matrix space, with the latter viewed
as an `R`-module (`R` is usually ``QQ`` or ``RDF`` for our purposes).

You should not construct any tensor products by calling the parent
directly. This is also why none of the classes are imported in the
global namespace. The come into play whenever you have vector or
matrix MIP linear expressions/constraints. The intended way to
construct them is implicitly by acting with vectors or matrices on
linear functions. For example::

    sage: mip.<x> = MixedIntegerLinearProgram('ppl')   # base ring is QQ
    sage: 3 + x[0] + 2*x[1]            # a linear function
    3 + x_0 + 2*x_1
    sage: x[0] * vector([3,4]) + 1     # vector linear function
    (1, 1) + (3, 4)*x_0
    sage: x[0] * matrix([[3,1],[4,0]]) + 1   # matrix linear function
    [1 + 3*x_0 x_0]
    [4*x_0     1  ]

Internally, all linear functions are stored as a dictionary whose

* keys are the index of the linear variable (and -1 for the constant
  term)

* values are the coefficient of that variable. That is, a number for
  linear functions, a vector for vector-valued functions, etc.

The entire dictionary can be accessed with the
:meth:`~sage.numerical.linear_tensor_element.LinearTensor.dict`
method. For convenience, you can also retrieve a single coefficient
with
:meth:`~sage.numerical.linear_tensor_element.LinearTensor.coefficient`. For
example::

    sage: mip.<b> = MixedIntegerLinearProgram()
    sage: f_scalar = (3 + b[7] + 2*b[9]);  f_scalar
    3 + x_0 + 2*x_1
    sage: f_scalar.dict()
    {-1: 3.0, 0: 1.0, 1: 2.0}
    sage: f_scalar.dict()[1]
    2.0
    sage: f_scalar.coefficient(b[9])
    2.0
    sage: f_scalar.coefficient(1)
    2.0

    sage: f_vector = b[7] * vector([3,4]) + 1;  f_vector
    (1.0, 1.0) + (3.0, 4.0)*x_0
    sage: f_vector.coefficient(-1)
    (1.0, 1.0)
    sage: f_vector.coefficient(b[7])
    (3.0, 4.0)
    sage: f_vector.coefficient(0)
    (3.0, 4.0)
    sage: f_vector.coefficient(1)
    (0.0, 0.0)

    sage: f_matrix = b[7] * matrix([[0,1], [2,0]]) + b[9] - 3;  f_matrix
    [-3 + x_1 x_0     ]
    [2*x_0    -3 + x_1]
    sage: f_matrix.coefficient(-1)
    [-3.0  0.0]
    [ 0.0 -3.0]
    sage: f_matrix.coefficient(0)
    [0.0 1.0]
    [2.0 0.0]
    sage: f_matrix.coefficient(1)
    [1.0 0.0]
    [0.0 1.0]

Just like :mod:`sage.numerical.linear_functions`, (in)equalities
become symbolic inequalities. See
:mod:`~sage.numerical.linear_tensor_constraints` for details.

.. NOTE::

    For brevity, we just use ``LinearTensor`` in class names. It is
    understood that this refers to the above tensor product
    construction.
"""

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from copy import copy

from sage.structure.parent import Parent
from sage.misc.cachefunc import cached_function
from sage.numerical.linear_functions import is_LinearFunction, LinearFunctionsParent_class
from sage.numerical.linear_tensor_element import LinearTensor


#*****************************************************************************
#
# Utility functions to test that something is a linear function / constraint
#
#*****************************************************************************

def is_LinearTensor(x):
    """
    Test whether ``x`` is a tensor product of linear functions with a
    free module.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: p = MixedIntegerLinearProgram()
        sage: x = p.new_variable(nonnegative=False)
        sage: from sage.numerical.linear_tensor import is_LinearTensor
        sage: is_LinearTensor(x[0] - 2*x[2])
        False
        sage: is_LinearTensor('a string')
        False
    """
    return isinstance(x, LinearTensor)


#*****************************************************************************
#
# Factory functions for the parents to ensure uniqueness
#
#*****************************************************************************

@cached_function
def LinearTensorParent(free_module_parent, linear_functions_parent):
    """
    Return the parent for the tensor product over the common ``base_ring``.

    The output is cached, so only a single parent is ever constructed
    for a given base ring.

    INPUT:
    
    - ``free_module_parent`` -- module. A free module, like vector or
      matrix space.

    - ``linear_functions_parent`` -- linear functions. The linear
      functions parent.

    OUTPUT:

    The parent of the tensor product of a free module and linear
    functions over a common base ring.

    EXAMPLES::

        sage: from sage.numerical.linear_functions import LinearFunctionsParent
        sage: from sage.numerical.linear_tensor import LinearTensorParent
        sage: LinearTensorParent(QQ^3, LinearFunctionsParent(QQ))
        Tensor product of Vector space of dimension 3 over Rational Field and Linear functions over Rational Field

        sage: LinearTensorParent(ZZ^3, LinearFunctionsParent(QQ))
        Traceback (most recent call last):
        ...
        ValueError: base rings must match
    """
    if free_module_parent.base_ring() != linear_functions_parent.base_ring():
        raise ValueError('base rings must match')
    if not isinstance(linear_functions_parent, LinearFunctionsParent_class):
        raise TypeError('linear_functions_parent must be a parent of linear functions')
    return LinearTensorParent_class(free_module_parent, linear_functions_parent)



#*****************************************************************************
#
# Parent of linear functions tensored with a free module
#
#*****************************************************************************

class LinearTensorParent_class(Parent):
    r"""
    The parent for all linear functions over a fixed base ring.

    .. warning::

        You should use :func:`LinearTensorParent` to construct
        instances of this class.

    INPUT/OUTPUT:

    See :func:`LinearTensorParent`

    EXAMPLES::

        sage: from sage.numerical.linear_tensor import LinearTensorParent_class
        sage: LinearTensorParent_class
        <class 'sage.numerical.linear_tensor.LinearTensorParent_class'>
    """
    Element = LinearTensor

    def __init__(self, free_module, linear_functions):
        """
        The Python constructor

        INPUT/OUTPUT:

        See :func:`LinearTensorParent`

        TESTS::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LinearFunctionsParent(RDF).tensor(RDF^2)
            Tensor product of Vector space of dimension 2 over Real Double
            Field and Linear functions over Real Double Field
        """
        self._free_module = free_module
        self._linear_functions = linear_functions
        base_ring = linear_functions.base_ring()
        from sage.categories.modules_with_basis import ModulesWithBasis
        Parent.__init__(self, base=base_ring, category=ModulesWithBasis(base_ring))

    def free_module(self):
        """
        Return the linear functions.

        See also :meth:`free_module`.

        OUTPUT:

        Parent of the linear functions, one of the factors in the
        tensor product construction.

        EXAMPLES::

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: lt = x[0] * vector(RDF, [1,2])
            sage: lt.parent().free_module()
            Vector space of dimension 2 over Real Double Field
            sage: lt.parent().free_module() is vector(RDF, [1,2]).parent()
            True
        """
        return self._free_module

    def is_vector_space(self):
        """
        Return whether the free module is a vector space.

        OUTPUT:

        Boolean. Whether the :meth:`free_module` factor in the tensor
        product is a vector space.

        EXAMPLES::

            sage: mip = MixedIntegerLinearProgram()
            sage: LF = mip.linear_functions_parent()
            sage: LF.tensor(RDF^2).is_vector_space()
            True
            sage: LF.tensor(RDF^(2,2)).is_vector_space()       
            False
        """
        from sage.modules.free_module import is_FreeModule
        return is_FreeModule(self.free_module())

    def is_matrix_space(self):
        """
        Return whether the free module is a matrix space.

        OUTPUT:

        Boolean. Whether the :meth:`free_module` factor in the tensor
        product is a matrix space.

        EXAMPLES::

            sage: mip = MixedIntegerLinearProgram()
            sage: LF = mip.linear_functions_parent()
            sage: LF.tensor(RDF^2).is_matrix_space()
            False
            sage: LF.tensor(RDF^(2,2)).is_matrix_space()
            True
        """
        from sage.matrix.matrix_space import is_MatrixSpace
        return is_MatrixSpace(self.free_module())

    def linear_functions(self):
        """
        Return the linear functions.

        See also :meth:`free_module`.

        OUTPUT:

        Parent of the linear functions, one of the factors in the
        tensor product construction.

        EXAMPLES::

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: lt = x[0] * vector([1,2])
            sage: lt.parent().linear_functions()
            Linear functions over Real Double Field
            sage: lt.parent().linear_functions() is mip.linear_functions_parent()
            True
        """
        return self._linear_functions

    def _repr_(self):
        """
        Return a string representation
        
        OUTPUT:

        String.

        EXAMPLES::

            sage: MixedIntegerLinearProgram().linear_functions_parent()
            Linear functions over Real Double Field
        """
        return 'Tensor product of {0} and {1}'.format(self.free_module(), self.linear_functions())

    def _convert_constant(self, m):
        """
        Convert ``m`` to a constant free module element.

        OUTPUT:

        A :meth:`free_module` element.

        EXAMPLES::

            sage: mip = MixedIntegerLinearProgram()
            sage: LF = mip.linear_functions_parent()
            sage: LF.tensor(RDF^2)._convert_constant(42)
            (42.0, 42.0)
            sage: LF.tensor(RDF^(2,2))._convert_constant(42)
            [42.0  0.0]
            [ 0.0 42.0]
        """
        M = self.free_module()
        m = M.base_ring()(m)
        if self.is_matrix_space():
            # Turn constants into diagonal matrices
            m_matrix = copy(M.zero_matrix())
            for i in range(min(M.ncols(), M.nrows())):
                m_matrix[i, i] = m
            m_matrix.set_immutable()
            return m_matrix
        elif self.is_vector_space():
            # Turn constants into vectors with all entries equal
            m_vector = M([m] * M.degree())
            return m_vector
        else:
            return M(m)
        
    def _element_constructor_(self, x):
        """
        Construct a :class:`LinearTensor` from ``x``.

        INPUT:

        - ``x`` -- anything that defines a
          :class:`~sage.numerical.linear_tensor_element.LinearTensor`. See
          examples.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: LT = p.linear_functions_parent().tensor(RDF^2)
            sage: LT._element_constructor_(123)
            (123.0, 123.0)

        Construct from dict with type conversion to RDF vector::

            sage: LT({1:[1, 2]})    # indirect doctest
            (1.0, 2.0)*x_1
            sage: type(_)
            <class 'sage.numerical.linear_tensor_element.LinearTensor'>

        Construct from scalar:

            sage: LT(123)    # indirect doctest
            (123.0, 123.0)
        
        Similar, over ``QQ`` and with matrices instead of vectors::

            sage: p_QQ = MixedIntegerLinearProgram(solver='ppl')
            sage: LT_QQ = p_QQ.linear_functions_parent().tensor(QQ^(2, 2))
            sage: LT_QQ({-1:[[1/2, 1/3], [2, 3]], 2:[[3/4, 1/4], [0, 0]]})
            [1/2 + 3/4*x_2 1/3 + 1/4*x_2]
            [2             3            ]
            sage: LT_QQ(42.1)
            [421/10 0     ]
            [0      421/10]

        Construct from a linear function::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LF_ZZ = LinearFunctionsParent(ZZ)
            sage: lf = LF_ZZ({-1:3, 1:2, 3:1})
            sage: LT(lf)
            (3.0, 3.0) + (2.0, 2.0)*x_1 + (1.0, 1.0)*x_3
        """
        M = self.free_module()
        if is_LinearTensor(x):
            x = x.dict()
        elif is_LinearFunction(x):
            x = dict([key, self._convert_constant(value)] for key, value in x.dict().items())
        elif isinstance(x, dict):
            x = dict([int(key), M(value)] for key, value in x.items())
        else:
            try:
                x = {-1: M(x)}
            except (TypeError, ValueError):
                x_R = M.base_ring()(x)
                x = {-1: self._convert_constant(x_R)}
        return self.element_class(self, x)

    def _coerce_map_from_(self, R):
        """
        Allow coercion of scalars into tensors.

        INPUT:

        - ``R`` -- a ring.

        OUTPUT:

        Boolean. Whether there is a coercion map.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: parent = p.linear_functions_parent()
            sage: parent.coerce(int(2))
            2
            sage: parent._coerce_map_from_(int)
            True

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: LF = mip.linear_functions_parent()
            sage: LT = LF.tensor(RDF^3)
            sage: LT.has_coerce_map_from(LF)
            True
        """
        if self.free_module().has_coerce_map_from(R):
            return True
        if self.linear_functions().has_coerce_map_from(R):
            return True
        return False

    def _an_element_(self):
        """
        Returns an element

        OUTPUT:

        A linear function tensored with a free module.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram().linear_functions_parent().tensor(RDF^2)
            sage: p._an_element_()
            (1.0, 0.0) + (5.0, 0.0)*x_2 + (7.0, 0.0)*x_5
            sage: p.an_element()   # indirect doctest
            (1.0, 0.0) + (5.0, 0.0)*x_2 + (7.0, 0.0)*x_5
        """
        m = self.free_module().an_element()
        return self._element_constructor_({-1:m, 2:5*m, 5:7*m})



