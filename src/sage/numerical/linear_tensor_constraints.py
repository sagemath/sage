"""
Constraints on Linear Functions Tensored with a Free Module

Here is an example of a vector-valued linear function::

    sage: mip.<x> = MixedIntegerLinearProgram('ppl')   # base ring is QQ
    sage: x[0] * vector([3,4]) + 1     # vector linear function
    (1, 1) + (3, 4)*x_0

Just like :mod:`~sage.numerical.linear_functions`, (in)equalities
become symbolic inequalities::

    sage: 3 + x[0] + 2*x[1] <= 10
    3 + x_0 + 2*x_1 <= 10
    sage: x[0] * vector([3,4]) + 1 <= 10
    (1, 1) + (3, 4)*x_0 <= (10, 10)
    sage: x[0] * matrix([[0,0,1],[0,1,0],[1,0,0]]) + x[1] * identity_matrix(3) >= 0
    [0 0 0]    [x_1 0         x_0]
    [0 0 0] <= [0   x_0 + x_1 0  ]
    [0 0 0]    [x_0 0         x_1]
"""
# ****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.misc.cachefunc import cached_function


#*****************************************************************************
#
# Utility functions to test that something is a linear function / constraint
#
#*****************************************************************************

def is_LinearTensorConstraint(x):
    """
    Test whether ``x`` is a constraint on module-valued linear functions.

    INPUT:

    - ``x`` -- anything

    OUTPUT: boolean

    EXAMPLES::

        sage: mip.<x> = MixedIntegerLinearProgram()
        sage: vector_ieq = (x[0] * vector([1,2]) <= x[1] * vector([2,3]))
        sage: from sage.numerical.linear_tensor_constraints import is_LinearTensorConstraint
        sage: is_LinearTensorConstraint(vector_ieq)
        doctest:warning...
        DeprecationWarning: The function is_LinearTensorConstraint is deprecated;
        use 'isinstance(..., LinearTensorConstraint)' instead.
        See https://github.com/sagemath/sage/issues/38184 for details.
        True
        sage: is_LinearTensorConstraint('a string')
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(38184,
                "The function is_LinearTensorConstraint is deprecated; "
                "use 'isinstance(..., LinearTensorConstraint)' instead.")
    return isinstance(x, LinearTensorConstraint)


#*****************************************************************************
#
# Factory functions for the parents to ensure uniqueness
#
#*****************************************************************************

@cached_function
def LinearTensorConstraintsParent(linear_functions_parent):
    """
    Return the parent for linear functions over ``base_ring``.

    The output is cached, so only a single parent is ever constructed
    for a given base ring.

    INPUT:

    - ``linear_functions_parent`` -- a
      :class:`~sage.numerical.linear_functions.LinearFunctionsParent_class`. The
      type of linear functions that the constraints are made out of.

    OUTPUT: the parent of the linear constraints with the given linear functions

    EXAMPLES::

        sage: from sage.numerical.linear_functions import LinearFunctionsParent
        sage: from sage.numerical.linear_tensor import LinearTensorParent
        sage: from sage.numerical.linear_tensor_constraints import \
        ....:     LinearTensorConstraintsParent, LinearTensorConstraintsParent
        sage: LF = LinearFunctionsParent(QQ)
        sage: LT = LinearTensorParent(QQ^2, LF)
        sage: LinearTensorConstraintsParent(LT)
        Linear constraints in the tensor product of Vector space of dimension 2
        over Rational Field and Linear functions over Rational Field
    """
    return LinearTensorConstraintsParent_class(linear_functions_parent)


# ****************************************************************************
#
# Elements of linear tensor constraints
#
# ****************************************************************************

class LinearTensorConstraint(Element):
    """
    Formal constraint involving two module-valued linear functions.

    .. NOTE::

        In the code, we use "linear tensor" as abbreviation for the
        tensor product (over the common base ring) of a :mod:`linear
        function <sage.numerical.linear_functions>` and a free module
        like a vector/matrix space.

    .. warning::

        This class has no reason to be instantiated by the user, and
        is meant to be used by instances of
        :class:`MixedIntegerLinearProgram`.

    INPUT:

    - ``parent`` -- the parent, a
      :class:`LinearTensorConstraintsParent_class`

    - ``lhs``, ``rhs`` -- two
      :class:`sage.numerical.linear_tensor_element.LinearTensor`. The
      left and right hand side of the constraint (in)equality.

    - ``equality`` -- boolean (default: ``False``); whether the
      constraint is an equality.  If ``False``, it is a ``<=``
      inequality.

    EXAMPLES::

        sage: mip.<b> = MixedIntegerLinearProgram()
        sage: (b[2]+2*b[3]) * vector([1,2]) <= b[8] * vector([2,3]) - 5
        (1.0, 2.0)*x_0 + (2.0, 4.0)*x_1 <= (-5.0, -5.0) + (2.0, 3.0)*x_2
    """

    def __init__(self, parent, lhs, rhs, equality):
        r"""
        Constructor for ``LinearTensorConstraint``.

        INPUT:

        See :class:`LinearTensorConstraint`.

        EXAMPLES::

            sage: mip.<b> = MixedIntegerLinearProgram()
            sage: b[2] * vector([1,2]) + 2*b[3] <= 0
            (1.0, 2.0)*x_0 + (2.0, 2.0)*x_1 <= (0.0, 0.0)
        """
        super().__init__(parent)
        self._lhs = lhs
        self._rhs = rhs
        self._equality = equality

    def is_equation(self):
        """
        Whether the constraint is a chained equation.

        OUTPUT: boolean

        EXAMPLES::

            sage: mip.<b> = MixedIntegerLinearProgram()
            sage: (b[0] * vector([1,2]) == 0).is_equation()
            True
            sage: (b[0] * vector([1,2]) >= 0).is_equation()
            False
        """
        return self._equality

    def is_less_or_equal(self):
        """
        Whether the constraint is a chained less-or_equal inequality.

        OUTPUT: boolean

        EXAMPLES::

            sage: mip.<b> = MixedIntegerLinearProgram()
            sage: (b[0] * vector([1,2]) == 0).is_less_or_equal()
            False
            sage: (b[0] * vector([1,2]) >= 0).is_less_or_equal()
            True
        """
        return not self._equality

    def lhs(self):
        """
        Return the left side of the (in)equality.

        OUTPUT:

        Instance of
        :class:`sage.numerical.linear_tensor_element.LinearTensor`. A
        linear function valued in a free module.

        EXAMPLES::

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: (x[0] * vector([1,2]) == 0).lhs()
            (1.0, 2.0)*x_0
        """
        return self._lhs

    def rhs(self):
        """
        Return the right side of the (in)equality.

        OUTPUT:

        Instance of
        :class:`sage.numerical.linear_tensor_element.LinearTensor`. A
        linear function valued in a free module.

        EXAMPLES::

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: (x[0] * vector([1,2]) == 0).rhs()
            (0.0, 0.0)
        """
        return self._rhs

    def _ascii_art_(self):
        """
        Return Ascii Art.

        OUTPUT:

        Ascii art of the constraint (in)equality.

        EXAMPLES::

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: ascii_art(x[0] * vector([1,2]) >= 0)
            (0.0, 0.0) <= (1.0, 2.0)*x_0
            sage: ascii_art(x[0] * matrix([[1,2],[3,4]]) >= 0)
            [0 0] <= [x_0   2*x_0]
            [0 0]    [3*x_0 4*x_0]
        """
        from sage.typeset.ascii_art import AsciiArt

        def matrix_art(m):
            lines = str(m).splitlines()
            return AsciiArt(lines, baseline=len(lines) // 2)
        comparator = AsciiArt([' == ' if self.is_equation() else ' <= '])
        return matrix_art(self.lhs()) + comparator + matrix_art(self.rhs())

    def _repr_(self):
        r"""
        Return a string representation of the constraint.

        OUTPUT: string

        EXAMPLES::

            sage: mip.<b> = MixedIntegerLinearProgram()
            sage: b[3] * vector([1,2]) <= (b[8] + 9) * vector([2,3])
            (1.0, 2.0)*x_0 <= (18.0, 27.0) + (2.0, 3.0)*x_1
            sage: b[3] * vector([1,2]) == (b[8] + 9) * vector([2,3])
            (1.0, 2.0)*x_0 == (18.0, 27.0) + (2.0, 3.0)*x_1
            sage: b[0] * identity_matrix(3) == 0
            [x_2 0   0  ]    [0 0 0]
            [0   x_2 0  ] == [0 0 0]
            [0   0   x_2]    [0 0 0]
        """
        if self.parent().linear_tensors().is_matrix_space():
            return str(self._ascii_art_())
        comparator = (' == ' if self.is_equation() else ' <= ')
        return str(self.lhs()) + comparator + str(self.rhs())


#*****************************************************************************
#
# Parent of linear constraints
#
#*****************************************************************************

class LinearTensorConstraintsParent_class(Parent):
    """
    Parent for :class:`LinearTensorConstraint`.

    .. warning::

        This class has no reason to be instantiated by the user, and
        is meant to be used by instances of
        :class:`MixedIntegerLinearProgram`. Also, use the
        :func:`LinearTensorConstraintsParent` factory function.

    INPUT/OUTPUT: see :func:`LinearTensorConstraintsParent`

    EXAMPLES::

        sage: p = MixedIntegerLinearProgram()
        sage: LT = p.linear_functions_parent().tensor(RDF^2);  LT
        Tensor product of Vector space of dimension 2 over Real Double
        Field and Linear functions over Real Double Field
        sage: from sage.numerical.linear_tensor_constraints import LinearTensorConstraintsParent
        sage: LTC = LinearTensorConstraintsParent(LT);  LTC
        Linear constraints in the tensor product of Vector space of
        dimension 2 over Real Double Field and Linear functions over
        Real Double Field
        sage: type(LTC)
        <class 'sage.numerical.linear_tensor_constraints.LinearTensorConstraintsParent_class'>
    """
    Element = LinearTensorConstraint

    def __init__(self, linear_tensor_parent):
        """
        The Python constructor.

        INPUT:

        - ``linear_tensor_parent`` -- instance of
          :class:`LinearTensorParent_class`

        TESTS::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LF = LinearFunctionsParent(RDF)
            sage: from sage.numerical.linear_tensor import LinearTensorParent
            sage: LT = LinearTensorParent(RDF^2, LF)
            sage: from sage.numerical.linear_tensor_constraints import LinearTensorConstraintsParent
            sage: LinearTensorConstraintsParent(LT)
            Linear constraints in the tensor product of Vector space of
            dimension 2 over Real Double Field and Linear functions over
            Real Double Field
        """
        Parent.__init__(self)
        self._LT = linear_tensor_parent
        self._LF = linear_tensor_parent.linear_functions()

    def linear_tensors(self):
        """
        Return the parent for the linear functions.

        OUTPUT: instance of :class:`sage.numerical.linear_tensor.LinearTensorParent_class`

        EXAMPLES::

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: ieq = (x[0] * vector([1,2]) >= 0)
            sage: ieq.parent().linear_tensors()
            Tensor product of Vector space of dimension 2 over Real Double
            Field and Linear functions over Real Double Field
        """
        return self._LT

    def linear_functions(self):
        """
        Return the parent for the linear functions.

        OUTPUT: instance of :class:`sage.numerical.linear_functions.LinearFunctionsParent_class`

        EXAMPLES::

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: ieq = (x[0] * vector([1,2]) >= 0)
            sage: ieq.parent().linear_functions()
            Linear functions over Real Double Field
        """
        return self._LF

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: ieq = (x[0] * vector([1,2]) >= 0)
            sage: ieq.parent()    # indirect doctests
            Linear constraints in the tensor product of Vector space of
            dimension 2 over Real Double Field and Linear functions over
            Real Double Field
        """
        return 'Linear constraints in the tensor product of {0} and {1}'.format(
            self.linear_tensors().free_module(), self.linear_functions())

    def _element_constructor_(self, left, right, equality):
        """
        Construct a :class:`LinearConstraint`.

        INPUT:

        - ``left`` -- a :class:`LinearTensor`, or something that can
          be converted into one, a list/tuple of
          :class:`LinearTensor`, or an existing
          :class:`LinearTensorConstraint`.

        - ``right`` -- a :class:`LinearTensor` or ``None``
          (default)

        - ``equality`` -- boolean; whether to
          construct an equation or a less-or-equal inequality

        OUTPUT: the :class:`LinearTensorConstraint` constructed from the input data

        EXAMPLES::

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: ieq = (x[0] * vector([1,2]) >= 0)
            sage: LTC = ieq.parent()
            sage: LTC._element_constructor_(1, 2, True)
            (1.0, 1.0) == (2.0, 2.0)
            sage: LTC(x[0], x[1], False)
            (1.0, 1.0)*x_0 <= (1.0, 1.0)*x_1
            sage: type(_)
            <class 'sage.numerical.linear_tensor_constraints.LinearTensorConstraintsParent_class.element_class'>
        """
        LT = self.linear_tensors()
        left = LT(left)
        right = LT(right)
        equality = bool(equality)
        return self.element_class(self, left, right, equality)

    def _an_element_(self):
        """
        Return an element.

        EXAMPLES::

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: ieq = (x[0] * vector([1,2]) >= 0)
            sage: ieq.parent().an_element()    # indirect doctest
            (0.0, 0.0) <= (1.0, 0.0) + (5.0, 0.0)*x_2 + (7.0, 0.0)*x_5
        """
        LT = self.linear_tensors()
        return LT.an_element() >= 0
