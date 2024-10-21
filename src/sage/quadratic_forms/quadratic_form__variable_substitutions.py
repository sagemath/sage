# sage_setup: distribution = sagemath-modules
"""
Variable Substitution, Multiplication, Division, Scaling
"""
# ****************************************************************************
#       Copyright (C) 2007 William Stein and Jonathan Hanke
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************


def swap_variables(self, r, s, in_place=False):
    r"""
    Switch the variables `x_r` and `x_s` in the quadratic form
    (replacing the original form if the ``in_place`` flag is True).

    INPUT:

    - ``r``, ``s`` -- integers `\geq 0`

    OUTPUT:

    a :class:`QuadraticForm` (by default, otherwise none)

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 4, range(1,11))
        sage: Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 2 3 4 ]
        [ * 5 6 7 ]
        [ * * 8 9 ]
        [ * * * 10 ]

        sage: Q.swap_variables(0,2)
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 8 6 3 9 ]
        [ * 5 2 7 ]
        [ * * 1 4 ]
        [ * * * 10 ]

        sage: Q.swap_variables(0,2).swap_variables(0,2)
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 2 3 4 ]
        [ * 5 6 7 ]
        [ * * 8 9 ]
        [ * * * 10 ]
    """
    if not in_place:
        Q = self.parent()(self.base_ring(), self.dim(), self.coefficients())
        Q.swap_variables(r, s, in_place=True)
        return Q

    # Switch diagonal elements
    tmp = self[r, r]
    self[r, r] = self[s, s]
    self[s, s] = tmp

    # Switch off-diagonal elements
    for i in range(self.dim()):
        if (i != r) and (i != s):
            tmp = self[r, i]
            self[r, i] = self[s, i]
            self[s, i] = tmp


def multiply_variable(self, c, i, in_place=False):
    r"""
    Replace the variables `x_i` by `c\cdot x_i` in the quadratic form
    (replacing the original form if the ``in_place`` flag is True).

    Here `c` must be an element of the base ring defining the
    quadratic form.

    INPUT:

    - ``c`` -- an element of ``self.base_ring()``

    - ``i`` -- integer `\geq 0`

    OUTPUT: a :class:`QuadraticForm` (by default, otherwise none)

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,9,5,7])
        sage: Q.multiply_variable(5, 0)
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 25 0 0 0 ]
        [ * 9 0 0 ]
        [ * * 5 0 ]
        [ * * * 7 ]
    """
    if not in_place:
        Q = self.parent()(self.base_ring(), self.dim(), self.coefficients())
        Q.multiply_variable(c, i, in_place=True)
        return Q

    # Stretch the diagonal element
    tmp = c * c * self[i, i]
    self[i, i] = tmp

    # Switch off-diagonal elements
    for k in range(self.dim()):
        if (k != i):
            tmp = c * self[k, i]
            self[k, i] = tmp


def divide_variable(self, c, i, in_place=False):
    r"""
    Replace the variables `x_i` by `(x_i)/c` in the quadratic form
    (replacing the original form if the ``in_place`` flag is True).

    Here `c` must be an element of the base ring defining the
    quadratic form, and the division must be defined in the base
    ring.

    INPUT:

    - ``c`` -- an element of ``self.base_ring()``

    - ``i`` -- integer `\geq 0`

    OUTPUT:

    a :class:`QuadraticForm` (by default, otherwise none)

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,9,5,7])
        sage: Q.divide_variable(3, 1)
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 0 0 0 ]
        [ * 1 0 0 ]
        [ * * 5 0 ]
        [ * * * 7 ]
    """
    if not in_place:
        Q = self.parent()(self.base_ring(), self.dim(), self.coefficients())
        Q.divide_variable(c, i, in_place=True)
        return Q

    # Stretch the diagonal element
    tmp = self[i, i] / (c * c)
    self[i, i] = tmp

    # Switch off-diagonal elements
    for k in range(self.dim()):
        if (k != i):
            tmp = self[k, i] / c
            self[k, i] = tmp


def scale_by_factor(self, c, change_value_ring_flag=False):
    """
    Scale the values of the quadratic form by the number `c`, if
    this is possible while still being defined over its base ring.

    If the flag is set to true, then this will alter the value ring
    to be the field of fractions of the original ring (if necessary).

    INPUT:

    - ``c`` -- a scalar in the fraction field of the value ring of the form

    OUTPUT: a quadratic form of the same dimension

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [3,9,18,27])
        sage: Q.scale_by_factor(3)
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 9 0 0 0 ]
        [ * 27 0 0 ]
        [ * * 54 0 ]
        [ * * * 81 ]

        sage: Q.scale_by_factor(1/3)
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 0 0 0 ]
        [ * 3 0 0 ]
        [ * * 6 0 ]
        [ * * * 9 ]
    """
    # Try to scale the coefficients while staying in the ring of values.
    new_coeff_list = [x * c for x in self.coefficients()]

    # Check if we can preserve the value ring and return result. -- USE THE BASE_RING FOR NOW...
    R = self.base_ring()
    try:
        list2 = [R(x) for x in new_coeff_list]
        Q = self.parent()(R, self.dim(), list2)
        return Q
    except ValueError:
        if not change_value_ring_flag:
            raise TypeError("we could not rescale the lattice in this way and preserve its defining ring")
        else:
            raise RuntimeError("this code is not tested by current doctests")
            F = R.fraction_field()
            list2 = [F(x) for x in new_coeff_list]
            Q = self.parent()(self.dim(), F, list2, R)  # DEFINE THIS!  IT WANTS TO SET THE EQUIVALENCE RING TO R, BUT WITH COEFFS IN F.
            # Q.set_equivalence_ring(R)
            return Q


def extract_variables(QF, var_indices):
    r"""
    Extract the variables (in order) whose indices are listed in
    ``var_indices``, to give a new quadratic form.

    INPUT:

    - ``var_indices`` -- list of integers `\geq 0`

    OUTPUT: a :class:`QuadraticForm`

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 4, range(10)); Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 0 1 2 3 ]
        [ * 4 5 6 ]
        [ * * 7 8 ]
        [ * * * 9 ]
        sage: Q.extract_variables([1,3])
        Quadratic form in 2 variables over Integer Ring with coefficients:
        [ 4 6 ]
        [ * 9 ]
    """
    m = len(var_indices)
    return QF.parent()(QF.base_ring(), m,
                       [QF[var_indices[i], var_indices[j]]
                        for i in range(m)
                        for j in range(i, m)])


def elementary_substitution(self, c, i, j, in_place=False):  # CHECK THIS!!!
    r"""
    Perform the substitution `x_i \longmapsto x_i + c\cdot x_j` (replacing the
    original form if the ``in_place`` flag is True).

    INPUT:

    - ``c`` -- an element of ``self.base_ring()``

    - ``i``, ``j`` -- integers `\geq 0`

    OUTPUT:

    a :class:`QuadraticForm` (by default, otherwise none)

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 4, range(1,11)); Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 2 3 4 ]
        [ * 5 6 7 ]
        [ * * 8 9 ]
        [ * * * 10 ]

        sage: Q.elementary_substitution(c=1, i=0, j=3)
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 2 3 6 ]
        [ * 5 6 9 ]
        [ * * 8 12 ]
        [ * * * 15 ]

    ::

        sage: R = QuadraticForm(ZZ, 4, range(1,11)); R
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 2 3 4 ]
        [ * 5 6 7 ]
        [ * * 8 9 ]
        [ * * * 10 ]

    ::

        sage: M = Matrix(ZZ, 4, 4, [1,0,0,1, 0,1,0,0, 0,0,1,0, 0,0,0,1]); M
        [1 0 0 1]
        [0 1 0 0]
        [0 0 1 0]
        [0 0 0 1]
        sage: R(M)
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 2 3 6 ]
        [ * 5 6 9 ]
        [ * * 8 12 ]
        [ * * * 15 ]
    """
    if not in_place:
        Q = self.parent()(self.base_ring(), self.dim(), self.coefficients())
        Q.elementary_substitution(c, i, j, True)
        return Q

    # Adjust the a_{k,j} coefficients
    ij_old = self[i, j]    # Store this since it's overwritten, but used in the a_{j,j} computation!
    for k in range(self.dim()):
        if (k != i) and (k != j):
            ans = self[j, k] + c * self[i, k]
            self[j, k] = ans
        elif (k == j):
            ans = self[j, k] + c * ij_old + c * c * self[i, i]
            self[j, k] = ans
        else:
            ans = self[j, k] + 2 * c * self[i, k]
            self[j, k] = ans


def add_symmetric(self, c, i, j, in_place=False):
    r"""
    Perform the substitution `x_j \longmapsto x_j + c\cdot x_i`, which has the
    effect (on associated matrices) of symmetrically adding
    `c` times the `j`-th row/column to the `i`-th row/column.

    NOTE: This is meant for compatibility with previous code,
    which implemented a matrix model for this class.  It is used
    in the method :meth:`local_normal_form`.

    INPUT:

    - ``c`` -- an element of ``self.base_ring()``

    - ``i``, ``j`` -- integers `\geq 0`

    OUTPUT:

    a :class:`QuadraticForm` (by default, otherwise none)

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 3, range(1,7)); Q
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 1 2 3 ]
        [ * 4 5 ]
        [ * * 6 ]
        sage: Q.add_symmetric(-1, 1, 0)
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 1 0 3 ]
        [ * 3 2 ]
        [ * * 6 ]
        sage: Q.add_symmetric(-3/2, 2, 0)     # ERROR: -3/2 isn't in the base ring ZZ
        Traceback (most recent call last):
        ...
        RuntimeError: this coefficient cannot be coerced
        to an element of the base ring for the quadratic form

    ::

        sage: Q = QuadraticForm(QQ, 3, range(1,7)); Q
        Quadratic form in 3 variables over Rational Field with coefficients:
        [ 1 2 3 ]
        [ * 4 5 ]
        [ * * 6 ]
        sage: Q.add_symmetric(-3/2, 2, 0)
        Quadratic form in 3 variables over Rational Field with coefficients:
        [ 1 2 0 ]
        [ * 4 2 ]
        [ * * 15/4 ]
    """
    return self.elementary_substitution(c, j, i, in_place)
