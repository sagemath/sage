# sage_setup: distribution = sagemath-pari
# sage.doctest: needs sage.libs.pari sage.modules
"""
Local Normal Form
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

import copy

from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.arith.misc import GCD, valuation, is_prime


def find_entry_with_minimal_scale_at_prime(self, p):
    r"""
    Find the entry of the quadratic form with minimal scale at the
    prime `p`, preferring diagonal entries in case of a tie.

    (I.e.  If
    we write the quadratic form as a symmetric matrix `M`, then this
    entry ``M[i,j]`` has the minimal valuation at the prime `p`.)

    .. NOTE::

        This answer is independent of the kind of matrix (Gram or
        Hessian) associated to the form.

    INPUT:

    - ``p`` -- a prime number > 0

    OUTPUT: a pair of integers `\geq 0`

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [6, 2, 20]); Q
        Quadratic form in 2 variables over Integer Ring with coefficients:
          [ 6 2 ]
          [ * 20 ]
        sage: Q.find_entry_with_minimal_scale_at_prime(2)
        (0, 1)
        sage: Q.find_entry_with_minimal_scale_at_prime(3)
        (1, 1)
        sage: Q.find_entry_with_minimal_scale_at_prime(5)
        (0, 0)
    """
    n = self.dim()
    min_val = Infinity
    ij_index = None
    val_2 = valuation(2, p)
    for d in range(n):           # d = difference j-i
        for e in range(n - d):   # e is the length of the diagonal with value d.

            # Compute the valuation of the entry
            if d == 0:
                tmp_val = valuation(self[e, e + d], p)
            else:
                tmp_val = valuation(self[e, e + d], p) - val_2

            # Check if it's any smaller than what we have
            if tmp_val < min_val:
                ij_index = (e, e + d)
                min_val = tmp_val

    # Return the result
    return ij_index


def local_normal_form(self, p):
    r"""
    Return a locally integrally equivalent quadratic form over
    the `p`-adic integers `\ZZ_p` which gives the Jordan decomposition.

    The Jordan components are written as sums of blocks of size `\leq 2`
    and are arranged by increasing scale, and then by increasing norm.
    This is equivalent to saying that we put the `1 \times 1` blocks before
    the `2 \times 2` blocks in each Jordan component.

    INPUT:

    - ``p`` -- a positive prime number

    OUTPUT: a quadratic form over `\ZZ`

    .. WARNING::

        Currently this only works for quadratic forms defined over `\ZZ`.

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [10,4,1])
        sage: Q.local_normal_form(5)
        Quadratic form in 2 variables over Integer Ring with coefficients:
          [ 1 0 ]
          [ * 6 ]

    ::

        sage: Q.local_normal_form(3)
        Quadratic form in 2 variables over Integer Ring with coefficients:
          [ 10 0 ]
          [ * 15 ]

        sage: Q.local_normal_form(2)
        Quadratic form in 2 variables over Integer Ring with coefficients:
          [ 1 0 ]
          [ * 6 ]
    """
    # Sanity Checks
    if self.base_ring() != ZZ:
        raise NotImplementedError("this currently only works for quadratic forms defined over ZZ")
    if not (p >= 2 and is_prime(p)):
        raise TypeError("p is not a positive prime number")

    # Some useful local variables
    Q = self.parent()(self.base_ring(), self.dim(), self.coefficients())

    # Prepare the final form to return
    Q_Jordan = copy.deepcopy(self)
    Q_Jordan.__init__(self.base_ring(), 0)

    while Q.dim() > 0:
        n = Q.dim()

        # Step 1: Find the minimally p-divisible matrix entry, preferring diagonals
        # -------------------------------------------------------------------------
        (min_i, min_j) = Q.find_entry_with_minimal_scale_at_prime(p)
        if min_i == min_j:
            min_val = valuation(2 * Q[min_i, min_j], p)
        else:
            min_val = valuation(Q[min_i, min_j], p)

        # Error if we still haven't seen nonzero coefficients!
        if min_val == Infinity:
            raise RuntimeError("the original matrix is degenerate")

        # Step 2: Arrange for the upper leftmost entry to have minimal valuation
        # ----------------------------------------------------------------------
        if min_i == min_j:
            block_size = 1
            Q.swap_variables(0, min_i, in_place=True)
        else:
            # Work in the upper-left 2x2 block, and replace it by its 2-adic equivalent form
            Q.swap_variables(0, min_i, in_place=True)
            Q.swap_variables(1, min_j, in_place=True)

            # 1x1 => make upper left the smallest
            if p != 2:
                block_size = 1
                Q.add_symmetric(1, 0, 1, in_place=True)
            # 2x2 => replace it with the appropriate 2x2 matrix
            else:
                block_size = 2

        # Step 3: Clear out the remaining entries
        #  ---------------------------------------
        min_scale = p ** min_val                             # This is the minimal valuation of the Hessian matrix entries.

        # Perform cancellation over Z by ensuring divisibility
        if block_size == 1:
            a = 2 * Q[0, 0]
            for j in range(block_size, n):
                b = Q[0, j]
                g = GCD(a, b)

                # Sanity Check:  a/g is a p-unit
                if valuation(g, p) != valuation(a, p):
                    raise RuntimeError("we have a problem with our rescaling not preserving p-integrality")

                Q.multiply_variable(ZZ(a / g), j, in_place=True)   # Ensures that the new b entry is divisible by a
                Q.add_symmetric(ZZ(-b / g), j, 0, in_place=True)  # Performs the cancellation

        elif block_size == 2:
            a1 = 2 * Q[0, 0]
            a2 = Q[0, 1]
            b1 = Q[1, 0]      # This is the same as a2
            b2 = 2 * Q[1, 1]

            big_det = a1 * b2 - a2 * b1
            small_det = big_det / (min_scale * min_scale)

            # Cancels out the rows/columns of the 2x2 block
            for j in range(block_size, n):
                a = Q[0, j]
                b = Q[1, j]

                # Ensures an integral result (scale jth row/column by big_det)
                Q.multiply_variable(big_det, j, in_place=True)

                # Performs the cancellation (by producing -big_det * jth row/column)
                Q.add_symmetric(ZZ(-(a * b2 - b * a2)), j, 0, in_place=True)
                Q.add_symmetric(ZZ(-(-a * b1 + b * a1)), j, 1, in_place=True)

                # Now remove the extra factor (non p-unit factor) in big_det we introduced above
                Q.divide_variable(ZZ(min_scale * min_scale), j, in_place=True)

            # Uses Cassels's proof to replace the remaining 2 x 2 block
            if (1 + small_det) % 8 == 0:
                Q[0, 0] = 0
                Q[1, 1] = 0
                Q[0, 1] = min_scale
            elif (5 + small_det) % 8 == 0:
                Q[0, 0] = min_scale
                Q[1, 1] = min_scale
                Q[0, 1] = min_scale
            else:
                raise RuntimeError("Error in LocalNormal: Impossible behavior for a 2x2 block! \n")

        # Check that the cancellation worked, extract the upper-left block, and trim Q to handle the next block.
        for i in range(block_size):
            for j in range(block_size, n):
                if Q[i, j] != 0:
                    raise RuntimeError(f"the cancellation did not work properly at entry ({i},{j})")
        Q_Jordan = Q_Jordan + Q.extract_variables(range(block_size))
        Q = Q.extract_variables(range(block_size, n))

    return Q_Jordan


def jordan_blocks_by_scale_and_unimodular(self, p, safe_flag=True):
    r"""
    Return a list of pairs `(s_i, L_i)` where `L_i` is a maximal
    `p^{s_i}`-unimodular Jordan component which is further decomposed into
    block diagonals of block size `\le 2`.

    For each `L_i` the `2 \times 2` blocks are listed after the `1 \times 1` blocks
    (which follows from the convention of the
    :meth:`local_normal_form` method).

    .. NOTE::

        The decomposition of each `L_i` into smaller blocks is not unique!

    The ``safe_flag`` argument allows us to select whether we want a copy of
    the output, or the original output.  By default ``safe_flag = True``, so we
    return a copy of the cached information.  If this is set to ``False``, then
    the routine is much faster but the return values are vulnerable to being
    corrupted by the user.

    INPUT:

    - ``p`` -- a prime number > 0

    OUTPUT:

    A list of pairs `(s_i, L_i)` where:

    - `s_i` is an integer,
    - `L_i` is a block-diagonal unimodular quadratic form over `\ZZ_p`.

    .. NOTE::

        These forms `L_i` are defined over the `p`-adic integers, but by a
        matrix over `\ZZ` (or `\QQ`?).

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,9,5,7])
        sage: Q.jordan_blocks_by_scale_and_unimodular(3)
        [(0, Quadratic form in 3 variables over Integer Ring with coefficients:
               [ 1 0 0 ]
               [ * 5 0 ]
               [ * * 7 ]),
         (2, Quadratic form in 1 variables over Integer Ring with coefficients:
               [ 1 ])]

    ::

        sage: Q2 = QuadraticForm(ZZ, 2, [1,1,1])
        sage: Q2.jordan_blocks_by_scale_and_unimodular(2)
        [(-1, Quadratic form in 2 variables over Integer Ring with coefficients:
                [ 2 2 ]
                [ * 2 ])]
        sage: Q = Q2 + Q2.scale_by_factor(2)
        sage: Q.jordan_blocks_by_scale_and_unimodular(2)
        [(-1, Quadratic form in 2 variables over Integer Ring with coefficients:
                [ 2 2 ]
                [ * 2 ]),
         (0, Quadratic form in 2 variables over Integer Ring with coefficients:
                [ 2 2 ]
                [ * 2 ])]
    """
    # Try to use the cached result
    try:
        if safe_flag:
            return copy.deepcopy(self.__jordan_blocks_by_scale_and_unimodular_dict[p])
        else:
            return self.__jordan_blocks_by_scale_and_unimodular_dict[p]
    except (KeyError, AttributeError):
        # Initialize the global dictionary if it doesn't exist
        if not hasattr(self, '__jordan_blocks_by_scale_and_unimodular_dict'):
            self.__jordan_blocks_by_scale_and_unimodular_dict = {}

    # Deal with zero dim'l forms
    if self.dim() == 0:
        return []

    # Find the Local Normal form of Q at p
    Q1 = self.local_normal_form(p)

    # Parse this into Jordan Blocks
    n = Q1.dim()
    tmp_Jordan_list = []
    i = 0
    start_ind = 0
    if n >= 2 and Q1[0, 1] != 0:
        start_scale = valuation(Q1[0, 1], p) - 1
    else:
        start_scale = valuation(Q1[0, 0], p)

    while i < n:

        # Determine the size of the current block
        if i == n - 1 or Q1[i, i + 1] == 0:
            block_size = 1
        else:
            block_size = 2

        # Determine the valuation of the current block
        if block_size == 1:
            block_scale = valuation(Q1[i, i], p)
        else:
            block_scale = valuation(Q1[i, i + 1], p) - 1

        # Process the previous block if the valuation increased
        if block_scale > start_scale:
            tmp_Jordan_list += [(start_scale, Q1.extract_variables(range(start_ind, i)).scale_by_factor(ZZ.one() / QQ(p)**start_scale))]
            start_ind = i
            start_scale = block_scale

        # Increment the index
        i += block_size

    # Add the last block
    tmp_Jordan_list += [(start_scale, Q1.extract_variables(range(start_ind, n)).scale_by_factor(ZZ.one() / QQ(p)**start_scale))]

    # Cache the result
    self.__jordan_blocks_by_scale_and_unimodular_dict[p] = tmp_Jordan_list

    # Return the result
    return tmp_Jordan_list


def jordan_blocks_in_unimodular_list_by_scale_power(self, p):
    r"""
    Return a list of Jordan components, whose component at index `i`
    should be scaled by the factor `p^i`.

    This is only defined for integer-valued quadratic forms
    (i.e., forms with base ring `\ZZ`), and the indexing only works
    correctly for `p=2` when the form has an integer Gram matrix.

    INPUT:

    - ``self`` -- a quadratic form over `\ZZ`, which has integer Gram matrix if `p = 2`
    - ``p`` -- a prime number > 0

    OUTPUT: list of `p`-unimodular quadratic forms

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 3, [2, -2, 0, 3, -5, 4])
        sage: Q.jordan_blocks_in_unimodular_list_by_scale_power(2)
        Traceback (most recent call last):
        ...
        TypeError: the given quadratic form has a Jordan component with a negative scale exponent

        sage: Q.scale_by_factor(2).jordan_blocks_in_unimodular_list_by_scale_power(2)
        [Quadratic form in 2 variables over Integer Ring with coefficients:
           [ 0 2 ]
           [ * 0 ],
         Quadratic form in 0 variables over Integer Ring with coefficients:
           ,
         Quadratic form in 1 variables over Integer Ring with coefficients:
           [ 345 ]]

        sage: Q.jordan_blocks_in_unimodular_list_by_scale_power(3)
        [Quadratic form in 2 variables over Integer Ring with coefficients:
           [ 2 0 ]
           [ * 10 ],
         Quadratic form in 1 variables over Integer Ring with coefficients:
           [ 2 ]]
    """
    # Sanity Check
    if self.base_ring() != ZZ:
        raise TypeError("this method only makes sense for integer-valued quadratic forms (i.e. defined over ZZ)")

    # Deal with zero dim'l forms
    if self.dim() == 0:
        return []

    # Find the Jordan Decomposition
    list_of_jordan_pairs = self.jordan_blocks_by_scale_and_unimodular(p)
    scale_list = [P[0] for P in list_of_jordan_pairs]
    s_max = max(scale_list)
    if min(scale_list) < 0:
        raise TypeError("the given quadratic form has a Jordan component with a negative scale exponent")

    # Make the new list of unimodular Jordan components
    zero_form = self.parent()(ZZ, 0)
    list_by_scale = [zero_form for _ in range(s_max + 1)]
    for P in list_of_jordan_pairs:
        list_by_scale[P[0]] = P[1]

    # Return the new list
    return list_by_scale
