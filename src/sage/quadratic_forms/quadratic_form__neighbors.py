"""
Neighbors
"""
from __future__ import print_function

from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from copy import deepcopy
from sage.quadratic_forms.extras import extend_to_primitive
from sage.matrix.constructor import matrix


# ############################################################################
# Routines used for understanding p-neighbors and computing classes in a genus
# ############################################################################


def find_primitive_p_divisible_vector__random(self, p):
    """
    Find a random `p`-primitive vector in `L/pL` whose value is `p`-divisible.

    .. note::

        Since there are about `p^{(n-2)}` of these lines, we have a `1/p`
        chance of randomly finding an appropriate vector.

    .. warning::

        If there are local obstructions for this to happen, then this algorithm
        will never terminate... =(  We should check for this too!

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [10,1,4])
        sage: Q.find_primitive_p_divisible_vector__random(5)    # random
        (1, 1)
        sage: Q.find_primitive_p_divisible_vector__random(5)    # random
        (1, 0)
        sage: Q.find_primitive_p_divisible_vector__random(5)    # random
        (2, 0)
        sage: Q.find_primitive_p_divisible_vector__random(5)    # random
        (2, 2)
        sage: Q.find_primitive_p_divisible_vector__random(5)    # random
        (3, 3)
        sage: Q.find_primitive_p_divisible_vector__random(5)    # random
        (3, 3)
        sage: Q.find_primitive_p_divisible_vector__random(5)    # random
        (2, 0)
    """
    n = self.dim()
    v = vector([ZZ.random_element(p) for _ in range(n)])

    # Repeatedly choose random vectors, and evaluate until the value
    # is p-divisible.
    while True:
        if (self(v) % p == 0) and (v != 0):
            return v
        else:
            v[ZZ.random_element(n)] = ZZ.random_element(p)
            # Replace a random entry and try again.


def find_primitive_p_divisible_vector__next(self, p, v=None):
    """
    Find the next `p`-primitive vector (up to scaling) in `L/pL` whose
    value is `p`-divisible, where the last vector returned was `v`.  For
    an initial call, no `v` needs to be passed.

    Returns vectors whose last non-zero entry is normalized to 0 or 1 (so no
    lines are counted repeatedly).  The ordering is by increasing the
    first non-normalized entry.  If we have tested all (lines of)
    vectors, then return None.

    OUTPUT:

    vector or None

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 2, [10,1,4])
        sage: v = Q.find_primitive_p_divisible_vector__next(5); v
        (1, 1)
        sage: v = Q.find_primitive_p_divisible_vector__next(5, v); v
        (1, 0)
        sage: v = Q.find_primitive_p_divisible_vector__next(5, v); v
    """
    # Initialize
    n = self.dim()
    if v is None:
        w = vector(ZZ, [0] * (n - 1) + [1])
    else:
        w = deepcopy(v)

    # Handle n = 1 separately.
    if n <= 1:
        raise NotImplementedError("Sorry -- Not implemented yet!")

    # Look for the last non-zero entry (which must be 1)
    nz = n - 1
    while w[nz] == 0:
        nz += -1

    # Test that the last non-zero entry is 1 (to detect tampering).
    if w[nz] != 1:
        print("Warning: The input vector to QuadraticForm.find_primitive_p_divisible_vector__next() is not normalized properly.")

    # Look for the next vector, until w == 0
    while True:

        # Look for the first non-maximal (non-normalized) entry
        ind = 0
        while (ind < nz) and (w[ind] == p - 1):
            ind += 1

        # Increment
        if ind < nz:
            w[ind] += 1
            for j in range(ind):
                w[j] = 0
        else:
            for j in range(ind + 1):    # Clear all entries
                w[j] = 0

            if nz != 0:
                # Move the non-zero normalized index over by one, or
                # return the zero vector
                w[nz - 1] = 1
                nz += -1

        # Test for zero vector
        if w == 0:
            return None

        # Test for p-divisibility
        if (self(w) % p == 0):
            return w

# -------------------------------------------------------------------------


def find_p_neighbor_from_vec(self, p, v):
    r"""
    Find the `p`-neighbor of this quadratic form associated to a given
    vector `v` satisfying:

    #. `Q(v) = 0  \pmod p`
    #. `v` is a non-singular point of the conic `Q(v) = 0 \pmod p`.

    Reference:  Gonzalo Tornaria's Thesis, Theorem 3.5, p34.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ,[1,1,1,1])
        sage: v = vector([0,2,1,1])
        sage: X = Q.find_p_neighbor_from_vec(3,v); X
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 3 10 0 -4 ]
        [ * 9 0 -6 ]
        [ * * 1 0 ]
        [ * * * 2 ]
    """
    R = self.base_ring()
    n = self.dim()
    B2 = self.matrix()

    # Find a (dual) vector w with B(v,w) != 0 (mod p)
    # We want the dot product with this to not be divisible by 2*p.
    v_dual = B2 * vector(v)

    y_ind = 0
    while ((y_ind < n) and (v_dual[y_ind] % p) == 0):
        # Check the dot product for the std basis vectors!
        y_ind += 1
    if y_ind == n:
        raise RuntimeError("Oops!  One of the standard basis vectors "
                           "should have worked.")
    w = vector(R, [R(i == y_ind) for i in range(n)])
    vw_prod = (v * self.matrix()).dot_product(w)

    # Lift the vector v to a vector v1 s.t. Q(v1) = 0 (mod p^2)
    s = self(v)
    if (s % p**2 != 0):
        al = (-s / (p * vw_prod)) % p
        v1 = v + p * al * w
        v1w_prod = (v1 * self.matrix()).dot_product(w)
    else:
        v1 = v
        v1w_prod = vw_prod

    # Construct a special p-divisible basis to use for the p-neighbor switch
    good_basis = extend_to_primitive([v1, w])
    for i in range(2, n):
        ith_prod = (good_basis[i] * self.matrix()).dot_product(v)
        c = (ith_prod / v1w_prod) % p
        good_basis[i] = good_basis[i] - c * w
        # Ensures that this extension has <v_i, v> = 0 (mod p)

    # Perform the p-neighbor switch
    good_basis[0] = vector([x / p for x in good_basis[0]])   # Divide v1 by p
    good_basis[1] = good_basis[1] * p                        # Multiply w by p

    # Return the associated quadratic form
    M = matrix(good_basis)
    QF = self.parent()
    return QF(R, M * self.matrix() * M.transpose())
