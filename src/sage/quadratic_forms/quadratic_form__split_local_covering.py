"""
Split Local Covering
"""
#######################################################################
# Routines that look for a split local covering for a given quadratic #
# form in 4 variables.                                                #
#######################################################################

from copy import deepcopy

from sage.quadratic_forms.extras import extend_to_primitive

import sage.rings.abc
from sage.rings.real_double import RDF
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.arith.misc import GCD


def cholesky_decomposition(self, bit_prec=53):
    r"""
    Give the Cholesky decomposition of this quadratic form `Q` as a real matrix
    of precision ``bit_prec``.

    RESTRICTIONS:

        `Q` must be given as a :class:`QuadraticForm` defined over `\ZZ`, `\QQ`, or some
        real field. If it is over some real field, then an error is raised if
        the precision given is not less than the defined precision of the real
        field defining the quadratic form!

    REFERENCE:

    - Cohen's "A Course in Computational Algebraic Number Theory" book, p 103.

    INPUT:

    - ``bit_prec`` -- a natural number (default: 53)

    OUTPUT: an upper triangular real matrix of precision ``bit_prec``


    .. TODO::

        If we only care about working over the real double field (``RDF``), then we
        can use the method :meth:`cholesky` present for square matrices over that.

    .. NOTE::

        There is a note in the original code reading

        ::

            Finds the Cholesky decomposition of a quadratic form -- as an upper-triangular matrix!
            (It's assumed to be global, hence twice the form it refers to.)  <-- Python revision asks:  Is this true?!? =|

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.cholesky_decomposition()
        [ 1.00000000000000 0.000000000000000 0.000000000000000]
        [0.000000000000000  1.00000000000000 0.000000000000000]
        [0.000000000000000 0.000000000000000  1.00000000000000]

    ::

        sage: Q = QuadraticForm(QQ, 3, range(1,7)); Q
        Quadratic form in 3 variables over Rational Field with coefficients:
        [ 1 2 3 ]
        [ * 4 5 ]
        [ * * 6 ]
        sage: Q.cholesky_decomposition()
        [ 1.00000000000000  1.00000000000000  1.50000000000000]
        [0.000000000000000  3.00000000000000 0.333333333333333]
        [0.000000000000000 0.000000000000000  3.41666666666667]
    """
    # Check that the precision passed is allowed.
    if isinstance(self.base_ring(), sage.rings.abc.RealField) and (self.base_ring().prec() < bit_prec):
        raise RuntimeError("the precision requested is greater than that of the given quadratic form")

    from sage.rings.real_mpfr import RealField

    # 1. Initialization
    n = self.dim()
    R = RealField(bit_prec)
    MS = MatrixSpace(R, n, n)
    Q = MS(R(0.5)) * MS(self.matrix())               # Initialize the real symmetric matrix A with the matrix for Q(x) = x^t * A * x

    # DIAGNOSTIC

    # 2. Loop on i
    for i in range(n):
        for j in range(i + 1, n):
            Q[j, i] = Q[i, j]             # Is this line redundant?
            Q[i, j] = Q[i, j] / Q[i, i]

        # 3. Main Loop
        for k in range(i + 1, n):
            for l in range(k, n):
                Q[k, l] = Q[k, l] - Q[k, i] * Q[i, l]

    # 4. Zero out the strictly lower-triangular entries
    for i in range(n):
        for j in range(i):
            Q[i, j] = 0

    return Q


def vectors_by_length(self, bound):
    r"""
    Return a list of short vectors together with their values.

    This is a naive algorithm which uses the Cholesky decomposition,
    but does not use the LLL-reduction algorithm.

    INPUT:

    - ``bound`` -- integer `\geq 0`

    OUTPUT:

    - a list ``L`` of length (``bound`` + 1) whose entry ``L[i]`` is a list of
      all vectors of length `i`.

    REFERENCES:

    This is a slightly modified version of Cohn's Algorithm
    2.7.5 in "A Course in Computational Number Theory", with the
    increment step moved around and slightly re-indexed to allow clean
    looping.

    .. NOTE::

        We could speed this up for very skew matrices by using LLL
        first, and then changing coordinates back, but for our purposes
        the simpler method is efficient enough.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1])
        sage: Q.vectors_by_length(5)                                                    # needs sage.symbolic
        [[[0, 0]],
         [[0, -1], [-1, 0]],
         [[-1, -1], [1, -1]],
         [],
         [[0, -2], [-2, 0]],
         [[-1, -2], [1, -2], [-2, -1], [2, -1]]]

    ::

        sage: Q1 = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q1.vectors_by_length(5)                                                   # needs sage.symbolic
        [[[0, 0, 0, 0]],
         [[-1, 0, 0, 0]],
         [],
         [[0, -1, 0, 0]],
         [[-1, -1, 0, 0], [1, -1, 0, 0], [-2, 0, 0, 0]],
         [[0, 0, -1, 0]]]

    ::

        sage: Q = QuadraticForm(ZZ, 4, [1,1,1,1, 1,0,0, 1,0, 1])
        sage: list(map(len, Q.vectors_by_length(2)))                                    # needs sage.symbolic
        [1, 12, 12]

    ::

        sage: Q = QuadraticForm(ZZ, 4, [1,-1,-1,-1, 1,0,0, 4,-3, 4])
        sage: list(map(len, Q.vectors_by_length(3)))                                    # needs sage.symbolic
        [1, 3, 0, 3]
    """
    # pari uses eps = 1e-6 ; nothing bad should happen if eps is too big
    # but if eps is too small, roundoff errors may knock off some
    # vectors of norm = bound (see #7100)
    eps = RDF(1e-6)
    bound = ZZ(max(bound, 0))  # bound is an integer
    Theta_Precision = bound + eps
    n = self.dim()

    # Make the vector of vectors which have a given value
    # (So theta_vec[i] will have all vectors v with Q(v) = i.)
    theta_vec = [[] for i in range(bound + 1)]

    # Initialize Q with zeros and Copy the Cholesky array into Q
    Q = self.cholesky_decomposition()

    # 1. Initialize
    T = n * [RDF(0)]    # Note: We index the entries as 0 --> n-1
    U = n * [RDF(0)]
    i = n - 1
    T[i] = RDF(Theta_Precision)
    U[i] = RDF(0)

    L = n * [0]
    x = n * [0]

    # 2. Compute bounds
    Z = (T[i] / Q[i][i]).sqrt(extend=False)
    L[i] = (Z - U[i]).floor()
    x[i] = (-Z - U[i]).ceil()

    done_flag = False
    Q_val = 0  # WARNING: Still need a good way of checking overflow for this value...

    # Big loop which runs through all vectors
    while not done_flag:

        # 3b. Main loop -- try to generate a complete vector x (when i=0)
        while (i > 0):
            T[i - 1] = T[i] - Q[i][i] * (x[i] + U[i]) * (x[i] + U[i])
            i = i - 1
            U[i] = 0
            for j in range(i + 1, n):
                U[i] = U[i] + Q[i][j] * x[j]

            # Now go back and compute the bounds...
            # 2. Compute bounds
            Z = (T[i] / Q[i][i]).sqrt(extend=False)
            L[i] = (Z - U[i]).floor()
            x[i] = (-Z - U[i]).ceil()

            # carry if we go out of bounds -- when Z is so small that
            # there aren't any integral vectors between the bounds
            # Note: this ensures T[i-1] >= 0 in the next iteration
            while (x[i] > L[i]):
                i += 1
                x[i] += 1

        # 4. Solution found (This happens when i = 0)
        Q_val_double = Theta_Precision - T[0] + Q[0][0] * (x[0] + U[0]) * (x[0] + U[0])
        Q_val = Q_val_double.round()

        # SANITY CHECK: Roundoff Error is < 0.001
        if abs(Q_val_double - Q_val) > 0.001:
            print(" x = ", x)
            print(" Float = ", Q_val_double, "   Long = ", Q_val)
            raise RuntimeError("The roundoff error is bigger than 0.001, so we should use more precision somewhere...")

        if Q_val <= bound:
            theta_vec[Q_val].append(deepcopy(x))

        # 5. Check if x = 0, for exit condition. =)
        j = 0
        done_flag = True
        while j < n:
            if x[j] != 0:
                done_flag = False
            j += 1

        # 3a. Increment (and carry if we go out of bounds)
        x[i] += 1
        while i < n - 1 and x[i] > L[i]:
            i += 1
            x[i] += 1

    return theta_vec


def complementary_subform_to_vector(self, v):
    r"""
    Find the `(n-1)`-dimensional quadratic form orthogonal to the vector `v`.

    .. NOTE::

        This is usually not a direct summand!

    .. NOTE::

        There is a minor difference in the cancellation
        code here (form the C++ version) since the notation ``Q[i,j]`` indexes
        coefficients of the quadratic polynomial here, not the symmetric
        matrix.  Also, it produces a better splitting now, for the full
        lattice (as opposed to a sublattice in the C++ code) since we
        now extend `v` to a unimodular matrix.

    INPUT:

    - ``v`` -- list of ``self.dim()`` integers

    OUTPUT: a :class:`QuadraticForm` over `\ZZ`

    EXAMPLES::

        sage: Q1 = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q1.complementary_subform_to_vector([1,0,0,0])
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 7 0 0 ]
        [ * 5 0 ]
        [ * * 3 ]

    ::

        sage: Q1.complementary_subform_to_vector([1,1,0,0])
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 7 0 0 ]
        [ * 5 0 ]
        [ * * 12 ]

    ::

        sage: Q1.complementary_subform_to_vector([1,1,1,1])
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 880 -480 -160 ]
        [ * 624 -96 ]
        [ * * 240 ]
    """
    n = self.dim()

    # Copy the quadratic form
    Q = deepcopy(self)

    # Find the first nonzero component of v, and call it nz  (Note: 0 <= nz < n)
    nz = 0
    while nz < n and v[nz] == 0:
        nz += 1

    # Abort if v is the zero vector
    if nz == n:
        raise TypeError("v cannot be the zero vector")

    # Make the change of basis matrix
    new_basis = extend_to_primitive(matrix(ZZ, n, 1, v))

    # Change Q (to Q1) to have v as its nz-th basis vector
    Q1 = Q(new_basis)

    # Pick out the value Q(v) of the vector
    d = Q1[0, 0]

    # For each row/column, perform elementary operations to cancel them out.
    for i in range(1, n):

        # Check if the (i,0)-entry is divisible by d,
        # and stretch its row/column if not.
        if Q1[i, 0] % d:
            Q1 = Q1.multiply_variable(d / GCD(d, Q1[i, 0] // 2), i)

        # Now perform the (symmetric) elementary operations to cancel out the (i,0) entries/
        Q1 = Q1.add_symmetric(-(Q1[i, 0] // 2) / (GCD(d, Q1[i, 0] // 2)), i, 0)

    # Check that we're done!
    done_flag = True
    for i in range(1, n):
        if Q1[0, i] != 0:
            done_flag = False

    if not done_flag:
        raise RuntimeError("There is a problem cancelling out the matrix entries! =O")

    # Return the complementary matrix
    return Q1.extract_variables(range(1, n))


def split_local_cover(self):
    r"""
    Try to find subform of the given (positive definite quaternary)
    quadratic form `Q` of the form

    .. MATH::

        d\cdot x^2 + T(y,z,w)

    where `d > 0` is as small as possible.

    This is done by exhaustive search on small vectors, and then
    comparing the local conditions of its sum with its complementary
    lattice and the original quadratic form `Q`.

    OUTPUT: a :class:`QuadraticForm` over `\ZZ`

    EXAMPLES::

        sage: Q1 = DiagonalQuadraticForm(ZZ, [7,5,3])
        sage: Q1.split_local_cover()                                                    # needs sage.symbolic
        Quadratic form in 3 variables over Integer Ring with coefficients:
        [ 3 0 0 ]
        [ * 5 0 ]
        [ * * 7 ]
    """
    from sage.quadratic_forms.quadratic_form import QuadraticForm

    # 0. If a split local cover already exists, then return it.
    if hasattr(self, "__split_local_cover"):
        if isinstance(self.__split_local_cover, QuadraticForm):  # Here the computation has been done.
            return self.__split_local_cover
        elif self.__split_local_cover in ZZ:    # Here it indexes the values already tried!
            current_length = self.__split_local_cover + 1
            Length_Max = current_length + 5
    else:
        current_length = 1
        Length_Max = 6

    # 1. Find a range of new vectors
    all_vectors = self.vectors_by_length(Length_Max)
    current_vectors = all_vectors[current_length]

    # Loop until we find a split local cover...
    while True:

        # 2. Check if any of the primitive ones produce a split local cover
        for v in current_vectors:
            Q = QuadraticForm(ZZ, 1, [current_length]) + self.complementary_subform_to_vector(v)
            if Q.local_representation_conditions() == self.local_representation_conditions():
                self.__split_local_cover = Q
                return Q

        # 3. Save what we have checked and get more vectors.
        self.__split_local_cover = current_length
        current_length += 1
        if current_length >= len(all_vectors):
            Length_Max += 5
            all_vectors = self.vectors_by_length(Length_Max)
        current_vectors = all_vectors[current_length]
