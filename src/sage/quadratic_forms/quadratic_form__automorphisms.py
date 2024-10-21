# sage_setup: distribution = sagemath-pari
# sage.doctest: needs sage.libs.pari sage.modules
"""
Automorphisms of Quadratic Forms
"""
# ****************************************************************************
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import Matrix
from sage.rings.integer_ring import ZZ

from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import vector
from sage.arith.misc import GCD


@cached_method
def basis_of_short_vectors(self, show_lengths=False):
    r"""
    Return a basis for `\ZZ^n` made of vectors with minimal lengths `Q(v)`.

    OUTPUT: a tuple of vectors, and optionally a tuple of values for each vector

    This uses :pari:`qfminim`.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.basis_of_short_vectors()
        ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1))
        sage: Q.basis_of_short_vectors(True)
        (((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)), (1, 3, 5, 7))

    The returned vectors are immutable::

        sage: v = Q.basis_of_short_vectors()[0]
        sage: v
        (1, 0, 0, 0)
        sage: v[0] = 0
        Traceback (most recent call last):
        ...
        ValueError: vector is immutable; please change a copy instead (use copy())
    """
    # Set an upper bound for the number of vectors to consider
    Max_number_of_vectors = 10000

    # Generate a PARI matrix for the associated Hessian matrix
    M_pari = self.__pari__()

    # Run through all possible minimal lengths to find a spanning set of vectors
    n = self.dim()
    M1 = Matrix([[0]])
    vec_len = 0
    while M1.rank() < n:
        vec_len += 1
        pari_mat = M_pari.qfminim(vec_len, Max_number_of_vectors)[2]
        number_of_vecs = ZZ(pari_mat.matsize()[1])
        vector_list = []
        for i in range(number_of_vecs):
            new_vec = vector([ZZ(x) for x in list(pari_mat[i])])
            vector_list.append(new_vec)

        # Make a matrix from the short vectors
        if vector_list:
            M1 = Matrix(vector_list)

    # Organize these vectors by length (and also introduce their negatives)
    max_len = vec_len // 2
    vector_list_by_length = [[] for _ in range(max_len + 1)]
    for v in vector_list:
        l = self(v)
        vector_list_by_length[l].append(v)
        vector_list_by_length[l].append(vector([-x for x in v]))

    # Make a matrix from the column vectors (in order of ascending length).
    sorted_list = [v for i in range(len(vector_list_by_length))
                   for v in vector_list_by_length[i]]
    sorted_matrix = Matrix(sorted_list).transpose()

    # Determine a basis of vectors of minimal length
    pivots = sorted_matrix.pivots()
    basis = tuple(sorted_matrix.column(i) for i in pivots)
    for v in basis:
        v.set_immutable()

    # Return the appropriate result
    if show_lengths:
        pivot_lengths = tuple(self(v) for v in basis)
        return basis, pivot_lengths
    else:
        return basis


def short_vector_list_up_to_length(self, len_bound, up_to_sign_flag=False):
    """
    Return a list of lists of short vectors `v`, sorted by length, with
    `Q(v) <` ``len_bound``.

    INPUT:

    - ``len_bound`` -- bound for the length of the vectors

    - ``up_to_sign_flag`` -- boolean (default: ``False``); if set to ``True``,
      then only one of the vectors of the pair `[v, -v]` is listed

    OUTPUT:

    A list of lists of vectors such that entry `[i]` contains all
    vectors of length `i`.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.short_vector_list_up_to_length(3)
        [[(0, 0, 0, 0)], [(1, 0, 0, 0), (-1, 0, 0, 0)], []]
        sage: Q.short_vector_list_up_to_length(4)
        [[(0, 0, 0, 0)],
         [(1, 0, 0, 0), (-1, 0, 0, 0)],
         [],
         [(0, 1, 0, 0), (0, -1, 0, 0)]]
        sage: Q.short_vector_list_up_to_length(5)
        [[(0, 0, 0, 0)],
         [(1, 0, 0, 0), (-1, 0, 0, 0)],
         [],
         [(0, 1, 0, 0), (0, -1, 0, 0)],
         [(1, 1, 0, 0),
          (-1, -1, 0, 0),
          (1, -1, 0, 0),
          (-1, 1, 0, 0),
          (2, 0, 0, 0),
          (-2, 0, 0, 0)]]
        sage: Q.short_vector_list_up_to_length(5, True)
        [[(0, 0, 0, 0)],
         [(1, 0, 0, 0)],
         [],
         [(0, 1, 0, 0)],
         [(1, 1, 0, 0), (1, -1, 0, 0), (2, 0, 0, 0)]]
        sage: m6 = matrix(6, [2, 1, 1, 1, -1, -1, 1, 2, 1, 1, -1, -1,
        ....:                 1, 1, 2, 0, -1, -1, 1, 1, 0, 2, 0, -1,
        ....:                 -1, -1, -1, 0, 2, 1, -1, -1, -1, -1, 1, 2])
        sage: Q = QuadraticForm(m6)
        sage: vs = Q.short_vector_list_up_to_length(8)
        sage: [len(vs[i]) for i in range(len(vs))]
        [1, 72, 270, 720, 936, 2160, 2214, 3600]

    The cases of ``len_bound < 2`` led to exception or infinite runtime before.

    ::

        sage: Q.short_vector_list_up_to_length(-1)
        []
        sage: Q.short_vector_list_up_to_length(0)
        []
        sage: Q.short_vector_list_up_to_length(1)
        [[(0, 0, 0, 0, 0, 0)]]

    In the case of quadratic forms that are not positive definite an error is raised.

    ::

        sage: QuadraticForm(matrix(2, [2, 0, 0, -2])).short_vector_list_up_to_length(3)
        Traceback (most recent call last):
        ...
        ValueError: Quadratic form must be positive definite in order to enumerate short vectors

    Check that PARI does not return vectors which are too long::

        sage: Q = QuadraticForm(matrix(2, [72, 12, 12, 120]))
        sage: len_bound_pari = 2*22953421 - 2; len_bound_pari
        45906840
        sage: vs = list(Q.__pari__().qfminim(len_bound_pari)[2])  # long time (18s on sage.math, 2014)
        sage: v = vs[0]; v  # long time
        [66, -623]~
        sage: v.Vec() * Q.__pari__() * v  # long time
        45902280
    """
    if not self.is_positive_definite():
        raise ValueError("Quadratic form must be positive definite "
                         "in order to enumerate short vectors")

    from sage.libs.pari.all import pari

    if len_bound <= 0:
        return []

    # Free module in which the vectors live
    V = FreeModule(ZZ, self.dim())

    # Adjust length for PARI. We need to subtract 1 because PARI returns
    # returns vectors of length less than or equal to b, but we want
    # strictly less. We need to double because the matrix is doubled.
    len_bound_pari = 2 * (len_bound - 1)

    # Call PARI's qfminim()
    parilist = self.__pari__().qfminim(len_bound_pari)[2].Vec()

    # List of lengths
    parilens = pari(r"(M,v) -> vector(#v, i, (v[i]~ * M * v[i])\2)")(self, parilist)

    # Sort the vectors into lists by their length
    vec_sorted_list = [[] for i in range(len_bound)]
    for i in range(len(parilist)):
        length = int(parilens[i])
        # In certain trivial cases, PARI can sometimes return longer
        # vectors than requested.
        if length < len_bound:
            sagevec = V(list(parilist[i]))
            vec_sorted_list[length].append(sagevec)
            if not up_to_sign_flag:
                vec_sorted_list[length].append(-sagevec)

    # Add the zero vector by hand
    vec_sorted_list[0].append(V.zero_vector())

    return vec_sorted_list


def short_primitive_vector_list_up_to_length(self, len_bound, up_to_sign_flag=False):
    r"""
    Return a list of lists of short primitive vectors `v`, sorted by length, with
    `Q(v) <` ``len_bound``.  The list in output `[i]` indexes all vectors of
    length `i`.  If the ``up_to_sign_flag`` is set to ``True``, then only one of
    the vectors of the pair `[v, -v]` is listed.

    .. NOTE::

        This processes the PARI/GP output to always give elements of type `\ZZ`.

    OUTPUT: list of lists of vectors

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.short_vector_list_up_to_length(5, True)
        [[(0, 0, 0, 0)],
         [(1, 0, 0, 0)],
         [],
         [(0, 1, 0, 0)],
         [(1, 1, 0, 0), (1, -1, 0, 0), (2, 0, 0, 0)]]
        sage: Q.short_primitive_vector_list_up_to_length(5, True)
        [[], [(1, 0, 0, 0)], [], [(0, 1, 0, 0)], [(1, 1, 0, 0), (1, -1, 0, 0)]]
    """
    # Get a list of short vectors
    full_vec_list = self.short_vector_list_up_to_length(len_bound, up_to_sign_flag)

    # Make a new list of the primitive vectors
    prim_vec_list = [[v for v in L if GCD(v) == 1]
                     for L in full_vec_list]

    # Return the list of primitive vectors
    return prim_vec_list


def _compute_automorphisms(self):
    """
    Call PARI to compute the automorphism group of the quadratic form.

    This uses :pari:`qfauto`.

    OUTPUT: none, this just caches the result

    TESTS::

        sage: DiagonalQuadraticForm(ZZ, [-1,1,1])._compute_automorphisms()
        Traceback (most recent call last):
        ...
        ValueError: not a definite form in QuadraticForm.automorphisms()
        sage: DiagonalQuadraticForm(GF(5), [1,1,1])._compute_automorphisms()
        Traceback (most recent call last):
        ...
        NotImplementedError: computing the automorphism group of a quadratic form is only supported over ZZ
    """
    if self.base_ring() is not ZZ:
        raise NotImplementedError("computing the automorphism group of a quadratic form is only supported over ZZ")
    if not self.is_definite():
        raise ValueError("not a definite form in QuadraticForm.automorphisms()")

    if hasattr(self, "__automorphisms_pari"):
        return

    A = self.__pari__().qfauto()
    self.__number_of_automorphisms = A[0]
    self.__automorphisms_pari = A[1]


def automorphism_group(self):
    """
    Return the group of automorphisms of the quadratic form.

    OUTPUT: a :class:`MatrixGroup`

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.automorphism_group()
        Matrix group over Rational Field with 3 generators (
        [ 0  0  1]  [1 0 0]  [ 1  0  0]
        [-1  0  0]  [0 0 1]  [ 0 -1  0]
        [ 0  1  0], [0 1 0], [ 0  0  1]
        )

    ::

        sage: DiagonalQuadraticForm(ZZ, [1,3,5,7]).automorphism_group()
        Matrix group over Rational Field with 4 generators (
        [-1  0  0  0]  [ 1  0  0  0]  [ 1  0  0  0]  [ 1  0  0  0]
        [ 0 -1  0  0]  [ 0 -1  0  0]  [ 0  1  0  0]  [ 0  1  0  0]
        [ 0  0 -1  0]  [ 0  0  1  0]  [ 0  0 -1  0]  [ 0  0  1  0]
        [ 0  0  0 -1], [ 0  0  0  1], [ 0  0  0  1], [ 0  0  0 -1]
        )

    The smallest possible automorphism group has order two, since we
    can always change all signs::

        sage: Q = QuadraticForm(ZZ, 3, [2, 1, 2, 2, 1, 3])
        sage: Q.automorphism_group()
        Matrix group over Rational Field with 1 generators (
        [-1  0  0]
        [ 0 -1  0]
        [ 0  0 -1]
        )
    """
    self._compute_automorphisms()

    from sage.matrix.matrix_space import MatrixSpace
    from sage.groups.matrix_gps.finitely_generated import MatrixGroup
    MS = MatrixSpace(self.base_ring().fraction_field(), self.dim(), self.dim())
    gens = [MS(x.sage()) for x in self.__automorphisms_pari]
    return MatrixGroup(gens)


def automorphisms(self):
    """
    Return the list of the automorphisms of the quadratic form.

    OUTPUT: list of matrices

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.number_of_automorphisms()
        48
        sage: 2^3 * factorial(3)
        48
        sage: len(Q.automorphisms())                                                    # needs sage.libs.gap
        48

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.number_of_automorphisms()
        16
        sage: aut = Q.automorphisms()                                                   # needs sage.libs.gap
        sage: len(aut)                                                                  # needs sage.libs.gap
        16
        sage: all(Q(M) == Q for M in aut)                                               # needs sage.libs.gap
        True

        sage: Q = QuadraticForm(ZZ, 3, [2, 1, 2, 2, 1, 3])
        sage: sorted(Q.automorphisms())                                                 # needs sage.libs.gap
        [
        [-1  0  0]  [1 0 0]
        [ 0 -1  0]  [0 1 0]
        [ 0  0 -1], [0 0 1]
        ]
    """
    return [x.matrix() for x in self.automorphism_group()]


def number_of_automorphisms(self):
    r"""
    Return the number of automorphisms (of det `1` and `-1`) of
    the quadratic form.

    OUTPUT: integer `\geq 2`

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 3, [1, 0, 0, 1, 0, 1], unsafe_initialization=True)
        sage: Q.number_of_automorphisms()
        48

    ::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.number_of_automorphisms()
        384
        sage: 2^4 * factorial(4)
        384
    """
    try:
        return self.__number_of_automorphisms
    except AttributeError:
        self._compute_automorphisms()
        return self.__number_of_automorphisms


def set_number_of_automorphisms(self, num_autos):
    r"""
    Set the number of automorphisms to be the value given.  No error
    checking is performed, to this may lead to erroneous results.

    The fact that this result was set externally is recorded in the
    internal list of external initializations, accessible by the
    method :meth:`list_external_initializations`.

    OUTPUT: none

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1, 1, 1])
        sage: Q.list_external_initializations()
        []
        sage: Q.set_number_of_automorphisms(-3)
        sage: Q.number_of_automorphisms()
        -3
        sage: Q.list_external_initializations()
        ['number_of_automorphisms']
    """
    self.__number_of_automorphisms = num_autos
    text = 'number_of_automorphisms'
    if text not in self._external_initialization_list:
        self._external_initialization_list.append(text)
