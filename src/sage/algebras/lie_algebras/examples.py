"""
Examples of Lie Algebras

There are the following examples of Lie algebras:

- A rather comprehensive family of 3-dimensional Lie algebras
- The Lie algebra of affine transformations of the line
- All abelian Lie algebras on free modules
- The Lie algebra of upper triangular matrices
- The Lie algebra of strictly upper triangular matrices
- The symplectic derivation Lie algebra
- The rank two Heisenberg Virasoro algebra

See also
:class:`sage.algebras.lie_algebras.virasoro.LieAlgebraRegularVectorFields`
and
:class:`sage.algebras.lie_algebras.virasoro.VirasoroAlgebra` for
other examples.

AUTHORS:

- Travis Scrimshaw (07-15-2013): Initial implementation
"""
# ****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.algebras.lie_algebras.virasoro import VirasoroAlgebra
from sage.algebras.lie_algebras.rank_two_heisenberg_virasoro import RankTwoHeisenbergVirasoro
from sage.algebras.lie_algebras.symplectic_derivation import SymplecticDerivationLieAlgebra as SymplecticDerivation
from sage.algebras.lie_algebras.onsager import OnsagerAlgebra
from sage.algebras.lie_algebras.onsager import OnsagerAlgebraACE as AlternatingCentralExtensionOnsagerAlgebra
from sage.algebras.lie_algebras.affine_lie_algebra import AffineLieAlgebra as Affine
from sage.algebras.lie_algebras.classical_lie_algebra import gl
from sage.algebras.lie_algebras.classical_lie_algebra import ClassicalMatrixLieAlgebra as ClassicalMatrix


# the next 6 lines are here to silent pyflakes and lgtm warnings
assert VirasoroAlgebra
assert RankTwoHeisenbergVirasoro
assert OnsagerAlgebra
assert SymplecticDerivation
assert Affine
assert gl
assert ClassicalMatrix


def three_dimensional(R, a, b, c, d, names=['X', 'Y', 'Z']):
    r"""
    The 3-dimensional Lie algebra over a given commutative ring `R`
    with basis `\{X, Y, Z\}` subject to the relations:

    .. MATH::

        [X, Y] = aZ + dY, \quad [Y, Z] = bX, \quad [Z, X] = cY + dZ

    where `a,b,c,d \in R`.

    This is always a well-defined 3-dimensional Lie algebra, as can
    be easily proven by computation.

    EXAMPLES::

        sage: L = lie_algebras.three_dimensional(QQ, 4, 1, -1, 2)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): 2*Y + 4*Z, ('X', 'Z'): Y - 2*Z, ('Y', 'Z'): X}
        sage: TestSuite(L).run()
        sage: L = lie_algebras.three_dimensional(QQ, 1, 0, 0, 0)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): Z}
        sage: L = lie_algebras.three_dimensional(QQ, 0, 0, -1, -1)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): -Y, ('X', 'Z'): Y + Z}
        sage: L = lie_algebras.three_dimensional(QQ, 0, 1, 0, 0)
        sage: L.structure_coefficients()
        Finite family {('Y', 'Z'): X}
        sage: lie_algebras.three_dimensional(QQ, 0, 0, 0, 0)
        Abelian Lie algebra on 3 generators (X, Y, Z) over Rational Field
        sage: Q.<a,b,c,d> = PolynomialRing(QQ)
        sage: L = lie_algebras.three_dimensional(Q, a, b, c, d)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): d*Y + a*Z, ('X', 'Z'): (-c)*Y + (-d)*Z, ('Y', 'Z'): b*X}
        sage: TestSuite(L).run()
    """
    if isinstance(names, str):
        names = names.split(',')
    X = names[0]
    Y = names[1]
    Z = names[2]
    from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
    s_coeff = {(X,Y): {Z:a, Y:d}, (Y,Z): {X:b}, (Z,X): {Y:c, Z:d}}
    return LieAlgebraWithStructureCoefficients(R, s_coeff, tuple(names))

def cross_product(R, names=['X', 'Y', 'Z']):
    r"""
    The Lie algebra of `\RR^3` defined by the usual cross product
    `\times`.

    EXAMPLES::

        sage: L = lie_algebras.cross_product(QQ)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): Z, ('X', 'Z'): -Y, ('Y', 'Z'): X}
        sage: TestSuite(L).run()
    """
    L = three_dimensional(R, 1, 1, 1, 0, names=names)
    L.rename("Lie algebra of RR^3 under cross product over {}".format(R))
    return L


def three_dimensional_by_rank(R, n, a=None, names=['X', 'Y', 'Z']):
    r"""
    Return a 3-dimensional Lie algebra of rank ``n``, where `0 \leq n \leq 3`.

    Here, the *rank* of a Lie algebra `L` is defined as the dimension
    of its derived subalgebra `[L, L]`. (We are assuming that `R` is
    a field of characteristic `0`; otherwise the Lie algebras
    constructed by this function are still well-defined but no longer
    might have the correct ranks.) This is not to be confused with
    the other standard definition of a rank (namely, as the
    dimension of a Cartan subalgebra, when `L` is semisimple).

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the rank
    - ``a`` -- the deformation parameter (used for `n = 2`); this should
      be a nonzero element of `R` in order for the resulting Lie
      algebra to actually have the right rank(?)
    - ``names`` -- (optional) the generator names

    EXAMPLES::

        sage: lie_algebras.three_dimensional_by_rank(QQ, 0)
        Abelian Lie algebra on 3 generators (X, Y, Z) over Rational Field
        sage: L = lie_algebras.three_dimensional_by_rank(QQ, 1)
        sage: L.structure_coefficients()
        Finite family {('Y', 'Z'): X}
        sage: L = lie_algebras.three_dimensional_by_rank(QQ, 2, 4)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): Y, ('X', 'Z'): Y + Z}
        sage: L = lie_algebras.three_dimensional_by_rank(QQ, 2, 0)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): Y}
        sage: lie_algebras.three_dimensional_by_rank(QQ, 3)
        sl2 over Rational Field
    """
    if isinstance(names, str):
        names = names.split(',')
    names = tuple(names)

    if n == 0:
        from sage.algebras.lie_algebras.abelian import AbelianLieAlgebra
        return AbelianLieAlgebra(R, names=names)

    if n == 1:
        L = three_dimensional(R, 0, 1, 0, 0, names=names) # Strictly upper triangular matrices
        L.rename("Lie algebra of 3x3 strictly upper triangular matrices over {}".format(R))
        return L

    if n == 2:
        if a is None:
            raise ValueError("The parameter 'a' must be specified")
        X = names[0]
        Y = names[1]
        Z = names[2]
        from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
        if a == 0:
            s_coeff = {(X,Y): {Y:R.one()}, (X,Z): {Y:R(a)}}
            # Why use R(a) here if R == 0 ? Also this has rank 1.
            L = LieAlgebraWithStructureCoefficients(R, s_coeff, tuple(names))
            L.rename("Degenerate Lie algebra of dimension 3 and rank 2 over {}".format(R))
        else:
            s_coeff = {(X,Y): {Y:R.one()}, (X,Z): {Y:R.one(), Z:R.one()}}
            # a doesn't appear here :/
            L = LieAlgebraWithStructureCoefficients(R, s_coeff, tuple(names))
            L.rename("Lie algebra of dimension 3 and rank 2 with parameter {} over {}".format(a, R))
        return L

    if n == 3:
        #return sl(R, 2)
        from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
        E = names[0]
        F = names[1]
        H = names[2]
        s_coeff = { (E,F): {H:R.one()}, (H,E): {E:R(2)}, (H,F): {F:R(-2)} }
        L = LieAlgebraWithStructureCoefficients(R, s_coeff, tuple(names))
        L.rename("sl2 over {}".format(R))
        return L

    raise ValueError("Invalid rank")

def affine_transformations_line(R, names=['X', 'Y'], representation='bracket'):
    """
    The Lie algebra of affine transformations of the line.

    EXAMPLES::

        sage: L = lie_algebras.affine_transformations_line(QQ)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): Y}
        sage: X, Y = L.lie_algebra_generators()
        sage: L[X, Y] == Y
        True
        sage: TestSuite(L).run()
        sage: L = lie_algebras.affine_transformations_line(QQ, representation="matrix")
        sage: X, Y = L.lie_algebra_generators()
        sage: L[X, Y] == Y
        True
        sage: TestSuite(L).run()
    """
    if isinstance(names, str):
        names = names.split(',')
    names = tuple(names)
    if representation == 'matrix':
        from sage.matrix.matrix_space import MatrixSpace
        MS = MatrixSpace(R, 2, sparse=True)
        one = R.one()
        gens = tuple(MS({(0,i):one}) for i in range(2))
        from sage.algebras.lie_algebras.lie_algebra import LieAlgebraFromAssociative
        return LieAlgebraFromAssociative(MS, gens, names=names)
    X = names[0]
    Y = names[1]
    from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
    s_coeff = {(X,Y): {Y:R.one()}}
    L = LieAlgebraWithStructureCoefficients(R, s_coeff, names=names)
    L.rename("Lie algebra of affine transformations of a line over {}".format(R))
    return L

def abelian(R, names=None, index_set=None):
    """
    Return the abelian Lie algebra generated by ``names``.

    EXAMPLES::

        sage: lie_algebras.abelian(QQ, 'x, y, z')
        Abelian Lie algebra on 3 generators (x, y, z) over Rational Field
    """
    if isinstance(names, str):
        names = names.split(',')
    elif isinstance(names, (list, tuple)):
        names = tuple(names)
    elif names is not None:
        if index_set is not None:
            raise ValueError("invalid generator names")
        index_set = names
        names = None
    from sage.rings.infinity import infinity
    if (index_set is not None
            and not isinstance(index_set, (list, tuple))
            and index_set.cardinality() == infinity):
        from sage.algebras.lie_algebras.abelian import InfiniteDimensionalAbelianLieAlgebra
        return InfiniteDimensionalAbelianLieAlgebra(R, index_set=index_set)
    from sage.algebras.lie_algebras.abelian import AbelianLieAlgebra
    return AbelianLieAlgebra(R, names=names, index_set=index_set)

def Heisenberg(R, n, representation="structure"):
    """
    Return the rank ``n`` Heisenberg algebra in the given representation.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the rank (a nonnegative integer or infinity)
    - ``representation`` -- (default: "structure") can be one of the following:

      - ``"structure"`` -- using structure coefficients
      - ``"matrix"`` -- using matrices

    EXAMPLES::

        sage: lie_algebras.Heisenberg(QQ, 3)
        Heisenberg algebra of rank 3 over Rational Field
    """
    from sage.rings.infinity import infinity
    if n == infinity:
        from sage.algebras.lie_algebras.heisenberg import InfiniteHeisenbergAlgebra
        return InfiniteHeisenbergAlgebra(R)
    if representation == "matrix":
        from sage.algebras.lie_algebras.heisenberg import HeisenbergAlgebra_matrix
        return HeisenbergAlgebra_matrix(R, n)
    from sage.algebras.lie_algebras.heisenberg import HeisenbergAlgebra
    return HeisenbergAlgebra(R, n)

def regular_vector_fields(R):
    r"""
    Return the Lie algebra of regular vector fields on `\CC^{\times}`.

    This is also known as the Witt (Lie) algebra.

    .. SEEALSO::

        :class:`~sage.algebras.lie_algebras.virasoro.LieAlgebraRegularVectorFields`

    EXAMPLES::

        sage: lie_algebras.regular_vector_fields(QQ)
        The Lie algebra of regular vector fields over Rational Field
    """
    from sage.algebras.lie_algebras.virasoro import LieAlgebraRegularVectorFields
    return LieAlgebraRegularVectorFields(R)

witt = regular_vector_fields

def pwitt(R, p):
    r"""
    Return the `p`-Witt Lie algebra over `R`.

    INPUT:

    - ``R`` -- the base ring
    - ``p`` -- a positive integer that is `0` in ``R``

    EXAMPLES::

        sage: lie_algebras.pwitt(GF(5), 5)
        The 5-Witt Lie algebra over Finite Field of size 5
    """
    from sage.algebras.lie_algebras.virasoro import WittLieAlgebra_charp
    return WittLieAlgebra_charp(R, p)

def upper_triangular_matrices(R, n):
    r"""
    Return the Lie algebra `\mathfrak{b}_k` of `k \times k` upper
    triangular matrices.

    .. TODO::

        This implementation does not know it is finite-dimensional and
        does not know its basis.

    EXAMPLES::

        sage: L = lie_algebras.upper_triangular_matrices(QQ, 4); L
        Lie algebra of 4-dimensional upper triangular matrices over Rational Field
        sage: TestSuite(L).run()
        sage: n0, n1, n2, t0, t1, t2, t3 = L.lie_algebra_generators()
        sage: L[n2, t2] == -n2
        True

    TESTS::

        sage: L = lie_algebras.upper_triangular_matrices(QQ, 1); L
        Lie algebra of 1-dimensional upper triangular matrices over Rational Field
        sage: TestSuite(L).run()
        sage: L = lie_algebras.upper_triangular_matrices(QQ, 0); L
        Lie algebra of 0-dimensional upper triangular matrices over Rational Field
        sage: TestSuite(L).run()
    """
    from sage.matrix.matrix_space import MatrixSpace
    from sage.algebras.lie_algebras.lie_algebra import LieAlgebraFromAssociative
    MS = MatrixSpace(R, n, sparse=True)
    one = R.one()
    names = tuple('n{}'.format(i) for i in range(n-1))
    names += tuple('t{}'.format(i) for i in range(n))
    gens = [MS({(i,i+1):one}) for i in range(n-1)]
    gens += [MS({(i,i):one}) for i in range(n)]
    L = LieAlgebraFromAssociative(MS, gens, names=names)
    L.rename("Lie algebra of {}-dimensional upper triangular matrices over {}".format(n, L.base_ring()))
    return L

def strictly_upper_triangular_matrices(R, n):
    r"""
    Return the Lie algebra `\mathfrak{n}_k` of strictly `k \times k` upper
    triangular matrices.

    .. TODO::

        This implementation does not know it is finite-dimensional and
        does not know its basis.

    EXAMPLES::

        sage: L = lie_algebras.strictly_upper_triangular_matrices(QQ, 4); L
        Lie algebra of 4-dimensional strictly upper triangular matrices over Rational Field
        sage: TestSuite(L).run()
        sage: n0, n1, n2 = L.lie_algebra_generators()
        sage: L[n2, n1]
        [ 0  0  0  0]
        [ 0  0  0 -1]
        [ 0  0  0  0]
        [ 0  0  0  0]

    TESTS::

        sage: L = lie_algebras.strictly_upper_triangular_matrices(QQ, 1); L
        Lie algebra of 1-dimensional strictly upper triangular matrices over Rational Field
        sage: TestSuite(L).run()
        sage: L = lie_algebras.strictly_upper_triangular_matrices(QQ, 0); L
        Lie algebra of 0-dimensional strictly upper triangular matrices over Rational Field
        sage: TestSuite(L).run()
    """
    from sage.matrix.matrix_space import MatrixSpace
    from sage.algebras.lie_algebras.lie_algebra import LieAlgebraFromAssociative
    MS = MatrixSpace(R, n, sparse=True)
    one = R.one()
    names = tuple('n{}'.format(i) for i in range(n-1))
    gens = tuple(MS({(i,i+1):one}) for i in range(n-1))
    L = LieAlgebraFromAssociative(MS, gens, names=names)
    L.rename("Lie algebra of {}-dimensional strictly upper triangular matrices over {}".format(n, L.base_ring()))
    return L

#####################################################################
## Classical Lie algebras


def sl(R, n, representation='bracket'):
    r"""
    The Lie algebra `\mathfrak{sl}_n`.

    The Lie algebra `\mathfrak{sl}_n` is the type `A_{n-1}` Lie algebra
    and is finite dimensional. As a matrix Lie algebra, it is given by
    the set of all `n \times n` matrices with trace 0.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the size of the matrix
    - ``representation`` -- (default: ``'bracket'``) can be one of
      the following:

      * ``'bracket'`` - use brackets and the Chevalley basis
      * ``'matrix'`` - use matrices

    EXAMPLES:

    We first construct `\mathfrak{sl}_2` using the Chevalley basis::

        sage: sl2 = lie_algebras.sl(QQ, 2); sl2
        Lie algebra of ['A', 1] in the Chevalley basis
        sage: E,F,H = sl2.gens()
        sage: E.bracket(F) == H
        True
        sage: H.bracket(E) == 2*E
        True
        sage: H.bracket(F) == -2*F
        True

    We now construct `\mathfrak{sl}_2` as a matrix Lie algebra::

        sage: sl2 = lie_algebras.sl(QQ, 2, representation='matrix')
        sage: E,F,H = sl2.gens()
        sage: E.bracket(F) == H
        True
        sage: H.bracket(E) == 2*E
        True
        sage: H.bracket(F) == -2*F
        True
    """
    if representation == 'bracket':
        from sage.algebras.lie_algebras.classical_lie_algebra import LieAlgebraChevalleyBasis
        return LieAlgebraChevalleyBasis(R, ['A', n-1])
    if representation == 'matrix':
        from sage.algebras.lie_algebras.classical_lie_algebra import sl as sl_matrix
        return sl_matrix(R, n)
    raise ValueError("invalid representation")


def su(R, n, representation='matrix'):
    r"""
    The Lie algebra `\mathfrak{su}_n`.

    The Lie algebra `\mathfrak{su}_n` is the compact real form of the
    type `A_{n-1}` Lie algebra and is finite-dimensional. As a matrix
    Lie algebra, it is given by the set of all `n \times n` skew-Hermitian
    matrices with trace 0.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the size of the matrix
    - ``representation`` -- (default: ``'matrix'``) can be one of
      the following:

      * ``'bracket'`` - use brackets and the Chevalley basis
      * ``'matrix'`` - use matrices

    EXAMPLES:

    We construct `\mathfrak{su}_2`, where the default is as a
    matrix Lie algebra::

        sage: su2 = lie_algebras.su(QQ, 2)
        sage: E,H,F = su2.basis()
        sage: E.bracket(F) == 2*H
        True
        sage: H.bracket(E) == 2*F
        True
        sage: H.bracket(F) == -2*E
        True

    Since `\mathfrak{su}_n` is the same as the type `A_{n-1}` Lie algebra,
    the bracket is the same as :func:`sl`::

        sage: su2 = lie_algebras.su(QQ, 2, representation='bracket')
        sage: su2 is lie_algebras.sl(QQ, 2, representation='bracket')
        True
    """
    if representation == 'bracket':
        from sage.algebras.lie_algebras.classical_lie_algebra import LieAlgebraChevalleyBasis
        return LieAlgebraChevalleyBasis(R, ['A', n-1])
    if representation == 'matrix':
        from sage.algebras.lie_algebras.classical_lie_algebra import MatrixCompactRealForm
        from sage.combinat.root_system.cartan_type import CartanType
        return MatrixCompactRealForm(R, CartanType(['A', n-1]))
    raise ValueError("invalid representation")

def so(R, n, representation='bracket'):
    r"""
    The Lie algebra `\mathfrak{so}_n`.

    The Lie algebra `\mathfrak{so}_n` is the type `B_k` Lie algebra
    if `n = 2k - 1` or the type `D_k` Lie algebra if `n = 2k`, and in
    either case is finite dimensional. As a matrix Lie algebra, it
    is given by the set of all real anti-symmetric `n \times n` matrices.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the size of the matrix
    - ``representation`` -- (default: ``'bracket'``) can be one of
      the following:

      * ``'bracket'`` - use brackets and the Chevalley basis
      * ``'matrix'`` - use matrices

    EXAMPLES:

    We first construct `\mathfrak{so}_5` using the Chevalley basis::

        sage: so5 = lie_algebras.so(QQ, 5); so5
        Lie algebra of ['B', 2] in the Chevalley basis
        sage: E1,E2, F1,F2, H1,H2 = so5.gens()
        sage: so5([E1, [E1, E2]])
        0
        sage: X = so5([E2, [E2, E1]]); X
        -2*E[alpha[1] + 2*alpha[2]]
        sage: H1.bracket(X)
        0
        sage: H2.bracket(X)
        -4*E[alpha[1] + 2*alpha[2]]
        sage: so5([H1, [E1, E2]])
        -E[alpha[1] + alpha[2]]
        sage: so5([H2, [E1, E2]])
        0

    We do the same construction of `\mathfrak{so}_4` using the Chevalley
    basis::

        sage: so4 = lie_algebras.so(QQ, 4); so4
        Lie algebra of ['D', 2] in the Chevalley basis
        sage: E1,E2, F1,F2, H1,H2 = so4.gens()
        sage: H1.bracket(E1)
        2*E[alpha[1]]
        sage: H2.bracket(E1) == so4.zero()
        True
        sage: E1.bracket(E2) == so4.zero()
        True

    We now construct `\mathfrak{so}_4` as a matrix Lie algebra::

        sage: sl2 = lie_algebras.sl(QQ, 2, representation='matrix')
        sage: E1,E2, F1,F2, H1,H2 = so4.gens()
        sage: H2.bracket(E1) == so4.zero()
        True
        sage: E1.bracket(E2) == so4.zero()
        True
    """
    if representation == 'bracket':
        from sage.algebras.lie_algebras.classical_lie_algebra import LieAlgebraChevalleyBasis
        if n % 2 == 0:
            return LieAlgebraChevalleyBasis(R, ['D', n//2])
        else:
            return LieAlgebraChevalleyBasis(R, ['B', (n-1)//2])
    if representation == 'matrix':
        from sage.algebras.lie_algebras.classical_lie_algebra import so as so_matrix
        return so_matrix(R, n)
    raise ValueError("invalid representation")


def sp(R, n, representation='bracket'):
    r"""
    The Lie algebra `\mathfrak{sp}_n`.

    The Lie algebra `\mathfrak{sp}_n` where `n = 2k` is the type `C_k`
    Lie algebra and is finite dimensional. As a matrix Lie algebra, it
    is given by the set of all matrices `X` that satisfy the equation:

    .. MATH::

        X^T M - M X = 0

    where

    .. MATH::

        M = \begin{pmatrix}
        0 & I_k \\
        -I_k & 0
        \end{pmatrix}.

    This is the Lie algebra of type `C_k`.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the size of the matrix
    - ``representation`` -- (default: ``'bracket'``) can be one of
      the following:

      * ``'bracket'`` - use brackets and the Chevalley basis
      * ``'matrix'`` - use matrices

    EXAMPLES:

    We first construct `\mathfrak{sp}_4` using the Chevalley basis::

        sage: sp4 = lie_algebras.sp(QQ, 4); sp4
        Lie algebra of ['C', 2] in the Chevalley basis
        sage: E1,E2, F1,F2, H1,H2 = sp4.gens()
        sage: sp4([E2, [E2, E1]])
        0
        sage: X = sp4([E1, [E1, E2]]); X
        2*E[2*alpha[1] + alpha[2]]
        sage: H1.bracket(X)
        4*E[2*alpha[1] + alpha[2]]
        sage: H2.bracket(X)
        0
        sage: sp4([H1, [E1, E2]])
        0
        sage: sp4([H2, [E1, E2]])
        -E[alpha[1] + alpha[2]]

    We now construct `\mathfrak{sp}_4` as a matrix Lie algebra::

        sage: sp4 = lie_algebras.sp(QQ, 4, representation='matrix'); sp4
        Symplectic Lie algebra of rank 4 over Rational Field
        sage: E1,E2, F1,F2, H1,H2 = sp4.gens()
        sage: H1.bracket(E1)
        [ 0  2  0  0]
        [ 0  0  0  0]
        [ 0  0  0  0]
        [ 0  0 -2  0]
        sage: sp4([E1, [E1, E2]])
        [0 0 2 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
    """
    if n % 2:
        raise ValueError("n must be even")
    if representation == 'bracket':
        from sage.algebras.lie_algebras.classical_lie_algebra import LieAlgebraChevalleyBasis
        return LieAlgebraChevalleyBasis(R, ['C', n//2])
    if representation == 'matrix':
        from sage.algebras.lie_algebras.classical_lie_algebra import sp as sp_matrix
        return sp_matrix(R, n)
    raise ValueError("invalid representation")
