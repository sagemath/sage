# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.modules
"""
Hard Lattice Generator

This module contains lattice related functions relevant in
cryptography.

Feel free to add more functionality.

AUTHORS:

- Richard Lindner <rlindner@cdc.informatik.tu-darmstadt.de>

- Michael Schneider <mischnei@cdc.informatik.tu-darmstadt.de>
"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.polynomial.polynomial_ring import PolynomialRing_general


def gen_lattice(type='modular', n=4, m=8, q=11, seed=None,
                quotient=None, dual=False, ntl=False, lattice=False):
    r"""
    This function generates different types of integral lattice bases
    of row vectors relevant in cryptography.

    Randomness can be set either with ``seed``, or by using
    :func:`sage.misc.randstate.set_random_seed`.

    INPUT:

    - ``type`` -- one of the following strings

      - ``'modular'`` -- default; a class of lattices for which
        asymptotic worst-case to average-case connections hold. For
        more refer to [Aj1996]_
      - ``'random'`` -- special case of modular (n=1); a dense class
        of lattice used for testing basis reduction algorithms
        proposed by Goldstein and Mayer [GM2002]_
      - ``'ideal'`` -- special case of modular; allows for a more
        compact representation proposed by [LM2006]_
      - ``'cyclotomic'`` -- special case of ideal; allows for
        efficient processing proposed by [LM2006]_

    - ``n`` -- determinant size, primal: `det(L) = q^n`, dual: `det(L) = q^{m-n}`.
      For ideal lattices this is also the degree of the quotient polynomial.

    - ``m`` -- lattice dimension, `L \subseteq Z^m`

    - ``q`` -- coefficient size, `q-Z^m \subseteq L`

    - ``seed`` -- randomness seed

    - ``quotient`` -- for the type ``'ideal'``, this determines the quotient
      polynomial. Ignored for all other types

    - ``dual`` -- set this flag if you want a basis for `q-dual(L)`, for example
      for Regev's LWE bases [Reg2005]_

    - ``ntl`` -- set this flag if you want the lattice basis in NTL readable
      format

    - ``lattice`` -- set this flag if you want a
      :class:`FreeModule_submodule_with_basis_integer` object instead
      of an integer matrix representing the basis

    OUTPUT: ``B`` a unique size-reduced triangular (primal: lower_left,
    dual: lower_right) basis of row vectors for the lattice in question

    EXAMPLES:

    Modular basis::

        sage: sage.crypto.gen_lattice(m=10, seed=42)
        [11  0  0  0  0  0  0  0  0  0]
        [ 0 11  0  0  0  0  0  0  0  0]
        [ 0  0 11  0  0  0  0  0  0  0]
        [ 0  0  0 11  0  0  0  0  0  0]
        [ 2  4  3  5  1  0  0  0  0  0]
        [ 1 -5 -4  2  0  1  0  0  0  0]
        [-4  3 -1  1  0  0  1  0  0  0]
        [-2 -3 -4 -1  0  0  0  1  0  0]
        [-5 -5  3  3  0  0  0  0  1  0]
        [-4 -3  2 -5  0  0  0  0  0  1]

    Random basis::

        sage: sage.crypto.gen_lattice(type='random', n=1, m=10, q=11^4, seed=42)
        [14641     0     0     0     0     0     0     0     0     0]
        [  431     1     0     0     0     0     0     0     0     0]
        [-4792     0     1     0     0     0     0     0     0     0]
        [ 1015     0     0     1     0     0     0     0     0     0]
        [-3086     0     0     0     1     0     0     0     0     0]
        [-5378     0     0     0     0     1     0     0     0     0]
        [ 4769     0     0     0     0     0     1     0     0     0]
        [-1159     0     0     0     0     0     0     1     0     0]
        [ 3082     0     0     0     0     0     0     0     1     0]
        [-4580     0     0     0     0     0     0     0     0     1]

    Ideal bases with quotient `x^n-1`, `m=2*n` are NTRU bases::

        sage: sage.crypto.gen_lattice(type='ideal', seed=42, quotient=x^4 - 1)          # needs sage.symbolic
        [11  0  0  0  0  0  0  0]
        [ 0 11  0  0  0  0  0  0]
        [ 0  0 11  0  0  0  0  0]
        [ 0  0  0 11  0  0  0  0]
        [-3 -3 -2  4  1  0  0  0]
        [ 4 -3 -3 -2  0  1  0  0]
        [-2  4 -3 -3  0  0  1  0]
        [-3 -2  4 -3  0  0  0  1]

    Ideal bases also work with polynomials::

        sage: R.<t> = PolynomialRing(ZZ)
        sage: sage.crypto.gen_lattice(type='ideal', seed=1234, quotient=t^4 - 1)        # needs sage.libs.pari
        [11  0  0  0  0  0  0  0]
        [ 0 11  0  0  0  0  0  0]
        [ 0  0 11  0  0  0  0  0]
        [ 0  0  0 11  0  0  0  0]
        [-3  4  1  4  1  0  0  0]
        [ 4 -3  4  1  0  1  0  0]
        [ 1  4 -3  4  0  0  1  0]
        [ 4  1  4 -3  0  0  0  1]

    Cyclotomic bases with n=2^k are SWIFFT bases::

        sage: sage.crypto.gen_lattice(type='cyclotomic', seed=42)                       # needs sage.libs.pari
        [11  0  0  0  0  0  0  0]
        [ 0 11  0  0  0  0  0  0]
        [ 0  0 11  0  0  0  0  0]
        [ 0  0  0 11  0  0  0  0]
        [-3 -3 -2  4  1  0  0  0]
        [-4 -3 -3 -2  0  1  0  0]
        [ 2 -4 -3 -3  0  0  1  0]
        [ 3  2 -4 -3  0  0  0  1]

    Dual modular bases are related to Regev's famous public-key
    encryption [Reg2005]_::

        sage: sage.crypto.gen_lattice(type='modular', m=10, seed=42, dual=True)
        [ 0  0  0  0  0  0  0  0  0 11]
        [ 0  0  0  0  0  0  0  0 11  0]
        [ 0  0  0  0  0  0  0 11  0  0]
        [ 0  0  0  0  0  0 11  0  0  0]
        [ 0  0  0  0  0 11  0  0  0  0]
        [ 0  0  0  0 11  0  0  0  0  0]
        [ 0  0  0  1 -5 -2 -1  1 -3  5]
        [ 0  0  1  0 -3  4  1  4 -3 -2]
        [ 0  1  0  0 -4  5 -3  3  5  3]
        [ 1  0  0  0 -2 -1  4  2  5  4]

    Relation of primal and dual bases::

        sage: B_primal = sage.crypto.gen_lattice(m=10, q=11, seed=42)
        sage: B_dual = sage.crypto.gen_lattice(m=10, q=11, seed=42, dual=True)
        sage: B_dual_alt = transpose(11*B_primal.inverse()).change_ring(ZZ)
        sage: B_dual_alt.hermite_form() == B_dual.hermite_form()
        True

    TESTS:

    Test some bad quotient polynomials::

        sage: sage.crypto.gen_lattice(type='ideal', seed=1234, quotient=cos(x))         # needs sage.symbolic
        Traceback (most recent call last):
        ...
        TypeError: self must be a numeric expression
        sage: sage.crypto.gen_lattice(type='ideal', seed=1234, quotient=x^23-1)         # needs sage.symbolic
        Traceback (most recent call last):
        ...
        ValueError: ideal basis requires n = quotient.degree()
        sage: R.<u,v> = ZZ[]
        sage: sage.crypto.gen_lattice(type='ideal', seed=1234, quotient=u+v)
        Traceback (most recent call last):
        ...
        TypeError: quotient should be a univariate polynomial

    We are testing output format choices::

        sage: sage.crypto.gen_lattice(m=10, q=11, seed=42)
        [11  0  0  0  0  0  0  0  0  0]
        [ 0 11  0  0  0  0  0  0  0  0]
        [ 0  0 11  0  0  0  0  0  0  0]
        [ 0  0  0 11  0  0  0  0  0  0]
        [ 2  4  3  5  1  0  0  0  0  0]
        [ 1 -5 -4  2  0  1  0  0  0  0]
        [-4  3 -1  1  0  0  1  0  0  0]
        [-2 -3 -4 -1  0  0  0  1  0  0]
        [-5 -5  3  3  0  0  0  0  1  0]
        [-4 -3  2 -5  0  0  0  0  0  1]

        sage: sage.crypto.gen_lattice(m=10, q=11, seed=42, ntl=True)
        [
        [11 0 0 0 0 0 0 0 0 0]
        [0 11 0 0 0 0 0 0 0 0]
        [0 0 11 0 0 0 0 0 0 0]
        [0 0 0 11 0 0 0 0 0 0]
        [2 4 3 5 1 0 0 0 0 0]
        [1 -5 -4 2 0 1 0 0 0 0]
        [-4 3 -1 1 0 0 1 0 0 0]
        [-2 -3 -4 -1 0 0 0 1 0 0]
        [-5 -5 3 3 0 0 0 0 1 0]
        [-4 -3 2 -5 0 0 0 0 0 1]
        ]

        sage: sage.crypto.gen_lattice(m=10, q=11, seed=42, lattice=True)                # needs fpylll
        Free module of degree 10 and rank 10 over Integer Ring
        User basis matrix:
        [ 0  0  1  1  0 -1 -1 -1  1  0]
        [-1  1  0  1  0  1  1  0  1  1]
        [-1  0  0  0 -1  1  1 -2  0  0]
        [-1 -1  0  1  1  0  0  1  1 -1]
        [ 1  0 -1  0  0  0 -2 -2  0  0]
        [ 2 -1  0  0  1  0  1  0  0 -1]
        [-1  1 -1  0  1 -1  1  0 -1 -2]
        [ 0  0 -1  3  0  0  0 -1 -1 -1]
        [ 0 -1  0 -1  2  0 -1  0  0  2]
        [ 0  1  1  0  1  1 -2  1 -1 -2]
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
    from sage.matrix.constructor import identity_matrix, block_matrix
    from sage.matrix.matrix_space import MatrixSpace
    from sage.rings.integer_ring import ZZ
    if seed is not None:
        from sage.misc.randstate import set_random_seed
        set_random_seed(seed)

    if type == 'random':
        if n != 1:
            raise ValueError('random bases require n = 1')

    ZZ_q = IntegerModRing(q)
    A = identity_matrix(ZZ_q, n)

    if type == 'random' or type == 'modular':
        R = MatrixSpace(ZZ_q, m-n, n)
        A = A.stack(R.random_element())

    elif type == 'ideal':
        if quotient is None:
            raise ValueError('ideal bases require a quotient polynomial')
        try:
            quotient = quotient.change_ring(ZZ_q)
        except (AttributeError, TypeError):
            quotient = quotient.polynomial(base_ring=ZZ_q)

        P = quotient.parent()
        # P should be a univariate polynomial ring over ZZ_q
        if not isinstance(P, PolynomialRing_general):
            raise TypeError("quotient should be a univariate polynomial")
        assert P.base_ring() is ZZ_q

        if quotient.degree() != n:
            raise ValueError('ideal basis requires n = quotient.degree()')
        R = P.quotient(quotient)
        for i in range(m//n):
            A = A.stack(R.random_element().matrix())

    elif type == 'cyclotomic':
        from sage.arith.misc import euler_phi
        from sage.misc.functional import cyclotomic_polynomial

        # we assume that n+1 <= min( euler_phi^{-1}(n) ) <= 2*n
        found = False
        for k in range(2*n,n,-1):
            if euler_phi(k) == n:
                found = True
                break
        if not found:
            raise ValueError("cyclotomic bases require that n "
                       "is an image of Euler's totient function")

        R = ZZ_q['x'].quotient(cyclotomic_polynomial(k, 'x'), 'x')
        for i in range(m//n):
            A = A.stack(R.random_element().matrix())

    # switch from representatives 0,...,(q-1) to (1-q)/2,....,(q-1)/2
    def minrep(a):
        if abs(a-q) < abs(a):
            return a-q
        else:
            return a
    A_prime = A[n:m].lift().apply_map(minrep)

    if not dual:
        B = block_matrix([[ZZ(q), ZZ.zero()], [A_prime, ZZ.one()] ],
                         subdivide=False)
    else:
        B = block_matrix([[ZZ.one(), -A_prime.transpose()],
            [ZZ.zero(), ZZ(q)]], subdivide=False)
        for i in range(m//2):
            B.swap_rows(i,m-i-1)

    if ntl and lattice:
        raise ValueError("Cannot specify ntl=True and lattice=True "
                         "at the same time")

    if ntl:
        return B._ntl_()
    elif lattice:
        from sage.modules.free_module_integer import IntegerLattice
        return IntegerLattice(B)
    else:
        return B
