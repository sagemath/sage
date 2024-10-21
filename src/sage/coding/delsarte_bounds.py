# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.numerical.mip
r"""
Delsarte (or linear programming) bounds

This module provides LP upper bounds for the parameters of codes, introduced in
by P. Delsarte in [De1973]_.

The exact LP solver PPL is used by default, ensuring that no
rounding/overflow problems occur.

AUTHORS:

- Dmitrii V. (Dima) Pasechnik (2012-10): initial implementation

- Dmitrii V. (Dima) Pasechnik (2015, 2021): minor fixes

- Charalampos Kokkalis (2021): Eberlein polynomials, general Q matrix codes
"""
# ****************************************************************************
#       Copyright (C) 2012, 2021 Dima Pasechnik <dimpase@gmail.com>
#       Copyright (C) 2021 Charalampos Kokkalis
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


def krawtchouk(n, q, l, x, check=True):
    r"""
    Compute `K^{n,q}_l(x)`, the Krawtchouk (a.k.a. Kravchuk) polynomial.

    See :wikipedia:`Kravchuk_polynomials`.

    It is defined by the generating function

    .. MATH::

        (1+(q-1)z)^{n-x}(1-z)^x=\sum_{l} K^{n,q}_l(x)z^l

    and is equal to

    .. MATH::

        K^{n,q}_l(x)=\sum_{j=0}^l (-1)^j (q-1)^{(l-j)} \binom{x}{j} \binom{n-x}{l-j}.

    INPUT:

    - ``n``, ``q``, ``x`` -- arbitrary numbers

    - ``l`` -- nonnegative integer

    - ``check`` -- check the input for correctness. ``True`` by
      default. Otherwise, pass it as it is. Use ``check=False`` at
      your own risk.

    .. SEEALSO::

        :class:`Symbolic Krawtchouk polynomials
        <sage.functions.orthogonal_polys.Func_krawtchouk>` `\tilde{K}_l(x; n, p)`
        which are related by

        .. MATH::

            (-q)^l K^{n,q^{-1}}_l(x) = \tilde{K}_l(x; n, 1-q).

    EXAMPLES::

        sage: codes.bounds.krawtchouk(24,2,5,4)
        2224
        sage: codes.bounds.krawtchouk(12300,4,5,6)
        567785569973042442072

    TESTS:

    Check that the bug reported on :issue:`19561` is fixed::

        sage: codes.bounds.krawtchouk(3,2,3,3)
        -1
        sage: codes.bounds.krawtchouk(int(3),int(2),int(3),int(3))
        -1
        sage: codes.bounds.krawtchouk(int(3),int(2),int(3),int(3),check=False)
        -1.0
        sage: codes.bounds.krawtchouk(24,2,5,4)
        2224

    Other unusual inputs::

        sage: codes.bounds.krawtchouk(sqrt(5),1-I*sqrt(3),3,55.3).n()                   # needs sage.symbolic
        211295.892797... + 1186.42763...*I
        sage: codes.bounds.krawtchouk(-5/2,7*I,3,-1/10)                                 # needs sage.symbolic
        480053/250*I - 357231/400
        sage: codes.bounds.krawtchouk(1,1,-1,1)
        Traceback (most recent call last):
        ...
        ValueError: l must be a nonnegative integer
        sage: codes.bounds.krawtchouk(1,1,3/2,1)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer
    """
    from sage.arith.misc import binomial
    from sage.arith.srange import srange
    # Use the expression in equation (55) of MacWilliams & Sloane, pg 151
    # We write jth term = some_factor * (j-1)th term
    if check:
        from sage.rings.integer_ring import ZZ
        l0 = ZZ(l)
        if l0 != l or l0 < 0:
            raise ValueError('l must be a nonnegative integer')
        l = l0
    kraw = jth_term = (q-1)**l * binomial(n, l)  # j=0
    for j in srange(1, l+1):
        jth_term *= -q*(l-j+1)*(x-j+1)/((q-1)*j*(n-j+1))
        kraw += jth_term
    return kraw


def eberlein(n, w, k, u, check=True):
    r"""
    Compute `E^{w,n}_k(x)`, the Eberlein polynomial.

    See :wikipedia:`Eberlein_polynomials`.

    It is defined as:

    .. MATH::

        E^{w,n}_k(u)=\sum_{j=0}^k (-1)^j \binom{u}{j} \binom{w-u}{k-j}
        \binom{n-w-u}{k-j},

    INPUT:

    - ``w``, ``k``, ``x`` -- arbitrary numbers

    - ``n`` -- nonnegative integer

    - ``check`` -- check the input for correctness. ``True`` by
      default. Otherwise, pass it as it is. Use ``check=False`` at
      your own risk.

    EXAMPLES::

        sage: codes.bounds.eberlein(24,10,2,6)
        -9

    TESTS:

    check normal inputs (various formats for arguments) ::

        sage: codes.bounds.eberlein(24,10,2,6)
        -9
        sage: codes.bounds.eberlein(int(24),int(10),int(2),int(6))
        -9
        sage: codes.bounds.eberlein(int(24),int(10),int(2),int(6),check=False)
        -9

    unusual inputs ::

        sage: codes.bounds.eberlein(-1,1,1,1)
        Traceback (most recent call last):
        ...
        ValueError: l must be a nonnegative integer
        sage: codes.bounds.eberlein(1,1,3/2,1)
        Traceback (most recent call last):
        ...
        TypeError: either m or x-m must be an integer
    """
    from sage.arith.misc import binomial
    from sage.arith.srange import srange

    if 2 * w > n:
        return eberlein(n, n - w, k, u)

    if check:
        from sage.rings.integer_ring import ZZ
        n0 = ZZ(n)
        if n0 != n or n0 < 0:
            raise ValueError('l must be a nonnegative integer')
        n = n0

    return sum([(-1)**j*binomial(u, j)*binomial(w-u, k-j)*binomial(n-w-u, k-j)
                for j in srange(k + 1)])


def _delsarte_LP_building(n, d, d_star, q, isinteger, solver, maxc=0):
    r"""
    LP builder - common for the two functions; not exported.

    EXAMPLES::

        sage: from sage.coding.delsarte_bounds import _delsarte_LP_building
        sage: _,p=_delsarte_LP_building(7, 3, 0, 2, False, "PPL")
        sage: p.show()
        Maximization:
          x_0 + x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7
        Constraints:
          constraint_0: 1 <= x_0 <= 1
          constraint_1: 0 <= x_1 <= 0
          constraint_2: 0 <= x_2 <= 0
          constraint_3: -7 x_0 - 5 x_1 - 3 x_2 - x_3 + x_4 + 3 x_5 + 5 x_6 + 7 x_7 <= 0
          constraint_4: -7 x_0 - 5 x_1 - 3 x_2 - x_3 + x_4 + 3 x_5 + 5 x_6 + 7 x_7 <= 0
          ...
          constraint_16: - x_0 + x_1 - x_2 + x_3 - x_4 + x_5 - x_6 + x_7 <= 0
        Variables:
          x_0 is a continuous variable (min=0, max=+oo)
          ...
          x_7 is a continuous variable (min=0, max=+oo)
    """
    from sage.numerical.mip import MixedIntegerLinearProgram

    p = MixedIntegerLinearProgram(maximization=True, solver=solver)
    A = p.new_variable(integer=isinteger, nonnegative=True)
    p.set_objective(sum([A[r] for r in range(n + 1)]))
    p.add_constraint(A[0] == 1)
    for i in range(1, d):
        p.add_constraint(A[i] == 0)
    for j in range(1, n + 1):
        rhs = sum([krawtchouk(n, q, j, r, check=False) * A[r]
                   for r in range(n + 1)])
        p.add_constraint(0 <= rhs)
        if j >= d_star:
            p.add_constraint(0 <= rhs)
        else:  # rhs is proportional to j-th weight of the dual code
            p.add_constraint(0 == rhs)

    if maxc > 0:
        p.add_constraint(sum([A[r] for r in range(n + 1)]), max=maxc)
    return A, p


def _delsarte_cwc_LP_building(n, d, w, solver, isinteger):
    r"""
    LP builder for Delsarte's LP for constant weight codes.

    It is used in :func:`delsarte_bound_constant_weight_code`; not exported.

    INPUT:

    - ``n`` -- the code length

    - ``w`` -- the weight of the code

    - ``d`` -- the (lower bound on) minimal distance of the code

    - ``solver`` -- the LP/ILP solver to be used. Defaults to
      ``PPL``. It is arbitrary precision, thus there will be no
      rounding errors. With other solvers (see
      :class:`MixedIntegerLinearProgram` for the list), you are on
      your own!

    - ``isinteger`` -- if ``True``, uses an integer programming solver
      (ILP), rather that an LP solver. Can be very slow if set to
      ``True``.

    EXAMPLES::

        sage: from sage.coding.delsarte_bounds import _delsarte_cwc_LP_building
        sage: _,p=_delsarte_cwc_LP_building(17, 4, 3, "PPL", False)
        sage: p.show()
        Maximization:
          x_0 + x_1 + 1
        Constraints:
          constraint_0: -1 <= 4/21 x_0 - 3/14 x_1
          constraint_1: -1 <= -23/273 x_0 + 3/91 x_1
          constraint_2: -1 <= 1/91 x_0 - 1/364 x_1
        Variables:
          x_0 is a continuous variable (min=0, max=+oo)
          x_1 is a continuous variable (min=0, max=+oo)
    """
    from sage.arith.misc import binomial
    from sage.numerical.mip import MixedIntegerLinearProgram

    p = MixedIntegerLinearProgram(maximization=True, solver=solver)
    A = p.new_variable(integer=isinteger, nonnegative=True)
    p.set_objective(sum([A[2*r] for r in range(d//2, w+1)]) + 1)

    def _q(k, i):
        mu_i = 1
        v_i = binomial(w, i)*binomial(n-w, i)
        return mu_i*eberlein(n, w, i, k)/v_i

    for k in range(1, w+1):
        p.add_constraint(sum([A[2*i]*_q(k, i) for i in range(d//2, w+1)]),
                         min=-1)

    return A, p


def delsarte_bound_constant_weight_code(n, d, w, return_data=False,
                                        solver='PPL', isinteger=False):
    r"""
    Find the Delsarte bound on a constant weight code.

    Find the Delsarte bound on a constant weight code of weight ``w``, length
    ``n``, lower bound on minimal distance ``d``.

    INPUT:

    - ``n`` -- the code length

    - ``d`` -- the (lower bound on) minimal distance of the code

    - ``w`` -- the weight of the code

    - ``return_data`` -- if ``True``, return a triple
      ``(W,LP,bound)``, where ``W`` is a weights vector, and ``LP``
      the Delsarte upper bound LP; both of them are Sage LP data.
      ``W`` need not be a weight distribution of a code.

    - ``solver`` -- the LP/ILP solver to be used. Defaults to
      ``PPL``. It is arbitrary precision, thus there will be no
      rounding errors. With other solvers (see
      :class:`MixedIntegerLinearProgram` for the list), you are on
      your own!

    - ``isinteger`` -- if ``True``, uses an integer programming solver
      (ILP), rather that an LP solver. Can be very slow if set to
      ``True``.

    EXAMPLES:

    The bound on the size of codes of length 17, weight 3, and minimal distance 4::

       sage: codes.bounds.delsarte_bound_constant_weight_code(17, 4, 3)
       45
       sage: a, p, val = codes.bounds.delsarte_bound_constant_weight_code(17, 4, 3, return_data=True)
       sage: [j for i,j in p.get_values(a).items()]
       [21, 70/3]

    The stricter bound (using ILP) on codes of length 17, weight 3, and minimal
    distance 4::

       sage: codes.bounds.delsarte_bound_constant_weight_code(17, 4, 3, isinteger=True)
       43
    """
    from sage.numerical.mip import MIPSolverException

    if d < 4:
        raise ValueError("Violated constraint d>=4 for "
                         "Binary Constant Weight Codes")

    if d >= 2*w or 2*w > n:
        raise ValueError("Violated constraint d<2w<=n for "
                         "Binary Constant Weight Codes")

    # minimum distance is even => if there is an odd lower bound on d we can
    # increase it by 1
    if d % 2:
        d += 1

    A, p = _delsarte_cwc_LP_building(n, d, w, solver, isinteger)
    try:
        bd = p.solve()
    except MIPSolverException as exc:
        print(f"Solver exception: {exc}")
        return (A, p, False) if return_data else False

    return (A, p, bd) if return_data else int(bd)


def delsarte_bound_hamming_space(n, d, q, return_data=False,
                                 solver='PPL', isinteger=False):
    r"""
    Find the Delsarte bound on codes in ``H_q^n`` of minimal distance ``d``.

    Find the Delsarte bound [De1973]_ on the size of codes in
    the Hamming space ``H_q^n`` of minimal distance ``d``.

    INPUT:

    - ``n`` -- the code length

    - ``d`` -- the (lower bound on) minimal distance of the code

    - ``q`` -- the size of the alphabet

    - ``return_data`` -- if ``True``, return a triple
      ``(W,LP,bound)``, where ``W`` is a weights vector, and ``LP``
      the Delsarte upper bound LP; both of them are Sage LP data.
      ``W`` need not be a weight distribution of a code.

    - ``solver`` -- the LP/ILP solver to be used. Defaults to
      ``PPL``. It is arbitrary precision, thus there will be no
      rounding errors. With other solvers (see
      :class:`MixedIntegerLinearProgram` for the list), you are on
      your own!

    - ``isinteger`` -- if ``True``, uses an integer programming solver
      (ILP), rather that an LP solver. Can be very slow if set to
      ``True``.

    EXAMPLES:

    The bound on the size of the `\GF{2}`-codes of length 11 and minimal distance 6::

       sage: codes.bounds.delsarte_bound_hamming_space(11, 6, 2)
       12
       sage: a, p, val = codes.bounds.delsarte_bound_hamming_space(11, 6, 2, return_data=True)
       sage: [j for i,j in p.get_values(a).items()]
       [1, 0, 0, 0, 0, 0, 11, 0, 0, 0, 0, 0]

    The bound on the size of the `\GF{2}`-codes of length 24 and minimal distance
    8, i.e. parameters of the extended binary Golay code::

       sage: a,p,x = codes.bounds.delsarte_bound_hamming_space(24,8,2,return_data=True)
       sage: x
       4096
       sage: [j for i,j in p.get_values(a).items()]
       [1, 0, 0, 0, 0, 0, 0, 0, 759, 0, 0, 0, 2576, 0, 0, 0, 759, 0, 0, 0, 0, 0, 0, 0, 1]

    The bound on the size of `\GF{4}`-codes of length 11 and minimal distance 3::

       sage: codes.bounds.delsarte_bound_hamming_space(11,3,4)
       327680/3

    An improvement of a known upper bound (150) from https://www.win.tue.nl/~aeb/codes/binary-1.html ::

       sage: a,p,x = codes.bounds.delsarte_bound_hamming_space(23,10,2,return_data=True,isinteger=True); x # long time
       148
       sage: [j for i,j in p.get_values(a).items()]                                      # long time
       [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 95, 0, 2, 0, 36, 0, 14, 0, 0, 0, 0, 0, 0, 0]

    Note that a usual LP, without integer variables, won't do the trick ::

       sage: codes.bounds.delsarte_bound_hamming_space(23,10,2).n(20)
       151.86

    Such an input is invalid::

       sage: codes.bounds.delsarte_bound_hamming_space(11,3,-4)
       Solver exception: PPL : There is no feasible solution
       False
    """
    from sage.numerical.mip import MIPSolverException
    A, p = _delsarte_LP_building(n, d, 0, q, isinteger, solver)
    try:
        bd = p.solve()
    except MIPSolverException as exc:
        print(f"Solver exception: {exc}")
        return (A, p, False) if return_data else False

    return (A, p, bd) if return_data else bd


def delsarte_bound_additive_hamming_space(n, d, q, d_star=1, q_base=0, return_data=False,
                                          solver='PPL', isinteger=False):
    r"""
    Find a modified Delsarte bound on additive codes in Hamming space `H_q^n` of minimal distance `d`.

    Find the Delsarte LP bound on ``F_{q_base}``-dimension of additive
    codes in Hamming space `H_q^n` of minimal distance ``d`` with
    minimal distance of the dual code at least ``d_star``.  If
    ``q_base`` is set to nonzero, then ``q`` is a power of
    ``q_base``, and the code is, formally, linear over
    ``F_{q_base}``. Otherwise it is assumed that ``q_base==q``.

    INPUT:

    - ``n`` -- the code length

    - ``d`` -- the (lower bound on) minimal distance of the code

    - ``q`` -- the size of the alphabet

    - ``d_star`` -- the (lower bound on) minimal distance of the dual code;
      only makes sense for additive codes

    - ``q_base`` -- if ``0``, the code is assumed to be linear. Otherwise,
      ``q=q_base^m`` and the code is linear over ``F_{q_base}``

    - ``return_data`` -- if ``True``, return a triple ``(W,LP,bound)``,
      where ``W`` is a weights vector, and ``LP`` the Delsarte bound
      LP; both of them are Sage LP data.  ``W`` need not be a weight
      distribution of a code, or, if ``isinteger==False``, even have
      integer entries.

    - ``solver`` -- the LP/ILP solver to be used. Defaults to ``'PPL'``. It is
      arbitrary precision, thus there will be no rounding errors. With
      other solvers (see :class:`MixedIntegerLinearProgram` for the
      list), you are on your own!

    - ``isinteger`` -- if ``True``, uses an integer programming solver (ILP),
      rather that an LP solver (can be very slow if set to ``True``)

    EXAMPLES:

    The bound on dimension of linear `\GF{2}`-codes of length 11 and minimal distance 6::

        sage: codes.bounds.delsarte_bound_additive_hamming_space(11, 6, 2)
        3
        sage: a,p,val = codes.bounds.delsarte_bound_additive_hamming_space(\
                             11, 6, 2, return_data=True)
        sage: [j for i,j in p.get_values(a).items()]
        [1, 0, 0, 0, 0, 0, 5, 2, 0, 0, 0, 0]

    The bound on the dimension of linear `\GF{4}`-codes of length 11 and minimal distance 3::

        sage: codes.bounds.delsarte_bound_additive_hamming_space(11,3,4)
        8

    The bound on the `\GF{2}`-dimension of additive `\GF{4}`-codes of length 11 and minimal
    distance 3::

        sage: codes.bounds.delsarte_bound_additive_hamming_space(11,3,4,q_base=2)
        16

    Such a ``d_star`` is not possible::

        sage: codes.bounds.delsarte_bound_additive_hamming_space(11,3,4,d_star=9)
        Solver exception: PPL : There is no feasible solution
        False

    TESTS::

        sage: a,p,x = codes.bounds.delsarte_bound_additive_hamming_space(\
                         19,15,7,return_data=True,isinteger=True)
        sage: [j for i,j in p.get_values(a).items()]
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 307, 0, 0, 1, 34]
        sage: codes.bounds.delsarte_bound_additive_hamming_space(19,15,7,solver='glpk')
        3
        sage: codes.bounds.delsarte_bound_additive_hamming_space(\
                 19,15,7, isinteger=True, solver='glpk')
        3
    """
    from sage.numerical.mip import MIPSolverException
    if q_base == 0:
        q_base = q

    kk = 0
    while q_base**kk < q:
        kk += 1

    if q_base**kk != q:
        print("Wrong q_base=", q_base, " for q=", q, kk)
        return False

    # this implementation assumes that our LP solver to be unable to do a hot
    # restart with an adjusted constraint

    m = kk*n  # this is to emulate repeat/until block
    bd = q**n+1

    while q_base**m < bd:
        # need to solve the LP repeatedly, as this is a new constraint!
        # we might become infeasible. More precisely, after rounding down
        # to the closest value of q_base^m, the LP, with the constraint that
        # the objective function is at most q_base^m,
        A, p = _delsarte_LP_building(n, d, d_star, q, isinteger,
                                     solver, q_base**m)
        try:
            bd = p.solve()
        except MIPSolverException as exc:
            print("Solver exception:", exc)
            return (A, p, False) if return_data else False
    # rounding the bound down to the nearest power of q_base, for q=q_base^m
    # bd_r = roundres(log(bd, base=q_base))
        m = -1
        while q_base**(m+1) < bd:
            m += 1
        if q_base**(m+1) == bd:
            m += 1

    return (A, p, m) if return_data else m


def _delsarte_Q_LP_building(q, d, solver, isinteger):
    r"""
    LP builder for Delsarte's LP for codes given Q matrix.

    LP builder for Delsarte's LP for codes, given Q matrix;
    used in :func:`delsarte_bound_Q_matrix`; not exported.

    INPUT:

    - ``q`` -- the Q matrix

    - ``d`` -- the (lower bound on) minimal distance of the code

    - ``solver`` -- the LP/ILP solver to be used. Defaults to
      ``PPL``. It is arbitrary precision, thus there will be no
      rounding errors. With other solvers (see
      :class:`MixedIntegerLinearProgram` for the list), you are on
      your own!

    - ``isinteger`` -- if ``True``, uses an integer programming solver
      (ILP), rather that an LP solver. Can be very slow if set to
      ``True``.

    EXAMPLES::

        sage: from sage.coding.delsarte_bounds import _delsarte_Q_LP_building
        sage: q = Matrix([[codes.bounds.krawtchouk(6,2,i,j) for j in range(7)] for i in range(7)])
        sage: _, p = _delsarte_Q_LP_building(q, 2, "PPL", False)
        sage: p.show()
        Maximization:
          x_0 + x_1 + x_2 + x_3 + x_4 + x_5 + x_6
        <BLANKLINE>
        Constraints:
          constraint_0: 1 <= x_0 <= 1
          constraint_1: 0 <= x_1 <= 0
          constraint_2: 0 <= 6 x_0 + 4 x_1 + 2 x_2 - 2 x_4 - 4 x_5 - 6 x_6
          constraint_3: 0 <= 15 x_0 + 5 x_1 - x_2 - 3 x_3 - x_4 + 5 x_5 + 15 x_6
          constraint_4: 0 <= 20 x_0 - 4 x_2 + 4 x_4 - 20 x_6
          constraint_5: 0 <= 15 x_0 - 5 x_1 - x_2 + 3 x_3 - x_4 - 5 x_5 + 15 x_6
          constraint_6: 0 <= 6 x_0 - 4 x_1 + 2 x_2 - 2 x_4 + 4 x_5 - 6 x_6
          constraint_7: 0 <= x_0 - x_1 + x_2 - x_3 + x_4 - x_5 + x_6
        Variables:
          x_0 is a continuous variable (min=0, max=+oo)
          x_1 is a continuous variable (min=0, max=+oo)
          x_2 is a continuous variable (min=0, max=+oo)
          x_3 is a continuous variable (min=0, max=+oo)
          x_4 is a continuous variable (min=0, max=+oo)
          x_5 is a continuous variable (min=0, max=+oo)
          x_6 is a continuous variable (min=0, max=+oo)
    """
    from sage.numerical.mip import MixedIntegerLinearProgram

    n, _ = q.dimensions()  # Q is a square matrix

    p = MixedIntegerLinearProgram(maximization=True, solver=solver)
    A = p.new_variable(integer=isinteger, nonnegative=True)
    p.set_objective(sum([A[i] for i in range(n)]))

    p.add_constraint(A[0] == 1)

    try:
        for i in range(1, d):
            p.add_constraint(A[i] == 0)
    except TypeError:
        for i in d:
            p.add_constraint(A[i] == 0)

    for k in range(1, n):
        p.add_constraint(sum([q[k][i] * A[i] for i in range(n)]), min=0)

    return A, p


def delsarte_bound_Q_matrix(q, d, return_data=False,
                            solver='PPL', isinteger=False):
    r"""
    Delsarte bound on a code with Q matrix ``q`` and lower bound on min. dist. ``d``.

    Find the Delsarte bound on a code with Q matrix ``q`` and lower bound on
    minimal distance ``d``.

    INPUT:

    - ``q`` -- the Q matrix

    - ``d`` -- the (lower bound on) minimal distance of the code

    - ``return_data`` -- if ``True``, return a triple
      ``(W,LP,bound)``, where ``W`` is a weights vector, and ``LP``
      the Delsarte upper bound LP; both of them are Sage LP data.
      ``W`` need not be a weight distribution of a code.

    - ``solver`` -- the LP/ILP solver to be used. Defaults to
      ``PPL``. It is arbitrary precision, thus there will be no
      rounding errors. With other solvers (see
      :class:`MixedIntegerLinearProgram` for the list), you are on
      your own!

    - ``isinteger`` -- if ``True``, uses an integer programming solver
      (ILP), rather that an LP solver. Can be very slow if set to
      ``True``.

    EXAMPLES:

    The bound on dimension of linear `\GF{2}`-codes of length 10 and minimal distance 6::

        sage: q_matrix = Matrix([[codes.bounds.krawtchouk(10,2,i,j) for i in range(11)]
        ....:                    for j in range(11)])
        sage: codes.bounds.delsarte_bound_Q_matrix(q_matrix, 6)
        2

        sage: a,p,val = codes.bounds.delsarte_bound_Q_matrix(q_matrix, 6, return_data=True)
        sage: [j for i,j in p.get_values(a).items()]
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

    TESTS:

    Cases for using Hamming scheme Q matrix::

        sage: q_matrix = Matrix([[codes.bounds.krawtchouk(10,2,i,j) for i in range(11)] for j in range(11)])
        sage: codes.bounds.delsarte_bound_Q_matrix(q_matrix, 6)
        2

        sage: a,p,val = codes.bounds.delsarte_bound_Q_matrix(q_matrix, 6, return_data=True)
        sage: [j for i,j in p.get_values(a).items()]
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    """
    from sage.numerical.mip import MIPSolverException
    from sage.structure.element import Matrix

    if not isinstance(q, Matrix):
        raise ValueError("Input to delsarte_bound_Q_matrix "
                         "should be a sage Matrix()")

    A, p = _delsarte_Q_LP_building(q, d, solver, isinteger)
    try:
        bd = p.solve()
    except MIPSolverException as exc:
        print(f"Solver exception: {exc}")
        return (A, p, False) if return_data else False

    return (A, p, bd) if return_data else bd
