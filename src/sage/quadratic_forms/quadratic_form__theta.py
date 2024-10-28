"""
Theta Series of Quadratic Forms

AUTHORS:

- Jonathan Hanke: initial code, theta series of degree 1

- Gonzalo Tornaria (2009-02-22): fixes and doctests

- Gonzalo Tornaria (2010-03-23): theta series of degree 2
"""
from copy import deepcopy

from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.integer_ring import ZZ


def theta_series(self, Max=10, var_str='q', safe_flag=True):
    r"""
    Compute the theta series as a power series in the variable given
    in ``var_str`` (which defaults to ``'q'``), up to the specified precision
    `O(q^{Max})`.

    This uses the PARI/GP function :pari:`qfrep`, wrapped by the
    theta_by_pari() method. This caches the result for future
    computations.

    The ``safe_flag`` allows us to select whether we want a copy of the
    output, or the original output.  It is only meaningful when a
    vector is returned, otherwise a copy is automatically made in
    creating the power series.  By default ``safe_flag`` = ``True``, so we
    return a copy of the cached information.  If this is set to ``False``,
    then the routine is much faster but the return values are
    vulnerable to being corrupted by the user.

    .. TODO::

        Allow the option ``Max='mod_form'`` to give enough coefficients
        to ensure we determine the theta series as a modular form.  This
        is related to the Sturm bound, but we will need to be careful about
        this (particularly for half-integral weights!).

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.theta_series()                                                          # needs sage.libs.pari
        1 + 2*q + 2*q^3 + 6*q^4 + 2*q^5 + 4*q^6 + 6*q^7 + 8*q^8 + 14*q^9 + O(q^10)

        sage: Q.theta_series(25)                                                        # needs sage.libs.pari
        1 + 2*q + 2*q^3 + 6*q^4 + 2*q^5 + 4*q^6 + 6*q^7 + 8*q^8 + 14*q^9 + 4*q^10
        + 12*q^11 + 18*q^12 + 12*q^13 + 12*q^14 + 8*q^15 + 34*q^16 + 12*q^17 + 8*q^18
        + 32*q^19 + 10*q^20 + 28*q^21 + 16*q^23 + 44*q^24 + O(q^25)
    """
    # Sanity Check: Max is an integer or an allowed string:
    try:
        M = ZZ(Max)
    except TypeError:
        M = -1

    if Max not in ['mod_form'] and not M >= 0:
        raise TypeError("Max = {Max} is not an integer >= 0 or an allowed string")

    if Max == 'mod_form':
        raise NotImplementedError("we have to figure out the correct number of Fourier coefficients to use")
        # return self.theta_by_pari(sturm_bound(self.level(), self.dim() / ZZ(2)) + 1, var_str, safe_flag)
    else:
        return self.theta_by_pari(M, var_str, safe_flag)


# ---- Compute the theta function by using the PARI/GP routine qfrep  ----

def theta_by_pari(self, Max, var_str='q', safe_flag=True):
    r"""
    Use PARI/GP to compute the theta function as a power series (or
    vector) up to the precision `O(q^{Max})`.  This also caches the result
    for future computations.

    If ``var_str`` = ``''``, then we return a vector `v` where ``v[i]`` counts the
    number of vectors of length `i`.

    The ``safe_flag`` allows us to select whether we want a copy of the
    output, or the original output.  It is only meaningful when a
    vector is returned, otherwise a copy is automatically made in
    creating the power series.  By default ``safe_flag=True``, so we
    return a copy of the cached information.  If this is set to ``False``,
    then the routine is much faster but the return values are
    vulnerable to being corrupted by the user.

    INPUT:

    - ``Max`` -- integer `\geq 0`
    - ``var_str`` -- string

    OUTPUT: a power series or a vector

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Prec = 100
        sage: compute = Q.theta_by_pari(Prec, '')                                       # needs sage.libs.pari
        sage: exact = [1] + [8 * sum([d  for d in divisors(i)  if d % 4 != 0])          # needs sage.libs.pari
        ....:                for i in range(1, Prec)]
        sage: compute == exact                                                          # needs sage.libs.pari
        True
    """
    # Try to use the cached result if it's enough precision
    if hasattr(self, '__theta_vec') and len(self.__theta_vec) >= Max:
        theta_vec = self.__theta_vec[:Max]
    else:
        theta_vec = self.representation_number_list(Max)
        # Cache the theta vector
        self.__theta_vec = theta_vec

    # Return the answer
    if not var_str:
        if safe_flag:
            return deepcopy(theta_vec)         # We must make a copy here to insure the integrity of the cached version!
        else:
            return theta_vec
    else:
        return PowerSeriesRing(ZZ, var_str)(theta_vec, Max)


# -- Compute the theta function by using an explicit Cholesky decomposition --


########################################################################
# Routines to compute the Fourier expansion of the theta function of Q #
# (to a given precision) via a Cholesky decomposition over RR.         #
#                                                                      #
# The Cholesky code was taken from:                                    #
# ~/Documents/290_Project/C/Ver13.2__3-5-2007/Matrix_mpz/Matrix_mpz.cc #
########################################################################


def theta_by_cholesky(self, q_prec):
    r"""
    Uses the real Cholesky decomposition to compute (the `q`-expansion of) the
    theta function of the quadratic form as a power series in `q` with terms
    correct up to the power `q^{\text{q\_prec}}`. (So its error is `O(q^
    {\text{q\_prec} + 1})`.)

    REFERENCE:

    Cohen's "A Course in Computational Algebraic Number Theory" book, p 102.

    EXAMPLES:

    Check the sum of 4 squares form against Jacobi's formula::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Theta = Q.theta_by_cholesky(10)
        sage: Theta
        1 + 8*q + 24*q^2 + 32*q^3 + 24*q^4 + 48*q^5 + 96*q^6
         + 64*q^7 + 24*q^8 + 104*q^9 + 144*q^10
        sage: Expected = [1] + [8*sum(d for d in divisors(n) if d%4)
        ....:                   for n in range(1, 11)]
        sage: Expected
        [1, 8, 24, 32, 24, 48, 96, 64, 24, 104, 144]
        sage: Theta.list() == Expected
        True

    Check the form `x^2 + 3y^2 + 5z^2 + 7w^2` represents everything except 2 and 22.::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Theta = Q.theta_by_cholesky(50)
        sage: Theta_list = Theta.list()
        sage: [m  for m in range(len(Theta_list))  if Theta_list[m] == 0]
        [2, 22]
    """
    # RAISE AN ERROR -- This routine is deprecated!
    # raise NotImplementedError("This routine is deprecated.  Try theta_series(), which uses theta_by_pari().")

    from sage.arith.misc import integer_ceil as ceil, integer_floor as floor
    from sage.misc.functional import sqrt
    from sage.rings.real_mpfr import RealField

    n = self.dim()
    theta = [0] * (q_prec + 1)
    PS = PowerSeriesRing(ZZ, 'q')

    bit_prec = 53                  # TO DO: Set this precision to reflect the appropriate roundoff
    Cholesky = self.cholesky_decomposition(bit_prec)     # error estimate, to be confident through our desired q-precision.
    Q = Cholesky      # <----  REDUNDANT!!!
    R = RealField(bit_prec)
    half = R(0.5)

    # 1. Initialize
    i = n - 1
    T = [R.zero()] * n
    U = [R.zero()] * n
    T[i] = R(q_prec)
    U[i] = 0
    L = [0] * n
    x = [0] * n

    # 2. Compute bounds
    # Z = sqrt(T[i] / Q[i,i])      # IMPORTANT NOTE: sqrt preserves the precision of the real number it's given... which is not so good... =|
    # L[i] = floor(Z - U[i])       # Note: This is a Sage Integer
    # x[i] = ceil(-Z - U[i]) - 1   # Note: This is a Sage Integer too

    done_flag = False
    from_step4_flag = False
    from_step3_flag = True        # We start by pretending this, since then we get to run through 2 and 3a once. =)

    # Big loop which runs through all vectors
    while not done_flag:

        # Loop through until we get to i=1 (so we defined a vector x)
        while from_step3_flag or from_step4_flag:              # IMPORTANT WARNING:  This replaces a do...while loop, so it may have to be adjusted!

            # Go to directly to step 3 if we're coming from step 4, otherwise perform step 2.
            if from_step4_flag:
                from_step4_flag = False
            else:
                # 2. Compute bounds
                from_step3_flag = False
                Z = sqrt(T[i] / Q[i, i])
                L[i] = floor(Z - U[i])
                x[i] = ceil(-Z - U[i]) - 1

            # 3a. Main loop
            x[i] += 1
            while x[i] > L[i]:
                i += 1
                x[i] += 1

            # 3b. Main loop
            if i > 0:
                from_step3_flag = True
                T[i - 1] = T[i] - Q[i, i] * (x[i] + U[i]) * (x[i] + U[i])
                i += -1
                U[i] = 0
                for j in range(i + 1, n):
                    U[i] += Q[i, j] * x[j]

        # 4. Solution found (This happens when i=0)
        from_step4_flag = True
        Q_val_double = q_prec - T[0] + Q[0, 0] * (x[0] + U[0]) * (x[0] + U[0])
        Q_val = floor(Q_val_double + half)        # Note: This rounds the value up, since the "round" function returns a float, but floor returns integer.

        # OPTIONAL SAFETY CHECK:
        eps = 0.000000001
        if abs(Q_val_double - Q_val) > eps:
            raise RuntimeError("Oh No! We have a problem with the floating point precision... \n"
                               + " Q_val_double = " + str(Q_val_double) + "\n"
                               + " Q_val = " + str(Q_val) + "\n"
                               + " x = " + str(x) + "\n")

        if Q_val <= q_prec:
            theta[Q_val] += 2

        # 5. Check if x = 0, for exit condition. =)
        done_flag = all(x[j] == 0 for j in range(n))

    # Set the value: theta[0] = 1
    theta[0] = 1

    # Return the series, truncated to the desired q-precision
    return PS(theta)


def theta_series_degree_2(Q, prec) -> dict:
    r"""
    Compute the theta series of degree 2 for the quadratic form `Q`.

    INPUT:

    - ``prec`` -- integer

    OUTPUT: dictionary, where:

    - keys are `{\rm GL}_2(\ZZ)`-reduced binary quadratic forms (given as triples of
      coefficients)
    - values are coefficients

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: Q2 = QuadraticForm(ZZ, 4, [1,1,1,1, 1,0,0, 1,0, 1])
        sage: S = Q2.theta_series_degree_2(10)
        sage: S[(0,0,2)]
        24
        sage: S[(1,0,1)]
        144
        sage: S[(1,1,1)]
        192

    AUTHORS:

    - Gonzalo Tornaria (2010-03-23)

    REFERENCE:

    - Raum, Ryan, Skoruppa, Tornaria, 'On Formal Siegel Modular Forms'
      (preprint)
    """
    from sage.arith.misc import integer_floor as floor
    from sage.misc.functional import sqrt
    from sage.misc.timing import cputime
    from sage.misc.verbose import verbose

    if Q.base_ring() != ZZ:
        raise TypeError("the quadratic form must be integral")
    if not Q.is_positive_definite():
        raise ValueError("the quadratic form must be positive definite")
    try:
        X = ZZ(prec - 1)    # maximum discriminant
    except TypeError:
        raise TypeError("prec is not an integer")

    if X < -1:
        raise ValueError("prec must be >= 0")

    if X == -1:
        return {}

    V = ZZ ** Q.dim()
    H = Q.Hessian_matrix()

    t = cputime()
    maxi = (X + 1) // 4
    v_list = (Q.vectors_by_length(maxi))        # assume a>0
    v_list = [[V(c) for c in vs] for vs in v_list]  # coerce vectors into V
    verbose("Computed vectors_by_length", t)

    # Deal with the singular part
    coeffs = {(0, 0, 0): ZZ.one()}
    for i in range(1, maxi + 1):
        coeffs[(0, 0, i)] = ZZ(2) * len(v_list[i])

    # Now deal with the non-singular part
    a_max = int(floor(sqrt(X / 3)))
    for a in range(1, a_max + 1):
        t = cputime()
        c_max = int(floor((a * a + X) / (4 * a)))
        for c in range(a, c_max + 1):
            for v1 in v_list[a]:
                v1_H = v1 * H

                def B_v1(v):
                    return v1_H * v2
                for v2 in v_list[c]:
                    b = abs(B_v1(v2))
                    if b <= a and 4 * a * c - b * b <= X:
                        qf = (a, b, c)
                        count = ZZ(4) if b == 0 else ZZ(2)
                        coeffs[qf] = coeffs.get(qf, ZZ.zero()) + count
        verbose("done a = %d" % a, t)

    return coeffs
