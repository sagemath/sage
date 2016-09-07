r"""
Joint Spectral Radius

This modules contains algorithms in conjunction with the
:wikipedia:`joint spectral radius <Joint_spectral_radius>`.


Various
=======

REFERENCES:

.. [Dumas2013] Philippe Dumas,
   *Joint spectral radius, dilation equations, and asymptotic behavior
   of radix-rational sequences*,
   Linear Algebra and its Applications 438 (2013), 2107-2126,
   :doi:`10.1016/j.laa.2012.10.013`.

.. [Gripenberg1996] Gustaf Gripenberg,
   *Computing the joint spectral radius*,
   Linear Algebra and its Applications 234 (1996), 43-60,
   :doi:`10.1016/0024-3795(94)00082-4`.

AUTHORS:

- Daniel Krenn (2016)

ACKNOWLEDGEMENT:

- Daniel Krenn is supported by the
  Austrian Science Fund (FWF): P 24644-N26.


Functions
=========
"""
#*****************************************************************************
#       Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


def joint_spectral_radius(S, delta=None, norm=None,
                          ring=None, max_iterations=None):
    r"""
    Return a lower and upper bound for the joint spectral radius
    of the given matrices.

    INPUT:

    - ``S`` -- an tuple or other iterable of matrices.

    - ``delta`` -- (default: ``0.1``) a real specifying the size
      of the returned interval.

    - ``norm`` -- (default: the `2`-norm) a function mapping a matrix ``M``
      to ``norm(M)``

    - ``ring`` -- (default: :mod:`RIF <sage.rings.real_mpfi>`) a
      real interval field.

    - ``max_iterations`` -- (default: ``1000``) a positive integer.

    OUTPUT:

    A :mod:`real interval field element <sage.rings.real_mpfi>`.

    ALGORITHM:

    This function implements the algorithm presented in [Gripenberg1996]_.

    EXAMPLES::

        sage: import logging
        sage: logging.basicConfig()
        sage: logger = logging.getLogger('sage.matrix.joint_spectral_radius')
        sage: logger.setLevel(logging.INFO)

    Example of [Gripenberg1996]_, Section 4 (result is between `0.6596789`
    and `0.6596924`)::

        sage: joint_spectral_radius(  # long time
        ....:     (Matrix([[3, 0],  [1, 3]]) / 5,
        ....:      Matrix([[3, -3], [0, -1]]) / 5),
        ....:     delta=RIF(0.0001))
        INFO:...:lower bound: 0.659678908955284?
        INFO:...:upper bound: 0.659778908955284?
        INFO:...:iterations m: 48
        INFO:...:max size of T_m: 21
        0.6597?

        sage: logger.setLevel(logging.DEBUG)
        sage: joint_spectral_radius(  # long time
        ....:     (Matrix([[3, 0],  [1, 3]]) / 5,
        ....:      Matrix([[3, -3], [0, -1]]) / 5),
        ....:     delta=RIF(0.01))
        DEBUG:...:m=2, alpha=0.6000000000000000?, beta=0.70827625302982189?, len(T)=4
        DEBUG:...:m=3, alpha=0.6000000000000000?, beta=0.706734263510838?, len(T)=6
        DEBUG:...:m=4, alpha=0.6000000000000000?, beta=0.704390998023197?, len(T)=8
        ...
        DEBUG:...:m=16, alpha=0.659678908955284?, beta=0.6707696905811751?, len(T)=1
        DEBUG:...:m=17, alpha=0.659678908955284?, beta=0.669678908955284?, len(T)=0
        INFO:...:lower bound: 0.659678908955284?
        INFO:...:upper bound: 0.669678908955284?
        INFO:...:iterations m: 17
        INFO:...:max size of T_m: 15
        0.66?

    ::

        sage: logger.setLevel(logging.WARNING)

    [Dumas2013]_, Example 3 (result is `1`)::

        sage: A0 = Matrix([[1, 1, 1], [0, 0, 0],  [0, 0, 0]])
        sage: A1 = Matrix([[1, 0, 0], [0, 1, -1], [0, 0, 0]])
        sage: A2 = Matrix([[0, 0, 0], [1, 1, 0],  [0, 0, 1]])
        sage: A3 = Matrix([[0, 0, 0], [0, 0, 0],  [1, -1, 1]])
        sage: joint_spectral_radius((A0, A1, A2, A3),  # long time
        ....:     delta=RIF(0.1)).str(style='brackets')
        '[1.0000000000000000 .. 1.1000000000000001]'

    [Dumas2013]_, Example 4 (result is `1`)::

        sage: A0 = Matrix([[1, 1/2, 0], [0, 1/2, 0], [0, 1/2, 1]])
        sage: A1 = Matrix([[1/2, 0, 0], [1/2, 1, 0], [1/2, 0, 1]])
        sage: joint_spectral_radius((A0, A1), delta=RIF(0.2)).str(style='brackets')
        '[1.0000000000000000 .. 1.2000000000000002]'

    [Dumas2013]_, Example 5 (result is `3\sqrt{2}=4.24\dots`)::

        sage: A0 = Matrix([[1, 0], [0, 1]])
        sage: A1 = Matrix([[3, -3], [3, 3]])
        sage: joint_spectral_radius((A0, A1), delta=RIF(0.2)).str(style='brackets')
        '[4.2426406871192847 .. 4.2426406871192848]'

    [Dumas2013]_, Example 6 (result is `1`)::

        sage: B0 = Matrix([[1/2, 0], [1/2, 1]])
        sage: B1 = Matrix([[1/2, 0], [-1/2, 1]])
        sage: Z = zero_matrix(2, 2)
        sage: A0 = block_matrix([[B1, Z], [Z, B0]])
        sage: A1 = block_matrix([[Z, Z], [B0, B1]])
        sage: joint_spectral_radius((A0, A1), delta=RIF(0.2)).str(style='brackets')
        '[1.0000000000000000 .. 1.2000000000000002]'

    TESTS::

        sage: A0 = Matrix([[1, 1, 1], [0, 0, 0],  [0, 0, 0]])
        sage: A1 = Matrix([[1, 0, 0], [0, 1, -1], [0, 0, 0]])
        sage: A2 = Matrix([[0, 0, 0], [1, 1, 0],  [0, 0, 1]])
        sage: A3 = Matrix([[0, 0, 0], [0, 0, 0],  [1, -1, 1]])
        sage: joint_spectral_radius((A0, A1, A2, A3),
        ....:     delta=RIF(0.1), max_iterations=5).str(style='brackets')
        WARNING:sage.matrix.joint_spectral_radius:Stopping since
        maximum number of iterations reached.
        '[1.0000000000000000 .. 1.1472026904398774]'
    """
    from sage.arith.srange import srange
    from sage.misc.misc_c import prod

    import logging
    logger = logging.getLogger(__name__)

    S = tuple(S)
    if ring is None:
        from sage.rings.real_mpfi import RIF
        ring = RIF
    R = ring
    if delta is None:
        delta = R(0.1)
    if norm is None:
        norm = lambda M: M.norm(2)
    Rnorm = lambda M: R(norm(M))
    if max_iterations is None:
        max_iterations = 1000

    def rho(M):
        return max(R(abs(v)) for v in M.eigenvalues())

    def pp(X, j):
        return Rnorm(prod(X[:j]))**(1/R(j))

    def p(X):
        return min(pp(X, j) for j in srange(1, len(X)+1))

    T = tuple(((M,), Rnorm(M)) for M in S)
    alpha = max(rho(M) for M in S)
    beta = max(pX for X, pX in T)  # pX equals Rnorm(M) here
    ell = 0

    for m in srange(2, max_iterations):

        prepreT = tuple((X + (M,), pX)
                        for X, pX in T for M in S)
        preT = tuple((X, min(pX, pp(X, len(X)))) for X, pX in prepreT)
        T = tuple((X, pX) for X, pX in preT if pX > alpha + delta)

        ell = max(ell, len(T))
        if T:
            alpha = max(alpha,
                        max(rho(prod(Y))**(1/R(m)) for Y, pY in T))
            beta = min(beta,
                       max(alpha + delta,
                           max(pY for Y, pY in T)))
        else:
            beta = min(beta, alpha + delta)

        logger.debug('m=%s, alpha=%s, beta=%s, len(T)=%s',
                     m, alpha, beta, len(T))

        if not T:
            break

    else:
        logger.warn('Stopping since maximum number of iterations reached.')

    logger.info('lower bound: %s', alpha)
    logger.info('upper bound: %s', beta)
    logger.info('iterations m: %s', m)
    logger.info('max size of T_m: %s', ell)

    return R(alpha, beta)
