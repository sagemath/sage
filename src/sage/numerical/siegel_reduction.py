r"""
Siegel reduction of period and Riemann matrices
===============================================

We consider `g\times 2g` complex matrices partitioned as `(\Omega_1 | \Omega_2)`
such that the associated `g\times g` *Riemann matrix* `\Omega=\Omega_1^{-1}\Omega_2`
satisfies:

* `\Omega` is symmetric,

* the imaginary part of `\Omega` is positive-definite.

This property is preserved under the right-action of `Sp(2g,\ZZ)`.
The notion of being *Siegel reduced* for Riemann matrices (see [DHBvHS2004]_)
can be extended to period matrices by defining a period matrix to be Siegel reduced
if the associated Riemann matrix is.

This module implements a routine to compute a Siegel reduced form, together
with the transformation matrix.

EXAMPLES::

    sage: from sage.numerical.siegel_reduction import siegel_reduction
    sage: from sage.schemes.riemann_surfaces.riemann_surface import numerical_inverse
    sage: CC = ComplexField(20)
    sage: P = matrix(CC,2,4,[1,3,1+5*I,12+10*I,0,1,I,4+3*I])
    sage: Phat, Gamma = siegel_reduction(P)
    sage: Phat
    [  1.0000   3.0000 5.0000*I 10.000*I]
    [ 0.00000   1.0000 1.0000*I 3.0000*I]
    sage: Gamma
    [ 1  0 -1  0]
    [ 0  1  0 -4]
    [ 0  0  1  0]
    [ 0  0  0  1]
    sage: numerical_inverse(Phat[:,:2])*Phat[:,2:]
    [2.0000*I 1.0000*I]
    [1.0000*I 3.0000*I]

We can also pass in a Riemann matrix::

    sage: Omega = numerical_inverse(P[:,:2])*P[:,2:]
    sage: Omega_hat , Gamma2 = siegel_reduction(Omega)
    sage: Phat[:,:2]^(-1)*Phat[:,2:] == Omega_hat
    True
    sage: Gamma == Gamma2
    True

AUTHORS:

 - Nils Bruin, Sohrab Ganjian (2021-09-08): initial version
"""

# ****************************************************************************
#       Copyright (C) 2021 Nils Bruin <nbruin@sfu.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix, block_matrix
from sage.schemes.riemann_surfaces.riemann_surface import numerical_inverse
from sage.numerical.riemann_theta import (
    cholesky_decomposition,
    imag_func,
    real_func,
    round_func,
)
from sage.libs.pari import pari


def _siegel_big_period_matrix(big_omega):
    """
    Return a Siegel reduced matrix and the transformation matrix gamma of a period matrix of a Riemann surface.

    INPUT:

    - ``big_omega`` -- gx2g Riemann period matrix of a Riemann surface

    OUTPUT: a tuple of

    - ``big_omega_hat`` -- gx2g Siegel reduced Riemann period matrix

    - ``gamma_matrix`` -- 2gx2g transformation matrix such that big_omega_hat = big_omega * gamma_matrix

    EXAMPLES:

    An example from a genus 2 curve::

        sage: from sage.numerical.siegel_reduction import _siegel_big_period_matrix
        sage: R.<X,Y>=QQ[]
        sage: C = Curve(Y^2-(X^6+X+1))
        sage: RS = C.riemann_surface()
        sage: PM = RS.period_matrix()
        sage: M, G = _siegel_big_period_matrix(PM)

    An example from a genus 5 curve::

        sage: R.<X,Y>=QQ[]
        sage: C = Curve(Y^2-(X^10+3))
        sage: RS = C.riemann_surface()
        sage: PM = RS.period_matrix()
        sage: M, G = _siegel_big_period_matrix(PM)

    REFERENCES:

    .. [DHBvHS2004]_
    .. [AC2019]_

    AUTHORS:

    - Nils Bruin, Sohrab Ganjian (2021-08-19): initial verision
    """
    big_omega_hat = big_omega
    CC = big_omega_hat.base_ring()
    half = CC(0.5)
    g = big_omega.nrows()
    I_g = matrix.identity(g)
    gamma_matrix = matrix.identity(2 * g)
    zero_matrix = matrix.zero(g)

    A = diagonal_matrix([0] + [1 for i in range(g - 1)])
    B = diagonal_matrix([-1] + [0 for i in range(g - 1)])
    AB_block_matrix = block_matrix(2, 2, [A, B, -B, A])

    while True:
        M_siegel = numerical_inverse(big_omega_hat[:, :g]) * big_omega_hat[:, g:]
        M_siegel = half * (M_siegel + M_siegel.transpose())

        T = cholesky_decomposition(M_siegel.apply_map(imag_func))
        U = pari.qflll(T).sage()
        TU = T * U
        TU_norm = [v.norm() for v in TU.columns()]
        min_index = TU_norm.index(min(TU_norm))

        if min_index != 0:
            temp = U[:, min_index]
            U[:, min_index] = U[:, 0]
            U[:, 0] = temp

        M_siegel = U.transpose() * M_siegel * U
        gamma_matrix = gamma_matrix * block_matrix(
            2, 2, [U.transpose().inverse(), zero_matrix, zero_matrix, U]
        )
        X = M_siegel.apply_map(real_func)
        x = X.apply_map(round_func)

        M_siegel = M_siegel - x
        gamma_matrix = gamma_matrix * block_matrix(2, 2, [I_g, -x, zero_matrix, I_g])

        if 1 <= M_siegel[0, 0].abs():
            big_omega_hat = big_omega * gamma_matrix
            break

        gamma_matrix = gamma_matrix * AB_block_matrix
        big_omega_hat = big_omega * gamma_matrix

    return big_omega_hat, gamma_matrix


def siegel_reduction(M):
    """
    Return a Siegel reduced matrix, together with the transformation matrix.

    INPUT:

    The input can be either a gxg Riemman matrix or a gx2g period matrix of a Riemann surface.

    - ``M`` -- gxg Riemann matrix or gx2g period matrix

    OUTPUT:

    The outputs are matrices omega_hat and gamma_matrix. Depedning on the size of the input, omega_hat
    can either be a gxg Siegel reduced Riemann Matrix or a gx2g Siegel reduced period matrix of a Riemann Surface.
    The former happens when the input matrix is gxg, and the latter occurs when the input is gx2g.

    - ``omega_hat`` -- gxg Siegel reduced Riemann matrix or gx2g Siegel reduced period matrix of a Riemann Surface

    - ``gamma_matrix`` -- 2gx2g transformation matrix

    EXAMPLES::

        sage: from sage.numerical.siegel_reduction import siegel_reduction
        sage: omega = (-1/(2*CC.pi()*CC.gen())) * Matrix(CC, [[111.207, 96.616],[96.616, 83.943]])
        sage: M, G = siegel_reduction(omega)

    An example from a genus 5 curve::

        sage: R.<X,Y>=QQ[]
        sage: C = Curve(Y^2-(X^10+3))
        sage: RS = C.riemann_surface()
        sage: RM = RS.riemann_matrix()
        sage: M, G = siegel_reduction(RM)

    REFERENCES:

    .. [DHBvHS2004]_
    .. [AC2019]_

    AUTHORS:

    - Nils Bruin, Sohrab Ganjian (2021-08-19): initial verision
    """
    g = M.nrows()
    n = M.ncols()

    if g == n:
        g = M.nrows()
        CC = M.base_ring()
        big_omega = matrix.identity(CC, g).augment(M)
        big_omega_hat, gamma_matrix = _siegel_big_period_matrix(big_omega)
        omega = numerical_inverse(big_omega_hat[:, :g]) * big_omega_hat[:, g:]
        return omega, gamma_matrix

    if 2 * g == n:
        return _siegel_big_period_matrix(M)

    raise ValueError(
        "Input matrix should either be a gxg Riemann matrix or a gx2g big period matrix"
    )
