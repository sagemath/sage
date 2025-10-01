r"""
Library of cythonized methods

AUTHORS:

- Harald Schilly (2011-01-16): initial version (#9623) partially based on work by Lauri Ruotsalainen
"""
# ***************************************************************************
#       Copyright (C) 2011 Harald Schilly <harald.schilly@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************


cpdef julia(ff_j, z, int iterations):
    """
    Helper function for the Julia Fractal interact example.

    INPUT:

    - ``ff_j`` -- fast callable for the inner iteration
    - ``z`` -- complex number
    - ``iterations`` -- number of loops

    TESTS::

        sage: # needs sage.symbolic
        sage: from sage.interacts.library_cython import julia
        sage: z = var('z')
        sage: c_real, c_imag = 1, 1
        sage: f = symbolic_expression(z**2 + c_real + c_imag * CDF.gen()).function(z)
        sage: ff_m = fast_callable(f, vars=[z], domain=CDF)
        sage: julia(ff_m, CDF(1,1), 3)
        1.0 + 3.0*I
    """
    for i in range(iterations):
        z = ff_j(z)
        if z.abs() > 2:
            break
    return z


cpdef mandel(ff_m, z, int iterations):
    """
    Helper function for the Mandelbrot Fractal interact example.

    INPUT:

    - ``ff_m`` -- fast callable for the inner iteration
    - ``z`` -- complex number
    - ``iterations`` -- number of loops

    TESTS::

        sage: # needs sage.symbolic
        sage: from sage.interacts.library_cython import mandel
        sage: z, c = var('z, c')
        sage: f = symbolic_expression(z**2 + c).function(z,c)
        sage: ff_m = fast_callable(f, vars=[z,c], domain=CDF)
        sage: mandel(ff_m, CDF(1,1), 3)
        1.0 + 3.0*I
    """
    c = z
    for i in range(iterations):
        z = ff_m(z, c)
        if z.abs() > 2:
            break
    return z


cpdef cellular(rule, int N):
    """
    Cythonized helper function for the cellular_automata fractal.

    Yields a matrix showing the evolution of a Wolfram's cellular automaton.
    Based on work by Pablo Angulo.
    http://wiki.sagemath.org/interact/misc#CellularAutomata

    INPUT:

    - ``rule`` -- determines how a cell's value is updated, depending
      on its neighbors
    - ``N`` -- number of iterations

    TESTS::

        sage: from sage.interacts.library_cython import cellular
        sage: rule = [1, 0, 1, 0, 0, 1, 1, 0]
        sage: N = 3
        sage: print(cellular(rule, N))
        [[0 0 0 1 0 0 0 0]
         [1 1 0 1 0 1 0 0]
         [0 1 1 1 1 1 0 0]]
    """
    from numpy import zeros
    cdef int j, k, l
    M = zeros((N, 2*N+2), dtype=int)
    M[0, N] = 1

    for j in range(1, N):
        for k in range(0, 2*N):
            l = 4 * M[j-1, k-1] + 2 * M[j-1, k] + M[j-1, k+1]
            M[j, k] = rule[l]
    return M
