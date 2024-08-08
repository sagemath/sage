r"""
Interface to LattE integrale programs
"""
# ****************************************************************************
#       Copyright (C) 2017 Vincent Delecroix <vincent.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.cpython.string import str_to_bytes, bytes_to_str

from subprocess import Popen, PIPE
from sage.rings.integer import Integer
from sage.features.latte import Latte_count, Latte_integrate


def count(arg, ehrhart_polynomial=False, multivariate_generating_function=False, raw_output=False, verbose=False, **kwds):
    r"""
    Call to the program count from LattE integrale.

    INPUT:

    - ``arg`` -- a cdd or LattE description string

    - ``ehrhart_polynomial``, ``multivariate_generating_function`` -- to
      compute Ehrhart polynomial or multivariate generating function instead of
      just counting points

    - ``raw_output`` -- if ``True`` then return directly the output string from LattE

    - For all other options of the count program, consult the LattE manual

    OUTPUT:

    Either a string (if ``raw_output`` if set to ``True``) or an integer (when
    counting points), or a polynomial (if ``ehrhart_polynomial`` is set to
    ``True``) or a multivariate THING (if ``multivariate_generating_function``
    is set to ``True``)

    EXAMPLES::

        sage: from sage.interfaces.latte import count    # optional - latte_int
        sage: P = 2 * polytopes.cube()

    Counting integer points from either the H or V representation::

        sage: count(P.cdd_Hrepresentation(), cdd=True)   # optional - latte_int
        125
        sage: count(P.cdd_Vrepresentation(), cdd=True)   # optional - latte_int
        125

    Ehrhart polynomial::

        sage: count(P.cdd_Hrepresentation(), cdd=True, ehrhart_polynomial=True)  # optional - latte_int
        64*t^3 + 48*t^2 + 12*t + 1

    Multivariate generating function currently only work with ``raw_output=True``::

        sage: opts = {'cdd': True,
        ....:         'multivariate_generating_function': True,
        ....:         'raw_output': True}
        sage: cddin = P.cdd_Hrepresentation()
        sage: print(count(cddin, **opts))  # optional - latte_int
        x[0]^2*x[1]^(-2)*x[2]^(-2)/((1-x[1])*(1-x[2])*(1-x[0]^(-1)))
         + x[0]^(-2)*x[1]^(-2)*x[2]^(-2)/((1-x[1])*(1-x[2])*(1-x[0]))
         + x[0]^2*x[1]^(-2)*x[2]^2/((1-x[1])*(1-x[2]^(-1))*(1-x[0]^(-1)))
         + x[0]^(-2)*x[1]^(-2)*x[2]^2/((1-x[1])*(1-x[0])*(1-x[2]^(-1)))
         + x[0]^2*x[1]^2*x[2]^(-2)/((1-x[2])*(1-x[1]^(-1))*(1-x[0]^(-1)))
         + x[0]^(-2)*x[1]^2*x[2]^(-2)/((1-x[2])*(1-x[0])*(1-x[1]^(-1)))
         + x[0]^2*x[1]^2*x[2]^2/((1-x[2]^(-1))*(1-x[1]^(-1))*(1-x[0]^(-1)))
         + x[0]^(-2)*x[1]^2*x[2]^2/((1-x[0])*(1-x[2]^(-1))*(1-x[1]^(-1)))

    TESTS:

    Testing raw output::

        sage: from sage.interfaces.latte import count   # optional - latte_int
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: count(cddin, cdd=True, raw_output=True)  # optional - latte_int
        '19'
        sage: count(cddin, cdd=True, raw_output=True, ehrhart_polynomial=True) # optional - latte_int
        ' + 1 * t^0 + 10/3 * t^1 + 8 * t^2 + 20/3 * t^3'
        sage: count(cddin, cdd=True, raw_output=True, multivariate_generating_function=True) # optional - latte_int
        'x[0]^(-1)*x[1]^(-1)/((1-x[0]*x[2])*(1-x[0]^(-1)*x[1])*...x[0]^(-1)*x[2]^(-1)))\n'

    Testing the ``verbose`` option::

        sage: n = count(cddin, cdd=True, verbose=True, raw_output=True)  # optional - latte_int
        This is LattE integrale ...
        ...
        Invocation: ...count '--redundancy-check=none' --cdd /dev/stdin
        ...
        Total Unimodular Cones: ...
        Maximum number of simplicial cones in memory at once: ...
        <BLANKLINE>
        ****  The number of lattice points is:   ****
        Total time: ... sec

    Trivial input for which LattE's preprocessor does all the work::

        sage: P = Polyhedron(vertices=[[0,0,0]])
        sage: cddin = P.cdd_Hrepresentation()
        sage: count(cddin, cdd=True, raw_output=False)  # optional - latte_int
        1

    Testing the runtime error::

        sage: P = Polyhedron(rays=[[0,1], [1,0]])
        sage: cddin = P.cdd_Hrepresentation()
        sage: count(cddin, cdd=True, raw_output=False)  # optional - latte_int
        Traceback (most recent call last):
        ...
        RuntimeError: LattE integrale program failed (exit code 1):
        This is LattE integrale ...
        ...
        The polyhedron is unbounded.
    """
    arg = str_to_bytes(arg)

    args = [Latte_count().absolute_filename()]
    if ehrhart_polynomial and multivariate_generating_function:
        raise ValueError
    if ehrhart_polynomial:
        args.append('--ehrhart-polynomial')
    elif multivariate_generating_function:
        args.append('--multivariate-generating-function')

    if 'redundancy_check' not in kwds:
        args.append('--redundancy-check=none')

    for key, value in kwds.items():
        if value is None or value is False:
            continue

        key = key.replace('_', '-')
        if value is True:
            args.append(f'--{key}')
        else:
            args.append(f'--{key}={value}')

    if multivariate_generating_function:
        from sage.misc.temporary_file import tmp_filename
        filename = tmp_filename()
        with open(filename, 'w') as f:
            f.write(bytes_to_str(arg))
        args += [filename]
    else:
        args += ['/dev/stdin']

    # The cwd argument is needed because latte
    # always produces diagnostic output files.
    import tempfile
    tempd = tempfile.TemporaryDirectory()

    latte_proc = Popen(args,
                       stdin=PIPE, stdout=PIPE,
                       stderr=(None if verbose else PIPE),
                       cwd=tempd.name)

    ans, err = latte_proc.communicate(arg)
    if err:
        err = bytes_to_str(err)
    ret_code = latte_proc.poll()
    if ret_code:
        if err is None:
            err = ", see error message above"
        else:
            err = ":\n" + err
        raise RuntimeError("LattE integrale program failed (exit code {})".format(ret_code) + err.strip())

    ans = bytes_to_str(ans)

    # There's an error handler below that uses the numOfLatticePoints
    # file created by latte, so we can't cleanup() the temporary
    # directory here. Instead we have to clean it up before the
    # (several) return statements.
    if ehrhart_polynomial:
        ans = ans.splitlines()[-2]
        if raw_output:
            tempd.cleanup()
            return ans
        else:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            from sage.rings.rational_field import QQ
            R = PolynomialRing(QQ, 't')
            tempd.cleanup()
            return R(ans)
    elif multivariate_generating_function:
        with open(filename + '.rat') as f:
            ans = f.read()
        if raw_output:
            tempd.cleanup()
            return ans
        else:
            raise NotImplementedError("there is no Sage object to handle multivariate series from LattE, use raw_output=True")
    else:
        if ans:  # Sometimes (when LattE's preproc does the work), no output appears on stdout.
            ans = ans.splitlines()[-1]
        if not ans:
            # opening a file is slow (30e-6s), so we read the file
            # numOfLatticePoints only in case of a IndexError above
            with open(tempd.name + '/numOfLatticePoints', 'r') as f:
                ans = f.read()

        if raw_output:
            tempd.cleanup()
            return ans
        else:
            tempd.cleanup()
            return Integer(ans)


def integrate(arg, polynomial=None, algorithm='triangulate', raw_output=False, verbose=False, **kwds):
    r"""
    Call to the function integrate from LattE integrale.

    INPUT:

    - ``arg`` -- a cdd or LattE description string

    - ``polynomial`` -- multivariate polynomial or valid LattE polynomial description string
      If given, the valuation parameter of LattE is set to integrate, and is set to volume otherwise

    - ``algorithm`` -- (default: ``'triangulate'``) the integration method; use 'triangulate' for
      polytope triangulation or 'cone-decompose' for tangent cone decomposition method

    - ``raw_output`` -- if ``True`` then return directly the output string from LattE

    - ``verbose`` -- if ``True`` then return directly verbose output from LattE

    - For all other options of the integrate program, consult the LattE manual

    OUTPUT: either a string (if ``raw_output`` if set to ``True``) or a rational

    EXAMPLES::

        sage: from sage.interfaces.latte import integrate   # optional - latte_int
        sage: P = 2 * polytopes.cube()
        sage: x, y, z = polygen(QQ, 'x, y, z')

    Integrating over a polynomial over a polytope in either the H or V representation::

        sage: integrate(P.cdd_Hrepresentation(), x^2*y^2*z^2, cdd=True)   # optional - latte_int
        4096/27
        sage: integrate(P.cdd_Vrepresentation(), x^2*y^2*z^2, cdd=True)   # optional - latte_int
        4096/27

    Computing the volume of a polytope in either the H or V representation::

        sage: integrate(P.cdd_Hrepresentation(), cdd=True)   # optional - latte_int
        64
        sage: integrate(P.cdd_Vrepresentation(), cdd=True)   # optional - latte_int
        64

    Polynomials given as a string in LattE description are also accepted::

        sage: integrate(P.cdd_Hrepresentation(), '[[1,[2,2,2]]]', cdd=True)   # optional - latte_int
        4096/27

    TESTS:

    Testing raw output::

        sage: from sage.interfaces.latte import integrate   # optional - latte_int
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: x, y, z = polygen(QQ, 'x, y, z')
        sage: f = 3*x^2*y^4*z^6 + 7*y^3*z^5
        sage: integrate(cddin, f, cdd=True, raw_output=True)  # optional - latte_int
        '629/47775'

    Testing the ``verbose`` option to integrate over a polytope::

        sage: ans = integrate(cddin, f, cdd=True, verbose=True, raw_output=True)  # optional - latte_int
        This is LattE integrale ...
        ...
        Invocation: ...integrate --valuation=integrate --triangulate --redundancy-check=none --cdd --monomials=... /dev/stdin
        ...

    Testing triangulate algorithm::

        sage: from sage.interfaces.latte import integrate   # optional - latte_int
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: integrate(cddin, algorithm='triangulate', cdd=True)  # optional - latte_int
        20/3

    Testing convex decomposition algorithm::

        sage: from sage.interfaces.latte import integrate   # optional - latte_int
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: integrate(cddin, algorithm='cone-decompose', cdd=True)  # optional - latte_int
        20/3

    Testing raw output::

        sage: from sage.interfaces.latte import integrate   # optional - latte_int
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: integrate(cddin, cdd=True, raw_output=True)  # optional - latte_int
        '20/3'

    Testing polynomial given as a string in LattE description::

        sage: from sage.interfaces.latte import integrate   # optional - latte_int
        sage: P = polytopes.cuboctahedron()
        sage: integrate(P.cdd_Hrepresentation(), '[[3,[2,4,6]],[7,[0, 3, 5]]]', cdd=True)   # optional - latte_int
        629/47775

    Testing the ``verbose`` option to compute the volume of a polytope::

        sage: from sage.interfaces.latte import integrate   # optional - latte_int
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: ans = integrate(cddin, cdd=True, raw_output=True, verbose=True)  # optional - latte_int
        This is LattE integrale ...
        ...
        Invocation: ...integrate --valuation=volume --triangulate --redundancy-check=none --cdd /dev/stdin
        ...

    Testing the runtime error::

        sage: P = Polyhedron(rays=[[1,0],[0,1]])
        sage: P._volume_latte()  # optional - latte_int
        Traceback (most recent call last):
        ...
        RuntimeError: LattE integrale program failed (exit code -6):
        This is LattE integrale ...
        ...
        determinant: nonsquare matrix
    """
    arg = str_to_bytes(arg)

    from sage.rings.rational import Rational

    args = [Latte_integrate().absolute_filename()]

    got_polynomial = bool(polynomial is not None)

    if got_polynomial:
        args.append('--valuation=integrate')
    else:
        args.append('--valuation=volume')

    if algorithm == 'triangulate':
        args.append('--triangulate')
    elif algorithm == 'cone-decompose':
        args.append('--cone-decompose')

    if 'redundancy_check' not in kwds:
        args.append('--redundancy-check=none')

    for key, value in kwds.items():
        if value is None or value is False:
            continue

        key = key.replace('_', '-')
        if value is True:
            args.append(f'--{key}')
        else:
            args.append(f'--{key}={value}')

    if got_polynomial:
        if not isinstance(polynomial, str):
            # transform polynomial to LattE description
            monomials_list = to_latte_polynomial(polynomial)
        else:
            monomials_list = str(polynomial)

        from sage.misc.temporary_file import tmp_filename
        filename_polynomial = tmp_filename()

        with open(filename_polynomial, 'w') as f:
            f.write(monomials_list)
            args += ['--monomials=' + filename_polynomial]

    args += ['/dev/stdin']

    # The cwd argument is needed because latte
    # always produces diagnostic output files.
    import tempfile
    tempd = tempfile.TemporaryDirectory()

    latte_proc = Popen(args,
                       stdin=PIPE, stdout=PIPE,
                       stderr=(None if verbose else PIPE),
                       cwd=tempd.name)

    ans, err = latte_proc.communicate(arg)
    if err:
        err = bytes_to_str(err)
    ret_code = latte_proc.poll()
    if ret_code:
        if err is None:
            err = ", see error message above"
        else:
            err = ":\n" + err
        raise RuntimeError("LattE integrale program failed (exit code {})".format(ret_code) + err.strip())

    ans = bytes_to_str(ans)
    ans = ans.splitlines()
    ans = ans[-5].split()
    assert ans[0] == 'Answer:'
    ans = ans[1]

    tempd.cleanup()
    if raw_output:
        return ans
    else:
        return Rational(ans)


def to_latte_polynomial(polynomial):
    r"""
    Helper function to transform a polynomial to its LattE description.

    INPUT:

    - ``polynomial`` -- a multivariate polynomial

    OUTPUT: string that describes the monomials list and exponent vectors

    TESTS:

    Testing a polynomial in three variables::

        sage: from sage.interfaces.latte import to_latte_polynomial
        sage: x, y, z = polygens(QQ, 'x, y, z')
        sage: f = 3*x^2*y^4*z^6 + 7*y^3*z^5
        sage: to_latte_polynomial(f)
        '[[3, [2, 4, 6]], [7, [0, 3, 5]]]'

        sage: to_latte_polynomial(x.parent().zero())
        '[]'

    Testing a univariate polynomial::

        sage: x = polygen(QQ, 'x')
        sage: to_latte_polynomial((x-1)^2)
        '[[1, [0]], [-2, [1]], [1, [2]]]'

        sage: to_latte_polynomial(x.parent().zero())
        '[]'
    """
    if polynomial == 0:
        return str([])

    from sage.rings.polynomial.polydict import ETuple

    coefficients_list = polynomial.coefficients()

    # transform list of exponents into a list of lists.
    # this branch handles the multivariate/univariate case
    if isinstance(polynomial.exponents()[0], ETuple):
        exponents_list = [list(exponent_vector_i) for exponent_vector_i in polynomial.exponents()]
    else:
        exponents_list = [[exponent_vector_i] for exponent_vector_i in polynomial.exponents()]

    # assuming that the order in coefficients() and exponents() methods match
    monomials_list = [list(monomial_i)
                      for monomial_i
                      in zip(coefficients_list, exponents_list)]

    return str(monomials_list)
