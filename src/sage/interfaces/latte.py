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
    Call to the program count from LattE integrale

    INPUT:

    - ``arg`` -- a cdd or LattE description string

    - ``ehrhart_polynomial``, ``multivariate_generating_function``  -- to
      compute Ehrhart polynomial or multivariate generating function instead of
      just counting points

    - ``raw_output`` -- if ``True`` then return directly the output string from LattE

    - For all other options of the count program, consult the LattE manual

    - ``name`` -- -- (default: ``'y'``) a string

      The variable names of the Laurent polynomial ring of the multivariate_generating_function

    - ``Factorization_sort`` (default: ``False``) and
      ``Factorization_simplify`` (default: ``False``) -- booleans

      These are passed on to
      :class:`sage.structure.factorization.Factorization` when creating
      the result.

    OUTPUT:

    Either a string (if ``raw_output`` if set to ``True``) or an integer (when
    counting points), or a polynomial (if ``ehrhart_polynomial`` is set to
    ``True``) or a tuple of ``Factorization`` objects
    :class:`~sage.structure.factorization.Factorization` whose factors are Laurent polynomials
    (if ``multivariate_generating_function`` is set to ``True``)

    EXAMPLES::

        sage: from sage.interfaces.latte import count
        sage: P = 2 * polytopes.cube()

    Counting integer points from either the H or V representation::

        sage: count(P.cdd_Hrepresentation(), cdd=True)      # optional - latte_int
        125
        sage: count(P.cdd_Vrepresentation(), cdd=True)      # optional - latte_int
        125

    Ehrhart polynomial::

        sage: count(P.cdd_Hrepresentation(), cdd=True,      # optional - latte_int
        ....:       ehrhart_polynomial=True)
        64*t^3 + 48*t^2 + 12*t + 1

    Returning a string of the multivariate generating function when ``raw_output=True``.
    Returning the summands of the multivariate generating function in a tuple of ``Factorization`` objects
    with the same format as
    :meth:`sage.geometry.polyhedron.generating_function.generating_function_of_integral_points`
    does when ``result_as_tuple=True``::

        sage: opts = {'cdd': True,
        ....:         'multivariate_generating_function': True,
        ....:         'raw_output': True}
        sage: cddin = P.cdd_Hrepresentation()
        sage: print(count(cddin, **opts))                                       # optional - latte_int
        x[0]^2*x[1]^(-2)*x[2]^(-2)/((1-x[1])*(1-x[2])*(1-x[0]^(-1)))
         + x[0]^(-2)*x[1]^(-2)*x[2]^(-2)/((1-x[1])*(1-x[2])*(1-x[0]))
         + x[0]^2*x[1]^(-2)*x[2]^2/((1-x[1])*(1-x[2]^(-1))*(1-x[0]^(-1)))
         + x[0]^(-2)*x[1]^(-2)*x[2]^2/((1-x[1])*(1-x[0])*(1-x[2]^(-1)))
         + x[0]^2*x[1]^2*x[2]^(-2)/((1-x[2])*(1-x[1]^(-1))*(1-x[0]^(-1)))
         + x[0]^(-2)*x[1]^2*x[2]^(-2)/((1-x[2])*(1-x[0])*(1-x[1]^(-1)))
         + x[0]^2*x[1]^2*x[2]^2/((1-x[2]^(-1))*(1-x[1]^(-1))*(1-x[0]^(-1)))
         + x[0]^(-2)*x[1]^2*x[2]^2/((1-x[0])*(1-x[2]^(-1))*(1-x[1]^(-1)))
        sage: count(cddin, cdd=True, multivariate_generating_function=True)     # optional - latte_int
        ((y0^2*y1^-2*y2^-2) * (-y1 + 1)^-1 * (-y2 + 1)^-1 * (1 - y0^-1)^-1,
         (y0^-2*y1^-2*y2^-2) * (-y1 + 1)^-1 * (-y2 + 1)^-1 * (-y0 + 1)^-1,
         (y0^2*y1^-2*y2^2) * (-y1 + 1)^-1 * (1 - y2^-1)^-1 * (1 - y0^-1)^-1,
         (y0^-2*y1^-2*y2^2) * (-y1 + 1)^-1 * (-y0 + 1)^-1 * (1 - y2^-1)^-1,
         (y0^2*y1^2*y2^-2) * (-y2 + 1)^-1 * (1 - y1^-1)^-1 * (1 - y0^-1)^-1,
         (y0^-2*y1^2*y2^-2) * (-y2 + 1)^-1 * (-y0 + 1)^-1 * (1 - y1^-1)^-1,
         y0^2*y1^2*y2^2 * (1 - y2^-1)^-1 * (1 - y1^-1)^-1 * (1 - y0^-1)^-1,
         (y0^-2*y1^2*y2^2) * (-y0 + 1)^-1 * (1 - y2^-1)^-1 * (1 - y1^-1)^-1)

    TESTS:

    Testing raw output::

        sage: from sage.interfaces.latte import count
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: count(cddin, cdd=True, raw_output=True)                           # optional - latte_int
        '19'
        sage: count(cddin, cdd=True, raw_output=True, ehrhart_polynomial=True)  # optional - latte_int
        ' + 1 * t^0 + 10/3 * t^1 + 8 * t^2 + 20/3 * t^3'
        sage: count(cddin, cdd=True, raw_output=True,                           # optional - latte_int
        ....:       multivariate_generating_function=True)
        'x[0]^(-1)*x[1]^(-1)/((1-x[0]*x[2])*(1-x[0]^(-1)*x[1])*(1-x[2]^(-1)))\n + x[0]^(-1)*x[1]^(-1)/((1-x[2])*(1-x[0]^(-1)*x[1])*(1-x[0]*x[2]^(-1)))\n + ... + x[0]*x[1]/((1-x[0]^(-1)*x[2])*(1-x[0]*x[1]^(-1))*(1-x[2]^(-1)))\n + x[0]*x[1]/((1-x[2])*(1-x[0]*x[1]^(-1))*(1-x[0]^(-1)*x[2]^(-1)))\n'

    Testing multivariate generating function::

        sage: from sage.interfaces.latte import count                           # optional - latte_int
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: count(cddin, cdd=True, multivariate_generating_function=True)     # optional - latte_int
        ((y0^-1*y1^-1) * (-y0*y2 + 1)^-1 * (1 - y0^-1*y1)^-1 * (1 - y2^-1)^-1,
         (y0^-1*y1^-1) * (-y2 + 1)^-1 * (1 - y0^-1*y1)^-1 * (-y0*y2^-1 + 1)^-1,
        ...
         y0*y1 * (1 - y0^-1*y2)^-1 * (-y0*y1^-1 + 1)^-1 * (1 - y2^-1)^-1,
         y0*y1 * (-y2 + 1)^-1 * (-y0*y1^-1 + 1)^-1 * (1 - y0^-1*y2^-1)^-1)

        sage: P = Polyhedron(rays=[[0,1], [1,0]])
        sage: cddin = P.cdd_Hrepresentation()
        sage: count(cddin, cdd=True, raw_output=True,                           # optional - latte_int
        ....:       multivariate_generating_function=True)
        '1/((1-x[1])*(1-x[0]))\n'
        sage: count(cddin, cdd=True, multivariate_generating_function=True)     # optional - latte_int
        (1 * (-y1 + 1)^-1 * (-y0 + 1)^-1,)

    Testing the ``verbose`` option::

        sage: from sage.interfaces.latte import count
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: n = count(cddin, cdd=True, verbose=True, raw_output=True)         # optional - latte_int
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
        sage: count(cddin, cdd=True, raw_output=False)                          # optional - latte_int
        1

    Testing the runtime error::

        sage: P = Polyhedron(rays=[[0,1], [1,0]])
        sage: cddin = P.cdd_Hrepresentation()
        sage: count(cddin, cdd=True, raw_output=False)                          # optional - latte_int
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

    mgf_kwds = {}
    for key,value in kwds.items():
        if key in ["name", "Factorization_sort", "Factorization_simplify", "sort_factors"]:
            mgf_kwds[key] = value
            continue
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
            tempd.cleanup()
            return str_to_multivariate_generating_function(ans, **mgf_kwds)
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


def str_to_multivariate_generating_function(raw_output_str, name=None, **kwds):
    r"""
    Helper function for :func:`count` if ``multivariate_generating_function`` is set to ``True``
    which preprocess the raw output string to a tuple of summands.

    TESTS:

        sage: from sage.interfaces.latte import count, str_to_multivariate_generating_function
        sage: P = Polyhedron(ieqs=[(0, 1, 0, 0), (0, -1, 1, 0)], eqns=[(0, -1, -1, 2)])
        sage: cddin = P.cdd_Hrepresentation()
        sage: raw_output_str = count(cddin, cdd=True, raw_output=True,
        ....:                        multivariate_generating_function=True)
        sage: str_to_multivariate_generating_function(raw_output_str, name='xi')
        (1 * (-xi0*xi1*xi2 + 1)^-1 * (-xi1^2*xi2 + 1)^-1,)
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    import re

    indices_regex = re.compile(r'(?<=\[)(\d*)(?=\])')
    indices = range(max(ZZ(index) for index in indices_regex.findall(raw_output_str)) + 1)
    if name is None:
        name = 'y'

    B = LaurentPolynomialRing(ZZ,
                        tuple(name + str(k) for k in indices),
                        len(indices))
    raw_output_list = raw_output_str[:-1].split('\n + ')
    return tuple(_str_to_multivariate_generating_function(a, B, **kwds) for a in raw_output_list)


def _str_to_multivariate_generating_function(summand, B=None,
                                             Factorization_sort=False, Factorization_simplify=False,
                                             sort_factors=False):
    r"""
    Helper function for :func:`str_to_multivariate_generating_function`
    which convert a summand string to a ``Factorization`` object.

    Each summand is in the format
    .. MATH::

        (\sum_{k\in K}y^{s_0}*y^{v_k}})\prod_{j\in J}(1 - y^{v_j})^{-1}.

    TESTS::

        sage: from sage.interfaces.latte import _str_to_multivariate_generating_function
        sage: B = LaurentPolynomialRing(ZZ, 'y', 3)
        sage: _str_to_multivariate_generating_function(
        ....:     '(-1)*x[0]^(-1)*x[2]/((1-x[0]^(-1)*x[1]^(-1))*(1-x[0]^(-1)*x[2]^(-1))*(1-x[0]))', B)
        (-y0^-1*y2) * (1 - y0^-1*y1^-1)^-1 * (1 - y0^-1*y2^-1)^-1 * (-y0 + 1)^-1
        sage: _str_to_multivariate_generating_function('(-1)/((1-x[0]*x[1]*x[2])*(1-x[1]^2*x[2]))\n', B)
        (-1) * (-y0*y1*y2 + 1)^-1 * (-y1^2*y2 + 1)^-1
        sage: _str_to_multivariate_generating_function(
        ....:    '((-1)*x[0]*x[2]^2 + x[1]^(-2)*x[2])/((1-x[0]*x[1]*x[2])*(1-x[1]^2*x[2]))', B)
        (-y0*y2^2 + y1^-2*y2) * (-y0*y1*y2 + 1)^-1 * (-y1^2*y2 + 1)^-1
    """
    from sage.rings.integer_ring import ZZ
    from sage.structure.factorization import Factorization
    import re

    numerator_str, denominator_str = summand.split('/')

    gen_regex = re.compile(r'(?<=\[)(\d*)(?=\])')
    exponent_regex = re.compile(r'([\d|-]+)')

    def str_to_laurent_monomial(monomial_str, B):
        result = 1
        for gen_str in monomial_str.split('*'):
            gen_exponent = gen_str.split('^')
            if len(gen_exponent) == 1:
                result *= B.gens()[ZZ(gen_regex.findall(gen_exponent[0])[0])]
            else:
                result *= B.gens()[ZZ(gen_regex.findall(gen_exponent[0])[0])] ** ZZ(exponent_regex.findall(gen_exponent[1])[0])
        return result

    def str_to_coef_times_laurent_monomial(monomial_str, B):
        if 'x' not in monomial_str:
            if '*' in monomial_str:
                aa, bb = monomial_str.split('*')
                return ZZ(aa.replace('(','').replace(')',''))*ZZ(bb.replace('(','').replace(')',''))
            else:
                return ZZ(monomial_str.replace('(','').replace(')',''))
        elif 'x' in monomial_str.split('*',1)[0]:
            return str_to_laurent_monomial(monomial_str, B)
        else:
            return ZZ(monomial_str.split('*',1)[0].replace('(','').replace(')',''))*str_to_laurent_monomial(monomial_str.split('*',1)[1],B)

    numerator = sum(str_to_coef_times_laurent_monomial(a, B) for a in numerator_str.split('+'))

    term_regex = re.compile(r'(?<=1-)(.+?)(?=$|\)$|\)\*\()')
    terms = (str_to_laurent_monomial(a, B) for a in term_regex.findall(denominator_str[1:-1]))
    if sort_factors:
        def key(t):
            D = t.dict().popitem()[0]
            return (-sum(abs(d) for d in D), D)
        terms = sorted(terms, key=key, reverse=True)
    return Factorization([(numerator, 1)] +
                         [(1-t, -1) for t in terms],
                         sort=Factorization_sort,
                         simplify=Factorization_simplify)


def integrate(arg, polynomial=None, algorithm='triangulate', raw_output=False, verbose=False, **kwds):
    r"""
    Call to the function integrate from LattE integrale.

    INPUT:

    - ``arg`` -- a cdd or LattE description string.

    - ``polynomial`` -- multivariate polynomial or valid LattE polynomial description string.
      If given, the valuation parameter of LattE is set to integrate, and is set to volume otherwise.

    - ``algorithm`` -- (default: 'triangulate') the integration method. Use 'triangulate' for
      polytope triangulation or 'cone-decompose' for tangent cone decomposition method.

    - ``raw_output`` -- if ``True`` then return directly the output string from LattE.

    - ``verbose`` -- if ``True`` then return directly verbose output from LattE.

    - For all other options of the integrate program, consult the LattE manual.

    OUTPUT:

    Either a string (if ``raw_output`` if set to ``True``) or a rational.

    EXAMPLES::

        sage: from sage.interfaces.latte import integrate
        sage: P = 2 * polytopes.cube()
        sage: x, y, z = polygen(QQ, 'x, y, z')

    Integrating over a polynomial over a polytope in either the H or V representation::

        sage: integrate(P.cdd_Hrepresentation(), x^2*y^2*z^2, cdd=True)         # optional - latte_int
        4096/27
        sage: integrate(P.cdd_Vrepresentation(), x^2*y^2*z^2, cdd=True)         # optional - latte_int
        4096/27

    Computing the volume of a polytope in either the H or V representation::

        sage: integrate(P.cdd_Hrepresentation(), cdd=True)                      # optional - latte_int
        64
        sage: integrate(P.cdd_Vrepresentation(), cdd=True)                      # optional - latte_int
        64

    Polynomials given as a string in LattE description are also accepted::

        sage: integrate(P.cdd_Hrepresentation(), '[[1,[2,2,2]]]', cdd=True)     # optional - latte_int
        4096/27

    TESTS:

    Testing raw output::

        sage: from sage.interfaces.latte import integrate
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: x, y, z = polygen(QQ, 'x, y, z')
        sage: f = 3*x^2*y^4*z^6 + 7*y^3*z^5
        sage: integrate(cddin, f, cdd=True, raw_output=True)                    # optional - latte_int
        '629/47775'

    Testing the ``verbose`` option to integrate over a polytope::

        sage: ans = integrate(cddin, f, cdd=True, verbose=True,                 # optional - latte_int
        ....:                 raw_output=True)
        This is LattE integrale ...
        ...
        Invocation: ...integrate --valuation=integrate --triangulate --redundancy-check=none --cdd --monomials=... /dev/stdin
        ...

    Testing triangulate algorithm::

        sage: from sage.interfaces.latte import integrate
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: integrate(cddin, algorithm='triangulate', cdd=True)               # optional - latte_int
        20/3

    Testing convex decomposition algorithm::

        sage: from sage.interfaces.latte import integrate
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: integrate(cddin, algorithm='cone-decompose', cdd=True)            # optional - latte_int
        20/3

    Testing raw output::

        sage: from sage.interfaces.latte import integrate
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: integrate(cddin, cdd=True, raw_output=True)                       # optional - latte_int
        '20/3'

    Testing polynomial given as a string in LattE description::

        sage: from sage.interfaces.latte import integrate
        sage: P = polytopes.cuboctahedron()
        sage: integrate(P.cdd_Hrepresentation(),                                # optional - latte_int
        ....:           '[[3,[2,4,6]],[7,[0, 3, 5]]]', cdd=True)
        629/47775

    Testing the ``verbose`` option to compute the volume of a polytope::

        sage: from sage.interfaces.latte import integrate
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: ans = integrate(cddin, cdd=True, raw_output=True, verbose=True)   # optional - latte_int
        This is LattE integrale ...
        ...
        Invocation: ...integrate --valuation=volume --triangulate --redundancy-check=none --cdd /dev/stdin
        ...

    Testing the runtime error::

        sage: P = Polyhedron(rays=[[1,0],[0,1]])
        sage: P._volume_latte()                                                 # optional - latte_int
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

    - ``polynomial`` -- a multivariate polynomial.

    OUTPUT:

    A string that describes the monomials list and exponent vectors.

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
