# sage_setup: distribution = sagemath-symbolics
"""
Symbolic sparse matrices

EXAMPLES::

    sage: matrix(SR, 2, 2, range(4), sparse=True)
    [0 1]
    [2 3]
    sage: matrix(SR, 2, 2, var('t'), sparse=True)
    [t 0]
    [0 t]

Arithmetic::

    sage: -matrix(SR, 2, range(4), sparse=True)
    [ 0 -1]
    [-2 -3]
    sage: m = matrix(SR, 2, [1..4], sparse=True); sqrt(2)*m
    [  sqrt(2) 2*sqrt(2)]
    [3*sqrt(2) 4*sqrt(2)]
    sage: m = matrix(SR, 4, [1..4^2], sparse=True)
    sage: m * m
    [ 90 100 110 120]
    [202 228 254 280]
    [314 356 398 440]
    [426 484 542 600]

    sage: m = matrix(SR, 3, [1, 2, 3], sparse=True); m
    [1]
    [2]
    [3]
    sage: m.transpose() * m
    [14]

Computing inverses::

    sage: M = matrix(SR, 2, var('a,b,c,d'), sparse=True)
    sage: ~M
    [1/a - b*c/(a^2*(b*c/a - d))           b/(a*(b*c/a - d))]
    [          c/(a*(b*c/a - d))              -1/(b*c/a - d)]
    sage: (~M*M).simplify_rational()
    [1 0]
    [0 1]
    sage: M = matrix(SR, 3, 3, range(9), sparse=True) - var('t')
    sage: (~M * M).simplify_rational()
    [1 0 0]
    [0 1 0]
    [0 0 1]

    sage: matrix(SR, 1, 1, 1, sparse=True).inverse()
    [1]
    sage: matrix(SR, 0, 0, sparse=True).inverse()
    []
    sage: matrix(SR, 3, 0, sparse=True).inverse()
    Traceback (most recent call last):
    ...
    ArithmeticError: self must be a square matrix

Transposition::

    sage: m = matrix(SR, 2, [sqrt(2), -1, pi, e^2], sparse=True)
    sage: m.transpose()
    [sqrt(2)      pi]
    [     -1     e^2]

``.T`` is a convenient shortcut for the transpose::

    sage: m.T
    [sqrt(2)      pi]
    [     -1     e^2]

Test pickling::

    sage: m = matrix(SR, 2, [sqrt(2), 3, pi, e], sparse=True); m
    [sqrt(2)       3]
    [     pi       e]
    sage: TestSuite(m).run()

Comparison::

    sage: m = matrix(SR, 2, [sqrt(2), 3, pi, e], sparse=True)
    sage: m == m
    True
    sage: m != 3
    True
    sage: m = matrix(SR,2,[1..4], sparse=True); n = m^2
    sage: (exp(m+n) - exp(m)*exp(n)).simplify_rational() == 0       # indirect test
    True


Determinant::

    sage: M = matrix(SR, 2, 2, [x,2,3,4], sparse=True)
    sage: M.determinant()
    4*x - 6
    sage: M = matrix(SR, 3,3,range(9), sparse=True)
    sage: M.det()
    0
    sage: t = var('t')
    sage: M = matrix(SR, 2, 2, [cos(t), sin(t), -sin(t), cos(t)], sparse=True)
    sage: M.det()
    cos(t)^2 + sin(t)^2
    sage: M = matrix([[sqrt(x),0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]], sparse=True)
    sage: det(M)
    sqrt(x)

Permanents::

    sage: M = matrix(SR, 2, 2, [x,2,3,4], sparse=True)
    sage: M.permanent()
    4*x + 6

Rank::

    sage: M = matrix(SR, 5, 5, range(25), sparse=True)
    sage: M.rank()
    2
    sage: M = matrix(SR, 5, 5, range(25), sparse=True) - var('t')
    sage: M.rank()
    5

    .. warning::

        :meth:`rank` may return the wrong answer if it cannot determine that a
        matrix element that is equivalent to zero is indeed so.

Copying symbolic matrices::

    sage: m = matrix(SR, 2, [sqrt(2), 3, pi, e], sparse=True)
    sage: n = copy(m)
    sage: n[0,0] = sin(1)
    sage: m
    [sqrt(2)       3]
    [     pi       e]
    sage: n
    [sin(1)      3]
    [    pi      e]

Conversion to Maxima::

    sage: m = matrix(SR, 2, [sqrt(2), 3, pi, e], sparse=True)
    sage: m._maxima_()
    matrix([sqrt(2),3],[%pi,%e])

TESTS:

Check that :issue:`12778` is fixed::

    sage: M = Matrix([[1, 0.9, 1/5, x^2], [2, 1.9, 2/5, x^3], [3, 2.9, 3/5, x^4]], sparse=True); M
    [                1 0.900000000000000               1/5               x^2]
    [                2  1.90000000000000               2/5               x^3]
    [                3  2.90000000000000               3/5               x^4]
    sage: parent(M)
    Full MatrixSpace of 3 by 4 sparse matrices over Symbolic Ring

Check that :issue:`35653` is fixed::

    sage: diagonal_matrix([x]).inverse()
    [1/x]
    sage: M = MatrixSpace(SR,2,2,sparse=True)
    sage: M([[x,0],[0,x]]).inverse()
    [1/x   0]
    [  0 1/x]
"""
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.factorization import Factorization

from sage.matrix.matrix_generic_sparse cimport Matrix_generic_sparse
from sage.matrix.constructor import matrix

cdef maxima

from sage.calculus.calculus import maxima

cdef class Matrix_symbolic_sparse(Matrix_generic_sparse):
    def echelonize(self, **kwds):
        """
        Echelonize using the classical algorithm.

        TESTS::

            sage: m = matrix([[cos(pi/5), sin(pi/5)], [-sin(pi/5), cos(pi/5)]], sparse=True)
            sage: m.echelonize(); m
            [1 0]
            [0 1]
        """
        return super().echelonize(algorithm='classical', **kwds)

    def eigenvalues(self, extend=True):
        """
        Compute the eigenvalues by solving the characteristic
        polynomial in maxima.

        The argument ``extend`` is ignored but kept for compatibility with
        other matrix classes.

        EXAMPLES::

            sage: a=matrix(SR,[[1,2],[3,4]], sparse=True)
            sage: a.eigenvalues()
            [-1/2*sqrt(33) + 5/2, 1/2*sqrt(33) + 5/2]

        TESTS:

        Check for :issue:`31700`::

            sage: m = matrix([[cos(pi/5), sin(pi/5)], [-sin(pi/5), cos(pi/5)]], sparse=True)
            sage: t = linear_transformation(m)
            sage: t.eigenvalues()
            [1/4*sqrt(5) - 1/4*sqrt(2*sqrt(5) - 10) + 1/4,
             1/4*sqrt(5) + 1/4*sqrt(2*sqrt(5) - 10) + 1/4]
        """
        maxima_evals = self._maxima_(maxima).eigenvalues()._sage_()
        if not len(maxima_evals):
            raise ArithmeticError("could not determine eigenvalues exactly using symbolic matrices; try using a different type of matrix via self.change_ring(), if possible")
        return sum([[ev] * int(mult) for ev, mult in zip(*maxima_evals)], [])

    def eigenvectors_left(self, other=None):
        r"""
        Compute the left eigenvectors of a matrix.

        INPUT:

        - ``other`` -- a square matrix `B` (default: ``None``) in a generalized
          eigenvalue problem; if ``None``, an ordinary eigenvalue problem is
          solved (currently supported only if the base ring of ``self`` is
          ``RDF`` or ``CDF``)

        OUTPUT:

        For each distinct eigenvalue, returns a list of the form (e,V,n)
        where e is the eigenvalue, V is a list of eigenvectors forming a
        basis for the corresponding left eigenspace, and n is the
        algebraic multiplicity of the eigenvalue.

        EXAMPLES::

            sage: A = matrix(SR,3,3,range(9), sparse=True); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: es = A.eigenvectors_left(); es
            [(-3*sqrt(6) + 6, [(1, -1/5*sqrt(6) + 4/5, -2/5*sqrt(6) + 3/5)], 1),
             (3*sqrt(6) + 6, [(1, 1/5*sqrt(6) + 4/5, 2/5*sqrt(6) + 3/5)], 1),
             (0, [(1, -2, 1)], 1)]
            sage: eval, [evec], mult = es[0]
            sage: delta = eval*evec - evec*A
            sage: abs(abs(delta)) < 1e-10
            3/5*sqrt(((2*sqrt(6) - 3)*(sqrt(6) - 2) + 7*sqrt(6) - 18)^2 + ((sqrt(6) - 2)*(sqrt(6) - 4) + 6*sqrt(6) - 14)^2) < (1.00000000000000e-10)
            sage: abs(abs(delta)).n() < 1e-10
            True

        ::

            sage: A = matrix(SR, 2, 2, var('a,b,c,d'), sparse=True)
            sage: A.eigenvectors_left()
            [(1/2*a + 1/2*d - 1/2*sqrt(a^2 + 4*b*c - 2*a*d + d^2), [(1, -1/2*(a - d + sqrt(a^2 + 4*b*c - 2*a*d + d^2))/c)], 1), (1/2*a + 1/2*d + 1/2*sqrt(a^2 + 4*b*c - 2*a*d + d^2), [(1, -1/2*(a - d - sqrt(a^2 + 4*b*c - 2*a*d + d^2))/c)], 1)]
            sage: es = A.eigenvectors_left(); es
            [(1/2*a + 1/2*d - 1/2*sqrt(a^2 + 4*b*c - 2*a*d + d^2), [(1, -1/2*(a - d + sqrt(a^2 + 4*b*c - 2*a*d + d^2))/c)], 1), (1/2*a + 1/2*d + 1/2*sqrt(a^2 + 4*b*c - 2*a*d + d^2), [(1, -1/2*(a - d - sqrt(a^2 + 4*b*c - 2*a*d + d^2))/c)], 1)]
            sage: eval, [evec], mult = es[0]
            sage: delta = eval*evec - evec*A
            sage: delta.apply_map(lambda x: x.full_simplify())
            (0, 0)

        This routine calls Maxima and can struggle with even small matrices
        with a few variables, such as a `3\times 3` matrix with three variables.
        However, if the entries are integers or rationals it can produce exact
        values in a reasonable time.  These examples create 0-1 matrices from
        the adjacency matrices of graphs and illustrate how the format and type
        of the results differ when the base ring changes.  First for matrices
        over the rational numbers, then the same matrix but viewed as a symbolic
        matrix. ::

            sage: G=graphs.CycleGraph(5)
            sage: am = G.adjacency_matrix(sparse=True)
            sage: spectrum = am.eigenvectors_left()
            sage: qqbar_evalue = spectrum[2][0]
            sage: type(qqbar_evalue)
            <class 'sage.rings.qqbar.AlgebraicNumber'>
            sage: qqbar_evalue
            0.618033988749895?

            sage: am = G.adjacency_matrix(sparse=True).change_ring(SR)
            sage: spectrum = am.eigenvectors_left()
            sage: symbolic_evalue = spectrum[2][0]
            sage: type(symbolic_evalue)
            <class 'sage.symbolic.expression.Expression'>
            sage: symbolic_evalue
            1/2*sqrt(5) - 1/2

            sage: bool(qqbar_evalue == symbolic_evalue)
            True

        A slightly larger matrix with a "nice" spectrum. ::

            sage: G = graphs.CycleGraph(6)
            sage: am = G.adjacency_matrix(sparse=True).change_ring(SR)
            sage: am.eigenvectors_left()
            [(-1, [(1, 0, -1, 1, 0, -1), (0, 1, -1, 0, 1, -1)], 2), (1, [(1, 0, -1, -1, 0, 1), (0, 1, 1, 0, -1, -1)], 2), (-2, [(1, -1, 1, -1, 1, -1)], 1), (2, [(1, 1, 1, 1, 1, 1)], 1)]

        TESTS::

            sage: A = matrix(SR, [[1, 2], [3, 4]], sparse=True)
            sage: B = matrix(SR, [[1, 1], [0, 1]], sparse=True)
            sage: A.eigenvectors_left(B)
            Traceback (most recent call last):
            ...
            NotImplementedError: generalized eigenvector decomposition is
            implemented for RDF and CDF, but not for Symbolic Ring

        Check that :issue:`23332` is fixed::

            sage: matrix([[x, x^2], [1, 0]], sparse=True).eigenvectors_left()
            [(-1/2*x*(sqrt(5) - 1), [(1, -1/2*x*(sqrt(5) + 1))], 1),
             (1/2*x*(sqrt(5) + 1), [(1, 1/2*x*(sqrt(5) - 1))], 1)]
        """
        if other is not None:
            raise NotImplementedError('generalized eigenvector decomposition '
                                      'is implemented for RDF and CDF, but '
                                      'not for %s' % self.base_ring())

        from sage.modules.free_module_element import vector
        from sage.rings.integer_ring import ZZ

        [evals, mults], evecs = self.transpose()._maxima_(maxima).eigenvectors()._sage_()
        result = []
        for e, evec, m in zip(evals, evecs, mults):
            result.append((e, [vector(v) for v in evec], ZZ(m)))

        return result

    def eigenvectors_right(self, other=None):
        r"""
        Compute the right eigenvectors of a matrix.

        INPUT:

        - ``other`` -- a square matrix `B` (default: ``None``) in a generalized
          eigenvalue problem; if ``None``, an ordinary eigenvalue problem is
          solved (currently supported only if the base ring of ``self`` is
          ``RDF`` or ``CDF``)

        OUTPUT:

        For each distinct eigenvalue, returns a list of the form (e,V,n)
        where e is the eigenvalue, V is a list of eigenvectors forming a
        basis for the corresponding right eigenspace, and n is the
        algebraic multiplicity of the eigenvalue.

        EXAMPLES::

            sage: A = matrix(SR,2,2,range(4), sparse=True); A
            [0 1]
            [2 3]
            sage: right = A.eigenvectors_right(); right
            [(-1/2*sqrt(17) + 3/2, [(1, -1/2*sqrt(17) + 3/2)], 1), (1/2*sqrt(17) + 3/2, [(1, 1/2*sqrt(17) + 3/2)], 1)]

        The right eigenvectors are nothing but the left eigenvectors of the
        transpose matrix::

            sage: left  = A.transpose().eigenvectors_left(); left
            [(-1/2*sqrt(17) + 3/2, [(1, -1/2*sqrt(17) + 3/2)], 1), (1/2*sqrt(17) + 3/2, [(1, 1/2*sqrt(17) + 3/2)], 1)]
            sage: right[0][1] == left[0][1]
            True

        TESTS::

            sage: A = matrix(SR, [[1, 2], [3, 4]], sparse=True)
            sage: B = matrix(SR, [[1, 1], [0, 1]], sparse=True)
            sage: A.eigenvectors_right(B)
            Traceback (most recent call last):
            ...
            NotImplementedError: generalized eigenvector decomposition is
            implemented for RDF and CDF, but not for Symbolic Ring

        Check that :issue:`23332` is fixed::

            sage: matrix([[x, x^2], [1, 0]], sparse=True).eigenvectors_right()
            [(-1/2*x*(sqrt(5) - 1), [(1, -1/2*(sqrt(5) + 1)/x)], 1),
             (1/2*x*(sqrt(5) + 1), [(1, 1/2*(sqrt(5) - 1)/x)], 1)]
        """
        return self.transpose().eigenvectors_left(other=other)

    def exp(self):
        r"""
        Return the matrix exponential of this matrix `X`, which is the matrix

        .. MATH::

           e^X = \sum_{k=0}^{\infty} \frac{X^k}{k!}.

        This function depends on maxima's matrix exponentiation
        function, which does not deal well with floating point
        numbers.  If the matrix has floating point numbers, they will
        be rounded automatically to rational numbers during the
        computation.

        EXAMPLES::

            sage: m = matrix(SR,2, [0,x,x,0], sparse=True); m
            [0 x]
            [x 0]
            sage: m.exp()
            [1/2*(e^(2*x) + 1)*e^(-x) 1/2*(e^(2*x) - 1)*e^(-x)]
            [1/2*(e^(2*x) - 1)*e^(-x) 1/2*(e^(2*x) + 1)*e^(-x)]
            sage: exp(m)
            [1/2*(e^(2*x) + 1)*e^(-x) 1/2*(e^(2*x) - 1)*e^(-x)]
            [1/2*(e^(2*x) - 1)*e^(-x) 1/2*(e^(2*x) + 1)*e^(-x)]

        Exponentiation works on 0x0 and 1x1 matrices, but the 1x1 example
        requires a patched version of maxima (:issue:`32898`) for now::

            sage: m = matrix(SR,0,[], sparse=True); m
            []
            sage: m.exp()
            []
            sage: m = matrix(SR,1,[2], sparse=True); m
            [2]
            sage: m.exp()  # not tested, requires patched maxima
            [e^2]

        Commuting matrices `m, n` have the property that
        `e^{m+n} = e^m e^n` (but non-commuting matrices need not)::

            sage: m = matrix(SR,2,[1..4], sparse=True); n = m^2
            sage: m*n
            [ 37  54]
            [ 81 118]
            sage: n*m
            [ 37  54]
            [ 81 118]

            sage: a = exp(m+n) - exp(m)*exp(n)
            sage: a.simplify_rational() == 0
            True

        The input matrix must be square::

            sage: m = matrix(SR,2,3,[1..6], sparse=True); exp(m)
            Traceback (most recent call last):
            ...
            ValueError: exp only defined on square matrices

        In this example we take the symbolic answer and make it
        numerical at the end::

            sage: exp(matrix(SR, [[1.2, 5.6], [3,4]], sparse=True)).change_ring(RDF)  # rel tol 1e-15
            [ 346.5574872980695  661.7345909344504]
            [354.50067371488416  677.4247827652946]

        Another example involving the reversed identity matrix, which
        we clumsily create::

            sage: m = identity_matrix(SR,4, sparse=True)
            sage: m = matrix(list(reversed(m.rows())), sparse=True) * x
            sage: exp(m)
            [1/2*(e^(2*x) + 1)*e^(-x)                        0                        0 1/2*(e^(2*x) - 1)*e^(-x)]
            [                       0 1/2*(e^(2*x) + 1)*e^(-x) 1/2*(e^(2*x) - 1)*e^(-x)                        0]
            [                       0 1/2*(e^(2*x) - 1)*e^(-x) 1/2*(e^(2*x) + 1)*e^(-x)                        0]
            [1/2*(e^(2*x) - 1)*e^(-x)                        0                        0 1/2*(e^(2*x) + 1)*e^(-x)]
        """
        if not self.is_square():
            raise ValueError("exp only defined on square matrices")
        if self.nrows() == 0:
            return self
        # Maxima's matrixexp function chokes on floating point numbers
        # so we automatically convert floats to rationals by passing
        # keepfloat: false
        m = self._maxima_(maxima)
        z = maxima('matrixexp(%s), keepfloat: false' % m.name())
        if self.nrows() == 1:
            # We do the following, because Maxima stupidly exp's 1x1
            # matrices into non-matrices!
            z = maxima('matrix([%s])' % z.name())

        return z._sage_()

    def charpoly(self, var='x', algorithm=None):
        r"""
        Compute the characteristic polynomial of ``self``, using maxima.

        .. NOTE::

            The characteristic polynomial is defined as `\det(xI-A)`.

        INPUT:

        - ``var`` -- (default: ``'x'``) name of variable of charpoly

        EXAMPLES::

            sage: M = matrix(SR, 2, 2, var('a,b,c,d'), sparse=True)
            sage: M.charpoly('t')
            t^2 + (-a - d)*t - b*c + a*d
            sage: matrix(SR, 5, [1..5^2], sparse=True).charpoly()
            x^5 - 65*x^4 - 250*x^3

        TESTS:

        The cached polynomial should be independent of the ``var``
        argument (:issue:`12292`). We check (indirectly) that the
        second call uses the cached value by noting that its result is
        not cached::

            sage: M = MatrixSpace(SR, 2, sparse=True)
            sage: A = M(range(0, 2^2))
            sage: type(A)
            <class 'sage.matrix.matrix_symbolic_sparse.Matrix_symbolic_sparse'>
            sage: A.charpoly('x')
            x^2 - 3*x - 2
            sage: A.charpoly('y')
            y^2 - 3*y - 2
            sage: A._cache['charpoly']
            x^2 - 3*x - 2

        Ensure the variable name of the polynomial does not conflict
        with variables used within the matrix (:issue:`14403`)::

            sage: Matrix(SR, [[sqrt(x), x],[1,x]], sparse=True).charpoly().list()
            [x^(3/2) - x, -x - sqrt(x), 1]

        Test that :issue:`13711` is fixed::

            sage: matrix([[sqrt(2), -1], [pi, e^2]], sparse=True).charpoly()
            x^2 + (-sqrt(2) - e^2)*x + pi + sqrt(2)*e^2

        Test that :issue:`26427` is fixed::

            sage: M = matrix(SR, 7, 7, SR.var('a', 49), sparse=True)
            sage: M.charpoly().degree() # long time
            7
        """
        cache_key = 'charpoly'
        cp = self.fetch(cache_key)
        if cp is not None:
            return cp.change_variable_name(var)
        from sage.symbolic.ring import SR

        # We must not use a variable name already present in the matrix
        vname = 'do_not_use_this_name_in_a_matrix_youll_compute_a_charpoly_of'
        vsym = SR(vname)

        cp = self._maxima_(maxima).charpoly(vname)._sage_().expand()
        cp = [cp.coefficient(vsym, i) for i in range(self.nrows() + 1)]
        cp = SR[var](cp)

        # Maxima has the definition det(matrix-xI) instead of
        # det(xI-matrix), which is what Sage uses elsewhere.  We
        # correct for the discrepancy.
        if self.nrows() % 2 == 1:
            cp = -cp

        self.cache(cache_key, cp)
        return cp

    def minpoly(self, var='x'):
        """
        Return the minimal polynomial of ``self``.

        EXAMPLES::

            sage: M = Matrix.identity(SR, 2, sparse=True)
            sage: M.minpoly()
            x - 1

            sage: t = var('t')
            sage: m = matrix(2, [1, 2, 4, t], sparse=True)
            sage: m.minimal_polynomial()
            x^2 + (-t - 1)*x + t - 8

        TESTS:

        Check that the variable `x` can occur in the matrix::

            sage: m = matrix([[x]], sparse=True)
            sage: m.minimal_polynomial('y')
            y - x
        """
        mp = self.fetch('minpoly')
        if mp is None:
            mp = self._maxima_lib_().jordan().minimalPoly().expand()
            d = mp.hipow('x')
            mp = [mp.coeff('x', i) for i in range(int(d) + 1)]
            mp = PolynomialRing(self.base_ring(), 'x')(mp)
            self.cache('minpoly', mp)
        return mp.change_variable_name(var)

    def fcp(self, var='x'):
        """
        Return the factorization of the characteristic polynomial of ``self``.

        INPUT:

        - ``var`` -- (default: ``'x'``) name of variable of charpoly

        EXAMPLES::

            sage: a = matrix(SR,[[1,2],[3,4]], sparse=True)
            sage: a.fcp()
            x^2 - 5*x - 2
            sage: [i for i in a.fcp()]
            [(x^2 - 5*x - 2, 1)]
            sage: a = matrix(SR,[[1,0],[0,2]], sparse=True)
            sage: a.fcp()
            (x - 2) * (x - 1)
            sage: [i for i in a.fcp()]
            [(x - 2, 1), (x - 1, 1)]
            sage: a = matrix(SR, 5, [1..5^2], sparse=True)
            sage: a.fcp()
            (x^2 - 65*x - 250) * x^3
            sage: list(a.fcp())
            [(x^2 - 65*x - 250, 1), (x, 3)]
        """
        from sage.symbolic.ring import SR
        sub_dict = {var: SR.var(var)}
        return Factorization(self.charpoly(var).subs(**sub_dict).factor_list())

    def jordan_form(self, subdivide=True, transformation=False):
        """
        Return a Jordan normal form of ``self``.

        INPUT:

        - ``self`` -- a square matrix

        - ``subdivide`` -- boolean (default: ``True``)

        - ``transformation`` -- boolean (default: ``False``)

        OUTPUT:

        If ``transformation`` is ``False``, only a Jordan normal form
        (unique up to the ordering of the Jordan blocks) is returned.
        Otherwise, a pair ``(J, P)`` is returned, where ``J`` is a
        Jordan normal form and ``P`` is an invertible matrix such that
        ``self`` equals ``P * J * P^(-1)``.

        If ``subdivide`` is ``True``, the Jordan blocks in the
        returned matrix ``J`` are indicated by a subdivision in
        the sense of :meth:`~sage.matrix.matrix2.subdivide`.

        EXAMPLES:

        We start with some examples of diagonalisable matrices::

            sage: a,b,c,d = var('a,b,c,d')
            sage: matrix([a], sparse=True).jordan_form()
            [a]
            sage: matrix([[a, 0], [1, d]], sparse=True).jordan_form(subdivide=True)
            [d|0]
            [-+-]
            [0|a]
            sage: matrix([[a, 0], [1, d]], sparse=True).jordan_form(subdivide=False)
            [d 0]
            [0 a]
            sage: matrix([[a, x, x], [0, b, x], [0, 0, c]], sparse=True).jordan_form()
            [c|0|0]
            [-+-+-]
            [0|b|0]
            [-+-+-]
            [0|0|a]

        In the following examples, we compute Jordan forms of some
        non-diagonalisable matrices::

            sage: matrix([[a, a], [0, a]], sparse=True).jordan_form()
            [a 1]
            [0 a]
            sage: matrix([[a, 0, b], [0, c, 0], [0, 0, a]], sparse=True).jordan_form()
            [c|0 0]
            [-+---]
            [0|a 1]
            [0|0 a]

        The following examples illustrate the ``transformation`` flag.
        Note that symbolic expressions may need to be simplified to
        make consistency checks succeed::

            sage: A = matrix([[x - a*c, a^2], [-c^2, x + a*c]], sparse=True)
            sage: J, P = A.jordan_form(transformation=True)
            sage: J, P
            (
            [x 1]  [-a*c    1]
            [0 x], [-c^2    0]
            )
            sage: A1 = P * J * ~P; A1
            [             -a*c + x (a*c - x)*a/c + a*x/c]
            [                 -c^2               a*c + x]
            sage: A1.simplify_rational() == A
            True

            sage: B = matrix([[a, b, c], [0, a, d], [0, 0, a]], sparse=True)
            sage: J, T = B.jordan_form(transformation=True)
            sage: J, T
            (
            [a 1 0]  [b*d   c   0]
            [0 a 1]  [  0   d   0]
            [0 0 a], [  0   0   1]
            )
            sage: (B * T).simplify_rational() == T * J
            True

        Finally, some examples involving square roots::

            sage: matrix([[a, -b], [b, a]], sparse=True).jordan_form()
            [a - I*b|      0]
            [-------+-------]
            [      0|a + I*b]
            sage: matrix([[a, b], [c, d]], sparse=True).jordan_form(subdivide=False)
            [1/2*a + 1/2*d - 1/2*sqrt(a^2 + 4*b*c - 2*a*d + d^2)                                                   0]
            [                                                  0 1/2*a + 1/2*d + 1/2*sqrt(a^2 + 4*b*c - 2*a*d + d^2)]
        """
        A = self._maxima_lib_()
        jordan_info = A.jordan()
        J = matrix(jordan_info.dispJordan()._sage_(), sparse=True)
        if subdivide:
            v = [x[1] for x in jordan_info]
            w = [sum(v[0:i]) for i in range(1, len(v))]
            J.subdivide(w, w)
        if transformation:
            P = A.diag_mode_matrix(jordan_info)._sage_()
            return J, matrix(P, sparse=True)
        else:
            return J

    def simplify(self):
        """
        Simplify ``self``.

        EXAMPLES::

            sage: var('x,y,z')
            (x, y, z)
            sage: m = matrix([[z, (x+y)/(x+y)], [x^2, y^2+2]], sparse=True); m
            [      z       1]
            [    x^2 y^2 + 2]
            sage: m.simplify()
            [      z       1]
            [    x^2 y^2 + 2]
        """
        return self.parent()([x.simplify() for x in self.list()])

    def simplify_trig(self):
        """
        EXAMPLES::

            sage: theta = var('theta')
            sage: M = matrix(SR, 2, 2, [cos(theta), sin(theta), -sin(theta), cos(theta)], sparse=True)
            sage: ~M
            [1/cos(theta) - sin(theta)^2/((sin(theta)^2/cos(theta) + cos(theta))*cos(theta)^2)                   -sin(theta)/((sin(theta)^2/cos(theta) + cos(theta))*cos(theta))]
            [                   sin(theta)/((sin(theta)^2/cos(theta) + cos(theta))*cos(theta))                                          1/(sin(theta)^2/cos(theta) + cos(theta))]
            sage: (~M).simplify_trig()
            [ cos(theta) -sin(theta)]
            [ sin(theta)  cos(theta)]
        """
        return self._maxima_(maxima).trigexpand().trigsimp()._sage_()

    def simplify_rational(self):
        """
        EXAMPLES::

            sage: M = matrix(SR, 3, 3, range(9), sparse=True) - var('t')
            sage: (~M*M)[0,0]
            t*(3*(2/t + (6/t + 7)/((t - 3/t - 4)*t))*(2/t + (6/t + 5)/((t - 3/t
            - 4)*t))/(t - (6/t + 7)*(6/t + 5)/(t - 3/t - 4) - 12/t - 8) + 1/t +
            3/((t - 3/t - 4)*t^2)) - 6*(2/t + (6/t + 5)/((t - 3/t - 4)*t))/(t -
            (6/t + 7)*(6/t + 5)/(t - 3/t - 4) - 12/t - 8) - 3*(6/t + 7)*(2/t +
            (6/t + 5)/((t - 3/t - 4)*t))/((t - (6/t + 7)*(6/t + 5)/(t - 3/t -
            4) - 12/t - 8)*(t - 3/t - 4)) - 3/((t - 3/t - 4)*t)
            sage: expand((~M*M)[0,0])
            1
            sage: (~M * M).simplify_rational()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return self._maxima_(maxima).fullratsimp()._sage_()

    def simplify_full(self):
        """
        Simplify a symbolic matrix by calling
        :meth:`Expression.simplify_full()` componentwise.

        INPUT:

        - ``self`` -- the matrix whose entries we should simplify

        OUTPUT: a copy of ``self`` with all of its entries simplified

        EXAMPLES:

        Symbolic matrices will have their entries simplified::

            sage: a,n,k = SR.var('a,n,k')
            sage: f1 = sin(x)^2 + cos(x)^2
            sage: f2 = sin(x/(x^2 + x))
            sage: f3 = binomial(n,k)*factorial(k)*factorial(n-k)
            sage: f4 = x*sin(2)/(x^a)
            sage: A = matrix(SR, [[f1,f2],[f3,f4]], sparse=True)
            sage: A.simplify_full()
            [                1    sin(1/(x + 1))]
            [     factorial(n) x^(-a + 1)*sin(2)]
        """
        M = self.parent()
        return M([expr.simplify_full() for expr in self])

    def canonicalize_radical(self):
        r"""
        Choose a canonical branch of each entry of ``self`` by calling
        :meth:`Expression.canonicalize_radical()` componentwise.

        EXAMPLES::

            sage: var('x','y')
            (x, y)
            sage: l1 = [sqrt(2)*sqrt(3)*sqrt(6) , log(x*y)]
            sage: l2 = [sin(x/(x^2 + x)) , 1]
            sage: m = matrix([l1, l2], sparse=True)
            sage: m
            [sqrt(6)*sqrt(3)*sqrt(2)                log(x*y)]
            [       sin(x/(x^2 + x))                       1]
            sage: m.canonicalize_radical()
            [              6 log(x) + log(y)]
            [ sin(1/(x + 1))               1]
        """
        M = self.parent()
        return M([expr.canonicalize_radical() for expr in self])

    def factor(self):
        """
        Operate point-wise on each element.

        EXAMPLES::

            sage: M = matrix(SR, 2, 2, x^2 - 2*x + 1, sparse=True); M
            [x^2 - 2*x + 1             0]
            [            0 x^2 - 2*x + 1]
            sage: M.factor()
            [(x - 1)^2         0]
            [        0 (x - 1)^2]
        """
        return matrix(self._maxima_(maxima).factor()._sage_(), sparse=True)

    def expand(self):
        """
        Operate point-wise on each element.

        EXAMPLES::

            sage: M = matrix(2, 2, range(4)) - var('x')
            sage: M*M
            [      x^2 + 2      -2*x + 3]
            [     -4*x + 6 (x - 3)^2 + 2]
            sage: (M*M).expand()
            [       x^2 + 2       -2*x + 3]
            [      -4*x + 6 x^2 - 6*x + 11]
        """
        from sage.misc.call import attrcall
        return self.apply_map(attrcall('expand'))

    def variables(self):
        """
        Return the variables of ``self``.

        EXAMPLES::

            sage: var('a,b,c,x,y')
            (a, b, c, x, y)
            sage: m = matrix([[x, x+2], [x^2, x^2+2]], sparse=True); m
            [      x   x + 2]
            [    x^2 x^2 + 2]
            sage: m.variables()
            (x,)
            sage: m = matrix([[a, b+c], [x^2, y^2+2]], sparse=True); m
            [      a   b + c]
            [    x^2 y^2 + 2]
            sage: m.variables()
            (a, b, c, x, y)
        """
        vars = set(sum([op.variables() for op in self.list()], ()))
        return tuple(sorted(vars, key=repr))

    def arguments(self):
        """
        Return a tuple of the arguments that ``self`` can take.

        EXAMPLES::

            sage: var('x,y,z')
            (x, y, z)
            sage: M = MatrixSpace(SR,2,2, sparse=True)
            sage: M(x).arguments()
            (x,)
            sage: M(x+sin(x)).arguments()
            (x,)
        """
        return self.variables()

    def number_of_arguments(self):
        """
        Return the number of arguments that ``self`` can take.

        EXAMPLES::

            sage: var('a,b,c,x,y')
            (a, b, c, x, y)
            sage: m = matrix([[a, (x+y)/(x+y)], [x^2, y^2+2]], sparse=True); m
            [      a       1]
            [    x^2 y^2 + 2]
            sage: m.number_of_arguments()
            3
        """
        return len(self.variables())

    def __call__(self, *args, **kwargs):
        """
        EXAMPLES::

            sage: var('x,y,z')
            (x, y, z)
            sage: M = MatrixSpace(SR,2,2, sparse=True)
            sage: h = M(sin(x)+cos(x))
            sage: h
            [cos(x) + sin(x)               0]
            [              0 cos(x) + sin(x)]
            sage: h(x=1)
            [cos(1) + sin(1)               0]
            [              0 cos(1) + sin(1)]
            sage: h(x=x)
            [cos(x) + sin(x)               0]
            [              0 cos(x) + sin(x)]

            sage: h = M((sin(x)+cos(x)).function(x))
            sage: h
            [cos(x) + sin(x)               0]
            [              0 cos(x) + sin(x)]

            sage: f = M([0,x,y,z]); f
            [0 x]
            [y z]
            sage: f.arguments()
            (x, y, z)
            sage: f()
            [0 x]
            [y z]
            sage: f(x=1)
            [0 1]
            [y z]
            sage: f(x=1,y=2)
            [0 1]
            [2 z]
            sage: f(x=1,y=2,z=3)
            [0 1]
            [2 3]
            sage: f({x:1,y:2,z:3})
            [0 1]
            [2 3]

        TESTS::

            sage: f(1, x=2)
            Traceback (most recent call last):
            ...
            ValueError: args and kwargs cannot both be specified
            sage: f(x=1,y=2,z=3,t=4)
            [0 1]
            [2 3]

            sage: h(1)
            Traceback (most recent call last):
            ...
            ValueError: use named arguments, like EXPR(x=..., y=...)
        """
        if kwargs and args:
            raise ValueError("args and kwargs cannot both be specified")

        if args:
            if len(args) == 1 and isinstance(args[0], dict):
                kwargs = {repr(x): vx for x, vx in args[0].iteritems()}
            else:
                raise ValueError('use named arguments, like EXPR(x=..., y=...)')

        new_entries = []
        for entry in self.list():
            try:
                new_entries.append(entry(**kwargs))
            except ValueError:
                new_entries.append(entry)

        return self.parent(new_entries)

    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j) except -1:
        r"""
        Return 1 if the entry ``(i, j)`` is zero, otherwise 0.

        EXAMPLES::

            sage: M = matrix(SR, [[0,1,0],[0,0,0]], sparse=True)
            sage: M.zero_pattern_matrix()  # indirect doctest
            [1 0 1]
            [1 1 1]
        """
        entry = self.get_unsafe(i, j)
        # See if we can avoid the full proof machinery that the entry is 0
        if entry.is_trivial_zero():
            return 1
        if entry:
            return 0
        else:
            return 1

    def function(self, *args):
        """
        Return a matrix over a callable symbolic expression ring.

        EXAMPLES::

            sage: x, y = var('x,y')
            sage: v = matrix([[x,y],[x*sin(y), 0]], sparse=True)
            sage: w = v.function([x,y]); w
            [       (x, y) |--> x        (x, y) |--> y]
            [(x, y) |--> x*sin(y)        (x, y) |--> 0]
            sage: w.parent()
            Full MatrixSpace of 2 by 2 sparse matrices over Callable function ring with arguments (x, y)
        """
        from sage.symbolic.callable import CallableSymbolicExpressionRing
        return matrix(CallableSymbolicExpressionRing(args),
                      self.nrows(), self.ncols(), self.list(), sparse=True)
