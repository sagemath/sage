# sage_setup: distribution = sagemath-singular
"""
Dense matrices over multivariate polynomials over fields

This implementation inherits from Matrix_generic_dense, i.e. it is not
optimized for speed only some methods were added.

AUTHOR:

* Martin Albrecht <malb@informatik.uni-bremen.de>
"""

# ****************************************************************************
#       Copyright (C) 2013 Martin Albrecht <malb@informatik.uni-bremen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.categories.fields import Fields
from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.matrix.matrix2 cimport Matrix

from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular
from sage.rings.polynomial.polynomial_singular_interface import can_convert_to_singular

from sage.libs.singular.function import singular_function


cdef class Matrix_mpolynomial_dense(Matrix_generic_dense):
    """
    Dense matrix over a multivariate polynomial ring over a field.
    """
    def echelon_form(self, algorithm='row_reduction', **kwds):
        """
        Return an echelon form of ``self`` using chosen algorithm.

        By default only a usual row reduction with no divisions or column swaps
        is returned.

        If Gauss-Bareiss algorithm is chosen, column swaps are recorded and can
        be retrieved via :meth:`swapped_columns`.

        INPUT:

        - ``algorithm`` -- string, which algorithm to use (default:
          'row_reduction'). Valid options are:

            - ``'row_reduction'`` (default) -- reduce as far as
              possible, only divide by constant entries

            - ``'frac'`` -- reduced echelon form over fraction field

            - ``'bareiss'`` -- fraction free Gauss-Bareiss algorithm
              with column swaps

        OUTPUT:

        The row echelon form of A depending on the chosen algorithm,
        as an immutable matrix.  Note that ``self`` is *not* changed
        by this command. Use ``A.echelonize()`` to change `A` in
        place.

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(GF(127), 2)
            sage: A = matrix(P, 2, 2, [1, x, 1, y])
            sage: A
            [1 x]
            [1 y]
            sage: A.echelon_form()
            [     1      x]
            [     0 -x + y]


        The reduced row echelon form over the fraction field is as follows::

            sage: A.echelon_form('frac') # over fraction field
            [1 0]
            [0 1]

        Alternatively, the Gauss-Bareiss algorithm may be chosen::

            sage: E = A.echelon_form('bareiss'); E
            [    1     y]
            [    0 x - y]

        After the application of the Gauss-Bareiss algorithm the swapped columns
        may inspected::

            sage: E.swapped_columns(), E.pivots()
            ((0, 1), (0, 1))
            sage: A.swapped_columns(), A.pivots()
            (None, (0, 1))

        Another approach is to row reduce as far as possible::

            sage: A.echelon_form('row_reduction')
            [     1      x]
            [     0 -x + y]
        """
        x = self.fetch('echelon_form_'+algorithm)
        if x is not None: return x

        if algorithm == "frac":
            E = self.matrix_over_field()
            E.echelonize(**kwds)
        else:
            E = self.__copy__()
            E.echelonize(algorithm=algorithm, **kwds)

        E.set_immutable()  # so we can cache the echelon form.
        self.cache('echelon_form_'+algorithm, E)

        if algorithm == "frac":
            self.cache('pivots', E.pivots())
        elif algorithm == "bareiss":
            l = E.swapped_columns()
            self.cache('pivots', tuple(sorted(l)))
        elif algorithm == "row_reduction":
            pass

        return E

    def pivots(self):
        """
        Return the pivot column positions of this matrix as a list of integers.

        This returns a list, of the position of the first nonzero entry in each
        row of the echelon form.

        OUTPUT: list of Python ints

        EXAMPLES::

            sage: matrix([PolynomialRing(GF(2), 2, 'x').gen()]).pivots()
            (0,)
            sage: K = QQ['x,y']
            sage: x, y = K.gens()
            sage: m = matrix(K, [(-x, 1, y, x - y), (-x*y, y, y^2 - 1, x*y - y^2 + x), (-x*y + x, y - 1, y^2 - y - 2, x*y - y^2 + x + y)])
            sage: m.pivots()
            (0, 2)
            sage: m.rank()
            2
        """
        x = self.fetch('pivots')
        if x is not None:
            return x

        self.echelon_form('frac')
        x = self.fetch('pivots')

        if x is None:
            raise RuntimeError("BUG: matrix pivots should have been set but weren't, matrix parent = '%s'"%self.parent())
        return x

    def echelonize(self, algorithm='row_reduction', **kwds):
        """
        Transform ``self`` into a matrix in echelon form over the same base ring as
        ``self``.

        If Gauss-Bareiss algorithm is chosen, column swaps are recorded and can
        be retrieved via :meth:`swapped_columns`.

        INPUT:

        - ``algorithm`` -- string, which algorithm to use. Valid options are:

            - ``'row_reduction'`` -- reduce as far as possible, only
              divide by constant entries

            - ``'bareiss'`` -- fraction free Gauss-Bareiss algorithm
              with column swaps

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQ, 2)
            sage: A = matrix(P, 2, 2, [1/2, x, 1, 3/4*y+1])
            sage: A
            [      1/2         x]
            [        1 3/4*y + 1]

            sage: B = copy(A)
            sage: B.echelonize('bareiss'); B
            [              1       3/4*y + 1]
            [              0 x - 3/8*y - 1/2]

            sage: B = copy(A)
            sage: B.echelonize('row_reduction'); B
            [               1              2*x]
            [               0 -2*x + 3/4*y + 1]

            sage: P.<x,y> = PolynomialRing(QQ, 2)
            sage: A = matrix(P,2,3,[2,x,0,3,y,1]); A
            [2 x 0]
            [3 y 1]

            sage: E = A.echelon_form('bareiss'); E
            [1 3 y]
            [0 2 x]
            sage: E.swapped_columns()
            (2, 0, 1)
            sage: A.pivots()
            (0, 1, 2)
        """
        self.check_mutability()

        if self._nrows == 0 or self._ncols == 0:
            self.cache('in_echelon_form_'+algorithm, True)
            self.cache('rank', 0)
            self.cache('pivots', tuple())
            return

        x = self.fetch('in_echelon_form_'+algorithm)
        if x is not None:
            return  # already known to be in echelon form

        if algorithm == 'bareiss':
            self._echelonize_gauss_bareiss()
        elif algorithm == 'row_reduction':
            self._echelonize_row_reduction()
        else:
            raise ValueError("Unknown algorithm '%s'" % algorithm)

    def _echelonize_gauss_bareiss(self):
        """
        Transform this matrix into a matrix in upper triangular form over the
        same base ring as ``self`` using the fraction free Gauss-Bareiss
        algorithm with column swaps.

        The performed column swaps can be accessed via
        :meth:`swapped_columns`.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: C = matrix(R, [[-6/5*x*y - y^2, -6*y^2 - 1/4*y],
            ....:                [  -1/3*x*y - 3, x*y - x]])
            sage: E = C.echelon_form('bareiss')     # indirect doctest
            sage: E
            [ -1/3*x*y - 3                                                          x*y - x]
            [            0 6/5*x^2*y^2 + 3*x*y^3 - 6/5*x^2*y - 11/12*x*y^2 + 18*y^2 + 3/4*y]
            sage: E.swapped_columns()
            (0, 1)

        ALGORITHM:

        Uses libSINGULAR or SINGULAR
        """
        cdef R = self.base_ring()

        x = self.fetch('in_echelon_form_bareiss')
        if x is not None:
            return  # already known to be in echelon form

        if isinstance(self.base_ring(), MPolynomialRing_libsingular):

            self.check_mutability()
            self.clear_cache()

            singular_bareiss = singular_function("bareiss")
            E, l = singular_bareiss(self.T)

            m = len(E)
            n = len(E[0])

            for r in range(m):
                for c in range(n):
                    self.set_unsafe(r, c, E[r][c])
                for c in range(n, self.ncols()):
                    self.set_unsafe(r, c, R._zero_element)
            for r in range(m, self.nrows()):
                for c in range(self.ncols()):
                    self.set_unsafe(r, c, R._zero_element)

            from sage.rings.integer_ring import ZZ
            l = [ZZ(e-1) for e in l]

            self.cache('in_echelon_form_bareiss',True)
            self.cache('rank', len(E))
            self.cache('pivots', tuple(range(len(E))))
            self.cache('swapped_columns', tuple(l))

        elif can_convert_to_singular(self.base_ring()):

            self.check_mutability()
            self.clear_cache()

            E,l = self.T._singular_().bareiss()._sage_(self.base_ring())

            # clear matrix
            for r from 0 <= r < self._nrows:
                for c from 0 <= c < self._ncols:
                    self.set_unsafe(r,c,R._zero_element)

            for r from 0 <= r < E.nrows():
                for c from 0 <= c < E.ncols():
                    self.set_unsafe(c,r, E[r,c])

            self.cache('in_echelon_form_bareiss',True)
            self.cache('rank', E.nrows())
            self.cache('pivots', tuple(range(E.nrows())))
            self.cache('swapped_columns', l)

        else:

            raise NotImplementedError("cannot apply Gauss-Bareiss algorithm over this base ring")

    def _echelonize_row_reduction(self):
        """
        Transform this matrix to a matrix in row reduced form as far as this is
        possible, i.e., only perform division by constant elements.

        EXAMPLES:

        If all entries are constant, then this method performs the same
        operations as ``echelon_form('default')`` ::

            sage: P.<x0,x1,y0,y1> = PolynomialRing(GF(127), 4)
            sage: A = Matrix(P, 4, 4, [-14,0,45,-55,-61,-16,0,0,0,0,25,-62,-22,0,52,0]); A
            [-14   0  45 -55]
            [-61 -16   0   0]
            [  0   0  25 -62]
            [-22   0  52   0]
            sage: E1 = A.echelon_form()
            sage: E2 = A.echelon_form('row_reduction')  # indirect doctest
            sage: E1 == E2
            True

        If no entries are constant, nothing happens::

            sage: P.<x0,x1,y0,y1> = PolynomialRing(GF(2), 4)
            sage: A = Matrix(P, 2, 2, [x0, y0, x0, y0]); A
            [x0 y0]
            [x0 y0]

            sage: B = copy(A)
            sage: A.echelonize('row_reduction') # modifies self
            sage: B == A
            True

        A more interesting example::

            sage: P.<x0,x1,y0,y1> = PolynomialRing(GF(2), 4)
            sage: l = [1, 1, 1, 1,     1,
            ....:      0, 1, 0, 1,    x0,
            ....:      0, 0, 1, 1,    x1,
            ....:      1, 1, 0, 0,    y0,
            ....:      0, 1, 0, 1,    y1,
            ....:      0, 1, 0, 0, x0*y0,
            ....:      0, 1, 0, 1, x0*y1,
            ....:      0, 0, 0, 0, x1*y0,
            ....:      0, 0, 0, 1, x1*y1]
            sage: A = Matrix(P, 9, 5, l)
            sage: B = A.__copy__()
            sage: B.echelonize('row_reduction'); B
            [                 1                  0                  0                  0     x0*y0 + x1 + 1]
            [                 0                  1                  0                  0              x0*y0]
            [                 0                  0                  1                  0    x0*y0 + x0 + x1]
            [                 0                  0                  0                  1         x0*y0 + x0]
            [                 0                  0                  0                  0            x0 + y1]
            [                 0                  0                  0                  0        x1 + y0 + 1]
            [                 0                  0                  0                  0         x0*y1 + x0]
            [                 0                  0                  0                  0              x1*y0]
            [                 0                  0                  0                  0 x0*y0 + x1*y1 + x0]

        This is the same result as SINGULAR's ``rowred`` command which
        returns::

            sage: E = A._singular_().rowred()._sage_(P)
            sage: E == B
            True

        ALGORITHM:

        Gaussian elimination with division limited to constant
        entries. Based on SINGULAR's rowred command.
        """
        from sage.matrix.constructor import matrix

        cdef int c, r, i, j, rc, start_row, nr, nc

        x = self.fetch('in_echelon_form_row_reduction')
        if x is not None:
            return  # already known to be in echelon form

        nr,nc = self.nrows(),self.ncols()
        F = self.base_ring().base_ring()
        cdef Matrix d = matrix(F, nr, nc)
        start_row = 0

        for r from 0 <= r < nr:
            for c from 0 <= c < nc:
                p = self.get_unsafe(r, c)
                if p.is_constant():
                    d.set_unsafe(r, c, p.constant_coefficient())

        for c from 0 <= c < nc:
            r = -1
            for rc from start_row <= rc < nr:
                if d.get_unsafe(rc, c):
                    r = rc
                    break
            if r!=-1:
                a_inverse = ~self.get_unsafe(r, c)
                self.rescale_row_c(r, a_inverse, c)
                self.swap_rows_c(r, start_row)

                for i from 0 <= i < nr:
                    if i != start_row:
                        minus_b = -self.get_unsafe(i, c)
                        self.add_multiple_of_row(i, start_row, minus_b, 0)

                start_row +=1

                d = d._parent(0)
                for i from start_row <= i < nr:
                    for j from c+1 <= j < nc:
                        if self.get_unsafe(i, j).is_constant():
                            d.set_unsafe(i, j, self.get_unsafe(i, j).constant_coefficient())

        self.cache('in_echelon_form_row_reduction', True)

    def swapped_columns(self):
        """
        Return which columns were swapped during the Gauss-Bareiss reduction.

        OUTPUT:

        Return a tuple representing the column swaps during the last application
        of the Gauss-Bareiss algorithm (see :meth:`echelon_form` for details).

        The tuple as length equal to the rank of ``self`` and the value at the
        `i`-th position indicates the source column which was put as the `i`-th
        column.

        If no Gauss-Bareiss reduction was performed yet, None is
        returned.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: C = random_matrix(R, 2, 2, terms=2)
            sage: while C.rank() != 2:
            ....:     C = random_matrix(R, 2, 2, terms=2)
            sage: C.swapped_columns()
            sage: E = C.echelon_form('bareiss')
            sage: sorted(E.swapped_columns())
            [0, 1]
        """
        return self.fetch('swapped_columns')

    def _fitting_ideal(self, i):
        r"""
        Return the `i`-th Fitting ideal of the matrix. This is the ideal generated
        by the `n - i` minors, where `n` is the number of columns.

        INPUT:

        - ``i`` -- integer

        OUTPUT: an ideal on the base ring

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: M = matrix(R, [[2*x-z, 0, y-z^2, 0], [0, z - y, z - x, 0],[z - y, x^2 - y, 0, z]])
            sage: M
            [ 2*x - z        0 -z^2 + y        0]
            [       0   -y + z   -x + z        0]
            [  -y + z  x^2 - y        0        z]
            sage: M.fitting_ideal(0)
            Ideal (0) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: M.fitting_ideal(1)
            Ideal (2*x^4 - 3*x^3*z + x^2*z^2 + y^2*z^2 - 2*y*z^3 + z^4 - 2*x^2*y - y^3 + 3*x*y*z + 2*y^2*z - 2*y*z^2, y*z^3 - z^4 - y^2*z + y*z^2, -2*x*y*z + 2*x*z^2 + y*z^2 - z^3, -2*x^2*z + 3*x*z^2 - z^3) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: M.fitting_ideal(2)
            Ideal (-x^3 + x^2*z + x*y - y*z, -y^2 + 2*y*z - z^2, -x^2*z^2 + x^2*y + y*z^2 - y^2, 2*x^3 - x^2*z - 2*x*y + y*z, -x*y + x*z + y*z - z^2, -y*z^2 + z^3 + y^2 - y*z, -2*x*y + 2*x*z + y*z - z^2, y*z^2 - z^3 - y^2 + y*z, 2*x^2 - 3*x*z + z^2, 2*x*z - z^2, -z^3 + y*z, -y*z + z^2, -x*z + z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: M.fitting_ideal(3)
            Ideal (2*x - z, -z^2 + y, -y + z, -x + z, -y + z, x^2 - y, z) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: M.fitting_ideal(4)
            Ideal (1) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: [R.ideal(M.minors(i)) == M._fitting_ideal(4 - i) for i in range(5)]
            [True, True, True, True, True]
        """
        minor = singular_function("minor")
        R = self.base_ring()
        for (nrow, row) in enumerate(self.rows()):
            if not row:
                N = self.delete_rows([nrow])
                return N._fitting_ideal(i)
            for (ncoef, coef) in enumerate(row):
                if all(coef.divides(f) for f in row):
                    N = self.__copy__()
                    for j in range(self.ncols()):
                        if j != ncoef:
                            N.add_multiple_of_column(j, ncoef, -R(self[nrow,j] / coef))
                    return N.fitting_ideal(i)
        for (ncolumn, column) in enumerate(self.columns()):
            if not column:
                N = self.delete_columns([ncolumn])
                return N._fitting_ideal(i-1)
            for (ncoef, coef) in enumerate(column):
                if all(coef.divides(f) for f in column):
                    N = self.__copy__()
                    for j in range(self.nrows()):
                        if j != ncoef:
                            N.add_multiple_of_row(j, ncoef, -R(self[j, ncolumn] / coef))
                    return N.fitting_ideal(i)
        rank = self.ncols() - i
        return R.ideal(minor(self, rank))

    def determinant(self, algorithm=None):
        """
        Return the determinant of this matrix.

        INPUT:

        - ``algorithm`` -- ignored

        EXAMPLES:

        We compute the determinant of the arbitrary `3x3` matrix::

            sage: R = PolynomialRing(QQ, 9, 'x')
            sage: A = matrix(R, 3, R.gens())
            sage: A
            [x0 x1 x2]
            [x3 x4 x5]
            [x6 x7 x8]
            sage: A.determinant()
            -x2*x4*x6 + x1*x5*x6 + x2*x3*x7 - x0*x5*x7 - x1*x3*x8 + x0*x4*x8

        We check if two implementations agree on the result::

            sage: R.<x,y> = QQ[]
            sage: C = matrix(R, [[-6/5*x*y - y^2, -6*y^2 - 1/4*y],
            ....:                [  -1/3*x*y - 3, x*y - x]])
            sage: C.determinant()
            -6/5*x^2*y^2 - 3*x*y^3 + 6/5*x^2*y + 11/12*x*y^2 - 18*y^2 - 3/4*y

            sage: C.change_ring(R.change_ring(QQbar)).det()
            (-6/5)*x^2*y^2 + (-3)*x*y^3 + 6/5*x^2*y + 11/12*x*y^2 + (-18)*y^2 + (-3/4)*y

        Finally, we check whether the Singular interface is working::

            sage: R.<x,y> = RR[]
            sage: C = matrix(R, [[0.368965517352886*y^2 + 0.425700773972636*x, -0.800362171389760*y^2 - 0.807635502485287],
            ....:                [0.706173539423122*y^2 - 0.915986060298440, 0.897165181570476*y + 0.107903328188376]])
            sage: C.determinant()
            0.565194587390682*y^4 + 0.33102301536914...*y^3 + 0.381923912175852*x*y - 0.122977163520282*y^2 + 0.0459345303240150*x - 0.739782862078649

        ALGORITHM: Calls Singular, libSingular or native implementation.

        TESTS::

            sage: R = PolynomialRing(QQ, 9, 'x')
            sage: matrix(R, 0, 0).det()
            1

            sage: R.<h,y> = QQ[]
            sage: m = matrix([[y,y,y,y]] * 4)  # larger than 3x3
            sage: m.charpoly()   # put charpoly in the cache
            x^4 - 4*y*x^3
            sage: m.det()
            0

        Check :issue:`23535` is fixed::

            sage: x = polygen(QQ)
            sage: K.<a,b> = NumberField([x^2 - 2, x^2 - 5])
            sage: R.<s,t> = K[]
            sage: m = matrix(R, 4, [y^i for i in range(4) for y in [a,b,s,t]])
            sage: m.det()
            (a - b)*s^3*t^2 + (-a + b)*s^2*t^3 + 3*s^3*t + (-3)*s*t^3 + (-5*a + 2*b)*s^3 + (2*a - 5*b)*s^2*t +
            (-2*a + 5*b)*s*t^2 + (5*a - 2*b)*t^3 + 3*b*a*s^2 + (-3*b*a)*t^2 + (10*a - 10*b)*s + (-10*a + 10*b)*t
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")

        d = self.fetch('det')

        if d is not None:
            return d

        if self._nrows == 0:
            d = self._coerce_element(1)
        elif self._nrows <= 3:
            from sage.matrix.matrix2 import Matrix
            d = Matrix.determinant(self)

        elif self.fetch('charpoly') is not None:
            # if charpoly known, then det is easy.
            D = self.fetch('charpoly')
            if D is not None:
                c = D[0]
                if self._nrows % 2:
                    c = -c
                d = self._coerce_element(c)
        else:
            R = self._base_ring

            if isinstance(R, MPolynomialRing_libsingular) and R.base_ring() in Fields():
                singular_det = singular_function("det")
                d = singular_det(self)

            elif can_convert_to_singular(self.base_ring()):
                d = R(self._singular_().det())

            else:
                from sage.matrix.matrix2 import Matrix
                d = Matrix.determinant(self)

        self.cache('det', d)
        return d
