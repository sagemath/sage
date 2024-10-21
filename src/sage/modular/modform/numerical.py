# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.flint
"""
Numerical computation of newforms
"""

# ****************************************************************************
#       Copyright (C) 2004-2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.fast_arith import prime_range
from sage.matrix.constructor import matrix
from sage.misc.verbose import verbose
from sage.misc.cachefunc import cached_method
from sage.misc.prandom import randint
from sage.modular.arithgroup.all import Gamma0
from sage.modular.modsym.all import ModularSymbols
from sage.modules.free_module_element import free_module_element as vector
from sage.rings.complex_double import CDF
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.structure.richcmp import richcmp_method, richcmp
from sage.structure.sage_object import SageObject
from sage.structure.sequence import Sequence

# This variable controls importing the SciPy library sparingly
scipy = None


@richcmp_method
class NumericalEigenforms(SageObject):
    """
    numerical_eigenforms(group, weight=2, eps=1e-20, delta=1e-2, tp=[2,3,5]).

    INPUT:

    - ``group`` -- a congruence subgroup of a Dirichlet character of
      order 1 or 2

    - ``weight`` -- integer >= 2

    - ``eps`` -- a small float; abs( ) < eps is what "equal to zero" is
      interpreted as for floating point numbers

    - ``delta`` -- a small-ish float; eigenvalues are considered distinct
      if their difference has absolute value at least delta

    - ``tp`` -- use the Hecke operators T_p for p in tp when searching
      for a random Hecke operator with distinct Hecke eigenvalues

    OUTPUT: a numerical eigenforms object, with the following useful methods:

    - :meth:`ap` -- return all eigenvalues of `T_p`

    - :meth:`eigenvalues` -- list of eigenvalues corresponding
      to the given list of primes, e.g.,::

          [[eigenvalues of T_2],
           [eigenvalues of T_3],
           [eigenvalues of T_5], ...]

    - :meth:`systems_of_eigenvalues` -- list of the systems of
      eigenvalues of eigenforms such that the chosen random linear
      combination of Hecke operators has multiplicity 1 eigenvalues.

    EXAMPLES::

        sage: n = numerical_eigenforms(23)
        sage: n == loads(dumps(n))
        True
        sage: n.ap(2)  # abs tol 1e-12
        [3.0, -1.6180339887498947, 0.6180339887498968]
        sage: n.systems_of_eigenvalues(7)  # abs tol 2e-12
        [
        [-1.6180339887498947, 2.2360679774997894, -3.2360679774997894],
        [0.6180339887498968, -2.236067977499788, 1.2360679774997936],
        [3.0, 4.0, 6.0]
        ]
        sage: n.systems_of_abs(7)  # abs tol 2e-12
        [
        [0.6180339887498943, 2.2360679774997894, 1.2360679774997887],
        [1.6180339887498947, 2.23606797749979, 3.2360679774997894],
        [3.0, 4.0, 6.0]
        ]
        sage: n.eigenvalues([2,3,5])  # rel tol 2e-12
        [[3.0, -1.6180339887498947, 0.6180339887498968],
         [4.0, 2.2360679774997894, -2.236067977499788],
         [6.0, -3.2360679774997894, 1.2360679774997936]]
    """
    def __init__(self, group, weight=2, eps=1e-20,
                 delta=1e-2, tp=[2,3,5]):
        """
        Create a new space of numerical eigenforms.

        EXAMPLES::

            sage: numerical_eigenforms(61) # indirect doctest
            Numerical Hecke eigenvalues for Congruence Subgroup Gamma0(61) of weight 2
        """
        if isinstance(group, (int, Integer)):
            group = Gamma0(Integer(group))
        self._group = group
        self._weight = Integer(weight)
        self._tp = tp
        if self._weight < 2:
            raise ValueError("weight must be at least 2")
        self._eps = eps
        self._delta = delta

    def __richcmp__(self, other, op):
        """
        Compare two spaces of numerical eigenforms.

        They are considered equal if and only if they come from the
        same space of modular symbols.

        EXAMPLES::

            sage: n = numerical_eigenforms(23)
            sage: n == loads(dumps(n))
            True
        """
        if not isinstance(other, NumericalEigenforms):
            return NotImplemented
        return richcmp(self.modular_symbols(), other.modular_symbols(), op)

    def level(self):
        """
        Return the level of this set of modular eigenforms.

        EXAMPLES::

            sage: n = numerical_eigenforms(61) ; n.level()
            61
        """
        return self._group.level()

    def weight(self):
        """
        Return the weight of this set of modular eigenforms.

        EXAMPLES::

            sage: n = numerical_eigenforms(61) ; n.weight()
            2
        """
        return self._weight

    def _repr_(self):
        """
        Print string representation of ``self``.

        EXAMPLES::

            sage: n = numerical_eigenforms(61) ; n
            Numerical Hecke eigenvalues for Congruence Subgroup Gamma0(61) of weight 2

            sage: n._repr_()
            'Numerical Hecke eigenvalues for Congruence Subgroup Gamma0(61) of weight 2'
        """
        return "Numerical Hecke eigenvalues for %s of weight %s" % (
            self._group, self._weight)

    @cached_method
    def modular_symbols(self):
        """
        Return the space of modular symbols used for computing this
        set of modular eigenforms.

        EXAMPLES::

            sage: n = numerical_eigenforms(61) ; n.modular_symbols()
            Modular Symbols space of dimension 5 for Gamma_0(61) of weight 2 with sign 1 over Rational Field
        """
        M = ModularSymbols(self._group,
                self._weight, sign=1)
        if M.base_ring() != QQ:
            raise ValueError("modular forms space must be defined over QQ")
        return M

    @cached_method
    def _eigenvectors(self):
        r"""
        Find numerical approximations to simultaneous eigenvectors in
        self.modular_symbols() for all T_p in self._tp.

        EXAMPLES::

            sage: n = numerical_eigenforms(61)
            sage: n._eigenvectors() # random order
            [              1.0    0.289473640239    0.176788851952    0.336707726757  2.4182243084e-16]
            [                0  -0.0702748344418    0.491416161212    0.155925712173    0.707106781187]
            [                0    0.413171180356    0.141163094698   0.0923242547901    0.707106781187]
            [                0    0.826342360711    0.282326189397     0.18464850958 6.79812569682e-16]
            [                0      0.2402380858    0.792225196393    0.905370774276 4.70805946682e-16]

        TESTS:

        This tests if this routine selects only eigenvectors with
        multiplicity one.  Two of the eigenvalues are
        (roughly) -92.21 and -90.30 so if we set ``eps = 2.0``
        then they should compare as equal, causing both eigenvectors
        to be absent from the matrix returned.  The remaining eigenvalues
        (ostensibly unique) are visible in the test, which should be
        independent of which eigenvectors are returned, but it does presume
        an ordering of these eigenvectors for the test to succeed.
        This exercises a correction in :issue:`8018`. ::

            sage: n = numerical_eigenforms(61, eps=2.0)
            sage: evectors = n._eigenvectors()
            sage: evalues = [(matrix((n._hecke_matrix*evectors).column(i))/matrix(evectors.column(i)))[0, 0]
            ....:            for i in range(evectors.ncols())]
            sage: diff = n._hecke_matrix*evectors - evectors*diagonal_matrix(evalues)
            sage: sum(abs(a) for a in diff.list()) < 1.0e-9
            True
        """
        verbose('Finding eigenvector basis')
        M = self.modular_symbols()

        tp = self._tp
        p = tp[0]
        t = M.T(p).matrix()
        for p in tp[1:]:
            t += randint(-50,50)*M.T(p).matrix()

        self._hecke_matrix = t

        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        evals,eig = scipy.linalg.eig(self._hecke_matrix.numpy(), right=True, left=False)
        B = matrix(eig)
        v = [CDF(evals[i]) for i in range(len(evals))]

        # Determine the eigenvectors with eigenvalues of multiplicity
        # one, with equality controlled by the value of eps
        # Keep just these eigenvectors
        eps = self._eps
        w = []
        for i in range(len(v)):
            e = v[i]
            uniq = True
            for j in range(len(v)):
                if uniq and i != j and abs(e-v[j]) < eps:
                    uniq = False
            if uniq:
                w.append(i)
        return B.matrix_from_columns(w)

    @cached_method
    def _easy_vector(self):
        """
        Return a very sparse vector v such that v times the eigenvector matrix
        has all entries nonzero.

        ALGORITHM:

        1. Choose row with the most nonzero entries.   (put 1 there)

        2. Consider submatrix of columns corresponding to zero entries
           in row chosen in 1.

        3. Find row of submatrix with most nonzero entries, and add
           appropriate multiple.  Repeat.

        EXAMPLES::

            sage: n = numerical_eigenforms(37)
            sage: n._easy_vector()                 # slightly random output
            (1.0, 1.0, 0)
            sage: n = numerical_eigenforms(43)
            sage: n._easy_vector()                 # slightly random output
            (1.0, 0, 1.0, 0)
            sage: n = numerical_eigenforms(125)
            sage: n._easy_vector()                 # slightly random output
            (0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        """
        E = self._eigenvectors()
        delta = self._delta
        x = (CDF**E.nrows()).zero_vector()
        if E.nrows() == 0:
            return x

        def best_row(M):
            """
            Find the best row among rows of M, i.e. the row
            with the most entries supported outside [-delta, delta].

            EXAMPLES::

                sage: numerical_eigenforms(61)._easy_vector() # indirect doctest
                (1.0, 0.0, 0.0, 0.0, 1.0)
            """
            R = M.rows()
            v = [len(support(r, delta)) for r in R]
            m = max(v)
            i = v.index(m)
            return i, R[i]

        i, e = best_row(E)

        x[i] = 1

        while True:
            s = set(support(e, delta))
            zp = [j for j in range(e.degree()) if j not in s]
            if not zp:
                break
            C = E.matrix_from_columns(zp)
            # best row
            i, f = best_row(C)
            x[i] += 1   # simplistic
            e = x * E

        self.__easy_vector = x
        return x

    @cached_method
    def _eigendata(self):
        """
        Return all eigendata for self._easy_vector().

        EXAMPLES::

            sage: numerical_eigenforms(61)._eigendata() # random order
            ((1.0, 0.668205013164, 0.219198805797, 0.49263343893, 0.707106781187), (1.0, 1.49654668896, 4.5620686498, 2.02990686579, 1.41421356237), [0, 1], (1.0, 1.0))
        """
        x = self._easy_vector()

        B = self._eigenvectors()

        def phi(y):
            """
            Take coefficients and a basis, and return that
            linear combination of basis vectors.

            EXAMPLES::

                sage: n = numerical_eigenforms(61) # indirect doctest
                sage: n._eigendata() # random order
                ((1.0, 0.668205013164, 0.219198805797, 0.49263343893, 0.707106781187), (1.0, 1.49654668896, 4.5620686498, 2.02990686579, 1.41421356237), [0, 1], (1.0, 1.0))
            """
            return y.element() * B

        phi_x = phi(x)
        V = phi_x.parent()
        phi_x_inv = V([a**(-1) for a in phi_x])
        eps = self._eps
        nzp = support(x, eps)
        x_nzp = vector(CDF, x.list_from_positions(nzp))
        self.__eigendata = (phi_x, phi_x_inv, nzp, x_nzp)
        return self.__eigendata

    @cached_method
    def ap(self, p):
        """
        Return a list of the eigenvalues of the Hecke operator `T_p`
        on all the computed eigenforms.  The eigenvalues match up
        between one prime and the next.

        INPUT:

        - ``p`` -- integer; a prime number

        OUTPUT: list of double precision complex numbers

        EXAMPLES::

            sage: n = numerical_eigenforms(11,4)
            sage: n.ap(2) # random order
            [9.0, 9.0, 2.73205080757, -0.732050807569]
            sage: n.ap(3) # random order
            [28.0, 28.0, -7.92820323028, 5.92820323028]
            sage: m = n.modular_symbols()
            sage: x = polygen(QQ, 'x')
            sage: m.T(2).charpoly('x').factor()
            (x - 9)^2 * (x^2 - 2*x - 2)
            sage: m.T(3).charpoly('x').factor()
            (x - 28)^2 * (x^2 + 2*x - 47)
        """
        p = Integer(p)
        if not p.is_prime():
            raise ValueError("p must be a prime")
        return Sequence(self.eigenvalues([p])[0], immutable=True)

    def eigenvalues(self, primes):
        """
        Return the eigenvalues of the Hecke operators corresponding
        to the primes in the input list of primes.   The eigenvalues
        match up between one prime and the next.

        INPUT:

        - ``primes`` -- list of primes

        OUTPUT: list of lists of eigenvalues

        EXAMPLES::

            sage: n = numerical_eigenforms(1,12)
            sage: n.eigenvalues([3,5,13])  # rel tol 2.4e-10
            [[177148.0, 252.00000000001896], [48828126.0, 4830.000000001376], [1792160394038.0, -577737.9999898539]]
        """
        primes = [Integer(p) for p in primes]
        for p in primes:
            if not p.is_prime():
                raise ValueError('each element of primes must be prime.')
        phi_x, phi_x_inv, nzp, x_nzp = self._eigendata()
        B = self._eigenvectors()

        def phi(y):
            """
            Take coefficients and a basis, and return that
            linear combination of basis vectors.

            EXAMPLES::

                sage: n = numerical_eigenforms(1,12)  # indirect doctest
                sage: n.eigenvalues([3,5,13])  # rel tol 2.4e-10
                [[177148.0, 252.00000000001896], [48828126.0, 4830.000000001376], [1792160394038.0, -577737.9999898539]]
            """
            return y.element() * B

        ans = []
        m = self.modular_symbols().ambient_module()
        for p in primes:
            t = m._compute_hecke_matrix_prime(p, nzp)
            w = phi(x_nzp*t)
            ans.append([w[i]*phi_x_inv[i] for i in range(w.degree())])
        return ans

    def systems_of_eigenvalues(self, bound):
        """
        Return all systems of eigenvalues for ``self`` for primes
        up to bound.

        EXAMPLES::

            sage: numerical_eigenforms(61).systems_of_eigenvalues(10)  # rel tol 1e-9
            [
            [-1.4811943040920152, 0.8060634335253695, 3.1563251746586642, 0.6751308705666477],
            [-1.0, -2.0000000000000027, -3.000000000000003, 1.0000000000000044],
            [0.3111078174659775, 2.903211925911551, -2.525427560843529, -3.214319743377552],
            [2.170086486626034, -1.7092753594369208, -1.63089761381512, -0.46081112718908984],
            [3.0, 4.0, 6.0, 8.0]
            ]
        """
        P = prime_range(bound)
        e = self.eigenvalues(P)
        v = Sequence([], cr=True)
        if len(e) == 0:
            return v
        for i in range(len(e[0])):
            v.append([e[j][i] for j in range(len(e))])
        v.sort()
        v.set_immutable()
        return v

    def systems_of_abs(self, bound):
        """
        Return the absolute values of all systems of eigenvalues for
        ``self`` for primes up to bound.

        EXAMPLES::

            sage: numerical_eigenforms(61).systems_of_abs(10)  # rel tol 1e-9
            [
            [0.3111078174659775, 2.903211925911551, 2.525427560843529, 3.214319743377552],
            [1.0, 2.0000000000000027, 3.000000000000003, 1.0000000000000044],
            [1.4811943040920152, 0.8060634335253695, 3.1563251746586642, 0.6751308705666477],
            [2.170086486626034, 1.7092753594369208, 1.63089761381512, 0.46081112718908984],
            [3.0, 4.0, 6.0, 8.0]
            ]
        """
        P = prime_range(bound)
        e = self.eigenvalues(P)
        v = Sequence([], cr=True)
        if len(e) == 0:
            return v
        for i in range(len(e[0])):
            v.append([abs(e[j][i]) for j in range(len(e))])
        v.sort()
        v.set_immutable()
        return v


def support(v, eps):
    """
    Given a vector `v` and a threshold eps, return all
    indices where `|v|` is larger than eps.

    EXAMPLES::

        sage: sage.modular.modform.numerical.support( numerical_eigenforms(61)._easy_vector(), 1.0 )
        []

        sage: sage.modular.modform.numerical.support( numerical_eigenforms(61)._easy_vector(), 0.5 )
        [0, 4]
    """
    return [i for i in range(v.degree()) if abs(v[i]) > eps]
