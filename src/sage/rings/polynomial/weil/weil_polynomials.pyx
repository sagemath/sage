## To enable OpenMP support, move the two lines starting with #distutils to the top.
# For best results, you may need to set the environment variable OMP_NUM_THREADS to a
# suitably large value (e.g., export OMP_NUM_THREADS=250).
#distutils: libraries = gomp
#distutils: extra_compile_args = -fopenmp
r"""
Iterator for Weil polynomials.

For `q` a prime power, a `q`-Weil polynomial is a monic polynomial with integer
coefficients whose complex roots all have absolute value `\sqrt{q}`. The class
``WeilPolynomials`` provides an iterable over a space of polynomials of this type.
It is possible to relax the monic condition by specifying one (or more) leading
coefficients. One may also impose certain congruence conditions. This can be
used to limit the Newton polygons of the resulting polynomials, or to lift
a polynomial specified by a congruence to a Weil polynomial.

For large jobs, one can set ``parallel=True`` to use ``OpenMP`` (if support was
enabled at compile time). Due to increased overhead, this is not recommended
for smaller problem sizes. To enable support, ensure that your compiler supports
``OpenMP`` and remove the appropriate ``#`` characters in the ``distutils`` commands below.
(You may also need to move those lines to the start of the file.)

AUTHOR:

-- Kiran S. Kedlaya (2007-05-28): initial version
-- (2015-08-29): switch from NTL to FLINT
-- (2017-10-03): consolidate Sage layer into .pyx file
                 define WeilPolynomials iterator
                 reverse convention for polynomials
                 pass multiprecision integers to/from C
-- (2019-02-02): update for Python3
                 improve parallel mode
-- (2019-12-19): final packaging for Sage (with help from David Roe)
-- (2026-02-10): extensive low-level optimization and reorganization
                 move some input sanitization out of C layer
                 improved algorithm for computing ranges
                 (using truncated Hausdorff moment condition)
                 more efficient application of Rolle condition (binary search)
                 choose number of processes based on available cores
                 import fmpz's directly as longs when possible
                 compute Hankel determinants via Dodgson condensation when possible
                 compute subresultants via Ducos method
                 more balanced work-splitting
                 more direct conversion of output to Sage polynomials
                 improved chunking in parallel mode
-- (2026-08-26): better answer conversion from FLINT to Sage
                 compute some Hankel determinants recursively
                 option to filter by linear constraints on Frobenius traces

A standalone version of this code can be found at
https://github.com/kedlaya/root-unitary
"""

# ****************************************************************************
#       Copyright (C) 2019-2026 Kiran S. Kedlaya <kskedl@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cython.parallel cimport prange
from cysignals.signals cimport sig_check
from cpython.mem cimport PyMem_Malloc, PyMem_Free

from sage.arith.misc import next_prime, primitive_root
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.functions.generalized import sgn

from sage.libs.gmp.mpz cimport mpz_set
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_vec cimport *
from sage.libs.flint.fmpz_poly cimport *
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint
from sage.structure.parent cimport Parent

cdef extern from "power_sums.c":
    ctypedef struct ps_static_data_t:
        pass

    ctypedef struct ps_dynamic_data_t:
        int flag        # State of the iterator (0 = inactive, 1 = running, 2 = found a solution, -1 = too many nodes)
        long node_count # Number of terminal nodes encountered
        fmpz *sympol    # Return value (a polynomial)

    int num_threads()
    int is_mpz(fmpz f)

    ps_static_data_t *ps_static_init(int d, const fmpz_t q, const fmpz_t lead, const fmpz *modlist, int num_constraints, const fmpz *constraints, long node_limit, int force_squarefree)
    ps_dynamic_data_t *ps_dynamic_init(int d, fmpz *coefflist)
    void ps_static_clear(ps_static_data_t *st_data)
    void ps_dynamic_clear(ps_dynamic_data_t *dy_data)
    void ps_cleanup(int n)

    int ps_dynamic_split(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, ps_dynamic_data_t *dy_data2) nogil
    int ps_next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int max_steps) nogil

cdef class dfs_manager:
    """
    Data structure to manage depth-first search.

    Such a structure is created and managed by an instance of ``WeilPolynomials_iter``.
    There is generally no need for a user to manipulate it directly.
    """
    cdef int d
    cdef int num_processes
    cdef long node_limit
    cdef ps_static_data_t *st_data
    cdef ps_dynamic_data_t **dy_data_buf
    cdef Parent ring

    def __cinit__(self, int d, q, coefflist, modlist, constraints,
                  long node_limit, int parallel, int force_squarefree, Parent ring):
        """
        Perform required C-level initialization (e.g., memory allocation).

        This is called automatically at object creation. It should never be called directly.

        TESTS:

        For some reason, Sage requires a dummy doctest to meet coverage requirements::

            sage: from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
            sage: w = WeilPolynomials(10, 1, sign=1, lead=[3,1,1])
            sage: it = iter(w)
            sage: it.process is not None # Verify object creation
            True
        """
        cdef fmpz_t temp_lead
        cdef fmpz_t temp_q
        cdef fmpz *temp_array
        cdef fmpz *temp_constraints
        cdef int i, j, num_constraints
        cdef Integer np

        self.d = d
        if parallel:
            np = (Integer(num_threads())**2).next_prime()
            while primitive_root(np) != 2:
                np = np.next_prime()
        else:
            np = Integer(1)
        self.num_processes = np
        self.ring = ring
        self.dy_data_buf = <ps_dynamic_data_t **>PyMem_Malloc(np * sizeof(ps_dynamic_data_t *))
        if not self.dy_data_buf:
            raise MemoryError()
        self.node_limit = node_limit
        fmpz_init(temp_lead)
        fmpz_set_mpz(temp_lead, Integer(coefflist[-1]).value)
        fmpz_init(temp_q)
        fmpz_set_mpz(temp_q, Integer(q).value)
        temp_array = _fmpz_vec_init(d+1)
        for i in range(d+1):
            fmpz_set_mpz(temp_array+i, Integer(modlist[i]).value)
        num_constraints = len(constraints)
        temp_constraints = _fmpz_vec_init((d+1) * num_constraints)
        for i in range(num_constraints):
            for j in range(min(d+1, len(constraints[i]))):
                fmpz_set_mpz(temp_constraints+(d+1)*i+j, Integer(constraints[i][j]).value)
        self.st_data = ps_static_init(d, temp_q, temp_lead,
                                         temp_array, num_constraints, temp_constraints,
                                         node_limit, force_squarefree)

        # Initialize processes, but assign work to only one process.
        # In parallel mode, other processes will get initialized later via work sharing.
        for i in range(d+1):
            fmpz_set_mpz(temp_array+i, Integer(coefflist[i]).value)
        self.dy_data_buf[0] = ps_dynamic_init(d, temp_array)
        for i in range(1, np):
            self.dy_data_buf[i] = ps_dynamic_init(d, NULL)

        fmpz_clear(temp_lead)
        fmpz_clear(temp_q)
        _fmpz_vec_clear(temp_array, d+1)
        _fmpz_vec_clear(temp_constraints, (d+1) * num_constraints)

    def __dealloc__(self):
        """
        Deallocate memory.
        """
        ps_static_clear(self.st_data)
        self.st_data = NULL
        if self.dy_data_buf != NULL:
            for i in range(self.num_processes):
                ps_dynamic_clear(self.dy_data_buf[i])
            PyMem_Free(self.dy_data_buf)
            self.dy_data_buf = NULL

    cpdef long node_count(self) noexcept:
        """
        Count nodes. This always gives 0 unless node_limit is specified.

        This method should not be called directly. Instead, use the ``node_count`` method
        of an instance of ``WeilPolynomials`` or ``WeilPolynomials_iter``.

        TESTS::

            sage: from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
            sage: w = WeilPolynomials(10, 1, sign=1, lead=[3,1,1], node_limit=10000000)
            sage: it = iter(w)
            sage: _ = next(it)
            sage: it.process.node_count()
            118
        """
        cdef long count = 0
        cdef int i
        for i in range(self.num_processes):
            count += self.dy_data_buf[i].node_count
        return count

    cdef Polynomial_integer_dense_flint extract_poly(self, int i):
        """
        Extract a polynomial computed by process i.
        """
        cdef int j, n = 2*self.d + 1
        cdef Polynomial_integer_dense_flint poly = self.ring.zero()
        cdef fmpz *sympol = self.dy_data_buf[i].sympol

        poly = poly._new() # Allocate a new polynomial in self.ring
        fmpz_poly_realloc(poly._poly, n)
        for j in range(n):
           fmpz_poly_set_coeff_fmpz(poly._poly, j, sympol+j)
        return poly

    cpdef object advance_exhaust(self, list ans, long ans_max=10000):
        """
        Advance the tree exhaustion.

        This method should not be called directly. Instead, use the iterator
        ``WeilPolynomials_iter`` or the iterable ``WeilPolynomials``.

        TESTS::

            sage: from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
            sage: w = WeilPolynomials(10,1,sign=1,lead=[3,1,1])
            sage: it = iter(w)
            sage: ans = []
            sage: ans = it.process.advance_exhaust(ans)
            sage: ans[0]
            3*x^10 + x^9 + x^8 - 5*x^7 + x^6 - 2*x^5 + x^4 - 5*x^3 + x^2 + x + 3
        """
        cdef int i, k = 1, flag
        cdef int t = 1, u = 0, np = self.num_processes, max_steps = 10000
        cdef long ans_count = 0

        if np == 1: # Serial mode
            while (t and not u and ans_count < ans_max):
                sig_check() # Check for interrupts
                t = ps_next_pol(self.st_data, self.dy_data_buf[0], max_steps)
                if t == 2: # Extract a solution
                    ans_count += 1
                    ans.append(self.extract_poly(0))
                elif t == -1: u = 1
        else: # Parallel mode
            while (t and not u and ans_count < ans_max):
                sig_check() # Check for interrupts
                t = 0
                for i in prange(np, nogil=True): # Step each process forward in parallel
                    t += ps_dynamic_split(self.st_data, self.dy_data_buf[i], self.dy_data_buf[(i+k) % np]) # Yield work
                    flag = ps_next_pol(self.st_data, self.dy_data_buf[i], max_steps)
                    t += flag
                    if flag == 2:
                        ans_count += 1
                        with gil: ans.append(self.extract_poly(i))
                    elif flag == -1: u += 1
                k = (k<<1) % np # Note that 2 is a primitive root mod np.
        if u:
            raise RuntimeError("Node limit ({0:%d}) exceeded".format(self.node_limit))
        return ans

class WeilPolynomials_iter():
    r"""
    Iterator created by WeilPolynomials.

    EXAMPLES::

        sage: from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
        sage: w = WeilPolynomials(10, 1, sign=1, lead=[3,1,1])
        sage: it = iter(w)
        sage: next(it)
        3*x^10 + x^9 + x^8 + 7*x^7 + 5*x^6 + 2*x^5 + 5*x^4 + 7*x^3 + x^2 + x + 3
        sage: w = WeilPolynomials(10, 1, sign=-1, lead=[3,1,1])
        sage: it = iter(w)
        sage: next(it)
        3*x^10 + x^9 + x^8 + 6*x^7 - 2*x^6 + 2*x^4 - 6*x^3 - x^2 - x - 3
    """
    def __init__(self, d, q, sign, lead, node_limit, parallel, squarefree, constraints=[], polring=None):
        r"""
        Create an iterator for Weil polynomials.

        EXAMPLES::

            sage: from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
            sage: w = WeilPolynomials(10, 1, sign=1, lead=[3,1,1])
            sage: it = iter(w)
            sage: next(it)
            3*x^10 + x^9 + x^8 + 7*x^7 + 5*x^6 + 2*x^5 + 5*x^4 + 7*x^3 + x^2 + x + 3
        """
        self.pol_internal = PolynomialRing(ZZ, name='x')
        self.pol = polring
        x = self.pol_internal.gen()
        d = Integer(d)
        if sign != 1 and sign != -1:
            raise ValueError("Invalid sign")
        if not q.is_integer() or q <= 0:
            raise ValueError("q must be a positive integer")
        if d == 0 and sign == -1: # No results
            self.process = None
            return
        if d % 2 == 0:
            if sign == 1:
                d2 = d//2
                num_cofactor = 0
                self.cofactor = 1
            else:
                d2 = d//2 - 1
                num_cofactor = 3
                self.cofactor = x**2 - q
        else:
            if not q.is_square():
                raise ValueError("Degree must be even if q is not a square")
            d2 = d//2
            if sign == 1:
                num_cofactor = 1
                self.cofactor = x + q.sqrt()
            else:
                num_cofactor = 2
                self.cofactor = x - q.sqrt()
        try:
            leadlist = list(lead)
        except TypeError:
            leadlist = [(lead, 0)]
        coefflist = []
        modlist = []
        for i in leadlist:
            try:
                (j, k) = i
            except TypeError:
                (j, k) = (i, 0)
            j = Integer(j)
            k = Integer(k)
            if len(modlist) == 0 and k != 0:
                raise ValueError("Leading coefficient must be specified exactly")
            if len(modlist) > 0 and ((k != 0 and modlist[-1]%k != 0) or (k == 0 and modlist[-1] != 0)):
                raise ValueError("Invalid moduli")
            coefflist.append(j)
            modlist.append(k)
        # Remove cofactor from initial coefficients
        if num_cofactor == 1: #cofactor x + sqrt(q)
            for i in range(1, len(coefflist)):
                coefflist[i] -= coefflist[i-1] * q.sqrt()
        elif num_cofactor == 2: #cofactor x - sqrt(q)
            for i in range(1, len(coefflist)):
                coefflist[i] += coefflist[i-1] * q.sqrt()
        elif num_cofactor == 3: #cofactor x^2 - q
            for i in range(2, len(coefflist)):
                coefflist[i] += coefflist[i-2]*q
        # Asymmetrize initial coefficients
        for i in range(len(coefflist)):
            for j in range(1, (len(coefflist)-i+1)//2):
                coefflist[i+2*j] -= (d2-i).binomial(j) * (q**j) * coefflist[i]
        for _ in range(d2+1-len(coefflist)):
            coefflist.append(0)
            modlist.append(1)
        coeffsign = sgn(coefflist[0])
        self.cofactor *= coeffsign
        coefflist = [x*coeffsign for x in reversed(coefflist)]
        if node_limit is None:
            node_limit = -1
        force_squarefree = Integer(squarefree)
        self.process = dfs_manager(d2, q, coefflist, modlist, constraints, node_limit,
                                   parallel, force_squarefree, self.pol_internal)
        self.q = q
        self.squarefree = squarefree
        self.ans = []
        self.num_cofactor = num_cofactor

    def __iter__(self):
        r"""
        Return the iterator (i.e. ``self``).

        EXAMPLES::

            sage: from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
            sage: w = WeilPolynomials(10, 1, sign=1, lead=[3,1,1])
            sage: it = iter(w)
            sage: it.__iter__() is it
            True
        """
        return self

    def __next__(self):
        r"""
        Step the iterator forward.

        EXAMPLES::

            sage: from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
            sage: w = WeilPolynomials(10, 1, sign=1, lead=[3,1,1])
            sage: it = iter(w)
            sage: next(it)
            3*x^10 + x^9 + x^8 + 7*x^7 + 5*x^6 + 2*x^5 + 5*x^4 + 7*x^3 + x^2 + x + 3
        """
        if self.process is None:
            raise StopIteration
        if len(self.ans) == 0:
            self.process.advance_exhaust(self.ans)
            if len(self.ans) == 0: # Iterator exhausted
                self.count = self.process.node_count()
                self.process = None
                raise StopIteration
        res = self.ans.pop()
        if self.pol is not None:
            res = self.pol(res)
        if self.num_cofactor: res *= self.cofactor
        return res

    def next(self): # For Python2 backward compatibility
        r"""
        Step the iterator forward.

        Included for Python2 backward compatibility.

        EXAMPLES::

            sage: from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
            sage: w = WeilPolynomials(10, 1, sign=1, lead=[3,1,1])
            sage: it = iter(w)
            sage: next(it)
            3*x^10 + x^9 + x^8 + 7*x^7 + 5*x^6 + 2*x^5 + 5*x^4 + 7*x^3 + x^2 + x + 3
        """
        return self.__next__()

    def node_count(self):
        r"""
        Return the number of terminal nodes found in the tree, excluding actual solutions.
        
        This always gives 0 unless node_limit is specified.

        EXAMPLES::

            sage: from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
            sage: w = WeilPolynomials(10, 1, sign=1, lead=[3,1,1], node_limit=10000000)
            sage: it = iter(w)
            sage: l = list(it)
            sage: it.node_count()
            118
        """
        if self.process is None:
            return self.count
        return self.process.node_count()


class WeilPolynomials():
    r"""
    Iterable for Weil polynomials, i.e., integer polynomials with all complex
    roots having a particular absolute value.

    Such polynomials `f` satisfy a functional equation

    .. MATH::

        T^d f(q/T) = s q^{d/2} f(T)

    where `d` is the degree of `f`, `s` is a sign and `q^{1/2}` is the absolute value
    of the roots of `f`.

    If parallel is False, then the order of values is descending lexicographical
    (i.e., polynomials with the largest coefficients of largest degrees sort first).

    If parallel is True, then the order of values is not specified. Beware that
    due to increased overhead, parallel execution may not yield a significant
    speedup for small problem sizes.

    INPUT:

    - ``d`` -- integer; the degree of the polynomials

    - ``q`` -- integer; the square of the complex absolute value of the roots

    - ``sign`` -- integer (default: `1`); the sign `s` of the functional equation

    - ``lead`` -- integer (default: `1`); list of integers or pairs of integers

        These are constraints on the leading coefficients of the generated polynomials.
        If pairs `(a, b)` of integers are given, they are treated as a constraint
        of the form `\equiv a \pmod{b}`; the moduli must be in decreasing order by
        divisibility, and the modulus of the leading coefficient must be 0.

    - ``node_limit`` -- integer (default: ``None``)

        If set, imposes an upper bound on the number of terminal nodes during the search
        (will raise a :exc:`RuntimeError` if exceeded).

    - ``parallel`` -- boolean (default: ``False``); whether to use multiple processes

        If set, will raise an error unless this file was compiled with OpenMP support
        (see instructions at the top of :mod:`sage.rings.polynomial.weil.weil_polynomials`).

    - ``squarefree`` -- boolean (default: ``False``)

        If set, only squarefree polynomials coprime to ``x^2 - q`` will be returned.

    - ``constraints`` -- list (default: empty)

        Each entry is a tuple expressing a lower bound constraint on the power sums `s_i`.
        The tuple `(c0, c1, ..., cn)` represents ``c_0 \leq c_1 s_1 + \cdots + c_n t_n``.

    - ``polring`` -- (optional) a polynomial ring in which to construct the results

        If not specified, answers are created in a polynomial ring over `ZZ` with variable `x`.

    EXAMPLES:

    Some simple cases::

        sage: from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
        sage: list(WeilPolynomials(2, 2))
        [x^2 + 2*x + 2, x^2 + x + 2, x^2 + 2, x^2 - x + 2, x^2 - 2*x + 2]
        sage: l = list(WeilPolynomials(4, 2))
        sage: l[0], l[-1]
        (x^4 + 4*x^3 + 8*x^2 + 8*x + 4, x^4 - 4*x^3 + 8*x^2 - 8*x + 4)
        sage: l = list(WeilPolynomials(3, 1, sign=-1))
        sage: l[0], l[-1]
        (x^3 + x^2 - x - 1, x^3 - 3*x^2 + 3*x - 1)

    By Kronecker's theorem, a monic integer polynomial has all roots of absolute
    value 1 if and only if it is a product of cyclotomic polynomials. For such a
    product to have positive sign of the functional equation, the factors `x-1`
    and `x+1` must each occur with even multiplicity. This code confirms
    Kronecker's theorem for polynomials of degree 6::

        sage: P.<x> = PolynomialRing(ZZ)
        sage: d = 6
        sage: ans1 = list(WeilPolynomials(d, 1, 1))
        sage: ans1.sort()
        sage: l = [(x-1)^2, (x+1)^2] + [cyclotomic_polynomial(n,x)
        ....:     for n in range(3, 2*d*d) if euler_phi(n) <= d]
        sage: w = WeightedIntegerVectors(d, [i.degree() for i in l])
        sage: ans2 = [prod(l[i]^v[i] for i in range(len(l))) for v in w]
        sage: ans2.sort()
        sage: print(ans1 == ans2)
        True

    Generating Weil polynomials with prescribed initial coefficients::

        sage: w = WeilPolynomials(10, 1, sign=1, lead=[3,1,1])
        sage: it = iter(w)
        sage: next(it)
        3*x^10 + x^9 + x^8 + 7*x^7 + 5*x^6 + 2*x^5 + 5*x^4 + 7*x^3 + x^2 + x + 3
        sage: w = WeilPolynomials(10, 1, sign=-1, lead=[3,1,1])
        sage: it = iter(w)
        sage: next(it)
        3*x^10 + x^9 + x^8 + 6*x^7 - 2*x^6 + 2*x^4 - 6*x^3 - x^2 - x - 3
        sage: list(WeilPolynomials(10, 2, lead=[1, -3, 5, -5, 5, -5]))
        [x^10 - 3*x^9 + 5*x^8 - 5*x^7 + 5*x^6 - 5*x^5 + 10*x^4 - 20*x^3 + 40*x^2 - 48*x + 32]

    Generating Weil polynomials with linear constraints on power sums::

        sage: w = WeilPolynomials(10, 1, sign=1, lead=[3,1])
        sage: l = list(w)
        sage: w2 = WeilPolynomials(10, 1, sign=1, lead=[3,1], constraints=[(3,1,0,1)])
        sage: l2 = list(w2)
        sage: l3 = [u for u in l if sum((i^3 + i) * j for i,j in u.roots(CC)) >= 2.99]
        sage: l2 == l3
        True
        sage: w3 = WeilPolynomials(10, 1, sign=1, lead=[3,1], constraints=[(2,0,0,1), (-3,2,1,-1)])
        sage: l4 = list(w3)
        sage: l5 = [u for u in l if sum(i^3 * j for i,j in u.roots(CC)) >= 1.99]
        sage: l6 = [u for u in l5 if sum((i^3 - i^2 - 2*i) * j for i,j in u.roots(CC)) <= 2.99 ]
        sage: l4 == l6
        True

    TESTS:

    Test restriction of initial coefficients::

        sage: w1 = WeilPolynomials(10, 1, sign=1, lead=3)
        sage: l1 = list(w1)
        sage: w2 = WeilPolynomials(10, 1, sign=1, lead=[3,1,1])
        sage: l2 = list(w2)
        sage: l3 = [i for i in l1 if i[1] == 1 and i[2] == 1]
        sage: l2 == l3
        True

        sage: w = WeilPolynomials(4, 2, lead=((1,0), (2,2)))
        sage: l = list(w)
        sage: l[0], l[-1]
        (x^4 + 4*x^3 + 8*x^2 + 8*x + 4, x^4 - 4*x^3 + 8*x^2 - 8*x + 4)
        sage: sorted(list(set(i[3] for i in l)))
        [-4, -2, 0, 2, 4]

    Test restriction to squarefree polynomials::

        sage: for (d,q,sign) in ((6,2,1), (6,4,-1), (5,4,-1)):
        ....:     w1 = WeilPolynomials(d, q, sign=sign)
        ....:     l1 = list(w1)
        ....:     w2 = WeilPolynomials(d, q, sign=sign, squarefree=True)
        ....:     l2 = list(w2)
        ....:     l3 = [i for i in l1 if i.is_squarefree()]
        ....:     print(l2 == l3)
        True
        True
        True

    Test that :issue:`29475` is resolved::

        sage: P.<x> = QQ[]
        sage: u = x^6 + x^5 + 6*x^4 - 2*x^3 + 66*x^2 + 121*x + 1331
        sage: u.is_weil_polynomial()
        True
        sage: u in WeilPolynomials(6, 11, 1, [1,1,6])
        True
        sage: u in WeilPolynomials(6, 11, 1, [(1,0),(1,11),(6,11)])
        True

    Test that :issue:`31809` is resolved::

        sage: foo = list(WeilPolynomials(12, 3, lead=(1,0,9,2,46), squarefree=False))
        sage: bar = list(WeilPolynomials(12, 3, lead=(1,0,9,2,46), squarefree=True))
        sage: bar == [f for f in foo if f.is_squarefree()]
        True

    Test that :issue:`32348` is resolved::

        sage: list(WeilPolynomials(10, 2, lead=(1,-3,5,-5,5,-5)))
        [x^10 - 3*x^9 + 5*x^8 - 5*x^7 + 5*x^6 - 5*x^5 + 10*x^4 - 20*x^3 + 40*x^2 - 48*x + 32]

    Test that :issue:`37860` is resolved::

        sage: list(WeilPolynomials(0, 1, sign=-1))
        []
    """
    def __init__(self, d, q, sign=1, lead=1, node_limit=None, parallel=False, squarefree=False, constraints=[], polring=None):
        r"""
        Initialize this iterable.

        EXAMPLES::

            sage: from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
            sage: w = WeilPolynomials(10, 1, sign=1, lead=[3,1,1])
            sage: w.__init__(10, 1, sign=1, lead=[3,1,-1]) # Change parameters before iterating
            sage: it = iter(w)
            sage: next(it) # Results reflect the changed parameters
            3*x^10 + x^9 - x^8 + 7*x^7 + 5*x^6 - 2*x^5 + 5*x^4 + 7*x^3 - x^2 + x + 3
        """
        if parallel and num_threads() == 1:
            raise RuntimeError("Parallel execution not supported")
        self.data = (d, q, sign, lead, node_limit, parallel, squarefree, constraints, polring)

    def __iter__(self):
        r"""
        Construct the associated iterator.

        EXAMPLES::

            sage: from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
            sage: w = WeilPolynomials(10, 1, sign=1, lead=[3,1,1])
            sage: it = w.__iter__()
            sage: next(it)
            3*x^10 + x^9 + x^8 + 7*x^7 + 5*x^6 + 2*x^5 + 5*x^4 + 7*x^3 + x^2 + x + 3
        """
        w = WeilPolynomials_iter(*self.data)
        self.w = w
        return w

    def node_count(self):
        r"""
        Return the number of terminal nodes found in the tree, excluding actual solutions.
        
        This always gives 0 unless node_limit is specified.

        EXAMPLES::

            sage: from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
            sage: w = WeilPolynomials(10, 1, sign=1, lead=[3,1,1], node_limit=10000000)
            sage: l = list(w)
            sage: w.node_count()
            118
        """
        return self.w.node_count()

