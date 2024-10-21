# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
S-Boxes and Their Algebraic Representations
"""
cimport cython
from cysignals.memory cimport check_allocarray, sig_free

from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport Element

from sage.combinat.integer_vector import IntegerVectors
from sage.crypto.boolean_function import BooleanFunction
from sage.crypto.boolean_function cimport hamming_weight, walsh_hadamard
from sage.matrix.constructor import matrix
from sage.matrix.matrix0 cimport Matrix
from sage.misc.cachefunc import cached_method
from sage.misc.functional import is_even
from sage.misc.misc_c import prod as mul
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.ideal import FieldIdeal, Ideal
from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


cdef Py_ssize_t _nterms(Py_ssize_t nvars, Py_ssize_t deg) noexcept:
    """
    Return the number of monomials possible up to a given
    degree.

    INPUT:

    - ``nvars`` -- number of variables

    - ``deg`` -- degree

    TESTS::

        sage: from sage.crypto.sbox import SBox
        sage: S = SBox(7,6,0,4,2,5,1,3)
        sage: F = S.polynomials(degree=3)  # indirect doctest                           # needs sage.libs.singular
    """
    cdef Py_ssize_t total = 1
    cdef Py_ssize_t divisor = 1
    cdef Py_ssize_t var_choices = 1

    cdef Py_ssize_t d
    for d in range(1, deg+1):
        var_choices *= (nvars - d + 1)
        divisor *= d
        total += var_choices // divisor

    return total


@cython.auto_pickle(True)
cdef class SBox(SageObject):
    r"""
    A substitution box or S-box is one of the basic components of
    symmetric key cryptography. In general, an S-box takes `m` input
    bits and transforms them into `n` output bits. This is called an
    `m \times n` S-box and is often implemented as a lookup table. These
    S-boxes are carefully chosen to resist linear and differential
    cryptanalysis [He2002]_.

    This module implements an S-box class which allows an algebraic
    treatment and determine various cryptographic properties.

    EXAMPLES:

    We consider the S-box of the block cipher PRESENT [BKLPPRSV2007]_::

        sage: from sage.crypto.sbox import SBox
        sage: S = SBox(12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2); S
        (12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2)
        sage: S(1)
        5

    Note that by default bits are interpreted in big endian
    order. This is not consistent with the rest of Sage, which has a
    strong bias towards little endian, but is consistent with most
    cryptographic literature::

        sage: S([0,0,0,1])
        [0, 1, 0, 1]

        sage: S = SBox(12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2, big_endian=False)
        sage: S(1)
        5
        sage: S([0,0,0,1])
        [1, 1, 0, 0]


    Now we construct an ``SBox`` object for the 4-bit small scale AES
    S-Box (cf. :mod:`sage.crypto.mq.sr`)::

        sage: sr = mq.SR(1,1,1,4, allow_zero_inversions=True)
        sage: S = SBox([sr.sub_byte(e) for e in list(sr.k)])
        sage: S
        (6, 5, 2, 9, 4, 7, 3, 12, 14, 15, 10, 0, 8, 1, 13, 11)

    AUTHORS:

    - Rusydi H. Makarim (2016-03-31) : added more functions to determine related cryptographic properties
    - Yann Laigle-Chapuy (2009-07-01): improve linear and difference matrix computation
    - Martin R. Albrecht (2008-03-12): initial implementation

    REFERENCES:

    - [He2002]_
    - [BKLPPRSV2007]_
    - [CDL2015]_
    """
    cdef list _S_list
    cdef object _ring
    cdef Py_ssize_t m
    cdef Py_ssize_t n
    cdef bint _big_endian
    cdef dict __dict__  # for cached_methods

    def __init__(self, *args, **kwargs):
        r"""
        Construct a substitution box (S-box) for a given lookup table `S`.

        INPUT:

        - ``S`` -- a finite iterable defining the S-box with integer or
          finite field elements

        - ``big_endian`` -- boolean (default: ``True``); controls whether bits
          shall be ordered in big endian order

        EXAMPLES:

        We construct a 3-bit S-box where e.g. the bits (0,0,1) are
        mapped to (1,1,1).::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3); S
            (7, 6, 0, 4, 2, 5, 1, 3)

            sage: S(0)
            7

        Construct S-box from univariate polynomial.::

            sage: R = PolynomialRing(GF(2**3), 'x')
            sage: inv = R.gen()**(2**3-2)
            sage: inv = SBox(inv); inv
            (0, 1, 5, 6, 7, 2, 3, 4)
            sage: inv.differential_uniformity()
            2

            sage: SBox(PolynomialRing(GF(3**3), 'x').gen())
            Traceback (most recent call last):
            ...
            TypeError: only polynomials over rings with characteristic 2 allowed

        TESTS::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox()
            Traceback (most recent call last):
            ...
            TypeError: no lookup table provided
            sage: S = SBox(1, 2, 3)
            Traceback (most recent call last):
            ...
            TypeError: lookup table length is not a power of 2
            sage: S = SBox(5, 6, 0, 3, 4, 2, 1, 2)
            sage: S.output_size()
            3
        """
        from sage.rings.polynomial.polynomial_element import Polynomial

        if "S" in kwargs:
            args = kwargs["S"]

        if len(args) == 1 and isinstance(args[0], Polynomial):
            # SBox defined via Univariate Polynomial, compute lookup table
            # by evaluating the polynomial on every base_ring element
            poly = args[0]
            R = poly.parent().base_ring()
            if R.characteristic() != 2:
                raise TypeError("only polynomials over rings with characteristic 2 allowed")
            S = [poly(v) for v in sorted(R)]
        elif len(args) == 1:  # iterables
            S = args[0]
        elif len(args) > 1:
            S = args
        else:
            raise TypeError("no lookup table provided")

        _S_list = []
        for e in S:
            if isinstance(e, Element) and isinstance(e.parent(), FiniteField):
                e = e.polynomial().change_ring(ZZ).subs(e.parent().characteristic())
            _S_list.append(e)
        S = _S_list

        if not ZZ(len(S)).is_power_of(2):
            raise TypeError("lookup table length is not a power of 2")
        self._S_list = S

        self.m = ZZ(len(S)).exact_log(2)
        self.n = ZZ(max(S)).nbits()
        self._big_endian = kwargs.get("big_endian", True)

        cdef Py_ssize_t i
        self._ring = PolynomialRing(
            GF(2),
            self.m + self.n,
            ["x%d" % i for i in range(self.m)] + ["y%d" % i for i in range(self.n)])

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: SBox(7,6,0,4,2,5,1,3)  # indirect doctest
            (7, 6, 0, 4, 2, 5, 1, 3)
        """
        return "(" + ", ".join(map(str, self)) + ")"

    def __len__(self):
        """
        Return the length of input bit strings.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: len(SBox(7,6,0,4,2,5,1,3))
            3
        """
        return self.m

    def __eq__(self, rhs):
        """
        S-boxes are considered to be equal if all construction
        parameters match.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: loads(dumps(S)) == S
            True
        """
        if not isinstance(rhs, SBox):
            raise NotImplemented

        cdef SBox other = <SBox> rhs
        return (self._S_list == other._S_list) and (self._big_endian == self._big_endian)

    def __ne__(self, other):
        """
        S-boxes are considered to be equal if all construction
        parameters match.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: S != S
            False
        """
        return not self.__eq__(other)

    cpdef list to_bits(self, x, n=None):
        """
        Return bitstring of length ``n`` for integer ``x``. The
        returned bitstring is guaranteed to have length ``n``.

        INPUT:

        - ``x`` -- integer

        - ``n`` -- bit length (optional)

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: S.to_bits(6)
            [1, 1, 0]

            sage: S.to_bits( S(6) )
            [0, 0, 1]

            sage: S( S.to_bits( 6 ) )
            [0, 0, 1]
        """
        if n is None and self.m == self.n:
            n = self.n

        F = GF(2)
        cdef list xs = [F(i) for i in ZZ(x).digits(base=2, padto=n)]

        if self._big_endian:
            xs.reverse()

        return xs

    def from_bits(self, x, n=None):
        """
        Return integer for bitstring ``x`` of length ``n``.

        INPUT:

        - ``x`` -- a bitstring

        - ``n`` -- bit length (optional)

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: S.from_bits( [1,1,0])
            6

            sage: S( S.from_bits( [1,1,0] ) )
            1
            sage: S.from_bits( S( [1,1,0] ) )
            1
        """
        if n is None and self.m == self.n:
            n = self.n

        if self._big_endian:
            x = list(reversed(x))

        return ZZ(self._rpad(x, n), 2)

    cdef list _rpad(self, list x, Py_ssize_t n=-1):
        """
        Right pads ``x`` such that ``len(x) == n``.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: S._rpad([1,1])  # not tested
            [1, 1, 0]
        """
        if n == -1 and self.m == self.n:
            n = self.n
        return x + [GF(2).zero()] * (n-len(x))

    def __call__(self, X):
        r"""
        Apply substitution to ``X``.

        If ``X`` is a list, it is interpreted as a sequence of bits
        depending on the bit order of this S-box.

        INPUT:

        - ``X`` -- either an integer, a tuple of `\GF{2}` elements of
          length ``len(self)`` or a finite field element in
          `\GF{2^n}`. As a last resort this function tries to convert
          ``X`` to an integer.

        EXAMPLES:

        We can call SBoxes with integers as inputs, this will
        return an integer::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(3, 0, 1, 3, 1, 0, 2, 2)
            sage: S(0)
            3

            sage: S([0,0,0])
            [1, 1]

            sage: S = SBox([7,6,0,4,2,5,1,3])
            sage: S(7)
            3

            sage: S[0]
            7

            sage: S(QQ(3))
            4

        Alternatively, we can call in with a list-like object, which
        will return a list::

            sage: S((0,2,3))
            [0, 1, 1]

            sage: S((0,0,1))
            [1, 1, 0]

        Calling it with a vector will return a vector::

            sage: S(vector(GF(2), [0,0,1]))
            (1, 1, 0)

            sage: type(vector(GF(2), [0,0,1])) == type(S(vector(GF(2), [0,0,1])))
            True

        An input from a finite field will be interpreted as given the
        coefficient vector as input::

            sage: k.<a> = GF(2^3)
            sage: S(a^2)  # interpreted as (0,0,1)
            a + 1
            sage: S([0, 0, 1])
            [1, 1, 0]
            sage: vector(a + 1)
            (1, 1, 0)

            sage: id = SBox(range(8))
            sage: all([x == id(x) for x in k])
            True

        Some examples for inputs that throw an :exc:`TypeError`::

            sage: S([1]*10^6)
            Traceback (most recent call last):
            ...
            TypeError: cannot apply SBox to provided element

            sage: S(1/2)
            Traceback (most recent call last):
            ...
            TypeError: cannot apply SBox to 1/2
        """
        # Handle integer inputs
        if isinstance(X, int):
            return self._S_list[<int> X]
        if isinstance(X, Integer):
            return self._S_list[<Integer> X]

        # Handle non-integer inputs: vectors, finite field elements to-integer-coercible elements
        # cdef int i
        if isinstance(X, Element):
            K = X.parent()
            if K.base_ring().characteristic() != 2:
                try:
                    X = ZZ(X)
                    return K(self._S_list[<Integer> X])
                except TypeError:
                    raise TypeError("cannot apply SBox to %s" % (X,))
            V = None
            try:
                V = K.vector_space(map=False)
            except AttributeError:
                try:
                    return self._S_list[ZZ(X)]
                except TypeError:
                    pass
            except TypeError:
                V = K.vector_space()
            # convert finite field element to vector
            if V is not None:
                X = V(X)
        else:
            K = None

        # At this point, we only handle X being a vector or list-like
        cdef list out
        cdef Py_ssize_t ell
        try:
            ell = len(X)
            X = list(X)
        except TypeError:
            # This will not pass the next test and ultimately raise an error
            ell = -1

        if ell == self.m:
            if self._big_endian:
                X = list(reversed(X))
            X = ZZ(X, 2)
            out = self.to_bits(self._S_list[X], self.n)
            if K is not None:
                return K(out)
            # NOTE: Parts of the code assume that when a list is passed
            #   in that a list is returned
            return out

        if len(str(X)) > 50:
            raise TypeError("cannot apply SBox to provided element")
        else:
            raise TypeError("cannot apply SBox to %s" % (X,))

    def __getitem__(self, X):
        """
        See :meth:`SBox.__call__`.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([7,6,0,4,2,5,1,3])
            sage: S[7]
            3
        """
        return self(X)

    def input_size(self):
        """
        Return the input size of this S-Box.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([0, 3, 2, 1, 1, 3, 2, 0])
            sage: S.input_size()
            3
        """
        return self.m

    def output_size(self):
        """
        Return the output size of this S-Box.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([0, 3, 2, 1, 1, 3, 2, 0])
            sage: S.output_size()
            2
        """
        return self.n

    def is_permutation(self):
        r"""
        Return ``True`` if this S-Box is a permutation.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: S.is_permutation()
            True

            sage: S = SBox(3,2,0,0,2,1,1,3)
            sage: S.is_permutation()
            False
        """
        if self.m != self.n:
            return False
        cdef Py_ssize_t m = self.m
        cdef Py_ssize_t i
        return len(set([self._S_list[i] for i in range(1 << m)])) == 1 << m

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: [e for e in S]
            [7, 6, 0, 4, 2, 5, 1, 3]
        """
        cdef Py_ssize_t i
        for i in range(1 << self.m):
            yield self._S_list[i]

    def derivative(self, u):
        r"""
        Return the derivative in direction of ``u``.

        INPUT:

        - ``u`` -- either an integer or a tuple/list of `\GF{2}` elements
          of length equal to ``m``


        The derivative of `F` in direction of `u` is defined as
        `x \mapsto F(x) + F(x + u)`.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: s = SBox(0,1,2,3)
            sage: s.derivative(1)
            (1, 1, 1, 1)
            sage: u = [1,0]
            sage: s.derivative(u)
            (1, 1, 1, 1)
            sage: v = vector(GF(2), [1,0])
            sage: s.derivative(v)
            (1, 1, 1, 1)
            sage: s.derivative(4)
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: from sage.crypto.sboxes import PRESENT
            sage: PRESENT.derivative(1).max_degree() < PRESENT.max_degree()             # needs sage.rings.polynomial.pbori
            True
        """
        from sage.structure.element import Vector
        nvars = self.m

        if isinstance(u, (tuple, list)):
            v = ZZ(u, base=2)
        elif isinstance(u, Vector):
            if u.base_ring() != GF(2):
                raise TypeError("base ring of input vector must be GF(2)")
            elif u.parent().dimension() != nvars:
                raise TypeError("input vector must be an element of a vector space with dimension %d" % (nvars,))
            v = ZZ(u.list(), base=2)
        else:
            v = u

        return SBox([self(x) ^ self(x ^ v)
                     for x in range(1 << self.input_size())])

    @cached_method
    def difference_distribution_table(self):
        """
        Return difference distribution table (DDT) ``A`` for this S-box.

        The rows of ``A`` encode the differences ``Delta I`` of the
        input and the columns encode the difference ``Delta O`` for
        the output. The bits are ordered according to the endianess of
        this S-box. The value at ``A[Delta I,Delta O]`` encodes how
        often ``Delta O`` is the actual output difference given
        ``Delta I`` as input difference.

        See [He2002]_ for an introduction to differential
        cryptanalysis.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: S.difference_distribution_table()
            [8 0 0 0 0 0 0 0]
            [0 2 2 0 2 0 0 2]
            [0 0 2 2 0 0 2 2]
            [0 2 0 2 2 0 2 0]
            [0 2 0 2 0 2 0 2]
            [0 0 2 2 2 2 0 0]
            [0 2 2 0 0 2 2 0]
            [0 0 0 0 2 2 2 2]
            sage: S = SBox(7,4,8,6)
            sage: S.difference_distribution_table()
            [4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 2 0 0 0 0 0 0 0 0 0 0 2 0]
            [0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 2]
            [0 2 0 0 0 0 0 0 0 0 0 0 2 0 0 0]

        TESTS::

        Testing square SBoxes::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: S.difference_distribution_table()
            [8 0 0 0 0 0 0 0]
            [0 2 2 0 2 0 0 2]
            [0 0 2 2 0 0 2 2]
            [0 2 0 2 2 0 2 0]
            [0 2 0 2 0 2 0 2]
            [0 0 2 2 2 2 0 0]
            [0 2 2 0 0 2 2 0]
            [0 0 0 0 2 2 2 2]

        Testing non-square SBoxes::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(8,8,8,8)
            sage: S.difference_distribution_table()
            [4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            sage: S = SBox(7,4,8,6)
            sage: S.difference_distribution_table()
            [4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 2 0 0 0 0 0 0 0 0 0 0 2 0]
            [0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 2]
            [0 2 0 0 0 0 0 0 0 0 0 0 2 0 0 0]
            sage: S = SBox(0,0,0,1,0,0,1,3)
            sage: S.difference_distribution_table()
            [8 0 0 0]
            [4 2 2 0]
            [2 4 0 2]
            [2 4 0 2]
            [4 2 2 0]
            [6 0 0 2]
            [2 4 0 2]
            [2 4 0 2]
        """
        cdef Py_ssize_t nrows = 1 << self.m
        cdef Py_ssize_t ncols = 1 << self.n
        cdef Py_ssize_t i, di

        cdef list L = [0]*(nrows*ncols)

        for i in range(nrows):
            si = self._S_list[i]
            for di in range(nrows):
                L[di*ncols + si ^ self._S_list[i ^ di]] += 1

        A = matrix(ZZ, nrows, ncols, L)
        A.set_immutable()

        return A

    def maximal_difference_probability_absolute(self):
        """
        Return the difference probability of the difference with the
        highest probability in absolute terms, i.e. how often it
        occurs in total.

        Equivalently, this is equal to the differential uniformity
        of this S-Box.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: S.maximal_difference_probability_absolute()
            2

        .. NOTE::

            This code is mainly called internally.
        """
        A = self.difference_distribution_table().__copy__()
        A[0, 0] = 0
        return max(map(abs, A.list()))

    differential_uniformity = maximal_difference_probability_absolute

    def maximal_difference_probability(self):
        r"""
        Return the difference probability of the difference with the
        highest probability in the range between 0.0 and 1.0
        indicating 0\% or 100\% respectively.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: S.maximal_difference_probability()
            0.25
        """
        return self.maximal_difference_probability_absolute() / (2.0**self.output_size())

    @cached_method
    def linear_approximation_table(self, scale='absolute_bias'):
        r"""
        Return linear approximation table (LAT) `A` for this S-box.

        The entry `A[\alpha,\beta]` corresponds to the probability
        `Pr[\alpha\cdot x = \beta\cdot S(x)]`, where `S` is this S-box
        mapping `n`-bit inputs to `m`-bit outputs.
        There are three typical notations for this probability used in
        the literature:

        - `Pr[\alpha\cdot x = \beta\cdot S(x)] = 1/2 + e(\alpha, \beta)`,
          where `e(\alpha, \beta)` is called the bias,
        - `2\cdot Pr[\alpha\cdot x = \beta\cdot S(x)] = 1 + c(\alpha, \beta)`,
          where `c(\alpha, \beta) = 2\cdot e(\alpha, \beta)` is the
          correlation, and
        - `2^{(m+1)}\cdot Pr[\alpha\cdot x = \beta\cdot S(x)] = 2^m +
          \hat{S}(\alpha, \beta)`, where `\hat{S}(\alpha, \beta)` is
          the Fourier coefficient of S.

        See [He2002]_ for an introduction to linear cryptanalysis.

        INPUT:

        - ``scale`` -- string to choose the scaling for the LAT, one of

          * "bias": elements are `e(\alpha, \beta)`
          * "correlation": elements are `c(\alpha, \beta)`
          * "absolute_bias": elements are `2^m\cdot e(\alpha, \beta)` (default)
          * "fourier_coefficient": elements are `\hat{S}(\alpha, \beta)`

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: lat_abs_bias = S.linear_approximation_table()
            sage: lat_abs_bias
            [ 4  0  0  0  0  0  0  0]
            [ 0  0  0  0  2  2  2 -2]
            [ 0  0 -2 -2 -2  2  0  0]
            [ 0  0 -2  2  0  0 -2 -2]
            [ 0  2  0  2 -2  0  2  0]
            [ 0 -2  0  2  0  2  0  2]
            [ 0 -2 -2  0  0 -2  2  0]
            [ 0 -2  2  0 -2  0  0 -2]

            sage: lat_abs_bias/(1 << S.input_size()) == S.linear_approximation_table(scale='bias')
            True

            sage: lat_abs_bias/(1 << (S.input_size()-1)) == S.linear_approximation_table(scale='correlation')
            True

            sage: lat_abs_bias*2 == S.linear_approximation_table(scale='fourier_coefficient')
            True

        According to this table the first bit of the input is equal
        to the third bit of the output 6 out of 8 times::

            sage: for i in srange(8): print(S.to_bits(i)[0] == S.to_bits(S(i))[2])
            False
            True
            True
            True
            False
            True
            True
            True
        """
        cdef Py_ssize_t m = self.m
        cdef Py_ssize_t n = self.n

        cdef Py_ssize_t nrows = 1 << m
        cdef Py_ssize_t ncols = 1 << n

        # directly compute the walsh_hadamard transform here, without
        # creating the BooleanFunction object
        cdef long* temp = <long*> check_allocarray(nrows*ncols, sizeof(long))
        cdef Py_ssize_t i, j

        for i in range(ncols):
            for j in range(nrows):
                temp[i*nrows + j] = 1 - (<int>(hamming_weight(i & self._S_list[j]) & 1) << 1)
            walsh_hadamard(&temp[i*nrows], m)

        cdef list L = [temp[i*nrows + j] for j in range(nrows) for i in range(ncols)]
        sig_free(temp)

        A = matrix(ZZ, nrows, ncols, L)

        if (scale is None) or (scale == "absolute_bias"):
            A /= 2
        elif scale == "bias":
            A /= 1 << (m+1)
        elif scale == "correlation":
            A /= 1 << m
        elif scale == "fourier_coefficient":
            pass
        else:
            raise ValueError("no such scaling for the LAT: %s" % scale)

        A.set_immutable()
        return A

    def maximal_linear_bias_absolute(self):
        r"""
        Return maximal linear bias, i.e. how often the linear
        approximation with the highest bias is true or false minus
        `2^{n-1}`.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: S.maximal_linear_bias_absolute()
            2
        """
        cdef list L = self.linear_approximation_table().list()
        L[0] = 0  # Set the upper left corner to be 0
        return max(map(abs, L))

    def maximal_linear_bias_relative(self):
        """
        Return maximal bias of all linear approximations of this
        S-box.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: S.maximal_linear_bias_relative()
            0.25
        """
        return self.maximal_linear_bias_absolute() / (2.0**self.m)

    def ring(self):
        """
        Create, return, and cache a polynomial ring for S-box polynomials.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: S.ring()
            Multivariate Polynomial Ring in x0, x1, x2, y0, y1, y2 over
             Finite Field of size 2
        """
        return self._ring

    def solutions(self, X=None, Y=None):
        """
        Return a dictionary of solutions to this S-box.

        INPUT:

        - ``X`` -- (optional) input variables

        - ``Y`` -- (optional) output variables

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([7,6,0,4,2,5,1,3])
            sage: F = S.polynomials()                                                   # needs sage.libs.singular
            sage: s = S.solutions()
            sage: any(f.subs(_s) for f in F for _s in s)                                # needs sage.libs.singular
            False
        """
        if X is None and Y is None:
            P = self.ring()
            gens = P.gens()
        else:
            P = X[0].parent()
            gens = X + Y

        cdef Py_ssize_t m = self.m

        cdef list solution, solutions = []
        cdef Py_ssize_t i
        for i in range(1 << m):
            solution = self.to_bits(i, m) + <list> self(self.to_bits(i, m))
            solutions.append(dict(zip(gens, solution)))

        return solutions

    def polynomials(self, X=None, Y=None, degree=2, groebner=False):
        """
        Return a list of polynomials satisfying this S-box.

        First, a simple linear fitting is performed for the given
        ``degree`` (cf. for example [BC2003]_). If ``groebner=True`` a
        Groebner basis is also computed for the result of that
        process.

        INPUT:

        - ``X`` -- (optional) input variables

        - ``Y`` -- (optional) output variables

        - ``degree`` -- (default: ``2``) integer > 0

        - ``groebner`` -- boolean (default: ``False``); calculate a reduced
          Groebner basis of the spanning polynomials to obtain more polynomials

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: P = S.ring()

        By default, this method returns an indirect representation::

            sage: S.polynomials()                                                       # needs sage.libs.singular
            [x0*x2 + x1 + y1 + 1,
             x0*x1 + x1 + x2 + y0 + y1 + y2 + 1,
             x0*y1 + x0 + x2 + y0 + y2,
             x0*y0 + x0*y2 + x1 + x2 + y0 + y1 + y2 + 1,
             x1*x2 + x0 + x1 + x2 + y2 + 1,
             x0*y0 + x1*y0 + x0 + x2 + y1 + y2,
             x0*y0 + x1*y1 + x1 + y1 + 1,
             x1*y2 + x1 + x2 + y0 + y1 + y2 + 1,
             x0*y0 + x2*y0 + x1 + x2 + y1 + 1,
             x2*y1 + x0 + y1 + y2,
             x2*y2 + x1 + y1 + 1,
             y0*y1 + x0 + x2 + y0 + y1 + y2,
             y0*y2 + x1 + x2 + y0 + y1 + 1,
             y1*y2 + x2 + y0]

        We can get a direct representation by computing a
        lexicographical Groebner basis with respect to the right
        variable ordering, i.e. a variable ordering where the output
        bits are greater than the input bits::

            sage: P.<y0,y1,y2,x0,x1,x2> = PolynomialRing(GF(2),6,order='lex')
            sage: S.polynomials([x0,x1,x2],[y0,y1,y2], groebner=True)                   # needs sage.libs.singular
            [y0 + x0*x1 + x0*x2 + x0 + x1*x2 + x1 + 1,
             y1 + x0*x2 + x1 + 1,
             y2 + x0 + x1*x2 + x1 + x2 + 1]

        TESTS:

        Check that :issue:`22453` is fixed::

            sage: from sage.crypto.sboxes import AES
            sage: aes_polys = AES.polynomials()  # long time
            sage: aes_polys[3]                   # long time
            x3*y0 + x5*y0 + x7*y0 + x6*y1 + x2*y2 + x3*y2 + x4*y2
                  + x2*y3 + x3*y3 + x5*y4 + x6*y4 + x3*y5 + x4*y5
                  + x4*y7 + x2 + x3 + y2 + y3 + y4 + 1
        """
        cdef Py_ssize_t m = self.m
        cdef Py_ssize_t n = self.n

        if X is None and Y is None:
            P = self.ring()
            X = P.gens()[:m]
            Y = P.gens()[m:]
        else:
            P = X[0].parent()

        gens = X + Y

        cdef Py_ssize_t i
        cdef list bits = []
        for i in range(1 << m):
            bits.append(self.to_bits(i, m) + <list> self(self.to_bits(i, m)))

        cdef Py_ssize_t ncols = (1 << m) + 1

        A = matrix(P, _nterms(m + n, degree), ncols)

        cdef list exponents = []
        cdef Py_ssize_t d
        for d in range(degree+1):
            exponents += IntegerVectors(d, max_length=m+n, min_length=m+n,
                                        min_part=0, max_part=1).list()

        row = 0
        for exponent in exponents:
            A[row, ncols-1] = mul([gens[i]**exponent[i] for i in range(len(exponent))])
            for col in range(1 << m):
                A[row, col] = mul([bits[col][i] for i in range(len(exponent)) if exponent[i]])
            row += 1

        rankSize = A.rank() - 1

        cdef Py_ssize_t c
        for c in range(ncols):
            A[0, c] = 1

        RR = A.echelon_form(algorithm='row_reduction')

        # extract spanning stet
        gens = (RR.column(ncols-1)[rankSize:]).list()

        if not groebner:
            return gens

        FI = set(FieldIdeal(P).gens())
        I = Ideal(gens + list(FI))
        gb = I.groebner_basis()

        gens = []
        for f in gb:
            # filter out field equations
            if f not in FI:
                gens.append(f)
        return gens

    def interpolation_polynomial(self, k=None):
        r"""
        Return a univariate polynomial over an extension field
        representing this S-box.

        If ``m`` is the input length of this S-box then the extension
        field is of degree ``m``.

        If the output length does not match the input length then a
        :exc:`TypeError` is raised.

        INPUT:

        - ``k`` -- (optional) an instance of `\GF{2^m}`

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: f = S.interpolation_polynomial()
            sage: f
            (a^2 + a + 1)*x^6 + a^2*x^5 + (a + 1)*x^4 + (a^2 + a)*x^3
             + x^2 + a*x + a^2 + a + 1

            sage: a = f.base_ring().gen()

            sage: f(0), S(0)
            (a^2 + a + 1, 7)

            sage: f(a^2 + 1), S(5)
            (a^2 + 1, 5)

        .. NOTE::

            The method-internal call to the S-box initially used a different
            endianess for handling finite field elements. This changed in
            :issue:`25633`, by calling the S-box directly.
        """
        if self.m != self.n:
            raise TypeError("Lagrange interpolation only supported if"
                            "self.input_size() == self.output_size()")

        cdef Py_ssize_t m = self.m
        if k is None:
            k = GF(2**m, 'a')

        cdef list l = []
        cdef int i
        for i in range(2**m):
            x = k(vector(self.to_bits(i, m)))
            l.append((x, self(x)))

        P = PolynomialRing(k, 'x')
        return P.lagrange_polynomial(l)

    def cnf(self, xi=None, yi=None, format=None):
        r"""
        Return a representation of this S-Box in conjunctive normal
        form.

        This function examines the truth tables for each output bit of
        the S-Box and thus has complexity `n * 2^m` for an `m \times n`
        S-Box.

        INPUT:

        - ``xi`` -- (default: ``1...m``) indices for the input variables

        - ``yi`` -- (default: ``m+1 ... m+n``) indices for the output variables

        - ``format`` -- (default: ``None``) output format, see below

        FORMATS:

        - ``None`` -- return a list of tuples of integers where each
          tuple represents a clause, the absolute value of an integer
          represents a variable and the sign of an integer indicates
          inversion

        - ``symbolic`` -- string that can be parsed by the
          ``SymbolicLogic`` package

        - ``dimacs`` -- string in DIMACS format which is the gold
          standard for SAT-solver input (cf. http://www.satlib.org/)

        - ``dimacs_headless`` -- string in DIMACS format, but without
          the header; this is useful for concatenation of outputs

        EXAMPLES:

        We give a very small example to explain the output format::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(1,2,0,3); S
            (1, 2, 0, 3)
            sage: cnf = S.cnf(); cnf
            [(1, 2, -3),  (1, 2, 4),
             (1, -2, 3),  (1, -2, -4),
             (-1, 2, -3), (-1, 2, -4),
             (-1, -2, 3), (-1, -2, 4)]

        This output completely describes the S-Box. For instance, we
        can check that ``S([0,1]) -> [1,0]`` satisfies every clause if
        the first input bit corresponds to the index ``1`` and the
        last output bit corresponds to the index ``3`` in the output.

        We can convert this representation to the DIMACS format::

            sage: print(S.cnf(format='dimacs'))
            p cnf 4 8
            1 2 -3 0
            1 2 4 0
            1 -2 3 0
            1 -2 -4 0
            -1 2 -3 0
            -1 2 -4 0
            -1 -2 3 0
            -1 -2 4 0

        For concatenation we can strip the header::

            sage: print(S.cnf(format='dimacs_headless'))
            1 2 -3 0
            1 2 4 0
            1 -2 3 0
            1 -2 -4 0
            -1 2 -3 0
            -1 2 -4 0
            -1 -2 3 0
            -1 -2 4 0

        This might be helpful in combination with the ``xi`` and
        ``yi`` parameter to assign indices manually::

            sage: print(S.cnf(xi=[10,20],yi=[30,40], format='dimacs_headless'))
            10 20 -30 0
            10 20 40 0
            10 -20 30 0
            10 -20 -40 0
            -10 20 -30 0
            -10 20 -40 0
            -10 -20 30 0
            -10 -20 40 0

        We can also return a string which is parse-able by the
        ``SymbolicLogic`` package::

            sage: log = SymbolicLogic()
            sage: s = log.statement(S.cnf(format='symbolic'))
            sage: log.truthtable(s)[1:]
            [['False', 'False', 'False', 'False', 'False'],
             ['False', 'False', 'False', 'True', 'False'],
             ['False', 'False', 'True', 'False', 'False'],
             ['False', 'False', 'True', 'True', 'True'],
             ['False', 'True', 'False', 'False', 'True'],
             ['False', 'True', 'False', 'True', 'True'],
             ['False', 'True', 'True', 'False', 'True'],
             ['False', 'True', 'True', 'True', 'True'],
             ['True', 'False', 'False', 'False', 'True'],
             ['True', 'False', 'False', 'True', 'True'],
             ['True', 'False', 'True', 'False', 'True'],
             ['True', 'False', 'True', 'True', 'True'],
             ['True', 'True', 'False', 'False', 'True'],
             ['True', 'True', 'False', 'True', 'True'],
             ['True', 'True', 'True', 'False', 'True'],
             ['True', 'True', 'True', 'True', 'True']]

        This function respects endianness of the S-Box::

            sage: S = SBox(1,2,0,3, big_endian=False); S
            (1, 2, 0, 3)
            sage: cnf = S.cnf(); cnf
            [(1, 2, -4), (1, 2, 3),
             (-1, 2, 4), (-1, 2, -3),
             (1, -2, -4), (1, -2, -3),
             (-1, -2, 4), (-1, -2, 3)]

        S-Boxes with m!=n also work:

            sage: o = list(range(8)) + list(range(8))
            sage: shuffle(o)
            sage: S = SBox(o)
            sage: S.is_permutation()
            False

            sage: len(S.cnf()) == 3*2^4
            True

        TESTS:

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(1,2,0,3, big_endian=False)
            sage: S.cnf([1000,1001,1002], [2000,2001,2002])
            Traceback (most recent call last):
            ...
            TypeError: first arg required to have length 2, got 3 instead
        """
        cdef Py_ssize_t m = self.m
        cdef Py_ssize_t n = self.n

        cdef Py_ssize_t i
        if xi is None:
            xi = [i+1 for i in range(m)]

        if yi is None:
            yi = [m+i+1 for i in range(n)]

        if len(xi) != m:
            raise TypeError("first arg required to have length %d, got %d instead" % (m, len(xi)))

        if len(yi) != n:
            raise TypeError("second arg required to have length %d, got %d instead" % (n, len(yi)))

        cdef list output_bits = list(range(n))
        if not self._big_endian:
            output_bits = list(reversed(output_bits))

        C = []  # the set of clauses
        cdef Py_ssize_t e, output_bit, v
        cdef list x, y
        for e in range(1 << m):
            x = self.to_bits(e, m)
            y = <list> self(x)  # evaluate at x
            for output_bit in output_bits:  # consider each bit
                clause = [(-1)**(int(v)) * i for v, i in zip(x, xi)]
                clause.append((-1)**(1-int(y[output_bit])) * yi[output_bit])
                C.append(tuple(clause))

        if format is None:
            return C
        elif format == 'symbolic':
            gd = self.ring().gens()
            formula = []
            for clause in C:
                clause = "|".join([str(gd[abs(v)-1]).replace("-", "~") for v in clause])
                formula.append("("+clause+")")
            return " & ".join(formula)

        elif format.startswith('dimacs'):
            if format == "dimacs_headless":
                header = ""
            else:
                header = "p cnf %d %d\n" % (m + n, len(C))
            values = " 0\n".join([" ".join(map(str, line)) for line in C])
            return header + values + " 0\n"
        else:
            raise ValueError("Format '%s' not supported" % (format,))

    def component_function(self, b):
        r"""
        Return a Boolean function corresponding to the component function
        `b \cdot S(x)`.

        If `S` is an `m \times n` S-Box, then `b \in \GF{2}^n` and
        `\cdot` denotes dot product of two vectors.

        INPUT:

        - ``b`` -- either an integer or a list/tuple/vector of `\GF{2}`
          elements of length ``self.output_size()``

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([7,6,0,4,2,5,1,3])
            sage: f3 = S.component_function(3)
            sage: f3.algebraic_normal_form()                                            # needs sage.rings.polynomial.pbori
            x0*x1 + x0*x2 + x0 + x2

            sage: f5 = S.component_function([1, 0, 1])
            sage: f5.algebraic_normal_form()                                            # needs sage.rings.polynomial.pbori
            x0*x2 + x0 + x1*x2

        TESTS::

            sage: from sage.crypto.sboxes import SBox
            sage: sb = SBox([0, 1, 2, 3, 0, 1, 2, 3])
            sage: sb.component_function([1, 0])
            Boolean function with 3 variables
        """
        cdef Py_ssize_t m = self.m
        cdef Py_ssize_t n = self.n

        try:
            b = list(b)
            if len(b) > n:
                raise ValueError("input (%s) is too long and would be truncated" % (b,))
            b = self.from_bits(b, n)
        except TypeError:
            try:
                b = ZZ(b)
            except TypeError:
                raise TypeError("cannot handle input argument %s" % (b,))

        cdef Py_ssize_t x
        ret = BooleanFunction([ZZ(b & self._S_list[x]).popcount() & 1 for x in range(1 << m)])

        return ret

    def nonlinearity(self):
        """
        Return the nonlinearity of this S-Box.

        The nonlinearity of an S-Box is defined as the minimum nonlinearity
        of all its component functions.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = mq.SR(1,4,4,8).sbox()
            sage: S.nonlinearity()
            112
        """
        return (1 << (self.m-1)) - self.maximal_linear_bias_absolute()

    def linearity(self):
        """
        Return the linearity of this S-Box.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = mq.SR(1, 4, 4, 8).sbox()
            sage: S.linearity()
            32
        """
        return self.maximal_linear_bias_absolute() << 1

    def is_apn(self):
        r"""
        Return ``True`` if this S-Box is an almost perfect nonlinear (APN)
        function.

        An `m \times m` S-Box `S` is called almost perfect nonlinear if for
        every nonzero `\alpha \in \GF{2}^m` and every
        `\beta \in \GF{2}^m`, the equation
        `S(x) \oplus S(x \oplus \alpha) = \beta` has 0 or 2 solutions.
        Equivalently, the differential uniformity of `S` is equal to 2.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([0,1,3,6,7,4,5,2])
            sage: S.is_apn()
            True
            sage: S.differential_uniformity()
            2
        """
        if self.input_size() != self.output_size():
            raise TypeError("APN function is only defined for"
                            "self.input_size() == self.output_size()")
        return self.differential_uniformity() == 2

    def differential_branch_number(self):
        r"""
        Return differential branch number of this S-Box.

        The differential branch number of an S-Box `S` is defined as

        .. MATH::

            \min_{v, w \neq v} \{ \mathrm{wt}(v \oplus w)
            + \mathrm{wt}(S(v) \oplus S(w)) \},

        where `\mathrm{wt}(x)` denotes the Hamming weight of vector `x`.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2])
            sage: S.differential_branch_number()
            3
        """
        cdef Py_ssize_t m = self.m
        cdef Py_ssize_t n = self.n
        cdef Py_ssize_t ret = (1 << m) + (1 << n)

        cdef Py_ssize_t a, b
        cdef Py_ssize_t x, y, w
        for a in range(1 << m):
            for b in range(1 << n):
                if a != b:
                    x = a ^ b
                    y = self._S_list[a] ^ self._S_list[b]
                    w = hamming_weight(x) + hamming_weight(y)
                    if w < ret:
                        ret = w
        return ret

    def linear_branch_number(self):
        r"""
        Return linear branch number of this S-Box.

        The linear branch number of an S-Box `S` is defined as

        .. MATH::

            \min_{\substack{\alpha \neq 0, \beta \\ \mathrm{LAM}(\alpha, \beta) \neq 0}}
                \{ \mathrm{wt}(\alpha) + \mathrm{wt}(\beta) \},

        where `\mathrm{LAM}(\alpha, \beta)` is the entry at row `\alpha` and
        column `\beta` of linear approximation matrix correspond to this
        S-Box. The `\mathrm{wt}(x)` denotes the Hamming weight of `x`.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2])
            sage: S.linear_branch_number()
            2

        TESTS::

            sage: f = SBox([0, 2, 0, 6, 2, 2, 3, 7])
            sage: f.linear_branch_number()
            1
        """
        cdef Py_ssize_t m = self.m
        cdef Py_ssize_t n = self.n
        cdef Matrix lat = <Matrix> self.linear_approximation_table()
        cdef Py_ssize_t ret = (1 << m) + (1 << n)

        cdef Py_ssize_t a, b, w
        for a in range(1 << m):
            for b in range(1, 1 << n):
                if lat.get_unsafe(a, b) != 0:
                    w = hamming_weight(a) + hamming_weight(b)
                    if w < ret:
                        ret = w
        return ret

    @cached_method
    def autocorrelation_table(self):
        r"""
        Return the autocorrelation table corresponding to this S-Box.

        for an `m \times n` S-Box `S`, its autocorrelation table entry at
        row `a \in \GF{2}^m` and column `b \in \GF{2}^n`
        (considering their integer representation) is defined as:

        .. MATH::

            \sum_{x \in \GF{2}^m} (-1)^{b \cdot S(x) \oplus
            b \cdot S(x \oplus a)}.

        Equivalently, the columns `b` of autocorrelation table correspond to
        the autocorrelation spectrum of component function `b \cdot S(x)`.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(7,6,0,4,2,5,1,3)
            sage: S.autocorrelation_table()                                             # needs sage.combinat
            [ 8  8  8  8  8  8  8  8]
            [ 8  0  0  0  0  0  0 -8]
            [ 8  0 -8  0  0  0  0  0]
            [ 8  0  0  0  0 -8  0  0]
            [ 8 -8  0  0  0  0  0  0]
            [ 8  0  0  0  0  0 -8  0]
            [ 8  0  0 -8  0  0  0  0]
            [ 8  0  0  0 -8  0  0  0]
        """
        from sage.combinat.matrices.hadamard_matrix import hadamard_matrix
        A = self.difference_distribution_table() * hadamard_matrix(1 << self.n)
        A.set_immutable()

        return A

    @cached_method
    def boomerang_connectivity_table(self):
        r"""
        Return the boomerang connectivity table (BCT) for this S-Box.

        Boomerang connectivity matrix of an invertible `m \times m`
        S-Box `S` is an `2^m \times 2^m` matrix with entry at row
        `\Delta_i \in \GF{2}^m` and column `\Delta_o \in \GF{2}^m`
        equal to

        .. MATH::

            |\{ x \in \GF{2}^m | S^{-1}( S(x) \oplus \Delta_o) \oplus
               S^{-1}( S(x \oplus \Delta_i) \oplus \Delta_o) = \Delta_i\}|.

        For more results concerning boomerang connectivity matrix,
        see [CHPSS18]_. The algorithm used here is the one from
        Dunkelman [Du2018]_.

        EXAMPLES::

            sage: from sage.crypto.sboxes import PRESENT
            sage: PRESENT.boomerang_connectivity_table()
            [16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16]
            [16  0  4  4  0 16  4  4  4  4  0  0  4  4  0  0]
            [16  0  0  6  0  4  6  0  0  0  2  0  2  2  2  0]
            [16  2  0  6  2  4  4  2  0  0  2  2  0  0  0  0]
            [16  0  0  0  0  4  2  2  0  6  2  0  6  0  2  0]
            [16  2  0  0  2  4  0  0  0  6  2  2  4  2  0  0]
            [16  4  2  0  4  0  2  0  2  0  0  4  2  0  4  8]
            [16  4  2  0  4  0  2  0  2  0  0  4  2  0  4  8]
            [16  4  0  2  4  0  0  2  0  2  0  4  0  2  4  8]
            [16  4  2  0  4  0  2  0  2  0  0  4  2  0  4  8]
            [16  0  2  2  0  4  0  0  6  0  2  0  0  6  2  0]
            [16  2  0  0  2  4  0  0  4  2  2  2  0  6  0  0]
            [16  0  6  0  0  4  0  6  2  2  2  0  0  0  2  0]
            [16  2  4  2  2  4  0  6  0  0  2  2  0  0  0  0]
            [16  0  2  2  0  0  2  2  2  2  0  0  2  2  0  0]
            [16  8  0  0  8  0  0  0  0  0  0  8  0  0  8 16]
        """
        from itertools import product

        cdef SBox Si = self.inverse()

        cdef Py_ssize_t nrows = 1 << self.m
        cdef Py_ssize_t ncols = 1 << self.n

        cdef list L = []

        cdef Py_ssize_t delta_in, x, i, j
        cdef list l, table, row
        for delta_in in range(ncols):
            table = [[] for _ in range(ncols)]
            for x in range(nrows):
                table[x ^ self._S_list[Si._S_list[x] ^ delta_in]].append(x)

            row = [0]*ncols
            for l in table:
                for i, j in product(l, l):
                    row[i ^ j] += 1
            L.extend(row)

        A = matrix(ZZ, nrows, ncols, L)
        A.set_immutable()
        return A

    def boomerang_uniformity(self):
        """
        Return the boomerang uniformity.

        The boomerang uniformity is defined as the highest entry in the
        boomerang connectivity table, ignoring the first row and column.

        EXAMPLES::

            sage: from sage.crypto.sboxes import AES
            sage: AES.boomerang_uniformity()
            6
        """
        bct = self.boomerang_connectivity_table()
        return max(bct.delete_rows([0]).delete_columns([0]).list())

    def linear_structures(self):
        r"""
        Return a list of 3-valued tuple `(b, \alpha, c)` such that `\alpha` is
        a `c`-linear structure of the component function `b \cdot S(x)`.

        A Boolean function `f : \GF{2}^m \mapsto \GF{2}` is said
        to have a `c`-linear structure if there exists a nonzero `\alpha` such
        that `f(x) \oplus f(x \oplus \alpha)` is a constant function `c`.

        An `m \times n` S-Box `S` has a linear structure if there exists a
        component function `b \cdot S(x)` that has a linear structure.

        The three valued tuple `(b, \alpha, c)` shows that `\alpha` is a
        `c`-linear structure of the component function `b \cdot S(x)`. This
        implies that for all output differences `\beta` of the S-Box
        correspond to input difference `\alpha`, we have `b \cdot \beta = c`.

        .. SEEALSO::

            :meth:`is_linear_structure`,
            :meth:`has_linear_structure`

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([0,1,3,6,7,4,5,2])
            sage: S.linear_structures()                                                 # needs sage.combinat
            [(1, 1, 1), (2, 2, 1), (3, 3, 1), (4, 4, 1),
             (5, 5, 1), (6, 6, 1), (7, 7, 1)]
        """
        cdef Py_ssize_t n = self.n
        cdef Py_ssize_t m = self.m
        cdef Matrix act = <Matrix> self.autocorrelation_table()
        cdef list ret = []

        cdef Py_ssize_t j, i, c
        for j in range(1, 1 << n):
            for i in range(1, 1 << m):
                if abs(act.get_unsafe(i, j)) == (1 << m):
                    c = ((1 - (act.get_unsafe(i, j) >> m)) >> 1)
                    ret.append((j, i, c))
        return ret

    def has_linear_structure(self):
        """
        Return ``True`` if there exists a nonzero component function of this
        S-Box that has a linear structure.

        .. SEEALSO::

            :meth:`is_linear_structure`,
            :meth:`linear_structures`.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2)
            sage: S.has_linear_structure()
            True
        """
        cdef Py_ssize_t i
        return any(self.component_function(i).has_linear_structure()
                   for i in range(1, 1 << self.n))

    def is_linear_structure(self, a, b):
        r"""
        Return ``True`` if `a` is a linear structure of the component function
        `b \cdot S(x)` where S is this `m \times n` S-Box.

        INPUT:

        - ``a`` -- either an integer or a tuple of `\GF{2}` elements of
          length equal to the input size of SBox
        - ``b`` -- either an integer or a tuple of `\GF{2}` elements of
          length equal to the output size of SBox

        .. SEEALSO::

            :meth:`linear_structures`,
            :meth:`has_linear_structure`

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2)
            sage: S.component_function(1).autocorrelation()
            (16, -16, 0, 0, 0, 0, 0, 0, -16, 16, 0, 0, 0, 0, 0, 0)
            sage: S.is_linear_structure(1, 1)
            True
            sage: S.is_linear_structure([1, 0, 0, 1], [0, 0, 0, 1])
            True
            sage: S.is_linear_structure([0, 1, 1, 1], 1)
            False
        """
        return self.component_function(b).is_linear_structure(a)

    def max_degree(self):
        """
        Return the maximal algebraic degree of all its component functions.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2])
            sage: S.max_degree()                                                        # needs sage.rings.polynomial.pbori
            3
        """
        ret = ZZ.zero()

        cdef Py_ssize_t i
        for i in range(self.n):
            deg_Si = self.component_function(1 << i).algebraic_degree()
            if deg_Si > ret:
                ret = deg_Si
        return ret

    def min_degree(self):
        """
        Return the minimal algebraic degree of all its component functions.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2])
            sage: S.min_degree()                                                        # needs sage.rings.polynomial.pbori
            2
        """
        ret = ZZ(self.m)

        cdef Py_ssize_t b
        for b in range(1, 1 << self.n):
            deg_bS = self.component_function(b).algebraic_degree()
            if deg_bS < ret:
                ret = deg_bS
        return ret

    def is_balanced(self):
        r"""
        Return ``True`` if this S-Box is balanced.

        An S-Box is balanced if all its component functions are balanced.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2])
            sage: S.is_balanced()
            True
        """
        cdef Py_ssize_t b
        for b in range(1, 1 << self.n):
            bS = self.component_function(b)
            if not bS.is_balanced():
                return False
        return True

    def is_almost_bent(self):
        r"""
        Return ``True`` if this S-Box is an almost bent (AB) function.

        An `m \times m` S-Box `S`, for `m` odd, is called almost bent if its
        nonlinearity is equal to `2^{m-1} - 2^{(m-1)/2}`.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([0,1,3,6,7,4,5,2])
            sage: S.is_almost_bent()
            True
        """
        if self.m != self.n:
            raise TypeError("almost bent function only exists for"
                            " self.input_size() == self.output_size()")

        cdef Py_ssize_t m = self.m

        if is_even(m):
            return False

        return self.nonlinearity() == 2**(m-1) - 2**((m-1)//2)

    def fixed_points(self):
        """
        Return a list of all fixed points of this S-Box.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([0,1,3,6,7,4,5,2])
            sage: S.fixed_points()
            [0, 1]
        """
        cdef Py_ssize_t i
        return [i for i in range(1 << self.m) if i == self._S_list[i]]

    def inverse(self):
        """
        Return the inverse of this S-Box.

        Note that the S-Box must be invertible, otherwise it will raise
        a :exc:`TypeError`.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([0, 1, 3, 6, 7, 4, 5, 2])
            sage: Sinv = S.inverse()
            sage: [Sinv(S(i)) for i in range(8)]
            [0, 1, 2, 3, 4, 5, 6, 7]
        """
        if not self.is_permutation():
            raise TypeError("S-Box must be a permutation")

        cdef Py_ssize_t i
        cdef list L = [self._S_list[i] for i in range(1 << self.m)]

        return SBox([L.index(i) for i in range(1 << self.m)],
                    big_endian=self._big_endian)

    def is_monomial_function(self):
        r"""
        Return ``True`` if this S-Box is a monomial/power function.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([0,1,3,6,7,4,5,2])
            sage: S.is_monomial_function()
            False
            sage: S_poly = S.interpolation_polynomial(); S_poly
            (a + 1)*x^6 + (a^2 + a + 1)*x^5 + (a^2 + a)*x^4
             + (a^2 + 1)*x^3 + a*x^2 + a*x

            sage: all([S(x) == S_poly(x) for x in S_poly.base_ring()])
            True

            sage: S = SBox(0,3,2,1)
            sage: S.interpolation_polynomial()
            x^2
            sage: S.is_monomial_function()
            True
        """
        return self.interpolation_polynomial().is_monomial()

    def is_plateaued(self):
        r"""
        Return ``True`` if this S-Box is plateaued, i.e. for all nonzero
        `b \in \GF{2}^n` the Boolean function `b \cdot S(x)`
        is plateaued.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox(0, 3, 1, 2, 4, 6, 7, 5)
            sage: S.is_plateaued()
            True
        """
        cdef Py_ssize_t b
        for b in range(1, 1 << self.n):
            bS = self.component_function(b)
            if not bS.is_plateaued():
                return False
        return True

    def is_bent(self):
        r"""
        Return ``True`` if this S-Box is bent, i.e. its nonlinearity
        is equal to `2^{m-1} - 2^{m/2 - 1}` where `m` is the input size
        of the S-Box.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: R.<x> = GF(2**2, 'a')[]
            sage: base = R.base_ring()
            sage: a = base.gen()
            sage: G = a * x^2 + 1
            sage: S = SBox([G(x * y**(14)) for x in sorted(base) for y in sorted(base)])
            sage: S.is_bent()
            True
            sage: S.nonlinearity()
            6
            sage: S.linear_approximation_table()
            [ 8 -2  2 -2]
            [ 0 -2  2 -2]
            [ 0 -2  2 -2]
            [ 0 -2  2 -2]
            [ 0 -2  2 -2]
            [ 0 -2 -2  2]
            [ 0  2  2  2]
            [ 0  2 -2 -2]
            [ 0 -2  2 -2]
            [ 0  2 -2 -2]
            [ 0 -2 -2  2]
            [ 0  2  2  2]
            [ 0 -2  2 -2]
            [ 0  2  2  2]
            [ 0  2 -2 -2]
            [ 0 -2 -2  2]
        """
        cdef Py_ssize_t m = self.m
        cdef Py_ssize_t n = self.n

        if not is_even(m) or n > m // 2:
            return False

        return self.nonlinearity() == 2**(m-1) - 2**(m//2 - 1)

    def is_involution(self):
        r"""
        Return ``True`` if this S-Box is an involution, i.e. the inverse S-Box
        is equal itself.

        EXAMPLES::

            sage: from sage.crypto.sbox import SBox
            sage: S = SBox([x**254 for x in sorted(GF(2**8))])
            sage: S.is_involution()
            True
        """
        return self == self.inverse()


cdef Py_ssize_t feistel_substitute(Py_ssize_t x, Py_ssize_t input_size, list sboxes) noexcept:
    """
    Compute a Feistel output using the given sboxes.

    INPUT:

    - ``x`` -- integer; the input to the Feistel construction
    - ``input_size`` -- integer; the bitsize of the Feistel construction
    - ``sboxes`` -- list of SBox; the sboxes applied in the Feistel construction

    EXAMPLES::

        sage: from sage.crypto.sbox import feistel_construction
        sage: from sage.crypto.sboxes import PRESENT as s
        sage: S = feistel_construction(s, s, s) # indirect doctest
    """
    cdef Py_ssize_t mask = (1 << input_size) - 1
    cdef Py_ssize_t xl = (x >> input_size) & mask
    cdef Py_ssize_t xr = x & mask

    cdef SBox sb
    for sb in sboxes:
        xl, xr = sb(xl) ^ xr, xl

    return (xl << input_size) | xr


cdef Py_ssize_t misty_substitute(Py_ssize_t x, Py_ssize_t input_size, list sboxes) noexcept:
    """
    Compute a Misty output using the given sboxes.

    INPUT:

    - ``x`` -- integer; the input to the Misty construction
    - ``input_size`` -- integer; the bitsize of the Misty construction
    - ``sboxes`` -- list of SBox; the sboxes applied in the Misty construction

    EXAMPLES::

        sage: from sage.crypto.sbox import misty_construction
        sage: from sage.crypto.sboxes import PRESENT as s
        sage: S = misty_construction(s, s, s) # indirect doctest
    """
    cdef Py_ssize_t mask = (1 << input_size) - 1
    cdef Py_ssize_t xl = (x >> input_size) & mask
    cdef Py_ssize_t xr = x & mask

    cdef SBox sb
    for sb in sboxes:
        xl, xr = sb(xr) ^ xl, xl

    return (xl << input_size) | xr


ctypedef Py_ssize_t (*_SBOX_CONSTR) (Py_ssize_t, Py_ssize_t, list) noexcept


cdef sbox_construction(_SBOX_CONSTR construction, list args):
    """
    Construct an Sbox from the given input sboxes that has a twice
    as big input size.

    INPUT:

    - ``args`` -- a finite iterable SBox objects

    EXAMPLES::

        sage: from sage.crypto.sbox import feistel_construction
        sage: from sage.crypto.sboxes import PRESENT as s
        sage: S = feistel_construction(s, s, s) # indirect doctest
    """
    if len(args) == 1:
        if isinstance(args[0], SBox):
            sboxes = [args[0]]
        else:
            sboxes = args[0]
    elif len(args) > 1:
        sboxes = args
    else:
        raise TypeError("no input provided")

    for sb in sboxes:
        if not isinstance(sb, SBox):
            raise TypeError("all inputs must be an instance of SBox object")

    cdef Py_ssize_t input_size = sboxes[0].input_size()
    cdef Py_ssize_t m = 2 * input_size

    cdef Py_ssize_t i
    return SBox([construction(i, input_size, sboxes) for i in range(1 << m)])


def feistel_construction(*args):
    r"""
    Return an S-Box constructed by Feistel structure using smaller S-Boxes in
    ``args``. The number of round in the construction is equal to the number of
    S-Boxes provided as input. For more results concerning the differential
    uniformity and the nonlinearity of S-Boxes constructed by Feistel structures
    see [CDL2015]_.

    INPUT:

    - ``args`` -- a finite iterable SBox objects

    EXAMPLES:

    Suppose we construct an `8 \times 8` S-Box with 3-round Feistel
    construction from the S-Box of PRESENT::

        sage: from sage.crypto.sbox import SBox
        sage: s = SBox(12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2)
        sage: from sage.crypto.sbox import feistel_construction
        sage: S = feistel_construction(s, s, s)

    The properties of the constructed S-Box can be easily examined::

        sage: S.nonlinearity()
        96
        sage: S.differential_branch_number()
        2
        sage: S.linear_branch_number()
        2
    """
    return sbox_construction(feistel_substitute, list(args))


def misty_construction(*args):
    r"""
    Return an S-Box constructed by MISTY structure using smaller
    S-Boxes in ``args``.

    The number of round in the construction is equal to the number
    of S-Boxes provided as input. For further result related to the
    nonlinearity and differential uniformity of the constructed S-Box
    one may consult [CDL2015]_.

    INPUT:

    - ``args`` -- a finite iterable SBox objects

    EXAMPLES:

    We construct an `8 \times 8` S-Box using 3-round MISTY structure with the
    following `4 \times 4` S-Boxes `S1, S2, S3` (see Example 2 in [CDL2015]_)::

        sage: from sage.crypto.sbox import SBox
        sage: S1 = SBox([0x4,0x0,0x1,0xF,0x2,0xB,0x6,0x7,0x3,0x9,0xA,0x5,0xC,0xD,0xE,0x8])
        sage: S2 = SBox([0x0,0x0,0x0,0x1,0x0,0xA,0x8,0x3,0x0,0x8,0x2,0xB,0x4,0x6,0xE,0xD])
        sage: S3 = SBox([0x0,0x7,0xB,0xD,0x4,0x1,0xB,0xF,0x1,0x2,0xC,0xE,0xD,0xC,0x5,0x5])
        sage: from sage.crypto.sbox import misty_construction
        sage: S = misty_construction(S1, S2, S3)

    The properties of the constructed S-Box can be easily examined::

        sage: S.differential_uniformity()
        8
        sage: S.linearity()
        64
    """
    return sbox_construction(misty_substitute, list(args))
