# -*- coding: utf-8 -*-
"""
Boolean functions

Those functions are used for example in LFSR based ciphers like
the filter generator or the combination generator.

This module allows to study properties linked to spectral analysis,
and also algebraic immunity.

EXAMPLES::

    sage: R.<x> = GF(2^8,'a')[]
    sage: from sage.crypto.boolean_function import BooleanFunction
    sage: B = BooleanFunction( x^254 ) # the Boolean function Tr(x^254)
    sage: B
    Boolean function with 8 variables
    sage: B.nonlinearity()
    112
    sage: B.algebraic_immunity()
    4

AUTHOR:

- Rusydi H. Makarim (2016-10-13): add functions related to linear structures
- Rusydi H. Makarim (2016-07-09): add is_plateaued()
- Yann Laigle-Chapuy (2010-02-26): add basic arithmetic
- Yann Laigle-Chapuy (2009-08-28): first implementation

"""

from cysignals.signals cimport sig_check
from libc.string cimport memcpy

from sage.structure.sage_object cimport SageObject
from sage.structure.richcmp cimport rich_to_bool
from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.polynomial.pbori.pbori import BooleanPolynomial
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.polynomial.polynomial_element import Polynomial

from sage.misc.superseded import deprecated_function_alias

from sage.data_structures.bitset_base cimport *

# for details about the implementation of hamming_weight (in .pxd),
# walsh_hadamard transform, reed_muller transform, and a lot
# more, see 'Matters computational' available on www.jjj.de.

cdef walsh_hadamard(long *f, int ldn):
    r"""
    The Walsh Hadamard transform is an orthogonal transform equivalent
    to a multidimensional discrete Fourier transform of size 2x2x...x2.
    It can be defined by the following formula:

    .. MATH:: W(j) = \sum_{i\in\{0,1\}^n} (-1)^{f(i)\oplus i \cdot j}

    EXAMPLES::

        sage: from sage.crypto.boolean_function import BooleanFunction
        sage: B = BooleanFunction([1,0,0,1])
        sage: B.walsh_hadamard_transform() # indirect doctest
        (0, 0, 0, -4)
    """
    cdef long n, ldm, m, mh, t1, t2, r, j, u, v
    n = 1 << ldn
    for ldm in range(1, ldn+1):
        m = (1 << ldm)
        mh = m // 2
        # If this is ``for r in range(0, n, m):``, then Cython generates horrible C code
        for 0 <= r < n by m:
            t1 = r
            t2 = r + mh
            for j in range(mh):
                sig_check()
                u = f[t1]
                v = f[t2]
                f[t1] = u + v
                f[t2] = u - v
                t1 += 1
                t2 += 1

cdef long yellow_code(unsigned long a):
    """
    The yellow-code is just a Reed Muller transform applied to a
    word.

    EXAMPLES::

        sage: from sage.crypto.boolean_function import BooleanFunction
        sage: R.<x,y,z> = BooleanPolynomialRing(3)
        sage: P = x*y
        sage: B = BooleanFunction( P )
        sage: B.truth_table() # indirect doctest
        (False, False, False, True, False, False, False, True)
    """
    cdef unsigned long s = (8*sizeof(unsigned long)) >> 1
    cdef unsigned long m = (~0UL) >> s
    cdef unsigned long r = a
    while(s):
        sig_check()
        r ^= (r&m) << s
        s >>= 1
        m ^= (m<<s)
    return r

cdef reed_muller(mp_limb_t* f, int ldn):
    r"""
    The Reed Muller transform (also known as binary Möbius transform)
    is an orthogonal transform. For a function `f` defined by

    .. MATH:: f(x) = \bigoplus_{I\subset\{1,\ldots,n\}} \left(a_I \prod_{i\in I} x_i\right)

    it allows to compute efficiently the ANF from the truth table and
    vice versa, using the formulae:

    .. MATH:: f(x) = \bigoplus_{support(x)\subset I} a_I
    .. MATH:: a_i  = \bigoplus_{I\subset support(x)} f(x)


    EXAMPLES::

        sage: from sage.crypto.boolean_function import BooleanFunction
        sage: R.<x,y,z> = BooleanPolynomialRing(3)
        sage: P = x*y
        sage: B = BooleanFunction( P )
        sage: B.truth_table() # indirect doctest
        (False, False, False, True, False, False, False, True)
    """
    cdef long n, ldm, m, mh, t1, t2, r, j
    n = 1 << ldn
    # intra word transform
    for r in range(n):
        f[r] = yellow_code(f[r])
    # inter word transform
    for ldm in range(1, ldn+1):
        m = 1 << ldm
        mh = m // 2
        # If this is ``for r in range(0, n, m):``, then Cython generates horrible C code
        for 0 <= r < n by m:
            t1 = r
            t2 = r + mh
            for j in range(mh):
                sig_check()
                f[t2] ^= f[t1]
                t1 += 1
                t2 += 1

cdef class BooleanFunction(SageObject):
    r"""
    This module implements Boolean functions represented as a truth table.

    We can construct a Boolean Function from either:

    - an integer - the result is the zero function with ``x`` variables;
    - a list - it is expected to be the truth table of the
      result. Therefore it must be of length a power of 2, and its
      elements are interpreted as Booleans;
    - a string - representing the truth table in hexadecimal;
    - a Boolean polynomial - the result is the corresponding Boolean function;
    - a polynomial P over an extension of GF(2) - the result is
      the Boolean function with truth table ``( Tr(P(x)) for x in
      GF(2^k) )``

    EXAMPLES:

    from the number of variables::

        sage: from sage.crypto.boolean_function import BooleanFunction
        sage: BooleanFunction(5)
        Boolean function with 5 variables

    from a truth table::

        sage: BooleanFunction([1,0,0,1])
        Boolean function with 2 variables

    note that elements can be of different types::

        sage: B = BooleanFunction([False, sqrt(2)])
        sage: B
        Boolean function with 1 variable
        sage: [b for b in B]
        [False, True]

    from a string::

        sage: BooleanFunction("111e")
        Boolean function with 4 variables

    from a :class:`sage.rings.polynomial.pbori.BooleanPolynomial`::

        sage: R.<x,y,z> = BooleanPolynomialRing(3)
        sage: P = x*y
        sage: BooleanFunction( P )
        Boolean function with 3 variables

    from a polynomial over a binary field::

        sage: R.<x> = GF(2^8,'a')[]
        sage: B = BooleanFunction( x^7 )
        sage: B
        Boolean function with 8 variables

    two failure cases::

        sage: BooleanFunction(sqrt(2))
        Traceback (most recent call last):
        ...
        TypeError: unable to init the Boolean function

        sage: BooleanFunction([1, 0, 1])
        Traceback (most recent call last):
        ...
        ValueError: the length of the truth table must be a power of 2
    """

    cdef bitset_t _truth_table
    cdef tuple _walsh_hadamard_transform
    cdef object _nvariables
    cdef object _nonlinearity
    cdef object _correlation_immunity
    cdef object _autocorrelation
    cdef object _absolute_indicator
    cdef object _sum_of_square_indicator

    def __cinit__(self, x):
        r"""
        Construct a Boolean Function.
        The input ``x`` can be either:

        - an integer - the result is the zero function with ``x`` variables;
        - a list - it is expected to be the truth table of the
          result. Therefore it must be of length a power of 2, and its
          elements are interpreted as Booleans;
        - a Boolean polynomial - the result is the corresponding Boolean function;
        - a polynomial P over an extension of GF(2) - the result is
          the Boolean function with truth table ``( Tr(P(x)) for x in
          GF(2^k) )``

        EXAMPLES:

        from the number of variables::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: BooleanFunction(5)
            Boolean function with 5 variables

        from a truth table::

            sage: BooleanFunction([1,0,0,1])
            Boolean function with 2 variables

        note that elements can be of different types::

            sage: B = BooleanFunction([False, sqrt(2)])
            sage: B
            Boolean function with 1 variable
            sage: [b for b in B]
            [False, True]

        from a :class:`sage.rings.polynomial.pbori.BooleanPolynomial`::

            sage: R.<x,y,z> = BooleanPolynomialRing(3)
            sage: P = x*y
            sage: BooleanFunction( P )
            Boolean function with 3 variables

        from a polynomial over a binary field::

            sage: R.<x> = GF(2^8,'a')[]
            sage: B = BooleanFunction( x^7 )
            sage: B
            Boolean function with 8 variables

        two failure cases::

            sage: BooleanFunction(sqrt(2))
            Traceback (most recent call last):
            ...
            TypeError: unable to init the Boolean function

            sage: BooleanFunction([1, 0, 1])
            Traceback (most recent call last):
            ...
            ValueError: the length of the truth table must be a power of 2
        """
        cdef mp_bitcnt_t i
        if isinstance(x, str):
            L = ZZ(len(x))
            if L.is_power_of(2):
                x = ZZ("0x"+x).digits(base=2,padto=4*L)
            else:
                raise ValueError("the length of the truth table must be a power of 2")
        from types import GeneratorType
        if isinstance(x, (list,tuple,GeneratorType)):
            # initialisation from a truth table

            # first, check the length
            L = ZZ(len(x))
            if L.is_power_of(2):
                self._nvariables = L.exact_log(2)
            else:
                raise ValueError("the length of the truth table must be a power of 2")

            # then, initialize our bitset
            bitset_init(self._truth_table, <mp_bitcnt_t> L)
            for i in range(L):
                bitset_set_to(self._truth_table, i, x[i])  # int(x[i])&1)

        elif isinstance(x, BooleanPolynomial):
            # initialisation from a Boolean polynomial
            self._nvariables = ZZ(x.parent().ngens())
            bitset_init(self._truth_table, <mp_bitcnt_t> (1<<self._nvariables))
            bitset_zero(self._truth_table)
            for m in x:
                i = sum([1<<k for k in m.iterindex()])
                bitset_set(self._truth_table, i)
            reed_muller(self._truth_table.bits,
                        ZZ(self._truth_table.limbs).exact_log(2))

        elif isinstance(x, (int, Integer)):
            # initialisation to the zero function
            self._nvariables = ZZ(x)
            bitset_init(self._truth_table, <mp_bitcnt_t> (1<<self._nvariables))
            bitset_zero(self._truth_table)

        elif isinstance(x, Polynomial):
            K = x.base_ring()
            if isinstance(K, FiniteField) and K.characteristic() == 2:
                self._nvariables = K.degree()
                bitset_init(self._truth_table, <mp_bitcnt_t> (1<<self._nvariables))
                bitset_zero(self._truth_table)
                try:
                    from sage.rings.finite_rings.finite_field_givaro import FiniteField_givaro
                except ImportError:
                    FiniteField_givaro = ()
                if isinstance(K, FiniteField_givaro):  # the ordering is not the same in this case
                    for u in K:
                        bitset_set_to(self._truth_table, ZZ(u._vector_().list(),2), (x(u)).trace())
                else:
                    for i,u in enumerate(K):
                        bitset_set_to(self._truth_table, i, (x(u)).trace())
        elif isinstance(x, BooleanFunction):
            self._nvariables = x.nvariables()
            bitset_init(self._truth_table, <mp_bitcnt_t> (1<<self._nvariables))
            bitset_copy(self._truth_table,(<BooleanFunction>x)._truth_table)
        else:
            raise TypeError("unable to init the Boolean function")

    def __dealloc__(self):
        bitset_free(self._truth_table)

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: BooleanFunction(4) #indirect doctest
            Boolean function with 4 variables
        """
        r = "Boolean function with " + self._nvariables.str() + " variable"
        if self._nvariables>1:
            r += "s"
        return r

    def __invert__(self):
        """
        Return the complement Boolean function of `self`.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction([0, 1, 1, 0, 1, 0, 0, 0])
            sage: (~B).truth_table(format='int')
            (1, 0, 0, 1, 0, 1, 1, 1)
        """
        cdef BooleanFunction res=BooleanFunction(self.nvariables())
        bitset_complement(res._truth_table, self._truth_table)
        return res

    def __add__(self, BooleanFunction other):
        """
        Return the element wise sum of ``self`` and ``other``,
        which must have the same number of variables.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: A = BooleanFunction([0, 1, 0, 1, 1, 0, 0, 1])
            sage: B = BooleanFunction([0, 1, 1, 0, 1, 0, 0, 0])
            sage: (A+B).truth_table(format='int')
            (0, 0, 1, 1, 0, 0, 0, 1)

        it also corresponds to the addition of algebraic normal forms::

            sage: S = A.algebraic_normal_form() + B.algebraic_normal_form()
            sage: (A+B).algebraic_normal_form() == S
            True

        TESTS::

            sage: A+BooleanFunction([0,1])
            Traceback (most recent call last):
            ...
            ValueError: the two Boolean functions must have the same number of variables
        """
        if self.nvariables() != other.nvariables():
            raise ValueError("the two Boolean functions must have the same number of variables")
        cdef BooleanFunction res = BooleanFunction(self)
        bitset_xor(res._truth_table, res._truth_table, other._truth_table)
        return res

    def __mul__(self, BooleanFunction other):
        """
        Return the elementwise multiplication of ``self`` and ``other``,
        which must have the same number of variables.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: A = BooleanFunction([0, 1, 0, 1, 1, 0, 0, 1])
            sage: B = BooleanFunction([0, 1, 1, 0, 1, 0, 0, 0])
            sage: (A*B).truth_table(format='int')
            (0, 1, 0, 0, 1, 0, 0, 0)

        it also corresponds to the multiplication of algebraic normal forms::

            sage: P = A.algebraic_normal_form() * B.algebraic_normal_form()
            sage: (A*B).algebraic_normal_form() == P
            True

        TESTS::

            sage: A*BooleanFunction([0,1])
            Traceback (most recent call last):
            ...
            ValueError: the two Boolean functions must have the same number of variables
        """
        if self.nvariables() != other.nvariables():
            raise ValueError("the two Boolean functions must have the same number of variables")
        cdef BooleanFunction res = BooleanFunction(self)
        bitset_and(res._truth_table, res._truth_table, other._truth_table)
        return res

    def __or__(BooleanFunction self, BooleanFunction other):
        """
        Return the concatenation of ``self`` and ``other``,
        which must have the same number of variables.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: A = BooleanFunction([0, 1, 0, 1])
            sage: B = BooleanFunction([0, 1, 1, 0])
            sage: (A|B).truth_table(format='int')
            (0, 1, 0, 1, 0, 1, 1, 0)

            sage: C = A.truth_table() + B.truth_table()
            sage: (A|B).truth_table(format='int') == C
            True

        TESTS::

            sage: A|BooleanFunction([0,1])
            Traceback (most recent call last):
            ...
            ValueError: the two Boolean functions must have the same number of variables
        """
        if (self._nvariables != other.nvariables()):
            raise ValueError("the two Boolean functions must have the same number of variables")

        cdef BooleanFunction res=BooleanFunction(self.nvariables()+1)
        cdef long i

        nb_limbs = self._truth_table.limbs
        if nb_limbs == 1:
            L = len(self)
            for i in range(L):
                res[i  ]=self[i]
                res[i+L]=other[i]
            return res

        memcpy(res._truth_table.bits,
               self._truth_table.bits, nb_limbs * sizeof(unsigned long))
        memcpy(&(res._truth_table.bits[nb_limbs]),
               other._truth_table.bits, nb_limbs * sizeof(unsigned long))

        return res

    def algebraic_normal_form(self):
        """
        Return the :class:`sage.rings.polynomial.pbori.BooleanPolynomial`
        corresponding to the algebraic normal form.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction([0,1,1,0,1,0,1,1])
            sage: P = B.algebraic_normal_form()
            sage: P
            x0*x1*x2 + x0 + x1*x2 + x1 + x2
            sage: [ P(*ZZ(i).digits(base=2,padto=3)) for i in range(8) ]
            [0, 1, 1, 0, 1, 0, 1, 1]
        """
        cdef bitset_t anf
        cdef mp_bitcnt_t inf, sup
        bitset_init(anf, <mp_bitcnt_t> (1<<self._nvariables))
        bitset_copy(anf, self._truth_table)
        reed_muller(anf.bits, ZZ(anf.limbs).exact_log(2))
        from sage.rings.polynomial.pbori.pbori import BooleanPolynomialRing
        R = BooleanPolynomialRing(self._nvariables,"x")
        G = R.gens()
        P = R(0)

        cdef long i, j, k
        for i in range(anf.limbs):
            if anf.bits[i]:
                inf = i*sizeof(long)*8
                sup = min((i+1)*sizeof(long)*8, (1<<self._nvariables))
                for j in range(inf, sup):
                    if bitset_in(anf,j):
                        m = R(1)
                        for k in range(self._nvariables):
                            if (j>>k)&1:
                                m *= G[k]
                        P += m
        bitset_free(anf)
        return P

    def nvariables(self):
        """
        The number of variables of this function.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: BooleanFunction(4).nvariables()
            4
        """
        return self._nvariables

    def truth_table(self, format='bin'):
        """
        The truth table of the Boolean function.

        INPUT: a string representing the desired format, can be either

        - 'bin' (default) : we return a tuple of Boolean values
        - 'int' : we return a tuple of 0 or 1 values
        - 'hex' : we return a string representing the truth_table in hexadecimal

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: R.<x,y,z> = BooleanPolynomialRing(3)
            sage: B = BooleanFunction( x*y*z + z + y + 1 )
            sage: B.truth_table()
            (True, True, False, False, False, False, True, False)
            sage: B.truth_table(format='int')
            (1, 1, 0, 0, 0, 0, 1, 0)
            sage: B.truth_table(format='hex')
            '43'

            sage: BooleanFunction('00ab').truth_table(format='hex')
            '00ab'

            sage: H = '0abbacadabbacad0'
            sage: len(H)
            16
            sage: T = BooleanFunction(H).truth_table(format='hex')
            sage: T == H
            True
            sage: H = H * 4
            sage: T = BooleanFunction(H).truth_table(format='hex')
            sage: T == H
            True
            sage: H = H * 4
            sage: T = BooleanFunction(H).truth_table(format='hex')
            sage: T == H
            True
            sage: len(T)
            256
            sage: B.truth_table(format='oct')
            Traceback (most recent call last):
            ...
            ValueError: unknown output format
        """
        if format == 'bin':
            return tuple(self)
        if format == 'int':
            return tuple(map(int,self))
        if format == 'hex':
            S = ZZ(self.truth_table(),2).str(16)
            S = "0"*((1<<(self._nvariables-2)) - len(S)) + S
            return S
        raise ValueError("unknown output format")

    def __len__(self):
        """
        Return the number of different input values.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: len(BooleanFunction(4))
            16
        """
        return 2**self._nvariables

    def __richcmp__(BooleanFunction self, other, int op):
        """
        Boolean functions are considered to be equal if the number of
        input variables is the same, and all the values are equal.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: b1 = BooleanFunction([0,1,1,0])
            sage: b2 = BooleanFunction([0,1,1,0])
            sage: b3 = BooleanFunction([0,1,1,1])
            sage: b4 = BooleanFunction([0,1])
            sage: b1 == b2
            True
            sage: b1 == b3
            False
            sage: b1 == b4
            False
        """
        if not isinstance(other, BooleanFunction):
            return NotImplemented
        o = <BooleanFunction> other
        return rich_to_bool(op, bitset_cmp(self._truth_table, o._truth_table))

    def __call__(self, x):
        """
        Return the value of the function for the given input.

        INPUT: either

        - a list -- then all elements are evaluated as Booleans
        - an integer -- then we consider its binary representation

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction([0,1,0,0])
            sage: B(1)
            1
            sage: B([1,0])
            1
            sage: B(4)
            Traceback (most recent call last):
            ...
            IndexError: index out of bound
        """
        if isinstance(x, (int, Integer)):
            if x >= self._truth_table.size:
                raise IndexError("index out of bound")
            return bitset_in(self._truth_table, <mp_bitcnt_t> x)
        elif isinstance(x, list):
            if len(x) != self._nvariables:
                raise ValueError("bad number of inputs")
            return self(ZZ([bool(_) for _ in x], 2))
        else:
            raise TypeError("cannot apply Boolean function to provided element")

    def __iter__(self):
        """
        Iterate through the value of the function.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction([0,1,1,0,1,0,1,0])
            sage: [int(b) for b in B]
            [0, 1, 1, 0, 1, 0, 1, 0]
        """
        return BooleanFunctionIterator(self)

    def _walsh_hadamard_transform_cached(self):
        """
        Return the cached Walsh Hadamard transform. *Unsafe*, no check.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction(3)
            sage: W = B.walsh_hadamard_transform()
            sage: B._walsh_hadamard_transform_cached() is W
            True
        """
        return self._walsh_hadamard_transform

    cpdef tuple walsh_hadamard_transform(self):
        r"""
        Compute the Walsh Hadamard transform `W` of the function `f`.

        .. MATH:: W(j) = \sum_{i\in\{0,1\}^n} (-1)^{f(i)\oplus i \cdot j}

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: R.<x> = GF(2^3,'a')[]
            sage: B = BooleanFunction( x^3 )
            sage: B.walsh_hadamard_transform()
            (0, -4, 0, 4, 0, 4, 0, 4)
        """
        cdef long *temp
        cdef mp_bitcnt_t i,n

        if self._walsh_hadamard_transform is None:
            n = self._truth_table.size
            temp = <long *>sig_malloc(sizeof(long)*n)

            for i in range(n):
                temp[i] = 1 - (bitset_in(self._truth_table, i) << 1)

            walsh_hadamard(temp, self._nvariables)
            self._walsh_hadamard_transform = tuple([temp[i] for i in range(n)])
            sig_free(temp)

        return self._walsh_hadamard_transform

    def absolute_walsh_spectrum(self):
        """
        Return the absolute Walsh spectrum fo the function.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("7969817CC5893BA6AC326E47619F5AD0")
            sage: sorted(B.absolute_walsh_spectrum().items())
            [(0, 64), (16, 64)]

            sage: B = BooleanFunction("0113077C165E76A8")
            sage: B.absolute_walsh_spectrum()
            {8: 64}
        """
        d = {}
        cdef long i
        for i in self.walsh_hadamard_transform():
            if abs(i) in d:
                d[abs(i)] += 1
            else:
                d[abs(i)] = 1
        return d

    def is_balanced(self):
        """
        Return ``True`` if the function takes the value ``True`` half of the time.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction(1)
            sage: B.is_balanced()
            False
            sage: B[0] = True
            sage: B.is_balanced()
            True
        """
        return self.walsh_hadamard_transform()[0] == 0

    def is_symmetric(self):
        """
        Return ``True`` if the function is symmetric, i.e. invariant under
        permutation of its input bits.

        Another way to see it is that the
        output depends only on the Hamming weight of the input.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction(5)
            sage: B[3] = 1
            sage: B.is_symmetric()
            False
            sage: V_B = [0, 1, 1, 0, 1, 0]
            sage: for i in srange(32): B[i] = V_B[i.popcount()]
            sage: B.is_symmetric()
            True
        """
        cdef mp_bitcnt_t i
        cdef list T = [self(2**i-1) for i in range(self._nvariables+1)]
        for i in range(1 << self._nvariables):
            sig_check()
            if T[hamming_weight(i)] != bitset_in(self._truth_table, i):
                return False
        return True

    def nonlinearity(self):
        """
        Return the nonlinearity of the function.

        This is the distance to the linear functions, or the number of
        output ones need to change to obtain a linear function.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction(5)
            sage: B[1] = B[3] = 1
            sage: B.nonlinearity()
            2
            sage: B = BooleanFunction("0113077C165E76A8")
            sage: B.nonlinearity()
            28
        """
        cdef long w
        if self._nonlinearity is None:
            self._nonlinearity = \
                ((1<<self._nvariables) - max(abs(w) for w in self.walsh_hadamard_transform())) >> 1
        return self._nonlinearity

    def is_bent(self):
        """
        Return ``True`` if the function is bent.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("0113077C165E76A8")
            sage: B.is_bent()
            True
        """
        if (self._nvariables & 1):
            return False
        return self.nonlinearity() == ((1<<self._nvariables)-(1<<(self._nvariables//2)))>>1

    def correlation_immunity(self):
        """
        Return the maximum value `m` such that the function is
        correlation immune of order `m`.

        A Boolean function is said to be correlation immune of order
        `m`, if the output of the function is statistically
        independent of the combination of any m of its inputs.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("7969817CC5893BA6AC326E47619F5AD0")
            sage: B.correlation_immunity()
            2
        """
        cdef long c, i
        if self._correlation_immunity is None:
            c = self._nvariables
            W = self.walsh_hadamard_transform()
            for i in range(len(W)):
                sig_check()
                if W[i]:
                    c = min(c, hamming_weight(i))
            self._correlation_immunity = ZZ(c-1)
        return self._correlation_immunity

    def resiliency_order(self):
        """
        Return the maximum value `m` such that the function is
        resilient of order `m`.

        A Boolean function is said to be resilient of order `m` if it
        is balanced and correlation immune of order `m`.

        If the function is not balanced, we return -1.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("077CE5A2F8831A5DF8831A5D077CE5A26996699669699696669999665AA5A55A")
            sage: B.resiliency_order()
            3
        """
        if not self.is_balanced():
            return -1
        return self.correlation_immunity()

    def autocorrelation(self):
        r"""
        Return the autocorrelation of the function, defined by

        .. MATH:: \Delta_f(j) = \sum_{i\in\{0,1\}^n} (-1)^{f(i)\oplus f(i\oplus j)}.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("03")
            sage: B.autocorrelation()
            (8, 8, 0, 0, 0, 0, 0, 0)
        """
        cdef long *temp
        cdef long i

        if self._autocorrelation is None:
            n = self._truth_table.size
            temp = <long *>sig_malloc(sizeof(long)*n)
            W = self.walsh_hadamard_transform()

            for i in range(n):
                sig_check()
                temp[i] = W[i]*W[i]

            walsh_hadamard(temp, self._nvariables)
            self._autocorrelation = tuple([temp[i] >> self._nvariables
                                           for i in range(n)])
            sig_free(temp)

        return self._autocorrelation

    def absolute_autocorrelation(self):
        """
        Return the absolute autocorrelation of the function.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("7969817CC5893BA6AC326E47619F5AD0")
            sage: sorted(B.absolute_autocorrelation().items())
            [(0, 33), (8, 58), (16, 28), (24, 6), (32, 2), (128, 1)]
        """
        d = {}
        cdef long i
        for i in self.autocorrelation():
            if abs(i) in d:
                d[abs(i)] += 1
            else:
                d[abs(i)] = 1
        return d

    def absolute_indicator(self):
        """
        Return the absolute indicator of the function.

        The absolute indicator is defined as the maximal absolute value of
        the autocorrelation.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("7969817CC5893BA6AC326E47619F5AD0")
            sage: B.absolute_indicator()
            32

        The old method's name contained a typo, it is deprecated::

            sage: B.absolut_indicator()
            doctest:warning
            ...
            DeprecationWarning: absolut_indicator is deprecated. Please use absolute_indicator instead.
            See https://github.com/sagemath/sage/issues/28001 for details.
            32
        """
        cdef long a
        if self._absolute_indicator is None:
            D = self.autocorrelation()
            self._absolute_indicator = max([abs(a) for a in D[1:]])
        return self._absolute_indicator

    absolut_indicator = deprecated_function_alias(28001, absolute_indicator)

    def sum_of_square_indicator(self):
        """
        Return the sum of square indicator of the function.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("7969817CC5893BA6AC326E47619F5AD0")
            sage: B.sum_of_square_indicator()
            32768
        """
        cdef long a
        if self._sum_of_square_indicator is None:
            D = self.autocorrelation()
            self._sum_of_square_indicator = sum(a**2 for a in D)
        return self._sum_of_square_indicator

    def annihilator(self,d, dim = False):
        r"""
        Return (if it exists) an annihilator of the boolean function of
        degree at most `d`, that is a Boolean polynomial `g` such that

        .. MATH::

            f(x)g(x) = 0 \forall x.

        INPUT:

        - ``d`` -- an integer;
        - ``dim`` -- a Boolean (default: False), if True, return also
          the dimension of the annihilator vector space.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: f = BooleanFunction("7969817CC5893BA6AC326E47619F5AD0")
            sage: f.annihilator(1) is None
            True
            sage: g = BooleanFunction( f.annihilator(3) )
            sage: set([ fi*g(i) for i,fi in enumerate(f) ])
            {0}
        """
        # NOTE: this is a toy implementation
        from sage.rings.polynomial.polynomial_ring_constructor import BooleanPolynomialRing_constructor
        R = BooleanPolynomialRing_constructor(self._nvariables,'x')
        G = R.gens()
        r = [R(1)]

        from sage.modules.free_module_element import vector
        s = vector(self.truth_table()).support()

        from sage.combinat.combination import Combinations
        from sage.misc.misc_c import prod

        from sage.matrix.constructor import Matrix
        from sage.arith.misc import binomial
        M = Matrix(GF(2), sum(binomial(self._nvariables,i) for i in range(d+1)), len(s))

        cdef long i
        for i in range(1, d+1):
            C = Combinations(self._nvariables, i)
            for c in C:
                sig_check()
                r.append(prod([G[i] for i in c]))

        cdef BooleanFunction t

        cdef long j
        cdef mp_bitcnt_t v

        for i,m in enumerate(r):
            t = BooleanFunction(m)
            for j,v in enumerate(s):
                sig_check()
                M[i,j] = bitset_in(t._truth_table,v)

        kg = M.kernel().gens()

        if kg:
            res = sum([kg[0][i]*ri for i,ri in enumerate(r)])
        else:
            res = None

        if dim:
            return res, len(kg)
        else:
            return res

    def algebraic_immunity(self, annihilator = False):
        """
        Return the algebraic immunity of the Boolean function.

        This is the smallest integer `i` such that there exists a non
        trivial annihilator for `self` or `~self`.

        INPUT:

        - ``annihilator`` -- boolean (default: ``False``); if ``True``,
          returns also an annihilator of minimal degree

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: R.<x0,x1,x2,x3,x4,x5> = BooleanPolynomialRing(6)
            sage: B = BooleanFunction(x0*x1 + x1*x2 + x2*x3 + x3*x4 + x4*x5)
            sage: B.algebraic_immunity(annihilator=True)
            (2, x0*x1 + x1*x2 + x2*x3 + x3*x4 + x4*x5 + 1)
            sage: B[0] +=1
            sage: B.algebraic_immunity()
            2

            sage: R.<x> = GF(2^8,'a')[]
            sage: B = BooleanFunction(x^31)
            sage: B.algebraic_immunity()
            4
        """
        f = self
        g = ~self
        cdef long i
        for i in range(self._nvariables):
            for fun in [f, g]:
                A = fun.annihilator(i)
                if A is not None:
                    if annihilator:
                        return i, A
                    else:
                        return i
        assert False, "you just found a bug!"

    def algebraic_degree(self):
        r"""
        Return the algebraic degree of this Boolean function.

        The algebraic degree of a Boolean function is defined as the degree
        of its algebraic normal form. Note that the degree of the constant
        zero function is defined to be equal to -1.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B.<x0, x1, x2, x3> = BooleanPolynomialRing()
            sage: f = BooleanFunction(x1*x2 + x1*x2*x3 + x1)
            sage: f.algebraic_degree()
            3
            sage: g = BooleanFunction([0, 0])
            sage: g.algebraic_degree()
            -1
        """
        return self.algebraic_normal_form().degree()

    def is_plateaued(self):
        r"""
        Return ``True`` if this function is plateaued, i.e. its Walsh transform
        takes at most three values `0` and `\pm \lambda`, where `\lambda` is some
        positive integer.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: R.<x0, x1, x2, x3> = BooleanPolynomialRing()
            sage: f = BooleanFunction(x0*x1 + x2 + x3)
            sage: f.walsh_hadamard_transform()
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, -8)
            sage: f.is_plateaued()
            True
        """
        W = self.absolute_walsh_spectrum()
        return (len(W) == 1) or (len(W) == 2 and 0 in W)

    def is_linear_structure(self, val):
        r"""
        Return ``True`` if ``val`` is a linear structure of this Boolean
        function.

        INPUT:

        - ``val`` -- either an integer or a tuple/list of `\GF{2}` elements
          of length equal to the number of variables

        .. SEEALSO::

            :meth:`has_linear_structure`,
            :meth:`linear_structures`.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: f = BooleanFunction([0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0])
            sage: f.is_linear_structure(1)
            True
            sage: l = [1, 0, 0, 1]
            sage: f.is_linear_structure(l)
            True
            sage: v = vector(GF(2), l)
            sage: f.is_linear_structure(v)
            True
            sage: f.is_linear_structure(7)
            False
            sage: f.is_linear_structure(20) #parameter is out of range
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: v = vector(GF(3), [1, 0, 1, 1])
            sage: f.is_linear_structure(v)
            Traceback (most recent call last):
            ...
            TypeError: base ring of input vector must be GF(2)
            sage: v = vector(GF(2), [1, 0, 1, 1, 1])
            sage: f.is_linear_structure(v)
            Traceback (most recent call last):
            ...
            TypeError: input vector must be an element of a vector space with dimension 4
            sage: f.is_linear_structure('X') #failure case
            Traceback (most recent call last):
            ...
            TypeError: cannot compute is_linear_structure() using parameter X
        """
        from sage.structure.element import is_Vector
        nvars = self._nvariables

        if isinstance(val, (tuple, list)):
            i = ZZ(val, base=2)
        elif is_Vector(val):
            if val.base_ring() != GF(2):
                raise TypeError("base ring of input vector must be GF(2)")
            elif val.parent().dimension() != nvars:
                raise TypeError("input vector must be an element of a vector"
                                " space with dimension %d" % (nvars,))
            i = ZZ(val.list(), base=2)
        else:
            i = val

        a = self.autocorrelation()
        try:
            return abs(a[i]) == 1 << nvars
        except IndexError:
            raise IndexError("index out of range")
        except TypeError:
            raise TypeError("cannot compute is_linear_structure() using parameter %s" % (val,))

    def has_linear_structure(self):
        r"""
        Return ``True`` if this function has a linear structure.

        An `n`-variable Boolean function `f` has a linear structure if
        there exists a nonzero `a \in \GF{2}^n` such that
        `f(x \oplus a) \oplus f(x)` is a constant function.

        .. SEEALSO::

            :meth:`is_linear_structure`,
            :meth:`linear_structures`.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: f = BooleanFunction([0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0])
            sage: f.has_linear_structure()
            True
            sage: f.autocorrelation()
            (16, -16, 0, 0, 0, 0, 0, 0, -16, 16, 0, 0, 0, 0, 0, 0)
            sage: g = BooleanFunction([0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1])
            sage: g.has_linear_structure()
            False
            sage: g.autocorrelation()
            (16, 4, 4, 4, 4, -4, -4, -4, -4, 4, -4, -4, -4, 4, -4, -4)
        """
        a = self.autocorrelation()
        nvars = self._nvariables
        cdef long i
        return any(abs(a[i]) == 1 << nvars for i in range(1, 1 << nvars))

    def linear_structures(self):
        r"""
        Return all linear structures of this Boolean function as a
        vector subspace of `\GF{2}^n`.

        .. SEEALSO::

            :meth:`is_linear_structure`,
            :meth:`has_linear_structure`.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: f = BooleanFunction([0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0])
            sage: LS = f.linear_structures()
            sage: LS.dimension()
            2
            sage: LS.basis_matrix()
            [1 0 0 0]
            [0 0 0 1]
            sage: LS.list()
            [(0, 0, 0, 0), (1, 0, 0, 0), (0, 0, 0, 1), (1, 0, 0, 1)]
        """
        from sage.modules.free_module import VectorSpace

        cdef long i
        nvars = self.nvariables()
        a = self.autocorrelation()
        l = [ZZ(i).digits(base=2, padto=nvars) for i in range(1<<nvars) if abs(a[i]) == 1<<nvars]
        V = VectorSpace(GF(2), nvars)
        return V.subspace(l)

    def derivative(self, u):
        r"""
        Return the derivative in direction of ``u``

        INPUT:

        - ``u`` -- either an integer or a tuple/list of `\GF{2}` elements
          of length equal to the number of variables


        The derivative of `f` in direction of `u` is defined as
        `x \mapsto f(x) + f(x + u)`.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: f = BooleanFunction([0,1,0,1,0,1,0,1])
            sage: f.derivative(1).algebraic_normal_form()
            1
            sage: u = [1,0,0]
            sage: f.derivative(u).algebraic_normal_form()
            1
            sage: v = vector(GF(2), u)
            sage: f.derivative(u).algebraic_normal_form()
            1
            sage: f.derivative(8).algebraic_normal_form()
            Traceback (most recent call last):
            ...
            IndexError: index out of bound
        """
        from sage.structure.element import is_Vector
        nvars = self._nvariables

        if isinstance(u, (tuple, list)):
            v = ZZ(u, base=2)
        elif is_Vector(u):
            if u.base_ring() != GF(2):
                raise TypeError("base ring of input vector must be GF(2)")
            elif u.parent().dimension() != nvars:
                raise TypeError("input vector must be an element of a vector space with dimension %d" % (nvars,))
            v = ZZ(u.list(), base=2)
        else:
            v = u

        return BooleanFunction([self(x) ^ self(x ^ v)
                                for x in range(1 << nvars)])

    def __setitem__(self, i, y):
        """
        Set a value of the function.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction([0,0,1,1])
            sage: B[0] = 1
            sage: B[2] = (3**17 == 9)
            sage: [b for b in B]
            [True, False, False, True]

        We take care to clear cached values::

            sage: W = B.walsh_hadamard_transform()
            sage: B[2] = 1
            sage: B._walsh_hadamard_transform_cached() is None
            True
        """
        self._clear_cache()
        bitset_set_to(self._truth_table, int(i), int(y)&1)

    def __getitem__(self, i):
        """
        Return the value of the function for the given input.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction([0,1,1,1])
            sage: [ int(B[i]) for i in range(len(B)) ]
            [0, 1, 1, 1]
        """
        return self(i)

    def _clear_cache(self):
        """
        Clear cached values.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction([0,1,1,0])
            sage: W = B.walsh_hadamard_transform()
            sage: B._walsh_hadamard_transform_cached() is None
            False
            sage: B._clear_cache()
            sage: B._walsh_hadamard_transform_cached() is None
            True
        """
        self._walsh_hadamard_transform = None
        self._nonlinearity = None
        self._correlation_immunity = None
        self._autocorrelation = None
        self._absolute_indicator = None
        self._sum_of_square_indicator = None

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction([0,1,1,0])
            sage: loads(dumps(B)) == B
            True
        """
        return unpickle_BooleanFunction, (self.truth_table(format='hex'),)


def unpickle_BooleanFunction(bool_list):
    """
    Specific function to unpickle Boolean functions.

    EXAMPLES::

        sage: from sage.crypto.boolean_function import BooleanFunction
        sage: B = BooleanFunction([0,1,1,0])
        sage: loads(dumps(B)) == B # indirect doctest
        True
    """
    return BooleanFunction(bool_list)


cdef class BooleanFunctionIterator:
    cdef long index, last
    cdef BooleanFunction f

    def __init__(self, f):
        """
        Iterator through the values of a Boolean function.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction(3)
            sage: type(B.__iter__())
            <class 'sage.crypto.boolean_function.BooleanFunctionIterator'>
        """
        self.f = f
        self.index = -1
        self.last = self.f._truth_table.size-1

    def __iter__(self):
        """
        Iterator through the values of a Boolean function.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction(1)
            sage: [b for b in B] # indirect doctest
            [False, False]
        """
        return self

    def __next__(self):
        """
        Next value.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction(1)
            sage: I = B.__iter__()
            sage: next(I)
            False
        """
        if self.index == self.last:
            raise StopIteration
        self.index += 1
        return bitset_in(self.f._truth_table, self.index)

##########################################
# Below we provide some constructions of #
# cryptographic Boolean function.        #
##########################################


def random_boolean_function(n):
    """
    Return a random Boolean function with `n` variables.

    EXAMPLES::

        sage: from sage.crypto.boolean_function import random_boolean_function
        sage: B = random_boolean_function(9)
        sage: B.nvariables()
        9
        sage: while not (210 < B.nonlinearity() < 220):
        ....:     B = random_boolean_function(9)
    """
    from sage.misc.randstate import current_randstate
    r = current_randstate().python_random()
    cdef BooleanFunction B = BooleanFunction(n)
    cdef bitset_t T
    cdef long i
    T[0] = B._truth_table[0]
    for i in range(T.limbs):
        sig_check()
        T.bits[i] = r.randrange(0,Integer(1)<<(sizeof(unsigned long)*8))
    return B
