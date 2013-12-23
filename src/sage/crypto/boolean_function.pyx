"""
Boolean functions

Those functions are used for example in LFSR based ciphers like
the filter generator or the combination generator.

This module allows to study properties linked to spectral analysis,
and also algebraic immunity.

EXAMPLES::

    sage: R.<x>=GF(2^8,'a')[]
    sage: from sage.crypto.boolean_function import BooleanFunction
    sage: B = BooleanFunction( x^254 ) # the Boolean function Tr(x^254)
    sage: B
    Boolean function with 8 variables
    sage: B.nonlinearity()
    112
    sage: B.algebraic_immunity()
    4

AUTHOR:

- Yann Laigle-Chapuy (2010-02-26): add basic arithmetic
- Yann Laigle-Chapuy (2009-08-28): first implementation

"""

from sage.structure.sage_object cimport SageObject

from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer
from sage.rings.finite_rings.constructor import GF
from sage.rings.polynomial.pbori import BooleanPolynomial
from sage.rings.finite_rings.constructor import is_FiniteField
from sage.rings.finite_rings.finite_field_givaro import FiniteField_givaro
from sage.rings.polynomial.polynomial_element import is_Polynomial

include "sage/misc/bitset.pxi"
from cpython.string cimport *

# for details about the implementation of hamming_weight_int,
# walsh_hadamard transform, reed_muller transform, and a lot
# more, see 'Matters computational' available on www.jjj.de.

cdef inline unsigned int hamming_weight_int(unsigned int x):
    # valid for 32bits
    x -=  (x>>1) & 0x55555555UL                        # 0-2 in 2 bits
    x  = ((x>>2) & 0x33333333UL) + (x & 0x33333333UL)  # 0-4 in 4 bits
    x  = ((x>>4) + x) & 0x0f0f0f0fUL                   # 0-8 in 8 bits
    x *= 0x01010101UL
    return x>>24

cdef walsh_hadamard(long *f, int ldn):
    r"""
    The Walsh Hadamard transform is an orthogonal transform equivalent
    to a multidimensional discrete Fourier transform of size 2x2x...x2.
    It can be defined by the following formula:

    .. math:: W(j) = \sum_{i\in\{0,1\}^n} (-1)^{f(i)\oplus i \cdot j}

    EXAMPLES::

        sage: from sage.crypto.boolean_function import BooleanFunction
        sage: B = BooleanFunction([1,0,0,1])
        sage: B.walsh_hadamard_transform() # indirect doctest
        (0, 0, 0, 4)
    """
    cdef long n, ldm, m, mh, t1, t2, r
    n = 1 << ldn
    for 1 <= ldm <= ldn:
        m  = (1<<ldm)
        mh = m//2
        for 0 <= r <n by m:
            t1 = r
            t2 = r+mh
            for 0 <= j < mh:
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
    cdef unsigned long s = (8*sizeof(unsigned long))>>1
    cdef unsigned long m = (~0UL) >> s
    cdef unsigned long r = a
    while(s):
        r ^= ( (r&m) << s )
        s >>= 1
        m ^= (m<<s)
    return r

cdef reed_muller(unsigned long *f, int ldn):
    r"""
    The Reed Muller transform (also known as binary Moebius transform)
    is an orthogonal transform. For a function `f` defined by

    .. math:: f(x) = \bigoplus_{I\subset\{1,\ldots,n\}} \left(a_I \prod_{i\in I} x_i\right)

    it allows to compute efficiently the ANF from the truth table and
    vice versa, using the formulae:

    .. math:: f(x) = \bigoplus_{support(x)\subset I} a_I
    .. math:: a_i  = \bigoplus_{I\subset support(x)} f(x)


    EXAMPLES::

        sage: from sage.crypto.boolean_function import BooleanFunction
        sage: R.<x,y,z> = BooleanPolynomialRing(3)
        sage: P = x*y
        sage: B = BooleanFunction( P )
        sage: B.truth_table() # indirect doctest
        (False, False, False, True, False, False, False, True)
    """
    cdef long n, ldm, m, mh, t1, t2, r
    n = 1 << ldn
    # intra word transform
    for 0 <= r < n:
        f[r] = yellow_code(f[r])
    # inter word transform
    for 1 <= ldm <= ldn:
        m  = (1<<ldm)
        mh = m//2
        for 0 <= r <n by m:
            t1 = r
            t2 = r+mh
            for 0 <= j < mh:
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
    cdef object _walsh_hadamard_transform
    cdef object _nvariables
    cdef object _nonlinearity
    cdef object _correlation_immunity
    cdef object _autocorrelation
    cdef object _absolut_indicator
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
        if PyString_Check(x):
            L = ZZ(len(x))
            if L.is_power_of(2):
                x = ZZ("0x"+x).digits(base=2,padto=4*L)
            else:
                raise ValueError, "the length of the truth table must be a power of 2"
        from types import GeneratorType
        if isinstance(x, (list,tuple,GeneratorType)):
        # initialisation from a truth table

            # first, check the length
            L = ZZ(len(x))
            if L.is_power_of(2):
                self._nvariables = L.exact_log(2)
            else:
                raise ValueError, "the length of the truth table must be a power of 2"

            # then, initialize our bitset
            bitset_init(self._truth_table, L)
            for 0<= i < L:
                bitset_set_to(self._truth_table, i, x[i])#int(x[i])&1)

        elif isinstance(x, BooleanPolynomial):
        # initialisation from a Boolean polynomial
            self._nvariables = ZZ(x.parent().ngens())
            bitset_init(self._truth_table, (1<<self._nvariables))
            bitset_zero(self._truth_table)
            for m in x:
                i = sum( [1<<k for k in m.iterindex()] )
                bitset_set(self._truth_table, i)
            reed_muller(self._truth_table.bits, ZZ(self._truth_table.limbs).exact_log(2) )

        elif isinstance(x, (int,long,Integer) ):
        # initialisation to the zero function
            self._nvariables = ZZ(x)
            bitset_init(self._truth_table,(1<<self._nvariables))
            bitset_zero(self._truth_table)

        elif is_Polynomial(x):
            K = x.base_ring()
            if is_FiniteField(K) and K.characteristic() == 2:
                self._nvariables = K.degree()
                bitset_init(self._truth_table,(1<<self._nvariables))
                bitset_zero(self._truth_table)
                if isinstance(K,FiniteField_givaro): #the ordering is not the same in this case
                    for u in K:
                        bitset_set_to(self._truth_table, ZZ(u._vector_().list(),2) , (x(u)).trace())
                else:
                    for i,u in enumerate(K):
                        bitset_set_to(self._truth_table, i , (x(u)).trace())
        elif isinstance(x, BooleanFunction):
            self._nvariables = x.nvariables()
            bitset_init(self._truth_table,(1<<self._nvariables))
            bitset_copy(self._truth_table,(<BooleanFunction>x)._truth_table)
        else:
            raise TypeError, "unable to init the Boolean function"

    def __dealloc__(self):
        bitset_free(self._truth_table)

    def _repr_(self):
        """
        EXAMPLE::

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

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B=BooleanFunction([0, 1, 1, 0, 1, 0, 0, 0])
            sage: (~B).truth_table(format='int')
            (1, 0, 0, 1, 0, 1, 1, 1)
        """
        cdef BooleanFunction res=BooleanFunction(self.nvariables())
        bitset_complement(res._truth_table, self._truth_table)
        return res

    def __add__(self, BooleanFunction other):
        """
        Return the element wise sum of `self`and `other` which must have the same number of variables.

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: A=BooleanFunction([0, 1, 0, 1, 1, 0, 0, 1])
            sage: B=BooleanFunction([0, 1, 1, 0, 1, 0, 0, 0])
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
        if (self.nvariables() != other.nvariables() ):
            raise ValueError("the two Boolean functions must have the same number of variables")
        cdef BooleanFunction res = BooleanFunction(self)
        bitset_xor(res._truth_table, res._truth_table, other._truth_table)
        return res

    def __mul__(self, BooleanFunction other):
        """
        Return the elementwise multiplication of `self`and `other` which must have the same number of variables.

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: A=BooleanFunction([0, 1, 0, 1, 1, 0, 0, 1])
            sage: B=BooleanFunction([0, 1, 1, 0, 1, 0, 0, 0])
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
        if (self.nvariables() != other.nvariables() ):
            raise ValueError("the two Boolean functions must have the same number of variables")
        cdef BooleanFunction res = BooleanFunction(self)
        bitset_and(res._truth_table, res._truth_table, other._truth_table)
        return res

    def __or__(BooleanFunction self, BooleanFunction other):
        """
        Return the concatenation of `self` and `other` which must have the same number of variables.

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: A=BooleanFunction([0, 1, 0, 1])
            sage: B=BooleanFunction([0, 1, 1, 0])
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

        nb_limbs = self._truth_table.limbs
        if nb_limbs == 1:
            L = len(self)
            for i in range(L):
                res[i  ]=self[i]
                res[i+L]=other[i]
            return res

        memcpy(res._truth_table.bits             , self._truth_table.bits, nb_limbs * sizeof(unsigned long))
        memcpy(&(res._truth_table.bits[nb_limbs]),other._truth_table.bits, nb_limbs * sizeof(unsigned long))

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
        bitset_init(anf, (1<<self._nvariables))
        bitset_copy(anf, self._truth_table)
        reed_muller(anf.bits, ZZ(anf.limbs).exact_log(2))
        from sage.rings.polynomial.pbori import BooleanPolynomialRing
        R = BooleanPolynomialRing(self._nvariables,"x")
        G = R.gens()
        P = R(0)
        for 0 <= i < anf.limbs:
            if anf.bits[i]:
                inf = i*sizeof(long)*8
                sup = min( (i+1)*sizeof(long)*8 , (1<<self._nvariables) )
                for inf <= j < sup:
                    if bitset_in(anf,j):
                        m = R(1)
                        for 0 <= k < self._nvariables:
                            if (j>>k)&1:
                                m *= G[k]
                        P+=m
        bitset_free(anf)
        return P

    def nvariables(self):
        """
        The number of variables of this function.

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: BooleanFunction(4).nvariables()
            4
        """
        return self._nvariables

    def truth_table(self,format='bin'):
        """
        The truth table of the Boolean function.

        INPUT: a string representing the desired format, can be either

        - 'bin' (default) : we return a tuple of Boolean values
        - 'int' : we return a tuple of 0 or 1 values
        - 'hex' : we return a string representing the truth_table in hexadecimal

        EXAMPLE::

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

            sage: B.truth_table(format='oct')
            Traceback (most recent call last):
            ...
            ValueError: unknown output format
        """
        if format is 'bin':
            return tuple(self)
        if format is 'int':
            return tuple(map(int,self))
        if format is 'hex':
            S = ""
            S = ZZ(self.truth_table(),2).str(16)
            S = "0"*((1<<(self._nvariables-2)) - len(S)) + S
            for 1 <= i < self._truth_table.limbs:
                if sizeof(long)==4:
                    t = "%04x"%self._truth_table.bits[i]
                if sizeof(long)==8:
                    t = "%08x"%self._truth_table.bits[i]
                S = t + S
            return S
        raise ValueError, "unknown output format"

    def __len__(self):
        """
        Return the number of different input values.

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: len(BooleanFunction(4))
            16
        """
        return 2**self._nvariables

    def __cmp__(self, other):
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
        cdef BooleanFunction o=other
        return bitset_cmp(self._truth_table, o._truth_table)

    def __call__(self, x):
        """
        Return the value of the function for the given input.

        INPUT: either

        - a list - then all elements are evaluated as Booleans
        - an integer - then we consider its binary representation

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction([0,1,0,0])
            sage: B(1)
            1
            sage: B([1,0])
            1
            sage: B(7)
            Traceback (most recent call last):
            ...
            IndexError: index out of bound

        """
        if isinstance(x, (int,long,Integer)):
            if x > self._truth_table.size:
                raise IndexError, "index out of bound"
            return bitset_in(self._truth_table,x)
        elif isinstance(x, list):
            if len(x) != self._nvariables:
                raise ValueError, "bad number of inputs"
            return self(ZZ(map(bool,x),2))
        else:
            raise TypeError, "cannot apply Boolean function to provided element"

    def __iter__(self):
        """
        Iterate through the value of the function.

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction([0,1,1,0,1,0,1,0])
            sage: [int(b) for b in B]
            [0, 1, 1, 0, 1, 0, 1, 0]

        """
        return BooleanFunctionIterator(self)

    def _walsh_hadamard_transform_cached(self):
        """
        Return the cached Walsh Hadamard transform. *Unsafe*, no check.

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction(3)
            sage: W = B.walsh_hadamard_transform()
            sage: B._walsh_hadamard_transform_cached() is W
            True
        """
        return self._walsh_hadamard_transform

    def walsh_hadamard_transform(self):
        r"""
        Compute the Walsh Hadamard transform `W` of the function `f`.

        .. math:: W(j) = \sum_{i\in\{0,1\}^n} (-1)^{f(i)\oplus i \cdot j}

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: R.<x> = GF(2^3,'a')[]
            sage: B = BooleanFunction( x^3 )
            sage: B.walsh_hadamard_transform()
            (0, 4, 0, -4, 0, -4, 0, -4)
        """
        cdef long *temp

        if self._walsh_hadamard_transform is None:
            n =  self._truth_table.size
            temp = <long *>sage_malloc(sizeof(long)*n)

            for 0<= i < n:
                temp[i] = (bitset_in(self._truth_table,i)<<1)-1

            walsh_hadamard(temp, self._nvariables)
            self._walsh_hadamard_transform = tuple( [temp[i] for i in xrange(n)] )
            sage_free(temp)

        return self._walsh_hadamard_transform

    def absolute_walsh_spectrum(self):
        """
        Return the absolute Walsh spectrum fo the function.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("7969817CC5893BA6AC326E47619F5AD0")
            sage: sorted([_ for _ in B.absolute_walsh_spectrum().iteritems() ])
            [(0, 64), (16, 64)]

            sage: B = BooleanFunction("0113077C165E76A8")
            sage: B.absolute_walsh_spectrum()
            {8: 64}
        """
        d = {}
        for i in self.walsh_hadamard_transform():
            if d.has_key(abs(i)):
                d[abs(i)] += 1
            else:
                d[abs(i)] = 1
        return d

    def is_balanced(self):
        """
        Return True if the function takes the value True half of the time.

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
        Return True if the function is symmetric, i.e. invariant under
        permutation of its input bits. Another way to see it is that the
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
        cdef list T = [ self(2**i-1) for i in range(self._nvariables+1) ]
        for i in xrange(2**self._nvariables):
            if T[ hamming_weight_int(i) ] != bitset_in(self._truth_table, i):
                return False
        return True

    def nonlinearity(self):
        """
        Return the nonlinearity of the function. This is the distance
        to the linear functions, or the number of output ones need to
        change to obtain a linear function.

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction(5)
            sage: B[1] = B[3] = 1
            sage: B.nonlinearity()
            2
            sage: B = BooleanFunction("0113077C165E76A8")
            sage: B.nonlinearity()
            28
        """
        if self._nonlinearity is None:
            self._nonlinearity = ( (1<<self._nvariables) - max( [abs(w) for w in self.walsh_hadamard_transform()] ) ) >> 1
        return self._nonlinearity

    def is_bent(self):
        """
        Return True if the function is bent.

        EXAMPLE::

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
        `m` , if the output of the function is statistically
        independent of the combination of any m of its inputs.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("7969817CC5893BA6AC326E47619F5AD0")
            sage: B.correlation_immunity()
            2
        """
        cdef int c
        if self._correlation_immunity is None:
            c = self._nvariables
            W = self.walsh_hadamard_transform()
            for 0 < i < len(W):
                if (W[i] != 0):
                    c = min( c , hamming_weight_int(i) )
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
        Return the autocorrelation fo the function, defined by

        .. math:: \Delta_f(j) = \sum_{i\in\{0,1\}^n} (-1)^{f(i)\oplus f(i\oplus j)}.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("03")
            sage: B.autocorrelation()
            (8, 8, 0, 0, 0, 0, 0, 0)
        """
        cdef long *temp

        if self._autocorrelation is None:
            n =  self._truth_table.size
            temp = <long *>sage_malloc(sizeof(long)*n)
            W = self.walsh_hadamard_transform()

            for 0<= i < n:
                temp[i] = W[i]*W[i]

            walsh_hadamard(temp, self._nvariables)
            self._autocorrelation = tuple( [temp[i]>>self._nvariables for i in xrange(n)] )
            sage_free(temp)

        return self._autocorrelation

    def absolute_autocorrelation(self):
        """
        Return the absolute autocorrelation of the function.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("7969817CC5893BA6AC326E47619F5AD0")
            sage: sorted([ _ for _ in B.absolute_autocorrelation().iteritems() ])
            [(0, 33), (8, 58), (16, 28), (24, 6), (32, 2), (128, 1)]
        """
        d = {}
        for i in self.autocorrelation():
            if d.has_key(abs(i)):
                d[abs(i)] += 1
            else:
                d[abs(i)] = 1
        return d

    def absolut_indicator(self):
        """
        Return the absolut indicator of the function. Ths is the maximal absolut
        value of the autocorrelation.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("7969817CC5893BA6AC326E47619F5AD0")
            sage: B.absolut_indicator()
            32
        """
        if self._absolut_indicator is None:
            D = self.autocorrelation()
            self._absolut_indicator = max([ abs(a) for a in D[1:] ])
        return self._absolut_indicator

    def sum_of_square_indicator(self):
        """
        Return the sum of square indicator of the function.

        EXAMPLES::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction("7969817CC5893BA6AC326E47619F5AD0")
            sage: B.sum_of_square_indicator()
            32768
        """
        if self._sum_of_square_indicator is None:
            D = self.autocorrelation()
            self._sum_of_square_indicator = sum([ a**2 for a in D ])
        return self._sum_of_square_indicator

    def annihilator(self,d, dim = False):
        r"""
        Return (if it exists) an annihilator of the boolean function of
        degree at most `d`, that is a Boolean polynomial `g` such that

        .. math::

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
            set([0])
        """
        # NOTE: this is a toy implementation
        from sage.rings.polynomial.pbori import BooleanPolynomialRing
        R = BooleanPolynomialRing(self._nvariables,'x')
        G = R.gens()
        r = [R(1)]

        from sage.modules.all import vector
        s = vector(self.truth_table()).support()

        from sage.combinat.combination import Combinations
        from sage.misc.misc import prod

        from sage.matrix.constructor import Matrix
        from sage.rings.arith import binomial
        M = Matrix(GF(2),sum([binomial(self._nvariables,i) for i in xrange(d+1)]),len(s))

        for i in xrange(1,d+1):
            C = Combinations(self._nvariables,i)
            for c in C:
                r.append(prod([G[i] for i in c]))

        cdef BooleanFunction t

        for i,m in enumerate(r):
            t = BooleanFunction(m)
            for j,v in enumerate(s):
                M[i,j] = bitset_in(t._truth_table,v)

        kg = M.kernel().gens()

        if len(kg)>0:
            res = sum([kg[0][i]*ri for i,ri in enumerate(r)])
        else:
            res = None

        if dim:
            return res,len(kg)
        else:
            return res

    def algebraic_immunity(self, annihilator = False):
        """
        Returns the algebraic immunity of the Boolean function. This is the smallest
        integer `i` such that there exists a non trivial annihilator for `self` or `~self`.

        INPUT:

        - annihilator - a Boolean (default: False), if True, returns also an annihilator of minimal degree.

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
        for i in xrange(self._nvariables):
            for fun in [f,g]:
                A = fun.annihilator(i)
                if A is not None:
                    if annihilator:
                        return i,A
                    else:
                        return i
        raise ValueError, "you just found a bug!"

    def __setitem__(self, i, y):
        """
        Set a value of the function.

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B=BooleanFunction([0,0,1,1])
            sage: B[0]=1
            sage: B[2]=(3**17 == 9)
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

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B=BooleanFunction([0,1,1,1])
            sage: [ int(B[i]) for i in range(len(B)) ]
            [0, 1, 1, 1]
        """
        return self.__call__(i)

    def _clear_cache(self):
        """
        Clear cached values.

        EXAMPLE::

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
        self._absolut_indicator = None
        self._sum_of_square_indicator = None

    def __reduce__(self):
        """
        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction([0,1,1,0])
            sage: loads(dumps(B)) == B
            True
        """
        return unpickle_BooleanFunction, (self.truth_table(format='hex'),)

def unpickle_BooleanFunction(bool_list):
    """
    Specific function to unpickle Boolean functions.

    EXAMPLE::

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

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction(3)
            sage: type(B.__iter__())
            <type 'sage.crypto.boolean_function.BooleanFunctionIterator'>
        """
        self.f = f
        self.index = -1
        self.last = self.f._truth_table.size-1

    def __iter__(self):
        """
        Iterator through the values of a Boolean function.

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction(1)
            sage: [b for b in B] # indirect doctest
            [False, False]
        """
        return self

    def __next__(self):
        """
        Next value.

        EXAMPLE::

            sage: from sage.crypto.boolean_function import BooleanFunction
            sage: B = BooleanFunction(1)
            sage: I = B.__iter__()
            sage: I.next()
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
    Returns a random Boolean function with `n` variables.

    EXAMPLE::

        sage: from sage.crypto.boolean_function import random_boolean_function
        sage: B = random_boolean_function(9)
        sage: B.nvariables()
        9
        sage: B.nonlinearity()
        217                     # 32-bit
        222                     # 64-bit
    """
    from sage.misc.randstate import current_randstate
    r = current_randstate().python_random()
    cdef BooleanFunction B = BooleanFunction(n)
    cdef bitset_t T
    T[0] = B._truth_table[0]
    for 0 <= i < T.limbs:
        T.bits[i] = r.randrange(0,Integer(1)<<(sizeof(unsigned long)*8))
    return B
