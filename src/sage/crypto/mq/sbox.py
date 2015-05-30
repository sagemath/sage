r"""
S-Boxes and Their Algebraic Representations
"""

from sage.combinat.integer_vector import IntegerVectors
from sage.matrix.constructor import Matrix
from sage.misc.misc_c import prod as mul
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.element_base import is_FiniteFieldElement
from sage.rings.finite_rings.constructor import FiniteField as GF
from sage.rings.ideal import FieldIdeal, Ideal
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.sage_object import SageObject

class SBox(SageObject):
    r"""
    A substitution box or S-box is one of the basic components of
    symmetric key cryptography. In general, an S-box takes ``m`` input
    bits and transforms them into ``n`` output bits. This is called an
    ``mxn`` S-box and is often implemented as a lookup table. These
    S-boxes are carefully chosen to resist linear and differential
    cryptanalysis [Heys02]_.

    This module implements an S-box class which allows an algebraic
    treatment.

    EXAMPLE:

    We consider the S-box of the block cipher PRESENT [PRESENT07]_::

        sage: S = mq.SBox(12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2); S
        (12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2)
        sage: S(1)
        5

    Note that by default bits are interpreted in big endian
    order. This is not consistent with the rest of Sage, which has a
    strong bias towards little endian, but is consistent with most
    cryptographic literature::

        sage: S([0,0,0,1])
        [0, 1, 0, 1]

        sage: S = mq.SBox(12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2, big_endian=False)
        sage: S(1)
        5
        sage: S([0,0,0,1])
        [1, 1, 0, 0]


    Now we construct an ``SBox`` object for the 4-bit small scale AES
    S-Box (cf. :mod:`sage.crypto.mq.sr`)::

        sage: sr = mq.SR(1,1,1,4, allow_zero_inversions=True)
        sage: S = mq.SBox([sr.sub_byte(e) for e in list(sr.k)])
        sage: S
        (6, 5, 2, 9, 4, 7, 3, 12, 14, 15, 10, 0, 8, 1, 13, 11)

    REFERENCES:

    .. [Heys02] H. Heys *A Tutorial on Linear and Differential
      Cryptanalysis* ; 2002' available at
      http://www.engr.mun.ca/~howard/PAPERS/ldc_tutorial.pdf

    .. [PRESENT07] A. Bogdanov, L. Knudsen, G. Leander, C. Paar,
      A. Poschmann, M. Robshaw, Y. Seurin, C. Vikkelsoe *PRESENT: An
      Ultra-Lightweight Block Cipher*; in Proceedings of CHES 2007;
      LNCS 7427; pp. 450-466; Springer Verlag 2007; available at
      http://www.crypto.rub.de/imperia/md/content/texte/publications/conferences/present_ches2007.pdf
    """

    def __init__(self, *args,  **kwargs):
        """
        Construct a substitution box (S-box) for a given lookup table
        `S`.

        INPUT:

        - ``S`` - a finite iterable defining the S-box with integer or
          finite field elements

        - ``big_endian`` - controls whether bits shall be ordered in
          big endian order (default: ``True``)

        EXAMPLE:

        We construct a 3-bit S-box where e.g. the bits (0,0,1) are
        mapped to (1,1,1).::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3); S
            (7, 6, 0, 4, 2, 5, 1, 3)

            sage: S(0)
            7

        TESTS::

            sage: S = mq.SBox()
            Traceback (most recent call last):
            ...
            TypeError: No lookup table provided.
            sage: S = mq.SBox(1, 2, 3)
            Traceback (most recent call last):
            ...
            TypeError: Lookup table length is not a power of 2.
            sage: S = mq.SBox(5, 6, 0, 3, 4, 2, 1, 2)
            sage: S.n
            3
        """
        if "S" in kwargs:
            S = kwargs["S"]
        elif len(args) == 1:
            S = args[0]
        elif len(args) > 1:
            S = args
        else:
            raise TypeError("No lookup table provided.")

        _S = []
        for e in S:
            if is_FiniteFieldElement(e):
                e = e.polynomial().change_ring(ZZ).subs( e.parent().characteristic() )
            _S.append(e)
        S = _S

        if not ZZ(len(S)).is_power_of(2):
            raise TypeError("Lookup table length is not a power of 2.")
        self._S = S

        self.m = ZZ(len(S)).exact_log(2)
        self.n = ZZ(max(S)).nbits()
        self._F = GF(2)
        self._big_endian = kwargs.get("big_endian",True)

    def _repr_(self):
        """
        EXAMPLE::

            sage: mq.SBox(7,6,0,4,2,5,1,3) #indirect doctest
            (7, 6, 0, 4, 2, 5, 1, 3)
        """
        return "(" + ", ".join(map(str,list(self))) + ")"

    def __len__(self):
        """
        Return the length of input bit strings.

        EXAMPLE::

            sage: len(mq.SBox(7,6,0,4,2,5,1,3))
            3
        """
        return self.m

    def __cmp__(self, other):
        """
        S-boxes are considered to be equal if all construction
        parameters match.

        EXAMPLE::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3)
            sage: loads(dumps(S)) == S
            True
        """
        return cmp((self._S,self._big_endian), (other._S,self._big_endian))

    def to_bits(self, x, n=None):
        """
        Return bitstring of length ``n`` for integer ``x``. The
        returned bitstring is guaranteed to have length ``n``.

        INPUT:

        - ``x`` - an integer

        - ``n`` - bit length (optional)

        EXAMPLE::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3)
            sage: S.to_bits(6)
            [1, 1, 0]

            sage: S.to_bits( S(6) )
            [0, 0, 1]

            sage: S( S.to_bits( 6 ) )
            [0, 0, 1]
        """
        if n is None and self.m == self.n:
            n = self.n

        if self._big_endian:
            swp = lambda x: list(reversed(x))
        else:
            swp = lambda x: x
        return swp(self._rpad( map(self._F,ZZ(x).digits(2)), n ))

    def from_bits(self, x, n=None):
        """
        Return integer for bitstring ``x`` of length ``n``.

        INPUT:

        - ``x`` - a bitstring

        - ``n`` - bit length (optional)

        EXAMPLE::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3)
            sage: S.from_bits( [1,1,0])
            6

            sage: S( S.from_bits( [1,1,0] ) )
            1
            sage: S.from_bits( S( [1,1,0] ) )
            1
        """
        if n is None and self.m == self.n:
            n = self.m

        if self._big_endian:
            swp = lambda x: list(reversed(x))
        else:
            swp = lambda x: x

        return ZZ( map(ZZ, self._rpad(swp(x), n)), 2)

    def _rpad(self,x, n=None):
        """
        Right pads ``x`` such that ``len(x) == n``.

        EXAMPLE::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3)
            sage: S._rpad([1,1])
            [1, 1, 0]
        """
        if n is None and self.m == self.n:
            n = self.n
        return  x + [self._F(0)]*(n-len(x))

    def __call__(self, X):
        """
        Apply substitution to ``X``.

        If ``X`` is a list, it is interpreted as a sequence of bits
        depending on the bit order of this S-box.

        INPUT:

        - ``X`` - either an integer, a tuple of `\GF{2}` elements of
          length ``len(self)`` or a finite field element in
          `\GF{2^n}`. As a last resort this function tries to convert
          ``X`` to an integer.

        EXAMPLE::

            sage: S = mq.SBox([7,6,0,4,2,5,1,3])
            sage: S(7)
            3

            sage: S((0,2,3))
            [0, 1, 1]

            sage: S[0]
            7

            sage: S[(0,0,1)]
            [1, 1, 0]

            sage: k.<a> = GF(2^3)
            sage: S(a^2)
            a

            sage: S(QQ(3))
            4

            sage: S([1]*10^6)
            Traceback (most recent call last):
            ...
            TypeError: Cannot apply SBox to provided element.

            sage: S(1/2)
            Traceback (most recent call last):
            ...
            TypeError: Cannot apply SBox to 1/2.

            sage: S = mq.SBox(3, 0, 1, 3, 1, 0, 2, 2)
            sage: S(0)
            3
            sage: S([0,0,0])
            [1, 1]
        """
        if isinstance(X, (int, long, Integer)):
            return self._S[ZZ(X)]

        try:
            from sage.modules.free_module_element import vector
            K = X.parent()
            if K.order() == 2**self.n:
                X = vector(X)
            else:
                raise TypeError
            if not self._big_endian:
                X = list(reversed(X))
            else:
                X = list(X)
            X = ZZ(map(ZZ,X),2)
            out =  self.to_bits(self._S[X], self.n)
            if self._big_endian:
                out = list(reversed(out))
            return K(vector(GF(2),out))
        except (AttributeError, TypeError):
            pass

        try:
            if len(X) == self.m:
                if self._big_endian:
                    X = list(reversed(X))
                X = ZZ(map(ZZ,X),2)
                out =  self._S[X]
                return self.to_bits(out,self.n)
        except TypeError:
            pass

        try:
            return self._S[ZZ(X)]
        except TypeError:
            pass

        if len(str(X)) > 50:
            raise TypeError("Cannot apply SBox to provided element.")
        else:
            raise TypeError("Cannot apply SBox to %s."%(X,))

    def __getitem__(self, X):
        """
        See  :meth:`SBox.__call__`.

        EXAMPLE::

            sage: S = mq.SBox([7,6,0,4,2,5,1,3])
            sage: S[7]
            3
        """
        return self(X)

    def is_permutation(self):
        """
        Return ``True`` if this S-Box is a permutation.

        EXAMPLE::

             sage: S = mq.SBox(7,6,0,4,2,5,1,3)
             sage: S.is_permutation()
             True

             sage: S = mq.SBox(3,2,0,0,2,1,1,3)
             sage: S.is_permutation()
             False
        """
        if self.m != self.n:
            return False
        return len(set([self(i) for i in range(2**self.m)])) == 2**self.m

    def __iter__(self):
        """
        EXAMPLE::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3)
            sage: [e for e in S]
            [7, 6, 0, 4, 2, 5, 1, 3]
        """
        for i in xrange(2**self.m):
            yield self(i)

    def difference_distribution_matrix(self):
        """
        Return difference distribution matrix ``A`` for this S-box.

        The rows of ``A`` encode the differences ``Delta I`` of the
        input and the columns encode the difference ``Delta O`` for
        the output. The bits are ordered according to the endianess of
        this S-box. The value at ``A[Delta I,Delta O]`` encodes how
        often ``Delta O`` is the actual output difference given
        ``Delta I`` as input difference.

        See [Heys02]_ for an introduction to differential
        cryptanalysis.

        EXAMPLE::

           sage: S = mq.SBox(7,6,0,4,2,5,1,3)
           sage: S.difference_distribution_matrix()
           [8 0 0 0 0 0 0 0]
           [0 2 2 0 2 0 0 2]
           [0 0 2 2 0 0 2 2]
           [0 2 0 2 2 0 2 0]
           [0 2 0 2 0 2 0 2]
           [0 0 2 2 2 2 0 0]
           [0 2 2 0 0 2 2 0]
           [0 0 0 0 2 2 2 2]
        """
        m = self.m
        n = self.n

        nrows = 1<<m
        ncols = 1<<n

        A = Matrix(ZZ, nrows, ncols)

        for i in range(nrows):
            si = self(i)
            for di in range(nrows):
                A[ di , si^self(i^di)] += 1
        return A

    def maximal_difference_probability_absolute(self):
        """
        Return the difference probability of the difference with the
        highest probability in absolute terms, i.e. how often it
        occurs in total.

        EXAMPLE::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3)
            sage: S.maximal_difference_probability_absolute()
            2

        .. note::

          This code is mainly called internally.
        """
        A = self.difference_distribution_matrix().__copy__()
        A[0,0] = 0
        return max(map(abs, A.list()))

    def maximal_difference_probability(self):
        r"""
        Return the difference probability of the difference with the
        highest probability in the range between 0.0 and 1.0
        indicating 0\% or 100\% respectively.

        EXAMPLE::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3)
            sage: S.maximal_difference_probability()
            0.25
        """
        return self.maximal_difference_probability_absolute()/(2.0**self.n)

    def linear_approximation_matrix(self):
        """
        Return linear approximation matrix ``A`` for this S-box.

        Let ``i_b`` be the ``b``-th bit of ``i`` and ``o_b`` the
        ``b``-th bit of ``o``. Then ``v = A[i,o]`` encodes the bias of
        the equation ``sum( i_b * x_i ) = sum( o_b * y_i )`` if
        ``x_i`` and ``y_i`` represent the input and output variables
        of the S-box.

        See [Heys02]_ for an introduction to linear cryptanalysis.

        EXAMPLE::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3)
            sage: S.linear_approximation_matrix()
            [ 4  0  0  0  0  0  0  0]
            [ 0  0  0  0  2  2  2 -2]
            [ 0  0 -2 -2 -2  2  0  0]
            [ 0  0 -2  2  0  0 -2 -2]
            [ 0  2  0  2 -2  0  2  0]
            [ 0 -2  0  2  0  2  0  2]
            [ 0 -2 -2  0  0 -2  2  0]
            [ 0 -2  2  0 -2  0  0 -2]

        According to this matrix the first bit of the input is equal
        to the third bit of the output 6 out of 8 times::

            sage: for i in srange(8): print S.to_bits(i)[0] == S.to_bits(S(i))[2]
            False
            True
            True
            True
            False
            True
            True
            True
        """
        try:
            return self._linear_approximation_matrix
        except AttributeError:
            pass

        m = self.m
        n = self.n

        nrows = 1<<m
        ncols = 1<<n

        from sage.crypto.boolean_function import BooleanFunction

        B = BooleanFunction(self.m)
        L = []
        for j in xrange(ncols):
            for i in xrange(nrows):
                B[i] = (self(i)&j).popcount()
            L.append(B.walsh_hadamard_transform())

        A = Matrix(ZZ, ncols, nrows, L)
        A = -A.transpose()/2

        self._linear_approximation_matrix = A
        return A

    def maximal_linear_bias_absolute(self):
        """
        Return maximal linear bias, i.e. how often the linear
        approximation with the highest bias is true or false minus
        `2^{n-1}`.

        EXAMPLE::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3)
            sage: S.maximal_linear_bias_absolute()
            2
        """
        A = self.linear_approximation_matrix().__copy__()
        A[0,0] = 0
        return max(map(abs, A.list()))

    def maximal_linear_bias_relative(self):
        """
        Return maximal bias of all linear approximations of this
        S-box.

        EXAMPLE::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3)
            sage: S.maximal_linear_bias_relative()
            0.25
        """
        return self.maximal_linear_bias_absolute()/(2.0**self.m)

    def ring(self):
        """
        Create, return and cache a polynomial ring for S-box
        polynomials.

        EXAMPLE::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3)
            sage: S.ring()
            Multivariate Polynomial Ring in x0, x1, x2, y0, y1, y2 over Finite Field of size 2
        """
        try:
            return self._ring
        except AttributeError:
            pass

        m = self.m
        n = self.n

        X = range(m)
        Y = range(n)
        self._ring = PolynomialRing(self._F, m+n, ["x%d"%i for i in X] + ["y%d"%i for i in Y])
        return self._ring

    def solutions(self, X=None, Y=None):
        """
        Return a dictionary of solutions to this S-box.

        INPUT:

        - ``X`` - input variables (default: ``None``)

        - ``Y`` - output variables (default: ``None``)

        EXAMPLE::

            sage: S = mq.SBox([7,6,0,4,2,5,1,3])
            sage: F = S.polynomials()
            sage: s = S.solutions()
            sage: any(f.subs(_s) for f in F for _s in s)
            False
        """
        if X is None and Y is None:
            P = self.ring()
            gens = P.gens()
        else:
            P = X[0].parent()
            gens = X + Y

        m = self.m
        n = self.n

        solutions = []
        for i in range(1<<m):
            solution = self.to_bits(i,m) + self( self.to_bits(i,m) )
            solutions.append( dict(zip(gens, solution)) )

        return solutions

    def polynomials(self, X=None, Y=None, degree=2, groebner=False):
        """
        Return a list of polynomials satisfying this S-box.

        First, a simple linear fitting is performed for the given
        ``degree`` (cf. for example [BC03]_). If ``groebner=True`` a
        Groebner basis is also computed for the result of that
        process.

        INPUT:

        - ``X`` - input variables

        - ``Y`` - output variables

        - ``degree`` - integer > 0 (default: ``2``)

        - ``groebner`` - calculate a reduced Groebner basis of the
          spanning polynomials to obtain more polynomials (default:
          ``False``)

        EXAMPLES::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3)
            sage: P = S.ring()

        By default, this method returns an indirect representation::

            sage: S.polynomials()
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
            sage: S.polynomials([x0,x1,x2],[y0,y1,y2], groebner=True)
            [y0 + x0*x1 + x0*x2 + x0 + x1*x2 + x1 + 1,
             y1 + x0*x2 + x1 + 1,
             y2 + x0 + x1*x2 + x1 + x2 + 1]

        REFERENCES:

        .. [BC03] A. Biryukov and C. D. Canniere *Block Ciphers and
          Systems of Quadratic Equations*; in Proceedings of Fast
          Software Encryption 2003; LNCS 2887; pp. 274-289,
          Springer-Verlag 2003.
        """
        def nterms(nvars, deg):
            """
            Return the number of monomials possible up to a given
            degree.

            INPUT:

            - ``nvars`` - number of variables

            - ``deg`` - degree

            TESTS::

                sage: S = mq.SBox(7,6,0,4,2,5,1,3)
                sage: F = S.polynomials(degree=3) # indirect doctest
            """
            total = 1
            divisor = 1
            var_choices = 1

            for d in xrange(1, deg+1):
                var_choices *= (nvars - d + 1)
                divisor *= d
                total += var_choices/divisor
            return total

        m = self.m
        n = self.n
        F = self._F

        if X is None and Y is None:
            P = self.ring()
            X = P.gens()[:m]
            Y = P.gens()[m:]
        else:
            P = X[0].parent()

        gens = X+Y

        bits = []
        for i in range(1<<m):
            bits.append( self.to_bits(i,m) + self(self.to_bits(i,m)) )

        ncols = (1<<m)+1

        A = Matrix(P, nterms(m+n, degree), ncols)

        exponents = []
        for d in range(degree+1):
            exponents += IntegerVectors(d, max_length=m+n, min_length=m+n, min_part=0, max_part=1).list()

        row = 0
        for exponent in exponents:
            A[row,ncols-1] = mul([gens[i]**exponent[i] for i in range(len(exponent))])
            for col in range(1<<m):
                A[row,col] = mul([bits[col][i] for i in range(len(exponent)) if exponent[i]])
            row +=1

        for c in range(ncols):
            A[0,c] = 1

        RR = A.echelon_form(algorithm='row_reduction')

        # extract spanning stet
        gens = (RR.column(ncols-1)[1<<m:]).list()

        if not groebner:
            return gens

        FI = set(FieldIdeal(P).gens())
        I = Ideal(gens + list(FI))
        gb = I.groebner_basis()

        gens = []
        for f in gb:
            if f not in FI: # filter out field equations
                gens.append(f)
        return gens

    def interpolation_polynomial(self, k=None):
        r"""
        Return a univariate polynomial over an extension field
        representing this S-box.

        If ``m`` is the input length of this S-box then the extension
        field is of degree ``m``.

        If the output length does not match the input length then a
        ``TypeError`` is raised.

        INPUT:

        - ``k`` - an instance of `\GF{2^m}` (default: ``None``)

        EXAMPLE::

            sage: S = mq.SBox(7,6,0,4,2,5,1,3)
            sage: f = S.interpolation_polynomial()
            sage: f
            x^6 + a*x^5 + (a + 1)*x^4 + (a^2 + a + 1)*x^3
              + (a^2 + 1)*x^2 + (a + 1)*x + a^2 + a + 1

            sage: a = f.base_ring().gen()

            sage: f(0), S(0)
            (a^2 + a + 1, 7)

            sage: f(a^2 + 1), S(5)
            (a^2 + 1, 5)
        """
        if self.m != self.n:
            raise TypeError("Lagrange interpolation only supported if self.m == self.n.")

        if k is None:
            k = GF(2**self.m,'a')
        l = []
        for i in xrange(2**self.m):
            i = self.to_bits(i, self.m)
            o = self(i)
            if self._big_endian:
                i = reversed(i)
                o = reversed(o)
            l.append( (k(vector(i)), k(vector(o))) )

        P = PolynomialRing(k,'x')
        return P.lagrange_polynomial(l)

    def cnf(self, xi=None, yi=None, format=None):
        """
        Return a representation of this S-Box in conjunctive normal
        form.

        This function examines the truth tables for each output bit of
        the S-Box and thus has complexity `n * 2^m` for an ``m x n``
        S-Box.

        INPUT:

        - ``xi`` - indices for the input variables (default: ``1...m``)

        - ``yi`` - indices for the output variables (default: ``m+1 ... m+n``)

        - ``format`` - output format, see below (default: ``None``)

        FORMATS:

        - ``None`` - return a list of tuples of integers where each
          tuple represents a clause, the absolute value of an integer
          represents a variable and the sign of an integer indicates
          inversion.

        - ``symbolic`` - a string that can be parsed by the
          ``SymbolicLogic`` package.

        - ``dimacs`` - a string in DIMACS format which is the gold
          standard for SAT-solver input (cf. http://www.satlib.org/).

        - ``dimacs_headless`` - a string in DIMACS format, but without
          the header. This is useful for concatenation of outputs.

        EXAMPLE:

        We give a very small example to explain the output format::

            sage: S = mq.SBox(1,2,0,3); S
            (1, 2, 0, 3)
            sage: cnf = S.cnf(); cnf
            [(1, 2, -3),  (1, 2, 4),
             (1, -2, 3),  (1, -2, -4),
             (-1, 2, -3), (-1, 2, -4),
             (-1, -2, 3), (-1, -2, 4)]

        This output completely describes the S-Box. For instance, we
        can check that ``S([0,1]) -> [1,0]`` satisfies every clause if
        the first input bit corresponds to the index ``1`` and the
        last output bit corresponds to the index ``3`` in the
        output.

        We can convert this representation to the DIMACS format::

            sage: print S.cnf(format='dimacs')
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

            sage: print S.cnf(format='dimacs_headless')
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

            sage: print S.cnf(xi=[10,20],yi=[30,40], format='dimacs_headless')
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

            sage: S = mq.SBox(1,2,0,3, big_endian=False); S
            (1, 2, 0, 3)
            sage: cnf = S.cnf(); cnf
            [(1, 2, -4), (1, 2, 3),
             (-1, 2, 4), (-1, 2, -3),
             (1, -2, -4), (1, -2, -3),
             (-1, -2, 4), (-1, -2, 3)]

        S-Boxes with m!=n also work:

            sage: o = range(8) + range(8)
            sage: shuffle(o)
            sage: S = mq.SBox(o)
            sage: S.is_permutation()
            False

            sage: len(S.cnf()) == 3*2^4
            True


        TESTS:

            sage: S = mq.SBox(1,2,0,3, big_endian=False)
            sage: S.cnf([1000,1001,1002], [2000,2001,2002])
            Traceback (most recent call last):
            ...
            TypeError: first arg required to have length 2, got 3 instead.
        """
        m, n = self.m, self.n

        if xi is None:
            xi = [i+1 for i in range(m)]

        if yi is None:
            yi = [m+i+1 for i in range(n)]

        if len(xi) != m:
            raise TypeError("first arg required to have length %d, got %d instead."%(m,len(xi)))

        if len(yi) != n:
            raise TypeError("second arg required to have length %d, got %d instead."%(n,len(yi)))

        output_bits = range(n)
        if not self._big_endian:
            output_bits = list(reversed(output_bits))

        C = [] # the set of clauses
        for e in xrange(2**m):
            x = self.to_bits(e, m)
            y = self(x) # evaluate at x
            for output_bit in output_bits: # consider each bit
                clause = [(-1)**(int(v)) * i for v,i in zip(x, xi)]
                clause.append( (-1)**(1-int(y[output_bit])) *  yi[output_bit] )
                C.append(tuple(clause))

        if format is None:
            return C
        elif format == 'symbolic':
            gd = self.ring().gens()
            formula = []
            for clause in C:
                clause = "|".join([str(gd[abs(v)-1]).replace("-","~") for v in clause])
                formula.append("("+clause+")")
            return " & ".join(formula)

        elif format.startswith('dimacs'):
            if format == "dimacs_headless":
                header = ""
            else:
                header = "p cnf %d %d\n"%(m+n,len(C))
            values = " 0\n".join([" ".join(map(str,line)) for line in C])
            return header + values + " 0\n"
        else:
            raise ValueError("Format '%s' not supported."%(format,))
