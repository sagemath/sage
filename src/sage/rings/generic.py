# sage_setup: distribution = sagemath-categories
r"""
Generic data structures and algorithms for rings

AUTHORS:

- Lorenz Panny (2022): :class:`ProductTree`, :func:`prod_with_derivative`
"""

from sage.misc.misc_c import prod


class ProductTree:
    r"""
    A simple binary product tree, i.e., a tree of ring elements in
    which every node equals the product of its children.
    (In particular, the *root* equals the product of all *leaves*.)

    Product trees are a very useful building block for fast computer
    algebra. For example, a quasilinear-time Discrete Fourier Transform
    (the famous *Fast* Fourier Transform) can be implemented as follows
    using the :meth:`remainders` method of this class::

        sage: # needs sage.rings.finite_rings
        sage: from sage.rings.generic import ProductTree
        sage: F = GF(65537)
        sage: a = F(1111)
        sage: assert a.multiplicative_order() == 1024
        sage: R.<x> = F[]
        sage: ms = [x - a^i for i in range(1024)]               # roots of unity
        sage: ys = [F.random_element() for _ in range(1024)]    # input vector
        sage: tree = ProductTree(ms)
        sage: zs = tree.remainders(R(ys))                       # compute FFT!
        sage: zs == [R(ys) % m for m in ms]
        True

    Similarly, the :meth:`interpolation` method can be used to implement
    the inverse Fast Fourier Transform::

        sage: tree.interpolation(zs).padded_list(len(ys)) == ys                         # needs sage.rings.finite_rings
        True

    This class encodes the tree as *layers*: Layer `0` is just a tuple
    of the leaves. Layer `i+1` is obtained from layer `i` by replacing
    each pair of two adjacent elements by their product, starting from
    the left. (If the length is odd, the unpaired element at the end is
    simply copied as is.) This iteration stops as soon as it yields a
    layer containing only a single element (the root).

    .. NOTE::

        Use this class if you need the :meth:`remainders` method.
        To compute just the product, :func:`prod` is likely faster.

    INPUT:

    - ``leaves`` -- an iterable of elements in a common ring

    EXAMPLES::

        sage: from sage.rings.generic import ProductTree
        sage: R.<x> = GF(101)[]
        sage: vs = [x - i for i in range(1,10)]
        sage: tree = ProductTree(vs)
        sage: tree.root()
        x^9 + 56*x^8 + 62*x^7 + 44*x^6 + 47*x^5 + 42*x^4 + 15*x^3 + 11*x^2 + 12*x + 13
        sage: tree.remainders(x^7 + x + 1)
        [3, 30, 70, 27, 58, 72, 98, 98, 23]
        sage: tree.remainders(x^100)
        [1, 1, 1, 1, 1, 1, 1, 1, 1]

    ::

        sage: # needs sage.libs.pari
        sage: vs = prime_range(100)
        sage: tree = ProductTree(vs)
        sage: tree.root().factor()
        2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47 * 53 * 59 * 61 * 67 * 71 * 73 * 79 * 83 * 89 * 97
        sage: tree.remainders(3599)
        [1, 2, 4, 1, 2, 11, 12, 8, 11, 3, 3, 10, 32, 30, 27, 48, 0, 0, 48, 49, 22, 44, 30, 39, 10]

    We can access the individual layers of the tree::

        sage: tree.layers                                                               # needs sage.libs.pari
        [(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97),
         (6, 35, 143, 323, 667, 1147, 1763, 2491, 3599, 4757, 5767, 7387, 97),
         (210, 46189, 765049, 4391633, 17120443, 42600829, 97),
         (9699690, 3359814435017, 729345064647247, 97),
         (32589158477190044730, 70746471270782959),
         (2305567963945518424753102147331756070,)]
    """
    def __init__(self, leaves):
        r"""
        Initialize a product tree having the given ring elements
        as its leaves.

        EXAMPLES::

            sage: from sage.rings.generic import ProductTree
            sage: vs = prime_range(100)                                                 # needs sage.libs.pari
            sage: tree = ProductTree(vs)                                                # needs sage.libs.pari
        """
        V = tuple(leaves)
        self.layers = [V]
        while len(V) > 1:
            V = tuple(prod(V[i:i+2]) for i in range(0,len(V),2))
            self.layers.append(V)

    def __len__(self):
        r"""
        Return the number of leaves of this product tree.

        EXAMPLES::

            sage: from sage.rings.generic import ProductTree
            sage: R.<x> = GF(101)[]
            sage: vs = [x - i for i in range(1,10)]
            sage: tree = ProductTree(vs)
            sage: len(tree)
            9
            sage: len(tree) == len(vs)
            True
            sage: len(tree.remainders(x^2))
            9
        """
        return len(self.layers[0])

    def __iter__(self):
        r"""
        Return an iterator over the leaves of this product tree.

        EXAMPLES::

            sage: from sage.rings.generic import ProductTree
            sage: R.<x> = GF(101)[]
            sage: vs = [x - i for i in range(1,10)]
            sage: tree = ProductTree(vs)
            sage: next(iter(tree)) == vs[0]
            True
            sage: list(tree) == vs
            True
        """
        return iter(self.layers[0])

    def root(self):
        r"""
        Return the value at the root of this product tree (i.e., the product of all leaves).

        EXAMPLES::

            sage: from sage.rings.generic import ProductTree
            sage: R.<x> = GF(101)[]
            sage: vs = [x - i for i in range(1,10)]
            sage: tree = ProductTree(vs)
            sage: tree.root()
            x^9 + 56*x^8 + 62*x^7 + 44*x^6 + 47*x^5 + 42*x^4 + 15*x^3 + 11*x^2 + 12*x + 13
            sage: tree.root() == prod(vs)
            True
        """
        assert len(self.layers[-1]) == 1
        return self.layers[-1][0]

    def leaves(self):
        r"""
        Return a tuple containing the leaves of this product tree.

        EXAMPLES::

            sage: from sage.rings.generic import ProductTree
            sage: R.<x> = GF(101)[]
            sage: vs = [x - i for i in range(1,10)]
            sage: tree = ProductTree(vs)
            sage: tree.leaves()
            (x + 100, x + 99, x + 98, ..., x + 93, x + 92)
            sage: tree.leaves() == tuple(vs)
            True
        """
        return self.layers[0]

    def remainders(self, x):
        r"""
        Given a value `x`, return a list of all remainders of `x`
        modulo the leaves of this product tree.

        The base ring must support the ``%`` operator for this
        method to work.

        INPUT:

        - ``x`` -- an element of the base ring of this product tree

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: from sage.rings.generic import ProductTree
            sage: vs = prime_range(100)
            sage: tree = ProductTree(vs)
            sage: n = 1085749272377676749812331719267
            sage: tree.remainders(n)
            [1, 1, 2, 1, 9, 1, 7, 15, 8, 20, 15, 6, 27, 11, 2, 6, 0, 25, 49, 5, 51, 4, 19, 74, 13]
            sage: [n % v for v in vs]
            [1, 1, 2, 1, 9, 1, 7, 15, 8, 20, 15, 6, 27, 11, 2, 6, 0, 25, 49, 5, 51, 4, 19, 74, 13]
        """
        X = [x]
        for V in reversed(self.layers):
            X = [X[i // 2] % V[i] for i in range(len(V))]
        return X

    _crt_bases = None

    def interpolation(self, xs):
        r"""
        Given a sequence ``xs`` of values, one per leaf, return a
        single element `x` which is congruent to the `i`\th value in
        ``xs`` modulo the `i`\th leaf, for all `i`.

        This is an explicit version of the Chinese remainder theorem;
        see also :meth:`CRT`. Using this product tree is faster for
        repeated calls since the required CRT bases are cached after
        the first run.

        The base ring must support the :func:`xgcd` function for this
        method to work.

        EXAMPLES::

            sage: from sage.rings.generic import ProductTree
            sage: vs = prime_range(100)
            sage: tree = ProductTree(vs)
            sage: tree.interpolation([1, 1, 2, 1, 9, 1, 7, 15, 8, 20, 15, 6, 27, 11, 2, 6, 0, 25, 49, 5, 51, 4, 19, 74, 13])
            1085749272377676749812331719267

        This method is faster than :func:`CRT` for repeated calls with
        the same moduli::

            sage: vs = prime_range(1000,2000)
            sage: rs = lambda: [randrange(1,100) for _ in vs]
            sage: tree = ProductTree(vs)
            sage: %timeit CRT(rs(), vs)             # not tested
            372 µs ± 3.34 µs per loop (mean ± std. dev. of 7 runs, 1,000 loops each)
            sage: %timeit tree.interpolation(rs())  # not tested
            146 µs ± 479 ns per loop (mean ± std. dev. of 7 runs, 10,000 loops each)
        """
        if self._crt_bases is None:
            from sage.arith.misc import CRT_basis
            self._crt_bases = []
            for V in self.layers[:-1]:
                B = tuple(CRT_basis(V[i:i+2]) for i in range(0, len(V), 2))
                self._crt_bases.append(B)
        if len(xs) != len(self.layers[0]):
            raise ValueError('number of given elements must equal the number of leaves')
        for basis, layer in zip(self._crt_bases, self.layers[1:]):
            xs = [sum(c*x for c, x in zip(cs, xs[2*i:2*i+2])) % mod
                  for i, (cs, mod) in enumerate(zip(basis, layer))]
        assert len(xs) == 1
        return xs[0]


def prod_with_derivative(pairs):
    r"""
    Given an iterable of pairs `(f, \partial f)` of ring elements,
    return the pair `(\prod f, \partial \prod f)`, assuming `\partial`
    is an operator obeying the standard product rule.

    This function is entirely algebraic, hence still works when the
    elements `f` and `\partial f` are all passed through some ring
    homomorphism first. One particularly useful instance of this is
    evaluating the derivative of a product of polynomials at a point
    without fully expanding the product; see the second example below.

    INPUT:

    - ``pairs`` -- an iterable of tuples `(f, \partial f)` of elements
      of a common ring

    ALGORITHM: Repeated application of the product rule.

    EXAMPLES::

        sage: from sage.rings.generic import prod_with_derivative
        sage: R.<x> = ZZ[]
        sage: fs = [x^2 + 2*x + 3, 4*x + 5, 6*x^7 + 8*x + 9]
        sage: prod(fs)
        24*x^10 + 78*x^9 + 132*x^8 + 90*x^7 + 32*x^4 + 140*x^3 + 293*x^2 + 318*x + 135
        sage: prod(fs).derivative()
        240*x^9 + 702*x^8 + 1056*x^7 + 630*x^6 + 128*x^3 + 420*x^2 + 586*x + 318
        sage: F, dF = prod_with_derivative((f, f.derivative()) for f in fs)
        sage: F
        24*x^10 + 78*x^9 + 132*x^8 + 90*x^7 + 32*x^4 + 140*x^3 + 293*x^2 + 318*x + 135
        sage: dF
        240*x^9 + 702*x^8 + 1056*x^7 + 630*x^6 + 128*x^3 + 420*x^2 + 586*x + 318

    The main reason for this function to exist is that it allows us to
    *evaluate* the derivative of a product of polynomials at a point
    `\alpha` without ever fully expanding the product *as a polynomial*::

        sage: alpha = 42
        sage: F(alpha)
        442943981574522759
        sage: dF(alpha)
        104645261461514994
        sage: us = [f(alpha) for f in fs]
        sage: vs = [f.derivative()(alpha) for f in fs]
        sage: prod_with_derivative(zip(us, vs))
        (442943981574522759, 104645261461514994)
    """
    class _aux:
        def __init__(self, f, df):
            self.f, self.df = f, df

        def __mul__(self, other):
            return _aux(self.f * other.f, self.df * other.f + self.f * other.df)

        def __iter__(self):
            yield self.f
            yield self.df

    return tuple(prod(_aux(*tup) for tup in pairs))
