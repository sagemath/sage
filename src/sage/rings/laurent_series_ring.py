"""
Laurent Series Rings

EXAMPLES::

    sage: R = LaurentSeriesRing(QQ, "x")
    sage: R.base_ring()
    Rational Field
    sage: S = LaurentSeriesRing(GF(17)['x'], 'y')
    sage: S
    Laurent Series Ring in y over Univariate Polynomial Ring in x over
    Finite Field of size 17
    sage: S.base_ring()
    Univariate Polynomial Ring in x over Finite Field of size 17

.. SEEALSO::

    * :func:`sage.misc.defaults.set_series_precision`
"""
from __future__ import print_function, absolute_import

import weakref

from . import laurent_series_ring_element

from . import polynomial
from . import ring

from sage.libs.pari.all import pari_gen
from sage.categories.fields import Fields
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationFields

laurent_series = {}
def LaurentSeriesRing(base_ring, name=None, names=None, default_prec=None, sparse=False):
    """
    EXAMPLES::

        sage: R = LaurentSeriesRing(QQ, 'x'); R
        Laurent Series Ring in x over Rational Field
        sage: x = R.0
        sage: g = 1 - x + x^2 - x^4 +O(x^8); g
        1 - x + x^2 - x^4 + O(x^8)
        sage: g = 10*x^(-3) + 2006 - 19*x + x^2 - x^4 +O(x^8); g
        10*x^-3 + 2006 - 19*x + x^2 - x^4 + O(x^8)

    You can also use more mathematical notation when the base is a
    field::

        sage: Frac(QQ[['x']])
        Laurent Series Ring in x over Rational Field
        sage: Frac(GF(5)['y'])
        Fraction Field of Univariate Polynomial Ring in y over Finite Field of size 5

    Here the fraction field is not just the Laurent series ring, so you
    can't use the ``Frac`` notation to make the Laurent
    series ring.

    ::

        sage: Frac(ZZ[['t']])
        Fraction Field of Power Series Ring in t over Integer Ring

    Laurent series rings are determined by their variable and the base
    ring, and are globally unique.

    ::

        sage: K = Qp(5, prec = 5)
        sage: L = Qp(5, prec = 200)
        sage: R.<x> = LaurentSeriesRing(K)
        sage: S.<y> = LaurentSeriesRing(L)
        sage: R is S
        False
        sage: T.<y> = LaurentSeriesRing(Qp(5,prec=200))
        sage: S is T
        True
        sage: W.<y> = LaurentSeriesRing(Qp(5,prec=199))
        sage: W is T
        False

    TESTS:

    Check if changing global series precision does it right (and
    that :trac:`17955` is fixed)::

        sage: set_series_precision(3)
        sage: R.<x> = LaurentSeriesRing(ZZ)
        sage: 1/(1 - 2*x)
        1 + 2*x + 4*x^2 + O(x^3)
        sage: set_series_precision(5)
        sage: R.<x> = LaurentSeriesRing(ZZ)
        sage: 1/(1 - 2*x)
        1 + 2*x + 4*x^2 + 8*x^3 + 16*x^4 + O(x^5)
        sage: set_series_precision(20)
    """
    if not names is None: name = names
    if name is None:
        raise TypeError("You must specify the name of the indeterminate of the Laurent series ring.")

    if default_prec is None:
        from sage.misc.defaults import series_precision
        default_prec = series_precision()

    global laurent_series
    key = (base_ring, name, default_prec, sparse)
    if key in laurent_series:
        x = laurent_series[key]()
        if x is not None: return x

    if isinstance(base_ring, ring.Field):
        R = LaurentSeriesRing_field(base_ring, name, default_prec, sparse)
    elif isinstance(base_ring, ring.IntegralDomain):
        R = LaurentSeriesRing_domain(base_ring, name, default_prec, sparse)
    elif isinstance(base_ring, ring.CommutativeRing):
        R = LaurentSeriesRing_generic(base_ring, name, default_prec, sparse)
    else:
        raise TypeError("base_ring must be a commutative ring")
    laurent_series[key] = weakref.ref(R)
    return R

def is_LaurentSeriesRing(x):
    """
    TESTS::

        sage: from sage.rings.laurent_series_ring import is_LaurentSeriesRing
        sage: K.<q> = LaurentSeriesRing(QQ)
        sage: is_LaurentSeriesRing(K)
        True
    """
    return isinstance(x, LaurentSeriesRing_generic)

class LaurentSeriesRing_generic(ring.CommutativeRing):
    """
    Univariate Laurent Series Ring

    EXAMPLES::

        sage: K = LaurentSeriesRing(CC, 'q')
        sage: K
        Laurent Series Ring in q over Complex Field with 53 bits of precision
        sage: loads(K.dumps()) == K
        True
        sage: P = QQ[['x']]
        sage: F = Frac(P)
        sage: TestSuite(F).run()

    When the base ring `k` is a field, the ring `k((x))` is a CDVF, that is
    a field equipped with a discrete valuation for which it is complete.
    The appropriate (sub)category is automatically set in this case::

        sage: k = GF(11)
        sage: R.<x> = k[[]]
        sage: F = Frac(R)
        sage: F.category()
        Category of complete discrete valuation fields
        sage: TestSuite(F).run()
    """

    def __init__(self, base_ring, name=None, default_prec=None, sparse=False, category=None):
        """
        Initialization

        EXAMPLES::

            sage: K.<q> = LaurentSeriesRing(QQ,default_prec=4); K
            Laurent Series Ring in q over Rational Field
            sage: 1 / (q-q^2)
            q^-1 + 1 + q + q^2 + O(q^3)
        """
        from .power_series_ring import PowerSeriesRing
        ring.CommutativeRing.__init__(self, base_ring, names=name,
                                                  category=getattr(self, '_default_category', Fields()))
        self._polynomial_ring = polynomial.polynomial_ring_constructor.PolynomialRing(self.base_ring(),
                                                                                      self.variable_name(),
                                                                                      sparse=sparse)
        self._power_series_ring = PowerSeriesRing(self.base_ring(),
                                                  self.variable_name(),
                                                  default_prec=default_prec,
                                                  sparse=sparse)

    def base_extend(self, R):
        """
        Returns the laurent series ring over R in the same variable as
        self, assuming there is a canonical coerce map from the base ring
        of self to R.

        EXAMPLES::

            sage: K.<x> = LaurentSeriesRing(QQ, default_prec=4)
            sage: K.base_extend(QQ['t'])
            Laurent Series Ring in x over Univariate Polynomial Ring in t over Rational Field
        """
        if R.has_coerce_map_from(self.base_ring()):
            return self.change_ring(R)
        else:
            raise TypeError("no valid base extension defined")

    def change_ring(self, R):
        """
        EXAMPLES::

            sage: K.<x> = LaurentSeriesRing(QQ, default_prec=4)
            sage: R = K.change_ring(ZZ); R
            Laurent Series Ring in x over Integer Ring
            sage: R.default_prec()
            4
        """
        return LaurentSeriesRing(R, self.variable_name(), default_prec=self.default_prec(), sparse=self.is_sparse())

    def is_sparse(self):
        """
        EXAMPLES::

            sage: K.<x> = LaurentSeriesRing(QQ, sparse=True)
            sage: K.is_sparse()
            True
        """
        return self.power_series_ring().is_sparse()

    def is_field(self, proof = True):
        """
        A Laurent series ring is a field if and only if the base ring is a field.

        TESTS::

            sage: LaurentSeriesRing(QQ,'t').is_field()
            True
            sage: LaurentSeriesRing(ZZ,'t').is_field()
            False
        """
        return self.base_ring().is_field()

    def is_dense(self):
        """
        EXAMPLES::

            sage: K.<x> = LaurentSeriesRing(QQ, sparse=True)
            sage: K.is_dense()
            False
        """
        return self.power_series_ring().is_dense()

    def __reduce__(self):
        """
        TESTS::

            sage: K.<q> = LaurentSeriesRing(QQ,default_prec=4)
            sage: L = loads(dumps(K))
            sage: L.default_prec()
            4
        """
        return self.__class__, (self.base_ring(), self.variable_name(), self.default_prec(), self.is_sparse())

    def _repr_(self):
        """
        EXAMPLES::

            sage: LaurentSeriesRing(QQ,'q') # indirect doctest
            Laurent Series Ring in q over Rational Field
            sage: LaurentSeriesRing(ZZ,'t',sparse=True)
            Sparse Laurent Series Ring in t over Integer Ring
        """
        s = "Laurent Series Ring in %s over %s"%(self.variable_name(), self.base_ring())
        if self.is_sparse():
            s = 'Sparse ' + s
        return s

    Element = laurent_series_ring_element.LaurentSeries

    def _element_constructor_(self, x, n=0):
        r"""
        Construct a Laurent series from `x`.

        INPUT:

        - ``x`` -- object that can be converted into a Laurent series

        - ``n`` -- (default: 0) multiply the result by `t^n`

        EXAMPLES::

            sage: R.<u> = LaurentSeriesRing(Qp(5, 10))
            sage: S.<t> = LaurentSeriesRing(RationalField())
            sage: R(t + t^2 + O(t^3))
            (1 + O(5^10))*u + (1 + O(5^10))*u^2 + O(u^3)

        Note that coercing an element into its own parent just produces
        that element again (since Laurent series are immutable)::

            sage: u is R(u)
            True

        Rational functions are accepted::

            sage: I = sqrt(-1)
            sage: K.<I> = QQ[I]
            sage: P.<t> = PolynomialRing(K)
            sage: L.<u> = LaurentSeriesRing(QQ[I])
            sage: L((t*I)/(t^3+I*2*t))
            1/2 + 1/4*I*u^2 - 1/8*u^4 - 1/16*I*u^6 + 1/32*u^8 +
            1/64*I*u^10 - 1/128*u^12 - 1/256*I*u^14 + 1/512*u^16 +
            1/1024*I*u^18 + O(u^20)

        ::

            sage: L(t*I) / L(t^3+I*2*t)
            1/2 + 1/4*I*u^2 - 1/8*u^4 - 1/16*I*u^6 + 1/32*u^8 +
            1/64*I*u^10 - 1/128*u^12 - 1/256*I*u^14 + 1/512*u^16 +
            1/1024*I*u^18 + O(u^20)

        TESTS:

        When converting from `R((z))` to `R((z))((w))`, the variable
        `z` is sent to `z` rather than to `w` (see :trac:`7085`)::

            sage: A.<z> = LaurentSeriesRing(QQ)
            sage: B.<w> = LaurentSeriesRing(A)
            sage: B(z)
            z
            sage: z/w
            z*w^-1

        Various conversions from PARI (see also :trac:`2508`)::

            sage: L.<q> = LaurentSeriesRing(QQ)
            sage: L.set_default_prec(10)
            doctest:...: DeprecationWarning: This method is deprecated.
            See http://trac.sagemath.org/16201 for details.
            sage: L(pari('1/x'))
            q^-1
            sage: L(pari('polchebyshev(5)'))
            5*q - 20*q^3 + 16*q^5
            sage: L(pari('polchebyshev(5) - 1/x^4'))
            -q^-4 + 5*q - 20*q^3 + 16*q^5
            sage: L(pari('1/polchebyshev(5)'))
            1/5*q^-1 + 4/5*q + 64/25*q^3 + 192/25*q^5 + 2816/125*q^7 + O(q^9)
            sage: L(pari('polchebyshev(5) + O(x^40)'))
            5*q - 20*q^3 + 16*q^5 + O(q^40)
            sage: L(pari('polchebyshev(5) - 1/x^4 + O(x^40)'))
            -q^-4 + 5*q - 20*q^3 + 16*q^5 + O(q^40)
            sage: L(pari('1/polchebyshev(5) + O(x^10)'))
            1/5*q^-1 + 4/5*q + 64/25*q^3 + 192/25*q^5 + 2816/125*q^7 + 8192/125*q^9 + O(q^10)
            sage: L(pari('1/polchebyshev(5) + O(x^10)'), -10)  # Multiply by q^-10
            1/5*q^-11 + 4/5*q^-9 + 64/25*q^-7 + 192/25*q^-5 + 2816/125*q^-3 + 8192/125*q^-1 + O(1)
            sage: L(pari('O(x^-10)'))
            O(q^-10)
        """
        from sage.rings.fraction_field_element import is_FractionFieldElement
        from sage.rings.polynomial.polynomial_element import is_Polynomial
        from sage.rings.polynomial.multi_polynomial_element import is_MPolynomial
        from sage.structure.element import parent

        P = parent(x)
        if isinstance(x, self.element_class) and n==0 and P is self:
            return x  # ok, since Laurent series are immutable (no need to make a copy)
        elif P is self.base_ring():
            # Convert x into a power series; if P is itself a Laurent
            # series ring A((t)), this prevents the implementation of
            # LaurentSeries.__init__() from effectively applying the
            # ring homomorphism A((t)) -> A((t))((u)) sending t to u
            # instead of the one sending t to t.  We cannot easily
            # tell LaurentSeries.__init__() to be more strict, because
            # A((t)) -> B((u)) is expected to send t to u if A admits
            # a coercion to B but A((t)) does not, and this condition
            # would be inefficient to check there.
            x = self.power_series_ring()(x)
        elif isinstance(x, pari_gen):
            t = x.type()
            if t == "t_RFRAC":   # Rational function
                x = self(self.polynomial_ring()(x.numerator())) / \
                    self(self.polynomial_ring()(x.denominator()))
                return (x << n)
            elif t == "t_SER":   # Laurent series
                n += x._valp()
                bigoh = n + x.length()
                x = self(self.polynomial_ring()(x.Vec()))
                return (x << n).add_bigoh(bigoh)
            else:  # General case, pretend to be a polynomial
                return self(self.polynomial_ring()(x)) << n
        elif is_FractionFieldElement(x) and \
             (x.base_ring() is self.base_ring() or x.base_ring() == self.base_ring()) and \
             (is_Polynomial(x.numerator()) or is_MPolynomial(x.numerator())):
            x = self(x.numerator())/self(x.denominator())
            return (x << n)
        return self.element_class(self, x, n)

    def _coerce_map_from_(self, P):
        """
        Return a coercion map from `P` to ``self``, or True, or None.

        The following rings admit a coercion map to the Laurent series
        ring `A((t))`:

        - any ring that admits a coercion map to `A` (including `A`
          itself);

        - any Laurent series ring, power series ring or polynomial
          ring in the variable `t` over a ring admitting a coercion
          map to `A`.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(ZZ)
            sage: S.<t> = PowerSeriesRing(QQ)
            sage: R.has_coerce_map_from(S) # indirect doctest
            False
            sage: R.has_coerce_map_from(R)
            True
            sage: R.<t> = LaurentSeriesRing(QQ['x'])
            sage: R.has_coerce_map_from(S)
            True
            sage: R.has_coerce_map_from(QQ['t'])
            True
            sage: R.has_coerce_map_from(ZZ['x']['t'])
            True
            sage: R.has_coerce_map_from(ZZ['t']['x'])
            False
            sage: R.has_coerce_map_from(ZZ['x'])
            True
        """
        A = self.base_ring()
        if A is P:
            return True
        f = A.coerce_map_from(P)
        if f is not None:
            return self.coerce_map_from(A) * f

        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        from sage.rings.power_series_ring import is_PowerSeriesRing
        if ((is_LaurentSeriesRing(P) or is_PowerSeriesRing(P) or is_PolynomialRing(P))
            and P.variable_name() == self.variable_name()
            and A.has_coerce_map_from(P.base_ring())):
            return True

    def __cmp__(self, other):
        """
        Compare this Laurent series ring to something else.

        Laurent series rings are considered equal if the base ring, variable
        names, and default truncation precision are the same.

        First the base rings are compared, then the variable names, then
        the default precision.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(ZZ)
            sage: S.<t> = LaurentSeriesRing(ZZ)
            sage: R is S
            True
            sage: R == S
            True
            sage: S.<t> = LaurentSeriesRing(ZZ, default_prec=10)
            sage: R == S
            False
            sage: LaurentSeriesRing(ZZ,'t') == LaurentSeriesRing(QQ,'t')
            False
            sage: LaurentSeriesRing(QQ,'t') == 5
            False
        """
        if not isinstance(other, LaurentSeriesRing_generic):
            return cmp(type(self),type(other))
        c = cmp(self.base_ring(), other.base_ring())
        if c: return c
        c = cmp(self.variable_name(), other.variable_name())
        if c: return c
        c = cmp(self.default_prec(), other.default_prec())
        if c: return c
        return 0


    def _is_valid_homomorphism_(self, codomain, im_gens):
        """
        EXAMPLES::

            sage: R.<x> = LaurentSeriesRing(GF(17))
            sage: S.<y> = LaurentSeriesRing(GF(19))
            sage: R.hom([y], S) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: images do not define a valid homomorphism
            sage: f = R.hom(x+x^3,R)
            sage: f(x^2)
            x^2 + 2*x^4 + x^6
        """
        ## NOTE: There are no ring homomorphisms from the ring of
        ## all formal power series to most rings, e.g, the p-adic
        ## field, since you can always (mathematically!) construct
        ## some power series that doesn't converge.
        ## Note that 0 is not a *ring* homomorphism.
        from .power_series_ring import is_PowerSeriesRing
        if is_PowerSeriesRing(codomain) or is_LaurentSeriesRing(codomain):
            return im_gens[0].valuation() > 0 and codomain.has_coerce_map_from(self.base_ring())
        return False

    def characteristic(self):
        """
        EXAMPLES::

            sage: R.<x> = LaurentSeriesRing(GF(17))
            sage: R.characteristic()
            17
        """
        return self.base_ring().characteristic()

    def residue_field(self):
        """
        Return the residue field of this Laurent series field
        if it is a complete discrete valuation field (i.e. if
        the base ring is a field, in which base it is also the
        residue field).

        EXAMPLES::

            sage: R.<x> = LaurentSeriesRing(GF(17))
            sage: R.residue_field()
            Finite Field of size 17

            sage: R.<x> = LaurentSeriesRing(ZZ)
            sage: R.residue_field()
            Traceback (most recent call last):
            ...
            TypeError: The base ring is not a field
        """
        if self.base_ring().is_field():
            return self.base_ring()
        else:
            raise TypeError("The base ring is not a field")

    def set_default_prec(self, n):
        """
        Set the default precision.

        This method is deprecated.

        TESTS::

            sage: R.<x> = LaurentSeriesRing(QQ)
            sage: R.set_default_prec(3)
            doctest:...: DeprecationWarning: This method is deprecated.
            See http://trac.sagemath.org/16201 for details.
            sage: 1/(x^5-x^7)
            x^-5 + x^-3 + O(x^-2)
        """
        from sage.misc.superseded import deprecation
        deprecation(16201, "This method is deprecated.")
        self.power_series_ring().set_default_prec(n)

    def default_prec(self):
        """
        Get the precision to which exact elements are truncated when
        necessary (most frequently when inverting).

        EXAMPLES::

            sage: R.<x> = LaurentSeriesRing(QQ, default_prec=5)
            sage: R.default_prec()
            5
        """
        return self.power_series_ring().default_prec()

    def is_exact(self):
        """
        Laurent series rings are inexact.

        EXAMPLES::

            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.is_exact()
            False
        """
        return False

    def gen(self, n=0):
        """
        EXAMPLES::

            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.gen()
            x
        """
        if n != 0:
            raise IndexError("Generator n not defined.")
        try:
            return self.__generator
        except AttributeError:
            self.__generator = laurent_series_ring_element.LaurentSeries(self, [0,1])
            return self.__generator

    def uniformizer(self):
        """
        Return a uniformizer of this Laurent series field if it is
        a discrete valuation field (i.e. if the base ring is actually
        a field). Otherwise, an error is raised.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: R.uniformizer()
            t

            sage: R.<t> = LaurentSeriesRing(ZZ)
            sage: R.uniformizer()
            Traceback (most recent call last):
            ...
            TypeError: The base ring is not a field
        """
        if self.base_ring().is_field():
            return self.gen()
        else:
            raise TypeError("The base ring is not a field")

    def ngens(self):
        """
        Laurent series rings are univariate.

        EXAMPLES::

            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.ngens()
            1
        """
        return 1

    def polynomial_ring(self):
        r"""
        If this is the Laurent series ring `R((t))`, return the
        polynomial ring `R[t]`.

        EXAMPLES::

            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.polynomial_ring()
            Univariate Polynomial Ring in x over Rational Field
        """
        return self._polynomial_ring

    def laurent_polynomial_ring(self):
        r"""
        If this is the Laurent series ring `R((t))`, return the Laurent
        polynomial ring `R[t,1/t]`.

        EXAMPLES::

            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.laurent_polynomial_ring()
            Univariate Laurent Polynomial Ring in x over Rational Field
        """
        try:
            return self.__laurent_polynomial_ring
        except AttributeError:
            self.__laurent_polynomial_ring = polynomial.laurent_polynomial_ring.LaurentPolynomialRing( \
                                         self.base_ring(), self.variable_name(), sparse=self.is_sparse())
            return self.__laurent_polynomial_ring

    def power_series_ring(self):
        r"""
        If this is the Laurent series ring `R((t))`, return the
        power series ring `R[[t]]`.

        EXAMPLES::

            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.power_series_ring()
            Power Series Ring in x over Rational Field
        """
        return self._power_series_ring

class LaurentSeriesRing_domain(LaurentSeriesRing_generic, ring.IntegralDomain):
    def __init__(self, base_ring, name=None, default_prec=None, sparse=False):
        """
        Initialization

        TESTS::

            sage: TestSuite(LaurentSeriesRing(ZZ,'t')).run()
        """
        LaurentSeriesRing_generic.__init__(self, base_ring, name, default_prec, sparse)

class LaurentSeriesRing_field(LaurentSeriesRing_generic, ring.Field):
    _default_category = CompleteDiscreteValuationFields()

    def __init__(self, base_ring, name=None, default_prec=None, sparse=False):
        """
        Initialization

        TESTS::

            sage: TestSuite(LaurentSeriesRing(QQ,'t')).run()
        """
        LaurentSeriesRing_generic.__init__(self, base_ring, name, default_prec, sparse)
