r"""
Power Series Rings

Power series rings are constructed in the standard Sage
fashion.  See also :doc:`multi_power_series_ring`.

EXAMPLES:

Construct rings and elements::

    sage: R.<t> = PowerSeriesRing(QQ)
    sage: R.random_element(6)  # random
    -4 - 1/2*t^2 - 1/95*t^3 + 1/2*t^4 - 12*t^5 + O(t^6)

::

    sage: R.<t,u,v> = PowerSeriesRing(QQ); R
    Multivariate Power Series Ring in t, u, v over Rational Field
    sage: p = -t + 1/2*t^3*u - 1/4*t^4*u + 2/3*v^5 + R.O(6); p
    -t + 1/2*t^3*u - 1/4*t^4*u + 2/3*v^5 + O(t, u, v)^6
    sage: p in R
    True

The default precision is specified at construction, but does not
bound the precision of created elements.

::

    sage: R.<t> = PowerSeriesRing(QQ, default_prec=5)
    sage: R.random_element(6)  # random
    1/2 - 1/4*t + 2/3*t^2 - 5/2*t^3 + 2/3*t^5 + O(t^6)

Construct univariate power series from a list of coefficients::

    sage: S = R([1, 3, 5, 7]); S
    1 + 3*t + 5*t^2 + 7*t^3

The default precision of a power series ring stays fixed and cannot be
changed. To work with different default precision, create a new power series
ring::

    sage: R.<x> = PowerSeriesRing(QQ, default_prec=10)
    sage: sin(x)
    x - 1/6*x^3 + 1/120*x^5 - 1/5040*x^7 + 1/362880*x^9 + O(x^10)
    sage: R.<x> = PowerSeriesRing(QQ, default_prec=15)
    sage: sin(x)
    x - 1/6*x^3 + 1/120*x^5 - 1/5040*x^7 + 1/362880*x^9 - 1/39916800*x^11 + 1/6227020800*x^13 + O(x^15)

An iterated example::

    sage: R.<t> = PowerSeriesRing(ZZ)
    sage: S.<t2> = PowerSeriesRing(R)
    sage: S
    Power Series Ring in t2 over Power Series Ring in t over Integer Ring
    sage: S.base_ring()
    Power Series Ring in t over Integer Ring

Sage can compute with power series over the symbolic ring.

::

    sage: K.<t> = PowerSeriesRing(SR, default_prec=5)
    sage: a, b, c = var('a,b,c')
    sage: f = a + b*t + c*t^2 + O(t^3)
    sage: f*f
    a^2 + 2*a*b*t + (b^2 + 2*a*c)*t^2 + O(t^3)
    sage: f = sqrt(2) + sqrt(3)*t + O(t^3)
    sage: f^2
    2 + 2*sqrt(3)*sqrt(2)*t + 3*t^2 + O(t^3)

Elements are first coerced to constants in ``base_ring``, then coerced
into the ``PowerSeriesRing``::

    sage: R.<t> = PowerSeriesRing(ZZ)
    sage: f = Mod(2, 3) * t; (f, f.parent())
    (2*t, Power Series Ring in t over Ring of integers modulo 3)

We make a sparse power series.

::

    sage: R.<x> = PowerSeriesRing(QQ, sparse=True); R
    Sparse Power Series Ring in x over Rational Field
    sage: f = 1 + x^1000000
    sage: g = f*f
    sage: g.degree()
    2000000

We make a sparse Laurent series from a power series generator::

    sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
    sage: latex(-2/3*(1/t^3) + 1/t + 3/5*t^2 + O(t^5))
    \frac{-\frac{2}{3}}{t^{3}} + \frac{1}{t} + \frac{3}{5}t^{2} + O(t^{5})
    sage: S = parent(1/t); S
    Sparse Laurent Series Ring in t over Rational Field

Choose another implementation of the attached polynomial ring::

    sage: R.<t> = PowerSeriesRing(ZZ)
    sage: type(t.polynomial())
    <... 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>
    sage: S.<s> = PowerSeriesRing(ZZ, implementation='NTL')
    sage: type(s.polynomial())
    <... 'sage.rings.polynomial.polynomial_integer_dense_ntl.Polynomial_integer_dense_ntl'>

AUTHORS:

- William Stein: the code
- Jeremy Cho (2006-05-17): some examples (above)
- Niles Johnson (2010-09): implement multivariate power series
- Simon King (2012-08): use category and coercion framework, :trac:`13412`

TESTS::

    sage: R.<t> = PowerSeriesRing(QQ)
    sage: R is loads(dumps(R))
    True
    sage: TestSuite(R).run()

::

    sage: R.<x> = PowerSeriesRing(QQ, sparse=True)
    sage: R is loads(dumps(R))
    True
    sage: TestSuite(R).run()

::

    sage: M = PowerSeriesRing(QQ, 't,u,v,w', default_prec=20)
    sage: M is loads(dumps(M))
    True
    sage: TestSuite(M).run()

"""

from . import power_series_poly
from . import power_series_mpoly
from .power_series_pari import PowerSeries_pari
from . import power_series_ring_element

from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.rings.polynomial.multi_polynomial_ring_base import is_MPolynomialRing
from .polynomial.polynomial_ring_constructor import PolynomialRing
from . import laurent_series_ring
from . import laurent_series_ring_element
from . import integer
from . import ring
from .infinity import infinity
import sage.misc.latex as latex
from sage.structure.nonexact import Nonexact

from sage.interfaces.magma import MagmaElement
from sage.rings.fraction_field_element import FractionFieldElement
from sage.misc.sage_eval import sage_eval

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.category_object import normalize_names
from sage.structure.element import parent
import sage.categories.commutative_rings as commutative_rings
_CommutativeRings = commutative_rings.CommutativeRings()
import sage.categories.integral_domains as integral_domains
_IntegralDomains = integral_domains.IntegralDomains()
import sage.categories.fields as fields
_Fields = fields.Fields()

from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationRings


def PowerSeriesRing(base_ring, name=None, arg2=None, names=None,
                    sparse=False, default_prec=None, order='negdeglex',
                    num_gens=None, implementation=None):
    r"""
    Create a univariate or multivariate power series ring over a given
    (commutative) base ring.

    INPUT:


    -  ``base_ring`` - a commutative ring

    -  ``name``, ``names`` - name(s) of the indeterminate

    - ``default_prec`` - the default precision used if an exact object must
       be changed to an approximate object in order to do an arithmetic
       operation.  If left as ``None``, it will be set to the global
       default (20) in the univariate case, and 12 in the multivariate case.

    -  ``sparse`` - (default: ``False``) whether power series
       are represented as sparse objects.

    - ``order`` - (default: ``negdeglex``) term ordering, for multivariate case

    - ``num_gens`` - number of generators, for multivariate case


    There is a unique power series ring over each base ring with given
    variable name. Two power series over the same base ring with
    different variable names are not equal or isomorphic.

    EXAMPLES (Univariate)::

        sage: R = PowerSeriesRing(QQ, 'x'); R
        Power Series Ring in x over Rational Field

    ::

        sage: S = PowerSeriesRing(QQ, 'y'); S
        Power Series Ring in y over Rational Field

    ::

        sage: R = PowerSeriesRing(QQ, 10)
        Traceback (most recent call last):
        ...
        ValueError: variable name '10' does not start with a letter

    ::

        sage: S = PowerSeriesRing(QQ, 'x', default_prec = 15); S
        Power Series Ring in x over Rational Field
        sage: S.default_prec()
        15

    EXAMPLES (Multivariate) See also :doc:`multi_power_series_ring`::

        sage: R = PowerSeriesRing(QQ, 't,u,v'); R
        Multivariate Power Series Ring in t, u, v over Rational Field

    ::

        sage: N = PowerSeriesRing(QQ,'w',num_gens=5); N
        Multivariate Power Series Ring in w0, w1, w2, w3, w4 over Rational Field

    Number of generators can be specified before variable name without using keyword::

        sage: M = PowerSeriesRing(QQ,4,'k'); M
        Multivariate Power Series Ring in k0, k1, k2, k3 over Rational Field

    Multivariate power series can be constructed using angle bracket or double square bracket notation::

        sage: R.<t,u,v> = PowerSeriesRing(QQ, 't,u,v'); R
        Multivariate Power Series Ring in t, u, v over Rational Field

        sage: ZZ[['s,t,u']]
        Multivariate Power Series Ring in s, t, u over Integer Ring

    Sparse multivariate power series ring::

        sage: M = PowerSeriesRing(QQ,4,'k',sparse=True); M
        Sparse Multivariate Power Series Ring in k0, k1, k2, k3 over
        Rational Field

    Power series ring over polynomial ring::

        sage: H = PowerSeriesRing(PolynomialRing(ZZ,3,'z'),4,'f'); H
        Multivariate Power Series Ring in f0, f1, f2, f3 over Multivariate
        Polynomial Ring in z0, z1, z2 over Integer Ring

    Power series ring over finite field::

        sage: S = PowerSeriesRing(GF(65537),'x,y'); S
        Multivariate Power Series Ring in x, y over Finite Field of size
        65537

    Power series ring with many variables::

        sage: R = PowerSeriesRing(ZZ, ['x%s'%p for p in primes(100)]); R
        Multivariate Power Series Ring in x2, x3, x5, x7, x11, x13, x17, x19,
        x23, x29, x31, x37, x41, x43, x47, x53, x59, x61, x67, x71, x73, x79,
        x83, x89, x97 over Integer Ring

    - Use :meth:`inject_variables` to make the variables available for
      interactive use.

      ::

        sage: R.inject_variables()
        Defining x2, x3, x5, x7, x11, x13, x17, x19, x23, x29, x31, x37,
        x41, x43, x47, x53, x59, x61, x67, x71, x73, x79, x83, x89, x97

        sage: f = x47 + 3*x11*x29 - x19 + R.O(3)
        sage: f in R
        True


    Variable ordering determines how series are displayed::

        sage: T.<a,b> = PowerSeriesRing(ZZ,order='deglex'); T
        Multivariate Power Series Ring in a, b over Integer Ring
        sage: T.term_order()
        Degree lexicographic term order
        sage: p = - 2*b^6 + a^5*b^2 + a^7 - b^2 - a*b^3 + T.O(9); p
        a^7 + a^5*b^2 - 2*b^6 - a*b^3 - b^2 + O(a, b)^9

        sage: U = PowerSeriesRing(ZZ,'a,b',order='negdeglex'); U
        Multivariate Power Series Ring in a, b over Integer Ring
        sage: U.term_order()
        Negative degree lexicographic term order
        sage: U(p)
        -b^2 - a*b^3 - 2*b^6 + a^7 + a^5*b^2 + O(a, b)^9


    TESTS::

        sage: N = PowerSeriesRing(QQ,'k',num_gens=5); N
        Multivariate Power Series Ring in k0, k1, k2, k3, k4 over Rational Field

    The following behavior of univariate power series ring will eventually
    be deprecated and then changed to return a multivariate power series
    ring::

        sage: N = PowerSeriesRing(QQ,'k',5); N
        Power Series Ring in k over Rational Field
        sage: N.default_prec()
        5
        sage: L.<m> = PowerSeriesRing(QQ,5); L
        Power Series Ring in m over Rational Field
        sage: L.default_prec()
        5

    By :trac:`14084`, a power series ring belongs to the category of integral
    domains, if the base ring does::

        sage: P = ZZ[['x']]
        sage: P.category()
        Category of integral domains
        sage: TestSuite(P).run()
        sage: M = ZZ[['x','y']]
        sage: M.category()
        Category of integral domains
        sage: TestSuite(M).run()

    Otherwise, it belongs to the category of commutative rings::

        sage: P = Integers(15)[['x']]
        sage: P.category()
        Category of commutative rings
        sage: TestSuite(P).run()
        sage: M = Integers(15)[['x','y']]
        sage: M.category()
        Category of commutative rings
        sage: TestSuite(M).run()

    .. SEEALSO::

        * :func:`sage.misc.defaults.set_series_precision`
    """
    #multivariate case:
    # examples for first case:
    # PowerSeriesRing(QQ,'x,y,z')
    # PowerSeriesRing(QQ,['x','y','z'])
    # PowerSeriesRing(QQ,['x','y','z'], 3)
    if names is None and name is not None:
        names = name
    if isinstance(names, (tuple, list)) and len(names) > 1 or (isinstance(names, str) and ',' in names):
        return _multi_variate(base_ring, num_gens=arg2, names=names,
                     order=order, default_prec=default_prec, sparse=sparse)
    # examples for second case:
    # PowerSeriesRing(QQ,3,'t')
    if arg2 is None and num_gens is not None:
        arg2 = names
        names = num_gens
    if (isinstance(arg2, str) and
            isinstance(names, (int, integer.Integer))):
        return _multi_variate(base_ring, num_gens=names, names=arg2,
                     order=order, default_prec=default_prec, sparse=sparse)


    # univariate case: the arguments to PowerSeriesRing used to be
    # (base_ring, name=None, default_prec=20, names=None, sparse=False),
    # and thus that is what the code below expects; this behavior is being
    # deprecated, and will eventually be removed.
    if default_prec is None and arg2 is None:
        from sage.misc.defaults import series_precision
        default_prec = series_precision()
    elif arg2 is not None:
        default_prec = arg2

    ## too many things (padics, elliptic curves) depend on this behavior,
    ## so no warning for now.
    ##

    # if isinstance(name, (int,integer.Integer)) or isinstance(arg2,(int,integer.Integer)):
    #     deprecation(trac_number, "This behavior of PowerSeriesRing is being deprecated in favor of constructing multivariate power series rings. (See Trac ticket #1956.)")

    # the following is the original, univariate-only code

    if isinstance(name, (int, integer.Integer)):
        default_prec = name
    if names is not None:
        name = names
    name = normalize_names(1, name)

    if name is None:
        raise TypeError("You must specify the name of the indeterminate of the Power series ring.")

    key = (base_ring, name, default_prec, sparse, implementation)
    if PowerSeriesRing_generic.__classcall__.is_in_cache(key):
        return PowerSeriesRing_generic(*key)

    if isinstance(name, (tuple, list)):
        assert len(name) == 1
        name = name[0]

    if not (name is None or isinstance(name, str)):
        raise TypeError("variable name must be a string or None")

    if base_ring in _Fields:
        R = PowerSeriesRing_over_field(base_ring, name, default_prec,
                                       sparse=sparse, implementation=implementation)
    elif base_ring in _IntegralDomains:
        R = PowerSeriesRing_domain(base_ring, name, default_prec,
                                   sparse=sparse, implementation=implementation)
    elif base_ring in _CommutativeRings:
        R = PowerSeriesRing_generic(base_ring, name, default_prec,
                                    sparse=sparse, implementation=implementation)
    else:
        raise TypeError("base_ring must be a commutative ring")
    return R

def _multi_variate(base_ring, num_gens=None, names=None,
                     order='negdeglex', default_prec=None, sparse=False):
    """
    Construct multivariate power series ring.
    """
    if names is None:
        raise TypeError("you must specify a variable name or names")

    if num_gens is None:
        if isinstance(names,str):
            num_gens = len(names.split(','))
        elif isinstance(names, (list, tuple)):
            num_gens = len(names)
        else:
            raise TypeError("variable names must be a string, tuple or list")
    names = normalize_names(num_gens, names)
    num_gens = len(names)
    if default_prec is None:
        default_prec = 12

    if base_ring not in commutative_rings.CommutativeRings():
        raise TypeError("base_ring must be a commutative ring")
    from sage.rings.multi_power_series_ring import MPowerSeriesRing_generic
    R = MPowerSeriesRing_generic(base_ring, num_gens, names,
                                 order=order, default_prec=default_prec, sparse=sparse)
    return R


def _single_variate():
    pass

def is_PowerSeriesRing(R):
    """
    Return True if this is a *univariate* power series ring.  This is in
    keeping with the behavior of ``is_PolynomialRing``
    versus ``is_MPolynomialRing``.

    EXAMPLES::

        sage: from sage.rings.power_series_ring import is_PowerSeriesRing
        sage: is_PowerSeriesRing(10)
        False
        sage: is_PowerSeriesRing(QQ[['x']])
        True
    """
    if isinstance(R, PowerSeriesRing_generic):
        return R.ngens() == 1
    else:
        return False

class PowerSeriesRing_generic(UniqueRepresentation, ring.CommutativeRing, Nonexact):
    """
    A power series ring.
    """

    def __init__(self, base_ring, name=None, default_prec=None, sparse=False,
                 implementation=None, category=None):
        """
        Initializes a power series ring.

        INPUT:


        -  ``base_ring`` - a commutative ring

        -  ``name`` - name of the indeterminate

        -  ``default_prec`` - the default precision

        -  ``sparse`` - whether or not power series are
           sparse

        - ``implementation`` -- either ``'poly'``, ``'mpoly'``, or
          ``'pari'``. Other values (for example ``'NTL'``) are passed to
          the attached polynomial ring.  The default is ``'pari'`` if
          the base field is a PARI finite field, and ``'poly'`` otherwise.

        If the base ring is a polynomial ring, then the option
        ``implementation='mpoly'`` causes computations to be done with
        multivariate polynomials instead of a univariate polynomial
        ring over the base ring.  Only use this for dense power series
        where you won't do too much arithmetic, but the arithmetic you
        do must be fast.  You must explicitly call
        ``f.do_truncation()`` on an element for it to truncate away
        higher order terms (this is called automatically before
        printing).

        EXAMPLES:

        This base class inherits from :class:`~sage.rings.ring.CommutativeRing`.
        Since :trac:`11900`, it is also initialised as such, and since :trac:`14084`
        it is actually initialised as an integral domain::

            sage: R.<x> = ZZ[[]]
            sage: R.category()
            Category of integral domains
            sage: TestSuite(R).run()

        When the base ring `k` is a field, the ring `k[[x]]` is not only a
        commutative ring, but also a complete discrete valuation ring (CDVR).
        The appropriate (sub)category is automatically set in this case::

            sage: k = GF(11)
            sage: R.<x> = k[[]]
            sage: R.category()
            Category of complete discrete valuation rings
            sage: TestSuite(R).run()

        It is checked that the default precision is non-negative
        (see :trac:`19409`)::

            sage: PowerSeriesRing(ZZ, 'x', default_prec=-5)
            Traceback (most recent call last):
            ...
            ValueError: default_prec (= -5) must be non-negative

        """
        from sage.rings.finite_rings.finite_field_pari_ffelt import FiniteField_pari_ffelt

        if implementation is None:
            if isinstance(base_ring, FiniteField_pari_ffelt):
                implementation = 'pari'
            else:
                implementation = 'poly'
            R = PolynomialRing(base_ring, name, sparse=sparse)
        elif implementation not in ['pari', 'mpoly']:     # see :trac:`28996`
            R = PolynomialRing(base_ring, name, sparse=sparse, implementation=implementation)
            implementation = 'poly'
        else:
            R = PolynomialRing(base_ring, name, sparse=sparse)

        self.__poly_ring = R
        self.__is_sparse = sparse
        if default_prec is None:
            from sage.misc.defaults import series_precision
            default_prec = series_precision()
        elif default_prec < 0:
            raise ValueError("default_prec (= %s) must be non-negative"
                             % default_prec)

        if implementation == 'poly':
            self.Element = power_series_poly.PowerSeries_poly
        elif implementation == 'mpoly':
            K = base_ring
            names = K.variable_names() + (name,)
            self.__mpoly_ring = PolynomialRing(K.base_ring(), names=names)
            assert is_MPolynomialRing(self.__mpoly_ring)
            self.Element = power_series_mpoly.PowerSeries_mpoly
        elif implementation == 'pari':
            self.Element = PowerSeries_pari
        else:
            raise ValueError('unknown power series implementation: %r' % implementation)

        ring.CommutativeRing.__init__(self, base_ring, names=name,
                                      category=getattr(self, '_default_category',
                                                       _CommutativeRings))
        Nonexact.__init__(self, default_prec)
        if self.Element is PowerSeries_pari:
            self.__generator = self.element_class(self, R.gen().__pari__())
        else:
            self.__generator = self.element_class(self, R.gen(), is_gen=True)

    def variable_names_recursive(self, depth=None):
        r"""
        Return the list of variable names of this and its base rings.

        EXAMPLES::

            sage: R = QQ[['x']][['y']][['z']]
            sage: R.variable_names_recursive()
            ('x', 'y', 'z')
            sage: R.variable_names_recursive(2)
            ('y', 'z')
        """
        if depth is None:
            from sage.rings.infinity import infinity
            depth = infinity

        if depth <= 0:
            all = ()
        elif depth == 1:
            all = self.variable_names()
        else:
            my_vars = self.variable_names()
            try:
               all = self.base_ring().variable_names_recursive(depth - len(my_vars)) + my_vars
            except AttributeError:
                all = my_vars
        if len(all) > depth:
            all = all[-depth:]
        return all

    def _repr_(self):
        """
        Print out a power series ring.

        EXAMPLES::

            sage: R = GF(17)[['y']]
            sage: R
            Power Series Ring in y over Finite Field of size 17
            sage: R.__repr__()
            'Power Series Ring in y over Finite Field of size 17'
            sage: R.rename('my power series ring')
            sage: R
            my power series ring
        """
        s = "Power Series Ring in %s over %s"%(self.variable_name(), self.base_ring())
        if self.is_sparse():
            s = 'Sparse ' + s
        return s

    def is_sparse(self):
        """
        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: t.is_sparse()
            False
            sage: R.<t> = PowerSeriesRing(ZZ, sparse=True)
            sage: t.is_sparse()
            True
        """
        return self.__is_sparse

    def is_dense(self):
        """
        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: t.is_dense()
            True
            sage: R.<t> = PowerSeriesRing(ZZ, sparse=True)
            sage: t.is_dense()
            False
        """
        return not self.__is_sparse

    def _latex_(self):
        r"""
        Display latex representation of this power series ring.

        EXAMPLES::

            sage: R = GF(17)[['y']]
            sage: latex(R)  # indirect doctest
            \Bold{F}_{17}[[y]]
            sage: R = GF(17)[['y12']]
            sage: latex(R)
            \Bold{F}_{17}[[y_{12}]]
        """
        return "%s[[%s]]"%(latex.latex(self.base_ring()), self.latex_variable_names()[0])

    def _coerce_map_from_(self, S):
        """
        A coercion from `S` exists, if `S` coerces into ``self``'s base ring,
        or if `S` is a univariate polynomial or power series ring with the
        same variable name as self, defined over a base ring that coerces into
        ``self``'s base ring.

        EXAMPLES::

            sage: A = GF(17)[['x']]
            sage: A.has_coerce_map_from(ZZ)  # indirect doctest
            True
            sage: A.has_coerce_map_from(ZZ['x'])
            True
            sage: A.has_coerce_map_from(ZZ['y'])
            False
            sage: A.has_coerce_map_from(ZZ[['x']])
            True

        """
        if self.base_ring().has_coerce_map_from(S):
            return True
        if (is_PolynomialRing(S) or is_PowerSeriesRing(S)) and self.base_ring().has_coerce_map_from(S.base_ring()) \
           and self.variable_names()==S.variable_names():
            return True

    def _element_constructor_(self, f, prec=infinity, check=True):
        """
        Coerce object to this power series ring.

        Returns a new instance unless the parent of f is self, in which
        case f is returned (since f is immutable).

        INPUT:


        -  ``f`` - object, e.g., a power series ring element

        -  ``prec`` - (default: infinity); truncation precision
           for coercion

        -  ``check`` - bool (default: True), whether to verify
           that the coefficients, etc., coerce in correctly.


        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R(t+O(t^5))    # indirect doctest
            t + O(t^5)
            sage: R(13)
            13
            sage: R(2/3)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: R([1,2,3])
            1 + 2*t + 3*t^2
            sage: S.<w> = PowerSeriesRing(QQ)
            sage: R(w + 3*w^2 + O(w^3))
            t + 3*t^2 + O(t^3)
            sage: x = polygen(QQ,'x')
            sage: R(x + x^2 + x^3 + x^5, 3)
            t + t^2 + O(t^3)
            sage: R(1/(1-x), prec=5)
            1 + t + t^2 + t^3 + t^4 + O(t^5)
            sage: R(1/x, 5)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Laurent Series Ring in t over
             Rational Field to Power Series Ring in t over Integer Ring

            sage: PowerSeriesRing(PowerSeriesRing(QQ,'x'),'y')(x)
            x
            sage: PowerSeriesRing(PowerSeriesRing(QQ,'y'),'x')(x)
            x
            sage: PowerSeriesRing(PowerSeriesRing(QQ,'t'),'y')(x)
            y
            sage: PowerSeriesRing(PowerSeriesRing(QQ,'t'),'y')(1/(1+x), 5)
            1 - y + y^2 - y^3 + y^4 + O(y^5)
            sage: PowerSeriesRing(PowerSeriesRing(QQ,'x',5),'y')(1/(1+x))
            1 - x + x^2 - x^3 + x^4 + O(x^5)
            sage: PowerSeriesRing(PowerSeriesRing(QQ,'y'),'x')(1/(1+x), 5)
            1 - x + x^2 - x^3 + x^4 + O(x^5)
            sage: PowerSeriesRing(PowerSeriesRing(QQ,'x'),'x')(x).coefficients()
            [x]

        Conversion from symbolic series::

            sage: x,y = var('x,y')
            sage: s=(1/(1-x)).series(x,3); s
            1 + 1*x + 1*x^2 + Order(x^3)
            sage: R.<x> = PowerSeriesRing(QQ)
            sage: R(s)
            1 + x + x^2 + O(x^3)
            sage: ex=(gamma(1-y)).series(y,3)
            sage: R.<y> = PowerSeriesRing(SR)
            sage: R(ex)
            1 + euler_gamma*y + (1/2*euler_gamma^2 + 1/12*pi^2)*y^2 + O(y^3)

        Laurent series with non-negative valuation are accepted (see
        :trac:`6431`)::

            sage: L.<q> = LaurentSeriesRing(QQ)
            sage: P = L.power_series_ring()
            sage: P(q)
            q
            sage: P(1/q)
            Traceback (most recent call last):
            ...
            TypeError: self is not a power series

        It is checked that the precision is non-negative
        (see :trac:`19409`)::

            sage: PowerSeriesRing(ZZ, 'x')(1, prec=-5)
            Traceback (most recent call last):
            ...
            ValueError: prec (= -5) must be non-negative

        """
        if prec is not infinity:
            prec = integer.Integer(prec)
            if prec < 0:
                raise ValueError("prec (= %s) must be non-negative" % prec)
        from sage.symbolic.expression import SymbolicSeries
        if isinstance(f, power_series_ring_element.PowerSeries) and f.parent() is self:
            if prec >= f.prec():
                return f
            f = f.truncate(prec)
        elif isinstance(f, laurent_series_ring_element.LaurentSeries) and f.parent().power_series_ring() is self:
            return self(f.power_series(), prec, check=check)
        elif isinstance(f, MagmaElement) and str(f.Type()) == 'RngSerPowElt':
            v = sage_eval(f.Eltseq())
            return self(v) * (self.gen(0)**f.Valuation())
        elif isinstance(f, FractionFieldElement):
            if self.base_ring().has_coerce_map_from(f.parent()):
                return self.element_class(self, [f], prec, check=check)
            else:
                num = self.element_class(self, f.numerator(), prec, check=check)
                den = self.element_class(self, f.denominator(), prec, check=check)
                return self.coerce(num/den)
        elif isinstance(f, SymbolicSeries):
            if str(f.default_variable()) == self.variable_name():
                return self.element_class(self, f.list(),
                                      f.degree(f.default_variable()), check=check)
            else:
                raise TypeError("Can only convert series into ring with same variable name.")
        return self.element_class(self, f, prec, check=check)

    def construction(self):
        """
        Return the functorial construction of self, namely, completion of
        the univariate polynomial ring with respect to the indeterminate
        (to a given precision).

        EXAMPLES::

            sage: R = PowerSeriesRing(ZZ, 'x')
            sage: c, S = R.construction(); S
            Univariate Polynomial Ring in x over Integer Ring
            sage: R == c(S)
            True
            sage: R = PowerSeriesRing(ZZ, 'x', sparse=True)
            sage: c, S = R.construction()
            sage: R == c(S)
            True

        """
        from sage.categories.pushout import CompletionFunctor
        if self.is_sparse():
            extras = {'sparse': True}
        else:
            extras = None
        return CompletionFunctor(self._names[0], self.default_prec(), extras), self._poly_ring()

    def _coerce_impl(self, x):
        """
        Return canonical coercion of x into self.

        Rings that canonically coerce to this power series ring R:

        - R itself

        - Any power series ring in the same variable whose base ring
          canonically coerces to the base ring of R.

        - Any ring that canonically coerces to the polynomial ring
          over the base ring of R.

        - Any ring that canonically coerces to the base ring of R

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R._coerce_(t + t^2)  # indirect doctest
            t + t^2
            sage: R._coerce_(1/t)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Laurent Series Ring in t over
             Rational Field to Power Series Ring in t over Integer Ring
            sage: R._coerce_(5)
            5
            sage: tt = PolynomialRing(ZZ,'t').gen()
            sage: R._coerce_(tt^2 + tt - 1)
            -1 + t + t^2
            sage: R._coerce_(1/2)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Power Series Ring in t over Integer Ring
            sage: S.<s> = PowerSeriesRing(ZZ)
            sage: R._coerce_(s)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Power Series Ring in s over Integer Ring to Power Series Ring in t over Integer Ring

        We illustrate canonical coercion between power series rings with
        compatible base rings::

            sage: R.<t> = PowerSeriesRing(GF(7)['w'])
            sage: S = PowerSeriesRing(ZZ, 't')
            sage: f = S([1,2,3,4]); f
            1 + 2*t + 3*t^2 + 4*t^3
            sage: g = R._coerce_(f); g
            1 + 2*t + 3*t^2 + 4*t^3
            sage: parent(g)
            Power Series Ring in t over Univariate Polynomial Ring in w over Finite Field of size 7
            sage: S._coerce_(g)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Power Series Ring in t over Univariate Polynomial Ring in w over Finite Field of size 7 to Power Series Ring in t over Integer Ring
        """
        try:
            P = x.parent()
            if is_PowerSeriesRing(P):
                if P.variable_name() == self.variable_name():
                    if self.has_coerce_map_from(P.base_ring()):
                        return self(x)
                    else:
                        raise TypeError("no natural map between bases of power series rings")

        except AttributeError:
            pass
        return self._coerce_try(x, [self.base_ring(), self.__poly_ring])



    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        r"""
        This gets called implicitly when one constructs a ring homomorphism
        from a power series ring.

        EXAMPLES::

            sage: S = RationalField(); R.<t>=PowerSeriesRing(S)
            sage: f = R.hom([0])
            sage: f(3)
            3
            sage: g = R.hom([t^2])
            sage: g(-1 + 3/5 * t)
            -1 + 3/5*t^2

        .. note::

           There are no ring homomorphisms from the ring of all formal
           power series to most rings, e.g, the p-adic field, since
           you can always (mathematically!) construct some power
           series that doesn't converge. Note that 0 is not a *ring*
           homomorphism.
        """
        if im_gens[0] == 0:
            return True   # this is allowed.
        if base_map is None and not codomain.has_coerce_map_from(self.base_ring()):
            return False
        from .laurent_series_ring import is_LaurentSeriesRing
        v = im_gens[0]
        if is_PowerSeriesRing(codomain) or is_LaurentSeriesRing(codomain):
            try:
                return v.valuation() > 0 or v.is_nilpotent()
            except NotImplementedError:
                return v.valuation() > 0
        try:
            return v.is_nilpotent()
        except NotImplementedError:
            pass
        return False

    def _poly_ring(self):
        """
        Return the underlying polynomial ring used to represent elements of
        this power series ring.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R._poly_ring()
            Univariate Polynomial Ring in t over Integer Ring
        """
        return self.__poly_ring

    def base_extend(self, R):
        """
        Return the power series ring over R in the same variable as self,
        assuming there is a canonical coerce map from the base ring of self
        to R.

        EXAMPLES::

            sage: R.<T> = GF(7)[[]]; R
            Power Series Ring in T over Finite Field of size 7
            sage: R.change_ring(ZZ)
            Power Series Ring in T over Integer Ring
            sage: R.base_extend(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: no base extension defined
        """
        if R.has_coerce_map_from(self.base_ring()):
            return self.change_ring(R)
        else:
            raise TypeError("no base extension defined")

    def change_ring(self, R):
        """
        Return the power series ring over R in the same variable as self.

        EXAMPLES::

            sage: R.<T> = QQ[[]]; R
            Power Series Ring in T over Rational Field
            sage: R.change_ring(GF(7))
            Power Series Ring in T over Finite Field of size 7
            sage: R.base_extend(GF(7))
            Traceback (most recent call last):
            ...
            TypeError: no base extension defined
            sage: R.base_extend(QuadraticField(3,'a'))
            Power Series Ring in T over Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878?
        """
        return PowerSeriesRing(R, name = self.variable_name(), default_prec = self.default_prec())

    def change_var(self, var):
        """
        Return the power series ring in variable ``var`` over the same base ring.

        EXAMPLES::

            sage: R.<T> = QQ[[]]; R
            Power Series Ring in T over Rational Field
            sage: R.change_var('D')
            Power Series Ring in D over Rational Field
        """
        return PowerSeriesRing(self.base_ring(), names = var, sparse=self.is_sparse())

    def is_exact(self):
        """
        Return False since the ring of power series over any ring is not
        exact.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R.is_exact()
            False
        """
        return False

    def gen(self, n=0):
        """
        Return the generator of this power series ring.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R.gen()
            t
            sage: R.gen(3)
            Traceback (most recent call last):
            ...
            IndexError: generator n>0 not defined
        """
        if n != 0:
            raise IndexError("generator n>0 not defined")
        return self.__generator

    def uniformizer(self):
        """
        Return a uniformizer of this power series ring if it is
        a discrete valuation ring (i.e., if the base ring is actually
        a field). Otherwise, an error is raised.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ)
            sage: R.uniformizer()
            t

            sage: R.<t> = PowerSeriesRing(ZZ)
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
        Return the number of generators of this power series ring.

        This is always 1.

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: R.ngens()
            1
        """
        return 1

    def random_element(self, prec=None, *args, **kwds):
        r"""
        Return a random power series.


        INPUT:

        -  ``prec`` - Integer specifying precision of output (default:
           default precision of self)

        -  ``*args, **kwds`` - Passed on to the ``random_element`` method for
           the base ring

        OUTPUT:

        -  Power series with precision ``prec`` whose coefficients are
           random elements from the base ring, randomized subject to the
           arguments ``*args`` and ``**kwds``


        ALGORITHM:

        Call the ``random_element`` method on the underlying polynomial
        ring.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ)
            sage: R.random_element(5)  # random
            -4 - 1/2*t^2 - 1/95*t^3 + 1/2*t^4 + O(t^5)
            sage: R.random_element(10)  # random
            -1/2 + 2*t - 2/7*t^2 - 25*t^3 - t^4 + 2*t^5 - 4*t^7 - 1/3*t^8 - t^9 + O(t^10)

        If given no argument, ``random_element`` uses default precision of self::

            sage: T = PowerSeriesRing(ZZ,'t')
            sage: T.default_prec()
            20
            sage: T.random_element()  # random
            4 + 2*t - t^2 - t^3 + 2*t^4 + t^5 + t^6 - 2*t^7 - t^8 - t^9 + t^11 - 6*t^12 + 2*t^14 + 2*t^16 - t^17 - 3*t^18 + O(t^20)
            sage: S = PowerSeriesRing(ZZ,'t', default_prec=4)
            sage: S.random_element()  # random
            2 - t - 5*t^2 + t^3 + O(t^4)


        Further arguments are passed to the underlying base ring (:trac:`9481`)::

            sage: SZ = PowerSeriesRing(ZZ,'v')
            sage: SQ = PowerSeriesRing(QQ,'v')
            sage: SR = PowerSeriesRing(RR,'v')

            sage: SZ.random_element(x=4, y=6)  # random
            4 + 5*v + 5*v^2 + 5*v^3 + 4*v^4 + 5*v^5 + 5*v^6 + 5*v^7 + 4*v^8 + 5*v^9 + 4*v^10 + 4*v^11 + 5*v^12 + 5*v^13 + 5*v^14 + 5*v^15 + 5*v^16 + 5*v^17 + 4*v^18 + 5*v^19 + O(v^20)
            sage: SZ.random_element(3, x=4, y=6)  # random
            5 + 4*v + 5*v^2 + O(v^3)
            sage: SQ.random_element(3, num_bound=3, den_bound=100)  # random
            1/87 - 3/70*v - 3/44*v^2 + O(v^3)
            sage: SR.random_element(3, max=10, min=-10)  # random
            2.85948321262904 - 9.73071330911226*v - 6.60414378519265*v^2 + O(v^3)

        """
        if prec is None:
            prec = self.default_prec()
        return self(self.__poly_ring.random_element(prec-1, *args, **kwds), prec)

    def __contains__(self, x):
        """
        Return True if x is an element of this power series ring or
        canonically coerces to this ring.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: t + t^2 in R
            True
            sage: 1/t in R
            False
            sage: 5 in R
            True
            sage: 1/3 in R
            False
            sage: S.<s> = PowerSeriesRing(ZZ)
            sage: s in R
            False
        """
        return self.has_coerce_map_from(parent(x))

    def is_field(self, proof = True):
        """
        Return False since the ring of power series over any ring is never
        a field.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R.is_field()
            False
        """
        return False

    def is_finite(self):
        """
        Return False since the ring of power series over any ring is never
        finite.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R.is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Return the characteristic of this power series ring, which is the
        same as the characteristic of the base ring of the power series
        ring.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R.characteristic()
            0
            sage: R.<w> = Integers(2^50)[[]]; R
            Power Series Ring in w over Ring of integers modulo 1125899906842624
            sage: R.characteristic()
            1125899906842624
        """
        return self.base_ring().characteristic()

    def residue_field(self):
        """
        Return the residue field of this power series ring.

        EXAMPLES::

            sage: R.<x> = PowerSeriesRing(GF(17))
            sage: R.residue_field()
            Finite Field of size 17
            sage: R.<x> = PowerSeriesRing(Zp(5))
            sage: R.residue_field()
            Finite Field of size 5
        """
        if self.base_ring().is_field():
            return self.base_ring()
        else:
            return self.base_ring().residue_field()

    def laurent_series_ring(self):
        """
        If this is the power series ring `R[[t]]`, return the
        Laurent series ring `R((t))`.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ,default_prec=5)
            sage: S = R.laurent_series_ring(); S
            Laurent Series Ring in t over Integer Ring
            sage: S.default_prec()
            5
            sage: f = 1+t; g=1/f; g
            1 - t + t^2 - t^3 + t^4 + O(t^5)
        """
        try:
            return self.__laurent_series_ring
        except AttributeError:
            self.__laurent_series_ring = laurent_series_ring.LaurentSeriesRing(
                                                 self.base_ring(), self.variable_name(), default_prec=self.default_prec(), sparse=self.is_sparse())
            return self.__laurent_series_ring


class PowerSeriesRing_domain(PowerSeriesRing_generic, ring.IntegralDomain):
    def fraction_field(self):
        """
        Return the Laurent series ring over the fraction field of the base
        ring.

        This is actually *not* the fraction field of this ring, but its
        completion with respect to the topology defined by the valuation.
        When we are working at finite precision, these two fields are
        indistinguishable; that is the reason why we allow ourselves to
        make this confusion here.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: R.fraction_field()
            Laurent Series Ring in t over Rational Field
            sage: Frac(R)
            Laurent Series Ring in t over Rational Field
        """
        laurent = self.laurent_series_ring()
        return laurent.change_ring(self.base_ring().fraction_field())


class PowerSeriesRing_over_field(PowerSeriesRing_domain):
    _default_category = CompleteDiscreteValuationRings()

    def fraction_field(self):
        """
        Return the fraction field of this power series ring, which is
        defined since this is over a field.

        This fraction field is just the Laurent series ring over the base
        field.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(GF(7))
            sage: R.fraction_field()
            Laurent Series Ring in t over Finite Field of size 7
            sage: Frac(R)
            Laurent Series Ring in t over Finite Field of size 7
        """
        return self.laurent_series_ring()


def unpickle_power_series_ring_v0(base_ring, name, default_prec, sparse):
    """
    Unpickle (deserialize) a univariate power series ring according to
    the given inputs.

    EXAMPLES::

        sage: P.<x> = PowerSeriesRing(QQ)
        sage: loads(dumps(P)) == P # indirect doctest
        True
    """
    return PowerSeriesRing(base_ring, name=name, default_prec = default_prec, sparse=sparse)
