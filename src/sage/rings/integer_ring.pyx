r"""
Ring `\ZZ` of Integers

The :class:`IntegerRing_class` represents the ring `\ZZ` of (arbitrary
precision) integers. Each integer is an instance of :class:`Integer`,
which is defined in a Pyrex extension module that wraps GMP integers
(the ``mpz_t`` type in GMP).

::

    sage: Z = IntegerRing(); Z
    Integer Ring
    sage: Z.characteristic()
    0
    sage: Z.is_field()
    False

There is a unique instance of the :class:`integer ring<IntegerRing_class>`.
To create an :class:`Integer`, coerce either a Python int, long, or a string. Various
other types will also coerce to the integers, when it makes sense.

::

    sage: a = Z(1234); a
    1234
    sage: b = Z(5678); b
    5678
    sage: type(a)
    <class 'sage.rings.integer.Integer'>
    sage: a + b
    6912
    sage: Z('94803849083985934859834583945394')
    94803849083985934859834583945394
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.long cimport *
from cpython.list cimport *
from cpython.object cimport Py_NE

from cysignals.signals cimport sig_check, sig_on, sig_off

from sage.libs.gmp.mpz cimport *
import sage.rings.infinity
import sage.rings.rational
import sage.rings.rational_field
import sage.rings.ideal
from sage.categories.basic import EuclideanDomains, DedekindDomains
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.noetherian_rings import NoetherianRings
from sage.rings.number_field.number_field_element_base import NumberFieldElement_base
from sage.structure.coerce cimport is_numpy_type
from sage.structure.element cimport parent
from sage.structure.parent cimport Parent
from sage.structure.richcmp cimport rich_to_bool

from sage.misc.misc_c import prod
from sage.misc.randstate cimport randstate, current_randstate, SAGE_RAND_MAX

cimport sage.rings.integer as integer

from sage.rings import ring

arith = None
cdef void late_import() noexcept:
    # A hack to avoid circular imports.
    global arith
    if arith is None:
        import sage.arith.all
        arith = sage.arith.all

cdef int number_of_integer_rings = 0

# Used by IntegerRing_class._randomize_mpz():
#   When the user requests an integer sampled from a
#   DiscreteGaussianDistributionIntegerSampler with parameter sigma, we store the pair
#   (sigma, sampler) here.  When the user requests an integer for the same
#   sigma, we do not recreate the sampler but take it from this "cache".
_prev_discrete_gaussian_integer_sampler = (None, None)


def is_IntegerRing(x):
    r"""
    Internal function: return ``True`` iff ``x`` is the ring `\ZZ` of integers.

    TESTS::

        sage: from sage.rings.integer_ring import is_IntegerRing
        sage: is_IntegerRing(ZZ)
        doctest:warning...
        DeprecationWarning: The function is_IntegerRing is deprecated;
        use 'isinstance(..., IntegerRing_class)' instead.
        See https://github.com/sagemath/sage/issues/38128 for details.
        True
        sage: is_IntegerRing(QQ)
        False
        sage: is_IntegerRing(parent(3))
        True
        sage: is_IntegerRing(parent(1/3))
        False
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38128,
                       "The function is_IntegerRing is deprecated; "
                       "use 'isinstance(..., IntegerRing_class)' instead.")
    return isinstance(x, IntegerRing_class)


cdef class IntegerRing_class(CommutativeRing):
    r"""
    The ring of integers.

    In order to introduce the ring `\ZZ` of integers, we illustrate
    creation, calling a few functions, and working with its elements.

    ::

        sage: Z = IntegerRing(); Z
        Integer Ring
        sage: Z.characteristic()
        0
        sage: Z.is_field()
        False
        sage: Z.category()
        Join of Category of Dedekind domains
            and Category of euclidean domains
            and Category of noetherian rings
            and Category of infinite enumerated sets
            and Category of metric spaces
        sage: Z(2^(2^5) + 1)
        4294967297

    One can give strings to create integers. Strings starting with
    ``0x`` are interpreted as hexadecimal, and strings starting with
    ``0o`` are interpreted as octal::

        sage: parent('37')
        <... 'str'>
        sage: parent(Z('37'))
        Integer Ring
        sage: Z('0x10')
        16
        sage: Z('0x1a')
        26
        sage: Z('0o20')
        16

    As an inverse to :meth:`~sage.rings.integer.Integer.digits`,
    lists of digits are accepted, provided that you give a base.
    The lists are interpreted in little-endian order, so that
    entry ``i`` of the list is the coefficient of ``base^i``::

        sage: Z([4,1,7], base=100)
        70104
        sage: Z([4,1,7], base=10)
        714
        sage: Z([3, 7], 10)
        73
        sage: Z([3, 7], 9)
        66
        sage: Z([], 10)
        0

    Alphanumeric strings can be used for bases 2..36; letters ``a`` to
    ``z`` represent numbers 10 to 36.  Letter case does not matter.
    ::

        sage: Z("sage", base=32)
        928270
        sage: Z("SAGE", base=32)
        928270
        sage: Z("Sage", base=32)
        928270
        sage: Z([14, 16, 10, 28], base=32)
        928270
        sage: 14 + 16*32 + 10*32^2 + 28*32^3
        928270

    We next illustrate basic arithmetic in `\ZZ`::

        sage: a = Z(1234); a
        1234
        sage: b = Z(5678); b
        5678
        sage: type(a)
        <class 'sage.rings.integer.Integer'>
        sage: a + b
        6912
        sage: b + a
        6912
        sage: a * b
        7006652
        sage: b * a
        7006652
        sage: a - b
        -4444
        sage: b - a
        4444

    When we divide two integers using ``/``, the result is automatically
    coerced to the field of rational numbers, even if the result is
    an integer.

    ::

        sage: a / b
        617/2839
        sage: type(a/b)
        <class 'sage.rings.rational.Rational'>
        sage: a/a
        1
        sage: type(a/a)
        <class 'sage.rings.rational.Rational'>

    For floor division, use the ``//`` operator instead::

        sage: a // b
        0
        sage: type(a//b)
        <class 'sage.rings.integer.Integer'>

    Next we illustrate arithmetic with automatic coercion. The types
    that coerce are: str, int, long, Integer.

    ::

        sage: a + 17
        1251
        sage: a * 374
        461516
        sage: 374 * a
        461516
        sage: a/19
        1234/19
        sage: 0 + Z(-64)
        -64

    Integers can be coerced::

        sage: a = Z(-64)
        sage: int(a)
        -64

    We can create integers from several types of objects::

        sage: Z(17/1)
        17
        sage: Z(Mod(19,23))
        19
        sage: Z(2 + 3*5 + O(5^3))                                                       # needs sage.rings.padics
        17

    Arbitrary numeric bases are supported; strings or list of integers
    are used to provide the digits (more details in
    :class:`IntegerRing_class`)::

        sage: Z("sage", base=32)
        928270
        sage: Z([14, 16, 10, 28], base=32)
        928270

    The :meth:`digits<~sage.rings.integer.Integer.digits>` method
    allows you to get the list of digits of an integer in a different
    basis (note that the digits are returned in little-endian order)::

        sage: b = Z([4,1,7], base=100)
        sage: b
        70104
        sage: b.digits(base=71)
        [27, 64, 13]

        sage: Z(15).digits(2)
        [1, 1, 1, 1]
        sage: Z(15).digits(3)
        [0, 2, 1]

    The :meth:`str<~sage.rings.integer.Integer.str>` method returns a
    string of the digits, using letters ``a`` to ``z`` to represent
    digits 10..36::

        sage: Z(928270).str(base=32)
        'sage'

    Note that :meth:`str<~sage.rings.integer.Integer.str>` only works
    with bases 2 through 36.

    TESTS::

        sage: TestSuite(ZZ).run()
        sage: list(ZZ)
        Traceback (most recent call last):
        ...
        NotImplementedError: len() of an infinite set

        sage: ZZ.is_finite()
        False
        sage: ZZ.cardinality()
        +Infinity
    """

    def __init__(self):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.rings.integer_ring import IntegerRing_class
            sage: A = IntegerRing_class()

        We check that ``ZZ`` is an infinite enumerated set
        (see :issue:`16239`)::

            sage: A in InfiniteEnumeratedSets()
            True
        """
        cat = (EuclideanDomains(), DedekindDomains(),
               InfiniteEnumeratedSets().Metric(), NoetherianRings())
        Parent.__init__(self, base=self, names=('x',), normalize=False,
                        category=cat)
        self._populate_coercion_lists_(init_no_parent=True,
                                       convert_method_name='_integer_')

    _element_constructor_ = integer.Integer

    def __reduce__(self):
        """
        For pickling.

        TESTS::

            sage: loads(dumps(ZZ)) is ZZ
            True
        """
        return IntegerRing, ()

    def __hash__(self):
        """
        Return the hash value of ``self``.

        TESTS::

            sage: from sage.rings.integer_ring import IntegerRing_class
            sage: A = IntegerRing_class()
            sage: A.__hash__()
            554590422
        """
        return 554590422

    def __richcmp__(left, right, int op):
        """
        Rich comparison of ``left`` and ``right``.

        TESTS::

            sage: from sage.rings.integer_ring import IntegerRing_class
            sage: ZZ == ZZ
            True
            sage: ZZ != QQ
            True
        """
        if left is right:
            return rich_to_bool(op, 0)

        if isinstance(right, IntegerRing_class):
            return rich_to_bool(op, 0)

        if isinstance(right, sage.rings.rational_field.RationalField):
            return rich_to_bool(op, -1)

        return op == Py_NE

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ZZ # indirect doctest
            Integer Ring
        """
        return "Integer Ring"

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        TESTS::

            sage: latex(ZZ) # indirect doctest
            \Bold{Z}
        """
        return "\\Bold{Z}"

    def __getitem__(self, x):
        r"""
        Return the ring `\ZZ[...]` obtained by adjoining to the integers one
        or several elements.

        EXAMPLES::

            sage: ZZ['x']
            Univariate Polynomial Ring in x over Integer Ring
            sage: ZZ['x,y']
            Multivariate Polynomial Ring in x, y over Integer Ring

            sage: # needs sage.rings.number_field sage.symbolic
            sage: ZZ[sqrt(2), sqrt(3)]
            Relative Order generated by [-sqrt3*sqrt2 + 3, 6*sqrt2 - 5*sqrt3, 11*sqrt2 - 9*sqrt3]
             in Number Field in sqrt2 with defining polynomial x^2 - 2 over its base field
            sage: R = ZZ[sqrt(5) + 1]; R
            Order of conductor 2 generated by a in Number Field in a
             with defining polynomial x^2 - 2*x - 4 with a = 3.236067977499790?
            sage: R.is_maximal()
            False
            sage: R = ZZ[(1 + sqrt(5))/2]; R
            Maximal Order generated by a in Number Field in a
             with defining polynomial x^2 - x - 1 with a = 1.618033988749895?
            sage: R.is_maximal()
            True
        """
        if x in self:
            return self

        if isinstance(x, NumberFieldElement_base):
            K, _ = parent(x).subfield(x)
            return K.order(K.gen())

        return CommutativeRing.__getitem__(self, x)

    def range(self, start, end=None, step=None):
        """
        Optimized range function for Sage integers.

        AUTHORS:

        - Robert Bradshaw (2007-09-20)

        EXAMPLES::

            sage: ZZ.range(10)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: ZZ.range(-5, 5)
            [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4]
            sage: ZZ.range(0, 50, 5)
            [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
            sage: ZZ.range(0, 50, -5)
            []
            sage: ZZ.range(50, 0, -5)
            [50, 45, 40, 35, 30, 25, 20, 15, 10, 5]
            sage: ZZ.range(50, 0, 5)
            []
            sage: ZZ.range(50, -1, -5)
            [50, 45, 40, 35, 30, 25, 20, 15, 10, 5, 0]

        It uses different code if the step doesn't fit in a long::

            sage: ZZ.range(0, 2^83, 2^80)
            [0, 1208925819614629174706176, 2417851639229258349412352,
             3626777458843887524118528, 4835703278458516698824704, 6044629098073145873530880,
             7253554917687775048237056, 8462480737302404222943232]

        Make sure :issue:`8818` is fixed::

            sage: ZZ.range(1r, 10r)
            [1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        if end is None:
            end = start
            start = Integer.__new__(Integer)
        if step is None:
            step = 1
        if type(step) is not int:
            if not isinstance(step, integer.Integer):
                step = integer.Integer(step)
            if mpz_fits_slong_p((<Integer>step).value):
                step = int(step)
        if not isinstance(start, integer.Integer):
            start = integer.Integer(start)
        if not isinstance(end, integer.Integer):
            end = integer.Integer(end)
        cdef integer.Integer a = <Integer>start
        cdef integer.Integer b = <Integer>end

        cdef int step_sign
        cdef long istep
        cdef integer.Integer zstep, last

        L = []
        if type(step) is int:
            istep = PyLong_AsLong(step)
            step_sign = istep
        else:
            zstep = <Integer>step
            step_sign = mpz_sgn(zstep.value)

        sig_on()
        while mpz_cmp(a.value, b.value)*step_sign < 0:
            last = a
            a = Integer.__new__(Integer)
            if type(step) is int: # count on branch prediction...
                if istep > 0:
                    mpz_add_ui(a.value, last.value, istep)
                else:
                    mpz_sub_ui(a.value, last.value, -istep)
            else:
                mpz_add(a.value, last.value, zstep.value)
            PyList_Append(L, last)
        sig_off()
        return L

    def __iter__(self):
        """
        Iterate over all integers: 0 1 -1 2 -2 3 -3 ...

        EXAMPLES::

            sage: for n in ZZ:
            ....:  if n < 3: print(n)
            ....:  else: break
            0
            1
            -1
            2
            -2
        """
        yield self(0)
        n = self(1)
        while True:
            sig_check()
            yield n
            yield -n
            n += 1

    cpdef _coerce_map_from_(self, S):
        r"""
        ``x`` canonically coerces to the integers `\ZZ` only if ``x``
        is an int, long or already an element of `\ZZ`.

        EXAMPLES::

            sage: ZZ.coerce(int(5)) # indirect doctest
            5
            sage: ZZ.coerce(GF(7)(2))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Finite Field of size 7 to Integer Ring

        The rational number ``3/1 == 3`` does not canonically coerce into the
        integers, since there is no canonical coercion map from the full
        field of rational numbers to the integers.

        ::

            sage: a = 3/1; parent(a)
            Rational Field
            sage: ZZ(a)
            3
            sage: ZZ.coerce(a)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Integer Ring

        Coercions are available from numpy integer types::

            sage: # needs numpy
            sage: import numpy
            sage: ZZ.coerce(numpy.int8('1'))
            1
            sage: ZZ.coerce(numpy.int32('32'))
            32
            sage: ZZ.coerce(numpy.int64('-12'))
            -12
            sage: ZZ.coerce(numpy.uint64('11'))
            11

        TESTS::

            sage: 5r + True
            6
            sage: 5 + True
            6

            sage: f = ZZ.coerce_map_from(int); f
            Native morphism:
              From: Set of Python objects of class 'int'
              To:   Integer Ring
            sage: f(4r)
            4
            sage: f(-7r)
            -7
        """
        if S is int:
            return sage.rings.integer.int_to_Z()
        elif S is bool:
            return True
        elif is_numpy_type(S):
            import numpy
            if issubclass(S, numpy.integer):
                return True
        return None

    def random_element(self, x=None, y=None, distribution=None):
        r"""
        Return a random integer.

        INPUT:

        - ``x``, ``y`` integers -- bounds for the result

        - ``distribution`` -- string:

          - ``'uniform'``
          - ``'mpz_rrandomb'``
          - ``'1/n'``
          - ``'gaussian'``

        OUTPUT: with no input, return a random integer

          If only one integer `x` is given, return an integer
          between 0 and `x-1`.

          If two integers are given, return an integer
          between `x` and `y-1` inclusive.

          If at least one bound is given, the default distribution is the
          uniform distribution; otherwise, it is the distribution described
          below.

          If the distribution ``'1/n'`` is specified, the bounds are ignored.

          If the distribution ``'mpz_rrandomb'`` is specified, the output is
          in the range from 0 to `2^x - 1`.

          If the distribution ``'gaussian'`` is specified, the output is
          sampled from a discrete Gaussian distribution with parameter
          `\sigma=x` and centered at zero. That is, the integer `v` is returned
          with probability proportional to `\mathrm{exp}(-v^2/(2\sigma^2))`.
          See :mod:`sage.stats.distributions.discrete_gaussian_integer` for
          details.  Note that if many samples from the same discrete Gaussian
          distribution are needed, it is faster to construct a
          :class:`sage.stats.distributions.discrete_gaussian_integer.DiscreteGaussianDistributionIntegerSampler`
          object which is then repeatedly queried.

        The default distribution for ``ZZ.random_element()`` is based on
        `X = \mbox{trunc}(4/(5R))`, where `R` is a random
        variable uniformly distributed between `-1` and `1`. This gives
        `\mbox{Pr}(X = 0) = 1/5`, and
        `\mbox{Pr}(X = n) = 2/(5|n|(|n|+1))` for
        `n \neq 0`. Most of the samples will be small; `-1`, `0`, and `1`
        occur with probability `1/5` each. But we also have a small but
        non-negligible proportion of "outliers";
        `\mbox{Pr}(|X| \geq n) = 4/(5n)`, so for instance, we
        expect that `|X| \geq 1000` on one in 1250 samples.

        We actually use an easy-to-compute truncation of the above
        distribution; the probabilities given above hold fairly well up to
        about `|n| = 10000`, but around `|n| = 30000` some
        values will never be returned at all, and we will never return
        anything greater than `2^{30}`.

        EXAMPLES::

            sage: ZZ.random_element().parent() is ZZ
            True

        The default uniform distribution is integers in `[-2, 2]`::

            sage: from collections import defaultdict
            sage: def add_samples(*args, **kwds):
            ....:     global dic, counter
            ....:     for _ in range(100):
            ....:         counter += 1
            ....:         dic[ZZ.random_element(*args, **kwds)] += 1

            sage: def prob(x):
            ....:     return 1/5
            sage: dic = defaultdict(Integer)
            sage: counter = 0.0
            sage: add_samples(distribution='uniform')
            sage: while any(abs(dic[i]/counter - prob(i)) > 0.01 for i in dic):
            ....:     add_samples(distribution='uniform')

        Here we use the distribution ``'1/n'``::

            sage: def prob(n):
            ....:     if n == 0:
            ....:         return 1/5
            ....:     return 2/(5*abs(n)*(abs(n) + 1))
            sage: dic = defaultdict(Integer)
            sage: counter = 0.0
            sage: add_samples(distribution='1/n')
            sage: while any(abs(dic[i]/counter - prob(i)) > 0.01 for i in dic):
            ....:     add_samples(distribution='1/n')

        If a range is given, the default distribution is uniform in that
        range::

            sage: -10 <= ZZ.random_element(-10, 10) < 10
            True
            sage: def prob(x):
            ....:     return 1/20
            sage: dic = defaultdict(Integer)
            sage: counter = 0.0
            sage: add_samples(-10, 10)
            sage: while any(abs(dic[i]/counter - prob(i)) > 0.01 for i in dic):
            ....:     add_samples(-10, 10)

            sage: 0 <= ZZ.random_element(5) < 5
            True
            sage: def prob(x):
            ....:     return 1/5
            sage: dic = defaultdict(Integer)
            sage: counter = 0.0
            sage: add_samples(5)
            sage: while any(abs(dic[i]/counter - prob(i)) > 0.01 for i in dic):
            ....:     add_samples(5)

            sage: while ZZ.random_element(10^50) < 10^49:
            ....:     pass

        Notice that the right endpoint is not included::

            sage: all(ZZ.random_element(-2, 2) < 2 for _ in range(100))
            True

        We return a sample from a discrete Gaussian distribution::

             sage: ZZ.random_element(11.0, distribution='gaussian').parent() is ZZ      # needs sage.modules
             True

        TESTS:

        Check that :issue:`32124` is fixed::

            sage: ZZ.random_element(5, -5, distribution='1/n').parent() is ZZ
            True
            sage: ZZ.random_element(5, -5, distribution='gaussian').parent() is ZZ      # needs sage.modules
            True
            sage: ZZ.random_element(5, -5, distribution='mpz_rrandomb').parent() is ZZ
            True

            sage: ZZ.random_element(-10, -5, distribution='mpz_rrandomb')
            Traceback (most recent call last):
            ...
            TypeError: x must be > 0
            sage: ZZ.random_element(-10, -5, distribution='gaussian')
            Traceback (most recent call last):
            ...
            TypeError: x must be > 0

        Checking error messages::

            sage: ZZ.random_element(-3)
            Traceback (most recent call last):
            ...
            TypeError: x must be > 0
            sage: ZZ.random_element(4, 2)
            Traceback (most recent call last):
            ...
            TypeError: x must be < y
        """
        cdef Integer z = Integer.__new__(Integer)
        if distribution == "1/n":
            x = None
            y = None
        elif distribution == "mpz_rrandomb" or distribution == "gaussian":
            y = None
        if x is not None and y is None and x <= 0:
            raise TypeError("x must be > 0")
        if x is not None and y is not None and x >= y:
            raise TypeError("x must be < y")
        self._randomize_mpz(z.value, x, y, distribution)
        return z

    cdef int _randomize_mpz(self, mpz_t value, x, y, distribution) except -1:
        """
        This is an internal function, called by random_element.

        INPUT:

        - ``value`` -- this is the variable in which the answer will be
          returned

        - ``x``, ``y``, ``distribution`` -- see :meth:`random_element`

        TESTS::

            sage: ZZ.random_element() # indirect doctest # random
            6
        """
        cdef integer.Integer r
        cdef integer.Integer n_max, n_min, n_width
        cdef randstate rstate = current_randstate()
        cdef int den = rstate.c_random()-SAGE_RAND_MAX/2
        if den == 0: den = 1
        if (distribution is None and x is None) or distribution == "1/n":
            mpz_set_si(value, (SAGE_RAND_MAX/5*2) / den)
        elif distribution is None or distribution == "uniform":
            if y is None:
                if x is None:
                    mpz_set_si(value, rstate.c_random()%5 - 2)
                else:
                    n_max = x if isinstance(x, integer.Integer) else self(x)
                    mpz_urandomm(value, rstate.gmp_state, n_max.value)
            else:
                n_min = x if isinstance(x, integer.Integer) else self(x)
                n_max = y if isinstance(y, integer.Integer) else self(y)
                n_width = n_max - n_min
                if mpz_sgn(n_width.value) <= 0:
                    n_min = self(-2)
                    n_width = self(5)
                mpz_urandomm(value, rstate.gmp_state, n_width.value)
                mpz_add(value, value, n_min.value)
        elif distribution == "mpz_rrandomb":
            if x is None:
                raise ValueError("must specify x to use 'distribution=mpz_rrandomb'")
            mpz_rrandomb(value, rstate.gmp_state, int(x))
        elif distribution == "gaussian":
            global _prev_discrete_gaussian_integer_sampler
            if x == _prev_discrete_gaussian_integer_sampler[0]:
                r = _prev_discrete_gaussian_integer_sampler[1]()
            else:
                from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
                D = DiscreteGaussianDistributionIntegerSampler(sigma=x, algorithm='uniform+logtable')
                r = D()
                _prev_discrete_gaussian_integer_sampler = (x, D)
            mpz_set(value, r.value)
        else:
            raise ValueError("Unknown distribution for the integers: %s" % distribution)

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        r"""
        Test whether the map from `\ZZ` to codomain, which takes the
        generator of `\ZZ` to ``im_gens[0]``, is a ring homomorphism.
        (This amounts to checking that 1 goes to 1.)

        TESTS::

            sage: ZZ._is_valid_homomorphism_(ZZ,[1])
            True
            sage: ZZ._is_valid_homomorphism_(ZZ,[2])
            False
            sage: ZZ._is_valid_homomorphism_(ZZ.quotient_ring(8), [ZZ.quotient_ring(8)(1)])
            True
        """
        if base_map is None:
            base_map = codomain.coerce_map_from(self)
            if base_map is None:
                return False
        try:
            return im_gens[0] == base_map(self.gen(0))
        except TypeError:
            return False

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: ZZ._repr_option('element_is_atomic')
            True
        """
        if key == 'element_is_atomic':
            return True
        return super()._repr_option(key)

    def is_field(self, proof = True):
        """
        Return ``False`` since the integers are not a field.

        EXAMPLES::

            sage: ZZ.is_field()
            False
        """
        return False

    def fraction_field(self):
        """
        Return the field of rational numbers - the fraction field of the
        integers.

        EXAMPLES::

            sage: ZZ.fraction_field()
            Rational Field
            sage: ZZ.fraction_field() == QQ
            True
        """
        return sage.rings.rational_field.Q

    def extension(self, poly, names, **kwds):
        """
        Return the order generated by the specified list of polynomials.

        INPUT:

        - ``poly`` -- list of one or more polynomials

        - ``names`` -- a parameter which will be passed to
          :func:`EquationOrder`

        - ``embedding`` -- a parameter which will be passed to
          :func:`EquationOrder`

        OUTPUT:

        - Given a single polynomial as input, return the order generated by a
          root of the polynomial in the field generated by a root of the
          polynomial.

          Given a list of polynomials as input, return the relative order
          generated by a root of the first polynomial in the list, over the
          order generated by the roots of the subsequent polynomials.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: ZZ.extension(x^2 - 5, 'a')                                            # needs sage.rings.number_field
            Order of conductor 2 generated by a in Number Field in a with defining polynomial x^2 - 5
            sage: ZZ.extension([x^2 + 1, x^2 + 2], 'a,b')                               # needs sage.rings.number_field
            Relative Order generated by [-b*a - 1, -3*a + 2*b] in Number Field in a
             with defining polynomial x^2 + 1 over its base field
        """
        from sage.rings.number_field.order import EquationOrder
        return EquationOrder(poly, names=names, **kwds)

    def quotient(self, I, names=None, **kwds):
        r"""
        Return the quotient of `\ZZ` by the ideal or integer ``I``.

        EXAMPLES::

            sage: ZZ.quo(6*ZZ)
            Ring of integers modulo 6
            sage: ZZ.quo(0*ZZ)
            Integer Ring
            sage: ZZ.quo(3)
            Ring of integers modulo 3
            sage: ZZ.quo(3*QQ)
            Traceback (most recent call last):
            ...
            TypeError: I must be an ideal of ZZ
        """
        if isinstance(I, sage.rings.integer.Integer):
            n = I
        elif isinstance(I, sage.rings.ideal.Ideal_generic):
            if not (I.ring() is self):
                raise TypeError("I must be an ideal of ZZ")
            n = I.gens()[0]
        else:
            raise TypeError("I must be an ideal of ZZ or an integer")
        if n == 0:
            return self
        return sage.rings.finite_rings.integer_mod_ring.IntegerModRing(n, **kwds)

    def residue_field(self, prime, check=True, names=None):
        r"""
        Return the residue field of the integers modulo the given prime, i.e.
        `\ZZ/p\ZZ`.

        INPUT:

        - ``prime`` -- a prime number

        - ``check`` -- boolean (default: ``True``); whether or not
          to check the primality of prime

        - ``names`` -- ignored (for compatibility with number fields)

        OUTPUT: the residue field at this prime

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: F = ZZ.residue_field(61); F
            Residue field of Integers modulo 61
            sage: pi = F.reduction_map(); pi
            Partially defined reduction map:
              From: Rational Field
              To:   Residue field of Integers modulo 61
            sage: pi(123/234)
            6
            sage: pi(1/61)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot reduce rational 1/61 modulo 61:
            it has negative valuation
            sage: lift = F.lift_map(); lift
            Lifting map:
              From: Residue field of Integers modulo 61
              To:   Integer Ring
            sage: lift(F(12345/67890))
            33
            sage: (12345/67890) % 61
            33

        Construction can be from a prime ideal instead of a prime::

            sage: ZZ.residue_field(ZZ.ideal(97))
            Residue field of Integers modulo 97

        TESTS::

            sage: ZZ.residue_field(ZZ.ideal(96))
            Traceback (most recent call last):
            ...
            TypeError: Principal ideal (96) of Integer Ring is not prime
            sage: ZZ.residue_field(96)
            Traceback (most recent call last):
            ...
            TypeError: 96 is not prime
        """
        if isinstance(prime, sage.rings.integer.Integer):
            p = self.ideal(prime)
        elif isinstance(prime, sage.rings.ideal.Ideal_generic):
            if not (prime.ring() is self):
                raise TypeError("%s is not an ideal of ZZ" % prime)
            p = prime
        else:
            raise TypeError("%s is neither an ideal of ZZ nor an integer" % prime)
        if check and not p.is_prime():
            raise TypeError("%s is not prime" % prime)
        from sage.rings.finite_rings.residue_field import ResidueField
        return ResidueField(p, names = None, check = check)

    def gens(self):
        """
        Return the tuple ``(1,)`` containing a single element, the additive
        generator of the integers, which is 1.

        EXAMPLES::

            sage: ZZ.gens(); ZZ.gens()[0]
            (1,)
            1
            sage: type(ZZ.gens()[0])
            <class 'sage.rings.integer.Integer'>
        """
        return (self(1), )

    def gen(self, n=0):
        """
        Return the additive generator of the integers, which is 1.

        INPUT:

        - ``n`` -- (default: 0) in a ring with more than one generator, the
          optional parameter `n` indicates which generator to return; since
          there is only one generator in this case, the only valid value for
          `n` is 0

        EXAMPLES::

            sage: ZZ.gen()
            1
            sage: type(ZZ.gen())
            <class 'sage.rings.integer.Integer'>
        """
        if n == 0:
            return self(1)
        else:
            raise IndexError("n must be 0")

    def ngens(self):
        """
        Return the number of additive generators of the ring, which is 1.

        EXAMPLES::

            sage: ZZ.ngens()
            1
            sage: len(ZZ.gens())
            1
        """
        return 1

    def degree(self):
        """
        Return the degree of the integers, which is 1.

        Here, degree refers to the rank of the ring as a module over the
        integers.

        EXAMPLES::

            sage: ZZ.degree()
            1
        """
        return 1

    def absolute_degree(self):
        """
        Return the absolute degree of the integers, which is 1.

        Here, absolute degree refers to the rank of the ring as a module
        over the integers.

        EXAMPLES::

            sage: ZZ.absolute_degree()
            1
        """
        return 1

    def characteristic(self):
        """
        Return the characteristic of the integers, which is 0.

        EXAMPLES::

            sage: ZZ.characteristic()
            0
        """
        return ZZ.zero()

    def krull_dimension(self):
        """
        Return the Krull dimension of the integers, which is 1.

        .. NOTE::

            This should rather be inherited from the category
            of ``DedekindDomains``.

        EXAMPLES::

            sage: ZZ.krull_dimension()
            1
        """
        return 1

    def completion(self, p, prec, extras = {}):
        r"""
        Return the metric completion of the integers at the prime `p`.

        INPUT:

        - ``p`` -- a prime (or ``infinity``)

        - ``prec`` -- the desired precision

        - ``extras`` -- any further parameters to pass to the method used to
          create the completion

        OUTPUT: the completion of `\ZZ` at `p`

        EXAMPLES::

            sage: ZZ.completion(infinity, 53)
            Integer Ring
            sage: ZZ.completion(5, 15, {'print_mode': 'bars'})                          # needs sage.rings.padics
            5-adic Ring with capped relative precision 15
        """
        if p == sage.rings.infinity.Infinity:
            return self
        else:
            from sage.rings.padics.factory import Zp
            return Zp(p, prec, **extras)

    def order(self):
        """
        Return the order (cardinality) of the integers, which is
        ``+Infinity``.

        EXAMPLES::

            sage: ZZ.order()
            +Infinity
        """
        return sage.rings.infinity.infinity

    def zeta(self, n=2):
        r"""
        Return a primitive ``n``-th root of unity in the integers, or raise an
        error if none exists.

        INPUT:

        - ``n`` -- (default: 2) a positive integer

        OUTPUT:

        an ``n``-th root of unity in `\ZZ`

        EXAMPLES::

            sage: ZZ.zeta()
            -1
            sage: ZZ.zeta(1)
            1
            sage: ZZ.zeta(3)
            Traceback (most recent call last):
            ...
            ValueError: no nth root of unity in integer ring
            sage: ZZ.zeta(0)
            Traceback (most recent call last):
            ...
            ValueError: n must be positive in zeta()
        """
        if n == 1:
            return sage.rings.integer.Integer(1)
        elif n == 2:
            return sage.rings.integer.Integer(-1)
        elif n < 1:
            raise ValueError("n must be positive in zeta()")
        else:
            raise ValueError("no nth root of unity in integer ring")

    def parameter(self):
        r"""
        Return an integer of degree 1 for the Euclidean property of `\ZZ`,
        namely 1.

        EXAMPLES::

            sage: ZZ.parameter()
            1
        """
        return self(1)

    def _roots_univariate_polynomial(self, p, ring=None, multiplicities=True, algorithm=None):
        r"""
        Return the roots of the univariate polynomial ``p``.

        INPUT:

        - ``p`` -- univariate integer polynomial

        - ``ring`` -- ring (default: ``None``); a ring containing `\ZZ` to
          compute the roots in; ``None`` is equivalent to ``ZZ``

        - ``multiplicities`` -- boolean (default: ``True``); whether to
          compute the multiplicities

        - ``algorithm`` -- ``'dense'``, ``'sparse'`` or ``None`` (default:
          ``None``); the algorithm to use

        OUTPUT:

        - If ``multiplicities=True``, the list of pairs `(r, n)` where
          `r` is a root and `n` the corresponding multiplicity;

        - If ``multiplicities=False``, the list of distincts roots with no
          information about the multiplicities.

        ALGORITHM:

        If ``algorithm`` is ``'dense'``, the roots are computed using
        :meth:`_roots_from_factorization`.

        If ``algorithm`` is ``'sparse'``, the roots are computed using the
        algorithm described in [CKS1999]_.

        If ``algorithm`` is ``None``, use the ``'dense'`` algorithm for
        polynomials of degree at most `100`, and ``'sparse'`` otherwise.

        .. NOTE::

            This is a helper method for
            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.roots`.

        TESTS::

            sage: # needs sage.libs.pari
            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: p = (x + 1)^23 * (x - 1)^23 * (x - 100) * (x + 5445)^5
            sage: ZZ._roots_univariate_polynomial(p)
            [(100, 1), (-5445, 5), (1, 23), (-1, 23)]
            sage: p *= (1 + x^3458645 - 76*x^3435423343 + x^45346567867756556)
            sage: ZZ._roots_univariate_polynomial(p)
            [(1, 23), (-1, 23), (100, 1), (-5445, 5)]
            sage: p *= x^156468451540687043504386074354036574634735074
            sage: ZZ._roots_univariate_polynomial(p)
            [(0, 156468451540687043504386074354036574634735074),
             (1, 23),
             (-1, 23),
             (100, 1),
             (-5445, 5)]
            sage: ZZ._roots_univariate_polynomial(p, multiplicities=False)
            [0, 1, -1, 100, -5445]

            sage: # needs sage.libs.pari
            sage: R.<x> = PolynomialRing(ZZ, sparse=False)
            sage: p = (x + 1)^23 * (x - 1)^23 * (x - 100) * (x + 5445)^5
            sage: ZZ._roots_univariate_polynomial(p)
            [(100, 1), (-5445, 5), (1, 23), (-1, 23)]
            sage: ZZ._roots_univariate_polynomial(p, multiplicities=False)
            [100, -5445, 1, -1]
            sage: ZZ._roots_univariate_polynomial(p, algorithm='sparse')
            [(100, 1), (-5445, 5), (1, 23), (-1, 23)]
            sage: ZZ._roots_univariate_polynomial(p, algorithm='dense')
            [(100, 1), (-5445, 5), (1, 23), (-1, 23)]
            sage: ZZ._roots_univariate_polynomial(p, algorithm='foobar')
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm 'foobar'

            sage: p = x^20 * p                                                          # needs sage.libs.pari
            sage: ZZ._roots_univariate_polynomial(p, algorithm='sparse')                # needs sage.libs.pari
            [(0, 20), (100, 1), (-5445, 5), (1, 23), (-1, 23)]
            sage: ZZ._roots_univariate_polynomial(p, algorithm='dense')                 # needs sage.libs.pari
            [(100, 1), (-5445, 5), (0, 20), (1, 23), (-1, 23)]
        """
        deg = p.degree()
        if deg < 0:
            raise ValueError("roots of 0 are not defined")

        # A specific algorithm is available only for integer roots of integer polynomials
        if ring is not self and ring is not None:
            raise NotImplementedError

        # Automatic choice of algorithm
        if algorithm is None:
            if deg > 100:
                algorithm = "sparse"
            else:
                algorithm = "dense"

        if algorithm != "dense" and algorithm != "sparse":
            raise ValueError("unknown algorithm '{}'".format(algorithm))

        # Check if the polynomial is a constant, in which case there are
        #   no roots. Note that the polynomial is not 0.
        if deg == 0:
            return []

        # The dense algorithm is to compute the roots from the factorization.
        if algorithm == "dense":
            cont = p.content_ideal().gen()
            if not cont.is_unit():
                p = p.map_coefficients(lambda c: c // cont)
            return p._roots_from_factorization(p.factor(), multiplicities)

        v = p.valuation()
        deg -= v
        cdef list roots

        # Root 0
        if v > 0:
            roots = [(self.zero(), v)] if multiplicities else [self.zero()]
            if deg == 0: # The shifted polynomial will be constant
                return roots
        else:
            roots = []

        p = p.shift(-v)
        cdef list e = p.exponents()
        cdef int i_min, i, j, k = len(e)

        # totally dense polynomial
        if k == 1 + deg:
            return roots + p._roots_from_factorization(p.factor(), multiplicities)

        cont = p.content_ideal().gen()
        if not cont.is_unit():
            p = p.map_coefficients(lambda c: c // cont)

        cdef list c = p.coefficients()

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(p.base_ring(), p.variable_name(), sparse=False)
        c_max_nbits = c[0].nbits()
        i_min = 0
        g = R.zero()

        # Looking for "large" gaps in the exponents
        # These gaps split the polynomial into lower degree components
        # Roots of modulus > 1 are common roots of the components
        for i in range(1, k):
            if e[i] - e[i-1] > c_max_nbits:
                g = g.gcd(R( {e[j] - e[i_min]: c[j] for j in range(i_min, i)} ))
                if g.is_one():
                    break
                i_min = i
                c_max_nbits = c[i].nbits()
            else:
                c_max_nbits = max(c[i].nbits(), c_max_nbits)

        # if no gap, directly return the roots of p
        if g.is_zero():
            roots.extend(p._roots_from_factorization(p.factor(), multiplicities))
            return roots

        g = g.gcd(R( {e[j] - e[i_min]: c[j] for j in range(i_min, k)} ))

        cdef list cc
        cdef list ee
        cdef int m1, m2
        cdef bint b1, b2
        # Computation of the roots of modulus 1, without multiplicities
        # 1 is root iff p(1) == 0 ; -1 iff p(-1) == 0
        if not multiplicities:
            if not sum(c):
                roots.append(self.one())
            s = 0
            for j in range(k):
                s += -c[j] if (e[j] % 2) else c[j]
            if not s:
                roots.append(-self.one())

        # Computation of the roots of modulus 1, with multiplicities
        # For the multiplicities, take the derivatives
        else:
            cc = c
            ee = e
            m1 = m2 = 0
            b1 = b2 = True

            for i in range(k):
                s1 = s2 = 0
                for j in range(k-i):
                    if b1: s1 += cc[j]
                    if b2: s2 += -cc[j] if (ee[j] % 2) else cc[j]
                if b1 and s1:
                    m1 = i
                    b1 = False
                if b2 and s2:
                    m2 = i
                    b2 = False
                # Stop asap
                if not (b1 or b2):
                    break

                # Sparse derivative, that is (p/x^v)' where v = p.val():
                ee = [ee[j] - ee[0] - 1 for j in range(1,k-i)]
                cc = [(ee[j] + 1) * cc[j+1] for j in range(k-i-1)]

            if m1 > 0:
                roots.append((self.one(), m1))
            if m2 > 0:
                roots.append((-self.one(), m2))

        # Add roots of modulus > 1 to `roots`:
        if multiplicities:
            roots.extend(r for r in g._roots_from_factorization(g.factor(), True)
                         if r[0].abs() > 1)
        else:
            roots.extend(r for r in g._roots_from_factorization(g.factor(), False)
                         if r.abs() > 1)

        return roots

    #################################
    #  Coercions to interfaces
    #################################

    def _gap_init_(self):
        """
        Return a GAP representation of ``self``.

        EXAMPLES::

            sage: gap(ZZ)  # indirect doctest                                           # needs sage.libs.gap
            Integers
        """
        return 'Integers'

    def _fricas_init_(self):
        """
        Return a FriCAS representation of ``self``.

        EXAMPLES::

            sage: fricas(ZZ)          # indirect doctest, optional - fricas
            Integer
        """
        return 'Integer'

    def _magma_init_(self, magma):
        """
        Return a magma representation of ``self``.

        EXAMPLES::

            sage: magma(ZZ)           # indirect doctest, optional - magma
            Integer Ring
        """
        return 'IntegerRing()'

    def _macaulay2_init_(self, macaulay2=None):
        """
        Return a macaulay2 representation of ``self``.

        EXAMPLES::

            sage: macaulay2(ZZ)       # optional - macaulay2
            ZZ
        """
        return "ZZ"

    def _polymake_init_(self):
        r"""
        Return the polymake representation of the integer ring.

        EXAMPLES::

            sage: polymake(ZZ)    # optional - jupymake # indirect doctest
            Integer
        """
        return '"Integer"'

    def _sympy_(self):
        r"""
        Return the SymPy set ``Integers``.

        EXAMPLES::

            sage: ZZ._sympy_()                                                          # needs sympy
            Integers
        """
        from sympy import Integers
        from sage.interfaces.sympy import sympy_init
        sympy_init()
        return Integers

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when
        evaluated.

        EXAMPLES::

            sage: sage_input(ZZ, verify=True)
            # Verified
            ZZ
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: ZZ._sage_input_(SageInputBuilder(), False)
            {atomic:ZZ}
        """
        return sib.name('ZZ')

    def valuation(self, p):
        r"""
        Return the discrete valuation with uniformizer ``p``.

        EXAMPLES::

            sage: v = ZZ.valuation(3); v                                                # needs sage.rings.padics
            3-adic valuation
            sage: v(3)                                                                  # needs sage.rings.padics
            1

        .. SEEALSO::

            :meth:`Order.valuation() <sage.rings.number_field.order.Order.valuation>`,
            :meth:`RationalField.valuation() <sage.rings.rational_field.RationalField.valuation>`
        """
        from sage.rings.padics.padic_valuation import pAdicValuation
        return pAdicValuation(self, p)

    def from_bytes(self, input_bytes, byteorder='big', is_signed=False):
        r"""
        Return the integer represented by the given array of bytes.

        Internally relies on the python ``int.from_bytes()`` method.

        INPUT:

        - ``input_bytes`` -- a bytes-like object or iterable producing bytes
        - ``byteorder`` -- string (default: ``'big'``); determines the byte order of
          ``input_bytes`` (can only be ``'big'`` or ``'little'``)
        - ``is_signed`` -- boolean (default: ``False``); determines whether to use two's
          compliment to represent the integer

        EXAMPLES::

            sage: ZZ.from_bytes(b'\x00\x10', byteorder='big')
            16
            sage: ZZ.from_bytes(b'\x00\x10', byteorder='little')
            4096
            sage: ZZ.from_bytes(b'\xfc\x00', byteorder='big', is_signed=True)
            -1024
            sage: ZZ.from_bytes(b'\xfc\x00', byteorder='big', is_signed=False)
            64512
            sage: ZZ.from_bytes([255, 0, 0], byteorder='big')
            16711680
            sage: type(_)
            <class 'sage.rings.integer.Integer'>
        """
        python_int = int.from_bytes(input_bytes, byteorder=byteorder, signed=is_signed)
        return self(python_int)

ZZ = IntegerRing_class()
Z = ZZ


def IntegerRing():
    """
    Return the integer ring.

    EXAMPLES::

        sage: IntegerRing()
        Integer Ring
        sage: ZZ==IntegerRing()
        True
    """
    return ZZ


def crt_basis(X, xgcd=None):
    r"""
    Compute and return a Chinese Remainder Theorem basis for the list ``X``
    of coprime integers.

    INPUT:

    - ``X`` -- list of Integers that are coprime in pairs

    - ``xgcd`` -- an optional parameter which is ignored

    OUTPUT:

    - ``E`` -- list of Integers such that ``E[i] = 1`` (mod ``X[i]``) and
      ``E[i] = 0`` (mod ``X[j]``) for all `j \neq i`

    For this explanation, let ``E[i]`` be denoted by `E_i`.

    The `E_i` have the property that if `A` is a list of objects, e.g.,
    integers, vectors, matrices, etc., where `A_i` is understood modulo
    `X_i`, then a CRT lift of `A` is simply

    .. MATH::

        \sum_i E_i A_i.

    ALGORITHM: To compute `E_i`, compute integers `s` and `t` such that

    .. MATH::

        s X_i + t \prod_{i \neq j} X_j = 1. (\*)

    Then

    .. MATH::

        E_i = t \prod_{i \neq j} X[j].

    Notice that equation
    (\*) implies that `E_i` is congruent to 1 modulo `X_i` and to 0
    modulo the other `X_j` for `j \neq i`.

    COMPLEXITY: We compute ``len(X)`` extended GCD's.

    EXAMPLES::

        sage: X = [11,20,31,51]
        sage: E = crt_basis([11,20,31,51])
        sage: E[0]%X[0], E[1]%X[0], E[2]%X[0], E[3]%X[0]
        (1, 0, 0, 0)
        sage: E[0]%X[1], E[1]%X[1], E[2]%X[1], E[3]%X[1]
        (0, 1, 0, 0)
        sage: E[0]%X[2], E[1]%X[2], E[2]%X[2], E[3]%X[2]
        (0, 0, 1, 0)
        sage: E[0]%X[3], E[1]%X[3], E[2]%X[3], E[3]%X[3]
        (0, 0, 0, 1)
    """
    if not isinstance(X, list):
        raise TypeError("X must be a list")
    if len(X) == 0:
        return []

    P = prod(X)

    Y = []
    # 2. Compute extended GCD's
    ONE = parent(X[0]).one()
    for i in range(len(X)):
        p = X[i]
        others = P // p
        g, _, t = p.xgcd(others)
        if g != ONE:
            raise ArithmeticError("the elements of the list X must be coprime in pairs")
        Y.append(t * others)
    return Y
