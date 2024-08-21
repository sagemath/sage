# sage_setup: distribution = sagemath-repl
r"""
Sage Input Formatting

This module provides the function :func:`sage_input` that takes an
arbitrary sage value and produces a sequence of commands that, if typed
at the ``sage:`` prompt, will recreate the value. If this is not
implemented for a particular value, then an exception is raised instead.
This might be useful in understanding a part of Sage, or for debugging.
For instance, if you have a value produced in a complicated way in the
middle of a debugging session, you could use :func:`sage_input` to find
a simple way to produce the same value. We attempt to produce commands
that are readable and idiomatic.::

    sage: sage_input(3)
    3
    sage: sage_input((polygen(RR) + RR(pi))^2, verify=True)                             # needs sage.symbolic
    # Verified
    R.<x> = RR[]
    x^2 + 6.2831853071795862*x + 9.869604401089358

With ``verify=True``, :func:`sage_input` also verifies the results, by
calling :func:`~sage.misc.sage_eval.sage_eval` on the result and
verifying that it is equal to the input.::

    sage: sage_input(GF(2)(1), verify=True)
    # Verified
    GF(2)(1)

We can generate code that works without the preparser, with
``preparse=False``; or we can generate code that will work whether or
not the preparser is enabled, with ``preparse=None``.  Generating code
with ``preparse=False`` may be useful to see how to create a certain
value in a Python or Cython source file.::

    sage: sage_input(5, verify=True)
    # Verified
    5
    sage: sage_input(5, preparse=False)
    ZZ(5)
    sage: sage_input(5, preparse=None)
    ZZ(5)
    sage: sage_input(5r, verify=True)
    # Verified
    5r
    sage: sage_input(5r, preparse=False)
    5
    sage: sage_input(5r, preparse=None)
    int(5)

Adding :func:`sage_input` support to your own classes is
straightforward.  You need to add a :func:`_sage_input_` method which
returns a :class:`SageInputExpression` (henceforth abbreviated as SIE)
which will reconstruct this instance of your class.

A ``_sage_input_`` method takes two parameters, conventionally named
``sib`` and ``coerced``.  The first argument is a
:class:`SageInputBuilder`; it has methods to build SIEs.  The second
argument, ``coerced``, is a boolean.  This is only useful if your class
is a subclass of :class:`Element` (although it is always present).  If
``coerced`` is ``False``, then your method must generate an expression
which will evaluate to a value of the correct type with the correct
parent.  If ``coerced`` is ``True``, then your method may generate an
expression of a type that has a canonical coercion to your type; and if
``coerced`` is 2, then your method may generate an expression of a type
that has a conversion to your type.

Let's work through some examples.  We'll build a sequence of functions
that would be acceptable as ``_sage_input_`` methods for the
:class:`~sage.rings.rational.Rational` class.

Here's the first and simplest version.::

    sage: def qq_sage_input_v1(self, sib, coerced):
    ....:     return sib(self.numerator())/sib(self.denominator())

We see that given a :class:`SageInputBuilder` ``sib``, you can construct
a SIE for a value ``v`` simply with ``sib(v)``, and you can construct a
SIE for a quotient with the division operator.  Of course, the other
operators also work, and so do function calls, method calls, subscripts,
etc.

We'll test with the following code, which you don't need to understand.
(It produces a list of 8 results, showing the formatted versions of -5/7
and 3, with the preparser either enabled or disabled and either with or
without an automatic coercion to QQ.)::

    sage: from sage.misc.sage_input import SageInputBuilder
    sage: def test_qq_formatter(fmt):
    ....:     results = []
    ....:     for v in [-5/7, QQ(3)]:
    ....:         for pp in [False, True]:
    ....:             for coerced in [False, True]:
    ....:                 sib = SageInputBuilder(preparse=pp)
    ....:                 results.append(sib.result(fmt(v, sib, coerced)))
    ....:     return results

    sage: test_qq_formatter(qq_sage_input_v1)
    [-ZZ(5)/ZZ(7), -ZZ(5)/ZZ(7), -5/7, -5/7, ZZ(3)/ZZ(1), ZZ(3)/ZZ(1), 3/1, 3/1]

Let's try for some shorter, perhaps nicer-looking output.  We'll start
by getting rid of the ``ZZ`` in the denominators; even without the
preparser, ``-ZZ(5)/7 == -ZZ(5)/ZZ(7)``.::

    sage: def qq_sage_input_v2(self, sib, coerced):
    ....:     return sib(self.numerator())/sib.int(self.denominator())

The ``int`` method on :class:`SageInputBuilder` returns a SIE for an
integer that is always represented in the simple way, without coercions.
(So, depending on the preparser mode, it might read in as an
:class:`~sage.rings.integer.Integer` or an ``int``.)::

    sage: test_qq_formatter(qq_sage_input_v2)
    [-ZZ(5)/7, -ZZ(5)/7, -5/7, -5/7, ZZ(3)/1, ZZ(3)/1, 3/1, 3/1]

Next let us get rid of the divisions by 1.  These are more complicated,
since if we are not careful we will get results in `\ZZ` instead of `\QQ`::

    sage: def qq_sage_input_v3(self, sib, coerced):
    ....:     if self.denominator() == 1:
    ....:         if coerced:
    ....:             return sib.int(self.numerator())
    ....:         else:
    ....:             return sib.name('QQ')(sib.int(self.numerator()))
    ....:     return sib(self.numerator())/sib.int(self.denominator())

We see that the \method{name} method gives an SIE representing a \sage
constant or function.::

    sage: test_qq_formatter(qq_sage_input_v3)
    [-ZZ(5)/7, -ZZ(5)/7, -5/7, -5/7, QQ(3), 3, QQ(3), 3]

This is the prettiest output we're going to get, but let's make one
further refinement.  Other :class:`_sage_input_` methods, like the one
for polynomials, analyze the structure of SIEs; they work better (give
prettier output) if negations are at the outside.  If the above code
were used for rationals, then ``sage_input(polygen(QQ) - 2/3)`` would
produce ``x + (-2/3)``; if we change to the following code, then we
would get ``x - 2/3`` instead.::

    sage: def qq_sage_input_v4(self, sib, coerced):
    ....:     num = self.numerator()
    ....:     neg = (num < 0)
    ....:     if neg: num = -num
    ....:     if self.denominator() == 1:
    ....:         if coerced:
    ....:             v = sib.int(num)
    ....:         else:
    ....:             v = sib.name('QQ')(sib.int(num))
    ....:     else:
    ....:         v = sib(num)/sib.int(self.denominator())
    ....:     if neg: v = -v
    ....:     return v

    sage: test_qq_formatter(qq_sage_input_v4)
    [-ZZ(5)/7, -ZZ(5)/7, -5/7, -5/7, QQ(3), 3, QQ(3), 3]

AUTHORS:

- Carl Witty (2008-04): new file

- Vincent Delecroix (2015-02): documentation formatting
"""

# ****************************************************************************
#       Copyright (C) 2008 Carl Witty <Carl.Witty@gmail.com>
#                     2015 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.real_mpfi', 'RealIntervalFieldElement')
lazy_import('sage.rings.complex_interval', 'ComplexIntervalFieldElement')


def sage_input(x, preparse=True, verify=False, allow_locals=False):
    r"""
    Return a sequence of commands that can be used to rebuild the object ``x``.

    INPUT:

    - ``x`` -- the value we want to find an input form for

    - ``preparse`` -- (default: ``True``) whether to generate code that requires
      the preparser.  With ``True``, generated code requires the preparser.
      With ``False``, generated code requires that the preparser not be used.
      With ``None``, generated code will work whether or not the preparser is
      used.

    - ``verify`` -- (default: ``False``) if ``True``, then the answer will be
      evaluated with :func:`sage_eval`, and an exception will be raised if the
      result is not equal to the original value.  (In fact, for ``verify=True``,
      :func:`sage_input` is effectively run three times, with ``preparse`` set
      to ``True``, ``False``, and ``None``, and all three results are checked.)
      This is particularly useful for doctests.

    - ``allow_locals`` -- (default: ``False``) if ``True``, then values that
      :func:`sage_input` cannot handle are returned in a dictionary, and the
      returned code assumes that this dictionary is passed as the ``locals``
      parameter of :func:`sage_eval`.  (Otherwise, if :func:`sage_input` cannot
      handle a value, an exception is raised.)

    EXAMPLES::

        sage: sage_input(GF(2)(1))
        GF(2)(1)
        sage: sage_input((GF(2)(0), GF(2)(1)), verify=True)
        # Verified
        GF_2 = GF(2)
        (GF_2(0), GF_2(1))

    When the preparser is enabled, we use the \sage generator syntax.::

        sage: K.<x> = GF(5)[]
        sage: sage_input(x^3 + 2*x, verify=True)
        # Verified
        R.<x> = GF(5)[]
        x^3 + 2*x
        sage: sage_input(x^3 + 2*x, preparse=False)
        R = GF(5)['x']
        x = R.gen()
        x**3 + 2*x

    The result of :func:`sage_input` is actually a pair of strings with a
    special ``__repr__`` method to print nicely.::

        sage: # needs sage.rings.real_mpfr sage.symbolic
        sage: r = sage_input(RealField(20)(pi), verify=True)
        sage: r
        # Verified
        RealField(20)(3.1415939)
        sage: isinstance(r, tuple)
        True
        sage: len(r)
        2
        sage: tuple(r)
        ('# Verified\n', 'RealField(20)(3.1415939)')

    We cannot find an input form for a function.::

        sage: sage_input((3, lambda x: x))
        Traceback (most recent call last):
        ...
        ValueError: cannot convert <function <lambda> at 0x...> to sage_input form

    But we can have :func:`sage_input` continue anyway, and return an input form
    for the rest of the expression, with ``allow_locals=True``.::

        sage: r = sage_input((3, lambda x: x), verify=True, allow_locals=True)
        sage: r
        LOCALS:
            _sil1: <function <lambda> at 0x...>
        # Verified
        (3, _sil1)
        sage: tuple(r)
        ('# Verified\n', '(3, _sil1)', {'_sil1': <function <lambda> at 0x...>})
    """
    if not verify:
        sib = SageInputBuilder(allow_locals=allow_locals, preparse=preparse)
        return sib.result(sib(x))

    # In verify mode, we actually compute and verify the answer with
    # all three settings of preparse.

    for pp in (True, False, None):
        sib = SageInputBuilder(allow_locals=allow_locals, preparse=pp)
        ans = sib.result(sib(x))
        verify_si_answer(x, ans, pp)
        if pp == preparse:
            ans_l = list(ans)
            ans_l[0] = '# Verified\n' + ans_l[0]
            final_answer = SageInputAnswer(*ans_l)

    return final_answer


class SageInputBuilder:
    r"""
    An instance of this class is passed to ``_sage_input_`` methods.
    It keeps track of the current state of the ``_sage_input_`` process,
    and contains many utility methods for building :class:`SageInputExpression`
    objects.

    In normal use, instances of :class:`SageInputBuilder` are created
    internally by :func:`sage_input`, but it may be useful to create
    an instance directly for testing or doctesting.

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder

    We can create a :class:`SageInputBuilder`, use it to create some
    :class:`SageInputExpression` s, and get a result.  (As mentioned
    above, this is only useful for testing or doctesting; normally
    you would just use :func:`sage_input`.)::

        sage: sib = SageInputBuilder()
        sage: sib.result((sib(3) + sib(4)) * (sib(5) + sib(6)))
        (3 + 4)*(5 + 6)
    """

    def __init__(self, allow_locals=False, preparse=True):
        r"""
        Initialize an instance of :class:`SageInputBuilder`.

        In normal use, instances of :class:`SageInputBuilder` are created
        internally by :func:`sage_input`, but it may be useful to create
        an instance directly for testing or doctesting.

        INPUT:

        - ``allow_locals`` -- (default: ``False``) if true, then values
          that cannot be converted to input form will be stored in
          a dictionary, which must be passed as the ``locals``
          when evaluating the result.

        - ``preparse`` -- (default: ``True``) if true, then the result
          will assume that the preparser is enabled.  If false, then
          the result will assume that the preparser is disabled.
          If ``None``, then the result will work whether or
          not the preparser is enabled.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: SageInputBuilder().preparse()
            True
            sage: SageInputBuilder(preparse=False).preparse()
            False
        """
        self._allow_locals = allow_locals
        self._preparse = preparse
        self._cached_types = set()
        self._cache = {}
        self._id_cache = {}
        self._parent_gens = {}
        self._next_local = 1
        self._locals = {}

    def __call__(self, x, coerced=False):
        r"""
        Try to convert an arbitrary value ``x`` into a
        :class:`SageInputExpression` (an SIE).

        We first check to see if an SIE has been cached for ``x``;
        if so, we return it.  If ``x`` is already an SIE, we return
        it unchanged.

        If ``x`` has a \method{_sage_input_} method, we call that
        method.

        Otherwise, if ``x`` is a value of some Python type that
        we know how to deal with, we convert it directly.

        Finally, for values we don't know how to convert, if
        ``self._allow_locals`` is true, we add it to a
        ``locals`` dictionary.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib(sib(3)))
            3

            sage: sib = SageInputBuilder()
            sage: sib.result(sib(GF(17)(5)))
            GF(17)(5)

        The argument ``coerced=True`` or ``coerced=2`` will get
        passed to the \method{_sage_input_} method of the argument.::

            sage: sib = SageInputBuilder()
            sage: sib.result(sib(GF(17)(5), True))
            5
            sage: sib.result(sib(RealField(200)(1.5), True))                            # needs sage.rings.real_mpfr
            1.5000000000000000000000000000000000000000000000000000000000000
            sage: sib.result(sib(RealField(200)(1.5), 2))                               # needs sage.rings.real_mpfr
            1.5

        Since :func:`sage_input` directly calls this method, all
        of the following are indirect doctests.::

            sage: sage_input(True)
            True
            sage: sage_input(-5r, verify=True)
            # Verified
            -5r
            sage: sage_input(7r, preparse=False, verify=True)
            # Verified
            7
            sage: sage_input(-11r, preparse=None, verify=True)
            # Verified
            -int(11)
            sage: sage_input(float(-infinity), preparse=True, verify=True)
            # Verified
            -float(infinity)
            sage: sage_input(float(NaN), preparse=True, verify=True)                    # needs sage.symbolic
            # Verified
            float(NaN)
            sage: sage_input(float(-pi), preparse=True, verify=True)                    # needs sage.symbolic
            # Verified
            float(-RR(3.1415926535897931))
            sage: sage_input(float(42), preparse=True, verify=True)                     # needs sage.rings.real_mpfr
            # Verified
            float(42)
            sage: sage_input("Hello, world\n", verify=True)
            # Verified
            'Hello, world\n'
            sage: sage_input("'", verify=True)
            # Verified
            "'"
            sage: sage_input('"', verify=True)
            # Verified
            '"'
            sage: sage_input(''' "'Hi,' she said." ''', verify=True)
            # Verified
            ' "\'Hi,\' she said." '
            sage: sage_input('Icky chars: \0\n\t\b\'\"\200\300\234', verify=True)
            # Verified
            'Icky chars: \x00\n\t\x08\'"\x80\xc0\x9c'
            sage: sage_input('unicode with spectral: \u1234\U00012345', verify=True)
            # Verified
            'unicode with spectral: \u1234\U00012345'
            sage: sage_input((2, 3.5, 'Hi'), verify=True)                               # needs sage.rings.real_mpfr
            # Verified
            (2, 3.5, 'Hi')
            sage: sage_input(lambda x: x)
            Traceback (most recent call last):
            ...
            ValueError: cannot convert <function <lambda> at 0x...> to sage_input form
            sage: sage_input(lambda x: x, allow_locals=True, verify=True)
            LOCALS:
              _sil1: <function <lambda> at 0x...>
            # Verified
            _sil1
        """
        # We want to look up x in our cache, to see if we've seen it before.
        # However, we don't want to assume that hashing x is always
        # efficient, so we only try the lookup if some value of the same
        # type as x has been cached.
        from sage.structure.all import parent

        if type(x) in self._cached_types:
            v = self._cache.get((parent(x), x))
            if v is not None:
                return v

        v = self._id_cache.get(id(x))
        if v is not None:
            return v[1]

        if isinstance(x, SageInputExpression):
            return x

        if hasattr(x, '_sage_input_'):
            return x._sage_input_(self, coerced)

        if x is None:
            return SIE_literal_stringrep(self, 'None')

        if isinstance(x, bool):
            return SIE_literal_stringrep(self, str(x))

        if isinstance(x, int):
            if self._preparse is True:
                if x < 0:
                    return -SIE_literal_stringrep(self, str(-x) + 'r')
                else:
                    return SIE_literal_stringrep(self, str(x) + 'r')
            elif self._preparse is False:
                return self.int(x)
            else:
                if x < 0:
                    return -self.name('int')(self.int(-x))
                else:
                    return self.name('int')(self.int(x))

        if isinstance(x, float):
            # floats could often have prettier output,
            # but I think they're rare enough in Sage that it's not
            # worth the effort.
            from math import inf
            if x == inf:
                return self.name('float')(self.name('infinity'))
            if x != x:
                return self.name('float')(self.name('NaN'))
            if x == -inf:
                return -self.name('float')(self.name('infinity'))
            if self._preparse is False and float(str(x)) == x:
                if x < 0:
                    return -SIE_literal_stringrep(self, str(-x))
                else:
                    return SIE_literal_stringrep(self, str(x))
            from sage.rings.real_mpfr import RR
            from sage.rings.integer_ring import ZZ
            rrx = RR(x)
            if rrx in ZZ and abs(rrx) < (1 << 53):
                return self.name('float')(self.int(ZZ(rrx)))
            return self.name('float')(RR(x))

        if isinstance(x, str):
            return SIE_literal_stringrep(self, repr(x))

        if isinstance(x, tuple):
            return SIE_tuple(self, [self(_) for _ in x], False)

        if isinstance(x, list):
            return SIE_tuple(self, [self(_) for _ in x], True)

        if isinstance(x, dict):
            return self.dict(x)

        if self._allow_locals:
            loc = self._next_local
            self._next_local += 1
            loc_name = '_sil%d' % loc
            self._locals[loc_name] = x
            return SIE_literal_stringrep(self, loc_name)
        else:
            raise ValueError("cannot convert {} to sage_input form".format(x))

    def preparse(self):
        r"""
        Check the preparse status.

        It returns ``True`` if the preparser will be enabled, ``False`` if it
        will be disabled, and ``None`` if the result must work whether or not
        the preparser is enabled.

        For example, this is useful in the \method{_sage_input_}
        methods of :class:`~sage.rings.integer.Integer` and :class:`RealNumber`; but most
        \method{_sage_input_} methods will not need to examine this.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: SageInputBuilder().preparse()
            True
            sage: SageInputBuilder(preparse=False).preparse()
            False
        """
        return self._preparse

    def int(self, n):
        r"""
        Return a raw SIE from the integer ``n``.

        As it is raw, it may read back as a Sage Integer or a Python int,
        depending on its size and whether the preparser is enabled.

        INPUT:

        - ``n`` -- a Sage Integer or a Python int

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.int(-3^50))
            -717897987691852588770249

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.int(-42r))
            -42
        """
        if n < 0:
            return -SIE_literal_stringrep(self, -n)
        else:
            return SIE_literal_stringrep(self, n)

    def float_str(self, n):
        r"""
        Given a string representing a floating-point number,
        produces a :class:`SageInputExpression` that formats as that
        string.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.float_str(repr(RR(e))))                                # needs sage.symbolic
            2.71828182845905
        """
        return SIE_literal_stringrep(self, n)

    def name(self, n):
        r"""
        Given a string representing a Python name,
        produces a :class:`SageInputExpression` for that name.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.name('pi') + sib.name('e'))
            pi + e
        """
        return SIE_literal_stringrep(self, n)

    def cache(self, x, sie, name):
        r"""
        INPUT:

        - ``x`` -- an arbitrary value

        - ``sie`` -- a :class:`SageInputExpression`

        - ``name`` -- a requested variable name

        Enters ``x`` and ``sie`` in a cache, so that subsequent calls
        ``self(x)`` will directly return ``sie``.  Also, marks the
        requested name of this ``sie`` to be ``name``.

        This should almost always be called as part of the
        \method{_sage_input_} method of a parent.  It may also be called
        on values of an arbitrary type, which may be useful if the values
        are both large and likely to be used multiple times in a single
        expression.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie42 = sib(GF(101)(42))
            sage: sib.cache(GF(101)(42), sie42, 'the_ultimate_answer')
            sage: sib.result(sib(GF(101)(42)) + sib(GF(101)(42)))
            the_ultimate_answer = GF(101)(42)
            the_ultimate_answer + the_ultimate_answer

        Note that we don't assign the result to a variable if the value
        is only used once.::

            sage: sib = SageInputBuilder()
            sage: sie42 = sib(GF(101)(42))
            sage: sib.cache(GF(101)(42), sie42, 'the_ultimate_answer')
            sage: sib.result(sib(GF(101)(42)) + sib(GF(101)(43)))
            GF_101 = GF(101)
            GF_101(42) + GF_101(43)
        """
        from sage.structure.all import parent

        self._cached_types.add(type(x))
        self._cache[(parent(x), x)] = sie
        sie._sie_preferred_varname = name

    def id_cache(self, x, sie, name):
        r"""
        INPUT:

        - ``x`` -- an arbitrary value

        - ``sie`` -- a :class:`SageInputExpression`

        - ``name`` -- a requested variable name

        Enters ``x`` and ``sie`` in a cache, so that subsequent calls
        ``self(x)`` will directly return ``sie``.  Also, marks the
        requested name of this ``sie`` to be ``name``.  Differs from
        the \method{cache} method in that the cache is keyed by
        ``id(x)`` instead of by ``x``.

        This may be called on values of an arbitrary type, which may
        be useful if the values are both large and likely to be used
        multiple times in a single expression; it should be preferred to
        \method{cache} if equality on the values is difficult or impossible
        to compute.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: x = polygen(ZZ)
            sage: sib = SageInputBuilder()
            sage: my_42 = 42*x
            sage: sie42 = sib(my_42)
            sage: sib.id_cache(my_42, sie42, 'the_ultimate_answer')
            sage: sib.result(sib(my_42) + sib(my_42))
            R.<x> = ZZ[]
            the_ultimate_answer = 42*x
            the_ultimate_answer + the_ultimate_answer

        Since id_cache keys off of object identity ("is"), the
        following does not trigger the cache.::

            sage: sib.result(sib(42*x) + sib(42*x))
            42*x + 42*x

        Note that we don't assign the result to a variable if the value
        is only used once.::

            sage: sib = SageInputBuilder()
            sage: my_42 = 42*x
            sage: sie42 = sib(my_42)
            sage: sib.id_cache(my_42, sie42, 'the_ultimate_answer')
            sage: sib.result(sib(my_42) + sib(43*x))
            R.<x> = ZZ[]
            42*x + 43*x
        """
        # If we just mapped id(x) -> sie, then it's possible that x could
        # be freed and another value allocated at the same position,
        # corrupting the cache.  But since we store x, that can't happen;
        # we don't even have to look at x when we read the cache.
        self._id_cache[id(x)] = (x, sie)
        sie._sie_preferred_varname = name

    def import_name(self, module, name, alt_name=None):
        r"""
        INPUT:

        - ``module``, ``name``, ``alt_name`` -- strings

        Creates an expression that will import a name from a module and
        then use that name.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: v1 = sib.import_name('sage.foo.bar', 'baz')
            sage: v2 = sib.import_name('sage.foo.bar', 'ZZ', 'not_the_real_ZZ')
            sage: sib.result(v1+v2)
            from sage.foo.bar import baz
            from sage.foo.bar import ZZ as not_the_real_ZZ
            baz + not_the_real_ZZ

        We adjust the names if there is a conflict.::

            sage: sib = SageInputBuilder()
            sage: v1 = sib.import_name('sage.foo', 'poly')
            sage: v2 = sib.import_name('sage.bar', 'poly')
            sage: sib.result(v1+v2)
            from sage.foo import poly as poly1
            from sage.bar import poly as poly2
            poly1 + poly2
        """
        return SIE_import_name(self, module, name, alt_name)

    def assign(self, e, val):
        r"""
        Construct a command that performs the assignment ``e=val``.

        Can only be used as an argument to the ``command`` method.

        INPUT:

        - ``e``, ``val`` -- SageInputExpression

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: circular = sib([None])
            sage: sib.command(circular, sib.assign(circular[0], circular))
            sage: sib.result(circular)
            si = [None]
            si[0] = si
            si
        """
        e = self(e)
        val = self(val)

        return SIE_assign(self, e, val)

    def command(self, v, cmd):
        r"""
        INPUT:

        - ``v``, ``cmd`` -- SageInputExpression

        Attaches a command to v, which will be executed before v is used.
        Multiple commands will be executed in the order added.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: incr_list = sib([])
            sage: sib.command(incr_list, incr_list.append(1))
            sage: sib.command(incr_list, incr_list.extend([2, 3]))
            sage: sib.result(incr_list)
            si = []
            si.append(1)
            si.extend([2, 3])
            si
        """
        v = self(v)
        cmd = self(cmd)

        v._sie_commands.append(cmd)

    def dict(self, entries):
        r"""
        Given a dictionary, or a list of (key, value) pairs,
        produces a :class:`SageInputExpression` representing
        the dictionary.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.dict({1:1, 2:5/2, 3:100/3}))
            {1:1, 2:5/2, 3:100/3}
            sage: sib.result(sib.dict([('hello', 'sunshine'), ('goodbye', 'rain')]))
            {'hello':'sunshine', 'goodbye':'rain'}
        """
        if isinstance(entries, dict):
            entries = list(entries.items())
        entries = [(self(key), self(val)) for (key, val) in entries]
        return SIE_dict(self, entries)

    def getattr(self, sie, attr):
        r"""
        Given a :class:`SageInputExpression` representing ``foo``
        and an attribute name bar, produce a :class:`SageInputExpression`
        representing ``foo.bar``.  Normally, you could just use
        attribute-access syntax, but that doesn't work if bar
        is some attribute that bypasses __getattr__ (such as if
        bar is '__getattr__' itself).

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.getattr(ZZ, '__getattr__')
            {getattr: {atomic:ZZ}.__getattr__}
            sage: sib.getattr(sib.name('foo'), '__new__')
            {getattr: {atomic:foo}.__new__}
        """
        return SIE_getattr(self, self(sie), attr)

    def empty_subscript(self, parent):
        r"""
        Given a :class:`SageInputExpression` representing ``foo``,
        produces a :class:`SageInputExpression` representing ``foo[]``.
        Since this is not legal Python syntax, it is useful only for
        producing the \sage generator syntax for a polynomial ring.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.empty_subscript(sib(2) + sib(3)))
            (2 + 3)[]

        The following calls this method indirectly.::

            sage: sage_input(polygen(ZZ['y']))
            R.<x> = ZZ['y'][]
            x
        """
        return SIE_subscript(self, parent, None)

    def use_variable(self, sie, name):
        r"""
        Marks the :class:`SageInputExpression` ``sie`` to use a variable
        even if it is only referenced once.  (If ``sie`` is the final
        top-level expression, though, it will not use a variable.)

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: e = sib.name('MatrixSpace')(ZZ, 10, 10)
            sage: sib.use_variable(e, 'MS')
            sage: sib.result(e.zero_matrix())
            MS = MatrixSpace(ZZ, 10, 10)
            MS.zero_matrix()

        Without the call to use_variable, we get this instead::

            sage: sib = SageInputBuilder()
            sage: e = sib.name('MatrixSpace')(ZZ, 10, 10)
            sage: sib.result(e.zero_matrix())
            MatrixSpace(ZZ, 10, 10).zero_matrix()

        And even with the call to use_variable, we don't use a variable here::

            sage: sib = SageInputBuilder()
            sage: e = sib.name('MatrixSpace')(ZZ, 10, 10)
            sage: sib.use_variable(e, 'MS')
            sage: sib.result(e)
            MatrixSpace(ZZ, 10, 10)
        """
        sie._sie_preferred_varname = name
        sie._sie_request_use_var = True

    def share(self, sie):
        r"""
        Mark the given expression as sharable, so that it will be replaced
        by a variable if it occurs multiple times in the expression.
        (Most non-single-token expressions are already sharable.)

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

        Without explicitly using .share(), string literals are not shared::

            sage: sib = SageInputBuilder()
            sage: e = sib('hello')
            sage: sib.result(sib((e, e)))
            ('hello', 'hello')

        See the difference if we use .share()::

            sage: sib = SageInputBuilder()
            sage: e = sib('hello')
            sage: sib.share(e)
            sage: sib.result(sib((e, e)))
            si = 'hello'
            (si, si)
        """
        sie._sie_share = True

    def parent_with_gens(self, parent, sie, gen_names, name, gens_syntax=None):
        r"""
        This method is used for parents with generators, to manage the
        \sage preparser generator syntax (like ``K.<x> = QQ[]``).

        The \method{_sage_input_} method of a parent class with
        generators should construct a :class:`SageInputExpression` for
        the parent, and then call this method with the parent itself,
        the constructed SIE, a sequence containing the names of the
        generators, and (optionally) another SIE to use if the \sage
        generator syntax is used; typically this will be the same as
        the first SIE except omitting a ``names`` parameter.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder


            sage: def test_setup(use_gens=True, preparse=True):
            ....:     sib = SageInputBuilder(preparse=preparse)
            ....:     gen_names=('foo', 'bar')
            ....:     parent = "some parent"
            ....:     normal_sie = sib.name('make_a_parent')(names=gen_names)
            ....:     if use_gens:
            ....:         gens_sie = sib.name('make_a_parent')()
            ....:     else:
            ....:         gens_sie = None
            ....:     name = 'the_thing'
            ....:     result = sib.parent_with_gens(parent, normal_sie,
            ....:                                   gen_names, name,
            ....:                                   gens_syntax=gens_sie)
            ....:     return sib, result

            sage: sib, par_sie = test_setup()
            sage: sib.result(par_sie)
            make_a_parent(names=('foo', 'bar'))

            sage: sib, par_sie = test_setup()
            sage: sib.result(sib(3) * sib.gen("some parent", 0))
            the_thing.<foo,bar> = make_a_parent()
            3*foo

            sage: sib, par_sie = test_setup(preparse=False)
            sage: sib.result(par_sie)
            make_a_parent(names=('foo', 'bar'))

            sage: sib, par_sie = test_setup(preparse=False)
            sage: sib.result(sib(3) * sib.gen("some parent", 0))
            the_thing = make_a_parent(names=('foo', 'bar'))
            foo,bar = the_thing.gens()
            ZZ(3)*foo

            sage: sib, par_sie = test_setup(use_gens=False)
            sage: sib.result(par_sie)
            make_a_parent(names=('foo', 'bar'))

            sage: sib, par_sie = test_setup(use_gens=False)
            sage: sib.result(sib(3) * sib.gen("some parent", 0))
            the_thing = make_a_parent(names=('foo', 'bar'))
            foo,bar = the_thing.gens()
            3*foo

            sage: sib, par_sie = test_setup()
            sage: sib.result(par_sie - sib.gen("some parent", 1))
            the_thing.<foo,bar> = make_a_parent()
            the_thing - bar
        """
        v = SIE_gens_constructor(self, sie, gen_names, gens_syntax=gens_syntax)
        self.cache(parent, v, name)
        gens = [SIE_gen(self, v, n) for n in gen_names]
        self._parent_gens[parent] = gens
        v._sie_gens = gens
        return v

    def gen(self, parent, n=0):
        r"""
        Given a parent, returns a :class:`SageInputExpression` for
        the `n`-th (default: 0) generator of the parent.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.gen(ZZ['y']))
            R.<y> = ZZ[]
            y
        """
        if parent not in self._parent_gens:
            self(parent)
            if parent not in self._parent_gens:
                raise ValueError("{} did not register generators for sage_input".format(parent))

        gens = self._parent_gens[parent]

        if n > len(gens):
            raise ValueError("{} registered only {} generators for sage_input".format(parent, len(gens)))

        return gens[n]

    def prod(self, factors, simplify=False):
        r"""
        Given a sequence, returns a :class:`SageInputExpression`
        for the product of the elements.

        With ``simplify=True``, performs some simplifications
        first.  If any element is formatted as a string ``'0'``,
        then that element is returned directly.  If any element is
        formatted as a string ``'1'``, then it is removed
        from the sequence (unless it is the only element in the sequence).
        And any negations are removed from the elements and moved to the
        outside of the product.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.prod([-1, 0, 1, -2]))
            -1*0*1*-2

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.prod([-1, 0, 1, 2], simplify=True))
            0

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.prod([-1, 2, -3, -4], simplify=True))
            -2*3*4

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.prod([-1, 1, -1, -1], simplify=True))
            -1

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.prod([1, 1, 1], simplify=True))
            1
        """
        neg = False
        factors = [self(factor) for factor in factors]
        if simplify:
            i = 0
            while i < len(factors):
                factor = factors[i]
                while isinstance(factor, SIE_unary) and factor._sie_op == '-':
                    neg = not neg
                    factor = factor._sie_operand
                    factors[i] = factor
                if isinstance(factor, SIE_literal_stringrep) and factor._sie_value == '0':
                    factors = [factor]
                    neg = False
                    break
                if isinstance(factor, SIE_literal_stringrep) and factor._sie_value == '1':
                    factors[i:i + 1] = []
                else:
                    i += 1
            if len(factors) == 0:
                factors.append(SIE_literal_stringrep(self, '1'))

        prod = factors[0]
        for factor in factors[1:]:
            prod = prod * factor
        if neg:
            prod = -prod
        return prod

    def sum(self, terms, simplify=False):
        r"""
        Given a sequence, returns a :class:`SageInputExpression`
        for the product of the elements.

        With ``simplify=True``, performs some simplifications
        first.  If any element is formatted as a string ``'0'``,
        then it is removed from the sequence (unless it is the only
        element in the sequence); and any instances of ``a + -b``
        are changed to ``a - b``.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.sum([-1, 0, 1, 0, -1]))
            -1 + 0 + 1 + 0 + -1

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.sum([-1, 0, 1, 0, -1], simplify=True))
            -1 + 1 - 1

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.sum([0, 0, 0], simplify=True))
            0
        """
        terms = [self(term) for term in terms]
        if simplify:
            i = 0
            while i < len(terms):
                term = terms[i]
                if isinstance(term, SIE_literal_stringrep) and term._sie_value == '0':
                    terms[i:i + 1] = []
                else:
                    i += 1
            if len(terms) == 0:
                terms.append(SIE_literal_stringrep(self, '0'))

        sum = terms[0]
        for term in terms[1:]:
            negate = False
            while simplify and isinstance(term, SIE_unary) and term._sie_op == '-':
                negate = not negate
                term = term._sie_operand
            if negate:
                sum = sum - term
            else:
                sum = sum + term
        return sum

    def result(self, e):
        r"""
        Given a :class:`SageInputExpression` constructed using ``self``,
        returns a tuple of a list of commands and an expression
        (and possibly a dictionary of local variables) suitable for
        :func:`sage_eval`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: r = sib.result(sib(6) * sib(7)); r
            6*7
            sage: tuple(r)
            ('', '6*7')
        """
        sif = SageInputFormatter()

        # Even if use_variable was called on e, don't automatically
        # use a variable for it.
        e._sie_request_use_var = False

        e._sie_prepare(sif)

        s = sif.format(e, 0)

        locals = self._locals
        if len(locals):
            return SageInputAnswer(sif._commands, sif.format(e, 0), locals)
        else:
            return SageInputAnswer(sif._commands, sif.format(e, 0))


# Python's precedence levels.  Hand-transcribed from section 5.14 of
# the Python 2 reference manual.  In the Python 3 reference manual
# this is section 6.16.
# See https://docs.python.org/3/reference/expressions.html
_prec_lambda = 2
_prec_or = 4
_prec_and = 6
_prec_not = 8
_prec_membership = 10
_prec_identity = 12
_prec_comparison = 14
_prec_bitor = 16
_prec_bitxor = 18
_prec_bitand = 20
_prec_shift = 22
_prec_addsub = 24
_prec_muldiv = 26
_prec_negate = 28
_prec_bitnot = 30
_prec_exponent = 32
_prec_attribute = 34
_prec_subscript = 36
_prec_slicing = 38
_prec_funcall = 40
_prec_atomic = 42


class SageInputExpression:
    r"""
    Subclasses of this class represent expressions for :func:`sage_input`.
    \sage classes should define a \method{_sage_input_} method, which
    will return an instance of :class:`SageInputExpression`, created using
    methods of :class:`SageInputBuilder`.

    To the extent possible, operations on :class:`SageInputExpression` objects
    construct a new :class:`SageInputExpression` representing that operation.
    That is, if ``a`` is a :class:`SageInputExpression`, then ``a + b``
    constructs a :class:`SageInputExpression` representing this sum.
    This also works for attribute access, function calls, subscripts, etc.
    Since arbitrary attribute accesses might be used to construct a new
    attribute-access expression, all internal attributes and methods
    have names that begin with ``_sie_`` to reduce the chance of
    collisions.

    It is expected that instances of this class will not be directly
    created outside this module; instead, instances will be created
    using methods of :class:`SageInputBuilder` and :class:`SageInputExpression`.

    Values of type :class:`SageInputExpression` print in a fairly ugly
    way, that reveals the internal structure of the expression tree.
    """

    def __init__(self, sib):
        r"""
        Initialize a :class:`SageInputExpression`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib(3) # indirect doctest
            sage: sie
            {atomic:3}
            sage: sie._sie_builder is sib
            True
        """
        self._sie_refcount = 0
        self._sie_builder = sib
        self._sie_context = None
        self._sie_preferred_varname = None
        self._sie_varname = None
        self._sie_request_use_var = False
        self._sie_use_var = False
        self._sie_requested_varname = False
        self._sie_commands = []

    def _sie_is_simple(self):
        r"""
        Return ``True`` if this :class:`SageInputExpression` is simple
        enough that duplicate uses are not worth caching.  Normally
        this will be true if the expression represents a single token.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sib.name('QQ')._sie_is_simple()
            True
            sage: sib(GF(2))._sie_is_simple()
            False
        """
        return False

    def _sie_referenced(self):
        r"""
        Return a list of the immediate subexpressions of this
        :class:`SageInputExpression`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: len(sib(GF(2))._sie_referenced())
            2
            sage: sib(5)._sie_referenced()
            []
        """
        return []

    def _sie_prepare(self, sif):
        r"""
        We traverse the entire expression DAG to prepare for printing.
        Here, we notice nodes with more than one parent, and mark them
        to replace with a variable (rather than generating the value
        multiple times).

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter
            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: pair = sib((GF(2), GF(2)))
            sage: single = sib(GF(2))
            sage: single._sie_refcount
            0
            sage: single._sie_use_var
            False
            sage: sib((GF(2), GF(2)))._sie_prepare(sif)
            sage: single._sie_refcount
            2
            sage: single._sie_use_var
            True
        """
        if self._sie_context is not sif:
            self._sie_context = sif
            self._sie_refcount = 0
        self._sie_refcount += 1
        if self._sie_request_use_var:
            self._sie_require_varname(sif)
            self._sie_use_var = True
        if not self._sie_is_simple():
            if self._sie_refcount == 2:
                self._sie_require_varname(sif)
                self._sie_use_var = True
        if self._sie_refcount == 1:
            for r in self._sie_referenced():
                r._sie_prepare(sif)
            for r in self._sie_commands:
                r._sie_prepare(sif)

    def _sie_require_varname(self, sif):
        r"""
        Mark this :class:`SageInputExpression` as requiring a variable name,
        and register it with a :class:`SageInputFormatter` (which will
        allocate a variable name at the end of the preparatory phase).

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter
            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib(3)
            sage: sie._sie_require_varname(sif)
            sage: sie._sie_requested_varname
            True
        """
        if not self._sie_requested_varname:
            sif.register_name(self._sie_preferred_varname)
            self._sie_requested_varname = True
            self._sie_generated = False

    def _sie_get_varname(self, sif):
        r"""
        Get the variable name that the :class:`SageInputFormatter` allocated
        for this :class:`SageInputExpression`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter
            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib(3)
            sage: sie._sie_require_varname(sif)
            sage: sie._sie_get_varname(sif)
            'si'
        """
        if self._sie_varname is None:
            self._sie_varname = sif.get_name(self._sie_preferred_varname)

        return self._sie_varname

    def _sie_is_negation(self):
        r"""
        Test whether a :class:`SageInputExpression` is a negation.

        Despite the obscure name, this is intended to be a public method.

        See the documentation for \method{SIE_unary._sie_is_negation}
        for useful examples.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter
            sage: sib = SageInputBuilder()
            sage: sie = sib.name('foo')
            sage: sie._sie_is_negation()
            False
        """
        return False

    def __call__(self, *args, **kwargs):
        r"""
        Given a :class:`SageInputExpression`, build a new
        :class:`SageInputExpression` representing a function call node
        (with ``self`` as the function).

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie(4)
            {call: {atomic:3}({atomic:4})}
        """
        args = [self._sie_builder(_) for _ in args]
        for k in kwargs:
            kwargs[k] = self._sie_builder(kwargs[k])
        return SIE_call(self._sie_builder, self, args, kwargs)

    def __getitem__(self, key):
        r"""
        Given a :class:`SageInputExpression`, build a new
        :class:`SageInputExpression` representing a subscript expression
        (with ``self`` as the value being subscripted).

        Currently, slices are not supported.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie[4]
            {subscr: {atomic:3}[{atomic:4}]}
            sage: sie[sib.name('x'), sib.name('y')]
            {subscr: {atomic:3}[{tuple: ({atomic:x}, {atomic:y})}]}
        """
        skey = self._sie_builder(key)
        return SIE_subscript(self._sie_builder, self, skey)

    def __getattr__(self, attr):
        r"""
        Given a :class:`SageInputExpression`, build a new
        :class:`SageInputExpression` representing an attribute access.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib.name('x')
            sage: sie.foo
            {getattr: {atomic:x}.foo}
            sage: sie.foo()
            {call: {getattr: {atomic:x}.foo}()}
        """
        return SIE_getattr(self._sie_builder, self, attr)

    def _rich_repr_(self, display_manager, **kwds):
        """
        Disable rich output.

        This is necessary because otherwise our :meth:`__getattr__`
        would be called.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib.name('x')
            sage: sie._rich_repr_(get_display_manager()) is None
            True
        """
        return None

    def __pow__(self, other):
        r"""
        Compute an expression tree for ``self ** other``.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie ^ 4
            {binop:** {atomic:3} {atomic:4}}
        """
        return self._sie_binop('**', other)

    def __mul__(self, other):
        r"""
        Compute an expression tree for ``self * other``.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie * 4
            {binop:* {atomic:3} {atomic:4}}
        """
        return self._sie_binop('*', other)

    def __truediv__(self, other):
        r"""
        Compute an expression tree for ``self / other``.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie / 4
            {binop:/ {atomic:3} {atomic:4}}
        """
        return self._sie_binop('/', other)

    def __add__(self, other):
        r"""
        Compute an expression tree for ``self + other``.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie + 4
            {binop:+ {atomic:3} {atomic:4}}
        """
        return self._sie_binop('+', other)

    def __sub__(self, other):
        r"""
        Compute an expression tree for ``self - other``.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie - 4
            {binop:- {atomic:3} {atomic:4}}
        """
        return self._sie_binop('-', other)

    def _sie_binop(self, op, other):
        r"""
        Compute an expression tree for ``self OP other``,
        where OP is a string representing a binary operator (such as
        '+' or '**').

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: v = sib.name('x')._sie_binop('%', sib.name('y'))
            sage: type(v)
            <class 'sage.misc.sage_input.SIE_binary'>
            sage: (v)._sie_op
            '%'
            sage: v
            {binop:% {atomic:x} {atomic:y}}
        """
        return SIE_binary(self._sie_builder, op, self, self._sie_builder(other))

    def __neg__(self):
        r"""
        Compute an expression tree for ``-self``.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: -sie
            {unop:- {atomic:3}}
        """
        return self._sie_unop('-')

    def __pos__(self):
        r"""
        Compute an expression tree for ``+self``.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: +sie
            {unop:+ {atomic:3}}
        """
        return self._sie_unop('+')

    def __invert__(self):
        r"""
        Compute an expression tree for ``~self``.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: ~sie
            {unop:~ {atomic:3}}
        """
        return self._sie_unop('~')

    def __abs__(self):
        r"""
        Compute an expression tree for ``abs(self)``.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: abs(sie)
            {call: {atomic:abs}({atomic:3})}
        """
        return self._sie_builder.name('abs')(self)

    def _sie_unop(self, op):
        r"""
        Compute an expression tree for ``OP self``,
        where OP is a string representing a unary operator (such as
        '-' or '~').

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: v = sie._sie_unop('~')
            sage: type(v)
            <class 'sage.misc.sage_input.SIE_unary'>
            sage: (v)._sie_op
            '~'
            sage: v
            {unop:~ {atomic:3}}
        """
        return SIE_unary(self._sie_builder, op, self)

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this expression, and the
        precedence of the top-level operator in the expression.

        EXAMPLES:

        Actually, all of these are examples of the \method{_sie_format}
        method on subclasses of :class:`SageInputExpression`;
        :class:`SageInputExpression` itself is an abstract base class
        (that cannot be instantiated).::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter
            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib(3)

            sage: for v in (sie, sie+7, sie/5):
            ....:     v._sie_prepare(sif)
            ....:     v._sie_format(sif)
            ('3', 42)
            ('3 + 7', 24)
            ('3/5', 26)
            sage: v = sib.assign(sib.name('foo').x, 3)
            sage: v._sie_prepare(sif)
            sage: v._sie_format(sif)
            Traceback (most recent call last):
            ...
            ValueError: Cannot format SIE_assign as expression
        """
        raise NotImplementedError

    def _sie_format_statement(self, sif):
        r"""
        Return the formatted string value of this expression, when
        used as a statement.

        On most :class:`SageInputExpression`s, this forwards directly
        to the \method{_sie_format} method.  However, on
        :class:`SageInputExpression`s that actually represent
        statements (such as :class:`SIE_assign`), this method
        has an implementation and \method{_sie_format} raises
        an error.  (This is to prevent accidental use of
        :class:`SIE_assign` as a value.)

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter
            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: v = sib(3)
            sage: v._sie_prepare(sif)
            sage: v._sie_format_statement(sif)
            '3'
            sage: v = sib.assign(sib.name('foo').x, 3)
            sage: v._sie_prepare(sif)
            sage: v._sie_format_statement(sif)
            'foo.x = 3'
        """
        result, prec = self._sie_format(sif)
        return result


class SIE_literal(SageInputExpression):
    r"""
    An abstract base class for ``literals`` (basically, values which
    consist of a single token).

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder, SIE_literal

        sage: sib = SageInputBuilder()
        sage: sie = sib(3)
        sage: sie
        {atomic:3}
        sage: isinstance(sie, SIE_literal)
        True
    """

    def _sie_is_simple(self):
        r"""
        Report that :class:`SIE_literal` values are not worth replacing by
        variables (for ``common subexpression elimination``) even if they
        occur multiple times in an expression.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie._sie_is_simple()
            True

            sage: sib.share(sie)
            sage: sie._sie_is_simple()
            False
            sage: sie._sie_share
            True
        """
        # Perhaps this should actually look at the formatted length of self,
        # and sometimes return false?  If some 50-digit integer occurs multiple
        # times in an expression, it might be better to do the replacement.
        return not self._sie_share


class SIE_literal_stringrep(SIE_literal):
    r"""
    Values in this class are leaves in a :func:`sage_input` expression
    tree.  Typically they represent a single token, and consist of the
    string representation of that token.  They are used for integer,
    floating-point, and string literals, and for name expressions.

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder, SIE_literal_stringrep

        sage: sib = SageInputBuilder()
        sage: isinstance(sib(3), SIE_literal_stringrep)
        True
        sage: isinstance(sib(3.14159, True), SIE_literal_stringrep)                     # needs sage.rings.real_mpfr
        True
        sage: isinstance(sib.name('pi'), SIE_literal_stringrep)
        True
        sage: isinstance(sib(False), SIE_literal_stringrep)
        True
        sage: sib(False)
        {atomic:False}
    """

    def __init__(self, sib, n):
        r"""
        Initialize a :class:`SIE_literal_stringrep` value.

        INPUT:

        - ``sib`` -- a :class:`SageInputBuilder`

        - ``n`` -- string; the value to be printed for this expression

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib(3)
            {atomic:3}
            sage: sib(3)._sie_value
            '3'
        """
        super().__init__(sib)
        self._sie_value = str(n)
        self._sie_share = False

    def __repr__(self):
        r"""
        Return a string representing this :class:`SIE_literal_stringrep`
        value.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib(3)
            {atomic:3}
            sage: sib("\n")
            {atomic:'\n'}
        """
        return "{atomic:%s}" % self._sie_value

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this expression, and an indication
        that it is ``atomic`` (never needs to be parenthesized).

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib(True)
            sage: sie._sie_prepare(sif)
            sage: sie._sie_format(sif)
            ('True', 42)
        """
        return self._sie_value, _prec_atomic


class SIE_call(SageInputExpression):
    r"""
    This class represents a function-call node in a :func:`sage_input`
    expression tree.

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sie = sib.name('GF')
        sage: sie(49)
        {call: {atomic:GF}({atomic:49})}
    """

    def __init__(self, sib, func, args, kwargs):
        r"""
        Initialize an instance of :class:`SIE_call`.

        INPUT:

        - ``sib`` -- a :class:`SageInputBuilder`

        - ``func`` -- a :class:`SageInputExpression` representing a function

        - ``args`` -- list of instances of :class:`SageInputExpression`
          representing the positional arguments

        - ``kwargs`` -- dictionary mapping strings to instances of
          :class:`SageInputExpression` representing the keyword arguments

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib('RealField')(53, rnd='RNDZ')
        """
        super().__init__(sib)
        self._sie_func = func
        self._sie_args = args
        self._sie_kwargs = kwargs

    def __repr__(self):
        r"""
        Return a string representing this :class:`SIE_call` value.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib('RealField')(53, rnd='RNDZ')
        """
        func = repr(self._sie_func)
        args = [repr(arg) for arg in self._sie_args]
        kwargs = sorted(k + '=' + repr(v)
                        for k, v in self._sie_kwargs.items())
        all_args = ', '.join(args + kwargs)
        return "{call: %s(%s)}" % (func, all_args)

    def _sie_referenced(self):
        r"""
        Return a list of the immediate subexpressions of this :class:`SIE_call`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib('RealField')(53, rnd='RNDZ')
            sage: sie._sie_referenced()
            [{atomic:53}, {atomic:'RealField'}, {atomic:'RNDZ'}]
        """
        refs = self._sie_args[:]
        refs.append(self._sie_func)
        refs.extend(self._sie_kwargs.values())
        return refs

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this expression, and an indication
        that it is a function call.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.name('RealField')(53, rnd='RNDZ')
            sage: sie._sie_prepare(sif)
            sage: sie._sie_format(sif)
            ("RealField(53, rnd='RNDZ')", 40)
        """
        func = sif.format(self._sie_func, _prec_attribute)
        args = [sif.format(arg, 0) for arg in self._sie_args]
        kwargs = sorted(k + '=' + sif.format(v, 0)
                        for k, v in self._sie_kwargs.items())
        all_args = ', '.join(args + kwargs)
        return ('%s(%s)' % (func, all_args), _prec_funcall)


class SIE_subscript(SageInputExpression):
    r"""
    This class represents a subscript node in a :func:`sage_input`
    expression tree.

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sie = sib.name('QQ')['x,y']
        sage: sie
        {subscr: {atomic:QQ}[{atomic:'x,y'}]}
    """

    def __init__(self, sib, coll, key):
        r"""
        Initialize an instance of :class:`SIE_subscript`.

        INPUT:

        - ``sib`` -- a :class:`SageInputBuilder`

        - ``coll`` -- a :class:`SageInputExpression` representing a collection

        - ``key`` -- a :class:`SageInputExpression` representing the subscript/key

        As a special case, ``key`` may be ``None``; this represents an
        empty subscript.  This is not legal Python syntax, but it is legal
        in the \sage preparser in examples like ``K.<x> = QQ[]``.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.name('QQ')['x']
            {subscr: {atomic:QQ}[{atomic:'x'}]}
            sage: sib.name('x')[1,2,3]
            {subscr: {atomic:x}[{tuple: ({atomic:1}, {atomic:2}, {atomic:3})}]}
            sage: sib.empty_subscript(sib.name('QQ'))
            {subscr: {atomic:QQ}[]}
        """
        super().__init__(sib)
        self._sie_coll = coll
        self._sie_key = key

    def __repr__(self):
        r"""
        Return a string representing this :class:`SIE_subscript` value.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.name('ZZ')['x,y']
            {subscr: {atomic:ZZ}[{atomic:'x,y'}]}
        """
        coll = repr(self._sie_coll)
        if self._sie_key is None:
            key = ''
        else:
            key = repr(self._sie_key)
        return "{subscr: %s[%s]}" % (coll, key)

    def _sie_referenced(self):
        r"""
        Return a list of the immediate subexpressions of this
        :class:`SIE_subscript`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib.name('GF')(5)['x,y']
            sage: sie._sie_referenced()
            [{call: {atomic:GF}({atomic:5})}, {atomic:'x,y'}]
        """
        refs = [self._sie_coll]
        if self._sie_key is not None:
            refs.append(self._sie_key)
        return refs

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this expression, and an
        indication that it is a subscript.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.name('QQ')['x']
            sage: sie._sie_prepare(sif)
            sage: sie._sie_format(sif)
            ("QQ['x']", 36)
        """
        coll = sif.format(self._sie_coll, _prec_attribute)
        if self._sie_key is None:
            key = ''
        else:
            key = sif.format(self._sie_key, 0)
        return '%s[%s]' % (coll, key), _prec_subscript


class SIE_getattr(SageInputExpression):
    r"""
    This class represents a getattr node in a :func:`sage_input`
    expression tree.

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sie = sib.name('CC').gen()
        sage: sie
        {call: {getattr: {atomic:CC}.gen}()}
    """
    def __init__(self, sib, obj, attr):
        r"""
        Initialize an instance of :class:`SIE_getattr`.

        INPUT:

        - ``sib`` -- a :class:`SageInputBuilder`

        - ``obj`` -- a :class:`SageInputExpression` representing an object

        - ``attr`` -- string; the attribute name

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.name('QQbar').zeta(5)
            {call: {getattr: {atomic:QQbar}.zeta}({atomic:5})}
        """
        super().__init__(sib)
        self._sie_obj = obj
        self._sie_attr = attr

    def __repr__(self):
        r"""
        Return a string representing this :class:`SIE_getattr` value.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.name('AA')(3).sqrt()
            {call: {getattr: {call: {atomic:AA}({atomic:3})}.sqrt}()}
        """
        obj = repr(self._sie_obj)
        return "{getattr: %s.%s}" % (obj, self._sie_attr)

    def _sie_referenced(self):
        r"""
        Return a list of the immediate subexpressions of this
        :class:`SIE_subscript`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib.name('CDF').gen
            sage: sie._sie_referenced()
            [{atomic:CDF}]
        """
        return [self._sie_obj]

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this expression, and an
        indication that it is an attribute reference.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.name('AA').common_polynomial
            sage: sie._sie_prepare(sif)
            sage: sie._sie_format(sif)
            ('AA.common_polynomial', 34)
        """
        obj = sif.format(self._sie_obj, _prec_exponent)
        return '%s.%s' % (obj, self._sie_attr), _prec_attribute


class SIE_tuple(SageInputExpression):
    r"""
    This class represents a tuple or list node in a :func:`sage_input`
    expression tree.

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sib((1, 'howdy'))
        {tuple: ({atomic:1}, {atomic:'howdy'})}
        sage: sib(["lists"])
        {list: ({atomic:'lists'})}
    """

    def __init__(self, sib, values, is_list):
        r"""
        Initialize an instance of :class:`SIE_tuple`.

        INPUT:

        - ``sib`` -- a :class:`SageInputBuilder`

        - ``values`` -- list of instances of :class:`SageInputExpression`
          representing the elements of this tuple

        - ``is_list`` -- is ``True`` if this class represents a list, False for a
          tuple

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib((3.5, -2))                                                        # needs sage.rings.real_mpfr
            {tuple: ({atomic:3.5}, {unop:- {atomic:2}})}
            sage: sib(["Hello", "world"])
            {list: ({atomic:'Hello'}, {atomic:'world'})}
        """
        super().__init__(sib)
        self._sie_values = values
        self._sie_is_list = is_list

    def __repr__(self):
        r"""
        Return a string representing this :class:`SIE_tuple` value.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib((2,3,5))
            {tuple: ({atomic:2}, {atomic:3}, {atomic:5})}
            sage: sib(["Hello", "world"])
            {list: ({atomic:'Hello'}, {atomic:'world'})}
        """
        kind = "list" if self._sie_is_list else "tuple"
        return "{%s: (%s)}" % \
            (kind, ', '.join(repr(v) for v in self._sie_values))

    def _sie_referenced(self):
        r"""
        Return a list of the immediate subexpressions of this
        :class:`SIE_tuple`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib((ZZ, GF(5)))
            sage: sie._sie_referenced()
            [{atomic:ZZ}, {call: {atomic:GF}({atomic:5})}]
        """
        return self._sie_values

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this tuple or list, and an
        indication that it is atomic (never needs to be parenthesized).

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: for v in ((), (1,), (1,2), [], [1], [1,2]):
            ....:      sie = sib(v)
            ....:      sie._sie_prepare(sif)
            ....:      sie._sie_format(sif)
            ('()', 42)
            ('(1,)', 42)
            ('(1, 2)', 42)
            ('[]', 42)
            ('[1]', 42)
            ('[1, 2]', 42)
        """
        values = [sif.format(val, 0) for val in self._sie_values]
        if self._sie_is_list:
            return '[%s]' % ', '.join(values), _prec_atomic
        else:
            if len(values) == 1:
                return '(%s,)' % values[0], _prec_atomic
            else:
                return '(%s)' % ', '.join(values), _prec_atomic


class SIE_dict(SageInputExpression):
    r"""
    This class represents a dict node in a :func:`sage_input`
    expression tree.

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sib.dict([('TeX', RR(pi)), ('Metafont', RR(e))])                          # needs sage.symbolic
        {dict: {{atomic:'TeX'}:{call: {atomic:RR}({atomic:3.1415926535897931})},
                {atomic:'Metafont'}:{call: {atomic:RR}({atomic:2.7182818284590451})}}}
        sage: sib.dict({-40:-40, 0:32, 100:212})
        {dict: {{unop:- {atomic:40}}:{unop:- {atomic:40}},
                {atomic:0}:{atomic:32}, {atomic:100}:{atomic:212}}}
    """

    def __init__(self, sib, entries):
        r"""
        Initialize an instance of :class:`SIE_dict`.

        INPUT:

        - ``sib`` -- a :class:`SageInputBuilder`

        - ``entries`` -- list of pairs of :class:`SageInputExpression`
          representing the entries of this dict

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.dict({'me':'good', 'you':'bad'})
            {dict: {{atomic:'me'}:{atomic:'good'}, {atomic:'you'}:{atomic:'bad'}}}
            sage: sib.dict([(10, 'PS2'), (12, 'PS2'), (13, 'PS3')])
            {dict: {{atomic:10}:{atomic:'PS2'}, {atomic:12}:{atomic:'PS2'}, {atomic:13}:{atomic:'PS3'}}}
        """
        super().__init__(sib)
        self._sie_entries = entries

    def __repr__(self):
        r"""
        Return a string representing this :class:`SIE_dict` value.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.dict({'keaton':'general', 'chan':'master'})
            {dict: {{atomic:'keaton'}:{atomic:'general'}, {atomic:'chan'}:{atomic:'master'}}}
        """
        return "{dict: {%s}}" % \
            ', '.join(repr(key) + ':' + repr(val)
                      for key, val in self._sie_entries)

    def _sie_referenced(self):
        r"""
        Return a list of the immediate subexpressions of this
        :class:`SIE_dict`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib.dict({1:'beguilement', 2:'legacy', 3:'passage'})
            sage: sie._sie_referenced()
            [{atomic:1}, {atomic:2}, {atomic:3}, {atomic:'beguilement'}, {atomic:'legacy'}, {atomic:'passage'}]
        """
        return [k for k, v in self._sie_entries] + [v for k, v in self._sie_entries]

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this dict, and an
        indication that it is atomic (never needs to be parenthesized).

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.dict({'carnivores':1, 'thinking':2, 'triumph':3})
            sage: sie._sie_prepare(sif)
            sage: sie._sie_format(sif)
            ("{'carnivores':1, 'thinking':2, 'triumph':3}", 42)
        """
        return "{%s}" % ', '.join(sif.format(k, 0) + ':' + sif.format(v, 0)
                                  for k, v in self._sie_entries), _prec_atomic


class SIE_binary(SageInputExpression):
    r"""
    This class represents an arithmetic expression with a binary operator
    and its two arguments, in a :func:`sage_input` expression tree.

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sib(3)+5
        {binop:+ {atomic:3} {atomic:5}}
    """

    def __init__(self, sib, op, lhs, rhs):
        r"""
        Initialize an instance of :class:`SIE_binary`.

        INPUT:

        - ``sib`` -- a :class:`SageInputBuilder`

        - ``op`` -- string representing a binary operator, such as '*' or '%'

        - ``lhs`` -- a :class:`SageInputExpression`

        - ``rhs`` -- a :class:`SageInputExpression`

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib(3)*5
            {binop:* {atomic:3} {atomic:5}}
        """
        super().__init__(sib)
        self._sie_op = op
        self._sie_operands = (lhs, rhs)

    def __repr__(self):
        r"""
        Return a string representing this :class:`SIE_binary` value.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib(7)/9
            {binop:/ {atomic:7} {atomic:9}}
        """
        return "{binop:%s %s %s}" % (self._sie_op, repr(self._sie_operands[0]), repr(self._sie_operands[1]))

    def _sie_referenced(self):
        r"""
        Return a tuple of the immediate subexpressions of this
        :class:`SIE_binary`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib.name('x') + 5
            sage: sie._sie_referenced()
            ({atomic:x}, {atomic:5})
        """
        return self._sie_operands

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this expression,
        and the precedence of the top-level operator in the expression.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: x = sib.name('x')
            sage: y = sib.name('y')
            sage: for v in (x+y, x*y, x**y):
            ....:     v._sie_prepare(sif)
            ....:     v._sie_format(sif)
            ('x + y', 24)
            ('x*y', 26)
            ('x^y', 32)

        Note that the printing for $x^y$ varies depending on whether the
        preparser is enabled.::

            sage: sibnp = SageInputBuilder(preparse=False)
            sage: sif = SageInputFormatter()
            sage: v = x**y
            sage: v._sie_prepare(sif)
            sage: v._sie_format(sif)
            ('x^y', 32)

        TESTS::

            sage: x = sib.name('x')
            sage: y = sib.name('y')
            sage: z = sib.name('z')
            sage: sib.result((x+y)+z)
            x + y + z
            sage: sib.result(x+(y+z))
            x + (y + z)
            sage: sib.result((x*y)*z)
            x*y*z
            sage: sib.result(x*(y*z))
            x*(y*z)
            sage: sib.result(x+(y*z))
            x + y*z
            sage: sib.result((x+y)*z)
            (x + y)*z
            sage: sib.result((x^y)^z)
            (x^y)^z
            sage: sib.result(x^(y^z))
            x^y^z
        """
        op = self._sie_op
        fop = op
        if op == '**':
            lhs = sif.format(self._sie_operands[0], _prec_exponent + 1)
            rhs = sif.format(self._sie_operands[1], _prec_exponent)
            if self._sie_builder.preparse():
                return '%s^%s' % (lhs, rhs), _prec_exponent
            else:
                return '%s**%s' % (lhs, rhs), _prec_exponent

        if op == '*':
            prec = _prec_muldiv
        elif op == '/':
            prec = _prec_muldiv
        elif op == '+':
            fop = ' + '
            prec = _prec_addsub
        elif op == '-':
            fop = ' - '
            prec = _prec_addsub
        else:
            raise ValueError('Unhandled op {} in SIE_binary'.format(op))

        lhs = sif.format(self._sie_operands[0], prec)
        rhs = sif.format(self._sie_operands[1], prec + 1)
        return '%s%s%s' % (lhs, fop, rhs), prec


class SIE_unary(SageInputExpression):
    r"""
    This class represents an arithmetic expression with a unary operator
    and its argument, in a :func:`sage_input` expression tree.

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: -sib(256)
        {unop:- {atomic:256}}
    """

    def __init__(self, sib, op, operand):
        r"""
        Initialize an instance of :class:`SIE_unary`.

        INPUT:

        - ``sib`` -- a :class:`SageInputBuilder`

        - ``op`` -- string representing a unary operator, such as '-'

        - ``operand`` -- a :class:`SageInputExpression`

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: -sib(3)
            {unop:- {atomic:3}}
        """
        super().__init__(sib)
        self._sie_op = op
        self._sie_operand = operand

    def __repr__(self):
        r"""
        Return a string representing this :class:`SIE_unary` value.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: -sib(15)
            {unop:- {atomic:15}}
        """
        return "{unop:%s %s}" % (self._sie_op, repr(self._sie_operand))

    def _sie_referenced(self):
        r"""
        Return a list of the immediate subexpressions of this
        :class:`SIE_unary`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = -sib.name('x')
            sage: sie._sie_referenced()
            [{atomic:x}]
        """
        return [self._sie_operand]

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this expression,
        and the precedence of the top-level operator in the expression.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: x = sib.name('x')
            sage: v = -x
            sage: v._sie_prepare(sif)
            sage: v._sie_format(sif)
            ('-x', 28)
            sage: v = ~x
            sage: v._sie_prepare(sif)
            sage: v._sie_format(sif)
            ('~x', 30)

        TESTS::

            sage: x = sib.name('x')
            sage: y = sib.name('y')
            sage: sib.result((-x)+y)
            -x + y
            sage: sib.result(x+(-y))
            x + -y
            sage: sib.result(-(x+y))
            -(x + y)
            sage: sib.result(-(-x))
            --x
            sage: sib.result(x-(-y))
            x - -y
            sage: sib.result((+x)+y)
            +x + y
            sage: sib.result(x+(+y))
            x + +y
            sage: sib.result(+(x+y))
            +(x + y)
            sage: sib.result(+(+x))
            ++x
            sage: sib.result(x+(+y))
            x + +y

        We assume that -(x*y) is always equal to (-x)*y.  Using this
        assumption, we print -(x*y) as -x*y, which parses as (-x)*y.::

            sage: sib.result(-(x*y))
            -x*y
            sage: sib.result((-x)*y)
            -x*y
            sage: sib.result(x*(-y))
            x*-y

        We do not do that for unary +, assuming that the user really
        means to express something by using unary +::

            sage: sib.result(+(x*y))
            +(x*y)
        """
        op = self._sie_op
        fop = op
        rprec = None
        if op == '-':
            # We print -(a*b) as -a*b, even though that will parse as
            # (-a)*b.
            prec = _prec_muldiv
            rprec = _prec_negate
        elif op == '+':
            prec = _prec_negate
        elif op == '~':
            prec = _prec_bitnot
        else:
            raise ValueError('Unhandled op {} in SIE_unary'.format(op))

        if rprec is None:
            rprec = prec

        return '%s%s' % (fop, sif.format(self._sie_operand, prec)), rprec

    def _sie_is_negation(self):
        r"""
        Test whether a :class:`SageInputExpression` is a negation.

        Despite the obscure name, this is intended to be a public method.

        This is used in the \method{_sage_input_} method for
        :class:`ComplexNumber`, so that ``sage_input(CC(-3))`` will
        produce ``-CC(3)`` instead of ``CC(-3)``.  (This is preferred
        so that you get ``x - CC(3)`` instead of ``x + CC(-3)``.)

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: x = sib.name('x')
            sage: v = -x

            sage: def mk_CC(b):
            ....:     if b._sie_is_negation():
            ....:         return -sib.name('CC')(b._sie_operand)
            ....:     else:
            ....:         return sib.name('CC')(b)

            sage: mk_CC(x)
            {call: {atomic:CC}({atomic:x})}
            sage: mk_CC(v)
            {unop:- {call: {atomic:CC}({atomic:x})}}
        """
        return self._sie_op == '-'


class SIE_gens_constructor(SageInputExpression):
    r"""
    This class represents an expression that can create a \sage parent
    with named generators, optionally using the \sage preparser
    generators syntax (like ``K.<x> = QQ[]``).

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: qq = sib.name('QQ')
        sage: sib.parent_with_gens("some parent", qq['x'],
        ....:                      ('x',), 'QQx',
        ....:                      gens_syntax=sib.empty_subscript(qq))
        {constr_parent: {subscr: {atomic:QQ}[{atomic:'x'}]} with gens: ('x',)}
    """

    def __init__(self, sib, constr, gen_names, gens_syntax=None):
        r"""
        Initialize an instance of :class:`SIE_gens_constructor`.

        INPUT:

        - ``sib`` -- a :class:`SageInputBuilder`

        - ``constr`` -- a :class:`SageInputExpression` for constructing this
          parent ``normally``

        - ``gen_names`` -- tuple of generator names

        - ``gens_syntax`` -- an optional :class:`SageInputExpression` for
          constructing this parent using the \sage preparser generators syntax

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: qq = sib.name('QQ')
            sage: sib.parent_with_gens("some parent", qq['x'],
            ....:                      ('x',), 'QQx',
            ....:                      gens_syntax=sib.empty_subscript(qq))
            {constr_parent: {subscr: {atomic:QQ}[{atomic:'x'}]} with gens: ('x',)}
        """
        super().__init__(sib)
        self._sie_constr = constr
        self._sie_gen_names = gen_names
        self._sie_gens = None  # will be overwritten from .parent_with_gens()
        self._sie_gens_constr = gens_syntax
        self._sie_assign_gens = False
        self._sie_generated = False

    def __repr__(self):
        r"""
        Return a string representing this :class:`SIE_gens_constructor` value.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: qq = sib.name('QQ')
            sage: sib.parent_with_gens("some parent", qq['x'],
            ....:                      ('x',), 'QQx',
            ....:                      gens_syntax=sib.empty_subscript(qq))
            {constr_parent: {subscr: {atomic:QQ}[{atomic:'x'}]} with gens: ('x',)}
        """
        return "{constr_parent: %s with gens: %s}" % (repr(self._sie_constr), self._sie_gen_names)

    def _sie_referenced(self):
        r"""
        Return a list of the immediate subexpressions of this
        :class:`SIE_gens_constructor`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: qq = sib.name('QQ')
            sage: gc = sib.parent_with_gens("some parent", qq['x'],
            ....:                           ('x',), 'QQx',
            ....:                           gens_syntax=sib.empty_subscript(qq))
            sage: gc._sie_referenced()
            [{subscr: {atomic:QQ}[{atomic:'x'}]}]
        """
        # This is used to determine if some expressions should be replaced
        # by variables (if the expression has more than one parent in
        # the expression DAG).  We assume that all expressions in
        # self._sie_gens_constr also occur in self._sie_constr.
        return [self._sie_constr]

    def _sie_gens_referenced(self, sif):
        r"""
        Mark that at least one of the generators in this
        :class:`SIE_gens_constructor` is used.  (This means we will actually
        construct all of the generators.)

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: qq = sib.name('QQ')
            sage: gc = sib.parent_with_gens("some parent", qq['x'],
            ....:                           ('x',), 'QQx',
            ....:                           gens_syntax=sib.empty_subscript(qq))
            sage: gc._sie_assign_gens
            False
            sage: gc._sie_gens_referenced(sif)
            sage: gc._sie_assign_gens
            True
        """
        self._sie_assign_gens = True
        self._sie_require_varname(sif)
        for gen in self._sie_gens:
            gen._sie_require_varname(sif)

    def _sie_add_command(self, sif):
        r"""
        Build commands to construct this parent and (if necessary)
        its associated generators.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: qq = sib.name('QQ')
            sage: gc = sib.parent_with_gens("some parent", qq['x'],
            ....:                           ('x',), 'QQx',
            ....:                           gens_syntax=sib.empty_subscript(qq))
            sage: gc._sie_gens_referenced(sif)
            sage: gc._sie_prepare(sif)
            sage: gc._sie_add_command(sif)
            sage: sif._commands
            'QQx.<x> = QQ[]\n'

        TESTS:

        There are several tricky cases here.

        We prefer the \sage preparser generators syntax::

            sage: sage_input(polygen(ZZ))
            R.<x> = ZZ[]
            x

        But of course we can't use that without the preparser::

            sage: sage_input(polygen(ZZ), preparse=False)
            R = ZZ['x']
            x = R.gen()
            x

        We also can't use the preparser syntax if there is a conflict
        between generator names.  For example, this works::

            sage: sage_input((polygen(ZZ), polygen(GF(17), 'y')))
            R1.<x> = ZZ[]
            R2.<y> = GF(17)[]
            (x, y)

        but this can't use the preparser syntax.::

            sage: sage_input((polygen(ZZ), polygen(GF(17))))
            R1 = ZZ['x']
            x1 = R1.gen()
            R2 = GF(17)['x']
            x2 = R2.gen()
            (x1, x2)

        If we never use the generators, then we don't bother with the
        preparser syntax.::

            sage: sage_input((ZZ['x'], ZZ['x'], GF(17)['y']))
            R = ZZ['x']
            (R, R, GF(17)['y'])
        """
        if not self._sie_generated:
            if self._sie_builder.preparse() and \
                    self._sie_gens_constr is not None and \
                    all(g._sie_got_preferred(sif) for g in self._sie_gens):
                s, _ = self._sie_gens_constr._sie_format(sif)
                sif._commands += '%s.<%s> = %s\n' % (self._sie_get_varname(sif), ','.join(self._sie_gen_names), s)
            else:
                s, _ = self._sie_constr._sie_format(sif)
                sif._commands += '%s = %s\n' % (self._sie_get_varname(sif), s)
                if self._sie_assign_gens:
                    if len(self._sie_gens) == 1:
                        sif._commands += '%s = %s.gen()\n' % (self._sie_gens[0]._sie_get_varname(sif), self._sie_get_varname(sif))
                    else:
                        sif._commands += '%s = %s.gens()\n' % (','.join(g._sie_get_varname(sif) for g in self._sie_gens), self._sie_get_varname(sif))
            self._sie_generated = True

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this parent-construction
        expression, and its precedence.

        As a side effect, if the generators of this parent are used,
        this adds commands to assign the generators to names.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: qq = sib.name('QQ')
            sage: gc = sib.parent_with_gens("some parent", qq['x'],
            ....:                           ('x',), 'QQx',
            ....:                           gens_syntax=sib.empty_subscript(qq))
            sage: gc._sie_gens_referenced(sif)
            sage: gc._sie_prepare(sif)
            sage: gc._sie_format(sif)
            ('QQx', 42)
            sage: sif._commands
            'QQx.<x> = QQ[]\n'
        """
        if self._sie_assign_gens:
            self._sie_add_command(sif)
            return self._sie_get_varname(sif), _prec_atomic

        return self._sie_constr._sie_format(sif)


class SIE_gen(SageInputExpression):
    r"""
    This class represents a named generator of a parent with named
    generators.

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sib.gen(ZZ['x'])
        {gen:x {constr_parent: {subscr: {atomic:ZZ}[{atomic:'x'}]} with gens: ('x',)}}
    """

    def __init__(self, sib, parent, name):
        r"""
        Initialize an instance of :class:`SIE_gen`.

        INPUT:

        - ``sib`` -- a :class:`SageInputBuilder`

        - ``parent`` -- a :class:`SIE_gens_constructor`

        - ``name`` -- string with the name of this generator

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.gen(ZZ['x']) # indirect doctest
            {gen:x {constr_parent: {subscr: {atomic:ZZ}[{atomic:'x'}]} with gens: ('x',)}}
        """
        super().__init__(sib)
        self._sie_parent = parent
        self._sie_preferred_varname = name

    def __repr__(self):
        r"""
        Return a string representing this :class:`SIE_gen` value.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.gen(ZZ['x']) # indirect doctest
            {gen:x {constr_parent: {subscr: {atomic:ZZ}[{atomic:'x'}]} with gens: ('x',)}}
        """
        return "{gen:%s %s}" % (self._sie_preferred_varname, repr(self._sie_parent))

    def _sie_is_simple(self):
        r"""
        Report that :class:`SIE_gen` values are single tokens.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.gen(ZZ['x'])._sie_is_simple()
            True
        """
        return True

    def _sie_prepare(self, sif):
        r"""
        We override the \method{_sie_prepare} method from
        :class:`SageInputExpression` to additionally mark the parent of this
        generator that the generator names must be assigned.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.gen(GF(13)['z'])
            sage: sie._sie_parent._sie_assign_gens
            False
            sage: sie._sie_prepare(sif)
            sage: sie._sie_parent._sie_assign_gens
            True
        """
        super()._sie_prepare(sif)
        self._sie_parent._sie_gens_referenced(sif)

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this named generator,
        and an indication that it is atomic.

        As a side effect, this generates commands to assign the generators
        of the parent to variables.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.gen(GF(41)['x'])
            sage: sie._sie_prepare(sif)
            sage: sie._sie_format(sif)
            ('x', 42)
            sage: sif._commands
            'R.<x> = GF(41)[]\n'
        """
        self._sie_parent._sie_add_command(sif)
        return self._sie_get_varname(sif), _prec_atomic

    def _sie_got_preferred(self, sif):
        r"""
        Check whether the :class:`SageInputFormatter` assigned us a
        variable name which is the same as the name of the generator
        name.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

        First we verify that if we use two generators with different
        names, then they get their preferred names.::

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: v = sib.gen(GF(2)['x']); w = sib.gen(GF(3)['y'])
            sage: v._sie_prepare(sif); w._sie_prepare(sif)
            sage: v._sie_got_preferred(sif)
            True
            sage: w._sie_got_preferred(sif)
            True

        Now, we repeat the experiment, except that the generators now
        have the same names.  In this case, the :class:`SageInputFormatter`
        will not use the generator name as the variable name, because
        of this conflict.::

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: v = sib.gen(GF(2)['x']); w = sib.gen(GF(3)['x'])
            sage: v._sie_prepare(sif); w._sie_prepare(sif)
            sage: v._sie_got_preferred(sif)
            False
            sage: w._sie_got_preferred(sif)
            False
        """
        return self._sie_get_varname(sif) == self._sie_preferred_varname


class SIE_import_name(SageInputExpression):
    r"""
    This class represents a name which has been imported from a module.

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sib.import_name('sage.rings.integer', 'make_integer')
        {import:sage.rings.integer/make_integer}
        sage: sib.import_name('sage.foo', 'happy', 'sad')
        {import:sage.foo/happy as sad}
    """

    def __init__(self, sib, module, name, alt_name=None):
        r"""
        Initialize an instance of :class:`SIE_import_name`.

        INPUT:

        - ``sib`` -- a :class:`SageInputBuilder`

        - ``module`` -- a module name

        - ``name`` -- an object name

        - ``alt_name`` -- an alternate object name, or None (the default)
                      to use name

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.import_name('sage.rings.integer', 'make_integer') # indirect doctest
            {import:sage.rings.integer/make_integer}
            sage: sib.import_name('sage.foo', 'happy', 'sad')
            {import:sage.foo/happy as sad}
        """
        super().__init__(sib)
        self._sie_formatted = False
        self._sie_module_name = module
        self._sie_object_name = name
        if alt_name is None:
            self._sie_preferred_varname = name
        else:
            self._sie_preferred_varname = alt_name

    def __repr__(self):
        r"""
        Return a string representing this :class:`SIE_import_name` value.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.import_name('sage.rings.integer', 'make_integer') # indirect doctest
            {import:sage.rings.integer/make_integer}
            sage: sib.import_name('sage.foo', 'happy', 'sad')
            {import:sage.foo/happy as sad}
        """
        return "{import:%s/%s%s}" % (self._sie_module_name, self._sie_object_name,
                                     "" if self._sie_object_name == self._sie_preferred_varname else " as %s" % self._sie_preferred_varname)

    def _sie_is_simple(self):
        r"""
        Report that :class:`SIE_import_name` values are single tokens.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.import_name('sage.rings.integer', 'make_integer')._sie_is_simple()
            True
        """
        return True

    def _sie_prepare(self, sif):
        r"""
        We override the \method{_sie_prepare} method from
        :class:`SageInputExpression` to request a variable name.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.import_name('sage.rings.integer', 'make_integer')
            sage: sie._sie_requested_varname
            False
            sage: sie._sie_prepare(sif)
            sage: sie._sie_requested_varname
            True
        """
        super()._sie_prepare(sif)
        self._sie_require_varname(sif)

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this import,
        and an indication that it is atomic.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: v1 = sib.import_name('sage.rings.integer', 'make_integer')
            sage: v2 = sib.import_name('sage.foo', 'happy', 'sad')
            sage: sie = v1(v2)
            sage: sie._sie_prepare(sif)
            sage: sie._sie_format(sif)
            ('make_integer(sad)', 40)
            sage: print(sif._commands)
            from sage.rings.integer import make_integer
            from sage.foo import happy as sad
        """
        name = self._sie_get_varname(sif)
        if self._sie_formatted:
            # Only run the import command once
            return name, _prec_atomic
        self._sie_formatted = True
        rename = ''
        if name != self._sie_object_name:
            rename = ' as ' + name
        sif._commands += 'from %s import %s%s\n' % (self._sie_module_name,
                                                    self._sie_object_name,
                                                    rename)
        return name, _prec_atomic


class SIE_assign(SageInputExpression):
    r"""
    This class represents an assignment command.

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sib.assign(sib.name('foo').x, sib.name('pi'))
        {assign: {getattr: {atomic:foo}.x} {atomic:pi}}
    """

    def __init__(self, sib, lhs, rhs):
        r"""
        Initialize an instance of :class:`SIE_assign`.

        INPUT:

        - ``sib`` -- a :class:`SageInputBuilder`

        - ``lhs`` -- the left-hand side of the assignment

        - ``rhs`` -- the right-hand side of the assignment

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.assign(sib.name('foo').x, sib.name('pi'))
            {assign: {getattr: {atomic:foo}.x} {atomic:pi}}
        """
        super().__init__(sib)
        self._sie_lhs = lhs
        self._sie_rhs = rhs

    def __repr__(self):
        r"""
        Return a string representing this :class:`SIE_assign` command.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.assign(sib.name('foo').x, sib.name('pi'))
            {assign: {getattr: {atomic:foo}.x} {atomic:pi}}
        """
        return "{assign: %s %s}" % (repr(self._sie_lhs), repr(self._sie_rhs))

    def _sie_referenced(self):
        r"""
        Return a list of the immediate subexpressions of this
        :class:`SIE_assign`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib.assign(sib.name('foo').x, sib.name('pi'))
            sage: sie._sie_referenced()
            [{getattr: {atomic:foo}.x}, {atomic:pi}]
        """
        return [self._sie_lhs, self._sie_rhs]

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this :class:`SIE_assign`
        as an expression.  Since an assignment is a statement, not
        an expression, always raises an error.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.assign(sib.name('foo').x, sib.name('pi'))
            sage: sie._sie_prepare(sif)
            sage: sie._sie_format(sif)
            Traceback (most recent call last):
            ...
            ValueError: Cannot format SIE_assign as expression
        """
        raise ValueError("Cannot format SIE_assign as expression")

    def _sie_format_statement(self, sif):
        r"""
        Return the formatted string of this :class:`SIE_assign`
        as a statement.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter
            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.assign(sib.name('foo').x, sib.name('pi'))
            sage: sie._sie_prepare(sif)
            sage: sie._sie_format_statement(sif)
            'foo.x = pi'
        """
        return '%s = %s' % (sif.format(self._sie_lhs, 0), sif.format(self._sie_rhs, 0))


class SageInputFormatter:
    r"""
    An instance of this class is used to keep track of variable names
    and a sequence of generated commands during the :func:`sage_input`
    formatting process.
    """

    def __init__(self):
        r"""
        Initialize an instance of :class:`SageInputFormatter`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputFormatter
            sage: sif = SageInputFormatter()
        """
        self._commands = ''
        self._names = set()
        self._dup_names = {}

    def format(self, e, prec):
        r"""
        Format a Sage input expression into a string.

        INPUT:

        - ``e`` -- a :class:`SageInputExpression`

        - ``prec`` -- integer representing a precedence level

        First, we check to see if ``e`` should be replaced by a variable.
        If so, we generate the command to assign the variable, and return
        the name of the variable.

        Otherwise, we format the expression by calling its \method{_sie_format}
        method, and add parentheses if necessary.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib(GF(5))

        Here we ``cheat`` by calling \method{_sie_prepare} twice, to make it
        use a variable.::

            sage: sie._sie_prepare(sif)
            sage: sie._sie_prepare(sif)
            sage: sif._commands
            ''
            sage: sif.format(sie, 0)
            'GF_5'
            sage: sif._commands
            'GF_5 = GF(5)\n'

        We demonstrate the use of commands, by showing how to construct
        code that will produce a random matrix::

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.name('matrix')(sib.name('ZZ'), 10, 10)
            sage: sib.command(sie, sie.randomize())
            sage: sie._sie_prepare(sif)
            sage: sif._commands
            ''
            sage: sif.format(sie, 0)
            'si'
            sage: sif._commands
            'si = matrix(ZZ, 10, 10)\nsi.randomize()\n'
        """
        if e._sie_use_var:
            if not e._sie_generated:
                s, _ = e._sie_format(self)
                # In complicated situations, this can get called
                # recursively...
                if not e._sie_generated:
                    self._commands += '%s = %s\n' % (e._sie_get_varname(self), s)
                    e._sie_generated = True

            formatted = e._sie_get_varname(self)
        else:
            s, iprec = e._sie_format(self)
            if iprec < prec:
                s = '(' + s + ')'
            formatted = s

        commands = e._sie_commands
        e._sie_commands = []

        for cmd in commands:
            s_cmd = cmd._sie_format_statement(self)
            self._commands += s_cmd + '\n'

        return formatted

    def register_name(self, name):
        r"""
        Register that some value would like to use a given name.
        If only one request for a name is received, then we will use the
        requested name; otherwise, we will add numbers to the end of the
        name to make it unique.

        If the input name is ``None``, then it is treated as a name of
        ``'si'``.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputFormatter

            sage: sif = SageInputFormatter()
            sage: sif._names, sif._dup_names
            (set(), {})
            sage: sif.register_name('x')
            sage: sif.register_name('y')
            sage: sif._names, sif._dup_names
            ({'x', 'y'}, {})
            sage: sif.register_name('x')
            sage: sif._names, sif._dup_names
            ({'x', 'y'}, {'x': 0})
        """
        if name is None:
            name = 'si'

        if name in self._names:
            self._dup_names[name] = 0
        else:
            self._names.add(name)

    def get_name(self, name):
        r"""
        Return a name corresponding to a given requested name.
        If only one request for a name is received, then we will use the
        requested name; otherwise, we will add numbers to the end of the
        name to make it unique.

        If the input name is ``None``, then it is treated as a name of
        ``'si'``.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputFormatter

            sage: sif = SageInputFormatter()
            sage: names = ('x', 'x', 'y', 'z')
            sage: for n in names: sif.register_name(n)
            sage: for n in names: sif.get_name(n)
            'x1'
            'x2'
            'y'
            'z'
        """
        if name is None:
            name = 'si'

        if name in self._dup_names:
            next = self._dup_names[name] + 1
            self._dup_names[name] = next
            return name + str(next)
        else:
            return name


def verify_same(a, b):
    r"""
    Verify that two Sage values are the same.  This is an extended equality
    test; it checks that the values are equal and that their parents are equal.
    (For values which are not Elements, the types are checked instead.)

    If the values are the same, we return ``None``; otherwise,
    we raise an exception.

    EXAMPLES::

        sage: from sage.misc.sage_input import verify_same
        sage: verify_same(1, 1)
        sage: verify_same(1, 2)
        Traceback (most recent call last):
        ...
        AssertionError: Expected 1 == 2
        sage: verify_same(1, 1r)
        Traceback (most recent call last):
        ...
        AttributeError: 'int' object has no attribute 'parent'...
        sage: verify_same(1r, 1)
        Traceback (most recent call last):
        ...
            assert(type(a) == type(b))
        AssertionError
        sage: verify_same(5, GF(7)(5))
        Traceback (most recent call last):
        ...
            assert(a.parent() == b.parent())
        AssertionError
    """
    from sage.structure.element import Element
    if isinstance(a, Element):
        assert a.parent() == b.parent()
    else:
        assert type(a) is type(b)
    if isinstance(a, (RealIntervalFieldElement, ComplexIntervalFieldElement)):
        assert a.endpoints() == b.endpoints(), "Expected %s == %s" % (a, b)
        return

    if not (a == b):
        # Verification failed => raise an AssertionError.
        #
        # There is an important exception: the IEEE-754 standard
        # recommends that NaN != NaN. So this comparison will fail for
        # any object involving NaN.
        #
        # If this case occurs, then a and b do not compare equal to
        # itself. In that case, we compare the string representations of
        # a and b.
        if not (a == a) and not (b == b):
            if repr(a) == repr(b):
                return  # Good!
        raise AssertionError("Expected %r == %r" % (a, b))


def verify_si_answer(x, answer, preparse):
    r"""
    Verify that evaluating ``answer`` gives a value equal to ``x``
    (with the same parent/type).  If ``preparse`` is ``True`` or
    ``False``, then we evaluate ``answer`` with the preparser
    enabled or disabled, respectively; if ``preparse`` is ``None``,
    then we evaluate ``answer`` both with the preparser enabled and
    disabled and check both results.

    On success, we return ``None``; on failure, we raise an exception.

    INPUT:

    - ``x`` -- an arbitrary Sage value

    - ``answer`` -- string, or a :class:`SageInputAnswer`

    - ``preparse`` -- ``True``, ``False``, or ``None``

    EXAMPLES::

        sage: from sage.misc.sage_input import verify_si_answer
        sage: verify_si_answer(1, '1', True)
        sage: verify_si_answer(1, '1', False)
        Traceback (most recent call last):
        ...
        AttributeError: 'int' object has no attribute 'parent'...
        sage: verify_si_answer(1, 'ZZ(1)', None)
    """
    from sage.misc.sage_eval import sage_eval
    if preparse is None:
        verify_same(x, sage_eval(answer, preparse=True))
        verify_same(x, sage_eval(answer, preparse=False))
    else:
        verify_same(x, sage_eval(answer, preparse=preparse))


class SageInputAnswer(tuple):
    r"""
    This class inherits from tuple, so it acts like a tuple when passed
    to :func:`sage_eval`; but it prints as a sequence of commands.

    EXAMPLES::

        sage: from sage.misc.sage_input import SageInputAnswer
        sage: v = SageInputAnswer('x = 22\n', 'x/7'); v
        x = 22
        x/7
        sage: isinstance(v, tuple)
        True
        sage: v[0]
        'x = 22\n'
        sage: v[1]
        'x/7'
        sage: len(v)
        2
        sage: v = SageInputAnswer('', 'sin(3.14)', {'sin': math.sin}); v
        LOCALS:
          sin: <built-in function sin>
        sin(3.14)
        sage: v[0]
        ''
        sage: v[1]
        'sin(3.14)'
        sage: v[2]
        {'sin': <built-in function sin>}
    """

    def __new__(cls, cmds, expr, locals=None):
        r"""
        Construct an instance of :class:`SageInputAnswer`.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputAnswer
            sage: v = SageInputAnswer('', 'sin(3.14)', {'sin': math.sin}); v
            LOCALS:
              sin: <built-in function sin>
            sin(3.14)
            sage: v[0]
            ''
            sage: v[1]
            'sin(3.14)'
            sage: v[2]
            {'sin': <built-in function sin>}
        """
        if locals:
            return tuple.__new__(cls, (cmds, expr, locals))
        else:
            return tuple.__new__(cls, (cmds, expr))

    def __repr__(self):
        r"""
        Return a string representation for a :class:`SageInputAnswer`,
        such that if you evaluate this :class:`SageInputAnswer` at the
        \sage command line, you get a result in a nice form ready to
        copy-and-paste.

        EXAMPLES::

            sage: from sage.misc.sage_input import SageInputAnswer
            sage: v = SageInputAnswer('', 'sin(3.14)', {'sin': math.sin}); v
            LOCALS:
              sin: <built-in function sin>
            sin(3.14)
        """
        if len(self) == 2:
            return self[0] + self[1]

        locals = self[2]
        locals_text = ''.join('  %s: %r\n' % (k, v)
                              for k, v in locals.items())
        return 'LOCALS:\n' + locals_text + self[0] + self[1]
