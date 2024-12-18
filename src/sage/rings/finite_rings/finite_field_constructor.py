# sage.doctest: needs sage.rings.finite_rings
r"""
Finite fields

Sage supports arithmetic in finite prime and extension fields.
Several implementation for prime fields are implemented natively in
Sage for several sizes of primes `p`. These implementations
are


- ``sage.rings.finite_rings.integer_mod.IntegerMod_int``,

- ``sage.rings.finite_rings.integer_mod.IntegerMod_int64``, and

- ``sage.rings.finite_rings.integer_mod.IntegerMod_gmp``.


Small extension fields of cardinality `< 2^{16}` are
implemented using tables of Zech logs via the Givaro C++ library
(``sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro``).
While this representation is very fast it is limited to finite
fields of small cardinality. Larger finite extension fields of
order `q >= 2^{16}` are internally represented as
polynomials over smaller finite prime fields. If the
characteristic of such a field is 2 then NTL is used internally to
represent the field
(``sage.rings.finite_rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e``).
In all other case the PARI C library is used
(``sage.rings.finite_rings.finite_field_pari_ffelt.FiniteField_pari_ffelt``).

However, this distinction is internal only and the user usually
does not have to worry about it because consistency across all
implementations is aimed for. In all extension field
implementations the user may either specify a minimal polynomial or
leave the choice to Sage.

For small finite fields the default choice are Conway polynomials.

The Conway polynomial `C_n` is the lexicographically first
monic irreducible, primitive polynomial of degree `n` over
`GF(p)` with the property that for a root `\alpha`
of `C_n` we have that
`\beta=
\alpha^{(p^n - 1)/(p^m - 1)}` is a root of
`C_m` for all `m` dividing `n`. Sage
contains a database of Conway polynomials which also can be queried
independently of finite field construction.

A pseudo-Conway polynomial satisfies all of the conditions required
of a Conway polynomial except the condition that it is lexicographically
first.  They are therefore not unique.  If no variable name is
specified for an extension field, Sage will fit the finite field
into a compatible lattice of field extensions defined by pseudo-Conway
polynomials. This lattice is stored in an
:class:`~sage.rings.algebraic_closure_finite_field.AlgebraicClosureFiniteField`
object; different algebraic closure objects can be created by using
a different ``prefix`` keyword to the finite field constructor.

Note that the computation of pseudo-Conway polynomials is expensive
when the degree is large and highly composite.  If a variable
name is specified then a random polynomial is used instead, which
will be much faster to find.

While Sage supports basic arithmetic in finite fields some more
advanced features for computing with finite fields are still not
implemented. For instance, Sage does not calculate embeddings of
finite fields yet.

EXAMPLES::

    sage: k = GF(5); type(k)
    <class 'sage.rings.finite_rings.finite_field_prime_modn.FiniteField_prime_modn_with_category'>

::

    sage: k = GF(5^2,'c'); type(k)                                                      # needs sage.libs.linbox
    <class 'sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro_with_category'>

One can also give the cardinality `q=p^n` as the tuple `(p,n)`::

    sage: k = GF((5, 2),'c'); k
    Finite Field in c of size 5^2

::

    sage: k = GF(2^16,'c'); type(k)                                                     # needs sage.libs.ntl
    <class 'sage.rings.finite_rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e_with_category'>

::

    sage: k = GF((3, 16),'c'); type(k)
    <class 'sage.rings.finite_rings.finite_field_pari_ffelt.FiniteField_pari_ffelt_with_category'>

Finite Fields support iteration, starting with 0.

::

    sage: k = GF(9, 'a')
    sage: for i,x in enumerate(k):  print("{} {}".format(i, x))
    0 0
    1 a
    2 a + 1
    3 2*a + 1
    4 2
    5 2*a
    6 2*a + 2
    7 a + 2
    8 1
    sage: for a in GF(5):
    ....:     print(a)
    0
    1
    2
    3
    4

We output the base rings of several finite fields.

::

    sage: k = GF(3); type(k)
    <class 'sage.rings.finite_rings.finite_field_prime_modn.FiniteField_prime_modn_with_category'>
    sage: k.base_ring()
    Finite Field of size 3

::

    sage: # needs sage.libs.linbox
    sage: k = GF(9,'alpha'); type(k)
    <class 'sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro_with_category'>
    sage: k.base_ring()
    Finite Field of size 3

::

    sage: k = GF((3, 40),'b'); type(k)
    <class 'sage.rings.finite_rings.finite_field_pari_ffelt.FiniteField_pari_ffelt_with_category'>
    sage: k.base_ring()
    Finite Field of size 3

Further examples::

    sage: GF(2).is_field()
    True
    sage: GF(next_prime(10^20)).is_field()
    True
    sage: GF(19^20,'a').is_field()
    True
    sage: GF(8,'a').is_field()
    True

AUTHORS:

- William Stein: initial version

- Robert Bradshaw: prime field implementation

- Martin Albrecht: Givaro and ntl.GF2E implementations
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections import defaultdict
from sage.structure.category_object import normalize_names, certify_names
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.integer import Integer

# the import below is just a redirection
from sage.rings.finite_rings.finite_field_base import is_FiniteField
assert is_FiniteField  # just to silent pyflakes

try:
    # We don't late import this because this means trouble with the Givaro library
    # On a Macbook Pro OSX 10.5.8, this manifests as a Bus Error on exiting Sage.
    # TODO: figure out why
    from .finite_field_givaro import FiniteField_givaro
except ImportError:
    FiniteField_givaro = None

try:
    from .finite_field_ntl_gf2e import FiniteField_ntl_gf2e
except ImportError:
    FiniteField_ntl_gf2e = None


from sage.structure.factory import UniqueFactory


class FiniteFieldFactory(UniqueFactory):
    """
    Return the globally unique finite field of given order with
    generator labeled by the given name and possibly with given
    modulus.

    INPUT:

    - ``order`` -- a prime power `p^n`, or a pair `(p,n)`

    - ``name`` -- string, optional.  Note that there can be a
      substantial speed penalty (in creating extension fields) when
      omitting the variable name, since doing so triggers the
      computation of pseudo-Conway polynomials in order to define a
      coherent lattice of extensions of the prime field.  The speed
      penalty grows with the size of extension degree and with
      the number of factors of the extension degree.
      Note that the second argument can also be the degree,
      where ``GF(p,n)`` is equivalent to ``GF((p,n))``.

    - ``modulus`` -- (optional) either a defining polynomial for the
      field, or a string specifying an algorithm to use to generate
      such a polynomial.  If ``modulus`` is a string, it is passed to
      :meth:`~sage.rings.polynomial.irreducible_element()` as the
      parameter ``algorithm``; see there for the permissible values of
      this parameter. In particular, you can specify
      ``modulus="primitive"`` to get a primitive polynomial.  You
      may not specify a modulus if you do not specify a variable name.

    - ``impl`` -- (optional) a string specifying the implementation of
      the finite field. Possible values are:

      - ``'modn'`` -- ring of integers modulo `p` (only for prime fields)

      - ``'givaro'`` -- Givaro, which uses Zech logs (only for fields
        of at most 65521 elements)

      - ``'ntl'`` -- NTL using GF2X (only in characteristic 2)

      - ``'pari'`` or ``'pari_ffelt'`` -- PARI's ``FFELT`` type (only
        for extension fields)

    - ``elem_cache`` -- (default: order < 500) cache all elements to
      avoid creation time; ignored unless ``impl='givaro'``

    - ``repr`` -- (default: ``'poly'``) ignored unless ``impl='givaro'``;
      controls the way elements are printed to the user:

      - 'log': repr is
        :meth:`~sage.rings.finite_rings.element_givaro.FiniteField_givaroElement.log_repr()`

      - 'int': repr is
        :meth:`~sage.rings.finite_rings.element_givaro.FiniteField_givaroElement.int_repr()`

      - 'poly': repr is
        :meth:`~sage.rings.finite_rings.element_givaro.FiniteField_givaroElement.poly_repr()`

    - ``check_irreducible`` -- verify that the polynomial modulus is
      irreducible

    - ``proof`` -- boolean (default: ``True``); if ``True``, use provable
      primality test; otherwise only use pseudoprimality test

    ALIAS: You can also use ``GF`` instead of ``FiniteField`` -- they
    are identical.

    EXAMPLES::

        sage: k.<a> = FiniteField(9); k
        Finite Field in a of size 3^2
        sage: parent(a)
        Finite Field in a of size 3^2
        sage: charpoly(a, 'y')
        y^2 + 2*y + 2

    We illustrate the proof flag.  The following example would hang
    for a very long time if we didn't use ``proof=False``.

    .. NOTE::

        Magma only supports ``proof=False`` for making finite fields,
        so falsely appears to be faster than Sage -- see :issue:`10975`.

    ::

        sage: k = FiniteField(10^1000 + 453, proof=False)
        sage: k = FiniteField((10^1000 + 453)^2, 'a', proof=False)      # long time -- about 5 seconds

    ::

        sage: F.<x> = GF(5)[]
        sage: K.<a> = GF(5**5, name='a', modulus=x^5 - x +1 )
        sage: f = K.modulus(); f
        x^5 + 4*x + 1
        sage: type(f)
         <class 'sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint'>

    By default, the given generator is not guaranteed to be primitive
    (a generator of the multiplicative group), use
    ``modulus="primitive"`` if you need this::

        sage: K.<a> = GF(5^45)
        sage: a.multiplicative_order()
        7105427357601001858711242675781
        sage: a.is_square()
        True
        sage: K.<b> = GF(5^45, modulus='primitive')
        sage: b.multiplicative_order()
        28421709430404007434844970703124

    The modulus must be irreducible::

        sage: K.<a> = GF(5**5, name='a', modulus=x^5 - x)
        Traceback (most recent call last):
        ...
        ValueError: finite field modulus must be irreducible but it is not

    You can't accidentally fool the constructor into thinking the
    modulus is irreducible when it is not, since it actually tests
    irreducibility modulo `p`.  Also, the modulus has to be of the
    right degree (this is always checked)::

        sage: F.<x> = QQ[]
        sage: factor(x^5 + 2)
        x^5 + 2
        sage: K.<a> = GF(5^5, modulus=x^5 + 2)
        Traceback (most recent call last):
        ...
        ValueError: finite field modulus must be irreducible but it is not
        sage: K.<a> = GF(5^5, modulus=x^3 + 3*x + 3, check_irreducible=False)
        Traceback (most recent call last):
        ...
        ValueError: the degree of the modulus does not equal the degree of the field

    Any type which can be converted to the polynomial ring `GF(p)[x]`
    is accepted as modulus::

        sage: K.<a> = GF(13^3, modulus=[1,0,0,2])
        sage: K.<a> = GF(13^10, modulus=pari("ffinit(13,10)"))
        sage: var('x')
        x
        sage: K.<a> = GF(13^2, modulus=x^2 - 2)
        sage: K.<a> = GF(13^2, modulus=sin(x))
        Traceback (most recent call last):
        ...
        TypeError: self must be a numeric expression

    If you wish to live dangerously, you can tell the constructor not
    to test irreducibility using ``check_irreducible=False``, but this
    can easily lead to crashes and hangs -- so do not do it unless you
    know that the modulus really is irreducible!

    ::

        sage: K.<a> = GF(5**2, name='a', modulus=x^2 + 2, check_irreducible=False)

    Even for prime fields, you can specify a modulus. This will not
    change how Sage computes in this field, but it will change the
    result of the :meth:`modulus` and :meth:`gen` methods::

        sage: k.<a> = GF(5, modulus='primitive')
        sage: k.modulus()
        x + 3
        sage: a
        2

    The order of a finite field must be a prime power::

        sage: GF(1)
        Traceback (most recent call last):
        ...
        ValueError: the order of a finite field must be at least 2
        sage: GF(100)
        Traceback (most recent call last):
        ...
        ValueError: the order of a finite field must be a prime power

    Finite fields with explicit random modulus are not cached::

        sage: k.<a> = GF(5**10, modulus='random')
        sage: n.<a> = GF(5**10, modulus='random')
        sage: while k.modulus() == n.modulus():
        ....:     n.<a> = GF(5**10, modulus='random')
        sage: n is k
        False
        sage: GF(5**10, 'a') is GF(5**10, 'a')
        True

    We check that various ways of creating the same finite field yield
    the same object, which is cached::

        sage: K = GF(7, 'a')
        sage: L = GF(7, 'b')
        sage: K is L           # name is ignored for prime fields
        True
        sage: K is GF(7, modulus=K.modulus())
        True
        sage: K = GF(4,'a'); K.modulus()
        x^2 + x + 1
        sage: L = GF(4,'a', K.modulus())
        sage: K is L
        True
        sage: M = GF(4,'a', K.modulus().change_variable_name('y'))
        sage: K is M
        True

    You may print finite field elements as integers. This currently
    only works if the order of field is `<2^{16}`, though::

        sage: k.<a> = GF(2^8, repr='int')
        sage: a
        2

    The following demonstrate coercions for finite fields using Conway
    polynomials::

        sage: k = GF(5^2); a = k.gen()
        sage: l = GF(5^5); b = l.gen()
        sage: a + b
        3*z10^5 + z10^4 + z10^2 + 3*z10 + 1

    Note that embeddings are compatible in lattices of such finite
    fields::

        sage: m = GF(5^3); c = m.gen()
        sage: (a+b)+c == a+(b+c)
        True
        sage: (a*b)*c == a*(b*c)
        True
        sage: from sage.categories.pushout import pushout
        sage: n = pushout(k, l)
        sage: o = pushout(l, m)
        sage: q = pushout(n, o)
        sage: q(o(b)) == q(n(b))
        True

    Another check that embeddings are defined properly::

        sage: k = GF(3**10)
        sage: l = GF(3**20)
        sage: l(k.gen()**10) == l(k.gen())**10
        True

    Using pseudo-Conway polynomials is slow for highly
    composite extension degrees::

        sage: k = GF(3^120)  # long time (about 3 seconds)
        sage: GF(3^40).gen().minimal_polynomial()(k.gen()^((3^120-1)/(3^40-1)))  # long time (because of previous line)
        0

    Before :issue:`17569`, the boolean keyword argument ``conway``
    was required when creating finite fields without a variable
    name.  This keyword argument is now removed (:issue:`21433`).
    You can still pass in ``prefix`` as an argument, which has the
    effect of changing the variable name of the algebraic closure::

        sage: K = GF(3^10, prefix='w'); L = GF(3^10); K is L
        False
        sage: K.variable_name(), L.variable_name()
        ('w10', 'z10')
        sage: list(K.polynomial()) == list(L.polynomial())
        True

    TESTS:

    Check that :issue:`16934` has been fixed::

        sage: k1.<a> = GF(17^14, impl='pari')
        sage: _ = a/2
        sage: k2.<a> = GF(17^14, impl='pari')
        sage: k1 is k2
        True

    Check that :issue:`21433` has been fixed::

        sage: K = GF(5^2)
        sage: L = GF(5^4)
        sage: from sage.categories.pushout import pushout
        sage: pushout(K,L) is L
        True

    Check that :issue:`25182` has been fixed::

        sage: GF(next_prime(2^63)^6)
        Finite Field in z6 of size 9223372036854775837^6

    Check that :issue:`31547` has been fixed::

        sage: q=2**152
        sage: GF(q,'a',modulus='primitive') == GF(q,'a',modulus='primitive')
        True

    Check that a number for ``name`` is interpreted as the degree instead::

        sage: GF(5, 2)
        Finite Field in z2 of size 5^2
        sage: GF((5, 2))
        Finite Field in z2 of size 5^2
        sage: GF(5, 2) and GF((5, 2))
        True
    """
    def __init__(self, *args, **kwds):
        """
        Initialization.

        EXAMPLES::

            sage: TestSuite(GF).run()
        """
        self._modulus_cache = defaultdict(dict)
        super().__init__(*args, **kwds)

    def create_key_and_extra_args(self, order, name=None, modulus=None, names=None,
                                  impl=None, proof=None,
                                  check_prime=True, check_irreducible=True,
                                  prefix=None, repr=None, elem_cache=None,
                                  **kwds):
        """
        EXAMPLES::

            sage: GF.create_key_and_extra_args(9, 'a')                                  # needs sage.libs.linbox
            ((9, ('a',), x^2 + 2*x + 2, 'givaro', 3, 2, True, None, 'poly', True, True, True), {})

        The order `q` can also be given as a pair `(p,n)`::

            sage: GF.create_key_and_extra_args((3, 2), 'a')                             # needs sage.libs.linbox
            ((9, ('a',), x^2 + 2*x + 2, 'givaro', 3, 2, True, None, 'poly', True, True, True), {})

        We do not take invalid keyword arguments and raise a value error
        to better ensure uniqueness::

            sage: GF.create_key_and_extra_args(9, 'a', foo='value')
            Traceback (most recent call last):
            ...
            TypeError: ...create_key_and_extra_args() got an unexpected keyword argument 'foo'

        Moreover, ``repr`` and ``elem_cache`` are ignored when not
        using givaro::

            sage: GF.create_key_and_extra_args(16, 'a', impl='ntl', repr='poly')        # needs sage.libs.ntl
            ((16, ('a',), x^4 + x + 1, 'ntl', 2, 4, True, None, None, None, True, True), {})
            sage: GF.create_key_and_extra_args(16, 'a', impl='ntl', elem_cache=False)   # needs sage.libs.ntl
            ((16, ('a',), x^4 + x + 1, 'ntl', 2, 4, True, None, None, None, True, True), {})
            sage: GF(16, impl='ntl') is GF(16, impl='ntl', repr='foo')                  # needs sage.libs.ntl
            True

        We handle extra arguments for the givaro finite field and
        create unique objects for their defaults::

            sage: GF(25, impl='givaro') is GF(25, impl='givaro', repr='poly')           # needs sage.libs.linbox
            True
            sage: GF(25, impl='givaro') is GF(25, impl='givaro', elem_cache=True)       # needs sage.libs.linbox
            True
            sage: GF(625, impl='givaro') is GF(625, impl='givaro', elem_cache=False)    # needs sage.libs.linbox
            True

        We explicitly take ``structure``, ``implementation`` and ``prec`` attributes
        for compatibility with :class:`~sage.categories.pushout.AlgebraicExtensionFunctor`
        but we ignore them as they are not used, see :issue:`21433`::

            sage: GF.create_key_and_extra_args(9, 'a', structure=None)                  # needs sage.libs.linbox
            ((9, ('a',), x^2 + 2*x + 2, 'givaro', 3, 2, True, None, 'poly', True, True, True), {})

        TESTS::

            sage: GF((6, 1), 'a')       # implicit doctest
            Traceback (most recent call last):
            ...
            ValueError: the order of a finite field must be a prime power

            sage: GF((9, 1), 'a')       # implicit doctest
            Traceback (most recent call last):
            ...
            ValueError: the order of a finite field must be a prime power

            sage: GF((5, 0), 'a')       # implicit doctest
            Traceback (most recent call last):
            ...
            ValueError: the order of a finite field must be a prime power

            sage: GF((3, 2, 1), 'a')    # implicit doctest
            Traceback (most recent call last):
            ...
            ValueError: wrong input for finite field constructor
            sage: GF(5, '^')
            sage: GF((5, 2), 1)
            sage: GF((5, 1), 3)
            sage: GF((5, 2), 3)
            sage: GF(25, 2)
            sage: GF((25, 2))
        """
        import sage.arith.all

        for key, val in kwds.items():
            if key not in ['structure', 'implementation', 'prec', 'embedding', 'latex_names']:
                raise TypeError("create_key_and_extra_args() got an unexpected keyword argument '%s'" % key)
            if not (val is None or isinstance(val, list) and all(c is None for c in val)):
                raise NotImplementedError("ring extension with prescribed %s is not implemented" % key)

        from sage.structure.proof.all import WithProof, arithmetic
        if proof is None:
            proof = arithmetic()
        with WithProof('arithmetic', proof):
            if isinstance(order, tuple):
                if len(order) != 2:
                    raise ValueError('wrong input for finite field constructor')
                p, n = map(Integer, order)
                if p < 2 or n < 1:
                    raise ValueError("the order of a finite field must be a prime power")
                order = p**n
            else:
                order = Integer(order)
                if order < 2:
                    raise ValueError("the order of a finite field must be at least 2")
                if isinstance(name, (int, Integer)):
                    p, n = order, name
                    name = None
                else:
                    p, n = order.perfect_power()
            # at this point, order = p**n
            # note that we haven't tested p for primality

            if n == 1:
                if impl is None:
                    impl = 'modn'
                if name is not None:
                    certify_names((name,))
                name = ('x',)  # Ignore name
                # Every polynomial of degree 1 is irreducible
                check_irreducible = False
            else:
                if names is not None:
                    name = names
                if name is None:
                    if prefix is None:
                        prefix = 'z'
                    name = prefix + str(n)
                    if modulus is not None:
                        raise ValueError("no modulus may be specified if variable name not given")
                    # Fpbar will have a strong reference, since algebraic_closure caches its results,
                    # and the coefficients of modulus lie in GF(p)
                    Fpbar = GF(p).algebraic_closure(prefix)
                    # This will give a Conway polynomial if p,n is small enough to be in the database
                    # and a pseudo-Conway polynomial if it's not.
                    modulus = Fpbar._get_polynomial(n)
                    check_irreducible = False
                name = normalize_names(1, name)

                if impl is None:
                    if order < zech_log_bound and FiniteField_givaro is not None:
                        impl = 'givaro'
                    elif p == 2 and FiniteField_ntl_gf2e is not None:
                        impl = 'ntl'
                    else:
                        impl = 'pari_ffelt'

            # Determine modulus.
            # For the 'modn' implementation, we use the following
            # optimization which we also need to avoid an infinite loop:
            # a modulus of None is a shorthand for x-1.
            if modulus is not None or impl != 'modn':
                from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                R = PolynomialRing(FiniteField(p), 'x')
                if modulus is None:
                    modulus = R.irreducible_element(n)
                if isinstance(modulus, str):
                    # A string specifies an algorithm to find a suitable modulus.
                    if modulus != "random" and modulus in self._modulus_cache[order]:
                        modulus = self._modulus_cache[order][modulus]
                    else:
                        self._modulus_cache[order][modulus] = modulus = R.irreducible_element(n, algorithm=modulus)
                else:
                    if isinstance(modulus, Polynomial):
                        modulus = modulus.change_variable_name('x')
                    modulus = R(modulus).monic()

                    if modulus.degree() != n:
                        raise ValueError("the degree of the modulus does not equal the degree of the field")
                # If modulus is x - 1 for impl="modn", set it to None
                if impl == 'modn' and modulus.list() == [-1,1]:
                    modulus = None
            if modulus is None:
                check_irreducible = False

            # Check extra arguments for givaro and setup their defaults
            # TODO: ntl takes a repr, but ignores it
            if impl == 'givaro':
                if repr is None:
                    repr = 'poly'
                if elem_cache is None:
                    elem_cache = (order < 500)
            else:
                # This has the effect of ignoring these keywords
                repr = None
                elem_cache = None

            return (order, name, modulus, impl, p, n, proof, prefix, repr, elem_cache, check_prime, check_irreducible), {}

    def create_object(self, version, key, **kwds):
        """
        EXAMPLES::

            sage: K = GF(19) # indirect doctest
            sage: TestSuite(K).run()

        We try to create finite fields with various implementations::

            sage: k = GF(2, impl='modn')
            sage: k = GF(2, impl='givaro')                                              # needs sage.libs.linbox
            sage: k = GF(2, impl='ntl')                                                 # needs sage.libs.ntl
            sage: k = GF(2, impl='pari')
            Traceback (most recent call last):
            ...
            ValueError: the degree must be at least 2
            sage: k = GF(2, impl='supercalifragilisticexpialidocious')
            Traceback (most recent call last):
            ...
            ValueError: no such finite field implementation: 'supercalifragilisticexpialidocious'
            sage: k.<a> = GF(2^15, impl='modn')
            Traceback (most recent call last):
            ...
            ValueError: the 'modn' implementation requires a prime order
            sage: k.<a> = GF(2^15, impl='givaro')                                       # needs sage.libs.linbox
            sage: k.<a> = GF(2^15, impl='ntl')                                          # needs sage.libs.ntl
            sage: k.<a> = GF(2^15, impl='pari')
            sage: k.<a> = GF(3^60, impl='modn')
            Traceback (most recent call last):
            ...
            ValueError: the 'modn' implementation requires a prime order
            sage: k.<a> = GF(3^60, impl='givaro')                                       # needs sage.libs.linbox
            Traceback (most recent call last):
            ...
            ValueError: q must be < 2^16
            sage: k.<a> = GF(3^60, impl='ntl')                                          # needs sage.libs.ntl
            Traceback (most recent call last):
            ...
            ValueError: q must be a 2-power
            sage: k.<a> = GF(3^60, impl='pari')
        """
        # IMPORTANT!  If you add a new class to the list of classes
        # that get cached by this factor object, then you *must* add
        # the following method to that class in order to fully support
        # pickling:
        #
        #     def __reduce__(self):   # and include good doctests, please!
        #         return self._factory_data[0].reduce_data(self)
        #
        # This is not in the base class for finite fields, since some finite
        # fields need not be created using this factory object, e.g., residue
        # class fields.

        if len(key) == 5:
            # for backward compatibility of pickles (see trac 10975).
            order, name, modulus, impl, _ = key
            p, n = Integer(order).factor()[0]
            proof = True
            prefix = kwds.get('prefix', None)
            # We can set the defaults here to be those for givaro
            #   as they are otherwise ignored
            repr = 'poly'
            elem_cache = (order < 500)
            check_prime = check_irreducible = False
        elif len(key) == 8:
            # For backward compatibility of pickles (see trac #21433)
            order, name, modulus, impl, _, p, n, proof = key
            prefix = kwds.get('prefix', None)
            # We can set the defaults here to be those for givaro
            #   as they are otherwise ignored
            repr = kwds.get('repr', 'poly')
            elem_cache = kwds.get('elem_cache', (order < 500))
            check_prime = check_irreducible = False
        elif len(key) == 10:
            order, name, modulus, impl, p, n, proof, prefix, repr, elem_cache = key
            check_prime = check_irreducible = False
        else:
            order, name, modulus, impl, p, n, proof, prefix, repr, elem_cache, check_prime, check_irreducible = key

        from sage.structure.proof.all import WithProof
        with WithProof('arithmetic', proof):
            if check_prime and not p.is_prime():
                raise ValueError("the order of a finite field must be a prime power")
            if check_irreducible and not modulus.is_irreducible():
                raise ValueError("finite field modulus must be irreducible but it is not")

        if impl == 'modn':
            if n != 1:
                raise ValueError("the 'modn' implementation requires a prime order")
            from .finite_field_prime_modn import FiniteField_prime_modn
            # Using a check option here is probably a worthwhile
            # compromise since this constructor is simple and used a
            # huge amount.
            K = FiniteField_prime_modn(order, check=False, modulus=modulus)
        else:
            # We have to do this with block so that the finite field
            # constructors below will use the proof flag that was
            # passed in when checking for primality, factoring, etc.
            # Otherwise, we would have to complicate all of their
            # constructors with check options.
            with WithProof('arithmetic', proof):
                if impl == 'givaro':
                    K = FiniteField_givaro(order, name, modulus, repr, elem_cache)
                elif impl == 'ntl':
                    K = FiniteField_ntl_gf2e(order, name, modulus)
                elif impl == 'pari_ffelt' or impl == 'pari':
                    from .finite_field_pari_ffelt import FiniteField_pari_ffelt
                    K = FiniteField_pari_ffelt(p, modulus, name)
                else:
                    raise ValueError("no such finite field implementation: %r" % impl)

            # Temporary; see create_key_and_extra_args() above.
            if prefix is not None:
                K._prefix = prefix

        return K


GF = FiniteField = FiniteFieldFactory("FiniteField")


def is_PrimeFiniteField(x):
    """
    Return ``True`` if ``x`` is a prime finite field.

    This function is deprecated.

    EXAMPLES::

        sage: from sage.rings.finite_rings.finite_field_constructor import is_PrimeFiniteField
        sage: is_PrimeFiniteField(QQ)
        doctest:...: DeprecationWarning: the function is_PrimeFiniteField is deprecated; use isinstance(x, sage.rings.finite_rings.finite_field_base.FiniteField) and x.is_prime_field() instead
        See https://github.com/sagemath/sage/issues/32664 for details.
        False
        sage: is_PrimeFiniteField(GF(7))
        True
        sage: is_PrimeFiniteField(GF(7^2, 'a'))
        False
        sage: is_PrimeFiniteField(GF(next_prime(10^90, proof=False)))
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(32664, "the function is_PrimeFiniteField is deprecated; use isinstance(x, sage.rings.finite_rings.finite_field_base.FiniteField) and x.is_prime_field() instead")

    from .finite_field_prime_modn import FiniteField_prime_modn
    from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_generic

    return isinstance(x, FiniteField_prime_modn) or \
           (isinstance(x, FiniteField_generic) and x.degree() == 1)


zech_log_bound = 2**16
