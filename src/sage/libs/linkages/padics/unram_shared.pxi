cimport cython


@cython.binding(True)
def frobenius_unram(self, arithmetic=True):
    """
    Return the image of this element under the Frobenius automorphism
    applied to its parent.

    INPUT:

    - ``self`` -- an element of an unramified extension
    - ``arithmetic`` -- whether to apply the arithmetic Frobenius (acting
      by raising to the `p`-th power on the residue field). If ``False`` is
      provided, the image of geometric Frobenius (raising to the `(1/p)`-th
      power on the residue field) will be returned instead.

    EXAMPLES::

        sage: R.<a> = Zq(5^4,3)
        sage: a.frobenius()
        (a^3 + a^2 + 3*a) + (3*a + 1)*5 + (2*a^3 + 2*a^2 + 2*a)*5^2 + O(5^3)
        sage: f = R.defining_polynomial()
        sage: f(a)
        O(5^3)
        sage: f(a.frobenius())
        O(5^3)
        sage: for i in range(4): a = a.frobenius()
        sage: a
        a + O(5^3)

        sage: R.<a> = Zq(5^4,3)
        sage: a.frobenius(arithmetic=False)
        (3*a^3 + 3*a^2 + a) + (a^3 + 4*a^2 + a + 4)*5 + (3*a^2 + 2*a + 3)*5^2 + O(5^3)

        sage: K.<a> = Qq(7^3,4)
        sage: b = (a+1)/7
        sage: c = b.frobenius(); c
        (3*a^2 + 5*a + 1)*7^-1 + (6*a^2 + 6*a + 6) + (4*a^2 + 3*a + 4)*7 + (6*a^2 + a + 6)*7^2 + O(7^3)
        sage: c.frobenius().frobenius()
        (a + 1)*7^-1 + O(7^3)

    An error will be raised if the parent of ``self`` is a ramified extension::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = Qp(5).extension(x^2 - 5)
        sage: a.frobenius()
        Traceback (most recent call last):
        ...
        NotImplementedError: Frobenius automorphism only implemented for unramified extensions

    TESTS:

    We check that :issue:`23575` is resolved::

        sage: x = R.random_element()
        sage: x.frobenius(arithmetic=false).frobenius() == x
        True
    """
    if self == 0:
        return self
    R = self.parent()
    p = R.prime()
    a = R.gen()
    frob_a = R._frob_gen(arithmetic)
    ppow = self.valuation()
    unit = self.unit_part()
    coefs = unit.expansion()
    ans = 0

    # Xavier's implementation based on Horner scheme
    for i in range(R.f()-1, -1, -1):
        update = 0
        for j in range(len(coefs)-1, -1, -1):
            update *= p
            try:
                update += coefs[j][i]
            except IndexError:
                pass
        ans *= frob_a
        ans += update
    return ans << ppow


@cython.binding(True)
def norm_unram(self, base = None):
    r"""
    Return the absolute or relative norm of this element.

    .. WARNING::

        This is not the `p`-adic absolute value.  This is a
        field theoretic norm down to a ground ring.  If you want the
        `p`-adic absolute value, use the ``abs()`` function instead.

    INPUT:

    - ``base`` -- a subfield of the parent `L` of this element; the norm is the
      relative norm from ``L`` to ``base``. Defaults to the absolute norm down
      to `\QQ_p` or `\ZZ_p`.

    EXAMPLES::

        sage: R = ZpCR(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: ((1+2*w)^5).norm()
        1 + 5^2 + O(5^5)
        sage: ((1+2*w)).norm()^5
        1 + 5^2 + O(5^5)

    TESTS::

        sage: R = ZpCA(5,5)
        sage: S.<x> = ZZ[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: ((1+2*w)^5).norm()
        1 + 5^2 + O(5^5)
        sage: ((1+2*w)).norm()^5
        1 + 5^2 + O(5^5)
        sage: R = ZpFM(5,5)
        sage: S.<x> = ZZ[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: ((1+2*w)^5).norm()
        1 + 5^2
        sage: ((1+2*w)).norm()^5
        1 + 5^2

    TESTS:

    Check that :issue:`11586` has been resolved::

        sage: R.<x> = QQ[]
        sage: f = x^2 + 3*x + 1
        sage: M.<a> = Qp(7).extension(f)
        sage: M(7).norm()
        7^2 + O(7^22)
        sage: b = 7*a + 35
        sage: b.norm()
        4*7^2 + 7^3 + O(7^22)
        sage: b*b.frobenius()
        4*7^2 + 7^3 + O(7^22)

    Check that :issue:`31845` is fixed::

        sage: R.<a> = Zq(4)
        sage: (a - a).norm()
        O(2^20)
    """
    if base is not None:
        if base is self.parent():
            return self
        else:
            raise NotImplementedError
    if self._is_exact_zero():
        return self.parent().ground_ring()(0)
    elif self._is_inexact_zero():
        return self.parent().ground_ring()(0, self.valuation())
    if self.valuation() == 0:
        return self.parent().ground_ring()(self.matrix_mod_pn().det())
    else:
        if self.parent().e() == 1:
            norm_of_uniformizer = self.parent().ground_ring().uniformizer_pow(self.parent().degree())
        else:
            norm_of_uniformizer = (-1)**self.parent().degree() * self.parent().defining_polynomial()[0]
        return self.parent().ground_ring()(self.unit_part().matrix_mod_pn().det()) * norm_of_uniformizer**self.valuation()


@cython.binding(True)
def trace_unram(self, base = None):
    r"""
    Return the absolute or relative trace of this element.

    If ``base`` is given then ``base`` must be a subfield of the
    parent `L` of ``self``, in which case the trace is the relative
    trace from `L` to ``base``.

    In all other cases, the trace is the absolute trace down to
    `\QQ_p` or `\ZZ_p`.

    EXAMPLES::

        sage: R = ZpCR(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = (2+3*w)^7
        sage: b = (6+w^3)^5
        sage: a.trace()
        3*5 + 2*5^2 + 3*5^3 + 2*5^4 + O(5^5)
        sage: a.trace() + b.trace()
        4*5 + 5^2 + 5^3 + 2*5^4 + O(5^5)
        sage: (a+b).trace()
        4*5 + 5^2 + 5^3 + 2*5^4 + O(5^5)

    TESTS::

        sage: R = ZpCA(5,5)
        sage: S.<x> = ZZ[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = (2+3*w)^7
        sage: b = (6+w^3)^5
        sage: a.trace()
        3*5 + 2*5^2 + 3*5^3 + 2*5^4 + O(5^5)
        sage: a.trace() + b.trace()
        4*5 + 5^2 + 5^3 + 2*5^4 + O(5^5)
        sage: (a+b).trace()
        4*5 + 5^2 + 5^3 + 2*5^4 + O(5^5)
        sage: R = ZpFM(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: a = (2+3*w)^7
        sage: b = (6+w^3)^5
        sage: a.trace()
        3*5 + 2*5^2 + 3*5^3 + 2*5^4
        sage: a.trace() + b.trace()
        4*5 + 5^2 + 5^3 + 2*5^4
        sage: (a+b).trace()
        4*5 + 5^2 + 5^3 + 2*5^4

    Check that :issue:`31845` is fixed::

        sage: R.<a> = Zq(4)
        sage: (a - a).trace()
        O(2^20)
    """
    if base is not None:
        if base is self.parent():
            return self
        else:
            raise NotImplementedError
    if self._is_exact_zero():
        return self.parent().ground_ring()(0)
    elif self._is_inexact_zero():
        return self.parent().ground_ring()(0, self.precision_absolute())
    if self.valuation() >= 0:
        return self.parent().ground_ring()(self.matrix_mod_pn().trace())
    else:
        shift = -self.valuation()
        return self.parent().ground_ring()((self * self.parent().prime() ** shift).matrix_mod_pn().trace()) / self.parent().prime()**shift
