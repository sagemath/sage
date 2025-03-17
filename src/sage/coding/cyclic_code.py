# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Cyclic code

Let `F` be a field. A `[n, k]` code `C` over `F` is called cyclic if every
cyclic shift of a codeword is also a codeword [Rot2006]_:

    .. MATH::

        \forall c \in C,
        c = (c_{0}, c_{1}, \dots , c_{n-1}) \in C
        \Rightarrow (c_{n-1}, c_{0}, \dots , c_{n-2}) \in C

Let `c = (c_0, c_1, \dots, c_{n-1})` be a codeword of `C`.
This codeword can be seen as a polynomial over `F_q[x]` as follows:
`\Sigma_{i=0}^{n-1} c_i x^i`.
There is a unique monic polynomial `g(x)` such that for every
`c(x) \in F_q[x]` of degree less than `n-1`, we have
`c(x) \in C \Leftrightarrow g(x) | c(x)`.
This polynomial is called the generator polynomial of `C`.

For now, only single-root cyclic codes (i.e. whose length `n` and field order
`q` are coprimes) are implemented.
"""

# *****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#                     2016 Julien Lavauzelle <julien.lavauzelle@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from .linear_code import (AbstractLinearCode,
                          LinearCodeSyndromeDecoder,
                          LinearCodeNearestNeighborDecoder)
from .encoder import Encoder
from .decoder import Decoder
from copy import copy
from sage.rings.integer import Integer
from sage.categories.homset import Hom
from sage.arith.misc import gcd
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as Zmod


def find_generator_polynomial(code, check=True):
    r"""
    Return a possible generator polynomial for ``code``.

    If the code is cyclic, the generator polynomial is the gcd of all the
    polynomial forms of the codewords. Conversely, if this gcd exactly
    generates the code ``code``, then ``code`` is cyclic.

    If ``check`` is set to ``True``, then it also checks that the code is
    indeed cyclic. Otherwise it doesn't.

    INPUT:

    - ``code`` -- a linear code

    - ``check`` -- whether the cyclicity should be checked

    OUTPUT:

    - the generator polynomial of ``code`` (if the code is cyclic).

    EXAMPLES::

        sage: from sage.coding.cyclic_code import find_generator_polynomial
        sage: C = codes.GeneralizedReedSolomonCode(GF(8, 'a').list()[1:], 4)
        sage: find_generator_polynomial(C)
        x^3 + (a^2 + 1)*x^2 + a*x + a^2 + 1
    """
    G = code.generator_matrix()
    F = code.base_ring()
    R = F['x']
    g = gcd(R(row.list()) for row in G)

    if check:
        n = code.length()
        k = code.dimension()
        if (g.degree() != n - k):
            raise ValueError("The code is not cyclic.")
        c = _to_complete_list(g, n)
        if any(vector(c[i:] + c[:i]) not in code for i in range(n)):
            raise ValueError("The code is not cyclic.")

    return g.monic()


def _to_complete_list(poly, length):
    r"""
    Return the vector of length exactly ``length`` corresponding to the
    coefficients of the provided polynomial. If needed, zeros are added.

    INPUT:

    - ``poly`` -- a polynomial

    - ``length`` -- integer

    OUTPUT: the list of coefficients

    EXAMPLES::

        sage: R = PolynomialRing(GF(2), 'X')
        sage: X = R.gen()
        sage: poly = X**4 + X + 1
        sage: sage.coding.cyclic_code._to_complete_list(poly, 7)
        [1, 1, 0, 0, 1, 0, 0]
    """
    L = poly.coefficients(sparse=False)
    return L + [poly.base_ring().zero()] * (length - len(L))


def bch_bound(n, D, arithmetic=False):
    r"""
    Return the BCH bound obtained for a cyclic code of length ``n`` and
    defining set ``D``.

    Consider a cyclic code `C`, with defining set `D`, length `n`, and minimum
    distance `d`. We have the following bound, called BCH bound, on `d`:
    `d \geq \delta + 1`, where `\delta` is the length of the longest arithmetic
    sequence (modulo `n`) of elements in `D`.

    That is, if `\exists c, \gcd(c,n) = 1` such that
    `\{l, l+c, \dots, l + (\delta - 1) \times c\} \subseteq D`,
    then `d \geq \delta + 1` [1]

    The BCH bound is often known in the particular case `c = 1`. The user can
    specify by setting ``arithmetic = False``.

    .. NOTE::

        As this is a specific use case of the BCH bound, it is *not* available
        in the global namespace.
        Call it by using ``sage.coding.cyclic_code.bch_bound``. You can also
        load it into the global namespace by typing
        ``from sage.coding.cyclic_code import bch_bound``.

    INPUT:

    - ``n`` -- integer

    - ``D`` -- list of integers

    - ``arithmetic`` -- (default: ``False``) if it is set to ``True``, then it
      computes the BCH bound using the longest arithmetic sequence definition

    OUTPUT:

    - ``(delta + 1, (l, c))`` -- such that ``delta + 1`` is the BCH bound, and
      ``l, c`` are the parameters of the longest arithmetic sequence
      (see below)

    EXAMPLES::

        sage: n = 15
        sage: D = [14,1,2,11,12]
        sage: sage.coding.cyclic_code.bch_bound(n, D)
        (3, (1, 1))

        sage: n = 15
        sage: D = [14,1,2,11,12]
        sage: sage.coding.cyclic_code.bch_bound(n, D, True)
        (4, (2, 12))
    """
    def longest_streak(step):
        max_len = 1
        max_offset = 0
        j = 0
        while j < n:
            h = j
            while isD[h * step % n]:
                h += 1
            if h - j > max_len:
                max_offset = j * step % n
                max_len = h - j
            j = h + 1
        return (max_len, max_offset)

    isD = [0] * n
    for d in D:
        try:
            isD[d] = 1
        except IndexError:
            raise ValueError("%s must contains integers between 0 and %s" %
                             (D, n - 1))
    if 0 not in isD:
        return (n + 1, (1, 0))

    if not arithmetic:
        one_len, offset = longest_streak(1)
        return (one_len + 1, (1, offset))
    else:
        n = Integer(n)
        longest_streak_list = [(longest_streak(step), step)
                               for step in n.coprime_integers(n // 2 + 1)
                               if step >= 1]
        (max_len, offset), step = max(longest_streak_list)
        return (max_len + 1, (step, offset))


class CyclicCode(AbstractLinearCode):
    r"""
    Representation of a cyclic code.

    We propose three different ways to create a new :class:`CyclicCode`, either by
    providing:

    - the generator polynomial and the length (1)
    - an existing linear code. In that case, a generator polynomial will be
      computed from the provided linear code's parameters (2)
    - (a subset of) the defining set of the cyclic code (3)

    For now, only single-root cyclic codes are implemented. That is, only
    cyclic codes such that its length `n` and field order `q` are coprimes.

    Depending on which behaviour you want, you need to specify the names of the
    arguments to :class:`CyclicCode`. See EXAMPLES section below for details.

    INPUT:

    - ``generator_pol`` -- (default: ``None``) the generator polynomial
      of ``self``; that is, the highest-degree monic polynomial which divides
      every polynomial representation of a codeword in ``self``

    - ``length`` -- (default: ``None``) the length of ``self``; it has to be
      bigger than the degree of ``generator_pol``

    - ``code`` -- (default: ``None``) a linear code

    - ``check`` -- boolean (default: ``False``);  whether the cyclicity of
      ``self`` must be checked while finding the generator polynomial. See
      :meth:`find_generator_polynomial` for details.

    - ``D`` -- (default: ``None``) a list of integers between ``0`` and
      ``length-1``, corresponding to (a subset of) the defining set of the code.
      Will be modified if it is not cyclotomic-closed.

    - ``field`` -- (default: ``None``) the base field of ``self``

    - ``primitive_root`` -- (default: ``None``) the primitive root of
      the splitting field which contains the roots of the generator polynomial.
      It has to be of multiplicative order ``length`` over this field.
      If the splitting field is not ``field``, it also have to be a polynomial
      in ``zx``, where ``x`` is the degree of the extension over the prime
      field. For instance, over ``GF(16)``, it must be a polynomial in ``z4``.

    EXAMPLES:

    We can construct a :class:`CyclicCode` object using three different methods.
    First (1), we provide a generator polynomial and a code length::

        sage: F.<x> = GF(2)[]
        sage: n = 7
        sage: g = x ** 3 + x + 1
        sage: C = codes.CyclicCode(generator_pol=g, length=n)
        sage: C
        [7, 4] Cyclic Code over GF(2)

    We can also provide a code (2). In that case, the program will try to
    extract a generator polynomial (see :meth:`find_generator_polynomial`
    for details)::

        sage: C = codes.GeneralizedReedSolomonCode(GF(8, 'a').list()[1:], 4)
        sage: Cc = codes.CyclicCode(code = C)
        sage: Cc
        [7, 4] Cyclic Code over GF(8)

    Finally, we can give (a subset of) a defining set for the code (3).
    In this case, the generator polynomial will be computed::

        sage: F = GF(16, 'a')
        sage: n = 15
        sage: Cc = codes.CyclicCode(length=n, field=F, D = [1,2])
        sage: Cc
        [15, 13] Cyclic Code over GF(16)
    """

    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, generator_pol=None, length=None, code=None, check=True,
                 D=None, field=None, primitive_root=None):
        r"""
        TESTS:

        If one provides a generator polynomial and a length, we check that
        the length is bigger than the degree of the polynomial::

            sage: F.<x> = GF(2)[]
            sage: n = 2
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            Traceback (most recent call last):
            ...
            ValueError: Only cyclic codes whose length and field order are coprimes are implemented.

        We also check that the polynomial is defined over a finite field::

            sage: F.<x> = RR[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            Traceback (most recent call last):
            ...
            ValueError: The generator polynomial must be defined over a finite field.

        And we check that the generator polynomial divides `x^{n} - 1`,
        where `n` is provided length::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 2 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            Traceback (most recent call last):
            ...
            ValueError: Provided polynomial must divide x^n - 1, where n is the provided length.

        In the case of a code is passed as argument, if it's not possible
        to extract a generator polynomial, an exception is raised::

            sage: G = matrix(GF(2), [[1, 1, 1], [0, 1, 1]])
            sage: C = codes.LinearCode(G)
            sage: Cc = codes.CyclicCode(code=C)
            Traceback (most recent call last):
            ...
            ValueError: The code is not cyclic.

        If the ``primitive_root`` does not lie in an extension of ``field``,
        or is not a primitive `n`-th root of unity, then
        an exception is raised::

            sage: F = GF(2)
            sage: n = 15
            sage: Dset = [1, 2, 4, 8]
            sage: alpha = GF(3).one()
            sage: Cc = codes.CyclicCode(D=Dset, field=F, length=n, primitive_root=alpha)
            Traceback (most recent call last):
            ...
            ValueError: primitive_root must belong to an extension of the base field
            sage: alpha = GF(16).one()
            sage: Cc = codes.CyclicCode(D=Dset, field=F, length=n, primitive_root=alpha)
            Traceback (most recent call last):
            ...
            ValueError: primitive_root must be a primitive n-th root of unity
            sage: alpha = GF(32).gen()
            sage: Cc = codes.CyclicCode(D=Dset, field=F, length=n, primitive_root=alpha)
            Traceback (most recent call last):
            ...
            ValueError: primitive_root must be a primitive n-th root of unity
        """
        # Case (1) : generator polynomial and length are provided.
        if (generator_pol is not None and length is not None and
                code is None and D is None and field is None and
                primitive_root is None):
            F = generator_pol.base_ring()
            if not F.is_finite() or not F.is_field():
                raise ValueError("The generator polynomial must be defined "
                                 "over a finite field.")
            q = F.cardinality()
            if not gcd(length, q) == 1:
                raise ValueError("Only cyclic codes whose length and field "
                                 "order are coprimes are implemented.")
            R = generator_pol.parent()
            deg = generator_pol.degree()
            if not isinstance(length, Integer):
                length = Integer(length)
            if not generator_pol.divides(R.gen() ** length - 1):
                raise ValueError("Provided polynomial must divide x^n - 1, "
                                 "where n is the provided length.")
            self._polynomial_ring = R
            self._dimension = length - deg
            if not generator_pol.is_monic():
                self._generator_polynomial = generator_pol.monic()
            else:
                self._generator_polynomial = generator_pol
            super().__init__(F, length, "Vector", "Syndrome")

        # Case (2) : a code is provided.
        elif (code is not None and
              generator_pol is None and length is None and D is None and
              field is None and primitive_root is None):
            if not isinstance(code, AbstractLinearCode):
                raise ValueError("code must be an AbstractLinearCode")
            F = code.base_ring()
            q = F.cardinality()
            n = code.length()
            if not gcd(n, q) == 1:
                raise ValueError("Only cyclic codes whose length and field "
                                 "order are coprimes are implemented.")
            g = find_generator_polynomial(code, check)
            self._polynomial_ring = g.parent()
            self._generator_polynomial = g
            self._dimension = code.dimension()
            super().__init__(code.base_ring(), n, "Vector", "Syndrome")

        # Case (3) : a defining set, a length and a field are provided
        elif (D is not None and length is not None and field is not None and
              generator_pol is None and code is None):
            F = field
            if not F.is_finite() or not F.is_field():
                raise ValueError("You must provide a finite field.")
            n = length
            q = F.cardinality()
            if not gcd(n, q) == 1:
                raise ValueError("Only cyclic codes whose length and field "
                                 "order are coprimes are implemented.")

            R = F['x']
            s = Zmod(n)(q).multiplicative_order()

            if primitive_root is not None:
                Fsplit = primitive_root.parent()
                try:
                    FE = Hom(F, Fsplit)[0]
                except Exception:
                    raise ValueError("primitive_root must belong to an "
                                     "extension of the base field")
                extension_degree = Fsplit.degree() // F.degree()
                if (extension_degree != s or
                        primitive_root.multiplicative_order() != n):
                    raise ValueError("primitive_root must be a primitive "
                                     "n-th root of unity")
                alpha = primitive_root
            else:
                Fsplit, FE = F.extension(Integer(s), map=True)
                alpha = Fsplit.zeta(n)

            Rsplit = Fsplit['xx']
            xx = Rsplit.gen()

            cosets = Zmod(n).cyclotomic_cosets(q, D)
            pows = [item for l in cosets for item in l]

            sec = FE.section()
            g = R.one()
            for J in cosets:
                pol = Rsplit.one()
                for j in J:
                    pol *= xx - alpha**j
                g *= R([sec(coeff) for coeff in pol])

            # we set class variables
            self._field_embedding = FE
            self._primitive_root = alpha
            self._defining_set = sorted(pows)
            self._polynomial_ring = R
            self._generator_polynomial = g
            self._dimension = n - g.degree()
            super().__init__(F, n, "Vector", "SurroundingBCH")

        else:
            raise AttributeError("You must provide either a code, or a list "
                                 "of powers and the length and the field, or "
                                 "a generator polynomial and the code length")

    def __contains__(self, word):
        r"""
        Return ``True`` if ``word`` belongs to ``self``, ``False`` otherwise.

        INPUT:

        - ``word`` -- the word to test

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: c = vector(GF(2), (1, 1, 1, 0, 0, 1, 0))
            sage: c in C
            True
        """
        g = self.generator_polynomial()
        R = self._polynomial_ring
        return (g.divides(R(word.list())) and word in self.ambient_space())

    def __eq__(self, other):
        r"""
        Test equality between CyclicCode objects.

        INPUT:

        - ``other`` -- the code to test

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C1 = codes.CyclicCode(generator_pol=g, length=n)
            sage: C2 = codes.CyclicCode(generator_pol=g, length=n)
            sage: C1 == C2
            True
        """
        if not isinstance(other, CyclicCode):
            return False
        else:
            R = self._polynomial_ring
            return (self.base_field() == other.base_field() and
                    self.length() == other.length() and
                    self.generator_polynomial() == R(other.generator_polynomial()))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: C
            [7, 4] Cyclic Code over GF(2)
        """
        return ("[%s, %s] Cyclic Code over GF(%s)"
                % (self.length(), self.dimension(),
                   self.base_field().cardinality()))

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: latex(C)
            [7, 4] \textnormal{ Cyclic Code over } \Bold{F}_{2}
        """
        return ("[%s, %s] \\textnormal{ Cyclic Code over } %s"
                % (self.length(), self.dimension(),
                   self.base_field()._latex_()))

    def generator_polynomial(self):
        r"""
        Return the generator polynomial of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: C.generator_polynomial()
            x^3 + x + 1
        """
        return self._generator_polynomial

    def field_embedding(self):
        r"""
        Return the base field embedding into the splitting field.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: C.field_embedding()
            Ring morphism:
              From: Finite Field of size 2
              To:   Finite Field in z3 of size 2^3
              Defn: 1 |--> 1
        """
        if not hasattr(self, "_field_embedding"):
            self.defining_set()
        return self._field_embedding

    def defining_set(self, primitive_root=None):
        r"""
        Return the set of exponents of the roots of ``self``'s generator
        polynomial over the extension field.

        Of course, it depends on the choice of the primitive root of
        the splitting field.

        INPUT:

        - ``primitive_root`` -- (optional) a primitive root of the extension
          field

        EXAMPLES:

        We provide a defining set at construction time::

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: C = codes.CyclicCode(length=n, field=F, D=[1,2])
            sage: C.defining_set()
            [1, 2]

        If the defining set was provided by the user, it might have been
        expanded at construction time. In this case, the expanded defining set
        will be returned::

            sage: C = codes.CyclicCode(length=13, field=F, D=[1, 2])
            sage: C.defining_set()
            [1, 2, 3, 5, 6, 9]

        If a generator polynomial was passed at construction time,
        the defining set is computed using this polynomial::

            sage: R.<x> = F[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: C.defining_set()
            [1, 2, 4]

        Both operations give the same result::

            sage: C1 = codes.CyclicCode(length=n, field=F, D=[1, 2, 4])
            sage: C1.generator_polynomial() == g
            True

        Another one, in a reversed order::

            sage: n = 13
            sage: C1 = codes.CyclicCode(length=n, field=F, D=[1, 2])
            sage: g = C1.generator_polynomial()
            sage: C2 = codes.CyclicCode(generator_pol=g, length=n)
            sage: C1.defining_set() == C2.defining_set()
            True
        """
        if (hasattr(self, "_defining_set") and
                (primitive_root is None or
                 primitive_root == self._primitive_root)):
            return self._defining_set
        else:
            F = self.base_field()
            n = self.length()
            q = F.cardinality()
            g = self.generator_polynomial()

            s = Zmod(n)(q).multiplicative_order()

            if primitive_root is None:
                Fsplit, FE = F.extension(Integer(s), map=True)
                alpha = Fsplit.zeta(n)
            else:
                try:
                    alpha = primitive_root
                    Fsplit = alpha.parent()
                    FE = Hom(Fsplit, F)[0]
                except ValueError:
                    raise ValueError("primitive_root does not belong to the "
                                     "right splitting field")
                if alpha.multiplicative_order() != n:
                    raise ValueError("primitive_root must have multiplicative "
                                     "order equal to the code length")

            Rsplit = Fsplit['xx']
            gsplit = Rsplit([FE(coeff) for coeff in g])
            roots = gsplit.roots(multiplicities=False)
            D = [root.log(alpha) for root in roots]

            self._field_embedding = FE
            self._primitive_root = alpha
            self._defining_set = sorted(D)
            return self._defining_set

    def primitive_root(self):
        r"""
        Return the primitive root of the splitting field that is used
        to build the defining set of the code.

        If it has not been specified by the user, it is set by default with the
        output of the ``zeta`` method of the splitting field.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: C.primitive_root()
            z3

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: a = F.gen()
            sage: Cc = codes.CyclicCode(length=n, field=F, D=[1,2],
            ....:                       primitive_root=a^2 + 1)
            sage: Cc.primitive_root()
            a^2 + 1
        """
        if hasattr(self, "_primitive_root"):
            return self._primitive_root
        else:
            self.defining_set()
            return self._primitive_root

    @cached_method
    def check_polynomial(self):
        r"""
        Return the check polynomial of ``self``.

        Let `C` be a cyclic code of length `n` and `g` its generator
        polynomial. The following: `h = \frac{x^n - 1}{g(x)}` is called `C`'s
        check polynomial.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: h = C.check_polynomial()
            sage: h == (x**n - 1)/C.generator_polynomial()
            True
        """
        R = self._polynomial_ring
        n = self.length()
        self._check_polynomial = (R.gen() ** n - 1) // self.generator_polynomial()
        return self._check_polynomial

    @cached_method
    def parity_check_matrix(self):
        r"""
        Return the parity check matrix of ``self``.

        The parity check matrix of a linear code `C` corresponds to the
        generator matrix of the dual code of `C`.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: C.parity_check_matrix()
            [1 0 1 1 1 0 0]
            [0 1 0 1 1 1 0]
            [0 0 1 0 1 1 1]
        """
        k = self.dimension()
        n = self.length()
        h = self.check_polynomial().reverse()
        l = _to_complete_list(h, n)
        M = matrix([l[-i:] + l[:-i] for i in range(n - k)])
        M.set_immutable()
        return M

    def bch_bound(self, arithmetic=False):
        r"""
        Return the BCH bound of ``self`` which is a bound on ``self``
        minimum distance.

        See :meth:`sage.coding.cyclic_code.bch_bound` for details.

        INPUT:

        - ``arithmetic`` -- (default: ``False``) if it is set to ``True``,
          then it computes the BCH bound using the longest arithmetic sequence
          definition

        OUTPUT:

        - ``(delta + 1, (l, c))`` -- such that ``delta + 1`` is the BCH bound,
          and ``l, c`` are the parameters of the largest arithmetic sequence

        EXAMPLES::

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: D = [14,1,2,11,12]
            sage: C = codes.CyclicCode(field=F, length=n, D = D)
            sage: C.bch_bound()
            (3, (1, 1))

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: D = [14,1,2,11,12]
            sage: C = codes.CyclicCode(field=F, length=n, D = D)
            sage: C.bch_bound(True)
            (4, (2, 12))
        """
        return bch_bound(self.length(), self.defining_set(), arithmetic)

    def surrounding_bch_code(self):
        r"""
        Return the surrounding BCH code of ``self``.

        EXAMPLES::

            sage: C = codes.CyclicCode(field=GF(2), length=63, D=[1, 7, 17])
            sage: C.dimension()
            45
            sage: CC = C.surrounding_bch_code()
            sage: CC
            [63, 51] BCH Code over GF(2) with designed distance 3
            sage: all(r in CC for r in C.generator_matrix())
            True
        """
        from .bch_code import BCHCode
        delta, params = self.bch_bound(arithmetic=True)
        return BCHCode(self.base_field(), self.length(), delta,
                       offset=params[1], jump_size=params[0])


class CyclicCodePolynomialEncoder(Encoder):
    r"""
    An encoder encoding polynomials into codewords.

    Let `C` be a cyclic code over some finite field `F`,
    and let `g` be its generator polynomial.

    This encoder encodes any polynomial `p \in F[x]_{<k}` by computing
    `c = p g` and returning the vector of its coefficients.

    INPUT:

    - ``code`` -- the associated code of this encoder

    EXAMPLES::

        sage: F.<x> = GF(2)[]
        sage: n = 7
        sage: g = x ** 3 + x + 1
        sage: C = codes.CyclicCode(generator_pol=g, length=n)
        sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
        sage: E
        Polynomial-style encoder for [7, 4] Cyclic Code over GF(2)
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E
            Polynomial-style encoder for [7, 4] Cyclic Code over GF(2)
        """
        if not isinstance(code, CyclicCode):
            raise ValueError("code has to be a CyclicCode")
        self._polynomial_ring = code._polynomial_ring
        super().__init__(code)

    def __eq__(self, other):
        r"""
        Test equality between CyclicCodePolynomialEncoder objects.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E1 = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E2 = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E1 == E2
            True
        """
        return (isinstance(other, CyclicCodePolynomialEncoder) and
                self.code() == other.code())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E
            Polynomial-style encoder for [7, 4] Cyclic Code over GF(2)
        """
        return "Polynomial-style encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: latex(E)
            \textnormal{Polynomial-style encoder for }[7, 4] \textnormal{ Cyclic Code over } \Bold{F}_{2}
        """
        return ("\\textnormal{Polynomial-style encoder for }%s" %
                self.code()._latex_())

    def encode(self, p):
        r"""
        Transform `p` into an element of the associated code of ``self``.

        INPUT:

        - ``p`` -- a polynomial from ``self`` message space

        OUTPUT: a codeword in associated code of ``self``

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: m = x ** 2 + 1
            sage: E.encode(m)
            (1, 1, 1, 0, 0, 1, 0)
        """
        C = self.code()
        k = C.dimension()
        n = C.length()
        if p.degree() >= k:
            raise ValueError("Degree of the message must be at most %s" % k - 1)
        res = _to_complete_list(p * C.generator_polynomial(), n)
        return vector(C.base_field(), res)

    def unencode_nocheck(self, c):
        r"""
        Return the message corresponding to ``c``.
        Does not check if ``c`` belongs to the code.

        INPUT:

        - ``c`` -- a vector with the same length as the code

        OUTPUT: an element of the message space

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: c = vector(GF(2), (1, 1, 1, 0, 0, 1, 0))
            sage: E.unencode_nocheck(c)
            x^2 + 1
        """
        R = self.message_space()
        g = self.code().generator_polynomial()
        p = R(c.list())
        return p // g

    def message_space(self):
        r"""
        Return the message space of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E.message_space()                                                     # needs sage.libs.ntl
            Univariate Polynomial Ring in x over Finite Field of size 2 (using GF2X)
        """
        return self._polynomial_ring


class CyclicCodeVectorEncoder(Encoder):
    r"""
    An encoder which can encode vectors into codewords.

    Let `C` be a cyclic code over some finite field `F`,
    and let `g` be its generator polynomial.

    Let `m = (m_1, m_2, \dots, m_k)` be a vector in `F^{k}`.
    This codeword can be seen as a polynomial over `F[x]`, as follows:
    `P_m = \Sigma_{i=0}^{k-1} m_i \times x^i`.

    To encode `m`, this encoder does the multiplication `P_m  g`.

    INPUT:

    - ``code`` -- the associated code of this encoder

    EXAMPLES::

        sage: F.<x> = GF(2)[]
        sage: n = 7
        sage: g = x ** 3 + x + 1
        sage: C = codes.CyclicCode(generator_pol=g, length=n)
        sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
        sage: E
        Vector-style encoder for [7, 4] Cyclic Code over GF(2)
    """

    def __init__(self, code):
        r"""

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E
            Vector-style encoder for [7, 4] Cyclic Code over GF(2)
        """
        if not isinstance(code, CyclicCode):
            raise ValueError("code has to be a CyclicCode")
        self._polynomial_ring = code._polynomial_ring
        super().__init__(code)

    def __eq__(self, other):
        r"""
        Test equality between CyclicCodeVectorEncoder objects.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E1 = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E2 = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E1 == E2
            True
        """
        return (isinstance(other, CyclicCodeVectorEncoder) and
                self.code() == other.code())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E
            Vector-style encoder for [7, 4] Cyclic Code over GF(2)
        """
        return "Vector-style encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: latex(E)
            \textnormal{Vector-style encoder for }[7, 4] \textnormal{ Cyclic Code over } \Bold{F}_{2}
        """
        return ("\\textnormal{Vector-style encoder for }%s" %
                self.code()._latex_())

    def encode(self, m):
        r"""
        Transform `m` into an element of the associated code of ``self``.

        INPUT:

        - ``m`` -- an element from ``self``'s message space

        OUTPUT: a codeword in the associated code of ``self``

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: m = vector(GF(2), (1, 0, 1, 0))
            sage: E.encode(m)
            (1, 1, 1, 0, 0, 1, 0)
        """
        if self.generator_matrix.cache is not None:
            return super().encode(m)

        k = self.code().dimension()
        n = self.code().length()
        F = self.code().base_field()
        R = self._polynomial_ring
        p = R(m.list())
        if p.degree() >= k:
            raise ValueError("Degree of the message must be at most %s" % k - 1)
        res = _to_complete_list(p * self.code().generator_polynomial(), n)
        return vector(F, res)

    def unencode_nocheck(self, c):
        r"""
        Return the message corresponding to ``c``.
        Does not check if ``c`` belongs to the code.

        INPUT:

        - ``c`` -- a vector with the same length as the code

        OUTPUT: an element of the message space

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: c = vector(GF(2), (1, 1, 1, 0, 0, 1, 0))
            sage: E.unencode_nocheck(c)
            (1, 0, 1, 0)
        """

        R = self._polynomial_ring
        g = self.code().generator_polynomial()
        p = R(c.list())
        l = _to_complete_list(p // g, self.message_space().dimension())
        return vector(self.code().base_field(), l)

    @cached_method
    def generator_matrix(self):
        r"""
        Return a generator matrix of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E.generator_matrix()
            [1 1 0 1 0 0 0]
            [0 1 1 0 1 0 0]
            [0 0 1 1 0 1 0]
            [0 0 0 1 1 0 1]
        """
        C = self.code()
        k = C.dimension()
        n = C.length()
        l = _to_complete_list(C.generator_polynomial(), n)
        M = matrix([l[-i:] + l[:-i] for i in range(k)])
        M.set_immutable()
        return M

    def message_space(self):
        r"""
        Return the message space of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol=g, length=n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E.message_space()
            Vector space of dimension 4 over Finite Field of size 2
        """
        return self.code().base_ring() ** self.code().dimension()


class CyclicCodeSurroundingBCHDecoder(Decoder):
    r"""
    A decoder which decodes through the surrounding BCH code of the cyclic
    code.

    INPUT:

    - ``code`` -- the associated code of this decoder

    - ``**kwargs`` -- all extra arguments are forwarded to the BCH decoder

    EXAMPLES::

        sage: C = codes.CyclicCode(field=GF(16), length=15, D=[14, 1, 2, 11, 12])
        sage: D = codes.decoders.CyclicCodeSurroundingBCHDecoder(C)
        sage: D
        Decoder through the surrounding BCH code of the [15, 10] Cyclic Code over GF(16)
    """
    def __init__(self, code, **kwargs):
        r"""

        EXAMPLES::

            sage: C = codes.CyclicCode(field=GF(16), length=15, D=[14, 1, 2, 11, 12])
            sage: D = codes.decoders.CyclicCodeSurroundingBCHDecoder(C)
            sage: D
            Decoder through the surrounding BCH code of the [15, 10] Cyclic Code over GF(16)
        """
        self._bch_code = code.surrounding_bch_code()
        self._bch_decoder = self._bch_code.decoder(**kwargs)
        self._decoder_type = copy(self._bch_decoder.decoder_type())
        super().__init__(code, code.ambient_space(), "Vector")

    def __eq__(self, other):
        r"""
        Test equality between CyclicCodeSurroundingBCHDecoder objects.

        EXAMPLES::

            sage: C = codes.CyclicCode(field=GF(16), length=15, D=[14, 1, 2, 11, 12])
            sage: D1 = C.decoder()
            sage: D2 = C.decoder()
            sage: D1 == D2
            True
        """
        return (isinstance(other, CyclicCodeSurroundingBCHDecoder) and
                self.code() == other.code() and
                self.bch_decoder() == other.bch_decoder())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.CyclicCode(field=GF(16), length=15, D=[14, 1, 2, 11, 12])
            sage: D = codes.decoders.CyclicCodeSurroundingBCHDecoder(C)
            sage: D
            Decoder through the surrounding BCH code of the [15, 10] Cyclic Code over GF(16)
        """
        return ("Decoder through the surrounding BCH code of the %s" %
                self.code())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.CyclicCode(field=GF(16), length=15, D=[14, 1, 2, 11, 12])
            sage: D = codes.decoders.CyclicCodeSurroundingBCHDecoder(C)
            sage: latex(D)
            \textnormal{Decoder through the surrounding BCH code of the }[15, 10] \textnormal{ Cyclic Code over } \Bold{F}_{2^{4}}
        """
        return ("\\textnormal{Decoder through the surrounding BCH code of "
                "the }%s" % self.code()._latex_())

    def bch_code(self):
        r"""
        Return the surrounding BCH code of
        :meth:`sage.coding.encoder.Encoder.code`.

        EXAMPLES::

            sage: C = codes.CyclicCode(field=GF(16), length=15, D=[14, 1, 2, 11, 12])
            sage: D = codes.decoders.CyclicCodeSurroundingBCHDecoder(C)
            sage: D.bch_code()
            [15, 12] BCH Code over GF(16) with designed distance 4
        """
        return self._bch_code

    def bch_decoder(self):
        r"""
        Return the decoder that will be used over the surrounding BCH code.

        EXAMPLES::

            sage: C = codes.CyclicCode(field=GF(16), length=15, D=[14, 1, 2, 11, 12])
            sage: D = codes.decoders.CyclicCodeSurroundingBCHDecoder(C)
            sage: D.bch_decoder()
            Decoder through the underlying GRS code of [15, 12] BCH Code
             over GF(16) with designed distance 4
        """
        return self._bch_decoder

    def decode_to_code(self, y):
        r"""
        Decodes ``r`` to an element in :meth:`sage.coding.encoder.Encoder.code`.

        EXAMPLES::

            sage: F = GF(16, 'a')
            sage: C = codes.CyclicCode(field=F, length=15, D=[14, 1, 2, 11, 12])
            sage: a = F.gen()
            sage: D = codes.decoders.CyclicCodeSurroundingBCHDecoder(C)
            sage: y = vector(F, [0, a^3, a^3 + a^2 + a, 1, a^2 + 1, a^3 + a^2 + 1,
            ....:                a^3 + a^2 + a, a^3 + a^2 + a, a^2 + a, a^2 + 1,
            ....:                a^2 + a + 1, a^3 + 1, a^2, a^3 + a, a^3 + a])
            sage: D.decode_to_code(y) in C
            True
        """
        return self.bch_code().decode_to_code(y)

    def decoding_radius(self):
        r"""
        Return maximal number of errors that ``self`` can decode.

        EXAMPLES::

            sage: C = codes.CyclicCode(field=GF(16), length=15, D=[14, 1, 2, 11, 12])
            sage: D = codes.decoders.CyclicCodeSurroundingBCHDecoder(C)
            sage: D.decoding_radius()
            1
        """
        return self._bch_decoder.decoding_radius()


# ###################### registration ##############################

CyclicCode._registered_encoders["Polynomial"] = CyclicCodePolynomialEncoder
CyclicCode._registered_encoders["Vector"] = CyclicCodeVectorEncoder
CyclicCode._registered_decoders["Syndrome"] = LinearCodeSyndromeDecoder
CyclicCode._registered_decoders["NearestNeighbor"] = LinearCodeNearestNeighborDecoder

CyclicCode._registered_decoders["SurroundingBCH"] = CyclicCodeSurroundingBCHDecoder
CyclicCodeSurroundingBCHDecoder._decoder_type = {"dynamic"}
