# -*- coding: utf-8 -*-
r"""
`p`-Selmer groups of number fields

This file contains code to compute `K(S,p)` where

- `K` is a number field
- `S` is a finite set of primes of `K`
- `p` is a prime number

For `m\ge2`, `K(S,m)` is defined to be the finite subgroup of
`K^*/(K^*)^m` consisting of elements represented by `a\in K^*` whose
valuation at all primes not in `S` is a multiple of `m`.  It fits in
the short exact sequence

.. MATH::

 1 \rightarrow O^*_{K,S}/(O^*_{K,S})^m \rightarrow K(S,m) \rightarrow Cl_{K,S}[m] \rightarrow 1

where `O^*_{K,S}` is the group of `S`-units of `K` and `Cl_{K,S}` the
`S`-class group.  When `m=p` is prime, `K(S,p)` is a
finite-dimensional vector space over `GF(p)`.  Its generators come
from three sources: units (modulo `p`'th powers); generators of the
`p`'th powers of ideals which are not principal but whose `p`'the
powers are principal; and generators coming from the prime ideals in
`S`.

The main function here is :meth:`pSelmerGroup`.  This will not
normally be used by users, who instead will access it through a method
of the NumberField class.


AUTHORS:

- John Cremona (2005-2021)

"""

# ****************************************************************************
#       Copyright (C) 2021 John Cremona <john.cremona@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.rational_field import QQ
from sage.misc.misc_c import prod

# A utility function to allow the same code to be used over QQ and
# over number fields:

def _ideal_generator(I):
    r"""
    Return the generator of a principal ideal.

    INPUT:

    - ``I`` (fractional ideal or integer) -- either a fractional ideal of a
      number field, which must be principal, or a rational integer.

    OUTPUT:

    A generator of `I` when `I` is a principal ideal, else `I` itself.

    EXAMPLES::

        sage: from sage.rings.number_field.selmer_group import _ideal_generator
        sage: _ideal_generator(5)
        5

        sage: K.<a> = QuadraticField(-11)
        sage: [_ideal_generator(K.prime_above(p)) for p in primes(25)]
        [2, 1/2*a - 1/2, -1/2*a - 3/2, 7, -a, 13, 17, 19, 1/2*a + 9/2]

    """
    try:
        return I.gens_reduced()[0]
    except AttributeError:
        return I.abs()

def _coords_in_C_p(I, C, p):
    r"""
    Return coordinates of the ideal ``I`` with respect to a basis of
    the ``p``-torsion of the ideal class group ``C``.

    INPUT:

    - ``I`` (ideal) -- a fractional ideal of a number field ``K``,
      whose ``p``'th power is principal.

    - ``C`` (class group) -- the ideal class group of ``K``.

    - ``p`` (prime) -- a prime number.

    OUTPUT:

    The coordinates of the ideal class `[I]` in the `p`-torsion
    subgroup `C[p]`. An error is raised if `I^p` is not principal.

    ALGORITHM:

    Find the coordinates of `[I]` with respect to generators of `C` as
    an abelian group, check that coordidates are 0 in cyclic factors
    of order prime to `p`, and return the list of `c/(n/p)` (mod `p`)
    for coordinates `c` for each cyclic factor of order `n` which is a
    multiple of `p`.

    EXAMPLES::

        sage: from sage.rings.number_field.selmer_group import _coords_in_C_p
        sage: K.<a> = QuadraticField(-5)
        sage: C = K.class_group()
        sage: C.order()
        2
        sage: P = K.prime_above(2)
        sage: C(P).order()
        2
        sage: _coords_in_C_p(P,C,2)
        [1]
        sage: _coords_in_C_p(P,C,3)
        Traceback (most recent call last):
        ...
        ValueError: The 3rd power of Fractional ideal (2, a + 1) is not principal

    """
    cyclic_orders = C.gens_orders()
    non_p_indices = [i           for i,n in enumerate(cyclic_orders) if not p.divides(n)]
    p_indices     = [(i, n // p) for i,n in enumerate(cyclic_orders) if p.divides(n)]

    coords = C(I).exponents()
    if all(coords[i] == 0 for i in non_p_indices) and all(coords[i] % n == 0 for i, n in p_indices):
        return [(coords[i] // n) % p for i, n in p_indices]
    raise ValueError("The {} power of {} is not principal".format(p.ordinal_str(),I))

def _coords_in_C_mod_p(I,C,p):
    r"""
    Return coordinates of the ideal ``I`` with respect to a basis of
    the ``p``-cotorsion of the ideal class group ``C``.

    INPUT:

    - ``I`` (ideal) -- a fractional ideal of a number field ``K``.

    - ``C`` (class group) -- the ideal class group of ``K``.

    - ``p`` (prime) -- a prime number.

    OUTPUT:

    The coordinates of the ideal class `[I]` in the `p`-cotorsion group `C/C^p`.

    ALGORITHM:

    Find the coordinates of `[I]` with respect to generators of `C` as
    an abelian group, and return the list of `c` (mod `p`)
    for coordinates `c` for each cyclic factor of order `n` which is a
    multiple of `p`.

    EXAMPLES::

        sage: from sage.rings.number_field.selmer_group import _coords_in_C_mod_p
        sage: K.<a> = QuadraticField(-5)
        sage: C = K.class_group()
        sage: [_coords_in_C_mod_p(K.prime_above(p), C, 2) for p in primes(25)]
        [[1], [1], [0], [1], [0], [0], [0], [0], [1]]

    An example where the class group has two primary components, one
    of which is not cyclic::

        sage: from sage.rings.number_field.selmer_group import _coords_in_C_mod_p
        sage: K.<a> = NumberField(x^2 - x + 58)
        sage: C = K.class_group()
        sage: C.gens_orders()
        (6, 2)
        sage: [_coords_in_C_mod_p(K.prime_above(p), C, 2) for p in primes(25)]
        [[1, 0], [1, 1], [1, 1], [0, 1], [1, 0], [0, 1], [0, 0], [0, 1], [0, 0]]
        sage: [_coords_in_C_mod_p(K.prime_above(p), C, 3) for p in primes(25)]
        [[2], [0], [1], [0], [0], [1], [0], [2], [0]]

    """
    cyclic_orders = C.gens_orders()
    p_indices = [i for i, n in enumerate(cyclic_orders) if p.divides(n)]
    coords = C(I).exponents()
    return [coords[i] % p for i in p_indices]

def _root_ideal(I, C, p):
    r"""
    Return a ``p``'th root of an ideal with respect to the class group.

    INPUT:

    - ``I`` (ideal) -- a fractional ideal of a number field ``K``,
      whose ideal class is a ``p``'th power.

    - ``C`` (class group) -- the ideal class group of ``K``.

    - ``p`` (prime) -- a prime number.

    OUTPUT:

    An ideal `J` such that `J^p` is in the same ideal class as `I`.

    EXAMPLES::

        sage: from sage.rings.number_field.selmer_group import _root_ideal
        sage: K.<a> = NumberField(x^2 - x + 58)
        sage: C = K.class_group()
        sage: cyclic_gens   = C.gens_ideals()
        sage: [C(I).order() for I in cyclic_gens]
        [6, 2]
        sage: C.gens_orders()
        (6, 2)
        sage: I = cyclic_gens[0]^2
        sage: J = _root_ideal(I, C, 2)
        sage: C(J^2) == C(I)
        True
        sage: I = cyclic_gens[0]^3
        sage: J = _root_ideal(I, C, 3)
        sage: C(J^3) == C(I)
        True

    """
    cyclic_orders = C.gens_orders()
    cyclic_gens   = C.gens_ideals()
    coords = C(I).exponents()

    # In the next line, e=(ci/p)%n should satisfy p*e=ci (mod n): we
    # are dividing the coordinate vector by p in the appropriate sense

    if not all(p.divides(ci) for ci, n in zip(coords, cyclic_orders) if p.divides(n)):
        raise ValueError("The ideal class of {} is not a {} power".format(I,p.ordinal_str()))

    w = [ci // p if p.divides(n) else (ci / p) % n for ci, n in zip(coords, cyclic_orders)]

    return prod([gen ** wi for wi, gen in zip(w, cyclic_gens)], C.number_field().ideal(1))

def coords_in_U_mod_p(u, U, p):
    r"""
    Return coordinates of a unit ``u`` with respect to a basis of the
    ``p``-cotorsion `U/U^p` of the unit group ``U``.

    INPUT:

    - ``u`` (algebraic unit) -- a unit in a number field ``K``.

    - ``U`` (unit group) -- the unit group of ``K``.

    - ``p`` (prime) -- a prime number.

    OUTPUT:

    The coordinates of the unit `u` in the `p`-cotorsion group `U/U^p`.

    ALGORITHM:

    Take the coordinate vector of `u` with respect to the generators
    of the unit group, drop the coordinate of the roots of unity
    factor if it is prime to `p`, and reduce the vector mod `p`.

    EXAMPLES::

        sage: from sage.rings.number_field.selmer_group import coords_in_U_mod_p
        sage: K.<a> = NumberField(x^4 - 5*x^2 + 1)
        sage: U = K.unit_group()
        sage: U
        Unit group with structure C2 x Z x Z x Z of Number Field in a with defining polynomial x^4 - 5*x^2 + 1
        sage: u0, u1, u2, u3 = U.gens_values()
        sage: u = u1*u2^2*u3^3
        sage: coords_in_U_mod_p(u,U,2)
        [0, 1, 0, 1]
        sage: coords_in_U_mod_p(u,U,3)
        [1, 2, 0]
        sage: u*=u0
        sage: coords_in_U_mod_p(u,U,2)
        [1, 1, 0, 1]
        sage: coords_in_U_mod_p(u,U,3)
        [1, 2, 0]

    """
    coords = U.log(u)
    start = 1 - int(p.divides(U.zeta_order())) # 0 or 1
    return  [c%p for c in coords[start:]]

def basis_for_p_cokernel(S, C, p):
    r"""
    Return a basis for the group of ideals supported on ``S`` (mod
    ``p``'th-powers) whose class in the class group ``C`` is a ``p``'th power,
    together with a function which takes the ``S``-exponents of such an
    ideal and returns its coordinates on this basis.

    INPUT:

    - ``S`` (list) -- a list of prime ideals in a number field ``K``.

    - ``C`` (class group) -- the ideal class group of ``K``.

    - ``p`` (prime) -- a prime number.

    OUTPUT:

    (tuple) (``b``, ``f``) where

    - ``b`` is a list of ideals which is a basis for the group of
      ideals supported on ``S`` (modulo ``p``'th powers) whose ideal
      class is a ``p``'th power;

    - ``f`` is a function which takes such an ideal and returns its
      coordinates with respect to this basis.

    EXAMPLES::

        sage: from sage.rings.number_field.selmer_group import basis_for_p_cokernel
        sage: K.<a> = NumberField(x^2 - x + 58)
        sage: S = K.ideal(30).support(); S
        [Fractional ideal (2, a),
        Fractional ideal (2, a + 1),
        Fractional ideal (3, a + 1),
        Fractional ideal (5, a + 1),
        Fractional ideal (5, a + 3)]
        sage: C = K.class_group()
        sage: C.gens_orders()
        (6, 2)
        sage: [C(P).exponents() for P in S]
        [(5, 0), (1, 0), (3, 1), (1, 1), (5, 1)]
        sage: b, f = basis_for_p_cokernel(S, C, 2); b
        [Fractional ideal (2), Fractional ideal (15, a + 13), Fractional ideal (5)]
        sage: b, f = basis_for_p_cokernel(S, C, 3); b
        [Fractional ideal (50, a + 18),
        Fractional ideal (10, a + 3),
        Fractional ideal (3, a + 1),
        Fractional ideal (5)]
        sage: b, f = basis_for_p_cokernel(S, C, 5); b
        [Fractional ideal (2, a),
        Fractional ideal (2, a + 1),
        Fractional ideal (3, a + 1),
        Fractional ideal (5, a + 1),
        Fractional ideal (5, a + 3)]

    """
    from sage.matrix.constructor import Matrix
    M = Matrix(GF(p), [_coords_in_C_mod_p(P, C, p) for P in S])
    k = M.left_kernel()
    bas = [prod([P ** bj.lift() for P, bj in zip(S, b.list())],
                C.number_field().ideal(1)) for b in k.basis()]
    return bas, k.coordinate_vector

# The main function


def pSelmerGroup(K, S, p, proof=None, debug=False):
    r"""
    Return the ``p``-Selmer group `K(S,p)` of the number field ``K``
    with respect to the prime ideals in ``S``

    INPUT:

    - ``K`` (number field) -- a number field, or `\QQ`.

    - ``S`` (list) -- a list of prime ideals in ``K``, or prime
      numbers when ``K`` is `\QQ`.

    - ``p`` (prime) -- a prime number.

    - ``proof`` - if True then compute the class group provably
      correctly. Default is True. Call :meth:`proof.number_field` to
      change this default globally.

    - ``debug`` (boolean, default ``False``) -- debug flag.

    OUTPUT:

    (tuple) ``KSp``, ``KSp_gens``, ``from_KSp``, ``to_KSp`` where

    - ``KSp`` is an abstract vector space over `GF(p)` isomorphic to `K(S,p)`;

    - ``KSp_gens`` is a list of elements of `K^*` generating `K(S,p)`;

    - ``from_KSp`` is a function from ``KSp`` to `K^*` implementing
      the isomorphism from the abstract `K(S,p)` to `K(S,p)` as a
      subgroup of `K^*/(K^*)^p`;

    - ``to_KSP`` is a partial function from `K^*` to ``KSp``, defined
      on elements `a` whose image in `K^*/(K^*)^p` lies in `K(S,p)`,
      mapping them via the inverse isomorphism to the abstract vector
      space ``KSp``.

    ALGORITHM:

    The list of generators of `K(S,p)` is the concatenation of three
    sublists, called ``alphalist``, ``betalist`` and ``ulist`` in the
    code.  Only ``alphalist`` depends on the primes in `S`.

    - ``ulist`` is a basis for `U/U^p` where `U` is the unit group.
      This is the list of fundamental units, including the generator
      of the group of roots of unity if its order is divisible by `p`.
      These have valuation `0` at all primes.

    - ``betalist`` is a list of the generators of the `p`'th powers of
      ideals which generate the `p`-torsion in the class group (so is
      empty if the class number is prime to `p`).  These have
      valuation divisible by `p` at all primes.

    - ``alphalist`` is a list of generators for each ideal `A` in a
      basis of those ideals supported on `S` (modulo `p`'th powers of
      ideals) which are `p`'th powers in the class group.  We find `B`
      such that `A/B^p` is principal and take a generator of it, for
      each `A` in a generating set.  As a special case, if all the
      ideals in `S` are principal then ``alphalist`` is a list of
      their generators.

    The map from the abstract space to `K^*` is easy: we just take the
    product of the generators to powers given by the coefficient
    vector.  No attempt is made to reduce the resulting product modulo
    `p`'th powers.

    The reverse map is more complicated.  Given `a\in K^*`:

    - write the principal ideal `(a)` in the form `AB^p` with `A`
      supported by `S` and `p`'th power free.  If this fails, then `a`
      does not represent an element of `K(S,p)` and an error is
      raised.

    - set `I_S` to be the group of ideals spanned by `S` mod `p`'th
      powers, and `I_{S,p}` the subgroup of `I_S` which maps to `0` in
      `C/C^p`.

    - Convert `A` to an element of `I_{S,p}`, hence find the
      coordinates of `a` with respect to the generators in
      ``alphalist``.

    - after dividing out by `A`, now `(a)=B^p` (with a different `a`
      and `B`).  Write the ideal class `[B]`, whose `p`'th power is
      trivial, in terms of the generators of `C[p]`; then `B=(b)B_1`,
      where the coefficients of `B_1` with respect to generators of
      `C[p]` give the coordinates of the result with respect to the
      generators in ``betalist``.

    - after dividing out by `B`, and by `b^p`, we now have `(a)=(1)`,
      so `a` is a unit, which can be expressed in terms of the unit
      generators.

    EXAMPLES:

    Over `\QQ` the the unit contribution is trivial unless `p=2` and
    the class group is trivial::

        sage: from sage.rings.number_field.selmer_group import pSelmerGroup
        sage: QS2, gens, fromQS2, toQS2 = pSelmerGroup(QQ, [2,3], 2)
        sage: QS2
        Vector space of dimension 3 over Finite Field of size 2
        sage: gens
        [2, 3, -1]
        sage: a = fromQS2([1,1,1]); a.factor()
        -1 * 2 * 3
        sage: toQS2(-6)
        (1, 1, 1)

        sage: QS3, gens, fromQS3, toQS3 = pSelmerGroup(QQ, [2,13], 3)
        sage: QS3
        Vector space of dimension 2 over Finite Field of size 3
        sage: gens
        [2, 13]
        sage: a = fromQS3([5,4]); a.factor()
        2^5 * 13^4
        sage: toQS3(a)
        (2, 1)
        sage: toQS3(a) == QS3([5,4])
        True

    A real quadratic field with class number 2, where the fundamental
    unit is a generator, and the class group provides another
    generator when `p=2`::

        sage: K.<a> = QuadraticField(-5)
        sage: K.class_number()
        2
        sage: P2 = K.ideal(2, -a+1)
        sage: P3 = K.ideal(3, a+1)
        sage: P5 = K.ideal(a)
        sage: KS2, gens, fromKS2, toKS2 = pSelmerGroup(K, [P2, P3, P5], 2)
        sage: KS2
        Vector space of dimension 4 over Finite Field of size 2
        sage: gens
        [a + 1, a, 2, -1]

    Each generator must have even valuation at primes not in `S`::

        sage: [K.ideal(g).factor() for g in gens]
        [(Fractional ideal (2, a + 1)) * (Fractional ideal (3, a + 1)),
        Fractional ideal (-a),
        (Fractional ideal (2, a + 1))^2,
        1]

        sage: toKS2(10)
        (0, 0, 1, 1)
        sage: fromKS2([0,0,1,1])
        -2
        sage: K(10/(-2)).is_square()
        True

        sage: KS3, gens, fromKS3, toKS3 = pSelmerGroup(K, [P2, P3, P5], 3)
        sage: KS3
        Vector space of dimension 3 over Finite Field of size 3
        sage: gens
        [1/2, 1/4*a + 1/4, a]

    The ``to`` and ``from`` maps are inverses of each other::

        sage: K.<a> = QuadraticField(-5)
        sage: S = K.ideal(30).support()
        sage: KS2, gens, fromKS2, toKS2 = pSelmerGroup(K, S, 2)
        sage: KS2
        Vector space of dimension 5 over Finite Field of size 2
        sage: assert all(toKS2(fromKS2(v))==v for v in KS2)
        sage: KS3, gens, fromKS3, toKS3 = pSelmerGroup(K, S, 3)
        sage: KS3
        Vector space of dimension 4 over Finite Field of size 3
        sage: assert all(toKS3(fromKS3(v))==v for v in KS3)
    """
    from sage.rings.number_field.number_field import proof_flag
    from sage.modules.free_module import VectorSpace
    from sage.sets.set import Set

    proof = proof_flag(proof)

    # Input check: p and all P in S must be prime.  Remove any repeats in S.

    S = list(Set(S))
    if not all(P.is_prime() for P in S):
        raise ValueError("elements of S must all be prime")
    if not p.is_prime():
        raise ValueError("p must be prime")

    F = GF(p)

    # Step 1. The unit contribution: all fundamental units, and also the
    # generating root of unity if its order is a multiple of p; we just
    # take generators of U/U^p.  These have valuation 0 everywhere.

    hK = 1 if K == QQ else K.class_number(proof=proof)
    C = K.class_group() if K == QQ else K.class_group(proof=proof)

    hKp = (hK%p == 0) # flag whether the class number is divisible by p

    if K == QQ:
        if p == 2:
            ulist = [QQ(-1)]
        else:
            ulist = []
    else:
        U = K.unit_group(proof=proof)
        ulist = U.gens_values()
        if U.zeta_order() % p:
            ulist = ulist[1:]

    if debug:
        print("{} generators in ulist = {}".format(len(ulist),ulist))

    # Step 2. The class group contribution: generators of the p'th
    # powers of ideals generating the p-torsion in the class group.
    # These have valuation divisible by p everywhere.

    if hKp:
        betalist = [_ideal_generator(c ** n)
                    for c, n in zip(C.gens_ideals(), C.gens_orders())
                    if n % p == 0]
    else:
        betalist = []

    if debug:
        print("{} generators in betalist = {}".format(len(betalist),betalist))

    # Step 3. The part depending on S: one generator for each ideal A
    # in a basis of those ideals supported on S (modulo p'th powers of
    # ideals) which is a p'th power in the class group.  We find B
    # such that A/B^p is principal and take a generator of that, for
    # each A in a generating set.

    # As a special case, when the class number is 1 we just take
    # generators of the primes in S.

    if hK > 1:
        T, f = basis_for_p_cokernel(S, C, p)
        alphalist = [_ideal_generator(I / _root_ideal(I, C, p) ** p) for I in T]
    else:
        f = lambda x:x
        alphalist = [_ideal_generator(P) for P in S]

    if debug:
        print("{} generators in alphalist = {}".format(len(alphalist), alphalist))

    # Now we have the generators of K(S,p), and define K(S,p) as an
    # abstract vector space:

    KSp_gens = alphalist + betalist + ulist
    KSp = VectorSpace(GF(p), len(KSp_gens))

    if debug:
        print("Generators of K(S,p) = {} (dimension {})".format(KSp_gens, len(KSp_gens)))

    # Now we define maps in each direction between the abstract space and K^*.

    # Define the easy map from KSp into K^*:

    def from_KSp(v):
        return prod([g ** vi for g, vi in zip(KSp_gens, v)], K(1))

    # Define the hard map from (a subgroup of) K^* to KSp:

    def to_KSp(a):
        # Check that a is in K(S,p):

        if not a:
            raise ValueError("argument {} should be nonzero".format(a))
        try:
            a = K(a)
        except ValueError:
            raise ValueError("argument {} should be in {}".format(a, K))

        if not all(P in S or a.valuation(P) % p == 0 for P in a.support()):
            raise ValueError("argument {} should have valuations divisible by {} at all primes in {}".format(a, p, S))

        # 1. (a) is a p'th power mod ideals in S, say (a)=AB^p, where
        # A is supported on S and is a linear combination of the
        # ideals T above.  Find the exponents of the P_i in S in A:

        S_vals = [F(a.valuation(P)) for P in S]
        avec = list(f(S_vals)) # coordinates of A w.r.t ideals in T (mod p'th powers)
        a1 = prod((alpha ** e for alpha, e in zip(alphalist,avec)), K(1))
        a /= a1
        if debug:
            print("alpha component is {} with coords {}".format(a1,avec))
            if K == QQ:
                print("continuing with quotient {} whose ideal should be a {}'th power: {}".format(a,p,a.factor()))
            else:
                print("continuing with quotient {} whose ideal should be a {}'th power: {}".format(a,p,K.ideal(a).factor()))

        # 2. Now (a) is a p'th power, say (a)=B^p.
        # Find B and the exponents of [B] w.r.t. basis of C[p]:

        supp = a.support()
        vals = [a.valuation(P) for P in supp]
        if debug:
            assert all(v % p == 0 for v in vals)
        one = K(1)    if K == QQ else K.ideal(1)
        aa  = a.abs() if K == QQ else K.ideal(a)
        B = prod((P ** (v // p) for P, v in zip(supp,vals)), one)
        if debug:
            assert B ** p == aa
            print("B={}".format(B))
            print("a={}".format(a))

        if hKp:
            bvec = _coords_in_C_p(B, C, p)
            a2 = prod((beta ** e for beta, e in zip(betalist, bvec)), K(1))
            a /= a2
            supp = a.support()
            vals = [a.valuation(P) for P in supp]
            if debug:
                assert all(v % p == 0 for v in vals)
            B = prod((P ** (v // p) for P, v in zip(supp, vals)), one)
            if debug:
                assert B ** p == aa
        else:
            bvec = []
            a2 = 1

        if debug:
            print("beta component is {} with coords {}".format(a2,bvec))
            print("continuing with quotient {} which should be a p'th power times a unit".format(a))

        # 3. Now (a) = (c)^p for some c, so a/c^p is a unit

        if K != QQ:
            assert B.is_principal()

        if debug:
            print("B={}".format(B))
        a3 = B if K==QQ else _ideal_generator(B)
        if debug:
            print("a3={}".format(a3))
        a /= a3 ** p
        if debug:
            print("dividing by {}th power of {}".format(p,a3))
            print("continuing with quotient {} which should be a unit".format(a))

        #4. Now a is a unit

        # NB not a.is_unit which is true for all a in K^*.  One could
        # also test K.ring_of_integers()(a).is_unit().

        if debug:
            if K == QQ:
                assert a.abs()==1
            else:
                assert K.ideal(a).is_one()

        if K == QQ:
            if p == 2:
                cvec = [1] if a == -1 else [0]
            else:
                cvec = []
        else:
            cvec = coords_in_U_mod_p(a,U,p)

        if debug:
            print("gamma component has coords {}".format(cvec))

        return KSp(avec + bvec + cvec)

    return KSp, KSp_gens, from_KSp, to_KSp
