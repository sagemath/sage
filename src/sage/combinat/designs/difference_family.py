r"""
Difference families

This module gathers everything related to difference families. One can build a
difference family (or check that it can be built) with :func:`difference_family`::

    sage: G,F = designs.difference_family(13,4,1)                                       # needs sage.libs.pari sage.modules

It defines the following functions:

{INDEX_OF_FUNCTIONS}

REFERENCES:

.. [BJL99-1] \T. Beth, D. Jungnickel, H. Lenz "Design theory Vol. I."
   Second edition. Encyclopedia of Mathematics and its Applications, 69. Cambridge
   University Press, (1999).

.. [BLJ99-2] \T. Beth, D. Jungnickel, H. Lenz "Design theory Vol. II."
   Second edition. Encyclopedia of Mathematics and its Applications, 78. Cambridge
   University Press, (1999).

.. [Bo39] \R. C. Bose, "On the construction of balanced incomplete block
   designs", Ann. Eugenics, 9 (1939), 353--399.

.. [Bu95] \M. Buratti "On simple radical difference families", J.
   Combinatorial Designs, 3 (1995) 161--168.

.. [Tu1965] \R. J. Turyn "Character sum and difference sets"
   Pacific J. Math. 15 (1965) 319--346.

.. [Tu1984] \R. J. Turyn "A special class of Williamson matrices and
   difference sets" J. Combinatorial Theory (A) 36 (1984) 111--115.

.. [Wi72] \R. M. Wilson "Cyclotomy and difference families in elementary Abelian
   groups", J. Number Theory, 4 (1972) 17--47.

Functions
---------
"""
# ****************************************************************************
#       Copyright (C) 2014 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import is_prime_power
from sage.misc.cachefunc import cached_function

from sage.categories.sets_cat import EmptySetError
import sage.arith.all as arith
from sage.misc.unknown import Unknown
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


def group_law(G):
    r"""
    Return a triple ``(identity, operation, inverse)`` that define the
    operations on the group ``G``.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import group_law
        sage: group_law(Zmod(3))
        (0, <built-in function add>, <built-in function neg>)
        sage: group_law(SymmetricGroup(5))                                              # needs sage.groups
        ((), <built-in function mul>, <built-in function inv>)
        sage: group_law(VectorSpace(QQ, 3))                                             # needs sage.modules
        ((0, 0, 0), <built-in function add>, <built-in function neg>)
    """
    import operator
    from sage.categories.groups import Groups
    from sage.categories.additive_groups import AdditiveGroups

    if G in Groups():            # multiplicative groups
        return (G.one(), operator.mul, operator.inv)
    elif G in AdditiveGroups():  # additive groups
        return (G.zero(), operator.add, operator.neg)
    else:
        raise ValueError("%s does not seem to be a group" % G)


def block_stabilizer(G, B):
    r"""
    Compute the left stabilizer of the block ``B`` under the action of ``G``.

    This function return the list of all `x\in G` such that `x\cdot B=B` (as a
    set).

    INPUT:

    - ``G`` -- a group (additive or multiplicative)
    - ``B`` -- a subset of ``G``

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import block_stabilizer

        sage: Z8 = Zmod(8)
        sage: block_stabilizer(Z8, [Z8(0),Z8(2),Z8(4),Z8(6)])
        [0, 2, 4, 6]
        sage: block_stabilizer(Z8, [Z8(0),Z8(2)])
        [0]

        sage: C = cartesian_product([Zmod(4),Zmod(3)])
        sage: block_stabilizer(C, [C((0,0)),C((2,0)),C((0,1)),C((2,1))])
        [(0, 0), (2, 0)]

        sage: b = list(map(Zmod(45),[1, 3, 7, 10, 22, 25, 30, 35, 37, 38, 44]))
        sage: block_stabilizer(Zmod(45),b)
        [0]
    """
    if not B:
        return list(G)
    identity, op, inv = group_law(G)
    b0 = inv(B[0])
    S = []
    for b in B:
        # fun: if we replace +(-b) with -b it completely fails!!
        bb0 = op(b, b0)  # bb0 = b-B[0]
        if all(op(bb0, c) in B for c in B):
            S.append(bb0)
    return S


def is_difference_family(G, D, v=None, k=None, l=None, verbose=False):
    r"""
    Check whether ``D`` forms a difference family in the group ``G``.

    INPUT:

    - ``G`` -- group of cardinality ``v``
    - ``D`` -- set of ``k``-subsets of ``G``
    - ``v``, ``k``, ``l`` -- (optional) parameters of the difference family
    - ``verbose`` -- boolean (default: ``False``); whether to print additional
      information

    .. SEEALSO::

        :func:`difference_family`

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import is_difference_family
        sage: G = Zmod(21)
        sage: D = [[0,1,4,14,16]]
        sage: is_difference_family(G, D, 21, 5)
        True

        sage: G = Zmod(41)
        sage: D = [[0,1,4,11,29],[0,2,8,17,21]]
        sage: is_difference_family(G, D, verbose=True)
        Too few:
          5 is obtained 0 times in blocks []
          14 is obtained 0 times in blocks []
          27 is obtained 0 times in blocks []
          36 is obtained 0 times in blocks []
        Too much:
          4 is obtained 2 times in blocks [0, 1]
          13 is obtained 2 times in blocks [0, 1]
          28 is obtained 2 times in blocks [0, 1]
          37 is obtained 2 times in blocks [0, 1]
        False
        sage: D = [[0,1,4,11,29],[0,2,8,17,22]]
        sage: is_difference_family(G, D)
        True

        sage: G = Zmod(61)
        sage: D = [[0,1,3,13,34],[0,4,9,23,45],[0,6,17,24,32]]
        sage: is_difference_family(G, D)
        True

        sage: # needs sage.modules
        sage: G = AdditiveAbelianGroup([3]*4)
        sage: a,b,c,d = G.gens()
        sage: D = [[d, -a+d, -c+d, a-b-d, b+c+d],
        ....:      [c, a+b-d, -b+c, a-b+d, a+b+c],
        ....:      [-a-b+c+d, a-b-c-d, -a+c-d, b-c+d, a+b],
        ....:      [-b-d, a+b+d, a-b+c-d, a-b+c, -b+c+d]]
        sage: is_difference_family(G, D)
        True

    The following example has a third block with a non-trivial stabilizer::

        sage: G = Zmod(15)
        sage: D = [[0,1,4],[0,2,9],[0,5,10]]
        sage: is_difference_family(G,D,verbose=True)
        It is a (15,3,1)-difference family
        True

    The function also supports multiplicative groups (non necessarily Abelian)::

        sage: # needs sage.groups
        sage: G = DihedralGroup(8)
        sage: x,y = G.gens()
        sage: i = G.one()
        sage: D1 = [[i,x,x^4], [i,x^2, y*x], [i,x^5,y], [i,x^6,y*x^2], [i,x^7,y*x^5]]
        sage: is_difference_family(G, D1, 16, 3, 2)
        True
        sage: from sage.combinat.designs.bibd import BIBD_from_difference_family
        sage: bibd = BIBD_from_difference_family(G, D1, lambd=2)

    TESTS::

        sage: # needs sage.rings.finite_rings
        sage: K = GF(3^2,'z')
        sage: z = K.gen()
        sage: D = [[1,z+1,2]]
        sage: _ = is_difference_family(K, D, verbose=True)
        the number of differences (=6) must be a multiple of v-1=8
        sage: _
        False
    """
    identity, mul, inv = group_law(G)

    Glist = list(G)

    D = [[G(_) for _ in d] for d in D]

    # Check v (and define it if needed)
    if v is None:
        v = len(Glist)
    else:
        if len(Glist) != v:
            if verbose:
                print("G must have cardinality v (=%d)" % int(v))
            return False

    # Check k (and define it if needed)
    if k is None:
        k = len(D[0])
    else:
        k = int(k)

    for d in D:
        if len(d) != k:
            if verbose:
                print("the block {} does not have length {}".format(d, k))
            return False

    # Check l (and define it if needed)
    #
    # - nb_diff: the number of pairs (with multiplicity) covered by the BIBD
    #            generated by the DF.
    #
    # - stab: the stabilizer of each set.
    nb_diff = 0
    stab = []
    for d in D:
        s = block_stabilizer(G,d)
        stab.append(s)
        nb_diff += k*(k-1) // len(s)
    if l is None:
        if nb_diff % (v-1) != 0:
            if verbose:
                print("the number of differences (={}) must be a multiple of v-1={}".format(nb_diff, v-1))
            return False
        l = nb_diff // (v-1)
    else:
        if nb_diff != l*(v-1):
            if verbose:
                print("the number of differences (={}) is not equal to l*(v-1) = {}".format(nb_diff, l*(v-1)))
            return False

    # Check that every x \in G-{0},occurs exactly l times as a difference
    counter = {g: 0 for g in Glist}
    where = {g: set() for g in Glist}
    del counter[identity]

    for i,d in enumerate(D):
        tmp_counter = {}
        for b in d:
            for c in d:
                if b == c:
                    continue
                gg = mul(b, inv(c))  # = b-c or bc^{-1}
                if gg not in tmp_counter:
                    tmp_counter[gg] = 0
                where[gg].add(i)
                tmp_counter[gg] += 1

        if sum(tmp_counter.values()) != k * (k - 1):
            if verbose:
                print("repeated element in the {}-th block {}".format(i, d))
            return False

        # Normalized number of occurrences added to counter
        stabi = len(stab[i])
        for gg in tmp_counter:
            counter[gg] += tmp_counter[gg]//stabi

    # Check the counter and report any error
    too_few = []
    too_much = []
    for g in Glist:
        if g == identity:
            continue
        if counter[g] < l:
            if verbose:
                too_few.append(g)
            else:
                return False
        if counter[g] > l:
            if verbose:
                too_much.append(g)
            else:
                return False

    if too_few:
        print("Too few:")
        for g in too_few:
            print("  {} is obtained {} times in blocks {}".format(
                        g, counter[g], sorted(where[g])))
    if too_much:
        print("Too much:")
        for g in too_much:
            print("  {} is obtained {} times in blocks {}".format(
                        g, counter[g], sorted(where[g])))
    if too_few or too_much:
        return False

    if verbose:
        print("It is a ({},{},{})-difference family".format(v, k, l))
    return True


def singer_difference_set(q, d):
    r"""
    Return a difference set associated to the set of hyperplanes in a projective
    space of dimension `d` over `GF(q)`.

    A Singer difference set has parameters:

    .. MATH::

        v = \frac{q^{d+1}-1}{q-1}, \quad
        k = \frac{q^d-1}{q-1}, \quad
        \lambda = \frac{q^{d-1}-1}{q-1}.

    The idea of the construction is as follows. One consider the finite field
    `GF(q^{d+1})` as a vector space of dimension `d+1` over `GF(q)`. The set of
    `GF(q)`-lines in `GF(q^{d+1})` is a projective plane and its set of
    hyperplanes form a balanced incomplete block design.

    Now, considering a multiplicative generator `z` of `GF(q^{d+1})`, we get a
    transitive action of a cyclic group on our projective plane from which it is
    possible to build a difference set.

    The construction is given in details in [Stinson2004]_, section 3.3.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import singer_difference_set, is_difference_family
        sage: G,D = singer_difference_set(3,2)                                          # needs sage.rings.finite_rings
        sage: is_difference_family(G, D, verbose=True)                                  # needs sage.rings.finite_rings
        It is a (13,4,1)-difference family
        True

        sage: G,D = singer_difference_set(4,2)                                          # needs sage.rings.finite_rings
        sage: is_difference_family(G, D, verbose=True)                                  # needs sage.rings.finite_rings
        It is a (21,5,1)-difference family
        True

        sage: G,D = singer_difference_set(3,3)                                          # needs sage.rings.finite_rings
        sage: is_difference_family(G, D, verbose=True)                                  # needs sage.rings.finite_rings
        It is a (40,13,4)-difference family
        True

        sage: G,D = singer_difference_set(9,3)                                          # needs sage.rings.finite_rings
        sage: is_difference_family(G, D, verbose=True)                                  # needs sage.rings.finite_rings
        It is a (820,91,10)-difference family
        True
    """
    q = Integer(q)
    assert q.is_prime_power()
    assert d >= 2

    from sage.rings.finite_rings.finite_field_constructor import GF
    from sage.rings.finite_rings.conway_polynomials import conway_polynomial
    from sage.rings.finite_rings.integer_mod_ring import Zmod

    # build a polynomial c over GF(q) such that GF(q)[x] / (c(x)) is a
    # GF(q**(d+1)) and such that x is a multiplicative generator.
    p,e = q.factor()[0]
    c = conway_polynomial(p, e*(d+1))
    if e != 1:
        # i.e. q is not a prime, so we factorize c over GF(q) and pick
        # one of its factor
        K = GF(q,'z')
        c = c.change_ring(K).factor()[0][0]
    else:
        K = GF(q)
    z = c.parent().gen()

    # Now we consider the GF(q)-subspace V spanned by (1,z,z^2,...,z^(d-1)) inside
    # GF(q^(d+1)). The multiplication by z is an automorphism of the
    # GF(q)-projective space built from GF(q^(d+1)). The difference family is
    # obtained by taking the integers i such that z^i belong to V.
    powers = [0]
    i = 1
    x = z
    k = (q**d-1)//(q-1)
    while len(powers) < k:
        if x.degree() <= (d-1):
            powers.append(i)
        x = (x*z).mod(c)
        i += 1

    return Zmod((q**(d+1)-1)//(q-1)), [powers]


def df_q_6_1(K, existence=False, check=True):
    r"""
    Return a `(q,6,1)`-difference family over the finite field `K`.

    The construction uses Theorem 11 of [Wi72]_.

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: from sage.combinat.designs.difference_family import is_difference_family, df_q_6_1
        sage: prime_powers = [v for v in range(31,500,30) if is_prime_power(v)]
        sage: parameters = [v for v in prime_powers
        ....:               if df_q_6_1(GF(v,'a'), existence=True) is True]
        sage: parameters
        [31, 151, 181, 211, 241, 271, 331, 361, 421]
        sage: for v in parameters:
        ....:     K = GF(v, 'a')
        ....:     df = df_q_6_1(K, check=True)
        ....:     assert is_difference_family(K, df, v, 6, 1)

    .. TODO::

        Do improvements due to Zhen and Wu 1999.
    """
    v = K.cardinality()
    x = K.multiplicative_generator()
    one = K.one()
    if v % 30 != 1:
        if existence:
            return False
        raise EmptySetError("k(k-1)=30 should divide (v-1)")
    t = (v-1) // 30  # number of blocks

    r = x**((v-1)//3)  # primitive cube root of unity
    r2 = r*r           # the other primitive cube root

    # we now compute the cosets of x**i
    xx = x**5
    to_coset = {x**i * xx**j: i for i in range(5) for j in range((v-1)/5)}

    for c in to_coset:  # the loop runs through all nonzero elements of K
        if c == one or c == r or c == r2:
            continue
        if len(set(to_coset[elt] for elt in (r-one, c*(r-one), c-one, c-r, c-r**2))) == 5:
            if existence:
                return True
            B = [one,r,r2,c,c*r,c*r2]
            D = [[xx**i * b for b in B] for i in range(t)]
            break
    else:
        if existence:
            return Unknown
        raise NotImplementedError("Wilson construction failed for v={}".format(v))

    if check and not is_difference_family(K, D, v, 6, 1):
        raise RuntimeError("Wilson 1972 construction failed! Please e-mail sage-devel@googlegroups.com")

    return D


def radical_difference_set(K, k, l=1, existence=False, check=True):
    r"""
    Return a difference set made of a cyclotomic coset in the finite field
    ``K`` and with parameters ``k`` and ``l``.

    Most of these difference sets appear in chapter VI.18.48 of the Handbook of
    combinatorial designs.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import radical_difference_set

        sage: D = radical_difference_set(GF(7), 3, 1); D                                # needs sage.rings.finite_rings
        [[1, 2, 4]]
        sage: sorted(x-y for x in D[0] for y in D[0] if x != y)                         # needs sage.rings.finite_rings
        [1, 2, 3, 4, 5, 6]

        sage: D = radical_difference_set(GF(16,'a'), 6, 2)                              # needs sage.rings.finite_rings
        sage: sorted(x-y for x in D[0] for y in D[0] if x != y)                         # needs sage.rings.finite_rings
        [1,
         1,
         a,
         a,
         a + 1,
         a + 1,
         a^2,
         a^2,
         ...
         a^3 + a^2 + a + 1,
         a^3 + a^2 + a + 1]

        sage: for k in range(2,50):                                                     # needs sage.rings.finite_rings
        ....:     for l in reversed(divisors(k*(k-1))):
        ....:         v = k*(k-1)//l + 1
        ....:         if is_prime_power(v) and radical_difference_set(GF(v,'a'),k,l,existence=True) is True:
        ....:             _ = radical_difference_set(GF(v,'a'),k,l)
        ....:             print("{:3} {:3} {:3}".format(v,k,l))
          3   2   1
          4   3   2
          7   3   1
          5   4   3
          7   4   2
         13   4   1
         11   5   2
          7   6   5
         11   6   3
         16   6   2
          8   7   6
          9   8   7
         19   9   4
         37   9   2
         73   9   1
         11  10   9
         19  10   5
         23  11   5
         13  12  11
         23  12   6
         27  13   6
         27  14   7
         16  15  14
         31  15   7
        ...
         41  40  39
         79  40  20
         83  41  20
         43  42  41
         83  42  21
         47  46  45
         49  48  47
        197  49  12
    """
    v = K.cardinality()

    if l*(v-1) != k*(k-1):
        if existence:
            return False
        raise EmptySetError("l*(v-1) is not equal to k*(k-1)")

    # trivial case
    if (v-1) == k:
        if existence:
            return True
        add_zero = False

    # q = 3 mod 4
    elif v % 4 == 3 and k == (v-1)//2:
        if existence:
            return True
        add_zero = False

    # q = 3 mod 4
    elif v % 4 == 3 and k == (v+1)//2:
        if existence:
            return True
        add_zero = True

    # q = 4t^2 + 1, t odd
    elif v % 8 == 5 and k == (v-1)//4 and arith.is_square((v-1)//4):
        if existence:
            return True
        add_zero = False

    # q = 4t^2 + 9, t odd
    elif v % 8 == 5 and k == (v+3)//4 and arith.is_square((v-9)//4):
        if existence:
            return True
        add_zero = True

    # exceptional case 1
    elif (v,k,l) == (16,6,2):
        if existence:
            return True
        add_zero = True

    # exceptional case 2
    elif (v,k,l) == (73,9,1):
        if existence:
            return True
        add_zero = False

    # are there more ??
    else:
        x = K.multiplicative_generator()
        D = K.cyclotomic_cosets(x**((v-1)//k), [K.one()])
        if is_difference_family(K, D, v, k, l):
            print("**  You found a new example of radical difference set **\n"
                  "**  for the parameters (v,k,l)=({},{},{}).            **\n"
                  "**  Please contact sage-devel@googlegroups.com        **\n".format(v, k, l))
            if existence:
                return True
            add_zero = False

        else:
            D = K.cyclotomic_cosets(x**((v-1)//(k-1)), [K.one()])
            D[0].insert(0,K.zero())
            if is_difference_family(K, D, v, k, l):
                print("**  You found a new example of radical difference set **\n"
                      "**  for the parameters (v,k,l)=({},{},{}).            **\n"
                      "**  Please contact sage-devel@googlegroups.com        **\n".format(v, k, l))
                if existence:
                    return True
                add_zero = True

            elif existence:
                return False
            else:
                raise EmptySetError("no radical difference set exist "
                        "for the parameters (v,k,l) = ({},{},{}".format(v,k,l))

    x = K.multiplicative_generator()
    if add_zero:
        r = x**((v-1)//(k-1))
        D = K.cyclotomic_cosets(r, [K.one()])
        D[0].insert(0, K.zero())
    else:
        r = x**((v-1)//k)
        D = K.cyclotomic_cosets(r, [K.one()])

    if check and not is_difference_family(K, D, v, k, l):
        raise RuntimeError("Sage tried to build a radical difference set with "
                "parameters ({},{},{}) but it seems that it failed! Please "
                "e-mail sage-devel@googlegroups.com".format(v,k,l))

    return D


def one_cyclic_tiling(A, n):
    r"""
    Given a subset ``A`` of the cyclic additive group `G = Z / nZ` return
    another subset `B` so that `A + B = G` and `|A| |B| = n` (i.e. any element
    of `G` is uniquely expressed as a sum `a+b` with `a` in `A` and `b` in `B`).

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import one_cyclic_tiling
        sage: tile = [0,2,4]
        sage: m = one_cyclic_tiling(tile,6); m
        [0, 3]
        sage: sorted((i+j)%6 for i in tile for j in m)
        [0, 1, 2, 3, 4, 5]

        sage: def print_tiling(tile, translat, n):
        ....:     for x in translat:
        ....:         print(''.join('X' if (i-x)%n in tile else '.' for i in range(n)))

        sage: tile = [0, 1, 2, 7]
        sage: m = one_cyclic_tiling(tile, 12)
        sage: print_tiling(tile, m, 12)
        XXX....X....
        ....XXX....X
        ...X....XXX.

        sage: tile = [0, 1, 5]
        sage: m = one_cyclic_tiling(tile, 12)
        sage: print_tiling(tile, m, 12)
        XX...X......
        ...XX...X...
        ......XX...X
        ..X......XX.

        sage: tile = [0, 2]
        sage: m = one_cyclic_tiling(tile, 8)
        sage: print_tiling(tile, m, 8)
        X.X.....
        ....X.X.
        .X.X....
        .....X.X

    ALGORITHM:

    Uses dancing links :mod:`sage.combinat.dlx`
    """
    # we first try a naive approach which correspond to what Wilson used in his
    # 1972 article
    n = int(n)
    d = len(A)
    if len(set(a % d for a in A)) == d:
        return [i*d for i in range(n//d)]

    # next, we consider an exhaustive search
    from sage.combinat.dlx import DLXMatrix

    rows = []
    for i in range(n):
        rows.append([i+1, [(i+a) % n+1 for a in A]])
    M = DLXMatrix(rows)
    for c in M:
        return [i-1 for i in c]


def one_radical_difference_family(K, k):
    r"""
    Search for a radical difference family on ``K`` using dancing links
    algorithm.

    For the definition of radical difference family, see
    :func:`radical_difference_family`. Here, we consider only radical difference
    family with `\lambda = 1`.

    INPUT:

    - ``K`` -- a finite field of cardinality `q`
    - ``k`` -- positive integer so that `k(k-1)` divides `q-1`

    OUTPUT: either a difference family or ``None`` if it does not exist

    ALGORITHM:

    The existence of a radical difference family is equivalent to a one
    dimensional tiling (or packing) problem in a cyclic group. This subsequent
    problem is solved by a call to the function :func:`one_cyclic_tiling`.

        Let `K^*` be the multiplicative group of the finite field `K`. A radical
        family has the form `\mathcal B = \{x_1 B, \ldots, x_k B\}`, where
        `B=\{x:x^{k}=1\}` (for `k` odd) or `B=\{x:x^{k-1}=1\}\cup \{0\}` (for
        `k` even). Equivalently, `K^*` decomposes as:

        .. MATH::

            K^* = \Delta (x_1 B) \cup \cdots \cup \Delta (x_k B)
            = x_1 \Delta B \cup \cdots \cup x_k \Delta B.

        We observe that `C=B\backslash 0` is a subgroup of the (cyclic) group
        `K^*`, that can thus be generated by some element `r`. Furthermore, we
        observe that `\Delta B` is always a union of cosets of `\pm C` (which is
        twice larger than `C`).

        .. MATH::

            \begin{array}{llll}
            (k\text{ odd} ) & \Delta B &= \{r^i-r^j:r^i\neq r^j\} &= \pm C\cdot \{r^i-1: 0 < i \leq m\}\\
            (k\text{ even}) & \Delta B &= \{r^i-r^j:r^i\neq r^j\}\cup C &= \pm C\cdot \{r^i-1: 0 < i < m\}\cup \pm C
            \end{array}

        where

        .. MATH::

            (k\text{ odd})\ m = (k-1)/2 \quad \text{and} \quad (k\text{ even})\ m = k/2.

        Consequently, `\mathcal B = \{x_1 B, \ldots, x_k B\}` is a radical
        difference family if and only if `\{x_1 (\Delta B/(\pm C)), \ldots, x_k
        (\Delta B/(\pm C))\}` is a partition of the cyclic group `K^*/(\pm C)`.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import (
        ....:    one_radical_difference_family,
        ....:    is_difference_family)

        sage: one_radical_difference_family(GF(13),4)                                   # needs sage.rings.finite_rings
        [[0, 1, 3, 9]]

    The parameters that appear in [Bu95]_::

        sage: df = one_radical_difference_family(GF(449), 8); df                        # needs sage.rings.finite_rings
        [[0, 1, 18, 25, 176, 324, 359, 444],
         [0, 9, 88, 162, 222, 225, 237, 404],
         [0, 11, 140, 198, 275, 357, 394, 421],
         [0, 40, 102, 249, 271, 305, 388, 441],
         [0, 49, 80, 93, 161, 204, 327, 433],
         [0, 70, 99, 197, 230, 362, 403, 435],
         [0, 121, 141, 193, 293, 331, 335, 382],
         [0, 191, 285, 295, 321, 371, 390, 392]]
        sage: is_difference_family(GF(449), df, 449, 8, 1)                              # needs sage.rings.finite_rings
        True
    """
    q = K.cardinality()
    x = K.multiplicative_generator()

    e = k*(k-1)
    if q % e != 1:
        raise ValueError("q%e is not 1")

    # We define A by (see the function's documentation):
    # ΔB = C.A
    if k % 2 == 1:
        m = (k-1) // 2
        r = x ** ((q-1) // k)     # k-th root of unity
        A = [r**i - 1 for i in range(1,m+1)]
    else:
        m = k // 2
        r = x ** ((q-1) // (k-1))  # (k-1)-th root of unity
        A = [r**i - 1 for i in range(1,m)]
        A.append(K.one())

    # instead of the complicated multiplicative group K^*/(±C) we use the
    # discrete logarithm to convert everything into the additive group Z/cZ
    c = m * (q-1) // e  # cardinal of ±C
    from sage.groups.generic import discrete_log
    logA = [discrete_log(a,x) % c for a in A]

    # if two elements of A are equal modulo c then no tiling is possible
    if len(set(logA)) != m:
        return None

    # brute force
    tiling = one_cyclic_tiling(logA, c)
    if tiling is None:
        return None

    D = K.cyclotomic_cosets(r, [x**i for i in tiling])
    if k % 2 == 0:
        for d in D:
            d.insert(K.zero(),0)
    return D


def radical_difference_family(K, k, l=1, existence=False, check=True):
    r"""
    Return a ``(v,k,l)``-radical difference family.

    Let fix an integer `k` and a prime power `q = t k(k-1) + 1`. Let `K` be a
    field of cardinality `q`. A `(q,k,1)`-difference family is *radical* if
    its base blocks are either: a coset of the `k`-th root of unity for `k` odd
    or a coset of `k-1`-th root of unity and `0` if `k` is even (the number `t`
    is the number of blocks of that difference family).

    The terminology comes from M. Buratti article [Bu95]_ but the first
    constructions go back to R. Wilson [Wi72]_.

    INPUT:

    - ``K`` -- a finite field
    - ``k`` -- positive integer; the size of the blocks
    - ``l`` -- integer (default: `1`); the `\lambda` parameter
    - ``existence`` -- if ``True``, then return either ``True`` if Sage knows
      how to build such design, ``Unknown`` if it does not and ``False`` if it
      knows that the design does not exist
    - ``check`` -- boolean (default: ``True``); if ``True`` then the result of
      the computation is checked before being returned. This should not be
      needed but ensures that the output is correct

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import radical_difference_family

        sage: radical_difference_family(GF(73), 9)                                      # needs sage.rings.finite_rings
        [[1, 2, 4, 8, 16, 32, 37, 55, 64]]

        sage: radical_difference_family(GF(281), 5)                                     # needs sage.rings.finite_rings
        [[1, 86, 90, 153, 232],
         [4, 50, 63, 79, 85],
         [5, 36, 149, 169, 203],
         [7, 40, 68, 219, 228],
         [9, 121, 212, 248, 253],
         [29, 81, 222, 246, 265],
         [31, 137, 167, 247, 261],
         [32, 70, 118, 119, 223],
         [39, 56, 66, 138, 263],
         [43, 45, 116, 141, 217],
         [98, 101, 109, 256, 279],
         [106, 124, 145, 201, 267],
         [111, 123, 155, 181, 273],
         [156, 209, 224, 264, 271]]

        sage: for k in range(5,10):                                                     # needs sage.rings.finite_rings
        ....:     print("k = {}".format(k))
        ....:     list_q = []
        ....:     for q in range(k*(k-1)+1, 2000, k*(k-1)):
        ....:          if is_prime_power(q):
        ....:              K = GF(q,'a')
        ....:              if radical_difference_family(K, k, existence=True) is True:
        ....:                  list_q.append(q)
        ....:                  _ = radical_difference_family(K,k)
        ....:     print(" ".join(str(p) for p in list_q))
        k = 5
        41 61 81 241 281 401 421 601 641 661 701 761 821 881 1181 1201 1301 1321
        1361 1381 1481 1601 1681 1801 1901
        k = 6
        181 211 241 631 691 1531 1831 1861
        k = 7
        337 421 463 883 1723
        k = 8
        449 1009
        k = 9
        73 1153 1873
    """
    v = K.cardinality()
    x = K.multiplicative_generator()
    e = k*(k-1)
    if (l*(v-1)) % e:
        raise ValueError("k (k-1) = {} should be a multiple of l (v-1) ={}".format(
                         k*(k-1), l*(v-1)))
    t = l*(v-1) // e  # number of blocks

    if t == 1:
        return radical_difference_set(K, k, l, existence=existence, check=check)

    elif l == (k-1):
        if existence:
            return True
        else:
            return K.cyclotomic_cosets(x**((v-1)//k))[1:]

    # all the other cases below concern the case l == 1
    elif l != 1:
        if existence:
            return Unknown
        raise NotImplementedError("No radical families implemented for l > 2")

    else:
        D = one_radical_difference_family(K,k)
        if D is None:
            if existence:
                return False
            raise EmptySetError("No such difference family")
        elif existence:
            return True

    if check and not is_difference_family(K, D, v, k, l):
        raise RuntimeError("radical_difference_family produced a wrong "
                           "difference family with parameters v={}, "
                           "k={}, l={}. Please contact "
                           "sage-devel@googlegroups.com".format(v,k,l))

    return D


def twin_prime_powers_difference_set(p, check=True):
    r"""
    Return a difference set on `GF(p) \times GF(p+2)`.

    The difference set is built from the following element of the Cartesian
    product of finite fields `GF(p) \times GF(p+2)`:

    - `(x,0)` with any `x`
    - `(x,y)` with `x` and `y` squares
    - `(x,y)` with `x` and `y` non-squares

    For more information see :wikipedia:`Difference_set`.

    INPUT:

    - ``check`` -- boolean (default: ``True``); if ``True``, then the result of
      the computation is checked before being returned. This should not be
      needed but ensures that the output is correct

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import twin_prime_powers_difference_set
        sage: G, D = twin_prime_powers_difference_set(3)
        sage: G
        The Cartesian product of (Finite Field of size 3, Finite Field of size 5)
        sage: D
        [[(1, 1), (1, 4), (2, 2), (2, 3), (0, 0), (1, 0), (2, 0)]]
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    from sage.categories.cartesian_product import cartesian_product
    from itertools import product
    Fp = FiniteField(p,'x')
    Fq = FiniteField(p+2,'x')
    Fpset = set(Fp)
    Fqset = set(Fq)
    Fp_squares = set(x**2 for x in Fpset)
    Fq_squares = set(x**2 for x in Fqset)

    # Pairs of squares, pairs of non-squares
    d = []
    d.extend(product(Fp_squares.difference([0]),Fq_squares.difference([0])))
    d.extend(product(Fpset.difference(Fp_squares),Fqset.difference(Fq_squares)))

    # All (x,0)
    d.extend((x,0) for x in Fpset)

    G = cartesian_product([Fp,Fq])

    if check and not is_difference_family(G, [d]):
        raise RuntimeError("twin_prime_powers_difference_set produced a wrong "
                           "difference set with p={}. Please contact "
                           "sage-devel@googlegroups.com".format(p))

    return G, [d]


def are_mcfarland_1973_parameters(v, k, lmbda, return_parameters=False):
    r"""
    Test whether ``(v,k,lmbda)`` is a triple that can be obtained from the
    construction from [McF1973]_.

    See :func:`mcfarland_1973_construction`.

    INPUT:

    - ``v``, ``k``, ``lmbda`` -- integers; parameters of the difference family
    - ``return_parameters`` -- boolean (default: ``False``); if ``True``, return a
      pair ``(True, (q, s))`` so that ``(q,s)`` can be used in the function
      :func:`mcfarland_1973_construction` to actually build a
      ``(v,k,lmbda)``-difference family. Or ``(False, None)`` if the
      construction is not possible

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: from sage.combinat.designs.difference_family import are_mcfarland_1973_parameters
        sage: are_mcfarland_1973_parameters(64, 28, 12)
        True
        sage: are_mcfarland_1973_parameters(64, 28, 12, return_parameters=True)
        (True, (2, 2))
        sage: are_mcfarland_1973_parameters(60, 13, 5)
        False
        sage: are_mcfarland_1973_parameters(98125, 19500, 3875)
        True
        sage: are_mcfarland_1973_parameters(98125, 19500, 3875, True)
        (True, (5, 3))

        sage: from sage.combinat.designs.difference_family import are_mcfarland_1973_parameters
        sage: for v in range(1, 100):                                                   # needs sage.rings.finite_rings
        ....:     for k in range(1,30):
        ....:         for l in range(1,15):
        ....:             if are_mcfarland_1973_parameters(v,k,l):
        ....:                 answer, (q,s) = are_mcfarland_1973_parameters(v,k,l,return_parameters=True)
        ....:                 print("{} {} {} {} {}".format(v,k,l,q,s))
        ....:                 assert answer is True
        ....:                 assert designs.difference_family(v,k,l,existence=True) is True
        ....:                 G,D = designs.difference_family(v,k,l)
        16 6 2 2 1
        45 12 3 3 1
        64 28 12 2 2
        96 20 4 4 1
    """
    if v <= k or k <= lmbda:
        return (False, None) if return_parameters else False
    k = ZZ(k)
    lmbda = ZZ(lmbda)
    qs, r = (k - lmbda).sqrtrem()  # sqrt(k-l) should be q^s
    if r or (qs*(qs-1)) % lmbda:
        return (False, None) if return_parameters else False

    q = qs*(qs-1) // lmbda + 1
    if (q <= 1 or
        v * (q-1) != qs*q * (qs*q+q-2) or
        k * (q-1) != qs * (qs*q-1)):
        return (False, None) if return_parameters else False

    # NOTE: below we compute the value of s so that qs = q^s. If the method
    # is_power_of of integers would be able to return the exponent, we could use
    # that... but currently this is not the case
    # see github issue #19792
    p1,a1 = qs.is_prime_power(get_data=True)
    p2,a2 = q.is_prime_power(get_data=True)

    if a1 == 0 or a2 == 0 or p1 != p2 or a1 % a2:
        return (False, None) if return_parameters else False

    return (True, (q, a1//a2)) if return_parameters else True


def mcfarland_1973_construction(q, s):
    r"""
    Return a difference set.

    The difference set returned has the following parameters

    .. MATH::

        v = \frac{q^{s+1}(q^{s+1}+q-2)}{q-1},
        k = \frac{q^s (q^{s+1}-1)}{q-1},
        \lambda = \frac{q^s(q^s-1)}{q-1}

    This construction is due to [McF1973]_.

    INPUT:

    - ``q``, ``s`` -- integers; parameters for the difference set (see the above
      formulas for the expression of ``v``, ``k``, ``l`` in terms of ``q`` and
      ``s``)

    .. SEEALSO::

        The function :func:`are_mcfarland_1973_parameters` makes the translation
        between the parameters `(q,s)` corresponding to a given triple
        `(v,k,\lambda)`.

    REFERENCES:

    .. [McF1973] Robert L. McFarland
       "A family of difference sets in non-cyclic groups"
       J. Combinatorial Theory (A) 15 (1973) 1--10.
       :doi:`10.1016/0097-3165(73)90031-9`

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import (
        ....:    mcfarland_1973_construction, is_difference_family)

        sage: G,D = mcfarland_1973_construction(3, 1)                                   # needs sage.modules
        sage: assert is_difference_family(G, D, 45, 12, 3)                              # needs sage.modules

        sage: G,D = mcfarland_1973_construction(2, 2)                                   # needs sage.modules
        sage: assert is_difference_family(G, D, 64, 28, 12)                             # needs sage.modules
    """
    from sage.rings.finite_rings.finite_field_constructor import GF
    from sage.modules.free_module import VectorSpace
    from sage.rings.finite_rings.integer_mod_ring import Zmod
    from sage.categories.cartesian_product import cartesian_product

    r = (q**(s+1)-1) // (q-1)
    F = GF(q,'a')
    V = VectorSpace(F, s+1)
    K = Zmod(r+1)

    G = cartesian_product([F]*(s+1) + [K])

    D = []
    for k, H in zip(K, V.subspaces(s)):
        for v in H:
            D.append(G(tuple(v) + (k,)))

    return G,[D]


def are_hadamard_difference_set_parameters(v, k, lmbda):
    r"""
    Check whether ``(v,k,lmbda)`` is of the form ``(4N^2, 2N^2 - N, N^2 - N)``.

    INPUT:

    - ``(v, k, lmbda)`` -- parameters of a difference set

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import are_hadamard_difference_set_parameters
        sage: are_hadamard_difference_set_parameters(36, 15, 6)
        True
        sage: are_hadamard_difference_set_parameters(60, 13, 5)
        False
    """
    N = k - 2*lmbda
    N2 = N*N
    return v == 4*N2 and k == 2*N2 - N and lmbda == N2 - N


@cached_function
def hadamard_difference_set_product_parameters(N):
    r"""
    Check whether a product construction is available for Hadamard difference
    set with parameter ``N``.

    This function looks for two integers `N_1` and `N_2` greater than `1`
    and so that `N = 2 N_1 N_2` and there exists Hadamard difference set with
    parameters `(4 N_i^2, 2N_i^2 - N_i, N_i^2 - N_i)`. If such pair exists,
    the output is the pair ``(N_1, N_2)`` otherwise it is ``None``.

    INPUT:

    - ``N`` -- positive integer

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import hadamard_difference_set_product_parameters
        sage: hadamard_difference_set_product_parameters(8)                             # needs sage.rings.finite_rings
        (2, 2)
    """
    if N % 2:
        return False

    for N1 in (N//2).divisors()[1:]:
        if 4*N1 > N:
            break
        v1 = 4*N1*N1
        k1 = 2*N1*N1 - N1
        l1 = N1*N1 - N1
        if not difference_family(v1, k1, l1, existence=True):
            continue
        N2 = N // (2*N1)
        v2 = 4*N2*N2
        k2 = 2*N2*N2 - N2
        l2 = N2*N2 - N2
        if not difference_family(v2, k2, l2, existence=True):
            continue

        return (N1,N2)

    return None


def hadamard_difference_set_product(G1, D1, G2, D2):
    r"""
    Make a product of two Hadamard difference sets.

    This product construction appears in [Tu1984]_.

    INPUT:

    - ``G1, D1``, ``G2, D2`` -- two Hadamard difference sets

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import hadamard_difference_set_product
        sage: from sage.combinat.designs.difference_family import is_difference_family

        sage: G1,D1 = designs.difference_family(16,6,2)                                 # needs sage.rings.finite_rings
        sage: G2,D2 = designs.difference_family(36,15,6)                                # needs sage.rings.finite_rings

        sage: G11,D11 = hadamard_difference_set_product(G1,D1,G1,D1)                    # needs sage.rings.finite_rings
        sage: assert is_difference_family(G11, D11, 256, 120, 56)                       # needs sage.rings.finite_rings
        sage: assert designs.difference_family(256, 120, 56, existence=True) is True    # needs sage.rings.finite_rings

        sage: G12,D12 = hadamard_difference_set_product(G1,D1,G2,D2)                    # needs sage.rings.finite_rings
        sage: assert is_difference_family(G12, D12, 576, 276, 132)                      # needs sage.rings.finite_rings
        sage: assert designs.difference_family(576, 276, 132, existence=True) is True   # needs sage.rings.finite_rings
    """
    from sage.categories.cartesian_product import cartesian_product

    G = cartesian_product([G1,G2])
    D1 = set(D1[0])
    D1c = set(s for s in G1 if s not in D1)
    D2 = set(D2[0])
    D2c = set(s for s in G2 if s not in D2)

    D = set().union((G((s1,s2)) for s1 in D1 for s2 in D2),
                    (G((s1,s2)) for s1 in D1c for s2 in D2c))

    return G, [[s for s in G if s not in D]]


def turyn_1965_3x3xK(k=4):
    r"""
    Return a difference set in either `C_3 \times C_3 \times C_4` or `C_3 \times
    C_3 \times C_2 \times C_2` with parameters `v=36`, `k=15`, `\lambda=6`.

    This example appears in [Tu1965]_.

    INPUT:

    - ``k`` -- either ``2`` (to get a difference set in `C_3 \times C_3 \times
      C_2 \times C_2`) or ``4`` (to get a difference set in `C_3 \times C_3
      \times C_3 \times C_4`)

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import turyn_1965_3x3xK
        sage: from sage.combinat.designs.difference_family import is_difference_family
        sage: G,D = turyn_1965_3x3xK(4)
        sage: assert is_difference_family(G, D, 36, 15, 6)
        sage: G,D = turyn_1965_3x3xK(2)
        sage: assert is_difference_family(G, D, 36, 15, 6)
    """
    from sage.categories.cartesian_product import cartesian_product
    from sage.rings.finite_rings.integer_mod_ring import Zmod

    if k == 2:
        G = cartesian_product([Zmod(3), Zmod(3), Zmod(2), Zmod(2)])
        K = [(0,0), (0,1), (1,0), (1,1)]
    elif k == 4:
        G = cartesian_product([Zmod(3), Zmod(3), Zmod(4)])
        K = [(0,), (1,), (2,), (3,)]
    else:
        raise ValueError("k must be 2 or 4")

    L = [[(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)],  # complement of y=0
         [(0,0),(1,1),(2,2)],                    # x-y=0
         [(0,0),(1,2),(2,1)],                    # x+y=0
         [(0,0),(0,1),(0,2)]]                    # x=0

    return G, [[G(v + k) for l, k in zip(L, K) for v in l]]


def _is_periodic_sequence(seq, period):
    r"""
    Check if the sequence is periodic with correct period.

    The sequence should have length at least twice the period, so that
    periodicity can be checked.

    INPUT:

    - ``seq`` -- the sequence to be tested (must have length at least twice the period)
    - ``period`` -- integer; the period that the sequence should have

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import _is_periodic_sequence
        sage: _is_periodic_sequence([0, 1, 2, 3, 0, 1, 2, 3], 4)
        True
        sage: _is_periodic_sequence([0, 1, 0, 1, 0, 1, 0, 1], 4)
        False
        sage: _is_periodic_sequence([0, 1, 1, 1, 0, 1, 2, 1], 4)
        False
    """
    assert len(seq) >= 2*period

    for per in range(1, period):
        first = seq[:per]
        periodic = True
        for j in range(1, len(seq)//per):
            if seq[j*per : (j+1)*per] != first:
                periodic = False
                break
        if periodic:
            return False
    if seq[:period] != seq[period : 2*period]:
        return False
    return True


def _create_m_sequence(q, n, check=True):
    r"""
    Create an m-sequence over GF(q) with period `q^n - 1`.

    Given a prime power `q`, the m-sequence is created as described by [Zie1959]_
    from a primitive function over the finite field `GF(q)`.

    Given a primitive function `f=c_0+c_1x+...+c_nx^n` over `K = GF(q)` of degree `n`,
    the recurrence is given by: `a_i = -c_0^{-1}(c_1a_{i-1} + ... + c_na{i-n})`.
    The first `n` elements will be `0, 0, ..., 0, 1` and these will give a maximal length recurrence sequence
    as shown in [Mit2008]_.

    INPUT:

    - ``q`` -- a prime power
    - ``n`` -- a nonnegative number
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the
      result is a sequence with correct period; detting it to ``False`` may
      speed up considerably the computation

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import _create_m_sequence
        sage: _create_m_sequence(3, 2)  # random                                        # needs sage.rings.finite_rings
        [1, 0, 1, 2, 2, 0, 2, 1]
        sage: _create_m_sequence(4, 2, check=False)  # random                           # needs sage.rings.finite_rings
        [1, 0, a, a + 1, a, a, 0, a + 1, 1, a + 1, a + 1, 0, 1, a, 1]
        sage: _create_m_sequence(6, 2)
        Traceback (most recent call last):
        ...
        ValueError: q must be a prime power
    """
    from sage.rings.finite_rings.finite_field_constructor import GF

    if not is_prime_power(q):
        raise ValueError('q must be a prime power')
    if n < 0:
        raise ValueError('n cannot be negative')

    K = GF(q, 'a')

    T = PolynomialRing(K, 'x')
    primitive = T.irreducible_element(n, algorithm='random')
    while not primitive.is_primitive():
        primitive = T.irreducible_element(n, algorithm='random')
    coeffs = primitive.coefficients()
    exps = primitive.exponents()

    period = q**n - 1
    seq_len = period*2 if check else period
    seq = [1] + [0]*(n-1)

    while len(seq) < seq_len:
        nxt = 0
        for i, coeff in zip(exps[1:], coeffs[1:]):
            nxt += coeff * seq[-i]
        seq.append(-coeffs[0].inverse() * nxt)

    if check:
        assert _is_periodic_sequence(seq, period)
    return seq[:period]


def _get_submodule_of_order(G, order):
    r"""
    Construct a submodule of the given order from group ``G``.

    This method tries to construct submodules from various elements of `G` until
    a submodule of the correct order is found.

    INPUT:

    - ``G`` -- an additive abelian group
    - ``order`` -- integer; the order of the desired submodule

    TESTS:

        sage: # needs sage.modules
        sage: from sage.combinat.designs.difference_family import _get_submodule_of_order
        sage: G = AdditiveAbelianGroup([48])
        sage: _get_submodule_of_order(G, 6).order()
        6
        sage: G = AdditiveAbelianGroup([13^2 - 1])
        sage: _get_submodule_of_order(G, 12).order()
        12
    """
    for el in G:
        H = G.submodule([el])
        if H.order() == order:
            return H
    return None


def relative_difference_set_from_m_sequence(q, N, check=True, return_group=False):
    r"""
    Construct `R((q^N-1)/(q-1), q-1, q^{N-1}, q^{N-2})` where ``q`` is a prime power and `N\ge 2`.

    The relative difference set is constructed over the set of additive integers modulo `q^N-1`,
    as described in Theorem 5.1 of [EB1966]_. Given an m-sequence `(a_i)` of period `q^N-1`, the
    set is: `R=\{i | 0 \le i \le q^{N-1}, a_i=1\}`.

    INPUT:

    - ``q`` -- a prime power
    - ``N`` -- a nonnegative number
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the
      result is a relative difference set before returning it
    - ``return_group`` -- boolean (default: ``False``); if ``True``, the function
      will also return the group from which the set is created

    OUTPUT:

    If ``return_group=False``, the function return only the relative difference
    set. Otherwise, it returns a tuple containing the group and the set.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import relative_difference_set_from_m_sequence
        sage: relative_difference_set_from_m_sequence(2, 4,               # random      # needs sage.modules sage.rings.finite_rings
        ....:                                         return_group=True)
        (Additive abelian group isomorphic to Z/15,
         [(0), (4), (5), (6), (7), (9), (11), (12)])
        sage: relative_difference_set_from_m_sequence(8, 2, check=False)  # random      # needs sage.modules sage.rings.finite_rings
        [(0), (6), (30), (40), (41), (44), (56), (61)]
        sage: relative_difference_set_from_m_sequence(6, 2)                             # needs sage.modules
        Traceback (most recent call last):
        ...
        ValueError: q must be a prime power

    TESTS::

        sage: from sage.combinat.designs.difference_family import is_relative_difference_set, _get_submodule_of_order
        sage: q, N = 5, 3
        sage: G, D = relative_difference_set_from_m_sequence(q, N, check=False,         # needs sage.modules sage.rings.finite_rings
        ....:                                                return_group=True)
        sage: H = _get_submodule_of_order(G, q-1)                                       # needs sage.modules sage.rings.finite_rings
        sage: is_relative_difference_set(D, G, H,                                       # needs sage.modules sage.rings.finite_rings
        ....:                            ((q^N-1)//(q-1), q-1, q^(N-1), q^(N-2)))
        True
        sage: q, N = 13, 2
        sage: G, D = relative_difference_set_from_m_sequence(q, N, check=False,         # needs sage.modules sage.rings.finite_rings
        ....:                                                return_group=True)
        sage: H = _get_submodule_of_order(G, q-1)                                       # needs sage.modules sage.rings.finite_rings
        sage: is_relative_difference_set(D, G, H,                                       # needs sage.modules sage.rings.finite_rings
        ....:                            ((q^N-1)//(q-1), q-1, q^(N-1), q^(N-2)))
        True
    """
    from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup

    if not is_prime_power(q):
        raise ValueError('q must be a prime power')
    if N < 2:
        raise ValueError('N must be at least 2')

    m_seq = _create_m_sequence(q, N, check=False)
    period = q**N - 1
    G = AdditiveAbelianGroup([period])

    set1 = [i for i in G if m_seq[i[0]] == 1]

    if check:
        H = _get_submodule_of_order(G, q-1)
        assert is_relative_difference_set(set1, G, H, (period // (q-1), q - 1, q**(N-1), q**(N-2)))

    if return_group:
        return G, set1
    return set1


def relative_difference_set_from_homomorphism(q, N, d, check=True, return_group=False):
    r"""
    Construct `R((q^N-1)/(q-1), n, q^{N-1}, q^{N-2}d)` where `nd = q-1`.

    Given a prime power `q`, a number `N \ge 2` and integers `d` such that `d | q-1` we create the
    relative difference set using the construction from Corollary 5.1.1 of [EB1966]_.

    INPUT:

    - ``q`` -- a prime power
    - ``N`` -- integer greater than 1
    - ``d`` -- integer which divides `q-1`
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the
      result is a relative difference set before returning it
    - ``return_group`` -- boolean (default: ``False``); if ``True``, the function
      will also return the group from which the set is created

    OUTPUT:

    If ``return_group=False``, the function return only the relative difference
    set. Otherwise, it returns a tuple containing the group and the set.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import relative_difference_set_from_homomorphism
        sage: relative_difference_set_from_homomorphism(7, 2, 3)        # random        # needs sage.modules sage.rings.finite_rings
        [(0), (3), (4), (2), (13), (7), (14)]
        sage: relative_difference_set_from_homomorphism(9, 2, 4,        # random        # needs sage.modules sage.rings.finite_rings
        ....:                                           check=False, return_group=True)
        (Additive abelian group isomorphic to Z/80,
         [(0), (4), (6), (13), (7), (12), (15), (8), (9)])
        sage: relative_difference_set_from_homomorphism(9, 2, 5)                        # needs sage.modules sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        ValueError: q-1 must be a multiple of d

    TESTS::

        sage: from sage.combinat.designs.difference_family import is_relative_difference_set, _get_submodule_of_order
        sage: q, N, d = 11, 2, 5
        sage: G, D = relative_difference_set_from_homomorphism(q, N, d, check=False,    # needs sage.modules sage.rings.finite_rings
        ....:                                                  return_group=True)
        sage: H = _get_submodule_of_order(G, (q-1)//d)                                  # needs sage.modules sage.rings.finite_rings
        sage: is_relative_difference_set(D, G, H,                                       # needs sage.modules sage.rings.finite_rings
        ....:                            ((q**N-1)//(q-1), (q-1)//d, q**(N-1), q**(N-2)*d))
        True
        sage: q, N, d = 9, 2, 4
        sage: G, D = relative_difference_set_from_homomorphism(q, N, d, check=False,    # needs sage.modules sage.rings.finite_rings
        ....:                                                  return_group=True)
        sage: H = _get_submodule_of_order(G, (q-1)//d)                                  # needs sage.modules sage.rings.finite_rings
        sage: is_relative_difference_set(D, G, H,                                       # needs sage.modules sage.rings.finite_rings
        ....:                            ((q**N-1)//(q-1), (q-1)//d, q**(N-1), q**(N-2)*d))
        True
    """
    from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup

    if not is_prime_power(q):
        raise ValueError('q must be a prime power')
    if N < 2:
        raise ValueError('N must be at least 2')
    if (q-1) % d != 0:
        raise ValueError('q-1 must be a multiple of d')

    G = AdditiveAbelianGroup([q**N - 1])
    K = _get_submodule_of_order(G, d)
    assert K is not None, 'Could not find kernel'

    G2 = G/K

    theta = G.hom([G2.gen(0)], G2)
    diff_set = relative_difference_set_from_m_sequence(q, N, check=False)
    second_diff_set = [theta(x) for x in diff_set]

    if check:
        H = _get_submodule_of_order(G2, (q-1) // d)
        assert is_relative_difference_set(second_diff_set, G2, H, ((q**N-1) // (q-1), (q-1) // d, q**(N-1), q**(N-2) * d))

    if return_group:
        return G2, second_diff_set
    return second_diff_set


def is_relative_difference_set(R, G, H, params, verbose=False):
    r"""
    Check if ``R`` is a difference set of ``G`` relative to ``H``, with the given parameters.

    This function checks that `G`, `H` and `R` have the orders specified in the parameters, and
    that `R` satisfies the definition of relative difference set (from [EB1966]_): the collection of
    differences `r-s`, `r,s \in R`, `r \neq s` contains only elements of `G` which are not in `H`, and contains
    every such element exactly `d` times.

    INPUT:

    - ``R`` -- list; the relative diffeence set of length `k`
    - ``G`` -- an additive abelian group of order `mn`
    - ``H`` -- list; a submodule of ``G`` of order `n`
    - ``params`` -- tuple in the form `(m, n, k, d)`
    - ``verbose`` -- boolean (default: ``False``); if ``True``, the function
      will be verbose when the sequences do not satisfy the constraints

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import _get_submodule_of_order, relative_difference_set_from_m_sequence, is_relative_difference_set
        sage: q, N = 5, 2
        sage: params = ((q^N-1) // (q-1), q - 1, q^(N-1), q^(N-2))
        sage: G, R = relative_difference_set_from_m_sequence(q, N, return_group=True)   # needs sage.libs.pari sage.modules
        sage: H = _get_submodule_of_order(G, q - 1)                                     # needs sage.libs.pari sage.modules
        sage: is_relative_difference_set(R, G, H, params)                               # needs sage.libs.pari sage.modules
        True

    If we pass the ``verbose`` argument, the function will explain why it failed::

        sage: R2 = [G[1], G[2], G[3], G[5], G[6]]                                       # needs sage.libs.pari sage.modules
        sage: is_relative_difference_set(R2, G, H, params, verbose=True)                # needs sage.libs.pari sage.modules
        There is a value in the difference set which is not repeated d times
        False
    """
    m, n, k, d = params
    if G.order() != m * n:
        if verbose:
            print('Incorrect order of G:', G.order())
        return False

    if H.order() != n:
        if verbose:
            print('Incorect order of H:', H.order())

    if len(R) != k:
        if verbose:
            print('Length of R not correct:', len(R))
        return False

    diff_set = {}
    for el1 in R:
        for el2 in R:
            if el1 != el2:
                idx = el1 - el2
                if idx not in diff_set:
                    diff_set[idx] = 0
                diff_set[idx] += 1
    values = [diff_set[x] for x in diff_set]
    if max(values) != d or min(values) != d:
        if verbose:
            print('There is a value in the difference set which is not repeated d times')
        return False

    for el in G:
        if el in H and el in diff_set:
            if verbose:
                print('An element of G is present in both the difference set and in H')
            return False
        if el not in H and el not in diff_set:
            if verbose:
                print('An element of G is not present in either one of H or the difference set')
            return False

    return True


def is_supplementary_difference_set(Ks, v=None, lmbda=None, G=None, verbose=False):
    r"""
    Check that the sets in ``Ks`` are `n-\{v; k_1, ..., k_n; \lambda \}` supplementary
    difference sets over group ``G`` of order ``v``.

    From the definition in [Spe1975]_: let  `S_1, S_2, ..., S_n` be `n` subsets of a group `G` of order `v`
    such that `|S_i| = k_i`. If, for each `g \in G`, `g \neq 0`, the total number of solutions of `a_i - a'_i = g`, with
    `a_i, a'_i \in S_i` is `\lambda`, then `S_1, S_2, ..., S_n` are `n-\{v; k_1, ..., k_n; \lambda\}` supplementary difference sets.

    One of the parameters ``v`` or ``G`` must always be specified. If ``G`` is not
    given, the function will use an ``AdditiveAbelianGroup`` of order ``v``.

    INPUT:

    - ``Ks`` -- list of sets to be checked
    - ``v`` -- integer; the parameter `v` of the supplementary difference sets
    - ``lmbda`` -- integer; the parameter `\lambda` of the supplementary difference sets
    - ``G`` -- a group of order `v`
    - ``verbose`` -- boolean (default: ``False``); if ``True``, the function will
      be verbose when the sets do not satisfy the constraints

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import supplementary_difference_set_from_rel_diff_set, is_supplementary_difference_set
        sage: G, [S1, S2, S3, S4] = supplementary_difference_set_from_rel_diff_set(17)  # needs sage.modules sage.rings.finite_rings
        sage: is_supplementary_difference_set([S1, S2, S3, S4], lmbda=16, G=G)          # needs sage.modules sage.rings.finite_rings
        True

    The parameter ``v`` can be given instead of ``G``::

        sage: is_supplementary_difference_set([S1, S2, S3, S4], v=16, lmbda=16)         # needs sage.modules sage.rings.finite_rings
        True
        sage: is_supplementary_difference_set([S1, S2, S3, S4], v=20, lmbda=16)         # needs sage.modules sage.rings.finite_rings
        False

    If ``verbose=True``, the function will be verbose::

        sage: is_supplementary_difference_set([S1, S2, S3, S4], lmbda=14, G=G,          # needs sage.modules sage.rings.finite_rings
        ....:                                 verbose=True)
        Number of pairs with difference (1) is 16, but lambda is 14
        False

    TESTS::

        sage: # needs sage.modules sage.rings.finite_rings
        sage: is_supplementary_difference_set([[1], [1]], lmbda=0, G=Zmod(3))
        True
        sage: is_supplementary_difference_set([S1, S2, S3, S4], v=17, lmbda=16, G=G)
        False
        sage: is_supplementary_difference_set([S1, S2, S3, S4], G=G)
        True
        sage: is_supplementary_difference_set([S1, S2, S3, S4], lmbda=16)
        Traceback (most recent call last):
        ...
        ValueError: one of G or v must be specified

    .. SEEALSO::

        :func:`supplementary_difference_set_from_rel_diff_set`
    """

    if G is None and v is None:
        raise ValueError('one of G or v must be specified')

    if G is None:
        from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup
        G = AdditiveAbelianGroup([v])

    if v is not None and G.order() != v:
        if verbose:
            print(f'G has order {G.order()}, but it should be {v}')
        return False

    differences_counter = {el: 0 for el in G}
    for K in Ks:
        for el1 in K:
            for el2 in K:
                diff = G(el1) - G(el2)
                differences_counter[diff] += 1

    for key, diff in differences_counter.items():
        if key == 0:
            continue
        if lmbda is None:
            lmbda = diff
        if diff != lmbda:
            if verbose:
                print(f'Number of pairs with difference {key} is {diff}, but lambda is {lmbda}')
            return False

    return True


def supplementary_difference_set_from_rel_diff_set(q, existence=False, check=True):
    r"""
    Construct `4-\{2v; v, v+1, v, v; 2v\}` supplementary difference sets where `q=2v+1`.

    The sets are created from relative difference sets as detailed in Theorem 3.3 of [Spe1975]_. this construction
    requires that `q` is an odd prime power and that there exists `s \ge 0` such that `(q-(2^{s+1}+1))/2^{s+1}` is
    an odd prime power.

    Note that the construction from [Spe1975]_ states that the resulting sets are `4-\{2v; v+1, v, v, v; 2v\}`
    supplementary difference sets. However, the implementation of that construction returns
    `4-\{2v; v, v+1, v, v; 2v\}` supplementary difference sets. This is not important, since the supplementary
    difference sets are not ordered.

    INPUT:

    - ``q`` -- an odd prime power
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check
      whether the supplementary difference sets can be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the sets
      are supplementary difference sets before returning them

    OUTPUT:

    If ``existence=False``, the function returns the 4 sets (containing integers),
    or raises an error if ``q`` does not satisfy the constraints.
    If ``existence=True``, the function returns a boolean representing whether
    supplementary difference sets can be constructed.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import supplementary_difference_set_from_rel_diff_set
        sage: supplementary_difference_set_from_rel_diff_set(17)  #random               # needs sage.libs.pari
        (Additive abelian group isomorphic to Z/16,
         [[(1), (5), (6), (7), (9), (13), (14), (15)],
          [(0), (2), (3), (5), (6), (10), (11), (13), (14)],
          [(0), (1), (2), (3), (5), (6), (7), (12)],
          [(0), (2), (3), (5), (6), (7), (9), (12)]])

    If ``existence=True``, the function returns a boolean::

        sage: supplementary_difference_set_from_rel_diff_set(7, existence=True)
        False
        sage: supplementary_difference_set_from_rel_diff_set(17, existence=True)
        True

    TESTS::

        sage: # needs sage.libs.pari
        sage: from sage.combinat.designs.difference_family import is_supplementary_difference_set
        sage: G, sets = supplementary_difference_set_from_rel_diff_set(17, check=False)
        sage: is_supplementary_difference_set(sets, lmbda=16, G=G)
        True
        sage: G, sets = supplementary_difference_set_from_rel_diff_set(9, check=False)
        sage: is_supplementary_difference_set(sets, lmbda=8, G=G)
        True
        sage: supplementary_difference_set_from_rel_diff_set(7)
        Traceback (most recent call last):
        ...
        ValueError: There is no s for which m-1 is an odd prime power
        sage: supplementary_difference_set_from_rel_diff_set(8)
        Traceback (most recent call last):
        ...
        ValueError: q must be an odd prime power
        sage: supplementary_difference_set_from_rel_diff_set(8, existence=True)
        False
        sage: supplementary_difference_set_from_rel_diff_set(7, existence=True)
        False
        sage: supplementary_difference_set_from_rel_diff_set(1, existence=True)
        False

    Check that the function works even when s > 1::

        sage: G, sets = supplementary_difference_set_from_rel_diff_set(353, check=False)  # long time, needs sage.libs.pari
        sage: is_supplementary_difference_set(sets, lmbda=352, G=G)     # long time, needs sage.libs.pari
        True

    .. SEEALSO::

        :func:`is_supplementary_difference_set`
    """
    s = 0
    m = -1

    while q > 2**(s+1) and (q-1) % 2**(s+1) == 0:
        prime_pow = (q-1)//2**(s+1) - 1
        if is_prime_power(prime_pow) and prime_pow % 2 == 1:
            m = (q - (2**(s+1) + 1)) // 2**(s+1) + 1
            break
        s += 1

    if existence:
        return is_prime_power(q) and q % 2 == 1 and m != -1

    if not is_prime_power(q) or q % 2 != 1:
        raise ValueError('q must be an odd prime power')
    if m == -1:
        raise ValueError('There is no s for which m-1 is an odd prime power')

    set1 = relative_difference_set_from_homomorphism(m - 1, 2, (m-2) // 2, check=False)

    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    P = PolynomialRing(ZZ, 'x')

    # Compute psi3, psi4
    hall = 0
    for d in set1:
        hall += P.monomial(d[0])

    def get_T(k):
        T = P.monomial(0) - 1
        for i in range(k):
            T += P.monomial(i)
        return T

    modulo = P.monomial(2*m) - 1

    diff = get_T(2*m) - (1+P.monomial(m))*hall
    diff = diff.mod(modulo)
    exp1, exp2 = diff.exponents()
    a = (exp1+exp2-m) // 2

    psi3 = (P.monomial(a) + hall).mod(modulo)
    psi4 = (P.monomial(a+m) + hall).mod(modulo)

    for i in range(s):
        m_start = 2**i * m
        psi3, psi4 = (psi3(P.monomial(2)) + P.monomial(1)*psi4(P.monomial(2))).mod(P.monomial(4*m_start)-1), \
                 (psi3(P.monomial(2)) + P.monomial(1)*(get_T(2*m_start)(P.monomial(2)) - psi4(P.monomial(2)))).mod(P.monomial(4*m_start)-1)

    # Construction of psi1, psi2
    G2, set2 = relative_difference_set_from_m_sequence(q, 2, check=False, return_group=True)
    s3 = get_fixed_relative_difference_set(G2, set2)

    phi_exps = []
    for i in range(len(s3)):
        for j in range(i+1, len(s3)):
            diff = s3[i] - s3[j]
            if diff % (q-1) == 0 and diff % (q**2-1) != 0:
                phi_exps.append(s3[i])

    exps1 = [(x+1)//2 for x in phi_exps if x % 2 == 1]
    exps2 = [x//2 for x in phi_exps if x % 2 == 0]

    theta1 = 0
    for exp in exps1:
        theta1 += P.monomial(exp)
    theta1 = theta1.mod(P.monomial(q-1)-1)

    theta2 = 0
    for exp in exps2:
        theta2 += P.monomial(exp)
    theta2 = theta2.mod(P.monomial(q-1) - 1)

    psi1 = ((1 + P.monomial((q-1)//2)) * theta1).mod(P.monomial(q-1) - 1)
    psi2 = (1 + (1 + P.monomial((q-1)//2)) * theta2).mod(P.monomial(q-1) - 1)

    from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup
    G = AdditiveAbelianGroup([q-1])
    K1 = [G[x] for x in psi1.exponents()]
    K2 = [G[x] for x in psi2.exponents()]
    K3 = [G[x] for x in psi3.exponents()]
    K4 = [G[x] for x in psi4.exponents()]

    if check:
        assert is_supplementary_difference_set([K1, K2, K3, K4], lmbda=q-1, G=G)

    return G, [K1, K2, K3, K4]


def supplementary_difference_set(q, existence=False, check=True):
    r"""
    Construct `4-\{2v; v, v+1, v, v; 2v\}` supplementary difference sets where `q=2v+1`.

    This is a deprecated version of :func:`supplementary_difference_set_from_rel_diff_set`,
    please use that instead.
    """
    from sage.misc.superseded import deprecation
    deprecation(35211, 'This function is deprecated, please use supplementary_difference_set_from_rel_diff_set instead.')

    if existence:
        return supplementary_difference_set_from_rel_diff_set(q, existence=True)
    _, s = supplementary_difference_set_from_rel_diff_set(q, check=check)
    return s


def get_fixed_relative_difference_set(G, rel_diff_set, as_elements=False):
    r"""
    Construct an equivalent relative difference set fixed by the size of the set.

    Given a relative difference set `R(q+1, q-1, q, 1)`, it is possible to find a translation
    of this set fixed by `q` (see Section 3 of [Spe1975]_). We say that a set is fixed by `t` if
    `\{td | d\in R\}= R`.

    In addition, the set returned by this function will contain the element `0`. This is needed in the
    construction of supplementary difference sets (see :func:`supplementary_difference_set_from_rel_diff_set`).

    INPUT:

    - ``G`` -- a group, of which ``rel_diff_set`` is a subset
    - ``rel_diff_set`` -- the relative difference set
    - ``as_elements`` -- boolean (default: ``False``); if ``True``, the list
      returned will contain elements of the abelian group (this may slow down
      the computation considerably)

    OUTPUT:

    By default, this function returns the set as a list of integers. However, if
    ``as_elements=True`` it will return the set as a list containing elements of
    the abelian group.
    If no such set can be found, the function will raise an error.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import relative_difference_set_from_m_sequence, get_fixed_relative_difference_set
        sage: G, s1 = relative_difference_set_from_m_sequence(5, 2, return_group=True)  # needs sage.libs.pari sage.modules
        sage: get_fixed_relative_difference_set(G, s1)  # random                        # needs sage.libs.pari sage.modules
        [2, 10, 19, 23, 0]

    If ``as_elements=True``, the result will contain elements of the group::

        sage: get_fixed_relative_difference_set(G, s1, as_elements=True)  # random      # needs sage.libs.pari sage.modules
        [(2), (10), (19), (23), (0)]

    TESTS::

        sage: # needs sage.libs.pari sage.modules
        sage: from sage.combinat.designs.difference_family import is_fixed_relative_difference_set
        sage: G, s1 = relative_difference_set_from_m_sequence(5, 2, return_group=True)
        sage: s2 = get_fixed_relative_difference_set(G, s1, as_elements=True)
        sage: is_fixed_relative_difference_set(s2, len(s2))
        True
        sage: G, s1 = relative_difference_set_from_m_sequence(9, 2, return_group=True)
        sage: s2 = get_fixed_relative_difference_set(G, s1, as_elements=True)
        sage: is_fixed_relative_difference_set(s2, len(s2))
        True
        sage: type(s2[0])
        <class 'sage.groups.additive_abelian.additive_abelian_group.AdditiveAbelianGroup_fixed_gens_with_category.element_class'>
        sage: s2 = get_fixed_relative_difference_set(G, s1)
        sage: type(s2[0])
        <class 'sage.rings.integer.Integer'>
    """
    q = len(rel_diff_set)

    s2 = None
    for el in G:
        fixed_set = [el+x for x in rel_diff_set]
        if is_fixed_relative_difference_set(fixed_set, q):
            s2 = fixed_set
            break
    assert s2 is not None, 'Cannot find fixed translation of the set'

    s3 = None
    for i in range(G.order()):
        temp = [((q+1)*i+x[0]) % G.order() for x in s2]
        if 0 in temp:
            s3 = temp
            break
    assert s3 is not None, 'Cannot find fixed set containing 0'

    if as_elements:
        return [G[x] for x in s3]
    return s3


def is_fixed_relative_difference_set(R, q):
    r"""
    Check if the relative difference set ``R`` is fixed by ``q``.

    A relative difference set  `R` is fixed by `q` if  `\{qd | d \in R\}= R` (see Section 3 of [Spe1975]_).

    INPUT:

    - ``R`` -- list containing elements of an abelian group; the relative
      difference set
    - ``q`` -- integer

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.combinat.designs.difference_family import relative_difference_set_from_m_sequence, get_fixed_relative_difference_set, is_fixed_relative_difference_set
        sage: G, s1 = relative_difference_set_from_m_sequence(7, 2, return_group=True)  # needs sage.libs.pari
        sage: s2 = get_fixed_relative_difference_set(G, s1, as_elements=True)           # needs sage.libs.pari
        sage: is_fixed_relative_difference_set(s2, len(s2))                             # needs sage.libs.pari
        True
        sage: G = AdditiveAbelianGroup([15])
        sage: s3 = [G[1], G[2], G[3], G[4]]
        sage: is_fixed_relative_difference_set(s3, len(s3))
        False

    If the relative difference set does not contain elements of the group, the method returns false::

        sage: G, s1 = relative_difference_set_from_m_sequence(7, 2, return_group=True)  # needs sage.libs.pari sage.modules
        sage: s2 = get_fixed_relative_difference_set(G, s1, as_elements=False)          # needs sage.libs.pari sage.modules
        sage: is_fixed_relative_difference_set(s2, len(s2))                             # needs sage.libs.pari sage.modules
        False
    """
    for el in R:
        if q * el not in R:
            return False
    return True


def skew_supplementary_difference_set_over_polynomial_ring(n, existence=False, check=True):
    r"""
    Construct skew supplementary difference sets over a polynomial ring of order ``n``.

    The skew supplementary difference sets for `n=81, 169` are taken from [Djo1994a]_.

    INPUT:

    - ``n`` -- integer; the parameter of the supplementary difference sets
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check
      whether the supplementary difference sets can be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the sets
      are supplementary difference sets with `S_1` skew before returning them;
      setting this parameter to ``False`` may speed up the computation considerably

    OUTPUT:

    If ``existence=False``, the function returns a Polynomial Ring of order ``n``
    and a list containing 4 sets, or raises an error if data for the given ``n``
    is not available.
    If ``existence=True``, the function returns a boolean representing whether
    skew supplementary difference sets can be constructed.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import skew_supplementary_difference_set_over_polynomial_ring
        sage: G, [S1, S2, S3, S4] = skew_supplementary_difference_set_over_polynomial_ring(81)      # needs sage.libs.pari

    If ``existence=True``, the function returns a boolean::

        sage: skew_supplementary_difference_set_over_polynomial_ring(81, existence=True)
        True
        sage: skew_supplementary_difference_set_over_polynomial_ring(17, existence=True)
        False

    TESTS::

        sage: skew_supplementary_difference_set_over_polynomial_ring(7)
        Traceback (most recent call last):
        ...
        NotImplementedError: skew SDS of order 7 not yet implemented
    """
    data = {
        81: (3, lambda x: x**4 - x**3 - 1, 16, 5,
             [1, 2, 4, 6, 8, 10, 12, 14], [1, 2, 3, 4, 10, 11, 13],
             [4, 5, 6, 8, 12, 13, 14], [2, 4, 5, 6, 7, 11, 12, 13, 15]),
        169: (13, lambda x: x**2 - 4*x + 6, 24, 7,
              [0, 2, 5, 7, 9, 10, 12, 15, 16, 18, 21, 22], [0, 1, 2, 7, 8, 9, 13, 14, 18, 20, 23],
              [1, 4, 6, 7, 9, 14, 16, 17, 20, 21, 23], [3, 5, 6, 9, 10, 12, 13, 14, 15, 17, 20])
    }

    if existence:
        return n in data

    if n not in data:
        raise NotImplementedError(f'skew SDS of order {n} not yet implemented')

    mod, poly, exp, order, ind1, ind2, ind3, ind4 = data[n]

    from sage.rings.finite_rings.integer_mod_ring import Zmod

    Z3 = Zmod(mod)
    R = ZZ['x']
    x = R.gen()
    F = Z3.extension(poly(x))

    H = [F.gen() ** (exp * i) for i in range(order)]

    cosets = []
    for i in range((n - 1) // (2 * order)):
        cosets.append([F.gen()**i * el for el in H])
        cosets.append([-F.gen()**i * el for el in H])

    def generate_set(index_set, cosets):
        return sum((cosets[idx] for idx in index_set), [])

    S1 = generate_set(ind1, cosets)
    S2 = generate_set(ind2, cosets)
    S3 = generate_set(ind3, cosets)
    S4 = generate_set(ind4, cosets)

    if check:
        lmbda = len(S1) + len(S2) + len(S3) + len(S4) - n
        assert is_supplementary_difference_set([S1, S2, S3, S4], lmbda=lmbda, G=F)
        assert _is_skew_set(F, S1)

    return F, [S1, S2, S3, S4]


def skew_supplementary_difference_set_with_paley_todd(n, existence=False, check=True):
    r"""
    Construct `4-\{n; n_1, n_2, n_3, n_4; \lambda\}` skew supplementary difference sets where `S_1` is the Paley-Todd difference set.

    The skew SDS returned have the property that `n_1 + n_2 + n_3 + n_4 = n + \lambda`.

    This construction is described in [DK2016]_. The function contains, for each
    value of `n`, a set `H` containing integers modulo `n`, and four sets `J, K, L`.
    Then, these are used to construct `(n; k_2, k_3, k_4; \lambda_2)` difference family,
    with `\lambda_2 = k_2 + k_3 + k_4 + (3n - 1) / 4`. Finally, these sets together
    with the Paley-Todd difference set form a skew supplementary difference set.

    INPUT:

    - ``n`` -- integer; the parameter of the supplementary difference set
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check
      whether the supplementary difference sets can be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the sets
      are supplementary difference sets with `S_1` skew before returning them;
      setting this parameter to ``False`` may speed up the computation considerably

    OUTPUT:

    If ``existence=False``, the function returns the group G of integers modulo
    ``n`` and a list containing 4 sets, or raises an error if data for the given
    ``n`` is not available.
    If ``existence=True``, the function returns a boolean representing whether
    skew supplementary difference sets can be constructed.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import skew_supplementary_difference_set_with_paley_todd
        sage: G, [S1, S2, S3, S4] = skew_supplementary_difference_set_with_paley_todd(239)

    If existence is ``True``, the function returns a boolean::

        sage: skew_supplementary_difference_set_with_paley_todd(239, existence=True)
        True
        sage: skew_supplementary_difference_set_with_paley_todd(17, existence=True)
        False

    TESTS::

        sage: skew_supplementary_difference_set_with_paley_todd(7)
        Traceback (most recent call last):
        ...
        NotImplementedError: data for skew SDS of order 7 not yet implemented
    """
    H_db = {
        239: [1, 10, 24, 44, 98, 100, 201],
    }

    indices = {
        239: [[1, 3, 5, 6, 15, 17, 19, 28, 34, 38, 39, 57, 58, 63, 85, 95, 107],
              [1, 3, 4, 5, 15, 16, 17, 18, 19, 21, 23, 29, 35, 45, 58, 63],
              [0, 1, 4, 6, 7, 8, 13, 16, 18, 34, 35, 45, 47, 58, 63, 95]],
    }

    if existence:
        return n in H_db

    if n not in H_db:
        raise NotImplementedError(f'data for skew SDS of order {n} not yet implemented')

    G = Zmod(n)
    H = {G(el) for el in H_db[n]}

    def generate_subset(indices, H):
        return list({el * idx for el in H for idx in indices})

    from sage.arith.misc import quadratic_residues

    S1 = [G(el) for el in quadratic_residues(n) if el != 0]
    S2 = generate_subset(indices[n][0], H)
    S3 = generate_subset(indices[n][1], H)
    S4 = generate_subset(indices[n][2], H)

    if check:
        lmbda = len(S1) + len(S2) + len(S3) + len(S4) - n
        assert is_supplementary_difference_set([S1, S2, S3, S4], lmbda=lmbda, G=G)
        assert _is_skew_set(G, S1)

    return G, [S1, S2, S3, S4]


def _construct_gs_difference_family_from_full(S1, S2, mu):
    r"""
    Construct a spin type Goethals-Seidel difference family given the first two
    sets.

    This construction is described in [Djo2024]_. Given the first two sets
    `S_1, S_2`, and the multiplier `\mu`, the last two sets are computed as
    `S_3 = \mu S_2`, `S_4 =\mu S_3`.

    The sets should contain elements of the additive group of integers modulo `n`.

    INPUT:

    - ``S1`` -- list of integer modulo `n`; the first set
    - ``S2`` -- list of integer modulo `n`; the second set

    OUTPUT:

    The Goethals-Seidel difference family ``(S_1, S_2, S_3, S_4)``.

    TESTS::

        sage: from sage.combinat.designs.difference_family import _construct_gs_difference_family_from_full
        sage: G = Zmod(9)
        sage: S1 = list(map(G, [0,3,6]))
        sage: S2 = list(map(G, [0, 1, 8]))
        sage: _construct_gs_difference_family_from_full(S1, S2, 4)
        [[0, 3, 6], [0, 1, 8], [0, 4, 5], [0, 7, 2]]
    """
    S3 = [mu * el for el in S2]
    S4 = [mu * el for el in S3]
    return [S1, S2, S3, S4]


def _construct_gs_difference_family_from_compact(rep1, rep2, H, mu):
    r"""
    Construct a spin type Goethals-Seidel difference family given a compact
    representation of the first two sets.

    This construction is described in [Djo2024]_. Given a subgroup `H` of a group
    `G`, and two sets of representatives `rep_1, rep_2` we can construct the first
    two sets: `S_1 = \bigcup_{x \in rep_1} x H` and  `S_2 = \bigcup_{x \in rep_2} x H`.
    The list ``H`` should contain elements of the additive group of integers
    modulo `n`.

    The other two sets are constructed using
    :func:`_construct_gs_difference_family_from_full`.

    INPUT:

    - ``rep1`` -- list of integers; the first set of representatives
    - ``rep2`` -- list of integers; the second set of representatives
    - ``H`` -- list of integers modulo `n`; the subgroup used to generate the
      first two sets

    OUTPUT:

    The Goethals-Seidel difference family ``(S_1, S_2, S_3, S_4)``.

    TESTS::

        sage: from sage.combinat.designs.difference_family import _construct_gs_difference_family_from_compact, is_supplementary_difference_set
        sage: G = Zmod(61)
        sage: H = list(map(G, [1, 9, 20, 34, 58]))
        sage: rep1 = [3, 4, 5, 6, 8, 10]
        sage: rep2 = [0, 8, 10, 13, 23, 26]
        sage: [S1, S2, S3, S4] = _construct_gs_difference_family_from_compact(rep1, rep2, H, 13)
        sage: is_supplementary_difference_set([S1, S2, S3, S4], G=G)
        True
    """
    S1 = set()
    S2 = set()
    for el in H:
        S1 = S1.union({x * el for x in rep1})
        S2 = S2.union({x * el for x in rep2})
    S1 = list(S1)
    S2 = list(S2)
    return _construct_gs_difference_family_from_full(S1, S2, mu)


def spin_goethals_seidel_difference_family(n, existence=False, check=True):
    r"""
    Construct a spin type Goethals-Seidel difference family with parameters
    `(n; k_1, k_2, k_3, k_4; \lambda)`.

    The construction is described in [Djo2024]_. This function contains, for
    each value of `n`, either a full representation of `S_1, S_2` together with
    the multiplier `\mu`, or a subgroup `H`, two sets of representatives, and the
    multiplier.
    This data is used to construct the difference family using the functions
    :func:`_construct_gs_difference_family_from_full` and
    :func:`_construct_gs_difference_family_from_compact`.

    Additionally, this function also checks if a (skew) difference family can be
    constructed using :func:`skew_spin_goethals_seidel_difference_family`.

    INPUT:

    - ``n`` -- integer; the parameter of the GS difference family
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check
      whether the difference family can be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the sets
      are a difference family before returning them;
      setting this parameter to ``False`` may speed up the computation considerably

    OUTPUT:

    If ``existence=False``, the function returns the group G of integers modulo
    ``n`` and a list containing 4 sets, or raises an error if data for the given
    ``n`` is not available.
    If ``existence=True``, the function returns a boolean representing whether
    the difference family can be constructed.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import spin_goethals_seidel_difference_family
        sage: G, [S1, S2, S3, S4] = spin_goethals_seidel_difference_family(73)

    If existence is ``True``, the function returns a boolean::

        sage: spin_goethals_seidel_difference_family(73, existence=True)
        True
        sage: spin_goethals_seidel_difference_family(5, existence=True)
        False

    TESTS::

        sage: from sage.combinat.designs.difference_family import is_supplementary_difference_set
        sage: G, [S1, S2, S3, S4] = spin_goethals_seidel_difference_family(9, check=False)
        sage: lmbda = len(S1) + len(S2) + len(S3) + len(S4) - 9
        sage: is_supplementary_difference_set([S1, S2, S3, S4], lmbda=lmbda, G=G)
        True
        sage: spin_goethals_seidel_difference_family(5)
        Traceback (most recent call last):
        ...
        NotImplementedError: Data for spin type Goethals Seidel family of order 5 not yet implemented
    """
    full_data = {
        7: ([0], [0, 1, 6], 2),
        9: ([0, 3, 6], [0, 1, 8], 4),
        13: ([0, 1, 4, 6],  [0, 4, 6, 7, 9], 3),
        19: ([4, 6, 9, 10, 13, 15], [0, 1, 5, 8, 9, 10, 11, 13], 7),
        21: ([1, 4, 5, 8, 10, 11, 12, 17, 19], [1, 3, 8, 9, 12, 13, 18, 20], 4),
        31: ([2, 4, 6, 12, 14, 16, 17, 19, 25, 26, 28, 29],
             [0, 3, 9, 11, 13, 14, 15, 16, 17, 18, 20, 22, 28], 5),
        37: ([0, 3, 4, 5, 7, 13, 18, 19, 24, 30, 32, 33, 34],
             [0, 1, 2, 3, 4, 6, 12, 13, 18, 19, 24, 25, 31, 33, 34, 35, 36], 10),
        39: ([1, 4, 6, 10, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 29, 33, 35, 38],
             [0, 2, 4, 6, 7, 10, 11, 14, 16, 19, 20, 22, 26, 32, 33, 38], 16),
        57: ([0, 4, 11, 12, 18, 19, 20, 25, 26, 27, 28, 29, 30, 31, 32, 37, 38, 39, 45, 46, 53],
             [1, 2, 5, 6, 7, 8, 9, 10, 12, 17, 19, 21, 22, 24, 25, 28, 30, 31, 34, 37, 39, 41, 42, 43, 44, 46, 53, 54],
             7)
    }
    compact_data = {
        73: ([1, 8, 64],
             [0, 9, 13, 18, 25, 26, 27, 35, 36, 43],
             [1, 2, 4, 9, 11, 14, 18, 21, 26, 34, 36, 43], 4),
        91: ([1, 16, 74],
             [0, 3, 4, 5, 8, 11, 19, 25, 27, 43, 45, 50, 55],
             [0, 1, 4, 5, 13, 14, 15, 25, 28, 33, 38, 43, 44, 49, 55], 9),
        93: ([1, 4, 16, 64, 70],
             [3, 10, 11, 14, 21, 23, 33, 34, 46],
             [3, 9, 11, 17, 23, 33, 34, 46, 62], 25),
        129: ([1, 4, 16, 64, 97, 121, 127],
              [1, 9, 10, 14, 19, 21, 23, 26, 27],
              [2, 5, 9, 10, 13, 18, 22, 27, 43, 86], 13),
        397: ([1, 16, 31, 99, 126, 167, 256, 273, 290, 333, 393],
              [3, 5, 9, 10, 11, 12, 18, 20, 21, 23, 29, 33, 36, 40, 44, 47, 61, 72],
              [2, 3, 6, 10, 17, 22, 24, 33, 34, 36, 40, 46, 47, 53, 58, 71, 72],
              34)
    }

    exist = n in full_data or n in compact_data or \
        skew_spin_goethals_seidel_difference_family(n, existence=True)
    if existence:
        return exist

    if not exist:
        raise NotImplementedError(f'Data for spin type Goethals Seidel family of order {n} not yet implemented')

    G = Zmod(n)
    if skew_spin_goethals_seidel_difference_family(n, existence=True):
        G, [S1, S2, S3, S4] = skew_spin_goethals_seidel_difference_family(n, check=False)
    elif n in full_data:
        S1, S2, mu = full_data[n]
        S1 = list(map(G, S1))
        S2 = list(map(G, S2))
        S1, S2, S3, S4 = _construct_gs_difference_family_from_full(S1, S2, mu)
    elif n in compact_data:
        H, rep1, rep2, mu = compact_data[n]
        H = list(map(G, H))
        S1, S2, S3, S4 = _construct_gs_difference_family_from_compact(rep1, rep2, H, mu)

    if check:
        lmbda = len(S1) + len(S2) + len(S3) + len(S4) - n
        assert is_supplementary_difference_set([S1, S2, S3, S4], lmbda=lmbda, G=G)
    return G, [S1, S2, S3, S4]


def skew_spin_goethals_seidel_difference_family(n, existence=False, check=True):
    r"""
    Construct skew spin type Goethals-Seidel difference family with parameters
    `(n; k_1, k_2, k_3, k_4; \lambda)`.

    The construction is described in [Djo2024]_. This function contains, for
    each value of `n`, either a full representation of `S_1, S_2` together with
    the multiplier `\mu`, or a subgroup `H`, two sets of representatives, and the
    multiplier.

    This data is used to construct the difference family using the functions
    :func:`_construct_gs_difference_family_from_full` and
    :func:`_construct_gs_difference_family_from_compact`.

    INPUT:

    - ``n`` -- integer; the parameter of the GS difference family
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check
      whether the skew difference family can be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the sets
      are a skew difference family before returning them;
      setting this parameter to ``False`` may speed up the computation considerably

    OUTPUT:

    If ``existence=False``, the function returns the group G of integers modulo
    ``n`` and a list containing 4 sets, or raises an error if data for the given
    ``n`` is not available.
    If ``existence=True``, the function returns a boolean representing whether
    the skew difference family can be constructed.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import skew_spin_goethals_seidel_difference_family
        sage: G, [S1, S2, S3, S4] = skew_spin_goethals_seidel_difference_family(61)

    If existence is ``True``, the function returns a boolean::

        sage: skew_spin_goethals_seidel_difference_family(61, existence=True)
        True
        sage: skew_spin_goethals_seidel_difference_family(5, existence=True)
        False

    TESTS::
        sage: from sage.combinat.designs.difference_family import is_supplementary_difference_set, _is_skew_set
        sage: G, [S1, S2, S3, S4] = skew_spin_goethals_seidel_difference_family(7, check=False)
        sage: lmbda = len(S1) + len(S2) + len(S3) + len(S4) - 7
        sage: is_supplementary_difference_set([S1, S2, S3, S4], lmbda=lmbda, G=G)
        True
        sage: _is_skew_set(G, S1)
        True
        sage: skew_spin_goethals_seidel_difference_family(5)
        Traceback (most recent call last):
        ...
        NotImplementedError: Data for skew spin type Goethals Seidel family of order 5 not yet implemented
    """
    full_data = {
        7: ([1, 2, 4], [1, 6], 2),
        19: ([1, 4, 5, 6, 7, 9, 11, 16, 17], [0, 1, 7, 8, 11, 12, 18], -2),
        37: ([2, 3, 4, 6, 8, 11, 15, 18, 20, 21, 23, 24, 25, 27, 28, 30, 32, 36],
             [0, 1, 2, 5, 9, 13, 14, 15, 22, 23, 24, 28, 32, 35, 36],
             10)
    }
    compact_data = {
        61: ([1, 9, 20, 34, 58],  [3, 4, 5, 6, 8, 10],  [0, 8, 10, 13, 23, 26], 13),
        127: ([1, 2, 4, 8, 16, 32, 64],
              [1, 3, 7, 9, 11, 19, 21, 23, 47],
              [0, 3, 7, 9, 11, 15, 29, 31, 55], 19),
        271: ([1, 28, 106, 125, 169, 178, 242, 248, 258],
              [1, 4, 5, 7, 8, 11, 14, 16, 19, 21, 22, 25, 31, 43, 44],
              [1, 2, 3, 5, 7, 8, 12, 19, 22, 27, 38, 42, 44, 51], 5),
        331: ([1, 74, 80, 85, 111, 120, 167, 180, 270, 274, 293],
              [5, 10, 11, 13, 16, 19, 20, 22, 32, 38, 53, 56, 64, 76, 101],
              [0, 4, 11, 16, 20, 28, 31, 37, 41, 49, 53, 56, 73, 88, 101], 31),
        397: ([1, 16, 31, 99, 126, 167, 256, 273, 290, 333, 393],
              [1, 6, 7, 8, 9, 10, 11, 12, 17, 18, 20, 21, 29, 34, 46, 47, 53, 106],
              [2, 11, 12, 17, 18, 20, 24, 27, 33, 34, 36, 40, 46, 47, 53, 58, 71],
              34),
        547: ([1, 46, 237, 261, 293, 350, 353, 375, 440, 475, 509, 517, 519],
              [1, 4, 5, 6, 10, 11, 13, 14, 17, 25, 29, 34, 35, 40, 49, 52, 55, 64, 69, 110, 123],
              [1, 4, 5, 11, 16, 17, 20, 26, 32, 33, 34, 41, 49, 52, 55, 64, 70, 80, 123, 207],
              40),
        631: ([1, 8, 43, 64, 79, 188, 228, 242, 279, 310, 339, 344, 512, 562, 587],
              [1, 2, 3, 4, 6, 7, 12, 13, 14, 17, 19, 21, 26, 27, 31, 38, 42, 52, 62, 76, 124],
              [0, 11, 13, 14, 18, 19, 21, 22, 29, 35, 39, 46, 62, 63, 65, 66, 67, 92, 117, 124, 187],
              2)
    }

    if existence:
        return n in full_data or n in compact_data

    if n not in full_data and n not in compact_data:
        raise NotImplementedError(f'Data for skew spin type Goethals Seidel family of order {n} not yet implemented')

    G = Zmod(n)
    if n in full_data:
        S1, S2, mu = full_data[n]
        S1 = list(map(G, S1))
        S2 = list(map(G, S2))
        S1, S2, S3, S4 = _construct_gs_difference_family_from_full(S1, S2, mu)
    if n in compact_data:
        H, rep1, rep2, mu = compact_data[n]
        H = list(map(G, H))
        S1, S2, S3, S4 = _construct_gs_difference_family_from_compact(rep1, rep2, H, mu)

    if check:
        lmbda = len(S1) + len(S2) + len(S3) + len(S4) - n
        assert is_supplementary_difference_set([S1, S2, S3, S4], lmbda=lmbda, G=G)
        assert _is_skew_set(G, S1)
    return G, [S1, S2, S3, S4]


def skew_supplementary_difference_set(n, existence=False, check=True, return_group=False):
    r"""
    Construct `4-\{n; n_1, n_2, n_3, n_4; \lambda\}` supplementary difference sets,
    where `S_1` is skew and `n_1 + n_2 + n_3 + n_4 = n+\lambda`.

    These sets are constructed from available data, as described in [Djo1994a]_. The set `S_1 \subset G` is
    always skew, i.e. `S_1 \cap (-S_1) = \emptyset` and `S_1 \cup (-S_1) = G \setminus \{0\}`.

    The data is taken from:

    * `n = 103, 151`: [Djo1994a]_
    * `n = 67, 113, 127, 157, 163, 181, 241`: [Djo1992a]_
    * `n = 37, 43`: [Djo1992b]_
    * `n = 39, 49, 65, 93, 121, 129, 133, 217, 219, 267`: [Djo1992c]_
    * `n = 97`: [Djo2008a]_
    * `n = 109, 145, 247`: [Djo2008b]_
    * `n = 73`: [Djo2023b]_
    * `n = 213, 631`: [DGK2014]_
    * `n = 331`: [DK2016]_

    Additional skew Supplementary difference sets are built using the function
    :func:`skew_supplementary_difference_set_over_polynomial_ring`, and
    :func:`skew_supplementary_difference_set_with_paley_todd`.

    INPUT:

    - ``n`` -- integer; the parameter of the supplementary difference set
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check
      whether the supplementary difference sets can be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the sets
      are supplementary difference sets with `S_1` skew before returning them;
      setting this parameter to ``False`` may speed up the computation considerably
    - ``return_group`` -- boolean (default: ``False``); if ``True``, the function
      will also return the group from which the sets are created

    OUTPUT:

    If ``existence=False``, the function returns a list containing 4 sets,
    or raises an error if data for the given ``n`` is not available. If
    ``return_group=True`` the function will additionally return the group from
    which the sets are created.
    If ``existence=True``, the function returns a boolean representing whether
    skew supplementary difference sets can be constructed.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import skew_supplementary_difference_set
        sage: [S1, S2, S3, S4] = skew_supplementary_difference_set(39)

    If ``return_group=True``, the function will also return the group::

        sage: G, [S1, S2, S3, S4] = skew_supplementary_difference_set(103, return_group=True)

    If ``existence=True``, the function returns a boolean::

        sage: skew_supplementary_difference_set(103, existence=True)
        True
        sage: skew_supplementary_difference_set(17, existence=True)
        False

    TESTS::

        sage: from sage.combinat.designs.difference_family import is_supplementary_difference_set, _is_skew_set
        sage: G, [S1, S2, S3, S4] = skew_supplementary_difference_set(113, check=False, return_group=True)
        sage: is_supplementary_difference_set([S1, S2, S3, S4], lmbda=len(S1)+len(S2)+len(S3)+len(S4)-113, G=G)
        True
        sage: _is_skew_set(G, S1)
        True
        sage: G, [S1, S2, S3, S4] = skew_supplementary_difference_set(67, check=False, return_group=True)
        sage: is_supplementary_difference_set([S1, S2, S3, S4], lmbda=len(S1)+len(S2)+len(S3)+len(S4)-67, G=G)
        True
        sage: _is_skew_set(G, S1)
        True
        sage: skew_supplementary_difference_set(17)
        Traceback (most recent call last):
        ...
        ValueError: Skew SDS of order 17 not yet implemented.
        sage: skew_supplementary_difference_set(17, existence=True)
        False
        sage: skew_supplementary_difference_set(127, existence=True)
        True
        sage: skew_supplementary_difference_set(81, existence=True)
        True

    .. NOTE::

        The data for `n=247` in [Djo2008b]_ contains a typo: the set `\alpha_2` should contain `223` instead of `233`.
        This can be verified by checking the resulting sets, which are given explicitly in the paper.
    """

    # If -1 is present in an index set, it means that {0} should be added to that set
    indices = {
        37: [[0, 3, 5, 7, 9, 10], [0, 5, 6, 7, 8],
             [1, 2, 6, 7, 9], [2, 6, 8, 9, 10]],
        39: [[1, 3, 5, 6, 8, 10, 12], [0, 1, 5, 8, 12, 13],
             [1, 3, 4, 7, 9, 12, 13], [0, 1, 2, 3, 7, 8]],
        43: [[1, 2, 4], [1, 2, 4], [0, 2, 3], [3, 4, -1]],
        49: [[1, 2, 5, 7, 8, 10, 13, 14], [4, 5, 6, 7, 10, 11],
             [0, 1, 2, 4, 6, 7, 12, 14], [1, 2, 3, 5, 6, 10, 12, 13, 14]],
        65: [[1, 3, 5, 6, 8, 10, 13, 14, 17, 18, 20, 22],
             [0, 3, 7, 10, 16, 17, 18, 20, 21],
             [2, 4, 6, 8, 9, 10, 14, 15, 16, 17, 18, 20],
             [5, 7, 8, 9, 11, 12, 13, 14, 16, 18, 19, 20, 21]],
        67:  [[0, 3, 5, 6, 9, 10, 13, 14, 17, 18, 20],
              [0, 2, 4, 9, 11, 12, 13, 16, 19, 21],
              [1, 3, 6, 10, 11, 13, 14, 16, 20, 21],
              [2, 4, 6, 8, 9, 11, 14, 17, 19]],
        73: [[4, 6, 8, 14], [8, 10, 12, 14], [4, 6, 10, 12], [-1, 0, 2, 10]],
        93: [[0, 3, 4, 6, 9, 10, 12, 14, 17, 18],
             [2, 3, 4, 5, 9, 13, 15, 18, 19],
             [1, 2, 3, 4, 5, 6, 7, 8, 16],
             [1, 4, 6, 11, 12, 13, 15, 16, 17, 18]],
        97: [[1, 2, 4, 6, 9, 11, 13, 14, 17, 18, 21, 23, 25, 27, 29, 30],
             [1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 23, 27, 29],
             [0, 1, 2, 5, 6, 12, 13, 15, 16, 20, 24, 25, 26, 29, 30, 31],
             [0, 2, 3, 4, 7, 8, 9, 11, 12, 13, 15, 16, 17, 18, 23, 28, 29]],
        103: [[1, 3, 4, 6, 8, 11, 12, 14, 17, 18, 20, 22, 25, 27, 28, 30, 32],
              [2, 9, 10, 12, 13, 14, 15, 16, 20, 21, 22, 23, 24, 26, 28, 29, 30],
              [0, 1, 2, 3, 4, 11, 12, 13, 16, 17, 19, 20, 21, 24, 25, 26, 28, 30, 31],
              [0, 1, 2, 3, 4, 5, 6, 13, 15, 18, 19, 20, 23, 24, 25, 26, 27, 28, 29, 31]],
        109: [[0, 2, 5, 7, 8, 10, 12, 15, 16, 19, 20, 23, 24, 26, 29, 30, 33, 34],
              [4, 5, 6, 7, 11, 15, 18, 19, 20, 22, 25, 30, 32, 33, 35],
              [0, 1, 5, 6, 9, 10, 11, 14, 17, 20, 24, 26, 27, 28, 29, 31, 32],
              [0, 3, 4, 6, 7, 9, 10, 12, 13, 22, 24, 25, 26, 27, 28, 29, 31, 33, 35]],
        113: [[0, 3, 4, 6, 8, 10, 13, 14],
              [1, 3, 8, 9, 10, 11, 12, 13],
              [0, 2, 3, 5, 6, 7, 12],
              [1, 2, 3, 5, 8, 9, 15]],
        121: [[0, 2, 4, 7, 8, 11, 13, 14, 16, 19, 20, 22],
              [0, 1, 4, 5, 8, 9, 10, 15, 17, 20, 23],
              [1, 2, 3, 7, 9, 16, 18, 19, 20, 21, 22, 23],
              [0, 2, 9, 10, 11, 12, 13, 14, 15, 17, 18, 21, 22, 23]],
        127: [[0, 3, 5, 7, 8, 10, 12, 14, 16],
              [0, 1, 3, 6, 7, 9, 10, 12, 14, 15],
              [0, 1, 3, 4, 5, 7, 8, 9, 15, 16],
              [1, 4, 5, 6, 9, 10, 13, 14, 15, 16]],
        129: [[1, 2, 4, 7, 9, 11, 12, 14, 16, 18],
              [0, 1, 2, 3, 9, 11, 14, 15, 19],
              [0, 1, 3, 6, 8, 10, 12, 16, 18, 19],
              [0, 3, 7, 8, 9, 10, 12, 14, 15, 17]],
        133: [[1, 2, 5, 6, 9, 11, 12, 14], [1, 4, 7, 9, 10, 12, 13, 15],
              [0, 5, 6, 8, 11, 12, 13, 15], [0, 1, 2, 5, 7, 8, 9, 13, 14, 15]],
        145: [[1, 2, 4, 7, 9, 10, 13, 14, 16, 19, 20, 22], [0, 2, 4, 7, 10, 11, 14, 18, 19, 20, 21, 22],
              [1, 3, 6, 9, 12, 13, 14, 17, 19, 20, 21, 22, 23], [2, 3, 5, 6, 7, 9, 12, 13, 15, 16, 19, 20, 21, 22, 23]],
        151: [[0, 3, 5, 6, 8, 11, 13, 14, 16, 19, 21, 23, 25, 27, 28],
              [2, 3, 6, 13, 16, 17, 20, 23, 25, 26, 27, 28, 29],
              [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 23, 24, 27, 28],
              [1, 4, 5, 10, 11, 12, 13, 14, 16, 18, 19, 22, 25, 26, 27, 28]],
        157: [[0, 2, 5, 7, 8, 11],
              [0, 4, 5, 6, 9, 11],
              [6, 7, 8, 9, 10, 11],
              [0, 5, 6, 7, 8, 10, 11]],
        163: [[0, 2, 5, 6, 9, 10, 13, 14, 17],
              [0, 1, 7, 10, 12, 15, 16, 17],
              [0, 1, 3, 5, 8, 13, 15, 16, 17],
              [3, 6, 7, 8, 11, 12, 13, 14, 16, 17]],
        181: [[0, 3, 5, 6, 8, 10, 13, 15, 16, 19],
              [4, 5, 7, 8, 11, 14, 15, 16, 18, 19],
              [0, 4, 10, 11, 13, 15, 16, 18, 19],
              [2, 4, 5, 7, 11, 13, 15, 17, 19]],
        213: [[1, 2, 5, 6, 9, 11, 12, 14, 16, 19, 20, 23, 24, 26, 29, 30],
              [3, 6, 8, 12, 13, 14, 15, 17, 20, 22, 23, 25, 26, 27, 28, 31],
              [2, 3, 5, 7, 9, 13, 16, 17, 19, 21, 23, 24, 27, 28, 29],
              [0, 5, 6, 9, 11, 13, 14, 17, 20, 22, 23, 26, 29, 31]],
        217: [[0, 3, 5, 7, 8, 11, 12, 14], [1, 3, 4, 7, 9, 11, 12, 15],
              [3, 4, 5, 6, 7, 9, 10, 14, 15], [1, 3, 4, 5, 7, 8, 11, 13, 14]],
        219: [[1, 3, 5, 6, 8, 11, 12, 15, 17, 18, 21, 22, 24],
              [2, 6, 8, 10, 11, 12, 13, 16, 19, 22, 23, 24],
              [0, 1, 5, 6, 10, 11, 13, 14, 17, 20, 21, 24, 25],
              [0, 2, 3, 4, 5, 6, 7, 11, 12, 13, 16, 20, 23]],
        241: [[0, 2, 4, 6, 8, 11, 12, 14],
              [1, 3, 4, 6, 7, 13, 14, 15],
              [6, 8, 9, 10, 12, 13, 14, 15],
              [3, 4, 5, 9, 10, 13, 14]],
        247: [[0, 2, 4, 7, 8, 10, 12, 15, 16, 18, 20, 23, 25, 27, 29],
              [0, 2, 7, 9, 11, 12, 14, 15, 16, 18, 20, 22, 26],
              [2, 3, 4, 12, 13, 14, 15, 16, 18, 20, 23, 24, 26, 27, 29],
              [0, 3, 4, 6, 10, 11, 12, 14, 18, 19, 20, 22, 25, 29]],
        267: [[0, 3, 4, 7, 8, 11, 13, 15, 16, 19, 21, 22, 25],
              [0, 1, 4, 5, 6, 8, 14, 15, 18, 21, 23],
              [0, 2, 4, 5, 7, 9, 10, 11, 14, 15, 16, 17, 25],
              [0, 1, 3, 4, 6, 14, 15, 16, 17, 18, 20, 22, 23, 25]],
        331: [[1, 2, 4, 7, 9, 10, 12, 15, 16, 18, 21, 22, 24, 26, 28],
              [-1, 0, 2, 6, 9, 11, 12, 14, 15, 17, 20, 21, 24, 25, 28],
              [-1, 0, 1, 5, 6, 7, 8, 9, 10, 12, 15, 18, 23, 28, 29],
              [-1, 0, 3, 7, 8, 10, 11, 12, 14, 16, 19, 20, 21, 26, 29]],
        631: [[0, 2, 4, 6, 9, 10, 12, 15, 16, 18, 20, 23, 24, 26, 29, 30, 32, 35, 36, 38, 40],
              [0, 1, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 16, 20, 23, 28, 29, 30, 32, 36, 38, 41],
              [0, 2, 3, 4, 6, 9, 10, 12, 14, 15, 16, 17, 18, 19, 20, 22, 24, 25, 26, 29, 34, 40],
              [0, 2, 4, 5, 6, 7, 8, 10, 15, 16, 18, 22, 23, 24, 26, 30, 31, 33, 35, 36, 37, 38]],
    }

    # If the element is a list, that is the coset.
    cosets_gens = {
        37: [1, 2, 3, 5, 6, 9],
        39: [1, 2, 3, 4, 6, 8, [13]],
        43: [1, 3, 7],
        49: [1, 2, 3, 4, 6, 7, 9, 12],
        65: [1, 2, 3, 5, 6, 7, 9, 10, 11, [13], 22, [26]],
        67: [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 17],
        73: [1, 3, 5, 9, 11, 13, 17, 25],
        93: [1, 2, 3, 5, 7, 9, 10, 14, 15, [31]],
        97: [1, 2, 3, 4, 5, 6, 7, 9, 10, 12, 13, 15, 18, 20, 23, 26],
        103: [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 15, 17, 19, 21, 23, 30],
        109: [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 23, 25, 30],
        113: [1, 2, 3, 5, 6, 9, 10, 13],
        121: [1, 2, 4, 5, 7, 8, 10, 11, 16, 17, 19, 20],
        127: [1, 3, 5, 7, 9, 11, 13, 19, 21],
        129: [1, 3, 5, 7, 9, 11, 13, 19, 21, [43]],
        133: [1, 2, 3, 6, 7, 9, 18, [19, 38, 76]],
        145: [1, [2, 17, 32, 72, 77, 127, 137], [3, 43, 48, 98, 108, 118, 133], [6, 51, 71, 86, 91, 96, 121],
              [7, 52, 82, 107, 112, 117, 132], [11, 21, 31, 46, 61, 101, 106], [14, 19, 69, 79, 89, 104, 119],
              [22, 42, 57, 62, 67, 92, 122], [5, 35, 80, 100, 115, 120, 125], [10, 15, 55, 70, 85, 95, 105],
              [29], [58]],
        151: [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 15, 22, 27, 29, 30],
        157: [1, 2, 3, 5, 9, 15],
        163: [1, 2, 3, 5, 6, 9, 10, 15, 18],
        181: [1, 2, 3, 4, 6, 7, 8, 12, 13, 24],
        213: [1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 15, 20, 22, 30, 43, 71],
        217: [1, 2, 4, 5, 7, 10, 19, [31, 62, 124]],
        219: [1, 2, 3, 5, 7, 9, 11, 15, 19, 22, 23, 33, [73]],
        241: [1, 2, 4, 5, 7, 13, 19, 35],
        247: [1, [2, 18, 31, 32, 41, 110, 122, 162, 223], [3, 27, 48, 165, 170, 183, 185, 211, 243],
              [5, 28, 45, 58, 80, 158, 187, 201, 226], [6, 54, 83, 93, 96, 119, 123, 175, 239], [7, 20, 63, 73, 112, 138, 163, 180, 232],
              [10, 56, 69, 90, 116, 127, 155, 160, 205], [11, 47, 99, 102, 111, 115, 150, 176, 177], [13, 52, 65, 78, 91, 117, 143, 208, 221],
              [14, 29, 40, 79, 113, 126, 146, 217, 224], [17, 25, 43, 49, 140, 142, 153, 194, 225], [19, 57, 171],
              [33, 34, 37, 50, 59, 86, 98, 141, 203], [35, 66, 68, 74, 100, 118, 159, 172, 196], [38, 95, 114]],
        267: [1, 2, 3, 5, 7, 9, 10, 13, 14, 15, 19, 39, [89]],
        331: [1, 2, 4, 5, 7, 8, 10, 13, 14, 16, 19, 20, 28, 32, 56],
        631: [1, 2, 3, 4, 5, 6, 7, 9, 12, 14, 17, 18, 19, 21, 23, 27, 31, 35, 38, 42, 62],
    }

    H_db = {
        37: [1, 10, -11],
        39: [1, 16, 22],
        43: [1, 4, 11, 16, 21, -2, -8],
        49: [1, 18, 30],
        65: [1, 16, 61],
        67: [1, 29, 37],
        73: [1, 2, 4, 8, 16, 32, 64, 55, 37],
        93: [1, 4, 16, 64, 70],
        97: [1, 35, 61],
        103: [1, 46, 56],
        109: [1, 45, 63],
        113: [1, 16, 28, 30, 49, 106, 109],
        121: [1, 3, 9, 27, 81],
        127: [1, 2, 4, 8, 16, 32, 64],
        129: [1, 4, 16, 64, 97, 121, 127],
        133: [1, 4, 16, 25, 64, 93, 100, 106, 123],
        145: [1, 16, 36, 81, 111, 136, 141],
        151: [1, 8, 19, 59, 64],
        157: [1, 14, 16, 39, 46, 67, 75, 93, 99, 101, 108, 130, 153],
        163: [1, 38, 40, 53, 58, 85, 104, 133, 140],
        181: [1, 39, 43, 48, 62, 65, 73, 80, 132],
        213: [1, 37, 91, 103, 172, 187, 190],
        217: [1, 8, 9, 25, 51, 64, 72, 78, 81, 142, 190, 191, 193, 200, 214],
        219: [1, 4, 16, 37, 55, 64, 148, 154, 178],
        241: [1, 15, 24, 54, 87, 91, 94, 98, 100, 119, 160, 183, 205, 225, 231],
        247: [1, 9, 16, 55, 61, 81, 139, 144, 235],
        267: [1, 4, 16, 64, 67, 91, 97, 121, 217, 223, 256],
        331: [1, 74, 80, 85, 111, 120, 167, 180, 270, 274, 293],
        631: [1, 8, 43, 64, 79, 188, 228, 242, 279, 310, 339, 344, 512, 562, 587],
    }

    G, S1, S2, S3, S4 = None, [], [], [], []

    if n in indices:
        if existence:
            return True
        G, [S1, S2, S3, S4] = _construction_supplementary_difference_set(n, H_db[n], indices[n], cosets_gens[n], check=False)
    elif skew_supplementary_difference_set_over_polynomial_ring(n, existence=True):
        if existence:
            return True
        G, [S1, S2, S3, S4] = skew_supplementary_difference_set_over_polynomial_ring(n, check=False)
    elif skew_supplementary_difference_set_with_paley_todd(n, existence=True):
        if existence:
            return True
        G, [S1, S2, S3, S4] = skew_supplementary_difference_set_with_paley_todd(n, check=False)
    elif skew_spin_goethals_seidel_difference_family(n, existence=True):
        if existence:
            return True
        G, [S1, S2, S3, S4] = skew_spin_goethals_seidel_difference_family(n, check=False)

    if existence:
        return False

    if G is None:
        raise ValueError(f'Skew SDS of order {n} not yet implemented.')

    if check:
        lmbda = len(S1) + len(S2) + len(S3) + len(S4) - n
        assert is_supplementary_difference_set([S1, S2, S3, S4], lmbda=lmbda, G=G)
        assert _is_skew_set(G, S1)

    if return_group:
        return G, [S1, S2, S3, S4]
    return [S1, S2, S3, S4]


def _construction_supplementary_difference_set(n, H, indices, cosets_gen, check=True):
    r"""
    Construct `4-\{n; n_1, n_2, n_3, n_4; \lambda\}` supplementary difference sets,
    where `n_1 + n_2 + n_3 + n_4 = n+\lambda`.

    This construction is described in [Djo1994a]_.

    Let H be a subgroup of Zmod(n) of order `t`. We construct the `2s = n/t` cosets
    `\alpha_0, .., \alpha_{2s-1}` by using the elements contained in a sequence `A`:
    `\alpha_{2i} = a_iH` and `\alpha_{2i+1} = -\alpha_{2i}`.

    Then, we use four indices sets `J_1, J_2, J_3, J_4` to construct the four
    supplementary difference sets: `S_i = \bigcup_{j\in J__i} \alpha_i`. Note that
    if `J_i` contains the value `-1`, this function will add `0` to the set `S_i`.

    To construct a coset `\alpha_{2i}` by listing it directly, replace the `2i`-th
    element of the list `A` with the desired set.

    INPUT:

    - ``n`` -- integer; the parameter of the supplementary difference set
    - ``H`` -- list of integers; the set `H` used to construct the cosets
    - ``indices`` -- list containing four lists of integers; the sets
      `J_1, J_2, J_3, J_4` described above
    - ``cosets_gen`` -- list containing integers or list of integers; the set `A`
      described above
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the sets
      are supplementary difference sets; setting this parameter to ``False`` may
      speed up the computation considerably

    OUTPUT:

    The function returns the ring of integers modulo ``n`` and a list containing
    the 4 sets.

    TESTS::

        sage: from sage.combinat.designs.difference_family import is_supplementary_difference_set, _construction_supplementary_difference_set
        sage: H = [1, 10, -11]
        sage: cosets_gen = [1, 2, 3, 5, 6, 9]
        sage: indices =  [[0, 3, 5, 7, 9, 10], [0, 5, 6, 7, 8], [1, 2, 6, 7, 9], [2, 6, 8, 9, 10]]
        sage: _construction_supplementary_difference_set(37, H, indices, cosets_gen)
        (Ring of integers modulo 37,
        [[32, 1, 33, 35, 34, 7, 9, 10, 12, 14, 16, 17, 18, 22, 24, 26, 29, 31],
         [32, 1, 33, 34, 5, 6, 7, 8, 10, 13, 18, 19, 23, 24, 26],
         [32, 2, 36, 5, 11, 13, 14, 15, 18, 19, 20, 24, 27, 29, 31],
         [2, 5, 6, 8, 9, 12, 13, 14, 15, 16, 19, 20, 23, 29, 31]])
        sage: H = [1, 16, 22]
        sage: cosets_gen = [1, 2, 3, 4, 6, 8, [13]]
        sage: indices = [[1, 3, 5, 6, 8, 10, 12], [0, 1, 5, 8, 12, 13], [1, 3, 4, 7, 9, 12, 13], [0, 1, 2, 3, 7, 8]]
        sage: _construction_supplementary_difference_set(39, H, indices, cosets_gen)
        (Ring of integers modulo 39,
        [[4, 6, 7, 8, 10, 11, 12, 13, 15, 17, 18, 20, 23, 25, 30, 34, 36, 37, 38],
         [1, 36, 38, 6, 12, 13, 15, 16, 17, 18, 22, 23, 26, 30],
         [3, 7, 9, 13, 14, 17, 21, 23, 24, 26, 27, 29, 33, 34, 35, 37, 38],
         [32, 1, 2, 34, 35, 5, 38, 37, 7, 6, 14, 15, 16, 17, 18, 22, 23, 29]])
        sage: H = [1, 4, 11, 16, 21, -2, -8]
        sage: cosets_gen = [1, 3, 7]
        sage: indices = [[1, 2, 4], [1, 2, 4], [0, 2, 3], [3, 4, -1]]
        sage: G, sets = _construction_supplementary_difference_set(43, H, indices, cosets_gen, check=False)
        sage: is_supplementary_difference_set(sets, lmbda=35, G=G)
        True

    .. SEEALSO::

        :func:`skew_supplementary_difference_set`
    """
    def generate_set(index_set, cosets):
        S = set()
        for idx in index_set:
            if idx == -1:
                S.add(Z(0))
            else:
                S = S.union(cosets[idx])
        return list(S)

    Z = Zmod(n)
    H = set(map(Z, H))

    cosets = []
    for el in cosets_gen:
        if isinstance(el, list):
            even_coset = {Z(x) for x in el}
        else:
            even_coset = {x*el for x in H}
        odd_coset = {-x for x in even_coset}
        cosets.append(even_coset)
        cosets.append(odd_coset)

    S1 = generate_set(indices[0], cosets)
    S2 = generate_set(indices[1], cosets)
    S3 = generate_set(indices[2], cosets)
    S4 = generate_set(indices[3], cosets)

    if check:
        lmbda = len(S1) + len(S2) + len(S3) + len(S4) - n
        assert is_supplementary_difference_set([S1, S2, S3, S4], lmbda=lmbda, G=Z)

    return Z, [S1, S2, S3, S4]


def supplementary_difference_set_hadamard(n, existence=False, check=True):
    r"""
    Construct `4-\{n; n_1, n_2, n_3, n_4; \lambda\}` supplementary difference sets,
    where `n_1 + n_2 + n_3 + n_4 = n+\lambda`.

    These sets are constructed from available data, as described in [Djo1994a]_.
    The data is taken from:

    * `n = 191`: [Djo2008c]_
    * `n = 239`: [Djo1994b]_
    * `n = 251`: [DGK2014]_

    Additional SDS are constructed using :func:`skew_supplementary_difference_set`.

    INPUT:

    - ``n`` -- integer; the parameter of the supplementary difference set
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check
      whether the supplementary difference sets can be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the sets
      are supplementary difference sets before returning them; Setting this
      parameter to ``False`` may speed up the computation considerably

    OUTPUT:

    If ``existence=False``, the function returns the ring of integers modulo
    ``n`` and a list containing the 4 sets, or raises an error if data for the
    given ``n`` is not available.
    If ``existence=True``, the function returns a boolean representing whether
    skew supplementary difference sets can be constructed.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import supplementary_difference_set_hadamard
        sage: G, [S1, S2, S3, S4] = supplementary_difference_set_hadamard(191)

    If ``existence=True``, the function returns a boolean::

        sage: supplementary_difference_set_hadamard(191, existence=True)
        True
        sage: supplementary_difference_set_hadamard(17, existence=True)
        False

    TESTS::

        sage: from sage.combinat.designs.difference_family import is_supplementary_difference_set
        sage: G, [S1, S2, S3, S4] = supplementary_difference_set_hadamard(191, check=False)
        sage: is_supplementary_difference_set([S1, S2, S3, S4], lmbda=len(S1)+len(S2)+len(S3)+len(S4)-191, G=G)
        True
        sage: G, [S1, S2, S3, S4] = supplementary_difference_set_hadamard(37, check=False)
        sage: is_supplementary_difference_set([S1, S2, S3, S4], lmbda=len(S1)+len(S2)+len(S3)+len(S4)-37, G=G)
        True
        sage: supplementary_difference_set_hadamard(11)
        Traceback (most recent call last):
        ...
        ValueError: SDS of order 11 not yet implemented.
        sage: supplementary_difference_set_hadamard(11, existence=True)
        False
        sage: supplementary_difference_set_hadamard(127, existence=True)
        True
    """

    indices = {
        191: [[1, 7, 9, 10, 11, 13, 17, 18, 25, 26, 30, 31, 33, 34, 35, 36, 37],
              [1, 4, 7, 9, 11, 12, 13, 14, 19, 21, 22, 23, 24, 25, 26, 29, 36, 37],
              [0, 3, 4, 5, 7, 8, 9, 16, 17, 19, 24, 25, 29, 30, 31, 33, 35, 37],
              [1, 3, 4, 5, 8, 11, 14, 18, 19, 20, 21, 23, 24, 25, 28, 29, 30, 32, 34, 35]],
        239: [[0, 1, 2, 3, 4, 5, 6, 7, 14, 18, 19, 21, 24, 25, 29, 30],
              [0, 1, 3, 7, 9, 12, 15, 18, 20, 22, 26, 28, 29, 30, 31, 32, 33],
              [2, 3, 4, 5, 8, 9, 10, 11, 13, 17, 19, 21, 22, 24, 27, 31, 32],
              [0, 1, 2, 3, 6, 7, 8, 11, 13, 15, 17, 18, 19, 22, 25, 26, 27, 32, 33]],
        251: [[2, 6, 8, 10, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 27, 28, 35, 36, 39, 41, 43, 44, 47, 48],
              [2, 5, 10, 11, 17, 18, 21, 23, 24, 25, 26, 28, 29, 30, 34, 35, 38, 39, 40, 41, 42, 43, 44, 49],
              [0, 2, 6, 7, 10, 11, 14, 15, 16, 18, 21, 22, 24, 26, 30, 35, 37, 38, 45, 46, 47, 48, 49],
              [1, 2, 3, 4, 8, 9, 12, 17, 21, 22, 27, 28, 29, 30, 33, 34, 39, 41, 42, 43, 46, 47, 48]],
    }

    cosets_gens = {
        191: [1, 2, 3, 4, 6, 8, 9, 11, 12, 13, 16, 17, 18, 19, 22, 32, 36, 38, 41],
        239: [1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 18, 21, 28, 35, 42],
        251: [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 14, 15, 17, 18, 19, 21, 28, 30, 33, 34, 35, 41, 43, 45, 68],
    }

    H_db = {
        191: [1, 39, 184, 109, 49],
        239: [1, 10, 24, 44, 98, 100, 201],
        251: [1, 20, 113, 149, 219],
    }

    if existence:
        return n in indices or skew_supplementary_difference_set(n, existence=True)

    sets = None
    G = None
    if n in indices:
        G, sets = _construction_supplementary_difference_set(n, H_db[n], indices[n], cosets_gens[n], check=False)
    elif skew_supplementary_difference_set(n, existence=True):
        G, sets = skew_supplementary_difference_set(n, check=False, return_group=True)
    elif spin_goethals_seidel_difference_family(n, existence=True):
        if existence:
            return True
        G, [S1, S2, S3, S4] = spin_goethals_seidel_difference_family(n, check=False)

    if sets is None:
        raise ValueError(f'SDS of order {n} not yet implemented.')

    S1, S2, S3, S4 = sets
    if check:
        lmbda = len(S1) + len(S2) + len(S3) + len(S4) - n
        assert is_supplementary_difference_set([S1, S2, S3, S4], lmbda=lmbda, G=G)

    return G, [S1, S2, S3, S4]


def _is_skew_set(G, S):
    r"""
    Check if ``S`` is a skew set over the group ``G``.

    From [Djo1994a]_, a set `S \subset G` (where `G` is a finite abelian group of order `n`) is of skew
    type if `S_1 \cap (-S_1) = \emptyset` and `S_1 \cup (-S_1) = G\setminus \{0\}`.

    INPUT:

    - ``G`` -- a group
    - ``S`` -- list containing elements of ``G``; the set to be checked

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import _is_skew_set
        sage: Z5 = Zmod(5)
        sage: _is_skew_set(Z5, [Z5(1), Z5(2)])
        True
        sage: _is_skew_set(Z5, [Z5(1), Z5(2), Z5(3)])
        False
        sage: _is_skew_set(Z5, [Z5(1)])
        False
    """
    for el in S:
        if -el in S:
            return False
    for el in G:
        if el == 0:
            continue
        if el not in S and -el not in S:
            return False
    return True


def are_complementary_difference_sets(G, A, B, verbose=False):
    r"""
    Check if ``A`` and ``B`` are complementary difference sets over the group ``G``.

    According to [Sze1971]_, two sets `A`, `B` of size `m` are complementary
    difference sets over a group `G` of size `2m+1` if:

    1. they are `2-\{2m+1; m, m; m-1\}` supplementary difference sets
    2. `A` is skew, i.e. `a \in A` implies `-a \not \in A`

    INPUT:

    - ``G`` -- a group of odd order
    - ``A`` -- set of elements of ``G``
    - ``B`` -- set of elements of ``G``
    - ``verbose`` -- boolean (default: ``False``); if ``True`` the function will
      be verbose when the sets do not satisfy the constraints

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import are_complementary_difference_sets
        sage: are_complementary_difference_sets(Zmod(7), [1, 2, 4], [1, 2, 4])
        True

    If ``verbose=True``, the function will be verbose::

        sage: are_complementary_difference_sets(Zmod(7), [1, 2, 5], [1, 2, 4], verbose=True)
        The sets are not supplementary difference sets with lambda = 2
        False

    TESTS::

        sage: are_complementary_difference_sets(Zmod(16), [1, 2, 4], [1, 2, 4])
        False
        sage: are_complementary_difference_sets(Zmod(7), [1, 2, 4], [1, 2, 3, 4])
        False
        sage: are_complementary_difference_sets(Zmod(19), [1, 4, 5, 6, 7, 9, 11, 16, 17], [1, 4, 5, 6, 7, 9, 11, 16, 17])
        True

    .. SEEALSO::

        :func:`is_supplementary_difference_set`
    """
    n = G.order()

    if n % 2 != 1:
        if verbose:
            print('G must have odd order')
        return False

    m = (n - 1) // 2

    if len(A) != m or len(B) != m:
        if verbose:
            print(f'A and B must have size {m}')
        return False

    if not is_supplementary_difference_set([A, B], lmbda=m-1, G=G):
        if verbose:
            print(f'The sets are not supplementary difference sets with lambda = {m-1}')
        return False

    if not _is_skew_set(G, A):
        if verbose:
            print('The set A is not skew')
        return False

    return True


def complementary_difference_setsI(n, check=True):
    r"""
    Construct complementary difference sets in a group of order `n \cong 3 \mod 4`, `n` a prime power.

    Let `G` be a Galois Field of order `n`, where `n` satisfies the requirements
    above. Let `A` be the set of nonzero quadratic elements in `G`, and `B = A`.
    Then `A` and `B` are complementary difference sets over a group of order `n`.
    This construction is described in [Sze1971]_.

    INPUT:

    - ``n`` -- integer; the order of the group `G`
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the sets
      are complementary difference sets before returning them

    OUTPUT:

    The function returns the Galois field of order ``n`` and the two sets, or raises
    an error if ``n`` does not satisfy the requirements of this construction.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import complementary_difference_setsI
        sage: complementary_difference_setsI(19)
        (Finite Field of size 19,
         [1, 4, 5, 6, 7, 9, 11, 16, 17],
         [1, 4, 5, 6, 7, 9, 11, 16, 17])

    TESTS::

        sage: from sage.combinat.designs.difference_family import are_complementary_difference_sets
        sage: G, A, B = complementary_difference_setsI(23, check=False)
        sage: are_complementary_difference_sets(G, A, B)
        True
        sage: complementary_difference_setsI(17)
        Traceback (most recent call last):
        ...
        ValueError: the parameter 17 is not valid
        sage: complementary_difference_setsI(15)                                        # needs sage.libs.pari
        Traceback (most recent call last):
        ...
        ValueError: the parameter 15 is not valid

    .. SEEALSO::

        :func:`are_complementary_difference_sets`
        :func:`complementary_difference_sets`
    """
    if not (n % 4 == 3 and is_prime_power(n)):
        raise ValueError(f'the parameter {n} is not valid')

    from sage.rings.finite_rings.finite_field_constructor import GF

    G = GF(n, 'a')
    A = list({x**2 for x in G if x**2 != 0})
    B = A
    if check:
        assert are_complementary_difference_sets(G, A, B)

    return G, A, B


def complementary_difference_setsII(n, check=True):
    r"""
    Construct complementary difference sets in a group of order `n = p^t`, where `p \cong 5 \mod 8` and `t \cong 1, 2, 3 \mod 4`.

    Consider a finite field `G` of order `n` and let `\rho` be the generator of
    the corresponding multiplicative group. Then, there are two different constructions,
    depending on whether `t` is even or odd.

    If `t \cong 2 \mod 4`, let `C_0` be the set of nonzero octic residues in `G`,
    and let `C_i = \rho^i C_0` for `1 \le i \le  7`.
    Then, `A = C_0 \cup C_1 \cup C_2 \cup C_3` and  `B = C_0 \cup C_1 \cup C_6 \cup C_7`.

    If `t` is odd, let `C_0` be the set of nonzero fourth powers in `G`, and let
    `C_i = \rho^i C_0` for `1 \le i \le  3`.
    Then, `A = C_0 \cup C_1` and  `B = C_0 \cup C_3`.

    For more details on this construction, see [Sze1971]_.

    INPUT:

    - ``n`` -- integer; the order of the group `G`
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the sets
      are complementary difference sets before returning them; setting this to
      ``False`` might speed up the computation for large values of ``n``

    OUTPUT:

    The function returns the Galois field of order ``n`` and the two sets, or raises
    an error if ``n`` does not satisfy the requirements of this construction.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import complementary_difference_setsII
        sage: complementary_difference_setsII(5)                                        # needs sage.libs.pari
        (Finite Field of size 5, [1, 2], [1, 3])

    TESTS::

        sage: # needs sage.libs.pari
        sage: from sage.combinat.designs.difference_family import are_complementary_difference_sets
        sage: G, A, B = complementary_difference_setsII(25, check=False)
        sage: are_complementary_difference_sets(G, A, B)
        True
        sage: G, A, B = complementary_difference_setsII(13, check=False)
        sage: are_complementary_difference_sets(G, A, B)
        True
        sage: complementary_difference_setsII(49)
        Traceback (most recent call last):
        ...
        ValueError: the parameter 49 is not valid
        sage: complementary_difference_setsII(15)
        Traceback (most recent call last):
        ...
        ValueError: the parameter 15 is not valid

    .. SEEALSO::

        :func:`are_complementary_difference_sets`
        :func:`complementary_difference_sets`
    """
    p, t = is_prime_power(n, get_data=True)
    if not (p % 8 == 5 and t > 0 and t % 4 in [1, 2, 3]):
        raise ValueError(f'the parameter {n} is not valid')

    from sage.rings.finite_rings.finite_field_constructor import GF
    G = GF(n, 'a')
    A, B = None, None

    if t % 2 == 0:
        rho = G.multiplicative_generator()
        C0 = list({el**8 for el in G if el != 0})
        C1, C2, C3, C6, C7 = ([rho**i * el for el in C0] for i in [1, 2, 3, 6, 7])
        A = C0 + C1 + C2 + C3
        B = C0 + C1 + C6 + C7
    else:
        rho = G.multiplicative_generator()
        C0 = list({el**4 for el in G if el**4 != 0})
        C1 = [rho * el for el in C0]
        C3 = [rho**3 * el for el in C0]
        A = C0 + C1
        B = C0 + C3

    if check:
        assert are_complementary_difference_sets(G, A, B)

    return G, A, B


def complementary_difference_setsIII(n, check=True):
    r"""
    Construct complementary difference sets in a group of order `n = 2m + 1`, where `4m + 3` is a prime power.

    Consider a finite field `G` of order `n` and let `\rho` be a primite element
    of this group. Now let `Q` be the set of nonzero quadratic residues in `G`,
    and let `A = \{ a | \rho^{2a} - 1 \in Q\}`, `B' = \{ b | -(\rho^{2b} + 1) \in Q\}`.
    Then `A` and `B = Q \setminus B'` are complementary difference sets over the ring
    of integers modulo `n`. For more details, see [Sz1969]_.

    INPUT:

    - ``n`` -- integer; the order of the group over which the sets are constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the sets
      are complementary difference sets before returning them; setting this to
      ``False`` might speed up the computation for large values of ``n``

    OUTPUT:

    The function returns the Galois field of order ``n`` and the two sets, or raises
    an error if ``n`` does not satisfy the requirements of this construction.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import complementary_difference_setsIII
        sage: complementary_difference_setsIII(11)                                      # needs sage.libs.pari
        (Ring of integers modulo 11, [1, 2, 5, 7, 8], [0, 1, 3, 8, 10])

    TESTS::

        sage: from sage.combinat.designs.difference_family import are_complementary_difference_sets
        sage: G, A, B = complementary_difference_setsIII(21, check=False)               # needs sage.libs.pari
        sage: are_complementary_difference_sets(G, A, B)                                # needs sage.libs.pari
        True
        sage: G, A, B = complementary_difference_setsIII(65, check=False)               # needs sage.libs.pari
        sage: are_complementary_difference_sets(G, A, B)                                # needs sage.libs.pari
        True
        sage: complementary_difference_setsIII(10)
        Traceback (most recent call last):
        ...
        ValueError: the parameter 10 is not valid
        sage: complementary_difference_setsIII(17)                                      # needs sage.libs.pari
        Traceback (most recent call last):
        ...
        ValueError: the parameter 17 is not valid

    .. SEEALSO::

        :func:`are_complementary_difference_sets`
        :func:`complementary_difference_sets`
    """
    m = (n - 1) // 2
    q = 4*m + 3
    if n % 2 != 1 or not is_prime_power(q):
        raise ValueError(f'the parameter {n} is not valid')

    from sage.rings.finite_rings.finite_field_constructor import GF
    G = Zmod(n)
    G2 = GF(q)
    rho = G2.primitive_element()

    Q = [rho ** (2*b) for b in range(1, n+1)]

    A = [G(a) for a in range(n) if rho**(2*a) - 1 in Q]
    B = [G(b) for b in range(n) if -rho**(2*b) - 1 not in Q]

    if check:
        assert are_complementary_difference_sets(G, A, B)

    return G, A, B


def complementary_difference_sets(n, existence=False, check=True):
    r"""
    Compute complementary difference sets over a group of order `n = 2m + 1`.

    According to [Sze1971]_, two sets `A`, `B` of size `m` are complementary
    difference sets over a group `G` of size `n = 2m + 1` if:

    1. they are `2-\{2m+1; m, m; m-1\}` supplementary difference sets
    2. `A` is skew, i.e. `a \in A` implies `-a \not \in A`

    This method tries to call :func:`complementary_difference_setsI`,
    :func:`complementary_difference_setsII` or :func:`complementary_difference_setsIII`
    if the parameter `n` satisfies the requirements of one of these functions.

    INPUT:

    - ``n`` -- integer; the order of the group over which the sets are constructed
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check
      whether the supplementary difference sets can be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the sets
      are complementary difference sets before returning them; setting this to
      ``False`` might speed up the computation for large values of ``n``

    OUTPUT:

    If ``existence=False``, the function returns group ``G`` and two complementary
    difference sets, or raises an error if data for the given ``n`` is not available.
    If ``existence=True``, the function returns a boolean representing whether
    complementary difference sets can be constructed for the given ``n``.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import complementary_difference_sets
        sage: complementary_difference_sets(15)                                         # needs sage.libs.pari
        (Ring of integers modulo 15, [1, 2, 4, 6, 7, 10, 12], [0, 1, 2, 6, 9, 13, 14])

    If ``existence=True``, the function returns a boolean::

        sage: complementary_difference_sets(15, existence=True)                         # needs sage.libs.pari
        True
        sage: complementary_difference_sets(16, existence=True)
        False

    TESTS::

        sage: from sage.combinat.designs.difference_family import are_complementary_difference_sets
        sage: G, A, B = complementary_difference_sets(29)                               # needs sage.libs.pari
        sage: are_complementary_difference_sets(G, A, B)                                # needs sage.libs.pari
        True
        sage: G, A, B = complementary_difference_sets(65)                               # needs sage.libs.pari
        sage: are_complementary_difference_sets(G, A, B)                                # needs sage.libs.pari
        True
        sage: complementary_difference_sets(10)
        Traceback (most recent call last):
        ...
        ValueError: the parameter n must be odd
        sage: complementary_difference_sets(17)                                         # needs sage.libs.pari
        Traceback (most recent call last):
        ...
        NotImplementedError: complementary difference sets of order 17 are not implemented yet

    .. SEEALSO::

        :func:`are_complementary_difference_sets`
    """
    if n % 2 == 0:
        if existence:
            return False
        raise ValueError('the parameter n must be odd')

    p, t = is_prime_power(n, get_data=True)
    G, A, B = None, None, None

    if n % 4 == 3 and t > 0:
        if existence:
            return True
        G, A, B = complementary_difference_setsI(n, check=False)
    elif p % 8 == 5 and t > 0 and t % 4 in [1, 2, 3]:
        if existence:
            return True
        G, A, B = complementary_difference_setsII(n, check=False)
    elif is_prime_power(2*n + 1):
        if existence:
            return True
        G, A, B = complementary_difference_setsIII(n, check=False)

    if existence:
        return False

    if G is None:
        raise NotImplementedError(f'complementary difference sets of order {n} are not implemented yet')

    if check:
        assert are_complementary_difference_sets(G, A, B)
    return G, A, B


def difference_family(v, k, l=1, existence=False, explain_construction=False, check=True):
    r"""
    Return a (``k``, ``l``)-difference family on an Abelian group of cardinality ``v``.

    Let `G` be a finite Abelian group. For a given subset `D` of `G`, we define
    `\Delta D` to be the multi-set of differences `\Delta D = \{x - y; x \in D,
    y \in D, x \not= y\}`. A `(G,k,\lambda)`-*difference family* is a collection
    of `k`-subsets of `G`, `D = \{D_1, D_2, \ldots, D_b\}` such that the union
    of the difference sets `\Delta D_i` for `i=1,...b`, seen as a multi-set,
    contains each element of `G \backslash \{0\}` exactly `\lambda`-times.

    When there is only one block, i.e. `\lambda(v - 1) = k(k-1)`, then a
    `(G,k,\lambda)`-difference family is also called a *difference set*.

    See also :wikipedia:`Difference_set`.

    If there is no such difference family, an :exc:`EmptySetError` is raised and
    if there is no construction at the moment :exc:`NotImplementedError`
    is raised.

    INPUT:

    - ``v``, ``k``, ``l`` -- parameters of the difference family. If ``l`` is
      not provided it is assumed to be ``1``
    - ``existence`` -- if ``True``, then return either ``True`` if Sage knows
      how to build such design, ``Unknown`` if it does not and ``False`` if it
      knows that the design does not exist
    - ``explain_construction`` -- instead of returning a difference family,
      returns a string that explains the construction used
    - ``check`` -- boolean (default: ``True``); if ``True``, then the result of
      the computation is checked before being returned. This should not be
      needed but ensures that the output is correct

    OUTPUT:

    A pair ``(G,D)`` made of a group `G` and a difference family `D` on that
    group. Or, if ``existence=True`` a troolean or if
    ``explain_construction=True`` a string.

    EXAMPLES::

        sage: G,D = designs.difference_family(73,4)                                     # needs sage.libs.pari
        sage: G                                                                         # needs sage.libs.pari
        Finite Field of size 73
        sage: D                                                                         # needs sage.libs.pari
        [[0, 1, 5, 18],
         [0, 3, 15, 54],
         [0, 9, 45, 16],
         [0, 27, 62, 48],
         [0, 8, 40, 71],
         [0, 24, 47, 67]]

        sage: print(designs.difference_family(73, 4, explain_construction=True))
        The database contains a (73,4)-evenly distributed set

        sage: # needs sage.libs.pari
        sage: G,D = designs.difference_family(15,7,3)
        sage: G
        Ring of integers modulo 15
        sage: D
        [[0, 1, 2, 4, 5, 8, 10]]
        sage: print(designs.difference_family(15,7,3,explain_construction=True))
        Singer difference set

        sage: # needs sage.libs.pari
        sage: print(designs.difference_family(91,10,1,explain_construction=True))
        Singer difference set
        sage: print(designs.difference_family(64,28,12, explain_construction=True))
        McFarland 1973 construction
        sage: print(designs.difference_family(576, 276, 132, explain_construction=True))
        Hadamard difference set product from N1=2 and N2=3

    For `k=6,7` we look at the set of small prime powers for which a
    construction is available::

        sage: def prime_power_mod(r, m):
        ....:     k = m+r
        ....:     while True:
        ....:         if is_prime_power(k):
        ....:             yield k
        ....:         k += m

        sage: # needs sage.libs.pari
        sage: from itertools import islice
        sage: l6 = {True: [], False: [], Unknown: []}
        sage: for q in islice(prime_power_mod(1,30), int(60)):
        ....:     l6[designs.difference_family(q,6,existence=True)].append(q)
        sage: l6[True]
        [31, 121, 151, 181, 211, ...,  3061, 3121, 3181]
        sage: l6[Unknown]
        [61]
        sage: l6[False]
        []

        sage: # needs sage.libs.pari
        sage: l7 = {True: [], False: [], Unknown: []}
        sage: for q in islice(prime_power_mod(1,42), int(60)):
        ....:     l7[designs.difference_family(q,7,existence=True)].append(q)
        sage: l7[True]
        [169, 337, 379, 421, 463, 547, 631, 673, 757, 841, 883, 967, ...,  4621, 4957, 5167]
        sage: l7[Unknown]
        [43, 127, 211, 2017, 2143, 2269, 2311, 2437, 2521, 2647, ..., 4999, 5041, 5209]
        sage: l7[False]
        []

    List available constructions::

        sage: for v in range(2,100):                                                    # needs sage.libs.pari
        ....:     constructions = []
        ....:     for k in range(2,10):
        ....:         for l in range(1,10):
        ....:             ret = designs.difference_family(v,k,l,existence=True)
        ....:             if ret is True:
        ....:                 constructions.append((k,l))
        ....:                 _ = designs.difference_family(v,k,l)
        ....:     if constructions:
        ....:         print("%2d: %s"%(v, ', '.join('(%d,%d)'%(k,l) for k,l in constructions)))
         3: (2,1)
         4: (3,2)
         5: (2,1), (4,3)
         6: (5,4)
         7: (2,1), (3,1), (3,2), (4,2), (6,5)
         8: (7,6)
         9: (2,1), (4,3), (8,7)
        10: (9,8)
        11: (2,1), (4,6), (5,2), (5,4), (6,3)
        13: (2,1), (3,1), (3,2), (4,1), (4,3), (5,5), (6,5)
        15: (3,1), (4,6), (5,6), (7,3), (7,6)
        16: (3,2), (5,4), (6,2)
        17: (2,1), (4,3), (5,5), (8,7)
        19: (2,1), (3,1), (3,2), (4,2), (6,5), (9,4), (9,8)
        21: (3,1), (4,3), (5,1), (6,3), (6,5)
        22: (4,2), (6,5), (7,4), (8,8)
        23: (2,1)
        25: (2,1), (3,1), (3,2), (4,1), (4,3), (6,5), (7,7), (8,7)
        27: (2,1), (3,1)
        28: (3,2), (6,5)
        29: (2,1), (4,3), (7,3), (7,6), (8,4), (8,6)
        31: (2,1), (3,1), (3,2), (4,2), (5,2), (5,4), (6,1), (6,5)
        33: (3,1), (5,5), (6,5)
        34: (4,2)
        35: (5,2)
        37: (2,1), (3,1), (3,2), (4,1), (4,3), (6,5), (9,2), (9,8)
        39: (3,1), (6,5)
        40: (3,2), (4,1)
        41: (2,1), (4,3), (5,1), (5,4), (6,3), (8,7)
        43: (2,1), (3,1), (3,2), (4,2), (6,5), (7,2), (7,3), (7,6), (8,4)
        45: (3,1), (5,1)
        46: (4,2), (6,2)
        47: (2,1)
        49: (2,1), (3,1), (3,2), (4,1), (4,3), (6,5), (8,7), (9,3)
        51: (3,1), (5,2), (6,3)
        52: (4,1)
        53: (2,1), (4,3)
        55: (3,1), (9,4)
        57: (3,1), (7,3), (8,1)
        59: (2,1)
        61: (2,1), (3,1), (3,2), (4,1), (4,3), (5,1), (5,4), (6,2), (6,3), (6,5)
        63: (3,1)
        64: (3,2), (4,1), (7,2), (7,6), (9,8)
        65: (5,1)
        67: (2,1), (3,1), (3,2), (6,5)
        69: (3,1)
        71: (2,1), (5,2), (5,4), (7,3), (7,6), (8,4)
        73: (2,1), (3,1), (3,2), (4,1), (4,3), (6,5), (8,7), (9,1), (9,8)
        75: (3,1), (5,2)
        76: (4,1)
        79: (2,1), (3,1), (3,2), (6,5)
        81: (2,1), (3,1), (4,3), (5,1), (5,4), (8,7)
        83: (2,1)
        85: (4,1), (7,2), (7,3), (8,2)
        89: (2,1), (4,3), (8,7)
        91: (6,1), (7,1)
        97: (2,1), (3,1), (3,2), (4,1), (4,3), (6,5), (8,7), (9,3)

    TESTS:

    Check more of the Wilson constructions from [Wi72]_::

        sage: Q5 = [241, 281,421,601,641, 661, 701, 821,881]
        sage: Q9 = [73, 1153, 1873, 2017]
        sage: Q15 = [76231]
        sage: Q4 = [13, 73, 97, 109, 181, 229, 241, 277, 337, 409, 421, 457]
        sage: Q8 = [1009, 3137, 3697]
        sage: for Q,k in [(Q4,4),(Q5,5),(Q8,8),(Q9,9),(Q15,15)]:                        # needs sage.libs.pari
        ....:     for q in Q:
        ....:         assert designs.difference_family(q,k,1,existence=True) is True
        ....:         _ = designs.difference_family(q,k,1)

    Check Singer difference sets::

        sage: sgp = lambda q,d: ((q**(d+1)-1)//(q-1), (q**d-1)//(q-1), (q**(d-1)-1)//(q-1))

        sage: for q in range(2,10):                                                     # needs sage.libs.pari
        ....:     if is_prime_power(q):
        ....:         for d in [2,3,4]:
        ....:           v,k,l = sgp(q,d)
        ....:           assert designs.difference_family(v,k,l,existence=True) is True
        ....:           _ = designs.difference_family(v,k,l)

    Check twin primes difference sets::

        sage: for p in [3,5,7,9,11]:                                                    # needs sage.libs.pari
        ....:     v = p*(p+2); k = (v-1)/2;  lmbda = (k-1)/2
        ....:     G,D = designs.difference_family(v,k,lmbda)

    Check Complementary difference sets::

        sage: for v in [15, 33, 35, 39, 51]:                                            # needs sage.libs.pari
        ....:     G, D = designs.difference_family(v, (v-1)//2, (v-1)//2-1)

    Check the database::

        sage: from sage.combinat.designs.database import DF,EDS
        sage: for v,k,l in DF:
        ....:     assert designs.difference_family(v,k,l,existence=True) is True
        ....:     df = designs.difference_family(v,k,l,check=True)

        sage: for k in EDS:                                                             # needs sage.libs.pari
        ....:     for v in EDS[k]:
        ....:         assert designs.difference_family(v,k,1,existence=True) is True
        ....:         df = designs.difference_family(v,k,1,check=True)

    Check the known Hadamard parameters::

        sage: for N in range(2,21):                                                     # needs sage.libs.pari
        ....:     v = 4*N^2; k = 2*N^2-N; l = N^2-N
        ....:     status = designs.difference_family(v,k,l,existence=True)
        ....:     print("{:2} {}".format(N,designs.difference_family(v,k,l,explain_construction=True) if status is True else status))
        2 McFarland 1973 construction
        3 Turyn 1965 construction
        4 McFarland 1973 construction
        5 False
        6 The database contains a (144,66,30)-difference family
        7 False
        8 McFarland 1973 construction
        9 Unknown
        10 Unknown
        11 False
        12 Hadamard difference set product from N1=2 and N2=3
        13 False
        14 Unknown
        15 Unknown
        16 McFarland 1973 construction
        17 False
        18 Hadamard difference set product from N1=3 and N2=3
        19 False
        20 Unknown

    Check a failing construction (:issue:`17528`)::

        sage: designs.difference_family(9,3)
        Traceback (most recent call last):
        ...
        NotImplementedError: No construction available for (9,3,1)-difference family

    Check that when ``existence=True`` we always obtain ``True``, ``False`` or ``Unknown``,
    and when ``explain_construction=True``, it is a string (see :issue:`24513`)::

        sage: designs.difference_family(3, 2, 1, existence=True)
        True
        sage: designs.difference_family(3, 2, 1, explain_construction=True)
        'Trivial difference family'

        sage: for _ in range(100):                                                      # needs sage.libs.pari
        ....:     v = randint(1, 30)
        ....:     k = randint(2, 30)
        ....:     l = randint(1, 30)
        ....:     res = designs.difference_family(v, k, l, existence=True)
        ....:     assert res is True or res is False or res is Unknown
        ....:     if res is True:
        ....:         assert isinstance(designs.difference_family(3, 2, 1, explain_construction=True), str)

    .. TODO::

        Implement recursive constructions from Buratti "Recursive for difference
        matrices and relative difference families" (1998) and Jungnickel
        "Composition theorems for difference families and regular planes" (1978)
    """
    from .block_design import are_hyperplanes_in_projective_geometry_parameters

    from .database import DF, EDS

    v = ZZ(v)
    k = ZZ(k)
    l = ZZ(l)

    if v < 0 or k < 0 or l < 0:
        if existence:
            return False
        raise EmptySetError("No difference family eixsts with negative parameters")

    if (v,k,l) in DF:
        if existence:
            return True
        elif explain_construction:
            return "The database contains a ({},{},{})-difference family".format(v,k,l)

        vv, blocks = next(iter(DF[v,k,l].items()))

        # Build the group
        from sage.rings.finite_rings.integer_mod_ring import Zmod
        if len(vv) == 1:
            G = Zmod(vv[0])
        else:
            from sage.categories.cartesian_product import cartesian_product
            G = cartesian_product([Zmod(i) for i in vv])

        df = [[G(i) for i in b] for b in blocks]

        if check and not is_difference_family(G, df, v=v, k=k, l=l):
            raise RuntimeError("There is an invalid ({},{},{})-difference "
                    "family in the database... Please contact "
                    "sage-devel@googlegroups.com".format(v,k,l))

        return G,df

    elif l == 1 and k in EDS and v in EDS[k]:
        if existence:
            return True
        elif explain_construction:
            return "The database contains a ({},{})-evenly distributed set".format(v,k)

        from sage.rings.finite_rings.finite_field_constructor import GF
        poly,B = EDS[k][v]
        if poly is None:  # q is prime
            K = G = GF(v)
        else:
            K = G = GF(v,'a',modulus=poly)

        B = [K(b) for b in B]
        e = k*(k-1)//2
        xe = G.multiplicative_generator()**e
        df = [[xe**j*b for b in B] for j in range((v-1)//(2*e))]
        if check and not is_difference_family(G, df, v=v, k=k, l=l):
            raise RuntimeError("There is an invalid ({},{})-evenly distributed "
                               "set in the database... Please contact "
                               "sage-devel@googlegroups.com".format(v, k))
        return G, df

    if k in [0,1]:
        # Then \Delta D_i is empty
        # So if G\{0} is empty is good, otherwise not
        if v == 1:
            if existence:
                return True
            from sage.rings.finite_rings.integer_mod_ring import Zmod
            l = [0] if k == 1 else []
            return Zmod(1),[l]

        if existence:
            return False
        raise EmptySetError("No difference family exists with k=1 and v!=1")

    e = k*(k-1)
    if (l*(v-1)) % e:
        if existence:
            return Unknown
        raise NotImplementedError("No construction available for ({},{},{})-difference family".format(v,k,l))

    # trivial construction
    if k == (v-1) and l == (v-2):
        if existence:
            return True
        elif explain_construction:
            return "Trivial difference family"

        from sage.rings.finite_rings.integer_mod_ring import Zmod
        G = Zmod(v)
        return G, [list(range(1, v))]

    factorization = arith.factor(v)
    if len(factorization) == 1:
        from sage.rings.finite_rings.finite_field_constructor import GF
        K = GF(v,'z')

    if are_mcfarland_1973_parameters(v,k,l):
        if existence:
            return True
        elif explain_construction:
            return "McFarland 1973 construction"
        else:
            _, (q,s) = are_mcfarland_1973_parameters(v,k,l,True)
            G,D = mcfarland_1973_construction(q,s)

    elif are_hyperplanes_in_projective_geometry_parameters(v,k,l):
        if existence:
            return True
        elif explain_construction:
            return "Singer difference set"
        else:
            _, (q,d) = are_hyperplanes_in_projective_geometry_parameters(v,k,l,True)
            G,D = singer_difference_set(q,d)

    elif are_hadamard_difference_set_parameters(v,k,l) and k-2*l == 3:
        if existence:
            return True
        elif explain_construction:
            return "Turyn 1965 construction"
        else:
            G,D = turyn_1965_3x3xK(4)

    elif are_hadamard_difference_set_parameters(v,k,l) and hadamard_difference_set_product_parameters(k-2*l):
        N1,N2 = hadamard_difference_set_product_parameters(k-2*l)
        if existence:
            return True
        elif explain_construction:
            return "Hadamard difference set product from N1={} and N2={}".format(N1,N2)
        else:
            v1 = 4*N1*N1
            v2 = 4*N2*N2
            k1 = 2*N1*N1 - N1
            k2 = 2*N2*N2 - N2
            l1 = N1*N1 - N1
            l2 = N2*N2 - N2
            G1, D1 = difference_family(v1,k1,l1)
            G2, D2 = difference_family(v2,k2,l2)
            G, D = hadamard_difference_set_product(G1,D1,G2,D2)

    elif are_hadamard_difference_set_parameters(v,k,l) and (k-2*l).is_prime():
        if existence:
            return False
        else:
            raise EmptySetError("by McFarland 1989 such difference family does not exist")

    elif len(factorization) == 1 and radical_difference_family(K, k, l, existence=True) is True:
        if existence:
            return True
        elif explain_construction:
            return "Radical difference family on a finite field"
        else:
            D = radical_difference_family(K,k,l)
            G = K

    elif (len(factorization) == 1
        and l == 1
        and k == 6
        and df_q_6_1(K, existence=True) is True):
        if existence:
            return True
        elif explain_construction:
            return "Wilson 1972 difference family made from the union of two cyclotomic cosets"
        else:
            D = df_q_6_1(K)
            G = K

    elif (k == (v-1)//2 and
          l == (k-1)//2 and
          len(factorization) == 2 and
          abs(pow(*factorization[0]) - pow(*factorization[1])) == 2):
        # Twin prime powers construction
        # i.e. v = p(p+2) where p and p+2 are prime powers
        #      k = (v-1)/2
        #      lambda = (k-1)/2 (ie 2l+1 = k)
        if existence:
            return True
        elif explain_construction:
            return "Twin prime powers difference family"
        else:
            p = pow(*factorization[0])
            q = pow(*factorization[1])
            if p > q:
                p,q = q,p
            G,D = twin_prime_powers_difference_set(p,check=False)

    elif (v-1)//2 == k and (v-1)//2-1 == l and complementary_difference_sets(v, existence=True):
        if existence:
            return True
        elif explain_construction:
            return "Complementary difference sets"
        else:
            G, A, B = complementary_difference_sets(v)
            D = [A, B]

    else:
        if existence:
            return Unknown
        raise NotImplementedError("No constructions for these parameters")

    if check and not is_difference_family(G,D,v=v,k=k,l=l,verbose=False):
        raise RuntimeError("There is a problem. Sage built the following "
                "difference family on G='{}' with parameters ({},{},{}):\n "
                "{}\nwhich seems to not be a difference family... "
                "Please contact sage-devel@googlegroups.com".format(G,v,k,l,D))

    return G, D


from sage.misc.rest_index_of_methods import gen_rest_table_index
import sys
__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=gen_rest_table_index(sys.modules[__name__]))
