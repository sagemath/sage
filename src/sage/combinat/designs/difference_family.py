# -*- coding: utf-8 -*-
r"""
Difference families

This module gathers everything related to difference families. One can build a
difference family (or check that it can be built) with :func:`difference_family`::

    sage: G,F = designs.difference_family(13,4,1)

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

from sage.misc.cachefunc import cached_function

from sage.categories.sets_cat import EmptySetError
import sage.arith.all as arith
from sage.misc.unknown import Unknown
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ


def group_law(G):
    r"""
    Return a triple ``(identity, operation, inverse)`` that define the
    operations on the group ``G``.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import group_law
        sage: group_law(Zmod(3))
        (0, <built-in function add>, <built-in function neg>)
        sage: group_law(SymmetricGroup(5))
        ((), <built-in function mul>, <built-in function inv>)
        sage: group_law(VectorSpace(QQ,3))
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
        raise ValueError("%s does not seem to be a group"%G)

def block_stabilizer(G, B):
    r"""
    Compute the left stabilizer of the block ``B`` under the action of ``G``.

    This function return the list of all `x\in G` such that `x\cdot B=B` (as a
    set).

    INPUT:

    - ``G`` -- a group (additive or multiplicative).

    - ``B`` -- a subset of ``G``.

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
        bb0 = op(b,b0) # bb0 = b-B[0]
        if all(op(bb0,c) in B for c in B):
            S.append(bb0)
    return S


def is_difference_family(G, D, v=None, k=None, l=None, verbose=False):
    r"""
    Check whether ``D`` forms a difference family in the group ``G``.

    INPUT:

    - ``G`` -- group of cardinality ``v``

    - ``D`` -- a set of ``k``-subsets of ``G``

    - ``v``, ``k`` and ``l`` -- optional parameters of the difference family

    - ``verbose`` - whether to print additional information

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

        sage: G = DihedralGroup(8)
        sage: x,y = G.gens()
        sage: i = G.one()
        sage: D1 = [[i,x,x^4], [i,x^2, y*x], [i,x^5,y], [i,x^6,y*x^2], [i,x^7,y*x^5]]
        sage: is_difference_family(G, D1, 16, 3, 2)
        True
        sage: from sage.combinat.designs.bibd import BIBD_from_difference_family
        sage: bibd = BIBD_from_difference_family(G,D1,lambd=2)

    TESTS::

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
    where   = {g: set() for g in Glist}
    del counter[identity]

    for i,d in enumerate(D):
        tmp_counter = {}
        for b in d:
            for c in d:
                if b == c:
                    continue
                gg = mul(b,inv(c)) # = b-c or bc^{-1}
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
    too_few  = []
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
        for g  in too_much:
            print("  {} is obtained {} times in blocks {}".format(
                        g, counter[g], sorted(where[g])))
    if too_few or too_much:
        return False

    if verbose:
        print("It is a ({},{},{})-difference family".format(v, k, l))
    return True

def singer_difference_set(q,d):
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
        sage: G,D = singer_difference_set(3,2)
        sage: is_difference_family(G,D,verbose=True)
        It is a (13,4,1)-difference family
        True

        sage: G,D = singer_difference_set(4,2)
        sage: is_difference_family(G,D,verbose=True)
        It is a (21,5,1)-difference family
        True

        sage: G,D = singer_difference_set(3,3)
        sage: is_difference_family(G,D,verbose=True)
        It is a (40,13,4)-difference family
        True

        sage: G,D = singer_difference_set(9,3)
        sage: is_difference_family(G,D,verbose=True)
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
    c = conway_polynomial(p,e*(d+1))
    if e != 1:  # i.e. q is not a prime, so we factorize c over GF(q) and pick
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

        sage: from sage.combinat.designs.difference_family import is_difference_family, df_q_6_1
        sage: prime_powers = [v for v in range(31,500,30) if is_prime_power(v)]
        sage: parameters = [v for v in prime_powers if df_q_6_1(GF(v,'a'), existence=True) is True]
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

    for c in to_coset: # the loop runs through all nonzero elements of K
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

        sage: D = radical_difference_set(GF(7), 3, 1); D
        [[1, 2, 4]]
        sage: sorted(x-y for x in D[0] for y in D[0] if x != y)
        [1, 2, 3, 4, 5, 6]

        sage: D = radical_difference_set(GF(16,'a'), 6, 2)
        sage: sorted(x-y for x in D[0] for y in D[0] if x != y)
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

        sage: for k in range(2,50):
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
    elif v%4 == 3 and k == (v-1)//2:
        if existence:
            return True
        add_zero = False

    # q = 3 mod 4
    elif v%4 == 3 and k == (v+1)//2:
        if existence:
            return True
        add_zero = True

    # q = 4t^2 + 1, t odd
    elif v%8 == 5 and k == (v-1)//4 and arith.is_square((v-1)//4):
        if existence:
            return True
        add_zero = False

    # q = 4t^2 + 9, t odd
    elif v%8 == 5 and k == (v+3)//4 and arith.is_square((v-9)//4):
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
            print("**  You found a new example of radical difference set **\n"\
                  "**  for the parameters (v,k,l)=({},{},{}).            **\n"\
                  "**  Please contact sage-devel@googlegroups.com        **\n".format(v, k, l))
            if existence:
                return True
            add_zero = False

        else:
            D = K.cyclotomic_cosets(x**((v-1)//(k-1)), [K.one()])
            D[0].insert(0,K.zero())
            if is_difference_family(K, D, v, k, l):
                print("**  You found a new example of radical difference set **\n"\
                      "**  for the parameters (v,k,l)=({},{},{}).            **\n"\
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

def one_cyclic_tiling(A,n):
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
    if len(set(a%d for a in A)) == d:
        return [i*d for i in range(n//d)]

    # next, we consider an exhaustive search
    from sage.combinat.dlx import DLXMatrix

    rows = []
    for i in range(n):
        rows.append([i+1, [(i+a)%n+1 for a in A]])
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

    - ``K`` -- a finite field of cardinality `q`.

    - ``k`` -- a positive integer so that `k(k-1)` divides `q-1`.

    OUTPUT:

    Either a difference family or ``None`` if it does not exist.

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

        sage: one_radical_difference_family(GF(13),4)
        [[0, 1, 3, 9]]

    The parameters that appear in [Bu95]_::

        sage: df = one_radical_difference_family(GF(449), 8); df
        [[0, 1, 18, 25, 176, 324, 359, 444],
         [0, 9, 88, 162, 222, 225, 237, 404],
         [0, 11, 140, 198, 275, 357, 394, 421],
         [0, 40, 102, 249, 271, 305, 388, 441],
         [0, 49, 80, 93, 161, 204, 327, 433],
         [0, 70, 99, 197, 230, 362, 403, 435],
         [0, 121, 141, 193, 293, 331, 335, 382],
         [0, 191, 285, 295, 321, 371, 390, 392]]
        sage: is_difference_family(GF(449), df, 449, 8, 1)
        True
    """
    q = K.cardinality()
    x = K.multiplicative_generator()

    e = k*(k-1)
    if q%e != 1:
        raise ValueError("q%e is not 1")

    # We define A by (see the function's documentation):
    # ΔB = C.A
    if k%2 == 1:
        m = (k-1) // 2
        r = x ** ((q-1) // k)     # k-th root of unity
        A = [r**i - 1 for i in range(1,m+1)]
    else:
        m = k // 2
        r = x ** ((q-1) // (k-1)) # (k-1)-th root of unity
        A = [r**i - 1 for i in range(1,m)]
        A.append(K.one())

    # instead of the complicated multiplicative group K^*/(±C) we use the
    # discrete logarithm to convert everything into the additive group Z/cZ
    c = m * (q-1) // e # cardinal of ±C
    from sage.groups.generic import discrete_log
    logA = [discrete_log(a,x)%c for a in A]

    # if two elements of A are equal modulo c then no tiling is possible
    if len(set(logA)) != m:
        return None

    # brute force
    tiling = one_cyclic_tiling(logA, c)
    if tiling is None:
        return None

    D = K.cyclotomic_cosets(r, [x**i for i in tiling])
    if k%2 == 0:
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

    - ``K`` - a finite field

    - ``k`` -- positive integer, the size of the blocks

    - ``l`` -- the `\lambda` parameter (default to `1`)

    - ``existence`` -- if ``True``, then return either ``True`` if Sage knows
      how to build such design, ``Unknown`` if it does not and ``False`` if it
      knows that the design does not exist.

    - ``check`` -- boolean (default: ``True``). If ``True`` then the result of
      the computation is checked before being returned. This should not be
      needed but ensures that the output is correct.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import radical_difference_family

        sage: radical_difference_family(GF(73),9)
        [[1, 2, 4, 8, 16, 32, 37, 55, 64]]

        sage: radical_difference_family(GF(281),5)
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

        sage: for k in range(5,10):
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

    - ``check`` -- boolean (default: ``True``). If ``True`` then the result of
      the computation is checked before being returned. This should not be
      needed but ensures that the output is correct.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import twin_prime_powers_difference_set
        sage: G,D = twin_prime_powers_difference_set(3)
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

    - ``v``, ``k``, ``lmbda`` - (integers) parameters of the difference family

    - ``return_parameters`` -- (boolean, default ``False``) if ``True`` return a
      pair ``(True, (q, s))`` so that ``(q,s)`` can be used in the function
      :func:`mcfarland_1973_construction` to actually build a
      ``(v,k,lmbda)``-difference family. Or ``(False, None)`` if the
      construction is not possible.

    EXAMPLES::

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
        sage: for v in range(1, 100):
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
        return (False,None) if return_parameters else False
    k = ZZ(k)
    lmbda = ZZ(lmbda)
    qs,r = (k - lmbda).sqrtrem() # sqrt(k-l) should be q^s
    if r or (qs*(qs-1))%lmbda:
        return (False,None) if return_parameters else False

    q = qs*(qs-1) // lmbda + 1
    if (q <= 1 or
        v * (q-1) != qs*q * (qs*q+q-2)  or
        k * (q-1)!= qs * (qs*q-1)):
        return (False,None) if return_parameters else False

    # NOTE: below we compute the value of s so that qs = q^s. If the method
    # is_power_of of integers would be able to return the exponent, we could use
    # that... but currently this is not the case
    # see trac ticket #19792
    p1,a1 = qs.is_prime_power(get_data=True)
    p2,a2 = q.is_prime_power(get_data=True)

    if a1 == 0 or a2 == 0 or p1 != p2 or a1%a2:
        return (False,None) if return_parameters else False

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

    - ``q``, ``s`` - (integers) parameters for the difference set (see the above
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

        sage: G,D = mcfarland_1973_construction(3, 1)
        sage: assert is_difference_family(G, D, 45, 12, 3)

        sage: G,D = mcfarland_1973_construction(2, 2)
        sage: assert is_difference_family(G, D, 64, 28, 12)
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
            D.append(G((tuple(v) + (k,))))

    return G,[D]

def are_hadamard_difference_set_parameters(v, k, lmbda):
    r"""
    Check whether ``(v,k,lmbda)`` is of the form ``(4N^2, 2N^2 - N, N^2 - N)``.

    INPUT:

    - ``(v,k,lmbda)`` -- parameters of a difference set

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

    This function looks for two integers `N_1` and `N_2`` greater than `1`
    and so that `N = 2 N_1 N_2` and there exists Hadamard difference set with
    parameters `(4 N_i^2, 2N_i^2 - N_i, N_i^2 - N_i)`. If such pair exists,
    the output is the pair ``(N_1, N_2)`` otherwise it is ``None``.

    INPUT:

    - ``N`` -- positive integer

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import hadamard_difference_set_product_parameters
        sage: hadamard_difference_set_product_parameters(8)
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

    - ``G1,D1``, ``G2,D2`` -- two Hadamard difference sets

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import hadamard_difference_set_product
        sage: from sage.combinat.designs.difference_family import is_difference_family

        sage: G1,D1 = designs.difference_family(16,6,2)
        sage: G2,D2 = designs.difference_family(36,15,6)

        sage: G11,D11 = hadamard_difference_set_product(G1,D1,G1,D1)
        sage: assert is_difference_family(G11, D11, 256, 120, 56)
        sage: assert designs.difference_family(256, 120, 56, existence=True) is True

        sage: G12,D12 = hadamard_difference_set_product(G1,D1,G2,D2)
        sage: assert is_difference_family(G12, D12, 576, 276, 132)
        sage: assert designs.difference_family(576, 276, 132, existence=True) is True
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
      \times C_3 \times C_4`).

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

    L = [[(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)], # complement of y=0
         [(0,0),(1,1),(2,2)],                   # x-y=0
         [(0,0),(1,2),(2,1)],                   # x+y=0
         [(0,0),(0,1),(0,2)]]                   # x=0

    return G, [[G(v + k) for l, k in zip(L, K) for v in l]]


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

    If there is no such difference family, an ``EmptySetError`` is raised and if
    there is no construction at the moment ``NotImplementedError`` is raised.

    INPUT:

    - ``v,k,l`` -- parameters of the difference family. If ``l`` is not provided
      it is assumed to be ``1``.

    - ``existence`` -- if ``True``, then return either ``True`` if Sage knows
      how to build such design, ``Unknown`` if it does not and ``False`` if it
      knows that the design does not exist.

    - ``explain_construction`` -- instead of returning a difference family,
      returns a string that explains the construction used.

    - ``check`` -- boolean (default: ``True``). If ``True`` then the result of
      the computation is checked before being returned. This should not be
      needed but ensures that the output is correct.

    OUTPUT:

    A pair ``(G,D)`` made of a group `G` and a difference family `D` on that
    group. Or, if ``existence`` is ``True`` a troolean or if
    ``explain_construction`` is ``True`` a string.

    EXAMPLES::

        sage: G,D = designs.difference_family(73,4)
        sage: G
        Finite Field of size 73
        sage: D
        [[0, 1, 5, 18],
         [0, 3, 15, 54],
         [0, 9, 45, 16],
         [0, 27, 62, 48],
         [0, 8, 40, 71],
         [0, 24, 47, 67]]

        sage: print(designs.difference_family(73, 4, explain_construction=True))
        The database contains a (73,4)-evenly distributed set

        sage: G,D = designs.difference_family(15,7,3)
        sage: G
        Ring of integers modulo 15
        sage: D
        [[0, 1, 2, 4, 5, 8, 10]]
        sage: print(designs.difference_family(15,7,3,explain_construction=True))
        Singer difference set

        sage: print(designs.difference_family(91,10,1,explain_construction=True))
        Singer difference set
        sage: print(designs.difference_family(64,28,12, explain_construction=True))
        McFarland 1973 construction
        sage: print(designs.difference_family(576, 276, 132, explain_construction=True))
        Hadamard difference set product from N1=2 and N2=3

    For `k=6,7` we look at the set of small prime powers for which a
    construction is available::

        sage: def prime_power_mod(r,m):
        ....:     k = m+r
        ....:     while True:
        ....:         if is_prime_power(k):
        ....:             yield k
        ....:         k += m

        sage: from itertools import islice
        sage: l6 = {True:[], False: [], Unknown: []}
        sage: for q in islice(prime_power_mod(1,30), int(60)):
        ....:     l6[designs.difference_family(q,6,existence=True)].append(q)
        sage: l6[True]
        [31, 121, 151, 181, 211, ...,  3061, 3121, 3181]
        sage: l6[Unknown]
        [61]
        sage: l6[False]
        []

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

        sage: for v in range(2,100):
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
        15: (3,1), (4,6), (5,6), (7,3)
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
        sage: for Q,k in [(Q4,4),(Q5,5),(Q8,8),(Q9,9),(Q15,15)]:
        ....:     for q in Q:
        ....:         assert designs.difference_family(q,k,1,existence=True) is True
        ....:         _ = designs.difference_family(q,k,1)

    Check Singer difference sets::

        sage: sgp = lambda q,d: ((q**(d+1)-1)//(q-1), (q**d-1)//(q-1), (q**(d-1)-1)//(q-1))

        sage: for q in range(2,10):
        ....:     if is_prime_power(q):
        ....:         for d in [2,3,4]:
        ....:           v,k,l = sgp(q,d)
        ....:           assert designs.difference_family(v,k,l,existence=True) is True
        ....:           _ = designs.difference_family(v,k,l)

    Check twin primes difference sets::

        sage: for p in [3,5,7,9,11]:
        ....:     v = p*(p+2); k = (v-1)/2;  lmbda = (k-1)/2
        ....:     G,D = designs.difference_family(v,k,lmbda)

    Check the database::

        sage: from sage.combinat.designs.database import DF,EDS
        sage: for v,k,l in DF:
        ....:     assert designs.difference_family(v,k,l,existence=True) is True
        ....:     df = designs.difference_family(v,k,l,check=True)

        sage: for k in EDS:
        ....:     for v in EDS[k]:
        ....:         assert designs.difference_family(v,k,1,existence=True) is True
        ....:         df = designs.difference_family(v,k,1,check=True)

    Check the known Hadamard parameters::

        sage: for N in range(2,21):
        ....:     v = 4*N^2; k = 2*N^2-N; l = N^2-N
        ....:     status = designs.difference_family(v,k,l,existence=True)
        ....:     print("{:2} {}".format(N,designs.difference_family(v,k,l,explain_construction=True) if status is True else status))
        2 McFarland 1973 construction
        3 Turyn 1965 construction
        4 McFarland 1973 construction
        5 False
        6 Unknown
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

    Check a failing construction (:trac:`17528`)::

        sage: designs.difference_family(9,3)
        Traceback (most recent call last):
        ...
        NotImplementedError: No construction available for (9,3,1)-difference family

    Check that when ``existence=True`` we always obtain ``True``, ``False`` or ``Unknown``,
    and when ``explain_construction=True``, it is a string (see :trac:`24513`)::

        sage: designs.difference_family(3, 2, 1, existence=True)
        True
        sage: designs.difference_family(3, 2, 1, explain_construction=True)
        'Trivial difference family'

        sage: for _ in range(100):
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
                     "sage-devel@googlegroups.com".format(v,k,l))
        return G,df

    if k in [0,1]:
        # Then \Delta D_i is empty
        # So if G\{0} is empty is good, otherwise not
        if v == 1:
            if existence:
                return True
            from sage.rings.finite_rings.integer_mod_ring import Zmod
            l = [0] if k ==1 else []
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
