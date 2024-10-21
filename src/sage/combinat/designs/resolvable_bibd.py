# sage_setup: distribution = sagemath-graphs
# sage.doctest: needs sage.rings.finite_rings
r"""
Resolvable Balanced Incomplete Block Design (RBIBD)

This module contains everything related to resolvable Balanced Incomplete Block
Designs. The constructions implemented here can be obtained through the
``designs.<tab>`` object::

    designs.resolvable_balanced_incomplete_block_design(15,3)

For Balanced Incomplete Block Design (BIBD) see the module :mod:`bibd
<sage.combinat.designs.bibd>`. A BIBD
is said to be *resolvable* if its blocks can be partitionned into parallel
classes, i.e.  partitions of its ground set.

The main function of this module is
:func:`resolvable_balanced_incomplete_block_design`, which calls all others.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`resolvable_balanced_incomplete_block_design` | Return a resolvable BIBD of parameters `v,k`.
    :func:`kirkman_triple_system` | Return a Kirkman Triple System on `v` points.
    :func:`v_4_1_rbibd` | Return a `(v,4,1)`-RBIBD
    :func:`PBD_4_7` | Return a `(v,\{4,7\})`-PBD
    :func:`PBD_4_7_from_Y` | Return a `(3v+1,\{4,7\})`-PBD from a `(v,\{4,5,7\},\NN-\{3,6,10\})`-GDD.

References:

.. [Stinson91] \D.R. Stinson,
   A survey of Kirkman triple systems and related designs,
   Volume 92, Issues 1-3, 17 November 1991, Pages 371-393,
   Discrete Mathematics,
   :doi:`10.1016/0012-365X(91)90294-C`

.. [RCW71] \D. K. Ray-Chaudhuri, R. M. Wilson,
   Solution of Kirkman's schoolgirl problem,
   Volume 19, Pages 187-203,
   Proceedings of Symposia in Pure Mathematics

.. [BJL99] \T. Beth, D. Jungnickel, H. Lenz,
   Design Theory 2ed.
   Cambridge University Press
   1999

Functions
---------
"""
from itertools import repeat
from sage.arith.misc import is_prime_power
from sage.combinat.designs.bibd import BalancedIncompleteBlockDesign
from sage.categories.sets_cat import EmptySetError
from .bibd import balanced_incomplete_block_design
from sage.misc.unknown import Unknown


def resolvable_balanced_incomplete_block_design(v,k,existence=False):
    r"""
    Return a resolvable BIBD of parameters `v,k`.

    A BIBD is said to be *resolvable* if its blocks can be partitionned into
    parallel classes, i.e. partitions of the ground set.

    INPUT:

    - ``v``, ``k`` -- integers

    - ``existence`` -- boolean; instead of building the design, return:

        - ``True`` -- meaning that Sage knows how to build the design

        - ``Unknown`` -- meaning that Sage does not know how to build the
          design, but that the design may exist (see :mod:`sage.misc.unknown`)

        - ``False`` -- meaning that the design does not exist

    .. SEEALSO::

        - :meth:`IncidenceStructure.is_resolvable`
        - :func:`~sage.combinat.designs.bibd.balanced_incomplete_block_design`

    EXAMPLES::

        sage: KTS15 = designs.resolvable_balanced_incomplete_block_design(15,3); KTS15
        (15,3,1)-Balanced Incomplete Block Design
        sage: KTS15.is_resolvable()
        True

    TESTS::

        sage: for v in range(40):
        ....:     for k in range(v):
        ....:         if designs.resolvable_balanced_incomplete_block_design(v,k,existence=True) is True:
        ....:             _ = designs.resolvable_balanced_incomplete_block_design(v,k)
    """
    # Trivial cases
    if v == 1 or k == v:
        return balanced_incomplete_block_design(v,k,existence=existence)

    # Non-existence of resolvable BIBD
    if (v < k or
        k < 2 or
        v % k != 0 or
        (v-1) % (k-1) != 0 or
        (v*(v-1)) % (k*(k-1)) != 0 or
        # From the Handbook of combinatorial designs:
        #
        # With lambda>1 the other exceptions is
        # (15,5,2)
        (k == 6 and v == 36) or
        # Fisher's inequality
        (v*(v-1))/(k*(k-1)) < v):
        if existence:
            return False
        raise EmptySetError("There exists no ({},{},{})-RBIBD".format(v,k,1))

    if k == 2:
        if existence:
            return True
        classes = [[[(c+i) % (v-1),(c+v-i) % (v-1)] for i in range(1, v//2)]
                   for c in range(v-1)]
        for i,classs in enumerate(classes):
            classs.append([v-1,i])

        B = BalancedIncompleteBlockDesign(v,
                                          sum(classes,[]),
                                          k=k,
                                          check=True,
                                          copy=False)
        B._classes = classes
        return B
    elif k == 3:
        return kirkman_triple_system(v,existence=existence)
    elif k == 4:
        return v_4_1_rbibd(v,existence=existence)
    else:
        if existence:
            return Unknown
        raise NotImplementedError("I don't know how to build a ({},{},1)-RBIBD!".format(v,3))


def kirkman_triple_system(v,existence=False):
    r"""
    Return a Kirkman Triple System on `v` points.

    A Kirkman Triple System `KTS(v)` is a resolvable Steiner Triple System. It
    exists if and only if `v\equiv 3\pmod{6}`.

    INPUT:

    - ``n`` -- integer

    - ``existence`` -- boolean (default: ``False``); whether to build the
      `KTS(n)` or only answer whether it exists

    .. SEEALSO::

        :meth:`IncidenceStructure.is_resolvable`

    EXAMPLES:

    A solution to Kirkmman's original problem::

        sage: kts = designs.kirkman_triple_system(15)
        sage: classes = kts.is_resolvable(1)[1]
        sage: names = '0123456789abcde'
        sage: def to_name(r_s_t):
        ....:     r, s, t = r_s_t
        ....:     return ' ' + names[r] + names[s] + names[t] + ' '
        sage: rows = ['   '.join(('Day {}'.format(i) for i in range(1,8)))]
        sage: rows.extend('   '.join(map(to_name,row)) for row in zip(*classes))
        sage: print('\n'.join(rows))
        Day 1   Day 2   Day 3   Day 4   Day 5   Day 6   Day 7
         07e     18e     29e     3ae     4be     5ce     6de
         139     24a     35b     46c     05d     167     028
         26b     03c     14d     257     368     049     15a
         458     569     06a     01b     12c     23d     347
         acd     7bd     78c     89d     79a     8ab     9bc

    TESTS::

        sage: for i in range(3,300,6):                                                  # needs sage.combinat
        ....:     _ = designs.kirkman_triple_system(i)
    """
    if v % 6 != 3:
        if existence:
            return False
        raise ValueError("There is no KTS({}) as v!=3 mod(6)".format(v))

    if existence:
        return False

    elif v == 3:
        return BalancedIncompleteBlockDesign(3,[[0,1,2]],k=3,lambd=1)

    elif v == 9:
        classes = [[[0, 1, 5], [2, 6, 7], [3, 4, 8]],
                   [[1, 6, 8], [3, 5, 7], [0, 2, 4]],
                   [[1, 4, 7], [0, 3, 6], [2, 5, 8]],
                   [[4, 5, 6], [0, 7, 8], [1, 2, 3]]]
        KTS = BalancedIncompleteBlockDesign(v,[tr for cl in classes for tr in cl],k=3,lambd=1,copy=False)
        KTS._classes = classes
        return KTS

    # Construction 1.1 from [Stinson91] (originally Theorem 6 from [RCW71])
    #
    # For all prime powers q=1 mod 6, there exists a KTS(2q+1)
    elif ((v-1)//2) % 6 == 1 and is_prime_power((v-1)//2):
        from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
        q = (v-1)//2
        K = GF(q,'x')
        a = K.primitive_element()
        t = (q - 1) // 6

        # m is the solution of a^m=(a^t+1)/2
        from sage.groups.generic import discrete_log
        m = discrete_log((a**t+1)/2, a)
        assert 2*a**m == a**t+1

        # First parallel class
        first_class = [[(0,1),(0,2),'inf']]
        b0 = K.one()
        b1 = a**t
        b2 = a**m
        first_class.extend([(b0*a**i,1),(b1*a**i,1),(b2*a**i,2)]
                            for i in list(range(t))+list(range(2*t,3*t))+list(range(4*t,5*t)))
        b0 = a**(m+t)
        b1 = a**(m+3*t)
        b2 = a**(m+5*t)
        first_class.extend([[(b0*a**i,2),(b1*a**i,2),(b2*a**i,2)]
                            for i in range(t)])

        # Action of K on the points
        action = lambda v,x : (v+x[0],x[1]) if len(x) == 2 else x

        # relabel to integer
        relabel = {(p,x): i+(x-1)*q
                   for i,p in enumerate(K)
                   for x in [1,2]}
        relabel['inf'] = 2*q

        classes = [[[relabel[action(p,x)] for x in tr] for tr in first_class]
                   for p in K]

        KTS = BalancedIncompleteBlockDesign(v,[tr for cl in classes for tr in cl],k=3,lambd=1,copy=False)

        KTS._classes = classes
        return KTS

    # Construction 1.2 from [Stinson91] (originally Theorem 5 from [RCW71])
    #
    # For all prime powers q=1 mod 6, there exists a KTS(3q)
    elif (v//3) % 6 == 1 and is_prime_power(v//3):
        from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
        q = v//3
        K = GF(q,'x')
        a = K.primitive_element()
        t = (q - 1) // 6
        A0 = [(0,0),(0,1),(0,2)]
        B = [[(a**i,j),(a**(i+2*t),j),(a**(i+4*t),j)] for j in range(3)
              for i in range(t)]
        A = [[(a**i,0),(a**(i+2*t),1),(a**(i+4*t),2)] for i in range(6*t)]

        # Action of K on the points
        action = lambda v,x: (v+x[0],x[1])

        # relabel to integer
        relabel = {(p,j): i+j*q
                   for i,p in enumerate(K)
                   for j in range(3)}

        B0 = [A0] + B + A[t:2*t] + A[3*t:4*t] + A[5*t:6*t]

        # Classes
        classes = [[[relabel[action(p,x)] for x in tr] for tr in B0]
                   for p in K]

        for i in list(range(t))+list(range(2*t,3*t))+list(range(4*t,5*t)):
            classes.append([[relabel[action(p,x)] for x in A[i]] for p in K])

        KTS = BalancedIncompleteBlockDesign(v,[tr for cl in classes for tr in cl],k=3,lambd=1,copy=False)
        KTS._classes = classes
        return KTS

    else:
        # This is Lemma IX.6.4 from [BJL99].
        #
        # This construction takes a (v,{4,7})-PBD. All points are doubled (x has
        # a copy x'), and an infinite point \infty is added.
        #
        # On all blocks of 2*4 points we "paste" a KTS(2*4+1) using the infinite
        # point, in such a way that all {x,x',infty} are set of the design. We
        # do the same for blocks with 2*7 points using a KTS(2*7+1).
        #
        # Note that the triples of points equal to {x,x',\infty} will be added
        # several times.
        #
        # As all those subdesigns are resolvable, each class of the KTS(n) is
        # obtained by considering a set {x,x',\infty} and all sets of all
        # parallel classes of the subdesign which contain this set.

        # We create the small KTS(n') we need, and relabel them such that
        # 01(n'-1),23(n'-1),... are blocks of the design.
        gdd4 = kirkman_triple_system(9)
        gdd7 = kirkman_triple_system(15)

        X = [B for B in gdd4 if 8 in B]
        for b in X:
            b.remove(8)
        X = sum(X, []) + [8]
        gdd4.relabel({v:i for i,v in enumerate(X)})
        gdd4 = gdd4.is_resolvable(True)[1]  # the relabeled classes

        X = [B for B in gdd7 if 14 in B]
        for b in X:
            b.remove(14)
        X = sum(X, []) + [14]
        gdd7.relabel({v:i for i,v in enumerate(X)})
        gdd7 = gdd7.is_resolvable(True)[1]  # the relabeled classes

        # The first parallel class contains 01(n'-1), the second contains
        # 23(n'-1), etc..
        # Then remove the blocks containing (n'-1)
        for B in gdd4:
            for i, b in enumerate(B):
                if 8 in b:
                    j = min(b)
                    del B[i]
                    B.insert(0, j)
                    break
        gdd4.sort()
        for B in gdd4:
            B.pop(0)

        for B in gdd7:
            for i, b in enumerate(B):
                if 14 in b:
                    j = min(b)
                    del B[i]
                    B.insert(0, j)
                    break
        gdd7.sort()
        for B in gdd7:
            B.pop(0)

        # Pasting the KTS(n') without {x,x',\infty} blocks
        classes = [[] for _ in repeat(None, (v - 1) // 2)]
        gdd = {4: gdd4, 7: gdd7}
        for B in PBD_4_7((v-1)//2,check=False):
            for i,classs in enumerate(gdd[len(B)]):
                classes[B[i]].extend([[2*B[x//2]+x % 2 for x in BB] for BB in classs])

        # The {x,x',\infty} blocks
        for i,classs in enumerate(classes):
            classs.append([2*i,2*i+1,v-1])

        KTS = BalancedIncompleteBlockDesign(v,
                                            blocks=[tr for cl in classes for tr in cl],
                                            k=3,
                                            lambd=1,
                                            check=True,
                                            copy=False)
        KTS._classes = classes
        assert KTS.is_resolvable()

        return KTS


def v_4_1_rbibd(v,existence=False):
    r"""
    Return a `(v,4,1)`-RBIBD.

    INPUT:

    - ``n`` -- integer

    - ``existence`` -- boolean (default: ``False``); whether to build the
      design or only answer whether it exists

    .. SEEALSO::

        - :meth:`IncidenceStructure.is_resolvable`
        - :func:`resolvable_balanced_incomplete_block_design`

    .. NOTE::

        A resolvable `(v,4,1)`-BIBD exists whenever `1\equiv 4\pmod(12)`. This
        function, however, only implements a construction of `(v,4,1)`-BIBD such
        that `v=3q+1\equiv 1\pmod{3}` where `q` is a prime power (see VII.7.5.a
        from [BJL99]_).

    EXAMPLES::

        sage: rBIBD = designs.resolvable_balanced_incomplete_block_design(28,4)
        sage: rBIBD.is_resolvable()
        True
        sage: rBIBD.is_t_design(return_parameters=True)
        (True, (2, 28, 4, 1))

    TESTS::

        sage: for q in prime_powers(2,30):  # indirect doctest
        ....:     if (3*q+1)%12 == 4:
        ....:         _ = designs.resolvable_balanced_incomplete_block_design(3*q+1,4)
    """
    # Volume 1, VII.7.5.a from [BJL99]_
    if v % 3 != 1 or not is_prime_power((v-1)//3):
        if existence:
            return Unknown
        raise NotImplementedError("I don't know how to build a ({},{},1)-RBIBD!".format(v,4))
    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
    q = (v-1)//3
    nn = (q-1)//4
    G = GF(q,'x')
    w = G.primitive_element()
    e = w**(nn)
    assert e**2 == -1

    first_class = [[(w**i,j),(-w**i,j),(e*w**i,j+1),(-e*w**i,j+1)]
                   for i in range(nn) for j in range(3)]

    first_class.append([(0,0),(0,1),(0,2),'inf'])

    label = {p:i for i,p in enumerate(G)}

    classes = [[[v-1 if x == 'inf' else (x[1] % 3)*q+label[x[0]+g] for x in S]
                for S in first_class]
               for g in G]

    BIBD = BalancedIncompleteBlockDesign(v,
                                         blocks=sum(classes,[]),
                                         k=4,
                                         check=True,
                                         copy=False)
    BIBD._classes = classes
    assert BIBD.is_resolvable()
    return BIBD


def PBD_4_7(v,check=True, existence=False):
    r"""
    Return a `(v,\{4,7\})`-PBD.

    For all `v` such that `n\equiv 1\pmod{3}` and `n\neq 10,19, 31` there exists
    a `(v,\{4,7\})`-PBD. This is proved in Proposition IX.4.5 from [BJL99]_,
    which this method implements.

    This construction of PBD is used by the construction of Kirkman Triple
    Systems.

    EXAMPLES::

        sage: from sage.combinat.designs.resolvable_bibd import PBD_4_7
        sage: PBD_4_7(22)
        Pairwise Balanced Design on 22 points with sets of sizes in [4, 7]

    TESTS:

    All values `\leq 300`::

        sage: for i in range(1,300,3):                                                  # needs sage.schemes
        ....:     if i not in [10,19,31]:
        ....:         assert PBD_4_7(i,existence=True) is True
        ....:         _ = PBD_4_7(i,check=True)
    """
    if v % 3 != 1 or v in [10,19,31]:
        if existence:
            return Unknown
        raise NotImplementedError
    if existence:
        return True

    from .group_divisible_designs import GroupDivisibleDesign
    from .group_divisible_designs import GDD_4_2
    from .bibd import PairwiseBalancedDesign
    from .bibd import balanced_incomplete_block_design

    if v == 22:
        # Beth/Jungnickel/Lenz: take KTS(15) and extend each of the 7 classes
        # with a new point. Make those new points a 7-set.
        KTS15 = kirkman_triple_system(15)
        blocks = [S+[i+15] for i,classs in enumerate(KTS15._classes) for S in classs]+[list(range(15, 22))]

    elif v == 34:
        # [BJL99] (p527,vol1), but originally Brouwer
        A = [(0,0),(1,1),(2,0),(4,1)]
        B = [(0,0),(1,0),(4,2)]
        C = [(0,0),(2,2),(5,0)]
        D = [(0,0),(0,1),(0,2)]

        A = [[(x+i,  y+j)     for x,y in A]        for i in range(9) for j in range(3)]
        B = [[(x+i,  y+i+j)   for x,y in B]+[27+j] for i in range(9) for j in range(3)]
        C = [[(x+i+j,y+2*i+j) for x,y in C]+[30+j] for i in range(9) for j in range(3)]
        D = [[(x+i,  y+i)     for x,y in D]+[33]   for i in range(9)]

        blocks = [[int(x) if not isinstance(x,tuple) else (x[1] % 3)*9+(x[0] % 9) for x in S]
                  for S in A+B+C+D+[list(range(27,34))]]
    elif v == 46:
        # [BJL99] (p527,vol1), but originally Brouwer
        A = [(1,0),(3,0),(9,0),(0,1)]
        B = [(2,0),(6,0),(5,0),(0,1)]
        C = [(0,0),(1,1),(4,2)]
        D = [(0,0),(2,1),(7,2)]
        E = [(0,0),(0,1),(0,2)]

        A = [[(x+i,  y+j) for x,y in A]        for i in range(13) for j in range(3)]
        B = [[(x+i,  y+j) for x,y in B]        for i in range(13) for j in range(3)]
        C = [[(x+i,  y+j) for x,y in C]+[39+j] for i in range(13) for j in range(3)]
        D = [[(x+i,  y+j) for x,y in D]+[42+j] for i in range(13) for j in range(3)]
        E = [[(x+i,  y+i) for x,y in E]+[45]   for i in range(13)]

        blocks = [[int(x) if not isinstance(x,tuple) else (x[1] % 3)*13+(x[0] % 13) for x in S]
                  for S in A+B+C+D+E+[list(range(39, 46))]]

    elif v == 58:
        # [BJL99] (p527,vol1), but originally Brouwer
        A = [(0,0),(1,0),(4,0),( 5,1)]
        B = [(0,0),(2,0),(8,0),(11,1)]
        C = [(0,0),(5,0),(2,1),(12,1)]
        D = [(0,0),(8,1),(7,2)]
        E = [(0,0),(6,1),(4,2)]
        F = [(0,0),(0,1),(0,2)]

        A = [[(x+i,  y+j) for x,y in A]        for i in range(17) for j in range(3)]
        B = [[(x+i,  y+j) for x,y in B]        for i in range(17) for j in range(3)]
        C = [[(x+i,  y+j) for x,y in C]        for i in range(17) for j in range(3)]
        D = [[(x+i,  y+j) for x,y in D]+[51+j] for i in range(17) for j in range(3)]
        E = [[(x+i,  y+j) for x,y in E]+[54+j] for i in range(17) for j in range(3)]
        F = [[(x+i,  y+i) for x,y in F]+[57]   for i in range(17)]

        blocks = [[int(x) if not isinstance(x,tuple) else (x[1] % 3)*17+(x[0] % 17) for x in S]
                  for S in A+B+C+D+E+F+[list(range(51,58))]]

    elif v == 70:
        # [BJL99] (p527,vol1), but originally Brouwer
        A = [(0,0),(1,0),(5,1),(13,1)]
        B = [(0,0),(4,0),(20,1),(10,1)]
        C = [(0,0),(16,0),(17,1),(19,1)]
        D = [(0,0),(2,1),(8,1),(11,1)]
        E = [(0,0),(3,2),(9,1)]
        F = [(0,0),(7,0),(14,1)]
        H = [(0,0),(0,1),(0,2)]

        A = [[(x+i,  y+j)       for x,y in A]        for i in range(21) for j in range(3)]
        B = [[(x+i,  y+j)       for x,y in B]        for i in range(21) for j in range(3)]
        C = [[(x+i,  y+j)       for x,y in C]        for i in range(21) for j in range(3)]
        D = [[(x+i,  y+j)       for x,y in D]        for i in range(21) for j in range(3)]
        E = [[(x+i,  y+j)       for x,y in E]+[63+j] for i in range(21) for j in range(3)]
        F = [[(x+3*i+j, y+ii+j) for x,y in F]+[66+j] for i in range( 7) for j in range(3) for ii in range(3)]
        H = [[(x+i,  y+i)       for x,y in H]+[69]   for i in range(21)]

        blocks = [[int(x) if not isinstance(x,tuple) else (x[1] % 3)*21+(x[0] % 21)
                   for x in S]
                  for S in A+B+C+D+E+F+H+[list(range(63,70))]]

    elif v == 82:
        # This construction is Theorem IX.3.16 from [BJL99] (p.627).
        #
        # A (15,{4},{3})-GDD from a (16,4)-BIBD
        from .group_divisible_designs import group_divisible_design
        from .orthogonal_arrays       import transversal_design
        GDD = group_divisible_design(3*5,K=[4],G=[3],check=False)
        TD = transversal_design(5,5)

        # A (75,{4},{15})-GDD
        GDD2 = [[3*B[x//3]+x % 3 for x in BB] for B in TD for BB in GDD]

        # We now complete the (75,{4},{15})-GDD into a (82,{4,7})-PBD. For this,
        # we add 7 new points that are added to all groups of size 15.
        #
        # On these groups a (15+7,{4,7})-PBD is pasted, in such a way that the 7
        # new points are a set of the final PBD
        PBD22 = PBD_4_7(15+7)
        S = next(SS for SS in PBD22 if len(SS) == 7)  # a set of size 7
        PBD22.relabel({v:i for i,v in enumerate([i for i in range(15+7) if i not in S] + S)})

        for B in PBD22:
            if B == S:
                continue
            for i in range(5):
                GDD2.append([x+i*15 if x < 15 else x+60 for x in B])

        GDD2.append(list(range(75,82)))
        blocks = GDD2

    elif v == 94:
        # IX.4.5.l from [BJL99].
        #
        # take 4 parallel lines from an affine plane of order 7, and a 5th
        # one. This is a (31,{4,5,7})-BIBD. And 94=3*31+1.
        from sage.combinat.designs.block_design import AffineGeometryDesign
        AF = AffineGeometryDesign(2,1,7)
        parall = []
        plus_one = None
        for S in AF:
            if all(x not in SS for SS in parall for x in S):
                parall.append(S)
            elif plus_one is None:
                plus_one = S
            if len(parall) == 4 and plus_one is not None:
                break
        X = set(sum(parall, plus_one))

        S_4_5_7 = [X.intersection(S) for S in AF]
        S_4_5_7 = [S for S in S_4_5_7 if len(S) > 1]
        S_4_5_7 = PairwiseBalancedDesign(X,
                                         blocks=S_4_5_7,
                                         K=[4,5,7],
                                         check=False)
        S_4_5_7.relabel()
        return PBD_4_7_from_Y(S_4_5_7,check=check)

    elif v == 127 or v == 142:
        # IX.4.5.o from [BJL99].
        #
        # Attach two or seven infinite points to a (40,4)-RBIBD to get a
        # (42,{4,5},{1,2,7})-GDD or a (47,{4,5},{1,2,7})-GDD
        points_to_add = 2 if v == 127 else 7
        rBIBD4 = v_4_1_rbibd(40)
        GDD = [S+[40+i] if i < points_to_add else S
               for i,classs in enumerate(rBIBD4._classes)
               for S in classs]
        if points_to_add == 7:
            GDD.append(list(range(40, 40 + points_to_add)))
            groups = [[x] for x in range(40+points_to_add)]
        else:
            groups = [[x] for x in range(40)]
            groups.append(list(range(40,40+points_to_add)))
        GDD = GroupDivisibleDesign(40+points_to_add,
                                   groups=groups,
                                   blocks=GDD,
                                   K=[2,4,5,7],
                                   check=False,
                                   copy=False)

        return PBD_4_7_from_Y(GDD,check=check)

    elif v % 6 == 1 and GDD_4_2((v - 1) // 6, existence=True) is True:
        # VII.5.17 from [BJL99]
        gdd = GDD_4_2((v - 1) // 6)
        return PBD_4_7_from_Y(gdd, check=check)

    elif v == 202:
        # IV.4.5.p from [BJL99]
        PBD = PBD_4_7(22,check=False)
        PBD = PBD_4_7_from_Y(PBD,check=False)
        return PBD_4_7_from_Y(PBD,check=check)

    elif balanced_incomplete_block_design(v,4,existence=True) is True:
        return balanced_incomplete_block_design(v,4)
    elif balanced_incomplete_block_design(v,7,existence=True) is True:
        return balanced_incomplete_block_design(v,7)
    else:
        from sage.combinat.designs.orthogonal_arrays import orthogonal_array
        # IX.4.5.m from [BJL99].
        #
        # This construction takes a TD(5,g) and truncates its last column to
        # size u: it yields a (4g+u,{4,5},{g,u})-GDD. If there exists a
        # (3g+1,{4,7})-PBD and a (3u+1,{4,7})-PBD, then we can apply the x->3x+1
        # construction on the truncated transversal design (which is a GDD).
        #
        # We write vv = 4g+u while satisfying the hypotheses.
        vv = (v - 1) // 3
        for g in range((vv + 5 - 1) // 5, vv // 4 + 1):
            u = vv-4*g
            if (orthogonal_array(5,g,existence=True) is True and
                PBD_4_7(3*g+1,existence=True) is True and
                PBD_4_7(3*u+1,existence=True) is True):
                from .orthogonal_arrays import transversal_design
                domain = set(range(vv))
                GDD = transversal_design(5,g)
                GDD = GroupDivisibleDesign(vv,
                                           groups=[[x for x in gr if x in domain] for gr in GDD.groups()],
                                           blocks=[[x for x in B if x in domain] for B in GDD],
                                           G=set([g,u]),
                                           K=[4,5],
                                           check=False)
                return PBD_4_7_from_Y(GDD,check=check)

    return PairwiseBalancedDesign(v,
                                  blocks=blocks,
                                  K=[4,7],
                                  check=check,
                                  copy=False)


def PBD_4_7_from_Y(gdd,check=True):
    r"""
    Return a `(3v+1,\{4,7\})`-PBD from a `(v,\{4,5,7\},\NN-\{3,6,10\})`-GDD.

    This implements Lemma IX.3.11 from [BJL99]_ (p.625). All points of the GDD
    are tripled, and a `+\infty` point is added to the design.

    - A group of size `s\in Y=\NN-\{3,6,10\}` becomes a set of size `3s`. Adding
      `\infty` to it gives it size `3s+1`, and this set is then replaced by a
      `(3s+1,\{4,7\})`-PBD.

    - A block of size `s\in\{4,5,7\}` becomes a `(3s,\{4,7\},\{3\})`-GDD.

    This lemma is part of the existence proof of `(v,\{4,7\})`-PBD as explained
    in IX.4.5 from [BJL99]_).

    INPUT:

    - ``gdd`` -- a `(v,\{4,5,7\},Y)`-GDD where `Y=\NN-\{3,6,10\}`

    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.

    EXAMPLES::

        sage: from sage.combinat.designs.resolvable_bibd import PBD_4_7_from_Y
        sage: PBD_4_7_from_Y(designs.transversal_design(7,8))                           # needs sage.schemes
        Pairwise Balanced Design on 169 points with sets of sizes in [4, 7]

    TESTS::

        sage: PBD_4_7_from_Y(designs.balanced_incomplete_block_design(10,10))           # needs sage.schemes
        Traceback (most recent call last):
        ...
        ValueError: The GDD should only contain blocks of size {4,5,7} but there are other: [10]
        sage: PBD_4_7_from_Y(designs.transversal_design(4,3))                           # needs sage.schemes
        Traceback (most recent call last):
        ...
        RuntimeError: A group has size 3 but I do not know how to build a (10,[4,7])-PBD
    """
    from .group_divisible_designs import group_divisible_design
    from .bibd import PairwiseBalancedDesign
    block_sizes = set(map(len,gdd._blocks))
    group_sizes = set(map(len,gdd._groups))
    if not block_sizes.issubset([4, 5, 7]):
        txt = list(block_sizes.difference([4, 5, 7]))
        raise ValueError("The GDD should only contain blocks of size {{4,5,7}} "
                         "but there are other: {}".format(txt))

    for gs in group_sizes:
        if PBD_4_7(3*gs+1,existence=True) is not True:
            raise RuntimeError("A group has size {} but I do not know how to "
                               "build a ({},[4,7])-PBD".format(gs,3*gs+1))

    GDD = {}  # the GDD we will need
    if 4 in block_sizes:
        # GDD[4] = GDD_from_BIBD(3*4,4)
        GDD[4] = group_divisible_design(3*4,K=[4],G=[3])
    if 5 in block_sizes:
        # GDD[5] = GDD_from_BIBD(3*5,4)
        GDD[5] = group_divisible_design(3*5,K=[4],G=[3])
    if 7 in block_sizes:
        # It is obtained from a PBD_4_7(22) by removing a point only contained
        # in sets of size 4
        GDD[7] = PBD_4_7(22)
        x = set(range(22)).difference(*[S for S in GDD[7] if len(S) != 4]).pop()
        relabel = sum((S for S in GDD[7] if x in S),[])  # the groups must be 012,345,...
        relabel = [xx for xx in relabel if xx != x]+[x]
        GDD[7].relabel({v:i for i,v in enumerate(relabel)})
        GDD[7] = [S for S in GDD[7] if 21 not in S]

    PBD = []

    # The blocks
    for B in gdd:
        for B_GDD in GDD[len(B)]:
            PBD.append([3*B[x//3]+(x % 3) for x in B_GDD])

    # The groups
    group_PBD = {gs:PBD_4_7(3*gs+1) for gs in group_sizes}
    for G in gdd.groups():
        gs = len(G)
        for B in group_PBD[gs]:
            PBD.append([3*G[x//3]+(x % 3) if x < 3*gs else 3*gdd.num_points()
                        for x in B])

    return PairwiseBalancedDesign(3*gdd.num_points()+1,
                                  blocks=PBD,
                                  K=[4,7],
                                  check=check,
                                  copy=False)
