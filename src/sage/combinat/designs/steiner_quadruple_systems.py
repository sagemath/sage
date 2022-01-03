r"""
Steiner Quadruple Systems

A Steiner Quadruple System on `n` points is a family `SQS_n \subset \binom {[n]}
4` of `4`-sets, such that any set `S\subset [n]` of size three is a subset of
exactly one member of `SQS_n`.

This module implements Haim Hanani's constructive proof that a Steiner Quadruple
System exists if and only if `n\equiv 2,4 \pmod 6`. Hanani's proof consists in 6
different constructions that build a large Steiner Quadruple System from a smaller
one, and though it does not give a very clear understanding of why it works (to say the
least)... it does !

The constructions have been implemented while reading two papers simultaneously,
for one of them sometimes provides the informations that the other one does
not. The first one is Haim Hanani's original paper [Han1960]_, and the other
one is a paper from Horan and Hurlbert which goes through all constructions
[HH2012]_.

It can be used through the ``designs`` object::

    sage: designs.steiner_quadruple_system(8)
    Incidence structure with 8 points and 14 blocks

AUTHORS:

- Nathann Cohen (May 2013, while listening to "*Le Blues Du Pauvre Delahaye*")

Index
-----

This module's main function is the following:

.. csv-table::
    :class: contentstable
    :widths: 15, 20, 65
    :delim: |

    | :func:`steiner_quadruple_system` | Return a Steiner Quadruple System on `n` points

This function redistributes its work among 6 constructions:

.. csv-table::
    :class: contentstable
    :widths: 15, 20, 65
    :delim: |

    Construction `1` | :func:`two_n`                | Return a Steiner Quadruple System on `2n` points
    Construction `2` | :func:`three_n_minus_two`    | Return a Steiner Quadruple System on `3n-2` points
    Construction `3` | :func:`three_n_minus_eight`  | Return a Steiner Quadruple System on `3n-8` points
    Construction `4` | :func:`three_n_minus_four`   | Return a Steiner Quadruple System on `3n-4` points
    Construction `5` | :func:`four_n_minus_six`     | Return a Steiner Quadruple System on `4n-6` points
    Construction `6` | :func:`twelve_n_minus_ten`   | Return a Steiner Quadruple System on `12n-10` points

It also defines two specific Steiner Quadruple Systems that the constructions
require, i.e. `SQS_{14}` and `SQS_{38}` as well as the systems of pairs
`P_{\alpha}(m)` and `\overline P_{\alpha}(m)` (see [Han1960]_).

Functions
---------
"""

from sage.misc.cachefunc import cached_function
from sage.combinat.designs.incidence_structures import IncidenceStructure

# Construction 1
def two_n(B):
    r"""
    Return a Steiner Quadruple System on `2n` points.

    INPUT:

    - ``B`` -- A Steiner Quadruple System on `n` points.

    EXAMPLES::

        sage: from sage.combinat.designs.steiner_quadruple_systems import two_n
        sage: for n in range(4, 30):
        ....:     if (n%6) in [2,4]:
        ....:         sqs = designs.steiner_quadruple_system(n)
        ....:         if not two_n(sqs).is_t_design(3,2*n,4,1):
        ....:             print("Something is wrong !")

    """
    n = B.num_points()
    Y = []

    # Line 1
    for x,y,z,t in B._blocks:
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    d = (a+b+c)%2
                    Y.append([x+a*n,y+b*n,z+c*n,t+d*n])

    # Line 2
    for j in range(n):
        for jj in range(j+1,n):
            Y.append([j,jj,n+j,n+jj])

    return IncidenceStructure(2*n,Y,check=False,copy=False)

# Construction 2
def three_n_minus_two(B):
    """
    Return a Steiner Quadruple System on `3n-2` points.

    INPUT:

    - ``B`` -- A Steiner Quadruple System on `n` points.

    EXAMPLES::

        sage: from sage.combinat.designs.steiner_quadruple_systems import three_n_minus_two
        sage: for n in range(4, 30):
        ....:     if (n%6) in [2,4]:
        ....:         sqs = designs.steiner_quadruple_system(n)
        ....:         if not three_n_minus_two(sqs).is_t_design(3,3*n-2,4,1):
        ....:             print("Something is wrong !")
    """
    n = B.num_points()
    A = n-1
    Y = []
    # relabel function
    r = lambda i,x : (i%3)*(n-1)+x
    for x,y,z,t in B._blocks:
        if t == A:
            # Line 2.
            for a in range(3):
                for b in range(3):
                    c = -(a+b)%3
                    Y.append([r(a,x),r(b,y),r(c,z),3*n-3])

            # Line 3.
            Y.extend([[r(i,x),r(i,y),r(i+1,z),r(i+2,z)] for i in range(3)])
            Y.extend([[r(i,x),r(i,z),r(i+1,y),r(i+2,y)] for i in range(3)])
            Y.extend([[r(i,y),r(i,z),r(i+1,x),r(i+2,x)] for i in range(3)])

        else:
            # Line 1.
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        d = -(a+b+c)%3
                        Y.append([r(a,x),r(b,y),r(c,z),r(d,t)])

    # Line 4.
    for j in range(n-1):
        for jj in range(j+1,n-1):
            Y.extend([[r(i,j),r(i,jj),r(i+1,j),r(i+1,jj)] for i in range(3)])

    # Line 5.
    for j in range(n-1):
        Y.append([r(0,j),r(1,j),r(2,j),3*n-3])

    return IncidenceStructure(3*n-2,Y,check=False,copy=False)

# Construction 3
def three_n_minus_eight(B):
    r"""
    Return a Steiner Quadruple System on `3n-8` points.

    INPUT:

    - ``B`` -- A Steiner Quadruple System on `n` points.

    EXAMPLES::

        sage: from sage.combinat.designs.steiner_quadruple_systems import three_n_minus_eight
        sage: for n in range(4, 30):
        ....:     if (n%12) == 2:
        ....:         sqs = designs.steiner_quadruple_system(n)
        ....:         if not three_n_minus_eight(sqs).is_t_design(3,3*n-8,4,1):
        ....:             print("Something is wrong !")

    """
    n = B.num_points()

    if (n%12) != 2:
        raise ValueError("n must be equal to 2 mod 12")

    B = relabel_system(B)
    r = lambda i,x : (i%3)*(n-4)+(x%(n-4))

    # Line 1.
    Y = [[x+2*(n-4) for x in B._blocks[-1]]]

    # Line 2.
    for s in B._blocks[:-1]:
        for i in range(3):
            Y.append([r(i,x) if x<= n-5 else x+2*(n-4) for x in s])


    # Line 3.
    for a in range(4):
        for aa in range(n-4):
            for aaa in range(n-4):
                aaaa = -(a+aa+aaa)%(n-4)
                Y.append([r(0,aa),r(1,aaa), r(2,aaaa),3*(n-4)+a])


    # Line 4.
    k = (n-14) // 12
    for i in range(3):
        for b in range(n-4):
            for bb in range(n-4):
                bbb = -(b+bb)%(n-4)
                for d in range(2*k+1):
                    Y.append([r(i+2,bbb), r(i, b+2*k+1+i*(4*k+2)-d) , r(i, b+2*k+2+i*(4*k+2)+d), r(i+1,bb)])

    # Line 5.
    for i in range(3):
        for alpha in range(4*k+2, 12*k+9):
            for ra,sa in P(alpha,6*k+5):
                for raa,saa in P(alpha,6*k+5):
                    Y.append([r(i,ra),r(i,sa),r(i+1,raa), r(i+1,saa)])

    return IncidenceStructure(3*n-8,Y,check=False,copy=False)

# Construction 4
def three_n_minus_four(B):
    r"""
    Return a Steiner Quadruple System on `3n-4` points.

    INPUT:

    - ``B`` -- A Steiner Quadruple System on `n` points where `n\equiv
      10\pmod{12}`.

    EXAMPLES::

        sage: from sage.combinat.designs.steiner_quadruple_systems import three_n_minus_four
        sage: for n in range(4, 30):
        ....:     if n%12 == 10:
        ....:         sqs = designs.steiner_quadruple_system(n)
        ....:         if not three_n_minus_four(sqs).is_t_design(3,3*n-4,4,1):
        ....:             print("Something is wrong !")

    """
    n = B.num_points()

    if n%12 != 10:
        raise ValueError("n must be equal to 10 mod 12")

    B = relabel_system(B)
    r = lambda i,x : (i%3)*(n-2)+(x%(n-2))

    # Line 1/2.
    Y = []
    for s in B._blocks:
        for i in range(3):
            Y.append([r(i,x) if x<= n-3 else x+2*(n-2) for x in s])

    # Line 3.
    for a in range(2):
        for aa in range(n-2):
            for aaa in range(n-2):
                aaaa = -(a+aa+aaa) % (n-2)
                Y.append([r(0,aa),r(1,aaa), r(2,aaaa),3*(n-2)+a])

    # Line 4.
    k = (n-10) // 12
    for i in range(3):
        for b in range(n-2):
            for bb in range(n-2):
                bbb = -(b+bb)%(n-2)
                for d in range(2*k+1):
                    Y.append([r(i+2,bbb), r(i, b+2*k+1+i*(4*k+2)-d) , r(i, b+2*k+2+i*(4*k+2)+d), r(i+1,bb)])

    # Line 5.
    from sage.graphs.graph_coloring import round_robin
    one_factorization = round_robin(2*(6*k+4)).edges()
    color_classes = [[] for j in range(2*(6*k+4)-1)]
    for u,v,l in one_factorization:
        color_classes[l].append((u,v))

    for i in range(3):
        for alpha in range(4*k+2, 12*k+6+1):
            for ra,sa in P(alpha, 6*k+4):
                for raa,saa in P(alpha, 6*k+4):
                    Y.append([r(i,ra),r(i,sa),r(i+1,raa), r(i+1,saa)])

    return IncidenceStructure(3*n-4,Y,check=False,copy=False)

# Construction 5
def four_n_minus_six(B):
    """
    Return a Steiner Quadruple System on `4n-6` points.

    INPUT:

    - ``B`` -- A Steiner Quadruple System on `n` points.

    EXAMPLES::

        sage: from sage.combinat.designs.steiner_quadruple_systems import four_n_minus_six
        sage: for n in range(4, 20):
        ....:     if (n%6) in [2,4]:
        ....:         sqs = designs.steiner_quadruple_system(n)
        ....:         if not four_n_minus_six(sqs).is_t_design(3,4*n-6,4,1):
        ....:             print("Something is wrong !")

    """
    n = B.num_points()
    f = n-2
    r = lambda i,ii,x : (2*(i%2)+(ii%2))*(n-2)+(x)%(n-2)

    # Line 1.
    Y = []
    for s in B._blocks:
        for i in range(2):
            for ii in range(2):
                Y.append([r(i,ii,x) if x<= n-3 else x+3*(n-2) for x in s])

    # Line 2/3/4/5
    k = f // 2
    for l in range(2):
        for eps in range(2):
            for c in range(k):
                for cc in range(k):
                    ccc = -(c+cc)%k
                    Y.append([4*(n-2)+l, r(0,0,2*c)  , r(0,1,2*cc-eps)  , r(1,eps,2*ccc+l)  ])
                    Y.append([4*(n-2)+l, r(0,0,2*c+1), r(0,1,2*cc-1-eps), r(1,eps,2*ccc+1-l)])
                    Y.append([4*(n-2)+l, r(1,0,2*c)  , r(1,1,2*cc-eps)  , r(0,eps,2*ccc+1-l)])
                    Y.append([4*(n-2)+l, r(1,0,2*c+1), r(1,1,2*cc-1-eps), r(0,eps,2*ccc+l)  ])

    # Line 6/7
    for h in range(2):
        for eps in range(2):
            for ccc in range(k):
                assert len(barP(ccc,k)) == k-1
                for rc,sc in barP(ccc,k):
                    for c in range(k):
                        cc = -(c+ccc)%k
                        Y.append([r(h,0,2*c+eps)  , r(h,1,2*cc-eps), r(h+1,0,rc), r(h+1,0,sc)])
                        Y.append([r(h,0,2*c-1+eps), r(h,1,2*cc-eps), r(h+1,1,rc), r(h+1,1,sc)])



    # Line 8/9
    for h in range(2):
        for eps in range(2):
            for ccc in range(k):
                for rc,sc in barP(k+ccc,k):
                    for c in range(k):
                        cc = -(c+ccc)%k
                        Y.append([r(h,0,2*c+eps)  , r(h,1,2*cc-eps), r(h+1,1,rc), r(h+1,1,sc)])
                        Y.append([r(h,0,2*c-1+eps), r(h,1,2*cc-eps), r(h+1,0,rc), r(h+1,0,sc)])


    # Line 10
    for h in range(2):
        for alpha in range(n-3):
            for ra,sa in P(alpha,k):
                for raa,saa in P(alpha,k):
                    Y.append([r(h,0,ra),r(h,0,sa),r(h,1,raa),r(h,1,saa)])

    return IncidenceStructure(4*n-6,Y,check=False,copy=False)

# Construction 6
def twelve_n_minus_ten(B):
    """
    Return a Steiner Quadruple System on `12n-6` points.

    INPUT:

    - ``B`` -- A Steiner Quadruple System on `n` points.

    EXAMPLES::

        sage: from sage.combinat.designs.steiner_quadruple_systems import twelve_n_minus_ten
        sage: for n in range(4, 15):
        ....:     if (n%6) in [2,4]:
        ....:         sqs = designs.steiner_quadruple_system(n)
        ....:         if not twelve_n_minus_ten(sqs).is_t_design(3,12*n-10,4,1):
        ....:             print("Something is wrong !")

    """
    n = B.num_points()
    B14 = steiner_quadruple_system(14)
    r = lambda i,x : i%(n-1)+(x%12)*(n-1)

    # Line 1.
    Y = []
    for s in B14._blocks:
        for i in range(n-1):
            Y.append([r(i,x) if x<= 11 else r(n-2,11)+x-11 for x in s])

    for s in B._blocks:
        if s[-1] == n-1:
            u,v,w,B = s
            dd = {0:u,1:v,2:w}
            d = lambda x:dd[x%3]
            for b in range(12):
                for bb in range(12):
                    bbb = -(b+bb)%12
                    for h in range(2):
                        # Line 2
                        Y.append([r(n-2,11)+1+h,r(u,b),r(v,bb),r(w,bbb+3*h)])

                    for i in range(3):
                        # Line 38.3
                        Y.append([r(d(i),b+4+i), r(d(i),b+7+i), r(d(i+1),bb), r(d(i+2),bbb)])

            for j in range(12):
                for eps in range(2):
                    for i in range(3):
                        # Line 38.4-38.7
                        Y.append([ r(d(i),j), r(d(i+1),j+6*eps  ), r(d(i+2),6*eps-2*j+1), r(d(i+2),6*eps-2*j-1)])
                        Y.append([ r(d(i),j), r(d(i+1),j+6*eps  ), r(d(i+2),6*eps-2*j+2), r(d(i+2),6*eps-2*j-2)])
                        Y.append([ r(d(i),j), r(d(i+1),j+6*eps-3), r(d(i+2),6*eps-2*j+1), r(d(i+2),6*eps-2*j+2)])
                        Y.append([ r(d(i),j), r(d(i+1),j+6*eps+3), r(d(i+2),6*eps-2*j-1), r(d(i+2),6*eps-2*j-2)])

            for j in range(6):
                for i in range(3):
                    for eps in range(2):
                        # Line 38.8
                        Y.append([ r(d(i),j), r(d(i),j+6), r(d(i+1),j+3*eps), r(d(i+1),j+6+3*eps)])

            for j in range(12):
                for i in range(3):
                    for eps in range(4):
                        # Line 38.11
                        Y.append([ r(d(i),j), r(d(i),j+1), r(d(i+1),j+3*eps), r(d(i+1),j+3*eps+1)])
                        # Line 38.12
                        Y.append([ r(d(i),j), r(d(i),j+2), r(d(i+1),j+3*eps), r(d(i+1),j+3*eps+2)])
                        # Line 38.13
                        Y.append([ r(d(i),j), r(d(i),j+4), r(d(i+1),j+3*eps), r(d(i+1),j+3*eps+4)])

            for alpha in [4,5]:
                for ra,sa in P(alpha,6):
                    for raa,saa in P(alpha,6):
                        for i in range(3):
                            for ii in range(i+1,3):
                                # Line 38.14
                                Y.append([ r(d(i),ra), r(d(i),sa), r(d(ii),raa), r(d(ii),saa)])

            for g in range(6):
                for eps in range(2):
                    for i in range(3):
                        for ii in range(3):
                            if i == ii:
                                continue
                            # Line 38.9
                            Y.append([ r(d(i),2*g+3*eps), r(d(i),2*g+6+3*eps), r(d(ii),2*g+1), r(d(ii),2*g+5)])
                            # Line 38.10
                            Y.append([ r(d(i),2*g+3*eps), r(d(i),2*g+6+3*eps), r(d(ii),2*g+2), r(d(ii),2*g+4)])

        else:
            x,y,z,t = s
            for a in range(12):
                for aa in range(12):
                    for aaa in range(12):
                        aaaa = -(a+aa+aaa)%12
                        # Line 3
                        Y.append([r(x,a), r(y,aa), r(z,aaa), r(t,aaaa)])
    return IncidenceStructure(12*n-10,Y,check=False,copy=False)

def relabel_system(B):
    r"""
    Relabels the set so that `\{n-4, n-3, n-2, n-1\}` is in `B`.

    INPUT:

    - ``B`` -- a list of 4-uples on `0,...,n-1`.

    EXAMPLES::

        sage: from sage.combinat.designs.steiner_quadruple_systems import relabel_system
        sage: SQS8 = designs.steiner_quadruple_system(8)
        sage: relabel_system(SQS8)
        Incidence structure with 8 points and 14 blocks
    """
    n = B.num_points()
    B0 = B._blocks[0]

    label = {
        B0[0] : n-4,
        B0[1] : n-3,
        B0[2] : n-2,
        B0[3] : n-1
        }

    def get_label(x):
        if x in label:
            return label[x]
        else:
            total = len(label)-4
            label[x] = total
            return total

    B = [[get_label(_) for _ in s] for s in B]
    return IncidenceStructure(n,B)

def P(alpha, m):
    r"""
    Return the collection of pairs `P_{\alpha}(m)`

    For more information on this system, see [Han1960]_.

    EXAMPLES::

        sage: from sage.combinat.designs.steiner_quadruple_systems import P
        sage: P(3,4)
        [(0, 5), (2, 7), (4, 1), (6, 3)]
    """
    if alpha >= 2*m-1:
        raise Exception
    if m%2==0:
        if alpha < m:
            if alpha%2 == 0:
                b = alpha // 2
                return [(2*a, (2*a + 2*b + 1)%(2*m)) for a in range(m)]
            else:
                b = (alpha-1) // 2
                return [(2*a, (2*a - 2*b - 1)%(2*m)) for a in range(m)]
        else:
            y = alpha - m
            pairs  = [(b,(2*y-b)%(2*m)) for b in range(y)]
            pairs += [(c,(2*m+2*y-c-2)%(2*m)) for c in range(2*y+1,m+y-1)]
            pairs += [(2*m+int(-1.5-.5*(-1)**y),y),(2*m+int(-1.5+.5*(-1)**y),m+y-1)]
            return pairs
    else:
        if alpha < m-1:
            if alpha % 2 == 0:
                b = alpha // 2
                return [(2*a,(2*a+2*b+1)%(2*m)) for a in range(m)]
            else:
                b = (alpha-1) // 2
                return [(2*a,(2*a-2*b-1)%(2*m)) for a in range(m)]
        else:
            y = alpha-m+1
            pairs  = [(b,2*y-b) for b in range(y)]
            pairs += [(c,2*m+2*y-c) for c in range(2*y+1,m+y)]
            pairs += [(y,m+y)]
            return pairs

def _missing_pair(n,l):
    r"""
    Return the smallest `(x,x+1)` that is not contained in `l`.

    EXAMPLES::

        sage: from sage.combinat.designs.steiner_quadruple_systems import _missing_pair
        sage: _missing_pair(6, [(0,1), (4,5)])
        (2, 3)
    """
    l = set(x for X in l for x in X)
    for x in range(n):
        if x not in l:
            break

    assert x not in l
    assert x + 1 not in l
    return (x, x + 1)


def barP(eps, m):
    r"""
    Return the collection of pairs `\overline P_{\alpha}(m)`

    For more information on this system, see [Han1960]_.

    EXAMPLES::

        sage: from sage.combinat.designs.steiner_quadruple_systems import barP
        sage: barP(3,4)
        [(0, 4), (3, 5), (1, 2)]
    """
    return barP_system(m)[eps]

@cached_function
def barP_system(m):
    r"""
    Return the 1-factorization of `K_{2m}` `\overline P(m)`

    For more information on this system, see [Han1960]_.

    EXAMPLES::

        sage: from sage.combinat.designs.steiner_quadruple_systems import barP_system
        sage: barP_system(3)
        [[(4, 3), (2, 5)],
        [(0, 5), (4, 1)],
        [(0, 2), (1, 3)],
        [(1, 5), (4, 2), (0, 3)],
        [(0, 4), (3, 5), (1, 2)],
        [(0, 1), (2, 3), (4, 5)]]
    """
    isequal = lambda e1,e2 : e1 == e2 or e1 == tuple(reversed(e2))
    pairs = []
    last = []

    if m % 2 == 0:
        # The first (shorter) collections of  pairs, obtained from P by removing
        # pairs. Those are added to 'last', a new list of pairs
        last = []
        for n in range(1, (m-2)//2+1):
            pairs.append([p for p in P(2*n,m) if not isequal(p,(2*n,(4*n+1)%(2*m)))])
            last.append((2*n,(4*n+1)%(2*m)))
            pairs.append([p for p in P(2*n-1,m) if not isequal(p,(2*m-2-2*n,2*m-1-4*n))])
            last.append((2*m-2-2*n,2*m-1-4*n))

        pairs.append([p for p in P(m,m) if not isequal(p,(2*m-2,0))])
        last.append((2*m-2,0))
        pairs.append([p for p in P(m+1,m) if not isequal(p,(2*m-1,1))])
        last.append((2*m-1,1))

        assert all(len(pp) == m-1 for pp in pairs)
        assert len(last) == m

        # Pairs of normal length

        pairs.append(P(0,m))
        pairs.append(P(m-1,m))

        for alpha in range(m+2,2*m-1):
            pairs.append(P(alpha,m))
        pairs.append(last)

        assert len(pairs) == 2*m

        # Now the points must be relabeled
        relabel = {}
        for n in range(1, (m-2)//2+1):
            relabel[2*n] = (4*n)%(2*m)
            relabel[4*n+1] = (4*n+1)%(2*m)
            relabel[2*m-2-2*n] = (4*n-2)%(2*m)
            relabel[2*m-1-4*n] = (4*n-1)%(2*m)

        relabel[2*m-2] = (1)%(2*m)
        relabel[0] = 0
        relabel[2*m-1] = 2*m-1
        relabel[1] = 2*m-2

    else:
        # The first (shorter) collections of  pairs, obtained from P by removing
        # pairs. Those are added to 'last', a new list of pairs

        last = []
        for n in range((m - 3) // 2 + 1):
            pairs.append([p for p in P(2*n,m) if not isequal(p,(2*n,(4*n+1)%(2*m)))])
            last.append((2*n,(4*n+1)%(2*m)))
            pairs.append([p for p in P(2*n+1,m) if not isequal(p,(2*m-2-2*n,2*m-3-4*n))])
            last.append((2*m-2-2*n,2*m-3-4*n))

        pairs.append([p for p in P(2*m-2,m) if not isequal(p,(m-1,2*m-1))])
        last.append((m-1,2*m-1))

        assert all(len(pp) == m-1 for pp in pairs)
        assert len(pairs) == m

        # Pairs of normal length

        for alpha in range(m-1,2*m-2):
            pairs.append(P(alpha,m))
        pairs.append(last)

        assert len(pairs) == 2*m

        # Now the points must be relabeled
        relabel = {}
        for n in range((m - 3) // 2 + 1):
            relabel[2*n] = (4*n)%(2*m)
            relabel[4*n+1] = (4*n+1)%(2*m)
            relabel[2*m-2-2*n] = (4*n+2)%(2*m)
            relabel[2*m-3-4*n] = (4*n+3)%(2*m)
        relabel[m-1] = (2*m-2)%(2*m)
        relabel[2*m-1] = 2*m-1

    assert len(relabel) == 2*m
    assert len(pairs) == 2*m

    # Relabeling the points

    pairs = [[(relabel[x],relabel[y]) for x,y in pp] for pp in pairs]

    # Pairs are sorted first according to their cardinality, then using the
    # number of the smallest point that they do NOT contain.
    pairs.sort(key=lambda x: _missing_pair(2*m+1,x))

    return pairs

@cached_function
def steiner_quadruple_system(n, check = False):
    r"""
    Return a Steiner Quadruple System on `n` points.

    INPUT:

    - ``n`` -- an integer such that `n\equiv 2,4\pmod 6`

    - ``check`` (boolean) -- whether to check that the system is a Steiner
      Quadruple System before returning it (`False` by default)

    EXAMPLES::

        sage: sqs4 = designs.steiner_quadruple_system(4)
        sage: sqs4
        Incidence structure with 4 points and 1 blocks
        sage: sqs4.is_t_design(3,4,4,1)
        True

        sage: sqs8 = designs.steiner_quadruple_system(8)
        sage: sqs8
        Incidence structure with 8 points and 14 blocks
        sage: sqs8.is_t_design(3,8,4,1)
        True

    TESTS::

        sage: for n in range(4, 100):                                      # long time
        ....:     if (n%6) in [2,4]:                                        # long time
        ....:         sqs = designs.steiner_quadruple_system(n, check=True) # long time
    """
    n = int(n)
    if not ((n%6) in [2, 4]):
        raise ValueError("n mod 6 must be equal to 2 or 4")
    elif n == 4:
        sqs = IncidenceStructure(4, [[0,1,2,3]], copy = False, check = False)
    elif n == 14:
        sqs = IncidenceStructure(14, _SQS14(), copy = False, check = False)
    elif n == 38:
        sqs = IncidenceStructure(38, _SQS38(), copy = False, check = False)
    elif n%12 in [4, 8]:
        nn =  n // 2
        sqs = two_n(steiner_quadruple_system(nn, check = False))
    elif n%18 in [4,10]:
        nn = (n+2) // 3
        sqs = three_n_minus_two(steiner_quadruple_system(nn, check = False))
    elif (n%36) == 34:
        nn = (n+8) // 3
        sqs = three_n_minus_eight(steiner_quadruple_system(nn, check = False))
    elif (n%36) == 26:
        nn = (n+4) // 3
        sqs = three_n_minus_four(steiner_quadruple_system(nn, check = False))
    elif n%24 in [2, 10]:
        nn = (n+6) // 4
        sqs = four_n_minus_six(steiner_quadruple_system(nn, check = False))
    elif n%72 in [14, 38]:
        nn = (n+10) // 12
        sqs = twelve_n_minus_ten(steiner_quadruple_system(nn, check = False))
    else:
        raise ValueError("This shouldn't happen !")

    if check and not sqs.is_t_design(3,n,4,1):
        raise RuntimeError("Something is very very wrong.")

    return sqs

def _SQS14():
    r"""
    Return a Steiner Quadruple System on 14 points.

    Obtained from the La Jolla Covering Repository.

    EXAMPLES::

        sage: from sage.combinat.designs.steiner_quadruple_systems import _SQS14
        sage: sqs14 = IncidenceStructure(_SQS14())
        sage: sqs14.is_t_design(3,14,4,1)
        True
    """
    return [[0, 1, 2, 5], [0, 1, 3, 6], [0, 1, 4, 13], [0, 1, 7, 10], [0, 1, 8, 9],
            [0, 1, 11, 12], [0, 2, 3, 4], [0, 2, 6, 12], [0, 2, 7, 9], [0, 2, 8, 11],
            [0, 2, 10, 13], [0, 3, 5, 13], [0, 3, 7, 11], [0, 3, 8, 10], [0, 3, 9, 12],
            [0, 4, 5, 9], [0, 4, 6, 11], [0, 4, 7, 8], [0, 4, 10, 12], [0, 5, 6, 8],
            [0, 5, 7, 12], [0, 5, 10, 11], [0, 6, 7, 13], [0, 6, 9, 10], [0, 8, 12, 13],
            [0, 9, 11, 13], [1, 2, 3, 13], [1, 2, 4, 12], [1, 2, 6, 9], [1, 2, 7, 11],
            [1, 2, 8, 10], [1, 3, 4, 5], [1, 3, 7, 8], [1, 3, 9, 11], [1, 3, 10, 12],
            [1, 4, 6, 10], [1, 4, 7, 9], [1, 4, 8, 11], [1, 5, 6, 11], [1, 5, 7, 13],
            [1, 5, 8, 12], [1, 5, 9, 10], [1, 6, 7, 12], [1, 6, 8, 13], [1, 9, 12, 13],
            [1, 10, 11, 13], [2, 3, 5, 11], [2, 3, 6, 7], [2, 3, 8, 12], [2, 3, 9, 10],
            [2, 4, 5, 13], [2, 4, 6, 8], [2, 4, 7, 10], [2, 4, 9, 11], [2, 5, 6, 10],
            [2, 5, 7, 8], [2, 5, 9, 12], [2, 6, 11, 13], [2, 7, 12, 13], [2, 8, 9, 13],
            [2, 10, 11, 12], [3, 4, 6, 9], [3, 4, 7, 12], [3, 4, 8, 13], [3, 4, 10, 11],
            [3, 5, 6, 12], [3, 5, 7, 10], [3, 5, 8, 9], [3, 6, 8, 11], [3, 6, 10, 13],
            [3, 7, 9, 13], [3, 11, 12, 13], [4, 5, 6, 7], [4, 5, 8, 10], [4, 5, 11, 12],
            [4, 6, 12, 13], [4, 7, 11, 13], [4, 8, 9, 12], [4, 9, 10, 13], [5, 6, 9, 13],
            [5, 7, 9, 11], [5, 8, 11, 13], [5, 10, 12, 13], [6, 7, 8, 9], [6, 7, 10, 11],
            [6, 8, 10, 12], [6, 9, 11, 12], [7, 8, 10, 13], [7, 8, 11, 12], [7, 9, 10, 12],
            [8, 9, 10, 11]]

def _SQS38():
    r"""
    Return a Steiner Quadruple System on 14 points.

    Obtained from the La Jolla Covering Repository.

    EXAMPLES::

        sage: from sage.combinat.designs.steiner_quadruple_systems import _SQS38
        sage: sqs38 = IncidenceStructure(_SQS38())
        sage: sqs38.is_t_design(3,38,4,1)
        True
    """
    # From the La Jolla Covering Repository
    return [[0, 1, 2, 14], [0, 1, 3, 34], [0, 1, 4, 31], [0, 1, 5, 27], [0, 1, 6, 17],
            [0, 1, 7, 12], [0, 1, 8, 36], [0, 1, 9, 10], [0, 1, 11, 18], [0, 1, 13, 37],
            [0, 1, 15, 35], [0, 1, 16, 22], [0, 1, 19, 33], [0, 1, 20, 25], [0, 1, 21, 23],
            [0, 1, 24, 32], [0, 1, 26, 28], [0, 1, 29, 30], [0, 2, 3, 10], [0, 2, 4, 9],
            [0, 2, 5, 28], [0, 2, 6, 15], [0, 2, 7, 36], [0, 2, 8, 23], [0, 2, 11, 22],
            [0, 2, 12, 13], [0, 2, 16, 25], [0, 2, 17, 18], [0, 2, 19, 30], [0, 2, 20, 35],
            [0, 2, 21, 29], [0, 2, 24, 34], [0, 2, 26, 31], [0, 2, 27, 32], [0, 2, 33, 37],
            [0, 3, 4, 18], [0, 3, 5, 23], [0, 3, 6, 32], [0, 3, 7, 19], [0, 3, 8, 20],
            [0, 3, 9, 17], [0, 3, 11, 25], [0, 3, 12, 24], [0, 3, 13, 27], [0, 3, 14, 31],
            [0, 3, 15, 22], [0, 3, 16, 28], [0, 3, 21, 33], [0, 3, 26, 36], [0, 3, 29, 35],
            [0, 3, 30, 37], [0, 4, 5, 7], [0, 4, 6, 28], [0, 4, 8, 25], [0, 4, 10, 30],
            [0, 4, 11, 20], [0, 4, 12, 32], [0, 4, 13, 36], [0, 4, 14, 29], [0, 4, 15, 27],
            [0, 4, 16, 35], [0, 4, 17, 22], [0, 4, 19, 23], [0, 4, 21, 34], [0, 4, 24, 33],
            [0, 4, 26, 37], [0, 5, 6, 24], [0, 5, 8, 26], [0, 5, 9, 29], [0, 5, 10, 20],
            [0, 5, 11, 13], [0, 5, 12, 14], [0, 5, 15, 33], [0, 5, 16, 37], [0, 5, 17, 35],
            [0, 5, 18, 19], [0, 5, 21, 25], [0, 5, 22, 30], [0, 5, 31, 32], [0, 5, 34, 36],
            [0, 6, 7, 30], [0, 6, 8, 33], [0, 6, 9, 12], [0, 6, 10, 18], [0, 6, 11, 37],
            [0, 6, 13, 31], [0, 6, 14, 35], [0, 6, 16, 29], [0, 6, 19, 25], [0, 6, 20, 27],
            [0, 6, 21, 36], [0, 6, 22, 23], [0, 6, 26, 34], [0, 7, 8, 11], [0, 7, 9, 33],
            [0, 7, 10, 21], [0, 7, 13, 20], [0, 7, 14, 22], [0, 7, 15, 31], [0, 7, 16, 34],
            [0, 7, 17, 29], [0, 7, 18, 24], [0, 7, 23, 26], [0, 7, 25, 32], [0, 7, 27, 28],
            [0, 7, 35, 37], [0, 8, 9, 37], [0, 8, 10, 27], [0, 8, 12, 18], [0, 8, 13, 30],
            [0, 8, 14, 15], [0, 8, 16, 21], [0, 8, 17, 19], [0, 8, 22, 35], [0, 8, 24, 31],
            [0, 8, 28, 34], [0, 8, 29, 32], [0, 9, 11, 30], [0, 9, 13, 23], [0, 9, 14, 18],
            [0, 9, 15, 25], [0, 9, 16, 26], [0, 9, 19, 28], [0, 9, 20, 36], [0, 9, 21, 35],
            [0, 9, 22, 24], [0, 9, 27, 31], [0, 9, 32, 34], [0, 10, 11, 36],
            [0, 10, 12, 15], [0, 10, 13, 26], [0, 10, 14, 16], [0, 10, 17, 37],
            [0, 10, 19, 29], [0, 10, 22, 31], [0, 10, 23, 32], [0, 10, 24, 35],
            [0, 10, 25, 34], [0, 10, 28, 33], [0, 11, 12, 16], [0, 11, 14, 24],
            [0, 11, 15, 26], [0, 11, 17, 31], [0, 11, 19, 21], [0, 11, 23, 34],
            [0, 11, 27, 29], [0, 11, 28, 35], [0, 11, 32, 33], [0, 12, 17, 20],
            [0, 12, 19, 35], [0, 12, 21, 28], [0, 12, 22, 25], [0, 12, 23, 27],
            [0, 12, 26, 29], [0, 12, 30, 33], [0, 12, 31, 34], [0, 12, 36, 37],
            [0, 13, 14, 33], [0, 13, 15, 29], [0, 13, 16, 24], [0, 13, 17, 21],
            [0, 13, 18, 34], [0, 13, 19, 32], [0, 13, 22, 28], [0, 13, 25, 35],
            [0, 14, 17, 26], [0, 14, 19, 20], [0, 14, 21, 32], [0, 14, 23, 36],
            [0, 14, 25, 28], [0, 14, 27, 30], [0, 14, 34, 37], [0, 15, 16, 36],
            [0, 15, 17, 23], [0, 15, 18, 20], [0, 15, 19, 34], [0, 15, 21, 37],
            [0, 15, 24, 28], [0, 15, 30, 32], [0, 16, 17, 32], [0, 16, 18, 27],
            [0, 16, 19, 31], [0, 16, 20, 33], [0, 16, 23, 30], [0, 17, 24, 27],
            [0, 17, 25, 33], [0, 17, 28, 36], [0, 17, 30, 34], [0, 18, 21, 26],
            [0, 18, 22, 29], [0, 18, 23, 28], [0, 18, 25, 31], [0, 18, 30, 35],
            [0, 18, 32, 37], [0, 18, 33, 36], [0, 19, 22, 26], [0, 19, 24, 37],
            [0, 19, 27, 36], [0, 20, 21, 31], [0, 20, 22, 37], [0, 20, 23, 24],
            [0, 20, 26, 30], [0, 20, 28, 32], [0, 20, 29, 34], [0, 21, 22, 27],
            [0, 21, 24, 30], [0, 22, 32, 36], [0, 22, 33, 34], [0, 23, 25, 29],
            [0, 23, 31, 37], [0, 23, 33, 35], [0, 24, 25, 26], [0, 24, 29, 36],
            [0, 25, 27, 37], [0, 25, 30, 36], [0, 26, 27, 33], [0, 26, 32, 35],
            [0, 27, 34, 35], [0, 28, 29, 37], [0, 28, 30, 31], [0, 29, 31, 33],
            [0, 31, 35, 36], [1, 2, 3, 15], [1, 2, 4, 35], [1, 2, 5, 32], [1, 2, 6, 28],
            [1, 2, 7, 18], [1, 2, 8, 13], [1, 2, 9, 37], [1, 2, 10, 11], [1, 2, 12, 19],
            [1, 2, 16, 36], [1, 2, 17, 23], [1, 2, 20, 34], [1, 2, 21, 26], [1, 2, 22, 24],
            [1, 2, 25, 33], [1, 2, 27, 29], [1, 2, 30, 31], [1, 3, 4, 11], [1, 3, 5, 10],
            [1, 3, 6, 29], [1, 3, 7, 16], [1, 3, 8, 37], [1, 3, 9, 24], [1, 3, 12, 23],
            [1, 3, 13, 14], [1, 3, 17, 26], [1, 3, 18, 19], [1, 3, 20, 31], [1, 3, 21, 36],
            [1, 3, 22, 30], [1, 3, 25, 35], [1, 3, 27, 32], [1, 3, 28, 33], [1, 4, 5, 19],
            [1, 4, 6, 24], [1, 4, 7, 33], [1, 4, 8, 20], [1, 4, 9, 21], [1, 4, 10, 18],
            [1, 4, 12, 26], [1, 4, 13, 25], [1, 4, 14, 28], [1, 4, 15, 32], [1, 4, 16, 23],
            [1, 4, 17, 29], [1, 4, 22, 34], [1, 4, 27, 37], [1, 4, 30, 36], [1, 5, 6, 8],
            [1, 5, 7, 29], [1, 5, 9, 26], [1, 5, 11, 31], [1, 5, 12, 21], [1, 5, 13, 33],
            [1, 5, 14, 37], [1, 5, 15, 30], [1, 5, 16, 28], [1, 5, 17, 36], [1, 5, 18, 23],
            [1, 5, 20, 24], [1, 5, 22, 35], [1, 5, 25, 34], [1, 6, 7, 25], [1, 6, 9, 27],
            [1, 6, 10, 30], [1, 6, 11, 21], [1, 6, 12, 14], [1, 6, 13, 15], [1, 6, 16, 34],
            [1, 6, 18, 36], [1, 6, 19, 20], [1, 6, 22, 26], [1, 6, 23, 31], [1, 6, 32, 33],
            [1, 6, 35, 37], [1, 7, 8, 31], [1, 7, 9, 34], [1, 7, 10, 13], [1, 7, 11, 19],
            [1, 7, 14, 32], [1, 7, 15, 36], [1, 7, 17, 30], [1, 7, 20, 26], [1, 7, 21, 28],
            [1, 7, 22, 37], [1, 7, 23, 24], [1, 7, 27, 35], [1, 8, 9, 12], [1, 8, 10, 34],
            [1, 8, 11, 22], [1, 8, 14, 21], [1, 8, 15, 23], [1, 8, 16, 32], [1, 8, 17, 35],
            [1, 8, 18, 30], [1, 8, 19, 25], [1, 8, 24, 27], [1, 8, 26, 33], [1, 8, 28, 29],
            [1, 9, 11, 28], [1, 9, 13, 19], [1, 9, 14, 31], [1, 9, 15, 16], [1, 9, 17, 22],
            [1, 9, 18, 20], [1, 9, 23, 36], [1, 9, 25, 32], [1, 9, 29, 35], [1, 9, 30, 33],
            [1, 10, 12, 31], [1, 10, 14, 24], [1, 10, 15, 19], [1, 10, 16, 26],
            [1, 10, 17, 27], [1, 10, 20, 29], [1, 10, 21, 37], [1, 10, 22, 36],
            [1, 10, 23, 25], [1, 10, 28, 32], [1, 10, 33, 35], [1, 11, 12, 37],
            [1, 11, 13, 16], [1, 11, 14, 27], [1, 11, 15, 17], [1, 11, 20, 30],
            [1, 11, 23, 32], [1, 11, 24, 33], [1, 11, 25, 36], [1, 11, 26, 35],
            [1, 11, 29, 34], [1, 12, 13, 17], [1, 12, 15, 25], [1, 12, 16, 27],
            [1, 12, 18, 32], [1, 12, 20, 22], [1, 12, 24, 35], [1, 12, 28, 30],
            [1, 12, 29, 36], [1, 12, 33, 34], [1, 13, 18, 21], [1, 13, 20, 36],
            [1, 13, 22, 29], [1, 13, 23, 26], [1, 13, 24, 28], [1, 13, 27, 30],
            [1, 13, 31, 34], [1, 13, 32, 35], [1, 14, 15, 34], [1, 14, 16, 30],
            [1, 14, 17, 25], [1, 14, 18, 22], [1, 14, 19, 35], [1, 14, 20, 33],
            [1, 14, 23, 29], [1, 14, 26, 36], [1, 15, 18, 27], [1, 15, 20, 21],
            [1, 15, 22, 33], [1, 15, 24, 37], [1, 15, 26, 29], [1, 15, 28, 31],
            [1, 16, 17, 37], [1, 16, 18, 24], [1, 16, 19, 21], [1, 16, 20, 35],
            [1, 16, 25, 29], [1, 16, 31, 33], [1, 17, 18, 33], [1, 17, 19, 28],
            [1, 17, 20, 32], [1, 17, 21, 34], [1, 17, 24, 31], [1, 18, 25, 28],
            [1, 18, 26, 34], [1, 18, 29, 37], [1, 18, 31, 35], [1, 19, 22, 27],
            [1, 19, 23, 30], [1, 19, 24, 29], [1, 19, 26, 32], [1, 19, 31, 36],
            [1, 19, 34, 37], [1, 20, 23, 27], [1, 20, 28, 37], [1, 21, 22, 32],
            [1, 21, 24, 25], [1, 21, 27, 31], [1, 21, 29, 33], [1, 21, 30, 35],
            [1, 22, 23, 28], [1, 22, 25, 31], [1, 23, 33, 37], [1, 23, 34, 35],
            [1, 24, 26, 30], [1, 24, 34, 36], [1, 25, 26, 27], [1, 25, 30, 37],
            [1, 26, 31, 37], [1, 27, 28, 34], [1, 27, 33, 36], [1, 28, 35, 36],
            [1, 29, 31, 32], [1, 30, 32, 34], [1, 32, 36, 37], [2, 3, 4, 16],
            [2, 3, 5, 36], [2, 3, 6, 33], [2, 3, 7, 29], [2, 3, 8, 19], [2, 3, 9, 14],
            [2, 3, 11, 12], [2, 3, 13, 20], [2, 3, 17, 37], [2, 3, 18, 24], [2, 3, 21, 35],
            [2, 3, 22, 27], [2, 3, 23, 25], [2, 3, 26, 34], [2, 3, 28, 30], [2, 3, 31, 32],
            [2, 4, 5, 12], [2, 4, 6, 11], [2, 4, 7, 30], [2, 4, 8, 17], [2, 4, 10, 25],
            [2, 4, 13, 24], [2, 4, 14, 15], [2, 4, 18, 27], [2, 4, 19, 20], [2, 4, 21, 32],
            [2, 4, 22, 37], [2, 4, 23, 31], [2, 4, 26, 36], [2, 4, 28, 33], [2, 4, 29, 34],
            [2, 5, 6, 20], [2, 5, 7, 25], [2, 5, 8, 34], [2, 5, 9, 21], [2, 5, 10, 22],
            [2, 5, 11, 19], [2, 5, 13, 27], [2, 5, 14, 26], [2, 5, 15, 29], [2, 5, 16, 33],
            [2, 5, 17, 24], [2, 5, 18, 30], [2, 5, 23, 35], [2, 5, 31, 37], [2, 6, 7, 9],
            [2, 6, 8, 30], [2, 6, 10, 27], [2, 6, 12, 32], [2, 6, 13, 22], [2, 6, 14, 34],
            [2, 6, 16, 31], [2, 6, 17, 29], [2, 6, 18, 37], [2, 6, 19, 24], [2, 6, 21, 25],
            [2, 6, 23, 36], [2, 6, 26, 35], [2, 7, 8, 26], [2, 7, 10, 28], [2, 7, 11, 31],
            [2, 7, 12, 22], [2, 7, 13, 15], [2, 7, 14, 16], [2, 7, 17, 35], [2, 7, 19, 37],
            [2, 7, 20, 21], [2, 7, 23, 27], [2, 7, 24, 32], [2, 7, 33, 34], [2, 8, 9, 32],
            [2, 8, 10, 35], [2, 8, 11, 14], [2, 8, 12, 20], [2, 8, 15, 33], [2, 8, 16, 37],
            [2, 8, 18, 31], [2, 8, 21, 27], [2, 8, 22, 29], [2, 8, 24, 25], [2, 8, 28, 36],
            [2, 9, 10, 13], [2, 9, 11, 35], [2, 9, 12, 23], [2, 9, 15, 22], [2, 9, 16, 24],
            [2, 9, 17, 33], [2, 9, 18, 36], [2, 9, 19, 31], [2, 9, 20, 26], [2, 9, 25, 28],
            [2, 9, 27, 34], [2, 9, 29, 30], [2, 10, 12, 29], [2, 10, 14, 20],
            [2, 10, 15, 32], [2, 10, 16, 17], [2, 10, 18, 23], [2, 10, 19, 21],
            [2, 10, 24, 37], [2, 10, 26, 33], [2, 10, 30, 36], [2, 10, 31, 34],
            [2, 11, 13, 32], [2, 11, 15, 25], [2, 11, 16, 20], [2, 11, 17, 27],
            [2, 11, 18, 28], [2, 11, 21, 30], [2, 11, 23, 37], [2, 11, 24, 26],
            [2, 11, 29, 33], [2, 11, 34, 36], [2, 12, 14, 17], [2, 12, 15, 28],
            [2, 12, 16, 18], [2, 12, 21, 31], [2, 12, 24, 33], [2, 12, 25, 34],
            [2, 12, 26, 37], [2, 12, 27, 36], [2, 12, 30, 35], [2, 13, 14, 18],
            [2, 13, 16, 26], [2, 13, 17, 28], [2, 13, 19, 33], [2, 13, 21, 23],
            [2, 13, 25, 36], [2, 13, 29, 31], [2, 13, 30, 37], [2, 13, 34, 35],
            [2, 14, 19, 22], [2, 14, 21, 37], [2, 14, 23, 30], [2, 14, 24, 27],
            [2, 14, 25, 29], [2, 14, 28, 31], [2, 14, 32, 35], [2, 14, 33, 36],
            [2, 15, 16, 35], [2, 15, 17, 31], [2, 15, 18, 26], [2, 15, 19, 23],
            [2, 15, 20, 36], [2, 15, 21, 34], [2, 15, 24, 30], [2, 15, 27, 37],
            [2, 16, 19, 28], [2, 16, 21, 22], [2, 16, 23, 34], [2, 16, 27, 30],
            [2, 16, 29, 32], [2, 17, 19, 25], [2, 17, 20, 22], [2, 17, 21, 36],
            [2, 17, 26, 30], [2, 17, 32, 34], [2, 18, 19, 34], [2, 18, 20, 29],
            [2, 18, 21, 33], [2, 18, 22, 35], [2, 18, 25, 32], [2, 19, 26, 29],
            [2, 19, 27, 35], [2, 19, 32, 36], [2, 20, 23, 28], [2, 20, 24, 31],
            [2, 20, 25, 30], [2, 20, 27, 33], [2, 20, 32, 37], [2, 21, 24, 28],
            [2, 22, 23, 33], [2, 22, 25, 26], [2, 22, 28, 32], [2, 22, 30, 34],
            [2, 22, 31, 36], [2, 23, 24, 29], [2, 23, 26, 32], [2, 24, 35, 36],
            [2, 25, 27, 31], [2, 25, 35, 37], [2, 26, 27, 28], [2, 28, 29, 35],
            [2, 28, 34, 37], [2, 29, 36, 37], [2, 30, 32, 33], [2, 31, 33, 35],
            [3, 4, 5, 17], [3, 4, 6, 37], [3, 4, 7, 34], [3, 4, 8, 30], [3, 4, 9, 20],
            [3, 4, 10, 15], [3, 4, 12, 13], [3, 4, 14, 21], [3, 4, 19, 25], [3, 4, 22, 36],
            [3, 4, 23, 28], [3, 4, 24, 26], [3, 4, 27, 35], [3, 4, 29, 31], [3, 4, 32, 33],
            [3, 5, 6, 13], [3, 5, 7, 12], [3, 5, 8, 31], [3, 5, 9, 18], [3, 5, 11, 26],
            [3, 5, 14, 25], [3, 5, 15, 16], [3, 5, 19, 28], [3, 5, 20, 21], [3, 5, 22, 33],
            [3, 5, 24, 32], [3, 5, 27, 37], [3, 5, 29, 34], [3, 5, 30, 35], [3, 6, 7, 21],
            [3, 6, 8, 26], [3, 6, 9, 35], [3, 6, 10, 22], [3, 6, 11, 23], [3, 6, 12, 20],
            [3, 6, 14, 28], [3, 6, 15, 27], [3, 6, 16, 30], [3, 6, 17, 34], [3, 6, 18, 25],
            [3, 6, 19, 31], [3, 6, 24, 36], [3, 7, 8, 10], [3, 7, 9, 31], [3, 7, 11, 28],
            [3, 7, 13, 33], [3, 7, 14, 23], [3, 7, 15, 35], [3, 7, 17, 32], [3, 7, 18, 30],
            [3, 7, 20, 25], [3, 7, 22, 26], [3, 7, 24, 37], [3, 7, 27, 36], [3, 8, 9, 27],
            [3, 8, 11, 29], [3, 8, 12, 32], [3, 8, 13, 23], [3, 8, 14, 16], [3, 8, 15, 17],
            [3, 8, 18, 36], [3, 8, 21, 22], [3, 8, 24, 28], [3, 8, 25, 33], [3, 8, 34, 35],
            [3, 9, 10, 33], [3, 9, 11, 36], [3, 9, 12, 15], [3, 9, 13, 21], [3, 9, 16, 34],
            [3, 9, 19, 32], [3, 9, 22, 28], [3, 9, 23, 30], [3, 9, 25, 26], [3, 9, 29, 37],
            [3, 10, 11, 14], [3, 10, 12, 36], [3, 10, 13, 24], [3, 10, 16, 23],
            [3, 10, 17, 25], [3, 10, 18, 34], [3, 10, 19, 37], [3, 10, 20, 32],
            [3, 10, 21, 27], [3, 10, 26, 29], [3, 10, 28, 35], [3, 10, 30, 31],
            [3, 11, 13, 30], [3, 11, 15, 21], [3, 11, 16, 33], [3, 11, 17, 18],
            [3, 11, 19, 24], [3, 11, 20, 22], [3, 11, 27, 34], [3, 11, 31, 37],
            [3, 11, 32, 35], [3, 12, 14, 33], [3, 12, 16, 26], [3, 12, 17, 21],
            [3, 12, 18, 28], [3, 12, 19, 29], [3, 12, 22, 31], [3, 12, 25, 27],
            [3, 12, 30, 34], [3, 12, 35, 37], [3, 13, 15, 18], [3, 13, 16, 29],
            [3, 13, 17, 19], [3, 13, 22, 32], [3, 13, 25, 34], [3, 13, 26, 35],
            [3, 13, 28, 37], [3, 13, 31, 36], [3, 14, 15, 19], [3, 14, 17, 27],
            [3, 14, 18, 29], [3, 14, 20, 34], [3, 14, 22, 24], [3, 14, 26, 37],
            [3, 14, 30, 32], [3, 14, 35, 36], [3, 15, 20, 23], [3, 15, 24, 31],
            [3, 15, 25, 28], [3, 15, 26, 30], [3, 15, 29, 32], [3, 15, 33, 36],
            [3, 15, 34, 37], [3, 16, 17, 36], [3, 16, 18, 32], [3, 16, 19, 27],
            [3, 16, 20, 24], [3, 16, 21, 37], [3, 16, 22, 35], [3, 16, 25, 31],
            [3, 17, 20, 29], [3, 17, 22, 23], [3, 17, 24, 35], [3, 17, 28, 31],
            [3, 17, 30, 33], [3, 18, 20, 26], [3, 18, 21, 23], [3, 18, 22, 37],
            [3, 18, 27, 31], [3, 18, 33, 35], [3, 19, 20, 35], [3, 19, 21, 30],
            [3, 19, 22, 34], [3, 19, 23, 36], [3, 19, 26, 33], [3, 20, 27, 30],
            [3, 20, 28, 36], [3, 20, 33, 37], [3, 21, 24, 29], [3, 21, 25, 32],
            [3, 21, 26, 31], [3, 21, 28, 34], [3, 22, 25, 29], [3, 23, 24, 34],
            [3, 23, 26, 27], [3, 23, 29, 33], [3, 23, 31, 35], [3, 23, 32, 37],
            [3, 24, 25, 30], [3, 24, 27, 33], [3, 25, 36, 37], [3, 26, 28, 32],
            [3, 27, 28, 29], [3, 29, 30, 36], [3, 31, 33, 34], [3, 32, 34, 36],
            [4, 5, 6, 18], [4, 5, 8, 35], [4, 5, 9, 31], [4, 5, 10, 21], [4, 5, 11, 16],
            [4, 5, 13, 14], [4, 5, 15, 22], [4, 5, 20, 26], [4, 5, 23, 37], [4, 5, 24, 29],
            [4, 5, 25, 27], [4, 5, 28, 36], [4, 5, 30, 32], [4, 5, 33, 34], [4, 6, 7, 14],
            [4, 6, 8, 13], [4, 6, 9, 32], [4, 6, 10, 19], [4, 6, 12, 27], [4, 6, 15, 26],
            [4, 6, 16, 17], [4, 6, 20, 29], [4, 6, 21, 22], [4, 6, 23, 34], [4, 6, 25, 33],
            [4, 6, 30, 35], [4, 6, 31, 36], [4, 7, 8, 22], [4, 7, 9, 27], [4, 7, 10, 36],
            [4, 7, 11, 23], [4, 7, 12, 24], [4, 7, 13, 21], [4, 7, 15, 29], [4, 7, 16, 28],
            [4, 7, 17, 31], [4, 7, 18, 35], [4, 7, 19, 26], [4, 7, 20, 32], [4, 7, 25, 37],
            [4, 8, 9, 11], [4, 8, 10, 32], [4, 8, 12, 29], [4, 8, 14, 34], [4, 8, 15, 24],
            [4, 8, 16, 36], [4, 8, 18, 33], [4, 8, 19, 31], [4, 8, 21, 26], [4, 8, 23, 27],
            [4, 8, 28, 37], [4, 9, 10, 28], [4, 9, 12, 30], [4, 9, 13, 33], [4, 9, 14, 24],
            [4, 9, 15, 17], [4, 9, 16, 18], [4, 9, 19, 37], [4, 9, 22, 23], [4, 9, 25, 29],
            [4, 9, 26, 34], [4, 9, 35, 36], [4, 10, 11, 34], [4, 10, 12, 37],
            [4, 10, 13, 16], [4, 10, 14, 22], [4, 10, 17, 35], [4, 10, 20, 33],
            [4, 10, 23, 29], [4, 10, 24, 31], [4, 10, 26, 27], [4, 11, 12, 15],
            [4, 11, 13, 37], [4, 11, 14, 25], [4, 11, 17, 24], [4, 11, 18, 26],
            [4, 11, 19, 35], [4, 11, 21, 33], [4, 11, 22, 28], [4, 11, 27, 30],
            [4, 11, 29, 36], [4, 11, 31, 32], [4, 12, 14, 31], [4, 12, 16, 22],
            [4, 12, 17, 34], [4, 12, 18, 19], [4, 12, 20, 25], [4, 12, 21, 23],
            [4, 12, 28, 35], [4, 12, 33, 36], [4, 13, 15, 34], [4, 13, 17, 27],
            [4, 13, 18, 22], [4, 13, 19, 29], [4, 13, 20, 30], [4, 13, 23, 32],
            [4, 13, 26, 28], [4, 13, 31, 35], [4, 14, 16, 19], [4, 14, 17, 30],
            [4, 14, 18, 20], [4, 14, 23, 33], [4, 14, 26, 35], [4, 14, 27, 36],
            [4, 14, 32, 37], [4, 15, 16, 20], [4, 15, 18, 28], [4, 15, 19, 30],
            [4, 15, 21, 35], [4, 15, 23, 25], [4, 15, 31, 33], [4, 15, 36, 37],
            [4, 16, 21, 24], [4, 16, 25, 32], [4, 16, 26, 29], [4, 16, 27, 31],
            [4, 16, 30, 33], [4, 16, 34, 37], [4, 17, 18, 37], [4, 17, 19, 33],
            [4, 17, 20, 28], [4, 17, 21, 25], [4, 17, 23, 36], [4, 17, 26, 32],
            [4, 18, 21, 30], [4, 18, 23, 24], [4, 18, 25, 36], [4, 18, 29, 32],
            [4, 18, 31, 34], [4, 19, 21, 27], [4, 19, 22, 24], [4, 19, 28, 32],
            [4, 19, 34, 36], [4, 20, 21, 36], [4, 20, 22, 31], [4, 20, 23, 35],
            [4, 20, 24, 37], [4, 20, 27, 34], [4, 21, 28, 31], [4, 21, 29, 37],
            [4, 22, 25, 30], [4, 22, 26, 33], [4, 22, 27, 32], [4, 22, 29, 35],
            [4, 23, 26, 30], [4, 24, 25, 35], [4, 24, 27, 28], [4, 24, 30, 34],
            [4, 24, 32, 36], [4, 25, 26, 31], [4, 25, 28, 34], [4, 27, 29, 33],
            [4, 28, 29, 30], [4, 30, 31, 37], [4, 32, 34, 35], [4, 33, 35, 37],
            [5, 6, 7, 19], [5, 6, 9, 36], [5, 6, 10, 32], [5, 6, 11, 22], [5, 6, 12, 17],
            [5, 6, 14, 15], [5, 6, 16, 23], [5, 6, 21, 27], [5, 6, 25, 30], [5, 6, 26, 28],
            [5, 6, 29, 37], [5, 6, 31, 33], [5, 6, 34, 35], [5, 7, 8, 15], [5, 7, 9, 14],
            [5, 7, 10, 33], [5, 7, 11, 20], [5, 7, 13, 28], [5, 7, 16, 27], [5, 7, 17, 18],
            [5, 7, 21, 30], [5, 7, 22, 23], [5, 7, 24, 35], [5, 7, 26, 34], [5, 7, 31, 36],
            [5, 7, 32, 37], [5, 8, 9, 23], [5, 8, 10, 28], [5, 8, 11, 37], [5, 8, 12, 24],
            [5, 8, 13, 25], [5, 8, 14, 22], [5, 8, 16, 30], [5, 8, 17, 29], [5, 8, 18, 32],
            [5, 8, 19, 36], [5, 8, 20, 27], [5, 8, 21, 33], [5, 9, 10, 12], [5, 9, 11, 33],
            [5, 9, 13, 30], [5, 9, 15, 35], [5, 9, 16, 25], [5, 9, 17, 37], [5, 9, 19, 34],
            [5, 9, 20, 32], [5, 9, 22, 27], [5, 9, 24, 28], [5, 10, 11, 29],
            [5, 10, 13, 31], [5, 10, 14, 34], [5, 10, 15, 25], [5, 10, 16, 18],
            [5, 10, 17, 19], [5, 10, 23, 24], [5, 10, 26, 30], [5, 10, 27, 35],
            [5, 10, 36, 37], [5, 11, 12, 35], [5, 11, 14, 17], [5, 11, 15, 23],
            [5, 11, 18, 36], [5, 11, 21, 34], [5, 11, 24, 30], [5, 11, 25, 32],
            [5, 11, 27, 28], [5, 12, 13, 16], [5, 12, 15, 26], [5, 12, 18, 25],
            [5, 12, 19, 27], [5, 12, 20, 36], [5, 12, 22, 34], [5, 12, 23, 29],
            [5, 12, 28, 31], [5, 12, 30, 37], [5, 12, 32, 33], [5, 13, 15, 32],
            [5, 13, 17, 23], [5, 13, 18, 35], [5, 13, 19, 20], [5, 13, 21, 26],
            [5, 13, 22, 24], [5, 13, 29, 36], [5, 13, 34, 37], [5, 14, 16, 35],
            [5, 14, 18, 28], [5, 14, 19, 23], [5, 14, 20, 30], [5, 14, 21, 31],
            [5, 14, 24, 33], [5, 14, 27, 29], [5, 14, 32, 36], [5, 15, 17, 20],
            [5, 15, 18, 31], [5, 15, 19, 21], [5, 15, 24, 34], [5, 15, 27, 36],
            [5, 15, 28, 37], [5, 16, 17, 21], [5, 16, 19, 29], [5, 16, 20, 31],
            [5, 16, 22, 36], [5, 16, 24, 26], [5, 16, 32, 34], [5, 17, 22, 25],
            [5, 17, 26, 33], [5, 17, 27, 30], [5, 17, 28, 32], [5, 17, 31, 34],
            [5, 18, 20, 34], [5, 18, 21, 29], [5, 18, 22, 26], [5, 18, 24, 37],
            [5, 18, 27, 33], [5, 19, 22, 31], [5, 19, 24, 25], [5, 19, 26, 37],
            [5, 19, 30, 33], [5, 19, 32, 35], [5, 20, 22, 28], [5, 20, 23, 25],
            [5, 20, 29, 33], [5, 20, 35, 37], [5, 21, 22, 37], [5, 21, 23, 32],
            [5, 21, 24, 36], [5, 21, 28, 35], [5, 22, 29, 32], [5, 23, 26, 31],
            [5, 23, 27, 34], [5, 23, 28, 33], [5, 23, 30, 36], [5, 24, 27, 31],
            [5, 25, 26, 36], [5, 25, 28, 29], [5, 25, 31, 35], [5, 25, 33, 37],
            [5, 26, 27, 32], [5, 26, 29, 35], [5, 28, 30, 34], [5, 29, 30, 31],
            [5, 33, 35, 36], [6, 7, 8, 20], [6, 7, 10, 37], [6, 7, 11, 33], [6, 7, 12, 23],
            [6, 7, 13, 18], [6, 7, 15, 16], [6, 7, 17, 24], [6, 7, 22, 28], [6, 7, 26, 31],
            [6, 7, 27, 29], [6, 7, 32, 34], [6, 7, 35, 36], [6, 8, 9, 16], [6, 8, 10, 15],
            [6, 8, 11, 34], [6, 8, 12, 21], [6, 8, 14, 29], [6, 8, 17, 28], [6, 8, 18, 19],
            [6, 8, 22, 31], [6, 8, 23, 24], [6, 8, 25, 36], [6, 8, 27, 35], [6, 8, 32, 37],
            [6, 9, 10, 24], [6, 9, 11, 29], [6, 9, 13, 25], [6, 9, 14, 26], [6, 9, 15, 23],
            [6, 9, 17, 31], [6, 9, 18, 30], [6, 9, 19, 33], [6, 9, 20, 37], [6, 9, 21, 28],
            [6, 9, 22, 34], [6, 10, 11, 13], [6, 10, 12, 34], [6, 10, 14, 31],
            [6, 10, 16, 36], [6, 10, 17, 26], [6, 10, 20, 35], [6, 10, 21, 33],
            [6, 10, 23, 28], [6, 10, 25, 29], [6, 11, 12, 30], [6, 11, 14, 32],
            [6, 11, 15, 35], [6, 11, 16, 26], [6, 11, 17, 19], [6, 11, 18, 20],
            [6, 11, 24, 25], [6, 11, 27, 31], [6, 11, 28, 36], [6, 12, 13, 36],
            [6, 12, 15, 18], [6, 12, 16, 24], [6, 12, 19, 37], [6, 12, 22, 35],
            [6, 12, 25, 31], [6, 12, 26, 33], [6, 12, 28, 29], [6, 13, 14, 17],
            [6, 13, 16, 27], [6, 13, 19, 26], [6, 13, 20, 28], [6, 13, 21, 37],
            [6, 13, 23, 35], [6, 13, 24, 30], [6, 13, 29, 32], [6, 13, 33, 34],
            [6, 14, 16, 33], [6, 14, 18, 24], [6, 14, 19, 36], [6, 14, 20, 21],
            [6, 14, 22, 27], [6, 14, 23, 25], [6, 14, 30, 37], [6, 15, 17, 36],
            [6, 15, 19, 29], [6, 15, 20, 24], [6, 15, 21, 31], [6, 15, 22, 32],
            [6, 15, 25, 34], [6, 15, 28, 30], [6, 15, 33, 37], [6, 16, 18, 21],
            [6, 16, 19, 32], [6, 16, 20, 22], [6, 16, 25, 35], [6, 16, 28, 37],
            [6, 17, 18, 22], [6, 17, 20, 30], [6, 17, 21, 32], [6, 17, 23, 37],
            [6, 17, 25, 27], [6, 17, 33, 35], [6, 18, 23, 26], [6, 18, 27, 34],
            [6, 18, 28, 31], [6, 18, 29, 33], [6, 18, 32, 35], [6, 19, 21, 35],
            [6, 19, 22, 30], [6, 19, 23, 27], [6, 19, 28, 34], [6, 20, 23, 32],
            [6, 20, 25, 26], [6, 20, 31, 34], [6, 20, 33, 36], [6, 21, 23, 29],
            [6, 21, 24, 26], [6, 21, 30, 34], [6, 22, 24, 33], [6, 22, 25, 37],
            [6, 22, 29, 36], [6, 23, 30, 33], [6, 24, 27, 32], [6, 24, 28, 35],
            [6, 24, 29, 34], [6, 24, 31, 37], [6, 25, 28, 32], [6, 26, 27, 37],
            [6, 26, 29, 30], [6, 26, 32, 36], [6, 27, 28, 33], [6, 27, 30, 36],
            [6, 29, 31, 35], [6, 30, 31, 32], [6, 34, 36, 37], [7, 8, 9, 21],
            [7, 8, 12, 34], [7, 8, 13, 24], [7, 8, 14, 19], [7, 8, 16, 17], [7, 8, 18, 25],
            [7, 8, 23, 29], [7, 8, 27, 32], [7, 8, 28, 30], [7, 8, 33, 35], [7, 8, 36, 37],
            [7, 9, 10, 17], [7, 9, 11, 16], [7, 9, 12, 35], [7, 9, 13, 22], [7, 9, 15, 30],
            [7, 9, 18, 29], [7, 9, 19, 20], [7, 9, 23, 32], [7, 9, 24, 25], [7, 9, 26, 37],
            [7, 9, 28, 36], [7, 10, 11, 25], [7, 10, 12, 30], [7, 10, 14, 26],
            [7, 10, 15, 27], [7, 10, 16, 24], [7, 10, 18, 32], [7, 10, 19, 31],
            [7, 10, 20, 34], [7, 10, 22, 29], [7, 10, 23, 35], [7, 11, 12, 14],
            [7, 11, 13, 35], [7, 11, 15, 32], [7, 11, 17, 37], [7, 11, 18, 27],
            [7, 11, 21, 36], [7, 11, 22, 34], [7, 11, 24, 29], [7, 11, 26, 30],
            [7, 12, 13, 31], [7, 12, 15, 33], [7, 12, 16, 36], [7, 12, 17, 27],
            [7, 12, 18, 20], [7, 12, 19, 21], [7, 12, 25, 26], [7, 12, 28, 32],
            [7, 12, 29, 37], [7, 13, 14, 37], [7, 13, 16, 19], [7, 13, 17, 25],
            [7, 13, 23, 36], [7, 13, 26, 32], [7, 13, 27, 34], [7, 13, 29, 30],
            [7, 14, 15, 18], [7, 14, 17, 28], [7, 14, 20, 27], [7, 14, 21, 29],
            [7, 14, 24, 36], [7, 14, 25, 31], [7, 14, 30, 33], [7, 14, 34, 35],
            [7, 15, 17, 34], [7, 15, 19, 25], [7, 15, 20, 37], [7, 15, 21, 22],
            [7, 15, 23, 28], [7, 15, 24, 26], [7, 16, 18, 37], [7, 16, 20, 30],
            [7, 16, 21, 25], [7, 16, 22, 32], [7, 16, 23, 33], [7, 16, 26, 35],
            [7, 16, 29, 31], [7, 17, 19, 22], [7, 17, 20, 33], [7, 17, 21, 23],
            [7, 17, 26, 36], [7, 18, 19, 23], [7, 18, 21, 31], [7, 18, 22, 33],
            [7, 18, 26, 28], [7, 18, 34, 36], [7, 19, 24, 27], [7, 19, 28, 35],
            [7, 19, 29, 32], [7, 19, 30, 34], [7, 19, 33, 36], [7, 20, 22, 36],
            [7, 20, 23, 31], [7, 20, 24, 28], [7, 20, 29, 35], [7, 21, 24, 33],
            [7, 21, 26, 27], [7, 21, 32, 35], [7, 21, 34, 37], [7, 22, 24, 30],
            [7, 22, 25, 27], [7, 22, 31, 35], [7, 23, 25, 34], [7, 23, 30, 37],
            [7, 24, 31, 34], [7, 25, 28, 33], [7, 25, 29, 36], [7, 25, 30, 35],
            [7, 26, 29, 33], [7, 27, 30, 31], [7, 27, 33, 37], [7, 28, 29, 34],
            [7, 28, 31, 37], [7, 30, 32, 36], [7, 31, 32, 33], [8, 9, 10, 22],
            [8, 9, 13, 35], [8, 9, 14, 25], [8, 9, 15, 20], [8, 9, 17, 18], [8, 9, 19, 26],
            [8, 9, 24, 30], [8, 9, 28, 33], [8, 9, 29, 31], [8, 9, 34, 36], [8, 10, 11, 18],
            [8, 10, 12, 17], [8, 10, 13, 36], [8, 10, 14, 23], [8, 10, 16, 31],
            [8, 10, 19, 30], [8, 10, 20, 21], [8, 10, 24, 33], [8, 10, 25, 26],
            [8, 10, 29, 37], [8, 11, 12, 26], [8, 11, 13, 31], [8, 11, 15, 27],
            [8, 11, 16, 28], [8, 11, 17, 25], [8, 11, 19, 33], [8, 11, 20, 32],
            [8, 11, 21, 35], [8, 11, 23, 30], [8, 11, 24, 36], [8, 12, 13, 15],
            [8, 12, 14, 36], [8, 12, 16, 33], [8, 12, 19, 28], [8, 12, 22, 37],
            [8, 12, 23, 35], [8, 12, 25, 30], [8, 12, 27, 31], [8, 13, 14, 32],
            [8, 13, 16, 34], [8, 13, 17, 37], [8, 13, 18, 28], [8, 13, 19, 21],
            [8, 13, 20, 22], [8, 13, 26, 27], [8, 13, 29, 33], [8, 14, 17, 20],
            [8, 14, 18, 26], [8, 14, 24, 37], [8, 14, 27, 33], [8, 14, 28, 35],
            [8, 14, 30, 31], [8, 15, 16, 19], [8, 15, 18, 29], [8, 15, 21, 28],
            [8, 15, 22, 30], [8, 15, 25, 37], [8, 15, 26, 32], [8, 15, 31, 34],
            [8, 15, 35, 36], [8, 16, 18, 35], [8, 16, 20, 26], [8, 16, 22, 23],
            [8, 16, 24, 29], [8, 16, 25, 27], [8, 17, 21, 31], [8, 17, 22, 26],
            [8, 17, 23, 33], [8, 17, 24, 34], [8, 17, 27, 36], [8, 17, 30, 32],
            [8, 18, 20, 23], [8, 18, 21, 34], [8, 18, 22, 24], [8, 18, 27, 37],
            [8, 19, 20, 24], [8, 19, 22, 32], [8, 19, 23, 34], [8, 19, 27, 29],
            [8, 19, 35, 37], [8, 20, 25, 28], [8, 20, 29, 36], [8, 20, 30, 33],
            [8, 20, 31, 35], [8, 20, 34, 37], [8, 21, 23, 37], [8, 21, 24, 32],
            [8, 21, 25, 29], [8, 21, 30, 36], [8, 22, 25, 34], [8, 22, 27, 28],
            [8, 22, 33, 36], [8, 23, 25, 31], [8, 23, 26, 28], [8, 23, 32, 36],
            [8, 24, 26, 35], [8, 25, 32, 35], [8, 26, 29, 34], [8, 26, 30, 37],
            [8, 26, 31, 36], [8, 27, 30, 34], [8, 28, 31, 32], [8, 29, 30, 35],
            [8, 31, 33, 37], [8, 32, 33, 34], [9, 10, 11, 23], [9, 10, 14, 36],
            [9, 10, 15, 26], [9, 10, 16, 21], [9, 10, 18, 19], [9, 10, 20, 27],
            [9, 10, 25, 31], [9, 10, 29, 34], [9, 10, 30, 32], [9, 10, 35, 37],
            [9, 11, 12, 19], [9, 11, 13, 18], [9, 11, 14, 37], [9, 11, 15, 24],
            [9, 11, 17, 32], [9, 11, 20, 31], [9, 11, 21, 22], [9, 11, 25, 34],
            [9, 11, 26, 27], [9, 12, 13, 27], [9, 12, 14, 32], [9, 12, 16, 28],
            [9, 12, 17, 29], [9, 12, 18, 26], [9, 12, 20, 34], [9, 12, 21, 33],
            [9, 12, 22, 36], [9, 12, 24, 31], [9, 12, 25, 37], [9, 13, 14, 16],
            [9, 13, 15, 37], [9, 13, 17, 34], [9, 13, 20, 29], [9, 13, 24, 36],
            [9, 13, 26, 31], [9, 13, 28, 32], [9, 14, 15, 33], [9, 14, 17, 35],
            [9, 14, 19, 29], [9, 14, 20, 22], [9, 14, 21, 23], [9, 14, 27, 28],
            [9, 14, 30, 34], [9, 15, 18, 21], [9, 15, 19, 27], [9, 15, 28, 34],
            [9, 15, 29, 36], [9, 15, 31, 32], [9, 16, 17, 20], [9, 16, 19, 30],
            [9, 16, 22, 29], [9, 16, 23, 31], [9, 16, 27, 33], [9, 16, 32, 35],
            [9, 16, 36, 37], [9, 17, 19, 36], [9, 17, 21, 27], [9, 17, 23, 24],
            [9, 17, 25, 30], [9, 17, 26, 28], [9, 18, 22, 32], [9, 18, 23, 27],
            [9, 18, 24, 34], [9, 18, 25, 35], [9, 18, 28, 37], [9, 18, 31, 33],
            [9, 19, 21, 24], [9, 19, 22, 35], [9, 19, 23, 25], [9, 20, 21, 25],
            [9, 20, 23, 33], [9, 20, 24, 35], [9, 20, 28, 30], [9, 21, 26, 29],
            [9, 21, 30, 37], [9, 21, 31, 34], [9, 21, 32, 36], [9, 22, 25, 33],
            [9, 22, 26, 30], [9, 22, 31, 37], [9, 23, 26, 35], [9, 23, 28, 29],
            [9, 23, 34, 37], [9, 24, 26, 32], [9, 24, 27, 29], [9, 24, 33, 37],
            [9, 25, 27, 36], [9, 26, 33, 36], [9, 27, 30, 35], [9, 27, 32, 37],
            [9, 28, 31, 35], [9, 29, 32, 33], [9, 30, 31, 36], [9, 33, 34, 35],
            [10, 11, 12, 24], [10, 11, 15, 37], [10, 11, 16, 27], [10, 11, 17, 22],
            [10, 11, 19, 20], [10, 11, 21, 28], [10, 11, 26, 32], [10, 11, 30, 35],
            [10, 11, 31, 33], [10, 12, 13, 20], [10, 12, 14, 19], [10, 12, 16, 25],
            [10, 12, 18, 33], [10, 12, 21, 32], [10, 12, 22, 23], [10, 12, 26, 35],
            [10, 12, 27, 28], [10, 13, 14, 28], [10, 13, 15, 33], [10, 13, 17, 29],
            [10, 13, 18, 30], [10, 13, 19, 27], [10, 13, 21, 35], [10, 13, 22, 34],
            [10, 13, 23, 37], [10, 13, 25, 32], [10, 14, 15, 17], [10, 14, 18, 35],
            [10, 14, 21, 30], [10, 14, 25, 37], [10, 14, 27, 32], [10, 14, 29, 33],
            [10, 15, 16, 34], [10, 15, 18, 36], [10, 15, 20, 30], [10, 15, 21, 23],
            [10, 15, 22, 24], [10, 15, 28, 29], [10, 15, 31, 35], [10, 16, 19, 22],
            [10, 16, 20, 28], [10, 16, 29, 35], [10, 16, 30, 37], [10, 16, 32, 33],
            [10, 17, 18, 21], [10, 17, 20, 31], [10, 17, 23, 30], [10, 17, 24, 32],
            [10, 17, 28, 34], [10, 17, 33, 36], [10, 18, 20, 37], [10, 18, 22, 28],
            [10, 18, 24, 25], [10, 18, 26, 31], [10, 18, 27, 29], [10, 19, 23, 33],
            [10, 19, 24, 28], [10, 19, 25, 35], [10, 19, 26, 36], [10, 19, 32, 34],
            [10, 20, 22, 25], [10, 20, 23, 36], [10, 20, 24, 26], [10, 21, 22, 26],
            [10, 21, 24, 34], [10, 21, 25, 36], [10, 21, 29, 31], [10, 22, 27, 30],
            [10, 22, 32, 35], [10, 22, 33, 37], [10, 23, 26, 34], [10, 23, 27, 31],
            [10, 24, 27, 36], [10, 24, 29, 30], [10, 25, 27, 33], [10, 25, 28, 30],
            [10, 26, 28, 37], [10, 27, 34, 37], [10, 28, 31, 36], [10, 29, 32, 36],
            [10, 30, 33, 34], [10, 31, 32, 37], [10, 34, 35, 36], [11, 12, 13, 25],
            [11, 12, 17, 28], [11, 12, 18, 23], [11, 12, 20, 21], [11, 12, 22, 29],
            [11, 12, 27, 33], [11, 12, 31, 36], [11, 12, 32, 34], [11, 13, 14, 21],
            [11, 13, 15, 20], [11, 13, 17, 26], [11, 13, 19, 34], [11, 13, 22, 33],
            [11, 13, 23, 24], [11, 13, 27, 36], [11, 13, 28, 29], [11, 14, 15, 29],
            [11, 14, 16, 34], [11, 14, 18, 30], [11, 14, 19, 31], [11, 14, 20, 28],
            [11, 14, 22, 36], [11, 14, 23, 35], [11, 14, 26, 33], [11, 15, 16, 18],
            [11, 15, 19, 36], [11, 15, 22, 31], [11, 15, 28, 33], [11, 15, 30, 34],
            [11, 16, 17, 35], [11, 16, 19, 37], [11, 16, 21, 31], [11, 16, 22, 24],
            [11, 16, 23, 25], [11, 16, 29, 30], [11, 16, 32, 36], [11, 17, 20, 23],
            [11, 17, 21, 29], [11, 17, 30, 36], [11, 17, 33, 34], [11, 18, 19, 22],
            [11, 18, 21, 32], [11, 18, 24, 31], [11, 18, 25, 33], [11, 18, 29, 35],
            [11, 18, 34, 37], [11, 19, 23, 29], [11, 19, 25, 26], [11, 19, 27, 32],
            [11, 19, 28, 30], [11, 20, 24, 34], [11, 20, 25, 29], [11, 20, 26, 36],
            [11, 20, 27, 37], [11, 20, 33, 35], [11, 21, 23, 26], [11, 21, 24, 37],
            [11, 21, 25, 27], [11, 22, 23, 27], [11, 22, 25, 35], [11, 22, 26, 37],
            [11, 22, 30, 32], [11, 23, 28, 31], [11, 23, 33, 36], [11, 24, 27, 35],
            [11, 24, 28, 32], [11, 25, 28, 37], [11, 25, 30, 31], [11, 26, 28, 34],
            [11, 26, 29, 31], [11, 29, 32, 37], [11, 30, 33, 37], [11, 31, 34, 35],
            [11, 35, 36, 37], [12, 13, 14, 26], [12, 13, 18, 29], [12, 13, 19, 24],
            [12, 13, 21, 22], [12, 13, 23, 30], [12, 13, 28, 34], [12, 13, 32, 37],
            [12, 13, 33, 35], [12, 14, 15, 22], [12, 14, 16, 21], [12, 14, 18, 27],
            [12, 14, 20, 35], [12, 14, 23, 34], [12, 14, 24, 25], [12, 14, 28, 37],
            [12, 14, 29, 30], [12, 15, 16, 30], [12, 15, 17, 35], [12, 15, 19, 31],
            [12, 15, 20, 32], [12, 15, 21, 29], [12, 15, 23, 37], [12, 15, 24, 36],
            [12, 15, 27, 34], [12, 16, 17, 19], [12, 16, 20, 37], [12, 16, 23, 32],
            [12, 16, 29, 34], [12, 16, 31, 35], [12, 17, 18, 36], [12, 17, 22, 32],
            [12, 17, 23, 25], [12, 17, 24, 26], [12, 17, 30, 31], [12, 17, 33, 37],
            [12, 18, 21, 24], [12, 18, 22, 30], [12, 18, 31, 37], [12, 18, 34, 35],
            [12, 19, 20, 23], [12, 19, 22, 33], [12, 19, 25, 32], [12, 19, 26, 34],
            [12, 19, 30, 36], [12, 20, 24, 30], [12, 20, 26, 27], [12, 20, 28, 33],
            [12, 20, 29, 31], [12, 21, 25, 35], [12, 21, 26, 30], [12, 21, 27, 37],
            [12, 21, 34, 36], [12, 22, 24, 27], [12, 22, 26, 28], [12, 23, 24, 28],
            [12, 23, 26, 36], [12, 23, 31, 33], [12, 24, 29, 32], [12, 24, 34, 37],
            [12, 25, 28, 36], [12, 25, 29, 33], [12, 26, 31, 32], [12, 27, 29, 35],
            [12, 27, 30, 32], [12, 32, 35, 36], [13, 14, 15, 27], [13, 14, 19, 30],
            [13, 14, 20, 25], [13, 14, 22, 23], [13, 14, 24, 31], [13, 14, 29, 35],
            [13, 14, 34, 36], [13, 15, 16, 23], [13, 15, 17, 22], [13, 15, 19, 28],
            [13, 15, 21, 36], [13, 15, 24, 35], [13, 15, 25, 26], [13, 15, 30, 31],
            [13, 16, 17, 31], [13, 16, 18, 36], [13, 16, 20, 32], [13, 16, 21, 33],
            [13, 16, 22, 30], [13, 16, 25, 37], [13, 16, 28, 35], [13, 17, 18, 20],
            [13, 17, 24, 33], [13, 17, 30, 35], [13, 17, 32, 36], [13, 18, 19, 37],
            [13, 18, 23, 33], [13, 18, 24, 26], [13, 18, 25, 27], [13, 18, 31, 32],
            [13, 19, 22, 25], [13, 19, 23, 31], [13, 19, 35, 36], [13, 20, 21, 24],
            [13, 20, 23, 34], [13, 20, 26, 33], [13, 20, 27, 35], [13, 20, 31, 37],
            [13, 21, 25, 31], [13, 21, 27, 28], [13, 21, 29, 34], [13, 21, 30, 32],
            [13, 22, 26, 36], [13, 22, 27, 31], [13, 22, 35, 37], [13, 23, 25, 28],
            [13, 23, 27, 29], [13, 24, 25, 29], [13, 24, 27, 37], [13, 24, 32, 34],
            [13, 25, 30, 33], [13, 26, 29, 37], [13, 26, 30, 34], [13, 27, 32, 33],
            [13, 28, 30, 36], [13, 28, 31, 33], [13, 33, 36, 37], [14, 15, 16, 28],
            [14, 15, 20, 31], [14, 15, 21, 26], [14, 15, 23, 24], [14, 15, 25, 32],
            [14, 15, 30, 36], [14, 15, 35, 37], [14, 16, 17, 24], [14, 16, 18, 23],
            [14, 16, 20, 29], [14, 16, 22, 37], [14, 16, 25, 36], [14, 16, 26, 27],
            [14, 16, 31, 32], [14, 17, 18, 32], [14, 17, 19, 37], [14, 17, 21, 33],
            [14, 17, 22, 34], [14, 17, 23, 31], [14, 17, 29, 36], [14, 18, 19, 21],
            [14, 18, 25, 34], [14, 18, 31, 36], [14, 18, 33, 37], [14, 19, 24, 34],
            [14, 19, 25, 27], [14, 19, 26, 28], [14, 19, 32, 33], [14, 20, 23, 26],
            [14, 20, 24, 32], [14, 20, 36, 37], [14, 21, 22, 25], [14, 21, 24, 35],
            [14, 21, 27, 34], [14, 21, 28, 36], [14, 22, 26, 32], [14, 22, 28, 29],
            [14, 22, 30, 35], [14, 22, 31, 33], [14, 23, 27, 37], [14, 23, 28, 32],
            [14, 24, 26, 29], [14, 24, 28, 30], [14, 25, 26, 30], [14, 25, 33, 35],
            [14, 26, 31, 34], [14, 27, 31, 35], [14, 28, 33, 34], [14, 29, 31, 37],
            [14, 29, 32, 34], [15, 16, 17, 29], [15, 16, 21, 32], [15, 16, 22, 27],
            [15, 16, 24, 25], [15, 16, 26, 33], [15, 16, 31, 37], [15, 17, 18, 25],
            [15, 17, 19, 24], [15, 17, 21, 30], [15, 17, 26, 37], [15, 17, 27, 28],
            [15, 17, 32, 33], [15, 18, 19, 33], [15, 18, 22, 34], [15, 18, 23, 35],
            [15, 18, 24, 32], [15, 18, 30, 37], [15, 19, 20, 22], [15, 19, 26, 35],
            [15, 19, 32, 37], [15, 20, 25, 35], [15, 20, 26, 28], [15, 20, 27, 29],
            [15, 20, 33, 34], [15, 21, 24, 27], [15, 21, 25, 33], [15, 22, 23, 26],
            [15, 22, 25, 36], [15, 22, 28, 35], [15, 22, 29, 37], [15, 23, 27, 33],
            [15, 23, 29, 30], [15, 23, 31, 36], [15, 23, 32, 34], [15, 24, 29, 33],
            [15, 25, 27, 30], [15, 25, 29, 31], [15, 26, 27, 31], [15, 26, 34, 36],
            [15, 27, 32, 35], [15, 28, 32, 36], [15, 29, 34, 35], [15, 30, 33, 35],
            [16, 17, 18, 30], [16, 17, 22, 33], [16, 17, 23, 28], [16, 17, 25, 26],
            [16, 17, 27, 34], [16, 18, 19, 26], [16, 18, 20, 25], [16, 18, 22, 31],
            [16, 18, 28, 29], [16, 18, 33, 34], [16, 19, 20, 34], [16, 19, 23, 35],
            [16, 19, 24, 36], [16, 19, 25, 33], [16, 20, 21, 23], [16, 20, 27, 36],
            [16, 21, 26, 36], [16, 21, 27, 29], [16, 21, 28, 30], [16, 21, 34, 35],
            [16, 22, 25, 28], [16, 22, 26, 34], [16, 23, 24, 27], [16, 23, 26, 37],
            [16, 23, 29, 36], [16, 24, 28, 34], [16, 24, 30, 31], [16, 24, 32, 37],
            [16, 24, 33, 35], [16, 25, 30, 34], [16, 26, 28, 31], [16, 26, 30, 32],
            [16, 27, 28, 32], [16, 27, 35, 37], [16, 28, 33, 36], [16, 29, 33, 37],
            [16, 30, 35, 36], [16, 31, 34, 36], [17, 18, 19, 31], [17, 18, 23, 34],
            [17, 18, 24, 29], [17, 18, 26, 27], [17, 18, 28, 35], [17, 19, 20, 27],
            [17, 19, 21, 26], [17, 19, 23, 32], [17, 19, 29, 30], [17, 19, 34, 35],
            [17, 20, 21, 35], [17, 20, 24, 36], [17, 20, 25, 37], [17, 20, 26, 34],
            [17, 21, 22, 24], [17, 21, 28, 37], [17, 22, 27, 37], [17, 22, 28, 30],
            [17, 22, 29, 31], [17, 22, 35, 36], [17, 23, 26, 29], [17, 23, 27, 35],
            [17, 24, 25, 28], [17, 24, 30, 37], [17, 25, 29, 35], [17, 25, 31, 32],
            [17, 25, 34, 36], [17, 26, 31, 35], [17, 27, 29, 32], [17, 27, 31, 33],
            [17, 28, 29, 33], [17, 29, 34, 37], [17, 31, 36, 37], [17, 32, 35, 37],
            [18, 19, 20, 32], [18, 19, 24, 35], [18, 19, 25, 30], [18, 19, 27, 28],
            [18, 19, 29, 36], [18, 20, 21, 28], [18, 20, 22, 27], [18, 20, 24, 33],
            [18, 20, 30, 31], [18, 20, 35, 36], [18, 21, 22, 36], [18, 21, 25, 37],
            [18, 21, 27, 35], [18, 22, 23, 25], [18, 23, 29, 31], [18, 23, 30, 32],
            [18, 23, 36, 37], [18, 24, 27, 30], [18, 24, 28, 36], [18, 25, 26, 29],
            [18, 26, 30, 36], [18, 26, 32, 33], [18, 26, 35, 37], [18, 27, 32, 36],
            [18, 28, 30, 33], [18, 28, 32, 34], [18, 29, 30, 34], [19, 20, 21, 33],
            [19, 20, 25, 36], [19, 20, 26, 31], [19, 20, 28, 29], [19, 20, 30, 37],
            [19, 21, 22, 29], [19, 21, 23, 28], [19, 21, 25, 34], [19, 21, 31, 32],
            [19, 21, 36, 37], [19, 22, 23, 37], [19, 22, 28, 36], [19, 23, 24, 26],
            [19, 24, 30, 32], [19, 24, 31, 33], [19, 25, 28, 31], [19, 25, 29, 37],
            [19, 26, 27, 30], [19, 27, 31, 37], [19, 27, 33, 34], [19, 28, 33, 37],
            [19, 29, 31, 34], [19, 29, 33, 35], [19, 30, 31, 35], [20, 21, 22, 34],
            [20, 21, 26, 37], [20, 21, 27, 32], [20, 21, 29, 30], [20, 22, 23, 30],
            [20, 22, 24, 29], [20, 22, 26, 35], [20, 22, 32, 33], [20, 23, 29, 37],
            [20, 24, 25, 27], [20, 25, 31, 33], [20, 25, 32, 34], [20, 26, 29, 32],
            [20, 27, 28, 31], [20, 28, 34, 35], [20, 30, 32, 35], [20, 30, 34, 36],
            [20, 31, 32, 36], [21, 22, 23, 35], [21, 22, 28, 33], [21, 22, 30, 31],
            [21, 23, 24, 31], [21, 23, 25, 30], [21, 23, 27, 36], [21, 23, 33, 34],
            [21, 25, 26, 28], [21, 26, 32, 34], [21, 26, 33, 35], [21, 27, 30, 33],
            [21, 28, 29, 32], [21, 29, 35, 36], [21, 31, 33, 36], [21, 31, 35, 37],
            [21, 32, 33, 37], [22, 23, 24, 36], [22, 23, 29, 34], [22, 23, 31, 32],
            [22, 24, 25, 32], [22, 24, 26, 31], [22, 24, 28, 37], [22, 24, 34, 35],
            [22, 26, 27, 29], [22, 27, 33, 35], [22, 27, 34, 36], [22, 28, 31, 34],
            [22, 29, 30, 33], [22, 30, 36, 37], [22, 32, 34, 37], [23, 24, 25, 37],
            [23, 24, 30, 35], [23, 24, 32, 33], [23, 25, 26, 33], [23, 25, 27, 32],
            [23, 25, 35, 36], [23, 27, 28, 30], [23, 28, 34, 36], [23, 28, 35, 37],
            [23, 29, 32, 35], [23, 30, 31, 34], [24, 25, 31, 36], [24, 25, 33, 34],
            [24, 26, 27, 34], [24, 26, 28, 33], [24, 26, 36, 37], [24, 28, 29, 31],
            [24, 29, 35, 37], [24, 30, 33, 36], [24, 31, 32, 35], [25, 26, 32, 37],
            [25, 26, 34, 35], [25, 27, 28, 35], [25, 27, 29, 34], [25, 29, 30, 32],
            [25, 31, 34, 37], [25, 32, 33, 36], [26, 27, 35, 36], [26, 28, 29, 36],
            [26, 28, 30, 35], [26, 30, 31, 33], [26, 33, 34, 37], [27, 28, 36, 37],
            [27, 29, 30, 37], [27, 29, 31, 36], [27, 31, 32, 34], [28, 30, 32, 37],
            [28, 32, 33, 35], [29, 33, 34, 36], [30, 34, 35, 37]]
