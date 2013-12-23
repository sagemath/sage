r"""
Linear code constructions

This file contains constructions of error-correcting codes which are
pure Python/Sage and not obtained from wrapping GUAVA functions.
The GUAVA wrappers are in guava.py.

All codes available here can be accessed through the ``codes`` object::

    sage: codes.HammingCode(3,GF(2))
    Linear code of length 7, dimension 4 over Finite Field of size 2

Let `F` be a finite field with `q` elements.
Here's a constructive definition of a cyclic code of length
`n`.

#. Pick a monic polynomial `g(x)\in F[x]` dividing
   `x^n-1`. This is called the generating polynomial of the
   code.

#. For each polynomial `p(x)\in F[x]`, compute
   `p(x)g(x)\ ({\rm mod}\ x^n-1)`. Denote the answer by
   `c_0+c_1x+...+c_{n-1}x^{n-1}`.

#. `{\bf c} =(c_0,c_1,...,c_{n-1})` is a codeword in
   `C`. Every codeword in `C` arises in this way
   (from some `p(x)`).

The polynomial notation for the code is to call
`c_0+c_1x+...+c_{n-1}x^{n-1}` the codeword (instead of
`(c_0,c_1,...,c_{n-1})`). The polynomial
`h(x)=(x^n-1)/g(x)` is called the check polynomial of
`C`.

Let `n` be a positive integer relatively prime to
`q` and let `\alpha` be a primitive
`n`-th root of unity. Each generator polynomial `g`
of a cyclic code `C` of length `n` has a
factorization of the form


.. math::

     g(x) = (x - \alpha^{k_1})...(x - \alpha^{k_r}),


where `\{k_1,...,k_r\} \subset \{0,...,n-1\}`. The
numbers `\alpha^{k_i}`, `1 \leq i \leq r`, are
called the zeros of the code `C`. Many families of cyclic
codes (such as BCH codes and the quadratic residue codes) are
defined using properties of the zeros of `C`.

- BCHCode - A 'Bose-Chaudhuri-Hockenghem code' (or BCH code for short)
  is the largest possible cyclic code of length n over field F=GF(q),
  whose generator polynomial has zeros (which contain the set)
  `Z = \{a^i\ |\ i \in C_b\cup ... C_{b+delta-2}\}`, where
  `a` is a primitive `n^{th}` root of unity in the
  splitting field `GF(q^m)`, `b` is an integer
  `0\leq b\leq n-delta+1` and `m` is the multiplicative
  order of `q` modulo `n`. The default here is `b=0`
  (unlike Guava, which has default `b=1`). Here `C_k` are
  the cyclotomic codes (see ``cyclotomic_cosets``).

- BinaryGolayCode, ExtendedBinaryGolayCode, TernaryGolayCode,
  ExtendedTernaryGolayCode the well-known"extremal" Golay codes,
  http://en.wikipedia.org/wiki/Golay_code

- cyclic codes - CyclicCodeFromGeneratingPolynomial (= CyclicCode),
  CyclicCodeFromCheckPolynomial,
  http://en.wikipedia.org/wiki/Cyclic_code

- DuadicCodeEvenPair, DuadicCodeOddPair: Constructs the "even
  (resp. odd) pair" of duadic codes associated to the "splitting" S1,
  S2 of n. This is a special type of cyclic code whose generator is
  determined by S1, S2. See chapter 6 in [HP]_.

- HammingCode - the well-known Hamming code,
  http://en.wikipedia.org/wiki/Hamming_code

- LinearCodeFromCheckMatrix - for specifying the code using the
  check matrix instead of the generator matrix.

- QuadraticResidueCodeEvenPair, QuadraticResidueCodeOddPair: Quadratic
  residue codes of a given odd prime length and base ring either don't
  exist at all or occur as 4-tuples - a pair of
  "odd-like" codes and a pair of "even-like" codes. If n 2 is prime
  then (Theorem 6.6.2 in [HP]_) a QR code exists over GF(q) iff q is a
  quadratic residue mod n. Here they are constructed as"even-like"
  duadic codes associated the splitting (Q,N) mod n, where Q is the
  set of non-zero quadratic residues and N is the non-residues.
  QuadraticResidueCode (a special case) and
  ExtendedQuadraticResidueCode are included as well.

- RandomLinearCode - Repeatedly applies Sage's random_element applied
  to the ambient MatrixSpace of the generator matrix until a full rank
  matrix is found.

- ReedSolomonCode - Given a finite field `F` of order `q`,
  let `n` and `k` be chosen such that
  `1 \leq k \leq n \leq q`.
  Pick `n` distinct elements of `F`, denoted
  `\{ x_1, x_2, ... , x_n \}`. Then, the codewords are obtained
  by evaluating every polynomial in `F[x]` of degree less than
  `k` at each `x_i`.

- ToricCode - Let `P` denote a list of lattice points in
  `\ZZ^d` and let `T` denote a listing of all
  points in `(F^x )^d`. Put `n=|T|` and let `k`
  denote the dimension of the vector space of functions
  `V = \mathrm{Span} \{x^e \ |\ e \in P\}`.
  The associated toric code `C` is
  the evaluation code which is the image of the evaluation map
  `eval_T : V \rightarrow F^n`, where `x^e` is the
  multi-index notation.

- WalshCode - a binary linear `[2^m,m,2^{m-1}]` code
  related to Hadamard matrices.
  http://en.wikipedia.org/wiki/Walsh_code

REFERENCES:

.. [HP] W. C. Huffman, V. Pless, Fundamentals of Error-Correcting
   Codes, Cambridge Univ. Press, 2003.

AUTHOR:

- David Joyner (2007-05): initial version

- " (2008-02): added cyclic codes, Hamming codes

- " (2008-03): added BCH code, LinearCodeFromCheckmatrix, ReedSolomonCode, WalshCode,
  DuadicCodeEvenPair, DuadicCodeOddPair, QR codes (even and odd)

- " (2008-09) fix for bug in BCHCode reported by F. Voloch

- " (2008-10) small docstring changes to WalshCode and walsh_matrix

Functions
---------

"""
############################################################################
## Copyright David Joyner, 2007. wdjoyner@gmail.com.
##  This is released under the GPL, version 2 or later (www.fsf.org).
#############################################################################


from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix
from sage.rings.finite_rings.constructor import FiniteField as GF
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.misc.misc import prod
from linear_code import LinearCodeFromVectorSpace, LinearCode
from sage.modules.free_module import span
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.structure.sequence import Sequence
from sage.rings.arith import GCD,LCM,divisors,quadratic_residues
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer import Integer
from sage.sets.set import Set
from sage.rings.finite_rings.integer_mod import Mod

############### utility functions ################


def cyclotomic_cosets(q, n, t = None):
    r"""
    INPUT: q,n,t positive integers (or t=None) Some type-checking of
    inputs is performed.

    OUTPUT: q-cyclotomic cosets mod n (or, if t is not None, the q-cyclotomic
    coset mod n containing t)

    Let q, n be relatively print positive integers and let
    `A = q^{ZZ}`. The group A acts on ZZ/nZZ by multiplication.
    The orbits of this action are "cyclotomic cosets", or more
    precisely "q-cyclotomic cosets mod n". Sometimes the smallest
    element of the coset is called the "coset leader". The algorithm
    will always return the cosets as sorted lists of lists, so the
    coset leader will always be the first element in the list.

    These cosets arise in the theory of duadic codes and minimal
    polynomials of finite fields. Fix a primitive element `z`
    of `GF(q^k)`. The minimal polynomial of `z^s` over
    `GF(q)` is given by

    .. math::

             M_s(x) = \prod_{i \in C_s} (x-z^i),


    where `C_s` is the q-cyclotomic coset mod n containing s,
    `n = q^k - 1`.

    EXAMPLES::

        sage: cyclotomic_cosets(2,11)
        [[0], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
        sage: cyclotomic_cosets(2,15)
        [[0], [1, 2, 4, 8], [3, 6, 9, 12], [5, 10], [7, 11, 13, 14]]
        sage: cyclotomic_cosets(2,15,5)
        [5, 10]
        sage: cyclotomic_cosets(3,16)
        [[0], [1, 3, 9, 11], [2, 6], [4, 12], [5, 7, 13, 15], [8], [10, 14]]
        sage: F.<z> = GF(2^4, "z")
        sage: P.<x> = PolynomialRing(F,"x")
        sage: a = z^5
        sage: a.minimal_polynomial()
        x^2 + x + 1
        sage: prod([x-z^i for i in [5, 10]])
        x^2 + x + 1
        sage: cyclotomic_cosets(3,2,0)
        [0]
        sage: cyclotomic_cosets(3,2,1)
        [1]
        sage: cyclotomic_cosets(3,2,2)
        [0]

    This last output looks strange but is correct, since the elements of
    the cosets are in ZZ/nZZ and 2 = 0 in ZZ/2ZZ.
    """
    from sage.misc.misc import srange
    if not(t==None) and type(t)<>Integer:
        raise TypeError,  "Optional input %s must None or an integer."%t
    if q<2 or n<2:
        raise TypeError,  "Inputs %s and %s must be > 1."%(q,n)
    if GCD(q,n) <> 1:
        raise TypeError,  "Inputs %s and %s must be relative prime."%(q,n)
    if t<>None and type(t)==Integer:
        S = Set([t*q**i%n for i in srange(n)])
        L = list(S)
        L.sort()
        return L
    ccs = Set([])
    ccs_list = [[0]]
    for s in range(1,n):
        if not(s in ccs):
            S = Set([s*q**i%n for i in srange(n)])
            L = list(S)
            L.sort()
            ccs = ccs.union(S)
            ccs_list.append(L)
    return ccs_list

def is_a_splitting(S1,S2,n):
    """
    INPUT: S1, S2 are disjoint sublists partitioning [1, 2, ..., n-1]
    n1 is an integer

    OUTPUT: a, b where a is True or False, depending on whether S1, S2
    form a "splitting" of n (ie, if there is a b1 such that b\*S1=S2
    (point-wise multiplication mod n), and b is a splitting (if a =
    True) or 0 (if a = False)

    Splittings are useful for computing idempotents in the quotient
    ring `Q = GF(q)[x]/(x^n-1)`. For

    EXAMPLES::

        sage: from sage.coding.code_constructions import is_a_splitting
        sage: n = 11; q = 3
        sage: C = cyclotomic_cosets(q,n); C
        [[0], [1, 3, 4, 5, 9], [2, 6, 7, 8, 10]]
        sage: S1 = C[1]
        sage: S2 = C[2]
        sage: is_a_splitting(S1,S2,11)
        (True, 2)
        sage: F = GF(q)
        sage: P.<x> = PolynomialRing(F,"x")
        sage: I = Ideal(P,[x^n-1])
        sage: Q.<x> = QuotientRing(P,I)
        sage: i1 = -sum([x^i for i in S1]); i1
        2*x^9 + 2*x^5 + 2*x^4 + 2*x^3 + 2*x
        sage: i2 = -sum([x^i for i in S2]); i2
        2*x^10 + 2*x^8 + 2*x^7 + 2*x^6 + 2*x^2
        sage: i1^2 == i1
        True
        sage: i2^2 == i2
        True
        sage: (1-i1)^2 == 1-i1
        True
        sage: (1-i2)^2 == 1-i2
        True

    We return to dealing with polynomials (rather than elements of
    quotient rings), so we can construct cyclic codes::

        sage: P.<x> = PolynomialRing(F,"x")
        sage: i1 = -sum([x^i for i in S1])
        sage: i2 = -sum([x^i for i in S2])
        sage: i1_sqrd = (i1^2).quo_rem(x^n-1)[1]
        sage: i1_sqrd  == i1
        True
        sage: i2_sqrd = (i2^2).quo_rem(x^n-1)[1]
        sage: i2_sqrd  == i2
        True
        sage: C1 = codes.CyclicCodeFromGeneratingPolynomial(n,i1)
        sage: C2 = codes.CyclicCodeFromGeneratingPolynomial(n,1-i2)
        sage: C1.dual_code() == C2
        True

    This is a special case of Theorem 6.4.3 in [HP]_.
    """
    if Set(S1).union(Set(S2)) <> Set(range(1,n)):
        raise TypeError, "Lists must partition [1,2,...,n-1]."
    if n<3:
        raise TypeError,  "Input %s must be > 2."%n
    for b in range(2,n):
        SS1 = Set([b*x%n for x in S1])
        SS2 = Set([b*x%n for x in S2])
        if SS1 == Set(S2) and SS2 == Set(S1):
            return True, b
    return False, 0


def lift2smallest_field(a):
    """
    INPUT: a is an element of a finite field GF(q)

    OUTPUT: the element b of the smallest subfield F of GF(q) for
    which F(b)=a.

    EXAMPLES::

        sage: from sage.coding.code_constructions import lift2smallest_field
        sage: FF.<z> = GF(3^4,"z")
        sage: a = z^10
        sage: lift2smallest_field(a)
        (2*z + 1, Finite Field in z of size 3^2)
        sage: a = z^40
        sage: lift2smallest_field(a)
        (2, Finite Field of size 3)

    AUTHORS:

    - John Cremona
    """
    FF = a.parent()
    k = FF.degree()
    if k == 1:
        return a, FF
    pol = a.minimal_polynomial()
    d = pol.degree()
    if d == k:
        return a, FF
    p = FF.characteristic()
    F = GF(p**d,"z")
    b = pol.roots(F,multiplicities=False)[0]
    return b, F

def lift2smallest_field2(a):
    """
    INPUT: a is an element of a finite field GF(q) OUTPUT: the element
    b of the smallest subfield F of GF(q) for which F(b)=a.

    EXAMPLES::

        sage: from sage.coding.code_constructions import lift2smallest_field2
        sage: FF.<z> = GF(3^4,"z")
        sage: a = z^40
        sage: lift2smallest_field2(a)
        (2, Finite Field of size 3)
        sage: FF.<z> = GF(2^4,"z")
        sage: a = z^15
        sage: lift2smallest_field2(a)
        (1, Finite Field of size 2)

    .. warning::

       Since coercion (the FF(b) step) has a bug in it, this
       *only works* in the case when you *know* F is a prime field.

    AUTHORS:

    - David Joyner
    """
    FF = a.parent()
    q = FF.order()
    if q.is_prime():
        return a,FF
    p = q.factor()[0][0]
    k = q.factor()[0][1]
    for d in divisors(k):
        F = GF(p**d,"zz")
        for b in F:
            if FF(b) == a:
                return b, F


def permutation_action(g,v):
    """
    Returns permutation of rows g\*v. Works on lists, matrices,
    sequences and vectors (by permuting coordinates). The code requires
    switching from i to i+1 (and back again) since the SymmetricGroup
    is, by convention, the symmetric group on the "letters" 1, 2, ...,
    n (not 0, 1, ..., n-1).

    EXAMPLES::

        sage: V = VectorSpace(GF(3),5)
        sage: v = V([0,1,2,0,1])
        sage: G = SymmetricGroup(5)
        sage: g = G([(1,2,3)])
        sage: permutation_action(g,v)
        (1, 2, 0, 0, 1)
        sage: g = G([()])
        sage: permutation_action(g,v)
        (0, 1, 2, 0, 1)
        sage: g = G([(1,2,3,4,5)])
        sage: permutation_action(g,v)
        (1, 2, 0, 1, 0)
        sage: L = Sequence([1,2,3,4,5])
        sage: permutation_action(g,L)
        [2, 3, 4, 5, 1]
        sage: MS = MatrixSpace(GF(3),3,7)
        sage: A = MS([[1,0,0,0,1,1,0],[0,1,0,1,0,1,0],[0,0,0,0,0,0,1]])
        sage: S5 = SymmetricGroup(5)
        sage: g = S5([(1,2,3)])
        sage: A
        [1 0 0 0 1 1 0]
        [0 1 0 1 0 1 0]
        [0 0 0 0 0 0 1]
        sage: permutation_action(g,A)
        [0 1 0 1 0 1 0]
        [0 0 0 0 0 0 1]
        [1 0 0 0 1 1 0]

    It also works on lists and is a "left action"::

        sage: v = [0,1,2,0,1]
        sage: G = SymmetricGroup(5)
        sage: g = G([(1,2,3)])
        sage: gv = permutation_action(g,v); gv
        [1, 2, 0, 0, 1]
        sage: permutation_action(g,v) == g(v)
        True
        sage: h = G([(3,4)])
        sage: gv = permutation_action(g,v)
        sage: hgv = permutation_action(h,gv)
        sage: hgv == permutation_action(h*g,v)
        True

    AUTHORS:

    - David Joyner, licensed under the GPL v2 or greater.
    """
    v_type_list = False
    if type(v) == list:
        v_type_list = True
        v = Sequence(v)
    V = v.parent()
    n = len(list(v))
    gv = []
    for i in range(n):
        gv.append(v[g(i+1)-1])
    if v_type_list:
        return gv
    return V(gv)

def walsh_matrix(m0):
    """
    This is the generator matrix of a Walsh code. The matrix of
    codewords correspond to a Hadamard matrix.

    EXAMPLES::

        sage: walsh_matrix(2)
        [0 0 1 1]
        [0 1 0 1]
        sage: walsh_matrix(3)
        [0 0 0 0 1 1 1 1]
        [0 0 1 1 0 0 1 1]
        [0 1 0 1 0 1 0 1]
        sage: C = LinearCode(walsh_matrix(4)); C
        Linear code of length 16, dimension 4 over Finite Field of size 2
        sage: C.spectrum()
        [1, 0, 0, 0, 0, 0, 0, 0, 15, 0, 0, 0, 0, 0, 0, 0, 0]

    This last code has minimum distance 8.

    REFERENCES:

    - http://en.wikipedia.org/wiki/Hadamard_matrix
    """
    m = int(m0)
    if m == 1:
        return matrix(GF(2), 1, 2, [ 0, 1])
    if m > 1:
        row2 = [x.list() for x in walsh_matrix(m-1).augment(walsh_matrix(m-1)).rows()]
        return matrix(GF(2), m, 2**m, [[0]*2**(m-1) + [1]*2**(m-1)] + row2)
    raise ValueError, "%s must be an integer > 0."%m0


##################### main constructions #####################



def BCHCode(n,delta,F,b=0):
    r"""
    A 'Bose-Chaudhuri-Hockenghem code' (or BCH code for short) is the
    largest possible cyclic code of length n over field F=GF(q), whose
    generator polynomial has zeros (which contain the set)
    `Z = \{a^{b},a^{b+1}, ..., a^{b+delta-2}\}`, where a is a
    primitive `n^{th}` root of unity in the splitting field
    `GF(q^m)`, b is an integer `0\leq b\leq n-delta+1`
    and m is the multiplicative order of q modulo n. (The integers
    `b,...,b+delta-2` typically lie in the range
    `1,...,n-1`.) The integer `delta \geq 1` is called
    the "designed distance". The length n of the code and the size q of
    the base field must be relatively prime. The generator polynomial
    is equal to the least common multiple of the minimal polynomials of
    the elements of the set `Z` above.

    Special cases are b=1 (resulting codes are called 'narrow-sense'
    BCH codes), and `n=q^m-1` (known as 'primitive' BCH
    codes).

    It may happen that several values of delta give rise to the same
    BCH code. The largest one is called the Bose distance of the code.
    The true minimum distance, d, of the code is greater than or equal
    to the Bose distance, so `d\geq delta`.

    EXAMPLES::

        sage: FF.<a> = GF(3^2,"a")
        sage: x = PolynomialRing(FF,"x").gen()
        sage: L = [b.minpoly() for b in [a,a^2,a^3]]; g = LCM(L)
        sage: f = x^(8)-1
        sage: g.divides(f)
        True
        sage: C = codes.CyclicCode(8,g); C
        Linear code of length 8, dimension 4 over Finite Field of size 3
        sage: C.minimum_distance()
        4
        sage: C = codes.BCHCode(8,3,GF(3),1); C
        Linear code of length 8, dimension 4 over Finite Field of size 3
        sage: C.minimum_distance()
        4
        sage: C = codes.BCHCode(8,3,GF(3)); C
        Linear code of length 8, dimension 5 over Finite Field of size 3
        sage: C.minimum_distance()
        3
        sage: C = codes.BCHCode(26, 5, GF(5), b=1); C
        Linear code of length 26, dimension 10 over Finite Field of size 5

    """
    from sage.misc.misc import srange
    q = F.order()
    R = IntegerModRing(n)
    m = R(q).multiplicative_order()
    FF = GF(q**m,"z")
    z = FF.gen()
    e = z.multiplicative_order()/n
    a = z**e # order n
    P = PolynomialRing(F,"x")
    x = P.gen()
    cosets = Set([])
    for i in srange(b,b+delta-1):
        cosets = cosets.union(Set(cyclotomic_cosets(q, n, i)))
    L0 = [a**j for j in cosets]
    L1 = [P(ai.minpoly()) for ai in L0]
    g = P(LCM(L1))
    #print cosets, "\n", g, "\n", (x**n-1).factor(), "\n", L1, "\n", g.divides(x**n-1)
    if not(g.divides(x**n-1)):
        raise ValueError, "BCH codes does not exist with the given input."
    return CyclicCodeFromGeneratingPolynomial(n,g)


def BinaryGolayCode():
    r"""
    BinaryGolayCode() returns a binary Golay code. This is a perfect
    [23,12,7] code. It is also (equivalent to) a cyclic code, with
    generator polynomial
    `g(x)=1+x^2+x^4+x^5+x^6+x^{10}+x^{11}`. Extending it yields
    the extended Golay code (see ExtendedBinaryGolayCode).

    EXAMPLE::

        sage: C = codes.BinaryGolayCode()
        sage: C
        Linear code of length 23, dimension 12 over Finite Field of size 2
        sage: C.minimum_distance()
        7
        sage: C.minimum_distance(algorithm='gap') # long time, check d=7
        7

    AUTHORS:

    - David Joyner (2007-05)
    """
    F = GF(2)
    B = [[1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
          [0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
          [0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
          [0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],\
          [0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],\
          [0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],\
          [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],\
          [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0],\
          [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0],\
          [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0],\
          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0],\
          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1]]
    # MS = MatrixSpace(F,12,23)
    # V = VectorSpace(F,23)
    V = span(B, F)
    return LinearCodeFromVectorSpace(V, d=7)


def CyclicCodeFromGeneratingPolynomial(n,g,ignore=True):
    r"""
    If g is a polynomial over GF(q) which divides `x^n-1` then
    this constructs the code "generated by g" (ie, the code associated
    with the principle ideal `gR` in the ring
    `R = GF(q)[x]/(x^n-1)` in the usual way).

    The option "ignore" says to ignore the condition that (a) the
    characteristic of the base field does not divide the length (the
    usual assumption in the theory of cyclic codes), and (b) `g`
    must divide `x^n-1`. If ignore=True, instead of returning
    an error, a code generated by `gcd(x^n-1,g)` is created.

    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(3),"x")
        sage: g = x-1
        sage: C = codes.CyclicCodeFromGeneratingPolynomial(4,g); C
        Linear code of length 4, dimension 3 over Finite Field of size 3
        sage: P.<x> = PolynomialRing(GF(4,"a"),"x")
        sage: g = x^3+1
        sage: C = codes.CyclicCodeFromGeneratingPolynomial(9,g); C
        Linear code of length 9, dimension 6 over Finite Field in a of size 2^2
        sage: P.<x> = PolynomialRing(GF(2),"x")
        sage: g = x^3+x+1
        sage: C = codes.CyclicCodeFromGeneratingPolynomial(7,g); C
        Linear code of length 7, dimension 4 over Finite Field of size 2
        sage: C.gen_mat()
        [1 1 0 1 0 0 0]
        [0 1 1 0 1 0 0]
        [0 0 1 1 0 1 0]
        [0 0 0 1 1 0 1]
        sage: g = x+1
        sage: C = codes.CyclicCodeFromGeneratingPolynomial(4,g); C
        Linear code of length 4, dimension 3 over Finite Field of size 2
        sage: C.gen_mat()
        [1 1 0 0]
        [0 1 1 0]
        [0 0 1 1]

    On the other hand, CyclicCodeFromPolynomial(4,x) will produce a
    ValueError including a traceback error message: "`x` must
    divide `x^4 - 1`". You will also get a ValueError if you
    type

    ::

        sage: P.<x> = PolynomialRing(GF(4,"a"),"x")
        sage: g = x^2+1

    followed by CyclicCodeFromGeneratingPolynomial(6,g). You will also
    get a ValueError if you type

    ::

        sage: P.<x> = PolynomialRing(GF(3),"x")
        sage: g = x^2-1
        sage: C = codes.CyclicCodeFromGeneratingPolynomial(5,g); C
        Linear code of length 5, dimension 4 over Finite Field of size 3

    followed by C = CyclicCodeFromGeneratingPolynomial(5,g,False), with
    a traceback message including "`x^2 + 2` must divide
    `x^5 - 1`".
    """
    P = g.parent()
    x = P.gen()
    F = g.base_ring()
    p = F.characteristic()
    if not(ignore) and p.divides(n):
        raise ValueError, 'The characteristic %s must not divide %s'%(p,n)
    if not(ignore) and not(g.divides(x**n-1)):
        raise ValueError, '%s must divide x^%s - 1'%(g,n)
    gn = GCD([g,x**n-1])
    d = gn.degree()
    coeffs = Sequence(gn.list())
    r1 = Sequence(coeffs+[0]*(n - d - 1))
    Sn = SymmetricGroup(n)
    s = Sn.gens()[0] # assumes 1st gen of S_n is (1,2,...,n)
    rows = [permutation_action(s**(-i),r1) for i in range(n-d)]
    MS = MatrixSpace(F,n-d,n)
    return LinearCode(MS(rows))

CyclicCode = CyclicCodeFromGeneratingPolynomial

def CyclicCodeFromCheckPolynomial(n,h,ignore=True):
    r"""
    If h is a polynomial over GF(q) which divides `x^n-1` then
    this constructs the code "generated by `g = (x^n-1)/h`"
    (ie, the code associated with the principle ideal `gR` in
    the ring `R = GF(q)[x]/(x^n-1)` in the usual way). The
    option "ignore" says to ignore the condition that the
    characteristic of the base field does not divide the length (the
    usual assumption in the theory of cyclic codes).

    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(3),"x")
        sage: C = codes.CyclicCodeFromCheckPolynomial(4,x + 1); C
        Linear code of length 4, dimension 1 over Finite Field of size 3
        sage: C = codes.CyclicCodeFromCheckPolynomial(4,x^3 + x^2 + x + 1); C
        Linear code of length 4, dimension 3 over Finite Field of size 3
        sage: C.gen_mat()
        [2 1 0 0]
        [0 2 1 0]
        [0 0 2 1]
    """
    P = h.parent()
    x = P.gen()
    d = h.degree()
    F = h.base_ring()
    p = F.characteristic()
    if not(ignore) and p.divides(n):
        raise ValueError, 'The characteristic %s must not divide %s'%(p,n)
    if not(h.divides(x**n-1)):
        raise ValueError, '%s must divide x^%s - 1'%(h,n)
    g = P((x**n-1)/h)
    return CyclicCodeFromGeneratingPolynomial(n,g)

def DuadicCodeEvenPair(F,S1,S2):
    r"""
    Constructs the "even pair" of duadic codes associated to the
    "splitting" (see the docstring for ``is_a_splitting``
    for the definition) S1, S2 of n.

    .. warning::

       Maybe the splitting should be associated to a sum of
       q-cyclotomic cosets mod n, where q is a *prime*.

    EXAMPLES::

        sage: from sage.coding.code_constructions import is_a_splitting
        sage: n = 11; q = 3
        sage: C = cyclotomic_cosets(q,n); C
        [[0], [1, 3, 4, 5, 9], [2, 6, 7, 8, 10]]
        sage: S1 = C[1]
        sage: S2 = C[2]
        sage: is_a_splitting(S1,S2,11)
        (True, 2)
        sage: codes.DuadicCodeEvenPair(GF(q),S1,S2)
        (Linear code of length 11, dimension 5 over Finite Field of size 3,
         Linear code of length 11, dimension 5 over Finite Field of size 3)
    """
    n = max(S1+S2)+1
    if not(is_a_splitting(S1,S2,n)):
        raise TypeError, "%s, %s must be a splitting of %s."%(S1,S2,n)
    q = F.order()
    k = Mod(q,n).multiplicative_order()
    FF = GF(q**k,"z")
    z = FF.gen()
    zeta = z**((q**k-1)/n)
    P1 = PolynomialRing(FF,"x")
    x = P1.gen()
    g1 = prod([x-zeta**i for i in S1+[0]])
    g2 = prod([x-zeta**i for i in S2+[0]])
    P2 = PolynomialRing(F,"x")
    x = P2.gen()
    gg1 = P2([lift2smallest_field(c)[0] for c in g1.coeffs()])
    gg2 = P2([lift2smallest_field(c)[0] for c in g2.coeffs()])
    C1 = CyclicCodeFromGeneratingPolynomial(n,gg1)
    C2 = CyclicCodeFromGeneratingPolynomial(n,gg2)
    return C1,C2

def DuadicCodeOddPair(F,S1,S2):
    """
    Constructs the "odd pair" of duadic codes associated to the
    "splitting" S1, S2 of n.

    .. warning::

       Maybe the splitting should be associated to a sum of
       q-cyclotomic cosets mod n, where q is a *prime*.

    EXAMPLES::

        sage: from sage.coding.code_constructions import is_a_splitting
        sage: n = 11; q = 3
        sage: C = cyclotomic_cosets(q,n); C
        [[0], [1, 3, 4, 5, 9], [2, 6, 7, 8, 10]]
        sage: S1 = C[1]
        sage: S2 = C[2]
        sage: is_a_splitting(S1,S2,11)
        (True, 2)
        sage: codes.DuadicCodeOddPair(GF(q),S1,S2)
        (Linear code of length 11, dimension 6 over Finite Field of size 3,
         Linear code of length 11, dimension 6 over Finite Field of size 3)

    This is consistent with Theorem 6.1.3 in [HP]_.
    """
    n = max(S1+S2)+1
    if not(is_a_splitting(S1,S2,n)):
        raise TypeError, "%s, %s must be a splitting of %s."%(S1,S2,n)
    q = F.order()
    k = Mod(q,n).multiplicative_order()
    FF = GF(q**k,"z")
    z = FF.gen()
    zeta = z**((q**k-1)/n)
    P1 = PolynomialRing(FF,"x")
    x = P1.gen()
    g1 = prod([x-zeta**i for i in S1+[0]])
    g2 = prod([x-zeta**i for i in S2+[0]])
    j = sum([x**i/n for i in range(n)])
    P2 = PolynomialRing(F,"x")
    x = P2.gen()
    coeffs1 = [lift2smallest_field(c)[0] for c in (g1+j).coeffs()]
    coeffs2 = [lift2smallest_field(c)[0] for c in (g2+j).coeffs()]
    gg1 = P2(coeffs1)
    gg2 = P2(coeffs2)
    C1 = CyclicCodeFromGeneratingPolynomial(n,gg1)
    C2 = CyclicCodeFromGeneratingPolynomial(n,gg2)
    return C1,C2


def ExtendedBinaryGolayCode():
    """
    ExtendedBinaryGolayCode() returns the extended binary Golay code.
    This is a perfect [24,12,8] code. This code is self-dual.

    EXAMPLES::

        sage: C = codes.ExtendedBinaryGolayCode()
        sage: C
        Linear code of length 24, dimension 12 over Finite Field of size 2
        sage: C.minimum_distance()
        8
        sage: C.minimum_distance(algorithm='gap') # long time, check d=8
        8

    AUTHORS:

    - David Joyner (2007-05)
    """
    B = [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1],\
         [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0],\
         [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1],\
         [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0],\
         [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1],\
         [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1],\
         [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1],\
         [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0],\
         [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0],\
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0],\
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1],\
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]]
    V = span(B, GF(2))
    return LinearCodeFromVectorSpace(V, d=8)
    # C = BinaryGolayCode()
    # return C.extended_code()


def ExtendedQuadraticResidueCode(n,F):
    r"""
    The extended quadratic residue code (or XQR code) is obtained from
    a QR code by adding a check bit to the last coordinate. (These
    codes have very remarkable properties such as large automorphism
    groups and duality properties - see [HP]_, Section 6.6.3-6.6.4.)

    INPUT:


    -  ``n`` - an odd prime

    -  ``F`` - a finite prime field F whose order must be a
       quadratic residue modulo n.


    OUTPUT: Returns an extended quadratic residue code.

    EXAMPLES::

        sage: C1 = codes.QuadraticResidueCode(7,GF(2))
        sage: C2 = C1.extended_code()
        sage: C3 = codes.ExtendedQuadraticResidueCode(7,GF(2)); C3
        Linear code of length 8, dimension 4 over Finite Field of size 2
        sage: C2 == C3
        True
        sage: C = codes.ExtendedQuadraticResidueCode(17,GF(2))
        sage: C
        Linear code of length 18, dimension 9 over Finite Field of size 2
        sage: C3 = codes.QuadraticResidueCodeOddPair(7,GF(2))[0]
        sage: C3x = C3.extended_code()
        sage: C4 = codes.ExtendedQuadraticResidueCode(7,GF(2))
        sage: C3x == C4
        True

    AUTHORS:

    - David Joyner (07-2006)
    """
    C = QuadraticResidueCodeOddPair(n,F)[0]
    return C.extended_code()

def ExtendedTernaryGolayCode():
    """
    ExtendedTernaryGolayCode returns a ternary Golay code. This is a
    self-dual perfect [12,6,6] code.

    EXAMPLES::

        sage: C = codes.ExtendedTernaryGolayCode()
        sage: C
        Linear code of length 12, dimension 6 over Finite Field of size 3
        sage: C.minimum_distance()
        6
        sage: C.minimum_distance(algorithm='gap') # long time, check d=6
        6

    AUTHORS:

    - David Joyner (11-2005)
    """
    B = [[1, 0, 0, 0, 0, 0, 2, 0, 1, 2, 1, 2],\
         [0, 1, 0, 0, 0, 0, 1, 2, 2, 2, 1, 0],\
         [0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1],\
         [0, 0, 0, 1, 0, 0, 1, 1, 0, 2, 2, 2],\
         [0, 0, 0, 0, 1, 0, 2, 1, 2, 2, 0, 1],\
         [0, 0, 0, 0, 0, 1, 0, 2, 1, 2, 2, 1]]
    V = span(B, GF(3))
    return LinearCodeFromVectorSpace(V, d=6)
    # C = TernaryGolayCode()
    # return C.extended_code()

def HammingCode(r,F):
    r"""
    Implements the Hamming codes.

    The `r^{th}` Hamming code over `F=GF(q)` is an
    `[n,k,d]` code with length `n=(q^r-1)/(q-1)`,
    dimension `k=(q^r-1)/(q-1) - r` and minimum distance
    `d=3`. The parity check matrix of a Hamming code has rows
    consisting of all nonzero vectors of length r in its columns,
    modulo a scalar factor so no parallel columns arise. A Hamming code
    is a single error-correcting code.

    INPUT:


    -  ``r`` - an integer 2

    -  ``F`` - a finite field.


    OUTPUT: Returns the r-th q-ary Hamming code.

    EXAMPLES::

        sage: codes.HammingCode(3,GF(2))
        Linear code of length 7, dimension 4 over Finite Field of size 2
        sage: C = codes.HammingCode(3,GF(3)); C
        Linear code of length 13, dimension 10 over Finite Field of size 3
        sage: C.minimum_distance()
        3
        sage: C.minimum_distance(algorithm='gap') # long time, check d=3
        3
        sage: C = codes.HammingCode(3,GF(4,'a')); C
        Linear code of length 21, dimension 18 over Finite Field in a of size 2^2

    While the ``codes`` object now gathers all code constructors,
    ``HammingCode`` is still available in the global namespace::

        sage: HammingCode(3,GF(2))
        doctest:1: DeprecationWarning: This method soon will not be available in that way anymore. To use it, you can now call it by typing codes.HammingCode
        See http://trac.sagemath.org/15445 for details.
        Linear code of length 7, dimension 4 over Finite Field of size 2

    """
    q = F.order()
    n =  (q**r-1)/(q-1)
    k = n-r
    MS = MatrixSpace(F,n,r)
    X = ProjectiveSpace(r-1,F)
    PFn = [list(p) for p in X.point_set(F).points(F)]
    H = MS(PFn).transpose()
    Cd = LinearCode(H)
    # Hamming code always has distance 3, so we provide the distance.
    return LinearCode(Cd.dual_code().gen_mat(), d=3)


def LinearCodeFromCheckMatrix(H):
    r"""
    A linear [n,k]-code C is uniquely determined by its generator
    matrix G and check matrix H. We have the following short exact
    sequence

    .. math::

        0 \rightarrow
        {\mathbf{F}}^k \stackrel{G}{\rightarrow}
        {\mathbf{F}}^n \stackrel{H}{\rightarrow}
        {\mathbf{F}}^{n-k} \rightarrow
        0.

    ("Short exact" means (a) the arrow `G` is injective, i.e.,
    `G` is a full-rank `k\times n` matrix, (b) the
    arrow `H` is surjective, and (c)
    `{\rm image}(G)={\rm kernel}(H)`.)

    EXAMPLES::

        sage: C = codes.HammingCode(3,GF(2))
        sage: H = C.check_mat(); H
        [1 0 1 0 1 0 1]
        [0 1 1 0 0 1 1]
        [0 0 0 1 1 1 1]
        sage: codes.LinearCodeFromCheckMatrix(H) == C
        True
        sage: C = codes.HammingCode(2,GF(3))
        sage: H = C.check_mat(); H
        [1 0 1 1]
        [0 1 1 2]
        sage: codes.LinearCodeFromCheckMatrix(H) == C
        True
        sage: C = codes.RandomLinearCode(10,5,GF(4,"a"))
        sage: H = C.check_mat()
        sage: codes.LinearCodeFromCheckMatrix(H) == C
        True
    """
    Cd = LinearCode(H)
    return Cd.dual_code()

def QuadraticResidueCode(n,F):
    r"""
    A quadratic residue code (or QR code) is a cyclic code whose
    generator polynomial is the product of the polynomials
    `x-\alpha^i` (`\alpha` is a primitive
    `n^{th}` root of unity; `i` ranges over the set of
    quadratic residues modulo `n`).

    See QuadraticResidueCodeEvenPair and QuadraticResidueCodeOddPair
    for a more general construction.

    INPUT:


    -  ``n`` - an odd prime

    -  ``F`` - a finite prime field F whose order must be a
       quadratic residue modulo n.


    OUTPUT: Returns a quadratic residue code.

    EXAMPLES::

        sage: C = codes.QuadraticResidueCode(7,GF(2))
        sage: C
        Linear code of length 7, dimension 4 over Finite Field of size 2
        sage: C = codes.QuadraticResidueCode(17,GF(2))
        sage: C
        Linear code of length 17, dimension 9 over Finite Field of size 2
        sage: C1 = codes.QuadraticResidueCodeOddPair(7,GF(2))[0]
        sage: C2 = codes.QuadraticResidueCode(7,GF(2))
        sage: C1 == C2
        True
        sage: C1 = codes.QuadraticResidueCodeOddPair(17,GF(2))[0]
        sage: C2 = codes.QuadraticResidueCode(17,GF(2))
        sage: C1 == C2
        True

    AUTHORS:

    - David Joyner (11-2005)
    """
    return QuadraticResidueCodeOddPair(n,F)[0]

def QuadraticResidueCodeEvenPair(n,F):
    """
    Quadratic residue codes of a given odd prime length and base ring
    either don't exist at all or occur as 4-tuples - a pair of
    "odd-like" codes and a pair of "even-like" codes. If n 2 is prime
    then (Theorem 6.6.2 in [HP]_) a QR code exists over GF(q) iff q is a
    quadratic residue mod n.

    They are constructed as "even-like" duadic codes associated the
    splitting (Q,N) mod n, where Q is the set of non-zero quadratic
    residues and N is the non-residues.

    EXAMPLES::

        sage: codes.QuadraticResidueCodeEvenPair(17,GF(13))
        (Linear code of length 17, dimension 8 over Finite Field of size 13,
         Linear code of length 17, dimension 8 over Finite Field of size 13)
        sage: codes.QuadraticResidueCodeEvenPair(17,GF(2))
        (Linear code of length 17, dimension 8 over Finite Field of size 2,
         Linear code of length 17, dimension 8 over Finite Field of size 2)
        sage: codes.QuadraticResidueCodeEvenPair(13,GF(9,"z"))
        (Linear code of length 13, dimension 6 over Finite Field in z of size 3^2,
         Linear code of length 13, dimension 6 over Finite Field in z of size 3^2)
        sage: C1 = codes.QuadraticResidueCodeEvenPair(7,GF(2))[0]
        sage: C1.is_self_orthogonal()
        True
        sage: C2 = codes.QuadraticResidueCodeEvenPair(7,GF(2))[1]
        sage: C2.is_self_orthogonal()
        True
        sage: C3 = codes.QuadraticResidueCodeOddPair(17,GF(2))[0]
        sage: C4 = codes.QuadraticResidueCodeEvenPair(17,GF(2))[1]
        sage: C3 == C4.dual_code()
        True

    This is consistent with Theorem 6.6.9 and Exercise 365 in [HP]_.
    """
    q = F.order()
    Q = quadratic_residues(n); Q.remove(0)  # non-zero quad residues
    N = range(1,n); tmp = [N.remove(x) for x in Q]  # non-zero quad non-residues
    if (n.is_prime() and n>2 and not(q in Q)):
        raise ValueError, "No quadratic residue code exists for these parameters."
    if not(is_a_splitting(Q,N,n)):
        raise TypeError, "No quadratic residue code exists for these parameters."
    return DuadicCodeEvenPair(F,Q,N)


def QuadraticResidueCodeOddPair(n,F):
    """
    Quadratic residue codes of a given odd prime length and base ring
    either don't exist at all or occur as 4-tuples - a pair of
    "odd-like" codes and a pair of "even-like" codes. If n 2 is prime
    then (Theorem 6.6.2 in [HP]_) a QR code exists over GF(q) iff q is a
    quadratic residue mod n.

    They are constructed as "odd-like" duadic codes associated the
    splitting (Q,N) mod n, where Q is the set of non-zero quadratic
    residues and N is the non-residues.

    EXAMPLES::

        sage: codes.QuadraticResidueCodeOddPair(17,GF(13))
        (Linear code of length 17, dimension 9 over Finite Field of size 13,
         Linear code of length 17, dimension 9 over Finite Field of size 13)
        sage: codes.QuadraticResidueCodeOddPair(17,GF(2))
        (Linear code of length 17, dimension 9 over Finite Field of size 2,
         Linear code of length 17, dimension 9 over Finite Field of size 2)
        sage: codes.QuadraticResidueCodeOddPair(13,GF(9,"z"))
        (Linear code of length 13, dimension 7 over Finite Field in z of size 3^2,
         Linear code of length 13, dimension 7 over Finite Field in z of size 3^2)
        sage: C1 = codes.QuadraticResidueCodeOddPair(17,GF(2))[1]
        sage: C1x = C1.extended_code()
        sage: C2 = codes.QuadraticResidueCodeOddPair(17,GF(2))[0]
        sage: C2x = C2.extended_code()
        sage: C2x.spectrum(); C1x.spectrum()
        [1, 0, 0, 0, 0, 0, 102, 0, 153, 0, 153, 0, 102, 0, 0, 0, 0, 0, 1]
        [1, 0, 0, 0, 0, 0, 102, 0, 153, 0, 153, 0, 102, 0, 0, 0, 0, 0, 1]
        sage: C2x == C1x.dual_code()
        True
        sage: C3 = codes.QuadraticResidueCodeOddPair(7,GF(2))[0]
        sage: C3x = C3.extended_code()
        sage: C3x.spectrum()
        [1, 0, 0, 0, 14, 0, 0, 0, 1]
        sage: C3x.is_self_dual()
        True

    This is consistent with Theorem 6.6.14 in [HP]_.
    """
    from sage.coding.code_constructions import is_a_splitting
    q = F.order()
    Q = quadratic_residues(n); Q.remove(0)  # non-zero quad residues
    N = range(1,n); tmp = [N.remove(x) for x in Q]  # non-zero quad non-residues
    if (n.is_prime() and n>2 and not(q in Q)):
        raise ValueError, "No quadratic residue code exists for these parameters."
    if not(is_a_splitting(Q,N,n)):
        raise TypeError, "No quadratic residue code exists for these parameters."
    return DuadicCodeOddPair(F,Q,N)


def RandomLinearCode(n,k,F):
    r"""
    The method used is to first construct a `k \times n`
    matrix using Sage's random_element method for the MatrixSpace
    class. The construction is probabilistic but should only fail
    extremely rarely.

    INPUT: Integers n,k, with `n>k`, and a finite field F

    OUTPUT: Returns a "random" linear code with length n, dimension k
    over field F.

    EXAMPLES::

        sage: C = codes.RandomLinearCode(30,15,GF(2))
        sage: C
        Linear code of length 30, dimension 15 over Finite Field of size 2
        sage: C = codes.RandomLinearCode(10,5,GF(4,'a'))
        sage: C
        Linear code of length 10, dimension 5 over Finite Field in a of size 2^2

    AUTHORS:

    - David Joyner (2007-05)
    """
    MS = MatrixSpace(F,k,n)
    for i in range(50):
        G = MS.random_element()
        if G.rank() == k:
            V = span(G.rows(), F)
            return LinearCodeFromVectorSpace(V)  # may not be in standard form
    MS1 = MatrixSpace(F,k,k)
    MS2 = MatrixSpace(F,k,n-k)
    Ik = MS1.identity_matrix()
    A = MS2.random_element()
    G = Ik.augment(A)
    return LinearCode(G)                          # in standard form


def ReedSolomonCode(n,k,F,pts = None):
    r"""
    Given a finite field `F` of order `q`, let
    `n` and `k` be chosen such that
    `1 \leq k \leq n \leq q`. Pick `n` distinct
    elements of `F`, denoted
    `\{ x_1, x_2, ... , x_n \}`. Then, the codewords are
    obtained by evaluating every polynomial in `F[x]` of degree
    less than `k` at each `x_i`:

    .. math::

       C = \left\{ \left( f(x_1), f(x_2), ..., f(x_n) \right), f \in F[x],
       {\rm deg}(f)<k \right\}.


    `C` is a `[n, k, n-k+1]` code. (In particular, `C` is MDS.)

    INPUT: n : the length k : the dimension F : the base ring pts :
    (optional) list of n points in F (if None then Sage picks n of them
    in the order given to the elements of F)

    EXAMPLES::

        sage: C = codes.ReedSolomonCode(6,4,GF(7)); C
        Linear code of length 6, dimension 4 over Finite Field of size 7
        sage: C.minimum_distance()
        3
        sage: C = codes.ReedSolomonCode(6,4,GF(8,"a")); C
        Linear code of length 6, dimension 4 over Finite Field in a of size 2^3
        sage: C.minimum_distance()
        3
        sage: C.minimum_distance(algorithm='gap') # long time, check d=n-k+1
        3
        sage: F.<a> = GF(3^2,"a")
        sage: pts = [0,1,a,a^2,2*a,2*a+1]
        sage: len(Set(pts)) == 6 # to make sure there are no duplicates
        True
        sage: C = codes.ReedSolomonCode(6,4,F,pts); C
        Linear code of length 6, dimension 4 over Finite Field in a of size 3^2
        sage: C.minimum_distance()
        3

    While the ``codes`` object now gathers all code constructors,
    ``ReedSolomonCode`` is still available in the global namespace::

        sage: ReedSolomonCode(6,4,GF(7))
        doctest:1: DeprecationWarning: This method soon will not be available in that way anymore. To use it, you can now call it by typing codes.ReedSolomonCode
        See http://trac.sagemath.org/15445 for details.
        Linear code of length 6, dimension 4 over Finite Field of size 7

    REFERENCES:

    - [W] http://en.wikipedia.org/wiki/Reed-Solomon
    """
    q = F.order()
    power = lambda x,n,F: (x==0 and n==0) and F(1) or F(x**n) # since 0^0 is undefined
    if n>q or k>n or k>q:
        raise ValueError, "RS codes does not exist with the given input."
    if not(pts == None) and not(len(pts)==n):
        raise ValueError, "You must provide exactly %s distinct points of %s"%(n,F)
    if (pts == None):
        pts = []
        i = 0
        for x in F:
            if i<n:
                pts.append(x)
                i = i+1
    MS = MatrixSpace(F, k, n)
    rowsG = []
    rowsG = [[power(x,j,F) for x in pts] for j in range(k)]
    G = MS(rowsG)
    return LinearCode(G, d=n-k+1)


def TernaryGolayCode():
    r"""
    TernaryGolayCode returns a ternary Golay code. This is a perfect
    [11,6,5] code. It is also equivalent to a cyclic code, with
    generator polynomial `g(x)=2+x^2+2x^3+x^4+x^5`.

    EXAMPLES::

        sage: C = codes.TernaryGolayCode()
        sage: C
        Linear code of length 11, dimension 6 over Finite Field of size 3
        sage: C.minimum_distance()
        5
        sage: C.minimum_distance(algorithm='gap') # long time, check d=5
        5

    AUTHORS:

    - David Joyner (2007-5)
    """
    F = GF(3)
    B = [[2, 0, 1, 2, 1, 1, 0, 0, 0, 0, 0],\
         [0, 2, 0, 1, 2, 1, 1, 0, 0, 0, 0],\
         [0, 0, 2, 0, 1, 2, 1, 1, 0, 0, 0],\
         [0, 0, 0, 2, 0, 1, 2, 1, 1, 0, 0],\
         [0, 0, 0, 0, 2, 0, 1, 2, 1, 1, 0],\
         [0, 0, 0, 0, 0, 2, 0, 1, 2, 1, 1]]
    V = span(B, F)
    return LinearCodeFromVectorSpace(V, d=5)

def ToricCode(P,F):
    r"""
    Let `P` denote a list of lattice points in
    `\ZZ^d` and let `T` denote the set of all
    points in `(F^x)^d` (ordered in some fixed way). Put
    `n=|T|` and let `k` denote the dimension of the
    vector space of functions `V = \mathrm{Span}\{x^e \ |\ e \in P\}`.
    The associated toric code `C` is the evaluation code which
    is the image of the evaluation map

    .. math::

        \mathrm{eval_T} : V \rightarrow F^n,


    where `x^e` is the multi-index notation
    (`x=(x_1,...,x_d)`, `e=(e_1,...,e_d)`, and
    `x^e = x_1^{e_1}...x_d^{e_d}`), where
    `eval_T (f(x)) = (f(t_1),...,f(t_n))`, and where
    `T=\{t_1,...,t_n\}`. This function returns the toric
    codes discussed in [J]_.

    INPUT:


    -  ``P`` - all the integer lattice points in a polytope
       defining the toric variety.

    -  ``F`` - a finite field.


    OUTPUT: Returns toric code with length n = , dimension k over field
    F.

    EXAMPLES::

         sage: C = codes.ToricCode([[0,0],[1,0],[2,0],[0,1],[1,1]],GF(7))
         sage: C
         Linear code of length 36, dimension 5 over Finite Field of size 7
         sage: C.minimum_distance()
         24
         sage: C = codes.ToricCode([[-2,-2],[-1,-2],[-1,-1],[-1,0],[0,-1],[0,0],[0,1],[1,-1],[1,0]],GF(5))
         sage: C
         Linear code of length 16, dimension 9 over Finite Field of size 5
         sage: C.minimum_distance()
         6
         sage: C = codes.ToricCode([ [0,0],[1,1],[1,2],[1,3],[1,4],[2,1],[2,2],[2,3],[3,1],[3,2],[4,1]],GF(8,"a"))
         sage: C
         Linear code of length 49, dimension 11 over Finite Field in a of size 2^3

    This is in fact a [49,11,28] code over GF(8). If you type next
    ``C.minimum_distance()`` and wait overnight (!), you
    should get 28.

    AUTHOR:

    - David Joyner (07-2006)

    REFERENCES:

    .. [J] D. Joyner, Toric codes over finite fields, Applicable
       Algebra in Engineering, Communication and Computing, 15, (2004), p. 63-79.
    """
    from sage.combinat.all import Tuples
    mset = [x for x in F if x!=0]
    d = len(P[0])
    pts = Tuples(mset,d).list()
    n = len(pts) # (q-1)^d
    k = len(P)
    e = P[0]
    B = []
    for e in P:
       tmpvar = [prod([t[i]**e[i] for i in range(d)]) for t in pts]
       B.append(tmpvar)
    # now B0 *should* be a full rank matrix
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(B))


def TrivialCode(F,n):
    MS = MatrixSpace(F,1,n)
    return LinearCode(MS(0))


def WalshCode(m):
    r"""
    Returns the binary Walsh code of length `2^m`. The matrix
    of codewords correspond to a Hadamard matrix. This is a (constant
    rate) binary linear `[2^m,m,2^{m-1}]` code.

    EXAMPLES::

        sage: C = codes.WalshCode(4); C
        Linear code of length 16, dimension 4 over Finite Field of size 2
        sage: C = codes.WalshCode(3); C
        Linear code of length 8, dimension 3 over Finite Field of size 2
        sage: C.spectrum()
        [1, 0, 0, 0, 7, 0, 0, 0, 0]
        sage: C.minimum_distance()
        4
        sage: C.minimum_distance(algorithm='gap') # check d=2^(m-1)
        4

    REFERENCES:

    - http://en.wikipedia.org/wiki/Hadamard_matrix

    - http://en.wikipedia.org/wiki/Walsh_code
    """
    return LinearCode(walsh_matrix(m), d=2**(m-1))
