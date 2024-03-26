# sage_setup: distribution = sagemath-repl
# sage.doctest: needs sage.libs.pari sage.modules
r"""
This file is meant to catch errors in the PARI/GP package which are not
caught by any other tests.

Check that :issue:`9876` has been fixed, this test comes from PARI's
self-test :pari:`rnfkummer` but was modified such that the answer is
canonical::

    sage: pari('K = bnfinit(y^4-52*y^2+26,1); pol = rnfkummer(bnrinit(K,3,1),Mat(5)); L = rnfinit(K, [pol, 10^6]); polredabs(polredbest(L.polabs))')  # long time
    x^20 - 112*x^18 + 5108*x^16 - 123460*x^14 + 1724337*x^12 - 14266996*x^10 + 69192270*x^8 - 188583712*x^6 + 260329852*x^4 - 141461008*x^2 + 19860776

Check that :issue:`10195` (PARI bug 1153) has been fixed::

    sage: print(gp.eval("mathnf([0,0,0,0,0,0,0,0,0,13;0,0,0,0,0,0,0,0,23,6;0,0,0,0,0,0,0,23,-4,-7;0,0,0,0,0,0,17,-3,5,-5;0,0,0,0,0,56,16,-16,-15,-17;0,0,0,0,57,24,-16,-25,2,-21;0,0,0,114,9,56,51,-52,25,-55;0,0,113,-31,-11,24,0,28,34,-16;0,50,3,2,16,-6,-2,7,-19,-21;118,43,51,23,37,-52,18,38,51,28],0)"))
    [787850171872400 32189386376004 356588299060422 742392731867995 282253457851430 665185047494955 664535243562463 744564809133574 113975061998590 527459013372200]
    [0 12 6 11 5 3 7 6 6 0]
    [0 0 3 1 2 1 1 0 0 0]
    [0 0 0 1 0 0 0 0 0 0]
    [0 0 0 0 1 0 0 0 0 0]
    [0 0 0 0 0 1 0 0 0 0]
    [0 0 0 0 0 0 1 0 0 0]
    [0 0 0 0 0 0 0 1 0 0]
    [0 0 0 0 0 0 0 0 1 0]
    [0 0 0 0 0 0 0 0 0 1]

Check that :issue:`11604` (PARI bug 1154) has been fixed::

    sage: A = Matrix(ZZ,4,4,[32982266684193100, 1368614777139719, 224591013270052693, 276460184982223238,1368614777139719,  56791380087354, 9319512049770279, 11471848267545007,224591013270052693, 9319512049770279,1529340971891522140, 1882541434053596358,276460184982223238, 11471848267545007, 1882541434053596358, 2317313350044091414])
    sage: pari(A).qfminim(2,0)
    [0, 0, [;]]

Check :issue:`13314`, the following should not give a
Segmentation Fault::

    sage: x = polygen(ComplexField(128))
    sage: p = x^84 + (16*x^4 - 1)^20 * (2^48*x^4 - 2049^4)
    sage: len(pari(p).polroots(precision=128))
    84

Check that the optional PARI databases work::

    sage: # optional - pari_elldata
    sage: gp.ellinit('"299998a1"')
    [1, 0, 1, 110, -3660, ...]
    sage: E = EllipticCurve("1728ba1")                                                  # needs sage.schemes
    sage: gp(E).ellidentify()                                                           # needs sage.schemes
    [["1728ba1", [0, 0, 0, -6, 6], [[1, 1]]], [1, 0, 0, 0]]

    sage: pari("ellmodulareqn(211)")                            # optional - pari_seadata
    [x^212 + (-y^7 + 5207*y^6 - 10241606*y^5 + 9430560101*y^4 - 4074860204015*y^3
      + 718868274900397*y^2 - 34897101275826114*y + 104096378056356968)*x^211...

The following requires the modular polynomials up to degree 223, while
only those up to degree 199 come standard in Sage::

    sage: p = next_prime(2^328)
    sage: E = EllipticCurve(GF(p), [6,1])                                               # needs sage.schemes
    sage: E.cardinality()  # long time (108s on sage.math, 2013), optional - pari_seadata, needs sage.schemes
    546812681195752981093125556779405341338292357723293496548601032930284335897180749997402596957976244

Create a number field with Galois group `A4`. Group `A4` corresponds to
transitive group `(12,3)` in GAP::

    sage: # optional - pari_galpol
    sage: R.<x> = PolynomialRing(ZZ)
    sage: pol = pari("galoisgetpol(12,3)[1]")
    sage: K.<a> = NumberField(R(pol))
    sage: factor(K.discriminant())
    163^8
    sage: [F.degree() for F,a,b in K.subfields()]
    [1, 3, 4, 4, 4, 4, 6, 6, 6, 12]
    sage: sorted([12/H.cardinality() for H in AlternatingGroup(4).subgroups()])
    [1, 3, 4, 4, 4, 4, 6, 6, 6, 12]
"""
