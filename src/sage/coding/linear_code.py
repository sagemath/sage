# -*- coding: utf-8 -*-
r"""
Linear Codes

VERSION: 1.2

Let `F` be a finite field.  Here, we will denote the finite field with `q`
elements by `\GF{q}`.  A subspace of `F^n` (with the standard basis) is
called a linear code of length `n`.  If its dimension is denoted `k` then we
typically store a basis of `C` as a `k \times n` matrix, with rows the basis
vectors.  It is called the generator matrix of `C`. The rows of the parity
check matrix of `C` are a basis for the code,

.. math::

      C^* = \{ v \in GF(q)^n\ |\ v\cdot c = 0,\ for \ all\ c \in C \},

called the dual space of `C`.

If `F=\GF{2}` then `C` is called a binary code.  If `F = \GF{q}` then `C` is
called a `q`-ary code.  The elements of a code `C` are called codewords.

The symmetric group `S_n` acts on `F^n` by permuting coordinates.  If an
element `p \in S_n` sends a code `C` of length `n` to itself (in other words,
every codeword of `C` is sent to some other codeword of `C`) then `p` is
called a permutation automorphism of `C`.  The (permutation) automorphism
group is denoted `Aut(C)`.

This file contains

#. LinearCode class definition; LinearCodeFromVectorspace conversion function,

#. The spectrum (weight distribution), covering_radius, minimum distance
   programs (calling Steve Linton's or CJ Tjhal's C programs),
   characteristic_function, and several implementations of the Duursma zeta
   function (sd_zeta_polynomial, zeta_polynomial, zeta_function,
   chinen_polynomial, for example),

#. interface with best_known_linear_code_www (interface with codetables.de
   since A. Brouwer's online tables have been disabled),
   bounds_minimum_distance which call tables in GUAVA (updated May 2006)
   created by Cen Tjhai instead of the online internet tables,

#. gen_mat, gen_mat_systematic, information_set, list, check_mat, decode,
   dual_code, extended_code, shortened, punctured, genus, binomial_moment,
   and divisor methods for LinearCode,

#. Boolean-valued functions such as "==", is_self_dual, is_self_orthogonal,
   is_subcode, is_permutation_automorphism, is_permutation_equivalent (which
   interfaces with Robert Miller's partition refinement code),

#. permutation methods: is_permutation_automorphism,
   permutation_automorphism_group, permuted_code, standard_form,
   module_composition_factors,

#. design-theoretic methods: assmus_mattson_designs (implementing
   Assmus-Mattson Theorem),

#. code constructions, such as HammingCode and ToricCode, are in a separate
   ``code_constructions.py`` module; in the separate ``guava.py`` module, you
   will find constructions, such as RandomLinearCodeGuava and
   BinaryReedMullerCode, wrapped from the corresponding GUAVA codes.

EXAMPLES::

    sage: MS = MatrixSpace(GF(2),4,7)
    sage: G = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
    sage: C = LinearCode(G)
    sage: C.basis()
    [(1, 1, 1, 0, 0, 0, 0),
     (1, 0, 0, 1, 1, 0, 0),
     (0, 1, 0, 1, 0, 1, 0),
     (1, 1, 0, 1, 0, 0, 1)]
    sage: c = C.basis()[1]
    sage: c in C
    True
    sage: c.nonzero_positions()
    [0, 3, 4]
    sage: c.support()
    [0, 3, 4]
    sage: c.parent()
    Vector space of dimension 7 over Finite Field of size 2

To be added:

#. More wrappers

#. GRS codes and special decoders.

#. `P^1` Goppa codes and group actions on `P^1` RR space codes.

REFERENCES:

- [HP] W. C. Huffman and V. Pless, Fundamentals of error-correcting codes,
  Cambridge Univ. Press, 2003.

- [Gu] GUAVA manual, http://www.gap-system.org/Packages/guava.html

AUTHORS:

- David Joyner (2005-11-22, 2006-12-03): initial version

- William Stein (2006-01-23): Inclusion in Sage

- David Joyner (2006-01-30, 2006-04): small fixes

- David Joyner (2006-07): added documentation, group-theoretical methods,
  ToricCode

- David Joyner (2006-08): hopeful latex fixes to documentation, added list and
  __iter__ methods to LinearCode and examples, added hamming_weight function,
  fixed random method to return a vector, TrivialCode, fixed subtle bug in
  dual_code, added galois_closure method, fixed mysterious bug in
  permutation_automorphism_group (GAP was over-using "G" somehow?)

- David Joyner (2006-08): hopeful latex fixes to documentation, added
  CyclicCode, best_known_linear_code, bounds_minimum_distance,
  assmus_mattson_designs (implementing Assmus-Mattson Theorem).

- David Joyner (2006-09): modified decode syntax, fixed bug in
  is_galois_closed, added LinearCode_from_vectorspace, extended_code,
  zeta_function

- Nick Alexander (2006-12-10): factor GUAVA code to guava.py

- David Joyner (2007-05): added methods punctured, shortened, divisor,
  characteristic_polynomial, binomial_moment, support for
  LinearCode. Completely rewritten zeta_function (old version is now
  zeta_function2) and a new function, LinearCodeFromVectorSpace.

- David Joyner (2007-11): added zeta_polynomial, weight_enumerator,
  chinen_polynomial; improved best_known_code; made some pythonic revisions;
  added is_equivalent (for binary codes)

- David Joyner (2008-01): fixed bug in decode reported by Harald Schilly,
  (with Mike Hansen) added some doctests.

- David Joyner (2008-02): translated standard_form, dual_code to Python.

- David Joyner (2008-03): translated punctured, shortened, extended_code,
  random (and renamed random to random_element), deleted zeta_function2,
  zeta_function3, added wrapper automorphism_group_binary_code to Robert
  Miller's code), added direct_sum_code, is_subcode, is_self_dual,
  is_self_orthogonal, redundancy_matrix, did some alphabetical reorganizing
  to make the file more readable. Fixed a bug in permutation_automorphism_group
  which caused it to crash.

- David Joyner (2008-03): fixed bugs in spectrum and zeta_polynomial, which
  misbehaved over non-prime base rings.

- David Joyner (2008-10): use CJ Tjhal's MinimumWeight if char = 2 or 3 for
  min_dist; add is_permutation_equivalent and improve
  permutation_automorphism_group using an interface with Robert Miller's code;
  added interface with Leon's code for the spectrum method.

- David Joyner (2009-02): added native decoding methods (see module_decoder.py)

- David Joyner (2009-05): removed dependence on Guava, allowing it to be an
  option. Fixed errors in some docstrings.

- Kwankyu Lee (2010-01): added methods gen_mat_systematic, information_set, and
  magma interface for linear codes.

- Niles Johnson (2010-08): :trac:`#3893`: ``random_element()`` should pass on ``*args`` and ``**kwds``.

- Thomas Feulner (2012-11): :trac:`13723`: deprecation of ``hamming_weight()``

TESTS::

    sage: MS = MatrixSpace(GF(2),4,7)
    sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
    sage: C  = LinearCode(G)
    sage: C == loads(dumps(C))
    True
"""

#******************************************************************************
#       Copyright (C) 2005 David Joyner <wdjoyner@gmail.com>
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#******************************************************************************

import urllib
import sage.modules.free_module as fm
import sage.modules.module as module
from sage.interfaces.all import gap
from sage.rings.finite_rings.constructor import FiniteField as GF
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector
from sage.rings.arith import GCD, rising_factorial, binomial
from sage.groups.all import SymmetricGroup
from sage.misc.misc import prod
from sage.misc.functional import log, is_even
from sage.rings.rational_field import QQ
from sage.structure.parent_gens import ParentWithGens
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.fraction_field import FractionField
from sage.rings.integer_ring import IntegerRing
from sage.rings.integer import Integer
from sage.combinat.set_partition import SetPartitions
from sage.misc.randstate import current_randstate
from sage.misc.decorators import rename_keyword

ZZ = IntegerRing()
VectorSpace = fm.VectorSpace

####################### coding theory functions ###############################

def hamming_weight(v):
    r"""
    Returns the Hamming weight of the vector ``v``, which is the number of
    non-zero entries.

    INPUT:

    - ``v`` - Vector

    OUTPUT:

    - Integer, the Hamming weight of ``v``

    EXAMPLES::

        sage: hamming_weight(vector(GF(2),[0,0,1]))
        doctest:1: DeprecationWarning: The global function hamming_weight(v) is
        deprecated, instead use v.hamming_weight().
        See http://trac.sagemath.org/13723 for details.
        1
        sage: hamming_weight(vector(GF(2),[0,0,0]))
        0
    """
    from sage.misc.superseded import deprecation
    deprecation(13723, "The global function hamming_weight(v) is deprecated, " +
                "instead use v.hamming_weight().")

    return len(v.nonzero_positions())

def code2leon(C):
    r"""
    Writes a file in Sage's temp directory representing the code C, returning
    the absolute path to the file. This is the Sage translation of the
    GuavaToLeon command in Guava's codefun.gi file.

    INPUT:

    - ``C`` - a linear code (over GF(p), p < 11)

    OUTPUT:

    - Absolute path to the file written

    EXAMPLES::

        sage: C = codes.HammingCode(3,GF(2)); C
        Linear code of length 7, dimension 4 over Finite Field of size 2
        sage: file_loc = sage.coding.linear_code.code2leon(C)
        sage: f = open(file_loc); print f.read()
        LIBRARY code;
        code=seq(2,4,7,seq(
        1,0,0,0,0,1,1,
        0,1,0,0,1,0,1,
        0,0,1,0,1,1,0,
        0,0,0,1,1,1,1
        ));
        FINISH;
        sage: f.close()

    """
    from sage.misc.misc import tmp_filename
    F = C.base_ring()
    p = F.order()  # must be prime and <11
    s = "LIBRARY code;\n"+"code=seq(%s,%s,%s,seq(\n"%(p,C.dimension(),C.length())
    Gr = [str(r)[1:-1].replace(" ","") for r in C.gen_mat().rows()]
    s += ",\n".join(Gr) + "\n));\nFINISH;"
    file_loc = tmp_filename()
    f = open(file_loc,"w")
    f.write(s)
    f.close()
    return file_loc

def wtdist_gap(Gmat, n, F):
    r"""
    INPUT:

    -  ``Gmat`` - String representing a GAP generator matrix G of a linear code

    -  ``n`` - Integer greater than 1, representing the number of columns of G
       (i.e., the length of the linear code)

    -  ``F`` - Finite field (in Sage), base field the code

    OUTPUT:

    -  Spectrum of the associated code

    EXAMPLES::

        sage: Gstr = 'Z(2)*[[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]]'
        sage: F = GF(2)
        sage: sage.coding.linear_code.wtdist_gap(Gstr, 7, F)
        [1, 0, 0, 7, 7, 0, 0, 1]

    Here ``Gstr`` is a generator matrix of the Hamming [7,4,3] binary code.

    ALGORITHM:

    Uses C programs written by Steve Linton in the kernel of GAP, so is fairly
    fast.

    AUTHORS:

    - David Joyner (2005-11)
    """
    G = gap(Gmat)
    q = F.order()
    k = gap(F)
    #C = G.GeneratorMatCode(k)
    #n = int(C.WordLength())
    z = 'Z(%s)*%s'%(q, [0]*n)     # GAP zero vector as a string
    _ = gap.eval("w:=DistancesDistributionMatFFEVecFFE("+Gmat+", GF("+str(q)+"),"+z+")")
    # for some reason, this commented code doesn't work:
    #dist0 = gap("DistancesDistributionMatFFEVecFFE("+Gmat+", GF("+str(q)+"),"+z+")")
    #v0 = dist0._matrix_(F)
    #print dist0,v0
    #d = G.DistancesDistributionMatFFEVecFFE(k, z)
    v = [eval(gap.eval("w["+str(i)+"]")) for i in range(1,n+2)] # because GAP returns vectors in compressed form
    return v

@rename_keyword(deprecation=11033, method="algorithm")
def min_wt_vec_gap(Gmat, n, k, F, algorithm=None):
    r"""
    Returns a minimum weight vector of the code generated by ``Gmat``.

    Uses C programs written by Steve Linton in the kernel of GAP, so is fairly
    fast. The option ``algorithm="guava"`` requires Guava. The default algorithm
    requires GAP but not Guava.

    INPUT:

    -  ``Gmat`` - String representing a GAP generator matrix G of a linear code
    -  n - Length of the code generated by G
    -  k - Dimension of the code generated by G
    -  F - Base field

    OUTPUT:

    -  Minimum weight vector of the code generated by ``Gmat``

    REMARKS:

    - The code in the default case allows one (for free) to also compute the
      message vector `m` such that `m\*G = v`, and the (minimum) distance, as
      a triple.  however, this output is not implemented.
    - The binary case can presumably be done much faster using Robert Miller's
      code (see the docstring for the spectrum method). This is also not (yet)
      implemented.

    EXAMPLES::

        sage: Gstr = "Z(2)*[[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]]"
        sage: sage.coding.linear_code.min_wt_vec_gap(Gstr,7,4,GF(2))
        (0, 1, 0, 1, 0, 1, 0)

        This output is different but still a minimum weight vector:

        sage: sage.coding.linear_code.min_wt_vec_gap(Gstr,7,4,GF(2),algorithm="guava")    # optional - gap_packages (Guava package)
        (0, 0, 1, 0, 1, 1, 0)

    Here ``Gstr`` is a generator matrix of the Hamming [7,4,3] binary code.

    AUTHORS:

    - David Joyner (11-2005)
    """
    current_randstate().set_seed_gap()
    if algorithm=="guava":
        gap.LoadPackage('"guava"')
        from sage.interfaces.gap import gfq_gap_to_sage
        gap.eval("G:="+Gmat)
        C = gap(Gmat).GeneratorMatCode(F)
        #n = int(C.length())
        cg = C.MinimumDistanceCodeword()
        c = [gfq_gap_to_sage(cg[j],F) for j in range(1,n+1)]
        V = VectorSpace(F,n)
        return V(c)
    qstr = str(F.order())
    zerovec = [0 for i in range(n)]
    zerovecstr = "Z("+qstr+")*"+str(zerovec)
    all = []
    gap.eval('Gmat:='+Gmat)
    for i in range(1,k+1):
        gap.eval("P:=AClosestVectorCombinationsMatFFEVecFFECoords(Gmat, GF("+qstr+"),"+zerovecstr+","+str(i)+",0); d:=WeightVecFFE(P[1])")
        v = gap("[P[1]]")
        m = gap("[P[2]]")
        dist = gap("d")
        #print v,m,dist
        #print [gap.eval("v["+str(i+1)+"]") for i in range(n)]
        all.append([v._matrix_(F), m._matrix_(F), int(dist)])
    ans = all[0]
    for x in all:
        if x[2]<ans[2] and x[2]>0:
            ans = x
    #print ans[0], ans[0].parent()
    return vector(F,[x for x in ans[0].rows()[0]]) # ugly 1xn matrix->vector coercion!

def best_known_linear_code(n, k, F):
    r"""
    Returns the best known (as of 11 May 2006) linear code of length ``n``,
    dimension ``k`` over field ``F``.  The function uses the tables described
    in ``bounds_minimum_distance`` to construct this code.

    This does not require an internet connection.

    EXAMPLES::

        sage: best_known_linear_code(10,5,GF(2))    # long time; optional - gap_packages (Guava package)
        Linear code of length 10, dimension 5 over Finite Field of size 2
        sage: gap.eval("C:=BestKnownLinearCode(10,5,GF(2))")     # long time; optional - gap_packages (Guava package)
        'a linear [10,5,4]2..4 shortened code'

    This means that best possible binary linear code of length 10 and
    dimension 5 is a code with minimum distance 4 and covering radius
    somewhere between 2 and 4. Use ``minimum_distance_why(10,5,GF(2))`` or
    ``print bounds_minimum_distance(10,5,GF(2))`` for further details.
    """
    q = F.order()
    C = gap("BestKnownLinearCode(%s,%s,GF(%s))"%(n,k,q))
    G = C.GeneratorMat()
    k = G.Length()
    n = G[1].Length()
    Gs = G._matrix_(F)
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(Gs))
    #return gap.eval("BestKnownLinearCode(%s,%s,GF(%s))"%(n,k,q))

def best_known_linear_code_www(n, k, F, verbose=False):
    r"""
    Explains the construction of the best known linear code over GF(q) with
    length n and dimension k, courtesy of the www page
    http://www.codetables.de/.

    INPUT:

    -  ``n`` - Integer, the length of the code

    -  ``k`` - Integer, the dimension of the code

    -  ``F`` - Finite field, of order 2, 3, 4, 5, 7, 8, or 9

    -  ``verbose`` - Bool (default: ``False``)

    OUTPUT:


    -  Text about why the bounds are as given

    EXAMPLES::

        sage: L = best_known_linear_code_www(72, 36, GF(2)) # optional - internet
        sage: print L                                       # optional - internet
        Construction of a linear code
        [72,36,15] over GF(2):
        [1]:  [73, 36, 16] Cyclic Linear Code over GF(2)
             CyclicCode of length 73 with generating polynomial x^37 + x^36 + x^34 +
        x^33 + x^32 + x^27 + x^25 + x^24 + x^22 + x^21 + x^19 + x^18 + x^15 + x^11 +
        x^10 + x^8 + x^7 + x^5 + x^3 + 1
        [2]:  [72, 36, 15] Linear Code over GF(2)
             Puncturing of [1] at 1
        last modified: 2002-03-20

    This function raises an ``IOError`` if an error occurs downloading data or
    parsing it. It raises a ``ValueError`` if the ``q`` input is invalid.

    AUTHORS:

    - Steven Sivek (2005-11-14)
    - David Joyner (2008-03)
    """
    q = F.order()
    if not q in [2, 3, 4, 5, 7, 8, 9]:
        raise ValueError, "q (=%s) must be in [2,3,4,5,7,8,9]"%q
    n = int(n)
    k = int(k)

    param = ("?q=%s&n=%s&k=%s"%(q,n,k)).replace('L','')

    url = "http://iaks-www.ira.uka.de/home/grassl/codetables/BKLC/BKLC.php"+param
    #url = "http://homepages.cwi.nl/htbin/aeb/lincodbd/"+param
    if verbose:
        print "Looking up the bounds at %s"%url
    f = urllib.urlopen(url)
    s = f.read()
    f.close()
    #print s

    i = s.find("<PRE>")
    j = s.find("</PRE>")
    if i == -1 or j == -1:
        raise IOError, "Error parsing data (missing pre tags)."
    text = s[i+5:j].strip()
    return text

def bounds_minimum_distance(n, k, F):
    r"""
    Calculates a lower and upper bound for the minimum distance of an optimal
    linear code with word length ``n`` and dimension ``k`` over the field
    ``F``.

    The function returns a record with the two bounds and an explanation for
    each bound. The function Display can be used to show the explanations.

    The values for the lower and upper bound are obtained from a table
    constructed by Cen Tjhai for GUAVA, derived from the table of
    Brouwer. (See http://www.win.tue.nl/ aeb/voorlincod.html or use the
    Sage function ``minimum_distance_why`` for the most recent data.)
    These tables contain lower and upper bounds for `q=2` (when ``n <= 257``),
    `q=3` (when ``n <= 243``), `q=4` (``n <= 256``). (Current as of
    11 May 2006.) For codes over other fields and for larger word lengths,
    trivial bounds are used.

    This does not require an internet connection. The format of the output is
    a little non-intuitive. Try ``bounds_minimum_distance(10,5,GF(2))`` for
    an example.

    This function requires optional GAP package (Guava).

    EXAMPLES::

        sage: print bounds_minimum_distance(10,5,GF(2)) # optional - gap_packages (Guava package)
        rec(
          construction :=
           [ <Operation "ShortenedCode">,
              [
                  [ <Operation "UUVCode">,
                      [
                          [ <Operation "DualCode">,
                              [ [ <Operation "RepetitionCode">, [ 8, 2 ] ] ] ],
                          [ <Operation "UUVCode">,
                              [
                                  [ <Operation "DualCode">,
                                      [ [ <Operation "RepetitionCode">, [ 4, 2 ] ] ] ]
                                    , [ <Operation "RepetitionCode">, [ 4, 2 ] ] ] ]
                         ] ], [ 1, 2, 3, 4, 5, 6 ] ] ],
          k := 5,
          lowerBound := 4,
          lowerBoundExplanation := ...
          n := 10,
          q := 2,
          references := rec(
               ),
          upperBound := 4,
          upperBoundExplanation := ... )
    """
    q = F.order()
    gap.eval("data := BoundsMinimumDistance(%s,%s,GF(%s))"%(n,k,q))
    Ldata = gap.eval("Display(data)")
    return Ldata

def self_orthogonal_binary_codes(n, k, b=2, parent=None, BC=None, equal=False,
    in_test=None):
    """
    Returns a Python iterator which generates a complete set of
    representatives of all permutation equivalence classes of
    self-orthogonal binary linear codes of length in ``[1..n]`` and
    dimension in ``[1..k]``.

    INPUT:

    -  ``n`` - Integer, maximal length

    -  ``k`` - Integer, maximal dimension

    -  ``b`` - Integer, requires that the generators all have weight divisible
       by ``b`` (if ``b=2``, all self-orthogonal codes are generated, and if
       ``b=4``, all doubly even codes are generated). Must be an even positive
       integer.

    -  ``parent`` - Used in recursion (default: ``None``)

    -  ``BC`` - Used in recursion (default: ``None``)

    -  ``equal`` - If ``True`` generates only [n, k] codes (default: ``False``)

    -  ``in_test`` - Used in recursion (default: ``None``)

    EXAMPLES:

    Generate all self-orthogonal codes of length up to 7 and dimension up
    to 3::

        sage: for B in self_orthogonal_binary_codes(7,3):
        ...    print B
        ...
        Linear code of length 2, dimension 1 over Finite Field of size 2
        Linear code of length 4, dimension 2 over Finite Field of size 2
        Linear code of length 6, dimension 3 over Finite Field of size 2
        Linear code of length 4, dimension 1 over Finite Field of size 2
        Linear code of length 6, dimension 2 over Finite Field of size 2
        Linear code of length 6, dimension 2 over Finite Field of size 2
        Linear code of length 7, dimension 3 over Finite Field of size 2
        Linear code of length 6, dimension 1 over Finite Field of size 2

    Generate all doubly-even codes of length up to 7 and dimension up
    to 3::

        sage: for B in self_orthogonal_binary_codes(7,3,4):
        ...    print B; print B.gen_mat()
        ...
        Linear code of length 4, dimension 1 over Finite Field of size 2
        [1 1 1 1]
        Linear code of length 6, dimension 2 over Finite Field of size 2
        [1 1 1 1 0 0]
        [0 1 0 1 1 1]
        Linear code of length 7, dimension 3 over Finite Field of size 2
        [1 0 1 1 0 1 0]
        [0 1 0 1 1 1 0]
        [0 0 1 0 1 1 1]

    Generate all doubly-even codes of length up to 7 and dimension up
    to 2::

        sage: for B in self_orthogonal_binary_codes(7,2,4):
        ...    print B; print B.gen_mat()
        Linear code of length 4, dimension 1 over Finite Field of size 2
        [1 1 1 1]
        Linear code of length 6, dimension 2 over Finite Field of size 2
        [1 1 1 1 0 0]
        [0 1 0 1 1 1]

    Generate all self-orthogonal codes of length equal to 8 and
    dimension equal to 4::

        sage: for B in self_orthogonal_binary_codes(8, 4, equal=True):
        ...     print B; print B.gen_mat()
        Linear code of length 8, dimension 4 over Finite Field of size 2
        [1 0 0 1 0 0 0 0]
        [0 1 0 0 1 0 0 0]
        [0 0 1 0 0 1 0 0]
        [0 0 0 0 0 0 1 1]
        Linear code of length 8, dimension 4 over Finite Field of size 2
        [1 0 0 1 1 0 1 0]
        [0 1 0 1 1 1 0 0]
        [0 0 1 0 1 1 1 0]
        [0 0 0 1 0 1 1 1]

    Since all the codes will be self-orthogonal, b must be divisible by
    2::

        sage: list(self_orthogonal_binary_codes(8, 4, 1, equal=True))
        Traceback (most recent call last):
        ...
        ValueError: b (1) must be a positive even integer.
    """
    d=int(b)
    if d!=b or d%2==1 or d <= 0:
        raise ValueError("b (%s) must be a positive even integer."%b)
    from binary_code import BinaryCode, BinaryCodeClassifier
    if k < 1 or n < 2:
        return
    if equal:
        in_test = lambda M : (M.ncols() - M.nrows()) <= (n-k)
        out_test = lambda C : (C.dimension() == k) and (C.length() == n)
    else:
        in_test = lambda M : True
        out_test = lambda C : True
    if BC is None:
        BC = BinaryCodeClassifier()
    if parent is None:
        for j in xrange(d, n+1, d):
            M = Matrix(GF(2), [[1]*j])
            if in_test(M):
                for N in self_orthogonal_binary_codes(n, k, d, M, BC, in_test=in_test):
                    if out_test(N): yield N
    else:
        C = LinearCode(parent)
        if out_test(C): yield C
        if k == parent.nrows():
            return
        for nn in xrange(parent.ncols()+1, n+1):
            if in_test(parent):
                for child in BC.generate_children(BinaryCode(parent), nn, d):
                    for N in self_orthogonal_binary_codes(n, k, d, child, BC, in_test=in_test):
                        if out_test(N): yield N

############################ linear codes python class ########################

class LinearCode(module.Module_old):
    r"""
    A class for linear codes over a finite field or finite ring. Each instance
    is a linear code determined by a generator matrix `G` (i.e., a
    `k \times n` matrix of (full) rank `k`, `k \leq n` over a finite field `F`.

    INPUT:

    -  ``G`` - a generator matrix over `F`. (``G`` can be defined over a
       finite ring but the matrices over that ring must have certain
       attributes, such as ``rank``.)

    - ``d`` - (Optional, default: ``None``) the minimum distance of the
      code. This is an optional parameter.

    .. note::
        The veracity of the minimum distance ``d``, if provided, is not
        checked.

    OUTPUT:

    The linear code of length `n` over `F` having `G` as a generator matrix.

    EXAMPLES::

        sage: MS = MatrixSpace(GF(2),4,7)
        sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
        sage: C  = LinearCode(G)
        sage: C
        Linear code of length 7, dimension 4 over Finite Field of size 2
        sage: C.base_ring()
        Finite Field of size 2
        sage: C.dimension()
        4
        sage: C.length()
        7
        sage: C.minimum_distance()
        3
        sage: C.spectrum()
        [1, 0, 0, 7, 7, 0, 0, 1]
        sage: C.weight_distribution()
        [1, 0, 0, 7, 7, 0, 0, 1]

    The minimum distance of the code, if known, can be provided as an
    optional parameter.::

        sage: C  = LinearCode(G, d=3)
        sage: C.minimum_distance()
        3

    Another example.::

        sage: MS = MatrixSpace(GF(5),4,7)
        sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
        sage: C  = LinearCode(G)
        sage: C
        Linear code of length 7, dimension 4 over Finite Field of size 5

    AUTHORS:

    - David Joyner (11-2005)
    """
    #    sage: C.minimum_distance_upper_bound()   # optional (net connection)
    #    3
    #    sage: C.minimum_distance_why()     # optional (net connection)
    #    Ub(7,4) = 3 follows by the Griesmer bound.
    def __init__(self, gen_mat, d=None):
        r"""
        See the docstring for :meth:`LinearCode`.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)    # indirect doctest
            sage: C
            Linear code of length 7, dimension 4 over Finite Field of size 2

        The minimum distance of the code, if known, can be provided as an
        optional parameter.::

            sage: C  = LinearCode(G, d=3)
            sage: C.minimum_distance()
            3
        """
        base_ring = gen_mat[0,0].parent()
        ParentWithGens.__init__(self, base_ring)
        self.__gens = gen_mat.rows()
        self.__gen_mat = gen_mat
        self.__length = len(gen_mat.row(0))
        self.__dim = gen_mat.rank()
        self.__distance = d

    def _repr_(self):
        r"""
        See the docstring for :meth:`LinearCode`.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)
            sage: C                     # indirect doctest
            Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        return "Linear code of length %s, dimension %s over %s"%(self.length(), self.dimension(), self.base_ring())

    def automorphism_group_binary_code(self):
        r"""
        This function is deprecated. Use permutation_automorphism_group instead.

        This only applies to linear binary codes and returns its (permutation)
        automorphism group. In other words, if the code `C` has length `n`
        then it returns the subgroup of the symmetric group `S_n`:

        .. math::

           \{ g \in S_n\ |\ g(c) \in C, \forall c\in C\},

        where `S_n` acts on `GF(2)^n` by permuting coordinates.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: G = C.automorphism_group_binary_code(); G
            doctest:...: DeprecationWarning: This function is deprecated...
            Permutation Group with generators [(4,5)(6,7), (4,6)(5,7), (2,3)(6,7), (2,4)(3,5), (1,2)(5,6)]
            sage: G.order()
            168
        """
        from sage.misc.superseded import deprecation
        deprecation(11033, "This function is deprecated.  Call the method"
                    +" permutation_automorphism_group instead.")
        #deprecated 4.7
        C = self
        F = C.base_ring()
        if F!=GF(2):
            raise NotImplementedError, "Only implemented for binary codes."
        return self.permutation_automorphism_group()

    def __iter__(self):
        """
        Return an iterator over the elements of this linear code.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: [list(c) for c in C if c.hamming_weight() < 4]
            [[0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 1, 1],
             [0, 1, 0, 0, 1, 0, 1], [0, 0, 1, 0, 1, 1, 0],
             [1, 1, 1, 0, 0, 0, 0], [1, 0, 0, 1, 1, 0, 0],
             [0, 1, 0, 1, 0, 1, 0], [0, 0, 1, 1, 0, 0, 1]]
        """
        from sage.modules.finite_submodule_iter import \
                                                FiniteFieldsubspace_iterator
        return FiniteFieldsubspace_iterator(self.gen_mat())

    def ambient_space(self):
        r"""
        Returns the ambient vector space of `self`.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.ambient_space()
            Vector space of dimension 7 over Finite Field of size 2
        """
        return VectorSpace(self.base_ring(),self.__length)

    def assmus_mattson_designs(self, t, mode=None):
        r"""
        Assmus and Mattson Theorem (section 8.4, page 303 of [HP]): Let
        `A_0, A_1, ..., A_n` be the weights of the codewords in a binary
        linear `[n , k, d]` code `C`, and let `A_0^*, A_1^*, ..., A_n^*` be
        the weights of the codewords in its dual `[n, n-k, d^*]` code `C^*`.
        Fix a `t`, `0<t<d`, and let

        .. math::

           s = |\{ i\ |\ A_i^* \not= 0, 0< i \leq n-t\}|.

        Assume `s\leq d-t`.

        1. If `A_i\not= 0` and `d\leq i\leq n`
           then `C_i = \{ c \in C\ |\ wt(c) = i\}` holds a simple t-design.

        2. If `A_i^*\not= 0` and `d*\leq i\leq n-t` then
           `C_i^* = \{ c \in C^*\ |\ wt(c) = i\}` holds a simple t-design.

        A block design is a pair `(X,B)`, where `X` is a non-empty finite set
        of `v>0` elements called points, and `B` is a non-empty finite
        multiset of size b whose elements are called blocks, such that each
        block is a non-empty finite multiset of `k` points. `A` design without
        repeated blocks is called a simple block design. If every subset of
        points of size `t` is contained in exactly `\lambda` blocks the block
        design is called a `t-(v,k,\lambda)` design (or simply a `t`-design
        when the parameters are not specified). When `\lambda=1` then the
        block design is called a `S(t,k,v)` Steiner system.

        In the Assmus and Mattson Theorem (1), `X` is the set `\{1,2,...,n\}`
        of coordinate locations and `B = \{supp(c)\ |\ c \in C_i\}` is the set
        of supports of the codewords of `C` of weight `i`. Therefore, the
        parameters of the `t`-design for `C_i` are

        ::

            t =       given
            v =       n
            k =       i   (k not to be confused with dim(C))
            b =       Ai
            lambda = b*binomial(k,t)/binomial(v,t) (by Theorem 8.1.6,
                                                       p 294, in [HP])

        Setting the ``mode="verbose"`` option prints out the values of the
        parameters.

        The first example below means that the binary [24,12,8]-code C has
        the property that the (support of the) codewords of weight 8 (resp.,
        12, 16) form a 5-design. Similarly for its dual code `C^*` (of course
        `C=C^*` in this case, so this info is extraneous). The test fails to
        produce 6-designs (ie, the hypotheses of the theorem fail to hold,
        not that the 6-designs definitely don't exist). The command
        assmus_mattson_designs(C,5,mode="verbose") returns the same value
        but prints out more detailed information.

        The second example below illustrates the blocks of the 5-(24, 8, 1)
        design (i.e., the S(5,8,24) Steiner system).

        EXAMPLES::

            sage: C = codes.ExtendedBinaryGolayCode()             #  example 1
            sage: C.assmus_mattson_designs(5)
            ['weights from C: ',
            [8, 12, 16, 24],
            'designs from C: ',
            [[5, (24, 8, 1)], [5, (24, 12, 48)], [5, (24, 16, 78)], [5, (24, 24, 1)]],
            'weights from C*: ',
            [8, 12, 16],
            'designs from C*: ',
            [[5, (24, 8, 1)], [5, (24, 12, 48)], [5, (24, 16, 78)]]]
            sage: C.assmus_mattson_designs(6)
            0
            sage: X = range(24)                           #  example 2
            sage: blocks = [c.support() for c in C if c.hamming_weight()==8]; len(blocks)  # long time computation
            759

        REFERENCE:

        - [HP] W. C. Huffman and V. Pless, Fundamentals of ECC,
          Cambridge Univ. Press, 2003.
        """
        C = self
        ans = []
        G = C.gen_mat()
        n = len(G.columns())
        Cp = C.dual_code()
        wts = C.spectrum()
        d = min([i for i in range(1,len(wts)) if wts[i]!=0])
        if t>=d:
            return 0
        nonzerowts = [i for i in range(len(wts)) if wts[i]!=0 and i<=n and i>=d]
        #print d,t,len(nonzerowts)
        if mode=="verbose":
            #print "\n"
            for w in nonzerowts:
                print "The weight ",w," codewords of C form a t-(v,k,lambda) design, where"
                print "      t = ",t," , v = ",n," , k = ",w," , lambda = ",wts[w]*binomial(w,t)/binomial(n,t)#,"\n"
                print "      There are ",wts[w]," blocks of this design."
        wtsp = Cp.spectrum()
        dp = min([i for i in range(1,len(wtsp)) if wtsp[i]!=0])
        nonzerowtsp = [i for i in range(len(wtsp)) if wtsp[i]!=0 and i<=n-t and i>=dp]
        s = len([i for i in range(1,n) if wtsp[i]!=0 and i<=n-t and i>0])
        if mode=="verbose":
            #print "\n"
            for w in nonzerowtsp:
                print "The weight ",w," codewords of C* form a t-(v,k,lambda) design, where"
                print "      t = ",t," , v = ",n," , k = ",w," , lambda = ",wtsp[w]*binomial(w,t)/binomial(n,t)#,"\n"
                print "      There are ",wts[w]," blocks of this design."
        if s<=d-t:
            des = [[t,(n,w,wts[w]*binomial(w,t)/binomial(n,t))] for w in nonzerowts]
            ans = ans + ["weights from C: ",nonzerowts,"designs from C: ",des]
            desp = [[t,(n,w,wtsp[w]*binomial(w,t)/binomial(n,t))] for w in nonzerowtsp]
            ans = ans + ["weights from C*: ",nonzerowtsp,"designs from C*: ",desp]
            return ans
        return 0

    def basis(self):
        r"""
        Returns a basis of `self`.

        EXAMPLES::

            sage: C = codes.HammingCode(3, GF(2))
            sage: C.basis()
            [(1, 0, 0, 0, 0, 1, 1), (0, 1, 0, 0, 1, 0, 1), (0, 0, 1, 0, 1, 1, 0), (0, 0, 0, 1, 1, 1, 1)]
        """
        return self.__gens

    # S. Pancratz, 19 Jan 2010:  In the doctests below, I removed the example
    # ``C.binomial_moment(3)``, which was also marked as ``#long``.  This way,
    # we shorten the doctests time while still maintaining a zero and a
    # non-zero example.
    def binomial_moment(self, i):
        r"""
        Returns the i-th binomial moment of the `[n,k,d]_q`-code `C`:

        .. math::

            B_i(C) = \sum_{S, |S|=i} \frac{q^{k_S}-1}{q-1}

        where `k_S` is the dimension of the shortened code `C_{J-S}`,
        `J=[1,2,...,n]`. (The normalized binomial moment is
        `b_i(C) = \binom(n,d+i)^{-1}B_{d+i}(C)`.) In other words, `C_{J-S}`
        is isomorphic to the subcode of C of codewords supported on S.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.binomial_moment(2)
            0
            sage: C.binomial_moment(4)    # long time
            35

        .. warning::

            This is slow.

        REFERENCE:

        - I. Duursma, "Combinatorics of the two-variable zeta function",
          Finite fields and applications, 109-136, Lecture Notes in
          Comput. Sci., 2948, Springer, Berlin, 2004.
        """
        n = self.length()
        k = self.dimension()
        d = self.minimum_distance()
        F = self.base_ring()
        q = F.order()
        J = range(1,n+1)
        Cp = self.dual_code()
        dp = Cp.minimum_distance()
        if i<d:
            return 0
        if i>n-dp and i<=n:
            return binomial(n,i)*(q**(i+k-n) -1)/(q-1)
        P = SetPartitions(J,2).list()
        b = QQ(0)
        for p in P:
            p = list(p)
            S = p[0]
            if len(S)==n-i:
                C_S = self.shortened(S)
                k_S = C_S.dimension()
                b = b + (q**(k_S) -1)/(q-1)
        return b

    def __contains__(self,v):
        r"""
        Returns True if `v` can be coerced into `self`. Otherwise, returns False.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: vector((1, 0, 0, 0, 0, 1, 1)) in C   # indirect doctest
            True
            sage: vector((1, 0, 0, 0, 2, 1, 1)) in C   # indirect doctest
            True
            sage: vector((1, 0, 0, 0, 0, 1/2, 1)) in C # indirect doctest
            False
        """
        A = self.ambient_space()
        C = A.subspace(self.gens())
        return C.__contains__(v)

    def characteristic(self):
        r"""
        Returns the characteristic of the base ring of `self`.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.characteristic()
            2
        """
        return (self.base_ring()).characteristic()

    def characteristic_polynomial(self):
        r"""
        Returns the characteristic polynomial of a linear code, as defined in
        van Lint's text [vL].

        EXAMPLES::

            sage: C = codes.ExtendedBinaryGolayCode()
            sage: C.characteristic_polynomial()
            -4/3*x^3 + 64*x^2 - 2816/3*x + 4096

        REFERENCES:

        - van Lint, Introduction to coding theory, 3rd ed., Springer-Verlag
          GTM, 86, 1999.
        """
        R = PolynomialRing(QQ,"x")
        x = R.gen()
        C = self
        Cd = C.dual_code()
        Sd = Cd.support()
        k = C.dimension()
        n = C.length()
        q = (C.base_ring()).order()
        return q**(n-k)*prod([1-x/j for j in Sd if j>0])

    def chinen_polynomial(self):
        """
        Returns the Chinen zeta polynomial of the code.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.chinen_polynomial()       # long time
            1/5*(2*sqrt(2)*t^3 + 2*sqrt(2)*t^2 + 2*t^2 + sqrt(2)*t + 2*t + 1)/(sqrt(2) + 1)
            sage: C = codes.TernaryGolayCode()
            sage: C.chinen_polynomial()       # long time
            1/7*(3*sqrt(3)*t^3 + 3*sqrt(3)*t^2 + 3*t^2 + sqrt(3)*t + 3*t + 1)/(sqrt(3) + 1)

        This last output agrees with the corresponding example given in
        Chinen's paper below.

        REFERENCES:

        - Chinen, K. "An abundance of invariant polynomials satisfying the
          Riemann hypothesis", April 2007 preprint.
        """
        from sage.functions.all import sqrt
        C = self
        n = C.length()
        RT = PolynomialRing(QQ,2,"Ts")
        T,s = FractionField(RT).gens()
        t = PolynomialRing(QQ,"t").gen()
        Cd = C.dual_code()
        k = C.dimension()
        q = (C.base_ring()).characteristic()
        d = C.minimum_distance()
        dperp = Cd.minimum_distance()
        if dperp > d:
            P = RT(C.zeta_polynomial())
            # Sage does not find dealing with sqrt(int) *as an algebraic object*
            # an easy thing to do. Some tricky gymnastics are used to
            # make Sage deal with objects over QQ(sqrt(q)) nicely.
            if is_even(n):
                Pd = q**(k-n/2)*RT(Cd.zeta_polynomial())*T**(dperp - d)
            if not(is_even(n)):
                Pd = s*q**(k-(n+1)/2)*RT(Cd.zeta_polynomial())*T**(dperp - d)
            CP = P+Pd
            f = CP/CP(1,s)
            return f(t,sqrt(q))
        if dperp < d:
            P = RT(C.zeta_polynomial())*T**(d - dperp)
            if is_even(n):
                Pd = q**(k-n/2)*RT(Cd.zeta_polynomial())
            if not(is_even(n)):
                Pd = s*q**(k-(n+1)/2)*RT(Cd.zeta_polynomial())
            CP = P+Pd
            f = CP/CP(1,s)
            return f(t,sqrt(q))
        if dperp == d:
            P = RT(C.zeta_polynomial())
            if is_even(n):
                Pd = q**(k-n/2)*RT(Cd.zeta_polynomial())
            if not(is_even(n)):
                Pd = s*q**(k-(n+1)/2)*RT(Cd.zeta_polynomial())
            CP = P+Pd
            f = CP/CP(1,s)
            return f(t,sqrt(q))

    def __cmp__(self, right):
        r"""
        Returns True if the generator matrices of `self` and `right` are
        equal.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G = MS([1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1])
            sage: G
            [1 0 0 0 0 1 1]
            [0 1 0 0 1 0 1]
            [0 0 1 0 1 1 0]
            [0 0 0 1 1 1 1]
            sage: D = LinearCode(G)
            sage: C == D
            True

            sage: Cperp = C.dual_code()
            sage: Cperpperp = Cperp.dual_code()
            sage: C == Cperpperp
            True
        """
        if not isinstance(right, LinearCode):
            return cmp(type(self), type(right))
        return cmp(self.__gen_mat, right.__gen_mat)

    def check_mat(self):
        r"""
        Returns the check matrix of ``self``.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: Cperp = C.dual_code()
            sage: C; Cperp
            Linear code of length 7, dimension 4 over Finite Field of size 2
            Linear code of length 7, dimension 3 over Finite Field of size 2
            sage: C.gen_mat()
             [1 0 0 0 0 1 1]
             [0 1 0 0 1 0 1]
             [0 0 1 0 1 1 0]
             [0 0 0 1 1 1 1]
            sage: C.check_mat()
             [1 0 1 0 1 0 1]
             [0 1 1 0 0 1 1]
             [0 0 0 1 1 1 1]
            sage: Cperp.check_mat()
             [1 0 0 0 0 1 1]
             [0 1 0 0 1 0 1]
             [0 0 1 0 1 1 0]
             [0 0 0 1 1 1 1]
            sage: Cperp.gen_mat()
             [1 0 1 0 1 0 1]
             [0 1 1 0 0 1 1]
             [0 0 0 1 1 1 1]
        """
        Cperp = self.dual_code()
        return Cperp.gen_mat()

    def covering_radius(self):
        r"""
        Wraps Guava's ``CoveringRadius`` command.

        The covering radius of a linear code `C` is the smallest number `r`
        with the property that each element `v` of the ambient vector space
        of `C` has at most a distance `r` to the code `C`. So for each
        vector `v` there must be an element `c` of `C` with `d(v,c) \leq  r`.
        A binary linear code with reasonable small covering radius is often
        referred to as a covering code.

        For example, if `C` is a perfect code, the covering radius is equal
        to `t`, the number of errors the code can correct, where `d = 2t+1`,
        with `d` the minimum distance of `C`.

        EXAMPLES::

            sage: C = codes.HammingCode(5,GF(2))
            sage: C.covering_radius()  # optional - gap_packages (Guava package)
            1
        """
        F = self.base_ring()
        G = self.gen_mat()
        gapG = gap(G)
        C = gapG.GeneratorMatCode(gap(F))
        r = C.CoveringRadius()
        try:
            return ZZ(r)
        except TypeError:
            raise RuntimeError("the covering radius of this code cannot be computed by Guava")

    @rename_keyword(deprecation=11033, method="algorithm")
    def decode(self, right, algorithm="syndrome"):
        r"""
        Decodes the received vector ``right`` to an element `c` in this code.

        Optional algorithms are "guava", "nearest neighbor" or "syndrome". The
        ``algorithm="guava"`` wraps GUAVA's ``Decodeword``.  Hamming codes have
        a special decoding algorithm; otherwise, ``"syndrome"`` decoding is
        used.

        INPUT:

        - ``right`` - Vector of length the length of this code
        - ``algorithm`` - Algorithm to use, one of ``"syndrome"``, ``"nearest
          neighbor"``, and ``"guava"`` (default: ``"syndrome"``)

        OUTPUT:

        - The codeword in this code closest to ``right``.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: MS = MatrixSpace(GF(2),1,7)
            sage: F = GF(2); a = F.gen()
            sage: v1 = [a,a,F(0),a,a,F(0),a]
            sage: C.decode(v1)
            (1, 1, 0, 1, 0, 0, 1)
            sage: C.decode(v1,algorithm="nearest neighbor")
            (1, 1, 0, 1, 0, 0, 1)
            sage: C.decode(v1,algorithm="guava")  # optional - gap_packages (Guava package)
            (1, 1, 0, 1, 0, 0, 1)
            sage: v2 = matrix([[a,a,F(0),a,a,F(0),a]])
            sage: C.decode(v2)
            (1, 1, 0, 1, 0, 0, 1)
            sage: v3 = vector([a,a,F(0),a,a,F(0),a])
            sage: c = C.decode(v3); c
            (1, 1, 0, 1, 0, 0, 1)
            sage: c in C
            True
            sage: C = codes.HammingCode(2,GF(5))
            sage: v = vector(GF(5),[1,0,0,2,1,0])
            sage: C.decode(v)
            (1, 0, 0, 2, 2, 0)
            sage: F.<a> = GF(4)
            sage: C = codes.HammingCode(2,F)
            sage: v = vector(F, [1,0,0,a,1])
            sage: C.decode(v)
            (a + 1, 0, 0, a, 1)
            sage: C.decode(v, algorithm="nearest neighbor")
            (a + 1, 0, 0, a, 1)
            sage: C.decode(v, algorithm="guava")  # optional - gap_packages (Guava package)
            (a + 1, 0, 0, a, 1)

        Does not work for very long codes since the syndrome table grows too
        large.
        """
        from decoder import decode
        if algorithm == 'syndrome' or algorithm == 'nearest neighbor':
            return decode(self,right)
        elif algorithm == 'guava':
            gap.LoadPackage('"guava"')
            code = gap.GeneratorMatCode(self.gen_mat(), self.base_ring())
            right = gap(list(right))
            right_word = gap.Codeword(right)
            result = gap.Decodeword(code, right_word)
            result = gap.VectorCodeword(result)
            from sage.interfaces.gap import gfq_gap_to_sage
            result = [gfq_gap_to_sage(v, self.base_ring()) for v in result]
            return self.ambient_space()(result)
        else:
            raise NotImplementedError, "Only 'syndrome','nearest neighbor','guava' are implemented."

    def divisor(self):
        r"""
        Returns the divisor of a code, which is the smallest integer `d_0 > 0`
        such that each `A_i > 0` iff `i` is divisible by `d_0`.

        EXAMPLES::

            sage: C = codes.ExtendedBinaryGolayCode()
            sage: C.divisor()   # Type II self-dual
            4
            sage: C = codes.QuadraticResidueCodeEvenPair(17,GF(2))[0]
            sage: C.divisor()
            2
        """
        C = self
        A = C.spectrum()
        n = C.length()
        V = VectorSpace(QQ,n+1)
        S = V(A).nonzero_positions()
        S0 = [S[i] for i in range(1,len(S))]
        #print S0
        if len(S)>1: return GCD(S0)
        return 1

    def dual_code(self):
        r"""
        This computes the dual code `Cd` of the code `C`,

        .. math::

            Cd = \{ v \in V\ |\ v\cdot c = 0,\ \forall c \in C \}.

        Does not call GAP.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.dual_code()
            Linear code of length 7, dimension 3 over Finite Field of size 2
            sage: C = codes.HammingCode(3,GF(4,'a'))
            sage: C.dual_code()
            Linear code of length 21, dimension 3 over Finite Field in a of size 2^2
        """
        G = self.gen_mat()
        H = G.transpose().kernel()
        V = H.ambient_vector_space()
        Cd = LinearCodeFromVectorSpace(V.span(H))
        return Cd
        #another way:
        #Gsf, p = standard_form(G)
        #k = len(G.rows())
        #n = len(G.columns())
        #MS = G.parent()
        #sG = G.matrix_from_columns(range(k,n))
        #Inmk = MatrixSpace(F,n-k,n-k).identity_matrix()
        #H = Inmk.augment(sG.transpose())
        #return LinearCode(H)

    def dimension(self):
        r"""
        Returns the dimension of this code.

        EXAMPLES::

            sage: G = matrix(GF(2),[[1,0,0],[1,1,0]])
            sage: C = LinearCode(G)
            sage: C.dimension()
            2
        """
        return self.__dim

    def direct_sum(self, other):
        """
        Returns the code given by the direct sum of the codes ``self`` and
        ``other``, which must be linear codes defined over the same base ring.

        EXAMPLES::

            sage: C1 = codes.HammingCode(3,GF(2))
            sage: C2 = C1.direct_sum(C1); C2
            Linear code of length 14, dimension 8 over Finite Field of size 2
            sage: C3 = C1.direct_sum(C2); C3
            Linear code of length 21, dimension 12 over Finite Field of size 2
        """
        C1 = self; C2 = other
        G1 = C1.gen_mat()
        G2 = C2.gen_mat()
        F = C1.base_ring()
        n1 = len(G1.columns())
        k1 = len(G1.rows())
        n2 = len(G2.columns())
        k2 = len(G2.rows())
        MS1 = MatrixSpace(F,k2,n1)
        MS2 = MatrixSpace(F,k1,n2)
        Z1 = MS1(0)
        Z2 = MS2(0)
        top = G1.augment(Z2)
        bottom = Z1.augment(G2)
        G = top.stack(bottom)
        return LinearCode(G)

    def __eq__(self, right):
        """
        Checks if ``self`` is equal to ``right``.

        EXAMPLES::

            sage: C1 = codes.HammingCode(3,GF(2))
            sage: C2 = codes.HammingCode(3,GF(2))
            sage: C1 == C2
            True
            sage: C2 = C1.extended_code()
            sage: C3 = C2.punctured([7])
            sage: C1 == C3
            True
        """
        slength = self.length()
        rlength = right.length()
        sdim = self.dimension()
        rdim = right.dimension()
        sF = self.base_ring()
        rF = right.base_ring()
        if slength != rlength:
            return False
        if sdim != rdim:
            return False
        if sF != rF:
            return False
        sbasis = self.gens()
        rbasis = right.gens()
        scheck = self.check_mat()
        rcheck = right.check_mat()
        for c in sbasis:
            if rcheck*c:
                return False
        for c in rbasis:
            if scheck*c:
                return False
        return True

    def extended_code(self):
        r"""
        If ``self`` is a linear code of length `n` defined over `F` then this
        returns the code of length `n+1` where the last digit `c_n` satisfies
        the check condition `c_0+...+c_n=0`. If ``self`` is an `[n,k,d]`
        binary code then the extended code `C^{\vee}` is an `[n+1,k,d^{\vee}]`
        code, where `d^=d` (if d is even) and `d^{\vee}=d+1` (if `d` is odd).

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(4,'a'))
            sage: C
            Linear code of length 21, dimension 18 over Finite Field in a of size 2^2
            sage: Cx = C.extended_code()
            sage: Cx
            Linear code of length 22, dimension 18 over Finite Field in a of size 2^2
        """
        G = self.gen_mat()
        F = self.base_ring()
        k = len(G.rows())
        MS1 = MatrixSpace(F,k,1)
        ck_sums = [-sum(G.rows()[i]) for i in range(k)]
        last_col = MS1(ck_sums)
        Gx = G.augment(last_col)
        return LinearCode(Gx)

    def galois_closure(self, F0):
        r"""
        If ``self`` is a linear code defined over `F` and `F_0` is a subfield
        with Galois group `G = Gal(F/F_0)` then this returns the `G`-module
        `C^-` containing `C`.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(4,'a'))
            sage: Cc = C.galois_closure(GF(2))
            sage: C; Cc
            Linear code of length 21, dimension 18 over Finite Field in a of size 2^2
            Linear code of length 21, dimension 20 over Finite Field in a of size 2^2
            sage: c = C.basis()[2]
            sage: V = VectorSpace(GF(4,'a'),21)
            sage: c2 = V([x^2 for x in c.list()])
            sage: c2 in C
            False
            sage: c2 in Cc
            True
        """
        G = self.gen_mat()
        F = self.base_ring()
        q = F.order()
        q0 = F0.order()
        a = log(q,q0)  # test if F/F0 is a field extension
        if not isinstance(a, Integer):
            raise ValueError,"Base field must be an extension of given field %s"%F0
        n = len(G.columns())
        k = len(G.rows())
        G0 = [[x**q0 for x in g.list()] for g in G.rows()]
        G1 = [[x for x in g.list()] for g in G.rows()]
        G2 = G0+G1
        MS = MatrixSpace(F,2*k,n)
        G3 = MS(G2)
        r = G3.rank()
        MS = MatrixSpace(F,r,n)
        Grref = G3.echelon_form()
        G = MS([Grref.row(i) for i in range(r)])
        return LinearCode(G)

    def __getitem__(self, i):
        r"""
        Returns the `i`-th codeword of this code.

        The implementation of this depends on the implementation of the
        :meth:`.__iter__` method.

        The implementation is as follows. Suppose that:

        - the primitive element of the base_ring of this code is `a`,
        - the prime subfield is `p`,
        - the field has order `p^m`,
        - the code has dimension `k`,
        - and the generator matrix is `G`.

        Then the :meth:`.__iter__` method returns the elements in this order:

        1. first, the following ordered list is returned:
           ``[i*a^0 * G[0] for i in range(p)]``
        2. Next, the following ordered list is returned:
           ``[i*a^0 * G[0] + a^1*G[0] for i in range(p)]``
        3. This continues till we get
           ``[(i*a^0 +(p-1)*a^1 +...+ (p-1)*a^(m-1))*G[0] for i in range(p)]``
        4. Then, we move to G[1]:
           ``[i*a^0 * G[0] + a^0*G[1] for i in range(p)]``,
         and so on.
         Hence the `i`-th element can be obtained by the p-adic expansion
         of `i` as ``[i_0, i_1, ...,i_{m-1}, i_m, i_{m+1}, ..., i_{km-1}].``
         The element that is generated is:

        .. math::

             \begin{aligned}
             & (i_0 a^0 + i_1 a^1 + \cdots + i_{m-1} a^{m-1}) G[0] + \\
             & (i_m a^0 + i_{m+1} a^1 + \cdots + i_{2m-1} a^{m-1}) G[1] + \\
             & \vdots\\
             & (i_{(k-1)m} a^0 + \cdots + i_{km-1} a^{m-1}) G[k-1]
             \end{aligned}

        EXAMPLES::

            sage: RS = codes.ReedSolomonCode(7, 3, GF(8, 'a'))
            sage: RS[24]
            (0, a^2 + a, a^2 + a + 1, a^2 + 1, 1, a, a^2)
            sage: RS[24] == RS.list()[24]
            True

        TESTS::

            sage: C = random_matrix(GF(25,'a'), 2, 7).row_space()
            sage: C = LinearCode(C.basis_matrix())
            sage: Clist = C.list()
            sage: all([C[i]==Clist[i] for i in xrange(len(C))])
            True

        Check that only the indices less than the size of the code are
        allowed::

            sage: C[25**2]
            Traceback (most recent call last):
            ...
            IndexError: The value of the index 'i' (=625) must be between
            0 and 'q^k -1' (=624), inclusive, where 'q' is the size of the
            base field and 'k' is the dimension of the code.

        """
        # IMPORTANT: If the __iter__() function implementation is changed
        # then the implementation here must also be changed so that
        # list(self)[i] and self[i] both return the same element.

        from sage.rings.padics.factory import Zp
        F = self.base_ring()
        maxindex = F.order()**self.dimension()-1
        if i < 0 or i > maxindex:
            raise IndexError("The value of the index 'i' (={}) must be between "
                             "0 and 'q^k -1' (={}), inclusive, where 'q' is "
                             "the size of the base field and 'k' is the "
                             "dimension of the code.".format(i, maxindex))

        a = F.primitive_element()
        m = F.degree()
        p = F.prime_subfield().order()
        A = [a**k for k in xrange(m)]
        G = self.gen_mat()
        N = self.dimension()*F.degree() # the total length of p-adic vector
        Z = Zp(p, N)
        ivec = Z(i).padded_list(N)

        codeword = 0
        row = 0
        for g in G:
            codeword += sum([ivec[j+row*m]*A[j] for j in xrange(m)])*g
            row += 1
        return codeword

    def gen_mat(self):
        r"""
        Return a generator matrix of this code.

        EXAMPLES::

            sage: C1 = codes.HammingCode(3,GF(2))
            sage: C1.gen_mat()
            [1 0 0 0 0 1 1]
            [0 1 0 0 1 0 1]
            [0 0 1 0 1 1 0]
            [0 0 0 1 1 1 1]
            sage: C2 = codes.HammingCode(2,GF(4,"a"))
            sage: C2.gen_mat()
            [    1     0     0 a + 1     a]
            [    0     1     0     1     1]
            [    0     0     1     a a + 1]
        """
        return self.__gen_mat

    def gen_mat_systematic(self):
        """
        Return a systematic generator matrix of the code.

        A generator matrix of a code is called systematic if it contains
        a set of columns forming an identity matrix.

        EXAMPLES::

            sage: G = matrix(GF(3),2,[1,-1,1,-1,1,1])
            sage: code = LinearCode(G)
            sage: code.gen_mat()
            [1 2 1]
            [2 1 1]
            sage: code.gen_mat_systematic()
            [1 2 0]
            [0 0 1]
        """
        return self.__gen_mat.echelon_form()

    def gens(self):
        r"""
        Returns the generators of this code as a list of vectors.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.gens()
             [(1, 0, 0, 0, 0, 1, 1), (0, 1, 0, 0, 1, 0, 1), (0, 0, 1, 0, 1, 1, 0), (0, 0, 0, 1, 1, 1, 1)]
        """
        return self.__gens

    def genus(self):
        r"""
        Returns the "Duursma genus" of the code, `\gamma_C = n+1-k-d`.

        EXAMPLES::

            sage: C1 = codes.HammingCode(3,GF(2)); C1
            Linear code of length 7, dimension 4 over Finite Field of size 2
            sage: C1.genus()
            1
            sage: C2 = codes.HammingCode(2,GF(4,"a")); C2
            Linear code of length 5, dimension 3 over Finite Field in a of size 2^2
            sage: C2.genus()
            0

        Since all Hamming codes have minimum distance 3, these computations
        agree with the definition, `n+1-k-d`.
        """
        d = self.minimum_distance()
        n = self.length()
        k = self.dimension()
        gammaC = n+1-k-d
        return gammaC

    def information_set(self):
        """
        Return an information set of the code.

        A set of column positions of a generator matrix of a code
        is called an information set if the corresponding columns
        form a square matrix of full rank.

        OUTPUT:

        - Information set of a systematic generator matrix of the code.

        EXAMPLES::

            sage: G = matrix(GF(3),2,[1,-1,0,-1,1,1])
            sage: code = LinearCode(G)
            sage: code.gen_mat_systematic()
            [1 2 0]
            [0 0 1]
            sage: code.information_set()
            (0, 2)
        """
        return self.__gen_mat.transpose().pivot_rows()

    def is_permutation_automorphism(self,g):
        r"""
        Returns `1` if `g` is an element of `S_n` (`n` = length of self) and
        if `g` is an automorphism of self.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(3))
            sage: g = SymmetricGroup(13).random_element()
            sage: C.is_permutation_automorphism(g)
            0
            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]])
            sage: C  = LinearCode(G)
            sage: S8 = SymmetricGroup(8)
            sage: g = S8("(2,3)")
            sage: C.is_permutation_automorphism(g)
            1
            sage: g = S8("(1,2,3,4)")
            sage: C.is_permutation_automorphism(g)
            0
        """
        basis = self.gen_mat().rows()
        H = self.check_mat()
        V = H.column_space()
        HGm = H*g.matrix()
        # raise TypeError, (type(H), type(V), type(basis[0]), type(Gmc))
        for c in basis:
            if HGm*c != V(0):
                return False
        return True

    @rename_keyword(deprecation=11033, method="algorithm")
    def is_permutation_equivalent(self,other,algorithm=None):
        """
        Returns ``True`` if ``self`` and ``other`` are permutation equivalent
        codes and ``False`` otherwise.

        The ``algorithm="verbose"`` option also returns a permutation (if
        ``True``) sending ``self`` to ``other``.

        Uses Robert Miller's double coset partition refinement work.

        EXAMPLES::

            sage: P.<x> = PolynomialRing(GF(2),"x")
            sage: g = x^3+x+1
            sage: C1 = codes.CyclicCodeFromGeneratingPolynomial(7,g); C1
            Linear code of length 7, dimension 4 over Finite Field of size 2
            sage: C2 = codes.HammingCode(3,GF(2)); C2
            Linear code of length 7, dimension 4 over Finite Field of size 2
            sage: C1.is_permutation_equivalent(C2)
            True
            sage: C1.is_permutation_equivalent(C2,algorithm="verbose")
            (True, (3,4)(5,7,6))
            sage: C1 = codes.RandomLinearCode(10,5,GF(2))
            sage: C2 = codes.RandomLinearCode(10,5,GF(3))
            sage: C1.is_permutation_equivalent(C2)
            False
        """
        from sage.groups.perm_gps.partn_ref.refinement_binary import NonlinearBinaryCodeStruct
        F = self.base_ring()
        F_o = other.base_ring()
        q = F.order()
        G = self.gen_mat()
        n = self.length()
        n_o = other.length()
        if F != F_o or n != n_o:
            return False
        k = len(G.rows())
        MS = MatrixSpace(F,q**k,n)
        CW1 = MS(self.list())
        CW2 = MS(other.list())
        B1 = NonlinearBinaryCodeStruct(CW1)
        B2 = NonlinearBinaryCodeStruct(CW2)
        ans = B1.is_isomorphic(B2)
        if ans!=False:
            if algorithm=="verbose":
                Sn = SymmetricGroup(n)
                return True, Sn([i+1 for i in ans])**(-1)
            return True
        return False

    def is_self_dual(self):
        """
        Returns ``True`` if the code is self-dual (in the usual Hamming inner
        product) and ``False`` otherwise.

        EXAMPLES::

            sage: C = codes.ExtendedBinaryGolayCode()
            sage: C.is_self_dual()
            True
            sage: C = codes.HammingCode(3,GF(2))
            sage: C.is_self_dual()
            False
        """
        return self == self.dual_code()


    def is_self_orthogonal(self):
        """
        Returns ``True`` if this code is self-orthogonal and ``False``
        otherwise.

        A code is self-orthogonal if it is a subcode of its dual.

        EXAMPLES::

            sage: C = codes.ExtendedBinaryGolayCode()
            sage: C.is_self_orthogonal()
            True
            sage: C = codes.HammingCode(3,GF(2))
            sage: C.is_self_orthogonal()
            False
            sage: C = QuasiQuadraticResidueCode(11)  # optional - gap_packages (Guava package)
            sage: C.is_self_orthogonal()             # optional - gap_packages (Guava package)
            True
        """
        return self.is_subcode(self.dual_code())

    def is_galois_closed(self):
        r"""
        Checks if ``self`` is equal to its Galois closure.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(4,"a"))
            sage: C.is_galois_closed()
            False
        """
        F = self.base_ring()
        p = F.characteristic()
        return self == self.galois_closure(GF(p))

    def is_subcode(self, other):
        """
        Returns ``True`` if ``self`` is a subcode of ``other``.

        EXAMPLES::

            sage: C1 = codes.HammingCode(3,GF(2))
            sage: G1 = C1.gen_mat()
            sage: G2 = G1.matrix_from_rows([0,1,2])
            sage: C2 = LinearCode(G2)
            sage: C2.is_subcode(C1)
            True
            sage: C1.is_subcode(C2)
            False
            sage: C3 = C1.extended_code()
            sage: C1.is_subcode(C3)
            False
            sage: C4 = C1.punctured([1])
            sage: C4.is_subcode(C1)
            False
            sage: C5 = C1.shortened([1])
            sage: C5.is_subcode(C1)
            False
            sage: C1 = codes.HammingCode(3,GF(9,"z"))
            sage: G1 = C1.gen_mat()
            sage: G2 = G1.matrix_from_rows([0,1,2])
            sage: C2 = LinearCode(G2)
            sage: C2.is_subcode(C1)
            True
        """
        G = self.gen_mat()
        for r in G.rows():
            if not(r in other):
                return False
        return True

    def __len__(self):
        r"""
        Return the size of this code.

        EXAMPLES::

            sage: C = codes.HammingCode(3, GF(2))
            sage: len(C)
            16
        """
        return self.base_ring().order()**self.dimension()

    def length(self):
        r"""
        Returns the length of this code.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.length()
            7
        """
        return self.__length

    def list(self):
        r"""
        Return a list of all elements of this linear code.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: Clist = C.list()
            sage: Clist[5]; Clist[5] in C
            (1, 0, 1, 0, 1, 0, 1)
            True
        """
        return list(self.__iter__())

    def _magma_init_(self, magma):
        r"""
        Retun a string representation in Magma of this linear code.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: Cm = magma(C)                 # optional - magma, indirect doctest
            sage: Cm.MinimumWeight()            # optional - magma
            3

        """
        G = magma(self.gen_mat())._ref()
        s = 'LinearCode(%s)' % G
        return s

    @rename_keyword(deprecation=11033, method="algorithm")
    def minimum_distance(self, algorithm=None):
        r"""
        Returns the minimum distance of this linear code.

        By default, this uses a GAP kernel function (in C and not part of
        Guava) written by Steve Linton.  If ``algorithm="guava"`` is set  and
        `q` is 2 or 3 then this uses a very fast program written in C written
        by CJ Tjhal. (This is much faster, except in some small examples.)

        Raises a ``ValueError`` in case there is no non-zero vector in this
        linear code.

        The minimum distance of the code is stored once it has been
        computed or provided during the initialization of :class:`LinearCode`.
        If ``algorithm`` is ``None`` and the stored value of minimum
        distance is found, then the stored value will be returned without
        recomputing the minimum distance again.

        INPUT:

        - ``algorithm`` - Method to be used, ``None``, ``"gap"``, or
          ``"guava"`` (default: ``None``).

        OUTPUT:

        - Integer, minimum distance of this code

        EXAMPLES::

            sage: MS = MatrixSpace(GF(3),4,7)
            sage: G = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.minimum_distance()
            3

        Once the minimum distance has been computed, it's value is stored.
        Hence the following command will return the value instantly,
        without further computations.::

            sage: C.minimum_distance()
            3

        If ``algorithm`` is provided, then the minimum distance will be
        recomputed even if there is a stored value from a previous run.::

            sage: C.minimum_distance(algorithm="gap")
            3
            sage: C.minimum_distance(algorithm="guava")  # optional - gap_packages (Guava package)
            3

        Another example.::

            sage: C = codes.HammingCode(2,GF(4,"a")); C
            Linear code of length 5, dimension 3 over Finite Field in a of size 2^2
            sage: C.minimum_distance()
            3

        TESTS::

            sage: C = codes.HammingCode(2,GF(4,"a"))
            sage: C.minimum_distance(algorithm='something')
            Traceback (most recent call last):
            ...
            ValueError: The algorithm argument must be one of None, 'gap' or 'guava'; got 'something'

        This shows that ticket :trac:`#6486` has been resolved::

            sage: G = matrix(GF(2),[[0,0,0]])
            sage: C = LinearCode(G)
            sage: C.minimum_distance()
            Traceback (most recent call last):
            ...
            ValueError: this linear code contains no non-zero vector
        """
        # Special code to handle the case where there is no non-zero vector.
        if self.dimension() == 0:
            raise ValueError("this linear code contains no non-zero vector")

        # If the minimum distance has already been computed or provided by
        # the user then simply return the stored value.
        # This is done only if algorithm is None.
        if self.__distance is not None and algorithm is None:
            return self.__distance
        if algorithm not in (None, "gap", "guava"):
            raise ValueError("The algorithm argument must be one of None, "
                        "'gap' or 'guava'; got '{0}'".format(algorithm))

        #sage: C.minimum_distance_upper_bound()  # optional (net connection)
        #5
        #    sage: C.minimum_distance_why()          # optional (net connection)
        #    Ub(10,5) = 5 follows by the Griesmer bound.
        F = self.base_ring()
        q = F.order()
        G = self.gen_mat()
        n = self.length()
        k = self.dimension()
        gapG = gap(G)
        if (q == 2 or q == 3) and algorithm=="guava":
            C = gapG.GeneratorMatCode(gap(F))
            d = C.MinimumWeight()
            #print "Running Guava's MinimumWeight ...\n"
            return ZZ(d)
        Gstr = "%s*Z(%s)^0"%(gapG, q)
        self.__distance = min_wt_vec_gap(Gstr,n,k,F).hamming_weight()
        return self.__distance

    def module_composition_factors(self, gp):
        r"""
        Prints the GAP record of the Meataxe composition factors module in
        Meataxe notation. This uses GAP but not Guava.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]])
            sage: C  = LinearCode(G)
            sage: gp = C.permutation_automorphism_group()

        Now type "C.module_composition_factors(gp)" to get the record printed.
        """
        F = self.base_ring()
        q = F.order()
        gens = gp.gens()
        G = self.gen_mat()
        n = len(G.columns())
        MS = MatrixSpace(F,n,n)
        mats = [] # initializing list of mats by which the gens act on self
        Fn = VectorSpace(F, n)
        W = Fn.subspace_with_basis(G.rows()) # this is self
        for g in gens:
            p = MS(g.matrix())
            m = [W.coordinate_vector(r*p) for r in G.rows()]
            mats.append(m)
        mats_str = str(gap([[list(r) for r in m] for m in mats]))
        gap.eval("M:=GModuleByMats("+mats_str+", GF("+str(q)+"))")
        print gap("MTX.CompositionFactors( M )")

    @rename_keyword(deprecation=11033, method="algorithm")
    def permutation_automorphism_group(self, algorithm="partition"):
        r"""
        If `C` is an `[n,k,d]` code over `F`, this function computes the
        subgroup `Aut(C) \subset S_n` of all permutation automorphisms of `C`.
        The binary case always uses the (default) partition refinement
        algorithm of Robert Miller.

        Note that if the base ring of `C` is `GF(2)` then this is the full
        automorphism group.

        INPUT:

        - ``algorithm`` - If ``"gap"`` then GAP's MatrixAutomorphism function
          (written by Thomas Breuer) is used. The implementation combines an
          idea of mine with an improvement suggested by Cary Huffman. If
          ``"gap+verbose"`` then code-theoretic data is printed out at
          several stages of the computation. If ``"partition"`` then the
          (default) partition refinement algorithm of Robert Miller is used.

        OUTPUT:

        - Permutation automorphism group

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]])
            sage: C  = LinearCode(G)
            sage: C
            Linear code of length 8, dimension 4 over Finite Field of size 2
            sage: G = C.permutation_automorphism_group()
            sage: G.order()
            144

        A less easy example involves showing that the permutation
        automorphism group of the extended ternary Golay code is the
        Mathieu group `M_{11}`.

        ::

            sage: C = codes.ExtendedTernaryGolayCode()
            sage: M11 = MathieuGroup(11)
            sage: M11.order()
            7920
            sage: G = C.permutation_automorphism_group()  # long time (6s on sage.math, 2011)
            sage: G.is_isomorphic(M11)                    # long time
            True

        Other examples::

            sage: C = codes.ExtendedBinaryGolayCode()
            sage: G = C.permutation_automorphism_group()
            sage: G.order()
            244823040
            sage: C = codes.HammingCode(5, GF(2))
            sage: G = C.permutation_automorphism_group()
            sage: G.order()
            9999360
            sage: C = codes.HammingCode(2,GF(3)); C
            Linear code of length 4, dimension 2 over Finite Field of size 3
            sage: C.permutation_automorphism_group(algorithm="partition")
            Permutation Group with generators [(1,3,4)]
            sage: C = codes.HammingCode(2,GF(4,"z")); C
            Linear code of length 5, dimension 3 over Finite Field in z of size 2^2
            sage: C.permutation_automorphism_group(algorithm="partition")
            Permutation Group with generators [(1,3)(4,5), (1,4)(3,5)]
            sage: C.permutation_automorphism_group(algorithm="gap")  # optional - gap_packages (Guava package)
            sage: C = codes.TernaryGolayCode()
            sage: C.permutation_automorphism_group(algorithm="gap")  # optional - gap_packages (Guava package)
            Permutation Group with generators [(3,4)(5,7)(6,9)(8,11), (3,5,8)(4,11,7)(6,9,10), (2,3)(4,6)(5,8)(7,10), (1,2)(4,11)(5,8)(9,10)]

        However, the option ``algorithm="gap+verbose"``, will print out::

            Minimum distance: 5 Weight distribution: [1, 0, 0, 0, 0, 132, 132,
            0, 330, 110, 0, 24]

            Using the 132 codewords of weight 5 Supergroup size: 39916800

        in addition to the output of
        ``C.permutation_automorphism_group(algorithm="gap")``.
        """
        F = self.base_ring()
        q = F.order()
        G = self.gen_mat() if 2*self.dimension() <= self.length() else self.dual_code().gen_mat()
        n = len(G.columns())
        k = len(G.rows())
        if "gap" in algorithm:
            gap.LoadPackage('"guava"')
            wts = self.spectrum()                                            # bottleneck 1
            nonzerowts = [i for i in range(len(wts)) if wts[i]!=0]
            Sn = SymmetricGroup(n)
            Gp = gap("SymmetricGroup(%s)"%n)               # initializing G in gap
            Gstr = str(gap(G))
            gap.eval("C:=GeneratorMatCode("+Gstr+",GF("+str(q)+"))")
            gap.eval("eltsC:=Elements(C)")
            if algorithm=="gap+verbose":
                print "\n Minimum distance: %s \n Weight distribution: \n %s"%(nonzerowts[1],wts)
            stop = 0                                          # only stop if all gens are autos
            for i in range(1,len(nonzerowts)):
                if stop == 1:
                    break
                wt = nonzerowts[i]
                if algorithm=="gap+verbose":
                    size = Gp.Size()
                    print "\n Using the %s codewords of weight %s \n Supergroup size: \n %s\n "%(wts[wt],wt,size)
                gap.eval("Cwt:=Filtered(eltsC,c->WeightCodeword(c)=%s)"%wt)   # bottleneck 2 (repeated
                gap.eval("matCwt:=List(Cwt,c->VectorCodeword(c))")            #        for each i until stop = 1)
                A = gap("MatrixAutomorphisms(matCwt)")
                #print "A = ",A, "\n Gp = ", Gp, "\n strGp = ", str(Gp)
                G2 = gap("Intersection2(%s,%s)"%(str(A).replace("\n",""),str(Gp).replace("\n",""))) #  bottleneck 3
                Gp = G2
                if Gp.Size()==1:
                    return PermutationGroup([()])
                autgp_gens = Gp.GeneratorsOfGroup()
                gens = [Sn(str(x).replace("\n","")) for x in autgp_gens]
                stop = 1                         # get ready to stop
                for x in gens:                   # if one of these gens is not an auto then don't stop
                    if not(self.is_permutation_automorphism(x)):
                        stop = 0
                        break
            G = PermutationGroup(gens)
            return G
        if algorithm=="partition":
            if q == 2:
                from sage.groups.perm_gps.partn_ref.refinement_binary import LinearBinaryCodeStruct
                B = LinearBinaryCodeStruct(G)
                autgp = B.automorphism_group()
                L = [[j+1 for j in gen] for gen in autgp[0]]
                AutGp = PermutationGroup(L)
            else:
                from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct
                from sage.matrix.constructor import matrix
                weights = {}
                for c in self:
                    wt = c.hamming_weight()
                    if wt not in weights:
                        weights[wt] = [c]
                    else:
                        weights[wt].append(c)
                weights.pop(0)
                AutGps = []
                for wt, words in weights.iteritems():
                    M = MatrixStruct(matrix(words))
                    autgp = M.automorphism_group()
                    L = [[j+1 for j in gen] for gen in autgp[0]]
                    G = PermutationGroup(L)
                    AutGps.append(G)
                if len(AutGps) > 0:
                    AutGp = AutGps[0]
                    for G in AutGps[1:]:
                        AutGp = AutGp.intersection(G)
                else:
                    return PermutationGroup([])
            return AutGp
        raise NotImplementedError("The only algorithms implemented currently are 'gap', 'gap+verbose', and 'partition'.")

    def permuted_code(self, p):
        r"""
        Returns the permuted code, which is equivalent to ``self`` via the
        column permutation ``p``.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: G = C.permutation_automorphism_group(); G
            Permutation Group with generators [(4,5)(6,7), (4,6)(5,7), (2,3)(6,7), (2,4)(3,5), (1,2)(5,6)]
            sage: g = G("(2,3)(6,7)")
            sage: Cg = C.permuted_code(g)
            sage: Cg
            Linear code of length 7, dimension 4 over Finite Field of size 2
            sage: C == Cg
            True
        """
        F = self.base_ring()
        G = self.gen_mat()
        n = len(G.columns())
        MS = MatrixSpace(F,n,n)
        Gp = G*MS(p.matrix().rows())
        return LinearCode(Gp)

    def punctured(self, L):
        r"""
        Returns the code punctured at the positions `L`,
        `L \subset \{1,2,...,n\}`. If this code `C` is of length `n` in
        GF(q) then the code `C^L` obtained from `C` by puncturing at the
        positions in `L` is the code of length `n-L` consisting of codewords
        of `C` which have their `i-th` coordinate deleted if `i \in L` and
        left alone if `i\notin L`:

        .. math::

            C^L = \{(c_{i_1},...,c_{i_N})\ |\ (c_1,...,c_n)\in C\},

        where `\{1,2,...,n\}-T = \{i_1,...,i_N\}`. In particular, if `L=\{j\}`
        then `C^L` is simply the code obtainen from `C` by deleting the `j-th`
        coordinate of each codeword. The code `C^L` is called the punctured
        code at `L`. The dimension of `C^L` can decrease if `|L|>d-1`.

        INPUT:

        - ``L`` - Subset of `\{1,...,n\}`, where `n` is the length of ``self``

        OUTPUT:

        - Linear code, the punctured code described above

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.punctured([1,2])
            Linear code of length 5, dimension 4 over Finite Field of size 2
        """
        G = self.gen_mat()
        GL = G.matrix_from_columns([i for i in range(G.ncols()) if i not in L])
        r = GL.rank()
        if r < GL.nrows():
            GL.echelonize()
            GL = GL[:r]
        return LinearCode(GL)

    def random_element(self, *args, **kwds):
        """
        Returns a random codeword; passes other positional and keyword
        arguments to ``random_element()`` method of vector space.

        OUTPUT:

        - Random element of the vector space of this code

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(4,'a'))
            sage: C.random_element() # random test
            (1, 0, 0, a + 1, 1, a, a, a + 1, a + 1, 1, 1, 0, a + 1, a, 0, a, a, 0, a, a, 1)

        Passes extra positional or keyword arguments through::

            sage: C.random_element(prob=.5, distribution='1/n') # random test
            (1, 0, a, 0, 0, 0, 0, a + 1, 0, 0, 0, 0, 0, 0, 0, 0, a + 1, a + 1, 1, 0, 0)
        """
        V = self.ambient_space()
        S = V.subspace(self.basis())
        return S.random_element(*args, **kwds)

    def redundancy_matrix(C):
        """
        If C is a linear [n,k,d] code then this function returns a
        `k \times (n-k)` matrix A such that G = (I,A) generates a code (in
        standard form) equivalent to C. If C is already in standard form and
        G = (I,A) is its generator matrix then this function simply returns
        that A.

        OUTPUT:

        - Matrix, the redundancy matrix

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.gen_mat()
             [1 0 0 0 0 1 1]
             [0 1 0 0 1 0 1]
             [0 0 1 0 1 1 0]
             [0 0 0 1 1 1 1]
            sage: C.redundancy_matrix()
             [0 1 1]
             [1 0 1]
             [1 1 0]
             [1 1 1]
            sage: C.standard_form()[0].gen_mat()
             [1 0 0 0 0 1 1]
             [0 1 0 0 1 0 1]
             [0 0 1 0 1 1 0]
             [0 0 0 1 1 1 1]
            sage: C = codes.HammingCode(2,GF(3))
            sage: C.gen_mat()
            [1 0 1 1]
            [0 1 1 2]
            sage: C.redundancy_matrix()
            [1 1]
            [1 2]
        """
        n = C.length()
        k = C.dimension()
        C1 = C.standard_form()[0]
        G1 = C1.gen_mat()
        return G1.matrix_from_columns(range(k,n))

    def sd_duursma_data(C, i):
        r"""
        Returns the Duursma data `v` and `m` of this formally s.d. code `C`
        and the type number `i` in (1,2,3,4).  Does *not* check if this code
        is actually sd.

        INPUT:

        - ``i`` - Type number

        OUTPUT:

        - Pair ``(v, m)`` as in Duursma [D]_

        REFERENCES:

        .. [D] I. Duursma, "Extremal weight enumerators and ultraspherical
           polynomials"

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),2,4)
            sage: G = MS([1,1,0,0,0,0,1,1])
            sage: C = LinearCode(G)
            sage: C == C.dual_code()  # checks that C is self dual
            True
            sage: for i in [1,2,3,4]: print C.sd_duursma_data(i)
            [2, -1]
            [2, -3]
            [2, -2]
            [2, -1]
        """
        n = C.length()
        d = C.minimum_distance()
        if i == 1:
            v = (n-4*d)/2 + 4
            m = d-3
        if i == 2:
            v = (n-6*d)/8 + 3
            m = d-5
        if i == 3:
            v = (n-4*d)/4 + 3
            m = d-4
        if i == 4:
            v = (n-3*d)/2 + 3
            m = d-3
        return [v,m]

    def sd_duursma_q(C,i,d0):
        r"""
        INPUT:

        -  ``C`` - sd code; does *not* check if `C` is actually an sd code
        -  ``i`` - Type number, one of 1,2,3,4
        -  ``d0`` - Divisor, the smallest integer such that each `A_i > 0` iff
           `i` is divisible by `d0`

        OUTPUT:

        - Coefficients `q_0, q_1, ...` of `q(T)` as in Duursma [D]_

        REFERENCES:

        - [D] - I. Duursma, "Extremal weight enumerators and ultraspherical
          polynomials"

        EXAMPLES::

            sage: C1 = codes.HammingCode(3,GF(2))
            sage: C2 = C1.extended_code(); C2
            Linear code of length 8, dimension 4 over Finite Field of size 2
            sage: C2.is_self_dual()
            True
            sage: C2.sd_duursma_q(1,1)
            2/5*T^2 + 2/5*T + 1/5
            sage: C2.sd_duursma_q(3,1)
            3/5*T^4 + 1/5*T^3 + 1/15*T^2 + 1/15*T + 1/15
        """
        q = (C.base_ring()).order()
        n = C.length()
        d = C.minimum_distance()
        d0 = C.divisor()
        if i==1 or i==2:
            if d>d0:
                c0 = QQ((n-d)*rising_factorial(d-d0,d0+1)*C.spectrum()[d])/rising_factorial(n-d0-1,d0+2)
            else:
                c0 = QQ((n-d)*C.spectrum()[d])/rising_factorial(n-d0-1,d0+2)
        if i==3 or i==4:
            if d>d0:
                c0 = rising_factorial(d-d0,d0+1)*C.spectrum()[d]/((q-1)*rising_factorial(n-d0,d0+1))
            else:
                c0 = C.spectrum()[d]/((q-1)*rising_factorial(n-d0,d0+1))
        v = ZZ(C.sd_duursma_data(i)[0])
        m = ZZ(C.sd_duursma_data(i)[1])
        if m<0 or v<0:
            raise ValueError("This case not implemented.")
        PR = PolynomialRing(QQ,"T")
        T = PR.gen()
        if i == 1:
            coefs = PR(c0*(1+3*T+2*T**2)**m*(2*T**2+2*T+1)**v).list()
            qc = [coefs[j]/binomial(4*m+2*v,m+j) for j in range(2*m+2*v+1)]
            q = PR(qc)
        if i == 2:
            F = ((T+1)**8+14*T**4*(T+1)**4+T**8)**v
            coefs = (c0*(1+T)**m*(1+4*T+6*T**2+4*T**3)**m*F).coeffs()
            qc = [coefs[j]/binomial(6*m+8*v,m+j) for j in range(4*m+8*v+1)]
            q = PR(qc)
        if i == 3:
            F = (3*T**2+4*T+1)**v*(1+3*T**2)**v
            # Note that: (3*T**2+4*T+1)(1+3*T**2)=(T+1)**4+8*T**3*(T+1)
            coefs = (c0*(1+3*T+3*T**2)**m*F).coeffs()
            qc = [coefs[j]/binomial(4*m+4*v,m+j) for j in range(2*m+4*v+1)]
            q = PR(qc)
        if i == 4:
            coefs = (c0*(1+2*T)**m*(4*T**2+2*T+1)**v).coeffs()
            qc = [coefs[j]/binomial(3*m+2*v,m+j) for j in range(m+2*v+1)]
            q = PR(qc)
        return q/q(1)

    def sd_zeta_polynomial(C, typ=1):
        r"""
        Returns the Duursma zeta function of a self-dual code using the
        construction in [D]_.

        INPUT:

        -  ``typ`` - Integer, type of this s.d. code; one of 1,2,3, or
           4 (default: 1)

        OUTPUT:

        -  Polynomial

        EXAMPLES::

            sage: C1 = codes.HammingCode(3,GF(2))
            sage: C2 = C1.extended_code(); C2
            Linear code of length 8, dimension 4 over Finite Field of size 2
            sage: C2.is_self_dual()
            True
            sage: C2.sd_zeta_polynomial()
            2/5*T^2 + 2/5*T + 1/5
            sage: C2.zeta_polynomial()
            2/5*T^2 + 2/5*T + 1/5
            sage: P = C2.sd_zeta_polynomial(); P(1)
            1
            sage: F.<z> = GF(4,"z")
            sage: MS = MatrixSpace(F, 3, 6)
            sage: G = MS([[1,0,0,1,z,z],[0,1,0,z,1,z],[0,0,1,z,z,1]])
            sage: C = LinearCode(G)  # the "hexacode"
            sage: C.sd_zeta_polynomial(4)
            1

        It is a general fact about Duursma zeta polynomials that `P(1) = 1`.

        REFERENCES:

        - [D] I. Duursma, "Extremal weight enumerators and ultraspherical
          polynomials"
        """
        d0 = C.divisor()
        P = C.sd_duursma_q(typ,d0)
        PR = P.parent()
        T = FractionField(PR).gen()
        if typ == 1:
            P0 = P
        if typ == 2:
            P0 = P/(1-2*T+2*T**2)
        if typ == 3:
            P0 = P/(1+3*T**2)
        if typ == 4:
            P0 = P/(1+2*T)
        return P0/P0(1)

    def shortened(self, L):
        r"""
        Returns the code shortened at the positions ``L``, where
        `L \subset \{1,2,...,n\}`.

        Consider the subcode `C(L)` consisting of all codewords `c\in C` which
        satisfy `c_i=0` for all `i\in L`. The punctured code `C(L)^L` is
        called the shortened code on `L` and is denoted `C_L`. The code
        constructed is actually only isomorphic to the shortened code defined
        in this way.

        By Theorem 1.5.7 in [HP], `C_L` is `((C^\perp)^L)^\perp`. This is used
        in the construction below.

        INPUT:

        - ``L`` - Subset of `\{1,...,n\}`, where `n` is the length of this code

        OUTPUT:

        - Linear code, the shortened code described above

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.shortened([1,2])
            Linear code of length 5, dimension 2 over Finite Field of size 2
        """
        Cd = self.dual_code()
        Cdp = Cd.punctured(L)
        return Cdp.dual_code()

    @rename_keyword(deprecation=11033, method="algorithm")
    def spectrum(self, algorithm=None):
        r"""
        Returns the spectrum of ``self`` as a list.

        The default algorithm uses a GAP kernel function (in C) written by
        Steve Linton.

        INPUT:

        - ``algorithm`` - ``None``, ``"gap"``, ``"leon"``, or ``"binary"``;
          defaults to ``"gap"`` except in the binary case.  If ``"gap"`` then
          uses the GAP function, if ``"leon"`` then uses Jeffrey Leon's
          software via Guava, and if ``"binary"`` then uses Sage native Cython
          code

        - List, the spectrum

        The optional algorithm (``"leon"``) may create a stack smashing error
        and a traceback but should return the correct answer. It appears to run
        much faster than the GAP algorithm in some small examples and much
        slower than the GAP algorithm in other larger examples.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G = MS([[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.spectrum()
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: F.<z> = GF(2^2,"z")
            sage: C = codes.HammingCode(2, F); C
            Linear code of length 5, dimension 3 over Finite Field in z of size 2^2
            sage: C.spectrum()
            [1, 0, 0, 30, 15, 18]
            sage: C = codes.HammingCode(3,GF(2)); C
            Linear code of length 7, dimension 4 over Finite Field of size 2
            sage: C.spectrum(algorithm="leon")   # optional - gap_packages (Guava package)
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C.spectrum(algorithm="gap")
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C.spectrum(algorithm="binary")
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C = codes.HammingCode(3,GF(3)); C
            Linear code of length 13, dimension 10 over Finite Field of size 3
            sage: C.spectrum() == C.spectrum(algorithm="leon")   # optional - gap_packages (Guava package)
            True
            sage: C = codes.HammingCode(2,GF(5)); C
            Linear code of length 6, dimension 4 over Finite Field of size 5
            sage: C.spectrum() == C.spectrum(algorithm="leon")   # optional - gap_packages (Guava package)
            True
            sage: C = codes.HammingCode(2,GF(7)); C
            Linear code of length 8, dimension 6 over Finite Field of size 7
            sage: C.spectrum() == C.spectrum(algorithm="leon")   # optional - gap_packages (Guava package)
            True

        """
        if algorithm is None:
            if self.base_ring().order() == 2:
                algorithm = "binary"
            else:
                algorithm = "gap"
        n = self.length()
        F = self.base_ring()
        G = self.gen_mat()
        if algorithm=="gap":
            Gstr = G._gap_init_()
            spec = wtdist_gap(Gstr,n,F)
            return spec
        elif algorithm=="binary":
            from sage.coding.binary_code import weight_dist
            return weight_dist(self.gen_mat())
        elif algorithm=="leon":
            if not(F.order() in [2,3,5,7]):
                raise NotImplementedError("The algorithm 'leon' is only implemented for q = 2,3,5,7.")
            # The GAP command DirectoriesPackageLibrary tells the location of the latest
            # version of the Guava libraries, so gives us the location of the Guava binaries too.
            guava_bin_dir = gap.eval('DirectoriesPackagePrograms("guava")[1]')
            guava_bin_dir = guava_bin_dir[guava_bin_dir.index('"') + 1:guava_bin_dir.rindex('"')]
            input = code2leon(self)
            from sage.misc.misc import tmp_filename
            output = tmp_filename()
            import os
            status = os.system(os.path.join(guava_bin_dir, 'wtdist')
                               + ' ' + input + "::code > " + output)
            if status != 0:
                raise RuntimeError("Problem calling Leon's wtdist program. Install gap_packages*.spkg and run './configure ../..; make'.")
            f = open(output); lines = f.readlines(); f.close()
            wts = [0]*(n+1)
            s = 0
            for L in lines:
                L = L.strip()
                if len(L) > 0:
                    o = ord(L[0])
                    if o >= 48 and o <= 57:
                        wt, num = L.split()
                        wts[eval(wt)] = eval(num)
            return wts
        else:
            raise NotImplementedError("The only algorithms implemented currently are 'gap', 'leon' and 'binary'.")

    def standard_form(self):
        r"""
        Returns the standard form of this linear code.

        An `[n,k]` linear code with generator matrix `G` in standard form is
        the row-reduced echelon form of `G` is `(I,A)`, where `I` denotes the
        `k \times k` identity matrix and `A` is a `k \times (n-k)` block. This
        method returns a pair `(C,p)` where `C` is a code permutation
        equivalent to ``self`` and `p` in `S_n`, with `n` the length of `C`,
        is the permutation sending ``self`` to `C`. This does not call GAP.

        Thanks to Frank Luebeck for (the GAP version of) this code.

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.gen_mat()
            [1 0 0 0 0 1 1]
            [0 1 0 0 1 0 1]
            [0 0 1 0 1 1 0]
            [0 0 0 1 1 1 1]
            sage: Cs,p = C.standard_form()
            sage: p
            ()
            sage: MS = MatrixSpace(GF(3),3,7)
            sage: G = MS([[1,0,0,0,1,1,0],[0,1,0,1,0,1,0],[0,0,0,0,0,0,1]])
            sage: C = LinearCode(G)
            sage: Cs, p = C.standard_form()
            sage: p
            (3,7)
            sage: Cs.gen_mat()
             [1 0 0 0 1 1 0]
             [0 1 0 1 0 1 0]
             [0 0 1 0 0 0 0]
        """
        from sage.coding.code_constructions import permutation_action as perm_action
        mat = self.gen_mat()
        MS = mat.parent()
        A = []
        k = len(mat.rows())
        M = mat.echelon_form()
        d = len(mat.columns())
        G = SymmetricGroup(d)
        perm = G([()])
        for i in range(1,k+1):
            r = M.rows()[i-1]
            j = r.nonzero_positions()[0]
            if (j < d and i <> j+1):
                perm = perm *G([(i,j+1)])
        if perm <> G([()]):
            for i in range(k):
                r = M.rows()[i]
                A.append(perm_action(perm,r))
        if perm == G([()]):
            A = M
        return LinearCode(MS(A)), perm

    def support(self):
        r"""
        Returns the set of indices `j` where `A_j` is nonzero, where
        spectrum(self) = `[A_0,A_1,...,A_n]`.

        OUTPUT:

        - List of integers

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.spectrum()
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C.support()
            [0, 3, 4, 7]
        """
        n = self.length()
        F = self.base_ring()
        V = VectorSpace(F,n+1)
        return V(self.spectrum()).support()

    def weight_enumerator(self, names="xy", name2=None):
        """
        Returns the weight enumerator of the code.

        INPUT:

        - ``names`` - String of length 2, containing two variable names
          (default: ``"xy"``). Alternatively, it can be a variable name or
          a string, or a tuple of variable names or strings.

        - ``name2`` - string or symbolic variable (default: ``None``).
          If ``name2`` is provided then it is assumed that ``names``
          contains only one variable.

        OUTPUT:

        - Polynomial over `\QQ`

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.weight_enumerator()
            x^7 + 7*x^4*y^3 + 7*x^3*y^4 + y^7
            sage: C.weight_enumerator(names="st")
            s^7 + 7*s^4*t^3 + 7*s^3*t^4 + t^7
            sage: (var1, var2) = var('var1, var2')
            sage: C.weight_enumerator((var1, var2))
            var1^7 + 7*var1^4*var2^3 + 7*var1^3*var2^4 + var2^7
            sage: C.weight_enumerator(var1, var2)
            var1^7 + 7*var1^4*var2^3 + 7*var1^3*var2^4 + var2^7

        """
        if name2 is not None:
            # We assume that actual variable names or strings are provided
            # for names if names2 is also provided. That is, names is not
            # a tuple or a list. Otherwise, PolynomialRing will return error
            names = (names, name2)
        spec = self.spectrum()
        n = self.length()
        R = PolynomialRing(QQ,2,names)
        x,y = R.gens()
        we = sum([spec[i]*x**(n-i)*y**i for i in range(n+1)])
        return we

    def zeta_polynomial(self, name="T"):
        r"""
        Returns the Duursma zeta polynomial of this code.

        Assumes that the minimum distances of this code and its dual are
        greater than 1.  Prints a warning to ``stdout`` otherwise.

        INPUT:

        - ``name`` - String, variable name (default: ``"T"``)

        OUTPUT:

        - Polynomial over `\QQ`

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.zeta_polynomial()
            2/5*T^2 + 2/5*T + 1/5
            sage: C = best_known_linear_code(6,3,GF(2))  # optional - gap_packages (Guava package)
            sage: C.minimum_distance()                   # optional - gap_packages (Guava package)
            3
            sage: C.zeta_polynomial()                    # optional - gap_packages (Guava package)
            2/5*T^2 + 2/5*T + 1/5
            sage: C = codes.HammingCode(4,GF(2))
            sage: C.zeta_polynomial()
            16/429*T^6 + 16/143*T^5 + 80/429*T^4 + 32/143*T^3 + 30/143*T^2 + 2/13*T + 1/13
            sage: F.<z> = GF(4,"z")
            sage: MS = MatrixSpace(F, 3, 6)
            sage: G = MS([[1,0,0,1,z,z],[0,1,0,z,1,z],[0,0,1,z,z,1]])
            sage: C = LinearCode(G)  # the "hexacode"
            sage: C.zeta_polynomial()
            1

        REFERENCES:

        - I. Duursma, "From weight enumerators to zeta functions", in
          Discrete Applied Mathematics, vol. 111, no. 1-2, pp. 55-73, 2001.
        """
        n = self.length()
        q = (self.base_ring()).order()
        d = self.minimum_distance()
        dperp = (self.dual_code()).minimum_distance()
        if d == 1 or dperp == 1:
            print "\n WARNING: There is no guarantee this function works when the minimum distance"
            print "            of the code or of the dual code equals 1.\n"
        RT = PolynomialRing(QQ,"%s"%name)
        R = PolynomialRing(QQ,3,"xy%s"%name)
        x,y,T = R.gens()
        we = self.weight_enumerator()
        A = R(we)
        B = A(x+y,y,T)-(x+y)**n
        Bs = B.coefficients()
        b = [Bs[i]/binomial(n,i+d) for i in range(len(Bs))]
        r = n-d-dperp+2
        #print B,Bs,b,r
        P_coeffs = []
        for i in range(len(b)):
           if i == 0:
              P_coeffs.append(b[0])
           if i == 1:
              P_coeffs.append(b[1] - (q+1)*b[0])
           if i>1:
              P_coeffs.append(b[i] - (q+1)*b[i-1] + q*b[i-2])
        #print P_coeffs
        P = sum([P_coeffs[i]*T**i for i in range(r+1)])
        return RT(P)/RT(P)(1)

    def zeta_function(self, name="T"):
        r"""
        Returns the Duursma zeta function of the code.

        INPUT:

        - ``name`` - String, variable name (default: ``"T"``)

        OUTPUT:

        - Element of `\QQ(T)`

        EXAMPLES::

            sage: C = codes.HammingCode(3,GF(2))
            sage: C.zeta_function()
            (2/5*T^2 + 2/5*T + 1/5)/(2*T^2 - 3*T + 1)
        """
        P =  self.zeta_polynomial()
        q = (self.base_ring()).characteristic()
        RT = PolynomialRing(QQ,"%s"%name)
        T = RT.gen()
        return P/((1-T)*(1-q*T))

    weight_distribution = spectrum

def LinearCodeFromVectorSpace(V, d=None):
    """
    Simply converts a vector subspace `V` of `GF(q)^n` into a `LinearCode`.

    INPUT:

    - ``V`` -- The vector space

    - ``d`` -- (Optional, default: ``None``) the minimum distance of the
      code, if known. This is an optional parameter.

    .. note::
        The veracity of the minimum distance ``d``, if provided, is not
        checked.

    EXAMPLES::

        sage: V = VectorSpace(GF(2), 8)
        sage: L = V.subspace([[1,1,1,1,0,0,0,0],[0,0,0,0,1,1,1,1]])
        sage: C = LinearCodeFromVectorSpace(L)
        sage: C.gen_mat()
        [1 1 1 1 0 0 0 0]
        [0 0 0 0 1 1 1 1]
        sage: C.minimum_distance()
        4

    Here, we provide the minimum distance of the code.::

        sage: C = LinearCodeFromVectorSpace(L, d=4)
        sage: C.minimum_distance()
        4
    """
    F = V.base_ring()
    B = V.basis()
    n = len(B[0].list())
    k = len(B)
    MS = MatrixSpace(F,k,n)
    G = MS([B[i].list() for i in range(k)])
    return LinearCode(G, d=d)

