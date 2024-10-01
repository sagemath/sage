# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Small Scale Variants of the AES (SR) Polynomial System Generator

Sage supports polynomial system generation for small scale (and full
scale) AES variants over `\GF{2}` and `\GF{2^e}`. Also, Sage supports
both the specification of SR as given in the papers [CMR2005]_ and
[CMR2006]_ and a variant of SR* which is equivalent to AES.

SR is a family of parameterizable variants of the AES suitable as a
framework for comparing different cryptanalytic techniques that can be
brought to bear on the AES. It is different from
:class:`Mini-AES <sage.crypto.block_cipher.miniaes.MiniAES>`, whose
purpose is as a teaching tool to help beginners understand the basic
structure and working of the full AES.

AUTHORS:

- Martin Albrecht (2008,2009-01): usability improvements

- Martin Albrecht (2007-09): initial version

- Niles Johnson (2010-08): (:issue:`3893`) ``random_element()`` should pass on ``*args`` and ``**kwds``.

EXAMPLES:

We construct SR(1,1,1,4) and study its properties.
::

    sage: sr = mq.SR(1, 1, 1, 4)

``n`` is the number of rounds, ``r`` the number of rows in the
state array, ``c`` the number of columns in the state array, and ``e`` the
degree of the underlying field.

::

    sage: sr.n, sr.r, sr.c, sr.e
    (1, 1, 1, 4)

By default variables are ordered reverse to as they appear, e.g.::

    sage: print(sr.R.repr_long())
    Polynomial Ring
      Base Ring : Finite Field in a of size 2^4
           Size : 20 Variables
       Block  0 : Ordering : deglex
                  Names    : k100, k101, k102, k103, x100, x101, x102, x103, w100, w101, w102, w103, s000, s001, s002, s003, k000, k001, k002, k003

However, this can be prevented by passing in ``reverse_variables=False`` to the constructor.

For SR(1, 1, 1, 4) the ``ShiftRows`` matrix isn't that interesting.::

    sage: sr.ShiftRows
    [1 0 0 0]
    [0 1 0 0]
    [0 0 1 0]
    [0 0 0 1]

Also, the ``MixColumns`` matrix is the identity matrix.::

    sage: sr.MixColumns
    [1 0 0 0]
    [0 1 0 0]
    [0 0 1 0]
    [0 0 0 1]

``Lin``, however, is not the identity matrix.::

    sage: sr.Lin
    [          a^2 + 1                 1         a^3 + a^2           a^2 + 1]
    [                a                 a                 1 a^3 + a^2 + a + 1]
    [          a^3 + a               a^2               a^2                 1]
    [                1               a^3             a + 1             a + 1]

``M`` and ``Mstar`` are identical for SR(1, 1, 1, 4)::

    sage: sr.M
    [          a^2 + 1                 1         a^3 + a^2           a^2 + 1]
    [                a                 a                 1 a^3 + a^2 + a + 1]
    [          a^3 + a               a^2               a^2                 1]
    [                1               a^3             a + 1             a + 1]

::

    sage: sr.Mstar
    [          a^2 + 1                 1         a^3 + a^2           a^2 + 1]
    [                a                 a                 1 a^3 + a^2 + a + 1]
    [          a^3 + a               a^2               a^2                 1]
    [                1               a^3             a + 1             a + 1]


However, for larger instances of SR Mstar is not equal to M::

    sage: sr = mq.SR(10,4,4,8)
    sage: sr.Mstar == ~sr.MixColumns * sr.M
    True

We can compute a Groebner basis for the ideals spanned by SR
instances to recover all solutions to the system.::

    sage: sr = mq.SR(1,1,1,4, gf2=True, polybori=True)
    sage: K = sr.base_ring()
    sage: a = K.gen()
    sage: K = [a]
    sage: P = [1]
    sage: F,s = sr.polynomial_system(P=P, K=K)                                          # needs sage.rings.polynomial.pbori
    sage: F.groebner_basis()                                                            # needs sage.rings.polynomial.pbori
    [k100, k101 + 1, k102, k103 + k003,
     x100 + 1, x101 + k003 + 1, x102 + k003 + 1,
     x103 + k003, w100, w101, w102 + 1, w103 + k003 + 1,
     s000 + 1, s001 + k003, s002 + k003, s003 + k003 + 1,
     k000, k001, k002 + 1]

Note that the order of ``k000``, ``k001``, ``k002`` and ``k003`` is
little endian. Thus the result ``k002 + 1, k001, k000`` indicates that
the key is either `a` or `a+1`. We can verify that both keys encrypt P
to the same ciphertext::

    sage: sr(P,[a])
    [0]
    sage: sr(P,[a+1])
    [0]

All solutions can easily be recovered using the variety function for ideals.::

   sage: I = F.ideal()                                                                  # needs sage.rings.polynomial.pbori
   sage: for V in I.variety():                                                          # needs sage.rings.polynomial.pbori sage.symbolic
   ....:    for k,v in sorted(V.items()):
   ....:       print("{} {}".format(k, v))
   ....:    print("\n")
   k003 0
   k002 1
   k001 0
   k000 0
   s003 1
   s002 0
   s001 0
   s000 1
   w103 1
   w102 1
   w101 0
   w100 0
   x103 0
   x102 1
   x101 1
   x100 1
   k103 0
   k102 0
   k101 1
   k100 0
   <BLANKLINE>
   k003 1
   k002 1
   k001 0
   k000 0
   s003 0
   s002 1
   s001 1
   s000 1
   w103 0
   w102 1
   w101 0
   w100 0
   x103 1
   x102 0
   x101 0
   x100 1
   k103 1
   k102 0
   k101 1
   k100 0

We can also verify the correctness of the variety by evaluating all
ideal generators on all points.::

   sage: for V in I.variety():                                                          # needs sage.rings.polynomial.pbori sage.symbolic
   ....:     for f in I.gens():
   ....:         if f.subs(V) != 0:
   ....:            print("epic fail")


Note that the S-Box object for SR can be constructed with a call to ``sr.sbox()``::

   sage: sr = mq.SR(1,1,1,4, gf2=True, polybori=True)
   sage: S = sr.sbox()

For example, we can now study the difference distribution table of ``S``::

   sage: S.difference_distribution_table()
   [16  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
   [ 0  2  2  2  2  0  0  0  2  0  0  0  2  4  0  0]
   [ 0  2  0  4  2  2  2  0  0  2  0  0  0  0  0  2]
   [ 0  2  4  0  0  2  0  0  2  2  0  2  0  0  2  0]
   [ 0  0  2  0  4  2  0  0  0  0  2  0  2  0  2  2]
   [ 0  0  0  2  0  0  0  2  4  2  0  0  2  0  2  2]
   [ 0  4  0  0  0  2  0  2  0  2  2  0  2  2  0  0]
   [ 0  2  0  0  0  0  2  0  0  0  0  2  4  2  2  2]
   [ 0  2  2  0  0  0  2  2  2  0  2  0  0  0  0  4]
   [ 0  0  2  2  0  0  0  0  0  2  2  4  0  2  0  2]
   [ 0  0  2  0  2  0  2  2  0  4  0  2  2  0  0  0]
   [ 0  0  0  0  2  0  2  0  2  2  4  0  0  2  2  0]
   [ 0  0  0  2  0  4  2  0  2  0  2  2  2  0  0  0]
   [ 0  0  0  0  2  2  0  4  2  0  0  2  0  2  0  2]
   [ 0  0  2  2  0  2  4  2  0  0  0  0  0  2  2  0]
   [ 0  2  0  2  2  0  0  2  0  0  2  2  0  0  4  0]

or use ``S`` to find alternative polynomial representations for the S-Box.::

   sage: S.polynomials(degree=3)                                                        # needs sage.libs.singular
   [x0*x1 + x1*x2 + x0*x3 + x0*y2 + x1 + y0 + y1 + 1,
    x0*x1 + x0*x2 + x0*y0 + x0*y1 + x0*y2 + x1 + x2 + y0 + y1 + y2,
    x0*x1 + x0*x2 + x0*x3 + x1*x3 + x0*y0 + x1*y0 + x0*y1 + x0*y3,
    x0*x1 + x0*x2 + x0*x3 + x1*x3 + x0*y0 + x1*y1 + x0*y3 + x1 + y0 + y1 + 1,
    x0*x1 + x0*x2 + x0*y2 + x1*y2 + x0*y3 + x0 + x1,
    x0*x3 + x1*x3 + x0*y1 + x0*y2 + x1*y3 + x0 + x1 + x2 + x3 + y0 + y1 + y3 + 1,
    x0*x1 + x1*x3 + x2*x3 + x0*y0 + x0*y2 + x0*y3 + x2 + y0 + y3,
    x0*x1 + x0*x2 + x0*x3 + x1*x3 + x2*y0 + x0*y2 + x0 + x2 + x3 + y3,
    x0*x3 + x1*x3 + x0*y0 + x2*y1 + x0*y2 + x3 + y3,
    x0*x1 + x0*x2 + x0*y0 + x0*y1 + x2*y2 + x0*y3 + x1 + y0 + y1 + 1,
    x0*x3 + x1*x3 + x0*y0 + x0*y1 + x0*y3 + x2*y3 + y0 + y3,
    x0*x1 + x0*x2 + x3*y0 + x0*y1 + x0*y3 + y0,
    x0*y0 + x0*y1 + x3*y1 + x0 + x2 + y0 + y3,
    x0*y0 + x3*y2 + y0,
    x0*x1 + x0*x2 + x0*x3 + x1*x3 + x0*y0 + x0*y2 + x3*y3 + y0,
    x0*x2 + x0*x3 + x0*y1 + y0*y1 + x0*y3 + x2 + x3 + y3,
    x0*x2 + x0*y0 + y0*y2 + x0*y3 + x0 + y0,
    x0*x1 + x0*x2 + x1*x3 + x0*y2 + y0*y3 + y0,
    x0*x1 + x0*y0 + y1*y2 + x0*y3 + x1 + x2 + y0 + 1,
    x0*x2 + x1*x3 + x0*y1 + x0*y2 + x0*y3 + y1*y3 + x0 + y0 + y3,
    x0*x1 + x0*x2 + x0*x3 + x0*y1 + x0*y2 + x0*y3 + y2*y3 + x0 + x1 + x2 + x3 + y1 + y3 + 1,
    x0*x1*x2 + x0*x3 + x0*y0 + x0*y1 + x0*y2 + x0,
    x0*x1*x3 + x0*x2 + x0*x3 + x0*y1 + x0*y3 + x0,
    x0*x1*y0 + x0*x1 + x0*y0 + x0,
    x0*x1*y1,
    x0*x1*y2 + x0*x2 + x0*y2 + x0*y3 + x0,
    x0*x1*y3 + x0*x1 + x0*x3 + x0*y0 + x0*y1 + x0*y2 + x0,
    x0*x2*x3 + x0*x1 + x0*x3 + x0*y1 + x0*y2 + x0*y3 + x0,
    x0*x2*y0 + x0*x1 + x0*x2 + x0*x3 + x0*y1 + x0*y2,
    x0*x2*y1 + x0*x2 + x0*x3 + x0*y0 + x0*y1 + x0*y2 + x0,
    x0*x2*y2 + x0*x2 + x0*y3 + x0,
    x0*x2*y3 + x0*x2 + x0*y3 + x0,
    x0*x3*y0 + x0*x1 + x0*x2 + x0*y0 + x0*y1 + x0*y3,
    x0*x3*y1 + x0*x2 + x0*y1 + x0*y3 + x0,
    x0*x3*y2,
    x0*x3*y3 + x0*x1 + x0*y1 + x0*y2 + x0*y3 + x0,
    x0*y0*y1 + x0*y1,
    x0*y0*y2 + x0*x2 + x0*y3 + x0,
    x0*y0*y3 + x0*x1 + x0*x3 + x0*y0 + x0*y1 + x0*y2 + x0*y3 + x0,
    x0*y1*y2 + x0*x2 + x0*y3 + x0,
    x0*y1*y3 + x0*x3 + x0*y0 + x0*y2 + x0*y3,
    x0*y2*y3 + x0*y2,
    x1*x2*x3 + x0*x1 + x1*x3 + x0*y0 + x0*y1 + x2 + x3 + y3,
    x1*x2*y0 + x0*x1 + x1*x3 + x0*y0 + x0*y1 + x2 + x3 + y3,
    x1*x2*y1 + x0*x1 + x1*x3 + x0*y0 + x1 + x2 + x3 + y0 + y1 + y3 + 1,
    x1*x2*y2 + x0*x1 + x0*y0 + x0*y1 + x0 + x1 + y0 + y1 + 1,
    x1*x2*y3 + x0*x1 + x1*x3 + x0*y0 + x1 + x2 + x3 + y0 + y1 + y3 + 1,
    x1*x3*y0 + x0*x1 + x0*x2 + x0*x3 + x1*x3 + x0*y0 + x0*y1 + x0*y3,
    x1*x3*y1 + x0*x2 + x0*x3 + x0*y3 + x2 + x3 + y3,
    x1*x3*y2 + x0*x2 + x0*x3 + x1*x3 + x0*y1 + x0*y3 + x0,
    x1*x3*y3 + x0*x1 + x0*x2 + x0*x3 + x0*y0 + x0*y1 + x0*y3,
    x1*y0*y1 + x0*x2 + x0*x3 + x0*y3 + x2 + x3 + y3,
    x1*y0*y2 + x0*x2 + x0*x3 + x1*x3 + x0*y1 + x0*y3 + x0,
    x1*y0*y3,
    x1*y1*y2 + x0*x1 + x0*x2 + x0*x3 + x1*x3 + x0*y0 + x0*y3 + x1 + y0 + y1 + 1,
    x1*y1*y3 + x0*x1 + x1*x3 + x0*y0 + x1 + x2 + x3 + y0 + y1 + y3 + 1,
    x1*y2*y3 + x0*x1 + x0*x2 + x1*x3 + x0*y0 + x0*y2 + x0*y3 + x0 + x1 + x2 + x3 + y0 + y1 + y3 + 1,
    x2*x3*y0 + x0*x1 + x0*x3 + x1*x3 + x0*y2 + x0*y3 + x2 + x3 + y3,
    x2*x3*y1 + x0*y1 + x0*y2 + x0*y3 + x3 + y0,
    x2*x3*y2 + x1*x3 + x0*y1 + x0 + x2 + x3 + y3,
    x2*x3*y3,
    x2*y0*y1 + x0*x2 + x0*x3 + x0*y0 + x0*y1 + x0*y2 + x0,
    x2*y0*y2 + x0*x2 + x1*x3 + x0*y1 + x0*y3 + x2 + x3 + y3,
    x2*y0*y3 + x0*x2 + x0*y3 + x0,
    x2*y1*y2 + x0*x1 + x0*x2 + x1*x3 + x0*y0 + x0*y3 + x0 + x1 + x2 + x3 + y0 + y1 + y3 + 1,
    x2*y1*y3 + x0*x3 + x1*x3 + x0*y0 + x0*y1 + x0*y3 + y0 + y3,
    x2*y2*y3 + x0*x1 + x0*x2 + x1*x3 + x0*y0 + x0*y3 + x0 + x1 + x2 + x3 + y0 + y1 + y3 + 1,
    x3*y0*y1 + x0*x3 + x0*y1 + x0 + x2 + x3 + y3,
    x3*y0*y2 + x0*y0 + y0,
    x3*y0*y3 + x1*x3 + x0*y1 + x0*y2 + x0*y3 + y0,
    x3*y1*y2 + x0*x2 + x0*x3 + x0*y3 + x2 + x3 + y3,
    x3*y1*y3 + x0*x2 + x0*x3 + x0*y0 + x0*y2 + x0,
    x3*y2*y3 + x0*x2 + x0*x3 + x1*x3 + x0*y0 + x0*y1 + x0*y3 + x0 + y0,
    y0*y1*y2 + x0*x3 + x0 + x2 + x3 + y3,
    y0*y1*y3 + x0*x3 + x0*y0 + x0*y2 + x0*y3,
    y0*y2*y3 + x0*x3 + x1*x3 + x0*y0 + x0*y1 + y0,
    y1*y2*y3 + x0*x1 + x0*x2 + x1*x3 + x0*y0 + x0*y3 + x0 + x1 + x2 + x3 + y0 + y1 + y3 + 1]

   sage: S.interpolation_polynomial()
   (a^2 + 1)*x^14 + a^2*x^13 + x^12 + a^2*x^11 + a*x^10 + (a^3 + a)*x^9 +
   (a^3 + 1)*x^7 + (a^3 + a^2 + a)*x^6 + a^2*x^5 + (a + 1)*x^4 + a^2*x^3 +
   (a^3 + a^2 + a)*x^2 + (a^3 + 1)*x + a^2 + a

The :class:`SR_gf2_2` gives an example how use alternative polynomial
representations of the S-Box for construction of polynomial systems.

TESTS::

    sage: sr == loads(dumps(sr))
    True

REFERENCES:

- [CMR2005]_

- [CMR2006]_

- [MR2002]_
"""

from sage.matrix.constructor import matrix, random_matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.flatten import flatten
from sage.misc.verbose import get_verbose
from sage.modules.vector_modn_dense import Vector_modn_dense
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
from sage.rings.polynomial.polynomial_ring_constructor import \
    BooleanPolynomialRing_constructor as BooleanPolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder
from sage.structure.element import Matrix

from .mpolynomialsystemgenerator import MPolynomialSystemGenerator


def SR(n=1, r=1, c=1, e=4, star=False, **kwargs):
    r"""
    Return a small scale variant of the AES polynomial system
    constructor subject to the following conditions:

    INPUT:

    - ``n`` -- the number of rounds (default: 1)
    - ``r`` -- the number of rows in the state array (default: 1)
    - ``c`` -- the number of columns in the state array (default: 1)
    - ``e`` -- the exponent of the finite extension field (default: 4)
    - ``star`` -- determines if SR\* or SR should be constructed (default: ``False``)
    - ``aes_mode`` -- as the SR key schedule specification differs
      slightly from the AES key schedule, this parameter controls
      which schedule to use (default: ``True``)
    - ``gf2`` -- generate polynomial systems over `\GF{2}` rather than
      over `\GF{2^e}` (default: ``False``)
    - ``polybori`` -- use the ``BooleanPolynomialRing`` as polynomial
      representation (default: ``True``, `\GF{2}` only)
    - ``order`` -- string to specify the term ordering of the
      variables (default: ``deglex``)
    - ``postfix`` -- string which is appended after the variable name
      (default: ``''``)
    - ``allow_zero_inversions`` -- boolean to control whether zero
      inversions raise an exception (default: ``False``)
    - ``correct_only`` -- only include correct inversion polynomials
      (default: ``False``, `\GF{2}` only)
    - ``biaffine_only`` -- only include bilinear and biaffine inversion
      polynomials (default: ``True``, `\GF{2}` only)

    EXAMPLES::

        sage: sr = mq.SR(1, 1, 1, 4)
        sage: ShiftRows = sr.shift_rows_matrix()
        sage: MixColumns = sr.mix_columns_matrix()
        sage: Lin = sr.lin_matrix()
        sage: M = MixColumns * ShiftRows * Lin
        sage: print(sr.hex_str_matrix(M))
         5 1 C 5
         2 2 1 F
         A 4 4 1
         1 8 3 3

    ::

        sage: sr = mq.SR(1, 2, 1, 4)
        sage: ShiftRows = sr.shift_rows_matrix()
        sage: MixColumns = sr.mix_columns_matrix()
        sage: Lin = sr.lin_matrix()
        sage: M = MixColumns * ShiftRows * Lin
        sage: print(sr.hex_str_matrix(M))
         F 3 7 F A 2 B A
         A A 5 6 8 8 4 9
         7 8 8 2 D C C 3
         4 6 C C 5 E F F
         A 2 B A F 3 7 F
         8 8 4 9 A A 5 6
         D C C 3 7 8 8 2
         5 E F F 4 6 C C

    ::

        sage: sr = mq.SR(1, 2, 2, 4)
        sage: ShiftRows = sr.shift_rows_matrix()
        sage: MixColumns = sr.mix_columns_matrix()
        sage: Lin = sr.lin_matrix()
        sage: M = MixColumns * ShiftRows * Lin
        sage: print(sr.hex_str_matrix(M))
         F 3 7 F 0 0 0 0 0 0 0 0 A 2 B A
         A A 5 6 0 0 0 0 0 0 0 0 8 8 4 9
         7 8 8 2 0 0 0 0 0 0 0 0 D C C 3
         4 6 C C 0 0 0 0 0 0 0 0 5 E F F
         A 2 B A 0 0 0 0 0 0 0 0 F 3 7 F
         8 8 4 9 0 0 0 0 0 0 0 0 A A 5 6
         D C C 3 0 0 0 0 0 0 0 0 7 8 8 2
         5 E F F 0 0 0 0 0 0 0 0 4 6 C C
         0 0 0 0 A 2 B A F 3 7 F 0 0 0 0
         0 0 0 0 8 8 4 9 A A 5 6 0 0 0 0
         0 0 0 0 D C C 3 7 8 8 2 0 0 0 0
         0 0 0 0 5 E F F 4 6 C C 0 0 0 0
         0 0 0 0 F 3 7 F A 2 B A 0 0 0 0
         0 0 0 0 A A 5 6 8 8 4 9 0 0 0 0
         0 0 0 0 7 8 8 2 D C C 3 0 0 0 0
         0 0 0 0 4 6 C C 5 E F F 0 0 0 0
    """
    if not kwargs.get("gf2", False):
        return SR_gf2n(n, r, c, e, star, **kwargs)
    else:
        return SR_gf2(n, r, c, e, star, **kwargs)


class SR_generic(MPolynomialSystemGenerator):
    def __init__(self, n=1, r=1, c=1, e=4, star=False, **kwargs):
        """
        Small Scale Variants of the AES.

        EXAMPLES::

            sage: sr = mq.SR(1, 1, 1, 4)
            sage: ShiftRows = sr.shift_rows_matrix()
            sage: MixColumns = sr.mix_columns_matrix()
            sage: Lin = sr.lin_matrix()
            sage: M = MixColumns * ShiftRows * Lin
            sage: print(sr.hex_str_matrix(M))
             5 1 C 5
             2 2 1 F
             A 4 4 1
             1 8 3 3
        """
        if n-1 not in range(10):
            raise TypeError("n must be between 1 and 10 (inclusive)")
        self._n = n

        if r not in (1, 2, 4):
            raise TypeError("r must be in (1, 2, 4)")
        self._r = r

        if c not in (1, 2, 4):
            raise TypeError("c must be in (1, 2, 4)")
        self._c = c

        if e not in (4, 8):
            raise TypeError("e must be either 4 or 8")
        self._e = e

        self._star = bool(star)

        self._base = self.base_ring()

        self._postfix = kwargs.get("postfix", "")
        self._order = kwargs.get("order", "deglex")
        self._aes_mode = kwargs.get("aes_mode", True)
        self._gf2 = kwargs.get("gf2", False)
        self._allow_zero_inversions = bool(kwargs.get("allow_zero_inversions", False))
        self._reverse_variables = bool(kwargs.get("reverse_variables", True))

        with AllowZeroInversionsContext(self):
            sub_byte_lookup = dict([(v, self.sub_byte(v)) for v in self._base])
        self._sub_byte_lookup = sub_byte_lookup

        if self._gf2:
            self._polybori = kwargs.get("polybori", True)

    def new_generator(self, **kwds):
        r"""
        Return a new ``SR`` instance equal to this instance
        except for the parameters passed explicitly to this function.

        INPUT:

        - ``**kwds`` -- see the ``SR`` constructor for accepted
          parameters

        EXAMPLES::

            sage: sr = mq.SR(2,1,1,4); sr
            SR(2,1,1,4)
            sage: sr.ring().base_ring()
            Finite Field in a of size 2^4

        ::

            sage: sr2 = sr.new_generator(gf2=True); sr2
            SR(2,1,1,4)
            sage: sr2.ring().base_ring()
            Finite Field of size 2
            sage: sr3 = sr2.new_generator(correct_only=True)
            sage: len(sr2.inversion_polynomials_single_sbox())
            20
            sage: len(sr3.inversion_polynomials_single_sbox())
            19
        """
        kwds.setdefault("n", self._n)
        kwds.setdefault("r", self._r)
        kwds.setdefault("c", self._c)
        kwds.setdefault("e", self._e)
        kwds.setdefault("star", self._star)
        kwds.setdefault("postfix", self._postfix)
        kwds.setdefault("order", self._order)
        kwds.setdefault("allow_zero_inversions", self._allow_zero_inversions)
        kwds.setdefault("aes_mode", self._aes_mode)
        kwds.setdefault("gf2", self._gf2)
        kwds.setdefault("reverse_variables", self._reverse_variables)

        try:
            polybori = self._polybori
        except AttributeError:
            polybori = False
        kwds.setdefault("polybori", polybori)

        try:
            correct_only = self._correct_only
        except AttributeError:
            correct_only = False
        kwds.setdefault("correct_only", correct_only)

        try:
            biaffine_only = self._biaffine_only
        except AttributeError:
            biaffine_only = False
        kwds.setdefault("biaffine_only", biaffine_only)

        if self._gf2 == kwds.get('gf2'):
            return self.__class__(**kwds)
        else:
            return SR(**kwds)

    def __getattr__(self, attr):
        """
        EXAMPLES::

            sage: sr = mq.SR(1, 2, 1, 4, gf2=True)
            sage: sr.Mstar
            [1 0 1 1 0 0 0 0]
            [1 1 0 1 0 0 0 0]
            [1 1 1 0 0 0 0 0]
            [0 1 1 1 0 0 0 0]
            [0 0 0 0 1 0 1 1]
            [0 0 0 0 1 1 0 1]
            [0 0 0 0 1 1 1 0]
            [0 0 0 0 0 1 1 1]
        """
        if attr == "e":
            return self._e
        elif attr == "c":
            return self._c
        elif attr == "n":
            return self._n
        elif attr == "r":
            return self._r

        elif attr == "R":
            self.R = self.ring()
            return self.R
        elif attr == "k":
            self.k = self.base_ring()
            return self.k

        elif attr == "Lin":
            self.Lin = self.lin_matrix()
            return self.Lin

        elif attr == "ShiftRows":
            self.ShiftRows = self.shift_rows_matrix()
            return self.ShiftRows

        elif attr == "MixColumns":
            self.MixColumns = self.mix_columns_matrix()
            return self.MixColumns

        elif attr == "M":
            self.M = self.MixColumns * self.ShiftRows * self.Lin
            return self.M

        elif attr == "Mstar":
            self.Mstar = self.ShiftRows * self.Lin
            return self.Mstar

        raise AttributeError("%s has no attribute %s" % (type(self), attr))

    def _repr_(self):
        """
        EXAMPLES::

            sage: sr = mq.SR(1, 2, 2, 4); sr #indirect doctest
            SR(1,2,2,4)
            sage: sr = mq.SR(1, 2, 2, 4, star=True); sr
            SR*(1,2,2,4)
        """
        if self._star:
            return "SR*(%d,%d,%d,%d)" % (self._n, self._r, self._c, self._e)
        else:
            return "SR(%d,%d,%d,%d)" % (self._n, self._r, self._c, self._e)

    def base_ring(self):
        r"""
        Return the base field of ``self`` as determined by
        ``self.e``.

        EXAMPLES::

            sage: sr = mq.SR(10, 2, 2, 4)
            sage: sr.base_ring().polynomial()
            a^4 + a + 1

        The Rijndael polynomial::

            sage: sr = mq.SR(10, 4, 4, 8)
            sage: sr.base_ring().polynomial()
            a^8 + a^4 + a^3 + a + 1
        """
        try:
            return self._base
        except AttributeError:
            if self._e == 4:
                self._base = GF(2**4, 'a', modulus=(1, 1, 0, 0, 1))
            elif self._e == 8:
                self._base = GF(2**8, 'a', modulus=(1, 1, 0, 1, 1, 0, 0, 0, 1))

            return self._base

    def __eq__(self, other):
        """
        Two generators are considered equal if they agree on all parameters
        passed to them during construction.

        EXAMPLES::

            sage: sr1 = mq.SR(2, 2, 2, 4)
            sage: sr1 == sr1
            True

            sage: sr2 = mq.SR(2, 2, 2, 4, gf2=True)
            sage: sr1 == sr2
            False
        """
        for name in ['n', 'r', 'c', 'e', '_postfix', '_order',
                     '_allow_zero_inversions', '_aes_mode', '_gf2', '_star']:
            lx = getattr(self, name)
            rx = getattr(other, name)
            if lx != rx:
                return False
        return True

    def __ne__(self, other):
        """
        Return whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: sr1 = mq.SR(2, 2, 2, 4)
            sage: sr1 != sr1
            False

            sage: sr2 = mq.SR(2, 2, 2, 4, gf2=True)
            sage: sr1 != sr2
            True
        """
        return not (self == other)

    def sub_bytes(self, d):
        r"""
        Perform the non-linear transform on ``d``.

        INPUT:

        - ``d`` -- state array or something coercible to a state array

        EXAMPLES::

            sage: sr = mq.SR(2, 1, 2, 8, gf2=True)
            sage: k = sr.base_ring()
            sage: A = Matrix(k, 1, 2 , [k(1), k.gen()])
            sage: sr.sub_bytes(A)
            [  a^6 + a^5 + a^4 + a^3 + a^2 a^6 + a^5 + a^4 + a^2 + a + 1]
        """
        d = self.state_array(d)
        return matrix(self.base_ring(), d.nrows(), d.ncols(), [self.sub_byte(b) for b in d.list()])

    def sub_byte(self, b):
        r"""
        Perform ``SubByte`` on a single byte/halfbyte ``b``.

        A ``ZeroDivision`` exception is raised if an attempt is made
        to perform an inversion on the zero element. This can be
        disabled by passing ``allow_zero_inversion=True`` to the
        constructor. A zero inversion can result in an inconsistent
        equation system.

        INPUT:

        - ``b`` -- an element in ``self.base_ring()``

        EXAMPLES:

        The S-Box table for `\GF{2^4}`::

            sage: sr = mq.SR(1, 1, 1, 4, allow_zero_inversions=True)
            sage: for e in sr.base_ring():
            ....:    print('% 20s % 20s'%(e, sr.sub_byte(e)))
                            0              a^2 + a
                            a              a^2 + 1
                          a^2                    a
                          a^3              a^3 + 1
                        a + 1                  a^2
                      a^2 + a          a^2 + a + 1
                    a^3 + a^2                a + 1
                  a^3 + a + 1            a^3 + a^2
                      a^2 + 1        a^3 + a^2 + a
                      a^3 + a    a^3 + a^2 + a + 1
                  a^2 + a + 1              a^3 + a
                a^3 + a^2 + a                    0
            a^3 + a^2 + a + 1                  a^3
                a^3 + a^2 + 1                    1
                      a^3 + 1        a^3 + a^2 + 1
                            1          a^3 + a + 1
        """
        if not b:
            if not self._allow_zero_inversions:
                raise ZeroDivisionError("A zero inversion occurred during an encryption or key schedule.")
            else:
                return self.sbox_constant()
        try:
            return self._sub_byte_lookup[b]
        except AttributeError:
            e = self.e
            k = self.k

            # inversion
            b = b ** ( 2**e - 2 )

            # GF(2) linear map
            if e == 4:
                if not hasattr(self, "_L"):
                    self._L = matrix(GF(2), 4, 4, [[1, 1, 1, 0],
                                                [0, 1, 1, 1],
                                                [1, 0, 1, 1],
                                                [1, 1, 0, 1]])

            elif e == 8:
                if not hasattr(self, "_L"):
                    self._L = matrix(GF(2), 8, 8, [[1, 0, 0, 0, 1, 1, 1, 1],
                                                [1, 1, 0, 0, 0, 1, 1, 1],
                                                [1, 1, 1, 0, 0, 0, 1, 1],
                                                [1, 1, 1, 1, 0, 0, 0, 1],
                                                [1, 1, 1, 1, 1, 0, 0, 0],
                                                [0, 1, 1, 1, 1, 1, 0, 0],
                                                [0, 0, 1, 1, 1, 1, 1, 0],
                                                [0, 0, 0, 1, 1, 1, 1, 1]])

            b = k(self._L * b._vector_())

            # constant addition
            if e == 4:
                b = b + k.from_integer(6)
            elif e == 8:
                b = b + k.from_integer(99)

            return b

    def sbox_constant(self):
        """
        Return the S-Box constant which is added after `L(x^{-1})` was
        performed. That is ``0x63`` if ``e == 8`` or ``0x6`` if ``e ==
        4``.

        EXAMPLES::

            sage: sr = mq.SR(10, 1, 1, 8)
            sage: sr.sbox_constant()
            a^6 + a^5 + a + 1
        """
        k = self.k
        if self.e == 4:
            return k.from_integer(6)
        elif self.e == 8:
            return k.from_integer(99)
        else:
            raise TypeError("sbox constant only defined for e in (4, 8)")

    def sbox(self, inversion_only=False):
        r"""
        Return an S-Box object for this SR instance.

        INPUT:

        - ``inversion_only`` -- do not include the `\GF{2}` affine map when
          computing the S-Box (default: ``False``)

        EXAMPLES::

            sage: sr = mq.SR(1,2,2,4, allow_zero_inversions=True)
            sage: S = sr.sbox(); S
            (6, 11, 5, 4, 2, 14, 7, 10, 9, 13, 15, 12, 3, 1, 0, 8)

            sage: sr.sub_byte(0)
            a^2 + a
            sage: sage_eval(str(sr.sub_byte(0)), {'a':2})
            6
            sage: S(0)
            6

            sage: sr.sub_byte(1)
            a^3 + a + 1
            sage: sage_eval(str(sr.sub_byte(1)), {'a':2})
            11
            sage: S(1)
            11

            sage: sr = mq.SR(1,2,2,8, allow_zero_inversions=True)
            sage: S = sr.sbox(); S
            (99, 124, 119, 123, 242, 107, 111, 197, 48, 1, 103, 43,
            254, 215, 171, 118, 202, 130, 201, 125, 250, 89, 71, 240,
            173, 212, 162, 175, 156, 164, 114, 192, 183, 253, 147, 38,
            54, 63, 247, 204, 52, 165, 229, 241, 113, 216, 49, 21, 4,
            199, 35, 195, 24, 150, 5, 154, 7, 18, 128, 226, 235, 39,
            178, 117, 9, 131, 44, 26, 27, 110, 90, 160, 82, 59, 214,
            179, 41, 227, 47, 132, 83, 209, 0, 237, 32, 252, 177, 91,
            106, 203, 190, 57, 74, 76, 88, 207, 208, 239, 170, 251,
            67, 77, 51, 133, 69, 249, 2, 127, 80, 60, 159, 168, 81,
            163, 64, 143, 146, 157, 56, 245, 188, 182, 218, 33, 16,
            255, 243, 210, 205, 12, 19, 236, 95, 151, 68, 23, 196,
            167, 126, 61, 100, 93, 25, 115, 96, 129, 79, 220, 34, 42,
            144, 136, 70, 238, 184, 20, 222, 94, 11, 219, 224, 50, 58,
            10, 73, 6, 36, 92, 194, 211, 172, 98, 145, 149, 228, 121,
            231, 200, 55, 109, 141, 213, 78, 169, 108, 86, 244, 234,
            101, 122, 174, 8, 186, 120, 37, 46, 28, 166, 180, 198,
            232, 221, 116, 31, 75, 189, 139, 138, 112, 62, 181, 102,
            72, 3, 246, 14, 97, 53, 87, 185, 134, 193, 29, 158, 225,
            248, 152, 17, 105, 217, 142, 148, 155, 30, 135, 233, 206,
            85, 40, 223, 140, 161, 137, 13, 191, 230, 66, 104, 65,
            153, 45, 15, 176, 84, 187, 22)

            sage: sr.sub_byte(0)
            a^6 + a^5 + a + 1

            sage: sage_eval(str(sr.sub_byte(0)), {'a':2})
            99
            sage: S(0)
            99

            sage: sr.sub_byte(1)
            a^6 + a^5 + a^4 + a^3 + a^2

            sage: sage_eval(str(sr.sub_byte(1)), {'a':2})
            124

            sage: S(1)
            124

            sage: sr = mq.SR(1,2,2,4, allow_zero_inversions=True)
            sage: S = sr.sbox(inversion_only=True); S
            (0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8)

            sage: S(0)
            0
            sage: S(1)
            1

            sage: S(sr.k.gen())
            a^3 + a + 1
        """
        from sage.crypto.sbox import SBox

        k = self.base_ring()
        if not inversion_only:
            with AllowZeroInversionsContext(self):
                S = [self.sub_byte(elem) for elem in sorted(k)]
            return SBox(S)
        else:
            e = self.e
            S = [elem ** (2**e - 2) for elem in sorted(k)]
            return SBox(S)

    def shift_rows(self, d):
        r"""
        Perform the ``ShiftRows`` operation on ``d``.

        INPUT:

        - ``d`` -- state array or something coercible to a state array

        EXAMPLES::

            sage: sr = mq.SR(10, 4, 4, 4)
            sage: E = sr.state_array() + 1; E
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

        ::

            sage: sr.shift_rows(E)
            [1 0 0 0]
            [1 0 0 0]
            [1 0 0 0]
            [1 0 0 0]
        """
        d = self.state_array(d)
        ret = []
        for i in range(d.nrows()):
            ret += list(d.row(i)[i % d.ncols():]) + list(d.row(i)[:i % d.ncols()])
        return matrix(self.base_ring(), self._r, self._c, ret)

    def mix_columns(self, d):
        r"""
        Perform the ``MixColumns`` operation on
        ``d``.

        INPUT:

        - ``d`` -- state array or something coercible to a state array

        EXAMPLES::

            sage: sr = mq.SR(10, 4, 4, 4)
            sage: E = sr.state_array() + 1; E
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

        ::

            sage: sr.mix_columns(E)
            [    a a + 1     1     1]
            [    1     a a + 1     1]
            [    1     1     a a + 1]
            [a + 1     1     1     a]
        """
        d = self.state_array(d)
        k = self.base_ring()
        a = k.gen()
        r = self._r
        if r == 1:
            M = matrix(self.base_ring(), 1, 1, [[1]])
        elif r == 2:
            M = matrix(self.base_ring(), 2, 2, [[a + 1, a],
                                              [a, a + 1]])

        elif r == 4:
            M = matrix(self.base_ring(), 4, 4, [[a, a+1, 1, 1],
                                              [1, a, a+1, 1],
                                              [1, 1, a, a+1],
                                              [a+1, 1, 1, a]])
        ret = []
        for column in d.columns():
            ret.append(M * column)
        # AES uses the column major ordering
        return matrix(k, d.ncols(), d.nrows(), ret).transpose()

    def add_round_key(self, d, key):
        r"""
        Perform the ``AddRoundKey`` operation on
        ``d`` using ``key``.

        INPUT:

        - ``d`` -- state array or something coercible to a state array

        - ``key`` -- state array or something coercible to a state array

        EXAMPLES::

            sage: sr = mq.SR(10, 4, 4, 4)
            sage: D = sr.random_state_array()
            sage: K = sr.random_state_array()
            sage: sr.add_round_key(D, K) == K + D
            True
        """
        d = self.state_array(d)
        key = self.state_array(key)

        return d+key

    def state_array(self, d=None):
        """
        Convert the parameter to a state array.

        INPUT:

        - ``d`` -- a matrix, a list, or a tuple (default: ``None``)

        EXAMPLES::

            sage: sr = mq.SR(2, 2, 2, 4)
            sage: k = sr.base_ring()
            sage: e1 = [k.from_integer(e) for e in range(2*2)]; e1
            [0, 1, a, a + 1]
            sage: e2 = sr.phi( Matrix(k, 2*2, 1, e1) )
            sage: sr.state_array(e1) # note the column major ordering
            [    0     a]
            [    1 a + 1]
            sage: sr.state_array(e2)
            [    0     a]
            [    1 a + 1]

        ::

            sage: sr.state_array()
            [0 0]
            [0 0]
        """
        r = self.r
        c = self.c
        e = self.e
        k = self.base_ring()

        if d is None:
            return matrix(k, r, c)

        if isinstance(d, Matrix):
            if d.nrows() == r*c*e:
                return matrix(k, c, r, self.antiphi(d).list()).transpose()
            elif d.ncols() == c and d.nrows() == r and d.base_ring() == k:
                return d

        if isinstance(d, tuple([list, tuple])):
            return matrix(k, c, r, d).transpose()

    def is_state_array(self, d):
        """
        Return ``True`` if ``d`` is a state array, i.e. has the correct
        dimensions and base field.

        EXAMPLES::

            sage: sr = mq.SR(2, 2, 4, 8)
            sage: k = sr.base_ring()
            sage: sr.is_state_array( matrix(k, 2, 4) )
            True

        ::

            sage: sr = mq.SR(2, 2, 4, 8)
            sage: k = sr.base_ring()
            sage: sr.is_state_array( matrix(k, 4, 4) )
            False
        """
        return isinstance(d, Matrix) and \
               d.nrows() == self.r and \
               d.ncols() == self.c and \
               d.base_ring() == self.base_ring()

    def random_state_array(self, *args, **kwds):
        r"""
        Return a random element in ``MatrixSpace(self.base_ring(),
        self.r, self.c)``.

        EXAMPLES::

            sage: sr = mq.SR(2, 2, 2, 4)
            sage: sr.random_state_array().parent()
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field in a of size 2^4
        """
        return random_matrix(self.base_ring(), self._r, self._c, *args, **kwds)

    def random_vector(self, *args, **kwds):
        r"""
        Return a random vector as it might appear in the algebraic
        expression of ``self``.

        EXAMPLES::

            sage: mq.SR(2, 2, 2, 4).random_vector().parent()
            Full MatrixSpace of 16 by 1 dense matrices over Finite Field in a of size 2^4

        .. NOTE::

           `\phi` was already applied to the result.
        """
        return self.vector(self.random_state_array(*args, **kwds))

    def random_element(self, elem_type='vector', *args, **kwds):
        """
        Return a random element for ``self``.  Other arguments and keywords are
        passed to random_* methods.

        INPUT:

        - ``elem_type`` -- either 'vector' or 'state array' (default: ``'vector'``)

        EXAMPLES::

            sage: sr = mq.SR()
            sage: sr.random_element().parent()
            Full MatrixSpace of 4 by 1 dense matrices over Finite Field in a of size 2^4
            sage: sr.random_element('state_array').parent()
            Full MatrixSpace of 1 by 1 dense matrices over Finite Field in a of size 2^4

        Passes extra positional or keyword arguments through::

            sage: sr.random_element(density=0)
            [0]
            [0]
            [0]
            [0]
        """
        if elem_type == "vector":
            return self.random_vector(*args, **kwds)
        elif elem_type == "state_array":
            return self.random_state_array(*args, **kwds)
        else:
            raise TypeError("parameter type not understood")

    def key_schedule(self, kj, i):
        """
        Return `k_i` for a given `i` and `k_j`
        with `j = i-1`.

        EXAMPLES::

            sage: sr = mq.SR(10, 4, 4, 8, star=True, allow_zero_inversions=True)
            sage: ki = sr.state_array()
            sage: for i in range(10):
            ....:     ki = sr.key_schedule(ki, i+1)
            sage: print(sr.hex_str_matrix(ki))
            B4 3E 23 6F
            EF 92 E9 8F
            5B E2 51 18
            CB 11 CF 8E
        """
        if i < 0:
            raise TypeError("i must be >= i")

        if i == 0:
            return kj

        r = self.r
        c = self.c
        F = self.base_ring()
        a = F.gen()
        SubByte = self.sub_byte

        rc = matrix(F, r, c, ([a**(i-1)] * c) + [F(0)]*((r-1)*c) )
        ki = matrix(F, r, c)

        if r == 1:
            s0 = SubByte(kj[0, c-1])

            if c > 1:
                for q in range(c):
                    ki[0, q] = s0 + sum([kj[0, t] for t in range(q+1) ])
            else:
                ki[0, 0] = s0

        elif r == 2:
            s0 = SubByte(kj[1, c-1])
            s1 = SubByte(kj[0, c-1])

            if c > 1:
                for q in range(c):
                    ki[0, q] = s0 + sum([ kj[0, t] for t in range(q+1) ])
                    ki[1, q] = s1 + sum([ kj[1, t] for t in range(q+1) ])
            else:
                ki[0, 0] = s0
                ki[1, 0] = s1

        elif r == 4:

            if self._aes_mode:
                s0 = SubByte(kj[1, c-1])
                s1 = SubByte(kj[2, c-1])
                s2 = SubByte(kj[3, c-1])
                s3 = SubByte(kj[0, c-1])
            else:
                s0 = SubByte(kj[3, c-1])
                s1 = SubByte(kj[2, c-1])
                s2 = SubByte(kj[1, c-1])
                s3 = SubByte(kj[0, c-1])

            if c > 1:
                for q in range(c):
                    ki[0, q] = s0 + sum([ kj[0, t] for t in range(q+1) ])
                    ki[1, q] = s1 + sum([ kj[1, t] for t in range(q+1) ])
                    ki[2, q] = s2 + sum([ kj[2, t] for t in range(q+1) ])
                    ki[3, q] = s3 + sum([ kj[3, t] for t in range(q+1) ])

            else:
                ki[0, 0] = s0
                ki[1, 0] = s1
                ki[2, 0] = s2
                ki[3, 0] = s3
        ki += rc

        return ki

    def __call__(self, P, K):
        r"""
        Encrypts the plaintext `P` using the key `K`.

        Both must be given as state arrays or coercible to state arrays.

        INPUT:

        - ``P`` -- plaintext as state array or something coercible to a
          qstate array

        - ``K`` -- key as state array or something coercible to a state
          array

        TESTS:

        The official AES test vectors::

            sage: sr = mq.SR(10, 4, 4, 8, star=True, allow_zero_inversions=True)
            sage: k = sr.base_ring()
            sage: plaintext = sr.state_array([k.from_integer(e) for e in range(16)])
            sage: key = sr.state_array([k.from_integer(e) for e in range(16)])
            sage: print(sr.hex_str_matrix( sr(plaintext, key) ))
            0A 41 F1 C6
            94 6E C3 53
            0B F0 94 EA
            B5 45 58 5A

        Brian Gladman's development vectors (dev_vec.txt)::

            sage: sr = mq.SR(10, 4, 4, 8, star=True, allow_zero_inversions=True, aes_mode=True)
            sage: k = sr.base_ring()
            sage: plain = '3243f6a8885a308d313198a2e0370734'
            sage: key = '2b7e151628aed2a6abf7158809cf4f3c'
            sage: from sage.misc.verbose import set_verbose
            sage: set_verbose(2)
            sage: cipher = sr(plain, key)
            R[01].start   193DE3BEA0F4E22B9AC68D2AE9F84808
            R[01].s_box   D42711AEE0BF98F1B8B45DE51E415230
            R[01].s_row   D4BF5D30E0B452AEB84111F11E2798E5
            R[01].m_col   046681E5E0CB199A48F8D37A2806264C
            R[01].k_sch   A0FAFE1788542CB123A339392A6C7605
            R[02].start   A49C7FF2689F352B6B5BEA43026A5049
            R[02].s_box   49DED28945DB96F17F39871A7702533B
            R[02].s_row   49DB873B453953897F02D2F177DE961A
            R[02].m_col   584DCAF11B4B5AACDBE7CAA81B6BB0E5
            R[02].k_sch   F2C295F27A96B9435935807A7359F67F
            R[03].start   AA8F5F0361DDE3EF82D24AD26832469A
            R[03].s_box   AC73CF7BEFC111DF13B5D6B545235AB8
            R[03].s_row   ACC1D6B8EFB55A7B1323CFDF457311B5
            R[03].m_col   75EC0993200B633353C0CF7CBB25D0DC
            R[03].k_sch   3D80477D4716FE3E1E237E446D7A883B
            R[04].start   486C4EEE671D9D0D4DE3B138D65F58E7
            R[04].s_box   52502F2885A45ED7E311C807F6CF6A94
            R[04].s_row   52A4C89485116A28E3CF2FD7F6505E07
            R[04].m_col   0FD6DAA9603138BF6FC0106B5EB31301
            R[04].k_sch   EF44A541A8525B7FB671253BDB0BAD00
            R[05].start   E0927FE8C86363C0D9B1355085B8BE01
            R[05].s_box   E14FD29BE8FBFBBA35C89653976CAE7C
            R[05].s_row   E1FB967CE8C8AE9B356CD2BA974FFB53
            R[05].m_col   25D1A9ADBD11D168B63A338E4C4CC0B0
            R[05].k_sch   D4D1C6F87C839D87CAF2B8BC11F915BC
            R[06].start   F1006F55C1924CEF7CC88B325DB5D50C
            R[06].s_box   A163A8FC784F29DF10E83D234CD503FE
            R[06].s_row   A14F3DFE78E803FC10D5A8DF4C632923
            R[06].m_col   4B868D6D2C4A8980339DF4E837D218D8
            R[06].k_sch   6D88A37A110B3EFDDBF98641CA0093FD
            R[07].start   260E2E173D41B77DE86472A9FDD28B25
            R[07].s_box   F7AB31F02783A9FF9B4340D354B53D3F
            R[07].s_row   F783403F27433DF09BB531FF54ABA9D3
            R[07].m_col   1415B5BF461615EC274656D7342AD843
            R[07].k_sch   4E54F70E5F5FC9F384A64FB24EA6DC4F
            R[08].start   5A4142B11949DC1FA3E019657A8C040C
            R[08].s_box   BE832CC8D43B86C00AE1D44DDA64F2FE
            R[08].s_row   BE3BD4FED4E1F2C80A642CC0DA83864D
            R[08].m_col   00512FD1B1C889FF54766DCDFA1B99EA
            R[08].k_sch   EAD27321B58DBAD2312BF5607F8D292F
            R[09].start   EA835CF00445332D655D98AD8596B0C5
            R[09].s_box   87EC4A8CF26EC3D84D4C46959790E7A6
            R[09].s_row   876E46A6F24CE78C4D904AD897ECC395
            R[09].m_col   473794ED40D4E4A5A3703AA64C9F42BC
            R[09].k_sch   AC7766F319FADC2128D12941575C006E
            R[10].s_box   E9098972CB31075F3D327D94AF2E2CB5
            R[10].s_row   E9317DB5CB322C723D2E895FAF090794
            R[10].k_sch   D014F9A8C9EE2589E13F0CC8B6630CA6
            R[10].output  3925841D02DC09FBDC118597196A0B32
            sage: set_verbose(0)
        """
        r, c, e = self.r, self.c, self.e
        F = self.base_ring()

        if isinstance(P, str):
            P = self.state_array([F.from_integer(ZZ(P[i: i + 2], 16)) for i in range(0, len(P), 2)])
        if isinstance(K, str):
            K = self.state_array([F.from_integer(ZZ(K[i: i + 2], 16)) for i in range(0, len(K), 2)])

        if self.is_state_array(P) and self.is_state_array(K):
            _type = self.state_array
        elif self.is_vector(P) and self.is_vector(K):
            _type = self.vector
        elif isinstance(P, (list,tuple)) and isinstance(K, (list,tuple)):
            if len(P) == len(K) == r*c:
                _type = self.state_array
            elif len(P) == len(K) == r*c*e:
                _type = self.vector
            else:
                raise TypeError("length %d or %d doesn't match either %d or %d" % (len(P),len(K),r*c,r*c*e))
        else:
            raise TypeError("plaintext or key parameter not understood")

        P = self.state_array(P)
        K = self.state_array(K)

        AddRoundKey = self.add_round_key
        SubBytes = self.sub_bytes
        MixColumns = self.mix_columns
        ShiftRows = self.shift_rows
        KeyExpansion = self.key_schedule

        P = AddRoundKey(P, K)

        for r in range(self._n-1):
            if get_verbose() >= 2:
                print("R[%02d].start   %s" % (r+1, self.hex_str_vector(P)))

            P = SubBytes(P)
            if get_verbose() >= 2:
                print("R[%02d].s_box   %s" % (r+1, self.hex_str_vector(P)))

            P = ShiftRows(P)
            if get_verbose() >= 2:
                print("R[%02d].s_row   %s" % (r+1, self.hex_str_vector(P)))

            P = MixColumns(P)
            if get_verbose() >= 2:
                print("R[%02d].m_col   %s" % (r+1, self.hex_str_vector(P)))

            K = KeyExpansion(K, r+1)
            if get_verbose() >= 2:
                print("R[%02d].k_sch   %s" % (r+1, self.hex_str_vector(K)))

            P = AddRoundKey(P, K)

        P = SubBytes(P)
        if get_verbose() >= 2:
            print("R[%02d].s_box   %s" % (self.n, self.hex_str_vector(P)))

        P = ShiftRows(P)
        if get_verbose() >= 2:
            print("R[%02d].s_row   %s" % (self.n, self.hex_str_vector(P)))

        if not self._star:
            P = MixColumns(P)
            if get_verbose() >= 2:
                print("R[%02d].m_col   %s" % (self.n, self.hex_str_vector(P)))

        K = KeyExpansion(K, self._n)
        if get_verbose() >= 2:
            print("R[%02d].k_sch   %s" % (self.n, self.hex_str_vector(K)))

        P = AddRoundKey(P, K)
        if get_verbose() >= 2:
            print("R[%02d].output  %s" % (self.n, self.hex_str_vector(P)))

        return _type(P)

    def hex_str(self, M, typ='matrix'):
        r"""
        Return a hex string for the provided AES state array/matrix.

        INPUT:

        - ``M`` -- state array

        - ``typ`` -- controls what to return, either 'matrix'
          or 'vector' (default: ``'matrix'``)

        EXAMPLES::

            sage: sr = mq.SR(2, 2, 2, 4)
            sage: k = sr.base_ring()
            sage: A = matrix(k, 2, 2, [1, k.gen(), 0, k.gen()^2])
            sage: sr.hex_str(A)
            ' 1 2 \n 0 4 \n'

        ::

            sage: sr.hex_str(A, typ='vector')
            '1024'
        """
        if typ == "matrix":
            return self.hex_str_matrix(M)
        elif typ == "vector":
            return self.hex_str_vector(M)
        else:
            raise TypeError("parameter type must either be 'matrix' or 'vector'")

    def hex_str_matrix(self, M):
        r"""
        Return a two-dimensional AES-like representation of the matrix M.

        That is, show the finite field elements as hex strings.

        INPUT:

        - ``M`` -- an AES state array

        EXAMPLES::

            sage: sr = mq.SR(2, 2, 2, 4)
            sage: k = sr.base_ring()
            sage: A = matrix(k, 2, 2, [1, k.gen(), 0, k.gen()^2])
            sage: sr.hex_str_matrix(A)
            ' 1 2 \n 0 4 \n'
        """
        e = M.base_ring().degree()
        st = [""]
        for x in range(M.nrows()):
            for y in range(M.ncols()):
                if e == 8:
                    st.append("%02X" % M[x, y].to_integer())
                else:
                    st.append("%X" % M[x, y].to_integer())
            st.append("\n")
        return " ".join(st)

    def hex_str_vector(self, M):
        """
        Return a one-dimensional AES-like representation of the matrix M.

        That is, show the finite field elements as hex strings.

        INPUT:

        - ``M`` -- an AES state array

        EXAMPLES::

            sage: sr = mq.SR(2, 2, 2, 4)
            sage: k = sr.base_ring()
            sage: A = matrix(k, 2, 2, [1, k.gen(), 0, k.gen()^2])
            sage: sr.hex_str_vector(A)
            '1024'
        """
        e = M.base_ring().degree()
        st = [""]
        for y in range(M.ncols()):
            for x in range(M.nrows()):
                if e == 8:
                    st.append("%02X" % M[x, y].to_integer())
                else:
                    st.append("%X" % M[x, y].to_integer())
            #st.append("\n")
        return "".join(st)

    def _insert_matrix_into_matrix(self, dst, src, row, col):
        """
        Insert matrix src into matrix dst starting at row and col.

        INPUT:

        - ``dst`` -- a matrix

        - ``src`` -- a matrix

        - ``row`` -- offset row

        - ``col`` -- offset columns

        EXAMPLES::

            sage: sr = mq.SR(10, 4, 4, 4)
            sage: a = sr.k.gen()
            sage: A = sr.state_array() + 1; A
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            sage: B = Matrix(sr.base_ring(), 2, 2, [0, a, a+1, a^2]); B
            [    0     a]
            [a + 1   a^2]
            sage: sr._insert_matrix_into_matrix(A, B, 1, 1)
            [    1     0     0     0]
            [    0     0     a     0]
            [    0 a + 1   a^2     0]
            [    0     0     0     1]
        """
        for i in range(src.nrows()):
            for j in range(src.ncols()):
                dst[row+i, col+j] = src[i, j]
        return dst

    def varformatstr(self, name, n=None, rc=None, e=None):
        r"""
        Return a format string which is understood by print et al.

        If a numerical value is omitted, the default value of ``self``
        is used.  The numerical values (``n``, ``rc``, ``e``) are used
        to determine the width of the respective fields in the format
        string.

        INPUT:

        - ``name`` -- name of the variable
        - ``n`` -- number of rounds (default: ``None``)
        - ``rc`` -- number of rows \* number of cols (default: ``None``)
        - ``e`` -- exponent of base field (default: ``None``)

        EXAMPLES::

            sage: sr = mq.SR(1, 2, 2, 4)
            sage: sr.varformatstr('x')
            'x%01d%01d%01d'
            sage: sr.varformatstr('x', n=1000)
            'x%03d%03d%03d'
        """
        if n is None:
            n = self.n
        if rc is None:
            rc = self.r * self.c
        if e is None:
            e = self.e

        l = str(max([  len(str(rc-1)), len(str(n-1)), len(str(e-1)) ] ))
        if name not in ("k", "s"):
            pf = self._postfix
        else:
            pf = ""
        format_string = name + pf + "%0" + l + "d" + "%0" + l + "d" + "%0" + l + "d"
        return format_string

    def varstr(self, name, nr, rc, e):
        """
        Return a string representing a variable for the small scale
        AES subject to the given constraints.

        INPUT:

        - ``name`` -- variable name
        - ``nr`` -- number of round to create variable strings for
        - ``rc`` -- row*column index in state array
        - ``e`` -- exponent of base field

        EXAMPLES::

            sage: sr = mq.SR(10, 1, 2, 4)
            sage: sr.varstr('x', 2, 1, 1)
            'x211'
        """
        format_string = self.varformatstr(name, self.n, self.r*self.c, self.e)
        return format_string % (nr, rc, e)

    def varstrs(self, name, nr, rc=None, e=None):
        """
        Return a list of strings representing variables in ``self``.

        INPUT:

        - ``name`` -- variable name
        - ``nr`` -- number of round to create variable strings for
        - ``rc`` -- number of rows * number of columns in the state array (default: ``None``)
        - ``e`` -- exponent of base field (default: ``None``)

        EXAMPLES::

            sage: sr = mq.SR(10, 1, 2, 4)
            sage: sr.varstrs('x', 2)
            ('x200', 'x201', 'x202', 'x203', 'x210', 'x211', 'x212', 'x213')
        """
        if rc is None:
            rc = self.r * self.c

        if e is None:
            e = self.e

        n = self._n

        format_string = self.varformatstr(name, n, rc, e)

        return tuple([format_string % (nr, rci, ei) for rci in range(rc) for ei in range(e)])

    def vars(self, name, nr, rc=None, e=None):
        """
        Return a list of variables in ``self``.

        INPUT:

        - ``name`` -- variable name
        - ``nr`` -- number of round to create variable strings for
        - ``rc`` -- number of rounds * number of columns in the state array (default: ``None``)
        - ``e`` -- exponent of base field (default: ``None``)

        EXAMPLES::

            sage: sr = mq.SR(10, 1, 2, 4)
            sage: sr.vars('x', 2)
            (x200, x201, x202, x203, x210, x211, x212, x213)
        """
        gd = self.variable_dict()
        return tuple([gd[s] for s in self.varstrs(name, nr, rc, e)])

    def variable_dict(self):
        """
        Return a dictionary to access variables in ``self.R`` by their
        names.

        EXAMPLES::

            sage: sr = mq.SR(1,1,1,4)
            sage: sr.variable_dict()
            {'k000': k000,
             'k001': k001,
             'k002': k002,
             'k003': k003,
             'k100': k100,
             'k101': k101,
             'k102': k102,
             'k103': k103,
             's000': s000,
             's001': s001,
             's002': s002,
             's003': s003,
             'w100': w100,
             'w101': w101,
             'w102': w102,
             'w103': w103,
             'x100': x100,
             'x101': x101,
             'x102': x102,
             'x103': x103}

            sage: sr = mq.SR(1,1,1,4,gf2=True)
            sage: sr.variable_dict()                                                    # needs sage.rings.polynomial.pbori
            {'k000': k000,
             'k001': k001,
             'k002': k002,
             'k003': k003,
             'k100': k100,
             'k101': k101,
             'k102': k102,
             'k103': k103,
             's000': s000,
             's001': s001,
             's002': s002,
             's003': s003,
             'w100': w100,
             'w101': w101,
             'w102': w102,
             'w103': w103,
             'x100': x100,
             'x101': x101,
             'x102': x102,
             'x103': x103}
        """
        try:
            R,gd = self._variable_dict
            if R is self.R:
                return gd
            else:
                pass
        except AttributeError:
            pass

        gd = self.R.gens_dict()
        self._variable_dict = self.R,gd
        return gd

    def block_order(self):
        """
        Return a block order for ``self`` where each round is a block.

        EXAMPLES::

            sage: sr = mq.SR(2, 1, 1, 4)
            sage: sr.block_order()
            Block term order with blocks:
            (Degree lexicographic term order of length 16,
             Degree lexicographic term order of length 16,
             Degree lexicographic term order of length 4)

        ::

            sage: P = sr.ring(order='block')
            sage: print(P.repr_long())
            Polynomial Ring
              Base Ring : Finite Field in a of size 2^4
                   Size : 36 Variables
               Block  0 : Ordering : deglex
                          Names    : k200, k201, k202, k203, x200, x201, x202, x203, w200, w201, w202, w203, s100, s101, s102, s103
               Block  1 : Ordering : deglex
                          Names    : k100, k101, k102, k103, x100, x101, x102, x103, w100, w101, w102, w103, s000, s001, s002, s003
               Block  2 : Ordering : deglex
                          Names    : k000, k001, k002, k003
        """
        r = self.r
        c = self.c
        e = self.e
        n = self.n

        T = None
        for _n in range(n):
            T = TermOrder('deglex', r*e + 3*r*c*e ) + T

        T += TermOrder('deglex', r*c*e)

        return T

    def ring(self, order=None, reverse_variables=None):
        r"""
        Construct a ring as a base ring for the polynomial system.

        By default, variables are ordered in the reverse of their natural
        ordering, i.e. the reverse of as they appear.

        INPUT:

        - ``order`` -- a monomial ordering (default: ``None``)
        - ``reverse_variables`` -- reverse rounds of variables (default: ``True``)

        The variable assignment is as follows:

        - `k_{i,j,l}` -- subkey round `i` word `j` conjugate/bit `l`
        - `s_{i,j,l}` -- subkey inverse round `i` word `j` conjugate/bit `l`
        - `w_{i,j,l}` -- inversion input round `i` word `j` conjugate/bit `l`
        - `x_{i,j,l}` -- inversion output round `i` word `j` conjugate/bit `l`


        Note that the variables are ordered in column major ordering
        in the state array and that the bits are ordered in little
        endian ordering.

        For example, if `x_{0,1,0}` is a variable over `\GF{2}` for
        `r=2` and `c=2` then refers to the *most* significant bit of
        the entry in the position (1,0) in the state array matrix.

        EXAMPLES::

            sage: sr = mq.SR(2, 1, 1, 4)
            sage: P = sr.ring(order='block')
            sage: print(P.repr_long())
            Polynomial Ring
              Base Ring : Finite Field in a of size 2^4
                   Size : 36 Variables
               Block  0 : Ordering : deglex
                          Names    : k200, k201, k202, k203, x200, x201, x202, x203, w200, w201, w202, w203, s100, s101, s102, s103
               Block  1 : Ordering : deglex
                          Names    : k100, k101, k102, k103, x100, x101, x102, x103, w100, w101, w102, w103, s000, s001, s002, s003
               Block  2 : Ordering : deglex
                          Names    : k000, k001, k002, k003
        """
        r = self.r
        c = self.c
        e = self.e
        n = self.n
        if not self._gf2:
            k = self.base_ring()
        else:
            k = GF(2)

        if order is not None:
            self._order = order
        if self._order == 'block':
            self._order = self.block_order()

        if reverse_variables is None:
            reverse_variables = self._reverse_variables

        if reverse_variables:
            process = reversed
        else:
            process = lambda x: x

        if reverse_variables:
            names = []
        else:
            names = self.varstrs("k", 0, r*c, e)

        for _n in process(list(range(n))):
            names += self.varstrs("k", _n+1, r*c, e)
            names += self.varstrs("x", _n+1, r*c, e)
            names += self.varstrs("w", _n+1, r*c, e)
            names += self.varstrs("s", _n, r, e)

        if reverse_variables:
            names += self.varstrs("k", 0, r*c, e)

        #from sage.rings.polynomial.pbori.pbori import BooleanPolynomialRing

        if self._gf2 and self._polybori:
            return BooleanPolynomialRing(2*n*r*c*e + (n+1)*r*c*e + n*r*e, names, order=self._order)
        else:
            return PolynomialRing(k, 2*n*r*c*e + (n+1)*r*c*e + n*r*e, names, order=self._order)

    def round_polynomials(self, i, plaintext=None, ciphertext=None):
        r"""
        Return list of polynomials for a given round `i`.

        If ``i == 0`` a plaintext must be provided, if ``i == n`` a
        ciphertext must be provided.

        INPUT:

        - ``i`` -- round number

        - ``plaintext`` -- plaintext (optional, mandatory in first round)

        - ``ciphertext`` -- ciphertext (optional, mandatory in last round)

        OUTPUT: tuple

        EXAMPLES::

            sage: sr = mq.SR(1, 1, 1, 4)
            sage: k = sr.base_ring()
            sage: p = [k.random_element() for _ in range(sr.r*sr.c)]
            sage: sr.round_polynomials(0, plaintext=p)
            (w100 + k000..., w101 + k001..., w102 + k002..., w103 + k003...)
        """
        r = self._r
        c = self._c
        e = self._e
        n = self._n
        R = self.R

        M = self.M

        _vars = self.vars

        if i == 0:
            w1 = matrix(R, r*c*e, 1, _vars("w", 1, r*c, e))
            k0 = matrix(R, r*c*e, 1, _vars("k", 0, r*c, e))
            if isinstance(plaintext, (tuple, list)) and len(plaintext) == r*c:
                plaintext = matrix(R, r*c*e, 1, self.phi(plaintext))
            return tuple((w1 + k0 + plaintext).list())

        elif i > 0 and i <= n:

            if self._star and i == n:
                M = self.Mstar

            xj = matrix(R, r*c*e, 1, _vars("x", i, r*c, e))
            ki = matrix(R, r*c*e, 1, _vars("k", i, r*c, e))
            rcon = matrix(R, r*c*e, 1, self.phi([self.sbox_constant()]*r*c))

            if i < n:
                wj = matrix(R, r*c*e, 1, _vars("w", i+1, r*c, e))
            if i == n:
                if isinstance(ciphertext, (tuple, list)) and len(ciphertext) == r*c:
                    ciphertext = matrix(R, r*c*e, 1, self.phi(ciphertext))
                wj = ciphertext

            lin = (wj + ki + M * xj + rcon).list()

            wi = matrix(R, r*c*e, 1, _vars("w", i, r*c, e))
            xi = matrix(R, r*c*e, 1, _vars("x", i, r*c, e))
            sbox = []
            sbox += self.inversion_polynomials(xi, wi, r*c*e)
            sbox += self.field_polynomials("x", i)
            sbox += self.field_polynomials("w", i)
            return tuple(lin + sbox)

    def key_schedule_polynomials(self, i):
        r"""
        Return polynomials for the `i`-th round of the key
        schedule.

        INPUT:

        - ``i`` -- round (`0 \leq i \leq n`)

        EXAMPLES::

            sage: sr = mq.SR(1, 1, 1, 4, gf2=True, polybori=False)

        The `0`-th subkey is the user provided key, so only conjugacy
        relations or field polynomials are added.::

            sage: sr.key_schedule_polynomials(0)
            (k000^2 + k000, k001^2 + k001, k002^2 + k002, k003^2 + k003)

        The 1st subkey is derived from the user provided key according to
        the key schedule which is non-linear.::

            sage: sr.key_schedule_polynomials(1)
            (k100 + s000 + s002 + s003,
             k101 + s000 + s001 + s003 + 1,
             k102 + s000 + s001 + s002 + 1,
             k103 + s001 + s002 + s003 + 1,
             k100^2 + k100, k101^2 + k101, k102^2 + k102, k103^2 + k103,
             s000^2 + s000, s001^2 + s001, s002^2 + s002, s003^2 + s003,
             s000*k000 + s000*k003 + s001*k002 + s002*k001 + s003*k000,
             s000*k000 + s000*k001 + s001*k000 + s001*k003 + s002*k002 + s003*k001,
             s000*k001 + s000*k002 + s001*k000 + s001*k001 + s002*k000 + s002*k003 + s003*k002,
             s000*k000 + s000*k001 + s000*k003 + s001*k001 + s002*k000 + s002*k002 + s003*k000 + k000,
             s000*k002 + s001*k000 + s001*k001 + s001*k003 + s002*k001 + s003*k000 + s003*k002 + k001,
             s000*k000 + s000*k001 + s000*k002 + s001*k002 + s002*k000 + s002*k001 + s002*k003 + s003*k001 + k002,
             s000*k001 + s001*k000 + s001*k002 + s002*k000 + s003*k001 + s003*k003 + k003,
             s000*k000 + s000*k002 + s000*k003 + s001*k000 + s001*k001 + s002*k002 + s003*k000 + s000,
             s000*k001 + s000*k003 + s001*k001 + s001*k002 + s002*k000 + s002*k003 + s003*k001 + s001,
             s000*k000 + s000*k002 + s001*k000 + s001*k002 + s001*k003 + s002*k000 + s002*k001 + s003*k002 + s002,
             s000*k001 + s000*k002 + s001*k000 + s001*k003 + s002*k001 + s003*k003 + s003,
             s000*k002 + s001*k001 + s002*k000 + s003*k003 + 1)
        """
        R = self.R
        r = self.r
        e = self.e
        c = self.c
        k = self.k
        a = k.gen()

        if i < 0:
            raise TypeError("i must by >= 0")

        if i == 0:
            return tuple(self.field_polynomials("k", i, r*c))
        else:
            L = self.lin_matrix(r)
            ki = matrix(R, r*c*e, 1, self.vars("k", i  , r*c, e))
            kj = matrix(R, r*c*e, 1, self.vars("k", i-1, r*c, e))
            si = matrix(R, r*e, 1, self.vars("s", i-1, r, e))

            rc = matrix(R, r*e, 1, self.phi([a**(i-1)] + [k(0)]*(r-1)) )
            d = matrix(R, r*e, 1, self.phi([self.sbox_constant()]*r) )

            sbox = []

            sbox += self.field_polynomials("k", i)
            sbox += self.field_polynomials("s", i-1, r)

            if r == 1:
                sbox += self.inversion_polynomials(kj[(c - 1)*e:(c - 1)*e + e], si[0:e], e)
            if r == 2:
                sbox += self.inversion_polynomials( kj[(2*c - 1)*e : (2*c - 1)*e + e] , si[0:1*e], e )
                sbox += self.inversion_polynomials( kj[(2*c - 2)*e : (2*c - 2)*e + e] , si[e:2*e], e )
            if r == 4:
                if self._aes_mode:
                    sbox += self.inversion_polynomials( kj[(4*c-3)*e  : (4*c-3)*e + e] , si[0*e : 1*e] , e )
                    sbox += self.inversion_polynomials( kj[(4*c-2)*e  : (4*c-2)*e + e] , si[1*e : 2*e] , e )
                    sbox += self.inversion_polynomials( kj[(4*c-1)*e  : (4*c-1)*e + e] , si[2*e : 3*e] , e )
                    sbox += self.inversion_polynomials( kj[(4*c-4)*e  : (4*c-4)*e + e] , si[3*e : 4*e] , e )
                else:
                    sbox += self.inversion_polynomials( kj[(4*c-1)*e  : (4*c-1)*e + e] , si[0*e : 1*e] , e )
                    sbox += self.inversion_polynomials( kj[(4*c-2)*e  : (4*c-2)*e + e] , si[1*e : 2*e] , e )
                    sbox += self.inversion_polynomials( kj[(4*c-3)*e  : (4*c-3)*e + e] , si[2*e : 3*e] , e )
                    sbox += self.inversion_polynomials( kj[(4*c-4)*e  : (4*c-4)*e + e] , si[3*e : 4*e] , e )

            si = L * si + d + rc
            Sum = matrix(R, r*e, 1)
            lin = []
            if c > 1:
                for q in range(c):
                    t = list(range(r*e*(q) , r*e*(q+1)))
                    Sum += kj.matrix_from_rows(t)
                    lin += (ki.matrix_from_rows(t) + si + Sum).list()

            else:
                lin += (ki + si).list()
            return tuple(lin + sbox)

    def polynomial_system(self, P=None, K=None, C=None):
        """
        Return a polynomial system for this small scale AES variant for a
        given plaintext-key pair.

        If neither ``P``, ``K`` nor ``C`` are provided, a random pair
        (``P``, ``K``) will be generated. If ``P`` and ``C`` are
        provided no ``K`` needs to be provided.

        INPUT:

        - ``P`` -- vector, list, or tuple (default: ``None``)
        - ``K`` -- vector, list, or tuple (default: ``None``)
        - ``C`` -- vector, list, or tuple (default: ``None``)

        EXAMPLES::

            sage: sr = mq.SR(1, 1, 1, 4, gf2=True, polybori=True)
            sage: P = sr.vector([0, 0, 1, 0])
            sage: K = sr.vector([1, 0, 0, 1])
            sage: F, s = sr.polynomial_system(P, K)                                     # needs sage.rings.polynomial.pbori

        This returns a polynomial system::

            sage: F                                                                     # needs sage.rings.polynomial.pbori
            Polynomial Sequence with 36 Polynomials in 20 Variables

        and a solution::

            sage: s  # random -- maybe we need a better doctest here?                   # needs sage.rings.polynomial.pbori
            {k000: 1, k001: 0, k003: 1, k002: 0}

        This solution is not the only solution that we can learn from the
        Groebner basis of the system.

        ::

            sage: F.groebner_basis()[-3:]                                               # needs sage.rings.polynomial.pbori
            [k000 + 1, k001, k003 + 1]

        In particular we have two solutions::

            sage: len(F.ideal().variety())                                              # needs sage.rings.polynomial.pbori
            2

        In the following example we provide ``C`` explicitly::

           sage: C = sr(P,K)
           sage: F,s = sr.polynomial_system(P=P, C=C)                                   # needs sage.rings.polynomial.pbori
           sage: F                                                                      # needs sage.rings.polynomial.pbori
           Polynomial Sequence with 36 Polynomials in 20 Variables

        Alternatively, we can use symbols for the ``P`` and
        ``C``. First, we have to create a polynomial ring::

            sage: # needs sage.rings.polynomial.pbori
            sage: sr = mq.SR(1, 1, 1, 4, gf2=True, polybori=True)
            sage: R = sr.R
            sage: vn = sr.varstrs("P",0,1,4) + R.variable_names() + sr.varstrs("C",0,1,4)
            sage: R = BooleanPolynomialRing(len(vn),vn)
            sage: sr.R = R


        Now, we can construct the purely symbolic equation system::

            sage: # needs sage.rings.polynomial.pbori
            sage: C = sr.vars("C",0); C
            (C000, C001, C002, C003)
            sage: P = sr.vars("P",0)
            sage: F,s = sr.polynomial_system(P=P,C=C)
            sage: F
            Polynomial Sequence with 36 Polynomials in 28 Variables
            sage: F.part(0)
            (P000 + w100 + k000, P001 + w101 + k001, P002 + w102 + k002, P003 + w103 + k003)
            sage: F.part(-2)
            (k100 + x100 + x102 + x103 + C000, k101 + x100 + x101 + x103 + C001 + 1, ...)

        We show that the (returned) key is a solution to the returned system::

            sage: sr = mq.SR(3,4,4,8, star=True, gf2=True, polybori=True)
            sage: while True:  # workaround (see :issue:`31891`)                         # needs sage.rings.polynomial.pbori
            ....:     try:
            ....:         F, s = sr.polynomial_system()
            ....:         break
            ....:     except ZeroDivisionError:
            ....:         pass
            sage: F.subs(s).groebner_basis()    # long time                             # needs sage.rings.polynomial.pbori
            Polynomial Sequence with 1248 Polynomials in 1248 Variables
        """
        plaintext = P
        key = K
        ciphertext = C

        system = []
        n = self._n

        data = []

        R = self.R
        r,c,e = self.r,self.c,self.e

        for d in (plaintext, key, ciphertext):
            if d is None:
                data.append( None )
            elif isinstance(d, (tuple, list)):
                if isinstance(d[0], int):
                    d = [GF(2)(_) for _ in d]
                if len(d) == r*c*e and (d[0].parent() is R or d[0].parent() == R):
                    data.append( matrix(R,r*c*e,1,d) )
                    continue
                try:
                    data.append( self.phi(self.state_array(d)) )
                except ValueError: # GF2 vectors maybe?
                    data.append( self.vector(d) )
            elif self.is_state_array(d):
                data.append( self.phi(d) )
            elif self.is_vector(d):
                data.append( d )
            else:
                data.append( False )

        plaintext, key, ciphertext = data

        if plaintext is False:
            raise TypeError("type %s of P not understood" % (type(plaintext)))
        elif plaintext is None:
            plaintext = self.random_element("vector")

        if key is None:
            key = self.random_element("vector")
        elif key is False and ciphertext is False:
            raise TypeError("type %s of K not understood" % (type(key)))

        if ciphertext is None:
            ciphertext = self(plaintext, key)
        elif ciphertext is False:
            raise TypeError("type %s of C not understood" % (type(ciphertext)))

        for i in range(n+1):
            system.append( self.round_polynomials(i, plaintext, ciphertext) )
            system.append( self.key_schedule_polynomials(i) )

        if key is not None:
            K = dict(zip(self.vars("k", 0), key.list()))
        else:
            K = None
        return PolynomialSequence(self.R, system), K


class SR_gf2n(SR_generic):
    r"""
    Small Scale Variants of the AES polynomial system constructor over
    `\GF{2^n}`.
    """
    def vector(self, d=None):
        """
        Construct a vector suitable for the algebraic representation of
        SR, i.e. BES.

        INPUT:

        - ``d`` -- values for vector, must be understood by ``self.phi``
          (default: ``None``)

        EXAMPLES::

            sage: sr = mq.SR()
            sage: sr
            SR(1,1,1,4)
            sage: k = sr.base_ring()
            sage: A = Matrix(k, 1, 1, [k.gen()])
            sage: sr.vector(A)
            [      a]
            [    a^2]
            [  a + 1]
            [a^2 + 1]
        """
        r = self.r
        c = self.c
        e = self.e
        k = self.base_ring()

        if d is None:
            return matrix(k, r*c*e, 1)
        elif d.ncols() == c and d.nrows() == r and d.base_ring() == k:
            return matrix(k, r*c*e, 1, self.phi(d).transpose().list())

    def is_vector(self, d):
        """
        Return ``True`` if ``d`` can be used as a vector for ``self``.

        EXAMPLES::

            sage: sr = mq.SR()
            sage: sr
            SR(1,1,1,4)
            sage: k = sr.base_ring()
            sage: A = Matrix(k, 1, 1, [k.gen()])
            sage: B = sr.vector(A)
            sage: sr.is_vector(A)
            False
            sage: sr.is_vector(B)
            True
        """
        return isinstance(d, Matrix) and \
               d.nrows() == self.r*self.c*self.e and \
               d.ncols() == 1 and \
               d.base_ring() == self.base_ring()

    def phi(self, l):
        r"""
        The operation `\phi` from [MR2002]_.

        Projects state arrays to their algebraic representation.

        INPUT:

        - ``l`` -- element to perform `\phi` on

        EXAMPLES::

            sage: sr = mq.SR(2, 1, 2, 4)
            sage: k = sr.base_ring()
            sage: A = matrix(k, 1, 2, [k.gen(), 0] )
            sage: sr.phi(A)
            [      a       0]
            [    a^2       0]
            [  a + 1       0]
            [a^2 + 1       0]
        """
        ret = []
        if isinstance(l, Matrix):
            for e in l.transpose().list():
                ret += [e**(2**i) for i in range(self.e)]
        else:
            for e in l:
                ret += [e**(2**i) for i in range(self.e)]
        if isinstance(l, list):
            return ret
        elif isinstance(l, tuple):
            return tuple(ret)
        elif isinstance(l, Matrix):
            return matrix(l.base_ring(), l.ncols(), l.nrows()*self.e, ret).transpose()
        else:
            raise TypeError

    def antiphi(self, l):
        r"""
        The operation `\phi^{-1}` from [MR2002]_ or the inverse of ``self.phi``.

        INPUT:

        - ``l`` -- a vector in the sense of :meth:`is_vector`

        EXAMPLES::

            sage: sr = mq.SR()
            sage: A = sr.random_state_array()
            sage: sr.antiphi(sr.phi(A)) == A
            True
        """
        if isinstance(l, Matrix):
            ret = l.transpose().list()[0:-1:self.e]
        else:
            ret = l[0:-1:self.e]

        if isinstance(l, list):
            return ret
        elif isinstance(l, tuple):
            return tuple(ret)
        elif isinstance(l, Matrix):
            return matrix(self.base_ring(), l.ncols(), l.nrows() // self.e,
                          ret).transpose()
        else:
            raise TypeError

    def shift_rows_matrix(self):
        """
        Return the ``ShiftRows`` matrix.

        EXAMPLES::

            sage: sr = mq.SR(1, 2, 2, 4)
            sage: s = sr.random_state_array()
            sage: r1 = sr.shift_rows(s)
            sage: r2 = sr.state_array( sr.shift_rows_matrix() * sr.vector(s) )
            sage: r1 == r2
            True
        """
        e = self.e
        r = self.r
        c = self.c
        k = self.base_ring()
        bs = r*c*e
        shift_rows = matrix(k, bs, bs)
        I = MatrixSpace(k, e, e)(1)
        for x in range(0, c):
            for y in range(0, r):
                _r = ((x*r)+y) * e
                _c = (((x*r)+((r+1)*y)) * e) % bs
                self._insert_matrix_into_matrix(shift_rows, I, _r, _c)

        return shift_rows

    def lin_matrix(self, length=None):
        """
        Return the ``Lin`` matrix.

        If no ``length`` is provided, the standard state space size is
        used. The key schedule calls this method with an explicit
        length argument because only ``self.r`` S-Box applications are
        performed in the key schedule.

        INPUT:

        - ``length`` -- length of state space (default: ``None``)

        EXAMPLES::

            sage: sr = mq.SR(1, 1, 1, 4)
            sage: sr.lin_matrix()
            [          a^2 + 1                 1         a^3 + a^2           a^2 + 1]
            [                a                 a                 1 a^3 + a^2 + a + 1]
            [          a^3 + a               a^2               a^2                 1]
            [                1               a^3             a + 1             a + 1]
        """
        r = self.r
        c = self.c
        e = self.e
        k = self.k

        if length is None:
            length = r*c

        lin = matrix(self.base_ring(), length*e, length*e)
        if e == 4:
            l = [k.from_integer(x) for x in (5, 1, 12, 5)]
            for k in range( 0, length ):
                for i in range(0, 4):
                    for j in range(0, 4):
                        lin[k*4+j, k*4+i] = l[(i-j) % 4] ** (2**j)
        elif e == 8:
            l = [k.from_integer(x) for x in (5, 9, 249, 37, 244, 1, 181, 143)]
            for k in range( 0, length ):
                for i in range(0, 8):
                    for j in range(0, 8):
                        lin[k*8+j, k*8+i] = l[(i-j) % 8] ** (2**j)

        return lin

    def mix_columns_matrix(self):
        """
        Return the ``MixColumns`` matrix.

        EXAMPLES::

            sage: sr = mq.SR(1, 2, 2, 4)
            sage: s = sr.random_state_array()
            sage: r1 = sr.mix_columns(s)
            sage: r2 = sr.state_array(sr.mix_columns_matrix() * sr.vector(s))
            sage: r1 == r2
            True
        """

        def D(b):
            """
            Return the `e x e` matrix `D` with `b^i` along the
            diagonal.

            EXAMPLES::

                sage: sr = mq.SR(1, 2, 1, 4)
                sage: sr.mix_columns_matrix() # indirect doctest
                [  a + 1       0       0       0       a       0       0       0]
                [      0 a^2 + 1       0       0       0     a^2       0       0]
                [      0       0       a       0       0       0   a + 1       0]
                [      0       0       0     a^2       0       0       0 a^2 + 1]
                [      a       0       0       0   a + 1       0       0       0]
                [      0     a^2       0       0       0 a^2 + 1       0       0]
                [      0       0   a + 1       0       0       0       a       0]
                [      0       0       0 a^2 + 1       0       0       0     a^2]
            """
            D = matrix(self.base_ring(), self._e, self._e)
            for i in range(self._e):
                D[i, i] = b**(2**i)
            return D

        r = self.r
        c = self.c
        e = self.e
        k = self.k
        a = k.gen()

        M = matrix(k, r*e, r*e)

        if r == 1:
            self._insert_matrix_into_matrix(M,   D(1), 0, 0)

        elif r == 2:
            self._insert_matrix_into_matrix(M, D(a+1), 0, 0)
            self._insert_matrix_into_matrix(M, D(a+1), e, e)
            self._insert_matrix_into_matrix(M,   D(a), e, 0)
            self._insert_matrix_into_matrix(M,   D(a), 0, e)

        elif r == 4:
            self._insert_matrix_into_matrix(M,   D(a),   0,   0)
            self._insert_matrix_into_matrix(M,   D(a),   e,   e)
            self._insert_matrix_into_matrix(M,   D(a), 2*e, 2*e)
            self._insert_matrix_into_matrix(M,   D(a), 3*e, 3*e)

            self._insert_matrix_into_matrix(M, D(a+1),   0,   e)
            self._insert_matrix_into_matrix(M, D(a+1),   e, 2*e)
            self._insert_matrix_into_matrix(M, D(a+1), 2*e, 3*e)
            self._insert_matrix_into_matrix(M, D(a+1), 3*e,   0)

            self._insert_matrix_into_matrix(M,   D(1),   0, 2*e)
            self._insert_matrix_into_matrix(M,   D(1),   e, 3*e)
            self._insert_matrix_into_matrix(M,   D(1), 2*e,   0)
            self._insert_matrix_into_matrix(M,   D(1), 3*e, 1*e)

            self._insert_matrix_into_matrix(M,   D(1),   0, 3*e)
            self._insert_matrix_into_matrix(M,   D(1),   e,   0)
            self._insert_matrix_into_matrix(M,   D(1), 2*e, 1*e)
            self._insert_matrix_into_matrix(M,   D(1), 3*e, 2*e)

        mix_columns = matrix(k, r*c*e, r*c*e)

        for i in range(c):
            self._insert_matrix_into_matrix(mix_columns, M, r*e*i, r*e*i)

        return mix_columns

    def inversion_polynomials(self, xi, wi, length):
        """
        Return polynomials to represent the inversion in the AES S-Box.

        INPUT:

        - ``xi`` -- output variables

        - ``wi`` -- input variables

        - ``length`` -- length of both lists

        EXAMPLES::

            sage: sr = mq.SR(1, 1, 1, 8)
            sage: R = sr.ring()
            sage: xi = Matrix(R, 8, 1, sr.vars('x', 1))
            sage: wi = Matrix(R, 8, 1, sr.vars('w', 1))
            sage: sr.inversion_polynomials(xi, wi, 8)
            [x100*w100 + 1,
            x101*w101 + 1,
            x102*w102 + 1,
            x103*w103 + 1,
            x104*w104 + 1,
            x105*w105 + 1,
            x106*w106 + 1,
            x107*w107 + 1]
        """
        return [xi[j, 0]*wi[j, 0] + 1 for j in range(length)]

    def field_polynomials(self, name, i, l=None):
        r"""
        Return list of conjugacy polynomials for a given round ``i``
        and name ``name``.

        INPUT:

        - ``name`` -- variable name
        - ``i`` -- round number
        - ``l`` -- r\*c (default: ``None``)

        EXAMPLES::

            sage: sr = mq.SR(3, 1, 1, 8)
            sage: sr.field_polynomials('x', 2)
            [x200^2 + x201,
            x201^2 + x202,
            x202^2 + x203,
            x203^2 + x204,
            x204^2 + x205,
            x205^2 + x206,
            x206^2 + x207,
            x207^2 + x200]
        """
        r = self._r
        c = self._c
        e = self._e

        if l is None:
            l = r*c

        _vars = self.vars(name, i, l, e)
        return [_vars[e*j+k]**2 - _vars[e*j+(k+1) % e] for j in range(l) for k in range(e)]


class SR_gf2(SR_generic):
    def __init__(self, n=1, r=1, c=1, e=4, star=False, **kwargs):
        r"""
        Small Scale Variants of the AES polynomial system constructor over
        `\GF{2}`. See help for SR.

        EXAMPLES::

            sage: sr = mq.SR(gf2=True)
            sage: sr
            SR(1,1,1,4)
        """
        SR_generic.__init__(self, n, r, c, e, star, **kwargs)
        self._correct_only = kwargs.get("correct_only", False)
        self._biaffine_only = kwargs.get("biaffine_only", True)

    def vector(self, d=None):
        """
        Construct a vector suitable for the algebraic representation of
        SR.

        INPUT:

        - ``d`` -- values for vector (default: ``None``)

        EXAMPLES::

            sage: sr = mq.SR(gf2=True)
            sage: sr
            SR(1,1,1,4)
            sage: k = sr.base_ring()
            sage: A = Matrix(k, 1, 1, [k.gen()])
            sage: sr.vector(A)
            [0]
            [0]
            [1]
            [0]
        """
        r = self.r
        c = self.c
        e = self.e
        k = GF(2)

        if d is None:
            return matrix(k, r*c*e, 1)
        elif isinstance(d, Matrix) and d.ncols() == c and d.nrows() == r and d.base_ring() == self.k:
            l = flatten([self.phi(x) for x in d.transpose().list()], (Vector_modn_dense,list,tuple))
            return matrix(k, r*c*e, 1, l)
        elif isinstance(d, (list, tuple)):
            if len(d) == self.r*self.c:
                l = flatten([self.phi(x) for x in d], (Vector_modn_dense,list,tuple))
                return matrix(k, r*c*e, 1, l)
            elif len(d) == self.r*self.c*self.e:
                return matrix(k, r*c*e, 1, d)
            else:
                raise TypeError
        else:
            raise TypeError

    def is_vector(self, d):
        """
        Return ``True`` if the given matrix satisfies the conditions
        for a vector as it appears in the algebraic expression of
        ``self``.

        INPUT:

        - ``d`` -- matrix

        EXAMPLES::

            sage: sr = mq.SR(gf2=True)
            sage: sr
            SR(1,1,1,4)
            sage: k = sr.base_ring()
            sage: A = Matrix(k, 1, 1, [k.gen()])
            sage: B = sr.vector(A)
            sage: sr.is_vector(A)
            False
            sage: sr.is_vector(B)
            True
        """
        return isinstance(d, Matrix) and \
               d.nrows() == self.r*self.c*self.e and \
               d.ncols() == 1 and \
               d.base_ring() == GF(2)

    def phi(self, l, diffusion_matrix=False):
        r"""
        The operation `\phi` from [MR2002]_.

        Given a list/matrix of elements in `\GF{2^e}`, return a
        matching list/matrix of elements in `\GF{2}`.

        INPUT:

        - ``l`` -- element to perform `\phi` on
        - ``diffusion_matrix`` -- if ``True``, the given matrix ``l`` is
          transformed to a matrix which performs the same operation
          over `\GF{2}` as ``l`` over `\GF{2^n}` (default: ``False``).

        EXAMPLES::

            sage: sr = mq.SR(2, 1, 2, 4, gf2=True)
            sage: k = sr.base_ring()
            sage: A = matrix(k, 1, 2, [k.gen(), 0] )
            sage: sr.phi(A)                                                             # needs sage.libs.gap
            [0 0]
            [0 0]
            [1 0]
            [0 0]
        """
        ret = []
        r, c, e = self.r, self.c, self.e

        # handle diffusion layer matrices first
        if isinstance(l, Matrix) and diffusion_matrix and \
           l.nrows() == r*c and l.ncols() == r*c and \
           l.base_ring() == self.k:
            B = matrix(GF(2), r*c*e, r*c*e)
            for x in range(r*c):
                for y in range(r*c):
                    T = self._mul_matrix(l[x, y])
                    self._insert_matrix_into_matrix(B, T, x*e, y*e)
            return B

        # ground field elements
        if l in self.k:
            return list(reversed(l._vector_()))

        # remaining matrices
        if isinstance(l, Matrix):
            for x in l.transpose().list():
                ret += list(reversed(x._vector_()))
        # or lists
        else:
            for x in l:
                ret += list(reversed(x._vector_()))

        if isinstance(l, list):
            return ret
        elif isinstance(l, tuple):
            return tuple(ret)
        elif isinstance(l, Matrix):
            return matrix(GF(2), l.ncols(), l.nrows()*self.e, ret).transpose()
        else:
            raise TypeError

    def antiphi(self, l):
        r"""
        The operation `\phi^{-1}` from [MR2002]_ or the inverse of ``self.phi``.

        INPUT:

        - ``l`` -- a vector in the sense of ``self.is_vector``

        EXAMPLES::

            sage: sr = mq.SR(gf2=True)
            sage: A = sr.random_state_array()
            sage: sr.antiphi(sr.phi(A)) == A                                            # needs sage.libs.gap
            True
        """
        e = self.e
        V = self.k.vector_space(map=False)

        if isinstance(l, Matrix):
            l2 = l.transpose().list()
        else:
            l2 = l

        ret = []
        for i in range(0, len(l2), e):
            ret.append( self.k(V(list(reversed(l2[i:i+e])))) )

        if isinstance(l, list):
            return ret
        elif isinstance(l, tuple):
            return tuple(ret)
        elif isinstance(l, Matrix):
            return matrix(self.base_ring(), self.r * self.c, 1, ret)
        else:
            raise TypeError

    def shift_rows_matrix(self):
        """
        Return the ``ShiftRows`` matrix.

        EXAMPLES::

            sage: sr = mq.SR(1, 2, 2, 4, gf2=True)
            sage: s = sr.random_state_array()
            sage: r1 = sr.shift_rows(s)
            sage: r2 = sr.state_array( sr.shift_rows_matrix() * sr.vector(s) )
            sage: r1 == r2
            True
        """
        r = self.r
        c = self.c
        k = self.k
        bs = r*c
        shift_rows = matrix(k, r*c, r*c)
        for x in range(0, c):
            for y in range(0, r):
                _r = ((x*r)+y)
                _c = ((x*r)+((r+1)*y)) % bs
                shift_rows[_r, _c] = 1
        return self.phi(shift_rows, diffusion_matrix=True)

    def mix_columns_matrix(self):
        """
        Return the ``MixColumns`` matrix.

        EXAMPLES::

            sage: sr = mq.SR(1, 2, 2, 4, gf2=True)
            sage: s = sr.random_state_array()
            sage: r1 = sr.mix_columns(s)
            sage: r2 = sr.state_array(sr.mix_columns_matrix() * sr.vector(s))
            sage: r1 == r2
            True
        """
        r = self.r
        c = self.c
        k = self.k
        a = k.gen()

        if r == 1:
            M = matrix(k, r, r, 1)

        elif r == 2:
            M = matrix(k, r, r, [a+1, a, a, a+1])

        elif r == 4:
            M = matrix(k, r, [a, a+1, 1, 1,
                             1, a, a+1, 1,
                             1, 1, a, a+1,
                             a+1, 1, 1, a])

        mix_columns = matrix(k, r*c, r*c)

        for i in range(c):
            self._insert_matrix_into_matrix(mix_columns, M, r*i, r*i)

        return self.phi(mix_columns, diffusion_matrix=True)

    def lin_matrix(self, length=None):
        """
        Return the ``Lin`` matrix.

        If no ``length`` is provided, the standard state space size is
        used. The key schedule calls this method with an explicit
        length argument because only ``self.r`` S-Box applications are
        performed in the key schedule.

        INPUT:

        - ``length`` -- length of state space (default: ``None``)

        EXAMPLES::

            sage: sr = mq.SR(1, 1, 1, 4, gf2=True)
            sage: sr.lin_matrix()
            [1 0 1 1]
            [1 1 0 1]
            [1 1 1 0]
            [0 1 1 1]
        """
        r, c, e = self.r, self.c, self.e

        if length is None:
            length = r*c

        if e == 8:
            Z = matrix(GF(2), 8, 8, [1, 0, 0, 0, 1, 1, 1, 1,
                                  1, 1, 0, 0, 0, 1, 1, 1,
                                  1, 1, 1, 0, 0, 0, 1, 1,
                                  1, 1, 1, 1, 0, 0, 0, 1,
                                  1, 1, 1, 1, 1, 0, 0, 0,
                                  0, 1, 1, 1, 1, 1, 0, 0,
                                  0, 0, 1, 1, 1, 1, 1, 0,
                                  0, 0, 0, 1, 1, 1, 1, 1])
        else:
            Z = matrix(GF(2), 4, 4, [1, 1, 1, 0,
                                  0, 1, 1, 1,
                                  1, 0, 1, 1,
                                  1, 1, 0, 1])

        Z = Z.transpose() # account for endianess mismatch

        lin = matrix(GF(2), length*e, length*e)

        for i in range(length):
            self._insert_matrix_into_matrix(lin, Z, i*e, i*e)
        return lin

    def _mul_matrix(self, x):
        r"""
        Given an element `x` in self.base_ring(), return a matrix
        which performs the same operation on a when interpreted over
        `\GF{2^e}` as `x` over `\GF{2^e}`.

        INPUT:

        - ``x`` -- an element in self.base_ring()

        EXAMPLES::

            sage: sr = mq.SR(gf2=True)
            sage: a = sr.k.gen()
            sage: A = sr._mul_matrix(a^2+1)
            sage: sr.antiphi( A *  sr.vector([a+1]) )
            [a^3 + a^2 + a + 1]

        ::

            sage: (a^2 + 1)*(a+1)
            a^3 + a^2 + a + 1
        """
        k = self.k
        e = self.e
        a = k.gen()

        columns = [list(reversed((x * a**i)._vector_()))
                   for i in reversed(range(e))]
        return matrix(GF(2), e, e, columns).transpose()

    def _square_matrix(self):
        """
        Return a matrix of dimension self.e x self.e which performs the
        squaring operation over `GF(2^n)` on vectors of length e.

        EXAMPLES::

            sage: sr = mq.SR(gf2=True)
            sage: a = sr.k.gen()
            sage: S = sr._square_matrix()
            sage: sr.antiphi( S *  sr.vector([a^3+1]) )
            [a^3 + a^2 + 1]

        ::

            sage: (a^3 + 1)^2
            a^3 + a^2 + 1
        """
        a = self.k.gen()
        e = self.e

        columns = []
        for i in reversed(range(e)):
            columns.append( list(reversed(((a**i)**2)._vector_())) )
        return matrix(GF(2), e , e, columns).transpose()

    def inversion_polynomials_single_sbox(self, x=None, w=None, biaffine_only=None, correct_only=None):
        """
        Return inversion polynomials of a single S-Box.

        INPUT:

        - ``xi`` -- output variables
        - ``wi`` -- input variables
        - ``length`` -- length of both lists

        EXAMPLES::

            sage: sr = mq.SR(1, 1, 1, 8, gf2=True)
            sage: len(sr.inversion_polynomials_single_sbox())
            24
            sage: len(sr.inversion_polynomials_single_sbox(correct_only=True))
            23
            sage: len(sr.inversion_polynomials_single_sbox(biaffine_only=False))
            40
            sage: len(sr.inversion_polynomials_single_sbox(biaffine_only=False, correct_only=True))
            39


            sage: sr = mq.SR(1, 1, 1, 8, gf2=True)
            sage: l0 = sr.inversion_polynomials_single_sbox(); len(l0)
            24
            sage: l1 = sr.inversion_polynomials_single_sbox(correct_only=True); len(l1)
            23
            sage: l2 = sr.inversion_polynomials_single_sbox(biaffine_only=False); len(l2)
            40
            sage: l3 = sr.inversion_polynomials_single_sbox(biaffine_only=False, correct_only=True); len(l3)
            39

            sage: set(l0) == set(sr._inversion_polynomials_single_sbox())
            True
            sage: set(l1) == set(sr._inversion_polynomials_single_sbox(correct_only=True))
            True
            sage: set(l2) == set(sr._inversion_polynomials_single_sbox(biaffine_only=False))
            True
            sage: set(l3) == set(sr._inversion_polynomials_single_sbox(biaffine_only=False, correct_only=True))
            True

            sage: sr = mq.SR(1, 1, 1, 4, gf2=True)
            sage: l0 = sr.inversion_polynomials_single_sbox(); len(l0)
            12
            sage: l1 = sr.inversion_polynomials_single_sbox(correct_only=True); len(l1)
            11
            sage: l2 = sr.inversion_polynomials_single_sbox(biaffine_only=False); len(l2)
            20
            sage: l3 = sr.inversion_polynomials_single_sbox(biaffine_only=False, correct_only=True); len(l3)
            19

            sage: set(l0) == set(sr._inversion_polynomials_single_sbox())
            True
            sage: set(l1) == set(sr._inversion_polynomials_single_sbox(correct_only=True))
            True
            sage: set(l2) == set(sr._inversion_polynomials_single_sbox(biaffine_only=False))
            True
            sage: set(l3) == set(sr._inversion_polynomials_single_sbox(biaffine_only=False, correct_only=True))
            True
        """
        e = self.e

        if biaffine_only is None:
            biaffine_only = self._biaffine_only
        if correct_only is None:
            correct_only = self._correct_only

        if x is None and w is None:
            # make sure it prints like in the book.
            names = ["w%d" % i for i in reversed(range(e))] + ["x%d" % i for i in reversed(range(e))]
            P = PolynomialRing(GF(2), e*2, names, order='lex')
            x = P.gens()[e:]
            w = P.gens()[:e]
        else:
            if isinstance(x, (tuple, list)):
                P = x[0].parent()
            elif isinstance(x, Matrix):
                P = x.base_ring()
            else:
                raise TypeError("x not understood")

            if isinstance(x, Matrix):
                x = x.column(0).list()
            if isinstance(w, Matrix):
                w = w.column(0).list()

        if e == 4:
            w3,w2,w1,w0 = w
            x3,x2,x1,x0 = x

            l = [w3*x3 + w3*x0 + w2*x1 + w1*x2 + w0*x3,
                 w3*x3 + w3*x2 + w2*x3 + w2*x0 + w1*x1 + w0*x2,
                 w3*x2 + w3*x1 + w2*x3 + w2*x2 + w1*x3 + w1*x0 + w0*x1,
                 w3*x3 + w3*x2 + w3*x0 + w2*x2 + w1*x3 + w1*x1 + w0*x3 + x3,
                 w3*x1 + w2*x3 + w2*x2 + w2*x0 + w1*x2 + w0*x3 + w0*x1 + x2,
                 w3*x3 + w3*x2 + w3*x1 + w2*x1 + w1*x3 + w1*x2 + w1*x0 + w0*x2 + x1,
                 w3*x2 + w2*x3 + w2*x1 + w1*x3 + w0*x2 + w0*x0 + x0,
                 w3*x3 + w3*x1 + w3*x0 + w3 + w2*x3 + w2*x2 + w1*x1 + w0*x3,
                 w3*x2 + w3*x0 + w2*x2 + w2*x1 + w2 + w1*x3 + w1*x0 + w0*x2,
                 w3*x3 + w3*x1 + w2*x3 + w2*x1 + w2*x0 + w1*x3 + w1*x2 + w1 + w0*x1,
                 w3*x2 + w3*x1 + w2*x3 + w2*x0 + w1*x2 + w0*x0 + w0]

            if not correct_only:
                l.append(w3*x1 + w2*x2 + w1*x3 + w0*x0 + 1)

            if not biaffine_only:
                l.extend([w3*x2 + w3*x1 + w3*x0 + w2*x3 + w2*x1 + w1*x3 + w1*x2 + w0*x3 + x3**2 + x3*x2 + x3*x1 + x2**2 + x1**2,
                          w3*x2 + w2*x2 + w2*x1 + w2*x0 + w1*x3 + w1*x1 + w0*x3 + w0*x2 + x3*x2 + x3*x1 + x3*x0 + x2**2 + x2*x1 + x2*x0 + x1*x0,
                          w3*x2 + w3*x1 + w2*x2 + w1*x2 + w1*x1 + w1*x0 + w0*x3 + w0*x1 + x3**2 + x3*x2 + x2*x0 + x1*x0,
                          w3*x3 + w3*x1 + w2*x3 + w2*x2 + w1*x3 + w0*x3 + w0*x2 + w0*x1 + w0*x0 + x3*x1 + x2*x1 + x2*x0 + x0**2,
                          w3**2 + w3*w2 + w3*w1 + w3*x2 + w3*x1 + w3*x0 + w2**2 + w2*x3 + w2*x1 + w1**2 + w1*x3 + w1*x2 + w0*x3,
                          w3*w2 + w3*w1 + w3*w0 + w3*x1 + w3*x0 + w2**2 + w2*w1 + w2*w0 + w2*x3 + w2*x2 + w2*x0 + w1*w0 + w1*x2 + w1*x1 + w0*x2,
                          w3**2 + w3*w2 + w3*x0 + w2*w0 + w2*x3 + w2*x2 + w2*x1 + w1*w0 + w1*x3 + w1*x1 + w1*x0 + w0*x1,
                          w3*w1 + w3*x3 + w3*x2 + w3*x1 + w3*x0 + w2*w1 + w2*w0 + w2*x2 + w2*x0 + w1*x3 + w1*x0 + w0**2 + w0*x0])
            return l

        else:
            w7,w6,w5,w4,w3,w2,w1,w0 = w
            x7,x6,x5,x4,x3,x2,x1,x0 = x

            l = [w7*x7 + w7*x5 + w7*x4 + w7*x0 + w6*x6 + w6*x5 + w6*x1 + w5*x7 + w5*x6 + w5*x2 + w4*x7 + w4*x3 + w3*x4 + w2*x5 + w1*x6 + w0*x7,
                 w7*x6 + w7*x4 + w7*x3 + w6*x7 + w6*x5 + w6*x4 + w6*x0 + w5*x6 + w5*x5 + w5*x1 + w4*x7 + w4*x6 + w4*x2 + w3*x7 + w3*x3 + w2*x4 + w1*x5 + w0*x6,
                 w7*x5 + w7*x3 + w7*x2 + w6*x6 + w6*x4 + w6*x3 + w5*x7 + w5*x5 + w5*x4 + w5*x0 + w4*x6 + w4*x5 + w4*x1 + w3*x7 + w3*x6 + w3*x2 + w2*x7 + w2*x3 + w1*x4 + w0*x5,
                 w7*x7 + w7*x4 + w7*x2 + w7*x1 + w6*x5 + w6*x3 + w6*x2 + w5*x6 + w5*x4 + w5*x3 + w4*x7 + w4*x5 + w4*x4 + w4*x0 + w3*x6 + w3*x5 + w3*x1 + w2*x7 + w2*x6 + w2*x2 + w1*x7 + w1*x3 + w0*x4,
                 w7*x7 + w7*x6 + w7*x5 + w7*x4 + w7*x3 + w7*x1 + w6*x7 + w6*x6 + w6*x5 + w6*x4 + w6*x2 + w5*x7 + w5*x6 + w5*x5 + w5*x3 + w4*x7 + w4*x6 + w4*x4 + w3*x7 + w3*x5 + w3*x0 + w2*x6 + w2*x1
                     + w1*x7 + w1*x2 + w0*x3,
                 w7*x6 + w7*x3 + w7*x2 + w6*x7 + w6*x4 + w6*x3 + w5*x5 + w5*x4 + w4*x6 + w4*x5 + w3*x7 + w3*x6 + w2*x7 + w2*x0 + w1*x1 + w0*x2,
                 w7*x7 + w7*x5 + w7*x2 + w7*x1 + w6*x6 + w6*x3 + w6*x2 + w5*x7 + w5*x4 + w5*x3 + w4*x5 + w4*x4 + w3*x6 + w3*x5 + w2*x7 + w2*x6 + w1*x7 + w1*x0 + w0*x1,
                 w7*x6 + w7*x5 + w7*x2 + w7*x0 + w6*x7 + w6*x4 + w6*x3 + w5*x7 + w5*x6 + w5*x3 + w5*x1 + w4*x5 + w4*x4 + w3*x7 + w3*x4 + w3*x2 + w2*x6 + w2*x5 + w1*x5 + w1*x3 + w0*x7 + w0*x6 + x7,
                 w7*x6 + w7*x3 + w7*x2 + w6*x6 + w6*x5 + w6*x2 + w6*x0 + w5*x7 + w5*x4 + w5*x3 + w4*x7 + w4*x6 + w4*x3 + w4*x1 + w3*x5 + w3*x4 + w2*x7 + w2*x4 + w2*x2 + w1*x6 + w1*x5 + w0*x5 + w0*x3
                     + x6,
                 w7*x7 + w7*x5 + w7*x4 + w7*x1 + w6*x6 + w6*x3 + w6*x2 + w5*x6 + w5*x5 + w5*x2 + w5*x0 + w4*x7 + w4*x4 + w4*x3 + w3*x7 + w3*x6 + w3*x3 + w3*x1 + w2*x5 + w2*x4 + w1*x7 + w1*x4 + w1*x2
                     + w0*x6 + w0*x5 + x5,
                 w7*x7 + w7*x5 + w7*x2 + w7*x1 + w6*x7 + w6*x5 + w6*x4 + w6*x1 + w5*x6 + w5*x3 + w5*x2 + w4*x6 + w4*x5 + w4*x2 + w4*x0 + w3*x7 + w3*x4 + w3*x3 + w2*x7 + w2*x6 + w2*x3 + w2*x1 + w1*x5
                     + w1*x4 + w0*x7 + w0*x4 + w0*x2 + x4,
                 w7*x5 + w7*x4 + w7*x3 + w7*x2 + w6*x5 + w6*x4 + w6*x3 + w6*x2 + w6*x1 + w5*x6 + w5*x5 + w5*x4 + w5*x3 + w4*x6 + w4*x5 + w4*x4 + w4*x3 + w4*x2 + w3*x7 + w3*x6 + w3*x5 + w3*x4 + w3*x0
                     + w2*x7 + w2*x6 + w2*x5 + w2*x4 + w2*x3 + w1*x7 + w1*x6 + w1*x5 + w1*x1 + w0*x7 + w0*x6 + w0*x5 + w0*x4 + x3,
                 w7*x7 + w7*x6 + w7*x5 + w7*x4 + w7*x3 + w7*x1 + w6*x7 + w6*x5 + w6*x2 + w5*x7 + w5*x6 + w5*x5 + w5*x4 + w5*x2 + w4*x6 + w4*x3 + w3*x7 + w3*x6 + w3*x5 + w3*x3 + w2*x7 + w2*x4 + w2*x0
                     + w1*x7 + w1*x6 + w1*x4 + w0*x5 + w0*x1 + x2,
                 w7*x6 + w7*x4 + w7*x1 + w6*x7 + w6*x6 + w6*x5 + w6*x4 + w6*x3 + w6*x1 + w5*x7 + w5*x5 + w5*x2 + w4*x7 + w4*x6 + w4*x5 + w4*x4 + w4*x2 + w3*x6 + w3*x3 + w2*x7 + w2*x6 + w2*x5 + w2*x3
                     + w1*x7 + w1*x4 + w1*x0 + w0*x7 + w0*x6 + w0*x4 + x1,
                 w7*x7 + w7*x4 + w7*x3 + w6*x7 + w6*x6 + w6*x3 + w6*x1 + w5*x5 + w5*x4 + w4*x7 + w4*x4 + w4*x2 + w3*x6 + w3*x5 + w2*x5 + w2*x3 + w1*x7 + w1*x6 + w0*x6 + w0*x4 + w0*x0 + x0,
                 w7*x6 + w7*x5 + w7*x3 + w7*x0 + w7 + w6*x7 + w6*x5 + w6*x2 + w6*x0 + w5*x7 + w5*x4 + w5*x2 + w5*x1 + w4*x6 + w4*x4 + w4*x3 + w3*x6 + w3*x5 + w3*x1 + w2*x7 + w2*x3 + w1*x5 + w0*x7,
                 w7*x5 + w7*x4 + w7*x2 + w6*x7 + w6*x6 + w6*x4 + w6*x1 + w6 + w5*x6 + w5*x3 + w5*x1 + w5*x0 + w4*x5 + w4*x3 + w4*x2 + w3*x7 + w3*x5 + w3*x4 + w3*x0 + w2*x7 + w2*x6 + w2*x2 + w1*x4
                     + w0*x6,
                 w7*x7 + w7*x4 + w7*x3 + w7*x1 + w6*x6 + w6*x5 + w6*x3 + w6*x0 + w5*x7 + w5*x5 + w5*x2 + w5*x0 + w5 + w4*x7 + w4*x4 + w4*x2 + w4*x1 + w3*x6 + w3*x4 + w3*x3 + w2*x6 + w2*x5 + w2*x1
                     + w1*x7 + w1*x3 + w0*x5,
                 w7*x7 + w7*x6 + w7*x3 + w7*x2 + w7*x0 + w6*x5 + w6*x4 + w6*x2 + w5*x7 + w5*x6 + w5*x4 + w5*x1 + w4*x6 + w4*x3 + w4*x1 + w4*x0 + w4 + w3*x5 + w3*x3 + w3*x2 + w2*x7 + w2*x5 + w2*x4
                     + w2*x0 + w1*x7 + w1*x6 + w1*x2 + w0*x4,
                 w7*x3 + w7*x2 + w7*x1 + w7*x0 + w6*x5 + w6*x4 + w6*x3 + w6*x2 + w6*x1 + w6*x0 + w5*x7 + w5*x6 + w5*x5 + w5*x4 + w5*x3 + w5*x2 + w5*x1 + w5*x0 + w4*x7 + w4*x6 + w4*x5 + w4*x4
                     + w4*x3 + w4*x2 + w4*x0 + w3*x7 + w3*x6 + w3*x5 + w3*x4 + w3*x2 + w3 + w2*x7 + w2*x6 + w2*x4 + w1*x6 + w1*x1 + w0*x3,
                 w7*x7 + w7*x6 + w7*x5 + w7*x3 + w7*x2 + w7*x1 + w6*x7 + w6*x5 + w6*x4 + w6*x3 + w6*x1 + w5*x7 + w5*x6 + w5*x5 + w5*x3 + w5*x0 + w4*x7 + w4*x5 + w4*x2 + w4*x1 + w3*x7 + w3*x4
                     + w3*x3 + w2*x6 + w2*x5 + w2 + w1*x7 + w1*x0 + w0*x2,
                 w7*x6 + w7*x5 + w7*x4 + w7*x2 + w7*x1 + w7*x0 + w6*x7 + w6*x6 + w6*x4 + w6*x3 + w6*x2 + w6*x0 + w5*x6 + w5*x5 + w5*x4 + w5*x2 + w4*x7 + w4*x6 + w4*x4 + w4*x1 + w4*x0 + w3*x6
                     + w3*x3 + w3*x2 + w2*x5 + w2*x4 + w1*x7 + w1*x6 + w1 + w0*x1,
                 w7*x7 + w7*x6 + w7*x4 + w7*x1 + w6*x6 + w6*x3 + w6*x1 + w6*x0 + w5*x5 + w5*x3 + w5*x2 + w4*x7 + w4*x5 + w4*x4 + w4*x0 + w3*x7 + w3*x6 + w3*x2 + w2*x4 + w1*x6 + w0*x0 + w0]

            if not correct_only:
                l.append(w7*x6 + w7*x5 + w7*x1 + w6*x7 + w6*x6 + w6*x2 + w5*x7 + w5*x3 + w4*x4 + w3*x5 + w2*x6 + w1*x7 + w0*x0 + 1)

            if not biaffine_only:
                l.extend([w7**2 + w7*w6 + w7*w3 + w7*w1 + w7*x7 + w7*x6 + w7*x5 + w7*x2 + w7*x1 + w7*x0 + w6**2 + w6*w0 + w6*x6 + w6*x5 + w6*x4 + w6*x3 + w6*x1 + w6*x0 + w5**2 + w5*w4 + w5*w3
                              + w5*w2 + w5*x7 + w5*x5 + w5*x4 + w5*x1 + w5*x0 + w4**2 + w4*w2 + w4*w0 + w4*x5 + w4*x4 + w4*x2 + w3*w2 + w3*x6 + w3*x3 + w3*x1 + w3*x0 + w2*x7 + w2*x5 + w2*x4
                              + w2*x0 + w1*x4 + w0**2 + w0*x0,
                          w7*x6 + w7*x4 + w7*x1 + w6*x7 + w6*x6 + w6*x5 + w6*x2 + w5*x7 + w5*x6 + w5*x5 + w5*x4 + w5*x3 + w5*x1 + w4*x5 + w4*x4 + w4*x3 + w4*x1 + w4*x0 + w3*x7 + w3*x5 + w3*x2
                              + w2*x7 + w2*x6 + w2*x3 + w1*x7 + w1*x6 + w1*x5 + w1*x4 + w1*x2 + w0*x6 + w0*x5 + w0*x4 + w0*x2 + w0*x1 + x7**2 + x7*x6 + x7*x5 + x7*x3 + x7*x1 + x7*x0 + x6*x2
                              + x6*x1 + x5*x4 + x5*x3 + x5*x2 + x5*x1 + x4*x3 + x4*x2 + x4*x1 + x3**2 + x3*x2 + x2*x1 + x2*x0,
                          w7*x5 + w7*x4 + w7*x3 + w7*x1 + w7*x0 + w6*x7 + w6*x5 + w6*x2 + w5*x7 + w5*x6 + w5*x3 + w4*x7 + w4*x6 + w4*x5 + w4*x4 + w4*x2 + w3*x6 + w3*x5 + w3*x4 + w3*x2 + w3*x1
                              + w2*x6 + w2*x3 + w1*x7 + w1*x4 + w0*x7 + w0*x6 + w0*x5 + w0*x3 + x7*x3 + x7*x2 + x6*x5 + x6*x4 + x6*x3 + x6*x2 + x6*x0 + x5*x4 + x5*x3 + x5*x2 + x4**2 + x4*x3
                              + x3*x2 + x3*x1,
                          w7*w3 + w7*w2 + w7*x6 + w7*x5 + w7*x4 + w7*x1 + w7*x0 + w6*w5 + w6*w4 + w6*w3 + w6*w2 + w6*w0 + w6*x5 + w6*x4 + w6*x3 + w6*x2 + w6*x0 + w5*w4 + w5*w3 + w5*w2 + w5*x7
                              + w5*x6 + w5*x4 + w5*x3 + w5*x0 + w4**2 + w4*w3 + w4*x7 + w4*x4 + w4*x3 + w4*x1 + w3*w2 + w3*w1 + w3*x7 + w3*x5 + w3*x2 + w3*x0 + w2*x6 + w2*x4 + w2*x3 + w1*x7
                              + w1*x3 + w0*x7,
                          w7*x5 + w7*x2 + w7*x1 + w6*x7 + w6*x6 + w6*x5 + w6*x4 + w6*x2 + w6*x1 + w5*x5 + w5*x3 + w5*x2 + w4*x3 + w4*x2 + w4*x1 + w3*x6 + w3*x3 + w3*x2 + w3*x0 + w2*x7 + w2*x6
                              + w2*x5 + w2*x3 + w2*x2 + w1*x6 + w1*x4 + w1*x3 + w0*x4 + w0*x3 + w0*x2 + x7*x5 + x7*x4 + x7*x1 + x7*x0 + x6*x0 + x5**2 + x5*x2 + x5*x1 + x5*x0 + x4**2 + x4*x0
                              + x3*x2 + x3*x0 + x1**2,
                          w7*w6 + w7*w5 + w7*w4 + w7*w3 + w7*x7 + w7*x5 + w7*x4 + w7*x3 + w7*x0 + w6**2 + w6*w5 + w6*w4 + w6*w2 + w6*w1 + w6*w0 + w6*x7 + w6*x4 + w6*x3 + w6*x2 + w6*x1 + w5*w4
                              + w5*w1 + w5*w0 + w5*x7 + w5*x6 + w5*x5 + w5*x3 + w5*x2 + w4*w2 + w4*w1 + w4*x7 + w4*x6 + w4*x3 + w4*x2 + w4*x0 + w3*w0 + w3*x7 + w3*x6 + w3*x4 + w3*x1 + w2**2
                              + w2*x5 + w2*x3 + w2*x2 + w1*x7 + w1*x6 + w1*x2 + w0*x6,
                          w7*w5 + w7*w4 + w7*w1 + w7*w0 + w7*x6 + w7*x2 + w6*w0 + w6*x6 + w6*x3 + w6*x2 + w6*x1 + w5**2 + w5*w2 + w5*w1 + w5*w0 + w5*x7 + w5*x6 + w5*x5 + w5*x2 + w4**2 + w4*w0
                              + w4*x6 + w4*x1 + w4*x0 + w3*w2 + w3*w0 + w3*x5 + w3*x4 + w3*x3 + w3*x2 + w3*x1 + w3*x0 + w2*x7 + w2*x6 + w2*x5 + w2*x4 + w2*x3 + w2*x2 + w2*x0 + w1**2 + w1*x7
                              + w1*x6 + w1*x4 + w0*x3,
                          w7*x7 + w7*x6 + w7*x5 + w7*x2 + w6*x7 + w6*x6 + w6*x5 + w6*x4 + w6*x3 + w6*x1 + w5*x5 + w5*x4 + w5*x3 + w5*x1 + w5*x0 + w4*x7 + w4*x5 + w4*x2 + w3*x7 + w3*x6 + w3*x3
                              + w2*x7 + w2*x6 + w2*x5 + w2*x4 + w2*x2 + w1*x6 + w1*x5 + w1*x4 + w1*x2 + w1*x1 + w0*x6 + w0*x3 + x7**2 + x7*x5 + x7*x3 + x6**2 + x6*x5 + x6*x2 + x6*x0 + x5**2
                              + x4**2 + x4*x3 + x4*x2 + x4*x1 + x3**2 + x3*x1 + x2*x1,
                          w7**2 + w7*w6 + w7*w5 + w7*w3 + w7*w1 + w7*w0 + w7*x6 + w7*x5 + w7*x3 + w7*x2 + w7*x1 + w6*w2 + w6*w1 + w6*x7 + w6*x6 + w6*x5 + w6*x2 + w6*x1 + w6*x0 + w5*w4 + w5*w3
                              + w5*w2 + w5*w1 + w5*x6 + w5*x5 + w5*x4 + w5*x3 + w5*x1 + w5*x0 + w4*w3 + w4*w2 + w4*w1 + w4*x7 + w4*x5 + w4*x4 + w4*x1 + w4*x0 + w3**2 + w3*w2 + w3*x5 + w3*x4
                              + w3*x2 + w2*w1 + w2*w0 + w2*x6 + w2*x3 + w2*x1 + w2*x0 + w1*x7 + w1*x5 + w1*x4 + w1*x0 + w0*x4,
                          w7*x7 + w7*x5 + w7*x2 + w6*x7 + w6*x6 + w6*x3 + w5*x7 + w5*x6 + w5*x5 + w5*x4 + w5*x2 + w4*x6 + w4*x5 + w4*x4 + w4*x2 + w4*x1 + w3*x6 + w3*x3 + w2*x7 + w2*x4 + w1*x7
                              + w1*x6 + w1*x5 + w1*x3 + w0*x7 + w0*x6 + w0*x5 + w0*x3 + w0*x2 + w0*x0 + x7**2 + x7*x6 + x7*x3 + x7*x1 + x6**2 + x6*x0 + x5**2 + x5*x4 + x5*x3 + x5*x2 + x4**2
                              + x4*x2 + x4*x0 + x3*x2 + x0**2,
                          w7*x7 + w7*x6 + w7*x5 + w7*x4 + w7*x3 + w7*x1 + w6*x5 + w6*x4 + w6*x3 + w6*x1 + w6*x0 + w5*x7 + w5*x5 + w5*x2 + w4*x7 + w4*x6 + w4*x3 + w3*x7 + w3*x6 + w3*x5 + w3*x4
                              + w3*x2 + w2*x6 + w2*x5 + w2*x4 + w2*x2 + w2*x1 + w1*x6 + w1*x3 + w0*x7 + w0*x4 + x7*x6 + x7*x5 + x7*x4 + x7*x3 + x6**2 + x6*x5 + x6*x4 + x6*x2 + x6*x1 + x6*x0
                              + x5*x4 + x5*x1 + x5*x0 + x4*x2 + x4*x1 + x3*x0 + x2**2,
                          w7*x5 + w7*x4 + w7*x3 + w7*x2 + w6*x7 + w6*x1 + w5*x5 + w5*x4 + w5*x3 + w5*x2 + w5*x1 + w4*x7 + w4*x6 + w4*x4 + w4*x3 + w3*x6 + w3*x5 + w3*x4 + w3*x3 + w2*x2 + w2*x0
                              + w1*x6 + w1*x5 + w1*x4 + w1*x3 + w1*x2 + w0*x7 + w0*x5 + w0*x4 + x7**2 + x7*x4 + x7*x2 + x6*x4 + x6*x3 + x6*x2 + x6*x1 + x5**2 + x5*x4 + x5*x3 + x5*x2 + x5*x0
                              + x4*x3 + x4*x2 + x4*x1 + x3**2 + x2*x0 + x1*x0,
                          w7*x6 + w7*x5 + w7*x3 + w7*x2 + w6*x5 + w6*x4 + w6*x3 + w6*x2 + w5*x7 + w5*x1 + w4*x5 + w4*x4 + w4*x3 + w4*x2 + w4*x1 + w3*x7 + w3*x6 + w3*x4 + w3*x3 + w2*x6 + w2*x5
                              + w2*x4 + w2*x3 + w1*x2 + w1*x0 + w0*x6 + w0*x5 + w0*x4 + w0*x3 + w0*x2 + x7*x5 + x7*x2 + x7*x0 + x6**2 + x6*x5 + x6*x2 + x6*x1 + x6*x0 + x5**2 + x5*x4 + x4**2
                              + x4*x2 + x4*x1 + x4*x0 + x3**2 + x3*x2 + x1*x0,
                          w7**2 + w7*w5 + w7*w3 + w7*x7 + w7*x6 + w7*x4 + w7*x3 + w7*x2 + w6**2 + w6*w5 + w6*w2 + w6*w0 + w6*x7 + w6*x6 + w6*x3 + w6*x2 + w6*x1 + w6*x0 + w5**2 + w5*x7 + w5*x6
                              + w5*x5 + w5*x4 + w5*x2 + w5*x1 + w4**2 + w4*w3 + w4*w2 + w4*w1 + w4*x6 + w4*x5 + w4*x2 + w4*x1 + w3**2 + w3*w1 + w3*x6 + w3*x5 + w3*x3 + w3*x0 + w2*w1 + w2*x7
                              + w2*x4 + w2*x2 + w2*x1 + w1*x6 + w1*x5 + w1*x1 + w0*x5,
                          w7*w5 + w7*w2 + w7*w0 + w7*x5 + w7*x3 + w6**2 + w6*w5 + w6*w2 + w6*w1 + w6*w0 + w6*x7 + w6*x3 + w6*x2 + w6*x0 + w5**2 + w5*w4 + w5*x7 + w5*x6 + w5*x4 + w5*x2 + w5*x0
                              + w4**2 + w4*w2 + w4*w1 + w4*w0 + w4*x6 + w4*x4 + w4*x3 + w4*x2 + w4*x0 + w3**2 + w3*w2 + w3*x7 + w3*x6 + w3*x4 + w3*x3 + w3*x2 + w3*x0 + w2*x7 + w2*x6 + w2*x4
                              + w2*x1 + w2*x0 + w1*w0 + w1*x5 + w1*x4 + w0*x1,
                          w7**2 + w7*w4 + w7*w2 + w7*x6 + w7*x4 + w7*x0 + w6*w4 + w6*w3 + w6*w2 + w6*w1 + w6*x4 + w6*x3 + w6*x1 + w5**2 + w5*w4 + w5*w3 + w5*w2 + w5*w0 + w5*x7 + w5*x5 + w5*x3
                              + w5*x1 + w5*x0 + w4*w3 + w4*w2 + w4*w1 + w4*x7 + w4*x5 + w4*x4 + w4*x3 + w4*x1 + w4*x0 + w3**2 + w3*x7 + w3*x5 + w3*x4 + w3*x3 + w3*x1 + w2*w0 + w2*x7 + w2*x5
                              + w2*x2 + w2*x1 + w1*w0 + w1*x6 + w1*x5 + w0*x2])

        return l

    def _inversion_polynomials_single_sbox(self, x=None, w=None, biaffine_only=None, correct_only=None):
        """
        Generate inversion polynomials of a single S-box.

        INPUT:

        - ``x`` -- output variables (default: ``None``)
        - ``w`` -- input variables  (default: ``None``)
        - ``biaffine_only`` -- only include biaffine polynomials (default: object default)
        - ``correct_only`` -- only include correct polynomials (default: object default)

        EXAMPLES::

            sage: sr = mq.SR(1, 1, 1, 8, gf2=True)
            sage: len(sr._inversion_polynomials_single_sbox())
            24
            sage: len(sr._inversion_polynomials_single_sbox(correct_only=True))
            23
            sage: len(sr._inversion_polynomials_single_sbox(biaffine_only=False))
            40
            sage: len(sr._inversion_polynomials_single_sbox(biaffine_only=False, correct_only=True))
            39
        """
        e = self.e

        if biaffine_only is None:
            biaffine_only = self._biaffine_only
        if correct_only is None:
            correct_only = self._correct_only

        if x is None and w is None:
            # make sure it prints like in the book.
            names = ["w%d" % i for i in reversed(range(e))] + ["x%d" % i for i in reversed(range(e))]
            P = PolynomialRing(GF(2), e*2, names, order='lex')
            x = matrix(P, e, 1, P.gens()[e:])
            w = matrix(P, e, 1, P.gens()[:e])
        else:
            if isinstance(x, (tuple, list)):
                P = x[0].parent()
            elif isinstance(x, Matrix):
                P = x.base_ring()
            else:
                raise TypeError("x not understood")

            if isinstance(x, (tuple, list)):
                x = matrix(P, e, 1, x)
            if isinstance(w, (tuple, list)):
                w = matrix(P, e, 1, w)

        T = self._mul_matrix(self.k.gen())
        o = matrix(P, e, 1, [0]*(e-1) + [1])

        columns = []
        for i in reversed(range(e)):
            columns.append((T**i * w).list())
        Cw = matrix(P, e, e, columns).transpose()

        columns = []
        for i in reversed(range(e)):
            columns.append((T**i * x).list())
        Cx = matrix(P, e, e, columns).transpose()

        S = self._square_matrix()

        l = []
        if correct_only:
            l.append( (Cw * x + o).list()[:-1] )
        else:
            l.append( (Cw * x + o).list() )
        l.append( (Cw * S * x + x).list() )
        l.append( (Cx * S * w + w).list() )
        if not biaffine_only:
            l.append( ((Cw * S**2 + Cx*S)*x).list() )
            l.append( ((Cx * S**2 + Cw*S)*w).list() )

        return sum(l, [])

    def inversion_polynomials(self, xi, wi, length):
        """
        Return polynomials to represent the inversion in the AES S-Box.

        INPUT:

        - ``xi`` -- output variables

        - ``wi`` -- input variables

        - ``length`` -- length of both lists

        EXAMPLES::

            sage: sr = mq.SR(1, 1, 1, 8, gf2=True)
            sage: xi = sr.vars('x', 1)                                                  # needs sage.rings.polynomial.pbori
            sage: wi = sr.vars('w', 1)                                                  # needs sage.rings.polynomial.pbori
            sage: sr.inversion_polynomials(xi, wi, len(xi))[:3]                         # needs sage.rings.polynomial.pbori
            [x100*w100 + x100*w102 + x100*w103 + x100*w107 + x101*w101 + x101*w102 + x101*w106 + x102*w100 + x102*w101 + x102*w105 + x103*w100 + x103*w104 + x104*w103 + x105*w102 + x106*w101 + x107*w100,
             x100*w101 + x100*w103 + x100*w104 + x101*w100 + x101*w102 + x101*w103 + x101*w107 + x102*w101 + x102*w102 + x102*w106 + x103*w100 + x103*w101 + x103*w105 + x104*w100 + x104*w104 + x105*w103 + x106*w102 + x107*w101,
             x100*w102 + x100*w104 + x100*w105 + x101*w101 + x101*w103 + x101*w104 + x102*w100 + x102*w102 + x102*w103 + x102*w107 + x103*w101 + x103*w102 + x103*w106 + x104*w100 + x104*w101 + x104*w105 + x105*w100 + x105*w104 + x106*w103 + x107*w102]
        """
        if isinstance(xi, Matrix):
            xi = xi.list()
        if isinstance(wi, Matrix):
            wi = wi.list()

        e = self.e
        l = []
        for j in range(0, length, e):
            l += self.inversion_polynomials_single_sbox(xi[j:j+e], wi[j:j+e])
        return l

    def field_polynomials(self, name, i, l=None):
        r"""
        Return list of field polynomials for a given round ``i`` and
        name ``name``.

        INPUT:

        - ``name`` -- variable name
        - ``i`` -- round number
        - ``l`` -- length of variable list (default: ``None`` = r\*c)

        EXAMPLES::

            sage: sr = mq.SR(3, 1, 1, 8, gf2=True, polybori=False)
            sage: sr.field_polynomials('x', 2)
            [x200^2 + x200, x201^2 + x201,
             x202^2 + x202, x203^2 + x203,
             x204^2 + x204, x205^2 + x205,
             x206^2 + x206, x207^2 + x207]

        ::

            sage: sr = mq.SR(3, 1, 1, 8, gf2=True, polybori=True)
            sage: sr.field_polynomials('x', 2)
            []
        """
        r = self._r
        c = self._c
        e = self._e

        if l is None:
            l = r*c

        if self._polybori:
            return []
        _vars = self.vars(name, i, l, e)
        return [_vars[e*j+k]**2 - _vars[e*j+k] for j in range(l) for k in range(e)]


class SR_gf2_2(SR_gf2):
    """
    This is an example how to customize the SR constructor.

    In this example, we replace the S-Box inversion polynomials by the
    polynomials generated by the S-Box class.
    """
    def inversion_polynomials_single_sbox(self, x=None, w=None, biaffine_only=None, correct_only=None, groebner=False):
        """
        Return inversion polynomials of a single S-Box.

        INPUT:

        - ``x`` -- output variables (default: ``None``)
        - ``w`` -- input variables  (default: ``None``)
        - ``biaffine_only`` -- ignored (always ``False``)
        - ``correct_only`` -- ignored (always ``True``)
        - ``groebner`` -- precompute the Groebner basis for this S-Box (default: ``False``)

        EXAMPLES::

            sage: from sage.crypto.mq.sr import SR_gf2_2
            sage: e = 4
            sage: sr = SR_gf2_2(1, 1, 1, e)
            sage: P = PolynomialRing(GF(2),['x%d'%i for i in range(e)] + ['w%d'%i for i in range(e)],order='lex')
            sage: X,W = P.gens()[:e],P.gens()[e:]
            sage: sr.inversion_polynomials_single_sbox(X, W, groebner=True)             # needs sage.libs.singular
            [x0 + w0*w1*w2 + w0*w1 + w0*w2 + w0*w3 + w0 + w1 + w2,
             x1 + w0*w1*w3 + w0*w3 + w0 + w1*w3 + w1 + w2*w3,
             x2 + w0*w2*w3 + w0*w2 + w0 + w1*w2 + w1*w3 + w2*w3,
             x3 + w0*w1*w2 + w0 + w1*w2*w3 + w1*w2 + w1*w3 + w1 + w2 + w3]

            sage: from sage.crypto.mq.sr import SR_gf2_2
            sage: e = 4
            sage: sr = SR_gf2_2(1, 1, 1, e)
            sage: sr.inversion_polynomials_single_sbox()                                # needs sage.libs.singular
            [w3*w1 + w3*w0 + w3*x2 + w3*x1 + w3 + w2*w1 + w1 + x3 + x2 + x1,
             w3*w2 + w3*w1 + w3*x3 + w2 + w1 + x3,
             w3*w2 + w3*w1 + w3*x2 + w3 + w2*x3 + x2 + x1,
             w3*w2 + w3*w1 + w3*x3 + w3*x2 + w3*x1 + w3 + w2*x2 + w0 + x3 + x2 + x1 + x0,
             w3*w2 + w3*w1 + w3*x1 + w3*x0 + w2*x1 + w0 + x3 + x0,
             w3*w2 + w3*w1 + w3*w0 + w3*x2 + w3*x1 + w2*w0 + w2*x0 + w0 + x3 + x2 + x1 + x0,
             w3*w2 + w3*x1 + w3 + w2*w0 + w1*w0 + w1 + x3 + x2,
             w3*w2 + w3*w1 + w3*x1 + w1*x3 + x3 + x2 + x1,
             w3*x3 + w3*x2 + w3*x0 + w3 + w1*x2 + w1 + w0 + x2 + x0,
             w3*w2 + w3*w1 + w3*x2 + w3*x1 + w1*x1 + w1 + w0 + x2 + x0,
             w3*w2 + w3*w1 + w3*w0 + w3*x3 + w3*x1 + w2*w0 + w1*x0 + x3 + x2,
             w3*w2 + w3*w1 + w3*x2 + w3*x1 + w3*x0 + w3 + w1 + w0*x3 + x3 + x2,
             w3*w2 + w3*w1 + w3*w0 + w3*x3 + w3 + w2*w0 + w1 + w0*x2 + x3 + x2,
             w3*w0 + w3*x2 + w2*w0 + w0*x1 + w0 + x3 + x1 + x0,
             w3*w0 + w3*x3 + w3*x0 + w2*w0 + w1 + w0*x0 + w0 + x3 + x2,
             w3*w2 + w3 + w1 + x3*x2 + x3 + x1,
             w3*w2 + w3*x3 + w1 + x3*x1 + x3 + x2,
             w3*w2 + w3*w0 + w3*x3 + w3*x2 + w3*x1 + w0 + x3*x0 + x1 + x0,
             w3*w2 + w3*w1 + w3*w0 + w3*x3 + w1 + w0 + x2*x1 + x2 + x0,
             w3*w2 + w2*w0 + w1 + x3 + x2*x0,
             w3*x3 + w3*x1 + w2*w0 + w1 + x3 + x2 + x1*x0 + x1]

        TESTS:

        Note that ``biaffine_only`` and ``correct_only`` are always
        ignored. The former is always false while the second is always
        true. They are only accepted for compatibility with the base
        class.

            sage: from sage.crypto.mq.sr import SR_gf2_2
            sage: e = 4
            sage: sr = SR_gf2_2(1, 1, 1, e)
            sage: l = sr.inversion_polynomials_single_sbox()                            # needs sage.libs.singular
            sage: l == sr.inversion_polynomials_single_sbox(biaffine_only=True, correct_only=False)                     # needs sage.libs.singular
            True
        """
        e = self.e
        if x is None and w is None:
            # make sure it prints like in the book.
            names = ["w%d" % i for i in reversed(range(e))] + ["x%d" % i for i in reversed(range(e))]
            P = PolynomialRing(GF(2), e*2, names, order='lex')
            x = P.gens()[e:]
            w = P.gens()[:e]

        S = self.sbox(inversion_only=True)
        F = S.polynomials(w, x, degree=e-2, groebner=groebner)
        return F


class AllowZeroInversionsContext:
    """
    Temporarily allow zero inversion.
    """
    def __init__(self, sr):
        """
        EXAMPLES::

            sage: from sage.crypto.mq.sr import AllowZeroInversionsContext
            sage: sr = mq.SR(1,2,2,4)
            sage: with AllowZeroInversionsContext(sr):
            ....:     sr.sub_byte(0)
            a^2 + a
        """
        self.sr = sr

    def __enter__(self):
        """
        EXAMPLES::

            sage: from sage.crypto.mq.sr import AllowZeroInversionsContext
            sage: sr = mq.SR(1,2,2,4)
            sage: sr.sub_byte(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: A zero inversion occurred during an encryption or key schedule.

            sage: with AllowZeroInversionsContext(sr):
            ....:     sr.sub_byte(0)
            a^2 + a
        """
        self.allow_zero_inversions = self.sr._allow_zero_inversions
        self.sr._allow_zero_inversions = True

    def __exit__(self, typ, value, tb):
        """
        EXAMPLES::

            sage: from sage.crypto.mq.sr import AllowZeroInversionsContext
            sage: sr = mq.SR(1,2,2,4)
            sage: with AllowZeroInversionsContext(sr):
            ....:     sr.sub_byte(0)
            a^2 + a
            sage: sr._allow_zero_inversions
            False
        """
        self.sr._allow_zero_inversions = self.allow_zero_inversions


def test_consistency(max_n=2, **kwargs):
    r"""
    Test all combinations of ``r``, ``c``, ``e`` and ``n`` in ``(1,
    2)`` for consistency of random encryptions and their polynomial
    systems. `\GF{2}` and `\GF{2^e}` systems are tested. This test takes
    a while.

    INPUT:

    - ``max_n`` -- maximal number of rounds to consider (default: 2)
    - ``kwargs`` -- are passed to the SR constructor

    TESTS:

    The following test called with ``max_n`` = 2 requires a LOT of RAM
    (much more than 2GB).  Since this might cause the doctest to fail
    on machines with "only" 2GB of RAM, we test ``max_n`` = 1, which
    has a more reasonable memory usage. ::

        sage: from sage.crypto.mq.sr import test_consistency
        sage: test_consistency(1)  # needs sage.rings.polynomial.pbori, long time (65s on sage.math, 2012)
        True
    """
    consistent = True
    for r in (1, 2, 4):
        for c in (1, 2, 4):
            for e in (4, 8):
                for n in range(1, max_n+1):
                    for gf2 in (True, False):
                        zero_division = True
                        while zero_division:
                            sr = SR(n, r, c, e, gf2=gf2, **kwargs)
                            try:
                                F, s = sr.polynomial_system()
                                F = F.subs(s)
                                consistent &= (F.groebner_basis()[0] != 1)
                                if not consistent:
                                    print(str(sr) + " is not consistent")
                                zero_division = False

                            except ZeroDivisionError:
                                pass
    return consistent
