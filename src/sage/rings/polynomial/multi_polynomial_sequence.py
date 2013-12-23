"""
Polynomial Sequences

We call a finite list of polynomials a ``Polynomial Sequence``.

Polynomial sequences in Sage can optionally be viewed as consisting of
various parts or sub-sequences. These kind of polynomial sequences
which naturally split into parts arise naturally for example in
algebraic cryptanalysis of symmetric cryptographic primitives. The
most prominent examples of these systems are: the small scale variants
of the AES [CMR05]_ (cf. :func:`sage.crypto.mq.sr.SR`) and Flurry/Curry [BPW06]_. By
default, a polynomial sequence has exactly one part.

AUTHORS:

- Martin Albrecht (2007ff): initial version
- Martin Albrecht (2009): refactoring, clean-up, new functions
- Martin Albrecht (2011): refactoring, moved to sage.rings.polynomial
- Alex Raichev (2011-06): added algebraic_dependence()
- Charles Bouillaguet (2013-1): added solve()

EXAMPLES:

As an example consider a small scale variant of the AES::

   sage: sr = mq.SR(2,1,2,4,gf2=True,polybori=True)
   sage: sr
   SR(2,1,2,4)

We can construct a polynomial sequence for a random plaintext-ciphertext
pair and study it::

   sage: set_random_seed(1)
   sage: F,s = sr.polynomial_system()
   sage: F
   Polynomial Sequence with 112 Polynomials in 64 Variables

   sage: r2 = F.part(2); r2
   (w200 + k100 + x100 + x102 + x103,
    w201 + k101 + x100 + x101 + x103 + 1,
    w202 + k102 + x100 + x101 + x102 + 1,
    w203 + k103 + x101 + x102 + x103,
    w210 + k110 + x110 + x112 + x113,
    w211 + k111 + x110 + x111 + x113 + 1,
    w212 + k112 + x110 + x111 + x112 + 1,
    w213 + k113 + x111 + x112 + x113,
    x100*w100 + x100*w103 + x101*w102 + x102*w101 + x103*w100,
    x100*w100 + x100*w101 + x101*w100 + x101*w103 + x102*w102 + x103*w101,
    x100*w101 + x100*w102 + x101*w100 + x101*w101 + x102*w100 + x102*w103 + x103*w102,
    x100*w100 + x100*w102 + x100*w103 + x101*w100 + x101*w101 + x102*w102 + x103*w100 + x100,
    x100*w101 + x100*w103 + x101*w101 + x101*w102 + x102*w100 + x102*w103 + x103*w101 + x101,
    x100*w100 + x100*w102 + x101*w100 + x101*w102 + x101*w103 + x102*w100 + x102*w101 + x103*w102 + x102,
    x100*w101 + x100*w102 + x101*w100 + x101*w103 + x102*w101 + x103*w103 + x103,
    x100*w100 + x100*w101 + x100*w103 + x101*w101 + x102*w100 + x102*w102 + x103*w100 + w100,
    x100*w102 + x101*w100 + x101*w101 + x101*w103 + x102*w101 + x103*w100 + x103*w102 + w101,
    x100*w100 + x100*w101 + x100*w102 + x101*w102 + x102*w100 + x102*w101 + x102*w103 + x103*w101 + w102,
    x100*w101 + x101*w100 + x101*w102 + x102*w100 + x103*w101 + x103*w103 + w103,
    x100*w102 + x101*w101 + x102*w100 + x103*w103 + 1,
    x110*w110 + x110*w113 + x111*w112 + x112*w111 + x113*w110,
    x110*w110 + x110*w111 + x111*w110 + x111*w113 + x112*w112 + x113*w111,
    x110*w111 + x110*w112 + x111*w110 + x111*w111 + x112*w110 + x112*w113 + x113*w112,
    x110*w110 + x110*w112 + x110*w113 + x111*w110 + x111*w111 + x112*w112 + x113*w110 + x110,
    x110*w111 + x110*w113 + x111*w111 + x111*w112 + x112*w110 + x112*w113 + x113*w111 + x111,
    x110*w110 + x110*w112 + x111*w110 + x111*w112 + x111*w113 + x112*w110 + x112*w111 + x113*w112 + x112,
    x110*w111 + x110*w112 + x111*w110 + x111*w113 + x112*w111 + x113*w113 + x113,
    x110*w110 + x110*w111 + x110*w113 + x111*w111 + x112*w110 + x112*w112 + x113*w110 + w110,
    x110*w112 + x111*w110 + x111*w111 + x111*w113 + x112*w111 + x113*w110 + x113*w112 + w111,
    x110*w110 + x110*w111 + x110*w112 + x111*w112 + x112*w110 + x112*w111 + x112*w113 + x113*w111 + w112,
    x110*w111 + x111*w110 + x111*w112 + x112*w110 + x113*w111 + x113*w113 + w113,
    x110*w112 + x111*w111 + x112*w110 + x113*w113 + 1)

We separate the system in independent subsystems::

    sage: C = Sequence(r2).connected_components(); C
    [[w213 + k113 + x111 + x112 + x113,
     w212 + k112 + x110 + x111 + x112 + 1,
     w211 + k111 + x110 + x111 + x113 + 1,
     w210 + k110 + x110 + x112 + x113,
     x110*w112 + x111*w111 + x112*w110 + x113*w113 + 1,
     x110*w112 + x111*w110 + x111*w111 + x111*w113 + x112*w111 + x113*w110 + x113*w112 + w111,
     x110*w111 + x111*w110 + x111*w112 + x112*w110 + x113*w111 + x113*w113 + w113,
     x110*w111 + x110*w113 + x111*w111 + x111*w112 + x112*w110 + x112*w113 + x113*w111 + x111,
     x110*w111 + x110*w112 + x111*w110 + x111*w113 + x112*w111 + x113*w113 + x113,
     x110*w111 + x110*w112 + x111*w110 + x111*w111 + x112*w110 + x112*w113 + x113*w112,
     x110*w110 + x110*w113 + x111*w112 + x112*w111 + x113*w110,
     x110*w110 + x110*w112 + x111*w110 + x111*w112 + x111*w113 + x112*w110 + x112*w111 + x113*w112 + x112,
     x110*w110 + x110*w112 + x110*w113 + x111*w110 + x111*w111 + x112*w112 + x113*w110 + x110,
     x110*w110 + x110*w111 + x111*w110 + x111*w113 + x112*w112 + x113*w111,
     x110*w110 + x110*w111 + x110*w113 + x111*w111 + x112*w110 + x112*w112 + x113*w110 + w110,
     x110*w110 + x110*w111 + x110*w112 + x111*w112 + x112*w110 + x112*w111 + x112*w113 + x113*w111 + w112],
    [w203 + k103 + x101 + x102 + x103,
    w202 + k102 + x100 + x101 + x102 + 1,
    w201 + k101 + x100 + x101 + x103 + 1,
    w200 + k100 + x100 + x102 + x103,
    x100*w102 + x101*w101 + x102*w100 + x103*w103 + 1,
    x100*w102 + x101*w100 + x101*w101 + x101*w103 + x102*w101 + x103*w100 + x103*w102 + w101,
    x100*w101 + x101*w100 + x101*w102 + x102*w100 + x103*w101 + x103*w103 + w103,
    x100*w101 + x100*w103 + x101*w101 + x101*w102 + x102*w100 + x102*w103 + x103*w101 + x101,
    x100*w101 + x100*w102 + x101*w100 + x101*w103 + x102*w101 + x103*w103 + x103, x100*w101 + x100*w102 + x101*w100 + x101*w101 + x102*w100 + x102*w103 + x103*w102,
    x100*w100 + x100*w103 + x101*w102 + x102*w101 + x103*w100,
    x100*w100 + x100*w102 + x101*w100 + x101*w102 + x101*w103 + x102*w100 + x102*w101 + x103*w102 + x102,
    x100*w100 + x100*w102 + x100*w103 + x101*w100 + x101*w101 + x102*w102 + x103*w100 + x100,
    x100*w100 + x100*w101 + x101*w100 + x101*w103 + x102*w102 + x103*w101,
    x100*w100 + x100*w101 + x100*w103 + x101*w101 + x102*w100 + x102*w102 + x103*w100 + w100,
    x100*w100 + x100*w101 + x100*w102 + x101*w102 + x102*w100 + x102*w101 + x102*w103 + x103*w101 + w102]]
    sage: C[0].groebner_basis()
    Polynomial Sequence with 30 Polynomials in 16 Variables

and compute the coefficient matrix::

    sage: A,v = Sequence(r2).coefficient_matrix()
    sage: A.rank()
    32

Using these building blocks we can implement a simple XL algorithm
easily::

    sage: sr = mq.SR(1,1,1,4, gf2=True, polybori=True, order='lex')
    sage: F,s = sr.polynomial_system()

    sage: monomials = [a*b for a in F.variables() for b in F.variables() if a<b]
    sage: len(monomials)
    190
    sage: F2 = Sequence(map(mul, cartesian_product_iterator((monomials, F))))
    sage: A,v = F2.coefficient_matrix(sparse=False)
    sage: A.echelonize()
    sage: A
    6840 x 4474 dense matrix over Finite Field of size 2 (type 'print A.str()' to see all of the entries)
    sage: A.rank()
    4056
    sage: A[4055]*v
    (k001*k003)

TEST::

    sage: P.<x,y> = PolynomialRing(QQ)
    sage: I = [[x^2 + y^2], [x^2 - y^2]]
    sage: F = Sequence(I, P)
    sage: loads(dumps(F)) == F
    True

.. NOTE::

   In many other computer algebra systems (cf. Singular) this class
   would be called ``Ideal`` but an ideal is a very distinct object
   from its generators and thus this is not an ideal in Sage.

.. [BPW06] J. Buchmann, A. Pychkine, R.-P. Weinmann
   *Block Ciphers Sensitive to Groebner Basis Attacks*
   in Topics in Cryptology -- CT RSA'06; LNCS 3860; pp. 313--331; Springer Verlag 2006;
   pre-print available at http://eprint.iacr.org/2005/200

Classes
-------
"""

from types import GeneratorType
from sage.misc.package import is_package_installed

from sage.structure.sequence import Sequence, Sequence_generic

from sage.rings.infinity import Infinity
from sage.rings.finite_rings.constructor import FiniteField as GF
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.quotient_ring import is_QuotientRing
from sage.rings.quotient_ring_element import QuotientRingElement
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.rings.polynomial.multi_polynomial import is_MPolynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.interfaces.singular import singular

def is_PolynomialSequence(F):
    """
    Return ``True`` if ``F`` is a ``PolynomialSequence``.

    INPUT:

    - ``F`` - anything

    EXAMPLE::

        sage: P.<x,y> = PolynomialRing(QQ)
        sage: I = [[x^2 + y^2], [x^2 - y^2]]
        sage: F = Sequence(I, P); F
        [x^2 + y^2, x^2 - y^2]

        sage: from sage.rings.polynomial.multi_polynomial_sequence import is_PolynomialSequence
        sage: is_PolynomialSequence(F)
        True
    """
    return isinstance(F,PolynomialSequence_generic)

def PolynomialSequence(arg1, arg2=None, immutable=False, cr=False, cr_str=None):
    """
    Construct a new polynomial sequence object.

    INPUT:

    - ``arg1`` - a multivariate polynomial ring, an ideal or a matrix

    - ``arg2`` - an iterable object of parts or polynomials
      (default:``None``)

      - ``immutable`` - if ``True`` the sequence is immutable (default: ``False``)

      - ``cr`` - print a line break after each element (default: ``False``)

      - ``cr_str`` - print a line break after each element if 'str' is
        called (default: ``None``)

    EXAMPLES::

        sage: P.<a,b,c,d> = PolynomialRing(GF(127),4)
        sage: I = sage.rings.ideal.Katsura(P)

    If a list of tuples is provided, those form the parts::

        sage: F = Sequence([I.gens(),I.gens()], I.ring()); F # indirect doctest
        [a + 2*b + 2*c + 2*d - 1,
         a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a,
         2*a*b + 2*b*c + 2*c*d - b,
         b^2 + 2*a*c + 2*b*d - c,
         a + 2*b + 2*c + 2*d - 1,
         a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a,
         2*a*b + 2*b*c + 2*c*d - b,
         b^2 + 2*a*c + 2*b*d - c]
        sage: F.nparts()
        2

    If an ideal is provided, the generators are used::

        sage: Sequence(I)
        [a + 2*b + 2*c + 2*d - 1,
         a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a,
         2*a*b + 2*b*c + 2*c*d - b,
         b^2 + 2*a*c + 2*b*d - c]

    If a list of polynomials is provided, the system has only one
    part::

        sage: F = Sequence(I.gens(), I.ring()); F
        [a + 2*b + 2*c + 2*d - 1,
         a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a,
         2*a*b + 2*b*c + 2*c*d - b,
         b^2 + 2*a*c + 2*b*d - c]
         sage: F.nparts()
         1
    """

    from sage.matrix.matrix import is_Matrix

    if is_MPolynomialRing(arg1) or is_QuotientRing(arg1):
        ring = arg1
        gens = arg2

    elif is_MPolynomialRing(arg2) or is_QuotientRing(arg2):
        ring = arg2
        gens = arg1

    elif is_Matrix(arg1) and arg2 is None:
        ring = arg1.base_ring()
        gens = arg1.list()

    elif isinstance(arg1, MPolynomialIdeal) and arg2 is None:
        ring = arg1.ring()
        gens = arg1.gens()

    elif isinstance(arg1, (list,tuple,GeneratorType)) and arg2 is None:
        gens = arg1

        try:
            e = iter(gens).next()
        except StopIteration:
            raise ValueError("Cannot determine ring from provided information.")

        if is_MPolynomial(e) or isinstance(e, QuotientRingElement):
            ring = e.parent()
        else:
            ring = iter(e).next().parent()
    else:
        raise TypeError("Cannot understand input.")

    try:
        e = iter(gens).next()

        if is_MPolynomial(e) or isinstance(e, QuotientRingElement):
            gens = tuple(gens)
            parts = (gens,)
            if not all(f.parent() is ring for f in gens):
                parts = ((ring(f) for f in gens),)
        else:
            parts = []
            _gens = []
            for part in gens:
                _part = []
                for gen in part:
                    if not gen.parent() is ring:
                        ring(gen)
                    _part.append(gen)
                _gens.extend(_part)
                parts.append(tuple(_part))
            gens = _gens
    except StopIteration:
        gens = tuple()
        parts = ((),)

    k = ring.base_ring()

    try: c = k.characteristic()
    except NotImplementedError: c = -1

    if c != 2:
        return PolynomialSequence_generic(parts, ring, immutable=immutable, cr=cr, cr_str=cr_str)
    elif k.degree() == 1:
        return PolynomialSequence_gf2(parts, ring, immutable=immutable, cr=cr, cr_str=cr_str)
    elif k.degree() > 1:
        return PolynomialSequence_gf2e(parts, ring, immutable=immutable, cr=cr, cr_str=cr_str)

class PolynomialSequence_generic(Sequence_generic):
    def __init__(self, parts, ring, immutable=False, cr=False, cr_str=None):
        """
        Construct a new system of multivariate polynomials.

        INPUT:

        - ``part`` - a list of lists with polynomials

        -  ``ring`` - a multivariate polynomial ring

        - ``immutable`` - if ``True`` the sequence is immutable (default: ``False``)

        - ``cr`` - print a line break after each element (default: ``False``)

        - ``cr_str`` - print a line break after each element if 'str'
          is called (default: ``None``)

        EXAMPLES::

            sage: P.<a,b,c,d> = PolynomialRing(GF(127),4)
            sage: I = sage.rings.ideal.Katsura(P)

            sage: Sequence([I.gens()], I.ring()) # indirect doctest
            [a + 2*b + 2*c + 2*d - 1, a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a, 2*a*b + 2*b*c + 2*c*d - b, b^2 + 2*a*c + 2*b*d - c]

        If an ideal is provided, the generators are used.::

            sage: Sequence(I)
            [a + 2*b + 2*c + 2*d - 1, a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a, 2*a*b + 2*b*c + 2*c*d - b, b^2 + 2*a*c + 2*b*d - c]

        If a list of polynomials is provided, the system has only one
        part.::

            sage: Sequence(I.gens(), I.ring())
            [a + 2*b + 2*c + 2*d - 1, a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a, 2*a*b + 2*b*c + 2*c*d - b, b^2 + 2*a*c + 2*b*d - c]
        """

        Sequence_generic.__init__(self, sum(parts,tuple()), ring, check=False, immutable=immutable,
                                  cr=cr, cr_str=cr_str, use_sage_types=True)
        self._ring = ring
        self._parts = parts

    def __copy__(self):
        """
        Return a copy of this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: copy(F) # indirect doctest
            Polynomial Sequence with 40 Polynomials in 20 Variables
            sage: type(F) == type(copy(F))
            True
        """
        return self.__class__(self._parts, self._ring, immutable=self.is_immutable())

    def ring(self):
        """
        Return the polynomial ring all elements live in.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: print F.ring().repr_long()
            Polynomial Ring
             Base Ring : Finite Field of size 2
                  Size : 20 Variables
              Block  0 : Ordering : deglex
                         Names    : k100, k101, k102, k103, x100, x101, x102, x103, w100, w101, w102, w103, s000, s001, s002, s003
              Block  1 : Ordering : deglex
                         Names    : k000, k001, k002, k003
        """
        return self._ring

    universe = ring

    def nparts(self):
        """
        Return number of parts of this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F.nparts()
            4
        """
        return len(self._parts)

    def parts(self):
        """
        Return a tuple of parts of this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: l = F.parts()
            sage: len(l)
            4
        """
        return tuple(self._parts)

    def part(self, i):
        """
        Return ``i``-th part of this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: R0 = F.part(1)
            sage: R0
            (k000^2 + k001, k001^2 + k002, k002^2 + k003, k003^2 + k000)
        """
        return self._parts[i]

    def ideal(self):
        """
        Return ideal spanned by the elements of this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: P = F.ring()
            sage: I = F.ideal()
            sage: I.elimination_ideal(P('s000*s001*s002*s003*w100*w101*w102*w103*x100*x101*x102*x103'))
            Ideal (k002 + (a^3 + a + 1)*k003 + (a^2 + 1),
                   k001 + (a^3)*k003, k000 + (a)*k003 + (a^2),
                   k103 + k003 + (a^2 + a + 1),
                   k102 + (a^3 + a + 1)*k003 + (a + 1),
                   k101 + (a^3)*k003 + (a^2 + a + 1),
                   k100 + (a)*k003 + (a),
                   k003^2 + (a)*k003 + (a^2))
            of Multivariate Polynomial Ring in k100, k101, k102, k103, x100, x101, x102, x103,
            w100, w101, w102, w103, s000, s001, s002, s003, k000, k001, k002, k003 over Finite Field in a of size 2^4
        """
        return self._ring.ideal(tuple(self))

    def groebner_basis(self, *args, **kwargs):
        """
        Compute and return a Groebner basis for the ideal spanned by
        the polynomials in this system.

        INPUT:

        - ``args`` - list of arguments passed to
          ``MPolynomialIdeal.groebner_basis`` call

        - ``kwargs`` - dictionary of arguments passed to
          ``MPolynomialIdeal.groebner_basis`` call

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: gb = F.groebner_basis()
            sage: Ideal(gb).basis_is_groebner()
            True
        """
        return self.ideal().groebner_basis(*args, **kwargs)

    def monomials(self):
        """
        Return an unordered tuple of monomials in this polynomial system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: len(F.monomials())
            49
        """
        M = set()
        for f in self:
            for m in f.monomials():
                M.add(m)
        return tuple(M)

    def nmonomials(self):
        """
        Return the number of monomials present in this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F.nmonomials()
            49
        """
        return len(self.monomials())

    def variables(self):
        """
        Return all variables present in this system. This tuple may or
        may not be equal to the generators of the ring of this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F.variables()[:10]
            (k003, k002, k001, k000, s003, s002, s001, s000, w103, w102)
        """
        V = set()
        for f in self:
            for v in f.variables():
                V.add(v)
        return tuple(sorted(V))

    def nvariables(self):
        """
        Return number of variables present in this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F.nvariables()
            20
        """
        return len(self.variables())

    def algebraic_dependence(self):
        r"""
        Returns the ideal of annihilating polynomials for the
        polynomials in ``self``, if those polynomials are algebraically
        dependent.
        Otherwise, returns the zero ideal.

        OUTPUT:

        If the polynomials `f_1,\ldots,f_r` in ``self`` are algebraically
        dependent, then the output is the ideal
        `\{F \in K[T_1,\ldots,T_r] : F(f_1,\ldots,f_r) = 0\}` of
        annihilating polynomials of `f_1,\ldots,f_r`.
        Here `K` is the coefficient ring of polynomial ring of `f_1,\ldots,f_r`
        and `T_1,\ldots,T_r` are new indeterminates.
        If `f_1,\ldots,f_r` are algebraically independent, then the output
        is the zero ideal in `K[T_1,\ldots,T_r]`.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = Sequence([x, x*y])
            sage: I = S.algebraic_dependence(); I
            Ideal (0) of Multivariate Polynomial Ring in T0, T1 over Rational Field

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = Sequence([x, (x^2 + y^2 - 1)^2, x*y - 2])
            sage: I = S.algebraic_dependence(); I
            Ideal (16 + 32*T2 - 8*T0^2 + 24*T2^2 - 8*T0^2*T2 + 8*T2^3 + 9*T0^4 - 2*T0^2*T2^2 + T2^4 - T0^4*T1 + 8*T0^4*T2 - 2*T0^6 + 2*T0^4*T2^2 + T0^8) of Multivariate Polynomial Ring in T0, T1, T2 over Rational Field
            sage: [F(S) for F in I.gens()]
            [0]

        ::

            sage: R.<x,y> = PolynomialRing(GF(7))
            sage: S = Sequence([x, (x^2 + y^2 - 1)^2, x*y - 2])
            sage: I = S.algebraic_dependence(); I
            Ideal (2 - 3*T2 - T0^2 + 3*T2^2 - T0^2*T2 + T2^3 + 2*T0^4 - 2*T0^2*T2^2 + T2^4 - T0^4*T1 + T0^4*T2 - 2*T0^6 + 2*T0^4*T2^2 + T0^8) of Multivariate Polynomial Ring in T0, T1, T2 over Finite Field of size 7
            sage: [F(S) for F in I.gens()]
            [0]

        .. NOTE::

            This function's code also works for sequences of polynomials from a
            univariate polynomial ring, but i don't know where in the Sage codebase
            to put it to use it to that effect.

        AUTHORS:

        - Alex Raichev (2011-06-22)
        """
        R = self.ring()
        K = R.base_ring()
        Xs = list(R.gens())
        r = len(self)
        d = len(Xs)

        # Expand R by r new variables.
        T = 'T'
        while T in [str(x) for x in Xs]:
            T = T+'T'
        Ts = [T + str(j) for j in range(r)]
        RR = PolynomialRing(K,d+r,tuple(Xs+Ts))
        Vs = list(RR.gens())
        Xs = Vs[0 :d]
        Ts = Vs[d:]

        J = RR.ideal([ Ts[j] - RR(self[j]) for j in range(r)])
        JJ = J.elimination_ideal(Xs)
        # By the elimination theorem, JJ is the kernel of the ring homorphism
        # `phi:K[\bar T] \to K[\bar X]` that fixes `K` and sends each
        # `T_i` to `f_i`.
        # So JJ is the ideal of annihilating polynomials of `f_1,\ldots,f_r`,
        # which is the zero ideal in case `f_1,\ldots,f_r` are algebraically
        # independent.

        # Coerce JJ into `K[T_1,\ldots,T_r]`.
        # Choosing the negdeglex order simply because i find it useful in my work.
        RRR = PolynomialRing(K,r,tuple(Ts),order='negdeglex')
        return RRR.ideal(JJ.gens())

    def coefficient_matrix(self, sparse=True):
        """
        Return tuple ``(A,v)`` where ``A`` is the coefficient matrix
        of this system and ``v`` the matching monomial vector.

        Thus value of ``A[i,j]`` corresponds the coefficient of the
        monomial ``v[j]`` in the ``i``-th polynomial in this system.

        Monomials are order w.r.t. the term ordering of
        ``self.ring()`` in reverse order, i.e. such that the smallest
        entry comes last.

        INPUT:

        - ``sparse`` - construct a sparse matrix (default: ``True``)

        EXAMPLE::

            sage: P.<a,b,c,d> = PolynomialRing(GF(127),4)
            sage: I = sage.rings.ideal.Katsura(P)
            sage: I.gens()
            [a + 2*b + 2*c + 2*d - 1,
             a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a,
             2*a*b + 2*b*c + 2*c*d - b,
             b^2 + 2*a*c + 2*b*d - c]

            sage: F = Sequence(I)
            sage: A,v = F.coefficient_matrix()
            sage: A
            [  0   0   0   0   0   0   0   0   0   1   2   2   2 126]
            [  1   0   2   0   0   2   0   0   2 126   0   0   0   0]
            [  0   2   0   0   2   0   0   2   0   0 126   0   0   0]
            [  0   0   1   2   0   0   2   0   0   0   0 126   0   0]

            sage: v
            [a^2]
            [a*b]
            [b^2]
            [a*c]
            [b*c]
            [c^2]
            [b*d]
            [c*d]
            [d^2]
            [  a]
            [  b]
            [  c]
            [  d]
            [  1]

            sage: A*v
            [        a + 2*b + 2*c + 2*d - 1]
            [a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a]
            [      2*a*b + 2*b*c + 2*c*d - b]
            [        b^2 + 2*a*c + 2*b*d - c]
        """
        R = self.ring()

        m = sorted(self.monomials(),reverse=True)
        nm = len(m)
        f = tuple(self)
        nf = len(f)

        #construct dictionary for fast lookups
        v = dict( zip( m , range(len(m)) ) )

        from sage.matrix.constructor import Matrix

        A = Matrix( R.base_ring(), nf, nm, sparse=sparse )

        for x in range( nf ):
            poly = f[x]
            for y in poly.monomials():
                A[ x , v[y] ] = poly.monomial_coefficient(y)

        return  A, Matrix(R,nm,1,m)

    def subs(self, *args, **kwargs):
        """
        Substitute variables for every polynomial in this system and
        return a new system. See ``MPolynomial.subs`` for calling
        convention.

        INPUT:

        -  ``args`` - arguments to be passed to ``MPolynomial.subs``
        -  ``kwargs`` - keyword arguments to be passed to ``MPolynomial.subs``

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system(); F
            Polynomial Sequence with 40 Polynomials in 20 Variables
            sage: F = F.subs(s); F
            Polynomial Sequence with 40 Polynomials in 16 Variables
        """
        return PolynomialSequence(self._ring, [tuple([f.subs(*args,**kwargs) for f in r]) for r in self._parts])

    def _singular_(self):
        """
        Return Singular ideal representation of this system.

        EXAMPLE::

            sage: P.<a,b,c,d> = PolynomialRing(GF(127))
            sage: I = sage.rings.ideal.Katsura(P)
            sage: F = Sequence(I); F
            [a + 2*b + 2*c + 2*d - 1,
             a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a,
             2*a*b + 2*b*c + 2*c*d - b,
             b^2 + 2*a*c + 2*b*d - c]
            sage: F._singular_()
            a+2*b+2*c+2*d-1,
            a^2+2*b^2+2*c^2+2*d^2-a,
            2*a*b+2*b*c+2*c*d-b,
            b^2+2*a*c+2*b*d-c
        """
        return singular.ideal(list(self))

    def _magma_init_(self, magma):
        """
        Return Magma ideal representation of the ideal spanned by this
        system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True)
            sage: F,s = sr.polynomial_system()
            sage: F.set_immutable()
            sage: magma(F)  # indirect doctest; optional - magma
            Ideal of Polynomial ring of rank 20 over GF(2)
            Order: Graded Reverse Lexicographical
            Variables: k100, k101, k102, k103, x100, x101, x102, x103, w100, w101, w102, w103, s000, s001, s002, s003, k000, k001, k002, k003
            Basis:
            [
            ...
            ]
        """
        P = magma(self.ring()).name()
        v = [x._magma_init_(magma) for x in list(self)]
        return 'ideal<%s|%s>'%(P, ','.join(v))

    def _repr_(self):
        """
        Return a string representation of this system.

        EXAMPLE::

            sage: P.<a,b,c,d> = PolynomialRing(GF(127))
            sage: I = sage.rings.ideal.Katsura(P)
            sage: F = Sequence(I); F # indirect doctest
            [a + 2*b + 2*c + 2*d - 1,
             a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a,
             2*a*b + 2*b*c + 2*c*d - b,
             b^2 + 2*a*c + 2*b*d - c]

        If the system contains 20 or more polynomials, a short summary
        is printed::

            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True)
            sage: F,s = sr.polynomial_system(); F
            Polynomial Sequence with 36 Polynomials in 20 Variables

        """
        if len(self) < 20:
            return Sequence_generic._repr_(self)
        else:
            return "Polynomial Sequence with %d Polynomials in %d Variables"%(len(self),self.nvariables())

    def __add__(self, right):
        """
        Add polynomial systems together, i.e. create a union of their
        polynomials.

        EXAMPLE::

            sage: P.<a,b,c,d> = PolynomialRing(GF(127))
            sage: I = sage.rings.ideal.Katsura(P)
            sage: F = Sequence(I)
            sage: F + [a^127 + a]
            [a + 2*b + 2*c + 2*d - 1,
             a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a,
             2*a*b + 2*b*c + 2*c*d - b,
             b^2 + 2*a*c + 2*b*d - c,
             a^127 + a]

            sage: F + P.ideal([a^127 + a])
            [a + 2*b + 2*c + 2*d - 1,
             a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a,
             2*a*b + 2*b*c + 2*c*d - b,
             b^2 + 2*a*c + 2*b*d - c,
             a^127 + a]

            sage: F + Sequence([a^127 + a], P)
            [a + 2*b + 2*c + 2*d - 1,
             a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a,
             2*a*b + 2*b*c + 2*c*d - b,
             b^2 + 2*a*c + 2*b*d - c,
             a^127 + a]
        """
        if is_PolynomialSequence(right) and right.ring() == self.ring():
            return PolynomialSequence(self.ring(), self.parts() + right.parts())

        elif isinstance(right,(tuple,list)) and all(map(lambda x: x.parent() == self.ring(), right)):
            return PolynomialSequence(self.ring(), self.parts() + (right,))

        elif isinstance(right,MPolynomialIdeal) and (right.ring() is self.ring() or right.ring() == self.ring()):
            return PolynomialSequence(self.ring(), self.parts() + (right.gens(),))

        else:
            raise TypeError("right must be a system over same ring as self.")

    def connection_graph(self):
        """
        Return the graph which has the variables of this system as
        vertices and edges between two variables if they appear in the
        same polynomial.

        EXAMPLE::

            sage: B.<x,y,z> = BooleanPolynomialRing()
            sage: F = Sequence([x*y + y + 1, z + 1])
            sage: F.connection_graph()
            Graph on 3 vertices
        """
        V = sorted(self.variables())
        from sage.graphs.graph import Graph
        g = Graph()
        g.add_vertices(sorted(V))
        for f in self:
            v = f.variables()
            a,tail = v[0],v[1:]
            for b in tail:
                g.add_edge((a,b))
        return g

    def connected_components(self):
        """
        Split the polynomial system in systems which do not share any
        variables.

        EXAMPLE:

        As an example consider one part of AES, which naturally
        splits into four subsystems which are independent::

            sage: sr = mq.SR(2,4,4,8,gf2=True,polybori=True)
            sage: F,s = sr.polynomial_system()
            sage: Fz = Sequence(F.part(2))
            sage: Fz.connected_components()
            [Polynomial Sequence with 128 Polynomials in 128 Variables,
             Polynomial Sequence with 128 Polynomials in 128 Variables,
             Polynomial Sequence with 128 Polynomials in 128 Variables,
             Polynomial Sequence with 128 Polynomials in 128 Variables]
        """
        g = self.connection_graph()
        C = sorted(g.connected_components())

        P = [[] for _ in xrange(len(C))]
        for f in self:
            for i,c in enumerate(C):
                if len(set(f.variables()).difference(c)) == 0:
                    P[i].append(f)
                    break
        P = sorted([PolynomialSequence(sorted(p)) for p in P])
        return P

    def _groebner_strategy(self):
        """
        Return the Singular Groebner Strategy object.

        This object allows to compute normal forms efficiently, since
        all conversion overhead is avoided.

        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(GF(127))
            sage: F = Sequence([x*y + z, y + z + 1])
            sage: F._groebner_strategy()
            Groebner Strategy for ideal generated by 2 elements over
            Multivariate Polynomial Ring in x, y, z over Finite Field of size 127
        """
        from sage.libs.singular.groebner_strategy import GroebnerStrategy
        return GroebnerStrategy(self.ideal())

    def __reduce__(self):
        """
        TESTS::

            sage: P.<x,y,z> = PolynomialRing(GF(127))
            sage: F = Sequence([x*y + z, y + z + 1])
            sage: loads(dumps(F)) == F # indirect doctest
            True
        """
        return PolynomialSequence, (self._ring, self._parts)

class PolynomialSequence_gf2(PolynomialSequence_generic):
    """
    Polynomial Sequences over `\mathbb{F}_2`.
    """
    def eliminate_linear_variables(self, maxlength=Infinity, skip=None, return_reductors=False):
        """
        Return a new system where linear leading variables are
        eliminated if the tail of the polynomial has length at most
        ``maxlength``.

        INPUT:

        - ``maxlength`` - an optional upper bound on the number of
          monomials by which a variable is replaced. If
          ``maxlength==+Infinity`` then no condition is checked.
          (default: +Infinity).

        - ``skip`` - an optional callable to skip eliminations. It
          must accept two parameters and return either ``True`` or
          ``False``. The two parameters are the leading term and the
          tail of a polynomial (default: ``None``).

        - ``return_reductors`` - if ``True`` the list of polynomials
          with linear leading terms which were used for reduction is
          also returned (default: ``False``).

        OUTPUT:

        When ``return_reductors==True``, then a pair of sequences of
        boolean polynomials are returned, along with the promises that:

          1. The union of the two sequences spans the
             same boolean ideal as the argument of the method

          2. The second sequence only contains linear polynomials, and
             it forms a reduced groebner basis (they all have pairwise
             distinct leading variables, and the leading variable of a
             polynomial does not occur anywhere in other polynomials).

          3. The leading variables of the second sequence do not occur
             anywhere in the first sequence (these variables have been
             eliminated).

        When ``return_reductors==False``, only the first sequence is
        returned.

        EXAMPLE::

            sage: B.<a,b,c,d> = BooleanPolynomialRing()
            sage: F = Sequence([c + d + b + 1, a + c + d, a*b + c, b*c*d + c])
            sage: F.eliminate_linear_variables() # everything vanishes
            []
            sage: F.eliminate_linear_variables(maxlength=2)
            [b + c + d + 1, b*c + b*d + c, b*c*d + c]
            sage: F.eliminate_linear_variables(skip=lambda lm,tail: str(lm)=='a')
            [a + c + d, a*c + a*d + a + c, c*d + c]

        The list of reductors can be requested by setting 'return_reductors' to ``True``::

            sage: B.<a,b,c,d> = BooleanPolynomialRing()
            sage: F = Sequence([a + b + d, a + b + c])
            sage: F,R = F.eliminate_linear_variables(return_reductors=True)
            sage: F
            []
            sage: R
            [a + b + d, c + d]

        TESTS:

        The function should really dispose of linear equations (:trac:`13968`)::

            sage: R.<x,y,z> = BooleanPolynomialRing()
            sage: S = Sequence([x+y+z+1, y+z])
            sage: S.eliminate_linear_variables(return_reductors=True)
            ([], [x + 1, y + z])


        The function should take care of linear variables created by previous
        substitution of linear variables ::

            sage: R.<x,y,z> = BooleanPolynomialRing()
            sage: S = Sequence([x*y*z+x*y+z*y+x*z, x+y+z+1, x+y])
            sage: S.eliminate_linear_variables(return_reductors=True)
            ([], [x + y, z + 1])

        .. NOTE::

            This is called "massaging" in [CBJ07]_.

        REFERENCES:

        .. [CBJ07] Gregory V. Bard, and Nicolas T. Courtois, and Chris Jefferson.
           *Efficient Methods for Conversion and Solution of Sparse Systems of Low-Degree
           Multivariate Polynomials over GF(2) via SAT-Solvers*.
           Cryptology ePrint Archive: Report 2007/024. available at
           http://eprint.iacr.org/2007/024
        """
        from sage.rings.polynomial.pbori import BooleanPolynomialRing
        from polybori import gauss_on_polys
        from polybori.ll import eliminate,ll_encode,ll_red_nf_redsb

        R = self.ring()

        if not isinstance(R, BooleanPolynomialRing):
            raise NotImplementedError("Only BooleanPolynomialRing's are supported.")

        F = self
        reductors = []

        if skip is None and maxlength==Infinity:
            # faster solution based on polybori.ll.eliminate
            while True:
                (this_step_reductors, _, higher) = eliminate(F)
                if this_step_reductors == []:
                   break
                reductors.extend( this_step_reductors )
                F = higher
        else:
            # slower, more flexible solution
            if skip is None:
                skip = lambda lm,tail: False

            while True:
                linear = []
                higher = []

                for f in F:
                    if f.degree() == 1 and len(f) <= maxlength + 1:
                        flm = f.lex_lead()
                        if skip(flm, f-flm):
                            higher.append(f)
                            continue
                        linear.append(f)
                    else:
                        higher.append(f)

                if not linear:
                    break

                linear = gauss_on_polys(linear)
                rb = ll_encode(linear)
                reductors.extend(linear)

                F = []
                for f in higher:
                    f = ll_red_nf_redsb(f, rb)
                    if f != 0:
                        F.append(f)

        ret = PolynomialSequence(R, higher)
        if return_reductors:
            reduced_reductors = gauss_on_polys(reductors)
            return ret, PolynomialSequence(R, reduced_reductors)
        else:
            return ret

    def _groebner_strategy(self):
        """
        Return the Singular Groebner Strategy object.

        This object allows to compute normal forms efficiently, since
        all conversion overhead is avoided.

        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(GF(2))
            sage: F = Sequence([x*y + z, y + z + 1])
            sage: F._groebner_strategy()
            Groebner Strategy for ideal generated by 2 elements over
            Multivariate Polynomial Ring in x, y, z over Finite Field of size 2

            sage: P.<x,y,z> = BooleanPolynomialRing()
            sage: F = Sequence([x*y + z, y + z + 1])
            sage: F._groebner_strategy()
            <sage.rings.polynomial.pbori.GroebnerStrategy object at 0x...>
        """
        from sage.rings.polynomial.pbori import BooleanPolynomialRing
        R = self.ring()

        if not isinstance(R, BooleanPolynomialRing):
            from sage.libs.singular.groebner_strategy import GroebnerStrategy
            return GroebnerStrategy(self.ideal())
        else:
            from sage.rings.polynomial.pbori import GroebnerStrategy
            g = GroebnerStrategy(R)
            for p in self:
                g.add_as_you_wish(p)
            g.reduction_strategy.opt_red_tail=True
            return g

    def solve(self, algorithm='polybori', n=1,  eliminate_linear_variables=True, verbose=False, **kwds):
        r"""
        Find solutions of this boolean polynomial system.

        This function provide a unified interface to several algorithms
        dedicated to solving systems of boolean equations. Depending on
        the particular nature of the system, some might be much faster
        than some others.

        INPUT:

        * ``self`` - a sequence of boolean polynomials

        * ``algorithm`` - the method to use. Possible values are
          ``polybori``, ``sat`` and ``exhaustive_search``. (default:
          ``polybori``, since it is always available)

        * ``n`` - number of solutions to return. If ``n == +Infinity``
          then all solutions are returned. If `n < \infty` then `n`
          solutions are returned if the equations have at least `n`
          solutions. Otherwise, all the solutions are
          returned. (default: ``1``)

        * ``eliminate_linear_variables`` - whether to eliminate
          variables that appear linearly. This reduces the number of
          variables (makes solving faster a priori), but is likely to
          make the equations denser (may make solving slower depending
          on the method).

        * ``verbose`` - whether to display progress and (potentially)
          useful information while the computation runs. (default:
          ``False``)

        EXAMPLES:

        Without argument, a single arbitrary solution is returned::

            sage: R.<x,y,z> = BooleanPolynomialRing()
            sage: S = Sequence([x*y+z, y*z+x, x+y+z+1])
            sage: sol = S.solve(); sol                       # random
            [{y: 1, z: 0, x: 0}]

        We check that it is actually a solution::

            sage: S.subs( sol[0] )
            [0, 0, 0]

        We obtain all solutions::

            sage: sols = S.solve(n=Infinity); sols           # random
            [{x: 0, y: 1, z: 0}, {x: 1, y: 1, z: 1}]
            sage: map( lambda x: S.subs(x), sols)
            [[0, 0, 0], [0, 0, 0]]

        We can force the use of exhaustive search if the optional
        package ``FES`` is present::

            sage: sol = S.solve(algorithm='exhaustive_search'); sol  # random, optional - FES
            [{x: 1, y: 1, z: 1}]
            sage: S.subs( sol[0] )
            [0, 0, 0]

        And we may use SAT-solvers if they are available::

            sage: sol = S.solve(algorithm='sat'); sol                     # random, optional - CryptoMiniSat
            [{y: 1, z: 0, x: 0}]
            sage: S.subs( sol[0] )
            [0, 0, 0]

        TESTS:

        Make sure that variables not occuring in the equations are no problem::

            sage: R.<x,y,z,t> = BooleanPolynomialRing()
            sage: S = Sequence([x*y+z, y*z+x, x+y+z+1])
            sage: sols = S.solve(n=Infinity)
            sage: map( lambda x: S.subs(x), sols)
            [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]

        Not eliminating linear variables::

            sage: sols = S.solve(n=Infinity, eliminate_linear_variables=False)
            sage: map( lambda x: S.subs(x), sols)
            [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]

        A tricky case where the linear equations are insatisfiable::

            sage: R.<x,y,z> = BooleanPolynomialRing()
            sage: S = Sequence([x*y*z+x*y+z*y+x*z, x+y+z+1, x+y+z])
            sage: S.solve()
            []

        """
        from sage.rings.polynomial.pbori import BooleanPolynomialRing
        from sage.modules.free_module import VectorSpace

        S = self
        R_origin = R_solving = self.ring()
        reductors = []

        if eliminate_linear_variables:
            T, reductors = self.eliminate_linear_variables(return_reductors=True)
            if T.variables() != ():
                R_solving = BooleanPolynomialRing( T.nvariables(), map(str, list(T.variables())) )
            S = PolynomialSequence( R_solving, [ R_solving(f) for f in T] )

        if S != []:
            if algorithm == "exhaustive_search":
                if not is_package_installed('fes'):
                    raise ValueError('algorithm=exhaustive_search requires the optional library FES. Run "install_package(\'fes\')" to install it.')
                from sage.libs.fes import exhaustive_search
                solutions = exhaustive_search(S, max_sols=n, verbose=verbose, **kwds)

            elif algorithm == "polybori":
                I = S.ideal()
                if verbose:
                    I.groebner_basis(full_prot=True, **kwds)
                else:
                    I.groebner_basis(**kwds)
                solutions = I.variety()
                if len(solutions) >= n:
                    solutions = solutions[:n]

            elif algorithm == "sat":
                from sage.sat.boolean_polynomials import solve as solve_sat
                if verbose:
                    solutions = solve_sat(S, n=n, s_verbosity=1, **kwds)
                else:
                    solutions = solve_sat(S, n=n, **kwds)
            else:
                raise ValueError("unknown 'algorithm' value")
        else:
            solutions = []

        if S.variables() == ():
            solved_variables = set()
        else:
            solved_variables = { R_origin(x).lm() for x in R_solving.gens() }
        eliminated_variables = { f.lex_lead() for f in reductors }
        leftover_variables = { x.lm() for x in R_origin.gens() } - solved_variables - eliminated_variables

        if leftover_variables != set():
            partial_solutions = solutions
            solutions = []
            for sol in partial_solutions:
                for v in VectorSpace( GF(2), len(leftover_variables) ):
                    new_solution = sol.copy()
                    for var,val in zip(leftover_variables, v):
                        new_solution[ var ] = val
                    solutions.append( new_solution )

        for r in reductors:
            for sol in solutions:
                sol[ r.lm() ] = r.subs(sol).constant_coefficient()

        return solutions




class PolynomialSequence_gf2e(PolynomialSequence_generic):
    """
    PolynomialSequence over `\mathbb{F}_{2^e}`, i.e extensions over
    GF(2).
    """

    def weil_restriction(self):
        """
        Project this polynomial system to `\mathbb{F}_2`.

        That is, compute the Weil restriction of scalars for the
        variety corresponding to this polynomial system and express it
        as a polynomial system over `\mathbb{F}_2`.

        EXAMPLE::

            sage: k.<a> = GF(2^2)
            sage: P.<x,y> = PolynomialRing(k,2)
            sage: a = P.base_ring().gen()
            sage: F = Sequence([x*y + 1, a*x + 1], P)
            sage: F2 = F.weil_restriction()
            sage: F2
            [x1*y0 + x0*y1 + x1*y1,
             x0*y0 + x1*y1 + 1,
             x0 + x1,
             x1 + 1,
             x0^2 + x0,
             x1^2 + x1,
             y0^2 + y0,
             y1^2 + y1]

        Another bigger example for a small scale AES::

           sage: sr = mq.SR(1,1,1,4,gf2=False)
           sage: F,s = sr.polynomial_system(); F
           Polynomial Sequence with 40 Polynomials in 20 Variables
           sage: F2 = F.weil_restriction(); F2
           Polynomial Sequence with 240 Polynomials in 80 Variables
        """
        from sage.rings.ideal import FieldIdeal
        J = self.ideal().weil_restriction()
        J += FieldIdeal(J.ring())
        return PolynomialSequence(J)

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override("sage.crypto.mq.mpolynomialsystem","MPolynomialSystem_generic", PolynomialSequence_generic)
register_unpickle_override("sage.crypto.mq.mpolynomialsystem","MPolynomialRoundSystem_generic", PolynomialSequence_generic)

