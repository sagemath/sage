r"""
Period lattices of elliptic curves and related functions

Let `E` be an elliptic curve defined over a number field `K`
(including `\QQ`).  We attach a period lattice (a discrete rank 2
subgroup of `\CC`) to each embedding of `K` into `\CC`.

In the case of real embeddings, the lattice is stable under complex
conjugation and is called a real lattice.  These have two types:
rectangular, (the real curve has two connected components and positive
discriminant) or non-rectangular (one connected component, negative
discriminant).

The periods are computed to arbitrary precision using the AGM (Gauss's
Arithmetic-Geometric Mean).

EXAMPLES::

    sage: x = polygen(ZZ, 'x')
    sage: K.<a> = NumberField(x^3 - 2)                                                  # needs sage.rings.number_field
    sage: E = EllipticCurve([0,1,0,a,a])                                                # needs sage.rings.number_field

First we try a real embedding::

    sage: emb = K.embeddings(RealField())[0]                                            # needs sage.rings.number_field
    sage: L = E.period_lattice(emb); L                                                  # needs sage.rings.number_field
    Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x^2 + a*x + a
     over Number Field in a with defining polynomial x^3 - 2
     with respect to the embedding Ring morphism:
       From: Number Field in a with defining polynomial x^3 - 2
       To:   Algebraic Real Field
       Defn: a |--> 1.259921049894873?

The first basis period is real::

    sage: L.basis()                                                                     # needs sage.rings.number_field
    (3.81452977217855, 1.90726488608927 + 1.34047785962440*I)
    sage: L.is_real()                                                                   # needs sage.rings.number_field
    True

For a basis `\omega_1,\omega_2` normalised so that `\omega_1/\omega_2`
is in the fundamental region of the upper half-plane, use the method
:meth:`~sage.schemes.elliptic_curves.period_lattice.PeriodLattice_ell.normalised_basis`
instead::

    sage: L.normalised_basis()                                                          # needs sage.rings.number_field
    (1.90726488608927 - 1.34047785962440*I, -1.90726488608927 - 1.34047785962440*I)

Next a complex embedding::

    sage: emb = K.embeddings(ComplexField())[0]                                         # needs sage.rings.number_field
    sage: L = E.period_lattice(emb); L                                                  # needs sage.rings.number_field
    Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x^2 + a*x + a
     over Number Field in a with defining polynomial x^3 - 2
     with respect to the embedding Ring morphism:
       From: Number Field in a with defining polynomial x^3 - 2
       To:   Algebraic Field
       Defn: a |--> -0.6299605249474365? - 1.091123635971722?*I

In this case, the basis `\omega_1`, `\omega_2` is always normalised so
that `\tau = \omega_1/\omega_2` is in the fundamental region in the
upper half plane::

    sage: # needs sage.rings.number_field
    sage: w1, w2 = L.basis(); w1, w2
    (-1.37588604166076 - 2.58560946624443*I, -2.10339907847356 + 0.428378776460622*I)
    sage: L.is_real()
    False
    sage: tau = w1/w2; tau
    0.387694505032876 + 1.30821088214407*I
    sage: L.normalised_basis()
    (-1.37588604166076 - 2.58560946624443*I, -2.10339907847356 + 0.428378776460622*I)

We test that bug :issue:`8415` (caused by a PARI bug fixed in v2.3.5) is OK::

    sage: # needs sage.rings.number_field
    sage: E = EllipticCurve('37a')
    sage: K.<a> = QuadraticField(-7)
    sage: EK = E.change_ring(K)
    sage: EK.period_lattice(K.complex_embeddings()[0])
    Period lattice associated to Elliptic Curve defined by y^2 + y = x^3 + (-1)*x
     over Number Field in a with defining polynomial x^2 + 7
      with a = 2.645751311064591?*I
     with respect to the embedding Ring morphism:
       From: Number Field in a with defining polynomial x^2 + 7
             with a = 2.645751311064591?*I
       To:   Algebraic Field
       Defn: a |--> -2.645751311064591?*I

REFERENCES:

- [CT2013]_

AUTHORS:

- ?: initial version.

- John Cremona:

  - Adapted to handle real embeddings of number fields, September 2008.

  - Added basis_matrix function, November 2008

  - Added support for complex embeddings, May 2009.

  - Added complex elliptic logs, March 2010; enhanced, October 2010.
"""

import sage.rings.abc

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.modules.free_module import FreeModule_generic_pid
from sage.rings.complex_mpfr import ComplexField, ComplexNumber, ComplexField_class
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.qqbar import AA, QQbar
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RealField, RealField_class, RealNumber
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.structure.richcmp import richcmp_method, richcmp, richcmp_not_equal

lazy_import('sage.libs.pari', 'pari')
lazy_import('sage.rings.number_field.number_field', 'refine_embedding')


class PeriodLattice(FreeModule_generic_pid):
    """
    The class for the period lattice of an algebraic variety.
    """
    pass


@richcmp_method
class PeriodLattice_ell(PeriodLattice):
    r"""
    The class for the period lattice of an elliptic curve.

    Currently supported are elliptic curves defined over `\QQ`, and
    elliptic curves defined over a number field with a real or complex
    embedding, where the lattice constructed depends on that
    embedding.
    """

    def __init__(self, E, embedding=None):
        r"""
        Initialises the period lattice by storing the elliptic curve and the embedding.

        INPUT:

        - ``E`` -- an elliptic curve

        - ``embedding`` -- (default: ``None``) an embedding of the base
          field `K` of ``E`` into a real or complex field.  If ``None``:

          - use the built-in coercion to `\RR` for `K=\QQ`;

          - use the first embedding into `\RR` given by
            ``K.embeddings(RealField())``, if there are any;

          - use the first embedding into `\CC` given by
            ``K.embeddings(ComplexField())``, if `K` is totally complex.

        .. NOTE::

            No periods are computed on creation of the lattice; see the
            methods :meth:`basis`, :meth:`normalised_basis` and
            :meth:`real_period` for precision setting.

        EXAMPLES:

        This function is not normally called directly, but will be
        called by the period_lattice() function of classes
        ell_number_field and ell_rational_field::

            sage: from sage.schemes.elliptic_curves.period_lattice import PeriodLattice_ell
            sage: E = EllipticCurve('37a')
            sage: PeriodLattice_ell(E)
            Period lattice associated to
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = PeriodLattice_ell(E, emb); L
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x^2 + a*x + a
             over Number Field in a with defining polynomial x^3 - 2
             with respect to the embedding Ring morphism:
               From: Number Field in a with defining polynomial x^3 - 2
               To:   Algebraic Real Field
               Defn: a |--> 1.259921049894873?

            sage: emb = K.embeddings(ComplexField())[0]                                 # needs sage.rings.number_field
            sage: L = PeriodLattice_ell(E, emb); L                                      # needs sage.rings.number_field
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x^2 + a*x + a
             over Number Field in a with defining polynomial x^3 - 2
             with respect to the embedding Ring morphism:
               From: Number Field in a with defining polynomial x^3 - 2
               To:   Algebraic Field
               Defn: a |--> -0.6299605249474365? - 1.091123635971722?*I

        TESTS::

            sage: # needs sage.rings.number_field
            sage: from sage.schemes.elliptic_curves.period_lattice import PeriodLattice_ell
            sage: K.<a> = NumberField(x^3 - 2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = PeriodLattice_ell(E,emb)
            sage: L == loads(dumps(L))
            True

        Elliptic curve over imaginary number field without ``embedding`` specified::

            sage: E = EllipticCurve(QQ[I], [5, -3*I])
            sage: L = PeriodLattice_ell(E, embedding=None)
            sage: L.elliptic_logarithm(E(I+1, I+2))  # abs tol 1e-15
            -0.773376784700140 - 0.177736018028666*I
            sage: L.elliptic_exponential(_)  # abs tol 1e-15
            (1.00000000000000 - 1.00000000000000*I : 2.00000000000000 - 1.00000000000000*I : 1.00000000000000)
        """
        # First we cache the elliptic curve with this period lattice:

        self.E = E

        # Next we cache the embedding into QQbar or AA which extends
        # the given embedding:

        K = E.base_field()
        self.is_approximate = isinstance(K, (RealField_class, ComplexField_class))
        if embedding is None:
            if K in (AA, QQbar):
                embedding = K.hom(QQbar)
                real = K == AA
            elif self.is_approximate:
                embedding = K.hom(K)
                real = isinstance(K, RealField_class)
            else:
                embs = K.embeddings(AA)
                real = len(embs) > 0
                if not real:
                    embs = K.embeddings(QQbar)
                embedding = embs[0]
        else:
            embedding = refine_embedding(embedding, Infinity)
            real = embedding(K.gen()).imag().is_zero()

        self.embedding = embedding

        # Next we compute and cache (in self.real_flag) the type of
        # the lattice: +1 for real rectangular, -1 for real
        # non-rectangular, 0 for non-real:

        self.real_flag = 0
        if real:
            self.real_flag = +1
            if embedding(E.discriminant()) < 0:
                self.real_flag = -1

        # The following algebraic data associated to E and the
        # embedding is cached:
        #
        # Ebar: the curve E base-changed to QQbar (or AA)
        # f2: the 2-division polynomial of Ebar
        # ei: the roots e1, e2, e3 of f2, as elements of QQbar (or AA)
        #
        # The ei are used both for period computation and elliptic
        # logarithms.

        if self.is_approximate:
            self.f2 = self.E.two_division_polynomial()
        else:
            self.Ebar = self.E.change_ring(self.embedding)
            self.f2 = self.Ebar.two_division_polynomial()
        if self.real_flag == 1: # positive discriminant
            self._ei = self.f2.roots(K if self.is_approximate else AA,multiplicities=False)
            self._ei.sort()  # e1 < e2 < e3
            e1, e2, e3 = self._ei
        elif self.real_flag == -1: # negative discriminant
            self._ei = self.f2.roots(ComplexField(K.precision()) if self.is_approximate else QQbar, multiplicities=False)
            self._ei = sorted(self._ei, key=lambda z: z.imag())
            e1, e3, e2 = self._ei # so e3 is real
            if not self.is_approximate:
                e3 = AA(e3)
            self._ei = [e1, e2, e3]
        else:
            self._ei = self.f2.roots(ComplexField(K.precision()) if self.is_approximate else QQbar, multiplicities=False)
            e1, e2, e3 = self._ei

        # The quantities sqrt(e_i-e_j) are cached (as elements of
        # QQbar) to be used in period computations:

        self._abc = (e3-e1).sqrt(), (e3-e2).sqrt(), (e2-e1).sqrt()

        PeriodLattice.__init__(self, base_ring=ZZ, rank=2, degree=1, sparse=False)

    def __richcmp__(self, other, op):
        r"""
        Comparison function for period lattices.

        TESTS::

            sage: # needs sage.rings.number_field
            sage: from sage.schemes.elliptic_curves.period_lattice import PeriodLattice_ell
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: embs = K.embeddings(ComplexField())
            sage: L1, L2, L3 = [PeriodLattice_ell(E, e) for e in embs]
            sage: L1 < L2 < L3
            True
        """
        if not isinstance(other, PeriodLattice_ell):
            return NotImplemented

        lx = self.E
        rx = other.E
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        a = self.E.base_field().gen()
        return richcmp(self.embedding(a), other.embedding(a), op)

    def __repr__(self):
        """
        Return the string representation of this period lattice.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice()
            Period lattice associated to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb); L
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x^2 + a*x + a over Number Field in a
             with defining polynomial x^3 - 2 with respect to the embedding Ring morphism:
               From: Number Field in a with defining polynomial x^3 - 2
               To:   Algebraic Real Field
               Defn: a |--> 1.259921049894873?
        """
        K = self.E.base_field()
        if K in (QQ, AA, QQbar) or isinstance(K, (RealField_class, ComplexField_class)):
            return "Period lattice associated to %s" % (self.E)
        return "Period lattice associated to %s with respect to the embedding %s" % (self.E, self.embedding)

    def __call__(self, P, prec=None):
        r"""
        Return the elliptic logarithm of a point `P`.

        INPUT:

        - ``P`` -- a point on the elliptic curve associated with this period
          lattice

        - ``prec`` -- (default: ``None``) precision in bits (default
          precision if ``None``)

        OUTPUT:

        (complex number) The elliptic logarithm of the point `P` with
        respect to this period lattice.  If `E` is the elliptic curve
        and `\sigma:K\to\CC` the embedding, then the returned value `z`
        is such that `z\pmod{L}` maps to `\sigma(P)` under the
        standard Weierstrass isomorphism from `\CC/L` to `\sigma(E)`.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: L = E.period_lattice()
            sage: E.discriminant() > 0
            True
            sage: L.real_flag
            1
            sage: P = E([-1,1])
            sage: P.is_on_identity_component ()
            False
            sage: L(P, prec=96)  # abs tol 1e-27
            0.4793482501902193161295330101 + 0.985868850775824102211203849*I
            sage: Q = E([3,5])
            sage: Q.is_on_identity_component()
            True
            sage: L(Q, prec=96)
            1.931128271542559442488585220

        Note that this is actually the inverse of the Weierstrass isomorphism::

            sage: L.elliptic_exponential(L(Q))
            (3.00000000000000 : 5.00000000000000 : 1.00000000000000)

        An example with negative discriminant, and a torsion point::

            sage: E = EllipticCurve('11a1')
            sage: L = E.period_lattice()
            sage: E.discriminant() < 0
            True
            sage: L.real_flag
            -1
            sage: P = E([16,-61])
            sage: L(P)
            0.253841860855911
            sage: L.real_period() / L(P)
            5.00000000000000
        """
        return self.elliptic_logarithm(P,prec)

    @cached_method
    def basis(self, prec=None, algorithm='sage'):
        r"""
        Return a basis for this period lattice as a 2-tuple.

        INPUT:

        - ``prec`` -- (default: ``None``) precision in bits (default
          precision if ``None``)

        - ``algorithm`` -- string (default: ``'sage'``); choice of
          implementation (for real embeddings only) between ``'sage'``
          (native Sage implementation) or ``'pari'`` (use the PARI
          library: only available for real embeddings)

        OUTPUT:

        (tuple of Complex) `(\omega_1,\omega_2)` where the lattice is
        `\ZZ\omega_1 + \ZZ\omega_2`.  If the lattice is real then
        `\omega_1` is real and positive, `\Im(\omega_2)>0` and
        `\Re(\omega_2/\omega_1)` is either `0` (for rectangular
        lattices) or `\frac{1}{2}` (for non-rectangular lattices).
        Otherwise, `\omega_1/\omega_2` is in the fundamental region of
        the upper half-plane.  If the latter normalisation is required
        for real lattices, use the method :meth:`normalised_basis`
        instead.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().basis()
            (2.99345864623196, 2.45138938198679*I)

        This shows that the issue reported at :issue:`3954` is fixed::

            sage: E = EllipticCurve('37a')
            sage: b1 = E.period_lattice().basis(prec=30)
            sage: b2 = E.period_lattice().basis(prec=30)
            sage: b1 == b2
            True

        This shows that the issue reported at :issue:`4064` is fixed::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().basis(prec=30)[0].parent()
            Real Field with 30 bits of precision
            sage: E.period_lattice().basis(prec=100)[0].parent()
            Real Field with 100 bits of precision

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.basis(64)
            (3.81452977217854509, 1.90726488608927255 + 1.34047785962440202*I)

            sage: # needs sage.rings.number_field
            sage: emb = K.embeddings(ComplexField())[0]
            sage: L = E.period_lattice(emb)
            sage: w1, w2 = L.basis(); w1, w2
            (-1.37588604166076 - 2.58560946624443*I, -2.10339907847356 + 0.428378776460622*I)
            sage: L.is_real()
            False
            sage: tau = w1/w2; tau
            0.387694505032876 + 1.30821088214407*I
        """
        # We divide into two cases: (1) Q, or a number field with a
        # real embedding; (2) a number field with a complex embedding.
        # In each case the periods are computed by a different
        # internal function.

        if self.is_real():
            return self._compute_periods_real(prec=prec, algorithm=algorithm)
        else:
            return self._compute_periods_complex(prec=prec)

    @cached_method
    def gens(self, prec=None, algorithm='sage'):
        r"""
        Return a basis for this period lattice as a 2-tuple.

        This is an alias for
        :meth:`~sage.schemes.elliptic_curves.period_lattice.PeriodLattice_ell.basis`.
        See the docstring there for a more in-depth explanation and further
        examples.

        INPUT:

        - ``prec`` -- (default: ``None``) precision in bits (default
          precision if ``None``)

        - ``algorithm`` -- string (default: ``'sage'``); choice of
          implementation (for real embeddings only) between ``'sage'``
          (native Sage implementation) or ``'pari'`` (use the PARI
          library: only available for real embeddings)

        OUTPUT:

        (tuple of Complex) `(\omega_1,\omega_2)` where the lattice is
        `\ZZ\omega_1 + \ZZ\omega_2`.  If the lattice is real then
        `\omega_1` is real and positive, `\Im(\omega_2)>0` and
        `\Re(\omega_2/\omega_1)` is either `0` (for rectangular
        lattices) or `\frac{1}{2}` (for non-rectangular lattices).
        Otherwise, `\omega_1/\omega_2` is in the fundamental region of
        the upper half-plane.  If the latter normalisation is required
        for real lattices, use the method :meth:`normalised_basis`
        instead.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().gens()
            (2.99345864623196, 2.45138938198679*I)

            sage: E.period_lattice().gens(prec=100)
            (2.9934586462319596298320099794, 2.4513893819867900608542248319*I)
        """
        return tuple(self.basis(prec=prec, algorithm=algorithm))

    @cached_method
    def normalised_basis(self, prec=None, algorithm='sage'):
        r"""
        Return a normalised basis for this period lattice as a 2-tuple.

        INPUT:

        - ``prec`` -- (default: ``None``) precision in bits (default
          precision if ``None``)

        - ``algorithm`` -- string (default: ``'sage'``); choice of
          implementation (for real embeddings only) between ``'sage'``
          (native Sage implementation) or ``'pari'`` (use the PARI
          library: only available for real embeddings)

        OUTPUT:

        (tuple of Complex) `(\omega_1,\omega_2)` where the lattice has
        the form `\ZZ\omega_1 + \ZZ\omega_2`.  The basis is normalised
        so that `\omega_1/\omega_2` is in the fundamental region of
        the upper half-plane.  For an alternative normalisation for
        real lattices (with the first period real), use the method
        :meth:`basis` instead.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().normalised_basis()
            (2.99345864623196, -2.45138938198679*I)

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.normalised_basis(64)
            (1.90726488608927255 - 1.34047785962440202*I,
            -1.90726488608927255 - 1.34047785962440202*I)

            sage: # needs sage.rings.number_field
            sage: emb = K.embeddings(ComplexField())[0]
            sage: L = E.period_lattice(emb)
            sage: w1, w2 = L.normalised_basis(); w1, w2
            (-1.37588604166076 - 2.58560946624443*I,
             -2.10339907847356 + 0.428378776460622*I)
            sage: L.is_real()
            False
            sage: tau = w1/w2; tau
            0.387694505032876 + 1.30821088214407*I
        """
        w1, w2 = self.basis(prec=prec, algorithm=algorithm)
        periods, _ = normalise_periods(w1, w2)
        return periods

    @cached_method
    def tau(self, prec=None, algorithm='sage'):
        r"""
        Return the upper half-plane parameter in the fundamental region.

        INPUT:

        - ``prec`` -- (default: ``None``) precision in bits (default
          precision if ``None``)

        - ``algorithm`` -- string (default: ``'sage'``); choice of
          implementation (for real embeddings only) between 'sage'
          (native Sage implementation) or 'pari' (use the PARI
          library: only available for real embeddings)

        OUTPUT:

        (Complex) `\tau = \omega_1/\omega_2` where the lattice has the
        form `\ZZ\omega_1 + \ZZ\omega_2`, normalised so that `\tau =
        \omega_1/\omega_2` is in the fundamental region of the upper
        half-plane.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: L = E.period_lattice()
            sage: L.tau()
            1.22112736076463*I

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: tau = L.tau(); tau
            -0.338718341018919 + 0.940887817679340*I
            sage: tau.abs()
            1.00000000000000
            sage: -0.5 <= tau.real() <= 0.5
            True

            sage: # needs sage.rings.number_field
            sage: emb = K.embeddings(ComplexField())[0]
            sage: L = E.period_lattice(emb)
            sage: tau = L.tau(); tau
            0.387694505032876 + 1.30821088214407*I
            sage: tau.abs()
            1.36444961115933
            sage: -0.5 <= tau.real() <= 0.5
            True
        """
        w1, w2 = self.normalised_basis(prec=prec, algorithm=algorithm)
        return w1/w2

    @cached_method
    def _compute_default_prec(self):
        r"""
        Internal function to compute the default precision to be used if nothing is passed in.
        """
        return self.E.base_field().precision() if self.is_approximate else RealField().precision()

    @cached_method
    def _compute_periods_real(self, prec=None, algorithm='sage'):
        r"""
        Internal function to compute the periods (real embedding case).

        INPUT:

        - ``prec`` -- integer or ``None`` (default); floating point
          precision (in bits). If ``None``, use the default precision.

        - ``algorithm`` string (default: ``'sage'``); choice of implementation between

          - ``'pari'`` -- use the PARI library

          - ``'sage'`` -- use a native Sage implementation (with the same underlying algorithm)

        OUTPUT:

        (tuple of Complex) `(\omega_1,\omega_2)` where the lattice has
        the form `\ZZ\omega_1 + \ZZ\omega_2`, `\omega_1` is real and
        `\omega_2/\omega_1` has real part either `0` or `frac{1}{2}`.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: embs = K.embeddings(CC)
            sage: Ls = [E.period_lattice(e) for e in embs]
            sage: [L.is_real() for L in Ls]
            [False, False, True]
            sage: Ls[2]._compute_periods_real(100)
            (3.8145297721785450936365098936,
            1.9072648860892725468182549468 + 1.3404778596244020196600112394*I)
            sage: Ls[2]._compute_periods_real(100, algorithm='pari')
            (3.8145297721785450936365098936,
            1.9072648860892725468182549468 - 1.3404778596244020196600112394*I)
        """
        if prec is None:
            prec = self._compute_default_prec()
        R = RealField(prec)
        C = ComplexField(prec)

        if algorithm == 'pari':
            ainvs = self.E.a_invariants()
            if self.E.base_field() is not QQ and not self.is_approximate:
                ainvs = [C(self.embedding(ai)).real() for ai in ainvs]

            # The precision for omega() is determined by ellinit()
            E_pari = pari.ellinit(ainvs, precision=prec)
            w1, w2 = E_pari.omega()
            return R(w1), C(w2)

        if algorithm != 'sage':
            raise ValueError("invalid value of 'algorithm' parameter")

        pi = R.pi()
        # Up to now everything has been exact in AA or QQbar (unless self.is_approximate),
        # but now we must go transcendental.  Only now is the desired precision used!
        if self.real_flag == 1: # positive discriminant
            a, b, c = (R(x) for x in self._abc)
            w1 = R(pi/a.agm(b))   # least real period
            w2 = C(0,pi/a.agm(c)) # least pure imaginary period
        else:
            a = C(self._abc[0])
            x, y, r = a.real().abs(), a.imag().abs(), a.abs()
            w1 = R(pi/r.agm(x)) # least real period
            w2 = R(pi/r.agm(y)) # least pure imaginary period /i
            w2 = C(w1,w2)/2

        return (w1,w2)

    @cached_method
    def _compute_periods_complex(self, prec=None, normalise=True):
        r"""
        Internal function to compute the periods (complex embedding case).

        INPUT:

        - ``prec`` -- integer or ``None`` (default); floating point precision
          (in bits); if ``None``, use the default precision

        - ``normalise`` -- boolean (default: ``True``); whether to normalise the
          basis after computation

        OUTPUT:

        (tuple of Complex) `(\omega_1,\omega_2)` where the lattice has
        the form `\ZZ\omega_1 + \ZZ\omega_2`.  If `normalise` is
        ``True``, the basis is normalised so that `(\omega_1/\omega_2)`
        is in the fundamental region of the upper half plane.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: embs = K.embeddings(CC)
            sage: Ls = [E.period_lattice(e) for e in embs]
            sage: [L.is_real() for L in Ls]
            [False, False, True]
            sage: L = Ls[0]
            sage: w1,w2 = L._compute_periods_complex(100); w1,w2
            (-1.3758860416607626645495991458 - 2.5856094662444337042877901304*I,
             -2.1033990784735587243397865076 + 0.42837877646062187766760569686*I)
            sage: tau = w1/w2; tau
            0.38769450503287609349437509561 + 1.3082108821440725664008561928*I
            sage: tau.real()
            0.38769450503287609349437509561
            sage: tau.abs()
            1.3644496111593345713923386773

        Without normalisation::

            sage: # needs sage.rings.number_field
            sage: w1,w2 = L._compute_periods_complex(normalise=False); w1,w2
            (2.10339907847356 - 0.428378776460622*I, 0.727513036812796 - 3.01398824270506*I)
            sage: tau = w1/w2; tau
            0.293483964608883 + 0.627038168678760*I
            sage: tau.real()
            0.293483964608883
            sage: tau.abs()          # > 1
            0.692321964451917
        """
        if prec is None:
            prec = self._compute_default_prec()
        C = ComplexField(prec)

        # Up to now everything has been exact in AA or QQbar (unless self.is_approximate),
        # but now we must go transcendental.  Only now is the desired precision used!
        pi = C.pi()
        a, b, c = (C(x) for x in self._abc)
        if (a+b).abs() < (a-b).abs():
            b = -b
        if (a+c).abs() < (a-c).abs():
            c = -c
        w1 = pi/a.agm(b)
        w2 = pi*C.gen()/a.agm(c)
        if (w1/w2).imag() < 0:
            w2 = -w2
        if normalise:
            w1w2, mat = normalise_periods(w1,w2)
            return w1w2
        return (w1,w2)

    def is_real(self):
        r"""
        Return ``True`` if this period lattice is real.

        EXAMPLES::

            sage: f = EllipticCurve('11a')
            sage: f.period_lattice().is_real()
            True

        ::

            sage: # needs sage.rings.number_field
            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve(K, [0,0,0,i,2*i])
            sage: emb = K.embeddings(ComplexField())[0]
            sage: L = E.period_lattice(emb)
            sage: L.is_real()
            False

        ::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)                                          # needs sage.rings.number_field
            sage: E = EllipticCurve([0,1,0,a,a])                                        # needs sage.rings.number_field
            sage: [E.period_lattice(emb).is_real() for emb in K.embeddings(CC)]         # needs sage.rings.number_field
            [False, False, True]

        ALGORITHM:

        The lattice is real if it is associated to a real embedding;
        such lattices are stable under conjugation.
        """
        return self.real_flag != 0

    def is_rectangular(self) -> bool:
        r"""
        Return ``True`` if this period lattice is rectangular.

        .. NOTE::

            Only defined for real lattices; a :exc:`RuntimeError`
            is raised for non-real lattices.

        EXAMPLES::

            sage: f = EllipticCurve('11a')
            sage: f.period_lattice().basis()
            (1.26920930427955, 0.634604652139777 + 1.45881661693850*I)
            sage: f.period_lattice().is_rectangular()
            False

        ::

            sage: f = EllipticCurve('37b')
            sage: f.period_lattice().basis()
            (1.08852159290423, 1.76761067023379*I)
            sage: f.period_lattice().is_rectangular()
            True

        ALGORITHM:

        The period lattice is rectangular precisely if the
        discriminant of the Weierstrass equation is positive, or
        equivalently if the number of real components is 2.
        """
        if self.is_real():
            return self.real_flag == +1
        raise RuntimeError("Not defined for non-real lattices.")

    def real_period(self, prec=None, algorithm='sage'):
        """
        Return the real period of this period lattice.

        INPUT:

        - ``prec`` -- integer or ``None`` (default); real precision in
          bits (default real precision if ``None``)

        - ``algorithm`` -- string (default: ``'sage'``); choice of
          implementation (for real embeddings only) between ``'sage'``
          (native Sage implementation) or ``'pari'`` (use the PARI
          library: only available for real embeddings)

        .. NOTE::

            Only defined for real lattices; a :exc:`RuntimeError`
            is raised for non-real lattices.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().real_period()
            2.99345864623196

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.real_period(64)
            3.81452977217854509
        """
        if self.is_real():
            return self.basis(prec,algorithm)[0]
        raise RuntimeError("Not defined for non-real lattices.")

    def omega(self, prec=None, bsd_normalise=False):
        r"""
        Return the real or complex volume of this period lattice.

        INPUT:

        - ``prec`` -- integer or ``None`` (default); real precision in
          bits (default real precision if ``None``)

        - ``bsd_normalise`` -- boolean (default: ``False``); flag to use
          BSD normalisation in the complex case

        OUTPUT:

        (real) For real lattices, this is the real period times the
        number of connected components.  For non-real lattices it is
        the complex area, or double the area if ``bsd_normalise`` is
        ``True``.

        .. NOTE::

            If the curve is given by a *global minimal* Weierstrass
            equation, then with ``bsd_normalise`` = ``True``, this
            gives the correct period in the BSD conjecture: the product
            of this quantity over all embeddings appears in the BSD
            formula. In general a correction factor is required to make
            allowance for the model.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().omega()
            5.98691729246392

        This is not a minimal model::

            sage: E = EllipticCurve([0, -432*6^2])
            sage: E.period_lattice().omega()
            0.486109385710056

        If you were to plug the above omega into the BSD conjecture,
        you would get an incorrect value, out by a factor of 2.  The
        following works though::

            sage: F = E.minimal_model()
            sage: F.period_lattice().omega()
            0.972218771420113

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.omega(64)
            3.81452977217854509

        A complex example (taken from J.E.Cremona and E.Whitley,
        *Periods of cusp forms and elliptic curves over imaginary
        quadratic fields*, Mathematics of Computation 62 No. 205
        (1994), 407-429).  See :issue:`29645` and :issue:`29782`::

            sage: # needs sage.rings.number_field
            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve([0,1-i,i,-i,0])
            sage: L = E.period_lattice(K.embeddings(CC)[0])
            sage: L.omega()
            8.80694160502647
            sage: L.omega(prec=200)
            8.8069416050264741493250743632295462227858630765392114070032
            sage: L.omega(bsd_normalise=True)
            17.6138832100529
        """
        if self.is_real():
            n_components = 2 if self.real_flag == 1 else 1
            return self.real_period(prec) * n_components
        else:
            bsd_factor = 2 if bsd_normalise else 1
            return self.complex_area(prec) * bsd_factor

    @cached_method
    def basis_matrix(self, prec=None, normalised=False):
        r"""
        Return the basis matrix of this period lattice.

        INPUT:

        - ``prec`` -- integer or ``None`` (default); real precision in
          bits (default real precision if ``None``)

        - ``normalised`` -- boolean (default: ``False``); if ``True`` and the
          embedding is real, use the normalised basis (see
          :meth:`normalised_basis`) instead of the default

        OUTPUT:

        A `2\times 2` real matrix whose rows are the lattice basis vectors,
        after identifying `\CC` with `\RR^2`.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().basis_matrix()
            [ 2.99345864623196 0.000000000000000]
            [0.000000000000000  2.45138938198679]

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.basis_matrix(64)
            [ 3.81452977217854509 0.000000000000000000]
            [ 1.90726488608927255  1.34047785962440202]

        See :issue:`4388`::

            sage: L = EllipticCurve('11a1').period_lattice()
            sage: L.basis_matrix()
            [ 1.26920930427955 0.000000000000000]
            [0.634604652139777  1.45881661693850]
            sage: L.basis_matrix(normalised=True)
            [0.634604652139777 -1.45881661693850]
            [-1.26920930427955 0.000000000000000]

        ::

            sage: L = EllipticCurve('389a1').period_lattice()
            sage: L.basis_matrix()
            [ 2.49021256085505 0.000000000000000]
            [0.000000000000000  1.97173770155165]
            sage: L.basis_matrix(normalised=True)
            [ 2.49021256085505 0.000000000000000]
            [0.000000000000000 -1.97173770155165]
        """
        from sage.matrix.constructor import Matrix

        if normalised:
            return Matrix([list(w) for w in self.normalised_basis(prec)])

        w1,w2 = self.basis(prec)
        if self.is_real():
            return Matrix([[w1,0],list(w2)])
        else:
            return Matrix([list(w) for w in (w1,w2)])

    def complex_area(self, prec=None):
        """
        Return the area of a fundamental domain for the period lattice
        of the elliptic curve.

        INPUT:

        - ``prec`` -- integer or ``None`` (default); real precision in
          bits (default real precision if ``None``)

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().complex_area()
            7.33813274078958

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: embs = K.embeddings(ComplexField())
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: [E.period_lattice(emb).is_real() for emb in K.embeddings(CC)]
            [False, False, True]
            sage: [E.period_lattice(emb).complex_area() for emb in embs]
            [6.02796894766694, 6.02796894766694, 5.11329270448345]
        """
        w1,w2 = self.basis(prec)
        return (w1*w2.conjugate()).imag().abs()

    def sigma(self, z, prec=None, flag=0):
        r"""
        Return the value of the Weierstrass sigma function for this elliptic curve  period lattice.

        INPUT:

        - ``z`` -- a complex number

        - ``prec`` -- (default: ``None``) real precision in bits
            (default real precision if ``None``)

        - ``flag`` --

          0: (default) ???;

          1: computes an arbitrary determination of log(sigma(z))

          2, 3: same using the product expansion instead of theta series. ???

        .. NOTE::

            The reason for the ???'s above, is that the PARI
            documentation for ellsigma is very vague.  Also this is
            only implemented for curves defined over `\QQ`.

        .. TODO::

            This function does not use any of the PeriodLattice functions
            and so should be moved to ell_rational_field.

        EXAMPLES::

            sage: EllipticCurve('389a1').period_lattice().sigma(CC(2,1))
            2.60912163570108 - 0.200865080824587*I
        """
        if prec is None:
            prec = self._compute_default_prec()
        try:
            return self.E.pari_curve().ellsigma(z, flag, precision=prec)
        except AttributeError:
            raise NotImplementedError("sigma function not yet implemented for period lattices of curves not defined over Q")

    def curve(self):
        r"""
        Return the elliptic curve associated with this period lattice.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: L = E.period_lattice()
            sage: L.curve() is E
            True

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(K.embeddings(RealField())[0])
            sage: L.curve() is E
            True

            sage: L = E.period_lattice(K.embeddings(ComplexField())[0])                 # needs sage.rings.number_field
            sage: L.curve() is E                                                        # needs sage.rings.number_field
            True
        """
        return self.E

    def ei(self):
        r"""
        Return the x-coordinates of the 2-division points of the elliptic curve associated
        with this period lattice, as elements of ``QQbar``.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: L = E.period_lattice()
            sage: L.ei()
            [-1.107159871688768?, 0.2695944364054446?, 0.8375654352833230?]

        In the following example, we should have one purely real 2-division point coordinate,
        and two conjugate purely imaginary coordinates.

        ::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(K.embeddings(RealField())[0])
            sage: x1,x2,x3 = L.ei()
            sage: abs(x1.real()) + abs(x2.real()) < 1e-14
            True
            sage: x1.imag(), x2.imag(), x3
            (-1.122462048309373?, 1.122462048309373?, -1.000000000000000?)

        ::

            sage: L = E.period_lattice(K.embeddings(ComplexField())[0])                 # needs sage.rings.number_field
            sage: L.ei()                                                                # needs sage.rings.number_field
            [-1.000000000000000? + 0.?e-1...*I,
             -0.9720806486198328? - 0.561231024154687?*I,
             0.9720806486198328? + 0.561231024154687?*I]
        """
        return self._ei

    def coordinates(self, z, rounding=None):
        r"""
        Return the coordinates of a complex number w.r.t. the lattice basis.

        INPUT:

        - ``z`` -- complex number

        - ``rounding`` -- (default: ``None``) whether and how to round the
          output (see below)

        OUTPUT:

        When ``rounding`` is ``None`` (the default), returns a tuple
        of reals `x`, `y` such that `z=xw_1+yw_2` where `w_1`, `w_2`
        are a basis for the lattice (normalised in the case of complex
        embeddings).

        When ``rounding`` is ``'round'``, returns a tuple of integers `n_1`,
        `n_2` which are the closest integers to the `x`, `y` defined
        above.  If `z` is in the lattice these are the coordinates of
        `z` with respect to the lattice basis.

        When ``rounding`` is ``'floor'``, returns a tuple of integers
        `n_1`, `n_2` which are the integer parts to the `x`, `y`
        defined above. These are used in :meth:`.reduce`

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: L = E.period_lattice()
            sage: w1, w2 = L.basis(prec=100)
            sage: P = E([-1,1])
            sage: zP = P.elliptic_logarithm(precision=100); zP  # abs tol 1e-28
            0.47934825019021931612953301006 + 0.98586885077582410221120384908*I
            sage: L.coordinates(zP)  # abs tol 1e-28
            (0.19249290511394227352563996419, 0.50000000000000000000000000000)
            sage: sum([x*w for x, w in zip(L.coordinates(zP), L.basis(prec=100))])  # abs tol 1e-28
            0.47934825019021931612953301006 + 0.98586885077582410221120384908*I

            sage: L.coordinates(12*w1 + 23*w2)  # abs tol 1e-28
            (12.000000000000000000000000000, 23.000000000000000000000000000)
            sage: L.coordinates(12*w1 + 23*w2, rounding='floor')
            (11, 22)
            sage: L.coordinates(12*w1 + 23*w2, rounding='round')
            (12, 23)
        """
        C = z.parent()
        if isinstance(C, sage.rings.abc.RealField):
            C = ComplexField(C.precision())
            z = C(z)
        else:
            if isinstance(C, sage.rings.abc.ComplexField):
                pass
            else:
                try:
                    C = ComplexField()
                    z = C(z)
                except TypeError:
                    raise TypeError("%s is not a complex number" % z)
        prec = C.precision()
        from sage.matrix.constructor import Matrix
        from sage.modules.free_module_element import vector
        if self.real_flag:
            w1, w2 = self.basis(prec)
            M = Matrix([[w1, 0], list(w2)])**(-1)
        else:
            w1, w2 = self.normalised_basis(prec)
            M = Matrix([list(w1), list(w2)])**(-1)
        u, v = vector(z) * M
        # Now z = u*w1+v*w2
        if rounding == 'round':
            return u.round(), v.round()
        if rounding == 'floor':
            return u.floor(), v.floor()
        return u, v

    def reduce(self, z):
        r"""
        Reduce a complex number modulo the lattice.

        INPUT:

        - ``z`` -- complex number

        OUTPUT:

        (complex) the reduction of `z` modulo the lattice, lying in
        the fundamental period parallelogram with respect to the
        lattice basis.  For curves defined over the reals (i.e. real
        embeddings) the output will be real when possible.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: L = E.period_lattice()
            sage: w1, w2 = L.basis(prec=100)
            sage: P = E([-1,1])
            sage: zP = P.elliptic_logarithm(precision=100); zP  # abs tol 1e-28
            0.47934825019021931612953301006 + 0.98586885077582410221120384908*I
            sage: z = zP + 10*w1 - 20*w2; z  # abs tol 1e-28
            25.381473858740770069343110929 - 38.448885180257139986236950114*I
            sage: L.reduce(z)  # abs tol 1e-28
            0.47934825019021931612953301006 + 0.98586885077582410221120384908*I
            sage: L.elliptic_logarithm(2*P)  # abs tol 1e-15
            0.958696500380439
            sage: L.reduce(L.elliptic_logarithm(2*P))  # abs tol 1e-15
            0.958696500380439
            sage: L.reduce(L.elliptic_logarithm(2*P) + 10*w1 - 20*w2)  # abs tol 1e-15
            0.958696500380444
        """
        C = z.parent()
        if isinstance(C, sage.rings.abc.RealField):
            C = ComplexField(C.precision())
            z = C(z)
        elif isinstance(C, sage.rings.abc.ComplexField):
            pass
        else:
            try:
                C = ComplexField()
                z = C(z)
            except TypeError:
                raise TypeError("%s is not a complex number" % z)
        prec = C.precision()
        if self.real_flag:
            w1, w2 = self.basis(prec)  # w1 real
        else:
            w1, w2 = self.normalised_basis(prec)
        u, v = self.coordinates(z, rounding='floor')
        z = z-u*w1-v*w2

        # Final adjustments for the real case.

        # NB We assume here that when the embedding is real then the
        # point is also real!

        if self.real_flag == 0:
            return z
        if self.real_flag == -1:
            k = (z.imag() / w2.imag()).round()
            z = z-k*w2
            return C(z.real(), 0)

        if ((2*z.imag()/w2.imag()).round()) % 2:
            return C(z.real(), w2.imag() / 2)
        else:
            return C(z.real(), 0)

    def e_log_RC(self, xP, yP, prec=None, reduce=True):
        r"""
        Return the elliptic logarithm of a real or complex point.

        INPUT:

        - ``xP``, ``yP`` (real or complex) -- Coordinates of a point on
          the embedded elliptic curve associated with this period
          lattice.

        - ``prec`` -- (default: ``None``) real precision in bits
          (default real precision if ``None``)

        - ``reduce`` -- boolean (default: ``True``); if ``True``, the result
          is reduced with respect to the period lattice basis

        OUTPUT:

        (complex number) The elliptic logarithm of the point `(xP,yP)`
        with respect to this period lattice.  If `E` is the elliptic
        curve and `\sigma:K\to\CC` the embedding, the returned
        value `z` is such that `z\pmod{L}` maps to `(xP,yP)=\sigma(P)`
        under the standard Weierstrass isomorphism from `\CC/L` to
        `\sigma(E)`.  If ``reduce`` is ``True``, the output is reduced
        so that it is in the fundamental period parallelogram with
        respect to the normalised lattice basis.

        ALGORITHM:

        Uses the complex AGM.  See [CT2013]_ for details.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: L = E.period_lattice()
            sage: P = E([-1,1])
            sage: xP, yP = [RR(c) for c in P.xy()]

        The elliptic log from the real coordinates::

            sage: L.e_log_RC(xP, yP)  # abs tol 1e-15
            0.479348250190219 + 0.985868850775824*I

        The same elliptic log from the algebraic point::

            sage: L(P)  # abs tol 1e-15
            0.479348250190219 + 0.985868850775824*I

        A number field example::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 2)
            sage: E = EllipticCurve([0,0,0,0,a])
            sage: v = K.real_places()[0]
            sage: L = E.period_lattice(v)
            sage: P = E.lift_x(1/3*a^2 + a + 5/3)
            sage: L(P)  # abs tol 1e-15
            3.51086196882538
            sage: xP, yP = [v(c) for c in P.xy()]
            sage: L.e_log_RC(xP, yP)  # abs tol 1e-15
            3.51086196882538

        Elliptic logs of real points which do not come from algebraic
        points::

            sage: # needs sage.rings.number_field
            sage: ER = EllipticCurve([v(ai) for ai in E.a_invariants()])
            sage: P = ER.lift_x(12.34)
            sage: xP, yP = P.xy()
            sage: xP, yP  # abs tol 1e-15
            (12.3400000000000, -43.3628968710567)
            sage: L.e_log_RC(xP, yP)  # abs tol 1e-15
            0.284656841192041
            sage: xP, yP = ER.lift_x(0).xy()  # abs tol 1e-15
            sage: L.e_log_RC(xP, yP)
            1.34921304541057

        Elliptic logs of complex points::

            sage: # needs sage.rings.number_field
            sage: v = K.complex_embeddings()[0]
            sage: L = E.period_lattice(v)
            sage: P = E.lift_x(1/3*a^2 + a + 5/3)
            sage: L(P)  # abs tol 1e-15
            1.68207104397706 - 1.87873661686704*I
            sage: xP, yP = [v(c) for c in P.xy()]
            sage: L.e_log_RC(xP, yP)  # abs tol 1e-15
            1.68207104397706 - 1.87873661686704*I
            sage: EC = EllipticCurve([v(ai) for ai in E.a_invariants()])
            sage: xP, yP = EC.lift_x(0).xy()
            sage: L.e_log_RC(xP, yP)  # abs tol 1e-15
            2.06711431204080 - 1.73451485683471*I
        """
        if prec is None:
            prec = self._compute_default_prec()
        # Note: using log2(prec) + 3 guard bits is usually enough.
        # To avoid computing a logarithm, we use 40 guard bits which
        # should be largely enough in practice.
        prec2 = prec + 40

        R = RealField(prec2)
        C = ComplexField(prec2)
        e1,e2,e3 = self._ei
        a1,a2,a3 = (self.embedding(a) for a in self.E.ainvs()[:3])

        wP = 2*yP+a1*xP+a3

        # We treat the case of 2-torsion points separately.  (Note
        # that Cohen's algorithm does not handle these properly.)

        if wP.is_zero():  # 2-torsion treated separately
            w1,w2 = self._compute_periods_complex(prec,normalise=False)
            if xP == e1:
                z = w2/2
            else:
                if xP == e3:
                    z = w1/2
                else:
                    z = (w1+w2)/2
            if reduce:
                z = self.reduce(z)
            return z

        # NB The first block of code works fine for real embeddings as
        # well as complex embeddings.  The special code for real
        # embeddings uses only real arithmetic in the iteration, and is
        # based on Cremona and Thongjunthug.

        # An older version, based on Cohen's Algorithm 7.4.8 also uses
        # only real arithmetic, and gives different normalisations,
        # but also causes problems (see #10026).  It is left in but
        # commented out below.

        if self.real_flag == 0:  # complex case

            a = C((e1-e3).sqrt())
            b = C((e1-e2).sqrt())
            if (a+b).abs() < (a-b).abs():
                b = -b
            r = C(((xP-e3)/(xP-e2)).sqrt())
            if r.real() < 0:
                r = -r
            t = -C(wP)/(2*r*(xP-e2))
            # eps controls the end of the loop. Since we aim at a target
            # precision of prec bits, eps = 2^(-prec) is enough.
            eps = R(1) >> prec
            while True:
                s = b*r+a
                a, b = (a+b)/2, (a*b).sqrt()
                if (a+b).abs() < (a-b).abs():
                    b = -b
                r = (a*(r+1)/s).sqrt()
                if (r.abs()-1).abs() < eps:
                    break
                if r.real() < 0:
                    r = -r
                t *= r
            z = ((a/t).arctan())/a
            z = ComplexField(prec)(z)
            if reduce:
                z = self.reduce(z)
            return z

        if self.real_flag == -1: # real, connected case
            z = C(self._abc[0]) # sqrt(e3-e1)
            a, y, b = z.real(), z.imag(), z.abs()
            uv = (xP-e1).sqrt()
            u, v = uv.real().abs(), uv.imag().abs()
            r = (u*a/(u*a+v*y)).sqrt()
            t = -r*R(wP)/(2*(u**2+v**2))
            on_egg = False
        else:                  # real, disconnected case
            a = R(e3-e1).sqrt()
            b = R(e3-e2).sqrt()
            if (a+b).abs() < (a-b).abs():
                b = -b
            on_egg = (xP < e3)
            if on_egg:
                r = a/R(e3-xP).sqrt()
                t = r*R(wP)/(2*R(xP-e1))
            else:
                r = R((xP-e1)/(xP-e2)).sqrt()
                t = -R(wP)/(2*r*R(xP-e2))

        # eps controls the end of the loop. Since we aim at a target
        # precision of prec bits, eps = 2^(-prec) is enough.
        eps = R(1) >> prec
        while True:
            s = b*r+a
            a, b = (a+b)/2, (a*b).sqrt()
            r = (a*(r+1)/s).sqrt()
            if (r-1).abs() < eps:
                break
            t *= r
        z = ((a / t).arctan()) / a
        if on_egg:
            w1, w2 = self._compute_periods_real(prec)
            z += w2 / 2
        z = ComplexField(prec)(z)
        if reduce:
            z = self.reduce(z)
        return z

    def elliptic_logarithm(self, P, prec=None, reduce=True):
        r"""
        Return the elliptic logarithm of a point.

        INPUT:

        - ``P`` -- point on the elliptic curve associated with this period
          lattice

        - ``prec`` -- (default: ``None``) real precision in bits
          (default real precision if ``None``)

        - ``reduce`` -- boolean (default: ``True``); if ``True``, the result
          is reduced with respect to the period lattice basis

        OUTPUT:

        (complex number) The elliptic logarithm of the point `P` with
        respect to this period lattice.  If `E` is the elliptic curve
        and `\sigma:K\to\CC` the embedding, the returned value `z`
        is such that `z\pmod{L}` maps to `\sigma(P)` under the
        standard Weierstrass isomorphism from `\CC/L` to `\sigma(E)`.
        If ``reduce`` is ``True``, the output is reduced so that it is
        in the fundamental period parallelogram with respect to the
        normalised lattice basis.

        ALGORITHM:

        Uses the complex AGM.  See [CT2013]_ for details.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: L = E.period_lattice()
            sage: E.discriminant() > 0
            True
            sage: L.real_flag
            1
            sage: P = E([-1,1])
            sage: P.is_on_identity_component ()
            False
            sage: L.elliptic_logarithm(P, prec=96)
            0.4793482501902193161295330101 + 0.9858688507758241022112038491*I
            sage: Q=E([3,5])
            sage: Q.is_on_identity_component()
            True
            sage: L.elliptic_logarithm(Q, prec=96)
            1.931128271542559442488585220

        Note that this is actually the inverse of the Weierstrass isomorphism::

            sage: L.elliptic_exponential(_)  # abs tol 1e-26
            (3.000000000000000000000000000 : 5.000000000000000000000000000 : 1.000000000000000000000000000)

        An example with negative discriminant, and a torsion point::

            sage: E = EllipticCurve('11a1')
            sage: L = E.period_lattice()
            sage: E.discriminant() < 0
            True
            sage: L.real_flag
            -1
            sage: P = E([16,-61])
            sage: L.elliptic_logarithm(P)
            0.253841860855911
            sage: L.real_period() / L.elliptic_logarithm(P)
            5.00000000000000

        An example where precision is problematic::

            sage: E = EllipticCurve([1, 0, 1, -85357462, 303528987048])  #18074g1
            sage: P = E([4458713781401/835903744, -64466909836503771/24167649046528, 1])
            sage: L = E.period_lattice()
            sage: L.ei()
            [5334.003952567705? - 1.964393150436?e-6*I,
             5334.003952567705? + 1.964393150436?e-6*I,
             -10668.25790513541?]
            sage: L.elliptic_logarithm(P,prec=100)
            0.27656204014107061464076203097

        Some complex examples, taken from the paper by Cremona and Thongjunthug::

            sage: # needs sage.rings.number_field
            sage: K.<i> = QuadraticField(-1)
            sage: a4 = 9*i - 10
            sage: a6 = 21 - i
            sage: E = EllipticCurve([0,0,0,a4,a6])
            sage: e1 = 3 - 2*i; e2 = 1 + i; e3 = -4 + i
            sage: emb = K.embeddings(CC)[1]
            sage: L = E.period_lattice(emb)
            sage: P = E(2 - i, 4 + 2*i)

        By default, the output is reduced with respect to the
        normalised lattice basis, so that its coordinates with respect
        to that basis lie in the interval [0,1)::

            sage: z = L.elliptic_logarithm(P, prec=100); z                              # needs sage.rings.number_field
            0.70448375537782208460499649302 - 0.79246725643650979858266018068*I
            sage: L.coordinates(z)                                                      # needs sage.rings.number_field
            (0.46247636364807931766105406092, 0.79497588726808704200760395829)

        Using ``reduce=False`` this step can be omitted.  In this case
        the coordinates are usually in the interval [-0.5,0.5), but
        this is not guaranteed.  This option is mainly for testing
        purposes::

            sage: z = L.elliptic_logarithm(P, prec=100, reduce=False); z                # needs sage.rings.number_field
            0.57002153834710752778063503023 + 0.46476340520469798857457031393*I
            sage: L.coordinates(z)                                                      # needs sage.rings.number_field
            (0.46247636364807931766105406092, -0.20502411273191295799239604171)

        The elliptic logs of the 2-torsion points are half-periods::

            sage: L.elliptic_logarithm(E(e1, 0), prec=100)                              # needs sage.rings.number_field
            0.64607575874356525952487867052 + 0.22379609053909448304176885364*I
            sage: L.elliptic_logarithm(E(e2, 0), prec=100)                              # needs sage.rings.number_field
            0.71330686725892253793705940192 - 0.40481924028150941053684639367*I
            sage: L.elliptic_logarithm(E(e3, 0), prec=100)                              # needs sage.rings.number_field
            0.067231108515357278412180731396 - 0.62861533082060389357861524731*I

        We check this by doubling and seeing that the resulting
        coordinates are integers::

            sage: L.coordinates(2*L.elliptic_logarithm(E(e1, 0), prec=100))             # needs sage.rings.number_field
            (1.0000000000000000000000000000, 0.00000000000000000000000000000)
            sage: L.coordinates(2*L.elliptic_logarithm(E(e2, 0), prec=100))             # needs sage.rings.number_field
            (1.0000000000000000000000000000, 1.0000000000000000000000000000)
            sage: L.coordinates(2*L.elliptic_logarithm(E(e3, 0), prec=100))             # needs sage.rings.number_field
            (0.00000000000000000000000000000, 1.0000000000000000000000000000)

        ::

            sage: # needs sage.rings.number_field
            sage: a4 = -78*i + 104
            sage: a6 = -216*i - 312
            sage: E = EllipticCurve([0,0,0,a4,a6])
            sage: emb = K.embeddings(CC)[1]
            sage: L = E.period_lattice(emb)
            sage: P = E(3 + 2*i, 14 - 7*i)
            sage: L.elliptic_logarithm(P)
            0.297147783912228 - 0.546125549639461*I
            sage: L.coordinates(L.elliptic_logarithm(P))
            (0.628653378040238, 0.371417754610223)
            sage: e1 = 1 + 3*i; e2 = -4 - 12*i; e3 = -e1 - e2
            sage: L.coordinates(L.elliptic_logarithm(E(e1, 0)))
            (0.500000000000000, 0.500000000000000)
            sage: L.coordinates(L.elliptic_logarithm(E(e2, 0)))
            (1.00000000000000, 0.500000000000000)
            sage: L.coordinates(L.elliptic_logarithm(E(e3, 0)))
            (0.500000000000000, 0.000000000000000)

        TESTS:

        See :issue:`10026` and :issue:`11767`::

            sage: # needs sage.rings.number_field
            sage: K.<w> = QuadraticField(2)
            sage: E = EllipticCurve([0, -1, 1, -3*w - 4, 3*w + 4])
            sage: T = E.simon_two_descent(lim1=20, lim3=5, limtriv=20)
            sage: P, Q = T[2]
            sage: embs = K.embeddings(CC)
            sage: Lambda = E.period_lattice(embs[0])
            sage: Lambda.elliptic_logarithm(P, 100)
            4.7100131126199672766973600998
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^2 + x + 5)
            sage: E = EllipticCurve(K, [0,0,1,-3,-5])
            sage: P = E([0,a])
            sage: Lambda = P.curve().period_lattice(K.embeddings(ComplexField(600))[0])
            sage: Lambda.elliptic_logarithm(P, prec=600)
            -0.842248166487739393375018008381693990800588864069506187033873183845246233548058477561706400464057832396643843146464236956684557207157300006542470428493573195030603817094900751609464 - 0.571366031453267388121279381354098224265947866751130917440598461117775339240176310729173301979590106474259885638797913383502735083088736326391919063211421189027226502851390118943491*I
            sage: K.<a> = QuadraticField(-5)
            sage: E = EllipticCurve([1,1,a,a,0])
            sage: P = E(0, 0)
            sage: L = P.curve().period_lattice(K.embeddings(ComplexField())[0])
            sage: L.elliptic_logarithm(P, prec=500)
            1.17058357737548897849026170185581196033579563441850967539191867385734983296504066660506637438866628981886518901958717288150400849746892393771983141354 - 1.13513899565966043682474529757126359416758251309237866586896869548539516543734207347695898664875799307727928332953834601460994992792519799260968053875*I
            sage: L.elliptic_logarithm(P, prec=1000)
            1.17058357737548897849026170185581196033579563441850967539191867385734983296504066660506637438866628981886518901958717288150400849746892393771983141354014895386251320571643977497740116710952913769943240797618468987304985625823413440999754037939123032233879499904283600304184828809773650066658885672885 - 1.13513899565966043682474529757126359416758251309237866586896869548539516543734207347695898664875799307727928332953834601460994992792519799260968053875387282656993476491590607092182964878750169490985439873220720963653658829712494879003124071110818175013453207439440032582917366703476398880865439217473*I

        Elliptic curve over ``QQbar``::

            sage: E = EllipticCurve(QQbar, [sqrt(2), I])
            sage: L = E.period_lattice()
            sage: P = E.lift_x(3)
            sage: L.elliptic_logarithm(P)
            -1.97657221097437 - 1.05021415535949*I
            sage: L.elliptic_exponential(_)  # abs tol 1e-15
            (3.00000000000000 + 9.20856947066460e-16*I : -5.59022723358798 - 0.0894418024719718*I : 1.00000000000000)
            sage: L.elliptic_logarithm(P, prec=100)  # abs tol 1e-15
            -3.4730631218714889933426781799 + 0.44627675553762761312098773197*I
            sage: L.elliptic_exponential(_)  # abs tol 1e-28
            (3.0000000000000000000000000000 - 1.4773628579202938936348512161e-30*I : -5.5902272335879800026836302686 - 0.089441802471969391005702381090*I : 1.0000000000000000000000000000)

        Real approximate field, negative discriminant. Note that the output precision uses the precision of the base field::

            sage: E = EllipticCurve(RealField(100), [1, 6])
            sage: L = E.period_lattice()
            sage: L.real_flag
            -1
            sage: P = E(3, 6)
            sage: L.elliptic_logarithm(P)
            2.4593388737550379526023682666
            sage: L.elliptic_exponential(_)
            (3.0000000000000000000000000000 : 5.9999999999999999999999999999 : 1.0000000000000000000000000000)

        Real approximate field, positive discriminant::

            sage: E = EllipticCurve(RealField(100), [-4, 3])
            sage: L = E.period_lattice()
            sage: L.real_flag
            1
            sage: P = E.lift_x(4)
            sage: L.elliptic_logarithm(P)
            0.51188849089267627141925354967
            sage: L.elliptic_exponential(_)
            (4.0000000000000000000000000000 : -7.1414284285428499979993998114 : 1.0000000000000000000000000000)

        Complex approximate field::

            sage: E = EllipticCurve(ComplexField(100), [I, 3*I+4])
            sage: L = E.period_lattice()
            sage: L.real_flag
            0
            sage: P = E.lift_x(4)
            sage: L.elliptic_logarithm(P)
            -1.1447032790074574712147458157 - 0.72429843602171875396186134806*I
            sage: L.elliptic_exponential(_)
            (4.0000000000000000000000000000 + 1.2025589033682610849950210280e-30*I : -8.2570982991257407680322611854 - 0.42387771989714340809597881586*I : 1.0000000000000000000000000000)
        """
        if P.curve() is not self.E:
            raise ValueError("Point is on the wrong curve")
        if prec is None:
            prec = self._compute_default_prec()
        if P.is_zero():
            return ComplexField(prec)(0)

        # Compute the real or complex coordinates of P:

        xP, yP = (self.embedding(coord) for coord in P.xy())

        # The real work is done over R or C now:

        return self.e_log_RC(xP, yP, prec, reduce=reduce)

    def elliptic_exponential(self, z, to_curve=True):
        r"""
        Return the elliptic exponential of a complex number.

        INPUT:

        - ``z`` -- complex number (viewed modulo this period lattice)

        - ``to_curve`` -- boolean (default: ``True``); see below

        OUTPUT:

        - If ``to_curve`` is False, a 2-tuple of real or complex
          numbers representing the point `(x,y) = (\wp(z),\wp'(z))`
          where `\wp` denotes the Weierstrass `\wp`-function with
          respect to this lattice.

        - If ``to_curve`` is True, the point `(X,Y) =
          (x-b_2/12,y-(a_1(x-b_2/12)-a_3)/2)` as a point in `E(\RR)`
          or `E(\CC)`, with `(x,y) = (\wp(z),\wp'(z))` as above, where
          `E` is the elliptic curve over `\RR` or `\CC` whose period
          lattice this is.

        - If the lattice is real and `z` is also real then the output
          is a pair of real numbers if ``to_curve`` is ``True``, or a
          point in `E(\RR)` if ``to_curve`` is ``False``.

        .. NOTE::

            The precision is taken from that of the input ``z``.

        EXAMPLES::

            sage: E = EllipticCurve([1,1,1,-8,6])
            sage: P = E(1, -2)
            sage: L = E.period_lattice()
            sage: z = L(P); z
            1.17044757240090
            sage: L.elliptic_exponential(z)
            (0.999999999999999 : -2.00000000000000 : 1.00000000000000)
            sage: _.curve()
            Elliptic Curve defined by y^2 + 1.00000000000000*x*y + 1.00000000000000*y
             = x^3 + 1.00000000000000*x^2 - 8.00000000000000*x + 6.00000000000000
             over Real Field with 53 bits of precision
            sage: L.elliptic_exponential(z,to_curve=False)
            (1.41666666666667, -2.00000000000000)
            sage: z = L(P, prec=201); z
            1.17044757240089592298992188482371493504472561677451007994189
            sage: L.elliptic_exponential(z)
            (1.00000000000000000000000000000000000000000000000000000000000
             : -2.00000000000000000000000000000000000000000000000000000000000
             : 1.00000000000000000000000000000000000000000000000000000000000)

        Examples over number fields::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^3 - 2)
            sage: embs = K.embeddings(CC)
            sage: E = EllipticCurve('37a')
            sage: EK = E.change_ring(K)
            sage: Li = [EK.period_lattice(e) for e in embs]
            sage: P = EK(-1, -1)
            sage: Q = EK(a - 1, 1 - a^2)
            sage: zi = [L.elliptic_logarithm(P) for L in Li]
            sage: [c.real() for c in Li[0].elliptic_exponential(zi[0])]
            [-1.00000000000000, -1.00000000000000, 1.00000000000000]
            sage: [c.real() for c in Li[0].elliptic_exponential(zi[1])]
            [-1.00000000000000, -1.00000000000000, 1.00000000000000]
            sage: [c.real() for c in Li[0].elliptic_exponential(zi[2])]
            [-1.00000000000000, -1.00000000000000, 1.00000000000000]

            sage: # needs sage.rings.number_field
            sage: zi = [L.elliptic_logarithm(Q) for L in Li]
            sage: Li[0].elliptic_exponential(zi[0])
            (-1.62996052494744 - 1.09112363597172*I
             : 1.79370052598410 - 1.37472963699860*I
             : 1.00000000000000)
            sage: [embs[0](c) for c in Q]
            [-1.62996052494744 - 1.09112363597172*I,
             1.79370052598410 - 1.37472963699860*I,
             1.00000000000000]
            sage: Li[1].elliptic_exponential(zi[1])
            (-1.62996052494744 + 1.09112363597172*I
             : 1.79370052598410 + 1.37472963699860*I
             : 1.00000000000000)
            sage: [embs[1](c) for c in Q]
            [-1.62996052494744 + 1.09112363597172*I,
             1.79370052598410 + 1.37472963699860*I,
             1.00000000000000]
            sage: [c.real() for c in Li[2].elliptic_exponential(zi[2])]
            [0.259921049894873, -0.587401051968199, 1.00000000000000]
            sage: [embs[2](c) for c in Q]
            [0.259921049894873, -0.587401051968200, 1.00000000000000]

        Test to show that :issue:`8820` is fixed::

            sage: # needs sage.rings.number_field
            sage: E = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-5)
            sage: L = E.change_ring(K).period_lattice(K.places()[0])
            sage: L.elliptic_exponential(CDF(.1,.1))  # abs tol 1e-15
            (0.0000142854026029 - 49.9960001066650*I
             : 249.520141250950 + 250.019855549131*I : 1.00000000000000)
            sage: L.elliptic_exponential(CDF(.1,.1), to_curve=False)
            (0.0000142854026029447 - 49.9960001066650*I,
             500.040282501900 + 500.039711098263*I)

        `z=0` is treated as a special case::

            sage: E = EllipticCurve([1,1,1,-8,6])
            sage: L = E.period_lattice()
            sage: L.elliptic_exponential(0)
            (0.000000000000000 : 1.00000000000000 : 0.000000000000000)
            sage: L.elliptic_exponential(0, to_curve=False)
            (+infinity, +infinity)

        ::

            sage: # needs sage.rings.number_field
            sage: E = EllipticCurve('37a')
            sage: K.<a> = QuadraticField(-5)
            sage: L = E.change_ring(K).period_lattice(K.places()[0])
            sage: P = L.elliptic_exponential(0); P
            (0.000000000000000 : 1.00000000000000 : 0.000000000000000)
            sage: P.parent()
            Abelian group of points on Elliptic Curve defined by
             y^2 + 1.00000000000000*y = x^3 + (-1.00000000000000)*x
             over Complex Field with 53 bits of precision

        Very small `z` are handled properly (see :issue:`8820`)::

            sage: # needs sage.rings.number_field
            sage: K.<a> = QuadraticField(-1)
            sage: E = EllipticCurve([0,0,0,a,0])
            sage: L = E.period_lattice(K.complex_embeddings()[0])
            sage: L.elliptic_exponential(1e-100)
            (0.000000000000000 : 1.00000000000000 : 0.000000000000000)

        The elliptic exponential of `z` is returned as `(0 : 1 : 0)` if
        the coordinates of `z` with respect to the period lattice are
        approximately integral::

            sage: (100/log(2.0,10))/0.8
            415.241011860920
            sage: L.elliptic_exponential((RealField(415)(1e-100))).is_zero()            # needs sage.rings.number_field
            True
            sage: L.elliptic_exponential((RealField(420)(1e-100))).is_zero()            # needs sage.rings.number_field
            False
        """
        C = z.parent()
        z_is_real = False
        if isinstance(C, sage.rings.abc.RealField):
            z_is_real = True
            C = ComplexField(C.precision())
            z = C(z)
        else:
            if isinstance(C, sage.rings.abc.ComplexField):
                z_is_real = z.is_real()
            else:
                try:
                    C = ComplexField()
                    z = C(z)
                    z_is_real = z.is_real()
                except TypeError:
                    raise TypeError("%s is not a complex number" % z)
        prec = C.precision()

        # test for the point at infinity:

        eps = (C(2)**(-0.8*prec)).real()  ## to test integrality w.r.t. lattice within 20%
        if all((t.round()-t).abs() < eps for t in self.coordinates(z)):
            K = z.parent()
            if to_curve:
                return self.curve().change_ring(K)(0)
            else:
                return K(Infinity), K(Infinity)

        # general number field code (including QQ):

        # We do not use PARI's ellztopoint function since it is only
        # defined for curves over the reals (note that PARI only
        # computes the period lattice basis in that case).  But Sage
        # can compute the period lattice basis over CC, and then
        # PARI's ellwp function works fine.

        # NB converting the PARI values to Sage values might land up
        # in real/complex fields of spuriously higher precision than
        # the input, since PARI's precision is in word-size chunks.
        # So we force the results back into the real/complex fields of
        # the same precision as the input.

        x, y = pari(self.basis(prec=prec)).ellwp(z, flag=1)
        x, y = (C(t) for t in (x, y))

        if self.real_flag and z_is_real:
            x = x.real()
            y = y.real()

        if to_curve:
            K = x.parent()
            if self.is_approximate:
                v = self.embedding
            else:
                v = refine_embedding(self.embedding, Infinity)
            a1, a2, a3, a4, a6 = (K(v(a)) for a in self.E.ainvs())
            b2 = K(v(self.E.b2()))
            x = x - b2 / 12
            y = (y - (a1 * x + a3)) / 2
            EK = EllipticCurve(K, [a1, a2, a3, a4, a6])
            return EK.point((x, y, K.one()), check=False)
        else:
            return (x, y)


def reduce_tau(tau):
    r"""
    Transform a point in the upper half plane to the fundamental region.

    INPUT:

    - ``tau`` -- complex number with positive imaginary part

    OUTPUT:

    (tuple) `(\tau',[a,b,c,d])` where `a,b,c,d` are integers such that

    - `ad-bc=1`;
    - `\tau'=(a\tau+b)/(c\tau+d)`;
    - `|\tau'|\ge1`;
    - `|\Re(\tau')|\le\frac{1}{2}`.

    EXAMPLES::

        sage: # needs sage.rings.real_mpfr sage.symbolic
        sage: from sage.schemes.elliptic_curves.period_lattice import reduce_tau
        sage: reduce_tau(CC(1.23,3.45))
        (0.230000000000000 + 3.45000000000000*I, [1, -1, 0, 1])
        sage: reduce_tau(CC(1.23,0.0345))
        (-0.463960069171512 + 1.35591888067914*I, [-5, 6, 4, -5])
        sage: reduce_tau(CC(1.23,0.0000345))
        (0.130000000001761 + 2.89855072463768*I, [13, -16, 100, -123])
    """
    assert tau.imag() > 0
    a, b = ZZ(1), ZZ(0)
    c, d = b, a
    k = tau.real().round()
    tau -= k
    a -= k*c
    b -= k*d
    while tau.abs() < 0.999:
        tau = -1/tau
        a, b, c, d = c, d, -a, -b
        k = tau.real().round()
        tau -= k
        a -= k*c
        b -= k*d
    assert a*d-b*c == 1
    assert tau.abs() >= 0.999 and tau.real().abs() <= 0.5
    return tau, [a,b,c,d]


def normalise_periods(w1, w2):
    r"""
    Normalise the period basis `(w_1,w_2)` so that `w_1/w_2` is in the fundamental region.

    INPUT:

    - ``w1``, ``w2`` -- two complex numbers with non-real ratio

    OUTPUT:

    (tuple) `((\omega_1',\omega_2'),[a,b,c,d])` where `a,b,c,d` are
    integers such that

    - `ad-bc=\pm1`;
    - `(\omega_1',\omega_2') = (a\omega_1+b\omega_2,c\omega_1+d\omega_2)`;
    - `\tau=\omega_1'/\omega_2'` is in the upper half plane;
    - `|\tau|\ge1` and `|\Re(\tau)|\le\frac{1}{2}`.

    EXAMPLES::

        sage: # needs sage.rings.real_mpfr sage.symbolic
        sage: from sage.schemes.elliptic_curves.period_lattice import reduce_tau, normalise_periods
        sage: w1 = CC(1.234, 3.456)
        sage: w2 = CC(1.234, 3.456000001)
        sage: w1/w2    # in lower half plane!
        0.999999999743367 - 9.16334785827644e-11*I
        sage: w1w2, abcd = normalise_periods(w1, w2)
        sage: a,b,c,d = abcd
        sage: w1w2 == (a*w1+b*w2, c*w1+d*w2)
        True
        sage: w1w2[0]/w1w2[1]
        1.23400010389203e9*I
        sage: a*d-b*c # note change of orientation
        -1
    """
    tau = w1/w2
    s = +1
    if tau.imag() < 0:
        w2 = -w2
        tau = -tau
        s = -1
    tau, abcd = reduce_tau(tau)
    a, b, c, d = abcd
    if s < 0:
        abcd = (a,-b,c,-d)
    return (a*w1+b*w2,c*w1+d*w2), abcd


def extended_agm_iteration(a, b, c):
    r"""
    Internal function for the extended AGM used in elliptic logarithm computation.

    INPUT:

    - ``a``, ``b``, ``c`` -- three real or complex numbers

    OUTPUT:

    (3-tuple) `(a_0,b_0,c_0)`, the limit of the iteration `(a,b,c) \mapsto ((a+b)/2,\sqrt{ab},(c+\sqrt{c^2+b^2-a^2})/2)`.

    EXAMPLES::

        sage: # needs sage.rings.real_mpfr
        sage: from sage.schemes.elliptic_curves.period_lattice import extended_agm_iteration
        sage: extended_agm_iteration(RR(1), RR(2), RR(3))
        (1.45679103104691, 1.45679103104691, 3.21245294970054)
        sage: extended_agm_iteration(CC(1,2), CC(2,3), CC(3,4))
        (1.46242448156430 + 2.47791311676267*I,
         1.46242448156430 + 2.47791311676267*I,
         3.22202144343535 + 4.28383734262540*I)

    TESTS::

        sage: extended_agm_iteration(1,2,3)                                             # needs sage.rings.real_mpfr
        Traceback (most recent call last):
        ...
        ValueError: values must be real or complex numbers
    """
    if not isinstance(a, (RealNumber,ComplexNumber)):
        raise ValueError("values must be real or complex numbers")
    eps = a.parent().one().real() >> (a.parent().precision() - 10)
    while True:
        a1 = (a + b) / 2
        b1 = (a * b).sqrt()
        delta = (b**2 - a**2) / c**2
        f = (1 + (1 + delta).sqrt()) / 2
        if (f.abs() - 1).abs() < eps:
            return a, b, c
        c *= f
        a, b = a1, b1
