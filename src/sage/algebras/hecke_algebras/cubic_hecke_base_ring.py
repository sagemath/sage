# sage.doctest: needs sage.libs.pari (for factorization)
r"""
Cubic Hecke Base Rings

This module contains special classes of polynomial rings
(:class:`CubicHeckeRingOfDefinition` and :class:`CubicHeckeExtensionRing`)
used in the context of
:class:`cubic Hecke algebras
<sage.algebras.hecke_algebras.cubic_hecke_algebra.CubicHeckeAlgebra>`.

AUTHORS:

- Sebastian Oehms May 2020: initial version
"""
# ###########################################################################
#       Copyright (C) 2020 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ###########################################################################
from sage.structure.category_object import normalize_names
from sage.structure.element import get_coercion_model
from sage.categories.action import Action
from sage.misc.verbose import verbose
from sage.misc.functional import cyclotomic_polynomial
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing_mpair
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.localization import Localization
from sage.algebras.splitting_algebra import solve_with_extension, SplittingAlgebra


# -----------------------------------------------------------------------------
# local helper functions
# -----------------------------------------------------------------------------
def normalize_names_markov(names, markov_trace_version):
    r"""
    Return a tuple of strings of variable names of length 3 resp. 4 (if
    ``markov_trace_version`` is ``True``) according to the given input names.

    INPUT:

    - ``names`` -- passed to :func:`~sage.structure.category_object.normalize_names`
    - ``markov_trace_version`` -- boolean; if set to ``True`` four names are
      expected the last of which corresponds to the writhe factor of the
      Markov trace

    EXAMPLES::

        sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
        sage: chbr.normalize_names_markov('a, b, c', False)
        ('a', 'b', 'c')
        sage: chbr.normalize_names_markov(('u', 'v', 'w', 's'), False)
        ('u', 'v', 'w')
    """
    if markov_trace_version:
        names = normalize_names(4, names)
    else:
        if isinstance(names, tuple):
            names = list(names)
        if isinstance(names, list) and len(names) > 3:
            names = normalize_names(3, names[0:3])
        else:
            names = normalize_names(3, names)
    return names


def register_ring_hom(ring_hom):
    r"""
    Register the given ring homomorphism as conversion map.

    EXAMPLES::

        sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
        sage: BR = chbr.CubicHeckeRingOfDefinition()
        sage: BR.create_specialization([E(5), E(7), E(3)])  # indirect doctest
        Universal Cyclotomic Field
        sage: _.convert_map_from(BR)
        Ring morphism:
          From: Multivariate Polynomial Ring in u, v, w
                  over Integer Ring localized at (w,)
          To:   Universal Cyclotomic Field
          Defn: u |--> E(5)
                v |--> E(7)
                w |--> E(3)
    """
    domain = ring_hom.domain()
    codomain = ring_hom.codomain()
    conversion_cached = codomain._is_conversion_cached(domain)

    if conversion_cached:
        test_map = codomain.convert_map_from(domain)
        try:
            if test_map != ring_hom:
                verbose('\nConversion:\n%s\n already exists and is different from:\n%s\n' % (test_map, ring_hom))
        except TypeError:
            verbose('\n Conversion:\n%s\n already exists and is not comparable to:\n%s\n' % (test_map, ring_hom))
    else:
        try:
            codomain.register_conversion(ring_hom)
        except ValueError:
            verbose('\nthe map:\n%s\ncannot be registered as conversion\n' % ring_hom)

    return


# -----------------------------------------------------------------------------
# class for the Galois Group action on the generic extension ring corresponding
# to the cubic equation
# -----------------------------------------------------------------------------
class GaloisGroupAction(Action):
    r"""
    Action on a multivariate polynomial ring by permuting the generators.

    EXAMPLES::

        sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
        sage: from operator import mul
        sage: R.<x, y, z> = ZZ[]
        sage: G = SymmetricGroup(3)
        sage: p = 5*x*y + 3*z**2
        sage: R._unset_coercions_used()
        sage: R.register_action(chbr.GaloisGroupAction(G, R, op=mul))
        sage: s = G([2,3,1])
        sage: s*p
        3*x^2 + 5*y*z
    """
    def _act_(self, perm, pol):
        r"""
        Application of the action.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: from operator import mul
            sage: R.<x, y> = QQ[]
            sage: G = SymmetricGroup(2)
            sage: A = chbr.GaloisGroupAction(G, R, op=mul)
            sage: p = ~5*x*y**2 + 3*x**2
            sage: s = G([2,1])
            sage: A._act_(s, p)
            1/5*x^2*y + 3*y^2
        """
        if not self.is_left():
            perm, pol = pol, perm
        pol_dict = {}
        for key, value in pol.monomial_coefficients().items():
            newkey = [0] * len(key)
            for pos, k in enumerate(key):
                newkey[perm(pos + 1) - 1] = k
            pol_dict[tuple(newkey)] = value
        return self.domain()(pol_dict)


###############################################################################
# EXTENSION RING
# -----------------------------------------------------------------------------
# Definition of the generic extension ring for the cubic Hecke algebra as
# Laurent polynomial ring in 3 indeterminates over the cyclotomic field of a
# third root of unity This is the most general ring over which the cubic Hecke
# algebra is semi-simple. In opposite to the generic base ring class, this class
# does not inherits from UniqueRepresentation since _test_pickling fails
# -----------------------------------------------------------------------------
class CubicHeckeExtensionRing(LaurentPolynomialRing_mpair):
    r"""
    The generic splitting algebra for the irreducible representations of
    the cubic Hecke algebra.

    This ring must contain three invertible indeterminates (representing
    the roots of the cubic equation) together with a third root of unity
    (needed for the 18-dimensional irreducibles of the cubic Hecke algebra
    on 4 strands).

    Therefore this ring is constructed as a multivariate Laurent polynomial
    ring in three indeterminates over a polynomial quotient ring over the
    integers with respect to the minimal polynomial of a third root of unity.

    The polynomial quotient ring is constructed as instance of
    :class:`SplittingAlgebra`.

    INPUT:

    - ``names`` -- (default: ``'u,v,w'``) string containing the names of the
      indeterminates separated by ``,`` or a triple of strings each of which
      are the names of one of the three indeterminates
    - ``order`` -- string (default: ``'degrevlex'``); the term order; see also
      :class:`~sage.rings.polynomial.laurent_polynomial_ring.LaurentPolynomialRing_mpair`
    - ``ring_of_definition`` -- (optional) a :class:`CubicHeckeRingOfDefinition`
      to specify the generic cubic Hecke base ring over which ``self`` may be
      realized as splitting ring via the ``as_splitting_algebra`` method
    - ``third_unity_root_name`` -- string (default: ``'e3'``); for setting the
      name of the third root if unity of ``self``
    - ``markov_trace_version`` -- boolean (default: ``False``); if this is
      set to ``True`` then ``self`` contains one invertible indeterminate in
      addition which is meant to represent the writhe factor of a Markov trace
      on the cubic Hecke algebra and which default name is ``s``

    EXAMPLES::

        sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
        sage: chbr.CubicHeckeExtensionRing('a, b, c')
        Multivariate Laurent Polynomial Ring in a, b, c
          over Splitting Algebra of x^2 + x + 1
            with roots [e3, -e3 - 1]
          over Integer Ring
        sage: _.an_element()
        b^2*c^-1 + e3*a
    """
    def __init__(self, names, order='degrevlex', ring_of_definition=None, third_unity_root_name='e3', markov_trace_version=False):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: ER = chbr.CubicHeckeExtensionRing('a, b, c')
            sage: TestSuite(ER).run()
        """
        # ----------------------------------------------------------------------
        # Setting connection with generic base ring (if given)
        # ----------------------------------------------------------------------
        self._ring_of_definition = None
        self._splitting_algebra = None

        if ring_of_definition is not None:
            if not isinstance(ring_of_definition, CubicHeckeRingOfDefinition):
                raise TypeError("generic base ring must be an instance of CubicHeckeRingOfDefinition")
            self._ring_of_definition = ring_of_definition

        # ----------------------------------------------------------------------
        # defining the base ring
        # note that we can't use ZZ.extension since it isn't possible to define
        # homomorphisms from orders in number fields, yet
        # ----------------------------------------------------------------------
        base_ring = SplittingAlgebra(cyclotomic_polynomial(3), [third_unity_root_name])

        # ----------------------------------------------------------------------
        # defining the ring itself
        # ----------------------------------------------------------------------
        self._names = normalize_names_markov(names, markov_trace_version)
        self._order = order

        pol_ring = PolynomialRing(base_ring, names=self._names, order=self._order, implementation=None)
        LaurentPolynomialRing_mpair.__init__(self, pol_ring)

        # ----------------------------------------------------------------------
        # setting Galois group action
        # ----------------------------------------------------------------------
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        from operator import mul
        self._galois_group = SymmetricGroup(3)
        galois_group_action = GaloisGroupAction(self._galois_group, self, op=mul)
        self._unset_coercions_used()
        self.register_action(galois_group_action)

        # ----------------------------------------------------------------------
        # Init of data used on demand
        # ----------------------------------------------------------------------
        self._mirror = None

    ############################################################################
    # overloaded inherited methods
    ############################################################################
    def construction(self):
        r"""
        Return ``None`` since this construction is not functorial.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: ER = chbr.CubicHeckeExtensionRing('a, b, c')
            sage: ER._test_category()   # indirect doctest
        """
        return None

    def __reduce__(self):
        r"""
        Used in pickling.

        TESTS::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: ER = chbr.CubicHeckeExtensionRing('a, b, c')
            sage: loads(dumps(ER)) == ER
            True
        """
        return CubicHeckeExtensionRing, (self._names, self._order, self._ring_of_definition)

    def _element_constructor_(self, x, mon=None):
        r"""
        Inherited element constructor overloaded to allow construction from
        ``GAP3`` ``MVP`` expressions.

        EXAMPLES::

            sage: # optional - gap3
            sage: CHA3 = algebras.CubicHecke(3)
            sage: GER = CHA3.extension_ring(generic=True)
            sage: sch7 = CHA3.chevie().SchurElements()[7]
            sage: GER(sch7)
            a*b*c^-2 + a^2*b^-1*c^-1 + a^-1*b^2*c^-1 + 2
            + a*b^-2*c + a^-2*b*c + a^-1*b^-1*c^2
            sage: rep4_gap3 = CHA3.chevie().Representations(4)
            sage: matrix(GER, rep4_gap3[1])
            [ b  0]
            [-b  c]
        """
        from sage.interfaces.gap3 import GAP3Element
        if isinstance(x, GAP3Element):
            return self._convert_from_gap3_mvp(x)
        return super()._element_constructor_(x, mon=mon)

    def _coerce_map_from_(self, R):
        r"""
        The rings that canonically coerce to ``self`` ar the ones from
        inheritance and the base ring of definition of the cubic Hecke algebra.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: BR = chbr.CubicHeckeRingOfDefinition()
            sage: ER = BR.extension_ring()
            sage: ER(BR.an_element())
            a*b + a*c + b*c + a*b^-1*c^-1 + 2*c^-1 + a^-1*b*c^-1 + 2*b^-1
            + 2*a^-1 + a^-1*b^-1*c
            sage: MBR = chbr.CubicHeckeRingOfDefinition(markov_trace_version=True)
            sage: MER = MBR.extension_ring()
            sage: MER(MBR.an_element())
            a*b*s^-1 + a*c*s^-1 + b*c*s^-1 + a*b^-1*c^-1 + 2*c^-1
            + a^-1*b*c^-1 + 2*b^-1 + 2*a^-1 + a^-1*b^-1*c
        """
        if isinstance(R, CubicHeckeRingOfDefinition):
            markov = R.markov_trace_version()
            a, b, c, *rem = self.gens()
            iu = a + b + c
            iv = a * b + a * c + b * c
            iw = a * b * c
            im_gens = [iu, iv, iw]
            if markov:
                if self.markov_trace_version():
                    im_gens += rem
                    # check of embedding fails in this case as long as the images of
                    # ``iu`` and ``iv`` need to be invertible (see comment in
                    # :meth:`__init__`).  # :class:`CubicHeckeRingOfDefinition`).
                    embedding_into_extension_ring = R.hom(im_gens, check=False)
            else:
                embedding_into_extension_ring = R.hom(im_gens)
            return embedding_into_extension_ring
        return super()._coerce_map_from_(R)

    def hom(self, im_gens, codomain=None, check=True, base_map=None):
        r"""
        Return a homomorphism of ``self``.

        INPUT:

        - ``im_gens`` -- tuple for the image of the generators of ``self``
        - ``codomain`` -- (optional) the codomain of the homomorphism

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: ER = chbr.CubicHeckeExtensionRing('a, b, c')
            sage: UCF = UniversalCyclotomicField()
            sage: map = ER.hom((UCF.gen(3),) + (UCF(3),UCF(4),UCF(5)))
            sage: ER.an_element()
            b^2*c^-1 + e3*a
            sage: map(_)
            -1/5*E(3) - 16/5*E(3)^2
        """
        gens = self.gens()
        num_gens = len(gens)

        if not isinstance(im_gens, (list, tuple)):
            im_gens = [im_gens]

        if len(im_gens) == num_gens + 1:
            e3, *im_remain = im_gens
            hom_cycl_gen = self.base_ring().hom([e3], codomain=e3.parent(), check=check, base_map=base_map)
            verbose("hom_cycl_gen %s" % hom_cycl_gen, level=2)
            return super().hom(im_remain, codomain=codomain, check=check, base_map=hom_cycl_gen)
        else:
            if base_map is None:
                raise ValueError('number of images must be four (including a '
                                 'third root of unity at first position) or a '
                                 'base_map (on %s) must be given' % self.base_ring())
            return super().hom(im_gens, codomain=codomain, check=check, base_map=base_map)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: ER = chbr.CubicHeckeExtensionRing('x, y, z')
            sage: ER.an_element()                             # indirect doctest
            y^2*z^-1 + e3*x
            sage: MER = chbr.CubicHeckeExtensionRing('x, y, z, s', markov_trace_version=True)
            sage: MER.an_element()                            # indirect doctest
            y^2*z^-1 + e3*x*s^-1
        """
        a, b, c, *rem = self.gens()
        e3 = self.cyclotomic_generator()
        s = self.one()
        if rem:
            s = rem[0]
        return b**2/c + a*e3/s

    ############################################################################
    # local methods
    ############################################################################
    def _is_markov_trace_version(self) -> bool:
        r"""
        Return whether ``self`` is the version containing the writhe parameter
        ``s`` for the Markov trace.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: ER = chbr.CubicHeckeExtensionRing('a, b, c')
            sage: ER._is_markov_trace_version()
            False
            sage: MER = chbr.CubicHeckeExtensionRing('a, b, c, s', markov_trace_version=True)
            sage: MER._is_markov_trace_version()
            True
        """
        return len(self.gens()) == 4

    # --------------------------------------------------------------------------
    # helper for element construction
    # --------------------------------------------------------------------------
    def _convert_from_gap3_mvp(self, mvp_expression):
        r"""
        Convert an element of ``GAP3`` interface containing Jean
        Michel's ``MVP`` (multivariate polynomials) to an element of ``self``.

        INPUT:

        - ``mvp_expression``  -- element of ``GAP3`` interface containing
          Jean Michel's ``MVP`` (multivariate polynomials)

        EXAMPLES::

            sage: # optional - gap3
            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: CHA3 = algebras.CubicHecke(3)
            sage: sch7 = CHA3.chevie().SchurElements()[7]
            sage: ER = chbr.CubicHeckeExtensionRing('a, b, c')
            sage: ER._convert_from_gap3_mvp(sch7)
            a*b*c^-2 + a^2*b^-1*c^-1 + a^-1*b^2*c^-1 + 2 + a*b^-2*c + a^-2*b*c + a^-1*b^-1*c^2
        """
        from sage.misc.sage_eval import sage_eval
        E3 = self.cyclotomic_generator()
        a, b, c, *rem = self.gens()
        na, nb, nc = self.variable_names()
        lc = {na: a, nb: b, nc: c, 'E3': E3}
        sage_expr = str(mvp_expression.FormatMaple())[1:-1]
        return sage_eval(sage_expr, locals=lc)

    ############################################################################
    # global methods
    ############################################################################
    def cyclotomic_generator(self):
        r"""
        Return the third root of unity as generator of the base ring
        of ``self``.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: ER  = chbr.CubicHeckeExtensionRing('a, b, c')
            sage: ER.cyclotomic_generator()
            e3
            sage: _**3 == 1
            True
        """
        return self(self.base_ring().gen())

    def conjugation(self):
        r"""
        Return an involution that performs *complex conjugation* with respect
        to base ring considered as order in the complex field.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: ER = chbr.CubicHeckeExtensionRing('x, y, z')
            sage: conj = ER.conjugation()
            sage: conj(ER.an_element())
            y^2*z^-1 + (-e3 - 1)*x
            sage: MER = chbr.CubicHeckeExtensionRing('x, y, z, s', markov_trace_version=True)
            sage: conj = MER.conjugation()
            sage: conj(MER.an_element())
            y^2*z^-1 + (-e3 - 1)*x*s^-1
        """
        e3 = self.cyclotomic_generator()
        return self.hom(tuple([e3**2] + list(self.gens())))

    def cubic_equation_galois_group(self):
        r"""
        Return the Galois group of the cubic equation, which is the permutation
        group on the three generators together with its action on ``self``.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: ER = chbr.CubicHeckeExtensionRing('a, b, c')
            sage: G = ER.cubic_equation_galois_group()
            sage: t = ER.an_element()
            sage: [(g ,g*t) for g in G]
            [((), b^2*c^-1 + e3*a),
            ((1,3,2), a^2*b^-1 + e3*c),
            ((1,2,3), e3*b + a^-1*c^2),
            ((2,3), e3*a + b^-1*c^2),
            ((1,3), a^-1*b^2 + e3*c),
            ((1,2), a^2*c^-1 + e3*b)]
        """
        return self._galois_group

    def mirror_involution(self):
        r"""
        Return the involution of ``self`` corresponding to the involution of
        the cubic Hecke algebra (with the same name).

        This means that it maps the generators of ``self`` to their inverses.

        .. NOTE::

           The mirror involution of the braid group does not factor through the
           cubic Hecke algebra over its base ring, but it does if it is
           considered as `\ZZ`-algebra. The base ring elements are transformed
           by this automorphism.

        OUTPUT: the involution as automorphism of ``self``

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: ER = chbr.CubicHeckeExtensionRing('p, q, r')
            sage: ER.mirror_involution()
            Ring endomorphism of Multivariate Laurent Polynomial Ring in p, q, r
                                 over Splitting Algebra of x^2 + x + 1
                                   with roots [e3, -e3 - 1]
                                 over Integer Ring
              Defn: p |--> p^-1
                    q |--> q^-1
                    r |--> r^-1
                    with map of base ring
            sage: _(ER.an_element())
            e3*p^-1 + q^-2*r

            sage: MER = chbr.CubicHeckeExtensionRing('p, q, r, s', markov_trace_version=True)
            sage: MER.mirror_involution()
            Ring endomorphism of Multivariate Laurent Polynomial Ring in p, q, r, s
              over Splitting Algebra of x^2 + x + 1
                with roots [e3, -e3 - 1] over Integer Ring
            Defn: p |--> p^-1
            q |--> q^-1
            r |--> r^-1
            s |--> s^-1
            with map of base ring
            sage: _(MER.an_element())
            e3*p^-1*s + q^-2*r
        """
        if self._mirror is None:
            e3 = self.base_ring().gen()
            if self._is_markov_trace_version():
                a, b, c, s = self.gens()
                self._mirror = self.hom([e3, ~a, ~b, ~c, ~s])
            else:
                a, b, c = self.gens()
                self._mirror = self.hom([e3, ~a, ~b, ~c])

        return self._mirror

    def create_specialization(self, im_cubic_equation_roots, im_writhe_parameter=None, var='T', third_unity_root_name='E3'):
        r"""
        Return an appropriate ring containing the elements from the list
        ``im_cubic_equation_roots`` defining a conversion map from ``self`` mapping
        the cubic equation roots of ``self`` to ``im_cubic_equation_roots``.

        INPUT:

        - ``im_cubic_equation_roots`` -- list or tuple of three ring elements
          such that there exists a ring homomorphism from the corresponding
          elements of ``self`` to them

        OUTPUT:

        A common parent containing the elements of ``im_cubic_equation_roots``
        together with their inverses.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: ER = chbr.CubicHeckeExtensionRing('a, b, c')
            sage: t = ER.an_element(); t
            b^2*c^-1 + e3*a
            sage: Sp1 = ER.create_specialization([E(5), E(7), E(3)]); Sp1
            Universal Cyclotomic Field
            sage: Sp1(t)
            -E(105)^11 - E(105)^16 - E(105)^26 - E(105)^37 - E(105)^41
            - E(105)^58 - E(105)^71 - E(105)^79 - E(105)^86 - E(105)^101
            sage: MER = chbr.CubicHeckeExtensionRing('a, b, c, s', markov_trace_version=True)
            sage: MER.create_specialization([E(5), E(7), E(3)], im_writhe_parameter=E(4))
            Universal Cyclotomic Field
            sage: a, b, c, s = MER.gens()
            sage: Sp1(MER(t)/s)
            E(420) + E(420)^29 + E(420)^89 + E(420)^149 + E(420)^169 + E(420)^209
            + E(420)^253 + E(420)^269 + E(420)^337 + E(420)^389

            sage: Z3 = CyclotomicField(3); E3=Z3.gen()
            sage: Sp2 = ER.create_specialization([E3, E3**2, Z3(1)])
            sage: Sp2(t)
            -1
            sage: MER.create_specialization([E3, E3**2, 1], im_writhe_parameter=2)
            Cyclotomic Field of order 3 and degree 2
            sage: Sp2(MER(t)*s)
            -2

            sage: Sp3 = ER.create_specialization([5, 7, 11])
            sage: Sp3(t)
            5*E3 + 49/11
        """
        # ----------------------------------------------------------------------
        # interpreting user given cubic equation roots and define the
        # corresponding specialized extension ring.
        # ----------------------------------------------------------------------

        if isinstance(im_cubic_equation_roots, tuple):
            im_cubic_equation_roots = list(im_cubic_equation_roots)

        if not isinstance(im_cubic_equation_roots, list):
            raise TypeError('cubic_equation_roots must be a list of three elements')

        if len(im_cubic_equation_roots) != 3:
            raise ValueError('there must be exactly three cubic_equation_roots')

        gens = self.gens()
        num_gens = len(gens)
        if im_writhe_parameter:
            if num_gens < 4:
                raise ValueError('im_writhe_parameter only possible for Markov-trace extension')
            im_gens = im_cubic_equation_roots + [im_writhe_parameter]
            a, b, c, s = im_gens
        else:
            if num_gens == 4:
                raise ValueError('im_writhe_parameter must be given for Markov-trace extension')
            im_gens = im_cubic_equation_roots
            a, b, c = im_gens

        image_ring = get_coercion_model().common_parent(*(im_gens))

        # ----------------------------------------------------------------------
        # make sure that all given cubic equation roots and their inverses
        # belong to image_ring
        # ----------------------------------------------------------------------
        try:
            image_ring = image_ring.localization(tuple(im_gens))
        except ValueError:
            pass

        im_gens = [image_ring(root) for root in im_gens]
        verbose('common parent of roots and inverses: %s' % (image_ring), level=2)

        image_ring_base = image_ring.base_ring()
        image_ring_map = None

        verbose('first choice: image_ring %s, image_ring_base %s' % (image_ring, image_ring_base), level=2)

        # ----------------------------------------------------------------------
        # make sure that a third root of unity belongs to image_ring
        # ----------------------------------------------------------------------

        E3 = None
        cp3 = cyclotomic_polynomial(3, var=var).change_ring(image_ring)
        cyclotomic_roots = solve_with_extension(cp3, [third_unity_root_name], var=var, flatten=True, warning=False)

        if len(cyclotomic_roots) > 0:
            E3 = cyclotomic_roots[0]
            verbose('third root of unity %s found in %s' % (E3, E3.parent()), level=2)

        if E3 is None:
            raise RuntimeError('cannot find a ring containing a third root of unity for the this choice of cubic roots!')

        hom_gens = [E3] + im_gens
        verbose('hom_gens %s' % hom_gens, level=2)

        image_ring = get_coercion_model().common_parent(*(hom_gens))
        verbose('common parent of roots and third root: %s' % image_ring, level=2)

        hom_gens = [image_ring(gen) for gen in hom_gens]

        image_ring_base = image_ring.base_ring()

        verbose('second choice: image_ring %s, image_ring_base %s' % (image_ring, image_ring_base), level=2)

        try:
            image_ring_map = self.hom(hom_gens, codomain=image_ring)
        except (ValueError, NotImplementedError):
            image_ring_map = self.hom(hom_gens, codomain=image_ring, check=False)
            verbose('check failed for embedding as ring morphism')

        verbose('specializing map defined %s' % image_ring_map, level=2)

        register_ring_hom(image_ring_map)
        return image_ring

    def as_splitting_algebra(self):
        r"""
        Return ``self`` as a :class:`SplittingAlgebra`; that is as an
        extension ring of the corresponding cubic Hecke algebra base ring
        (``self._ring_of_definition``, as a :class:`CubicHeckeRingOfDefinition`)
        splitting its cubic equation into linear factors, such that the roots
        are images of the generators of ``self``.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: GBR = chbr.CubicHeckeRingOfDefinition()
            sage: GER = GBR.extension_ring()
            sage: ER = GER.as_splitting_algebra(); ER
            Splitting Algebra of T^2 + T + 1 with roots [E3, -E3 - 1]
              over Splitting Algebra of h^3 - u*h^2 + v*h - w
                with roots [a, b, -b - a + u]
              over Multivariate Polynomial Ring in u, v, w
              over Integer Ring localized at (w,)
            sage: ER(GER.an_element())
            a*E3 + ((u/(-w))*a^2 + ((u^2 - v)/w)*a)*b + a - u
            sage: ER(GBR.an_element())
            (u^2 + v*w)/w

            sage: MBR = chbr.CubicHeckeRingOfDefinition(markov_trace_version=True)
            sage: MER = MBR.extension_ring()
            sage: ES = MER.as_splitting_algebra(); ES
            Splitting Algebra of T^2 + T + 1 with roots [E3, -E3 - 1]
              over Splitting Algebra of h^3 - u*h^2 + v*h - w
                with roots [a, b, -b - a + u]
              over Multivariate Polynomial Ring in u, v, w, s
              over Integer Ring localized at (s, w, v, u)
            sage: ES(MER.an_element())
            (((-1)/(-s))*a)*E3 + ((u/(-w))*a^2 + ((u^2 - v)/w)*a)*b + a - u
            sage: ES(MBR.an_element())
            (u^2*s + v*w)/(w*s)
        """
        if self._splitting_algebra is not None:
            verbose("End (short)", level=2)
            return self._splitting_algebra

        if self._ring_of_definition is None:
            verbose("constructing generic base ring", level=2)
            self._ring_of_definition = CubicHeckeRingOfDefinition()

        markov = self._is_markov_trace_version()

        BR = self._ring_of_definition
        root_names = list(self._names)
        var_s = None
        if markov:
            var_s = root_names.pop()  # s not needed as root
            a, b, c, s = self.gens()
        else:
            a, b, c = self.gens()

        root_names.pop()  # c not needed as root

        FSR = SplittingAlgebra(BR.cubic_equation(), root_names, warning=False)
        splitting_roots = FSR.splitting_roots()
        verbose('splitting roots %s' % splitting_roots, level=2)

        A, B, C = splitting_roots
        e3 = self.cyclotomic_generator()
        if var_s:
            fsr_s = FSR.scalar_base_ring().gens_dict()[var_s]
            S = self.create_specialization([A, B, C], im_writhe_parameter=fsr_s)
            # check of embedding fails in this case as long as the images of
            # ``iu`` and ``iv`` need to be invertible (see comment in
            # :meth:`__init__` of :class:`CubicHeckeRingOfDefinition`).
            map_back = S.hom([e3, b, a, a + b + c, a*b+a*c+b*c, a*b*c, s], check=False)
        else:
            S = self.create_specialization([A, B, C])
            map_back = S.hom([e3, b, a, a + b + c, a*b+a*c+b*c, a*b*c])
        self.register_coercion(map_back)
        self._splitting_algebra = S
        return self._splitting_algebra

    def field_embedding(self, characteristic=0):
        r"""
        Return a field embedding of ``self``.

        INPUT:

        - ``characteristic`` -- integer (default: `0`); the characteristic
          of the field

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: BR = chbr.CubicHeckeRingOfDefinition()
            sage: ER = BR.extension_ring()
            sage: ER.field_embedding()
            Ring morphism:
            From: Multivariate Laurent Polynomial Ring in a, b, c
                    over Splitting Algebra of x^2 + x + 1
                      with roots [e3, -e3 - 1]
                    over Integer Ring
            To:   Fraction Field of Multivariate Polynomial Ring in a, b, c
                    over Cyclotomic Field of order 3 and degree 2
            Defn: a |--> a
                  b |--> b
                  c |--> c
            with map of base ring

            sage: ER.field_embedding(characteristic=5)
            Ring morphism:
            From: Multivariate Laurent Polynomial Ring in a, b, c
                    over Splitting Algebra of x^2 + x + 1
                      with roots [e3, -e3 - 1]
                    over Integer Ring
            To:   Fraction Field of Multivariate Polynomial Ring in a, b, c
                    over Finite Field in a of size 5^2
            Defn: a |--> a
                  b |--> b
                  c |--> c
            with map of base ring

            sage: MER = ER.markov_trace_version()
            sage: MER.field_embedding()
            Ring morphism:
            From: Multivariate Laurent Polynomial Ring in a, b, c, s
                    over Splitting Algebra of x^2 + x + 1
                      with roots [e3, -e3 - 1]
                    over Integer Ring
            To:   Fraction Field of Multivariate Polynomial Ring in a, b, c, s
                    over Cyclotomic Field of order 3 and degree 2
            Defn: a |--> a
                  b |--> b
                  c |--> c
                  s |--> s
            with map of base ring
        """
        if characteristic == 0:
            from sage.rings.number_field.number_field import CyclotomicField
            C3 = CyclotomicField(3)
            E3 = C3.gen()
        else:
            if not ZZ(characteristic).is_prime():
                raise ValueError('characteristic must be a prime integer')
            from sage.rings.finite_rings.finite_field_constructor import GF
            from sage.misc.functional import cyclotomic_polynomial
            G = GF(characteristic)
            c3 = cyclotomic_polynomial(3).change_ring(G)
            C3 = c3.splitting_field('a')
            E3 = c3.change_ring(C3).roots()[0][0]

        B = self.base_ring()
        embBase = B.hom((E3,))
        if not C3.has_coerce_map_from(B):
            C3._unset_coercions_used()
            C3.register_coercion(embBase)
        P = C3[self.variable_names()]
        F = P.fraction_field()
        emb = self.hom((F(E3),) + F.gens())
        F = emb.codomain()
        if not F.has_coerce_map_from(self):
            F._unset_coercions_used()
            F.register_coercion(emb)
        return emb

    def markov_trace_version(self):
        r"""
        Return the Markov trace version of ``self``.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: ER = chbr.CubicHeckeExtensionRing('a, b, c')
            sage: ER.markov_trace_version()
            Multivariate Laurent Polynomial Ring in a, b, c, s
              over Splitting Algebra of x^2 + x + 1
                with roots [e3, -e3 - 1] over Integer Ring
        """
        if self._is_markov_trace_version():
            return self
        names = self.variable_names() + ('s',)
        return self.__class__(names=names, order=self._order, markov_trace_version=True)


################################################################################
# Ring of Definition
# ------------------------------------------------------------------------------
# Definition of the ring of definition for the cubic Hecke algebra as polynomial
# ring in 2 indeterminates over univariate Laurent polynomial ring over the
# integers. This is the most general ring over which the cubic Hecke algebra may
# be defined.
# ------------------------------------------------------------------------------
class CubicHeckeRingOfDefinition(Localization):
    r"""
    The *ring of definition* of the cubic Hecke algebra.

    It contains one invertible indeterminate (representing the product of the
    roots of the cubic equation) and two non invertible indeterminates.

    .. NOTE::

        We follow a suggestion by Ivan Marin in the name *ring of definition*.
        We avoid alternative names like *generic* or *universal* base ring
        as these have some issues. The first option could be misleading
        in view of the term *generic point* used in algebraic geometry, which
        would mean the function field in ``u, v, w``, here.

        The second option is problematic since the base ring itself is not a
        universal object. Rather, the universal object is the cubic Hecke algebra
        considered as a `\ZZ`-algebra including ``u, v, w`` as pairwise commuting
        indeterminates. From this point of view the base ring appears to be a
        subalgebra of this universal object generated by ``u, v, w``.

    INPUT:

    - ``names`` -- (default: ``'u,v,w'``) string containing the names of the
      indeterminates separated by ``,`` or a triple of strings each of which
      are the names of one of the three indeterminates
    - ``order`` -- string (default: ``'degrevlex'``); the term order; see also
      :class:`~sage.rings.polynomial.laurent_polynomial_ring.LaurentPolynomialRing_mpair`
    - ``ring_of_definition`` -- (optional) a :class:`CubicHeckeRingOfDefinition`
      to specify the generic cubic Hecke base ring over which ``self`` may be
      realized as splitting ring via the ``as_splitting_algebra`` method
    - ``markov_trace_version`` -- boolean (default: ``False``); if this is
      set to ``True`` then ``self`` contains one invertible indeterminate in
      addition which is meant to represent the writhe factor of a Markov trace
      on the cubic Hecke algebra and which default name is ``s``

    EXAMPLES::

        sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
        sage: BR = chbr.CubicHeckeRingOfDefinition()
        sage: u, v, w = BR.gens()
        sage: ele = 3*u*v-5*w**(-2)
        sage: ER = BR.extension_ring()
        sage: ER(ele)
        3*a^2*b + 3*a*b^2 + 3*a^2*c + 9*a*b*c + 3*b^2*c
        + 3*a*c^2 + 3*b*c^2 + (-5)*a^-2*b^-2*c^-2
        sage: phi1 = BR.hom( [4,3,1/1] )
        sage: phi1(ele)
        31

        sage: LL.<t> = LaurentPolynomialRing(ZZ)
        sage: phi2=BR.hom( [LL(4),LL(3),t] )
        sage: phi2(ele)
        -5*t^-2 + 36

        sage: BR.create_specialization( [E(5), E(7), E(3)] )
        Universal Cyclotomic Field
        sage: _(ele)
        -3*E(105) - 5*E(105)^2 - 5*E(105)^8 - 5*E(105)^11 - 5*E(105)^17
        - 5*E(105)^23 - 5*E(105)^26 - 5*E(105)^29 - 5*E(105)^32 - 5*E(105)^38
        - 5*E(105)^41 - 5*E(105)^44 - 5*E(105)^47 - 5*E(105)^53 - 5*E(105)^59
        - 5*E(105)^62 - 5*E(105)^68 - 8*E(105)^71 - 5*E(105)^74 - 5*E(105)^83
        - 5*E(105)^86 - 5*E(105)^89 - 5*E(105)^92 - 5*E(105)^101 - 5*E(105)^104
    """
    def __init__(self, names=('u', 'v', 'w', 's'), order='degrevlex', markov_trace_version=False):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: BR = chbr.CubicHeckeRingOfDefinition()
            sage: TestSuite(BR).run()
        """
        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------
        # Saving class-globals
        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------
        names = normalize_names_markov(names, markov_trace_version)
        if len(names) == 4:
            # invertible_names = names[2:4]  # s must be invertible, too
            #
            # because the formal Markov trace of the both basis elements
            # ``self.get_order()[568]`` (``KnotInfo.L8a7_1``) and
            # ``self.get_order()[596]`` (mirror image of ``KnotInfo.K9_46``)
            # have the indeterminate ``v`` in their denominator we need to have
            # all indeterminates invertible (``u`` as well for the mirror images)
            invertible_names = names
        else:
            invertible_names = names[2]

        self._order = order

        # ----------------------------------------------------------------------
        # ---------------------------------------------------------------------
        # Init of self
        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------
        base_ring = PolynomialRing(ZZ, names, order=order)
        Localization.__init__(self, base_ring, invertible_names)

        # ----------------------------------------------------------------------
        # Init of data used on demand
        # ----------------------------------------------------------------------
        self._mirror = None
        return

    ############################################################################
    # overloaded inherited methods
    ############################################################################
    def _defining_names(self):
        r"""
        Return the generators of ``self`` as the defining names.

        This method is cached in the parent class.
        This causes trouble if a second instance of ``self`` is used for another
        cubic Hecke algebra in the same session. To avoid this it is overloaded
        without ``cached_method`` decorator.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: BR = chbr.CubicHeckeRingOfDefinition()
            sage: BR._defining_names()
            (u, v, w)
        """
        return self.gens()

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: BR = chbr.CubicHeckeRingOfDefinition()
            sage: BR.an_element()                            # indirect doctest
            (u^2 + v*w)/w
            sage: MBR = chbr.CubicHeckeRingOfDefinition(markov_trace_version=True)
            sage: MBR.an_element()                           # indirect doctest
            (u^2*s + v*w)/(w*s)
        """
        u, v, w, *rem = self.gens()
        s = self.one()
        if rem:
            s = rem[0]
        return u**2/w+v/s

    ############################################################################
    # Local Methods
    ############################################################################
    def _is_markov_trace_version(self) -> bool:
        r"""
        Return whether ``self`` is the version containing the writhe parameter
        ``s`` for the Markov trace.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: BR = chbr.CubicHeckeRingOfDefinition()
            sage: BR._is_markov_trace_version()
            False
            sage: MBR = chbr.CubicHeckeRingOfDefinition(markov_trace_version=True)
            sage: MBR._is_markov_trace_version()
            True
        """
        return len(self.gens()) == 4

    ############################################################################
    # Global Methods
    ############################################################################
    def cubic_equation(self, var='h', as_coefficients=False):
        r"""
        Return the cubic equation over which the cubic Hecke algebra is defined.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: BR = chbr.CubicHeckeRingOfDefinition()
            sage: BR.cubic_equation()
            h^3 - u*h^2 + v*h - w
            sage: BR.cubic_equation(var='t')
            t^3 - u*t^2 + v*t - w
            sage: BR.cubic_equation(as_coefficients=True)
            [-w, v, -u, 1]
        """
        u, v, w, *rem = self.gens()
        cf = [-w, v, -u, 1]
        if as_coefficients:
            return cf
        P = PolynomialRing(self, var)

        return P(cf)

    def mirror_involution(self):
        r"""
        Return the involution of ``self`` corresponding to the involution of the
        cubic Hecke algebra (with the same name).

        This means that it maps the last generator of ``self`` to its inverse
        and both others to their product with the image of the former.

        From the cubic equation for a braid generator `\beta_i`:

        .. MATH::

            \beta_i^3 - u \beta_i^2 + v\beta_i -w = 0.

        One deduces the following cubic equation for `\beta_i^{-1}`:

        .. MATH::

            \beta_i^{-3} -\frac{v}{w} \beta_i^{-2} + \frac{u}{w}\beta_i^{-1}
            - \frac{1}{w} = 0.

        .. NOTE::

           The mirror involution of the braid group does not factor through
           the cubic Hecke algebra over its base ring, but it does if it is
           considered as `\ZZ`-algebra. The base ring elements are transformed
           by this automorphism.

        OUTPUT: the involution as automorphism of ``self``

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: BR = chbr.CubicHeckeRingOfDefinition()
            sage: BR.mirror_involution()
            Ring endomorphism of Multivariate Polynomial Ring in u, v, w
                                 over Integer Ring localized at (w,)
              Defn: u |--> v/w
                    v |--> u/w
                    w |--> 1/w
            sage: _(BR.an_element())
            (v^2 + u)/w

            sage: MBR = chbr.CubicHeckeRingOfDefinition(markov_trace_version=True)
            sage: MBR.mirror_involution()
            Ring endomorphism of Multivariate Polynomial Ring in u, v, w, s
                                 over Integer Ring localized at (s, w, v, u)
            Defn: u |--> v/w
            v |--> u/w
            w |--> 1/w
            s |--> 1/s
            sage: _(MBR.an_element())
            (v^2 + u*s)/w
        """
        if self._mirror is None:
            if self._is_markov_trace_version():
                u, v, w, s = self.gens()
                self._mirror = self.hom([v / w, u / w, ~w, ~s])
            else:
                u, v, w = self.gens()
                self._mirror = self.hom([v / w, u / w, ~w])
        return self._mirror

    def create_specialization(self, im_cubic_equation_parameters, im_writhe_parameter=None):
        r"""
        Return an appropriate Ring containing the elements from the list
        ``im_cubic_equation_parameters`` having a conversion map from ``self``
        mapping the cubic equation parameters of ``self`` to
        ``im_cubic_equation_parameters``.

        INPUT:

        - ``im_cubic_equation_parameters`` -- list or tuple of three ring
          elements such that there exists a ring homomorphism from the
          corresponding elements of ``self`` to them

        OUTPUT:

        A common parent containing the elements of ``im_cubic_equation_parameters``
        together with an inverse of the third element.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: BR = chbr.CubicHeckeRingOfDefinition()
            sage: t = BR.an_element(); t
            (u^2 + v*w)/w
            sage: Sp1 = BR.create_specialization([E(5), E(7), E(3)]); Sp1
            Universal Cyclotomic Field
            sage: Sp1(t)
            E(105) + E(105)^8 + E(105)^29 - E(105)^37 + E(105)^43 - E(105)^52
            + E(105)^64 - E(105)^67 + E(105)^71 - E(105)^82 + E(105)^92
            - E(105)^97

            sage: MBR = chbr.CubicHeckeRingOfDefinition(markov_trace_version=True)
            sage: MBR.create_specialization([E(5), E(7), E(3)], im_writhe_parameter=E(4))
            Universal Cyclotomic Field
            sage: u, v, w, s = MBR.gens()
            sage: Sp1(MBR(t)/s)
            E(420)^13 - E(420)^53 + E(420)^73 - E(420)^109 - E(420)^137
            - E(420)^221 + E(420)^253 - E(420)^277 + E(420)^313 - E(420)^361
            + E(420)^373 - E(420)^389

            sage: Z3 = CyclotomicField(3); E3=Z3.gen()
            sage: Sp2 = BR.create_specialization([E3, E3**2, 1]); Sp2
            Cyclotomic Field of order 3 and degree 2
            sage: Sp2(t)
            -2*zeta3 - 2
            sage: MBR.create_specialization([E3, E3**2, 1], im_writhe_parameter=2)
            Cyclotomic Field of order 3 and degree 2
            sage: Sp2(MBR(t)/s)
            -zeta3 - 1

            sage: Sp3 = BR.create_specialization([5, 7, 11]); Sp3
            Integer Ring localized at (11,)
            sage: Sp3(t)
            102/11
        """
        # ----------------------------------------------------------------------
        # setting the base_ring  according to the cubic_equation_parameters
        # ----------------------------------------------------------------------
        if isinstance(im_cubic_equation_parameters, tuple):
            im_cubic_equation_parameters = list(im_cubic_equation_parameters)

        if not isinstance(im_cubic_equation_parameters, list):
            raise TypeError('cubic_equation_parameters must be a list of three elements')

        if len(im_cubic_equation_parameters) != 3:
            raise ValueError('there must be exactly three cubic_equation_parameters')

        gens = self.gens()
        num_gens = len(gens)
        if im_writhe_parameter:
            if num_gens < 4:
                raise ValueError('im_writhe_parameter only possible for Markov-trace extension')
            im_gens = im_cubic_equation_parameters + [im_writhe_parameter]
            u, v, w, s = im_gens
        else:
            if num_gens == 4:
                raise ValueError('im_writhe_parameter must be given for Markov-trace extension')
            im_gens = im_cubic_equation_parameters
            u, v, w = im_gens

        image_ring = None
        image_ring_map = None
        image_ring_base = w.parent()

        # ----------------------------------------------------------------------
        # short exit on trivial invocation
        # ----------------------------------------------------------------------
        if image_ring_base is self and im_gens == gens:
            return self

        image_ring = get_coercion_model().common_parent(*(im_gens))

        # ----------------------------------------------------------------------
        # make sure that the inverse of w belongs to image_ring
        # ----------------------------------------------------------------------
        try:
            image_ring = image_ring.localization(w)
        except ValueError:
            pass

        im_gens = [image_ring(para) for para in im_gens]

        verbose('common parent of parameters and inverses: %s' % image_ring, level=2)

        try:
            image_ring_map = self.hom(im_gens, codomain=image_ring)
        except ValueError:
            image_ring_map = self.hom(im_gens, codomain=image_ring, check=False)
            verbose('Warning: check failed for embedding as ring morphism')

        register_ring_hom(image_ring_map)
        return image_ring

    # --------------------------------------------------------------------------
    # Definition of the generic extension ring for the cubic hecke algebra as
    # Laurent polynomial ring in 3 indeterminates over cyclotomic field of order
    # 3. The generic extension ring guarantees semisimplicity of the cubic Hecke
    # algebra.
    # --------------------------------------------------------------------------
    @cached_method
    def extension_ring(self, names=('a', 'b', 'c', 's')):
        r"""
        Return the generic extension ring attached to ``self``.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: BR = chbr.CubicHeckeRingOfDefinition()
            sage: BR.extension_ring()
            Multivariate Laurent Polynomial Ring in a, b, c
              over Splitting Algebra of x^2 + x + 1
                with roots [e3, -e3 - 1]
              over Integer Ring
        """
        markov = self._is_markov_trace_version()
        return CubicHeckeExtensionRing(names, ring_of_definition=self, markov_trace_version=markov)

    def markov_trace_version(self):
        r"""
        Return the extension of the ring of definition needed to treat the
        formal Markov traces.

        This appends an additional variable ``s`` to measure the writhe
        of knots and makes ``u`` and ``v`` invertible.

        EXAMPLES::

            sage: from sage.algebras.hecke_algebras import cubic_hecke_base_ring as chbr
            sage: GBR = chbr.CubicHeckeRingOfDefinition()
            sage: GBR.markov_trace_version()
            Multivariate Polynomial Ring in u, v, w, s
              over Integer Ring localized at (s, w, v, u)
        """
        if self._is_markov_trace_version():
            return self
        names = self.base_ring().variable_names() + ('s',)
        return self.__class__(names=names, order=self._order, markov_trace_version=True)

    def specialize_homfly(self):
        r"""
        Return a map to the two variable Laurent polynomial ring that is
        the parent of the HOMFLY-PT polynomial.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: CHA2 = algebras.CubicHecke(2)
            sage: K5_1 = KnotInfo.K5_1.link()
            sage: br = CHA2(K5_1.braid())
            sage: mt = br.formal_markov_trace()
            sage: MT = mt.base_ring()
            sage: f = MT.specialize_homfly(); f
            Composite map:
              From: Multivariate Polynomial Ring in u, v, w, s over Integer Ring
                    localized at (s, w, v, u)
              To:   Multivariate Laurent Polynomial Ring in L, M over Integer Ring
              Defn:   Ring morphism:
                      From: Multivariate Polynomial Ring in u, v, w, s
                            over Integer Ring localized at (s, w, v, u)
                      To:   Multivariate Polynomial Ring in L, M
                            over Integer Ring localized at (M, M - 1, L)
                      Defn: u |--> -M + 1
                            v |--> -M + 1
                            w |--> 1
                            s |--> L
                    then
                      Conversion map:
                      From: Multivariate Polynomial Ring in L, M
                            over Integer Ring localized at (M, M - 1, L)
                      To:   Multivariate Laurent Polynomial Ring in L, M
                            over Integer Ring
            sage: sup = mt.support()
            sage: h1 = sum(f(mt.coefficient(b)) * b.regular_homfly_polynomial() for b in sup)
            sage: L, M = f.codomain().gens()
            sage: h2 = K5_1.homfly_polynomial()
            sage: h1*L**(-5) == h2  # since the writhe of K5_1 is 5
            True
        """
        if not self._is_markov_trace_version():
            raise ValueError('Functionality available for Markov trace version, only')
        from sage.knots.link import Link
        H = Link([]).homfly_polynomial().parent()
        L, M = H.gens()
        HL = H.localization(1 - M)
        u = HL(1 - M)
        phi = self.hom((u, u, HL.one(), HL(L)))
        inc = H.convert_map_from(HL)
        return inc * phi

    def specialize_kauffman(self):
        r"""
        Return a map to the two variable Laurent polynomial ring that is
        the parent of the Kauffman polynomial.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: CHA2 = algebras.CubicHecke(2)
            sage: K5_1 = KnotInfo.K5_1.link()
            sage: br = CHA2(K5_1.braid())
            sage: mt = br.formal_markov_trace()
            sage: MT = mt.base_ring()
            sage: f = MT.specialize_kauffman(); f
            Composite map:
              From: Multivariate Polynomial Ring in u, v, w, s over Integer Ring
                    localized at (s, w, v, u)
              To:   Multivariate Laurent Polynomial Ring in a, z over Integer Ring
              Defn:   Ring morphism:
                      From: Multivariate Polynomial Ring in u, v, w, s
                            over Integer Ring localized at (s, w, v, u)
                      To:   Multivariate Polynomial Ring in a, z
                            over Integer Ring localized at (z, a, a + z, a*z + 1)
                      Defn: u |--> (a*z + 1)/a
                            v |--> (a + z)/a
                            w |--> 1/a
                            s |--> a
                    then
                      Conversion map:
                      From: Multivariate Polynomial Ring in a, z over Integer Ring
                            localized at (z, a, a + z, a*z + 1)
                      To:   Multivariate Laurent Polynomial Ring in a, z
                            over Integer Ring
            sage: sup = mt.support()
            sage: k1 = sum(f(mt.coefficient(b)) * b.regular_kauffman_polynomial() for b in sup)
            sage: a, z = f.codomain().gens()
            sage: k2 = KnotInfo.K5_1.kauffman_polynomial()
            sage: k1*a**(-5) == k2  # since the writhe of K5_1 is 5
            True
        """
        if not self._is_markov_trace_version():
            raise ValueError('Functionality available for Markov trace version, only')
        from sage.knots.knotinfo import KnotInfo
        K = KnotInfo.L2a1_1.kauffman_polynomial().parent()
        a, z = K.gens()
        ku = z * a + 1
        kv = z + a
        KL = K.localization((ku, kv))
        u = KL(ku / a)
        v = KL(kv / a)
        phi = self.hom((u, v, KL(~a), KL(a)))
        inc = K.convert_map_from(KL)
        return inc * phi

    def specialize_links_gould(self):
        r"""
        Return a map to the two variable Laurent polynomial ring that is
        the parent of the Links-Gould polynomial.

        EXAMPLES::

            sage: from sage.knots.knotinfo import KnotInfo
            sage: CHA2 = algebras.CubicHecke(2)
            sage: K5_1 = KnotInfo.K5_1.link()
            sage: br = CHA2(K5_1.braid())
            sage: mt = br.formal_markov_trace()
            sage: MT = mt.base_ring()
            sage: f = MT.specialize_links_gould(); f
            Composite map:
              From: Multivariate Polynomial Ring in u, v, w, s over Integer Ring
                    localized at (s, w, v, u)
              To:   Multivariate Laurent Polynomial Ring in t0, t1 over Integer Ring
              Defn:   Ring morphism:
                      From: Multivariate Polynomial Ring in u, v, w, s
                            over Integer Ring localized at (s, w, v, u)
                      To:   Multivariate Polynomial Ring in t0, t1 over Integer Ring
                            localized at (t1, t0, t0 + t1 - 1, t0*t1 - t0 - t1)
                      Defn: u |--> t0 + t1 - 1
                            v |--> t0*t1 - t0 - t1
                            w |--> -t0*t1
                            s |--> 1
                    then
                      Conversion map:
                      From: Multivariate Polynomial Ring in t0, t1 over Integer Ring
                            localized at (t1, t0, t0 + t1 - 1, t0*t1 - t0 - t1)
                      To:   Multivariate Laurent Polynomial Ring in t0, t1 over Integer Ring
            sage: sup = mt.support()
            sage: sum(f(mt.coefficient(b)) * b.links_gould_polynomial() for b in sup)
            -t0^4*t1 - t0^3*t1^2 - t0^2*t1^3 - t0*t1^4 + t0^4 + 2*t0^3*t1 + 2*t0^2*t1^2
            + 2*t0*t1^3 + t1^4 - t0^3 - 2*t0^2*t1 - 2*t0*t1^2 - t1^3 + t0^2 + 2*t0*t1
            + t1^2 - t0 - t1 + 1
        """
        if not self._is_markov_trace_version():
            raise ValueError('Functionality available for Markov trace version, only')
        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
        L = LaurentPolynomialRing(ZZ, 't0, t1')
        t0, t1 = L.gens()
        lu = t0 + t1 - 1
        lv = t0 * t1 - t0 - t1
        lw = -t0 * t1
        LL = L.localization((lu, lv))
        u = LL(lu)
        v = LL(lv)
        w = LL(lw)
        phi = self.hom((u, v, w, LL.one()))
        inc = L.convert_map_from(LL)
        return inc * phi
