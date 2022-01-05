"""
Graded Rings of Modular Forms

This module contains functions to find generators for the graded ring of
modular forms of given level.

AUTHORS:

- William Stein (2007-08-24): first version
"""
#*****************************************************************************
#       Copyright (C) 2007 William Stein
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.richcmp import richcmp_method, richcmp
from sage.rings.all import Integer, QQ, ZZ
from sage.misc.all import prod
from sage.misc.verbose import verbose
from sage.misc.cachefunc import cached_method
from sage.modular.arithgroup.all import Gamma0, is_CongruenceSubgroup
from .constructor                 import ModularForms
from sage.structure.sage_object  import SageObject
from random import shuffle

def _span_of_forms_in_weight(forms, weight, prec, stop_dim=None, use_random=False):
    r"""
    Utility function. Given a nonempty list of pairs ``(k,f)``, where `k` is an
    integer and `f` is a power series, and a weight l, return all weight l
    forms obtained by multiplying together the given forms.

    INPUT:

    - ``forms`` -- list of pairs `(k, f)` with k an integer and f a power
      series (all over the same base ring)
    - ``weight`` -- an integer
    - ``prec`` -- an integer (less than or equal to the precision of all the
      forms in ``forms``) -- precision to use in power series computations.
    - ``stop_dim`` -- an integer: stop as soon as we have enough forms to span
      a submodule of this rank (a saturated one if the base ring is `\ZZ`).
      Ignored if ``use_random`` is False.
    - ``use_random`` -- which algorithm to use. If True, tries random products
      of the generators of the appropriate weight until a large enough
      submodule is found (determined by ``stop_dim``). If False, just tries
      everything.

    Note that if the given forms do generate the whole space, then
    ``use_random=True`` will often be quicker (particularly if the weight is
    large); but if the forms don't generate, the randomized algorithm is no
    help and will actually be substantially slower, because it needs to do
    repeated echelon form calls to check if vectors are in a submodule, while
    the non-randomized algorithm just echelonizes one enormous matrix at the
    end.

    EXAMPLES::

        sage: import sage.modular.modform.find_generators as f
        sage: forms = [(4, 240*eisenstein_series_qexp(4,5)), (6,504*eisenstein_series_qexp(6,5))]
        sage: f._span_of_forms_in_weight(forms, 12, prec=5)
        Vector space of degree 5 and dimension 2 over Rational Field
        Basis matrix:
        [        1         0    196560  16773120 398034000]
        [        0         1       -24       252     -1472]
        sage: f._span_of_forms_in_weight(forms, 24, prec=5)
        Vector space of degree 5 and dimension 3 over Rational Field
        Basis matrix:
        [          1           0           0    52416000 39007332000]
        [          0           1           0      195660    12080128]
        [          0           0           1         -48        1080]
        sage: ModularForms(1, 24).q_echelon_basis(prec=5)
        [
        1 + 52416000*q^3 + 39007332000*q^4 + O(q^5),
        q + 195660*q^3 + 12080128*q^4 + O(q^5),
        q^2 - 48*q^3 + 1080*q^4 + O(q^5)
        ]

    Test the alternative randomized algorithm::

        sage: f._span_of_forms_in_weight(forms, 24, prec=5, use_random=True, stop_dim=3)
        Vector space of degree 5 and dimension 3 over Rational Field
        Basis matrix:
        [          1           0           0    52416000 39007332000]
        [          0           1           0      195660    12080128]
        [          0           0           1         -48        1080]
    """
    t = verbose('multiplying forms up to weight %s'%weight)
    # Algorithm: run through the monomials of the appropriate weight, and build
    # up the vector space they span.

    n = len(forms)
    R = forms[0][1].base_ring()
    V = R ** prec
    W = V.zero_submodule()
    shortforms = [f[1].truncate_powerseries(prec) for f in forms]

    # List of weights
    from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
    wts = list(WeightedIntegerVectors(weight, [f[0] for f in forms]))
    t = verbose("calculated weight list", t)
    N = len(wts)

    if use_random:
        if stop_dim is None:
            raise ValueError("stop_dim must be provided if use_random is True")
        shuffle(wts)

        for c in range(N):
            w = V(prod(shortforms[i]**wts[c][i] for i in range(n)).padded_list(prec))
            if w in W:
                continue
            W = V.span(list(W.gens()) + [w])
            if stop_dim and W.rank() == stop_dim:
                if R != ZZ or W.index_in_saturation() == 1:
                    verbose("Succeeded after %s of %s" % (c, N), t)
                    return W
        verbose("Nothing worked", t)
        return W
    else:
        G = [V(prod(forms[i][1]**c[i] for i in range(n)).padded_list(prec)) for c in wts]
        t = verbose('found %s candidates' % N, t)
        W = V.span(G)
        verbose('span has dimension %s' % W.rank(), t)
        return W

def find_generators(*args):
    r"""
    This function, which existed in earlier versions of Sage, has now been
    replaced by the :meth:`~ModularFormsRing.generators` method of
    ModularFormsRing objects.

    EXAMPLES::

        sage: from sage.modular.modform.find_generators import find_generators
        sage: find_generators()
        Traceback (most recent call last):
        ...
        NotImplementedError: find_generators has been removed -- use ModularFormsRing.generators()
    """
    raise NotImplementedError("find_generators has been removed -- use ModularFormsRing.generators()")

def basis_for_modform_space(*args):
    r"""
    This function, which existed in earlier versions of Sage, has now been
    replaced by the :meth:`~ModularFormsRing.q_expansion_basis` method of
    ModularFormsRing objects.

    EXAMPLES::

        sage: from sage.modular.modform.find_generators import basis_for_modform_space
        sage: basis_for_modform_space()
        Traceback (most recent call last):
        ...
        NotImplementedError: basis_for_modform_space has been removed -- use ModularFormsRing.q_expansion_basis()
    """
    raise NotImplementedError("basis_for_modform_space has been removed -- use ModularFormsRing.q_expansion_basis()")


@richcmp_method
class ModularFormsRing(SageObject):

    def __init__(self, group, base_ring=QQ):
        r"""
        The ring of modular forms (of weights 0 or at least 2) for a congruence
        subgroup of `{\rm SL}_2(\ZZ)`, with coefficients in a specified base ring.

        INPUT:

        - ``group`` -- a congruence subgroup of `{\rm SL}_2(\ZZ)`, or a
          positive integer `N` (interpreted as `\Gamma_0(N)`)

        - ``base_ring`` (ring, default: `\QQ`) -- a base ring, which should be
          `\QQ`, `\ZZ`, or the integers mod `p` for some prime `p`.

        EXAMPLES::

            sage: ModularFormsRing(Gamma1(13))
            Ring of modular forms for Congruence Subgroup Gamma1(13) with coefficients in Rational Field
            sage: m = ModularFormsRing(4); m
            Ring of modular forms for Congruence Subgroup Gamma0(4) with coefficients in Rational Field
            sage: m.modular_forms_of_weight(2)
            Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(4) of weight 2 over Rational Field
            sage: m.modular_forms_of_weight(10)
            Modular Forms space of dimension 6 for Congruence Subgroup Gamma0(4) of weight 10 over Rational Field
            sage: m == loads(dumps(m))
            True
            sage: m.generators()
            [(2, 1 + 24*q^2 + 24*q^4 + 96*q^6 + 24*q^8 + O(q^10)),
            (2, q + 4*q^3 + 6*q^5 + 8*q^7 + 13*q^9 + O(q^10))]
            sage: m.q_expansion_basis(2,10)
            [1 + 24*q^2 + 24*q^4 + 96*q^6 + 24*q^8 + O(q^10),
             q + 4*q^3 + 6*q^5 + 8*q^7 + 13*q^9 + O(q^10)]
            sage: m.q_expansion_basis(3,10)
            []
            sage: m.q_expansion_basis(10,10)
            [1 + 10560*q^6 + 3960*q^8 + O(q^10),
             q - 8056*q^7 - 30855*q^9 + O(q^10),
             q^2 - 796*q^6 - 8192*q^8 + O(q^10),
             q^3 + 66*q^7 + 832*q^9 + O(q^10),
             q^4 + 40*q^6 + 528*q^8 + O(q^10),
             q^5 + 20*q^7 + 190*q^9 + O(q^10)]

        TESTS:

        Check that :trac:`15037` is fixed::

            sage: ModularFormsRing(3.4)
            Traceback (most recent call last):
            ...
            ValueError: Group (=3.40000000000000) should be a congruence subgroup
            sage: ModularFormsRing(Gamma0(2), base_ring=PolynomialRing(ZZ,x))
            Traceback (most recent call last):
            ...
            ValueError: Base ring (=Univariate Polynomial Ring in x over Integer Ring) should be QQ, ZZ or a finite prime field
        """
        if isinstance(group, (int, Integer)):
            group = Gamma0(group)
        elif not is_CongruenceSubgroup(group):
            raise ValueError("Group (=%s) should be a congruence subgroup" % group)

        if base_ring != ZZ and not base_ring.is_prime_field():
            raise ValueError("Base ring (=%s) should be QQ, ZZ or a finite prime field" % base_ring)

        self.__group = group
        self.__base_ring = base_ring
        self.__cached_maxweight = ZZ(-1)
        self.__cached_gens = []
        self.__cached_cusp_maxweight = ZZ(-1)
        self.__cached_cusp_gens = []

    def group(self):
        r"""
        Return the congruence subgroup for which this is the ring of modular forms.

        EXAMPLES::

            sage: R = ModularFormsRing(Gamma1(13))
            sage: R.group() is Gamma1(13)
            True
        """
        return self.__group

    def base_ring(self):
        r"""
        Return the coefficient ring of this modular forms ring.

        EXAMPLES::

            sage: ModularFormsRing(Gamma1(13)).base_ring()
            Rational Field
            sage: ModularFormsRing(Gamma1(13), base_ring = ZZ).base_ring()
            Integer Ring
        """
        return self.__base_ring

    def __richcmp__(self, other, op):
        r"""
        Compare self to other.

        Rings are equal if and only if their groups and base rings are.

        EXAMPLES::

            sage: ModularFormsRing(3) == 3
            False
            sage: ModularFormsRing(Gamma0(3)) == ModularFormsRing(Gamma0(7))
            False
            sage: ModularFormsRing(Gamma0(3)) == ModularFormsRing(Gamma0(3))
            True
        """
        if not isinstance(other, ModularFormsRing):
            return NotImplemented

        return richcmp((self.group(), self.base_ring()),
                       (other.group(), other.base_ring()), op)

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: ModularFormsRing(Gamma0(13))._repr_()
            'Ring of modular forms for Congruence Subgroup Gamma0(13) with coefficients in Rational Field'
            sage: ModularFormsRing(Gamma1(13), base_ring=ZZ)._repr_()
            'Ring of modular forms for Congruence Subgroup Gamma1(13) with coefficients in Integer Ring'
        """
        return "Ring of modular forms for %s with coefficients in %s" % (self.group(), self.base_ring())

    def modular_forms_of_weight(self, weight):
        """
        Return the space of modular forms on this group of the given weight.

        EXAMPLES::

            sage: R = ModularFormsRing(13)
            sage: R.modular_forms_of_weight(10)
            Modular Forms space of dimension 11 for Congruence Subgroup Gamma0(13) of weight 10 over Rational Field
            sage: ModularFormsRing(Gamma1(13)).modular_forms_of_weight(3)
            Modular Forms space of dimension 20 for Congruence Subgroup Gamma1(13) of weight 3 over Rational Field
        """
        return ModularForms(self.group(), weight)

    def generators(self, maxweight=8, prec=10, start_gens=[], start_weight=2):
        r"""
        If `R` is the base ring of self, then this function calculates a set of
        modular forms which generate the `R`-algebra of all modular forms of
        weight up to ``maxweight`` with coefficients in `R`.

        INPUT:

        - ``maxweight`` (integer, default: 8) -- check up to this weight for
          generators

        - ``prec`` (integer, default: 10) -- return `q`-expansions to this
          precision

        - ``start_gens`` (list, default: ``[]``) -- list of pairs `(k, f)`, or
          triples `(k, f, F)`, where:

          - `k` is an integer,
          - `f` is the `q`-expansion of a modular form of weight `k`, as a power series over the base ring of self,
          - `F` (if provided) is a modular form object corresponding to F.

          If this list is nonempty, we find a minimal generating set containing
          these forms. If `F` is not supplied, then `f` needs to have
          sufficiently large precision (an error will be raised if this is not
          the case); otherwise, more terms will be calculated from the modular
          form object `F`.

        - ``start_weight`` (integer, default: 2) -- calculate the graded
          subalgebra of forms of weight at least ``start_weight``.

        OUTPUT:

        a list of pairs (k, f), where f is the q-expansion to precision
        ``prec`` of a modular form of weight k.

        .. SEEALSO::

            :meth:`gen_forms`, which does exactly the same thing, but returns
            Sage modular form objects rather than bare power series, and keeps
            track of a lifting to characteristic 0 when the base ring is a
            finite field.

        .. note::

            If called with the default values of ``start_gens`` (an empty list)
            and ``start_weight`` (2), the values will be cached for re-use on
            subsequent calls to this function. (This cache is shared with
            :meth:`gen_forms`). If called with non-default values for these
            parameters, caching will be disabled.

        EXAMPLES::

            sage: ModularFormsRing(SL2Z).generators()
            [(4, 1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + 60480*q^6 + 82560*q^7 + 140400*q^8 + 181680*q^9 + O(q^10)), (6, 1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 - 4058208*q^6 - 8471232*q^7 - 17047800*q^8 - 29883672*q^9 + O(q^10))]
            sage: s = ModularFormsRing(SL2Z).generators(maxweight=5, prec=3); s
            [(4, 1 + 240*q + 2160*q^2 + O(q^3))]
            sage: s[0][1].parent()
            Power Series Ring in q over Rational Field

            sage: ModularFormsRing(1).generators(prec=4)
            [(4, 1 + 240*q + 2160*q^2 + 6720*q^3 + O(q^4)), (6, 1 - 504*q - 16632*q^2 - 122976*q^3 + O(q^4))]
            sage: ModularFormsRing(2).generators(prec=12)
            [(2, 1 + 24*q + 24*q^2 + 96*q^3 + 24*q^4 + 144*q^5 + 96*q^6 + 192*q^7 + 24*q^8 + 312*q^9 + 144*q^10 + 288*q^11 + O(q^12)), (4, 1 + 240*q^2 + 2160*q^4 + 6720*q^6 + 17520*q^8 + 30240*q^10 + O(q^12))]
            sage: ModularFormsRing(4).generators(maxweight=2, prec=20)
            [(2, 1 + 24*q^2 + 24*q^4 + 96*q^6 + 24*q^8 + 144*q^10 + 96*q^12 + 192*q^14 + 24*q^16 + 312*q^18 + O(q^20)), (2, q + 4*q^3 + 6*q^5 + 8*q^7 + 13*q^9 + 12*q^11 + 14*q^13 + 24*q^15 + 18*q^17 + 20*q^19 + O(q^20))]

        Here we see that for ``\Gamma_0(11)`` taking a basis of forms in weights 2
        and 4 is enough to generate everything up to weight 12 (and probably
        everything else).::

            sage: v = ModularFormsRing(11).generators(maxweight=12)
            sage: len(v)
            3
            sage: [k for k, _ in v]
            [2, 2, 4]
            sage: dimension_modular_forms(11,2)
            2
            sage: dimension_modular_forms(11,4)
            4

        For congruence subgroups not containing -1, we miss out some forms since we
        can't calculate weight 1 forms at present, but we can still find generators
        for the ring of forms of weight `\ge 2`::

            sage: ModularFormsRing(Gamma1(4)).generators(prec=10, maxweight=10)
            [(2, 1 + 24*q^2 + 24*q^4 + 96*q^6 + 24*q^8 + O(q^10)),
            (2, q + 4*q^3 + 6*q^5 + 8*q^7 + 13*q^9 + O(q^10)),
            (3, 1 + 12*q^2 + 64*q^3 + 60*q^4 + 160*q^6 + 384*q^7 + 252*q^8 + O(q^10)),
            (3, q + 4*q^2 + 8*q^3 + 16*q^4 + 26*q^5 + 32*q^6 + 48*q^7 + 64*q^8 + 73*q^9 + O(q^10))]

        Using different base rings will change the generators::

            sage: ModularFormsRing(Gamma0(13)).generators(maxweight=12, prec=4)
            [(2, 1 + 2*q + 6*q^2 + 8*q^3 + O(q^4)), (4, 1 + O(q^4)), (4, q + O(q^4)), (4, q^2 + O(q^4)), (4, q^3 + O(q^4)), (6, 1 + O(q^4)), (6, q + O(q^4))]
            sage: ModularFormsRing(Gamma0(13),base_ring=ZZ).generators(maxweight=12, prec=4)
            [(2, 1 + 2*q + 6*q^2 + 8*q^3 + O(q^4)), (4, q + 4*q^2 + 10*q^3 + O(q^4)), (4, 2*q^2 + 5*q^3 + O(q^4)), (4, q^2 + O(q^4)), (4, -2*q^3 + O(q^4)), (6, O(q^4)), (6, O(q^4)), (12, O(q^4))]
            sage: [k for k,f in ModularFormsRing(1, QQ).generators(maxweight=12)]
            [4, 6]
            sage: [k for k,f in ModularFormsRing(1, ZZ).generators(maxweight=12)]
            [4, 6, 12]
            sage: [k for k,f in ModularFormsRing(1, Zmod(5)).generators(maxweight=12)]
            [4, 6]
            sage: [k for k,f in ModularFormsRing(1, Zmod(2)).generators(maxweight=12)]
            [4, 6, 12]

        An example where ``start_gens`` are specified::

            sage: M = ModularForms(11, 2); f = (M.0 + M.1).qexp(8)
            sage: ModularFormsRing(11).generators(start_gens = [(2, f)])
            Traceback (most recent call last):
            ...
            ValueError: Requested precision cannot be higher than precision of approximate starting generators!
            sage: f = (M.0 + M.1).qexp(10); f
            1 + 17/5*q + 26/5*q^2 + 43/5*q^3 + 94/5*q^4 + 77/5*q^5 + 154/5*q^6 + 86/5*q^7 + 36*q^8 + 146/5*q^9 + O(q^10)
            sage: ModularFormsRing(11).generators(start_gens = [(2, f)])
            [(2, 1 + 17/5*q + 26/5*q^2 + 43/5*q^3 + 94/5*q^4 + 77/5*q^5 + 154/5*q^6 + 86/5*q^7 + 36*q^8 + 146/5*q^9 + O(q^10)), (2, 1 + 12*q^2 + 12*q^3 + 12*q^4 + 12*q^5 + 24*q^6 + 24*q^7 + 36*q^8 + 36*q^9 + O(q^10)), (4, 1 + O(q^10))]
        """
        sgs = []
        for x in start_gens:
            if len(x) == 2:
                if x[1].prec() < prec:
                    raise ValueError("Requested precision cannot be higher than precision of approximate starting generators!")
                sgs.append((x[0], x[1], None))
            else:
                sgs.append(x)

        G = self._find_generators(maxweight, tuple(sgs), start_weight)

        ret = []
        # Returned generators may be a funny mixture of precisions if start_gens has been used.
        for k, f, F in G:
            if f.prec() < prec:
                f = F.qexp(prec).change_ring(self.base_ring())
            else:
                f = f.truncate_powerseries(prec)
            ret.append((k, f))

        return ret


    def gen_forms(self, maxweight=8, start_gens=[], start_weight=2):
        r"""
        This function calculates a list of modular forms generating this ring
        (as an algebra over the appropriate base ring). It differs from
        :meth:`generators` only in that it returns Sage modular form objects,
        rather than bare `q`-expansions; and if the base ring is a finite
        field, the modular forms returned will be forms in characteristic 0
        with integral `q`-expansions whose reductions modulo `p` generate the
        ring of modular forms mod `p`.

        INPUT:

        - ``maxweight`` (integer, default: 8) -- calculate forms generating all
          forms up to this weight.

        - ``start_gens`` (list, default: ``[]``) -- a list of modular forms. If
          this list is nonempty, we find a minimal generating set containing
          these forms.

        - ``start_weight`` (integer, default: 2) -- calculate the graded
          subalgebra of forms of weight at least ``start_weight``.

        .. note::

            If called with the default values of ``start_gens`` (an empty list)
            and ``start_weight`` (2), the values will be cached for re-use on
            subsequent calls to this function. (This cache is shared with
            :meth:`generators`). If called with non-default values for these
            parameters, caching will be disabled.

        EXAMPLES::

            sage: A = ModularFormsRing(Gamma0(11), Zmod(5)).gen_forms(); A
            [1 + 12*q^2 + 12*q^3 + 12*q^4 + 12*q^5 + O(q^6), q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6), q - 9*q^4 - 10*q^5 + O(q^6)]
            sage: A[0].parent()
            Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
        """
        sgs = tuple( (F.weight(), None, F) for F in start_gens )
        G = self._find_generators(maxweight, sgs, start_weight)
        return [F for k,f,F in G]

    def _find_generators(self, maxweight, start_gens, start_weight):
        r"""
        For internal use. This function is called by :meth:`generators` and
        :meth:`gen_forms`: it returns a list of triples `(k, f, F)` where `F`
        is a modular form of weight `k` and `f` is its `q`-expansion coerced
        into the base ring of self.

        INPUT:

        - maxweight: maximum weight to try
        - start_weight: minimum weight to try
        - start_gens: a sequence of tuples of the form `(k, f, F)`, where `F` is a
          modular form of weight `k` and `f` is its `q`-expansion coerced into
          ``self.base_ring()`. Either (but not both) of `f` and `F` may be
          None.

        OUTPUT:

        a list of tuples, formatted as with ``start_gens``.

        EXAMPLES::

            sage: R = ModularFormsRing(Gamma1(4))
            sage: R._find_generators(8, (), 2)
            [(2, 1 + 24*q^2 + 24*q^4 + 96*q^6 + 24*q^8 + O(q^9), 1 + 24*q^2 + 24*q^4 + O(q^6)), (2, q + 4*q^3 + 6*q^5 + 8*q^7 + O(q^9), q + 4*q^3 + 6*q^5 + O(q^6)), (3, 1 + 12*q^2 + 64*q^3 + 60*q^4 + 160*q^6 + 384*q^7 + 252*q^8 + O(q^9), 1 + 12*q^2 + 64*q^3 + 60*q^4 + O(q^6)), (3, q + 4*q^2 + 8*q^3 + 16*q^4 + 26*q^5 + 32*q^6 + 48*q^7 + 64*q^8 + O(q^9), q + 4*q^2 + 8*q^3 + 16*q^4 + 26*q^5 + O(q^6))]
        """
        default_params = (start_gens == () and start_weight == 2)

        if default_params and self.__cached_maxweight != -1:
            verbose("Already know generators up to weight %s -- using those" % self.__cached_maxweight)

            if self.__cached_maxweight >= maxweight:
                return [(k, f, F) for k, f, F in self.__cached_gens if k <= maxweight]

            start_gens = self.__cached_gens
            start_weight = self.__cached_maxweight + 1

        if self.group().is_even():
            increment = 2
        else:
            increment = 1

        working_prec = self.modular_forms_of_weight(maxweight).sturm_bound()

        # parse the list of start gens
        G = []
        for x in start_gens:
            k, f, F = x
            if F is None and f.prec() < working_prec:
                raise ValueError("Need start gens to precision at least %s" % working_prec)
            elif f is None or f.prec() < working_prec:
                f = F.qexp(working_prec).change_ring(self.base_ring())
            G.append((k, f, F))

        k = start_weight
        if increment == 2 and (k % 2) == 1:
            k += 1

        while k <= maxweight:

            if self.modular_forms_of_weight(k).dimension() == 0:
                k += increment
                continue

            verbose('Looking at k = %s' % k)
            M = self.modular_forms_of_weight(k)

            # 1. Multiply together all forms in G that give an element
            #    of M.
            if G:
                F = _span_of_forms_in_weight(G, k, M.sturm_bound(), None, False)
            else:
                F = (self.base_ring() ** M.sturm_bound()).zero_submodule()

            # 2. If the dimension of the span of the result is equal
            #    to the dimension of M, increment k.
            if F.rank() == M.dimension():
                if self.base_ring().is_field() or F.index_in_saturation() == 1:
                    # TODO: Do something clever if the submodule's of the right
                    # rank but not saturated -- avoid triggering needless
                    # modular symbol computations.
                    verbose('Nothing new in weight %s' % k)
                    k += increment
                    continue

            # 3. If the dimension is less, compute a basis for G, and
            #    try adding basis elements of M into G.

            verbose("Known generators span a subspace of dimension %s of space of dimension %s" % (F.dimension(), M.dimension()))
            if self.base_ring() == ZZ:
                verbose("saturation index is %s" % F.index_in_saturation())

            t = verbose("Computing more modular forms at weight %s" % k)
            kprec = M.sturm_bound()
            if self.base_ring() == QQ:
                B = M.q_echelon_basis(working_prec)
            else:
                B = M.q_integral_basis(working_prec)
            t = verbose("done computing forms", t)
            V = F.ambient_module().submodule_with_basis([f.padded_list(kprec) for f in B])
            Q = V / F
            for q in Q.gens():
                try:
                    qc = V.coordinates(Q.lift(q))
                except AttributeError:
                    # work around a silly free module bug
                    qc = V.coordinates(q.lift())
                qcZZ = [ZZ(_) for _ in qc] # lift to ZZ so we can define F
                f = sum([B[i] * qcZZ[i] for i in range(len(B))])
                F = M(f)
                G.append((k, f.change_ring(self.base_ring()), F))

            verbose('added %s new generators' % Q.ngens(), t)
            k += increment

        if default_params:
            self.__cached_maxweight = maxweight
            self.__cached_gens = G

        return G

    @cached_method
    def q_expansion_basis(self, weight, prec=None, use_random=True):
        r"""
        Calculate a basis of q-expansions for the space of modular forms of the
        given weight for this group, calculated using the ring generators given
        by ``find_generators``.

        INPUT:

        - ``weight`` (integer) -- the weight
        - ``prec`` (integer or ``None``, default: ``None``) -- power series
          precision. If ``None``, the precision defaults to the Sturm bound for
          the requested level and weight.
        - ``use_random`` (boolean, default: True) -- whether or not to use a
          randomized algorithm when building up the space of forms at the given
          weight from known generators of small weight.

        EXAMPLES::

            sage: m = ModularFormsRing(Gamma0(4))
            sage: m.q_expansion_basis(2,10)
            [1 + 24*q^2 + 24*q^4 + 96*q^6 + 24*q^8 + O(q^10),
            q + 4*q^3 + 6*q^5 + 8*q^7 + 13*q^9 + O(q^10)]
            sage: m.q_expansion_basis(3,10)
            []

            sage: X = ModularFormsRing(SL2Z)
            sage: X.q_expansion_basis(12, 10)
            [1 + 196560*q^2 + 16773120*q^3 + 398034000*q^4 + 4629381120*q^5 + 34417656000*q^6 + 187489935360*q^7 + 814879774800*q^8 + 2975551488000*q^9 + O(q^10),
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + 84480*q^8 - 113643*q^9 + O(q^10)]

        We calculate a basis of a massive modular forms space, in two ways.
        Using this module is about twice as fast as Sage's generic code. ::

            sage: A = ModularFormsRing(11).q_expansion_basis(30, prec=40) # long time (5s)
            sage: B = ModularForms(Gamma0(11), 30).q_echelon_basis(prec=40) # long time (9s)
            sage: A == B # long time
            True

        Check that absurdly small values of ``prec`` don't mess things up::

            sage: ModularFormsRing(11).q_expansion_basis(10, prec=5)
            [1 + O(q^5), q + O(q^5), q^2 + O(q^5), q^3 + O(q^5), q^4 + O(q^5), O(q^5), O(q^5), O(q^5), O(q^5), O(q^5)]
        """
        d = self.modular_forms_of_weight(weight).dimension()
        if d == 0:
            return []

        if prec is None:
            prec=self.modular_forms_of_weight(weight).sturm_bound()

        working_prec = max(prec, self.modular_forms_of_weight(weight).sturm_bound())

        gen_weight = min(6, weight)

        while True:
            verbose("Trying to generate the %s-dimensional space at weight %s using generators of weight up to %s" % (d, weight, gen_weight))
            G = self.generators(maxweight=gen_weight, prec=working_prec)
            V = _span_of_forms_in_weight(G, weight, prec=working_prec, use_random=use_random, stop_dim=d)
            if V.rank() == d and (self.base_ring().is_field() or V.index_in_saturation() == 1):
                break
            else:
                gen_weight += 1
                verbose("Need more generators: trying again with generators of weight up to %s" % gen_weight)

        R = G[0][1].parent()
        return [R(list(x), prec=prec) for x in V.gens()]

    def cuspidal_ideal_generators(self, maxweight=8, prec=None):
        r"""
        Calculate generators for the ideal of cuspidal forms in this ring, as a
        module over the whole ring.

        EXAMPLES::

            sage: ModularFormsRing(Gamma0(3)).cuspidal_ideal_generators(maxweight=12)
            [(6, q - 6*q^2 + 9*q^3 + 4*q^4 + O(q^5), q - 6*q^2 + 9*q^3 + 4*q^4 + 6*q^5 + O(q^6))]
            sage: [k for k,f,F in ModularFormsRing(13, base_ring=ZZ).cuspidal_ideal_generators(maxweight=14)]
            [4, 4, 4, 6, 6, 12]
        """
        working_prec = self.modular_forms_of_weight(maxweight).sturm_bound()

        if self.__cached_cusp_maxweight > -1:
            k = self.__cached_cusp_maxweight + 1
            verbose("Already calculated cusp gens up to weight %s -- using those" % (k-1))

            # we may need to increase the precision of the cached cusp
            # generators
            G =  []
            for j,f,F in self.__cached_cusp_gens:
                if f.prec() >= working_prec:
                    f = F.qexp(working_prec).change_ring(self.base_ring())
                G.append( (j,f,F) )
        else:
            k = 2
            G = []


        while k <= maxweight:
            t = verbose("Looking for cusp generators in weight %s" % k)

            kprec = self.modular_forms_of_weight(k).sturm_bound()

            flist = []

            for (j, f, F) in G:
                for g in self.q_expansion_basis(k - j, prec=kprec):
                    flist.append(g*f)
            A = self.base_ring() ** kprec
            W = A.span([A(f.padded_list(kprec)) for f in flist])

            S = self.modular_forms_of_weight(k).cuspidal_submodule()
            if (W.rank() == S.dimension()
                and (self.base_ring().is_field() or W.index_in_saturation() == 1)):
                    verbose("Nothing new in weight %s" % k, t)
                    k += 1
                    continue

            t = verbose("Known cusp generators span a submodule of dimension %s of space of dimension %s" % (W.rank(), S.dimension()), t)

            B = S.q_integral_basis(prec=working_prec)
            V = A.span([A(f.change_ring(self.base_ring()).padded_list(kprec)) for f in B])
            Q = V/W

            for q in Q.gens():
                try:
                    qc = V.coordinates(Q.lift(q))
                except AttributeError:
                    # work around a silly free module bug
                    qc = V.coordinates(q.lift())
                qcZZ = [ZZ(_) for _ in qc] # lift to ZZ so we can define F
                f = sum([B[i] * qcZZ[i] for i in range(len(B))])
                F = S(f)
                G.append((k, f.change_ring(self.base_ring()), F))

            verbose('added %s new generators' % Q.ngens(), t)
            k += 1

        self.__cached_cusp_maxweight = maxweight
        self.__cached_cusp_gens = G

        if prec is None:
            return G
        elif prec <= working_prec:
            return [ (k, f.truncate_powerseries(prec), F) for k,f,F in G]
        else:
            # user wants increased precision, so we may as well cache that
            Gnew = [ (k, F.qexp(prec).change_ring(self.base_ring()), F) for k,f,F in G]
            self.__cached_cusp_gens = Gnew
            return Gnew

    def cuspidal_submodule_q_expansion_basis(self, weight, prec=None):
        r"""
        Calculate a basis of `q`-expansions for the space of cusp forms of
        weight ``weight`` for this group.

        INPUT:

        - ``weight`` (integer) -- the weight
        - ``prec`` (integer or None) -- precision of `q`-expansions to return

        ALGORITHM: Uses the method :meth:`cuspidal_ideal_generators` to
        calculate generators of the ideal of cusp forms inside this ring. Then
        multiply these up to weight ``weight`` using the generators of the
        whole modular form space returned by :meth:`q_expansion_basis`.

        EXAMPLES::

            sage: R = ModularFormsRing(Gamma0(3))
            sage: R.cuspidal_submodule_q_expansion_basis(20)
            [q - 8532*q^6 - 88442*q^7 + O(q^8), q^2 + 207*q^6 + 24516*q^7 + O(q^8), q^3 + 456*q^6 + O(q^8), q^4 - 135*q^6 - 926*q^7 + O(q^8), q^5 + 18*q^6 + 135*q^7 + O(q^8)]

        We compute a basis of a space of very large weight, quickly (using this
        module) and slowly (using modular symbols), and verify that the answers
        are the same. ::

            sage: A = R.cuspidal_submodule_q_expansion_basis(80, prec=30)  # long time (1s on sage.math, 2013)
            sage: B = R.modular_forms_of_weight(80).cuspidal_submodule().q_expansion_basis(prec=30)  # long time (19s on sage.math, 2013)
            sage: A == B # long time
            True
        """
        d = self.modular_forms_of_weight(weight).cuspidal_submodule().dimension()
        if d == 0:
            return []

        minprec = self.modular_forms_of_weight(weight).sturm_bound()
        if prec is None:
            prec = working_prec = minprec
        else:
            working_prec = max(prec, minprec)

        gen_weight = min(6, weight)

        while True:
            verbose("Trying to generate the %s-dimensional cuspidal submodule at weight %s using generators of weight up to %s" % (d, weight, gen_weight))
            G = self.cuspidal_ideal_generators(maxweight=gen_weight, prec=working_prec)

            flist = []
            for (j, f, F) in G:
                for g in self.q_expansion_basis(weight - j, prec=working_prec):
                    flist.append(g*f)

            A = self.base_ring() ** working_prec
            W = A.span([A(f.padded_list(working_prec)) for f in flist])
            if W.rank() == d and (self.base_ring().is_field() or W.index_in_saturation() == 1):
                break
            else:
                gen_weight += 1
                verbose("Need more generators: trying again with generators of weight up to %s" % gen_weight)

        R = G[0][1].parent()
        return [R(list(x), prec=prec) for x in W.gens()]
