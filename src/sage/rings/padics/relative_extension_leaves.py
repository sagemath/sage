# sage_setup: distribution = sagemath-pari
"""
Relative extensions of `p`-adic rings

We represent general extensions of `p`-adic rings as a two-step extension:
first an unramified extension of Qp, followed by an Eisenstein extension
of the result.

This file contains the parent classes for such extensions.
"""

#*****************************************************************************
#       Copyright (C) 2017 David Roe <roed.math@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.morphism import Morphism
from sage.categories.homset import Hom
from .generic_nodes import pAdicFixedModRingGeneric, pAdicCappedAbsoluteRingGeneric, pAdicCappedRelativeRingGeneric, pAdicCappedRelativeFieldGeneric, pAdicFloatingPointRingGeneric, pAdicFloatingPointFieldGeneric
from .eisenstein_extension_generic import EisensteinExtensionGeneric
from .relative_ramified_FM import RelativeRamifiedFixedModElement
from .relative_ramified_CA import RelativeRamifiedCappedAbsoluteElement
from .relative_ramified_CR import RelativeRamifiedCappedRelativeElement
from .relative_ramified_FP import RelativeRamifiedFloatingPointElement
from .pow_computer_relative import PowComputer_relative_maker


class pAdicRelativeBaseringInjection(Morphism):
    """
    The injection of the unramified base into the two-step extension.

    INPUT:

    - ``R`` -- an unramified `p`-adic ring or field
    - ``S`` -- an eisenstein extension of ``R``

    EXAMPLES::

        sage: # needs sage.libs.ntl
        sage: K.<a> = Qq(125)
        sage: R.<x> = K[]
        sage: W.<w> = K.extension(x^3 + 15*a*x - 5*(1+a^2))
        sage: f = W.coerce_map_from(K); f
        Generic morphism:
          From: 5-adic Unramified Extension Field in a defined by x^3 + 3*x + 3
          To:   5-adic Eisenstein Extension Field in w defined by x^3 + 15*a*x - 5*a^2 - 5 over its base field
    """
    def __init__(self, R, S):
        """
        Initialization.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: K.<a> = Qq(125)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^3 + 15*a*x - 5*(1+a^2))
            sage: f = W.coerce_map_from(K)
            sage: type(f)
            <class 'sage.rings.padics.relative_extension_leaves.pAdicRelativeBaseringInjection'>
        """
        if not R.is_field() or S.is_field():
            Morphism.__init__(self, Hom(R, S))
        else:
            from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
            Morphism.__init__(self, Hom(R, S, SetsWithPartialMaps()))

    def _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: K.<a> = Qq(125,2)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^3 + 15*a*x - 5*(1+a^2))
            sage: f = W.coerce_map_from(K)
            sage: f(a+5)  # indirect doctest
            a + (4*a^2 + 4*a + 3)*w^3 + (a + 2)*w^4 + (2*a^2 + 4*a + 2)*w^5 + O(w^6)
        """
        if x.is_zero():
            return self.codomain()(0,x.precision_absolute())
        else:
            return self.codomain()([x])

    def _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both).

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: K.<a> = Qq(125,2)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^3 + 15*a*x - 5*(1+a^2))
            sage: f = W.coerce_map_from(K)
            sage: f(5*a,5)
            (4*a^2 + a + 3)*w^3 + (a^2 + 2*a)*w^4 + O(w^5)
            sage: f(5*a,8,2)  # indirect doctest
            (4*a^2 + a + 3)*w^3 + (a^2 + 2*a)*w^4 + O(w^5)
        """
        return self.codomain()([x], *args, **kwds)

    def section(self):
        """
        Map back to the base ring.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: K.<a> = Qq(125,2)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^3 + 15*a*x - 5*(1+a^2))
            sage: f = W.coerce_map_from(K)
            sage: g = f.section()
            sage: g(a + w - w)
            a + O(5^2)
        """
        return pAdicRelativeBaseringSection(self.codomain(), self.domain())


class pAdicRelativeBaseringSection(Morphism):
    """
    The map from a two-step extension back to its maximal unramified subextension.

    EXAMPLES::

        sage: # needs sage.libs.ntl
        sage: K.<a> = Qq(2^10)
        sage: R.<x> = K[]
        sage: W.<w> = K.extension(x^4 + 2*a*x^2 - 16*x - 6)
        sage: f = K.convert_map_from(W); f
        Generic morphism:
          From: 2-adic Eisenstein Extension Field in w defined by x^4 + 2*a*x^2 - 16*x - 6 over its base field
          To:   2-adic Unramified Extension Field in a defined by x^10 + x^6 + x^5 + x^3 + x^2 + x + 1
    """
    def __init__(self, S, R):
        """
        Initialization.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: K.<a> = Qq(2^10)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^4 + 2*a*x^2 - 16*x - 6*a)
            sage: f = K.convert_map_from(W); type(f)
            <class 'sage.rings.padics.relative_extension_leaves.pAdicRelativeBaseringSection'>
        """
        from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
        Morphism.__init__(self, Hom(S, R, SetsWithPartialMaps()))

    def _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: K.<a> = Qq(2^10)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^4 + 2*a*x^2 - 16*x - 6*a)
            sage: f = K.convert_map_from(W)
            sage: f(a + w - w)  # indirect doctest
            a + O(2^20)
            sage: f(w)
            Traceback (most recent call last):
            ...
            ValueError: Element not contained in base ring
        """
        f = x.polynomial()
        if f.degree() > 0:
            raise ValueError("Element not contained in base ring")
        return f[0]

    def _call_with_args(self, x, args=(), kwds={}):
        """
        Used when specifying absolute or relative precision.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: K.<a> = Qq(2^10)
            sage: R.<x> = K[]
            sage: W.<w> = K.extension(x^4 + 2*a*x^2 - 16*x - 6*a)
            sage: f = K.convert_map_from(W)
            sage: f(a, 5)  # indirect doctest
            a + O(2^5)
        """
        return self.codomain()(self._call_(x), *args, **kwds)


class RelativeRamifiedExtensionRingFixedMod(EisensteinExtensionGeneric, pAdicFixedModRingGeneric):
    """
    Two-step extension ring with fixed-mod precision.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: A.<a> = ZqFM(2^10)
        sage: R.<x> = A[]
        sage: W.<w> = A.extension(x^4 + 2*a*x^2 - 16*x - 6*a); W
        2-adic Eisenstein Extension Ring in w defined by x^4 + 2*a*x^2 - 16*x - 6*a over its base ring
        sage: w^4 + 2*a*w^2 - 16*w - 6*a == 0
        True
    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        """
        Initialization.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: A.<a> = ZqFM(5^4)
            sage: R.<x> = A[]
            sage: W.<w> = A.extension(x^3 - 25*(a+1)*x + 10*(a^2+2))
            sage: TestSuite(W).run(max_samples=16)      # long time
        """
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(prec=unram_prec+1)
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, False, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'fixed-mod')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedFixedModElement)
        from .relative_ramified_FM import pAdicCoercion_ZZ_FM, pAdicConvert_QQ_FM
        self.register_coercion(pAdicCoercion_ZZ_FM(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        self.register_conversion(pAdicConvert_QQ_FM(self))


class RelativeRamifiedExtensionRingCappedAbsolute(EisensteinExtensionGeneric, pAdicCappedAbsoluteRingGeneric):
    """
    Two-step extension ring with capped absolute precision.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: A.<a> = ZqCA(2^10)
        sage: R.<x> = A[]
        sage: W.<w> = A.extension(x^4 + 2*a*x^2 - 16*x - 6*a); W
        2-adic Eisenstein Extension Ring in w defined by x^4 + 2*a*x^2 - 16*x - 6*a over its base ring
        sage: w^4 + 2*a*w^2 - 16*w - 6*a == 0
        True
    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        """
        Initialization.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: A.<a> = ZqCA(5^4)
            sage: R.<x> = A[]
            sage: W.<w> = A.extension(x^3 - 25*(a+1)*x + 10*(a^2+2))
            sage: TestSuite(W).run(max_samples=16)      # long time
        """
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(show_prec=False, type='floating-point')
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, False, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'capped-abs')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedCappedAbsoluteElement)
        from .relative_ramified_CA import pAdicCoercion_ZZ_CA, pAdicConvert_QQ_CA
        self.register_coercion(pAdicCoercion_ZZ_CA(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        self.register_conversion(pAdicConvert_QQ_CA(self))


class RelativeRamifiedExtensionRingCappedRelative(EisensteinExtensionGeneric, pAdicCappedRelativeRingGeneric):
    """
    Two-step extension ring with capped relative precision.

    EXAMPLES::

        sage: # needs sage.libs.ntl
        sage: A.<a> = ZqCR(2^10)
        sage: R.<x> = A[]
        sage: W.<w> = A.extension(x^4 + 2*a*x^2 - 16*x - 6*a); W
        2-adic Eisenstein Extension Ring in w defined by x^4 + 2*a*x^2 - 16*x - 6*a over its base ring
        sage: w^4 + 2*a*w^2 - 16*w - 6*a == 0
        True
    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        """
        Initialization.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: A.<a> = ZqCR(5^4)
            sage: R.<x> = A[]
            sage: W.<w> = A.extension(x^3 - 25*(a+1)*x + 10*(a^2+2))
            sage: TestSuite(W).run(max_samples=16)      # long time
        """
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(show_prec=False, type='floating-point')
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, False, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'capped-rel')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedCappedRelativeElement)
        from .relative_ramified_CR import pAdicCoercion_ZZ_CR, pAdicConvert_QQ_CR
        self.register_coercion(pAdicCoercion_ZZ_CR(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        self.register_conversion(pAdicConvert_QQ_CR(self))


class RelativeRamifiedExtensionFieldCappedRelative(EisensteinExtensionGeneric, pAdicCappedRelativeFieldGeneric):
    """
    Two-step extension field with capped relative precision.

    EXAMPLES::

        sage: # needs sage.libs.ntl
        sage: A.<a> = QqCR(2^10)
        sage: R.<x> = A[]
        sage: W.<w> = A.extension(x^4 + 2*a*x^2 - 16*x - 6*a); W
        2-adic Eisenstein Extension Field in w defined by x^4 + 2*a*x^2 - 16*x - 6*a over its base field
        sage: w^4 + 2*a*w^2 - 16*w - 6*a == 0
        True
    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        """
        Initialization.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: A.<a> = QqCR(5^4)
            sage: R.<x> = A[]
            sage: W.<w> = A.extension(x^3 - 25*(a+1)*x + 10*(a^2+2))
            sage: TestSuite(W).run(max_samples=16)      # long time
        """
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring().change(show_prec=False, type='floating-point')
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, True, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'capped-rel')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedCappedRelativeElement)
        from .relative_ramified_CR import pAdicCoercion_ZZ_CR, pAdicCoercion_QQ_CR
        self.register_coercion(pAdicCoercion_ZZ_CR(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        # We also want to convert down to the ring of integers: this is used in teichmuller expansion
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring().integer_ring(), self))
        self.register_coercion(pAdicCoercion_QQ_CR(self))


class RelativeRamifiedExtensionRingFloatingPoint(EisensteinExtensionGeneric, pAdicFloatingPointRingGeneric):
    """
    Two-step extension ring with floating point precision.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: A.<a> = ZqFP(2^10)
        sage: R.<x> = A[]
        sage: W.<w> = A.extension(x^4 + 2*a*x^2 - 16*x - 6*a); W
        2-adic Eisenstein Extension Ring in w defined by x^4 + 2*a*x^2 - 16*x - 6*a over its base ring
        sage: w^4 + 2*a*w^2 - 16*w - 6*a == 0
        True
    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        """
        Initialization.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: A.<a> = ZqFP(5^4)
            sage: R.<x> = A[]
            sage: W.<w> = A.extension(x^3 - 25*(a+1)*x + 10*(a^2+2))
            sage: TestSuite(W).run(max_samples=16)      # long time
        """
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring()#.change(field=False, show_prec=False)
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, False, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'floating-point')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedFloatingPointElement)
        from .relative_ramified_FP import pAdicCoercion_ZZ_FP, pAdicConvert_QQ_FP
        self.register_coercion(pAdicCoercion_ZZ_FP(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        self.register_conversion(pAdicConvert_QQ_FP(self))


class RelativeRamifiedExtensionFieldFloatingPoint(EisensteinExtensionGeneric, pAdicFloatingPointFieldGeneric):
    """
    Two-step extension field with floating point precision.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: A.<a> = QqFP(2^10)
        sage: R.<x> = A[]
        sage: W.<w> = A.extension(x^4 + 2*a*x^2 - 16*x - 6*a); W
        2-adic Eisenstein Extension Field in w defined by x^4 + 2*a*x^2 - 16*x - 6*a over its base field
        sage: w^4 + 2*a*w^2 - 16*w - 6*a == 0
        True
    """
    def __init__(self, exact_modulus, approx_modulus, prec, print_mode, shift_seed, names, implementation):
        """
        Initialization.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: A.<a> = QqFP(5^4)
            sage: R.<x> = A[]
            sage: W.<w> = A.extension(x^3 - 25*(a+1)*x + 10*(a^2+2))
            sage: TestSuite(W).run(max_samples=16)      # long time
        """
        self._exact_modulus = exact_modulus
        unram_prec = (prec + approx_modulus.degree() - 1) // approx_modulus.degree()
        KFP = approx_modulus.base_ring()#.change(field=False, show_prec=False)
        self.prime_pow = PowComputer_relative_maker(approx_modulus.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, True, exact_modulus.change_ring(KFP), shift_seed.change_ring(KFP), 'floating-point')
        self._implementation = 'Polynomial'
        EisensteinExtensionGeneric.__init__(self, approx_modulus, prec, print_mode, names, RelativeRamifiedFloatingPointElement)
        from .relative_ramified_FP import pAdicCoercion_ZZ_FP, pAdicCoercion_QQ_FP
        self.register_coercion(pAdicCoercion_ZZ_FP(self))
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring(), self))
        # We also want to convert down to the ring of integers: this is used in teichmuller expansion
        self.register_coercion(pAdicRelativeBaseringInjection(approx_modulus.base_ring().integer_ring(), self))
        self.register_coercion(pAdicCoercion_QQ_FP(self))
