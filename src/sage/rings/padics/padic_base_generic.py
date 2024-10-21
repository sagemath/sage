# sage_setup: distribution = sagemath-pari
r"""
`p`-adic Base Generic

A superclass for implementations of `\ZZ_p` and `\QQ_p`.

AUTHORS:

- David Roe
"""

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from .padic_generic import pAdicGeneric
from .misc import precprint
from sage.rings.padics.pow_computer import PowComputer
from sage.rings.padics.padic_capped_relative_element import pAdicCoercion_ZZ_CR, pAdicCoercion_QQ_CR, pAdicConvert_QQ_CR
from sage.rings.padics.padic_capped_absolute_element import pAdicCoercion_ZZ_CA, pAdicConvert_QQ_CA
from sage.rings.padics.padic_fixed_mod_element import pAdicCoercion_ZZ_FM, pAdicConvert_QQ_FM
from sage.rings.padics.padic_floating_point_element import pAdicCoercion_ZZ_FP, pAdicCoercion_QQ_FP, pAdicConvert_QQ_FP


class pAdicBaseGeneric(pAdicGeneric):
    _implementation = 'GMP'

    def __init__(self, p, prec, print_mode, names, element_class):
        """
        Initialization.

        TESTS::

            sage: R = Zp(5)  # indirect doctest
        """
        if self.is_relaxed():
            from sage.rings.padics.pow_computer_flint import PowComputer_flint
            self.prime_pow = PowComputer_flint(p, 1, 1, 1, self.is_field())
        else:
            self.prime_pow = PowComputer(p, max(min(prec - 1, 30), 1), prec, self.is_field(), self._prec_type())
        pAdicGeneric.__init__(self, self, p, prec, print_mode, names, element_class)
        if self.is_field():
            if self.is_capped_relative():
                coerce_list = [pAdicCoercion_ZZ_CR(self), pAdicCoercion_QQ_CR(self)]
                convert_list = []
            elif self.is_floating_point():
                coerce_list = [pAdicCoercion_ZZ_FP(self), pAdicCoercion_QQ_FP(self)]
                convert_list = []
            elif self.is_lattice_prec():
                coerce_list = [QQ]
                convert_list = []
            elif self.is_relaxed():
                coerce_list = [QQ]
                convert_list = []
            else:
                raise RuntimeError
        elif self.is_capped_relative():
            coerce_list = [pAdicCoercion_ZZ_CR(self)]
            convert_list = [pAdicConvert_QQ_CR(self)]
        elif self.is_capped_absolute():
            coerce_list = [pAdicCoercion_ZZ_CA(self)]
            convert_list = [pAdicConvert_QQ_CA(self)]
        elif self.is_fixed_mod():
            coerce_list = [pAdicCoercion_ZZ_FM(self)]
            convert_list = [pAdicConvert_QQ_FM(self)]
        elif self.is_floating_point():
            coerce_list = [pAdicCoercion_ZZ_FP(self)]
            convert_list = [pAdicConvert_QQ_FP(self)]
        elif self.is_lattice_prec():
            coerce_list = [ZZ]
            convert_list = [QQ]
        elif self.is_relaxed():
            coerce_list = [ZZ]
            convert_list = [QQ]
        else:
            raise RuntimeError
        self.Element = element_class
        self._populate_coercion_lists_(coerce_list=coerce_list, convert_list=convert_list)

    def _repr_(self, do_latex=False):
        r"""
        Return a print representation of this `p`-adic ring or field.

        EXAMPLES::

            sage: K = Zp(17); K  # indirect doctest
            17-adic Ring with capped relative precision 20
            sage: latex(K)
            \Bold{Z}_{17}
            sage: K = ZpCA(17); K  # indirect doctest
            17-adic Ring with capped absolute precision 20
            sage: latex(K)
            \Bold{Z}_{17}
            sage: K = ZpFP(17); K  # indirect doctest
            17-adic Ring with floating precision 20
            sage: latex(K)
            \Bold{Z}_{17}
            sage: K = ZpFM(7); K
            7-adic Ring of fixed modulus 7^20
            sage: latex(K)  # indirect doctest
            \Bold{Z}_{7}
            sage: K = ZpLF(2); K   # indirect doctest
            doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
            See https://github.com/sagemath/sage/issues/23505 for details.
            2-adic Ring with lattice-float precision
            sage: latex(K)
            \Bold{Z}_{2}
            sage: K = Qp(17); K  # indirect doctest
            17-adic Field with capped relative precision 20
            sage: latex(K)
            \Bold{Q}_{17}
            sage: K = QpFP(17); K  # indirect doctest
            17-adic Field with floating precision 20
            sage: latex(K)
            \Bold{Q}_{17}
            sage: K = QpLC(2); K   # indirect doctest
            2-adic Field with lattice-cap precision
            sage: latex(K)
            \Bold{Q}_{2}
        """
        if do_latex:
            if self.is_field():
                s = r"\Bold{Q}_{%s}" % self.prime()
            else:
                s = r"\Bold{Z}_{%s}" % self.prime()
            if hasattr(self, '_label') and self._label:
                s = r"\verb'%s' (\simeq %s)" % (self._label, s)
        else:
            s = "Field " if self.is_field() else "Ring "
            s = "%s-adic " % self.prime() + s + precprint(self._prec_type(), self.precision_cap(), self.prime())
            if hasattr(self, '_label') and self._label:
                s += " (label: %s)" % self._label
        return s

    def exact_field(self):
        """
        Return the rational field.

        For compatibility with extensions of `p`-adics.

        EXAMPLES::

            sage: Zp(5).exact_field()
            Rational Field
        """
        from sage.rings.rational_field import QQ
        return QQ

    def exact_ring(self):
        """
        Return the integer ring.

        EXAMPLES::

            sage: Zp(5).exact_ring()
            Integer Ring
        """
        from sage.rings.integer_ring import ZZ
        return ZZ

    def is_isomorphic(self, ring):
        r"""
        Return whether ``self`` and ``ring`` are isomorphic, i.e. whether
        ``ring`` is an implementation of `\ZZ_p` for the same prime as ``self``.

        INPUT:

        - ``self`` -- a `p`-adic ring

        - ``ring`` -- a ring

        OUTPUT: ``boolean`` -- whether ``ring`` is an implementation of \ZZ_p`
        for the same prime as ``self``

        EXAMPLES::

            sage: R = Zp(5, 15, print_mode='digits'); S = Zp(5, 44, print_max_terms=4); R.is_isomorphic(S)
            True
        """
        return isinstance(ring, pAdicBaseGeneric) and self.prime() == ring.prime() and self.is_field() == ring.is_field()

    def gen(self, n=0):
        """
        Return the `n`-th generator of this extension.  For base
        rings/fields, we consider the generator to be the prime.

        EXAMPLES::

            sage: R = Zp(5); R.gen()
            5 + O(5^21)
        """
        if n != 0:
            raise IndexError("only one generator")
        return self(self.prime())

    def modulus(self, exact=False):
        r"""
        Return the polynomial defining this extension.

        For compatibility with extension fields; we define the modulus to be `x-1`.

        INPUT:

        - ``exact`` -- boolean (default: ``False``); whether to return a
          polynomial with integer entries

        EXAMPLES::

            sage: Zp(5).modulus(exact=True)
            x
        """
        return self.defining_polynomial(exact=exact)

    def absolute_discriminant(self):
        """
        Return the absolute discriminant of this `p`-adic ring.

        EXAMPLES::

            sage: Zp(5).absolute_discriminant()
            1
        """
        return 1

    def discriminant(self, K=None):
        """
        Return the discriminant of this `p`-adic ring over ``K``.

        INPUT:

        - ``self`` -- a `p`-adic ring

        - ``K`` -- a sub-ring of ``self`` or ``None`` (default: ``None``)

        OUTPUT:

        integer; the discriminant of this ring over ``K`` (or the absolute
        discriminant if ``K`` is ``None``)

        EXAMPLES::

            sage: Zp(5).discriminant()
            1
        """
        if (K is None or K is self):
            return 1
        else:
            raise ValueError("Ground Ring must be a subring of self")

    def is_abelian(self):
        """
        Return whether the Galois group is abelian, i.e. ``True``.
        #should this be automorphism group?

        EXAMPLES::

            sage: R = Zp(3, 10,'fixed-mod'); R.is_abelian()
            True
        """
        return True

    def is_normal(self):
        """
        Return whether or not this is a normal extension, i.e. ``True``.

        EXAMPLES::

            sage: R = Zp(3, 10,'fixed-mod'); R.is_normal()
            True
        """
        return True

    def uniformizer(self):
        """
        Return a uniformizer for this ring.

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod', 'series')
            sage: R.uniformizer()
            3
        """
        return self(self.prime_pow._prime())

    def uniformizer_pow(self, n):
        """
        Return the `n`-th power of the uniformizer of ``self`` (as
        an element of ``self``).

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.uniformizer_pow(5)
            5^5 + O(5^25)
            sage: R.uniformizer_pow(infinity)
            0
        """
        return self(self.prime_pow(n))

    def _uniformizer_print(self):
        """
        Return how the uniformizer is supposed to print.

        EXAMPLES::

            sage: R = Zp(5, names='pi'); R._uniformizer_print()
            'pi'
        """
        return self.variable_name()

    def has_pth_root(self):
        r"""
        Return whether or not `\ZZ_p` has a primitive `p`-th
        root of unity.

        EXAMPLES::

            sage: Zp(2).has_pth_root()
            True
            sage: Zp(17).has_pth_root()
            False
        """
        return (self.prime() == 2)

    def has_root_of_unity(self, n):
        r"""
        Return whether or not `\ZZ_p` has a primitive `n`-th
        root of unity.

        INPUT:

        - ``self`` -- a `p`-adic ring

        - ``n`` -- integer

        OUTPUT:

        boolean; whether ``self`` has primitive `n`-th root of unity

        EXAMPLES::

            sage: R=Zp(37)
            sage: R.has_root_of_unity(12)
            True
            sage: R.has_root_of_unity(11)
            False
        """
        if (self.prime() == 2):
            return n.divides(2)
        else:
            return n.divides(self.prime() - 1)

    def zeta(self, n=None):
        r"""
        Return a generator of the group of roots of unity.

        INPUT:

        - ``self`` -- a `p`-adic ring

        - ``n`` -- integer or ``None`` (default: ``None``)

        OUTPUT:

        ``element``; a generator of the `n`-th roots of unity, or a generator
        of the full group of roots of unity if ``n`` is ``None``

        EXAMPLES::

            sage: R = Zp(37,5)
            sage: R.zeta(12)
            8 + 24*37 + 37^2 + 29*37^3 + 23*37^4 + O(37^5)
        """
        if (self.prime() == 2):
            if (n is None) or (n == 2):
                return self(-1)
            if n == 1:
                return self(1)
            else:
                raise ValueError("No, %sth root of unity in self" % n)
        else:
            from sage.rings.finite_rings.finite_field_constructor import GF
            return self.teichmuller(GF(self.prime()).zeta(n).lift())

    def zeta_order(self):
        """
        Return the order of the group of roots of unity.

        EXAMPLES::

            sage: R = Zp(37); R.zeta_order()
            36
            sage: Zp(2).zeta_order()
            2
        """
        if (self.prime() == 2):
            return 2
        else:
            return self.prime() - 1

    def plot(self, max_points=2500, **args):
        r"""
        Create a visualization of this `p`-adic ring as a fractal
        similar to a generalization of the Sierpi\'nski
        triangle.

        The resulting image attempts to capture the
        algebraic and topological characteristics of `\ZZ_p`.

        INPUT:

        - ``max_points`` -- the maximum number or points to plot,
          which controls the depth of recursion (default: 2500)

        - ``**args`` -- color, size, etc. that are passed to the
          underlying point graphics objects

        REFERENCES:

        - Cuoco, A. ''Visualizing the `p`-adic Integers'', The
          American Mathematical Monthly, Vol. 98, No. 4 (Apr., 1991),
          pp. 355-364

        EXAMPLES::

            sage: Zp(3).plot()                                                          # needs sage.plot
            Graphics object consisting of 1 graphics primitive
            sage: Zp(5).plot(max_points=625)                                            # needs sage.plot
            Graphics object consisting of 1 graphics primitive
            sage: Zp(23).plot(rgbcolor=(1,0,0))                                         # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        if 'pointsize' not in args:
            args['pointsize'] = 1
        from sage.misc.mrange import cartesian_product_iterator
        from sage.rings.real_double import RDF
        from sage.plot.all import points
        p = self.prime()
        phi = 2*RDF.pi()/p
        V = RDF**2
        vs = [V([(phi*t).sin(), (phi*t).cos()]) for t in range(p)]
        all = []
        depth = max(RDF(max_points).log(p).floor(), 1)
        scale = min(RDF(1.5/p), 1/RDF(3))
        pts = [vs]*depth
        if depth == 1 and 23 < p < max_points:
            extras = int(max_points/p)
            if p/extras > 5:
                pts = [vs]*depth + [vs[::extras]]
        for digits in cartesian_product_iterator(pts):
            p = sum([v * scale**n for n, v in enumerate(digits)])
            all.append(tuple(p))
        g = points(all, **args)
        # Set default plotting options
        g.axes(False)
        g.set_aspect_ratio(1)
        return g
