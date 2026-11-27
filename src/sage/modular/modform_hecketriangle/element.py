# sage.doctest: needs sage.combinat sage.graphs
r"""
Elements of Hecke modular forms spaces

AUTHORS:

- Jonas Jermann (2013): initial version
"""

# ****************************************************************************
#       Copyright (C) 2013-2014 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .graded_ring_element import FormsRingElement


class FormsElement(FormsRingElement):
    """
    (Hecke) modular forms.
    """

    def __init__(self, parent, rat) -> None:
        r"""
        An element of a space of (Hecke) modular forms.

        INPUT:

        - ``parent`` -- a modular form space

        - ``rat`` -- a rational function which corresponds to a
          modular form in the modular form space

        OUTPUT:

        A (Hecke) modular form element corresponding to the given
        rational function with the given parent space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: (x,y,z,d)=var("x,y,z,d")
            sage: MF = ModularForms(n=5, k=20/3, ep=1)
            sage: MF.default_prec(3)
            sage: el = MF(x^5*d-y^2*d)
            sage: el
            q - 9/(200*d)*q^2 + O(q^3)
            sage: el.rat()
            x^5*d - y^2*d
            sage: el.parent()
            ModularForms(n=5, k=20/3, ep=1) over Integer Ring
            sage: el.rat().parent()
            Fraction Field of Multivariate Polynomial Ring in x, y, z, d over Integer Ring

            sage: subspace = MF.subspace([MF.gen(1)])
            sage: ss_el = subspace(x^5*d-y^2*d)
            sage: ss_el == el
            True
            sage: ss_el.parent()
            Subspace of dimension 1 of ModularForms(n=5, k=20/3, ep=1) over Integer Ring
        """
        super().__init__(parent, rat)

        if self.AT(["quasi"]) >= self._analytic_type:
            pass
        elif not (self.is_homogeneous() and
                  self._weight == parent.weight() and
                  self._ep == parent.ep()):
            raise ValueError("{} does not correspond to an element of {}.".format(rat, parent))

        from .subspace import SubSpaceForms
        if isinstance(parent, SubSpaceForms) and (parent._module is not None):
            try:
                self.coordinate_vector()
            except TypeError:
                raise ValueError(f"{rat} does not correspond to an element of {parent}")

    def _repr_(self) -> str:
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: (x,y,z,d)=var("x,y,z,d")
            sage: QuasiModularForms(n=5, k=10, ep=-1)(x^3*z^3-y^3)
            21/(20*d)*q - 4977/(16000*d^2)*q^2 + 297829/(12800000*d^3)*q^3 + 27209679/(20480000000*d^4)*q^4 + O(q^5)
            sage: QuasiModularForms(n=infinity, k=8, ep=1)(x*(x-y^2))
            64*q + 512*q^2 + 768*q^3 - 4096*q^4 + O(q^5)
        """

        return self._qexp_repr()

    def _latex_(self) -> str:
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: (x,y,z,d)=var("x,y,z,d")
            sage: latex(QuasiModularForms(n=5, k=10, ep=-1)(x^3*z^3-y^3))
            f_{\rho}^{3} E_{2}^{3} -  f_{i}^{3}
            sage: latex(QuasiModularForms(n=infinity, k=8, ep=1)(x*(x-y^2)))
            -E_{4} f_{i}^{2} + E_{4}^{2}
        """
        return super()._latex_()

    def coordinate_vector(self):
        r"""
        Return the coordinate vector of ``self`` with
        respect to ``self.parent().gens()``.

        .. NOTE::

            This uses the corresponding function of the
            parent. If the parent has not defined a coordinate
            vector function or a module for coordinate vectors
            then an exception is raised by the parent
            (default implementation).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=4, k=24, ep=-1)
            sage: MF.gen(0).coordinate_vector().parent()
            Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.gen(0).coordinate_vector()
            (1, 0, 0)
            sage: subspace = MF.subspace([MF.gen(0), MF.gen(2)])
            sage: subspace.gen(0).coordinate_vector().parent()
            Vector space of dimension 2 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: subspace.gen(0).coordinate_vector()
            (1, 0)
            sage: subspace.gen(0).coordinate_vector() == subspace.coordinate_vector(subspace.gen(0))
            True
        """
        return self.parent().coordinate_vector(self)

    def ambient_coordinate_vector(self):
        r"""
        Return the coordinate vector of ``self`` with
        respect to ``self.parent().ambient_space().gens()``.

        The returned coordinate vector is an element
        of ``self.parent().module()``.

        .. NOTE::

            This uses the corresponding function of the
            parent. If the parent has not defined a coordinate
            vector function or an ambient module for
            coordinate vectors then an exception is raised
            by the parent (default implementation).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=4, k=24, ep=-1)
            sage: MF.gen(0).ambient_coordinate_vector().parent()
            Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.gen(0).ambient_coordinate_vector()
            (1, 0, 0)
            sage: subspace = MF.subspace([MF.gen(0), MF.gen(2)])
            sage: subspace.gen(0).ambient_coordinate_vector().parent()
            Vector space of degree 3 and dimension 2 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            Basis matrix:
            [1 0 0]
            [0 0 1]
            sage: subspace.gen(0).ambient_coordinate_vector()
            (1, 0, 0)
            sage: subspace.gen(0).ambient_coordinate_vector() == subspace.ambient_coordinate_vector(subspace.gen(0))
            True
        """
        return self.parent().ambient_coordinate_vector(self)

    def lseries(self, num_prec=53, max_imaginary_part=0):
        r"""
        Return the `L`-series of ``self`` if ``self`` is modular and holomorphic.

        This relies on the implementation in PARI.

        INPUT:

        - ``num_prec`` -- integer (default: 53); giving the required precision

        - ``max_imaginary_part`` -- a real number (default: 0), indicating up
          to which imaginary part the `L`-series is going to be studied

        - ``max_asymp_coeffs`` -- ignored

        OUTPUT:

        An interface to PARI functions for computing `L`-series, namely
        the series given by the Fourier coefficients of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform.eis_series import eisenstein_series_lseries
            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: f = ModularForms(n=3, k=4).E4()/240
            sage: L = f.lseries(); L
            L-series associated to the modular form 1/240 + q + 9*q^2 + 28*q^3 + 73*q^4 + O(q^5)
            sage: L.conductor
            1
            sage: L(1).prec()
            53
            sage: L.check_functional_equation() < 2^(-50)
            True
            sage: L(1)
            -0.0304484570583...
            sage: abs(L(1) - eisenstein_series_lseries(4)(1)) < 2^(-53)
            True
            sage: L.derivative(1, 1)
            -0.0504570844798...
            sage: L.derivative(1, 2)/2
            -0.0350657360354...
            sage: L.taylor_series(1, 3)
            -0.0304484570583... - 0.0504570844798...*z - 0.0350657360354...*z^2 + O(z^3)
            sage: coeffs = f.q_expansion_vector(min_exp=0, max_exp=20, fix_d=True)
            sage: sum([coeffs[k] * ZZ(k)^(-10) for k in range(1,len(coeffs))]).n(53)
            1.00935215408...
            sage: L(10)
            1.00935215649...

            sage: f = ModularForms(n=6, k=4).E4()
            sage: L = f.lseries(num_prec=200)
            sage: L.conductor
            3
            sage: L.check_functional_equation() < 2^(-180)
            True
            sage: L(1)
            -2.92305187760575399490414692523085855811204642031749788...
            sage: L(1).prec()
            200
            sage: coeffs = f.q_expansion_vector(min_exp=0, max_exp=20, fix_d=True)
            sage: sum([coeffs[k] * ZZ(k)^(-10) for k in range(1,len(coeffs))]).n(53)
            24.2281438789...
            sage: L(10).n(53)
            24.2281439447...

            sage: f = ModularForms(n=8, k=6, ep=-1).E6()
            sage: L = f.lseries()
            sage: L.check_functional_equation() < 2^(-45)
            True
            sage: L.taylor_series(3, 3)
            0.000000000000... + 0.867197036668...*z + 0.261129628199...*z^2 + O(z^3)
            sage: coeffs = f.q_expansion_vector(min_exp=0, max_exp=20, fix_d=True)
            sage: sum([coeffs[k]*k^(-10) for k in range(1,len(coeffs))]).n(53)
            -13.0290002560...
            sage: L(10).n(53)
            -13.0290184579...

            sage: # long time
            sage: f = ModularForms(n=17, k=24).Delta()^2
            sage: L = f.lseries()
            sage: L.check_functional_equation() < 2^(-50)
            True
            sage: L.taylor_series(12, 3)
            0.000683924755280... - 0.000875942285963...*z + 0.000647618966023...*z^2 + O(z^3)
            sage: coeffs = f.q_expansion_vector(min_exp=0, max_exp=20, fix_d=True)
            sage: sum([coeffs[k]*k^(-30) for k in range(1,len(coeffs))]).n(53)
            9.31562890589...e-10
            sage: L(30).n(53)
            9.31562890589...e-10

            sage: f = ModularForms(n=infinity, k=2, ep=-1).f_i()
            sage: L = f.lseries()
            sage: L.check_functional_equation() < 2^(-50)
            True
            sage: L.taylor_series(1, 3)
            0.000000000000... + 5.76543616701...*z + 9.92776715593...*z^2 + O(z^3)
            sage: coeffs = f.q_expansion_vector(min_exp=0, max_exp=20, fix_d=True)
            sage: sum([coeffs[k] * ZZ(k)^(-10) for k in range(1,len(coeffs))]).n(53)
            -23.9781792831...
            sage: L(10).n(53)
            -23.9781792831...
        """
        from sage.rings.integer import Integer
        from sage.lfunctions.pari import lfun_generic, LFunction
        from sage.libs.pari import pari

        if (not (self.is_modular() and self.is_holomorphic())
                or self.weight() == 0):
            raise NotImplementedError("L-series are only implemented for "
                                      "non-trivial holomorphic modular forms")

        conductor = self.group().lam()**2
        if self.group().is_arithmetic():
            conductor = Integer(conductor)
        else:
            minpoly = conductor.minpoly()._pari_init_()
            # conductor is 4 cos(2 i pi / 2n)^2
            # how to pick up the correct real root ?
            pari_conductor = pari(f"()->polrootsreal({minpoly}, [0,4])[1]")
            # conductor in a number field with real embedding or in AA
            # TODO : pass this to PARI as a zero-arity enclosure ()->
            conductor = pari_conductor  # conductor.n(num_prec)

        if self.is_cuspidal():
            poles = []
            residues = []
        else:
            poles = 0   # automatic, instead of [weight]
            residues = [None]

        Lpari = lfun_generic(conductor=conductor,
                             gammaV=[0, 1],
                             weight=self.weight(),
                             eps=self.ep(),
                             poles=poles,
                             residues=residues)
        L = LFunction(Lpari, prec=num_prec, max_im=max_imaginary_part)
        nterms = Integer(L.cost())
        an_list = list(self.q_expansion_vector(min_exp=0,
                                               max_exp=nterms + 1,
                                               fix_d=True))
        Lpari.init_coeffs(an_list[1:])
        L = LFunction(Lpari, prec=num_prec, max_im=max_imaginary_part)
        is_good = L.check_functional_equation()
        assert is_good < 1e-10, is_good

        mf = "cusp" if self.is_cuspidal() else "modular"
        L.rename(f"L-series associated to the {mf} form {self}")
        return L
