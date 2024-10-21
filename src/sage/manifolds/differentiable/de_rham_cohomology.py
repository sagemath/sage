# sage_setup: distribution = sagemath-symbolics
r"""
De Rham Cohomology

Let `M` and `N` be differentiable manifolds and `\varphi\colon M \to N` be
a differentiable map. Then the associated de Rham complex is given by

.. MATH::

    0 \rightarrow \Omega^0(M,\varphi) \xrightarrow{\mathrm{d}_0}
        \Omega^1(M,\varphi) \xrightarrow{\mathrm{d}_1} \dots
        \xrightarrow{\mathrm{d}_{n-1}} \Omega^n(M,\varphi)
        \xrightarrow{\mathrm{d}_{n}} 0,

where `\Omega^k(M,\varphi)` is the module of differential forms of degree `k`,
and `d_k` is the associated exterior derivative. Then the `k`-*th de Rham
cohomology group* is given by

.. MATH::

    H^k_{\mathrm{dR}}(M, \varphi) =
        \left. \mathrm{ker}(\mathrm{d}_k) \middle/
        \mathrm{im}(\mathrm{d}_{k-1}) \right. ,

and the corresponding ring is obtained by

.. MATH::

    H^*_{\mathrm{dR}}(M, \varphi) = \bigoplus^n_{k=0} H^k_{\mathrm{dR}}(M, \varphi).

The de Rham cohomology ring is implemented via :class:`DeRhamCohomologyRing`.
Its elements, the cohomology classes, are represented by
:class:`DeRhamCohomologyClass`.

AUTHORS:

- Michael Jung (2021) : initial version
"""

#******************************************************************************
#       Copyright (C) 2021 Michael Jung <m.jung@vu.nl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.element import AlgebraElement
from sage.categories.algebras import Algebras
from .characteristic_cohomology_class import (CharacteristicCohomologyClassRing,
                                              CharacteristicCohomologyClassRingElement)

class DeRhamCohomologyClass(AlgebraElement):
    r"""
    Define a cohomology class in the de Rham cohomology ring.

    INPUT:

    - ``parent`` -- de Rham cohomology ring represented by an instance of
      :class:`DeRhamCohomologyRing`
    - ``representative`` -- a closed (mixed) differential form representing the
      cohomology class

    .. NOTE::

        The current implementation only provides basic features. Comparison via
        exact forms are not supported at the time being.

    EXAMPLES::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: C = M.de_rham_complex()
        sage: H = C.cohomology()
        sage: omega = M.diff_form(1, [1,1], name='omega')
        sage: u = H(omega); u
        [omega]

    Cohomology classes can be lifted to the algebra of mixed differential
    forms::

        sage: u.lift()
        Mixed differential form omega on the 2-dimensional differentiable
         manifold M

    However, comparison of two cohomology classes is limited the time being::

        sage: eta = M.diff_form(1, [1,1], name='eta')
        sage: H(eta) == u
        True
        sage: H.one() == u
        Traceback (most recent call last):
        ...
        NotImplementedError: comparison via exact forms is currently not supported
    """
    def __init__(self, parent, representative):
        r"""
        Construct an element of the de Rham cohomology ring.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: omega = M.diff_form(1, [1,1], name='omega', latex_name=r'\omega')
            sage: u = H(omega)
            sage: TestSuite(u).run(skip=['_test_eq', '_test_nonzero_equal'])  # equality not fully supported yet
        """
        super().__init__(parent=parent)
        self._representative = representative

    def _repr_(self):
        r"""
        Return a string representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: H.an_element()  # indirect doctest
            [one]
            sage: H.an_element()._repr_()
            '[one]'
        """
        name = self._representative._name
        if name is None:
            name = 'unnamed form'
        return f"[{name}]"

    def _latex_(self):
        r"""
        Return a LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', latex_name=r'\mathcal{M}')
            sage: X.<x,y> = M.chart()
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: omega = M.diff_form(1, [1,1], name='omega', latex_name=r'\omega')
            sage: u = H(omega)
            sage: latex(u)  # indirect doctest
            \left[\omega\right]
            sage: u._latex_()
            '\\left[\\omega\\right]'
        """
        latex_name = self._representative._latex_name
        if latex_name is None:
            latex_name = r'\mathrm{unnamed form}'
        return rf"\left[{latex_name}\right]"

    def representative(self):
        r"""
        Return a representative of ``self`` in the associated de Rham
        complex.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: omega = M.diff_form(2, name='omega')
            sage: omega[0,1] = x
            sage: omega.display()
            omega = x dx∧dy
            sage: u = H(omega); u
            [omega]
            sage: u.representative()
            Mixed differential form omega on the 2-dimensional differentiable
             manifold M
        """
        return self._representative

    lift = representative

    def _add_(self, other):
        r"""
        Addition of two cohomology classes.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: omega = M.diff_form(1, [1,1], name='omega')
            sage: eta = M.diff_form(1, [1,-1], name='eta')
            sage: H(omega) + H(eta)
            [omega+eta]
        """
        return self.parent()(self.representative() + other.representative())

    def cup(self, other):
        r"""
        Cup product of two cohomology classes.

        INPUT:

        - ``other`` -- another cohomology class in the de Rham cohomology

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: omega = M.diff_form(1, [1,1], name='omega')
            sage: eta = M.diff_form(1, [1,-1], name='eta')
            sage: H(omega).cup(H(eta))
            [omega∧eta]
        """
        return self * other

    def _mul_(self, other):
        r"""
        Cup product of two cohomology classes.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: omega = M.diff_form(1, [1,1], name='omega')
            sage: eta = M.diff_form(1, [1,-1], name='eta')
            sage: H(omega) * H(eta)
            [omega∧eta]
        """
        return self.parent()(self.representative().wedge(other.representative()))

    def _rmul_(self, scalar):
        r"""
        Multiplication with a scalar.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: omega = M.diff_form(1, [1,1], name='omega')
            sage: 1/2*H(omega)
            [1/2∧omega]
        """
        return self.parent(scalar * self.representative())

    def _sub_(self, other):
        r"""
        Subtraction of two cohomology classes.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: omega = M.diff_form(1, [1,1], name='omega')
            sage: eta = M.diff_form(1, [1,-1], name='eta')
            sage: H(omega) - H(eta)
            [omega-eta]
        """
        return self.parent()(self.representative() - other.representative())

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        .. WARNING::

            At current stage, the equality operator only checks whether the
            representatives are equal. No further checks are supported so far.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: omega = M.diff_form(1, [1,1], name='omega')
            sage: eta = M.diff_form(1, [1,1], name='eta')
            sage: H(omega) == H(eta)
            True
            sage: H(omega) == H.one()
            Traceback (most recent call last):
            ...
            NotImplementedError: comparison via exact forms is currently not supported
        """
        if self is other:
            return True
        if isinstance(other, type(self)):
            if self.representative() == other.representative():
                return True
        raise NotImplementedError('comparison via exact forms is currently not supported')

class DeRhamCohomologyRing(Parent, UniqueRepresentation):
    r"""
    The de Rham cohomology ring of a de Rham complex.

    This ring is naturally endowed with a multiplication induced by the wedge
    product, called *cup product*, see :meth:`DeRhamCohomologyClass.cup`.

    .. NOTE::

        The current implementation only provides basic features. Comparison via
        exact forms are not supported at the time being.

    INPUT:

    - ``de_rham_complex`` -- a de Rham complex, typically an instance of
      :class:`~sage.manifolds.differentiable.mixed_form_algebra.MixedFormAlgebra`

    EXAMPLES:

    We define the de Rham cohomology ring on a parallelizable manifold `M`::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: C = M.de_rham_complex()
        sage: H = C.cohomology(); H
        De Rham cohomology ring on the 2-dimensional differentiable manifold M

    Its elements are induced by closed differential forms on `M`::

        sage: beta = M.diff_form(1, [1,0], name='beta')
        sage: beta.display()
        beta = dx
        sage: d1 = C.differential(1)
        sage: d1(beta).display()
        dbeta = 0
        sage: b = H(beta); b
        [beta]

    Cohomology classes can be lifted to the algebra of mixed differential
    forms::

        sage: b.representative()
        Mixed differential form beta on the 2-dimensional differentiable
         manifold M

    The ring admits a zero and unit element::

        sage: H.zero()
        [zero]
        sage: H.one()
        [one]
    """
    def __init__(self, de_rham_complex):
        r"""
        Construct the de Rham cohomology ring.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology(); H
            De Rham cohomology ring on the 2-dimensional differentiable
             manifold M
            sage: TestSuite(H).run(skip=['_test_elements',
            ....:                        '_test_elements_eq_symmetric',
            ....:                       '_test_elements_eq_transitive',
            ....:                       '_test_elements_neq'])  # equality not fully supported yet
        """
        base_field = de_rham_complex.base_ring()
        Parent.__init__(self, base=base_field, category=Algebras(base_field))
        self._de_rham_complex = self._module = de_rham_complex
        self._manifold = de_rham_complex._domain

    Element = DeRhamCohomologyClass

    def _element_constructor_(self, x):
        r"""
        Construct an element of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: H._element_constructor_(C.one())
            [one]

        Non-cycle element::

            sage: omega = M.diff_form(1, name='omega')
            sage: omega[0] = y
            sage: omega.display()
            omega = y dx
            sage: H(omega)
            Traceback (most recent call last):
            ...
            ValueError: Mixed differential form omega on the 2-dimensional
             differentiable manifold M must be a closed form
        """
        if isinstance(x, CharacteristicCohomologyClassRingElement):
            x = x.representative()
        elif x not in self._module:
            raise TypeError(f"{x} must be an element of {self._module}")
        x = self._module(x)
        if x.derivative() != 0:
            raise ValueError(f"{x} must be a closed form")
        return self.element_class(self, x)

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: H.has_coerce_map_from(QQ)
            True

        ::

            sage: M = Manifold(4, 'M')
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: TM = M.tangent_bundle()
            sage: C = TM.characteristic_cohomology_class_ring(); C
            Algebra of characteristic cohomology classes of the Tangent bundle
             TM over the 4-dimensional differentiable manifold M
            sage: H.has_coerce_map_from(C)
            True
        """
        if isinstance(other, CharacteristicCohomologyClassRing):
            # TODO: we need to be careful if manifolds have boundary!
            return other._vbundle._base_space == self._manifold
        return super()._coerce_map_from_(other)

    def _repr_(self):
        r"""
        Return a string representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology(); H
            De Rham cohomology ring on the 2-dimensional differentiable
             manifold M
        """
        desc = "De Rham cohomology ring "
        if self._module._dest_map is self._manifold.identity_map():
            desc += "on the {}".format(self._manifold)
        else:
            desc += "along the {} mapped ".format(self._manifold)
            desc += "into the {} ".format(self._module._ambient_domain)
            if self._module._dest_map._name is None:
                dm_name = "unnamed map"
            else:
                dm_name = self._module._dest_map._name
            desc += "via " + dm_name
        return desc

    def _latex_(self):
        r"""
        Return a LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(3, 'M', latex_name=r'\mathcal{M}')
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: H._latex_()
            'H^*_{\\mathrm{dR}}\\left(\\mathcal{M}\\right)'
            sage: latex(H)  # indirect doctest
            H^*_{\mathrm{dR}}\left(\mathcal{M}\right)
        """
        latex_name = r"H^*_{\mathrm{dR}}\left(" + self._manifold._latex_name
        if self._module._dest_map is not self._manifold.identity_map():
            dm_latex_name = self._module._dest_map._latex_name
            if dm_latex_name is None:
                dm_latex_name = r"\mathrm{unnamed\; map}"
            latex_name += "," + dm_latex_name
        latex_name += r"\right)"
        return latex_name

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: H.an_element()
            [one]
        """
        return self.one()

    @cached_method
    def zero(self):
        r"""
        Return the zero element of ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: H.zero()
            [zero]
            sage: H.zero().representative()
            Mixed differential form zero on the 2-dimensional differentiable
             manifold M
        """
        return self.element_class(self, self._module.zero())

    @cached_method
    def one(self):
        r"""
        Return the one element of ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: C = M.de_rham_complex()
            sage: H = C.cohomology()
            sage: H.one()
            [one]
            sage: H.one().representative()
            Mixed differential form one on the 2-dimensional differentiable
             manifold M
        """
        return self.element_class(self, self._module.one())
