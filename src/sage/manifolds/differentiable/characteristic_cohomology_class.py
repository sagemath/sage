# sage_setup: distribution = sagemath-symbolics
r"""
Characteristic cohomology classes

A *characteristic class* `\kappa` is a natural transformation that
associates with each vector bundle `E \to M` a cohomology class
`\kappa(E) \in H^*(M;R)` such that for any continuous map `f\colon N \to M`
from another topological manifold `N`, the *naturality condition* is
satisfied:

.. MATH::

    f^*\kappa(E) = \kappa(f^* E) \in H^*(N;R)

The cohomology class `\kappa(E)` is called *characteristic cohomology class*.
Roughly speaking, characteristic cohomology classes measure the non-triviality
of vector bundles.

One way to obtain and compute characteristic classes in the de Rham cohomology
with coefficients in the ring `\CC` is via the so-called *Chern-Weil theory*
using the curvature of a differentiable vector bundle.

For that let `\nabla` be a connection on `E`, `e` a local frame on
`E` and `\Omega` be the corresponding curvature matrix
(see: :meth:`~sage.manifolds.differentiable.bundle_connection.BundleConnection.curvature_form`).

Namely, if `P: \mathrm{Mat}_{n \times n}(\CC) \to \CC` is an invariant
polynomial, the object

.. MATH::

    \left[ P \left( \Omega \right) \right] \in H^{2*}_{\mathrm{dR}}(M, \CC)

is well-defined, independent of the choice of `\nabla` (the proof can be
found in [Roe1988]_ pp. 31) and fulfills the naturality condition.
This is the foundation of the Chern-Weil theory and the following definitions.

.. NOTE::

    This documentation is rich of examples, but sparse in explanations. Please
    consult the references for more details.

AUTHORS:

- Michael Jung (2019) : initial version
- Michael Jung (2021) : new algorithm; complete refactoring

REFERENCES:

- [Mil1974]_
- [Roe1988]_

Contents
--------

We consider the following three types of classes:

- :ref:`additive`
- :ref:`multiplicative`
- :ref:`Pfaffian`

.. _additive:

Additive Classes
----------------

In the **complex** case, let `f` be a holomorphic function around zero. Then
we call

.. MATH::

    \left[\mathrm{tr}\left( f\left( \frac{\Omega}{2 \pi i} \right)
        \right)\right] \in H^{2*}_{\mathrm{dR}}(M, \CC)

the *additive characteristic class associated to* `f` of the complex vector
bundle `E`.

Important and predefined additive classes are:

- *Chern Character* with `f(x) = \exp(x)`

In the **real** case, let `g` be a holomorphic function around zero with
`g(0)=0`. Then we call

.. MATH::

    \left[\mathrm{tr}\left( \frac{1}{2} g\left( -\frac{\Omega^2}{4 \pi^2}
        \right) \right)\right] \in H^{4*}_{\mathrm{dR}}(M, \CC)

the *additive characteristic class associated to* `g` of the **real** vector
bundle `E`.

EXAMPLES:

Consider the **Chern character** on some 2-dimensional spacetime::

    sage: M = Manifold(2, 'M', structure='Lorentzian')
    sage: X.<t,x> = M.chart()
    sage: E = M.vector_bundle(1, 'E', field='complex'); E
    Differentiable complex vector bundle E -> M of rank 1 over the base space
     2-dimensional Lorentzian manifold M
    sage: e = E.local_frame('e')

Let us define the connection `\nabla^E` in terms of an electro-magnetic
potential `A(t)`::

    sage: nab = E.bundle_connection('nabla^E', latex_name=r'\nabla^E')
    sage: omega = M.one_form(name='omega')
    sage: A = function('A')
    sage: nab.set_connection_form(0, 0)[1] = I*A(t)
    sage: nab[0, 0].display()
    connection (0,0) of bundle connection nabla^E w.r.t. Local frame
     (E|_M, (e_0)) = I*A(t) dx
    sage: nab.set_immutable()

The Chern character is then given by::

    sage: ch = E.characteristic_cohomology_class('ChernChar'); ch
    Characteristic cohomology class ch(E) of the Differentiable complex vector
     bundle E -> M of rank 1 over the base space 2-dimensional Lorentzian
     manifold M

The corresponding characteristic form w.r.t. the bundle connection can be
obtained via :meth:`get_form`::

    sage: ch_form = ch.get_form(nab); ch_form.display_expansion()
    ch(E, nabla^E) = 1 + 1/2*d(A)/dt/pi dt∧dx

.. _multiplicative:

Multiplicative Classes
----------------------

In the **complex** case, let `f` be a holomorphic function around zero.
Then we call

.. MATH::

    \left[\det\left( f\left( \frac{\Omega}{2 \pi i} \right)
        \right)\right] \in H^{2*}_{\mathrm{dR}}(M, \CC)

the *multiplicative characteristic class associated to* `f` of the complex
vector bundle `E`.

Important and predefined multiplicative classes on complex vector bundles are:

- *Chern class* with `f(x) = 1+x`
- *Todd class* with `f(x) = \frac{x}{1-\exp(-x)}`

In the **real** case, let `g` be a holomorphic function around zero with
`g(0)=1`. Then we call

.. MATH::

    \left[\det\left( \sqrt{ g \left( -\frac{\Omega^2}{4 \pi^2} \right) } \right)
        \right] \in H^{4*}_{\mathrm{dR}}(M, \CC)

the *multiplicative characteristic class associated to* `g` on the **real**
vector bundle `E`.

Important and predefined multiplicative classes on real vector bundles are:

- *Pontryagin class* with `g(x) = 1+x`
- `\hat{A}` *class* with `g(x) = \frac{\sqrt{x}/2}{\sinh(\sqrt{x}/2)}`
- *Hirzebruch class* with `g(x) = \frac{\sqrt{x}}{\tanh(\sqrt{x})}`

EXAMPLES:

We consider the **Chern class** of the tautological line bundle `\gamma^1` over
`\CC\mathbf{P}^1`::

    sage: M = Manifold(2, 'CP^1', start_index=1)
    sage: U = M.open_subset('U')
    sage: c_cart.<x,y> = U.chart() # homogeneous coordinates in real terms
    sage: c_comp.<z, zbar> = U.chart(r'z:z zbar:\bar{z}') # complexification
    sage: cart_to_comp = c_cart.transition_map(c_comp, (x+I*y, x-I*y))
    sage: comp_to_cart = cart_to_comp.inverse()
    sage: E = M.vector_bundle(1, 'gamma^1', field='complex')
    sage: e = E.local_frame('e', domain=U)

To apply the Chern-Weil approach, we need a bundle connection in terms of a
connection one form. To achieve this, we take the connection induced from the
hermitian metric on the trivial bundle
`\CC^2 \times \CC\mathbf{P}^1 \supset \gamma^1`. In this the frame `e`
corresponds to the section `[z:1] \mapsto (z,1)` and its magnitude-squared
is given by `1+|z|^2`::

    sage: nab = E.bundle_connection('nabla')
    sage: omega = U.one_form(name='omega')
    sage: omega[c_comp.frame(),1,c_comp] = zbar/(1+z*zbar)
    sage: nab[e, 1, 1] = omega
    sage: nab.set_immutable()

Now, the Chern class can be constructed::

    sage: c = E.characteristic_cohomology_class('Chern'); c
    Characteristic cohomology class c(gamma^1) of the Differentiable complex
     vector bundle gamma^1 -> CP^1 of rank 1 over the base space 2-dimensional
     differentiable manifold CP^1
    sage: c_form = c.get_form(nab)
    sage: c_form.display_expansion(c_comp.frame(), chart=c_comp)
    c(gamma^1, nabla) = 1 + 1/2*I/(pi + pi*z^2*zbar^2 + 2*pi*z*zbar) dz∧dzbar

Since `U` and `\CC\mathbf{P}^1` differ only by a point and therefore a null
set, it is enough to integrate the top form over the domain `U`::

    sage: integrate(integrate(c_form[2][[1,2]].expr(c_cart), x, -infinity, infinity).full_simplify(),
    ....:           y, -infinity, infinity)
    1

The result shows that `c_1(\gamma^1)` generates the second integer
cohomology of `\CC\mathbf{P}^1`.

.. _Pfaffian:

Pfaffian Classes
----------------

Usually, there is no such thing as "Pfaffian classes" in literature. However,
using the matrix' Pfaffian and inspired by the aforementioned definitions,
such classes can be defined as follows.

Let `E` be a real vector bundle of rank `2n` and `f` an odd real function
being analytic at zero. Furthermore, let `\Omega` be skew-symmetric, which
certainly will be true if `\nabla` is metric and `e` is orthonormal. Then
we call

.. MATH::

    \left[\mathrm{Pf}\left( f\left( \frac{\Omega}{2 \pi} \right) \right)\right]
        \in H^{2n*}(M,\RR)

the *Pfaffian class associated to f*.

The most important Pfaffian class is the *Euler class* which is simply given by
`f(x)=x`.

EXAMPLES:

We consider the **Euler class** of `S^2`::

    sage: M.<x,y> = manifolds.Sphere(2, coordinates='stereographic')
    sage: TM = M.tangent_bundle()
    sage: e_class = TM.characteristic_cohomology_class('Euler'); e_class
    Characteristic cohomology class e(TS^2) of the Tangent bundle TS^2 over the
     2-sphere S^2 of radius 1 smoothly embedded in the Euclidean space E^3

To compute a particular representative of the Euler class, we need to determine
a connection, which is in this case given by the standard metric::

    sage: g = M.metric('g') # standard metric on S2, long time
    sage: nab = g.connection() # long time
    sage: nab.set_immutable() # long time

Now the representative of the Euler class with respect to the connection
`\nabla_g` induced by the standard metric can be computed::

    sage: e_class_form = e_class.get_form(nab) # long time
    sage: e_class_form.display_expansion() # long time
    e(TS^2, nabla_g) = 2/(pi + pi*x^4 + pi*y^4 + 2*pi*x^2 + 2*(pi + pi*x^2)*y^2) dx∧dy

Let us check whether this form represents the Euler class correctly::

    sage: # long time
    sage: expr = e_class_form[2][[1,2]].expr()
    sage: expr = integrate(expr, x, -infinity, infinity)
    sage: expr = expr.simplify_full()
    sage: integrate(expr, y, -infinity, infinity)
    2

As we can see, the integral coincides with the Euler characteristic of `S^2` so
that our form actually represents the Euler class appropriately.
"""

# *****************************************************************************
#       Copyright (C) 2021 Michael Jung <m.jung at vu.nl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.algebras.finite_gca import FiniteGCAlgebra
from sage.combinat.free_module import IndexedFreeModuleElement
from sage.misc.fast_methods import Singleton
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from .affine_connection import AffineConnection
from .bundle_connection import BundleConnection
from .levi_civita_connection import LeviCivitaConnection
from sage.symbolic.expression import Expression
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class CharacteristicCohomologyClassRingElement(IndexedFreeModuleElement):
    r"""
    Characteristic cohomology class.

    Let `E \to M` be a real/complex vector bundle of rank `k`. A characteristic
    cohomology class of `E` is generated by either

    - Chern classes if `E` is complex,
    - Pontryagin classes if `E` is real,
    - Pontryagin classes and the Euler class if `E` is real and orientable,

    via the Chern-Weil homomorphism.

    For details about the ring structure, see
    :class:`CharacteristicCohomologyClassRing`.

    To construct a characteristic cohomology class, please use
    :func:`CharacteristicCohomologyClass`.

    EXAMPLES::

        sage: M = Manifold(12, 'M')
        sage: TM = M.tangent_bundle()
        sage: CR = TM.characteristic_cohomology_class_ring()
        sage: p1, p2, p3 = CR.gens()
        sage: p1*p2+p3
        Characteristic cohomology class (p_1⌣p_2 + p_3)(TM) of the Tangent
         bundle TM over the 12-dimensional differentiable manifold M
        sage: A = TM.characteristic_cohomology_class('AHat'); A
        Characteristic cohomology class A^(TM) of the Tangent bundle TM over
         the 12-dimensional differentiable manifold M
        sage: A == 1 - p1/24 + (7*p1^2-4*p2)/5760 + (44*p1*p2-31*p1^3-16*p3)/967680
        True
    """
    def __init__(self, parent, x, name=None, latex_name=None):
        r"""
        Construct a characteristic cohomology class.

        TESTS::

            sage: M = Manifold(8, 'M')
            sage: TM = M.tangent_bundle()
            sage: p = TM.characteristic_cohomology_class('Pontryagin')
            sage: TestSuite(p).run()
        """
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._mixed_forms = {}  # dict. of characteristic forms of self
                                # (key: bundle connection)
        super().__init__(parent, x)

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(8, 'M')
            sage: TM = M.tangent_bundle()
            sage: p = TM.characteristic_cohomology_class('Pontryagin')
            sage: p._repr_()
            'Characteristic cohomology class p(TM) of the Tangent bundle TM
             over the 8-dimensional differentiable manifold M'
            sage: p  # indirect doctest
            Characteristic cohomology class p(TM) of the Tangent bundle TM over
             the 8-dimensional differentiable manifold M

        ::

            sage: x = var('x')
            sage: k = TM.characteristic_cohomology_class(1+x^2, class_type='multiplicative')
            sage: k._repr_()
            'Characteristic cohomology class (1 + p_1^2 - 2*p_2)(TM) of the
             Tangent bundle TM over the 8-dimensional differentiable manifold M'
            sage: k  # indirect doctest
            Characteristic cohomology class (1 + p_1^2 - 2*p_2)(TM) of the
             Tangent bundle TM over the 8-dimensional differentiable manifold M
        """
        if self._name is None:
            name = f'({super()._repr_()})'
        else:
            name = self._name
        vbundle = self.parent()._vbundle
        name = f'{name}({vbundle._name})'
        return f'Characteristic cohomology class {name} of the {vbundle}'

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(8, 'M')
            sage: TM = M.tangent_bundle()
            sage: p = TM.characteristic_cohomology_class('Pontryagin')
            sage: p._latex_()
            'p\\left(TM\\right)'
            sage: latex(p)  # indirect doctest
            p\left(TM\right)

        ::

            sage: x = var('x')
            sage: k = TM.characteristic_cohomology_class(1+x^2, class_type='multiplicative')
            sage: k._latex_()
            '\\left(1 + p_1^{2} - 2 p_2\\right)\\left(TM\\right)'
            sage: latex(k)
            \left(1 + p_1^{2} - 2 p_2\right)\left(TM\right)
        """
        if self._latex_name is None:
            latex = r'\left(' + super()._latex_() + r'\right)'
        else:
            latex = self._latex_name
        vbundle = self.parent()._vbundle
        latex += r'\left(' + vbundle._latex_name + r'\right)'
        return latex

    def get_form(self, nab):
        r"""
        Return the characteristic form of ``self``.

        INPUT:

        - ``nab`` -- get the characteristic form w.r.t. to the
          connection ``nab``

        OUTPUT:

        - an instance of
          :class:`sage.manifolds.differentiable.mixed_form.MixedForm`

        EXAMPLES:

        Trivial characteristic form on Euclidean space::

            sage: M = manifolds.EuclideanSpace(4)
            sage: TM = M.tangent_bundle()
            sage: one = TM.characteristic_cohomology_class_ring().one()
            sage: g = M.metric()
            sage: nab = g.connection()
            sage: nab.set_immutable()
            sage: one.get_form(nab)
            Mixed differential form one on the 4-dimensional Euclidean space E^4

        Pontryagin form on the 4-sphere::

            sage: M = manifolds.Sphere(4)
            sage: TM = M.tangent_bundle()
            sage: p = TM.characteristic_cohomology_class('Pontryagin'); p
            Characteristic cohomology class p(TS^4) of the Tangent bundle TS^4
             over the 4-sphere S^4 of radius 1 smoothly embedded in the
             5-dimensional Euclidean space E^5
            sage: g = M.metric() # long time
            sage: nab = g.connection() # long time
            sage: nab.set_immutable() # long time
            sage: p_form = p.get_form(nab); p_form # long time
            Mixed differential form p(TS^4, nabla_g) on the 4-sphere S^4 of
             radius 1 smoothly embedded in the 5-dimensional Euclidean space E^5
            sage: p_form.display_expansion() # long time
            p(TS^4, nabla_g) = 1
        """
        if nab not in self._mixed_forms:
            dom = nab._domain
            A = dom.mixed_form_algebra()

            # trivial cases
            if self == 1:
                self._mixed_forms[nab] = A(dom._one_scalar_field)
            elif self == 0:
                self._mixed_forms[nab] = A(dom._zero_scalar_field)
            else:  # non-trivial case
                from functools import reduce

                parent = self.parent()
                algorithm = parent._algorithm

                grading = parent.print_options()['sorting_key']
                res = [dom.diff_form_module(i).zero()
                       for i in range(dom._dim + 1)]
                for ind, c in self:
                    deg = grading(ind)
                    gen_pow = [algorithm.get_gen_pow(nab, i, ind[i])
                               for i in range(len(ind))]
                    res[deg] += c * reduce(lambda x, y: x.wedge(y), gen_pow)

                res = A(res)  # convert result into mixed form

                # preparse names
                vbundle = parent._vbundle
                if self._name is None:
                    name = f'({super()._repr_()})'
                else:
                    name = self._name
                if self._latex_name is None:
                    latex_name = r'\left(' + super()._latex_() + r'\right)'
                else:
                    latex_name = self._latex_name
                # appendix
                append_name = f'({vbundle._name}, {nab._name})'
                append_latex_name = r'\left(' + vbundle._latex_name
                append_latex_name += ', ' + nab._latex_name + r'\right)'

                # set names of components
                from sage.arith.misc import gcd

                step = gcd(parent._degrees)  # step size of (possibly) nonzero
                for i in range(dom._dim // step + 1):
                    # enumerate (possibly) nonzero components
                    comp_name = name + f'_{i}' + append_name
                    comp_latex_name = latex_name + r'_{' + str(i) + '}'
                    comp_latex_name += append_latex_name
                    res[step * i].set_name(name=comp_name,
                                           latex_name=comp_latex_name)

                # set global names
                res._name = name + append_name
                res._latex_name = latex_name + append_latex_name

                res.set_immutable()

                self._mixed_forms[nab] = res  # cache result in dict

        return self._mixed_forms[nab]

    def representative(self, nab=None):
        r"""
        Return any representative of ``self``.

        INPUT:

        - ``nab`` -- (default: ``None``) if stated, return the representative
          w.r.t. to the connection ``nab``; otherwise an arbitrary already
          computed representative will be chosen.

        OUTPUT:

        - an instance of
          :class:`sage.manifolds.differentiable.mixed_form.MixedForm`

        EXAMPLES:

        Define the 4-dimensional Euclidean space::

            sage: M = manifolds.EuclideanSpace(4)
            sage: TM = M.tangent_bundle()
            sage: one = TM.characteristic_cohomology_class_ring().one()

        No characteristic form has been computed so far, thus we get an error::

            sage: one.representative()
            Traceback (most recent call last):
            ...
            AttributeError: cannot pick a representative

        Get a characteristic form::

            sage: g = M.metric()
            sage: nab = g.connection()
            sage: nab.set_immutable()
            sage: one.get_form(nab)
            Mixed differential form one on the 4-dimensional Euclidean space E^4

        Now, the result is cached and `representative` returns a form::

            sage: one.representative()
            Mixed differential form one on the 4-dimensional Euclidean space E^4

        Alternatively, the option ``nab`` can be used to return the
        characteristic form w.r.t. a fixed connection::

            sage: one.representative(nab)
            Mixed differential form one on the 4-dimensional Euclidean space E^4

        .. SEEALSO::

            :meth:`CharacteristicCohomologyClassRingElement.get_form`
        """
        if nab is None:
            if not self._mixed_forms:
                raise AttributeError('cannot pick a representative')
            return next(iter(self._mixed_forms.values()))
        return self.get_form(nab)

    def set_name(self, name=None, latex_name=None):
        r"""
        Set (or change) the text name and LaTeX name of the characteristic
        class.

        INPUT:

        - ``name`` -- (string; default: ``None``) name given to the
          characteristic cohomology class
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote
          the characteristic cohomology class; if ``None`` while ``name`` is
          provided, the LaTeX symbol is set to ``name``

        EXAMPLES::

            sage: M = Manifold(8, 'M')
            sage: TM = M.tangent_bundle()
            sage: x = var('x')
            sage: k = TM.characteristic_cohomology_class(1+x^2,
            ....:                               class_type='multiplicative'); k
            Characteristic cohomology class (1 + p_1^2 - 2*p_2)(TM) of the
             Tangent bundle TM over the 8-dimensional differentiable manifold M
            sage: k.set_name(name='k', latex_name=r'\kappa')
            sage: k
            Characteristic cohomology class k(TM) of the Tangent bundle TM over
             the 8-dimensional differentiable manifold M
            sage: latex(k)
            \kappa\left(TM\right)
        """
        if name is not None:
            self._name = name
            if latex_name is None:
                self._latex_name = self._name
        if latex_name is not None:
            self._latex_name = latex_name


class CharacteristicCohomologyClassRing(FiniteGCAlgebra):
    r"""
    Characteristic cohomology class ring.

    Let `E \to M` be a real or complex vector bundle of rank `k` and `R` be a
    torsion-free subring of `\CC`.

    Let `BG` be the classifying space of the group `G`. As for vector bundles,
    we consider

    - `G = O(k)` if `E` is real,
    - `G = SO(k)` if `E` is real and oriented,
    - `G = U(k)` if `E` is complex.

    The cohomology ring `H^*(BG; R)` can be explicitly expressed for the
    aforementioned cases:

    .. MATH::

        H^*(BG; R) \cong \begin{cases}
            R[c_1, \ldots c_k] & \text{if } G = U(k), \\
            R[p_1, \ldots p_{\lfloor \frac{k}{2}\rfloor}] & \text{if } G = O(k), \\
            R[p_1, \ldots p_k, e] \big/ (p_k-e^2) & \text{if } G = SO(2k), \\
            R[p_1, \ldots p_k, e] & \text{if } G = SO(2k+1). \\
        \end{cases}

    The Chern-Weil homomorphism relates the generators in the de Rham
    cohomology as follows. If `\Omega` is a curvature form matrix on `E`, for
    the Chern classes we get

    .. MATH::

        \left[ \det\left( 1 + \frac{t \Omega}{2 \pi i} \right) \right] = 1 +
        \sum^k_{n=1} c_n(E) t^n,

    for the Pontryagin classes we have

    .. MATH::

        \left[ \det\left( 1 - \frac{t \Omega}{2 \pi} \right) \right] = 1 +
        \sum^{\lfloor\frac{k}{2} \rfloor}_{n=1} p_n(E) t^n,

    and for the Euler class we obtain

    .. MATH::

        \left[ \mathrm{Pf}\left(\frac{\Omega}{2 \pi} \right) \right] = e(E).

    Consequently, the cohomology ring `H^*(BG; R)` is mapped (not
    necessarily injectively) to a subring of `H^*_\mathrm{dR}(M, \CC)` via
    the Chern-Weil homomorphism. This implementation attempts to represent this
    subring.

    .. NOTE::

        Some generators might have torsion in `H^*(BG; R)` giving a zero
        element in the de Rham cohomology. Those generators are still
        considered in the implementation. Generators whose degree exceed the
        dimension of the base space, however, are ignored.

    INPUT:

    - ``base`` -- base ring
    - ``vbundle`` -- vector bundle

    EXAMPLES:

    Characteristic cohomology class ring over the tangent bundle of an
    8-dimensional manifold::

        sage: M = Manifold(8, 'M')
        sage: TM = M.tangent_bundle()
        sage: CR = TM.characteristic_cohomology_class_ring(); CR
        Algebra of characteristic cohomology classes of the Tangent bundle TM
         over the 8-dimensional differentiable manifold M
        sage: CR.gens()
        (Characteristic cohomology class (p_1)(TM) of the Tangent bundle TM over
         the 8-dimensional differentiable manifold M,
         Characteristic cohomology class (p_2)(TM) of the Tangent bundle TM
         over the 8-dimensional differentiable manifold M)

    The default base ring is `\QQ`::

        sage: CR.base_ring()
        Rational Field

    Characteristic cohomology class ring over a complex vector bundle::

        sage: M = Manifold(4, 'M')
        sage: E = M.vector_bundle(2, 'E', field='complex')
        sage: CR_E = E.characteristic_cohomology_class_ring(); CR_E
        Algebra of characteristic cohomology classes of the Differentiable
         complex vector bundle E -> M of rank 2 over the base space
         4-dimensional differentiable manifold M
        sage: CR_E.gens()
        (Characteristic cohomology class (c_1)(E) of the Differentiable complex
         vector bundle E -> M of rank 2 over the base space 4-dimensional
         differentiable manifold M,
         Characteristic cohomology class (c_2)(E) of the Differentiable
         complex vector bundle E -> M of rank 2 over the base space
         4-dimensional differentiable manifold M)

    Characteristic cohomology class ring over an oriented manifold::

        sage: S2 = manifolds.Sphere(2, coordinates='stereographic')
        sage: TS2 = S2.tangent_bundle()
        sage: S2.has_orientation()
        True
        sage: CR = TS2.characteristic_cohomology_class_ring()
        sage: CR.gens()
        (Characteristic cohomology class (e)(TS^2) of the Tangent bundle TS^2
         over the 2-sphere S^2 of radius 1 smoothly embedded in the Euclidean
         space E^3,)
    """
    Element = CharacteristicCohomologyClassRingElement

    def __init__(self, base, vbundle):
        r"""
        Construct a characteristic cohomology ring.

        TESTS::

            sage: M = Manifold(8, 'M')
            sage: TM = M.tangent_bundle()
            sage: CR = TM.characteristic_cohomology_class_ring()
            sage: TestSuite(CR).run()
        """
        self._vbundle = vbundle
        self._domain = vbundle._base_space
        dim = self._domain._dim
        rk = vbundle._rank
        if vbundle._field_type == 'complex':
            ran = min(rk, dim // 2)
            names = [f'c_{i}' for i in range(1, ran + 1)]
            degrees = [2 * i for i in range(1, ran + 1)]
            self._algorithm = ChernAlgorithm()
        elif vbundle._field_type == 'real':
            ran = min(rk // 2, dim // 4)
            names = [f'p_{i}' for i in range(1, ran + 1)]
            degrees = [4 * i for i in range(1, ran + 1)]
            self._algorithm = PontryaginAlgorithm()
            if vbundle.has_orientation():
                # add Euler class generator
                # Euler should be first entry; see `PontryaginEulerAlgorithm`
                names = ['e'] + names
                degrees = [rk] + degrees
                self._algorithm = PontryaginEulerAlgorithm()
                # TODO: add relation e^2=p_k for dim=2*k
        else:
            raise TypeError(f'Characteristic cohomology classes not supported '
                            f'for vector bundles with '
                            f'field type {vbundle._field_type}')

        if not names or not degrees:
            raise ValueError('cannot find any generators')

        names = tuple(names)  # hashable
        degrees = tuple(degrees)  # hashable
        super().__init__(base=base, names=names, degrees=degrees,
                         max_degree=dim, mul_symbol='⌣',
                         mul_latex_symbol=r'\smile')

    def _element_constructor_(self, x, **kwargs):
        r"""
        Convert ``x`` into ``self``.

        TESTS::

            sage: M = Manifold(8, 'M')
            sage: TM = M.tangent_bundle()
            sage: CR = TM.characteristic_cohomology_class_ring()
            sage: p = CR('Pontryagin'); p
            Characteristic cohomology class p(TM) of the Tangent bundle TM over
             the 8-dimensional differentiable manifold M
            sage: CR(p, name='pontr')
            Characteristic cohomology class pontr(TM) of the Tangent bundle
             TM over the 8-dimensional differentiable manifold M
        """
        if isinstance(x, (str, Expression)) or isinstance(x, Polynomial):
            return self._build_element(x, **kwargs)

        R = self.base_ring()

        if x in R:
            one_basis = self.one_basis()
            d = {one_basis: R(x)}
        elif isinstance(x, CharacteristicCohomologyClassRingElement):
            d = x._monomial_coefficients
        # x is an element of the basis enumerated set;
        # This is a very ugly way of testing this
        elif ((hasattr(self._indices, 'element_class') and
               isinstance(self._indices.element_class, type) and
               isinstance(x, self._indices.element_class)) or
              self.parent()(x) == self._indices):
            d = {x: R.one()}
        elif x in self._indices:
            d = {self._indices(x): R.one()}
        else:
            raise TypeError(f"do not know how to make x (= {x}) "
                            f"an element of self (={self})")
        name, latex_name = kwargs.get('name'), kwargs.get('latex_name')
        return self.element_class(self, d, name=name, latex_name=latex_name)

    @cached_method
    def _build_element(self, *args, **kwargs):
        r"""
        Construct a characteristic cohomology class.

        The result is cached.

        INPUT:

        - ``val`` -- the input data corresponding to the characteristic class
          using the Chern-Weil homomorphism; this argument can be either a
          symbolic expression, a polynomial or one of the following predefined
          classes:

          - ``'Chern'`` -- total Chern class,
          - ``'ChernChar'`` -- Chern character,
          - ``'Todd'`` -- Todd class,
          - ``'Pontryagin'`` -- total Pontryagin class,
          - ``'Hirzebruch'`` -- Hirzebruch class,
          - ``'AHat'`` -- `\hat{A}` class,
          - ``'Euler'`` -- Euler class.

        - ``name`` -- (default: ``None``) string representation given to the
          characteristic class; if ``None`` the default algebra representation or
          predefined name is used
        - ``latex_name`` -- (default: ``None``) LaTeX name given to the
          characteristic class; if ``None`` the value of ``name`` is used
        - ``class_type`` -- (default: ``None``) class type of the characteristic
          cohomology class; the following options are possible:

          - ``'multiplicative'`` -- returns a class of multiplicative type
          - ``'additive'`` -- returns a class of additive type
          - ``'Pfaffian'`` -- returns a class of Pfaffian type

          This argument must be stated if ``val`` is a polynomial or symbolic
          expression.

        EXAMPLES:

        Total Pontryagin class of an 8-dimensional manifold::

            sage: M = Manifold(8, 'M')
            sage: TM = M.tangent_bundle()
            sage: p = TM.characteristic_cohomology_class('Pontryagin'); p
            Characteristic cohomology class p(TM) of the Tangent bundle TM over the
             8-dimensional differentiable manifold M

        Define a multiplicative class (see :func:`multiplicative_sequence`)::

            sage: P.<x> = PolynomialRing(QQ)
            sage: f = 1 + x - x^2
            sage: f_class = TM.characteristic_cohomology_class(f, class_type='multiplicative'); f_class
            Characteristic cohomology class (1 + p_1 - p_1^2 + 3*p_2)(TM) of the
             Tangent bundle TM over the 8-dimensional differentiable manifold M

        Pass a symbolic expression, whose Taylor expansion at zero will be used::

            sage: M = Manifold(8, 'M')
            sage: TM = M.tangent_bundle()
            sage: x = var('x')
            sage: f = cos(x)
            sage: f_class = TM.characteristic_cohomology_class(f, class_type='multiplicative'); f_class
            Characteristic cohomology class (1 - 1/2*p_1^2 + p_2)(TM) of the Tangent
             bundle TM over the 8-dimensional differentiable manifold M
        """
        name, latex_name = kwargs.get('name'), kwargs.get('latex_name')
        base_ring = self.base_ring()
        class_type = kwargs.get('class_type')
        vbundle = self._vbundle
        val = args[0]
        dim = vbundle._base_space._dim

        # predefined classes accessible via class names
        if isinstance(val, str):
            from sage.arith.misc import factorial, bernoulli

            P = PolynomialRing(base_ring, 'x')
            x = P.gen()
            if val == 'Chern':
                if vbundle._field_type != 'complex':
                    raise ValueError(f'total Chern class not defined on {vbundle}')
                if name is None:
                    name = 'c'
                class_type = 'multiplicative'
                val = 1 + x
            elif val == 'Pontryagin':
                if vbundle._field_type != 'real':
                    raise ValueError(f'total Pontryagin class not defined on {vbundle}')
                if name is None:
                    name = 'p'
                class_type = 'multiplicative'
                val = 1 + x
            elif val == 'ChernChar':
                if vbundle._field_type != 'complex':
                    raise ValueError(f'Chern character not defined on {vbundle}')
                if name is None:
                    name = 'ch'
                if latex_name is None:
                    latex_name = r'\mathrm{ch}'
                class_type = 'additive'
                coeff = [1 / factorial(k) for k in
                         range(dim // 2 + 1)]  # exp(x)
                val = P(coeff)
            elif val == 'Todd':
                if vbundle._field_type != 'complex':
                    raise ValueError(f'Todd class not defined on {vbundle}')
                if name is None:
                    name = 'Td'
                if latex_name is None:
                    latex_name = r'\mathrm{Td}'
                class_type = 'multiplicative'
                val = 1 + x / 2
                for k in range(1, dim // 2 + 1):
                    val += (-1) ** (k + 1) / factorial(2 * k) * bernoulli(
                        2 * k) * x ** (2 * k)
            elif val == 'Hirzebruch':
                if vbundle._field_type != 'real':
                    raise ValueError(f'Hirzebruch class not defined on {vbundle}')
                if name is None:
                    name = 'L'
                if latex_name is None:
                    latex_name = 'L'
                class_type = 'multiplicative'
                coeff = [2 ** (2 * k) * bernoulli(2 * k) / factorial(2 * k)
                         for k in range(dim // 4 + 1)]
                val = P(coeff)
            elif val == 'AHat':
                if vbundle._field_type != 'real':
                    raise ValueError(f'AHat class not defined on {vbundle}')
                if name is None:
                    name = 'A^'
                if latex_name is None:
                    latex_name = r'\hat{A}'
                class_type = 'multiplicative'
                coeff = [- (2 ** (2 * k) - 2) / 2 ** (2 * k) * bernoulli(
                    2 * k) / factorial(2 * k)
                         for k in range(dim // 4 + 1)]
                val = P(coeff)
            elif val == 'Euler':
                if vbundle._field_type != 'real' or not vbundle.has_orientation():
                    raise ValueError(f'Euler class not defined on {vbundle}')
                if name is None:
                    name = 'e'
                class_type = 'Pfaffian'
                val = x
            else:
                raise ValueError(f'predefined class "{val}" unknown')

        # turn symbolic expression into a polynomial via Taylor expansion
        if isinstance(val, Expression):
            x = val.default_variable()
            P = PolynomialRing(base_ring, x)

            if vbundle._field_type == 'real':
                pow_range = dim // 4
            elif vbundle._field_type == 'complex':
                pow_range = dim // 2
            else:
                raise ValueError(f'field type of {vbundle} must be real or complex')

            val = P(val.taylor(x, 0, pow_range))

        # turn polynomial into a characteristic cohomology class via sequences
        if isinstance(val, Polynomial):
            if class_type is None:
                raise TypeError(f'class_type must be stated if {val} '
                                f'is a polynomial')
            n = self.ngens()
            s = 0  # shift; important in case of Euler class generator
            if self._algorithm is PontryaginEulerAlgorithm():
                s = 1  # skip Euler class
                n -= 1  # ignore Euler class

            if class_type == 'additive':
                sym = additive_sequence(val, vbundle._rank, n)
            elif class_type == 'multiplicative':
                sym = multiplicative_sequence(val, n)
            elif class_type == 'Pfaffian':
                P = val.parent()
                x = P.gen()
                val = (val(x) - val(-x)) / 2  # project to odd functions
                val = P([(-1) ** k * val[2 * k + 1] for k in range(n + 1)])
                sym = multiplicative_sequence(val, n)
            else:
                raise AttributeError('unkown class type')

            d = {}
            w_vec = self._weighted_vectors
            for p, c in sym:
                vec = [0] * self.ngens()
                if class_type == 'Pfaffian':
                    vec[0] = 1  # always multiply with e
                for i in p:
                    vec[i - 1 + s] += 1
                key = w_vec(vec)
                d[key] = c
            res = self._from_dict(d)
            res.set_name(name=name, latex_name=latex_name)
            return res

        # Nothing worked? Then something went wrong!
        raise ValueError(f'cannot convert {val} into an element of {self}')

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(8, 'M')
            sage: TM = M.tangent_bundle()
            sage: CR = TM.characteristic_cohomology_class_ring()
            sage: CR._repr_()
            'Algebra of characteristic cohomology classes of the Tangent bundle
             TM over the 8-dimensional differentiable manifold M'
            sage: CR  # indirect doctest
            Algebra of characteristic cohomology classes of the Tangent bundle
             TM over the 8-dimensional differentiable manifold M
        """
        vbundle = self._vbundle
        repr = f'Algebra of characteristic cohomology classes of the {vbundle}'
        return repr


# *****************************************************************************
# ALGORITHMS
# *****************************************************************************

def multiplicative_sequence(q, n=None):
    r"""
    Turn the polynomial ``q`` into its multiplicative sequence.

    Let `q` be a polynomial and `x_1, \ldots x_n` indeterminates. The
    *multiplicative sequence of* `q` is then given by the polynomials `K_j`

    .. MATH::

        \sum_{j=0}^n K_j(\sigma_1, \ldots, \sigma_j) z^j =
        \prod_{i=1}^{n} q(z \,x_i),

    where `\sigma_i` is the `i`-th elementary symmetric polynomial in the
    indeterminates `x_i`.

    INPUT:

    - ``q`` -- polynomial to turn into its multiplicative sequence
    - ``n`` -- (default: ``None``) the highest order `n` of the sequence;
      if ``None``, the order of ``q`` is assumed

    OUTPUT: a symmetric polynomial representing the multiplicative sequence

    EXAMPLES::

        sage: P.<x> = PolynomialRing(QQ)
        sage: from sage.manifolds.differentiable.characteristic_cohomology_class import multiplicative_sequence
        sage: f = 1 + x - x^2
        sage: sym = multiplicative_sequence(f); sym
        e[] + e[1] - e[1, 1] + 3*e[2]

    The maximal order of the result can be stated with ``n``::

        sage: sym_5 = multiplicative_sequence(f, n=5); sym_5
        e[] + e[1] - e[1, 1] + 3*e[2] - e[2, 1] + e[2, 2] + 4*e[3] - 3*e[3, 1]
         + e[3, 2] + 7*e[4] - 4*e[4, 1] + 11*e[5]
    """
    from sage.combinat.sf.sf import SymmetricFunctions
    from sage.combinat.partition import Partitions
    from sage.misc.misc_c import prod

    if n is None:
        n = q.degree()

    R = q.parent().base_ring()
    Sym = SymmetricFunctions(R)
    m = Sym.m()

    # Get the multiplicative sequence in the monomial basis:
    mon_pol = m._from_dict({p: prod(q[i] for i in p)
                            for k in range(n + 1)
                            for p in Partitions(k)})
    return Sym.e()(mon_pol)


def additive_sequence(q, k, n=None):
    r"""
    Turn the polynomial ``q`` into its additive sequence.

    Let `q` be a polynomial and `x_1, \ldots x_n` indeterminates. The
    *additive sequence of* `q` is then given by the polynomials `Q_j`

    .. MATH::

        \sum_{j=0}^n Q_j(\sigma_1, \ldots, \sigma_j) z^j =
        \sum_{i=1}^{k} q(z \,x_i),

    where `\sigma_i` is the `i`-th elementary symmetric polynomial in the
    indeterminates `x_i`.

    INPUT:

    - ``q`` -- polynomial to turn into its additive sequence
    - ``k`` -- maximal index `k` of the sum
    - ``n`` -- (default: ``None``) the highest order of the sequence `n`;
      if ``None``, the order of ``q`` is assumed

    OUTPUT: a symmetric polynomial representing the additive sequence

    EXAMPLES::

        sage: P.<x> = PolynomialRing(QQ)
        sage: from sage.manifolds.differentiable.characteristic_cohomology_class import additive_sequence
        sage: f = 1 + x - x^2
        sage: sym = additive_sequence(f, 2); sym
        2*e[] + e[1] - e[1, 1] + 2*e[2]

    The maximal order of the result can be stated with ``n``::

        sage: sym_1 = additive_sequence(f, 2, 1); sym_1
        2*e[] + e[1]
    """
    from sage.combinat.sf.sf import SymmetricFunctions
    from sage.combinat.partition import Partitions

    if n is None:
        n = q.degree()

    R = q.parent().base_ring()
    Sym = SymmetricFunctions(R)
    m = Sym.m()

    # Express the additive sequence in the monomial basis, the 0-th
    # order term must be treated separately; here comes ``rk`` into play:
    m_dict = {Partitions(0)([]): k * q[0]}
    m_dict.update({Partitions(k)([k]): q[k] for k in range(1, n + 1)})
    mon_pol = m._from_dict(m_dict)
    return Sym.e()(mon_pol)


def fast_wedge_power(form, n):
    r"""
    Return the wedge product power of `form` using a square-and-wedge
    algorithm.

    INPUT:

    - ``form`` -- a differential form
    - ``n`` -- nonnegative integer

    EXAMPLES::

        sage: M = Manifold(4, 'M')
        sage: X.<x,y,z,t> = M.chart()
        sage: omega = M.diff_form(2, name='omega')
        sage: omega[0,1] = t*y^2 + 2*x
        sage: omega[0,3] = z - 2*y
        sage: from sage.manifolds.differentiable.characteristic_cohomology_class import fast_wedge_power
        sage: fast_wedge_power(omega, 0)
        Scalar field 1 on the 4-dimensional differentiable manifold M
        sage: fast_wedge_power(omega, 1)
        2-form omega on the 4-dimensional differentiable manifold M
        sage: fast_wedge_power(omega, 2)
        4-form omega∧omega on the 4-dimensional differentiable manifold M
    """
    if n == 0:
        return form._domain._one_scalar_field
    elif n < 0:
        raise ValueError("'n' must be nonnegative")
    val = form
    while not (n & 1):
        val = val.wedge(val)
        n >>= 1

    # Now multiply together the correct factors form^(2^i)
    res = val
    n >>= 1
    while n:
        val = val.wedge(val)
        if n & 1:
            res = val.wedge(res)
        n >>= 1

    return res


class Algorithm_generic(SageObject):
    r"""
    Abstract algorithm class to compute the characteristic forms of the
    generators.

    EXAMPLES::

        sage: from sage.manifolds.differentiable.characteristic_cohomology_class import Algorithm_generic
        sage: algorithm = Algorithm_generic()
        sage: algorithm.get
        Cached version of <function Algorithm_generic.get at 0x...>
        sage: algorithm.get_local
        Traceback (most recent call last):
        ...
        NotImplementedError: <abstract method get_local at 0x...>
        sage: algorithm.get_gen_pow
        Cached version of <function Algorithm_generic.get_gen_pow at 0x...>
    """
    @cached_method
    def get(self, nab):
        r"""
        Return the global characteristic forms of the generators w.r.t. a given
        connection.

        The result is cached.

        OUTPUT:

        - a list containing the generator's global characteristic forms as
          instances of
          :class:`sage.manifolds.differentiable.diff_form.DiffForm`

        EXAMPLES::

            sage: M = manifolds.EuclideanSpace(4)
            sage: TM = M.tangent_bundle()
            sage: g = M.metric()
            sage: nab = g.connection()
            sage: nab.set_immutable()

        ::

            sage: p = TM.characteristic_cohomology_class('Pontryagin')
            sage: p_form = p.get_form(nab); p_form # long time
            Mixed differential form p(TE^4, nabla_g) on the 4-dimensional
             Euclidean space E^4
            sage: p_form.display_expansion() # long time
            p(TE^4, nabla_g) = 1
        """
        if isinstance(nab, AffineConnection):
            vbundle = nab._domain.tangent_bundle()
        elif isinstance(nab, BundleConnection):
            vbundle = nab._vbundle
        else:
            raise TypeError(f'{nab} must be a connection')
        dom = nab._domain
        res = []  # will be specified within first iteration
        for frame in dom._get_min_covering(nab._coefficients):
            cmat = [[nab.curvature_form(i, j, frame) for j in vbundle.irange()]
                    for i in vbundle.irange()]
            res_loc = self.get_local(cmat)
            if not res:
                # until now, degrees of generators were unknown
                res = [dom.diff_form(loc_form.degree())
                       for loc_form in res_loc]
            for form, loc_form in zip(res, res_loc):
                form.set_restriction(loc_form)
            # TODO: make `res` immutable?
        return res

    @abstract_method
    def get_local(self, cmat):
        r"""
        Abstract method to get the local forms of the generators w.r.t. a given
        curvature form matrix ``cmat``.

        OUTPUT: list containing the generator's local characteristic forms

        ALGORITHM:

        The inherited class determines the algorithm.

        EXAMPLES:

        4-dimensional Euclidean space::

            sage: M = manifolds.EuclideanSpace(4)
            sage: TM = M.tangent_bundle()
            sage: g = M.metric()
            sage: nab = g.connection()
            sage: e = M.frames()[0]  # select standard frame
            sage: cmat = [ [nab.curvature_form(i, j, e)  # long time
            ....:           for j in TM.irange()]
            ....:         for i in TM.irange()]

        Import the algorithm::

            sage: from sage.manifolds.differentiable.characteristic_cohomology_class import PontryaginAlgorithm
            sage: algorithm = PontryaginAlgorithm()
            sage: [p1] = algorithm.get_local(cmat) # long time
            sage: p1.display() # long time
            0

        A concrete implementation is given by a
        :class:`sage.misc.fast_methods.Singleton`::

            sage: algorithm is PontryaginAlgorithm()
            True
        """
        pass

    @cached_method
    def get_gen_pow(self, nab, i, n):
        r"""
        Return the `n`-th power of the `i`-th generator's characteristic form
        w.r.t ``nab``.

        The result is cached.

        EXAMPLES::

            sage: M = manifolds.EuclideanSpace(8)
            sage: TM = M.tangent_bundle()
            sage: g = M.metric()
            sage: nab = g.connection()
            sage: nab.set_immutable()

        ::

            sage: A = TM.characteristic_cohomology_class('AHat')
            sage: A_form = A.get_form(nab); A_form  # long time
            Mixed differential form A^(TE^8, nabla_g) on the 8-dimensional
             Euclidean space E^8
            sage: A_form.display_expansion()  # long time
            A^(TE^8, nabla_g) = 1
        """
        if n == 0:
            return nab._domain._one_scalar_field  # no computation necessary
        return fast_wedge_power(self.get(nab)[i], n)


class ChernAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Chern forms.

    EXAMPLES:

    Define a complex line bundle over a 2-dimensional manifold::

        sage: M = Manifold(2, 'M', structure='Lorentzian')
        sage: X.<t,x> = M.chart()
        sage: E = M.vector_bundle(1, 'E', field='complex'); E
        Differentiable complex vector bundle E -> M of rank 1 over the base space
         2-dimensional Lorentzian manifold M
        sage: e = E.local_frame('e')
        sage: nab = E.bundle_connection('nabla^E', latex_name=r'\nabla^E')
        sage: omega = M.one_form(name='omega')
        sage: A = function('A')
        sage: nab.set_connection_form(0, 0)[1] = I*A(t)
        sage: nab[0, 0].display()
        connection (0,0) of bundle connection nabla^E w.r.t. Local frame
         (E|_M, (e_0)) = I*A(t) dx
        sage: nab.set_immutable()

    Import the algorithm and apply ``nab`` to it::

        sage: from sage.manifolds.differentiable.characteristic_cohomology_class import ChernAlgorithm
        sage: algorithm = ChernAlgorithm()
        sage: algorithm.get(nab)
        [2-form on the 2-dimensional Lorentzian manifold M]
        sage: algorithm.get(nab)[0].display()
        1/2*d(A)/dt/pi dt∧dx

    Check some equalities::

        sage: cmat = [[nab.curvature_form(0, 0, e)]]
        sage: algorithm.get(nab)[0] == algorithm.get_local(cmat)[0]  # bundle trivial
        True
        sage: algorithm.get_gen_pow(nab, 0, 1) == algorithm.get(nab)[0]
        True
    """
    def get_local(self, cmat):
        r"""
        Return the local Chern forms w.r.t. a given curvature form matrix.

        OUTPUT:

        - a list containing the local characteristic Chern forms as
          instances of
          :class:`sage.manifolds.differentiable.diff_form.DiffForm`

        ALGORITHM::

            The algorithm is based on the Faddeev-LeVerrier algorithm for the
            characteristic polynomial.

        EXAMPLES:

        Define a complex line bundle over a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='Lorentzian')
            sage: X.<t,x> = M.chart()
            sage: E = M.vector_bundle(1, 'E', field='complex'); E
            Differentiable complex vector bundle E -> M of rank 1 over the base
             space 2-dimensional Lorentzian manifold M
            sage: e = E.local_frame('e')
            sage: nab = E.bundle_connection('nabla^E', latex_name=r'\nabla^E')
            sage: omega = M.one_form(name='omega')
            sage: A = function('A')
            sage: nab.set_connection_form(0, 0)[1] = I*A(t)
            sage: nab[0, 0].display()
            connection (0,0) of bundle connection nabla^E w.r.t. Local frame
             (E|_M, (e_0)) = I*A(t) dx
            sage: cmat = [[nab.curvature_form(i, j, e) for j in E.irange()]
            ....:         for i in E.irange()]

        Import the algorithm and apply ``cmat`` to it::

            sage: from sage.manifolds.differentiable.characteristic_cohomology_class import ChernAlgorithm
            sage: algorithm = ChernAlgorithm()
            sage: algorithm.get_local(cmat)
            [2-form on the 2-dimensional Lorentzian manifold M]
        """
        from sage.symbolic.constants import pi, I

        dom = cmat[0][0]._domain
        rk = len(cmat)
        dim = dom._dim
        ran = min(rk, dim // 2)
        if ran < 1:
            return []  # nothing to compute
        fac = I / (2 * pi)
        res = []
        m = cmat
        for k in range(1, ran):
            c = -sum(m[i][i] for i in range(rk)) / k
            res.append(fac * c)
            for i in range(rk):
                m[i][i] += c
            fac *= I / (2 * pi)
            m = [[sum(cmat[i][l].wedge(m[l][j]) for l in range(rk))
                  for j in range(rk)] for i in range(rk)]
        res.append(-fac * sum(m[i][i] for i in range(rk)) / ran)
        return res


class PontryaginAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Pontryagin forms.

    EXAMPLES:

    5-dimensional Euclidean space::

        sage: M = manifolds.EuclideanSpace(5)
        sage: g = M.metric()
        sage: nab = g.connection()
        sage: nab.set_immutable()

    Import the algorithm::

        sage: from sage.manifolds.differentiable.characteristic_cohomology_class import PontryaginAlgorithm
        sage: algorithm = PontryaginAlgorithm()
        sage: [p1] = algorithm.get(nab) # long time
        sage: p1.display() # long time
        0
    """
    def get_local(self, cmat):
        r"""
        Return the local Pontryagin forms w.r.t. a given curvature form matrix.

        OUTPUT: list containing the local characteristic Pontryagin forms

        ALGORITHM::

            The algorithm is based on the Faddeev-LeVerrier algorithm for the
            characteristic polynomial.

        EXAMPLES:

        5-dimensional Euclidean space::

            sage: M = manifolds.EuclideanSpace(5)
            sage: TM = M.tangent_bundle()
            sage: g = M.metric()
            sage: nab = g.connection()
            sage: e = M.frames()[0]  # select standard frame
            sage: cmat = [ [nab.curvature_form(i, j, e) # long time
            ....:           for j in TM.irange()]
            ....:         for i in TM.irange()]

        Import the algorithm::

            sage: from sage.manifolds.differentiable.characteristic_cohomology_class import PontryaginAlgorithm
            sage: algorithm = PontryaginAlgorithm()
            sage: [p1] = algorithm.get_local(cmat) # long time
            sage: p1.display() # long time
            0
        """
        from sage.symbolic.constants import pi

        dom = cmat[0][0]._domain
        rk = len(cmat)
        dim = dom._dim
        ran = min(rk // 2, dim // 4)
        if ran < 1:
            return []  # nothing to compute
        fac = 1 / (2 * pi) ** 2
        res = []
        m = cmat2 = [[sum(cmat[i][l].wedge(cmat[l][j])
                          for l in range(rk))
                      for j in range(rk)] for i in range(rk)]
        for k in range(1, ran):
            c = -sum(m[i][i] for i in range(rk)) / (2 * k)
            res.append(fac * c)
            for i in range(rk):
                m[i][i] += c
            fac *= 1 / (2 * pi) ** 2
            m = [[sum(cmat2[i][l].wedge(m[l][j]) for l in range(rk))
                  for j in range(rk)] for i in range(rk)]
        res.append(-fac * sum(m[i][i] for i in range(rk)) / (2 * ran))
        return res


class EulerAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Euler forms.

    EXAMPLES:

    Consider the 2-dimensional Euclidean space::

        sage: M = manifolds.EuclideanSpace(2)
        sage: g = M.metric()
        sage: nab = g.connection()
        sage: nab.set_immutable()

    Import the algorithm and apply ``nab`` to it::

        sage: from sage.manifolds.differentiable.characteristic_cohomology_class import EulerAlgorithm
        sage: algorithm = EulerAlgorithm()
        sage: algorithm.get(nab)
        [2-form on the Euclidean plane E^2]
        sage: algorithm.get(nab)[0].display()
        0
    """
    @cached_method
    def get(self, nab):
        r"""
        Return the global characteristic forms of the generators w.r.t. a given
        connection.

        INPUT:

        - ``nab`` -- a metric connection `\nabla`

        OUTPUT: list containing the global characteristic Euler form

        ALGORITHM:

        Assume that `\nabla` is compatible with the metric `g`, and let
        `(s_1, \ldots, s_n)` be any oriented frame. Denote by
        `G_s = (g(s_i, s_j))_{ij}` the metric tensor and let
        `\Omega_s` be the curvature form matrix of `\nabla` w.r.t. `s`.
        Then, we get:

        .. MATH::

            \left(G_s \cdot \Omega_s \right)_{ij} = g\!\left(R(.,.)s_i, s_j\right),

        where `R` is the Riemannian curvature tensor w.r.t. `\nabla`.

        The characteristic Euler form is now obtained by the expression

        .. MATH::

            \frac{1}{\sqrt{\left|\det(G_s)\right|}} \
                \mathrm{Pf}\!\left(G_s \cdot \frac{\Omega_s}{2 \pi}\right).

        EXAMPLES:

        Consider the 2-sphere::

            sage: M.<x,y> = manifolds.Sphere(2, coordinates='stereographic')
            sage: g = M.metric() # long time
            sage: nab = g.connection() # long time
            sage: nab.set_immutable() # long time

        Import the algorithm and apply ``nab`` to it::

            sage: from sage.manifolds.differentiable.characteristic_cohomology_class import EulerAlgorithm
            sage: algorithm = EulerAlgorithm()
            sage: algorithm.get(nab) # long time
            [2-form on the 2-sphere S^2 of radius 1 smoothly embedded in the
             Euclidean space E^3]
            sage: algorithm.get(nab)[0].display() # long time
            2/(pi + pi*x^4 + pi*y^4 + 2*pi*x^2 + 2*(pi + pi*x^2)*y^2) dx∧dy

        REFERENCES:

        - [Che1944]_
        - [Baer2020]_
        """
        if not isinstance(nab, LeviCivitaConnection):
            raise TypeError('Euler forms are currently only supported for '
                            'Levi-Civita connections')
        dom = nab._domain
        vbundle = dom.tangent_bundle()
        rk = vbundle._rank
        if not vbundle.has_orientation():
            raise ValueError('Euler forms can only be defined for orientable '
                             'vector bundles')
        if rk % 2 != 0:
            raise ValueError('Euler forms are currently only supported for '
                             'vector bundles with even rank')
        res = dom.diff_form(rk)
        g = nab._metric
        for frame in dom._get_min_covering(vbundle.orientation()):
            # (G_s * Ω_s)_ij = g(R(.,.)s_i, s_j)
            gcmat = [[sum(g[[frame, i, j]] * nab.curvature_form(j, k, frame)
                          for j in vbundle.irange())
                      for k in vbundle.irange()] for i in vbundle.irange()]
            [res_loc] = self.get_local(gcmat)  # Pf(G_s * Ω_s) mod const.
            # e = 1 / sqrt(|det(G_s)|) * Pf(G_s * Ω_s) mod const.
            det = g.det(frame)
            if det.is_trivial_zero():
                raise ValueError(f'metric {g} must be non-degenerate')
            sqrt_det = det.abs().sqrt()
            res.set_restriction(res_loc / sqrt_det)  # local Euler form
            # TODO: make `res` immutable?
        return [res]

    def get_local(self, cmat):
        r"""
        Return the normalized Pfaffian w.r.t. a given curvature form matrix.

        The normalization is given by the factor
        `\left(\frac{1}{2 \pi}\right)^{\frac{k}{2}}`, where `k` is the
        dimension of the curvature matrix.

        OUTPUT: list containing the normalized Pfaffian of a given curvature form

        .. NOTE::

            In general, the output does *not* represent the local
            characteristic Euler form. The result is only guaranteed to be the
            local Euler form if ``cmat`` is given w.r.t. an orthonormal
            oriented frame. See :meth:`get` for details.

        ALGORITHM::

            The algorithm is based on the Bär-Faddeev-LeVerrier algorithm for
            the Pfaffian.

        EXAMPLES:

        Consider the 2-sphere::

            sage: M.<th,phi> = manifolds.Sphere(2)  # use spherical coordinates
            sage: TM = M.tangent_bundle()
            sage: g = M.metric()
            sage: nab = g.connection()
            sage: e = M.frames()[0]  # select frame (opposite orientation!)
            sage: cmat = [[nab.curvature_form(i, j, e) for j in TM.irange()]
            ....:         for i in TM.irange()]

        We need some preprocessing because the frame is not orthonormal::

            sage: gcmat = [[sum(g[[e, i, j]] * nab.curvature_form(j, k, e)
            ....:               for j in TM.irange())
            ....:           for k in TM.irange()] for i in TM.irange()]

        Now, ``gcmat`` is guaranteed to be skew-symmetric and can be applied
        to :meth:`get_local`. Then, the result must be normalized::

            sage: from sage.manifolds.differentiable.characteristic_cohomology_class import EulerAlgorithm
            sage: algorithm = EulerAlgorithm()
            sage: euler = -algorithm.get_local(gcmat)[0] / sqrt(g.det(frame=e))
            sage: euler.display()
            1/2*sin(th)/pi dth∧dphi
        """
        from sage.symbolic.constants import pi

        rk = len(cmat)
        ran = rk // 2
        m = a = [cmat[i].copy() for i in range(rk)]
        for i in range(0, rk, 2):
            m[i], m[i + 1] = m[i + 1], m[i]  # swap entries
            for k in range(rk):
                m[k][i + 1] = -m[k][i + 1]
        for k in range(1, ran):
            e = -sum(m[i][i] for i in range(rk)) / (2 * k)
            for i in range(rk):
                m[i][i] += e
            m = [[sum(a[i][l].wedge(m[l][j]) for l in range(rk))
                  for j in range(rk)] for i in range(rk)]
        e = -sum(m[i][i] for i in range(rk)) / (2 * ran)  # Pfaffian mod sign
        e *= (-1 / (2 * pi)) ** ran  # normalize
        return [e]


class PontryaginEulerAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Euler and Pontryagin forms.

    EXAMPLES:

    6-dimensional Euclidean space::

        sage: M = manifolds.EuclideanSpace(6)
        sage: g = M.metric()
        sage: nab = g.connection()
        sage: nab.set_immutable()

    Import the algorithm::

        sage: from sage.manifolds.differentiable.characteristic_cohomology_class import PontryaginEulerAlgorithm
        sage: algorithm = PontryaginEulerAlgorithm()
        sage: e_form, p1_form = algorithm.get(nab)  # long time
        sage: e_form.display()  # long time
        0
        sage: p1_form.display()  # long time
        0
    """

    @cached_method
    def get(self, nab):
        r"""
        Return the global characteristic forms of the generators w.r.t. a given
        connection.

        OUTPUT:

        - a list containing the global Euler form in the first entry, and the
          global Pontryagin forms in the remaining entries.

        EXAMPLES:

        4-dimensional Euclidean space::

            sage: M = manifolds.EuclideanSpace(4)
            sage: g = M.metric()
            sage: nab = g.connection()
            sage: nab.set_immutable()

        Import the algorithm::

            sage: from sage.manifolds.differentiable.characteristic_cohomology_class import PontryaginEulerAlgorithm
            sage: algorithm = PontryaginEulerAlgorithm()
            sage: algorithm.get(nab) # long time
            [4-form on the 4-dimensional Euclidean space E^4,
             4-form on the 4-dimensional Euclidean space E^4]
        """
        return EulerAlgorithm().get(nab) + PontryaginAlgorithm().get(nab)

    def get_local(self, cmat):
        r"""
        Return the local Euler and Pontryagin forms w.r.t. a given curvature
        form matrix.

        .. NOTE::

            Similar as for :class:`EulerAlgorithm`, the first entry only
            represents the Euler form if the curvature form matrix is chosen
            w.r.t. an orthonormal oriented frame. See
            :meth:`EulerAlgorithm.get_local` for details.

        OUTPUT:

        - a list containing the local Euler form in the first entry, and the
          local Pontryagin forms in the remaining entries.

        EXAMPLES:

        6-dimensional Euclidean space::

            sage: M = manifolds.EuclideanSpace(6)
            sage: TM = M.tangent_bundle()
            sage: g = M.metric()
            sage: nab = g.connection()
            sage: e = M.frames()[0]  # select the standard frame
            sage: cmat = [ [nab.curvature_form(i, j, e)  # long time
            ....:           for j in TM.irange()]
            ....:         for i in TM.irange() ]

        Import the algorithm::

            sage: from sage.manifolds.differentiable.characteristic_cohomology_class import PontryaginEulerAlgorithm
            sage: algorithm = PontryaginEulerAlgorithm()
            sage: e, p1 = algorithm.get_local(cmat)  # long time
            sage: e.display()  # long time
            0
            sage: p1.display()  # long time
            0
        """
        res = EulerAlgorithm().get_local(cmat)  # first entry is Euler class
        res += PontryaginAlgorithm().get_local(cmat)  # rest Pontryagin
        return res

    @cached_method
    def get_gen_pow(self, nab, i, n):
        r"""
        Return the `n`-th power of the `i`-th generator w.r.t ``nab``.

        The result is cached.

        EXAMPLES:

        4-dimensional Euclidean space::

            sage: M = manifolds.EuclideanSpace(4)
            sage: TM = M.tangent_bundle()
            sage: g = M.metric()
            sage: nab = g.connection()
            sage: nab.set_immutable()

        Import the algorithm::

            sage: from sage.manifolds.differentiable.characteristic_cohomology_class import PontryaginEulerAlgorithm
            sage: algorithm = PontryaginEulerAlgorithm()
            sage: e = algorithm.get_gen_pow(nab, 0, 1)  # Euler, long time
            sage: e.display() # long time
            0
            sage: p1_pow2 = algorithm.get_gen_pow(nab, 1, 2)  # 1st Pontryagin squared, long time
            sage: p1_pow2 # long time
            8-form zero on the 4-dimensional Euclidean space E^4
        """
        if n == 0:
            return nab._domain._one_scalar_field  # no computation necessary
        if i == 0:
            return fast_wedge_power(EulerAlgorithm().get(nab)[0], n)
        return fast_wedge_power(PontryaginAlgorithm().get(nab)[i - 1], n)
