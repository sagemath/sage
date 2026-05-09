r"""
Introduction to Drinfeld modular forms
======================================

This tutorial outlines the definitions, the notations, and the
implementation of Drinfeld modular forms in SageMath. We assume that the
reader has basic knowledge of classical modular forms, as we will often
make analogies to this setting. We also assume little knowledge of
Drinfeld modules; for this topic, the interested reader can consult the
SageMath reference manual
:ref:`Drinfeld modules <sage.rings.function_field.drinfeld_modules.drinfeld_module>`.

.. RUBRIC:: Preliminary notations

Let `q` be a prime power and let `A` be the ring of functions of
`\mathbb{P}^1/\mathbb{F}_q` which are regular outside a closed point
`\infty`. This ring is the polynomial ring `\mathbb{F}_q[T]`. We denote
by `K := \mathbb{F}_q(T)` rational function field. We endow `K` with the
`1/T`-adic valuation and let `K_{\infty} := \mathbb{F}_q((1/T))` be the
completion of `K`. Next, we define `\mathbb{C}_{\infty}` to be the
completion of an algebraic closure of `K_{\infty}`. Lastly, we denote
by `\tau : x\mapsto x^q` the `q`-Frobenius.

.. NOTE::

    The above construction of `\mathbb{C}_{\infty}` is the same as the
    construction of `\mathbb{C}_p` in the case of `p`-adic numbers
    (see :wikipedia:`P-adic_number#Algebraic_closure`).

In SageMath, we create the rational function field by first creating a
univariate polynomial ring over `\mathbb{F}_q` and, following this, by
constructing its field of fractions::

    sage: A = GF(3)['T']
    sage: K.<T> = Frac(A)
    sage: K
    Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3
    sage: K.base()  # returns A
    Univariate Polynomial Ring in T over Finite Field of size 3

.. RUBRIC:: Drinfeld period domain and action of `\mathrm{GL}_r(K_{\infty})`

In the classical setting, the domain of any modular form is the complex
upper half plane `\mathcal{H}:=\{w\in \mathbb{C} : \mathrm{im}(w)>0\}`.
The analogue of this plane in the function field setting is the
*Drinfeld period domain of rank* `r > 1` and it is defined by

.. MATH::

    \Omega^r(\mathbb{C}_{\infty}) :=
    \mathbb{P}^{r-1}(\mathbb{C}_{\infty})
    \setminus \{K_{\infty}\text{-rational hyperplanes}\}.

This space is a rigid analytic space and, after fixing an arbitrary
nonzero constant `\xi` in `\mathbb{C}_{\infty}`, we identify its
elements with the set of column vectors
`(w_1,\ldots, w_{r-1}, w_{r})^{\mathrm{T}}` in `\mathbb{C}_{\infty}^r`
such that the `w_i` are `K_{\infty}`-linearly independant and
`w_r = \xi`. Note that `\xi` is unspecified, but the reader can assume
that `\xi = 1` without any loss of significant information. Its value
can be interesting simply for normalization purposes.

We define a left action of `\mathrm{GL}_r(K_{\infty})` on
`\Omega^r(\mathbb{C}_{\infty})` by setting

.. MATH::

    \gamma(w) := j(\gamma, w)^{-1}\gamma w

where `j(\gamma, w) := \xi^{-1} \cdot (\text{last entry of }\gamma w)`.

.. RUBRIC:: Universal Drinfeld module over `\Omega^r(\mathbb{C}_{\infty})`

For any `w = (w_1, \ldots, w_{r-1}, \xi)` in
`\Omega^r(\mathbb{C}_{\infty})` we have a corresponding discrete
`A`-module `\Lambda^w` which is free of rank `r`:

.. MATH::

    \Lambda^w := Aw_1 \oplus \cdots \oplus Aw_{r-1} \oplus A\xi.

An important result is that we have analytic uniformization which is the
analogue of complex uniformization for elliptic curves. In our setting,
elliptic curves are replaced by Drinfeld modules. In short, there exists
a corresponding Drinfeld module

.. MATH::

    \phi^w : T \mapsto T + g_1(w)\tau + \cdots
    + g_{r - 1}(w)\tau^{r-1} + g_{r}(w)\tau^{r}.

such that the exponential of `\phi^w` induces an isomorphism (of abelian
group) between the additive group `\mathbb{C}_{\infty}` and the quotient
`\mathbb{C}_{\infty} / \Lambda^w`. Background material on Drinfeld
modules and their analytic uniformization can be found in section 4.3
and 4.6 of [Gos1998]_.

The Drinfeld module
`\phi^w : A \to \mathbb{C}_{\infty} \{\tau\}` is called the
*universal Drinfeld* `\mathbb{F}_q[T]`-*module over*
`\Omega^r(\mathbb{C}_{\infty})` and its coefficients
`g_i : \Omega^r(\mathbb{C}_{\infty}) \to \mathbb{C}_{\infty}`
are rigid analytic functions satisfying the *invariance property*:

.. MATH::

    g_i(\gamma(w)) = j(\gamma, w)^{1 - q^i} g_i(w),
    ~\forall \gamma\in \mathrm{GL}_r(A)

where `g_{r}(w)` never vanishes. Moreover, these coefficients `g_i`
admits an expansion at infinity, analogous to `q`-expansion principle
for classical modular forms. The functions `g_i` are known as the
*coefficients forms at* `T`. More generally, the coefficients of the
image `\phi^w_a` for any `a\in A` are called the *coefficient forms at*
`a` and they are an algebraic combination of the coefficient forms at
`T`.

In the rank two case, the expansion at infinity is of the form

.. MATH::

    g_i(w) = \sum_{i = 0}^{\infty} a_n(g_i)u(w)^i

where `u(w) := e(w)^{-1}` and `e` is the exponential function of the
Carlitz module `\rho:T\mapsto T + \tau`. The analytic parameter `u` is
called the *parameter at infinity*.

A *Drinfeld modular form* of rank `r`, weight `k`, type `m` for
`\mathrm{GL}_r(A)` is a rigid analytic function

.. MATH::

    f:\Omega^r(\mathbb{C}_{\infty}) \to \mathbb{C}_{\infty}

such that

* `f(\gamma(w)) = \mathrm{det}(\gamma)^m j(\gamma, w)^k f(w)` for all
  `\gamma` in `\mathrm{GL}_r(A)` and
  `w\in \Omega^r(\mathbb{C}_{\infty})`;

* `f` is holomorphic at infinity.

Without diving into the details, we mention that the second condition is
similar to the classical case. More specifically, in the rank two
situation, the expansion of `f` is given by a power series in `u`
`f = \sum_{n\geq 0} a_n(f) u^n` where `a_n(f)\in \mathbb{C}_{\infty}`.

Lastly, we also mention that the integer `m` only depends on its class
modulo `q-1`.

Note that all the above theory is covered in much greater details in
part I of [BRP2018]_.

.. RUBRIC:: Ring of Drinfeld modular forms

Letting `M_k^{r, m}(\mathrm{GL}_r(A))` denote the space of rank `r`,
weight `k\in (q - 1)\mathbb{Z}` and type `m` Drinfeld modular forms,
we define

.. MATH::

    M^{r, 0}(\mathrm{GL}_r(A)) :=
    \bigoplus_{k\in ZZ} M_k^{r, 0}(\mathrm{GL}_r(A))

to be the graded ring of all Drinfeld modular forms of type 0. The
graduation is given by the weight of a modular form. Similarly, we let
`M^{r}(\mathrm{GL}_r(A)) \supset M^{r, 0}(\mathrm{GL}_r(A))` be the ring
of all Drinfeld modular forms of rank `r` and arbitrary type. By
theorem 17.5 in part III of [BRP2018]_, we have

.. MATH::

    M^{r, 0}(\mathrm{GL}_r(A))
    = \mathbb{C}_{\infty}[g_1,\ldots, g_{r-1}, g_{r}].

and

.. MATH::

    M^{r}(\mathrm{GL}_r(A))
    = \mathbb{C}_{\infty}[g_1,\ldots, g_{r-1}, h_{r}].

where `h_r` is a weight `(q^r - 1)/(q - 1)` modular forms of type `1`
which is a `(q-1)`-root of `g_r` sometimes known as *Gekeler's* `h`
*function*, see theorem 3.8 of [Gek2017]_ for the precise definition of
this function.

.. RUBRIC:: SageMath implementation

In SageMath, we model the ring of type 0 Drinfeld modular forms over `K`
as a finitely generated ring in the coefficients forms `g_i`:

.. MATH::

    K[g_1,\ldots, g_{r-1}, g_r].

Hence, any ring element is seen as a formal algebraic combination of the
coefficient forms `g_i` over `K`. Likewise, the ring of arbitrary
type forms is generated by `g_1\ldots, g_{r-1}, h_r`.

To create the ring of type zero and rank `r` Drinfeld modular forms, one
uses the class
:class:`~sage.modular.drinfeld_modform.ring.DrinfeldModularForms`::

    sage: A = GF(3)['T']
    sage: K.<T> = Frac(A)
    sage: M = DrinfeldModularForms(K, 3)  # rank 3
    sage: M
    Ring of Drinfeld modular forms of rank 3 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3

To create the ring of arbitrary types modular forms, one passes the
keyword argument ``has_type=True``::

    sage: M = DrinfeldModularForms(K, 4, has_type=True)
    sage: M.gens()
    (g1, g2, g3, h4)
    sage: h4 = M.3
    sage: h4.weight()
    40

For more information about the functionalities of the implementation,
one should consult the documentation of the main classes:

- Parent class:
  :class:`~sage.modular.drinfeld_modform.ring.DrinfeldModularForms`

- Element class:
  :class:`~sage.modular.drinfeld_modform.element.DrinfeldModularFormsElement`

.. RUBRIC:: References

A good introduction to Drinfeld modular forms of rank 2, see Gekeler's
paper [Gek1988]_. See also [BRP2018]_ for a detailed exposition of the
arbitrary rank theory.
"""
