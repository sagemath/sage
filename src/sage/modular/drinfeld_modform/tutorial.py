r"""
Introduction to Drinfeld modular forms
======================================

This tutorial outlines the definitions, the notations, and the
implementation of Drinfeld modular forms in SageMath. We assume that the
reader has basic knowledge of classical modular forms, as we
will often make analogies to this setting. We also assume little
knowledge of Drinfeld modules. The interested reader can consult the
Sage reference manual about this topic :ref:`Drinfeld modules`.

.. RUBRIC:: Preliminary notations

Let `q` be a prime power and let `A` be the ring of functions of
`\mathbb{P}^1/\mathbb{F}_q` which are regular outside `\infty`. This
ring is the polynomial ring `\mathbb{F}_q[T]`. We denote by
`K := \mathbb{F}_q(T)` its fraction field. We endow `K` with the
`1/T`-adic valuation and let `K_{\infty} := \mathbb{F}_q((1/T))` be the
completion of `K`. Next, we define `\mathbb{C}_{\infty}` to
be the completion of an algebraic closure of `K_{\infty}`. Lastly, We
denote by `\tau : x\mapsto x^q` the `q`-Frobenius.

In SageMath, we create the rational function field by first creating a
univariate polynomial ring over `\mathbb{F}_q` and then constructing its
field of fractions::

    sage: A = GF(3)['T']
    sage: K.<T> = Frac(A)
    sage: K
    Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3
    sage: K.base()  # returns A
    Univariate Polynomial Ring in T over Finite Field of size 3

.. RUBRIC:: Drinfeld period domain

In the classical setting, the domain of any modular forms is the complex
upper half plane `\mathcal{H}:=\{w\in \mathbb{C} : \mathrm{im}(w)>0\}`.
The analogue of this plane in the function field setting is the
*Drinfeld period domain of rank* `r > 1` and it is defined by

.. MATH::

    \Omega^r(\mathbb{C}_{\infty}) :=
    \mathbb{P}^{r-1}(\mathbb{C}_{\infty})
    \setminus \{K_{\infty}\text{-rational hyperplanes}\}.

This space is a rigid analytic space and we identify its elements
with the set of column vectors
`(w_1,\ldots, w_{r-1}, w_{r})^{\mathrm{T}}` in `\mathbb{C}_{\infty}^r`
such that the `w_i` are `K_{\infty}`-linearly independant and
`w_r = \xi`, a nonzero constant in `\mathbb{C}_{\infty}`.

We define a left action of `\mathrm{GL}_r(K_{\infty})` on
`\Omega^r(\mathbb{C}_{\infty})` by setting

.. MATH::

    \gamma(w) := j(\gamma, w)^{-1}\gamma w

where `j(\gamma, w) := \xi^{-1} \cdot (\text{last entry of }\gamma w)`.

.. RUBRIC:: Universal Drinfeld module over `\Omega^r(\mathbb{C}_{\infty})`

For any `w = (w_1, \ldots, w_{r-1}, \xi)` in
`\Omega^r(\mathbb{C}_{\infty})` we have a corresponding discrete
`A`-module which is free of rank `r`:

.. MATH::

    \Lambda^w := Aw_1 \oplus \cdots \oplus Aw_{r-1} \oplus A\xi.

By analytic uniformization, there exists a corresponding Drinfeld module

.. MATH::

    \phi^w : T \mapsto T + g_1(w)\tau + \cdots
    + g_{r - 1}(w)\tau^{r-1} + g_{r}(w)\tau^{r}.

This Drinfeld module is called the *universal
Drinfeld `\mathbb{F}_q[T]`-module over `\Omega^r(\mathbb{C}_{\infty})`*
and its coefficients
`g_i : \Omega^r(\mathbb{C}_{\infty}) \rightarrow \mathbb{C}_{\infty}`
are rigid analytic functions satisfying the invariance property:

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

In the rank two case, this expansion at
infinity is of the form

.. MATH::

    g_i(w) = \sum_{i = 0}^{\infty} a_n(g_i)u(w)^i

where `u(w) := e(w)^{-1}` and `e` is the exponential of the Carlitz
module `\rho:T\mapsto T + \tau`. The analytic parameter `u` is called
the *parameter at infinity*.

A *Drinfeld modular form* of rank `r`, weight `k`, type `m` for
`\mathrm{GL}_r(A)` is a rigid analytic function
`f:\Omega^r(\mathbb{C}_{\infty}) \rightarrow \mathbb{C}_{\infty}`
such that

* `f(\gamma(w)) = \mathrm{det}(\gamma)^m j(\gamma, w)^k f(w)` for all
  `\gamma` in `\mathrm{GL}_r(A)` and
  `w\in \Omega^r(\mathbb{C}_{\infty})`;

* `f` is holomorphic at infinity.

The second condition is similar to the classical case. In the rank two
situation, this expansion is simply given by
`f = \sum_{n\geq 0} a_n(f) u^n` where `a_n(f)\in \mathbb{C}_{\infty}`.

The reader is refered to
part I of [BRP2018]_ for more information about the analytic theory of
Drinfeld modular form of arbitrary rank.

.. RUBRIC:: Ring of Drinfeld modular forms

Letting `M_k^{r, m}(\mathrm{GL}_r(A))` denote the space of rank `r`,
weight `k\in (q - 1)\mathbb{Z}` and type `m~(\mathrm{mod}~q-1)`
Drinfeld modular forms, we define

.. MATH::

    M^{r, 0}(\mathrm{GL}_r(A)) :=
    \bigoplus_{k\in ZZ} M_k^{r}(\mathrm{GL}_r(A))

to be the graded ring of all Drinfeld modular forms of type 0. The
graduation is given by the weight of a modular form. Similarly, we let
`M^{r}(\mathrm{GL}_r(A)) \supset M^{r, 0}(\mathrm{GL}_r(A))` be the ring
of all Drinfeld modular forms of arbitrary type. By
theorem 17.5 in part III of [BRP2018]_, we have

.. MATH::

    M^{r, 0}(\mathrm{GL}_r(A))
    = \mathbb{C}_{\infty}[g_1,\ldots, g_{r-1}, g_{r}].

and

.. MATH::

    M^{r, 0}(\mathrm{GL}_r(A))
    = \mathbb{C}_{\infty}[g_1,\ldots, g_{r-1}, h_{r}].

where `h_r` is a weight `(q^r - 1)/(q - 1)` modular forms of type `1`
which is a `(q-1)`-root of `g_r`.

To create the ring of type zero forms in SageMath we use the class
`DrinfeldModularForms`::

    sage: A = GF(3)['T']
    sage: K.<T> = Frac(A)
    sage: M = DrinfeldModularForms(K, 3)  # rank 3
    sage: M
    Ring of Drinfeld modular forms of rank 3 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3

This ring is generated by the coefficient forms that we see purely as
symbolic objects, each of which we attach its respective weight::

    sage: M.gens()
    [g1, g2, g3]
    sage: M.inject_variables()
    Defining g1, g2, g3
    sage: g1.weight()
    2
    sage: g2.weight()
    8
    sage: g3.weight()
    26

To create the ring of arbitrary types modular forms, one passes the
keyword argument ``has_type=True``::

    sage: M = DrinfeldModularForms(K, 4, has_type=True)
    sage: M.gens()
    [g1, g2, g3, h4]
    sage: h4 = M.3
    sage: h4.weight()
    40

.. RUBRIC:: Coefficients forms

There are bound methods for accessing the coefficients forms at any
`a`::

    sage: A = GF(2)['T']
    sage: K.<T> = Frac(A)
    sage: M = DrinfeldModularForms(K, 3)
    sage: M.coefficient_forms()
    [g1, g2, g3]
    sage: M.coefficient_forms(T^2)
    [(T^2 + T)*g1,
     g1^3 + (T^4 + T)*g2,
     g1^4*g2 + g1*g2^2 + (T^8 + T)*g3,
     g1^8*g3 + g1*g3^2 + g2^5,
     g2^8*g3 + g2*g3^4,
     g3^9]

To access a single coefficient form, use the method
``coefficient_form``::

    sage: A = GF(5)['T']
    sage: K.<T> = Frac(A)
    sage: M = DrinfeldModularForms(K, 3, has_type=True)
    sage: M.coefficient_form(1)
    g1
    sage: M.coefficient_form(2, T^2 + 1)
    g1^6 + (T^25 + T)*g2
    sage: M.coefficient_form(3)  # h3^(q - 1) = g3
    h3^4

.. RUBRIC:: References

See [BRP2018]_, [Gek1988]_.
"""
