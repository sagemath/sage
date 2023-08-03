Integrable Highest Weight Representations of Affine Lie algebras
================================================================

.. linkall

In this section `\mathfrak{g}` can be an arbitrary Kac-Moody Lie Algebra
made with a symmetrizable, indecomposable Cartan matrix.

Suppose that `V` is a representation with a weight space decomposition
as in :ref:`roots_and_weights`.  Let `\alpha` be a real root, and let
`\mathfrak{g}_\alpha` be the corresponding root space, that is,
the one-dimensional weight space for `\alpha` in the adjoint
representation of `\mathfrak{g}` on itself. Then `-\alpha` is also
a root. The two one-dimensional spaces `\mathfrak{g}_\alpha` and
`\mathfrak{g}_{-\alpha}` generate a Lie algebra isomorphic to
`\mathfrak{sl}_2`. The module `V` is called *integrable* if for each
such `\alpha` the representation of `\mathfrak{sl}_2` obtained this
way integrates to a representation of the Lie group `\operatorname{SL}_2`.
Since this group contains an element that stabilizes `\mathfrak{h}`
and induces the corresponding simple reflection on the weight lattice,
integrability implies that the weight multiplicities are invariant
under the action of the Weyl group.

If the Kac-Moody Lie algebra `\mathfrak{g}` is finite-dimensional
then the integrable highest weight representations are
just the irreducible finite-dimensional ones. For a general
Kac-Moody Lie algebra the integrable highest weight representations
are the analogs of the finite-dimensional ones, that is,
with :class:`WeylCharacterRing` elements. Their
theory has many aspects in common with the finite-dimensional
representations of finite-dimensional simple Lie algebras,
such as the parametrization by dominant weights, and
generalizations of the Weyl denominator and character
formulas, due to Macdonald and Kac respectively.

If `\Lambda` is a dominant weight, then the irreducible
highest weight module `L(\Lambda)` defined in :ref:`roots_and_weights`
is integrable. Moreover every highest weight integrable representation
arises this way, so these representations are in bijection with the
cone `P^+` of dominant weights.

The affine case
---------------

Now we assume that `\mathfrak{g}` is affine. The integrable
highest weight representations (and their crystals) are
extremely interesting. Integrable highest weight representations of
`\mathfrak{g}` arise in a variety of contexts, from string
theory to the modular representation theory of the symmetric
groups, and the theory of modular forms. The representation `L(\Lambda_0)`
is particularly ubiquitous and is called the *basic representation*.

Therefore in [KMPS]_ (published in 1990) tabulated data for
many of these representations. They wrote

    We present a vast quantity of numerical data in tabular form, this
    being the only source for such information. The computations are tedious
    and not particularly straightforward when it is necessary to carry them
    out individually. We hope the appearance of this book will spur interest
    in a field that has become, in barely 20 years, deeply rewarding and
    full of promise for the future. It would indeed be gratifying if these
    tables were to appear to the scientists of 2040 as obsolete as the
    dust-gathering compilations of transcendental functions appear for us
    today because of their availability on every pocket calculator.

As we will explain, Sage can reproduce the contents of these tables. 
Moreover the tables in [KMPS]_ are limited to the untwisted types,
but Sage also implements the twisted types.

Although Sage can reproduce the tables in the second volume of [KMPS]_,
the work remains very useful. The first volume is a down-to-earth
and very helpful exposition of the theory of integrable representations
of affine Lie algebras with explicit examples and explanations of the
connections with mathematical physics and vertex operators.

.. _support_integrable:

The support of an integrable highest weight representation
----------------------------------------------------------

Let `\Lambda \in P^+` and let `V = L(\Lambda)` be the integrable
representation with highest weight `\Lambda`. If `\mu` is another
weight, let `\operatorname{mult}(\mu)` denote the multiplicity of
the weight `\mu` in `L(\Lambda)`. Define the *support* of the
representation `\operatorname{supp}(V)` to be the set of `\mu`
such that `\operatorname{mult}(\mu) > 0`.

If `\operatorname{mult}(\mu) > 0` then `\lambda-\mu` is a linear
combination of the simple roots with nonnegative integer coefficients.
Moreover `\operatorname{supp}(V)` is contained in the paraboloid

.. MATH::

    (\Lambda+\rho | \Lambda+\rho) - (\mu+\rho | \mu+\rho) \geq 0

where `(\, | \,)` is the invariant inner product on the weight
lattice and `\rho` is the Weyl vector (:ref:`untwisted_affine`).
Moreover if `\mu \in \operatorname{supp}(V)` then `\Lambda - \mu`
is an element of the root lattice `Q` ([Kac]_, Propositions 11.3 and 11.4).
    
We organize the weight multiplicities into sequences called
*string functions* or *strings* as follows. By [Kac]_, Proposition 11.3
or Corollary 11.9, for fixed `\mu` the function
`\operatorname{mult}(\mu - k\delta)` of `k` is an increasing sequence.
We adjust `\mu` by a multiple of `\delta` to the beginning
of the positive part of the sequence. Thus we define `\mu` to be
*maximal* if `\operatorname{mult}(\mu) \neq 0` but
`\operatorname{mult}(\mu + \delta) = 0`.

Since `\delta` is fixed under the action of the affine Weyl group, and
since the weight multiplicities are Weyl group invariant, the function
`k \mapsto \operatorname{mult}(\mu - k \delta)` is unchanged if `\mu`
is replaced by `w(\mu)` for some Weyl group element `w`. Now every
Weyl orbit contains a dominant weight.  Therefore in enumerating the
string we may assume that the weight `\mu` is dominant. There are only
a finite number of dominant maximal weights. Thus there are only a
finite number of such strings to be computed.

Modular Forms
-------------

Remarkably, [KacPeterson]_ showed that each string is the set of
Fourier coefficients of a weakly holomorphic modular form; see also
[Kac]_ Chapters 12 and 13. Here *weakly holomorphic modular* means
that the form is allowed to have poles at cusps.

To this end we define the *modular characteristic*

.. MATH::

    m_\Lambda = \frac{|\Lambda+\rho|^2}{2(k+h^\vee)} - \frac{|\rho|^2}{2h^\vee}.

Here `k = (\Lambda | \delta)` is the *level* of the representation and
`h^\vee` is the dual Coxeter number (:ref:`coxeternumber`).
If `\mu` is a weight, define

.. MATH::

    m_{\Lambda,\mu} = m_\Lambda - \frac{|\mu|^2}{2k}.

Let `\Lambda` be a weight, which we may assume maximal. Then Kac
and Peterson defined the *string function*

.. MATH::

    c_\mu^\Lambda = q^{m_{\Lambda,\mu}}
        \sum_{n\in\ZZ} \operatorname{mult}(\mu - n\delta) q^n.

Although these do arise as partition functions in string theory, the
term "string" here does not refer to physical strings.

The string function `c_\mu^\Lambda` is a weakly holomorphic modular
form, possibly of half-integral weight. See [Kac]_, Corollary 13.10,
or [KacPeterson]_. It can have poles at infinity, but multiplying
`c_\mu^\Lambda` by `\eta(\tau)^{\dim\,\mathfrak{g}^\circ}` gives
a holomorphic modular form (for some level). Here `\eta` is the
Dedekind eta function:

.. MATH::

    \eta(\tau) = q^{1/24} \prod_{k=1}^\infty(1-q^k),
        \qquad q = e^{2\pi i \tau}.

The weight of this modular form `\eta(\tau)^{\dim\,\mathfrak{g}^\circ}
c^\Lambda_\lambda` is the number of positive roots of `\mathfrak{g}^\circ`.

Sage methods for integrable representations
-------------------------------------------

In this section we will show how to use Sage to compute with
integrable highest weight representations of affine Lie algebras.
For further documentation, see the reference manual
:class:`~sage.combinat.root_system.integrable_representations.IntegrableRepresentation`.

In the following example, we work with the integrable representation
with highest weight `2 \Lambda_0` for `\widehat{\mathfrak{sl}}_2`,
that is, `A_1^{(1)}`. First we create a dominant weight in
the extended weight lattice, then create the ``IntegrableRepresentation``
class. We compute the strings. There are two, since there are two
dominant maximal weights. One of them is the highest weight `2\Lambda_0`,
and the other is `2\Lambda_1 - \delta`::

    sage: L = RootSystem("A1~").weight_lattice(extended=True)
    sage: Lambda = L.fundamental_weights()
    sage: delta = L.null_root()
    sage: W = L.weyl_group(prefix="s")
    sage: s0, s1 = W.simple_reflections()
    sage: V = IntegrableRepresentation(2*Lambda[0])
    sage: V.strings()
    {2*Lambda[0]: [1, 1, 3, 5, 10, 16, 28, 43, 70, 105, 161, 236],
     2*Lambda[1] - delta: [1, 2, 4, 7, 13, 21, 35, 55, 86, 130, 196, 287]}
    sage: mw1, mw2 = V.dominant_maximal_weights(); mw1, mw2
    (2*Lambda[0], 2*Lambda[1] - delta)

We see there are two dominant maximal weights, `2 \Lambda_0` and
`2 \Lambda_1 - \delta`. We obtain every maximal weight from these
by applying Weyl group elements. These lie inside the paraboloid
described in :ref:`support_integrable`. Here are a few more
maximal weights::

    sage: pairs = [(s0*s1*s0, mw1), (s0*s1, mw2), (s0, mw1), (W.one(), mw2),
    ....:          (W.one(), mw1), (s1, mw2), (s1*s0, mw1), (s1*s0*s1, mw2)]
    sage: [w.action(mw) for (w, mw) in pairs]
    [-6*Lambda[0] + 8*Lambda[1] - 8*delta,
     -4*Lambda[0] + 6*Lambda[1] - 5*delta,
     -2*Lambda[0] + 4*Lambda[1] - 2*delta,
     2*Lambda[1] - delta,
     2*Lambda[0],
     4*Lambda[0] - 2*Lambda[1] - delta,
     6*Lambda[0] - 4*Lambda[1] - 2*delta,
     8*Lambda[0] - 6*Lambda[1] - 5*delta]

We confirm that the string function for one in the Weyl orbit
is the same as that for ``mw2``, calculated above::

    sage: s1.action(mw2)
    4*Lambda[0] - 2*Lambda[1] - delta
    sage: [V.mult(s0.action(mw2)-k*delta) for k in [0..10]]
    [1, 2, 4, 7, 13, 21, 35, 55, 86, 130, 196]

String functions of integrable representations often appear
in the Online Encyclopedia of Integer Sequences::

    sage: [oeis(x) for x in V.strings().values()]    # optional - internet
    [0: A233758: Bisection of A006950 (the even part).,
     0: A233759: Bisection of A006950 (the odd part).]

Reading what the OEIS tells us about the sequence :oeis:`A006950`,
we learn that the two strings are the odd and even parts of the series

.. MATH::

   \prod_{k=1}^\infty \frac{1+q^{2k-1}}{1-q^{2k}}
   = \prod_{k=1}^\infty \frac{1-q^{2k}}{(1-q^k)(1-q^{4k})}
   = q^{1/8} \frac{\eta(2\tau)}{\eta(\tau)\eta(4\tau)}

This is *not* a modular form because of the factor `q^{1/8}` in
front of the ratio of eta functions.

Let us confirm what the Online Encyclopedia tells us by computing
the above product::

    sage: PS.<q> = PowerSeriesRing(QQ)
    sage: prod([(1+q^(2*k-1))/(1-q^(2*k)) for k in [1..20]])
    1 + q + q^2 + 2*q^3 + 3*q^4 + 4*q^5 + 5*q^6 + 7*q^7 + 10*q^8
     + 13*q^9 + 16*q^10 + 21*q^11 + 28*q^12 + 35*q^13 + 43*q^14
     + 55*q^15 + 70*q^16 + 86*q^17 + 105*q^18 + 130*q^19 + O(q^20)

We see the values of the two strings interspersed in this
product, with the `2 \Lambda_0` string values in the even
positions and the `2 \Lambda_1 - \delta` values in the odd positions.

To compute `c^{2\Lambda_0}_\lambda`, which is guaranteed to be
a modular form, we must compute the modular characteristics.
We are interested in the cases where `\lambda` is one of the
two dominant maximal weights::

     sage: [V.modular_characteristic(x) for x in [2*Lambda[0], 2*Lambda[1]-delta]]
     [-1/16, 7/16]

This gives us the string functions

.. MATH::

    \begin{aligned}
    c^{2\Lambda_0}_{2\Lambda_0} & = q^{-1/16}(1+q+3q^2+5q^3+10q^4+16q^5+\cdots),\\
    c^{2\Lambda_0}_{2\Lambda_1-\delta} & = q^{7/16}(1+2q+4q^2+7q^3+13q^4+21q^5+\cdots).
    \end{aligned}

These are both weakly holomorphic modular forms. Any linear combination
of these two is also a weakly holomorphic modular form. For example we
may replace `\tau` by `\tau/2` in our previous identity and get

.. MATH::

    c^{2\Lambda_0}_{2\Lambda_0} + c^{2\Lambda_0}_{2\Lambda_1-\delta}
    = \frac{\eta(\tau)}{\eta(\tau/2)\eta(2\tau)}.

Many more examples may be found in [KacPeterson]_ and [KMPS]_.

Let `V` be the integrable highest weight representation with highest
weight `\Lambda`. If `\mu` is in the support of `V` then `\Lambda - \mu`
is of the form `\sum_i n_i\alpha_i` where `\alpha_i` are the simple roots.
Sage employs an internal representation of the weights as tuples
`(n_0, n_1, \ldots)`. You can convert weights to and from this
notation as follows::

    sage: L = RootSystem(['E',6,2]).weight_lattice(extended=True)
    sage: Lambda = L.fundamental_weights()
    sage: delta = L.null_root()
    sage: V = IntegrableRepresentation(Lambda[0])
    sage: V.strings()
    {Lambda[0]: [1, 2, 7, 14, 35, 66, 140, 252, 485, 840, 1512, 2534]}
    sage: V.to_weight((1,2,0,1,0))
    Lambda[0] - 3*Lambda[1] + 4*Lambda[2] - 2*Lambda[3] + Lambda[4] - delta
    sage: V.from_weight(Lambda[0] - 3*Lambda[1] + 4*Lambda[2] - 2*Lambda[3] + Lambda[4] - delta)
    (1, 2, 0, 1, 0)
    sage: V.from_weight(Lambda[0]-delta)
    (1, 2, 3, 2, 1)

In reporting the strings, one may set the optional parameter depth to
get more or fewer values. In certain cases even the first coefficient
of the string is significant.  See [JayneMisra2014]_ and [KimLeeOh2017]_.

Catalan numbers (:oeis:`A000108`)::

    sage: P = RootSystem(['A',12,1]).weight_lattice(extended=true)
    sage: Lambda = P.fundamental_weights()
    sage: IntegrableRepresentation(2*Lambda[0]).strings(depth=1) # long time
    {2*Lambda[0]: [1],
     Lambda[1] + Lambda[12] - delta: [1],
     Lambda[2] + Lambda[11] - 2*delta: [2],
     Lambda[3] + Lambda[10] - 3*delta: [5],
     Lambda[4] + Lambda[9] - 4*delta: [14],
     Lambda[5] + Lambda[8] - 5*delta: [42],
     Lambda[6] + Lambda[7] - 6*delta: [132]}

Catalan triangle numbers (:oeis:`A000245`)::

    sage: sorted(IntegrableRepresentation(Lambda[0]+Lambda[2]).strings(depth=1).values()) # long time
    [[1], [3], [9], [12], [28], [90], [297]]

Central binomial coefficients (:oeis:`A001700`, :oeis:`A128015`)::

    sage: P = RootSystem(['B',8,1]).weight_lattice(extended=true)
    sage: Lambda = P.fundamental_weights()
    sage: sorted(IntegrableRepresentation(Lambda[0]+Lambda[1]).strings(depth=1).values()) # long time
    [[1], [1], [1], [3], [3], [10], [10], [35], [35], [126]]
