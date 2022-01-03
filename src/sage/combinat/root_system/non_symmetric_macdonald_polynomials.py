r"""
Nonsymmetric Macdonald polynomials

AUTHORS:

- Anne Schilling and Nicolas M. Thiéry (2013): initial version

ACKNOWLEDGEMENTS:

The initial version of this code (together with :class:`.root_lattice_realization_algebras.Algebras`
and :class:`.hecke_algebra_representation.HeckeAlgebraRepresentation`) was written by
Anne Schilling and Nicolas M. Thiery during the ICERM Semester Program on "Automorphic Forms,
Combinatorial Representation Theory and Multiple Dirichlet Series" (January 28, 2013 - May 3, 2013)
with the help of Dan Bump, Ben Brubaker, Bogdan Ion, Dan Orr, Arun Ram, Siddhartha Sahi, and Mark Shimozono.
Special thanks go to Bogdan Ion and Mark Shimozono for their patient explanations and hand computations
to check the code.
"""
# ****************************************************************************
#       Copyright (C) 2013 Nicolas M. Thiery <nthiery at users.sf.net>
#                          Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.integer_ring import ZZ
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.root_system.hecke_algebra_representation import CherednikOperatorsEigenvectors


class NonSymmetricMacdonaldPolynomials(CherednikOperatorsEigenvectors):
    r"""
    Nonsymmetric Macdonald polynomials

    INPUT:

    - ``KL`` -- an affine Cartan type or the group algebra of a
      realization of the affine weight lattice
    - ``q``, ``q1``, ``q2`` -- parameters in the base ring of the group algebra (default: ``q``, ``q1``, ``q2``)
    - ``normalized`` -- a boolean (default: ``True``)
      whether to normalize the result to have leading coefficient 1

    This implementation covers all reduced affine root systems.
    The polynomials are constructed recursively by the application
    of intertwining operators.

    .. TODO::

        - Non-reduced case (Koornwinder polynomials).
        - Non-equal parameters for the affine Hecke algebra.
        - Choice of convention (dominant/anti-dominant, ...).
        - More uniform implementation of the `T_0^\vee` operator.
        - Optimizations, in particular in the calculation of the
          eigenvalues for the recursion.

    EXAMPLES:

    We construct the family of nonsymmetric Macdonald polynomials in
    three variables in type `A`::

        sage: E = NonSymmetricMacdonaldPolynomials(["A",2,1])

    They are constructed as elements of the group algebra of the
    classical weight lattice `L_0` (or one of its realizations, such as
    the ambient space, which is used here) and indexed by elements of `L_0`::

        sage: L0 = E.keys(); L0
        Ambient space of the Root system of type ['A', 2]

    Here is the nonsymmetric Macdonald polynomial with leading term
    `[2,0,1]`::

        sage: E[L0([2,0,1])]
        ((-q*q1-q*q2)/(-q*q1-q2))*B[(1, 1, 1)] + ((-q1-q2)/(-q*q1-q2))*B[(2, 1, 0)] + B[(2, 0, 1)]

    It can be seen as a polynomial (or in general a Laurent
    polynomial) by interpreting each term as an exponent vector. The
    parameter `q` is the exponential of the null (co)root, whereas
    `q_1` and `q_2` are the two eigenvalues of each generator
    `T_i` of the affine Hecke algebra (see the background section for
    details).

    By setting `q_1=t`, `q_2=-1` and using the
    :meth:`.root_lattice_realization_algebras.Algebras.ElementMethods.expand`
    method, we recover the nonsymmetric Macdonald polynomial as
    computed by [HHL06]_'s combinatorial formula::

        sage: K = QQ['q,t'].fraction_field()
        sage: q,t = K.gens()
        sage: E = NonSymmetricMacdonaldPolynomials(["A",2,1], q=q, q1=t, q2=-1)
        sage: vars = K['x0,x1,x2'].gens()
        sage: E[L0([2,0,1])].expand(vars)
        (t - 1)/(q*t - 1)*x0^2*x1 + x0^2*x2 + (q*t - q)/(q*t - 1)*x0*x1*x2

        sage: from sage.combinat.sf.ns_macdonald import E
        sage: E([2,0,1])
        (t - 1)/(q*t - 1)*x0^2*x1 + x0^2*x2 + (q*t - q)/(q*t - 1)*x0*x1*x2

    Here is a type `G_2^{(1)}` nonsymmetric Macdonald polynomial::

        sage: E = NonSymmetricMacdonaldPolynomials(["G",2,1])
        sage: L0 = E.keys()
        sage: omega = L0.fundamental_weights()
        sage: E[ omega[2]-omega[1] ]
        ((-q*q1^3*q2-q*q1^2*q2^2)/(q*q1^4-q2^4))*B[(0, 0, 0)] + B[(1, -1, 0)] + ((-q1*q2^3-q2^4)/(q*q1^4-q2^4))*B[(1, 0, -1)]

    Many more examples are given after the background section.

    .. SEEALSO::

        - :func:`sage.combinat.sf.ns_macdonald.E`
        - :meth:`SymmetricFunctions.macdonald`

    .. RUBRIC:: Background

    .. RUBRIC:: The polynomial module

    The nonsymmetric Macdonald polynomials are a distinguished basis of the "polynomial" module
    of the affine Hecke algebra. Given::

        - a ground ring `K`, which contains the input parameters `q, q_1, q_2`
        - an affine root system, specified by a Cartan type `C`
        - a realization `L` of the weight lattice of type `C`

    the polynomial module is the group algebra `K[L_0]` of the classical
    weight lattice `L_0` with coefficients in `K`. It is isomorphic to the
    Laurent polynomial ring over `K` generated by the formal exponentials
    of any basis of `L_0`.

    In our running example `L` is the ambient space of type `C_2^{(1)}`::

        sage: K = QQ['q,q1,q2'].fraction_field()
        sage: q, q1, q2 = K.gens()
        sage: C = CartanType(["C",2,1])
        sage: L = RootSystem(C).ambient_space(); L
        Ambient space of the Root system of type ['C', 2, 1]

        sage: L.simple_roots()
        Finite family {0: -2*e[0] + e['delta'], 1: e[0] - e[1], 2: 2*e[1]}
        sage: omega = L.fundamental_weights(); omega
        Finite family {0: e['deltacheck'], 1: e[0] + e['deltacheck'], 2: e[0] + e[1] + e['deltacheck']}

        sage: L0 = L.classical(); L0
        Ambient space of the Root system of type ['C', 2]
        sage: KL0 = L0.algebra(K); KL0
        Algebra of the Ambient space of the Root system of type ['C', 2]
        over Fraction Field of Multivariate Polynomial Ring in q, q1, q2 over Rational Field

    .. RUBRIC:: Affine Hecke algebra

    The affine Hecke algebra is generated by elements `T_i` for ``i`` in the
    set of affine Dynkin nodes. They satisfy the same braid relations
    as the simple reflections `s_i` of the affine Weyl group.
    The `T_i` satisfy the quadratic relation

    .. MATH::

        (T_i-q_1)\circ(T_i-q_2) = 0,

    where `q_1` and `q_2` are the input parameters.  Some of the
    representation theory requires that `q_1` and `q_2` satisfy
    additional relations; typically one uses the specializations
    `q_1=u` and `q_2=-1/u` or `q_1=t` and `q_2=-1`). This can be
    achieved by constructing an appropriate field and passing `q_1`
    and `q_2` appropriately; see the examples. In principle, the
    parameter(s) could further depend on ``i``; this is not yet
    implemented but the code has been designed in such a way that this
    feature is easy to add.

    .. RUBRIC:: Demazure-Lusztig operators

    The ``i``-th Demazure-Lusztig operator is an operator on `K[L]`
    which interpolates between the reflection `s_i` and the Demazure operator `\pi_i`
    (see :meth:`.root_lattice_realization.RootLatticeRealization.Algebras.ParentMethods.demazure_lusztig_operators`).::

        sage: KL = L.algebra(K); KL
        Algebra of the Ambient space of the Root system of type ['C', 2, 1]
        over Fraction Field of Multivariate Polynomial Ring in q, q1, q2 over Rational Field
        sage: T = KL.demazure_lusztig_operators(q1, q2)
        sage: x = KL.monomial(omega[1]); x
        B[e[0] + e['deltacheck']]
        sage: T[2](x)
        q1*B[e[0] + e['deltacheck']]
        sage: T[1](x)
        (q1+q2)*B[e[0] + e['deltacheck']] + q1*B[e[1] + e['deltacheck']]
        sage: T[0](x)
        q1*B[e[0] + e['deltacheck']]

    The affine Hecke algebra acts on `K[L]` by letting the generators `T_i` act by
    the Demazure-Lusztig operators. The class
    :class:`sage.combinat.root_system.hecke_algebra_representation.HeckeAlgebraRepresentation`
    implements some simple generic features for representations of affine Hecke algebras
    defined by the action of their `T`-generators.::

        sage: T
        A representation of the (q1, q2)-Hecke algebra of type ['C', 2, 1] on Algebra of the Ambient space of the Root system of type ['C', 2, 1] over Fraction Field of Multivariate Polynomial Ring in q, q1, q2 over Rational Field
        sage: type(T)
        <class 'sage.combinat.root_system.hecke_algebra_representation.HeckeAlgebraRepresentation'>
        sage: T._test_relations()                 # long time (1.3s)

    Here we construct the operator `q_1 T_2^{-1}\circ T_1^{-1}T_0`
    from a signed reduced word::

        sage: T.Tw([0,1,2],[1,-1,-1], q1^2)
        Generic endomorphism of Algebra of the Ambient space of the Root system of type ['C', 2, 1]
        over Fraction Field of Multivariate Polynomial Ring in q, q1, q2 over Rational Field

    (note the reversal of the word). Inverses are computed using the
    quadratic relation.

    .. RUBRIC:: Cherednik operators

    The affine Hecke algebra contains elements `Y_\lambda` indexed by
    the coroot lattice. Their action on `K[L]` is implemented in Sage::

        sage: Y = T.Y(); Y
        Lazy family (...)_{i in Coroot lattice of the Root system of type ['C', 2, 1]}
        sage: alphacheck = Y.keys().simple_roots()
        sage: Y1 = Y[alphacheck[1]]
        sage: Y1(x)
        ((q1^2+2*q1*q2+q2^2)/(-q1*q2))*B[e[0] + e['deltacheck']]
        + ((-q1^2-2*q1*q2-q2^2)/(-q2^2))*B[-e[1] + e['deltacheck']]
        + ((-q1^2-q1*q2)/(-q2^2))*B[2*e[0] - e[1] - e['delta']
        + e['deltacheck']] + ((q1^3+q1^2*q2)/(-q2^3))*B[e[0] - e['delta']
        + e['deltacheck']] + ((q1^3+q1^2*q2)/(-q2^3))*B[e[0] - 2*e[1] - e['delta']
        + e['deltacheck']] + ((q1+q2)/(-q2))*B[e[1] + e['deltacheck']]
        + ((q1^3+2*q1^2*q2+q1*q2^2)/(-q2^3))*B[-e[1] - e['delta'] + e['deltacheck']]
        + ((q1^3+q1^2*q2)/(-q2^3))*B[2*e[0] - e[1] - 2*e['delta'] + e['deltacheck']]
        + ((q1^3+2*q1^2*q2+q1*q2^2)/(-q2^3))*B[-e[0] - e['delta'] + e['deltacheck']]
        + ((q1^3+2*q1^2*q2+q1*q2^2)/(-q2^3))*B[e[0] - 2*e['delta'] + e['deltacheck']]
        + ((q1^3+q1^2*q2)/(-q2^3))*B[3*e[0] - 3*e['delta'] + e['deltacheck']]
        + ((q1^3+q1^2*q2)/(-q2^3))*B[-e[0] - 2*e[1] - e['delta'] + e['deltacheck']]
        + ((q1^3+q1^2*q2)/(-q2^3))*B[e[0] - 2*e[1] - 2*e['delta'] + e['deltacheck']]
        + (q1^3/(-q2^3))*B[3*e[0] - 2*e[1] - 3*e['delta'] + e['deltacheck']]

    The Cherednik operators span a Laurent polynomial ring inside the
    affine Hecke algebra; namely `\lambda\mapsto Y_\lambda` is a group
    isomorphism from the classical root lattice (viewed additively) to
    the affine Hecke algebra (viewed multiplicatively). In practice,
    `Y_\lambda` is constructed by computing combinatorially its signed
    reduced word (and an overall scalar factor) using the periodic
    orientation of the alcove model in the coweight lattice (see
    :meth:`.hecke_algebra_representation.HeckeAlgebraRepresentation.Y_lambdacheck`)::

        sage: Lcheck = L.root_system.coweight_lattice()
        sage: w = Lcheck.reduced_word_of_translation(Lcheck(alphacheck[1])); w
        [0, 2, 1, 0, 2, 1]
        sage: Lcheck.signs_of_alcovewalk(w)
        [1, -1, 1, -1, 1, 1]

    .. RUBRIC:: Level zero representation of the affine Hecke algebra

    The action of the affine Hecke algebra on `K[L]` induces
    an action on `K[L_0]`: the action of `T_i` on `X^\lambda` for `\lambda` a
    classical weight in `L_0` is obtained by embedding the weight at
    level zero in the affine weight lattice (see
    :meth:`.weight_lattice_realizations.WeightLatticeRealizations.ParentMethods.embed_at_level`)
    applying the Demazure-Lusztig operator there, and projecting from `K[L]\to K[L_0]`
    mapping the exponential of `\delta` to `q` (see
    :meth:`.root_lattice_realization_algebras.Algebras.ParentMethods.q_project`). This is implemented in
    :meth:`.root_lattice_realization_algebras.Algebras.ParentMethods.demazure_lusztig_operators_on_classical`::

        sage: T = KL.demazure_lusztig_operators_on_classical(q, q1,q2)
        sage: omega = L0.fundamental_weights()
        sage: x = KL0.monomial(omega[1])
        sage: T[0](x)
        (-q*q2)*B[(-1, 0)]

    For classical nodes these are the usual Demazure-Lusztig operators::

        sage: T[1](x)
        (q1+q2)*B[(1, 0)] + q1*B[(0, 1)]

    .. RUBRIC:: Nonsymmetric Macdonald polynomials

    We can now finally define the nonsymmetric Macdonald polynomials.
    Because the Cherednik operators commute (and there is no radical),
    they can be simultaneously diagonalized; namely, `K[L_0]` admits a
    `K`-basis of joint eigenvectors for the `Y_\lambda`.  For `\mu \in
    L_0`, the nonsymmetric Macdonald polynomial `E_\mu` is the unique
    eigenvector of the family of Cherednik operators `Y_\lambda`
    having `\mu` as leading term::

        sage: E = NonSymmetricMacdonaldPolynomials(KL, q, q1, q2); E
        The family of the Macdonald polynomials of type ['C', 2, 1] with parameters q, q1, q2

    Or for short::

        sage: E = NonSymmetricMacdonaldPolynomials(C)

    .. RUBRIC:: Recursive construction of the nonsymmetric Macdonald polynomials

    The generators `T_i` of the affine Hecke algebra almost skew
    commute with the Cherednik operators. More precisely, one
    can deform them into the so-called intertwining operators:

    .. MATH:: \tau_i = T_i - (q_1+q_2) \frac{Y_i^{a-1}}{1-Y_i^a}\,.

    (where `a=1` except for `i=0` in type `BC` where `a=a_0=2`) which
    satisfy the following skew commutation relations:

    .. MATH:: \tau_i Y_\lambda = \tau_i Y_{s_i\lambda} \,.

    If `s_i \mu \ne \mu`, applying `\tau_i` on an eigenvector `E_\mu`
    produces a new eigenvector (essentially `E_{s_i\mu}`) with a
    distinct eigenvalue.  It follows that the eigenvectors indexed by
    an affine Weyl orbit of weights, may be recursively computed from
    a single weight in the orbit.

    In the case at hand, there is a little complication: namely, the
    simple reflections `s_i` acting at level 0 do not act transitively
    on classical weights; in fact the orbits for the classical Weyl
    group and for the affine Weyl group are the same. Thus, one can
    construct the nonsymmetric Macdonald polynomials for all weights
    from those for the classical dominant weights, but one is lacking
    a creation operator to construct the nonsymmetric Macdonald
    polynomials for dominant weights.

    .. RUBRIC:: Twisted Demazure-Lusztig operators

    To compensate for this, one needs to consider another affinization
    of the action of the classical Demazure-Lusztig operators
    `T_1,\dots,T_n`, which gives rise to the double affine Hecke algebra.
    Following Cherednik, one adds another operator `T_0^\vee` implemented in:
    :meth:`.root_lattice_realization_algebras.Algebras.ParentMethods.T0_check_on_basis`.
    See also:
    :meth:`.root_lattice_realization_algebras.Algebras.ParentMethods.twisted_demazure_lusztig_operators`.

    Depending on the type (untwisted or not), this is a representation
    of the affine Hecke algebra for another affinization of the
    classical Cartan type. The corresponding action of the affine Weyl
    group -- which is used to compute the recursion on `\mu` -- occurs
    in the corresponding weight lattice realization::

        sage: E.L()
        Ambient space of the Root system of type ['C', 2, 1]
        sage: E.L_prime()
        Coambient space of the Root system of type ['B', 2, 1]
        sage: E.L_prime().classical()
        Ambient space of the Root system of type ['C', 2]

    See :meth:`L_prime` and
    :meth:`.cartan_type.CartanType_affine.other_affinization`.

    REFERENCES:

    .. [HaimanICM] \M. Haiman, Cherednik algebras, Macdonald polynomials and combinatorics,
       Proceedings of the International Congress of Mathematicians,
       Madrid 2006, Vol. III, 843-872.

    .. [HHL06] \J. Haglund, M. Haiman and N. Loehr,
       A combinatorial formula for nonsymmetric Macdonald polynomials,
       Amer. J. Math. 130, No. 2 (2008), 359-383.

    .. [LNSSS12] \C. Lenart, S. Naito, D. Sagaki, A. Schilling, M. Shimozono,
       A uniform model for Kirillov-Reshetikhin crystals I: Lifting
       the parabolic quantum Bruhat graph, preprint :arXiv:`1211.2042`
       [math.QA]

    .. RUBRIC:: More examples

    We show how to create the nonsymmetric Macdonald polynomials in
    two different ways and check that they are the same::

        sage: K = QQ['q,u'].fraction_field()
        sage: q, u = K.gens()
        sage: E = NonSymmetricMacdonaldPolynomials(['D',3,1], q, u, -1/u)
        sage: omega = E.keys().fundamental_weights()
        sage: E[omega[1]+omega[3]]
        ((-q*u^2+q)/(-q*u^4+1))*B[(1/2, -1/2, 1/2)] + ((-q*u^2+q)/(-q*u^4+1))*B[(1/2, 1/2, -1/2)] + B[(3/2, 1/2, 1/2)]

        sage: KL = RootSystem(["D",3,1]).ambient_space().algebra(K)
        sage: P = NonSymmetricMacdonaldPolynomials(KL, q, u, -1/u)
        sage: E[omega[1]+omega[3]] == P[omega[1]+omega[3]]
        True
        sage: E[E.keys()((0,1,-1))]
        ((-q*u^2+q)/(-q*u^2+1))*B[(0, 0, 0)] + ((-u^2+1)/(-q*u^2+1))*B[(1, 1, 0)]
        + ((-u^2+1)/(-q*u^2+1))*B[(1, 0, -1)] + B[(0, 1, -1)]

    In type `A`, there is also a combinatorial implementation of the
    nonsymmetric Macdonald polynomials in terms of augmented diagram
    fillings as in [HHL06]_. See
    :func:`sage.combinat.sf.ns_macdonald.E`.  First we check that
    these polynomials are indeed eigenvectors of the Cherednik
    operators::

        sage: K = QQ['q,t'].fraction_field()
        sage: q,t = K.gens()
        sage: q1 = t; q2 = -1
        sage: KL = RootSystem(["A",2,1]).ambient_space().algebra(K)
        sage: KL0 = KL.classical()
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q, q1, q2)
        sage: omega = E.keys().fundamental_weights()
        sage: w = omega[1]
        sage: import sage.combinat.sf.ns_macdonald as NS
        sage: p = NS.E([1,0,0]); p
        x0
        sage: pp = KL0.from_polynomial(p)
        sage: E.eigenvalues(KL0.from_polynomial(p))
        [t, (-1)/(-q*t^2), t]

        sage: def eig(l): return E.eigenvalues(KL0.from_polynomial(NS.E(l)))

        sage: eig([1,0,0])
        [t, (-1)/(-q*t^2), t]
        sage: eig([2,0,0])
        [q*t, (-1)/(-q^2*t^2), t]
        sage: eig([3,0,0])
        [q^2*t, (-1)/(-q^3*t^2), t]
        sage: eig([2,0,4])
        [(-1)/(-q^3*t), 1/(q^2*t), q^4*t^2]

    Next we check explicitly that they agree with the current implementation::

        sage: K = QQ['q','t'].fraction_field()
        sage: q,t = K.gens()
        sage: KL = RootSystem(["A",1,1]).ambient_lattice().algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q, t, -1)
        sage: L0 = E.keys()
        sage: KL0 = KL.classical()
        sage: P = K['x0,x1']
        sage: def EE(weight): return E[L0(weight)].expand(P.gens())
        sage: import sage.combinat.sf.ns_macdonald as NS
        sage: EE([0,0])
        1
        sage: NS.E([0,0])
        1
        sage: EE([1,0])
        x0
        sage: NS.E([1,0])
        x0
        sage: EE([0,1])
        (t - 1)/(q*t - 1)*x0 + x1
        sage: NS.E([0,1])
        (t - 1)/(q*t - 1)*x0 + x1

        sage: NS.E([2,0])
        x0^2 + (q*t - q)/(q*t - 1)*x0*x1
        sage: EE([2,0])
        x0^2 + (q*t - q)/(q*t - 1)*x0*x1

    The same, directly in the ambient lattice with several shifts::

        sage: E[L0([2,0])]
        ((-q*t+q)/(-q*t+1))*B[(1, 1)] + B[(2, 0)]
        sage: E[L0([1,-1])]
        ((-q*t+q)/(-q*t+1))*B[(0, 0)] + B[(1, -1)]
        sage: E[L0([0,-2])]
        ((-q*t+q)/(-q*t+1))*B[(-1, -1)] + B[(0, -2)]

    Systematic checks with Sage's implementation of [HHL06]_::

        sage: assert all(EE([x,y]) == NS.E([x,y]) for d in range(5) for x,y in IntegerVectors(d,2))

    With the current implementation, we can compute nonsymmetric
    Macdonald polynomials for any type, for example for type `E_6^{(1)}`::

        sage: K=QQ['q,u'].fraction_field()
        sage: q, u = K.gens()
        sage: KL = RootSystem(["E",6,1]).weight_space(extended=True).algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q,u,-1/u)
        sage: L0 = E.keys()

        sage: E[L0.fundamental_weight(1).weyl_action([2,4,3,2,1])]
        ((-u^2+1)/(-q*u^16+1))*B[-Lambda[1] + Lambda[3]] + ((-u^2+1)/(-q*u^16+1))*B[Lambda[1]]
        + B[-Lambda[2] + Lambda[5]] + ((-u^2+1)/(-q*u^16+1))*B[Lambda[2] - Lambda[4] + Lambda[5]]
        + ((-u^2+1)/(-q*u^16+1))*B[-Lambda[3] + Lambda[4]]

        sage: E[L0.fundamental_weight(2).weyl_action([2,5,3,4,2])]  # long time (6s)
        ((-q^2*u^20+q^2*u^18+q*u^2-q)/(-q^2*u^32+2*q*u^16-1))*B[0]
        + B[Lambda[1] - Lambda[3] + Lambda[4] - Lambda[5] + Lambda[6]]
        + ((-u^2+1)/(-q*u^16+1))*B[Lambda[1] - Lambda[3] + Lambda[5]]
        + ((-q*u^20+q*u^18+u^2-1)/(-q^2*u^32+2*q*u^16-1))*B[-Lambda[2] + Lambda[4]]
        + ((-q*u^20+q*u^18+u^2-1)/(-q^2*u^32+2*q*u^16-1))*B[Lambda[2]]
        + ((u^4-2*u^2+1)/(q^2*u^32-2*q*u^16+1))*B[Lambda[3] - Lambda[4] + Lambda[5]]
        + ((-u^2+1)/(-q*u^16+1))*B[Lambda[3] - Lambda[5] + Lambda[6]]

        sage: E[L0.fundamental_weight(1)+L0.fundamental_weight(6)]  # long time (13s)
        ((q^2*u^10-q^2*u^8-q^2*u^2+q^2)/(q^2*u^26-q*u^16-q*u^10+1))*B[0]
        + ((-q*u^2+q)/(-q*u^10+1))*B[Lambda[1] - Lambda[2] + Lambda[6]]
        + ((-q*u^2+q)/(-q*u^10+1))*B[Lambda[1] + Lambda[2] - Lambda[4] + Lambda[6]]
        + ((-q*u^2+q)/(-q*u^10+1))*B[Lambda[1] - Lambda[3] + Lambda[4] - Lambda[5] + Lambda[6]]
        + ((-q*u^2+q)/(-q*u^10+1))*B[Lambda[1] - Lambda[3] + Lambda[5]] + B[Lambda[1] + Lambda[6]]
        + ((-q*u^2+q)/(-q*u^10+1))*B[-Lambda[2] + Lambda[4]] + ((-q*u^2+q)/(-q*u^10+1))*B[Lambda[2]]
        + ((-q*u^2+q)/(-q*u^10+1))*B[Lambda[3] - Lambda[4] + Lambda[5]]
        + ((-q*u^2+q)/(-q*u^10+1))*B[Lambda[3] - Lambda[5] + Lambda[6]]

    We test various other types::

        sage: K=QQ['q,u'].fraction_field()
        sage: q, u = K.gens()
        sage: KL = RootSystem(["A",5,2]).ambient_space().algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL, q, u, -1/u)
        sage: L0 = E.keys()
        sage: E[L0.fundamental_weight(2)]
        ((-q*u^2+q)/(-q*u^8+1))*B[(0, 0, 0)] + B[(1, 1, 0)]
        sage: E[L0((0,-1,1))]                                       # long time (1.5s)
        ((-q^2*u^10+q^2*u^8-q*u^6+q*u^4+q*u^2+u^2-q-1)/(-q^3*u^12+q^2*u^8+q*u^4-1))*B[(0, 0, 0)]
        + ((-u^2+1)/(-q*u^4+1))*B[(1, -1, 0)]
        + ((u^6-u^4-u^2+1)/(q^3*u^12-q^2*u^8-q*u^4+1))*B[(1, 1, 0)]
        + ((u^4-2*u^2+1)/(q^3*u^12-q^2*u^8-q*u^4+1))*B[(1, 0, -1)]
        + ((q^2*u^12-q^2*u^10-u^2+1)/(q^3*u^12-q^2*u^8-q*u^4+1))*B[(1, 0, 1)] + B[(0, -1, 1)]
        + ((-u^2+1)/(-q^2*u^8+1))*B[(0, 1, -1)] + ((-u^2+1)/(-q^2*u^8+1))*B[(0, 1, 1)]

        sage: K=QQ['q,u'].fraction_field()
        sage: q, u = K.gens()
        sage: KL = RootSystem(["E",6,2]).ambient_space().algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q,u,-1/u)
        sage: L0 = E.keys()
        sage: E[L0.fundamental_weight(4)]                           # long time (5s)
        ((-q^3*u^20+q^3*u^18+q^2*u^2-q^2)/(-q^3*u^28+q^2*u^22+q*u^6-1))*B[(0, 0, 0, 0)]
        + ((-q*u^2+q)/(-q*u^6+1))*B[(1/2, 1/2, -1/2, -1/2)] + ((-q*u^2+q)/(-q*u^6+1))*B[(1/2, 1/2, -1/2, 1/2)]
        + ((-q*u^2+q)/(-q*u^6+1))*B[(1/2, 1/2, 1/2, -1/2)] + ((-q*u^2+q)/(-q*u^6+1))*B[(1/2, 1/2, 1/2, 1/2)]
        + ((q*u^2-q)/(q*u^6-1))*B[(1, 0, 0, 0)] + B[(1, 1, 0, 0)] + ((-q*u^2+q)/(-q*u^6+1))*B[(0, 1, 0, 0)]
        sage: E[L0((1,-1,0,0))]                                     # long time (23s)
        ((q^3*u^18-q^3*u^16+q*u^4-q^2*u^2-2*q*u^2+q^2+q)/(q^3*u^18-q^2*u^12-q*u^6+1))*B[(0, 0, 0, 0)]
        + ((-q^3*u^18+q^3*u^16+q*u^2-q)/(-q^3*u^18+q^2*u^12+q*u^6-1))*B[(1/2, -1/2, -1/2, -1/2)]
        + ((-q^3*u^18+q^3*u^16+q*u^2-q)/(-q^3*u^18+q^2*u^12+q*u^6-1))*B[(1/2, -1/2, -1/2, 1/2)]
        + ((q^3*u^18-q^3*u^16-q*u^2+q)/(q^3*u^18-q^2*u^12-q*u^6+1))*B[(1/2, -1/2, 1/2, -1/2)]
        + ((q^3*u^18-q^3*u^16-q*u^2+q)/(q^3*u^18-q^2*u^12-q*u^6+1))*B[(1/2, -1/2, 1/2, 1/2)]
        + ((q*u^8-q*u^6-q*u^2+q)/(q^3*u^18-q^2*u^12-q*u^6+1))*B[(1/2, 1/2, -1/2, -1/2)]
        + ((q*u^8-q*u^6-q*u^2+q)/(q^3*u^18-q^2*u^12-q*u^6+1))*B[(1/2, 1/2, -1/2, 1/2)]
        + ((-q*u^8+q*u^6+q*u^2-q)/(-q^3*u^18+q^2*u^12+q*u^6-1))*B[(1/2, 1/2, 1/2, -1/2)]
        + ((-q*u^8+q*u^6+q*u^2-q)/(-q^3*u^18+q^2*u^12+q*u^6-1))*B[(1/2, 1/2, 1/2, 1/2)]
        + ((-q^2*u^18+q^2*u^16-q*u^8+q*u^6+q*u^2+u^2-q-1)/(-q^3*u^18+q^2*u^12+q*u^6-1))*B[(1, 0, 0, 0)]
        + B[(1, -1, 0, 0)] + ((-u^2+1)/(-q^2*u^12+1))*B[(1, 1, 0, 0)] + ((-u^2+1)/(-q^2*u^12+1))*B[(1, 0, -1, 0)]
        + ((u^2-1)/(q^2*u^12-1))*B[(1, 0, 1, 0)] + ((-u^2+1)/(-q^2*u^12+1))*B[(1, 0, 0, -1)]
        + ((-u^2+1)/(-q^2*u^12+1))*B[(1, 0, 0, 1)] + ((-q*u^2+q)/(-q*u^6+1))*B[(0, -1, 0, 0)]
        + ((-q*u^4+2*q*u^2-q)/(-q^3*u^18+q^2*u^12+q*u^6-1))*B[(0, 1, 0, 0)]
        + ((-q*u^4+2*q*u^2-q)/(-q^3*u^18+q^2*u^12+q*u^6-1))*B[(0, 0, -1, 0)]
        + ((-q*u^4+2*q*u^2-q)/(-q^3*u^18+q^2*u^12+q*u^6-1))*B[(0, 0, 1, 0)]
        + ((-q*u^4+2*q*u^2-q)/(-q^3*u^18+q^2*u^12+q*u^6-1))*B[(0, 0, 0, -1)]
        + ((-q*u^4+2*q*u^2-q)/(-q^3*u^18+q^2*u^12+q*u^6-1))*B[(0, 0, 0, 1)]

    Next we test a twisted type (checked against Maple computation by
    Bogdan Ion for `q_1=t^2` and `q_2=-1`)::

        sage: E = NonSymmetricMacdonaldPolynomials(["A",5,2])
        sage: omega = E.keys()

        sage: E[omega[1]]
        B[(1, 0, 0)]

        sage: E[-omega[1]]
        B[(-1, 0, 0)] + ((q*q1^6+q*q1^5*q2+q1*q2^5+q2^6)/(q^3*q1^6+q^2*q1^5*q2+q*q1*q2^5+q2^6))*B[(1, 0, 0)] + ((q1+q2)/(q*q1+q2))*B[(0, -1, 0)] + ((q1+q2)/(q*q1+q2))*B[(0, 1, 0)] + ((q1+q2)/(q*q1+q2))*B[(0, 0, -1)] + ((q1+q2)/(q*q1+q2))*B[(0, 0, 1)]

        sage: E[omega[2]]
        ((-q1*q2^3-q2^4)/(q*q1^4-q2^4))*B[(1, 0, 0)] + B[(0, 1, 0)]

        sage: E[-omega[2]]
        ((q^2*q1^7+q^2*q1^6*q2-q1*q2^6-q2^7)/(q^3*q1^7-q^2*q1^5*q2^2+q*q1^2*q2^5-q2^7))*B[(1, 0, 0)] + B[(0, -1, 0)]
        + ((q*q1^5*q2^2+q*q1^4*q2^3-q1*q2^6-q2^7)/(q^3*q1^7-q^2*q1^5*q2^2+q*q1^2*q2^5-q2^7))*B[(0, 1, 0)]
        + ((-q1*q2-q2^2)/(q*q1^2-q2^2))*B[(0, 0, -1)] + ((q1*q2+q2^2)/(-q*q1^2+q2^2))*B[(0, 0, 1)]

        sage: E[-omega[1]-omega[2]]
        ((q^3*q1^6+q^3*q1^5*q2+2*q^2*q1^6+3*q^2*q1^5*q2-q^2*q1^4*q2^2-2*q^2*q1^3*q2^3-q*q1^5*q2-2*q*q1^4*q2^2+q*q1^3*q2^3+2*q*q1^2*q2^4-q*q1*q2^5-q*q2^6+q1^3*q2^3+q1^2*q2^4-2*q1*q2^5-2*q2^6)/(q^4*q1^6+q^3*q1^5*q2-q^3*q1^4*q2^2+q*q1^2*q2^4-q*q1*q2^5-q2^6))*B[(0, 0, 0)] + B[(-1, -1, 0)] + ((q*q1^4+q*q1^3*q2+q1*q2^3+q2^4)/(q^3*q1^4+q^2*q1^3*q2+q*q1*q2^3+q2^4))*B[(-1, 1, 0)] + ((q1+q2)/(q*q1+q2))*B[(-1, 0, -1)] + ((-q1-q2)/(-q*q1-q2))*B[(-1, 0, 1)] + ((q*q1^4+q*q1^3*q2+q1*q2^3+q2^4)/(q^3*q1^4+q^2*q1^3*q2+q*q1*q2^3+q2^4))*B[(1, -1, 0)] + ((q^2*q1^6+q^2*q1^5*q2+q*q1^5*q2-q*q1^3*q2^3-q1^5*q2-q1^4*q2^2+q1^3*q2^3+q1^2*q2^4-q1*q2^5-q2^6)/(q^4*q1^6+q^3*q1^5*q2-q^3*q1^4*q2^2+q*q1^2*q2^4-q*q1*q2^5-q2^6))*B[(1, 1, 0)] + ((q*q1^4+2*q*q1^3*q2+q*q1^2*q2^2-q1^3*q2-q1^2*q2^2+q1*q2^3+q2^4)/(q^3*q1^4+q^2*q1^3*q2+q*q1*q2^3+q2^4))*B[(1, 0, -1)] + ((q*q1^4+2*q*q1^3*q2+q*q1^2*q2^2-q1^3*q2-q1^2*q2^2+q1*q2^3+q2^4)/(q^3*q1^4+q^2*q1^3*q2+q*q1*q2^3+q2^4))*B[(1, 0, 1)] + ((q1+q2)/(q*q1+q2))*B[(0, -1, -1)] + ((q1+q2)/(q*q1+q2))*B[(0, -1, 1)] + ((q*q1^4+2*q*q1^3*q2+q*q1^2*q2^2-q1^3*q2-q1^2*q2^2+q1*q2^3+q2^4)/(q^3*q1^4+q^2*q1^3*q2+q*q1*q2^3+q2^4))*B[(0, 1, -1)] + ((q*q1^4+2*q*q1^3*q2+q*q1^2*q2^2-q1^3*q2-q1^2*q2^2+q1*q2^3+q2^4)/(q^3*q1^4+q^2*q1^3*q2+q*q1*q2^3+q2^4))*B[(0, 1, 1)]

        sage: E[omega[1]-omega[2]]
        ((q^3*q1^7+q^3*q1^6*q2-q*q1*q2^6-q*q2^7)/(q^3*q1^7-q^2*q1^5*q2^2+q*q1^2*q2^5-q2^7))*B[(0, 0, 0)] + B[(1, -1, 0)]
        + ((q*q1^5*q2^2+q*q1^4*q2^3-q1*q2^6-q2^7)/(q^3*q1^7-q^2*q1^5*q2^2+q*q1^2*q2^5-q2^7))*B[(1, 1, 0)] + ((-q1*q2-q2^2)/(q*q1^2-q2^2))*B[(1, 0, -1)]
        + ((q1*q2+q2^2)/(-q*q1^2+q2^2))*B[(1, 0, 1)]

        sage: E[omega[3]]
        ((-q1*q2^2-q2^3)/(-q*q1^3-q2^3))*B[(1, 0, 0)] + ((-q1*q2^2-q2^3)/(-q*q1^3-q2^3))*B[(0, 1, 0)] + B[(0, 0, 1)]

        sage: E[-omega[3]]
        ((q*q1^4*q2+q*q1^3*q2^2-q1*q2^4-q2^5)/(-q^2*q1^5-q2^5))*B[(1, 0, 0)] + ((q*q1^4*q2+q*q1^3*q2^2-q1*q2^4-q2^5)/(-q^2*q1^5-q2^5))*B[(0, 1, 0)]
        + B[(0, 0, -1)] + ((-q1*q2^4-q2^5)/(-q^2*q1^5-q2^5))*B[(0, 0, 1)]

    .. RUBRIC:: Comparison with the energy function of crystals

    Next we test that the nonsymmetric Macdonald polynomials at `t=0`
    match with the one-dimensional configuration sums involving
    Kirillov-Reshetikhin crystals for various types. See
    [LNSSS12]_::

        sage: K = QQ['q,t'].fraction_field()
        sage: q,t = K.gens()
        sage: KL = RootSystem(["A",5,2]).ambient_space().algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL, q, t, -1)
        sage: omega = E.keys().fundamental_weights()
        sage: E[-omega[1]].map_coefficients(lambda x:x.subs(t=0))
        B[(-1, 0, 0)] + B[(1, 0, 0)] + B[(0, -1, 0)] + B[(0, 1, 0)] + B[(0, 0, -1)] + B[(0, 0, 1)]
        sage: E[-omega[2]].map_coefficients(lambda x:x.subs(t=0))   # long time (3s)
        (q+2)*B[(0, 0, 0)] + B[(-1, -1, 0)] + B[(-1, 1, 0)] + B[(-1, 0, -1)]
        + B[(-1, 0, 1)] + B[(1, -1, 0)] + B[(1, 1, 0)] + B[(1, 0, -1)] + B[(1, 0, 1)]
        + B[(0, -1, -1)] + B[(0, -1, 1)] + B[(0, 1, -1)] + B[(0, 1, 1)]

    ::

        sage: KL = RootSystem(["C",3,1]).ambient_space().algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q, t,-1)
        sage: omega = E.keys().fundamental_weights()
        sage: E[-omega[2]].map_coefficients(lambda x:x.subs(t=0))   # long time (5s)
        2*B[(0, 0, 0)] + B[(-1, -1, 0)] + B[(-1, 1, 0)] + B[(-1, 0, -1)]
        + B[(-1, 0, 1)] + B[(1, -1, 0)] + B[(1, 1, 0)] + B[(1, 0, -1)] + B[(1, 0, 1)]
        + B[(0, -1, -1)] + B[(0, -1, 1)] + B[(0, 1, -1)] + B[(0, 1, 1)]

    ::

        sage: R = RootSystem(['C',3,1])
        sage: KL = R.weight_lattice(extended=True).algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q, t,-1)
        sage: omega = E.keys().fundamental_weights()
        sage: La = R.weight_space().basis()
        sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1])
        sage: E[-2*omega[1]].map_coefficients(lambda x:x.subs(t=0)) == LS.one_dimensional_configuration_sum(q) # long time (15s)
        True
        sage: LS = crystals.ProjectedLevelZeroLSPaths(La[1]+La[2])
        sage: E[-omega[1]-omega[2]].map_coefficients(lambda x:x.subs(t=0)) == LS.one_dimensional_configuration_sum(q) # long time (45s)
        True

    ::

        sage: R = RootSystem(['C',2,1])
        sage: KL = R.weight_lattice(extended=True).algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q, t,-1)
        sage: omega = E.keys().fundamental_weights()
        sage: La = R.weight_space().basis()
        sage: for d in range(1,3):                                  # long time (10s)
        ....:     for x,y in IntegerVectors(d,2):
        ....:         weight = x*La[1]+y*La[2]
        ....:         weight0 = -x*omega[1]-y*omega[2]
        ....:         LS = crystals.ProjectedLevelZeroLSPaths(weight)
        ....:         assert E[weight0].map_coefficients(lambda x:x.subs(t=0)) == LS.one_dimensional_configuration_sum(q)

    ::

        sage: R = RootSystem(['B',3,1])
        sage: KL = R.weight_lattice(extended=True).algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q, t,-1)
        sage: omega = E.keys().fundamental_weights()
        sage: La = R.weight_space().basis()
        sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1])
        sage: E[-2*omega[1]].map_coefficients(lambda x:x.subs(t=0)) == LS.one_dimensional_configuration_sum(q) # long time (23s)
        True
        sage: B = crystals.KirillovReshetikhin(['B',3,1],1,1)
        sage: T = crystals.TensorProduct(B,B)
        sage: T.one_dimensional_configuration_sum(q) == LS.one_dimensional_configuration_sum(q) # long time (2s)
        True

    ::

        sage: R = RootSystem(['BC',3,2])
        sage: KL = R.weight_lattice(extended=True).algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q, t,-1)
        sage: omega = E.keys().fundamental_weights()
        sage: La = R.weight_space().basis()
        sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1])
        sage: E[-2*omega[1]].map_coefficients(lambda x:x.subs(t=0)) == LS.one_dimensional_configuration_sum(q) # long time (21s)
        True

    ::

        sage: R = RootSystem(CartanType(['BC',3,2]).dual())
        sage: KL = R.weight_space(extended=True).algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q, t,-1)
        sage: omega = E.keys().fundamental_weights()
        sage: La = R.weight_space().basis()
        sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1])
        sage: g = E[-2*omega[1]].map_coefficients(lambda x:x.subs(t=0)) # long time (30s)
        sage: f = LS.one_dimensional_configuration_sum(q)           # long time (1.5s)
        sage: P = g.support()[0].parent()                           # long time
        sage: B = P.algebra(q.parent())                             # long time
        sage: sum(p[1]*B(P(p[0])) for p in f) == g                  # long time
        True

    ::

        sage: C = CartanType(['G',2,1])
        sage: R = RootSystem(C.dual())
        sage: K = QQ['q,t'].fraction_field()
        sage: q,t = K.gens()
        sage: KL = R.weight_lattice(extended=True).algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL, q, t,-1)
        sage: omega = E.keys().fundamental_weights()
        sage: La = R.weight_space().basis()
        sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1])
        sage: E[-2*omega[1]].map_coefficients(lambda x:x.subs(t=0)) == LS.one_dimensional_configuration_sum(q) # not tested, long time (20s)
        True
        sage: LS = crystals.ProjectedLevelZeroLSPaths(La[1]+La[2])
        sage: E[-omega[1]-omega[2]].map_coefficients(lambda x:x.subs(t=0)) == LS.one_dimensional_configuration_sum(q) # not tested, long time (23s)
        True

    The next test breaks if the energy is not scaled by the
    translation factor for dual type `G_2^{(1)}`::

        sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1]+La[2])
        sage: E[-2*omega[1]-omega[2]].map_coefficients(lambda x:x.subs(t=0)) == LS.one_dimensional_configuration_sum(q) # not tested, very long time (100s)
        True

        sage: R = RootSystem(['D',4,1])
        sage: KL = R.weight_lattice(extended=True).algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL, q, t,-1)
        sage: omega = E.keys().fundamental_weights()
        sage: La = R.weight_space().basis()
        sage: for d in range(1,2):                                  # long time (41s)
        ....:     for a,b,c,d in IntegerVectors(d,4):
        ....:         weight = a*La[1]+b*La[2]+c*La[3]+d*La[4]
        ....:         weight0 = -a*omega[1]-b*omega[2]-c*omega[3]-d*omega[4]
        ....:         LS = crystals.ProjectedLevelZeroLSPaths(weight)
        ....:         assert E[weight0].map_coefficients(lambda x:x.subs(t=0)) == LS.one_dimensional_configuration_sum(q)

    TESTS:

    Calculations checked with Bogdan Ion 2013/04/18::

        sage: K = QQ['q,t'].fraction_field()
        sage: q,t=K.gens()
        sage: E = NonSymmetricMacdonaldPolynomials(["B",2,1], q=q,q1=t,q2=-1/t)
        sage: L0 = E.keys()
        sage: omega = L0.fundamental_weights()

        sage: E[omega[1]]
        ((-q*t^4+q*t^2)/(-q*t^6+1))*B[(0, 0)] + B[(1, 0)]
        sage: E[omega[2]]
        B[(1/2, 1/2)]
        sage: E[-omega[1]]
        ((-q^2*t^8+q^2*t^6-q*t^6+2*q*t^4-q*t^2+t^2-1)/(-q^3*t^8+q^2*t^6+q*t^2-1))*B[(0, 0)] + B[(-1, 0)] + ((-q*t^8+q*t^6+t^2-1)/(-q^3*t^8+q^2*t^6+q*t^2-1))*B[(1, 0)] + ((-t^2+1)/(-q*t^2+1))*B[(0, -1)] + ((t^2-1)/(q*t^2-1))*B[(0, 1)]
        sage: E[L0([0,1])]
        ((-q*t^4+q*t^2)/(-q*t^4+1))*B[(0, 0)] + ((-t^2+1)/(-q*t^4+1))*B[(1, 0)] + B[(0, 1)]
        sage: E[L0([1,1])]
        ((q*t^2-q)/(q*t^2-1))*B[(0, 0)] + ((-q*t^2+q)/(-q*t^2+1))*B[(1, 0)] + B[(1, 1)] + ((-q*t^2+q)/(-q*t^2+1))*B[(0, 1)]

        sage: E = NonSymmetricMacdonaldPolynomials(["A",2,1], q=q,q1=t,q2=-1/t)
        sage: L0 = E.keys()
        sage: factor(E[L0([-1,0,1])][L0.zero()])
        (t - 1) * (t + 1) * (q*t^2 - 1)^-3 * (q*t^2 + 1)^-1 * (q^3*t^6 + 2*q^2*t^6 - 3*q^2*t^4 - 2*q*t^2 - t^2 + q + 2)

    Checking step by step calculations in type `BC` with Bogdan Ion 2013/04/18::

        sage: K = QQ['q,t'].fraction_field()
        sage: q,t=K.gens()
        sage: E = NonSymmetricMacdonaldPolynomials(["BC",1,2], q=q,q1=t,q2=-1/t)
        sage: KL0 = E.domain()
        sage: L0 = E.keys()
        sage: omega = L0.fundamental_weights()
        sage: e = L0.basis()
        sage: E._T_Y[1] ( KL0.monomial(e[0]) )
        1/t*B[(-1)]
        sage: E._T_Y[0] ( KL0.monomial(L0.zero()) )
        t*B[(0)]
        sage: E._T_Y[0] ( KL0.monomial(-e[0]))
        ((-t^2+1)/(q*t))*B[(0)] + 1/(q^2*t)*B[(1)]

        sage: Y = E.Y()
        sage: alphacheck = Y.keys().simple_roots()
        sage: Y0 = Y[alphacheck[0]]
        sage: Y1 = Y[alphacheck[1]]
        sage: Y0
        Generic endomorphism of Algebra of the Ambient space of the Root system of type ['C', 1]
        over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field
        sage: Y0.word, Y0.signs, Y0.scalar
        ((0, 1), (-1, -1), 1/q)
        sage: Y1.word, Y1.signs, Y1.scalar
        ((1, 0), (1, 1), 1)

        sage: T0_check = E._T[0]

    Comparing with Bogdan Ion's hand calculations for type `BC`, 2013/05/13:

    .. TODO:: add his notes in latex

    ::

        sage: K = QQ['q,q1,q2'].fraction_field()
        sage: q,q1,q2=K.gens()
        sage: L = RootSystem(["A",4,2]).ambient_space()
        sage: L.cartan_type()
        ['BC', 2, 2]
        sage: L.null_root()
        2*e['delta']
        sage: L.simple_roots()
        Finite family {0: -e[0] + e['delta'], 1: e[0] - e[1], 2: 2*e[1]}
        sage: KL = L.algebra(K)
        sage: KL0 = KL.classical()
        sage: L0 = L.classical()
        sage: L0.cartan_type()
        ['C', 2]

        sage: E = NonSymmetricMacdonaldPolynomials(KL, q=q,q1=q1,q2=q2)
        sage: E.keys()
        Ambient space of the Root system of type ['C', 2]
        sage: E.keys().simple_roots()
        Finite family {1: (1, -1), 2: (0, 2)}
        sage: omega = E.keys().fundamental_weights()

        sage: E[0*omega[1]]
        B[(0, 0)]
        sage: E[omega[1]]
        ((-q*q1*q2^3-q*q2^4)/(q^2*q1^4-q2^4))*B[(0, 0)] + B[(1, 0)]

        sage: E[2*omega[2]]      # long time # not checked against Bogdan's notes, but a good self-consistency test
        ((-q^12*q1^6-q^12*q1^5*q2+2*q^10*q1^5*q2+5*q^10*q1^4*q2^2+3*q^10*q1^3*q2^3+2*q^8*q1^5*q2+4*q^8*q1^4*q2^2+q^8*q1^3*q2^3-q^8*q1^2*q2^4+q^8*q1*q2^5+q^8*q2^6-q^6*q1^3*q2^3+q^6*q1^2*q2^4+4*q^6*q1*q2^5+2*q^6*q2^6+q^4*q1^3*q2^3+3*q^4*q1^2*q2^4+4*q^4*q1*q2^5+2*q^4*q2^6)/(-q^12*q1^6-q^10*q1^5*q2-q^8*q1^3*q2^3+q^6*q1^4*q2^2-q^6*q1^2*q2^4+q^4*q1^3*q2^3+q^2*q1*q2^5+q2^6))*B[(0, 0)] + ((q^7*q1^2*q2+2*q^7*q1*q2^2+q^7*q2^3+q^5*q1^2*q2+2*q^5*q1*q2^2+q^5*q2^3)/(-q^8*q1^3-q^6*q1^2*q2+q^2*q1*q2^2+q2^3))*B[(-1, 0)] + ((-q^6*q1*q2-q^6*q2^2)/(q^6*q1^2-q2^2))*B[(-1, -1)] + ((q^6*q1^2*q2+2*q^6*q1*q2^2+q^6*q2^3+q^4*q1^2*q2+2*q^4*q1*q2^2+q^4*q2^3)/(-q^8*q1^3-q^6*q1^2*q2+q^2*q1*q2^2+q2^3))*B[(-1, 1)] + ((-q^3*q1*q2-q^3*q2^2)/(q^6*q1^2-q2^2))*B[(-1, 2)] + ((q^7*q1^3+q^7*q1^2*q2-q^7*q1*q2^2-q^7*q2^3-2*q^5*q1^2*q2-4*q^5*q1*q2^2-2*q^5*q2^3-2*q^3*q1^2*q2-4*q^3*q1*q2^2-2*q^3*q2^3)/(q^8*q1^3+q^6*q1^2*q2-q^2*q1*q2^2-q2^3))*B[(1, 0)] + ((q^6*q1^2*q2+2*q^6*q1*q2^2+q^6*q2^3+q^4*q1^2*q2+2*q^4*q1*q2^2+q^4*q2^3)/(-q^8*q1^3-q^6*q1^2*q2+q^2*q1*q2^2+q2^3))*B[(1, -1)] + ((q^8*q1^3+q^8*q1^2*q2+q^6*q1^3+q^6*q1^2*q2-q^6*q1*q2^2-q^6*q2^3-2*q^4*q1^2*q2-4*q^4*q1*q2^2-2*q^4*q2^3-q^2*q1^2*q2-3*q^2*q1*q2^2-2*q^2*q2^3)/(q^8*q1^3+q^6*q1^2*q2-q^2*q1*q2^2-q2^3))*B[(1, 1)] + ((q^5*q1^2+q^5*q1*q2-q^3*q1*q2-q^3*q2^2-q*q1*q2-q*q2^2)/(q^6*q1^2-q2^2))*B[(1, 2)] + ((-q^6*q1^2-q^6*q1*q2+q^4*q1*q2+q^4*q2^2+q^2*q1*q2+q^2*q2^2)/(-q^6*q1^2+q2^2))*B[(2, 0)] + ((-q^3*q1*q2-q^3*q2^2)/(q^6*q1^2-q2^2))*B[(2, -1)] + ((-q^5*q1^2-q^5*q1*q2+q^3*q1*q2+q^3*q2^2+q*q1*q2+q*q2^2)/(-q^6*q1^2+q2^2))*B[(2, 1)] + B[(2, 2)] + ((q^7*q1^2*q2+2*q^7*q1*q2^2+q^7*q2^3+q^5*q1^2*q2+2*q^5*q1*q2^2+q^5*q2^3)/(-q^8*q1^3-q^6*q1^2*q2+q^2*q1*q2^2+q2^3))*B[(0, -1)] + ((q^7*q1^3+q^7*q1^2*q2-q^7*q1*q2^2-q^7*q2^3-2*q^5*q1^2*q2-4*q^5*q1*q2^2-2*q^5*q2^3-2*q^3*q1^2*q2-4*q^3*q1*q2^2-2*q^3*q2^3)/(q^8*q1^3+q^6*q1^2*q2-q^2*q1*q2^2-q2^3))*B[(0, 1)] + ((q^6*q1^2+q^6*q1*q2-q^4*q1*q2-q^4*q2^2-q^2*q1*q2-q^2*q2^2)/(q^6*q1^2-q2^2))*B[(0, 2)]
        sage: E.recursion(2*omega[2])
        [0, 1, 0, 2, 1, 0, 2, 1, 0]

    Some tests that the `T` s are implemented properly by hand
    defining the `Y` s in terms of them::

        sage: T = E._T_Y
        sage: Ye1     = T.Tw((1,2,1,0), scalar = (-1/(q1*q2))^2)
        sage: Ye2     = T.Tw((2,1,0,1), signs = (1,1,1,-1), scalar = (-1/(q1*q2)))
        sage: Yalpha0 = T.Tw((0,1,2,1), signs = (-1,-1,-1,-1), scalar = q^-1*(-q1*q2)^2)
        sage: Yalpha1 = T.Tw((1,2,0,1,2,0), signs=(1,1,-1,1,-1,1), scalar = -1/(q1*q2))
        sage: Yalpha2 = T.Tw((2,1,0,1,2,1,0,1), signs = (1,1,1,-1,1,1,1,-1), scalar = (1/(q1*q2))^2)

        sage: Ye1(KL0.one())
        q1^2/q2^2*B[(0, 0)]
        sage: Ye2(KL0.one())
        ((-q1)/q2)*B[(0, 0)]

        sage: Yalpha0(KL0.one())
        q2^2/(q*q1^2)*B[(0, 0)]
        sage: Yalpha1(KL0.one())
        ((-q1)/q2)*B[(0, 0)]
        sage: Yalpha2(KL0.one())
        q1^2/q2^2*B[(0, 0)]

    Testing the `Y` s directly::

        sage: Y = E.Y()
        sage: Y.keys()
        Coroot lattice of the Root system of type ['BC', 2, 2]
        sage: alpha = Y.keys().simple_roots()
        sage: L(alpha[0])
        -2*e[0] + e['deltacheck']
        sage: L(alpha[1])
        e[0] - e[1]
        sage: L(alpha[2])
        e[1]

        sage: Y[alpha[0]].word
        (0, 1, 2, 1)
        sage: Y[alpha[0]].signs
        (-1, -1, -1, -1)
        sage: Y[alpha[0]].scalar # mind that Sage's q is the usual q^{1/2}
        q1^2*q2^2/q
        sage: Y[alpha[0]](KL0.one())
        q2^2/(q*q1^2)*B[(0, 0)]

        sage: Y[alpha[1]].word
        (1, 2, 0, 1, 2, 0)
        sage: Y[alpha[1]].signs
        (1, 1, -1, 1, -1, 1)
        sage: Y[alpha[1]].scalar
        1/(-q1*q2)

        sage: Y[alpha[2]].word    # Bogdan says it should be the square of that; do we need to take translation factors into account or not?
        (2, 1, 0, 1)
        sage: Y[alpha[2]].signs
        (1, 1, 1, -1)
        sage: Y[alpha[2]].scalar
        1/(-q1*q2)

    Checking the provided nonsymmetric Macdonald polynomial::

        sage: E10 = KL0.monomial(L0((1,0))) + KL0( q*(1-(-q1/q2)) / (1-q^2*(-q1/q2)^4) )
        sage: E10 == E[omega[1]]
        True
        sage: E.eigenvalues(E10)  # not checked
        [q*q1^2/q2^2, q2^3/(-q^2*q1^3), q1/(-q2)]

    Checking T0check::

        sage: T0check_on_basis = KL.T0_check_on_basis(q1,q2, convention="dominant")
        sage: T0check_on_basis.phi # note: this is in fact a0 phi
        (2, 0)
        sage: T0check_on_basis.v   # what to match it with?
        (1,)
        sage: T0check_on_basis.j   # what to match it with?
        2
        sage: T0check_on_basis(KL0.basis().keys().zero())
        ((-q1^2)/q2)*B[(1, 0)]

        sage: T0check = E._T[0]
        sage: T0check(KL0.one())
        ((-q1^2)/q2)*B[(1, 0)]


    Systematic tests of nonsymmetric Macdonald polynomials in type
    `A_1^{(1)}`, in the weight lattice. Each time, we specify the
    eigenvalues for the action of `Y_{\alpha_0}`, and `Y_{\alpha_1}`::

        sage: K = QQ['q','t'].fraction_field()
        sage: q,t = K.gens()
        sage: KL = RootSystem(["A",1,1]).weight_lattice(extended=True).algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q, t, -1)
        sage: omega = E.keys().fundamental_weights()

        sage: x = E[0*omega[1]]; x
        B[0]
        sage: E.eigenvalues(x)
        [1/(q*t), t]
        sage: x.is_one()
        True
        sage: x.parent()
        Algebra of the Weight lattice of the Root system of type ['A', 1]
        over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field

        sage: E[omega[1]]
        B[Lambda[1]]
        sage: E.eigenvalues(_)
        [t, 1/(q*t)]
        sage: E[2*omega[1]]
        ((-q*t+q)/(-q*t+1))*B[0] + B[2*Lambda[1]]
        sage: E.eigenvalues(_)
        [q*t, 1/(q^2*t)]
        sage: E[3*omega[1]]
        ((-q^2*t+q^2)/(-q^2*t+1))*B[-Lambda[1]] + ((-q^2*t+q^2-q*t+q)/(-q^2*t+1))*B[Lambda[1]] + B[3*Lambda[1]]
        sage: E.eigenvalues(_)
        [q^2*t, 1/(q^3*t)]
        sage: E[4*omega[1]]
        ((q^5*t^2-q^5*t+q^4*t^2-2*q^4*t+q^3*t^2+q^4-2*q^3*t+q^3-q^2*t+q^2)/(q^5*t^2-q^3*t-q^2*t+1))*B[0] + ((-q^3*t+q^3)/(-q^3*t+1))*B[-2*Lambda[1]] + ((-q^3*t+q^3-q^2*t+q^2-q*t+q)/(-q^3*t+1))*B[2*Lambda[1]] + B[4*Lambda[1]]
        sage: E.eigenvalues(_)
        [q^3*t, 1/(q^4*t)]
        sage: E[6*omega[1]]
        ((-q^12*t^3+q^12*t^2-q^11*t^3+2*q^11*t^2-2*q^10*t^3-q^11*t+4*q^10*t^2-2*q^9*t^3-2*q^10*t+5*q^9*t^2-2*q^8*t^3-4*q^9*t+6*q^8*t^2-q^7*t^3+q^9-5*q^8*t+5*q^7*t^2-q^6*t^3+q^8-6*q^7*t+4*q^6*t^2+2*q^7-5*q^6*t+2*q^5*t^2+2*q^6-4*q^5*t+q^4*t^2+2*q^5-2*q^4*t+q^4-q^3*t+q^3)/(-q^12*t^3+q^9*t^2+q^8*t^2+q^7*t^2-q^5*t-q^4*t-q^3*t+1))*B[0] + ((-q^5*t+q^5)/(-q^5*t+1))*B[-4*Lambda[1]] + ((q^9*t^2-q^9*t+q^8*t^2-2*q^8*t+q^7*t^2+q^8-2*q^7*t+q^6*t^2+q^7-2*q^6*t+q^5*t^2+q^6-2*q^5*t+q^5-q^4*t+q^4)/(q^9*t^2-q^5*t-q^4*t+1))*B[-2*Lambda[1]] + ((q^9*t^2-q^9*t+q^8*t^2-2*q^8*t+2*q^7*t^2+q^8-3*q^7*t+2*q^6*t^2+q^7-4*q^6*t+2*q^5*t^2+2*q^6-4*q^5*t+q^4*t^2+2*q^5-3*q^4*t+q^3*t^2+2*q^4-2*q^3*t+q^3-q^2*t+q^2)/(q^9*t^2-q^5*t-q^4*t+1))*B[2*Lambda[1]] + ((q^5*t-q^5+q^4*t-q^4+q^3*t-q^3+q^2*t-q^2+q*t-q)/(q^5*t-1))*B[4*Lambda[1]] + B[6*Lambda[1]]
        sage: E.eigenvalues(_)
        [q^5*t, 1/(q^6*t)]
        sage: E[-omega[1]]
        B[-Lambda[1]] + ((-t+1)/(-q*t+1))*B[Lambda[1]]
        sage: E.eigenvalues(_)
        [(-1)/(-q^2*t), q*t]

    As expected, `e^{-\omega}` is not an eigenvector::

        sage: E.eigenvalues(KL.classical().monomial(-omega[1]))
        Traceback (most recent call last):
        ...
        AssertionError

    We proceed by comparing against the examples from the appendix of
    [HHL06]_ in type `A_2^{(1)}`::

        sage: K = QQ['q','t'].fraction_field()
        sage: q,t = K.gens()
        sage: KL = RootSystem(["A",2,1]).ambient_space().algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q, t, -1)
        sage: L0 = E.keys()
        sage: omega = L0.fundamental_weights()
        sage: P = K['x0,x1,x2']
        sage: def EE(weight): return E[L0(weight)].expand(P.gens())

        sage: EE([0,0,0])
        1
        sage: EE([1,0,0])
        x0
        sage: EE([0,1,0])
        (t - 1)/(q*t^2 - 1)*x0 + x1
        sage: EE([0,0,1])
        (t - 1)/(q*t - 1)*x0 + (t - 1)/(q*t - 1)*x1 + x2
        sage: EE([1,1,0])
        x0*x1
        sage: EE([1,0,1])
        (t - 1)/(q*t^2 - 1)*x0*x1 + x0*x2
        sage: EE([0,1,1])
        (t - 1)/(q*t - 1)*x0*x1 + (t - 1)/(q*t - 1)*x0*x2 + x1*x2
        sage: EE([2,0,0])
        x0^2 + (q*t - q)/(q*t - 1)*x0*x1 + (q*t - q)/(q*t - 1)*x0*x2

        sage: EE([0,2,0])
        (t - 1)/(q^2*t^2 - 1)*x0^2 + (q^2*t^3 - q^2*t^2 + q*t^2 - 2*q*t + q - t + 1)/(q^3*t^3 - q^2*t^2 - q*t + 1)*x0*x1 + x1^2 + (q*t^2 - 2*q*t + q)/(q^3*t^3 - q^2*t^2 - q*t + 1)*x0*x2 + (q*t - q)/(q*t - 1)*x1*x2

    Systematic checks with Sage's implementation of [HHL06]_::

        sage: import sage.combinat.sf.ns_macdonald as NS
        sage: assert all(EE([x,y,z]) == NS.E([x,y,z]) for d in range(5) for x,y,z in IntegerVectors(d,3)) # long time (9s)

    We check that we get eigenvectors for generic `q_1`, `q_2`::

        sage: K = QQ['q,q1,q2'].fraction_field()
        sage: q,q1,q2 = K.gens()
        sage: KL = RootSystem(["A",2,1]).ambient_space().algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q, q1, q2)
        sage: L0 = E.keys()
        sage: omega = L0.fundamental_weights()
        sage: E[2*omega[2]]
        ((-q*q1-q*q2)/(-q*q1-q2))*B[(1, 2, 1)] + ((-q*q1-q*q2)/(-q*q1-q2))*B[(2, 1, 1)] + B[(2, 2, 0)]
        sage: for d in range(4):                                    # long time (9s)
        ....:     for weight in IntegerVectors(d,3).map(list).map(L0):
        ....:         eigenvalues = E.eigenvalues(E[L0(weight)])

    Some type `C` calculations::

        sage: K = QQ['q','t'].fraction_field()
        sage: q, t = K.gens()
        sage: KL = RootSystem(["C",2,1]).ambient_space().algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q, t, -1)
        sage: L0 = E.keys()
        sage: omega = L0.fundamental_weights()
        sage: E[0*omega[1]]
        B[(0, 0)]
        sage: E.eigenvalues(_)  # checked for i=0 with previous calculation
        [1/(q*t^3), t, t]
        sage: E[omega[1]]
        B[(1, 0)]
        sage: E.eigenvalues(_)  # not checked
        [t, 1/(q*t^3), t]

        sage: E[-omega[1]]          # consistent with before refactoring
        B[(-1, 0)] + ((-t+1)/(-q*t+1))*B[(1, 0)] + ((-t+1)/(-q*t+1))*B[(0, -1)] + ((t-1)/(q*t-1))*B[(0, 1)]
        sage: E.eigenvalues(_)  # not checked
        [(-1)/(-q^2*t^3), q*t, t]
        sage: E[-omega[1]+omega[2]] # consistent with before refactoring
        ((-t+1)/(-q*t^3+1))*B[(1, 0)] + B[(0, 1)]
        sage: E.eigenvalues(_)  # not checked
        [t, q*t^3, (-1)/(-q*t^2)]
        sage: E[omega[1]-omega[2]]  # consistent with before refactoring
        ((-t+1)/(-q*t^2+1))*B[(1, 0)] + B[(0, -1)] + ((-t+1)/(-q*t^2+1))*B[(0, 1)]
        sage: E.eigenvalues(_)  # not checked
        [1/(q^2*t^3), 1/(q*t), q*t^2]

        sage: E[-omega[2]]
        ((-q^2*t^4+q^2*t^3-q*t^3+2*q*t^2-q*t+t-1)/(-q^3*t^4+q^2*t^3+q*t-1))*B[(0, 0)] + B[(-1, -1)] + ((-t+1)/(-q*t+1))*B[(-1, 1)] + ((t-1)/(q*t-1))*B[(1, -1)] + ((-q*t^4+q*t^3+t-1)/(-q^3*t^4+q^2*t^3+q*t-1))*B[(1, 1)]
        sage: E.eigenvalues(_)  # not checked                       # long time (1s)
        [1/(q^3*t^3), t, q*t]
        sage: E[-omega[2]].map_coefficients(lambda c: c.subs(t=0))     # checking against crystals
        B[(0, 0)] + B[(-1, -1)] + B[(-1, 1)] + B[(1, -1)] + B[(1, 1)]

        sage: E[2*omega[2]]
        ((-q^6*t^7+q^6*t^6-q^5*t^6+2*q^5*t^5-q^4*t^5-q^5*t^3+3*q^4*t^4-3*q^4*t^3+q^3*t^4+q^4*t^2-2*q^3*t^2+q^3*t-q^2*t+q^2)/(-q^6*t^7+q^5*t^6+q^4*t^4+q^3*t^4-q^3*t^3-q^2*t^3-q*t+1))*B[(0, 0)] + ((-q^3*t^2+q^3*t)/(-q^3*t^3+1))*B[(-1, -1)] + ((-q^3*t^3+2*q^3*t^2-q^3*t)/(-q^4*t^4+q^3*t^3+q*t-1))*B[(-1, 1)] + ((-q^3*t^3+2*q^3*t^2-q^3*t)/(-q^4*t^4+q^3*t^3+q*t-1))*B[(1, -1)] + ((-q^4*t^4+q^4*t^3-q^3*t^3+2*q^3*t^2-q^2*t^3-q^3*t+2*q^2*t^2-q^2*t+q*t-q)/(-q^4*t^4+q^3*t^3+q*t-1))*B[(1, 1)] + ((q*t-q)/(q*t-1))*B[(2, 0)] + B[(2, 2)] + ((-q*t+q)/(-q*t+1))*B[(0, 2)]
        sage: E.eigenvalues(_)  # not checked
        [q^3*t^3, t, (-1)/(-q^2*t^2)]

    The following computations were calculated by hand::

        sage: KL0 = KL.classical()
        sage: E11 = KL0.sum_of_terms([[L0([1,1]), 1], [L0([0,0]), (-q*t^2 + q*t)/(1-q*t^3)]])
        sage: E11 == E[omega[2]]
        True
        sage: E.eigenvalues(E11)
        [q*t^3, t, (-1)/(-q*t^2)]

        sage: E1m1 = KL0.sum_of_terms([[L0([1,-1]), 1], [L0([1,1]), (1-t)/(1-q*t^2)], [L0([0,0]), q*t*(1-t)/(1-q*t^2)] ])
        sage: E1m1 == E[2*omega[1]-omega[2]]
        True
        sage: E.eigenvalues(E1m1)
        [1/(q*t), 1/(q^2*t^3), q*t^2]

    Now we present an example for a twisted affine root system. The
    results are eigenvectors::

        sage: K = QQ['q','t'].fraction_field()
        sage: q, t = K.gens()
        sage: KL = RootSystem("C2~*").ambient_space().algebra(K)
        sage: E = NonSymmetricMacdonaldPolynomials(KL,q, t, -1)
        sage: omega = E.keys().fundamental_weights()
        sage: E[0*omega[1]]
        B[(0, 0)]
        sage: E.eigenvalues(_)
        [1/(q*t^2), t, t]
        sage: E[omega[1]]
        ((-q*t+q)/(-q*t^2+1))*B[(0, 0)] + B[(1, 0)]
        sage: E.eigenvalues(_)
        [q*t^2, 1/(q^2*t^3), t]

        sage: E[-omega[1]]
        ((-q*t+q-t+1)/(-q^2*t+1))*B[(0, 0)] + B[(-1, 0)] + ((-t+1)/(-q^2*t+1))*B[(1, 0)] + ((-t+1)/(-q^2*t+1))*B[(0, -1)] + ((t-1)/(q^2*t-1))*B[(0, 1)]
        sage: E.eigenvalues(_)
        [(-1)/(-q^3*t^2), q^2*t, t]
        sage: E[-omega[1]+omega[2]]
        B[(-1/2, 1/2)] + ((-t+1)/(-q^2*t^3+1))*B[(1/2, -1/2)] + ((-q*t^3+q*t^2-t+1)/(-q^2*t^3+1))*B[(1/2, 1/2)]
        sage: E.eigenvalues(_)
        [(-1)/(-q^2*t^2), q^2*t^3, (-1)/(-q*t)]
        sage: E[omega[1]-omega[2]]
        B[(1/2, -1/2)] + ((-t+1)/(-q*t^2+1))*B[(1/2, 1/2)]
        sage: E.eigenvalues(_)
        [t, 1/(q^2*t^3), q*t^2]

    Type BC, comparison with calculations with Maple by Bogdan Ion::

        sage: K = QQ['q','t'].fraction_field()
        sage: q,t = K.gens()
        sage: def to_SR(x): return x.expand([SR.var('x%s'%i) for i in range(1,x.parent().basis().keys().dimension()+1)]).subs(q=SR.var('q'), t=SR.var('t'))
        sage: var('x1,x2,x3')
        (x1, x2, x3)

        sage: E = NonSymmetricMacdonaldPolynomials(["BC",2,2], q=q, q1=t^2,q2=-1)
        sage: omega=E.keys().fundamental_weights()
        sage: expected = (t-1)*(t+1)*(2+q^4+2*q^2-2*t^2-2*q^2*t^2-t^4*q^2-q^4*t^4+t^4-3*q^6*t^6-2*q^4*t^6+2*q^6*t^8+2*q^4*t^8+t^10*q^8)*q^4/((q^2*t^3-1)*(q^2*t^3+1)*(t*q-1)*(t*q+1)*(t^2*q^3+1)*(t^2*q^3-1))+(t-1)^2*(t+1)^2*(2*q^2+q^4+2+q^4*t^2)*q^3*x1/((t^2*q^3+1)*(t^2*q^3-1)*(t*q-1)*(t*q+1))+(t-1)^2*(t+1)^2*(q^2+1)*q^5/((t^2*q^3+1)*(t^2*q^3-1)*(t*q-1)*(t*q+1)*x1)+(t-1)^2*(t+1)^2*(q^2+1)*q^4*x2/((t^2*q^3+1)*(t^2*q^3-1)*(t*q-1)*(t*q+1)*x1)+(t-1)^2*(t+1)^2*(2*q^2+q^4+2+q^4*t^2)*q^3*x2/((t^2*q^3+1)*(t^2*q^3-1)*(t*q-1)*(t*q+1))+(t-1)^2*(t+1)^2*(q^2+1)*q^5/((t^2*q^3+1)*(t^2*q^3-1)*(t*q-1)*(t*q+1)*x2)+x1^2*x2^2+(t-1)*(t+1)*(-2*q^2-q^4-2+2*q^2*t^2+t^2+q^6*t^4+q^4*t^4)*q^2*x2*x1/((t^2*q^3+1)*(t^2*q^3-1)*(t*q-1)*(t*q+1))+(t-1)*(t+1)*(q^2+1+q^4*t^2)*q*x2^2*x1/((t^2*q^3-1)*(t^2*q^3+1))+(t-1)*(t+1)*q^3*x1^2/((t^2*q^3-1)*(t^2*q^3+1)*x2)+(t-1)*(t+1)*(q^2+1+q^4*t^2)*q*x2*x1^2/((t^2*q^3-1)*(t^2*q^3+1))+(t-1)*(t+1)*q^6/((t^2*q^3+1)*(t^2*q^3-1)*x1*x2)+(t-1)*(t+1)*(q^2+1+q^4*t^2)*q^2*x1^2/((t^2*q^3-1)*(t^2*q^3+1))+(t-1)*(t+1)*(q^2+1+q^4*t^2)*q^2*x2^2/((t^2*q^3-1)*(t^2*q^3+1))+(t-1)*(t+1)*q^3*x2^2/((t^2*q^3-1)*(t^2*q^3+1)*x1)+(t-1)^2*(t+1)^2*(q^2+1)*q^4*x1/((t^2*q^3+1)*(t^2*q^3-1)*(t*q-1)*(t*q+1)*x2)
        sage: to_SR(E[2*omega[2]]) - expected                       # long time (3.5s)
        0

        sage: E = NonSymmetricMacdonaldPolynomials(["BC",3,2], q=q, q1=t^2,q2=-1)
        sage: omega=E.keys().fundamental_weights()
        sage: mu = -3*omega[1] + 3*omega[2] - omega[3]; mu
        (-1, 2, -1)
        sage: expected = (t-1)^2*(t+1)^2*(3*q^2+q^4+1+t^2*q^4+q^2*t^2-3*t^4*q^2-5*t^6*q^4+2*t^8*q^4-4*t^8*q^6-q^8*t^10+2*t^10*q^6-2*q^8*t^12+t^14*q^8-t^14*q^10+q^10*t^16+q^8*t^16+q^10*t^18+t^18*q^12)*x2*x1/((q^3*t^5+1)*(q^3*t^5-1)*(t*q-1)*(t*q+1)*(t^3*q^2+1)*(t^3*q^2-1)*(t^2*q-1)*(t^2*q+1))+(t-1)^2*(t+1)^2*(q^2*t^6+2*t^6*q^4-q^4*t^4+t^4*q^2-q^2*t^2+t^2-2-q^2)*q^2*x1/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x2)+(t-1)^2*(t+1)^2*(-q^2-1+t^4*q^2-q^4*t^4+2*t^6*q^4)*x1^2/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1))+(t+1)*(t-1)*x2^2*x3/((t*q-1)*(t*q+1)*x1)+(t-1)^2*(t+1)^2*(3*q^2+q^4+2+t^2*q^4+2*q^2*t^2-4*t^4*q^2+q^4*t^4-6*t^6*q^4+t^8*q^4-4*t^8*q^6-q^8*t^10+t^10*q^6-3*q^8*t^12-2*t^14*q^10+2*t^14*q^8+2*q^10*t^16+q^8*t^16+t^18*q^12+2*q^10*t^18)*q*x2/((q^3*t^5+1)*(q^3*t^5-1)*(t*q-1)*(t*q+1)*(t^3*q^2+1)*(t^3*q^2-1)*(t^2*q-1)*(t^2*q+1))+(t-1)^2*(t+1)^2*(1+q^4+2*q^2+t^2*q^4-3*t^4*q^2+q^2*t^6-5*t^6*q^4+3*t^8*q^4-4*t^8*q^6+2*t^10*q^6-q^8*t^12-t^14*q^10+t^14*q^8+q^10*t^16+t^18*q^12)*x3*x1/((q^3*t^5+1)*(q^3*t^5-1)*(t*q-1)*(t*q+1)*(t^3*q^2+1)*(t^3*q^2-1)*(t^2*q-1)*(t^2*q+1))+(t-1)^2*(t+1)^2*(2*q^2+1+q^4+t^2*q^4-t^2+q^2*t^2-4*t^4*q^2+q^4*t^4+q^2*t^6-5*t^6*q^4+3*t^8*q^4-4*t^8*q^6+2*t^10*q^6+q^6*t^12-2*q^8*t^12-2*t^14*q^10+2*t^14*q^8+q^10*t^16+t^18*q^12)*q*x3/((q^3*t^5+1)*(q^3*t^5-1)*(t*q-1)*(t*q+1)*(t^3*q^2+1)*(t^3*q^2-1)*(t^2*q-1)*(t^2*q+1))+(t-1)^2*(t+1)^2*(1+t^2+t^4*q^2)*q*x3*x2^2/((t*q-1)*(t*q+1)*(t^3*q^2+1)*(t^3*q^2-1))+(t-1)^2*(t+1)^2*(-q^2-2-q^2*t^2+t^4-q^4*t^4-t^4*q^2+3*q^2*t^6-t^6*q^4-t^8*q^6+t^8*q^4+t^10*q^4+2*q^6*t^12-q^8*t^12+t^14*q^8)*q*x3*x2*x1/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1))+(t-1)*(t+1)*x1^2/((q^3*t^5-1)*(q^3*t^5+1)*x3*x2)+(t-1)*(t+1)*(-q^2-1+t^4*q^2-q^4*t^4+2*t^6*q^4)*x2^2/((t*q-1)*(t*q+1)*(t^3*q^2+1)*(t^3*q^2-1))+(t-1)*(t+1)*(t^3*q-1)*(t^3*q+1)*x3*x2^2*x1/((t*q-1)*(t*q+1)*(t^3*q^2+1)*(t^3*q^2-1))+(t-1)^2*(t+1)^2*(q^2+1)*q*x1/((t*q+1)*(t*q-1)*(q^3*t^5+1)*(q^3*t^5-1)*x3*x2)+(t-1)^2*(t+1)^2*(t^3*q-1)*(t^3*q+1)*x3*x2*x1^2/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1))+(t-1)^2*(t+1)^2*q^3*x3/((t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x1*x2)+(t-1)*(t+1)*(-1-q^2+q^2*t^2+t^10*q^6)*q*x2/((t*q+1)*(t*q-1)*(q^3*t^5+1)*(q^3*t^5-1)*x3*x1)+x2^2/(x1*x3)+(t-1)*(t+1)*q*x2^2/((t*q-1)*(t*q+1)*x3)+(t-1)^3*(t+1)^3*(1+t^2+t^4*q^2)*q*x2*x1^2/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1))+(t-1)^2*(t+1)^2*q*x1^2/((t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x3)+(t-1)^2*(t+1)^2*(q^2*t^6+2*t^6*q^4-q^4*t^4+t^4*q^2-q^2*t^2+t^2-2-q^2)*q^3/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x2)+(t-1)*(t+1)*(q^2+2-t^2+q^4*t^4-t^4*q^2-3*t^6*q^4+t^8*q^4-2*t^10*q^6-q^8*t^12+q^6*t^12+q^8*t^16+q^10*t^16)*q^2*x2/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x1)+(t-1)^2*(t+1)^2*(q^2+1)*q^2/((t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x3*x2)+(t-1)*(t+1)*(1+q^4+2*q^2-2*q^2*t^2+t^4*q^6-q^4*t^4-3*q^6*t^6-t^6*q^4+2*t^8*q^6-t^10*q^6-q^8*t^10-t^14*q^10+t^14*q^8+2*q^10*t^16)*x2/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x3)+(t-1)^2*(t+1)^2*(-q^2-2-q^2*t^2-q^4*t^4+2*t^6*q^4+t^10*q^6+q^8*t^12+t^14*q^8)*q^3/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x1)+(t-1)^2*(t+1)^2*(-1-q^2-q^2*t^2+t^2+t^4*q^2-q^4*t^4+2*t^6*q^4)*q^2*x3/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x2)+(t-1)*(t+1)*q*x2^2/((t*q-1)*(t*q+1)*x1)+(t-1)^2*(t+1)^2*(1+t^2+t^4*q^2)*q*x2^2*x1/((t*q-1)*(t*q+1)*(t^3*q^2+1)*(t^3*q^2-1))+(t-1)^2*(t+1)^2*q*x1^2/((t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x2)+(t-1)^2*(t+1)^2*(-1-q^4-2*q^2-t^2*q^4-q^2*t^2+t^4*q^2-t^4*q^6-2*q^4*t^4+3*t^6*q^4-q^6*t^6-t^8*q^8+t^8*q^6+2*t^10*q^6-q^10*t^12+3*q^8*t^12+2*t^14*q^10)*x3*x2/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1))+(t-1)*(t+1)*(q^2+1-t^2+q^4*t^4-t^4*q^2+q^2*t^6-3*t^6*q^4+t^8*q^4-t^10*q^6+q^6*t^12-q^8*t^12+q^10*t^16)*q^2*x3/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x1)+(t-1)*(t+1)*(-1-q^2+q^2*t^2+t^10*q^6)*q^2/((t*q-1)*(t*q+1)*(q^3*t^5+1)*(q^3*t^5-1)*x1*x3)+(t-1)*(t+1)*(1+q^4+2*q^2-3*q^2*t^2+t^4*q^6-q^4*t^4-3*q^6*t^6-t^6*q^4+t^8*q^4+2*t^8*q^6-t^10*q^6+t^14*q^8-t^14*q^10+q^10*t^16)*x1/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x3)+(t-1)^2*(t+1)^2*(3*q^2+q^4+2+q^2*t^2-t^2+t^2*q^4-6*t^4*q^2+q^4*t^4-7*t^6*q^4+q^2*t^6+3*t^8*q^4-4*t^8*q^6+t^10*q^4+3*t^10*q^6-q^8*t^12-t^14*q^10+t^14*q^8+q^8*t^16+q^10*t^18)*q*x1/((q^3*t^5+1)*(q^3*t^5-1)*(t*q-1)*(t*q+1)*(t^3*q^2+1)*(t^3*q^2-1)*(t^2*q-1)*(t^2*q+1))+(t-1)^2*(t+1)^2*(-q^2-2-q^2*t^2-q^4*t^4+2*t^6*q^4+t^10*q^6+q^6*t^12+t^14*q^8)*q*x2*x1/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x3)+(t+1)*(t-1)*x2^2*x1/((t*q-1)*(t*q+1)*x3)+(t-1)^3*(t+1)^3*(1+t^2+t^4*q^2)*q*x3*x1^2/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1))+(t-1)*(t+1)*q^3/((q^3*t^5+1)*(q^3*t^5-1)*x1*x2*x3)+(t-1)^2*(t+1)^2*(3+3*q^2+q^4+2*q^2*t^2-t^2+t^2*q^4-6*t^4*q^2+q^4*t^4-8*t^6*q^4+q^2*t^6+2*t^8*q^4-4*t^8*q^6+t^10*q^4+2*t^10*q^6-2*q^8*t^12-t^14*q^10+t^14*q^8+q^8*t^16+q^10*t^16+2*q^10*t^18)*q^2/((q^3*t^5+1)*(q^3*t^5-1)*(t*q-1)*(t*q+1)*(t^3*q^2+1)*(t^3*q^2-1)*(t^2*q-1)*(t^2*q+1))+(t-1)^2*(t+1)^2*(-q^4-2*q^2-1-t^2*q^4-t^4*q^6+2*q^6*t^6+t^6*q^4+t^10*q^6+q^8*t^12+t^14*q^10)*q/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x3)+(t-1)^2*(t+1)^2*(-1-q^2-q^2*t^2+t^2+t^4*q^2-q^4*t^4+2*t^6*q^4)*q*x3*x1/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x2)+(t-1)^2*(t+1)^2*x2*x1^2/((t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x3)+(t-1)^2*(t+1)^2*x3*x1^2/((t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x2)+(t-1)^2*(t+1)^2*q^4/((t*q+1)*(t*q-1)*(q^3*t^5+1)*(q^3*t^5-1)*x1*x2)+(t-1)^2*(t+1)^2*(-q^2-1-q^2*t^2-q^4*t^4+t^6*q^4+t^10*q^6+q^8*t^12+t^14*q^10)*q*x3*x2/((t^3*q^2-1)*(t^3*q^2+1)*(t*q+1)*(t*q-1)*(q^3*t^5-1)*(q^3*t^5+1)*x1)
        sage: to_SR(E[mu]) - expected                               # long time (20s)
        0

        sage: E = NonSymmetricMacdonaldPolynomials(["BC",1,2], q=q, q1=t^2,q2=-1)
        sage: omega=E.keys().fundamental_weights()
        sage: mu = -4*omega[1]; mu
        (-4)
        sage: expected = (t-1)*(t+1)*(-1+q^2*t^2-q^2-3*q^10-7*q^26*t^8+5*t^2*q^6-q^16-3*q^4+4*t^10*q^30-4*t^6*q^22-10*q^20*t^6+2*q^32*t^10-3*q^6-4*q^8+q^34*t^10-4*t^8*q^24-2*q^12-q^14+2*q^22*t^10+4*q^26*t^10+4*q^28*t^10+t^6*q^30-2*q^32*t^8-2*t^8*q^22+2*q^24*t^10-q^20*t^2-2*t^6*q^12+t^8*q^14+2*t^4*q^24-4*t^8*q^30+2*t^8*q^20-9*t^6*q^16+3*q^26*t^6+q^28*t^6+3*t^2*q^4+2*q^18*t^8-6*t^6*q^14+4*t^4*q^22-2*q^24*t^6+3*t^2*q^12+7*t^4*q^20-t^2*q^16+11*q^18*t^4-2*t^2*q^18+9*q^16*t^4-t^4*q^6+6*q^8*t^2+5*q^10*t^2-6*q^28*t^8+q^12*t^4+8*t^4*q^14-10*t^6*q^18-q^4*t^4+q^16*t^8-2*t^4*q^8)/((t*q^4-1)*(t*q^4+1)*(q^7*t^2-1)*(q^7*t^2+1)*(t*q^3-1)*(t*q^3+1)*(q^5*t^2+1)*(q^5*t^2-1))+(q^2+1)*(q^4+1)*(t-1)*(t+1)*(-1+q^2*t^2-q^2+t^2*q^6-q^4+t^6*q^22+3*q^10*t^4+t^2-q^8-2*t^8*q^24+q^22*t^10+q^26*t^10-2*t^8*q^22+q^24*t^10-4*t^6*q^12-2*t^8*q^20-3*t^6*q^16+2*t^2*q^4-t^6*q^10-2*t^6*q^14+t^8*q^12-t^2*q^12+2*q^16*t^4+q^8*t^2-q^10*t^2+3*q^12*t^4+2*t^4*q^14+t^6*q^18-2*q^4*t^4+q^16*t^8+q^20*t^10)*q*x1/((t*q^4-1)*(t*q^4+1)*(q^7*t^2-1)*(q^7*t^2+1)*(t*q^3-1)*(t*q^3+1)*(q^5*t^2+1)*(q^5*t^2-1))+(q^2+1)*(q^4+1)*(t-1)*(t+1)*(1+q^8+q^4+q^2-q^8*t^2-2*t^2*q^4-t^2*q^6+t^2*q^12-t^2+t^4*q^6-2*q^16*t^4-t^4*q^14-2*q^12*t^4+t^6*q^12+t^6*q^16+t^6*q^18+t^6*q^14)*q/((t*q^4-1)*(t*q^4+1)*(q^7*t^2-1)*(q^7*t^2+1)*(t*q^3-1)*(t*q^3+1)*x1)+(t-1)*(t+1)*(-1-q^2-q^6-q^4-q^8+t^2*q^4-t^2*q^14+t^2*q^6-q^10*t^2+q^8*t^2-t^2*q^12+q^12*t^4+q^10*t^4+q^16*t^4+2*t^4*q^14)*(q^4+1)/((q^7*t^2+1)*(q^7*t^2-1)*(t*q^4-1)*(t*q^4+1)*x1^2)+(t-1)*(t+1)*(q^4+1)*(q^2+1)*q/((t*q^4-1)*(t*q^4+1)*x1^3)+(q^4+1)*(t-1)*(t+1)*(1+q^6+q^8+q^2+q^4-q^2*t^2-3*t^2*q^4+q^10*t^2+t^2*q^12-2*t^2*q^6-q^8*t^2-2*q^16*t^4+q^4*t^4+t^4*q^6-q^10*t^4-2*q^12*t^4-2*t^4*q^14+t^6*q^12+t^6*q^18+2*t^6*q^16+t^6*q^14)*x1^2/((t*q^4-1)*(t*q^4+1)*(q^7*t^2-1)*(q^7*t^2+1)*(t*q^3-1)*(t*q^3+1))+(t-1)*(t+1)*(-1-t^2*q^6+t^2+t^4*q^8)*(q^4+1)*(q^2+1)*q*x1^3/((q^7*t^2+1)*(q^7*t^2-1)*(t*q^4-1)*(t*q^4+1))+1/x1^4+(t-1)*(t+1)*x1^4/((t*q^4-1)*(t*q^4+1))
        sage: to_SR(E[mu]) - expected
        0

    Type `BC` dual, comparison with hand calculations by Bogdan Ion::

        sage: K = QQ['q,q1,q2'].fraction_field()
        sage: q,q1,q2 = K.gens()
        sage: ct = CartanType(["BC",2,2]).dual()
        sage: E = NonSymmetricMacdonaldPolynomials(ct, q=q, q1=q1, q2=q2)
        sage: KL = E.domain(); KL
        Algebra of the Ambient space of the Root system of type ['B', 2]
        over Fraction Field of Multivariate Polynomial Ring in q, q1, q2 over Rational Field
        sage: alpha = E.keys().simple_roots(); alpha
        Finite family {1: (1, -1), 2: (0, 1)}
        sage: omega=E.keys().fundamental_weights(); omega
        Finite family {1: (1, 0), 2: (1/2, 1/2)}
        sage: epsilon = E.keys().basis(); epsilon
        Finite family {0: (1, 0), 1: (0, 1)}

    Note: Sage's `q` is the usual `q^2`::

        sage: E.L().null_root()
        e['delta']
        sage: E.L().null_coroot()
        2*e['deltacheck']

    Some eigenvectors::

        sage: E[0*omega[1]]
        B[(0, 0)]
        sage: E[omega[1]]
        ((-q^2*q1^3*q2-q^2*q1^2*q2^2)/(q^2*q1^4-q2^4))*B[(0, 0)] + B[(1, 0)]
        sage: Eomega1 = KL.one() * (q^2*(-q1/q2)^2*(1-(-q1/q2))) / (1-q^2*(-q1/q2)^4) + KL.monomial(omega[1])
        sage: E[omega[1]] == Eomega1
        True

    Checking the `Y` s::

        sage: Y = E.Y()
        sage: alphacheck = Y.keys().simple_roots()
        sage: Y0 = Y[alphacheck[0]]
        sage: Y1 = Y[alphacheck[1]]
        sage: Y2 = Y[alphacheck[2]]

        sage: Y0.word, Y0.signs, Y0.scalar
        ((0, 1, 2, 1, 0, 1, 2, 1), (-1, -1, -1, -1, -1, -1, -1, -1), q1^4*q2^4/q^2)
        sage: Y1.word, Y1.signs, Y1.scalar
        ((1, 2, 0, 1, 2, 0), (1, 1, -1, 1, -1, 1), 1/(-q1*q2))
        sage: Y2.word, Y2.signs, Y2.scalar
        ((2, 1, 0, 1), (1, 1, 1, -1), 1/(-q1*q2))

        sage: E.eigenvalues(0*omega[1])
        [q2^4/(q^2*q1^4), q1/(-q2), q1/(-q2)]

    Checking the `T` and `T^{-1}` s::

        sage: T = E._T_Y
        sage: Tinv0 = T.Tw_inverse([0])
        sage: Tinv1 = T.Tw_inverse([1])
        sage: Tinv2 = T.Tw_inverse([2])

        sage: for x in [0*epsilon[0], -epsilon[0], -epsilon[1], epsilon[0], epsilon[1]]:
        ....:     x = KL.monomial(x)
        ....:     assert Tinv0(T[0](x)) == x and T[0](Tinv0(x)) == x
        ....:     assert Tinv1(T[1](x)) == x and T[1](Tinv1(x)) == x
        ....:     assert Tinv2(T[2](x)) == x and T[2](Tinv2(x)) == x

        sage: start = E[omega[1]]; start
        ((-q^2*q1^3*q2-q^2*q1^2*q2^2)/(q^2*q1^4-q2^4))*B[(0, 0)] + B[(1, 0)]
        sage: Tinv1(Tinv2(Tinv1(Tinv0(Tinv1(Tinv2(Tinv1(Tinv0(start)))))))) * (q1*q2)^4/q^2 == Y0(start)
        True
        sage: Y0(start) == q^2*q1^4/q2^4 * start
        True

    Checking the relation between the `Y` s::

        sage: q^2 * Y0(Y1(Y1(Y2(Y2(start))))) == start
        True
        sage: for x in [0*epsilon[0], -epsilon[0], -epsilon[1], epsilon[0], epsilon[1]]:
        ....:     x = KL.monomial(x)
        ....:     assert q^2 * Y0(Y1(Y1(Y2(Y2(start))))) == start

    """

    @staticmethod
    def __classcall__(cls, KL, q='q', q1='q1', q2='q2', normalized=True):
        r"""
        EXAMPLES::

            sage: NonSymmetricMacdonaldPolynomials(["B", 2, 1])
            The family of the Macdonald polynomials of type ['B', 2, 1] with parameters q, q1, q2
        """
        from sage.combinat.root_system.cartan_type import CartanType
        K = None
        #if KL in Algebras:
        if isinstance(KL, CombinatorialFreeModule): # temporary work around C3 issue ...
            K = KL.base_ring()
        else:
            if q == 'q':
                from sage.rings.rational_field import QQ
                K = QQ['q','q1','q2'].fraction_field()
            else:
                K = q.parent()
            KL = CartanType(KL).root_system().ambient_space().algebra(K)
        q = K(q)
        q1 = K(q1)
        q2 = K(q2)
        return super(NonSymmetricMacdonaldPolynomials, cls).__classcall__(cls, KL, q, q1, q2, normalized)

    def __init__(self, KL, q, q1, q2, normalized):
        r"""
        Initializes the nonsymmetric Macdonald polynomial class.

        INPUT:

        - ``KL`` -- algebra over weight space
        - ``q``, ``q1``, ``q2`` -- parameters
        - ``normalized`` -- a boolean (default: True)
           whether to normalize the result to have leading coefficient 1

        EXAMPLES::

            sage: K = QQ['q,q1,q2'].fraction_field()
            sage: q, q1, q2 = K.gens()
            sage: KL = RootSystem(["A",1,1]).weight_space(extended = True).algebra(K)
            sage: NonSymmetricMacdonaldPolynomials(KL,q, q1, q2)
            The family of the Macdonald polynomials of type ['A', 1, 1] with parameters q, q1, q2

            sage: KL = RootSystem(["A",1,1]).ambient_space().algebra(K)
            sage: NonSymmetricMacdonaldPolynomials(KL,q, q1, q2)
            The family of the Macdonald polynomials of type ['A', 1, 1] with parameters q, q1, q2

            sage: KL = RootSystem(["A",1,1]).weight_space().algebra(K)
            sage: NonSymmetricMacdonaldPolynomials(KL,q, q1, q2)
            Traceback (most recent call last):
            ...
            AssertionError: The weight lattice needs to be extended!

        """
        # TODO: check all the choices!
        self._KL = KL
        self._L = KL.basis().keys()
        assert self._L.is_extended(), "The weight lattice needs to be extended!"
        self._q = q
        self._q1 = q1
        self._q2 = q2
        assert self.L_prime().classical() is self.L().classical()
        T   = KL.twisted_demazure_lusztig_operators     (   q1, q2, convention="dominant")
        T_Y = KL.demazure_lusztig_operators_on_classical(q, q1, q2, convention="dominant")
        CherednikOperatorsEigenvectors.__init__(self, T, T_Y, normalized = normalized)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: NonSymmetricMacdonaldPolynomials(["B", 2, 1])
            The family of the Macdonald polynomials of type ['B', 2, 1] with parameters q, q1, q2
        """
        return "The family of the Macdonald polynomials of type %s with parameters %s, %s, %s"%(self.cartan_type(),self._q, self._q1, self._q2)

    # This is redundant with the cartan_type method of
    # CherednikOperatorsEigenvectors, but we need it very early in the
    # initialization, before self._T_Y is set ...
    @cached_method
    def cartan_type(self):
        r"""
        Return Cartan type of ``self``.

        EXAMPLES::

            sage: NonSymmetricMacdonaldPolynomials(["B", 2, 1]).cartan_type()
            ['B', 2, 1]
        """
        return self._L.cartan_type()

    def L(self):
        r"""
        Return the affinization of the classical weight space.

        EXAMPLES::

            sage: NonSymmetricMacdonaldPolynomials(["B", 2, 1]).L()
            Ambient space of the Root system of type ['B', 2, 1]
        """
        return self._L

    @cached_method
    def L_check(self):
        r"""
        Return the other affinization of the classical weight space.

        .. TODO:: should this just return `L` in the simply laced case?

        EXAMPLES::

            sage: NonSymmetricMacdonaldPolynomials(["B", 2, 1]).L_check()
            Coambient space of the Root system of type ['C', 2, 1]
            sage: NonSymmetricMacdonaldPolynomials(["B", 2, 1]).L_check().classical()
            Ambient space of the Root system of type ['B', 2]
        """
        from sage.combinat.root_system.weight_space import WeightSpace
        from sage.combinat.root_system.type_affine import AmbientSpace
        L = self.L()
        other_affine_root_system = self.cartan_type().classical().dual().affine().root_system()
        if isinstance(L, WeightSpace): # TODO: make a nicer test
            return other_affine_root_system.coweight_space(L.base_ring(), extended=True)
        else:
            assert isinstance(L, AmbientSpace)
            return other_affine_root_system.coambient_space(L.base_ring())

    @cached_method
    def L_prime(self):
        r"""
        The affine space where classical weights are lifted for the recursion.

        Also the parent of `\rho'`.

        EXAMPLES:

        In the twisted case, this is the affinization of the classical
        ambient space::

            sage: NonSymmetricMacdonaldPolynomials("B2~*").L()
            Ambient space of the Root system of type ['B', 2, 1]^*
            sage: NonSymmetricMacdonaldPolynomials("B2~*").L().classical()
            Ambient space of the Root system of type ['C', 2]

            sage: NonSymmetricMacdonaldPolynomials("B2~*").L_prime()
            Ambient space of the Root system of type ['B', 2, 1]^*
            sage: NonSymmetricMacdonaldPolynomials("B2~*").L_prime().classical()
            Ambient space of the Root system of type ['C', 2]

        In the untwisted case, this is the other affinization of the
        classical ambient space::

            sage: NonSymmetricMacdonaldPolynomials("B2~").L()
            Ambient space of the Root system of type ['B', 2, 1]
            sage: NonSymmetricMacdonaldPolynomials("B2~").L().classical()
            Ambient space of the Root system of type ['B', 2]

            sage: NonSymmetricMacdonaldPolynomials("B2~").L_prime()
            Coambient space of the Root system of type ['C', 2, 1]
            sage: NonSymmetricMacdonaldPolynomials("B2~").L_prime().classical()
            Ambient space of the Root system of type ['B', 2]

        For simply laced, the two affinizations coincide::

            sage: NonSymmetricMacdonaldPolynomials("A2~").L()
            Ambient space of the Root system of type ['A', 2, 1]
            sage: NonSymmetricMacdonaldPolynomials("A2~").L().classical()
            Ambient space of the Root system of type ['A', 2]

            sage: NonSymmetricMacdonaldPolynomials("A2~").L_prime()
            Coambient space of the Root system of type ['A', 2, 1]
            sage: NonSymmetricMacdonaldPolynomials("A2~").L_prime().classical()
            Ambient space of the Root system of type ['A', 2]

        .. NOTE:: do we want the coambient space of type `A_2^{(1)}` instead?

        For type BC::

            sage: NonSymmetricMacdonaldPolynomials(["BC",3,2]).L_prime()
            Ambient space of the Root system of type ['BC', 3, 2]
        """
        ct = self.cartan_type()
        if ct.is_untwisted_affine():
            return self.L_check()
        else:
            return self.L()

    @cached_method
    def L0(self):
        r"""
        Return the space indexing the monomials of the nonsymmetric Macdonald polynomials.

        EXAMPLES::

            sage: NonSymmetricMacdonaldPolynomials("B2~").L0()
            Ambient space of the Root system of type ['B', 2]
            sage: NonSymmetricMacdonaldPolynomials("B2~*").L0()
            Ambient space of the Root system of type ['C', 2]
        """
        return self.L().classical()

    @cached_method
    def KL0(self):
        r"""
        Return the group algebra where the nonsymmetric Macdonald polynomials live.

        EXAMPLES::

            sage: NonSymmetricMacdonaldPolynomials("B2~").KL0()
            Algebra of the Ambient space of the Root system of type ['B', 2]
            over Fraction Field of Multivariate Polynomial Ring in q, q1, q2 over Rational Field
            sage: NonSymmetricMacdonaldPolynomials("B2~*").KL0()
            Algebra of the Ambient space of the Root system of type ['C', 2]
            over Fraction Field of Multivariate Polynomial Ring in q, q1, q2 over Rational Field

        """
        return self._KL.classical()

    @lazy_attribute
    def Q_to_Qcheck(self):
        r"""
        The reindexing of the index set of the Y's by the coroot lattice.

        EXAMPLES::

            sage: E = NonSymmetricMacdonaldPolynomials("C2~")
            sage: alphacheck = E.Y().keys().simple_roots()
            sage: E.Q_to_Qcheck(alphacheck[0])
            alphacheck[0] - alphacheck[2]
            sage: E.Q_to_Qcheck(alphacheck[1])
            alphacheck[1]
            sage: E.Q_to_Qcheck(alphacheck[2])
            alphacheck[2]

            sage: x = alphacheck[1] + 2*alphacheck[2]
            sage: x.parent()
            Root lattice of the Root system of type ['B', 2, 1]
            sage: E.Q_to_Qcheck(x)
            alphacheck[1] + 2*alphacheck[2]
            sage: _.parent()
            Coroot lattice of the Root system of type ['C', 2, 1]
        """
        #assert self.cartan_type().is_untwisted_affine()
        Qcheck = self._T_Y.Y().keys()
        Q = Qcheck.cartan_type().other_affinization().root_system().root_lattice()
        assert Q.classical() is Qcheck.classical()
        return Q.module_morphism(Qcheck.simple_roots_tilde().__getitem__, codomain=Qcheck)

    def Y(self):
        r"""
        Return the family of `Y` operators whose eigenvectors are the nonsymmetric Macdonald polynomials.

        EXAMPLES::

            sage: NonSymmetricMacdonaldPolynomials("C2~").Y()
            Lazy family (<lambda>(i))_{i in Root lattice of the Root system of type ['B', 2, 1]}
            sage: _.keys().classical()
            Root lattice of the Root system of type ['B', 2]

            sage: NonSymmetricMacdonaldPolynomials("C2~*").Y()
            Lazy family (<...Y_lambdacheck...>(i))_{i in Coroot lattice of the Root system of type ['C', 2, 1]^*}
            sage: _.keys().classical()
            Root lattice of the Root system of type ['C', 2]

            sage: NonSymmetricMacdonaldPolynomials(["BC", 3, 2]).Y()
            Lazy family (<...Y_lambdacheck...>(i))_{i in Coroot lattice of the Root system of type ['BC', 3, 2]}
            sage: _.keys().classical()
            Root lattice of the Root system of type ['B', 3]
        """
        from sage.sets.family import Family
        Y = self._T_Y.Y()
        ct = self.cartan_type()
        # TODO: improve test
        if ct.dual().is_untwisted_affine() or ct.type() == "BC":
            return Y
        Q = self.Q_to_Qcheck.domain()
        return Family(Q, lambda lambdacheck: Y[self.Q_to_Qcheck(lambdacheck)])

    def affine_lift(self, mu):
        r"""
        Return the affinization of `\mu` in `L'`.

        INPUT:

        - ``mu`` -- a classical weight `\mu`

        .. SEEALSO::

            - :meth:`.hecke_algebra_representation.CherednikOperatorsEigenvectors.affine_lift`
            - :meth:`affine_retract`
            - :meth:`L_prime`

        EXAMPLES:

        In the untwisted case, this is the other affinization at level 1::

            sage: E = NonSymmetricMacdonaldPolynomials("B2~")
            sage: L0 = E.keys(); L0
            Ambient space of the Root system of type ['B', 2]
            sage: omega = L0.fundamental_weights()
            sage: E.affine_lift(omega[1])
            e[0] + e['deltacheck']
            sage: E.affine_lift(omega[1]).parent()
            Coambient space of the Root system of type ['C', 2, 1]

        In the twisted case, this is the usual affinization at level 1::

            sage: E = NonSymmetricMacdonaldPolynomials("B2~*")
            sage: L0 = E.keys(); L0
            Ambient space of the Root system of type ['C', 2]
            sage: omega = L0.fundamental_weights()
            sage: E.affine_lift(omega[1])
            e[0] + e['deltacheck']
            sage: E.affine_lift(omega[1]).parent()
            Ambient space of the Root system of type ['B', 2, 1]^*
        """
        return self.L_prime().embed_at_level(mu, 1)

    def twist(self, mu, i):
        r"""
        Act by `s_i` on the affine weight `\mu`.

        This calls ``simple_reflection``; which is semantically the
        same as the default implementation.

        EXAMPLES::

            sage: W = WeylGroup(["B",3])
            sage: W.element_class._repr_=lambda x: "".join(str(i) for i in x.reduced_word())
            sage: K = QQ['q1,q2']
            sage: q1, q2 = K.gens()
            sage: KW = W.algebra(K)
            sage: T = KW.demazure_lusztig_operators(q1, q2, affine=True)
            sage: E = T.Y_eigenvectors()
            sage: w = W.an_element(); w
            123
            sage: E.twist(w,1)
            1231
        """
        return mu.simple_reflection(i)

    def affine_retract(self, mu):
        r"""
        Retract the affine weight `\mu` into a classical weight.

        INPUT:

        - ``mu`` -- an affine weight `\mu` in `L'`

        .. SEEALSO::

            - :meth:`.hecke_algebra_representation.HeckeAlgebraRepresentation.affine_retract`
            - :meth:`affine_lift`
            - :meth:`L_prime`

        EXAMPLES::

            sage: E = NonSymmetricMacdonaldPolynomials("B2~")
            sage: L0 = E.keys(); L0
            Ambient space of the Root system of type ['B', 2]
            sage: omega = L0.fundamental_weights()
            sage: E.affine_lift(omega[1])
            e[0] + e['deltacheck']
            sage: E.affine_retract(E.affine_lift(omega[1]))
            (1, 0)
        """
        assert mu in self.L_prime()
        return self.L0()(mu)

    def __getitem__(self, mu):
        r"""
        Return the nonsymmetric Macdonald polynomial `E_\mu`.

        INPUT:

        - ``mu`` -- a weight `\mu` that lifts to a level 0 element of the affine weight lattice

        This methods simply checks the weight and calls
        :meth:`.hecke_algebra_representation.CherednikOperatorsEigenvectors.__getitem__`.

        .. NOTE::

            Any element of the finite weight lattice lifts to a level
            0 element of the affine weight lattice.
            Exception: `\omega_n` in type `BC_n` dual.

        EXAMPLES::

            sage: ct = CartanType(["BC",2,2]).dual()
            sage: E = NonSymmetricMacdonaldPolynomials(ct)
            sage: omega = E.keys().fundamental_weights()
            sage: omega[2]
            (1/2, 1/2)
            sage: E[omega[2]]
            Traceback (most recent call last):
            ...
            ValueError: 1/2*e[0] + 1/2*e[1] does not lift to a level 0 element of the affine weight lattice
            sage: E[2*omega[2]]
            ((q^2*q1^2+q^2*q1*q2)/(q^2*q1^2-q2^2))*B[(0, 0)] + ((-q^2*q1^2-q^2*q1*q2)/(-q^2*q1^2+q2^2))*B[(1, 0)] + B[(1, 1)] + ((-q^2*q1^2-q^2*q1*q2)/(-q^2*q1^2+q2^2))*B[(0, 1)]
        """
        muaff = self._L.embed_at_level(mu, 0)
        if not all(muaff.scalar(coroot) in ZZ for coroot in self._L.simple_coroots()):
            raise ValueError("%s does not lift to a level 0 element of the affine weight lattice"%muaff)
        return super(NonSymmetricMacdonaldPolynomials, self).__getitem__(mu)


    @cached_method
    def rho_prime(self): # Should be rho_prime_check
        r"""
        Return the level 0 sum of the classical fundamental weights in `L'`.

        .. SEEALSO:: :meth:`L_prime`

        EXAMPLES:

        Untwisted case::

            sage: NonSymmetricMacdonaldPolynomials("B2~").rho_prime() # CHECKME
            3/2*e[0] + 1/2*e[1]
            sage: NonSymmetricMacdonaldPolynomials("B2~").rho_prime().parent()
            Coambient space of the Root system of type ['C', 2, 1]

        Twisted case::

            sage: NonSymmetricMacdonaldPolynomials("B2~*").rho_prime() # CHECKME
            2*e[0] + e[1]
            sage: NonSymmetricMacdonaldPolynomials("B2~*").rho_prime().parent()
            Ambient space of the Root system of type ['B', 2, 1]^*
        """
        return self.L_prime().rho_classical()

    def eigenvalue_experimental(self, mu, l):
        r"""
        Return the eigenvalue of `Y^{\lambda^\vee}` acting on the macdonald polynomial `E_\mu`.

        INPUT:

        - ``mu`` -- the index `\mu` of an eigenvector
        - `l` -- an index `\lambda^\vee` of some `Y`

        .. NOTE::

            - This method is currently not used; most tests below even
              test the naive method. They are left here as a basis for
              a future implementation.

            - This is actually equivariant, as long as `s_i` does not
              fix `\lambda`.

            - This method is only really needed for
              `\lambda^\vee=\alpha^\vee_i` with `i=0,...,n`.

        See Corollary 6.11 of [Haiman06]_.

        EXAMPLES::

            sage: K = QQ['q,t'].fraction_field()
            sage: q,t = K.gens()
            sage: q1 = t
            sage: q2 = -1
            sage: KL = RootSystem(["A",1,1]).ambient_space().algebra(K)
            sage: E = NonSymmetricMacdonaldPolynomials(KL,q, q1, q2)
            sage: L0 = E.keys()
            sage: E.eigenvalues(L0([0,0])) # Checked by hand by Mark and Arun
            [1/(q*t), t]
            sage: alpha = E.Y().keys().simple_roots()
            sage: E.eigenvalue_experimental(L0([0,0]), alpha[0]) # todo: not implemented
            1/(q*t)
            sage: E.eigenvalue_experimental(L0([0,0]), alpha[1])
            t

        Some examples of eigenvalues (not mathematically checked!!!)::

            sage: E.eigenvalues(L0([1,0]))
            [t, 1/(q*t)]
            sage: E.eigenvalues(L0([0,1]))
            [1/(q^2*t), q*t]
            sage: E.eigenvalues(L0([1,1]))
            [1/(q*t), t]
            sage: E.eigenvalues(L0([2,1]))
            [t, 1/(q*t)]
            sage: E.eigenvalues(L0([-1,1]))
            [(-1)/(-q^3*t), q^2*t]
            sage: E.eigenvalues(L0([-2,1]))
            [(-1)/(-q^4*t), q^3*t]
            sage: E.eigenvalues(L0([-2,0]))
            [(-1)/(-q^3*t), q^2*t]

        Some type `B` examples::

            sage: K = QQ['q,t'].fraction_field()
            sage: q,t = K.gens()
            sage: q1 = t
            sage: q2 = -1
            sage: L = RootSystem(["B",2,1]).ambient_space()
            sage: KL = L.algebra(K)
            sage: E = NonSymmetricMacdonaldPolynomials(KL,q, q1, q2)
            sage: L0 = E.keys()
            sage: alpha = L.simple_coroots()
            sage: E.eigenvalue(L0((0,0)), alpha[0]) # not checked # not tested
            q/t
            sage: E.eigenvalue(L0((1,0)), alpha[1]) # What Mark got by hand # not tested
            q
            sage: E.eigenvalue(L0((1,0)), alpha[2]) # not checked # not tested
            t
            sage: E.eigenvalue(L0((1,0)), alpha[0]) # not checked # not tested
            1

            sage: L = RootSystem("B2~*").ambient_space()
            sage: KL = L.algebra(K)
            sage: E = NonSymmetricMacdonaldPolynomials(KL,q, q1, q2)
            sage: L0 = E.keys()
            sage: alpha = L.simple_coroots()
            sage: E.eigenvalue(L0((0,0)), alpha[0]) # assuming Mark's calculation is correct, one should get # not tested
            1/(q*t^2)

        The expected value can more or less be read off from equation
        (37), Corollary 6.15 of [Haiman06]_

        .. TODO::

            - Use proposition 6.9 of [Haiman06]_ to check the action
              of the `Y` s on monomials.

            - Generalize to any `q_1`, `q_2`.

            - Check claim by Mark: all scalar products should occur in
              the finite weight lattice, with alpha 0 being the
              appropriate projection of the affine alpha 0. Question:
              can this be emulated by being at level 0?
        """
        assert self.Y().keys().is_parent_of(l)
        L_prime = self.L_prime()
        L0 = L_prime.classical()
        I0 = L0.index_set()
        assert L0.is_parent_of(mu)
        # Should we view mu as a translation, and ask for its alcove walk?
        muaff = self.affine_lift(mu) # embeds mu at level 1 in L_prime
        w = reversed(mu.reduced_word(I0, positive=False)) # the reduced word for w_\mu, Prop. 6.9 of [Haiman06]_
        # mu should be scaled to make sure it implements a translation
        #w = reversed(L.reduced_word_of_translation(L(mu)))
        #x = L.embed_at_level(L0.rho(),1)
        #x = L.rho() / L.rho().level()
        x = self.rho_prime()
        l = self.L_prime().coroot_lattice()(l) # there might need to be a `nu` here
        for i in w:
            x = x.simple_reflection(i)
        q1,q2 = self.hecke_parameters(1) # TODO: clean up
        t = -q2/q1  # TODO: generalize for any eigenvalue
        # In type BC, maybe this should be q^...*a[0]
        return self._q**(-muaff.scalar(l)) * t**(-x.scalar(l))

    def seed(self, mu):
        r"""
        Return `E_\mu` for `\mu` minuscule, i.e. in the fundamental alcove.

        INPUT:

        - ``mu`` -- the index `\mu` of an eigenvector

        EXAMPLES::

            sage: E = NonSymmetricMacdonaldPolynomials(["A",2,1])
            sage: omega = E.keys().fundamental_weights()
            sage: E.seed(omega[1])
            B[(1, 0, 0)]
        """
        return self.KL0().monomial(mu)

    def symmetric_macdonald_polynomial(self, mu):
        r"""
        Return the symmetric Macdonald polynomial indexed by `\mu`.

        INPUT:

        - ``mu`` -- a dominant weight `\mu`

        .. WARNING::

            The result is Weyl-symmetric only for Hecke parameters of
            the form `q_1=v` and `q_2=-1/v`.  In general the value of
            `v` below, should be the square root of `-q_1/q_2`, but the
            use of `q_1=t` and `q_2=-1` results in nonintegral powers of `t`.

        EXAMPLES::

            sage: K = QQ['q,v,t'].fraction_field()
            sage: q,v,t = K.gens()
            sage: E = NonSymmetricMacdonaldPolynomials(['A',2,1], q, v, -1/v)
            sage: om = E.L0().fundamental_weights()
            sage: E.symmetric_macdonald_polynomial(om[2])
            B[(1, 1, 0)] + B[(1, 0, 1)] + B[(0, 1, 1)]
            sage: E.symmetric_macdonald_polynomial(2*om[1])
            ((q*v^2+v^2-q-1)/(q*v^2-1))*B[(1, 1, 0)] + ((q*v^2+v^2-q-1)/(q*v^2-1))*B[(1, 0, 1)] + B[(2, 0, 0)] + ((q*v^2+v^2-q-1)/(q*v^2-1))*B[(0, 1, 1)] + B[(0, 2, 0)] + B[(0, 0, 2)]
            sage: f = E.symmetric_macdonald_polynomial(E.L0()((2,1,0))); f
            ((2*q*v^4+v^4-q*v^2+v^2-q-2)/(q*v^4-1))*B[(1, 1, 1)] + B[(1, 2, 0)] + B[(1, 0, 2)] + B[(2, 1, 0)] + B[(2, 0, 1)] + B[(0, 1, 2)] + B[(0, 2, 1)]

        We compare with the type `A` Macdonald polynomials
        coming from symmetric functions::

            sage: P = SymmetricFunctions(K).macdonald().P()
            sage: g = P[2,1].expand(3); g
            x0^2*x1 + x0*x1^2 + x0^2*x2 + (2*q*t^2 - q*t - q  + t^2 + t - 2)/(q*t^2 - 1)*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2
            sage: fe = f.expand(g.parent().gens()); fe
            x0^2*x1 + x0*x1^2 + x0^2*x2 + (2*q*v^4 - q*v^2 - q + v^4 + v^2 - 2)/(q*v^4 - 1)*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2
            sage: g.map_coefficients(lambda x: x.subs(t=v*v)) == fe
            True

            sage: E = NonSymmetricMacdonaldPolynomials(['C',3,1], q, v, -1/v)
            sage: om = E.L0().fundamental_weights()
            sage: E.symmetric_macdonald_polynomial(om[1]+om[2])
            B[(-2, -1, 0)] + B[(-2, 1, 0)] + B[(-2, 0, -1)] + B[(-2, 0, 1)] + ((4*q^3*v^14+2*q^2*v^14-2*q^3*v^12+2*q^2*v^12-2*q^3*v^10+q*v^12-5*q^2*v^10-5*q*v^4+q^2*v^2-2*v^4+2*q*v^2-2*v^2+2*q+4)/(q^3*v^14-q^2*v^10-q*v^4+1))*B[(-1, 0, 0)] + B[(-1, -2, 0)] + ((2*q*v^4+v^4-q*v^2+v^2-q-2)/(q*v^4-1))*B[(-1, -1, -1)] + ((2*q*v^4+v^4-q*v^2+v^2-q-2)/(q*v^4-1))*B[(-1, -1, 1)] + ((2*q*v^4+v^4-q*v^2+v^2-q-2)/(q*v^4-1))*B[(-1, 1, -1)] + ((2*q*v^4+v^4-q*v^2+v^2-q-2)/(q*v^4-1))*B[(-1, 1, 1)] + B[(-1, 2, 0)] + B[(-1, 0, -2)] + B[(-1, 0, 2)] + ((4*q^3*v^14+2*q^2*v^14-2*q^3*v^12+2*q^2*v^12-2*q^3*v^10+q*v^12-5*q^2*v^10-5*q*v^4+q^2*v^2-2*v^4+2*q*v^2-2*v^2+2*q+4)/(q^3*v^14-q^2*v^10-q*v^4+1))*B[(1, 0, 0)] + B[(1, -2, 0)] + ((2*q*v^4+v^4-q*v^2+v^2-q-2)/(q*v^4-1))*B[(1, -1, -1)] + ((2*q*v^4+v^4-q*v^2+v^2-q-2)/(q*v^4-1))*B[(1, -1, 1)] + ((2*q*v^4+v^4-q*v^2+v^2-q-2)/(q*v^4-1))*B[(1, 1, -1)] + ((2*q*v^4+v^4-q*v^2+v^2-q-2)/(q*v^4-1))*B[(1, 1, 1)] + B[(1, 2, 0)] + B[(1, 0, -2)] + B[(1, 0, 2)] + B[(2, -1, 0)] + B[(2, 1, 0)] + B[(2, 0, -1)] + B[(2, 0, 1)] + B[(0, -2, -1)] + B[(0, -2, 1)] + ((-4*q^3*v^14-2*q^2*v^14+2*q^3*v^12-2*q^2*v^12+2*q^3*v^10-q*v^12+5*q^2*v^10+5*q*v^4-q^2*v^2+2*v^4-2*q*v^2+2*v^2-2*q-4)/(-q^3*v^14+q^2*v^10+q*v^4-1))*B[(0, -1, 0)] + B[(0, -1, -2)] + B[(0, -1, 2)] + ((-4*q^3*v^14-2*q^2*v^14+2*q^3*v^12-2*q^2*v^12+2*q^3*v^10-q*v^12+5*q^2*v^10+5*q*v^4-q^2*v^2+2*v^4-2*q*v^2+2*v^2-2*q-4)/(-q^3*v^14+q^2*v^10+q*v^4-1))*B[(0, 1, 0)] + B[(0, 1, -2)] + B[(0, 1, 2)] + B[(0, 2, -1)] + B[(0, 2, 1)] + ((4*q^3*v^14+2*q^2*v^14-2*q^3*v^12+2*q^2*v^12-2*q^3*v^10+q*v^12-5*q^2*v^10-5*q*v^4+q^2*v^2-2*v^4+2*q*v^2-2*v^2+2*q+4)/(q^3*v^14-q^2*v^10-q*v^4+1))*B[(0, 0, -1)] + ((4*q^3*v^14+2*q^2*v^14-2*q^3*v^12+2*q^2*v^12-2*q^3*v^10+q*v^12-5*q^2*v^10-5*q*v^4+q^2*v^2-2*v^4+2*q*v^2-2*v^2+2*q+4)/(q^3*v^14-q^2*v^10-q*v^4+1))*B[(0, 0, 1)]

        An example for type `G`::

            sage: E = NonSymmetricMacdonaldPolynomials(['G',2,1], q, v, -1/v)
            sage: om = E.L0().fundamental_weights()
            sage: E.symmetric_macdonald_polynomial(2*om[1])
            ((3*q^6*v^22+3*q^5*v^22-3*q^6*v^20+q^4*v^22-4*q^5*v^20+q^4*v^18-q^5*v^16+q^3*v^18-2*q^4*v^16+q^5*v^14-q^3*v^16+q^4*v^14-4*q^4*v^12+q^2*v^14+q^5*v^10-8*q^3*v^12+4*q^4*v^10-4*q^2*v^12+8*q^3*v^10-q*v^12-q^4*v^8+4*q^2*v^10-q^2*v^8+q^3*v^6-q*v^8+2*q^2*v^6-q^3*v^4+q*v^6-q^2*v^4+4*q*v^2-q^2+3*v^2-3*q-3)/(q^6*v^22-q^5*v^20-q^4*v^12-q^3*v^12+q^3*v^10+q^2*v^10+q*v^2-1))*B[(0, 0, 0)] + ((q*v^2+v^2-q-1)/(q*v^2-1))*B[(-2, 1, 1)] + B[(-2, 2, 0)] + B[(-2, 0, 2)] + ((-q*v^2-v^2+q+1)/(-q*v^2+1))*B[(-1, -1, 2)] + ((2*q^4*v^12+2*q^3*v^12-2*q^4*v^10-2*q^3*v^10+q^2*v^8-q^3*v^6+q*v^8-2*q^2*v^6+q^3*v^4-q*v^6+q^2*v^4-2*q*v^2-2*v^2+2*q+2)/(q^4*v^12-q^3*v^10-q*v^2+1))*B[(-1, 1, 0)] + ((-q*v^2-v^2+q+1)/(-q*v^2+1))*B[(-1, 2, -1)] + ((2*q^4*v^12+2*q^3*v^12-2*q^4*v^10-2*q^3*v^10+q^2*v^8-q^3*v^6+q*v^8-2*q^2*v^6+q^3*v^4-q*v^6+q^2*v^4-2*q*v^2-2*v^2+2*q+2)/(q^4*v^12-q^3*v^10-q*v^2+1))*B[(-1, 0, 1)] + ((-q*v^2-v^2+q+1)/(-q*v^2+1))*B[(1, -2, 1)] + ((-2*q^4*v^12-2*q^3*v^12+2*q^4*v^10+2*q^3*v^10-q^2*v^8+q^3*v^6-q*v^8+2*q^2*v^6-q^3*v^4+q*v^6-q^2*v^4+2*q*v^2+2*v^2-2*q-2)/(-q^4*v^12+q^3*v^10+q*v^2-1))*B[(1, -1, 0)] + ((-q*v^2-v^2+q+1)/(-q*v^2+1))*B[(1, 1, -2)] + ((-2*q^4*v^12-2*q^3*v^12+2*q^4*v^10+2*q^3*v^10-q^2*v^8+q^3*v^6-q*v^8+2*q^2*v^6-q^3*v^4+q*v^6-q^2*v^4+2*q*v^2+2*v^2-2*q-2)/(-q^4*v^12+q^3*v^10+q*v^2-1))*B[(1, 0, -1)] + B[(2, -2, 0)] + ((q*v^2+v^2-q-1)/(q*v^2-1))*B[(2, -1, -1)] + B[(2, 0, -2)] + B[(0, -2, 2)] + ((-2*q^4*v^12-2*q^3*v^12+2*q^4*v^10+2*q^3*v^10-q^2*v^8+q^3*v^6-q*v^8+2*q^2*v^6-q^3*v^4+q*v^6-q^2*v^4+2*q*v^2+2*v^2-2*q-2)/(-q^4*v^12+q^3*v^10+q*v^2-1))*B[(0, -1, 1)] + ((2*q^4*v^12+2*q^3*v^12-2*q^4*v^10-2*q^3*v^10+q^2*v^8-q^3*v^6+q*v^8-2*q^2*v^6+q^3*v^4-q*v^6+q^2*v^4-2*q*v^2-2*v^2+2*q+2)/(q^4*v^12-q^3*v^10-q*v^2+1))*B[(0, 1, -1)] + B[(0, 2, -2)]

        """
        if self.cartan_type().classical() != mu.parent().cartan_type() or not mu.is_dominant():
            raise ValueError("%s must be a dominant weight for the classical subrootsystem of %s" % (mu, self.cartan_type()))
        v = self._q1
        KL0 = self.KL0()
        s = KL0.zero()
        # efficiently compute the finite Hecke symmetrization of the
        # nonsymmetric Macdonald polynomial of the dominant weight mu
        # by searching the Weyl orbit of mu and remembering
        Torbit = {}
        for c in mu._orbit_iter():
            i = c.first_descent()
            if i is None:
                Torbit[c] = self[mu] # the nonsymmetric Macdonald polynomial of mu
            else:
                Torbit[c] = v * self._T.Tw([i])(Torbit[c.simple_reflection(i)])
            s = s + Torbit[c]
        return s

