# sage.doctest: needs sage.combinat sage.graphs sage.modules
r"""
Iwahori-Hecke Algebras

AUTHORS:

- Daniel Bump, Nicolas Thiery (2010): Initial version

- Brant Jones, Travis Scrimshaw, Andrew Mathas (2013):
  Moved into the category framework and implemented the
  Kazhdan-Lusztig `C` and `C^{\prime}` bases

- Chase Meadors, Tianyuan Xu (2021):
  Implemented direct computation of products in the
  `C^{\prime}` basis using du Cloux's Coxeter3 package
"""
# ****************************************************************************
#  Copyright (C) 2013 Brant Jones <brant at math.jmu.edu>
#                     Daniel Bump <bump at match.stanford.edu>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from functools import cmp_to_key
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.realizations import Realizations, Category_realization_of_parent
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.finite_dimensional_algebras_with_basis import FiniteDimensionalAlgebrasWithBasis
from sage.categories.coxeter_groups import CoxeterGroups
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.arith.misc import is_square
from sage.combinat.root_system.coxeter_group import CoxeterGroup
from sage.sets.family import Family
from sage.combinat.free_module import CombinatorialFreeModule


def normalized_laurent_polynomial(R, p):
    r"""
    Return a normalized version of the (Laurent polynomial) ``p`` in the
    ring ``R``.

    Various ring operations in ``sage`` return an element of the field of
    fractions of the parent ring even though the element is "known" to belong to
    the base ring. This function is a hack to recover from this. This occurs
    somewhat haphazardly with Laurent polynomial rings::

        sage: R.<q> = LaurentPolynomialRing(ZZ)
        sage: [type(c) for c in (q**-1).coefficients()]
        [<class 'sage.rings.integer.Integer'>]

    It also happens in any ring when dividing by units::

        sage: type(3/1)
        <class 'sage.rings.rational.Rational'>
        sage: type(-1/-1)
        <class 'sage.rings.rational.Rational'>

    This function is a variation on a suggested workaround of Nils Bruin.

    EXAMPLES::

        sage: from sage.algebras.iwahori_hecke_algebra import normalized_laurent_polynomial
        sage: type(normalized_laurent_polynomial(ZZ, 3/1))
        <class 'sage.rings.integer.Integer'>
        sage: R.<q> = LaurentPolynomialRing(ZZ)
        sage: [type(c) for c in normalized_laurent_polynomial(R, q**-1).coefficients()]
        [<class 'sage.rings.integer.Integer'>]
        sage: R.<u,v> = LaurentPolynomialRing(ZZ,2)
        sage: p = normalized_laurent_polynomial(R, 2*u**-1*v**-1 + u*v)
        sage: ui = normalized_laurent_polynomial(R, u^-1)
        sage: vi = normalized_laurent_polynomial(R, v^-1)
        sage: p(ui, vi)
        2*u*v + u^-1*v^-1
        sage: q = u+v+ui
        sage: q(ui, vi)
        u + v^-1 + u^-1
    """
    try:
        return R({k: R._base(c) for k, c in p.monomial_coefficients().items()})
    except (AttributeError, TypeError):
        return R(p)


def index_cmp(x, y):
    """
    Compare two term indices ``x`` and ``y`` by Bruhat order, then by word
    length, and then by the generic comparison.

    EXAMPLES::

        sage: from sage.algebras.iwahori_hecke_algebra import index_cmp
        sage: W = WeylGroup(['A',2,1])
        sage: x = W.from_reduced_word([0,1])
        sage: y = W.from_reduced_word([0,2,1])
        sage: x.bruhat_le(y)
        True
        sage: index_cmp(x, y)
        1
    """
    if x.bruhat_le(y) or x.length() < y.length():
        return 1
    if y.bruhat_le(x) or x.length() > y.length():
        return -1
    # fallback case, in order to define a total order
    if x < y:
        return -1
    if x > y:
        return 1
    return 0


sorting_key = cmp_to_key(index_cmp)


class IwahoriHeckeAlgebra(Parent, UniqueRepresentation):
    r"""
    The Iwahori-Hecke algebra of the Coxeter group ``W``
    with the specified parameters.

    INPUT:

    - ``W`` -- a Coxeter group or Cartan type
    - ``q1`` -- a parameter
    - ``q2`` -- (default: ``-1``) another parameter
    - ``base_ring`` -- (default: ``q1.parent()``) a ring containing ``q1``
      and ``q2``

    The Iwahori-Hecke algebra [Iwa1964]_ is a deformation of the group algebra of
    a Weyl group or, more generally, a Coxeter group. These algebras are
    defined by generators and relations and they depend on a deformation
    parameter `q`. Taking `q = 1`, as in the following example, gives a ring
    isomorphic to the group algebra of the corresponding Coxeter group.

    Let `(W, S)` be a Coxeter system and let `R` be a commutative ring
    containing elements `q_1` and `q_2`. Then the *Iwahori-Hecke algebra*
    `H = H_{q_1,q_2}(W,S)` of `(W,S)` with parameters `q_1` and `q_2` is the
    unital associative algebra with generators `\{T_s \mid s\in S\}` and
    relations:

    .. MATH::

        \begin{aligned}
            (T_s - q_1)(T_s - q_2) &= 0\\
            T_r T_s T_r \cdots &= T_s T_r T_s \cdots,
        \end{aligned}

    where the number of terms on either side of the second relations (the braid
    relations) is the order of `rs` in the Coxeter group `W`, for `r,s \in S`.

    Iwahori-Hecke algebras are fundamental in many areas of mathematics,
    ranging from the representation theory of Lie groups and quantum groups,
    to knot theory and statistical mechanics. For more information see,
    for example, [KL1979]_, [HKP2010]_, [Jon1987]_ and
    :wikipedia:`Iwahori-Hecke_algebra`.

    .. RUBRIC:: Bases

    A reduced expression for an element `w \in W` is any minimal length
    word `w = s_1 \cdots s_k`, with `s_i \in S`. If `w = s_1 \cdots s_k` is a
    reduced expression for `w` then Matsumoto's Monoid Lemma implies that
    `T_w = T_{s_1} \cdots T_{s_k}` depends on `w` and not on the choice of
    reduced expressions. Moreover, `\{ T_w \mid w\in W \}` is a basis for the
    Iwahori-Hecke algebra `H` and

    .. MATH::

        T_s T_w = \begin{cases}
           T_{sw}, &                     \text{if } \ell(sw) = \ell(w)+1,\\
           (q_1+q_2)T_w -q_1q_2 T_{sw}, & \text{if } \ell(sw) = \ell(w)-1.
        \end{cases}

    The `T`-basis of `H` is implemented for any choice of parameters
    ``q_1`` and ``q_2``::

        sage: R.<u,v> = LaurentPolynomialRing(ZZ,2)
        sage: H = IwahoriHeckeAlgebra('A3', u,v)
        sage: T = H.T()
        sage: T[1]
        T[1]
        sage: T[1,2,1] + T[2]
        T[1,2,1] + T[2]
        sage: T[1] * T[1,2,1]
        (u+v)*T[1,2,1] + (-u*v)*T[2,1]
        sage: T[1]^-1
        (-u^-1*v^-1)*T[1] + (v^-1+u^-1)

    Working over the Laurent polynomial ring `Z[q^{\pm 1/2}]` Kazhdan and
    Lusztig proved that there exist two distinguished bases
    `\{ C^{\prime}_w \mid w \in W \}` and `\{ C_w \mid w \in W \}` of `H`
    which are uniquely determined by the properties that they are invariant
    under the bar involution on `H` and have triangular transitions matrices
    with polynomial entries of a certain form with the `T`-basis;
    see [KL1979]_ for a precise statement.

    It turns out that the Kazhdan-Lusztig bases can be defined (by
    specialization) in `H` whenever `-q_1 q_2` is a square in the base ring.
    The Kazhdan-Lusztig bases are implemented inside `H` whenever `-q_1 q_2`
    has a square root::

        sage: H = IwahoriHeckeAlgebra('A3', u^2, -v^2)
        sage: T = H.T(); Cp = H.Cp(); C = H.C()
        sage: T(Cp[1])
        (u^-1*v^-1)*T[1] + (u^-1*v)
        sage: T(C[1])
        (u^-1*v^-1)*T[1] + (-u*v^-1)
        sage: Cp(C[1])
        Cp[1] + (-u*v^-1-u^-1*v)
        sage: elt = Cp[2]*Cp[3]+C[1]; elt
        Cp[2,3] + Cp[1] + (-u*v^-1-u^-1*v)
        sage: c = C(elt); c
        C[2,3] + C[1] + (u*v^-1+u^-1*v)*C[3] + (u*v^-1+u^-1*v)*C[2] + (u^2*v^-2+2+u^-2*v^2)
        sage: t = T(c); t
        (u^-2*v^-2)*T[2,3] + (u^-1*v^-1)*T[1] + (u^-2)*T[3] + (u^-2)*T[2] + (-u*v^-1+u^-2*v^2)
        sage: Cp(t)
        Cp[2,3] + Cp[1] + (-u*v^-1-u^-1*v)
        sage: Cp(c)
        Cp[2,3] + Cp[1] + (-u*v^-1-u^-1*v)

    The conversions to and from the Kazhdan-Lusztig bases are done behind the
    scenes whenever the Kazhdan-Lusztig bases are well-defined. Once a suitable
    Iwahori-Hecke algebra is defined they will work without further
    intervention.

    For example, with the "standard parameters", so that
    `(T_r-q^2)(T_r+1) = 0`::

        sage: R.<q> = LaurentPolynomialRing(ZZ)
        sage: H = IwahoriHeckeAlgebra('A3', q^2)
        sage: T=H.T(); Cp=H.Cp(); C=H.C()
        sage: C(T[1])
        q*C[1] + q^2
        sage: elt = Cp(T[1,2,1]); elt
        q^3*Cp[1,2,1] - q^2*Cp[2,1] - q^2*Cp[1,2] + q*Cp[1] + q*Cp[2] - 1
        sage: C(elt)
        q^3*C[1,2,1] + q^4*C[2,1] + q^4*C[1,2] + q^5*C[1] + q^5*C[2] + q^6

    With the "normalized presentation", so that `(T_r-q)(T_r+q^{-1}) = 0`::

        sage: R.<q> = LaurentPolynomialRing(ZZ)
        sage: H = IwahoriHeckeAlgebra('A3', q, -q^-1)
        sage: T=H.T(); Cp=H.Cp(); C=H.C()
        sage: C(T[1])
        C[1] + q
        sage: elt = Cp(T[1,2,1]); elt
        Cp[1,2,1] - (q^-1)*Cp[2,1] - (q^-1)*Cp[1,2] + (q^-2)*Cp[1] + (q^-2)*Cp[2] - (q^-3)
        sage: C(elt)
        C[1,2,1] + q*C[2,1] + q*C[1,2] + q^2*C[1] + q^2*C[2] + q^3

    In the group algebra, so that `(T_r-1)(T_r+1) = 0`::

        sage: H = IwahoriHeckeAlgebra('A3', 1)
        sage: T=H.T(); Cp=H.Cp(); C=H.C()
        sage: C(T[1])
        C[1] + 1
        sage: Cp(T[1,2,1])
        Cp[1,2,1] - Cp[2,1] - Cp[1,2] + Cp[1] + Cp[2] - 1
        sage: C(_)
        C[1,2,1] + C[2,1] + C[1,2] + C[1] + C[2] + 1

    On the other hand, if the Kazhdan-Lusztig bases are not well-defined (when
    `-q_1 q_2` is not a square), attempting to use the Kazhdan-Lusztig bases
    triggers an error::

        sage: R.<q>=LaurentPolynomialRing(ZZ)
        sage: H = IwahoriHeckeAlgebra('A3', q)
        sage: C=H.C()
        Traceback (most recent call last):
        ...
        ValueError: The Kazhdan_Lusztig bases are defined only when -q_1*q_2 is a square

    We give an example in affine type::

        sage: R.<v> = LaurentPolynomialRing(ZZ)
        sage: H = IwahoriHeckeAlgebra(['A',2,1], v^2)
        sage: T=H.T(); Cp=H.Cp(); C=H.C()
        sage: C(T[1,0,2])
        v^3*C[1,0,2] + v^4*C[1,0] + v^4*C[0,2] + v^4*C[1,2]
         + v^5*C[0] + v^5*C[2] + v^5*C[1] + v^6
        sage: Cp(T[1,0,2])
        v^3*Cp[1,0,2] - v^2*Cp[1,0] - v^2*Cp[0,2] - v^2*Cp[1,2]
         + v*Cp[0] + v*Cp[2] + v*Cp[1] - 1
        sage: T(C[1,0,2])
        (v^-3)*T[1,0,2] - (v^-1)*T[1,0] - (v^-1)*T[0,2] - (v^-1)*T[1,2]
         + v*T[0] + v*T[2] + v*T[1] - v^3
        sage: T(Cp[1,0,2])
        (v^-3)*T[1,0,2] + (v^-3)*T[1,0] + (v^-3)*T[0,2] + (v^-3)*T[1,2]
         + (v^-3)*T[0] + (v^-3)*T[2] + (v^-3)*T[1] + (v^-3)

    EXAMPLES:

    We start by creating a Iwahori-Hecke algebra together with the three bases
    for these algebras that are currently supported::

        sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
        sage: H = IwahoriHeckeAlgebra('A3', v**2)
        sage: T = H.T()
        sage: C = H.C()
        sage: Cp = H.Cp()

    It is also possible to define these three bases quickly using
    the :meth:`inject_shorthands` method.

    Next we create our generators for the `T`-basis and do some basic
    computations and conversions between the bases::

        sage: T1,T2,T3 = T.algebra_generators()
        sage: T1 == T[1]
        True
        sage: T1*T2 == T[1,2]
        True
        sage: T1 + T2
        T[1] + T[2]
        sage: T1*T1
        -(1-v^2)*T[1] + v^2
        sage: (T1 + T2)*T3 + T1*T1 - (v + v^-1)*T2
        T[3,1] + T[2,3] - (1-v^2)*T[1] - (v^-1+v)*T[2] + v^2
        sage: Cp(T1)
        v*Cp[1] - 1
        sage: Cp((v^1 - 1)*T1*T2 - T3)
        -(v^2-v^3)*Cp[1,2] + (v-v^2)*Cp[1] - v*Cp[3] + (v-v^2)*Cp[2] + v
        sage: C(T1)
        v*C[1] + v^2
        sage: p = C(T2*T3 - v*T1); p
        v^2*C[2,3] - v^2*C[1] + v^3*C[3] + v^3*C[2] - (v^3-v^4)
        sage: Cp(p)
        v^2*Cp[2,3] - v^2*Cp[1] - v*Cp[3] - v*Cp[2] + (1+v)
        sage: Cp(T2*T3 - v*T1)
        v^2*Cp[2,3] - v^2*Cp[1] - v*Cp[3] - v*Cp[2] + (1+v)

    In addition to explicitly creating generators, we have two shortcuts to
    basis elements. The first is by using elements of the underlying Coxeter
    group, the other is by using reduced words::

        sage: s1,s2,s3 = H.coxeter_group().gens()
        sage: T[s1*s2*s1*s3] == T[1,2,1,3]
        True
        sage: T[1,2,1,3] == T1*T2*T1*T3
        True

    TESTS:

    We check the defining properties of the bases::

        sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
        sage: H = IwahoriHeckeAlgebra('A3', v**2)
        sage: W = H.coxeter_group()
        sage: T = H.T()
        sage: C = H.C()
        sage: Cp = H.Cp()
        sage: T(Cp[1])
        (v^-1)*T[1] + (v^-1)
        sage: T(C[1])
        (v^-1)*T[1] - v
        sage: C(Cp[1])
        C[1] + (v^-1+v)
        sage: Cp(C[1])
        Cp[1] - (v^-1+v)
        sage: all(C[x] == C[x].bar() for x in W) # long time
        True
        sage: all(Cp[x] == Cp[x].bar() for x in W) # long time
        True
        sage: all(T(C[x]).bar() == T(C[x]) for x in W) # long time
        True
        sage: all(T(Cp[x]).bar() == T(Cp[x]) for x in W) # long time
        True
        sage: KL = KazhdanLusztigPolynomial(W, v)
        sage: term = lambda x,y: (-1)^y.length() * v^(-2*y.length()) * KL.P(y, x).substitute(v=v^-2)*T[y]
        sage: all(T(C[x]) == (-v)^x.length()*sum(term(x,y) for y in W) for x in W) # long time
        True
        sage: all(T(Cp[x]) == v^-x.length()*sum(KL.P(y,x).substitute(v=v^2)*T[y] for y in W) for x in W) # long time
        True

    We check conversion between the bases for type `B_2` as well as some of
    the defining properties::

        sage: H = IwahoriHeckeAlgebra(['B',2], v**2)
        sage: W = H.coxeter_group()
        sage: T = H.T()
        sage: C = H.C()
        sage: Cp = H.Cp()
        sage: all(T[x] == T(C(T[x])) for x in W) # long time
        True
        sage: all(T[x] == T(Cp(T[x])) for x in W) # long time
        True
        sage: all(C[x] == C(T(C[x])) for x in W) # long time
        True
        sage: all(C[x] == C(Cp(C[x])) for x in W) # long time
        True
        sage: all(Cp[x] == Cp(T(Cp[x])) for x in W) # long time
        True
        sage: all(Cp[x] == Cp(C(Cp[x])) for x in W) # long time
        True
        sage: all(T(C[x]).bar() == T(C[x]) for x in W) # long time
        True
        sage: all(T(Cp[x]).bar() == T(Cp[x]) for x in W) # long time
        True
        sage: KL = KazhdanLusztigPolynomial(W, v)
        sage: term = lambda x,y: (-1)^y.length() * v^(-2*y.length()) * KL.P(y, x).substitute(v=v^-2)*T[y]
        sage: all(T(C[x]) == (-v)^x.length()*sum(term(x,y) for y in W) for x in W) # long time
        True
        sage: all(T(Cp[x]) == v^-x.length()*sum(KL.P(y,x).substitute(v=v^2)*T[y] for y in W) for x in W) # long time
        True

    .. TODO::

        Implement multi-parameter Iwahori-Hecke algebras together with their
        Kazhdan-Lusztig bases. That is, Iwahori-Hecke algebras with (possibly)
        different parameters for each conjugacy class of simple reflections
        in the underlying Coxeter group.

    .. TODO::

        When given "generic parameters" we should return the generic
        Iwahori-Hecke algebra with these parameters and allow the user to
        work inside this algebra rather than doing calculations behind the
        scenes in a copy of the generic Iwahori-Hecke algebra. The main
        problem is that it is not clear how to recognise when the
        parameters are "generic".
    """
    @staticmethod
    def __classcall_private__(cls, W, q1, q2=-1, base_ring=None):
        r"""
        TESTS::

            sage: H = IwahoriHeckeAlgebra("A2", 1)
            sage: W = CoxeterGroup("A2")
            sage: H.coxeter_group() == W
            True
            sage: H.cartan_type() == CartanType("A2")
            True
            sage: H._q2 == -1
            True
            sage: H2 = IwahoriHeckeAlgebra(W, QQ(1), base_ring=ZZ)
            sage: H is H2
            True
        """
        if W not in CoxeterGroups():
            W = CoxeterGroup(W)
        if base_ring is None:
            base_ring = q1.parent()
        else:
            q1 = base_ring(q1)
        q2 = base_ring(q2)
        return super().__classcall__(cls, W, q1, q2, base_ring)

    def __init__(self, W, q1, q2, base_ring):
        r"""
        Initialize and return the two parameter Iwahori-Hecke algebra ``self``.

        EXAMPLES::

            sage: R.<q1,q2> = QQ[]
            sage: H = IwahoriHeckeAlgebra("A2", q1, q2=q2, base_ring=Frac(R))
            sage: TestSuite(H).run()
        """
        self._W = W
        self._coxeter_type = W.coxeter_type()

        self._q1 = q1
        self._q2 = q2

        # Used when multiplying generators: minor speed-up as it avoids the
        # need to constantly add and multiply the parameters when applying the
        # quadratic relation: T^2 = (q1+q2)T - q1*q2
        self._q_sum = q1 + q2
        self._q_prod = -q1 * q2

        # If -q1*q2 is a square then it makes sense to talk of he Kazhdan-Lusztig
        # basis of the Iwhaori-Hecke algebra. In this case we set
        # self._root=\sqrt{q1*q2}. The Kazhdan-Lusztig bases will be computed in
        # the generic case behind the scenes and then specialized to this # algebra.
        is_Square, root = is_square(self._q_prod, root=True)
        if is_Square:
            # Attach the generic Hecke algebra and the basis change maps
            self._root = root
            self._generic_iwahori_hecke_algebra = IwahoriHeckeAlgebra_nonstandard(W)
            self._shorthands = ['C', 'Cp', 'T']
        else:
            # Can we actually remove the bases C and Cp in this case?
            self._root = None
            self._shorthands = ['T']

        # if 2 is a unit in the base ring then add th A and B bases
        try:
            base_ring(base_ring.one() / 2)
            self._shorthands.extend(['A', 'B'])
        except (TypeError, ZeroDivisionError):
            pass

        if W.is_finite():
            self._category = FiniteDimensionalAlgebrasWithBasis(base_ring)
        else:
            self._category = AlgebrasWithBasis(base_ring)
        Parent.__init__(self, base=base_ring, category=self._category.WithRealizations())

        self._is_generic = False  # needed for initialisation of _KLHeckeBasis

        # The following is used by the bar involution = self._bar_on_coefficients
        try:
            self._inverse_base_ring_generators = {g: self.base_ring()(g) ** -1
                    for g in self.base_ring().variable_names()}
        except TypeError:
            self._inverse_base_ring_generators = {}

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: R.<q1,q2> = QQ[]
            sage: IwahoriHeckeAlgebra("A2", q1**2, q2**2, base_ring=Frac(R))
            Iwahori-Hecke algebra of type A2 in q1^2,q2^2 over Fraction Field of Multivariate Polynomial Ring in q1, q2 over Rational Field
        """
        try:
            ct = self._coxeter_type._repr_(compact=True)
        except TypeError:
            ct = repr(self._coxeter_type)
        return "Iwahori-Hecke algebra of type {} in {},{} over {}".format(
            ct, self._q1, self._q2, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: R.<q1,q2> = QQ[]
            sage: H = IwahoriHeckeAlgebra("A2", q1**2, q2**2, base_ring=Frac(R))
            sage: latex(H)
            \mathcal{H}_{q_{1}^{2},q_{2}^{2}}\left(A_{2},
             \mathrm{Frac}(\Bold{Q}[q_{1}, q_{2}])\right)
            sage: R.<q> = LaurentPolynomialRing(ZZ)
            sage: H = IwahoriHeckeAlgebra("A2", q)
            sage: latex(H)
            \mathcal{H}_{q,-1}\left(A_{2}, \Bold{Z}[q^{\pm 1}]\right)
        """
        from sage.misc.latex import latex
        return "\\mathcal{{H}}_{{{},{}}}\\left({}, {}\\right)".format(latex(self._q1),
                latex(self._q2), latex(self._coxeter_type), latex(self.base_ring()))

    def _bar_on_coefficients(self, c):
        r"""
        Given a Laurent polynomial ``c`` return the Laurent polynomial obtained
        by applying the (generic) bar involution to ``c`` .

        This is the ring homomorphism of Laurent polynomials in
        `\ZZ[u,u^{-1},v,v^{-1}]` which sends `u` to `u^{-1}` and `v`
        to `v^{-1}`.

        EXAMPLES::

            sage: R.<q>=LaurentPolynomialRing(ZZ)
            sage: H = IwahoriHeckeAlgebra("A3",q^2)
            sage: H._bar_on_coefficients(q)
            q^-1
        """
        return normalized_laurent_polynomial(self._base, c).substitute(**self._inverse_base_ring_generators)

    def coxeter_type(self):
        r"""
        Return the Coxeter type of ``self``.

        EXAMPLES::

            sage: IwahoriHeckeAlgebra("D4", 1).coxeter_type()
            Coxeter type of ['D', 4]
        """
        return self._coxeter_type

    def cartan_type(self):
        r"""
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: IwahoriHeckeAlgebra("D4", 1).cartan_type()
            ['D', 4]
        """
        try:
            return self._coxeter_type.cartan_type()
        except AttributeError:
            return None

    def coxeter_group(self):
        r"""
        Return the Coxeter group of ``self``.

        EXAMPLES::

            sage: IwahoriHeckeAlgebra("B2", 1).coxeter_group()
            Finite Coxeter group over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? with Coxeter matrix:
            [1 4]
            [4 1]
        """
        return self._W

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the `T`-basis).

        EXAMPLES::

            sage: H = IwahoriHeckeAlgebra("B2", 1)
            sage: H.a_realization()
            Iwahori-Hecke algebra of type B2 in 1,-1 over Integer Ring in the T-basis
        """
        return self.T()

    def q1(self):
        """
        Return the parameter `q_1` of ``self``.

        EXAMPLES::

            sage: H = IwahoriHeckeAlgebra("B2", 1)
            sage: H.q1()
            1
        """
        return self._q1

    def q2(self):
        """
        Return the parameter `q_2` of ``self``.

        EXAMPLES::

            sage: H = IwahoriHeckeAlgebra("B2", 1)
            sage: H.q2()
            -1
        """
        return self._q2

    class _BasesCategory(Category_realization_of_parent):
        r"""
        The category of bases of a Iwahori-Hecke algebra.
        """
        def __init__(self, base):
            r"""
            Initialize the bases of a Iwahori-Hecke algebra.

            INPUT:

            - ``base`` -- a Iwahori-Hecke algebra

            TESTS::

                sage: H = IwahoriHeckeAlgebra("B2", 1)
                sage: bases = H._BasesCategory()
                sage: H.T() in bases
                True
            """
            Category_realization_of_parent.__init__(self, base)

        def super_categories(self):
            r"""
            The super categories of ``self``.

            EXAMPLES::

                sage: H = IwahoriHeckeAlgebra("B2", 1)
                sage: bases = H._BasesCategory()
                sage: bases.super_categories()
                [Category of realizations of Iwahori-Hecke algebra of type B2 in 1,-1 over Integer Ring,
                 Category of finite dimensional algebras with basis over Integer Ring]
            """
            return [Realizations(self.base()), self.base()._category]

        def _repr_(self):
            r"""
            Return the representation of ``self``.

            EXAMPLES::

                sage: H = IwahoriHeckeAlgebra("B2", 1)
                sage: H._BasesCategory()
                Category of bases of Iwahori-Hecke algebra of type B2 in 1,-1 over Integer Ring
            """
            return "Category of bases of %s" % self.base()

        class ParentMethods:
            r"""
            This class collects code common to all the various bases. In most
            cases, these are just default implementations that will get
            specialized in a basis.
            """
            def _repr_(self):
                """
                Text representation of this basis of Iwahori-Hecke algebra.

                EXAMPLES::

                    sage: H = IwahoriHeckeAlgebra("B2", 1)
                    sage: H.T()
                    Iwahori-Hecke algebra of type B2 in 1,-1 over Integer Ring in the T-basis
                    sage: H.C()
                    Iwahori-Hecke algebra of type B2 in 1,-1 over Integer Ring in the C-basis
                    sage: H.Cp()
                    Iwahori-Hecke algebra of type B2 in 1,-1 over Integer Ring in the Cp-basis
                """
                return "%s in the %s-basis" % (self.realization_of(),
                                               self._basis_name)

            def __getitem__(self, i):
                """
                Return the basis element indexed by ``i``.

                INPUT:

                - ``i`` -- either an element of the Coxeter group or a
                  reduced word

                .. WARNING::

                    If ``i`` is not a reduced expression then the basis element
                    indexed by the corresponding element of the algebra is
                    returned rather than the corresponding product of the
                    generators::

                        sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
                        sage: T = IwahoriHeckeAlgebra('A3', v**2).T()
                        sage: T[1,1] == T[1] * T[1]
                        False

                EXAMPLES::

                    sage: H = IwahoriHeckeAlgebra("B2", 1)
                    sage: T = H.T()
                    sage: G = H.coxeter_group()
                    sage: T[G.one()]
                    1
                    sage: T[G.simple_reflection(1)]
                    T[1]
                    sage: T[G.from_reduced_word([1,2,1])]
                    T[1,2,1]
                    sage: T[[]]
                    1
                    sage: T[1]
                    T[1]
                    sage: T[1,2,1]
                    T[1,2,1]
                """
                W = self.realization_of().coxeter_group()
                if i in ZZ:
                    return self(W.simple_reflection(i))
                if i in W:
                    return self(i)
                if i == []:
                    return self.one()
                return self(W.from_reduced_word(i))

            def is_field(self, proof=True):
                """
                Return whether this Iwahori-Hecke algebra is a field.

                EXAMPLES::

                    sage: T = IwahoriHeckeAlgebra("B2", 1).T()
                    sage: T.is_field()
                    False
                """
                return False

            def is_commutative(self) -> bool:
                """
                Return whether this Iwahori-Hecke algebra is commutative.

                EXAMPLES::

                    sage: T = IwahoriHeckeAlgebra("B2", 1).T()
                    sage: T.is_commutative()
                    False
                """
                return self.base_ring().is_commutative() \
                    and self.realization_of().coxeter_group().is_commutative()

            @cached_method
            def one_basis(self):
                r"""
                Return the identity element in the Weyl group, as per
                ``AlgebrasWithBasis.ParentMethods.one_basis``.

                EXAMPLES::

                    sage: H = IwahoriHeckeAlgebra("B2", 1)
                    sage: H.T().one_basis()
                    [1 0]
                    [0 1]
                """
                return self.realization_of().coxeter_group().one()

            def index_set(self):
                r"""
                Return the index set of ``self``.

                EXAMPLES::

                    sage: IwahoriHeckeAlgebra("B2", 1).T().index_set()
                    (1, 2)
                """
                return self.realization_of().coxeter_group().index_set()

            @cached_method
            def algebra_generators(self):
                r"""
                Return the generators.

                They do not have order two but satisfy a quadratic relation.
                They coincide with the simple reflections in the Coxeter group
                when `q_1 = 1` and `q_2 = -1`. In this special case,
                the Iwahori-Hecke algebra is identified with the group algebra
                of the Coxeter group.

                EXAMPLES:

                In the standard basis::

                    sage: R.<q> = QQ[]
                    sage: H = IwahoriHeckeAlgebra("A3", q).T()
                    sage: T = H.algebra_generators(); T
                    Finite family {1: T[1], 2: T[2], 3: T[3]}
                    sage: T.list()
                    [T[1], T[2], T[3]]
                    sage: [T[i] for i in [1,2,3]]
                    [T[1], T[2], T[3]]
                    sage: T1,T2,T3 = H.algebra_generators()
                    sage: T1
                    T[1]
                    sage: H = IwahoriHeckeAlgebra(['A',2,1], q).T()
                    sage: T = H.algebra_generators(); T
                    Finite family {0: T[0], 1: T[1], 2: T[2]}
                    sage: T.list()
                    [T[0], T[1], T[2]]
                    sage: [T[i] for i in [0,1,2]]
                    [T[0], T[1], T[2]]
                    sage: [T0, T1, T2] = H.algebra_generators()
                    sage: T0
                    T[0]

                In the Kazhdan-Lusztig basis::

                    sage: R = LaurentPolynomialRing(QQ, 'v')
                    sage: v = R.gen(0)
                    sage: H = IwahoriHeckeAlgebra('A5', v**2)
                    sage: C = H.C()
                    sage: C.algebra_generators()
                    Finite family {1: C[1], 2: C[2], 3: C[3], 4: C[4], 5: C[5]}
                    sage: C.algebra_generators().list()
                    [C[1], C[2], C[3], C[4], C[5]]
                """
                return self.basis().keys().simple_reflections().map(self.monomial)

            def algebra_generator(self, i):
                r"""
                Return the `i`-th generator of ``self``.

                EXAMPLES:

                In the standard basis::

                    sage: R.<q>=QQ[]
                    sage: H = IwahoriHeckeAlgebra("A3", q).T()
                    sage: [H.algebra_generator(i) for i in H.index_set()]
                    [T[1], T[2], T[3]]

                In the Kazhdan-Lusztig basis::

                    sage: R = LaurentPolynomialRing(QQ, 'v')
                    sage: v = R.gen(0)
                    sage: H = IwahoriHeckeAlgebra('A5', v**2)
                    sage: C = H.C()
                    sage: [C.algebra_generator(i) for i in H.coxeter_group().index_set()]
                    [C[1], C[2], C[3], C[4], C[5]]
                """
                return self.algebra_generators()[i]

            @abstract_method(optional=True)
            def bar_on_basis(self, w):
                """
                Return the bar involution on the basis element of ``self``
                indexed by ``w``.

                EXAMPLES::

                    sage: R.<v> = LaurentPolynomialRing(QQ)
                    sage: H = IwahoriHeckeAlgebra('A3', v**2)
                    sage: W = H.coxeter_group()
                    sage: s1,s2,s3 = W.simple_reflections()
                    sage: Cp = H.Cp()
                    sage: Cp.bar_on_basis(s1*s2*s1*s3)
                    Cp[1,2,3,1]
                """

            @abstract_method(optional=True)
            def hash_involution_on_basis(self, w):
                """
                Return the bar involution on the basis element of ``self``
                indexed by ``w``.

                EXAMPLES::

                    sage: R.<v> = LaurentPolynomialRing(QQ)
                    sage: H = IwahoriHeckeAlgebra('A3', v**2)
                    sage: W = H.coxeter_group()
                    sage: s1,s2,s3 = W.simple_reflections()
                    sage: Cp = H.Cp()
                    sage: C = H.C()
                    sage: C(Cp.hash_involution_on_basis(s1*s2*s1*s3))
                    C[1,2,3,1]
                """

        class ElementMethods:
            def bar(self):
                r"""
                Return the bar involution of ``self``.

                The bar involution `\overline{\phantom{x}}` is an antilinear
                `\ZZ`-algebra involution defined by the identity on `\ZZ`,
                sending `q^{1/2} \mapsto q^{-1/2}`, and `\overline{T_w} =
                T_{w^{-1}}^{-1}`.

                REFERENCES:

                - :wikipedia:`Iwahori-Hecke_algebra#Canonical_basis`

                EXAMPLES:

                We first test on a single generator::

                    sage: R.<q> = LaurentPolynomialRing(QQ)
                    sage: H = IwahoriHeckeAlgebra('A3', q)
                    sage: T = H.T()
                    sage: T1,T2,T3 = T.algebra_generators()
                    sage: T1.bar()
                    (q^-1)*T[1] + (q^-1-1)
                    sage: T1.bar().bar() == T1
                    True

                Next on a multiple of generators::

                    sage: b = (T1*T2*T1).bar(); b
                    (q^-3)*T[1,2,1] + (q^-3-q^-2)*T[2,1] + (q^-3-q^-2)*T[1,2]
                     + (q^-3-2*q^-2+q^-1)*T[1] + (q^-3-2*q^-2+q^-1)*T[2]
                     + (q^-3-2*q^-2+2*q^-1-1)
                    sage: b.bar() == T1*T2*T1
                    True

                A sum::

                    sage: s = T1 + T2
                    sage: b = s.bar(); b
                    (q^-1)*T[1] + (q^-1)*T[2] + (2*q^-1-2)
                    sage: b.bar() == s
                    True

                A more complicated example::

                    sage: p = T1*T2 + (1-q+q^-1)*T3 - q^3*T1*T3
                    sage: p.bar()
                    -(q^-5)*T[3,1] + (q^-2)*T[1,2]
                     - (q^-5-q^-4-q^-2+q^-1)*T[1]
                     - (q^-5-q^-4+q^-2-q^-1-1)*T[3]
                     + (q^-2-q^-1)*T[2]
                     - (q^-5-2*q^-4+q^-3-1+q)
                    sage: p.bar().bar() == p
                    True

                This also works for arbitrary ``q1`` and ``q2``::

                    sage: R.<q1,q2> = LaurentPolynomialRing(QQ)
                    sage: H = IwahoriHeckeAlgebra('A3', q1, q2=-q2)
                    sage: T = H.T()
                    sage: T1,T2,T3 = T.algebra_generators()
                    sage: p = T1*T3 + T2
                    sage: p.bar()
                    (q1^-2*q2^-2)*T[3,1]
                     + (-q1^-1*q2^-2+q1^-2*q2^-1)*T[1]
                     + (-q1^-1*q2^-2+q1^-2*q2^-1)*T[3]
                     + (q1^-1*q2^-1)*T[2]
                     + (-q2^-1+q1^-1+q2^-2-2*q1^-1*q2^-1+q1^-2)
                    sage: p.bar().bar() == p
                    True

                Next we have an example in the `C` basis::

                    sage: R.<v> = LaurentPolynomialRing(QQ)
                    sage: H = IwahoriHeckeAlgebra('A3', v**2)
                    sage: C = H.C()
                    sage: p = C[1]*C[3] + C[2]
                    sage: p.bar()
                    C[3,1] + C[2]
                    sage: p.bar().bar() == p
                    True

                For the `C^{\prime}` basis as well::

                    sage: R.<v> = LaurentPolynomialRing(QQ)
                    sage: H = IwahoriHeckeAlgebra('A3', v**2)
                    sage: Cp = H.Cp()
                    sage: p = Cp[1]*Cp[3] + Cp[2]
                    sage: p.bar()
                    Cp[3,1] + Cp[2]

                TESTS:

                We check that doing the computations explicitly in the `T`
                basis gives the same results and with bar invariant
                coefficients::

                    sage: R.<v> = LaurentPolynomialRing(QQ)
                    sage: H = IwahoriHeckeAlgebra('A3', v**2)
                    sage: Cp = H.Cp()
                    sage: T = H.T()
                    sage: Cp(T(Cp[1,2,1])) == Cp[1,2,1]
                    True
                    sage: p = 4*Cp[1]*Cp[3] + (v^2 + v^-2 - 2)*Cp[2]
                    sage: Cp(T(p).bar()) == p
                    True
                """
                B = self.parent()
                if B.bar_on_basis is NotImplemented:
                    T = B.realization_of().T()
                    return B(T(self).bar())
                H = B.realization_of()
                return sum(H._bar_on_coefficients(c) * B.bar_on_basis(w)
                           for (w, c) in self)

            def hash_involution(self):
                r"""
                Return the hash involution of ``self``.

                The hash involution `\alpha` is a `\ZZ`-algebra
                involution of the Iwahori-Hecke algebra determined by
                `q^{1/2} \mapsto q^{-1/2}`, and `T_w \mapsto
                (-q_1 q_2)^{-\ell(w)} T_w`, for `w` an element of the
                corresponding Coxeter group.

                This map is defined in [KL1979]_ and it is used to
                change between the `C` and `C^{\prime}` bases because
                `\alpha(C_w) = (-1)^{\ell(w)} C_w'`.

                EXAMPLES::

                    sage: R.<v> = LaurentPolynomialRing(QQ)
                    sage: H = IwahoriHeckeAlgebra('A3', v**2)
                    sage: T = H.T()
                    sage: T1,T2,T3 = T.algebra_generators()
                    sage: elt = T1.hash_involution(); elt
                    -(v^-2)*T[1]
                    sage: elt.hash_involution()
                    T[1]
                    sage: elt = T1*T2 + (v^3 - v^-1 + 2)*T3*T1*T2*T3
                    sage: elt.hash_involution()
                    (v^-11+2*v^-8-v^-7)*T[1,2,3,2] + (v^-4)*T[1,2]
                    sage: elt.hash_involution().hash_involution() == elt
                    True

                With the Kazhdan-Lusztig `C^{\prime}` basis::

                    sage: Cp = H.Cp()
                    sage: p = Cp[1]*Cp[3] + Cp[2]
                    sage: q = p.hash_involution(); q
                    Cp[3,1] - (v^-1+v)*Cp[1] - (v^-1+v)*Cp[3] - Cp[2] + (v^-2+v^-1+2+v+v^2)
                    sage: q.hash_involution() == p
                    True

                With the Kazhdan-Lusztig `C` basis::

                    sage: C = H.C()
                    sage: p = C[1]*C[3] + C[2]
                    sage: q = p.hash_involution(); q
                    C[3,1] + (v^-1+v)*C[1] + (v^-1+v)*C[3] - C[2] + (v^-2-v^-1+2-v+v^2)
                    sage: q.hash_involution() == p
                    True
                """
                basis = self.parent()
                if basis.hash_involution_on_basis is NotImplemented:
                    T = basis.realization_of().T()
                    return basis(T(self).hash_involution())

                H = basis.realization_of()
                return basis(sum(H._bar_on_coefficients(c) * basis.hash_involution_on_basis(w) for (w, c) in self))

            def goldman_involution(self):
                r"""
                Return the Goldman involution of ``self``.

                The Goldman involution is the algebra involution of the
                Iwahori-Hecke algebra determined by

                .. MATH::

                    T_w \mapsto (-q_1 q_2)^{\ell(w)} T_{w^{-1}}^{-1},

                where `w` is
                an element of the corresponding Coxeter group. The main point
                here is that `q_1 q_2 T_s^{-1} = -T_s + q_1 + q_2`, for
                each simple reflection `s`.

                This map is defined in [Iwa1964]_. The *alternating Hecke algebra*
                is the fixed-point subalgebra the Iwahori-Hecke algebra under
                this involution.

                EXAMPLES::

                    sage: R.<v> = LaurentPolynomialRing(QQ)
                    sage: H = IwahoriHeckeAlgebra('A3', v**2)
                    sage: T = H.T()
                    sage: T[1].goldman_involution()
                    -T[1] - (1-v^2)
                    sage: T[2].goldman_involution()
                    -T[2] - (1-v^2)
                    sage: T[1,2].goldman_involution()
                    T[1,2] + (1-v^2)*T[1] + (1-v^2)*T[2] + (1-2*v^2+v^4)
                    sage: elt=T[1,2]+ v*T[1]
                    sage: elt.goldman_involution()
                    T[1,2] + (1-v-v^2)*T[1] + (1-v^2)*T[2] + (1-v-2*v^2+v^3+v^4)
                    sage: elt.goldman_involution().goldman_involution() == elt
                    True
                    sage: H.A()(elt).goldman_involution()==elt.goldman_involution()
                    True

                With different parameters::

                    sage: H = IwahoriHeckeAlgebra('A3', q1=v, q2=-v^-1)
                    sage: T = H.T()
                    sage: T[1].goldman_involution()
                    -T[1] - (v^-1-v)
                    sage: T[2].goldman_involution()
                    -T[2] - (v^-1-v)
                    sage: T[1,2].goldman_involution()
                    T[1,2] + (v^-1-v)*T[1] + (v^-1-v)*T[2] + (v^-2-2+v^2)
                    sage: elt=T[1,2]+ v*T[1]
                    sage: elt.goldman_involution()
                    T[1,2] + (v^-1-2*v)*T[1] + (v^-1-v)*T[2] + (v^-2-3+2*v^2)
                    sage: elt.goldman_involution().goldman_involution() == elt
                    True
                    sage: H.A()(elt).goldman_involution()==elt.goldman_involution()
                    True

                With the `A` basis::

                    sage: A = H.A()
                    sage: p = A[1,3] + A[2]
                    sage: q = p.goldman_involution(); q
                    A[3,1] - A[2]
                    sage: q.goldman_involution() == p
                    True

                TESTS::

                    sage: all(h.goldman_involution().goldman_involution() == h for h in T.basis()) # long time
                    True
                """
                basis = self.parent()
                if hasattr(basis, 'goldman_involution_on_basis'):
                    return basis.sum(c * basis.goldman_involution_on_basis(w) for (w, c) in self)

                T = basis.realization_of().T()
                return basis(T(self).goldman_involution())

            def specialize_to(self, new_hecke, num_vars=2):
                r"""
                Return the element in the Iwahori-Hecke algebra ``new_hecke``
                with respect to the same basis which is obtained from ``self``
                by specializing the generic parameters in this algebra to the
                parameters of ``new_hecke``.

                INPUT:

                - ``new_hecke`` -- the Hecke algebra in specialized parameters
                - ``num_vars`` -- the number of variables to specialize

                .. WARNING::

                    This is not always defined. In particular, the number of
                    generators must match ``num_vars``

                EXAMPLES::

                    sage: R.<a,b> = LaurentPolynomialRing(ZZ)
                    sage: H = IwahoriHeckeAlgebra("A3", a^2, -b^2)
                    sage: T = H.T()
                    sage: elt = T[1,2,1] + 3*T[1] - a*b*T[3]
                    sage: S.<q> = LaurentPolynomialRing(ZZ)
                    sage: HS = IwahoriHeckeAlgebra("A3", q^2, -1)
                    sage: selt = elt.specialize_to(HS); selt
                    T[1,2,1] + 3*T[1] + q^2*T[3]
                    sage: GA = IwahoriHeckeAlgebra("A3", 1, -1)
                    sage: elt.specialize_to(GA)
                    T[1,2,1] + 3*T[1] + T[3]

                We need to specify that we are only specializing
                one argument::

                    sage: selt.specialize_to(GA)
                    Traceback (most recent call last):
                    ...
                    TypeError: Wrong number of arguments
                    sage: selt.specialize_to(GA, 1)
                    T[1,2,1] + 3*T[1] + T[3]
                """
                q1 = new_hecke._q1
                q2 = new_hecke._q2
                new_basis = getattr(new_hecke, self.parent()._basis_name)()

                # is there an easier way than this to convert the
                # coefficients to the correct base ring for new_hecke?
                if num_vars == 2:
                    args = (q1, q2)
                elif num_vars == 1:
                    args = (q1,)
                else:
                    return new_basis._from_dict({w: new_hecke._base(c)
                                                 for w, c in self})

                return new_basis._from_dict({w: new_hecke._base(c(args))
                                             for w, c in self})

    class _Basis(CombinatorialFreeModule, BindableClass):
        r"""
        Technical methods (i.e., not mathematical) that are inherited by each
        basis of the algebra. These methods cannot be defined in the category.
        """
        def __init__(self, algebra, prefix=None):
            r"""
            Initialises a basis class for the Iwahori-Hecke algebra ``algebra``.
            Optionally, a ``prefix`` can be set which is used when printing the
            basis elements. The prefix defaults to ``self._basis-name``, which
            is the name of the basis class.

            EXAMPLES::

                sage: H = IwahoriHeckeAlgebra("G2",1)
                sage: t = H.T(prefix='t')
                sage: t[1]
                t[1]
            """
            if prefix is None:
                self._prefix = self._basis_name
            else:
                self._prefix = prefix

            CombinatorialFreeModule.__init__(self,
                                             algebra.base_ring(),
                                             algebra._W,
                                             category=algebra._BasesCategory(),
                                             sorting_key=sorting_key,
                                             prefix=self._prefix)

        # This **must** match the name of the class in order for
        #   specialize_to() to work
        _basis_name = 'B'

        def _repr_term(self, t):
            r"""
            Return the string representation of the term indexed by ``t``.

            EXAMPLES::

                sage: R.<q> = QQ[]
                sage: H = IwahoriHeckeAlgebra("A3", q)
                sage: W = H.coxeter_group()
                sage: H.T()._repr_term(W.from_reduced_word([1,2,3]))
                'T[1,2,3]'
            """
            redword = t.reduced_word()
            if not redword:
                return "1"
            return self._print_options['prefix'] + '[%s]' % ','.join('%d' % i for i in redword)

        def _latex_term(self, t):
            r"""
            Return latex for the term indexed by ``t``.

            EXAMPLES::

                sage: R.<v> = QQ[]
                sage: H = IwahoriHeckeAlgebra("A3", q1=v**2, q2=-1)
                sage: W = H.coxeter_group()
                sage: H.T()._latex_term(W.from_reduced_word([1,2,3]))
                'T_{1}T_{2}T_{3}'
            """
            redword = t.reduced_word()
            if not redword:
                return '1'
            return ''.join("%s_{%d}" % (self._print_options['prefix'], i)
                           for i in redword)

        def product_on_basis(self, w1, w2):
            r"""
            Return the product of the two basis elements indexed by ``w1`` and
            ``w2``.

            The computation is actually done by converting to the
            `T`-basis, multiplying and then converting back.

            EXAMPLES::

                sage: R = LaurentPolynomialRing(QQ, 'v')
                sage: v = R.gen(0)
                sage: H = IwahoriHeckeAlgebra('A2', v**2)
                sage: s1,s2 = H.coxeter_group().simple_reflections()
                sage: [H.Cp().product_on_basis(s1,x) for x in [s1,s2]]
                [(v^-1+v)*Cp[1], Cp[1,2]]
                sage: [H.C().product_on_basis(s1,x) for x in [s1,s2]]
                [-(v^-1+v)*C[1], C[1,2]]
            """
            return self(self.to_T_basis(w1) * self.to_T_basis(w2))

    class T(_Basis):
        r"""
        The standard basis of Iwahori-Hecke algebra.

        For every simple reflection `s_i` of the Coxeter group, there is a
        corresponding generator `T_i` of Iwahori-Hecke algebra. These
        are subject to the relations:

        .. MATH::

            (T_i - q_1) (T_i - q_2) = 0

        together with the braid relations:

        .. MATH::

            T_i T_j T_i \cdots = T_j T_i T_j \cdots,

        where the number of terms on each of the two sides is the order of
        `s_i s_j` in the Coxeter group.

        Weyl group elements form a basis of Iwahori-Hecke algebra `H`
        with the property that if `w_1` and `w_2` are Coxeter group elements
        such that `\ell(w_1 w_2) = \ell(w_1) + \ell(w_2)` then
        `T_{w_1 w_2} = T_{w_1} T_{w_2}`.

        With the default value `q_2 = -1` and with `q_1 = q` the
        generating relation may be written
        `T_i^2 = (q-1) \cdot T_i + q \cdot 1` as in [Iwa1964]_.

        EXAMPLES::

            sage: H = IwahoriHeckeAlgebra("A3", 1)
            sage: T = H.T()
            sage: T1,T2,T3 = T.algebra_generators()
            sage: T1*T2*T3*T1*T2*T1 == T3*T2*T1*T3*T2*T3
            True
            sage: w0 = T(H.coxeter_group().long_element())
            sage: w0
            T[1,2,3,1,2,1]
            sage: T = H.T(prefix='s')
            sage: T.an_element()
            s[1,2,3] + 2*s[1] + 3*s[2] + 1

        TESTS::

            sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
            sage: H = IwahoriHeckeAlgebra('A3', v**2)
            sage: W = H.coxeter_group()
            sage: T = H.T()
            sage: C = H.C()
            sage: Cp = H.Cp()
            sage: all(T(C(T[x])) == T[x] for x in W) # long time
            True
            sage: all(T(Cp(T[x])) == T[x] for x in W) # long time
            True

        We check a property of the bar involution and `R`-polynomials::

            sage: KL = KazhdanLusztigPolynomial(W, v)
            sage: all(T[x].bar() == sum(v^(-2*y.length()) * KL.R(y, x).substitute(v=v^-2) * T[y] for y in W) for x in W) # long time
            True
        """
        _basis_name = "T"   # this is used, for example, by specialize_to and is the default prefix

        def inverse_generator(self, i):
            r"""
            Return the inverse of the `i`-th generator, if it exists.

            This method is only available if the Iwahori-Hecke algebra
            parameters ``q1`` and ``q2`` are both invertible.  In this case,
            the algebra generators are also invertible and this method
            returns the inverse of ``self.algebra_generator(i)``.

            EXAMPLES::

                sage: P.<q1, q2>=QQ[]
                sage: F = Frac(P)
                sage: H = IwahoriHeckeAlgebra("A2", q1, q2=q2, base_ring=F).T()
                sage: H.base_ring()
                Fraction Field of Multivariate Polynomial Ring in q1, q2 over Rational Field
                sage: H.inverse_generator(1)
                -1/(q1*q2)*T[1] + ((q1+q2)/(q1*q2))
                sage: H = IwahoriHeckeAlgebra("A2", q1, base_ring=F).T()
                sage: H.inverse_generator(2)
                -(1/(-q1))*T[2] + ((q1-1)/(-q1))
                sage: P1.<r1, r2> = LaurentPolynomialRing(QQ)
                sage: H1 = IwahoriHeckeAlgebra("B2", r1, q2=r2, base_ring=P1).T()
                sage: H1.base_ring()
                Multivariate Laurent Polynomial Ring in r1, r2 over Rational Field
                sage: H1.inverse_generator(2)
                (-r1^-1*r2^-1)*T[2] + (r2^-1+r1^-1)
                sage: H2 = IwahoriHeckeAlgebra("C2", r1, base_ring=P1).T()
                sage: H2.inverse_generator(2)
                (r1^-1)*T[2] + (-1+r1^-1)
            """
            A = self.realization_of()
            try:
                # This currently works better than ~(self._q1) if
                # self.base_ring() is a Laurent polynomial ring since it
                # avoids accidental coercion into a field of fractions.
                i1 = normalized_laurent_polynomial(A._base, A._q1 ** -1)
                i2 = normalized_laurent_polynomial(A._base, A._q2 ** -1)
            except Exception:
                raise ValueError("%s and %s must be invertible" % (A._q1, A._q2))
            return (-i1*i2)*self.algebra_generator(i)+(i1+i2)

        @cached_method
        def inverse_generators(self):
            r"""
            Return the inverses of all the generators, if they exist.

            This method is only available if ``q1`` and ``q2`` are invertible.
            In that case, the algebra generators are also invertible.

            EXAMPLES::

                sage: P.<q> = PolynomialRing(QQ)
                sage: F = Frac(P)
                sage: H = IwahoriHeckeAlgebra("A2", q, base_ring=F).T()
                sage: T1,T2 = H.algebra_generators()
                sage: U1,U2 = H.inverse_generators()
                sage: U1*T1,T1*U1
                (1, 1)
                sage: P1.<q> = LaurentPolynomialRing(QQ)
                sage: H1 = IwahoriHeckeAlgebra("A2", q, base_ring=P1).T(prefix='V')
                sage: V1,V2 = H1.algebra_generators()
                sage: W1,W2 = H1.inverse_generators()
                sage: [W1,W2]
                [(q^-1)*V[1] + (q^-1-1), (q^-1)*V[2] + (q^-1-1)]
                sage: V1*W1, W2*V2
                (1, 1)
            """
            return Family(self.index_set(), self.inverse_generator)

        def product_on_basis(self, w1, w2):
            r"""
            Return `T_{w_1} T_{w_2}`, where `w_1` and `w_2` are words in the
            Coxeter group.

            EXAMPLES::

                sage: R.<q> = QQ[]; H = IwahoriHeckeAlgebra("A2", q)
                sage: T = H.T()
                sage: s1,s2 = H.coxeter_group().simple_reflections()
                sage: [T.product_on_basis(s1,x) for x in [s1,s2]]
                [(q-1)*T[1] + q, T[1,2]]
            """
            result = self.monomial(w1)
            for i in w2.reduced_word():
                result = self.product_by_generator(result, i)
            return result

        def product_by_generator_on_basis(self, w, i, side='right'):
            r"""
            Return the product `T_w T_i` (resp. `T_i T_w`) if ``side`` is
            ``'right'`` (resp. ``'left'``).

            If the quadratic relation is `(T_i-u)(T_i-v) = 0`, then we have

            .. MATH::

                T_w T_i = \begin{cases}
                T_{ws_i} & \text{if } \ell(ws_i) = \ell(w) + 1, \\
                (u+v) T_{ws_i} - uv T_w & \text{if } \ell(w s_i) = \ell(w)-1.
                \end{cases}

            The left action is similar.

            INPUT:

            - ``w`` -- an element of the Coxeter group
            - ``i`` -- an element of the index set
            - ``side`` -- ``'right'`` (default) or ``'left'``

            EXAMPLES::

                sage: R.<q> = QQ[]; H = IwahoriHeckeAlgebra("A2", q)
                sage: T = H.T()
                sage: s1,s2 = H.coxeter_group().simple_reflections()
                sage: [T.product_by_generator_on_basis(w, 1) for w in [s1,s2,s1*s2]]
                [(q-1)*T[1] + q, T[2,1], T[1,2,1]]
                sage: [T.product_by_generator_on_basis(w, 1, side='left') for w in [s1,s2,s1*s2]]
                [(q-1)*T[1] + q, T[1,2], (q-1)*T[1,2] + q*T[2]]
            """
            wi = w.apply_simple_reflection(i, side=side)
            A = self.realization_of()
            if w.has_descent(i, side=side):
                # 10% faster than a plain addition on the example of #12528
                return self.sum_of_terms(((w, A._q_sum), (wi, A._q_prod)),
                                         distinct=True)
            else:
                return self.monomial(wi)

        def product_by_generator(self, x, i, side='right'):
            r"""
            Return `T_i \cdot x`, where `T_i` is the `i`-th generator. This is
            coded individually for use in ``x._mul_()``.

            EXAMPLES::

                sage: R.<q> = QQ[]; H = IwahoriHeckeAlgebra("A2", q).T()
                sage: T1, T2 = H.algebra_generators()
                sage: [H.product_by_generator(x, 1) for x in [T1,T2]]
                [(q-1)*T[1] + q, T[2,1]]
                sage: [H.product_by_generator(x, 1, side = "left") for x in [T1,T2]]
                [(q-1)*T[1] + q, T[1,2]]
            """
            return self.linear_combination((self.product_by_generator_on_basis(w, i, side), c)
                                           for (w, c) in x)

        def to_C_basis(self, w):
            r"""
            Return `T_w` as a linear combination of `C`-basis elements.

            EXAMPLES::

                sage: R = LaurentPolynomialRing(QQ, 'v')
                sage: v = R.gen(0)
                sage: H = IwahoriHeckeAlgebra('A2', v**2)
                sage: s1,s2 = H.coxeter_group().simple_reflections()
                sage: T = H.T()
                sage: C = H.C()
                sage: T.to_C_basis(s1)
                v*T[1] + v^2
                sage: C(T(s1))
                v*C[1] + v^2
                sage: C(v^-1*T(s1) - v)
                C[1]
                sage: C(T(s1*s2)+T(s1)+T(s2)+1)
                v^2*C[1,2] + (v+v^3)*C[1] + (v+v^3)*C[2] + (1+2*v^2+v^4)
                sage: C(T(s1*s2*s1))
                v^3*C[1,2,1] + v^4*C[2,1] + v^4*C[1,2] + v^5*C[1] + v^5*C[2] + v^6
            """
            H = self.realization_of()
            generic_T = H._generic_iwahori_hecke_algebra.T()
            return generic_T.to_C_basis(w).specialize_to(H)

        def to_Cp_basis(self, w):
            r"""
            Return `T_w` as a linear combination of `C^{\prime}`-basis
            elements.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ)
                sage: H = IwahoriHeckeAlgebra('A2', v**2)
                sage: s1,s2 = H.coxeter_group().simple_reflections()
                sage: T = H.T()
                sage: Cp = H.Cp()
                sage: T.to_Cp_basis(s1)
                v*Cp[1] - 1
                sage: Cp(T(s1))
                v*Cp[1] - 1
                sage: Cp(T(s1)+1)
                v*Cp[1]
                sage: Cp(T(s1*s2)+T(s1)+T(s2)+1)
                v^2*Cp[1,2]
                sage: Cp(T(s1*s2*s1))
                v^3*Cp[1,2,1] - v^2*Cp[2,1] - v^2*Cp[1,2] + v*Cp[1] + v*Cp[2] - 1
            """
            H = self.realization_of()
            generic_T = H._generic_iwahori_hecke_algebra.T()
            return generic_T.to_Cp_basis(w).specialize_to(H)

        def bar_on_basis(self, w):
            """
            Return the bar involution of `T_w`, which is `T^{-1}_{w^-1}`.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ)
                sage: H = IwahoriHeckeAlgebra('A3', v**2)
                sage: W = H.coxeter_group()
                sage: s1,s2,s3 = W.simple_reflections()
                sage: T = H.T()
                sage: b = T.bar_on_basis(s1*s2*s3); b
                (v^-6)*T[1,2,3] + (v^-6-v^-4)*T[3,1]
                 + (v^-6-v^-4)*T[1,2] + (v^-6-v^-4)*T[2,3]
                 + (v^-6-2*v^-4+v^-2)*T[1] + (v^-6-2*v^-4+v^-2)*T[3]
                 + (v^-6-2*v^-4+v^-2)*T[2] + (v^-6-3*v^-4+3*v^-2-1)
                sage: b.bar()
                T[1,2,3]
            """
            return self.monomial(w.inverse()).inverse()

        def hash_involution_on_basis(self, w):
            r"""
            Return the hash involution on the basis element ``self[w]``.

            The hash involution `\alpha` is a `\ZZ`-algebra
            involution of the Iwahori-Hecke algebra determined by
            `q^{1/2} \mapsto q^{-1/2}`, and `T_w \mapsto
            (-q_1 q_2)^{-\ell(w)} T_w`, for `w` an element of the
            corresponding Coxeter group.

            This map is defined in [KL1979]_ and it is used to change between
            the `C` and `C^{\prime}` bases because
            `\alpha(C_w) = (-1)^{\ell(w)}C^{\prime}_w`.

            This function is not intended to be called directly. Instead, use
            :meth:`hash_involution`.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
                sage: H = IwahoriHeckeAlgebra('A3', v**2)
                sage: T=H.T()
                sage: s=H.coxeter_group().simple_reflection(1)
                sage: T.hash_involution_on_basis(s)
                -(v^-2)*T[1]
                sage: T[s].hash_involution()
                -(v^-2)*T[1]
                sage: h = T[1]*T[2] + (v^3 - v^-1 + 2)*T[3,1,2,3]
                sage: h.hash_involution()
                (v^-11+2*v^-8-v^-7)*T[1,2,3,2] + (v^-4)*T[1,2]
                sage: h.hash_involution().hash_involution() == h
                True
            """
            H = self.realization_of()
            return (-H._q_prod)**(-w.length())*self.monomial(w)

        def goldman_involution_on_basis(self, w):
            r"""
            Return the Goldman involution to the basis element
            indexed by ``w``.

            The goldman involution is the algebra involution of the
            Iwahori-Hecke algebra determined by

            .. MATH::

                T_w \mapsto (-q_1 q_2)^{\ell(w)} T_{w^{-1}}^{-1},

            where `w` is an element of the corresponding Coxeter group.

            This map is defined in [Iwa1964]_ and it is used to define the
            alternating subalgebra of the Iwahori-Hecke algebra, which is the
            fixed-point subalgebra of the Goldman involution.

            This function is not intended to be called directly. Instead, use
            :meth:`goldman_involution`.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
                sage: H = IwahoriHeckeAlgebra('A3', v**2)
                sage: T=H.T()
                sage: s=H.coxeter_group().simple_reflection(1)
                sage: T.goldman_involution_on_basis(s)
                -T[1] - (1-v^2)
                sage: T[s].goldman_involution()
                -T[1] - (1-v^2)
                sage: h = T[1]*T[2] + (v^3 - v^-1 + 2)*T[3,1,2,3]
                sage: h.goldman_involution()
                -(v^-1-2-v^3)*T[1,2,3,2]
                 - (v^-1-2-v+2*v^2-v^3+v^5)*T[3,1,2]
                 - (v^-1-2-v+2*v^2-v^3+v^5)*T[1,2,3]
                 - (v^-1-2-v+2*v^2-v^3+v^5)*T[2,3,2]
                 - (v^-1-2-2*v+4*v^2-2*v^4+2*v^5-v^7)*T[3,1]
                 - (v^-1-3-2*v+4*v^2-2*v^4+2*v^5-v^7)*T[1,2]
                 - (v^-1-2-2*v+4*v^2-2*v^4+2*v^5-v^7)*T[3,2]
                 - (v^-1-2-2*v+4*v^2-2*v^4+2*v^5-v^7)*T[2,3]
                - (v^-1-3-2*v+5*v^2+v^3-4*v^4+v^5+2*v^6-2*v^7+v^9)*T[1]
                 - (v^-1-2-3*v+6*v^2+2*v^3-6*v^4+2*v^5+2*v^6-3*v^7+v^9)*T[3]
                 - (v^-1-3-3*v+7*v^2+2*v^3-6*v^4+2*v^5+2*v^6-3*v^7+v^9)*T[2]
                 - (v^-1-3-3*v+8*v^2+3*v^3-9*v^4+6*v^6-3*v^7-2*v^8+3*v^9-v^11)
                sage: h.goldman_involution().goldman_involution() == h
                True
            """
            H = self.realization_of()
            return (-H._q_prod)**w.length() * self.monomial(w.inverse()).inverse()

        class Element(CombinatorialFreeModule.Element):
            r"""
            A class for elements of an Iwahori-Hecke algebra in the `T` basis.

            TESTS::

                sage: R.<q> = QQ[]
                sage: H = IwahoriHeckeAlgebra("B3",q).T()
                sage: T1,T2,T3 = H.algebra_generators()
                sage: T1+2*T2*T3
                2*T[2,3] + T[1]
                sage: T1*T1
                (q-1)*T[1] + q

                sage: R.<q1,q2> = QQ[]
                sage: H = IwahoriHeckeAlgebra("A2", q1, q2=q2).T(prefix='x')
                sage: sum(H.algebra_generators())^2
                x[2,1] + x[1,2] + (q1+q2)*x[1] + (q1+q2)*x[2] + (-2*q1*q2)

                sage: H = IwahoriHeckeAlgebra("A2", q1, q2=q2).T(prefix='t')
                sage: t1,t2 = H.algebra_generators()
                sage: (t1-t2)^3
                (q1^2-q1*q2+q2^2)*t[1] + (-q1^2+q1*q2-q2^2)*t[2]

                sage: R.<q> = QQ[]
                sage: H = IwahoriHeckeAlgebra("G2", q).T()
                sage: [T1, T2] = H.algebra_generators()
                sage: T1*T2*T1*T2*T1*T2 == T2*T1*T2*T1*T2*T1
                True
                sage: T1*T2*T1 == T2*T1*T2
                False

                sage: H = IwahoriHeckeAlgebra("A2", 1).T()
                sage: [T1,T2] = H.algebra_generators()
                sage: T1+T2
                T[1] + T[2]

                sage: -(T1+T2)
                -T[1] - T[2]
                sage: 1-T1
                -T[1] + 1

                sage: T1.parent()
                Iwahori-Hecke algebra of type A2 in 1,-1 over Integer Ring in the T-basis
            """
            def __invert__(self):
                r"""
                Return the inverse if ``self`` is a basis element.

                An element is a basis element if it is `T_w` where `w` is in
                the Weyl group. The base ring must be a field or Laurent
                polynomial ring. Other elements of the ring have inverses but
                the inverse method is only implemented for the basis elements.

                EXAMPLES::

                    sage: R.<q> = LaurentPolynomialRing(QQ)
                    sage: H = IwahoriHeckeAlgebra("A2", q).T()
                    sage: [T1,T2] = H.algebra_generators()
                    sage: x = (T1*T2).inverse(); x   # indirect doctest
                    (q^-2)*T[2,1] + (q^-2-q^-1)*T[1] + (q^-2-q^-1)*T[2] + (q^-2-2*q^-1+1)
                    sage: x*T1*T2
                    1

                TESTS:

                We check some alternative forms of input for inverting
                an element::

                    sage: R.<q> = LaurentPolynomialRing(QQ)
                    sage: H = IwahoriHeckeAlgebra("A2", q).T()
                    sage: T1,T2 = H.algebra_generators()
                    sage: ~(T1*T2)
                    (q^-2)*T[2,1] + (q^-2-q^-1)*T[1] + (q^-2-q^-1)*T[2] + (q^-2-2*q^-1+1)
                    sage: (T1*T2)^(-1)
                    (q^-2)*T[2,1] + (q^-2-q^-1)*T[1] + (q^-2-q^-1)*T[2] + (q^-2-2*q^-1+1)
                """
                if len(self) != 1:
                    raise NotImplementedError("inverse only implemented for basis elements (monomials in the generators)" % self)
                H = self.parent()
                w = self.support_of_term()

                return H.prod(H.inverse_generator(i) for i in reversed(w.reduced_word()))

    standard = T

    class _KLHeckeBasis(_Basis):
        r"""
        Abstract class for the common methods for the Kazhdan-Lusztig `C` and
        `C^{\prime}` bases.
        """
        def __init__(self, IHAlgebra, prefix=None):
            r"""
            Initialize the Kazhdan-Lusztig basis of the Iwahori-Hecke
            algebra ``IHAlgebra``.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ)
                sage: H = IwahoriHeckeAlgebra('A3', v**2)
                sage: Cp = H.Cp()
                sage: C = H.C()
            """
            if IHAlgebra._root is None:
                raise ValueError('The Kazhdan_Lusztig bases are defined '
                                 'only when -q_1*q_2 is a square')

            if IHAlgebra._is_generic:
                klbasis = IwahoriHeckeAlgebra_nonstandard._KLHeckeBasis
            else:
                klbasis = IwahoriHeckeAlgebra._KLHeckeBasis
            super(klbasis, self).__init__(IHAlgebra, prefix)

            # Define conversion from the KL-basis to the T-basis via
            # specialization from the generic Hecke algebra
            self.module_morphism(self.to_T_basis, codomain=IHAlgebra.T(), category=self.category()
                                 ).register_as_coercion()

            # ...and from the T_basis to the KL-basis.
            T = IHAlgebra.T()
            T.module_morphism(getattr(T, 'to_{}_basis'.format(self._basis_name)),
                              codomain=self, category=self.category()
                              ).register_as_coercion()

        def bar_on_basis(self, w):
            r"""
            Return the bar involution on the Kazhdan-Lusztig basis element
            indexed by ``w``. By definition, all Kazhdan-Lusztig basis elements
            are fixed by the bar involution.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ)
                sage: H = IwahoriHeckeAlgebra('A3', v**2)
                sage: W = H.coxeter_group()
                sage: s1,s2,s3 = W.simple_reflections()
                sage: Cp = H.Cp()
                sage: Cp.bar_on_basis(s1*s2*s1*s3)
                Cp[1,2,3,1]
            """
            return self.monomial(w)

        def to_T_basis(self, w):
            r"""
            Return the Kazhdan-Lusztig basis element ``self[w]`` as a linear
            combination of ``T``-basis elements.

            EXAMPLES::

                sage: H=IwahoriHeckeAlgebra("A3",1); Cp=H.Cp(); C=H.C()
                sage: s=H.coxeter_group().simple_reflection(1)
                sage: C.to_T_basis(s)
                T[1] - 1
                sage: Cp.to_T_basis(s)
                T[1] + 1
            """
            H = self.realization_of()
            generic_KL = getattr(H._generic_iwahori_hecke_algebra, self._basis_name)()
            return generic_KL.to_T_basis(w).specialize_to(H)

    class Cp(_KLHeckeBasis):
        r"""
        The `C^{\prime}` Kazhdan-Lusztig basis of Iwahori-Hecke algebra.

        Assuming the standard quadratic relations of `(T_r-q)(T_r+1)=0`, for
        every element `w` in the Coxeter group, there is a unique element
        `C^{\prime}_w` in the Iwahori-Hecke algebra which is uniquely determined
        by the two properties:

        .. MATH::

            \begin{aligned}
                \overline{ C^{\prime}_w } &= C^{\prime}_w, \\
                C^{\prime}_w &= q^{-\ell(w)/2}
                    \sum_{v \leq w} P_{v,w}(q) T_v,
            \end{aligned}

        where `\leq` is the Bruhat order on the underlying Coxeter group and
        `P_{v,w}(q) \in \ZZ[q,q^{-1}]` are polynomials in `\ZZ[q]` such that
        `P_{w,w}(q) = 1` and if `v < w` then `\deg P_{v,w}(q) \leq
        \frac{1}{2}(\ell(w)-\ell(v)-1)`.

        More generally, if the quadratic relations are of the form
        (T_s-q_1)(T_s-q_2)=0` and `\sqrt{-q_1q_2}` exists then, for a simple
        reflection `s`, the corresponding Kazhdan-Lusztig basis element is:

        .. MATH::

            C^{\prime}_s = (-q_1 q_2)^{-1/2} (T_s + 1).

        See [KL1979]_ for more details.

        If the optional ``coxeter3`` package is available and the
        Iwahori--Hecke algebra was initialized in the "standard" presentation
        where `\{q_1,q_2\} = \{v^2,1\}` as sets or the "normalized"
        presentation where `\{q_1,q_2\} = \{v,-v^{-1}\}` as sets, the function
        :func:`product_on_basis` in this class computes products in the
        `C^{\prime}`-basis directly in the basis itself, using ``coxeter3`` to
        calculate certain `\mu`-coefficients quickly. If the above conditions
        are not all met, the function computes such products indirectly, by
        converting elements to the `T`-basis, computing products there, and
        converting back. The indirect method can be prohibitively slow for
        more complex calculations; the direct method is faster.

        EXAMPLES::

            sage: R = LaurentPolynomialRing(QQ, 'v')
            sage: v = R.gen(0)
            sage: H = IwahoriHeckeAlgebra('A5', v**2)
            sage: W = H.coxeter_group()
            sage: s1,s2,s3,s4,s5 = W.simple_reflections()
            sage: T = H.T()
            sage: Cp = H.Cp()
            sage: T(s1)**2
            -(1-v^2)*T[1] + v^2
            sage: T(Cp(s1))
            (v^-1)*T[1] + (v^-1)
            sage: T(Cp(s1)*Cp(s2)*Cp(s1))
            (v^-3)*T[1,2,1] + (v^-3)*T[2,1] + (v^-3)*T[1,2]
             + (v^-3+v^-1)*T[1] + (v^-3)*T[2] + (v^-3+v^-1)

        ::

            sage: R = LaurentPolynomialRing(QQ, 'v')
            sage: v = R.gen(0)
            sage: H = IwahoriHeckeAlgebra('A3', v**2)
            sage: W = H.coxeter_group()
            sage: s1,s2,s3 = W.simple_reflections()
            sage: Cp = H.Cp()
            sage: Cp(s1*s2*s1)
            Cp[1,2,1]
            sage: Cp(s1)**2
            (v^-1+v)*Cp[1]
            sage: Cp(s1)*Cp(s2)*Cp(s1)
            Cp[1,2,1] + Cp[1]
            sage: Cp(s1)*Cp(s2)*Cp(s3)*Cp(s1)*Cp(s2)    # long time
            Cp[1,2,3,1,2] + Cp[1,2,1] + Cp[3,1,2]

        In the following product computations, whether ``coxeter3`` is
        installed makes a big difference: without ``coxeter3`` the product in
        type `H_4` takes about 5 seconds to compute and the product in type
        `A_9` seems infeasible, while with ``coxeter3`` both the computations
        are instant::

            sage: H = IwahoriHeckeAlgebra('H4', v**2)   # optional - coxeter3
            sage: Cp = H.Cp()                           # optional - coxeter3
            sage: Cp[3,4,3]*Cp[3,4,3,4]*Cp[1,2,3,4]     # optional - coxeter3
            (v^-2+2+v^2)*Cp[3,4,3,4,1,2,3,4,2]
            + (v^-2+2+v^2)*Cp[3,4,3,4,3,1,2]
            + (v^-3+3*v^-1+3*v+v^3)*Cp[3,4,3,4,3,1]
            + (v^-1+v)*Cp[3,4,1,2,3,4]
            + (v^-1+v)*Cp[3,4,1,2]

            sage: H = IwahoriHeckeAlgebra('A9', v**2)   # optional - coxeter3
            sage: Cp = H.Cp()                           # optional - coxeter3
            sage: Cp[1,2,1,8,9,8]*Cp[1,2,3,7,8,9]       # optional - coxeter3
            (v^-2+2+v^2)*Cp[7,8,9,7,8,7,1,2,3,1]
            + (v^-2+2+v^2)*Cp[8,9,8,7,1,2,3,1]
            + (v^-3+3*v^-1+3*v+v^3)*Cp[8,9,8,1,2,3,1]

        To use ``coxeter3`` for product computations most efficiently, we
        recommend creating the Iwahori-Hecke algebra from a Coxeter group
        implemented with ``coxeter3`` to avoid unnecessary conversions, as in
        the following example with the same product computed in the last one::

            sage: # optional - coxeter3
            sage: R = LaurentPolynomialRing(QQ, 'v')
            sage: v = R.gen(0)
            sage: W = CoxeterGroup('A9', implementation='coxeter3')
            sage: H = IwahoriHeckeAlgebra(W, v**2)
            sage: Cp = H.Cp()
            sage: Cp[1,2,1,8,9,8]*Cp[1,2,3,7,8,9]
            (v^-2+2+v^2)*Cp[1,2,1,3,7,8,7,9,8,7]
            + (v^-2+2+v^2)*Cp[1,2,1,3,8,9,8,7]
            + (v^-3+3*v^-1+3*v+v^3)*Cp[1,2,1,3,8,9,8]

        TESTS::

            sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
            sage: H = IwahoriHeckeAlgebra('A3', v**2)
            sage: W = H.coxeter_group()
            sage: T = H.T()
            sage: C = H.C()
            sage: Cp = H.Cp()
            sage: all(Cp(T(Cp[x])) == Cp[x] for x in W) # long time
            True
            sage: all(Cp(C(Cp[x])) == Cp[x] for x in W) # long time
            True
        """
        _basis_name = 'Cp'   # this is used, for example, by specialize_to and is the default prefix

        def __init__(self, IHAlgebra, prefix=None):
            r"""
            TESTS::

                sage: # optional - coxeter3
                sage: R.<v> = LaurentPolynomialRing(QQ)
                sage: W = CoxeterGroup('A3', implementation='coxeter3')
                sage: H = IwahoriHeckeAlgebra(W, v**2)
                sage: Cp = H.Cp()
                sage: Cp._delta == v + ~v
                True
                sage: Cp._W_Coxeter3 == H._W
                True
                sage: H = IwahoriHeckeAlgebra(W, QQ(1))
                sage: Cp = H.Cp()
                sage: Cp._W_Coxeter3 is None
                True
            """
            super().__init__(IHAlgebra, prefix)

            self._W_Coxeter3 = None

            # See if we meet the conditions to use the direct product_on_basis algorithm.
            # To use v + ~v as the value delta, we need the standard or
            #   normalized presentations of the Hecke algebra.
            v = IHAlgebra.base_ring().gen(0)
            parameters = {IHAlgebra.q1(), IHAlgebra.q2()}
            if v == IHAlgebra.base_ring().one() or (parameters != {v**2, -1} and parameters != {v, -1/v}):
                return

            # check if products can be computed directly using ``coxeter3``
            try:
                from sage.libs.coxeter3.coxeter_group import CoxeterGroup as Coxeter3Group
            except ImportError:
                return

            self._delta = v + ~v
            if isinstance(IHAlgebra._W, Coxeter3Group):
                self._W_Coxeter3 = IHAlgebra._W
            else:
                self._W_Coxeter3 = CoxeterGroup(IHAlgebra._W.coxeter_type(), implementation='coxeter3')

        def hash_involution_on_basis(self, w):
            r"""
            Return the effect of applying the hash involution to the basis
            element ``self[w]``.

            This function is not intended to be called directly. Instead, use
            :meth:`hash_involution`.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
                sage: H = IwahoriHeckeAlgebra('A3', v**2)
                sage: Cp=H.Cp()
                sage: s=H.coxeter_group().simple_reflection(1)
                sage: Cp.hash_involution_on_basis(s)
                -Cp[1] + (v^-1+v)
                sage: Cp[s].hash_involution()
                -Cp[1] + (v^-1+v)
            """
            return (-1)**w.length() * self(self.realization_of().C().monomial(w))

        def product_on_basis(self, w1, w2):
            r"""
            Return the expansion of `C^{\prime}_{w_1} \cdot C^{\prime}_{w_2}`
            in the `C^{\prime}`-basis.

            If ``coxeter3`` is installed and the Iwahori--Hecke algebra is in
            the standard or normalized presentation, the product is computed
            directly using the method described in ALGORITHM. If not, the
            product is computed indirectly by converting the factors to the
            `T`-basis, computing the product there, and converting back.

            The following formulas for products of the forms `C^{\prime}_s
            \cdot C^{\prime}_w`  and `C^{\prime}_w \cdot C^{\prime}_s`, where
            `s` is a generator of the Coxeter group and `w` an arbitrary
            element, are key to the direct computation method. The formulas are
            valid for both the standard and normalized presentation of the
            Hecke algebra.

            .. MATH::

                C^{\prime}_s \cdot C^{\prime}_w = \begin{cases}
                  (q+q^{-1})C^{\prime}_{w},   & \text{if } \ell(sw) = \ell(w)-1,\\
                  C^{\prime}_{sw}+\sum_{v\leq w, sv \leq v} \mu(v,w)C^{\prime}_v,
                                              & \text{if } \ell(sw) = \ell(w)+1.
                \end{cases}

                \qquad\qquad

                C^{\prime}_w \cdot C^{\prime}_s = \begin{cases}
                  (q+q^{-1})C^{\prime}_{w},   & \text{if } \ell(ws) = \ell(w)-1,\\
                  C^{\prime}_{ws}+\sum_{v\leq w, vs \leq v} \mu(v,w)C^{\prime}_v,
                                              & \text{if } \ell(ws) = \ell(w)+1.
                \end{cases}

            In the above, `\leq` is the Bruhat order on the Coxeter group and
            `\mu(v,w)` is the "leading coefficient of Kazhdan-Lusztig
            polynomials"; see [KL1979]_ and [Lus2013]_ for more details. The
            method designates the computation of the `\mu`-coefficients to
            Sage's interface to Fokko du Cloux's ``coxeter3`` package, which is
            why the method requires the creation of the Coxeter group using the
            ``'coxeter3'`` implementation.

            ALGORITHM:

            The direct algorithm for computing `C^{\prime}_x \cdot
            C^{\prime}_y` runs in two steps as follows.

            If `\ell(x) \leq \ell(y)`, we first decompose `C^{\prime}_x` into
            a polynomial in the generators `C^{\prime}_s (s\in S)` and then
            multiply that polynomial with `C^{\prime}_y`. If `\ell(x) >
            \ell(y)`, we decompose `C^{\prime}_y` into a polynomial in
            `C^{\prime}_s (s\in S)` and multiply that polynomial with
            `C^{\prime}_x`. The second step (multiplication) is done by
            repeatedly applying the formulas displayed earlier directly. The
            first step (decomposition) is done by induction on the Bruhat order
            as follows: for every element `u\in W` with length `\ell(u)>1`,
            pick a left descent `s` of `u` and write `u=sw` (so `w=su`), then
            note that

            .. MATH::

                C^{\prime}_u = C^{\prime}_s \cdot C^{\prime}_{w}
                               - \sum_{v \le u; sv < v} \mu(v,w) C^{\prime}_v

            by the earlier formulas, where the element `w` and all elements
            `v`'s on the right side are lower than `u` in the Bruhat order;
            this allows us to finish the computation by decomposing the lower
            order terms `C^{\prime}_w` and each `C^{\prime}_v`. For example,
            for `u=121, s=1, w=21` in type `A_3` we have `C^{\prime}_{121} =
            C^{\prime}_1  C^{\prime}_{21} - C^{\prime}_1`, where the lower
            order term `C^{\prime}_{21}` further decomposes into `C^{\prime}_2
            C^{\prime}_1`, therefore

            .. MATH::

                C^{\prime}_{121} = C^{\prime}_1 C^{\prime}_2 C^{\prime}_1
                                   - C^{\prime}_1.

            We note that the base cases `\ell(x)=1` or `\ell(x)=0` of the above
            induction occur when `x` is itself a Coxeter generator `s` or the
            group identity, respectively. The decomposition is trivial in these
            cases (we have `C^{\prime}_x=C^{\prime}_s` or `C^{\prime}_x=1`, the
            unit of the Hecke algebra).

            EXAMPLES::

                sage: # optional - coxeter3
                sage: R.<v> = LaurentPolynomialRing(ZZ, 'v')
                sage: W = CoxeterGroup('A3', implementation='coxeter3')
                sage: H = IwahoriHeckeAlgebra(W, v**2); Cp=H.Cp()
                sage: Cp.product_on_basis(W([1,2,1]), W([3,1]))
                (v^-1+v)*Cp[1,2,1,3]
                sage: Cp.product_on_basis(W([1,2,1]), W([3,1,2]))
                (v^-1+v)*Cp[1,2,1,3,2] + (v^-1+v)*Cp[1,2,1]
            """
            if self._W_Coxeter3 is None:
                # We do not meet the conditions to use the direct product
                # algorithm; fall back to conversion to/from the T-basis.
                return super().product_on_basis(w1, w2)

            # If self._W_Coxeter3 is not the underlying Coxeter group, we need
            # to convert elements first for this algorithm.
            if self._W_Coxeter3 != self.realization_of()._W:
                w1 = self._W_Coxeter3.from_reduced_word(w1.reduced_word())
                w2 = self._W_Coxeter3.from_reduced_word(w2.reduced_word())

            # Decomposition: write one of C'_{w1} and C'_{w2} as a polynomial in the
            # generators C'_{s}.
            if len(w1) <= len(w2):
                side = 'left'
                gen_expression = self._decompose_into_generators(w1)
                other_element = self.monomial(w2)
            else:
                side = 'right'
                gen_expression = self._decompose_into_generators(w2)
                other_element = self.monomial(w1)

            # Multiplication: multiply the generators in each term of the above
            # polynomial onto other_element and add that summand onto result.
            result = self.zero()
            for (p, coeff) in gen_expression.items():
                summand = coeff * other_element
                if side == 'right':
                    for s in p:
                        summand = self._product_with_generator(s, summand, side)
                else:
                    for s in reversed(p):
                        summand = self._product_with_generator(s, summand, side)
                result += summand

            # Again, if self._W_Coxeter3 is not the underlying Coxeter group,
            # we need to convert the result. Specifically, make sure basis
            # elements appearing therein are actually indexed by elements of
            # the original underlying Coxeter group.
            if self._W_Coxeter3 != self.realization_of()._W:
                _W = self.realization_of()._W
                result = self._from_dict({_W.from_reduced_word(w.reduced_word()): c
                                          for (w, c) in result}, remove_zeros=False)

            return result

        def _product_with_generator_on_basis(self, s, w, side='left'):
            r"""
            Compute the product of `C^{\prime}_s` and `C^{\prime}_w`, putting
            `C^{\prime}_s` on the given ``side``.

            INPUT:

            - ``s`` -- integer in ``self.index_set()``
            - ``w`` -- a word in ``self.coxeter_group()``
            - ``side`` -- string; ``'left'`` or ``'right'``

            EXAMPLES::

                sage: # optional - coxeter3
                sage: R.<v> = LaurentPolynomialRing(ZZ, 'v')
                sage: W = CoxeterGroup('A3', implementation='coxeter3')
                sage: H = IwahoriHeckeAlgebra(W, v**2); Cp = H.Cp()
                sage: Cp._product_with_generator_on_basis(1, W([2,1]), 'left')
                Cp[1,2,1] + Cp[1]
                sage: Cp._product_with_generator_on_basis(1, W([2,1]), 'right')
                (v^-1+v)*Cp[2,1]
                sage: Cp._product_with_generator_on_basis(2, W([1,3,2,1,3]), 'right')
                Cp[1,2,1,3,2,1] + Cp[1,2,3,2] + Cp[1,3,2,1]
            """
            # use the product formula described in the class' documentation
            if w.has_descent(s, side=side):
                return self._from_dict({w: self._delta}, remove_zeros=False)

            element = {}
            between = self._W_Coxeter3.bruhat_interval([], w)
            R = self.base_ring()
            for x in between:
                # Get (coxeter3-implemented) group element corresponding to x
                x_elt = self._W_Coxeter3(x)
                if x_elt.has_descent(s, side=side):
                    # Compute mu-coefficient via coxeter3
                    coeff = R(x.mu_coefficient(w))
                    if coeff:
                        element[x_elt] = coeff
            longer_word = self._W_Coxeter3([s]) * w if side == 'left' else w * self._W_Coxeter3([s])
            element[longer_word] = R.one()
            return self._from_dict(element, remove_zeros=False)

        def _product_with_generator(self, s, x, side='left'):
            r"""
            Compute the product of `C^{\prime}_s` with any linear
            combination of `C^{\prime}`-basis elements.

            INPUT:

            - ``s`` -- integer in ``self.index_set()``
            - ``x`` -- any element of ``self``
            - ``side`` -- string; ``'left'`` or ``'right'``

            EXAMPLES::

                sage: # optional - coxeter3
                sage: R.<v> = LaurentPolynomialRing(ZZ, 'v')
                sage: W = CoxeterGroup('A3', implementation='coxeter3')
                sage: H = IwahoriHeckeAlgebra(W, v**2); Cp = H.Cp()
                sage: Cp._product_with_generator(1, Cp[1]+Cp[2], 'left')
                Cp[1,2] + (v^-1+v)*Cp[1]
                sage: Cp._product_with_generator(1, Cp[1]+Cp[2], 'right')
                Cp[2,1] + (v^-1+v)*Cp[1]
            """
            return self.linear_combination((self._product_with_generator_on_basis(s, w, side), coeff) for (w, coeff) in x)

        def _decompose_into_generators(self, u):
            r"""
            Decompose `C^{\prime}_u` into a polynomial in the KL generators
            `C^{\prime}_s`; see the ALGORITHM section of
            :func:`product_on_basis`.

            OUTPUT:

            A dictionary keyed by tuples with integer values. Each entry
            represents a term, where the tuple represents a monomial term
            in the KL generators and the value represents the coefficient
            of that term. For example, an item `(1,2): 3` stands for
            `3 \cdot C^{\prime}_1 C^{\prime}_2`.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(ZZ, 'v')             # optional - coxeter3
                sage: W = CoxeterGroup('A3', implementation='coxeter3')  # optional - coxeter3
                sage: H = IwahoriHeckeAlgebra(W, v**2); Cp=H.Cp()        # optional - coxeter3

            When `u` is itself a generator `s`, the decomposition is trivial::

                sage: Cp._decompose_into_generators(W([1]))              # optional - coxeter3
                {(1,): 1}

            Another example, where `C^{\prime}_u` happens to be a monomial
            (e.g., `C'_{21}  = C'_2 C'_1`)::

                sage: Cp._decompose_into_generators(W([2,1]))            # optional - coxeter3
                {(2, 1): 1}

            In more general situations the sum is a polynomial (e.g.,
            `C'_{121} = C'_1 C'_2 C'_1 - C'_1)`::

                sage: Cp._decompose_into_generators(W([1,2,1]))          # optional - coxeter3
                {(1,): -1, (1, 2, 1): 1}
                sage: Cp._decompose_into_generators(W([1,2,3,1,2]))      # optional - coxeter3
                {(1,): 1, (1, 2, 1): -1, (1, 2, 1, 3, 2): 1, (1, 3, 2): -1}

            TESTS::

                sage: Cp._decompose_into_generators(W([]))               # optional - coxeter3
                {(): 1}
            """
            # l(y) = 0 or 1
            if not len(u):
                return {(): 1}
            if len(u) == 1:
                return {(u[0],): 1}

            # l(y) > 1, use the recursive method described in product_on_basis
            s = u[0]
            w = u[1:]  # so CpC_s * CpC_w = CpC_u + lower order terms

            # get the lower order terms ("sum_term")
            sum_term = {}
            between = self._W_Coxeter3.bruhat_interval([], w)
            R = self.base_ring()
            for v in between:
                # Get (coxeter3-implemented) group element corresponding to v
                v_elt = self._W_Coxeter3(v)
                if v_elt.has_left_descent(s):
                    # Compute mu-coefficient via coxeter3
                    coeff = R(v.mu_coefficient(w))
                    if coeff:
                        sum_term[v_elt] = coeff

            # recursion: decompose C'_s * C'_w and the lower order terms
            result = {(s,) + gens: coeff for (gens, coeff) in self._decompose_into_generators(w).items()}
            zero = R.zero()
            for (z, c1) in sum_term.items():
                # Subtract off each term from sum_term.
                for (gens, c2) in self._decompose_into_generators(z).items():
                    result[gens] = result.get(gens, zero) - c1 * c2

            return result

    C_prime = Cp

    class C(_KLHeckeBasis):
        r"""
        The Kazhdan-Lusztig `C`-basis of Iwahori-Hecke algebra.

        Assuming the standard quadratic relations of `(T_r-q)(T_r+1)=0`, for
        every element `w` in the Coxeter group, there is a unique element
        `C_w` in the Iwahori-Hecke algebra which is uniquely determined
        by the two properties:

        .. MATH::

            \begin{aligned}
                \overline{C_w} &= C_w \\
                C_w &= (-1)^{\ell(w)} q^{\ell(w)/2}
                    \sum_{v \leq w} (-q)^{-\ell(v)}\overline{P_{v,w}(q)} T_v
            \end{aligned}

        where `\leq` is the Bruhat order on the underlying Coxeter group and
        `P_{v,w}(q)\in\ZZ[q,q^{-1}]` are polynomials in `\ZZ[q]` such that
        `P_{w,w}(q) = 1` and if `v < w` then
        `\deg P_{v,w}(q) \leq \frac{1}{2}(\ell(w) - \ell(v) - 1)`.
        This is related to the `C^{\prime}` Kazhdan-Lusztig basis by `C_i =
        -\alpha(C_i^{\prime})` where `\alpha` is the `\ZZ`-linear Hecke
        involution defined by `q^{1/2} \mapsto q^{-1/2}` and `\alpha(T_i) =
        -(q_1 q_2)^{-1/2} T_i`.

        More generally, if the quadratic relations are of the form
        (T_s-q_1)(T_s-q_2)=0` and `\sqrt{-q_1q_2}` exists then, for a simple
        reflection `s`, the corresponding Kazhdan-Lusztig basis element is:

        .. MATH::

            C_s = (-q_1 q_2)^{1/2} (1 - (-q_1 q_2)^{-1/2} T_s).

        See [KL1979]_ for more details.

        EXAMPLES::

            sage: R.<v> = LaurentPolynomialRing(QQ)
            sage: H = IwahoriHeckeAlgebra('A5', v**2)
            sage: W = H.coxeter_group()
            sage: s1,s2,s3,s4,s5 = W.simple_reflections()
            sage: T = H.T()
            sage: C = H.C()
            sage: T(s1)**2
            -(1-v^2)*T[1] + v^2
            sage: T(C(s1))
            (v^-1)*T[1] - v
            sage: T(C(s1)*C(s2)*C(s1))
            (v^-3)*T[1,2,1] - (v^-1)*T[2,1] - (v^-1)*T[1,2]
             + (v^-1+v)*T[1] + v*T[2] - (v+v^3)

        ::

            sage: R.<v> = LaurentPolynomialRing(QQ)
            sage: H = IwahoriHeckeAlgebra('A3', v**2)
            sage: W = H.coxeter_group()
            sage: s1,s2,s3 = W.simple_reflections()
            sage: C = H.C()
            sage: C(s1*s2*s1)
            C[1,2,1]
            sage: C(s1)**2
            -(v^-1+v)*C[1]
            sage: C(s1)*C(s2)*C(s1)
            C[1,2,1] + C[1]

        TESTS::

            sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
            sage: H = IwahoriHeckeAlgebra('A3', v**2)
            sage: W = H.coxeter_group()
            sage: T = H.T()
            sage: C = H.C()
            sage: Cp = H.Cp()
            sage: all(C(T(C[x])) == C[x] for x in W) # long time
            True
            sage: all(C(Cp(C[x])) == C[x] for x in W) # long time
            True

        Check the defining property between `C` and `C^{\prime}`::

            sage: T(C[1])
            (v^-1)*T[1] - v
            sage: -T(Cp[1]).hash_involution()
            (v^-1)*T[1] - v
            sage: T(Cp[1] + Cp[2]).hash_involution()
            -(v^-1)*T[1] - (v^-1)*T[2] + 2*v
            sage: -T(C[1] + C[2])
            -(v^-1)*T[1] - (v^-1)*T[2] + 2*v
            sage: Cp(-C[1].hash_involution())
            Cp[1]
            sage: Cp(-C[1,2,3].hash_involution())
            Cp[1,2,3]
            sage: Cp(C[1,2,1,3].hash_involution())
            Cp[1,2,3,1]
            sage: all(C((-1)**x.length()*Cp[x].hash_involution()) == C[x] for x in W) # long time
            True
        """
        _basis_name = "C"   # this is used, for example, by specialize_to and is the default prefix

        def hash_involution_on_basis(self, w):
            r"""
            Return the effect of applying the hash involution to the basis
            element ``self[w]``.

            This function is not intended to be called directly. Instead, use
            :meth:`hash_involution`.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
                sage: H = IwahoriHeckeAlgebra('A3', v**2)
                sage: C=H.C()
                sage: s=H.coxeter_group().simple_reflection(1)
                sage: C.hash_involution_on_basis(s)
                -C[1] - (v^-1+v)
                sage: C[s].hash_involution()
                -C[1] - (v^-1+v)
            """
            return (-1)**w.length() * self(self.realization_of().Cp().monomial(w))

    class A(_Basis):
        r"""
        The `A`-basis of an Iwahori-Hecke algebra.

        The `A`-basis of the Iwahori-Hecke algebra is the simplest basis
        that is invariant under the Goldman involution `\#`, up to sign.
        For `w` in the underlying Coxeter group define:

        .. MATH::

            A_w = T_w + (-1)^{\ell(w)}T_w^{\#}
                = T_w + (-1)^{\ell(w)}T_{w^{-1}}^{-1}

        This gives a basis of the Iwahori-Hecke algebra whenever 2 is a unit
        in the base ring. The `A`-basis induces a `\ZZ / 2\ZZ`-grading
        on the Iwahori-Hecke algebra.

        The `A`-basis is a basis only when `2` is invertible. An error
        is raised whenever `2` is not a unit in the base ring.

        EXAMPLES::

            sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
            sage: H = IwahoriHeckeAlgebra('A3', v**2)
            sage: A=H.A(); T=H.T()
            sage: T(A[1])
            T[1] + (1/2-1/2*v^2)
            sage: T(A[1,2])
            T[1,2] + (1/2-1/2*v^2)*T[1] + (1/2-1/2*v^2)*T[2] + (1/2-v^2+1/2*v^4)
            sage: A[1]*A[2]
            A[1,2] - (1/4-1/2*v^2+1/4*v^4)

        TESTS::

            sage: R.<v> = LaurentPolynomialRing(GF(2), 'v')
            sage: H = IwahoriHeckeAlgebra('A3', v**2)
            sage: H.A()
            Traceback (most recent call last):
            ...
            TypeError: the A-basis is defined only when 2 is invertible
        """
        _basis_name = "A"

        def __init__(self, IHAlgebra, prefix=None):
            r"""
            Initialize the `A`-basis of the Iwahori-Hecke algebra ``IHAlgebra``.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ)
                sage: H = IwahoriHeckeAlgebra('A3', v**2)
                sage: A = H.A()
            """
            R = IHAlgebra.base_ring()
            try:
                R(R.one()/2)
            except (TypeError, ZeroDivisionError):
                raise TypeError('the A-basis is defined only when 2 is invertible')

            super().__init__(IHAlgebra, prefix)

            # Define and register coercions from the A basis to the T basis and back again
            from_A_to_T = self.module_morphism(self.to_T_basis, codomain=IHAlgebra.T(),
                                               triangular='lower', key=sorting_key,
                                               category=self.category())
            from_A_to_T.register_as_coercion()
            from_T_to_A = ~from_A_to_T
            from_T_to_A.register_as_coercion()

        def to_T_basis(self, w):
            r"""
            Return the `A`-basis element ``self[w]`` as a linear
            combination of `T`-basis elements.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ)
                sage: H = IwahoriHeckeAlgebra('A3', v**2); A=H.A(); T=H.T()
                sage: s=H.coxeter_group().simple_reflection(1)
                sage: A.to_T_basis(s)
                T[1] + (1/2-1/2*v^2)
                sage: T(A[1,2])
                T[1,2] + (1/2-1/2*v^2)*T[1] + (1/2-1/2*v^2)*T[2] + (1/2-v^2+1/2*v^4)
                sage: A(T[1,2])
                A[1,2] - (1/2-1/2*v^2)*A[1] - (1/2-1/2*v^2)*A[2]
            """
            T = self.realization_of().T()
            return (T.monomial(w) + (-1)**w.length()*T.goldman_involution_on_basis(w)) / 2

        def goldman_involution_on_basis(self, w):
            r"""
            Return the effect of applying the Goldman involution to the basis
            element ``self[w]``.

            This function is not intended to be called directly. Instead, use
            :meth:`goldman_involution`.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
                sage: H = IwahoriHeckeAlgebra('A3', v**2)
                sage: A=H.A()
                sage: s=H.coxeter_group().simple_reflection(1)
                sage: A.goldman_involution_on_basis(s)
                -A[1]
                sage: A[1,2].goldman_involution()
                A[1,2]
            """
            return (-1)**w.length() * self.monomial(w)

    class B(_Basis):
        r"""
        The `B`-basis of an Iwahori-Hecke algebra.

        The `B`-basis is the unique basis of the Iwahori-Hecke algebra that
        is invariant under the Goldman involution, up to sign, and invariant
        under the Kazhdan-Lusztig bar involution. In the generic case, the
        `B`-basis becomes the group basis of the group algebra of the Coxeter
        group the `B`-basis upon setting the Hecke parameters equal to `1`.
        If `w` is an element of the corresponding Coxeter group then
        the `B`-basis element `B_w` is uniquely determined by the conditions
        that `B_w^{\#} = (-1)^{\ell(w)} B_w`, where `\#` is the
        :meth:`Goldman involution <goldman_involution>` and

        .. MATH::

            B_w = T_w + \sum_{v<w}b_{vw}(q) T_v

        where `b_{vw}(q) \neq 0` only if `v < w` in the Bruhat order and
        `\ell(v) \not\equiv \ell(w) \pmod 2`.

        This gives a basis of the Iwahori-Hecke algebra whenever `2` is a
        unit in the base ring. The `B`-basis induces a `\ZZ / 2 \ZZ`-grading
        on the Iwahori-Hecke algebra. The `B`-basis elements are also
        invariant under the Kazhdan-Lusztig bar involution and hence
        are related to the Kazhdan-Lusztig bases.

        The `B`-basis is a basis only when `2` is invertible. An error
        is raised whenever `2` is not a unit in the base ring.

        EXAMPLES::

            sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
            sage: H = IwahoriHeckeAlgebra('A3', v**2)
            sage: A=H.A(); T=H.T(); Cp=H.Cp()
            sage: T(A[1])
            T[1] + (1/2-1/2*v^2)
            sage: T(A[1,2])
            T[1,2] + (1/2-1/2*v^2)*T[1] + (1/2-1/2*v^2)*T[2] + (1/2-v^2+1/2*v^4)
            sage: A[1]*A[2]
            A[1,2] - (1/4-1/2*v^2+1/4*v^4)
            sage: Cp(A[1]*A[2])
            v^2*Cp[1,2] - (1/2*v+1/2*v^3)*Cp[1] - (1/2*v+1/2*v^3)*Cp[2]
             + (1/4+1/2*v^2+1/4*v^4)
            sage: Cp(A[1])
            v*Cp[1] - (1/2+1/2*v^2)
            sage: Cp(A[1,2])
            v^2*Cp[1,2] - (1/2*v+1/2*v^3)*Cp[1]
             - (1/2*v+1/2*v^3)*Cp[2] + (1/2+1/2*v^4)
            sage: Cp(A[1,2,1])
            v^3*Cp[1,2,1] - (1/2*v^2+1/2*v^4)*Cp[2,1]
             - (1/2*v^2+1/2*v^4)*Cp[1,2] + (1/2*v+1/2*v^5)*Cp[1]
             + (1/2*v+1/2*v^5)*Cp[2] - (1/2+1/2*v^6)

        TESTS::

            sage: R.<v> = LaurentPolynomialRing(ZZ, 'v')
            sage: H = IwahoriHeckeAlgebra('A3', v**2)
            sage: H.B()
            Traceback (most recent call last):
            ...
            TypeError: the B-basis is defined only when 2 is invertible
        """
        _basis_name = "B"

        def __init__(self, IHAlgebra, prefix=None):
            r"""
            Initialize the `B`-basis of the Iwahori-Hecke algebra ``IHAlgebra``.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ)
                sage: H = IwahoriHeckeAlgebra('A3', v**2)
                sage: B = H.B()
            """
            R = IHAlgebra.base_ring()
            try:
                R(R.one()/2)
            except (TypeError, ZeroDivisionError):
                raise TypeError('the B-basis is defined only when 2 is invertible')

            super().__init__(IHAlgebra, prefix)

            # Define and register coercions from the B basis to the T basis and back again
            from_B_to_T = self.module_morphism(self.to_T_basis, codomain=IHAlgebra.T(),
                                               triangular='lower', key=sorting_key,
                                               category=self.category())
            from_B_to_T.register_as_coercion()
            from_T_to_B = ~from_B_to_T
            from_T_to_B.register_as_coercion()

        @cached_method
        def to_T_basis(self, w):
            r"""
            Return the `B`-basis element ``self[w]`` as a linear
            combination of `T`-basis elements.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ)
                sage: H = IwahoriHeckeAlgebra('A3', v**2); B=H.B(); T=H.T()
                sage: s=H.coxeter_group().simple_reflection(1)
                sage: B.to_T_basis(s)
                T[1] + (1/2-1/2*v^2)
                sage: T(B[1,2])
                T[1,2] + (1/2-1/2*v^2)*T[1] + (1/2-1/2*v^2)*T[2]
                sage: B(T[1,2])
                B[1,2] - (1/2-1/2*v^2)*B[1] - (1/2-1/2*v^2)*B[2] + (1/2-v^2+1/2*v^4)
            """
            T = self.realization_of().T()
            Bw = T(self.realization_of().A()[w])
            odd = [v for v in Bw.support()
                   if v != w and not (v.length() - w.length()) % 2]
            return Bw - T.sum(Bw.coefficient(v) * self.to_T_basis(v)
                              for v in odd)

        def goldman_involution_on_basis(self, w):
            r"""
            Return the Goldman involution to the basis element
            indexed by ``w``.

            This function is not intended to be called directly. Instead, use
            :meth:`goldman_involution`.

            EXAMPLES::

                sage: R.<v> = LaurentPolynomialRing(QQ, 'v')
                sage: H = IwahoriHeckeAlgebra('A3', v**2)
                sage: B=H.B()
                sage: s=H.coxeter_group().simple_reflection(1)
                sage: B.goldman_involution_on_basis(s)
                -B[1]
                sage: B[1,2].goldman_involution()
                B[1,2]
            """
            return (-1)**w.length() * self.monomial(w)


# The IwahoriHeckeAlgebra_nonstandard class must have the same basis classes as
# the IwahoriHeckeAlgebra class with the same name and they should inherit from
# the respective basis class
class IwahoriHeckeAlgebra_nonstandard(IwahoriHeckeAlgebra):
    r"""
    This is a class which is used behind the scenes by
    :class:`IwahoriHeckeAlgebra` to compute the Kazhdan-Lusztig bases. It is
    not meant to be used directly. It implements the slightly idiosyncratic
    (but convenient) Iwahori-Hecke algebra with two parameters which is
    defined over the Laurent polynomial ring `\ZZ[u,u^{-1},v,v^{-1}]` in
    two variables and has quadratic relations:

    .. MATH::

        (T_r - u)(T_r + v^2/u) = 0.

    The point of these relations is that the product of the two parameters is
    `v^2` which is a square in `\ZZ[u,u^{-1},v,v^{-1}]`. Consequently, the
    Kazhdan-Lusztig bases are defined for this algebra.

    More generally, if we have a Iwahori-Hecke algebra with two parameters
    which has quadratic relations of the form:

    .. MATH::

        (T_r - q_1)(T_r - q_2) = 0

    where `-q_1 q_2` is a square then the Kazhdan-Lusztig bases are
    well-defined for this algebra.  Moreover, these bases be computed by
    specialization from the generic Iwahori-Hecke algebra using the
    specialization which sends `u \mapsto q_1` and `v \mapsto \sqrt{-q_1 q_2}`,
    so that `v^2 / u \mapsto -q_2`.

    For example, if `q_1 = q = Q^2` and `q_2 = -1` then `u \mapsto q` and
    `v \mapsto \sqrt{q} = Q`; this is the standard presentation of the
    Iwahori-Hecke algebra with `(T_r - q)(T_r + 1) = 0`. On the other hand,
    when `q_1 = q` and `q_2 = -q^{-1}` then `u \mapsto q` and `v \mapsto 1`.
    This is the normalized presentation with `(T_r - v)(T_r + v^{-1}) = 0`.

    .. WARNING::

        This class uses non-standard parameters for the Iwahori-Hecke algebra
        and are related to the standard parameters by an outer automorphism
        that is non-trivial on the `T`-basis.
    """
    @staticmethod
    def __classcall_private__(cls, W):
        r"""
        TESTS::

            sage: H1 = sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra_nonstandard("A2")
            sage: W = CoxeterGroup("A2")
            sage: H2 = sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra_nonstandard(W)
            sage: H1 is H2
            True
        """
        if W not in CoxeterGroups():
            W = CoxeterGroup(W)
        return super().__classcall__(cls, W)

    def __init__(self, W):
        r"""
        EXAMPLES::

            sage: H = sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra_nonstandard("A2")
            sage: TestSuite(H).run()
        """
        self._W = W
        self._coxeter_type = W.coxeter_type()

        base_ring = LaurentPolynomialRing(ZZ, 'u,v')
        u, v = base_ring.gens()

        # We don't want to call IwahoriHeckeAlgebra.__init__ because this would
        # try and attach a generic Hecke algebra to this algebra leading to
        # an infinite loop.
        self._q1 = u
        self._q2 = normalized_laurent_polynomial(base_ring, -v**2*u**-1)
        self._root = v

        # Used when multiplying generators: minor speed-up as it avoids the
        # need to constantly add and multiply the parameters when applying the
        # quadratic relation: T^2 = (q1+q2)T - q1*q2
        self._q_sum = normalized_laurent_polynomial(base_ring, self._q1+self._q2)
        self._q_prod = normalized_laurent_polynomial(base_ring, -self._q1*self._q2)

        self.u_inv = normalized_laurent_polynomial(base_ring, u**-1)
        self.v_inv = normalized_laurent_polynomial(base_ring, v**-1)

        self._shorthands = ['C', 'Cp', 'T']

        if W.is_finite():
            self._category = FiniteDimensionalAlgebrasWithBasis(base_ring)
        else:
            self._category = AlgebrasWithBasis(base_ring)
        Parent.__init__(self, base=base_ring, category=self._category.WithRealizations())
        self._is_generic = True  # needed for initialising _KLHeckeBasis

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra_nonstandard("A2")
            A generic Iwahori-Hecke algebra of type A2 in u,-u^-1*v^2 over
             Multivariate Laurent Polynomial Ring in u, v over Integer Ring
        """
        try:
            ct = self._coxeter_type._repr_(compact=True)
        except TypeError:
            ct = repr(self._coxeter_type)
        return "A generic Iwahori-Hecke algebra of type {} in {},{} over {}".format(
                ct, self._q1, self._q2, self.base_ring())

    def _bar_on_coefficients(self, c):
        r"""
        Given a Laurent polynomial ``c`` return the Laurent polynomial obtained
        by applying the (generic) bar involution to ``c`` .

        This is the ring homomorphism of Laurent polynomials in
        `\ZZ[u,u^{-1},v,v^{-1}]` which sends `u` to `u^{-1}` and `v`
        to `v^{-1}`.

        EXAMPLES::

            sage: R.<q>=LaurentPolynomialRing(ZZ)
            sage: H=IwahoriHeckeAlgebra("A3",q^2)
            sage: GH=H._generic_iwahori_hecke_algebra
            sage: GH._bar_on_coefficients(GH.u_inv)
            u
            sage: GH._bar_on_coefficients(GH.v_inv)
            v
        """
        return normalized_laurent_polynomial(self._base, c)(self.u_inv, self.v_inv)

    class _BasesCategory(IwahoriHeckeAlgebra._BasesCategory):
        """
        Category of bases for a generic Iwahori-Hecke algebra.
        """
        def super_categories(self):
            r"""
            The super categories of ``self``.

            EXAMPLES::

                sage: H = sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra_nonstandard("B2")
                sage: H._BasesCategory().super_categories()
                [Category of bases of A generic Iwahori-Hecke algebra of type B2 in u,-u^-1*v^2 over
                  Multivariate Laurent Polynomial Ring in u, v over Integer Ring]
            """
            return [IwahoriHeckeAlgebra._BasesCategory(self.base())]

        class ElementMethods:
            def specialize_to(self, new_hecke):
                r"""
                Return the element in the Iwahori-Hecke algebra ``new_hecke``
                with respect to the same basis which is obtained from ``self``
                by specializing the generic parameters in this algebra to the
                parameters of ``new_hecke``.

                The generic Iwahori-Hecke algebra is defined over
                `\ZZ[u^\pm, v^\pm]` and has parameters ``u`` and
                ``-v^2/u``. The specialization map sends ``u`` to
                ``new_hecke._q1`` and ``v`` to ``new_hecke._root`` which is
                thesquare root of ``-new_hecke._q1*new_hecke._q2``, so
                `-v^2/u` is sent to ``new_hecke._q2``.

                This function is not intended to be called directly. Rather it
                is called behind the scenes to convert between the
                Kazhdan-Lusztig and standard bases of the Iwahori-Hecke
                algebras.

                EXAMPLES::

                    sage: R.<a,b>=LaurentPolynomialRing(ZZ,2)
                    sage: H=IwahoriHeckeAlgebra("A3",a^2,-b^2)
                    sage: GH=H._generic_iwahori_hecke_algebra
                    sage: GH.T()(GH.C()[1])
                    (v^-1)*T[1] + (-u*v^-1)
                    sage: ( GH.T()(GH.C()[1]) ).specialize_to(H)
                    (a^-1*b^-1)*T[1] + (-a*b^-1)
                    sage: GH.C()( GH.T()[1] )
                    v*C[1] + u
                    sage: GH.C()( GH.T()[1] ).specialize_to(H)
                    a*b*C[1] + a^2
                    sage: H.C()( H.T()[1] )
                    a*b*C[1] + a^2
                """
                hecke = self.parent().realization_of()
                q1 = new_hecke._q1
                root = new_hecke._root
                # is there an easier way than this to convert the
                # coefficients to the correct base ring for new_hecke?

                def new_coeff(c):
                    return new_hecke._base(normalized_laurent_polynomial(hecke._base, c)(q1, root))
                new_basis = getattr(new_hecke, self.parent()._basis_name)()
                return new_basis._from_dict({w: new_coeff(c) for w, c in self})

    class T(IwahoriHeckeAlgebra.T):
        r"""
        The `T`-basis for the generic Iwahori-Hecke algebra.
        """
        @cached_method
        def to_Cp_basis(self, w):
            r"""
            Return `T_w` as a linear combination of `C^{\prime}`-basis
            elements.

            EXAMPLES::

                sage: H = sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra_nonstandard("A2")
                sage: s1,s2 = H.coxeter_group().simple_reflections()
                sage: T = H.T()
                sage: Cp = H.Cp()
                sage: T.to_Cp_basis(s1)
                v*Cp[1] + (-u^-1*v^2)
                sage: Cp(T(s1))
                v*Cp[1] + (-u^-1*v^2)
                sage: Cp(T(s1)+1)
                v*Cp[1] + (-u^-1*v^2+1)
                sage: Cp(T(s1*s2)+T(s1)+T(s2)+1)
                v^2*Cp[1,2] + (-u^-1*v^3+v)*Cp[1] + (-u^-1*v^3+v)*Cp[2]
                 + (u^-2*v^4-2*u^-1*v^2+1)
                sage: Cp(T(s1*s2*s1))
                v^3*Cp[1,2,1] + (-u^-1*v^4)*Cp[2,1] + (-u^-1*v^4)*Cp[1,2]
                 + (u^-2*v^5)*Cp[1] + (u^-2*v^5)*Cp[2] + (-u^-3*v^6)
            """
            A = self.realization_of()
            Cp = A.Cp()

            if w == A._W.one():  # the identity element of the Coxeter group
                return Cp.one()

            T0 = self.zero()
            inp = self.monomial(w)
            result = Cp.zero()
            while inp != T0:
                (x, c) = inp.trailing_item(key=sorting_key)
                inp = inp - c * A._root**x.length() * Cp.to_T_basis(x)
                result = result + c * A._root**x.length() * Cp.monomial(x)

            return result

        @cached_method
        def to_C_basis(self, w):
            r"""
            Return `T_w` as a linear combination of `C`-basis elements.

            To compute this we piggy back off the `C^{\prime}`-basis
            conversion using the observation that the hash involution sends
            `T_w` to `(-q_1 q_1)^{\ell(w)} T_w` and `C_w` to
            `(-1)^{\ell(w)} C^{\prime}_w`. Therefore, if

            .. MATH::

                T_w = \sum_v a_{vw} C^{\prime}_v

            then

            .. MATH::

                T_w = (-q_1 q_2)^{\ell(w)} \Big( \sum_v a_{vw} C^{\prime}_v
                                                 \Big)^\#
                    = \sum_v (-1)^{\ell(v)} \overline{a_{vw}} C_v

            Note that we cannot just apply :meth:`hash_involution` here because
            this involution always returns the answer with respect to the
            same basis.

            EXAMPLES::

                sage: H = sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra_nonstandard("A2")
                sage: s1,s2 = H.coxeter_group().simple_reflections()
                sage: T = H.T()
                sage: C = H.C()
                sage: T.to_C_basis(s1)
                v*T[1] + u
                sage: C(T(s1))
                v*C[1] + u
                sage: C(T( C[1] ))
                C[1]
                sage: C(T(s1*s2)+T(s1)+T(s2)+1)
                v^2*C[1,2] + (u*v+v)*C[1] + (u*v+v)*C[2] + (u^2+2*u+1)
                sage: C(T(s1*s2*s1))
                v^3*C[1,2,1] + u*v^2*C[2,1] + u*v^2*C[1,2] + u^2*v*C[1] + u^2*v*C[2] + u^3
            """
            H = self.realization_of()
            q_w = (-H._q_prod)**w.length()
            return self.sum_of_terms((v, (-1)**v.length()*q_w*H._bar_on_coefficients(c))
                                     for (v, c) in self.to_Cp_basis(w))

    class Cp(IwahoriHeckeAlgebra.Cp):
        r"""
        The Kazhdan-Lusztig `C^{\prime}`-basis for the generic Iwahori-Hecke
        algebra.
        """
        @cached_method
        def to_T_basis(self, w):
            r"""
            Return `C^{\prime}_w` as a linear combination of `T`-basis
            elements.

            EXAMPLES::

                sage: H = sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra_nonstandard("A3")
                sage: s1,s2,s3 = H.coxeter_group().simple_reflections()
                sage: T = H.T()
                sage: Cp = H.Cp()
                sage: Cp.to_T_basis(s1)
                (v^-1)*T[1] + (u^-1*v)
                sage: Cp.to_T_basis(s1*s2)
                (v^-2)*T[1,2] + (u^-1)*T[1] + (u^-1)*T[2] + (u^-2*v^2)
                sage: Cp.to_T_basis(s1*s2*s1)
                (v^-3)*T[1,2,1] + (u^-1*v^-1)*T[2,1] + (u^-1*v^-1)*T[1,2]
                 + (u^-2*v)*T[1] + (u^-2*v)*T[2] + (u^-3*v^3)
                sage: T(Cp(s1*s2*s1))
                (v^-3)*T[1,2,1] + (u^-1*v^-1)*T[2,1] + (u^-1*v^-1)*T[1,2]
                 + (u^-2*v)*T[1] + (u^-2*v)*T[2] + (u^-3*v^3)
                sage: T(Cp(s2*s1*s3*s2))
                (v^-4)*T[2,3,1,2] + (u^-1*v^-2)*T[2,3,1] + (u^-1*v^-2)*T[1,2,1]
                 + (u^-1*v^-2)*T[3,1,2] + (u^-1*v^-2)*T[2,3,2] + (u^-2)*T[2,1]
                 + (u^-2)*T[3,1] + (u^-2)*T[1,2] + (u^-2)*T[3,2]
                 + (u^-2)*T[2,3] + (u^-3*v^2)*T[1] + (u^-3*v^2)*T[3]
                 + (u^-1+u^-3*v^2)*T[2] + (u^-2*v^2+u^-4*v^4)
            """
            A = self.realization_of()
            T = A.T()
            Ts = T.algebra_generators()

            if w == A._W.one():  # the identity element of the Coxeter group
                return T.one()

            s = w.first_descent()
            ws = w.apply_simple_reflection(s)

            cpw_s = self.to_T_basis(ws) * A.v_inv * (Ts[s] - A._q2 * T.one())

            i = 1

            def key_func(x):
                return sorting_key(x.leading_support())

            while i < len(cpw_s):
                (x, c) = sorted(cpw_s.terms(), key=key_func)[i].leading_item()
                mu = normalized_laurent_polynomial(A._base, c)[0, -x.length()]    # the coefficient of v^-len(x)
                if mu != 0:
                    cpw_s -= mu * self.to_T_basis(x)
                else:
                    i += 1

            return cpw_s

    C_prime = Cp

    class C(IwahoriHeckeAlgebra.C):
        r"""
        The Kazhdan-Lusztig `C`-basis for the generic Iwahori-Hecke algebra.
        """
        @cached_method
        def to_T_basis(self, w):
            r"""
            Return `C_w` as a linear combination of `T`-basis elements.

            EXAMPLES::

                sage: H = sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra_nonstandard("A3")
                sage: s1,s2,s3 = H.coxeter_group().simple_reflections()
                sage: T = H.T()
                sage: C = H.C()
                sage: C.to_T_basis(s1)
                (v^-1)*T[1] + (-u*v^-1)
                sage: C.to_T_basis(s1*s2)
                (v^-2)*T[1,2] + (-u*v^-2)*T[1] + (-u*v^-2)*T[2] + (u^2*v^-2)
                sage: C.to_T_basis(s1*s2*s1)
                (v^-3)*T[1,2,1] + (-u*v^-3)*T[2,1] + (-u*v^-3)*T[1,2]
                 + (u^2*v^-3)*T[1] + (u^2*v^-3)*T[2] + (-u^3*v^-3)
                sage: T(C(s1*s2*s1))
                (v^-3)*T[1,2,1] + (-u*v^-3)*T[2,1] + (-u*v^-3)*T[1,2]
                 + (u^2*v^-3)*T[1] + (u^2*v^-3)*T[2] + (-u^3*v^-3)
                sage: T(C(s2*s1*s3*s2))
                (v^-4)*T[2,3,1,2] + (-u*v^-4)*T[2,3,1] + (-u*v^-4)*T[1,2,1]
                 + (-u*v^-4)*T[3,1,2] + (-u*v^-4)*T[2,3,2] + (u^2*v^-4)*T[2,1]
                 + (u^2*v^-4)*T[3,1] + (u^2*v^-4)*T[1,2] + (u^2*v^-4)*T[3,2]
                 + (u^2*v^-4)*T[2,3] + (-u^3*v^-4)*T[1] + (-u^3*v^-4)*T[3]
                 + (-u^3*v^-4-u*v^-2)*T[2] + (u^4*v^-4+u^2*v^-2)
            """
            # Treat our index as an index for the C'-basis, convert to the T-basis and
            # then apply the Hecke involution to the result. This gives the
            # desired result because C_w = (-1)^{len(w)) \tau( C_w' ), where
            # \tau is the Hecke involution.
            return (-1)**w.length()*self.realization_of().Cp().to_T_basis(w).hash_involution()


from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.algebras.iwahori_hecke_algebra',
                           'IwahoriHeckeAlgebraT', IwahoriHeckeAlgebra)
