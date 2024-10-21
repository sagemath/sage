# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.graphs sage.groups sage.modules
"""
Weyl Character Rings
"""
# ****************************************************************************
#  Copyright (C) 2011 Daniel Bump <bump at match.stanford.edu>
#                     Nicolas Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sage.combinat.root_system.branching_rules
from sage.categories.algebras import Algebras
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.root_system import RootSystem
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
from sage.misc.functional import is_even
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ


class WeylCharacterRing(CombinatorialFreeModule):
    r"""
    A class for rings of Weyl characters.

    Let `K` be a compact Lie group, which we assume is semisimple and
    simply-connected. Its complexified Lie algebra `L` is the Lie algebra of a
    complex analytic Lie group `G`. The following three categories are
    equivalent: finite-dimensional representations of `K`; finite-dimensional
    representations of `L`; and finite-dimensional analytic representations of
    `G`. In every case, there is a parametrization of the irreducible
    representations by their highest weight vectors. For this theory of Weyl,
    see (for example):

    * Adams, *Lectures on Lie groups*
    * Broecker and Tom Dieck, *Representations of Compact Lie groups*
    * Bump, *Lie Groups*
    * Fulton and Harris, *Representation Theory*
    * Goodman and Wallach, *Representations and Invariants of the Classical Groups*
    * Hall, *Lie Groups, Lie Algebras and Representations*
    * Humphreys, *Introduction to Lie Algebras and their representations*
    * Procesi, *Lie Groups*
    * Samelson, *Notes on Lie Algebras*
    * Varadarajan, *Lie groups, Lie algebras, and their representations*
    * Zhelobenko, *Compact Lie Groups and their Representations*.

    Computations that you can do with these include computing their
    weight multiplicities, products (thus decomposing the tensor
    product of a representation into irreducibles) and branching
    rules (restriction to a smaller group).

    There is associated with `K`, `L` or `G` as above a lattice, the weight
    lattice, whose elements (called weights) are characters of a Cartan
    subgroup or subalgebra. There is an action of the Weyl group `W` on
    the lattice, and elements of a fixed fundamental domain for `W`, the
    positive Weyl chamber, are called dominant. There is for each
    representation a unique highest dominant weight that occurs with
    nonzero multiplicity with respect to a certain partial order, and
    it is called the highest weight vector.

    EXAMPLES::

        sage: L = RootSystem("A2").ambient_space()
        sage: [fw1,fw2] = L.fundamental_weights()
        sage: R = WeylCharacterRing(['A',2], prefix='R')
        sage: [R(1),R(fw1),R(fw2)]
        [R(0,0,0), R(1,0,0), R(1,1,0)]

    Here ``R(1)``, ``R(fw1)``, and ``R(fw2)`` are irreducible representations
    with highest weight vectors `0`, `\Lambda_1`, and `\Lambda_2` respectively
    (the first two fundamental weights).

    For type `A` (also `G_2`, `F_4`, `E_6` and `E_7`) we will take as the
    weight lattice not the weight lattice of the semisimple group, but for a
    larger one. For type `A`, this means we are concerned with the
    representation theory of `K = U(n)` or `G = GL(n, \CC)` rather than `SU(n)`
    or `SU(n, \CC)`. This is useful since the representation theory of `GL(n)`
    is ubiquitous, and also since we may then represent the fundamental
    weights (in :mod:`sage.combinat.root_system.root_system`) by vectors
    with integer entries. If you are only interested in `SL(3)`, say, use
    ``WeylCharacterRing(['A',2])`` as above but be aware that ``R([a,b,c])``
    and ``R([a+1,b+1,c+1])`` represent the same character of `SL(3)` since
    ``R([1,1,1])`` is the determinant.

    For more information, see the thematic tutorial *Lie Methods and
    Related Combinatorics in Sage*, available at:

    https://doc.sagemath.org/html/en/thematic_tutorials/lie.html
    """
    @staticmethod
    def __classcall__(cls, ct, base_ring=ZZ, prefix=None, style='lattice', k=None, conjugate=False, cyclotomic_order=None, fusion_labels=None, inject_variables=False):
        """
        TESTS::

            sage: R = WeylCharacterRing("G2", style='coroots')
            sage: R.cartan_type() is CartanType("G2")
            True
            sage: R.base_ring() is ZZ
            True
        """
        ct = CartanType(ct)
        if prefix is None:
            if ct.is_atomic():
                prefix = ct[0] + str(ct[1])
            else:
                prefix = repr(ct)
        return super().__classcall__(cls, ct, base_ring=base_ring, prefix=prefix, style=style, k=k, conjugate=conjugate, cyclotomic_order=cyclotomic_order, fusion_labels=fusion_labels, inject_variables=inject_variables)

    def __init__(self, ct, base_ring=ZZ, prefix=None, style='lattice', k=None, conjugate=False, cyclotomic_order=None, fusion_labels=None, inject_variables=False):
        """
        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: TestSuite(A2).run()
        """
        ct = CartanType(ct)
        self._cartan_type = ct
        self._rank = ct.rank()
        self._base_ring = base_ring
        self._space = RootSystem(self._cartan_type).ambient_space()
        self._origin = self._space.zero()
        if prefix is None:
            if ct.is_atomic():
                prefix = ct[0] + str(ct[1])
            else:
                prefix = repr(ct)
        self._prefix = prefix
        self._style = style
        self._fusion_labels = None
        self._field = None
        self._basecoer = None
        self._k = k
        if k is not None:
            self._k = Integer(k)
        if ct.is_irreducible():
            self._opposition = ct.opposition_automorphism()
            self._highest = self._space.highest_root()
            self._hip = self._highest.inner_product(self._highest)
        if style == "coroots":
            self._word = self._space.weyl_group().long_element().reduced_word()

        # Set the basis
        if k is not None:
            self._prefix += str(k)
            fw = self._space.fundamental_weights()

            def next_level(wt):
                return [wt + la for la in fw if self.level(wt + la) <= k]
            B = list(RecursivelyEnumeratedSet([self._space.zero()], next_level))
            B = [self._space.from_vector_notation(wt, style='coroots') for wt in B]
        else:
            B = self._space

        cat = AlgebrasWithBasis(base_ring).Commutative()
        if k is None:
            cat = cat.Subobjects().Graded()
        else:
            cat = cat.FiniteDimensional()
        CombinatorialFreeModule.__init__(self, base_ring, B, category=cat)

        # Register the embedding of self into ambient as a coercion
        self.lift.register_as_coercion()
        # Register the partial inverse as a conversion
        self.register_conversion(self.retract)

        # Record properties of the FusionRing
        # mg = square of long to short root lengths
        # nf = normalizing factor for the inner product
        # fg = order of the fundamental group (except for Type B)
        if k is not None:
            if ct[0] in ['A', 'D', 'E']:
                self._m_g = 1
            elif ct[0] in ['B', 'C', 'F']:
                self._m_g = 2
            else:
                self._m_g = 3
            if ct[0] in ['B', 'F']:
                self._nf = 2
            else:
                self._nf = 1
            self._h_check = ct.dual_coxeter_number()
            self._l = self._m_g * (self._k + self._h_check)
            if conjugate:
                self._conj = -1
            else:
                self._conj = 1
            if ct[0] == 'A':
                self._fg = ct[1] + 1
            elif ct[0] == 'E' and ct[1] == 6:
                self._fg = 3
            elif ct[0] == 'E' and ct[1] == 7:
                self._fg = 2
            elif ct[0] == 'D':
                self._fg = 2
            else:
                self._fg = 1
            if cyclotomic_order is None:
                self._cyclotomic_order = self._fg * self._l
            else:
                self._cyclotomic_order = cyclotomic_order
            self._fusion_labels = fusion_labels
            if fusion_labels:
                self.fusion_labels(labels=fusion_labels, inject_variables=inject_variables)

    @cached_method
    def ambient(self):
        """
        Return the weight ring of ``self``.

        EXAMPLES::

            sage: WeylCharacterRing("A2").ambient()
            The Weight ring attached to The Weyl Character Ring of Type A2 with Integer Ring coefficients
        """
        return WeightRing(self)

    # Eventually, one could want to put the cache_method here rather
    # than on _irr_weights. Or just to merge this method and _irr_weights
    def lift_on_basis(self, irr):
        """
        Expand the basis element indexed by the weight ``irr`` into the
        weight ring of ``self``.

        INPUT:

        - ``irr`` -- a dominant weight

        This is used to implement :meth:`lift`.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: v = A2._space([2,1,0]); v
            (2, 1, 0)
            sage: A2.lift_on_basis(v)
            2*a2(1,1,1) + a2(1,2,0) + a2(1,0,2) + a2(2,1,0) + a2(2,0,1) + a2(0,1,2) + a2(0,2,1)

        This is consistent with the analogous calculation with symmetric
        Schur functions::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s[2,1].expand(3)
            x0^2*x1 + x0*x1^2 + x0^2*x2 + 2*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2
        """
        return self.ambient()._from_dict(self._irr_weights(irr))

    def demazure_character(self, hwv, word, debug=False):
        r"""
        Compute the Demazure character.

        INPUT:

        - ``hwv`` -- a (usually dominant) weight
        - ``word`` -- a Weyl group word

        Produces the Demazure character with highest weight ``hwv`` and
        ``word`` as an element of the weight ring. Only available if
        ``style="coroots"``. The Demazure operators are also available as
        methods of :class:`WeightRing` elements, and as methods of crystals.
        Given a
        :class:`~sage.combinat.crystals.tensor_product.CrystalOfTableaux`
        with given highest weight vector, the Demazure method on the
        crystal will give the equivalent of this method, except that
        the Demazure character of the crystal is given as a sum of
        monomials instead of an element of the :class:`WeightRing`.

        See :meth:`WeightRing.Element.demazure` and
        :meth:`sage.categories.classical_crystals.ClassicalCrystals.ParentMethods.demazure_character`

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2",style='coroots')
            sage: h = sum(A2.fundamental_weights()); h
            (2, 1, 0)
            sage: A2.demazure_character(h,word=[1,2])
            a2(0,0) + a2(-2,1) + a2(2,-1) + a2(1,1) + a2(-1,2)
            sage: A2.demazure_character((1,1),word=[1,2])
            a2(0,0) + a2(-2,1) + a2(2,-1) + a2(1,1) + a2(-1,2)
        """
        if self._style != "coroots":
            raise ValueError('demazure method unavailable: use style="coroots"')
        hwv = self._space.from_vector_notation(hwv, style='coroots')
        return self.ambient()._from_dict(self._demazure_weights(hwv, word=word, debug=debug))

    @lazy_attribute
    def lift(self):
        """
        The embedding of ``self`` into its weight ring.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: A2.lift
            Generic morphism:
              From: The Weyl Character Ring of Type A2 with Integer Ring coefficients
              To:   The Weight ring attached to The Weyl Character Ring of Type A2 with Integer Ring coefficients

        ::

            sage: x = -A2(2,1,1) - A2(2,2,0) + A2(3,1,0)
            sage: A2.lift(x)
            a2(1,3,0) + a2(1,0,3) + a2(3,1,0) + a2(3,0,1) + a2(0,1,3) + a2(0,3,1)

        As a shortcut, you may also do::

            sage: x.lift()
            a2(1,3,0) + a2(1,0,3) + a2(3,1,0) + a2(3,0,1) + a2(0,1,3) + a2(0,3,1)

        Or even::

            sage: a2 = WeightRing(A2)
            sage: a2(x)
            a2(1,3,0) + a2(1,0,3) + a2(3,1,0) + a2(3,0,1) + a2(0,1,3) + a2(0,3,1)
        """
        return self.module_morphism(self.lift_on_basis,
                                    codomain=self.ambient(),
                                    category=AlgebrasWithBasis(self.base_ring()))

    def _retract(self, chi):
        """
        Construct a Weyl character from an invariant element of the weight ring.

        INPUT:

        - ``chi`` -- a linear combination of weights which
          shall be invariant under the action of the Weyl group

        OUTPUT: the corresponding Weyl character

        Please use instead the morphism :meth:`retract` which is
        implemented using this method.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: a2 = WeightRing(A2)

        ::

            sage: v = A2._space([3,1,0]); v
            (3, 1, 0)
            sage: chi = a2.sum_of_monomials(v.orbit()); chi
            a2(1,3,0) + a2(1,0,3) + a2(3,1,0) + a2(3,0,1) + a2(0,1,3) + a2(0,3,1)
            sage: A2._retract(chi)
            -A2(2,1,1) - A2(2,2,0) + A2(3,1,0)
        """
        return self.char_from_weights(dict(chi))

    @lazy_attribute
    def retract(self):
        """
        The partial inverse map from the weight ring into ``self``.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: a2 = WeightRing(A2)
            sage: A2.retract
            Generic morphism:
              From: The Weight ring attached to The Weyl Character Ring of Type A2 with Integer Ring coefficients
              To:   The Weyl Character Ring of Type A2 with Integer Ring coefficients

        ::

            sage: v = A2._space([3,1,0]); v
            (3, 1, 0)
            sage: chi = a2.sum_of_monomials(v.orbit()); chi
            a2(1,3,0) + a2(1,0,3) + a2(3,1,0) + a2(3,0,1) + a2(0,1,3) + a2(0,3,1)
            sage: A2.retract(chi)
            -A2(2,1,1) - A2(2,2,0) + A2(3,1,0)

        The input should be invariant::

            sage: A2.retract(a2.monomial(v))
            Traceback (most recent call last):
            ...
            ValueError: multiplicity dictionary may not be Weyl group invariant

        As a shortcut, you may use conversion::

            sage: A2(chi)
            -A2(2,1,1) - A2(2,2,0) + A2(3,1,0)
            sage: A2(a2.monomial(v))
            Traceback (most recent call last):
            ...
            ValueError: multiplicity dictionary may not be Weyl group invariant
        """
        from sage.categories.homset import Hom
        from sage.categories.morphism import SetMorphism
        category = Algebras(self.base_ring())
        return SetMorphism(Hom(self.ambient(), self, category), self._retract)

    def _repr_(self):
        """
        EXAMPLES::

            sage: WeylCharacterRing("A3")
            The Weyl Character Ring of Type A3 with Integer Ring coefficients
        """
        if self._k is None:
            return "The Weyl Character Ring of Type {} with {} coefficients".format(self._cartan_type._repr_(compact=True), self._base_ring)
        else:
            return "The Fusion Ring of Type {} and level {} with {} coefficients".format(self._cartan_type._repr_(compact=True), self._k, self._base_ring)

    def __call__(self, *args):
        """
        Construct an element of ``self``.

        The input can either be an object that can be coerced or
        converted into ``self`` (an element of ``self``, of the base
        ring, of the weight ring), or a dominant weight. In the later
        case, the basis element indexed by that weight is returned.

        To specify the weight, you may give it explicitly. Alternatively,
        you may give a tuple of integers. Normally these are the
        components of the vector in the standard realization of
        the weight lattice as a vector space. Alternatively, if
        the ring is constructed with ``style = "coroots"``, you may
        specify the weight by giving a set of integers, one for each
        fundamental weight; the weight is then the linear combination
        of the fundamental weights with these coefficients.

        As a syntactical shorthand, for tuples of length at least two,
        the parenthesis may be omitted.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: [A2(x) for x in [-2,-1,0,1,2]]
            [-2*A2(0,0,0), -A2(0,0,0), 0, A2(0,0,0), 2*A2(0,0,0)]
            sage: [A2(2,1,0), A2([2,1,0]), A2(2,1,0)== A2([2,1,0])]
            [A2(2,1,0), A2(2,1,0), True]
            sage: A2([2,1,0]) == A2(2,1,0)
            True
            sage: l = -2*A2(0,0,0) - A2(1,0,0) + A2(2,0,0) + 2*A2(3,0,0)
            sage: [l in A2, A2(l) == l]
            [True, True]
            sage: P.<q> = QQ[]
            sage: A2 = WeylCharacterRing(['A',2], base_ring = P)
            sage: [A2(x) for x in [-2,-1,0,1,2,-2*q,-q,q,2*q,(1-q)]]
            [-2*A2(0,0,0), -A2(0,0,0), 0, A2(0,0,0), 2*A2(0,0,0), -2*q*A2(0,0,0), -q*A2(0,0,0),
            q*A2(0,0,0), 2*q*A2(0,0,0), (-q+1)*A2(0,0,0)]
            sage: R.<q> = ZZ[]
            sage: A2 = WeylCharacterRing(['A',2], base_ring = R, style='coroots')
            sage: q*A2(1)
            q*A2(0,0)
            sage: [A2(x) for x in [-2,-1,0,1,2,-2*q,-q,q,2*q,(1-q)]]
            [-2*A2(0,0), -A2(0,0), 0, A2(0,0), 2*A2(0,0), -2*q*A2(0,0), -q*A2(0,0), q*A2(0,0), 2*q*A2(0,0), (-q+1)*A2(0,0)]
        """
        # The purpose of this __call__ method is only to handle the
        # syntactical shorthand; otherwise it just delegates the work
        # to the coercion model, which itself will call
        # _element_constructor_ if the input is made of exactly one
        # object which can't be coerced into self
        if len(args) > 1:
            args = (args,)
        return super().__call__(*args)

    def _element_constructor_(self, weight):
        """
        Construct a monomial from a dominant weight.

        INPUT:

        - ``weight`` -- an element of the weight space, or a tuple

        This method is responsible for constructing an appropriate
        dominant weight from ``weight``, and then return the monomial
        indexed by that weight. See :meth:`__call__` and
        :meth:`sage.combinat.root_system.ambient_space.AmbientSpace.from_vector`.

        TESTS::

            sage: A2 = WeylCharacterRing("A2")
            sage: A2._element_constructor_([2,1,0])
            A2(2,1,0)
        """
        weight = self._space.from_vector_notation(weight, style=self._style)
        if not weight.is_dominant_weight():
            raise ValueError("{} is not a dominant element of the weight lattice".format(weight))
        if self._k is not None:
            if self.level(weight) > self._k:
                raise ValueError("{} has level greater than {}".format(weight, self._k))
        return self.monomial(weight)

    def product_on_basis(self, a, b):
        r"""
        Compute the tensor product of two irreducible representations ``a``
        and ``b``.

        EXAMPLES::

            sage: D4 = WeylCharacterRing(['D',4])
            sage: spin_plus = D4(1/2,1/2,1/2,1/2)
            sage: spin_minus = D4(1/2,1/2,1/2,-1/2)
            sage: spin_plus * spin_minus # indirect doctest
            D4(1,0,0,0) + D4(1,1,1,0)
            sage: spin_minus * spin_plus
            D4(1,0,0,0) + D4(1,1,1,0)

        Uses the Brauer-Klimyk method.
        """
        # The method is asymmetrical, and as a rule of thumb
        # it is fastest to switch the factors so that the
        # smaller character is the one that is decomposed
        # into weights.
        if sum(a.coefficients()) > sum(b.coefficients()):
            a, b = b, a
        return self._product_helper(self._irr_weights(a), b)

    def _product_helper(self, d1, b):
        """
        Helper function for :meth:`product_on_basis`.

        INPUT:

        - ``d1`` -- dictionary of weight multiplicities
        - ``b`` -- a dominant weight

        If ``d1`` is the dictionary of weight multiplicities of a character,
        returns the product of that character by the irreducible character
        with highest weight ``b``.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: r = A2(1,0,0)
            sage: [A2._product_helper(r.weight_multiplicities(),x) for x in A2.space().fundamental_weights()]
            [A2(1,1,0) + A2(2,0,0), A2(1,1,1) + A2(2,1,0)]
        """
        d = {}
        for k in d1:
            [epsilon, g] = self.dot_reduce(b + k)
            if epsilon == 1:
                d[g] = d.get(g, 0) + d1[k]
            elif epsilon == -1:
                d[g] = d.get(g, 0) - d1[k]
        return self._from_dict(d, coerce=True)

    def dot_reduce(self, a):
        r"""
        Auxiliary function for :meth:`product_on_basis`.

        Return a pair `[\epsilon, b]` where `b` is a dominant weight and
        `\epsilon` is 0, 1 or -1. To describe `b`, let `w` be an element of
        the Weyl group such that `w(a + \rho)` is dominant. If
        `w(a + \rho) - \rho` is dominant, then `\epsilon` is the sign of
        `w` and `b` is `w(a + \rho) - \rho`. Otherwise, `\epsilon` is zero.

        INPUT:

        - ``a`` -- a weight

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: weights = sorted(A2(2,1,0).weight_multiplicities().keys(), key=str); weights
            [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 1, 1), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
            sage: [A2.dot_reduce(x) for x in weights]
            [[0, (0, 0, 0)], [-1, (1, 1, 1)], [-1, (1, 1, 1)], [1, (1, 1, 1)], [0, (0, 0, 0)], [0, (0, 0, 0)], [1, (2, 1, 0)]]
        """
        alphacheck = self._space.simple_coroots()
        alpha = self._space.simple_roots()
        [epsilon, ret] = [1, a]
        done = False
        while not done:
            done = True
            for i in self._space.index_set():
                c = ret.inner_product(alphacheck[i])
                if c == -1:
                    return [0, self._space.zero()]
                elif c < -1:
                    epsilon = -epsilon
                    ret -= (1 + c) * alpha[i]
                    done = False
                    break
            if self._k is not None:
                l = self.level(ret)
                k = self._k
                if l > k:
                    if l == k + 1:
                        return [0, self._space.zero()]
                    else:
                        epsilon = -epsilon
                        ret = self.affine_reflect(ret, k + 1)
                        done = False
        return [epsilon, ret]

    def affine_reflect(self, wt, k=0):
        r"""
        Return the reflection of wt in the hyperplane `\theta`.

        Optionally, this also shifts by a multiple `k` of `\theta`.

        INPUT:

        - ``wt`` -- a weight
        - ``k`` -- (optional) a positive integer

        EXAMPLES::

            sage: B22 = FusionRing("B2",2)
            sage: fw = B22.fundamental_weights(); fw
            Finite family {1: (1, 0), 2: (1/2, 1/2)}
            sage: [B22.affine_reflect(x,2) for x in fw]
            [(2, 1), (3/2, 3/2)]
        """
        coef = ZZ(2 * wt.inner_product(self._highest) / self._hip)
        return wt + (k - coef) * self._highest

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: WeylCharacterRing("A3").some_elements()
            [A3(1,0,0,0), A3(1,1,0,0), A3(1,1,1,0)]
        """
        return [self.monomial(x) for x in self.fundamental_weights()]

    def one_basis(self):
        """
        Return the index of 1 in ``self``.

        EXAMPLES::

            sage: WeylCharacterRing("A3").one_basis()
            (0, 0, 0, 0)
            sage: WeylCharacterRing("A3").one()
            A3(0,0,0,0)
        """
        return self._space.zero()

    @cached_method
    def _irr_weights(self, hwv):
        """
        Compute the weights of an irreducible as a dictionary.

        Given a dominant weight ``hwv``, this produces a dictionary of
        weight multiplicities for the irreducible representation
        with highest weight vector ``hwv``. This method is cached
        for efficiency.

        INPUT:

        - ``hwv`` -- a dominant weight

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: v = A2.fundamental_weights()[1]; v
            (1, 0, 0)
            sage: A2._irr_weights(v)
            {(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1}
        """
        if self._style == "coroots":
            return self._demazure_weights(hwv)
        else:
            return irreducible_character_freudenthal(hwv)

    def _demazure_weights(self, hwv, word='long', debug=False):
        """
        Compute the weights of a Demazure character.

        This method duplicates the functionality of :meth:`_irr_weights`, under
        the assumption that ``style = "coroots"``, but allows an optional
        parameter ``word``. (This is not allowed in :meth:`_irr_weights` since
        it would interfere with the ``@cached_method``.) Produces the
        dictionary of weights for the irreducible character with highest
        weight ``hwv`` when ``word`` is omitted, or for the Demazure character
        if ``word`` is included.

        INPUT:

        - ``hwv`` -- a dominant weight

        EXAMPLES::

            sage: B2 = WeylCharacterRing("B2", style='coroots')
            sage: [B2._demazure_weights(v, word=[1,2]) for v in B2.fundamental_weights()]
            [{(1, 0): 1, (0, 1): 1}, {(-1/2, 1/2): 1, (1/2, -1/2): 1, (1/2, 1/2): 1}]
        """
        alphacheck = self._space.simple_coroots()
        dd = {}
        h = tuple(int(hwv.inner_product(alphacheck[j]))
                  for j in self._space.index_set())
        dd[h] = 1
        return self._demazure_helper(dd, word=word, debug=debug)

    def _demazure_helper(self, dd, word='long', debug=False):
        r"""
        Assumes ``style = "coroots"``. If the optional parameter ``word`` is
        specified, produces a Demazure character (defaults to the long Weyl
        group element.

        INPUT:

        - ``dd`` -- dictionary of weights

        - ``word`` -- (optional) a Weyl group reduced word

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2",style='coroots')
            sage: dd = {}; dd[(1,1)]=int(1)
            sage: A2._demazure_helper(dd,word=[1,2])
            {(0, 0, 0): 1, (-1, 1, 0): 1, (1, -1, 0): 1, (1, 0, -1): 1, (0, 1, -1): 1}
        """
        if self._style != "coroots":
            raise ValueError('_demazure_helper method unavailable: use style="coroots"')
        index_set = self._space.index_set()
        alphacheck = self._space.simple_coroots()
        alpha = self._space.simple_roots()
        r = self.rank()
        cm = {}
        supp = []
        for i in index_set:
            temp = []
            cm[i] = [0] * r
            for ind, j in enumerate(index_set):
                cm[i][ind] = int(alpha[i].inner_product(alphacheck[j]))
                if cm[i][ind]:
                    temp.append(ind)
            supp.append(temp)
            if debug:
                print("cm[%s]=%s" % (i, cm[i]))
        accum = dd
        if word == "long":
            word = self._word
        for i in reversed(word):
            if debug:
                print("i=%s" % i)
            next = {}
            for v in accum:
                coroot = v[i - 1]
                if debug:
                    print("   v=%s, coroot=%s" % (v, coroot))
                if coroot >= 0:
                    mu = v
                    for j in range(coroot + 1):
                        next[mu] = next.get(mu, 0) + accum[v]
                        if debug:
                            print("     mu=%s, next[mu]=%s" % (mu, next[mu]))
                        mu = list(mu)
                        for k in supp[i - 1]:
                            mu[k] -= cm[i][k]
                        mu = tuple(mu)
                else:
                    mu = v
                    for j in range(-1 - coroot):
                        mu = list(mu)
                        for k in supp[i - 1]:
                            mu[k] += cm[i][k]
                        mu = tuple(mu)
                        next[mu] = next.get(mu, 0) - accum[v]
                        if debug:
                            print("     mu=%s, next[mu]=%s" % (mu, next[mu]))
            accum = {}
            for v in next:
                accum[v] = next[v]
        ret = {}
        for v in accum:
            if accum[v]:
                ret[self._space.from_vector_notation(v, style='coroots')] = accum[v]
        return ret

    @cached_method
    def _weight_multiplicities(self, x):
        """
        Produce weight multiplicities for the (possibly reducible)
        WeylCharacter ``x``.

        EXAMPLES::

            sage: B2 = WeylCharacterRing("B2",style='coroots')
            sage: chi = 2*B2(1,0)
            sage: B2._weight_multiplicities(chi)
            {(0, 0): 2, (-1, 0): 2, (1, 0): 2, (0, -1): 2, (0, 1): 2}
        """
        d = {}
        m = x._monomial_coefficients
        for k in m:
            c = m[k]
            d1 = self._irr_weights(k)
            for l in d1:
                if l in d:
                    d[l] += c * d1[l]
                else:
                    d[l] = c * d1[l]
        for k in list(d):
            if d[k] == 0:
                del d[k]
            else:
                d[k] = self._base_ring(d[k])
        return d

    def base_ring(self):
        """
        Return the base ring of ``self``.

        EXAMPLES::

            sage: R = WeylCharacterRing(['A',3], base_ring = CC); R.base_ring()
            Complex Field with 53 bits of precision
        """
        return self._base_ring

    def irr_repr(self, hwv):
        """
        Return a string representing the irreducible character with highest
        weight vector ``hwv``.

        EXAMPLES::

            sage: B3 = WeylCharacterRing("B3")
            sage: [B3.irr_repr(v) for v in B3.fundamental_weights()]
            ['B3(1,0,0)', 'B3(1,1,0)', 'B3(1/2,1/2,1/2)']
            sage: B3 = WeylCharacterRing("B3", style='coroots')
            sage: [B3.irr_repr(v) for v in B3.fundamental_weights()]
            ['B3(1,0,0)', 'B3(0,1,0)', 'B3(0,0,1)']
        """
        return self._prefix + self._wt_repr(hwv)

    def level(self, wt):
        """
        Return the level of the weight, defined to be the value of
        the weight on the coroot associated with the highest root.

        EXAMPLES::

            sage: R = FusionRing("F4",2); [R.level(x) for x in R.fundamental_weights()]
            [2, 3, 2, 1]
            sage: [CartanType("F4~").dual().a()[x] for x in [1..4]]
            [2, 3, 2, 1]
        """
        return ZZ(2 * wt.inner_product(self._highest) / self._hip)

    def _dual_helper(self, wt):
        """
        If `w_0` is the long Weyl group element and `wt` is an
        element of the weight lattice, this returns `-w_0(wt)`.

        EXAMPLES::

            sage: A3=WeylCharacterRing("A3")
            sage: [A3._dual_helper(x) for x in A3.fundamental_weights()]
            [(0, 0, 0, -1), (0, 0, -1, -1), (0, -1, -1, -1)]
        """
        if self.cartan_type()[0] == 'A':  # handled separately for GL(n) compatibility
            return self.space()([-x for x in reversed(wt.to_vector().list())])
        ret = 0
        alphacheck = self._space.simple_coroots()
        fw = self._space.fundamental_weights()
        for i in self._space.index_set():
            ret += wt.inner_product(alphacheck[i]) * fw[self._opposition[i]]
        return ret

    def _wt_repr(self, wt):
        """
        Produce a representation of a vector in either coweight or
        lattice notation (following the appendices in Bourbaki, Lie Groups and
        Lie Algebras, Chapters 4,5,6), depending on whether the parent
        :class:`WeylCharacterRing` is created with ``style="coweights"``
        or not.

        EXAMPLES::

            sage: [fw1,fw2]=RootSystem("G2").ambient_space().fundamental_weights(); fw1,fw2
            ((1, 0, -1), (2, -1, -1))
            sage: [WeylCharacterRing("G2")._wt_repr(v) for v in [fw1,fw2]]
            ['(1,0,-1)', '(2,-1,-1)']
            sage: [WeylCharacterRing("G2",style='coroots')._wt_repr(v) for v in [fw1,fw2]]
            ['(1,0)', '(0,1)']
        """
        if self._style == "lattice":
            vec = wt.to_vector()
        elif self._style == "coroots":
            vec = [wt.inner_product(x) for x in self.simple_coroots()]
        else:
            raise ValueError("unknown style")
        hstring = str(vec[0])
        for i in range(1, len(vec)):
            hstring = hstring + "," + str(vec[i])
        return "(" + hstring + ")"

    def _repr_term(self, t):
        """
        Representation of the monomial corresponding to a weight ``t``.

        EXAMPLES::

            sage: G2 = WeylCharacterRing("G2") # indirect doctest
            sage: [G2._repr_term(x) for x in G2.fundamental_weights()]
            ['G2(1,0,-1)', 'G2(2,-1,-1)']
        """
        if self._fusion_labels is not None:
            t = tuple([t.inner_product(x) for x in self.simple_coroots()])
            return self._fusion_labels[t]
        else:
            return self.irr_repr(t)

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: WeylCharacterRing("A2").cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def fundamental_weights(self):
        """
        Return the fundamental weights.

        EXAMPLES::

            sage: WeylCharacterRing("G2").fundamental_weights()
            Finite family {1: (1, 0, -1), 2: (2, -1, -1)}
        """
        return self._space.fundamental_weights()

    def simple_roots(self):
        """
        Return the simple roots.

        EXAMPLES::

            sage: WeylCharacterRing("G2").simple_roots()
            Finite family {1: (0, 1, -1), 2: (1, -2, 1)}
        """
        return self._space.simple_roots()

    def simple_coroots(self):
        """
        Return the simple coroots.

        EXAMPLES::

            sage: WeylCharacterRing("G2").simple_coroots()
            Finite family {1: (0, 1, -1), 2: (1/3, -2/3, 1/3)}
        """
        return self._space.simple_coroots()

    def highest_root(self):
        """
        Return the highest root.

        EXAMPLES::

            sage: WeylCharacterRing("G2").highest_root()
             (2, -1, -1)
        """
        return self._space.highest_root()

    def positive_roots(self):
        """
        Return the positive roots.

        EXAMPLES::

            sage: WeylCharacterRing("G2").positive_roots()
            [(0, 1, -1), (1, -2, 1), (1, -1, 0), (1, 0, -1), (1, 1, -2), (2, -1, -1)]
        """
        return self._space.positive_roots()

    def dynkin_diagram(self):
        """
        Return the Dynkin diagram of ``self``.

        EXAMPLES::

            sage: WeylCharacterRing("E7").dynkin_diagram()
                    O 2
                    |
                    |
            O---O---O---O---O---O
            1   3   4   5   6   7
            E7
        """
        return self.space().dynkin_diagram()

    def extended_dynkin_diagram(self):
        """
        Return the extended Dynkin diagram, which is the Dynkin diagram
        of the corresponding untwisted affine type.

        EXAMPLES::

            sage: WeylCharacterRing("E7").extended_dynkin_diagram()
                        O 2
                        |
                        |
            O---O---O---O---O---O---O
            0   1   3   4   5   6   7
            E7~
        """
        return self.cartan_type().affine().dynkin_diagram()

    def rank(self):
        """
        Return the rank.

        EXAMPLES::

            sage: WeylCharacterRing("G2").rank()
            2
        """
        return self._rank

    def space(self):
        """
        Return the weight space associated to ``self``.

        EXAMPLES::

            sage: WeylCharacterRing(['E',8]).space()
            Ambient space of the Root system of type ['E', 8]
        """
        return self._space

    def char_from_weights(self, mdict):
        """
        Construct a Weyl character from an invariant linear combination
        of weights.

        INPUT:

        - ``mdict`` -- dictionary mapping weights to coefficients,
          and representing a linear combination of weights which
          shall be invariant under the action of the Weyl group

        OUTPUT: the corresponding Weyl character

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: v = A2._space([3,1,0]); v
            (3, 1, 0)
            sage: d = dict([(x,1) for x in v.orbit()]); d
            {(1, 3, 0): 1,
             (1, 0, 3): 1,
             (3, 1, 0): 1,
             (3, 0, 1): 1,
             (0, 1, 3): 1,
             (0, 3, 1): 1}
            sage: A2.char_from_weights(d)
            -A2(2,1,1) - A2(2,2,0) + A2(3,1,0)
        """
        return self._from_dict(self._char_from_weights(mdict), coerce=True)

    def _char_from_weights(self, mdict):
        """
        Helper method for :meth:`char_from_weights`.

        INPUT:

        - ``mdict`` -- dictionary of weight multiplicities

        The output of this method is a dictionary whose keys are dominant
        weights that is the same as the :meth:`monomial_coefficients` method
        of ``self.char_from_weights()``.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: v = A2._space([3,1,0])
            sage: d = dict([(x,1) for x in v.orbit()])
            sage: A2._char_from_weights(d)
            {(2, 1, 1): -1, (2, 2, 0): -1, (3, 1, 0): 1}
        """
        hdict = {}
        ddict = mdict.copy()
        while ddict:
            highest = max((x.inner_product(self._space.rho()), x) for x in ddict)[1]
            if not highest.is_dominant():
                raise ValueError("multiplicity dictionary may not be Weyl group invariant")
            sdict = self._irr_weights(highest)
            c = ddict[highest]
            if highest in hdict:
                hdict[highest] += c
            else:
                hdict[highest] = c
            for k in sdict:
                if k in ddict:
                    if ddict[k] == c*sdict[k]:
                        del ddict[k]
                    else:
                        ddict[k] = ddict[k]-c*sdict[k]
                else:
                    ddict[k] = -c*sdict[k]
        return hdict

    def adjoint_representation(self):
        """
        Return the adjoint representation as an element of the WeylCharacterRing.

        EXAMPLES::

            sage: G2 = WeylCharacterRing("G2",style='coroots')
            sage: G2.adjoint_representation()
            G2(0,1)
        """
        return self(self.highest_root())

    def maximal_subgroups(self):
        r"""
        This method is only available if the Cartan type of
        ``self`` is irreducible and of rank no greater than 8.
        This method produces a list of the maximal subgroups
        of ``self``, up to (possibly outer) automorphisms. Each line
        in the output gives the Cartan type of a maximal subgroup
        followed by a command that creates the branching rule.

        EXAMPLES::

            sage: WeylCharacterRing("E6").maximal_subgroups()
            D5:branching_rule("E6","D5","levi")
            C4:branching_rule("E6","C4","symmetric")
            F4:branching_rule("E6","F4","symmetric")
            A2:branching_rule("E6","A2","miscellaneous")
            G2:branching_rule("E6","G2","miscellaneous")
            A2xG2:branching_rule("E6","A2xG2","miscellaneous")
            A1xA5:branching_rule("E6","A1xA5","extended")
            A2xA2xA2:branching_rule("E6","A2xA2xA2","extended")

        Note that there are other embeddings of (for example
        `A_2` into `E_6` as nonmaximal subgroups. These
        embeddings may be constructed by composing branching
        rules through various subgroups.

        Once you know which maximal subgroup you are interested
        in, to create the branching rule, you may either
        paste the command to the right of the colon from the
        above output onto the command line, or alternatively
        invoke the related method :meth:`maximal_subgroup`::

            sage: branching_rule("E6","G2","miscellaneous")
            miscellaneous branching rule E6 => G2
            sage: WeylCharacterRing("E6").maximal_subgroup("G2")
            miscellaneous branching rule E6 => G2

        It is believed that the list of maximal subgroups is complete, except that some
        subgroups may be not be invariant under outer automorphisms. It is reasonable
        to want a list of maximal subgroups that is complete up to conjugation,
        but to obtain such a list you may have to apply outer automorphisms.
        The group of outer automorphisms modulo inner automorphisms is isomorphic
        to the group of symmetries of the Dynkin diagram, and these are available
        as branching rules. The following example shows that while
        a branching rule from `D_4` to `A_1\times C_2` is supplied,
        another different one may be obtained by composing it with the
        triality automorphism of `D_4`::

            sage: [D4,A1xC2]=[WeylCharacterRing(x,style='coroots') for x in ["D4","A1xC2"]]
            sage: fw = D4.fundamental_weights()
            sage: b = D4.maximal_subgroup("A1xC2")
            sage: [D4(fw).branch(A1xC2,rule=b) for fw in D4.fundamental_weights()]
            [A1xC2(1,1,0),
            A1xC2(2,0,0) + A1xC2(2,0,1) + A1xC2(0,2,0),
            A1xC2(1,1,0),
            A1xC2(2,0,0) + A1xC2(0,0,1)]
            sage: b1 = branching_rule("D4","D4","triality")*b
            sage: [D4(fw).branch(A1xC2,rule=b1) for fw in D4.fundamental_weights()]
            [A1xC2(1,1,0),
            A1xC2(2,0,0) + A1xC2(2,0,1) + A1xC2(0,2,0),
            A1xC2(2,0,0) + A1xC2(0,0,1),
            A1xC2(1,1,0)]
        """
        return sage.combinat.root_system.branching_rules.maximal_subgroups(self.cartan_type())

    def maximal_subgroup(self, ct):
        """
        Return a branching rule or a list of branching rules.

        INPUT:

        - ``ct`` -- the Cartan type of a maximal subgroup of ``self``

        In rare cases where there is
        more than one maximal subgroup (up to outer automorphisms)
        with the given Cartan type, the function returns a list of
        branching rules.

        EXAMPLES::

            sage: WeylCharacterRing("E7").maximal_subgroup("A2")
            miscellaneous branching rule E7 => A2
            sage: WeylCharacterRing("E7").maximal_subgroup("A1")
            [iii branching rule E7 => A1, iv branching rule E7 => A1]

        For more information, see the related method :meth:`maximal_subgroups`.
        """
        return sage.combinat.root_system.branching_rules.maximal_subgroups(self.cartan_type(), mode='get_rule')[ct]

    class Element(CombinatorialFreeModule.Element):
        """
        A class for Weyl characters.
        """

        def cartan_type(self):
            """
            Return the Cartan type of ``self``.

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2")
                sage: A2([1,0,0]).cartan_type()
                ['A', 2]
            """
            return self.parent()._cartan_type

        def degree(self):
            """
            Return the degree of ``self``.

            This is the dimension of the associated module.

            EXAMPLES::

                sage: B3 = WeylCharacterRing(['B',3])
                sage: [B3(x).degree() for x in B3.fundamental_weights()]
                [7, 21, 8]
            """
            L = self.parent()._space
            return sum(L.weyl_dimension(k) * c for k, c in self)

        def branch(self, S, rule='default'):
            """
            Return the restriction of the character to the subalgebra.

            If no rule is specified, we will try to specify one.

            INPUT:

            - ``S`` -- a Weyl character ring for a Lie subgroup or subalgebra

            - ``rule`` -- a branching rule

            See :func:`~sage.combinat.root_system.branching_rules.branch_weyl_character`
            for more information about branching rules.

            EXAMPLES::

                sage: B3 = WeylCharacterRing(['B',3])
                sage: A2 = WeylCharacterRing(['A',2])
                sage: [B3(w).branch(A2,rule='levi') for w in B3.fundamental_weights()]
                [A2(0,0,0) + A2(1,0,0) + A2(0,0,-1),
                A2(0,0,0) + A2(1,0,0) + A2(1,1,0) + A2(1,0,-1) + A2(0,-1,-1) + A2(0,0,-1),
                A2(-1/2,-1/2,-1/2) + A2(1/2,-1/2,-1/2) + A2(1/2,1/2,-1/2) + A2(1/2,1/2,1/2)]
            """
            return sage.combinat.root_system.branching_rules.branch_weyl_character(self, self.parent(), S, rule=rule)

        def dual(self):
            """
            The involution that replaces a representation with
            its contragredient. (For Fusion rings, this is the
            conjugation map.)

            EXAMPLES::

                sage: A3 = WeylCharacterRing("A3", style='coroots')
                sage: A3(1,0,0)^2
                A3(0,1,0) + A3(2,0,0)
                sage: (A3(1,0,0)^2).dual()
                A3(0,1,0) + A3(0,0,2)
            """
            if not self.parent().cartan_type().is_irreducible():
                raise NotImplementedError("dual method is not implemented for reducible types")
            d = self.monomial_coefficients()
            WCR = self.parent()
            return sum(d[k] * WCR._element_constructor_(self.parent()._dual_helper(k))
                       for k in d)

        def highest_weight(self):
            """
            Return the parametrizing dominant weight
            of an irreducible character.

            This method is only available for basis elements.

            EXAMPLES::

                sage: G2 = WeylCharacterRing("G2", style='coroots')
                sage: [x.highest_weight() for x in [G2(1,0),G2(0,1)]]
                [(1, 0, -1), (2, -1, -1)]
            """
            if len(self.monomial_coefficients()) != 1:
                raise ValueError("fusion weight is valid for basis elements only")
            return self.leading_support()

        def __pow__(self, n):
            """
            Return the `n`-th power of ``self``.

            We override the method in :mod:`sage.monoids.monoids` since
            using the Brauer-Klimyk algorithm, it is more efficient to
            compute ``a*(a*(a*a))`` than ``(a*a)*(a*a)``.

            EXAMPLES::

                sage: B4 = WeylCharacterRing("B4",style='coroots')
                sage: spin = B4(0,0,0,1)
                sage: [spin^k for k in [0,1,3]]
                [B4(0,0,0,0), B4(0,0,0,1), 5*B4(0,0,0,1) + 4*B4(1,0,0,1) + 3*B4(0,1,0,1) + 2*B4(0,0,1,1) + B4(0,0,0,3)]
                sage: spin^-1
                Traceback (most recent call last):
                ...
                ValueError: cannot invert self (= B4(0,0,0,1))
                sage: x = 2 * B4.one(); x
                2*B4(0,0,0,0)
                sage: x^-3
                1/8*B4(0,0,0,0)
            """
            n = ZZ(n)
            if not n:
                return self.parent().one()
            if n < 0:
                self = ~self
                n = -n

            res = self
            for i in range(n - 1):
                res = self * res
            return res

        def is_irreducible(self):
            """
            Return whether ``self`` is an irreducible character.

            EXAMPLES::

                sage: B3 = WeylCharacterRing(['B',3])
                sage: [B3(x).is_irreducible() for x in B3.fundamental_weights()]
                [True, True, True]
                sage: sum(B3(x) for x in B3.fundamental_weights()).is_irreducible()
                False
            """
            return self.coefficients() == [1]

        @cached_method
        def symmetric_power(self, k):
            r"""
            Return the `k`-th symmetric power of ``self``.

            INPUT:

            - ``k`` -- nonnegative integer

            The algorithm is based on the
            identity `k h_k = \sum_{r=1}^k p_k h_{k-r}` relating the power-sum
            and complete symmetric polynomials. Applying this to the
            eigenvalues of an element of the parent Lie group in the
            representation ``self``, the `h_k` become symmetric powers and
            the `p_k` become Adams operations, giving an efficient recursive
            implementation.

            EXAMPLES::

                sage: B3 = WeylCharacterRing("B3",style='coroots')
                sage: spin = B3(0,0,1)
                sage: spin.symmetric_power(6)
                B3(0,0,0) + B3(0,0,2) + B3(0,0,4) + B3(0,0,6)
            """
            par = self.parent()
            if k == 0:
                return par.one()
            if k == 1:
                return self
            ret = par.zero()
            for r in range(1, k + 1):
                adam_r = self._adams_operator_helper(r)
                ret += par.linear_combination((par._product_helper(adam_r, l), c)
                                              for l, c in self.symmetric_power(k - r))
            m = ret.weight_multiplicities()
            dd = {key: val / k for key, val in m.items()}
            return self.parent().char_from_weights(dd)

        @cached_method
        def exterior_power(self, k):
            r"""
            Return the `k`-th exterior power of ``self``.

            INPUT:

            - ``k`` -- nonnegative integer

            The algorithm is based on the
            identity `k e_k = \sum_{r=1}^k (-1)^{k-1} p_k e_{k-r}` relating the
            power-sum and elementary symmetric polynomials. Applying this to
            the eigenvalues of an element of the parent Lie group in the
            representation ``self``, the `e_k` become exterior powers and
            the `p_k` become Adams operations, giving an efficient recursive
            implementation.

            EXAMPLES::

                sage: B3 = WeylCharacterRing("B3",style='coroots')
                sage: spin = B3(0,0,1)
                sage: spin.exterior_power(6)
                B3(1,0,0) + B3(0,1,0)
            """
            par = self.parent()
            if k == 0:
                return par.one()
            if k == 1:
                return self
            ret = par.zero()
            for r in range(1, k + 1):
                adam_r = self._adams_operator_helper(r)
                if is_even(r):
                    ret -= par.linear_combination((par._product_helper(adam_r, l), c) for (l, c) in self.exterior_power(k-r))
                else:
                    ret += par.linear_combination((par._product_helper(adam_r, l), c) for (l, c) in self.exterior_power(k-r))
            dd = {}
            m = ret.weight_multiplicities()
            for l in m:
                dd[l] = m[l]/k
            return self.parent().char_from_weights(dd)

        def adams_operator(self, r):
            """
            Return the `r`-th Adams operation of ``self``.

            INPUT:

            - ``r`` -- positive integer

            This is a virtual character,
            whose weights are the weights of ``self``, each multiplied by `r`.

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2")
                sage: A2(1,1,0).adams_operator(3)
                A2(2,2,2) - A2(3,2,1) + A2(3,3,0)
            """
            return self.parent().char_from_weights(self._adams_operator_helper(r))

        adams_operation = adams_operator

        def _adams_operator_helper(self, r):
            """
            Helper function for Adams operations.

            INPUT:

            - ``r`` -- positive integer

            Return the dictionary of weight multiplicities for the Adams
            operation, needed for internal use by symmetric and exterior powers.

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2")
                sage: A2(1,1,0)._adams_operator_helper(3)
                {(3, 3, 0): 1, (3, 0, 3): 1, (0, 3, 3): 1}
            """
            d = self.weight_multiplicities()
            return {r * key: val for key, val in d.items()}

        def symmetric_square(self):
            """
            Return the symmetric square of the character.

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2",style='coroots')
                sage: A2(1,0).symmetric_square()
                A2(2,0)
            """
            # Conceptually, this converts self to the weight ring,
            # computes its square there, and converts the result back.
            #
            # This implementation uses that this is a squaring (and not
            # a generic product) in the weight ring to optimize by
            # running only through pairs of weights instead of couples.
            c = self.weight_multiplicities()
            ckeys = list(c)
            d = {}
            for j in range(len(ckeys)):
                for i in range(j+1):
                    ci = ckeys[i]
                    cj = ckeys[j]
                    t = ci + cj
                    if i < j:
                        coef = c[ci]*c[cj]
                    else:
                        coef = c[ci]*(c[ci]+1)/2
                    if t in d:
                        d[t] += coef
                    else:
                        d[t] = coef
            for k in list(d):
                if d[k] == 0:
                    del d[k]
            return self.parent().char_from_weights(d)

        def exterior_square(self):
            """
            Return the exterior square of the character.

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2",style='coroots')
                sage: A2(1,0).exterior_square()
                A2(0,1)
            """
            c = self.weight_multiplicities()
            ckeys = list(c)
            d = {}
            for j in range(len(ckeys)):
                for i in range(j+1):
                    ci = ckeys[i]
                    cj = ckeys[j]
                    t = ci + cj
                    if i < j:
                        coef = c[ci]*c[cj]
                    else:
                        coef = c[ci]*(c[ci]-1)/2
                    if t in d:
                        d[t] += coef
                    else:
                        d[t] = coef
            for k in list(d):
                if d[k] == 0:
                    del d[k]
            return self.parent().char_from_weights(d)

        def frobenius_schur_indicator(self):
            r"""
            Return:

            - `1` if the representation is real (orthogonal)

            - `-1` if the representation is quaternionic (symplectic)

            - `0` if the representation is complex (not self dual)

            The Frobenius-Schur indicator of a character `\chi`
            of a compact group `G` is the Haar integral over the
            group of `\chi(g^2)`. Its value is 1, -1 or 0. This
            method computes it for irreducible characters of
            compact Lie groups by checking whether the symmetric
            and exterior square characters contain the trivial
            character.

            .. TODO::

                Try to compute this directly without actually calculating
                the full symmetric and exterior squares.

            EXAMPLES::

                sage: B2 = WeylCharacterRing("B2",style='coroots')
                sage: B2(1,0).frobenius_schur_indicator()
                 1
                sage: B2(0,1).frobenius_schur_indicator()
                 -1
            """
            if not self.is_irreducible():
                raise ValueError("Frobenius-Schur indicator is only valid for irreducible characters")
            z = self.parent()._space.zero()
            if self.symmetric_square().coefficient(z) != 0:
                return 1
            if self.exterior_square().coefficient(z) != 0:
                return -1
            return 0

        def weight_multiplicities(self):
            """
            Return the dictionary of weight multiplicities for the Weyl
            character ``self``.

            The character does not have to be irreducible.

            EXAMPLES::

                sage: B2 = WeylCharacterRing("B2",style='coroots')
                sage: B2(0,1).weight_multiplicities()
                {(-1/2, -1/2): 1, (-1/2, 1/2): 1, (1/2, -1/2): 1, (1/2, 1/2): 1}
            """
            return self.parent()._weight_multiplicities(self)

        def inner_product(self, other):
            """
            Compute the inner product with another character.

            The irreducible characters are an orthonormal basis with respect
            to the usual inner product of characters, interpreted as functions
            on a compact Lie group, by Schur orthogonality.

            INPUT:

            - ``other`` -- another character

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2")
                sage: [f1,f2] = A2.fundamental_weights()
                sage: r1 = A2(f1)*A2(f2); r1
                A2(1,1,1) + A2(2,1,0)
                sage: r2 = A2(f1)^3; r2
                A2(1,1,1) + 2*A2(2,1,0) + A2(3,0,0)
                sage: r1.inner_product(r2)
                3
            """
            return sum(self.coefficient(x) * other.coefficient(x)
                       for x in self.monomial_coefficients())

        def invariant_degree(self):
            """
            Return the multiplicity of the trivial representation in ``self``.

            Multiplicities of other irreducibles may be obtained
            using :meth:`multiplicity`.

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2",style='coroots')
                sage: rep = A2(1,0)^2*A2(0,1)^2; rep
                2*A2(0,0) + A2(0,3) + 4*A2(1,1) + A2(3,0) + A2(2,2)
                sage: rep.invariant_degree()
                2
            """
            return self.coefficient(self.parent().space()(0))

        def multiplicity(self, other):
            """
            Return the multiplicity of the irreducible ``other`` in ``self``.

            INPUT:

            - ``other`` -- an irreducible character

            EXAMPLES::

                sage: B2 = WeylCharacterRing("B2",style='coroots')
                sage: rep = B2(1,1)^2; rep
                B2(0,0) + B2(1,0) + 2*B2(0,2) + B2(2,0) + 2*B2(1,2) + B2(0,4) + B2(3,0) + B2(2,2)
                sage: rep.multiplicity(B2(0,2))
                2
            """
            if not other.is_irreducible():
                raise ValueError("{} is not irreducible".format(other))
            return self.coefficient(other.support()[0])


def irreducible_character_freudenthal(hwv, debug=False):
    r"""
    Return the dictionary of multiplicities for the irreducible
    character with highest weight `\lambda`.

    The weight multiplicities are computed by the Freudenthal multiplicity
    formula. The algorithm is based on recursion relation that is stated,
    for example, in Humphrey's book on Lie Algebras. The multiplicities are
    invariant under the Weyl group, so to compute them it would be sufficient
    to compute them for the weights in the positive Weyl chamber. However
    after some testing it was found to be faster to compute every
    weight using the recursion, since the use of the Weyl group is
    expensive in its current implementation.

    INPUT:

    - ``hwv`` -- a dominant weight in a weight lattice

    - ``L`` -- the ambient space

    EXAMPLES::

        sage: WeylCharacterRing("A2")(2,1,0).weight_multiplicities() # indirect doctest
        {(1, 1, 1): 2, (1, 2, 0): 1, (1, 0, 2): 1, (2, 1, 0): 1,
         (2, 0, 1): 1, (0, 1, 2): 1, (0, 2, 1): 1}
    """
    L = hwv.parent()
    rho = L.rho()
    mdict = {}
    current_layer = {hwv: 1}

    simple_roots = L.simple_roots()
    positive_roots = L.positive_roots()

    while current_layer:
        next_layer = {}
        for mu in current_layer:
            if current_layer[mu] != 0:
                mdict[mu] = current_layer[mu]
                for alpha in simple_roots:
                    next_layer[mu - alpha] = None
        if debug:
            print(next_layer)

        for mu in next_layer:
            if next_layer[mu] is None:
                accum = 0
                for alpha in positive_roots:
                    mu_plus_i_alpha = mu + alpha
                    while mu_plus_i_alpha in mdict:
                        accum += mdict[mu_plus_i_alpha] * (mu_plus_i_alpha).inner_product(alpha)
                        mu_plus_i_alpha += alpha
                if accum == 0:
                    next_layer[mu] = 0
                else:
                    hwv_plus_rho = hwv + rho
                    mu_plus_rho = mu + rho
                    next_layer[mu] = ZZ(2 * accum) / ZZ((hwv_plus_rho).inner_product(hwv_plus_rho) - (mu_plus_rho).inner_product(mu_plus_rho))
        current_layer = next_layer
    return mdict


class WeightRing(CombinatorialFreeModule):
    """
    The weight ring, which is the group algebra over a weight lattice.

    A Weyl character may be regarded as an element of the weight ring.
    In fact, an element of the weight ring is an element of the
    :class:`Weyl character ring <WeylCharacterRing>` if and only if it is
    invariant under the action of the Weyl group.

    The advantage of the weight ring over the Weyl character ring
    is that one may conduct calculations in the weight ring that
    involve sums of weights that are not Weyl group invariant.

    EXAMPLES::

        sage: A2 = WeylCharacterRing(['A',2])
        sage: a2 = WeightRing(A2)
        sage: wd = prod(a2(x/2)-a2(-x/2) for x in a2.space().positive_roots()); wd
        a2(-1,1,0) - a2(-1,0,1) - a2(1,-1,0) + a2(1,0,-1) + a2(0,-1,1) - a2(0,1,-1)
        sage: chi = A2([5,3,0]); chi
        A2(5,3,0)
        sage: a2(chi)
        a2(1,2,5) + 2*a2(1,3,4) + 2*a2(1,4,3) + a2(1,5,2) + a2(2,1,5)
        + 2*a2(2,2,4) + 3*a2(2,3,3) + 2*a2(2,4,2) + a2(2,5,1) + 2*a2(3,1,4)
        + 3*a2(3,2,3) + 3*a2(3,3,2) + 2*a2(3,4,1) + a2(3,5,0) + a2(3,0,5)
        + 2*a2(4,1,3) + 2*a2(4,2,2) + 2*a2(4,3,1) + a2(4,4,0) + a2(4,0,4)
        + a2(5,1,2) + a2(5,2,1) + a2(5,3,0) + a2(5,0,3) + a2(0,3,5)
        + a2(0,4,4) + a2(0,5,3)
        sage: a2(chi)*wd
        -a2(-1,3,6) + a2(-1,6,3) + a2(3,-1,6) - a2(3,6,-1) - a2(6,-1,3) + a2(6,3,-1)
        sage: sum((-1)^w.length()*a2([6,3,-1]).weyl_group_action(w) for w in a2.space().weyl_group())
        -a2(-1,3,6) + a2(-1,6,3) + a2(3,-1,6) - a2(3,6,-1) - a2(6,-1,3) + a2(6,3,-1)
        sage: a2(chi)*wd == sum((-1)^w.length()*a2([6,3,-1]).weyl_group_action(w) for w in a2.space().weyl_group())
        True
    """
    @staticmethod
    def __classcall__(cls, parent, prefix=None):
        """
        TESTS::

            sage: A3 = WeylCharacterRing("A3", style='coroots')
            sage: a3 = WeightRing(A3)
            sage: a3.cartan_type(), a3.base_ring(), a3.parent()
            (['A', 3], Integer Ring, The Weyl Character Ring of Type A3 with Integer Ring coefficients)
        """
        return super().__classcall__(cls, parent, prefix=prefix)

    def __init__(self, parent, prefix):
        """
        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: a2 = WeightRing(A2)
            sage: TestSuite(a2).run()

        TESTS::

            sage: A1xA1 = WeylCharacterRing("A1xA1")
            sage: a1xa1 = WeightRing(A1xA1)
            sage: TestSuite(a1xa1).run()
            sage: a1xa1.an_element()
            a1xa1(2,2,3,0)
        """
        self._parent = parent
        self._style = parent._style
        self._prefix = prefix
        self._space = parent._space
        self._cartan_type = parent._cartan_type
        self._rank = parent._rank
        self._origin = parent._origin
        self._base_ring = parent._base_ring
        if prefix is None:
            # TODO: refactor this fragile logic into CartanType's
            if self._parent._prefix.replace('x', '_').isupper():
                # The 'x' workaround above is to support reducible Cartan types like 'A1xB2'
                prefix = self._parent._prefix.lower()
            elif self._parent._prefix.islower():
                prefix = self._parent._prefix.upper()
            else:
                # TODO: this only works for irreducible Cartan types!
                prefix = (self._cartan_type[0].lower() + str(self._rank))
        self._prefix = prefix
        category = AlgebrasWithBasis(self._base_ring).Commutative()
        CombinatorialFreeModule.__init__(self, self._base_ring, self._space, category=category)

    def _repr_(self):
        """
        EXAMPLES::

            sage: P.<q>=QQ[]
            sage: G2 = WeylCharacterRing(['G',2], base_ring = P)
            sage: WeightRing(G2) # indirect doctest
            The Weight ring attached to The Weyl Character Ring of Type G2 with Univariate Polynomial Ring in q over Rational Field coefficients
        """
        return "The Weight ring attached to %s" % self._parent

    def __call__(self, *args):
        """
        Construct an element of ``self``.

        The input can either be an object that can be coerced or
        converted into ``self`` (an element of ``self``, of the base
        ring, of the weight ring), or a dominant weight. In the later
        case, the basis element indexed by that weight is returned.

        To specify the weight, you may give it explicitly. Alternatively,
        you may give a tuple of integers. Normally these are the
        components of the vector in the standard realization of
        the weight lattice as a vector space. Alternatively, if
        the ring is constructed with style='coroots', you may
        specify the weight by giving a set of integers, one for each
        fundamental weight; the weight is then the linear combination
        of the fundamental weights with these coefficients.

        As a syntactical shorthand, for tuples of length at least two,
        the parenthesis may be omitted.

        EXAMPLES::

            sage: a2 = WeightRing(WeylCharacterRing(['A',2]))
            sage: a2(-1)
            -a2(0,0,0)
        """
        # The purpose of this __call__ method is only to handle the
        # syntactical shorthand; otherwise it just delegates the work
        # to the coercion model, which itself will call
        # _element_constructor_ if the input is made of exactly one
        # object which can't be coerced into self
        if len(args) > 1:
            args = (args,)
        return super().__call__(*args)

    def _element_constructor_(self, weight):
        """
        Construct a monomial from a weight.

        INPUT:

        - ``weight`` -- an element of the weight space, or a tuple

        This method is responsible for constructing an appropriate
        weight from the data in ``weight``, and then return the
        monomial indexed by that weight. See :meth:`__call__` and
        :meth:`sage.combinat.root_system.ambient_space.AmbientSpace.from_vector`.

        TESTS::

            sage: A2 = WeylCharacterRing("A2")
            sage: A2._element_constructor_([2,1,0])
            A2(2,1,0)
        """
        weight = self._space.from_vector_notation(weight, style=self._style)
        return self.monomial(weight)

    def product_on_basis(self, a, b):
        """
        Return the product of basis elements indexed by ``a`` and ``b``.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: a2 = WeightRing(A2)
            sage: a2(1,0,0) * a2(0,1,0) # indirect doctest
            a2(1,1,0)
        """
        return self(a + b)

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: A3 = WeylCharacterRing("A3")
            sage: a3 = WeightRing(A3)
            sage: a3.some_elements()
            [a3(1,0,0,0), a3(1,1,0,0), a3(1,1,1,0)]
        """
        return [self.monomial(x) for x in self.fundamental_weights()]

    def one_basis(self):
        """
        Return the index of `1`.

        EXAMPLES::

            sage: A3 = WeylCharacterRing("A3")
            sage: WeightRing(A3).one_basis()
            (0, 0, 0, 0)
            sage: WeightRing(A3).one()
            a3(0,0,0,0)
        """
        return self._space.zero()

    def parent(self):
        """
        Return the parent Weyl character ring.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: a2 = WeightRing(A2)
            sage: a2.parent()
            The Weyl Character Ring of Type A2 with Integer Ring coefficients
            sage: a2.parent() == A2
            True
        """
        return self._parent

    def weyl_character_ring(self):
        """
        Return the parent Weyl Character Ring.

        A synonym for ``self.parent()``.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: a2 = WeightRing(A2)
            sage: a2.weyl_character_ring()
            The Weyl Character Ring of Type A2 with Integer Ring coefficients
        """
        return self._parent

    def cartan_type(self):
        """
        Return the Cartan type.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: WeightRing(A2).cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def space(self):
        """
        Return the weight space realization associated to ``self``.

        EXAMPLES::

            sage: E8 = WeylCharacterRing(['E',8])
            sage: e8 = WeightRing(E8)
            sage: e8.space()
            Ambient space of the Root system of type ['E', 8]
        """
        return self._space

    def fundamental_weights(self):
        """
        Return the fundamental weights.

        EXAMPLES::

            sage: WeightRing(WeylCharacterRing("G2")).fundamental_weights()
            Finite family {1: (1, 0, -1), 2: (2, -1, -1)}
        """
        return self._space.fundamental_weights()

    def simple_roots(self):
        """
        Return the simple roots.

        EXAMPLES::

            sage: WeightRing(WeylCharacterRing("G2")).simple_roots()
            Finite family {1: (0, 1, -1), 2: (1, -2, 1)}
        """
        return self._space.simple_roots()

    def positive_roots(self):
        """
        Return the positive roots.

        EXAMPLES::

            sage: WeightRing(WeylCharacterRing("G2")).positive_roots()
            [(0, 1, -1), (1, -2, 1), (1, -1, 0), (1, 0, -1), (1, 1, -2), (2, -1, -1)]
        """
        return self._space.positive_roots()

    def wt_repr(self, wt):
        r"""
        Return a string representing the irreducible character with
        highest weight vector ``wt``.

        Uses coroot notation if the associated
        Weyl character ring is defined with ``style="coroots"``.

        EXAMPLES::

            sage: G2 = WeylCharacterRing("G2")
            sage: [G2.ambient().wt_repr(x) for x in G2.fundamental_weights()]
            ['g2(1,0,-1)', 'g2(2,-1,-1)']
            sage: G2 = WeylCharacterRing("G2",style='coroots')
            sage: [G2.ambient().wt_repr(x) for x in G2.fundamental_weights()]
            ['g2(1,0)', 'g2(0,1)']
        """
        return self._prefix + self.parent()._wt_repr(wt)

    def _repr_term(self, t):
        """
        Representation of the monomial corresponding to a weight ``t``.

        EXAMPLES::

            sage: G2 = WeylCharacterRing("G2")
            sage: g2 = WeightRing(G2)
            sage: [g2(x) for x in g2.fundamental_weights()] # indirect doctest
            [g2(1,0,-1), g2(2,-1,-1)]
        """
        return self.wt_repr(t)

    class Element(CombinatorialFreeModule.Element):
        """
        A class for weight ring elements.
        """

        def cartan_type(self):
            """
            Return the Cartan type.

            EXAMPLES::

                sage: A2 = WeylCharacterRing("A2")
                sage: a2 = WeightRing(A2)
                sage: a2([0,1,0]).cartan_type()
                ['A', 2]
            """
            return self.parent()._cartan_type

        def weyl_group_action(self, w):
            """
            Return the action of the Weyl group element ``w`` on ``self``.

            EXAMPLES::

                sage: G2 = WeylCharacterRing(['G',2])
                sage: g2 = WeightRing(G2)
                sage: L = g2.space()
                sage: [fw1, fw2] = L.fundamental_weights()
                sage: sum(g2(fw2).weyl_group_action(w) for w in L.weyl_group())
                2*g2(-2,1,1) + 2*g2(-1,-1,2) + 2*g2(-1,2,-1) + 2*g2(1,-2,1) + 2*g2(1,1,-2) + 2*g2(2,-1,-1)
            """
            return self.map_support(w.action)

        def character(self):
            """
            Assuming that ``self`` is invariant under the Weyl group, this will
            express it as a linear combination of characters. If ``self`` is
            not Weyl group invariant, this method will not terminate.

            EXAMPLES::

                sage: A2 = WeylCharacterRing(['A',2])
                sage: a2 = WeightRing(A2)
                sage: W = a2.space().weyl_group()
                sage: mu = a2(2,1,0)
                sage: nu = sum(mu.weyl_group_action(w) for w in W) ; nu
                a2(1,2,0) + a2(1,0,2) + a2(2,1,0) + a2(2,0,1) + a2(0,1,2) + a2(0,2,1)
                sage: nu.character()
                -2*A2(1,1,1) + A2(2,1,0)
            """
            return self.parent().parent().char_from_weights(self.monomial_coefficients())

        def scale(self, k):
            """
            Multiply a weight by `k`.

            The operation is extended by linearity to the weight ring.

            INPUT:

            - ``k`` -- nonzero integer

            EXAMPLES::

                sage: g2 = WeylCharacterRing("G2",style='coroots').ambient()
                sage: g2(2,3).scale(2)
                g2(4,6)
            """
            if k == 0:
                raise ValueError("parameter must be nonzero")
            d1 = self.monomial_coefficients()
            d2 = {k * mu: coeff for mu, coeff in d1.items()}
            return self.parent()._from_dict(d2)

        def shift(self, mu):
            r"""
            Add `\mu` to any weight.

            Extended by linearity to the weight ring.

            INPUT:

            - ``mu`` -- a weight

            EXAMPLES::

                sage: g2 = WeylCharacterRing("G2",style='coroots').ambient()
                sage: [g2(1,2).shift(fw) for fw in g2.fundamental_weights()]
                [g2(2,2), g2(1,3)]
            """
            d1 = self.monomial_coefficients()
            d2 = {mu + nu: val for nu, val in d1.items()}
            return self.parent()._from_dict(d2)

        def demazure(self, w, debug=False):
            r"""
            Return the result of applying the Demazure operator `\partial_w`
            to ``self``.

            INPUT:

            - ``w`` -- a Weyl group element, or its reduced word

            If `w = s_i` is a simple reflection, the operation `\partial_w`
            sends the weight `\lambda` to

            .. MATH::

                \frac{\lambda - s_i \cdot \lambda + \alpha_i}{1 + \alpha_i},

            where the numerator is divisible the denominator in the weight
            ring. This is extended by multiplicativity to all `w` in the
            Weyl group.

            EXAMPLES::

                sage: B2 = WeylCharacterRing("B2",style='coroots')
                sage: b2 = WeightRing(B2)
                sage: b2(1,0).demazure([1])
                b2(1,0) + b2(-1,2)
                sage: b2(1,0).demazure([2])
                b2(1,0)
                sage: r = b2(1,0).demazure([1,2]); r
                b2(1,0) + b2(-1,2)
                sage: r.demazure([1])
                b2(1,0) + b2(-1,2)
                sage: r.demazure([2])
                b2(0,0) + b2(1,0) + b2(1,-2) + b2(-1,2)
            """
            if isinstance(w, list):
                word = w
            else:
                word = w.reduced_word()
            d1 = self.monomial_coefficients()
            d = {}
            alphacheck = self.parent()._space.simple_coroots()
            for v in d1:
                d[tuple(v.inner_product(alphacheck[j])
                        for j in self.parent().space().index_set())] = d1[v]
            return self.parent()._from_dict(self.parent().parent()._demazure_helper(d, word, debug=debug))

        def demazure_lusztig(self, i, v):
            r"""
            Return the result of applying the Demazure-Lusztig operator
            `T_i` to ``self``.

            INPUT:

            - ``i`` -- an element of the index set (or a reduced word or
              Weyl group element)
            - ``v`` -- an element of the base ring

            If `R` is the parent WeightRing, the Demazure-Lusztig operator
            `T_i` is the linear map `R \to R` that sends (for a weight
            `\lambda`) `R(\lambda)` to

            .. MATH::

                (R(\alpha_i)-1)^{-1} \bigl(R(\lambda) - R(s_i\lambda)
                - v(R(\lambda) - R(\alpha_i + s_i \lambda)) \bigr)

            where the numerator is divisible by the denominator in `R`.
            The Demazure-Lusztig operators give a representation of the
            Iwahori--Hecke algebra associated to the Weyl group. See

            * Lusztig, Equivariant `K`-theory and representations of Hecke
              algebras, Proc. Amer. Math. Soc. 94 (1985), no. 2, 337-342.
            * Cherednik, *Nonsymmetric Macdonald polynomials*. IMRN 10,
              483-515 (1995).

            In the examples, we confirm the braid and quadratic relations
            for type `B_2`.

            EXAMPLES::

                sage: P.<v> = PolynomialRing(QQ)
                sage: B2 = WeylCharacterRing("B2",style='coroots',base_ring=P); b2 = B2.ambient()
                sage: def T1(f): return f.demazure_lusztig(1,v)
                sage: def T2(f): return f.demazure_lusztig(2,v)
                sage: T1(T2(T1(T2(b2(1,-1)))))
                (v^2-v)*b2(0,-1) + v^2*b2(-1,1)
                sage: [T1(T1(f))==(v-1)*T1(f)+v*f for f in [b2(0,0), b2(1,0), b2(2,3)]]
                [True, True, True]
                sage: [T1(T2(T1(T2(b2(i,j))))) == T2(T1(T2(T1(b2(i,j))))) for i in [-2..2] for j in [-1,1]]
                [True, True, True, True, True, True, True, True, True, True]

            Instead of an index `i` one may use a reduced word or
            Weyl group element::

                sage: b2(1,0).demazure_lusztig([2,1],v)==T2(T1(b2(1,0)))
                True
                sage: W = B2.space().weyl_group(prefix='s')
                sage: [s1,s2]=W.simple_reflections()
                sage: b2(1,0).demazure_lusztig(s2*s1,v)==T2(T1(b2(1,0)))
                True
            """
            if i in self.parent().space().index_set():
                rho = self.parent().space().from_vector_notation(self.parent().space().rho(), style='coroots')
                inv = self.scale(-1)
                return (-inv.shift(-rho).demazure([i]).shift(rho) + v * inv.demazure([i])).scale(-1)
            elif isinstance(i, list):
                if not i:
                    return self
                elif len(i) == 1:
                    return self.demazure_lusztig(i[0], v)
                else:
                    return self.demazure_lusztig(i[1:], v).demazure_lusztig(i[:1], v)
            else:
                try:
                    return self.demazure_lusztig(i.reduced_word(), v)
                except Exception:
                    raise ValueError("unknown index {}".format(i))
