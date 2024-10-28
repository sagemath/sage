# sage.doctest: needs sage.combinat sage.modules
"""
Rational Cherednik Algebras
"""
# ****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.algebras import Algebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.combinat.root_system.root_system import RootSystem
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ


class RationalCherednikAlgebra(CombinatorialFreeModule):
    r"""
    A rational Cherednik algebra.

    Let `k` be a field. Let `W` be a complex reflection group acting on
    a vector space `\mathfrak{h}` (over `k`). Let `\mathfrak{h}^*` denote
    the corresponding dual vector space. Let `\cdot` denote the
    natural action of `w` on `\mathfrak{h}` and `\mathfrak{h}^*`. Let
    `\mathcal{S}` denote the set of reflections of `W` and  `\alpha_s`
    and `\alpha_s^{\vee}` are the associated root and coroot of `s`. Let
    `c = (c_s)_{s \in W}` such that `c_s = c_{tst^{-1}}` for all `t \in W`.

    The *rational Cherednik algebra* is the `k`-algebra
    `H_{c,t}(W) = T(\mathfrak{h} \oplus \mathfrak{h}^*) \otimes kW` with
    parameters `c, t \in k` that is subject to the relations:

    .. MATH::

        \begin{aligned}
        w \alpha & = (w \cdot \alpha) w,
        \\ \alpha^{\vee} w & = w (w^{-1} \cdot \alpha^{\vee}),
        \\ \alpha \alpha^{\vee} & = \alpha^{\vee} \alpha
        + t \langle \alpha^{\vee}, \alpha \rangle
        + \sum_{s \in \mathcal{S}} c_s \frac{\langle \alpha^{\vee},
        \alpha_s \rangle \langle \alpha^{\vee}_s, \alpha \rangle}{
        \langle \alpha^{\vee}, \alpha \rangle} s,
        \end{aligned}

    where `w \in W` and `\alpha \in \mathfrak{h}` and
    `\alpha^{\vee} \in \mathfrak{h}^*`.

    INPUT:

    - ``ct`` -- a finite Cartan type
    - ``c`` -- the parameters `c_s` given as an element or a tuple, where
      the first entry is the one for the long roots and (for
      non-simply-laced types) the second is for the short roots
    - ``t`` -- the parameter `t`
    - ``base_ring`` -- (optional) the base ring
    - ``prefix`` -- (default: ``('a', 's', 'ac')``) the prefixes

    .. TODO::

        Implement a version for complex reflection groups.

    REFERENCES:

    - [GGOR2003]_
    - [EM2001]_
    """
    @staticmethod
    def __classcall_private__(cls, ct, c=1, t=None, base_ring=None, prefix=('a', 's', 'ac')):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: R1 = algebras.RationalCherednik(['B',2], 1, 1, QQ)
            sage: R2 = algebras.RationalCherednik(CartanType(['B',2]), [1,1], 1, QQ, ('a', 's', 'ac'))
            sage: R1 is R2
            True
        """
        ct = CartanType(ct)
        if not ct.is_finite():
            raise ValueError("the Cartan type must be finite")
        if base_ring is None:
            if t is None:
                base_ring = QQ
            else:
                base_ring = t.parent()
        if t is None:
            t = base_ring.one()
        else:
            t = base_ring(t)

        # Normalize the parameter c
        if isinstance(c, (tuple, list)):
            if ct.is_simply_laced():
                if len(c) != 1:
                    raise ValueError("1 parameter c_s must be given for simply-laced types")
                c = (base_ring(c[0]),)
            else:
                if len(c) != 2:
                    raise ValueError("2 parameters c_s must be given for non-simply-laced types")
                c = (base_ring(c[0]), base_ring(c[1]))
        else:
            c = base_ring(c)
            if ct.is_simply_laced():
                c = (c,)
            else:
                c = (c, c)

        return super().__classcall__(cls, ct, c, t, base_ring, tuple(prefix))

    def __init__(self, ct, c, t, base_ring, prefix):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: k = QQ['c,t']
            sage: R = algebras.RationalCherednik(['A',2], k.gen(0), k.gen(1))
            sage: TestSuite(R).run()  # long time
        """
        self._c = c
        self._t = t
        self._cartan_type = ct
        self._weyl = RootSystem(ct).root_lattice().weyl_group(prefix=prefix[1])
        self._hd = IndexedFreeAbelianMonoid(ct.index_set(), prefix=prefix[0],
                                            bracket=False)
        self._h = IndexedFreeAbelianMonoid(ct.index_set(), prefix=prefix[2],
                                           bracket=False)
        indices = DisjointUnionEnumeratedSets([self._hd, self._weyl, self._h])
        CombinatorialFreeModule.__init__(self, base_ring, indices,
                                         category=Algebras(base_ring).WithBasis().Graded(),
                                         sorting_key=self._genkey)

    def _genkey(self, t):
        r"""
        Construct a key for comparison for a term indexed by ``t``.

        The key we create is the tuple in the following order:

        - overall degree
        - length of the Weyl group element
        - the Weyl group element
        - the element of `\mathfrak{h}`
        - the element of `\mathfrak{h}^*`

        EXAMPLES::

            sage: R = algebras.RationalCherednik(['A',2], 1, 1, QQ)
            sage: R.an_element()**2 # indirect doctest
            9*ac1^2 + 10*I + 6*a1*ac1 + 6*s1 + 3/2*s2 + 3/2*s1*s2*s1 + a1^2
        """
        return (self.degree_on_basis(t), t[1].length(), t[1], str(t[0]), str(t[2]))

    @lazy_attribute
    def _reflections(self):
        """
        A dictionary of reflections to a pair of the associated root
        and coroot.

        EXAMPLES::

            sage: R = algebras.RationalCherednik(['B',2], [1,2], 1, QQ)
            sage: [R._reflections[k] for k in sorted(R._reflections, key=str)]
            [(alpha[1], alphacheck[1], 1),
             (alpha[1] + alpha[2], 2*alphacheck[1] + alphacheck[2], 2),
             (alpha[2], alphacheck[2], 2),
             (alpha[1] + 2*alpha[2], alphacheck[1] + alphacheck[2], 1)]
        """
        d = {}
        for r in RootSystem(self._cartan_type).root_lattice().positive_roots():
            s = self._weyl.from_reduced_word(r.associated_reflection())
            if r.is_short_root():
                c = self._c[1]
            else:
                c = self._c[0]
            d[s] = (r, r.associated_coroot(), c)
        return d

    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: RationalCherednikAlgebra(['A',4], 2, 1, QQ)
            Rational Cherednik Algebra of type ['A', 4] with c=2 and t=1
             over Rational Field
            sage: algebras.RationalCherednik(['B',2], [1,2], 1, QQ)
            Rational Cherednik Algebra of type ['B', 2] with c_L=1 and c_S=2
             and t=1 over Rational Field
        """
        ret = "Rational Cherednik Algebra of type {} with ".format(self._cartan_type)
        if self._cartan_type.is_simply_laced():
            ret += "c={}".format(self._c[0])
        else:
            ret += "c_L={} and c_S={}".format(*self._c)
        return ret + " and t={} over {}".format(self._t, self.base_ring())

    def _repr_term(self, t):
        """
        Return a string representation of the term indexed by ``t``.

        EXAMPLES::

            sage: R = algebras.RationalCherednik(['A',2], 1, 1, QQ)
            sage: R.an_element() # indirect doctest
            3*ac1 + 2*s1 + a1
            sage: R.one() # indirect doctest
            I
        """
        r = []
        if t[0] != self._hd.one():
            r.append(t[0])
        if t[1] != self._weyl.one():
            r.append(t[1])
        if t[2] != self._h.one():
            r.append(t[2])
        if not r:
            return 'I'
        return '*'.join(repr(x) for x in r)

    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: R = algebras.RationalCherednik(['A',2], 1, 1, QQ)
            sage: list(R.algebra_generators())
            [a1, a2, s1, s2, ac1, ac2]
        """
        keys = ['a'+str(i) for i in self._cartan_type.index_set()]
        keys += ['s'+str(i) for i in self._cartan_type.index_set()]
        keys += ['ac'+str(i) for i in self._cartan_type.index_set()]

        def gen_map(k):
            if k[0] == 's':
                i = int(k[1:])
                return self.monomial((self._hd.one(),
                                       self._weyl.group_generators()[i],
                                       self._h.one()))
            if k[1] == 'c':
                i = int(k[2:])
                return self.monomial((self._hd.one(),
                                       self._weyl.one(),
                                       self._h.monoid_generators()[i]))

            i = int(k[1:])
            return self.monomial((self._hd.monoid_generators()[i],
                                   self._weyl.one(),
                                   self._h.one()))
        return Family(keys, gen_map)

    @cached_method
    def one_basis(self):
        """
        Return the index of the element `1`.

        EXAMPLES::

            sage: R = algebras.RationalCherednik(['A',2], 1, 1, QQ)
            sage: R.one_basis()
            (1, 1, 1)
        """
        return (self._hd.one(), self._weyl.one(), self._h.one())

    def product_on_basis(self, left, right):
        r"""
        Return ``left`` multiplied by ``right`` in ``self``.

        EXAMPLES::

            sage: R = algebras.RationalCherednik(['A',2], 1, 1, QQ)
            sage: a2 = R.algebra_generators()['a2']
            sage: ac1 = R.algebra_generators()['ac1']
            sage: a2 * ac1  # indirect doctest
            a2*ac1
            sage: ac1 * a2
            -I + a2*ac1 - s1 - s2 + 1/2*s1*s2*s1
            sage: x = R.an_element()
            sage: [y * x for y in R.some_elements()]
            [0,
             3*ac1 + 2*s1 + a1,
             9*ac1^2 + 10*I + 6*a1*ac1 + 6*s1 + 3/2*s2 + 3/2*s1*s2*s1 + a1^2,
             3*a1*ac1 + 2*a1*s1 + a1^2,
             3*a2*ac1 + 2*a2*s1 + a1*a2,
             3*s1*ac1 + 2*I - a1*s1,
             3*s2*ac1 + 2*s2*s1 + a1*s2 + a2*s2,
             3*ac1^2 - 2*s1*ac1 + 2*I + a1*ac1 + 2*s1 + 1/2*s2 + 1/2*s1*s2*s1,
             3*ac1*ac2 + 2*s1*ac1 + 2*s1*ac2 - I + a1*ac2 - s1 - s2 + 1/2*s1*s2*s1]
            sage: [x * y for y in R.some_elements()]
            [0,
             3*ac1 + 2*s1 + a1,
             9*ac1^2 + 10*I + 6*a1*ac1 + 6*s1 + 3/2*s2 + 3/2*s1*s2*s1 + a1^2,
             6*I + 3*a1*ac1 + 6*s1 + 3/2*s2 + 3/2*s1*s2*s1 - 2*a1*s1 + a1^2,
             -3*I + 3*a2*ac1 - 3*s1 - 3*s2 + 3/2*s1*s2*s1 + 2*a1*s1 + 2*a2*s1 + a1*a2,
             -3*s1*ac1 + 2*I + a1*s1,
             3*s2*ac1 + 3*s2*ac2 + 2*s1*s2 + a1*s2,
             3*ac1^2 + 2*s1*ac1 + a1*ac1,
             3*ac1*ac2 + 2*s1*ac2 + a1*ac2]
        """
        # Make copies of the internal dictionaries
        dl = dict(left[2]._monomial)
        dr = dict(right[0]._monomial)

        # If there is nothing to commute
        if not dl and not dr:
            return self.monomial((left[0], left[1] * right[1], right[2]))

        R = self.base_ring()
        I = self._cartan_type.index_set()
        P = PolynomialRing(R, 'x', len(I))
        G = P.gens()
        gens_dict = {a:G[i] for i,a in enumerate(I)}
        Q = RootSystem(self._cartan_type).root_lattice()
        alpha = Q.simple_roots()
        alphacheck = Q.simple_coroots()

        def commute_w_hd(w, al): # al is given as a dictionary
            ret = P.one()
            for k in al:
                x = sum(c * gens_dict[i] for i, c in alpha[k].weyl_action(w))
                ret *= x**al[k]
            ret = ret.monomial_coefficients()
            for k in ret:
                yield (self._hd({I[i]: e for i, e in enumerate(k) if e != 0}), ret[k])

        # Do Lac Ra if they are both non-trivial
        if dl and dr:
            il = next(iter(dl.keys()))
            ir = next(iter(dr.keys()))

            # Compute the commutator
            terms = self._product_coroot_root(il, ir)

            # remove the generator from the elements
            dl[il] -= 1
            if dl[il] == 0:
                del dl[il]
            dr[ir] -= 1
            if dr[ir] == 0:
                del dr[ir]

            # We now commute right roots past the left reflections: s Ra = Ra' s
            cur = self._from_dict({(hd, s*right[1], right[2]): c * cc
                                    for s,c in terms
                                    for hd, cc in commute_w_hd(s, dr)})
            cur = self.monomial((left[0], left[1], self._h(dl))) * cur

            # Add back in the commuted h and hd elements
            rem = self.monomial((left[0], left[1], self._h(dl)))
            rem = rem * self.monomial((self._hd({ir:1}), self._weyl.one(),
                                        self._h({il:1})))
            rem = rem * self.monomial((self._hd(dr), right[1], right[2]))

            return cur + rem

        if dl:
            # We have La Ls Lac Rs Rac,
            #   so we must commute Lac Rs = Rs Lac'
            #   and obtain La (Ls Rs) (Lac' Rac)
            ret = P.one()
            for k in dl:
                x = sum(c * gens_dict[i]
                        for i,c in alphacheck[k].weyl_action(right[1].reduced_word(),
                                                             inverse=True))
                ret *= x**dl[k]
            ret = ret.monomial_coefficients()
            w = left[1]*right[1]
            return self._from_dict({(left[0], w,
                                      self._h({I[i]: e for i,e in enumerate(k)
                                               if e != 0}) * right[2]
                                     ): ret[k]
                                     for k in ret})

        # Otherwise dr is non-trivial and we have La Ls Ra Rs Rac,
        #   so we must commute Ls Ra = Ra' Ls
        w = left[1]*right[1]
        return self._from_dict({(left[0] * hd, w, right[2]): c
                                 for hd, c in commute_w_hd(left[1], dr)})

    @cached_method
    def _product_coroot_root(self, i, j):
        r"""
        Return the product `\alpha^{\vee}_i \alpha_j`.

        EXAMPLES::

            sage: k = QQ['c,t']
            sage: R = algebras.RationalCherednik(['A',3], k.gen(0), k.gen(1))
            sage: sorted(R._product_coroot_root(1, 1))
            [(s1, 2*c),
             (s1*s2*s1, 1/2*c),
             (s1*s2*s3*s2*s1, 1/2*c),
             (1, 2*t),
             (s3, 0),
             (s2, 1/2*c),
             (s2*s3*s2, 1/2*c)]

            sage: sorted(R._product_coroot_root(1, 2))
            [(s1, -c),
             (s1*s2*s1, 1/2*c),
             (s1*s2*s3*s2*s1, 0),
             (1, -t),
             (s3, 0),
             (s2, -c),
             (s2*s3*s2, -1/2*c)]

            sage: sorted(R._product_coroot_root(1, 3))
            [(s1, 0),
             (s1*s2*s1, -1/2*c),
             (s1*s2*s3*s2*s1, 1/2*c),
             (1, 0),
             (s3, 0),
             (s2, 1/2*c),
             (s2*s3*s2, -1/2*c)]
        """
        Q = RootSystem(self._cartan_type).root_lattice()
        ac = Q.simple_coroot(i)
        al = Q.simple_root(j)

        R = self.base_ring()
        terms = [(self._weyl.one(), self._t * R(ac.scalar(al)))]
        for s in self._reflections:
            # p[0] is the root, p[1] is the coroot, p[2] the value c_s
            pr, pc, c = self._reflections[s]
            terms.append((s, c * R(ac.scalar(pr) * pc.scalar(al)
                                    / pc.scalar(pr))))
        return tuple(terms)

    def degree_on_basis(self, m):
        """
        Return the degree on the monomial indexed by ``m``.

        EXAMPLES::

            sage: R = algebras.RationalCherednik(['A',2], 1, 1, QQ)
            sage: [R.degree_on_basis(g.leading_support())
            ....:  for g in R.algebra_generators()]
            [1, 1, 0, 0, -1, -1]
        """
        return m[0].length() - m[2].length()

    @cached_method
    def trivial_idempotent(self):
        r"""
        Return the trivial idempotent of ``self``.

        Let `e = |W|^{-1} \sum_{w \in W} w` is the trivial idempotent.
        Thus `e^2 = e` and `eW = We`. The trivial idempotent is used
        in the construction of the spherical Cherednik algebra from
        the rational Cherednik algebra by `U_{c,t}(W) = e H_{c,t}(W) e`.

        EXAMPLES::

            sage: R = algebras.RationalCherednik(['A',2], 1, 1, QQ)
            sage: R.trivial_idempotent()
            1/6*I + 1/6*s1 + 1/6*s2 + 1/6*s2*s1 + 1/6*s1*s2 + 1/6*s1*s2*s1
        """
        coeff = self.base_ring()(~self._weyl.cardinality())
        hd_one = self._hd.one() # root - a
        h_one = self._h.one() # coroot - ac
        return self._from_dict({(hd_one, w, h_one): coeff for w in self._weyl},
                               remove_zeros=False)

    @cached_method
    def deformed_euler(self):
        """
        Return the element `eu_k`.

        EXAMPLES::

            sage: R = algebras.RationalCherednik(['A',2], 1, 1, QQ)
            sage: R.deformed_euler()
            2*I + 2/3*a1*ac1 + 1/3*a1*ac2 + 1/3*a2*ac1 + 2/3*a2*ac2
             + s1 + s2 + s1*s2*s1
        """
        I = self._cartan_type.index_set()
        G = self.algebra_generators()
        cm = ~CartanMatrix(self._cartan_type)
        n = len(I)
        ac = [G['ac'+str(i)] for i in I]
        la = [sum(cm[i,j]*G['a'+str(I[i])] for i in range(n)) for j in range(n)]
        return self.sum(ac[i]*la[i] for i in range(n))

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: R = algebras.RationalCherednik(['A',2], 1, 1, QQ)
            sage: R.an_element()
            3*ac1 + 2*s1 + a1
        """
        G = self.algebra_generators()
        i = str(self._cartan_type.index_set()[0])
        return G['a'+i] + 2*G['s'+i] + 3*G['ac'+i]

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: R = algebras.RationalCherednik(['A',2], 1, 1, QQ)
            sage: R.some_elements()
            [0, I, 3*ac1 + 2*s1 + a1, a1, a2, s1, s2, ac1, ac2]
        """
        ret = [self.zero(), self.one(), self.an_element()]
        ret += list(self.algebra_generators())
        return ret
