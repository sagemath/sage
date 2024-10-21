# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Affine factorization crystal of type `A`
"""
#*****************************************************************************
#  Copyright (C) 2014 Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import lazy_import
from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.crystals import CrystalMorphism
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.homset import Hom
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.rsk import RSK

lazy_import('sage.combinat.root_system.weyl_group', 'WeylGroup')


class AffineFactorizationCrystal(UniqueRepresentation, Parent):
    r"""
    The crystal on affine factorizations with a cut-point, as introduced
    by [MS2015]_.

    INPUT:

    - ``w`` -- an element in an (affine) Weyl group or a skew shape of `k`-bounded partitions (if `k` was specified)

    - ``n`` -- the number of factors in the factorization

    - ``x`` -- (default: ``None``) the cut point; if not specified it is determined as the minimal missing residue in ``w``

    - ``k`` -- (default: ``None``) positive integer, specifies that ``w`` is `k`-bounded or a `k+1`-core when specified

    EXAMPLES::

        sage: W = WeylGroup(['A',3,1], prefix='s')
        sage: w = W.from_reduced_word([2,3,2,1])
        sage: B = crystals.AffineFactorization(w,3); B
        Crystal on affine factorizations of type A2 associated to s2*s3*s2*s1
        sage: B.list()
        [(1, s2, s3*s2*s1),
         (1, s3*s2, s3*s1),
         (1, s3*s2*s1, s3),
         (s3, s2, s3*s1),
         (s3, s2*s1, s3),
         (s3*s2, s1, s3),
         (s3*s2*s1, 1, s3),
         (s3*s2*s1, s3, 1),
         (s3*s2, 1, s3*s1),
         (s3*s2, s3, s1),
         (s3*s2, s3*s1, 1),
         (s2, 1, s3*s2*s1),
         (s2, s3, s2*s1),
         (s2, s3*s2, s1),
         (s2, s3*s2*s1, 1)]

    We can also access the crystal by specifying a skew shape in terms of `k`-bounded partitions::

        sage: crystals.AffineFactorization([[3,1,1],[1]], 3, k=3)
        Crystal on affine factorizations of type A2 associated to s2*s3*s2*s1

    We can compute the highest weight elements::

        sage: hw = [w for w in B if w.is_highest_weight()]
        sage: hw
        [(1, s2, s3*s2*s1)]
        sage: hw[0].weight()
        (3, 1, 0)

    And show that this crystal is isomorphic to the tableau model of the same weight::

        sage: C = crystals.Tableaux(['A',2],shape=[3,1])
        sage: GC = C.digraph()
        sage: GB = B.digraph()
        sage: GC.is_isomorphic(GB, edge_labels=True)
        True

    The crystal operators themselves move elements between adjacent factors::

        sage: b = hw[0];b
        (1, s2, s3*s2*s1)
        sage: b.f(1)
        (1, s3*s2, s3*s1)

    The cut point `x` is not supposed to occur in the reduced words for `w`::

        sage: B = crystals.AffineFactorization([[3,2],[2]],4,x=0,k=3)
        Traceback (most recent call last):
        ...
        ValueError: x cannot be in reduced word of s0*s3*s2
    """
    @staticmethod
    def __classcall_private__(cls, w, n, x=None, k=None):
        r"""
        Classcall to mend the input.

        TESTS::

            sage: A = crystals.AffineFactorization([[3,1],[1]], 4, k=3); A
            Crystal on affine factorizations of type A3 associated to s3*s2*s1
            sage: AC = crystals.AffineFactorization([Core([4,1],4),Core([1],4)], 4, k=3)
            sage: AC is A
            True
        """
        if k is not None:
            from sage.combinat.core import Core
            from sage.combinat.partition import Partition
            W = WeylGroup(['A',k,1], prefix='s')
            if isinstance(w[0], Core):
                w = [w[0].to_bounded_partition(), w[1].to_bounded_partition()]
            else:
                w = [Partition(w[0]), Partition(w[1])]
            w0 = W.from_reduced_word(w[0].from_kbounded_to_reduced_word(k))
            w1 = W.from_reduced_word(w[1].from_kbounded_to_reduced_word(k))
            w = w0*(w1.inverse())
        return super().__classcall__(cls, w, n, x)

    def __init__(self, w, n, x=None):
        r"""
        EXAMPLES::

            sage: B = crystals.AffineFactorization([[3,2],[2]],4,x=0,k=3)
            Traceback (most recent call last):
            ...
            ValueError: x cannot be in reduced word of s0*s3*s2

            sage: B = crystals.AffineFactorization([[3,2],[2]],4,k=3)
            sage: B.x
            1
            sage: B.w
            s0*s3*s2
            sage: B.k
            3
            sage: B.n
            4

        TESTS::

            sage: W = WeylGroup(['A',3,1], prefix='s')
            sage: w = W.from_reduced_word([2,3,2,1])
            sage: B = crystals.AffineFactorization(w,3)
            sage: TestSuite(B).run()  # long time
        """
        Parent.__init__(self, category=ClassicalCrystals())
        self.n = n
        self.k = w.parent().n-1
        self.w = w
        cartan_type = CartanType(['A',n-1])
        self._cartan_type = cartan_type
        from sage.combinat.sf.sf import SymmetricFunctions
        from sage.rings.rational_field import QQ
        Sym = SymmetricFunctions(QQ)
        s = Sym.schur()
        support = s(w.stanley_symmetric_function()).support()
        support = [ [0]*(n-len(mu))+[mu[len(mu)-i-1] for i in range(len(mu))] for mu in support]
        generators = [tuple(p) for mu in support for p in affine_factorizations(w,n,mu)]
        #generators = [tuple(p) for p in affine_factorizations(w, n)]
        self.module_generators = [self(t) for t in generators]
        if x is None:
            if generators:
                x = min( set(range(self.k+1)).difference(set(
                            sum([i.reduced_word() for i in generators[0]],[]))))
            else:
                x = 0
        if x in set(w.reduced_word()):
            raise ValueError("x cannot be in reduced word of {}".format(w))
        self.x = x

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: W = WeylGroup(['A',3,1], prefix='s')
            sage: w = W.from_reduced_word([3,2,1])
            sage: crystals.AffineFactorization(w,4)
            Crystal on affine factorizations of type A3 associated to s3*s2*s1

            sage: crystals.AffineFactorization([[3,1],[1]], 4, k=3)
            Crystal on affine factorizations of type A3 associated to s3*s2*s1
        """
        return "Crystal on affine factorizations of type A{} associated to {}".format(self.n-1, self.w)

    # temporary workaround while an_element is overridden by Parent
    _an_element_ = EnumeratedSets.ParentMethods._an_element_

    @lazy_attribute
    def _tableaux_isomorphism(self):
        """
        Return the isomorphism from ``self`` to the tableaux model.

        EXAMPLES::

            sage: W = WeylGroup(['A',3,1], prefix='s')
            sage: w = W.from_reduced_word([3,2,1])
            sage: B = crystals.AffineFactorization(w,4)
            sage: B._tableaux_isomorphism
            ['A', 3] Crystal morphism:
              From: Crystal on affine factorizations of type A3 associated to s3*s2*s1
              To:   The crystal of tableaux of type ['A', 3] and shape(s) [[3]]

            sage: W = WeylGroup(['A',3,1], prefix='s')
            sage: w = W.from_reduced_word([2,1,3,2])
            sage: B = crystals.AffineFactorization(w,3)
            sage: B._tableaux_isomorphism
            ['A', 2] Crystal morphism:
              From: Crystal on affine factorizations of type A2 associated to s2*s3*s1*s2
              To:   The crystal of tableaux of type ['A', 2] and shape(s) [[2, 2]]
        """
        # Constructing the tableaux crystal
        from sage.combinat.crystals.tensor_product import CrystalOfTableaux

        def mg_to_shape(mg):
            l = list(mg.weight().to_vector())
            while l and l[-1] == 0:
                l.pop()
            return l
        sh = [mg_to_shape(mg) for mg in self.highest_weight_vectors()]
        C = CrystalOfTableaux(self.cartan_type(), shapes=sh)
        phi = FactorizationToTableaux(Hom(self, C, category=self.category()))
        phi.register_as_coercion()
        return phi

    class Element(ElementWrapper):

        def e(self, i):
            r"""
            Return the action of `e_i` on ``self``.

            EXAMPLES::

                sage: B = crystals.AffineFactorization([[3,1],[1]], 4, k=3)
                sage: W = B.w.parent()
                sage: t = B((W.one(),W.one(),W.from_reduced_word([3]),W.from_reduced_word([2,1]))); t
                (1, 1, s3, s2*s1)
                sage: t.e(1)
                (1, 1, 1, s3*s2*s1)
            """
            if i not in self.index_set():
                raise ValueError("i must be in the index set")
            b = self.bracketing(i)
            if not b[0]:
                return None
            W = self.parent().w.parent()
            x = self.parent().x
            k = self.parent().k
            n = self.parent().n
            a = min(b[0])
            left = [j for j in (self.value[n-i-1]).reduced_word() if j != (a+x) % (k+1)]
            right = [(j-x) % (k+1) for j in (self.value[n-i]).reduced_word()]
            m = max([j for j in range(a) if (j+x) % (k+1) not in left])
            right += [m+1]
            right.sort(reverse=True)
            right = [(j+x) % (k+1) for j in right]
            t = [self.value[j] for j in range(n-i-1)] + [W.from_reduced_word(left)] + [W.from_reduced_word(right)] + [self.value[j] for j in range(n-i+1,n)]
            return self.parent()(tuple(t))

        def f(self, i):
            r"""
            Return the action of `f_i` on ``self``.

            EXAMPLES::

                sage: B = crystals.AffineFactorization([[3,1],[1]], 4, k=3)
                sage: W = B.w.parent()
                sage: t = B((W.one(),W.one(),W.from_reduced_word([3]),W.from_reduced_word([2,1]))); t
                (1, 1, s3, s2*s1)
                sage: t.f(2)
                (1, s3, 1, s2*s1)
                sage: t.f(1)
                (1, 1, s3*s2, s1)
            """
            if i not in self.index_set():
                raise ValueError("i must be in the index set")
            b = self.bracketing(i)
            if not b[1]:
                return None
            W = self.parent().w.parent()
            x = self.parent().x
            k = self.parent().k
            n = self.parent().n
            a = max(b[1])
            right = [j for j in (self.value[n-i]).reduced_word() if j != (a+x) % (k+1)]
            left = [(j-x) % (k+1) for j in (self.value[n-i-1]).reduced_word()]
            m = min([j for j in range(a+1,k+2) if (j+x) % (k+1) not in right])
            left += [m-1]
            left.sort(reverse=True)
            left = [(j+x) % (k+1) for j in left]
            t = [self.value[j] for j in range(n-i-1)] + [W.from_reduced_word(left)] + [W.from_reduced_word(right)] + [self.value[j] for j in range(n-i+1,n)]
            return self.parent()(tuple(t))

        def bracketing(self, i):
            r"""
            Removes all bracketed letters between `i`-th and `i+1`-th entry.

            EXAMPLES::

                sage: B = crystals.AffineFactorization([[3,1],[1]], 3, k=3, x=4)
                sage: W = B.w.parent()
                sage: t = B((W.one(),W.from_reduced_word([3]),W.from_reduced_word([2,1]))); t
                (1, s3, s2*s1)
                sage: t.bracketing(1)
                [[3], [2, 1]]
            """
            n = self.parent().n
            x = self.parent().x
            k = self.parent().k
            right = (self.value[n-i]).reduced_word()
            left = (self.value[n-i-1]).reduced_word()
            right_n = [(j-x) % (k+1) for j in right]
            left_n = [(j-x) % (k+1) for j in left]
            left_unbracketed = []
            while left_n:
                m = max(left_n)
                left_n.remove(m)
                l = [j for j in right_n if j > m]
                if l:
                    right_n.remove(min(l))
                else:
                    left_unbracketed += [m]
            return [list(left_unbracketed), list(right_n)]

        def to_tableau(self):
            """
            Return the tableau representation of ``self``.

            Uses the recording tableau of a minor variation of
            Edelman-Greene insertion. See Theorem 4.11 in [MS2015]_.

            EXAMPLES::

                sage: W = WeylGroup(['A',3,1], prefix='s')
                sage: w = W.from_reduced_word([2,1,3,2])
                sage: B = crystals.AffineFactorization(w,3)
                sage: for x in B:
                ....:     x
                ....:     x.to_tableau().pp()
                (1, s2*s1, s3*s2)
                  1  1
                  2  2
                (s2, s1, s3*s2)
                  1  1
                  2  3
                (s2, s3*s1, s2)
                  1  2
                  2  3
                (s2*s1, 1, s3*s2)
                  1  1
                  3  3
                (s2*s1, s3, s2)
                  1  2
                  3  3
                (s2*s1, s3*s2, 1)
                  2  2
                  3  3
            """
            return self.parent()._tableaux_isomorphism(self)


def affine_factorizations(w, l, weight=None):
    r"""
    Return all factorizations of `w` into `l` factors or of weight ``weight``.

    INPUT:

    - ``w`` -- an (affine) permutation or element of the (affine) Weyl group

    - ``l`` -- nonnegative integer

    - ``weight`` -- (default: ``None``) tuple of nonnegative integers
      specifying the length of the factors

    EXAMPLES::

       sage: W = WeylGroup(['A',3,1], prefix='s')
       sage: w = W.from_reduced_word([3,2,3,1,0,1])
       sage: from sage.combinat.crystals.affine_factorization import affine_factorizations
       sage: affine_factorizations(w,4)
       [[s2, s3, s0, s2*s1*s0],
       [s2, s3, s2*s0, s1*s0],
       [s2, s3, s2*s1*s0, s1],
       [s2, s3*s2, s0, s1*s0],
       [s2, s3*s2, s1*s0, s1],
       [s2, s3*s2*s1, s0, s1],
       [s3*s2, s3, s0, s1*s0],
       [s3*s2, s3, s1*s0, s1],
       [s3*s2, s3*s1, s0, s1],
       [s3*s2*s1, s3, s0, s1]]

       sage: W = WeylGroup(['A',2], prefix='s')
       sage: w0 = W.long_element()
       sage: affine_factorizations(w0,3)
       [[1, s1, s2*s1],
       [1, s2*s1, s2],
       [s1, 1, s2*s1],
       [s1, s2, s1],
       [s1, s2*s1, 1],
       [s2, s1, s2],
       [s2*s1, 1, s2],
       [s2*s1, s2, 1]]
       sage: affine_factorizations(w0,3,(0,1,2))
       [[1, s1, s2*s1]]
       sage: affine_factorizations(w0,3,(1,1,1))
       [[s1, s2, s1], [s2, s1, s2]]
       sage: W = WeylGroup(['A',3], prefix='s')
       sage: w0 = W.long_element()
       sage: affine_factorizations(w0,6,(1,1,1,1,1,1))  # long time
       [[s1, s2, s1, s3, s2, s1],
       [s1, s2, s3, s1, s2, s1],
       [s1, s2, s3, s2, s1, s2],
       [s1, s3, s2, s1, s3, s2],
       [s1, s3, s2, s3, s1, s2],
       [s2, s1, s2, s3, s2, s1],
       [s2, s1, s3, s2, s1, s3],
       [s2, s1, s3, s2, s3, s1],
       [s2, s3, s1, s2, s1, s3],
       [s2, s3, s1, s2, s3, s1],
       [s2, s3, s2, s1, s2, s3],
       [s3, s1, s2, s1, s3, s2],
       [s3, s1, s2, s3, s1, s2],
       [s3, s2, s1, s2, s3, s2],
       [s3, s2, s1, s3, s2, s3],
       [s3, s2, s3, s1, s2, s3]]
       sage: affine_factorizations(w0,6,(0,0,0,1,2,3))
       [[1, 1, 1, s1, s2*s1, s3*s2*s1]]
    """
    if weight is None:
        if l == 0:
            if w.is_one():
                return [[]]
            else:
                return []
        else:
            return [[u]+p for (u,v) in w.left_pieri_factorizations() for p in affine_factorizations(v,l-1) ]
    else:
        if l != len(weight):
            return []
        if l == 0:
            if w.is_one():
                return [[]]
            else:
                return []
        else:
            return [[u]+p for (u,v) in w.left_pieri_factorizations(max_length=weight[0]) if u.length() == weight[0]
                    for p in affine_factorizations(v,l-1,weight[1:]) ]

#####################################################################
## Crystal isomorphisms


class FactorizationToTableaux(CrystalMorphism):
    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        TESTS::

            sage: W = WeylGroup(['A',3,1], prefix='s')
            sage: w = W.from_reduced_word([2,1,3,2])
            sage: B = crystals.AffineFactorization(w,3)
            sage: phi = B._tableaux_isomorphism
            sage: [phi(b) for b in B]
            [[[1, 1], [2, 2]],
             [[1, 1], [2, 3]],
             [[1, 2], [2, 3]],
             [[1, 1], [3, 3]],
             [[1, 2], [3, 3]],
             [[2, 2], [3, 3]]]
        """
        p = []
        q = []
        for i, factor in enumerate(reversed(x.value)):
            word = factor.reduced_word()
            p += [i + 1] * len(word)
            # We sort for those pesky commutative elements
            # The word is most likely in reverse order to begin with
            q += sorted(reversed(word))
        C = self.codomain()
        return C(RSK(p, q, insertion=RSK.rules.EG)[1])

    def is_isomorphism(self) -> bool:
        """
        Return ``True`` as this is an isomorphism.

        EXAMPLES::

            sage: W = WeylGroup(['A',3,1], prefix='s')
            sage: w = W.from_reduced_word([2,1,3,2])
            sage: B = crystals.AffineFactorization(w,3)
            sage: phi = B._tableaux_isomorphism
            sage: phi.is_isomorphism()
            True

        TESTS::

            sage: W = WeylGroup(['A',4,1], prefix='s')
            sage: w = W.from_reduced_word([2,1,3,2,4,3,2,1])
            sage: B = crystals.AffineFactorization(w, 4)         # long time
            sage: phi = B._tableaux_isomorphism                  # long time
            sage: all(phi(b).e(i) == phi(b.e(i)) and             # long time
            ....:     phi(b).f(i) == phi(b.f(i))
            ....:     for b in B for i in B.index_set())
            True
            sage: set(phi(b) for b in B) == set(phi.codomain())  # long time
            True
        """
        return True

    is_embedding = is_isomorphism
    is_surjective = is_isomorphism
