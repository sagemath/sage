r"""
Schur algebras for `GL_n`

This file implements:

- Schur algebras for `GL_n` over an arbitrary field.

- The canonical action of the Schur algebra on a tensor power of the standard
  representation.

- Using the above to calculate the characters of irreducible `GL_n` modules.

AUTHORS:

- Eric Webster (2010-07-01): implement Schur algebra

- Hugh Thomas (2011-05-08): implement action of Schur algebra and characters
  of irreducible modules
"""

# ****************************************************************************
#       Copyright (C) 2010 Eric Webster
#       Copyright (C) 2011 Hugh Thomas <hugh.ross.thomas@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools

from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.rings import Rings
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModule_Tensor
from sage.combinat.integer_lists import IntegerListsLex
from sage.combinat.partition import Partitions, Partition
from sage.combinat.permutation import Permutations
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
from sage.combinat.tableau import SemistandardTableaux
from sage.arith.all import binomial
from sage.matrix.constructor import Matrix
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


def _schur_I_nr_representatives(n, r):
    r"""
    Internal function called by :func:`schur_representation_indices`,
    which generates all weakly increasing tuples of length ``r`` in the
    alphabet ``1, 2, ..., n``.

    EXAMPLES::

        sage: from sage.algebras.schur_algebra import _schur_I_nr_representatives
        sage: _schur_I_nr_representatives(2, 4)
        ((1, 1, 1, 1), (1, 1, 1, 2), (1, 1, 2, 2), (1, 2, 2, 2), (2, 2, 2, 2))
    """
    if r == 0:
        return ()

    index = []
    element = [1]
    while element:
        if element[-1] > n:
            element.pop()
            if element:
                element[-1] += 1
            continue

        if len(element) == r:
            index.append(tuple(element))
            element[-1] += 1
            continue

        element.append(element[-1])

    return tuple(index)


def schur_representative_indices(n, r):
    r"""
    Return a set which functions as a basis of `S_K(n,r)`.

    More specifically, the basis for `S_K(n,r)` consists of
    equivalence classes of pairs of tuples of length ``r`` on the alphabet
    `\{1, \dots, n\}`, where the equivalence relation is simultaneous
    permutation of the two tuples.  We can therefore fix a
    representative for each equivalence class in which the entries of
    the first tuple weakly increase, and the entries of the second tuple
    whose corresponding values in the first tuple are equal, also
    weakly increase.

    EXAMPLES::

        sage: from sage.algebras.schur_algebra import schur_representative_indices
        sage: schur_representative_indices(2, 2)
        [((1, 1), (1, 1)), ((1, 1), (1, 2)),
         ((1, 1), (2, 2)), ((1, 2), (1, 1)),
         ((1, 2), (1, 2)), ((1, 2), (2, 1)),
         ((1, 2), (2, 2)), ((2, 2), (1, 1)),
         ((2, 2), (1, 2)), ((2, 2), (2, 2))]
    """
    basis = []
    I_nr_repr = _schur_I_nr_representatives(n, r)
    for e in I_nr_repr:
        j = 0
        k = 0
        I1 = []
        l = len(e)
        while k < l:
            if e[k] != e[j]:
                I2 = []
                if j == 0:
                    I1 = _schur_I_nr_representatives(n, k)
                else:
                    I2 = _schur_I_nr_representatives(n, k - j)
                    I = []
                    for m1 in range(len(I1)):
                        for m2 in range(len(I2)):
                            I.append(I1[m1] + I2[m2])
                    I1 = I
                j = k
            elif k == l - 1:
                I2 = []
                k += 1
                if j == 0:
                    I1 = _schur_I_nr_representatives(n, k)
                else:
                    I2 = _schur_I_nr_representatives(n, k - j)
                    I = []
                    for m1 in range(len(I1)):
                        for m2 in range(len(I2)):
                            I.append(I1[m1] + I2[m2])
                    I1 = I
            else:
                k += 1

        for v in I1:
            basis.append((tuple(e), tuple(v)))

    return basis


def schur_representative_from_index(i0, i1):
    r"""
    Simultaneously reorder a pair of tuples to obtain the equivalent
    element of the distinguished basis of the Schur algebra.

    .. SEEALSO::

        :func:`schur_representative_indices`

    INPUT:

    - A pair of tuples of length `r` with elements in `\{1,\dots,n\}`

    OUTPUT:

    - The corresponding pair of tuples ordered correctly.

    EXAMPLES::

        sage: from sage.algebras.schur_algebra import schur_representative_from_index
        sage: schur_representative_from_index([2,1,2,2], [1,3,0,0])
        ((1, 2, 2, 2), (3, 0, 0, 1))
    """
    w = []
    for i, val in enumerate(i0):
        w.append((val, i1[i]))
    w.sort()
    i0 = []
    i1 = []
    for pair in w:
        i0.append(pair[0])
        i1.append(pair[1])
    return (tuple(i0), tuple(i1))


class SchurAlgebra(CombinatorialFreeModule):
    r"""
    A Schur algebra.

    Let `R` be a commutative ring, `n` be a positive integer, and `r`
    be a non-negative integer. Define `A_R(n,r)` to be the set of
    homogeneous polynomials of degree `r` in `n^2` variables `x_{ij}`.
    Therefore we can write `R[x_{ij}] = \bigoplus_{r \geq 0} A_R(n,r)`,
    and `R[x_{ij}]` is known to be a bialgebra with coproduct given by
    `\Delta(x_{ij}) = \sum_l x_{il} \otimes x_{lj}` and counit
    `\varepsilon(x_{ij}) = \delta_{ij}`. Therefore `A_R(n,r)` is a
    subcoalgebra of `R[x_{ij}]`. The *Schur algebra* `S_R(n,r)` is the
    linear dual to `A_R(n,r)`, that is `S_R(n,r) := \hom(A_R(n,r), R)`,
    and `S_R(n,r)` obtains its algebra structure naturally by dualizing
    the comultiplication of `A_R(n,r)`.

    Let `V = R^n`. One of the most important properties of the Schur
    algebra `S_R(n, r)` is that it is isomorphic to the endomorphisms
    of `V^{\otimes r}` which commute with the natural action of `S_r`.

    EXAMPLES::

        sage: S = SchurAlgebra(ZZ, 2, 2); S
        Schur algebra (2, 2) over Integer Ring

    REFERENCES:

    - [Gr2007]_
    - :wikipedia:`Schur_algebra`
    """
    def __init__(self, R, n, r):
        """
        Initialize ``self``.

        TESTS::

            sage: S = SchurAlgebra(ZZ, 2, 2)
            sage: TestSuite(S).run()

        ::

            sage: SchurAlgebra(ZZ, -2, 2)
            Traceback (most recent call last):
            ...
            ValueError: n (=-2) must be a positive integer
            sage: SchurAlgebra(ZZ, 2, -2)
            Traceback (most recent call last):
            ...
            ValueError: r (=-2) must be a non-negative integer
            sage: SchurAlgebra('niet', 2, 2)
            Traceback (most recent call last):
            ...
            ValueError: R (=niet) must be a commutative ring
        """
        if n not in ZZ or n <= 0:
            raise ValueError("n (={}) must be a positive integer".format(n))
        if r not in ZZ or r < 0:
            raise ValueError("r (={}) must be a non-negative integer".format(r))
        if R not in Rings.Commutative():
            raise ValueError("R (={}) must be a commutative ring".format(R))

        self._n = n
        self._r = r

        CombinatorialFreeModule.__init__(self, R,
                                         schur_representative_indices(n, r),
                                         prefix='S', bracket=False,
                                         category=AlgebrasWithBasis(R).FiniteDimensional())

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: SchurAlgebra(ZZ, 4, 4)
            Schur algebra (4, 4) over Integer Ring
        """
        msg = "Schur algebra ({}, {}) over {}"
        return msg.format(self._n, self._r, self.base_ring())

    @cached_method
    def one(self):
        """
        Return the element `1` of ``self``.

        EXAMPLES::

            sage: S = SchurAlgebra(ZZ, 2, 2)
            sage: e = S.one(); e
            S((1, 1), (1, 1)) + S((1, 2), (1, 2)) + S((2, 2), (2, 2))

            sage: x = S.an_element()
            sage: x * e == x
            True
            sage: all(e * x == x for x in S.basis())
            True

            sage: S = SchurAlgebra(ZZ, 4, 4)
            sage: e = S.one()
            sage: x = S.an_element()
            sage: x * e == x
            True
        """
        tt = IntegerListsLex(length=self._r, min_part=1, max_part=self._n,
                             min_slope=0)
        words = [tuple(u) for u in tt]
        return self.sum(self._monomial((w, w)) for w in words)

    def product_on_basis(self, e_ij, e_kl):
        r"""
        Return the product of basis elements.

        EXAMPLES::

            sage: S = SchurAlgebra(QQ, 2, 3)
            sage: B = S.basis()

        If we multiply two basis elements `x` and `y`, such that
        `x[1]` and `y[0]` are not permutations of each other, the
        result is zero::

            sage: S.product_on_basis(((1, 1, 1), (1, 1, 2)), ((1, 2, 2), (1, 1, 2)))
            0

        If we multiply a basis element `x` by a basis element which
        consists of the same tuple repeated twice (on either side),
        the result is either zero (if the previous case applies) or `x`::

            sage: ww = B[((1, 2, 2), (1, 2, 2))]
            sage: x = B[((1, 2, 2), (1, 1, 2))]
            sage: ww * x
            S((1, 2, 2), (1, 1, 2))

        An arbitrary product, on the other hand, may have multiplicities::

            sage: x = B[((1, 1, 1), (1, 1, 2))]
            sage: y = B[((1, 1, 2), (1, 2, 2))]
            sage: x * y
            2*S((1, 1, 1), (1, 2, 2))
        """
        j = e_ij[1]

        i = e_ij[0]
        l = e_kl[1]

        l = sorted(l)

        # Find basis elements (p,q) such that p ~ i and q ~ l
        e_pq = []
        for v in self.basis().keys():
            if v[0] == i and sorted(v[1]) == l:
                e_pq.append(v)

        b = self.basis()
        product = self.zero()

        # Find s in I(n,r) such that (p,s) ~ (i,j) and (s,q) ~ (k,l)
        for e in e_pq:
            Z_ijklpq = self.base_ring().zero()
            for s in Permutations([xx for xx in j]):
                if (schur_representative_from_index(e[0], s) == e_ij
                        and schur_representative_from_index(s, e[1]) == e_kl):
                    Z_ijklpq += self.base_ring().one()
            product += Z_ijklpq * b[e]

        return product

    def dimension(self):
        r"""
        Return the dimension of ``self``.

        The dimension of the Schur algebra `S_R(n, r)` is

        .. MATH::

            \dim S_R(n,r) = \binom{n^2+r-1}{r}.

        EXAMPLES::

            sage: S = SchurAlgebra(QQ, 4, 2)
            sage: S.dimension()
            136
            sage: S = SchurAlgebra(QQ, 2, 4)
            sage: S.dimension()
            35
        """
        return binomial(self._n ** 2 + self._r - 1, self._r)


class SchurTensorModule(CombinatorialFreeModule_Tensor):
    r"""
    The space `V^{\otimes r}` where `V = R^n` equipped with a left action
    of the Schur algebra `S_R(n,r)` and a right action of the symmetric
    group `S_r`.

    Let `R` be a commutative ring and `V = R^n`. We consider the module
    `V^{\otimes r}` equipped with a natural right action of the symmetric
    group `S_r` given by

    .. MATH::

        (v_1 \otimes v_2 \otimes \cdots \otimes v_n) \sigma
        = v_{\sigma(1)} \otimes v_{\sigma(2)} \otimes \cdots
        \otimes v_{\sigma(n)}.

    The Schur algebra `S_R(n,r)` is naturally isomorphic to the
    endomorphisms of `V^{\otimes r}` which commutes with the `S_r` action.
    We get the natural left action of `S_R(n,r)` by this isomorphism.

    EXAMPLES::

        sage: T = SchurTensorModule(QQ, 2, 3); T
        The 3-fold tensor product of a free module of dimension 2
         over Rational Field
        sage: A = SchurAlgebra(QQ, 2, 3)
        sage: P = Permutations(3)
        sage: t = T.an_element(); t
        2*B[1] # B[1] # B[1] + 2*B[1] # B[1] # B[2] + 3*B[1] # B[2] # B[1]
        sage: a = A.an_element(); a
        2*S((1, 1, 1), (1, 1, 1)) + 2*S((1, 1, 1), (1, 1, 2))
         + 3*S((1, 1, 1), (1, 2, 2))
        sage: p = P.an_element(); p
        [3, 1, 2]
        sage: y = a * t; y
        14*B[1] # B[1] # B[1]
        sage: y * p
        14*B[1] # B[1] # B[1]
        sage: z = t * p; z
        2*B[1] # B[1] # B[1] + 3*B[1] # B[1] # B[2] + 2*B[2] # B[1] # B[1]
        sage: a * z
        14*B[1] # B[1] # B[1]

    We check the commuting action property::

        sage: all( (bA * bT) * p == bA * (bT * p)
        ....:      for bT in T.basis() for bA in A.basis() for p in P)
        True
    """
    def __init__(self, R, n, r):
        """
        Initialize ``self``.

        TESTS::

            sage: T = SchurTensorModule(QQ, 2, 3)
            sage: TestSuite(T).run()
        """
        C = CombinatorialFreeModule(R, list(range(1, n + 1)))
        self._n = n
        self._r = r
        self._sga = SymmetricGroupAlgebra(R, r)
        self._schur = SchurAlgebra(R, n, r)
        cat = ModulesWithBasis(R).TensorProducts().FiniteDimensional()
        CombinatorialFreeModule_Tensor.__init__(self, tuple([C] * r), category=cat)
        g = self._schur.module_morphism(self._monomial_product, codomain=self)
        self._schur_action = self.module_morphism(g, codomain=self, position=1)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: SchurTensorModule(QQ, 2, 3)
            The 3-fold tensor product of a free module of dimension 2
            over Rational Field
        """
        msg = "The {}-fold tensor product of a free module of dimension {}"
        msg += " over {}"
        return msg.format(self._r, self._n, self.base_ring())

    def _monomial_product(self, xi, v):
        """
        Result of acting by the basis element ``xi`` of the corresponding
        Schur algebra on the basis element ``v`` of ``self``.

        EXAMPLES::

            sage: T = SchurTensorModule(QQ, 2, 3)
            sage: xi = T._schur.basis().keys()[4]; xi
            ((1, 1, 2), (1, 1, 1))
            sage: T._monomial_product(xi, (1, 1, 1))
            B[1] # B[1] # B[2] + B[1] # B[2] # B[1] + B[2] # B[1] # B[1]
        """
        ret = []
        for i in itertools.product(list(range(1, self._n + 1)), repeat=self._r):
            if schur_representative_from_index(i, v) == xi:
                ret.append(tuple(i))
        return self.sum_of_monomials(ret)

    class Element(CombinatorialFreeModule_Tensor.Element):
        def _acted_upon_(self, elt, self_on_left=False):
            """
            Return the action of ``elt`` on ``self``.

            We add the *left* action of the Schur algebra, and the *right*
            actions of the symmetric group algebra and the symmetric group.

            EXAMPLES::

                sage: T = SchurTensorModule(QQ, 2, 4)
                sage: x = T.an_element()
                sage: A = SchurAlgebra(QQ, 2, 4)
                sage: y = A.an_element()
                sage: y * x
                14*B[1] # B[1] # B[1] # B[1]
                sage: x * y
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: ...

            ::

                sage: SGA = SymmetricGroupAlgebra(QQ, 4)
                sage: y = SGA.an_element()
                sage: x * y
                14*B[1] # B[1] # B[1] # B[1] + 17*B[1] # B[1] # B[1] # B[2]
                 + 7*B[1] # B[1] # B[2] # B[1] + 9*B[1] # B[2] # B[1] # B[1]
                 + 2*B[2] # B[1] # B[1] # B[1]
                sage: y * x
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: ...

            ::

                sage: S = Permutations(4)
                sage: y = S.an_element()
                sage: x * y
                2*B[1] # B[1] # B[1] # B[1] + 3*B[1] # B[1] # B[1] # B[2]
                 + 2*B[2] # B[1] # B[1] # B[1]
            """
            P = self.parent()
            if self_on_left:
                if elt in P._sga:
                    return P.sum_of_terms((tuple([m[i - 1] for i in me]),
                                           c * ce)
                                          for m, c in self for me, ce in elt)

                if elt in P._sga._indices:
                    return P.sum_of_terms((tuple([m[i - 1] for i in elt]), c)
                                          for m, c in self)

            elif elt in P._schur:  # self_on_left is False
                return P._schur_action(elt, self)
            return super(SchurTensorModule.Element, self)._acted_upon_(elt, self_on_left)


def GL_irreducible_character(n, mu, KK):
    r"""
    Return the character of the irreducible module indexed by ``mu``
    of `GL(n)` over the field ``KK``.

    INPUT:

    - ``n`` -- a positive integer
    - ``mu`` -- a partition of at most ``n`` parts
    - ``KK`` -- a field

    OUTPUT:

    a symmetric function which should be interpreted in ``n``
    variables to be meaningful as a character

    EXAMPLES:

    Over `\QQ`, the irreducible character for `\mu` is the Schur
    function associated to `\mu`, plus garbage terms (Schur
    functions associated to partitions with more than `n` parts)::

        sage: from sage.algebras.schur_algebra import GL_irreducible_character
        sage: sbasis = SymmetricFunctions(QQ).s()
        sage: z = GL_irreducible_character(2, [2], QQ)
        sage: sbasis(z)
        s[2]

        sage: z = GL_irreducible_character(4, [3, 2], QQ)
        sage: sbasis(z)
        -5*s[1, 1, 1, 1, 1] + s[3, 2]

    Over a Galois field, the irreducible character for `\mu` will
    in general be smaller.

    In characteristic `p`, for a one-part partition `(r)`, where
    `r = a_0 + p a_1 + p^2 a_2 + \dots`, the result is (see [Gr2007]_,
    after 5.5d) the product of `h[a_0], h[a_1]( pbasis[p]), h[a_2]
    ( pbasis[p^2]), \dots,` which is consistent with the following ::

        sage: from sage.algebras.schur_algebra import GL_irreducible_character
        sage: GL_irreducible_character(2, [7], GF(3))
        m[4, 3] + m[6, 1] + m[7]
    """
    mbasis = SymmetricFunctions(QQ).m()
    r = sum(mu)
    M = SchurTensorModule(KK, n, r)
    A = M._schur
    SGA = M._sga

    #make ST the superstandard tableau of shape mu
    from sage.combinat.tableau import from_shape_and_word
    ST = from_shape_and_word(mu, list(range(1, r + 1)), convention='English')

    #make ell the reading word of the highest weight tableau of shape mu
    ell = [i + 1 for i, l in enumerate(mu) for dummy in range(l)]

    e = M.basis()[tuple(ell)]  # the element e_l

    # This is the notation `\{X\}` from just before (5.3a) of [Gr2007]_.
    S = SGA._indices
    BracC = SGA._from_dict({S(x.tuple()): x.sign() for x in ST.column_stabilizer()},
                           remove_zeros=False)
    f = e * BracC  # M.action_by_symmetric_group_algebra(e, BracC)

    # [Green, Theorem 5.3b] says that a basis of the Carter-Lusztig
    # module V_\mu is given by taking this f, and multiplying by all
    # xi_{i,ell} with ell as above and i semistandard.

    carter_lusztig = []
    for T in SemistandardTableaux(mu, max_entry=n):
        i = tuple(flatten(T))
        schur_rep = schur_representative_from_index(i, tuple(ell))
        y = A.basis()[schur_rep] * e  # M.action_by_Schur_alg(A.basis()[schur_rep], e)
        carter_lusztig.append(y.to_vector())

    #Therefore, we now have carter_lusztig as a list giving the basis
    #of `V_\mu`

    #We want to think of expressing this character as a sum of monomial
    #symmetric functions.

    #We will determine a basis element for each m_\lambda in the
    #character, and we want to keep track of them by \lambda.

    #That means that we only want to pick out the basis elements above for
    #those semistandard words whose content is a partition.

    contents = Partitions(r, max_length=n).list()
    # all partitions of r, length at most n

    # JJ will consist of a list for each element of `contents`,
    # recording the list
    # of semistandard tableaux words with that content

    # graded_basis will consist of the corresponding basis element
    graded_basis = []
    JJ = []
    for i in range(len(contents)):
        graded_basis.append([])
        JJ.append([])
    for T in SemistandardTableaux(mu, max_entry=n):
        i = tuple(flatten(T))
        # Get the content of T
        con = [0] * n
        for a in i:
            con[a - 1] += 1
        try:
            P = Partition(con)
            P_index = contents.index(P)
            JJ[P_index].append(i)
            schur_rep = schur_representative_from_index(i, tuple(ell))
            x = A.basis()[schur_rep] * f  # M.action_by_Schur_alg(A.basis()[schur_rep], f)
            graded_basis[P_index].append(x.to_vector())
        except ValueError:
            pass

    #There is an inner product on the Carter-Lusztig module V_\mu; its
    #maximal submodule is exactly the kernel of the inner product.

    #Now, for each possible partition content, we look at the graded piece of
    #that degree, and we record how these elements pair with each of the
    #elements of carter_lusztig.

    #The kernel of this pairing is the part of this graded piece which is
    #not in the irreducible module for \mu.

    length = len(carter_lusztig)

    phi = mbasis.zero()
    for aa in range(len(contents)):
        mat = []
        for kk in range(len(JJ[aa])):
            temp = []
            for j in range(length):
                temp.append(graded_basis[aa][kk].inner_product(carter_lusztig[j]))
            mat.append(temp)
        angle = Matrix(mat)
        phi += (len(JJ[aa]) - angle.nullity()) * mbasis(contents[aa])
    return phi
