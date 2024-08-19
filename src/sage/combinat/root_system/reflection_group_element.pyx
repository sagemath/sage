r"""
Reflection group elements

.. SEEALSO::

    :mod:`sage.combinat.root_system.reflection_group_complex`,
    :mod:`sage.combinat.root_system.reflection_group_real`

AUTHORS:

- Christian Stump (initial version 2011--2015)
- Travis Scrimshaw (14-03-2017): moved element code
"""
# ***************************************************************************
#       Copyright (C) 2011-2016 Christian Stump <christian.stump at gmail.com>
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************
import re

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc_c import prod
from sage.arith.functions import lcm
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.combinat.root_system.reflection_group_c import reduced_word_c, reduce_in_coset
from sage.matrix.constructor import Matrix
from sage.matrix.special import identity_matrix

TUPLE = re.compile(r'(?:\([0-9,]*\))+')


cdef class ComplexReflectionGroupElement(PermutationGroupElement):
    """
    An element in a complex reflection group.
    """
    def __hash__(self):
        r"""
        Return a hash for this reflection group element.

        This hash stores both the element as a reduced word and the parent group.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',5])                      # optional - gap3
            sage: W_hash = set(hash(w) for w in W)                  # optional - gap3
            sage: len(W_hash) == W.cardinality()                    # optional - gap3
            True

        TESTS:

        Check that types B and C are hashed differently, see :issue:`29726`::

            sage: WB = ReflectionGroup(['B',5])                     # optional - gap3
            sage: WC = ReflectionGroup(['C',5])                     # optional - gap3

            sage: WB_hash = set(hash(w) for w in WB)                # optional - gap3
            sage: WC_hash = set(hash(w) for w in WC)                # optional - gap3

            sage: len(WB_hash) == WB.cardinality()                  # optional - gap3
            True

            sage: len(WC_hash) == WC.cardinality()                  # optional - gap3
            True

            sage: WB_hash.intersection(WC_hash)                     # optional - gap3
            set()

        Check that :issue:`34912` is fixed::

            sage: # optional - gap3
            sage: G4 = ReflectionGroup(4)
            sage: g0, g1 = G4.gens()
            sage: elt = g0^2 * g1 * g0^2 * g1
            sage: elt
            (1,12)(2,24)(3,19)(4,22)(5,17)(6,20)(7,23)(8,9)(10,21)(11,13)(14,18)(15,16)
            sage: y = (elt * G4.gen(1)) * G4.gen(1) * G4.gen(1)
            sage: elt == y
            True
            sage: hash(elt) == hash(y)
            True
        """
        return hash(self._parent) | super().__hash__()

    def reduced_word(self):
        r"""
        Return a word in the simple reflections to obtain ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((5,1,1), index_set=['a'], hyperplane_index_set=['x'], reflection_index_set=['A','B','C','D']) # optional - gap3
            sage: [w.reduced_word() for w in W]                     # optional - gap3
            [[], ['a'], ['a', 'a'], ['a', 'a', 'a'], ['a', 'a', 'a', 'a']]

        .. SEEALSO:: :meth:`reduced_word_in_reflections`
        """
        I = self._parent._index_set
        return [I[i] for i in self._reduced_word]

    @lazy_attribute
    def _reduced_word(self):
        r"""
        Compute a reduced word and stores it into ``self._reduced_word``.

        TESTS::

            sage: W = ReflectionGroup((5,1,1))                      # optional - gap3
            sage: w = W.an_element()                                # optional - gap3
            sage: w._reduced_word                                   # optional - gap3
            [0]
        """
        W = self._parent
        gens = [W.simple_reflection(j) for j in W._index_set]
        return _gap_factorization(self, gens)

    # @cached_in_parent_method
    def reduced_word_in_reflections(self):
        r"""
        Return a word in the reflections to obtain ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup((5,1,1), index_set=['a'], reflection_index_set=['A','B','C','D']) # optional - gap3
            sage: [w.reduced_word_in_reflections() for w in W]      # optional - gap3
            [[], ['A'], ['B'], ['C'], ['D']]

        .. SEEALSO:: :meth:`reduced_word`
        """
        if self.is_one():
            return []

        W = self._parent
        gens = [W.reflection(j) for j in W._reflection_index_set]
        word = _gap_factorization(self, gens)
        return [self._parent._reflection_index_set[i] for i in word]

    def length(self):
        r"""
        Return the length of ``self`` in generating reflections.

        This is the minimal numbers of generating reflections needed
        to obtain ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(4)                            # optional - gap3
            sage: for w in W.iteration_tracking_words():            # optional - gap3
            ....:     print("{} {}".format(w.reduced_word(), w.length()))
            [] 0
            [1] 1
            [2] 1
            [1, 1] 2
            [1, 2] 2
            [2, 1] 2
            [2, 2] 2
            [1, 1, 2] 3
            [1, 2, 1] 3
            [1, 2, 2] 3
            [2, 1, 1] 3
            [2, 2, 1] 3
            [1, 1, 2, 1] 4
            [1, 1, 2, 2] 4
            [1, 2, 1, 1] 4
            [1, 2, 2, 1] 4
            [2, 1, 1, 2] 4
            [2, 2, 1, 1] 4
            [1, 1, 2, 1, 1] 5
            [1, 1, 2, 2, 1] 5
            [1, 2, 1, 1, 2] 5
            [1, 2, 2, 1, 1] 5
            [1, 1, 2, 1, 1, 2] 6
            [1, 1, 2, 2, 1, 1] 6

            sage: data = {w: (len(w.reduced_word()), w.length())    # optional - gap3
            ....:         for w in W.iteration_tracking_words()}
            sage: for w in W:                                       # optional - gap3
            ....:     assert data[w] == (w.length(), w.length()), w
        """
        return ZZ(len(self.reduced_word()))

    # @cached_in_parent_method
    def to_matrix(self, on_space='primal'):
        r"""
        Return ``self`` as a matrix acting on the underlying vector
        space.

        - ``on_space`` -- (default: ``'primal'``) whether
          to act as the reflection representation on the given
          basis, or to act on the dual reflection representation
          on the dual basis

        EXAMPLES::

            sage: W = ReflectionGroup((3,1,2))                                          # optional - gap3
            sage: data = {w: [w.to_matrix(), w.to_matrix(on_space='dual')] for w in W}  # optional - gap3
            sage: for w in W.iteration_tracking_words():                                # optional - gap3
            ....:     w.reduced_word()
            ....:     mats = [w.to_matrix(), w.to_matrix(on_space='dual')]
            ....:     mats
            ....:     assert data[w] == mats
            []
            [
            [1 0]  [1 0]
            [0 1], [0 1]
            ]
            [1]
            [
            [E(3)    0]  [E(3)^2      0]
            [   0    1], [     0      1]
            ]
            [2]
            [
            [0 1]  [0 1]
            [1 0], [1 0]
            ]
            [1, 1]
            [
            [E(3)^2      0]  [E(3)    0]
            [     0      1], [   0    1]
            ]
            [1, 2]
            [
            [   0 E(3)]  [     0 E(3)^2]
            [   1    0], [     1      0]
            ]
            [2, 1]
            [
            [   0    1]  [     0      1]
            [E(3)    0], [E(3)^2      0]
            ]
            [1, 1, 2]
            [
            [     0 E(3)^2]  [   0 E(3)]
            [     1      0], [   1    0]
            ]
            [1, 2, 1]
            [
            [   0 E(3)]  [     0 E(3)^2]
            [E(3)    0], [E(3)^2      0]
            ]
            [2, 1, 1]
            [
            [     0      1]  [   0    1]
            [E(3)^2      0], [E(3)    0]
            ]
            [2, 1, 2]
            [
            [   1    0]  [     1      0]
            [   0 E(3)], [     0 E(3)^2]
            ]
            [1, 1, 2, 1]
            [
            [     0 E(3)^2]  [     0   E(3)]
            [  E(3)      0], [E(3)^2      0]
            ]
            [1, 2, 1, 1]
            [
            [     0   E(3)]  [     0 E(3)^2]
            [E(3)^2      0], [  E(3)      0]
            ]
            [1, 2, 1, 2]
            [
            [E(3)    0]  [E(3)^2      0]
            [   0 E(3)], [     0 E(3)^2]
            ]
            [2, 1, 1, 2]
            [
            [     1      0]  [   1    0]
            [     0 E(3)^2], [   0 E(3)]
            ]
            [1, 1, 2, 1, 1]
            [
            [     0 E(3)^2]  [   0 E(3)]
            [E(3)^2      0], [E(3)    0]
            ]
            [1, 1, 2, 1, 2]
            [
            [E(3)^2      0]  [  E(3)      0]
            [     0   E(3)], [     0 E(3)^2]
            ]
            [1, 2, 1, 1, 2]
            [
            [  E(3)      0]  [E(3)^2      0]
            [     0 E(3)^2], [     0   E(3)]
            ]
            [1, 1, 2, 1, 1, 2]
            [
            [E(3)^2      0]  [E(3)    0]
            [     0 E(3)^2], [   0 E(3)]
            ]
        """
        W = self._parent
        if W._reflection_representation is None:
            mat = self.canonical_matrix()
        else:
            refl_repr = W._reflection_representation
            id_mat = identity_matrix(QQ, refl_repr[W.index_set()[0]].nrows())
            mat = prod((refl_repr[i] for i in self.reduced_word()), id_mat)

        if on_space == "primal":
            pass
        elif on_space == "dual":
            mat = mat.inverse().transpose()
        else:
            raise ValueError('on_space must be "primal" or "dual"')

        mat.set_immutable()
        return mat

    matrix = to_matrix

    def canonical_matrix(self):
        r"""
        Return the matrix of ``self`` in the canonical faithful representation.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: W = WeylGroup(['A',2], prefix='s', implementation='permutation')
            sage: for w in W:
            ....:     w.reduced_word()
            ....:     w.canonical_matrix()
            []
            [1 0]
            [0 1]
            [2]
            [ 1  1]
            [ 0 -1]
            [1]
            [-1  0]
            [ 1  1]
            [1, 2]
            [-1 -1]
            [ 1  0]
            [2, 1]
            [ 0  1]
            [-1 -1]
            [1, 2, 1]
            [ 0 -1]
            [-1  0]
        """
        W = self._parent
        Phi = W.roots()
        cdef list inds = [W._index_set_inverse[i] for i in W.independent_roots().keys()]
        mat = W.base_change_matrix() * Matrix([Phi[self.perm[i]] for i in inds])
        mat.set_immutable()
        return mat

    cpdef action(self, vec, on_space='primal'):
        r"""
        Return the image of ``vec`` under the action of ``self``.

        INPUT:

        - ``vec`` -- vector in the basis given by the simple root

        - ``on_space`` -- (default: ``'primal'``) whether
          to act as the reflection representation on the given
          basis, or to act on the dual reflection representation
          on the dual basis

        EXAMPLES::

            sage: W = ReflectionGroup((3,1,2))                      # optional - gap3
            sage: w = W.from_reduced_word([1, 2, 1, 1, 2])          # optional - gap3
            sage: for alpha in W.independent_roots():               # optional - gap3
            ....:     print("%s -> %s"%(alpha,w.action(alpha)))
            (1, 0) -> (E(3), 0)
            (-1, 1) -> (-E(3), E(3)^2)
        """
        mat = self.matrix(on_space=on_space)
        return vec * mat

    cpdef _act_on_(self, vec, bint self_on_left):
        r"""
        Define the action of ``self`` as a linear transformation
        on the vector space, in the basis given by the simple
        roots.

        - ``vec`` -- the vector (an iterable) to act on

        - ``self_on_left`` -- whether the action of ``self`` is on
          the left or on the right

        EXAMPLES::

            sage: W = ReflectionGroup((3,1,2))                      # optional - gap3
            sage: w = W.from_reduced_word([1, 2, 1, 1, 2])          # optional - gap3
            sage: for alpha in W.independent_roots():               # optional - gap3
            ....:     print("%s -> %s"%(alpha, w * alpha))
            (1, 0) -> (E(3), 0)
            (-1, 1) -> (-E(3), E(3)^2)
        """
        if not self_on_left:
            return (~self).action(vec)
        return self.action(vec)

    cpdef action_on_root_indices(self, i):
        """
        Return the right action on the set of roots.

        INPUT:

        - ``i`` -- index of the root to act on

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',3])
            sage: w = W.w0
            sage: N = len(W.roots())
            sage: [w.action_on_root_indices(i) for i in range(N)]
            [8, 7, 6, 10, 9, 11, 2, 1, 0, 4, 3, 5]

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2], reflection_index_set=['A','B','C'])
            sage: w = W.w0
            sage: N = len(W.roots())
            sage: [w.action_on_root_indices(i) for i in range(N)]
            [4, 3, 5, 1, 0, 2]

        TESTS::

            sage: W = ReflectionGroup(4)                 # optional - gap3
            sage: N = len(W.roots())                     # optional - gap3
            sage: all(sorted([w.action_on_root_indices(i) for i in range(N)]) == list(range(N)) for w in W)   # optional - gap3
            True
        """
        return self.perm[i]

    def action_on_root(self, root):
        r"""
        Return the root obtained by applying ``self`` to ``root``
        on the right.

        INPUT:

        - ``root`` -- the root to act on

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                      # optional - gap3
            sage: for w in W:                                       # optional - gap3
            ....:     print("%s %s"%(w.reduced_word(),
            ....:           [w.action_on_root(beta,side='left') for beta in W.positive_roots()]))
            [] [(1, 0), (0, 1), (1, 1)]
            [2] [(1, 1), (0, -1), (1, 0)]
            [1] [(-1, 0), (1, 1), (0, 1)]
            [1, 2] [(0, 1), (-1, -1), (-1, 0)]
            [2, 1] [(-1, -1), (1, 0), (0, -1)]
            [1, 2, 1] [(0, -1), (-1, 0), (-1, -1)]

        TESTS::

            sage: W = ReflectionGroup(4); Phi = sorted(W.roots())   # optional - gap3
            sage: all(sorted([w.action_on_root(beta) for beta in Phi]) == Phi for w in W)   # optional - gap3
            True
        """
        Phi = self._parent.roots()
        return Phi[self.action_on_root_indices(Phi.index(root))]

    def to_permutation_of_roots(self):
        r"""
        Return ``self`` as a permutation of the roots with indexing
        starting at `1`.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                      # optional - gap3
            sage: for w in W:       # optional - gap3
            ....:     perm = w.to_permutation_of_roots()
            ....:     print("{} {}".format(perm, perm == w))
            () True
            (1,3)(2,5)(4,6) True
            (1,4)(2,3)(5,6) True
            (1,6,2)(3,5,4) True
            (1,2,6)(3,4,5) True
            (1,5)(2,4)(3,6) True
        """
        W = self._parent
        return PermutationGroupElement(self, W)

    # @cached_in_parent_method
    def fix_space(self):
        r"""
        Return the fix space of ``self``.

        This is the sub vector space of the underlying vector space
        on which ``self`` acts trivially.

        EXAMPLES::

            sage: W = ReflectionGroup((1,1,3))                      # optional - gap3
            sage: for w in W:                                       # optional - gap3
            ....:     w.reduced_word()
            ....:     w.fix_space()
            []
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
            [2]
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 0]
            [1]
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [0 1]
            [1, 2]
            Vector space of degree 2 and dimension 0 over Rational Field
            Basis matrix:
            []
            [2, 1]
            Vector space of degree 2 and dimension 0 over Rational Field
            Basis matrix:
            []
            [1, 2, 1]
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [ 1 -1]

            sage: W = ReflectionGroup(23)                 # optional - gap3
            sage: W.gen(0).fix_space()                    # optional - gap3
            Vector space of degree 3 and dimension 2 over Universal Cyclotomic Field
            Basis matrix:
            [0 1 0]
            [0 0 1]
        """
        I = identity_matrix(QQ, self._parent.rank())
        return (self.to_matrix() - I).right_kernel()

    # @cached_in_parent_method
    def reflection_eigenvalues(self, is_class_representative=False):
        r"""
        Return the reflection eigenvalues of ``self``.

        INPUT:

        - ``is_class_representative`` -- boolean (default: ``False``); whether
          to first replace ``self`` by the representative of its
          conjugacy class

        EXAMPLES::

            sage: W = ReflectionGroup(4)                            # optional - gap3
            sage: for w in W.iteration_tracking_words():            # optional - gap3
            ....:      w.reflection_eigenvalues()
            [0, 0]
            [1/3, 0]
            [1/3, 0]
            [2/3, 0]
            [1/6, 1/2]
            [1/6, 1/2]
            [2/3, 0]
            [1/4, 3/4]
            [1/4, 3/4]
            [1/4, 3/4]
            [1/4, 3/4]
            [1/4, 3/4]
            [1/3, 0]
            [1/2, 5/6]
            [1/3, 0]
            [1/2, 5/6]
            [1/2, 5/6]
            [1/2, 5/6]
            [1/6, 1/2]
            [2/3, 0]
            [1/6, 1/2]
            [2/3, 0]
            [1/2, 1/2]
            [1/4, 3/4]

            sage: data = {w: w.reflection_eigenvalues() for w in W}  # optional - gap3
            sage: all(w.reflection_eigenvalues() == data[w] for w in W.iteration_tracking_words())  # optional - gap3
            True
        """
        return self._parent.reflection_eigenvalues(self, is_class_representative=is_class_representative)

    # @cached_in_parent_method
    def galois_conjugates(self):
        r"""
        Return all Galois conjugates of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(4)                            # optional - gap3
            sage: for w in W.iteration_tracking_words():            # optional - gap3
            ....:     print(w.galois_conjugates())
            [[1 0]
             [0 1]]
            [[   1    0]
             [   0 E(3)], [     1      0]
             [     0 E(3)^2]]
            [[ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
             [ 4/3*E(3) + 2/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2],
             [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
             [ 2/3*E(3) + 4/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]]
            [[     1      0]
             [     0 E(3)^2], [   1    0]
             [   0 E(3)]]
            [[ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
             [-2/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2],
             [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
             [ 2/3*E(3) - 2/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]]
            [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
             [ 4/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2],
             [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
             [ 2/3*E(3) + 4/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]]
            [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
             [ 2/3*E(3) + 4/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2],
             [ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
             [ 4/3*E(3) + 2/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]]
            [[ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
             [-2/3*E(3) - 4/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2],
             [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
             [-4/3*E(3) - 2/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]]
            [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
             [-2/3*E(3) + 2/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2],
             [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
             [ 2/3*E(3) - 2/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]]
            [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
             [-4/3*E(3) - 2/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2],
             [ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
             [-2/3*E(3) - 4/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]]
            [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
             [ 4/3*E(3) + 2/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2],
             [-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
             [ 2/3*E(3) + 4/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]]
            [[-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
             [ 2/3*E(3) + 4/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2],
             [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
             [ 4/3*E(3) + 2/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]]
            [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
             [-2/3*E(3) - 4/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2],
             [-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
             [-4/3*E(3) - 2/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]]
            [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]
             [ 2/3*E(3) - 2/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2],
             [ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
             [-2/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]]
            [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
             [-2/3*E(3) + 2/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2],
             [-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
             [ 2/3*E(3) - 2/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]]
            [[-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
             [-4/3*E(3) - 2/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2],
             [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
             [-2/3*E(3) - 4/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]]
            [[   -1     0]
             [    0 -E(3)], [     -1       0]
             [      0 -E(3)^2]]
            [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
             [ 2/3*E(3) + 4/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2],
             [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
             [ 4/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]]
            [[ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
             [-2/3*E(3) - 4/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2],
             [-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
             [-4/3*E(3) - 2/3*E(3)^2  1/3*E(3) + 2/3*E(3)^2]]
            [[-1/3*E(3) + 1/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2]
             [ 2/3*E(3) - 2/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2],
             [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]
             [-2/3*E(3) + 2/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]]
            [[     -1       0]
             [      0 -E(3)^2], [   -1     0]
             [    0 -E(3)]]
            [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
             [-4/3*E(3) - 2/3*E(3)^2 -2/3*E(3) - 1/3*E(3)^2],
             [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
             [-2/3*E(3) - 4/3*E(3)^2 -1/3*E(3) - 2/3*E(3)^2]]
            [[-1  0]
             [ 0 -1]]
            [[-1/3*E(3) + 1/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2]
             [ 2/3*E(3) - 2/3*E(3)^2  1/3*E(3) - 1/3*E(3)^2],
             [ 1/3*E(3) - 1/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]
             [-2/3*E(3) + 2/3*E(3)^2 -1/3*E(3) + 1/3*E(3)^2]]

            sage: data = {w: w.galois_conjugates() for w in W}      # optional - gap3
            sage: all(w.galois_conjugates() == data[w] for w in W.iteration_tracking_words())  # optional - gap3
            True
        """
        rk = self._parent.rank()
        M = self.to_matrix().list()
        m = lcm([x.conductor() if hasattr(x, "conductor") else 1 for x in M])
        cdef list M_gals = [x.galois_conjugates(m) if hasattr(x, "galois_conjugates") else [x] for x in M]
        cdef list conjugates = []
        cdef int i
        for i in range(len(M_gals[0])):
            conjugates.append(Matrix(rk, [X[i] for X in M_gals]))
        return conjugates

cdef class RealReflectionGroupElement(ComplexReflectionGroupElement):
    @lazy_attribute
    def _reduced_word(self):
        r"""
        Compute a reduced word and stores it into ``self._reduced_word``.
        The words are in ``range(n)`` and not in the index set.

        TESTS::

            sage: W = ReflectionGroup(['A',2])                      # optional - gap3
            sage: [w._reduced_word for w in W]                      # optional - gap3
            [[], [1], [0], [0, 1], [1, 0], [0, 1, 0]]
        """
        return reduced_word_c(self._parent, self)

    def reduced_word_in_reflections(self):
        r"""
        Return a word in the reflections to obtain ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2], index_set=['a','b'], reflection_index_set=['A','B','C']) # optional - gap3
            sage: [(w.reduced_word(), w.reduced_word_in_reflections()) for w in W]  # optional - gap3
            [([], []),
             (['b'], ['B']),
             (['a'], ['A']),
             (['a', 'b'], ['A', 'B']),
             (['b', 'a'], ['A', 'C']),
             (['a', 'b', 'a'], ['C'])]

        .. SEEALSO:: :meth:`reduced_word`
        """
        if self.is_one():
            return []

        W = self._parent
        r = self.reflection_length()
        R = W.reflections()
        I = W.reflection_index_set()
        cdef list word = []
        cdef RealReflectionGroupElement w
        while r > 0:
            for i in I:
                w = <RealReflectionGroupElement>(R[i]._mul_(self))
                if w.reflection_length() < r:
                    word.append(i)
                    r -= 1
                    self = w
                    break
        return word

    def length(self):
        r"""
        Return the length of ``self`` in generating reflections.

        This is the minimal numbers of generating reflections needed
        to obtain ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                      # optional - gap3
            sage: for w in W:                                       # optional - gap3
            ....:     print("%s %s"%(w.reduced_word(), w.length()))
            [] 0
            [2] 1
            [1] 1
            [1, 2] 2
            [2, 1] 2
            [1, 2, 1] 3
        """
        return ZZ(len(self._reduced_word))

    cpdef bint has_left_descent(self, i) noexcept:
        r"""
        Return whether ``i`` is a left descent of ``self``.

        This is done by testing whether ``i`` is mapped by ``self``
        to a negative root.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(["A",3])
            sage: s = W.simple_reflections()
            sage: (s[1]*s[2]).has_left_descent(1)
            True
            sage: (s[1]*s[2]).has_left_descent(2)
            False
        """
        W = self._parent
        # we also check == because 0-based indexing
        return self.perm[W._index_set_inverse[i]] >= W.number_of_reflections()

    cpdef bint has_descent(self, i, side='left', positive=False) noexcept:
        r"""
        Return whether `i` is a descent (or ascent) of ``self``.

        This is done by testing whether `i` is mapped by ``self``
        to a negative root.

        INPUT:

        - ``i`` -- an index of a simple reflection
        - ``side`` -- (default: ``'right'``) ``'left'`` or ``'right'``
        - ``positive`` -- boolean (default: ``False``)

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(["A",3])
            sage: s = W.simple_reflections()
            sage: (s[1]*s[2]).has_descent(1)
            True
            sage: (s[1]*s[2]).has_descent(2)
            False
        """
        if not isinstance(positive, bool):
            raise TypeError("%s is not a boolean" % (bool))

        if i not in self._parent.index_set():
            raise ValueError("the given index %s is not in the index set" % i)

        negative = not positive

        if side == 'left':
            return self.has_left_descent(i) is negative
        elif side == 'right':
            return self.has_right_descent(i) is negative
        else:
            raise ValueError('side must be "left" or "right"')

    def coset_representative(self, index_set, side='right'):
        """
        Return the unique shortest element of the Coxeter group
        `W` which is in the same left (resp. right) coset as
        ``self``, with respect to the parabolic subgroup `W_I`.

        INPUT:

        - ``index_set`` -- a subset (or iterable) of the nodes of the index set
        - ``side`` -- (default: ``right``) ``'left'`` or ``'right'``

        EXAMPLES::

            sage: # needs sage.graphs
            sage: W = CoxeterGroup(['A',4], implementation='permutation')
            sage: s = W.simple_reflections()
            sage: w = s[2] * s[1] * s[3]
            sage: w.coset_representative([]).reduced_word()
            [2, 1, 3]
            sage: w.coset_representative([1]).reduced_word()
            [2, 3]
            sage: w.coset_representative([1,2]).reduced_word()
            [2, 3]
            sage: w.coset_representative([1,3]                 ).reduced_word()
            [2]
            sage: w.coset_representative([2,3]                 ).reduced_word()
            [2, 1]
            sage: w.coset_representative([1,2,3]               ).reduced_word()
            []
            sage: w.coset_representative([],      side = 'left').reduced_word()
            [2, 1, 3]
            sage: w.coset_representative([1],     side = 'left').reduced_word()
            [2, 1, 3]
            sage: w.coset_representative([1,2],   side = 'left').reduced_word()
            [3]
            sage: w.coset_representative([1,3],   side = 'left').reduced_word()
            [2, 1, 3]
            sage: w.coset_representative([2,3],   side = 'left').reduced_word()
            [1]
            sage: w.coset_representative([1,2,3], side = 'left').reduced_word()
            []
        """
        S = tuple(self._parent.simple_reflections())
        N = self._parent.number_of_reflections()
        I = self._parent._index_set_inverse
        return reduce_in_coset(self, S, [I[i] for i in index_set], N, side=='left')

    def to_matrix(self, side='right', on_space='primal'):
        r"""
        Return ``self`` as a matrix acting on the underlying vector
        space.

        - ``side`` -- (default: ``'right'``) whether the
          action of ``self`` is on the ``'left'`` or on the ``'right'``

        - ``on_space`` -- (default: ``'primal'``) whether
          to act as the reflection representation on the given
          basis, or to act on the dual reflection representation
          on the dual basis

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])           # optional - gap3
            sage: for w in W:                            # optional - gap3
            ....:     w.reduced_word()
            ....:     [w.to_matrix(), w.to_matrix(on_space='dual')]
            []
            [
            [1 0]  [1 0]
            [0 1], [0 1]
            ]
            [2]
            [
            [ 1  1]  [ 1  0]
            [ 0 -1], [ 1 -1]
            ]
            [1]
            [
            [-1  0]  [-1  1]
            [ 1  1], [ 0  1]
            ]
            [1, 2]
            [
            [-1 -1]  [ 0 -1]
            [ 1  0], [ 1 -1]
            ]
            [2, 1]
            [
            [ 0  1]  [-1  1]
            [-1 -1], [-1  0]
            ]
            [1, 2, 1]
            [
            [ 0 -1]  [ 0 -1]
            [-1  0], [-1  0]
            ]

        TESTS::

            sage: W = ReflectionGroup(['F',4])           # optional - gap3
            sage: all(w.to_matrix(side='left') == W.from_reduced_word(reversed(w.reduced_word())).to_matrix(side='right').transpose() for w in W) # optional - gap3
            True
            sage: all(w.to_matrix(side='right') == W.from_reduced_word(reversed(w.reduced_word())).to_matrix(side='left').transpose() for w in W) # optional - gap3
            True
        """
        W = self._parent
        cdef RealReflectionGroupElement w
        if W._reflection_representation is None:
            if side == "left":
                w = <RealReflectionGroupElement>(~self)
            elif side == "right":
                w = <RealReflectionGroupElement>(self)
            else:
                raise ValueError('side must be "left" or "right"')
            mat = w.canonical_matrix()
        else:
            refl_repr = W._reflection_representation
            id_mat = identity_matrix(QQ, refl_repr[W.index_set()[0]].nrows())
            mat = prod((refl_repr[i] for i in self.reduced_word()), id_mat)

        if on_space == "primal":
            if side == "left":
                mat = mat.transpose()
        elif on_space == "dual":
            if side == "left":
                mat = mat.inverse()
            else:
                mat = mat.inverse().transpose()
        else:
            raise ValueError('on_space must be "primal" or "dual"')

        mat.set_immutable()
        return mat

    matrix = to_matrix

    cpdef action(self, vec, side='right', on_space='primal'):
        r"""
        Return the image of ``vec`` under the action of ``self``.

        INPUT:

        - ``vec`` -- vector in the basis given by the simple root

        - ``side`` -- (default: ``'right'``) whether the
          action of ``self`` is on the ``'left'`` or on the ``'right'``

        - ``on_space`` -- (default: ``'primal'``) whether
          to act as the reflection representation on the given
          basis, or to act on the dual reflection representation
          on the dual basis

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                      # optional - gap3
            sage: for w in W:                                       # optional - gap3
            ....:     print("%s %s"%(w.reduced_word(),
            ....:           [w.action(weight,side='left') for weight in W.fundamental_weights()]))
            [] [(2/3, 1/3), (1/3, 2/3)]
            [2] [(2/3, 1/3), (1/3, -1/3)]
            [1] [(-1/3, 1/3), (1/3, 2/3)]
            [1, 2] [(-1/3, 1/3), (-2/3, -1/3)]
            [2, 1] [(-1/3, -2/3), (1/3, -1/3)]
            [1, 2, 1] [(-1/3, -2/3), (-2/3, -1/3)]

        TESTS::

            sage: W = ReflectionGroup(['B',3])                      # optional - gap3
            sage: all(w.action(alpha,side='right') == w.action_on_root(alpha,side='right')  # optional - gap3
            ....:     for w in W for alpha in W.simple_roots())
            True
            sage: all(w.action(alpha,side='left') == w.action_on_root(alpha,side='left')  #optional - gap3
            ....:     for w in W for alpha in W.simple_roots())
            True
        """
        W = self._parent
        n = W.rank()
        Phi = W.roots()
        cdef RealReflectionGroupElement w
        if side == "right":
            w = <RealReflectionGroupElement> self
        elif side == "left":
            w = <RealReflectionGroupElement>(~self)
        else:
            raise ValueError('side must be "left" or "right"')
        cdef int j
        ret = Phi[0].parent().zero()
        if on_space == "primal":
            for j in range(n):
                ret += vec[j] * Phi[w.perm[j]]
            return ret
        elif on_space == "dual":
            w = <RealReflectionGroupElement>(~w)
            for j in range(n):
                ret += Phi[w.perm[j]] * vec[j]
            return ret
        else:
            raise ValueError('on_space must be "primal" or "dual"')

    cpdef _act_on_(self, vec, bint self_on_left):
        r"""
        Give the action of ``self`` as a linear transformation on
        the vector space, in the basis given by the simple roots.

        INPUT:

        - ``vec`` -- the vector (an iterable) to act on

        - ``self_on_left`` -- whether the action of ``self`` is on
          the left or on the right

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])           # optional - gap3
            sage: w = W.from_reduced_word([1,2])         # optional - gap3
            sage: for root in W.positive_roots():        # optional - gap3
            ....:     print("%s -> %s"%(root, w*root))
            (1, 0) -> (0, 1)
            (0, 1) -> (-1, -1)
            (1, 1) -> (-1, 0)

            sage: for root in W.positive_roots():        # optional - gap3
            ....:     print("%s -> %s"%(root, root*w))
            (1, 0) -> (-1, -1)
            (0, 1) -> (1, 0)
            (1, 1) -> (0, -1)
        """
        if self_on_left:
            return self.action(vec, side='left')
        return self.action(vec, side='right')

    cpdef action_on_root_indices(self, i, side='right'):
        """
        Return the action on the set of roots.

        INPUT:

        - ``i`` -- index of the root to act on

        - ``side`` -- (default: ``'right'``) whether the
          action is on the left or on the right

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',3])
            sage: w = W.w0
            sage: N = len(W.roots())
            sage: [w.action_on_root_indices(i,side='left') for i in range(N)]
            [8, 7, 6, 10, 9, 11, 2, 1, 0, 4, 3, 5]

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2], reflection_index_set=['A','B','C'])
            sage: w = W.w0
            sage: N = len(W.roots())
            sage: [w.action_on_root_indices(i,side='left') for i in range(N)]
            [4, 3, 5, 1, 0, 2]
        """
        cdef RealReflectionGroupElement w
        if side == "right":
            w = self
        elif side == "left":
            w = <RealReflectionGroupElement>(~self)
        else:
            raise ValueError('side must be "left" or "right"')
        return w.perm[i]

    def action_on_root(self, root, side='right'):
        r"""
        Return the root obtained by applying ``self`` to ``root``.

        INPUT:

        - ``root`` -- the root to act on

        - ``side`` -- (default: ``'right'``) whether the
          action is on the left or on the right

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])           # optional - gap3
            sage: for w in W:                            # optional - gap3
            ....:     print("%s %s"%(w.reduced_word(),
            ....:           [w.action_on_root(beta,side='left') for beta in W.positive_roots()]))
            [] [(1, 0), (0, 1), (1, 1)]
            [2] [(1, 1), (0, -1), (1, 0)]
            [1] [(-1, 0), (1, 1), (0, 1)]
            [1, 2] [(0, 1), (-1, -1), (-1, 0)]
            [2, 1] [(-1, -1), (1, 0), (0, -1)]
            [1, 2, 1] [(0, -1), (-1, 0), (-1, -1)]

            sage: W = ReflectionGroup(['A',2])           # optional - gap3
            sage: for w in W:                            # optional - gap3
            ....:     print("%s %s"%(w.reduced_word(),
            ....:           [w.action_on_root(beta,side='right') for beta in W.positive_roots()]))
            [] [(1, 0), (0, 1), (1, 1)]
            [2] [(1, 1), (0, -1), (1, 0)]
            [1] [(-1, 0), (1, 1), (0, 1)]
            [1, 2] [(-1, -1), (1, 0), (0, -1)]
            [2, 1] [(0, 1), (-1, -1), (-1, 0)]
            [1, 2, 1] [(0, -1), (-1, 0), (-1, -1)]
        """
        Phi = self._parent.roots()
        return Phi[self.action_on_root_indices(Phi.index(root), side=side)]

    def inversion_set(self, side='right'):
        r"""
        Return the inversion set of ``self``.

        This is the set `\{\beta \in \Phi^+ : s(\beta) \in \Phi^-\}`,
        where `s` is ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                      # optional - gap3
            sage: for w in W:                                       # optional - gap3
            ....:     print("%s %s"%(w.reduced_word(), w.inversion_set()))
            [] []
            [2] [(0, 1)]
            [1] [(1, 0)]
            [1, 2] [(1, 0), (1, 1)]
            [2, 1] [(0, 1), (1, 1)]
            [1, 2, 1] [(1, 0), (0, 1), (1, 1)]

            sage: W.from_reduced_word([1,2]).inversion_set(side='left') # optional - gap3
            [(0, 1), (1, 1)]
        """
        N = self._parent.number_of_reflections()
        Phi = self._parent.roots()
        cdef int i
        if side == "left":
            self = <RealReflectionGroupElement>(~self)
        elif side != "right":
            raise ValueError('side must be "left" or "right"')
        return [Phi[i] for i in range(N) if self.perm[i] >= N]


def _gap_factorization(w, gens):
    r"""
    Return a factorization of ``w`` using the generators ``gens``.

    .. WARNING::

        This is only available through GAP3 and Chevie.

    EXAMPLES::

        sage: from sage.combinat.root_system.reflection_group_element import _gap_factorization
        sage: W = ReflectionGroup((1,1,3))                              # optional - gap3
        sage: gens = [W.simple_reflection(i) for i in W.index_set()]    # optional - gap3
        sage: [_gap_factorization(w,gens) for w in W]                   # optional - gap3
        [[], [1], [0], [0, 1], [1, 0], [0, 1, 0]]
    """
    from sage.interfaces.gap3 import gap3

    gap3.execute('W := GroupWithGenerators(%s)' % str(gens))
    gap3.execute(_gap_factorization_code)
    fac = gap3('MinimalWord(W,%s)' % str(w)).sage()
    return [i-1 for i in fac]


_gap_factorization_code = r"""
# MinimalWord(G,w)
# given a permutation group G find some expression of minimal length in the
# generators of G and their inverses of the element w (an inverse is
# represented by a negative index).
# To speed up  later calls to  the same function  the fields G.base, G.words,
# G.nbwordslength are kept.
MinimalWord:=function(G,w)
  local decode,i,p,g,h,n,bag,nbe,nbf,new,gens,inds;
  # to save space elements of G are represented as image of the base, and
  # words are represented as: index of previous elt, last generator applied;
  if not IsBound(G.base) then
    StabChain(G);g:=G; G.base:=[];
    while IsBound(g.orbit) do Add(G.base,g.orbit[1]); g:=g.stabilizer; od;
  fi;
  w:=OnTuples(G.base,w);
  if not IsBound(G.words) then
    G.words:=[G.base]; G.lastmult:=[[0,0]];
    G.nbwordslength:=[1];
  fi;
  gens:=ShallowCopy(G.generators);inds:=[1..Length(gens)];
  #  for g in G.generators do
  #    if g<>g^-1 then Add(gens,g^-1);Add(inds,-Position(gens,g));fi;
  #  od;
  bag:=Set(G.words);
  nbe:=0;nbf:=0;
  decode:=function(i)local w;w:=[];
    while i<>1 do Add(w,G.lastmult[i][2]); i:=G.lastmult[i][1];od;
    return Reversed(w);
  end;
  while true do
    if w in bag then return decode(Position(G.words,w));fi;
    new:=Length(G.words);
    for g in [1..Length(gens)] do
      for h in [1+Sum(G.nbwordslength{[1..Length(G.nbwordslength)-1]})..new] do
         n:=OnTuples(G.words[h],gens[g]);
         if n in bag then
           nbe:=nbe+1;# if nbe mod 500=1 then Print(".\c");fi;
         else
           nbf:=nbf+1;# if nbf mod 500=1 then Print("*\c");fi;
       Add(G.words,n);Add(G.lastmult,[h,inds[g]]);AddSet(bag,n);
         fi;
       od;
    od;
    Add(G.nbwordslength,Length(G.words)-new);
    Print("\n",G.nbwordslength[Length(G.nbwordslength)]," elements of length ",
      Length(G.nbwordslength)-1);
  od;
end;"""


def _gap_return(S, coerce_obj='self'):
    r"""
    Return the string ``S`` after a few modifications are done.

    This is a stupid internal function to take GAP output as a string,
    replace a few things, to then turn it into a Sage object.

    TESTS::

        sage: from sage.combinat.root_system.reflection_group_complex import _gap_return
        sage: _gap_return("[ (), (1,4)(2,3)(5,6), (1,6,2)(3,5,4) ]")    # optional - gap3
        ['()', '(1,4)(2,3)(5,6)', '(1,6,2)(3,5,4)']
    """
    S = S.replace(' ', '').replace('\n', '')
    return TUPLE.findall(S)
