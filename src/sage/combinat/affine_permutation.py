# sage.doctest: needs sage.combinat sage.groups
r"""
Affine Permutations
"""

# ****************************************************************************
#       Copyright (C) 2013 Tom Denton <sdenton4@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from itertools import repeat

from sage.categories.affine_weyl_groups import AffineWeylGroups
from sage.combinat.composition import Composition
from sage.combinat.partition import Partition
from sage.combinat.root_system.cartan_type import CartanType
from sage.misc.cachefunc import cached_method
from sage.misc.constant_function import ConstantFunction
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.misc.prandom import randint
from sage.rings.integer import Integer
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

lazy_import('sage.combinat.root_system.weyl_group', 'WeylGroup')
lazy_import('sage.groups.perm_gps.permgroup_named', 'SymmetricGroup')


class AffinePermutation(ClonableArray):
    r"""
    An affine permutation, represented in the window notation, and
    considered as a bijection from `\ZZ` to `\ZZ`.

    EXAMPLES::

        sage: A = AffinePermutationGroup(['A',7,1])
        sage: p = A([3, -1, 0, 6, 5, 4, 10, 9])
        sage: p
        Type A affine permutation with window [3, -1, 0, 6, 5, 4, 10, 9]
    """

    def __init__(self, parent, lst, check=True):
        r"""
        Initialize ``self``.

        INPUT:

        - ``parent`` -- the parent affine permutation group

        - ``lst`` -- list giving the base window of the affine permutation

        - ``check`` -- whether to test if the affine permutation is valid

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p = A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p
            Type A affine permutation with window [3, -1, 0, 6, 5, 4, 10, 9]
            sage: TestSuite(p).run()

        TESTS:

        Check that :issue:`26436` is fixed::

            sage: A = AffinePermutationGroup(['A',3,1])
            sage: p = A([-3/1,2/1,3/1,8/1])
            sage: q = ~p
            sage: q * p
            Type A affine permutation with window [1, 2, 3, 4]
        """
        if check:
            lst = [Integer(val) for val in lst]
        self.k = parent.k
        self.n = self.k + 1
        # This N doesn't matter for type A, but comes up in all other types.
        if parent.cartan_type()[0] == 'A':
            self.N = self.n
        elif parent.cartan_type()[0] in ['B', 'C', 'D']:
            self.N = 2 * self.k + 1
        elif parent.cartan_type()[0] == 'G':
            self.N = 6
        else:
            raise NotImplementedError('unsupported Cartan type')
        ClonableArray.__init__(self, parent, lst, check)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p = A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p
            Type A affine permutation with window [3, -1, 0, 6, 5, 4, 10, 9]
        """
        return ("Type " + self.parent().cartan_type().letter
                + " affine permutation with window " + str(list(self)))

    def __rmul__(self, q):
        r"""
        Given ``self`` and `q`, returns ``self*q``.

        INPUT:

        - ``q`` -- an element of ``self.parent()``

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p = A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: q = A([0, 2, 3, 4, 5, 6, 7, 9])
            sage: p.__rmul__(q)
            Type A affine permutation with window [1, -1, 0, 6, 5, 4, 10, 11]
        """
        l = [self.value(q.value(i)) for i in range(1,len(self)+1)]
        return type(self)(self.parent(), l, check=False)

    def __lmul__(self, q):
        r"""
        Given ``self`` and `q`, returns ``q*self``.

        INPUT:

        - ``q`` -- an element of ``self.parent()``

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p = A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: q = A([0,2,3,4,5,6,7,9])
            sage: p.__lmul__(q)
            Type A affine permutation with window [3, -1, 1, 6, 5, 4, 10, 8]
        """
        #if self.parent().right_to_left:
        #    self,q=q,self
        #... product rule
        l = [q.value(self.value(i)) for i in range(1,len(self)+1)]
        return type(self)(self.parent(), l, check=False)

    def __mul__(self, q):
        r"""
        Given ``self`` and `q`, returns ``self*q``.

        INPUT:

        - ``q`` -- an element of ``self.parent()``

        EXAMPLES::

            sage: p = AffinePermutationGroup(['A',7,1])([3, -1, 0, 6, 5, 4, 10, 9])
            sage: s1 = AffinePermutationGroup(['A',7,1]).one().apply_simple_reflection(1)
            sage: p * s1
            Type A affine permutation with window [-1, 3, 0, 6, 5, 4, 10, 9]
            sage: p.apply_simple_reflection(1, 'right')
            Type A affine permutation with window [-1, 3, 0, 6, 5, 4, 10, 9]
        """
        return self.__rmul__(q)

    @cached_method
    def __invert__(self):
        r"""
        Return the inverse affine permutation.

        EXAMPLES::

            sage: p = AffinePermutationGroup(['A',7,1])([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.inverse()  # indirect doctest
            Type A affine permutation with window [0, -1, 1, 6, 5, 4, 10, 11]
        """
        inv = [self.position(i) for i in range(1, len(self) + 1)]
        return type(self)(self.parent(), inv, check=False)

    def apply_simple_reflection(self, i, side='right'):
        r"""
        Apply a simple reflection.

        INPUT:

        - ``i`` -- integer
        - ``side`` -- (default: ``'right'``) determines whether to apply the
          reflection on the ``'right'`` or ``'left'``

        EXAMPLES::

            sage: p = AffinePermutationGroup(['A',7,1])([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.apply_simple_reflection(3)
            Type A affine permutation with window [3, -1, 6, 0, 5, 4, 10, 9]
            sage: p.apply_simple_reflection(11)
            Type A affine permutation with window [3, -1, 6, 0, 5, 4, 10, 9]
            sage: p.apply_simple_reflection(3, 'left')
            Type A affine permutation with window [4, -1, 0, 6, 5, 3, 10, 9]
            sage: p.apply_simple_reflection(11, 'left')
            Type A affine permutation with window [4, -1, 0, 6, 5, 3, 10, 9]
        """
        if side == 'right':
            return self.apply_simple_reflection_right(i)
        if side == 'left':
            return self.apply_simple_reflection_left(i)

    def __call__(self, i):
        r"""
        Return the image of the integer ``i`` under this permutation.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.value(1) #indirect doctest
            3
            sage: p.value(9)
            11
        """
        return self.value(i)

    def is_i_grassmannian(self, i=0, side='right') -> bool:
        r"""
        Test whether ``self`` is `i`-grassmannian, i.e., either is the
        identity or has ``i`` as the sole descent.

        INPUT:

        - ``i`` -- an element of the index set
        - ``side`` -- determines the side on which to check the descents

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.is_i_grassmannian()
            False
            sage: q=A.from_word([3,2,1,0])
            sage: q.is_i_grassmannian()
            True
            sage: q=A.from_word([2,3,4,5])
            sage: q.is_i_grassmannian(5)
            True
            sage: q.is_i_grassmannian(2, side='left')
            True
        """
        return self == self.parent().one() or self.descents(side) == [i]

    def index_set(self) -> tuple:
        r"""
        Index set of the affine permutation group.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: A.index_set()
            (0, 1, 2, 3, 4, 5, 6, 7)
        """
        return tuple(range(self.k + 1))

    def lower_covers(self, side='right'):
        r"""
        Return lower covers of ``self``.

        The set of affine permutations of one less length related by
        multiplication by a simple transposition on the indicated side.
        These are the elements that ``self`` covers in weak order.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.lower_covers()
            [Type A affine permutation with window [-1, 3, 0, 6, 5, 4, 10, 9],
             Type A affine permutation with window [3, -1, 0, 5, 6, 4, 10, 9],
             Type A affine permutation with window [3, -1, 0, 6, 4, 5, 10, 9],
             Type A affine permutation with window [3, -1, 0, 6, 5, 4, 9, 10]]
        """
        S = self.descents(side)
        return [self.apply_simple_reflection(i, side) for i in S]

    def is_one(self) -> bool:
        r"""
        Test whether the affine permutation is the identity.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.is_one()
            False
            sage: q=A.one()
            sage: q.is_one()
            True
        """
        return self == self.parent().one()

    def reduced_word(self):
        r"""
        Return a reduced word for the affine permutation.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.reduced_word()
            [0, 7, 4, 1, 0, 7, 5, 4, 2, 1]
        """
        #This is about 25% faster than the default algorithm.
        x = self
        i = 0
        word = []
        while not x.is_one():
            if x.has_descent(i):
                x = x.apply_simple_reflection_right(i)
                word.append(i)
            i = (i+1) % (self.k+1)
        word.reverse()
        return word

    def signature(self):
        r"""
        Signature of the affine permutation, `(-1)^l`, where `l` is the
        length of the permutation.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.signature()
            1
        """
        return (-1)**self.length()

    @cached_method
    def to_weyl_group_element(self):
        r"""
        The affine Weyl group element corresponding to the affine permutation.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.to_weyl_group_element()
            [ 0 -1  0  1  0  0  1  0]
            [ 1 -1  0  1  0  0  1 -1]
            [ 1 -1  0  1  0  0  0  0]
            [ 0  0  0  1  0  0  0  0]
            [ 0  0  0  1  0 -1  1  0]
            [ 0  0  0  1 -1  0  1  0]
            [ 0  0  0  0  0  0  1  0]
            [ 0 -1  1  0  0  0  1  0]
        """
        W = self.parent().weyl_group()
        return W.from_reduced_word(self.reduced_word())

    def grassmannian_quotient(self, i=0, side='right'):
        r"""
        Return the Grassmannian quotient.

        Factors ``self`` into a unique product of a Grassmannian and a
        finite-type element.  Returns a tuple containing the Grassmannian
        and finite elements, in order according to ``side``.

        INPUT:

        - ``i`` -- (default: 0) an element of the index set; the descent
          checked for

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: gq=p.grassmannian_quotient()
            sage: gq
            (Type A affine permutation with window [-1, 0, 3, 4, 5, 6, 9, 10],
             Type A affine permutation with window [3, 1, 2, 6, 5, 4, 8, 7])
            sage: gq[0].is_i_grassmannian()
            True
            sage: 0 not in gq[1].reduced_word()
            True
            sage: prod(gq)==p
            True

            sage: gqLeft=p.grassmannian_quotient(side='left')
            sage: 0 not in gqLeft[0].reduced_word()
            True
            sage: gqLeft[1].is_i_grassmannian(side='left')
            True
            sage: prod(gqLeft)==p
            True
        """
        fin = self.parent().one()
        gr = self
        D = gr.descents(side=side)
        while not (D == [i] or D == []):
            m = D[0]
            if m == i:
                m = D[1]
            if side == 'right':
                fin = fin.apply_simple_reflection(m, side='left')
                gr = gr.apply_simple_reflection(m, side='right')
            else:
                fin = fin.apply_simple_reflection(m, side='right')
                gr = gr.apply_simple_reflection(m, side='left')
            D = gr.descents(side=side)
        if side == 'right':
            return (gr, fin)
        else:
            return (fin, gr)


class AffinePermutationTypeA(AffinePermutation):
    #----------------------
    #Type-specific methods.
    #(Methods existing in all types, but with type-specific definition.)
    #----------------------
    def check(self):
        r"""
        Check that ``self`` is an affine permutation.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p = A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p
            Type A affine permutation with window [3, -1, 0, 6, 5, 4, 10, 9]
            sage: q = A([1,2,3])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: length of list must be k+1=8
            sage: q = A([1,2,3,4,5,6,7,0])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: window does not sum to 36
            sage: q = A([1,1,3,4,5,6,7,9])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: entries must have distinct residues
        """
        if not self:
            return
        k = self.parent().k
        # Type A
        if len(self) != k + 1:
            raise ValueError(f"length of list must be k+1={k + 1}")
        sigma = (k + 2).binomial(2)
        if sigma != sum(self):
            raise ValueError(f"window does not sum to {sigma}")
        l = sorted(i % (k + 1) for i in self)
        if any(i != j for i, j in enumerate(l)):
            raise ValueError("entries must have distinct residues")

    def value(self, i, base_window=False):
        r"""
        Return the image of the integer ``i`` under this permutation.

        INPUT:

        - ``base_window`` -- boolean; indicating whether ``i`` is in the
          base window; if ``True``, will run a bit faster, but the method
          will screw up if ``i`` is not actually in the index set

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p = A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.value(1)
            3
            sage: p.value(9)
            11
        """
        if base_window:
            self[i-1]
        window = (i-1) // (self.k+1)
        return self[(i-1) % (self.k+1)] + window*(self.k+1)

    def position(self, i):
        r"""
        Find the position ``j`` such the ``self.value(j) == i``.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p = A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.position(3)
            1
            sage: p.position(11)
            9
        """
        for r in range(self.k+1):
            if self[r] % (self.k+1) == i % (self.k+1):
                #i sits in position i, but some number of windows away.
                diff = (i-self[r]) // (self.k+1)
                return r + diff*(self.k+1) + 1
        return False

    def apply_simple_reflection_right(self, i):
        r"""
        Apply the simple reflection to positions `i`, `i+1`.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: p = AffinePermutationGroup(['A',7,1])([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.apply_simple_reflection_right(3)
            Type A affine permutation with window [3, -1, 6, 0, 5, 4, 10, 9]
            sage: p.apply_simple_reflection_right(11)
            Type A affine permutation with window [3, -1, 6, 0, 5, 4, 10, 9]
        """
        j = i % (self.k+1)
        #Cloning is currently kinda broken, in that caches don't clear which
        #leads to strangeness with the cloned object.
        #The clone approach is quite a bit (2x) faster, though, so this should
        #switch once the caching situation is fixed.
        #with self.clone(check=False) as l:
        l = self[:]
        if j == 0:
            a = l[0]
            l[0] = l[-1] - (self.k+1)
            l[-1] = a + (self.k+1)
        else:
            a = l[j-1]
            l[j-1] = l[j]
            l[j] = a
        #return l
        return type(self)(self.parent(), l, check=False)

    def apply_simple_reflection_left(self, i):
        r"""
        Apply the simple reflection to the values `i`, `i+1`.

        EXAMPLES::

            sage: p = AffinePermutationGroup(['A',7,1])([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.apply_simple_reflection_left(3)
            Type A affine permutation with window [4, -1, 0, 6, 5, 3, 10, 9]
            sage: p.apply_simple_reflection_left(11)
            Type A affine permutation with window [4, -1, 0, 6, 5, 3, 10, 9]
        """
        #Here are a couple other methods we tried out, but turned out
        #to be slower than the current implementation.
        #1) This one was very bad:
        #   return self.inverse().apply_simple_reflection_right(i).inverse()
        #2) Also bad, though not quite so bad:
        #   return (self.parent().simple_reflection(i))*self
        i = i % (self.k+1)
        #Cloning is currently kinda broken, in that caches don't clear which
        #leads to strangeness with the cloned object.
        #The clone approach is quite a bit faster, though, so this should switch
        #once the caching situation is fixed.
        #with self.clone(check=False) as l:
        l = []
        if i != self.k:
            for m in range(self.k + 1):
                res = self[m] % (self.k + 1)
                if res == i:
                    l.append(self[m] + 1)
                elif res == i + 1:
                    l.append(self[m] - 1)
                else:
                    l.append(self[m])
        if i == self.k:
            for m in range(self.k + 1):
                res = self[m] % (self.k + 1)
                if res == i:
                    l.append(self[m] + 1)
                elif res == 0:
                    l.append(self[m] - 1)
                else:
                    l.append(self[m])
        return type(self)(self.parent(), l, check=False)

    def has_right_descent(self, i) -> bool:
        r"""
        Determine whether there is a descent at ``i``.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: p = AffinePermutationGroup(['A',7,1])([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.has_right_descent(1)
            True
            sage: p.has_right_descent(9)
            True
            sage: p.has_right_descent(0)
            False
        """
        return self.value(i) > self.value(i+1)

    def has_left_descent(self, i) -> bool:
        r"""
        Determine whether there is a descent at ``i``.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: p = AffinePermutationGroup(['A',7,1])([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.has_left_descent(1)
            True
            sage: p.has_left_descent(9)
            True
            sage: p.has_left_descent(0)
            True
        """
        # This is much faster than the default method of taking the inverse and
        # then finding right descents...
        return self.position(i) > self.position(i + 1)

    def to_type_a(self):
        r"""
        Return an embedding of ``self`` into the affine permutation group of
        type `A`.  (For type `A`, just returns ``self``.)

        EXAMPLES::

            sage: p = AffinePermutationGroup(['A',7,1])([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.to_type_a() is p
            True
        """
        return self

    #----------------------
    #Type-A-specific methods.
    #Only available in Type A.
    #----------------------

    def flip_automorphism(self):
        r"""
        The Dynkin diagram automorphism which fixes `s_0` and reverses all
        other indices.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.flip_automorphism()
            Type A affine permutation with window [0, -1, 5, 4, 3, 9, 10, 6]
        """
        #Note: There should be a more combinatorial (ie, faster) way to do this.
        w = [(self.k+1-i) % (self.k+1) for i in self.reduced_word()]
        return self.parent().from_word(w)

    def promotion(self):
        r"""
        The Dynkin diagram automorphism which sends `s_i` to `s_{i+1}`.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p = A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.promotion()
            Type A affine permutation with window [2, 4, 0, 1, 7, 6, 5, 11]
        """
        l = [self[-1] - self.k]
        l.extend(self[i] + 1 for i in range(self.k))
        return type(self)(self.parent(), l)

    def maximal_cyclic_factor(self, typ='decreasing', side='right', verbose=False):
        r"""
        For an affine permutation `x`, find the unique maximal subset `A`
        of the index set such that `x = yd_A` is a reduced product.

        INPUT:

        - ``typ`` -- ``'increasing'`` or ``'decreasing'``
          (default: ``'decreasing'``); chooses whether to find increasing
          or decreasing sets

        - ``side`` -- ``'right'`` or ``'left'`` (default: ``'right'``); chooses
          whether to find maximal sets starting from the left or the right

        - ``verbose`` -- boolean;  if ``True``, outputs information about how
          the cyclically increasing element was found

        EXAMPLES::

            sage: p = AffinePermutationGroup(['A',7,1])([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.maximal_cyclic_factor()
            [7, 5, 4, 2, 1]
            sage: p.maximal_cyclic_factor(side='left')
            [1, 0, 7, 5, 4]
            sage: p.maximal_cyclic_factor('increasing','right')
            [4, 5, 7, 0, 1]
            sage: p.maximal_cyclic_factor('increasing','left')
            [0, 1, 2, 4, 5]
        """
        k = self.k
        if side[0] == 'r':
            descents = self.descents(side='right')
            side = 'right'
        else:
            descents = self.descents(side='left')
            side = 'left'
        #for now, assume side is 'right')
        best_T = []
        for i in descents:
            y = self.clone().apply_simple_reflection(i,side)
            T = [i]
            j = i
            for _ in range(1, self.k):
                if (typ[0],side[0]) == ('d', 'r'):
                    j = (j+1) % (k+1)
                if (typ[0],side[0]) == ('i', 'r'):
                    j = (j-1) % (k+1)
                if (typ[0],side[0]) == ('d', 'l'):
                    j = (j-1) % (k+1)
                if (typ[0],side[0]) == ('i', 'l'):
                    j = (j+1) % (k+1)
                if y.has_descent(j, side):
                    y = y.apply_simple_reflection(j,side)
                    T.append(j % (k+1))
            if verbose:
                print(i, T)
            if len(T) > len(best_T):
                best_T = T
        #if (typ[0],side[0])==('i','r'): best_T.reverse()
        #if (typ[0],side[0])==('d','l'): best_T.reverse()
        #if typ[0]=='d': best_T.reverse()
        if side[0] == 'r':
            best_T.reverse()
        return best_T

    def maximal_cyclic_decomposition(self, typ='decreasing', side='right', verbose=False):
        r"""
        Find the unique maximal decomposition of ``self`` into cyclically
        decreasing/increasing elements.

        INPUT:

        - ``typ`` -- ``'increasing'`` or ``'decreasing'``
          (default: ``'decreasing'``); chooses whether to find increasing
          or decreasing sets

        - ``side`` -- ``'right'`` or ``'left'`` (default: ``'right'``) chooses
          whether to find maximal sets starting from the left or the right

        - ``verbose`` -- boolean (default: ``False``); print extra information
          while finding the decomposition

        EXAMPLES::

            sage: p = AffinePermutationGroup(['A',7,1])([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.maximal_cyclic_decomposition()
            [[0, 7], [4, 1, 0], [7, 5, 4, 2, 1]]
            sage: p.maximal_cyclic_decomposition(side='left')
            [[1, 0, 7, 5, 4], [1, 0, 5], [2, 1]]
            sage: p.maximal_cyclic_decomposition(typ='increasing', side='right')
            [[1], [5, 0, 1, 2], [4, 5, 7, 0, 1]]
            sage: p.maximal_cyclic_decomposition(typ='increasing', side='left')
            [[0, 1, 2, 4, 5], [4, 7, 0, 1], [7]]

        TESTS::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: S=p.maximal_cyclic_decomposition()
            sage: p==prod(A.from_word(l) for l in S)
            True
            sage: S=p.maximal_cyclic_decomposition(typ='increasing', side='left')
            sage: p==prod(A.from_word(l) for l in S)
            True
            sage: S=p.maximal_cyclic_decomposition(typ='increasing', side='right')
            sage: p==prod(A.from_word(l) for l in S)
            True
            sage: S=p.maximal_cyclic_decomposition(typ='decreasing', side='right')
            sage: p==prod(A.from_word(l) for l in S)
            True
        """
        y = self.clone()
        listy = []
        if verbose:
            print('length of x:', self.length())
        while not y.is_one():
            S = y.maximal_cyclic_factor(typ, side, verbose)
            listy.append(S[:])
            if side[0] == 'r':
                S.reverse()
            for i in S:
                if side[0] == 'r':
                    y = y.apply_simple_reflection_right(i)
                else:
                    y = y.apply_simple_reflection_left(i)
            if verbose:
                print(S, y.length())
        if side[0] == 'r':
            listy.reverse()
        return listy

    def to_lehmer_code(self, typ='decreasing', side='right'):
        r"""
        Return the affine Lehmer code.

        There are four such codes; the options ``typ`` and ``side`` determine
        which code is generated.  The codes generated are the shape of the
        maximal cyclic decompositions of ``self`` according to the given
        ``typ`` and ``side`` options.

        INPUT:

        - ``typ`` -- ``'increasing'`` or ``'decreasing'``
          (default: ``'decreasing'``); chooses whether to find increasing
          or decreasing sets

        - ``side`` -- ``'right'`` or ``'left'`` (default: ``'right'``) chooses
          whether to find maximal sets starting from the left or the right

        EXAMPLES::

            sage: import itertools
            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: orders = ('increasing','decreasing')
            sage: sides = ('left','right')
            sage: for o,s in itertools.product(orders, sides):
            ....:   p.to_lehmer_code(o,s)
            [2, 3, 2, 0, 1, 2, 0, 0]
            [2, 2, 0, 0, 2, 1, 0, 3]
            [3, 1, 0, 0, 2, 1, 0, 3]
            [0, 3, 3, 0, 1, 2, 0, 1]
            sage: for a in itertools.product(orders, sides):
            ....:   A.from_lehmer_code(p.to_lehmer_code(a[0],a[1]), a[0],a[1])==p
            True
            True
            True
            True
        """
        code = [0 for i in range(self.k+1)]
        if typ[0] == 'i' and side[0] == 'r':
            #Find number of positions to the right of position i with smaller
            #value than the number in position i.
            for i in range(self.k+1):
                a = self(i)
                for j in range(i+1, i+self.k+1):
                    b = self(j)
                    if b < a:
                        code[i] += (a-b) // (self.k+1) + 1
        elif typ[0] == 'd' and side[0] == 'r':
            #Find number of positions to the left of position i with larger
            #value than the number in position i.  Then cyclically shift
            #the resulting vector.
            for i in range(self.k+1):
                a = self(i)
                for j in range(i-self.k, i):
                    b = self(j)
                    # A small rotation is necessary for the reduced word from
                    # the Lehmer code to match the element.
                    if a < b:
                        code[i-1] += ((b-a)//(self.k+1)+1)
        elif typ[0] == 'i' and side[0] == 'l':
            #Find number of positions to the right of i smaller than i, then
            #cyclically shift the resulting vector.
            for i in range(self.k+1):
                pos = self.position(i)
                for j in range(pos+1, pos+self.k+1):
                    b = self(j)
                    #A small rotation is necessary for the reduced word from
                    #the lehmer code to match the element.
                    if b < i:
                        code[i-1] += (i-b) // (self.k+1) + 1
        elif typ[0] == 'd' and side[0] == 'l':
            #Find number of positions to the left of i larger than i.
            for i in range(self.k+1):
                pos = self.position(i)
                for j in range(pos-self.k, pos):
                    b = self(j)
                    if b > i:
                        code[i] += (b-i) // (self.k+1) + 1
        return Composition(code)

    def is_fully_commutative(self) -> bool:
        r"""
        Determine whether ``self`` is fully commutative.

        This means that it has no reduced word with a braid.

        This uses a specific algorithm.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p = A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.is_fully_commutative()
            False
            sage: q = A([-3, -2, 0, 7, 9, 2, 11, 12])
            sage: q.is_fully_commutative()
            True
        """
        if self == self.parent().one():
            return True
        c = self.to_lehmer_code()
        firstnonzero = None
        m = -1
        for i in range(self.n):
            if c[i] > 0:
                if firstnonzero is None:
                    firstnonzero = i
                if m != -1 and c[i] - (i - m) >= c[m]:
                    return False
                m = i
        # now check m (the last nonzero) against first nonzero.
        d = self.n - (m - firstnonzero)
        return not c[firstnonzero] - d >= c[m]

    def to_bounded_partition(self, typ='decreasing', side='right'):
        r"""
        Return the `k`-bounded partition associated to the dominant element
        obtained by sorting the Lehmer code.

        INPUT:

        - ``typ`` -- ``'increasing'`` or ``'decreasing'`` (default: ``'decreasing'``);
          chooses whether to find increasing or decreasing sets

        - ``side`` -- ``'right'`` or ``'left'`` (default: ``'right'``); chooses
          whether to find maximal sets starting from the left or the right

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',2,1])
            sage: p=A.from_lehmer_code([4,1,0])
            sage: p.to_bounded_partition()
            [2, 1, 1, 1]
        """
        c = sorted(self.to_lehmer_code(typ, side))
        c.reverse()
        return Partition(c).conjugate()

    def to_core(self, typ='decreasing', side='right'):
        r"""
        Return the core associated to the dominant element obtained by sorting
        the Lehmer code.

        INPUT:

        - ``typ`` -- ``'increasing'`` or ``'decreasing'`` (default: ``'decreasing'``.)

        - ``side`` -- ``'right'`` or ``'left'`` (default: ``'right'``.)  Chooses whether to
          find maximal sets starting from the left or the right

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',2,1])
            sage: p=A.from_lehmer_code([4,1,0])
            sage: p.to_bounded_partition()
            [2, 1, 1, 1]
            sage: p.to_core()
            [4, 2, 1, 1]
        """
        return self.to_bounded_partition(typ,side).to_core(self.k)

    def to_dominant(self, typ='decreasing', side='right'):
        r"""
        Finds the Lehmer code and then sorts it.  Returns the affine permutation
        with the given sorted Lehmer code; this element is 0-dominant.

        INPUT:

        - ``typ`` -- ``'increasing'`` or ``'decreasing'``
          (default: ``'decreasing'``) chooses whether to find increasing
          or decreasing sets

        - ``side`` -- ``'right'`` or ``'left'`` (default: ``'right'``) chooses
          whether to find maximal sets starting from the left or the right

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.to_dominant()
            Type A affine permutation with window [-2, -1, 1, 3, 4, 8, 10, 13]
            sage: p.to_dominant(typ='increasing', side='left')
            Type A affine permutation with window [3, 4, -1, 5, 0, 9, 6, 10]
        """
        if self.is_i_grassmannian(side=side):
            return self
        c = sorted(self.to_lehmer_code(typ,side))
        c.reverse()
        return self.parent().from_lehmer_code(c, typ, side)

    def tableau_of_word(self, w, typ='decreasing', side='right', alpha=None):
        r"""
        Finds a tableau on the Lehmer code of ``self`` corresponding
        to the given reduced word.

        For a full description of this algorithm, see [Den2012]_.

        INPUT:

        - ``w`` -- a reduced word for ``self``
        - ``typ`` -- ``'increasing'`` or ``'decreasing'``; the type of
          Lehmer code used
        - ``side`` -- ``'right'`` or ``'left'``
        - ``alpha`` -- a content vector; ``w`` should be of type ``alpha``;
          specifying ``alpha`` produces semistandard tableaux

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.tableau_of_word(p.reduced_word())
            [[], [1, 6, 9], [2, 7, 10], [], [3], [4, 8], [], [5]]
            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: w=p.reduced_word()
            sage: w
            [0, 7, 4, 1, 0, 7, 5, 4, 2, 1]
            sage: alpha=[5,3,2]
            sage: p.tableau_of_word(p.reduced_word(), alpha=alpha)
            [[], [1, 2, 3], [1, 2, 3], [], [1], [1, 2], [], [1]]
            sage: p.tableau_of_word(p.reduced_word(), side='left')
            [[1, 4, 9], [6], [], [], [3, 7], [8], [], [2, 5, 10]]
            sage: p.tableau_of_word(p.reduced_word(), typ='increasing', side='right')
            [[9, 10], [1, 2], [], [], [3, 4], [8], [], [5, 6, 7]]
            sage: p.tableau_of_word(p.reduced_word(), typ='increasing', side='left')
            [[1, 2], [4, 5, 6], [9, 10], [], [3], [7, 8], [], []]
        """
        g = self.parent().simple_reflections()
        # Check w is reduced....:should probably throw an exception otherwise.
        x0 = prod(g[i] for i in w)
        if x0.length() != len(w):
            raise ValueError("word was not reduced")
        if alpha is None:
            alpha = Composition([1 for i in w])
        else:
            if sum(alpha) != len(w):
                raise ValueError("size of alpha must match length of w")
            alpha = Composition(alpha)
        # TODO: We should probably check that w is of type alpha! probably a different function.
        # Now we actually build the recording tableau.
        tab = [[] for _ in repeat(None, self.k + 1)]
        label = 1
        al_index = 0
        j = 0
        x = self.parent().one()
        cx = x.to_lehmer_code(typ, side)
        n = len(w)-1
        for i in range(len(w)):
            if side[0] == 'r':
                #y=g[w[n-i]]*x
                y = x.apply_simple_reflection_left(w[n-i])
            else:
                y = x.apply_simple_reflection_right(w[i])
            cy = y.to_lehmer_code(typ, side)
            for r in range(self.k+1):
                if cy[r] > cx[r]:
                    tab[r].append(label)
                    j += 1
                    if j == alpha[al_index]:
                        al_index += 1
                        j = 0
                        label += 1
                    break
            x = y
            cx = cy
        return tab

#-------------------------------------------------------------------------------


class AffinePermutationTypeC(AffinePermutation):
    #----------------------
    #Type-specific methods.
    #(Methods existing in all types, but with type-specific definition.)
    #----------------------
    def check(self):
        r"""
        Check that ``self`` is an affine permutation.

        EXAMPLES::

            sage: C = AffinePermutationGroup(['C',4,1])
            sage: x = C([-1,5,3,7])
            sage: x
            Type C affine permutation with window [-1, 5, 3, 7]
        """
        if not self:
            return
        k = self.parent().k
        if len(self) != k:
            raise ValueError("length of list must be k=" + str(k))
        reslist = []
        for i in self:
            r = i % self.N
            if r == 0:
                raise ValueError("entries may not have residue 0 mod 2k+1")
            if not (r not in reslist and self.N-r not in reslist):
                raise ValueError("entries must have distinct residues")
            reslist.append(r)

    def value(self, i):
        r"""
        Return the image of the integer ``i`` under this permutation.

        EXAMPLES::

            sage: C = AffinePermutationGroup(['C',4,1])
            sage: x = C.one()
            sage: [x.value(i) for i in range(-10,10)] == list(range(-10,10))
            True
        """
        N = 2*self.k + 1
        window = i // N
        index = i % N
        if index == 0:
            return i
        if index <= self.k:
            return self[index-1]+window*N
        if index > self.k:
            return -(self[N-index-1]-N)+window*N

    def position(self, i):
        r"""
        Find the position `j` such the ``self.value(j)=i``.

        EXAMPLES::

            sage: C = AffinePermutationGroup(['C',4,1])
            sage: x = C.one()
            sage: [x.position(i) for i in range(-10,10)] == list(range(-10,10))
            True
        """
        N = 2*self.k + 1
        index = i % N
        if index == 0:
            return i
        for r in range(len(self)):
            if self[r] % N == index:
                #i sits in position i, but some number of windows away.
                diff = (i-self[r]) // N
                return r + diff*N + 1
            if self[r] % N == N - index:
                #then we sit some number of windows from position -r.
                diff = (i+self[r]) // N
                return -r + diff*N - 1
        return False

    def apply_simple_reflection_right(self, i):
        r"""
        Apply the simple reflection indexed by ``i`` on positions.

        EXAMPLES::

            sage: C = AffinePermutationGroup(['C',4,1])
            sage: x=C([-1,5,3,7])
            sage: for i in C.index_set(): x.apply_simple_reflection_right(i)
            Type C affine permutation with window [1, 5, 3, 7]
            Type C affine permutation with window [5, -1, 3, 7]
            Type C affine permutation with window [-1, 3, 5, 7]
            Type C affine permutation with window [-1, 5, 7, 3]
            Type C affine permutation with window [-1, 5, 3, 2]
        """
        if i not in self.parent().index_set():
            raise ValueError('index not in index set')
        j = i
        l = self[:]
        if j != 0 and j != self.k:
            a = l[j-1]
            l[j-1] = l[j]
            l[j] = a
        elif j == 0:
            l[0] = -l[0]
        elif j == self.k:
            l[self.k-1] = self(self.k+1)
        #return l
        return type(self)(self.parent(), l, check=False)

    def apply_simple_reflection_left(self, i):
        r"""
        Apply the simple reflection indexed by ``i`` on values.

        EXAMPLES::

            sage: C = AffinePermutationGroup(['C',4,1])
            sage: x = C([-1,5,3,7])
            sage: for i in C.index_set(): x.apply_simple_reflection_left(i)
            Type C affine permutation with window [1, 5, 3, 7]
            Type C affine permutation with window [-2, 5, 3, 8]
            Type C affine permutation with window [-1, 5, 2, 6]
            Type C affine permutation with window [-1, 6, 4, 7]
            Type C affine permutation with window [-1, 4, 3, 7]
        """
        if i not in self.parent().index_set():
            raise ValueError('index not in index set')
        j = self.N - i
        l = []
        if i != self.k and i != 0:
            for m in range(self.k):
                res = self[m] % self.N
                if res == i:
                    l.append(self[m] + 1)
                elif res == i + 1:
                    l.append(self[m] - 1)
                elif res == j:
                    l.append(self[m] - 1)
                elif res == j - 1:
                    l.append(self[m] + 1)
                else:
                    l.append(self[m])
        elif i == 0:
            for m in range(self.k):
                res = self[m] % self.N
                if res == 1:
                    l.append(self[m] - 2)
                elif res == self.N - 1:
                    l.append(self[m] + 2)
                else:
                    l.append(self[m])
        elif i == self.k:
            for m in range(self.k):
                res = self[m] % self.N
                if res == i:
                    l.append(self[m] + 1)
                elif res == j:
                    l.append(self[m] - 1)
                else:
                    l.append(self[m])
        return type(self)(self.parent(), l, check=False)

    def has_right_descent(self, i) -> bool:
        r"""
        Determine whether there is a descent at index ``i``.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: C = AffinePermutationGroup(['C',4,1])
            sage: x = C([-1,5,3,7])
            sage: for i in C.index_set(): x.has_right_descent(i)
            True
            False
            True
            False
            True
        """
        return self.value(i) > self.value(i+1)

    def has_left_descent(self, i) -> bool:
        r"""
        Determine whether there is a descent at ``i``.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: C = AffinePermutationGroup(['C',4,1])
            sage: x = C([-1,5,3,7])
            sage: for i in C.index_set(): x.has_left_descent(i)
            True
            False
            True
            False
            True
        """
        # This is much faster than the default method of taking the inverse and
        # then finding right descents...
        return self.position(i) > self.position(i + 1)

    def to_type_a(self):
        r"""
        Return an embedding of ``self`` into the affine permutation group of
        type `A`.

        EXAMPLES::

            sage: C = AffinePermutationGroup(['C',4,1])
            sage: x = C([-1,5,3,7])
            sage: x.to_type_a()
            Type A affine permutation with window [-1, 5, 3, 7, 2, 6, 4, 10, 9]
        """
        A = AffinePermutationGroup(['A', self.N-1, 1])
        return A([self.value(i) for i in range(1, self.N+1)])


class AffinePermutationTypeB(AffinePermutationTypeC):
    #----------------------
    #Type-specific methods.
    #(Methods existing in all types, but with type-specific definition.)
    #----------------------
    def check(self):
        r"""
        Check that ``self`` is an affine permutation.

        EXAMPLES::

            sage: B = AffinePermutationGroup(['B',4,1])
            sage: x = B([-5,1,6,-2])
            sage: x
            Type B affine permutation with window [-5, 1, 6, -2]
        """
        if not self:
            return
        k = self.parent().k
        # Check window length.
        if len(self) != k:
            raise ValueError("length of list must be k=" + str(k))
        # Check for repeated residues.
        reslist = []
        for i in self:
            r = i % self.N
            if r == 0:
                raise ValueError("entries may not have residue 0 mod 2k+1")
            if not (r not in reslist and self.N - r not in reslist):
                raise ValueError("entries must have distinct residues")
            reslist.append(r)
        # Check that we have an even number of 'small' elements right of the zeroth entry.
        s = sum(-i // self.N+1 for i in (self.value(j) for j in range(1,self.N+1)) if i < 0)
        if s % 2:
            raise ValueError("type B affine permutations have an even number of "
                             "entries less than 0 to the right of the 0th position")

    def apply_simple_reflection_right(self, i):
        r"""
        Apply the simple reflection indexed by ``i`` on positions.

        EXAMPLES::

            sage: B = AffinePermutationGroup(['B',4,1])
            sage: p=B([-5,1,6,-2])
            sage: p.apply_simple_reflection_right(1)
            Type B affine permutation with window [1, -5, 6, -2]
            sage: p.apply_simple_reflection_right(0)
            Type B affine permutation with window [-1, 5, 6, -2]
            sage: p.apply_simple_reflection_right(4)
            Type B affine permutation with window [-5, 1, 6, 11]
        """
        if i not in self.parent().index_set():
            raise ValueError('index not in index set')
        j = i
        l = self[:]
        if j != 0 and j != self.k:
            #just swap l[j], l[j-1]
            (l[j-1], l[j]) = (l[j], l[j-1])
        elif j == 0:
            l[0] = -self(2)
            l[1] = -self(1)
        elif j == self.k:
            l[self.k-1] = self(self.k+1)
        #return l
        return type(self)(self.parent(), l, check=False)

    def apply_simple_reflection_left(self, i):
        r"""
        Apply the simple reflection indexed by ``i`` on values.

        EXAMPLES::

            sage: B = AffinePermutationGroup(['B',4,1])
            sage: p=B([-5,1,6,-2])
            sage: p.apply_simple_reflection_left(0)
            Type B affine permutation with window [-5, -2, 6, 1]
            sage: p.apply_simple_reflection_left(2)
            Type B affine permutation with window [-5, 1, 7, -3]
            sage: p.apply_simple_reflection_left(4)
            Type B affine permutation with window [-4, 1, 6, -2]
        """
        if i not in self.parent().index_set():
            raise ValueError('index not in index set')
        j = self.N - i
        l = []
        if i != self.k and i != 0:
            for m in range(self.k):
                res = self[m] % self.N
                if res == i:
                    l.append(self[m]+1)
                elif res == i + 1:
                    l.append(self[m]-1)
                elif res == j:
                    l.append(self[m]-1)
                elif res == j - 1:
                    l.append(self[m]+1)
                else:
                    l.append(self[m])
        elif i == 0:
            for m in range(self.k):
                res = self[m] % self.N
                if res == 1:
                    l.append(self[m]-3)
                elif res == self.N - 2:
                    l.append(self[m]+3)
                elif res == 2:
                    l.append(self[m]-3)
                elif res == self.N - 1:
                    l.append(self[m]+3)
                else:
                    l.append(self[m])
        elif i == self.k:
            for m in range(self.k):
                res = self[m] % self.N
                if res == i:
                    l.append(self[m] + 1)
                elif res == j:
                    l.append(self[m] - 1)
                else:
                    l.append(self[m])
        return type(self)(self.parent(), l, check=False)

    def has_right_descent(self, i) -> bool:
        r"""
        Determines whether there is a descent at index ``i``.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: B = AffinePermutationGroup(['B',4,1])
            sage: p = B([-5,1,6,-2])
            sage: [p.has_right_descent(i) for i in B.index_set()]
            [True, False, False, True, False]
        """
        if i == 0:
            return self.value(-2) > self.value(1)
        return self.value(i) > self.value(i+1)

    def has_left_descent(self, i) -> bool:
        r"""
        Determines whether there is a descent at ``i``.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: B = AffinePermutationGroup(['B',4,1])
            sage: p=B([-5,1,6,-2])
            sage: [p.has_left_descent(i) for i in B.index_set()]
            [True, True, False, False, True]
        """
        if i == 0:
            return self.position(-2) > self.position(1)
        return self.position(i) > self.position(i+1)


class AffinePermutationTypeD(AffinePermutationTypeC):
    #----------------------
    #Type-specific methods.
    #(Methods existing in all types, but with type-specific definition.)
    #----------------------
    def check(self):
        r"""
        Check that ``self`` is an affine permutation.

        EXAMPLES::

            sage: D = AffinePermutationGroup(['D',4,1])
            sage: p = D([1,-6,5,-2])
            sage: p
            Type D affine permutation with window [1, -6, 5, -2]
        """
        if not self:
            return
        k = self.parent().k
        # Check window length.
        if len(self) != k:
            raise ValueError("length of list must be k=" + str(k))
        #Check for repeated residues.
        reslist = []
        for i in self:
            r = i % self.N
            if r == 0:
                raise ValueError("entries may not have residue 0 mod 2k+1")
            if not (r not in reslist and self.N-r not in reslist):
                raise ValueError("entries must have distinct residues")
            reslist.append(r)
        # Check that we have an even number of 'big' elements left of the kth entry.
        s = sum(i // self.N + 1 - (i % self.N <= self.k)
                for i in (self.value(j) for j in range(-self.k,self.k+1)) if i > self.k)
        if s % 2:
            raise ValueError("type D affine permutations have an even number of entries"
                             " greater than x.k weakly to the left of the x.k position")
        # Check that we have an even number of 'small' elements right of the zeroth entry.
        s = sum(-i // self.N+1 for i in (self.value(j) for j in range(1,self.N+1)) if i < 0)
        if s % 2:
            raise ValueError("type D affine permutations have an even number of entries"
                             " less than 0 to the right of the 0th position")

    def apply_simple_reflection_right(self, i):
        r"""
        Apply the simple reflection indexed by ``i`` on positions.

        EXAMPLES::

            sage: D = AffinePermutationGroup(['D',4,1])
            sage: p=D([1,-6,5,-2])
            sage: p.apply_simple_reflection_right(0)
            Type D affine permutation with window [6, -1, 5, -2]
            sage: p.apply_simple_reflection_right(1)
            Type D affine permutation with window [-6, 1, 5, -2]
            sage: p.apply_simple_reflection_right(4)
            Type D affine permutation with window [1, -6, 11, 4]
        """
        if i not in self.parent().index_set():
            raise ValueError('index not in index set')
        j = i
        l = self[:]
        if j != 0 and j != self.k:
            a = l[j-1]
            l[j-1] = l[j]
            l[j] = a
        elif j == 0:
            c = l[0]
            l[0] = -l[1]
            l[1] = -c
        elif j == self.k:
            l[self.k-2] = self(self.k+1)
            l[self.k-1] = self(self.k+2)
        #return l
        return type(self)(self.parent(), l, check=False)

    def apply_simple_reflection_left(self, i):
        r"""
        Apply simple reflection indexed by ``i`` on values.

        EXAMPLES::

            sage: D = AffinePermutationGroup(['D',4,1])
            sage: p=D([1,-6,5,-2])
            sage: p.apply_simple_reflection_left(0)
            Type D affine permutation with window [-2, -6, 5, 1]
            sage: p.apply_simple_reflection_left(1)
            Type D affine permutation with window [2, -6, 5, -1]
            sage: p.apply_simple_reflection_left(4)
            Type D affine permutation with window [1, -4, 3, -2]
        """
        if i not in self.parent().index_set():
            raise ValueError('index not in index set')
        j = self.N - i
        l = []
        if i and i != self.k:
            for m in range(self.k):
                res = self[m] % self.N
                if res == i:
                    l.append(self[m]+1)
                elif res == i+1:
                    l.append(self[m]-1)
                elif res == j:
                    l.append(self[m]-1)
                elif res == j-1:
                    l.append(self[m]+1)
                else:
                    l.append(self[m])
        elif i == 0:
            for m in range(self.k):
                res = self[m] % self.N
                if res == 1:
                    l.append(self[m]-3)
                elif res == self.N-2:
                    l.append(self[m]+3)
                elif res == 2:
                    l.append(self[m]-3)
                elif res == self.N-1:
                    l.append(self[m]+3)
                else:
                    l.append(self[m])
        elif i == self.k:
            for m in range(self.k):
                res = self[m] % self.N
                if res == self.k:
                    l.append(self[m]+2)
                elif res == self.k+2:
                    l.append(self[m]-2)
                elif res == self.k-1:
                    l.append(self[m]+2)
                elif res == self.k+1:
                    l.append(self[m]-2)
                else:
                    l.append(self[m])
        return type(self)(self.parent(), l, check=False)

    def has_right_descent(self, i) -> bool:
        r"""
        Determine whether there is a descent at index ``i``.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: D = AffinePermutationGroup(['D',4,1])
            sage: p=D([1,-6,5,-2])
            sage: [p.has_right_descent(i) for i in D.index_set()]
            [True, True, False, True, False]
        """
        if i == 0:
            return self.value(-2) > self.value(1)
        if i == self.k:
            return self.value(i) > self.value(i+2)
        return self.value(i) > self.value(i+1)

    def has_left_descent(self, i) -> bool:
        r"""
        Determine whether there is a descent at ``i``.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: D = AffinePermutationGroup(['D',4,1])
            sage: p=D([1,-6,5,-2])
            sage: [p.has_left_descent(i) for i in D.index_set()]
            [True, True, False, True, True]
        """
        if i == 0:
            return self.position(-2) > self.position(1)
        if i == self.k:
            return self.position(i) > self.position(i+2)
        return self.position(i) > self.position(i+1)


class AffinePermutationTypeG(AffinePermutation):
    #----------------------
    #Type-specific methods.
    #(Methods existing in all types, but with type-specific definition.)
    #----------------------
    def check(self):
        r"""
        Check that ``self`` is an affine permutation.

        EXAMPLES::

            sage: G = AffinePermutationGroup(['G',2,1])
            sage: p = G([2, 10, -5, 12, -3, 5])
            sage: p
            Type G affine permutation with window [2, 10, -5, 12, -3, 5]
        """
        if not self:
            return
        if not len(self) == 6:
            raise ValueError("length of list must be 6")
        #Check that we have an even number of 'big' elements left of the 7th entry.
        s = sum(i//6 - (i % 6 == 0) for i in self if i > 6)
        if s % 2:
            raise ValueError("type G affine permutations have an even number of"
                             " entries greater than 6 to the left of the 7th position")
        #Check that we have an even number of 'small' elements right of the zeroth entry.
        s = sum(-i//6 + 1 for i in self if i <= 0)
        if s % 2:
            raise ValueError("type G affine permutations have an even number of"
                             " entries less than 0 to the right of the 0th position")

    def value(self, i, base_window=False):
        r"""
        Return the image of the integer ``i`` under this permutation.

        INPUT:

        - ``base_window`` -- boolean indicating whether ``i`` is between 1 and
          `k+1`; if ``True``, will run a bit faster, but the method will screw
          up if ``i`` is not actually in the index set

        EXAMPLES::

            sage: G = AffinePermutationGroup(['G',2,1])
            sage: p=G([2, 10, -5, 12, -3, 5])
            sage: [p.value(i) for i in [1..12]]
            [2, 10, -5, 12, -3, 5, 8, 16, 1, 18, 3, 11]
        """
        N = 6
        if base_window:
            self[i-1]
        window = (i-1) // N
        return self[(i-1) % N] + window*(N)

    def position(self, i):
        r"""
        Find the position ``j`` such the ``self.value(j) == i``.

        EXAMPLES::

            sage: G = AffinePermutationGroup(['G',2,1])
            sage: p = G([2, 10, -5, 12, -3, 5])
            sage: [p.position(i) for i in p]
            [1, 2, 3, 4, 5, 6]
        """
        N = 6
        for r in range(N):
            if self[r] % N == i % N:
                #i sits in position i, but some number of windows away.
                diff = (i-self[r]) // N
                return r + diff*N + 1
        return False

    def apply_simple_reflection_right(self, i):
        r"""
        Apply the simple reflection indexed by ``i`` on positions.

        EXAMPLES::

            sage: G = AffinePermutationGroup(['G',2,1])
            sage: p = G([2, 10, -5, 12, -3, 5])
            sage: p.apply_simple_reflection_right(0)
            Type G affine permutation with window [-9, -1, -5, 12, 8, 16]
            sage: p.apply_simple_reflection_right(1)
            Type G affine permutation with window [10, 2, 12, -5, 5, -3]
            sage: p.apply_simple_reflection_right(2)
            Type G affine permutation with window [2, -5, 10, -3, 12, 5]
        """
        if i not in self.parent().index_set():
            raise ValueError('index not in index set')
        j = i
        l = self[:]
        if j == 1:
            l[0] = self(2)
            l[1] = self(1)
            l[2] = self(4)
            l[3] = self(3)
            l[4] = self(6)
            l[5] = self(5)
        elif j == 2:
            l[1] = self(3)
            l[2] = self(2)
            l[3] = self(5)
            l[4] = self(4)
        elif j == 0:
            l[0] = self(-1)
            l[1] = self(0)
            l[4] = self(7)
            l[5] = self(8)
        #return l
        return type(self)(self.parent(), l, check=False)

    def apply_simple_reflection_left(self, i):
        r"""
        Apply simple reflection indexed by `i` on values.

        EXAMPLES::

            sage: G = AffinePermutationGroup(['G',2,1])
            sage: p=G([2, 10, -5, 12, -3, 5])
            sage: p.apply_simple_reflection_left(0)
            Type G affine permutation with window [0, 10, -7, 14, -3, 7]
            sage: p.apply_simple_reflection_left(1)
            Type G affine permutation with window [1, 9, -4, 11, -2, 6]
            sage: p.apply_simple_reflection_left(2)
            Type G affine permutation with window [3, 11, -5, 12, -4, 4]
        """
        if i not in self.parent().index_set():
            raise ValueError('index not in index set')
        l = []
        if i == 1:
            for m in range(6):
                res = self[m] % 6
                if res == 1 or res == 3 or res == 5:
                    l.append(self[m]+1)
                elif res == 2 or res == 4 or res == 0:
                    l.append(self[m]-1)
                else:
                    l.append(self[m])
        elif i == 2:
            for m in range(6):
                res = self[m] % 6
                if res == 2 or res == 4:
                    l.append(self[m]+1)
                elif res == 3 or res == 5:
                    l.append(self[m]-1)
                else:
                    l.append(self[m])
        elif i == 0:
            for m in range(6):
                res = self[m] % 6
                if res == 1 or res == 2:
                    l.append(self[m]-2)
                elif res == 5 or res == 0:
                    l.append(self[m]+2)
                else:
                    l.append(self[m])
        return type(self)(self.parent(), l, check=False)

    def has_right_descent(self, i) -> bool:
        r"""
        Determines whether there is a descent at index `i`.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: G = AffinePermutationGroup(['G',2,1])
            sage: p = G([2, 10, -5, 12, -3, 5])
            sage: [p.has_right_descent(i) for i in G.index_set()]
            [False, False, True]
        """
        if i not in self.parent().index_set():
            raise ValueError('index not in index set')
        if i == 0:
            return self.value(0) > self.value(2)
        return self.value(i) > self.value(i+1)

    def has_left_descent(self, i) -> bool:
        r"""
        Determines whether there is a descent at ``i``.

        INPUT:

        - ``i`` -- integer

        EXAMPLES::

            sage: G = AffinePermutationGroup(['G',2,1])
            sage: p = G([2, 10, -5, 12, -3, 5])
            sage: [p.has_left_descent(i) for i in G.index_set()]
            [False, True, False]
        """
        if i not in self.parent().index_set():
            raise ValueError('index not in index set')
        if i == 0:
            return self.position(0) > self.position(2)
        return self.position(i) > self.position(i+1)

    def to_type_a(self):
        r"""
        Return an embedding of ``self`` into the affine permutation group of
        type A.

        EXAMPLES::

            sage: G = AffinePermutationGroup(['G',2,1])
            sage: p = G([2, 10, -5, 12, -3, 5])
            sage: p.to_type_a()
            Type A affine permutation with window [2, 10, -5, 12, -3, 5]
        """
        A = AffinePermutationGroup(['A', 5, 1])
        return A(self)


#-------------------------------------------------------------------------
#    Class of all affine permutations.
#-------------------------------------------------------------------------

def AffinePermutationGroup(cartan_type):
    r"""
    Wrapper function for specific affine permutation groups.

    These are combinatorial implementations of the affine Weyl groups of
    types `A`, `B`, `C`, `D`, and `G` as permutations of the set of all integers.
    the basic algorithms are derived from [BB2005]_ and [Eri1995]_.

    EXAMPLES::

        sage: ct = CartanType(['A',7,1])
        sage: A = AffinePermutationGroup(ct)
        sage: A
        The group of affine permutations of type ['A', 7, 1]

    We define an element of ``A``::

        sage: p = A([3, -1, 0, 6, 5, 4, 10, 9])
        sage: p
        Type A affine permutation with window [3, -1, 0, 6, 5, 4, 10, 9]

    We find the value `p(1)`, considering `p` as a bijection on the integers.
    This is the same as calling the
    :meth:`~sage.combinat.affine_permutation.AffinePermutation.value` method::

        sage: p.value(1)
        3
        sage: p(1) == p.value(1)
        True

    We can also find the position of the integer 3 in `p` considered as a
    sequence, equivalent to finding `p^{-1}(3)`::

        sage: p.position(3)
        1
        sage: (p^-1)(3)
        1

    Since the affine permutation group is a group, we demonstrate its
    group properties::

        sage: A.one()
        Type A affine permutation with window [1, 2, 3, 4, 5, 6, 7, 8]

        sage: q = A([0, 2, 3, 4, 5, 6, 7, 9])
        sage: p * q
        Type A affine permutation with window [1, -1, 0, 6, 5, 4, 10, 11]
        sage: q * p
        Type A affine permutation with window [3, -1, 1, 6, 5, 4, 10, 8]

        sage: p^-1
        Type A affine permutation with window [0, -1, 1, 6, 5, 4, 10, 11]
        sage: p^-1 * p == A.one()
        True
        sage: p * p^-1 == A.one()
        True

    If we decide we prefer the Weyl Group implementation of the affine Weyl
    group, we can easily get it::

        sage: p.to_weyl_group_element()
        [ 0 -1  0  1  0  0  1  0]
        [ 1 -1  0  1  0  0  1 -1]
        [ 1 -1  0  1  0  0  0  0]
        [ 0  0  0  1  0  0  0  0]
        [ 0  0  0  1  0 -1  1  0]
        [ 0  0  0  1 -1  0  1  0]
        [ 0  0  0  0  0  0  1  0]
        [ 0 -1  1  0  0  0  1  0]

    We can find a reduced word and do all of the other things one expects in
    a Coxeter group::

        sage: p.has_right_descent(1)
        True
        sage: p.apply_simple_reflection(1)
        Type A affine permutation with window [-1, 3, 0, 6, 5, 4, 10, 9]
        sage: p.apply_simple_reflection(0)
        Type A affine permutation with window [1, -1, 0, 6, 5, 4, 10, 11]
        sage: p.reduced_word()
        [0, 7, 4, 1, 0, 7, 5, 4, 2, 1]
        sage: p.length()
        10

    The following methods are particular to type `A`.
    We can check if the element is fully commutative::

        sage: p.is_fully_commutative()
        False
        sage: q.is_fully_commutative()
        True

    We can also compute the affine Lehmer code of the permutation,
    a weak composition with `k + 1` entries::

        sage: p.to_lehmer_code()
        [0, 3, 3, 0, 1, 2, 0, 1]

    Once we have the Lehmer code, we can obtain a `k`-bounded partition by
    sorting the Lehmer code, and then reading the row lengths.
    There is a unique 0-Grassmanian (dominant) affine permutation associated
    to this `k`-bounded partition, and a `k`-core as well. ::

        sage: p.to_bounded_partition()
        [5, 3, 2]
        sage: p.to_dominant()
        Type A affine permutation with window [-2, -1, 1, 3, 4, 8, 10, 13]
        sage: p.to_core()
        [5, 3, 2]

    Finally, we can take a reduced word for `p` and insert it to find a
    standard composition tableau associated uniquely to that word::

        sage: p.tableau_of_word(p.reduced_word())
        [[], [1, 6, 9], [2, 7, 10], [], [3], [4, 8], [], [5]]

    We can also form affine permutation groups in types `B`, `C`, `D`,
    and `G`::

        sage: B = AffinePermutationGroup(['B',4,1])
        sage: B.an_element()
        Type B affine permutation with window [-1, 3, 4, 11]

        sage: C = AffinePermutationGroup(['C',4,1])
        sage: C.an_element()
        Type C affine permutation with window [2, 3, 4, 10]

        sage: D = AffinePermutationGroup(['D',4,1])
        sage: D.an_element()
        Type D affine permutation with window [-1, 3, 11, 5]

        sage: G = AffinePermutationGroup(['G',2,1])
        sage: G.an_element()
        Type G affine permutation with window [0, 4, -1, 8, 3, 7]
    """
    ct = CartanType(cartan_type)
    if ct.letter == 'A':
        return AffinePermutationGroupTypeA(ct)
    if ct.letter == 'B':
        return AffinePermutationGroupTypeB(ct)
    if ct.letter == 'C':
        return AffinePermutationGroupTypeC(ct)
    if ct.letter == 'D':
        return AffinePermutationGroupTypeD(ct)
    if ct.letter == 'G':
        return AffinePermutationGroupTypeG(ct)
    raise NotImplementedError('Cartan type provided is not implemented as an affine permutation group')


class AffinePermutationGroupGeneric(UniqueRepresentation, Parent):
    """
    The generic affine permutation group class, in which we define all type-free
    methods for the specific affine permutation groups.
    """

    # ----------------------
    # Type-free methods.
    # ----------------------

    def __init__(self, cartan_type):
        r"""
        TESTS::

            sage: AffinePermutationGroup(['A',7,1])
            The group of affine permutations of type ['A', 7, 1]
        """
        Parent.__init__(self, category=AffineWeylGroups())
        ct = CartanType(cartan_type)
        self.k = Integer(ct.n)
        self.n = ct.rank()
        # This N doesn't matter for type A, but comes up in all other types.
        if ct.letter == 'A':
            self.N = self.k + 1
        elif ct.letter == 'B' or ct.letter == 'C' or ct.letter == 'D':
            self.N = 2 * self.k + 1
        elif ct.letter == 'G':
            self.N = 6
        self._cartan_type = ct

    def _element_constructor_(self, *args, **keywords):
        r"""
        TESTS::

            sage: AffinePermutationGroup(['A',7,1])([3, -1, 0, 6, 5, 4, 10, 9])
            Type A affine permutation with window [3, -1, 0, 6, 5, 4, 10, 9]
        """
        return self.element_class(self, *args, **keywords)

    def _repr_(self) -> str:
        r"""
        TESTS::

            sage: AffinePermutationGroup(['A',7,1])
            The group of affine permutations of type ['A', 7, 1]
        """
        return "The group of affine permutations of type " + str(self.cartan_type())

    def _test_enumeration(self, n=4, **options):
        r"""
        Test that ``self`` has same number of elements of length ``n`` as the
        Weyl Group implementation.

        Combined with ``self._test_coxeter_relations`` this shows isomorphism
        up to length ``n``.

        TESTS::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: A._test_enumeration(3)
        """
        tester = self._tester(**options)
        n1 = len(list(self.elements_of_length(n)))
        W = self.weyl_group()
        I = W.weak_order_ideal(ConstantFunction(True), side='right')
        n2 = len(list(I.elements_of_depth_iterator(n)))
        tester.assertEqual(n1, n2, "number of (ranked) elements of affine"
                                   " permutation group disagrees with Weyl group")

    def weyl_group(self):
        r"""
        Return the Weyl Group of the same type as ``self``.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: A.weyl_group()
            Weyl Group of type ['A', 7, 1] (as a matrix group acting on the root space)
        """
        return WeylGroup(self._cartan_type)

    def classical(self):
        r"""
        Return the finite permutation group.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: A.classical()
            Symmetric group of order 8! as a permutation group
        """
        if self._cartan_type.letter == 'A':
            return SymmetricGroup(self.k+1)
        return WeylGroup(self._cartan_type.classical())

    def cartan_type(self):
        r"""
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: AffinePermutationGroup(['A',7,1]).cartan_type()
            ['A', 7, 1]
        """
        return self._cartan_type

    def cartan_matrix(self):
        r"""
        Return the Cartan matrix of ``self``.

        EXAMPLES::

            sage: AffinePermutationGroup(['A',7,1]).cartan_matrix()
            [ 2 -1  0  0  0  0  0 -1]
            [-1  2 -1  0  0  0  0  0]
            [ 0 -1  2 -1  0  0  0  0]
            [ 0  0 -1  2 -1  0  0  0]
            [ 0  0  0 -1  2 -1  0  0]
            [ 0  0  0  0 -1  2 -1  0]
            [ 0  0  0  0  0 -1  2 -1]
            [-1  0  0  0  0  0 -1  2]
        """
        return self.cartan_type().cartan_matrix()

    def is_crystallographic(self):
        r"""
        Tells whether the affine permutation group is crystallographic.

        EXAMPLES::

            sage: AffinePermutationGroup(['A',7,1]).is_crystallographic()
            True
        """
        return self.cartan_type().is_crystallographic()

    def index_set(self):
        r"""
        EXAMPLES::

            sage: AffinePermutationGroup(['A',7,1]).index_set()
            (0, 1, 2, 3, 4, 5, 6, 7)
        """
        return self.cartan_type().index_set()

    _index_set = index_set

    def reflection_index_set(self):
        r"""
        EXAMPLES::

            sage: AffinePermutationGroup(['A',7,1]).reflection_index_set()
            (0, 1, 2, 3, 4, 5, 6, 7)
        """
        return self.cartan_type().index_set()

    def rank(self):
        r"""
        Rank of the affine permutation group, equal to `k+1`.

        EXAMPLES::

            sage: AffinePermutationGroup(['A',7,1]).rank()
            8
        """
        return self.k + 1

    def random_element(self, n=None):
        r"""
        Return a random affine permutation of length ``n``.

        If ``n`` is not specified, then ``n`` is chosen as a random
        nonnegative integer in `[0, 1000]`.

        Starts at the identity, then chooses an upper cover at random.
        Not very uniform: actually constructs a uniformly random reduced word
        of length `n`.  Thus we most likely get elements with lots of reduced
        words!

        For the actual code, see
        :meth:`sage.categories.coxeter_group.random_element_of_length`.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: A.random_element() # random
            Type A affine permutation with window [-12, 16, 19, -1, -2, 10, -3, 9]
            sage: p = A.random_element(10)
            sage: p.length() == 10
            True
        """
        if n is None:
            n = randint(0, 1000)
        return self.random_element_of_length(n)

    def from_word(self, w):
        r"""
        Build an affine permutation from a given word.
        Note: Already in category as ``from_reduced_word``, but this is less
        typing!

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: A.from_word([0, 7, 4, 1, 0, 7, 5, 4, 2, 1])
            Type A affine permutation with window [3, -1, 0, 6, 5, 4, 10, 9]
        """
        return self.from_reduced_word(w)

    @cached_method
    def _an_element_(self):
        r"""
        Return a Coxeter element.

        EXAMPLES::

            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p=A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: A.from_word([0, 7, 4, 1, 0, 7, 5, 4, 2, 1])
            Type A affine permutation with window [3, -1, 0, 6, 5, 4, 10, 9]
        """
        return self.from_reduced_word(self.index_set())


class AffinePermutationGroupTypeA(AffinePermutationGroupGeneric):
    #------------------------
    #Type-specific methods.
    #(Methods in all types, but with specific definition.)
    #------------------------

    @cached_method
    def one(self):
        r"""
        Return the identity element.

        EXAMPLES::

            sage: AffinePermutationGroup(['A',7,1]).one()
            Type A affine permutation with window [1, 2, 3, 4, 5, 6, 7, 8]

        TESTS::

            sage: A = AffinePermutationGroup(['A',5,1])
            sage: A==loads(dumps(A))
            True
            sage: TestSuite(A).run()
        """
        return self(list(range(1, self.k + 2)))

    #------------------------
    #Type-unique methods.
    #(Methods which do not exist in all types.)
    #------------------------
    def from_lehmer_code(self, C, typ='decreasing', side='right'):
        r"""
        Return the affine permutation with the supplied Lehmer code (a weak
        composition with `k+1` parts, at least one of which is 0).

        INPUT:

        - ``typ`` -- ``'increasing'`` or ``'decreasing'``
          (default: ``'decreasing'``); type of product
        - ``side`` -- ``'right'`` or ``'left'`` (default: ``'right'``);
          whether the decomposition is from the right or left

        EXAMPLES::

            sage: import itertools
            sage: A = AffinePermutationGroup(['A',7,1])
            sage: p = A([3, -1, 0, 6, 5, 4, 10, 9])
            sage: p.to_lehmer_code()
            [0, 3, 3, 0, 1, 2, 0, 1]
            sage: A.from_lehmer_code(p.to_lehmer_code()) == p
            True
            sage: orders = ('increasing','decreasing')
            sage: sides = ('left','right')
            sage: all(A.from_lehmer_code(p.to_lehmer_code(o,s),o,s) == p
            ....:     for o,s in itertools.product(orders,sides))
            True
        """
        if len(C) - 1 != self.k:
            raise ValueError("composition must have {} entries".format(self.k+1))
        if 0 not in C:
            raise ValueError("composition must contain a zero entry")
        k = self.k
        #Find a zero entry in C.
        for r in range(self.k+1):
            if C[r] == 0:
                break
        D = list(C)
        #The s0 and t0 are +-1, dependent on typ and side.
        if (typ[0],side[0]) == ('d','r'):
            (t0,s0) = (-1, 1)
        if (typ[0],side[0]) == ('i','r'):
            (t0,s0) = ( 1, 1)
        if (typ[0],side[0]) == ('d','l'):
            (t0,s0) = (-1,-1)
        if (typ[0],side[0]) == ('i','l'):
            (t0,s0) = ( 1,-1)
        row = 0
        #Method is to build a reduced word from the composition.
        #We create a list of cyclically in/decreasing words appearing in
        #the decomposition corresponding to the composition C,
        #and then build the element.
        listy = []
        while sum(D) > 0:
            l = ['x'] * (self.k + 1)
            ll = []
            #read off a row of C.
            for j in range(self.k+1):
                pos = (r + s0*t0*j) % (k+1)
                residue = (r + s0*t0*(row + j)) % (k+1)
                if D[pos] != 0:
                    ll.append(residue)
                    l[pos] = [residue]
                    D[pos] -= 1
            if side[0] == 'l':
                ll.reverse()
            listy.append(ll)
            row += 1
        if side[0] == 'r':
            listy.reverse()
        x = self.one()
        for ll in listy:
            for i in ll:
                x = x.apply_simple_reflection_right(i)
        return x

    Element = AffinePermutationTypeA


class AffinePermutationGroupTypeC(AffinePermutationGroupGeneric):
    #------------------------
    #Type-specific methods.
    #(Methods in all types, but with specific definition.)
    #------------------------

    @cached_method
    def one(self):
        r"""
        Return the identity element.

        EXAMPLES::

            sage: ct=CartanType(['C',4,1])
            sage: C = AffinePermutationGroup(ct)
            sage: C.one()
            Type C affine permutation with window [1, 2, 3, 4]
            sage: C.one()*C.one()==C.one()
            True

        TESTS::

            sage: C = AffinePermutationGroup(['C',4,1])
            sage: C==loads(dumps(C))
            True
            sage: TestSuite(C).run()
        """
        return self(list(range(1, self.k + 1)))

    Element = AffinePermutationTypeC


class AffinePermutationGroupTypeB(AffinePermutationGroupTypeC):
    #------------------------
    #Type-specific methods.
    #(Methods in all types, but with specific definition.)
    #------------------------
    Element = AffinePermutationTypeB


class AffinePermutationGroupTypeD(AffinePermutationGroupTypeC):
    #------------------------
    #Type-specific methods.
    #(Methods in all types, but with specific definition.)
    #------------------------
    Element = AffinePermutationTypeD


class AffinePermutationGroupTypeG(AffinePermutationGroupGeneric):
    #------------------------
    #Type-specific methods.
    #(Methods in all types, but with specific definition.)
    #------------------------
    @cached_method
    def one(self):
        r"""
        Return the identity element.

        EXAMPLES::

            sage: AffinePermutationGroup(['G',2,1]).one()
            Type G affine permutation with window [1, 2, 3, 4, 5, 6]

        TESTS::

            sage: G = AffinePermutationGroup(['G',2,1])
            sage: G==loads(dumps(G))
            True
            sage: TestSuite(G).run()
        """
        return self([1,2,3,4,5,6])

    Element = AffinePermutationTypeG
