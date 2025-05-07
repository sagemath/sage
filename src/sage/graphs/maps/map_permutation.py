from sage.all import Permutation
from sage.graphs.maps.map_error import InvalidMapPermutationArgument
import numpy as np


class MapPermutation:
    """
    A class to represent permutation used in different maps
    """

    def __init__(self, lst, trust=False) -> None:
        """
        Init the MapPermutation

        INPUT:
        - lst a list representing the permutation or a list of tuples representing 
        the cycle of the permutationor an integer representing the size 
        of the permutation(In This case it will return the identify of size lst) 
        or a Permutation.

        -trust a parameter telling whether or not
        to pass some test to verify if lst represent 
        a valid permutation.

        EXAMPLES::

            sage: MapPermutation(3)
            [1, 2, 3]
            sage: MapPermutation([2,3,1])
            [2, 3, 1]
            sage: MapPermutation([(2,3) , (1,4)])
            [4, 3, 2, 1]

        NOTE::
            O(n) where n is the size of the permutation
        """
        if isinstance(lst, Permutation):
            self._init_from_permutation(lst, trust=trust)
            return
        try:
            if lst == int(lst) or lst > 0:
                self._init_from_number(lst)
                return
        except Exception as _:
            pass

        try:
            if type(lst[0]) == type((42,)):
                self._init_from_cycle_list(lst, trust=trust)
                return
            self._init_from_list(lst)
        except Exception as _:
            raise InvalidMapPermutationArgument()

    def _init_from_cycle_list(self, lst, trust=False):
        """
        Init the permutation from a list of cycles

        INPUT:

            -lst a list of tuples representing the cycle of self
            -trust a boolean indicating if check should be passed

        EXAMPLES::

            sage: lst = [(4,3),(2,1)]
            sage: MapPermutation(lst)._init_from_cycle_list(lst)

        .. NOTE::
            O(m),Used internally not intended to be used by the user
        """
        cnt = np.max(np.array(list(map(lambda x: np.max(np.array(x)), lst))))
        self._tab = np.arange(1, cnt+1, dtype=int)
        self._rtab = np.arange(1, cnt+1, dtype=int)
        for e in lst:
            prev = 0
            for i in e:
                if i <= 0:
                    raise InvalidMapPermutationArgument()

                if prev != 0:
                    self._rtab[i-1] = prev
                    self._tab[prev-1] = i
                prev = i

            self._tab[e[-1]-1] = e[0]
            self._rtab[e[0]-1] = e[-1]

        if trust:
            return

        if (self._tab == 0).sum() != 0:
            raise InvalidMapPermutationArgument()

    def _init_from_number(self, n):
        """
        INPUT:
            n the size of the identify permutation

        EXAMPLES::

            sage: MapPermutation(3)._init_from_number(3)

        .. NOTE::
            O(n),Used internally not intended to be used by the user
        """
        self._tab = np.arange(1, n+1)
        self._rtab = np.arange(1, n+1)

    def _init_from_permutation(self, perm, trust=False):
        """
        INPUT:
            -perm a Permutation object
            -trust a boolean indicating if check should be
            passed

        EXAMPLES::

            sage: MapPermutation(3)._init_from_number(3)

        .. NOTE::
            O(n),Used internally not intended to be used by the user
        """

        self._init_from_list(list(perm), trust=trust)

    def _init_from_list(self, list, trust=False):
        """
        Init the permutation from a list of index

        INPUT:
            -lst a list of index
            -trust a boolean indicating if check should 
            be passed

        EXAMPLES::

            sage: lst = [(4,3),(2,1)]
            sage: MapPermutation(lst)._init_from_cycle_list(lst)

        .. NOTE::
            O(m),Used internally not intended to be used by the user

        """
        self._tab = np.array(list)
        self._rtab = np.zeros(len(list), dtype=int)
        self._rtab[self._tab-1] = np.arange(1, len(list)+1)
        if trust:
            return
        if not np.issubdtype(self._tab.dtype, np.integer) or ((self._tab > len(self._tab)) + (self._tab <= 0)).sum() != 0 or len(np.unique(self._tab)) != len(self._tab):
            raise InvalidMapPermutationArgument()

    def size(self):
        """
        OUTPUT: 
            The size of self

        EXAMPLES::

            sage: MapPermutation([(3,32),(6,7)]).size()
            32

        .. NOTE::
            O(1)
        """
        return len(self._tab)

    def apply(self, i):
        """
        INPUT:
            i an index

        OUTPUT:
            self(i)

        EXAMPLES::

            sage: MapPermutation([(3,22),(22,33),(7,14)]).apply(3)
            22

        .. NOTE::
            O(1)
        """
        if i > self.size():
            return i
        return self._tab[i-1]

    def inverseApply(self, i):
        """
        INPUT:
            i an index

        OUTPUT:
            j such that self(j) = i

        EXAMPLES::

            sage: MapPermutation([(3,22),(22,33),(7,14,8)]).inverseApply(7)
            8

        .. NOTE::
            O(1)
        """

        if i > self.size():
            return i
        return self._rtab[i-1]

    def __repr__(self) -> str:
        """
        OUTPUT:
            The string representation of self

        EXAMPLES::

            sage: MapPermutation([(3,9),(2,1),(7,14,8)])
            [2, 1, 9, 4, 5, 6, 14, 7, 3, 10, 11, 12, 13, 8] 

        """
        return str(list(self))

    def pretty_repr(self):
        """
        OUTPUT:

            A more pretty string representation of self

        EXAMPLES::

            sage: MapPermutation([(3,9),(2,1),(7,14,8)]).pretty_repr()
            'MapPermutation: [(1, 2), (3, 9), (4,), (5,), (6,), (7, 14, 8), (10,), (11,), (12,), (13,)]'

        """
        return f"MapPermutation: {self.to_cycles()}"

    def pretty_print(self):
        """
        Print self in a more pretty form

        EXAMPLES::

            sage: MapPermutation([(3,9),(2,1),(7,14,8)]).pretty_print()
            MapPermutation: [(1, 2), (3, 9), (4,), (5,), (6,), (7, 14, 8), (10,), (11,), (12,), (13,)]

        """
        print(self.pretty_repr())

    def to_cycles(self):
        """
        OUTPUT:
            A list of tuple representing the cycle decomposition of self

        .. NOTE::
            O(n)
        """
        check = np.zeros(self.size()+1)
        cycles = []
        for i in range(1, self.size()+1):
            if check[i]:
                continue
            check[i] = 1
            cycle = [i]
            curIndex = self(i)
            while curIndex != i:
                check[curIndex] = 1
                cycle.append(curIndex)
                curIndex = self(curIndex)
            cycles.append(tuple(cycle))
        return cycles

    def inverse(self):
        """
        This function calculate  the inverse of self

        OUTPUT:
            - The inverse of self

        EXAMPLES::

            sage: MapPermutation([(3,9),(2,1),(7,14,8)]).inverse()
            [2, 1, 9, 4, 5, 6, 8, 14, 3, 10, 11, 12, 13, 7]

        .. NOTE::
            O(n),where n is the number of element of the permutation
        """
        return MapPermutation(list(self._rtab))

    def __call__(self, i):
        """
        OUTPUT:

            self(i)

        EXAMPLES::

            sage: MapPermutation([(3,9),(2,1),(7,14,8)])(3)
            9 

        """
        return self.apply(i)

    def __getitem__(self, id):
        """
        INPUT:

            id

        OUTPUT:
            id th item

        EXAMPLES::
            sage: MapPermutation([(3,9),(2,1),(7,14,8)]).__getitem__(3)
            9 
        """
        return self(id)

    def __iter__(self):
        """

        OUTPUT:
            self, with self.index to 0

        EXAMPLES::

            sage: MapPermutation([(3,9),(2,1),(7,14,8)]).__iter__()
            [2, 1, 9, 4, 5, 6, 14, 7, 3, 10, 11, 12, 13, 8]

        """
        self.index = 0
        return self

    def __next__(self):
        """

        OUTPUT:
            The next value in during an iteration or raise StopIteration

        EXAMPLES::

            sage: MapPermutation([(3,9),(2,1),(7,14,8)]).__iter__().__next__()
            2

        """

        if self.index < len(self):
            result = self(self.index + 1)
            self.index += 1
            return result
        else:
            raise StopIteration

    def __len__(self):
        """

        OUTPUT:
            The size of self

        EXAMPLES::

            sage: len(MapPermutation([(3,9),(2,1),(7,14,8)]))
            14 

        """
        return self.size()

    def __eq__(self, other):
        """
        INPUT:
            other, an object

        OUTPUT:
            A boolean indicating other and self are equal

        EXAMPLES::
            sage: m = MapPermutation([(3,9),(2,1),(7,14,8)])
            sage: m == m
            True
        """
        if isinstance(other, MapPermutation):
            return list(other) == list(self)
        return False

    def left_action_product(self, rperm):
        """
        This function calculate self*perm where * is the composition operation between permutation

        INPUT:

            -rperm: Another MapPermutation

        OUTPUT:
            -A MapPermutation of size max(rperm.size(),self.size()) representing the composition self*rperm

        EXAMPLES::

            sage: m = MapPermutation([(3,9),(2,1),(7,14,8)])
            sage: t = MapPermutation([(1,2)])
            sage: m.left_action_product(t)
            [1, 2, 9, 4, 5, 6, 14, 7, 3, 10, 11, 12, 13, 8]

        ..NOTE::
            O(max(n,t)) ,where n is the size of self and t the size of rperm
        """
        outSize = max(self.size(), rperm.size())

        outList = [self(rperm(i)) for i in range(1, outSize + 1)]

        return MapPermutation(outList)

    def number_of_fixed_points(self):
        """
        OUTPUT:
            The number of fixed point ( we only consider i such that i<=self.size())

        EXAMPLES::

            sage: MapPermutation([(3,9),(2,1),(7,14,8)]).number_of_fixed_points()
            7

        .. NOTE::
            O(n) where n is the size of self
        """
        return np.sum(self._tab == np.arange(1, self.size()+1))

    def right_action_product(self, lperm):
        """
        This function calculate lperm*self where * is the composition operation between permutation

        INPUT:
            -lperm: Another map permutation

        OUTPUT:
            -A MapPermutation of size max(lperm.size(),self.size()) representing the composition lperm*self

        EXAMPLES::

            sage: m = MapPermutation([(3,9),(2,1),(7,14,8)])
            sage: t = MapPermutation([(1,2)])
            sage: m.right_action_product(t)
            [1, 2, 9, 4, 5, 6, 14, 7, 3, 10, 11, 12, 13, 8]

        .. NOTE::
            O(max(n,t)),where n is the size of self and t the size of lperm
        """
        return lperm.left_action_product(self)

    def __mul__(self, b):
        """
        This function calculate lperm*self where * is the composition operation between permutation

        INPUT:
            -b: Another map permutation

        OUTPUT:
            -A MapPermutation of size max(b.size(),self.size()) representing the composition self*b

        EXAMPLES::

            sage: m = MapPermutation([(3,9),(2,1),(7,14,8)])
            sage: t = MapPermutation([(1,2,4)])
            sage: m*t
            [1, 4, 9, 2, 5, 6, 14, 7, 3, 10, 11, 12, 13, 8]

        .. NOTE::
            O(max(n,t)),where n is the size of self and t the size of  b
        """

        return self.left_action_product(b)
