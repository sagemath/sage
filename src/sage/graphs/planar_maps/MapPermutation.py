from sage.all import Permutation
from sage.graphs.planar_maps.MapError import InvalidMapPermutationArgument
import numpy as np


class MapPermutation:

    def __init__(self, lst, trust=False) -> None:
        if isinstance(lst, Permutation):
            self._init_from_permutation(lst, trust=trust)
            return
        try:
            if lst == int(lst) or lst > 0:
                self._init_from_number(lst)
                return
        except:
            pass

        try:
            if type(lst[0]) == type((42,)):
                self._init_from_cycle_list(lst, trust=trust)
                return
            self._init_from_list(lst)
        except:
            raise InvalidMapPermutationArgument()

    def _init_from_cycle_list(self, lst, trust=False):
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
        self._tab = np.arange(1, n+1)
        self._rtab = np.arange(1, n+1)

    def _init_from_permutation(self, perm, trust=False):
        self._init_from_list(list(perm), trust=trust)

    def _init_from_list(self, list, trust=False):
        self._tab = np.array(list)
        self._rtab = np.zeros(len(list), dtype=int)
        self._rtab[self._tab-1] = np.arange(1, len(list)+1)
        if trust:
            return
        if not np.issubdtype(self._tab.dtype, np.integer) or ((self._tab > len(self._tab)) + (self._tab <= 0)).sum() != 0 or len(np.unique(self._tab)) != len(self._tab):
            raise InvalidMapPermutationArgument()

    def size(self):
        """
        Returns: 
            The size of self
        ---
        O(1)
        """
        return len(self._tab)

    def apply(self, i):
        """
        Args:
            i
        ---
        Returns:
            self(i)
        ---
        O(1)
        """
        if i > self.size():
            return i
        return self._tab[i-1]

    def inverseApply(self, i):
        """
        Args:
            i an index
        ----
        Returns:
            j such that self(j) = i
        ----
        O(1)
        """

        if i > self.size():
            return i
        return self._rtab[i-1]

    def __repr__(self) -> str:
        return str(list(self))

    def pretty_repr(self):
        return f"MapPermutation: {self.to_cycles()}"

    def pretty_print(self):
        """
        Print self in a more pretty form
        """
        print(self.pretty_repr())

    def to_cycles(self):
        """
        Returns:
            A list of tuple representing the cycle decomposition of self
        ---
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
        -------
        Returns:
            - The inverse of self
        -------
        O(n)
        where n is the number of element of the permutation
        -------
        """
        return MapPermutation(list(self._rtab))

    def __call__(self, i):
        return self.apply(i)

    def __getitem__(self, id):
        return self(id)

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        if self.index < len(self):
            result = self(self.index + 1)
            self.index += 1
            return result
        else:
            raise StopIteration

    def __len__(self):
        return self.size()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return list(other) == list(self)
        return False

    def left_action_product(self, rperm):
        """
        This function calculate self*perm where * is the composition operation between permutation
        -------
        Args:
            -rperm: Another MapPermutation
        Returns:
            -A MapPermutation of size max(rperm.size(),self.size()) representing the composition self*rperm
        -------
        O(max(n,t))
        where n is the size of self and t the size of rperm
        -------
        """
        outSize = max(self.size(), rperm.size())

        outList = [self(rperm(i)) for i in range(1, outSize + 1)]

        return MapPermutation(outList)

    def number_of_fixed_points(self):
        """
        Returns: the number of fixed point ( we only consider i such that i<=self.size())
        ----
        O(n) where n is the size of self
        """
        return np.sum(self._tab == np.arange(1, self.size()+1))

    def right_action_product(self, lperm):
        """
        This function calculate lperm*self where * is the composition operation between permutation
        -------
        Args:
            -lperm: Another map permutation
        Returns:
            -A MapPermutation of size max(lperm.size(),self.size()) representing the composition lperm*self
        -------
        O(max(n,t))
        where n is the size of self and t the size of lperm
        -------
        """
        return lperm.left_action_product(self)

    def __mul__(self, b):
        return self.left_action_product(b)
