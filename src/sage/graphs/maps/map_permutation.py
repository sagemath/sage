from typing import Any
from sage.all import Permutation
from sage.graphs.maps.map_error import InvalidMapPermutationArgumentError
from sage.rings.integer import Integer
import numpy as np


class MapPermutation:
    """A class to efficiently represent Permutations used in different maps."""

    def __init__(self, lst: int | Permutation | list[int] | list[tuple[int, ...]], trust=False):
        """
        Init the MapPermutation

        INPUT:
        - ``lst`` -- int | Permutation | list[int] | list[tuple[int, ...]] ; an integer (interpreted as the size), or a Sage permutation object, or a list representing the permutation, or a list of tuples representing the cycles of the permutation.
        - ``trust`` -- bool ;  a parameter telling whether or not
        to pass some test to verify if lst represent
        a valid permutation.

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: MapPermutation(3)
            [1, 2, 3]
            sage: MapPermutation([2,3,1])
            [2, 3, 1]
            sage: MapPermutation([(2,3) , (1,4)])
            [4, 3, 2, 1]

        NOTE:
            O(n) where n is the size of the permutation
        """

        def isInt(x) -> bool:
            try:
                _ = int(x)
            except Exception:
                return False
            else:
                return True

        if isinstance(lst, Permutation):
            self._init_from_permutation(lst, trust=trust)
        elif isInt(lst):
            self._init_from_number(int(lst))
        elif type(lst) in (list, tuple) and len(lst) >= 1 and isInt(lst[0]):        # we need to include ints, Sage Integers, numpy ints...
            self._init_from_list(list(map(int, lst)))
        elif type(lst) in (list, tuple) and len(lst) >= 1 and type(lst[0]) in (list, tuple):
            self._init_from_cycle_list([tuple(map(int, t)) for t in lst], trust=trust)
        else:
            raise InvalidMapPermutationArgumentError()

    def _init_from_cycle_list(self, lst: list[tuple[int, ...]], trust=False) -> None:
        """
        Init the permutation from a list of cycles

        INPUT:

            -``lst`` -- list[tuple[int, ...]] ;  a list of tuples of ints representing the cycle of self
            - ``trust`` -- bool ; a boolean indicating if check should be passed

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: lst = [(4,3),(2,1)]
            sage: MapPermutation(lst)._init_from_cycle_list(lst)

        NOTE:
            O(m),Used internally not intended to be used by the user
        """
        cnt = max(map(max, lst))
        self._tab = np.arange(1, cnt+1, dtype=int)
        self._rtab = np.arange(1, cnt+1, dtype=int)
        for e in lst:
            prev = 0
            for i in e:
                if i <= 0:
                    raise InvalidMapPermutationArgumentError()

                if prev != 0:
                    self._rtab[i-1] = prev
                    self._tab[prev-1] = i
                prev = i

            self._tab[e[-1]-1] = e[0]
            self._rtab[e[0]-1] = e[-1]

        if trust:
            return

        if (self._tab == 0).sum() != 0:
            raise InvalidMapPermutationArgumentError()

    def _init_from_number(self, n: int) -> None:
        """
        INPUT:
            - ``n`` -- int ; the size of the identify permutation

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: MapPermutation(3)._init_from_number(3)

        NOTE:
            O(n),Used internally not intended to be used by the user
        """
        self._tab = np.arange(1, n+1)
        self._rtab = np.arange(1, n+1)

    def _init_from_permutation(self, perm: Permutation, trust=False) -> None:
        """
        INPUT:
        - ``perm`` -- Permutation; a Permutation object
        - ``trust`` -- bool; a boolean indicating if check should be
        passed

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: MapPermutation(3)._init_from_number(3)

        NOTE:
            O(n),Used internally not intended to be used by the user
        """

        self._init_from_list(list(perm), trust=trust)

    def _init_from_list(self, list: list[int], trust=False) -> None:
        """
        Init the permutation from a list of index

        INPUT:
            -``lst`` -- List[int] ; a list of index
            - ``trust`` -- bool ; a boolean indicating if check should
            be passed

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: lst = [(4,3),(2,1)]
            sage: MapPermutation(lst)._init_from_cycle_list(lst)

        NOTE:
            O(m),Used internally not intended to be used by the user

        """
        self._tab = np.array(list)
        self._rtab = np.zeros(len(list), dtype=int)
        self._rtab[self._tab-1] = np.arange(1, len(list)+1)
        if trust:
            return
        if not np.issubdtype(self._tab.dtype, np.integer) or ((self._tab > len(self._tab)) + (self._tab <= 0)).sum() != 0 or len(np.unique(self._tab)) != len(self._tab):
            raise InvalidMapPermutationArgumentError()
        
    def to_list(self) -> list[int]:
        """
        Return the permutation self, as a list of integers.
        
        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: MapPermutation([(2,3,4),(6,7)]).to_list()
            [1, 3, 4, 2, 5, 7, 6]
            
        NOTE:
            O(n)
        """
        return [int(self(i)) for i in range(1, self.size() + 1)]

    def size(self) -> int:
        """
        OUTPUT:
            The size of self

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: MapPermutation([(3,32),(6,7)]).size()
            32

        NOTE:
            O(1)
        """
        return len(self._tab)

    def apply(self, i: int) -> int:
        """
        INPUT:
            - ``i`` -- int ; An index

        OUTPUT:
            self(i)

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: int(MapPermutation([(3,22),(22,33),(7,14)]).apply(3))
            22

        NOTE:
            O(1)
        """
        if i > self.size():
            return i
        return self._tab[i-1]

    def inverseApply(self, i: int) -> int:
        """
        INPUT:
            - ``i`` -- int ; An index

        OUTPUT:
            j such that self(j) = i

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: int(MapPermutation([(3,22),(22,33),(7,14,8)]).inverseApply(7))
            8

        NOTE:
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

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: MapPermutation([(3,9),(2,1),(7,14,8)])
            [2, 1, 9, 4, 5, 6, 14, 7, 3, 10, 11, 12, 13, 8]

        """
        return "[" + ", ".join(map(str, self)) + "]"  # Permet d'afficher de la même manière les int et les np.int64
        # return str(list(self))

    def pretty_repr(self) -> str:
        """
        OUTPUT:

            A more pretty string representation of self

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: MapPermutation([(3,9),(2,1),(7,14,8)]).pretty_repr()
            'MapPermutation: [(1, 2), (3, 9), (4,), (5,), (6,), (7, 14, 8), (10,), (11,), (12,), (13,)]'

        """
        return f"MapPermutation: {self.to_cycles()}"

    def pretty_print(self) -> None:
        """
        Print self in a more pretty form

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: MapPermutation([(3,9),(2,1),(7,14,8)]).pretty_print()
            MapPermutation: [(1, 2), (3, 9), (4,), (5,), (6,), (7, 14, 8), (10,), (11,), (12,), (13,)]

        """
        print(self.pretty_repr())

    def to_cycles(self) -> list[tuple[int, ...]]:
        """
        OUTPUT:
            A list of tuple representing the cycle decomposition of self

        NOTE:
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
                cycle.append(int(curIndex))
                curIndex = self(curIndex)
            cycles.append(tuple(cycle))
        return cycles

    def inverse(self) -> "MapPermutation":
        """
        This function calculate  the inverse of self

        OUTPUT:
            - The inverse of self

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: MapPermutation([(3,9),(2,1),(7,14,8)]).inverse()
            [2, 1, 9, 4, 5, 6, 8, 14, 3, 10, 11, 12, 13, 7]

        NOTE:
            O(n),where n is the number of element of the permutation
        """
        return MapPermutation(list(self._rtab))

    def __call__(self, i: int) -> int:
        """
        OUTPUT:

            self(i)

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: int(MapPermutation([(3,9),(2,1),(7,14,8)])(3))
            9

        """
        return self.apply(i)

    def __getitem__(self, id: int) -> int:
        """
        INPUT:

            -``id`` -- int

        OUTPUT:
            id th item

        EXAMPLES::
            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: int(MapPermutation([(3,9),(2,1),(7,14,8)]).__getitem__(3))
            9
        """
        return self(id)

    def __iter__(self) -> "MapPermutation":
        """

        OUTPUT:
            self, with self.index to 0

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: MapPermutation([(3,9),(2,1),(7,14,8)]).__iter__()
            [2, 1, 9, 4, 5, 6, 14, 7, 3, 10, 11, 12, 13, 8]

        """
        self.index = 0
        return self

    def __next__(self) -> int:
        """

        OUTPUT:
            The next value in during an iteration or raise StopIteration

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: int(MapPermutation([(3,9),(2,1),(7,14,8)]).__iter__().__next__())
            2

        """

        if self.index < len(self):
            result = self(self.index + 1)
            self.index += 1
            return result
        else:
            raise StopIteration

    def __len__(self) -> int:
        """

        OUTPUT:
            The size of self

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: len(MapPermutation([(3,9),(2,1),(7,14,8)]))
            14

        """
        return self.size()

    def __eq__(self, other: Any) -> bool:
        """
        INPUT:
            - ``other`` -- MapPermutation

        OUTPUT:
            A boolean indicating other and self are equal

        EXAMPLES::
            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: m = MapPermutation([(3,9),(2,1),(7,14,8)])
            sage: m == m
            True
        """
        if isinstance(other, MapPermutation):
            return list(other) == list(self)
        return False

    def left_action_product(self, rperm: "MapPermutation") -> "MapPermutation":
        """
        This function calculate self*perm where * is the composition operation between permutation

        INPUT:

            - ``rperm`` -- MapPermutation

        OUTPUT:
            -A MapPermutation of size max(rperm.size(),self.size()) representing the composition self*rperm

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: m = MapPermutation([(3,9),(2,1),(7,14,8)])
            sage: t = MapPermutation([(1,2)])
            sage: m.left_action_product(t)
            [1, 2, 9, 4, 5, 6, 14, 7, 3, 10, 11, 12, 13, 8]

        NOTE:
            O(max(n,t)) ,where n is the size of self and t the size of rperm
        """
        outSize = max(self.size(), rperm.size())

        outList = [self(rperm(i)) for i in range(1, outSize + 1)]

        return MapPermutation(outList)

    def number_of_fixed_points(self) -> int:
        """
        OUTPUT:
            The number of fixed point ( we only consider i such that i<=self.size())

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: int(MapPermutation([(3,9),(2,1),(7,14,8)]).number_of_fixed_points())
            7

        NOTE:
            O(n) where n is the size of self
        """
        return np.sum(self._tab == np.arange(1, self.size()+1))

    def right_action_product(self, lperm: "MapPermutation") -> "MapPermutation":
        """
        This function calculate lperm*self where * is the composition operation between permutation

        INPUT:
            -``lperm`` -- MapPermutation

        OUTPUT:
            -A MapPermutation of size max(lperm.size(),self.size()) representing the composition lperm*self

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: m = MapPermutation([(3,9),(2,1),(7,14,8)])
            sage: t = MapPermutation([(1,2)])
            sage: m.right_action_product(t)
            [1, 2, 9, 4, 5, 6, 14, 7, 3, 10, 11, 12, 13, 8]

        NOTE:
            O(max(n,t)),where n is the size of self and t the size of lperm
        """
        return lperm.left_action_product(self)

    def __mul__(self, b: "MapPermutation") -> "MapPermutation":
        """
        This function calculate lperm*self where * is the composition operation between permutation

        INPUT:
            - ``b`` -- MapPermutation

        OUTPUT:
            -A MapPermutation of size max(b.size(),self.size()) representing the composition self*b

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: m = MapPermutation([(3,9),(2,1),(7,14,8)])
            sage: t = MapPermutation([(1,2,4)])
            sage: m*t
            [1, 4, 9, 2, 5, 6, 14, 7, 3, 10, 11, 12, 13, 8]

        NOTE:
            O(max(n,t)),where n is the size of self and t the size of  b
        """

        return self.left_action_product(b)
