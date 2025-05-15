from sage.graphs.maps.cyclic_chained_list import CyclicChainedList
from sage.all import Permutation
from sage.graphs.maps.map_permutation import MapPermutation
from sage.graphs.maps.map_error import NotImplementedError


class PrimitiveRotatingPermutation(MapPermutation):

    """

    This class represent a more primitive version of rotating permutation useful in PrimitiveMutableLabelledMap

    WARNING: We take as a convention for this class that if i is bigger than the size of self,
    then self(i) = i.
      """

    def __init__(self, lst):
        """
        This function initiate the primitive rotating permutation, lst can be  a Permutation or a list of int or list of tuple representing the cycle of
        the permutation or a MapPermutation or an integer representing the size of the permutation(in this case self will represent the identity permutation of size lst).

        INPUT:

        - ``lst`` -- List[int] | List[Tuples] | int | Permutation | MapPermutation ; a list representing the permutation or a list of tuples representing 
        the cycle of the permutationor an integer representing the size 
        of the permutation(In This case it will return the identify of size lst) 
        or a Permutation or a MapPermutation.


        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: PrimitiveRotatingPermutation([(1,3,4),(7,8,10),(9,6)])
            [3, 2, 4, 1, 5, 9, 8, 10, 6, 7]

        NOTE:
            O(n) where n is the size of the permutation

        """
        if isinstance(lst, Permutation) or isinstance(lst, MapPermutation):
            self.__init__(list(lst))
            return

        # Important will stay
        # Associate to a node its image
        self._permCycle = {}
        # Size of the permutation
        self._n = 0
        # Number of cycles in the permutation
        self._numCycles = 0

        # Number of fixed point in the permutation
        self._numberOfFixedPoint = 0

        try:
            if lst == int(lst) and lst > 0:
                # If lst is an integer we just set our permutation to be the
                # identity
                self._n = lst
                self._numCycles = self._n
                return
        except BaseException:
            pass

        mx = 0
        seen = []
        try:
            # We're directly using the cycle representation to initialise the
            # permutation
            if isinstance(lst[0], type((42,))):
                for l in lst:
                    for i in l:
                        mx = max(i, mx)
                        if i != int(i) or i <= 0:
                            raise ValueError(
                                f"Invalid argument: {i} isn't a strictly positive integer in the list given")
                seen = [False for i in range(mx + 1)]
                cnt = 0
                for l in lst:
                    k = 0
                    while k < len(l):
                        i = l[k]
                        if seen[i]:
                            raise ValueError(
                                f"Invalid argument: {i} appears at least two times in list given it cannot be a permutation.")
                        seen[i] = True
                        k += 1
                        while k < len(l) and l[k] == i:
                            k += 1
                            continue
                        cnt += 1
                self._numCycles += mx - cnt
                self._numberOfFixedPoint += mx - cnt
                for l in lst:
                    prevNode = None
                    for i in l:
                        newNode = CyclicChainedList(i)
                        self._permCycle[i] = newNode
                        if prevNode is not None:
                            prevNode.insertAfter(newNode)
                        prevNode = newNode
                    self._numberOfFixedPoint += len(l) == 1
                    self._numCycles += 1

            else:
                mx = len(lst)
                for i in lst:
                    if i != int(i) or i <= 0:
                        raise ValueError(
                            f"Invalid argument : {i} isn't a strictly positive integer in the list given")
                    if i > len(lst):
                        raise ValueError(
                            f"{i} is bigger than the size of the given list")
                seen = [False for i in range(mx + 1)]
                for i in lst:
                    if seen[i]:
                        raise ValueError(
                            f"Invalid argument: {i} appears at least two time in the list given it cannot be a permutation..")
                    seen[i] = True

                seen = [False for i in range(mx + 1)]
                for i in range(1, mx + 1):
                    if seen[i]:
                        continue
                    prevNode = None
                    curElement = i
                    cnt = 0
                    while not seen[curElement]:
                        newNode = CyclicChainedList(curElement)
                        self._permCycle[curElement] = newNode
                        if prevNode is not None:
                            prevNode.insertAfter(newNode)
                        prevNode = newNode
                        seen[curElement] = True
                        curElement = lst[curElement - 1]
                        cnt += 1
                    self._numCycles += 1
                    self._numberOfFixedPoint += cnt == 1
        except ValueError as e:
            raise e
        except BaseException:
            raise ValueError("Invalid argument: The argument given must be Permutation or MapPermutation or a non empty list of integers representing the permutation or a non empty list of tuples representing the cycles of the permutations or a positive integer.")
        self._n = mx

    # OK

    def size(self):
        """
        OUTPUT:
            The size of the permutation

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm.size()
            8

        NOTE:
            O(1)
        """

        return self._n

    # OK
    def deleteLastKIndex(self, k):
        """
        This function will delete the last k index from self

        INPUT:
        -``k`` -- int ; the number of node to delete

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4),(7,8,10),(9,6)])
            sage: rperm
            [3, 2, 4, 1, 5, 9, 8, 10, 6, 7]
            sage: rperm.deleteLastKIndex(3)
            sage: rperm
            [3, 2, 4, 1, 5, 6, 7]

        NOTE:
            O(k) 
        """

        if k > self.size():
            raise ValueError(
                f"Cannot delete {k} last element in a PrimitiveRotatingPermutation of size {self.size()}")
        for _ in range(k):
            self.delete(self._n)

    # OK

    def delete(self, index):
        """
        This will delete index of the corresponding cycle note that after this operation if we note the original
        size of self as n, the which contained index will count one less element,
        self will be of size n-1 and if n != index the element numbered n will relabeled as index.
        For instance if self is the permutation(1, 2, 3)(4, 5) and we delete 2 it will become(1, 3)(4, 2),
        If n = 1 an error or index is not a strictly positive integer <= n an error will be raised.

        INPUT:
        - ``index`` -- int ;

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4),(7,8,10),(9,6)])
            sage: rperm
            [3, 2, 4, 1, 5, 9, 8, 10, 6, 7]
            sage: rperm.delete(10)
            sage: rperm
            [3, 2, 4, 1, 5, 9, 8, 7, 6]

        NOTE: 
            O(1),index must be an strictly positive integer and self.size() >= 2 otherwise an error will be raised
        """

        if self.size() == 1:
            raise ValueError(
                "Cannot delete an element from a Permutation of size 1")
        if index != int(index) or index <= 0 or index > self.size():
            raise ValueError(
                "{index} isn't a strictly positive integer <= self.size()")

        u = self(index)
        v = self(u)

        nPrev = self.size()
        node = self.getNode(index)

        node.remove()

        if v == index and u != index:
            self._numberOfFixedPoint += 1

        if u == index:
            self._numberOfFixedPoint -= 1
            self._numCycles -= 1

        self._n -= 1

        self._permCycle.pop(index)

        if nPrev != index:
            try:
                self._permCycle[nPrev].val = index

                self._permCycle[index] = self._permCycle[nPrev]

                self._permCycle.pop(nPrev)
            except BaseException:
                pass

    # OK
    def inverseApply(self, i):
        """
        This function apply  the inverse self on i, we take as a convention i if i is an integer > self.size(), self.inverseApply(i) = i

        INPUT:
        - ``i`` -- int

        OUTPUT:
            j such that self(j) = i


        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4),(7,8,10),(9,6)])
            sage: rperm.inverseApply(10)
            8

        NOTE:
            O(1)
        """

        if i != int(i) or i <= 0:
            raise ValueError("{i} isn't a positive integer")
        try:
            return self._permCycle[i].prev.val
        except BaseException:
            return i

    # OK

    def swapIndex(self, index, otherIndex):
        """
        This function swap the index role in the permutation

        INPUT:
        - ``index`` -- int ;<= self.size()
        - ``otherIndex`` -- int ; <= self.size()

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4),(7,8,10),(9,6)])
            sage: rperm
            [3, 2, 4, 1, 5, 9, 8, 10, 6, 7]
            sage: rperm.swapIndex(1,10)
            sage: rperm
            [7, 2, 4, 10, 5, 9, 8, 1, 6, 3]

        NOTE:
            O(1)
        """

        nodeIndex = self.getNode(index)
        nodeOther = self.getNode(otherIndex)
        self._permCycle[otherIndex] = nodeIndex
        self._permCycle[index] = nodeOther
        self._permCycle[otherIndex].val = otherIndex
        self._permCycle[index].val = index

    def cutDelete(self, startIndex, endIndex):
        """
        This will cut the cycle in two part startIndex...endIndex and the rest , and than will delete startIndex and endIndex.

        INPUT: 
        - ``startIndex`` -- int ; on the same cycle as endIndex such that (startIndex,endIndex)=(self.size(),self.size()-1)  
        - ``endIndex`` -- int ;  on the same cycle as startIndex such that (startIndex,endIndex)=(self.size(),self.size()-1)

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4),(2,7,11,8,10),(9,6)])
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4), (2, 7, 11, 8, 10), (5,), (6, 9)]
            sage: rperm.cutDelete(10,11)
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4), (2, 7), (5,), (6, 9), (8,)]

        NOTE:
            O(1)
        """

        assert startIndex != endIndex

        tempNewIndex = self._n+1
        tempNewIndexOther = self._n+2
        self.cutAdd(startIndex, endIndex, tempNewIndex, tempNewIndexOther)

        self.deleteLastKIndex(4)

    # OK

    def cutAdd(self, startIndex, endIndex, newIndexStart, newIndexEnd):
        """
        This implement a special operation.In a nutshell it cut a cycle and add two index in each cycle,
        let denote A = startIndex, B = endIndex, C = newIndexStart, D = newIndexEnd and say the cycle is of the form F -> A -> S -> .. -> T -> B -> R -> ... -> F
        than the situation will be the following after a call to this function, A -> S -> ... -> T -> D -> A and F -> C -> B -> R -> ... -> F

        INPUT:
        - ``startIndex`` -- int ; on same cycle as ``endIndex``
        - ``endIndex ``  -- int ; on same cycle as ``startIndex``
        - ``newIndexStart`` -- int; such that {newIndexEnd, newIndexStart} = {n+1, n+2} 
        - ``newIndexEnd`` -- int : such that {newIndexEnd, newIndexStart} = {n+1, n+2} 

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4),(2,7,11,8,10),(9,6)])
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4), (2, 7, 11, 8, 10), (5,), (6, 9)]
            sage: rperm.cutAdd(1,4,12,13)
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 13), (2, 7, 11, 8, 10), (4, 12), (5,), (6, 9)]  

        NOTE:
            O(1)
        """

        if newIndexEnd == newIndexStart:
            raise ValueError(
                f"{newIndexEnd} and {newIndexStart} must be different")
        if newIndexStart <= self.size() or newIndexStart > self.size() + 2:
            raise ValueError(
                f"{newIndexStart} must be  >{self.size()} and <= {self.size() + 2}")
        if newIndexEnd <= self.size() or newIndexEnd > self.size() + 2:
            raise ValueError(
                f"{newIndexEnd} must be  >{self.size()} and <= {self.size() + 2}")
        if startIndex == endIndex:
            self.addBefore(startIndex)
            self.stretch(1)
            return

        nodeStartIndex = self.getNode(startIndex)
        nodeEndIndex = self.getNode(endIndex)

        # Updating scalar attribute accordingly
        self._numCycles += 1
        self._n += 2

        nodeNewIndexStart = self.getNode(newIndexStart)
        nodeNewIndexEnd = self.getNode(newIndexEnd)

        # NodeNewIndexStart processing
        tmpNode = self.getNode(self.inverseApply(startIndex))
        nodeNewIndexStart.prev = tmpNode
        tmpNode.nxt = nodeNewIndexStart

        nodeNewIndexStart.nxt = nodeEndIndex

        # NodeNewIndexEnd processing
        tmpNode = self.getNode(self.inverseApply(endIndex))
        nodeNewIndexEnd.prev = tmpNode
        tmpNode.nxt = nodeNewIndexEnd

        nodeNewIndexEnd.nxt = nodeStartIndex

        nodeEndIndex.prev = nodeNewIndexStart
        nodeStartIndex.prev = nodeNewIndexEnd

    # OK

    def labelToTheEnd(self, listIndexes):
        """
        This is a helper function  it just move all of the element in listIndexes to the last indices

        INPUT:
        - ``listIndexes`` -- List[int]

        OUTPUT:

            A map giving a correspondence between the old index and the new
            if it was changed.

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4),(2,7,11,8,10),(9,6)])
            sage: rperm
            [3, 7, 4, 1, 5, 9, 11, 10, 6, 2, 8]
            sage: rperm.labelToTheEnd([3,2,11])
            {2: 9, 3: 10}
            sage: rperm
            [10, 6, 9, 1, 5, 2, 11, 3, 7, 4, 8]

        NOTE:
            O(len(listIndexes))

        """

        for index in listIndexes:
            if index != int(index) or index <= 0 or index > self.size():
                raise ValueError(
                    f"In labelToTheEnd : {index} isn't a strictly positive integer <= {self.size()}")

        indexMap = set()
        for index in listIndexes:
            indexMap.add(index)

        indexCandidate = set()
        for j in range(len(indexMap)):
            indexCandidate.add(self.size() - j)

        for index in list(indexCandidate):
            if index in indexMap:
                indexCandidate.remove(index)
                indexMap.remove(index)

        corresOut = {}
        for index in list(indexMap):
            if index not in indexMap:
                continue
            corresIndex = indexCandidate.pop()
            corresOut[index] = corresIndex
            indexMap.remove(index)
            self.swapIndex(index, corresIndex)
        return corresOut
    # OK

    def bruteAddCycles(self, cycles):
        """
        Another helper function that add cyclein cycles, this one assumed is more dangerous than addCycles
        cause it assumed that the cycles are well formed thus the term brute

        INPUT:
        - ``cycles`` -- List[Tuple] ; list of cycles as tuple

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4)])
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4), (2,)]
            sage: rperm.stretch(3)
            sage: rperm.bruteAddCycles([(5,6,2)] )
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4), (2, 5, 6), (7,)]

        NOTE:
            O(len(cycles))
        """

        for c in cycles:
            for i in range(len(c) - 1):
                self.addAfterGeneral(c[i], c[i + 1])

    # OK
    def addCycles(self, cycles):
        """
        Another helper function it will raise an error if element of the cycles
        are not > self.size() and <= self.size()+len(cycles), the cycle must be well formed

        INPUT:
        - ``cycles`` -- List[Tuple] ; list of cycles

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4)])
            sage: rperm.addCycles([(5,6)])
            sage: rperm
            [3, 2, 4, 1, 6, 5]

        NOTE:
            O(len(cycles))
        """

        testSet = set()
        N = 0
        for c in cycles:
            N += len(c)

        for c in cycles:
            for e in c:
                testSet.add(e)
                if e <= self.size() or e <= 0 or e != int(e) or e > self.size() + N:
                    raise ValueError("{cycles} isn't valid")
        if len(testSet) != N:
            raise ValueError("{cycles} isn't valid")

        self.stretch(N)

        self.bruteAddCycles(cycles)

    # OK
    def isValidIndex(self, index):
        """
        Check if index is a integer > 0 and <=self.size()
        otherwise raise an Error

        INPUT:
        - ``index`` -- int

        EXAMPLES:: 

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4)])
            sage: try :
            ....:     rperm.isValidIndex(100)
            ....: except:
            ....:     pass
            ....:

        NOTE:
            O(1)
        """

        if index <= 0 or index != int(index) or index > self.size():
            raise ValueError(f"{index} isn't valid")
    # OK

    def addAfterGeneral(self, index, otherIndex):
        """
        This is a more general version of addAfter it only assumed that otherIndex is a fixed point
        and will add it after index in its cycle
        INPUT:
        - ``index`` -- int
        -- ``otherIndex`` -- int

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4)])
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4), (2,)]
            sage: rperm.addAfterGeneral(4,2)
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4, 2)]

        NOTE:
            O(1)
        """

        self.isValidIndex(index)
        self.isValidIndex(otherIndex)
        if index == otherIndex:
            return
        if not self(otherIndex) == otherIndex:
            raise ValueError(
                f"Can only add after fixed point {otherIndex} isn't one")

        self._numberOfFixedPoint -= self(index) == index

        node = self.getNode(index)

        newNode = self.getNode(otherIndex)

        self._numberOfFixedPoint -= 1
        self._numCycles -= 1

        node.insertAfter(newNode)

    # OK
    def addBeforeGeneral(self, index, otherIndex):
        """
        More general version of addBeforeit only assumed that otherIndex is a fixed point
        and will add it before index in its cycle
        INPUT:
        - ``index``  -- int
        - ``otherIndex`` -- int

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4)])
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4), (2,)]
            sage: rperm.addBeforeGeneral(3,2)
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 2, 3, 4)]

        NOTE:
            O(1)
        """

        self.isValidIndex(index)
        indexPrev = self.inverseApply(index)
        self.addAfterGeneral(indexPrev, otherIndex)

    # OK
    def mergeDelete(self, index, otherIndex):
        """
        Assuming that index and otherIndex are not in the same cycle it will do the
        following first index and otherIndex will be sent to self.size() self.size()-1 they will be deleted and given
        that before we add: U -> ... -> V -> index -> R -> U and F -> ... -> T -> otherIndex -> Q -> F, we will have after
        U -> ... -> V -> Q -> F -> ... -> T -> R -> U.

        index, otherIndex are two node not on the same cycle, WARNING: if it isn't the case there is no guarantee and no error will be raised

        INPUT:
        - ``index`` -- int
        - ``otherIndex`` -- int

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4), (2, 5, 7, 8), (6,)]
            sage: rperm.mergeDelete(3,7)
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 2, 5, 4), (6,)]

        NOTE:
            O(1)
        """

        backUpNumberOfFixedPoint = self.number_of_fixed_points()

        self.labelToTheEnd([index, otherIndex])

        if self(
                self._n) == self._n or self(self._n - 1) == self._n-1:
            self.deleteLastKIndex(2)
            return

        beforeIndex = self.inverseApply(self._n)
        afterIndex = self.apply(self._n - 1)

        self.deleteLastKIndex(2)

        nodeBefore = self.getNode(beforeIndex)
        nodeAfter = self.getNode(afterIndex)

        if nodeBefore.nxt == nodeBefore:
            if nodeAfter.nxt == nodeAfter:
                nodeBefore.nxt = nodeAfter
                nodeBefore.prev = nodeAfter
                nodeAfter.nxt = nodeBefore
                nodeAfter.prev = nodeBefore
            else:
                tmpNode = nodeAfter.prev
                tmpNode.nxt = nodeBefore
                nodeBefore.prev = tmpNode
                nodeAfter.prev = nodeBefore
                nodeBefore.nxt = nodeAfter

        else:
            if nodeAfter.nxt == nodeAfter:
                tmpNode = nodeBefore.nxt
                tmpNode.prev = nodeAfter
                nodeAfter.nxt = tmpNode
                nodeAfter.prev = nodeBefore
                nodeBefore.nxt = nodeAfter
            else:
                tmpNodeBefore = nodeBefore.nxt
                tmpNodeAfter = nodeAfter.prev
                tmpNodeBefore.prev = tmpNodeAfter
                tmpNodeAfter.nxt = tmpNodeBefore

                nodeBefore.nxt = nodeAfter
                nodeAfter.prev = nodeBefore

        self._numCycles -= 1
        self._numberOfFixedPoint = backUpNumberOfFixedPoint

    # OK
    def getNode(self, index):
        """
        This function will return the node associated to index
        and if it doesn't exit it will create one note that if index > self.size()
        it will raise an error.

        INPUT:
        - ``index`` -- int

        OUTPUT:
            The node associated to index

        EXAMPLES::
            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm.getNode(3) != rperm.getNode(4)
            True

        NOTE:
            O(1)
        """

        if index != int(index) or index <= 0 or index > self.size():
            raise ValueError(
                "{index} isn't a strictly positive integer <= self.size()")

        try:
            node = self._permCycle[index]
        except BaseException:
            node = CyclicChainedList(index)
            self._permCycle[index] = node
        return node

    # OK

    def stretch(self, m):
        """
        This function will increase the size of the permutation by m,all the new index will
        be fixed point

        INPUT:
        - ``m`` -- int; ``m``  >= 0 

        EXAMPLES:: 

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm
            [3, 5, 4, 1, 7, 6, 8, 2]
            sage: rperm.stretch(5)
            sage: rperm
            [3, 5, 4, 1, 7, 6, 8, 2, 9, 10, 11, 12, 13] 

        NOTE:
            O(1)
        """

        self._n += m
        self._numberOfFixedPoint += m
        self._numCycles += m
    # OK

    def addAfter(self, index):
        """
        Let denote n=self.size() given that  n>=index>=1, this will increase the size of self by one and add
        the new element n+1 on the cycle of index after index.You should note that if index>self.size() this will raise an error.

        INPUT:
        - ``index`` -- int

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4), (2, 5, 7, 8), (6,)]
            sage: rperm.addAfter(6)
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4), (2, 5, 7, 8), (6, 9)]

        NOTE:
            O(1)
        """

        if index != int(index) or index <= 0 or index > self.size():
            raise ValueError(
                f"{index} isn't a strictly positive integer <= {self.size()}")

        self._numberOfFixedPoint -= self(index) == index
        nPrev = self.size()

        self.stretch(1)

        node = self.getNode(index)

        newNode = self.getNode(nPrev + 1)

        self._numberOfFixedPoint -= 1
        self._numCycles -= 1

        node.insertAfter(newNode)
    # OK

    def addBefore(self, index):
        """
        Let denote n=self.size() given that  n>=index>=1, this will increase the size of self by one and add
        the new element n+1 on the cycle of index before index.You should note that if index>self.size() this will raise an error.

        INPUT:
        - ``index`` -- int

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4), (2, 5, 7, 8), (6,)]
            sage: rperm.addBefore(3)
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 9, 3, 4), (2, 5, 7, 8), (6,)]

        NOTE:
            O(1)
        """

        if index != int(index) or index <= 0 or index > self.size():
            raise ValueError(
                "{index} isn't a strictly positive integer <= self.size()")

        node = self.getNode(index)
        prevIndex = node.prev.val

        self.addAfter(prevIndex)

    def numberInCycle(self, index):
        """
        Not implemented for PrimitiveRotatingPermutation

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: try:
            ....:     rperm.numberInCycle(1)
            ....: except:
            ....:     print("OK")
            ....:
            OK


        """
        raise NotImplementedError(self)

    # OK
    def numberOfCycles(self):
        """
        OUTPUT: 
            the number of cycle of self

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4), (2, 5, 7, 8), (6,)]
            sage: rperm.numberOfCycles()
            3

        NOTE:
            O(1)
        """

        return self._numCycles

    # OK
    def sameCycle(self, i, j):
        """
        Not implemented for PrimitiveRotatingPermutation

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: try:
            ....:     rperm.sameCycle(1,2)
            ....: except:
            ....:     print("OK")
            ....:
            OK

        """

        raise NotImplementedError(self)

    # OK

    def __repr__(self):
        """

        OUTPUT:
        A string representation of self

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm
            [3, 5, 4, 1, 7, 6, 8, 2]

        """
        return "[" + ", ".join(map(str, self)) + "]"  # Permet d'afficher de la même manière les int et les np.int64
        # return str(list(self))

    # OK

    def pretty_repr(self):
        """
        OUTPUT:
            Return a string representation of self in a more pretty form

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm.pretty_repr()
            'Primitive Rotating permutation: [(1, 3, 4), (2, 5, 7, 8), (6,)]'
        """

        return f"Primitive Rotating permutation: {self.to_cycles()}"

    # OK
    def pretty_print(self):
        """
        Print self in a more pretty form

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm.pretty_print()
            Primitive Rotating permutation: [(1, 3, 4), (2, 5, 7, 8), (6,)] 

        """

        print(self.pretty_repr())

    # OK
    def to_cycles(self):
        """
        This method calculate a list of tuple representing the cycle of self
        OUTPUT:
            - lst a list of tuples representing the cycles of self given in increasing order of their minimum elements

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm.to_cycles()
            [(1, 3, 4), (2, 5, 7, 8), (6,)]

        NOTE:
            O(n),where n is the number of element of self
        """

        seen = [False for i in range(self.size() + 1)]
        cycles = []
        for i in range(1, self.size() + 1):
            if seen[i]:
                continue
            try:
                node = self._permCycle[i]
                cycle = node.getValList()
                cycles.append(tuple(map(int, cycle)))
                for j in cycle:
                    seen[j] = True
            except BaseException:
                cycles.append((int(i),))

        return cycles

    # OK
    def inverse(self):
        """
        This function calculate  the inverse of self
        OUTPUT:
            - The inverse of self

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm.inverse()*rperm
            [1, 2, 3, 4, 5, 6, 7, 8]

        NOTE:
            O(n),where n is the number of element of the permutation
        """

        cycles = self.to_cycles()
        return MapPermutation([tuple(reversed(e)) for e in cycles])

    # OK
    def apply(self, i):
        """
        This function apply self on i , we take as a convention i if i is an integer > self.size() , self.apply(i) = i
        INPUT:
        - ``i`` -- int

        OUTPUT:
            self(i)

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm.apply(7)
            8 

        NOTE:
            O(1)
        """

        if i != int(i) or i <= 0:
            raise ValueError("{i} isn't a positive integer")
        try:
            return self._permCycle[i].nxt.val
        except BaseException:
            return i

    # OK

    def number_of_fixed_points(self):
        """
        OUTPUT: 
            the number of fixed point ( we only consider i such that i<=self.size())

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm.number_of_fixed_points()
            1

        NOTE:
            O(1)
        """
        return self._numberOfFixedPoint

    def __eq__(self, other):
        """
        INPUT:
        - ``other`` -- MapPermutation 

        OUTPUT:
            A boolean indicating if self and other are equal

        EXAMPLES::

            sage: from sage.graphs.maps.primitive_rotating_permutation import PrimitiveRotatingPermutation
            sage: rperm = PrimitiveRotatingPermutation([(1,3,4), (7,8,2,5)])
            sage: rperm == rperm
            True 

        """

        if isinstance(other, MapPermutation):
            return list(other) == list(self)
        return False
