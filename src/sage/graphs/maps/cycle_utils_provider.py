from sage.graphs.maps.splay_tree import SplayTree, SplayNode


class CycleUtilsProvider:
    """
    This class is an abstraction of some utils used by RotatingPermutation
    It is mainly use for abstracting some operations
    """

    def __init__(self, cycles):
        """
        Init the CycleUtilsProvider

        EXAMPLES:
            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])

        .. NOTE::
            O(nlog(n)) where n is the sum of size of the cycles
        """
        # It should contain only non fixed point index as key
        self.nodeMap = {}

        # For every cycle we will create a splay tree for it
        # Note that the value contained in the tree must be integer
        # For methods to work
        for c in cycles:
            # If it is a fixed point continue
            if len(c) == 1:
                continue
            # Else create a splay for it
            splayTree = SplayTree(range(len(c)))
            for i in range(len(c)):
                e = c[i]
                self.nodeMap[e] = splayTree.getNode(i)
                self.nodeMap[e].index = e
    # OK
    # Return the number of element in the same cycle as index
    # O(log(n))

    def numberInCycle(self, index):
        """
        INPUT:
            index

        OUTPUT:
            The number of element in the same cycle as index

        EXAMPLES::

            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider.numberInCycle(5)
            2

        .. NOTE::
            O(log(n))
        """
        if self.isFixedPoint(index):
            return 1
        return self.getSplayTree(index).size()
    # OK
    # Return whether i and j are in the same cycle
    # O(log(n))

    def sameCycle(self, i, j):
        """
        INPUT:
            i,j indexes

        OUTPUT:
            A boolean indicating if i and j are on the same cycle


        EXAMPLES::

            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider.sameCycle(5,6)
            False
            sage: provider.sameCycle(7,9)
            True

        .. NOTE::
            O(log(n))
        """
        if self.isFixedPoint(i) or self.isFixedPoint(j):
            return False
        rootI = self.nodeMap[i].getRoot()
        rootJ = self.nodeMap[j].getRoot()
        self.nodeMap[i].SafeSplay()
        return rootI == rootJ

    # If otherIndex is a fixed node it will add it after index
    # Otherwise it will raise an Error
    # O(log(n))
    def addAfter(self, index, otherIndex):
        """
        Add otherIndex in the cycle of index after index

        INPUT:
            index,otherIndex , otherIndex != index and otherIndex is a fixed point

        EXAMPLES::

            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider.addAfter(1,33)
            sage: provider.getCycleList(1)
            [1, 33, 5]

        .. NOTE::
            O(log(n))
        """
        assert index != otherIndex
        if not self.isFixedPoint(otherIndex):
            raise ValueError(f"{otherIndex} isn't a fixed point ")
        splayTree = None
        value = 0
        if index not in self.nodeMap:
            splayTree = SplayTree()
            splayTree.insert(index)
            self.nodeMap[index] = splayTree.getNode(index)
            value = index
        else:
            splayTree = self.nodeMap[index].getSplayTree()
            value = splayTree.root.value + splayTree.root.offset

        leftSplayTree, rightSplayTree = splayTree.split(value + 1 / 2)
        leftSplayTree.insert(value + 1)

        self.nodeMap[otherIndex] = leftSplayTree.getNode(value + 1)

        leftSplayTree.shift(-1)
        rightSplayTree.merge(leftSplayTree)

        self.nodeMap[otherIndex].index = otherIndex
        self.nodeMap[index].index = index

    # Return a boolean indicating if two index in  listIndexes
    # are in the sameCycle efficiently (here O( len(listIndexes)log(n)))
    def checkTwoInTheSameCycle(self, listIndexes):
        """
        INPUT:
            listIndexes a list of indexes

        OUTPUT:
            A boolean indicating if there are two indexes in listIndexes in the same cycle

        EXAMPLES::
            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider.checkTwoInTheSameCycle([33,1,7])
            False
            sage: provider.checkTwoInTheSameCycle([8,1,7])
            True 

        .. NOTE::
            O(len(listIndexes)log(n))
        """
        mapIndexes = set(listIndexes)
        for e in list(mapIndexes):
            if self.isFixedPoint(e):
                try:
                    mapIndexes.remove(e)
                except:
                    pass

        splayTreeMap = set()
        for e in mapIndexes:

            splayTree = self.getSplayTree(e)

            if splayTree in splayTreeMap:
                return True

            splayTreeMap.add(splayTree)

        return False

    # If otherIndex is a fixed node it will add  it before index
    # O(log(n))

    def addBefore(self, index, otherIndex):
        """
        Add otherIndex in the cycle of index before index
        INPUT:
            index,otherIndex , otherIndex != index and otherIndex is a fixed point

        EXAMPLES::
            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider.addBefore(1,33)
            sage: provider.getCycleList(1)
            [33, 1, 5]

        .. NOTE::
            O(log(n))
        """

        assert index != otherIndex
        if not self.isFixedPoint(otherIndex):
            raise ValueError(f"{otherIndex} isn't a fixed point ")

        self._safeIndex(index)

        self.nodeMap[otherIndex] = self.nodeMap[index]

        self.nodeMap[otherIndex].index = otherIndex

        self.nodeMap.pop(index)

        self.addAfter(otherIndex, index)

    # O(log(n))
    # We don't have node for fixed point index
    # Hence it can be dangerous to access the node map of an index when it is a fixed point
    # This will create a node if it doesn't exist for index and the corresponding splaytree
    # Useful during operations to not have to have a different logic for fixed point that won't be anymore after the operations
    def _safeIndex(self, index):
        """
        Used internaly,because we don't have node for fixed point index , it may create a different logic for them
        this will temporarly create a node during this operations , the use must be careful that after having used 
        this temporary node , that it is deleted from node map or isn't anymore a fixed point

        EXAMPLES::
            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider._safeIndex(12) 

        .. NOTE::
            O(log(m))
        """
        if not self.isFixedPoint(index):
            return
        splayTree = SplayTree()
        splayTree.insert(index)
        self.nodeMap[index] = splayTree.getNode(index)
        self.nodeMap[index].index = index

    # Def swapLabel of index and otherIndex but keep their topology
    # OK
    # O(log(n))
    def swapIndex(self, index, otherIndex):
        """
        This will only swap the label of the node associated to index and otherIndex while keeping
        a relabelling of some sort
        INPUt:
            index,otherIndex

        EXAMPLES::
            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider.getCycleList(7)
            [7, 8, 9, 11]
            sage: provider.swapIndex(9,8)
            sage: provider.getCycleList(7)
            [7, 9, 8, 11]
        .. NOTE::
            O(log(n))
        """
        if index not in self.nodeMap and otherIndex not in self.nodeMap:
            return
        if index not in self.nodeMap and otherIndex in self.nodeMap:
            self.nodeMap[index] = self.nodeMap[otherIndex]
            self.nodeMap[index].index = index
            self.nodeMap.pop(otherIndex)
            return
        if index in self.nodeMap and otherIndex not in self.nodeMap:
            self.nodeMap[otherIndex] = self.nodeMap[index]
            self.nodeMap[otherIndex].index = otherIndex
            self.nodeMap.pop(index)
            return
        tmp = self.nodeMap[index]
        self.nodeMap[index] = self.nodeMap[otherIndex]
        self.nodeMap[otherIndex] = tmp

        self.nodeMap[index].index = index
        self.nodeMap[otherIndex].index = otherIndex

    # OK
    # O(sizeOfCycle)
    def getCycleList(self, index):
        """
        INPUT:
            index
        OUTPUT:
            A list of all index in the cycle of index

        EXAMPLES::

            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider.getCycleList(7)
            [7, 8, 9, 11] 

        .. NOTE::
            O(t+log(n)) where t is the size of the cycle of index
        """
        if self.isFixedPoint(index):
            return [index]
        splayTree = self.getSplayTree(index)
        return splayTree.indexList()

    # Detach the node if it isn't a fixed point and make it a fixed point
    # O(log(n))
    # OK
    def detach(self, index):
        """
        This will make index a fixed point

        INPUT: 
            index

        EXAMPLES::

            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider.getCycleList(7)
            [7, 8, 9, 11]
            sage: provider.detach(9)
            sage: provider.getCycleList(7)
            [7, 8, 11]
            sage: provider.isFixedPoint(9)
            True

        .. NOTE::
            O(log(n))
        """
        if self.isFixedPoint(index):
            return
        # Delete the key of index from it splay tree
        splayTree = self.getSplayTree(index)
        value = splayTree.root.value + splayTree.root.offset
        splayTree.delete(value)

        # Delete index from the nodeMap
        # We can only have no fixed point inside it
        self.nodeMap.pop(index)

        # If there is only one element left in the old cycle of index
        # We must also make him fixed point
        if splayTree.size() == 1:
            node = splayTree.root
            self.detach(node.index)

    # O(1)
    # OK
    def isFixedPoint(self, index):
        """
        INPUT:
            index
        OUTPUT: 
            A boolean indicating if index is a fixed point or not 

        EXAMPLES::
            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider.getCycleList(7)
            [7, 8, 9, 11]
            sage: provider.detach(9)
            sage: provider.getCycleList(7)
            [7, 8, 11]
            sage: provider.isFixedPoint(9)
            True
            sage: provider.isFixedPoint(11)
            False

        .. NOTE::
            O(1)
        """
        return index not in self.nodeMap

    # OK
    # O(log(n))
    def getValue(self, index):
        """
        This will return the key associated to index in the splay tree corresponding to his cycle,
        while making sure that the node asscociated to index is the root of the tree
        INPUT:
            index not a fixed point

        OUTPUT:
            The value associated to index in the splay tree corresponding to his cycle 

        EXAMPLES::

            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider.getValue(11)
            3

        .. NOTE::
            O(log(n))
        """
        if self.isFixedPoint(index):
            raise ValueError(
                f"A fixed point doesn't have a value: {index} is a fixed point")
        splayTree = self.getSplayTree(index)
        return splayTree.root.value + splayTree.root.offset

    # OK
    # O(log(n))
    def getSplayTree(self, index):
        """
        Returns the splay tree associated to index while making sure that the node associated to index become the root
        index must not be a fixed point otherwise an error will be raised

        INPUT:
            index not a fixed point

        OUTPUT:
            The splay associated to index 

        EXAMPLES::

            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider.getSplayTree(1) == provider.getSplayTree(5)
            True 

        .. NOTE::

            O(log(n))

        """
        if self.isFixedPoint(index):
            raise ValueError(
                f"A fixed point doesn't have a splay tree : {index} is a fixed point")
        return self.nodeMap[index].getSplayTree()
    # OK
    # O(log(n))

    def makeMin(self, index):
        """
        This will make index the max and the root of the splay tree, while keeping the same cycle topology
        INPUT:
            index

        EXAMPLES::

            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider.getSplayTree(1).min() == provider.getValue(5)
            False
            sage: provider.makeMin(5)
            sage: provider.getSplayTree(1).min() == provider.getValue(5)
            True 

        .. NOTE::
            O(log(n))

        """
        if self.isFixedPoint(index):
            return
        value = self.getValue(index)
        splayTree = self.getSplayTree(index)

        # Note that .min on a splay tree make the min become the root
        if splayTree.min() == value:
            return

        left, right = splayTree.split(value - 1 / 2)
        left.shift(right.max() + 1 - left.min())
        left.merge(right)
    # OK
    # O(log(n))

    def makeMax(self, index):
        """
        This will make index the max and the root of the splay tree, while keeping the same cycle topology

        INPUT:
            index

        EXAMPLES::

            sage: provider = CycleUtilsProvider([(1,5),(7,8,9,11)])
            sage: provider.getSplayTree(7).max() == provider.getValue(7)
            False
            sage: provider.makeMax(7)
            sage: provider.getSplayTree(7).max() == provider.getValue(7)
            True

        .. NOTE::
            O(log(n))
        """
        if self.isFixedPoint(index):
            return
        value = self.getValue(index)
        splayTree = self.getSplayTree(index)

        # Note that .max on a splay tree make the max become the root
        if splayTree.max() == value:
            return
        left, right = splayTree.split(value + 1 / 2)
        left.shift(right.max() + 1 - left.min())
        left.merge(right)

    # OK
    # O(log(n))
    def merge(self, beforeIndex, afterIndex):
        """
        beforeIndex and afterIndex must be on different cycle , it will merge both cycle in the following manner C = ...->beforeIndex and D = afterIndex-> ... 
        it will create the following ...->beforeIndex->afterIndex-> ...
        INPUT:
            beforeIndex,afterIndex two indexes not on the same cycle, an error will be raised otherwise

        EXAMPLES::

            sage: provider = CycleUtilsProvider([(4242,1,5,42),(7,8,9,11)])
            sage: provider.merge(1,8)
            sage: provider.getCycleList(1)
            [5, 42, 4242, 1, 8, 9, 11, 7]

        .. NOTE::
            O(log(n))
        """
        if self.sameCycle(beforeIndex, afterIndex):
            raise ValueError(
                f"Cannot merge detache two indexes in the same cycle  : {beforeIndex} and {afterIndex}")

        self._safeIndex(beforeIndex)
        self._safeIndex(afterIndex)

        self.makeMax(beforeIndex)
        self.makeMin(afterIndex)

        leftTree = self.getSplayTree(beforeIndex)
        rightSplayTree = self.getSplayTree(afterIndex)

        rightSplayTree.shift(leftTree.max() + 1 - rightSplayTree.min())
        rightSplayTree.merge(leftTree)

    # OK
    # O(log(n))
    def cut(self, startIndex, endIndex):
        """
        startIndex and endIndex must be on the same cycle , it will cut their cycle into two part one startIndex...endIndex ,
        and the rest
        INPUT:
            startIndex and endIndex they must be on the same cycle, an error will be raise otherwise

        EXAMPLES::

            sage: provider = CycleUtilsProvider([(4242,1,5,424242,42),(7,8,9,11)])
            sage: provider.cut(5,424242)
            sage: provider.getCycleList(1)
            [4242, 1, 42]
            sage: provider.getCycleList(5)
            [5, 424242]

        .. NOTE::
            O(log(n))
        """
        if not self.sameCycle(startIndex, endIndex):
            raise ValueError(
                f"Cut can only be called for index in the same cycle {startIndex} and {endIndex} are not")
        if startIndex == endIndex:
            self.detach(startIndex)
            return
        startValue = self.getValue(startIndex)
        endValue = self.getValue(endIndex)
        splayTree = self.getSplayTree(startIndex)

        rest = None
        if startValue < endValue:
            left, otherRight = splayTree.split(endValue + 1 / 2)
            otherLeft, _ = left.split(startValue - 1 / 2)
            rest = otherRight.merge(otherLeft)
        else:
            otherLeft, right = splayTree.split(startValue - 1 / 2)
            left, otherRight = otherLeft.split(endValue + 1 / 2)
            left.shift(right.max() + 1 - left.min())
            left.merge(right)
            rest = otherRight

        if rest.size() == 1:
            self.nodeMap.pop(rest.indexList()[0])
