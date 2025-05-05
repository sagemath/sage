# This file contain an implementation of a custom
# Variant of SplayTree used in CycleUtilsProviderMain
# Be careful when modifying this file to make sure that it is coherent with what is done in CycleUtilsProvider

# OK
def swapRight(node, otherNode):
    """
    Exchange the right child of node and otherNode

    INPUT:
        node,otherNode two SplayNode

    EXAMPLES::

        sage: node = SplayNode()
        sage: otherNode = SplayNode()
        sage: otherNode.right = SplayNode(4,otherNode)
        sage: node.right = SplayNode(3,node)
        sage: swapRight(node,otherNode)
        sage: node.right.value
        4
        sage: otherNode.right.value
        3
        sage: node.right.parent == node and otherNode.right.parent == otherNode
        True

    .. NOTE::
        O(1)

    """
    tmp = node.right
    node.right = otherNode.right
    otherNode.right = tmp

    if node.right is not None:
        node.right.parent = node

    if otherNode.right is not None:
        otherNode.right.parent = otherNode


# OK
def swapLeft(node, otherNode):
    """
    Swap the left child of node and otherNode
    INPUT:
        node,otherNode two SplayNode

    EXAMPLES::

        sage: node = SplayNode()
        sage: otherNode = SplayNode()
        sage: node.left = SplayNode(3,node)
        sage: otherNode.left = SplayNode(4,otherNode)
        sage: swapLeft(node,otherNode)
        sage: node.left.value
        4
        sage: otherNode.left.value
        3
        sage: node.left.parent == node and otherNode.left.parent == otherNode
        True

    .. NOTE:: 
        O(1)

    """
    tmp = node.left
    node.left = otherNode.left
    otherNode.left = tmp

    if node.left is not None:
        node.left.parent = node

    if otherNode.left is not None:
        otherNode.left.parent = otherNode


# OK
def SwapNonTopologicalExceptIndex(node, otherNode):
    """
    Swap the following attribut cnt,value,offset,splayTree of node and otherNode

    INPUT:

        node,otherNode two SplayNode

    EXAMPLES::
        sage: node = SplayNode()
        sage: otherNode = SplayNode()
        sage: node.offset = 7
        sage: otherNode.offset = 8
        sage: SwapNonTopologicalExceptIndex(node,otherNode)
        sage: node.offset
        8
        sage: otherNode.offset
        7

    .. NOTE::
        O(1)
    """
    # Value
    tmp = node.value
    node.value = otherNode.value
    otherNode.value = tmp

    # Offset
    tmp = node.offset
    node.offset = otherNode.offset
    otherNode.offset = tmp

    # cnt
    tmp = node.cnt
    node.cnt = otherNode.cnt
    otherNode.cnt = tmp

    # SplayTree
    tmp = node.splayTree
    node.splayTree = otherNode.splayTree
    otherNode.splayTree = tmp


# OK
def makeParentKnow(node):
    """
    Make the parent know(by changing is left or right attribute) that node is his new children

    INPUT:
        node a SplayNode

    EXAMPLES::

        sage: node = SplayNode(4)
        sage: child = SplayNode(3,node)
        sage: node.left,node.right
        (None, None)
        sage: makeParentKnow(child)
        sage: node.left == child
        True 

    .. NOTE::
        O(1)
    """
    if node.parent is None:
        return
    parent = node.parent
    if isLeftChild(parent, node):
        parent.left = node
    else:
        parent.right = node


# OK
def swapNodeButNotIndex(node, otherNode):
    """
    Swap every attribute except index

    INPUT:
        node,otherNode two SplayNode

    EXAMPLES::
        sage: node = SplayNode(4)
        sage: node.index =54
        sage: otherNode = SplayNode(12)
        sage: otherNode.index = 42
        sage: swapNodeButNotIndex(node,otherNode)
        sage: node.value
        12 
        sage: otherNode.value
        4
        sage: otherNode.index
        42
        sage: node.index
        54

    .. NOTE::
        O(1)
    """
    # OK
    if node == otherNode:
        return

    if node.parent == otherNode.parent:
        swapLeft(node, otherNode)
        swapRight(node, otherNode)
        SwapNonTopologicalExceptIndex(node, otherNode)
        makeParentKnow(node)
        makeParentKnow(otherNode)
        return

    # OK
    if node.parent == otherNode:
        swapNodeButNotIndex(otherNode, node)
        return

    # OK
    if otherNode.parent == node:
        swapLeft(node, otherNode)
        swapRight(node, otherNode)

        otherNode.parent = node.parent

        SwapNonTopologicalExceptIndex(node, otherNode)

        makeParentKnow(otherNode)

        node.parent = otherNode
        return

    # Left
    swapLeft(node, otherNode)
    # Right
    swapRight(node, otherNode)
    # Non topological field(i.e != left,right,parent) and != index
    SwapNonTopologicalExceptIndex(node, otherNode)

    # Parent
    tmp = node.parent
    node.parent = otherNode.parent
    otherNode.parent = tmp

    makeParentKnow(node)
    makeParentKnow(otherNode)


def numberOfElement(node):
    """
    Returns node.cnt or 0 if node is None

    EXAMPLES::

        sage: numberOfElement(None)
        0 

    .. NOTE::
        O(1)
    """
    return 0 if node is None else node.cnt


def valueToTheLeft(parentValue, value):
    """
    INPUT:
        parentValue ,value

    OUTPUT:
        Returns a boolean indciating if parentValue>value

    EXAMPLES::
        sage: from SplayTree import valueToTheLeft
        sage: valueToTheLeft(7,8)
        False

    .. NOTE::
        O(1)
    """
    return parentValue > value


def valueToTheRight(parentValue, value):
    """
    INPUT:
        parentValue,value

    OUTPUT::
        Returns a boolean indciating if parentValue<value

    EXAMPLES::

        sage: valueToTheRight(7,8)
        True

    .. NOTE::
        O(1)
    """
    return parentValue < value


def isLeftChild(parentNode, node):
    """

    INPUT:
        parentNode,node

    OUTPUT:
        Returns a boolean indicating if node.value<parentNode.value


    EXAMPLES::

        sage: node = SplayNode(4)
        sage: child = SplayNode(3,node)
        sage: makeParentKnow(child)
        sage: isRightChild(node,child)
        False
        sage: isLeftChild(node,child)
        True

    .. NOTE::
        O(1)
    """
    return valueToTheLeft(parentNode.value, node.value + node.offset)


def isRightChild(parentNode, node):
    """
    INPUT:
        parentNode,node 

    OUTPUT:

        Returns a boolean indicating if node.value>parentNode.value

    EXAMPLES::

        sage: node = SplayNode(4)
        sage: child = SplayNode(3,node)
        sage: makeParentKnow(child)
        sage: isRightChild(node,child)
        False
        sage: isLeftChild(node,child)
        True


    .. NOTE::
        O(1)
    """
    return valueToTheRight(parentNode.value, node.value + node.offset)


class SplayNode:
    """
    This class is an internal class used in SplayTree and CycleUtilsProvider,
    most operation are O(log(n)) amortized where n is the size of the tree iff there is
    a splay after each of them , which is done most of the time outside of this class by the splay tree
    on which it is attached.
    """

    def __init__(self, value=None, parent=None) -> None:
        """
        Init the splaynode 

        INPUT:
            value,parent both default to None

        EXAMPLES::
            sage: node = SplayNode(4)

        .. NOTE::
            O(1)
        """
        # The value contained in the node note that
        # It doesn't represent the key associated to node
        # Which is the sum of the value and offset
        self.value = value

        # Left child
        self.left = None
        # Right child
        self.right = None

        # Parent if it is None it means this is the
        # Root
        self.parent = parent

        # The number of node in the subtree of node
        self.cnt = 0

        # The offset contained in the node
        self.offset = 0

        # The splay tree on which the node is attached
        self.splayTree = None

        # The index which is associated
        self.index = None

    def SafeSplay(self):
        """
        This will splay self while being sure that after self.splayTree point to the correct value
        and oldRoot.splayTree is None.

        EXAMPLES::

            sage: node = SplayNode(4)
            sage: node.SafeSplay()

        .. NOTE::
            O(log(n))
        """
        oldRoot = self.getRoot()
        if oldRoot != self:
            self.splayTree = oldRoot.splayTree
            oldRoot.splayTree = None
        self.splay()

    def getSplayTree(self):
        """
        OUTPUT::

            This return a reference to the splay tree of self , and make self the root of the splay tree or None
            if it doesn't exist.

        EXAMPLES::

            sage: node = SplayNode(4)
            sage: node.getSplayTree() is None
            True 

        .. NOTE::
            O(log(n))
        """
        # It is guaranteed after colling this function that self is the root
        root = self.getRoot()
        self.splay()
        splayTree = root.splayTree
        if splayTree is not None:
            splayTree.changeRoot(self)
        return splayTree

    def getRoot(self):
        """
        OUTPUT::
            The root of the tree on which self is attached

        EXAMPLES::

            sage: child = SplayNode(3,node)
            sage: node = SplayNode(4)
            sage: child = SplayNode(3,node)
            sage: makeParentKnow(child)
            sage: child.getRoot() == node
            True

        .. NOTE::
            O(log(n)), don't forget to splay on self
            to guarantee the amortized complexity
        """
        node = self
        while not node.isRoot():
            node = node.parent

        return node

    def height(self):
        """

        OUTPUT:
            The height of the subtree of self

        EXAMPLES::

            sage: node = SplayNode(4)
            sage: child = SplayNode(3,node)
            sage: makeParentKnow(child)
            sage: child.height()
            0
            sage: node.height()
            1

        .. NOTE:
            O(n),written recursively hence if you have
            too small recursion limit it can crash

        """
        if self.isEmpty():
            return -1
        maxChildHeight = -1
        if self.right is not None:
            maxChildHeight = self.right.height()
        if self.left is not None:
            maxChildHeight = max(self.left.height(), maxChildHeight)
        return 1 + maxChildHeight

    def sortedList(self, offset=0):
        """

        INPUT:
            offset (default to 0) the offset to add
            to each value of the subtree

        OUTPUT:

            Returns the sorted list of key in the subtree 

        EXAMPLES::

            sage: node = SplayNode(4)
            sage: child = SplayNode(3,node)
            sage: makeParentKnow(child)
            sage: node.sortedList()
            [3, 4]

        .. NOTE::
            O(n)
        """
        if self.isEmpty():
            return []
        out = [self.value + offset + self.offset]
        if self.right is not None:
            out = out + self.right.sortedList(self.offset + offset)
        if self.left is not None:
            out = self.left.sortedList(offset + self.offset) + out
        return out

    def indexList(self):
        """

        OUTPUT:

        Returns the list of index in the subtree associated to self in the order of their key


        EXAMPLES::

            sage: node = SplayNode(4)
            sage: child = SplayNode(3,node)
            sage: makeParentKnow(child)
            sage: node.index = 4
            sage: child.index = 75
            sage: node.indexList()
            [75, 4] 

        ..NOTE::
            O(n)

        """
        if self.isEmpty():
            return []
        out = [self.index]
        if self.right is not None:
            out = out + self.right.indexList()
        if self.left is not None:
            out = self.left.indexList() + out
        return out

    def isEmpty(self):
        """

        OUTPUT:
        A boolean indcating if self is empty

        EXAMPLES::

            sage: node = SplayNode()
            sage: node.isEmpty()
            True


        .. NOTE::
            O(1)
        """
        return self.value is None

    def findSmallestGreater(self, value):
        """
        INPUT:

            value

        OUTPUT:
        A couple node,offset such that if there is smallest greater key in self that value node.offset+offset is this
        value.

        EXAMPLES::

            sage: node = SplayNode(22)
            sage: child = SplayNode(3,node)
            sage: makeParentKnow(child)
            sage: node.findSmallestGreater(1)[0] == child
            True

        .. NOTE::
            O(log(n)),to keep the O(log(m)) amortized value you must splay node.
        """

        if self.isEmpty():
            return self, 0

        node = self
        offset = self.offset
        bestNode = node
        bestOffset = self.offset
        while True:
            if value == node.value + offset:
                return node, offset
            if valueToTheLeft(node.value + offset, value):
                bestNode = node
                bestOffset = offset
                if node.left is None:
                    return bestNode, bestOffset
                node = node.left
            else:
                if node.right is None:
                    return bestNode, bestOffset
                node = node.right
            offset += node.offset

    def findBiggestSmaller(self, value):
        """
        INPUT:
            value

        OUTPUT:
        A couple offset,node such that if it exist node.value+offset is the biggest value in
        self subtree such that <= value, otherwise the smallest bigger

        EXAMPLES::

            sage: node = SplayNode(22)
            sage: child = SplayNode(3,node)
            sage: makeParentKnow(child)
            sage: node.findBiggestSmaller(33)[0] == node
            True 

        .. NOTE::
        O(log(n)),you must splay on node to guarantee the amortized O(log(n)) complexity
        """
        if self.isEmpty():
            return self, 0

        node = self
        offset = self.offset
        bestNode = node
        bestOffset = self.offset
        while True:
            if value == node.value + offset:
                return node, offset
            if valueToTheLeft(node.value + offset, value):
                if node.left is None:
                    return bestNode, bestOffset
                node = node.left
            else:
                bestNode = node
                bestOffset = offset
                if node.right is None:
                    return bestNode, bestOffset
                node = node.right
            offset += node.offset

    def insert(self, newValue):
        """


        Insert newValue inside the tree 

        INPUT:
            newValue

        OUTPUT:
            return b,node,offset such that b=True if it was a real new value otherwise False
            and node.value+offset = value

        EXAMPLES::

            sage: node = SplayNode(22) 
            sage: isNew,newNode,offset = node.insert(34)
            sage: isNew
            True
            sage: newNode.value + offset == 34
            True

        .. NOTE::
        O(log(n))
        """
        if self.isEmpty():
            self.value = newValue
            self.addCountUpward(1)
            return True, self, 0
        node = self
        offset = self.offset
        while (True):
            if node.value + offset == newValue:
                return False, node, offset
            if valueToTheLeft(node.value + offset, newValue):
                if node.left is None:
                    node.left = SplayNode(newValue - offset, node)
                    node.left.addCountUpward(1)
                    return True, node.left, offset
                else:
                    node = node.left
            else:
                if node.right is None:
                    node.right = SplayNode(newValue - offset, node)
                    node.right.addCountUpward(1)
                    return True, node.right, offset
                else:
                    node = node.right
            offset += node.offset

    def addCountUpward(self, toAdd):
        """

        Add toAdd to the cnt attribute of all the node in the path toward the root from self

        INPUT:

            toAdd

        EXAMPLES::

            sage: node = SplayNode(22)
            sage: child = SplayNode(3,node)
            sage: cchild = SplayNode(1,child)
            sage: makeParentKnow(child)
            sage: makeParentKnow(cchild)
            sage: child.cnt = 1
            sage: node.cnt = 2
            sage: cchild.cnt
            0
            sage: child.addCountUpward(4)
            sage: child.cnt,node.cnt
            (5, 6)
            sage: cchild.cnt
            0

        .. NOTE::
        O(log(n)),not intended to be used by the user, you should splay on node after calling it.
        """
        node = self
        node.cnt += toAdd
        while not node.isRoot():
            node = node.parent
            node.cnt += toAdd

    def find(self, value):
        """
        Find a node such that node.value+node.offset() == value, return the offset and the node,
        if it doesn't exist it will return the node,offset the smallest real value i.e node.value+offset < value

        INPUT:
            value

        OUTPUT:
            node,offset corresponding to the above description

        EXAMPLES::

            sage: node = SplayNode(22)
            sage: node.find(22)[0] == node
            True 

        .. NOTE:
        O(log(n)),should splay on node to guarantee amortized log complexity
        """
        if self.isEmpty():
            return self, 0
        node = self
        offset = self.offset
        while True:
            if value == node.value + offset:
                return node, offset
            if valueToTheLeft(node.value + offset, value):
                if node.left is None:
                    return node, offset
                node = node.left
            else:
                if node.right is None:
                    return node, offset
                node = node.right
            offset += node.offset

    def min(self):
        """
        OUTPUT:
            node,offset such that node.value+offset is the minimal key

        EXAMPLES::

            sage: node = SplayNode(22)
            sage: child = SplayNode(3,node)
            sage: makeParentKnow(child)
            sage: node.min()[0] == child
            True

        .. NOTE::
        O(log(n)), splay on node to guarantee amortized log time
        """
        minimum = self
        offset = self.offset
        while minimum.left is not None:
            minimum = minimum.left
            offset += minimum.offset
        return minimum, offset

    def delete(self, value):
        """
        Delete value from self subtree and return the nearest node to the deleted one

        OUTPUT:
            node corresponding to the above description

        EXAMPLES::

            sage: node = SplayNode(22)
            sage: child = SplayNode(3,node)
            sage: makeParentKnow(child)
            sage: node.sortedList()
            [3, 22]
            sage: node.delete(3) == node
            True
            sage: node.sortedList()
            [22]

        ..NOTE::
        O(log(n)) , splay on the node returned
        """
        if self.isEmpty():
            return self

        node, offset = self.find(value)
        if node.value + offset != value:
            return node

        if node.left is None and node.right is None:
            if node.isRoot():
                node.value = None
                node.addCountUpward(-1)
                return node
            node.parent.addCountUpward(-1)
            if isLeftChild(node.parent, node):
                node.parent.left = None
                return node.parent
            node.parent.right = None
            return node.parent

        if node.left is None and node.right is not None:
            node.right.parent = node.parent
            if node.parent is None:
                return node.right

            if (isLeftChild(node.parent, node)):
                node.parent.left = node.right
            else:
                node.parent.right = node.right

             # Offset
            node.right.offset += node.offset

            node.parent.addCountUpward(-1)
            return node.right

        if node.right is None and node.left is not None:
            node.left.parent = node.parent
            if node.parent is None:
                return node.left
            if isLeftChild(node.parent, node):
                node.parent.left = node.left
            else:
                node.parent.right = node.left

            # Offset
            node.left.offset += node.offset

            node.parent.addCountUpward(-1)
            return node.left

        rep, repOffset = node.right.min()

        if isLeftChild(rep.parent, rep):
            rep.parent.left = rep.right
        else:
            rep.parent.right = rep.right

        if rep.right is not None:
            rep.right.parent = rep.parent
            rep.right.offset += rep.offset

        node.value = rep.value + repOffset
        rep.parent.addCountUpward(-1)

        # This line is mainly here to not break the association
        # in CycleUtilsProvider between index <-> node
        swapNodeButNotIndex(node, rep)

        return node.parent

    def max(self):
        """
        OUTPUT:
        Returns node,offset such that node.value+offset is the 
        maximum key in the subtree of self

        EXAMPLES::

            sage: node = SplayNode(22)
            sage: child = SplayNode(3,node)
            sage: makeParentKnow(child)
            sage: node.max()[0] == node
            True

        ..NOTE:
        O(log(n)),splay on the returned node to guarantee log amortized time
        """
        maximum = self
        offset = self.offset
        while maximum.right is not None:
            maximum = maximum.right
            offset += maximum.offset
        return maximum, offset

    def _isBst(self, offset=0):
        """ 
        Verify that self is a binary search tree and raise an error otherwise

        OUTPUT:
        returns the min and max of the subtree of self

        EXAMPLES::

            sage: node = SplayNode(22)
            sage: node._isBst()
            (22, 22)

        .. NOTE::
        O(n),written recursively hence if a too small recursion limit may crash
        """
        if self.isEmpty():
            return None, None
        min = self.value + self.offset + offset
        max = self.value + self.offset + offset
        if self.left is not None:
            minLeft, maxLeft = self.left._isBst(self.offset + offset)
            assert valueToTheLeft(self.value + offset + self.offset, maxLeft)
            min = minLeft
        if self.right is not None:
            minRight, maxRight = self.right._isBst()
            assert valueToTheRight(self.value + offset + self.offset, minRight)
            max = maxRight
        return min, max

    def isBst(self):
        """
        Raise an error if self isn't a bst

        EXAMPLES::

            sage: try:
            ....:     node.isBst()
            ....: except:
            ....:     print("No")
            ....:

        ..NOTE:
        O(n), written recursively if too small recursion limit it may crash
        """
        self._isBst()

    def rightRotation(self):
        """
        Do a right rotation while updating the attribute accordingly

        OUTPUT::
            The node which is after the call the father of self 
            if self has a left child otherwise self

        EXAMPLES::

            sage: node = SplayNode(22)
            sage: child = SplayNode(3,node)
            sage: cchild = SplayNode(1,child)
            sage: makeParentKnow(child)
            sage: makeParentKnow(cchild)
            sage: node.rightRotation() == child
            True
            sage: child.isRoot()
            True 

        .. NOTE::
        O(1)
        """
        if self.isEmpty():
            return self

        oldLeft = self.left

        if oldLeft is None:
            return self

        # Parent
        if self.parent is not None:
            parent = self.parent
            if isLeftChild(parent, self):
                parent.left = oldLeft
            else:
                parent.right = oldLeft

        if oldLeft.right is not None:
            oldLeft.right.parent = self

        oldLeft.parent = self.parent
        self.parent = oldLeft

        # Offset
        x = oldLeft.offset
        y = self.offset

        oldLeft.offset = x + y
        self.offset = -x
        if oldLeft.right is not None:
            oldLeft.right.offset += x

        # Cnt
        totCnt = self.cnt

        self.cnt = totCnt - oldLeft.cnt + numberOfElement(oldLeft.right)

        oldLeft.cnt = totCnt

        # Child

        self.left = oldLeft.right
        oldLeft.right = self

        return oldLeft

    def leftRotation(self):
        """
        Do a left rotation while updating attribute accordingly

        OUTPUT:
            The node which is after the call the father of self
            if self has a left child otherwise self

        EXAMPLES::

            False
            sage: node = SplayNode(22)
            sage: child = SplayNode(33,node)
            sage: cchild = SplayNode(44,child)
            sage: makeParentKnow(cchild)
            sage: makeParentKnow(child)
            sage: node.leftRotation() == child
            True
            sage: child.isRoot()
            True

        .. NOTE::
        O(1)
        """
        if self.isEmpty():
            return self

        oldRight = self.right
        if oldRight is None:
            return self

        # Parent
        if self.parent is not None:
            parent = self.parent
            if isLeftChild(parent, self):
                parent.left = oldRight
            else:
                parent.right = oldRight

        if oldRight.left is not None:
            oldRight.left.parent = self

        oldRight.parent = self.parent
        self.parent = oldRight

        # Offset
        y = self.offset
        x = oldRight.offset
        oldRight.offset = x + y
        self.offset = -x
        if oldRight.left is not None:
            oldRight.left.offset += x

        # Cnt
        totCnt = self.cnt

        self.cnt = totCnt - oldRight.cnt + numberOfElement(oldRight.left)

        oldRight.cnt = totCnt

        # Child
        self.right = oldRight.left
        oldRight.left = self

        return oldRight

    def isRoot(self):
        """

        OUTPUT:
        Returns a boolean indcating if self is the root of the tree or not

        EXAMPLES::
            sage: node = SplayNode(22)
            sage: child = SplayNode(33,node)
            sage: cchild = SplayNode(44,child)
            sage: makeParentKnow(cchild)
            sage: makeParentKnow(child)
            sage: node.leftRotation() == child
            True
            sage: child.isRoot()
            True

        .. NOTE::
        O(1)
        """
        return self.parent is None

    def getOffset(self):
        """
        OUTPUT:
        Returns the sum of offset from the root to self 

        EXAMPLES::

            sage: node = SplayNode(22)
            sage: child = SplayNode(33,node)
            sage: cchild = SplayNode(44,child)
            sage: makeParentKnow(child)
            sage: makeParentKnow(cchild)
            sage: node.offset = 5
            sage: child.offset = 2
            sage: cchild.offset = 1
            sage: cchild.getOffset()
            8

        .. NOTE::
        O(log(m))
        """
        offset = self.offset
        node = self
        while not node.isRoot():
            node = node.parent
            offset += node.offset
        return offset

    def splay(self):
        """
        This applied splay the node toward the root of the tree,
        note that it isn't safe mainly cause it doesn't update the splayTree attribute of the root
        You should use SafeSplay instead

        OUTPUT:
            self


        EXAMPLES::

            sage: node = SplayNode(22)
            sage: child.offset = 2
            sage: cchild = SplayNode(44,child)
            sage: makeParentKnow(cchild)
            sage: makeParentKnow(child)
            sage: cchild.isRoot()
            True

        .. NOTE::
        O(log(m))
        """
        while not self.isRoot():
            parent = self.parent

            if parent.isRoot():
                if isLeftChild(parent, self):
                    parent.rightRotation()
                else:
                    parent.leftRotation()
                continue

            grandparent = parent.parent

            if isLeftChild(grandparent, parent) and isLeftChild(parent, self):
                grandparent.rightRotation()
                parent.rightRotation()

            elif isRightChild(grandparent, parent) and isRightChild(parent, self):
                grandparent.leftRotation()
                parent.leftRotation()

            elif isRightChild(grandparent, parent) and isLeftChild(parent, self):
                parent.rightRotation()
                grandparent.leftRotation()

            else:
                parent.leftRotation()
                grandparent.rightRotation()
        return self


class SplayTree():
    def __init__(self, lst=[], root=None):
        """
        This class implement a special version of splay tree where it is possible to schift all the
        value contained in the splay tree by a constant.

        Initialises the splay tree

        INPUT:
            -lst: A list of element to add too self
            -root: A possible Node on which to root the splayTree if set to None it will be set to
                Node()

        EXAMPLES::

            sage: sp = SplayTree([0,1,2])

        .. NOTE::
        O(n * log(max(n,m)) ) where n is the length of lst and m is the number of element in the subtree
            of root
        """

        self.root = root
        if self.root is None:
            self.root = SplayNode()
        self.root.splayTree = self
        self.valid = True
        self.insertList(lst)

    def changeRoot(self, root):
        """
        Change the root of the splay tree to root

        INPUT: 
            root


        EXAMPLES::

            sage: sp = SplayTree([0,1,2])
            sage: oldRoot = sp.root
            sage: sp.changeRoot(SplayNode())
            sage: sp.root != oldRoot
            True

        O(log(m)) where m is the size of the tree of root,not intended to be used directly
        """
        if not root.isRoot():
            raise ValueError("The argument must be the root of his tree")
        self.root.splayTree = None
        self.root = root
        self.root.splayTree = self

    def shift(self, toAdd):
        """
        This will shift all the value contained in the tree by toAdd
        INPUT:
            toAdd

        EXAMPLES::

            sage: sp = SplayTree([0,1,2])
            sage: sp.toList()
            [0, 1, 2]
            sage: sp.shift(42)
            sage: sp.toList()
            [42, 43, 44]

        .. NOTE::
            O(1)
        """
        self.checkValid()
        self.root.offset += toAdd

    def checkValid(self):
        """
        Check whether self is a valid SplayTree if not raise an Error

        EXAMPLES::

            sage: sp = SplayTree([0,1,2])
            sage: sp.checkValid()

        .. NOTE::
            O(1)
        """
        if self.valid:
            return
        raise ValueError("This isn't a valid instance of a splayTree anymore")

    def split(self, value):
        """
        This will split self into two splayTree a,b such that a contains all the element <= value
        and b all the element >value.You should note that a == self after the operation.
        INPUT:
            value : The value on which to split
        OUTPUT:
            (a,b) such that a contains all the element <= value and b all the element > value

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.toList()
            [0, 1, 2, 7, 22, 33]
            sage: a,b = sp.split(7)
            sage: a.toList()
            [0, 1, 2, 7]
            sage: b.toList()
            [22, 33]

        ..NOTE::
        O(log(n)) where n is the size of self
        """
        self.checkValid()
        if self.isEmpty():
            return self, SplayTree()

        if self.findBiggestSmaller(value) is None:
            return SplayTree(), self
        right = self.root.right
        self.root.right = None
        self.root.cnt -= numberOfElement(right)

        if right is not None:
            right.parent = None
            right.offset += self.root.offset
        return self, SplayTree(root=right)

    def indexList(self):
        """

        OUTPUT::
            index list in the order of the key 

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.indexList()
            [None, None, None, None, None, None]

        """
        self.checkValid()
        return self.root.indexList()

    def toList(self):
        """
        This will return a sorted list of all the key contains in self

        OUTPUT:  
            A sorted list of all the element of self

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.toList()
            [0, 1, 2, 7, 22, 33]

        .. NOTE::
        O(n)
        """
        self.checkValid()
        return self.root.sortedList()

    def merge(self, otherSplayTree):
        """
        This function return a splay tree containing all the element of self and otherSplayTree ,
        you need to have otherSplayTree.max()<self.min() or it will raise an Error, moreover after
        applied self isn't anymore a valid SplayTree and otherSplayTree point to the merged tree.

        INPUT:
            otherSplayTree : The otherSplayTree such that for otherSplayTree.max()<self.min()
        OUTPUT:
            The merge splay tree containing all the element of self and otherSplayTree

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: qp = SplayTree([35,44,42])
            sage: qp.toList()
            [35, 42, 44]
            sage: sp.toList()
            [0, 1, 2, 7, 22, 33]
            sage: qp.merge(sp) == sp
            True
            sage: sp.toList()
            [0, 1, 2, 7, 22, 33, 35, 42, 44]

        .. NOTE::
        O(log(n)+log(m)) where n is the size of self and m is the size of otherSplayTree
        """
        self.checkValid()

        if otherSplayTree.isEmpty():
            otherSplayTree.root = self.root
            self.root.splayTree = otherSplayTree
            self.valid = False
            return otherSplayTree

        if self.isEmpty():
            self.valid = False
            return otherSplayTree

        if self.min() <= otherSplayTree.max():
            raise ValueError(
                "All the element of the left splay tree must be < than the element of self")
        self.valid = False
        assert otherSplayTree.root.right is None

        otherSplayTree.root.right = self.root
        self.root.parent = otherSplayTree.root
        otherSplayTree.root.cnt += self.root.cnt
        self.root.offset -= otherSplayTree.root.offset

        return otherSplayTree

    def _height(self):
        """
        This function returns the height of self

        OUTPUT: 
            the height of self

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp._height()
            3 

        .. NOTE::
        O(n) where n = self.size(),written recursively hence with too small rrecusion
        limit may crash.
        """
        self.checkValid()
        return self.root.height()

    def isEmpty(self):
        """
        OUTPUT: 
            A boolean indicating if self is Empty

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.isEmpty()
            False 

        .. NOTE::
        O(1)
        """
        self.checkValid()
        return self.root.isEmpty()

    def _isBst(self):
        """

        OUTPUT: 
        A boolean indicating if self is a valid binary search tree

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp._isBst()
            True 


        .. NOTE::
        O(n) where n =  self.size(), written recursively hence with too much small
        recursion limit it may crash
        """

        self.checkValid()
        try:
            self.root.isBst()
            return True
        except BaseException:
            return False

    def _isRoot(self):
        """
        OUTPUT: 
            A boolean indicating if self.root is really the root

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp._isRoot()
            True

        .. NOTE::
        O(1)
        """
        self.checkValid()
        return self.root.isRoot()

    def getNode(self, value):
        """
        The node corresponding to value in the tree if value is in the tree otherwise it return None,
        note that node.value+node.getOffSet() == value and not node.value == value, moreover it will make if it
        exists the node the root

        INPUT: 
            value

        Returns:  
            The node corresponding to value in the tree if value is in the tree otherwise
            it return None

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.getNode(3) !=  sp.getNode(2)
            True

        ..NOTE::
        O(log(n))
        """
        self.checkValid()
        if self.isEmpty():
            return None
        node, _ = self.root.find(value)
        self.changeRoot(node.splay())
        if self.root.value + self.root.offset != value:
            return None
        return self.root

    def insert(self, newValue):
        """
        Insert newValue inside a tree, return a boolean True if newValue was inserted
        and False if newValue was already present in self

        INPUT:
            value

        OUTPUT: 
            True if newValue was inserted and False if newValue was already present in the tree

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.toList()
            [0, 1, 2, 7, 22, 33]
            sage: sp.insert(55)
            True
            sage: sp.toList()
            [0, 1, 2, 7, 22, 33, 55]

        .. NOTE::
        O(log(n)) where n = self.size()
        """

        self.checkValid()
        isNew, node, _ = self.root.insert(newValue)
        self.changeRoot(node.splay())
        return isNew

    def delete(self, value):
        """
        Delete value inside self

        INPUT: value

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.toList()
            [0, 1, 2, 7, 22, 33]
            sage: sp.delete(7)
            sage: sp.toList()
            [0, 1, 2, 22, 33 


        .. NOTE::
        O(log(n)) where n = self.size()
        """
        self.checkValid()
        self.changeRoot(self.root.delete(value).splay())

    def min(self):
        """
        This method return the minimum of self, and make the corresponding node the root
        OUTPUT: 
            the minimum of self if self isn't empty and None otherwise

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.min()
            0

        .. NOTE::
        O(log(n)) where n = self.size()
        """
        self.checkValid()
        if self.isEmpty():
            return None
        node, _ = self.root.min()
        self.changeRoot(node.splay())

        return self.root.value + self.root.offset

    def max(self):
        """
        This method return the maximum of self, and make the corresponding node 
        the root


        OUTPUT: 
            the maximum of self if self isn't empty and None otherwise

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.max()
            33 

        .. NOTE::
        O(log(n)) where n = self.size()
        """

        self.checkValid()
        if self.isEmpty():
            return None
        node, _ = self.root.max()
        self.changeRoot(node.splay())
        return self.root.value + self.root.offset

    def find(self, value):
        """
        This return a boolean indicating if value is in self, and if it exists
        make the corresponding node the root

        INPUT: 
            value
        OUTPUT: 
            A boolean indicating if value is in self

        EXAMPLES:: 

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.find(7)
            True
            sage: sp.find(8)
            False


        .. NOTE::
        O(log(n)) where n = self.size()
        """

        self.checkValid()
        if self.isEmpty():
            return False
        node, _ = self.root.find(value)
        self.changeRoot(node.splay())
        return self.root.value + self.root.offset == value

    def insertList(self, list):
        """
        This will insert the element of list inside self

        INPUT: 
            list, a list of element to add

        EXAMPLES:: 

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.toList()
            [0, 1, 2, 7, 22, 33]
            sage: sp.insertList([42,4242,424242])
            sage: sp.toList()
            [0, 1, 2, 7, 22, 33, 42, 4242, 424242]

        .. NOTE::
        O(m*log(n+m)) where n = self.size() and m = len(list)
        """

        self.checkValid()
        for e in list:
            self.insert(e)

    def findSmallestGreater(self, value):
        """
        This will return the smallest element of self >= value or None if it doesn't exist,
        if the value exist it will be made the root otherwise the greatest smaller will be made the root
        INPUT: 
            value
        OUTPUT: 
            the smallest element of self >= value or None if it  doesn't exist

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.toList()
            [0, 1, 2, 7, 22, 33]
            sage: sp.findSmallestGreater(12)
            22

        .. NOTE::
        O(log(n)) where n = self.size()
        """

        self.checkValid()
        if self.isEmpty():
            return None
        node, _ = self.root.findSmallestGreater(value)
        self.changeRoot(node.splay())
        if self.root.value + self.root.offset < value:
            return None
        return self.root.value + self.root.offset

    def size(self):
        """
        OUTPUT: 
            Size of self

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.toList()
            [0, 1, 2, 7, 22, 33]
            sage: sp.size()
            6

        .. NOTE::
        O(1)
        """
        self.checkValid()
        return self.root.cnt

    def findBiggestSmaller(self, value):
        """
        This will return the smallest element of self <= value or None if it doesn't exist, 
        if the value exist it will be made the root otherwise the smallest bigger will be made the root


        INPUT: 
            value
        OUTPUT: 
            the biggest element of self <= value or None if doesn't exist

        EXAMPLES::

            sage: sp = SplayTree([0,7,1,33,2,22])
            sage: sp.toList()
            [0, 1, 2, 7, 22, 33]
            sage: sp.findBiggestSmaller(8)
            7

        .. NOTE::
        O(log(n)) where n = self.size()
        """

        self.checkValid()
        if self.isEmpty():
            return None
        node, _ = self.root.findBiggestSmaller(value)
        self.changeRoot(node.splay())
        if self.root.value + self.root.offset > value:
            return None
        return self.root.value + self.root.offset
