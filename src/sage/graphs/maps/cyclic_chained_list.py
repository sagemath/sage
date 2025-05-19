"""Define the internal CyclicChainedList class."""


class CyclicChainedList:
    """
    This is an internal class representing cyclic chained list node used in the class RotatingPermutation.

    Note that you must be careful when directly manipulating the object to not create more than one cycle. All
    the basic implemented methods are guaranteed to not alter this invariant.
    """

    def __init__(self, val: int):
        """
        INPUT:

        - ``val`` -- int; the value contained in the node

        EXAMPLES::

            sage: from sage.graphs.maps.cyclic_chained_list import CyclicChainedList
            sage: CyclicChainedList(3)
            NodeOfCyclicChainedList(val=3)
        """
        self.nxt = self
        self.prev = self
        self.val = val

    def __repr__(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.graphs.maps.cyclic_chained_list import CyclicChainedList
            sage: CyclicChainedList(3)
            NodeOfCyclicChainedList(val=3)
        """
        return f"NodeOfCyclicChainedList(val={self.val})"

    def getValList(self) -> list:
        """
        OUTPUT:

        a list containing the value of each node in the same cycle as self

        EXAMPLES::

            sage: from sage.graphs.maps.cyclic_chained_list import CyclicChainedList
            sage: node = CyclicChainedList(3)
            sage: node.insertAfter(CyclicChainedList(4))
            sage: node.getValList()
            [3, 4]

        .. NOTE::

            O(1)
        """
        cycleNode = []
        cycleNode.append(self)
        curNode = self.nxt
        while curNode != self:
            cycleNode.append(curNode)
            curNode = curNode.nxt

        return [node.val for node in cycleNode]

    def insertAfter(self, otherNode: "CyclicChainedList") -> None:
        """
        This function will insert otherNode after self i.e if self is in a cycle  A->self->A where  A is the rest of the cycle and that otherNode is in a cycle of the form otherNode ->R ->otherNode where R represent the rest of the cycle this function will fuse them into a cycle
        of the form A->self->otherNode->R->A

        INPUT:

        otherNode, another node of a CyclicChainedList
        he must be alone in his cycle

        EXAMPLES::

            sage: from sage.graphs.maps.cyclic_chained_list import CyclicChainedList
            sage: node = CyclicChainedList(3)
            sage: node.insertAfter(CyclicChainedList(4))
            sage: node.getValList()
            [3, 4]

        .. NOTE::

            O(1)
        """
        oldNxt = self.nxt
        oldPrevOther = otherNode.prev
        self.nxt = otherNode
        otherNode.prev = self
        oldPrevOther.nxt = oldNxt
        oldNxt.prev = oldPrevOther

    def insertBefore(self, otherNode: "CyclicChainedList") -> None:
        """
        This function will insert otherNode before self i.e if self is in a cycle  A->self->A where  A is the rest of the cycle and that otherNode is
        in a cycle of the form otherNode ->R ->otherNode where R represent the rest of the cycle this function will fuse them into a cycle
        of the form A->otherNode->R->self->A

        INPUT:

        otherNode,another node of a CyclicChainedList
        he must be alone in his cycle

        EXAMPLES::

            sage: from sage.graphs.maps.cyclic_chained_list import CyclicChainedList
            sage: node = CyclicChainedList(3)
            sage: node.insertAfter(CyclicChainedList(4))
            sage: node.insertBefore(CyclicChainedList(5))
            sage: node.getValList()
            [3, 4, 5]

        .. NOTE::

            O(1)
        """
        self.prev.insertAfter(otherNode)

    def remove(self) -> None:
        """
        If ``self`` is not an isolated node (i.e such that self.prev = self.nxt = self) it will remove self from the where it is present
        otherwise it will not do anything

        EXAMPLES::

            sage: from sage.graphs.maps.cyclic_chained_list import CyclicChainedList
            sage: node = CyclicChainedList(3)
            sage: node.insertAfter(CyclicChainedList(4))
            sage: node.insertBefore(CyclicChainedList(5))
            sage: node.getValList()
            [3, 4, 5]
            sage: otherNode = node.nxt
            sage: otherNode.remove()
            sage: node.getValList()
            [3, 5]

        .. NOTE::

            O(1)
        """
        if self.prev == self:
            return
        prev = self.prev
        nxt = self.nxt

        self.nxt = self
        self.prev = self
        prev.nxt = nxt
        nxt.prev = prev
