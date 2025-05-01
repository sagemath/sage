class CyclicChainedList:
    """
    This is an internal class represent cyclic chained list node used in the class RotatingPermutation
    Note that you must be careful when directly manipulating the object to not create more than one cycle , all
    the basic implemented method are guaranteed to not alter this invariant.
    """

    def __init__(self, val):
        """
        val: the value contained in the node
        """
        self.nxt = self
        self.prev = self
        self.val = val

    def getValList(self):
        """
        Return: a list containing the value of each node in the same cycle as self
        O(1)
        """
        cycleNode = []
        cycleNode.append(self)
        curNode = self.nxt
        while curNode != self:
            cycleNode.append(curNode)
            curNode = curNode.nxt

        return [node.val for node in cycleNode]

    def insertAfter(self, otherNode):
        """
        This function will insert otherNode after self i.e if self is in a cycle  A->self->A where  A is the rest of the cycle and that otherNode is
        in a cycle of the form otherNode ->R ->otherNode where R represent the rest of the cycle this function will fuse them into a cycle
        of the form A->self->otherNode->R->A
        O(1)
        """
        oldNxt = self.nxt
        oldPrevOther = otherNode.prev
        self.nxt = otherNode
        otherNode.prev = self
        oldPrevOther.nxt = oldNxt
        oldNxt.prev = oldPrevOther

    def insertBefore(self, otherNode):
        """
        This function will insert otherNode before self i.e if self is in a cycle  A->self->A where  A is the rest of the cycle and that otherNode is
        in a cycle of the form otherNode ->R ->otherNode where R represent the rest of the cycle this function will fuse them into a cycle
        of the form A->otherNode->R->self->A
        O(1)
        """
        self.prev.insertAfter(otherNode)

    def remove(self):
        """
        If self isn't an isolated node(i.e such that self.prev = self.nxt = self) it will remove self from the where it is present
        otherwise it won't do anything
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
