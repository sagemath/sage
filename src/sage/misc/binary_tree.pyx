# sage_setup: distribution = sagemath-categories
"""
Binary trees

This implements a binary tree in Cython.

AUTHORS:

- Tom Boothby (2007-02-15).  Initial version free for any use (public domain).
"""
from cysignals.memory cimport sig_malloc, sig_free

from cpython.ref cimport PyObject, Py_INCREF, Py_XDECREF

cdef binary_tree_node *BinaryTreeNode(int key, object value) noexcept:
    cdef binary_tree_node *t
    t = <binary_tree_node *>sig_malloc(sizeof(binary_tree_node))
    t.key = key
    t.left = NULL
    t.right = NULL
    Py_INCREF(value)
    t.value = <void *>value
    return t

cdef void free_binary_tree_node(binary_tree_node *self) noexcept:
    Py_XDECREF(<PyObject *>self.value)
    sig_free(self)

cdef inline void binary_tree_dealloc(binary_tree_node *self) noexcept:
    if self != NULL:
        binary_tree_dealloc(self.left)
        binary_tree_dealloc(self.right)
        free_binary_tree_node(self)


cdef void binary_tree_insert(binary_tree_node *self, int key, object value) noexcept:
    if self.key == key:
        return
    elif self.key > key:
        if self.left == NULL:
            self.left = BinaryTreeNode(key, value)
        else:
            binary_tree_insert(self.left, key, value)
    else:
        if self.right == NULL:
            self.right = BinaryTreeNode(key, value)
        else:
            binary_tree_insert(self.right, key, value)

cdef object binary_tree_get(binary_tree_node *self, int key):
    if self.key == key:
        return <object>self.value
    elif self.key > key:
        if self.left == NULL:
            return None
        else:
            return binary_tree_get(self.left, key)
    else:
        if self.right == NULL:
            return None
        else:
            return binary_tree_get(self.right, key)

cdef object binary_tree_delete(binary_tree_node *self, int key):
    cdef object t
    if self.key > key:
        if self.left == NULL:
            return None
        elif self.left.key == key:
            t = <object>self.left.value
            self.left = binary_tree_left_excise(self.left)
            return t
        else:
            return binary_tree_delete(self.left, key)
    else:
        if self.right == NULL:
            return None
        elif self.right.key == key:
            t = <object>self.right.value
            self.right = binary_tree_right_excise(self.right)
            return t
        else:
            return binary_tree_delete(self.right, key)

cdef binary_tree_node *binary_tree_left_excise(binary_tree_node *self) noexcept:
    cdef binary_tree_node *left
    cdef binary_tree_node *cur
    if self.left == NULL:
        left = self.right
    elif self.right == NULL:
        left = self.left
    else:
        left = self.left
        cur = self.right
        while cur.right != NULL:
            cur = cur.right
        cur.right = self.left.right
    free_binary_tree_node(self)
    return left


cdef binary_tree_node *binary_tree_right_excise(binary_tree_node *self) noexcept:
    cdef binary_tree_node *right
    cdef binary_tree_node *cur
    if self.right == NULL:
        right = self.left
    elif self.left == NULL:
        right = self.right
    else:
        right = self.right
        cur = self.left
        while cur.left != NULL:
            cur = cur.left
        cur.left = self.right.left
    free_binary_tree_node(self)
    return right


cdef binary_tree_node *binary_tree_head_excise(binary_tree_node *self) noexcept:
    cdef binary_tree_node *cur
    cdef int right
    # We have a pointer we're about to free.  Chances are, we'll never
    # see this pointer again.  Thus, its least significant bit is
    # "random" enough to resist bias.
    right = (<int>self)&1
    if self.right == NULL:
        return self.left
    if self.left == NULL:
        return self.right
    if right:
        # move right branch to left, return left
        cur = self.left
        while cur.right != NULL:
            cur = cur.right
        cur.right = self.right
        cur = self.left
    else:
        # move left branch to right, return right
        cur = self.right
        while cur.left != NULL:
            cur = cur.left
        cur.left = self.left
        cur = self.right
    free_binary_tree_node(self)
    return cur


cdef int LIST_PREORDER, LIST_POSTORDER, LIST_INORDER, LIST_KEYS, LIST_VALUES
LIST_PREORDER  = 1
LIST_INORDER = 2
LIST_POSTORDER = 4
LIST_KEYS = 8
LIST_VALUES = 16

cdef object binary_tree_list(binary_tree_node *cur, int behavior):
    if behavior & LIST_KEYS:
        item = int(cur.key)
    else:
        item = <object>cur.value

    if behavior & LIST_PREORDER:
        arry = [item]
    else:
        arry = []

    if cur.left != NULL:
        arry.extend(binary_tree_list(cur.left, behavior))

    if behavior & LIST_INORDER:
        arry.append(item)

    if cur.right != NULL:
        arry.extend(binary_tree_list(cur.right, behavior))

    if behavior & LIST_POSTORDER:
        arry.append(item)

    return arry


cdef class BinaryTree:
    """
    A simple binary tree with integer keys.
    """
    def __cinit__(BinaryTree self):
        self.head = NULL

    def __dealloc__(BinaryTree self):
        """
        TESTS:

        We test that :issue:`18897` is fixed::

            sage: def test():
            ....:     from sage.rings.polynomial.polynomial_compiled import CompiledPolynomialFunction
            ....:     import gc
            ....:     from collections import Counter
            ....:     gc.collect()
            ....:     pre={id(c) for c in gc.get_objects()}
            ....:     L = [-1, 9, -22, 21, -8, 1]
            ....:     for _ in range(100):
            ....:         CompiledPolynomialFunction(L)  # this creates and deallocs a binary tree
            ....:     gc.collect()
            ....:     post=Counter(type(o) for o in gc.get_objects() if id(o) not in pre)
            ....:     return [(k,v) for (k,v) in post.items() if v>10]
            sage: test()   # indirect doctest
            []
        """
        binary_tree_dealloc(self.head)

    def insert(BinaryTree self, object key, object value=None):
        """
        Insert a key-value pair into the BinaryTree.

        Duplicate keys are ignored.

        The first parameter, key, should be an int, or coercible (one-to-one)
        into an int.

        EXAMPLES::

            sage: from sage.misc.binary_tree import BinaryTree
            sage: t = BinaryTree()
            sage: t.insert(1)
            sage: t.insert(0)
            sage: t.insert(2)
            sage: t.insert(0,1)
            sage: t.get(0)
            0
        """
        cdef int ckey
        if value is None:
            value = key
        ckey = int(key)
        if self.head is NULL:
            self.head = BinaryTreeNode(ckey, value)
        else:
            binary_tree_insert(self.head, ckey, value)

    def delete(BinaryTree self, int key):
        """
        Remove a the node corresponding to key, and return the value
        associated with it.

        EXAMPLES::

            sage: from sage.misc.binary_tree import BinaryTree
            sage: t = BinaryTree()
            sage: t.insert(3,3)
            sage: t.insert(1,1)
            sage: t.insert(2,2)
            sage: t.insert(0,0)
            sage: t.insert(5,5)
            sage: t.insert(6,6)
            sage: t.insert(4,4)
            sage: t.delete(0)
            0
            sage: t.delete(3)
            3
            sage: t.delete(5)
            5
            sage: t.delete(2)
            2
            sage: t.delete(6)
            6
            sage: t.delete(1)
            1
            sage: t.delete(0)
            sage: t.get_max()
            4
            sage: t.get_min()
            4
        """
        cdef object r
        if self.head == NULL:
            return None
        elif self.head.key == key:
            r = <object>self.head.value
            self.head = binary_tree_head_excise(self.head)
            return r
        else:
            return binary_tree_delete(self.head, key)

    def get(BinaryTree self, int key):
        """
        Return the value associated with the key given.

        EXAMPLES::

            sage: from sage.misc.binary_tree import BinaryTree
            sage: t = BinaryTree()
            sage: t.insert(0, Matrix([[0,0], [1,1]]))                                   # needs sage.modules
            sage: t.insert(0, 1)
            sage: t.get(0)                                                              # needs sage.modules
            [0 0]
            [1 1]
        """
        if self.head == NULL:
            return None
        return binary_tree_get(self.head, key)

    def contains(BinaryTree self, int key):
        """
        Return whether a node with the given key exists in the tree.

        EXAMPLES::

            sage: from sage.misc.binary_tree import BinaryTree
            sage: t = BinaryTree()
            sage: t.contains(1)
            False
            sage: t.insert(1, 1)
            sage: t.contains(1)
            True
        """
        if self.head == NULL:
            return False
        return binary_tree_get(self.head, key) is not None

    def get_max(BinaryTree self):
        """
        Return the value of the node with the maximal key value.
        """
        cdef binary_tree_node *cur
        if self.head == NULL:
            return None
        cur = self.head
        while cur.right != NULL:
            cur = cur.right
        return <object>cur.value

    def get_min(BinaryTree self):
        """
        Return the value of the node with the minimal key value.
        """
        cdef binary_tree_node *cur
        if self.head == NULL:
            return None
        cur = self.head
        while cur.left != NULL:
            cur = cur.left
        return <object>cur.value

    def pop_max(BinaryTree self):
        """
        Return the value of the node with the maximal key value,
        and remove that node from the tree.

        EXAMPLES::

            sage: from sage.misc.binary_tree import BinaryTree
            sage: t = BinaryTree()
            sage: t.insert(4,'e')
            sage: t.insert(2,'c')
            sage: t.insert(0,'a')
            sage: t.insert(1,'b')
            sage: t.insert(3,'d')
            sage: t.insert(5,'f')
            sage: while not t.is_empty():
            ....:     print(t.pop_max())
            f
            e
            d
            c
            b
            a
        """
        cdef binary_tree_node *cur
        cdef object max
        if self.head == NULL:
            return None
        if self.head.right == NULL:
            max = <object>self.head.value
            cur = self.head.left
            free_binary_tree_node(self.head)
            self.head = cur
            return max
        cur = self.head
        while cur.right.right != NULL:
            cur = cur.right
        max = <object>cur.right.value
        cur.right = binary_tree_right_excise(cur.right)
        return max

    def pop_min(BinaryTree self):
        """
        Return the value of the node with the minimal key value,
        and remove that node from the tree.

        EXAMPLES::

            sage: from sage.misc.binary_tree import BinaryTree
            sage: t = BinaryTree()
            sage: t.insert(4,'e')
            sage: t.insert(2,'c')
            sage: t.insert(0,'a')
            sage: t.insert(1,'b')
            sage: t.insert(3,'d')
            sage: t.insert(5,'f')
            sage: while not t.is_empty():
            ....:     print(t.pop_min())
            a
            b
            c
            d
            e
            f
        """
        cdef binary_tree_node *cur
        cdef object min
        if self.head == NULL:
            return None
        if self.head.left == NULL:
            min = <object>self.head.value
            cur = self.head.right
            free_binary_tree_node(self.head)
            self.head = cur
            return min
        cur = self.head
        while cur.left.left != NULL:
            cur = cur.left
        min = <object>cur.left.value
        cur.left = binary_tree_left_excise(cur.left)
        return min

    def is_empty(BinaryTree self):
        """
        Return whether the tree has no nodes.

        EXAMPLES::

            sage: from sage.misc.binary_tree import BinaryTree
            sage: t = BinaryTree()
            sage: t.is_empty()
            True
            sage: t.insert(0,0)
            sage: t.is_empty()
            False
        """
        return self.head == NULL

    def keys(BinaryTree self, order='inorder'):
        """
        Return the keys sorted according to "order" parameter.

        The order can be one of
        "inorder", "preorder", or "postorder"
        """
        if self.head == NULL:
            return []

        if order == "postorder":
            o = LIST_POSTORDER
        elif order == "inorder":
            o = LIST_INORDER
        else:
            o = LIST_PREORDER

        return binary_tree_list(self.head, LIST_KEYS + o)

    def values(BinaryTree self, order='inorder'):
        """
        Return the keys sorted according to "order" parameter.

        The order can be one of
        "inorder", "preorder", or "postorder"
        """
        if self.head == NULL:
            return []

        if order == "postorder":
            o = LIST_POSTORDER
        elif order == "inorder":
            o = LIST_INORDER
        else:
            o = LIST_PREORDER

        return binary_tree_list(self.head, LIST_VALUES + o)

    def _headkey_(BinaryTree self):
        """
        Used by the stress tester.  Don't think a user would care.

        Email tom if you care what the headkey is.
        """
        if self.head == NULL:
            return 0
        return self.head.key


class Test:
    def random(self):
        self.binary_tree()

    def binary_tree(self, values = 100, cycles = 100000):
        """
        Perform a sequence of random operations, given random inputs
        to stress test the binary tree structure.

        This was useful during development to find memory leaks /
        segfaults.  Cycles should be at least 100 times as large as
        values, or the delete, contains, and get methods won't hit
        very often.

        INPUT:

        - ``values`` -- number of possible values to use

        - ``cycles`` -- number of operations to perform

        TESTS::

            sage: sage.misc.binary_tree.Test().random()
        """
        from sage.misc.prandom import randint
        t = BinaryTree()
        for i in range(cycles):
            r = randint(0, 8)
            s = randint(0, values)
            if r == 1:
                t.insert(s)
            elif r == 2:
                t.delete(t._headkey_())
            elif r == 3:
                t.get(s)
            elif r == 4:
                t.contains(s)
            elif r == 5:
                t.get_max()
            elif r == 6:
                t.get_min()
            elif r == 7:
                t.pop_min()
            elif r == 8:
                t.pop_max()
            else:
                t.delete(s)
