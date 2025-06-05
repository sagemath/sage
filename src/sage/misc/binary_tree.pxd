cdef struct binary_tree_node:
    int key
    binary_tree_node *left
    binary_tree_node *right
    void *value

# cdef binary_tree_node *BinaryTreeNode(int, object)
# cdef void free_binary_tree_node(binary_tree_node *)
# cdef void binary_tree_dealloc(binary_tree_node *)
# cdef void binary_tree_insert(binary_tree_node *self, int, object)
# cdef object binary_tree_get(binary_tree_node *, int)
# cdef object binary_tree_delete(binary_tree_node *, int)
# cdef binary_tree_node *binary_tree_left_excise(binary_tree_node *)
# cdef binary_tree_node *binary_tree_right_excise(binary_tree_node *)
# cdef binary_tree_node *binary_tree_head_excise(binary_tree_node *)
# cdef object binary_tree_list(binary_tree_node *, int)


# cdef int LIST_PREORDER, LIST_POSTORDER, LIST_INORDER, LIST_KEYS, LIST_VALUES
# LIST_PREORDER  = 1
# LIST_INORDER = 2
# LIST_POSTORDER = 4
# LIST_KEYS = 8
# LIST_VALUES = 16


cdef class BinaryTree:
    cdef binary_tree_node *head
