from typing import Optional, List

class binary_tree_node:
    key: int
    left: Optional['binary_tree_node']
    right: Optional['binary_tree_node']
    value: object

def BinaryTreeNode(key: int, value: object) -> binary_tree_node:
    ...

def free_binary_tree_node(node: binary_tree_node) -> None:
    ...

def binary_tree_dealloc(node: Optional[binary_tree_node]) -> None:
    ...

def binary_tree_insert(node: binary_tree_node, key: int, value: object) -> None:
    ...

def binary_tree_get(node: Optional[binary_tree_node], key: int) -> Optional[object]:
    ...

def binary_tree_delete(node: binary_tree_node, key: int) -> Optional[object]:
    ...

def binary_tree_left_excise(node: binary_tree_node) -> Optional[binary_tree_node]:
    ...

def binary_tree_right_excise(node: binary_tree_node) -> Optional[binary_tree_node]:
    ...

def binary_tree_head_excise(node: binary_tree_node) -> Optional[binary_tree_node]:
    ...

LIST_PREORDER: int
LIST_POSTORDER: int
LIST_INORDER: int
LIST_KEYS: int
LIST_VALUES: int

def binary_tree_list(node: binary_tree_node, behavior: int) -> List[object]:
    ...

class BinaryTree:
    head: Optional[binary_tree_node]

    def __cinit__(self) -> None:
        ...

    def __dealloc__(self) -> None:
        ...

    def insert(self, key: object, value: Optional[object] = None) -> None:
        ...

    def delete(self, key: int) -> Optional[object]:
        ...

    def get(self, key: int) -> Optional[object]:
        ...

    def contains(self, key: int) -> bool:
        ...

    def get_max(self) -> Optional[object]:
        ...

    def get_min(self) -> Optional[object]:
        ...

    def pop_max(self) -> Optional[object]:
        ...

    def pop_min(self) -> Optional[object]:
        ...

    def is_empty(self) -> bool:
        ...

    def keys(self, order: str = 'inorder') -> List[int]:
        ...

    def values(self, order: str = 'inorder') -> List[object]:
        ...

    def _headkey_(self) -> int:
        ...
