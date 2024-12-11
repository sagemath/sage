/*
 * Pairing heap
 *
 * Implements a pairing heap data structure as described in [1]. See also [2]
 * for more details.
 *
 * This implementation is templated by the type TI of items and the type TV of
 * the value associated with an item. The type TI must be either a standard type
 * (int, size_t, etc.) or a type equipped with a has function as supported by
 * std::unordered_map. The top of the heap is the item with smallest value,
 * i.e., this is a min heap data structure. The number of items in the heap is
 * not fixed. It supports the following operations:
 *
 * - empty(): return true if the heap is empty, and false otherwise.
 *
 * - push(item, value): push an item to the heap with specified value.
 *
 * - top(): access the pair (item, value) at the top of the heap, i.e., with
 *   smallest value in time O(1).
 *   This operation assumes that the heap is not empty.
 *
 * - top_item(): access the item at the top of the heap in time O(1).
 *   This operation assumes that the heap is not empty.
 *
 * - top_value(): access the value of the item at the top of the heap in O(1).
 *   This operation assumes that the heap is not empty.
 *
 * - pop(): remove top item from the heap in amortize time O(log(n))
 *
 * - decrease(item, new_value): change the value associated with the item to the
 *   specified value ``new_value`` in time o(log(n)). The new value must be
 *   smaller than the previous one. Otherwise the structure of the heap is no
 *   longer guaranteed.
 *   If the item is not already in the heap, this method calls method ``push``.
 *
 * - contains(item): check whether specified item is in the heap in time O(1).
 *
 * - value(item): return the value associated with the item in the heap.
 *   This operation assumes that the item is already in the heap.
 *
 * References:
 *
 * [1] M. L. Fredman, R. Sedgewick, D. D. Sleator, and R. E. Tarjan.
 *     "The pairing heap: a new form of self-adjusting heap".
 *     Algorithmica. 1 (1): 111-129, 1986. doi:10.1007/BF01840439.
 *
 * [2] https://en.wikipedia.org/wiki/Pairing_heap
 *
 * Author:
 * - David Coudert <david.coudert@inria.fr>
 *
 */

#ifndef PAIRING_HEAP_H
#define PAIRING_HEAP_H

#include <unordered_map>
#include <stdexcept>


namespace pairing_heap {

  template<
    typename TV,  // type of values
    typename T    // type of the child class
    >
  struct PairingHeapNodeBase {
    public:

      bool operator<=(PairingHeapNodeBase const& other) const {
	return static_cast<T const*>(this)->le_implem(static_cast<T const&>(other));
      }

      // Pair list of heaps and return pointer to the top of resulting heap
      static T *_pair(T *p) {
	if (p == nullptr) {
	  return nullptr;
	}

	/*
	 * Move toward the end of the list, counting elements along the way.
	 * This is done in order to:
	 * - know whether the list has odd or even number of nodes
	 * - speed up going-back through the list
	 */
	size_t children = 1;
	T *it = p;
	while (it->next != nullptr) {
	  it = it->next;
	  children++;
	}

	T *result;

	if (children % 2 == 1) {
	  T *a = it;
	  it = it->prev;
	  a->prev = a->next = nullptr;
	  result = a;
	} else {
	  T *a = it;
	  T *b = it->prev;
	  it = it->prev->prev;
	  a->prev = a->next = b->prev = b->next = nullptr;
	  result = _merge(a, b);
	}

	for (size_t i = 0; i < (children - 1) / 2; i++) {
	  T *a = it;
	  T *b = it->prev;
	  it = it->prev->prev;
	  a->prev = a->next = b->prev = b->next = nullptr;
	  result = _merge(_merge(a, b), result);
	}

	return result;
      } // end _pair


      // Merge 2 heaps and return pointer to the top of resulting heap
      static T *_merge(T *a, T *b) {
	if (*a <= *b)  { // Use comparison method of PairingHeapNodeBase
	  _link(a, b);
	  return a;
	} else {
	  _link(b, a);
	  return b;
	}
      } // end _merge


      // Make b a child of a
      static void _link(T *a, T *b) {
	if (a->child != nullptr) {
	  b->next = a->child;
	  a->child->prev = b;
	}
	b->prev = a;
	a->child = b;
      } // end _link


      // Remove p from its parent children list
      static void _unlink(T *p) {
	if (p->prev->child == p) {
	  p->prev->child = p->next;
	} else {
	  p->prev->next = p->next;
	}
	if (p->next != nullptr) {
	  p->next->prev = p->prev;
	}
	p->prev = nullptr;
	p->next = nullptr;
      } // end _unlink


      TV value;  // value associated to the node
      T * prev;  // Previous sibling of the node or parent
      T * next;  // Next sibling of the node
      T * child; // First child of the node

    protected:
      // Only derived class can build a PairingHeapNodeBase
      explicit PairingHeapNodeBase(const TV &some_value)
	: value{some_value}, prev{nullptr}, next{nullptr}, child{nullptr} {
	}
  }; // end struct PairingHeapNodeBase


  template<
    typename TI,  // type of items stored in the node
    typename TV   // type of values associated with the stored item
                  // Assumes TV is a comparable type
    >
  class PairingHeapNode
    : public PairingHeapNodeBase<TV, PairingHeapNode<TI, TV>> {

  public:
    PairingHeapNode(TI const& some_item, TV const& some_value)
      : Base_(some_value), item(some_item) {
    }

    bool le_implem(PairingHeapNode const& other) const {
      return this->value <= other.value;
    }

    TI item;  // item contained in the node

  private:
    using Base_ = PairingHeapNodeBase<TV, PairingHeapNode<TI, TV>>;
  };


  class PairingHeapNodePy
    : public PairingHeapNodeBase<PyObject *, PairingHeapNodePy> {
  public:
    PairingHeapNodePy(PyObject *some_value)
      : Base_(some_value) {
    }

    bool le_implem(PairingHeapNodePy const& other) const {
      return PyObject_RichCompareBool(this->value, other.value, Py_LE);
    }

  private:
    using Base_ = PairingHeapNodeBase<PyObject *, PairingHeapNodePy>;
  };



  template<
    typename TI,  // type of items stored in the node
    typename TV   // type of values associated with the stored item
                  // Assume TV is a comparable type
    >
  class PairingHeap
  {
  public:
    using HeapNodeType = PairingHeapNode<TI, TV>;

    // Constructor
    explicit PairingHeap()
      : root(nullptr) {
    }

    // Copy constructor
    PairingHeap(PairingHeap<TI, TV> const *other)
      : root(nullptr) {
      for (auto const& it: other->nodes) {
	push(it.first, it.second->value);
      }
    }

    // Destructor
    virtual ~PairingHeap() {
      for (auto const& it: nodes) {
	delete it.second;
      }
    }

    // Return true if the heap is empty, else false
    bool empty() const {
      return root == nullptr;
    }

    // Return true if the heap is not empty, else false
    explicit operator bool() const {
      return root != nullptr;
    }

    // Insert an item into the heap with specified value (priority)
    void push(const TI &some_item, const TV &some_value) {
      if (nodes.find(some_item) != nodes.end()) {
	throw std::invalid_argument("item already in the heap");
      }
      PairingHeapNode<TI, TV> *p = new PairingHeapNode<TI, TV>(some_item, some_value);
      nodes[some_item] = p;
      root = root == nullptr ? p : HeapNodeType::_merge(root, p);
    }

    // Return the top pair (item, value) of the heap
    std::pair<TI, TV> top() const {
      if (root ==  nullptr) {
	throw std::domain_error("trying to access the top of an empty heap");
      }
      return std::make_pair(root->item, root->value);
    }

    // Return the top item of the heap
    TI top_item() const {
      if (root == nullptr) {
	throw std::domain_error("trying to access the top of an empty heap");
      }
      return root->item;
    }

    // Return the top value of the heap
    TV top_value() const {
      if (root == nullptr) {
	throw std::domain_error("trying to access the top of an empty heap");
      }
      return root->value;
    }

    // Remove the top element from the heap. Do nothing if empty
    void pop() {
      if (root != nullptr) {
	PairingHeapNode<TI, TV> *p = root->child;
	nodes.erase(root->item);
	delete root;
	root = HeapNodeType::_pair(p);
      }
    }

    // Decrease the value of specified item
    // If the item is not in the heap, push it
    void decrease(const TI &some_item, const TV &new_value) {
      if (contains(some_item)) {
	PairingHeapNode<TI, TV> *p = nodes[some_item];
	if (p->value <= new_value) {
	  throw std::invalid_argument("the new value must be less than the current value");
	}
	p->value = new_value;
	if (p->prev != nullptr) {
	  HeapNodeType::_unlink(p);
	  root = HeapNodeType::_merge(root, p);
	}
      } else {
	push(some_item, new_value);
      }
    }

    // Check if specified item is in the heap
    bool contains(TI const& some_item) const {
      return nodes.find(some_item) != nodes.end();
    }

    // Return the value associated with the item
    TV value(const TI &some_item) const {
      auto it = nodes.find(some_item);
      if (it == nodes.end()) {
	throw std::invalid_argument("the specified item is not in the heap");
      }
      return it->second->value;
    }

    // Return the number of items in the heap
    size_t size() const {
      return nodes.size();
    }

  private:

    // Pointer to the top of the heap
    PairingHeapNode<TI, TV> *root;

    // Map used to access stored items
    std::unordered_map<TI, PairingHeapNode<TI, TV> *> nodes;

  }; // end class PairingHeap

} // end namespace pairing_heap

#endif
