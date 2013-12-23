# -*- coding: utf-8 -*-
r"""
Posets

This module implements finite partially ordered sets. It defines:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :class:`FinitePoset` | A class for finite posets
    :class:`FinitePosets_n` | A class for finite posets up to isomorphism (i.e. unlabeled posets)
    :meth:`Poset` | Construct a finite poset from various forms of input data.
    :meth:`is_poset` | Tests whether a directed graph is acyclic and transitively reduced.

**List of Poset methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FinitePoset.antichains_iterator` | Returns an iterator over the antichains of the poset.
    :meth:`~FinitePoset.antichains` | Returns the antichains of the poset.
    :meth:`~FinitePoset.bottom` | Returns the bottom element of the poset, if it exists.
    :meth:`~FinitePoset.cardinality` | Returns the number of elements in the poset.
    :meth:`~FinitePoset.chains` | Returns all the chains of ``self``
    :meth:`~FinitePoset.chain_polytope` | Returns the chain polytope of the poset.
    :meth:`~FinitePoset.chain_polynomial` | Returns the chain polynomial of the poset.
    :meth:`~FinitePoset.closed_interval` | Returns a list of the elements `z` such that `x \le z \le y`.
    :meth:`~FinitePoset.compare_elements` | Compare `x` and `y` in the poset.
    :meth:`~FinitePoset.comparability_graph` | Returns the comparability graph of the poset.
    :meth:`~FinitePoset.cover_relations_iterator` | Returns an iterator for the cover relations of the poset.
    :meth:`~FinitePoset.cover_relations` | Returns the list of pairs [u,v] which are cover relations
    :meth:`~FinitePoset.covers` | Returns True if y covers x and False otherwise.
    :meth:`~FinitePoset.coxeter_transformation` | Returns the matrix of the Auslander-Reiten translation acting on the Grothendieck group of the derived category of modules
    :meth:`~FinitePoset.dual` | Returns the dual poset of the given poset.
    :meth:`~FinitePoset.evacuation` | Computes evacuation on the linear extension associated to the poset ``self``.
    :meth:`~FinitePoset.f_polynomial` | Returns the f-polynomial of a bounded poset.
    :meth:`~FinitePoset.flag_f_polynomial` | Returns the flag f-polynomial of a bounded and ranked poset.
    :meth:`~FinitePoset.flag_h_polynomial` | Returns the flag h-polynomial of a bounded and ranked poset.
    :meth:`~FinitePoset.frank_network` | Returns Frank's network (a DiGraph along with a cost function on its edges) associated to ``self``.
    :meth:`~FinitePoset.graphviz_string` | Returns a representation in the DOT language, ready to render in graphviz.
    :meth:`~FinitePoset.greene_shape` | Computes the Greene-Kleitman partition aka Greene shape of the poset ``self``.
    :meth:`~FinitePoset.h_polynomial` | Returns the h-polynomial of a bounded poset.
    :meth:`~FinitePoset.has_bottom` | Returns True if the poset has a unique minimal element.
    :meth:`~FinitePoset.hasse_diagram` | Returns the Hasse diagram of ``self`` as a Sage :class:`DiGraph`.
    :meth:`~FinitePoset.has_top` | Returns True if the poset contains a unique maximal element, and False otherwise.
    :meth:`~FinitePoset.incomparability_graph` | Returns the incomparability graph of the poset.
    :meth:`~FinitePoset.interval` | Returns a list of the elements `z` such that `x \le z \le y`.
    :meth:`~FinitePoset.is_bounded` | Returns True if the poset contains a unique maximal element and a unique minimal element, and False otherwise.
    :meth:`~FinitePoset.is_chain` | Returns True if the poset is totally ordered, and False otherwise.
    :meth:`~FinitePoset.is_EL_labelling` | Returns whether ``f`` is an EL labelling of ``self``
    :meth:`~FinitePoset.is_gequal` | Returns ``True`` if `x` is greater than or equal to `y` in the poset, and ``False`` otherwise.
    :meth:`~FinitePoset.is_graded` | Returns whether this poset is graded.
    :meth:`~FinitePoset.is_greater_than` | Returns ``True`` if `x` is greater than but not equal to `y` in the poset, and ``False`` otherwise.
    :meth:`~FinitePoset.is_isomorphic` | Returns True if both posets are isomorphic.
    :meth:`~FinitePoset.is_join_semilattice` | Returns True is the poset has a join operation, and False otherwise.
    :meth:`~FinitePoset.is_lequal` | Returns ``True`` if `x` is less than or equal to `y` in the poset, and ``False`` otherwise.
    :meth:`~FinitePoset.is_less_than` | Returns ``True`` if `x` is less than but not equal to `y` in the poset, and ``False`` otherwise.
    :meth:`~FinitePoset.is_linear_extension` | Returns whether ``l`` is a linear extension of ``self``
    :meth:`~FinitePoset.is_meet_semilattice` | Returns True if self has a meet operation, and False otherwise.
    :meth:`~FinitePoset.join_matrix` | Returns a matrix whose ``(i,j)`` entry is ``k``, where ``self.linear_extension()[k]`` is the join (least upper bound) of ``self.linear_extension()[i]`` and ``self.linear_extension()[j]``.
    :meth:`~FinitePoset.is_incomparable_chain_free` | Returns whether the poset is `(m+n)`-free.
    :meth:`~FinitePoset.is_ranked` | Returns whether this poset is ranked.
    :meth:`~FinitePoset.is_slender` | Returns whether the poset ``self`` is slender or not.
    :meth:`~FinitePoset.lequal_matrix` | Computes the matrix whose ``(i,j)`` entry is 1 if ``self.linear_extension()[i] < self.linear_extension()[j]`` and 0 otherwise
    :meth:`~FinitePoset.level_sets` | Returns a list l such that l[i+1] is the set of minimal elements of the poset obtained by removing the elements in l[0], l[1], ..., l[i].
    :meth:`~FinitePoset.linear_extension` | Returns a linear extension of this poset.
    :meth:`~FinitePoset.linear_extensions` | Returns the enumerated set of all the linear extensions of this poset
    :meth:`~FinitePoset.list` | List the elements of the poset. This just returns the result of :meth:`linear_extension`.
    :meth:`~FinitePoset.lower_covers_iterator` | Returns an iterator for the lower covers of the element y. An lower cover of y is an element x such that y x is a cover relation.
    :meth:`~FinitePoset.lower_covers` | Returns a list of lower covers of the element y. An lower cover of y is an element x such that y x is a cover relation.
    :meth:`~FinitePoset.maximal_chains` | Returns all maximal chains of this poset.  Each chain is listed in increasing order.
    :meth:`~FinitePoset.maximal_elements` | Returns a list of the maximal elements of the poset.
    :meth:`~FinitePoset.meet_matrix` | Returns a matrix whose ``(i,j)`` entry is ``k``, where ``self.linear_extension()[k]`` is the meet (greatest lower bound) of ``self.linear_extension()[i]`` and ``self.linear_extension()[j]``.
    :meth:`~FinitePoset.minimal_elements` | Returns a list of the minimal elements of the poset.
    :meth:`~FinitePoset.mobius_function_matrix` | Returns a matrix whose ``(i,j)`` entry is the value of the Mobius function evaluated at ``self.linear_extension()[i]`` and ``self.linear_extension()[j]``.
    :meth:`~FinitePoset.mobius_function` | Returns the value of the Mobius function of the poset on the elements x and y.
    :meth:`~FinitePoset.open_interval` | Returns a list of the elements `z` such that `x < z < y`. The order is that induced by the ordering in
    :meth:`~FinitePoset.order_complex` | Returns the order complex associated to this poset.
    :meth:`~FinitePoset.order_filter` | Returns the order filter generated by a list of elements.
    :meth:`~FinitePoset.order_ideal` | Returns the order ideal generated by a list of elements.
    :meth:`~FinitePoset.order_polynomial` | Returns the order polynomial of the poset.
    :meth:`~FinitePoset.order_polytope` | Returns the order polytope of the poset.
    :meth:`~FinitePoset.p_partition_enumerator` | Returns a `P`-partition enumerator of the poset.
    :meth:`~FinitePoset.plot` | Returns a Graphic object corresponding the Hasse diagram of the poset.
    :meth:`~FinitePoset.product` | Returns the cartesian product of ``self`` and ``other``.
    :meth:`~FinitePoset.promotion` | Computes the (extended) promotion on the linear extension of the poset ``self``
    :meth:`~FinitePoset.random_subposet` | Returns a random subposet that contains each element with probability p.
    :meth:`~FinitePoset.rank_function` | Returns a rank function of the poset, if it exists.
    :meth:`~FinitePoset.rank` | Returns the rank of an element, or the rank of the poset if element is None.
    :meth:`~FinitePoset.relabel` | Returns a copy of this poset with its elements relabelled
    :meth:`~FinitePoset.relations_iterator` | Returns an iterator for all the relations of the poset.
    :meth:`~FinitePoset.relations` | Returns a list of all relations of the poset.
    :meth:`~FinitePoset.show` | Shows the Graphics object corresponding the Hasse diagram of the poset.
    :meth:`~FinitePoset.subposet` | Returns the poset containing elements with partial order induced by that of self.
    :meth:`~FinitePoset.top` | Returns the top element of the poset, if it exists.
    :meth:`~FinitePoset.unwrap` | Unwraps an element of this poset
    :meth:`~FinitePoset.upper_covers_iterator` | Returns an iterator for the upper covers of the element y. An upper cover of y is an element x such that y x is a cover relation.
    :meth:`~FinitePoset.upper_covers` | Returns a list of upper covers of the element y. An upper cover of y is an element x such that y x is a cover relation.
    :meth:`~FinitePoset.with_linear_extension` | Returns a copy of ``self`` with a different default linear extension
    :meth:`~FinitePoset.zeta_polynomial` | Returns the zeta polynomial of the poset.

Classes and functions
---------------------
"""
#*****************************************************************************
#       Copyright (C) 2008 Peter Jipsen <jipsen@chapman.edu>,
#                          Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import random
import copy
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc_c import prod
from sage.misc.superseded import deprecated_function_alias
from sage.categories.category import Category
from sage.categories.sets_cat import Sets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.posets import Posets
from sage.categories.finite_posets import FinitePosets
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring import polygen
from sage.graphs.digraph import DiGraph
from sage.graphs.digraph_generators import digraphs
from sage.combinat.posets.hasse_diagram import HasseDiagram
from sage.combinat.posets.elements import PosetElement
from sage.combinat.combinatorial_map import combinatorial_map


def Poset(data=None, element_labels=None, cover_relations=False, linear_extension=False, category = None, facade = None, key = None):
    r"""
    Construct a finite poset from various forms of input data.

    INPUT:

    - ``data`` -- different input are accepted by this constructor:

        1. A two-element list or tuple `(E, R)`, where `E` is a
           collection of elements of the poset and `R` is a collection
           of relations `x<=y`, each represented as a two-element
           lists/tuples/iterables such as [x,y]. The poset is then the
           transitive closure of the provided relations. If
           ``cover_relations=True``, then `R` is assumed to contain
           exactly the cover relations of the poset. If `E` is empty,
           then `E` is taken to be the set of elements appearing in
           the relations `R`.

        2. A two-element list or tuple `(E, f)`, where `E` is the set
           of elements of the poset and `f` is a function such that,
           for any pair `x,y` of elements of `E`, `f(x,y)` returns
           whether `x <= y`. If ``cover_relations=True``, then
           `f(x,y)` should return whether `x` is covered by `y`.

        3. A dictionary, list or tuple of upper covers: ``data[x]`` is
           a list of the elements that cover the element `x` in the
           poset.

           .. WARNING::

              If data is a list or tuple of length `2`, then it is
              handled by the above case..

        4. An acyclic, loop-free and multi-edge free ``DiGraph``. If
           ``cover_relations`` is ``True``, then the edges of the
           digraph are assumed to correspond to the cover relations of
           the poset. Otherwise, the cover relations are computed.

        5. A previously constructed poset (the poset itself is returned).

    - ``element_labels`` -- (default: None); an optional list or
      dictionary of objects that label the poset elements.

    - ``cover_relations`` -- a boolean (default: False); whether the
      data can be assumed to describe a directed acyclic graph whose
      arrows are cover relations; otherwise, the cover relations are
      first computed.

    - ``linear_extension`` -- a boolean (default: False); whether to
      use the provided list of elements as default linear extension
      for the poset; otherwise a linear extension is computed.

    - ``facade`` -- a boolean or ``None`` (default); whether the
      :meth:`Poset`'s elements should be wrapped to make them aware of the Poset
      they belong to.

      * If ``facade = True``, the :meth:`Poset`'s elements are exactly those
        given as input.

      * If ``facade = False``, the :meth:`Poset`'s elements will become
        :class:`~sage.combinat.posets.posets.PosetElement` objects.

      * If ``facade = None`` (default) the expected behaviour is the behaviour
        of ``facade = True``, unless the opposite can be deduced from the
        context (i.e. for instance if a :meth:`Poset` is built from another
        :meth:`Poset`, itself built with ``facade = False``)

    OUTPUT:

        ``FinitePoset`` -- an instance of the :class:`FinitePoset`` class.

    If ``category`` is specified, then the poset is created in this
    category instead of :class:`FinitePosets`.

    .. seealso:: :class:`Posets`, :class:`~sage.categories.posets.Posets`, :class:`FinitePosets`

    EXAMPLES:

    1. Elements and cover relations::

          sage: elms = [1,2,3,4,5,6,7]
          sage: rels = [[1,2],[3,4],[4,5],[2,5]]
          sage: Poset((elms, rels), cover_relations = True, facade = False)
          Finite poset containing 7 elements

       Elements and non-cover relations::

          sage: elms = [1,2,3,4]
          sage: rels = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
          sage: P = Poset( [elms,rels] ,cover_relations=False); P
          Finite poset containing 4 elements
          sage: P.cover_relations()
          [[1, 2], [2, 3], [3, 4]]

    2. Elements and function: the standard permutations of [1, 2, 3, 4]
       with the Bruhat order::

          sage: elms = Permutations(4)
          sage: fcn = lambda p,q : p.bruhat_lequal(q)
          sage: Poset((elms, fcn))
          Finite poset containing 24 elements

       With a function that identifies the cover relations: the set
       partitions of {1, 2, 3} ordered by refinement::

          sage: elms = SetPartitions(3)
          sage: def fcn(A, B):
          ...     if len(A) != len(B)+1:
          ...         return False
          ...     for a in A:
          ...         if not any(set(a).issubset(b) for b in B):
          ...             return False
          ...     return True
          sage: Poset((elms, fcn), cover_relations=True)
          Finite poset containing 5 elements

    3. A dictionary of upper covers::

          sage: Poset({'a':['b','c'], 'b':['d'], 'c':['d'], 'd':[]})
          Finite poset containing 4 elements

       A list of upper covers::

          sage: Poset([[1,2],[4],[3],[4],[]])
          Finite poset containing 5 elements

       A list of upper covers and a dictionary of labels::

          sage: elm_labs = {0:"a",1:"b",2:"c",3:"d",4:"e"}
          sage: P = Poset([[1,2],[4],[3],[4],[]],elm_labs, facade = False)
          sage: P.list()
          [a, b, c, d, e]

       .. warning::

         The special case where the argument data is a list or tuple of
         length 2 is handled by the above cases. So you cannot use this
         method to input a 2-element poset.

    4. An acyclic DiGraph.

       ::

          sage: dag = DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
          sage: Poset(dag)
          Finite poset containing 6 elements

       Any directed acyclic graph without loops or multiple edges, as long
       as ``cover_relations=False``::

          sage: dig = DiGraph({0:[2,3], 1:[3,4,5], 2:[5], 3:[5], 4:[5]})
          sage: dig.allows_multiple_edges()
          False
          sage: dig.allows_loops()
          False
          sage: dig.transitive_reduction() == dig
          False
          sage: Poset(dig, cover_relations=False)
          Finite poset containing 6 elements
          sage: Poset(dig, cover_relations=True)
          Traceback (most recent call last):
          ...
          ValueError: Hasse diagram is not transitively reduced.

    .. rubric:: Default Linear extension

    Every poset `P` obtained with ``Poset`` comes equipped with a
    default linear extension, which is also used for enumerating
    its elements. By default, this linear extension is computed,
    and has no particular significance::

        sage: P = Poset((divisors(12), attrcall("divides")))
        sage: P.list()
        [1, 2, 4, 3, 6, 12]
        sage: P.linear_extension()
        [1, 2, 4, 3, 6, 12]

    You may enforce a specific linear extension using the
    ``linear_extension`` option::

        sage: P = Poset((divisors(12), attrcall("divides")), linear_extension=True)
        sage: P.list()
        [1, 2, 3, 4, 6, 12]
        sage: P.linear_extension()
        [1, 2, 3, 4, 6, 12]

    Depending on popular request, ``Poset`` might eventually get
    modified to always use the provided list of elements as
    default linear extension, when it is one.

    .. seealso:: :meth:`FinitePoset.linear_extensions`

    .. rubric:: Facade posets

    When ``facade = False``, the elements of a poset are wrapped so as to make
    them aware that they belong to that poset::

        sage: P = Poset(DiGraph({'d':['c','b'],'c':['a'],'b':['a']}), facade = False)
        sage: d,c,b,a = list(P)
        sage: a.parent() is P
        True

    This allows for comparing elements according to `P`::

        sage: c < a
        True

    However, this may have surprising effects::

        sage: my_elements = ['a','b','c','d']
        sage: any(x in my_elements for x in P)
        False

    and can be anoying when one wants to manipulate the elements of
    the poset::

        sage: a + b
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand type(s) for +: 'FinitePoset_with_category.element_class' and 'FinitePoset_with_category.element_class'
        sage: a.element + b.element
        'ac'

    By default, facade posets are constructed instead::

        sage: P = Poset(DiGraph({'d':['c','b'],'c':['a'],'b':['a']}))

    In this example, the elements of the poset remain plain strings::

        sage: d,c,b,a = list(P)
        sage: type(a)
        <type 'str'>

    Of course, those strings are not aware of `P`. So to compare two
    such strings, one needs to query `P`::

        sage: a < b
        True
        sage: P.lt(a,b)
        False

    which models the usual mathematical notation `a <_P b`.

    Most operations seem to still work, but at this point there is no
    guarantee whatsoever::

        sage: P.list()
        ['d', 'b', 'c', 'a']
        sage: P.principal_order_ideal('a')
        ['d', 'b', 'c', 'a']
        sage: P.principal_order_ideal('b')
        ['d', 'b']
        sage: P.principal_order_ideal('d')
        ['d']
        sage: TestSuite(P).run()

    .. warning::

        :class:`DiGraph` is used to construct the poset, and the
        vertices of a :class:`DiGraph` are converted to plain Python
        :class:`int`'s if they are :class:`Integer`'s::

            sage: G = DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: type(G.vertices()[0])
            <type 'int'>

        This is worked around by systematically converting back the
        vertices of a poset to :class:`Integer`'s if they are
        :class:`int`'s::

            sage: P = Poset((divisors(15), attrcall("divides")), facade = False)
            sage: type(P.an_element().element)
            <type 'sage.rings.integer.Integer'>

            sage: P = Poset((divisors(15), attrcall("divides")), facade=True)
            sage: type(P.an_element())
            <type 'sage.rings.integer.Integer'>

        This may be abusive::

            sage: P = Poset((range(5), operator.le), facade = True)
            sage: P.an_element().parent()
            Integer Ring

    .. rubric:: Unique representation

    As most parents, :class:`Poset` have unique representation (see
    :class:`UniqueRepresentation`. Namely if two posets are created
    from two equal data, then they are not only equal but actually
    identical::

        sage: data1 = [[1,2],[3],[3]]
        sage: data2 = [[1,2],[3],[3]]
        sage: P1 = Poset(data1)
        sage: P2 = Poset(data2)
        sage: P1 == P2
        True
        sage: P1 is P2
        True

    In situations where this behaviour is not desired, one can use the
    ``key`` option::

        sage: P1 = Poset(data1, key = "foo")
        sage: P2 = Poset(data2, key = "bar")
        sage: P1 is P2
        False
        sage: P1 == P2
        False

    ``key`` can be any hashable value and is passed down to
    :class:`UniqueRepresentation`. It is otherwise ignored by the
    poset constructor.

    TESTS::

        sage: P = Poset([[1,2],[3],[3]])
        sage: type(hash(P))
        <type 'int'>

    Bad input::

        sage: Poset([1,2,3], lambda x,y : x<y)
        Traceback (most recent call last):
        ...
        ValueError: elements_label should be a dict or a list if different from None. (Did you intend data to be equal to a pair ?)
    """
    # Avoiding some errors from the user when data should be a pair
    if (element_labels is not None and
        not isinstance(element_labels, dict) and
        not isinstance(element_labels, list)):
        raise ValueError("elements_label should be a dict or a list if "+
                         "different from None. (Did you intend data to be "+
                         "equal to a pair ?)")

    #Convert data to a DiGraph
    elements = None
    D = {}
    if isinstance(data, FinitePoset):
        if element_labels is None and category is None and facade is None:
            return data
        else:
            return FinitePoset(data, data._elements, category = category, facade = facade)
    elif data is None: # type 0
        D = DiGraph()
    elif isinstance(data, DiGraph): # type 4
        D = copy.deepcopy(data)
    elif isinstance(data, dict): # type 3: dictionary of upper covers
        D = DiGraph(data)
    elif isinstance(data,(list,tuple)): # types 1, 2, 3 (list/tuple)
        if len(data) == 2: # types 1 or 2
            if callable(data[1]): # type 2
                elements, function = data
                relations = []
                for x in elements:
                    for y in elements:
                        if function(x,y) is True:
                            relations.append([x,y])
            else: # type 1
                elements, relations = data
                # check that relations are relations
                for r in relations:
                    try:
                        u, v = r
                    except ValueError:
                        raise TypeError, "not a list of relations"
            D = DiGraph()
            D.add_vertices(elements)
            D.add_edges(relations)
        elif len(data) > 2:
            # type 3, list/tuple of upper covers
            D = DiGraph(dict([[Integer(i),data[i]] for i in range(len(data))]))
        else:
            raise ValueError, "not valid poset data."

    # DEBUG: At this point D should be a DiGraph.
    if not isinstance(D,DiGraph):
        raise TypeError, "BUG: D should be a digraph."

    # Determine cover relations, if necessary.
    if cover_relations is False:
        D = D.transitive_reduction()

    # Check that the digraph does not contain loops, multiple edges
    # and is transitively reduced.
    if D.has_loops():
        raise ValueError, "Hasse diagram contains loops."
    elif D.has_multiple_edges():
        raise ValueError, "Hasse diagram contains multiple edges."
    elif cover_relations is True and not D.is_transitively_reduced():
        raise ValueError, "Hasse diagram is not transitively reduced."

    if linear_extension and elements is not None:
        lin_ext = list(elements)
    else:
        # Compute a linear extension of the poset (a topological sort).
        try:
            lin_ext = D.topological_sort()
        except StandardError:
            raise ValueError, "Hasse diagram contains cycles."

    # Relabel using the linear_extension.
    # So range(len(D)) becomes a linear extension of the poset.
    rdict = dict([[lin_ext[i],i] for i in range(len(lin_ext))])
    D.relabel(rdict)

    # Set element labels.
    if element_labels is None:
        elements = lin_ext
        # Work around the fact that, currently, when a DiGraph is
        # created with Integer's as vertices, those vertices are
        # converted to plain int's. This is a bit abusive.
        elements = [ Integer(i) if isinstance(i,int) else i for i in elements ]
    else:
        elements = [element_labels[z] for z in lin_ext]

    return FinitePoset(D,elements, category = category, facade = facade, key = key)

class FinitePoset(UniqueRepresentation, Parent):
    r"""
    Constructs a (finite) `n`-element poset from a set of elements and a
    directed acyclic graph or poset.

    INPUT:

    - ``hasse_diagram`` -- an instance of this class (``FinitePoset``),
      or a digraph that is transitively-reduced, acyclic, loop-free,
      multiedge-free, and with vertices indexed by ``range(n)``. We also
      assume that ``range(n)`` is a linear extension of the poset.

    - ``elements`` - an optional list of elements, with ``element[i]``
      corresponding to vertex ``i``. If ``elements`` is ``None``, then it is
      set to be the vertex set of the digraph.

    - ``category`` -- :class:`FinitePosets`, or a subcategory thereof.

    - ``facade`` -- a boolean or ``None`` (default); whether the
      :class:`~sage.combinat.posets.posets.FinitePoset`'s elements should be
      wrapped to make them aware of the Poset they belong to.

      * If ``facade = True``, the
        :class:`~sage.combinat.posets.posets.FinitePoset`'s elements are exactly
        those given as input.

      * If ``facade = False``, the
        :class:`~sage.combinat.posets.posets.FinitePoset`'s elements will become
        :class:`~sage.combinat.posets.posets.PosetElement` objects.

      * If ``facade = None`` (default) the expected behaviour is the behaviour
        of ``facade = True``, unless the opposite can be deduced from the
        context (i.e. for instance if a
        :class:`~sage.combinat.posets.posets.FinitePoset` is built from another
        :class:`~sage.combinat.posets.posets.FinitePoset`, itself built with
        ``facade = False``)

    - ``key`` -- any hashable value (default: ``None``).

    EXAMPLES::

        sage: uc = [[2,3], [], [1], [1], [1], [3,4]]
        sage: from sage.combinat.posets.posets import FinitePoset
        sage: P = FinitePoset(DiGraph(dict([[i,uc[i]] for i in range(len(uc))])), facade = False); P
        Finite poset containing 6 elements
        sage: P.cover_relations()
        [[0, 2], [0, 3], [2, 1], [3, 1], [4, 1], [5, 3], [5, 4]]
        sage: TestSuite(P).run()
        sage: P.category()
        Category of finite posets
        sage: P.__class__
        <class 'sage.combinat.posets.posets.FinitePoset_with_category'>

        sage: Q = sage.combinat.posets.posets.FinitePoset(P, facade = False); Q
        Finite poset containing 6 elements

        sage: Q is P
        True

    We keep the same underlying hasse diagram, but change the elements::

        sage: Q = sage.combinat.posets.posets.FinitePoset(P, elements=[1,2,3,4,5,6], facade = False); Q
        Finite poset containing 6 elements
        sage: Q.cover_relations()
        [[1, 3], [1, 4], [3, 2], [4, 2], [5, 2], [6, 4], [6, 5]]

    We test the facade argument::

        sage: P = Poset(DiGraph({'a':['b'],'b':['c'],'c':['d']}), facade = False)
        sage: P.category()
        Category of finite posets
        sage: parent(P[0]) is P
        True

        sage: Q = Poset(DiGraph({'a':['b'],'b':['c'],'c':['d']}), facade = True)
        sage: Q.category()
        Category of facade finite posets
        sage: parent(Q[0]) is str
        True
        sage: TestSuite(Q).run(skip = ['_test_an_element']) # is_parent_of is not yet implemented

    Changing a non facade poset to a facade poset::

        sage: PQ = Poset(P, facade = True)
        sage: PQ.category()
        Category of facade finite posets
        sage: parent(PQ[0]) is str
        True
        sage: PQ is Q
        True

    Changing a facade poset to a non facade poset::

        sage: QP = Poset(Q, facade = False)
        sage: QP.category()
        Category of finite posets
        sage: parent(QP[0]) is QP
        True

    .. note::

       A class that inherits from this class needs to define
       ``Element``. This is the class of the elements that the inheriting
       class contains. For example, for this class, ``FinitePoset``,
       ``Element`` is ``PosetElement``.  It can also define ``_dual_class`` which
       is the class of dual posets of this
       class. E.g. ``FiniteMeetSemilattice._dual_class`` is
       ``FiniteJoinSemilattice``.

    TESTS:

    Equality is derived from :class:`UniqueRepresentation`. We check that this
    gives consistent results::

        sage: P = Poset([[1,2],[3],[3]])
        sage: P == P
        True
        sage: Q = Poset([[1,2],[],[1]])
        sage: Q == P
        False
        sage: p1, p2 = Posets(2).list()
        sage: p2 == p1, p1 != p2
        (False, True)
        sage: [[p1.__eq__(p2) for p1 in Posets(2)] for p2 in Posets(2)]
        [[True, False], [False, True]]
        sage: [[p2.__eq__(p1) for p1 in Posets(2)] for p2 in Posets(2)]
        [[True, False], [False, True]]
        sage: [[p2 == p1 for p1 in Posets(3)] for p2 in Posets(3)]
        [[True, False, False, False, False], [False, True, False, False, False], [False, False, True, False, False], [False, False, False, True, False], [False, False, False, False, True]]

        sage: [[p1.__ne__(p2) for p1 in Posets(2)] for p2 in Posets(2)]
        [[False, True], [True, False]]
        sage: P = Poset([[1,2,4],[3],[3]])
        sage: Q = Poset([[1,2],[],[1],[4]])
        sage: P != Q
        True
        sage: P != P
        False
        sage: Q != Q
        False
        sage: [[p1.__ne__(p2) for p1 in Posets(2)] for p2 in Posets(2)]
        [[False, True], [True, False]]
    """

    @staticmethod
    def __classcall__(cls, hasse_diagram, elements = None, category = None, facade = None, key = None):
        """
        Normalizes the arguments passed to the constructor

        TESTS::

            sage: P = sage.combinat.posets.posets.FinitePoset(DiGraph())
            sage: type(P)
            <class 'sage.combinat.posets.posets.FinitePoset_with_category'>
            sage: TestSuite(P).run()

        See also the extensive tests in the class documentation
        """
        assert isinstance(hasse_diagram, (FinitePoset, DiGraph))
        if isinstance(hasse_diagram, FinitePoset):
            if elements is None:
                elements = hasse_diagram._elements
            if category is None:
                category = hasse_diagram.category()
                if facade is False and category.is_subcategory(Sets().Facades()):
                    # We need to remove Sets().Facades() from the category
                    # This is fragile ...
                    from sage.categories.category import JoinCategory
                    assert isinstance(category, JoinCategory)
                    categories = list(category.super_categories())
                    categories.remove(Sets().Facades())
                    category = Category.join(categories)
            if facade is None:
                facade = hasse_diagram in Sets().Facades()
            hasse_diagram = hasse_diagram._hasse_diagram
        else:
            hasse_diagram = HasseDiagram(hasse_diagram)
            if elements is None:
                elements = hasse_diagram.vertices()
            if facade is None:
                facade = True
        elements = tuple(elements)
        category = Category.join([FinitePosets().or_subcategory(category), FiniteEnumeratedSets()])
        return super(FinitePoset, cls).__classcall__(cls, hasse_diagram = hasse_diagram, elements = elements,
                                                     category = category, facade = facade, key = key)

    def __init__(self, hasse_diagram, elements, category, facade, key):
        """
        EXAMPLES::

            sage: P = Poset(DiGraph({'a':['b'],'b':['c'],'c':['d']}), facade = False)
            sage: type(P)
            <class 'sage.combinat.posets.posets.FinitePoset_with_category'>

        The internal data structure currently consists of:

        - the Hasse diagram of the poset, represented by a DiGraph
          with vertices labelled 0,...,n-1 according to a linear
          extension of the poset (that is if `i \mapsto j` is an edge
          then `i<j`), together with some extra methods (see
          :class:`sage.combinat.posets.hasse_diagram.HasseDiagram`)::

            sage: P._hasse_diagram
            Hasse diagram of a poset containing 4 elements
            sage: P._hasse_diagram.cover_relations()
            [(0, 1), (1, 2), (2, 3)]

        - a tuple of the original elements, not wrapped as elements of
          ``self`` (but see also ``P._list``)::

            sage: P._elements
            ('a', 'b', 'c', 'd')

          ``P._elements[i]`` gives the element of ``P`` corresponding
          to the vertex ``i``

        - a dictionary mapping back elements to vertices::

            sage: P._element_to_vertex_dict
            {'a': 0, 'c': 2, 'b': 1, 'd': 3}

        - and a boolean stating whether the poset is a facade poset::

            sage: P._is_facade
            False

        This internal data structure is subject to change at any
        point. Do not break encapsulation!

        TESTS::

            sage: TestSuite(P).run()

        See also the extensive tests in the class documentation.
        """
        Parent.__init__(self, category = category, facade = facade)
        self._hasse_diagram = hasse_diagram
        self._elements = elements
        self._element_to_vertex_dict = dict( (elements[i], i) for i in range(len(elements)) )
        self._is_facade = facade

    @lazy_attribute
    def _list(self):
        """
        The list of the elements of ``self``, each wrapped to have
        ``self`` as parent

        EXAMPLES::

            sage: P = Poset(DiGraph({'a':['b'],'b':['c'],'c':['d']}), facade = False)
            sage: L = P._list; L
            (a, b, c, d)
            sage: type(L[0])
            <class 'sage.combinat.posets.elements.FinitePoset_with_category.element_class'>
            sage: L[0].parent() is P
            True

        Constructing them once for all makes future conversions
        between the vertex id and the element faster. This also
        ensures unique representation of the elements of this poset,
        which could be used to later speed up certain operations
        (equality test, ...)
        """
        if self._is_facade:
            return self._elements
        else:
            return tuple(self.element_class(self, self._elements[vertex], vertex)
                         for vertex in range(len(self._elements)))

    # This defines the type (class) of elements of poset.
    Element = PosetElement

    def _element_to_vertex(self, element):
        """
        Given an element of the poset (wrapped or not), returns the
        corresponding vertex of the Hasse diagram.

        EXAMPLES::

            sage: P = Poset((divisors(15), attrcall("divides")))
            sage: list(P)
            [1, 3, 5, 15]
            sage: x = P(5)
            sage: P._element_to_vertex(x)
            2

        The same with a non-wrapped element of `P`::

            sage: P._element_to_vertex(5)
            2

        TESTS::

            sage: P = Poset((divisors(15), attrcall("divides")), facade = False)

        Testing for wrapped elements::

            sage: all(P._vertex_to_element(P._element_to_vertex(x)) is x for x in P)
            True

        Testing for non-wrapped elements::

            sage: all(P._vertex_to_element(P._element_to_vertex(x)) is P(x) for x in divisors(15))
            True

        Testing for non-wrapped elements for a facade poset::

            sage: P = Poset((divisors(15), attrcall("divides")), facade = True)
            sage: all(P._vertex_to_element(P._element_to_vertex(x)) is x for x in P)
            True
        """
        if isinstance(element, self.element_class) and element.parent() is self:
            return element.vertex
        else:
            try:
                return self._element_to_vertex_dict[element]
            except KeyError:
                raise ValueError, "element (=%s) not in poset"%element

    def _vertex_to_element(self, vertex):
        """
        Return the element of ``self`` corresponding to the vertex
        ``vertex`` of the Hasse diagram.

        It is wrapped if ``self`` is not a facade poset.

        EXAMPLES::

            sage: P = Poset((divisors(15), attrcall("divides")), facade = False)
            sage: x = P._vertex_to_element(2)
            sage: x
            5
            sage: x.parent() is P
            True

            sage: P = Poset((divisors(15), attrcall("divides")), facade = True)
            sage: x = P._vertex_to_element(2)
            sage: x
            5
            sage: x.parent() is ZZ
            True
        """
        return self._list[vertex]

    def unwrap(self, element):
        """
        Return the element ``element`` of the poset ``self`` in
        unwrapped form.

        INPUT:

        - ``element`` -- an element of ``self``

        EXAMPLES::

            sage: P = Poset((divisors(15), attrcall("divides")), facade = False)
            sage: x = P.an_element(); x
            1
            sage: x.parent()
            Finite poset containing 4 elements
            sage: P.unwrap(x)
            1
            sage: P.unwrap(x).parent()
            Integer Ring

        For a non facade poset, this is equivalent to using the
        ``.element`` attribute::

            sage: P.unwrap(x) is x.element
            True

        For a facade poset, this does nothing::

            sage: P = Poset((divisors(15), attrcall("divides")), facade=True)
            sage: x = P.an_element()
            sage: P.unwrap(x) is x
            True

        This method is useful in code where we don't know if ``P`` is
        a facade poset or not.
        """
        if self._is_facade:
            return element
        else:
            return element.element

    def __contains__(self, x):
        r"""
        Returns True if x is an element of the poset.

        TESTS::

            sage: from sage.combinat.posets.posets import FinitePoset
            sage: P5 = FinitePoset(DiGraph({(5,):[(4,1),(3,2)], \
                    (4,1):[(3,1,1),(2,2,1)], \
                    (3,2):[(3,1,1),(2,2,1)], \
                    (3,1,1):[(2,1,1,1)], \
                    (2,2,1):[(2,1,1,1)], \
                    (2,1,1,1):[(1,1,1,1,1)], \
                    (1,1,1,1,1):[]}))
            sage: x = P5.list()[3]
            sage: x in P5
            True

        For the sake of speed, an element with the right class and
        parent is assumed to be in this parent. This can possibly be
        counterfeited by feeding garbage to the constructor::

            sage: x = P5.element_class(P5, "a", 5)
            sage: x in P5
            True
        """
        if isinstance(x, self.element_class):
            return x.parent() is self
        return x in self._element_to_vertex_dict

    is_parent_of = __contains__

    def _element_constructor_(self, element):
        """
        Constructs an element of ``self``

        EXAMPLES::

            sage: from sage.combinat.posets.posets import FinitePoset
            sage: P = FinitePoset(DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]}), facade = False)
            sage: P(5)
            5
            sage: Q = FinitePoset(DiGraph({5:[2,3], 1:[3,4], 2:[0], 3:[0], 4:[0]}), facade = False)
            sage: Q(5)
            5

        Accessing the n-th element of ``self`` as ``P(i)`` is deprecated::

            sage: P(5) == P(-1)
            doctest:...: DeprecationWarning: Accessing the i-th element of a poset as P(i) is deprecated. Please use P[i]
            See http://trac.sagemath.org/13109 for details.
            True
            sage: Q(5) == Q(-1)
            True
            sage: R = FinitePoset(DiGraph({'a':['b','c'], 'b':['d'], 'c':['d'], 'd':[]}), facade = False)
            sage: R(0)
            a
            sage: R('a') == R(0)
            True
            sage: R('d') == R(-1)
            True

        Please use instead ``P[i]``::

            sage: R('a') == R[0]
            True
            sage: R('d') == R[-1]
            True

        TESTS::

            sage: P = Poset(DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]}), facade = False)
            sage: all(P(x) is x for x in P)
            True
            sage: P = Poset((divisors(15), attrcall("divides")), facade = True)
            sage: all(P(x) is x for x in P)
            True
        """
        try:
            return self._list[self._element_to_vertex_dict[element]]
        except KeyError:
            if isinstance(element,Integer):
                from sage.misc.superseded import deprecation
                deprecation(13109, "Accessing the i-th element of a poset as P(i) is deprecated. Please use P[i]")
                if element > -1:
                    return self.element_class(self, \
                        self._elements[element], element)
                else:
                    return self.element_class(self, \
                        self._elements[element], self.cardinality()+element)
            else:
                raise ValueError, "%s is not an element of this poset"%type(element)

    def __call__(self, element):
        """
        Creates elements of this poset

        This overrides the generic call method for all parents
        :meth:`Parent.__call__`, as a work around to allow for facade
        posets over plain Python objects (e.g. Python's
        int's). Indeed, the default __call__ method expects the result
        of :meth:`_element_constructor_` to be a Sage element (see
        :meth:`sage.structure.coerce_maps.DefaultConvertMap_unique._call_`)::

            sage: P = Poset(DiGraph({'d':['c','b'],'c':['a'],'b':['a']}),
            ...             facade = True)
            sage: P('a')              # indirect doctest
            'a'
            sage: TestSuite(P).run()
            sage: P = Poset(((False, True), operator.eq), facade = True)
            sage: P(True)
            1
        """
        if self._is_facade and element in self._element_to_vertex_dict:
            return element
        return super(FinitePoset, self).__call__(element)

    def hasse_diagram(self, wrapped = True):
        r"""
        Return the Hasse diagram of ``self`` as a Sage :class:`DiGraph`. If
        ``dot2tex`` is installed, then this sets the Hasse diagram's latex
        options to use the ``dot2tex`` formatting.

        .. TODO::

            Should the vertices of the diagram have the poset as parent?

        EXAMPLES::

            sage: Q = Poset({5:[2,3], 1:[3,4], 2:[0], 3:[0], 4:[0]}, facade = False)
            sage: Q.hasse_diagram()
            Digraph on 6 vertices

            sage: P = Poset({'a':['b'],'b':['d'],'c':['d'],'d':['f'],'e':['f'],'f':[]}, facade = False)
            sage: H = P.hasse_diagram()
            sage: P.cover_relations()
            [[e, f], [c, d], [a, b], [b, d], [d, f]]
            sage: H.edges()
            [(a, b, None), (c, d, None), (b, d, None), (e, f, None), (d, f, None)]

            sage: P = Poset((divisors(15), attrcall("divides")), facade = False)
            sage: H = P.hasse_diagram()
            sage: H.vertices()
            [1, 5, 3, 15]
            sage: H.edges()
            [(1, 3, None), (1, 5, None), (5, 15, None), (3, 15, None)]
            sage: H.set_latex_options(format = "dot2tex")   # optional - dot2tex
            sage: view(H, tight_page=True) # optional - dot2tex
        """
        G = DiGraph(self._hasse_diagram).relabel(self._list, inplace = False)
        from sage.graphs.dot2tex_utils import have_dot2tex
        if have_dot2tex():
            G.set_latex_options(format='dot2tex',
                                prog='dot',
                                layout='acyclic',
                                edge_options=lambda x: {'backward':True})
        return G

    def _latex_(self):
        r"""
        Returns a latex method for the poset.

        EXAMPLES::

            sage: P = Poset(([1,2], [[1,2]]), cover_relations = True)
            sage: print P._latex_() #optional - dot2tex graphviz
            \begin{tikzpicture}
            %
            \useasboundingbox (0,0) rectangle (5.0cm,5.0cm);
            %
            \definecolor{cv0}{rgb}{0.0,0.0,0.0}
            \definecolor{cfv0}{rgb}{1.0,1.0,1.0}
            \definecolor{clv0}{rgb}{0.0,0.0,0.0}
            \definecolor{cv1}{rgb}{0.0,0.0,0.0}
            \definecolor{cfv1}{rgb}{1.0,1.0,1.0}
            \definecolor{clv1}{rgb}{0.0,0.0,0.0}
            \definecolor{cv0v1}{rgb}{0.0,0.0,0.0}
            %
            \Vertex[style={minimum size=1.0cm,draw=cv0,fill=cfv0,text=clv0,shape=circle},LabelOut=false,L=\hbox{$1$},x=0.0cm,y=0.0cm]{v0}
            \Vertex[style={minimum size=1.0cm,draw=cv1,fill=cfv1,text=clv1,shape=circle},LabelOut=false,L=\hbox{$2$},x=5.0cm,y=5.0cm]{v1}
            %
            \Edge[lw=0.1cm,style={post, bend right,color=cv0v1,},](v0)(v1)
            %
            \end{tikzpicture}
        """
        return self.hasse_diagram()._latex_()

    def _repr_(self):
        r"""
        Returns a string representation of the poset.

        TESTS::

            sage: partitions_of_five = {(5,):[(4,1),(3,2)], \
                    (4,1):[(3,1,1),(2,2,1)], \
                    (3,2):[(3,1,1),(2,2,1)], \
                    (3,1,1):[(2,1,1,1)], \
                    (2,2,1):[(2,1,1,1)], \
                    (2,1,1,1):[(1,1,1,1,1)], \
                    (1,1,1,1,1):[]}
            sage: P5 = Poset(partitions_of_five)
            sage: P5._repr_()
            'Finite poset containing 7 elements'
        """
        return "Finite poset containing %s elements"%self._hasse_diagram.order()

    def __iter__(self):
        """
        Iterates through the elements of a linear extension of the poset.

        EXAMPLES::

            sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: sorted(D.__iter__())
            [0, 1, 2, 3, 4]
        """
        return iter(self._list)

    def linear_extension(self, linear_extension=None, check=True):
        """
        Returns a linear extension of this poset.

        INPUT:

        - ``linear_extension`` -- a list of the elements of ``self`` (default: ``None``)
        - ``check`` -- a boolean (default: True);
          whether to check that ``linear_extension`` is indeed a
          linear extension of ``self``.

        .. seealso:: :meth:`is_linear_extension`, :meth:`linear_extensions`

        EXAMPLES::

            sage: P = Poset((divisors(15), attrcall("divides")), facade = True)

        Without optional argument, the default linear extension of the
        poset is returned, as a plain list::

            sage: P.linear_extension()
            [1, 3, 5, 15]

        Otherwise, a full-featured linear extension is constructed
        as an element of ``P.linear_extensions()``::

            sage: l = P.linear_extension([1,5,3,15]); l
            [1, 5, 3, 15]
            sage: type(l)
            <class 'sage.combinat.posets.linear_extensions.LinearExtensionsOfPoset_with_category.element_class'>
            sage: l.parent()
            The set of all linear extensions of Finite poset containing 4 elements

        By default, the linear extension is checked for correctness::

            sage: l = P.linear_extension([1,3,15,5])
            Traceback (most recent call last):
            ...
            ValueError: [1, 3, 15, 5] is not a linear extension of Finite poset containing 4 elements

        This can be disabled (at your own risks!) with::

            sage: P.linear_extension([1,3,15,5], check=False)
            [1, 3, 15, 5]

        .. todo::

            - Is it acceptable to have those two features for a single method?

            - In particular, we miss a short idiom to get the default
              linear extension
        """
        if linear_extension is not None:
            return self.linear_extensions()(linear_extension, check=check)
        else:
            # TODO: do we care whether this is a list or tuple?
            return list(self._list)

    @cached_method
    def linear_extensions(self, facade=False):
        """
        Returns the enumerated set of all the linear extensions of this poset

        INPUT:

        - ``facade`` -- a boolean (default: ``False``);
          whether to return the linear extensions as plain lists

          .. warning::

            The ``facade`` option is not yet fully functional::

                sage: P = Poset((divisors(12), attrcall("divides")), linear_extension=True)
                sage: L = P.linear_extensions(facade=True); L
                The set of all linear extensions of Finite poset containing 6 elements
                sage: L([1, 2, 3, 4, 6, 12])
                Traceback (most recent call last):
                ...
                TypeError: Cannot convert list to sage.structure.element.Element

        .. seealso:: :meth:`linear_extension`, :meth:`is_linear_extension`

        EXAMPLES::

            sage: P = Poset((divisors(12), attrcall("divides")), linear_extension=True)
            sage: P.list()
            [1, 2, 3, 4, 6, 12]
            sage: L = P.linear_extensions(); L
            The set of all linear extensions of Finite poset containing 6 elements
            sage: l = L.an_element(); l
            [1, 2, 3, 4, 6, 12]
            sage: L.cardinality()
            5
            sage: L.list()
            [[1, 2, 3, 4, 6, 12], [1, 2, 3, 6, 4, 12], [1, 2, 4, 3, 6, 12], [1, 3, 2, 4, 6, 12], [1, 3, 2, 6, 4, 12]]

        Each element is aware that it is a linear extension of `P`::

            sage: type(l.parent())
            <class 'sage.combinat.posets.linear_extensions.LinearExtensionsOfPoset_with_category'>

        With ``facade=True``, the elements of ``L`` are plain lists instead::

            sage: L = P.linear_extensions(facade=True)
            sage: l = L.an_element()
            sage: type(l)
            <type 'list'>

        .. warning::

            In Sage <= 4.8, this function used to return a plain list
            of lists. To recover the previous functionality, please
            use::

                sage: L = list(P.linear_extensions(facade=True)); L
                [[1, 2, 3, 4, 6, 12], [1, 2, 3, 6, 4, 12], [1, 2, 4, 3, 6, 12], [1, 3, 2, 4, 6, 12], [1, 3, 2, 6, 4, 12]]
                sage: type(L[0])
                <type 'list'>

        TESTS::

            sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: list(D.linear_extensions())
            [[0, 1, 2, 3, 4], [0, 1, 2, 4, 3], [0, 2, 1, 3, 4], [0, 2, 1, 4, 3], [0, 2, 4, 1, 3]]

        """
        from linear_extensions import LinearExtensionsOfPoset
        return LinearExtensionsOfPoset(self, facade = facade)

    def is_linear_extension(self, l):
        """
        Returns whether ``l`` is a linear extension of ``self``

        INPUT:

        - ``l`` -- a list (or iterable) containing all of the elements of ``self`` exactly once

        .. seealso:: :meth:`linear_extension`, :meth:`linear_extensions`

        EXAMPLES::

            sage: P = Poset((divisors(12), attrcall("divides")), facade=True, linear_extension=True)
            sage: P.list()
            [1, 2, 3, 4, 6, 12]
            sage: P.is_linear_extension([1, 2, 4, 3, 6, 12])
            True
            sage: P.is_linear_extension([1, 2, 4, 6, 3, 12])
            False

            sage: [p for p in Permutations(list(P)) if P.is_linear_extension(p)]
            [[1, 2, 3, 4, 6, 12], [1, 2, 3, 6, 4, 12], [1, 2, 4, 3, 6, 12], [1, 3, 2, 4, 6, 12], [1, 3, 2, 6, 4, 12]]
            sage: list(P.linear_extensions())
            [[1, 2, 3, 4, 6, 12], [1, 2, 3, 6, 4, 12], [1, 2, 4, 3, 6, 12], [1, 3, 2, 4, 6, 12], [1, 3, 2, 6, 4, 12]]

        .. note:: this is used and systematically tested in :class:`~sage.combinat.posets.linear_extensions.LinearExtensionsOfPosets`

        TESTS:

        Check that :trac:`15313` is fixed::

            sage: P = Poset((divisors(12), attrcall("divides")), facade=True, linear_extension=True)
            sage: P.is_linear_extension([1,2,4,3,6,12,1337])
            False
            sage: P.is_linear_extension([1,2,4,3,6,666,12,1337])
            False
            sage: P = Poset(DiGraph(5))
            sage: P.is_linear_extension(['David', 'McNeil', 'La', 'Lamentable', 'Aventure', 'de', 'Simon', 'Wiesenthal'])
            False
        """
        index = { x:i for (i,x) in enumerate(l) }
        return (len(l) == self.cardinality() and
                all(x in index for x in self) and
                all(index[i] < index[j] for (i,j) in self.cover_relations()))

    def list(self):
        """
        List the elements of the poset. This just returns the result
        of :meth:`linear_extension`.

        EXAMPLES::

            sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] }, facade = False)
            sage: D.list()
            [0, 1, 2, 3, 4]
            sage: type(D.list()[0])
            <class 'sage.combinat.posets.elements.FinitePoset_with_category.element_class'>
        """
        return list(self.linear_extension())

    def plot(self, label_elements=True, element_labels=None,
             vertex_size=300, vertex_colors=None,
             layout='acyclic',
             **kwds):
        """
        Returns a Graphic object for the Hasse diagram of the poset.

        The poset is increasing from bottom to top.

        By default, the vertices are labelled.

        If the poset is ranked, the plot uses the rank function for
        the heights of the vertices.

        INPUT:

        - ``label_elements`` (default: ``True``) - whether to display element labels

        - ``element_labels`` (default: ``None``) - a dictionary of element labels

        EXAMPLES::

            sage: D = Poset({ 1:[2,3], 2:[4], 3:[4,5] })
            sage: D.plot(label_elements=False)
            sage: D.plot()
            sage: type(D.plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: elm_labs = {1:'a', 2:'b', 3:'c', 4:'d', 5:'e'}
            sage: D.plot(element_labels=elm_labs)

        Plot of the empy poset::

            sage: P = Poset({})
            sage: P.plot()

        Plot of a ranked poset::

            sage: P = Poset(DiGraph('E@ACA@?'))
            sage: P.is_ranked()
            True
            sage: P.plot()

        TESTS:

        We check that ``label_elements`` and ``element_labels`` are honored::

            sage: def get_plot_labels(P): return sorted(t.string for t in P if isinstance(t, sage.plot.text.Text))
            sage: P1 = Poset({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: P2 = Poset({ 0:[1,2], 1:[3], 2:[3,4] }, facade=True)
            sage: get_plot_labels(P1.plot(label_elements=False))
            []
            sage: get_plot_labels(P1.plot(label_elements=True))
            ['0', '1', '2', '3', '4']
            sage: element_labels = {0:'a', 1:'b', 2:'c', 3:'d', 4:'e'}
            sage: get_plot_labels(P1.plot(element_labels=element_labels))
            ['a', 'b', 'c', 'd', 'e']
            sage: get_plot_labels(P2.plot(element_labels=element_labels))
            ['a', 'b', 'c', 'd', 'e']

        """
        from collections import defaultdict
        graph = self.hasse_diagram()
        rank_function = self.rank_function()
        if rank_function:
            heights = defaultdict(list)
        else:
            heights = None
        # if relabelling is needed
        if label_elements and element_labels is not None:
            relabelling = dict((self(element), label)
                               for (element, label) in element_labels.items())
            graph = graph.relabel(relabelling, inplace = False)
            if rank_function: # use the rank function to set the heights
                for i in self:
                    heights[rank_function(i)].append(relabelling[i])
        else: # otherwise
            if rank_function: # use the rank function to set the heights
                for i in self:
                    heights[rank_function(i)].append(i)
        return graph.plot(vertex_labels=label_elements,
                          vertex_size=vertex_size,
                          vertex_colors=vertex_colors,
                          layout=layout,
                          heights=heights,
                          **kwds)

    def show(self, label_elements=True, element_labels=None,
            vertex_size=300, vertex_colors=None, layout='acyclic', **kwds):
        """
        Shows the Graphics object corresponding the Hasse diagram of the
        poset. Optionally, it is labelled.

        INPUT:

        -  ``label_elements`` - whether to display element
           labels

        -  ``element_labels`` - a dictionary of element
           labels

        EXAMPLES::

            sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: D.plot(label_elements=False)
            sage: D.show()
            sage: elm_labs = {0:'a', 1:'b', 2:'c', 3:'d', 4:'e'}
            sage: D.show(element_labels=elm_labs)
        """
        self.plot(label_elements=label_elements, element_labels=element_labels,
            vertex_size=vertex_size, vertex_colors=vertex_colors, layout=layout).show(**kwds)

    @combinatorial_map(name="to graph")
    def to_graph(self):
        """
        Return the graph of ``self`` corresponding to forgetting the
        poset structure.

        EXAMPLES::

            sage: P = Poset({0:[1,2],1:[3],2:[3],3:[]})
            sage: P.to_graph()
            Graph on 4 vertices
            sage: P = Poset()
            sage: G = P.to_graph(); G
            Graph on 0 vertices

        Check that it is hashable::

            sage: hash(G) == hash(G)
            True
        """
        from sage.graphs.graph import Graph
        G = Graph(self.hasse_diagram())
        G._immutable = True
        return G

    def level_sets(self):
        """
        Return a list ``l`` such that ``l[i]`` is the set of minimal
        elements of the poset obtained from ``self`` by removing the
        elements in ``l[0], l[1], ..., l[i-1]``. (In particular,
        ``l[0]`` is the set of minimal elements of ``self``.)

        EXAMPLES::

            sage: P = Poset({0:[1,2],1:[3],2:[3],3:[]})
            sage: [len(x) for x in P.level_sets()]
            [1, 2, 1]

        ::

            sage: Q = Poset({0:[1,2], 1:[3], 2:[4], 3:[4]})
            sage: [len(x) for x in Q.level_sets()]
            [1, 2, 1, 1]
        """
        return [map(self._vertex_to_element, level) for level in
                self._hasse_diagram.level_sets()]

    def cover_relations(self,element=None):
        """
        Returns the list of pairs [u,v] of elements of the poset such that
        u v is a cover relation (that is, u v and there does not exist z
        such that u z v).

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: Q.cover_relations()
            [[1, 2], [0, 2], [2, 3], [3, 4]]
        """
        return [c for c in self.cover_relations_iterator()]

    def cover_relations_iterator(self):
        """
        Returns an iterator for the cover relations of the poset.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: type(Q.cover_relations_iterator())
            <type 'generator'>
            sage: [z for z in Q.cover_relations_iterator()]
            [[1, 2], [0, 2], [2, 3], [3, 4]]
        """
        for u,v,l in self._hasse_diagram.edge_iterator():
            yield map(self._vertex_to_element,(u,v))

    def relations(self):
        r"""
        Returns a list of all relations of the poset.

        OUTPUT:

        A list of pairs (each pair is a list), where the first element
        of the pair is less than or equal to the second element.

        Pairs are produced in a rough sort of lexicographic order,
        where earlier elements are from lower levels of the poset.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: Q.relations()
            [[1, 1], [1, 2], [1, 3], [1, 4], [0, 0], [0, 2], [0, 3], [0, 4], [2, 2], [2, 3], [2, 4], [3, 3], [3, 4], [4, 4]]

        AUTHOR:

        - Rob Beezer (2011-05-04)
        """
        return list(self.relations_iterator())

    def relations_iterator(self):
        r"""
        Returns an iterator for all the relations of the poset.

        OUTPUT:

        A generator that produces pairs (each pair is a list), where the
        first element of the pair is less than or equal to the second element.

        Pairs are produced in a rough sort of lexicographic order,
        where earlier elements are from lower levels of the poset.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: type(Q.relations_iterator())
            <type 'generator'>
            sage: [z for z in Q.relations_iterator()]
            [[1, 1], [1, 2], [1, 3], [1, 4], [0, 0], [0, 2], [0, 3], [0, 4], [2, 2], [2, 3], [2, 4], [3, 3], [3, 4], [4, 4]]

        AUTHOR:

        - Rob Beezer (2011-05-04)
        """
        # Relies on vertices the fact that _elements correspond to the rows and
        # columns of the lequal matrix
        leq_mat = self.lequal_matrix()
        n = leq_mat.nrows()
        elements = self._elements
        for i in range(n):
            for j in range(i, n):
                if leq_mat[i,j]:
                    yield [elements[i], elements[j]]

    def is_incomparable_chain_free(self, m, n = None):
        r"""
        Returns ``True`` if the poset is `(m+n)`-free (that is, there is no pair
        of incomparable chains of lengths `m` and `n`), and ``False`` if not.

        If ``m`` is a tuple of pairs of chain lengths, returns ``True`` if the poset
        does not contain a pair of incomparable chains whose lengths comprise
        one of the chain pairs, and ``False`` if not.
        A poset is `(m+n)`-free if it contains no induced subposet that is
        isomorphic to the poset consisting of two disjoint chains of lengths
        `m` and `n`.  See, for example, Exercise 15 in Chapter 3 of
        [EnumComb1]_.

        INPUT:

        - ``m`` - tuple of pairs of nonnegative integers
        - ``m``, ``n`` - nonnegative integers

        EXAMPLES::

            sage: P = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: P.is_incomparable_chain_free(1, 1)
            False
            sage: P.is_incomparable_chain_free(2, 1)
            True

        ::

            sage: P = Poset(((0, 1, 2, 3, 4), ((0, 1), (1, 2), (0, 3), (4, 2))))
            sage: P.is_incomparable_chain_free(((3, 1), (2, 2)))
            True

        ::

            sage: P = Poset((("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"), (("d", "a"), ("e", "a"), ("f", "a"), ("g", "a"), ("h", "b"), ("f", "b"), ("h", "c"), ("g", "c"), ("h", "d"), ("i", "d"), ("h", "e"), ("i", "e"), ("j", "f"), ("i", "f"), ("j", "g"), ("i", "g"), ("j", "h"))))
            sage: P.is_incomparable_chain_free(3, 1)
            True
            sage: P.is_incomparable_chain_free(2, 2)
            False

        ::

            sage: [len([p for p in Posets(n) if p.is_incomparable_chain_free(((3, 1), (2, 2)))]) for n in range(6)]
            [1, 1, 2, 5, 14, 42]

        TESTS::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: Q.is_incomparable_chain_free(2, 20/10)
            True
            sage: Q.is_incomparable_chain_free(2, pi)
            Traceback (most recent call last):
            ...
            TypeError: 2 and pi must be integers.
            sage: Q.is_incomparable_chain_free(2, -1)
            Traceback (most recent call last):
            ...
            ValueError: 2 and -1 must be nonnegative integers.
            sage: P = Poset(((0, 1, 2, 3, 4), ((0, 1), (1, 2), (0, 3), (4, 2))))
            sage: P.is_incomparable_chain_free((3, 1))
            Traceback (most recent call last):
            ...
            TypeError: (3, 1) is not a tuple of tuples.
            sage: P.is_incomparable_chain_free([3, 1], [2, 2])
            Traceback (most recent call last):
            ...
            TypeError: [3, 1] and [2, 2] must be integers.
            sage: P.is_incomparable_chain_free([[3, 1], [2, 2]])
            True
            sage: P.is_incomparable_chain_free(([3, 1], [2, 2]))
            True
            sage: P.is_incomparable_chain_free([3, 1], 2)
            Traceback (most recent call last):
            ...
            TypeError: [3, 1] and 2 must be integers.
            sage: P.is_incomparable_chain_free(([3, 1], [2, 2, 2]))
            Traceback (most recent call last):
            ...
            ValueError: '([3, 1], [2, 2, 2])' is not a tuple of length-2 tuples.

        AUTHOR:

        - Eric Rowland (2013-05-28)

        REFERENCES:

        .. [EnumComb1] Richard P. Stanley,
           *Enumerative Combinatorics, volume 1*,
           Second Edition,
           Cambridge University Press (2011).
           http://math.mit.edu/~rstan/ec/ec1/
        """
        if n is None:
            try:
                chain_pairs = [tuple(chain_pair) for chain_pair in m]
            except TypeError:
                raise TypeError('%s is not a tuple of tuples.' % str(tuple(m)))
            if not all(len(chain_pair) is 2 for chain_pair in chain_pairs):
                raise ValueError('%r is not a tuple of length-2 tuples.' % str(tuple(m)))
            return all(self.is_incomparable_chain_free(*chain_pair) for chain_pair in chain_pairs)
        try:
            m, n = Integer(m), Integer(n)
        except TypeError:
            raise TypeError('%s and %s must be integers.' % (m, n))
        if m < 0 or n < 0:
            raise ValueError("%s and %s must be nonnegative integers." % (m, n))
        twochains = digraphs.TransitiveTournament(m) + digraphs.TransitiveTournament(n)
        return self.hasse_diagram().transitive_closure().subgraph_search(twochains, induced = True) is None

    def is_lequal(self, x, y):
        """
        Returns ``True`` if `x` is less than or equal to `y` in the poset, and
        ``False`` otherwise.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = Q(0),Q(1),Q(4)
            sage: Q.is_lequal(x,y)
            False
            sage: Q.is_lequal(y,x)
            False
            sage: Q.is_lequal(x,z)
            True
            sage: Q.is_lequal(y,z)
            True
            sage: Q.is_lequal(z,z)
            True
        """
        r = self.compare_elements(x,y)
        return r == 0 or r == -1
    le = is_lequal

    def is_less_than(self, x, y):
        """
        Returns ``True`` if `x` is less than but not equal to `y` in the poset,
        and ``False`` otherwise.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = Q(0),Q(1),Q(4)
            sage: Q.is_less_than(x,y)
            False
            sage: Q.is_less_than(y,x)
            False
            sage: Q.is_less_than(x,z)
            True
            sage: Q.is_less_than(y,z)
            True
            sage: Q.is_less_than(z,z)
            False
        """
        return self.compare_elements(x,y) == -1
    lt = is_less_than

    def is_gequal(self, x, y):
        """
        Returns ``True`` if `x` is greater than or equal to `y` in the poset,
        and ``False`` otherwise.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = Q(0),Q(1),Q(4)
            sage: Q.is_gequal(x,y)
            False
            sage: Q.is_gequal(y,x)
            False
            sage: Q.is_gequal(x,z)
            False
            sage: Q.is_gequal(z,x)
            True
            sage: Q.is_gequal(z,y)
            True
            sage: Q.is_gequal(z,z)
            True
        """
        r = self.compare_elements(x,y)
        return r == 0 or r == 1
    ge = is_gequal

    def is_greater_than(self, x, y):
        """
        Returns ``True`` if `x` is greater than but not equal to `y` in the
        poset, and ``False`` otherwise.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = Q(0),Q(1),Q(4)
            sage: Q.is_greater_than(x,y)
            False
            sage: Q.is_greater_than(y,x)
            False
            sage: Q.is_greater_than(x,z)
            False
            sage: Q.is_greater_than(z,x)
            True
            sage: Q.is_greater_than(z,y)
            True
            sage: Q.is_greater_than(z,z)
            False
        """
        return self.compare_elements(x,y) == 1
    gt = is_greater_than

    def compare_elements(self, x, y):
        r"""
        Compare `x` and `y` in the poset.

        If ``x`` = ``y``, then ``0`` is returned;
        if ``x`` < ``y``, then ``-1`` is returned;
        if ``x`` > ``y``, then ``1`` is returned;
        and if ``x`` and ``y`` are not comparable,
        then ``None`` is returned.

        EXAMPLES::

            sage: P = Poset([[1,2],[4],[3],[4],[]])
            sage: P.compare_elements(0,0)
            0
            sage: P.compare_elements(0,4)
            -1
            sage: P.compare_elements(4,0)
            1
            sage: P.compare_elements(1,2)

        """
        i, j = map(self._element_to_vertex,(x,y))
        if i == j:
            return 0
        elif self._hasse_diagram.is_less_than(i, j):
            return -1
        elif self._hasse_diagram.is_less_than(j, i):
            return  1
        else:
            return None

    def minimal_elements(self):
        """
        Returns a list of the minimal elements of the poset.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P(0) in P.minimal_elements()
            True
            sage: P(1) in P.minimal_elements()
            True
            sage: P(2) in P.minimal_elements()
            True
        """
        return map(self._vertex_to_element, self._hasse_diagram.minimal_elements())

    def maximal_elements(self):
        """
        Returns a list of the maximal elements of the poset.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.maximal_elements()
            [4]
        """
        return map(self._vertex_to_element, self._hasse_diagram.maximal_elements())

    def bottom(self):
        """
        Returns the bottom element of the poset, if it exists.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.bottom() is None
            True
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.bottom()
            0
        """
        hasse_bot = self._hasse_diagram.bottom()
        if hasse_bot is None:
            return None
        else:
            return self._vertex_to_element(hasse_bot)

    def has_bottom(self):
        """
        Returns True if the poset has a unique minimal element.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.has_bottom()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.has_bottom()
            True
        """
        return self._hasse_diagram.has_bottom()

    def top(self):
        """
        Returns the top element of the poset, if it exists.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.top() is None
            True
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.top()
            1

        TESTS::

            sage: R = Poset([[0],[]])
            sage: R.list()
            [0]
            sage: R.top() #Trac #10776
            0

        """
        hasse_top = self._hasse_diagram.top()
        if hasse_top is None:
            return None
        else:
            return self._vertex_to_element(hasse_top)

    def has_top(self):
        """
        Returns True if the poset contains a unique maximal element, and
        False otherwise.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.has_top()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.has_top()
            True
        """
        return self._hasse_diagram.has_top()

    def is_bounded(self):
        """
        Return ``True`` if the poset ``self`` is bounded, and ``False``
        otherwise.

        We call a poset bounded if it contains a unique maximal element
        and a unique minimal element.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.is_bounded()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.is_bounded()
            True
        """
        return self._hasse_diagram.is_bounded()

    def is_chain(self):
        """
        Returns True if the poset is totally ordered, and False otherwise.

        EXAMPLES::

            sage: L = Poset({0:[1],1:[2],2:[3],3:[4]})
            sage: L.is_chain()
            True
            sage: V = Poset({0:[1,2]})
            sage: V.is_chain()
            False
        """
        return self._hasse_diagram.is_chain()

    def is_chain_of_poset(self, o, ordered=False):
        """
        Return whether an iterable ``o`` is a chain of ``self``,
        including a check for ``o`` being ordered from smallest
        to largest element if the keyword ``ordered`` is set to
        ``True``.

        INPUT:

        - ``o`` -- an iterable (e. g., list, set, or tuple)
          containing some elements of ``self``

        - ``ordered`` -- a Boolean (default: ``False``) which
          decides whether the notion of a chain includes being
          ordered

        OUTPUT:

        If ``ordered`` is set to ``False``, the truth value of
        the following assertion is returned: The subset of ``self``
        formed by the elements of ``o`` is a chain in ``self``.

        If ``ordered`` is set to ``True``, the truth value of
        the following assertion is returned: Every element of the
        list ``o`` is (strictly!) smaller than its successor in
        ``self``. (This makes no sense if ``ordered`` is a set.)

        EXAMPLES::

            sage: P = Poset((divisors(12), attrcall("divides")))
            sage: sorted(P.list())
            [1, 2, 3, 4, 6, 12]
            sage: P.is_chain_of_poset([2, 4])
            True
            sage: P.is_chain_of_poset([12, 6])
            True
            sage: P.is_chain_of_poset([12, 6], ordered=True)
            False
            sage: P.is_chain_of_poset([6, 12], ordered=True)
            True
            sage: P.is_chain_of_poset(())
            True
            sage: P.is_chain_of_poset((), ordered=True)
            True
            sage: P.is_chain_of_poset((3, 4, 12))
            False
            sage: P.is_chain_of_poset((3, 6, 12, 1))
            True
            sage: P.is_chain_of_poset((3, 6, 12, 1), ordered=True)
            False
            sage: P.is_chain_of_poset((3, 6, 12), ordered=True)
            True
            sage: P.is_chain_of_poset((1, 1, 3))
            True
            sage: P.is_chain_of_poset((1, 1, 3), ordered=True)
            False
            sage: P.is_chain_of_poset((1, 3), ordered=True)
            True
            sage: P.is_chain_of_poset((6, 1, 1, 3))
            True
            sage: P.is_chain_of_poset((2, 1, 1, 3))
            False
        """
        if ordered:
            sorted_o = o
            return all(self.lt(a, b) for a, b in zip(sorted_o, sorted_o[1:]))
        else:
            # _element_to_vertex can be assumed to be a linear extension
            # of the poset according to the documentation of class
            # HasseDiagram.
            sorted_o = sorted(o, key=self._element_to_vertex)
            return all(self.le(a, b) for a, b in zip(sorted_o, sorted_o[1:]))

    def is_EL_labelling(self, f, return_raising_chains=False):
        r"""
        Returns ``True`` if ``f`` is an EL labelling of ``self``.

        A labelling `f` of the edges of the Hasse diagram of a poset
        is called an EL labelling (edge lexicographic labelling) if
        for any two elements `u` and `v` with `u \leq v`,

            - there is a unique `f`-raising chain from `u` to `v` in
              the Hasse diagram, and this chain is lexicographically
              first among all chains from `u` to `v`.

        For more details, see [Bj1980]_.

        INPUT:

        - ``f`` -- a function taking two elements ``a`` and ``b`` in
          ``self`` such that ``b`` covers ``a`` and returning elements
          in a totally ordered set.

        - ``return_raising_chains`` (optional; default:``False``) if
          ``True``, returns the set of all raising chains in ``self``,
          if possible.

        EXAMPLES:

        Let us consider a Boolean poset::

            sage: P = Poset([[(0,0),(0,1),(1,0),(1,1)],[[(0,0),(0,1)],[(0,0),(1,0)],[(0,1),(1,1)],[(1,0),(1,1)]]],facade=True)
            sage: label = lambda a,b: min( i for i in [0,1] if a[i] != b[i] )
            sage: P.is_EL_labelling(label)
            True
            sage: P.is_EL_labelling(label,return_raising_chains=True)
            {((0, 0), (0, 1)): [1], ((0, 0), (1, 0)): [0], ((0, 1), (1, 1)): [0], ((1, 0), (1, 1)): [1], ((0, 0), (1, 1)): [0, 1]}

        REFERENCES:

            .. [Bj1980] Anders Björner,
               *Shellable and Cohen-Macaulay partially ordered sets*,
               Trans. Amer. Math. Soc. 260 (1980), 159-183,
               :doi:`10.1090/S0002-9947-1980-0570784-2`
        """
        label_dict = { (a,b):f(a,b) for a,b in self.cover_relations_iterator() }
        if return_raising_chains:
            raising_chains = {}
        for a,b in self.interval_iterator():
            P = self.subposet(self.interval(a,b))
            max_chains = sorted( [ [ label_dict[(chain[i],chain[i+1])] for i in range(len(chain)-1) ] for chain in P.maximal_chains() ] )
            if max_chains[0] != sorted(max_chains[0]) or any( max_chains[i] == sorted(max_chains[i]) for i in range(1,len(max_chains)) ):
                return False
            elif return_raising_chains:
                raising_chains[(a,b)] = max_chains[0]
        if return_raising_chains:
            return raising_chains
        else:
            return True

    def rank_function(self):
        r"""
        Return the (normalized) rank function of the poset,
        if it exists.

        A *rank function* of a poset `P` is a function `r`
        that maps elements of `P` to integers and satisfies:
        `r(x) = r(y) + 1` if `x` covers `y`. The function `r`
        is normalized such that its minimum value on every
        connected component of the Hasse diagram of `P` is
        `0`. This determines the function `r` uniquely (when
        it exists).

        OUTPUT:

        - a lambda function, if the poset admits a rank function
        - ``None``, if the poset does not admit a rank function

        EXAMPLES::

            sage: P = Poset(([1,2,3,4],[[1,4],[2,3],[3,4]]), facade=True)
            sage: P.rank_function() is not None
            True
            sage: P = Poset(([1,2,3,4,5],[[1,2],[2,3],[3,4],[1,5],[5,4]]), facade=True)
            sage: P.rank_function() is not None
            False
            sage: P = Poset(([1,2,3,4,5,6,7,8],[[1,4],[2,3],[3,4],[5,7],[6,7]]), facade=True)
            sage: f = P.rank_function(); f is not None
            True
            sage: f(5)
            0
            sage: f(2)
            0

        TESTS::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: r = P.rank_function()
            sage: for u,v in P.cover_relations_iterator():
            ....:     if r(v) != r(u) + 1:
            ....:         print "Bug in rank_function!"

        ::

            sage: Q = Poset([[1,2],[4],[3],[4],[]])
            sage: Q.rank_function() is None
            True
        """
        hasse_rf = self._hasse_diagram.rank_function()
        if hasse_rf is None:
            return None
        else:
            return lambda z: hasse_rf(self._element_to_vertex(z))

    def rank(self, element=None):
        r"""
        Return the rank of an element ``element`` in the poset ``self``,
        or the rank of the poset if ``element`` is ``None``.

        (The rank of a poset is the length of the longest chain of
        elements of the poset.)

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]], facade = False)
            sage: P.rank(5)
            2
            sage: P.rank()
            3
            sage: Q = Poset([[1,2],[3],[],[]])

            sage: P = Posets.SymmetricGroupBruhatOrderPoset(4)
            sage: [(v,P.rank(v)) for v in P]
            [('1234', 0),
             ('1324', 1),
            ...
             ('4231', 5),
             ('4321', 6)]
        """
        if element is None:
            return len(self.level_sets())-1
        elif self.is_ranked():
            return self.rank_function()(element)
        else:
            raise ValueError, "Poset is not ranked."

    def is_ranked(self):
        r"""
        Returns whether this poset is ranked.

        A poset is *ranked* if it admits a rank function. For more information
        about the rank function, see :meth:`~sage.combinat.posets.hasse_diagram.HasseDiagram.rank_function`.

        .. SEEALSO:: :meth:`is_graded`.

        EXAMPLES::

            sage: P = Poset([[1],[2],[3],[4],[]])
            sage: P.is_ranked()
            True
            sage: Q = Poset([[1,5],[2,6],[3],[4],[],[6,3],[4]])
            sage: Q.is_ranked()
            False
            sage: P = Poset( ([1,2,3,4],[[1,2],[2,4],[3,4]] ))
            sage: P.is_ranked()
            True
        """
        return bool(self.rank_function())

    def is_graded(self):
        r"""
        Returns whether this poset is graded.

        A poset is *graded* if all its maximal chains have the same length.
        There are various competing definitions for graded posets (see
        :wikipedia:`Graded_poset`). This definition is from section 3.1 of
        Richard Stanley's *Enumerative Combinatorics, Vol. 1* [EnumComb1]_.

        Note that every graded poset is ranked, but the converse is not
        true.

        .. SEEALSO:: :meth:`is_ranked`

        .. TODO::

            The current algorithm could be improvable. See
            :trac:`13223`.

        EXAMPLES::

            sage: P = Poset([[1],[2],[3],[4],[]])
            sage: P.is_graded()
            True
            sage: Q = Poset([[1,5],[2,6],[3],[4],[],[6,3],[4]])
            sage: Q.is_graded()
            False
            sage: P = Poset( ([1,2,3,4],[[1,2],[2,4],[3,4]] ))
            sage: P.is_graded()
            False
            sage: P = Poset({1: [2, 3], 4: [5]})
            sage: P.is_graded()
            True
            sage: P = Poset({1: [2, 3], 3: [4]})
            sage: P.is_graded()
            False
            sage: P = Poset({1: [2, 3], 4: []})
            sage: P.is_graded()
            False
            sage: P = Posets.BooleanLattice(4)
            sage: P.is_graded()
            True
            sage: P = RootSystem(['D',4]).root_poset()
            sage: P.is_graded()
            True
            sage: P = Poset({})
            sage: P.is_graded()
            True

        TESTS:

        Here we test that the empty poset is graded::

            sage: Poset([[],[]]).is_graded()
            True
        """
        # The code below is not really optimized, but beats looking at
        # every maximal chain...
        if len(self) <= 2:
            return True
        # Let's work with the Hasse diagram in order to avoid some
        # indirection (the output doesn't depend on the vertex labels).
        hasse = self._hasse_diagram
        rf = hasse.rank_function()
        if rf is None:
            return False    # because every graded poset is ranked.
        if not all(rf(i) == 0 for i in hasse.minimal_elements()):
            return False
        maxes = hasse.maximal_elements()
        rank = rf(maxes[0])
        return all(rf(i) == rank for i in maxes)

    def covers(self,x,y):
        """
        Returns True if y covers x and False otherwise.

        EXAMPLES::

            sage: Q = Poset([[1,5],[2,6],[3],[4],[],[6,3],[4]])
            sage: Q.covers(Q(1),Q(6))
            True
            sage: Q.covers(Q(1),Q(4))
            False
        """
        return self._hasse_diagram.has_edge(*map(self._element_to_vertex,(x,y)))

    def upper_covers_iterator(self,y):
        """
        Returns an iterator for the upper covers of the element y. An upper
        cover of y is an element x such that y x is a cover relation.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: type(Q.upper_covers_iterator(0))
            <type 'generator'>
        """
        for x in self._hasse_diagram.neighbor_out_iterator(self._element_to_vertex(y)):
            yield self._vertex_to_element(x)

    def upper_covers(self,y):
        """
        Returns a list of upper covers of the element y. An upper cover of
        y is an element x such that y x is a cover relation.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: map(Q.upper_covers,Q.list())
            [[2], [2], [3], [4], []]
        """
        return [x for x in self.upper_covers_iterator(y)]

    def lower_covers_iterator(self,y):
        """
        Returns an iterator for the lower covers of the element y. An lower
        cover of y is an element x such that y x is a cover relation.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: type(Q.lower_covers_iterator(0))
            <type 'generator'>
        """
        for x in self._hasse_diagram.neighbor_in_iterator(self._element_to_vertex(y)):
            yield self._vertex_to_element(x)

    def lower_covers(self,y):
        """
        Returns a list of lower covers of the element y. An lower cover of
        y is an element x such that y x is a cover relation.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: map(Q.lower_covers,Q.list())
            [[], [], [1, 0], [2], [3]]
        """
        return [x for x in self.lower_covers_iterator(y)]

    def cardinality(self):
        """
        Returns the number of elements in the poset.

        EXAMPLES::

            sage: Poset([[1,2,3],[4],[4],[4],[]]).cardinality()
            5
        """
        return Integer(self._hasse_diagram.order())

    size = deprecated_function_alias(8735, cardinality)

    def mobius_function(self,x,y):
        r"""
        Returns the value of the Mobius function of the poset on the
        elements x and y.

        EXAMPLES::

            sage: P = Poset([[1,2,3],[4],[4],[4],[]])
            sage: P.mobius_function(P(0),P(4))
            2
            sage: sum([P.mobius_function(P(0),v) for v in P])
            0
            sage: sum([abs(P.mobius_function(P(0),v)) \
            ...        for v in P])
            6
            sage: for u,v in P.cover_relations_iterator():
            ...    if P.mobius_function(u,v) != -1:
            ...        print "Bug in mobius_function!"

        ::

            sage: Q = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: Q.mobius_function(Q(0),Q(7))
            0
            sage: Q.mobius_function(Q(0),Q(5))
            0
            sage: Q.mobius_function(Q(2),Q(7))
            2
            sage: Q.mobius_function(Q(3),Q(3))
            1
            sage: sum([Q.mobius_function(Q(0),v) for v in Q])
            0
        """
        i,j = map(self._element_to_vertex,(x,y))
        return self._hasse_diagram.mobius_function(i,j)

    def mobius_function_matrix(self, ring = ZZ, sparse = False):
        r"""
        Returns a matrix whose ``(i,j)`` entry is the value of the Mobius
        function evaluated at ``self.linear_extension()[i]`` and
        ``self.linear_extension()[j]``.

        INPUT:

        - ``ring`` -- the ring of coefficients (default: ``ZZ``)
        - ``sparse`` -- whether the returned matrix is sparse or not
          (default: ``True``)

        EXAMPLES::

            sage: P = Poset([[4,2,3],[],[1],[1],[1]])
            sage: x,y = (P.linear_extension()[0],P.linear_extension()[1])
            sage: P.mobius_function(x,y)
            -1
            sage: M = P.mobius_function_matrix(); M
            [ 1 -1 -1 -1  2]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]
            [ 0  0  0  0  1]
            sage: M[0,4]
            2
            sage: M[0,1]
            -1

        We now demonstrate the usage of the optional parameters::

            sage: P.mobius_function_matrix(ring=QQ, sparse=False).parent()
            Full MatrixSpace of 5 by 5 dense matrices over Rational Field
        """
        M = self._hasse_diagram.mobius_function_matrix()
        if ring is not ZZ:
            M = M.change_ring(ring)
        if not sparse:
            M = M.dense_matrix()
        return M


    def lequal_matrix(self, ring = ZZ, sparse = False):
        """
        Computes the matrix whose ``(i,j)`` entry is 1 if
        ``self.linear_extension()[i] < self.linear_extension()[j]`` and 0
        otherwise.

        INPUT:

        - ``ring`` -- the ring of coefficients (default: ``ZZ``)
        - ``sparse`` -- whether the returned matrix is sparse or not
          (default: ``True``)

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]], facade = False)
            sage: LEQM = P.lequal_matrix(); LEQM
            [1 1 1 1 1 1 1 1]
            [0 1 0 1 0 0 0 1]
            [0 0 1 1 1 0 1 1]
            [0 0 0 1 0 0 0 1]
            [0 0 0 0 1 0 0 1]
            [0 0 0 0 0 1 1 1]
            [0 0 0 0 0 0 1 1]
            [0 0 0 0 0 0 0 1]
            sage: LEQM[1,3]
            1
            sage: P.linear_extension()[1] < P.linear_extension()[3]
            True
            sage: LEQM[2,5]
            0
            sage: P.linear_extension()[2] < P.linear_extension()[5]
            False

        We now demonstrate the usage of the optional parameters::

            sage: P.lequal_matrix(ring=QQ, sparse=False).parent()
            Full MatrixSpace of 8 by 8 dense matrices over Rational Field
        """
        M = self._hasse_diagram.lequal_matrix()
        if ring is not ZZ:
            M = M.change_ring(ring)
        if not sparse:
            M = M.dense_matrix()
        return M

    def coxeter_transformation(self):
        r"""
        Returns the matrix of the Auslander-Reiten translation acting
        on the Grothendieck group of the derived category of modules
        on the poset ``self``, in the basis of simple modules. This matrix is
        usually called the Coxeter transformation.

        EXAMPLES::

            sage: Posets.PentagonPoset().coxeter_transformation()
            [ 0  0  0  0 -1]
            [ 0  0  0  1 -1]
            [ 0  1  0  0 -1]
            [-1  1  1  0 -1]
            [-1  1  0  1 -1]

        TESTS::

            sage: M = Posets.PentagonPoset().coxeter_transformation()
            sage: M**8 == 1
            True
        """
        return self._hasse_diagram.coxeter_transformation()

    def meet_matrix(self):
        """
        Returns a matrix whose ``(i,j)`` entry is ``k``, where
        ``self.linear_extension()[k]`` is the meet (greatest lower bound) of
        ``self.linear_extension()[i]`` and ``self.linear_extension()[j]``.

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]], facade = False)
            sage: M = P.meet_matrix(); M
            [0 0 0 0 0 0 0 0]
            [0 1 0 1 0 0 0 1]
            [0 0 2 2 2 0 2 2]
            [0 1 2 3 2 0 2 3]
            [0 0 2 2 4 0 2 4]
            [0 0 0 0 0 5 5 5]
            [0 0 2 2 2 5 6 6]
            [0 1 2 3 4 5 6 7]
            sage: M[P(4).vertex,P(3).vertex] == P(0).vertex
            True
            sage: M[P(5).vertex,P(2).vertex] == P(2).vertex
            True
            sage: M[P(5).vertex,P(2).vertex] == P(5).vertex
            False
        """
        return self._hasse_diagram.meet_matrix()

    def is_meet_semilattice(self):
        r"""
        Returns True if self has a meet operation, and False otherwise.

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]], facade = False)
            sage: P.is_meet_semilattice()
            True

            sage: P = Poset([[1,2],[3],[3],[]])
            sage: P.is_meet_semilattice()
            True

            sage: P = Poset({0:[2,3],1:[2,3]})
            sage: P.is_meet_semilattice()
            False
        """
        return self._hasse_diagram.is_meet_semilattice()

    def join_matrix(self):
        """
        Returns a matrix whose ``(i,j)`` entry is ``k``, where
        ``self.linear_extension()[k]`` is the join (least upper bound) of
        ``self.linear_extension()[i]`` and ``self.linear_extension()[j]``.

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]], facade = False)
            sage: J = P.join_matrix(); J
            [0 1 2 3 4 5 6 7]
            [1 1 3 3 7 7 7 7]
            [2 3 2 3 4 6 6 7]
            [3 3 3 3 7 7 7 7]
            [4 7 4 7 4 7 7 7]
            [5 7 6 7 7 5 6 7]
            [6 7 6 7 7 6 6 7]
            [7 7 7 7 7 7 7 7]
            sage: J[P(4).vertex,P(3).vertex] == P(7).vertex
            True
            sage: J[P(5).vertex,P(2).vertex] == P(5).vertex
            True
            sage: J[P(5).vertex,P(2).vertex] == P(2).vertex
            False
        """
        return self._hasse_diagram.join_matrix()

    def is_join_semilattice(self):
        """
        Returns True is the poset has a join operation, and False
        otherwise.

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: P.is_join_semilattice()
            True

            sage: P = Poset([[1,2],[3],[3],[]])
            sage: P.is_join_semilattice()
            True

            sage: P = Poset({0:[2,3],1:[2,3]})
            sage: P.is_join_semilattice()
            False
        """
        return self._hasse_diagram.is_join_semilattice()

    def is_isomorphic(self,other):
        """
        Returns True if both posets are isomorphic.

        EXAMPLES::

            sage: P = Poset(([1,2,3],[[1,3],[2,3]]))
            sage: Q = Poset(([4,5,6],[[4,6],[5,6]]))
            sage: P.is_isomorphic( Q )
            True
        """
        if hasattr(other,'hasse_diagram'):
            return self.hasse_diagram().is_isomorphic( other.hasse_diagram() )
        else:
            raise ValueError, 'The input is not a finite poset.'

    import __builtin__ # Caveat: list is overridden by the method list above!!!
    def antichains(self, element_constructor = __builtin__.list):
        """
        Returns the antichains of the poset.

        INPUT:

         - ``element_constructor`` -- a function taking an iterable as
           argument (default: list)

        OUTPUT: an enumerated set

        An *antichain* of a poset is a collection of elements of the
        poset that are pairwise incomparable.

        EXAMPLES::

            sage: A = Posets.PentagonPoset().antichains(); A
            Set of antichains of Finite lattice containing 5 elements
            sage: list(A)
            [[], [0], [1], [1, 2], [1, 3], [2], [3], [4]]
            sage: A.cardinality()
            8
            sage: A[3]
            [1, 2]
            sage: list(Posets.AntichainPoset(3).antichains())
            [[], [2], [2, 1], [2, 1, 0], [2, 0], [1], [1, 0], [0]]
            sage: list(Posets.ChainPoset(3).antichains())
            [[], [0], [1], [2]]

        To get the antichains of a given size one can currently use::

            sage: list(A.elements_of_depth_iterator(2))
            [[1, 2], [1, 3]]

        Eventually the following syntax will be accepted::

            sage: A.subset(size = 2) # todo: not implemented

        To get the antichains as, say, sets, one may use the
        ``element_constructor`` option::

            sage: list(Posets.ChainPoset(3).antichains(element_constructor = set))
            [set([]), set([0]), set([1]), set([2])]

        .. note:: Internally, this uses
            :class:`sage.combinat.subsets_pairwise.PairwiseCompatibleSubsets`
            and :class:`SearchForest`. At this point, iterating
            through this set is about twice slower than using
            :meth:`antichains_iterator` (tested on
            ``posets.AntichainPoset(15)``). The algorithm is the same
            (depth first search through the tree), but
            :meth:`antichains_iterator` manually inlines things which
            apparently avoids some infrastructure overhead.

            On the other hand, this returns a full featured enumerated
            set, with containment testing, etc.

        """
        vertex_to_element = self._vertex_to_element
        def f(antichain):
            return element_constructor(vertex_to_element(x) for x in antichain)
        result = self._hasse_diagram.antichains(element_class = f)
        result.rename("Set of antichains of %s"%self)
        return result

    def antichains_iterator(self):
        """
        Returns an iterator over the antichains of the poset.

        EXAMPLES::

            sage: Posets.PentagonPoset().antichains_iterator()
            <generator object antichains_iterator at ...>

        .. seealso:: :meth:`antichains`
        """
        vertex_to_element = self._vertex_to_element
        for antichain in self._hasse_diagram.antichains_iterator():
            yield map(vertex_to_element, antichain)

    def chains(self, element_constructor=__builtin__.list, exclude=None):
        """
        Return all the chains of ``self``.

        INPUT:

        - ``element_constructor`` -- a function taking an iterable as
           argument (default: ``list``)

        - ``exclude`` -- elements of the poset to be excluded
          (default: ``None``)

        OUTPUT:

        The enumerated set of all chains of ``self``, each of which
        is given as an ``element_constructor``.

        A *chain* of a poset is a set of elements of the poset
        that are pairwise comparable.

        EXAMPLES::

            sage: A = Posets.PentagonPoset().chains(); A
            Set of chains of Finite lattice containing 5 elements
            sage: list(A)
            [[], [0], [0, 1], [0, 1, 4], [0, 2], [0, 2, 3], [0, 2, 3, 4], [0, 2, 4], [0, 3], [0, 3, 4], [0, 4], [1], [1, 4], [2], [2, 3], [2, 3, 4], [2, 4], [3], [3, 4], [4]]

        To get the chains of a given size one can currently use::

            sage: list(A.elements_of_depth_iterator(2))
            [[0, 1], [0, 2], [0, 3], [0, 4], [1, 4], [2, 3], [2, 4], [3, 4]]


        For bounded posets, one can exclude the bounds as follows::

            sage: P = Posets.DiamondPoset(5)
            sage: list(P.chains(exclude=[0, 4]))
            [[], [1], [2], [3]]

        Another example of exclusion of vertices::

            sage: P = Poset({1: [2, 3], 2: [4], 3: [4, 5]})
            sage: list(P.chains(element_constructor=tuple, exclude=[3]))
            [(), (1,), (1, 2), (1, 2, 4), (1, 4), (1, 5), (2,), (2, 4), (4,), (5,)]

        Eventually the following syntax will be accepted::

            sage: A.subset(size = 2) # todo: not implemented

        .. seealso:: :meth:`maximal_chains`, :meth:`antichains`
        """
        vertex_to_element = self._vertex_to_element
        def f(chain):
            return element_constructor(vertex_to_element(x) for x in chain)
        if not(exclude is None):
            exclude = [self._element_to_vertex(x) for x in exclude]
        result = self._hasse_diagram.chains(element_class = f,
                                            exclude=exclude)
        result.rename("Set of chains of %s" % self)
        return result

    def product(self,other):
        """
        Returns the cartesian product of ``self`` and ``other``.

        EXAMPLES::

            sage: P = Posets.ChainPoset(3)
            sage: Q = Posets.ChainPoset(4)
            sage: PQ = P.product(Q) ; PQ
            Finite poset containing 12 elements
            sage: len(PQ.hasse_diagram().edges())
            17
            sage: Q.product(P).is_isomorphic(PQ)
            True

            sage: P = Posets.BooleanLattice(2)
            sage: Q = P.product(P)
            sage: Q.is_isomorphic(Posets.BooleanLattice(4))
            True
        """
        return Poset(self._hasse_diagram.cartesian_product(other._hasse_diagram),cover_relations=True)

    def interval_iterator(self):
        """
        Returns an iterator over all pairs `x<y` in ``self``.

        EXAMPLES::

            sage: list(Posets.PentagonPoset().interval_iterator())
            [[0, 1], [0, 2], [0, 3], [0, 4], [1, 4], [2, 3], [2, 4], [3, 4]]

        .. seealso:: :meth:`maximal_chains`, :meth:`chains`
        """
        return self.chains().elements_of_depth_iterator(2)

    def dual(self):
        """
        Returns the dual poset of the given poset.

        EXAMPLE::

            sage: P = Poset(([1,2,3],[[1,2],[1,3]]))
            sage: P.cover_relations()
            [[1, 2], [1, 3]]
            sage: Q = P.dual()
            sage: Q.cover_relations()
            [[3, 1], [2, 1]]

            sage: P = LatticePoset([[1,2],[3],[3]], facade = True)
            sage: P.cover_relations()
            [[0, 1], [0, 2], [1, 3], [2, 3]]
            sage: Q = P.dual()
            sage: Q.cover_relations()
            [[3, 2], [3, 1], [2, 0], [1, 0]]
            sage: Q.category()
            Category of facade finite lattice posets
            sage: Q.__class__
            <class 'sage.combinat.posets.lattices.FiniteLatticePoset_with_category'>

            sage: P = MeetSemilattice([[1,2],[3],[3]])
            sage: P.dual().__class__
            <class 'sage.combinat.posets.lattices.FiniteJoinSemilattice_with_category'>
            sage: P = JoinSemilattice([[1,2],[3],[3]])
            sage: P.dual().__class__
            <class 'sage.combinat.posets.lattices.FiniteMeetSemilattice_with_category'>
        """
        n = self.cardinality()
        dual_hasse_digraph = DiGraph(self._hasse_diagram).reverse()
        dual_hasse_digraph.relabel(lambda i: n-i-1)

        return self._dual_class(dual_hasse_digraph,
                                elements = reversed(self._elements),
                                category = self.category(),
                                facade = self._is_facade)

    def relabel(self, relabelling):
        r"""
        Returns a copy of this poset with its elements relabelled

        INPUT:

        - ``relabelling`` -- a function or dictionnary

          This function should map each (non-wrapped) element of
          ``self`` to some distinct object.

        EXAMPLES::

            sage: P = Poset((divisors(12), attrcall("divides")), linear_extension=True, facade = False)
            sage: P.list()
            [1, 2, 3, 4, 6, 12]
            sage: P.cover_relations()
            [[1, 2], [1, 3], [2, 4], [2, 6], [3, 6], [4, 12], [6, 12]]
            sage: Q = P.relabel(lambda x: 12/x)
            sage: Q.list()
            [12, 6, 4, 3, 2, 1]
            sage: Q.cover_relations()
            [[12, 6], [12, 4], [6, 3], [6, 2], [4, 2], [3, 1], [2, 1]]

        Here we relabel the elements of a poset by {0,1,2, ...}, using
        a dictionary::

            sage: P = Poset((divisors(12), attrcall("divides")), linear_extension=True, facade = False)
            sage: relabelling = {c.element:i for (i,c) in enumerate(P)}; relabelling
            {1: 0, 2: 1, 3: 2, 4: 3, 6: 4, 12: 5}
            sage: Q = P.relabel(relabelling)
            sage: Q.list()
            [0, 1, 2, 3, 4, 5]
            sage: Q.cover_relations()
            [[0, 1], [0, 2], [1, 3], [1, 4], [2, 4], [3, 5], [4, 5]]

        Mind the ``c.element``; this is because the relabelling is
        applied to the elements of the poset without the wrapping.
        Thanks to this convention, the same relabelling function can
        be used both for facade or non facade posets::

            sage: P = Poset((divisors(12), attrcall("divides")), facade = True, linear_extension=True)
            sage: P.list()
            [1, 2, 3, 4, 6, 12]
            sage: Q = P.relabel(lambda x: 12/x)
            sage: Q.list()
            [12, 6, 4, 3, 2, 1]
            sage: Q.cover_relations()
            [[12, 6], [12, 4], [6, 3], [6, 2], [4, 2], [3, 1], [2, 1]]

        .. note::

            As can be seen in the above examples, the default linear
            extension of ``Q`` is that of ``P`` after relabelling. In
            particular, ``P`` and ``Q`` share the same internal Hasse
            diagram.
        """
        from sage.misc.stopgap import stopgap
        stopgap("Relabelling posets is known to break equality between posets (P == Q)", 14019)

        assert not isinstance(relabelling, (tuple, list)), "relabelling by tuple or list not yet defined"
        if isinstance(relabelling, dict):
            relabelling = relabelling.__getitem__
        elements = tuple(relabelling(x) for x in self._elements)
        return FinitePoset(self._hasse_diagram,
                           elements = elements,
                           category=self.category(),
                           facade=self._is_facade)

    def canonical_label(self):
        """
        Return the unique poset on the labels `\{0, \ldots, n-1\}` (where `n`
        is the number of elements in ``self``) that is isomorphic to ``self``
        and invariant in the isomorphism class.

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.canonical_label()`

        EXAMPLES::

            sage: P = Poset((divisors(12), attrcall("divides")), linear_extension=True, facade = False)
            sage: P.list()
            [1, 2, 3, 4, 6, 12]
            sage: P.cover_relations()
            [[1, 2], [1, 3], [2, 4], [2, 6], [3, 6], [4, 12], [6, 12]]
            sage: Q = P.canonical_label()
            sage: Q.list()
            [0, 1, 2, 3, 4, 5]
            sage: Q.cover_relations()
            [[0, 2], [0, 3], [1, 5], [2, 4], [3, 1], [3, 4], [4, 5]]

        As a facade::

            sage: P = Poset((divisors(12), attrcall("divides")), facade = True, linear_extension=True)
            sage: P.list()
            [1, 2, 3, 4, 6, 12]
            sage: P.cover_relations()
            [[1, 2], [1, 3], [2, 4], [2, 6], [3, 6], [4, 12], [6, 12]]
            sage: Q = P.canonical_label()
            sage: Q.list()
            [0, 1, 2, 3, 4, 5]
            sage: Q.cover_relations()
            [[0, 2], [0, 3], [1, 5], [2, 4], [3, 1], [3, 4], [4, 5]]
        """
        return FinitePoset(DiGraph(self._hasse_diagram).canonical_label(),
                           elements=range(len(self._elements)),
                           category=self.category(),
                           facade=self._is_facade)

    def with_linear_extension(self, linear_extension):
        """
        Returns a copy of ``self`` with a different default linear extension

        EXAMPLES::

            sage: P = Poset((divisors(12), attrcall("divides")), linear_extension=True)
            sage: P.cover_relations()
            [[1, 2], [1, 3], [2, 4], [2, 6], [3, 6], [4, 12], [6, 12]]
            sage: list(P)
            [1, 2, 3, 4, 6, 12]
            sage: Q = P.with_linear_extension([1,3,6,2,4,12])
            sage: list(Q)
            [1, 3, 6, 2, 4, 12]
            sage: Q.cover_relations()
            [[1, 3], [1, 2], [3, 6], [6, 12], [2, 6], [2, 4], [4, 12]]

        TESTS:

        We check that we can pass in a list of elements of P instead::

            sage: Q = P.with_linear_extension(map(P, [1,3,6,2,4,12]))
            sage: list(Q)
            [1, 3, 6, 2, 4, 12]
            sage: Q.cover_relations()
            [[1, 3], [1, 2], [3, 6], [6, 12], [2, 6], [2, 4], [4, 12]]

        We check that this works for facade posets too::

            sage: P = Poset((divisors(12), attrcall("divides")), facade=True)
            sage: Q = P.with_linear_extension([1,3,6,2,4,12])
            sage: list(Q)
            [1, 3, 6, 2, 4, 12]
            sage: Q.cover_relations()
            [[1, 3], [1, 2], [3, 6], [6, 12], [2, 6], [2, 4], [4, 12]]
            sage: sorted(Q.cover_relations()) == sorted(P.cover_relations())
            True

        .. note::

            With the current implementation, this requires relabelling
            the internal Dynkin diagram which is `O(n+m)`, where `n`
            is the number of elements and `m` the number of cover
            relations.

        """
        new_vertices = [ self._element_to_vertex(element) for element in linear_extension]
        new_elements = [ self._elements[i] for i in new_vertices ]
        vertex_relabelling = dict(zip(new_vertices, range(len(new_vertices))))
        return FinitePoset(self._hasse_diagram.relabel(vertex_relabelling, inplace=False),
                           elements = new_elements,
                           category=self.category(),
                           facade=self._is_facade)

    def graphviz_string(self,graph_string="graph",edge_string="--"):
        r"""
        Returns a representation in the DOT language, ready to render in
        graphviz.

        REFERENCES:

        - http://www.graphviz.org/doc/info/lang.html

        EXAMPLES::

            sage: P = Poset({'a':['b'],'b':['d'],'c':['d'],'d':['f'],'e':['f'],'f':[]})
            sage: print P.graphviz_string()
            graph {
            "f";"d";"b";"a";"c";"e";
            "f"--"e";"d"--"c";"b"--"a";"d"--"b";"f"--"d";
            }
        """
        s = '%s {\n' % graph_string
        for v in reversed(self.list()):
            s+= '"%s";'%v
        s+= '\n'
        for u, v in self.cover_relations_iterator():
            s+= '"%s"%s"%s";' % (v, edge_string, u)
        s+= "\n}"
        return s

    def subposet(self, elements):
        """
        Returns the poset containing elements with partial order induced by
        that of self.

        EXAMPLES::

            sage: P = Poset({"a":["c","d"], "b":["d","e"], "c":["f"], "d":["f"], "e":["f"]}, facade = False)
            sage: Q = P.subposet(["a","b","f"]); Q
            Finite poset containing 3 elements
            sage: Q.cover_relations()
            [[b, f], [a, f]]

        A subposet of a facade poset is again a facade poset::

            sage: P = Poset({"a":["c","d"], "b":["d","e"], "c":["f"], "d":["f"], "e":["f"]}, facade=True)
            sage: Q = P.subposet(["a","b","f"]); Q
            Finite poset containing 3 elements
            sage: Q.cover_relations()
            [['b', 'f'], ['a', 'f']]

        One may specified wrapped elements or not::

            sage: P = Poset({"a":["c","d"], "b":["d","e"], "c":["f"], "d":["f"], "e":["f"]}, facade = False)
            sage: Q = P.subposet([P("a"),P("b"),P("f")]); Q
            Finite poset containing 3 elements
            sage: Q.cover_relations()
            [[b, f], [a, f]]

            sage: B = posets.BooleanLattice(2)
            sage: above = B.principal_order_filter(0)
            sage: Q = B.subposet(above)
            sage: above_new = Q.principal_order_filter(Q.list()[0])
            sage: Q.subposet(above_new)
            Finite poset containing 4 elements

        TESTS::

            sage: P.subposet(("a","b","f"))
            Finite poset containing 3 elements
            sage: P.subposet(["a","b","x"])
            Traceback (most recent call last):
            ...
            ValueError: <type 'str'> is not an element of this poset
            sage: P.subposet(3)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
        """
        # Type checking is performed by the following line:
        elements  = [self(e) for e in elements]
        relations = []
        for u in elements:
            for v in elements:
                if self.is_less_than(u,v):
                    relations.append([u,v])
        if not self._is_facade:
            elements = [e.element for e in elements]
            relations = [[u.element,v.element] for u,v in relations]
        return Poset((elements, relations), cover_relations=False, facade=self._is_facade)

    def random_subposet(self, p):
        """
        Returns a random subposet that contains each element with
        probability p.

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: Q = P.random_subposet(.25)
        """
        elements = []
        p = float(p)
        for v in self:
            if random.random() <= p:
                elements.append(v)
        return self.subposet(elements)

    def order_filter(self,elements):
        """
        Returns the order filter generated by the elements of an
        iterable ``elements``.

        `I` is an order filter if, for any `x` in `I` and `y` such that
        `y \ge x`, then `y` is in `I`.

        EXAMPLES::

            sage: B = Posets.BooleanLattice(4)
            sage: B.order_filter([3,8])
            [3, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        """
        vertices = sorted(map(self._element_to_vertex,elements))
        of = self._hasse_diagram.order_filter(vertices)
        return map(self._vertex_to_element,of)

    def order_ideal(self,elements):
        """
        Returns the order ideal generated by the elements of an
        iterable ``elements``.

        `I` is an order ideal if, for any `x` in `I` and `y` such that
        `y \le x`, then `y` is in `I`.

        EXAMPLES::

            sage: B = Posets.BooleanLattice(4)
            sage: B.order_ideal([7,10])
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]
            sage: B.order_ideal(iter(range(4, 9)))
            [0, 1, 2, 3, 4, 5, 6, 7, 8]
        """
        vertices = map(self._element_to_vertex,elements)
        oi = self._hasse_diagram.order_ideal(vertices)
        return map(self._vertex_to_element,oi)

    def interval(self, x, y):
        """
        Returns a list of the elements `z` such that `x \le z \le y`.
        The order is that induced by the ordering in
        ``self.linear_extension()``.

        INPUT:


        -  ``x`` - any element of the poset

        -  ``y`` - any element of the poset


        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: P = Poset(dag)
            sage: I = set(map(P,[2,5,6,4,7]))
            sage: I == set(P.interval(2,7))
            True

        ::

            sage: dg = DiGraph({"a":["b","c"], "b":["d"], "c":["d"]})
            sage: P = Poset(dg, facade = False)
            sage: P.interval("a","d")
            [a, c, b, d]
        """
        return map(self._vertex_to_element,self._hasse_diagram.interval(
                self._element_to_vertex(x),self._element_to_vertex(y)))

    def closed_interval(self, x, y):
        """
        Returns a list of the elements `z` such that `x \le z \le y`.
        The order is that induced by the ordering in
        ``self.linear_extension()``.

        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: P = Poset(dag)
            sage: I = set(map(P,[2,5,6,4,7]))
            sage: I == set(P.closed_interval(2,7))
            True
        """
        return self.interval(x,y)

    def open_interval(self, x, y):
        """
        Returns a list of the elements `z` such that `x < z < y`. The
        order is that induced by the ordering in
        ``self.linear_extension()``.

        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: P = Poset(dag)
            sage: I = set(map(P,[5,6,4]))
            sage: I == set(P.open_interval(2,7))
            True

        ::

            sage: dg = DiGraph({"a":["b","c"], "b":["d"], "c":["d"]})
            sage: P = Poset(dg, facade = False)
            sage: P.open_interval("a","d")
            [c, b]
        """
        return map(self._vertex_to_element,self._hasse_diagram.open_interval(
                self._element_to_vertex(x),self._element_to_vertex(y)))

    def comparability_graph(self):
        r"""
        Returns the comparability graph of ``self``.

        See :wikipedia:`Comparability_graph`

        .. SEEALSO:: :meth:`incomparability_graph`, :mod:`sage.graphs.comparability`

        EXAMPLES::

            sage: p = posets.ChainPoset(4)
            sage: p.comparability_graph().is_isomorphic(graphs.CompleteGraph(4))
            True

            sage: p = posets.DiamondPoset(5)
            sage: g = p.comparability_graph(); g
            Comparability graph on 5 vertices
            sage: g.size()
            7
        """
        G = self.hasse_diagram().transitive_closure().to_undirected()
        G.rename('Comparability graph on %s vertices' % self.cardinality())
        return G

    def incomparability_graph(self):
        r"""
        Returns the incomparability graph of ``self``.

        This is the complement of the comparability graph.

        .. SEEALSO:: :meth:`comparability_graph`, :mod:`sage.graphs.comparability`

        EXAMPLES::

            sage: p = posets.ChainPoset(4)
            sage: p.incomparability_graph().size()
            0

            sage: p = posets.DiamondPoset(5)
            sage: g = p.incomparability_graph(); g
            Incomparability graph on 5 vertices
            sage: g.size()
            3
        """
        G = self.comparability_graph().complement()
        G.rename('Incomparability graph on %s vertices' % self.cardinality())
        return G

    def maximal_chains(self, partial=None):
        """
        Returns all maximal chains of this poset.

        Each chain is listed in increasing order.

        INPUT:

        -  ``partial`` - list (optional).  If present, find all maximal
           chains starting with the elements in partial.

        Returns list of the maximal chains of this poset.

        This is used in constructing the order complex for the poset.

        EXAMPLES::

            sage: P = Posets.BooleanLattice(3)
            sage: P.maximal_chains()
            [[0, 1, 3, 7], [0, 1, 5, 7], [0, 2, 3, 7], [0, 2, 6, 7], [0, 4, 5, 7], [0, 4, 6, 7]]
            sage: P.maximal_chains(partial=[0,2])
            [[0, 2, 3, 7], [0, 2, 6, 7]]
            sage: Q = Posets.ChainPoset(6)
            sage: Q.maximal_chains()
            [[0, 1, 2, 3, 4, 5]]
        """
        if partial is None or len(partial) == 0:
            start = self.minimal_elements()
            partial = []
        else:
            start = self.upper_covers(partial[-1])
        if len(start) == 0:
            return [partial]
        if len(start) == 1:
            return self.maximal_chains(partial=partial + start)
        parts = [partial + [x] for x in start]
        answer = []
        for new in parts:
            answer += self.maximal_chains(partial=new)
        return answer

    def order_complex(self, on_ints=False):
        """
        Returns the order complex associated to this poset.

        The order complex is the simplicial complex with vertices equal
        to the elements of the poset, and faces given by the chains.

        INPUT:

        - ``on_ints`` -- a boolean (default: False)

        EXAMPLES::

            sage: P = Posets.BooleanLattice(3)
            sage: S = P.order_complex(); S
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6, 7) and 6 facets
            sage: S.f_vector()
            [1, 8, 19, 18, 6]
            sage: S.homology()      # S is contractible
            {0: 0, 1: 0, 2: 0, 3: 0}
            sage: Q = P.subposet([1,2,3,4,5,6])
            sage: Q.order_complex().homology()    # a circle
            {0: 0, 1: Z}

            sage: P = Poset((divisors(15), attrcall("divides")), facade = True)
            sage: P.order_complex()
            Simplicial complex with vertex set (1, 3, 5, 15) and facets {(1, 3, 15), (1, 5, 15)}

        If ``on_ints``, then the elements of the poset are labelled
        `0,1,\dots` in the chain complex::

            sage: P.order_complex(on_ints=True)
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 3)}
        """
        from sage.homology.simplicial_complex import SimplicialComplex
        L = self.list()
        if on_ints:
            iso = dict( [ (L[i],i) for i in range(len(L)) ] )

        facets = []
        for f in self.maximal_chains():
            # TODO: factor out the logic for on_ints / facade / ...
            # We will want to do similar things elsewhere
            if on_ints:
                facets.append([iso[a]    for a in f])
            elif self._is_facade:
                facets.append([a         for a in f])
            else:
                facets.append([a.element for a in f])

        return SimplicialComplex(facets)

    def order_polytope(self):
        r"""
        Return the order polytope of the poset ``self``.

        The order polytope of a finite poset `P` is defined as the subset
        of `\RR^P` consisting of all maps `x : P \to \RR` satisfying

        .. MATH::

            0 \leq x(p) \leq 1 \mbox{ for all } p \in P,

        and

        .. MATH::

            x(p) \leq x(q) \mbox{ for all } p, q \in P
            \mbox{ satisfying } p < q.

        This polytope was defined and studied in [St1986]_.

        EXAMPLES::

            sage: P = posets.AntichainPoset(3)
            sage: Q = P.order_polytope();Q
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 8 vertices
            sage: P = posets.PentagonPoset()
            sage: Q = P.order_polytope();Q
            A 5-dimensional polyhedron in QQ^5 defined as the convex hull of 8 vertices

            sage: P = Poset([[1,2,3],[[1,2],[1,3]]])
            sage: Q = P.order_polytope()
            sage: Q.contains((1,0,0))
            False
            sage: Q.contains((0,1,1))
            True

        REFERENCES:

        .. [St1986] Richard Stanley. *Two poset polytopes*,
           Discrete Comput. Geom. (1986), :doi:`10.1007/BF02187680`
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        ineqs = [[0] + [ZZ(j==v)-ZZ(j==u) for j in self]
                 for u,v,w in self.hasse_diagram().edges()]
        for i in self.maximal_elements():
            ineqs += [[1] + [-ZZ(j==i) for j in self]]
        for i in self.minimal_elements():
            ineqs += [[0] + [ZZ(j==i) for j in self]]
        return Polyhedron(ieqs=ineqs)

    def chain_polytope(self):
        r"""
        Return the chain polytope of the poset ``self``.

        The chain polytope of a finite poset `P` is defined as the subset
        of `\RR^P` consisting of all maps `x : P \to \RR` satisfying

        .. MATH::

            x(p) \geq 0 \mbox{ for all } p \in P,

        and

        .. MATH::

            x(p_1) + x(p_2) + \ldots + x(p_k) \leq 1
            \mbox{ for all chains } p_1 < p_2 < \ldots < p_k
            \mbox{ in } P.

        This polytope was defined and studied in [St1986]_.

        EXAMPLES::

            sage: P = posets.AntichainPoset(3)
            sage: Q = P.chain_polytope();Q
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 8 vertices
            sage: P = posets.PentagonPoset()
            sage: Q = P.chain_polytope();Q
            A 5-dimensional polyhedron in QQ^5 defined as the convex hull of 8 vertices
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        ineqs=[[1]+[-ZZ(j in chain) for j in self] for chain in self.maximal_chains()]
        for i in self:
            ineqs+=[[0]+[ZZ(j==i) for j in self]]
        return Polyhedron(ieqs=ineqs)

    def zeta_polynomial(self):
        r"""
        Return the zeta polynomial of the poset ``self``.

        The zeta polynomial of a poset is the unique polynomial `Z(q)`
        such that for every integer `m > 1`, `Z(m)` is the number of
        weakly increasing sequences `x_1 \leq x_2 \leq \dots \leq x_{m-1}`
        of elements of the poset.

        The polynomial `Z(q)` is integral-valued, but generally doesn't
        have integer coefficients. It can be computed as

        .. MATH::

            Z(q) = \sum_{k \geq 1} \dbinom{q-2}{k-1} c_k,

        where `c_k` is the number of all chains of length `k` in the
        poset.

        For more information, see section 3.12 of [EnumComb1]_.

        In particular, `Z(2)` is the number of vertices and `Z(3)` is
        the number of intervals.

        EXAMPLES::

            sage: Posets.ChainPoset(2).zeta_polynomial()
            q
            sage: Posets.ChainPoset(3).zeta_polynomial()
            1/2*q^2 + 1/2*q
            sage: P = posets.PentagonPoset()
            sage: P.zeta_polynomial()
            1/6*q^3 + q^2 - 1/6*q
            sage: P = Posets.DiamondPoset(5)
            sage: P.zeta_polynomial()
            3/2*q^2 - 1/2*q

        TESTS:

        Checking the simplest cases::

            sage: Poset({}).zeta_polynomial()
            0
            sage: Poset({1: []}).zeta_polynomial()
            1
            sage: Poset({1: [], 2: []}).zeta_polynomial()
            2
        """
        q = polygen(QQ, 'q')
        g = sum(q**len(ch) for ch in self._hasse_diagram.chains())
        n = g.degree()
        f = g[max(n, 1)]
        while n > 1:
            f = (q - n)*f
            n = n - 1
            f = g[n] + f/n
        return f

    def f_polynomial(self):
        r"""
        Return the `f`-polynomial of a bounded poset ``self``.

        This is the `f`-polynomial of the order complex of the poset
        minus its bounds.

        The coefficient of `q^i` is the number of chains of
        `i+1` elements containing both bounds of the poset.

        .. SEEALSO::

            :meth:`is_bounded`, :meth:`h_polynomial`, :meth:`order_complex`,
            :meth:`sage.homology.cell_complex.GenericCellComplex.f_vector`

        .. WARNING::

            This is slightly different from the ``fPolynomial``
            method in Macaulay2.

        EXAMPLES::

            sage: P = Posets.DiamondPoset(5)
            sage: P.f_polynomial()
            3*q^2 + q
            sage: P = Poset({1:[2,3],2:[4],3:[5],4:[6],5:[7],6:[7]})
            sage: P.f_polynomial()
            q^4 + 4*q^3 + 5*q^2 + q

            sage: P = Poset({2: []})
            sage: P.f_polynomial()
            1
        """
        q = polygen(ZZ, 'q')
        hasse = self._hasse_diagram
        if len(hasse) == 1:
            return q.parent().one()
        maxi = hasse.top()
        mini = hasse.bottom()
        if (mini is None) or (maxi is None):
            raise TypeError('the poset is not bounded')
        return sum(q**(len(ch)+1) for ch in hasse.chains(exclude=[mini,
                                                                  maxi]))

    def h_polynomial(self):
        r"""
        Return the `h`-polynomial of a bounded poset ``self``.

        This is the `h`-polynomial of the order complex of the poset
        minus its bounds.

        This is related to the `f`-polynomial by a simple change
        of variables:

        .. MATH::

            h(q) = (1-q)^{\deg f} f \left( \frac{q}{1-q} \right),

        where `f` and `h` denote the `f`-polynomial and the
        `h`-polynomial, respectively.

        See :wikipedia:`h-vector`.

        .. SEEALSO::

            :meth:`is_bounded`, :meth:`f_polynomial`, :meth:`order_complex`,
            :meth:`sage.homology.simplicial_complex.SimplicialComplex.h_vector`

        .. WARNING::

            This is slightly different from the ``hPolynomial``
            method in Macaulay2.

        EXAMPLES::

            sage: P = Posets.AntichainPoset(3).order_ideals_lattice()
            sage: P.h_polynomial()
            q^3 + 4*q^2 + q
            sage: P = Posets.DiamondPoset(5)
            sage: P.h_polynomial()
            2*q^2 + q
            sage: P = Poset({1: []})
            sage: P.h_polynomial()
            1
        """
        q = polygen(ZZ, 'q')
        hasse = self._hasse_diagram
        if len(hasse) == 1:
            return q.parent().one()
        maxi = hasse.top()
        mini = hasse.bottom()
        if (mini is None) or (maxi is None):
            raise TypeError('the poset is not bounded')
        f = sum(q**(len(ch)) for ch in hasse.chains(exclude=[mini, maxi]))
        d = f.degree()
        f = (1-q)**d * q * f(q=q/(1-q))
        return q.parent(f)

    def flag_f_polynomial(self):
        r"""
        Return the flag `f`-polynomial of a bounded and ranked poset
        ``self``.

        This is the sum, over all chains containing both bounds,
        of a monomial encoding the ranks of the elements of the chain.

        More precisely, if `P` is a bounded ranked poset, then the
        flag `f`-polynomial of `P` is defined as the polynomial

        .. MATH::

            \sum_{\substack{p_0 < p_1 < \ldots < p_k, \\
                            p_0 = \min P, \ p_k = \max P}}
            x_{\rho(p_1)} x_{\rho(p_2)} \cdots x_{\rho(p_k)}
            \in \ZZ[x_1, x_2, \cdots, x_n]

        where `\min P` and `\max P` are (respectively) the minimum and
        the maximum of `P`, where `\rho` is the rank function of `P`
        (normalized to satisfy `\rho(\min P) = 0`), and where
        `n` is the rank of `\max P`. (Note that the indeterminate
        `x_0` doesn't actually appear in the polynomial.)

        For technical reasons, the polynomial is returned in the
        slightly larger ring `\ZZ[x_0, x_1, x_2, \cdots, x_{n+1}]` by
        this method.

        See :wikipedia:`h-vector`.

        .. SEEALSO:: :meth:`is_bounded`, :meth:`flag_h_polynomial`

        EXAMPLES::

            sage: P = Posets.DiamondPoset(5)
            sage: P.flag_f_polynomial()
            3*x1*x2 + x2

            sage: P = Poset({1:[2,3],2:[4],3:[5],4:[6],5:[6]})
            sage: fl = P.flag_f_polynomial(); fl
            2*x1*x2*x3 + 2*x1*x3 + 2*x2*x3 + x3
            sage: q = polygen(ZZ,'q')
            sage: fl(q,q,q,q) == P.f_polynomial()
            True

            sage: P = Poset({1:[2,3,4],2:[5],3:[5],4:[5],5:[6]})
            sage: P.flag_f_polynomial()
            3*x1*x2*x3 + 3*x1*x3 + x2*x3 + x3

            sage: P = Poset({2: [3]})
            sage: P.flag_f_polynomial()
            x1

            sage: P = Poset({2: []})
            sage: P.flag_f_polynomial()
            1
        """
        hasse = self._hasse_diagram
        maxi = hasse.top()
        mini = hasse.bottom()
        if (mini is None) or (maxi is None):
            raise TypeError('the poset is not bounded')
        rk = hasse.rank_function()
        if rk is None:
            raise TypeError('the poset should be ranked')
        n = rk(maxi)
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        if n == 0:
            return PolynomialRing(ZZ, 'x', 1).one()
        anneau = PolynomialRing(ZZ, 'x', n+1)
        x = anneau.gens()
        return x[n] * sum(prod(x[rk(i)] for i in ch) for ch in hasse.chains(exclude=[mini, maxi]))

    def flag_h_polynomial(self):
        r"""
        Return the flag `h`-polynomial of a bounded and ranked poset
        ``self``.

        If `P` is a bounded ranked poset whose maximal element has
        rank `n` (where the minimal element is set to have rank `0`),
        then the flag `h`-polynomial of `P` is defined as the
        polynomial

        .. MATH::

            \prod_{k=1}^n (1-x_k) \cdot f \left(\frac{x_1}{1-x_1},
            \frac{x_2}{1-x_2}, \cdots, \frac{x_n}{1-x_n}\right)
            \in \ZZ[x_1, x_2, \cdots, x_n],

        where `f` is the flag `f`-polynomial of `P` (see
        :meth:`flag_f_polynomial`).

        For technical reasons, the polynomial is returned in the
        slightly larger ring `\QQ[x_0, x_1, x_2, \cdots, x_{n+1}]` by
        this method.

        See :wikipedia:`h-vector`.

        .. SEEALSO:: :meth:`is_bounded`, :meth:`flag_f_polynomial`

        EXAMPLES::

            sage: P = Posets.DiamondPoset(5)
            sage: P.flag_h_polynomial()
            2*x1*x2 + x2

            sage: P = Poset({1:[2,3],2:[4],3:[5],4:[6],5:[6]})
            sage: fl = P.flag_h_polynomial(); fl
            -x1*x2*x3 + x1*x3 + x2*x3 + x3
            sage: q = polygen(ZZ,'q')
            sage: fl(q,q,q,q) == P.h_polynomial()
            True

            sage: P = Poset({1:[2,3,4],2:[5],3:[5],4:[5],5:[6]})
            sage: P.flag_h_polynomial()
            2*x1*x3 + x3

            sage: P = posets.ChainPoset(4)
            sage: P.flag_h_polynomial()
            x3

            sage: P = Poset({2: [3]})
            sage: P.flag_h_polynomial()
            x1

            sage: P = Poset({2: []})
            sage: P.flag_h_polynomial()
            1
        """
        hasse = self._hasse_diagram
        maxi = hasse.top()
        mini = hasse.bottom()
        if (mini is None) or (maxi is None):
            raise TypeError('the poset is not bounded')
        rk = hasse.rank_function()
        if rk is None:
            raise TypeError('the poset should be ranked')
        n = rk(maxi)
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        if n == 0:
            return PolynomialRing(QQ, 'x', 1).one()
        anneau = PolynomialRing(QQ, 'x', n+1)
        x = anneau.gens()
        return prod(1-x[k] for k in range(1, n)) * x[n] \
               * sum(prod(x[rk(i)]/(1-x[rk(i)]) for i in ch)
                     for ch in hasse.chains(exclude=[mini, maxi]))

    def characteristic_polynomial(self):
        r"""
        Return the characteristic polynomial of a graded poset ``self``.

        If `P` is a graded poset with rank `n` and a unique minimal
        element `\hat{0}`, then the characteristic polynomial of
        `P` is defined to be

        .. MATH::

            \sum_{x \in P} \mu(\hat{0}, x) q^{n-\rho(x)} \in \ZZ[q],

        where `\rho` is the rank function, and `\mu` is the Moebius
        function of `P`.

        See section 3.10 of [EnumComb1]_.

        EXAMPLES::

            sage: P = Posets.DiamondPoset(5)
            sage: P.characteristic_polynomial()
            q^2 - 3*q + 2
            sage: P = Poset({1:[2,3],2:[4],3:[5],4:[6],5:[6],6:[7]})
            sage: P.characteristic_polynomial()
            q^4 - 2*q^3 + q
            sage: P = Poset({1: []})
            sage: P.characteristic_polynomial()
            1
        """
        hasse = self._hasse_diagram
        rk = hasse.rank_function()
        if rk is None:
            raise TypeError('the poset should be ranked')
        n = rk(hasse.maximal_elements()[0])
        x0 = hasse.minimal_elements()[0]
        q = polygen(ZZ, 'q')
        return sum(hasse.mobius_function(x0, x) * q**(n - rk(x)) for x in hasse)

    def chain_polynomial(self):
        """
        Return the chain polynomial of ``self``.

        The coefficient of `q^k` is the number of chains of length `k`
        in ``self``. The length of a chain is the number of elements.

        .. WARNING::

            This is not what has been called the chain polynomial
            in [St1986]_. The latter is identical with the order
            polynomial (:meth:`order_polynomial`).

        EXAMPLES::

            sage: P = Posets.ChainPoset(3)
            sage: t = P.chain_polynomial(); t
            q^3 + 3*q^2 + 3*q + 1
            sage: t(1) == len(list(P.chains()))
            True

            sage: P = Posets.BooleanLattice(3)
            sage: P.chain_polynomial()
            6*q^4 + 18*q^3 + 19*q^2 + 8*q + 1

            sage: P = Posets.AntichainPoset(5)
            sage: P.chain_polynomial()
            5*q + 1

            sage: P = Poset({})
            sage: P.chain_polynomial()
            1
            sage: parent(P.chain_polynomial())
            Univariate Polynomial Ring in q over Integer Ring

            sage: R = Poset({1: []})
            sage: R.chain_polynomial()
            q + 1
        """
        hasse = self._hasse_diagram
        q = polygen(ZZ, 'q')
        one = q.parent().one()
        hasse_size = hasse.cardinality()
        chain_polys = [0]*hasse_size
        # chain_polys[i] will be the generating function for the
        # chains with topmost vertex i (in the labelling of the
        # Hasse diagram).
        for i in range(hasse_size):
            chain_polys[i] = q + sum(q*chain_polys[j]
                                     for j in hasse.principal_order_ideal(i))
        return one + sum(chain_polys)

    def order_polynomial(self):
        """
        Return the order polynomial of ``self``.

        The order polynomial `\Omega_P(q)` of a poset `P` is defined
        as the unique polynomial `S` such that for each integer
        `m \geq 1`, `S(m)` is the number of order-preserving maps
        from `P` to `\{1,\ldots,m\}`.

        See sections 3.12 and 3.15 of [EnumComb1]_, and also
        [St1986]_.

        .. SEEALSO:: :meth:`order_polytope`

        EXAMPLES::

            sage: P = Posets.AntichainPoset(3)
            sage: P.order_polynomial()
            q^3

            sage: P = Posets.ChainPoset(3)
            sage: f = P.order_polynomial(); f
            1/6*q^3 + 1/2*q^2 + 1/3*q
            sage: [f(i) for i in range(4)]
            [0, 1, 4, 10]
        """
        return self.order_ideals_lattice().zeta_polynomial()

    def promotion(self, i=1):
        r"""
        Computes the (extended) promotion on the linear extension of the poset ``self``

        INPUT:

        - ``i`` -- an integer between `1` and `n` (default: `1`)

        OUTPUT:

        - an isomorphic poset, with the same default linear extension

        The extended promotion is defined on a poset ``self`` of size
        `n` by applying the promotion operator `\tau_i \tau_{i+1}
        \cdots \tau_{n-1}` to the default linear extension `\pi` of ``self``
        (see :meth:`~sage.combinat.posets.linear_extensions.LinearExtensionOfPoset.promotion`),
        and relabelling ``self`` accordingly. For more details see [St2009]_.

        When the vertices of the poset ``self`` are labelled by
        `\{1,2,\ldots,n\}`, the linear extension is the identity, and
        `i=1`, the above algorithm corresponds to the promotion
        operator on posets defined by Schützenberger as
        follows. Remove `1` from ``self`` and replace it by the
        minimum `j` of all labels covering `1` in the poset. Then,
        remove `j` and replace it by the minimum of all labels
        covering `j`, and so on.  This process ends when a label is a
        local maximum. Place the label `n+1` at this vertex.  Finally,
        decrease all labels by `1`.

        REFERENCES:

            .. [St2009] Richard Stanley,
               *Promotion and evacuation*,
               Electron. J. Combin. 16 (2009), no. 2, Special volume in honor of Anders Björner,
               Research Paper 9, 24 pp.

        EXAMPLES::

            sage: P = Poset(([1,2], [[1,2]]), facade = False)
            sage: P.promotion()
            Finite poset containing 2 elements
            sage: P == P.promotion()
            True

            sage: P = Poset(([1,2,3,4,5,6,7], [[1,2],[1,4],[2,3],[2,5],[3,6],[4,7],[5,6]]))
            sage: P.list()
            [1, 2, 3, 5, 6, 4, 7]
            sage: Q = P.promotion(4); Q
            Finite poset containing 7 elements
            sage: Q.cover_relations()
            [[1, 2], [1, 6], [2, 3], [2, 5], [3, 7], [5, 7], [6, 4]]

        Note that if one wants to obtain the promotion defined by
        Schützenberger's algorithm directly on the poset, one needs
        to make sure the linear extension is the identity::

            sage: P = P.with_linear_extension([1,2,3,4,5,6,7])
            sage: P.list()
            [1, 2, 3, 4, 5, 6, 7]
            sage: Q = P.promotion(4); Q
            Finite poset containing 7 elements
            sage: Q.cover_relations()
            [[1, 2], [1, 6], [2, 3], [2, 4], [3, 5], [4, 5], [6, 7]]
            sage: Q = P.promotion()
            sage: Q.cover_relations()
            [[1, 2], [1, 3], [2, 4], [2, 5], [3, 6], [4, 7], [5, 7]]

        Here is an example for a poset not labelled by `\{1,2,\ldots,n\}`::

            sage: P = Poset((divisors(30), attrcall("divides")), linear_extension = True)
            sage: P.list()
            [1, 2, 3, 5, 6, 10, 15, 30]
            sage: P.cover_relations()
            [[1, 2], [1, 3], [1, 5], [2, 6], [2, 10], [3, 6], [3, 15], [5, 10], [5, 15], [6, 30], [10, 30], [15, 30]]
            sage: Q = P.promotion(4); Q
            Finite poset containing 8 elements
            sage: Q.cover_relations()
            [[1, 2], [1, 3], [1, 6], [2, 5], [2, 15], [3, 5], [3, 10], [5, 30], [6, 10], [6, 15], [10, 30], [15, 30]]

        .. seealso::

            - :meth:`linear_extension`
            - :meth:`with_linear_extension` and the ``linear_extension`` option of :func:`Poset`
            - :meth:`~sage.combinat.posets.linear_extensions.LinearExtensionOfPoset.promotion`
            - :meth:`evacuation`

        AUTHOR:

        - Anne Schilling (2012-02-18)
        """
        return self.linear_extension(self.linear_extension()).promotion(i).to_poset()

    def evacuation(self):
        r"""
        Compute evacuation on the linear extension associated to the poset ``self``.

        OUTPUT:

        - an isomorphic poset, with the same default linear extension

        Evacuation is defined on a poset ``self`` of size `n` by
        applying the evacuation operator
        `(\tau_1 \cdots \tau_{n-1}) (\tau_1 \cdots \tau_{n-2}) \cdots (\tau_1)`,
        to the default linear extension `\pi` of ``self``
        (see :meth:`~sage.combinat.posets.linear_extensions.LinearExtensionOfPoset.evacuation`),
        and relabelling ``self`` accordingly. For more details see [Stan2009]_.

        .. SEEALSO::

            - :meth:`linear_extension`
            - :meth:`with_linear_extension` and the ``linear_extension`` option of :func:`Poset`
            - :meth:`~sage.combinat.posets.linear_extensions.LinearExtensionOfPoset.evacuation`
            - :meth:`promotion`

        REFERENCES:

        .. [Stan2009] Richard Stanley,
           *Promotion and evacuation*,
           Electron. J. Combin. 16 (2009), no. 2, Special volume in honor of
           Anders Björner,
           Research Paper 9, 24 pp.

        EXAMPLES::

            sage: P = Poset(([1,2], [[1,2]]), facade = False)
            sage: P.evacuation()
            Finite poset containing 2 elements
            sage: P.evacuation() == P
            True

            sage: P = Poset(([1,2,3,4,5,6,7], [[1,2],[1,4],[2,3],[2,5],[3,6],[4,7],[5,6]]), linear_extension = True, facade = False)
            sage: P.list()
            [1, 2, 3, 4, 5, 6, 7]
            sage: Q = P.evacuation(); Q
            Finite poset containing 7 elements
            sage: Q.cover_relations()
            [[1, 2], [1, 3], [2, 5], [3, 4], [3, 6], [4, 7], [6, 7]]

        Note that the results depend on the linear extension associated to the poset::

            sage: P = Poset(([1,2,3,4,5,6,7], [[1,2],[1,4],[2,3],[2,5],[3,6],[4,7],[5,6]]))
            sage: P.list()
            [1, 2, 3, 5, 6, 4, 7]
            sage: Q = P.evacuation(); Q
            Finite poset containing 7 elements
            sage: Q.cover_relations()
            [[1, 2], [1, 5], [2, 3], [5, 6], [5, 4], [6, 7], [4, 7]]

        Here is an example of a poset where the vertices are not labelled by `\{1,2,\ldots,n\}`::

            sage: P = Poset((divisors(15), attrcall("divides")), linear_extension = True)
            sage: P.list()
            [1, 3, 5, 15]
            sage: Q = P.evacuation(); Q
            Finite poset containing 4 elements
            sage: Q.cover_relations()
            [[1, 3], [1, 5], [3, 15], [5, 15]]

        AUTHOR:

        - Anne Schilling (2012-02-18)
        """
        return self.linear_extension(self.linear_extension()).evacuation().to_poset()

    def is_slender(self):
        r"""
        Return whether the poset ``self`` is slender or not.

        It is assumed for this method that ``self`` is a finite graded poset.
        A finite poset `P` is called slender if every rank 2 interval contains
        three or four elements. See [Stan2009]_.

        EXAMPLES::

            sage: P = Poset(([1,2,3,4],[[1,2],[1,3],[2,4],[3,4]]), facade = True)
            sage: P.is_slender()
            True
            sage: P = Poset(([1,2,3,4,5],[[1,2],[1,3],[1,4],[2,5],[3,5],[4,5]]), facade = True)
            sage: P.is_slender()
            False

            sage: W = WeylGroup(['A',2])
            sage: G = W.bruhat_poset()
            sage: G.is_slender()
            True
            sage: W = WeylGroup(['A',3])
            sage: G = W.bruhat_poset()
            sage: G.is_slender()
            True
        """
        for x in self:
            d = {}
            S = self.upper_covers(x)
            Y = [ c for y in S for c in self.upper_covers(y) ]
            for y in Y:
                d[y] = d.get(y,0) + 1
            if not all( d[y]<3 for y in d.keys() ):
                return False
        return True

    def frank_network(self):
        r"""
        Computes Frank's network of the poset ``self``. This is defined in
        Section 8 of [BF1999]_.

        OUTPUT:

        A pair `(G, e)`, where `G` is Frank's network of `P` encoded as a
        :class:`DiGraph`, and `e` is the cost function on its edges encoded
        as a dictionary (indexed by these edges, which in turn are encoded
        as tuples of 2 vertices).

        .. NOTE::

            Frank's network of `P` is a certain directed graph with `2|P| + 2`
            vertices, defined in Section 8 of [BF1999]_. Its set of vertices
            consists of two vertices `(0, p)` and `(1, p)` for each element
            `p` of `P`, as well as two vertices `(-1, 0)` and `(2, 0)`.
            (These notations are not the ones used in [BF1999]_; see the table
            below for their relation.) The edges are:

            - for each `p` in `P`, an edge from `(-1, 0)` to `(0, p)`;

            - for each `p` in `P`, an edge from `(1, p)` to `(2, 0)`;

            - for each `p` and `q` in `P` such that `x \geq y`, an edge from
              `(0, p)` to `(1, q)`.

            We make this digraph into a network in the sense of flow theory as
            follows: The vertex `(-1, 0)` is considered as the source of this
            network, and the vertex `(2, 0)` as the sink. The cost function is
            defined to be `1` on the edge from `(0, p)` to `(1, p)` for each
            `p \in P`, and to be `0` on every other edge. The capacity is `1`
            on each edge. Here is how to translate this notations into that
            used in [BF1999]_::

              our notations                    [BF1999]
                 (-1, 0)                          s
                 (0, p)                          x_p
                 (1, p)                          y_p
                 (2, 0)                           t
                  a[e]                           a(e)

        REFERENCES:

        .. [BF1999] Thomas Britz, Sergey Fomin,
           *Finite posets and Ferrers shapes*,
           Advances in Mathematics 158, pp. 86-127 (2001),
           :arxiv:`math/9912126` (the arXiv version has less errors).

        EXAMPLES::

            sage: ps = [[16,12,14,-13],[[12,14],[14,-13],[12,16],[16,-13]]]
            sage: G, e = Poset(ps).frank_network()
            sage: G.edges()
            [((-1, 0), (0, -13), None), ((-1, 0), (0, 12), None), ((-1, 0), (0, 14), None), ((-1, 0), (0, 16), None), ((0, -13), (1, -13), None), ((0, -13), (1, 12), None), ((0, -13), (1, 14), None), ((0, -13), (1, 16), None), ((0, 12), (1, 12), None), ((0, 14), (1, 12), None), ((0, 14), (1, 14), None), ((0, 16), (1, 12), None), ((0, 16), (1, 16), None), ((1, -13), (2, 0), None), ((1, 12), (2, 0), None), ((1, 14), (2, 0), None), ((1, 16), (2, 0), None)]
            sage: e
            {((-1, 0), (0, 14)): 0, ((0, -13), (1, 12)): 0, ((-1, 0), (0, -13)): 0, ((0, 16), (1, 12)): 0, ((1, 16), (2, 0)): 0, ((0, -13), (1, 16)): 0, ((1, -13), (2, 0)): 0, ((0, -13), (1, -13)): 1, ((0, -13), (1, 14)): 0, ((-1, 0), (0, 16)): 0, ((0, 12), (1, 12)): 1, ((-1, 0), (0, 12)): 0, ((1, 14), (2, 0)): 0, ((1, 12), (2, 0)): 0, ((0, 14), (1, 12)): 0, ((0, 16), (1, 16)): 1, ((0, 14), (1, 14)): 1}
            sage: qs = [[1,2,3,4,5,6,7,8,9],[[1,3],[3,4],[5,7],[1,9],[2,3]]]
            sage: Poset(qs).frank_network()
            (Digraph on 20 vertices, {((0, 3), (1, 1)): 0, ((1, 8), (2, 0)): 0, ((-1, 0), (0, 3)): 0, ((0, 6), (1, 6)): 1, ((1, 9), (2, 0)): 0, ((0, 9), (1, 9)): 1, ((1, 7), (2, 0)): 0, ((0, 3), (1, 2)): 0, ((0, 3), (1, 3)): 1, ((0, 4), (1, 4)): 1, ((1, 2), (2, 0)): 0, ((0, 4), (1, 3)): 0, ((-1, 0), (0, 5)): 0, ((-1, 0), (0, 8)): 0, ((1, 3), (2, 0)): 0, ((0, 1), (1, 1)): 1, ((1, 1), (2, 0)): 0, ((0, 8), (1, 8)): 1, ((0, 4), (1, 1)): 0, ((1, 4), (2, 0)): 0, ((0, 2), (1, 2)): 1, ((-1, 0), (0, 1)): 0, ((0, 7), (1, 7)): 1, ((-1, 0), (0, 2)): 0, ((0, 7), (1, 5)): 0, ((0, 9), (1, 1)): 0, ((0, 5), (1, 5)): 1, ((-1, 0), (0, 9)): 0, ((-1, 0), (0, 7)): 0, ((0, 4), (1, 2)): 0, ((-1, 0), (0, 6)): 0, ((-1, 0), (0, 4)): 0, ((1, 6), (2, 0)): 0, ((1, 5), (2, 0)): 0})

        AUTHOR:

        - Darij Grinberg (2013-05-09)
        """
        from sage.graphs.digraph import DiGraph
        P0 = [(0, i) for i in self]
        pdict = { (-1, 0): P0, (2, 0): [] }
        for i in self:
            pdict[(0, i)] = [(1, j) for j in self if self.ge(i, j)]
            pdict[(1, i)] = [(2, 0)]
        G = DiGraph(pdict)
        a = { (u, v): 0 for (u, v, l) in G.edge_iterator() }
        for i in self:
            a[((0, i), (1, i))] = 1
        return (G, a)

    def greene_shape(self):
        r"""
        Return the Greene-Kleitman partition of ``self``.

        The Greene-Kleitman partition of a finite poset `P` is the partition
        `(c_1 - c_0, c_2 - c_1, c_3 - c_2, \ldots)`, where `c_k` is the
        maximum cardinality of a union of `k` chains of `P`. Equivalently,
        this is the conjugate of the partition `(a_1 - a_0, a_2 - a_1, a_3 -
        a_2, \ldots)`, where `a_k` is the maximum cardinality of a union of
        `k` antichains of `P`.

        See many sources, e. g., [BF1999]_, for proofs of this equivalence.

        EXAMPLES::

            sage: P = Poset([[3,2,1],[[3,1],[2,1]]])
            sage: P.greene_shape()
            [2, 1]
            sage: P = Poset([[1,2,3,4],[[1,4],[2,4],[4,3]]])
            sage: P.greene_shape()
            [3, 1]
            sage: P = Poset([[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22],[[1,4],[2,4],[4,3]]])
            sage: P.greene_shape()
            [3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: P = Poset([[],[]])
            sage: P.greene_shape()
            []

        AUTHOR:

        - Darij Grinberg (2013-05-09)
        """
        from sage.combinat.partition import Partition
        (G, a) = self.frank_network()
        n = len(self)
        chron = _ford_fulkerson_chronicle(G, (-1, 0), (2, 0), a)
        size = 0
        ps = []
        part = 0
        (pold, vold) = (0, 0)
        while size != n:
            (p, v) = chron.next()
            if v > vold:
                size += p
                if part > 0:
                    ps.append(part)
            elif p > pold:
                part += 1
            (pold, vold) = (p, v)
        ps.reverse()
        return Partition(ps)

    def p_partition_enumerator(self, tup, R, check=False):
        r"""
        Return a `P`-partition enumerator of ``self``.

        Given a total order `\prec` on the vertices of a poset `P`, a
        `P`-partition enumerator is the quasisymmetric function
        `\sum_f \prod_{p \in P} x_{f(p)}`, where the first sum is taken over
        all `P`-partitions `f`.

        A `P`-partition is a function `f : P \to \{1,2,3,...\}` satisfying
        the following properties for any two elements `i` and `j` of `P`:

        - if `i \prec j` then `f(i) \leq f(j)`,

        - if `j \prec i` then `f(j) < f(i)`.

        INPUT:

        - ``tup`` -- A tuple of elements of `P` representing a total order
          (this does not have to be a linear extension)

        - ``R`` -- A commutative ring

        OUTPUT:

        The `P`-partition enumerator of ``self`` according to ``tup`` in the
        algebra `QSym` over the base ring `R`.

        EXAMPLES::

            sage: P = Poset([[1,2,3,4],[[1,4],[2,4],[4,3]]])
            sage: FP = P.p_partition_enumerator((3,1,2,4), QQ, check=True); FP
            2*M[1, 1, 1, 1] + 2*M[1, 2, 1] + M[2, 1, 1] + M[3, 1]

            sage: expansion = FP.expand(5)
            sage: xs = expansion.parent().gens()
            sage: expansion == sum([xs[a]*xs[b]*xs[c]*xs[d] for a in range(5) for b in range(5) for c in range(5) for d in range(5) if a <= b and c <= b and b < d])
            True

            sage: P = Poset([[],[]])
            sage: FP = P.p_partition_enumerator((), QQ, check=True); FP
            M[]
        """
        if check:
            if sorted(self.list()) != sorted(tup):
                raise ValueError("the elements of tup are not those of P")
        from sage.combinat.composition import Composition
        from sage.combinat.ncsf_qsym.qsym import QuasiSymmetricFunctions
        QR = QuasiSymmetricFunctions(R)
        n = len(tup)
        res = QR.zero()
        tupdict = dict(zip(tup, range(n)))
        for lin in self.linear_extensions(facade=True):
            descents = [i + 1 for i in xrange(n-1) if tupdict[lin[i]] > tupdict[lin[i+1]]]
            res += QR.Fundamental()(Composition(from_subset=(descents, n)))
        return res

FinitePoset._dual_class = FinitePoset

##### Posets #####

class FinitePosets_n(UniqueRepresentation, Parent):
    r"""
    The finite enumerated set of all posets on `n` vertices, up to an isomorphism.

    EXAMPLES::

        sage: P = Posets(3)
        sage: P.cardinality()
        5
        sage: for p in P: print p.cover_relations()
        []
        [[1, 2]]
        [[0, 1], [0, 2]]
        [[0, 1], [1, 2]]
        [[0, 2], [1, 2]]
    """

    def __init__(self, n):
        r"""
        EXAMPLES::

            sage: P = Posets(3); P
            Posets containing 3 vertices
            sage: P.category()
            Category of finite enumerated sets
            sage: P.__class__
            <class 'sage.combinat.posets.posets.FinitePosets_n_with_category'>
            sage: TestSuite(P).run()
        """
        Parent.__init__(self, category = FiniteEnumeratedSets())
        self._n = n

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: P = Posets(3)
            sage: P._repr_()
            'Posets containing 3 vertices'
        """
        return "Posets containing %s vertices" % self._n

    def __contains__(self, P):
        """
        EXAMPLES::

            sage: posets.PentagonPoset() in Posets(5)
            True
            sage: posets.PentagonPoset() in Posets(3)
            False
            sage: 1 in Posets(3)
            False
        """
        return P in FinitePosets() and P.cardinality() == self._n

    def __iter__(self):
        """
        Returns an iterator of representatives of the isomorphism classes
        of finite posets of a given size.

        .. note::

           This uses the DiGraph iterator as a backend to construct
           transitively-reduced, acyclic digraphs.

        EXAMPLES::
            sage: P = Posets(2)
            sage: list(P)
            [Finite poset containing 2 elements, Finite poset containing 2 elements]
        """
        from sage.graphs.digraph_generators import DiGraphGenerators
        for dig in DiGraphGenerators()(self._n, is_poset):
            # We need to relabel the digraph since range(self._n) must be a linear
            # extension. Too bad we need to compute this again. TODO: Fix this.
            label_dict = dict(zip(dig.topological_sort(),range(dig.order())))
            yield FinitePoset(dig.relabel(label_dict,inplace=False))

    def cardinality(self, from_iterator=False):
        r"""
        Return the cardinality of this object.

        .. note::

            By default, this returns pre-computed values obtained from
            the On-Line Encyclopedia of Integer Sequences (:oeis:`A000112`).
            To override this, pass the argument ``from_iterator=True``.

        EXAMPLES::

            sage: P = Posets(3)
            sage: P.cardinality()
            5
            sage: P.cardinality(from_iterator=True)
            5
        """
        # Obtained from The On-Line Encyclopedia of Integer Sequences;
        # this is sequence number A000112.
        known_values = [1, 1, 2, 5, 16, 63, 318, 2045, 16999, 183231,
                2567284, 46749427, 1104891746, 33823827452, 1338193159771,
                68275077901156, 4483130665195087]
        if not from_iterator and self._n < len(known_values):
            return Integer(known_values[self._n])
        else:
            return super(FinitePosets_n, self).cardinality()

# For backward compatibility of pickles of the former Posets()
Posets_all = Posets

##### Miscellaneous functions #####

def is_poset(dig):
    r"""
    Tests whether a directed graph is acyclic and transitively
    reduced.

    EXAMPLES::

        sage: from sage.combinat.posets.posets import is_poset
        sage: dig = DiGraph({0:[2,3], 1:[3,4,5], 2:[5], 3:[5], 4:[5]})
        sage: is_poset(dig)
        False
        sage: is_poset(dig.transitive_reduction())
        True
    """
    return dig.is_directed_acyclic() and dig.is_transitively_reduced()

def _ford_fulkerson_chronicle(G, s, t, a):
    r"""
    Iterate through the Ford-Fulkerson algorithm for an acyclic directed
    graph with all edge capacities equal to `1`. This is an auxiliary algorithm
    for use by the :meth:`FinitePoset.greene_shape` method of finite posets,
    and is lacking some of the functionality that a general Ford-Fulkerson
    algorithm implementation should have.

    INPUT:

    - ``G`` -- an acyclic directed graph

    - ``s`` -- a vertex of `G` as the source

    - ``t`` -- a vertex of `G` as the sink

    - ``a`` -- a cost function (on the set of edges of ``G``) encoded as
      a dictionary. The keys of this dictionary are encoded as pairs
      of vertices.

    OUTPUT:

    An iterator which iterates through the values `(p, v)` during the
    application of the Ford-Fulkerson algorithm applied to the graph
    `G` with source `s`, sink `t`, cost function `a` and capacity `1`
    on each edge. Here, `p` denotes the value of the potential, and `v`
    denotes the value of the flow at every moment during the execution
    of the algorithm. The algorithm starts at `(p, v) = (0, 0)`.
    Every time ``next()`` is called, the iterator performs one step of
    the algorithm (incrementing either `p` or `v`) and yields the
    resulting pair `(p, v)`. Note that `(0, 0)` is never yielded.
    The iterator goes on for eternity, since the stopping condition
    is not implemented. This is OK for use in the ``greene_partition``
    function, since that one knows when to stop.

    The notation used here is that of Section 7 of [BF1999]_.

    .. WARNING::

        This method is tailor-made for its use in the
        :meth:`FinitePoset.greene_shape()` method of a finite poset. It's not
        very useful in general. First of all, as said above, the iterator
        does not know when to halt. Second, `G` needs to be acyclic for it
        to correctly work. This must be amended if this method is ever to be
        used outside the Greene-Kleitman partition construction. For the
        Greene-Kleitman partition, this is a non-issue since Frank's network
        is always acyclic.

    EXAMPLES::

        sage: from sage.combinat.posets.posets import _ford_fulkerson_chronicle
        sage: G = DiGraph({1: [3,6,7], 2: [4], 3: [7], 4: [], 6: [7,8], 7: [9], 8: [9,12], 9: [], 10: [], 12: []})
        sage: s = 1
        sage: t = 9
        sage: (1, 6, None) in G.edges()
        True
        sage: (1, 6) in G.edges()
        False
        sage: a = {(1, 6): 4, (2, 4): 0, (1, 3): 4, (1, 7): 1, (3, 7): 6, (7, 9): 1, (6, 7): 3, (6, 8): 1, (8, 9): 0, (8, 12): 2}
        sage: ffc = _ford_fulkerson_chronicle(G, s, t, a)
        sage: ffc.next()
        (1, 0)
        sage: ffc.next()
        (2, 0)
        sage: ffc.next()
        (2, 1)
        sage: ffc.next()
        (3, 1)
        sage: ffc.next()
        (4, 1)
        sage: ffc.next()
        (5, 1)
        sage: ffc.next()
        (5, 2)
        sage: ffc.next()
        (6, 2)
        sage: ffc.next()
        (7, 2)
        sage: ffc.next()
        (8, 2)
        sage: ffc.next()
        (9, 2)
        sage: ffc.next()
        (10, 2)
        sage: ffc.next()
        (11, 2)
    """
    from sage.graphs.digraph import DiGraph

    # pi: potential function as a dictionary.
    pi = { v: 0 for v in G.vertex_iterator() }
    # p: value of the potential pi.
    p = 0

    # f: flow function as a dictionary.
    f = { (u, v): 0 for (u, v, l) in G.edge_iterator() }
    # val: value of the flow f. (Can't call it v due to Python's asinine
    # handling of for loops.)
    val = 0

    # capacity: capacity function as a dictionary. Here, just the
    # indicator function of the set of arcs of G.
    capacity = { (u, v): 1 for (u, v, l) in G.edge_iterator() }

    while True:

        # Step MC1 in Britz-Fomin, Algorithm 7.2.

        # Gprime: directed graph G' from Britz-Fomin, Section 7.
        Gprime = DiGraph()
        Gprime.add_vertices(G.vertices())
        for (u,v,l) in G.edge_iterator():
            if pi[v] - pi[u] == a[(u, v)]:
                if f[(u, v)] < capacity[(u, v)]:
                    Gprime.add_edge(u, v)
                elif f[(u, v)] > 0:
                    Gprime.add_edge(v, u)

        # X: list of vertices of G' reachable from s, along with
        # the shortest paths from s to them.
        X = Gprime.shortest_paths(s)
        if t in X.keys():
            # Step MC2a in Britz-Fomin, Algorithm 7.2.
            shortest_path = X[t]
            shortest_path_in_edges = zip(shortest_path[:-1],shortest_path[1:])
            for (u, v) in shortest_path_in_edges:
                if v in G.neighbors_out(u):
                    f[(u, v)] += 1
                else:
                    f[(v, u)] -= 1
            val += 1
        else:
            # Step MC2b in Britz-Fomin, Algorithm 7.2.
            for v in G.vertex_iterator():
                if not X.has_key(v):
                    pi[v] += 1
            p += 1

        yield (p, val)

