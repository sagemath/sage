r"""
Strong and weak tableaux

There are two types of `k`-tableaux: strong `k`-tableaux and weak `k`-tableaux.
Standard weak `k`-tableaux correspond to saturated chains in the weak order,
whereas standard strong `k`-tableaux correspond to saturated chains in the strong Bruhat order.
For semistandard tableaux, the notion of weak and strong horizontal strip is necessary.
More information can be found in [LLMS2006]_ .

.. SEEALSO:: :meth:`sage.combinat.k_tableau.StrongTableau`, :meth:`sage.combinat.k_tableau.WeakTableau`

Authors:

- Anne Schilling and Mike Zabrocki (2013): initial version
- Avi Dalal and Nate Gallup (2013): implementation of `k`-charge
"""
#*****************************************************************************
#       Copyright (C) 2013 Anne Schilling <anne at math.ucdavis.edu>
#                          Mike Zabrocki  <zabrocki at mathstat.yorku.ca>
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
#****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.parent import Parent
from sage.structure.list_clone import ClonableList
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.combinat.skew_tableau import SkewTableau, SemistandardSkewTableaux
from sage.combinat.partition import Partition, Partitions
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.combinat.core import Core
from sage.rings.integer_ring import ZZ
from sage.functions.generalized import sgn
from sage.misc.flatten import flatten
from sage.combinat.skew_partition import SkewPartition
from sage.combinat.tableau import Tableaux
from sage.combinat.composition import Composition
import copy


def WeakTableau(t, k, inner_shape = [], representation = "core"):
    r"""
    This is the dispatcher method for the element class of weak `k`-tableaux.

    Standard weak `k`-tableaux correspond to saturated chains in the weak order.
    There are three formulations of weak tableaux, one in terms of cores, one in terms
    of `k`-bounded partitions, and one in terms of factorizations of affine Grassmannian
    elements. For semistandard weak `k`-tableaux, all letters of the same value have to
    satisfy the conditions of a horizontal strip. In the affine Grassmannian formulation this
    means that all factors are cyclically decreasing elements. For more information, see
    for example [LLMSSZ2013]_.

    INPUT:

    - ``t`` -- a weak `k`-tableau in the specified representation:

      - for the 'core' representation ``t`` is a list of lists where each subtableaux
        should have a `k+1`-core shape; ``None`` is allowed as an entry for skew weak
        `k`-tableaux
      - for the 'bounded' representation ``t`` is a list of lists where each subtableaux
        should have a `k`-bounded shape; ``None`` is allowed as an entry for skew weak
        `k`-tableaux
      - for the 'factorized_permutation' representation ``t`` is either a list of
        cyclically decreasing Weyl group elements or a list of reduced words of cyclically
        decreasing Weyl group elements; to indicate a skew tableau in this representation,
        ``inner_shape`` should be the inner shape as a `(k+1)`-core

    - ``k`` -- positive integer

    - ``inner_shape`` -- this entry is only relevant for the 'factorized_permutation'
      representation and specifies the inner shape in case the tableau is skew
      (default: ``[]``)

    - ``representation`` -- 'core', 'bounded', or 'factorized_permutation'
      (default: 'core')

    EXAMPLES:

    Here is an example of a weak 3-tableau in core representation::

        sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
        sage: t.shape()
        [5, 2, 1]
        sage: t.weight()
        (2, 2, 2)
        sage: type(t)
        <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>

    And now we give a skew weak 3-tableau in core representation::

        sage: ts = WeakTableau([[None, 1, 1, 2, 2], [None, 2], [1]], 3)
        sage: ts.shape()
        ([5, 2, 1], [1, 1])
        sage: ts.weight()
        (2, 2)
        sage: type(ts)
        <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>

    Next we create the analogue of the first example in bounded representation::

        sage: tb = WeakTableau([[1,1,2],[2,3],[3]], 3, representation="bounded")
        sage: tb.shape()
        [3, 2, 1]
        sage: tb.weight()
        (2, 2, 2)
        sage: type(tb)
        <class 'sage.combinat.k_tableau.WeakTableaux_bounded_with_category.element_class'>
        sage: tb.to_core_tableau()
        [[1, 1, 2, 2, 3], [2, 3], [3]]
        sage: t == tb.to_core_tableau()
        True

    And the analogue of the skew example in bounded representation::

        sage: tbs = WeakTableau([[None, 1, 2], [None, 2], [1]], 3, representation = "bounded")
        sage: tbs.shape()
        ([3, 2, 1], [1, 1])
        sage: tbs.weight()
        (2, 2)
        sage: tbs.to_core_tableau()
        [[None, 1, 1, 2, 2], [None, 2], [1]]
        sage: ts.to_bounded_tableau() == tbs
        True

    Finally we do the same examples for the factorized permutation representation::

        sage: tf = WeakTableau([[2,0],[3,2],[1,0]], 3, representation = "factorized_permutation")
        sage: tf.shape()
        [5, 2, 1]
        sage: tf.weight()
        (2, 2, 2)
        sage: type(tf)
        <class 'sage.combinat.k_tableau.WeakTableaux_factorized_permutation_with_category.element_class'>
        sage: tf.to_core_tableau() == t
        True

        sage: tfs = WeakTableau([[0,3],[2,1]], 3, inner_shape = [1,1], representation = 'factorized_permutation')
        sage: tfs.shape()
        ([5, 2, 1], [1, 1])
        sage: tfs.weight()
        (2, 2)
        sage: type(tfs)
        <class 'sage.combinat.k_tableau.WeakTableaux_factorized_permutation_with_category.element_class'>
        sage: tfs.to_core_tableau()
        [[None, 1, 1, 2, 2], [None, 2], [1]]

    Another way to pass from one representation to another is as follows::

        sage: ts
        [[None, 1, 1, 2, 2], [None, 2], [1]]
        sage: ts.parent()._representation
        'core'
        sage: ts.representation('bounded')
        [[None, 1, 2], [None, 2], [1]]

    To test whether a given semistandard tableau is a weak `k`-tableau in the bounded representation,
    one can ask::

        sage: t = Tableau([[1,1,2],[2,3],[3]])
        sage: t.is_k_tableau(3)
        True
        sage: t = SkewTableau([[None, 1, 2], [None, 2], [1]])
        sage: t.is_k_tableau(3)
        True
        sage: t = SkewTableau([[None, 1, 1], [None, 2], [2]])
        sage: t.is_k_tableau(3)
        False

    TESTS::

        sage: t = WeakTableau([[2,0],[3,2],[1,0]], 3, representation = "bla")
        Traceback (most recent call last):
        ...
        NotImplementedError: The representation option needs to be 'core', 'bounded', or 'factorized_permutation'
    """
    if representation == "core":
        return WeakTableau_core(t, k)
    elif representation == "bounded":
        return WeakTableau_bounded(t, k)
    elif representation == "factorized_permutation":
        return WeakTableau_factorized_permutation(t, k, inner_shape = inner_shape)
    else:
        raise NotImplementedError("The representation option needs to be 'core', 'bounded', or 'factorized_permutation'")

def WeakTableaux(k, shape , weight, representation = "core"):
    r"""
    This is the dispatcher method for the parent class of weak `k`-tableaux.

    INPUT:

    - ``k`` -- positive integer
    - ``shape`` -- shape of the weak `k`-tableaux; for the 'core' and
      'factorized_permutation' representation, the shape is inputted as a `(k+1)`-core;
      for the 'bounded' representation, the shape is inputted as a `k`-bounded partition;
      for skew tableaux, the shape is inputted as a tuple of the outer and inner shape
    - ``weight`` -- the weight of the weak `k`-tableaux as a list or tuple
    - ``representation`` -- 'core', 'bounded', or 'factorized_permutation' (default: 'core')

    EXAMPLES::

        sage: T = WeakTableaux(3, [5,2,1], [1,1,1,1,1,1])
        sage: T.list()
        [[[1, 3, 4, 5, 6], [2, 6], [4]],
        [[1, 2, 4, 5, 6], [3, 6], [4]],
        [[1, 2, 3, 4, 6], [4, 6], [5]],
        [[1, 2, 3, 4, 5], [4, 5], [6]]]
        sage: T.cardinality()
        4

        sage: T = WeakTableaux(3, [[5,2,1], [2]], [1,1,1,1])
        sage: T.list()
        [[[None, None, 2, 3, 4], [1, 4], [2]],
        [[None, None, 1, 2, 4], [2, 4], [3]],
        [[None, None, 1, 2, 3], [2, 3], [4]]]

        sage: T = WeakTableaux(3, [3,2,1], [1,1,1,1,1,1], representation = 'bounded')
        sage: T.list()
        [[[1, 3, 5], [2, 6], [4]],
        [[1, 2, 5], [3, 6], [4]],
        [[1, 2, 3], [4, 6], [5]],
        [[1, 2, 3], [4, 5], [6]]]

        sage: T = WeakTableaux(3, [[3,2,1], [2]], [1,1,1,1], representation = 'bounded')
        sage: T.list()
        [[[None, None, 3], [1, 4], [2]],
        [[None, None, 1], [2, 4], [3]],
        [[None, None, 1], [2, 3], [4]]]

        sage: T = WeakTableaux(3, [5,2,1], [1,1,1,1,1,1], representation = 'factorized_permutation')
        sage: T.list()
        [[s0, s3, s2, s1, s3, s0],
        [s0, s3, s2, s3, s1, s0],
        [s0, s2, s3, s2, s1, s0],
        [s2, s0, s3, s2, s1, s0]]

        sage: T = WeakTableaux(3, [[5,2,1], [2]], [1,1,1,1], representation = 'factorized_permutation')
        sage: T.list()
        [[s0, s3, s2, s3], [s0, s2, s3, s2], [s2, s0, s3, s2]]
    """
    if representation == "core":
        return WeakTableaux_core(k, shape, weight)
    elif representation == "bounded":
        return WeakTableaux_bounded(k, shape, weight)
    elif representation == "factorized_permutation":
        return WeakTableaux_factorized_permutation(k, shape, weight)
    else:
        raise NotImplementedError("The representation option needs to be 'core', 'bounded', or 'factorized_permutation'")

#Abstract class for the elements of weak tableau
class WeakTableau_abstract(ClonableList,
        metaclass=InheritComparisonClasscallMetaclass):
    r"""
    Abstract class for the various element classes of WeakTableau.
    """
    def shape(self):
        r"""
        Return the shape of ``self``.

        When the tableau is straight, the outer shape is returned.
        When the tableau is skew, the tuple of the outer and inner shape is returned.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.shape()
            [5, 2, 1]
            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.shape()
            ([5, 2, 1], [2])

            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: t.shape()
            [3, 2, 1]
            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.shape()
            ([3, 2, 1], [2])

            sage: t = WeakTableau([[2],[0,3],[2,1,0]], 3, representation = 'factorized_permutation')
            sage: t.shape()
            [5, 2, 1]
            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.shape()
            ([5, 2, 1], [2])
        """
        return self.parent().shape()

    def weight(self):
        r"""
        Return the weight of ``self``.

        The weight is a tuple whose `i`-th entry is the number of labels `i` in the
        bounded representation of ``self``.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.weight()
            (2, 2, 2)
            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.weight()
            (1, 1, 1, 1)
            sage: t = WeakTableau([[None,2,3],[3]],2)
            sage: t.weight()
            (0, 1, 1)

            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: t.weight()
            (3, 2, 1)
            sage: t = WeakTableau([[1,1,2],[2,3],[3]], 3, representation = 'bounded')
            sage: t.weight()
            (2, 2, 2)
            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.weight()
            (1, 1, 1, 1)

            sage: t = WeakTableau([[2],[0,3],[2,1,0]], 3, representation = 'factorized_permutation')
            sage: t.weight()
            (3, 2, 1)
            sage: t = WeakTableau([[2,0],[3,2],[1,0]], 3, representation = 'factorized_permutation')
            sage: t.weight()
            (2, 2, 2)
            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.weight()
            (2, 2)
        """
        return self.parent()._weight

    def size(self):
        r"""
        Return the size of the shape of ``self``.

        In the bounded representation, the size of the shape is the number of boxes in the
        outer shape minus the number of boxes in the inner shape. For the core and
        factorized permutation representation, the size is the length of the outer shape
        minus the length of the inner shape.

        .. SEEALSO:: :meth:`sage.combinat.core.Core.length`

        EXAMPLES::

            sage: t = WeakTableau([[None, 1, 1, 2, 2], [None, 2], [1]], 3)
            sage: t.shape()
            ([5, 2, 1], [1, 1])
            sage: t.size()
            4
            sage: t = WeakTableau([[1,1,2],[2,3],[3]], 3, representation="bounded")
            sage: t.shape()
            [3, 2, 1]
            sage: t.size()
            6
        """
        return self.parent().size()

    def intermediate_shapes(self):
        r"""
        Return the intermediate shapes of ``self``.

        A (skew) tableau with letters `1,2,\ldots,\ell` can be viewed as a sequence of shapes,
        where the `i`-th shape is given by the shape of the subtableau on letters `1,2,\ldots,i`.
        The output is the list of these shapes.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            sage: t.intermediate_shapes()
            [[], [2], [4, 1], [5, 2, 1]]

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.intermediate_shapes()
            [[2], [2, 1], [3, 1, 1], [4, 1, 1], [5, 2, 1]]

            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: t.intermediate_shapes()
            [[], [3], [3, 2], [3, 2, 1]]

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.intermediate_shapes()
            [[2], [3], [3, 1], [3, 1, 1], [3, 2, 1]]

            sage: t = WeakTableau([[0],[3],[2],[3]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.intermediate_shapes()
            [[2], [2, 1], [3, 1, 1], [4, 1, 1], [5, 2, 1]]
        """
        if self.parent()._representation in ['core', 'bounded']:
            return intermediate_shapes(self)
        else:
            return intermediate_shapes(self.to_core_tableau())

    def pp(self):
        r"""
        Return a pretty print string of the tableau.

        EXAMPLES::

            sage: t = WeakTableau([[None, 1, 1, 2, 2], [None, 2], [1]], 3)
            sage: t.pp()
            .  1  1  2  2
            .  2
            1
            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.pp()
            [s2*s0, s3*s2]
        """
        if self.parent()._representation in ['core', 'bounded']:
            print(self._repr_diagram())
        else:
            print(self)

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        EXAMPLES::

            sage: T = WeakTableaux(3, [5,2,1], [1,1,1,1,1,1], representation='core')
            sage: t = T[0]
            sage: hash(t) == hash(t)
            True
            sage: T = WeakTableaux(3, [2,2,1], [1,1,1,1,1], representation='bounded')
            sage: t = T[0]
            sage: hash(t) == hash(t)
            True
            sage: T = WeakTableaux(3, [5,2,1], [1,1,1,1,1,1], representation='factorized_permutation')
            sage: t = T[0]
            sage: hash(t) == hash(t)
            True
        """
        if self.parent()._representation in ['core', 'bounded']:
            return hash(tuple(tuple(x) for x in self)) + hash(self.parent().k)
        else:
            return super(WeakTableau_abstract, self).__hash__()

    def _latex_(self):
        r"""
        Return a latex method for the tableau.

        EXAMPLES::

            sage: t = WeakTableau([[None, 1, 1, 2, 2], [None, 2], [1]], 3)
            sage: latex(t)
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{5}c}\cline{1-5}
            \lr{}&\lr{1}&\lr{1}&\lr{2}&\lr{2}\\\cline{1-5}
            \lr{}&\lr{2}\\\cline{1-2}
            \lr{1}\\\cline{1-1}
            \end{array}$}
            }

            sage: t = WeakTableau([[0,3],[2,1]], 3, inner_shape = [1,1], representation = 'factorized_permutation')
            sage: latex(t)
            [s_{0}s_{3},s_{2}s_{1}]
        """
        def chi(x):
            if x is None:
                return ""
            if x in ZZ:
                return x
            return "%s"%x
        if self.parent()._representation in ['core', 'bounded']:
            t = [[chi(x) for x in row] for row in self]
            from .output import tex_from_array
            return tex_from_array(t)
        else:
            return "["+"".join(self[i]._latex_()+',' for i in range(len(self)-1))+self[len(self)-1]._latex_()+"]"

    def representation(self, representation = 'core'):
        r"""
        Return the analogue of ``self`` in the specified representation.

        INPUT:

        - ``representation`` -- 'core', 'bounded', or 'factorized_permutation' (default: 'core')

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            sage: t.parent()._representation
            'core'
            sage: t.representation('bounded')
            [[1, 1, 2, 4], [2, 3, 5], [3, 4], [5, 6], [6], [7]]
            sage: t.representation('factorized_permutation')
            [s0, s3*s1, s2*s1, s0*s4, s3*s0, s4*s2, s1*s0]

            sage: tb = WeakTableau([[1, 1, 2, 4], [2, 3, 5], [3, 4], [5, 6], [6], [7]], 4, representation = 'bounded')
            sage: tb.parent()._representation
            'bounded'
            sage: tb.representation('core') == t
            True
            sage: tb.representation('factorized_permutation')
            [s0, s3*s1, s2*s1, s0*s4, s3*s0, s4*s2, s1*s0]

            sage: tp = WeakTableau([[0],[3,1],[2,1],[0,4],[3,0],[4,2],[1,0]], 4, representation = 'factorized_permutation')
            sage: tp.parent()._representation
            'factorized_permutation'
            sage: tp.representation('core') == t
            True
            sage: tp.representation('bounded') == tb
            True
        """
        t = self
        if self.parent()._representation in ['bounded', 'factorized_permutation']:
            t = t.to_core_tableau()
        if representation == 'core':
            return t
        elif representation == 'bounded':
            return t.to_bounded_tableau()
        elif representation == 'factorized_permutation':
            return t.to_factorized_permutation_tableau()
        else:
            raise ValueError("The representation must be one of 'core', 'bounded', or 'factorized_permutation'")

#Abstract class for the parents of weak tableaux
class WeakTableaux_abstract(UniqueRepresentation, Parent):
    r"""
    Abstract class for the various parent classes of WeakTableaux.
    """
    def shape(self):
        r"""
        Return the shape of the tableaux of ``self``.

        When ``self`` is the class of straight tableaux, the outer shape is returned.
        When ``self`` is the class of skew tableaux, the tuple of the outer and inner
        shape is returned.

        Note that in the 'core' and 'factorized_permutation' representation, the shapes
        are `(k+1)`-cores.  In the 'bounded' representation, the shapes are `k`-bounded
        partitions.

        If the user wants to access the skew shape (even if the inner shape is empty),
        please use ``self._shape``.

        EXAMPLES::

            sage: T = WeakTableaux(3, [5,2,2], [2,2,2,1])
            sage: T.shape()
            [5, 2, 2]
            sage: T._shape
            ([5, 2, 2], [])
            sage: T = WeakTableaux(3, [[5,2,2], [1]], [2,1,2,1])
            sage: T.shape()
            ([5, 2, 2], [1])

            sage: T = WeakTableaux(3, [3,2,2], [2,2,2,1], representation = 'bounded')
            sage: T.shape()
            [3, 2, 2]
            sage: T._shape
            ([3, 2, 2], [])
            sage: T = WeakTableaux(3, [[3,2,2], [1]], [2,1,2,1], representation = 'bounded')
            sage: T.shape()
            ([3, 2, 2], [1])

            sage: T = WeakTableaux(3, [4,1], [2,2], representation = 'factorized_permutation')
            sage: T.shape()
            [4, 1]
            sage: T._shape
            ([4, 1], [])
            sage: T = WeakTableaux(4, [[6,2,1], [2]], [2,1,1,1], representation = 'factorized_permutation')
            sage: T.shape()
            ([6, 2, 1], [2])
        """
        if self._skew:
            return (self._outer_shape, self._inner_shape)
        return self._outer_shape

    def size(self):
        r"""
        Return the size of the shape.

        In the bounded representation, the size of the shape is the number of boxes in the
        outer shape minus the number of boxes in the inner shape. For the core and
        factorized permutation representation, the size is the length of the outer shape
        minus the length of the inner shape.

        EXAMPLES::

            sage: T = WeakTableaux(3, [5,2,1], [1,1,1,1,1,1])
            sage: T.size()
            6
            sage: T = WeakTableaux(3, [3,2,1], [1,1,1,1,1,1], representation = 'bounded')
            sage: T.size()
            6
            sage: T = WeakTableaux(4, [[6,2,1], [2]], [2,1,1,1], 'factorized_permutation')
            sage: T.size()
            5
        """
        if self._representation == 'bounded':
            return self._outer_shape.size() - self._inner_shape.size()
        else:
            return self._outer_shape.length() - self._inner_shape.length()

    def representation(self, representation = 'core'):
        r"""
        Return the analogue of ``self`` in the specified representation.

        INPUT:

        - ``representation`` -- 'core', 'bounded', or 'factorized_permutation' (default: 'core')

        EXAMPLES::

            sage: T = WeakTableaux(3, [5,2,1], [1,1,1,1,1,1])
            sage: T._representation
            'core'
            sage: T.representation('bounded')
            Bounded weak 3-Tableaux of (skew) 3-bounded shape [3, 2, 1] and weight (1, 1, 1, 1, 1, 1)
            sage: T.representation('factorized_permutation')
            Factorized permutation (skew) weak 3-Tableaux of shape [5, 2, 1] and weight (1, 1, 1, 1, 1, 1)

            sage: T = WeakTableaux(3, [3,2,1], [1,1,1,1,1,1], representation = 'bounded')
            sage: T._representation
            'bounded'
            sage: T.representation('core')
            Core weak 3-Tableaux of (skew) core shape [5, 2, 1] and weight (1, 1, 1, 1, 1, 1)
            sage: T.representation('bounded')
            Bounded weak 3-Tableaux of (skew) 3-bounded shape [3, 2, 1] and weight (1, 1, 1, 1, 1, 1)
            sage: T.representation('bounded') == T
            True
            sage: T.representation('factorized_permutation')
            Factorized permutation (skew) weak 3-Tableaux of shape [5, 2, 1] and weight (1, 1, 1, 1, 1, 1)
            sage: T.representation('factorized_permutation') == T
            False

            sage: T = WeakTableaux(3, [5,2,1], [1,1,1,1,1,1], representation = 'factorized_permutation')
            sage: T._representation
            'factorized_permutation'
            sage: T.representation('core')
            Core weak 3-Tableaux of (skew) core shape [5, 2, 1] and weight (1, 1, 1, 1, 1, 1)
            sage: T.representation('bounded')
            Bounded weak 3-Tableaux of (skew) 3-bounded shape [3, 2, 1] and weight (1, 1, 1, 1, 1, 1)
            sage: T.representation('factorized_permutation')
            Factorized permutation (skew) weak 3-Tableaux of shape [5, 2, 1] and weight (1, 1, 1, 1, 1, 1)
        """
        outer_shape = self._outer_shape
        inner_shape = self._inner_shape
        weight = self._weight
        if (self._representation in ['core', 'factorized_permutation']) and representation == 'bounded':
            outer_shape = outer_shape.to_bounded_partition()
            inner_shape = inner_shape.to_bounded_partition()
        if self._representation == 'bounded' and (representation in ['core', 'factorized_permutation']):
            outer_shape = outer_shape.to_core(self.k)
            inner_shape = inner_shape.to_core(self.k)
        return WeakTableaux(self.k, [outer_shape, inner_shape], weight, representation = representation)


#Weak Tableaux in terms of cores
class WeakTableau_core(WeakTableau_abstract):
    r"""
    A (skew) weak `k`-tableau represented in terms of `(k+1)`-cores.
    """
    @staticmethod
    def __classcall_private__(cls, t, k):
        r"""
        Implements the shortcut ``WeakTableau_core(t, k)`` to ``WeakTableaux_core(k, shape , weight)(t)``
        where ``shape`` is the shape of the tableau and ``weight`` is its weight.

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableau_core
            sage: t = WeakTableau_core([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.check()
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>
            sage: TestSuite(t).run()
            sage: t.parent()._skew
            False

            sage: t = WeakTableau_core([[None, None, 1, 1, 2], [1, 2], [2]],3)
            sage: t.check()
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>
            sage: TestSuite(t).run()
            sage: t.parent()._skew
            True
        """
        if isinstance(t, cls):
            return t
        tab = SkewTableau(list(t))
        outer = Core(tab.outer_shape(),k+1)
        inner = Core(tab.inner_shape(),k+1)
        weight = WeakTableau_bounded.from_core_tableau(t,k).weight()
        return WeakTableaux_core(k, [outer, inner], weight)(t)

    def __init__(self, parent, t):
        r"""
        Initialization of weak `k`-tableau ``t`` in core representation.

        INPUT:

        - ``t`` -- weak tableau in core representation; the input is supposed to be a list
          of lists specifying the rows of the tableau;
          ``None`` is allowed as an entry for skew weak `k`-tableaux

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableau_core, WeakTableaux_core
            sage: T = WeakTableaux_core(3,[5,2,1],[2,2,2])
            sage: t = T([[1, 1, 2, 2, 3], [2, 3], [3]]); t
            [[1, 1, 2, 2, 3], [2, 3], [3]]
            sage: c = WeakTableau_core([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            sage: T = WeakTableaux_core(3,[5,2,1],[2,2,2])
            sage: t = T([[1, 1, 2, 2, 3], [2, 3], [3]]); t
            [[1, 1, 2, 2, 3], [2, 3], [3]]
            sage: c == t
            True
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>
            sage: t.parent()
            Core weak 3-Tableaux of (skew) core shape [5, 2, 1] and weight (2, 2, 2)
            sage: TestSuite(t).run()

            sage: t = WeakTableau_core([[None, None, 1, 1, 2], [1, 2], [2]],3);  t
            [[None, None, 1, 1, 2], [1, 2], [2]]
            sage: t.weight()
            (2, 2)
            sage: t.shape()
            ([5, 2, 1], [2])
            sage: TestSuite(t).run()
        """
        self.k = parent.k
        ClonableList.__init__(self, parent, t)

    def _repr_diagram(self):
        r"""
        Return a string representation of ``self`` as a diagram.

        EXAMPLES::

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: print(t._repr_diagram())
            .  .  2  3  4
            1  4
            2
        """
        t = SkewTableau(list(self))
        return t._repr_diagram()

    def shape_core(self):
        r"""
        Return the shape of ``self`` as a `(k+1)`-core.

        When the tableau is straight, the outer shape is returned as a core.  When the
        tableau is skew, the tuple of the outer and inner shape is returned as cores.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            sage: t.shape_core()
            [5, 2, 1]

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.shape_core()
            ([5, 2, 1], [2])
        """
        return self.shape()

    def shape_bounded(self):
        r"""
        Return the shape of ``self`` as a `k`-bounded partition.

        When the tableau is straight, the outer shape is returned as a `k`-bounded
        partition.  When the tableau is skew, the tuple of the outer and inner shape is
        returned as `k`-bounded partitions.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            sage: t.shape_bounded()
            [3, 2, 1]

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.shape_bounded()
            ([3, 2, 1], [2])
        """
        if self.parent()._skew:
            return tuple([r.to_bounded_partition() for r in self.shape_core()])
        return self.shape_core().to_bounded_partition()

    def check(self):
        r"""
        Check that ``self`` is a valid weak `k`-tableau.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2], [2]], 2)
            sage: t.check()
            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.check()

        TESTS::

            sage: T = WeakTableaux(2, [3,1], [1,1,1,1])
            sage: t = T([[1,2,3],[3]])
            Traceback (most recent call last):
            ...
            ValueError: The weight of the parent does not agree with the weight of the tableau!

            sage: t = WeakTableau([[1, 2, 2], [1]], 2)
            Traceback (most recent call last):
            ...
            ValueError: The tableau is not semistandard!
        """
        if not self.parent()._weight == WeakTableau_bounded.from_core_tableau(self,self.k).weight():
            raise ValueError("The weight of the parent does not agree with the weight of the tableau!")
        t = SkewTableau(list(self))
        if t not in SemistandardSkewTableaux():
            raise ValueError("The tableau is not semistandard!")
        outer = Core(t.outer_shape(),self.k+1)
        inner = Core(t.inner_shape(),self.k+1)
        if self.parent()._outer_shape != outer:
            raise ValueError("The outer shape of the parent does not agree with the outer shape of the tableau!")
        if self.parent()._inner_shape != inner:
            raise ValueError("The inner shape of the parent does not agree with the inner shape of the tableau!")
        self.to_bounded_tableau().check()

    def to_bounded_tableau(self):
        r"""
        Return the bounded representation of the weak `k`-tableau ``self``.

        Each restricted subtableau of the output is a `k`-bounded partition.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: c = t.to_bounded_tableau(); c
            [[1, 1, 2], [2, 3], [3]]
            sage: type(c)
            <class 'sage.combinat.k_tableau.WeakTableaux_bounded_with_category.element_class'>

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.to_bounded_tableau()
            [[None, None, 3], [1, 4], [2]]
            sage: t.to_bounded_tableau().to_core_tableau() == t
            True
        """
        shapes = [ Core(p,self.k+1).to_bounded_partition() for p in self.intermediate_shapes() ]
        if self.parent()._skew:
            l = [[None]*i for i in shapes[0]]
        else:
            l = []
        for i in range(1,len(shapes)):
            p = shapes[i]
            if len(l) < len(p):
                l += [[]]
            l_new = []
            for j in range(len(l)):
                l_new += [l[j] + [i]*(p[j]-len(l[j]))]
            l = l_new
        return WeakTableau_bounded(l, self.k)

    def to_factorized_permutation_tableau(self):
        r"""
        Return the factorized permutation representation of the weak `k`-tableau ``self``.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: c = t.to_factorized_permutation_tableau(); c
            [s2*s0, s3*s2, s1*s0]
            sage: type(c)
            <class 'sage.combinat.k_tableau.WeakTableaux_factorized_permutation_with_category.element_class'>
            sage: c.to_core_tableau() == t
            True

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: c = t.to_factorized_permutation_tableau(); c
            [s0, s3, s2, s3]
            sage: c._inner_shape
            [2]
            sage: c.to_core_tableau() == t
            True

        TESTS::

            sage: t = WeakTableau([], 4)
            sage: c = t.to_factorized_permutation_tableau(); c
            [1]
            sage: c._inner_shape
            []
            sage: c.to_core_tableau() == t
            True
        """
        shapes = [ Core(p,self.k+1).to_grassmannian() for p in self.intermediate_shapes() ]
        perms = [ shapes[i]*(shapes[i-1].inverse()) for i in range(len(shapes)-1,0,-1)]
        return WeakTableau_factorized_permutation(perms, self.k, inner_shape = self.parent()._inner_shape)

    def residues_of_entries(self, v):
        r"""
        Return a list of residues of cells of weak `k`-tableau ``self`` labeled by ``v``.

        INPUT:

        - ``v`` -- a label of a cell in ``self``

        OUTPUT:

        - a list of residues

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            sage: t.residues_of_entries(1)
            [0, 1]

            sage: t = WeakTableau([[None, None, 1, 1, 4], [1, 4], [3]], 3)
            sage: t.residues_of_entries(1)
            [2, 3]
        """
        S = set((j - i) % (self.k+1)
                for i in range(len(self))
                for j in range(len(self[i]))
                if self[i][j] == v)
        return sorted(S)

    def dictionary_of_coordinates_at_residues(self, v):
        r"""
        Return a dictionary assigning to all residues of ``self`` with label ``v`` a list
        of cells with the given residue.

        INPUT:

        - ``v`` -- a label of a cell in ``self``

        OUTPUT:

        - dictionary assigning coordinates in ``self`` to residues

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            sage: t.dictionary_of_coordinates_at_residues(3)
            {0: [(0, 4), (1, 1)], 2: [(2, 0)]}

            sage: t = WeakTableau([[None, None, 1, 1, 4], [1, 4], [3]], 3)
            sage: t.dictionary_of_coordinates_at_residues(1)
            {2: [(0, 2)], 3: [(0, 3), (1, 0)]}

            sage: t = WeakTableau([], 3)
            sage: t.dictionary_of_coordinates_at_residues(1)
            {}
        """
        d = {}
        for r in self.residues_of_entries(v):
            d[r] = []
            for i in range(len(self)):
                for j in range(len(self[i])):
                    if self[i][j]==v and (j - i)%(self.k+1) == r:
                        d[r]+=[(i,j)]
        return d

    def list_of_standard_cells(self):
        r"""
        Return a list of lists of the coordinates of the standard cells of ``self``.

        INPUT:

        - ``self`` -- a weak `k`-tableau in core representation with partition weight

        OUTPUT:

        - a list of lists of coordinates

        .. WARNING::

            This method currently only works for straight weak tableaux with partition
            weight.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.list_of_standard_cells()
            [[(0, 1), (1, 0), (2, 0)], [(0, 0), (0, 2), (1, 1)]]
            sage: t = WeakTableau([[1, 1, 1, 2], [2, 2, 3]], 5)
            sage: t.list_of_standard_cells()
            [[(0, 2), (1, 1), (1, 2)], [(0, 1), (1, 0)], [(0, 0), (0, 3)]]
            sage: t = WeakTableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            sage: t.list_of_standard_cells()
            [[(0, 1), (1, 0), (2, 0), (0, 5), (3, 0), (4, 0), (5, 0)], [(0, 0), (0, 2), (1, 1), (2, 1), (1, 2), (3, 1)]]

        TESTS::

            sage: t = WeakTableau([],3)
            sage: t.list_of_standard_cells()
            []

            sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            sage: t.list_of_standard_cells()
            Traceback (most recent call last):
            ...
            ValueError: This method only works for straight tableaux!

            sage: t = WeakTableau([[1,2],[2]], 3)
            sage: t.list_of_standard_cells()
            Traceback (most recent call last):
            ...
            ValueError: This method only works for weak tableaux with partition weight!
        """
        if self.parent()._skew:
            raise ValueError("This method only works for straight tableaux!")
        if self.weight() not in Partitions(sum(self.weight())):
            raise ValueError("This method only works for weak tableaux with partition weight!")
        if not self:
            return []
        mu = Partition(self.weight()).conjugate()
        already_used = []
        out = []
        for i in range(self[0].count(1)):
            standard_cells = [(0,self[0].count(1) - i - 1)]
            r = self[0].count(1) - i - 1
            for v in range(1,mu[i]):
                D = self.dictionary_of_coordinates_at_residues(v+1)
                new_D = {a: b for (a, b) in D.items()
                         if all(x not in already_used for x in b)}
                r = (r - min([self.k+1 - (x-r)%(self.k+1) for x in new_D]))%(self.k+1)
                standard_cells.append(new_D[r][-1])
                already_used += new_D[r]
            out.append(standard_cells)
        return out

    def k_charge(self, algorithm = "I"):
        r"""
        Return the `k`-charge of ``self``.

        INPUT:

        - ``algorithm`` -- (default: "I") if "I", computes `k`-charge using the `I`
          algorithm, otherwise uses the `J`-algorithm

        OUTPUT:

        - a nonnegative integer

        For the definition of `k`-charge and the various algorithms to compute it see
        Section 3.3 of [LLMSSZ2013]_.

        .. SEEALSO:: :meth:`k_charge_I` and :meth:`k_charge_J`

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.k_charge()
            2
            sage: t = WeakTableau([[1, 3, 4, 5, 6], [2, 6], [4]], 3)
            sage: t.k_charge()
            8
            sage: t = WeakTableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            sage: t.k_charge()
            12

        TESTS::

            sage: T = WeakTableaux(4, [13,9,5,3,2,1,1], [4,3,3,2,2,1,1,1])
            sage: T.cardinality()
            6
            sage: all(t.k_charge_I() == t.k_charge_J() for t in T)
            True
        """
        if algorithm == "I":
            return self.k_charge_I()
        return self.k_charge_J()

    def k_charge_I(self):
        r"""
        Return the `k`-charge of ``self`` using the `I`-algorithm.

        For the definition of `k`-charge and the `I`-algorithm see Section 3.3 of [LLMSSZ2013]_.

        OUTPUT:

        - a nonnegative integer

        .. SEEALSO:: :meth:`k_charge` and :meth:`k_charge_J`

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.k_charge_I()
            2
            sage: t = WeakTableau([[1, 3, 4, 5, 6], [2, 6], [4]], 3)
            sage: t.k_charge_I()
            8
            sage: t = WeakTableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            sage: t.k_charge_I()
            12

        TESTS::

            sage: t = WeakTableau([[None, None, 1, 1, 4], [1, 4], [3]], 3)
            sage: t.k_charge_I()
            Traceback (most recent call last):
            ...
            ValueError: k-charge is not defined for skew weak tableaux
        """
        if self.parent()._skew:
            raise ValueError("k-charge is not defined for skew weak tableaux")
        stt = self.list_of_standard_cells()
        kch = 0
        for sw in stt:
            Ii = 0
            for r in range(len(sw)-1):
                if sw[r][1] < sw[r+1][1]:
                    Ii += 1 + abs(self.parent().diag(sw[r+1],sw[r]))
                else:
                    Ii += - abs(self.parent().diag(sw[r],sw[r+1]))
                kch += Ii
        return kch

    def k_charge_J(self):
        r"""
        Return the `k`-charge of ``self`` using the `J`-algorithm.

        For the definition of `k`-charge and the `J`-algorithm see Section 3.3 of [LLMSSZ2013]_.

        OUTPUT:

        - a nonnegative integer

        .. SEEALSO:: :meth:`k_charge` and :meth:`k_charge_I`

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: t.k_charge_J()
            2
            sage: t = WeakTableau([[1, 3, 4, 5, 6], [2, 6], [4]], 3)
            sage: t.k_charge_J()
            8
            sage: t = WeakTableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            sage: t.k_charge_J()
            12

        TESTS::

            sage: t = WeakTableau([[None, None, 1, 1, 4], [1, 4], [3]], 3)
            sage: t.k_charge_I()
            Traceback (most recent call last):
            ...
            ValueError: k-charge is not defined for skew weak tableaux

            sage: t = WeakTableau([[1, 1, 2, 3], [2, 4, 4], [3, 6], [5]], 4, representation='bounded')
            sage: t.k_charge() == t.k_charge(algorithm = 'J')
            True
        """
        if self.parent()._skew:
            raise ValueError("k-charge is not defined for skew weak tableaux")
        stt = self.list_of_standard_cells()
        kch = 0
        for sw in stt:
            Ji = 0
            for i in range(len(sw)-1):
                c = (self._height_of_restricted_subword(sw,i+2)+1,0)
                cdi = self.parent().circular_distance((-c[0])%(self.k+1),(sw[i][1]-sw[i][0])%(self.k+1))
                cdi1 = self.parent().circular_distance((-c[0])%(self.k+1),(sw[i+1][1]-sw[i+1][0])%(self.k+1))
                if (cdi > cdi1):
                    Ji += 1
                kch += Ji + self.parent().diag(sw[i+1],c)
        return kch

    def _height_of_restricted_subword(self, sw, r):
        r"""
        Return the row of the highest addable cell of the subtableau of ``self`` with letters `\le r`
        (excluding letters `r` in standard subwords before ``sw``).

        Restrict the weak `k`-tableau ``self`` to letters `\le r` and remove all letters
        `r` that appeared in a previous standard subword selected by
        :meth:`list_of_standard_cells`.

        INPUT:

        - ``sw`` -- one of the subwords of standard cells of ``self``
        - ``r`` -- nonnegative integer

        OUTPUT:

        - a nonnegative integer

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            sage: s = t.list_of_standard_cells()[0]; s
            [(0, 1), (1, 0), (2, 0)]
            sage: t._height_of_restricted_subword(s,2)
            1

            sage: t = WeakTableau([[1, 3, 4, 5, 6], [2, 6], [4]], 3)
            sage: s = t.list_of_standard_cells()[0]; s
            [(0, 0), (1, 0), (0, 1), (2, 0), (0, 3), (1, 1)]
            sage: t._height_of_restricted_subword(s,4)
            2

            sage: t = WeakTableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            sage: s = t.list_of_standard_cells()[0]; s
            [(0, 1), (1, 0), (2, 0), (0, 5), (3, 0), (4, 0), (5, 0)]
            sage: t._height_of_restricted_subword(s,6)
            4
        """
        R = [v for v in self.shape().to_partition().cells() if self[v[0]][v[1]] < r]
        L = [v for v in sw if self[v[0]][v[1]] <= r]
        return max(v[0] for v in L+R)

class WeakTableaux_core(WeakTableaux_abstract):
    r"""
    The class of (skew) weak `k`-tableaux in the core representation of shape ``shape``
    (as `k+1`-core) and weight ``weight``.

    INPUT:

    - ``k`` -- positive integer
    - ``shape`` -- the shape of the `k`-tableaux represented as a `(k+1)`-core; if the
      tableaux are skew, the shape is a tuple of the outer and inner shape (both as
      `(k+1)`-cores)
    - ``weight`` -- the weight of the `k`-tableaux

    EXAMPLES::

        sage: T = WeakTableaux(3, [4,1], [2,2])
        sage: T.list()
        [[[1, 1, 2, 2], [2]]]

        sage: T = WeakTableaux(3, [[5,2,1], [2]], [1,1,1,1])
        sage: T.list()
        [[[None, None, 2, 3, 4], [1, 4], [2]],
        [[None, None, 1, 2, 4], [2, 4], [3]],
        [[None, None, 1, 2, 3], [2, 3], [4]]]
    """

    @staticmethod
    def __classcall_private__(cls, k, shape, weight):
        r"""
        Straighten arguments before unique representation.

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_core
            sage: T = WeakTableaux_core(3, [2,1], [1,1,1])
            sage: TestSuite(T).run()
            sage: T = WeakTableaux_core(3, [[5,2,1], [2]], [1,1,1,1])
            sage: TestSuite(T).run()
        """
        if shape == [] or shape[0] in ZZ:
            shape = (Core(shape, k+1), Core([],k+1))
        else:
            shape = tuple([Core(r,k+1) for r in shape])
        return super(WeakTableaux_core, cls).__classcall__(cls, k, shape, tuple(weight))

    def __init__(self, k, shape, weight):
        r"""
        Initializes the parent class of (skew) weak `k`-tableaux in core representation.

        INPUT:

        - ``k`` -- positive integer
        - ``outer_shape`` -- the outer shape of the `k`-tableaux represented as a
          `(k+1)`-core
        - ``weight`` -- the weight of the `k`-tableaux
        - ``inner_shape`` -- the inner shape of the skew `k`-tableaux represented as a
          `(k+1)`-core;  for straight tableaux the inner shape does not need to be
          specified (default: [])

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_core
            sage: T = WeakTableaux_core(3, [4,1], [2,2])
            sage: TestSuite(T).run()
            sage: T = WeakTableaux_core(3, [[5,2,1], [2]], [1,1,1,1])
            sage: TestSuite(T).run()
        """
        self.k = k
        self._skew = bool(shape[1])
        self._outer_shape = shape[0]
        self._inner_shape = shape[1]
        self._shape = (self._outer_shape, self._inner_shape)
        self._weight = weight
        self._representation = 'core'
        Parent.__init__(self, category = FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_core
            sage: repr(WeakTableaux_core(3, [2,1], [1,1,1]))
            'Core weak 3-Tableaux of (skew) core shape [2, 1] and weight (1, 1, 1)'
            sage: repr(WeakTableaux_core(3, [[5,2,1], [2]], [1,1,1,1]))
            'Core weak 3-Tableaux of (skew) core shape ([5, 2, 1], [2]) and weight (1, 1, 1, 1)'
        """
        return "Core weak %s-Tableaux of (skew) core shape %s and weight %s"%(self.k, self.shape(), self._weight)

    def __iter__(self):
        r"""
        TESTS::

            sage: T = WeakTableaux(3, [4,1], [2,2])
            sage: T.list()
            [[[1, 1, 2, 2], [2]]]
            sage: T = WeakTableaux(3, [5,2,2], [2,2,2,1])
            sage: T.list()
            [[[1, 1, 3, 3, 4], [2, 2], [3, 3]], [[1, 1, 2, 2, 3], [2, 3], [3, 4]]]
            sage: T = WeakTableaux(3, [[5,2,2], [1]], [2,1,2,1])
            sage: T.list()
            [[[None, 1, 3, 3, 4], [1, 2], [3, 3]],
            [[None, 1, 2, 3, 3], [1, 3], [2, 4]],
            [[None, 1, 1, 2, 3], [2, 3], [3, 4]]]
        """
        for t in WeakTableaux_bounded(self.k, [self._outer_shape.to_bounded_partition(), self._inner_shape.to_bounded_partition()], self._weight):
            yield t.to_core_tableau()

    def diag(self, c, ha):
        r"""
        Return the number of diagonals strictly between cells ``c`` and ``ha`` of the same residue as ``c``.

        INPUT:

        - ``c`` -- a cell in the lattice
        - ``ha`` -- another cell in the lattice with bigger row and smaller column than `c`

        OUTPUT:

        - a nonnegative integer

        EXAMPLES::

            sage: T = WeakTableaux(4, [5,2,2], [2,2,2,1])
            sage: T.diag((1,2),(4,0))
            0
        """
        return divmod((c[1]-c[0])-(ha[1]-ha[0])-1, self.k+1)[0]

    def circular_distance(self, cr, r):
        r"""
        Return the shortest counterclockwise distance between ``cr`` and ``r`` modulo `k+1`.

        INPUT:

        - ``cr``, ``r`` -- nonnegative integers between `0` and `k`

        OUTPUT:

        - a positive integer

        EXAMPLES::

            sage: T = WeakTableaux(10, [], [])
            sage: T.circular_distance(8, 6)
            2
            sage: T.circular_distance(8, 8)
            0
            sage: T.circular_distance(8, 9)
            10
        """
        return self.k - ((r+self.k-cr)%(self.k+1))

    Element = WeakTableau_core


#Weak tableaux in terms of `k`-bounded partitions
class WeakTableau_bounded(WeakTableau_abstract):
    r"""
    A (skew) weak `k`-tableau represented in terms of `k`-bounded partitions.
    """
    @staticmethod
    def __classcall_private__(cls, t, k):
        r"""
        Implements the shortcut ``WeakTableau_bounded(t, k)`` to ``WeakTableaux_bounded(k, shape, weight)(t)``
        where ``shape`` is the shape of the tableau and ``weight`` is its weight.

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableau_bounded
            sage: t = WeakTableau_bounded([[1,1,2],[2,3],[3]],3)
            sage: t.check()
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_bounded_with_category.element_class'>
            sage: TestSuite(t).run()
            sage: t.parent()._skew
            False

            sage: t = WeakTableau_bounded([[None, None, 1], [1, 2], [2]], 3)
            sage: t.check()
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_bounded_with_category.element_class'>
            sage: TestSuite(t).run()
            sage: t.parent()._skew
            True
        """
        if isinstance(t, cls):
            return t
        tab = SkewTableau(list(t))
        outer = tab.outer_shape()
        inner = tab.inner_shape()
        weight = tuple(tab.weight())
        if outer.conjugate().length() > k:
            raise ValueError("The shape of %s is not %s-bounded"%(t, k))
        return WeakTableaux_bounded(k, [outer, inner], weight)(t)

    def __init__(self, parent, t):
        r"""
        Initialization of (skew) weak `k`-tableau ``t`` in `k`-bounded representation.

        INPUT:

        - ``t`` -- weak tableau in `k`-bounded representation; the input is supposed to be
          a list of iterables specifying the rows of the tableau;  ``None`` is allowed as an
          entry for skew weak `k`-tableaux

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableau_bounded, WeakTableaux_bounded
            sage: c = WeakTableau_bounded([[1,1,2],[2,3],[3]],3)
            sage: T = WeakTableaux_bounded(3,[3,2,1],[2,2,2])
            sage: t = T([[1,1,2],[2,3],[3]]); t
            [[1, 1, 2], [2, 3], [3]]
            sage: c == t
            True
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_bounded_with_category.element_class'>
            sage: t.parent()
            Bounded weak 3-Tableaux of (skew) 3-bounded shape [3, 2, 1] and weight (2, 2, 2)
            sage: TestSuite(t).run()

            sage: t = WeakTableau_bounded([[None, None, 1], [2, 4], [3]], 3)
            sage: t.shape()
            ([3, 2, 1], [2])
            sage: t.weight()
            (1, 1, 1, 1)
            sage: TestSuite(t).run()

            sage: t = T([[1,1,3],[2,2],[3]])
            Traceback (most recent call last):
            ...
            ValueError: This is not a proper weak 3-tableau
        """
        k = parent.k
        self.k = k
        if parent._outer_shape.conjugate().length() > k:
            raise ValueError("%s is not a %s-bounded tableau"%(t, k))
        ClonableList.__init__(self, parent, [list(r) for r in t])

    def _repr_diagram(self):
        r"""
        Return a string representation of ``self`` as a diagram.

        EXAMPLES::

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: print(t._repr_diagram())
            .  .  1
            2  4
            3
            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: print(t._repr_diagram())
            1  1  1
            2  2
            3
        """
        t = SkewTableau(list(self))
        return t._repr_diagram()

    def shape_core(self):
        r"""
        Return the shape of ``self`` as `(k+1)`-core.

        When the tableau is straight, the outer shape is returned as a `(k+1)`-core.
        When the tableau is skew, the tuple of the outer and inner shape is returned as
        `(k+1)`-cores.

        EXAMPLES::

            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: t.shape_core()
            [5, 2, 1]

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.shape_core()
            ([5, 2, 1], [2])
        """
        if self.parent()._skew:
            return tuple([r.to_core(self.k) for r in self.shape_bounded()])
        return self.shape_bounded().to_core(self.k)

    def shape_bounded(self):
        r"""
        Return the shape of ``self`` as `k`-bounded partition.

        When the tableau is straight, the outer shape is returned as a `k`-bounded
        partition.  When the tableau is skew, the tuple of the outer and inner shape is
        returned as `k`-bounded partitions.

        EXAMPLES::

            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: t.shape_bounded()
            [3, 2, 1]

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.shape_bounded()
            ([3, 2, 1], [2])
        """
        return self.shape()

    def check(self):
        r"""
        Check that ``self`` is a valid weak `k`-tableau.

        EXAMPLES::

            sage: t = WeakTableau([[1,1],[2]], 2, representation = 'bounded')
            sage: t.check()

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.check()

        TESTS::

            sage: t = WeakTableau([[1,1,3],[2,2],[3]], 3, representation = 'bounded')
            Traceback (most recent call last):
            ...
            ValueError: This is not a proper weak 3-tableau

            sage: T = WeakTableaux(3, [3,1], [2,1], representation = 'bounded')
            sage: t = T([[None, 1,1], [2]])
            Traceback (most recent call last):
            ...
            ValueError: The inner shape of the parent does not agree with the inner shape of the tableau!

            sage: t = WeakTableau([[1,1],[1]], 3, representation = 'bounded')
            Traceback (most recent call last):
            ...
            ValueError: The tableaux is not semistandard!
        """
        t = SkewTableau(list(self))
        if t not in SemistandardSkewTableaux():
            raise ValueError("The tableaux is not semistandard!")
        if not self.parent()._weight == tuple(t.weight()):
            raise ValueError("The weight of the parent does not agree with the weight of the tableau!")
        outer = t.outer_shape()
        inner = t.inner_shape()
        if self.parent()._outer_shape != outer:
            raise ValueError("The outer shape of the parent does not agree with the outer shape of the tableau!")
        if self.parent()._inner_shape != inner:
            raise ValueError("The inner shape of the parent does not agree with the inner shape of the tableau!")
        if not t.is_k_tableau(self.k):
            raise ValueError("This is not a proper weak %s-tableau"%(self.k))

    def _is_k_tableau(self):
        r"""
        Checks whether ``self`` is a valid weak `k`-tableau.

        EXAMPLES::

            sage: t = WeakTableau([[1,1,1],[2,2],[3]], 3, representation = 'bounded')
            sage: t._is_k_tableau()
            True

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t._is_k_tableau()
            True
        """
        shapes = self.intermediate_shapes()
        kshapes = [ la.k_conjugate(self.k) for la in shapes ]
        return all( kshapes[i+1].contains(kshapes[i]) for i in range(len(shapes)-1) )

    def to_core_tableau(self):
        r"""
        Return the weak `k`-tableau ``self`` where the shape of each restricted tableau is a `(k+1)`-core.

        EXAMPLES::

            sage: t = WeakTableau([[1,1,2,4],[2,3,5],[3,4],[5,6],[6],[7]], 4, representation = 'bounded')
            sage: c = t.to_core_tableau(); c
            [[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]]
            sage: type(c)
            <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>
            sage: t = WeakTableau([], 4, representation = 'bounded')
            sage: t.to_core_tableau()
            []

            sage: from sage.combinat.k_tableau import WeakTableau_bounded
            sage: t = WeakTableau([[1,1,2],[2,3],[3]], 3, representation = 'bounded')
            sage: WeakTableau_bounded.from_core_tableau(t.to_core_tableau(),3)
            [[1, 1, 2], [2, 3], [3]]
            sage: t == WeakTableau_bounded.from_core_tableau(t.to_core_tableau(),3)
            True

            sage: t = WeakTableau([[None, None, 1], [2, 4], [3]], 3, representation = 'bounded')
            sage: t.to_core_tableau()
            [[None, None, 1, 2, 4], [2, 4], [3]]
            sage: t == WeakTableau_bounded.from_core_tableau(t.to_core_tableau(),3)
            True
        """
        shapes = [ p.to_core(self.k) for p in self.intermediate_shapes() ]
        if self.parent()._skew:
            l = [[None]*i for i in shapes[0]]
        else:
            l = []
        for i in range(1,len(shapes)):
            p = shapes[i]
            if len(l) < len(p):
                l += [[]]
            l_new = []
            for j in range(len(l)):
                l_new += [l[j] + [i]*(p[j]-len(l[j]))]
            l = l_new
        return WeakTableau_core(l, self.k)

    @classmethod
    def from_core_tableau(cls, t, k):
        r"""
        Construct weak `k`-bounded tableau from in `k`-core tableau.

        EXAMPLES::

            sage: from sage.combinat.k_tableau import WeakTableau_bounded
            sage: WeakTableau_bounded.from_core_tableau([[1, 1, 2, 2, 3], [2, 3], [3]], 3)
            [[1, 1, 2], [2, 3], [3]]

            sage: WeakTableau_bounded.from_core_tableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
            [[None, None, 3], [1, 4], [2]]

            sage: WeakTableau_bounded.from_core_tableau([[None,2,3],[3]], 2)
            [[None, 2], [3]]
        """
        t = SkewTableau(list(t))
        shapes = [ Core(p, k+1).to_bounded_partition() for p in intermediate_shapes(t) ]#.to_chain() ]
        if t.inner_shape() == Partition([]):
            l = []
        else:
            l = [[None]*i for i in shapes[0]]
        for i in range(1,len(shapes)):
            p = shapes[i]
            if len(l) < len(p):
                l += [[]]
            l_new = []
            for j in range(len(l)):
                l_new += [l[j] + [i]*(p[j]-len(l[j]))]
            l = l_new
        return cls(l, k)

    def k_charge(self, algorithm = 'I'):
        r"""
        Return the `k`-charge of ``self``.

        INPUT:

        - ``algorithm`` -- (default: "I") if "I", computes `k`-charge using the `I`
          algorithm, otherwise uses the `J`-algorithm

        OUTPUT:

        - a nonnegative integer

        For the definition of `k`-charge and the various algorithms to compute it see Section 3.3 of [LLMSSZ2013]_.

        EXAMPLES::

            sage: t = WeakTableau([[1, 1, 2], [2, 3], [3]], 3, representation = 'bounded')
            sage: t.k_charge()
            2
            sage: t = WeakTableau([[1, 3, 5], [2, 6], [4]], 3, representation = 'bounded')
            sage: t.k_charge()
            8
            sage: t = WeakTableau([[1, 1, 2, 4], [2, 3, 5], [3, 4], [5, 6], [6], [7]], 4, representation = 'bounded')
            sage: t.k_charge()
            12
        """
        return self.to_core_tableau().k_charge(algorithm = algorithm)


class WeakTableaux_bounded(WeakTableaux_abstract):
    r"""
    The class of (skew) weak `k`-tableaux in the bounded representation of shape ``shape``
    (as `k`-bounded partition or tuple of `k`-bounded partitions in the skew case) and
    weight ``weight``.

    INPUT:

    - ``k`` -- positive integer
    - ``shape`` -- the shape of the `k`-tableaux represented as a `k`-bounded partition;
      if the tableaux are skew, the shape is a tuple of the outer and inner shape each
      represented as a `k`-bounded partition
    - ``weight`` -- the weight of the `k`-tableaux

    EXAMPLES::

        sage: T = WeakTableaux(3, [3,1], [2,2], representation = 'bounded')
        sage: T.list()
        [[[1, 1, 2], [2]]]

        sage: T = WeakTableaux(3, [[3,2,1], [2]], [1,1,1,1], representation = 'bounded')
        sage: T.list()
        [[[None, None, 3], [1, 4], [2]],
        [[None, None, 1], [2, 4], [3]],
        [[None, None, 1], [2, 3], [4]]]
    """
    @staticmethod
    def __classcall_private__(cls, k, shape, weight):
        r"""
        Straighten arguments before unique representation.

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_bounded
            sage: T = WeakTableaux_bounded(3, [2,1], [1,1,1])
            sage: TestSuite(T).run()
            sage: T = WeakTableaux_bounded(3, [[3,2,1], [2]], [1,1,1,1])
            sage: TestSuite(T).run()
        """
        if shape == [] or shape[0] in ZZ:
            shape = (Partition(shape), Partition([]))
        else:
            shape = tuple([Partition(r) for r in shape])
        return super(WeakTableaux_bounded, cls).__classcall__(cls, k, shape, tuple(weight))

    def __init__(self, k, shape, weight):
        r"""
        Initializes the parent class of (skew) weak `k`-tableaux in bounded representation.

        INPUT:

        - ``k`` -- positive integer
        - ``shape`` -- the shape of the `k`-tableaux represented as a `k`-bounded
          partition; if the tableaux are skew, the shape is a tuple of the outer and inner
          shape each represented as a `k`-bounded partition
        - ``weight`` -- the weight of the `k`-tableaux

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_bounded
            sage: T = WeakTableaux_bounded(3, [3,1], [2,2])
            sage: TestSuite(T).run()
            sage: T = WeakTableaux_bounded(3, [[3,2,1], [2]], [1,1,1,1])
            sage: TestSuite(T).run()
        """
        self.k = k
        self._skew = bool(shape[1])
        self._outer_shape = Partition(shape[0])
        self._inner_shape = Partition(shape[1])
        self._shape = (self._outer_shape, self._inner_shape)
        self._weight = tuple(weight)
        self._representation = 'bounded'
        Parent.__init__(self, category = FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_bounded
            sage: repr(WeakTableaux_bounded(3, [2,1], [1,1,1]))
            'Bounded weak 3-Tableaux of (skew) 3-bounded shape [2, 1] and weight (1, 1, 1)'
            sage: repr(WeakTableaux_bounded(3, [[3,2,1], [2]], [1,1,1,1]))
            'Bounded weak 3-Tableaux of (skew) 3-bounded shape ([3, 2, 1], [2]) and weight (1, 1, 1, 1)'
        """
        return "Bounded weak %s-Tableaux of (skew) %s-bounded shape %s and weight %s"%(self.k, self.k, self.shape(), self._weight)

    def __iter__(self):
        r"""
        TESTS::

            sage: T = WeakTableaux(3, [3,1], [2,2], representation = 'bounded')
            sage: T.list()
            [[[1, 1, 2], [2]]]
            sage: T = WeakTableaux(3, [3,2,2], [2,2,2,1], representation = 'bounded')
            sage: T.list()
            [[[1, 1, 4], [2, 2], [3, 3]], [[1, 1, 2], [2, 3], [3, 4]]]
            sage: T = WeakTableaux(3, [[3,2,2], [1]], [2,1,2,1], representation = 'bounded')
            sage: T.list()
            [[[None, 1, 4], [1, 2], [3, 3]],
            [[None, 1, 3], [1, 3], [2, 4]],
            [[None, 1, 1], [2, 3], [3, 4]]]
        """
        for t in SemistandardSkewTableaux([self._outer_shape, self._inner_shape], self._weight):
            if t.is_k_tableau(self.k):
                yield self(t)

    Element = WeakTableau_bounded

#Weak tableaux in terms of factorized permutations
class WeakTableau_factorized_permutation(WeakTableau_abstract):
    r"""
    A weak (skew) `k`-tableau represented in terms of factorizations of affine
    permutations into cyclically decreasing elements.
    """
    @staticmethod
    def straighten_input(t, k):
        r"""
        Straightens input.

        INPUT:

        - ``t`` -- a list of reduced words or a list of elements in the Weyl group of type
          `A_k^{(1)}`
        - ``k`` -- a positive integer

        EXAMPLES::

            sage: from sage.combinat.k_tableau import WeakTableau_factorized_permutation
            sage: WeakTableau_factorized_permutation.straighten_input([[2,0],[3,2],[1,0]], 3)
            (s2*s0, s3*s2, s1*s0)
            sage: W = WeylGroup(['A',4,1])
            sage: WeakTableau_factorized_permutation.straighten_input([W.an_element(),W.an_element()], 4)
            (s0*s1*s2*s3*s4, s0*s1*s2*s3*s4)

        TESTS::

            sage: WeakTableau_factorized_permutation.straighten_input([W.an_element(),W.an_element()], 3)
            Traceback (most recent call last):
            ...
            ValueError: inconsistent number of rows: should be 4 but got 5
        """
        W = WeylGroup(['A', k, 1], prefix='s')
        if len(t) > 0:
            if isinstance(t[0], (list, tuple)):
                w_tuple = tuple(W.from_reduced_word(p) for p in t)
            else:
                w_tuple = tuple(W(r) for r in t)
        else:
            w_tuple = tuple([W.one()])
        return w_tuple

    @staticmethod
    def __classcall_private__(cls, t, k, inner_shape = []):
        r"""
        Implements the shortcut ``WeakTableau_factorized_permutation(t, k)`` to
        ``WeakTableaux_factorized_permutation(k, shape, weight)(t)``
        where ``shape`` is the shape of the tableau as a `(k+1)`-core (or a tuple of
        `(k+1)`-cores if the tableau is skew) and ``weight`` is its weight.

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableau_factorized_permutation
            sage: t = WeakTableau_factorized_permutation([[2,0],[3,2],[1,0]], 3)
            sage: t.check()
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_factorized_permutation_with_category.element_class'>
            sage: TestSuite(t).run()

            sage: t = WeakTableau_factorized_permutation([[0,3],[2,1]], 3, inner_shape = [1,1])
            sage: t.check()
            sage: TestSuite(t).run()

            sage: t = WeakTableau_factorized_permutation([], 3); t
            [1]
            sage: t.check()
            sage: TestSuite(t).run()
        """
        if isinstance(t, cls):
            return t
        W = WeylGroup(['A', k, 1], prefix='s')
        w = cls.straighten_input(t, k)
        weight = tuple(w[i].length() for i in range(len(w) - 1, -1, -1))
        inner_shape = Core(inner_shape, k + 1)
        outer_shape = (W.prod(w)*W(inner_shape.to_grassmannian())).affine_grassmannian_to_core()
        return WeakTableaux_factorized_permutation(k, [outer_shape, inner_shape], weight)(w)

    def __init__(self, parent, t):
        r"""
        Initialization of (skew) weak `k`-tableau ``t`` in factorized permutation representation.

        INPUT:

        - ``t`` -- (skew) weak tableau in factorized permutation representation; the input
          can either be a list of reduced words of cyclically decreasing elements, or a
          list of cyclically decreasing elements;  when the tableau is skew, the inner
          shape needs to be specified as a `(k+1)`-core

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableau_factorized_permutation, WeakTableaux_factorized_permutation
            sage: c = WeakTableau_factorized_permutation([[2,0],[3,2],[1,0]], 3)
            sage: T = WeakTableaux_factorized_permutation(3, [5,2,1],[2,2,2])
            sage: t = T([[2,0],[3,2],[1,0]]); t
            [s2*s0, s3*s2, s1*s0]
            sage: c == t
            True
            sage: type(t)
            <class 'sage.combinat.k_tableau.WeakTableaux_factorized_permutation_with_category.element_class'>
            sage: t.parent()
            Factorized permutation (skew) weak 3-Tableaux of shape [5, 2, 1] and weight (2, 2, 2)
            sage: TestSuite(t).run()

            sage: t = WeakTableau_factorized_permutation([[2,0],[3,2]], 3, inner_shape = [2]); t
            [s2*s0, s3*s2]
            sage: t._inner_shape
            [2]
            sage: t.weight()
            (2, 2)
            sage: t.shape()
            ([5, 2, 1], [2])
            sage: TestSuite(t).run()

            sage: t = T([[3,0],[0,3],[1,0]])
            Traceback (most recent call last):
            ...
            ValueError: The outer shape of the parent does not agree with the outer shape of the tableau!

            sage: t = WeakTableau_factorized_permutation([], 3); t
            [1]
            sage: t.parent()._outer_shape
            []
            sage: t.parent()._weight
            (0,)
        """
        self.k = parent.k
        self._inner_shape = parent._inner_shape
        ClonableList.__init__(self, parent, self.straighten_input(t, parent.k))

    def shape_core(self):
        r"""
        Return the shape of ``self`` as a `(k+1)`-core.

        When the tableau is straight, the outer shape is returned as a core.
        When the tableau is skew, the tuple of the outer and inner shape is returned as
        cores.

        EXAMPLES::

            sage: t = WeakTableau([[2],[0,3],[2,1,0]], 3, representation = 'factorized_permutation')
            sage: t.shape_core()
            [5, 2, 1]

            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.shape()
            ([5, 2, 1], [2])
        """
        return self.shape()

    def shape_bounded(self):
        r"""
        Return the shape of ``self`` as a `k`-bounded partition.

        When the tableau is straight, the outer shape is returned as a `k`-bounded
        partition.  When the tableau is skew, the tuple of the outer and inner shape is
        returned as `k`-bounded partitions.

        EXAMPLES::

            sage: t = WeakTableau([[2],[0,3],[2,1,0]], 3, representation = 'factorized_permutation')
            sage: t.shape_bounded()
            [3, 2, 1]

            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.shape_bounded()
            ([3, 2, 1], [2])
        """
        if self.parent()._skew:
            return tuple([r.to_bounded_partition() for r in self.shape_core()])
        return self.shape_core().to_bounded_partition()

    def check(self):
        r"""
        Check that ``self`` is a valid weak `k`-tableau.

        EXAMPLES::

            sage: t = WeakTableau([[2],[0,3],[2,1,0]], 3, representation = 'factorized_permutation')
            sage: t.check()

        TESTS::

            sage: t = WeakTableau([[2,0],[3,2]], 3, representation = 'factorized_permutation')
            Traceback (most recent call last):
            ...
            ValueError: this only works on type 'A' affine Grassmannian elements

            sage: T = WeakTableaux(3, [4,1], [2,1], representation = 'factorized_permutation')
            sage: t = T([[2],[1],[0]])
            Traceback (most recent call last):
            ...
            ValueError: The weight of the parent does not agree with the weight of the tableau!
        """
        weight = tuple(self[i].length() for i in range(len(self) - 1, -1, -1))
        if not self.parent()._weight == weight:
            raise ValueError("The weight of the parent does not agree "
                             "with the weight of the tableau!")
        W = self[0].parent()
        outer = (W.prod(self)*W((self._inner_shape).to_grassmannian())).affine_grassmannian_to_core()
        if self.parent()._outer_shape != outer:
            raise ValueError("The outer shape of the parent does not agree with the outer shape of the tableau!")
        if not self._is_k_tableau():
            raise ValueError("This is not a proper weak %s-tableau"%(self.k))

    def _is_k_tableau(self):
        r"""
        Checks whether ``self`` is a valid weak `k`-tableau.

        EXAMPLES::

            sage: t = WeakTableau([[2],[0,3],[2,1,0]], 3, representation = 'factorized_permutation')
            sage: t._is_k_tableau()
            True

            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t._is_k_tableau()
            True
        """
        W = self[0].parent()
        if (W.prod(self)*W(self.parent()._inner_shape.to_grassmannian())).is_affine_grassmannian():
            return all( r.is_pieri_factor() for r in self )
        return False

    def to_core_tableau(self):
        r"""
        Return the weak `k`-tableau ``self`` where the shape of each restricted tableau is a `(k+1)`-core.

        EXAMPLES::

            sage: t = WeakTableau([[0], [3,1], [2,1], [0,4], [3,0], [4,2], [1,0]], 4, representation = 'factorized_permutation'); t
            [s0, s3*s1, s2*s1, s0*s4, s3*s0, s4*s2, s1*s0]
            sage: c = t.to_core_tableau(); c
            [[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]]
            sage: type(c)
            <class 'sage.combinat.k_tableau.WeakTableaux_core_with_category.element_class'>
            sage: t = WeakTableau([[]], 4, representation = 'factorized_permutation'); t
            [1]
            sage: t.to_core_tableau()
            []

            sage: from sage.combinat.k_tableau import WeakTableau_factorized_permutation
            sage: t = WeakTableau([[2,0],[3,2],[1,0]], 3, representation = 'factorized_permutation')
            sage: WeakTableau_factorized_permutation.from_core_tableau(t.to_core_tableau(), 3)
            [s2*s0, s3*s2, s1*s0]
            sage: t == WeakTableau_factorized_permutation.from_core_tableau(t.to_core_tableau(), 3)
            True

            sage: t = WeakTableau([[2,0],[3,2]], 3, inner_shape = [2], representation = 'factorized_permutation')
            sage: t.to_core_tableau()
            [[None, None, 1, 1, 2], [1, 2], [2]]
            sage: t == WeakTableau_factorized_permutation.from_core_tableau(t.to_core_tableau(), 3)
            True
        """
        W = self[0].parent()
        factor = W(self._inner_shape.to_grassmannian())
        shapes = [factor]
        for i in range(len(self)-1,-1,-1):
            factor = self[i]*factor
            shapes += [factor.affine_grassmannian_to_core()]
        if self.parent()._skew:
            l = [[None]*i for i in self._inner_shape]
        else:
            l = []
        for i in range(1,len(shapes)):
            p = shapes[i]
            if len(l) < len(p):
                l += [[]]
            l_new = []
            for j in range(len(l)):
                l_new += [l[j] + [i]*(p[j]-len(l[j]))]
            l = l_new
        return WeakTableau_core(l, self.k)

    @classmethod
    def from_core_tableau(cls, t, k):
        r"""
        Construct weak factorized affine permutation tableau from a `k`-core tableau.

        EXAMPLES::

            sage: from sage.combinat.k_tableau import WeakTableau_factorized_permutation
            sage: WeakTableau_factorized_permutation.from_core_tableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
            [s2*s0, s3*s2, s1*s0]
            sage: WeakTableau_factorized_permutation.from_core_tableau([[1, 1, 2, 3, 4, 4, 5, 5, 6], [2, 3, 5, 5, 6], [3, 4, 7], [5, 6], [6], [7]], 4)
            [s0, s3*s1, s2*s1, s0*s4, s3*s0, s4*s2, s1*s0]
            sage: WeakTableau_factorized_permutation.from_core_tableau([[None, 1, 1, 2, 2], [None, 2], [1]], 3)
            [s0*s3, s2*s1]
        """
        t = SkewTableau(list(t))
        shapes = [ Core(p, k+1).to_grassmannian() for p in intermediate_shapes(t) ] #t.to_chain() ]
        perms = [ shapes[i]*(shapes[i-1].inverse()) for i in range(len(shapes)-1,0,-1)]
        return cls(perms, k, inner_shape = t.inner_shape())

    def k_charge(self, algorithm = 'I'):
        r"""
        Return the `k`-charge of ``self``.

        OUTPUT:

        - a nonnegative integer

        EXAMPLES::

            sage: t = WeakTableau([[2,0],[3,2],[1,0]], 3, representation = 'factorized_permutation')
            sage: t.k_charge()
            2
            sage: t = WeakTableau([[0],[3],[2],[1],[3],[0]], 3, representation = 'factorized_permutation')
            sage: t.k_charge()
            8
            sage: t = WeakTableau([[0],[3,1],[2,1],[0,4],[3,0],[4,2],[1,0]], 4, representation = 'factorized_permutation')
            sage: t.k_charge()
            12
        """
        return self.to_core_tableau().k_charge(algorithm = algorithm)


class WeakTableaux_factorized_permutation(WeakTableaux_abstract):
    r"""
    The class of (skew) weak `k`-tableaux in the factorized permutation representation of shape ``shape`` (as `k+1`-core
    or tuple of `(k+1)`-cores in the skew case) and weight ``weight``.

    INPUT:

    - ``k`` -- positive integer
    - ``shape`` -- the shape of the `k`-tableaux represented as a `(k+1)`-core;
      in the skew case the shape is a tuple of the outer and inner shape both as `(k+1)`-cores
    - ``weight`` -- the weight of the `k`-tableaux

    EXAMPLES::

        sage: T = WeakTableaux(3, [4,1], [2,2], representation = 'factorized_permutation')
        sage: T.list()
        [[s3*s2, s1*s0]]

        sage: T = WeakTableaux(4, [[6,2,1], [2]], [2,1,1,1], representation = 'factorized_permutation')
        sage: T.list()
        [[s0, s4, s3, s4*s2], [s0, s3, s4, s3*s2], [s3, s0, s4, s3*s2]]
    """
    @staticmethod
    def __classcall_private__(cls, k, shape, weight):
        r"""
        Straighten arguments before unique representation.

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_factorized_permutation
            sage: T = WeakTableaux_factorized_permutation(3, [2,1], [1,1,1])
            sage: TestSuite(T).run()
            sage: T = WeakTableaux_factorized_permutation(4, [[6,2,1], [2]], [2,1,1,1])
            sage: TestSuite(T).run() # long time
        """
        if shape == [] or shape[0] in ZZ:
            shape = (Core(shape, k+1), Core([],k+1))
        else:
            shape = tuple([Core(r,k+1) for r in shape])
        return super(WeakTableaux_factorized_permutation, cls).__classcall__(cls, k, shape, tuple(weight))

    def __init__(self, k, shape, weight):
        r"""
        Initializes the parent class of weak `k`-tableaux in factorized permutation representation.

        INPUT:

        - ``k`` -- positive integer
        - ``shape`` -- the shape of the `k`-tableaux represented as a `(k+1)`-core;
          in the skew case the shape is a tuple of the outer and inner shape both as
          `(k+1)`-cores
        - ``weight`` -- the weight of the `k`-tableaux

        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_factorized_permutation
            sage: T = WeakTableaux_factorized_permutation(3, [4,1], [2,2])
            sage: TestSuite(T).run()
            sage: T = WeakTableaux_factorized_permutation(4, [[6,2,1], [2]], [2,1,1,1])
            sage: TestSuite(T).run() # long time
        """
        self.k = k
        self._skew = bool(shape[1])
        self._outer_shape = Core(shape[0], k+1)
        self._inner_shape = Core(shape[1], k+1)
        self._shape = (self._outer_shape, self._inner_shape)
        self._weight = weight
        self._representation = 'factorized_permutation'
        Parent.__init__(self, category = FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.k_tableau import WeakTableaux_factorized_permutation
            sage: repr(WeakTableaux_factorized_permutation(3, [2,1], [1,1,1]))
            'Factorized permutation (skew) weak 3-Tableaux of shape [2, 1] and weight (1, 1, 1)'
            sage: repr(WeakTableaux_factorized_permutation(4, [[6,2,1], [2]], [2,1,1,1]))
            'Factorized permutation (skew) weak 4-Tableaux of shape ([6, 2, 1], [2]) and weight (2, 1, 1, 1)'
        """
        return "Factorized permutation (skew) weak %s-Tableaux of shape %s and weight %s"%(self.k, self.shape(), self._weight)

    def __iter__(self):
        r"""
        TESTS::

            sage: T = WeakTableaux(3, [4,1], [2,2], representation = 'factorized_permutation')
            sage: T.list()
            [[s3*s2, s1*s0]]
            sage: T = WeakTableaux(3, [5,2,2], [2,2,2,1], representation = 'factorized_permutation')
            sage: T.list()
            [[s0, s3*s2, s0*s3, s1*s0], [s3, s2*s0, s3*s2, s1*s0]]
            sage: T = WeakTableaux(4, [[6,2,1], [2]], [2,1,1,1], representation = 'factorized_permutation')
            sage: T.list()
            [[s0, s4, s3, s4*s2], [s0, s3, s4, s3*s2], [s3, s0, s4, s3*s2]]
        """
        for t in WeakTableaux_core(self.k, self.shape(), self._weight):
            yield WeakTableau_factorized_permutation.from_core_tableau(t, self.k)

    Element = WeakTableau_factorized_permutation


######## END weak tableaux BEGIN strong tableaux

class StrongTableau(ClonableList, metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A (standard) strong `k`-tableau is a (saturated) chain in Bruhat order.

    Combinatorially, it is a sequence of embedded `k+1`-cores (subject to some conditions)
    together with a set of markings.

    A strong cover in terms of cores corresponds to certain translated ribbons. A marking
    corresponds to the choice of one of the translated ribbons, which is indicated by
    marking the head (southeast most cell in French notation) of the chosen ribbon.  For
    more information, see [LLMS2006]_ and [LLMSSZ2013]_.

    In Sage, a strong `k`-tableau is created by specifying `k`, a standard strong
    tableau together with its markings, and a weight `\mu`. Here the standard tableau is
    represented by a sequence of `k+1`-cores

    .. MATH::

        \lambda^{(0)} \subseteq \lambda^{(1)} \subseteq \cdots \subseteq \lambda^{(m)}

    where each of the `\lambda^{(i)}` is a `k+1`-core.  The standard tableau is a filling
    of the diagram for the core `\lambda^{(m)}/\lambda^{(0)}` where a strong cover
    is represented by letters `\pm i` in the skew shape `\lambda^{(i)}/\lambda^{(i-1)}`.
    Each skew `(k+1)`-core `\lambda^{(i)}/\lambda^{(i-1)}` is a ribbon or multiple
    copies of the same ribbon which are separated by `k+1` diagonals.  Precisely one of
    the copies of the ribbons will be marked in the largest diagonal of the connected
    component (the 'head' of the ribbon).  The marked cells are indicated by negative
    signs.

    The strong tableau is stored as a standard strong marked tableau (referred to as the
    standard part of the strong tableau) and a vector representing the weight.

    EXAMPLES::

        sage: StrongTableau( [[-1, -2, -3], [3]], 2, [3] )
        [[-1, -1, -1], [1]]
        sage: StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3, [2,2,3,1])
        [[-1, -1, -2, -3], [-2, 3, -3, 4], [2, 3], [-3, -4]]

    Alternatively, the strong `k`-tableau can also be entered directly in semistandard
    format and then the standard tableau and the weight are computed and stored::

        sage: T = StrongTableau([[-1,-1,-1],[1]], 2); T
        [[-1, -1, -1], [1]]
        sage: T.to_standard_list()
        [[-1, -2, -3], [3]]
        sage: T.weight()
        (3,)
        sage: T = StrongTableau([[-1, -1, -2, -3], [-2, 3, -3, 4], [2, 3], [-3, -4]], 3); T
        [[-1, -1, -2, -3], [-2, 3, -3, 4], [2, 3], [-3, -4]]
        sage: T.to_standard_list()
        [[-1, -2, -4, -7], [-3, 6, -6, 8], [4, 7], [-5, -8]]
        sage: T.weight()
        (2, 2, 3, 1)
    """

    def __init__(self, parent, T):
        """
        INPUT:

        - ``parent`` -- an instance of ``StrongTableaux``
        - ``T`` -- standard marked strong (possibly skew) `k`-tableau or a semistandard
          marked strong (possibly skew) `k`-tableau with inner cells represented by
          ``None``

        EXAMPLES::

            sage: T = StrongTableau( [[-1, -2, -3]], 3 ); T
            [[-1, -2, -3]]
            sage: T
            [[-1, -2, -3]]
            sage: T.weight()
            (1, 1, 1)
            sage: T.size()
            3
            sage: T.parent()
            Set of strong 3-tableaux of shape [3] and of weight (1, 1, 1)
            sage: StrongTableau( [[-1, -2, -3], [3]], 2 )
            [[-1, -2, -3], [3]]
            sage: StrongTableau( [[-1, -1, 2], [-2]], 2 )
            [[-1, -1, 2], [-2]]
            sage: T = StrongTableau( [[-1, -2, 3], [-3]], 2, weight=[2,1] ); T
            [[-1, -1, 2], [-2]]
            sage: T = StrongTableau( [[-1, -2, 3], [-3]], 2, weight=[0,2,1] ); T
            [[-2, -2, 3], [-3]]
            sage: T.weight()
            (0, 2, 1)
            sage: T.size()
            3
            sage: T.parent()
            Set of strong 2-tableaux of shape [3, 1] and of weight (0, 2, 1)
            sage: StrongTableau( [[-1, -2, 3], [-3]], 2, weight=[1,2] )
            Traceback (most recent call last):
            ...
            ValueError: The weight=(1, 2) and the markings on the standard tableau=[[-1, -2, 3], [-3]] do not agree.
            sage: StrongTableau( [[None, None, -2, -4], [None, None], [-1, -3], [2, 4], [-5], [5], [5], [5]], 4 )
            [[None, None, -2, -4], [None, None], [-1, -3], [2, 4], [-5], [5], [5], [5]]
            sage: StrongTableau( [[None, None, -2, -4], [None, None], [-1, -3], [2, 4], [-5], [5], [5], [5]], 4, weight=[2,2,1] )
            [[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]]
            sage: StrongTableau( [[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4)
            [[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]]

        TESTS::

            sage: T = StrongTableau([], 3); T
            []
            sage: T.weight()
            ()
            sage: T.parent()
            Set of strong 3-tableaux of shape [] and of weight ()
            sage: T = StrongTableau( [[None, None], [None, None]], 4, weight=() ); T
            [[None, None], [None, None]]
            sage: T.size()
            0
        """
        self.k = parent.k
        self._tableau = T
        ClonableList.__init__(self, parent, T)

    @staticmethod
    def __classcall_private__(cls, T, k, weight=None):
        r"""
        Straighten input and implement the shortcut ``StrongTableau(T, k, weight=None)``
        to ``StrongTableaux(k, shape, weight)(T)``.

        TESTS::

            sage: t = StrongTableau( [[-1, -2, -3]], 3 )
            sage: t.parent()
            Set of strong 3-tableaux of shape [3] and of weight (1, 1, 1)
            sage: TestSuite(t).run()
            sage: t = StrongTableau( [[-1, -2, 3], [-3]], 2, weight=[2,1] )
            sage: TestSuite(t).run()
            sage: StrongTableau([[-1,-1,-1]], 3)
            [[-1, -1, -1]]
            sage: StrongTableau([[None, None, None], [None]], 2)
            [[None, None, None], [None]]

            sage: StrongTableau([[-1, -2, -2], [1]], 2)
            Traceback (most recent call last):
            ...
            ValueError: Unable to parse strong marked tableau : [[-1, -2, -2], [1]]

            sage: StrongTableau([[-1,-1,-1,-1]], 3)
            Traceback (most recent call last):
            ...
            ValueError: [4] is not a 4-core

            sage: StrongTableau([[-1, -2], [2]], 3)
            Traceback (most recent call last):
            ...
            ValueError: The marks in [[-1, -2], [2]] are not correctly placed.

            sage: StrongTableau([[None, None, None], [None]], 3)
            Traceback (most recent call last):
            ...
            ValueError: [3, 1] is not a 4-core

            sage: StrongTableau([[None, -1, 2], [-2]], 2, [2])
            Traceback (most recent call last):
            ...
            ValueError: The weight=(2,) and the markings on the standard tableau=[[None, -1, 2], [-2]] do not agree.
        """
        if isinstance(T, cls):
            return T
        outer_shape = Core([len(t) for t in T], k + 1)
        loop = (row.count(None) for row in T)
        inner_shape = Core([x for x in loop if x], k + 1)
        Te = [v for row in T for v in row if v is not None] + [0]
        count_marks = tuple(Te.count(-(i+1)) for i in range(-min(Te)))
        if not all(v == 1 for v in count_marks):
            # if T is not standard -> turn into standard
            if weight is not None and tuple(weight) != count_marks:
                raise ValueError("Weight = %s and tableau = %s do not agree" % (weight, T))
            tijseq = StrongTableaux.marked_CST_to_transposition_sequence(T, k)
            if tijseq is None or len(tijseq)<sum(list(count_marks)):
                raise ValueError("Unable to parse strong marked tableau : %s"%T)
            T = StrongTableaux.transpositions_to_standard_strong( tijseq, k, [[None]*r for r in inner_shape] ) # build from scratch
            T = T.set_weight( count_marks )
            return T
        else:
            if weight is not None:
                count_marks = tuple(weight) # in the case that it is standard + weight
        return StrongTableaux.__classcall__(StrongTableaux, k, (outer_shape, inner_shape), count_marks)(T)

    def check(self):
        r"""
        Check that ``self`` is a valid strong `k`-tableau.

        This function verifies that the outer and inner shape of the parent class is equal to
        the outer and inner shape of the tableau, that the tableau portion of ``self`` is
        a valid standard tableau, that the marks are placed correctly and that the size
        and weight agree.

        EXAMPLES::

            sage: T = StrongTableau([[-1, -1, -2], [2]], 2)
            sage: T.check()
            sage: T = StrongTableau([[None, None, 2, -4, -4], [-1, 4], [-2]], 3)
            sage: T.check()

        TESTS::

            sage: ST = StrongTableaux(2, [3,1], [1,1,1,1])
            sage: ST([[-1,-2,3],[-3]])
            Traceback (most recent call last):
            ...
            ValueError: The size of the tableau [[-1, -2, 3], [-3]] and weight (1, 1, 1, 1) do not match
            sage: ST([[-1,-3],[-2],[3]])
            Traceback (most recent call last):
            ...
            ValueError: The outer shape of the parent does not agree with the outer shape of the tableau!

            sage: StrongTableau([[-1, -2, 2], [1]], 2)
            Traceback (most recent call last):
            ...
            ValueError: The marks in [[-1, -2, 2], [1]] are not correctly placed.

            sage: StrongTableau([[-1, -2, 3], [3]], 2)
            Traceback (most recent call last):
            ...
            ValueError: The marks in [[-1, -2, 3], [3]] are not correctly placed.

            sage: StrongTableau([[-1,-2,-4,7],[-3,6,-6,8],[4,-7],[-5,-8]], 3, [2,2,3,1])
            Traceback (most recent call last):
            ...
            ValueError: The weight=(2, 2, 3, 1) and the markings on the standard tableau=[[-1, -2, -4, 7], [-3, 6, -6, 8], [4, -7], [-5, -8]] do not agree.
        """
        T = SkewTableau(self.to_standard_list())
        outer = Core(T.outer_shape(),self.k+1)
        inner = Core(T.inner_shape(),self.k+1)
        if self.parent()._outer_shape != outer:
            raise ValueError("The outer shape of the parent does not agree with the outer shape of the tableau!")
        if self.parent()._inner_shape != inner:
            raise ValueError("The inner shape of the parent does not agree with the inner shape of the tableau!")
        if not self._is_valid_marked():
            raise ValueError("The marks in %s are not correctly placed."%(self.to_standard_list()))
        if not self._is_valid_standard():
            raise ValueError("At least one shape in %s is not a valid %s-core."%(self.to_standard_list(), self.k+1))
        if not self.outer_shape().length()-self.inner_shape().length()==self.size():
            raise ValueError("The size of the tableau %s and weight %s do not match"%(self.to_standard_list(),self.weight()))
        if not self.is_column_strict_with_weight( self.weight() ):
            raise ValueError("The weight=%s and the markings on the standard tableau=%s do not agree."%(self.weight(),self.to_standard_list()))

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        EXAMPLES::

            sage: t = StrongTableau([[-1, -1, -2], [2]], 2)
            sage: hash(t) == hash(t)
            True
        """
        return hash(tuple(tuple(x) for x in self)) + hash(self.parent().k)

    def _is_valid_marked( self ):
        r"""
        Check the validity of marks of a potential tableau ``self``.

        This method is called by method :meth:`check` and is not meant to be
        accessed by the user.

        This method first checks that there is one marked cell for the size of the
        tableau.  Then, for each marked cell, it verifies that the cell below and to the
        right is not a positive value.

        In other words this checks that the marked cells are at the head of the connected
        components.  This function verifies that the markings of ``self`` are
        consistent with a strong marked standard tableau.

        INPUT:

        - ``self`` -- a list of lists representing a potential *standard* marked tableau

        OUTPUT:

        - a boolean, ``True`` if the marks are properly placed in the tableau

        EXAMPLES::

            sage: all( T._is_valid_marked() for T in StrongTableaux.standard_marked_iterator(3, 6))
            True
            sage: StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3)._is_valid_marked()
            True
            sage: StrongTableau([[-1, -2, 3], [3]], 2)
            Traceback (most recent call last):
            ...
            ValueError: The marks in [[-1, -2, 3], [3]] are not correctly placed.

        Marking in the wrong place::

            sage: StrongTableau([[None, None, -4, 5, -5], [None, None], [-1, -3], [2], [-2], [2], [3]], 4)
            Traceback (most recent call last):
            ...
            ValueError: The marks in [[None, None, -4, 5, -5], [None, None], [-1, -3], [2], [-2], [2], [3]] are not correctly placed.

        No marking on a 2::

            sage: StrongTableau([[None, None, -4, 5, -5], [None, None], [-1, -3], [2], [2], [2], [3]], 4)
            Traceback (most recent call last):
            ...
            ValueError: Unable to parse strong marked tableau : [[None, None, -4, 5, -5], [None, None], [-1, -3], [2], [2], [2], [3]]

        TESTS::

            sage: StrongTableau([[None, None, None], [None]], 2)._is_valid_marked()
            True
            sage: StrongTableau([], 4)._is_valid_marked()
            True
        """
        T = self.to_standard_list()
        size = Core([len(t) for t in T], self.k+1).length()
        inner_size = Core([y for y in (len([x for x in row if x is None]) for row in T) if y > 0], self.k+1).length()
        if len(set(v for v in flatten(list(T)) if v in ZZ and v < 0)) != size - inner_size:
            return False # TT does not have exactly self.size() marked cells
        for i in range(len(T)):
            for j in range(len(T[i])):
                v = T[i][j]
                if v is not None and v<0 and ((i!=0 and T[i-1][j]==abs(v)) or (j<len(T[i])-1 and T[i][j+1]==abs(v))):
                    return False
        return True

    def _is_valid_standard( self ):
        r"""
        Test if ``self`` has a valid strong (un)marked standard part of the tableau.

        This method is called by method :meth:`check` and is not meant to be
        accessed by the user.

        This methods returns ``True`` if every intermediate shape (restricted to values
        less than or equal to `i` for each `i`) is a `k+1`-core and that the length
        of the `i+1`-restricted core is the length of the `i`-restricted core plus 1.

        OUTPUT:

        - a boolean, ``True`` means the standard strong marked tableau is valid

        EXAMPLES::

            sage: all( T._is_valid_standard() for T in StrongTableaux.standard_marked_iterator(4, 6))
            True

        Inner shape is not a 3-core::

            sage: StrongTableau([[None, None, None], [-1]], 2)
            Traceback (most recent call last):
            ...
            ValueError: [3] is not a 3-core

        Restrict to 1 and 2 is not a 5-core::

            sage: StrongTableau([[None, None, -4, 5, -5], [None, None], [-1, -3], [-2], [2], [3], [3]], 4)
            Traceback (most recent call last):
            ...
            ValueError: At least one shape in [[None, None, -4, 5, -5], [None, None], [-1, -3], [-2], [2], [3], [3]] is not a valid 5-core.

        TESTS::

            sage: StrongTableau([[None, None, None], [None]], 2)._is_valid_standard()
            True
            sage: StrongTableau([], 4)._is_valid_standard()
            True
        """
        Tshapes = intermediate_shapes(self.to_unmarked_standard_list())
        if not all(Partition(la).is_core(self.k + 1) for la in Tshapes):
            return False
        Tsizes = [Core(lam, self.k + 1).length() for lam in Tshapes]
        return all(Tsizes[i] == Tsizes[i+1]-1 for i in range(len(Tsizes)-1))

    def is_column_strict_with_weight( self, mu ):
        """
        Test if ``self`` is a column strict tableau with respect to the weight ``mu``.

        INPUT:

        - ``mu`` -- a vector of weights

        OUTPUT:

        - a boolean, ``True`` means the underlying column strict strong marked tableau is valid

        EXAMPLES::

            sage: StrongTableau([[-1, -2, -3], [3]], 2).is_column_strict_with_weight([3])
            True
            sage: StrongTableau([[-1, -2, 3], [-3]], 2).is_column_strict_with_weight([3])
            False

        TESTS::

            sage: StrongTableau([[None, None, None], [None]], 2).is_column_strict_with_weight([])
            True
            sage: StrongTableau([], 4).is_column_strict_with_weight([])
            True
        """
        ss = 0
        for i in range(len(mu)):
            for j in range(mu[i]-1):
                # the markings should move from left to right
                if self.content_of_marked_head( ss+j+1 ) >= self.content_of_marked_head( ss+j+2 ):
                    return False
            ss += mu[i]
        return True

    def _repr_diagram(self):
        r"""
        Return a string representing the pretty print of the tableau.

        EXAMPLES::

            sage: StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3, [2,2,3,1])._repr_diagram()
            ' -1 -1 -2 -3\n -2  3 -3  4\n  2  3\n -3 -4'
            sage: StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4)._repr_diagram()
            '  .  . -1 -2\n  .  .\n -1 -2\n  1  2\n -3\n  3\n  3\n  3'
            sage: StrongTableau([], 4)._repr_diagram()
            ''
        """
        return SkewTableau(self.to_list())._repr_diagram()

    def _repr_list(self):
        r"""
        Return a string representing the list of lists of the tableau.

        EXAMPLES::

            sage: StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3, [2,2,3,1])._repr_list()
            '[[-1, -1, -2, -3], [-2, 3, -3, 4], [2, 3], [-3, -4]]'
            sage: StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4)._repr_list()
            '[[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]]'
            sage: StrongTableau([], 4)._repr_list()
            '[]'
        """
        return repr(self.to_list())

    def _repr_compact(self):
        """
        Return a compact string representation of ``self``.

        EXAMPLES::

            sage: StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3, [2,2,3,1])._repr_compact()
            '-1,-1,-2,-3/-2,3,-3,4/2,3/-3,-4'
            sage: StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4)._repr_compact()
            '.,.,-1,-2/.,./-1,-2/1,2/-3/3/3/3'
            sage: StrongTableau([],4)._repr_compact()
            '-'
        """
        return SkewTableau(self.to_list())._repr_compact()

    def _repr_(self):
        r"""
        Return a representation of ``self``.

        To display a strong marked tableau we display the semistandard version.

        EXAMPLES::

            sage: StrongTableau( [[-1, -2, -3]], 3 )
            [[-1, -2, -3]]
            sage: StrongTableau( [[-1, -2, -3]], 3 , weight=[3])
            [[-1, -1, -1]]
            sage: StrongTableau( [], 3 )
            []
            sage: T = StrongTableau([[-1,-2,3],[-3]],2)
            sage: T
            [[-1, -2, 3], [-3]]
            sage: Tableaux.options(display="diagram")
            sage: T
             -1 -2  3
             -3
            sage: Tableaux.options(convention="French")
            sage: T
             -3
             -1 -2  3
            sage: Tableaux.options(display="compact")
            sage: T
            -1,-2,3/-3
            sage: Tableaux.options(display="list",convention="English")
        """
        return self.parent().options._dispatch(self, '_repr_', 'display')

    def cell_of_marked_head(self, v):
        r"""
        Return location of marked head labeled by ``v`` in the standard part of ``self``.

        Return the coordinates of the ``v``-th marked cell in the strong standard tableau
        ``self``.  If there is no mark, then the value returned is `(0, r)` where `r` is
        the length of the first row.

        INPUT:

        - ``v`` -- an integer representing the label in the standard tableau

        OUTPUT:

        - a pair of the coordinates of the marked cell with entry ``v``

        EXAMPLES::

            sage: T = StrongTableau([[-1, -3, 4, -5], [-2], [-4]], 3)
            sage: [ T.cell_of_marked_head(i) for i in range(1,7)]
            [(0, 0), (1, 0), (0, 1), (2, 0), (0, 3), (0, 4)]
            sage: T = StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4)
            sage: [ T.cell_of_marked_head(i) for i in range(1,7)]
            [(2, 0), (0, 2), (2, 1), (0, 3), (4, 0), (0, 4)]

        TESTS::

            sage: StrongTableau([],4).cell_of_marked_head(4)
            (0, 0)
        """
        T = self.to_standard_list()
        if T==[]:
            return (0,0)
        for i in range(len(T)):
            for j in range(len(T[i])):
                if T[i][j]==-v:
                    return (i,j)
        return (0,len(T[0]))

    def content_of_marked_head(self, v):
        r"""
        Return the diagonal of the marked label ``v`` in the standard part of ``self``.

        Return the content (the `j-i` coordinate of the cell) of the ``v``-th marked cell
        in the strong standard tableau ``self``.  If there is no mark, then the value
        returned is the size of first row.

        INPUT:

        - ``v`` -- an integer representing the label in the standard tableau

        OUTPUT:

        - an integer representing the residue of the location of the mark

        EXAMPLES::

            sage: [ StrongTableau([[-1, -3, 4, -5], [-2], [-4]], 3).content_of_marked_head(i) for i in range(1,7)]
            [0, -1, 1, -2, 3, 4]
            sage: T = StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4)
            sage: [ T.content_of_marked_head(i) for i in range(1,7)]
            [-2, 2, -1, 3, -4, 4]

        TESTS::

            sage: StrongTableau([],4).content_of_marked_head(4)
            0
        """
        c = self.cell_of_marked_head(v)
        return c[1]-c[0]

    def cells_of_marked_ribbon(self, v):
        r"""
        Return a list of all cells the marked ribbon labeled by ``v`` in the standard part of ``self``.

        Return the list of coordinates of the cells which are in the marked
        ribbon with label ``v`` in the standard part of the tableau.  Note that
        the result is independent of the weight of the tableau.

        The cells are listed from largest content (where the mark is located)
        to the smallest.  Hence, the first entry in this list will be the marked cell.

        INPUT:

        - ``v`` -- the entry of the standard tableau

        OUTPUT:

        - a list of pairs representing the coordinates of the cells of
          the marked ribbon

        EXAMPLES::

            sage: T = StrongTableau([[-1, -1, -2, -2, 3], [2, -3], [-3]],3)
            sage: T.to_standard_list()
            [[-1, -2, -3, -4, 6], [4, -6], [-5]]
            sage: T.cells_of_marked_ribbon(1)
            [(0, 0)]
            sage: T.cells_of_marked_ribbon(4)
            [(0, 3)]
            sage: T = StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3)
            sage: T.cells_of_marked_ribbon(6)
            [(1, 2), (1, 1)]
            sage: T.cells_of_marked_ribbon(9)
            []
            sage: T = StrongTableau([[None, None, -1, -1, 3], [1, -3], [-3]],3)
            sage: T.to_standard_list()
            [[None, None, -1, -2, 4], [2, -4], [-3]]
            sage: T.cells_of_marked_ribbon(1)
            [(0, 2)]

        TESTS::

            sage: StrongTableau([],3).cells_of_marked_ribbon(1)
            []
        """
        d = self.content_of_marked_head(v)
        T = SkewTableau(self.to_unmarked_standard_list())
        cells = []
        while d is not None:
            adt=[c for c in T.cells_by_content(d) if T[c[0]][c[1]]==v]
            if adt==[]:
                d = None
            else:
                d -= 1
                cells += adt
        return cells

    def cell_of_highest_head( self, v ):
        """
        Return the cell of the highest head of label ``v`` in the standard part of ``self``.

        Return the cell where the head of the ribbon in the highest row is located
        in the underlying standard tableau.  If there is no cell with entry ``v`` then
        the cell returned is `(0, r)` where `r` is the length of the first row.

        This cell is calculated by iterating through the diagonals of the tableau.

        INPUT:

        - ``v`` -- an integer indicating the label in the standard tableau

        OUTPUT:

        - a pair of integers indicating the coordinates of the head of the highest
          ribbon with label ``v``

        EXAMPLES::

            sage: T = StrongTableau([[-1,2,-3],[-2,3],[3]], 1)
            sage: [T.cell_of_highest_head(v) for v in range(1,5)]
            [(0, 0), (1, 0), (2, 0), (0, 3)]
            sage: T = StrongTableau([[None,None,-3,4],[3,-4]],2)
            sage: [T.cell_of_highest_head(v) for v in range(1,5)]
            [(1, 0), (1, 1), (0, 4), (0, 4)]

        TESTS::

            sage: StrongTableau([],2).cell_of_highest_head(1)
            (0, 0)
        """
        Tlist = SkewTableau(self.to_standard_list())
        if Tlist==[]:
            return (0, 0)
        r = len(Tlist[0])
        dout = (0, r)
        for d in range(-len(Tlist),r+1):
            for c in Tlist.cells_by_content(d):
                if nabs(Tlist[c[0]][c[1]])==v:
                    dout = c
            if dout!=(0, r) and dout[1]-dout[0]!=d:
                return dout
        return dout

    def content_of_highest_head( self, v ):
        r"""
        Return the diagonal of the highest head of the cells labeled ``v`` in the standard part of ``self``.

        Return the content of the cell of the head in the highest row of all ribbons labeled by ``v`` of
        the underlying standard tableau.  If there is no cell with entry ``v`` then
        the value returned is the length of the first row.

        INPUT:

        - ``v`` -- an integer representing the label in the standard tableau

        OUTPUT:

        - an integer representing the content of the head of the highest
          ribbon with label ``v``

        EXAMPLES::

            sage: [StrongTableau([[-1,2,-3],[-2,3],[3]], 1).content_of_highest_head(v) for v in range(1,5)]
            [0, -1, -2, 3]

        TESTS::

            sage: StrongTableau([], 4).content_of_highest_head(1)
            0
            sage: StrongTableau([[-1,-1]], 4).content_of_highest_head(3)
            2
        """
        c = self.cell_of_highest_head(v)
        return c[1]-c[0]

    def cells_head_dictionary(self):
        r"""
        Return a dictionary with the locations of the heads of all markings.

        Return a dictionary of values and lists of cells where the heads with the values
        are located.

        OUTPUT:

        - a dictionary with keys the entries in the tableau and values are the coordinates
          of the heads with those entries

        EXAMPLES::

            sage: T = StrongTableau([[-1,-2,-4,7],[-3,6,-6,8],[4,-7],[-5,-8]], 3)
            sage: T.cells_head_dictionary()
            {1: [(0, 0)],
             2: [(0, 1)],
             3: [(1, 0)],
             4: [(2, 0), (0, 2)],
             5: [(3, 0)],
             6: [(1, 2)],
             7: [(2, 1), (0, 3)],
             8: [(3, 1), (1, 3)]}
            sage: T = StrongTableau([[None, 4, -4, -6, -7, 8, 8, -8], [None, -5, 8, 8, 8], [-3, 6]],3)
            sage: T.cells_head_dictionary()
            {1: [(2, 0)],
             2: [(0, 2)],
             3: [(1, 1)],
             4: [(2, 1), (0, 3)],
             5: [(0, 4)],
             6: [(1, 4), (0, 7)]}
             sage: StrongTableau([[None, None], [None, -1]], 4).cells_head_dictionary()
             {1: [(1, 1)]}

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).cells_head_dictionary()
            {}
            sage: StrongTableau([],4).cells_head_dictionary()
            {}
        """
        return StrongTableaux.cells_head_dictionary(self.to_unmarked_standard_list())

    def cells_of_heads(self, v):
        r"""
        Return a list of cells of the heads with label ``v`` in the standard part of ``self``.

        A list of cells which are heads of the ribbons with label ``v`` in the
        standard part of the tableau ``self``.  If there is no cell labelled by ``v`` then return the empty
        list.

        INPUT:

        - ``v`` -- an integer label

        OUTPUT:

        - a list of pairs of integers of the coordinates of the heads of the ribbons
          with label ``v``

        EXAMPLES::

            sage: T = StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4)
            sage: T.cells_of_heads(1)
            [(2, 0)]
            sage: T.cells_of_heads(2)
            [(3, 0), (0, 2)]
            sage: T.cells_of_heads(3)
            [(2, 1)]
            sage: T.cells_of_heads(4)
            [(3, 1), (0, 3)]
            sage: T.cells_of_heads(5)
            [(4, 0)]
            sage: T.cells_of_heads(6)
            []

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).cells_of_heads(1)
            []
            sage: StrongTableau([],4).cells_of_heads(1)
            []
        """
        dout = self.cells_head_dictionary()
        if v in dout:
            return dout[v]
        else:
            return []

    def contents_of_heads(self, v):
        r"""
        A list of contents of the cells which are heads of the ribbons with label ``v``.

        If there is no cell labelled by ``v`` then return the empty list.

        INPUT:

        - ``v`` -- an integer label

        OUTPUT:

        - a list of integers of the content of the heads of the ribbons with label ``v``

        EXAMPLES::

            sage: T = StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4)
            sage: T.contents_of_heads(1)
            [-2]
            sage: T.contents_of_heads(2)
            [-3, 2]
            sage: T.contents_of_heads(3)
            [-1]
            sage: T.contents_of_heads(4)
            [-2, 3]
            sage: T.contents_of_heads(5)
            [-4]
            sage: T.contents_of_heads(6)
            []

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).contents_of_heads(1)
            []
            sage: StrongTableau([],4).contents_of_heads(1)
            []
        """
        return [c[1]-c[0] for c in self.cells_of_heads(v)]

    def entries_by_content(self, diag):
        r"""
        Return the entries on the diagonal of ``self``.

        Return the entries in the tableau that are in the cells `(i,j)` with
        `j-i` equal to ``diag`` (that is, with content equal to ``diag``).

        INPUT:

        - ``diag`` -- an integer indicating the diagonal

        OUTPUT:

        - a list (perhaps empty) of labels on the diagonal ``diag``

        EXAMPLES::

            sage: T = StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4)
            sage: T.entries_by_content(0)
            []
            sage: T.entries_by_content(1)
            []
            sage: T.entries_by_content(2)
            [-1]
            sage: T.entries_by_content(-2)
            [-1, 2]

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).entries_by_content(1)
            []
            sage: StrongTableau([],4).entries_by_content(1)
            []
        """
        return SkewTableau(self.to_list()).entries_by_content(diag)

    def entries_by_content_standard(self, diag):
        r"""
        Return the entries on the diagonal of the standard part of ``self``.

        Return the entries in the tableau that are in the cells `(i,j)` with
        `j-i` equal to ``diag`` (that is, with content equal to ``diag``) in the
        standard tableau.

        INPUT:

        - ``diag`` -- an integer indicating the diagonal

        OUTPUT:

        - a list (perhaps empty) of labels on the diagonal ``diag``

        EXAMPLES::

            sage: T = StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4)
            sage: T.entries_by_content_standard(0)
            []
            sage: T.entries_by_content_standard(1)
            []
            sage: T.entries_by_content_standard(2)
            [-2]
            sage: T.entries_by_content_standard(-2)
            [-1, 4]

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).entries_by_content_standard(1)
            []
            sage: StrongTableau([],4).entries_by_content_standard(1)
            []
        """
        return SkewTableau(self.to_standard_list()).entries_by_content(diag)

    def ribbons_above_marked(self, v):
        r"""
        Number of ribbons of label ``v`` higher than the marked ribbon in the standard part.

        Return the number of copies of the ribbon with label ``v``  in the standard part
        of ``self`` which are in a higher row than the marked ribbon.  Note that the result
        is independent of the weight of the tableau.

        INPUT:

        - ``v`` -- the entry of the standard tableau

        OUTPUT:

        - an integer representing the number of copies of the ribbon above the marked
          ribbon

        EXAMPLES::

            sage: T = StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3)
            sage: T.ribbons_above_marked(4)
            1
            sage: T.ribbons_above_marked(6)
            0
            sage: T.ribbons_above_marked(9)
            0
            sage: StrongTableau([[-1,-2,-3,-4],[2,3,4],[3,4],[4]], 1).ribbons_above_marked(4)
            3

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).ribbons_above_marked(1)
            0
            sage: StrongTableau([],4).ribbons_above_marked(1)
            0
        """
        d = self.content_of_marked_head(v)
        count = 0
        for i in range(self.k+1, len(self.to_standard_list())+d, self.k+1):
            count += int(v in self.entries_by_content_standard(d-i))
        return count

    def height_of_ribbon(self, v):
        r"""
        The number of rows occupied by one of the ribbons with label ``v``.

        The number of rows occupied by the marked ribbon with label ``v``
        (and by consequence the number of rows occupied by any ribbon with the same label)
        in the standard part of ``self``.

        INPUT:

        - ``v`` -- the label of the standard marked tableau

        OUTPUT:

        - a non-negative integer representing the number of rows
          occupied by the ribbon which is marked

        EXAMPLES::

            sage: T = StrongTableau([[-1, -1, -2, -2, 3], [2, -3], [-3]],3)
            sage: T.to_standard_list()
            [[-1, -2, -3, -4, 6], [4, -6], [-5]]
            sage: T.height_of_ribbon(1)
            1
            sage: T.height_of_ribbon(4)
            1
            sage: T = StrongTableau([[None,None,1,-2],[None,-3,4,-5],[-1,3],[-4,5]], 3)
            sage: T.height_of_ribbon(3)
            2
            sage: T.height_of_ribbon(6)
            0

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).height_of_ribbon(1)
            0
            sage: StrongTableau([],4).height_of_ribbon(1)
            0
        """
        return len(set(c[0] for c in self.cells_of_marked_ribbon(v)))

    def number_of_connected_components(self, v):
        r"""
        Number of connected components of ribbons with label ``v`` in the standard part.

        The number of connected components is calculated by finding the number of cells
        with label ``v`` in the standard part of the tableau and dividing by the number
        of cells in the ribbon.

        INPUT:

        - ``v`` -- the label of the standard marked tableau

        OUTPUT:

        - a non-negative integer representing the number of connected
          components

        EXAMPLES::

            sage: T = StrongTableau([[-1, -1, -2, -2, 3], [2, -3], [-3]],3)
            sage: T.to_standard_list()
            [[-1, -2, -3, -4, 6], [4, -6], [-5]]
            sage: T.number_of_connected_components(1)
            1
            sage: T.number_of_connected_components(4)
            2
            sage: T = StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3)
            sage: T.number_of_connected_components(6)
            1
            sage: T.number_of_connected_components(9)
            0

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).number_of_connected_components(1)
            0
            sage: StrongTableau([],4).number_of_connected_components(1)
            0
        """
        sz = len(self.cells_of_marked_ribbon(v))
        if sz==0:
            return 0
        T = self.to_standard_list()
        nocells = len([i for i in range(len(T)) for j in range(len(T[i])) if T[i][j]==v])+1
        return ZZ(nocells/sz)

    def intermediate_shapes(self):
        r"""
        Return the intermediate shapes of ``self``.

        A (skew) tableau with letters `1, 2, \ldots, \ell` can be viewed as a sequence of
        shapes, where the `i`-th shape is given by the shape of the subtableau on letters
        `1, 2, \ldots, i`.

        The output is the list of these shapes.  The marked cells are ignored so to
        recover the strong tableau one would need the intermediate shapes and the
        :meth:`content_of_marked_head` for each pair of adjacent shapes in the list.

        OUTPUT:

        - a list of lists of integers representing `k+1`-cores

        EXAMPLES::

            sage: T = StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3, [2,2,3,1])
            sage: T.intermediate_shapes()
            [[], [2], [3, 1, 1], [4, 3, 2, 1], [4, 4, 2, 2]]
            sage: T = StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4)
            sage: T.intermediate_shapes()
            [[2, 2], [3, 2, 1, 1], [4, 2, 2, 2], [4, 2, 2, 2, 1, 1, 1, 1]]

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).intermediate_shapes()
            [[2, 1]]
            sage: StrongTableau([],4).intermediate_shapes()
            [[]]
        """
        return intermediate_shapes(self.to_unmarked_list())

    def pp( self ):
        r"""
        Print the strong tableau ``self`` in pretty print format.

        EXAMPLES::

            sage: T = StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3, [2,2,3,1])
            sage: T.pp()
            -1 -1 -2 -3
            -2  3 -3  4
             2  3
            -3 -4
            sage: T = StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4)
            sage: T.pp()
              .  . -1 -2
              .  .
             -1 -2
              1  2
             -3
              3
              3
              3
            sage: Tableaux.options(convention="French")
            sage: T.pp()
              3
              3
              3
             -3
              1  2
             -1 -2
              .  .
              .  . -1 -2
            sage: Tableaux.options(convention="English")
        """
        print(self._repr_diagram())

    def outer_shape( self ):
        r"""
        Return the outer shape of ``self``.

        This method returns the outer shape of ``self`` as viewed as a ``Core``.
        The outer shape of a strong tableau is always a `(k+1)`-core.

        OUTPUT:

        - a `(k+1)`-core

        EXAMPLES::

            sage: StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4).outer_shape()
            [4, 2, 2, 2, 1, 1, 1, 1]
            sage: StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3, [2,2,3,1]).outer_shape()
            [4, 4, 2, 2]

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).outer_shape()
            [2, 1]
            sage: StrongTableau([],4).outer_shape()
            []
        """
        return self.parent().outer_shape()

    def inner_shape( self ):
        r"""
        Return the inner shape of ``self``.

        If ``self`` is a strong skew tableau, then this method returns the inner shape
        (the shape of the cells labelled with ``None``).
        If ``self`` is not skew, then the inner shape is empty.

        OUTPUT:

        - a `(k+1)`-core

        EXAMPLES::

            sage: StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4).inner_shape()
            [2, 2]
            sage: StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3, [2,2,3,1]).inner_shape()
            []

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).inner_shape()
            [2, 1]
            sage: StrongTableau([],4).inner_shape()
            []
        """
        return self.parent().inner_shape()

    def shape( self ):
        r"""
        Return the shape of ``self``.

        If ``self`` is a skew tableau then return a pair of `k+1`-cores consisting of the
        outer and the inner shape.  If ``self`` is strong tableau with no inner shape then
        return a `k+1`-core.

        INPUT:

        - ``form`` - optional argument to indicate 'inner', 'outer' or 'skew' (default : 'outer')

        OUTPUT:

        - a `k+1`-core or a pair of `k+1`-cores if form is not 'inner' or 'outer'

        EXAMPLES::

            sage: T = StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4)
            sage: T.shape()
            ([4, 2, 2, 2, 1, 1, 1, 1], [2, 2])
            sage: StrongTableau([[-1, -2, 3], [-3]], 2).shape()
            [3, 1]
            sage: type(StrongTableau([[-1, -2, 3], [-3]], 2).shape())
            <class 'sage.combinat.core.Cores_length_with_category.element_class'>

        TESTS::

            sage: StrongTableau([[None, None, None], [None]], 2).shape()
            ([3, 1], [3, 1])
            sage: StrongTableau([],4).shape()
            []
        """
        return self.parent().shape()

    def weight( self ):
        r"""
        Return the weight of the tableau.

        The weight is a list of non-negative integers indicating the number of 1s,
        number of 2s, number of 3s, etc.

        OUTPUT:

        - a list of non-negative integers

        EXAMPLES::

            sage: T = StrongTableau([[-1, -2, -3, 4], [-4], [-5]], 3); T.weight()
            (1, 1, 1, 1, 1)
            sage: T.set_weight([3,1,1]).weight()
            (3, 1, 1)
            sage: StrongTableau([[-1,-1,-2,-3],[-2,3,-3,4],[2,3],[-3,-4]], 3).weight()
            (2, 2, 3, 1)

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).weight()
            ()
            sage: StrongTableau([],4).weight()
            ()
        """
        return self.parent()._weight

    def size( self ):
        """
        Return the size of the strong tableau.

        The size of the strong tableau is the sum of the entries in the
        :meth:`weight`.  It will also be equal to the length of the
        outer shape (as a `k+1`-core) minus the length of the inner shape.

        .. SEEALSO:: :meth:`sage.combinat.core.Core.length`

        OUTPUT:

        - a non-negative integer

        EXAMPLES::

            sage: StrongTableau([[-1, -2, -3, 4], [-4], [-5]], 3).size()
            5
            sage: StrongTableau([[None, None, -1, 2], [-2], [-3]], 3).size()
            3

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).size()
            0
            sage: StrongTableau([],4).size()
            0
        """
        return sum(self.weight())

    def to_list( self ):
        """
        Return the marked column strict (possibly skew) tableau as a list of lists.

        OUTPUT:

        - a list of lists of integers or ``None``

        EXAMPLES::

            sage: StrongTableau([[-1, -2, -3, 4], [-4], [-5]], 3).set_weight([2,1,1,1]).to_list()
            [[-1, -1, -2, 3], [-3], [-4]]
            sage: StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4).to_list()
            [[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]]
            sage: StrongTableau([[-1, -2, -3, 4], [-4], [-5]], 3, [3,1,1]).to_list()
            [[-1, -1, -1, 2], [-2], [-3]]

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).to_list()
            [[None, None], [None]]
            sage: StrongTableau([],4).to_list()
            []
        """
        def f(v):
            # f is a function which maps v or -v to the weight value corresponding to the partition mu
            if v is None:
                return None
            else:
                return sgn(v)*min([i for i in range(len(self.weight())+1) if sum(self.weight()[:i])>=abs(v)])
        return [[f(v) for v in row] for row in self.to_standard_list()]

    def to_unmarked_list( self ):
        """
        Return the tableau as a list of lists with markings removed.

        Return the list of lists of the rows of the tableau where the markings have been
        removed.

        OUTPUT:

        - a list of lists of integers or ``None``

        EXAMPLES::

            sage: T = StrongTableau( [[-1, -2, -3, 4], [-4], [-5]], 3, [3,1,1])
            sage: T.to_unmarked_list()
            [[1, 1, 1, 2], [2], [3]]
            sage: TT = T.set_weight([2,1,1,1])
            sage: TT.to_unmarked_list()
            [[1, 1, 2, 3], [3], [4]]
            sage: StrongTableau( [[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4).to_unmarked_list()
            [[None, None, 1, 2], [None, None], [1, 2], [1, 2], [3], [3], [3], [3]]

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).to_unmarked_list()
            [[None, None], [None]]
            sage: StrongTableau([],4).to_unmarked_list()
            []
        """
        return [[nabs(v) for v in row] for row in self.to_list()]

    def to_standard_list(self):
        """
        Return the underlying standard strong tableau as a list of lists.

        Internally, for a strong tableau the standard strong tableau and its weight
        is stored separately.  This method returns the underlying standard part.

        OUTPUT:

        - a list of lists of integers or ``None``

        EXAMPLES::

            sage: StrongTableau([[-1, -2, -3, 4], [-4], [-5]], 3, [3,1,1]).to_standard_list()
            [[-1, -2, -3, 4], [-4], [-5]]
            sage: StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4).to_standard_list()
            [[None, None, -2, -4], [None, None], [-1, -3], [2, 4], [-5], [5], [5], [5]]

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).to_standard_list()
            [[None, None], [None]]
            sage: StrongTableau([],4).to_standard_list()
            []
        """
        return self._tableau

    def to_standard_tableau(self):
        """
        Return the underlying standard strong tableau as a ``StrongTableau`` object.

        Internally, for a strong tableau the standard strong tableau and its weight
        is stored separately.  This method returns the underlying standard part as a
        ``StrongTableau``.

        OUTPUT:

        - a strong tableau with standard weight

        EXAMPLES::

            sage: T = StrongTableau([[-1, -2, -3, 4], [-4], [-5]], 3, [3,1,1])
            sage: T.to_standard_tableau()
            [[-1, -2, -3, 4], [-4], [-5]]
            sage: T.to_standard_tableau() == T.to_standard_list()
            False
            sage: StrongTableau([[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4).to_standard_tableau()
            [[None, None, -2, -4], [None, None], [-1, -3], [2, 4], [-5], [5], [5], [5]]

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).to_standard_tableau()
            [[None, None], [None]]
            sage: StrongTableau([],4).to_standard_tableau()
            []
        """
        return StrongTableau(self._tableau, self.k)

    def to_unmarked_standard_list( self ):
        """
        Return the standard part of the tableau as a list of lists with markings removed.

        Return the list of lists of the rows of the tableau where the markings have been
        removed.

        OUTPUT:

        - a list of lists of integers or ``None``

        EXAMPLES::

            sage: StrongTableau( [[-1, -2, -3, 4], [-4], [-5]], 3, [3,1,1]).to_unmarked_standard_list()
            [[1, 2, 3, 4], [4], [5]]
            sage: StrongTableau( [[None, None, -1, -2], [None, None], [-1, -2], [1, 2], [-3], [3], [3], [3]], 4).to_unmarked_standard_list()
            [[None, None, 2, 4], [None, None], [1, 3], [2, 4], [5], [5], [5], [5]]

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).to_unmarked_standard_list()
            [[None, None], [None]]
            sage: StrongTableau([],4).to_unmarked_standard_list()
            []
        """
        return [[nabs(_) for _ in x] for x in self.to_standard_list()]

    def _latex_(self):
        r"""
        Return a latex method for the tableau.

        EXAMPLES::

            sage: T = StrongTableau( [[None, -1, -2, 3], [2, -3]], 2, weight=[2,1] )
            sage: Tableaux.options(convention = "English")
            sage: latex(T)
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\cline{1-4}
            \lr{}&\lr{1^\ast}&\lr{1^\ast}&\lr{2}\\\cline{1-4}
            \lr{1}&\lr{2^\ast}\\\cline{1-2}
            \end{array}$}
            }
            sage: Tableaux.options(convention = "French")
            sage: latex(T)
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[t]{*{4}c}\cline{1-2}
            \lr{1}&\lr{2^\ast}\\\cline{1-4}
            \lr{}&\lr{1^\ast}&\lr{1^\ast}&\lr{2}\\\cline{1-4}
            \end{array}$}
            }
        """
        def chi(x):
            if x is None:
                return ""
            if x in ZZ:
                s = "%s"%abs(x)
                if x<0:
                    s += "^\\ast"
                return s
            return "%s"%x
        T = [[chi(x) for x in row] for row in self.to_list()]
        from .output import tex_from_array
        return tex_from_array(T)

    def restrict( self, r ):
        r"""
        Restrict the standard part of the tableau to the labels `1, 2, \ldots, r`.

        Return the tableau consisting of the labels of the standard part of ``self``
        restricted to the labels of `1` through ``r``.  The result is another
        ``StrongTableau`` object.

        INPUT:

        - ``r`` -- an integer

        OUTPUT:

        - A strong tableau

        EXAMPLES::

            sage: T = StrongTableau([[None, None, -4, 5, -5], [None, None], [-1, -3], [-2], [2], [2], [3]], 4, weight=[1,1,1,1,1])
            sage: T.restrict(3)
            [[None, None], [None, None], [-1, -3], [-2], [2], [2], [3]]
            sage: TT = T.restrict(0)
            sage: TT
            [[None, None], [None, None]]
            sage: TT == StrongTableau( [[None, None], [None, None]], 4 )
            True
            sage: T.restrict(5) == T
            True

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).restrict(1)
            [[None, None], [None]]
            sage: StrongTableau([],4).restrict(1)
            []
        """
        rr = sum(self.weight()[:r])
        rest_tab = [y for y in ([x for x in row if x is None or abs(x)<=rr] for row in self.to_standard_list()) if y]
        new_parent = StrongTableaux( self.k, (Core([len(x) for x in rest_tab], self.k+1), self.inner_shape()), self.weight()[:r] )
        return new_parent(rest_tab)

    def set_weight( self, mu ):
        """
        Sets a new weight ``mu`` for ``self``.

        This method first tests if the underlying standard tableau is column-strict with
        respect to the weight ``mu``.  If it is, then it changes the weight and returns
        the tableau; otherwise it raises an error.

        INPUT:

        - ``mu`` -- a list of non-negative integers representing the new weight

        EXAMPLES::

            sage: StrongTableau( [[-1, -2, -3], [3]], 2 ).set_weight( [3] )
            [[-1, -1, -1], [1]]
            sage: StrongTableau( [[-1, -2, -3], [3]], 2 ).set_weight( [0,3] )
            [[-2, -2, -2], [2]]
            sage: StrongTableau( [[-1, -2, 3], [-3]], 2 ).set_weight( [2, 0, 1] )
            [[-1, -1, 3], [-3]]
            sage: StrongTableau( [[-1, -2, 3], [-3]], 2 ).set_weight( [3] )
            Traceback (most recent call last):
            ...
            ValueError: [[-1, -2, 3], [-3]] is not a semistandard strong tableau with respect to the partition [3]

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).set_weight([])
            [[None, None], [None]]
            sage: StrongTableau([],4).set_weight([])
            []
        """
        if sum(mu)!=self.size() or self.is_column_strict_with_weight( mu ):
            return StrongTableaux.__classcall__(StrongTableaux, self.k, (self.outer_shape(), self.inner_shape()), tuple(mu))(self.to_standard_list())
        else:
            raise ValueError("%s is not a semistandard strong tableau with respect to the partition %s"%(self,mu))

    def left_action( self, tij ):
        r"""
        Action of transposition ``tij`` on ``self`` by adding marked ribbons.

        Computes the left action of the transposition ``tij`` on the tableau.
        If ``tij`` acting on the element of the affine Grassmannian raises the length by 1,
        then this function will add a cell to the standard tableau.

        INPUT:

        - ``tij`` -- a transposition represented as a pair `(i, j)`.

        OUTPUT:

        - ``self`` after it has been modified by the action of the transposition ``tij``

        EXAMPLES::

            sage: StrongTableau( [[None, -1, -2, -3], [3], [-4]], 3, weight=[1,1,1,1] ).left_action([0,1])
            [[None, -1, -2, -3, 5], [3, -5], [-4]]
            sage: StrongTableau( [[None, -1, -2, -3], [3], [-4]], 3, weight=[1,1,1,1] ).left_action([4,5])
            [[None, -1, -2, -3, -5], [3, 5], [-4]]
            sage: T = StrongTableau( [[None, -1, -2, -3], [3], [-4]], 3, weight=[1,1,1,1] )
            sage: T.left_action([-3,-2])
            [[None, -1, -2, -3], [3], [-4], [-5]]
            sage: T = StrongTableau( [[None, -1, -2, -3], [3], [-4]], 3, weight=[3,1] )
            sage: T.left_action([-3,-2])
            [[None, -1, -1, -1], [1], [-2], [-3]]
            sage: T
            [[None, -1, -1, -1], [1], [-2]]
            sage: T.check()
            sage: T.weight()
            (3, 1)

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).left_action([-2,-1])
            [[None, None], [None], [-1]]
            sage: StrongTableau([],4).left_action([0,1])
            [[-1]]
        """
        T = StrongTableaux._left_action_list(copy.deepcopy( self.to_standard_list() ), tij, self.size()+1, self.k)
        return StrongTableau( T, self.k, self.weight()+(1,) )

    def follows_tableau( self ):
        r"""
        Return a list of strong marked tableaux with length one longer than ``self``.

        Return list of all strong tableaux obtained from ``self`` by extending to a core
        which follows the shape of ``self`` in the strong order.

        OUTPUT:

        - a list of strong tableaux which follow ``self`` in strong order

        EXAMPLES::

            sage: T = StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3, [2,2,3,1])
            sage: T.follows_tableau()
            [[[-1, -1, -2, -3, 5, 5, -5], [-2, 3, -3, 4], [2, 3], [-3, -4]],
             [[-1, -1, -2, -3, 5], [-2, 3, -3, 4], [2, 3, 5], [-3, -4], [-5]],
             [[-1, -1, -2, -3, 5], [-2, 3, -3, 4], [2, 3, -5], [-3, -4], [5]],
             [[-1, -1, -2, -3, -5], [-2, 3, -3, 4], [2, 3, 5], [-3, -4], [5]],
             [[-1, -1, -2, -3], [-2, 3, -3, 4], [2, 3], [-3, -4], [-5], [5], [5]]]
            sage: StrongTableau([[-1,-2],[-3,-4]],3).follows_tableau()
            [[[-1, -2, 5, 5, -5], [-3, -4]], [[-1, -2, 5], [-3, -4], [-5]],
             [[-1, -2, -5], [-3, -4], [5]], [[-1, -2], [-3, -4], [-5], [5], [5]]]

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).follows_tableau()
            [[[None, None, -1], [None]], [[None, None], [None, -1]], [[None, None], [None], [-1]]]
            sage: StrongTableau([],4).follows_tableau()
            [[[-1]]]
        """
        v = self.size()+1
        out = []
        for T in StrongTableaux.follows_tableau_unsigned_standard( self.to_standard_list(), self.k ):
            for m in StrongTableaux.cells_head_dictionary(T)[v]:
                TT = copy.deepcopy(T)
                TT[m[0]][m[1]] = -v
                out.append(StrongTableau(TT, self.k, self.weight()+(1,)))
        return out

    def spin_of_ribbon( self, v ):
        r"""
        Return the spin of the ribbon with label ``v`` in the standard part of ``self``.

        The spin of a ribbon is an integer statistic.  It is the sum of `(h-1) r` plus
        the number of connected components above the marked one where `h` is the height
        of the marked ribbon and `r` is the number of connected components.

        .. SEEALSO:: :meth:`height_of_ribbon`, :meth:`number_of_connected_components`,
          :meth:`ribbons_above_marked`

        INPUT:

        - ``v`` -- a label of the standard part of the tableau

        OUTPUT:

        - an integer value representing the spin of the ribbon with label ``v``.

        EXAMPLES::

            sage: T = StrongTableau([[-1,-2,5,6],[-3,-4,-7,8],[-5,-6],[7,-8]], 3)
            sage: [T.spin_of_ribbon(v) for v in range(1,9)]
            [0, 0, 0, 0, 0, 0, 1, 0]
            sage: T = StrongTableau([[None,None,-1,-3],[-2,3,-3,4],[2,3],[-3,-4]], 3)
            sage: [T.spin_of_ribbon(v) for v in range(1,7)]
            [0, 1, 0, 0, 1, 0]

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).spin_of_ribbon(1)
            0
            sage: StrongTableau([],4).spin_of_ribbon(1)
            0
        """
        return (self.height_of_ribbon(v)-1)*self.number_of_connected_components(v)+self.ribbons_above_marked(v)

    def spin( self ):
        r"""
        Return the spin statistic of the tableau ``self``.

        The spin is an integer statistic on a strong marked tableau.  It is
        the sum of `(h-1) r` plus the number of connected components above the
        marked one where `h` is the height of the marked ribbon and `r` is
        the number of connected components.

        .. SEEALSO:: :meth:`height_of_ribbon`, :meth:`number_of_connected_components`,
          :meth:`ribbons_above_marked`

        The `k`-Schur functions with a parameter `t` can be defined as

        .. MATH::

            s^{(k)}_\lambda[X; t] = \sum_T t^{spin(T)} m_{weight(T)}[X]

        where the sum is over all column strict marked strong `k`-tableaux
        of shape `\lambda` and partition content.

        OUTPUT:

        - an integer value representing the spin.

        EXAMPLES::

            sage: StrongTableau([[-1,-2,5,6],[-3,-4,-7,8],[-5,-6],[7,-8]], 3, [2,2,3,1]).spin()
            1
            sage: StrongTableau([[-1,-2,-4,-7],[-3,6,-6,8],[4,7],[-5,-8]], 3, [2,2,3,1]).spin()
            2
            sage: StrongTableau([[None,None,-1,-3],[-2,3,-3,4],[2,3],[-3,-4]], 3).spin()
            2
            sage: ks3 = SymmetricFunctions(QQ['t'].fraction_field()).kschur(3)
            sage: t = ks3.realization_of().t
            sage: m = ks3.ambient().realization_of().m()
            sage: myks221 = sum(sum(t**T.spin() for T in StrongTableaux(3,[3,2,1],weight=mu))*m(mu) for mu in Partitions(5, max_part=3))
            sage: myks221 == m(ks3[2,2,1])
            True
            sage: h = ks3.ambient().realization_of().h()
            sage: Core([4,4,2,2],4).to_bounded_partition()
            [2, 2, 2, 2]
            sage: ks3[2,2,2,2].lift().scalar(h[3,3,2]) == sum( t**T.spin() for T in StrongTableaux(3, [4,4,2,2], weight=[3,3,2]) )
            True

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).spin()
            0
            sage: StrongTableau([],4).spin()
            0
        """
        return sum(self.spin_of_ribbon(v) for v in range(1,self.size()+1))

    def to_transposition_sequence( self ):
        """
        Return a list of transpositions corresponding to ``self``.

        Given a strong column strict tableau ``self`` returns the list of transpositions
        which when applied to the left of an empty tableau gives the corresponding strong
        standard tableau.

        OUTPUT:

        - a list of pairs of values ``[i,j]`` representing the transpositions `t_{ij}`

        EXAMPLES::

            sage: T = StrongTableau([[-1, -1, -1], [1]],2)
            sage: T.to_transposition_sequence()
            [[2, 3], [1, 2], [0, 1]]
            sage: T = StrongTableau([[-1, -1, 2], [-2]],2)
            sage: T.to_transposition_sequence()
            [[-1, 0], [1, 2], [0, 1]]
            sage: T = StrongTableau([[None, -1, 2, -3], [-2, 3]],2)
            sage: T.to_transposition_sequence()
            [[3, 4], [-1, 0], [1, 2]]

        TESTS::

            sage: StrongTableau([[None, None], [None]], 4).to_transposition_sequence()
            []
            sage: StrongTableau([],4).to_transposition_sequence()
            []
        """
        return StrongTableaux.marked_CST_to_transposition_sequence( self.to_standard_list(), self.k )

class StrongTableaux(UniqueRepresentation, Parent):

    def __init__( self, k, shape, weight ):
        r"""
        TESTS::

            sage: strongT = StrongTableaux(2, [3,1], weight=[2,1])
            sage: TestSuite(strongT).run()

            sage: strongT = StrongTableaux(0, [2,2], weight=[2,2])
            Traceback (most recent call last):
            ...
            ValueError: The input k has to be a positive integer
        """
        self._outer_shape = shape[0]
        self._inner_shape = shape[1]
        self.k = k
        if weight is None:
            self._weight = (1,)*(self._outer_shape.length()-self._inner_shape.length())
        else:
            self._weight = weight
        Parent.__init__(self, category = FiniteEnumeratedSets())

    @staticmethod
    def __classcall_private__(cls, k, shape, weight=None):
        r"""
        Straighten arguments before unique representation.

        TESTS::

            sage: ST3 = StrongTableaux(3, [2,2], weight=[1,1,1,1])
            sage: TestSuite(ST3).run()
        """
        if k<=0:
            raise ValueError("The input k has to be a positive integer")
        if shape==[] or shape[0] in ZZ:
            outer_shape = Core(shape,k+1)
            inner_shape = Core([],k+1)
        else:
            outer_shape = Core(shape[0],k+1)
            inner_shape = Core(shape[1],k+1)
        if weight is not None:
            weight = tuple(weight)
        return super(StrongTableaux, cls).__classcall__(cls, k, (outer_shape, inner_shape), weight)

    def _repr_( self ):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: StrongTableaux(3, [2,2], weight=[1,1,1,1])
            Set of strong 3-tableaux of shape [2, 2] and of weight (1, 1, 1, 1)
            sage: StrongTableaux(3, [2,2])
            Set of strong 3-tableaux of shape [2, 2] and of weight (1, 1, 1, 1)
            sage: StrongTableaux(3, [[2,2],[1]], weight=[0,0,2,1])
            Set of strong 3-tableaux of shape [[2, 2], [1]] and of weight (0, 0, 2, 1)
            sage: StrongTableaux(3, [[],[]], weight=[])
            Set of strong 3-tableaux of shape [] and of weight ()
       """
        if self._inner_shape==Core([],self.k+1):
            s = "Set of strong %s-tableaux"%self.k
            s +=" of shape %s"%self._outer_shape
        else:
            s = "Set of strong %s-tableaux"%self.k
            s +=" of shape [%s, %s]"%(self._outer_shape, self._inner_shape)
        s +="%sand of weight %s"%(" ",self._weight)
        return s

    options = Tableaux.options

    def an_element(self):
        r"""
        Return the first generated element of the class of ``StrongTableaux``.

        EXAMPLES::

            sage: ST = StrongTableaux(3, [3], weight=[3])
            sage: ST.an_element()
            [[-1, -1, -1]]
        """
        return next(iter(self))

    def outer_shape(self):
        r"""
        Return the outer shape of the class of strong tableaux.

        OUTPUT:

        - a `k+1`-core

        EXAMPLES::

            sage: StrongTableaux( 2, [3,1] ).outer_shape()
            [3, 1]
            sage: type(StrongTableaux( 2, [3,1] ).outer_shape())
            <class 'sage.combinat.core.Cores_length_with_category.element_class'>
            sage: StrongTableaux( 4, [[2,1], [1]] ).outer_shape()
            [2, 1]
        """
        return self._outer_shape

    def inner_shape(self):
        r"""
        Return the inner shape of the class of strong tableaux.

        OUTPUT:

        - a `k+1`-core

        EXAMPLES::

            sage: StrongTableaux( 2, [3,1] ).inner_shape()
            []
            sage: type(StrongTableaux( 2, [3,1] ).inner_shape())
            <class 'sage.combinat.core.Cores_length_with_category.element_class'>
            sage: StrongTableaux( 4, [[2,1], [1]] ).inner_shape()
            [1]
        """
        return self._inner_shape

    def shape(self):
        r"""
        Return the shape of ``self``.

        If the ``self`` has an inner shape return a pair consisting of an inner and
        an outer shape.  If the inner shape is empty then return only the outer shape.

        OUTPUT:

        - a `k+1`-core or a pair of `k+1`-cores

        EXAMPLES::

            sage: StrongTableaux( 2, [3,1] ).shape()
            [3, 1]
            sage: type(StrongTableaux( 2, [3,1] ).shape())
            <class 'sage.combinat.core.Cores_length_with_category.element_class'>
            sage: StrongTableaux( 4, [[2,1], [1]] ).shape()
            ([2, 1], [1])
        """
        if self._inner_shape:
            return (self._outer_shape, self._inner_shape)
        return self._outer_shape

    def __iter__(self):
        r"""
        TESTS::

            sage: ST = StrongTableaux(3, [4,1], weight=[2,2])
            sage: ST.list()
            [[[-1, -1, -2, -2], [2]], [[-1, -1, 2, -2], [-2]]]
            sage: ST = StrongTableaux(3, [5,2,2], weight=[2,2,2,1])
            sage: ST.cardinality()
            14
            sage: StrongTableaux(3, [5,2,2], weight=[3,3,1]).list()
            [[[-1, -1, -1, -2, -2], [-2, 2], [2, -3]], [[-1, -1, -1, 2, -2], [-2, -2], [2, -3]], [[-1, -1, -1, -2, -3], [-2, -2], [2, 2]]]
            sage: StrongTableaux(3, [4,1,1]).cardinality()
            10
            sage: StrongTableaux(3, [5,2,2], weight=[6,1]).list() # there are no strong column strict tableaux of shape [5,2,2] and weight (6,1)
            []
            sage: StrongTableaux(3, [[5,2,2], [3,1,1]], weight=[2,1]).list()
            [[[None, None, None, -1, -1], [None, 1], [None, -2]],
             [[None, None, None, 1, -1], [None, -1], [None, -2]],
             [[None, None, None, -1, -2], [None, -1], [None, 1]]]
            sage: StrongTableaux(2, [[4,3,3,2,2,1,1], [2,1,1]], weight=[1,1,1,1]).cardinality()
            150
            sage: StrongTableaux(2, [[7,5,3,1], [2,1,1]], weight=[2,2]).cardinality()
            18
            sage: StrongTableaux(2, [[3,1],[3,1]]).list()
            [[[None, None, None], [None]]]
            sage: StrongTableaux(4, []).list()
            [[]]
        """
        size = sum(self._weight)
        if size==0:
            yield self([[None]*(row) for row in self._inner_shape])
        else:
            for unT in StrongTableaux.standard_unmarked_iterator( self.k, size, self._outer_shape, self._inner_shape ):
                for T in StrongTableaux.marked_given_unmarked_and_weight_iterator( unT, self.k, self._weight ):
                    yield T

    @classmethod
    def standard_unmarked_iterator( cls, k, size, outer_shape=None, inner_shape=[] ):
        r"""
        An iterator for standard unmarked strong tableaux.

        An iterator which generates all unmarked tableaux of a given ``size`` which are
        contained in ``outer_shape`` and which contain the ``inner_shape``.

        These are built recursively by building all standard marked strong tableaux of
        size ``size`` `-1` and adding all possible covers.

        If ``outer_shape`` is ``None`` then there is no restriction on the shape of the
        tableaux which are created.

        INPUT:

        - ``k``, ``size`` - a positive integers
        - ``outer_shape`` - a list representing a `k+1`-core (default: ``None``)
        - ``inner_shape`` - a list representing a `k+1`-core (default: [])

        OUTPUT:

        - an iterator which lists all standard strong unmarked tableaux with ``size``
          cells and which are contained in ``outer_shape`` and contain ``inner_shape``

        EXAMPLES::

            sage: list(StrongTableaux.standard_unmarked_iterator(2, 3))
            [[[1, 2, 3], [3]], [[1, 2], [3], [3]], [[1, 3, 3], [2]], [[1, 3], [2], [3]]]
            sage: list(StrongTableaux.standard_unmarked_iterator(2, 1, inner_shape=[1,1]))
            [[[None, 1, 1], [None]], [[None, 1], [None], [1]]]
            sage: len(list(StrongTableaux.standard_unmarked_iterator(4,4)))
            10
            sage: len(list(StrongTableaux.standard_unmarked_iterator(4,6)))
            98
            sage: len(list(StrongTableaux.standard_unmarked_iterator(4,4, inner_shape=[2,2])))
            92
            sage: len(list(StrongTableaux.standard_unmarked_iterator(4,4, outer_shape=[5,2,2,1], inner_shape=[2,2])))
            10

        TESTS::

            sage: list(StrongTableaux.standard_unmarked_iterator(2,0, outer_shape=[3,1], inner_shape=[3,1]))
            [[[None, None, None], [None]]]
            sage: list(StrongTableaux.standard_unmarked_iterator(4,0, outer_shape=[]))
            [[]]
        """
        if size==0:
            if outer_shape is None or Core(outer_shape,k+1).contains(inner_shape):
                yield [[None]*(inner_shape[i]) for i in range(len(inner_shape))]
        else:
            for T in cls.standard_unmarked_iterator(k, size-1, outer_shape, inner_shape):
                for TT in cls.follows_tableau_unsigned_standard(T, k):
                    if outer_shape is None or Core(outer_shape, k+1).contains([len(_) for _ in TT]):
                        yield TT

    @classmethod
    def marked_given_unmarked_and_weight_iterator(cls, unmarkedT, k, weight):
        r"""
        An iterator generating strong marked tableaux from an unmarked strong tableau.

        Iterator which lists all marked tableaux of weight ``weight`` such that the
        standard unmarked part of the tableau is equal to ``unmarkedT``.

        INPUT:

        - ``unmarkedT`` - a list of lists representing a strong unmarked tableau
        - ``k`` - a positive integer
        - ``weight`` - a list of non-negative integers indicating the weight

        OUTPUT:

        - an iterator that returns ``StrongTableau`` objects

        EXAMPLES::

            sage: ST = StrongTableaux.marked_given_unmarked_and_weight_iterator([[1,2,3],[3]], 2, [3])
            sage: list(ST)
            [[[-1, -1, -1], [1]]]
            sage: ST = StrongTableaux.marked_given_unmarked_and_weight_iterator([[1,2,3],[3]], 2, [0,3])
            sage: list(ST)
            [[[-2, -2, -2], [2]]]
            sage: ST = StrongTableaux.marked_given_unmarked_and_weight_iterator([[1,2,3],[3]], 2, [1,2])
            sage: list(ST)
            [[[-1, -2, -2], [2]]]
            sage: ST = StrongTableaux.marked_given_unmarked_and_weight_iterator([[1,2,3],[3]], 2, [2,1])
            sage: list(ST)
            [[[-1, -1, 2], [-2]], [[-1, -1, -2], [2]]]
            sage: ST = StrongTableaux.marked_given_unmarked_and_weight_iterator([[None, None, 1, 2, 4], [2, 4], [3]], 3, [3,1])
            sage: list(ST)
            []
            sage: ST = StrongTableaux.marked_given_unmarked_and_weight_iterator([[None, None, 1, 2, 4], [2, 4], [3]], 3, [2,2])
            sage: list(ST)
            [[[None, None, -1, -1, 2], [1, -2], [-2]],
             [[None, None, -1, -1, -2], [1, 2], [-2]]]

        TESTS::

            sage: list(StrongTableaux.marked_given_unmarked_and_weight_iterator([[None, None, None],[None]], 2, []))
            [[[None, None, None], [None]]]
            sage: list(StrongTableaux.marked_given_unmarked_and_weight_iterator([], 4, weight=[]))
            [[]]
        """
        td = StrongTableaux.cells_head_dictionary(unmarkedT)
        if td == {}:  # the tableau is empty
            yield StrongTableau(unmarkedT, k, [])
        else:
            import itertools
            dsc = Composition(weight).descents()
            for m in itertools.product(*[td[key] for key in sorted(td)]):
                if all(((m[i][1]-m[i][0]<m[i+1][1]-m[i+1][0]) or (i in dsc)) for i in range(len(m)-1)):
                   yield StrongTableaux.add_marking(unmarkedT, m, k, weight)

    @classmethod
    def add_marking( cls, unmarkedT, marking, k, weight ):
        r"""
        Add markings to a partially marked strong tableau.

        Given an partially marked standard tableau and a list of cells where the marks
        should be placed along with a ``weight``, return the semi-standard marked strong
        tableau.  The marking should complete the marking so that the result is a
        strong standard marked tableau.

        INPUT:

        - ``unmarkedT`` - a list of lists which is a partially marked strong `k`-tableau
        - ``marking`` - a list of pairs of coordinates where cells are to be marked
        - ``k`` - a positive integer
        - ``weight`` - a tuple of the weight of the output tableau

        OUTPUT:

        - a ``StrongTableau`` object

        EXAMPLES::

            sage: StrongTableaux.add_marking([[None,1,2],[2]], [(0,1), (1,0)], 2, [1,1])
            [[None, -1, 2], [-2]]
            sage: StrongTableaux.add_marking([[None,1,2],[2]], [(0,1), (1,0)], 2, [2])
            Traceback (most recent call last):
            ...
            ValueError: The weight=(2,) and the markings on the standard tableau=[[None, -1, 2], [-2]] do not agree.
            sage: StrongTableaux.add_marking([[None,1,2],[2]], [(0,1), (0,2)], 2, [2])
            [[None, -1, -1], [1]]

        TESTS::

            sage: StrongTableaux.add_marking([[None,None,None],[None]], [], 2, [])
            [[None, None, None], [None]]
            sage: StrongTableaux.add_marking([], [], 2, [])
            []
        """
        def msgn(c,v):
            if c in marking:
                return -v
            else:
                return v
        return StrongTableau([[msgn((i,j),unmarkedT[i][j]) for j in range(len(unmarkedT[i]))] for i in range(len(unmarkedT))], k, weight )

    @classmethod
    def _left_action_list( cls, Tlist, tij, v, k ):
        r"""
        Act by the transposition ``tij`` if it increases the size of the tableau by 1.

        This method modifies the tableau ``Tlist`` instead of returning a copy.

        INPUT:

        - ``Tlist`` - a partial standard strong `k`-tableau as a list of lists
        - ``tij`` - a pair of integers representing a transposition
        - ``v`` - the label to add to the tableau
        - ``k`` - a positive integer

        OUTPUT:

        - a list of lists, in particular, it is ``Tlist``

        EXAMPLES::

            sage: StrongTableaux._left_action_list( [[None]], [1,2], 10, 2 )
            [[None, -10]]
            sage: StrongTableaux._left_action_list( [[None]], [1,2], 10, 1 )
            [[None, -10], [10]]
            sage: StrongTableaux._left_action_list( [[None]], [2,3], 10, 1 )
            Traceback (most recent call last):
            ...
            ValueError: [2, 3] is not a single step up in the strong lattice
            sage: StrongTableaux._left_action_list( [[None]], [3,4], 10, 1 )
            [[None, 10], [10]]
            sage: T = StrongTableaux._left_action_list( [[None]], [1,2], 10, 2 )
            sage: StrongTableaux._left_action_list( T, [2,3], 4, 2 )
            [[None, -10, -4], [4]]
            sage: T
            [[None, -10, -4], [4]]
        """
        innershape = Core([len(_) for _ in Tlist], k+1)
        outershape = innershape.affine_symmetric_group_action(tij, transposition=True)
        if outershape.length()==innershape.length()+1:
            for c in SkewPartition([outershape.to_partition(),innershape.to_partition()]).cells():
                while c[0]>=len(Tlist):
                    Tlist.append([])
                Tlist[c[0]].append( v )
                if len(Tlist[c[0]])-c[0]==tij[1]:
                    Tlist[c[0]][-1] = -Tlist[c[0]][-1]  #mark the cell that is on the j-1 diagonal
            return Tlist
        else:
            raise ValueError("%s is not a single step up in the strong lattice"%tij)

    @classmethod
    def follows_tableau_unsigned_standard( cls, Tlist, k ):
        r"""
        Return a list of strong tableaux one longer in length than ``Tlist``.

        Return list of all standard strong tableaux obtained from ``Tlist`` by extending to
        a core which follows the shape of ``Tlist`` in the strong order.  It does not put
        the markings on the last entry that it adds but it does keep the markings on all
        entries smaller.  The objects returned are not ``StrongTableau`` objects (and
        cannot be) because the last entry will not properly marked.

        INPUT:

        - ``Tlist`` -- a filling of a `k+1`-core as a list of lists
        - ``k`` - an integer

        OUTPUT:

        - a list of strong tableaux which follow ``Tlist`` in strong order

        EXAMPLES::

            sage: StrongTableaux.follows_tableau_unsigned_standard([[-1, -1, -2, -3], [-2, 3, -3, 4], [2, 3], [-3, -4]], 3)
            [[[-1, -1, -2, -3, 5, 5, 5], [-2, 3, -3, 4], [2, 3], [-3, -4]],
             [[-1, -1, -2, -3, 5], [-2, 3, -3, 4], [2, 3, 5], [-3, -4], [5]],
             [[-1, -1, -2, -3], [-2, 3, -3, 4], [2, 3], [-3, -4], [5], [5], [5]]]
            sage: StrongTableaux.follows_tableau_unsigned_standard([[None,-1],[-2,-3]],3)
            [[[None, -1, 4, 4, 4], [-2, -3]], [[None, -1, 4], [-2, -3], [4]],
             [[None, -1], [-2, -3], [4], [4], [4]]]

        TESTS::

            sage: StrongTableaux.follows_tableau_unsigned_standard([[None, None, None], [None]], 2)
            [[[None, None, None, 1], [None, 1]], [[None, None, None], [None], [1]]]
            sage: StrongTableaux.follows_tableau_unsigned_standard([], 4)
            [[[1]]]
        """
        v = max([0]+[abs(v) for rows in Tlist for v in rows if v is not None])+1
        out = []
        sh = Core([len(_) for _ in Tlist], k+1)
        for ga in sh.strong_covers():
            T = copy.deepcopy(Tlist)
            T += [[] for i in range(len(ga)-len(T))]
            for c in SkewPartition([ga.to_partition(), sh.to_partition()]).cells():
                T[c[0]] += [v]
            out.append(T)
        return out

    @classmethod
    def standard_marked_iterator( cls, k, size, outer_shape=None, inner_shape=[] ):
        r"""
        An iterator for generating standard strong marked tableaux.

        An iterator which generates all standard marked `k`-tableaux of a given ``size``
        which are contained in ``outer_shape`` and contain the ``inner_shape``.
        If ``outer_shape`` is ``None`` then there is no restriction on the shape of the
        tableaux which are created.

        INPUT:

        - ``k`` - a positive integer
        - ``size`` - a positive integer
        - ``outer_shape`` - a list which is a `k+1`-core (default: ``None``)
        - ``inner_shape`` - a list which is a `k+1`-core (default: [])

        OUTPUT:

        - an iterator which returns the standard marked tableaux with ``size`` cells
          and that are contained in ``outer_shape`` and contain ``inner_shape``

        EXAMPLES::

            sage: list(StrongTableaux.standard_marked_iterator(2, 3))
            [[[-1, -2, 3], [-3]], [[-1, -2, -3], [3]], [[-1, -2], [-3], [3]], [[-1, 3, -3], [-2]], [[-1, 3], [-2], [-3]], [[-1, -3], [-2], [3]]]
            sage: list(StrongTableaux.standard_marked_iterator(2, 1, inner_shape=[1,1]))
            [[[None, 1, -1], [None]], [[None, 1], [None], [-1]], [[None, -1], [None], [1]]]
            sage: len(list(StrongTableaux.standard_marked_iterator(4,4)))
            10
            sage: len(list(StrongTableaux.standard_marked_iterator(4,6)))
            140
            sage: len(list(StrongTableaux.standard_marked_iterator(4,4, inner_shape=[2,2])))
            200
            sage: len(list(StrongTableaux.standard_marked_iterator(4,4, outer_shape=[5,2,2,1], inner_shape=[2,2])))
            24

        TESTS::

            sage: list(StrongTableaux.standard_marked_iterator(2,0,inner_shape=[3,1]))
            [[[None, None, None], [None]]]
            sage: list(StrongTableaux.standard_marked_iterator(4,0))
            [[]]
        """
        for T in cls.standard_unmarked_iterator( k, size, outer_shape, inner_shape ):
            for TT in cls.marked_given_unmarked_and_weight_iterator( T, k, [1]*(size) ):
                yield TT

    @classmethod
    def cells_head_dictionary( cls, T ):
        r"""
        Return a dictionary with the locations of the heads of all markings.

        Return a dictionary of values and lists of cells where the heads with the values
        are located in a strong standard unmarked tableau ``T``.

        INPUT:

        - ``T`` -- a strong standard unmarked tableau as a list of lists

        OUTPUT:

        - a dictionary with keys the entries in the tableau and values are the coordinates
          of the heads with those entries

        EXAMPLES::

            sage: StrongTableaux.cells_head_dictionary([[1,2,4,7],[3,6,6,8],[4,7],[5,8]])
            {1: [(0, 0)],
             2: [(0, 1)],
             3: [(1, 0)],
             4: [(2, 0), (0, 2)],
             5: [(3, 0)],
             6: [(1, 2)],
             7: [(2, 1), (0, 3)],
             8: [(3, 1), (1, 3)]}
            sage: StrongTableaux.cells_head_dictionary([[None, 2, 2, 4, 5, 6, 6, 6], [None, 3, 6, 6, 6], [1, 4]])
            {1: [(2, 0)],
             2: [(0, 2)],
             3: [(1, 1)],
             4: [(2, 1), (0, 3)],
             5: [(0, 4)],
             6: [(1, 4), (0, 7)]}

        TESTS::

             sage: StrongTableaux.cells_head_dictionary([[None, None, None],[None]])
             {}
             sage: StrongTableaux.cells_head_dictionary([])
             {}
        """
        if T == []:
            return {}
        ST = SkewTableau(T)
        dout = {}
        for i in range(-len(T), len(T[0])):
            nextv = ST.entries_by_content(i + 1)
            for c in ST.cells_by_content(i):
                v = T[c[0]][c[1]]
                if v not in nextv:
                    if v in dout:
                        dout[v] += [c]
                    else:
                        dout[v] = [c]
        return dout

    @classmethod
    def marked_CST_to_transposition_sequence(self, T, k):
        """
        Return a list of transpositions corresponding to ``T``.

        Given a strong column strict tableau ``T`` returns the list of transpositions
        which when applied to the left of an empty tableau gives the corresponding strong
        standard tableau.

        INPUT:

        - ``T`` -- a non-empty column strict tableau as a list of lists
        - ``k`` -- a positive integer

        OUTPUT:

        - a list of pairs of values ``[i,j]`` representing the transpositions `t_{ij}`

        EXAMPLES::

            sage: CST_to_trans = StrongTableaux.marked_CST_to_transposition_sequence
            sage: CST_to_trans([[-1, -1, -1], [1]], 2)
            [[2, 3], [1, 2], [0, 1]]
            sage: CST_to_trans([], 2)
            []
            sage: CST_to_trans([[-2, -2, -2], [2]], 2)
            [[2, 3], [1, 2], [0, 1]]
            sage: CST_to_trans([[-1, -2, -2, -2, -2], [-2, 2], [2]], 3)
            [[4, 5], [3, 4], [2, 3], [1, 2], [-1, 0], [0, 1]]
            sage: CST_to_trans([[-1, -2, -5, 5, -5, 5, -5], [-3, -4, 5, 5], [5]],3)
            [[5, 7], [3, 5], [2, 3], [0, 1], [-1, 0], [1, 2], [0, 1]]
            sage: CST_to_trans([[-1, -2, -3, 4, -7], [-4, -6], [-5, 6]],3)
            [[4, 5], [-1, 1], [-2, -1], [-1, 0], [2, 3], [1, 2], [0, 1]]

        TESTS::

            sage: StrongTableaux.marked_CST_to_transposition_sequence([[None, None, None], [None]], 2)
            []
            sage: StrongTableaux.marked_CST_to_transposition_sequence([], 4)
            []
        """
        LL = list(T)
        if not LL or all(v is None for v in sum(LL,[])):
            return []
        marks = [v for row in T for v in row if v is not None and v<0] + [0]
        m = -min(marks) # the largest marked cell
        transeq = [] # start with the empty list and append on the right
        sh = Core([len(_) for _ in T], k+1)
        j = max([ c-r for r,row in enumerate(LL) for c,val in enumerate(row)
                  if val == -m ])
        P = sh.to_partition()
        for l in range(k):
            msh = sh.affine_symmetric_group_action([j-l,j+1], transposition=True)
            mP = msh.to_partition()
            # my worry here is that the affine symmetric group action might apply an invalid
            # transposition but get something of the right length anyway.  How do I test if it is applying
            # a valid or invalid transposition?
            if msh.length() == sh.length() - 1:
                # if applying t_{j-l,j+1} reduces the size of the shape by 1
                valcells = [] # values in all the cells except content j
                regcells = [] # values in the cells with content j
                valid = True
                for (x,y) in SkewPartition([P, mP]).cells():
                    if y-x != j:
                        if LL[x][y] != m:
                            valid = False
                            break
                        valcells.append(LL[x][y])
                    else:
                        regcells.append(LL[x][y])
                if valid and regcells == [-m]:
                    # if all labels that are not content j are v and the label
                    # with content j = -m
                    mcells = mP.cells()
                    MM = [[LL[a][b] for b in range(len(LL[a])) if (a,b) in mcells]
                          for a in range(len(mP))]
                    transeq = self.marked_CST_to_transposition_sequence(MM, k)
                    if transeq is not None:
                        return [[j-l, j+1]] + transeq

    @classmethod
    def transpositions_to_standard_strong( self, transeq, k, emptyTableau=[] ):
        """
        Return a strong tableau corresponding to a sequence of transpositions.

        This method returns the action by left multiplication on the empty strong tableau
        by transpositions specified by ``transeq``.

        INPUT:

        - ``transeq`` -- a sequence of transpositions `t_{ij}` (a list of pairs).
        - ``emptyTableau`` -- (default: ``[]``) an empty list or a skew strong tableau
          possibly consisting of ``None`` entries

        OUTPUT:

        - a ``StrongTableau`` object

        EXAMPLES::

            sage: StrongTableaux.transpositions_to_standard_strong([[0,1]], 2)
            [[-1]]
            sage: StrongTableaux.transpositions_to_standard_strong([[-2,-1], [2,3]], 2, [[None, None]])
            [[None, None, -1], [1], [-2]]
            sage: StrongTableaux.transpositions_to_standard_strong([[2, 3], [1, 2], [0, 1]], 2)
            [[-1, -2, -3], [3]]
            sage: StrongTableaux.transpositions_to_standard_strong([[-1, 0], [1, 2], [0, 1]], 2)
            [[-1, -2, 3], [-3]]
            sage: StrongTableaux.transpositions_to_standard_strong([[3, 4], [-1, 0], [1, 2]], 2, [[None]])
            [[None, -1, 2, -3], [-2, 3]]

        TESTS::

            sage: StrongTableaux.transpositions_to_standard_strong([], 2, [[None, None, None], [None]])
            [[None, None, None], [None]]
            sage: StrongTableaux.transpositions_to_standard_strong([], 4, [])
            []
        """
        out = copy.deepcopy(emptyTableau)
        for i in range(1,len(transeq)+1):
            out = StrongTableaux._left_action_list(out, transeq[-i], i, k)
        return StrongTableau(out, k, weight = (1,)*len(transeq))

    Element = StrongTableau

#### common or global functions related to weak/strong tableaux

def nabs(v):
    r"""
    Return the absolute value of ``v`` or ``None``.

    INPUT:

    - ``v`` -- either an integer or ``None``

    OUTPUT:

    - either a non-negative integer or ``None``

    EXAMPLES::

        sage: from sage.combinat.k_tableau import nabs
        sage: nabs(None)
        sage: nabs(-3)
        3
        sage: nabs(None)
    """
    if v is None:
        return v
    else:
        return abs(v)

def intermediate_shapes(t):
    r"""
    Return the intermediate shapes of tableau ``t``.

    A (skew) tableau with letters `1, 2,\ldots, \ell` can be viewed as a sequence of
    shapes, where the `i`-th shape is given by the shape of the subtableau on letters
    `1, 2, \ldots, i`.  The output is the list of these shapes.

    OUTPUT:

    - a list of lists representing partitions

    EXAMPLES::

        sage: from sage.combinat.k_tableau import intermediate_shapes
        sage: t = WeakTableau([[1, 1, 2, 2, 3], [2, 3], [3]],3)
        sage: intermediate_shapes(t)
        [[], [2], [4, 1], [5, 2, 1]]

        sage: t = WeakTableau([[None, None, 2, 3, 4], [1, 4], [2]], 3)
        sage: intermediate_shapes(t)
        [[2], [2, 1], [3, 1, 1], [4, 1, 1], [5, 2, 1]]
    """
    shapes = []
    t = SkewTableau(list(t))
    for i in range(len(t.weight())+1):
        shapes += [ t.restrict(i).outer_shape()]
    return shapes
