r"""
Composition Tableaux

AUTHORS:

- Chris Berg, Jeff Ferreira (2012-9): Initial version
"""

from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.family import Family
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.composition import Composition, Compositions
from sage.combinat.partition import Partition
from sage.combinat.combinat import CombinatorialElement
from sage.rings.integer import Integer
from sage.combinat.backtrack import GenericBacktracker
import copy


class CompositionTableau(CombinatorialElement, metaclass=ClasscallMetaclass):
    r"""
    A composition tableau.

    A *composition tableau* `t` of shape `I = (I_1, \ldots, I_{\ell})` is an
    array of boxes in rows,  `I_i` boxes in row `i`, filled with positive
    integers such that:

    1) the entries in the rows of `t` weakly decrease from left to right,
    2) the left-most column of `t` strictly increase from top to bottom.
    3) Add zero entries to the rows of `t` until the resulting array is
       rectangular of  shape `\ell \times m`. For `1 \leq i < j \leq \ell,
       2 \leq k \leq m` and `(t(j,k) \neq 0`, and also if `t(j,k) \geq t(i,k))`
       implies `t(j,k) > t(i,k-1).`

    INPUT:

    - ``t`` -- A list of lists

    EXAMPLES::

        sage: CompositionTableau([[1],[2,2]])
        [[1], [2, 2]]
        sage: CompositionTableau([[1],[3,2],[4,4]])
        [[1], [3, 2], [4, 4]]
        sage: CompositionTableau([])
        []
    """
    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a composition tableau is only ever constructed as
        an ``element_class`` call of an appropriate parent.

        TESTS::

            sage: t = CompositionTableau([[1],[2,2]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Composition Tableaux
            sage: t.category()
            Category of elements of Composition Tableaux
        """
        if isinstance(t, CompositionTableau):
            return t
        return CompositionTableaux_all().element_class(CompositionTableaux_all(), t)

    def __init__(self, parent, t):
        r"""
        Initialize a composition tableau.

        TESTS::

            sage: t = CompositionTableaux()([[1],[2,2]])
            sage: s = CompositionTableaux(3)([[1],[2,2]])
            sage: s == t
            True
            sage: t.parent()
            Composition Tableaux
            sage: s.parent()
            Composition Tableaux of size 3 and maximum entry 3
            sage: r = CompositionTableaux()(s)
            sage: r.parent()
            Composition Tableaux
        """
        if isinstance(t, CompositionTableau):
            CombinatorialElement.__init__(self, parent, t._list)
            return

        # CombinatorialObject verifies that t is a list
        # We must verify t is a list of lists
        if not all(isinstance(row, list) for row in t):
            raise ValueError("A composition tableau must be a list of lists.")

        if not [len(_) for _ in t] in Compositions():
            raise ValueError("A composition tableau must be a list of non-empty lists.")

        # Verify rows weakly decrease from left to right
        for row in t:
            if any(row[i] < row[i+1] for i in range(len(row)-1)):
                raise ValueError("Rows must weakly decrease from left to right.")

        # Verify leftmost column strictly increases from top to bottom
        first_col = [row[0] for row in t if t!=[[]]]
        if any(first_col[i] >= first_col[i+1] for i in range(len(t)-1)):
            raise ValueError("Leftmost column must strictly increase from top to bottom.")

        # Verify triple condition
        l = len(t)
        m = max([len(_) for _ in t]+[0])
        TT = [row+[0]*(m-len(row)) for row in t]
        for i in range(l):
            for j in range(i+1,l):
                for k in range(1,m):
                    if TT[j][k] and TT[i][k] <= TT[j][k] <= TT[i][k-1]:
                        raise ValueError("Triple condition must be satisfied.")

        CombinatorialElement.__init__(self, parent, t)

    def _repr_diagram(self):
        r"""
        Return a string representation of ``self`` as an array.

        EXAMPLES::

            sage: t = CompositionTableau([[1],[3,2],[4,4]])
            sage: print(t._repr_diagram())
              1
              3  2
              4  4
        """
        return '\n'.join(("".join(("%3s" % str(x) for x in row))
                          for row in self))

    def __call__(self, *cell):
        r"""
        Return the value in the corresponding cell of ``self``.

        EXAMPLES::

            sage: t = CompositionTableau([[1],[3,2],[4,4]])
            sage: t(1,1)
            2
            sage: t(2,0)
            4
            sage: t(2,2)
            Traceback (most recent call last):
            ...
            IndexError: The cell (2,2) is not contained in [[1], [3, 2], [4, 4]]
        """
        try:
            i, j = cell
        except ValueError:
            i, j = cell[0]

        try:
            return self[i][j]
        except IndexError:
            raise IndexError("The cell (%d,%d) is not contained in %s"%(i,j,self))

    def pp(self):
        r"""
        Return a pretty print string of ``self``.

        EXAMPLES::

            sage: CompositionTableau([[1],[3,2],[4,4]]).pp()
            1
            3  2
            4  4
        """
        print(self._repr_diagram())

    def size(self):
        r"""
        Return the number of boxes in ``self``.

        EXAMPLES::

            sage: CompositionTableau([[1],[3,2],[4,4]]).size()
            5
        """
        return sum(len(row) for row in self)

    def weight(self):
        r"""
        Return a composition where entry `i` is the number of times that `i` appears in
        ``self``.

        EXAMPLES::

            sage: CompositionTableau([[1],[3,2],[4,4]]).weight()
            [1, 1, 1, 2, 0]
        """
        w = {i: 0 for i in range(1, self.size() + 1)}
        for row in self:
            for i in row:
                w[i] += 1
        return Composition([w[i] for i in range(1, self.size()+1)])

    def descent_set(self):
        r"""
        Return the set of all `i` that do *not* have `i+1` appearing strictly
        to the left of `i` in ``self``.

        EXAMPLES::

            sage: CompositionTableau([[1],[3,2],[4,4]]).descent_set()
            [1, 3]
        """
        cols = {}
        for row in self:
            for (col,i) in enumerate(row):
                cols[i] = col
        des_set = sorted([i for i in cols if i+1 in cols and cols[i+1] >= cols[i]])
        return des_set

    def descent_composition(self):
        r"""
        Return the composition corresponding to the set of all `i` that do
        not have `i+1` appearing strictly to the left of `i` in ``self``.

        EXAMPLES::

            sage: CompositionTableau([[1],[3,2],[4,4]]).descent_composition()
            [1, 2, 2]
        """
        return Composition(from_subset=(self.descent_set(), self.size()))

    def shape_composition(self):
        r"""
        Return a Composition object which is the shape of ``self``.

        EXAMPLES::

            sage: CompositionTableau([[1,1],[3,2],[4,4,3]]).shape_composition()
            [2, 2, 3]
            sage: CompositionTableau([[2,1],[3],[4]]).shape_composition()
            [2, 1, 1]
        """
        return Composition([len(row) for row in self])

    def shape_partition(self):
        r"""
        Return a partition which is the shape of ``self`` sorted into weakly
        decreasing order.

        EXAMPLES::

            sage: CompositionTableau([[1,1],[3,2],[4,4,3]]).shape_partition()
            [3, 2, 2]
            sage: CompositionTableau([[2,1],[3],[4]]).shape_partition()
            [2, 1, 1]
        """
        return Partition(sorted([len(row) for row in self], reverse=True))

    def is_standard(self):
        r"""
        Return ``True`` if ``self`` is a standard composition tableau and
        ``False`` otherwise.

        EXAMPLES::

            sage: CompositionTableau([[1,1],[3,2],[4,4,3]]).is_standard()
            False
            sage: CompositionTableau([[2,1],[3],[4]]).is_standard()
            True
        """
        entries = sum(self,[])
        return sorted(entries) == list(range(1, self.size() + 1))

class CompositionTableaux(UniqueRepresentation, Parent):
    r"""
    Composition tableaux.

    INPUT:

    Keyword arguments:

    - ``size`` -- the size of the composition tableaux
    - ``shape`` -- the shape of the composition tableaux
    - ``max_entry`` -- the maximum entry for the composition tableaux

    Positional arguments:

    - The first argument is interpreted as ``size`` or ``shape`` depending on
      whether it is an integer or a composition.

    EXAMPLES::

        sage: CT = CompositionTableaux(3); CT
        Composition Tableaux of size 3 and maximum entry 3
        sage: list(CT)
        [[[1], [2], [3]],
         [[1], [2, 2]],
         [[1], [3, 2]],
         [[1], [3, 3]],
         [[2], [3, 3]],
         [[1, 1], [2]],
         [[1, 1], [3]],
         [[2, 1], [3]],
         [[2, 2], [3]],
         [[1, 1, 1]],
         [[2, 1, 1]],
         [[2, 2, 1]],
         [[2, 2, 2]],
         [[3, 1, 1]],
         [[3, 2, 1]],
         [[3, 2, 2]],
         [[3, 3, 1]],
         [[3, 3, 2]],
         [[3, 3, 3]]]

        sage: CT = CompositionTableaux([1,2,1]); CT
        Composition tableaux of shape [1, 2, 1] and maximum entry 4
        sage: list(CT)
        [[[1], [2, 2], [3]],
         [[1], [2, 2], [4]],
         [[1], [3, 2], [4]],
         [[1], [3, 3], [4]],
         [[2], [3, 3], [4]]]

        sage: CT = CompositionTableaux(shape=[1,2,1],max_entry=3); CT
        Composition tableaux of shape [1, 2, 1] and maximum entry 3
        sage: list(CT)
        [[[1], [2, 2], [3]]]

        sage: CT = CompositionTableaux(2,max_entry=3); CT
        Composition Tableaux of size 2 and maximum entry 3
        sage: list(CT)
        [[[1], [2]],
         [[1], [3]],
         [[2], [3]],
         [[1, 1]],
         [[2, 1]],
         [[2, 2]],
         [[3, 1]],
         [[3, 2]],
         [[3, 3]]]

        sage: CT = CompositionTableaux(0); CT
        Composition Tableaux of size 0 and maximum entry 0
        sage: list(CT)
        [[]]
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`CompositionTableaux` for
        more information.

        TESTS::

            sage: CT = CompositionTableaux(3); CT
            Composition Tableaux of size 3 and maximum entry 3
            sage: CT = CompositionTableaux(size=3); CT
            Composition Tableaux of size 3 and maximum entry 3
            sage: CT = CompositionTableaux([1,2]); CT
            Composition tableaux of shape [1, 2] and maximum entry 3
            sage: CT = CompositionTableaux(shape=[1,2]); CT
            Composition tableaux of shape [1, 2] and maximum entry 3
            sage: CT = CompositionTableaux(shape=[]); CT
            Composition tableaux of shape [] and maximum entry 0
            sage: CT = CompositionTableaux(0); CT
            Composition Tableaux of size 0 and maximum entry 0
            sage: CT = CompositionTableaux(max_entry=3); CT
            Composition tableaux with maximum entry 3
            sage: CT = CompositionTableaux([1,2],max_entry=3); CT
            Composition tableaux of shape [1, 2] and maximum entry 3
            sage: CT = CompositionTableaux(size=2,shape=[1,2]); CT
            Traceback (most recent call last):
            ...
            ValueError: size and shape are different sizes
        """
        # Process keyword arguments first
        n = kwargs.get('n', None)
        size = kwargs.get('size', n)

        comp = kwargs.get('comp', None)
        shape = kwargs.get('shape', comp)

        max_entry = kwargs.get('max_entry', None)

        # Process positional arguments
        if args:
            # The first arg could be either a size or a shape
            if isinstance(args[0], (int, Integer)):
                if size is not None:
                    raise ValueError("size was specified more than once")
                else:
                    size = args[0]
            else:
                if shape is not None:
                    raise ValueError("the shape was specified more than once")
                shape = args[0]

        # Consistency checks
        if size is not None:
            if not isinstance(size, (int, Integer)):
                raise ValueError("size must be an integer")
            elif size < 0:
                raise ValueError("size must be non-negative")

        if shape is not None:
            # use in (and not isinstance) below so that lists can be used as
            # shorthand
            if shape not in Compositions():
                raise ValueError("shape must be a composition")
            if any(i == 0 for i in shape):
                raise ValueError("shape must have non-zero parts")
            shape = Composition(shape)

        if (size is not None) and (shape is not None):
            if sum(shape) != size:
                raise ValueError("size and shape are different sizes")

        if max_entry is not None:
            if not isinstance(max_entry, (int, Integer)):
                raise ValueError("max_entry must be an integer")
            elif max_entry <= 0:
                raise ValueError("max_entry must be positive")

        # Dispatch to appropriate class
        if (shape is not None):
            return CompositionTableaux_shape(shape, max_entry)

        if (size is not None):
            return CompositionTableaux_size(size, max_entry)

        return CompositionTableaux_all(max_entry)

    def __init__(self, **kwds):
        r"""
        Initialize ``self``.

        TESTS::

            sage: CT = CompositionTableaux()
            sage: TestSuite(CT).run()
        """
        if 'max_entry' in kwds:
            self.max_entry = kwds['max_entry']
            kwds.pop('max_entry')
        else:
            self.max_entry = None
        super(CompositionTableaux, self).__init__(**kwds)

    Element = CompositionTableau

    def _element_constructor_(self, t):
        r"""
        Construct an object from ``t`` as an element of ``self``, if
        possible.

        INPUT:

        - ``t`` -- data which can be interpreted as a composition tableau

        OUTPUT:

        - The corresponding CompositionTableau object

        TESTS::

            sage: CT = CompositionTableaux(3)
            sage: CT([[1],[2,2]]).parent() is CT
            True
            sage: CT([[1],[1,2]])
            Traceback (most recent call last):
            ...
            ValueError: [[1], [1, 2]] is not an element of Composition Tableaux of size 3 and maximum entry 3.
        """
        if t not in self:
            raise ValueError("%s is not an element of %s." % (t, self))

        return self.element_class(self, t)

    def __contains__(self, T):
        r"""
        Return ``True`` if ``T`` can be interpreted as
        :class:`CompositionTableau`.

        TESTS::

            sage: [[1],[2,2]] in CompositionTableaux(3)
            True
            sage: [[1],[2,2]] in CompositionTableaux(shape=[1,2])
            True
            sage: CompositionTableau([[1],[2,2]]) in CompositionTableaux()
            True
            sage: [[1],[2,2],[2]] in CompositionTableaux()
            False
        """
        if isinstance(T, CompositionTableau):
            return True

        # leftmost column of T strictly increases from top to bottom
        first_col = [row[0] for row in T]
        if any(first_col[i] >= first_col[i+1] for i in range(len(T)-1)):
            return False
        # rows of T weakly decrease from left to right
        for row in T:
            if any(row[i] < row[i+1] for i in range(len(row)-1)):
                return False
        # for 1 <= i < j <= len(comp), for 2 <= k <= m,
        #   T[j,k] \neq 0 and T[j,k] >= T[i,k] ==> T[j,k] > T[i,k-1]
        l = len(T)
        m = max([len(_) for _ in T]+[0])
        TT = [row+[0]*(m-len(row)) for row in T]
        for i in range(l):
            for j in range(i+1,l):
                for k in range(1,m):
                    if TT[j][k] != 0 and TT[j][k] >= TT[i][k] and TT[j][k] <= TT[i][k-1]:
                        return False
        return True

class CompositionTableaux_all(CompositionTableaux, DisjointUnionEnumeratedSets):
    r"""
    All composition tableaux.
    """
    def __init__(self, max_entry=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: CT = CompositionTableaux()
            sage: TestSuite(CT).run()
        """
        self.max_entry = max_entry
        CT_n = lambda n: CompositionTableaux_size(n, max_entry)
        DisjointUnionEnumeratedSets.__init__(self,
                    Family(NonNegativeIntegers(), CT_n),
                    facade=True, keepkey = False)

    def _repr_(self):
        r"""
        TESTS::

            sage: CompositionTableaux(3)
            Composition Tableaux of size 3 and maximum entry 3

            sage: CompositionTableaux()
            Composition Tableaux
        """
        if self.max_entry is not None:
            return "Composition tableaux with maximum entry %s"%str(self.max_entry)
        return "Composition Tableaux"

    def an_element(self):
        r"""
        Return a particular element of ``self``.

        EXAMPLES::

            sage: CT = CompositionTableaux()
            sage: CT.an_element()
            [[1, 1], [2]]
        """
        return self.element_class(self, [[1, 1], [2]])

class CompositionTableaux_size(CompositionTableaux):
    r"""
    Composition tableaux of a fixed size `n`.

    INPUT:

    - ``n`` -- a nonnegative integer.
    - ``max_entry`` -- a nonnegative integer. This keyword argument defaults to ``n``.

    OUTPUT:

    - The class of composition tableaux of size ``n``.
    """
    def __init__(self, n, max_entry=None):
        r"""
        Initializes the class of composition tableaux of size ``n``.

        TESTS::

            sage: CT = CompositionTableaux(4)
            sage: TestSuite(CT).run()
        """
        if max_entry is None:
            max_entry = n
        super(CompositionTableaux_size, self).__init__(max_entry=max_entry,
                category=FiniteEnumeratedSets())
        self.size = n

    def __contains__(self,x):
        r"""
        TESTS::

            sage: [[1],[2,2]] in CompositionTableaux(3)
            True
            sage: [[1],[2,2]] in CompositionTableaux(4)
            False
        """
        return CompositionTableaux.__contains__(self, x) and sum(map(len,x)) == self.size

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: [t for t in CompositionTableaux(3)]
            [[[1], [2], [3]],
             [[1], [2, 2]],
             [[1], [3, 2]],
             [[1], [3, 3]],
             [[2], [3, 3]],
             [[1, 1], [2]],
             [[1, 1], [3]],
             [[2, 1], [3]],
             [[2, 2], [3]],
             [[1, 1, 1]],
             [[2, 1, 1]],
             [[2, 2, 1]],
             [[2, 2, 2]],
             [[3, 1, 1]],
             [[3, 2, 1]],
             [[3, 2, 2]],
             [[3, 3, 1]],
             [[3, 3, 2]],
             [[3, 3, 3]]]

            sage: CompositionTableaux(3)[0].parent() is CompositionTableaux(3)
            True
        """
        for comp in Compositions(self.size):
            for T in CompositionTableaux_shape(comp,self.max_entry):
                yield self.element_class(self, T)

    def _repr_(self):
        r"""
        TESTS::

            sage: CompositionTableaux(3)
            Composition Tableaux of size 3 and maximum entry 3
        """
        return "Composition Tableaux of size %s and maximum entry %s"%(str(self.size), str(self.max_entry))

    def _an_element_(self):
        r"""
        Return a particular element of ``self``.

        EXAMPLES::

            sage: CT = CompositionTableaux(4)
            sage: CT.an_element()
            [[1, 1, 1], [2]]
            sage: CompositionTableaux(0).an_element()
            []
            sage: CompositionTableaux(1).an_element()
            [[1]]
        """
        if self.size == 0:
            return self.element_class(self, [])
        if self.size == 1:
            return self.element_class(self,[[1]])

        return self.element_class(self, [[1]*(self.size-1),[2]])

class CompositionTableaux_shape(CompositionTableaux):
    r"""
    Composition tableaux of a fixed shape ``comp`` with a given max entry.

    INPUT:

    - ``comp`` -- a composition.
    - ``max_entry`` -- a nonnegative integer. This keyword argument defaults
      to the size of ``comp``.
    """
    def  __init__(self, comp, max_entry=None):
        """
        Initialize ``sefl``.

        TESTS::

            sage: CT = CompositionTableaux([1,2])
            sage: TestSuite(CT).run()

            sage: CT = CompositionTableaux([1,2], max_entry=4)
            sage: TestSuite(CT).run()
        """
        if max_entry is None:
            max_entry = sum(comp)
        super(CompositionTableaux_shape, self).__init__(max_entry = max_entry,
              category = FiniteEnumeratedSets())
        self.shape = comp

    def __iter__(self):
        r"""
        An iterator for composition tableaux of a given shape.

        EXAMPLES::

            sage: [t for t in CompositionTableaux([1,2])]
            [[[1], [2, 2]], [[1], [3, 2]], [[1], [3, 3]], [[2], [3, 3]]]
            sage: [t for t in CompositionTableaux([1,2],max_entry=4)]
            [[[1], [2, 2]],
             [[1], [3, 2]],
             [[1], [3, 3]],
             [[1], [4, 2]],
             [[1], [4, 3]],
             [[1], [4, 4]],
             [[2], [3, 3]],
             [[2], [4, 3]],
             [[2], [4, 4]],
             [[3], [4, 4]]]
        """
        if sum(self.shape) == 0:
            yield CompositionTableau([])
        else:
            for z in CompositionTableauxBacktracker(self.shape, self.max_entry):
                yield CompositionTableau(z)

    def __contains__(self, x):
        r"""
        TESTS::

            sage: [[2],[4,3]] in CompositionTableaux([1,2])
            True
            sage: [[2],[3,2]] in CompositionTableaux([1,2])
            False
        """
        return CompositionTableaux.__contains__(self, x) and [len(_) for _ in x] == self.shape

    def _repr_(self):
        r"""
        TESTS::

            sage: CompositionTableaux([1,2,1])
            Composition tableaux of shape [1, 2, 1] and maximum entry 4
            sage: CompositionTableaux([1,2,1],max_entry=3)
            Composition tableaux of shape [1, 2, 1] and maximum entry 3
        """
        return "Composition tableaux of shape %s and maximum entry %s" % (str(self.shape), str(self.max_entry))

    def an_element(self):
        r"""
        Return a particular element of :class:`CompositionTableaux_shape`.

        EXAMPLES::

            sage: CT = CompositionTableaux([1,2,1])
            sage: CT.an_element()
            [[1], [2, 2], [3]]
        """
        if self.shape == []:
            return self.element_class(self, [])
        t = [[i]*len for (i,len) in enumerate(self.shape,start=1)]
        return self.element_class(self, t)


class CompositionTableauxBacktracker(GenericBacktracker):
    r"""
    A backtracker class for generating sets of composition tableaux.
    """
    def __init__(self, shape, max_entry=None):
        """
        EXAMPLES::

            sage: from sage.combinat.composition_tableau import CompositionTableauxBacktracker
            sage: n = CompositionTableauxBacktracker([1,3,2])
            sage: n._ending_position
            (2, 1)
            sage: n._initial_state
            (0, 0)
        """
        self._shape = shape
        self._n = sum(shape)
        self._initial_data = [ [None]*s for s in shape ]
        if max_entry is None:
            max_entry=sum(shape)
        self.max_entry=max_entry

        # The ending position will be at the lowest box which is farthest right
        ending_row = len(shape)-1
        ending_col = shape[-1]-1
        self._ending_position = (ending_row, ending_col)

        # Get the highest box that is farthest left
        starting_row = 0
        starting_col = 0

        GenericBacktracker.__init__(self, self._initial_data, (starting_row, starting_col))

    def _rec(self, obj, state):
        r"""
        EXAMPLES::

            sage: from sage.combinat.composition_tableau import CompositionTableauxBacktracker
            sage: n = CompositionTableauxBacktracker([1,3,2])
            sage: obj = [ [None], [None, None, None], [None, None] ]
            sage: list(n._rec(obj, n._initial_state))
            [([[1], [None, None, None], [None, None]], (1, 0), False),
             ([[2], [None, None, None], [None, None]], (1, 0), False),
             ([[3], [None, None, None], [None, None]], (1, 0), False),
             ([[4], [None, None, None], [None, None]], (1, 0), False),
             ([[5], [None, None, None], [None, None]], (1, 0), False),
             ([[6], [None, None, None], [None, None]], (1, 0), False)]
        """
        # Append zeros to a copy of obj
        obj_copy = copy.deepcopy(obj)
        N = max(len(u) for u in obj_copy)
        for a in range(len(obj_copy)):
            Na = len(obj_copy[a])
            obj_copy[a] += [0] * (N - Na)

        # We need to set the i,j^th entry.
        i, j = state

        # Get the next state
        new_state = self.get_next_pos(i, j)
        yld = True if new_state is None else False

        for k in range(1,self.max_entry +1):
            #We check to make sure that k does not violate the rule weak decrease in rows
            if j!=0 and obj[i][j-1] < k:
                continue

            #We check to make sure that k does not violate strict increase in first column
            if j == 0 and i != 0 and k <= obj[i-1][j]:
                continue

            #We check to make sure that k does not violate the Triple Rule
            if j != 0 and i != 0 and any(k == obj_copy[m][j] for m in range(i)):
                continue
            if j != 0 and i != 0 and any(obj_copy[m][j] < k and k <= obj_copy[m][j-1] for m in range(i)):
                continue

            #Fill in the in the i,j box with k
            obj[i][j] = k
            obj_copy[i][j] = k

            #Yield the object
            yield copy.deepcopy(obj), new_state, yld

    def get_next_pos(self, ii, jj):
        r"""
        EXAMPLES::

            sage: from sage.combinat.composition_tableau import CompositionTableauxBacktracker
            sage: T = CompositionTableau([[2,1],[5,4,3,2],[6],[7,7,6]])
            sage: n = CompositionTableauxBacktracker(T.shape_composition())
            sage: n.get_next_pos(1,1)
            (1, 2)
        """
        if (ii, jj) == self._ending_position:
            return None

        for j in range(jj+1, self._shape[ii]):
            if self._shape[ii] >= j:
                return ii, j

        return ii+1, 0

