# sage.doctest: optional - gap3
r"""
Finite real reflection groups
-------------------------------

Let `V` be a finite-dimensional real vector space. A reflection of
`V` is an operator `r \in \operatorname{GL}(V)` that has order `2`
and fixes pointwise a hyperplane in `V`.
In the present implementation, finite real reflection groups are
tied with a root system.

Finite real reflection groups with root systems have been classified
according to finite Cartan-Killing types.
For more definitions and classification types of finite complex
reflection groups, see :wikipedia:`Complex_reflection_group`.

The point of entry to work with reflection groups is :func:`~sage.combinat.root_system.reflection_group_real.ReflectionGroup`
which can be used with finite Cartan-Killing types::

    sage: ReflectionGroup(['A',2])
    Irreducible real reflection group of rank 2 and type A2
    sage: ReflectionGroup(['F',4])
    Irreducible real reflection group of rank 4 and type F4
    sage: ReflectionGroup(['H',3])
    Irreducible real reflection group of rank 3 and type H3

AUTHORS:

- Christian Stump (initial version 2011--2015)

.. WARNING::

    Uses the GAP3 package *Chevie* which is available as an
    experimental package (installed by ``sage -i gap3``) or to
    download by hand from `Jean Michel's website
    <http://webusers.imj-prg.fr/~jean.michel/gap3/>`_.
"""
# ****************************************************************************
#       Copyright (C) 2011-2016 Christian Stump <christian.stump at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_function, cached_method, cached_in_parent_method
from sage.combinat.root_system.cartan_type import CartanType, CartanType_abstract
from sage.interfaces.gap3 import gap3
from sage.rings.integer_ring import ZZ
from sage.combinat.root_system.reflection_group_complex import ComplexReflectionGroup, IrreducibleComplexReflectionGroup
from sage.combinat.root_system.reflection_group_element import RealReflectionGroupElement


def ReflectionGroup(*args, **kwds):
    r"""
    Construct a finite (complex or real) reflection group as a Sage
    permutation group by fetching the permutation representation of the
    generators from chevie's database.

    INPUT:

    Can be one or multiple of the following:

    - a triple `(r, p, n)` with `p` divides `r`, which denotes the group
      `G(r, p, n)`

    - an integer between `4` and `37`, which denotes an exceptional
      irreducible complex reflection group

    - a finite Cartan-Killing type

    EXAMPLES:

    Finite reflection groups can be constructed from

    Cartan-Killing classification types::

        sage: W = ReflectionGroup(['A',3]); W
         Irreducible real reflection group of rank 3 and type A3

        sage: W = ReflectionGroup(['H',4]); W
         Irreducible real reflection group of rank 4 and type H4

        sage: W = ReflectionGroup(['I',5]); W
         Irreducible real reflection group of rank 2 and type I2(5)

    the complex infinite family `G(r,p,n)` with `p` divides `r`::

        sage: W = ReflectionGroup((1,1,4)); W
        Irreducible real reflection group of rank 3 and type A3

        sage: W = ReflectionGroup((2,1,3)); W
        Irreducible real reflection group of rank 3 and type B3

    Chevalley-Shepard-Todd exceptional classification types::

        sage: W = ReflectionGroup(23); W
        Irreducible real reflection group of rank 3 and type H3

    Cartan types and matrices::

        sage: ReflectionGroup(CartanType(['A',2]))
        Irreducible real reflection group of rank 2 and type A2

        sage: ReflectionGroup(CartanType((['A',2],['A',2])))
        Reducible real reflection group of rank 4 and type A2 x A2

        sage: C = CartanMatrix(['A',2])
        sage: ReflectionGroup(C)
        Irreducible real reflection group of rank 2 and type A2

    multiples of the above::

        sage: W = ReflectionGroup(['A',2],['B',2]); W
        Reducible real reflection group of rank 4 and type A2 x B2

        sage: W = ReflectionGroup(['A',2],4); W
        Reducible complex reflection group of rank 4 and type A2 x ST4

        sage: W = ReflectionGroup((4,2,2),4); W
        Reducible complex reflection group of rank 4 and type G(4,2,2) x ST4
    """
    if not is_chevie_available():
        raise ImportError("the GAP3 package 'chevie' is needed to work with (complex) reflection groups")

    from sage.interfaces.gap3 import gap3
    gap3.load_package("chevie")

    error_msg = "the input data (%s) is not valid for reflection groups"

    W_types = []
    is_complex = False
    for arg in args:
        # preparsing
        if isinstance(arg, list):
            X = tuple(arg)
        else:
            X = arg

        # precheck for valid input data
        if not (isinstance(X, (CartanType_abstract,tuple)) or (X in ZZ and 4 <= X <= 37)):
            raise ValueError(error_msg % X)

        # transforming two reducible types and an irreducible type
        if isinstance(X, CartanType_abstract):
            if not X.is_finite():
                raise ValueError(error_msg % X)
            if hasattr(X,"cartan_type"):
                X = X.cartan_type()
            if X.is_irreducible():
                W_types.extend([(X.letter, X.n)])
            else:
                W_types.extend([(x.letter, x.n) for x in X.component_types()])

        elif X == (2,2,2) or X == ('I',2):
            W_types.extend([('A',1), ('A',1)])

        elif X == (2,2,3):
            W_types.extend([('A', 3)])

        else:
            W_types.append(X)

    # converting the real types given as complex types
    # and then checking for real vs complex
    for i,W_type in enumerate(W_types):
        if W_type in ZZ:
            if W_type == 23:
                W_types[i] = ('H', 3)
            elif W_type == 28:
                W_types[i] = ('F', 4)
            elif W_type == 30:
                W_types[i] = ('H', 4)
            elif W_type == 35:
                W_types[i] = ('E', 6)
            elif W_type == 36:
                W_types[i] = ('E', 7)
            elif W_type == 37:
                W_types[i] = ('E', 8)
        if isinstance(W_type,tuple) and len(W_type) == 3:
            if W_type[0] == W_type[1] == 1:
                W_types[i] = ('A', W_type[2]-1)
            elif W_type[0] == 2 and W_type[1] == 1:
                W_types[i] = ('B', W_type[2])
            elif W_type[0] == W_type[1] == 2:
                W_types[i] = ('D', W_type[2])
            elif W_type[0] == W_type[1] and W_type[2] == 2:
                W_types[i] = ('I', W_type[0])

        W_type = W_types[i]
        # check for real vs complex
        if W_type in ZZ or (isinstance(W_type, tuple) and len(W_type) == 3):
            is_complex = True

    for index_set_kwd in ['index_set', 'hyperplane_index_set', 'reflection_index_set']:
        index_set = kwds.get(index_set_kwd, None)
        if index_set is not None:
            if isinstance(index_set, (list, tuple)):
                kwds[index_set_kwd] = tuple(index_set)
            else:
                raise ValueError('the keyword %s must be a list or tuple' % index_set_kwd)

    if len(W_types) == 1:
        if is_complex is True:
            cls = IrreducibleComplexReflectionGroup
        else:
            cls = IrreducibleRealReflectionGroup
    else:
        if is_complex is True:
            cls = ComplexReflectionGroup
        else:
            cls = RealReflectionGroup
    return cls(tuple(W_types),
               index_set=kwds.get('index_set', None),
               hyperplane_index_set=kwds.get('hyperplane_index_set', None),
               reflection_index_set=kwds.get('reflection_index_set', None))


@cached_function
def is_chevie_available():
    r"""
    Test whether the GAP3 Chevie package is available.

    EXAMPLES::

        sage: # needs sage.groups
        sage: from sage.combinat.root_system.reflection_group_real import is_chevie_available
        sage: is_chevie_available() # random
        False
        sage: is_chevie_available() in [True, False]
        True
    """
    try:
        from sage.interfaces.gap3 import gap3
        gap3._start()
        gap3.load_package("chevie")
        return True
    except Exception:
        return False

#####################################################################
## Classes


class RealReflectionGroup(ComplexReflectionGroup):
    """
    A real reflection group given as a permutation group.

    .. SEEALSO::

        :func:`ReflectionGroup`
    """

    def __init__(self, W_types, index_set=None, hyperplane_index_set=None, reflection_index_set=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: W = ReflectionGroup(['A',3])
            sage: TestSuite(W).run()
        """
        W_types = tuple([tuple(W_type) if isinstance(W_type, (list,tuple)) else W_type
                         for W_type in W_types])
        cartan_types = []
        for W_type in W_types:
            W_type = CartanType(W_type)
            if not W_type.is_finite() or not W_type.is_irreducible():
                raise ValueError("the given Cartan type of a component is not irreducible and finite")
            cartan_types.append( W_type )
        if len(W_types) == 1:
            cls = IrreducibleComplexReflectionGroup
        else:
            cls = ComplexReflectionGroup
        cls.__init__(self, W_types, index_set=index_set,
                                    hyperplane_index_set=hyperplane_index_set,
                                    reflection_index_set=reflection_index_set)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3],['B',2],['I',5],['I',6])
            sage: W._repr_()
            'Reducible real reflection group of rank 9 and type A3 x B2 x I2(5) x G2'
        """
        type_str = ''
        for W_type in self._type:
            type_str += self._irrcomp_repr_(W_type)
            type_str += ' x '
        type_str = type_str[:-3]
        return 'Reducible real reflection group of rank %s and type %s' % (self._rank, type_str)

    def iteration(self, algorithm='breadth', tracking_words=True):
        r"""
        Return an iterator going through all elements in ``self``.

        INPUT:

        - ``algorithm`` -- (default: ``'breadth'``) must be one of
          the following:

          * ``'breadth'`` -- iterate over in a linear extension of the
            weak order
          * ``'depth'`` -- iterate by a depth-first-search
          * ``'parabolic'`` -- iterate by using parabolic subgroups

        - ``tracking_words`` -- boolean (default: ``True``); whether or not to
          keep track of the reduced words and store them in ``_reduced_word``

        .. NOTE::

            The fastest iteration is the parabolic iteration and the
            depth first algorithm without tracking words is second.
            In particular, ``'depth'`` is ~1.5x faster than ``'breadth'``.

        .. NOTE::

            The ``'parabolic'`` iteration does not track words and requires
            keeping the subgroup corresponding to `I \setminus \{i\}` in
            memory (for each `i` individually).

        EXAMPLES::

            sage: W = ReflectionGroup(["B",2])

            sage: for w in W.iteration("breadth",True):
            ....:     print("%s %s"%(w, w._reduced_word))
            () []
            (1,3)(2,6)(5,7) [1]
            (1,5)(2,4)(6,8) [0]
            (1,7,5,3)(2,4,6,8) [0, 1]
            (1,3,5,7)(2,8,6,4) [1, 0]
            (2,8)(3,7)(4,6) [1, 0, 1]
            (1,7)(3,5)(4,8) [0, 1, 0]
            (1,5)(2,6)(3,7)(4,8) [0, 1, 0, 1]

            sage: for w in W.iteration("depth", False): w
            ()
            (1,3)(2,6)(5,7)
            (1,5)(2,4)(6,8)
            (1,3,5,7)(2,8,6,4)
            (1,7)(3,5)(4,8)
            (1,7,5,3)(2,4,6,8)
            (2,8)(3,7)(4,6)
            (1,5)(2,6)(3,7)(4,8)
        """
        from sage.combinat.root_system.reflection_group_c import Iterator
        return iter(Iterator(self, N=self.number_of_reflections(),
                             algorithm=algorithm, tracking_words=tracking_words))

    def __iter__(self):
        r"""
        Return an iterator going through all elements in ``self``.

        For options and faster iteration see :meth:`iteration`.

        EXAMPLES::

            sage: W = ReflectionGroup(["B",2])

            sage: for w in W: print("%s %s"%(w, w._reduced_word))
            () []
            (1,3)(2,6)(5,7) [1]
            (1,5)(2,4)(6,8) [0]
            (1,7,5,3)(2,4,6,8) [0, 1]
            (1,3,5,7)(2,8,6,4) [1, 0]
            (2,8)(3,7)(4,6) [1, 0, 1]
            (1,7)(3,5)(4,8) [0, 1, 0]
            (1,5)(2,6)(3,7)(4,8) [0, 1, 0, 1]
        """
        return self.iteration(algorithm='breadth', tracking_words=True)

    @cached_method
    def bipartite_index_set(self):
        r"""
        Return the bipartite index set of a real reflection group.

        EXAMPLES::

            sage: W = ReflectionGroup(["A",5])
            sage: W.bipartite_index_set()
            [[1, 3, 5], [2, 4]]

            sage: W = ReflectionGroup(["A",5],index_set=['a','b','c','d','e'])
            sage: W.bipartite_index_set()
            [['a', 'c', 'e'], ['b', 'd']]
        """
        L, R = self._gap_group.BipartiteDecomposition().sage()
        inv = self._index_set_inverse
        L = [i for i in self._index_set if inv[i] + 1 in L]
        R = [i for i in self._index_set if inv[i] + 1 in R]
        return [L, R]

    def cartan_type(self):
        r"""
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3])
            sage: W.cartan_type()
            ['A', 3]

            sage: W = ReflectionGroup(['A',3], ['B',3])
            sage: W.cartan_type()
            A3xB3 relabelled by {1: 3, 2: 2, 3: 1}

        TESTS:

        Check that dihedral types are handled properly::

            sage: W = ReflectionGroup(['I',3]); W
            Irreducible real reflection group of rank 2 and type A2

            sage: W = ReflectionGroup(['I',4]); W
            Irreducible real reflection group of rank 2 and type C2

            sage: W = ReflectionGroup(['I',5]); W
            Irreducible real reflection group of rank 2 and type I2(5)
        """
        if len(self._type) == 1:
            ct = self._type[0]
            if ct['series'] == "I":
                C = CartanType([ct['series'], ct['bond']])
            else:
                C = CartanType([ct['series'], ct['rank']])
            CG = C.coxeter_diagram()
            G = self.coxeter_diagram()
            return C.relabel(CG.is_isomorphic(G, edge_labels=True, certificate=True)[1])
        else:
            return CartanType([W.cartan_type() for W in self.irreducible_components()])

    def positive_roots(self):
        r"""
        Return the positive roots of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3], ['B',2])
            sage: W.positive_roots()
            [(1, 0, 0, 0, 0),
             (0, 1, 0, 0, 0),
             (0, 0, 1, 0, 0),
             (0, 0, 0, 1, 0),
             (0, 0, 0, 0, 1),
             (1, 1, 0, 0, 0),
             (0, 1, 1, 0, 0),
             (0, 0, 0, 1, 1),
             (1, 1, 1, 0, 0),
             (0, 0, 0, 2, 1)]

            sage: W = ReflectionGroup(['A',3])
            sage: W.positive_roots()
            [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1), (1, 1, 1)]
        """
        return self.roots()[:self.number_of_reflections()]

    def almost_positive_roots(self):
        r"""
        Return the almost positive roots of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3], ['B',2])
            sage: W.almost_positive_roots()
            [(-1, 0, 0, 0, 0),
             (0, -1, 0, 0, 0),
             (0, 0, -1, 0, 0),
             (0, 0, 0, -1, 0),
             (0, 0, 0, 0, -1),
             (1, 0, 0, 0, 0),
             (0, 1, 0, 0, 0),
             (0, 0, 1, 0, 0),
             (0, 0, 0, 1, 0),
             (0, 0, 0, 0, 1),
             (1, 1, 0, 0, 0),
             (0, 1, 1, 0, 0),
             (0, 0, 0, 1, 1),
             (1, 1, 1, 0, 0),
             (0, 0, 0, 2, 1)]

            sage: W = ReflectionGroup(['A',3])
            sage: W.almost_positive_roots()
            [(-1, 0, 0),
             (0, -1, 0),
             (0, 0, -1),
             (1, 0, 0),
             (0, 1, 0),
             (0, 0, 1),
             (1, 1, 0),
             (0, 1, 1),
             (1, 1, 1)]
        """
        return [-beta for beta in self.simple_roots()] + self.positive_roots()

    def root_to_reflection(self, root):
        r"""
        Return the reflection along the given ``root``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])
            sage: for beta in W.roots(): W.root_to_reflection(beta)
            (1,4)(2,3)(5,6)
            (1,3)(2,5)(4,6)
            (1,5)(2,4)(3,6)
            (1,4)(2,3)(5,6)
            (1,3)(2,5)(4,6)
            (1,5)(2,4)(3,6)
        """
        Phi = self.roots()
        R = self.reflections()
        i = Phi.index(root) + 1
        j = Phi.index(-root) + 1
        for r in R:
            if r(i) == j:
                return r
        raise AssertionError("there is a bug in root_to_reflection")

    def reflection_to_positive_root(self, r):
        r"""
        Return the positive root orthogonal to the given reflection.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])
            sage: for r in W.reflections():
            ....:     print(W.reflection_to_positive_root(r))
            (1, 0)
            (0, 1)
            (1, 1)
        """
        Phi = self.roots()
        N = len(Phi) // 2
        for i in range(1, N + 1):
            if r(i) == i + N:
                return Phi[i - 1]
        raise AssertionError("there is a bug in reflection_to_positive_root")

    @cached_method
    def fundamental_weights(self):
        r"""
        Return the fundamental weights of ``self`` in terms of the simple roots.

        The fundamental weights are defined by
        `s_j(\omega_i) = \omega_i - \delta_{i=j}\alpha_j`
        for the simple reflection `s_j` with corresponding simple
        roots `\alpha_j`.

        In other words, the transpose Cartan matrix sends the weight
        basis to the root basis. Observe again that the action here is
        defined as a right action, see the example below.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3], ['B',2])
            sage: W.fundamental_weights()
            Finite family {1: (3/4, 1/2, 1/4, 0, 0), 2: (1/2, 1, 1/2, 0, 0), 3: (1/4, 1/2, 3/4, 0, 0), 4: (0, 0, 0, 1, 1/2), 5: (0, 0, 0, 1, 1)}

            sage: W = ReflectionGroup(['A',3])
            sage: W.fundamental_weights()
            Finite family {1: (3/4, 1/2, 1/4), 2: (1/2, 1, 1/2), 3: (1/4, 1/2, 3/4)}

            sage: W = ReflectionGroup(['A',3])
            sage: S = W.simple_reflections()
            sage: N = W.fundamental_weights()
            sage: for i in W.index_set():
            ....:     for j in W.index_set():
            ....:         print("{} {} {} {}".format(i, j, N[i], N[i]*S[j].to_matrix()))
            1 1 (3/4, 1/2, 1/4) (-1/4, 1/2, 1/4)
            1 2 (3/4, 1/2, 1/4) (3/4, 1/2, 1/4)
            1 3 (3/4, 1/2, 1/4) (3/4, 1/2, 1/4)
            2 1 (1/2, 1, 1/2) (1/2, 1, 1/2)
            2 2 (1/2, 1, 1/2) (1/2, 0, 1/2)
            2 3 (1/2, 1, 1/2) (1/2, 1, 1/2)
            3 1 (1/4, 1/2, 3/4) (1/4, 1/2, 3/4)
            3 2 (1/4, 1/2, 3/4) (1/4, 1/2, 3/4)
            3 3 (1/4, 1/2, 3/4) (1/4, 1/2, -1/4)
        """
        from sage.sets.family import Family
        m = self.cartan_matrix().transpose().inverse()
        Delta = tuple(self.simple_roots())
        zero = Delta[0].parent().zero()
        weights = [sum([m[i, j] * sj for j, sj in enumerate(Delta)], zero)
                   for i in range(len(Delta))]
        for weight in weights:
            weight.set_immutable()
        return Family({ind:weights[i] for i, ind in enumerate(self._index_set)})

    def fundamental_weight(self, i):
        r"""
        Return the fundamental weight with index ``i``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3])
            sage: [ W.fundamental_weight(i) for i in W.index_set() ]
            [(3/4, 1/2, 1/4), (1/2, 1, 1/2), (1/4, 1/2, 3/4)]
        """
        return self.fundamental_weights()[i]

    @cached_method
    def coxeter_diagram(self):
        """
        Return the Coxeter diagram associated to ``self``.

        EXAMPLES::

            sage: G = ReflectionGroup(['B',3])
            sage: G.coxeter_diagram().edges(labels=True, sort=True)
            [(1, 2, 4), (2, 3, 3)]
        """
        from sage.graphs.graph import Graph
        from itertools import combinations

        V = self.index_set()
        S = self.simple_reflections()
        E = []
        for i, j in combinations(V, 2):
            o = (S[i] * S[j]).order()
            if o >= 3:
                E.append((i, j, o))
        return Graph([V, E], format='vertices_and_edges', immutable=True)

    @cached_method
    def coxeter_matrix(self):
        """
        Return the Coxeter matrix associated to ``self``.

        EXAMPLES::

            sage: G = ReflectionGroup(['A',3])
            sage: G.coxeter_matrix()
            [1 3 2]
            [3 1 3]
            [2 3 1]
        """
        return self.cartan_type().coxeter_matrix()

    @cached_method
    def right_coset_representatives(self, J):
        r"""
        Return the right coset representatives of ``self`` for the
        parabolic subgroup generated by the simple reflections in ``J``.

        EXAMPLES::

            sage: W = ReflectionGroup(["A",3])
            sage: for J in Subsets([1,2,3]): W.right_coset_representatives(J)
            [(), (2,5)(3,9)(4,6)(8,11)(10,12), (1,4)(2,8)(3,5)(7,10)(9,11),
             (1,7)(2,4)(5,6)(8,10)(11,12), (1,2,10)(3,6,5)(4,7,8)(9,12,11),
             (1,4,6)(2,3,11)(5,8,9)(7,10,12), (1,6,4)(2,11,3)(5,9,8)(7,12,10),
             (1,7)(2,6)(3,9)(4,5)(8,12)(10,11),
             (1,10,2)(3,5,6)(4,8,7)(9,11,12), (1,2,3,12)(4,5,10,11)(6,7,8,9),
             (1,5,9,10)(2,12,8,6)(3,4,7,11), (1,6)(2,9)(3,8)(5,11)(7,12),
             (1,8)(2,7)(3,6)(4,10)(9,12), (1,10,9,5)(2,6,8,12)(3,11,7,4),
             (1,12,3,2)(4,11,10,5)(6,9,8,7), (1,3)(2,12)(4,10)(5,11)(6,8)(7,9),
             (1,5,12)(2,9,4)(3,10,8)(6,7,11), (1,8,11)(2,5,7)(3,12,4)(6,10,9),
             (1,11,8)(2,7,5)(3,4,12)(6,9,10), (1,12,5)(2,4,9)(3,8,10)(6,11,7),
             (1,3,7,9)(2,11,6,10)(4,8,5,12), (1,9,7,3)(2,10,6,11)(4,12,5,8),
             (1,11)(3,10)(4,9)(5,7)(6,12), (1,9)(2,8)(3,7)(4,11)(5,10)(6,12)]
            [(), (2,5)(3,9)(4,6)(8,11)(10,12), (1,4)(2,8)(3,5)(7,10)(9,11),
             (1,2,10)(3,6,5)(4,7,8)(9,12,11), (1,4,6)(2,3,11)(5,8,9)(7,10,12),
             (1,6,4)(2,11,3)(5,9,8)(7,12,10), (1,2,3,12)(4,5,10,11)(6,7,8,9),
             (1,5,9,10)(2,12,8,6)(3,4,7,11), (1,6)(2,9)(3,8)(5,11)(7,12),
             (1,3)(2,12)(4,10)(5,11)(6,8)(7,9),
             (1,5,12)(2,9,4)(3,10,8)(6,7,11), (1,3,7,9)(2,11,6,10)(4,8,5,12)]
            [(), (2,5)(3,9)(4,6)(8,11)(10,12), (1,7)(2,4)(5,6)(8,10)(11,12),
             (1,4,6)(2,3,11)(5,8,9)(7,10,12),
             (1,7)(2,6)(3,9)(4,5)(8,12)(10,11),
             (1,10,2)(3,5,6)(4,8,7)(9,11,12), (1,2,3,12)(4,5,10,11)(6,7,8,9),
             (1,10,9,5)(2,6,8,12)(3,11,7,4), (1,12,3,2)(4,11,10,5)(6,9,8,7),
             (1,8,11)(2,5,7)(3,12,4)(6,10,9), (1,12,5)(2,4,9)(3,8,10)(6,11,7),
             (1,11)(3,10)(4,9)(5,7)(6,12)]
            [(), (1,4)(2,8)(3,5)(7,10)(9,11), (1,7)(2,4)(5,6)(8,10)(11,12),
             (1,2,10)(3,6,5)(4,7,8)(9,12,11), (1,6,4)(2,11,3)(5,9,8)(7,12,10),
             (1,10,2)(3,5,6)(4,8,7)(9,11,12), (1,5,9,10)(2,12,8,6)(3,4,7,11),
             (1,8)(2,7)(3,6)(4,10)(9,12), (1,12,3,2)(4,11,10,5)(6,9,8,7),
             (1,3)(2,12)(4,10)(5,11)(6,8)(7,9),
             (1,11,8)(2,7,5)(3,4,12)(6,9,10), (1,9,7,3)(2,10,6,11)(4,12,5,8)]
            [(), (2,5)(3,9)(4,6)(8,11)(10,12), (1,4,6)(2,3,11)(5,8,9)(7,10,12),
             (1,2,3,12)(4,5,10,11)(6,7,8,9)]
            [(), (1,4)(2,8)(3,5)(7,10)(9,11), (1,2,10)(3,6,5)(4,7,8)(9,12,11),
             (1,6,4)(2,11,3)(5,9,8)(7,12,10), (1,5,9,10)(2,12,8,6)(3,4,7,11),
             (1,3)(2,12)(4,10)(5,11)(6,8)(7,9)]
            [(), (1,7)(2,4)(5,6)(8,10)(11,12), (1,10,2)(3,5,6)(4,8,7)(9,11,12),
             (1,12,3,2)(4,11,10,5)(6,9,8,7)]
            [()]
        """
        from sage.combinat.root_system.reflection_group_element import _gap_return
        J_inv = [self._index_set_inverse[j] + 1 for j in J]
        S = str(gap3('ReducedRightCosetRepresentatives(%s,ReflectionSubgroup(%s,%s))' % (self._gap_group._name, self._gap_group._name, J_inv)))
        return [self(w, check=False) for w in _gap_return(S)]

    def simple_root_index(self, i):
        r"""
        Return the index of the simple root `\alpha_i`.

        This is the position of `\alpha_i` in the list of simple roots.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',3])
            sage: [W.simple_root_index(i) for i in W.index_set()]
            [0, 1, 2]
        """
        return self._index_set_inverse[i]

    def bruhat_cone(self, x, y, side='upper', backend='cdd'):
        r"""
        Return the (upper or lower) Bruhat cone associated to the interval ``[x,y]``.

        To a cover relation `v \prec w` in strong Bruhat order you can assign a positive
        root `\beta` given by the unique reflection `s_\beta` such that `s_\beta v = w`.

        The upper Bruhat cone of the interval `[x,y]` is the non-empty, polyhedral cone generated
        by the roots corresponding to `x \prec a` for all atoms `a` in the interval.
        The lower Bruhat cone of the interval `[x,y]` is the non-empty, polyhedral cone generated
        by the roots corresponding to `c \prec y` for all coatoms `c` in the interval.

        INPUT:

        - ``x`` -- an element in the group `W`
        - ``y`` -- an element in the group `W`
        - ``side`` -- (default: ``'upper'``) must be one of the following:

          * ``'upper'`` -- return the upper Bruhat cone of the interval [``x``, ``y``]
          * ``'lower'`` -- return the lower Bruhat cone of the interval [``x``, ``y``]

        - ``backend`` -- string (default: ``'cdd'``); the backend to use to create the polyhedron

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])
            sage: x = W.from_reduced_word([1])
            sage: y = W.w0
            sage: W.bruhat_cone(x, y)
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays

            sage: W = ReflectionGroup(['E',6])
            sage: x = W.one()
            sage: y = W.w0
            sage: W.bruhat_cone(x, y, side='lower')
            A 6-dimensional polyhedron in QQ^6 defined as the convex hull of 1 vertex and 6 rays

        TESTS::

            sage: W = ReflectionGroup(['A',2])
            sage: x = W.one()
            sage: y = W.w0
            sage: W.bruhat_cone(x, y, side='nonsense')
            Traceback (most recent call last):
            ...
            ValueError: side must be either 'upper' or 'lower'

        REFERENCES:

        - [Dy1994]_
        - [JS2021]_
        """
        if side == 'upper':
            roots = [self.reflection_to_positive_root(x * r * x.inverse())
                     for z, r in x.bruhat_upper_covers_reflections()
                     if z.bruhat_le(y)]
        elif side == 'lower':
            roots = [self.reflection_to_positive_root(y * r * y.inverse())
                     for z, r in y.bruhat_lower_covers_reflections()
                     if x.bruhat_le(z)]
        else:
            raise ValueError("side must be either 'upper' or 'lower'")

        from sage.geometry.polyhedron.constructor import Polyhedron
        if self.is_crystallographic():
            return Polyhedron(vertices=[[0] * self.rank()],
                              rays=roots,
                              ambient_dim=self.rank(),
                              backend=backend)
        if backend == 'cdd':
            from warnings import warn
            warn("Using floating point numbers for roots of unity. This might cause numerical errors!")
            from sage.rings.real_double import RDF as base_ring
        else:
            from sage.rings.qqbar import AA as base_ring

        return Polyhedron(vertices=[[0] * self.rank()],
                          rays=roots,
                          ambient_dim=self.rank(),
                          base_ring=base_ring,
                          backend=backend)

    class Element(RealReflectionGroupElement, ComplexReflectionGroup.Element):

        @cached_in_parent_method
        def right_coset_representatives(self):
            r"""
            Return the right coset representatives of ``self``.

            EXAMPLES::

                sage: W = ReflectionGroup(['A',2])
                sage: for w in W:
                ....:     rcr = w.right_coset_representatives()
                ....:     print("%s %s"%(w.reduced_word(),
                ....:                    [v.reduced_word() for v in rcr]))
                [] [[], [2], [1], [2, 1], [1, 2], [1, 2, 1]]
                [2] [[], [2], [1]]
                [1] [[], [1], [1, 2]]
                [1, 2] [[]]
                [2, 1] [[]]
                [1, 2, 1] [[], [2], [2, 1]]
            """
            W = self.parent()
            T = W.reflections()
            T_fix = [i + 1 for i in T.keys()
                     if self.fix_space().is_subspace(T[i].fix_space())]
            S = str(gap3('ReducedRightCosetRepresentatives(%s,ReflectionSubgroup(%s,%s))' % (W._gap_group._name, W._gap_group._name, T_fix)))
            from sage.combinat.root_system.reflection_group_element import _gap_return
            return [W(w, check=False) for w in _gap_return(S)]

        def left_coset_representatives(self):
            r"""
            Return the left coset representatives of ``self``.

            .. SEEALSO:: :meth:`right_coset_representatives`

            EXAMPLES::

                sage: W = ReflectionGroup(['A',2])
                sage: for w in W:
                ....:     lcr = w.left_coset_representatives()
                ....:     print("%s %s"%(w.reduced_word(),
                ....:                    [v.reduced_word() for v in lcr]))
                [] [[], [2], [1], [1, 2], [2, 1], [1, 2, 1]]
                [2] [[], [2], [1]]
                [1] [[], [1], [2, 1]]
                [1, 2] [[]]
                [2, 1] [[]]
                [1, 2, 1] [[], [2], [1, 2]]
            """
            return [ (~w) for w in self.right_coset_representatives() ]


class IrreducibleRealReflectionGroup(RealReflectionGroup, IrreducibleComplexReflectionGroup):
    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: for i in [2..7]: ReflectionGroup(["I", i])
            Reducible real reflection group of rank 2 and type A1 x A1
            Irreducible real reflection group of rank 2 and type A2
            Irreducible real reflection group of rank 2 and type C2
            Irreducible real reflection group of rank 2 and type I2(5)
            Irreducible real reflection group of rank 2 and type G2
            Irreducible real reflection group of rank 2 and type I2(7)
        """
        type_str = self._irrcomp_repr_(self._type[0])
        return 'Irreducible real reflection group of rank %s and type %s' % (self._rank,type_str)

    class Element(RealReflectionGroup.Element, IrreducibleComplexReflectionGroup.Element):
        pass
