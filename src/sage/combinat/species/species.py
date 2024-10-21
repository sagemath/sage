# sage_setup: distribution = sagemath-combinat
"""
Combinatorial Species

This file defines the main classes for working with combinatorial
species, operations on them, as well as some implementations of
basic species required for other constructions.

This code is based on the work of Ralf Hemmecke and Martin Rubey's
Aldor-Combinat, which can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/aldor/combinat/index.html.
In particular, the relevant section for this file can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/AldorCombinat/combinatse8.html.

Weighted Species:

As a first application of weighted species, we count unlabeled
ordered trees by total number of nodes and number of internal
nodes. To achieve this, we assign a weight of `1` to the
leaves and of `q` to internal nodes::

    sage: q = QQ['q'].gen()
    sage: leaf = species.SingletonSpecies()
    sage: internal_node = species.SingletonSpecies(weight=q)
    sage: L = species.LinearOrderSpecies(min=1)
    sage: T = species.CombinatorialSpecies(min=1)
    sage: T.define(leaf + internal_node*L(T))
    sage: T.isotype_generating_series()[0:6]                                            # needs sage.modules
    [0, 1, q, q^2 + q, q^3 + 3*q^2 + q, q^4 + 6*q^3 + 6*q^2 + q]

Consider the following::

    sage: T.isotype_generating_series().coefficient(4)                                  # needs sage.modules
    q^3 + 3*q^2 + q

This means that, among the trees on `4` nodes, one has a
single internal node, three have two internal nodes, and one has
three internal nodes.
"""
# ****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from .generating_series import OrdinaryGeneratingSeriesRing, ExponentialGeneratingSeriesRing, CycleIndexSeriesRing
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.combinat.species.misc import accept_size
from sage.combinat.species.structure import StructuresWrapper, IsotypesWrapper
from functools import reduce


class GenericCombinatorialSpecies(SageObject):
    def __init__(self, min=None, max=None, weight=None):
        """
        TESTS::

            sage: P = species.PermutationSpecies(size=3)
            sage: P._weight
            1
            sage: P._min
            3
            sage: P._max
            4
        """
        self._weight = weight if weight is not None else QQ(1)
        self._min = min
        self._max = max

    def __hash__(self):
        """
        Return a hash of the unique info tuple.

        EXAMPLES::

            sage: hash(species.SetSpecies()) #random
            -152204909943771174
        """
        return hash(self._unique_info())

    def _unique_info(self):
        """
        Return a tuple which should uniquely identify the species.

        EXAMPLES::

            sage: species.SetSpecies()._unique_info()
            (<class 'sage.combinat.species.set_species.SetSpecies'>, None, None, 1)
            sage: species.SingletonSpecies()._unique_info()
            (<class 'sage.combinat.species.characteristic_species.SingletonSpecies'>,
             None,
             None,
             1)

        ::

            sage: X = species.SingletonSpecies()
            sage: Y = X + X
            sage: Y._unique_info()
            (<class 'sage.combinat.species.sum_species.SumSpecies'>,
             None,
             None,
             1,
             Singleton species,
             Singleton species)
        """
        info = (self.__class__, self._min, self._max, self._weight)
        if hasattr(self, "_state_info") and self._state_info:
            info += tuple(self._state_info)
        return info

    def __eq__(self, x):
        """
        Test equality between two species.

        EXAMPLES::

            sage: X = species.SingletonSpecies()
            sage: X + X == X + X
            True
            sage: X == X
            True
            sage: X == species.EmptySetSpecies()
            False
            sage: X == X*X
            False

        ::

            sage: X = species.SingletonSpecies()
            sage: E = species.EmptySetSpecies()
            sage: L = CombinatorialSpecies()
            sage: L.define(E+X*L)
            sage: K = CombinatorialSpecies()
            sage: K.define(E+X*L)
            sage: L == K
            True
        """
        if not isinstance(x, GenericCombinatorialSpecies):
            return False
        return self._unique_info() == x._unique_info()

    def __ne__(self, other):
        """
        Check whether ``self`` and ``other`` are not equal.

        EXAMPLES::

            sage: X = species.SingletonSpecies()
            sage: X + X == X + X
            True
            sage: X != X
            False
            sage: X != species.EmptySetSpecies()
            True
            sage: X != X*X
            True

            sage: X = species.SingletonSpecies()
            sage: E = species.EmptySetSpecies()
            sage: L = CombinatorialSpecies()
            sage: L.define(E+X*L)
            sage: K = CombinatorialSpecies()
            sage: K.define(E+X*L)
            sage: L != K
            False
        """
        return not (self == other)

    def __getstate__(self):
        r"""
        This is used during the pickling process and returns a dictionary
        of the data needed to create this object during the unpickling
        process. It returns an (\*args, \*\*kwds) tuple which is to be
        passed into the constructor for the class of this species. Any
        subclass should define a ``_state_info`` list for any arguments which
        need to be passed in the constructor.

        EXAMPLES::

            sage: C = species.CharacteristicSpecies(5)
            sage: args, kwds = C.__getstate__()
            sage: args
            {0: 5}
            sage: sorted(kwds.items())
            [('max', None), ('min', None), ('weight', 1)]
        """
        kwds = {'weight': self._weight, 'min': self._min, 'max': self._max}
        try:
            return (dict(enumerate(self._state_info)), kwds)
        except AttributeError:
            return ({}, kwds)

    def __setstate__(self, state):
        """
        This is used during unpickling to recreate this object from the
        data provided by the ``__getstate__`` method.

        TESTS::

            sage: C2 = species.CharacteristicSpecies(2)
            sage: C4 = species.CharacteristicSpecies(4)
            sage: C2
            Characteristic species of order 2
            sage: C2.__setstate__(C4.__getstate__()); C2
            Characteristic species of order 4
        """
        args_dict, kwds = state
        self.__class__.__init__(self, *[args_dict[i] for i in range(len(args_dict))], **kwds)

    def weighted(self, weight):
        """
        Return a version of this species with the specified weight.

        EXAMPLES::

            sage: t = ZZ['t'].gen()
            sage: C = species.CycleSpecies(); C
            Cyclic permutation species
            sage: C.weighted(t)
            Cyclic permutation species with weight=t
        """
        args_dict, kwds = self.__getstate__()
        kwds.update({'weight': weight})
        return self.__class__(*[args_dict[i] for i in range(len(args_dict))], **kwds)

    def __repr__(self):
        """
        Return a string representation of this species.

        EXAMPLES::

            sage: CombinatorialSpecies()
            Combinatorial species

        ::

            sage: species.SetSpecies()
            Set species
            sage: species.SetSpecies(min=1)
            Set species with min=1
            sage: species.SetSpecies(min=1, max=4)
            Set species with min=1, max=4
            sage: t = ZZ['t'].gen()
            sage: species.SetSpecies(min=1, max=4, weight=t)
            Set species with min=1, max=4, weight=t
        """
        if hasattr(self, "_name"):
            name = self._name if isinstance(self._name, str) else self._name()
        else:
            name = "Combinatorial species"

        options = []

        if self._min is not None:
            options.append('min=%s' % self._min)
        if self._max is not None:
            options.append('max=%s' % self._max)
        if self._weight != 1:
            options.append('weight=%s' % self._weight)

        if options:
            name += " with " + ", ".join(options)

        return name

    def __add__(self, g):
        """
        Return the sum of ``self`` and ``g``.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: F = P + P; F
            Sum of (Permutation species) and (Permutation species)
            sage: F.structures([1,2]).list()
            [[1, 2], [2, 1], [1, 2], [2, 1]]
        """
        from .sum_species import SumSpecies
        if not isinstance(g, GenericCombinatorialSpecies):
            raise TypeError("g must be a combinatorial species")
        return SumSpecies(self, g)

    sum = __add__

    def __mul__(self, g):
        """
        Return the product of ``self`` and ``g``.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: F = P * P; F
            Product of (Permutation species) and (Permutation species)
        """
        from .product_species import ProductSpecies
        if not isinstance(g, GenericCombinatorialSpecies):
            raise TypeError("g must be a combinatorial species")
        return ProductSpecies(self, g)

    product = __mul__

    def __call__(self, g):
        """
        EXAMPLES::

            sage: S = species.SetSpecies()
            sage: S(S)
            Composition of (Set species) and (Set species)
        """
        from .composition_species import CompositionSpecies
        if not isinstance(g, GenericCombinatorialSpecies):
            raise TypeError("g must be a combinatorial species")
        return CompositionSpecies(self, g)

    composition = __call__

    def functorial_composition(self, g):
        """
        Return the functorial composition of ``self`` with ``g``.

        EXAMPLES::

            sage: E = species.SetSpecies()
            sage: E2 = E.restricted(min=2, max=3)
            sage: WP = species.SubsetSpecies()
            sage: P2 = E2*E
            sage: G = WP.functorial_composition(P2)
            sage: G.isotype_generating_series()[0:5]                                    # needs sage.modules
            [1, 1, 2, 4, 11]
        """
        from .functorial_composition_species import FunctorialCompositionSpecies
        if not isinstance(g, GenericCombinatorialSpecies):
            raise TypeError("g must be a combinatorial species")
        return FunctorialCompositionSpecies(self, g)

    @accept_size
    def restricted(self, min=None, max=None):
        """
        Return the restriction of the species.

        INPUT:

        - ``min`` -- (optional) integer

        - ``max`` -- (optional) integer

        EXAMPLES::

            sage: S = species.SetSpecies().restricted(min=3); S
            Set species with min=3
            sage: S.structures([1,2]).list()
            []
            sage: S.generating_series()[0:5]
            [0, 0, 0, 1/6, 1/24]
        """
        kwargs = {'min': self._min if min is None else min,
                  'max': self._max if max is None else max,
                  'weight': self._weight}
        return self.__class__(**kwargs)

    def structures(self, labels, structure_class=None):
        """
        EXAMPLES::

            sage: F = CombinatorialSpecies()
            sage: F.structures([1,2,3]).list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return StructuresWrapper(self, labels, structure_class)

    def isotypes(self, labels, structure_class=None):
        """
        EXAMPLES::

            sage: F = CombinatorialSpecies()
            sage: F.isotypes([1,2,3]).list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return IsotypesWrapper(self, labels, structure_class=structure_class)

    def _check(self, n=5):
        """
        Return ``True`` if the number of structures and isomorphism types
        generated is the same as the number found from the generating
        series.

        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: P._check()                                                            # needs sage.libs.flint
            True
        """
        st = self.structures(range(n))
        it = self.isotypes(range(n))

        try:
            return (len(st.list()) == st.cardinality() and
                    len(it.list()) == it.cardinality())
        except NotImplementedError:
            return False

    def __pow__(self, n):
        r"""
        Return this species to the power `n`.

        This uses a binary exponentiation algorithm to perform the
        powering.

        EXAMPLES::

            sage: One = species.EmptySetSpecies()
            sage: X = species.SingletonSpecies()
            sage: X^2
            Product of (Singleton species) and (Singleton species)
            sage: X^5
            Product of (Singleton species) and (Product of (Product of
            (Singleton species) and (Singleton species)) and (Product
            of (Singleton species) and (Singleton species)))

            sage: (X^2).generating_series()[0:4]
            [0, 0, 1, 0]
            sage: (X^3).generating_series()[0:4]
            [0, 0, 0, 1]
            sage: ((One+X)^3).generating_series()[0:4]
            [1, 3, 3, 1]
            sage: ((One+X)^7).generating_series()[0:8]
            [1, 7, 21, 35, 35, 21, 7, 1]

            sage: x = QQ[['x']].gen()
            sage: coeffs = ((1+x+x+x**2)**25+O(x**10)).padded_list()
            sage: T = ((One+X+X+X^2)^25)
            sage: T.generating_series()[0:10] == coeffs
            True
            sage: X^1 is X
            True
            sage: A = X^32
            sage: A.digraph()                                                           # needs sage.graphs
            Multi-digraph on 6 vertices

        TESTS::

            sage: X**(-1)
            Traceback (most recent call last):
            ...
            ValueError: only positive exponents are currently supported
        """
        from sage.rings.integer import Integer
        import operator
        n = Integer(n)
        if n <= 0:
            raise ValueError("only positive exponents are currently supported")
        digits = n.digits(2)
        squares = [self]
        for i in range(len(digits) - 1):
            squares.append(squares[-1] * squares[-1])
        return reduce(operator.mul, (s for i, s in zip(digits, squares)
                                     if i != 0))

    def _get_series(self, series_ring_class, prefix, base_ring=None):
        """
        Return the generating / isotype generating / cycle index series
        ring. The purpose of this method is to restrict the result of
        _series_helper to ``self._min`` and ``self._max``.

        EXAMPLES::

            sage: P = species.PermutationSpecies(min=2, max=4)
            sage: P.generating_series()[0:8] #indirect doctest
            [0, 0, 1, 1, 0, 0, 0, 0]
        """
        series = self._series_helper(series_ring_class, prefix, base_ring=base_ring)
        # We need to restrict the series based on the min
        # and max of this species.  Note that if min and max
        # are both None (as in the default case), then the restrict
        # method will just return series.
        if self._min is None and self._max is None:
            return series
        return series.parent()(lambda n: series[n],
                               valuation=self._min, degree=self._max)

    def _series_helper(self, series_ring_class, prefix, base_ring=None):
        """
        This code handles much of the common work involved in getting the
        generating series for this species (such as determining the
        correct base ring to pass down to the subclass, determining which
        method on the subclass to call to get the series object, etc.)

        INPUT:

        - ``series_ring_class`` -- a class for the series ring such as
          ``ExponentialGeneratingSeriesRing``, etc.

        - ``prefix`` -- the string prefix associated with the generating series
          such as "cis" for the cycle index series. This prefix appears in the
          methods that are implemented in the subclass.

        - ``base_ring`` -- the ring in which the coefficients of the generating
          series live. If it is not specified, then it is determined by the
          weight of the species.

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import OrdinaryGeneratingSeriesRing
            sage: S = species.SetSpecies()
            sage: itgs = S._series_helper(OrdinaryGeneratingSeriesRing, "itgs")
            sage: itgs[:3]
            [1, 1, 1]

        ::

            sage: itgs = S._series_helper(OrdinaryGeneratingSeriesRing, "itgs", base_ring=RDF)
            sage: itgs[:3]
            [1.0, 1.0, 1.0]
        """
        prefix = "_" + prefix

        # Get the base ring
        if base_ring is None:
            base_ring = self.weight_ring()
        else:
            # The specified base ring must have maps from both
            # the rational numbers and the weight ring
            if not base_ring.has_coerce_map_from(QQ):
                raise ValueError("specified base ring does not contain the rationals")
            if not base_ring.has_coerce_map_from(self.weight_ring()):
                raise ValueError("specified base ring is incompatible with the weight ring of self")

        series_ring = series_ring_class(base_ring)

        # Try to return things like self._gs(base_ring)
        # This is used when the subclass wants to just
        # handle creating the generating series itself;
        # for example, returning the exponential of a
        # generating series.
        try:
            return getattr(self, prefix)(series_ring, base_ring)
        except AttributeError:
            pass

        # Try to return things like self._gs_callable(base_ring).
        # This is used when the subclass just provides a callable
        # for the coefficients of the generating series.  Optionally,
        # the subclass can specify the order of the series.
        try:
            callable = getattr(self, prefix + "_callable")
            try:
                return series_ring(lambda n: callable(base_ring, n), valuation=self._order())
            except AttributeError:
                return series_ring(lambda n: callable(base_ring, n))
        except AttributeError:
            pass

        # Try to use things like self._gs_term(base_ring).
        # This is used when the generating series is just a single
        # term.
        try:
            return series_ring(getattr(self, prefix + "_term")(base_ring),
                                    self._order())
        except AttributeError:
            pass

        # Try to use things like self._gs_list(base_ring).
        # This is used when the coefficients of the generating series
        # can be given by a finite list with the last coefficient repeating.
        # The generating series with all ones coefficients is generated this
        # way.
        try:
            return series_ring(lambda n: getattr(self, prefix + "_list")(base_ring, n))
        except AttributeError:
            pass

        raise NotImplementedError

    @cached_method
    def generating_series(self, base_ring=None):
        r"""
        Return the generating series for this species.

        This is an exponential generating series so the `n`-th
        coefficient of the series corresponds to the number of labeled
        structures with `n` labels divided by `n!`.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: g = P.generating_series()
            sage: g[:4]
            [1, 1, 1, 1]
            sage: g.counts(4)
            [1, 1, 2, 6]
            sage: P.structures([1,2,3]).list()
            [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
            sage: len(_)
            6
        """
        return self._get_series(ExponentialGeneratingSeriesRing, "gs", base_ring)

    @cached_method
    def isotype_generating_series(self, base_ring=None):
        r"""
        Return the isotype generating series for this species.

        The `n`-th coefficient of this series corresponds to the number
        of isomorphism types for the structures on `n` labels.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: g = P.isotype_generating_series()
            sage: g[0:4]                                                                # needs sage.libs.flint
            [1, 1, 2, 3]
            sage: g.counts(4)                                                           # needs sage.libs.flint
            [1, 1, 2, 3]
            sage: P.isotypes([1,2,3]).list()                                            # needs sage.libs.flint
            [[2, 3, 1], [2, 1, 3], [1, 2, 3]]
            sage: len(_)                                                                # needs sage.libs.flint
            3
        """
        return self._get_series(OrdinaryGeneratingSeriesRing, "itgs", base_ring)

    @cached_method
    def cycle_index_series(self, base_ring=None):
        r"""
        Return the cycle index series for this species.

        The cycle index series is a sequence of symmetric functions.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: g = P.cycle_index_series()                                            # needs sage.modules
            sage: g[0:4]                                                                # needs sage.modules
            [p[], p[1], p[1, 1] + p[2], p[1, 1, 1] + p[2, 1] + p[3]]
        """
        return self._get_series(CycleIndexSeriesRing, "cis", base_ring)

    def is_weighted(self):
        """
        Return ``True`` if this species has a nontrivial weighting associated
        with it.

        EXAMPLES::

            sage: C = species.CycleSpecies()
            sage: C.is_weighted()
            False
        """
        return self._weight != 1

    def weight_ring(self):
        """
        Return the ring in which the weights of this species occur.

        By default, this is just the field of rational numbers.

        EXAMPLES::

            sage: species.SetSpecies().weight_ring()
            Rational Field
        """
        if self.is_weighted():
            return self._weight.parent()
        else:
            return QQ

    def _common_parent(self, parents):
        """
        Return a parent that contains all the parents
        in the given list of parents.

        EXAMPLES::

            sage: C = species.CombinatorialSpecies()
            sage: C._common_parent([QQ, ZZ['t']])
            Univariate Polynomial Ring in t over Rational Field
        """
        assert parents
        from sage.structure.element import get_coercion_model
        cm = get_coercion_model()

        common = parents[0]
        for p in parents[1:]:
            common = cm.explain(common, p, verbosity=0)
            if common is None:
                raise ValueError("unable to find a common parent")
        return common

    def digraph(self):
        """
        Return a directed graph where the vertices are the individual
        species that make up this one.

        EXAMPLES::

            sage: X = species.SingletonSpecies()
            sage: B = species.CombinatorialSpecies()
            sage: B.define(X+B*B)
            sage: g = B.digraph(); g                                                    # needs sage.graphs
            Multi-digraph on 4 vertices

            sage: sorted(g, key=str)                                                    # needs sage.graphs
            [Combinatorial species,
             Product of (Combinatorial species) and (Combinatorial species),
             Singleton species,
             Sum of (Singleton species) and
              (Product of (Combinatorial species) and (Combinatorial species))]

            sage: d = {sp: i for i, sp in enumerate(g)}                                 # needs sage.graphs
            sage: g.relabel(d)                                                          # needs sage.graphs
            sage: g.canonical_label().edges(sort=True)                                  # needs sage.graphs
            [(0, 3, None), (2, 0, None), (2, 0, None), (3, 1, None), (3, 2, None)]
        """
        from sage.graphs.digraph import DiGraph
        d = DiGraph(multiedges=True)
        self._add_to_digraph(d)
        return d

    def _add_to_digraph(self, d):
        """
        Add this species as a vertex to the digraph ``d`` along with any
        'children' of this species. For example, sum species would add
        itself as a vertex and an edge between itself and each of its
        summands.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: d = DiGraph(multiedges=True)
            sage: X = species.SingletonSpecies()
            sage: X._add_to_digraph(d); d
            Multi-digraph on 1 vertex
            sage: (X+X)._add_to_digraph(d); d
            Multi-digraph on 2 vertices
            sage: d.edges(sort=True)
            [(Sum of (Singleton species) and (Singleton species), Singleton species, None),
             (Sum of (Singleton species) and (Singleton species), Singleton species, None)]
        """
        d.add_vertex(self)

        if not hasattr(self, "_state_info"):
            return

        for child in self._state_info:
            if not isinstance(child, GenericCombinatorialSpecies):
                continue
            d.add_edge(self, child)
            child._add_to_digraph(d)

    def algebraic_equation_system(self):
        """
        Return a system of algebraic equations satisfied by this species.

        The nodes are numbered in the order that they appear as vertices of
        the associated digraph.

        EXAMPLES::

            sage: B = species.BinaryTreeSpecies()
            sage: B.algebraic_equation_system()                                         # needs sage.graphs
            [-node3^2 + node1, -node1 + node3 + (-z)]

        ::

            sage: sorted(B.digraph().vertex_iterator(), key=str)                        # needs sage.graphs
            [Combinatorial species with min=1,
             Product of (Combinatorial species with min=1)
                    and (Combinatorial species with min=1),
             Singleton species,
             Sum of (Singleton species)
                and (Product of (Combinatorial species with min=1)
                and (Combinatorial species with min=1))]

        ::

            sage: B.algebraic_equation_system()[0].parent()                             # needs sage.graphs
            Multivariate Polynomial Ring in node0, node1, node2, node3 over
             Fraction Field of Univariate Polynomial Ring in z over Rational Field
        """
        d = self.digraph()

        Qz = QQ['z'].fraction_field()

        # Generate the variable names and the corresponding polynomial rings
        var_names = ["node%s" % i for i in range(d.num_verts())]
        R = Qz[", ".join(var_names)]
        R_gens_dict = R.gens_dict()

        # A dictionary mapping the nodes to variables
        vertices = sorted(d.vertex_iterator(), key=str)
        var_mapping = {node: R_gens_dict[name]
                       for node, name in zip(vertices, var_names)}
        var_mapping['z'] = Qz.gen()

        eqns = []
        subs = {}
        for species in vertices:
            try:
                eqn = species._equation(var_mapping)
                if eqn in Qz or eqn in R.gens():
                    subs[var_mapping[species]] = eqn
                else:
                    eqns.append(var_mapping[species] - eqn)
            except AttributeError:
                raise NotImplementedError
        return [eq.subs(subs) for eq in eqns]
