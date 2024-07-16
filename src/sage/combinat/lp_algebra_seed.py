r"""
Laurent Phenomenon Algebra Seeds

This class implements seeds and their mutations for Lam and Pylyavskyy's *Laurent
phenomenon algebras* (LP algebras).

Fixing a unique factorization domain `A`, a pair `(\mathbf{x}, \mathbf{f})` is
said to be an *LP seed* if `\mathbf{x}=\{x_1, ..., x_n\}` is a transcendence
basis for the field of rational functions in `n` independent variables over
`\text{Frac}(A)`, and `\mathbf{f} = \{f_1, ..., f_n\}` is a collection of
irreducible polynomials over `A` encoding the exchange relations.

One can view LP seeds and their corresponding LP algebras as a vast
generalisation of Fomin and Zelevinsky's cluster algebras. This module provides
basic functionality for investigating their properties.

AUTHORS:

- Oliver Daisey (2023-03-20): initial version
"""

# ****************************************************************************
#       Copyright (C) 2023 Oliver Daisey <oliver.j.daisey at durham.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from sage.arith.misc import factor, gcd
from sage.graphs.graph import Graph
from sage.rings.fraction_field import FractionField
from sage.rings.infinity import infinity
from sage.rings.integer_ring import IntegerRing_class, ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import RationalField
from sage.structure.sage_object import SageObject
from sage.misc.latex import latex
from sage.symbolic.ring import SR


class LPASeed(SageObject):
    r"""
    Initialise a Laurent phenomenon algebra seed.

    INPUT:

    - ``data`` -- can be one of the following:

        * dict - dictionary of initial variable names to their
            corresponding exchange polynomials
        * LPASeed object

    - ``coefficients`` -- tuple of symbolic variables (default: ``()``); the
        labels of, if any, the coefficients of the exchange polynomials. If no
        coefficients are provided, the module attempts to detect them from the
        input data.

    - ``base_ring`` -- unique factorisation domain (default: ``ZZ``);
        the ring which we take the exchange polynomials over; currently
        only supports ``ZZ`` or ``QQ``

    EXAMPLES:

    This example initialises a linear Laurent phenomenon algebra in two
    variables::

        sage: var('x1,x2')
        (x1, x2)
        sage: S = LPASeed({x1: 1 + x2, x2: 1 + x1})
        sage: S
        A seed with cluster variables [x1, x2] and exchange polynomials [x2 + 1, x1 + 1]

    We add some coefficients to get the generic linear LP algebra in three
    variables::

        sage: var('x1,x2,x3,a0,a2,a3,b0,b1,b3,c0,c1,c2')
        (x1, x2, x3, a0, a2, a3, b0, b1, b3, c0, c1, c2)
        sage: S = LPASeed({x1: a0 + a2*x2 + a3*x3, x2: b0 + b1*x1 + b3*x3, x3: c0 + c1*x1 + c2*x2}, coefficients=[a0,a2,a3,b0,b1,b3,c0,c1,c2],base_ring=ZZ)
        sage: S
        A seed with cluster variables [x1, x2, x3] and exchange polynomials [x2*a2 + x3*a3 + a0, x1*b1 + x3*b3 + b0, x1*c1 + x2*c2 + c0]

    More complicated polynomials are allowed, as long as they are
    irreducible::

        sage: var('x1,x2,x3')
        (x1, x2, x3)
        sage: S = LPASeed({x1: 1 + x2*x3^2 + 4*x3^3, x2: 2 - x1^2, x3: 4 + x1^3*x2^2 - 3*x1})
        sage: S
        A seed with cluster variables [x1, x2, x3] and exchange polynomials [x2*x3^2 + 4*x3^3 + 1, -x1^2 + 2, x1^3*x2^2 - 3*x1 + 4]

    Nonirreducible polynomials will raise an exception::

        sage: var('x1, x2')
        (x1, x2)
        sage: S = LPASeed({x1: 4 - x2^2, x2: 1 + x1})
        Traceback (most recent call last):
        ...
        ValueError: (LP2) fail: -x2^2 + 4 is not irreducible over Integer Ring

    Different base rings are allowed::

        sage: var('x1,x2')
        (x1, x2)
        sage: S = LPASeed({x1: 1 + x2, x2: 1 + x1}, base_ring=QQ)
        sage: S
        A seed with cluster variables [x1, x2] and exchange polynomials [x2 + 1, x1 + 1]
    """

    def __init__(self, data, coefficients=(), base_ring=ZZ):
        r"""
        Initialize an LP seed.

        EXAMPLES::

            sage: var('x1, x2, x3')
            (x1, x2, x3)
            sage: S = LPASeed({x1: 1 + x2 + x3, x2: 1 + x1 + x3, x3: 1 + x1 + x2})
            sage: TestSuite(S).run()
        """
        # unpack supplied data

        # make a copy
        if isinstance(data, LPASeed):

            self._names = data._names
            self._base_ring = data._base_ring
            self._coefficients = data._coefficients
            self._ambient_field = data._ambient_field
            self._polynomial_ring = data._polynomial_ring
            self._rank = data._rank
            self._exchange_polys = copy(data._exchange_polys)
            self._laurent_polys = copy(data._laurent_polys)
            self._cluster_vars = copy(data._cluster_vars)
            self._mutation_sequence = copy(data._mutation_sequence)
            return

        # we assume we have dictionary of variable/poly pairs
        if not isinstance(data, dict):
            raise TypeError("the input must be a dict or a LPASeed")

        names = tuple(data)  # ensure names are immutable

        # if we are not supplied any coefficients, check input data
        if not coefficients:
            for name in names:
                # variables present in this polynomial
                name_poly_vars = set(SR(data[name]).variables())
                # add the variables that are not already names
                coefficients += tuple(name_poly_vars.difference(set(names)))

        coefficients = tuple(set(coefficients))  # coefficients are immutable and have no repeats
        exchange_polys = list(data.values())

        self._names = names

        # construct the ambient ring from input data

        # first we get generators for ambient ring
        variables = names + coefficients

        # now attempt to build polynomial ring from this data
        if not isinstance(base_ring, (IntegerRing_class, RationalField)):
            raise NotImplementedError("only ZZ or QQ is supported")

        self._base_ring = base_ring
        self._polynomial_ring = PolynomialRing(self._base_ring,
                                               names=variables)

        self._ambient_field = FractionField(self._polynomial_ring)

        self._rank = len(self._names)

        # we get initial cluster variables by casting initial variables as
        # rational functions
        self._cluster_vars = [self._ambient_field(self._names[i])
                              for i in range(self._rank)]

        # take what input data we were given and try to use it to
        # construct polynomials
        self._exchange_polys = []
        for i in range(self._rank):
            self._exchange_polys.append(
                self._polynomial_ring(exchange_polys[i]))

        self._coefficients = coefficients
        self._check_seed()

        # ensure mutation sequence / seed list initialises as empty lists
        self._mutation_sequence = []

        # begin with correct Laurent polynomials
        self._laurent_polys = self._exchange_polys.copy()
        self._compute_laurent()

    def _check_seed(self):
        r"""
        Perform some mathematical checks on this seed.

        This includes checking that each exchange polynomial does not depend on
        its corresponding cluster variable, that each exchange polynomial is
        irreducible, and finally that no cluster variable divides any exchange
        polynomial.

        TESTS::

            sage: var('x1, x2')
            (x1, x2)
            sage: LPASeed({x1: 1})
            Traceback (most recent call last):
            ...
            ValueError: (LP2) fail: 1 is not irreducible over Integer Ring

            sage: LPASeed({x1: 0})
            Traceback (most recent call last):
            ...
            ValueError: (LP2) fail: exchange polynomial is zero

            sage: LPASeed({x1: 1 + 2*x2, x2: 3 + 4*x1})
            A seed with cluster variables [x1, x2] and exchange polynomials [2*x2 + 1, 4*x1 + 3]
        """
        # (LP1) check polynomials do not depend on their cluster variable

        for i in range(self._rank):
            currPoly = self._exchange_polys[i]
            varList = currPoly.variables()
            if self._names[i] in varList:
                raise ValueError("(LP1) fail: %s depends on %s" %
                                 (self._exchange_polys[i], self._names[i]))

        # (LP2) check exchange polynomials are irreducible
        # and not divisible by any cluster variable

        for f in self._exchange_polys:
            # check irreducibility
            if f == 0:
                raise ValueError("(LP2) fail: exchange polynomial is zero")
            L = list(factor(f))
            if len(L) > 1 or f.is_unit():
                raise ValueError("(LP2) fail: %s is not irreducible over %s"
                                 % (f, self._base_ring))

            # check not divisible by any variable
            for var in self._names:
                if f % self._polynomial_ring(var) == 0:
                    raise ValueError("(LP2) fail: %s divides %s" % (var, f))

    # a private method to compute the exchange Laurent polynomials of `self`

    def _compute_laurent(self):
        r"""
        Compute the exchange Laurent polynomials of ``self``.

        TESTS:

        This should leave the seed invariant::

            sage: var('x1, x2, a')
            (x1, x2, a)
            sage: S = LPASeed({x1: a, x2: a}, coefficients=[a])
            sage: S._compute_laurent()
            sage: S.laurent_polys()
            [a/x2, a/x1]
        """
        # work with copies as we perform substitutions
        exchange_polys = copy(self._exchange_polys)
        laurent_polys = copy(self._exchange_polys)

        for i in range(self._rank):

            # this is the exchange polynomial we are
            # finding laurent polynomial for
            current_poly = exchange_polys[i]

            # we build the laurent polynomial from scratch by iterating
            # through the other variables, dividing by appropriate variable
            # if necessary
            for j in range(self._rank):
                if i != j:

                    # this is the polynomial we want to check
                    # divisibility for
                    sub_poly = current_poly.subs(
                        **{str(self._names[j]): exchange_polys[j]})

                    # calculate maximal power of exchange_polys[j] that
                    # divides sub_poly
                    counter = 0
                    while True:
                        (q, r) = sub_poly.quo_rem(
                            exchange_polys[j] ** (counter + 1))
                        if r != 0:
                            break
                        counter = counter + 1

                    # divide exchange polynomial by this maximal power
                    laurent_polys[i] = laurent_polys[i] / (
                        self._polynomial_ring(self._names[j]) ** counter)

        # after performing all substitutions, set internal laurent polynomials
        self._laurent_polys = laurent_polys

    # mutates self at ith variable

    def mutate(self, i, inplace=True):
        r"""
        Mutate this seed at the ``i`` th index.

        INPUT:

        - ``i`` -- integer or iterable of integers; the index/indices to mutate
          ``self`` at, where we index from 0

        - ``inplace`` -- boolean (default: ``True``); whether to mutate the
          current instance of the seed, or return a new ``LPASeed`` object

        EXAMPLES:

        We mutate a rank-two Laurent phenomenon algebra at the first index::

            sage: var('x1,x2,a0,a2,b0,b1')
            (x1, x2, a0, a2, b0, b1)
            sage: S = LPASeed({x1: a0 + a2*x2, x2: b0 + b1*x1}, coefficients=[a0,a2,b0,b1])
            sage: S
            A seed with cluster variables [x1, x2] and exchange polynomials [x2*a2 + a0, x1*b1 + b0]
            sage: S.mutate(0, inplace=True)
            A seed with cluster variables [(x2*a2 + a0)/x1, x2] and exchange polynomials [x2*a2 + a0, x1*b0 + a0*b1]

        Mutating at the same index is an involution::

            sage: var('x1,x2,x3')
            (x1, x2, x3)
            sage: S = LPASeed({x1: 1 + x2*x3, x2: 1 + x1^2 + x3^2, x3: 1 + x1 + x2})
            sage: T = S.mutate([2,2], inplace=False)
            sage: T == S
            True
        """
        # input preprocessing

        if not inplace:
            self = LPASeed(self)

        # is our input data iterable?
        if hasattr(i, '__iter__'):
            for index in i:
                if index not in range(self._rank):
                    raise IndexError('iterable %s contains a nonvalid index %s'
                                     % (i, index))
                self.mutate(index)
            return self

        elif i not in range(self._rank):
            raise IndexError('did not pass in a valid integer index')

        # Mutate cluster variables:

        xi = self._names[i]  # the variable we are mutating at

        # get the new cluster variable by applying exchange relation
        # need to cast names as polynomials for subs method
        d = {str(self._names[j]): self._cluster_vars[j] for j in range(
            self._rank)}
        exchange_laurent_poly = self._laurent_polys[i].subs(**d)
        self._cluster_vars[i] = exchange_laurent_poly / self._cluster_vars[i]

        # Mutate exchange polynomials:

        for j in range(self._rank):
            if j != i and xi in self._exchange_polys[j].variables():

                xj = self._names[j]  # the variable corresponding to this poly

                # MUTATION ALGORITHM:

                # SUBSTITUTION:
                h = self._laurent_polys[i].subs(**{str(xj): 0})
                G = self._exchange_polys[j].subs(**{str(xi): (h / xi)})
                G = self._polynomial_ring(G.numerator())

                # CANCELLATION:
                G_factors = list(G.factor())
                H = 1  # this will be G with all common factors with h removed
                for factor in G_factors:
                    if gcd(self._ambient_field(h).numerator(), factor[0]) == 1:
                        H = H * (factor[0]) ** (factor[1])

                # NORMALISATION:
                self._exchange_polys[j] = H.numerator()

        # after completing all transformations, update the sequence of mutations
        self._mutation_sequence.append(i)

        # make sure we have up-to-date laurent polynomials
        self._compute_laurent()

        return self

    def is_mutation_equivalent(self, other_seed):
        r"""
        Return whether this seed and ``other_seed`` are mutation equivalent.

        INPUT:

        - ``other_seed`` -- ``LPASeed`` object; the seed we wish to compare to

        OUTPUT:

        - ``True`` if the two seeds are mutation equivalent, and ``False``
          otherwise

        EXAMPLES:

        A seed is mutation equivalent to any of its mutations::

            sage: var('x1,x2,x3')
            (x1, x2, x3)
            sage: S = LPASeed({x1: 1 + x2 + x3, x2: 1 + x1 + x3, x3: 1 + x1 + x2})
            sage: T = S.mutate(0, inplace=False)
            sage: S.is_mutation_equivalent(T)
            True
        """
        if not isinstance(other_seed, LPASeed):
            raise ValueError('%s is not a seed!' % (other_seed))

        for i in range(self._rank):
            seed_test = LPASeed(self)
            seed_test.mutate(i)
            if seed_test == other_seed:
                return True
        return False

    def mutation_class_iter(self, depth=infinity, verbose=False,
                            return_paths=False, algorithm='BFS'):
        r"""
        Return an iterator for the mutation class of this seed.

        INPUT:

        - ``depth`` -- Integer (default: ``infinity``); only return seeds at
          most ``depth`` mutations away from the initial seed

        - ``verbose`` -- Boolean (default: ``False``); if ``True``, the
          current depth of recursion for the chosen algorithm is shown while
          computing

        - ``return_paths`` -- Boolean (default: ``False``); if ``True``, a
          path of mutations from ``self`` to the given seed is returned as well

        - ``algorithm`` -- String (default: ``'BFS'``); the search algorithm to
          find new seeds. Currently supported options::

          * 'BFS' - breadth-first search
          * 'DFS' - depth-first search

        EXAMPLES:

        We iterate over the mutation class for a rank two seed::

            sage: var('x1,x2')
            (x1, x2)
            sage: S = LPASeed({x1: 1 + x2, x2: 1 + x1})
            sage: t = S.mutation_class_iter()
            sage: for seed in t: print(seed.cluster())
            [x1, x2]
            [(x2 + 1)/x1, x2]
            [x1, (x1 + 1)/x2]
            [(x2 + 1)/x1, (x1 + x2 + 1)/(x1*x2)]
            [(x1 + x2 + 1)/(x1*x2), (x1 + 1)/x2]

        Non finite-type works if we specify a fixed depth, but seeds can get big
        rather quickly::

            sage: var('x1,x2')
            (x1, x2)
            sage: S = LPASeed({x1: 1 + x2 + x2^2, x2: 1 + x1 + x1^2})
            sage: t = S.mutation_class_iter(depth=15,algorithm='DFS')
            sage: for seed in t: print(seed.cluster()[0].denominator()) # long time
            1
            x1
            1
            x1*x2^2
            x1*x2^2
            x1^3*x2^4
            x1^3*x2^4
            x1^5*x2^6
            x1^5*x2^6
            x1^7*x2^8
            x1^7*x2^8
            x1^9*x2^10
            x1^9*x2^10
            x1^11*x2^12
            x1^11*x2^12
            x1^13*x2^14

        We can print computational statistics when computing examples::

            sage: var('x1,x2,x3,a0,a2,a3,b0,b1,b3,c0,c1,c2')
            (x1, x2, x3, a0, a2, a3, b0, b1, b3, c0, c1, c2)
            sage: S = LPASeed({x1: a0 + a2*x2 + a3*x3, x2: b0 + b1*x1 + b3*x3, x3: c0 + c1*x1 + c2*x2})
            sage: t = S.mutation_class_iter(verbose=True)
            sage: for seed in t: None
            depth: 0     found: 1
            depth: 1     found: 4
            depth: 2     found: 10
            depth: 3     found: 16
        """

        if algorithm == 'BFS':

            return self._mutation_class_iter_bfs(
                depth=depth, verbose=verbose, return_paths=return_paths)

        elif algorithm == 'DFS':

            return self._mutation_class_iter_dfs(
                depth=depth, verbose=verbose, return_paths=return_paths)

        else:

            raise ValueError('nonsupported search algorithm: %s' % (algorithm))

    def _mutation_class_iter_bfs(self, depth=infinity, verbose=False,
                                 return_paths=False):

        # initialise

        n = self._rank
        seeds_found = [self]
        current_depth = 0
        seeds_to_check = [self]

        # If we are showing depth, show some statistics
        if verbose:
            dc = str(current_depth)
            dc += ' ' * (5 - len(dc))
            nr = str(len(seeds_found))
            nr += ' ' * (10 - len(nr))
            print("depth: %s found: %s" % (dc, nr))

        if return_paths:
            yield (self, [])
        else:
            yield self

        new_seeds_found = True

        while new_seeds_found and current_depth < depth:

            current_depth += 1
            new_seeds = []  # reset new seed list

            for seed in seeds_to_check:

                # we do not need to check the index we last mutated at
                if not seed._mutation_sequence:
                    last_index = None
                else:
                    last_index = seed._mutation_sequence[-1]

                for i in range(n):
                    if i == last_index:
                        continue

                    seed2 = seed.mutate(i, inplace=False)
                    if seed2 not in seeds_found:
                        new_seeds += [seed2]
                        seeds_found += [seed2]
                        if return_paths:
                            yield (seed2, seed2.mutation_sequence())
                        else:
                            yield seed2

            seeds_to_check = new_seeds

            if new_seeds and verbose:
                dc = str(current_depth)
                dc += ' ' * (5 - len(dc))
                nr = str(len(seeds_found))
                nr += ' ' * (10 - len(nr))
                print("depth: %s found: %s" % (dc, nr))

    def _mutation_class_iter_dfs(self, depth=infinity, verbose=False,
                                 return_paths=False):

        # initialise

        n = self._rank
        seeds_found = [self]
        current_depth = 0
        seeds_to_check = [self]

        # If we are showing depth, show some statistics
        if verbose:
            dc = str(current_depth)
            dc += ' ' * (5 - len(dc))
            nr = str(len(seeds_found))
            nr += ' ' * (10 - len(nr))
            print("depth: %s found: %s" % (dc, nr))

        if return_paths:
            yield (self, [])
        else:
            yield self

        new_seeds_found = False

        while seeds_to_check:

            current_depth += 1
            seed = seeds_to_check.pop()

            # are we still within depth constraint?

            if current_depth < depth:

                # we do not need to check the index we last mutated at
                if not seed._mutation_sequence:
                    last_index = None
                else:
                    last_index = seed._mutation_sequence[-1]

                for i in range(n):
                    if i == last_index:
                        continue

                    seed2 = seed.mutate(i, inplace=False)

                    if seed2 not in seeds_found:

                        new_seeds_found = True
                        seeds_found.append(seed2)
                        seeds_to_check.append(seed2)

                        if return_paths:
                            yield (seed2, seed2.mutation_sequence())
                        else:
                            yield seed2

            if not new_seeds_found:

                current_depth -= 1

            elif verbose:
                dc = str(current_depth)
                dc += ' ' * (5 - len(dc))
                nr = str(len(seeds_found))
                nr += ' ' * (10 - len(nr))
                print("depth: %s found: %s" % (dc, nr))

    def mutation_class(self, depth=infinity, verbose=False,
                       return_paths=False, algorithm='BFS'):
        r"""
        Return the mutation class of ``self`` with respect to
        certain constraints.

        .. SEEALSO::

            :meth:`mutation_class_iter`

        INPUT:

        - ``depth`` -- (default: ``infinity``) integer, only seeds with
          distance at most ``depth`` from ``self`` are returned
        - ``verbose`` -- (default: ``False``) if ``True``, the actual depth
          of the mutation is shown
        - ``return_paths`` -- (default: ``False``) if ``True``, a path of
          mutation sequences from ``self`` to the given seed is returned as well
        - ``algorithm`` -- string (default: ``'BFS'``); the search algorithm to
          find new seeds; currently supported options:

          * 'BFS' - breadth-first search
          * 'DFS' - depth-first search

        EXAMPLES:

        .. SEEALSO::

            For further examples see :meth:`mutation_class_iter`.

        We validate the possible sizes of mutation classes in rank two::

            sage: var('x1,x2,A,B,C,D,E,F');
            (x1, x2, A, B, C, D, E, F)
            sage: S = LPASeed({x1: C, x2: C}, coefficients=[C])
            sage: len(S.mutation_class())
            3
            sage: S = LPASeed({x1: A, x2: C + D*x1}, coefficients=[A,C,D])
            sage: len(S.mutation_class())
            4
            sage: S = LPASeed({x1: A + B*x2, x2: C + D*x1}, coefficients=[A,B,C,D])
            sage: len(S.mutation_class())
            5
            sage: S = LPASeed({x1: A + B*x2 + C*x2^2, x2: D + E*x1}, coefficients=[A,B,C,D,E])
            sage: len(S.mutation_class())
            6
            sage: S = LPASeed({x1: A + B*x2 + C*x2^2 + D*x2^3, x2: E + F*x1}, coefficients=[A,B,C,D,E,F])
            sage: len(S.mutation_class())
            8
        """
        return [S for S in self.mutation_class_iter(depth=depth,
                                                    verbose=verbose,
                                                    return_paths=return_paths,
                                                    algorithm=algorithm)]

    def cluster_class_iter(self, depth=infinity, verbose=False,
                           algorithm='BFS'):
        r"""
        Iterator for the cluster class of ``self`` with respect to certain
        constraints.

        .. SEEALSO::

            :meth:`mutation_class_iter`

        INPUT:

        - ``depth`` -- (default: ``infinity``) integer, only clusters with
          distance at most ``depth`` from ``self`` are returned
        - ``verbose`` -- (default: ``False``) if ``True``, the actual depth
          of the mutation is shown
        - ``return_paths`` -- (default: ``False``) if ``True``, a path of
          mutation sequences from ``self`` to the given seed is returned as well
        - ``algorithm`` -- string (default: ``'BFS'``); the search algorithm to
          find new seeds; currently supported options:

          * 'BFS' - breadth-first search
          * 'DFS' - depth-first search

        EXAMPLES:

        We check a classic example::

            sage: var('a,f,C')
            (a, f, C)
            sage: S = LPASeed({a: f + C, f: a + C}, coefficients=[C])
            sage: t = S.cluster_class_iter()
            sage: for cluster in t: print(cluster)
            [a, f]
            [(f + C)/a, f]
            [a, (a + C)/f]
            [(f + C)/a, (a + f + C)/(a*f)]
            [(a + f + C)/(a*f), (a + C)/f]

        .. SEEALSO::

            For further examples see :meth:`mutation_class_iter`.
        """
        mc_iter = self.mutation_class_iter(depth=depth, verbose=verbose,
                                           algorithm=algorithm)
        for c in mc_iter:
            yield c.cluster()

    def cluster_class(self, depth=infinity, verbose=False,
                      algorithm='BFS'):
        r"""
        Return the cluster class of ``self`` with respect to certain
        constraints.

        .. SEEALSO::

            :meth:`mutation_class_iter`

        INPUT:

        - ``depth`` -- (default: ``infinity``) integer, only clusters with
          distance at most ``depth`` from ``self`` are returned
        - ``verbose`` -- (default: ``False``) if ``True``, the actual depth
          of the mutation is shown
        - ``return_paths`` -- (default: ``False``) if ``True``, a path of
          mutation sequences from ``self`` to the given seed is returned as well
        - ``algorithm`` -- string (default: ``'BFS'``); the search algorithm to
          find new seeds; currently supported options:

          * 'BFS' - breadth-first search
          * 'DFS' - depth-first search

        .. SEEALSO::

            For further examples see :meth:`cluster_class_iter`.

        TESTS::

            sage: var('x1, x2, x3')
            (x1, x2, x3)
            sage: LPASeed({x1: 2},base_ring=ZZ).cluster_class()
            [[x1], [2/x1]]
        """
        return [c for c in self.cluster_class_iter(depth=depth,
                                                   verbose=verbose,
                                                   algorithm=algorithm)]

    def variable_class_iter(self, depth=infinity, algorithm='BFS'):
        r"""
        Return an iterator for all cluster variables in the mutation class of
        ``self`` in seeds at most ``depth`` away from ``self``.

        INPUT:

        - ``depth`` -- (default:``infinity``) integer, only seeds with
          distance at most ``depth`` from ``self`` are returned
        - ``algorithm`` -- string (default: ``'BFS'``); the search algorithm to
          find new seeds; currently supported options:

          * 'BFS' - breadth-first search
          * 'DFS' - depth-first search

        EXAMPLES:

        We define a simple iterator for the denominators of seeds in the
        mutation class::

            sage: var('x1, x2, x3, a0, a2, a3, b0, b1, b3, c0, c1, c2')
            (x1, x2, x3, a0, a2, a3, b0, b1, b3, c0, c1, c2)
            sage: S = LPASeed({x1: a0 + a2*x2 + a3*x3, x2: b0 + b1*x1 + b3*x3, x3: c0 + c1*x1 + c2*x2})
            sage: t = S.variable_class_iter()
            sage: for variable in t: print(variable.denominator()) # long time
            1
            1
            1
            x1
            x2
            x3
            x1*x2
            x1*x3
            x2*x3
            x1*x2*x3
        """
        mut_iter = self.mutation_class_iter(depth=depth, verbose=False,
                                            algorithm=algorithm)
        var_class = set()

        for seed in mut_iter:
            for x in seed.cluster():
                if x not in var_class:
                    var_class.add(x)
                    yield x

    def variable_class(self, depth=infinity):
        r"""
        Return all cluster variables in the mutation class of ``self``. These
        are exactly the generators for the LP algebra generated by this seed.

        INPUT:

        - ``depth`` -- (default:``infinity``) integer, only seeds with distance
          at most depth from ``self`` are returned

        .. SEEALSO::

            For further examples see :meth:`variable_class_iter`.

        We find the generators for various LP algebras::

            sage: var('x1,x2,x3')
            (x1, x2, x3)
            sage: S = LPASeed({x1: 1 + x2, x2: 1 + x1})
            sage: S.variable_class()
            [(x1 + x2 + 1)/(x1*x2), (x2 + 1)/x1, (x1 + 1)/x2, x2, x1]
            sage: S = LPASeed({x1: 1 + x2 + x3, x2: 1 + x1 + x3, x3: 1 + x1 + x2})
            sage: S.variable_class()
            [(x1 + x2 + x3 + 1)/(x1*x2*x3), (x2 + x3 + 1)/x1, (x1 + x3 + 1)/x2, (x1 + x2 + 1)/x3, x3, x2, x1]
        """
        var_iter = self.variable_class_iter(depth=depth)
        return sorted(var_iter)

    def is_equivalent(self, other):
        r"""
        Return whether ``self`` and ``other`` are equivalent as LP seeds.
        Two seeds are equivalent if and only if there is a permutation of the
        cluster variables of one seed to get the cluster variables of the other
        seed, up to unit multipliers. Note we also overload equality to
        equivalence.

        INPUT:

        - ``other`` -- ``LPASeed``; the seed which we are comparing ``self`` to

        EXAMPLES:

        Mutating this rank two example five times yields an equivalent seed::

            sage: var('x1,x2')
            (x1, x2)
            sage: S = LPASeed({x1: 1 + x2, x2: 1 + x1})
            sage: S.is_equivalent(S.mutate([0,1,0,1,0], inplace=False))
            True
        """
        if not isinstance(other, LPASeed):
            raise ValueError('%s is not a seed!' % (other))

        if self._rank != other._rank:
            return False

        n = self._rank
        L = []
        for i in range(n):
            for j in range(n):
                x = other.cluster()[j]
                y = self.cluster()[i]
                t = self._ambient_field(x / y)
                if t.numerator().is_unit() and t.denominator().is_unit():
                    L.append(j)
        return len(L) == n

    def exchange_graph(self):
        r"""
        Return the exchange graph of ``self``.

        EXAMPLES:

        We work out the exchange graph for a rank-two example::

            sage: var('x1, x2')
            (x1, x2)
            sage: LPASeed({x1: 1 + x2, x2: 1 + x1}).exchange_graph()
            Graph on 5 vertices

        .. PLOT::

            var('x1, x2')
            G = LPASeed({x1: 1 + x2, x2: 1 + x1}).exchange_graph()
            sphinx_plot(G)

        A rank three example::

            sage: var('x1, x2, x3')
            (x1, x2, x3)
            sage: LPASeed({x1: 1 + x2 + x3, x2: 1 + x1 + x3, x3: 1 + x1 + x2}).exchange_graph()
            Graph on 10 vertices

        .. PLOT::

            var('x1, x2, x3')
            G = LPASeed({x1: 1 + x2 + x3, x2: 1 + x1 + x3, x3: 1 + x1 + x2}).exchange_graph()
            sphinx_plot(G)
        """
        # initialise

        covers = []
        n = self.rank()
        stack = [self]
        known_seeds = []

        # do a depth-first search
        while stack:
            i = stack.pop()
            for k in range(n):
                j = i.mutate(k, inplace=False)
                # need good convenient method of representing seeds on graph
                covers.append((frozenset(i.cluster()), frozenset(j.cluster())))
                if j not in known_seeds:
                    known_seeds += [j]
                    stack.append(j)
        G = Graph(covers)
        G.relabel()
        return G

    def cluster(self):
        r"""
        Return the cluster variables of ``self``.

        EXAMPLES:

        Get the cluster variables after performing a mutation::

            sage: var('x1, x2')
            (x1, x2)
            sage: S = LPASeed({x1: 1 + x2, x2: 1 + x1})
            sage: S.mutate(0)
            A seed with cluster variables [(x2 + 1)/x1, x2] and exchange polynomials [x2 + 1, x1 + 1]
            sage: S.cluster()
            [(x2 + 1)/x1, x2]
        """
        return list(self._ambient_field(x) for x in self._cluster_vars)

    def exchange_polys(self):
        r"""
        Return the exchange polynomials of ``self``.

        EXAMPLES:

        Get the exchange polynomials after performing a mutation::

            sage: var('x1, x2')
            (x1, x2)
            sage: S = LPASeed({x1: 3 + 4*x2, x2: 5 + 6*x1})
            sage: S.mutate(0)
            A seed with cluster variables [(4*x2 + 3)/x1, x2] and exchange polynomials [4*x2 + 3, 5*x1 + 18]
            sage: S.exchange_polys()
            [4*x2 + 3, 5*x1 + 18]
        """
        return list(self._exchange_polys)

    def laurent_polys(self):
        r"""
        Return the exchange Laurent polynomials of ``self``.

        EXAMPLES:

        This seed has non-trivial exchange Laurent polynomials::

            sage: var('x1, x2, x3')
            (x1, x2, x3)
            sage: S = LPASeed({x1: 1 + x2 + x3, x2: 1 + x3, x3: 1 + x2})
            sage: S.laurent_polys()
            [(x2 + x3 + 1)/(x2*x3), x3 + 1, x2 + 1]
        """
        return list(self._ambient_field(f) for f in self._laurent_polys)

    def rank(self):
        r"""
        Return the rank of ``self``. This is the number of
        cluster variables in ``self``.

        EXAMPLES:

        The rank is the number of cluster variables in the seed::

            sage: var('x1,x2,x3,x4')
            (x1, x2, x3, x4)
            sage: LPASeed({x1: 1 + x2, x2: 1 + x1, x3: 1 + x4, x4: 1 + x1}).rank()
            4
        """
        return self._rank

    def randomly_mutate(self, depth, inplace=True):
        r"""
        Randomly mutates this seed at ``depth`` indices. Useful for working out
        if a seed produces a finite type LP algebra.

        INPUT:

        - ``depth`` -- integer, the number of random mutations to perform.

        EXAMPLES:

        Check that randomly mutating a seed keeps it within the mutation class::

            sage: var('x1,x2,x3')
            (x1, x2, x3)
            sage: S = LPASeed({x1: 1 + x2 + x3, x2: 1 + x1 + x3, x3: 1 + x1 + x2})
            sage: T = LPASeed(S)
            sage: S.randomly_mutate(5) # random
            sage: S in T.mutation_class()
            True

        """
        import random
        indices = [random.randint(0, self._rank - 1) for _ in range(depth)]

        self = self.mutate(indices, inplace)

        return self

    def mutation_sequence(self):
        r"""
        Return the list of indices we have mutated ``self`` at.

        EXAMPLES:

        We look at the mutation sequences computed by the mutation class
        iterator::

            sage: var('x1, x2, x3')
            (x1, x2, x3)
            sage: S = LPASeed({x1: 1 + x2 + x3, x2: 1 + x1 + x3, x3: 1 + x1 + x2})
            sage: t = S.mutation_class_iter(algorithm='BFS')
            sage: for seed in t: print(seed.mutation_sequence())
            []
            [0]
            [1]
            [2]
            [0, 1]
            [0, 2]
            [1, 0]
            [1, 2]
            [2, 0]
            [2, 1]
            sage: t = S.mutation_class_iter(algorithm='DFS')
            sage: for seed in t: print(seed.mutation_sequence())
            []
            [0]
            [1]
            [2]
            [2, 0]
            [2, 1]
            [2, 1, 2]
            [2, 1, 2, 0]
            [2, 1, 2, 0, 2]
            [2, 1, 2, 0, 2, 0]
        """
        return _remove_repeat_indices(self._mutation_sequence)

    # equality is currently defined as equivalence (might change later)

    def __eq__(self, other):
        r"""
        Check if ``self`` and ``other`` are equal. This means checking that they
        are equivalent.

        TESTS::

            sage: var('x1')
            x1
            sage: S = LPASeed({x1:2},base_ring=ZZ)
            sage: T = S.mutate(0, inplace=False)
            sage: T==S
            False
        """

        # we check if cluster variables are is_equivalent of each other,
        # up to unit multiple. Since clusters determine the seed,
        # this is sufficient
        return self.is_equivalent(other)

    def _copy_(self):
        r"""
        Return a copy of ``self``.

        TESTS::

            sage: var('x1,x2')
            (x1, x2)
            sage: S = LPASeed({x1: 1 + x2, x2: 1 + x1})
            sage: T1 = S._copy_()
            sage: T1.mutate(0, inplace=True)
            A seed with cluster variables [(x2 + 1)/x1, x2] and exchange polynomials [x2 + 1, x1 + 1]
            sage: T2 = S.mutate(0, inplace=False)
            sage: T1==T2
            True
        """
        return LPASeed(self)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: var('x1, x2')
            (x1, x2)
            sage: print(LPASeed({x1: 1 + x2, x2: 1 + x1}))
            A seed with cluster variables [x1, x2] and exchange polynomials [x2 + 1, x1 + 1]
        """
        return ("A seed with cluster variables {0}"
                " and exchange polynomials {1}").format(self.cluster(),
                                                        self.exchange_polys())

    def __hash__(self):
        """
        Return a hash of ``self``.

        EXAMPLES:

        We check that two seeds with the same hash are equal::

            sage: var('x1,x2')
            (x1, x2)
            sage: S = LPASeed({x1: 1 + x2, x2: 1 + x1})
            sage: T = LPASeed(S)
            sage: S.mutate([0,1,0])
            A seed with cluster variables [(x1 + 1)/x2, (x1 + x2 + 1)/(x1*x2)] and exchange polynomials [x2 + 1, x1 + 1]
            sage: T.mutate([1,0])
            A seed with cluster variables [(x1 + x2 + 1)/(x1*x2), (x1 + 1)/x2] and exchange polynomials [x2 + 1, x1 + 1]
            sage: hash(S) == hash(T)
            True
            sage: S == T
            True
        """
        return hash(frozenset(self.cluster()))

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of this seed.

        TESTS::

            sage: var('x1, x2')
            (x1, x2)
            sage: latex(LPASeed({x1: 1 + x2, x2: 1 + x1}))
            \left(\left(x_{1}, x_{2} + 1\right), \left(x_{2}, x_{1} + 1\right)\right)
        """
        cvars = self.cluster()
        polys = self.exchange_polys()
        cvar_poly_pairs = list()
        for i in range(self.rank()):
            cvar_poly_pairs.append(latex(tuple([cvars[i], polys[i]])))

        return latex(tuple(cvar_poly_pairs))

    def are_laurent_polys_trivial(self):
        r"""
        Return whether this seed has Laurent polynomials with nontrivial denominators.
        """
        for poly in self.laurent_polys():

            if poly.denominator() != 1:

                return False

        return True

    def is_mutation_infinite(self, give_reason=True):
        r"""
        Perform some heuristic checks on the mutation class of this seed to work out if it is mutation infinite or not.

        Completely probabilistic and may false positive; do not use in research setting.
        Warning: This could take a long time for some seeds!
        """

        MAX_DEGREE = 20  # when to decide a seed mutation class is too large based on degrees of exchange polynomials
        counter = 0
        for seed in self.mutation_class_iter():
            counter += 1
            if not seed.passes_rank_two_check():  # type: ignore
                if give_reason:
                    print("Fails rank two check after mutating at indices %s" %
                          (seed.mutation_sequence()))
                return True

            # check degrees of exchange polynomials

            for poly in seed.exchange_polys():  # type: ignore

                for var in poly.variables():

                    if poly.degree(var) > MAX_DEGREE:
                        if give_reason:
                            print("Fails polynomial degree finiteness check after mutating at indices %s" % (
                                seed.mutation_sequence()))
                        return True

        # if we managed to iterate through everything, then it must be finite

        if give_reason:
            print("Mutation finite, has %s seeds in mutation class" % (counter))
        return False

    def passes_rank_two_check(self):
        r"""
        Check that, when restricting to any two variables, that the product of the degrees is strictly less than 4.
        """

        N = self.rank()
        polys = self.exchange_polys()

        for i in range(N):
            for j in range(i, N):
                if polys[i].degree(self._polynomial_ring(self._names[j])) * polys[j].degree(self._polynomial_ring(self._names[i])) >= 4:
                    return False

        return True

    def get_laurent_poly_denominators(self):
        r"""
        Return a list containing the denominators of the Laurent polynomials for this seed.
        """
        return [poly.denominator() for poly in self.laurent_polys()]

    @staticmethod
    def create_generic_seed(vars, polys):

        for i in range(len(polys)):

            coefficients = []
            M = polys[i].monomials()
            new_coeffs = list(var('a_%d%d' % (i, j) for j in range(len(M))))
            coefficients += new_coeffs
            new_poly = 0
            for j in range(len(M)):
                new_poly += new_coeffs[j]*M[j]

            polys[i] = new_poly

        return LPASeed({k: v for k, v in zip(vars, polys)}, coefficients=coefficients)


def _remove_repeat_indices(L):
    r"""
    Remove mutations twice at the same index from the list of indices ``L``.

    This will result in a list with no consecutive entries equal.

    TESTS::

        sage: var('x1, x2, x3, x4, x5')
        (x1, x2, x3, x4, x5)
        sage: S = LPASeed({x1:x2+x3,x2:x1+x4,x3:x1+x2,x4:x1+x3,x5:x1+x4})
        sage: S.mutate([0,1,1,3,2,4,4,2,3])
        A seed with cluster variables [(x2 + x3)/x1, x2, x3, x4, x5] and exchange polynomials [x2 + x3, x1*x4 + x3, x1 + 1, x1*x3 + x2 + x3, x1*x4 + x2 + x3]
        sage: S.mutation_sequence()
        [0]
    """
    G = []
    flag = False  # we raise the flag iff an adjacent duplicate exists
    index = 0
    while index < len(L):
        if index == len(L) - 1 or L[index] != L[index + 1]:
            G.append(L[index])
            index += 1
        else:
            flag = True
            index += 2  # skip over next element

    if flag:
        return _remove_repeat_indices(G)
    else:
        return G
