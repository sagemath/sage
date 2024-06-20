r"""
Implementation of flag algebras, with a class for combinatorial theories


A combinatorial theory is any theory with universal axioms only, 
(therefore the elements satisfy a heredetary property). This 
implementation allows the construction of any such theory, and 
can perform flag algebraic computations on them. The theory of
flag algebras is from [Raz2007]_

To find out more about flags, how to create and manipulate them,
see :mod:`sage.algebras.flag`. This docstring is for combinatorial
theories and combinatorial optimization problems using flag algebras.

The examples will use ::

    sage: from sage.algebras.flag_algebras import *

To create a `CombinatorialTheory` object we need to first know the signature.
For example graphs have one relational symbol for the edges, of arity 2.
We could provide for example `edges=2` as a parameter to the constructor, 
showing this fact.

In addition we need two important functions:

A generator function, showing the program how to generate elements of 
the theory for a fixed size. For example, to graphs with ordered vertices (
or equivalently 0-1 symmetric matrices with 0 diagonal) can be generated
with the following code ::
    
    sage: def test_generator_ov_graph(n):
    ....:    full = list(itertools.combinations(range(n), int(2)))
    ....:    for ii in range(binomial(n, 2)+1):
    ....:        for xx in itertools.combinations(full, int(ii)):
    ....:            yield {'edges': xx}
    

This function takes `n` integer as input and returns the possible
ways one can construct such graphs with ordered vertices. Note
the elements returned are dictionaries. The key is `'edges'` and
the value is the actual list of edges. `'edges'` key is important to
match the one defined in the signature.

The second important function is the identifier function. This takes
as input an element of the theory, with some of the points marked,
and returns something that uniquely identifies this element, up to
automorphism. The automorphisms must respect the marked points. For
ordered vertex graphs, since there is no automorphism of the vertices
giving the exact same order this is easy, we can just return the edges
sorted with the marked points and the total number of points. ::
    
    sage: def test_identify_ov_graph(n, ftype_points, edges):
    ....:    return (n, tuple(ftype_points), \
    ....:    tuple(sorted(list(edges))))
    ....: 

For theories where members can have more automorphisms, this might be
harder. Then to create the theory, provide a name, the two functions and
the signature ::

    sage: TestOVGraphTheory = CombinatorialTheory('TestOVGraph', \
    ....: test_generator_ov_graph, test_identify_ov_graph, edges=2)

The following common theories are already implemented:
-GraphTheory
-ThreeGraphTheory
-DiGraphTheory
-TournamentTheory
-PermutationTheory
-OVGraphTheory (graphs with ordered vertices)
-OEGraphTheory (graphs with ordered edges)
-RamseyGraphTheory (see [LiPf2021]_ for explanation)

The rest of this docstring will use `GraphTheory` since the number 
of structures is realtively small there. `Flag` docstring shows
ways to create and calculate with flags. To create a flag we need to
provide the vertex size and a list of elements for each signature. 
To create an edge flag for graphs, use ::

    sage: e = GraphTheory(2, edges=[[0, 1]])

To create a triangle `K_3` we can write ::
    
    sage: k3 = GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]])

If we have a theory, say GraphTheory, we can simply exclude structures
with :func:`exclude`. This takes a list of flags or a single flag :: 

    sage: GraphTheory.exclude(k3)

Excluding structures overwrites the theory, and in the future will only
consider members without the excluded structures. So with the above, 
excluding a triangle will make GraphTheory not generate `k3`, or any larger
flag with induced `k3` in it. We can check this by generating flags of size
`4` ::

    sage: GraphTheory.generate_flags(4)
    (Flag on 4 points, ftype from [] with edges=[],
     Flag on 4 points, ftype from [] with edges=[[0, 3]],
     Flag on 4 points, ftype from [] with edges=[[0, 3], [1, 3]],
     Flag on 4 points, ftype from [] with edges=[[0, 3], [1, 3], [2, 3]],
     Flag on 4 points, ftype from [] with edges=[[0, 2], [1, 3]],
     Flag on 4 points, ftype from [] with edges=[[0, 2], [0, 3], [1, 3]],
     Flag on 4 points, ftype from [] with edges=[[0, 2], [0, 3], [1, 2], [1, 3]])

Excluding structures overwrites the previously excluded structures. Calling 
`exclude` without any arguments excludes nothing, giving back the original
theory.

To optimize the density of a flag (or linear combination of flags) in a theory,
we can call :func:`optimize_problem`. To try to find the maximum number of 
edges `e` in `k3` free graphs we can write ::
    
    sage: x = GraphTheory.optimize_problem(e, 3)
    ...
    Optimal solution found.
    sage: abs(x-0.5)<1e-6
    True
    
The second parameter, `optimize_problem(e, 3)` indicates the maximum size the
program expands the flags. We can reset the excluded graphs, and try to minimize 
the density of triangles and empty triples ::
    
    sage: GraphTheory.exclude()
    sage: e3 = GraphTheory(3)
    sage: x = GraphTheory.optimize_problem(e3+k3, 3, maximize=False)
    ...
    sage: abs(x-0.25)<1e-6
    True

The `solve_problem` function requires csdpy, an sdp solver. But it is possible
to get these relations directly, by expressing the inequalities 
as a sum of squares ::

    sage: GraphTheory.exclude(k3)
    sage: pe = GraphTheory(2, ftype=[0]) - 1/2
    sage: 1/2 >= pe.mul_project(pe) * 2 + e
    True

:func:`mul_project` is short for multiplication and projection, the 
non-negativity of self multiplication is preserved after projection so
`pe.mul_project(pe)` is non-negative on all large enough structures 
giving that `1/2` is larger than `e` in all large enough structures.

The following longer example shows that the density of K^3_4 is always
less than 3/8 in K^3_5-free hypergraphs. It uses the ThreeGraphTheory
object to deal with 3-uniform hypergraphs and hand-picked squares.
These values come from [Bod2023]_::
    
    sage: em5 = ThreeGraphTheory(5, edges=[])
    sage: ThreeGraphTheory.exclude(em5)
    
    sage: em4 = ThreeGraphTheory(4, edges=[])
    
    sage: em3p1 = ThreeGraphTheory(3, ftype_points=[0], edges=[])
    sage: sq1 = (em3p1 - 3/4).mul_project(em3p1 - 3/4) # todo: not implemented
    
    sage: la4p2 = ThreeGraphTheory(4, ftype_points=[0, 1], edges=[[0, 2, 3]])
    sage: lb4p2 = ThreeGraphTheory(4, ftype_points=[0, 1], edges=[[1, 2, 3]])
    sage: sq2 = (la4p2 - lb4p2).mul_project(la4p2 - lb4p2) # todo: not implemented
    
    sage: ma4p3 = ThreeGraphTheory(4, ftype_points=[0, 1, 2], edges=[[0, 1, 3]])
    sage: mb4p3 = ThreeGraphTheory(4, ftype_points=[0, 1, 2], edges=[[0, 2, 3]])
    sage: mc4p3 = ThreeGraphTheory(4, ftype_points=[0, 1, 2], edges=[[1, 2, 3]])
    sage: sq3 = (ma4p3 + mb4p3 + mc4p3 - 1/2).mul_project(ma4p3 + mb4p3 + mc4p3 - 1/2) # todo: not implemented
     
    sage: em4p3 = ThreeGraphTheory(4, ftype_points=[0, 1, 2])
    sage: sq4 = (em4p3 - 1/2).mul_project(em4p3 - 1/2) # todo: not implemented
    
    sage: n5q4 = ThreeGraphTheory(5, ftype_points=[0, 1, 2, 3], edges=[[0, 1, 2]])
    sage: sq5 = (n5q4 - 1/2).mul_project(n5q4 - 1/2) # todo: not implemented
    
    sage: oa5t4 = ThreeGraphTheory(5, ftype_points=[0, 1, 2, 3], edges=[[0, 1, 4]])
    sage: ob5t4 = ThreeGraphTheory(5, ftype_points=[0, 1, 2, 3], edges=[[2, 3, 4]])
    sage: sq6 = (oa5t4 - ob5t4).mul_project(oa5t4 - ob5t4) # todo: not implemented
     
    sage: sos = 2/3 * sq1 + 1/6 * sq2 + 13/12 * sq3 + 11/12 * sq4 + 2 * sq5 + 1/2 * sq6 # todo: not implemented
    sage: 3/8 >= sos + em4 # todo: not implemented
    True


.. SEEALSO::
    :func:`CombinatorialTheory.__init__`
    :func:`CombinatorialTheory.exclude`
    :func:`CombinatorialTheory.optimize_problem`
    :func:`CombinatorialTheory.generate_flags`

AUTHORS:

- Levente Bodnar (Dec 2023): Initial version

"""

# ****************************************************************************
#       Copyright (C) 2023 LEVENTE BODNAR <bodnalev at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from functools import lru_cache
import itertools

from sage.structure.richcmp import richcmp
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import CommutativeAlgebraElement
from sage.rings.ring import CommutativeAlgebra
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.algebras.flag import Flag

from sage.categories.sets_cat import Sets

from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix

from sage.misc.prandom import randint
from sage.arith.misc import falling_factorial, binomial

from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph
from sage.misc.lazy_import import lazy_import
lazy_import("sage.graphs.graph_generators", "graphs")
lazy_import("sage.graphs.digraph_generators", "digraphs")
lazy_import("sage.graphs.hypergraph_generators", "hypergraphs")

import pickle
import os

class CombinatorialTheory(Parent, UniqueRepresentation):
    
    Element = Flag
    
    def __init__(self, name, generator, identifier, **signature):
        r"""
        Initialize a Combinatorial Theory
        
        A combinatorial theory is any theory with universal axioms only, 
        (therefore the elements satisfy a heredetary property).
        See the file docstring for more information.

        INPUT:

        - ``name`` -- string; name of the Theory
        - ``generator`` -- function; generates elements 
            of the theory. For a given input ``n`` 
            returns a list of elements of the theory
            in a dictionary format (for each
            value in the signature, one dictionary 
            entry describing the blocks corresponding to
            that signature)
        - ``identifier`` -- function; given a structure
            with the matching signature of this theory,
            returns a unique identifier, such that
            automorphic structures return the same 
            value.
        - ``**signature`` -- named integers; the signature
            of the theory, for each name a corresponding number
            giving the arity of that symbol

        OUTPUT: A CombinatorialTheory object

        EXAMPLES::

        This example shows how to create the theory for graphs 
        with ordered vertices (or equivalently 0-1 matrices)::
            
            sage: from sage.algebras.flag_algebras import *
            sage: def test_generator_ov_graph(n):
            ....:    full = list(itertools.combinations(range(n), int(2)))
            ....:    for ii in range(binomial(n, 2)+1):
            ....:        for xx in itertools.combinations(full, int(ii)):
            ....:            yield {'edges': xx}
            ....: 
            sage: def test_identify_ov_graph(n, ftype_points, edges):
            ....:    return (n, tuple(ftype_points), \
            ....:    tuple(sorted(list(edges))))
            ....: 
            sage: TestOVGraphTheory = CombinatorialTheory('TestOVGraph', \
            ....: test_generator_ov_graph, test_identify_ov_graph, edges=2)
            sage: TestOVGraphTheory
            Theory for TestOVGraph

        .. NOTE::

            There are pre-constructed CombinatorialTheory objects
            in sage.algebras.flag_algebras for the following:
            -GraphTheory
            -ThreeGraphTheory
            -DiGraphTheory
            -TournamentTheory
            -PermutationTheory
            -OVGraphTheory (graphs with ordered vertices)
            -OEGraphTheory (graphs with ordered edges)
            -RamseyGraphTheory (see [LiPf2021]_ for explanation)
        """
        self._signature = signature
        self._excluded = []
        self._generator = generator
        self._identifier = identifier
        self._name = name
        Parent.__init__(self, category=(Sets(), ))
        self._populate_coercion_lists_()
    
    def _compress(self, numbers):
        table = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_"
        
        p = 1007016245527451
        rem = 0
        for ii in numbers:
            rem = (rem*16 + ii) % p
        ret = ""
        if rem==0:
            ret = table[0]
        while rem > 0:
            ret += table[rem%64]
            rem //= 64
        return ret
    
    def clear(self):
        self._identify.cache_clear()
        ns = "calcs/" + self._name + "."
        for xx in os.listdir("calcs/"):
            if xx.startswith(ns):
                os.remove("calcs/"+xx)
    
    def _save(self, ind, ret, is_table):
        ns = "calcs/" + self._name + "."
        if is_table:
            excluded, n1, n2, large_ftype, ftype_inj = ind
            numsind = [0]
            for xx in excluded:
                numsind += xx.raw_numbers()
            numsind += [n1, n2] + large_ftype.raw_numbers() + list(ftype_inj)
        else:
            excluded, n, ftype = ind
            numsind = [1]
            for xx in excluded:
                numsind += xx.raw_numbers()
            numsind += [n] + ftype.raw_numbers()
        save_name = ns + self._compress(numsind)
        os.makedirs(os.path.dirname(save_name), exist_ok=True)
        with open(save_name, 'wb') as file:
            pickle.dump(ret, file)
    
    def _load(self, ind, is_table):
        ns = "calcs/" + self._name + "."
        
        if is_table:
            excluded, n1, n2, large_ftype, ftype_inj = ind
            numsind = [0]
            for xx in excluded:
                numsind += xx.raw_numbers()
            numsind += [n1, n2] + large_ftype.raw_numbers() + list(ftype_inj)
        else:
            excluded, n, ftype = ind
            numsind = [1]
            for xx in excluded:
                numsind += xx.raw_numbers()
            numsind += [n] + ftype.raw_numbers()
        load_name = ns + self._compress(numsind)
        
        if os.path.isfile(load_name):
            with open(load_name, 'rb') as file:
                ret = pickle.load(file)
                if is_table:
                    return ret
                else:
                    for xx in ret:
                        xx._set_parent(self)
                    return ret
        return None
    
    def _repr_(self):
        r"""
        Give a nice string representation of the theory object

        OUTPUT: The representation as a string

        EXAMPLES::

            sage: from sage.algebras.flag_algebras import *
            sage: print(GraphTheory)
            Theory for Graph
        """
        return 'Theory for {}'.format(self._name)
    
    def _element_constructor_(self, n, **kwds):
        r"""
        Construct elements of this theory

        INPUT:

        - ``n`` -- integer; number of points of the flag
        - ``**kwds`` -- can contain ftype_points, listing
            the points that will form part of the ftype;
            and can contain the blocks for each signature.
            If they are not included, they are assumed to 
            be empty lists.

        OUTPUT: A Flag with the given parameters

        EXAMPLES::

        Create an empty graph on 3 vertices ::

            sage: from sage.algebras.flag_algebras import *
            sage: GraphTheory(3)
            Flag on 3 points, ftype from [] with edges=[]
        
        Create an edge with one point marked as an ftype ::
        
            sage: GraphTheory(2, ftype_points=[0], edges=[[0, 1]])
            Flag on 2 points, ftype from [0] with edges=[[0, 1]]
            
        Create a RamseyGraphTheory flag, a fully colored
        triangle (useful for calculating R(K_3), see 
        :func:`solve_problem`) ::
            
            sage: RamseyGraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]])
            Flag on 3 points, ftype from [] with edges=[[0, 1], [0, 2], [1, 2]], edges_marked=[]

        .. NOTE::

            Different input parameters can result in equal objects, for 
            example the following two graphs are automorphic::
            sage: b1 = [[0, 1], [0, 2], [0, 4], [1, 3], [2, 4]]
            sage: b2 = [[0, 4], [1, 2], [1, 3], [2, 3], [3, 4]]
            sage: g1 = GraphTheory(5, edges=b1)
            sage: g2 = GraphTheory(5, edges=b2)
            sage: g1==g2
            True

        .. SEEALSO::

            :func:`__init__` of :class:`Flag`
        """
        return self.element_class(self, n, **kwds)
    
    def empty_element(self):
        r"""
        Returns the empty element, ``n``=0 and no blocks

        OUTPUT: The empty element of the CombinatorialTheory

        EXAMPLES::

            sage: from sage.algebras.flag_algebras import *
            sage: GraphTheory.empty_element()
            Ftype on 0 points with edges=[]

        .. NOTE::

            Since the underlying vertex set (empty set)
            is the same as the ftype point set, this is
            an ftype

        .. SEEALSO::

            :func:`empty`
        """
        return self.element_class(self, 0)
    
    empty = empty_element
    
    def _an_element_(self, n=0, ftype=None):
        r"""
        Returns a random element

        INPUT:

        - ``n`` -- integer (default: `0`); size of the element
        - ``ftype`` -- Flag (default: `None`); ftype of the element
            if not provided then returns an element with empty ftype

        OUTPUT: A Flag with matching parameters
        """
        if ftype==None:
            ftype = self.empty_element()
        if n==None or n==ftype.size():
            return ftype
        ls = self.generate_flags(n, ftype)
        return ls[randint(0, len(ls)-1)]
    
    def some_elements(self):
        r"""
        Returns a list of elements
        """
        pt = self.element_class(self, 1, ftype=[0])
        return [self._an_element_(), 
                self._an_element_(n=2), 
                self._an_element_(ftype=pt), 
                self._an_element_(n=2, ftype=pt)]
    
    def signature(self):
        return self._signature
    
    def identify(self, n, ftype_points, **blocks):
        r"""
        The function used to test for equality.

        INPUT:

        - ``n`` -- integer; size of the flag
        - ``ftype_points`` -- list; the points of the ftype
        - ``**blocks`` -- the blocks for each signature

        OUTPUT: The identifier of the structure defined by the
            ``identifier`` function in the __init__

        .. SEEALSO::

            :func:`Flag.unique`
        """
        blocks = {key:tuple([tuple(xx) for xx in blocks[key]]) 
                  for key in blocks}
        return self._identify(n, tuple(ftype_points), **blocks)
    
    @lru_cache(maxsize=None)
    def _identify(self, n, ftype_points, **blocks):
        r"""
        The hidden _identify, the inputs are in a tuple form
        and is cached for some speed
        """
        return self._identifier(n, ftype_points, **blocks)
    
    def exclude(self, flags=[]):
        r"""
        Exclude some induced flags from the theory
        
        This allows creation of CombinatorialTheory -s with excluded
        flags. The flags are not allowed to appear as an induced
        substructure in any of the generated flags later.

        INPUT:

        - ``flags`` -- list of flags or a flag (default: `[]`); 
            The list of flags to exclude, flags are treated as
            a singleton list

        EXAMPLES::

        How to create triangle-free graphs ::

            sage: from sage.algebras.flag_algebras import *
            sage: triangle = GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]])
            sage: GraphTheory.exclude(triangle)
        
        There are 14 graphs on 5 vertices without triangles ::
        
            sage: len(GraphTheory.generate_flags(5))
            14

        .. NOTE::

            Calling :func:`exclude` again will overwrite the list
            of excluded structures. So calling exclude() again, gives
            back the original theory

        TESTS::

            sage: from sage.algebras.flag_algebras import *
            sage: ThreeGraphTheory.exclude(ThreeGraphTheory(4))
            sage: len(ThreeGraphTheory.generate_flags(5))
            23
            sage: TournamentTheory.exclude(TournamentTheory(3, edges=[[0, 1], [1, 2], [2, 0]]))
            sage: TournamentTheory.generate_flags(5)
            (Flag on 5 points, ftype from [] with edges=[[1, 0], [2, 0], [2, 1], [3, 0], [3, 1], [3, 2], [4, 0], [4, 1], [4, 2], [4, 3]],)
            
        """
        if type(flags)==Flag:
            self._excluded = [flags]
        else:
            self._excluded = flags
    
    def _check_excluded(self, elms):
        r"""
        Helper to check the excluded structures in generation
        """
        flg = elms[0]
        for xx in elms[1]:
            if xx<= flg:
                return False
        return True
    
    def optimize_problem(self, target_element, target_size, \
                         ftypes=None, maximize=True, certificate=False, \
                         positives=None):
        r"""
        Try to maximize or minimize the value of `target_element`
        
        The algorithm calculates the multiplication tables and 
        sends the SDP problem to cvxopt.
        
        INPUT:

        - ``target_element`` -- Flag or FlagAlgebraElement; 
            the target whose density this function tries to
            maximize or minimize in large structures.
        - ``target_size`` -- integer; The program evaluates
            flags and the relations between them up to this
            size.
        - ``ftypes`` -- list of ftypes (default: `None`); 
            the list of ftypes the program will consider,
            if not provided (or `None`) then generates all
            possible ftypes.
        - ``maximize`` -- boolean (default: `True`); 
        - ``certificate`` -- boolean (default: `False`);
        - ``positives`` -- list of flag algebra elements, 
            optimizer will assume those are positive, can
            have different types

        OUTPUT: A bound for the optimization problem. If 
            certificate is requested then returns the entire
            output of the solver as the second argument.

        EXAMPLES::
        
        
        Mantel's theorem, calculate the maximum density
        of edges in triangle-free graphs ::

            sage: from sage.algebras.flag_algebras import *
            sage: GraphTheory.exclude(GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]]))
            sage: x = GraphTheory.optimize_problem(GraphTheory(2, edges=[[0, 1]]), 3)
            ...
            Optimal solution found.
            sage: abs(x-0.5)<1e-6
            True
        
        Generalized Turan problem, for example maximum density of K_3
        in K_5 -free graphs. The complement is calculated, optimizing empty
        E_3 in E_5 -free graphs. They are equivalent

            sage: GraphTheory.exclude(GraphTheory(5))
            sage: x = GraphTheory.optimize_problem(GraphTheory(3), 5) # long time (5 second)
            ...
            Optimal solution found.
        
        
        Generalized Turan problem for hypergraphs, maximum density of K^3_4
        in K^3_5 -free 3-graphs. Again, the complement is calculated. This is 
        long and should not be normally tested

            sage: ThreeGraphTheory.exclude(ThreeGraphTheory(5))
            sage: ThreeGraphTheory.optimize_problem(ThreeGraphTheory(4), 6) # todo: not implemented
            ...
            Optimal solution found.
        
        
        The minimum number of transitive tournaments is attained at 
        a random tournament [CoRa2015]_::

            sage: x = TournamentTheory.optimize_problem( \
            ....: TournamentTheory(3, edges=[[0, 1], [0, 2], [1, 2]]), \
            ....: 3, maximize=False)
            ...
            Optimal solution found.
        
        
        Ramsey's theorem, the K_5 is the largest 2-colorable complete
        graph without monochromatic K_3. (see [LiPf2021]_) ::

            sage: RamseyGraphTheory.exclude(RamseyGraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]]))
            sage: x = RamseyGraphTheory.optimize_problem(RamseyGraphTheory(2), 4, maximize=False)
            ...
            Optimal solution found.
            sage: abs(x-0.2)<1e-6
            True
        
        .. NOTE::

            If `target_size` is too large, the calculation might take a 
            really long time. On a personal computer the following 
            maximum target_size values are recommended:
            -GraphTheory: 8
            -ThreeGraphTheory: 6
            -DiGraphTheory: 5
            -TournamentTheory: 9
            -PermutationTheory: 8
            -RamseyGraphTheory: 8
            -OVGraphTheory: 5
            -OEGraphTheory: 4
        """
        import time
        
        from csdpy import solve_sdp
        import numpy as np
        from tqdm import tqdm
        import sys
        
        current = time.time()
        #calculate constraints from positive vectors
        if positives == None:
            constraints_flags = []
            constraints_vals = []
        else:
            constraints_flags = []
            for ii in range(len(positives)):
                fv = positives[ii]
                if isinstance(fv, Flag):
                    continue
                d = target_size - fv.size()
                k = fv.ftype().size()
                terms = fv.afae().parent().generate_flags(k+d)
                constraints_flags += [fv.mul_project(xx) for xx in terms]
            constraints_vals = [0]*len(constraints_flags)
        
        #calculate ftypes
        if ftypes is None:
            flags = [flag for kk in range(2-target_size%2, target_size-1, 2) 
                      for flag in self.generate_flags(kk)]
            ftypes = [flag.subflag([], ftype_points=list(range(flag.size()))) \
                      for flag in flags]

        print("Ftypes constructed in {:.2f}s".format(time.time() - current), flush=True); current = time.time()
        block_sizes = [len(self.generate_flags((target_size + \
                       ftype.size())//2, ftype)) for ftype in ftypes]
        constraints = len(self.generate_flags(target_size))

        alg = FlagAlgebra(QQ, self)
        
        one_vector = alg(target_size, [1]*len(self.generate_flags(target_size)))
        constraints_flags.extend([one_vector, one_vector*(-1)])
        constraints_vals.extend([1, -1])

        block_sizes.extend([-constraints, -2*len(constraints_vals)])
        block_num = len(block_sizes)
        mat_inds = []
        mat_vals = []
        print("Block sizes done in {:.2f}s".format(time.time() - current), flush=True); current = time.time()
        print("Block sizes are {}".format(block_sizes), flush=True)
        print("Calculating product matrices for {} ftypes and {} structures".format(len(ftypes), constraints), flush=True)
        for ii, ftype in (pbar := tqdm(enumerate(ftypes), file=sys.stdout)):
            ns = (target_size + ftype.size())//2
            fls = self.generate_flags(ns, ftype)
            table = self.mul_project_table(ns, ns, ftype, [])
            for gg, mm in enumerate(table):
                dd = mm._dict()
                if len(dd)>0:
                    inds, values = zip(*mm._dict().items())
                    iinds, jinds = zip(*inds)
                    for cc in range(len(iinds)):
                        if iinds[cc]>=jinds[cc]:
                            mat_inds.extend([gg+1, ii+1, iinds[cc]+1, 
                                             jinds[cc]+1])
                            mat_vals.append(values[cc])
            pbar.set_description("{} is complete".format(ftype))
        
        print("Table calculation done in {:.2f}s".format(time.time() - current), flush=True); current = time.time()
        if maximize:
            avals = (target_element*(-1)<<(target_size - \
                                           target_element.size())).values()
        else:
            avals = (target_element<<(target_size - \
                                      target_element.size())).values()

        for ii in range(len(constraints_vals)):
            mat_inds.extend([0, block_num, 1+ii, 1+ii])
            mat_vals.append(constraints_vals[ii])
        
        constraints_flags_vec = [(xx<<(target_size-xx.size())).values() for xx in constraints_flags]
        for gg in range(constraints):
            mat_inds.extend([gg+1, block_num-1, gg+1, gg+1])
            mat_vals.append(1)
            for ii in range(len(constraints_flags_vec)):
                mat_inds.extend([gg+1, block_num, ii+1, ii+1])
                mat_vals.append(constraints_flags_vec[ii][gg])
        print("Target and constraint calculation done in {:.2f}s\n".format(time.time() - current), flush=True); current = time.time()
        
        sdp_result = solve_sdp(block_sizes, list(avals), 
                               mat_inds, mat_vals)
        if maximize:
            ret = max(-sdp_result['primal'], -sdp_result['dual'])
        else:
            ret = min(sdp_result['primal'], sdp_result['dual'])
        print("Result is {}".format(ret), flush=True)
        if certificate:
            ralg = FlagAlgebra(RR, self)
            vec = ralg(target_size, sdp_result['y'])
            ret = (ret, vec, sdp_result)
        return ret
    
    def _gfe(self, excluded, n, ftype):
        r"""
        Cached version of generate flags excluded

        .. SEEALSO::

            :func:`generate_flags`
        """
        if ftype==None:
            ftype = self.empty()
        ind = (excluded, n, ftype)
        loaded = self._load(ind, False)
        if loaded != None:
            return loaded
        
        import multiprocessing as mp
        
        if ftype.size()==0: #just generate empty elements
            if n==0:
                ret = (self.empty_element(), )
                self._save(ind, ret, False)
                return ret
            if len(excluded)==0: #just return the output of the generator
                ret = tuple([self.element_class(self, n, **xx) for xx in self._generator(n)])
                self._save(ind, ret, False)
                return ret
            
            #otherwise check each generated for the excluded values
            slist = [(xx, excluded) for xx in self._gfe(tuple(), n, None)]
            pool = mp.Pool(mp.cpu_count()-1)
            canincl = pool.map(self._check_excluded, slist)
            pool.close(); pool.join()
            ret = tuple([slist[ii][0] for ii in range(len(slist)) if canincl[ii]])
            self._save(ind, ret, False)
            return ret
        
        #generate flags by first getting the empty structures then finding the flags
        empstrs = self._gfe(excluded, n, None)
        pool = mp.Pool(mp.cpu_count()-1)
        pares = pool.map(ftype._ftypes_inside, empstrs)
        pool.close(); pool.join()
        ret = []
        for coll in pares:
            for xx in coll:
                if xx not in ret:
                    ret.append(xx)
        ret = tuple(ret)
        self._save(ind, ret, False)
        return ret
    
    def generate_flags(self, n, ftype=None):
        r"""
        Returns the list of flags with a given size and ftype

        INPUT:

        - ``n`` -- integer; the size of the returned structures
        - ``ftype`` -- Flag; the ftype of the returned structures

        OUTPUT: List of all flags with given size and ftype

        EXAMPLES::

        There are 4 graphs on 3 vertices. Flags with empty
        ftype correspond to elements of the theory ::

            sage: from sage.algebras.flag_algebras import *
            sage: len(GraphTheory.generate_flags(3))
            4
        
        There are 6 graph flags with one vertex ftype. The
        "cherry" ([[0, 1], [0, 2]]) and the complement can be 
        marked two different ways to a flag
        
            sage: len(GraphTheory.generate_flags(3, GraphTheory(1, ftype_points=[0])))
            6
        
        .. NOTE::

            See the notes on :func:`optimize_problem`. A large `n` can
            result in large number of structures.
        """
        if ftype==None:
            ftype = self.empty()
        else:
            if not ftype.is_ftype():
                raise ValueError('{} is not an Ftype'.format(ftype))
            if ftype not in self:
                raise ValueError('{} is not a part of this theory'.format(ftype))
        ftype_size = ftype.size()
        if n<ftype_size:
            return tuple()
        elif n==ftype_size:
            return (ftype, )
        return self._gfe(tuple(self._excluded), n, ftype)
    
    def mul_project_table(self, n1, n2, large_ftype, ftype_inj=None):
        r"""
        Returns the multiplication projection table
        
        `ftype_inj` specifies a projection, an embedding of a
        smaller ftype inside the `large_ftype`. Two flag vectors
        with `n1` and `n2` underlying flag size set and `large_ftype`
        ftype can be multiplied together and projected with this table.
        The result is a tuple of sparse matrices corresponding to the 
        coefficients in the list of flags with the projected ftype

        INPUT:

        - ``n1`` -- integer; Vertex size of first flag vector
        - ``n2`` -- integer; Vertex size of second flag vector
        - ``large_ftype`` -- Flag; The ftype of the flag vectors
        - ``ftype_inj`` -- list (default: `None`); The injection
            from the smaller ftype to the `large_ftype`. If `None` 
            then the calculation is performed without any projection
            (so `ftype_inj` is the identity bijection)

        OUTPUT: A tuple of sparse matrices

        .. NOTE::

            This just transforms the input to a standard form
            and then calls :func:`_mpte` which is cached for speed. 
            This table is used in all sorts of operations, flag algebra
            multiplication and projection. But the main use happens in 
            optimize_problem, where these tables form the semidefinite
            constraint.

        .. SEEALSO::

            :func:`_mpte`
            :func:`optimize_problem`
            :func:`FlagAlgebraElement.mul_project`
            :func:`FlagAlgebraElement._mul_`
            :func:`FlagAlgebraElement.project`

        TESTS::

            sage: from sage.algebras.flag_algebras import *
            sage: table = GraphTheory.mul_project_table(2, 2, GraphTheory(1, ftype_points=[0]), [])
            sage: table[1][0, 0]
            1/3
            
            sage: table = RamseyGraphTheory.mul_project_table(3, 3, RamseyGraphTheory(2, ftype_points=[0, 1]), [])
            sage: table[3][1, 1]
            1/6
        """
        large_size = large_ftype.size()
        if ftype_inj==None:
            ftype_inj = tuple(range(large_size))
        else:
            ftype_inj = tuple(ftype_inj)
            if len([ii for ii in ftype_inj if (ii not in range(large_size))])!=0:
                raise ValueError('ftype_inj must map into the points of {}'.format(large_ftype))
            if len(set(ftype_inj)) != len(ftype_inj):
                raise ValueError('ftype_inj must be injective (no repeated elements)')
        return self._mpte(tuple(self._excluded), n1, n2, large_ftype, ftype_inj)
    
    def _mpte(self, excluded, n1, n2, large_ftype, ftype_inj):
        r"""
        The (hidden) cached version of :func:`mul_project_table`
        """
        ind = (excluded, n1, n2, large_ftype, ftype_inj)
        loaded = self._load(ind, True)
        if loaded != None:
            return loaded
        
        from sage.matrix.args import MatrixArgs
        import multiprocessing as mp
        ftype_inj = list(ftype_inj)
        large_size = large_ftype.size()
        small_ftype = large_ftype.subflag([], ftype_points=ftype_inj)
        small_size = small_ftype.size()
        ftype_remap = ftype_inj + [ii for ii in range(large_size) if (ii not in ftype_inj)]
        N = n1 + n2 - large_size
        
        Nflgs = self._gfe(excluded, N, small_ftype)
        n1flgs = self._gfe(excluded, n1, large_ftype)
        n2flgs = self._gfe(excluded, n2, large_ftype)
        
        slist = tuple((flg, n1, n1flgs, n2, n2flgs, ftype_remap, large_ftype, small_ftype) for flg in Nflgs)
        
        pool = mp.Pool(mp.cpu_count()-1)
        mats = pool.map(self._density_wrapper, slist)
        pool.close(); pool.join()
        
        norm = falling_factorial(N - small_size, large_size - small_size) 
        norm *= binomial(N - large_size, n1 - large_size)
        
        ret = tuple([MatrixArgs(QQ, mat[0], mat[1], entries=mat[2]).matrix()/norm for mat in mats])
        
        self._save(ind, ret, True)
        return ret
    
    def _density_wrapper(self, ar):
        r"""
        Helper function used in the parallelization of calculating densities
        """
        return ar[0].densities(ar[1], ar[2], ar[3], ar[4], ar[5], ar[6], ar[7])

class FlagAlgebraElement(CommutativeAlgebraElement):
    def __init__(self, parent, n, values):
        r"""
        Initialize a Flag Algebra Element
        
        INPUT:

        - ``parent`` -- FlagAlgebra; The parent :class:`FlagAlgebra`
        - ``n`` -- integer; the size of the flags
        - ``values`` -- vector; the values representing the
            linear combination of :class:`Flag` elements

        OUTPUT: A FlagAlgebraElement object with the given initial parameters

        .. NOTE::

            It is recommended to use the FlagAlgebra element constructor

        .. SEEALSO::

            :func:`FlagAlgebra._element_constructor_`
        """
        self._flags = parent.generate_flags(n)
        if len(values)!=len(self._flags):
            raise ValueError('The coefficients must have the same length as the number of flags')
        self._n = n
        base = parent.base()
        self._values = vector(base, values)
        self._ftype = parent.ftype()
        CommutativeAlgebraElement.__init__(self, parent)
    
    def ftype(self):
        r"""
        Return the ftype of the parent FlagAlgebra
        
        The algebra is defined over elements of a CombinatorialTheory, all
        having the same ftype.

        OUTPUT: The ftype of the parent FlagAlgebra. A :class:`Flag` element

        EXAMPLES::

        The ftype of a :class:`Flag` is the same as the ftype of the 
        :class:`FlagAlgebraElement` we can construct from it ::
        
            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3)
            sage: g.ftype()
            Ftype on 0 points with edges=[]
            sage: g.ftype()==g.afae().ftype()
            True
        
        .. SEEALSO::

            :func:`ftype` in :class:`FlagAlgebra`
            :func:`ftype` in :class:`Flag`
        """
        return self._ftype
    
    def size(self):
        r"""
        Return the size of the vertex set of flags in this element
        
        OUTPUT: The size of each flag is :func:`flags`. 

        TESTS::

            sage: from sage.algebras.flag_algebras import *
            sage: FG = FlagAlgebra(QQ, GraphTheory)
            sage: FGElem = FG._an_element_()
            sage: FGElem.size() == FGElem.flags()[0].size()
            True
        """
        return self._n
    
    vertex_number = size
    
    def flags(self):
        r"""
        Returns the list of flags, corresponding to each base element
        
        The flags returned are the list of flags with size equal to 
        `self.size()` and ftype to `self.ftype()`. Their number is 
        the same as the length of `self.values()`. 

        OUTPUT: The list of flags

        EXAMPLES::

        3 vertex graphs with empty ftype ::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3)
            sage: g.afae()
            Flag Algebra Element over Rational Field
            1 - Flag on 3 points, ftype from [] with edges=[]
            0 - Flag on 3 points, ftype from [] with edges=[[0, 2]]
            0 - Flag on 3 points, ftype from [] with edges=[[0, 2], [1, 2]]
            0 - Flag on 3 points, ftype from [] with edges=[[0, 1], [0, 2], [1, 2]]
            sage: g.afae().flags()
            (Flag on 3 points, ftype from [] with edges=[],
             Flag on 3 points, ftype from [] with edges=[[0, 2]],
             Flag on 3 points, ftype from [] with edges=[[0, 2], [1, 2]],
             Flag on 3 points, ftype from [] with edges=[[0, 1], [0, 2], [1, 2]])

        .. NOTE::

            This is the same as 
            `self.theory().generate_flags(self.size(), self.ftype())`

        .. SEEALSO::

            :func:`CombinatorialTheory.generate_flags`
            :func:`size`
            :func:`ftype`
            :func:`values`
            :func:`Flag.afae`

        TESTS::

            sage: from sage.algebras.flag_algebras import *
            sage: g.afae().flags() == g.theory().generate_flags(g.size(), g.ftype())
            True
        """
        return self._flags
    
    def values(self):
        r"""
        Returns the vector of values, corresponding to each element 
        in :func:`flags`
        
        OUTPUT: A vector

        EXAMPLES::

        A flag transformed to a flag algebra element has 
        all zeroes except one entry, itself ::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3)
            sage: g.afae().values()
            (1, 0, 0, 0)

        .. SEEALSO::

            :func:`flags`
            :func:`_vector_`
            :func:`Flag.afae`
        """
        return self._values
    
    def _vector_(self, R):
        r"""
        Returns the vector of values, corresponding to each element 
        in :func:`flags` in a given base.
        
        OUTPUT: A vector

        EXAMPLES::

        A flag transformed to a flag algebra element has 
        all zeroes except one entry, itself ::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3)
            sage: g.afae()._vector_(QQ['x'])
            (1, 0, 0, 0)

        .. SEEALSO::

            :func:`values()`
            :func:`Flag.afae`
        """
        return vector(R, self._values)
    
    def __len__(self):
        r"""
        Returns the length, the number of elements
        in :func:`flags` or :func:`values` (which is the same)

        EXAMPLES::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3)
            sage: len(g.afae())
            4

        .. SEEALSO::

            :func:`flags`
            :func:`values`
            :func:`Flag.afae`
            :func:`__iter__`
        """
        return len(self._values)
    
    def __iter__(self):
        r"""
        Iterates over the elements of self
        
        It yields (number, flag) tuples, indicating the
        coefficient of the flag.

        EXAMPLES::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3)
            sage: for x in g.afae():
            ....:   print(x)
            (1, Flag on 3 points, ftype from [] with edges=[])
            (0, Flag on 3 points, ftype from [] with edges=[[0, 2]])
            (0, Flag on 3 points, ftype from [] with edges=[[0, 2], [1, 2]])
            (0, Flag on 3 points, ftype from [] with edges=[[0, 1], [0, 2], [1, 2]])


        .. SEEALSO::

            :func:`flags`
            :func:`values`
            :func:`Flag.afae`
            :func:`__len__`
        """
        for ii in range(len(self)):
            yield (self._values[ii], self._flags[ii])
    
    def _repr_(self):
        r"""
        Give a string representation
        
        Lists the flags and the corresponding coefficients,
        each on a separate line. If the list is too long 
        then only shows nonzero entries.


        EXAMPLES::

        Short list, so display all ::

            sage: from sage.algebras.flag_algebras import *
            sage: gf = GraphTheory(3).afae()
            sage: gf
            Flag Algebra Element over Rational Field
            1 - Flag on 3 points, ftype from [] with edges=[]
            0 - Flag on 3 points, ftype from [] with edges=[[0, 2]]
            0 - Flag on 3 points, ftype from [] with edges=[[0, 2], [1, 2]]
            0 - Flag on 3 points, ftype from [] with edges=[[0, 1], [0, 2], [1, 2]]
            
        Long list, only the nonzero entries are displayed 
        (for some reason this is failing the doctest, but should show 
        `Flag Algebra Element over Rational Field
        1 - Flag on 5 points, ftype from [] with edges=[]
        1 - Flag on 5 points, ftype from [] with edges=[[0, 3], [1, 4]]
        
        ` I can't figure out why)::
        
            sage: g1 = GraphTheory(5)
            sage: g2 = GraphTheory(5, edges=[[0, 1], [3, 4]])
            sage: g1+g2
            ...
            
        .. SEEALSO::

            :func:`Flag._repr_`
        """
        sttrl = ['Flag Algebra Element over {}'.format(self.parent().base())]
        strs = [str(xx) for xx in self.values()]
        maxstrlen = max([len(xx) for xx in strs])
        flgs = self.flags()
        for ii in range(len(self)):
            if len(self)<20 or abs(self.values()[ii])>=1e-8:
                sttrl.append(('{:<'+str(maxstrlen)+'} - {}').format(strs[ii], str(flgs[ii])))
        return "\n".join(sttrl)
    
    def as_flag_algebra_element(self):
        r"""
        Returns self.
        
        Only here to allow calling this function on
        both flags and flag algebra elements

        .. SEEALSO::

            :func:`Flag.afae`
        """
        return self
    
    afae = as_flag_algebra_element
    
    def _add_(self, other):
        r"""
        Adds to FlagAlgebraElements together

        OUTPUT: The sum

        EXAMPLES::

        The smaller size is shifted to match the larger ::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3).afae()
            sage: e = GraphTheory(2).afae()
            sage: e+g
            Flag Algebra Element over Rational Field
            2   - Flag on 3 points, ftype from [] with edges=[]
            2/3 - Flag on 3 points, ftype from [] with edges=[[0, 2]]
            1/3 - Flag on 3 points, ftype from [] with edges=[[0, 2], [1, 2]]
            0   - Flag on 3 points, ftype from [] with edges=[[0, 1], [0, 2], [1, 2]]

        .. NOTE::

            The result's size will match the size of the larger component

        .. SEEALSO::

            :func:`Flag._add_`
            :func:`__lshift__`
            :func:`_sub_`
        """
        nm = max(self.size(), other.size())
        vals = (self<<(nm-self.size())).values() + (other<<(nm-other.size())).values()
        return self.__class__(self.parent(), nm, vals)
    
    def _sub_(self, other):
        r"""
        Subtract a FlagAlgebraElement from this

        EXAMPLES::

        This also shifts the smaller flag to match the larger ::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3).afae()
            sage: e = GraphTheory(2).afae()
            sage: e-g
            Flag Algebra Element over Rational Field
            0   - Flag on 3 points, ftype from [] with edges=[]
            2/3 - Flag on 3 points, ftype from [] with edges=[[0, 2]]
            1/3 - Flag on 3 points, ftype from [] with edges=[[0, 2], [1, 2]]
            0   - Flag on 3 points, ftype from [] with edges=[[0, 1], [0, 2], [1, 2]]

        .. SEEALSO::

            :func:`Flag._sub_`
            :func:`__lshift__`
            :func:`_add_`
        """
        nm = max(self.size(), other.size())
        vals = (self<<(nm-self.size())).values() - (other<<(nm-other.size())).values()
        return self.__class__(self.parent(), nm, vals)
    
    def _mul_(self, other):
        r"""
        Multiplies two elements together
        
        The result will have size 
        `self.size() + other.size() - self.ftype().size()`

        EXAMPLES::

        Two empty edges multiplied together has size 4 ::

            sage: from sage.algebras.flag_algebras import *
            sage: e = GraphTheory(2).afae()
            sage: (e*e).size()
            4
            
        But if pointed (size of ftype is 1), then the size is 3 ::
            
            sage: pe = GraphTheory(2, ftype=[0])
            sage: pe*pe
            Flag Algebra Element over Rational Field
            1 - Flag on 3 points, ftype from [0] with edges=[]
            0 - Flag on 3 points, ftype from [0] with edges=[[0, 2]]
            1 - Flag on 3 points, ftype from [1] with edges=[[0, 2]]
            0 - Flag on 3 points, ftype from [0] with edges=[[0, 2], [1, 2]]
            0 - Flag on 3 points, ftype from [2] with edges=[[0, 2], [1, 2]]
            0 - Flag on 3 points, ftype from [0] with edges=[[0, 1], [0, 2], [1, 2]]
            
        Can also multiply with constants:
            
            sage: pe*3
            Flag Algebra Element over Rational Field
            3 - Flag on 2 points, ftype from [0] with edges=[]
            0 - Flag on 2 points, ftype from [0] with edges=[[0, 1]]
        
        .. NOTE::
            
            Multiplying and then projecting to a smaller ftype
            can be performed more efficiently with :func:`mul_project`
        
        .. SEEALSO::

            :func:`Flag._mul_`
            :func:`mul_project`
        """
        if self.size()==self.ftype().size():
            vals = other.values() * self.values()[0]
            return self.__class__(self.parent(), other.size(), vals)
        if other.size()==self.ftype().size():
            vals = self.values() * other.values()[0]
            return self.__class__(self.parent(), self.size(), vals)
        table = self.parent().mpt(self.size(), other.size())
        N = self.size() + other.size() - self.ftype().size()
        vals = [self.values() * mat * other.values() for mat in table]
        return self.__class__(self.parent(), N, vals)
    
    def __truediv__(self, other):
        r"""
        Divide by a scalar

        INPUT:

        - ``other`` -- number; any number such that `1` can be divided with that

        OUTPUT: The `FlagAlgebraElement` resulting from the division

        EXAMPLES::

        If 1 can be divided by that, then the division is allowed ::

            sage: from sage.algebras.flag_algebras import *
            sage: var('x')
            x
            sage: g = GraphTheory(2)
            sage: g.afae()/x
            Flag Algebra Element over Symbolic Ring
            1/x - Flag on 2 points, ftype from [] with edges=[]
            0   - Flag on 2 points, ftype from [] with edges=[[0, 1]]
        
        .. NOTE::
            
            This is the linear extension of :func:`Flag.__truediv__`
        
        .. SEEALSO::

            :func:`Flag.afae`
            :func:`Flag.__truediv__`
        """
        return self * (1/other)
    
    def __lshift__(self, amount):
        r"""
        `FlagAlgebraElement`, equal to this, with size is shifted by the amount
        
        The result will have size equal to 
        `self.size() + amount`, but the elements will be equal
        
        EXAMPLES::

        Edge shifted to size `3` ::

            sage: from sage.algebras.flag_algebras import *
            sage: edge = GraphTheory(2, edges=[[0, 1]])
            sage: (edge.afae()<<1).values()
            (0, 1/3, 2/3, 1)

        .. NOTE::
            
            This is the linear extension of :func:`Flag.__lshift__`

        .. SEEALSO::

            :func:`Flag.__lshift__`
            :func:`Flag.afae()`
        """
        if amount<0:
            raise ValueError('Can not have negative shifts')
        if amount==0:
            return self
        table = self.parent().mpt(self.size(), amount + self.ftype().size())
        vals = [sum(self.values() * mat) for mat in table]
        return self.__class__(self.parent(), self.size() + amount, vals)
    
    def project(self, ftype_inj=tuple()):
        r"""
        Project this to a smaller ftype
        

        INPUT:

        - ``ftype_inj`` -- tuple (default: (, )); the injection of the
            projected ftype inside the larger ftype

        OUTPUT: the `FlagAlgebraElement` resulting from the projection

        EXAMPLES::

        If the center of a cherry is flagged, then the projection has
        coefficient 1/3 ::

            sage: from sage.algebras.flag_algebras import *
            sage: p_cherry = GraphTheory(3, edges=[[0, 1], [0, 2]], ftype_points=[0])
            sage: p_cherry.afae().project().values()
            (0, 0, 1/3, 0)

        .. NOTE::
            
            This is the linear extension of :func:`Flag.project`

            If `ftype_inj==tuple(range(self.ftype().size()))` then this
            does nothing.

        .. SEEALSO::

            :func:`Flag.project`
        """
        return self.mul_project(1, ftype_inj)
    
    def mul_project(self, other, ftype_inj=tuple()):
        r"""
        Multiply self with other, and the project the result.

        INPUT:

        - ``ftype_inj`` -- tuple (default: (, )); the injection of the
            projected ftype inside the larger ftype

        OUTPUT: the `FlagAlgebraElement` resulting from the multiplication
            and projection

        EXAMPLES::

        Pointed edge multiplied with itself and projected ::

            sage: from sage.algebras.flag_algebras import *
            sage: felem = GraphTheory(2, edges=[[0, 1]], ftype_points=[0]).afae()
            sage: felem.mul_project(felem).values()
            (0, 0, 1/3, 1)

        .. NOTE::
            
            This is the bilinear extension of :func:`Flag.mul_project`

            If `ftype_inj==tuple(range(self.ftype().size()))` then this
            is the same as usual multiplication.

        .. SEEALSO::

            :func:`_mul_`
            :func:`project`
            :func:`Flag.mul_project`
        """
        other = self.parent(other)
        ftype_inj = tuple(ftype_inj)
        new_ftype = self.ftype().subflag([], ftype_points=ftype_inj)
        table = self.parent().mpt(self.size(), other.size(), ftype_inj)
        N = self.size() + other.size() - self.ftype().size()
        vals = [self.values() * mat * other.values() for mat in table]
        
        TargetAlgebra = FlagAlgebra(self.parent().base(), self.parent().combinatorial_theory(), new_ftype)
        return TargetAlgebra(N, vals)
    
    def density(self, other):
        r"""
        The density of self in other.
        
        Randomly choosing self.size() points in other, the
        probability of getting self.

        EXAMPLES::

        Density of an edge in the cherry graph is 2/3 ::

            sage: from sage.algebras.flag_algebras import *
            sage: cherry = GraphTheory(3, edges=[[0, 1], [0, 2]]).afae()
            sage: edge = GraphTheory(2, edges=[[0, 1]]).afae()
            sage: cherry.density(edge)
            2/3

        .. NOTE::
            
            This is the bilinear extension of :func:`Flag.mul_project`
        
        .. SEEALSO::
        
            :func:`Flag.density`
        """
        diff = self.size() - other.size()
        if diff < 0:
            raise ValueError('The target can not be larger')
        return self.values() * (other<<diff).values()
    
    def _richcmp_(self, other, op):
        r"""
        Compares the elements.
        
        Since the parent agrees, the ftype too. They are shifted to the
        same size and the values compared elementwise.

        EXAMPLES::

        Trivial example `g <= 2*g` ::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3).afae()
            sage: g <= 2*g
            True
            
        Since `g` has zero coefficients ::
        
            sage: g < 2*g
            False
            
        Shifting the size gives equal elements ::
        
            sage: g == (g<<1)
            True
            
        More sophisticated example, Mantel's theorem.
        Using the fact that for large enough structures
        `x.mul_project(x)` is always positive, we get that 
        `1/2 >= GraphTheory(2)`, so K_3 free graphs can not
        contain more than 1/2 density of K_2 ::
            
            sage: GraphTheory.exclude(GraphTheory(3))
            sage: f = GraphTheory(2, ftype=[0]) - 1/2
            sage: 1/2 >= f.mul_project(f)*2 + GraphTheory(2)
            True

        .. NOTE::
            
            When comparing Flags with FlagAlgebraElements, it
            will return False, unless the Flag is transformed
            to a FlagAlgebraElement.

        .. SEEALSO::

            :func:`mul_project`
        """
        nm = max(self.size(), other.size())
        v1 = (self<<(nm-self.size())).values()
        v2 = (other<<(nm-other.size())).values()
        return all([richcmp(v1[ii], v2[ii], op) for ii in range(len(v1))])

class FlagAlgebra(CommutativeAlgebra, UniqueRepresentation):
    def __init__(self, base, theory, ftype=None):
        r"""
        Initialize a FlagAlgebra

        INPUT:

        - ``base`` -- Ring; The base ring, this FlagAlgebra is constructed over
            This must contain the rationals `QQ`
        - ``theory`` -- CombinatorialTheory; The combinatorial theory
            this flag algebra is based on
        - ``ftype`` -- Flag (default=`None`); The ftype of the elements
            in this FlagAlgebra. The default `None` gives the empty type

        OUTPUT: The resulting FlagAlgebra

        EXAMPLES::

        Create the FlagAlgebra for GraphTheory (without any ftype) ::

            sage: from sage.algebras.flag_algebras import *
            sage: GraphFlagAlgebra = FlagAlgebra(QQ, GraphTheory)
            sage: GraphFlagAlgebra
            Flag Algebra with Ftype on 0 points with edges=[] over Rational Field
        
        Create the FlagAlgebra for TournamentTheory with point ftype ::

            sage: FlagAlgebra(QQ, TournamentTheory, TournamentTheory(1, ftype_points=[0]))
            Flag Algebra with Ftype on 1 points with edges=[] over Rational Field
        """
        if ftype==None:
            ftype = theory.empty_element()
        else:
            if not ftype.is_ftype():
                raise ValueError('{} is not an Ftype'.format(ftype))
            if ftype.theory() != theory:
                raise ValueError('{} is not a part of {}'.format(ftype, theory))
        if not base.has_coerce_map_from(QQ):
            raise ValueError('The base must contain the rationals')
        
        self._theory = theory
        self._ftype = ftype
        CommutativeAlgebra.__init__(self, base)
    
    Element = FlagAlgebraElement
    
    def _element_constructor_(self, *args, **kwds):
        r"""
        Constructs a FlagAlgebraElement with the given parameters
        
        If a single value is provided it constructs the constant in the Algebra
        If a Flag is provided it constructs the element whose only `1` value is
        that Flag.
        If a different FlagAlgebraElement is provided, then checks if the
        `theory` and the `ftype` agrees, then tries to coerce the values to the
        base of self.
        Otherwise uses the constructor of FlagAlgebraElement, which accepts a
        size value and a list of coefficients, whose length must be precisely the
        number of flags.

        EXAMPLES::

        Construct from a constant ::

            sage: from sage.algebras.flag_algebras import *
            sage: FA = FlagAlgebra(QQ, GraphTheory)
            sage: FA(3)
            Flag Algebra Element over Rational Field
            3 - Ftype on 0 points with edges=[]
        
        Construct from a flag ::
        
            sage: g = GraphTheory(2)
            sage: el = FA(g)
            sage: el
            Flag Algebra Element over Rational Field
            1 - Flag on 2 points, ftype from [] with edges=[]
            0 - Flag on 2 points, ftype from [] with edges=[[0, 1]]
            
        Construct from a FlagAlgebraElement with smaller base ::
        
            sage: FAX = FlagAlgebra(QQ['x'], GraphTheory)
            sage: FAX(el)
            Flag Algebra Element over Univariate Polynomial Ring in x over Rational Field
            1 - Flag on 2 points, ftype from [] with edges=[]
            0 - Flag on 2 points, ftype from [] with edges=[[0, 1]]
            
        Constructing the element directly from coefficients ::
            
            sage: FA(2, [3, 4])
            Flag Algebra Element over Rational Field
            3 - Flag on 2 points, ftype from [] with edges=[]
            4 - Flag on 2 points, ftype from [] with edges=[[0, 1]]

        .. SEEALSO::

            :func:`FlagAlgebraElement.__init__`
        """
        if len(args)==1:
            v = args[0]
            base = self.base()
            if isinstance(v, Flag):
                if v.ftype()==self.ftype():
                    flags = self.generate_flags(v.size())
                    vec = vector(base, [(1 if xx==v else 0) for xx in flags])
                    return self.element_class(self, v.size(), vec)
            elif isinstance(v, FlagAlgebraElement):
                if v.ftype()==self.ftype():
                    if self.base()==v.parent().base():
                        return v
                    elif self.base().has_coerce_map_from(v.parent().base()):
                        vals = vector(self.base(), v.values())
                        return self.element_class(self, v.size(), vals)
            elif v in base:
                return self.element_class(self, self.ftype().size(), vector(base, [v]))
            raise ValueError('Can\'t construct an element from {}'.format(v))
        return self.element_class(self, *args, **kwds)
    
    def _coerce_map_from_(self, S):
        r"""
        Checks if it can be coerced from S
        """
        if self.base().has_coerce_map_from(S):
            return True
        if S==self.theory():
            return True
        if isinstance(S, FlagAlgebra):
            if S.ftype()==self.ftype() and self.base().has_coerce_map_from(S.base()):
                return True
        return False
    
    def _pushout_(self, S):
        r"""
        Constructs the pushout FlagAlgebra
        """
        if S.has_coerce_map_from(self.base()):
            return FlagAlgebra(S, self.theory(), self.ftype())
        return None
    
    def _repr_(self):
        r"""
        Returns a short text representation

        EXAMPLES::

            sage: from sage.algebras.flag_algebras import *
            sage: FlagAlgebra(QQ, GraphTheory)
            Flag Algebra with Ftype on 0 points with edges=[] over Rational Field

        .. SEEALSO::

            :func:`Flag._repr_`
        """
        return 'Flag Algebra with {} over {}'.format(self.ftype(), self.base())
    
    def ftype(self):
        r"""
        Returns the ftype of this FlagAlgebra.

        EXAMPLES::

        Without specifying anything in the constructor, the ftype
        is empty ::

            sage: from sage.algebras.flag_algebras import *
            sage: FA = FlagAlgebra(QQ, GraphTheory)
            sage: FA.ftype()
            Ftype on 0 points with edges=[]

        .. NOTE::

            This is the same ftype as the ftype of the elements, 
            and the ftype of the flags in those elements. 

        .. SEEALSO::

            :func:`Flag.ftype`
            :func:`FlagAlgebraElement.ftype`
        """
        return self._ftype
    
    def combinatorial_theory(self):
        r"""
        Returns the :class:`CombinatorialTheory` object, whose
        flags form the basis of this FlagAlgebra

        EXAMPLES::

        This is the same as provided in the constructor ::

            sage: from sage.algebras.flag_algebras import *
            sage: FA = FlagAlgebra(QQ, GraphTheory)
            sage: FA.theory()
            Theory for Graph

        .. SEEALSO::

            :func:`__init__`
            :class:`CombinatorialTheory`
        """
        return self._theory
    
    theory = combinatorial_theory
    
    def base_ring(self):
        r"""
        Returns the base_ring
        
        Same as `base_ring` of the `base` provided in the constructor

        EXAMPLES::

            sage: from sage.algebras.flag_algebras import *
            sage: FA = FlagAlgebra(QQ, GraphTheory)
            sage: FA.base_ring()
            Rational Field

        .. SEEALSO::

            :func:`base`
        """
        return self.base().base_ring()
    
    def characteristic(self):
        r"""
        Returns the characteristic
        
        Same as `characteristic` of the `base` provided in the constructor

        EXAMPLES::

            sage: from sage.algebras.flag_algebras import *
            sage: FA = FlagAlgebra(QQ, GraphTheory)
            sage: FA.characteristic()
            0

        .. SEEALSO::

            :func:`base`
        """
        return self.base().characteristic()
    
    def generate_flags(self, n):
        r"""
        Generates flags of a given size, and `ftype` of `self`

        .. NOTE::

            Same as `CombinatorialTheory.generate_flags` with `ftype`
            from `self`

        .. SEEALSO::

            :func:`CombinatorialTheory.generate_flags`
        """
        return self.theory().generate_flags(n, self.ftype())
    
    def _an_element_(self):
        r"""
        Returns an element
        """
        a = self.base().an_element()
        f = self.combinatorial_theory()._an_element_(n=self.ftype().size() + 1, ftype=self.ftype())
        return self(f)*a
    
    def some_elements(self):
        r"""
        Returns a small list of elements
        """
        return [self.an_element(),self(self.base().an_element())]
    
    def mul_project_table(self, n1, n2, ftype_inj=None):
        r"""
        Returns the multiplication projection table
        
        This is the same as :func:`CombinatorialTheory.mul_project_table`
        with self.ftype() substituted in.

        .. SEEALSO::

            :func:`CombinatorialTheory.mul_project_table`
        """
        return self.theory().mul_project_table(n1, n2, self.ftype(), ftype_inj)
    
    mpt = mul_project_table


#************************************************
#
#   Implementation of a few (common) theories
#
#************************************************


def _generator_graph(n):
    r"""
    Given `n` integer, generates the graphs of size `n` using nauty
    and returns them in a dictionary required for Flag constructors
    """
    for xx in graphs.nauty_geng(str(n)):
        yield {'edges': tuple(xx.edges(labels=None))}

def _identify_graph(n, ftype_points, edges):
    r"""
    Creates a unique identifier for a graph using canonical labelings
    """
    partition = [[ii] for ii in ftype_points] + [list(set(range(n)).difference(set(ftype_points))), ]
    g = Graph([list(range(n)), edges], format='vertices_and_edges')
    blocks = tuple(g.canonical_label(partition=partition).edges(labels=None, sort=True))
    ftype_points = tuple(range(len(ftype_points)))
    return (n, ftype_points, blocks)

def _generator_threegraph(n):
    r"""
    Given `n` integer, generates the threegraphs of size `n` using nauty
    and returns them in a dictionary required for Flag constructors.
    
    Can also be modified to return hypergraphs for any size, but
    larger edge sizes result in too large theories usually.
    """
    r = 3
    for ee in range(binomial(n, r)+1):
        for xx in hypergraphs.nauty(ee, n, uniform=r):
            yield {'edges': xx}

def _identify_hypergraph(n, ftype_points, edges):
    r"""
    Identifies hypergraphs by creating a canonical label for the adjacency
    bipartite graph.
    """
    g = Graph([list(range(n+len(edges))), [(i+n,x) for i,b in enumerate(edges) for x in b]], 
              format='vertices_and_edges')
    partt = [[ii] for ii in ftype_points] + \
            [[ii for ii in range(n) if ii not in ftype_points]] + \
            [list(range(n,n+len(edges)))]
    blocks = tuple(g.canonical_label(partition=partt).edges(labels=None, sort=True))
    ftype_points = tuple(range(len(ftype_points)))
    return (n, ftype_points, blocks)

def _generator_digraph(n):
    r"""
    Given `n` integer, generates the digraphs of size `n` using nauty
    and returns them in a dictionary required for Flag constructors
    """
    gen = graphs.nauty_geng(str(n))
    for xx in digraphs.nauty_directg(gen):
        yield {'edges': tuple(xx.edges(labels=None))}

def _identify_digraph(n, ftype_points, edges):
    r"""
    Creates a unique identifier for a graph using canonical labelings
    """
    partition = [[ii] for ii in ftype_points] + [list(set(range(n)).difference(set(ftype_points))), ]
    g = DiGraph([list(range(n)), edges], format='vertices_and_edges')
    blocks = tuple(g.canonical_label(partition=partition).edges(labels=None, sort=True))
    ftype_points = tuple(range(len(ftype_points)))
    return (n, ftype_points, blocks)

def _generator_tournament(n):
    r"""
    Given `n` integer, generates the tournaments of size `n` using nauty
    and returns them in a dictionary required for Flag constructors
    """
    for xx in digraphs.tournaments_nauty(n):
        yield {'edges': tuple(xx.edges(labels=None))}

def _generator_permutation(n):
    r"""
    Given `n` integer, generates the permutations of objects `n`
    and returns the ordering as a binary relation in dictionary
    form required for Flag constructors
    """
    for perm in itertools.permutations(range(n)):
        yield {'edges': tuple(itertools.combinations(perm, r=2))}

def _identify_permutation(n, ftype_points, edges):
    r"""
    Returns a unique representation of this permutation
    """
    return (ftype_points, tuple(sorted(edges)))

def _identify_oe_graph(n, ftype_points, edges):
    r"""
    Identifies ordered edge graphs by creating a canonical label for 
    the adjacency bipartite graph.
    """
    g = Graph([list(range(n+len(edges))), [(i+n,x) for i,b in enumerate(edges) for x in b]], 
              format='vertices_and_edges')
    partt = [[ii] for ii in ftype_points] + \
            [[ii] for ii in range(n, n+len(edges))] + \
            [[ii for ii in range(n) if ii not in ftype_points]]
    blocks = tuple(g.canonical_label(partition=partt).edges(labels=None, sort=True))
    ftype_points = tuple(range(len(ftype_points)))
    return (n, ftype_points, blocks)

def _generator_oe_graph(n):
    r"""
    Given `n` integer, generates the graphs on `n` vertices
    with different (non-isomorphic) edge orderings.
    """
    for xx in graphs.nauty_geng(str(n)):
        unordered = tuple(xx.edges(labels=None))
        unique = []
        for perm in itertools.permutations(unordered):
            rel = _identify_oe_graph(n, [], perm)
            if rel not in unique:
                unique.append(rel)
                yield {'edges': perm}

def _generator_ov_graph(n):
    r"""
    Given `n` integer, generates the graphs on `n` vertices
    with different (non-isomorphic) vertex orderings.
    
    This is the same set as the boolean symmetric 
    nxn matrices with 0s on the diagonal.
    """
    full = list(itertools.combinations(range(n), int(2)))
    for ii in range(binomial(n, 2)+1):
        for xx in itertools.combinations(full, int(ii)):
            yield {'edges': xx}

def _identify_ov_graph(n, ftype_points, edges):
    r"""
    Returns a unique representation for this ordered
    vertex graph
    """
    return (n, tuple(ftype_points), tuple(sorted(list(edges))))



def _ramsey_graphs_from_graph(xx, n, excess_edges, two_edge_triples):
    r"""
    From a graph creates all the Ramsey graphs using nauty multig.
    
    A Ramsey graph here is essentially a multigraph where the edge
    colorings are represented by multi edges.
    """
    from sage.features.nauty import NautyExecutable
    import subprocess, select
    import shlex
    directg_path = NautyExecutable("multig").absolute_filename()
    edges = xx.edges(labels=None)
    le = len(edges)
    options=" -T -m2 -e{}".format(le+excess_edges)
    sub = subprocess.Popen(
        shlex.quote(directg_path) + ' {0}'.format(options),
        shell=True,
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='latin-1'
    )
    sub.stdin.write(xx.graph6_string())
    sub.stdin.close()
    unique = []
    for ll in sub.stdout:
        if ll and ll[0]==str(n):
            edges_marked = []
            seq = ll[:-1].split(" ")
            for ii in range(1, len(seq)//3):
                try:
                    ed = (int(seq[ii*3]), int(seq[ii*3 + 1]))
                except:
                    print("error on line: \n", ll, "\nhappened at graph: ", xx)
                    return
                if seq[ii*3 + 2]=="2":
                    edges_marked.append(ed)
            two_good = True
            for trp in two_edge_triples:
                count = sum([(trp[0], trp[1]) in edges_marked, 
                             (trp[0], trp[2]) in edges_marked, 
                             (trp[1], trp[2]) in edges_marked])
                if count == 1:
                    two_good = False
                    break
            if two_good:
                yield {'edges': edges, 'edges_marked': edges_marked}
    for line in sub.stderr:
        pass
    sub.wait()

def _generator_ramsey_graph(n):
    r"""
    Given `n` integer, generates the graphs on `n` vertices
    with some of the edges marked. 
    
    This is color-blind, so coloring the complement of the graph
    results in the same graph.
    """
    for xx in graphs.nauty_geng(str(n)):
        edges = xx.edges(labels=None)
        good = True
        two_edge_triples = []
        for trp in itertools.combinations(range(n), int(3)):
            count = sum([(trp[0], trp[1]) in edges, 
                         (trp[0], trp[2]) in edges, 
                         (trp[1], trp[2]) in edges])
            if count == 1:
                good = False
                break
            if count == 2:
                two_edge_triples.append(trp)
        if good: 
            le = len(edges)
            for ii in range(le//2 + 1):
                for rr in _ramsey_graphs_from_graph(xx, n, ii, two_edge_triples):
                    yield rr

def _identify_ramsey_graph(n, ftype_points, edges, edges_marked):
    r"""
    Creates a canonical label for the Ramsey graph and returns
    it as the identifier.
    """
    if len(edges_marked) == len(edges)/2:
        edges_marked_alter = [xx for xx in edges if xx not in edges_marked]
        g = Graph([list(range(n+len(edges)+len(edges_marked))), 
                   [(i+n, x) for i,b in enumerate(edges) for x in b] + 
                  [(i+n+len(edges), x) for i, b in enumerate(edges_marked) for x in b]], 
                  format='vertices_and_edges')
        partt = [[ii] for ii in ftype_points] + \
                [[ii for ii in range(n) if ii not in ftype_points]] + \
                [list(range(n,n+len(edges)))] + \
                [list(range(n+len(edges), n+len(edges)+len(edges_marked)))]
        blocks = tuple(g.canonical_label(partition=partt).edges(labels=None, sort=True))
        
        g_alter = Graph([list(range(n+len(edges)+len(edges_marked_alter))), 
                   [(i+n, x) for i,b in enumerate(edges) for x in b] + 
                  [(i+n+len(edges), x) for i, b in enumerate(edges_marked_alter) for x in b]], 
                  format='vertices_and_edges')
        partt_alter = [[ii] for ii in ftype_points] + \
                [[ii for ii in range(n) if ii not in ftype_points]] + \
                [list(range(n,n+len(edges)))] + \
                [list(range(n+len(edges), n+len(edges)+len(edges_marked_alter)))]
        blocks_alter = tuple(g_alter.canonical_label(partition=partt_alter).edges(labels=None, sort=True))
        return (n, tuple(range(len(ftype_points))), min(blocks, blocks_alter))
    else:
        if len(edges_marked) > len(edges)/2:
            edges_marked = [xx for xx in edges if xx not in edges_marked]
        g = Graph([list(range(n+len(edges)+len(edges_marked))), 
                   [(i+n, x) for i,b in enumerate(edges) for x in b] + 
                  [(i+n+len(edges), x) for i, b in enumerate(edges_marked) for x in b]], 
                  format='vertices_and_edges')
        partt = [[ii] for ii in ftype_points] + \
                [[ii for ii in range(n) if ii not in ftype_points]] + \
                [list(range(n,n+len(edges)))] + \
                [list(range(n+len(edges), n+len(edges)+len(edges_marked)))]
        blocks = tuple(g.canonical_label(partition=partt).edges(labels=None, sort=True))
        return (n, tuple(range(len(ftype_points))), blocks)

GraphTheory = CombinatorialTheory('Graph', 
                                  _generator_graph, 
                                  _identify_graph, 
                                  edges=2)

ThreeGraphTheory = CombinatorialTheory('3-Graph', 
                                       _generator_threegraph, 
                                       _identify_hypergraph, 
                                       edges=3)

DiGraphTheory = CombinatorialTheory('DiGraph', 
                                    _generator_digraph, 
                                    _identify_digraph, 
                                    edges=2)

TournamentTheory = CombinatorialTheory('Tournament', 
                                       _generator_tournament, 
                                       _identify_digraph, 
                                       edges=2)
#Note: TournamentTheory is equivalent to 
#DiGraphTheory.exclude(DiGraphTheory(2, edges=[[0, 1], [1, 0]]))

PermutationTheory = CombinatorialTheory('Permutation', 
                                        _generator_permutation, 
                                        _identify_permutation, 
                                        edges=2)

OEGraphTheory = CombinatorialTheory('OEdgeGraph', 
                                    _generator_oe_graph, 
                                    _identify_oe_graph, 
                                    edges=2)

OVGraphTheory = CombinatorialTheory('OEdgeGraph', 
                                    _generator_ov_graph, 
                                    _identify_ov_graph, 
                                    edges=2)

RamseyGraphTheory = CombinatorialTheory('RamseyGraph', 
                                        _generator_ramsey_graph, 
                                        _identify_ramsey_graph, 
                                        edges=2, 
                                        edges_marked=2)


