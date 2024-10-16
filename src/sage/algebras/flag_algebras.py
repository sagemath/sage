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
from sage.structure.element import Element
from sage.rings.rational_field import QQ
from sage.rings.semirings.non_negative_integer_semiring import NN
from sage.rings.integer import Integer
from sage.rings.real_mpfr import RR
from sage.algebras.flag import Flag

from sage.categories.sets_cat import Sets

from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix

from sage.misc.prandom import randint
from sage.arith.misc import falling_factorial, binomial
from sage.misc.functional import log, round

from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph
from sage.misc.lazy_import import lazy_import
lazy_import("sage.graphs.graph_generators", "graphs")
lazy_import("sage.graphs.digraph_generators", "digraphs")
lazy_import("sage.graphs.hypergraph_generators", "hypergraphs")

import pickle
import os
from tqdm import tqdm

class CombinatorialTheory(Parent, UniqueRepresentation):
    
    Element = Flag
    
    def __init__(self, name, generator, identifier, size_combine=None, **signature):
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
        if size_combine==None:
            self._sizes = NN
            self._size_combine = None
        else:
            self._size_combine = size_combine
            self._sizes = [ii for ii in range(100) if size_combine(0, ii, 0) == ii]
        self._excluded = []
        self._cache = {}
        self._generator = generator
        self._identifier = identifier
        self._name = name
        Parent.__init__(self, category=(Sets(), ))
        self._populate_coercion_lists_()
    
    def _compress(self, numbers):
        r"""
        Compresses the list of numbers to a short string.
        
        This is used when naming files in saving and loading, to store the parameters of this theory.

        EXAMPLES::

            sage: from sage.algebras.flag_algebras import *
            sage: GraphTheory._compress([1, 1, 10, 4, 2, 11, 15])
            '_kKgb'
        """
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
        r"""
        Clears the saved data generated by this theory, and the cached functions.
        """
        self._identify.cache_clear()
        for xx in os.listdir(self._calcs_dir()):
            if xx.startswith(self._name):
                os.remove(self._calcs_dir()+xx)
    
    def _save(self, ind, ret, is_table):
        r"""
        Saves the calculated data to persistent storage.
        """
        ns = self._calcs_dir() + self._name + "."
        
        if is_table:
            excluded, n1, n2, N, large_ftype, ftype_inj = ind
            numsind = [0]
            for xx in excluded:
                numsind += xx.raw_numbers()
            numsind += [n1, n2, N] + large_ftype.raw_numbers() + list(ftype_inj)
            if len(ret)<3000:
                self._cache[ind] = ret
        else:
            excluded, n, ftype = ind
            numsind = [1]
            for xx in excluded:
                numsind += xx.raw_numbers()
            numsind += [n] + ftype.raw_numbers()
            self._cache[ind] = ret
        save_name = ns + self._compress(numsind)
        os.makedirs(os.path.dirname(save_name), exist_ok=True)
        with open(save_name, 'wb') as file:
            pickle.dump(ret, file)
    
    def _other_save(self, name, data):
        save_name = self._calcs_dir() + name
        with open(save_name, 'wb') as file:
            pickle.dump(ret, file)
    
    def _other_load(self, name):
        save_name = self._calcs_dir() + name
        with open(save_name, 'rb') as file:
            pickle.load(ret, file)
    
    def _load(self, ind, is_table):
        r"""
        Tries to load persistent data.
        """
        ns = self._calcs_dir() + self._name + "."
        
        if ind in self._cache.keys():
            return self._cache[ind]
        
        if is_table:
            excluded, n1, n2, N, large_ftype, ftype_inj = ind
            numsind = [0]
            for xx in excluded:
                numsind += xx.raw_numbers()
            numsind += [n1, n2, N] + large_ftype.raw_numbers() + list(ftype_inj)
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
                if not is_table:
                    for xx in ret:
                        xx._set_parent(self)
                self._cache[ind] = ret
                return ret
        return None
    
    def _calcs_dir(self):
        if os.path.isdir("../calcs"):
            return "../calcs/"
        return "calcs/"
    
    def _certs_dir(self):
        if os.path.isdir("../certs"):
            return "../certs/"
        return "certs/"
    
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
        if n==0:
            #this is needed for sum to work
            alg = FlagAlgebra(QQ, self)
            return alg(0)
        if n in self.sizes():
            return self.element_class(self, n, **kwds)
        raise ValueError("For theory {}, size {} is not allowed.".format(self._name, n))
    
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
        res = [self._an_element_()]
        if 1 in self.sizes():
            res.append(self.element_class(self, 1, ftype=[0]))
        if 2 in self.sizes():
            res.append(self._an_element_(n=2))
        return res
    
    def flag_compact_repr(self, flag):
        blocks = flag.blocks()
        ret = ["n:{}".format(flag.size())]
        if len(flag.ftype_points())!=0:
            ret.append("t:"+"".join(map(str, flag.ftype_points())))
        for name in self._signature.keys():
            desc = name + ":"
            arity = self._signature[name]
            if arity==1:
                desc += "".join([str(xx[0]) for xx in blocks[name]])
            else:
                desc += ",".join(["".join(map(str, ed)) for ed in blocks[name]])
            ret.append(desc)
        return "; ".join(ret)
    
    def blowup_construction(self, target_size, pattern_size, symbolic=False, symmetric=True, unordered=False, **kwargs):
        r"""
        Returns a blowup construction, based on a given pattern

        INPUT:

        - ``target_size`` -- integer; size of the resulting FlagAlgebraElement
        - ``pattern_size`` -- integer; size of the pattern of the blowup
        - ``symbolic`` -- boolean (default: `False`); if the resulting 
            construction has part sizes symbolic variables
        - ``symmetric`` -- boolean (default: `True`); if the construction is
            symmetric. Speeds up calculation
        - ``unordered`` -- boolean (default: `False`); if the construction's parts are
            unordered. Slows down calculation
        - ``**kwargs`` -- the parameters of the pattern, one for each signature
            element.

        OUTPUT: A FlagAlgebraElement with values corresponding to the one resulting from a
            blowup construction
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(QQ, pattern_size, "X")
        gs = R.gens()
        res = 0

        if symmetric:
            terms = ((sum(gs))**target_size).dict()
            iterator = tqdm(terms)
            for exps in iterator:
                verts = []
                for ind, exp in enumerate(exps):
                    verts += [ind]*exp
                coeff = terms[exps]/(pattern_size**target_size)
                if symbolic:
                    coeff = terms[exps]
                    for ind, exp in enumerate(exps):
                        coeff *= gs[ind]**exp
                blocks = {}
                for rel in kwargs:
                    if rel not in self.signature():
                        continue
                    reledges = kwargs[rel]
                    bladd = []
                    for edge in reledges:
                        clusters = [[ii for ii in range(target_size) if verts[ii]==ee] for ee in edge]
                        bladd += list(set([tuple(sorted(xx)) for xx in itertools.product(*clusters) if len(set(xx))==len(edge)]))
                    blocks[rel] = bladd
                res += self(target_size, **blocks).afae()*coeff
        else:
            rep = int(target_size if not unordered else target_size - 1)
            iterator = tqdm(itertools.product(range(pattern_size), repeat=rep))
            for verts in iterator:
                if unordered:
                    verts = [0] + list(verts)

                coeff = 1/(pattern_size**rep)
                if symbolic:
                    coeff = 1
                    for ind in verts:
                        coeff *= gs[ind]

                blocks = {}
                for rel in kwargs:
                    if rel not in self.signature():
                        continue
                    reledges = kwargs[rel]
                    bladd = []
                    for edge in reledges:
                        clusters = [[ii for ii in range(target_size) if verts[ii]==ee] for ee in edge]
                        bladd += list(set([tuple(sorted(xx)) for xx in itertools.product(*clusters) if len(set(xx))==len(edge)]))
                    blocks[rel] = bladd
                res += self(target_size, **blocks).afae() * coeff
        return res
    
    def signature(self):
        return self._signature
    
    def sizes(self):
        return self._sizes
    
    def size_combine(self, k, n1, n2):
        if k<0 or n1<0 or n2<0:
            raise ValueError("Can't have negative size.")
        if n1<k or n2<k:
            raise ValueError("Can't have larger ftype size than flag size.")
        ret = n1+n2-k
        if self._size_combine != None:
            ret = self._size_combine(k, n1, n2)
        if ret==None:
            raise ValueError("Size combination is not allowed.")
        if ret<0:
            raise ValueError("The size combination resulted in a negative value.")
        return ret
    
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
        if n not in self.sizes():
            return None
        blocks = {key:tuple([tuple(xx) for xx in blocks[key]]) 
                  for key in blocks}
        if len(ftype_points)!=0 and not hasattr(ftype_points[0], "__getitem__"):
            ftype_points = [(ii, ) for ii in ftype_points]
        else:
            ftype_points = [tuple(xx) for xx in ftype_points]
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
            if flags.unique() != None:
                self._excluded = [flags]
            else:
                self._excluded = []
        else:
            self._excluded = [xx for xx in flags if xx.unique() != None]
    
    def _check_excluded(self, elms):
        r"""
        Helper to check the excluded structures in generation
        """
        flg = elms[0]
        for xx in elms[1]:
            if xx <= flg:
                return False
        return True
    
    def _adjust_table_phi(self, table_constructor, phi_vectors_exact, test=False):
        r"""
        Helper to modify a table constructor, incorporating extra data from
        constructions (phi_vectors_exact)
        """
        if len(phi_vectors_exact)==0:
            return table_constructor

        for param in table_constructor.keys():
            ns, ftype, target_size = param
            table = self.mul_project_table(ns, ns, ftype, ftype_inj=[], target_size=target_size)
            Zs = [[None for _ in range(len(phi_vectors_exact))] for _ in range(len(table_constructor[param]))]
            for gg, morig in enumerate(table):
                for ii, base in enumerate(table_constructor[param]):
                    mat = base * morig * base.T
                    for phind, phi_vector_exact in enumerate(phi_vectors_exact):
                        if Zs[ii][phind]==None:
                            Zs[ii][phind] = mat*phi_vector_exact[gg]
                        else:
                            Zs[ii][phind] += mat*phi_vector_exact[gg]

            new_bases = []
            for ii, Zgroup in enumerate(Zs):
                Z = None
                for Zjj in Zgroup:
                    if test and (not Zjj.is_positive_semidefinite()):
                        print("Construction based Z matrix for {} is not semidef: {}".format(ftype, min(Zjj.eigenvalues())))
                    if Z==None:
                        Z = Zjj
                    else:
                        Z.augment(Zjj)
                Zk = Z.kernel()
                Zkern = Zk.basis_matrix()
                if Zkern.nrows()>0:
                    new_bases.append(matrix(QQ, Zkern * table_constructor[param][ii], sparse=True))
            table_constructor[param] = new_bases

        return table_constructor

    def _print_eigenvalues(self, table_constructor, sdp_result):
        r"""
        Helper to quickly print the eigenvalues of each X matrix from the result
        of an SDP.
        """
        block_index = 0
        for params in table_constructor.keys():
            ns, ftype, target_size = params
            table = self.mul_project_table(ns, ns, ftype, ftype_inj=[], target_size=target_size)

            for plus_index, base in enumerate(table_constructor[params]):
                X_approx = matrix(sdp_result['X'][block_index + plus_index])
                X_eigenvalues = X_approx.eigenvalues()
                print("{} index {} has eigenvalues {}\n\n".format(ftype, plus_index, X_eigenvalues))
            block_index += len(table_constructor[params])
    
    def _tables_to_sdp_data(self, table_constructor, prev_data=None):
        r"""
        Helper to transform the data from the multiplication tables to an SDP input
        """
        if prev_data==None:
            block_sizes = []
            target = []
            mat_inds = []
            mat_vals = []
        else:
            block_sizes, target, mat_inds, mat_vals = prev_data
        block_index = len(block_sizes) + 1
        for params in table_constructor.keys():
            ns, ftype, target_size = params
            table = self.mul_project_table(ns, ns, ftype, ftype_inj=[], target_size=target_size)
            block_sizes += [base.nrows() for base in table_constructor[params]]

            #only loop through the table once
            for gg, morig in enumerate(table):
                #for each base change create the entries
                for plus_index, base in enumerate(table_constructor[params]):
                    mm = base * morig * base.T
                    dd = mm._dict()
                    if len(dd)>0:
                        inds, values = zip(*mm._dict().items())
                        iinds, jinds = zip(*inds)
                        for cc in range(len(iinds)):
                            if iinds[cc]>=jinds[cc]:
                                mat_inds.extend([gg+1, block_index + plus_index, iinds[cc]+1, jinds[cc]+1])
                                mat_vals.append(values[cc])
            block_index += len(table_constructor[params])
        return block_sizes, target, mat_inds, mat_vals
    
    def _constraints_to_sdp_data(self, constraints_data, prev_data=None):
        r"""
        Helper to transform the data from the constraints to an SDP input data
        """
        if prev_data==None:
            block_sizes = []
            target = []
            mat_inds = []
            mat_vals = []
        else:
            block_sizes, target, mat_inds, mat_vals = prev_data
        flag_num, constraints_vals, constraints_flags_vec, one_vector = constraints_data
        block_index = len(block_sizes) + 1
        constr_num = len(constraints_vals)
        for ii in range(constr_num):
            mat_inds.extend([0, block_index+1, 1+ii, 1+ii])
            mat_vals.append(constraints_vals[ii])

        for gg in range(flag_num):
            mat_inds.extend([gg+1, block_index, gg+1, gg+1])
            mat_vals.append(1)
            for ii in range(constr_num):
                mat_inds.extend([gg+1, block_index+1, ii+1, ii+1])
                mat_vals.append(constraints_flags_vec[ii][gg])
        block_sizes += [-flag_num, -constr_num]
        return block_sizes, target, mat_inds, mat_vals
    
    def _target_to_sdp_data(self, target, prev_data=None):
        r"""
        Helper to transform the target to an SDP input data
        """
        if prev_data==None:
            return [], list(target), [], []
        prev_data[1] = list(target)
        return prev_data
    
    def _get_relevant_ftypes(self, target_size):
        r"""
        Returns the ftypes useful for optimizing up to `target_size`
        """
        plausible_sizes = []
        for fs in self.sizes():
            if fs>=target_size:
                break
            if fs==0:
                continue
            plausible_sizes.append(fs)
        ftype_pairs = []
        for fs, ns in itertools.combinations(plausible_sizes, r=int(2)):
            if self.size_combine(fs, ns, ns) <= target_size:
                kk = ns-fs
                found = False
                for ii, (bfs, bns) in enumerate(ftype_pairs):
                    if bns-bfs==kk:
                        found = True
                        if ns>bns:
                            ftype_pairs[ii]=(fs, ns)
                        break
                if not found:
                    ftype_pairs.append((fs, ns))

        ftype_data = []
        for fs, ns in ftype_pairs:
            ftype_flags = self.generate_flags(fs)
            ftypes = [flag.subflag([], ftype_points=list(range(fs))) for flag in ftype_flags]
            for xx in ftypes:
                ftype_data.append((ns, xx, target_size))
        ftype_data.sort()
        return ftype_data
    
    def _create_table_constructor(self, ftype_data, target_size):
        r"""
        Table constructor is a dictionary that holds the data to construct
        all the multiplication tables. 
        
        For each ftype and base change it provides the data to create the multiplication table.
        Also pre-computes the multiplication tables if they are not calculated yet.
        """
        
        sym_asym_mats = [self.sym_asym_bases(dat[0], dat[1]) for dat in ftype_data]

        table_constructor = {}
        for ii, dat in (pbar := tqdm(enumerate(ftype_data))):
            ns, ftype, target_size = dat
            #pre-calculate the table here
            table = self.mul_project_table(ns, ns, ftype, ftype_inj=[], target_size=target_size)
            if table==None:
                pbar.set_description("Structures with size {} and {} had singular multiplication table!".format(ns, ftype))
                continue
            sym_base, asym_base = sym_asym_mats[ii]
            bases = []
            if sym_base.nrows()!=0:
                bases.append(sym_base)
            if asym_base.nrows()!=0:
                bases.append(asym_base)
            table_constructor[dat] = bases
            pbar.set_description("Done with mult table for {}".format(ftype))
        return table_constructor
    
    def _create_constraints_data(self, positives, target_element, target_size):
        r"""
        Creates the data that holds the linear constraints
        """
        
        
        base_flags = self.generate_flags(target_size)
        
        if positives == None:
            positives_list_exact = []
            constraints_vals = []
        else:
            positives_list_exact = []
            for ii in (pbar:= tqdm(range(len(positives)))):
                fv = positives[ii]
                if isinstance(fv, Flag):
                    continue
                kf = fv.ftype().size()
                nf = fv.size()
                if self._size_combine == None:
                    df = target_size - nf + kf
                else:
                    df = -1
                    for xx in self.sizes():
                        if self._size_combine(kf, nf, xx)==target_size:
                            df = xx
                            break
                mult_table = self.mul_project_table(nf, df, fv.ftype(), ftype_inj=[], target_size=target_size)
                fvvals = fv.values()
                m = matrix(QQ, [vector(fvvals*mat) for mat in mult_table])
                positives_list_exact += list(m.T)
                pbar.set_description("Done with positivity constraint {}".format(ii))
            constraints_vals = [0]*len(positives_list_exact)
        
        # The one vector is also calculated here and is a linear constraint
        if target_element.ftype().size()==0:
            one_vector = vector([1]*len(base_flags))
        else:
            one_vector = (target_element.ftype().project()<<(target_size - target_element.ftype().size())).values()
        positives_list_exact.extend([one_vector, one_vector*(-1)])
        constraints_vals.extend([1, -1])
        
        return len(base_flags), constraints_vals, positives_list_exact, one_vector
    
    
    def _round_sdp_solution_no_phi(self, sdp_result, sdp_data, table_constructor, \
                                   constraints_data, denom=1024):
        import numpy as np
        from numpy import linalg as LA
        from sage.functions.other import ceil

        #unpack variables

        block_sizes, target_list_exact, mat_inds, mat_vals = sdp_data
        target_vector_exact = vector(target_list_exact)
        flags_num, constraints_vals, positives_list_exact, one_vector = constraints_data
        positives_matrix_exact = matrix(QQ, len(positives_list_exact), flags_num, positives_list_exact)
        
        one_vector_exact = positives_matrix_exact.rows()[-2] # find the one_vector from the equality constraint
        positives_matrix_exact = positives_matrix_exact[:-2, :] # remove the equality constraints

        flags_num = -block_sizes[-2] # same as |F_n|

        X_matrices_approx = sdp_result['X'][:-2]
        X_matrices_rounded = []
        print("Rounding X matrices")
        for X in tqdm(X_matrices_approx):

            Xr = _round_matrix(X, method=0, denom=denom)
            Xnp = np.array(Xr)
            eigenvalues, eigenvectors = LA.eig(Xnp)
            emin = min(eigenvalues)
            if emin<0:
                eminr = ceil(-emin*denom)/denom
                Xr = matrix(QQ, Xr) + diagonal_matrix(QQ, [eminr]*len(X), sparse=True)
            X_matrices_rounded.append(Xr)
        X_matrices_flat = [vector(_flatten_matrix(X.rows(), doubled=False)) for X in (X_matrices_rounded)]

        e_vector_approx = sdp_result['X'][-1][:-2]
        e_vector_rounded = vector(_round_list(e_vector_approx, force_pos=True, method=0, denom=denom))
        
        phi_vector_approx = sdp_result['y']
        phi_vector_rounded = vector(_round_list(phi_vector_approx, force_pos=True, method=0, denom=denom))

        slacks = target_vector_exact - positives_matrix_exact.T*e_vector_rounded
        block_index = 0
        print("Calculating resulting bound")
        for params in tqdm(table_constructor.keys()):
            ns, ftype, target_size = params
            table = self.mul_project_table(ns, ns, ftype, ftype_inj=[], target_size=target_size)
            for gg, morig in enumerate(table):
                for plus_index, base in enumerate(table_constructor[params]):
                    block_dim = block_sizes[block_index + plus_index]
                    X_flat = X_matrices_flat[block_index + plus_index]
                    M = base * morig * base.T
                    M_flat_vector_exact = vector(_flatten_matrix(M.rows(), doubled=True))
                    slacks[gg] -= M_flat_vector_exact*X_flat
            block_index += len(table_constructor[params])
        # scale back slacks with the one vector, the minimum is the final result
        result = min([slacks[ii]/oveii for ii, oveii in enumerate(one_vector_exact) if oveii!=0])
        # pad the slacks, so it is all positive where it counts
        slacks -= result*one_vector_exact
        
        print("The rounded result is {}".format(result))
        
        return result, X_matrices_rounded, e_vector_rounded, slacks, [phi_vector_rounded]
    
    def _round_sdp_solution_phi(self, sdp_result, sdp_data, table_constructor, \
                            constraints_data, phi_vectors_exact, denom=1024):
        r"""
        Round the SDP results output to get something exact.
        """
        
        #unpack variables
        block_sizes, target_list_exact, mat_inds, mat_vals = sdp_data
        target_vector_exact = vector(target_list_exact)
        flags_num, constraints_vals, positives_list_exact, one_vector = constraints_data
        positives_matrix_exact = matrix(QQ, len(positives_list_exact), flags_num, positives_list_exact)
        
        no_constr = len(phi_vectors_exact)==0
        phi_vector_exact = vector([0]*positives_matrix_exact.ncols()) if no_constr else phi_vectors_exact[0]
        
        one_vector_exact = positives_matrix_exact.rows()[-2] # find the one_vector from the equality constraint
        positives_matrix_exact = positives_matrix_exact[:-2, :] # remove the equality constraints

        flags_num = -block_sizes[-2] # same as |F_n|

        c_vector_approx = vector(sdp_result['X'][-2]) # dim: |F_n|, c vector, primal slack for flags
        c_vector_rounded = vector(_round_list(c_vector_approx, method=0, denom=denom)) # as above but rounded

        # The F (FF) flag indecies where the c vector is zero/nonzero
        c_zero_inds = [FF for FF, xx in enumerate(c_vector_approx) if (abs(xx)<1e-6 or phi_vector_exact[FF]!=0)]
        c_nonzero_inds = [FF for FF in range(flags_num) if FF not in c_zero_inds]



        positives_num = -block_sizes[-1] - 2 # same as m, number of positive constraints (-2 for the equality)

        phi_pos_vector_exact = positives_matrix_exact*phi_vector_exact # dim: m, witness that phi is positive

        e_vector_approx = vector(sdp_result['X'][-1][:-2]) # dim: m, the e vector, primal slack for positivitives
        e_vector_rounded = vector(_round_list(e_vector_approx, method=0, denom=denom)) # as above but rounded

        # The f (ff) positivity constraints where the e vector is zero/nonzero
        e_zero_inds = [ff for ff, xx in enumerate(e_vector_approx) if (abs(xx)<1e-6 or phi_pos_vector_exact[ff]!=0)]
        e_nonzero_inds = [ff for ff in range(positives_num) if ff not in e_zero_inds]



        bound_exact = target_vector_exact*phi_vector_exact 
        # the constraints for the flags that are exact
        corrected_target_relevant_exact = vector([target_vector_exact[FF] - bound_exact for FF in c_zero_inds])
        # the d^f_F matrix, but only the relevant parts for the rounding
        # so F where c_F = 0 and f where e_f != 0
        positives_matrix_relevant_exact = matrix(QQ, len(e_nonzero_inds), len(c_zero_inds), \
                                                 [[positives_matrix_exact[ff][FF] for FF in c_zero_inds] for ff in e_nonzero_inds])
        # the e vector, but only the nonzero entries
        e_nonzero_list_rounded = [e_vector_rounded[ff] for ff in e_nonzero_inds]
        
        
        
        # 
        # Flatten the matrices relevant for the rounding
        # 
        # M table transforms to a matrix, (with nondiagonal entries doubled)
        # only the FF index matrices corresponding with tight constraints are used
        # 
        # X transforms to a vector
        # only the semidefinite blocks are used
        # 

        # The relevant entries of M flattened to a matrix this will be indexed by 
        # c_zero_inds and the triples from the types
        
        print("Rounding X matrices")
        
        M_flat_relevant_matrix_exact = matrix(QQ, len(c_zero_inds), 0, 0, sparse=True)
        X_flat_vector_rounded = [] # The rounded X values flattened to a list
        block_index = 0
        block_info = []
        for params in table_constructor.keys():
            ns, ftype, target_size = params
            table = self.mul_project_table(ns, ns, ftype, ftype_inj=[], target_size=target_size)

            for plus_index, base in enumerate(table_constructor[params]):
                block_info.append([ftype, base])
                X_approx = sdp_result['X'][block_index + plus_index]
                X_flat_vector_rounded += _round_list(_flatten_matrix(X_approx), method=0, denom=denom)

                M_extra = []

                for FF in c_zero_inds:
                    M_FF = table[FF]
                    M_extra.append(_flatten_matrix((base * M_FF * base.T).rows(), doubled=True))

                M_flat_relevant_matrix_exact = M_flat_relevant_matrix_exact.augment(matrix(M_extra))
            block_index += len(table_constructor[params])


        # 
        # Append the relevant M matrix and the X with the additional values from
        # the positivity constraints. 
        #
        # Then correct the x vector values
        # 

        M_matrix_final = M_flat_relevant_matrix_exact.augment(positives_matrix_relevant_exact.T)
        x_vector_final = vector(X_flat_vector_rounded+e_nonzero_list_rounded)


        # Correct the values of the x vector, based on the minimal L_2 norm
        x_vector_corr = x_vector_final - M_matrix_final.T * \
        (M_matrix_final * M_matrix_final.T).pseudoinverse() * \
        (M_matrix_final*x_vector_final - corrected_target_relevant_exact) 
        
        #
        # Recover the X matrices and e vector from the corrected x
        #

        e_nonzero_vector_corr = x_vector_corr[-len(e_nonzero_inds):]
        e_vector_corr = vector(QQ, positives_num, dict(zip(e_nonzero_inds, e_nonzero_vector_corr)))
        if len(e_vector_corr)>0 and min(e_vector_corr)<0:
            print("Linear coefficient is negative: {}".format(min(e_vector_corr)))
            return None
        
        X_final = []
        slacks = target_vector_exact - positives_matrix_exact.T*e_vector_corr
        block_index = 0
        print("Calculating resulting bound")
        for params in tqdm(table_constructor.keys()):
            ns, ftype, target_size = params
            table = self.mul_project_table(ns, ns, ftype, ftype_inj=[], target_size=target_size)
            for plus_index, base in enumerate(table_constructor[params]):
                block_dim = block_sizes[block_index + plus_index]
                X_ii_small, x_vector_corr = _unflatten_matrix(x_vector_corr, block_dim)
                X_ii_small = matrix(X_ii_small)
                
                # verify semidefiniteness
                if not X_ii_small.is_positive_semidefinite():
                    print("Rounded X matrix {} is not semidefinite: {}".format(block_index+plus_index, min(X_ii_small.eigenvalues())))
                    return None
                
                # update slacks
                for gg, morig in enumerate(table):
                    M = base * morig * base.T
                    M_flat_vector_exact = vector(_flatten_matrix(M.rows(), doubled=True))
                    slacks[gg] -= M_flat_vector_exact*vector(_flatten_matrix(X_ii_small.rows(), doubled=False))
                
                X_final.append(X_ii_small)
            block_index += len(table_constructor[params])
        
        # scale back slacks with the one vector, the minimum is the final result
        result = min([slacks[ii]/oveii for ii, oveii in enumerate(one_vector_exact) if oveii!=0])
        # pad the slacks, so it is all positive where it counts
        slacks -= result*one_vector_exact
        
        return result, X_final, e_vector_corr, slacks, phi_vectors_exact
    
    
    def _fix_X_bases(self, table_constructor, X_original):
        r"""
        Transforms the X matrices to a base that agrees with the original list of flags
        
        Basically undoes the sym/asym changes and the reductions by the constructions' kernels.
        """
        if X_original==None:
            return None
        X_flats = []
        block_index = 0
        for params in table_constructor.keys():
            ns, ftype, target_size = params
            X_ii = None
            for plus_index, base in enumerate(table_constructor[params]):
                if X_ii == None:
                    X_ii = base.T * X_original[block_index + plus_index] * base
                else:
                    X_ii += base.T * X_original[block_index + plus_index] * base
            block_index += len(table_constructor[params])
            X_flats.append(vector(QQ, _flatten_matrix(X_ii.rows())))
        return X_flats
    
    def _fix_X_bases_pld(self, table_constructor, X_original):
        r"""
        Transforms the X matrices to a base that agrees with the original list of flags
        
        Basically undoes the sym/asym changes and the reductions by the constructions' kernels.
        Also changes the X matrices to the P L D L.T P.T form
        """
        if X_original==None:
            return None
        P_dicts = []
        L_mats = []
        D_vecs = []
        block_index = 0
        for params in table_constructor.keys():
            ns, ftype, target_size = params
            X_ii = None
            for plus_index, base in enumerate(table_constructor[params]):
                if X_ii == None:
                    X_ii = base.T * X_original[block_index + plus_index] * base
                else:
                    X_ii += base.T * X_original[block_index + plus_index] * base
            block_index += len(table_constructor[params])
            P, L, D = X_ii.block_ldlt()
            # last check that D has only diagonals
            if not all([xx[0]==xx[1] for xx in D._dict().keys()]):
                print("LDLT factoring failed")
                return None
            P_dicts.append(P._dict())
            L_mats.append(vector(_flatten_matrix(L.T.rows())))
            D_vecs.append(vector(D.diagonal()))
        return P_dicts, L_mats, D_vecs
    
    def _format_optimizer_output(self, table_constructor, mult=1, sdp_output=None, rounding_output=None, file="default"):
        r"""
        Formats the outputs to a nice certificate
        
        The result contains: the final bound, the X matrices, the linear coefficients, the slacks,
                            a guess or exact construction, the list of base flags, the list of used (typed) flags
        """
        
        is_json = file.endswith(".json")
        
        target_size = 0
        typed_flags = {}
        for params in table_constructor.keys():
            ns, ftype, target_size = params
            if is_json:
                safekey = "target:" + str(ns) + "; " + self.flag_compact_repr(ftype)
                typed_flags[safekey] = [self.flag_compact_repr(xx) for xx in self.generate_flags(ns, ftype)]
            else:
                typed_flags[(ns, ftype)] = self.generate_flags(ns, ftype)
        
        if is_json:
            base_flags = [self.flag_compact_repr(xx) for xx in self.generate_flags(target_size)]
        else:
            base_flags = self.generate_flags(target_size)
        
        result = None
        X_original = None
        e_vector = None
        slacks = None
        phi_vecs = None
        
        if sdp_output!=None:
            result = sdp_output['primal']
            X_original = [matrix(dat) for dat in sdp_output['X'][:-2]]
            e_vector = vector(sdp_output['X'][-1])
            slacks = vector(sdp_output['X'][-2])
            phi_vecs = [vector(sdp_output['y'])]
        elif rounding_output!=None:
            result, X_original, e_vector, slacks, phi_vecs = rounding_output
        
        result *= mult
        #Ps, Ls, Ds = self._fix_X_bases_pld(table_constructor, X_original)
        Xs = self._fix_X_bases(table_constructor, X_original)
        
        cert_dict = {"result": result, 
                     "X matrices": Xs,
                     "e vector": e_vector,
                     "slack vector": slacks,
                     "phi vectors": phi_vecs,
                     "base flags": base_flags,
                     "typed flags": typed_flags
                    }
        
        if file!="default":
            if is_json:
                import json
                file = self._certs_dir() + file
                with open(file, "w") as f:
                    json.dump(cert_dict, f, indent=2, default=str)
            else:
                import pickle
                if not file.endswith(".pickle"):
                    file += ".pickle"
                file = self._certs_dir() + file
                with open(file, "wb") as f:
                    pickle.dump(cert_dict, f)
        
        return cert_dict
    
    def optimize_problem(self, target_element, target_size, maximize=True, positives=None, \
                         construction=None, file=None, exact=False, denom=1024):
        r"""
        Try to maximize or minimize the value of `target_element`
        
        The algorithm calculates the multiplication tables and 
        sends the SDP problem to CSDPY.
        
        INPUT:
 
        - ``target_element`` -- Flag or FlagAlgebraElement; 
            the target whose density this function tries to
            maximize or minimize in large structures.
        - ``target_size`` -- integer; The program evaluates
            flags and the relations between them up to this
            size.
        - ``maximize`` -- boolean (default: `True`); 
        - ``file`` -- file to save the certificate 
            (default: `None`); Use None to not create 
            certificate
        - ``positives`` -- list of flag algebra elements, 
            optimizer will assume those are positive, can
            have different types
        - ``construction`` -- a list or a single element of 
            `FlagAlgebraElement`s; to consider in the kernel
            `None` by default, in this case the construction
            is guessed from a preliminary run.
        - ``exact`` -- boolean; to round the result or not

        OUTPUT: A bound for the optimization problem. If 
            certificate is requested then returns the entire
            output of the solver as the second argument.
        
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
        from csdpy import solve_sdp
        import sys
        import io
        import time

        #
        # Initial setup
        #

        if target_size not in self.sizes():
            raise ValueError("For theory {}, size {} is not allowed.".format(self._name, target_size))

        base_flags = self.generate_flags(target_size)
        print("Base flags generated, their number is {}".format(len(base_flags)))
        mult = -1 if maximize else 1
        target_vector_exact = (target_element.project()*(mult)<<(target_size - target_element.size())).values()
        sdp_data = self._target_to_sdp_data(target_vector_exact)
        
        #
        # Create the relevant ftypes
        #
        
        ftype_data = self._get_relevant_ftypes(target_size)
        print("The relevant ftypes are constructed, their number is {}".format(len(ftype_data)))
        flags = [self.generate_flags(dat[0], dat[1]) for dat in ftype_data]
        flag_sizes = [len(xx) for xx in flags]
        print("Block sizes before symmetric/asymmetric change is applied: {}".format(flag_sizes))
        
        #
        # Create the table constructor and add it to sdp_data
        #
        
        table_constructor = self._create_table_constructor(ftype_data, target_size)
        sdp_data = self._tables_to_sdp_data(table_constructor, prev_data=sdp_data)
        print("Tables finished")

        #
        # Add constraints data and add it to sdp_data
        #
        
        constraints_data = self._create_constraints_data(positives, target_element, target_size)
        sdp_data = self._constraints_to_sdp_data(constraints_data, prev_data=sdp_data)
        print("Constraints finished")
        
        #
        # If construction is None or [] then run the optimizer without any construction
        #
        
        if construction==None or construction==[]:
            print("Running sdp without construction. Used block sizes are {}".format(sdp_data[0]))
            
            time.sleep(float(0.1))
            initial_sol = solve_sdp(*sdp_data)
            time.sleep(float(0.1))
            
            # Format the result and return it if floating point values are fine
            if (not exact):
                if file==None:
                    return initial_sol['primal'] * mult
                return self._format_optimizer_output(table_constructor, mult=mult, sdp_output=initial_sol, file=file)
            
            # Guess the construction in this case
            if construction==None:
                one_vector = constraints_data[-1]
                phi_vector_original = initial_sol['y']
                phi_vector_rounded, error_coeff = _round_adaptive(initial_sol['y'], one_vector)
                if error_coeff<1e-6:
                    alg = FlagAlgebra(QQ, self)
                    construction = alg(target_size, phi_vector_rounded)
                    phipr = str(construction)
                    print("The initial run gave an accurate looking construction")
                    if len(phipr)<1000:
                        print("Rounded construction vector is: \n{}".format(phipr))
                else:
                    print("The initial run didn't provide an accurate construction")
                    construction = []
            
            # If nothing was provided or the guess failed, then round the current solution
            if construction==[]:
                rounding_output = self._round_sdp_solution_no_phi(initial_sol, sdp_data, table_constructor, \
                                                                 constraints_data, denom=denom)
                if file==None:
                    return rounding_output[0] * mult
                return self._format_optimizer_output(table_constructor, mult=mult, rounding_output=rounding_output, file=file)
        
        
        #
        # Run the optimizer (again if a construction was guessed) with the construction
        #
        
        if isinstance(construction, FlagAlgebraElement):
            phi_vectors_exact = [construction.values()]
        else:
            phi_vectors_exact = [xx.values() for xx in construction]

        #
        # Adjust the table to consider the kernel from y_rounded
        #

        print("Adjusting table with kernels from construction")
        table_constructor = self._adjust_table_phi(table_constructor, phi_vectors_exact)
        sdp_data = self._target_to_sdp_data(target_vector_exact)
        sdp_data = self._tables_to_sdp_data(table_constructor, prev_data=sdp_data)
        sdp_data = self._constraints_to_sdp_data(constraints_data, prev_data=sdp_data)
        
        #
        # Then run the optimizer
        #
        
        print("Running SDP after kernel correction. Used block sizes are {}".format(sdp_data[0]))
        time.sleep(float(0.1))
        final_sol = solve_sdp(*sdp_data)
        time.sleep(float(0.1))

        # Quickly deal with the case when no rounding is needed
        if (not exact):
            if file==None:
                return final_sol['primal'] * mult
            return self._format_optimizer_output(table_constructor, mult=mult, sdp_output=final_sol, file=file)
        
        
        print("Starting the rounding of the result")
        rounding_output = self._round_sdp_solution_phi(final_sol, sdp_data, table_constructor, \
                                                                         constraints_data, phi_vectors_exact, denom=denom)
        if rounding_output==None:
            print("Rounding based on construction was unsuccessful")
            rounding_output = self._round_sdp_solution_no_phi(final_sol, sdp_data, table_constructor, \
                                                            constraints_data, denom=denom)
        
        print("Final rounded bound is {}".format(rounding_output[0]*mult))
        
        if file==None:
            return rounding_output[0] * mult
        return self._format_optimizer_output(table_constructor, mult=mult, rounding_output=rounding_output, file=file)
    
    
    optimize = optimize_problem
    
    
    def verify_certificate(self, file_or_cert, target_element, target_size, maximize=True, positives=None, construction=None):
        r"""
        Verifies the certificate provided by the optimizer written to `file`
        """
        import pickle
        
        #
        # Load the certificate file (only pickle is supported)
        #
        
        if isinstance(file_or_cert, str):
            file = file_or_cert
            if not file.endswith(".pickle"):
                file += ".pickle"
            file = self._certs_dir() + file
            with open(file, 'rb') as file:
                certificate = pickle.load(file)
        else:
            certificate = file_or_cert
        
        
        #
        # Checking eigenvalues and positivity constraints
        #
        
        res = certificate["result"]
        e_values = vector(certificate["e vector"])

        if len(e_values)>0 and min(e_values)<0:
            print("Solution is not valid!")
            print("Linear constraint's coefficient is negative {}".format(min(e_values)))
            return -1
        
        X_flats = certificate["X matrices"]
        #Ps = certificate["P dictionaries"]
        #Ds = certificate["D vectors"]
        #Ls = certificate["L matrices"]
        print("Checking X matrices")
        for ii, Xf in tqdm(enumerate(X_flats)):
            X = matrix(QQ, _unflatten_matrix(Xf)[0])
            if not (X.is_positive_semidefinite()):
                print("Solution is not valid!")
                print("Matrix {} is not semidefinite".format(ii))
                return -1
            #P = matrix(QQ, len(Ddiag), len(Ddiag), Ps[ii], sparse=True)
            #D = diagonal_matrix(QQ, Ddiag, sparse=True)
            #Larr, _ = _unflatten_matrix(Ls[ii], dim=len(Ddiag), doubled=False, upper=True)
            #L = matrix(QQ, Larr).T
            #PL = P*L
            #X = PL * D * PL.T
            #X_flats.append(vector(QQ, _flatten_matrix(X.rows())))

        print("Solution matrices are all positive semidefinite, linear coefficients are all non-negative")

        #
        # Initial setup
        #

        mult = -1 if maximize else 1
        base_flags = certificate["base flags"]
        target_vector_exact = (target_element.project()*(mult)<<(target_size - target_element.size())).values()
        if target_element.ftype().size()==0:
            one_vector = vector([1]*len(base_flags))
        else:
            one_vector = (target_element.ftype().project()<<(target_size - target_element.ftype().size())).values()
        
        ftype_data = list(certificate["typed flags"].keys())

        #
        # Create the semidefinite matrix data
        #

        table_list = []
        print("Calculating multiplication tables")
        for ii, dat in tqdm(enumerate(ftype_data)):
            ns, ftype = dat
            #calculate the table
            table = self.mul_project_table(ns, ns, ftype, ftype_inj=[], target_size=target_size)
            if table!=None:
                table_list.append(table)

        #
        # Create the data from linear constraints
        #

        positives_list_exact = []
        if positives != None:
            for ii, fv in enumerate(positives):
                if isinstance(fv, Flag):
                    continue
                nf = fv.size()
                df = target_size + fv.ftype().size() - nf
                mult_table = self.mul_project_table(nf, df, fv.ftype(), ftype_inj=[], target_size=target_size)
                fvvals = fv.values()
                m = matrix(QQ, [vector(fvvals*mat) for mat in mult_table])
                positives_list_exact += list(m.T)
                print("Done with positivity constraint {}".format(ii))
        positives_matrix_exact = matrix(QQ, len(positives_list_exact), len(base_flags), positives_list_exact)

        print("Done calculating linear constraints")

        #
        # Calculate the bound the solution provides
        #
        
        print("Calculating the bound provided by the certificate")
        
        slacks = target_vector_exact - positives_matrix_exact.T*e_values
        for ii, table in tqdm(enumerate(table_list)):
            for gg, mat_gg in enumerate(table):
                mat_flat = vector(_flatten_matrix(mat_gg.rows(), doubled=True))
                slacks[gg] -= mat_flat * X_flats[ii]
        result = min([slacks[ii]/oveii for ii, oveii in enumerate(one_vector) if oveii!=0])
        result *= mult
        print("The solution is valid, it proves the bound {}".format(result))

        return result
    
    verify = verify_certificate
    
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
        if n not in self.sizes():
            raise ValueError("For theory {}, size {} is not allowed.".format(self._name, n))
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
    
    generate = generate_flags
    
    def mul_project_table(self, n1, n2, large_ftype, ftype_inj=None, target_size=None):
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
        if target_size==None:
            target_size = self.size_combine(large_size, n1, n2)
        elif target_size not in self.sizes():
            raise ValueError("For theory {}, size {} is not allowed.".format(self._name, target_size))
        return self._mpte(tuple(self._excluded), target_size, n1, n2, large_ftype, ftype_inj)
    
    mpt = mul_project_table
    table = mul_project_table
    
    def _mpte(self, excluded, N, n1, n2, large_ftype, ftype_inj):
        r"""
        The (hidden) cached version of :func:`mul_project_table`
        """
        ind = (excluded, N, n1, n2, large_ftype, ftype_inj)
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
        
        Nflgs = self._gfe(excluded, N, small_ftype)
        n1flgs = self._gfe(excluded, n1, large_ftype)
        n2flgs = self._gfe(excluded, n2, large_ftype)
        
        slist = tuple((flg, n1, n1flgs, n2, n2flgs, ftype_remap, large_ftype, small_ftype) for flg in Nflgs)
        
        pool = mp.Pool(mp.cpu_count()-1)
        mats = pool.map(self._density_wrapper, slist)
        pool.close(); pool.join()
        
        if all([mat[4]==0 for mat in mats]):
            #This is a degenerate mult table
            return None
        
        ret = tuple([ \
            MatrixArgs(QQ, mat[0], mat[1], entries=mat[2]).matrix()/max(1, QQ(mat[4]*mat[3]/(max(1, mat[5])))) \
            for mat in mats])
        
        self._save(ind, ret, True)
        return ret
    
    def sym_asym_bases(self, n, ftype=None):
        r"""
        Generate the change of base matrices for the symmetric
        and the asymmetric subspaces
        """

        flags = self.generate_flags(n, ftype)
        uniques = []
        sym_base = []
        asym_base = []
        for xx in flags:
            xxid = self.identify(n, [xx.ftype_points()], **xx.blocks())
            if xxid not in uniques:
                uniques.append(xxid)
                sym_base.append(xx.afae())
            else:
                sym_ind = uniques.index(xxid)
                asym_base.append(sym_base[sym_ind] - xx.afae())
                sym_base[sym_ind] += xx
        m_sym = matrix(len(sym_base), len(flags), [xx.values() for xx in sym_base], sparse=True)
        m_asym = matrix(len(asym_base), len(flags), [xx.values() for xx in asym_base], sparse=True)
        return m_sym, m_asym
    
    def _density_wrapper(self, ar):
        r"""
        Helper function used in the parallelization of calculating densities
        """
        return ar[0].densities(ar[1], ar[2], ar[3], ar[4], ar[5], ar[6], ar[7])


def _flatten_matrix(mat, doubled=False):
    r"""
    Flatten a symmetric matrix, optionally double non-diagonal elements
    """
    res = []
    factor = 2 if doubled else 1
    for ii in range(len(mat)):
        res.append(mat[ii][ii])
        res += [factor*mat[ii][jj] for jj in range(ii+1, len(mat))]
    return res

def _unflatten_matrix(ls, dim=-1, doubled=False, upper=False):
    r"""
    Unflatten a symmetric matrix, optionally correct for the doubled non-diagonal elements
    """
    if dim==-1:
        dim = Integer(round((1/2) * ((8*len(ls)+1)**(1/2) - 1) ))
    mat = [[0]*dim for ii in range(dim)]
    factor = 2 if doubled else 1
    index = 0
    for ii in range(dim):
        # Fill the diagonal element
        mat[ii][ii] = ls[index]
        index += 1
        # Fill the off-diagonal elements
        for jj in range(ii + 1, dim):
            mat[ii][jj] = ls[index] / factor
            if not upper:
                mat[jj][ii] = ls[index] / factor
            index += 1
    return matrix(mat), ls[index:]

def _round(value, method=1, quotient_bound=7, denom_bound=9, denom=1024):
    r"""
    Helper function, to round a number using either 
    method=0 - simple fixed denominator
    method=1 - continued fractions
    """
    if method==0:
        return QQ(round(value*denom)/denom)
    else:
        from sage.rings.continued_fraction import continued_fraction
        cf = continued_fraction(value)
        for ii, xx in enumerate(cf.quotients()):
            if xx>=2**quotient_bound or cf.denominator(ii)>2**(denom_bound):
                if ii>0:
                    return cf.convergent(ii-1)
                return 0
        return cf.value()

def _round_list(ls, force_pos=False, method=1, quotient_bound=7, denom_bound=9, denom=1024):
    r"""
    Helper function, to round a list
    """
    if force_pos:
        return [max(_round(xx, method, quotient_bound, denom_bound, denom), 0) for xx in ls]
    else:
        return [_round(xx, method, quotient_bound, denom_bound, denom) for xx in ls]

def _round_matrix(mat, method=1, quotient_bound=7, denom_bound=9, denom=1024):
    r"""
    Helper function, to round a matrix
    """
    return matrix(QQ, [_round_list(xx, False, method, quotient_bound, denom_bound, denom) for xx in mat])

def _round_ldl(mat, method=1, quotient_bound=7, denom_bound=9, denom=1024):
    r"""
    Helper function, to round a matrix using ldl decomposition
    """
    mat_ldl = matrix(mat).block_ldlt()
    P = matrix(QQ, mat_ldl[0])
    L = matrix(QQ, _round_matrix(mat_ldl[1], method, quotient_bound, denom_bound, denom))
    D = diagonal_matrix(QQ, _round_list(mat_ldl[2].diagonal(), True, method, quotient_bound, denom_bound, denom))
    pl = P*L
    return pl*D*pl.T

def _round_adaptive(ls, onevec, denom=1024):
    r"""
    Adaptive rounding based on continued fraction and preserving an inner product
    with `onevec`
    
    If the continued fraction rounding fails fall back to a simple denominator method
    """
    best_vec = None
    best_error = 1000
    best_lcm = 1000000000
    
    orig = vector(ls)
    for resol1 in range(5, 20):
        resol2 = round(resol1*1.5)
        rls = vector([_round(xx, quotient_bound=resol1, denom_bound=resol2) for xx in orig])
        ip = rls*onevec
        if ip != 0 and abs(ip - 1)<best_error:
            if ip.as_integer_ratio()[1] > best_lcm**1.5 and ip != 1:
                continue
            best_vec = rls/ip
            best_error = abs(ip - 1)
            best_lcm = ip.as_integer_ratio()[1]
    if best_vec==None:
        rvec = vector(QQ, _round_list(ls, True, method=0, denom=denom))
        best_vec = rvec/(rvec*onevec)
    return best_vec, ((best_vec-orig)/(len(orig)**0.5)).norm()

def _fraction_print(val, thr=20):
    r"""
    Print out a fraction. If it is too big, then print an approximation
    indicating with a ? that it is not an exact value.
    """
    if len(str(val))>thr:
        return str(val.n())+"?"
    else:
        return str(val)

class FlagAlgebraElement(Element):
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
        self._values = vector(base, values, sparse=True)
        self._ftype = parent.ftype()
        Element.__init__(self, parent)
    
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
        return vector(R, self._values, sparse=True)
    
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
            if len(self)<10:
                sttrl.append(('{:<'+str(maxstrlen)+'} - {}').format(strs[ii], str(flgs[ii])))
            else:
                include = True
                try: 
                    include = abs(float(self.values()[ii]))>=1e-8
                except: 
                    include = self.values()[ii]!=0
                if include:
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
        N = self.parent().theory().size_combine(self.ftype().size(), self.size(), other.size())
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
        ressize = amount + self.size()
        table = self.parent().mpt(self.size(), self.ftype().size(), target_size=ressize)
        vals = [sum(self.values() * mat) for mat in table]
        return self.__class__(self.parent(), ressize, vals)
    
    def __getitem__(self, flag):
        if (not isinstance(flag, Flag)) or (not flag.ftype()==self.ftype()) or (not self.size()==flag.size()):
            raise TypeError("Indecies must be Flags with matching ftype and size, not {}".format(str(type(flag))))
        return self.values()[self.flags().index(flag)]
    
    def __setitem__(self, flag, value):
        if (not isinstance(flag, Flag)) or (not flag.ftype()==self.ftype()) or (not self.size()==flag.size()):
            raise TypeError("Indecies must be Flags with matching ftype and size, not {}".format(str(type(flag))))
        self.values()[self.flags().index(flag)] = value
    
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
    
    def mul_project(self, other, ftype_inj=tuple(), target_size=None):
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
        if new_ftype==None or new_ftype.unique()==None:
            raise ValueError("The ftype injection maps to an invalid ftype.")
        N = self.parent().theory().size_combine(self.ftype().size(), self.size(), other.size())
        if target_size!=None:
            if target_size<N:
                raise ValueError("Target size is smaller than minimum allowed size for this operation.")
            N = target_size
        table = self.parent().mpt(self.size(), other.size(), ftype_inj=ftype_inj, target_size=N)
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
    
    def set_sum(self, target=1):
        r"""
        Make the symbolic variables sum to the given target.
        """
        valvec = self.values()
        if len(valvec)==0:
            return self
        par = valvec[0].parent()
        try:
            gs = par.gens()
        except:
            return self
        repl = {gs[-1]: target - sum(gs[:-1])}
        nvec = vector([xx.subs(repl) for xx in valvec])
        return self.parent()(self.size(), nvec)

    def subs(self, args):
        r"""
        Symbolic substitution.
        """
        valvec = self.values()
        if len(valvec)==0:
            return self
        par = valvec[0].parent()
        try:
            gs = par.gens()
        except:
            return self
        repl = {gs[ii]:args[ii] for ii in range(min(len(args), len(gs)))}
        nvec = vector([xx.subs(repl) for xx in valvec])
        retalg = FlagAlgebra(QQ, self.parent().theory())
        return retalg(self.size(), nvec)

    def derivative(self, times):
        r"""
        Symbolic differentiation.
        """
        valvec = self.values()
        if len(valvec)==0:
            return self
        par = valvec[0].parent()
        try:
            gs = par.gens()
        except:
            return self
        rvec = []
        for ff, vv in enumerate(valvec):
            if vv==0:
                rvec.append(0)
                continue
            aval = vv
            for ii in range(min(len(times), len(gs))):
                aval = aval.derivative(gs[ii], times[ii])
            rvec.append(aval)
        return self.parent()(self.size(), vector(rvec))

    def derivatives(self, point):
        r"""
        Returns all symbolic derivatives evaluated at a given point.
        """
        times = []
        valvec = self.values()
        if len(valvec)==0:
            return self
        par = valvec[0].parent()
        try:
            gs = par.gens()
        except:
            return self
        for xx in self.values():
            if xx==0:
                continue
            dd = xx.dict()
            for kk in dd.keys():
                kk = tuple(kk)
                if kk in times:
                    continue
                for ll in itertools.product(*[range(ii+1) for ii in kk]):
                    if ll not in times:
                        times.append(ll)
        res = []
        for xx in times:
            der = self.derivative(xx).subs(point)
            minnz = 1000000
            for xx in der.values():
                if int(xx)!=0:
                    minnz = min(abs(xx), minnz)
            if minnz != 1000000:
                res.append(der/minnz)
        return res
    
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

class FlagAlgebra(Parent, UniqueRepresentation):
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
        Parent.__init__(self, base)
    
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
                    vec = vector(base, [(1 if xx==v else 0) for xx in flags], sparse=True)
                    return self.element_class(self, v.size(), vec)
            elif isinstance(v, FlagAlgebraElement):
                if v.ftype()==self.ftype():
                    if self.base()==v.parent().base():
                        return v
                    elif self.base().has_coerce_map_from(v.parent().base()):
                        vals = vector(self.base(), v.values(), sparse=True)
                        return self.element_class(self, v.size(), vals)
            elif v in base:
                return self.element_class(self, self.ftype().size(), vector(base, [v], sparse=True))
            raise ValueError('Can\'t construct an element from {}'.format(v))
        return self.element_class(self, *args, **kwds)
    
    def sizes(self):
        return self.theory().sizes()
    
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
        f = self.combinatorial_theory()._an_element_(n=self.ftype().size(), ftype=self.ftype())
        return self(f)*a
    
    def some_elements(self):
        r"""
        Returns a small list of elements
        """
        return [self.an_element(),self(self.base().an_element())]
    
    def mul_project_table(self, n1, n2, ftype_inj=None, target_size=None):
        r"""
        Returns the multiplication projection table
        
        This is the same as :func:`CombinatorialTheory.mul_project_table`
        with self.ftype() substituted in.

        .. SEEALSO::

            :func:`CombinatorialTheory.mul_project_table`
        """
        return self.theory().mul_project_table(n1, n2, self.ftype(), ftype_inj, target_size)
    
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
    ftype_union = [jj for ff in ftype_points for jj in ff]
    partition = list(ftype_points) + [[ii for ii in range(n) if ii not in ftype_union]]
    g = Graph([list(range(n)), edges], format='vertices_and_edges')
    blocks = tuple(g.canonical_label(partition=partition).edges(labels=None, sort=True))
    return (n, tuple([len(xx) for xx in ftype_points]), blocks)

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
    ftype_union = [jj for ff in ftype_points for jj in ff]
    partt = list(ftype_points) + \
            [[ii for ii in range(n) if ii not in ftype_union]] + \
            [list(range(n,n+len(edges)))]
    blocks = tuple(g.canonical_label(partition=partt).edges(labels=None, sort=True))
    return (n, tuple([len(xx) for xx in ftype_points]), blocks)

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
    partition = list(ftype_points) + [list(set(range(n)).difference(set(ftype_points))), ]
    g = DiGraph([list(range(n)), edges], format='vertices_and_edges')
    blocks = tuple(g.canonical_label(partition=partition).edges(labels=None, sort=True))
    return (n, tuple([len(xx) for xx in ftype_points]), blocks)

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
    ftype_union = [jj for ff in ftype_points for jj in ff]
    partt = list(ftype_points) + \
            [[ii] for ii in range(n, n+len(edges))] + \
            [[ii for ii in range(n) if ii not in ftype_union]]
    blocks = tuple(g.canonical_label(partition=partt).edges(labels=None, sort=True))
    return (n, tuple([len(xx) for xx in ftype_points]), blocks)

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
    ftype_union = [jj for ff in ftype_points for jj in ff]
    if len(edges_marked) == len(edges)/2:
        edges_marked_alter = [xx for xx in edges if xx not in edges_marked]
        g = Graph([list(range(n+len(edges)+len(edges_marked))), 
                   [(i+n, x) for i,b in enumerate(edges) for x in b] + 
                  [(i+n+len(edges), x) for i, b in enumerate(edges_marked) for x in b]], 
                  format='vertices_and_edges')
        partt = list(ftype_points) + \
                [[ii for ii in range(n) if ii not in ftype_union]] + \
                [list(range(n,n+len(edges)))] + \
                [list(range(n+len(edges), n+len(edges)+len(edges_marked)))]
        blocks = tuple(g.canonical_label(partition=partt).edges(labels=None, sort=True))
        
        g_alter = Graph([list(range(n+len(edges)+len(edges_marked_alter))), 
                   [(i+n, x) for i,b in enumerate(edges) for x in b] + 
                  [(i+n+len(edges), x) for i, b in enumerate(edges_marked_alter) for x in b]], 
                  format='vertices_and_edges')
        partt_alter = list(ftype_points) + \
                [[ii for ii in range(n) if ii not in ftype_union]] + \
                [list(range(n,n+len(edges)))] + \
                [list(range(n+len(edges), n+len(edges)+len(edges_marked_alter)))]
        blocks_alter = tuple(g_alter.canonical_label(partition=partt_alter).edges(labels=None, sort=True))
        return (n, tuple([len(xx) for xx in ftype_points]), min(blocks, blocks_alter))
    else:
        if len(edges_marked) > len(edges)/2:
            edges_marked = [xx for xx in edges if xx not in edges_marked]
        g = Graph([list(range(n+len(edges)+len(edges_marked))), 
                   [(i+n, x) for i,b in enumerate(edges) for x in b] + 
                  [(i+n+len(edges), x) for i, b in enumerate(edges_marked) for x in b]], 
                  format='vertices_and_edges')
        partt = list(ftype_points) + \
                [[ii for ii in range(n) if ii not in ftype_union]] + \
                [list(range(n,n+len(edges)))] + \
                [list(range(n+len(edges), n+len(edges)+len(edges_marked)))]
        blocks = tuple(g.canonical_label(partition=partt).edges(labels=None, sort=True))
        return (n, tuple([len(xx) for xx in ftype_points]), blocks)

def _cube_graphs(d, edge_num):
    from sage.features.nauty import NautyExecutable
    import subprocess, select
    import shlex
    directg_path = NautyExecutable("multig").absolute_filename()
    dcube = graphs.CubeGraph(d, embedding=0)
    le = len(dcube.edges(labels=None))
    options=" -T -m2 -e{}".format(le+edge_num)
    sub = subprocess.Popen(
        shlex.quote(directg_path) + ' {0}'.format(options),
        shell=True,
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='latin-1'
    )
    sub.stdin.write(dcube.graph6_string())
    sub.stdin.close()

    for line in sub.stdout:
        if line and line[0]==str(2**d) or line[:2]==str(2**d):
            edges = []
            seq = line[:-1].split(" ")
            for ii in range(1, len(seq)//3):
                try:
                    ed = (int(seq[ii*3]), int(seq[ii*3 + 1]))
                except:
                    print("error on line: \n", line)
                    return
                edges.append(ed)
                if seq[ii*3 + 2]=="2":
                    edges.append((ed[0], ed[1]))
            yield {'edges': edges}
    for line in sub.stderr:
        pass
    sub.wait()

def _generator_cube_graphs(n):
    if n==0:
        yield {'edges': []}
        return
    if log(n, 2) in NN:
        d = log(n, 2)
        for enum in range(2**(n-1)*n + 1):
            for xx in _cube_graphs(d, enum):
                yield xx
    else:
        return None

def _identify_cube_graphs(n, ftype_points, edges):
    if n==0:
        return (), ()
    if not (log(n, 2) in NN):
        return None
    d = log(n, 2)
    if len(set(edges))!=2**(d-1) * d:
        return None
    ftype_union = [jj for ff in ftype_points for jj in ff]
    partt = list(ftype_points) + [[ff for ff in range(n) if ff not in ftype_union]]
    g = Graph([list(range(n)), edges], format='vertices_and_edges', multiedges=True, vertex_labels=False)
    blocks = g.canonical_label(partition=partt).edges(labels=False, sort=True)
    return tuple(blocks), tuple([len(xx) for xx in ftype_points])

def _cube_points(d, point_num, edges):
    unique = []
    n = 2**d
    for ps in itertools.combinations(range(n), point_num):
        points = [[ii] for ii in ps]
        id = _identify_cube_points(n, [], edges, points)
        if id not in unique:
            unique.append(id)
            yield {'edges': edges, 'points': points}

def _generator_cube_points(n):
    if n==0:
        yield {'edges':[], 'points':[]}
        return
    if log(n, 2) in NN:
        d = log(n, 2)
        dcube = graphs.CubeGraph(d, embedding=0)
        edges = Graph(dcube.graph6_string(), format="graph6").edges(labels=None)
        for pnum in range(n+1):
            for xx in _cube_points(d, pnum, edges):
                yield xx
    else:
        return None

def _identify_cube_points(n, ftype_points, edges, points):
    if n==0:
        return (), ()
    if not (log(n, 2) in NN):
        return None
    d = log(n, 2)
    if len(set(edges))!=2**(d-1) * d:
        return None
    ftype_union = [jj for ff in ftype_points for jj in ff]
    g_parts = list(ftype_points) + \
              [[ii for ii in range(n) if ii not in ftype_union]] + [[n]]
    g_verts = list(range(n+1))
    g_edges = list(edges) + [(ii[0], n) for ii in points]
    g = Graph([g_verts, g_edges], format='vertices_and_edges')
    blocks = tuple(g.canonical_label(partition=g_parts).edges(labels=None, sort=True))
    return (n, tuple([len(xx) for xx in ftype_points]), blocks)

def _cube_size_combine(k, n1, n2):
    if k==0 and n1==0 and n2==0:
        return 0
    if k==0:
        if n1>0 and n2>0:
            return None
    n1 = max(n1, 1)
    n2 = max(n2, 1)
    k = max(k, 1)
    if log(n1, 2) not in NN or log(n2, 2) not in NN or log(k, 2) not in NN:
        return None
    return QQ(n1*n2/k)

def colored_identify(k, order_partition, n, ftype_points, **kwargs):
    r"""
    A general identifier code, it works on any edge arity and color number,
    and supports color symmetries using the `order_partition`
    """
    is_graph = (k==2)
    color_number = sum(len(xx) for xx in order_partition)
    edges = kwargs["edges"]
    Cs = [[cx[0] for cx in kwargs["C{}".format(ii)]] for ii in range(color_number)]
    ftype_union = [jj for ff in ftype_points for jj in ff]
    g_parts = list(ftype_points) + \
              [[ii for ii in range(n) if ii not in ftype_union]]
    ppadd = 0 if is_graph else len(edges)
    g_verts = list(range(n+ppadd+color_number))
    g_parts.append(list(range(n, n+ppadd)))

    g_parts += [[n+ppadd+ii for ii in partition_j] for partition_j in order_partition]
    
    if is_graph:
        g_edges = list(edges)
        for ii in range(color_number):
            g_edges += [(xx, n+ii) for xx in Cs[ii]]
    else:
        g_edges = [(i+n,x) for i,b in enumerate(edges) for x in b]
        for ii in range(color_number):
            g_edges += [(xx, n+len(edges)+ii) for xx in Cs[ii]]
    g = Graph([g_verts, g_edges], format='vertices_and_edges')
    blocks = tuple(g.canonical_label(partition=g_parts).edges(labels=None, sort=True))
    return (n, tuple([len(xx) for xx in ftype_points]), blocks)

def colored_generate(k, order_partition, n):
    r"""
    A general generator code, it works on any edge arity and color number,
    and supports color symmetries using the `order_partition`
    """
    color_number = sum(len(xx) for xx in order_partition)
    if k==2:
        BT = GraphTheory
    if k==3:
        BT = ThreeGraphTheory
    for xx in BT.generate_flags(n):
        unique = []
        edges = xx.blocks()['edges']
        
        for yy in itertools.product(range(color_number), repeat=int(n)):
            yy = list(yy)
            Cs = {"C{}".format(cc):[[ii] for ii, oo in enumerate(yy) if oo==cc] for cc in range(color_number)}
            iden = colored_identify(k==2, order_partition, n, [], edges=edges, **Cs)
            if iden not in unique:
                unique.append(iden)
                Cs["edges"] = edges
                yield Cs

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
#DiGraphTheory.exclude([DiGraphTheory(2, edges=[[0, 1], [1, 0]]), DiGraphTheory(2)])

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

# should only be used up to size 8.
# note that for the next size, even with only 15 edges (out of the 32) there are 1.4M nonisomorphic graphs
# len(list(_cube_graphs(4, 15)))
# 1479300
HypercubeGraphTheory = CombinatorialTheory('HypercubeGraph', 
                                           _generator_cube_graphs, 
                                           _identify_cube_graphs, 
                                           size_combine=_cube_size_combine, 
                                           edges=2)

# should only be used up to size 16.
HypercubeVertexTheory = CombinatorialTheory('HypercubeVertex', 
                                         _generator_cube_points, 
                                         _identify_cube_points, 
                                         size_combine=_cube_size_combine, 
                                         edges=2, points=1)