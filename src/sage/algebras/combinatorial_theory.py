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

    sage: GraphTheory.generate(4)
    (Flag on 4 points, ftype from () with edges=(),
     Flag on 4 points, ftype from () with edges=(01),
     Flag on 4 points, ftype from () with edges=(01 03),
     Flag on 4 points, ftype from () with edges=(02 13),
     Flag on 4 points, ftype from () with edges=(01 02 03),
     Flag on 4 points, ftype from () with edges=(01 02 13),
     Flag on 4 points, ftype from () with edges=(02 03 12 13))

Calling :func:`exclude` adds new structures to the list of flags we don't want
to generate, and we don't want to appear in larger strctures as an induced
subflag. To reset the theory and don't exclude anything, call :func:`reset`.

To optimize the density of a flag (or linear combination of flags) in a theory,
we can call :func:`optimize`. To try to find the maximum number of 
edges `e` in `k3` free graphs we can write ::
    
    sage: x = GraphTheory.optimize(e, 3)
    ...
    Success: SDP solved
    ...
    sage: abs(x-0.5)<1e-6
    True
    
The second parameter, `optimize_problem(e, 3)` indicates the maximum size the
program expands the flags. We can reset the excluded graphs, and try to minimize 
the density of triangles and empty triples ::
    
    sage: GraphTheory.reset()
    sage: e3 = GraphTheory(3)
    sage: x = GraphTheory.optimize(e3+k3, 3, maximize=False)
    ...
    sage: abs(x-0.25)<1e-6
    True

The :func:`optimize` function requires csdpy, an sdp solver. But it is possible
to get these relations directly, by expressing the inequalities 
as a sum of squares ::

    sage: GraphTheory.reset()
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

- Levente Bodnar (2023-2025): Main development

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

import itertools

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.all import QQ, NN, Integer, infinity
from sage.algebras.flag import Flag, Pattern, inductive_generator, overlap_generator
from sage.algebras.flag_algebras import FlagAlgebra, FlagAlgebraElement

from sage.categories.sets_cat import Sets
from sage.all import vector, matrix, diagonal_matrix

from sage.misc.prandom import randint
from sage.arith.misc import falling_factorial, binomial, factorial
from sage.misc.functional import round
from functools import lru_cache

import hashlib
import pickle
import os
from tqdm import tqdm

def combine(name, *theories, symmetries=False):
    if not isinstance(name, str):
        raise ValueError("Name must be a string")
    
    #Sanity checks
    if len(theories)==0:
        raise ValueError("At least one theory is expected!")
    if len(theories)==1:
        import warnings
        warnings.warn("Warning, only one theory was provided." + 
                      "It will be returned with the same name.")
        return theories[0]

    #Check if we can use symmetry, and the resulting groups
    can_symmetry = True
    next_group = 0
    result_signature = {}
    result_symmetry = []
    result_excluded = []

    for theory in theories:
        next_group_increment = 0
        if len(theory._signature.keys())!=1:
            can_symmetry = False
        for kk in theory._signature:
            if kk in result_signature:
                raise ValueError("The relation names must be different!")
            tkk = dict(theory._signature[kk])
            if can_symmetry:
                for ll in result_signature:
                    tll = result_signature[ll]
                    if tll["arity"]!=tkk["arity"] or \
                    tll["ordered"]!=tkk["ordered"]:
                        can_symmetry = False
            next_group_increment = max(next_group_increment, tkk["group"]+1)
            tkk["group"] += next_group
            result_signature[kk] = tkk
        result_symmetry += list(theory._symmetries)
        result_excluded += list(theory._excluded)
        next_group += next_group_increment

    if not can_symmetry:
        if len(theories)!=2:
            raise ValueError("Can't combine more than 2 theories with " + 
                             "different parameters.")
        if symmetries is not False:
            import warnings
            warnings.warn("Combined theories have different parameters, " + 
                          "symmetries will be ignored.")
            symmetries = False

    if symmetries is not False:
        #Make everything in the same group
        for xx in result_signature:
            result_signature[xx]["group"] = 0
        if symmetries is True:
            #This case symmetry is trivial for the entire group
            result_symmetry = [(len(theories), len(theories), tuple())]
        else:
            #This case symmetry is as provided by the parameter
            m = 0
            formatted_sym = []
            for edge in symmetries:
                m = max(m, edge[0], edge[1])
                formatted_sym.append(tuple(sorted(list(edge))))
            formatted_sym = tuple(sorted(formatted_sym))
            result_symmetry = [(len(theories), m+1, formatted_sym)]
    #Note that otherwise we use each symmetry from the combined pieces
    theory_data = {
        "signature": result_signature, 
        "symmetries": result_symmetry
        }
    ser_data = _serialize_data(theory_data)
    if len(theories)==2 and not symmetries:
        ser_data = (tuple(theories), ser_data)
    else:
        ser_data = (None, ser_data)
    ret_theory = CombinatorialTheory(name, _from_data=ser_data)
    return ret_theory

def _serialize_data(data):
    #For the signature
    signature = data["signature"]
    sered_signature = []
    for xx in signature:
        ll = tuple(signature[xx].values())
        sered_signature.append((xx, ll))
    sered_signature = tuple(sered_signature)
    return (sered_signature, tuple(data["symmetries"]))

def test_generate():
    def test_theory(TT, nstart, nend, vals):
        print("\nTesting {}".format(str(TT)))
        for ii,jj in enumerate(range(nstart, nend+1)):
            print("Size {}, the number is {} (should be {})".format(
                jj, 
                len(TT.generate(jj)), 
                vals[ii]
                ))
    
    CG = combine("CGraph", Color0, GraphTheory)
    Cs = combine("Cs", Color0, Color1, symmetries=True)
    CG.clear()
    Cs.clear()
    test_theory(Color0, 5, 10, [6, 7, 8, 9, 10, 11])
    test_theory(Cs, 3, 8, [13, 22, 34, 50, 70, 95])
    test_theory(GraphTheory, 3, 7, [4, 11, 34, 156, 1044])
    test_theory(ThreeGraphTheory, 3, 6, [2, 5, 34, 2136])
    test_theory(DiGraphTheory, 2, 5, [3, 16, 218, 9608])
    test_theory(DiThreeGraphTheory, 3, 3, [16])
    test_theory(CG, 2, 6, [6, 20, 90, 544, 5096])
    Cs.exclude([Cs(1), Cs(1, C0=[[0]], C1=[[0]])])
    GraphTheory.exclude(GraphTheory(3))
    CGp = combine("CGsym", GraphTheory, Cs)
    CGp.clear()
    test_theory(Cs, 4, 8, [3, 3, 4, 4, 5])
    test_theory(GraphTheory, 3, 8, [3, 7, 14, 38, 107, 410])
    test_theory(CGp, 2, 5, [4, 8, 32, 106])
    Css = combine("Colors3Sym", Color0, Color1, Color2, symmetries=True)
    Css.clear()
    pe = Css(1)
    p0 = Css.p(1, C2=[0], C1=[0])
    Css.exclude([pe, p0])
    test_theory(Css, 3, 8, [3, 4, 5, 7, 8, 10])
    Cyc4 = combine("Cyclic4", Color0, Color1, Color2, Color3, 
                   symmetries=CyclicSymmetry(4))
    Cyc4.clear()
    Cyc4.exclude([
        Cyc4(1),
        Cyc4.p(1, C0=[0], C1=[0]),
        Cyc4.p(1, C0=[0], C2=[0])
    ])
    GraphTheory.reset()
    T4 = combine("Cyclic4Graph", Cyc4, GraphTheory)
    Cyc6 = combine("Cyclic6", Color0, Color1, Color2, Color3, Color4, Color5, 
                   symmetries=CyclicSymmetry(6))
    Cyc6.clear()
    Cyc6.exclude([
        Cyc6(1),
        Cyc6.p(1, C0=[0], C1=[0]),
        Cyc6.p(1, C0=[0], C2=[0]),
        Cyc6.p(1, C0=[0], C3=[0])
    ])
    T6 = combine("Cyclic6Graph", Cyc6, GraphTheory)
    T4.clear()
    T6.clear()
    test_theory(Cyc4, 4, 9, [10, 14, 22, 30, 43, 55])
    test_theory(Cyc6, 3, 7, [10, 22, 42, 80, 132])
    test_theory(T4, 2, 5, [6, 30, 260, 3052])
    test_theory(T6, 2, 4, [8, 62, 754])

def clear_all_calculations(theory_name=None):
    calcs_dir = os.path.join(os.getenv('HOME'), '.sage', 'calcs')
    if not os.path.exists(calcs_dir):
        return
    for xx in os.listdir(calcs_dir):
        if theory_name==None or str(xx).startswith(theory_name):
            file_path = os.path.join(calcs_dir, xx)
            os.remove(file_path)

def show_all_calculations(theory_name=None):
    calcs_dir = os.path.join(os.getenv('HOME'), '.sage', 'calcs')
    if not os.path.exists(calcs_dir):
        return
    for xx in os.listdir(calcs_dir):
        file_path = os.path.join(calcs_dir, xx)
        file_theory = str(xx).split(".")[0]
        if theory_name==None:
            with open(file_path , "rb") as file:
                data = pickle.load(file)
            if data != None:
                print(file_theory + ": " + str(data["key"][:2]))
        elif str(xx).startswith(theory_name):
            with open(file_path , "rb") as file:
                data = pickle.load(file)
            if data != None:
                print(data["key"][:2])

class CombinatorialTheory(Parent, UniqueRepresentation):
    
    Element = Flag
    
    def __init__(self, name, relation_name="edges", arity=2, 
                 is_ordered=False, _from_data=None):
        r"""
        Initialize a Combinatorial Theory
        
        A combinatorial theory is any theory with universal axioms only, 
        (therefore the elements satisfy a heredetary property).
        See the file docstring for more information.

        INPUT:

        - ``name`` -- string; name of the Theory
        - ``relation_name`` -- string; name of the relation
        - ``arity`` -- integer; arity of the relation
        - ``is_ordered`` -- boolean; if the values are ordered
        - ``_from_data`` -- list; only used internally

        OUTPUT: A CombinatorialTheory object
        """
        self._name = name
        
        if _from_data != None:
            self._sources = _from_data[0]
            _from_data = _from_data[1]
            
            sered_signature = _from_data[0]
            self._signature = {}
            max_group = -1
            for ll in sered_signature:
                key = ll[0]
                val = {
                    "arity": ll[1][0],
                    "ordered": ll[1][1],
                    "group": ll[1][2]
                }
                max_group = max(max_group, val["group"])
                self._signature[key] = val
            self._symmetries = _from_data[1]
            if len(self._symmetries) != max_group+1:
                print(self._symmetries)
                print(self._signature)
                raise ValueError("Provided data has different symmetry " + 
                                 "set size than group number")
        else:
            if arity < 1 or (arity not in NN):
                raise ValueError("Arity must be nonzero positive integer!")
            self._signature = {relation_name: {
                "arity": arity,
                "ordered": is_ordered,
                "group": 0
            }}
            self._sources = None
            self._symmetries = ((1, 1, tuple()), )
        self._excluded = tuple()
        self._no_question = False
        Parent.__init__(self, category=(Sets(), ))
        self._populate_coercion_lists_()
    
    def _repr_(self):
        r"""
        Give a nice string representation of the theory object

        OUTPUT: The representation as a string

        EXAMPLES::

            sage: print(GraphTheory)
            Theory for Graph
        """
        return 'Theory for {}'.format(self._name)

    def signature(self):
        r"""
        Returns the signature data for this theory

        OUTPUT: A dictionary containing the signature data

        EXAMPLES::

            sage: GraphTheory.signature()
            {'edges': {'arity': 2, 'group': 0, 'ordered': False}}
        """
        return self._signature

    def symmetries(self):
        r"""
        Returns the symmetry data for this theory

        OUTPUT: A tuple containing the signature data

        EXAMPLES::

            sage: GraphTheory.symmetries()
            ((1, 1, ()),)
        """
        return self._symmetries

    #Parent methods
    def _element_constructor_(self, n, **kwds):
        r"""
        Construct elements of this theory

        INPUT:

        - ``n`` -- the size of the flag
        - ``**kwds`` -- can contain ftype_points, listing
            the points that will form part of the ftype;
            and can contain the blocks for each signature.
            If they are not included, they are assumed to 
            be empty lists.

        OUTPUT: A Flag with the given parameters

        EXAMPLES::

        Create an empty graph on 3 vertices ::

            sage: GraphTheory(3)
            Flag on 3 points, ftype from () with edges=()
        
        Create an edge with one point marked as an ftype ::
        
            sage: GraphTheory(2, ftype_points=[0], edges=[[0, 1]])
            Flag on 2 points, ftype from (0,) with edges=(01)

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

        if isinstance(n, Flag) or isinstance(n, Pattern):
            if n.parent()==self:
                return n
            n = n.as_pattern()
            return self.pattern(n.size(), ftype=n.ftype_points(), 
                                **n.as_pattern().blocks())

        ftype_points = tuple()
        if 'ftype_points' in kwds:
            try:
                ftype_points = tuple(kwds['ftype_points'])
            except:
                raise ValueError("The provided ftype_points must be iterable")
        elif 'ftype' in kwds:
            try:
                ftype_points = tuple(kwds['ftype'])
            except:
                raise ValueError("The provided ftype must be iterable")
        
        blocks = {}
        for xx in self._signature.keys():
            blocks[xx] = tuple()
            unary = (self._signature[xx]["arity"]==1)
            if xx in kwds:
                try:
                    blocks[xx] = tuple(kwds[xx])
                except:
                    raise ValueError("The provided {} must be iterable".format(xx))
                if unary:
                    if len(blocks[xx])>0:
                        try:
                            tuple(blocks[xx][0])
                        except:
                            blocks[xx] = tuple([[aa] for aa in blocks[xx]])
                        
        return self.element_class(self, n, ftype_points, **blocks)
    
    def empty_element(self):
        r"""
        Returns the empty element, ``n``=0 and no blocks

        OUTPUT: The empty element of the CombinatorialTheory

        EXAMPLES::

            sage: GraphTheory.empty_element()
            Ftype on 0 points with edges=()

        .. NOTE::
            This has an alias called :func:`empty`
            Since the underlying vertex set (empty set)
            is the same as the ftype point set, this is
            an ftype

        .. SEEALSO::

            :func:`empty`
        """
        blocks = {}
        for xx in self._signature:
            blocks[xx] = tuple()
        return self.element_class(self, 0, tuple(), **blocks)
    
    empty = empty_element
    
    def pattern(self, n, **kwds):
        r"""
        Construct patterns for this theory

        INPUT:

        - ``n`` -- the size of the flag
        - ``**kwds`` -- can contain ftype_points, listing
            the points that will form part of the ftype;
            and can contain the blocks for each signature.
            If they are not included, they are assumed to 
            be empty lists. Can also contain missing relations
            for each signature entry.

        OUTPUT: A Pattern with the given parameters

        EXAMPLES::

        Create a pattern on 3 vertices with one edge required
        and one edge missing ::

            sage: GraphTheory(3, edges=[[0, 1]], edges_m=[[1, 2]])
            Flag on 3 points, ftype from () with edges=(01)

        .. NOTE::
            Also has alias :func:`P`, :func:`p`, :func:`Pattern`

        .. SEEALSO::

            :func:`__init__` of :class:`Pattern`
        """
        ftype_points = tuple()
        if 'ftype_points' in kwds:
            try:
                ftype_points = tuple(kwds['ftype_points'])
            except:
                raise ValueError("The provided ftype_points must be iterable")
        elif 'ftype' in kwds:
            try:
                ftype_points = tuple(kwds['ftype'])
            except:
                raise ValueError("The provided ftype must be iterable")
        if len(ftype_points)==n:
            return self._element_constructor_(n, **kwds)
        
        blocks = {}
        for xx in self._signature.keys():
            blocks[xx] = tuple()
            blocks[xx+"_m"] = tuple()
            unary = self._signature[xx]["arity"]==1


            if xx in kwds:
                try:
                    blocks[xx] = tuple(kwds[xx])
                except:
                    raise ValueError(
                        "The provided {} must be iterable".format(xx)
                        )
                if unary:
                    if len(blocks[xx])>0:
                        try:
                            tuple(blocks[xx][0])
                        except:
                            blocks[xx] = tuple([[aa] for aa in blocks[xx]])
            
            for xx_missing in [xx+"_m", xx+"_missing", xx+"_miss"]:
                if xx_missing in kwds:
                    blocks[xx+"_m"] = kwds[xx_missing]


                    try:
                        blocks[xx+"_m"] = tuple(kwds[xx_missing])
                    except:
                        raise ValueError(
                            "The provided {} must be iterable".format(xx_missing)
                            )
                    if unary:
                        if len(blocks[xx])>0:
                            try:
                                tuple(blocks[xx][0])
                            except:
                                blocks[xx+"_m"] = tuple(
                                    [[aa] for aa in blocks[xx+"_m"]]
                                    )
        return Pattern(self, n, ftype_points, **blocks)
    
    p = pattern
    P = pattern
    Pattern = pattern

    def _an_element_(self, n=0, ftype=None):
        r"""
        Returns a random element

        INPUT:

        - ``n`` -- integer (default: `0`); size of the element
        - ``ftype`` -- Flag (default: `None`); ftype of the element
            if not provided then returns an element with empty ftype

        OUTPUT: A Flag with matching parameters

        EXAMPLES::

            sage: GraphTheory._an_element_()
            Ftype on 0 points with edges=()
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

        EXAMPLES::

            sage: GraphTheory.some_elements()
            [Ftype on 0 points with edges=()]
        """
        res = [self._an_element_()]
        return res

    #Persistend data management
    def _calcs_dir(self):
        r"""
        Returns the path where the calculations are stored.

        EXAMPLES::

            sage: GraphTheory._calcs_dir()
            '/home/bodnalev/.sage/calcs'
        """
        calcs_dir = os.path.join(os.getenv('HOME'), '.sage', 'calcs')
        if not os.path.exists(calcs_dir):
            os.makedirs(calcs_dir)
        return calcs_dir

    def _save(self, data, key=None, path=None, name=None):
        r"""
        Saves a calculation to persistent memory.
        """
        if name==None:
            if key==None:
                raise ValueError(
                    "Either the key or the name must be provided!"
                    )
            serialized_key = pickle.dumps((self, key))
            hashed_key = hashlib.sha256(serialized_key).hexdigest()
            file_name = self._name + "." + hashed_key
        else:
            file_name = name

        if path==None:
            file_path = os.path.join(self._calcs_dir(), file_name)
        elif path=="":
            file_path = file_name
        else:
            if not os.path.exists(path):
                os.makedirs(path)
            file_path = os.path.join(path, file_name)
        save_object = {'key': key, 'data': data}
        with open(file_path, "wb") as file:
            pickle.dump(save_object, file)

    def _load(self, key=None, path=None, name=None):
        r"""
        Tries to load a calculation from persistent memory.
        """
        if key!=None:
            serialized_key = pickle.dumps((self, key))
            hashed_key = hashlib.sha256(serialized_key).hexdigest()
            file_name = self._name + "." + hashed_key
        if name!=None:
            file_name = name

        if path==None:
            file_path = os.path.join(self._calcs_dir(), file_name)
        else:
            if not os.path.exists(path):
                os.makedirs(path)
            file_path = os.path.join(path, file_name)
        
        if not os.path.exists(file_path):
            return None

        with open(file_path, "rb") as file:
            save_object = pickle.load(file)
        
        if key!=None and save_object != None and save_object['key'] != key:
            import warnings
            warnings.warn("Hash collision or corrupted data!")
            return None
        
        return save_object['data']

    def show_files(self):
        r"""
        Shows the persistent files saved from this theory.
        """
        for xx in os.listdir(self._calcs_dir()):
            if xx.startswith(self._name + "."):
                file_path = os.path.join(self._calcs_dir(), xx)
                with open(file_path , "rb") as file:
                    data = pickle.load(file)
                if data != None:
                    print(data["key"][:2])

    def clear(self):
        r"""
        Clears all calculation from the persistent memory.
        """
        for xx in os.listdir(self._calcs_dir()):
            if xx.startswith(self._name + "."):
                file_path = os.path.join(self._calcs_dir(), xx)
                os.remove(file_path)

    def _serialize(self, excluded=None):
        r"""
        Serializes this theory. Note this contains information about 
        the structures excluded from this theory.

        EXAMPLES::

            sage: GraphTheory._serialize()
            {'excluded': ((3, (), ((0, 1), (0, 2), (1, 2))),),
             'name': 'Graph',
             'signature': {'edges': {'arity': 2, 'group': 0, 'ordered': False}},
             'sources': None,
             'symmetries': ((1, 1, ()),)}
        """
        if excluded==None:
            excluded = self.get_total_excluded(100000)
        else:
            excluded = tuple(excluded)
        sourceser = None
        if self._sources != None:
            sourceser = (
                self._sources[0]._serialize(),
                self._sources[1]._serialize()
            )
        return {
            "name": self._name,
            "signature": self._signature,
            "symmetries": self._symmetries,
            "sources": sourceser,
            "excluded": tuple([xx._serialize() for xx in excluded])
        }

    #Optimizing and rounding

    def blowup_construction(self, target_size, pattern_size, 
                            symbolic=False, symmetric=True, unordered=False, 
                            **kwargs):
        r"""
        Returns a blowup construction, based on a given pattern

        INPUT:

        - ``target_size`` -- integer; size of the resulting FlagAlgebraElement
        - ``pattern_size`` -- integer; size of the pattern of the blowup
        - ``symbolic`` -- boolean (default: `False`); if the resulting 
            construction has part sizes symbolic variables
        - ``symmetric`` -- boolean (default: `True`); if the construction is
            symmetric. Speeds up calculation
        - ``unordered`` -- boolean (default: `False`); if the construction's 
            parts are unordered. Slows down calculation
        - ``**kwargs`` -- the parameters of the pattern, one for each signature
            element.

        OUTPUT: A FlagAlgebraElement with values corresponding to the one 
            resulting from a blowup construction
        
        EXAMPLES::

            sage: GraphTheory.blowup_construction(3, 2, edges=[[0, 1]])
            ...
            Flag Algebra Element over Rational Field
            1/4 - Flag on 3 points, ftype from () with edges=()
            0   - Flag on 3 points, ftype from () with edges=(01)
            3/4 - Flag on 3 points, ftype from () with edges=(01 02)
            
            sage: GraphTheory.reset()
            sage: GraphTheory.blowup_construction(3, 2, edges=[[0, 1], [1, 1]], symbolic=True)
            ...
            Flag Algebra Element over Multivariate Polynomial Ring in X0, X1 over Rational Field
            X0^3             - Flag on 3 points, ftype from () with edges=()
            0                - Flag on 3 points, ftype from () with edges=(01)
            3*X0^2*X1        - Flag on 3 points, ftype from () with edges=(01 02)
            3*X0*X1^2 + X1^3 - Flag on 3 points, ftype from () with edges=(01 02 12)
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
                        clusters = [
                            [ii for ii in range(target_size) if verts[ii]==ee] \
                                for ee in edge
                            ]
                        bladd += list( \
                            set([tuple(sorted(xx)) \
                            for xx in itertools.product(*clusters) \
                            if len(set(xx))==len(edge)]) \
                                    )
                    blocks[rel] = bladd
                try:
                    res += self(target_size, **blocks).afae()*coeff
                except:
                    raise ValueError(
                        "The construction contains excluded structures: " + 
                        self(target_size, **blocks)
                        )
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
                        clusters = [
                            [ii for ii in range(target_size) if verts[ii]==ee] \
                            for ee in edge
                            ]
                        bladd += list( \
                            set([tuple(sorted(xx)) \
                            for xx in itertools.product(*clusters) \
                            if len(set(xx))==len(edge)]) \
                                    )
                    blocks[rel] = bladd
                try:
                    res += self(target_size, **blocks).afae() * coeff
                except:
                    raise ValueError(
                        "The construction contains excluded structures: " + 
                        self(target_size, **blocks)
                        )
        return res

    def _adjust_table_phi(self, table_constructor, phi_vectors_exact, 
                          test=False, ring=QQ):
        r"""
        Helper to modify a table constructor, incorporating extra data from
        constructions (phi_vectors_exact)
        """
        if len(phi_vectors_exact)==0:
            return table_constructor

        for param in table_constructor.keys():
            ns, ftype, target_size = param
            table = self.mul_project_table(ns, ns, ftype, ftype_inj=[], 
                                           target_size=target_size)
            Zs = [
                [None for _ in range(len(phi_vectors_exact))] \
                    for _ in range(len(table_constructor[param]))
                ]
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
                        print("Construction based Z matrix for " + 
                              "{} is not semidef: {}".format(
                                  ftype, min(Zjj.eigenvalues())
                                  ))
                    if Z==None:
                        Z = Zjj
                    else:
                        Z.augment(Zjj)
                Zk = Z.kernel()
                Zkern = Zk.basis_matrix()
                if Zkern.nrows()>0:
                    new_bases.append(
                        matrix(ring, 
                               Zkern * table_constructor[param][ii], 
                               sparse=True)
                               )
            table_constructor[param] = new_bases

        return table_constructor

    def _print_eigenvalues(self, table_constructor, sdp_result):
        r"""
        Helper to quickly print the eigenvalues of each X matrix from the result
        of an SDP.
        """
        block_index = 0
        for params in table_constructor.keys():
            ftype = params[1]
            for plus_index in range(len(table_constructor[params])):
                X_approx = matrix(sdp_result['X'][block_index + plus_index])
                X_eigenvalues = X_approx.eigenvalues()
                print("{} index {} has eigenvalues {}\n\n".format(
                    ftype, plus_index, X_eigenvalues
                    ))
            block_index += len(table_constructor[params])
    
    def _tables_to_sdp_data(self, table_constructor, prev_data=None):
        r"""
        Helper to transform the data from the multiplication 
        tables to an SDP input
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
            table = self.mul_project_table(
                ns, ns, ftype, ftype_inj=[], target_size=target_size
                )
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
                                mat_inds.extend(
                                    [gg+1, block_index + plus_index, 
                                     iinds[cc]+1, jinds[cc]+1]
                                    )
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
        flag_num, constraints_vals, constraints_flags_vec, one_vector = \
            constraints_data
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
    
    def _make_sdp_data_integer(self, sdp_data):
        from sage.arith.functions import lcm
        block_sizes, target, mat_inds, mat_vals = sdp_data

        mat_vals_factor = 1
        for xx in mat_vals:
            if xx!=0:
                mat_vals_factor = lcm(mat_vals_factor, QQ(xx).denominator())
        mat_vals = [Integer(xx*mat_vals_factor) for xx in mat_vals]
        
        target_factor = 1
        for xx in target:
            if xx!=0:
                target_factor = lcm(target_factor, QQ(xx).denominator())
        target = [Integer(xx*target_factor) for xx in target]

        return (block_sizes, target, mat_inds, mat_vals)

    def _get_relevant_ftypes(self, target_size):
        r"""
        Returns the ftypes useful for optimizing up to `target_size`
        """
        plausible_sizes = list(range(1, target_size))
        ftype_pairs = []
        for fs, ns in itertools.combinations(plausible_sizes, r=int(2)):
            if ns+ns-fs <= target_size:
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
            ftypes = [flag.subflag([], ftype_points=list(range(fs))) \
                    for flag in ftype_flags]
            for xx in ftypes:
                ftype_data.append((ns, xx, target_size))
        ftype_data.sort()
        return ftype_data
    
    def _create_table_constructor(self, ftype_data, target_size):
        r"""
        Table constructor is a dictionary that holds the data to construct
        all the multiplication tables. 
        
        For each ftype and base change it provides the data to create 
        the multiplication table. Also pre-computes the multiplication 
        tables if they are not calculated yet.
        """
        
        sym_asym_mats = [
            self.sym_asym_bases(dat[0], dat[1]) for dat in ftype_data
            ]

        table_constructor = {}
        for ii, dat in (pbar := tqdm(enumerate(ftype_data))):
            ns, ftype, target_size = dat
            #pre-calculate the table here
            table = self.mul_project_table(
                ns, ns, ftype, ftype_inj=[], target_size=target_size
                )
            if table==None:
                pbar.set_description(
                    "{} ({}) had singular table!".format(ftype, ns)
                    )
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
                df = target_size - nf + kf
                mult_table = self.mul_project_table(
                    nf, df, fv.ftype(), ftype_inj=[], target_size=target_size
                    )
                fvvals = fv.values()
                m = matrix(QQ, [vector(fvvals*mat) for mat in mult_table])
                positives_list_exact += list(m.T)
                pbar.set_description(
                    "Done with positivity constraint {}".format(ii)
                    )
            constraints_vals = [0]*len(positives_list_exact)
        
        # The one vector is also calculated here and is a linear constraint
        if target_element.ftype().size()==0:
            one_vector = vector([1]*len(base_flags))
        else:
            one_vector = (target_element.ftype().project()<<(
                target_size - target_element.ftype().size()
                )).values()
        positives_list_exact.extend([one_vector, one_vector*(-1)])
        constraints_vals.extend([1, -1])
        
        return len(base_flags), constraints_vals, \
            positives_list_exact, one_vector
    
    def _round_sdp_solution_no_phi(self, sdp_result, sdp_data, 
                                   table_constructor, constraints_data, 
                                   denom=1024):
        import numpy as np
        from numpy import linalg as LA
        from sage.functions.other import ceil

        #unpack variables

        block_sizes, target_list_exact, mat_inds, mat_vals = sdp_data
        target_vector_exact = vector(target_list_exact)
        flags_num, constraints_vals, positives_list_exact, one_vector = \
            constraints_data
        positives_matrix_exact = matrix(
            QQ, len(positives_list_exact), flags_num, positives_list_exact
            )
        
        # find the one_vector from the equality constraint
        one_vector_exact = positives_matrix_exact.rows()[-2] 
        # remove the equality constraints
        positives_matrix_exact = positives_matrix_exact[:-2, :] 

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
                Xr = matrix(QQ, Xr) + \
                diagonal_matrix(QQ, [eminr]*len(X), sparse=True)
            X_matrices_rounded.append(Xr)
        X_matrices_flat = [
            vector(_flatten_matrix(X.rows(), doubled=False)) \
                for X in (X_matrices_rounded)
            ]

        e_vector_approx = sdp_result['X'][-1][:-2]
        e_vector_rounded = vector(QQ, 
            _round_list(e_vector_approx, force_pos=True, method=0, denom=denom)
            )
        
        phi_vector_approx = sdp_result['y']
        phi_vector_rounded = vector(QQ, 
            _round_list(phi_vector_approx, force_pos=True, method=0, denom=denom)
            )

        slacks = target_vector_exact - positives_matrix_exact.T*e_vector_rounded
        block_index = 0
        print("Calculating resulting bound")
        for params in tqdm(table_constructor.keys()):
            ns, ftype, target_size = params
            table = self.mul_project_table(
                ns, ns, ftype, ftype_inj=[], target_size=target_size
                )
            for gg, morig in enumerate(table):
                for plus_index, base in enumerate(table_constructor[params]):
                    block_dim = block_sizes[block_index + plus_index]
                    X_flat = X_matrices_flat[block_index + plus_index]
                    M = base * morig * base.T
                    M_flat_vector_exact = vector(QQ, 
                        _flatten_matrix(M.rows(), doubled=True)
                        )
                    slacks[gg] -= M_flat_vector_exact*X_flat
            block_index += len(table_constructor[params])
        # scale back slacks with the one vector, the minimum is the final result
        result = min(
            [slacks[ii]/oveii for ii, oveii in \
             enumerate(one_vector_exact) if oveii!=0]
            )
        # pad the slacks, so it is all positive where it counts
        slacks -= result*one_vector_exact
        
        print("The rounded result is {}".format(result))
        
        return result, X_matrices_rounded, e_vector_rounded, \
            slacks, [phi_vector_rounded]
    
    def _round_sdp_solution_phi(self, sdp_result, sdp_data, 
                                table_constructor, constraints_data, 
                                phi_vectors_exact, denom=1024, ring=QQ):
        r"""
        Round the SDP results output to get something exact.
        """
        
        #unpack variables
        block_sizes, target_list_exact, mat_inds, mat_vals = sdp_data
        target_vector_exact = vector(ring, target_list_exact)
        flags_num, constraints_vals, positives_list_exact, one_vector = \
            constraints_data
        positives_matrix_exact = matrix(
            ring, len(positives_list_exact), flags_num, positives_list_exact
            )
        
        no_constr = len(phi_vectors_exact)==0
        phi_vector_exact = vector(
            ring, 
            [0]*positives_matrix_exact.ncols()
            ) if no_constr else phi_vectors_exact[0]
        
        # find the one_vector from the equality constraint
        one_vector_exact = positives_matrix_exact.rows()[-2] 
        # remove the equality constraints
        positives_matrix_exact = positives_matrix_exact[:-2, :] 
        # same as |F_n|
        flags_num = -block_sizes[-2] 

        # dim: |F_n|, c vector, primal slack for flags
        c_vector_approx = vector(sdp_result['X'][-2]) 
        # as above but rounded, in case it is needed
        # c_vector_rounded = vector(
        #     QQ, 
        #     _round_list(c_vector_approx, method=0, denom=denom)
        #     ) 

        # The F (FF) flag indecies where the c vector is zero/nonzero
        c_zero_inds = [
            FF for FF, xx in enumerate(c_vector_approx) if \
                (abs(xx)<1e-6 or phi_vector_exact[FF]!=0)
                ]
        # c_nonzero_inds = [FF for FF in range(flags_num) \ 
        # if FF not in c_zero_inds]


        # same as m, number of positive constraints (-2 for the equality)
        positives_num = -block_sizes[-1] - 2 

        # dim: m, witness that phi is positive
        phi_pos_vector_exact = positives_matrix_exact*phi_vector_exact 

        # dim: m, the e vector, primal slack for positivitives
        e_vector_approx = vector(sdp_result['X'][-1][:-2])
        # as above but rounded
        e_vector_rounded = vector(QQ, 
                                  _round_list(e_vector_approx, method=0, denom=denom)
                                  ) 

        # The f (ff) positivity constraints where the e vector is zero/nonzero
        e_zero_inds = [
            ff for ff, xx in enumerate(e_vector_approx) if \
                (abs(xx)<1e-6 or phi_pos_vector_exact[ff]!=0)
                ]
        e_nonzero_inds = [
            ff for ff in range(positives_num) if ff not in e_zero_inds
            ]



        bound_exact = target_vector_exact*phi_vector_exact 
        # the constraints for the flags that are exact
        corrected_target_relevant_exact = vector(
            ring, 
            [target_vector_exact[FF] - bound_exact for FF in c_zero_inds]
            )
        # the d^f_F matrix, but only the relevant parts for the rounding
        # so F where c_F = 0 and f where e_f != 0
        positives_matrix_relevant_exact = matrix(
            ring, len(e_nonzero_inds), len(c_zero_inds), 
            [[positives_matrix_exact[ff][FF] for FF in c_zero_inds] \
             for ff in e_nonzero_inds]
             )
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
        
        M_flat_relevant_matrix_exact = matrix(
            ring, len(c_zero_inds), 0, 0, sparse=True
            )
        # The rounded X values flattened to a list
        X_flat_vector_rounded = [] 
        block_index = 0
        block_info = []
        for params in table_constructor.keys():
            ns, ftype, target_size = params
            table = self.mul_project_table(
                ns, ns, ftype, ftype_inj=[], target_size=target_size
                )

            for plus_index, base in enumerate(table_constructor[params]):
                block_info.append([ftype, base])
                X_approx = sdp_result['X'][block_index + plus_index]
                X_flat_vector_rounded += \
                    _round_list(
                        _flatten_matrix(X_approx), method=0, denom=denom
                        )

                M_extra = []

                for FF in c_zero_inds:
                    M_FF = table[FF]
                    M_extra.append(
                        _flatten_matrix(
                            (base * M_FF * base.T).rows(), doubled=True
                            )
                        )

                M_flat_relevant_matrix_exact = \
                    M_flat_relevant_matrix_exact.augment(
                        matrix(ring, M_extra)
                        )
            block_index += len(table_constructor[params])


        # 
        # Append the relevant M matrix and the X with the additional values from
        # the positivity constraints. 
        #
        # Then correct the x vector values
        # 

        M_matrix_final = M_flat_relevant_matrix_exact.augment(
            positives_matrix_relevant_exact.T
            )
        x_vector_final = vector(
            ring, 
            X_flat_vector_rounded+e_nonzero_list_rounded
            )


        # Correct the values of the x vector, based on the minimal L_2 norm
        x_vector_corr = x_vector_final - M_matrix_final.T * \
        (M_matrix_final * M_matrix_final.T).pseudoinverse() * \
        (M_matrix_final*x_vector_final - corrected_target_relevant_exact) 
        
        #
        # Recover the X matrices and e vector from the corrected x
        #

        e_nonzero_vector_corr = x_vector_corr[-len(e_nonzero_inds):]
        if len(e_nonzero_vector_corr)>0 and min(e_nonzero_vector_corr)<0:
            print("Linear coefficient is negative: {}".format(
                min(e_nonzero_vector_corr)
                ))
            e_nonzero_vector_corr = [max(xx, 0) for xx in e_nonzero_vector_corr]
        e_vector_dict = dict(zip(e_nonzero_inds, e_nonzero_vector_corr))
        e_vector_corr = vector(ring, positives_num, e_vector_dict)
        
        
        X_final = []
        slacks = target_vector_exact - positives_matrix_exact.T*e_vector_corr
        block_index = 0
        print("Calculating resulting bound")
        for params in tqdm(table_constructor.keys()):
            ns, ftype, target_size = params
            table = self.mul_project_table(
                ns, ns, ftype, ftype_inj=[], target_size=target_size
                )
            for plus_index, base in enumerate(table_constructor[params]):
                block_dim = block_sizes[block_index + plus_index]
                X_ii_small, x_vector_corr = _unflatten_matrix(
                    x_vector_corr, block_dim
                    )
                X_ii_small = matrix(ring, X_ii_small)
                
                # verify semidefiniteness
                if not X_ii_small.is_positive_semidefinite():
                    print("Rounded X matrix "+ 
                          "{} is not semidefinite: {}".format(
                              block_index+plus_index, 
                              min(X_ii_small.eigenvalues())
                              ))
                    return None
                
                # update slacks
                for gg, morig in enumerate(table):
                    M = base * morig * base.T
                    M_flat_vector_exact = vector(
                        _flatten_matrix(M.rows(), doubled=True)
                        )
                    slacks[gg] -= M_flat_vector_exact*vector(
                        _flatten_matrix(X_ii_small.rows(), doubled=False)
                        )
                
                X_final.append(X_ii_small)
            block_index += len(table_constructor[params])
        
        # scale back slacks with the one vector, the minimum is the final result
        result = min([slacks[ii]/oveii \
                      for ii, oveii in enumerate(one_vector_exact) if \
                        oveii!=0])
        # pad the slacks, so it is all positive where it counts
        slacks -= result*one_vector_exact
        
        return result, X_final, e_vector_corr, slacks, phi_vectors_exact
    
    def _fix_X_bases(self, table_constructor, X_original):
        r"""
        Transforms the X matrices to a base that agrees with the original 
        list of flags
        
        Basically undoes the sym/asym changes and the reductions by the 
        constructions' kernels.
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
        Transforms the X matrices to a base that agrees with the original 
        list of flags
        
        Basically undoes the sym/asym changes and the reductions by the 
        constructions' kernels.
        Also changes the X matrices to the P L D L.T P.T form
        """
        if X_original==None:
            return None
        P_dicts = []
        L_mats = []
        D_vecs = []
        block_index = 0
        for params in table_constructor.keys():
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
    
    def _format_optimizer_output(self, table_constructor, mult=1, 
                                 sdp_output=None, rounding_output=None, 
                                 file=None):
        r"""
        Formats the outputs to a nice certificate
        
        The result contains: the final bound, the X matrices, the linear 
        coefficients, the slacks, a guess or exact construction, the 
        list of base flags, the list of used (typed) flags
        """
        
        target_size = 0
        typed_flags = {}
        for params in table_constructor.keys():
            ns, ftype, target_size = params
            typed_flags[(ns, ftype)] = self.generate_flags(ns, ftype)
        
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
        
        if file!=None and file!="":
            if not file.endswith(".pickle"):
                file += ".pickle"
            with open(file, "wb") as file_handle:
                pickle.dump(cert_dict, file_handle)
        if file=="notebook":
            return cert_dict
        return result
    
    def optimize_problem(self, target_element, target_size, maximize=True, 
                         positives=None, construction=None, file=None, 
                         exact=False, denom=1024, ring=QQ):
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
        """
        from csdpy import solve_sdp
        import time

        #
        # Initial setup
        #
        base_flags = self.generate_flags(target_size)
        print("Base flags generated, their number is {}".format(
            len(base_flags)
            ))
        mult = -1 if maximize else 1
        target_vector_exact = (
            target_element.project()*(mult) << \
                (target_size - target_element.size())
                ).values()
        sdp_data = self._target_to_sdp_data(target_vector_exact)
        
        #
        # Create the relevant ftypes
        #
        
        ftype_data = self._get_relevant_ftypes(target_size)
        print("The relevant ftypes are constructed, their " + 
              "number is {}".format(len(ftype_data)))
        flags = [self.generate_flags(dat[0], dat[1]) for dat in ftype_data]
        flag_sizes = [len(xx) for xx in flags]
        print("Block sizes before symmetric/asymmetric change is" + 
              " applied: {}".format(flag_sizes))
        
        #
        # Create the table constructor and add it to sdp_data
        #
        
        table_constructor = self._create_table_constructor(
            ftype_data, target_size
            )
        sdp_data = self._tables_to_sdp_data(
            table_constructor, prev_data=sdp_data
            )
        print("Tables finished")

        #
        # Add constraints data and add it to sdp_data
        #
        
        constraints_data = self._create_constraints_data(
            positives, target_element, target_size
            )
        sdp_data = self._constraints_to_sdp_data(
            constraints_data, prev_data=sdp_data
            )
        print("Constraints finished")
        
        #
        # If construction is None or [] then run the optimizer 
        # without any construction
        #
        
        if construction==None or construction==[]:
            print("Running sdp without construction. " + 
                  "Used block sizes are {}".format(sdp_data[0]))
            
            time.sleep(float(0.1))
            initial_sol = solve_sdp(*sdp_data)
            time.sleep(float(0.1))

            # Format the result and return it if floating point values are fine
            if (not exact):
                if file==None:
                    return initial_sol['primal'] * mult
                return self._format_optimizer_output(
                    table_constructor, 
                    mult=mult, 
                    sdp_output=initial_sol, 
                    file=file
                    )
            
            # Guess the construction in this case
            if construction==None:
                one_vector = constraints_data[-1]
                #phi_vector_original = initial_sol['y']
                phi_vector_rounded, error_coeff = _round_adaptive(
                    initial_sol['y'], one_vector
                    )
                if error_coeff<1e-6:
                    alg = FlagAlgebra(self, QQ)
                    construction = alg(target_size, phi_vector_rounded)
                    phipr = str(construction)
                    print("The initial run gave an accurate "+
                          "looking construction")
                    if len(phipr)<1000:
                        print("Rounded construction vector "+
                              "is: \n{}".format(phipr))
                else:
                    print("The initial run didn't provide an "+
                          "accurate construction")
                    construction = []
            
            # If nothing was provided or the guess failed, 
            # then round the current solution
            if construction==[]:
                rounding_output = self._round_sdp_solution_no_phi(
                    initial_sol, sdp_data, table_constructor, 
                    constraints_data, denom=denom)
                if file==None:
                    return rounding_output[0] * mult
                return self._format_optimizer_output(
                    table_constructor, 
                    mult=mult, 
                    rounding_output=rounding_output, 
                    file=file
                    )
        
        
        #
        # Run the optimizer (again if a construction was guessed) 
        # with the construction
        #
        
        if isinstance(construction, FlagAlgebraElement):
            phi_vectors_exact = [construction.values()]
        else:
            phi_vectors_exact = [xx.values() for xx in construction]

        #
        # Adjust the table to consider the kernel from y_rounded
        #

        print("Adjusting table with kernels from construction")
        table_constructor = self._adjust_table_phi(
            table_constructor, phi_vectors_exact, ring=ring
            )
        sdp_data = self._target_to_sdp_data(target_vector_exact)
        sdp_data = self._tables_to_sdp_data(
            table_constructor, prev_data=sdp_data
            )
        sdp_data = self._constraints_to_sdp_data(
            constraints_data, prev_data=sdp_data
            )
        
        #
        # Then run the optimizer
        #
        
        print("Running SDP after kernel correction. "+
              "Used block sizes are {}".format(sdp_data[0]))
        time.sleep(float(0.1))
        final_sol = solve_sdp(*sdp_data)
        time.sleep(float(0.1))

        # Quickly deal with the case when no rounding is needed
        if (not exact):
            if file==None:
                return final_sol['primal'] * mult
            return self._format_optimizer_output(
                table_constructor, 
                mult=mult, 
                sdp_output=final_sol, 
                file=file
                )
        
        
        print("Starting the rounding of the result")
        rounding_output = self._round_sdp_solution_phi(
            final_sol, sdp_data, table_constructor, 
            constraints_data, phi_vectors_exact, 
            denom=denom, ring=ring
            )
        if rounding_output==None:
            print("Rounding based on construction was unsuccessful")
            rounding_output = self._round_sdp_solution_no_phi(
                final_sol, sdp_data, table_constructor, 
                constraints_data, denom=denom
                )
        
        print("Final rounded bound is {}".format(rounding_output[0]*mult))
        
        if file==None:
            return rounding_output[0] * mult
        return self._format_optimizer_output(
            table_constructor, 
            mult=mult, 
            rounding_output=rounding_output, 
            file=file
            )
    
    optimize = optimize_problem
    
    def external_optimize(self, target_element, target_size, maximize=True, 
                          positives=None, construction=None, file=None, 
                          specific_ftype=None):
        if (not isinstance(file, str)) or file=="":
            raise ValueError("File name is invalid.")
        if not file.endswith(".dat-s"):
            if file.endswith(".dat"):
                file += "-s"
            else:
                file += ".dat-s"
        #
        # Initial setup
        #
        base_flags = self.generate_flags(target_size)
        print("Base flags generated, their number "+
              "is {}".format(len(base_flags)))
        mult = -1 if maximize else 1
        target_vector_exact = (
            target_element.project()*(mult)<<
            (target_size - target_element.size())
            ).values()
        sdp_data = self._target_to_sdp_data(target_vector_exact)
        
        #
        # Create the relevant ftypes
        #
        
        if specific_ftype==None:
            ftype_data = self._get_relevant_ftypes(target_size)
        else:
            ftype_data = specific_ftype
        print("The relevant ftypes are constructed, "+
              "their number is {}".format(len(ftype_data)))
        flags = [self.generate_flags(dat[0], dat[1]) for dat in ftype_data]
        flag_sizes = [len(xx) for xx in flags]
        print("Block sizes before symmetric/asymmetric change " + 
              "is applied: {}".format(flag_sizes))
        
        #
        # Create the table constructor and add it to sdp_data
        #
        
        table_constructor = self._create_table_constructor(
            ftype_data, target_size
            )
        if not (construction==None or construction==[]):
            if isinstance(construction, FlagAlgebraElement):
                phi_vectors_exact = [construction.values()]
            else:
                phi_vectors_exact = [xx.values() for xx in construction]
            print("Adjusting table with kernels from construction")
            table_constructor = self._adjust_table_phi(
                table_constructor, phi_vectors_exact
                )
        sdp_data = self._tables_to_sdp_data(
            table_constructor, prev_data=sdp_data
            )
        print("Tables finished")

        #
        # Create constraints data and add it to sdp_data
        #
        
        constraints_data = self._create_constraints_data(
            positives, target_element, target_size
            )
        sdp_data = self._constraints_to_sdp_data(
            constraints_data, prev_data=sdp_data
            )
        print("Constraints finished")


        #
        # Make sdp data integer and write it to a file
        #
        sdp_data = self._make_sdp_data_integer(sdp_data)
        
        with open(file, "a") as file:
            block_sizes, target, mat_inds, mat_vals = sdp_data

            file.write("{}\n{}\n".format(len(target), len(block_sizes)))
            file.write(" ".join(map(str, block_sizes)) + "\n")
            for xx in target:
                file.write("%.1e " % xx)
            file.write("\n")

            for ii in range(len(mat_vals)):
                file.write("{} {} {} {} {}\n".format(
                    mat_inds[ii*4 + 0],
                    mat_inds[ii*4 + 1],
                    mat_inds[ii*4 + 2],
                    mat_inds[ii*4 + 3],
                    "%.1e" % mat_vals[ii]
                ))

    def verify_certificate(self, file_or_cert, target_element, target_size, 
                           maximize=True, positives=None, construction=None):
        r"""
        Verifies the certificate provided by the optimizer 
        written to `file`
        """
        
        #
        # Load the certificate file (only pickle is supported)
        #
        
        if isinstance(file_or_cert, str):
            file = file_or_cert
            if not file.endswith(".pickle"):
                file += ".pickle"
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
            print("Linear constraint's coefficient is negative ", 
                  min(e_values))
            return -1
        
        X_flats = certificate["X matrices"]
        # Ps = certificate["P dictionaries"]
        # Ds = certificate["D vectors"]
        # Ls = certificate["L matrices"]
        print("Checking X matrices")
        for ii, Xf in tqdm(enumerate(X_flats)):
            X = matrix(QQ, _unflatten_matrix(Xf)[0])
            if not (X.is_positive_semidefinite()):
                print("Solution is not valid!")
                print("Matrix {} is not semidefinite".format(ii))
                return -1
            # P = matrix(QQ, len(Ddiag), len(Ddiag), Ps[ii], sparse=True)
            # D = diagonal_matrix(QQ, Ddiag, sparse=True)
            # Larr, _ = _unflatten_matrix(
            #     Ls[ii], dim=len(Ddiag), 
            #     doubled=False, upper=True
            #     )
            # L = matrix(QQ, Larr).T
            # PL = P*L
            # X = PL * D * PL.T
            # X_flats.append(vector(QQ, _flatten_matrix(X.rows())))

        print("Solution matrices are all positive semidefinite, " + 
              "linear coefficients are all non-negative")

        #
        # Initial setup
        #

        mult = -1 if maximize else 1
        base_flags = certificate["base flags"]
        target_vector_exact = (
            target_element.project()*(mult)<<\
                (target_size - target_element.size())
                ).values()
        if target_element.ftype().size()==0:
            one_vector = vector([1]*len(base_flags))
        else:
            one_vector = (
                target_element.ftype().project()<<\
                    (target_size - target_element.ftype().size())
                    ).values()
        
        ftype_data = list(certificate["typed flags"].keys())

        #
        # Create the semidefinite matrix data
        #

        table_list = []
        print("Calculating multiplication tables")
        for ii, dat in tqdm(enumerate(ftype_data)):
            ns, ftype = dat
            #calculate the table
            table = self.mul_project_table(
                ns, ns, ftype, ftype_inj=[], target_size=target_size
                )
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
                mult_table = self.mul_project_table(
                    nf, df, fv.ftype(), ftype_inj=[], target_size=target_size
                    )
                fvvals = fv.values()
                m = matrix(QQ, [vector(fvvals*mat) for mat in mult_table])
                positives_list_exact += list(m.T)
                print("Done with positivity constraint {}".format(ii))
        else:
            e_values = vector(QQ, 0)
        positives_matrix_exact = matrix(
            QQ, 
            len(positives_list_exact), len(base_flags), 
            positives_list_exact
            )

        print("Done calculating linear constraints")

        #
        # Calculate the bound the solution provides
        #
        
        print("Calculating the bound provided by the certificate")
        
        slacks = target_vector_exact - positives_matrix_exact.T*e_values
        for ii, table in tqdm(enumerate(table_list)):
            for gg, mat_gg in enumerate(table):
                mat_flat = vector(
                    _flatten_matrix(mat_gg.rows(), doubled=True)
                    )
                slacks[gg] -= mat_flat * X_flats[ii]
        result = min(
            [slacks[ii]/oveii for ii, oveii in enumerate(one_vector) \
             if oveii!=0]
            )
        result *= mult
        print("The solution is valid, it proves "+
              "the bound {}".format(result))

        return result
    
    verify = verify_certificate
    
    
    #Generating flags
    def _guess_number(self, n):
        if n==0:
            return 1
        excluded = self.get_total_excluded(n)
        key = ("generate", n, 
               self._serialize(excluded), self.empty()._serialize())
        loaded = self._load(key=key)
        if loaded != None:
            return len(loaded)
        if self._sources==None:
            max_arity = -1
            for xx in self._signature:
                max_arity = max(max_arity, self._signature[xx]["arity"])
            if max_arity==1 or n<max_arity:
                return 1
            
            check_bits = 0
            for xx in self._signature:
                arity =self._signature[xx]["arity"]
                if arity != 1:
                    factor = 1
                    if self._signature[xx]["ordered"]:
                        factor = factorial(arity)
                    check_bits += factor*(binomial(n-2, arity-2))
            
            prev_guess = len(self.generate(n-1))
            sign_perm = len(self._signature_perms())
            return binomial(prev_guess + 1, 2) * (2**check_bits) * sign_perm
        else:
            t0, t1 = self._sources
            guess0 = t0._guess_number(n)
            guess1 = t1._guess_number(n)
            return guess0*guess1*8
    
    def generate_with_overlaps(self, theory0, theory1, n):
        ls0 = theory0.generate(n)
        ls1 = theory1.generate(n)
        return overlap_generator(n, self, ls0, ls1, tuple())

    def no_question(self, val=None):
        if val==None:
            val = True
        self._no_question = val

    def get_total_excluded(self, n):
        if self._sources == None:
            ret = [xx for xx in self._excluded if xx.size()<=n]
        else:
            ret = [xx for xx in self._excluded if xx.size()<=n]
            ret += list(self._sources[0].get_total_excluded(n))
            ret += list(self._sources[1].get_total_excluded(n))
        return tuple(ret)

    def generate_flags(self, n, ftype=None, run_bound=500000):
        r"""
        Returns the list of flags with a given size and ftype

        INPUT:

        - ``n`` -- integer; the size of the returned structures
        - ``ftype`` -- Flag; the ftype of the returned structures

        OUTPUT: List of all flags with given size and ftype
        
        EXAMPLES::

        There are 4 graphs on 3 vertices. Flags with empty
        ftype correspond to elements of the theory ::

            
            sage: len(GraphTheory.generate_flags(3))
            4
        
        .. NOTE::

            :func:`generate` is an alias for this.
            See the notes on :func:`optimize_problem`. A large `n` can
            result in large number of structures.
        """
        
        #Handling edge cases
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

        
        #Trying to load
        excluded = self.get_total_excluded(n)
        key = ("generate", n, 
               self._serialize(excluded), ftype._serialize())
        loaded = self._load(key=key)
        if loaded != None:
            return loaded

        if ftype.size()==0:
            def just_generate():
                if self._sources != None:
                    t0, t1 = self._sources
                    small0 = t0.generate(n)
                    small1 = t1.generate(n)
                    small_excl = tuple(
                        [xx for xx in self._excluded if xx.size()<=n]
                        )
                    ret = overlap_generator(n, self, 
                                            small0, small1, 
                                            small_excl)
                else:
                    prev = self.generate_flags(n-1, run_bound=run_bound)
                    ret = inductive_generator(n, self, prev, excluded)
                return ret

            # No ftype generation needed, just generate inductively
            if run_bound==infinity or n<=3 or self._no_question:
                ret = just_generate()
            else:
                guess = self._guess_number(n)
                if guess < run_bound:
                    ret = just_generate()
                else:
                    confirm = input("This might take a while: " + 
                                    "{}. Continue? y/n\n".format(guess))
                    if "y" in confirm.lower():
                        return self.generate_flags(n, run_bound=infinity)
                    else:
                        raise RuntimeError("Calculation interrupted")
        else:
            # First generate the structures without ftypes then find them
            empstrs = self.generate_flags(n)
            guess = len(empstrs) * falling_factorial(n, ftype.size())
            if guess < run_bound:
                ret = self._find_ftypes(empstrs, ftype)
            else:
                confirm = input("This might take a while: " + 
                                "{}. Continue? y/n\n".format(guess))
                if "y" in confirm.lower():
                    return self.generate_flags(n, ftype, run_bound=infinity)
                else:
                    raise RuntimeError("Calculation interrupted")
        self._save(ret, key)
        return ret
    
    def _find_ftypes(self, empstrs, ftype):
        import multiprocessing as mp
        pool = mp.Pool(mp.cpu_count()-1)
        pares = pool.map(ftype.ftypes_inside, empstrs)
        pool.close(); pool.join()
        return tuple(itertools.chain.from_iterable(pares))
    
    generate = generate_flags

    def exclude(self, structs=None, force=False):
        #Set up structs to contain all we want to exclude
        if structs==None:
            structs = []
        elif type(structs)==Flag or type(structs)==Pattern:
            structs = [structs]
        if not force:
            structs += list(self._excluded)
        
        #Make structs sorted, so it is as small as possible
        structs.sort(key=lambda x : x.size())
        self._excluded = tuple()

        for xx in structs:
            if isinstance(xx, Pattern):
                if xx.theory()!=self:
                    xx = self.Pattern(xx.size(), **xx.blocks())
                extension = xx.compatible_flags()
                self._excluded = tuple(list(self._excluded) + extension)
            elif isinstance(xx, Flag):
                if xx.theory()!=self:
                    xxpat = xx.as_pattern()
                    selfpat = self.Pattern(xxpat.size(), **xxpat.blocks())
                    extension = selfpat.compatible_flags()
                    self._excluded = tuple(list(self._excluded) + extension)
                else:
                    if xx in self.generate(xx.size()):
                        self._excluded = tuple(list(self._excluded) + [xx])

    def reset_exclude(self):
        self.exclude(force=True)

    reset = reset_exclude

    def match_pattern(self, pattern):
        if pattern is Flag:
            return [pattern]
        ss = pattern
        if len(ss.ftype_points())!=0:
            ss = ss.subpattern()
        return [
            xx for xx in self.generate_flags(ss.size(), ss.ftype()) \
                if ss.is_compatible(xx)
                ]
    
    match = match_pattern

    @lru_cache(maxsize=None)
    def _signature_perms(self):
        
        terms = self._signature.keys()
        groups = [self._signature[xx]["group"] for xx in terms]

        grouped_terms = {}
        for group, term in zip(groups, terms):
            grouped_terms.setdefault(group, []).append(term)

        group_permutations = {
            group: list(itertools.permutations(terms)) \
                for group, terms in grouped_terms.items()
            }

        grouped_permutations = itertools.product(
            *(group_permutations[group] for group \
              in sorted(grouped_terms.keys()))
            )

        all_permutations = []
        for grouped_perm in grouped_permutations:
            flat_perm = []
            for perm in grouped_perm:
                flat_perm.extend(perm)
            all_permutations.append(tuple(flat_perm))
        return all_permutations
    
    #Generating tables

    def mul_project_table(self, n1, n2, large_ftype, ftype_inj=None, 
                          target_size=None):
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

            sage: table = GraphTheory.mul_project_table(2, 2, GraphTheory(1, ftype_points=[0]), [])
            sage: table[1][0, 0]
            1/3
        """

        #Sanity checks
        large_size = large_ftype.size()
        if ftype_inj==None:
            ftype_inj = tuple(range(large_size))
        else:
            ftype_inj = tuple(ftype_inj)
            checklist = [ii for ii in ftype_inj if \
                         (ii not in range(large_size))]
            if len(checklist)!=0:
                raise ValueError("ftype_inj must map into the " + 
                                 "points of {}".format(large_ftype))
            if len(set(ftype_inj)) != len(ftype_inj):
                raise ValueError("ftype_inj must be injective " + 
                                 "(no repeated elements)")
        if target_size==None:
            target_size = n1+n2 - large_size

        #Trying to load
        excluded = self.get_total_excluded(target_size)
        key = ("table", (n1, n2, target_size), 
               large_ftype._serialize(), tuple(ftype_inj), 
               self._serialize(excluded))
        loaded = self._load(key=key)
        if loaded != None:
            return loaded
        
        N = target_size

        from sage.matrix.args import MatrixArgs
        import multiprocessing as mp
        ftype_inj = list(ftype_inj)
        large_size = large_ftype.size()
        small_ftype = large_ftype.subflag([], ftype_points=ftype_inj)
        small_size = small_ftype.size()
        ftype_remap = ftype_inj + [
            ii for ii in range(large_size) if (ii not in ftype_inj)
            ]
        
        Nflgs = self.generate(N, small_ftype)
        n1flgs = self.generate(n1, large_ftype)
        n2flgs = self.generate(n2, large_ftype)
        
        slist = tuple((
            flg, n1, n1flgs, n2, n2flgs, 
            ftype_remap, large_ftype, small_ftype
            ) for flg in Nflgs
            )
        
        pool = mp.Pool(mp.cpu_count()-1)
        mats = pool.map(self._density_wrapper, slist)
        pool.close(); pool.join()

        norm = falling_factorial(N - small_size, large_size - small_size) 
        norm *= binomial(N - large_size, n1 - large_size)
        ret = tuple([
            MatrixArgs(QQ, 
                       mat[0], mat[1], 
                       entries=mat[2]
                       ).matrix()/norm for mat in mats])
        
        self._save(ret, key=key)
        return ret
    
    mpt = mul_project_table
    table = mul_project_table
    
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
            xxid = xx.unique(weak=True)
            if xxid not in uniques:
                uniques.append(xxid)
                sym_base.append(xx.afae())
            else:
                sym_ind = uniques.index(xxid)
                asym_base.append(sym_base[sym_ind] - xx.afae())
                sym_base[sym_ind] += xx
        m_sym = matrix(
            len(sym_base), len(flags), 
            [xx.values() for xx in sym_base], sparse=True
            )
        m_asym = matrix(
            len(asym_base), len(flags), 
            [xx.values() for xx in asym_base], sparse=True
            )
        return m_sym, m_asym
    
    def _density_wrapper(self, ar):
        r"""
        Helper function used in the parallelization of calculating densities
        """
        return ar[0].densities(*ar[1:])

Theory = CombinatorialTheory

#Pre-defined theories
GraphTheory = Theory("Graph")
DiGraphTheory = Theory("DiGraph", arity=2, is_ordered=True)
ThreeGraphTheory = Theory("ThreeGraph", arity=3)
DiThreeGraphTheory = Theory("DiThreeGraph", arity=3, is_ordered=True)
FourGraphTheory = Theory("FourGraph", arity=4)
Color0 = Theory("Color0", relation_name="C0", arity=1)
Color1 = Theory("Color1", relation_name="C1", arity=1)
Color2 = Theory("Color2", relation_name="C2", arity=1)
Color3 = Theory("Color3", relation_name="C3", arity=1)
Color4 = Theory("Color4", relation_name="C4", arity=1)
Color5 = Theory("Color5", relation_name="C5", arity=1)
Color6 = Theory("Color6", relation_name="C6", arity=1)
Color7 = Theory("Color7", relation_name="C7", arity=1)

#Pre-defined symmetries
def CyclicSymmetry(n):
    succ = list(range(1, n)) + [0]
    ret = []
    for ii in range(n):
        ret.append([ii, succ[ii]])
        ret.append([ii, n + ii])
        ret.append([ii, 2*n + ii])
        ret.append([n + ii, 2*n + succ[ii]])
        ret.append([ii, 2*n + succ[ii]])
    return ret
K4mSymmetry = [
    [0, 6], [1, 6], [4, 6],
    [0, 7], [1, 7], [5, 7],
    [0, 8], [2, 8], [3, 8],
    [0, 9], [2, 9], [5, 9],
    [0, 10], [3, 10], [4, 10],
    [1, 11], [2, 11], [3, 11],
    [1, 12], [2, 12], [4, 12],
    [1, 13], [3, 13], [5, 13],
    [2, 14], [4, 14], [5, 14],
    [3, 15], [4, 15], [5, 15]
]
FullSymmetry = True
NoSymmetry = False

#Primitive rounding methods
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
    Unflatten a symmetric matrix, optionally correct for the doubled 
    non-diagonal elements
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

def _round(value, method=1, quotient_bound=7, denom_bound=9, 
           denom=1024):
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

def _round_list(ls, force_pos=False, method=1, quotient_bound=7, 
                denom_bound=9, denom=1024):
    r"""
    Helper function, to round a list
    """
    if force_pos:
        return [max(
            _round(xx, method, quotient_bound, denom_bound, denom), 0
            ) for xx in ls]
    else:
        return [_round(xx, method, quotient_bound, 
                       denom_bound, denom) for xx in ls]


def _round_matrix(mat, method=1, quotient_bound=7, denom_bound=9, 
                  denom=1024):
    r"""
    Helper function, to round a matrix
    """
    try:
        return matrix(QQ, [_round_list(xx, False, 
                                       method, quotient_bound, 
                                       denom_bound, denom
                                       ) for xx in mat])
    except:
        #This happens when a semidef constraint turns out to be just linear
        return diagonal_matrix(QQ, _round_list(mat, True, 
                                               method, quotient_bound, 
                                               denom_bound, denom))

def _round_adaptive(ls, onevec, denom=1024):
    r"""
    Adaptive rounding based on continued fraction and preserving 
    an inner product with `onevec`
    
    If the continued fraction rounding fails fall back to a simple 
    denominator method
    """
    best_vec = None
    best_error = 1000
    best_lcm = 1000000000
    
    orig = vector(ls)
    for resol1 in range(5, 20):
        resol2 = round(resol1*1.5)
        rls = vector([_round(xx, quotient_bound=resol1, denom_bound=resol2) \
                      for xx in orig])
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
