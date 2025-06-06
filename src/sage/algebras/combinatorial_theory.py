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
from sage.all import QQ, NN, Integer, ZZ, infinity, RR, RDF, RealField
from sage.algebras.flag import BuiltFlag, ExoticFlag, Pattern, inductive_generator, overlap_generator
from sage.algebras.flag_algebras import FlagAlgebra, FlagAlgebraElement

from sage.categories.sets_cat import Sets
from sage.all import vector, matrix, diagonal_matrix

from sage.misc.prandom import randint
from sage.arith.misc import falling_factorial, binomial, factorial
from sage.misc.functional import round
from sage.functions.other import ceil
from functools import lru_cache

import hashlib
import pickle
import os
from tqdm import tqdm

# This should really be in the doctests
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

# Primitive rounding methods
def _flatten_matrix(mat, doubled=False):
    r"""
    Flatten a symmetric matrix, optionally double non-diagonal elements
    """
    res = []
    factor = 2 if doubled else 1
    try:
        for ii in range(len(mat)):
            res.append(mat[ii][ii])
            res += [factor*mat[ii][jj] for jj in range(ii+1, len(mat))]
    except:
        for ii in range(len(mat)):
            res.append(mat[ii])
            res += [0 for jj in range(ii+1, len(mat))]
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

def _remove_kernel(mat, factor=1024, threshold=1e-4):
    d = mat.nrows()
    M_scaled = matrix(ZZ, [[round(xx*factor) for xx in vv] for vv in mat])
    M_augmented = matrix.identity(ZZ, d).augment(M_scaled)
    LLL_reduced = M_augmented.LLL()
    LLL_coeffs = LLL_reduced[:, :d]
    norm_test = LLL_coeffs * mat
    #print("norms are: ", " ".join([str(int(log(rr.norm(1)/d, 10))) for rr in norm_test]))
    kernel_base = [LLL_coeffs[ii] for ii,rr in enumerate(norm_test) if rr.norm(1)/d < threshold]
    #print("resulting kernel base is:\n", kernel_base)
    if len(kernel_base)==0:
        return mat, matrix.identity(d, sparse=True)
    K = matrix(ZZ, kernel_base).stack(matrix.identity(d))
    image_space = matrix(K.gram_schmidt()[0][len(kernel_base):, :], sparse=True)
    kernel_removed_mat = image_space * mat * image_space.T
    norm_factor = image_space * image_space.T
    mat_recover = image_space.T * norm_factor.inverse()
    return kernel_removed_mat, mat_recover

def custom_psd_test(mat):
    dim = mat.nrows()
    for ii in range(dim):
        pmin = mat[:ii+1, :ii+1]
        if pmin.det()<0:
            return False
    return True

class _CombinatorialTheory(Parent, UniqueRepresentation):
    def __init__(self, name):
        self._name = name
        self._excluded = tuple()
        self._printlevel = 1
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
        """
        return self._symmetries

    def sizes(self):
        return NN

    # Persistend data management
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
        Saves data to persistent storage.

        The file name is determined based on the provided arguments. If ``name`` is not
        given, a hashed key is generated from ``(self, key)`` and appended to ``self._name``.
        The file is saved in the directory given by ``path`` if provided, or in the directory
        returned by ``self._calcs_dir()`` if ``path`` is None. If ``path`` is an empty string,
        the file is saved in the current working directory.

        INPUT:

            - ``data`` -- any serializable object; the data to be saved.
            - ``key`` -- optional; used to generate the file name if ``name`` is not provided.
            - ``path`` -- optional; the directory path where the file will be saved.
            - ``name`` -- optional; explicit file name. Overrides key-based file naming if provided.

        OUTPUT:
            None

        EXAMPLES:

            sage: obj = SomeClass()  # Assume SomeClass defines _name and _calcs_dir appropriately.
            sage: obj._save("example data", key="sample")
            sage: file_name = obj._name + "." + hashlib.sha256(pickle.dumps((obj, "sample"))).hexdigest()
            sage: file_path = os.path.join(obj._calcs_dir(), file_name)
            sage: os.path.exists(file_path)
            True

        TESTS:

            sage: obj._save("data", name="test_file.pkl")
            sage: os.path.exists("test_file.pkl")
            True

            sage: obj._save("data")  # Neither key nor name provided.
            Traceback (most recent call last):
            ...
            ValueError: Either the key or the name must be provided!
        """
        if name == None:
            if key == None:
                raise ValueError("Either the key or the name must be provided!")
            serialized_key = pickle.dumps((self, key))
            hashed_key = hashlib.sha256(serialized_key).hexdigest()
            file_name = self._name + "." + hashed_key
        else:
            file_name = name

        if path == None:
            file_path = os.path.join(self._calcs_dir(), file_name)
        elif path == "":
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
        Loads data from persistent storage.

        The file name is determined by the provided ``key`` or explicitly by ``name``.
        If ``name`` is not provided, a hash is computed from ``(self, key)`` and appended to
        ``self._name``. The file is then sought in the directory specified by ``path`` (if given)
        or in the directory returned by ``self._calcs_dir()`` (if ``path`` is None). If the file
        does not exist, the function returns ``None``. Additionally, if a key is provided and the
        loaded data's key does not match the provided key, a warning is issued and ``None`` is returned.

        INPUT:

            - ``key`` -- optional; used to generate the file name if ``name`` is not provided.
            - ``path`` -- optional; directory path where the file is expected to be.
            - ``name`` -- optional; explicit file name. Overrides key-based file naming if provided.

        OUTPUT:
            The data stored in the file if found and valid; otherwise, ``None``.

        EXAMPLES:

            sage: obj = SomeClass()  # Assume SomeClass defines _name and _calcs_dir appropriately.
            sage: # Attempt to load data using a key.
            sage: data = obj._load(key="sample_key")
            sage: data  # Should return the saved data or None if not found.

        TESTS:

            sage: # Test: loading a non-existent file returns None.
            sage: obj._load(name="nonexistent_file.pkl")
            None

            sage: # Test: key mismatch triggers a warning and returns None.
            sage: obj._load(key="incorrect_key", name="existing_file.pkl")
            Traceback (most recent call last):
            ...
            Warning: Hash collision or corrupted data!
        """
        if key != None:
            serialized_key = pickle.dumps((self, key))
            hashed_key = hashlib.sha256(serialized_key).hexdigest()
            file_name = self._name + "." + hashed_key
        if name != None:
            file_name = name

        if path == None:
            file_path = os.path.join(self._calcs_dir(), file_name)
        else:
            if not os.path.exists(path):
                os.makedirs(path)
            file_path = os.path.join(path, file_name)
        
        if not os.path.exists(file_path):
            return None

        with open(file_path, "rb") as file:
            save_object = pickle.load(file)
        
        if key != None and save_object != None and save_object['key'] != key:
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
    
    def printlevel(self, val):
        self._printlevel = val

    # Optimizing and rounding

    def blowup_construction(self, target_size, parts, **kwargs):
        #
        # Initial setup, parsing parameters, setting up args
        #
        from sage.all import get_coercion_model, PolynomialRing
        from sage.all import multinomial_coefficients
        
        coercion = get_coercion_model()
        base_field = QQ
        if target_size < 0:
            raise ValueError("Target size must be non-negative.")
        
        # Parsing the parts
        part_vars_names = []
        part_weights_raw = []
        num_parts = 0
        if isinstance(parts, Integer):
            num_parts = parts
            if parts <= 0:
                raise ValueError("Number of parts must be positive")
            part_weights_raw = [1 / parts] * parts
        elif isinstance(parts, (list, tuple)):
            num_parts = len(parts)
            for p_spec in parts:
                part_weights_raw.append(p_spec)
                if isinstance(p_spec, str):
                    part_vars_names.append(p_spec)
                else:
                    base_field = coercion.common_parent(
                        base_field, p_spec.parent()
                        )
        else:
            raise TypeError(
                "The provided parts must be an integer or a list/tuple"
                )

        # Parsing the relations
        _temp_prob_vars = set()
        current_signature_map_for_scan = self.signature()
        for rel_name, definition in kwargs.items():
            if rel_name not in current_signature_map_for_scan:
                continue
            if isinstance(definition, dict):
                for _, prob_spec_scan in definition.items():
                    if isinstance(prob_spec_scan, str):
                        _temp_prob_vars.add(prob_spec_scan)
                    else:
                        base_field = coercion.common_parent(
                            base_field, prob_spec_scan
                            )

        # Create base ring
        prob_vars_names_list = sorted(list(_temp_prob_vars))
        # For now, just merge the params, maybe throw error if there's overlap
        all_var_names = sorted(list(set(part_vars_names + prob_vars_names_list)))
        if all_var_names:
            R = PolynomialRing(base_field, names=all_var_names)
            var_map = {
                name: R.gen(i) for i, name in enumerate(all_var_names)
                }
        else:
            R = base_field
            var_map = {}
        part_weights_R = [
            var_map[p_raw] if isinstance(p_raw, str) else R(p_raw) 
            for p_raw in part_weights_raw
            ]
        
        # Separate relations to deterministic and probabilistic
        deterministic_blowup_relations_R = set()
        probabilistic_blowup_relations_R = {}
        current_signature_map = self.signature()
        for rel_name, definition in kwargs.items():
            if rel_name not in current_signature_map:
                continue
            sig_details = current_signature_map[rel_name]
            isordered = sig_details["ordered"]
            if isinstance(definition, (list, tuple)):
                for edge in definition:
                    if isordered:
                        key = (rel_name, tuple(edge))
                    else:
                        key = (rel_name, tuple(sorted(edge)))
                    deterministic_blowup_relations_R.add(key)
            elif isinstance(definition, dict):
                for edge, prob in definition.items():
                    if isordered:
                        key = (rel_name, tuple(edge))
                    else:
                        key = (rel_name, tuple(sorted(edge)))
                    try:
                        prob_R_val = R(prob)
                    except:
                        prob_R_val = var_map[prob]
                    if prob_R_val == 1:
                        deterministic_blowup_relations_R.add(key)
                    elif prob_R_val != 0:
                        probabilistic_blowup_relations_R[key] = prob_R_val
            else:
                raise TypeError("Relations must be lists or dictionaries")

        #
        # Main calculation
        #
        
        # Helper to get probabilistic part
        def calculate_contribution(verts_assignment):
            vertex_indices = list(range(target_size))
            base_relations_for_this_outcome = {}
            potential_probabilistic_specs = []

            # Organize probabilistic relations
            for rel_name_sig, sig_details in current_signature_map.items():
                arity_sig = sig_details["arity"]
                is_ordered_sig = sig_details["ordered"]
                if arity_sig > target_size:
                    continue
                if is_ordered_sig:
                    vert_iterator = itertools.permutations(vertex_indices, arity_sig)
                else:
                    vert_iterator = itertools.combinations(vertex_indices, arity_sig)
                
                for v_tuple in vert_iterator:
                    parts_v_tuple = tuple(verts_assignment[v_idx] for v_idx in v_tuple)
                    key_in_blowup_def = (rel_name_sig, parts_v_tuple)
                    if key_in_blowup_def in deterministic_blowup_relations_R:
                        if rel_name_sig not in base_relations_for_this_outcome:
                            base_relations_for_this_outcome[rel_name_sig] = []
                        base_relations_for_this_outcome[rel_name_sig].append(list(v_tuple))
                    elif key_in_blowup_def in probabilistic_blowup_relations_R:
                        prob_for_rel = probabilistic_blowup_relations_R[key_in_blowup_def]
                        potential_probabilistic_specs.append({
                            "name": rel_name_sig, "verts": v_tuple, "prob": prob_for_rel
                        })
            
            num_potential_probabilistic = len(potential_probabilistic_specs)
            if num_potential_probabilistic == 0:
                structure = self(target_size, **base_relations_for_this_outcome)
                return structure.afae()
            
            ret_prob = 0
            for i in range(1 << num_potential_probabilistic):
                current_outcome_prob_R = 1
                relations_for_this_specific_outcome = {
                    k: list(v) for k, v in base_relations_for_this_outcome.items()
                }
                for j in range(num_potential_probabilistic):
                    spec = potential_probabilistic_specs[j]
                    is_present_in_outcome = (i >> j) & 1
                    if is_present_in_outcome:
                        current_outcome_prob_R *= spec["prob"]
                        rel_name, vert_list = spec["name"], list(spec["verts"])
                        if rel_name not in relations_for_this_specific_outcome:
                            relations_for_this_specific_outcome[rel_name] = []
                        relations_for_this_specific_outcome[rel_name].append(vert_list)
                    else:
                        current_outcome_prob_R *= (1 - spec["prob"])
                if current_outcome_prob_R == 0:
                    continue
                structure = self(target_size, **relations_for_this_specific_outcome)
                ret_prob += structure.afae() * current_outcome_prob_R
            return ret_prob

        # Get blowup components
        res = 0
        for exps_counts, mult_factor in multinomial_coefficients(
            num_parts, target_size
            ).items():
            verts_assignment_pattern = [
                part_idx for part_idx, count_in_part in enumerate(exps_counts) 
                for _ in range(count_in_part)
            ]
            weight_product_part_R = 1
            for i_part, count_in_part_val in enumerate(exps_counts):
                if count_in_part_val > 0:
                    weight_product_part_R *= \
                    (part_weights_R[i_part] ** count_in_part_val)
            current_coeff_R = mult_factor * weight_product_part_R
            term_contribution = calculate_contribution(verts_assignment_pattern)
            res += term_contribution * current_coeff_R
        return res

    def get_Z_matrices(self, phi_vector, table_constructor, test=True):
        Zs = []
        for param in table_constructor.keys():
            ns, ftype, target_size = param
            table = self.mul_project_table(ns, ns, ftype, ftype_inj=[], 
                                           target_size=target_size)
            Zm = [None for _ in range(len(table_constructor[param]))]
            for gg, morig in enumerate(table):
                if phi_vector[gg]==0:
                    continue
                for ii, base in enumerate(table_constructor[param]):
                    mat = base * morig * base.T
                    if Zm[ii]==None:
                        Zm[ii] = mat*phi_vector[gg]
                    else:
                        Zm[ii] += mat*phi_vector[gg]
            Zs.append(Zm)
            for Zii in Zm:
                if test and (not Zii.is_positive_semidefinite()):
                    self.fprint("Construction based Z matrix for " + 
                                "{} is not semidef: {}".format(
                                    ftype, min(Zii.eigenvalues())
                                    ))
        return Zs

    def _adjust_table_phi(self, table_constructor, phi_vectors_exact, 
                          test=False):
        r"""
        Helper to modify a table constructor, incorporating extra data from
        constructions (phi_vectors_exact)
        """
        if len(phi_vectors_exact)==0:
            return table_constructor

        to_pop = []
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
                        self.fprint("Construction based Z matrix for " + 
                              "{} is not semidef: {}".format(
                                  ftype, min(Zjj.eigenvalues())
                                  ))
                    if Z==None:
                        Z = Zjj
                    else:
                        Z.augment(Zjj)
                Zk = Z.kernel()
                Zkern = Zk.basis_matrix()
                #print("removing kernel:\n", Z.image().basis_matrix())
                if Zkern.nrows()>0:
                    new_bases.append(
                        matrix(Zkern * table_constructor[param][ii], 
                               sparse=True)
                               )
            if len(new_bases)!=0:
                table_constructor[param] = new_bases
            else:
                to_pop.append(param)
        
        for param in to_pop:
            table_constructor.pop(param, None)
        return table_constructor
    
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
        flag_num, constraints_vals, constraints_flags_vec, _, __ = \
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
        if self._printlevel>0:
            iterator = tqdm(enumerate(ftype_data))
            pbar = iterator
        else:
            iterator = enumerate(ftype_data)
            pbar = None
        for ii, dat in iterator:
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
            if not (pbar is None):
                pbar.set_description("Done with mult table for {}".format(ftype))
        return table_constructor
    
    def _create_constraints_data(self, positives, target_element, target_size):
        r"""
        Creates the data that holds the linear constraints
        """
        
        base_flags = self.generate_flags(target_size)
        factor_flags = []
        
        if positives == None:
            positives_list_exact = []
            constraints_vals = []
        else:
            positives_list_exact = []
            if self._printlevel>0:
                iterator = tqdm(range(len(positives)))
                pbar = iterator
            else:
                iterator = range(len(positives))
                pbar = None
            for ii in iterator:
                fv = positives[ii]
                if isinstance(fv, ExoticFlag) or isinstance(fv, BuiltFlag):
                    continue
                kf = fv.ftype().size()
                nf = fv.size()
                df = target_size - nf + kf
                mult_table = self.mul_project_table(
                    nf, df, fv.ftype(), ftype_inj=[], target_size=target_size
                    )
                factor_flags.append(self.generate(df, fv.ftype()))
                fvvals = fv.values()
                m = matrix([vector(fvvals*mat) for mat in mult_table])
                positives_list_exact += list(m.T)
                if not (pbar is None):
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
            positives_list_exact, one_vector, factor_flags
    
    def _round_sdp_solution_no_phi(self, sdp_result, sdp_data, 
                                   table_constructor, constraints_data, 
                                   **params):
        r"""
        Rounds the SDP solution without adjusting the phi vector.

        This function processes an SDP solution by rounding its matrices and slack variables.
        It performs the following steps:
        
        - Rounds each matrix in the SDP result using a specified denominator (default: 1024),
            correcting for any negative eigenvalues.
        - Flattens the rounded matrices and computes slack variables based on the provided
            table constructor and constraints data.
        - Determines the final result by scaling the slack variables with a one-vector extracted
            from the constraints.
        
        INPUT:

            - ``sdp_result`` -- dictionary;
            - ``sdp_data`` -- tuple;
            - ``table_constructor`` -- dictionary;
            - ``constraints_data`` -- tuple;
            - ``**kwargs`` -- keyword arguments for rounding parameters:
                * ``denom`` (integer, default: 1024): the denominator used for rounding.

        EXAMPLES:

            sage: 

        TESTS:

            sage: 
        """
        
        import numpy as np
        from numpy import linalg as LA
        from sage.functions.other import ceil
        from sage.all import coercion_model

        # set up parameters
        denom = params.get("denom", 1024)

        #unpack variables

        block_sizes, target_list_exact, mat_inds, mat_vals = sdp_data
        target_vector_exact = vector(target_list_exact)
        flags_num, ___, positives_list_exact, _, __ = \
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
        self.fprint("Rounding X matrices")
        if self._printlevel>0:
            iterator = tqdm(X_matrices_approx)
        else:
            iterator = X_matrices_approx
        for X in iterator:
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
        self.fprint("Calculating resulting bound")
        if self._printlevel > 0:
            iterator = tqdm(table_constructor.keys())
        else:
            iterator = table_constructor.keys()
        for params in iterator:
            ns, ftype, target_size = params
            table = self.mul_project_table(
                ns, ns, ftype, ftype_inj=[], target_size=target_size
                )
            for gg, morig in enumerate(table):
                for plus_index, base in enumerate(table_constructor[params]):
                    block_dim = block_sizes[block_index + plus_index]
                    X_flat = X_matrices_flat[block_index + plus_index]
                    M = base * morig * base.T
                    M_flat_vector_exact = vector( 
                        _flatten_matrix(M.rows(), doubled=True)
                        )
                    prod = M_flat_vector_exact*X_flat
                    slacks = slacks - vector(prod.parent(), len(slacks), {gg:prod})
            block_index += len(table_constructor[params])
        # scale back slacks with the one vector, the minimum is the final result
        result = min(
            [slacks[ii]/oveii for ii, oveii in \
             enumerate(one_vector_exact) if oveii!=0]
            )
        # pad the slacks, so it is all positive where it counts
        slacks -= result*one_vector_exact
        
        self.fprint("The rounded result is {}".format(result))
        
        return result, X_matrices_rounded, e_vector_rounded, \
            slacks, [phi_vector_rounded]
    
    def _round_sdp_solution_no_phi_alter(self, sdp_result, sdp_data, 
                                         table_constructor, constraints_data, 
                                         **params):
        r"""
        Round the SDP results output to get something exact.
        """
        import gc
        import time

        # set up parameters
        denom = params.get("denom", 1024)
        slack_threshold = params.get("slack_threshold", 1e-9)
        linear_threshold = params.get("linear_threshold", 1e-6)
        kernel_threshold = params.get("kernel_threshold", 1e-4)
        kernel_denom = params.get("kernel_denom", 1024)
        
        # unpack variables
        block_sizes, target_list_exact, _, __ = sdp_data
        target_vector_exact = vector(target_list_exact)
        flags_num, _, positives_list_exact, __, ___ = constraints_data
        _ = None; __ = None; ___ = None; gc.collect()
        positives_matrix_exact = matrix(
            len(positives_list_exact), flags_num, positives_list_exact
            )
        
        # find the one_vector from the equality constraint
        one_vector_exact = positives_matrix_exact.rows()[-2] 
        # remove the equality constraints
        positives_matrix_exact = positives_matrix_exact[:-2, :]

        # dim: |F_n|, c vector, primal slack for flags
        c_vector_approx = vector(sdp_result['X'][-2]) 

        c_zero_inds = [
            FF for FF, xx in enumerate(c_vector_approx) if 
            (xx<=slack_threshold)
            ]

        # same as m, number of positive constraints (-2 for the equality)
        positives_num = -block_sizes[-1] - 2 
        # dim: m, the e vector, primal slack for positivitives
        e_vector_approx = vector(sdp_result['X'][-1][:-2])
        # as above but rounded
        e_vector_rounded = vector(QQ, 
                                  _round_list(e_vector_approx, method=0, denom=denom)
                                  ) 

        # The f (ff) positivity constraints where the e vector is zero/nonzero
        e_zero_inds = [
            ff for ff, xx in enumerate(e_vector_approx) if \
                (xx<linear_threshold)
                ]
        e_nonzero_inds = [
            ff for ff in range(positives_num) if ff not in e_zero_inds
            ]



        bound_exact = _round(sdp_result['primal'], method=1)
        # the constraints for the flags that are exact
        corrected_target_relevant_exact = vector(
            [target_vector_exact[FF] - bound_exact for FF in c_zero_inds]
            )
        # the d^f_F matrix, but only the relevant parts for the rounding
        # so F where c_F = 0 and f where e_f != 0
        positives_matrix_relevant_exact = matrix(
            len(e_nonzero_inds), len(c_zero_inds), 
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
        
        self.fprint("Flattening X matrices")
        start_time = time.time()
        M_flat_relevant_matrix_exact = matrix(
            len(c_zero_inds), 0, 0, sparse=True
            )
        X_flat_list = []
        block_index = 0
        X_recover_bases = []
        X_sizes_corrected = []

        for params in table_constructor.keys():
            ns, ftype, target_size = params
            table = self.mul_project_table(
                ns, ns, ftype, ftype_inj=[], target_size=target_size
                )

            for plus_index, base in enumerate(table_constructor[params]):
                
                X_approx = matrix(sdp_result['X'][block_index + plus_index])
                X_kernel_removed, recover_base = _remove_kernel(X_approx, kernel_denom, kernel_threshold)
                X_recover_bases.append(recover_base)
                X_sizes_corrected.append(X_kernel_removed.nrows())
                X_rounded_flattened = _round_list(
                    _flatten_matrix(X_kernel_removed.rows()), method=0, denom=denom
                    )
                X_flat_list.extend(X_rounded_flattened)
                
                M_extra = []

                X_adjusted_base = recover_base.T * base
                for FF in c_zero_inds:
                    M_FF = table[FF]
                    M_extra.append(
                        _flatten_matrix(
                            (X_adjusted_base * M_FF * X_adjusted_base.T).rows(), doubled=True
                            )
                        )

                M_flat_relevant_matrix_exact = \
                    M_flat_relevant_matrix_exact.T.stack(
                        matrix(M_extra).T
                        ).T
                
                gc.collect()

            block_index += len(table_constructor[params])

        # 
        # Append the relevant M matrix and the X with the additional values from
        # the positivity constraints. Then correct the x vector values
        # 

        M_matrix_final = M_flat_relevant_matrix_exact.T.stack(
            positives_matrix_relevant_exact
            ).T
        x_vector_final = vector(
            X_flat_list + e_nonzero_list_rounded
            )
        self.fprint("This took {}s".format(time.time() - start_time))
        start_time = time.time()
        self.fprint("Correcting flat X matrices")
        self.fprint("Dimensions: ", M_matrix_final.dimensions())
        # Correct the values of the x vector, based on the minimal L_2 norm
        residual = M_matrix_final * x_vector_final - corrected_target_relevant_exact
        M_self_prod = matrix(M_matrix_final * M_matrix_final.T, sparse=False)
        x_vector_corr = x_vector_final - M_matrix_final.T * (
            M_self_prod.pseudoinverse() * residual
        )
        self.fprint("This took {}s".format(time.time() - start_time))
        start_time = time.time()

        self.fprint("Unflattening X matrices")

        #
        # Recover the X matrices and e vector from the corrected x
        #
        
        e_nonzero_vector_corr = x_vector_corr[-len(e_nonzero_inds):]
        if len(e_nonzero_vector_corr)>0 and min(e_nonzero_vector_corr)<0:
            self.fprint("Linear coefficient is negative: {}".format(
                min(e_nonzero_vector_corr)
                ))
            e_nonzero_vector_corr = [max(xx, 0) for xx in e_nonzero_vector_corr]
        e_vector_dict = dict(zip(e_nonzero_inds, e_nonzero_vector_corr))
        e_vector_base = vector(e_vector_dict.values()).base_ring()
        e_vector_corr = vector(e_vector_base, positives_num, e_vector_dict)
        self.fprint("This took {}s".format(time.time() - start_time))
        start_time = time.time()

        X_final = []
        slacks = target_vector_exact - positives_matrix_exact.T*e_vector_corr
        block_index = 0
        self.fprint("Calculating resulting bound")
        if self._printlevel > 0:
            iterator = tqdm(table_constructor.keys())
        else:
            iterator = table_constructor.keys()
        for params in iterator:
            ns, ftype, target_size = params
            table = self.mul_project_table(
                ns, ns, ftype, ftype_inj=[], target_size=target_size
                )
            for plus_index, base in enumerate(table_constructor[params]):
                block_dim = X_sizes_corrected[block_index + plus_index]
                X_ii_raw, x_vector_corr = _unflatten_matrix(
                    x_vector_corr, block_dim
                    )
                X_ii_raw = matrix(X_ii_raw)
                recover_base = X_recover_bases[block_index + plus_index]
                X_ii_small = recover_base * X_ii_raw * recover_base.T
                
                # verify semidefiniteness
                if not X_ii_small.is_positive_semidefinite():
                    self.fprint("Rounded X matrix "+ 
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
                    prod = M_flat_vector_exact*vector(
                        _flatten_matrix(X_ii_small.rows(), doubled=False)
                        )
                    slacks = slacks - vector(prod.parent(), len(slacks), {gg:prod})
                
                X_final.append(X_ii_small)
            block_index += len(table_constructor[params])

        # scale back slacks with the one vector, the minimum is the final result
        result = min([slacks[ii]/oveii \
                    for ii, oveii in enumerate(one_vector_exact) if \
                        oveii!=0])
        # pad the slacks, so it is all positive where it counts
        slacks -= result*one_vector_exact
        self.fprint("This took {}s".format(time.time() - start_time))
        start_time = time.time()

        return result, X_final, e_vector_corr, slacks, []

    def _round_sdp_solution_phi(self, sdp_result, sdp_data, 
                                table_constructor, constraints_data, 
                                phi_vectors_exact, **params):
        r"""
        Round the SDP results output to get something exact.
        """
        import gc
        import time

        # set up parameters
        denom = params.get("denom", 1024)
        slack_threshold = params.get("slack_threshold", 1e-9)
        linear_threshold = params.get("linear_threshold", 1e-6)
        kernel_threshold = params.get("kernel_threshold", 1e-4)
        kernel_denom = params.get("kernel_denom", 1024)
        
        # unpack variables
        block_sizes, target_list_exact, _, __ = sdp_data
        target_vector_exact = vector(target_list_exact)
        flags_num, _, positives_list_exact, __, ___ = constraints_data
        _ = None; __ = None; ___ = None; gc.collect()
        positives_matrix_exact = matrix(
            len(positives_list_exact), flags_num, positives_list_exact
            )
        
        no_constr = len(phi_vectors_exact)==0
        phi_vector_exact = vector(
            [0]*positives_matrix_exact.ncols()
            ) if no_constr else phi_vectors_exact[0]
        
        # find the one_vector from the equality constraint
        one_vector_exact = positives_matrix_exact.rows()[-2] 
        # remove the equality constraints
        positives_matrix_exact = positives_matrix_exact[:-2, :]

        # dim: |F_n|, c vector, primal slack for flags
        c_vector_approx = vector(sdp_result['X'][-2]) 

        c_zero_inds = [
            FF for FF, xx in enumerate(c_vector_approx) if 
            (abs(xx)<=slack_threshold or phi_vector_exact[FF]!=0)
            ]

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
                (abs(xx)<linear_threshold or phi_pos_vector_exact[ff]!=0)
                ]
        e_nonzero_inds = [
            ff for ff in range(positives_num) if ff not in e_zero_inds
            ]



        bound_exact = target_vector_exact*phi_vector_exact 
        # the constraints for the flags that are exact
        corrected_target_relevant_exact = vector( 
            [target_vector_exact[FF] - bound_exact for FF in c_zero_inds]
            )
        # the d^f_F matrix, but only the relevant parts for the rounding
        # so F where c_F = 0 and f where e_f != 0
        positives_matrix_relevant_exact = matrix(
            len(e_nonzero_inds), len(c_zero_inds), 
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
        
        self.fprint("Flattening X matrices")
        start_time = time.time()
        M_flat_relevant_matrix_exact = matrix(
            len(c_zero_inds), 0, 0, sparse=True
            )
        X_flat_list = []
        block_index = 0
        X_recover_bases = []
        X_sizes_corrected = []

        for params in table_constructor.keys():
            ns, ftype, target_size = params
            table = self.mul_project_table(
                ns, ns, ftype, ftype_inj=[], target_size=target_size
                )

            for plus_index, base in enumerate(table_constructor[params]):
                
                X_approx = matrix(sdp_result['X'][block_index + plus_index])
                X_kernel_removed, recover_base = _remove_kernel(X_approx, kernel_denom, kernel_threshold)
                X_recover_bases.append(recover_base)
                X_sizes_corrected.append(X_kernel_removed.nrows())
                X_rounded_flattened = _round_list(
                    _flatten_matrix(X_kernel_removed.rows()), method=0, denom=denom
                    )
                X_flat_list.extend(X_rounded_flattened)
                
                M_extra = []

                X_adjusted_base = recover_base.T * base
                for FF in c_zero_inds:
                    M_FF = table[FF]
                    M_extra.append(
                        _flatten_matrix(
                            (X_adjusted_base * M_FF * X_adjusted_base.T).rows(), doubled=True
                            )
                        )

                M_flat_relevant_matrix_exact = \
                    M_flat_relevant_matrix_exact.T.stack(
                        matrix(M_extra).T
                        ).T
                
                gc.collect()

            block_index += len(table_constructor[params])

        # 
        # Append the relevant M matrix and the X with the additional values from
        # the positivity constraints. Then correct the x vector values
        # 

        M_matrix_final = M_flat_relevant_matrix_exact.T.stack(
            positives_matrix_relevant_exact
            ).T
        x_vector_final = vector(
            X_flat_list + e_nonzero_list_rounded
            )
        self.fprint("This took {}s".format(time.time() - start_time))
        start_time = time.time()
        self.fprint("Correcting flat X matrices")
        self.fprint("Dimensions: ", M_matrix_final.dimensions())
        # Correct the values of the x vector, based on the minimal L_2 norm
        residual = M_matrix_final * x_vector_final - corrected_target_relevant_exact
        M_self_prod = matrix(M_matrix_final * M_matrix_final.T, sparse=False)
        x_vector_corr = x_vector_final - M_matrix_final.T * (
            M_self_prod.pseudoinverse() * residual
        )
        self.fprint("This took {}s".format(time.time() - start_time))
        start_time = time.time()

        self.fprint("Unflattening X matrices")

        #
        # Recover the X matrices and e vector from the corrected x
        #
        
        e_nonzero_vector_corr = x_vector_corr[-len(e_nonzero_inds):]
        if len(e_nonzero_vector_corr)>0 and min(e_nonzero_vector_corr)<0:
            self.fprint("Linear coefficient is negative: {}".format(
                min(e_nonzero_vector_corr)
                ))
            e_nonzero_vector_corr = [max(xx, 0) for xx in e_nonzero_vector_corr]
        e_vector_dict = dict(zip(e_nonzero_inds, e_nonzero_vector_corr))
        e_vector_base = vector(e_vector_dict.values()).base_ring()
        e_vector_corr = vector(e_vector_base, positives_num, e_vector_dict)
        self.fprint("This took {}s".format(time.time() - start_time))
        start_time = time.time()

        X_final = []
        slacks = target_vector_exact - positives_matrix_exact.T*e_vector_corr
        block_index = 0
        self.fprint("Calculating resulting bound")
        if self._printlevel > 0:
            iterator = tqdm(table_constructor.keys())
        else:
            iterator = table_constructor.keys()
        for params in iterator:
            ns, ftype, target_size = params
            table = self.mul_project_table(
                ns, ns, ftype, ftype_inj=[], target_size=target_size
                )
            for plus_index, base in enumerate(table_constructor[params]):
                block_dim = X_sizes_corrected[block_index + plus_index]
                X_ii_raw, x_vector_corr = _unflatten_matrix(
                    x_vector_corr, block_dim
                    )
                X_ii_raw = matrix(X_ii_raw)
                recover_base = X_recover_bases[block_index + plus_index]
                X_ii_small = recover_base * X_ii_raw * recover_base.T
                
                # verify semidefiniteness
                invalid = False
                try:
                    if not X_ii_small.is_positive_semidefinite():
                        self.fprint("Rounded X matrix "+ 
                            "{} is not semidefinite: {}".format(
                                block_index+plus_index, 
                                min(X_ii_small.eigenvalues())
                                ))
                        invalid = True
                except:
                    RFF = RealField(prec=100)
                    if not matrix(RFF, X_ii_small).is_positive_semidefinite():
                        self.fprint("Rounded X matrix "+ 
                            "{} is not semidefinite: {}".format(
                                block_index+plus_index, 
                                min(matrix(RFF, X_ii_small).eigenvalues())
                                ))
                        invalid = True
                if invalid:
                    return None
                
                # update slacks
                for gg, morig in enumerate(table):
                    M = base * morig * base.T
                    M_flat_vector_exact = vector(
                        _flatten_matrix(M.rows(), doubled=True)
                        )
                    prod = M_flat_vector_exact*vector(
                        _flatten_matrix(X_ii_small.rows(), doubled=False)
                        )
                    slacks = slacks - vector(prod.parent(), len(slacks), {gg:prod})
                
                X_final.append(X_ii_small)
            block_index += len(table_constructor[params])

        # scale back slacks with the one vector, the minimum is the final result
        result = min([slacks[ii]/oveii \
                    for ii, oveii in enumerate(one_vector_exact) if \
                        oveii!=0])
        # pad the slacks, so it is all positive where it counts
        slacks -= result*one_vector_exact
        self.fprint("This took {}s".format(time.time() - start_time))
        start_time = time.time()

        return result, X_final, e_vector_corr, slacks, phi_vectors_exact
    
    def _fix_X_bases(self, table_constructor, X_original):
        r"""
        Transforms the X matrices to a base that agrees with the original 
        list of flags
        
        Basically undoes the sym/asym changes and the reductions by the 
        constructions' kernels.
        """
        from fractions import Fraction
        if X_original==None:
            return None
        X_flats = []
        block_index = 0
        for params in table_constructor.keys():
            X_ii = None
            for plus_index, base in enumerate(table_constructor[params]):
                if X_ii == None:
                    X_ii = base.T * X_original[block_index + plus_index] * base
                else:
                    X_ii += base.T * X_original[block_index + plus_index] * base
            block_index += len(table_constructor[params])
            X_flat = _flatten_matrix(X_ii.rows())
            if QQ.has_coerce_map_from(X_ii.base_ring()):
                X_flat = [
                    Fraction(int(num.numerator()), int(num.denominator())) for 
                    num in X_flat
                    ]
            elif X_ii.base_ring()==RDF or X_ii.base_ring()==RR:
                X_flat = [float(num) for num in X_flat]
            else:
                X_flat = vector(X_flat)
            X_flats.append(X_flat)
        return X_flats
    
    def _format_optimizer_output(self, table_constructor, mult=1, 
                                 sdp_output=None, rounding_output=None, 
                                 file=None, target_size=0, **misc):
        r"""
        Formats the outputs to a nice certificate
        
        The result contains: the final bound, the X matrices, the linear 
        coefficients, the slacks, a guess or exact construction, the 
        list of base flags, the list of used (typed) flags
        """
        
        from fractions import Fraction

        typed_flags = {}
        for params in table_constructor.keys():
            ns, ftype, _target_size = params
            typed_flags[(ns, ftype._pythonize())] = [
                flg._pythonize() for flg in self.generate_flags(ns, ftype)
            ]
        if target_size==0:
            raise ValueError("Empty type constraints!")

        base_flags = [
            flg._pythonize() for flg in self.generate_flags(target_size)
        ]

        factor_flags = misc.get("factor_flags", None)
        if factor_flags != None:
            python_factor_flags = [
                [xx._pythonize() for xx in factors] for 
                factors in factor_flags 
            ]
        else:
            python_factor_flags = []

        def pythonize(dim, data):
            if data==None:
                return None
            if dim==0:
                try:
                    num, denom = (QQ(data)).as_integer_ratio()
                    return Fraction(int(num), int(denom))
                except:
                    return data
            return [pythonize(dim-1, xx) for xx in data]

        result = None
        X_original = None
        e_vector = None
        slacks = None
        phi_vecs = None
        target = pythonize(1, misc.get("target", None))
        positives = pythonize(2, misc.get("positives", None))
        if positives!=None:
            positives = positives[:-2]
        maximize = mult==-1

        if sdp_output!=None:
            #these should be regular float (and stay like that)
            result = sdp_output['primal']
            oresult = result
            e_vector = sdp_output['X'][-1][:-2]
            slacks = sdp_output['X'][-2]
            phi_vecs = [sdp_output['y']]
            #except this, but this is transformed back
            X_original = [matrix(xx) for xx in sdp_output['X'][:-2]]
        elif rounding_output!=None:
            oresult, X_original, e_vector, slacks, phi_vecs = rounding_output
            result = pythonize(0, oresult)
            e_vector = pythonize(1, e_vector)
            slacks = pythonize(1, slacks)
            phi_vecs = pythonize(2, phi_vecs)
        result *= int(mult)
        oresult *= mult
        Xs = self._fix_X_bases(table_constructor, X_original)
        
        cert_dict = {"result": result, 
                     "X matrices": Xs,
                     "e vector": e_vector,
                     "slack vector": slacks,
                     "phi vectors": phi_vecs,
                     "base flags": base_flags,
                     "typed flags": typed_flags,
                     "target": target,
                     "positives": positives,
                     "factor flags": python_factor_flags,
                     "maximize": maximize,
                     "target size": int(target_size) 
                    }
        
        if file!=None and file!="" and file!="notebook":
            if not file.endswith(".pickle"):
                file += ".pickle"
            with open(file, "wb") as file_handle:
                pickle.dump(cert_dict, file_handle)
        if file=="notebook":
            return cert_dict
        return oresult
    
    def _handle_sdp_params(self, **params):
        sdp_params=["axtol", "atytol", "objtol", "pinftol", "dinftol", "maxiter", 
                    "minstepfrac", "maxstepfrac", "minstepp", "minstepd", "usexzgap", 
                    "tweakgap", "affine", "printlevel", "perturbobj", "fastmode"]
        
        if "precision" in params and params["precision"]!=None:
            precision = params["precision"]
            if "axtol" not in params:
                params["axtol"] = precision
            if "atytol" not in params:
                params["atytol"] = precision
            if "objtol" not in params:
                params["objtol"] = precision
            if "pinftol" not in params:
                params["pinftol"] = 1/precision
            if "dinftol" not in params:
                params["dinftol"] = 1/precision
            if "minstepp" not in params:
                params["minstepp"] = precision
            if "minstepd" not in params:
                params["minstepd"] = precision
        
        if self._printlevel != 1:
            if params=={}:
                params = {"printlevel": self._printlevel}
            else:
                if "printlevel" not in params:
                    params["printlevel"] = self._printlevel
        if params!={}:
            with open("param.csdp", "w") as paramsfile:
                for key, value in params.items():
                    if str(key) not in sdp_params:
                        continue
                    paramsfile.write(f"{key}={value}\n")
            if "printlevel" in params:
                self._printlevel = params["printlevel"]
            else:
                self._printlevel = 1

    def solve_sdp(self, target_element, target_size, construction, 
                  maximize=True, positives=None, file=None, 
                  specific_ftype=None, **params):
        r"""
        TODO Docstring
        """
        from csdpy import solve_sdp
        import time

        #
        # Initial setup
        #

        self._handle_sdp_params(**params)
        base_flags = self.generate_flags(target_size)
        self.fprint("Base flags generated, their number is {}".format(
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
        
        if specific_ftype==None:
            ftype_data = self._get_relevant_ftypes(target_size)
        else:
            ftype_data = specific_ftype
        self.fprint("The relevant ftypes are constructed, their " + 
              "number is {}".format(len(ftype_data)))
        flags = [self.generate_flags(dat[0], dat[1]) for dat in ftype_data]
        flag_sizes = [len(xx) for xx in flags]
        self.fprint("Block sizes before symmetric/asymmetric change is" + 
              " applied: {}".format(flag_sizes))
        
        #
        # Create the table constructor and adjust it based on construction
        #
        
        table_constructor = self._create_table_constructor(
            ftype_data, target_size
            )
        if isinstance(construction, FlagAlgebraElement):
            phi_vectors_exact = [construction.values()]
        else:
            phi_vectors_exact = [xx.values() for xx in construction]
        self.fprint("Adjusting table with kernels from construction")
        table_constructor = self._adjust_table_phi(
            table_constructor, phi_vectors_exact
            )
        sdp_data = self._tables_to_sdp_data(
            table_constructor, prev_data=sdp_data
            )
        self.fprint("Tables finished")

        #
        # Add constraints data and add it to sdp_data
        #
        
        constraints_data = self._create_constraints_data(
            positives, target_element, target_size
            )
        sdp_data = self._constraints_to_sdp_data(
            constraints_data, prev_data=sdp_data
            )
        self.fprint("Constraints finished")
        
        #
        # Then run the optimizer
        #
        
        self.fprint("Running SDP. Used block sizes are {}".format(sdp_data[0]))
        time.sleep(float(0.1))
        final_sol = solve_sdp(*sdp_data)
        time.sleep(float(0.1))

        if file!=None:
            if not file.endswith(".pickle"):
                file += ".pickle"
            with open(file, "wb") as file_handle:
                save_data = (
                    final_sol, 
                    (sdp_data[0], sdp_data[1], None, None), 
                    table_constructor, 
                    (constraints_data[0], None, constraints_data[2], None, constraints_data[4]), 
                    phi_vectors_exact, 
                    mult, target_size)
                pickle.dump(save_data, file_handle)

        return final_sol["primal"]*mult

    def round_solution(self, sdp_output_file, certificate_file=None, **params):
        r"""
        TODO Docstring
        """
        if not sdp_output_file.endswith(".pickle"):
            sdp_output_file += ".pickle"
        with open(sdp_output_file, "rb") as file_handle:
            save_data = pickle.load(file_handle)
        #
        # Unpack the data
        #
        
        sdp_result, sdp_data, table_constructor, \
            constraints_data, phi_vectors_exact, \
            mult, target_size = save_data
        

        #
        # Perform the rounding
        #

        self.fprint("Starting the rounding of the result")
        rounding_output = self._round_sdp_solution_phi( \
            sdp_result,
            sdp_data, 
            table_constructor,
            constraints_data,
            phi_vectors_exact,
            **params
            )
        if rounding_output==None:
            print("Rounding was unsuccessful!")
            return
        value = rounding_output[0]*mult
        self.fprint("Final rounded bound is {}".format(value))
        
        if certificate_file==None:
            return value
        return self._format_optimizer_output(
            table_constructor, 
            mult=mult,
            rounding_output=rounding_output,
            file=certificate_file,
            target=vector(sdp_data[1])*mult,
            positives=constraints_data[2],
            factor_flags=constraints_data[4],
            target_size=target_size
            )

    def optimize_problem(self, target_element, target_size, maximize=True, 
                         positives=None, construction=None, file=None, 
                         exact=False, specific_ftype=None, 
                         **params):
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

        self._handle_sdp_params(**params)
        
        constr_print_limit = params.get("constr_print_limit", 1000)
        constr_error_threshold = params.get("constr_error_threshold", 1e-6)
        
        if target_size not in self.sizes():
            raise ValueError("For theory {}, size {} is not allowed.".format(self._name, target_size))

        base_flags = self.generate_flags(target_size)
        self.fprint("Base flags generated, their number is {}".format(
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
        
        if specific_ftype==None:
            ftype_data = self._get_relevant_ftypes(target_size)
        else:
            ftype_data = specific_ftype
        self.fprint("The relevant ftypes are constructed, their " + 
              "number is {}".format(len(ftype_data)))
        flags = [self.generate_flags(dat[0], dat[1]) for dat in ftype_data]
        flag_sizes = [len(xx) for xx in flags]
        self.fprint("Block sizes before symmetric/asymmetric change is" + 
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
        self.fprint("Tables finished")

        #
        # Add constraints data and add it to sdp_data
        #
        
        constraints_data = self._create_constraints_data(
            positives, target_element, target_size
            )
        sdp_data = self._constraints_to_sdp_data(
            constraints_data, prev_data=sdp_data
            )
        self.fprint("Constraints finished")
        
        #
        # Helper for returning the data, writing to file, and some cleanup
        #

        def help_return(value, sdpo=None, roundo=None):
            try:
                os.remove("param.csdp")
            except OSError:
                pass
            if file==None:
                return value
            return self._format_optimizer_output(
                table_constructor, 
                mult=mult, 
                sdp_output=sdpo, 
                rounding_output=roundo,
                file=file,
                target=target_vector_exact*mult,
                positives=constraints_data[2],
                factor_flags=constraints_data[4],
                target_size=target_size
                )

        #
        # If construction is None or [] then run the optimizer 
        # without any construction
        #
        
        if construction==None or construction==[] or construction is False:
            self.fprint("Running sdp without construction. " + 
                  "Used block sizes are {}".format(sdp_data[0]))
            
            time.sleep(float(0.1))
            initial_sol = solve_sdp(*sdp_data)
            time.sleep(float(0.1))

            # Format the result and return it if floating point values are fine
            if (not exact):
                return help_return(initial_sol['primal'] * mult, 
                                   sdpo=initial_sol)
                
            # Guess the construction in this case
            if construction==None:
                one_vector = constraints_data[3]
                phi_vector_rounded, error_coeff = _round_adaptive(
                    initial_sol['y'], one_vector
                    )
                if error_coeff<constr_error_threshold:
                    alg = FlagAlgebra(self, QQ)
                    construction = alg(target_size, phi_vector_rounded)
                    phipr = str(construction)
                    self.fprint("The initial run gave an accurate "+
                          "looking construction")
                    if len(phipr)<constr_print_limit:
                        self.fprint("Rounded construction vector "+
                              "is: \n{}".format(phipr))
                else:
                    self.fprint("The initial run didn't provide an "+
                          "accurate construction")
                    construction = []
            # If nothing was provided or the guess failed, 
            # then round the current solution
            if construction==[]:
                try:
                    rounding_output = self._round_sdp_solution_no_phi_alter(
                        initial_sol, sdp_data, table_constructor, 
                        constraints_data, **params)
                except:
                    rounding_output = None
                    construction = False
                
                if rounding_output!=None:
                    return help_return(rounding_output[0] * mult, roundo=rounding_output)
                else:
                    construction = False
            
            if construction==False:
                rounding_output = self._round_sdp_solution_no_phi(
                    initial_sol, sdp_data, table_constructor, 
                    constraints_data, **params)
                return help_return(rounding_output[0] * mult, roundo=rounding_output)
        
        
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

        self.fprint("Adjusting table with kernels from construction")
        table_constructor = self._adjust_table_phi(
            table_constructor, phi_vectors_exact
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
        
        self.fprint("Running SDP after kernel correction. "+
              "Used block sizes are {}".format(sdp_data[0]))
        time.sleep(float(0.1))
        final_sol = solve_sdp(*sdp_data)
        time.sleep(float(0.1))

        # Quickly deal with the case when no rounding is needed
        if (not exact):
            return help_return(final_sol['primal'] * mult, sdpo=final_sol)
        
        self.fprint("Starting the rounding of the result")
        rounding_output = self._round_sdp_solution_phi(
            final_sol, sdp_data, table_constructor, 
            constraints_data, phi_vectors_exact, 
            **params
            )
        if rounding_output==None:
            self.fprint("Rounding based on construction was unsuccessful")
            rounding_output = self._round_sdp_solution_no_phi(
                final_sol, sdp_data, table_constructor, 
                constraints_data, **params
                )
        
        self.fprint("Final rounded bound is {}".format(rounding_output[0]*mult))
        
        return help_return(rounding_output[0] * mult, roundo=rounding_output)
    
    optimize = optimize_problem
    
    def fprint(self, *args):
        if self._printlevel != 0:
            print(*args)

    def external_optimize(self, target_element, target_size, maximize=True, 
                          positives=None, construction=None, file=None, 
                          specific_ftype=None, **params):
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
        self.fprint("Base flags generated, their number "+
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
        self.fprint("The relevant ftypes are constructed, "+
              "their number is {}".format(len(ftype_data)))
        flags = [self.generate_flags(dat[0], dat[1]) for dat in ftype_data]
        flag_sizes = [len(xx) for xx in flags]
        self.fprint("Block sizes before symmetric/asymmetric change " + 
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
            self.fprint("Adjusting table with kernels from construction")
            table_constructor = self._adjust_table_phi(
                table_constructor, phi_vectors_exact
                )
        sdp_data = self._tables_to_sdp_data(
            table_constructor, prev_data=sdp_data
            )
        self.fprint("Tables finished")

        #
        # Create constraints data and add it to sdp_data
        #
        
        constraints_data = self._create_constraints_data(
            positives, target_element, target_size
            )
        sdp_data = self._constraints_to_sdp_data(
            constraints_data, prev_data=sdp_data
            )
        self.fprint("Constraints finished")
        #
        # Make sdp data integer and write it to a file
        #
        precision = params.get("precision", 20)
        R = RealField(prec=round(precision*log(10, 2)), sci_not=True)
        
        with open(file, "w") as file:
            block_sizes, target, mat_inds, mat_vals = sdp_data
            file.write("{}\n{}\n".format(len(target), len(block_sizes)))
            file.write(" ".join(map(str, block_sizes)) + "\n")
            for xx in target:
                file.write(str(R(xx)) + " ")
            file.write("\n")

            for ii in range(len(mat_vals)):
                file.write("{} {} {} {} {}\n".format(
                    mat_inds[ii*4 + 0],
                    mat_inds[ii*4 + 1],
                    mat_inds[ii*4 + 2],
                    mat_inds[ii*4 + 3],
                    str(R(mat_vals[ii]))
                ))

    def construction_from_certificate(self, file_or_cert):
        if isinstance(file_or_cert, str):
            file = file_or_cert
            if not file.endswith(".pickle"):
                file += ".pickle"
            with open(file, 'rb') as file:
                certificate = pickle.load(file)
        else:
            certificate = file_or_cert
        conss = certificate["phi vectors"]
        if "target size" not in certificate:
            raise ValueError("The certificate contains no target size")
        tsize = certificate["target size"]
        if len(conss)==0:
            raise ValueError("The certificate contains no constructions")
        cons = conss[0]
        if len(cons)==0:
            return 0
        if isinstance(cons[0], float):
            R = RR
            vals = cons
        else:
            R = QQ
            vals = []
            for ii, xx in enumerate(cons):
                try:
                    vals.append(R(xx))
                except:
                    R = xx.parent()
                    for jj in range(ii):
                        vals[jj] = R(vals[jj])
        FA = FlagAlgebra(self, R)
        return FA(tsize, vals)

    def verify_certificate(self, file_or_cert, target_element=None, 
                           target_size=None, maximize=True, positives=None, 
                           **params):
        r"""
        Verifies the certificate provided by the optimizer 
        written to `file`
        """
        
        if self._printlevel != 1:
            if "printlevel" in params:
                self._printlevel = params["printlevel"]

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
        
        from fractions import Fraction
        def to_sage(dim, data):
            if dim==0:
                if isinstance(data, Fraction):
                    return QQ(data)
                if isinstance(data, float):
                    return RR(data)
                return data
            return [to_sage(dim-1, xx) for xx in data]

        #
        # Checking eigenvalues and positivity constraints
        #
        
        e_values = to_sage(1, certificate["e vector"])

        if len(e_values)>0 and min(e_values)<0:
            print("Solution is not valid!")
            self.fprint("Linear constraint's coefficient is negative ", 
                  min(e_values))
            return -1
        
        X_flats = to_sage(2, certificate["X matrices"])
        self.fprint("Checking X matrices")
        if self._printlevel > 0:
            iterator = tqdm(enumerate(X_flats))
        else:
            iterator = enumerate(X_flats)
        for ii, Xf in iterator:
            X = matrix(_unflatten_matrix(Xf)[0])
            try:
                if not (X.is_positive_semidefinite()):
                    print("Solution is not valid!")
                    self.fprint("Matrix {} is not semidefinite".format(ii))
                    return -1
            except:
                RFF = RealField(prec=100)
                if not matrix(RFF, X).is_positive_semidefinite():
                    print("Solution is not valid!")
                    self.fprint("Matrix {} is not semidefinite".format(ii))
                    return -1

        self.fprint("Solution matrices are all positive semidefinite, " + 
              "linear coefficients are all non-negative")

        #
        # Initial setup
        #

        if target_size==None:
            if "target size" in certificate:
                target_size = to_sage(0, certificate["target size"])
            else:
                raise ValueError("Target size must be specified "+
                                 "if it is not part of the certificate!")

        mult = -1 if maximize else 1
        base_flags = certificate["base flags"]
        one_vector = vector([1]*len(base_flags))
        if target_element==None:
            if "target" in certificate:
                target_vector_exact = vector(to_sage(1, certificate["target"]))
            else:
                raise ValueError("Target must be specified "+
                                 "if it is not part of the certificate!")
        else:
            target_vector_exact = (
                target_element.project()<<\
                    (target_size - target_element.size())
                    ).values()
            if not (target_element.ftype().afae()==1):
                one_vector = (
                    target_element.ftype().project()<<\
                        (target_size - target_element.ftype().size())
                        ).values()
        target_vector_exact *= mult
        ftype_data = list(certificate["typed flags"].keys())
        
        #
        # Create the semidefinite matrix data
        #

        table_list = []
        self.fprint("Calculating multiplication tables")
        if self._printlevel > 0:
            iterator = tqdm(enumerate(ftype_data))
        else:
            iterator = enumerate(ftype_data)
        for ii, dat in iterator:
            ns, ftype = dat
            ftype = self._element_constructor_(ftype[0], ftype=ftype[1], **dict(ftype[2]))
            #calculate the table
            table = self.mul_project_table(
                ns, ns, ftype, ftype_inj=[], target_size=target_size
                )
            if table!=None:
                table_list.append(table)

        #
        # Create the data from linear constraints
        #

        if positives != None:
            positives_list_exact = []
            for ii, fv in enumerate(positives):
                if isinstance(fv, BuiltFlag) or isinstance(fv, ExoticFlag):
                    continue
                nf = fv.size()
                df = target_size + fv.ftype().size() - nf
                mult_table = self.mul_project_table(
                    nf, df, fv.ftype(), ftype_inj=[], target_size=target_size
                    )
                fvvals = fv.values()
                m = matrix([vector(fvvals*mat) for mat in mult_table])
                positives_list_exact += list(m.T)
                self.fprint("Done with positivity constraint {}".format(ii))
            positives_matrix_exact = matrix(
                len(positives_list_exact), len(base_flags), 
                positives_list_exact
                )
            e_values = vector(e_values[:len(positives_list_exact)])
        else:
            if "positives" in certificate:
                posls = to_sage(2, certificate["positives"])[:-2]
                positives_matrix_exact = matrix(len(posls), len(base_flags), posls)
                e_values = vector(e_values[:len(posls)])
            else:
                positives_matrix_exact = matrix(0, len(base_flags), [])
                e_values = vector(QQ, [])

        self.fprint("Done calculating linear constraints")

        #
        # Calculate the bound the solution provides
        #
        
        self.fprint("Calculating the bound provided by the certificate")
        
        slacks = target_vector_exact - positives_matrix_exact.T*e_values
        if self._printlevel > 0:
            iterator = tqdm(enumerate(table_list))
        else:
            iterator = enumerate(table_list)
        for ii, table in iterator:
            slvecdel = []
            for mat_gg in table:
                mat_flat = vector(
                    _flatten_matrix(mat_gg.rows(), doubled=True)
                    )
                slvecdel.append(mat_flat * vector(X_flats[ii]))
            slacks -= vector(slvecdel)
        
        result = min(
            [slacks[ii]/oveii for ii, oveii in enumerate(one_vector) \
             if oveii!=0]
            )
        result *= mult
        print("The solution is valid, it proves "+
              "the bound {}".format(result))

        return result
    
    verify = verify_certificate
    

    # Generating flags

    def _find_ftypes(self, empstrs, ftype):
        import multiprocessing as mp
        pool = mp.Pool(mp.cpu_count()-1)
        pares = pool.map(ftype.ftypes_inside, empstrs)
        pool.close(); pool.join()
        return tuple(itertools.chain.from_iterable(pares))

    # Generating tables

    def sym_asym_bases(self, n, ftype=None):
        r"""
        Generate the change of base matrices for the symmetric
        and the asymmetric subspaces
        """

        flags = self.generate_flags(n, ftype)
        uniques = []
        sym_base = []
        asym_base = []

        #Correct method
        sym_base_lasts = []
        for xx in flags:
            xxid = xx.unique(weak=True)[0]
            if xxid not in uniques:
                uniques.append(xxid)
                sym_base.append(xx.afae())
                sym_base_lasts.append(xx.afae())
            else:
                sym_ind = uniques.index(xxid)
                asym_base.append(sym_base_lasts[sym_ind] - xx)
                sym_base[sym_ind] += xx
                sym_base_lasts[sym_ind] = xx

        #Old worse method (but sometimes helps rounding)
        # for xx in flags:
        #     xxid = xx.unique(weak=True)[0]
        #     if xxid not in uniques:
        #         uniques.append(xxid)
        #         sym_base.append(xx.afae())
        #     else:
        #         sym_ind = uniques.index(xxid)
        #         asym_base.append(sym_base[sym_ind] - xx.afae())
        #         sym_base[sym_ind] += xx
        
        
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

class BuiltTheory(_CombinatorialTheory):
    
    Element = BuiltFlag
    
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
        self._no_question = True
        _CombinatorialTheory.__init__(self, name)
    
    # Parent methods
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

        if isinstance(n, BuiltFlag) or isinstance(n, Pattern):
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
    
    # Optimizing and rounding

    def _get_relevant_ftypes(self, target_size):
        r"""
        Returns the relevant ftypes for optimizing up to ``target_size``.

        INPUT:

            - ``target_size`` -- integer; the target size for which the ftypes are generated.

        OUTPUT:
            list; a list of tuples where each tuple is of the form (ns, ftype, target_size).

        EXAMPLES:

            sage: GraphTheory._get_relevant_ftypes(3)
            asd

        TESTS:

            sage: GraphTheory._get_relevant_ftypes(-1)
            asd
        """
        if target_size < 0:
            raise ValueError("Target size should be non-negative!")
        
        plausible_sizes = list(range(1, target_size))
        ftype_pairs = []
        for fs, ns in itertools.combinations(plausible_sizes, r=int(2)):
            if ns + ns - fs <= target_size:
                kk = ns - fs
                found = False
                for ii, (bfs, bns) in enumerate(ftype_pairs):
                    if bns - bfs == kk:
                        found = True
                        if ns > bns:
                            ftype_pairs[ii] = (fs, ns)
                        break
                if not found:
                    ftype_pairs.append((fs, ns))

        ftype_data = []
        for fs, ns in ftype_pairs:
            ftype_flags = self.generate_flags(fs)
            ftypes = [flag.subflag([], ftype_points=list(range(fs)))
                    for flag in ftype_flags]
            for xx in ftypes:
                ftype_data.append((ns, xx, target_size))
        ftype_data.sort()
        return ftype_data

    # Generating flags
    
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
    
    def no_question(self, val=None):
        if val==None:
            val = True
        self._no_question = val

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
            if run_bound==infinity or guess < run_bound or n<=3 or self._no_question:
                ret = self._find_ftypes(empstrs, ftype)
            else:
                confirm = input("This might take a while: " + 
                                "{}. Continue? y/n\n".format(guess))
                if "y" in confirm.lower():
                    ret = self.generate_flags(n, ftype, run_bound=infinity)
                else:
                    raise RuntimeError("Calculation interrupted")
        self._save(ret, key)
        return ret
    
    generate = generate_flags

    def exclude(self, structs=None, force=False):
        #Set up structs to contain all we want to exclude
        if structs==None:
            structs = []
        elif type(structs)==BuiltFlag or type(structs)==Pattern:
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
            elif isinstance(xx, BuiltFlag):
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

    def get_total_excluded(self, n=100):
        if self._sources == None:
            ret = [xx for xx in self._excluded if xx.size()<=n]
        else:
            ret = [xx for xx in self._excluded if xx.size()<=n]
            ret += list(self._sources[0].get_total_excluded(n))
            ret += list(self._sources[1].get_total_excluded(n))
        return tuple(ret)

    def match_pattern(self, pattern):
        if pattern is BuiltFlag:
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
    
    
    # Generating tables

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

class ExoticTheory(_CombinatorialTheory):
    
    Element = ExoticFlag
    
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
            self._sizes = [ii for ii in range(1000) if size_combine(0, ii, 0) == ii]
        self._generator = generator
        self._identifier = identifier
        self._sources = None
        self._symmetries = None
        _CombinatorialTheory.__init__(self, name)
    
    # Parent methods

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

        if isinstance(n, ExoticFlag):
            return n

        if n not in self.sizes():
            ValueError("For theory {}, size {} is not allowed.".format(self._name, n))

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
            if xx in kwds:
                try:
                    blocks[xx] = tuple(kwds[xx])
                except:
                    raise ValueError("The provided {} must be iterable".format(xx))
                
        return self.element_class(self, n, ftype_points, **blocks)

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
        return self._element_constructor_(0)
    
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
    
    def pattern(self, *args, **kwds):
        raise NotImplementedError("Patterns are not implemented for" + 
                                  str(self))
    
    p = pattern
    P = pattern
    Pattern = pattern

    # Optimizing and rounding

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

    # Generating flags

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

    def exclude(self, structs=None, force=False):
        #Set up structs to contain all we want to exclude
        if structs==None:
            structs = []
        elif type(structs)==ExoticFlag:
            structs = [structs]
        if not force:
            structs += list(self._excluded)
        
        #Make structs sorted, so it is as small as possible
        structs.sort(key=lambda x : x.size())
        self._excluded = tuple()

        for xx in structs:
            if isinstance(xx, ExoticFlag):
                if xx.theory()!=self:
                    print("{} is from a different theory".format(xx))
                else:
                    if xx in self.generate(xx.size()):
                        self._excluded = tuple(list(self._excluded) + [xx])
            else:
                print("Excluding objects with type {} is not supported".format(type(xx)))

    def reset_exclude(self):
        self.exclude(force=True)

    reset = reset_exclude

    def get_total_excluded(self, n):
        ret = [xx for xx in self._excluded if xx.size()<=n]
        return tuple(ret)
    
    def _check_excluded(self, elms):
        r"""
        Helper to check the excluded structures in generation
        """
        flg = elms[0]
        for xx in elms[1]:
            if xx <= flg:
                return False
        return True
    
    def _gfe(self, excluded, n, ftype):
        r"""
        Cached version of generate flags excluded

        .. SEEALSO::

            :func:`generate_flags`
        """
        if ftype==None:
            ftype = self.empty()
        
        key = ("generate", n, 
               self._serialize(excluded), ftype._serialize())
        loaded = self._load(key=key)
        if loaded != None:
            return loaded
        
        import multiprocessing as mp
        
        if ftype.size()==0: #just generate empty elements
            if n==0:
                ret = (self.empty_element(), )
            elif len(excluded)==0: #just return the output of the generator
                ret = tuple([self.element_class(self, n, tuple(), **xx) for xx in self._generator(n)])
            else:
                #otherwise check each generated for the excluded values
                slist = [(xx, excluded) for xx in self._gfe(tuple(), n, None)]
                pool = mp.Pool(mp.cpu_count()-1)
                canincl = pool.map(self._check_excluded, slist)
                pool.close(); pool.join()
                ret = tuple([slist[ii][0] for ii in range(len(slist)) if canincl[ii]])
        else:
            #generate flags by first getting the empty structures then finding the flags
            empstrs = self._gfe(excluded, n, None)
            pool = mp.Pool(mp.cpu_count()-1)
            pares = pool.map(ftype.ftypes_inside, empstrs)
            pool.close(); pool.join()
            ret = []
            for coll in pares:
                for xx in coll:
                    if xx not in ret:
                        ret.append(xx)
            ret = tuple(ret)
        self._save(ret, key=key)
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

    # Generating tables
    
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
        #Trying to load
        key = ("table", (n1, n2, N), 
               large_ftype._serialize(), tuple(ftype_inj), 
               self._serialize(excluded))
        loaded = self._load(key=key)
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
        
        self._save(ret, key=key)
        return ret
    

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
        if not isinstance(theory, BuiltTheory):
            raise ValueError("Only built theories are allowed!")
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
    
    def _serialize_data(data):
        #For the signature
        signature = data["signature"]
        sered_signature = []
        for xx in signature:
            ll = tuple(signature[xx].values())
            sered_signature.append((xx, ll))
        sered_signature = tuple(sered_signature)
        return (sered_signature, tuple(data["symmetries"]))
    
    ser_data = _serialize_data(theory_data)
    if len(theories)==2 and not symmetries:
        ser_data = (tuple(theories), ser_data)
    else:
        ser_data = (None, ser_data)
    ret_theory = BuiltTheory(name, _from_data=ser_data)
    return ret_theory

# Pre-defined theories
GraphTheory = BuiltTheory("Graph")
DiGraphTheory = BuiltTheory("DiGraph", arity=2, is_ordered=True)
ThreeGraphTheory = BuiltTheory("ThreeGraph", arity=3)
DiThreeGraphTheory = BuiltTheory("DiThreeGraph", arity=3, is_ordered=True)
FourGraphTheory = BuiltTheory("FourGraph", arity=4)
Color0 = BuiltTheory("Color0", relation_name="C0", arity=1)
Color1 = BuiltTheory("Color1", relation_name="C1", arity=1)
Color2 = BuiltTheory("Color2", relation_name="C2", arity=1)
Color3 = BuiltTheory("Color3", relation_name="C3", arity=1)
Color4 = BuiltTheory("Color4", relation_name="C4", arity=1)
Color5 = BuiltTheory("Color5", relation_name="C5", arity=1)
Color6 = BuiltTheory("Color6", relation_name="C6", arity=1)
Color7 = BuiltTheory("Color7", relation_name="C7", arity=1)

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

from sage.misc.functional import log
from sage.graphs.graph import Graph

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
    from sage.graphs.graph_generators import graphs
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

def _cube_graphs(d, edge_num):
    from sage.features.nauty import NautyExecutable
    import subprocess, select
    import shlex
    from sage.graphs.graph_generators import graphs
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
    from sage.graphs.graph_generators import graphs
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

PermutationTheory = ExoticTheory('Permutation', 
                                        _generator_permutation, 
                                        _identify_permutation, 
                                        edges=2)

OEGraphTheory = ExoticTheory('OEdgeGraph', 
                                    _generator_oe_graph, 
                                    _identify_oe_graph, 
                                    edges=2)

OVGraphTheory = ExoticTheory('OVertexGraph', 
                                    _generator_ov_graph, 
                                    _identify_ov_graph, 
                                    edges=2)

# should only be used up to size 8.
HypercubeGraphTheory = ExoticTheory('HypercubeGraph', 
                                     _generator_cube_graphs, 
                                     _identify_cube_graphs, 
                                     size_combine=_cube_size_combine, 
                                     edges=2)

# should only be used up to size 16.
HypercubeVertexTheory = ExoticTheory('HypercubeVertex', 
                                         _generator_cube_points, 
                                         _identify_cube_points, 
                                         size_combine=_cube_size_combine, 
                                         edges=2, points=1)

def Theory(name, *args, **kwdargs):
    if name=="Permutation":
        return PermutationTheory
    if name=="OEdgeGraph":
        return OEGraphTheory
    if name=="OVertexGraph":
        return OVGraphTheory
    if name=="HypercubeGraph":
        return HypercubeGraphTheory
    if name=="HypercubeVertex":
        return HypercubeVertexTheory
    return BuiltTheory(name, *args, **kwdargs)