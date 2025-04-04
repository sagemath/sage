r"""
Implementation of Flag, elements of :class:`CombinatorialTheory`

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
from sage.all import QQ
from cysignals.signals cimport sig_check
from sage.structure.element cimport Element
from sage.structure.coerce cimport coercion_model
from blisspy cimport canonical_form_from_edge_list, automorphism_group_gens_from_edge_list
from tqdm import tqdm

# Elementary block operations
cdef tuple _subblock_helper(tuple points, tuple block, bint inverse = False):
    if len(block)==0:
        return tuple()
    cdef set points_set = set(points)
    cdef dict points_index = {p: i for i, p in enumerate(points)}
    if inverse:
        points_index = {i: p for i, p in enumerate(points)}
    cdef list ret = []
    cdef bint gd
    cdef int yy
    for xx in block:
        gd = True
        for yy in xx:
            if yy not in points_set:
                gd = False
                break
        if gd:
            ret.append(tuple([points_index[yy] for yy in xx]))
    return tuple(ret)

cdef dict _merge_blocks(dict block0, dict block1, tuple only_include):
    cdef dict merged = {}
    cdef str key
    cdef tuple tp
    cdef int xx
    for key in block0:
        merged[key] = tuple(list(block0[key]) + [
            tp for tp in block1[key] if all([(xx in tp) for xx in only_include])
        ])
    return merged

cdef dict _perm_blocks(dict blocks, tuple perm, bint inverse = False):
    cdef dict ret
    cdef str xx
    ret = {
        xx: _subblock_helper(perm, blocks[xx], inverse) 
        for xx in blocks.keys()
    }
    return ret

cdef dict _standardize_blocks(dict blocks, dict signature, bint pattern):
    cdef dict ret = {}
    cdef str xx, kk
    for xx in blocks:
        if not pattern:
            if signature[xx]["ordered"]:
                ret[xx] = tuple(sorted([tuple(yy) for yy in blocks[xx]]))
            else:
                ret[xx] = tuple(sorted([tuple(sorted(yy)) for yy in blocks[xx]]))
        else:
            kk = xx
            if xx not in signature:
                kk = xx[:-2]
            if signature[kk]["ordered"]:
                ret[xx] = tuple(sorted([tuple(yy) for yy in blocks[xx]]))
            else:
                ret[xx] = tuple(sorted([tuple(sorted(yy)) for yy in blocks[xx]]))
    return ret

cdef dict _perm_signature(dict blocks, tuple perm):
    cdef dict ret = {}
    cdef int ii
    cdef str xx
    for ii, xx in enumerate(blocks.keys()):
        ret[xx] = blocks[perm[ii]]
    return ret

cdef dict _perm_pattern_signature(dict blocks, tuple perm, dict orig_sign):
    cdef dict ret = {}
    cdef int ii
    cdef str xx
    for ii, xx in enumerate(orig_sign.keys()):
        ret[xx] = blocks[perm[ii]]
        ret[xx+"_m"] = blocks[perm[ii]+"_m"]
    return ret

# Group operations

cdef set _generate_group(tuple generators, int n):
    cdef set group = set()
    cdef list to_check = [tuple(range(n))]
    cdef tuple perm, gen, new_perm
    cdef int i
    while to_check:
        perm = to_check.pop()
        if perm in group:
            continue
        group.add(perm)
        for gen in generators:
            new_perm = tuple(gen[perm[i]] for i in range(n))
            if new_perm not in group:
                to_check.append(new_perm)
    return group

cdef list _compute_coset_reps(set G, set H, int n):
    cdef list coset_reps = []
    cdef tuple g, c, h_tuple
    cdef int i
    cdef bint in_existing_coset
    cdef list h
    for g in G:
        in_existing_coset = False
        for c in coset_reps:
            h = [0] * n
            for i in range(n):
                h[c[i]] = g[i]
            h_tuple = tuple(h)
            if h_tuple in H:
                in_existing_coset = True
                break
        if not in_existing_coset:
            coset_reps.append(g)
    return coset_reps

cdef inline tuple _compose2(tuple a, tuple b, tuple c, int n):
    cdef int ii
    cdef list out_list = [0]*n
    for ii in range(n):
        out_list[ii] = a[b[c[ii]]]
    return tuple(out_list)

cdef inline tuple _compose(tuple a, tuple b, int n):
    cdef int ii
    cdef list out_list = [0]*n
    for ii in range(n):
        out_list[ii] = a[b[ii]]
    return tuple(out_list)

cdef tuple _invert_c(tuple p, int n):
    cdef int ii
    cdef list inv_list = [0]*n
    for ii in range(n):
        inv_list[p[ii]] = ii
    return tuple(inv_list)

cdef set _generate_coset(set G, tuple sigma, int n):
    cdef set coset = set()
    cdef tuple g
    for g in G:
        coset.add(_compose(g, sigma, n))
    return coset

cdef dict _left_coset_representative_map(set G0, int n):
    cdef dict rep_map = {}
    cdef set visited = set()
    cdef tuple sigma
    cdef set coset
    for sigma in itertools.permutations(range(n)):
        if sigma in visited:
            continue
        rep_map[sigma] = sigma
        coset = _generate_coset(G0, sigma, n)
        for c in coset:
            visited.add(c)
            rep_map[c] = sigma
    return rep_map

cdef set _find_double_coset_reps(set G0, set G1, dict rep_map, int n):
    cdef set T = set(rep_map.values())
    cdef set double_coset_reps = set()
    cdef tuple t, r, g1
    cdef tuple g1_inv, r_inv, candidate_g0
    cdef bint is_new
    
    for t in T:
        is_new = True
        for r in double_coset_reps:
            for g1 in G1:
                g1_inv = _invert_c(g1, n)
                r_inv = _invert_c(r, n)
                candidate_g0 = _compose2(t, g1_inv, r_inv, n)
                if candidate_g0 in G0:
                    is_new = False
                    break
            if not is_new:
                break
        if is_new:
            double_coset_reps.add(t)
    return double_coset_reps

# Helpers for the generators

cdef int _get_max_arity(dict signature, int n=1000):
    cdef int max_arity = 0
    cdef int curr_arity
    cdef str xx
    for xx in signature:
        curr_arity = signature[xx]['arity']
        if curr_arity>n:
            continue
        max_arity = max(max_arity, curr_arity)
    return max_arity

cpdef tuple _get_single_extensions(dict relation, int n, int fix):
    cdef int arity = relation['arity']
    if n<fix or fix>arity or n<arity:
        return tuple([tuple()])
    cdef bint ordered = relation['ordered']
    cdef tuple extension_tuples = tuple(itertools.combinations(range(fix, n), r=arity-fix))
    cdef tuple ord_base
    cdef tuple xx
    cdef int r

    if ordered:
        ord_base = tuple(
            itertools.chain.from_iterable([
                itertools.permutations(
                    list(range(fix))+list(xx)
                    ) for xx in extension_tuples
                ])
            )
    else:
        ord_base = tuple([
            tuple(
                list(range(fix))+list(xx)
                ) for xx in extension_tuples
        ])
    
    return tuple(itertools.chain.from_iterable(
        [itertools.combinations(ord_base, r) for r in range(len(ord_base) + 1)]
    ))

cdef tuple _get_all_extensions(dict signature, int n, int fix):
    cdef str xx

    cdef list terms = [
        _get_single_extensions(signature[xx], n, fix) for xx in signature
    ]
    
    cdef tuple poss
    cdef list ret = []
    for poss in itertools.product(*terms):
        ret.append({xx: poss[ii] for ii,xx in enumerate(signature)})
    return tuple(ret)

cdef bint _excluded_compatible(int n, BuiltFlag flag, tuple excluded, int max_signature):
    cdef list base_points = list(range(max_signature))
    cdef BuiltFlag fexclii
    cdef Pattern pexclii
    cdef tuple extra_points
    cdef int sii
    for ii in excluded:
        fexclii = ii
        sii = fexclii.size()
        if sii<max_signature:
            continue 
        for extra_points in itertools.combinations(range(max_signature, n), sii-max_signature):
            if flag.subflag(points=base_points+list(extra_points))==fexclii:
                return False
    return True

cpdef tuple inductive_generator(int n, theory, tuple smaller_structures, tuple excluded):
    
    cdef dict signature = theory.signature()
    cdef int max_arity = _get_max_arity(signature, n)
    #Handle the trivial case, when only the empty structure is possible
    cdef dict final_overlap
    if max_arity==0:
        final_overlap = {}
        for xx in signature:
            final_overlap[xx] = tuple()
        return tuple([BuiltFlag(theory, n, tuple(), **final_overlap)])

    #Handle the unary case
    cdef BuiltFlag F, final_flag
    cdef dict F_canonical_blocks

    cdef tuple perm_0 = tuple([n-1] + list(range(n-1)))
    cdef tuple extensions
    cdef dict ext
    cdef dict ret = {}
    cdef tuple patt

    if max_arity==1:
        extensions = _get_all_extensions(signature, n, 1)
        for F in smaller_structures:
            F_canonical_blocks = _perm_blocks(F._blocks, perm_0)
            for ext in extensions:
                sig_check()
                final_overlap = _merge_blocks(F_canonical_blocks, ext, tuple())
                final_flag = BuiltFlag(theory, n, tuple(), **final_overlap)
                if _excluded_compatible(n, final_flag, excluded, 1):
                    patt = final_flag._relation_list()
                    if patt not in ret:
                        ret[patt] = [final_flag]
                    elif final_flag not in ret[patt]:
                        ret[patt].append(final_flag)
        combined_tuple = tuple(itertools.chain.from_iterable(
            [ret[key] for key in sorted(ret.keys())]
        ))
        return combined_tuple
    
    extensions = _get_all_extensions(signature, n, 2)
    cdef list signature_perms = theory._signature_perms()

    #For the pre-calculation
    cdef dict subf_classes, key_lookup
    cdef int missing_point
    cdef list included_points
    cdef int ii
    cdef BuiltFlag Fs, Fs_canonical, F_canonical
    cdef tuple Fs_relabel, F_relabel
    cdef list Fs_cosets
    subf_classes = {}
    key_lookup = {}
    cdef set already_checked = set()
    cdef BuiltFlag pointed_flag
    
    cdef dict tad
    cdef list to_add
    cdef bint gd
    cdef dict F_can_perm_blocks
    cdef tuple sign_perm
    
    #Pre calculate the cosets and Fs classes
    #print("generate coset precalc")
    #for F in tqdm(smaller_structures):
    for F in smaller_structures:
        for missing_point in range(n-1):
            sig_check()
            pointed_flag = F.subflag(ftype_points=[missing_point])
            if pointed_flag in already_checked:
                continue
            already_checked.add(pointed_flag)

            included_points = [ii for ii in range(n-1) if ii!=missing_point]
            Fs = F.subflag(points=included_points)
            Fs_relabel = Fs.unique()[1]
            Fs_canonical = Fs.subflag(points=Fs_relabel)

            if Fs_canonical in key_lookup.keys():
                Fs_canonical = key_lookup[Fs_canonical]
            else:
                key_lookup[Fs_canonical] = Fs_canonical
                subf_classes[Fs_canonical] = []

            F_relabel = tuple([missing_point] + [included_points[ii] for ii in Fs_relabel])
            F_canonical = F.subflag(points=F_relabel)
            Fs_cosets = F_canonical._find_coset_representatives(Fs_canonical, (0, ))
            F_canonical_blocks = F_canonical._blocks

            if len(signature_perms)<=1:
                subf_classes[Fs_canonical].append((F_canonical_blocks, Fs_cosets))
            else:
                to_add = []
                gd = True
                for sign_perm in signature_perms:
                    if Fs._blocks==_perm_signature(Fs._blocks, sign_perm):
                        F_can_perm_blocks = _perm_signature(F_canonical_blocks, sign_perm)
                        for tad in to_add:
                            if tad==F_can_perm_blocks:
                                gd = False
                                break
                        if gd:
                            to_add.append(F_can_perm_blocks)
                subf_classes[Fs_canonical] += [(tad, Fs_cosets) for tad in to_add]

    #For the merging
    cdef dict G_canonical_blocks, FG_overlap, G_canonical_blocks_permed
    cdef tuple dummy_1
    cdef tuple G_perm
    cdef list G_perm_list

    #Check ways to combine
    #print("generate checking pairs")
    #for Fs_canonical in tqdm(subf_classes):
    for Fs_canonical in subf_classes:
        for ii, (F_canonical_blocks, _) in enumerate(subf_classes[Fs_canonical]):
            F_canonical_blocks = _perm_blocks(F_canonical_blocks, perm_0)
            for G_canonical_blocks, Gs_cosets in subf_classes[Fs_canonical][ii:]:
                for G_perm in Gs_cosets:
                    G_perm_list = list(G_perm)
                    G_perm_list.insert(1, n-1)
                    G_canonical_blocks_permed = _perm_blocks(G_canonical_blocks, tuple(G_perm_list))
                    FG_overlap = _merge_blocks(F_canonical_blocks, G_canonical_blocks_permed, (0, ))
                    for ext in extensions:
                        sig_check()
                        final_overlap = _merge_blocks(FG_overlap, ext, tuple())
                        final_flag = BuiltFlag(theory, n, tuple(), **final_overlap)
                        if _excluded_compatible(n, final_flag, excluded, 2):
                            patt = final_flag._relation_list()
                            if patt not in ret:
                                ret[patt] = [final_flag]
                            elif final_flag not in ret[patt]:
                                ret[patt].append(final_flag)
    
    combined_tuple = tuple(itertools.chain.from_iterable(
        [ret[key] for key in sorted(ret.keys())]
    ))
    return combined_tuple

cpdef tuple overlap_generator(int n, theoryR, tuple small0, tuple small1, tuple excluded):
    cdef list ret = []

    cdef int ii, jj
    cdef BuiltFlag fl0, fl1, final_flag

    cdef set allperm = set(itertools.permutations(range(n)))

    cdef dict fl0_reps
    cdef set aut0, aut1
    cdef tuple perm, patt
    cdef dict fl0blocks, fl1blocks, overlap, fl0_permed

    for ii, fl0 in enumerate(small0):
        aut0 = fl0.automorphisms()
        fl0blocks = fl0._blocks
        fl0_reps = _left_coset_representative_map(aut0, n)
        for jj, fl1 in enumerate(small1):
            aut1 = fl1.automorphisms()
            fl1blocks = fl1._blocks
            for perm in _find_double_coset_reps(aut0, aut1, fl0_reps, n):
                fl0_permed = _perm_blocks(fl0blocks, perm)
                overlap = fl0_permed | fl1blocks
                final_flag = BuiltFlag(theoryR, n, tuple(), **overlap)
                if _excluded_compatible(n, final_flag, excluded, 0):
                    # patt = (fl0, fl1)
                    # if patt not in ret:
                    #     ret[patt] = [final_flag]
                    # elif final_flag not in ret[patt]:
                    #     ret[patt].append(final_flag)
                    ret.append(final_flag)
    # combined_tuple = tuple(itertools.chain.from_iterable(
    #     [ret[key] for key in sorted(ret.keys())]
    # ))
    return tuple(ret)

cdef class _Flag(Element):
    cdef int _n
    cdef int _ftype_size
    
    cdef tuple _ftype_points
    cdef tuple _not_ftype_points
    cdef dict _blocks
    cdef tuple _unique

    def __init__(self, theory, n, ftype):
        self._n = int(n)
        self._ftype_points = tuple(ftype)
        self._ftype_size = len(self._ftype_points)
        self._not_ftype_points = None
        Element.__init__(self, theory)

    def _repr_(self):
        blocks = self._blocks
        blocks_reps = []
        for xx in blocks.keys():
            brx = xx + '=('
            bsb = []
            for ee in blocks[xx]:
                bsb.append("".join(map(str, ee)))
            brx += ' '.join(bsb)
            brx += ')'
            blocks_reps.append(brx)
        strblocks = ', '.join(blocks_reps)
        if self.is_ftype():
            return 'Ftype on {} points with {}'.format(self.size(), strblocks)
        return 'Flag on {} points, ftype from {} with {}'.format(self.size(), self.ftype_points(), strblocks)
    
    def _serialize(self):
        ret = [self.size(), self.ftype_points()]
        blocks = self._blocks
        for kk in blocks:
            ret.append(blocks[kk])
        return tuple(ret)

    def _pythonize(self):
        blret = {}
        blocks = self._blocks
        for kk in blocks:
            blret[kk] = tuple([
                tuple([int(pp) for pp in edge]) for edge in blocks[kk]
            ])
        return tuple([
            int(self.size()), 
            tuple([int(ii) for ii in self.ftype_points()]), 
            tuple(blret.items())
            ])

    # Basic properties

    def combinatorial_theory(self):
        r"""
        Returns the combinatorial theory this flag is a member of
        
        This is the same as the parent.

        .. SEEALSO::

            :func:`theory`
            :func:`parent`
        """
        return self.parent()
    
    theory = combinatorial_theory
    
    cpdef int size(self):
        r"""
        Returns the size of the vertex set of this Flag.

        OUTPUT: integer, the number of vertices.

        EXAMPLES::

        This is the size parameter in the `Flag` initialization ::

            sage: GraphTheory(4).size()
            4
        """
        return self._n
    
    vertex_number = size

    cpdef blocks(self, key=None):
        r"""
        Returns the blocks

        INPUT:

        - ``key`` -- 

        OUTPUT: 
        """
        if key==None:
            return self._blocks
        else:
            return self._blocks[key]
    
    cpdef tuple ftype_points(self):
        r"""
        The points of the ftype inside self.
        
        This gives an injection of ftype into self

        OUTPUT: list of integers

        EXAMPLES::

            
            sage: two_pointed_triangle = GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]], ftype=[0, 1])
            sage: two_pointed_triangle.ftype_points()
            (0, 1)

        .. SEEALSO::

            :func:`__init__`
        """
        return self._ftype_points

    cpdef tuple not_ftype_points(self):
        r"""
        This is a helper function, caches the points that are not
        part of the ftype.
        """
        if self._not_ftype_points != None:
            return self._not_ftype_points
        cdef int ii
        self._not_ftype_points = tuple([ii for ii in range(self.size()) if ii not in self._ftype_points])
        return self._not_ftype_points

    cpdef bint is_ftype(self):
        r"""
        Returns `True` if this flag is an ftype.

        .. SEEALSO::

            :func:`_repr_`
        """
        return self._n == self._ftype_size

    # Isomorphisms and related stuff

    cpdef tuple unique(self, bint weak = False):
        return ((), ())
    
    cpdef bint weak_equal(self, _Flag other):
        if self.theory() != other.theory():
            return False
        if self._ftype_size != other._ftype_size:
            return False
        cdef tuple sun = self.unique(weak=True)
        cdef tuple oun = other.unique(weak=True)
        return sun[0] == oun[0]
    
    cpdef bint normal_equal(self, _Flag other):
        if self.theory() != other.theory():
            return False
        if self._ftype_size != other._ftype_size:
            return False
        cdef tuple sun = self.unique(weak=False)
        cdef tuple oun = other.unique(weak=False)
        return sun[0] == oun[0]
    
    cpdef bint strong_equal(self, _Flag other):
        if self.theory() != other.theory():
            return False
        if self.ftype_points() != other.ftype_points():
            return False
        return self._blocks == other._blocks

    # Flag algebra compatibility
    def afae(self):
        from sage.algebras.flag_algebras import FlagAlgebra
        alg = FlagAlgebra(self.theory(), QQ, self.ftype())
        return alg(self)

    def custom_coerce(self, other):
        from sage.algebras.flag_algebras import FlagAlgebra, FlagAlgebraElement
        if isinstance(other, _Flag) or isinstance(other, Pattern):
            if self.ftype()!=other.ftype():
                raise ValueError("The ftypes must agree.")
            alg = FlagAlgebra(self.theory(), QQ, self.ftype())
            return (alg(self), alg(other))
        elif isinstance(other, FlagAlgebraElement):
            if self.ftype()!=other.ftype():
                raise ValueError("The ftypes must agree.")
            base = other.base()
            alg = other.parent()
            return (alg(self), other)
        else:
            base = coercion_model.common_parent(QQ, other.parent())
            alg = FlagAlgebra(self.theory(), base, self.ftype())
            return (alg(self), alg(base(other)))

    def _add_(self, other):
        sf, of = self.custom_coerce(other)
        return sf._add_(of)

    def _sub_(self, other):
        sf, of = self.custom_coerce(other)
        return sf._sub_(of)

    def _neg_(self):
        from sage.all import Integer
        return self.__mul__(Integer(-1))

    def __mul__(self, other):
        sf, of = self.custom_coerce(other)
        return sf._mul_(of)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        if other==0:
            return 1
        if other==1:
            return self
        else:
            ret = self
            for _ in range(other-1):
                ret *= self
            return ret

    def __lshift__(self, amount):
        r"""
        `FlagAlgebraElement`, equal to this, with size is shifted by the amount

        EXAMPLES::

        Edge shifted to size `3` ::

            sage: edge = GraphTheory(2, edges=[[0, 1]])
            sage: (edge>>1).values()
            (0, 1/3, 2/3, 1)

        .. SEEALSO::

            :func:`FlagAlgebraElement.__lshift__`
        """
        return self.afae().__lshift__(amount)
    
    def __rshift__(self, amount):
        r"""
        `FlagAlgebraElement`, averaged to an `amount` smaller size

        EXAMPLES::

        Cherry averaged to size `2` ::
            
            sage: cherry = GraphTheory(3, edges=[[0, 1], [0, 2]])
            sage: (cherry<<1).values()
            (0, 1/3, 2/3, 1)

        .. SEEALSO::

            :func:`FlagAlgebraElement.__rshift__`
        """
        return self.afae().__rshift__(amount)

    def __truediv__(self, other):
        r"""
        Divide by a scalar

        INPUT:

        - ``other`` -- number; any number such that `1` can be divided with that

        OUTPUT: The `FlagAlgebraElement` resulting from the division

        EXAMPLES::

        Divide by `2` ::

            
            sage: g = GraphTheory(3)
            sage: (g/2).values()
            (1/2, 0, 0, 0)
            
        Even for `x` symbolic `1/x` is defined, so the division is understood ::
            sage: var('x')
            x
            sage: g = GraphTheory(2)
            sage: g/x
            Flag Algebra Element over Symbolic Ring
            1/x - Flag on 2 points, ftype from () with edges=()
            0   - Flag on 2 points, ftype from () with edges=(01)
        
        .. NOTE::

            Dividing by `Flag` or `FlagAlgebraElement` is not allowed, only
            numbers such that the division is defined in some extension
            of the rationals.

        .. SEEALSO::

            :func:`FlagAlgebraElement.__truediv__`
        """
        return self.afae().__truediv__(other)
    
    def __eq__(self, other):
        r"""
        Compare two flags for == (equality)
        
        This is the isomorphism defined by the identifiers,
        respecting the types.

        .. SEEALSO::

            :func:`unique`
            :func:`theory`
            :func:`CombinatorialTheory.identify`
        """
        if type(other)!=type(self):
            return False
        if self.theory()!=other.theory():
            return False
        return self.normal_equal(other)
    
    def __lt__(self, other):
        r"""
        Compare two flags for < (proper induced inclusion)
        
        Returns true if self appears as a proper induced structure 
        inside other.

        .. SEEALSO::

            :func:`__le__`
        """
        if type(other)!=type(self):
            return False
        if self.theory()!=other.theory():
            return False
        if self.size()>=other.size():
            return False
        if self.ftype() != other.ftype():
            return False
        for subp in itertools.combinations(other.not_ftype_points(), self.size()-self.ftype().size()):
            sig_check()
            osub = other.subflag(subp)
            if self==osub:
                return True
        return False
    
    def __le__(self, other):
        r"""
        Compare two flags for <= (induced inclusion)
        
        Returns true if self appears as an induced structure inside
        other.

        EXAMPLES::

        Edge appears in a 4 star ::

            
            sage: star = GraphTheory(4, edges=[[0, 1], [0, 2], [0, 3]])
            sage: edge = GraphTheory(2, edges=[[0, 1]])
            sage: edge <= star
            True
            
        The ftypes must agree ::
        
            sage: p_edge = GraphTheory(2, edges=[[0, 1]], ftype_points=[0])
            sage: p_edge <= star
            False
        
        But when ftypes agree, the inclusion must respect it ::
            
            sage: pstar = star.subflag(ftype_points=[0])
            sage: sub1 = GraphTheory(3, ftype=[0], edges=[[0, 1], [0, 2]])
            sage: sub1 <= pstar
            True
            sage: sub2 = GraphTheory(3, ftype=[1], edges=[[0, 1], [0, 2]])
            sage: sub2 <= pstar
            False

        .. SEEALSO::

            :func:`__lt__`
            :func:`__eq__`
            :func:`unique`
        """
        return self==other or self<other
    
    def __hash__(self):
        r"""
        A hash based on the unique identifier
        so this is compatible with `__eq__`.
        """
        return hash(self.unique()[0])
    
    def __getstate__(self):
        r"""
        Saves this flag to a dictionary
        """
        dd = {'theory': self.theory(),
              'n': self._n, 
              'ftype_points': self._ftype_points, 
              'blocks':self._blocks, 
              'unique':self._unique}
        return dd
    
    def __setstate__(self, dd):
        r"""
        Loads this flag from a dictionary
        """
        self._n = dd['n']
        self._ftype_points = dd['ftype_points']
        self._ftype_size = len(self._ftype_points)
        self._not_ftype_points = None
        self._blocks = dd['blocks']
        self._unique = dd['unique']
    
    def project(self, ftype_inj=tuple()):
        r"""
        Project this `Flag` to a smaller ftype
        

        INPUT:

        - ``ftype_inj`` -- tuple (default: (, )); the injection of the
            projected ftype inside the larger ftype

        OUTPUT: the `FlagAlgebraElement` resulting from the projection

        EXAMPLES::

        If the center of a cherry is flagged, then the projection has
        coefficient 1/3 ::

            
            sage: p_cherry = GraphTheory(3, edges=[[0, 1], [0, 2]], ftype_points=[0])
            sage: p_cherry.project().values()
            (0, 0, 1/3, 0)

        .. NOTE::

            If `ftype_inj==tuple(range(self.ftype().size()))` then this
            does nothing.

        .. SEEALSO::

            :func:`FlagAlgebraElement.project`
        """
        return self.afae().project(ftype_inj)
    
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

            
            sage: p_edge = GraphTheory(2, edges=[[0, 1]], ftype_points=[0])
            sage: p_edge.mul_project(p_edge).values()
            (0, 0, 1/3, 1)

        .. NOTE::

            If `ftype_inj==tuple(range(self.ftype().size()))` then this
            is the same as usual multiplication.

        .. SEEALSO::

            :func:`_mul_`
            :func:`project`
            :func:`FlagAlgebraElement.mul_project`
        """
        return self.afae().mul_project(other, ftype_inj)
    
    def density(self, other):
        r"""
        The density of self in other.
        
        Randomly choosing self.size() points in other, the
        probability of getting self.

        EXAMPLES::

        Density of an edge in the cherry graph is 2/3 ::

            
            sage: cherry = GraphTheory(3, edges=[[0, 1], [0, 2]])
            sage: edge = GraphTheory(2, edges=[[0, 1]])
            sage: cherry.density(edge)
            2/3
        
        .. SEEALSO::
        
            :func:`FlagAlgebraElement.density`
        """
        safae = self.afae()
        oafae = safae.parent(other)
        return self.afae().density(oafae)

cdef class BuiltFlag(_Flag):
    cdef set _automorphisms
    cdef tuple _weak_unique
    cdef BuiltFlag _ftype
    
    def __init__(self, theory, n, ftype, **params):
        self._blocks = _standardize_blocks(params, theory._signature, False)
        self._automorphisms = None
        self._weak_unique = None
        self._unique = None
        self._ftype = None
        _Flag.__init__(self, theory, n, ftype)
    
    cpdef tuple _relation_list(self):
        cdef dict ret = {}
        cdef dict signature = self.theory().signature()
        cdef int group
        for xx in signature:
            group = signature[xx]["group"]
            if group in ret:
                ret[group].append(len(self._blocks[xx]))
            else:
                ret[group] = [len(self._blocks[xx])]
        return tuple([tuple(sorted(xx)) for xx in ret.values()])
    
    cpdef Pattern as_pattern(self):
        if self.ftype().size() != 0:
            import warnings
            warnings.warn("Transforming Flag to Pattern with nonempty ftype. Ftype will be ignored!")
        cdef dict pat_blocks = dict(self._blocks)
        cdef str xx
        cdef dict rel
        cdef tuple sx
        cdef list mx
        cdef int n = self.size()
        for xx in self.theory().signature():
            rel = self.theory().signature()[xx]
            sx = self._blocks[xx]
            mx = [tup for tup in _generate_all_relations(n, rel) if tup not in sx]
            pat_blocks[xx+"_m"] = tuple(mx)
        return Pattern(self.theory(), n, tuple(), **pat_blocks)

    # Basic properties
    
    cpdef dict signature(self):
        return self.theory().signature()

    cpdef BuiltFlag ftype(self):
        r"""
        Returns the ftype of this `Flag`

        EXAMPLES::

        Ftype of a pointed triangle is just a point ::

            
            sage: pointed_triangle = GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]], ftype=[0])
            sage: pointed_triangle.ftype()
            Ftype on 1 points with edges=()
        
        And with two points it is ::
        
            sage: two_pointed_triangle = GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]], ftype=[0, 1])
            sage: two_pointed_triangle.ftype()
            Ftype on 2 points with edges=(01)

        .. NOTE::

            This is essentially the subflag, but the order of points matter. The result is saved
            for speed.

        .. SEEALSO::

            :func:`subflag`
        """
        if self._ftype==None:
            if self.is_ftype():
                self._ftype = self
            self._ftype = self.subflag([])
        return self._ftype
    
    # Isomorphisms and related stuff

    cpdef tuple unique(self, bint weak = False):
        if weak and self._weak_unique != None:
            return self._weak_unique
        if (not weak) and self._unique != None:
            return self._unique

        cdef tuple symmetry_graph_data = self._symmetry_graph(weak)

        cdef tuple result = canonical_form_from_edge_list(
            symmetry_graph_data[0],
            symmetry_graph_data[1],
            symmetry_graph_data[2],
            symmetry_graph_data[3],
            symmetry_graph_data[4],
            symmetry_graph_data[5]
        )
        cdef list new_edges = sorted(result[0])
        cdef tuple uniret = tuple([len(self.ftype_points()), weak] + new_edges)
        cdef tuple relabel =tuple([result[1][i] for i in range(self.size())])
        relabel = tuple([relabel.index(ii) for ii in range(self.size())])
        cdef tuple ret = (uniret, relabel)
        if weak:
            self._weak_unique = ret
        else:
            self._unique = ret
        return ret
    
    cdef tuple _symmetry_graph(self, bint weak = False):
        #Should return:
        #-vertex number, 
        #-edges as two lists,
        #-edge label number
        #-edge labels
        #-partition

        cdef dict blocks = self._blocks
        cdef int next_vertex = self.size()
        cdef dict signature = self.signature()

        #Data for relations
        cdef dict rel_info
        cdef str rel_name
        cdef int group
        cdef bint ordered

        #
        #Creating lookup tables for the vertices/groups
        #
        cdef dict unary_relation_vertices = {}
        cdef dict tuple_vertices = {}
        cdef dict group_vertices = {}
        cdef dict groups = {}
        cdef list group_relation_vertices
        cdef tuple t

        cdef dict all_symmetry_vertices_remap = {}

        for rel_name in signature:
            rel_info = signature[rel_name]
            arity = rel_info['arity']
            group = rel_info['group']

            #Layer 1 vertices for the relations
            if arity == 1:
                unary_relation_vertices[rel_name] = next_vertex
                next_vertex += 1
            else:
                occurrences = blocks[rel_name]
                tuple_vertices[rel_name] = {}
                for t in occurrences:
                    tuple_vertices[rel_name][t] = next_vertex
                    next_vertex += 1
            
            #Layer 2 vertices
            group_vertices[rel_name] = next_vertex
            #Creating groups
            if group not in groups:
                #This is the first time we see this group number
                #Initiaize the group list with it's name
                groups[group] = [rel_name]
                all_symmetry_vertices_remap[group] = [next_vertex]
            else:
                #This group number was seen before so just update the lists
                groups[group].append(rel_name)
                all_symmetry_vertices_remap[group].append(next_vertex)
            next_vertex += 1
        

        #Adding the remaining layer 2 vertices for the symmetry graphs
        cdef tuple symmetries = self.theory()._symmetries
        cdef int n_sym, m_sym, ii
        cdef tuple edges_sym
        for group in groups:
            n_sym, m_sym, edges_sym = symmetries[group]
            all_symmetry_vertices_remap[group] += list(range(next_vertex, next_vertex+m_sym-n_sym))
            next_vertex += m_sym - n_sym

        #Creating the partition, first the vertices from layer 0
        cdef list partition

        if weak:
            partition = [self.not_ftype_points(), self.ftype_points()]
        else:
            partition = [self.not_ftype_points()]
            partition += [[ii] for ii in self.ftype_points()]

        #For groups also add the symmetry edges
        cdef list Vout = []
        cdef list Vin = []
        cdef tuple edge
        cdef list group_symmetry_remap
        for group in groups:
            #Layer 2 partition
            partition.append([group_vertices[rel_name] for rel_name in groups[group]])
            n_sym, m_sym, edges_sym = symmetries[group]
            group_symmetry_remap = all_symmetry_vertices_remap[group]
            partition.append(group_symmetry_remap[n_sym:])
            for edge in edges_sym:
                Vin.append(group_symmetry_remap[edge[0]])
                Vout.append(group_symmetry_remap[edge[1]])

            #Layer 1 partition
            group_relation_vertices = []
            for rel_name in groups[group]:
                if rel_name in unary_relation_vertices:
                    group_relation_vertices.append(unary_relation_vertices[rel_name])
                else:
                    group_relation_vertices += list(tuple_vertices[rel_name].values())
            partition.append(group_relation_vertices)

        # Build the edge lists for the rest
        cdef list labels = None
        cdef int edge_label_num = 1
        cdef list conns
        
        for rel_name in unary_relation_vertices:
            #This will find connections to unary_relation_vertices[rel_name] in Layer 1

            #From layer 0
            conns = [t[0] for t in blocks[rel_name]]
            #From layer 2
            conns.append(group_vertices[rel_name])
            Vout += conns
            Vin += [unary_relation_vertices[rel_name]] * len(conns)

        for rel_name in tuple_vertices:
            #Same but for every element in unary_relation_vertices[rel_name]
            rel_info = signature[rel_name]
            ordered = rel_info['ordered']
            arity = rel_info['arity']

            #If this is the first time we realize things must be ordered
            if ordered and arity!=1:
                if labels==None:
                    labels = [0]*len(Vin)
                edge_label_num = max(edge_label_num, arity)
            
            for block in tuple_vertices[rel_name]:
                conns = [group_vertices[rel_name]] + list(block)
                Vout += conns
                Vin += [tuple_vertices[rel_name][block]] * len(conns)
                if ordered and arity!=1:
                    labels.append(0)
                    labels += list(range(arity))
                elif labels!=None:
                    labels += [0] * len(conns)
        
        cdef tuple ret = (next_vertex, Vout, Vin, edge_label_num, labels, partition)
        return ret
    
    def sage_symmetry_graph(self):
        symmetry_graph_data = self._symmetry_graph()
        Vnr, Vout, Vin, Lnr, labels, partition = symmetry_graph_data
        vert_inds = []
        for xx in range(Vnr):
            for ii, par in enumerate(partition):
                if xx in par:
                    vert_inds.append(ii)
        vertices = list(zip(range(Vnr), vert_inds))
        Vout = [vertices[ii] for ii in Vout]
        Vin = [vertices[ii] for ii in Vin]
        if labels==None:
            labels = [0]*len(Vin)
        edges = list(zip(Vin, Vout, labels))
        from sage.graphs.graph import Graph
        return Graph([vertices, edges], format="vertices_and_edges")

    cpdef tuple automorphism_generators(self):
        if self._automorphisms!=None:
            return self._automorphisms
        cdef tuple symmetry_graph_data = self._symmetry_graph()
        cdef tuple result = automorphism_group_gens_from_edge_list(
            symmetry_graph_data[0],
            symmetry_graph_data[1],
            symmetry_graph_data[2],
            symmetry_graph_data[3],
            symmetry_graph_data[4],
            symmetry_graph_data[5]
        )
        return result

    cpdef set automorphisms(self):
        if self._automorphisms!=None:
            return self._automorphisms
        cdef tuple symmetry_graph_data = self._symmetry_graph()
        cdef tuple result = automorphism_group_gens_from_edge_list(
            symmetry_graph_data[0],
            symmetry_graph_data[1],
            symmetry_graph_data[2],
            symmetry_graph_data[3],
            symmetry_graph_data[4],
            symmetry_graph_data[5]
        )
        try:
            self._automorphisms = _generate_group(result, len(self.not_ftype_points()))
        except:
            print("The error occurred at ", self)
            raise RuntimeError("Stop here")
        return self._automorphisms

    cdef list _find_coset_representatives(self, BuiltFlag subflag, tuple points_missing):
        cdef set G_self, G_small, G_subf
        cdef list extended_perm
        cdef int n = self.size()
        cdef int k = subflag.size()
        cdef int ii

        G_self = self.automorphisms()
        G_small = subflag.automorphisms()
        G_subf = set()
        cdef list remap = [ii for ii in range(n) if ii not in points_missing]
        for perm in G_small:
            extended_perm = [remap[perm[ii]] for ii in range(k)]
            for ii in points_missing:
                extended_perm.insert(ii, ii)
            G_subf.add(tuple(extended_perm))
        return _compute_coset_reps(G_subf, G_self & G_subf, self.size())
    
    cpdef list nonequal_permutations(self):
        if len(self.ftype_points())==0:
            return self._typeless_nonequal_permutations()
        else:
            return self._typed_nonequal_permutations()
    
    cpdef list _typeless_nonequal_permutations(self):
        cdef list ssc = self.signature_changes()
        cdef set G_self = self.automorphisms()
        cdef set G_all = set(itertools.permutations(range(self.size())))
        cdef list cosreps = _compute_coset_reps(G_all, G_self, self.size())
        cdef BuiltFlag xx
        cdef tuple perm
        return [xx.subflag(points=perm) for perm in cosreps for xx in ssc]

    cpdef list _typed_nonequal_permutations(self):
        cdef list ssc = self.signature_changes()
        cdef list ret = []
        cdef BuiltFlag xx, yy, zz
        cdef bint gd
        for zz in ssc:
            for perm in itertools.permutations(range(self.size())):
                xx = zz.subflag(points=perm)
                gd = True
                for yy in ret:
                    if yy.strong_equal(xx):
                        gd = False
                        break
                if gd:
                    ret.append(xx)
        return ret

    cpdef list _generate_overlaps(self, other_theory, result_theory):
        if self.ftype().size() != 0:
            raise ValueError("Ftype must be empty")
        cdef list result = []
        cdef list result_partial = []
        cdef tuple other_flags = other_theory.generate(self.size())
        cdef BuiltFlag other, other_permed, overlap
        
        for other in other_flags:
            result_partial = []
            for other_permed in other.nonequal_permutations():
                overlap = BuiltFlag(result_theory, self.size(), **self._blocks, **other_permed._blocks)
                if overlap not in result_partial:
                    result_partial.append(overlap)
            result += result_partial
        return result

    
    # Core loops

    cpdef BuiltFlag subflag(self, points=None, ftype_points=None):
        r"""
        Returns the induced subflag.
        """
        cdef int ii
        if ftype_points==None:
            ftype_points = self._ftype_points
        if points==None:
            points = tuple(range(self._n))
        else:
            points = tuple(list(points) + [ii for ii in ftype_points if ii not in points])
        if len(points)==0:
            return self.theory().empty()
        ftype_points = tuple(ftype_points)
        if points==tuple(range(self.size())) and ftype_points==self._ftype_points:
            return self
        blocks = {xx: _subblock_helper(points, self._blocks[xx]) for xx in self._blocks.keys()}
        new_ftype_points = [points.index(ii) for ii in ftype_points]
        return BuiltFlag(self.theory(), len(points), new_ftype_points, **blocks)

    cpdef list ftypes_inside(self, target):
        r"""
        Returns the possible ways self ftype appears in target

        INPUT:

        - ``target`` -- Flag; the flag where we are looking for copies of self

        OUTPUT: list of Flags with ftype matching as self, not necessarily unique
        """
        cdef list ret = []
        cdef list lrp = list(range(target.size()))
        cdef tuple ftype_points
        cdef BuiltFlag newflag
        for ftype_points in itertools.permutations(range(target.size()), self._n):
            sig_check()
            if target.subflag(ftype_points, ftype_points)==self:
                newflag = target.subflag(lrp, ftype_points)
                if newflag not in ret:
                    ret.append(newflag)
        return ret

    cpdef list signature_changes(self):
        cdef list ret = []
        cdef list perms = self.theory()._signature_perms()
        cdef tuple perm
        cdef dict nblocks
        cdef str xx
        cdef int ii
        if len(perms)==1:
            return [self]
        for perm in perms:
            nblocks = {}
            for ii, xx in enumerate(self._blocks.keys()):
                nblocks[perm[ii]] = self._blocks[xx]
            ret.append(BuiltFlag(self.theory(), self.size(), self._ftype_points, **nblocks))
        return ret
    
    cpdef densities(self, int n1, tuple n1flgs, int n2, tuple n2flgs, \
    list ftype_remap, BuiltFlag large_ftype, BuiltFlag small_ftype):
        r"""
        Returns the density matrix, indexed by the entries of `n1flgs` and `n2flgs`
        
        The matrix returned has entry `(i, j)` corresponding to the possibilities of
        `n1flgs[i]` and `n2flgs[j]` inside self, projected to the small ftype. 
        
        This is the same as counting the ways we can choose `n1` and `n2` points
        inside `self`, such that the two sets cover the entire `self.size()` point
        set, and calculating the probability that the overlap induces an ftype
        isomorphic to `large_ftype` and the points sets are isomorphic to 
        `n1flgs[i]` and `n2flgs[j]`.

        INPUT:

        - ``n1`` -- integer; the size of the first flag list
        - ``n1flgs`` -- tuple of flags; the first flag tuple (each of size `n1`)
        - ``n2`` -- integer; the size of the second flag list
        - ``n2flgs`` -- tuple of flags; the second flag tuple (each of size `n2`)
        - ``ftype_remap`` -- list; shows how to remap `small_ftype` into `large_ftype`
        - ``large_ftype`` -- ftype; the ftype of the overlap
        - ``small_ftype`` -- ftype; the ftype of self

        OUTPUT: a sparse matrix corresponding with the counts

        .. SEEALSO::

            :func:`CombinatorialTheory.mul_project_table`
            :func:`FlagAlgebra.mul_project_table`
            :func:`FlagAlgebraElement.mul_project`
        """
        cdef int N = self._n
        cdef int small_size = small_ftype.size()
        cdef int large_size = large_ftype.size()
        cdef int ctr = 0
        
        cdef dict ret = {}
        cdef tuple small_points = self._ftype_points
        cdef int ii, vii, n1_ind, n2_ind
        cdef tuple difference
        cdef list large_points, remaining_points, not_large_points
        cdef BuiltFlag n1_subf, n2_subf, ind_large_ftype

        for difference in itertools.permutations(self.not_ftype_points(), large_size - small_size):
            sig_check()
            large_points = [0]*len(ftype_remap)
            for ii in range(len(ftype_remap)):
                vii = ftype_remap[ii]
                if vii<small_size:
                    large_points[ii] = small_points[vii]
                else:
                    large_points[ii] = difference[vii-small_size]
            ind_large_ftype = self.subflag([], ftype_points=large_points)
            if ind_large_ftype==large_ftype:
                not_large_points = [ii for ii in range(N) if ii not in large_points]
                for n1_extra_points in itertools.combinations(not_large_points, n1 - large_size):
                    n1_subf = self.subflag(n1_extra_points, ftype_points=large_points)
                    try:
                        n1_ind = n1flgs.index(n1_subf)
                    except ValueError:
                        raise ValueError("Could not find \n", n1_subf, "\nin the list of ", \
                                         n1, " sized flags with ", large_ftype, \
                                         ".\nThis can happen if the generator and identifier ",\
                                         "(from the current CombinatorialTheory) is incompatible, ",\
                                         "or if the theory is not heredetary")
                    
                    remaining_points = [ii for ii in not_large_points if ii not in n1_extra_points]
                    for n2_extra_points in itertools.combinations(remaining_points, n2 - large_size):
                        n2_subf = self.subflag(n2_extra_points, ftype_points=large_points)
                        try:
                            n2_ind = n2flgs.index(n2_subf)
                        except:
                            raise ValueError("Could not find \n", n2_subf, "\nin the list of ", \
                                             n2, " sized flags with ", large_ftype, \
                                             ".\nThis can happen if the generator and identifier ",\
                                             "(from the current CombinatorialTheory) is incompatible, ",\
                                             "or if the theory is not heredetary")
                        try:
                            ret[(n1_ind, n2_ind)] += 1
                        except:
                            ret[(n1_ind, n2_ind)] = 1
        return (len(n1flgs), len(n2flgs), ret)

    # Flag algebra compatibility

    def __getstate__(self):
        r"""
        Saves this flag to a dictionary
        """
        dd = _Flag.__getstate__(self)
        dd['weak_unique'] = self._weak_unique
        dd['automorphisms'] = self._automorphisms
        return dd
    
    def __setstate__(self, dd):
        r"""
        Loads this flag from a dictionary
        """
        _Flag.__setstate__(self, dd)
        self._set_parent(dd['theory'])
        self._weak_unique = dd['weak_unique']
        self._automorphisms = dd['automorphisms']
        self._ftype = None

cdef class ExoticFlag(_Flag):

    cdef ExoticFlag _ftype
    
    def __init__(self, theory, n, ftype, **params):
        r"""
        Initialize a :class:`Flag` element

        INPUT:

        - ``theory`` -- :class:`CombinatorialTheory`; the underlying theory, 
            which is also the parent
        - ``n`` -- integer; the size (number of vertices) of the flag
        - ``**params`` -- optional parameters; can contain the points of the
            ftype as a list of points under the name "ftype" of "ftype_points", 
            if not provided then empty ftype is assumed. Can also contain the 
            blocks for each element in the signature, if not provided then an 
            empty block set is assumed.

        OUTPUT: The resulting flag

        EXAMPLES::

        Create a simple GraphTheory triangle ::
            
            sage: from sage.algebras.flag_algebras import *
            sage: from sage.algebras.flag import Flag
            sage: Flag(GraphTheory, 3, edges=[[0, 1], [0, 2], [1, 2]])
            Flag on 3 points, ftype from [] with edges=[[0, 1], [0, 2], [1, 2]]
        
        Create a DiGraphTheory edge out from a pointed ftype ::

            sage: Flag(DiGraphTheory, 2, edges=[[0, 1]], ftype=[0])
            Flag on 2 points, ftype from [0] with edges=[[0, 1]]
        
        .. NOTE::

            It is recommended to create flags using the parent's element constructor

        .. SEEALSO::

            :func:`CombinatorialTheory._element_constructor_`
        """
        self._blocks = {}
        for xx in theory._signature.keys():
            if xx in params:
                self._blocks[xx] = tuple([tuple(yy) for yy in params[xx]])
            else:
                self._blocks[xx] = tuple()
        self._unique = ()
        self._ftype = None
        _Flag.__init__(self, theory, n, ftype)
    
    # Basic properties

    cpdef ExoticFlag ftype(self):
        r"""
        Returns the ftype of this `Flag`

        EXAMPLES::

        Ftype of a pointed triangle is just a point ::

            
            sage: pointed_triangle = GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]], ftype=[0])
            sage: pointed_triangle.ftype()
            Ftype on 1 points with edges=()
        
        And with two points it is ::
        
            sage: two_pointed_triangle = GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]], ftype=[0, 1])
            sage: two_pointed_triangle.ftype()
            Ftype on 2 points with edges=(01)

        .. NOTE::

            This is essentially the subflag, but the order of points matter. The result is saved
            for speed.

        .. SEEALSO::

            :func:`subflag`
        """
        if self._ftype==None:
            if self.is_ftype():
                self._ftype = self
            self._ftype = self.subflag([])
        return self._ftype
    

    # Isomorphisms and related stuff

    cpdef tuple unique(self, bint weak = False):
        r"""
        This returns a unique identifier that can equate isomorphic
        objects

        EXAMPLES::

        Isomorphic graphs have the same :func:`unique` value ::

            sage: from sage.algebras.flag_algebras import *
            sage: b1 = [[0, 1], [0, 2], [0, 4], [1, 3], [2, 4]]
            sage: b2 = [[0, 4], [1, 2], [1, 3], [2, 3], [3, 4]]
            sage: g1 = GraphTheory(5, edges=b1)
            sage: g2 = GraphTheory(5, edges=b2)
            sage: g1.unique() == g2.unique()
            True
            
        .. NOTE::

            The value returned here depends on the values of
            the parent `CombinatorialTheory`

        .. SEEALSO::

            :func:`theory`
            :func:`CombinatorialTheory.identify`
            :func:`__eq__`
        """
        if weak:
            return (self.theory().identify(self._n, [self._ftype_points], 
            **self._blocks), None)
        if self._unique==():
            self._unique = self.theory().identify(
                self._n, self._ftype_points, **self._blocks)
        return (self._unique, None)
    
    # Core loops

    cpdef ExoticFlag subflag(self, points=None, ftype_points=None):
        r"""
        Returns the induced subflag.
        
        The resulting sublaf contains the union of points and ftype_points
        and has ftype constructed from ftype_points. 

        INPUT:

        - ``points`` -- list (default: `None`); the points inducing the subflag.
            If not provided (or `None`) then this is the entire vertex set, so
            only the ftype changes
        - ``ftype_points`` list (default: `None`); the points inducing the ftype
            of the subflag. If not provided (or `None`) then the original ftype
            point set is used, so the result of the ftype will be the same

        OUTPUT: The induced sub Flag

        EXAMPLES::

        Same ftype ::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3, edges=[[0, 1]], ftype=[0])
            sage: g.subflag([0, 2])
            Flag on 2 points, ftype from [0] with edges=[]
            
        Only change ftype ::
            
            sage: g.subflag(ftype_points=[0, 1])
            Flag on 3 points, ftype from [0, 1] with edges=[[0, 1]]

        .. NOTE::

            As the ftype points can be chosen, the result can have different
            ftype as self.

        TESTS::

            sage: g.subflag()==g
            True
        """
        if ftype_points==None:
            ftype_points = self._ftype_points
        if points==None:
            points = tuple(range(self._n))
        else:
            points = tuple([ii for ii in range(self._n) if (ii in points or ii in ftype_points)])
        if len(points)==self._n and ftype_points==self._ftype_points:
            return self
        blocks = {xx: _subblock_helper(points, self._blocks[xx]) for xx in self._blocks.keys()}
        new_ftype_points = [points.index(ii) for ii in ftype_points]
        return ExoticFlag(self.parent(), len(points), ftype=new_ftype_points, **blocks)

    cpdef list ftypes_inside(self, target):
        r"""
        Returns the possible ways self ftype appears in target

        INPUT:

        - ``target`` -- Flag; the flag where we are looking for copies of self

        OUTPUT: list of Flags with ftype matching as self, not necessarily unique
        """
        cdef list ret = []
        cdef list lrp = list(range(target.size()))
        cdef tuple ftype_points
        cdef ExoticFlag newflag
        for ftype_points in itertools.permutations(range(target.size()), self._n):
            sig_check()
            if target.subflag(ftype_points, ftype_points)==self:
                newflag = target.subflag(lrp, ftype_points)
                if newflag not in ret:
                    ret.append(newflag)
        return ret

    cpdef densities(self, n1, n1flgs, n2, n2flgs, ftype_remap, large_ftype, small_ftype):
        r"""
        Returns the density matrix, indexed by the entries of `n1flgs` and `n2flgs`
        
        The matrix returned has entry `(i, j)` corresponding to the possibilities of
        `n1flgs[i]` and `n2flgs[j]` inside self, projected to the small ftype. 
        
        This is the same as counting the ways we can choose `n1` and `n2` points
        inside `self`, such that the two sets cover the entire `self.size()` point
        set, and calculating the probability that the overlap induces an ftype
        isomorphic to `large_ftype` and the points sets are isomorphic to 
        `n1flgs[i]` and `n2flgs[j]`.

        INPUT:

        - ``n1`` -- integer; the size of the first flag list
        - ``n1flgs`` -- list of flags; the first flag list (each of size `n1`)
        - ``n2`` -- integer; the size of the second flag list
        - ``n2flgs`` -- list of flags; the second flag list (each of size `n2`)
        - ``ftype_remap`` -- list; shows how to remap `small_ftype` into `large_ftype`
        - ``large_ftype`` -- ftype; the ftype of the overlap
        - ``small_ftype`` -- ftype; the ftype of self

        OUTPUT: a sparse matrix corresponding with the counts

        .. SEEALSO::

            :func:`CombinatorialTheory.mul_project_table`
            :func:`FlagAlgebra.mul_project_table`
            :func:`FlagAlgebraElement.mul_project`
        """
        cdef int N = self._n
        cdef int small_size = small_ftype.size()
        cdef int large_size = large_ftype.size()
        cdef int ctr = 0
        cdef bint chk = 0
        cdef int valid_ftypes = 0
        cdef int correct_ftypes = 0
        cdef int valid_flag_pairs = 0
        
        ret = {}
        small_points = self._ftype_points
        for difference in itertools.permutations(self.not_ftype_points(), large_size - small_size):
            sig_check()
            large_points = [0]*len(ftype_remap)
            for ii in range(len(ftype_remap)):
                vii = ftype_remap[ii]
                if vii<small_size:
                    large_points[ii] = small_points[vii]
                else:
                    large_points[ii] = difference[vii-small_size]
            ind_large_ftype = self.subflag([], ftype_points=large_points)
            if ind_large_ftype.unique()[0]==None:
                continue
            
            valid_ftypes += 1
            if ind_large_ftype==large_ftype:
                correct_ftypes += 1
                not_large_points = [ii for ii in range(N) if ii not in large_points]
                for n1_extra_points in itertools.combinations(not_large_points, n1 - large_size):
                    n1_subf = self.subflag(n1_extra_points, ftype_points=large_points)
                    if n1_subf.unique()[0]==None:
                        continue
                    try:
                        n1_ind = n1flgs.index(n1_subf)
                    except ValueError:
                        raise ValueError("Could not find \n", n1_subf, "\nin the list of ", \
                                         n1, " sized flags with ", large_ftype, \
                                         ".\nThis can happen if the generator and identifier ",\
                                         "(from the current CombinatorialTheory) is incompatible, ",\
                                         "or if the theory is not heredetary")
                    
                    remaining_points = [ii for ii in not_large_points if ii not in n1_extra_points]
                    for n2_extra_points in itertools.combinations(remaining_points, n2 - large_size):
                        n2_subf = self.subflag(n2_extra_points, ftype_points=large_points)
                        if n2_subf.unique()[0]==None:
                            continue
                        valid_flag_pairs += 1
                        try:
                            n2_ind = n2flgs.index(n2_subf)
                        except:
                            raise ValueError("Could not find \n", n2_subf, "\nin the list of ", \
                                             n2, " sized flags with ", large_ftype, \
                                             ".\nThis can happen if the generator and identifier ",\
                                             "(from the current CombinatorialTheory) is incompatible, ",\
                                             "or if the theory is not heredetary")
                        try:
                            ret[(n1_ind, n2_ind)] += 1
                        except:
                            ret[(n1_ind, n2_ind)] = 1
        return (len(n1flgs), len(n2flgs), ret, valid_ftypes, valid_flag_pairs, correct_ftypes)

    # Flag algebra compatibility
    
    def __setstate__(self, dd):
        r"""
        Loads this flag from a dictionary
        """
        _Flag.__setstate__(self, dd)
        self._set_parent(dd['theory'])
        self._ftype = None

cdef tuple _generate_all_relations(int n, dict relation):
    cdef int arity = relation["arity"]
    cdef bint is_ordered = relation["ordered"]
    if is_ordered:
        return tuple(itertools.permutations(range(n), r=arity))
    return tuple(itertools.combinations(range(n), arity))

cdef bint _block_refinement(tuple block_flag, tuple block_pattern, tuple missing_pattern):
    cdef tuple xx
    for xx in block_pattern:
        if xx not in block_flag:
            return False
    for xx in missing_pattern:
        if xx in block_flag:
            return False
    return True

cdef dict _pattern_overlap_fix(tuple ftype_points, dict blocks, dict signature):
    cdef str xx
    cdef tuple blxx, msxx, tp
    cdef dict fix = {}
    for xx in signature:
        blxx = blocks[xx]
        msxx = blocks[xx+"_m"]
        #Check they are disjoint
        if fix=={}:
            for tp in blxx:
                if tp in msxx:
                    fix = blocks
                    fix[xx+"_m"] = tuple([tp for tp in msxx if tp not in blxx])
                    break
        else:
            fix[xx+"_m"] = tuple([tp for tp in msxx if tp not in blxx])
    return fix

cdef dict _pattern_ftype_fix(tuple ftype_points, dict blocks, dict signature):
    cdef str xx
    cdef tuple blxx, msxx, tp, alrel, fttp
    cdef dict fix = {}
    for xx in signature:
        blxx = blocks[xx]
        msxx = blocks[xx+"_m"]
        #Check ftype has no undefined edges
        alrel = _generate_all_relations(len(ftype_points), signature[xx])
        if fix=={}:
            for tp in alrel:
                fttp = tuple([ftype_points[ii] for ii in tp])
                if (fttp not in blxx) and (fttp not in msxx):
                    fix = blocks
                    newmsxx = list(blocks[xx+"_m"])
                    for tp in alrel:
                        fttp = tuple([ftype_points[ii] for ii in tp])
                        if (fttp not in blxx) and (fttp not in msxx):
                            newmsxx.append(fttp)
                    fix[xx+"_m"] = tuple(newmsxx)
                    break
        else:
            newmsxx = list(blocks[xx+"_m"])
            for tp in alrel:
                fttp = tuple([ftype_points[ii] for ii in tp])
                if (fttp not in blxx) and (fttp not in msxx):
                    newmsxx.append(fttp)
            fix[xx+"_m"] = tuple(newmsxx)
    return fix

cdef class Pattern(_Flag):
    
    cdef BuiltFlag _ftype
    
    def __init__(self, theory, n, ftype, **params):
        _Flag.__init__(self, theory, n, ftype)
        cdef dict blocks = _standardize_blocks(params, theory._signature, True)
        cdef dict blfix = _pattern_overlap_fix(self._ftype_points, blocks, theory._signature)
        if blfix!={}:
            import warnings
            warnings.warn("""The pattern is initialized with required relations and 
            missing relations overlapping. The required relations will take priority.""")
            blocks = _standardize_blocks(blfix, theory._signature, True)
        blfix = _pattern_ftype_fix(self._ftype_points, blocks, theory._signature)
        if blfix!={}:
            import warnings
            warnings.warn("""The pattern is initialized with optional relations inside 
            the ftype. Those relations will be missing in the resulting pattern.""")
            blocks = _standardize_blocks(blfix, theory._signature, True)
        self._blocks = blocks
        self._ftype = None
    
    def _repr_(self):
        blocks = self._blocks
        blocks_reps = []
        for xx in blocks.keys():
            brx = xx + '=('
            bsb = []
            for ee in blocks[xx]:
                bsb.append("".join(map(str, ee)))
            brx += ' '.join(bsb)
            brx += ')'
            blocks_reps.append(brx)
        strblocks = ', '.join(blocks_reps)
        return 'Pattern on {} points, ftype from {} with {}'.format(self.size(), self.ftype_points(), strblocks)
    
    __str__ = _repr_
    
    def _serialize(self):
        ret = ["pattern", self.size(), self.ftype_points()]
        blocks = self._blocks
        for kk in blocks:
            ret.append(blocks[kk])
        return tuple(ret)
    
    #Basic properties

    cpdef BuiltFlag ftype(self):
        if self._ftype==None:
            blocks = {xx: _subblock_helper(self._ftype_points, self._blocks[xx]) for xx in self.theory().signature().keys()}
            self._ftype = BuiltFlag(self.theory(), len(self._ftype_points), self._ftype_points, **blocks)
        return self._ftype
    
    cpdef Pattern as_pattern(self):
        return self

    #Pattern methods

    cpdef bint is_compatible(self, BuiltFlag other):
        if self.size() != other.size():
            return False
        if self.theory() != other.theory():
            return False
        if self.ftype() != other.ftype():
            return False
        
        cdef dict signature = self.theory().signature()

        cdef dict oblocks = _perm_blocks(
            other._blocks, 
            tuple(itertools.chain(other.not_ftype_points(), other.ftype_points()))
        )
        oblocks = _standardize_blocks(oblocks, signature, False)
        
        cdef dict sblocks = _perm_blocks(
            self._blocks,
            tuple(itertools.chain(self.not_ftype_points(), self.ftype_points()))
        )

        cdef tuple rems = tuple(range(len(self.not_ftype_points()), self.size()))
        cdef dict tsblocks, chsblocks

        cdef list signperms = self.theory()._signature_perms()

        for notftype_perm in itertools.permutations(range(len(self.not_ftype_points()))):
            tsblocks = _perm_blocks(sblocks, tuple(itertools.chain(notftype_perm, rems)))
            for sperm in signperms:
                chsblocks = _perm_pattern_signature(tsblocks, sperm, self.theory().signature())
                chsblocks = _standardize_blocks(chsblocks, signature, True)
                if all([
                _block_refinement(oblocks[xx], chsblocks[xx], chsblocks[xx+"_m"]) 
                for xx in self.theory().signature().keys()
                ]):
                    return True
        return False

    cpdef subpattern(self, points=None, ftype_points=None):
        if ftype_points==None:
            ftype_points = self._ftype_points
        
        if points==None:
            points = list(ftype_points) + [ii for ii in range(self._n) if ii not in ftype_points]
        else:
            points = list(ftype_points) + [ii for ii in points if ii not in ftype_points]
        if set(ftype_points)!=set(self._ftype_points):
            raise ValueError("Subflag for patterns is not defined with different ftype!")
        blocks = {xx: _subblock_helper(points, self._blocks[xx]) for xx in self._blocks.keys()}
        new_ftype_points = [points.index(ii) for ii in ftype_points]
        return Pattern(self.theory(), len(points), new_ftype_points, **blocks)

    #Compatibility with flag algebras

    cpdef list compatible_flags(self):
        aflags = self.theory().generate(self.size(), self.ftype())
        ret = [xx for xx in aflags if self.is_compatible(xx)]
        return ret

    def as_flag_algebra_element(self, base=QQ):
        from sage.algebras.flag_algebras import FlagAlgebra
        from sage.all import vector

        targ_alg = FlagAlgebra(self.theory(), base=base, ftype=self.ftype())
        aflags = self.theory().generate(self.size(), self.ftype())
        targ_vec = {}
        for ii,xx in enumerate(aflags):
            if self.is_compatible(xx):
                targ_vec[ii]=1
        return targ_alg(self.size(), vector(base, len(aflags), targ_vec))
    
    afae = as_flag_algebra_element
