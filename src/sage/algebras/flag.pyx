r"""
Implementation of Flag, elements of :class:`CombinatorialTheory`

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

import itertools
from sage.rings.rational_field import QQ
from cysignals.signals cimport sig_check
from sage.structure.element cimport Element
from blisspy cimport canonical_form_from_edge_list, automorphism_group_gens_from_edge_list

# Elementary block operations
cdef tuple _subblock_helper(tuple points, tuple block):
    if len(block)==0:
        return tuple()
    cdef set points_set = set(points)
    cdef dict points_index = {p: i for i, p in enumerate(points)}
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

cdef dict _perm_blocks(dict blocks, tuple perm):
    cdef dict ret
    cdef str xx
    ret = {
        xx: _subblock_helper(perm, blocks[xx]) 
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


# Automorphism operations

cdef set _generate_group(tuple generators, int n):
    cdef set group = set()
    cdef list to_check = [tuple(range(n))]
    cdef tuple perm, gen, new_perm
    cdef int i
    try:
        while to_check:
            perm = to_check.pop()
            if perm in group:
                continue
            group.add(perm)
            for gen in generators:
                new_perm = tuple(gen[perm[i]] for i in range(n))
                if new_perm not in group:
                    to_check.append(new_perm)
    except:
        print("error occurred with running this on ", generators, " n=", n)
        raise RuntimeError("Stop here")
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

# Helpers for the inductive generator

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

cdef bint _excluded_compatible(int n, Flag flag, tuple excluded, int max_signature):
    cdef list base_points = list(range(max_signature))
    cdef Flag fexclii
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

cpdef tuple inductive_generator(int n, theory, tuple smaller_structures, dict signature, tuple excluded):
    cdef int max_arity = _get_max_arity(signature, n)
    #Handle the trivial case, when only the empty structure is possible
    cdef dict final_overlap
    if max_arity==0:
        final_overlap = {}
        for xx in signature:
            final_overlap[xx] = tuple()
        return tuple([Flag(theory, n, tuple(), **final_overlap)])
    


    #Handle the unary case
    cdef Flag F, final_flag
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
                final_flag = Flag(theory, n, tuple(), **final_overlap)
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

    #For the pre-calculation
    cdef dict subf_classes, key_lookup
    cdef int missing_point
    cdef list included_points
    cdef int ii
    cdef Flag Fs, Fs_canonical, F_canonical
    cdef tuple Fs_relabel, F_relabel
    cdef list Fs_cosets
    subf_classes = {}
    key_lookup = {}
    cdef set already_checked = set()
    cdef Flag pointed_flag
    
    #Pre calculate the cosets and Fs classes
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
            subf_classes[Fs_canonical].append((F_canonical._blocks, Fs_cosets))

    #For the merging
    cdef dict G_canonical_blocks, FG_overlap, G_canonical_blocks_permed
    cdef tuple dummy_1
    cdef tuple G_perm
    cdef list G_perm_list

    #Check ways to combine
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
                        final_flag = Flag(theory, n, tuple(), **final_overlap)
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

cdef class Flag(Element):
    
    cdef int _n
    cdef int _ftype_size
    
    cdef tuple _ftype_points
    cdef tuple _not_ftype_points
    cdef dict _blocks
    cdef tuple _unique
    cdef set _automorphisms
    cdef tuple _weak_unique
    
    cdef Flag _ftype
    
    def __init__(self, theory, n, ftype, **params):
        self._n = int(n)
        self._ftype_points = tuple(ftype)
        self._ftype_size = len(self._ftype_points)
        self._not_ftype_points = None
        self._blocks = _standardize_blocks(params, theory._signature, False)
        self._unique = None
        self._automorphisms = None
        self._weak_unique = None
        self._ftype = None
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

    cpdef tuple _relation_list(self):
        cdef dict ret = {}
        cdef dict signature = self.parent().signature()
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
        for xx in self.parent().signature():
            rel = self.parent().signature()[xx]
            sx = self._blocks[xx]
            mx = [tup for tup in _generate_all_relations(n, rel) if tup not in sx]
            pat_blocks[xx+"_m"] = tuple(mx)
        return Pattern(self.theory(), n, tuple(), **pat_blocks)

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

            sage: from sage.algebras.flag_algebras import *
            sage: GraphTheory(4).size()
            4
        """
        return self._n
    
    vertex_number = size

    cpdef dict signature(self):
        return self.theory().signature()

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

            sage: from sage.algebras.flag_algebras import *
            sage: two_pointed_triangle = GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]], ftype=[0, 1])
            sage: two_pointed_triangle.ftype_points()
            [0, 1]

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

    cpdef Flag ftype(self):
        r"""
        Returns the ftype of this `Flag`

        EXAMPLES::

        Ftype of a pointed triangle is just a point ::

            sage: from sage.algebras.flag_algebras import *
            sage: pointed_triangle = GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]], ftype=[0])
            sage: pointed_triangle.ftype()
            Ftype on 1 points with edges=[]
        
        And with two points it is ::
        
            sage: two_pointed_triangle = GraphTheory(3, edges=[[0, 1], [0, 2], [1, 2]], ftype=[0, 1])
            sage: two_pointed_triangle.ftype()
            Ftype on 2 points with edges=[[0, 1]]

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
    
    cpdef bint is_ftype(self):
        r"""
        Returns `True` if this flag is an ftype.

        .. SEEALSO::

            :func:`_repr_`
        """
        return self._n == self._ftype_size


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
        cdef tuple symmetries = self.parent()._symmetries

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
            
            #Creating groups
            if group not in groups:
                groups[group] = [rel_name]
            else:
                groups[group].append(rel_name)
            #Layer 2 vertices
            group_vertices[rel_name] = next_vertex
            next_vertex += 1

        #Creating the partition, first the vertices from layer 0
        cdef list partition

        if weak:
            partition = [self.ftype_points(), self.not_ftype_points()]
        else:
            partition = [[ii] for ii in self.ftype_points()]
            partition.append(self.not_ftype_points())

        for group in groups:
            #Layer 2 partition
            partition.append([group_vertices[rel_name] for rel_name in groups[group]])

            #Layer 1 partition
            group_relation_vertices = []
            for rel_name in groups[group]:
                if rel_name in unary_relation_vertices:
                    group_relation_vertices.append(unary_relation_vertices[rel_name])
                else:
                    group_relation_vertices += list(tuple_vertices[rel_name].values())
            partition.append(group_relation_vertices)

        # Build the edge lists
        cdef list Vout = []
        cdef list Vin = []
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
        #print(ret)
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

    cpdef bint weak_equal(self, Flag other):
        if self.theory() != other.theory():
            return False
        if self._ftype_size != other._ftype_size:
            return False
        cdef tuple sun = self.unique(weak=True)
        cdef tuple oun = other.unique(weak=True)
        return sun[0] == oun[0]
    
    cpdef bint normal_equal(self, Flag other):
        if self.theory() != other.theory():
            return False
        if self._ftype_size != other._ftype_size:
            return False
        cdef tuple sun = self.unique(weak=False)
        cdef tuple oun = other.unique(weak=False)
        return sun[0] == oun[0]
    
    cpdef bint strong_equal(self, Flag other):
        if self.theory() != other.theory():
            return False
        if self.ftype_points() != other.ftype_points():
            return False
        return self._blocks == other._blocks
    
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
            self._automorphisms = _generate_group(result, self.size())
        except:
            print("The error occurred at ", self)
            raise RuntimeError("Stop here")
        return self._automorphisms

    cdef list _find_coset_representatives(self, Flag subflag, tuple points_missing):
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


    # Core loops

    cpdef Flag subflag(self, points=None, ftype_points=None):
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
            return self.parent().empty()
        ftype_points = tuple(ftype_points)
        if points==tuple(range(self.size())) and ftype_points==self._ftype_points:
            return self
        blocks = {xx: _subblock_helper(points, self._blocks[xx]) for xx in self._blocks.keys()}
        new_ftype_points = [points.index(ii) for ii in ftype_points]
        return Flag(self.parent(), len(points), new_ftype_points, **blocks)

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
        cdef Flag newflag
        for ftype_points in itertools.permutations(range(target.size()), self._n):
            sig_check()
            if target.subflag(ftype_points, ftype_points)==self:
                newflag = target.subflag(lrp, ftype_points)
                if newflag not in ret:
                    ret.append(newflag)
        return ret

    cpdef tuple nonequal_permutations(self, include_signature=False):
        cdef list ret = []
        cdef tuple perm
        cdef Flag newflag
        cdef Flag xx
        cdef bint can_add
        for perm in itertools.permutations(range(self.size())):
            sig_check()
            newflag = self.subflag(points=perm)
            can_add = True
            for xx in ret:
                if xx.strong_equal(newflag):
                    can_add = False
                    break
            if can_add:
                if include_signature:
                    ret += newflag.signature_changes()
                else:
                    ret.append(newflag)
        return tuple(ret)
    
    cpdef list signature_changes(self):
        cdef list ret = []
        cdef tuple perms = self.parent()._signature_perms()
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
            ret.append(Flag(self.parent(), self.size(), self._ftype_points, **nblocks))
        return ret
    
    cpdef densities(self, int n1, tuple n1flgs, int n2, tuple n2flgs, \
    list ftype_remap, Flag large_ftype, Flag small_ftype):
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
        cdef Flag n1_subf, n2_subf, ind_large_ftype

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

    def as_flag_algebra_element(self, basis=QQ):
        r"""
        Transforms this `Flag` to a `FlagAlgebraElement` over a given basis

        INPUT:

        - ``basis`` -- Ring (default: `QQ`); the base of
            the FlagAlgebra where the target will live.

        OUTPUT: A `FlagAlgebraElement` representing this `Flag`

        .. SEEALSO::

            :class:`FlagAlgebra`
            :func:`FlagAlgebra._element_constructor_`
            :class:`FlagAlgebraElement`
        """
        from sage.algebras.flag_algebras import FlagAlgebra
        targ_alg = FlagAlgebra(basis, self.theory(), self.ftype())
        return targ_alg(self)
    
    afae = as_flag_algebra_element
    
    def as_operand(self):
        r"""
        Turns this `Flag` into a `FlagAlgebraElement` so operations can be performed on it

        .. SEEALSO::

            :func:`as_flag_algebra_element`
        """
        return self.afae(QQ)

    def _add_(self, other):
        r"""
        Add two Flags together
        
        The flags must have the same ftype. Different sizes are 
            all shifted to the larger one.

        OUTPUT: The :class:`FlagAlgebraElement` object, 
            which is the sum of the two parameters

        EXAMPLES::

        Adding to self is 2*self ::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3)
            sage: g+g==2*g
            True
        
        Adding two distinct elements with the same size gives a vector 
        with exactly two `1` entries ::

            sage: h = GraphTheory(3, edges=[[0, 1]])
            sage: (g+h).values()
            (1, 1, 0, 0)
        
        Adding with different size the smaller flag
        is shifted to have the same size ::
        
            sage: e = GraphTheory(2)
            sage: (e+h).values()
            (1, 5/3, 1/3, 0)

        .. SEEALSO::

            :func:`FlagAlgebraElement._add_`
            :func:`__lshift__`

        """
        if self.ftype()!=other.ftype():
            raise TypeError("The terms must have the same ftype")
        return self.afae()._add_(other.afae())
    
    def _sub_(self, other):
        r"""
        Subtract a Flag from `self`
        
        The flags must have the same ftype. Different sizes are 
            all shifted to the larger one.

        EXAMPLES::
            
            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(2)
            sage: h = GraphTheory(3, edges=[[0, 1]])
            sage: (g-h).values()
            (1, -1/3, 1/3, 0)

        .. SEEALSO::

            :func:`_add_`
            :func:`__lshift__`
            :func:`FlagAlgebraElement._sub_`

        """
        if self.ftype()!=other.ftype():
            raise TypeError("The terms must have the same ftype")
        return self.afae()._sub_(other.afae())
    
    def _mul_(self, other):
        r"""
        Multiply two flags together.
        
        The flags must have the same ftype. The result
        will have the same ftype and size 
        `self.size() + other.size() - self.ftype().size()`

        OUTPUT: The :class:`FlagAlgebraElement` object, 
            which is the product of the two parameters

        EXAMPLES::

        Pointed edge multiplied by itself ::

            sage: from sage.algebras.flag_algebras import *
            sage: pe = GraphTheory(2, edges=[[0, 1]], ftype=[0])
            sage: (pe*pe).values()
            (0, 0, 0, 0, 1, 1)

        .. SEEALSO::

            :func:`FlagAlgebraElement._mul_`
            :func:`mul_project`
            :func:`CombinatorialTheory.mul_project_table`

        TESTS::

            sage: sum((pe*pe*pe*pe).values())
            11
            sage: e = GraphTheory(2)
            sage: (e*e).values()
            (1, 2/3, 1/3, 0, 2/3, 1/3, 0, 0, 1/3, 0, 0)
        """
        if self.ftype()!=other.ftype():
            raise TypeError("The terms must have the same ftype")
        return self.afae()._mul_(other.afae())
    
    def __lshift__(self, amount):
        r"""
        `FlagAlgebraElement`, equal to this, with size is shifted by the amount

        EXAMPLES::

        Edge shifted to size `3` ::

            sage: from sage.algebras.flag_algebras import *
            sage: edge = GraphTheory(2, edges=[[0, 1]])
            sage: (edge<<1).values()
            (0, 1/3, 2/3, 1)

        .. SEEALSO::

            :func:`FlagAlgebraElement.__lshift__`
        """
        return self.afae().__lshift__(amount)
    
    def __truediv__(self, other):
        r"""
        Divide by a scalar

        INPUT:

        - ``other`` -- number; any number such that `1` can be divided with that

        OUTPUT: The `FlagAlgebraElement` resulting from the division

        EXAMPLES::

        Divide by `2` ::

            sage: from sage.algebras.flag_algebras import *
            sage: g = GraphTheory(3)
            sage: (g/2).values()
            (1/2, 0, 0, 0)
            
        Even for `x` symbolic `1/x` is defined, so the division is understood ::
            sage: var('x')
            x
            sage: g = GraphTheory(2)
            sage: g/x
            Flag Algebra Element over Symbolic Ring
            1/x - Flag on 2 points, ftype from [] with edges=[]
            0   - Flag on 2 points, ftype from [] with edges=[[0, 1]]
        
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
        if self.parent()!=other.parent():
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
        if self.parent()!=other.parent():
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

            sage: from sage.algebras.flag_algebras import *
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
              'unique':self._unique,
              'weak_unique':self._weak_unique}
        return dd
    
    def __setstate__(self, dd):
        r"""
        Loads this flag from a dictionary
        """
        self._set_parent(dd['theory'])
        self._n = dd['n']
        self._ftype_points = dd['ftype_points']
        self._ftype_size = len(self._ftype_points)
        self._not_ftype_points = None
        self._blocks = dd['blocks']
        self._unique = dd['unique']
        self._weak_unique = dd['weak_unique']
    
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

            sage: from sage.algebras.flag_algebras import *
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

            sage: from sage.algebras.flag_algebras import *
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

            sage: from sage.algebras.flag_algebras import *
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

cdef class Pattern(Element):
    
    cdef int _n
    cdef int _ftype_size
    
    cdef tuple _ftype_points
    cdef tuple _not_ftype_points
    cdef dict _blocks
    
    cdef Flag _ftype
    
    def __init__(self, theory, n, ftype, **params):
        self._n = int(n)
        self._ftype_points = tuple(ftype)
        self._ftype_size = len(self._ftype_points)
        self._not_ftype_points = None
        self._blocks = _standardize_blocks(params, theory._signature, True)
        self._ftype = None
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
        return 'Pattern on {} points, ftype from {} with {}'.format(self.size(), self.ftype_points(), strblocks)
    
    __str__ = _repr_
    
    def _serialize(self):
        ret = ["pattern", self.size(), self.ftype_points()]
        blocks = self._blocks
        for kk in blocks:
            ret.append(blocks[kk])
        return tuple(ret)
    
    #Basic properties

    def combinatorial_theory(self):
        return self.parent()
    
    theory = combinatorial_theory
    

    def size(self):
        return self._n
    
    vertex_number = size

    cpdef blocks(self, str key=None, bint missing=True):
        if key!=None:
            if missing:
                key += "_m"
            return self._blocks[key]
        if missing:
            return self._blocks
        return {kk:self._blocks[kk] for kk in self.theory().signature().keys()}

    cpdef ftype(self):
        if self._ftype==None:
            blocks = {xx: _subblock_helper(self._ftype_points, self._blocks[xx]) for xx in self.parent().signature().keys()}
            self._ftype = Flag(self.parent(), len(self._ftype_points), self._ftype_points, **blocks)
        return self._ftype
    
    cpdef ftype_points(self):
        return self._ftype_points
    
    cpdef not_ftype_points(self):
        if self._not_ftype_points != None:
            return self._not_ftype_points
        self._not_ftype_points = [ii for ii in range(self.size()) if ii not in self._ftype_points]
        return self._not_ftype_points
    
    def is_ftype(self):
        return False

    #Pattern methods

    cpdef bint is_compatible(self, Flag other):
        if self.size() != other.size():
            return False
        if self.theory() != other.theory():
            return False
        if self.ftype() != other.ftype():
            return False
        
        cdef dict sb = self._blocks
        cdef dict ob
        for operm in other.nonequal_permutations():
            ob = operm.blocks()
            if all([_block_refinement(ob[xx], sb[xx], sb[xx+"_m"]) for xx in self.parent().signature().keys()]):
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
        return Pattern(self.parent(), len(points), new_ftype_points, **blocks)

    #Compatibility with flag algebras

    cpdef list compatible_flags(self):
        aflags = self.parent().generate(self.size())
        ret = [xx for xx in aflags if self.is_compatible(xx)]
        return ret

    def as_flag_algebra_element(self, basis=QQ):
        from sage.algebras.flag_algebras import FlagAlgebra
        from sage.modules.free_module_element import vector

        targ_alg = FlagAlgebra(basis, self.theory(), self.ftype())
        aflags = self.parent().generate(self.size())
        targ_vec = {}
        for ii,xx in enumerate(aflags):
            if self.is_compatible(xx):
                targ_vec[ii]=1
        return targ_alg(self.size(), vector(basis, len(aflags), targ_vec))
    
    afae = as_flag_algebra_element
    
    def as_operand(self):
        return self.afae(QQ)
    
    cpdef _add_(self, other):
        if self.ftype()!=other.ftype():
            raise TypeError("The terms must have the same ftype")
        return self.afae()._add_(other.afae())
    
    cpdef _sub_(self, other):
        if self.ftype()!=other.ftype():
            raise TypeError("The terms must have the same ftype")
        return self.afae()._sub_(other.afae())
    
    cpdef _mul_(self, other):
        if self.ftype()!=other.ftype():
            raise TypeError("The terms must have the same ftype")
        return self.afae()._mul_(other.afae())
    
    def __lshift__(self, amount):
        return self.afae().__lshift__(amount)
    
    def __truediv__(self, other):
        return self.afae().__truediv__(other)
    
    def project(self, ftype_inj=tuple()):
        return self.afae().project(ftype_inj)
    
    def mul_project(self, other, ftype_inj=tuple()):
        return self.afae().mul_project(other, ftype_inj)
    
    def density(self, other):
        safae = self.afae()
        oafae = safae.parent(other)
        return self.afae().density(other)
    