[0;31mType:[0m            module
[0;31mString form:[0m     <module 'sage.combinat.species.library' from '/home/mainak/sage/src/sage/combinat/species/library.py'>
[0;31mFile:[0m            ~/sage/src/sage/combinat/species/library.py
[0;31mSource:[0m         
[0;34m"""[0m
[0;34mExamples of Combinatorial Species[0m
[0;34m"""[0m[0;34m[0m
[0;34m[0m[0;31m# ****************************************************************************[0m[0;34m[0m
[0;34m[0m[0;31m#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,[0m[0;34m[0m
[0;34m[0m[0;31m#[0m[0;34m[0m
[0;34m[0m[0;31m#  Distributed under the terms of the GNU General Public License (GPL)[0m[0;34m[0m
[0;34m[0m[0;31m#[0m[0;34m[0m
[0;34m[0m[0;31m#    This code is distributed in the hope that it will be useful,[0m[0;34m[0m
[0;34m[0m[0;31m#    but WITHOUT ANY WARRANTY; without even the implied warranty of[0m[0;34m[0m
[0;34m[0m[0;31m#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU[0m[0;34m[0m
[0;34m[0m[0;31m#    General Public License for more details.[0m[0;34m[0m
[0;34m[0m[0;31m#[0m[0;34m[0m
[0;34m[0m[0;31m#  The full text of the GPL is available at:[0m[0;34m[0m
[0;34m[0m[0;31m#[0m[0;34m[0m
[0;34m[0m[0;31m#                  https://www.gnu.org/licenses/[0m[0;34m[0m
[0;34m[0m[0;31m# ****************************************************************************[0m[0;34m[0m
[0;34m[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0;34m.[0m[0mset_species[0m [0;32mimport[0m [0mSetSpecies[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0;34m.[0m[0mpartition_species[0m [0;32mimport[0m [0mPartitionSpecies[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0;34m.[0m[0msubset_species[0m [0;32mimport[0m [0mSubsetSpecies[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0;34m.[0m[0mrecursive_species[0m [0;32mimport[0m [0mCombinatorialSpecies[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0;34m.[0m[0mcharacteristic_species[0m [0;32mimport[0m [0mCharacteristicSpecies[0m[0;34m,[0m [0mSingletonSpecies[0m[0;34m,[0m [0mEmptySetSpecies[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0;34m.[0m[0mcycle_species[0m [0;32mimport[0m [0mCycleSpecies[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0;34m.[0m[0mlinear_order_species[0m [0;32mimport[0m [0mLinearOrderSpecies[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0;34m.[0m[0mpermutation_species[0m [0;32mimport[0m [0mPermutationSpecies[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0;34m.[0m[0mempty_species[0m [0;32mimport[0m [0mEmptySpecies[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0;34m.[0m[0msum_species[0m [0;32mimport[0m [0mSumSpecies[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0;34m.[0m[0mproduct_species[0m [0;32mimport[0m [0mProductSpecies[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0;34m.[0m[0mcomposition_species[0m [0;32mimport[0m [0mCompositionSpecies[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0;34m.[0m[0mfunctorial_composition_species[0m [0;32mimport[0m [0mFunctorialCompositionSpecies[0m[0;34m[0m
[0;34m[0m[0;34m[0m
[0;34m[0m[0;32mfrom[0m [0msage[0m[0;34m.[0m[0mmisc[0m[0;34m.[0m[0mcachefunc[0m [0;32mimport[0m [0mcached_function[0m[0;34m[0m
[0;34m[0m[0;34m[0m
[0;34m[0m[0;34m[0m
[0;34m[0m[0;34m@[0m[0mcached_function[0m[0;34m[0m
[0;34m[0m[0;32mdef[0m [0mSimpleGraphSpecies[0m[0;34m([0m[0;34m)[0m[0;34m:[0m[0;34m[0m
[0;34m[0m    [0;34m"""[0m
[0;34m    Return the species of simple graphs.[0m
[0;34m[0m
[0;34m    EXAMPLES::[0m
[0;34m[0m
[0;34m        sage: S = species.SimpleGraphSpecies()[0m
[0;34m        sage: S.generating_series().counts(10)[0m
[0;34m        [1, 1, 2, 8, 64, 1024, 32768, 2097152, 268435456, 68719476736][0m
[0;34m        sage: S.cycle_index_series()[:5][0m
[0;34m        [p[],[0m
[0;34m         p[1],[0m
[0;34m         p[1, 1] + p[2],[0m
[0;34m         4/3*p[1, 1, 1] + 2*p[2, 1] + 2/3*p[3],[0m
[0;34m         8/3*p[1, 1, 1, 1] + 4*p[2, 1, 1] + 2*p[2, 2] + 4/3*p[3, 1] + p[4]][0m
[0;34m        sage: S.isotype_generating_series()[:6][0m
[0;34m        [1, 1, 2, 4, 11, 34][0m
[0;34m[0m
[0;34m    TESTS::[0m
[0;34m[0m
[0;34m        sage: seq = S.isotype_generating_series().counts(6)[1:][0m
[0;34m        sage: oeis(seq)[0]                              # optional -- internet[0m
[0;34m        A000088: Number of graphs on n unlabeled nodes.[0m
[0;34m[0m
[0;34m    ::[0m
[0;34m[0m
[0;34m        sage: seq = S.generating_series().counts(10)[1:][0m
[0;34m        sage: oeis(seq)[0]                              # optional -- internet[0m
[0;34m        A006125: a(n) = 2^(n*(n-1)/2).[0m
[0;34m    """[0m[0;34m[0m
[0;34m[0m    [0mE[0m [0;34m=[0m [0mSetSpecies[0m[0;34m([0m[0;34m)[0m[0;34m[0m
[0;34m[0m    [0mE2[0m [0;34m=[0m [0mSetSpecies[0m[0;34m([0m[0msize[0m[0;34m=[0m[0;36m2[0m[0;34m)[0m[0;34m[0m
[0;34m[0m    [0mWP[0m [0;34m=[0m [0mSubsetSpecies[0m[0;34m([0m[0;34m)[0m[0;34m[0m
[0;34m[0m    [0mP2[0m [0;34m=[0m [0mE2[0m [0;34m*[0m [0mE[0m[0;34m[0m
[0;34m[0m    [0;32mreturn[0m [0mWP[0m[0;34m.[0m[0mfunctorial_composition[0m[0;34m([0m[0mP2[0m[0;34m)[0m[0;34m[0m
[0;34m[0m[0;34m[0m
[0;34m[0m[0;34m[0m
[0;34m[0m[0;34m@[0m[0mcached_function[0m[0;34m[0m
[0;34m[0m[0;32mdef[0m [0mBinaryTreeSpecies[0m[0;34m([0m[0;34m)[0m[0;34m:[0m[0;34m[0m
[0;34m[0m    [0;34mr"""[0m
[0;34m    Return the species of binary trees on `n` leaves.[0m
[0;34m[0m
[0;34m    The species of binary trees `B` is defined by `B = X + B \cdot B`,[0m
[0;34m    where `X` is the singleton species.[0m
[0;34m[0m
[0;34m    EXAMPLES::[0m
[0;34m[0m
[0;34m        sage: B = species.BinaryTreeSpecies()[0m
[0;34m        sage: B.generating_series().counts(10)[0m
[0;34m        [0, 1, 2, 12, 120, 1680, 30240, 665280, 17297280, 518918400][0m
[0;34m        sage: B.isotype_generating_series().counts(10)[0m
[0;34m        [0, 1, 1, 2, 5, 14, 42, 132, 429, 1430][0m
[0;34m        sage: B._check()[0m
[0;34m        True[0m
[0;34m[0m
[0;34m    ::[0m
[0;34m[0m
[0;34m        sage: B = species.BinaryTreeSpecies()[0m
[0;34m        sage: a = B.structures([1,2,3,4,5])[187]; a[0m
[0;34m        2*((5*3)*(4*1))[0m
[0;34m        sage: a.automorphism_group()[0m
[0;34m        Permutation Group with generators [()][0m
[0;34m[0m
[0;34m    TESTS::[0m
[0;34m[0m
[0;34m        sage: seq = B.isotype_generating_series().counts(10)[1:][0m
[0;34m        sage: oeis(seq)[0]                              # optional -- internet[0m
[0;34m        A000108: Catalan numbers: ...[0m
[0;34m    """[0m[0;34m[0m
[0;34m[0m    [0mB[0m [0;34m=[0m [0mCombinatorialSpecies[0m[0;34m([0m[0mmin[0m[0;34m=[0m[0;36m1[0m[0;34m)[0m[0;34m[0m
[0;34m[0m    [0mX[0m [0;34m=[0m [0mSingletonSpecies[0m[0;34m([0m[0;34m)[0m[0;34m[0m
[0;34m[0m    [0mB[0m[0;34m.[0m[0mdefine[0m[0;34m([0m[0mX[0m [0;34m+[0m [0mB[0m [0;34m*[0m [0mB[0m[0;34m)[0m[0;34m[0m
[0;34m[0m    [0;32mreturn[0m [0mB[0m[0;34m[0m
[0;34m[0m[0;34m[0m
[0;34m[0m[0;34m[0m
[0;34m[0m[0;34m@[0m[0mcached_function[0m[0;34m[0m
[0;34m[0m[0;32mdef[0m [0mBinaryForestSpecies[0m[0;34m([0m[0;34m)[0m[0;34m:[0m[0;34m[0m
[0;34m[0m    [0;34m"""[0m
[0;34m    Return the species of binary forests.[0m
[0;34m[0m
[0;34m    Binary forests are defined as sets of binary trees.[0m
[0;34m[0m
[0;34m    EXAMPLES::[0m
[0;34m[0m
[0;34m        sage: F = species.BinaryForestSpecies()[0m
[0;34m        sage: F.generating_series().counts(10)[0m
[0;34m        [1, 1, 3, 19, 193, 2721, 49171, 1084483, 28245729, 848456353][0m
[0;34m        sage: F.isotype_generating_series().counts(10)[0m
[0;34m        [1, 1, 2, 4, 10, 26, 77, 235, 758, 2504][0m
[0;34m        sage: F.cycle_index_series()[:7][0m
[0;34m        [p[],[0m
[0;34m         p[1],[0m
[0;34m         3/2*p[1, 1] + 1/2*p[2],[0m
[0;34m         19/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3],[0m
[0;34m         193/24*p[1, 1, 1, 1] + 3/4*p[2, 1, 1] + 5/8*p[2, 2] + 1/3*p[3, 1] + 1/4*p[4],[0m
[0;34m         907/40*p[1, 1, 1, 1, 1] + 19/12*p[2, 1, 1, 1] + 5/8*p[2, 2, 1] + 1/2*p[3, 1, 1] + 1/6*p[3, 2] + 1/4*p[4, 1] + 1/5*p[5],[0m
[0;34m         49171/720*p[1, 1, 1, 1, 1, 1] + 193/48*p[2, 1, 1, 1, 1] + 15/16*p[2, 2, 1, 1] + 61/48*p[2, 2, 2] + 19/18*p[3, 1, 1, 1] + 1/6*p[3, 2, 1] + 7/18*p[3, 3] + 3/8*p[4, 1, 1] + 1/8*p[4, 2] + 1/5*p[5, 1] + 1/6*p[6]][0m
[0;34m[0m
[0;34m    TESTS::[0m
[0;34m[0m
[0;34m        sage: seq = F.isotype_generating_series().counts(10)[1:][0m
[0;34m        sage: oeis(seq)[0]                              # optional -- internet[0m
[0;34m        A052854: Number of forests of ordered trees on n total nodes.[0m
[0;34m    """[0m[0;34m[0m
[0;34m[0m    [0mB[0m [0;34m=[0m [0mBinaryTreeSpecies[0m[0;34m([0m[0;34m)[0m[0;34m[0m
[0;34m[0m    [0mS[0m [0;34m=[0m [0mSetSpecies[0m[0;34m([0m[0;34m)[0m[0;34m[0m
[0;34m[0m    [0mF[0m [0;34m=[0m [0mS[0m[0;34m([0m[0mB[0m[0;34m)[0m[0;34m[0m
[0;34m[0m    [0;32mreturn[0m [0mF[0m[0;34m[0m
[0;34m[0m[0;34m[0m
[0;34m[0m[0;34m[0m
[0;34m[0m[0;32mdel[0m [0mcached_function[0m  [0;31m# so it doesn't get picked up by tab completion[0m[0;34m[0m[0;34m[0m[0m
[0;31mClass docstring:[0m
Create a module object.

The name must be a string; the optional doc argument can have any
type.
[0;31mInit docstring:[0m  Initialize self.  See help(type(self)) for accurate signature.
