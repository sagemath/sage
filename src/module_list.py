#!/usr/bin/env python

import os
from distutils.extension import Extension
from sage.env import SAGE_LOCAL

SAGE_INC = os.path.join(SAGE_LOCAL, 'include')

#########################################################
### BLAS setup
#########################################################

## Choose cblas library -- note -- make sure to update sage/misc/cython.py
## if you change this!!
if os.environ.has_key('SAGE_BLAS'):
    BLAS=os.environ['SAGE_BLAS']
    BLAS2=os.environ['SAGE_BLAS']
elif os.path.exists('%s/lib/libatlas.so'%os.environ['SAGE_LOCAL']):
    BLAS='cblas'
    BLAS2='atlas'
elif os.path.exists('/usr/lib/libcblas.dylib') or \
     os.path.exists('/usr/lib/libcblas.so'):
    BLAS='cblas'
    BLAS2='cblas'
elif os.path.exists('/usr/lib/libblas.dll.a'):
    BLAS='gslcblas'
    BLAS2='gslcblas'
else:
    # This is very slow  (?), but *guaranteed* to be available.
    BLAS='gslcblas'
    BLAS2='gslcblas'


#########################################################
### Commonly used definitions
#########################################################

flint_depends = [SAGE_INC + '/flint/flint.h']
singular_depends = [SAGE_INC + '/libsingular.h']
givaro_depends = [SAGE_INC + '/givaro/givconfig.h']

singular_incs = [SAGE_INC + '/singular', SAGE_INC + '/factory']

#########################################################
### M4RI flags
#########################################################

import ast
m4ri_extra_compile_args = ["-std=c99"]
for line in open(SAGE_INC + "/m4ri/m4ri_config.h"):
    if not line.startswith("#define __M4RI_SIMD_CFLAGS"):
        continue
    m4ri_sse2_cflags = ast.literal_eval(line[len("#define __M4RI_SIMD_CFLAGS"):].strip())
    m4ri_extra_compile_args.extend( [flag.strip() for flag in m4ri_sse2_cflags.split(" ") if flag.strip()] )
    break

singular_libs = ['singular', 'flint', 'ntl', 'gmpxx', 'gmp', 'readline', 'm']

#########################################################
### Givaro flags
#########################################################

givaro_extra_compile_args =['-D__STDC_LIMIT_MACROS']

#########################################################
### PolyBoRi settings
#########################################################

polybori_extra_compile_args = []
polybori_major_version = '0.8'

#############################################################
### List of modules
###
### Note that the list of modules is sorted alphabetically
### by extension name. Please keep this list sorted when
### adding new modules!
###
#############################################################

from sage_setup.optional_extension import OptionalExtension
UNAME = os.uname()

def uname_specific(name, value, alternative):
    if name in UNAME[0]:
        return value
    else:
        return alternative


ext_modules = [

    ################################
    ##
    ## sage.algebras
    ##
    ################################

    Extension('sage.algebras.quatalg.quaternion_algebra_element',
               sources = ['sage/algebras/quatalg/quaternion_algebra_element.pyx'],
               language='c++',
               libraries = ["flint", "gmp", "gmpxx", "m", "stdc++", "ntl"],
               depends = flint_depends),

    Extension('sage.algebras.letterplace.free_algebra_letterplace',
              sources = ['sage/algebras/letterplace/free_algebra_letterplace.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs,
              depends = singular_depends),

    Extension('sage.algebras.letterplace.free_algebra_element_letterplace',
              sources = ['sage/algebras/letterplace/free_algebra_element_letterplace.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs,
              depends = singular_depends),

    Extension('sage.algebras.letterplace.letterplace_ideal',
              sources = ['sage/algebras/letterplace/letterplace_ideal.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs,
              depends = singular_depends),

    Extension('sage.algebras.quatalg.quaternion_algebra_cython',
               sources = ['sage/algebras/quatalg/quaternion_algebra_cython.pyx'],
               language='c++',
               depends = flint_depends,
               libraries = ["flint", "gmp", "gmpxx", "m", "stdc++", "ntl"]),

    ################################
    ##
    ## sage.calculus
    ##
    ################################

    Extension('sage.calculus.var',
              sources = ['sage/calculus/var.pyx']),

    Extension('sage.calculus.riemann',
              sources = ['sage/calculus/riemann.pyx']),

    Extension('sage.calculus.interpolators',
              sources = ['sage/calculus/interpolators.pyx']),


    ################################
    ##
    ## sage.categories
    ##
    ################################

    Extension('*', ['sage/categories/**/*.pyx']),

    ################################
    ##
    ## sage.coding
    ##
    ################################

    Extension('sage.coding.codecan.codecan',
              sources = ['sage/coding/codecan/codecan.pyx'],
              libraries = ['gmp', 'flint'],
              include_dirs = ['sage/groups/perm_gps/partn_ref2/'],
              depends = flint_depends),

    Extension('*', ['sage/coding/**/*.pyx']),

    ################################
    ##
    ## sage.combinat
    ##
    ################################

    Extension('sage.combinat.expnums',
              sources = ['sage/combinat/expnums.pyx'],
              libraries = ['gmp']),

    Extension('sage.combinat.matrices.dancing_links',
              sources = ['sage/combinat/matrices/dancing_links.pyx'],
              language='c++'),

    Extension('sage.structure.list_clone',
              sources=['sage/structure/list_clone.pyx']),

    Extension('sage.structure.list_clone_demo',
              sources=['sage/structure/list_clone_demo.pyx']),

    Extension('sage.structure.list_clone_timings_cy',
              sources=['sage/structure/list_clone_timings_cy.pyx']),

    Extension('sage.sets.finite_set_map_cy',
              sources=['sage/sets/finite_set_map_cy.pyx']),

    Extension('sage.combinat.partitions',
              sources = ['sage/combinat/partitions.pyx',
                         'sage/combinat/partitions_c.cc'],
              libraries = ['gmp', 'mpfr'],
              depends = ['sage/combinat/partitions_c.h'],
              language='c++'),

    Extension('sage.combinat.words.word_datatypes',
            sources=['sage/combinat/words/word_datatypes.pyx']),

    Extension('sage.combinat.words.word_char',
            sources=['sage/combinat/words/word_char.pyx']),

    Extension('sage.combinat.permutation_cython',
              sources=['sage/combinat/permutation_cython.pyx']),

    Extension('sage.combinat.dict_addition',
              sources=['sage/combinat/dict_addition.pyx']),

    Extension('sage.combinat.debruijn_sequence',
              sources=['sage/combinat/debruijn_sequence.pyx']),

    Extension('sage.combinat.degree_sequences',
              sources = ['sage/combinat/degree_sequences.pyx']),

    Extension('sage.combinat.combinat_cython',
              sources=['sage/combinat/combinat_cython.pyx'],
              libraries=['gmp']),

    Extension('sage.combinat.enumeration_mod_permgroup',
              sources=['sage/combinat/enumeration_mod_permgroup.pyx']),

    Extension('sage.combinat.q_bernoulli',
              sources = ['sage/combinat/q_bernoulli.pyx']),

    Extension('sage.combinat.crystals.letters',
              sources=['sage/combinat/crystals/letters.pyx']),

    Extension('sage.combinat.designs.subhypergraph_search',
              sources=['sage/combinat/designs/subhypergraph_search.pyx']),

    Extension('sage.combinat.designs.designs_pyx',
              sources=['sage/combinat/designs/designs_pyx.pyx'],
              libraries=['gmp']),

    Extension('sage.combinat.designs.orthogonal_arrays_find_recursive',
              sources=['sage/combinat/designs/orthogonal_arrays_find_recursive.pyx']),

    Extension('sage.combinat.designs.evenly_distributed_sets',
              sources=['sage/combinat/designs/evenly_distributed_sets.pyx']),

    ################################
    ##
    ## sage.crypto
    ##
    ################################

    Extension('sage.crypto.boolean_function',
              sources = ['sage/crypto/boolean_function.pyx'],
              libraries=['gmp']),


    ################################
    ##
    ## sage.data_structures
    ##
    ################################

    Extension('sage.data_structures.bounded_integer_sequences',
              sources = ['sage/data_structures/bounded_integer_sequences.pyx'],
              libraries = ['gmp']),

    Extension('sage.data_structures.bitset',
              sources = ['sage/data_structures/bitset.pyx'],
              libraries = ['gmp']),

    ################################
    ##
    ## sage.ext
    ##
    ################################

    Extension('*', ['sage/ext/*.pyx']),

    ################################
    ##
    ## sage.finance
    ##
    ################################

    Extension('sage.finance.fractal',
              sources = ['sage/finance/fractal.pyx']),

    Extension('sage.finance.markov_multifractal_cython',
              sources = ['sage/finance/markov_multifractal_cython.pyx']),

    Extension('sage.finance.option',
              sources = ['sage/finance/option.pyx']),

    Extension('sage.finance.time_series',
              sources = ['sage/finance/time_series.pyx']),


    ################################
    ##
    ## sage.functions
    ##
    ################################

    Extension('sage.functions.prime_pi',
        sources = ['sage/functions/prime_pi.pyx'],
        libraries = ['pari', 'gmp'],
        extra_compile_args = ['-std=c99']),

    ################################
    ##
    ## sage.games
    ##
    ################################

    Extension('sage.games.sudoku_backtrack',
              sources = ['sage/games/sudoku_backtrack.pyx']),

    ################################
    ##
    ## sage.geometry
    ##
    ################################

    Extension('sage.geometry.point_collection',
              sources = ['sage/geometry/point_collection.pyx']),

    Extension('sage.geometry.toric_lattice_element',
              sources = ['sage/geometry/toric_lattice_element.pyx'],
              libraries=['gmp']),

    Extension('sage.geometry.integral_points',
              sources = ['sage/geometry/integral_points.pyx']),

    Extension('sage.geometry.triangulation.base',
              sources = ['sage/geometry/triangulation/base.pyx',
                         'sage/geometry/triangulation/functions.cc',
                         'sage/geometry/triangulation/data.cc',
                         'sage/geometry/triangulation/triangulations.cc'],
              depends = ['sage/geometry/triangulation/functions.h',
                         'sage/geometry/triangulation/data.h',
                         'sage/geometry/triangulation/triangulations.h'],
              language="c++"),

    ################################
    ##
    ## sage.graphs
    ##
    ################################

    Extension('sage.graphs.asteroidal_triples',
              sources = ['sage/graphs/asteroidal_triples.pyx']),

    Extension('sage.graphs.chrompoly',
              sources = ['sage/graphs/chrompoly.pyx'],
              libraries = ['gmp']),

    Extension('sage.graphs.cliquer',
              sources = ['sage/graphs/cliquer.pyx', 'sage/graphs/cliquer/cl.c'],
              libraries = ['cliquer']),

    Extension('sage.graphs.centrality',
              sources = ['sage/graphs/centrality.pyx']),

    Extension('sage.graphs.independent_sets',
              sources = ['sage/graphs/independent_sets.pyx'],
              libraries=['gmp']),

    Extension('sage.graphs.graph_decompositions.vertex_separation',
              sources = ['sage/graphs/graph_decompositions/vertex_separation.pyx']),

    Extension('sage.graphs.graph_decompositions.graph_products',
              sources = ['sage/graphs/graph_decompositions/graph_products.pyx']),

    Extension('sage.graphs.convexity_properties',
              sources = ['sage/graphs/convexity_properties.pyx'],
              libraries = ['gmp']),

    Extension('sage.graphs.comparability',
              sources = ['sage/graphs/comparability.pyx']),

    Extension('sage.graphs.generic_graph_pyx',
              sources = ['sage/graphs/generic_graph_pyx.pyx'],
              libraries = ['gmp']),

    Extension('sage.graphs.graph_generators_pyx',
              sources = ['sage/graphs/graph_generators_pyx.pyx']),

    Extension('sage.graphs.distances_all_pairs',
              sources = ['sage/graphs/distances_all_pairs.pyx'],
              libraries = ['gmp']),

    Extension('sage.graphs.base.graph_backends',
              sources = ['sage/graphs/base/graph_backends.pyx']),

    Extension('sage.graphs.base.static_dense_graph',
              sources = ['sage/graphs/base/static_dense_graph.pyx'],
              libraries = ['gmp']),

    Extension('sage.graphs.base.static_sparse_graph',
              sources = ['sage/graphs/base/static_sparse_graph.pyx'],
              libraries = ['gmp']),

    Extension('sage.graphs.base.static_sparse_backend',
              sources = ['sage/graphs/base/static_sparse_backend.pyx']),

    Extension('sage.graphs.weakly_chordal',
              sources = ['sage/graphs/weakly_chordal.pyx']),

    Extension('sage.graphs.matchpoly',
              sources = ['sage/graphs/matchpoly.pyx'],
              libraries = ['gmp', 'flint'],
              extra_compile_args = ['-std=c99'],
              depends = flint_depends),

    OptionalExtension("sage.graphs.mcqd",
              ["sage/graphs/mcqd.pyx"],
              language = "c++",
              package = 'mcqd'),

    OptionalExtension("sage.graphs.bliss",
              ["sage/graphs/bliss.pyx"],
              language = "c++",
              libraries = ['bliss'],
              package = 'bliss'),

    OptionalExtension('sage.graphs.modular_decomposition',
              sources = ['sage/graphs/modular_decomposition.pyx'],
              libraries = ['modulardecomposition'],
              package = 'modular_decomposition'),

    Extension('sage.graphs.planarity',
              sources = ['sage/graphs/planarity.pyx',
                         'sage/graphs/planarity_c/graphColorVertices.c',
                         'sage/graphs/planarity_c/graphColorVertices_Extensions.c',
                         'sage/graphs/planarity_c/graphDrawPlanar.c',
                         'sage/graphs/planarity_c/graphDrawPlanar_Extensions.c',
                         'sage/graphs/planarity_c/graphEmbed.c',
                         'sage/graphs/planarity_c/graphExtensions.c',
                         'sage/graphs/planarity_c/graphIO.c',
                         'sage/graphs/planarity_c/graphIsolator.c',
                         'sage/graphs/planarity_c/graphK23Search.c',
                         'sage/graphs/planarity_c/graphK23Search_Extensions.c',
                         'sage/graphs/planarity_c/graphK33Search.c',
                         'sage/graphs/planarity_c/graphK33Search_Extensions.c',
                         'sage/graphs/planarity_c/graphK4Search.c',
                         'sage/graphs/planarity_c/graphK4Search_Extensions.c',
                         'sage/graphs/planarity_c/graphNonplanar.c',
                         'sage/graphs/planarity_c/graphOuterplanarObstruction.c',
                         'sage/graphs/planarity_c/graphPreprocess.c',
                         'sage/graphs/planarity_c/graphTests.c',
                         'sage/graphs/planarity_c/graphUtils.c',
                         'sage/graphs/planarity_c/listcoll.c',
                         'sage/graphs/planarity_c/planarity.c',
                         'sage/graphs/planarity_c/planarityCommandLine.c',
                         'sage/graphs/planarity_c/planarityRandomGraphs.c',
                         'sage/graphs/planarity_c/planaritySpecificGraph.c',
                         'sage/graphs/planarity_c/planarityUtils.c',
                         'sage/graphs/planarity_c/stack.c'],
              depends = ['sage/graphs/planarity_c/appconst.h',
                         'sage/graphs/planarity_c/graphColorVertices.h',
                         'sage/graphs/planarity_c/graphColorVertices.private.h',
                         'sage/graphs/planarity_c/graphDrawPlanar.h',
                         'sage/graphs/planarity_c/graphDrawPlanar.private.h',
                         'sage/graphs/planarity_c/graphExtensions.h',
                         'sage/graphs/planarity_c/graphExtensions.private.h',
                         'sage/graphs/planarity_c/graphFunctionTable.h',
                         'sage/graphs/planarity_c/graph.h',
                         'sage/graphs/planarity_c/graphK23Search.h',
                         'sage/graphs/planarity_c/graphK23Search.private.h',
                         'sage/graphs/planarity_c/graphK33Search.h',
                         'sage/graphs/planarity_c/graphK33Search.private.h',
                         'sage/graphs/planarity_c/graphK4Search.h',
                         'sage/graphs/planarity_c/graphK4Search.private.h',
                         'sage/graphs/planarity_c/graphStructures.h',
                         'sage/graphs/planarity_c/listcoll.h',
                         'sage/graphs/planarity_c/planarity.h',
                         'sage/graphs/planarity_c/platformTime.h',
                         'sage/graphs/planarity_c/stack.h']),

    Extension('sage.graphs.graph_decompositions.rankwidth',
              sources = ['sage/graphs/graph_decompositions/rankwidth.pyx',
                         'sage/graphs/graph_decompositions/rankwidth_c/rw.c']),

    Extension('sage.graphs.graph_decompositions.bandwidth',
              sources = ['sage/graphs/graph_decompositions/bandwidth.pyx']),

    Extension('sage.graphs.spanning_tree',
              sources = ['sage/graphs/spanning_tree.pyx']),

    Extension('sage.graphs.trees',
              sources = ['sage/graphs/trees.pyx']),

    Extension('sage.graphs.genus',
              sources = ['sage/graphs/genus.pyx']),

    Extension('sage.graphs.hyperbolicity',
              sources = ['sage/graphs/hyperbolicity.pyx'],
              libraries = ['gmp']),

    ################################
    ##
    ## sage.graphs.base
    ##
    ################################

    Extension('sage.graphs.base.c_graph',
              sources = ['sage/graphs/base/c_graph.pyx'],
              libraries=['gmp']),

    Extension('sage.graphs.base.sparse_graph',
              sources = ['sage/graphs/base/sparse_graph.pyx'],
              libraries=['gmp']),

    Extension('sage.graphs.base.dense_graph',
              sources = ['sage/graphs/base/dense_graph.pyx'],
              libraries=['gmp']),

    ################################
    ##
    ## sage.groups
    ##
    ################################

    Extension('sage.groups.group',
              sources = ['sage/groups/group.pyx']),

    Extension('sage.groups.old',
              sources = ['sage/groups/old.pyx']),

    Extension('sage.groups.libgap_wrapper',
              sources = ['sage/groups/libgap_wrapper.pyx']),

    Extension('sage.groups.perm_gps.permgroup_element',
              sources = ['sage/groups/perm_gps/permgroup_element.pyx']),

    Extension('sage.groups.semimonomial_transformations.semimonomial_transformation',
              sources = ['sage/groups/semimonomial_transformations/semimonomial_transformation.pyx']),

    ###################################
    ##
    ## sage.groups.perm_gps.partn_ref
    ##
    ###################################

    Extension('sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label',
              sources = ['sage/groups/perm_gps/partn_ref/automorphism_group_canonical_label.pyx'],
              libraries = ['gmp', 'flint'],
              extra_compile_args = ['-std=c99'],
              depends = flint_depends),

    Extension('sage.groups.perm_gps.partn_ref.canonical_augmentation',
              sources = ['sage/groups/perm_gps/partn_ref/canonical_augmentation.pyx'],
              libraries = ['gmp', 'flint'],
              extra_compile_args = ['-std=c99'],
              depends = flint_depends),

    Extension('sage.groups.perm_gps.partn_ref.double_coset',
              sources = ['sage/groups/perm_gps/partn_ref/double_coset.pyx'],
              libraries = ['gmp', 'flint'],
              extra_compile_args = ['-std=c99'],
              depends = flint_depends),

    Extension('sage.groups.perm_gps.partn_ref.refinement_binary',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_binary.pyx'],
              libraries = ['gmp', 'flint'],
              extra_compile_args = ['-std=c99'],
              depends = flint_depends),

    Extension('sage.groups.perm_gps.partn_ref.refinement_graphs',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_graphs.pyx'],
              libraries = ['gmp', 'flint'],
              extra_compile_args = ['-std=c99'],
              depends = flint_depends),

    Extension('sage.groups.perm_gps.partn_ref.refinement_lists',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_lists.pyx'],
              libraries = ['gmp', 'flint'],
              extra_compile_args = ['-std=c99'],
              depends = flint_depends),

    Extension('sage.groups.perm_gps.partn_ref.refinement_matrices',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_matrices.pyx'],
              libraries = ['gmp', 'flint'],
              extra_compile_args = ['-std=c99'],
              depends = flint_depends),

    Extension('sage.groups.perm_gps.partn_ref.refinement_python',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_python.pyx'],
              libraries = ['gmp', 'flint'],
              extra_compile_args = ['-std=c99'],
              depends = flint_depends),

    Extension('sage.groups.perm_gps.partn_ref.refinement_sets',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_sets.pyx'],
              libraries = ['gmp', 'flint'],
              extra_compile_args = ['-std=c99'],
              depends = flint_depends),

    ###################################
    ##
    ## sage.groups.perm_gps.partn_ref2
    ##
    ###################################

    Extension('sage.groups.perm_gps.partn_ref2.refinement_generic',
              sources = ['sage/groups/perm_gps/partn_ref2/refinement_generic.pyx'],
              libraries = ["flint", "gmp", "gmpxx", "stdc++"],
              extra_compile_args=["-std=c99"],
              depends = flint_depends +
                  ['sage/groups/perm_gps/partn_ref2/refinement_generic.h']),

    ################################
    ##
    ## sage.gsl
    ##
    ################################

    Extension('sage.gsl.callback',
              sources = ['sage/gsl/callback.pyx'],
              libraries = ['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.dwt',
              sources = ['sage/gsl/dwt.pyx'],
              libraries=['gsl',BLAS],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.fft',
              sources = ['sage/gsl/fft.pyx'],
              libraries = ['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.gsl_array',
              sources = ['sage/gsl/gsl_array.pyx'],
              libraries=['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.integration',
              sources = ['sage/gsl/integration.pyx'],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')],
              libraries=['gsl',BLAS, BLAS2]),

    Extension('sage.gsl.interpolation',
              sources = ['sage/gsl/interpolation.pyx'],
              libraries = ['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.ode',
              sources = ['sage/gsl/ode.pyx'],
              libraries=['gsl',BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.gsl.probability_distribution',
              sources = ['sage/gsl/probability_distribution.pyx'],
              libraries=['gsl', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    ################################
    ##
    ## sage.interacts
    ##
    ################################

    Extension('sage.interacts.library_cython',
              sources = ['sage/interacts/library_cython.pyx'],
              libraries = []),

    ################################
    ##
    ## sage.libs
    ##
    ################################

    OptionalExtension('sage.libs.coxeter3.coxeter',
              sources = ['sage/libs/coxeter3/coxeter.pyx'],
              include_dirs = [os.path.join(SAGE_INC, 'coxeter')],
              language="c++",
              libraries = ['coxeter3'],
              package = 'coxeter3'),

    Extension('sage.libs.ecl',
              sources = ["sage/libs/ecl.pyx"],
              libraries = ["ecl", "gmp"],
              include_dirs = [SAGE_INC + '/ecl'],
              depends = [SAGE_INC + '/ecl/ecl.h']),

    OptionalExtension("sage.libs.fes",
             ["sage/libs/fes.pyx"],
             language = "c",
             libraries = ['fes'],
             package = 'fes'),

    Extension('sage.libs.flint.flint',
              sources = ["sage/libs/flint/flint.pyx"],
              libraries = ["flint", "gmp", "gmpxx", "m", "stdc++"],
              extra_compile_args = ["-std=c99", "-D_XPG6"],
              depends = flint_depends),

    Extension('sage.libs.flint.fmpz_poly',
              sources = ["sage/libs/flint/fmpz_poly.pyx"],
              libraries = ["flint", "gmp", "gmpxx", "m", "stdc++"],
              extra_compile_args = ["-std=c99", "-D_XPG6"],
              depends = flint_depends),

    Extension('sage.libs.flint.arith',
              sources = ["sage/libs/flint/arith.pyx"],
              libraries = ["flint", "gmp", "gmpxx", "m", "stdc++"],
              extra_compile_args = ["-std=c99", "-D_XPG6"],
              depends = flint_depends),

    Extension('sage.libs.fplll.fplll',
              sources = ['sage/libs/fplll/fplll.pyx'],
              libraries = ['gmp', 'mpfr', 'stdc++', 'fplll'],
              language="c++",
              include_dirs = [SAGE_INC + '/fplll'],
              extra_compile_args=["-DFPLLL_V3_COMPAT"],
              depends = [SAGE_INC + "/fplll/fplll.h"] + flint_depends),

    Extension('sage.libs.gmp.pylong',
              sources = ['sage/libs/gmp/pylong.pyx'],
              libraries = ['gmp']),

    Extension('sage.libs.gmp.rational_reconstruction',
              sources = ['sage/libs/gmp/rational_reconstruction.pyx'],
              libraries = ['gmp']),

    Extension('sage.libs.linbox.linbox',
              sources = ['sage/libs/linbox/linbox.pyx'],
              # For this to work on cygwin, linboxsage *must* be
              # before ntl.
              libraries = ['linboxsage', 'ntl', 'iml', 'linbox',
                           'stdc++', 'givaro', 'mpfr', 'gmp', 'gmpxx', BLAS, BLAS2],
              language = 'c++',
              extra_compile_args = givaro_extra_compile_args,
              depends = givaro_depends),

    Extension('sage.libs.lcalc.lcalc_Lfunction',
              sources = ['sage/libs/lcalc/lcalc_Lfunction.pyx'],
              libraries = ['m', 'ntl', 'mpfr', 'gmp', 'gmpxx',
                           'Lfunction'],
              include_dirs = [SAGE_INC + "/libLfunction"],
              extra_compile_args=["-O3", "-ffast-math"],
              language = 'c++'),

    Extension('sage.libs.libecm',
              sources = ['sage/libs/libecm.pyx'],
              libraries = ['ecm', 'gmp'],
              extra_link_args = uname_specific("Linux", ["-Wl,-z,noexecstack"],
                                                        []),
              depends = [SAGE_INC + "/ecm.h"]),

    Extension('sage.libs.lrcalc.lrcalc',
              sources = ["sage/libs/lrcalc/lrcalc.pyx"],
              include_dirs = [SAGE_INC + '/lrcalc/'],
              libraries = ["lrcalc"]),

    Extension('sage.libs.mwrank.mwrank',
              sources = ["sage/libs/mwrank/mwrank.pyx",
                         "sage/libs/mwrank/wrap.cc"],
              define_macros = [("NTL_ALL",None)],
              depends = ["sage/libs/mwrank/wrap.h"] +
                        [ SAGE_INC + "/eclib/" + h for h in
                          ["curve.h","egr.h","descent.h","points.h","isogs.h",
                            "marith.h","htconst.h","interface.h"]
                        ],
              libraries = ["ec",
                           "ntl", "pari", "gmp", "gmpxx", "stdc++", "m"]),

    Extension('sage.libs.pari.gen',
              sources = ["sage/libs/pari/gen.pyx"],
              libraries = ['pari', 'gmp']),

    Extension('sage.libs.pari.handle_error',
              sources = ["sage/libs/pari/handle_error.pyx"],
              libraries = ['pari', 'gmp']),

    Extension('sage.libs.pari.pari_instance',
              sources = ["sage/libs/pari/pari_instance.pyx"],
              extra_compile_args = ["-std=c99", "-D_XPG6"],
              libraries = [ 'pari', 'gmp', 'flint'],
              depends = flint_depends),

    Extension('sage.libs.ppl',
              sources = ['sage/libs/ppl.pyx', 'sage/libs/ppl_shim.cc'],
              libraries = ['ppl', 'gmpxx', 'gmp', 'm'],
              language="c++",
              depends = [SAGE_INC + "/ppl.hh"]),

    Extension('sage.libs.ratpoints',
              sources = ["sage/libs/ratpoints.pyx"],
              depends = [SAGE_INC + '/ratpoints.h'],
              libraries = ["ratpoints", "gmp"]),

    Extension('sage.libs.readline',
              sources = ['sage/libs/readline.pyx'],
              libraries = ['readline']),

    Extension('sage.libs.singular.singular',
              sources = ['sage/libs/singular/singular.pyx'],
              libraries = ['givaro'] + singular_libs,
              language="c++",
              include_dirs = singular_incs,
              depends = singular_depends + givaro_depends,
              extra_compile_args = givaro_extra_compile_args),

    Extension('sage.libs.singular.polynomial',
              sources = ['sage/libs/singular/polynomial.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs,
              depends = singular_depends),

    Extension('sage.libs.singular.ring',
              sources = ['sage/libs/singular/ring.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs,
              depends = singular_depends),

    Extension('sage.libs.singular.groebner_strategy',
              sources = ['sage/libs/singular/groebner_strategy.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs,
              depends = singular_depends),

    Extension('sage.libs.singular.function',
              sources = ['sage/libs/singular/function.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs,
              depends = singular_depends,
              extra_compile_args = givaro_extra_compile_args),

    Extension('sage.libs.singular.option',
              sources = ['sage/libs/singular/option.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs,
              depends = singular_depends),

    Extension('sage.libs.symmetrica.symmetrica',
              sources = ["sage/libs/symmetrica/%s"%s for s in ["symmetrica.pyx"]],
              include_dirs = ['/usr/include/malloc/'],
              libraries = ["symmetrica"],
              depends = [SAGE_INC + "/symmetrica/def.h"]),

    Extension('sage.libs.mpmath.utils',
              sources = ["sage/libs/mpmath/utils.pyx"],
              libraries = ['mpfr', 'gmp']),

    Extension('sage.libs.mpmath.ext_impl',
              sources = ["sage/libs/mpmath/ext_impl.pyx"],
              libraries = ['mpfr', 'gmp']),

    Extension('sage.libs.mpmath.ext_main',
              sources = ["sage/libs/mpmath/ext_main.pyx"],
              libraries = ['gmp']),

    Extension('sage.libs.mpmath.ext_libmp',
              sources = ["sage/libs/mpmath/ext_libmp.pyx"],
              libraries = ['gmp']),

    ################################
    ##
    ## sage.libs.gap
    ##
    ################################

    Extension('sage.libs.gap.util',
              sources = ["sage/libs/gap/util.pyx"],
              libraries = ['gmp', 'gap', 'm']),

    Extension('sage.libs.gap.element',
              sources = ["sage/libs/gap/element.pyx"],
              libraries = ['gmp', 'gap', 'm']),

    Extension('sage.libs.gap.libgap',
              sources = ["sage/libs/gap/libgap.pyx"],
              libraries = ['gmp', 'gap', 'm']),

    ###################################
    ##
    ## sage.libs.cremona
    ##
    ###################################

    Extension('sage.libs.cremona.homspace',
              sources = ["sage/libs/cremona/homspace.pyx"],
              libraries = ['ec', 'ntl', 'pari',
                           'gmpxx', 'gmp', 'm'],
              language='c++',
              define_macros = [("NTL_ALL",None)],
              depends = [ SAGE_INC + "/eclib/" + h for h in
                          ["interface.h","bigrat.h","rat.h","curve.h",
                           "moddata.h","symb.h","cusp.h","homspace.h","mat.h"]
                        ]),

    Extension('sage.libs.cremona.mat',
              sources = ["sage/libs/cremona/mat.pyx"],
              libraries = ['ec', 'ntl', 'pari',
                           'gmpxx', 'gmp', 'm'],
              language='c++',
              define_macros = [("NTL_ALL",None)],
              depends = [ SAGE_INC + "/eclib/" + h for h in
                          ["interface.h","bigrat.h","rat.h","curve.h",
                           "moddata.h","symb.h","cusp.h","homspace.h","mat.h"]
                        ]),

    Extension('sage.libs.cremona.newforms',
              sources = ["sage/libs/cremona/newforms.pyx"],
              libraries = ['ec', 'ntl', 'pari',
                           'gmpxx', 'gmp', 'm'],
              language='c++',
              define_macros = [("NTL_ALL",None)],
              depends = [ SAGE_INC + "/eclib/" + h for h in
                          ["interface.h","bigrat.h","rat.h","curve.h",
                           "moddata.h","symb.h","cusp.h","xsplit.h","method.h",
                           "oldforms.h","homspace.h","cperiods.h","newforms.h"]
                        ]),

    ###################################
    ##
    ## sage.libs.ntl
    ##
    ###################################

    Extension('sage.libs.ntl.error',
              sources = ["sage/libs/ntl/error.pyx"],
              libraries = ["ntl", "gmp", "gmpxx"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2',
              sources = ["sage/libs/ntl/ntl_GF2.pyx"],
              libraries = ["ntl", "gmp", "gmpxx"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2E',
              sources = ["sage/libs/ntl/ntl_GF2E.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2EContext',
              sources = ["sage/libs/ntl/ntl_GF2EContext.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2EX',
              sources = ["sage/libs/ntl/ntl_GF2EX.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_GF2X',
              sources = ["sage/libs/ntl/ntl_GF2X.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_lzz_p',
              sources = ["sage/libs/ntl/ntl_lzz_p.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_lzz_pContext',
              sources = ["sage/libs/ntl/ntl_lzz_pContext.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_lzz_pX',
              sources = ["sage/libs/ntl/ntl_lzz_pX.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_mat_GF2',
              sources = ["sage/libs/ntl/ntl_mat_GF2.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_mat_GF2E',
              sources = ["sage/libs/ntl/ntl_mat_GF2E.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_mat_ZZ',
              sources = ["sage/libs/ntl/ntl_mat_ZZ.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ',
              sources = ["sage/libs/ntl/ntl_ZZ.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZX',
              sources = ["sage/libs/ntl/ntl_ZZX.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_p',
              sources = ["sage/libs/ntl/ntl_ZZ_p.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pContext',
              sources = ["sage/libs/ntl/ntl_ZZ_pContext.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pE',
              sources = ["sage/libs/ntl/ntl_ZZ_pE.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pEContext',
              sources = ["sage/libs/ntl/ntl_ZZ_pEContext.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pEX',
              sources = ["sage/libs/ntl/ntl_ZZ_pEX.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.libs.ntl.ntl_ZZ_pX',
              sources = ["sage/libs/ntl/ntl_ZZ_pX.pyx"],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    ################################
    ##
    ## sage.matrix
    ##
    ################################

    Extension('sage.matrix.action',
              sources = ['sage/matrix/action.pyx']),

    Extension('sage.matrix.echelon_matrix',
              sources = ['sage/matrix/echelon_matrix.pyx']),

    Extension('sage.matrix.change_ring',
              sources = ['sage/matrix/change_ring.pyx'],
              libraries=[BLAS, BLAS2, 'gmp']),

    Extension('sage.matrix.matrix',
              sources = ['sage/matrix/matrix.pyx']),

    Extension('sage.matrix.matrix0',
              sources = ['sage/matrix/matrix0.pyx']),

    Extension('sage.matrix.matrix1',
              sources = ['sage/matrix/matrix1.pyx']),

    Extension('sage.matrix.matrix2',
              sources = ['sage/matrix/matrix2.pyx']),

    Extension('sage.matrix.matrix_complex_double_dense',
              sources = ['sage/matrix/matrix_complex_double_dense.pyx'],
              libraries=[BLAS, BLAS2]),

    Extension('sage.matrix.matrix_cyclo_dense',
              sources = ['sage/matrix/matrix_cyclo_dense.pyx'],
              language = "c++",
              libraries=['ntl', 'gmp']),

    Extension('sage.matrix.matrix_dense',
              sources = ['sage/matrix/matrix_dense.pyx']),

    Extension('sage.matrix.matrix_double_dense',
              sources = ['sage/matrix/matrix_double_dense.pyx'],
              libraries=[BLAS, BLAS2]),

    Extension('sage.matrix.matrix_generic_dense',
              sources = ['sage/matrix/matrix_generic_dense.pyx']),

    Extension('sage.matrix.matrix_generic_sparse',
              sources = ['sage/matrix/matrix_generic_sparse.pyx']),

    Extension('sage.matrix.matrix_integer_dense',
              sources = ['sage/matrix/matrix_integer_dense.pyx'],
              extra_compile_args = ['-std=c99'] + m4ri_extra_compile_args,
              # order matters for cygwin!!
              libraries = ['iml', 'pari', 'ntl', 'gmp', 'm', 'flint', BLAS, BLAS2],
              depends = [SAGE_INC + '/m4ri/m4ri.h'] + flint_depends),

    Extension('sage.matrix.matrix_integer_sparse',
              sources = ['sage/matrix/matrix_integer_sparse.pyx'],
              libraries = ['gmp']),

    Extension('sage.matrix.matrix_mod2_dense',
              sources = ['sage/matrix/matrix_mod2_dense.pyx'],
              libraries = ['gmp','m4ri', 'gd', 'png12', 'z'],
              extra_compile_args = m4ri_extra_compile_args,
              depends = [SAGE_INC + "/png.h", SAGE_INC + "/m4ri/m4ri.h"]),

    Extension('sage.matrix.matrix_mod2e_dense',
              sources = ['sage/matrix/matrix_mod2e_dense.pyx'],
              libraries = ['m4rie', 'm4ri', 'm'],
              depends = [SAGE_INC + "/m4rie/m4rie.h"],
              include_dirs = [SAGE_INC + '/m4rie'],
              extra_compile_args = m4ri_extra_compile_args),

    Extension('sage.matrix.matrix_modn_dense_float',
              sources = ['sage/matrix/matrix_modn_dense_float.pyx'],
              language="c++",
              libraries = ['linbox', 'givaro', 'mpfr', 'gmpxx', 'gmp', BLAS, BLAS2],
              extra_compile_args = ['-DDISABLE_COMMENTATOR'] + givaro_extra_compile_args),

    Extension('sage.matrix.matrix_modn_dense_double',
              sources = ['sage/matrix/matrix_modn_dense_double.pyx'],
              language="c++",
              libraries = ['linbox', 'givaro', 'mpfr', 'gmpxx', 'gmp', BLAS, BLAS2],
              extra_compile_args = ["-D_XPG6", "-DDISABLE_COMMENTATOR"]
                    + m4ri_extra_compile_args + givaro_extra_compile_args),

    Extension('sage.matrix.matrix_modn_sparse',
              sources = ['sage/matrix/matrix_modn_sparse.pyx'],
              libraries = ['gmp']),

    Extension('sage.matrix.matrix_mpolynomial_dense',
              sources = ['sage/matrix/matrix_mpolynomial_dense.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs,
              depends = singular_depends),

    Extension('sage.matrix.matrix_rational_dense',
              sources = ['sage/matrix/matrix_rational_dense.pyx'],
              extra_compile_args = ["-std=c99", "-D_XPG6"] + m4ri_extra_compile_args,
              # order matters for cygwin!!
              libraries = ['iml', 'pari', 'ntl', 'gmp', 'm', 'flint', BLAS, BLAS2],
              depends = [SAGE_INC + '/m4ri/m4ri.h'] + flint_depends),

    Extension('sage.matrix.matrix_rational_sparse',
              sources = ['sage/matrix/matrix_rational_sparse.pyx'],
              libraries = ['gmp']),

    Extension('sage.matrix.matrix_real_double_dense',
              sources = ['sage/matrix/matrix_real_double_dense.pyx'],
              libraries=[BLAS, BLAS2]),

    Extension('sage.matrix.matrix_sparse',
              sources = ['sage/matrix/matrix_sparse.pyx']),

    Extension('sage.matrix.matrix_symbolic_dense',
              sources = ['sage/matrix/matrix_symbolic_dense.pyx']),

    Extension('sage.matrix.matrix_window',
              sources = ['sage/matrix/matrix_window.pyx']),

    Extension('sage.matrix.misc',
              sources = ['sage/matrix/misc.pyx'],
              libraries=['mpfr','gmp']),

    Extension('sage.matrix.strassen',
              sources = ['sage/matrix/strassen.pyx']),

    ################################
    ##
    ## sage.matroids
    ##
    ################################

    Extension('sage.matroids.matroid',
            ['sage/matroids/matroid.pyx']),

    Extension('sage.matroids.extension',
            ['sage/matroids/extension.pyx'],
              libraries = ['gmp']),

    Extension('sage.matroids.set_system',
            ['sage/matroids/set_system.pyx'],
              libraries = ['gmp']),

    Extension('sage.matroids.lean_matrix',
            ['sage/matroids/lean_matrix.pyx'],
              libraries = ['gmp']),

    Extension('sage.matroids.basis_exchange_matroid',
            ['sage/matroids/basis_exchange_matroid.pyx'],
              libraries = ['gmp']),

    Extension('sage.matroids.basis_matroid',
            ['sage/matroids/basis_matroid.pyx'],
              libraries = ['gmp']),

    Extension('sage.matroids.linear_matroid',
            ['sage/matroids/linear_matroid.pyx'],
              libraries = ['gmp']),

    Extension('sage.matroids.circuit_closures_matroid',
            ['sage/matroids/circuit_closures_matroid.pyx']),

    Extension('sage.matroids.unpickling',
            ['sage/matroids/unpickling.pyx'],
              libraries = ['gmp']),

    ################################
    ##
    ## sage.media
    ##
    ################################

    Extension('sage.media.channels',
              sources = ['sage/media/channels.pyx']),

    ################################
    ##
    ## sage.misc
    ##
    ################################

    Extension('*', ['sage/misc/*.pyx']),

    # Only include darwin_utilities on OS_X >= 10.5
    OptionalExtension('sage.misc.darwin_utilities',
        sources = ['sage/misc/darwin_memory_usage.c',
                   'sage/misc/darwin_utilities.pyx'],
        depends = ['sage/misc/darwin_memory_usage.h'],
        condition = (UNAME[0] == "Darwin" and not UNAME[2].startswith('8.'))),

    ################################
    ##
    ## sage.modular
    ##
    ################################

    Extension('sage.modular.arithgroup.congroup',
              sources = ['sage/modular/arithgroup/congroup.pyx']),

    Extension('sage.modular.arithgroup.farey_symbol',
              sources = ['sage/modular/arithgroup/farey_symbol.pyx',
                         'sage/modular/arithgroup/farey.cpp',
                         'sage/modular/arithgroup/sl2z.cpp'],
              libraries = ['gmpxx', 'gmp'],
              language = 'c++'),

    Extension('sage.modular.arithgroup.arithgroup_element',
              sources = ['sage/modular/arithgroup/arithgroup_element.pyx']),

    Extension('sage.modular.modform.eis_series_cython',
              sources = ['sage/modular/modform/eis_series_cython.pyx'],
              libraries = ["gmp", "flint"],
              extra_compile_args = ['-std=c99'],
              depends = flint_depends),

    Extension('sage.modular.modsym.apply',
              sources = ['sage/modular/modsym/apply.pyx'],
              libraries = ["flint", "gmp", "gmpxx", "m", "stdc++"],
              extra_compile_args=["-std=c99", "-D_XPG6"],
              depends = flint_depends),

    Extension('sage.modular.modsym.manin_symbol',
              sources = ['sage/modular/modsym/manin_symbol.pyx']),

    Extension('sage.modular.modsym.relation_matrix_pyx',
              sources = ['sage/modular/modsym/relation_matrix_pyx.pyx']),

    Extension('sage.modular.modsym.heilbronn',
              sources = ['sage/modular/modsym/heilbronn.pyx'],
              libraries = ["flint", "gmp", "gmpxx", "m", "stdc++"],
              extra_compile_args=["-std=c99", "-D_XPG6"],
              depends = flint_depends),

    Extension('sage.modular.modsym.p1list',
              sources = ['sage/modular/modsym/p1list.pyx'],
              libraries = ['gmp']),

    ################################
    ##
    ## sage.modules
    ##
    ################################

    Extension('sage.modules.finite_submodule_iter',
              sources = ['sage/modules/finite_submodule_iter.pyx']),

    Extension('sage.modules.free_module_element',
              sources = ['sage/modules/free_module_element.pyx']),

    Extension('sage.modules.module',
              sources = ['sage/modules/module.pyx']),

    Extension('sage.modules.vector_complex_double_dense',
              ['sage/modules/vector_complex_double_dense.pyx'],
              libraries = [BLAS, BLAS2]),

    Extension('sage.modules.vector_double_dense',
              ['sage/modules/vector_double_dense.pyx'],
              libraries = [BLAS, BLAS2]),

    Extension('sage.modules.vector_integer_dense',
              sources = ['sage/modules/vector_integer_dense.pyx'],
              libraries = ['gmp']),

    Extension('sage.modules.vector_modn_dense',
              extra_compile_args = ['-std=c99'],
              sources = ['sage/modules/vector_modn_dense.pyx']),

    Extension('sage.modules.vector_mod2_dense',
              sources = ['sage/modules/vector_mod2_dense.pyx'],
              libraries = ['gmp','m4ri', 'png12', 'gd'],
              extra_compile_args = m4ri_extra_compile_args,
              depends = [SAGE_INC + "/png.h", SAGE_INC + "/m4ri/m4ri.h"]),

    Extension('sage.modules.vector_rational_dense',
              sources = ['sage/modules/vector_rational_dense.pyx'],
              libraries = ['gmp']),

    Extension('sage.modules.vector_real_double_dense',
              ['sage/modules/vector_real_double_dense.pyx'],
              libraries = [BLAS, BLAS2]),

    # Extension('sage.modules.vector_rational_sparse',
    #           sources = ['sage/modules/vector_rational_sparse.pyx'],
    #           libraries = ['gmp']),

    ################################
    ##
    ## sage.numerical
    ##
    ################################


    Extension("sage.numerical.mip",
              ["sage/numerical/mip.pyx"],
              libraries=["stdc++"]),

    Extension("sage.numerical.linear_functions",
              ["sage/numerical/linear_functions.pyx"],
              libraries=["stdc++"]),

    Extension("sage.numerical.linear_tensor_element",
              ["sage/numerical/linear_tensor_element.pyx"],
              include_dirs=[SAGE_INC],
              libraries=["stdc++"]),

    Extension("sage.numerical.backends.generic_backend",
              ["sage/numerical/backends/generic_backend.pyx"],
              libraries=["stdc++"]),

    Extension("sage.numerical.backends.glpk_backend",
              ["sage/numerical/backends/glpk_backend.pyx"],
              language = 'c++',
              libraries=["stdc++", "glpk", "gmp", "z"]),

    Extension("sage.numerical.backends.ppl_backend",
              ["sage/numerical/backends/ppl_backend.pyx"],
              libraries=["stdc++"]),

    Extension("sage.numerical.backends.cvxopt_backend",
              ["sage/numerical/backends/cvxopt_backend.pyx"],
              libraries=["stdc++"]),

    Extension("sage.numerical.backends.glpk_graph_backend",
              ["sage/numerical/backends/glpk_graph_backend.pyx"],
              language = 'c++',
              libraries=["stdc++", "glpk", "gmp", "z"]),

    OptionalExtension("sage.numerical.backends.gurobi_backend",
              ["sage/numerical/backends/gurobi_backend.pyx"],
              libraries = ["stdc++", "gurobi"],
              condition = os.path.isfile(SAGE_INC + "/gurobi_c.h") and
                  os.path.isfile(SAGE_LOCAL + "/lib/libgurobi.so")),

    OptionalExtension("sage.numerical.backends.cplex_backend",
              ["sage/numerical/backends/cplex_backend.pyx"],
              libraries = ["stdc++", "cplex"],
              condition = os.path.isfile(SAGE_INC + "/cplex.h") and
                  os.path.isfile(SAGE_LOCAL + "/lib/libcplex.a")),

    OptionalExtension("sage.numerical.backends.coin_backend",
              ["sage/numerical/backends/coin_backend.pyx"],
              language = 'c++',
              libraries = ["Cbc", "CbcSolver", "Cgl", "Clp", "CoinUtils", "OsiCbc", "OsiClp", "Osi", "lapack"],
              package = 'cbc'),

    ################################
    ##
    ## sage.plot
    ##
    ################################

    Extension('sage.plot.complex_plot',
              sources = ['sage/plot/complex_plot.pyx']),

    Extension('sage.plot.plot3d.base',
              sources = ['sage/plot/plot3d/base.pyx'],
              extra_compile_args=["-std=c99"]),

    Extension('sage.plot.plot3d.implicit_surface',
              sources = ['sage/plot/plot3d/implicit_surface.pyx'],
              libraries = ['gsl']),

    Extension('sage.plot.plot3d.index_face_set',
              sources = ['sage/plot/plot3d/index_face_set.pyx'],
              extra_compile_args=["-std=c99"]),

    Extension('sage.plot.plot3d.parametric_surface',
              sources = ['sage/plot/plot3d/parametric_surface.pyx']),

    Extension('sage.plot.plot3d.shapes',
              sources = ['sage/plot/plot3d/shapes.pyx']),

    Extension('sage.plot.plot3d.transform',
              sources = ['sage/plot/plot3d/transform.pyx']),

    ################################
    ##
    ## sage.quadratic_forms
    ##
    ################################

    Extension('sage.quadratic_forms.count_local_2',
              sources = ['sage/quadratic_forms/count_local_2.pyx'],
              libraries = ['gmp']),

    Extension('sage.quadratic_forms.quadratic_form__evaluate',
              sources = ['sage/quadratic_forms/quadratic_form__evaluate.pyx']),


    Extension('sage.quadratic_forms.ternary',
              sources = ['sage/quadratic_forms/ternary.pyx']),

    ################################
    ##
    ## sage.repl
    ##
    ################################

    Extension('sage.repl.inputhook',
              sources = ['sage/repl/inputhook.pyx']),

    Extension('sage.repl.readline_extra_commands',
              sources = ['sage/repl/readline_extra_commands.pyx'],
              libraries = ['readline']),

    ################################
    ##
    ## sage.rings
    ##
    ################################

    Extension('sage.rings.sum_of_squares',
              sources = ['sage/rings/sum_of_squares.pyx'],
              libraries = ['m']),

    Extension('sage.rings.bernmm',
              sources = ['sage/rings/bernmm.pyx',
                         'sage/rings/bernmm/bern_modp.cpp',
                         'sage/rings/bernmm/bern_modp_util.cpp',
                         'sage/rings/bernmm/bern_rat.cpp'],
              libraries = ['ntl', 'gmp', 'stdc++', 'pthread'],
              depends = ['sage/rings/bernmm/bern_modp.h',
                         'sage/rings/bernmm/bern_modp_util.h',
                         'sage/rings/bernmm/bern_rat.h'],
              language = 'c++',
              define_macros=[('USE_THREADS', '1'),
                             ('THREAD_STACK_SIZE', '4096')]),

    Extension('sage.rings.bernoulli_mod_p',
              sources = ['sage/rings/bernoulli_mod_p.pyx'],
              libraries=['ntl'],
              language = 'c++',
              include_dirs = ['sage/libs/ntl/']),

    Extension('sage.rings.complex_double',
              sources = ['sage/rings/complex_double.pyx'],
              extra_compile_args=["-std=c99", "-D_XPG6"],
              libraries = (['gsl', BLAS, BLAS2, 'pari', 'gmp', 'm'])),

    Extension('sage.rings.complex_interval',
              sources = ['sage/rings/complex_interval.pyx'],
              libraries = ['mpfi', 'mpfr', 'gmp']),

    Extension('sage.rings.complex_number',
              sources = ['sage/rings/complex_number.pyx'],
              libraries = ['mpfr', 'gmp']),

    Extension('sage.rings.integer',
              sources = ['sage/rings/integer.pyx'],
              libraries=['ntl', 'pari', 'flint', 'gmp'],
              depends = flint_depends),

    Extension('sage.rings.integer_ring',
              sources = ['sage/rings/integer_ring.pyx'],
              libraries=['ntl', 'gmp']),

    Extension('sage.rings.factorint',
              sources = ['sage/rings/factorint.pyx'],
              libraries=['gmp']),

    Extension('sage.rings.fast_arith',
              sources = ['sage/rings/fast_arith.pyx'],
              libraries=['pari','gmp']),

    Extension('sage.rings.fraction_field_element',
              sources = ['sage/rings/fraction_field_element.pyx']),

    Extension('sage.rings.fraction_field_FpT',
              sources = ['sage/rings/fraction_field_FpT.pyx'],
              libraries = ["flint", "gmp", "gmpxx", "ntl", "zn_poly"],
              language = 'c++',
              depends = flint_depends),

    Extension('sage.rings.laurent_series_ring_element',
              sources = ['sage/rings/laurent_series_ring_element.pyx']),

    Extension('sage.rings.morphism',
              sources = ['sage/rings/morphism.pyx']),

    Extension('sage.rings.complex_mpc',
              sources = ['sage/rings/complex_mpc.pyx'],
              libraries = ['mpc', 'mpfr', 'gmp']),

    Extension('sage.rings.noncommutative_ideals',
              sources = ['sage/rings/noncommutative_ideals.pyx']),

    Extension('sage.rings.power_series_mpoly',
              sources = ['sage/rings/power_series_mpoly.pyx']),

    Extension('sage.rings.power_series_poly',
              sources = ['sage/rings/power_series_poly.pyx']),

    Extension('sage.rings.power_series_ring_element',
              sources = ['sage/rings/power_series_ring_element.pyx']),

    Extension('sage.rings.rational',
              sources = ['sage/rings/rational.pyx'],
              libraries=['ntl', 'gmp']),

    Extension('sage.rings.real_double',
              sources = ['sage/rings/real_double.pyx'],
              libraries = ['gsl', 'gmp', BLAS, BLAS2],
              define_macros=[('GSL_DISABLE_DEPRECATED','1')]),

    Extension('sage.rings.real_interval_absolute',
              sources = ['sage/rings/real_interval_absolute.pyx'],
              libraries = ['gmp']),

    OptionalExtension("sage.rings.real_arb",
                      ["sage/rings/real_arb.pyx"],
                      libraries = ['arb', 'mpfi', 'mpfr'],
                      include_dirs = [SAGE_INC + '/flint'],
                      depends = flint_depends,
                      package = 'arb'),

    Extension('sage.rings.real_lazy',
              sources = ['sage/rings/real_lazy.pyx']),

    Extension('sage.rings.real_mpfi',
              sources = ['sage/rings/real_mpfi.pyx'],
              libraries = ['mpfi', 'mpfr', 'gmp']),

    Extension('sage.rings.real_mpfr',
              sources = ['sage/rings/real_mpfr.pyx'],
              libraries = ['mpfr', 'pari', 'gmp']),

    Extension('sage.rings.finite_rings.residue_field',
              sources = ['sage/rings/finite_rings/residue_field.pyx']),

    Extension('sage.rings.ring',
              sources = ['sage/rings/ring.pyx']),

    ################################
    ##
    ## sage.rings.finite_rings
    ##
    ################################

    Extension('sage.rings.finite_rings.finite_field_base',
              sources = ['sage/rings/finite_rings/finite_field_base.pyx']),

    Extension('sage.rings.finite_rings.element_base',
              sources = ['sage/rings/finite_rings/element_base.pyx']),

    Extension('sage.rings.finite_rings.integer_mod',
              sources = ['sage/rings/finite_rings/integer_mod.pyx'],
              libraries = ['gmp']),

    Extension('sage.rings.finite_rings.element_givaro',
              sources = ["sage/rings/finite_rings/element_givaro.pyx"],
              # this order is needed to compile under windows.
              libraries = ['givaro', 'ntl', 'pari', 'gmpxx', 'gmp', 'm'],
              language='c++',
              extra_compile_args = givaro_extra_compile_args),

    Extension('sage.rings.finite_rings.element_ntl_gf2e',
              sources = ['sage/rings/finite_rings/element_ntl_gf2e.pyx'],
              libraries = ['ntl', 'pari', 'gmp'],
              language = 'c++'),

    Extension('sage.rings.finite_rings.element_pari_ffelt',
              sources = ['sage/rings/finite_rings/element_pari_ffelt.pyx'],
              libraries = ['pari', 'gmp']),

    Extension('sage.rings.finite_rings.hom_finite_field',
              sources = ["sage/rings/finite_rings/hom_finite_field.pyx"]),

    Extension('sage.rings.finite_rings.hom_prime_finite_field',
              sources = ["sage/rings/finite_rings/hom_prime_finite_field.pyx"]),

    Extension('sage.rings.finite_rings.hom_finite_field_givaro',
              sources = ["sage/rings/finite_rings/hom_finite_field_givaro.pyx"],
              # this order is needed to compile under windows.
              libraries = ['givaro', 'ntl', 'gmpxx', 'gmp', 'm'],
              language='c++',
              extra_compile_args = givaro_extra_compile_args),

    ################################
    ##
    ## sage.rings.function_field
    ##
    ################################

    Extension('sage.rings.function_field.function_field_element',
              sources = ['sage/rings/function_field/function_field_element.pyx']),

    ################################
    ##
    ## sage.rings.number_field
    ##
    ################################

    Extension('sage.rings.number_field.number_field_base',
              sources = ['sage/rings/number_field/number_field_base.pyx']),

    Extension('sage.rings.number_field.number_field_element',
              sources = ['sage/rings/number_field/number_field_element.pyx'],
              libraries=['ntl','gmp'],
              language = 'c++'),

    Extension('sage.rings.number_field.number_field_element_quadratic',
              sources = ['sage/rings/number_field/number_field_element_quadratic.pyx'],
              libraries=['ntl', 'gmp'],
              language = 'c++'),

    Extension('sage.rings.number_field.number_field_morphisms',
              sources = ['sage/rings/number_field/number_field_morphisms.pyx']),

    Extension('sage.rings.number_field.totallyreal',
              sources = ['sage/rings/number_field/totallyreal.pyx'],
              libraries = ['pari', 'gmp']),

    Extension('sage.rings.number_field.totallyreal_data',
              sources = ['sage/rings/number_field/totallyreal_data.pyx'],
              libraries = ['gmp']),

    ################################
    ##
    ## sage.rings.padics
    ##
    ################################

    Extension('sage.rings.padics.morphism',
              sources = ['sage/rings/padics/morphism.pyx']),

    Extension('sage.rings.padics.common_conversion',
              sources = ['sage/rings/padics/common_conversion.pyx'],
              libraries=['gmp']),

    Extension('sage.rings.padics.local_generic_element',
              sources = ['sage/rings/padics/local_generic_element.pyx']),

    Extension('sage.rings.padics.padic_capped_absolute_element',
              sources = ['sage/rings/padics/padic_capped_absolute_element.pyx'],
              libraries=['gmp']),

    Extension('sage.rings.padics.padic_capped_relative_element',
              sources = ['sage/rings/padics/padic_capped_relative_element.pyx'],
              libraries=['gmp']),

    Extension('sage.rings.padics.padic_ext_element',
              sources = ['sage/rings/padics/padic_ext_element.pyx'],
              libraries=['ntl', 'gmp', 'gmpxx', 'm'],
              language='c++'),

    Extension('sage.rings.padics.padic_fixed_mod_element',
              sources = ['sage/rings/padics/padic_fixed_mod_element.pyx'],
              libraries=['gmp']),

    Extension('sage.rings.padics.padic_generic_element',
              sources = ['sage/rings/padics/padic_generic_element.pyx'],
              libraries=['gmp']),

    Extension('sage.rings.padics.padic_printing',
              sources = ['sage/rings/padics/padic_printing.pyx'],
              libraries=['gmp', 'ntl', 'gmpxx', 'm'],
              language='c++'),

    Extension('sage.rings.padics.padic_ZZ_pX_CA_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_CA_element.pyx'],
              libraries = ['ntl', 'gmp', 'gmpxx','m'],
              language='c++'),

    Extension('sage.rings.padics.padic_ZZ_pX_CR_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_CR_element.pyx'],
              libraries=['ntl', 'gmp', 'gmpxx','m'],
              language='c++'),

    Extension('sage.rings.padics.padic_ZZ_pX_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_element.pyx'],
              libraries=['ntl', 'gmp', 'gmpxx', 'm'],
              language='c++'),

    Extension('sage.rings.padics.padic_ZZ_pX_FM_element',
              sources = ['sage/rings/padics/padic_ZZ_pX_FM_element.pyx'],
              libraries=['ntl', 'gmp', 'gmpxx', 'm'],
              language='c++'),

    Extension('sage.rings.padics.pow_computer',
              sources = ['sage/rings/padics/pow_computer.pyx'],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    Extension('sage.rings.padics.pow_computer_ext',
              sources = ['sage/rings/padics/pow_computer_ext.pyx'],
              libraries = ["ntl", "gmp", "gmpxx", "m"],
              language='c++'),

    ################################
    ##
    ## sage.rings.polynomial
    ##
    ################################

    Extension('sage.rings.polynomial.cyclotomic',
              sources = ['sage/rings/polynomial/cyclotomic.pyx']),

    Extension('sage.rings.polynomial.laurent_polynomial',
              sources = ['sage/rings/polynomial/laurent_polynomial.pyx']),

    Extension('sage.rings.polynomial.multi_polynomial',
              sources = ['sage/rings/polynomial/multi_polynomial.pyx']),

    Extension('sage.rings.polynomial.multi_polynomial_ideal_libsingular',
              sources = ['sage/rings/polynomial/multi_polynomial_ideal_libsingular.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs,
              depends = singular_depends),

    Extension('sage.rings.polynomial.plural',
              sources = ['sage/rings/polynomial/plural.pyx'],
              libraries = ['m', 'readline', 'singular', 'givaro', 'gmpxx', 'gmp'],
              language="c++",
              include_dirs = singular_incs,
              depends = [SAGE_INC + "/libsingular.h"],
              extra_compile_args = givaro_extra_compile_args),

    Extension('sage.rings.polynomial.multi_polynomial_libsingular',
              sources = ['sage/rings/polynomial/multi_polynomial_libsingular.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs,
              depends = singular_depends),

    Extension('sage.rings.polynomial.multi_polynomial_ring_generic',
              sources = ['sage/rings/polynomial/multi_polynomial_ring_generic.pyx']),

    Extension('sage.rings.polynomial.polynomial_number_field',
              sources = ['sage/rings/polynomial/polynomial_number_field.pyx']),

    Extension('sage.rings.polynomial.polydict',
              sources = ['sage/rings/polynomial/polydict.pyx']),

    Extension('sage.rings.polynomial.polynomial_compiled',
               sources = ['sage/rings/polynomial/polynomial_compiled.pyx']),

    Extension('sage.rings.polynomial.polynomial_element',
              sources = ['sage/rings/polynomial/polynomial_element.pyx']),

    Extension('sage.rings.polynomial.polynomial_gf2x',
              sources = ['sage/rings/polynomial/polynomial_gf2x.pyx'],
              libraries = ['ntl', 'stdc++', 'gmp'],
              extra_compile_args = m4ri_extra_compile_args,
              language = 'c++',
              depends = [SAGE_INC + '/m4ri/m4ri.h'],
              include_dirs = ['sage/libs/ntl/']),

    Extension('sage.rings.polynomial.polynomial_zz_pex',
              sources = ['sage/rings/polynomial/polynomial_zz_pex.pyx'],
              libraries = ['ntl', 'stdc++', 'gmp'],
              language = 'c++',
              include_dirs = ['sage/libs/ntl/']),

    Extension('sage.rings.polynomial.polynomial_zmod_flint',
              sources = ['sage/rings/polynomial/polynomial_zmod_flint.pyx'],
              libraries = ["flint", "gmp", "gmpxx", "ntl", "zn_poly"],
              language = 'c++',
              depends = flint_depends),

    Extension('sage.rings.polynomial.polynomial_integer_dense_flint',
              sources = ['sage/rings/polynomial/polynomial_integer_dense_flint.pyx'],
              language = 'c++',
              libraries = ["flint", "ntl", "gmpxx", "gmp"],
              depends = flint_depends),

    Extension('sage.rings.polynomial.polynomial_integer_dense_ntl',
              sources = ['sage/rings/polynomial/polynomial_integer_dense_ntl.pyx'],
              libraries = ['ntl', 'stdc++', 'gmp'],
              language = 'c++',
              include_dirs = ['sage/libs/ntl/']),

    Extension('sage.rings.polynomial.polynomial_rational_flint',
              sources = ['sage/rings/polynomial/polynomial_rational_flint.pyx'],
              libraries = ["flint", "ntl", "gmpxx", "gmp"],
              language = 'c++',
              depends = flint_depends),

    Extension('sage.rings.polynomial.polynomial_modn_dense_ntl',
              sources = ['sage/rings/polynomial/polynomial_modn_dense_ntl.pyx'],
              libraries = ['ntl', 'stdc++', 'gmp'],
              language = 'c++',
              include_dirs = ['sage/libs/ntl/']),

    Extension('sage.rings.polynomial.polynomial_ring_homomorphism',
              sources = ['sage/rings/polynomial/polynomial_ring_homomorphism.pyx']),

    Extension('sage.rings.polynomial.pbori',
              sources = ['sage/rings/polynomial/pbori.pyx'],
              libraries=['polybori-' + polybori_major_version,
                         'polybori_groebner-' + polybori_major_version, 'm4ri', 'gd', 'png12'],
              include_dirs = [SAGE_INC, "sage/libs/polybori"],
              depends = [SAGE_INC + "/polybori/" + hd + ".h" for hd in ["polybori", "config"] ] + \
                        [SAGE_INC + '/m4ri/m4ri.h'],
              extra_compile_args = polybori_extra_compile_args + m4ri_extra_compile_args,
              language = 'c++'),

    Extension('sage.rings.polynomial.polynomial_real_mpfr_dense',
              sources = ['sage/rings/polynomial/polynomial_real_mpfr_dense.pyx'],
              libraries = ['mpfr', 'gmp']),

    Extension('sage.rings.polynomial.real_roots',
              sources = ['sage/rings/polynomial/real_roots.pyx'],
              libraries=['mpfr', 'gmp']),

    Extension('sage.rings.polynomial.symmetric_reduction',
              sources = ['sage/rings/polynomial/symmetric_reduction.pyx']),

    ################################
    ##
    ## sage.rings.semirings
    ##
    ################################

    Extension('sage.rings.semirings.tropical_semiring',
              sources = ['sage/rings/semirings/tropical_semiring.pyx']),

    ################################
    ##
    ## sage.schemes
    ##
    ################################

    Extension('sage.schemes.elliptic_curves.descent_two_isogeny',
              sources = ['sage/schemes/elliptic_curves/descent_two_isogeny.pyx'],
              extra_compile_args=["-std=c99"],
              depends = [SAGE_INC + '/ratpoints.h',
                         SAGE_INC + '/gmp.h'] +
                         flint_depends,
              libraries = ['flint', 'gmp', 'ratpoints']),

    Extension('sage.schemes.elliptic_curves.period_lattice_region',
              sources = ['sage/schemes/elliptic_curves/period_lattice_region.pyx']),

    Extension('sage.schemes.hyperelliptic_curves.hypellfrob',
              sources = ['sage/schemes/hyperelliptic_curves/hypellfrob.pyx',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/hypellfrob.cpp',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_ntl.cpp',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_zn_poly.cpp'],
              libraries = ['ntl', 'stdc++', 'gmp', 'zn_poly'],
              depends = ['sage/schemes/hyperelliptic_curves/hypellfrob/hypellfrob.h',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_ntl.h',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_zn_poly.h'],
              language = 'c++',
              include_dirs = ['sage/libs/ntl/',
                              'sage/schemes/hyperelliptic_curves/hypellfrob/']),

    Extension('sage.schemes.toric.divisor_class',
              sources = ['sage/schemes/toric/divisor_class.pyx'],
              libraries = ['gmp']),

    ################################
    ##
    ## sage.sets
    ##
    ################################

    Extension('sage.sets.disjoint_set',
              sources = ['sage/sets/disjoint_set.pyx'],
              libraries = ['gmp', 'flint'],
              extra_compile_args = ['-std=c99'],
              depends = flint_depends),

    Extension('sage.sets.recursively_enumerated_set',
              sources = ['sage/sets/recursively_enumerated_set.pyx']),

    ################################
    ##
    ## sage.stats
    ##
    ################################

    Extension('sage.stats.hmm.util',
              sources = ['sage/stats/hmm/util.pyx']),

    Extension('sage.stats.hmm.distributions',
              sources = ['sage/stats/hmm/distributions.pyx']),

    Extension('sage.stats.hmm.hmm',
              sources = ['sage/stats/hmm/hmm.pyx']),

    Extension('sage.stats.hmm.chmm',
              sources = ['sage/stats/hmm/chmm.pyx'],
              extra_compile_args=["-std=c99"]),

    Extension('sage.stats.intlist',
              sources = ['sage/stats/intlist.pyx']),

    Extension('sage.stats.distributions.discrete_gaussian_integer',
              sources = ['sage/stats/distributions/discrete_gaussian_integer.pyx', 'sage/stats/distributions/dgs_gauss_mp.c', 'sage/stats/distributions/dgs_gauss_dp.c', 'sage/stats/distributions/dgs_bern.c'],
              depends = ['sage/stats/distributions/dgs_gauss.h', 'sage/stats/distributions/dgs_bern.h', 'sage/stats/distributions/dgs_misc.h'],
              libraries = ['gmp', 'mpfr'],
              extra_compile_args=["-std=c99", "-D_XOPEN_SOURCE=600"],
          ),

    ################################
    ##
    ## sage.structure
    ##
    ################################

    Extension('sage.structure.category_object',
              sources = ['sage/structure/category_object.pyx']),

    Extension('sage.structure.coerce',
              sources = ['sage/structure/coerce.pyx']),

    Extension('sage.structure.coerce_actions',
              sources = ['sage/structure/coerce_actions.pyx']),

    Extension('sage.structure.coerce_dict',
              sources = ['sage/structure/coerce_dict.pyx']),

    Extension('sage.structure.coerce_maps',
              sources = ['sage/structure/coerce_maps.pyx']),

    Extension('sage.structure.debug_options',
              sources=['sage/structure/debug_options.pyx']),

    # Compile this with -Os because it works around a bug with
    # GCC-4.7.3 + Cython 0.19 on Itanium, see Trac #14452. Moreover, it
    # actually results in faster code than -O3.
    Extension('sage.structure.element',
              sources = ['sage/structure/element.pyx'],
              extra_compile_args=["-Os"]),

    Extension('sage.structure.element_wrapper',
              sources = ['sage/structure/element_wrapper.pyx']),

    Extension('sage.structure.factory',
              sources = ['sage/structure/factory.pyx']),

    Extension('sage.structure.generators',
              sources = ['sage/structure/generators.pyx']),

    Extension('sage.structure.mutability',
              sources = ['sage/structure/mutability.pyx']),

    Extension('sage.structure.misc',
              sources = ['sage/structure/misc.pyx']),

    Extension('sage.structure.parent',
              sources = ['sage/structure/parent.pyx']),

    Extension('sage.structure.parent_base',
              sources = ['sage/structure/parent_base.pyx']),

    Extension('sage.structure.parent_gens',
              sources = ['sage/structure/parent_gens.pyx']),

    Extension('sage.structure.parent_old',
              sources = ['sage/structure/parent_old.pyx']),

    Extension('sage.structure.sage_object',
              sources = ['sage/structure/sage_object.pyx']),

    ################################
    ##
    ## sage.symbolic
    ##
    ################################

    Extension('sage.symbolic.function',
              sources = ['sage/symbolic/function.pyx']),

    Extension('sage.symbolic.ring',
              sources = ['sage/symbolic/ring.pyx']),

    Extension('*', ['sage/symbolic/*.pyx']),

    ################################
    ##
    ## sage.tests
    ##
    ################################

    Extension('sage.tests.interrupt',
              sources = ['sage/tests/interrupt.pyx', 'sage/tests/c_lib.c']),

    Extension('sage.tests.stl_vector',
              sources = ['sage/tests/stl_vector.pyx'],
              libraries = ['gmp'],
              language = 'c++'),

    Extension('sage.tests.cython',
              sources = ['sage/tests/cython.pyx']),

    ################################
    ##
    ## sage.sat
    ##
    ################################

    OptionalExtension("sage.sat.solvers.cryptominisat.cryptominisat",
              ["sage/sat/solvers/cryptominisat/cryptominisat.pyx"],
              include_dirs = [SAGE_INC, SAGE_INC+"/cmsat"],
              language = "c++",
              libraries = ['cryptominisat', 'z'],
              package = 'cryptominisat'),

    OptionalExtension("sage.sat.solvers.cryptominisat.solverconf",
              ["sage/sat/solvers/cryptominisat/solverconf.pyx", "sage/sat/solvers/cryptominisat/solverconf_helper.cpp"],
              include_dirs = [SAGE_INC, SAGE_INC+"/cmsat"],
              language = "c++",
              libraries = ['cryptominisat', 'z'],
              package = 'cryptominisat'),

    Extension('sage.sat.solvers.satsolver',
              sources = ['sage/sat/solvers/satsolver.pyx']),

    ################################
    ##
    ## sage.schemes
    ##
    ################################

    Extension('sage.schemes.projective.projective_morphism_helper',
              sources = ['sage/schemes/projective/projective_morphism_helper.pyx']),
]

# Add auto-generated modules
import sage_setup.autogen.interpreters
ext_modules += sage_setup.autogen.interpreters.modules
