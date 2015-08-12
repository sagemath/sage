import os
from glob import glob
from distutils.extension import Extension
from sage.env import SAGE_LOCAL

SAGE_INC = os.path.join(SAGE_LOCAL, 'include')

#########################################################
### BLAS setup
#########################################################

## Choose cblas library -- note -- make sure to update sage/misc/cython.py
## if you change this!!
if 'SAGE_BLAS' in os.environ:
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
### Commonly used definitions and aliases
#########################################################

singular_incs = [SAGE_INC + '/singular', SAGE_INC + '/factory']

aliases = dict(
        GSL_LIBRARIES=['gsl', BLAS, BLAS2],
        INTERRUPT_DEPENDS=glob("sage/ext/interrupt/*.h"),
        )

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

#########################################################
### Library order
#########################################################

# This list defines the *order* of linking libraries. Cython allows
# defining libraries using "# distutils: libraries = LIB". However, if
# there are multiple libraries, the order is undefined so we need to
# manually reorder the libraries according to this list. The order is
# important in particular for Cygwin. Any libraries which are not
# listed here will be added at the end of the list (without changing
# their relative order). There is one exception: stdc++ is always put
# at the very end of the list.
library_order_list = [
    "singular", "ec", "ecm",
    "linboxsage", "ntl", "iml", "linbox", "givaro",
    "gsl", "pari", "flint", "ratpoints", "ecl", "glpk", "ppl",
    "arb", "mpfi", "mpfr", "mpc", "gmp", "gmpxx",
    "polybori-" + polybori_major_version,
    "polybori_groebner-" + polybori_major_version,
    "m4rie", "m4ri",
    "zn_poly", "gap",
    "gd", "png12",
    "m", "readline", "Lfunction",
    BLAS, BLAS2,
    "cryptominisat", "fplll", "z"]

# Make a dict with library:order pairs, where the order are negative
# integers sorted according to library_order_list. When sorting,
# unlisted libraries have order 0, so they appear after the libraries
# in library_order_list.
n = len(library_order_list)
library_order = {}
for i in range(n):
    lib = library_order_list[i]
    library_order[lib] = i-n

library_order["stdc++"] = 1000

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
               libraries = ["flint", "gmp", "gmpxx", "m", "ntl"]),

    Extension('sage.algebras.letterplace.free_algebra_letterplace',
              sources = ['sage/algebras/letterplace/free_algebra_letterplace.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs),

    Extension('sage.algebras.letterplace.free_algebra_element_letterplace',
              sources = ['sage/algebras/letterplace/free_algebra_element_letterplace.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs),

    Extension('sage.algebras.letterplace.letterplace_ideal',
              sources = ['sage/algebras/letterplace/letterplace_ideal.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs),

    Extension('sage.algebras.quatalg.quaternion_algebra_cython',
               sources = ['sage/algebras/quatalg/quaternion_algebra_cython.pyx'],
               language='c++',
               libraries = ["flint", "gmp", "gmpxx", "m", "ntl"]),

    ################################
    ##
    ## sage.calculus
    ##
    ################################

    Extension('*', ['sage/calculus/*.pyx']),

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
              libraries = ['flint']),

    Extension('*', ['sage/coding/**/*.pyx']),

    ################################
    ##
    ## sage.combinat
    ##
    ################################

    Extension('*', ['sage/combinat/**/*.pyx']),

    ################################
    ##
    ## sage.crypto
    ##
    ################################

    Extension('*', ['sage/crypto/*.pyx']),

    ################################
    ##
    ## sage.data_structures
    ##
    ################################

    Extension('*', ['sage/data_structures/*.pyx']),

    ################################
    ##
    ## sage.ext
    ##
    ################################

    Extension('*', ['sage/ext/**/*.pyx']),

    ################################
    ##
    ## sage.finance
    ##
    ################################

    Extension('*', ['sage/finance/*.pyx']),

    ################################
    ##
    ## sage.functions
    ##
    ################################

    Extension('sage.functions.prime_pi',
        sources = ['sage/functions/prime_pi.pyx'],
        extra_compile_args = ['-std=c99']),

    ################################
    ##
    ## sage.games
    ##
    ################################

    Extension('*', ['sage/games/*.pyx']),

    ################################
    ##
    ## sage.geometry
    ##
    ################################

    Extension('sage.geometry.point_collection',
              sources = ['sage/geometry/point_collection.pyx']),

    Extension('sage.geometry.toric_lattice_element',
              sources = ['sage/geometry/toric_lattice_element.pyx']),

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
              sources = ['sage/graphs/chrompoly.pyx']),

    Extension('sage.graphs.cliquer',
              sources = ['sage/graphs/cliquer.pyx']),

    Extension('sage.graphs.centrality',
              sources = ['sage/graphs/centrality.pyx']),

    Extension('sage.graphs.independent_sets',
              sources = ['sage/graphs/independent_sets.pyx']),

    Extension('sage.graphs.graph_decompositions.fast_digraph',
              sources = ['sage/graphs/graph_decompositions/fast_digraph.pyx']),

    Extension('sage.graphs.graph_decompositions.vertex_separation',
              sources = ['sage/graphs/graph_decompositions/vertex_separation.pyx']),

    Extension('sage.graphs.graph_decompositions.graph_products',
              sources = ['sage/graphs/graph_decompositions/graph_products.pyx']),

    Extension('sage.graphs.convexity_properties',
              sources = ['sage/graphs/convexity_properties.pyx']),

    Extension('sage.graphs.comparability',
              sources = ['sage/graphs/comparability.pyx']),

    Extension('sage.graphs.generic_graph_pyx',
              sources = ['sage/graphs/generic_graph_pyx.pyx']),

    Extension('sage.graphs.graph_generators_pyx',
              sources = ['sage/graphs/graph_generators_pyx.pyx']),

    Extension('sage.graphs.distances_all_pairs',
              sources = ['sage/graphs/distances_all_pairs.pyx']),

    Extension('sage.graphs.base.graph_backends',
              sources = ['sage/graphs/base/graph_backends.pyx']),

    Extension('sage.graphs.base.static_dense_graph',
              sources = ['sage/graphs/base/static_dense_graph.pyx']),

    Extension('sage.graphs.base.static_sparse_graph',
              sources = ['sage/graphs/base/static_sparse_graph.pyx']),

    Extension('sage.graphs.base.static_sparse_backend',
              sources = ['sage/graphs/base/static_sparse_backend.pyx']),

    Extension('sage.graphs.weakly_chordal',
              sources = ['sage/graphs/weakly_chordal.pyx']),

    Extension('sage.graphs.matchpoly',
              sources = ['sage/graphs/matchpoly.pyx'],
              libraries = ['flint'],
              extra_compile_args = ['-std=c99']),

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
              sources = ['sage/graphs/planarity.pyx'],
              libraries=['planarity']),

    Extension('sage.graphs.strongly_regular_db',
              sources = ['sage/graphs/strongly_regular_db.pyx']),

    Extension('sage.graphs.graph_decompositions.rankwidth',
              sources = ['sage/graphs/graph_decompositions/rankwidth.pyx'],
              libraries=['rw']),

    Extension('sage.graphs.graph_decompositions.bandwidth',
              sources = ['sage/graphs/graph_decompositions/bandwidth.pyx']),

    Extension('sage.graphs.graph_decompositions.cutwidth',
              sources = ['sage/graphs/graph_decompositions/cutwidth.pyx']),

    Extension('sage.graphs.spanning_tree',
              sources = ['sage/graphs/spanning_tree.pyx']),

    Extension('sage.graphs.trees',
              sources = ['sage/graphs/trees.pyx']),

    Extension('sage.graphs.genus',
              sources = ['sage/graphs/genus.pyx']),

    Extension('sage.graphs.hyperbolicity',
              sources = ['sage/graphs/hyperbolicity.pyx']),

    Extension('sage.graphs.base.c_graph',
              sources = ['sage/graphs/base/c_graph.pyx']),

    Extension('sage.graphs.base.sparse_graph',
              sources = ['sage/graphs/base/sparse_graph.pyx']),

    Extension('sage.graphs.base.dense_graph',
              sources = ['sage/graphs/base/dense_graph.pyx']),

    Extension('sage.graphs.base.boost_graph',
              sources = ['sage/graphs/base/boost_graph.pyx'],
              language = 'c++'),

    ################################
    ##
    ## sage.groups
    ##
    ################################

    Extension('*', ['sage/groups/*.pyx']),

    Extension('sage.groups.semimonomial_transformations.semimonomial_transformation',
              sources = ['sage/groups/semimonomial_transformations/semimonomial_transformation.pyx']),

    ###################################
    ##
    ## sage.groups.perm_gps
    ##
    ###################################

    Extension('sage.groups.perm_gps.permgroup_element',
              sources = ['sage/groups/perm_gps/permgroup_element.pyx']),

    Extension('sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label',
              sources = ['sage/groups/perm_gps/partn_ref/automorphism_group_canonical_label.pyx'],
              libraries = ['flint'],
              extra_compile_args = ['-std=c99']),

    Extension('sage.groups.perm_gps.partn_ref.canonical_augmentation',
              sources = ['sage/groups/perm_gps/partn_ref/canonical_augmentation.pyx'],
              libraries = ['flint'],
              extra_compile_args = ['-std=c99']),

    Extension('sage.groups.perm_gps.partn_ref.double_coset',
              sources = ['sage/groups/perm_gps/partn_ref/double_coset.pyx'],
              libraries = ['flint'],
              extra_compile_args = ['-std=c99']),

    Extension('sage.groups.perm_gps.partn_ref.refinement_binary',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_binary.pyx'],
              libraries = ['flint'],
              extra_compile_args = ['-std=c99']),

    Extension('sage.groups.perm_gps.partn_ref.refinement_graphs',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_graphs.pyx'],
              libraries = ['flint'],
              extra_compile_args = ['-std=c99']),

    Extension('sage.groups.perm_gps.partn_ref.refinement_lists',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_lists.pyx'],
              libraries = ['flint'],
              extra_compile_args = ['-std=c99']),

    Extension('sage.groups.perm_gps.partn_ref.refinement_matrices',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_matrices.pyx'],
              libraries = ['flint'],
              extra_compile_args = ['-std=c99']),

    Extension('sage.groups.perm_gps.partn_ref.refinement_python',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_python.pyx'],
              libraries = ['flint'],
              extra_compile_args = ['-std=c99']),

    Extension('sage.groups.perm_gps.partn_ref.refinement_sets',
              sources = ['sage/groups/perm_gps/partn_ref/refinement_sets.pyx'],
              libraries = ['flint'],
              extra_compile_args = ['-std=c99']),

    Extension('sage.groups.perm_gps.partn_ref2.refinement_generic',
              sources = ['sage/groups/perm_gps/partn_ref2/refinement_generic.pyx'],
              libraries = ["flint", "gmp", "gmpxx", "stdc++"],
              extra_compile_args=["-std=c99"],
              depends = ['sage/groups/perm_gps/partn_ref2/refinement_generic.h']),

    ################################
    ##
    ## sage.gsl
    ##
    ################################

    Extension('*', ['sage/gsl/*.pyx']),

    ################################
    ##
    ## sage.interacts
    ##
    ################################

    Extension('*', ['sage/interacts/*.pyx']),

    ################################
    ##
    ## sage.interfaces
    ##
    ################################

    Extension('*', ['sage/interfaces/*.pyx']),

    ################################
    ##
    ## sage.lfunctions
    ##
    ################################

    Extension('sage.lfunctions.zero_sums',
              sources = ['sage/lfunctions/zero_sums.pyx'],
              libraries = ["m","flint"]),

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
              libraries = ["ecl"],
              depends = [SAGE_INC + '/ecl/ecl.h']),

    OptionalExtension("sage.libs.fes",
             ["sage/libs/fes.pyx"],
             language = "c",
             libraries = ['fes'],
             package = 'fes'),

    Extension('sage.libs.flint.flint',
              sources = ["sage/libs/flint/flint.pyx"],
              libraries = ["flint", "gmp", "gmpxx", "m", "stdc++"],
              extra_compile_args = ["-std=c99", "-D_XPG6"]),

    Extension('sage.libs.flint.fmpz_poly',
              sources = ["sage/libs/flint/fmpz_poly.pyx"],
              libraries = ["flint", "gmp", "gmpxx", "m", "stdc++"],
              extra_compile_args = ["-std=c99", "-D_XPG6"]),

    Extension('sage.libs.flint.arith',
              sources = ["sage/libs/flint/arith.pyx"],
              libraries = ["flint", "gmp", "gmpxx", "m", "stdc++"],
              extra_compile_args = ["-std=c99", "-D_XPG6"]),

    Extension('sage.libs.fplll.fplll',
              sources = ['sage/libs/fplll/fplll.pyx'],
              libraries = ['gmp', 'mpfr', 'fplll'],
              language="c++",
              include_dirs = [SAGE_INC + '/fplll'],
              extra_compile_args=["-DFPLLL_V3_COMPAT"],
              depends = [SAGE_INC + "/fplll/fplll.h"]),

    Extension('sage.libs.gmp.pylong',
              sources = ['sage/libs/gmp/pylong.pyx']),

    Extension('sage.libs.gmp.rational_reconstruction',
              sources = ['sage/libs/gmp/rational_reconstruction.pyx']),

    Extension('sage.libs.linbox.linbox',
              sources = ['sage/libs/linbox/linbox.pyx'],
              libraries = ['linboxsage', 'ntl', 'iml', 'linbox',
                           'givaro', 'mpfr', 'gmp', 'gmpxx', BLAS, BLAS2],
              language = 'c++',
              extra_compile_args = givaro_extra_compile_args,
              depends = [os.path.join(SAGE_INC, 'givaro', 'givconfig.h')]),

    Extension('sage.libs.lcalc.lcalc_Lfunction',
              sources = ['sage/libs/lcalc/lcalc_Lfunction.pyx'],
              libraries = ['m', 'ntl', 'mpfr', 'gmp', 'gmpxx',
                           'Lfunction'],
              extra_compile_args=["-O3", "-ffast-math"],
              language = 'c++'),

    Extension('sage.libs.libecm',
              sources = ['sage/libs/libecm.pyx'],
              libraries = ['ecm'],
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
              libraries = ["ec", "pari",
                           "ntl", "gmp", "gmpxx", "stdc++", "m"]),

    Extension('sage.libs.pari.closure',
              sources = ["sage/libs/pari/closure.pyx"],
              libraries = ['pari', 'gmp']),

    Extension('sage.libs.pari.gen',
              sources = ["sage/libs/pari/gen.pyx"]),

    Extension('sage.libs.pari.handle_error',
              sources = ["sage/libs/pari/handle_error.pyx"]),

    Extension('sage.libs.pari.pari_instance',
              sources = ["sage/libs/pari/pari_instance.pyx"],
              extra_compile_args = ["-std=c99", "-D_XPG6"],
              libraries = ['flint']),

    Extension('sage.libs.ppl',
              sources = ['sage/libs/ppl.pyx', 'sage/libs/ppl_shim.cc'],
              libraries = ['ppl', 'gmpxx', 'gmp', 'm'],
              language="c++",
              depends = [SAGE_INC + "/ppl.hh"]),

    Extension('sage.libs.ratpoints',
              sources = ["sage/libs/ratpoints.pyx"],
              depends = [SAGE_INC + '/ratpoints.h'],
              libraries = ["ratpoints"]),

    Extension('sage.libs.readline',
              sources = ['sage/libs/readline.pyx'],
              libraries = ['readline']),

    Extension('sage.libs.singular.singular',
              sources = ['sage/libs/singular/singular.pyx'],
              libraries = ['givaro'] + singular_libs,
              language="c++",
              include_dirs = singular_incs,
              extra_compile_args = givaro_extra_compile_args),

    Extension('sage.libs.singular.polynomial',
              sources = ['sage/libs/singular/polynomial.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs),

    Extension('sage.libs.singular.ring',
              sources = ['sage/libs/singular/ring.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs),

    Extension('sage.libs.singular.groebner_strategy',
              sources = ['sage/libs/singular/groebner_strategy.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs),

    Extension('sage.libs.singular.function',
              sources = ['sage/libs/singular/function.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs,
              extra_compile_args = givaro_extra_compile_args),

    Extension('sage.libs.singular.option',
              sources = ['sage/libs/singular/option.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs),

    Extension('sage.libs.symmetrica.symmetrica',
              sources = ["sage/libs/symmetrica/symmetrica.pyx"],
              libraries = ["symmetrica"],
              depends = [SAGE_INC + "/symmetrica/def.h"]),

    Extension('sage.libs.mpmath.utils',
              sources = ["sage/libs/mpmath/utils.pyx"],
              libraries = ['mpfr']),

    Extension('sage.libs.mpmath.ext_impl',
              sources = ["sage/libs/mpmath/ext_impl.pyx"],
              libraries = ['mpfr']),

    Extension('sage.libs.mpmath.ext_main',
              sources = ["sage/libs/mpmath/ext_main.pyx"]),

    Extension('sage.libs.mpmath.ext_libmp',
              sources = ["sage/libs/mpmath/ext_libmp.pyx"]),

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

    Extension('sage.libs.ntl.convert',
              sources = ["sage/libs/ntl/convert.pyx"],
              libraries = ["ntl", "gmp", "gmpxx"],
              language='c++'),

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
              sources = ['sage/matrix/change_ring.pyx']),

    Extension('sage.matrix.matrix',
              sources = ['sage/matrix/matrix.pyx']),

    Extension('sage.matrix.matrix0',
              sources = ['sage/matrix/matrix0.pyx']),

    Extension('sage.matrix.matrix1',
              sources = ['sage/matrix/matrix1.pyx']),

    Extension('sage.matrix.matrix2',
              sources = ['sage/matrix/matrix2.pyx']),

    Extension('sage.matrix.matrix_complex_double_dense',
              sources = ['sage/matrix/matrix_complex_double_dense.pyx']),

    Extension('sage.matrix.matrix_cyclo_dense',
              sources = ['sage/matrix/matrix_cyclo_dense.pyx'],
              language = "c++",
              libraries=['ntl']),

    Extension('sage.matrix.matrix_dense',
              sources = ['sage/matrix/matrix_dense.pyx']),

    Extension('sage.matrix.matrix_double_dense',
              sources = ['sage/matrix/matrix_double_dense.pyx']),

    Extension('sage.matrix.matrix_generic_dense',
              sources = ['sage/matrix/matrix_generic_dense.pyx']),

    Extension('sage.matrix.matrix_generic_sparse',
              sources = ['sage/matrix/matrix_generic_sparse.pyx']),

    Extension('sage.matrix.matrix_integer_dense',
              sources = ['sage/matrix/matrix_integer_dense.pyx'],
              extra_compile_args = ['-std=c99'] + m4ri_extra_compile_args,
              libraries = ['iml', 'ntl', 'gmp', 'm', 'flint', BLAS, BLAS2],
              depends = [SAGE_INC + '/m4ri/m4ri.h']),

    Extension('sage.matrix.matrix_integer_sparse',
              sources = ['sage/matrix/matrix_integer_sparse.pyx']),

    Extension('sage.matrix.matrix_mod2_dense',
              sources = ['sage/matrix/matrix_mod2_dense.pyx'],
              libraries = ['m4ri', 'gd', 'png12', 'z'],
              extra_compile_args = m4ri_extra_compile_args,
              depends = [SAGE_INC + "/png.h", SAGE_INC + "/m4ri/m4ri.h"]),

    Extension('sage.matrix.matrix_mod2e_dense',
              sources = ['sage/matrix/matrix_mod2e_dense.pyx'],
              libraries = ['m4rie', 'm4ri', 'm'],
              depends = [SAGE_INC + "/m4rie/m4rie.h"],
              extra_compile_args = m4ri_extra_compile_args),

    Extension('sage.matrix.matrix_modn_dense_float',
              sources = ['sage/matrix/matrix_modn_dense_float.pyx'],
              language="c++",
              libraries = ['ntl', 'linbox', 'givaro', 'mpfr', 'gmpxx', 'gmp', BLAS, BLAS2],
              extra_compile_args = ['-DDISABLE_COMMENTATOR'] + givaro_extra_compile_args),

    Extension('sage.matrix.matrix_modn_dense_double',
              sources = ['sage/matrix/matrix_modn_dense_double.pyx'],
              language="c++",
              libraries = ['ntl', 'linbox', 'givaro', 'mpfr', 'gmpxx', 'gmp', BLAS, BLAS2],
              extra_compile_args = ["-D_XPG6", "-DDISABLE_COMMENTATOR"]
                    + m4ri_extra_compile_args + givaro_extra_compile_args),

    Extension('sage.matrix.matrix_modn_sparse',
              sources = ['sage/matrix/matrix_modn_sparse.pyx']),

    Extension('sage.matrix.matrix_mpolynomial_dense',
              sources = ['sage/matrix/matrix_mpolynomial_dense.pyx'],
              libraries = singular_libs,
              language="c++",
              include_dirs = singular_incs),

    Extension('sage.matrix.matrix_rational_dense',
              sources = ['sage/matrix/matrix_rational_dense.pyx'],
              extra_compile_args = ["-std=c99", "-D_XPG6"] + m4ri_extra_compile_args,
              libraries = ['iml', 'ntl', 'm', 'flint', BLAS, BLAS2],
              depends = [SAGE_INC + '/m4ri/m4ri.h']),

    Extension('sage.matrix.matrix_rational_sparse',
              sources = ['sage/matrix/matrix_rational_sparse.pyx']),

    Extension('sage.matrix.matrix_real_double_dense',
              sources = ['sage/matrix/matrix_real_double_dense.pyx']),

    Extension('sage.matrix.matrix_sparse',
              sources = ['sage/matrix/matrix_sparse.pyx']),

    Extension('sage.matrix.matrix_symbolic_dense',
              sources = ['sage/matrix/matrix_symbolic_dense.pyx']),

    Extension('sage.matrix.matrix_window',
              sources = ['sage/matrix/matrix_window.pyx']),

    Extension('sage.matrix.misc',
              sources = ['sage/matrix/misc.pyx'],
              libraries=['mpfr']),

    Extension('sage.matrix.strassen',
              sources = ['sage/matrix/strassen.pyx']),

    ################################
    ##
    ## sage.matroids
    ##
    ################################

    Extension('*', ['sage/matroids/*.pyx']),

    ################################
    ##
    ## sage.media
    ##
    ################################

    Extension('*', ['sage/media/*.pyx']),

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
              libraries = ["flint"],
              extra_compile_args = ['-std=c99']),

    Extension('sage.modular.modform.l_series_gross_zagier_coeffs',
              sources = ['sage/modular/modform/l_series_gross_zagier_coeffs.pyx']),

    Extension('sage.modular.modsym.apply',
              sources = ['sage/modular/modsym/apply.pyx'],
              libraries = ["flint", "gmp", "gmpxx", "m", "stdc++"],
              extra_compile_args=["-std=c99", "-D_XPG6"]),

    Extension('sage.modular.modsym.manin_symbol',
              sources = ['sage/modular/modsym/manin_symbol.pyx']),

    Extension('sage.modular.modsym.relation_matrix_pyx',
              sources = ['sage/modular/modsym/relation_matrix_pyx.pyx']),

    Extension('sage.modular.modsym.heilbronn',
              sources = ['sage/modular/modsym/heilbronn.pyx'],
              libraries = ["flint", "gmp", "gmpxx", "m", "stdc++"],
              extra_compile_args=["-std=c99", "-D_XPG6"]),

    Extension('sage.modular.modsym.p1list',
              sources = ['sage/modular/modsym/p1list.pyx']),

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
              ['sage/modules/vector_complex_double_dense.pyx']),

    Extension('sage.modules.vector_double_dense',
              ['sage/modules/vector_double_dense.pyx']),

    Extension('sage.modules.vector_integer_dense',
              sources = ['sage/modules/vector_integer_dense.pyx']),

    Extension('sage.modules.vector_modn_dense',
              extra_compile_args = ['-std=c99'],
              sources = ['sage/modules/vector_modn_dense.pyx']),

    Extension('sage.modules.vector_mod2_dense',
              sources = ['sage/modules/vector_mod2_dense.pyx'],
              libraries = ['m4ri', 'png12', 'gd'],
              extra_compile_args = m4ri_extra_compile_args,
              depends = [SAGE_INC + "/png.h", SAGE_INC + "/m4ri/m4ri.h"]),

    Extension('sage.modules.vector_rational_dense',
              sources = ['sage/modules/vector_rational_dense.pyx']),

    Extension('sage.modules.vector_real_double_dense',
              ['sage/modules/vector_real_double_dense.pyx']),

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
              libraries=["stdc++"]),

    Extension("sage.numerical.backends.generic_backend",
              ["sage/numerical/backends/generic_backend.pyx"],
              libraries=["stdc++"]),

    Extension("sage.numerical.backends.glpk_backend",
              ["sage/numerical/backends/glpk_backend.pyx"],
              language = 'c++',
              libraries=["gmp", "glpk", "z"]),

    Extension("sage.numerical.backends.ppl_backend",
              ["sage/numerical/backends/ppl_backend.pyx"],
              libraries=["stdc++"]),

    Extension("sage.numerical.backends.cvxopt_backend",
              ["sage/numerical/backends/cvxopt_backend.pyx"],
              libraries=["stdc++"]),

    Extension("sage.numerical.backends.glpk_graph_backend",
              ["sage/numerical/backends/glpk_graph_backend.pyx"],
              language = 'c++',
              libraries=["gmp", "glpk", "z"]),

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
    ## sage.parallel
    ##
    ################################

    Extension('*', ['sage/parallel/**/*.pyx']),

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
              sources = ['sage/plot/plot3d/implicit_surface.pyx']),

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

    Extension('*', ['sage/quadratic_forms/*.pyx']),

    ###############################
    ##
    ## sage.quivers
    ##
    ###############################

    Extension('*', ['sage/quivers/*.pyx']),

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
              libraries = ['ntl', 'pthread'],
              depends = ['sage/rings/bernmm/bern_modp.h',
                         'sage/rings/bernmm/bern_modp_util.h',
                         'sage/rings/bernmm/bern_rat.h'],
              language = 'c++',
              define_macros=[('USE_THREADS', '1'),
                             ('THREAD_STACK_SIZE', '4096')]),

    Extension('sage.rings.bernoulli_mod_p',
              sources = ['sage/rings/bernoulli_mod_p.pyx'],
              libraries=['ntl'],
              language = 'c++'),

    OptionalExtension("sage.rings.complex_ball_acb",
                      ["sage/rings/complex_ball_acb.pyx"],
                      libraries=['arb', 'mpfi', 'mpfr'],
                      include_dirs=[SAGE_INC + '/flint'],
                      package='arb'),

    Extension('sage.rings.complex_double',
              sources = ['sage/rings/complex_double.pyx'],
              extra_compile_args=["-std=c99", "-D_XPG6"],
              libraries = (['m'])),

    Extension('sage.rings.complex_interval',
              sources = ['sage/rings/complex_interval.pyx'],
              libraries = ['gmp', 'mpfi', 'mpfr']),

    Extension('sage.rings.complex_number',
              sources = ['sage/rings/complex_number.pyx'],
              libraries = ['gmp', 'mpfr']),

    Extension('sage.rings.integer',
              sources = ['sage/rings/integer.pyx'],
              libraries=['ntl', 'flint']),

    Extension('sage.rings.integer_ring',
              sources = ['sage/rings/integer_ring.pyx'],
              libraries=['ntl']),

    Extension('sage.rings.factorint',
              sources = ['sage/rings/factorint.pyx']),

    Extension('sage.rings.fast_arith',
              sources = ['sage/rings/fast_arith.pyx']),

    Extension('sage.rings.fraction_field_element',
              sources = ['sage/rings/fraction_field_element.pyx']),

    Extension('sage.rings.fraction_field_FpT',
              sources = ['sage/rings/fraction_field_FpT.pyx'],
              libraries = ["flint", "gmp", "gmpxx", "ntl", "zn_poly"],
              language = 'c++'),

    Extension('sage.rings.laurent_series_ring_element',
              sources = ['sage/rings/laurent_series_ring_element.pyx']),

    Extension('sage.rings.morphism',
              sources = ['sage/rings/morphism.pyx']),

    Extension('sage.rings.complex_mpc',
              sources = ['sage/rings/complex_mpc.pyx'],
              libraries = ['gmp', 'mpc', 'mpfr']),

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
              libraries=['ntl']),

    Extension('sage.rings.real_double',
              sources = ['sage/rings/real_double.pyx']),

    Extension('sage.rings.real_interval_absolute',
              sources = ['sage/rings/real_interval_absolute.pyx']),

    OptionalExtension("sage.rings.real_arb",
                      ["sage/rings/real_arb.pyx"],
                      libraries = ['arb', 'mpfi', 'mpfr'],
                      include_dirs = [SAGE_INC + '/flint'],
                      package = 'arb'),

    Extension('sage.rings.real_lazy',
              sources = ['sage/rings/real_lazy.pyx']),

    Extension('sage.rings.real_mpfi',
              sources = ['sage/rings/real_mpfi.pyx'],
              libraries = ['mpfi', 'mpfr']),

    Extension('sage.rings.real_mpfr',
              sources = ['sage/rings/real_mpfr.pyx'],
              libraries = ['mpfr']),

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
              sources = ['sage/rings/finite_rings/integer_mod.pyx']),

    Extension('sage.rings.finite_rings.element_givaro',
              sources = ["sage/rings/finite_rings/element_givaro.pyx"],
              libraries = ['givaro', 'ntl', 'gmpxx', 'gmp', 'm'],
              language='c++',
              extra_compile_args = givaro_extra_compile_args),

    Extension('sage.rings.finite_rings.element_ntl_gf2e',
              sources = ['sage/rings/finite_rings/element_ntl_gf2e.pyx'],
              libraries = ['ntl'],
              language = 'c++'),

    Extension('sage.rings.finite_rings.element_pari_ffelt',
              sources = ['sage/rings/finite_rings/element_pari_ffelt.pyx']),

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
              libraries=['ntl'],
              language = 'c++'),

    Extension('sage.rings.number_field.number_field_element_quadratic',
              sources = ['sage/rings/number_field/number_field_element_quadratic.pyx'],
              libraries=['ntl'],
              language = 'c++'),

    Extension('sage.rings.number_field.number_field_morphisms',
              sources = ['sage/rings/number_field/number_field_morphisms.pyx']),

    Extension('sage.rings.number_field.totallyreal',
              sources = ['sage/rings/number_field/totallyreal.pyx']),

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
              sources = ['sage/rings/padics/common_conversion.pyx']),

    Extension('sage.rings.padics.local_generic_element',
              sources = ['sage/rings/padics/local_generic_element.pyx']),

    Extension('sage.rings.padics.padic_capped_absolute_element',
              sources = ['sage/rings/padics/padic_capped_absolute_element.pyx']),

    Extension('sage.rings.padics.padic_capped_relative_element',
              sources = ['sage/rings/padics/padic_capped_relative_element.pyx']),

    Extension('sage.rings.padics.padic_ext_element',
              sources = ['sage/rings/padics/padic_ext_element.pyx'],
              libraries=['ntl', 'gmp', 'gmpxx', 'm'],
              language='c++'),

    Extension('sage.rings.padics.padic_fixed_mod_element',
              sources = ['sage/rings/padics/padic_fixed_mod_element.pyx']),

    Extension('sage.rings.padics.padic_generic_element',
              sources = ['sage/rings/padics/padic_generic_element.pyx']),

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
              include_dirs = singular_incs),

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
              include_dirs = singular_incs),

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
              libraries = ['gmp', 'ntl'],
              extra_compile_args = m4ri_extra_compile_args,
              language = 'c++',
              depends = [SAGE_INC + '/m4ri/m4ri.h']),

    Extension('sage.rings.polynomial.polynomial_zz_pex',
              sources = ['sage/rings/polynomial/polynomial_zz_pex.pyx'],
              libraries = ['ntl'],
              language = 'c++'),

    Extension('sage.rings.polynomial.polynomial_zmod_flint',
              sources = ['sage/rings/polynomial/polynomial_zmod_flint.pyx'],
              libraries = ["flint", "gmp", "gmpxx", "ntl", "zn_poly"],
              language = 'c++'),

    Extension('sage.rings.polynomial.polynomial_integer_dense_flint',
              sources = ['sage/rings/polynomial/polynomial_integer_dense_flint.pyx'],
              language = 'c++',
              libraries = ["flint", "ntl", "gmpxx", "gmp"]),

    Extension('sage.rings.polynomial.polynomial_integer_dense_ntl',
              sources = ['sage/rings/polynomial/polynomial_integer_dense_ntl.pyx'],
              libraries = ['ntl'],
              language = 'c++'),

    Extension('sage.rings.polynomial.polynomial_rational_flint',
              sources = ['sage/rings/polynomial/polynomial_rational_flint.pyx'],
              libraries = ["flint", "ntl", "gmpxx", "gmp"],
              language = 'c++'),

    Extension('sage.rings.polynomial.polynomial_modn_dense_ntl',
              sources = ['sage/rings/polynomial/polynomial_modn_dense_ntl.pyx'],
              libraries = ['ntl'],
              language = 'c++'),

    Extension('sage.rings.polynomial.polynomial_ring_homomorphism',
              sources = ['sage/rings/polynomial/polynomial_ring_homomorphism.pyx']),

    Extension('sage.rings.polynomial.pbori',
              sources = ['sage/rings/polynomial/pbori.pyx'],
              libraries=['polybori-' + polybori_major_version,
                         'polybori_groebner-' + polybori_major_version, 'm4ri', 'gd', 'png12'],
              depends = [SAGE_INC + "/polybori/" + hd + ".h" for hd in ["polybori", "config"] ] + \
                        [SAGE_INC + '/m4ri/m4ri.h'],
              extra_compile_args = polybori_extra_compile_args + m4ri_extra_compile_args,
              language = 'c++'),

    Extension('sage.rings.polynomial.polynomial_real_mpfr_dense',
              sources = ['sage/rings/polynomial/polynomial_real_mpfr_dense.pyx'],
              libraries = ['gmp', 'mpfr']),

    Extension('sage.rings.polynomial.real_roots',
              sources = ['sage/rings/polynomial/real_roots.pyx'],
              libraries=['mpfr']),

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
    ## sage.sat
    ##
    ################################

    OptionalExtension("sage.sat.solvers.cryptominisat.cryptominisat",
              sources = ["sage/sat/solvers/cryptominisat/cryptominisat.pyx"],
              include_dirs = [os.path.join(SAGE_INC, "cmsat")],
              language = "c++",
              libraries = ['cryptominisat', 'z'],
              package = 'cryptominisat'),

    OptionalExtension("sage.sat.solvers.cryptominisat.solverconf",
              sources = ["sage/sat/solvers/cryptominisat/solverconf.pyx", "sage/sat/solvers/cryptominisat/solverconf_helper.cpp"],
              include_dirs = [os.path.join(SAGE_INC, "cmsat")],
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

    Extension('sage.schemes.elliptic_curves.descent_two_isogeny',
              sources = ['sage/schemes/elliptic_curves/descent_two_isogeny.pyx'],
              extra_compile_args=["-std=c99"],
              depends = [SAGE_INC + '/ratpoints.h',
                         SAGE_INC + '/gmp.h'],
              libraries = ['flint', 'ratpoints']),

    Extension('sage.schemes.elliptic_curves.period_lattice_region',
              sources = ['sage/schemes/elliptic_curves/period_lattice_region.pyx']),

    Extension('sage.schemes.hyperelliptic_curves.hypellfrob',
              sources = ['sage/schemes/hyperelliptic_curves/hypellfrob.pyx',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/hypellfrob.cpp',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_ntl.cpp',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_zn_poly.cpp'],
              libraries = ['gmp', 'ntl', 'zn_poly'],
              depends = ['sage/schemes/hyperelliptic_curves/hypellfrob/hypellfrob.h',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_ntl.h',
                         'sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_zn_poly.h'],
              language = 'c++',
              include_dirs = ['sage/libs/ntl/',
                              'sage/schemes/hyperelliptic_curves/hypellfrob/']),

    Extension('sage.schemes.projective.projective_morphism_helper',
              sources = ['sage/schemes/projective/projective_morphism_helper.pyx']),

    Extension('sage.schemes.toric.divisor_class',
              sources = ['sage/schemes/toric/divisor_class.pyx']),

    ################################
    ##
    ## sage.sets
    ##
    ################################

    Extension('sage.sets.disjoint_set',
              sources = ['sage/sets/disjoint_set.pyx'],
              libraries = ['flint'],
              extra_compile_args = ['-std=c99']),

    Extension('sage.sets.finite_set_map_cy',
              sources=['sage/sets/finite_set_map_cy.pyx']),

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
              libraries = ['mpfr'],
              extra_compile_args=["-std=c99", "-D_XOPEN_SOURCE=600"],
          ),

    ################################
    ##
    ## sage.structure
    ##
    ################################

    # Compile this with -Os because it works around a bug with
    # GCC-4.7.3 + Cython 0.19 on Itanium, see Trac #14452. Moreover, it
    # actually results in faster code than -O3.
    Extension('sage.structure.element',
              sources = ['sage/structure/element.pyx'],
              extra_compile_args=["-Os"]),

    Extension('*', ['sage/structure/*.pyx']),

    ################################
    ##
    ## sage.symbolic
    ##
    ################################

    Extension('*', ['sage/symbolic/*.pyx']),

    ################################
    ##
    ## sage.tests
    ##
    ################################

    Extension('sage.tests.stl_vector',
              sources = ['sage/tests/stl_vector.pyx'],
              language = 'c++'),

    Extension('sage.tests.cython',
              sources = ['sage/tests/cython.pyx']),
]

# Add auto-generated modules
import sage_setup.autogen.interpreters
ext_modules += sage_setup.autogen.interpreters.modules
