#########################################################
### Library order
#########################################################

import os

# BEGIN copied from module_list.py (but #29706 removes the original).
# TODO: When #29706 is merged, simplify this module using the expanded cython_aliases.

import pkgconfig

# CBLAS can be one of multiple implementations
cblas_pc = pkgconfig.parse('cblas')
cblas_libs = cblas_pc['libraries']
cblas_library_dirs = cblas_pc['library_dirs']
cblas_include_dirs = cblas_pc['include_dirs']

# TODO: Remove Cygwin hack by installing a suitable cblas.pc
if os.path.exists('/usr/lib/libblas.dll.a'):
    cblas_libs = ['gslcblas']

# LAPACK can be one of multiple implementations
lapack_pc = pkgconfig.parse('lapack')
lapack_libs = lapack_pc['libraries']
lapack_library_dirs = lapack_pc['library_dirs']
lapack_include_dirs = lapack_pc['include_dirs']

# GD image library
gd_pc = pkgconfig.parse('gdlib')
gd_libs = gd_pc['libraries']
gd_library_dirs = gd_pc['library_dirs']
gd_include_dirs = gd_pc['include_dirs']

# PNG image library
png_pc = pkgconfig.parse('libpng')
png_libs = png_pc['libraries']
png_library_dirs = png_pc['library_dirs']
png_include_dirs = png_pc['include_dirs']

# zlib
try:
    zlib_pc = pkgconfig.parse('zlib')
except pkgconfig.PackageNotFoundError:
    from collections import defaultdict
    zlib_pc = defaultdict(list, {'libraries': ['z']})
zlib_libs = zlib_pc['libraries']
zlib_library_dirs = zlib_pc['library_dirs']
zlib_include_dirs = zlib_pc['include_dirs']

#########################################################
### M4RI flags
#########################################################

m4ri_pc = pkgconfig.parse('m4ri')
m4ri_libs = m4ri_pc['libraries']
m4ri_library_dirs = m4ri_pc['library_dirs']
m4ri_include_dirs = m4ri_pc['include_dirs']

m4ri_extra_compile_args = pkgconfig.cflags('m4ri').split()
try:
    m4ri_extra_compile_args.remove("-pedantic")
except ValueError:
    pass

# END copied from module_list.py (but #29706 removes the original).


# This list defines the *order* of linking libraries. A library should
# be put *before* any library it links to. Cython allows
# defining libraries using "# distutils: libraries = LIB". However, if
# there are multiple libraries, the order is undefined so we need to
# manually reorder the libraries according to this list. The order is
# important in particular for Cygwin. Any libraries which are not
# listed here will be added at the end of the list (without changing
# their relative order).
from sage.env import cython_aliases
aliases = cython_aliases()

arb_dylib_name = aliases["ARB_LIBRARY"]
library_order_list = aliases["SINGULAR_LIBRARIES"] + [
    "giac", "intl", "curl",
    "ec", "ecm"
] + aliases["LINBOX_LIBRARIES"] + aliases["FFLASFFPACK_LIBRARIES"] + aliases["GSL_LIBRARIES"] + [
    "pari", "flint", "ratpoints", "ecl", "glpk", "ppl",
    arb_dylib_name, "mpfi", "mpfr", "mpc", "ntl", "gmp", "gmpxx",
    "brial",
    "brial_groebner",
    "m4rie",
] + m4ri_libs + [
    "zn_poly", "gap",
] + gd_libs + png_libs + [
    "m", "readline", "Lfunction" ,
] + cblas_libs + zlib_libs

# Make a dict with library:order pairs, where the order are negative
# integers sorted according to library_order_list. When sorting,
# unlisted libraries have order 0, so they appear after the libraries
# in library_order_list.
n = len(library_order_list)
library_order = {}
for i in range(n):
    lib = library_order_list[i]
    library_order[lib] = i-n
