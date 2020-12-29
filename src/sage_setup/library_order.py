#########################################################
### Library order
#########################################################

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
] + aliases["M4RI_LIBRARIES"] + [
    "zn_poly", "gap",
] + aliases["GDLIB_LIBRARIES"] + aliases["LIBPNG_LIBRARIES"] + [
    "m", "readline", "Lfunction" ,
] + aliases["CBLAS_LIBRARIES"] + aliases["ZLIB_LIBRARIES"]

# Make a dict with library:order pairs, where the order are negative
# integers sorted according to library_order_list. When sorting,
# unlisted libraries have order 0, so they appear after the libraries
# in library_order_list.
n = len(library_order_list)
library_order = {}
for i in range(n):
    lib = library_order_list[i]
    library_order[lib] = i-n
