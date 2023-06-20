include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-polyhedra
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Convex polyhedra in arbitrary dimension, mixed integer linear optimization
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_pplpy
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_sagemath_modules
    SPKG_INSTALL_REQUIRES_sagemath_glpk

[options.extras_require]
test        = SPKG_INSTALL_REQUIRES_sagemath_repl

# polyhedral libraries
4ti2        = # FIXME
cddlib      = # FIXME
latte_int   = # FIXME
normaliz    = SPKG_INSTALL_REQUIRES_pynormaliz
polymake    = SPKG_INSTALL_REQUIRES_jupymake
ppl         = # no extra required
topcom      = # FIXME

# optimization libraries
cbc         = sagemath-polyhedra[cbc_sage]
cbc_sage    = SPKG_INSTALL_REQUIRES_sage_numerical_backends_coin
coin        = sagemath-polyhedra[cbc_sage]
coin_sage   = sagemath-polyhedra[cbc_sage]
cplex       = sagemath-polyhedra[cplex_sage]
cplex_sage  = SPKG_INSTALL_REQUIRES_sage_numerical_backends_cplex
cvxopt      = sagemath-polyhedra[cvxopt_sage]
cvxopt_sage = SPKG_INSTALL_REQUIRES_cvxopt
cvxpy       = SPKG_INSTALL_REQUIRES_cvxpy
glpk        = sagemath-polyhedra[glpk_sage]
glpk_sage   = # no extra required
gurobi      = sagemath-polyhedra[gurobi_sage]
gurobi_sage = SPKG_INSTALL_REQUIRES_sage_numerical_backends_gurobi
scip        = SPKG_INSTALL_REQUIRES_pyscipopt

# supported rings
QQ  =
ZZ  =
RDF = # FIXME: cddlib
NumberField = # FIXME

# features
