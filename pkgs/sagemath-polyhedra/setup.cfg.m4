# -*- conf-unix -*-
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
    esyscmd(`sage-get-system-packages install-requires \
        gmpy2          \
        cysignals      \
        pplpy          \
        memory_allocator \
        sagemath_modules \
        sagemath_glpk \
        | sed "2,\$s/^/    /;"')dnl

[options.extras_require]
test = esyscmd(`sage-get-system-packages install-requires sagemath_repl')

# polyhedral libraries
4ti2        = # FIXME
cddlib      = # FIXME
latte_int   = # FIXME
normaliz    = esyscmd(`sage-get-system-packages install-requires pynormaliz')
polymake    = esyscmd(`sage-get-system-packages install-requires jupymake')
ppl         =                   # no extra required
topcom      = # FIXME

# optimization libraries
cbc         = sagemath-polyhedra[cbc_sage]
cbc_sage    = esyscmd(`sage-get-system-packages install-requires sage_numerical_backends_coin')
coin        = sagemath-polyhedra[cbc_sage]
coin_sage   = sagemath-polyhedra[cbc_sage]
cplex       = sagemath-polyhedra[cplex_sage]
cplex_sage  = esyscmd(`sage-get-system-packages install-requires sage_numerical_backends_cplex')
cvxopt      = sagemath-polyhedra[cvxopt_sage]
cvxopt_sage = esyscmd(`sage-get-system-packages install-requires cvxopt')
cvxpy       = esyscmd(`sage-get-system-packages install-requires cvxpy')
glpk        = sagemath-polyhedra[glpk_sage]
glpk_sage   =           # no extra required
gurobi      = sagemath-polyhedra[gurobi_sage]
gurobi_sage = esyscmd(`sage-get-system-packages install-requires sage_numerical_backends_gurobi')
scip        = esyscmd(`sage-get-system-packages install-requires pyscipopt')

# supported rings
QQ =
ZZ =
RDF = # FIXME: cddlib
NumberField = # FIXME

# features
