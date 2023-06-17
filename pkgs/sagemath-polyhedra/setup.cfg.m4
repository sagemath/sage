# -*- conf-unix -*-
[metadata]
name = sagemath-polyhedra
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Convex polyhedra in arbitrary dimension
long_description = file: README.rst
long_description_content_type = text/x-rst
license = GNU General Public License (GPL) v2 or later
author = The Sage Developers
author_email = sage-support@googlegroups.com
url = https://www.sagemath.org

classifiers =
    Development Status :: 6 - Mature
    Intended Audience :: Education
    Intended Audience :: Science/Research
    License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
    Programming Language :: Python :: 3 :: Only
    Programming Language :: Python :: 3.8
    Programming Language :: Python :: 3.9
    Programming Language :: Python :: 3.10
    Programming Language :: Python :: 3.11
    Programming Language :: Python :: Implementation :: CPython
    Topic :: Scientific/Engineering :: Mathematics

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
normaliz    = esyscmd(`sage-get-system-packages install-requires pynormaliz')
polymake    = esyscmd(`sage-get-system-packages install-requires jupymake')
ppl =                   # no extra required

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
