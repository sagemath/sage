include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_meson_python
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_pkgconfig
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_modules
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
]
build-backend = "mesonpy"

[project]
name = "sagemath-polyhedra"
description = "Sage: Open Source Mathematics Software: Convex polyhedra in arbitrary dimension, mixed integer linear optimization"
dependencies = [
    SPKG_INSTALL_REQUIRES_sage_conf
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_pplpy
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_glpk
    SPKG_INSTALL_REQUIRES_sagemath_modules
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
test        = [SPKG_INSTALL_REQUIRES_sagemath_repl]

# general libraries
flint       = [SPKG_INSTALL_REQUIRES_sagemath_flint]
fpylll      = [SPKG_INSTALL_REQUIRES_fpylll]
linbox      = [SPKG_INSTALL_REQUIRES_sagemath_linbox]
pari        = [SPKG_INSTALL_REQUIRES_sagemath_pari]

# polyhedral libraries
4ti2        = []  # FIXME
cddlib      = []  # FIXME
latte       = ["sagemath-polyhedra[latte_int]"]  # alias
latte_int   = []  # FIXME
normaliz    = [SPKG_INSTALL_REQUIRES_pynormaliz]
polymake    = [SPKG_INSTALL_REQUIRES_jupymake]
ppl         = []  # no extra required
topcom      = []  # FIXME

# optimization libraries
cbc         = ["sagemath-polyhedra[cbc_sage]"]
cbc_sage    = [SPKG_INSTALL_REQUIRES_sage_numerical_backends_coin]
coin        = ["sagemath-polyhedra[cbc_sage]"]
coin_sage   = ["sagemath-polyhedra[cbc_sage]"]
cplex       = ["sagemath-polyhedra[cplex_sage]"]
cplex_sage  = [SPKG_INSTALL_REQUIRES_sage_numerical_backends_cplex]
cvxopt      = ["sagemath-polyhedra[cvxopt_sage]"]
cvxopt_sage = [SPKG_INSTALL_REQUIRES_cvxopt]
cvxpy       = [SPKG_INSTALL_REQUIRES_cvxpy]
glpk        = ["sagemath-polyhedra[glpk_sage]"]
glpk_sage   = []  # no extra required
gurobi      = ["sagemath-polyhedra[gurobi_sage]"]
gurobi_sage = [SPKG_INSTALL_REQUIRES_sage_numerical_backends_gurobi]
scip        = [SPKG_INSTALL_REQUIRES_pyscipopt]

# supported rings
QQ          = []
ZZ          = []
RDF         = []  # FIXME: cddlib
NumberField = ["sagemath-polyhedra[flint]"]

# features
graphs      = [SPKG_INSTALL_REQUIRES_sagemath_graphs]
groups      = [SPKG_INSTALL_REQUIRES_sagemath_groups]
plot        = [SPKG_INSTALL_REQUIRES_sagemath_plot]
posets      = [SPKG_INSTALL_REQUIRES_sagemath_graphs]
toric       = [SPKG_INSTALL_REQUIRES_sagemath_graphs]

# the whole package
standard    = ["sagemath-polyhedra[flint,fpylll,linbox,glpk,graphs,groups,pari,plot,RDF]"]

[tool.setuptools]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}

[external]
# External dependencies in the format proposed by https://peps.python.org/pep-0725
build-requires = [
  "virtual:compiler/c",
  "virtual:compiler/cpp",
  "pkg:generic/pkg-config",
]

host-requires = [
  "pkg:generic/gmp",
  "pkg:generic/mpc",
  "pkg:generic/mpfr",
]

dependencies = [
]
