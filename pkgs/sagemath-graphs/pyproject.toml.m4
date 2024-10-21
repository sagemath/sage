include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_meson_python
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_pkgconfig
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
]
build-backend = "mesonpy"

[project]
name = "sagemath-graphs"
description = "Sage: Open Source Mathematics Software: Graphs, posets, hypergraphs, designs, abstract complexes, combinatorial polyhedra, abelian sandpiles, quivers"
dependencies = [
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_sagemath_categories
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
test        = [SPKG_INSTALL_REQUIRES_sagemath_repl]

# libraries
bliss       = [SPKG_INSTALL_REQUIRES_sagemath_bliss]
gap         = [SPKG_INSTALL_REQUIRES_sagemath_gap]
igraph      = [SPKG_INSTALL_REQUIRES_python_igraph]
mcqd        = [SPKG_INSTALL_REQUIRES_sagemath_mcqd]
networkx    = [SPKG_INSTALL_REQUIRES_networkx]
pari        = [SPKG_INSTALL_REQUIRES_sagemath_pari]
tdlib       = [SPKG_INSTALL_REQUIRES_sagemath_tdlib]

# features
combinat    = [SPKG_INSTALL_REQUIRES_sagemath_combinat]
editor      = [SPKG_INSTALL_REQUIRES_phitigra]
homology    = [SPKG_INSTALL_REQUIRES_sagemath_modules]
mip         = [SPKG_INSTALL_REQUIRES_sagemath_polyhedra]
modules     = [SPKG_INSTALL_REQUIRES_sagemath_modules]
plot        = [SPKG_INSTALL_REQUIRES_sagemath_plot]
polyhedra   = [SPKG_INSTALL_REQUIRES_sagemath_polyhedra]
repl        = [SPKG_INSTALL_REQUIRES_sagemath_repl]
sat         = [SPKG_INSTALL_REQUIRES_sagemath_combinat]

standard    = ["sagemath-graphs[combinat,databases,mip,modules,plot,polyhedra,repl]"]

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
  "pkg:generic/boost",
  "pkg:generic/gmp",
  "pkg:generic/mpc",
  "pkg:generic/mpfr",
]

dependencies = [
]
