include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
#
# Note we include numpy here to build some modules that cimport numpy,
# but it is not part of the install-requires.
requires = [
    "sage_setup[autogen]",
    SPKG_INSTALL_REQUIRES_meson_python
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_pkgconfig
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_mpmath
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_numpy
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
]
build-backend = "mesonpy"

[project]
name = "sagemath-modules"
description = "Sage: Open Source Mathematics Software: Vectors, matrices, tensors, vector spaces, affine spaces, modules and algebras, additive groups, quadratic forms, homology, coding theory, matroids"
dependencies = [
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_mpmath
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
test    = [SPKG_INSTALL_REQUIRES_sagemath_repl]

# extras by packages
flint   = [SPKG_INSTALL_REQUIRES_sagemath_flint]
gap     = [SPKG_INSTALL_REQUIRES_sagemath_gap]
gsl     = []  # No extra needed
linbox  = [SPKG_INSTALL_REQUIRES_sagemath_linbox]
m4ri    = ["sagemath-modules[linbox]"]
m4rie   = ["sagemath-modules[linbox]"]
meataxe = [SPKG_INSTALL_REQUIRES_sagemath_meataxe]
mpfi    = ["sagemath-modules[flint]"]
mpfr    = []  # No extra needed
mpmath  = []  # No extra needed
ntl     = [SPKG_INSTALL_REQUIRES_sagemath_ntl]
numpy   = [SPKG_INSTALL_REQUIRES_numpy]
pari    = [SPKG_INSTALL_REQUIRES_sagemath_pari]

# extras by rings
RDF     = ["sagemath-modules[numpy]"]
CDF     = ["sagemath-modules[numpy]"]
RR      = []  # No extra needed
CC      = []  # No extra needed
RIF     = ["sagemath-modules[flint]"]
CIF     = ["sagemath-modules[flint]"]
RBF     = ["sagemath-modules[flint]"]
CBF     = ["sagemath-modules[flint]"]
GF      = ["sagemath-modules[pari]"]
GF2     = ["sagemath-modules[m4ri]"]
GF2e    = ["sagemath-modules[m4rie]"]
GF2n    = ["sagemath-modules[m4rie]"]
GFpn    = ["sagemath-modules[meataxe]"]
QQbar   = ["sagemath-modules[NumberField]"]
AA      = ["sagemath-modules[NumberField]"]
UCF     = ["sagemath-modules[gap]"]
Zp      = ["sagemath-modules[pari]"]
Qp      = ["sagemath-modules[Zp]"]
Zq      = ["sagemath-modules[Zp]"]
Qq      = ["sagemath-modules[Zp]"]
SR      = [SPKG_INSTALL_REQUIRES_sagemath_symbolics]
FiniteField     = ["sagemath-modules[GF]"]
NumberField     = ["sagemath-modules[flint]"]
QuadraticField  = ["sagemath-modules[NumberField]"]
CyclotomicField = ["sagemath-modules[NumberField]"]

# extras by features
invariant   = [SPKG_INSTALL_REQUIRES_sagemath_groups]
combinat    = [SPKG_INSTALL_REQUIRES_sagemath_combinat]
padics      = ["sagemath-modules[Zp]"]

# the whole package
standard    = ["sagemath-modules[invariant,combinat,padics,NumberField,FiniteField,m4ri,m4rie,flint,linbox,numpy,mpfi,ntl,pari]"]

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
  "pkg:generic/gsl",
  "pkg:generic/mpc",
  "pkg:generic/mpfr",
]

dependencies = [
]
