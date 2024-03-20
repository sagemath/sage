include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_meson_python
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_pkgconfig
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_modules
    SPKG_INSTALL_REQUIRES_scipy
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
]
build-backend = "mesonpy"

[project]
name = "sagemath-schemes"
description = "Sage: Open Source Mathematics Software: Schemes, varieties, elliptic curves, algebraic Riemann surfaces, modular forms, arithmetic dynamics"
dependencies = [
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_scipy
    SPKG_INSTALL_REQUIRES_sagemath_modules
    SPKG_INSTALL_REQUIRES_sagemath_singular
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
test    = [SPKG_INSTALL_REQUIRES_sagemath_repl]

# extras by packages (same as sagemath-modules)
flint   = [SPKG_INSTALL_REQUIRES_sagemath_flint]
linbox  = []  # FIXME
m4ri    = []  # FIXME
m4rie   = []  # FIXME
meataxe = [SPKG_INSTALL_REQUIRES_sagemath_meataxe]
mpfi    = []  # FIXME
ntl     = [SPKG_INSTALL_REQUIRES_sagemath_ntl]
numpy   = [SPKG_INSTALL_REQUIRES_numpy]
pari    = [SPKG_INSTALL_REQUIRES_sagemath_pari]

# extras by packages (specific to sagemath-schemes)

singular = []  # no extra needed

# extras by rings; same as in sagemath-modules
RDF     = ["sagemath-schemes[numpy]"]
CDF     = ["sagemath-schemes[numpy]"]
RR      = []  # no extra needed
CC      = []  # no extra needed
RIF     = []
CIF     = []
RBF     = ["sagemath-schemes[flint]"]
CBF     = ["sagemath-schemes[flint]"]
GF      = ["sagemath-schemes[pari]"]
GF2     = ["sagemath-schemes[m4ri]"]
GF2e    = ["sagemath-schemes[m4rie]"]
GF2n    = ["sagemath-schemes[m4rie]"]
GFpn    = ["sagemath-schemes[meataxe]"]
QQbar   = ["sagemath-schemes[NumberField]"]
AA      = ["sagemath-schemes[NumberField]"]
UCF     = ["sagemath-schemes[NumberField]"]
Zp      = []  # FIXME
Qp      = ["sagemath-schemes[Zp]"]
Zq      = ["sagemath-schemes[Zp]"]
Qq      = ["sagemath-schemes[Zp]"]
FiniteField     = ["sagemath-schemes[GF]"]
NumberField     = []  # FIXME
QuadraticField  = ["sagemath-schemes[NumberField]"]
CyclotomicField = ["sagemath-schemes[NumberField]"]

# extras by features
toric           = [SPKG_INSTALL_REQUIRES_sagemath_polyhedra
                   SPKG_INSTALL_REQUIRES_sagemath_graphs]
padics          = ["sagemath-schemes[Zp]"]

# the whole package
standard        = ["sagemath-schemes[toric,padics,NumberField,FiniteField,flint,linbox,mpfi,ntl,numpy,pari,singular]"]

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
