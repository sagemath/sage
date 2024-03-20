include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
# sage_conf is needed for library name of the flint library.
requires = [
    SPKG_INSTALL_REQUIRES_meson_python
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_sage_conf
    SPKG_INSTALL_REQUIRES_numpy
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_modules
    SPKG_INSTALL_REQUIRES_sagemath_ntl
    SPKG_INSTALL_REQUIRES_sagemath_pari
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cypari
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_pkgconfig
]
build-backend = "mesonpy"

[project]
name = "sagemath-flint"
description = "Sage: Open Source Mathematics Software: Fast computations with FLINT and arb"
dependencies = [
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_ntl
    SPKG_INSTALL_REQUIRES_cypari
    SPKG_INSTALL_REQUIRES_numpy
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
test = [
     SPKG_INSTALL_REQUIRES_sagemath_repl
     SPKG_INSTALL_REQUIRES_sagemath_modules
]

[tool.setuptools]
packages = [
    "sage.libs.arb",
    "sage.libs.flint",
    "sage.libs.mpfi",
    "sage.matrix",
    "sage.rings",
    "sage.rings.number_field",
    "sage.rings.padics",
    "sage.rings.polynomial",
]
include-package-data = false

[tool.setuptools.package-data]
"sage.libs.arb" = ["*.pxd"]
"sage.libs.flint" = ["*.pxd"]
"sage.matrix" = ["matrix_complex_ball_dense.pxd"]
"sage.rings" = ["*_arb.pxd"]
"sage.rings.number_field" = ["number_field_element_quadratic.pxd"]
"sage.rings.padics" = ["*_flint_*.pxd"]
"sage.rings.polynomial" = [
    "*_flint.pxd",
    "*_arb.pxd",
]

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
  "pkg:generic/flint",
  "pkg:generic/gmp",
  "pkg:generic/mpc",
  "pkg:generic/mpfr",
]

dependencies = [
]
