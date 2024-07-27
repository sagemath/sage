include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_meson_python
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_sage_conf
    SPKG_INSTALL_REQUIRES_pkgconfig
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_modules
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_numpy
    SPKG_INSTALL_REQUIRES_mpmath
    SPKG_INSTALL_REQUIRES_cypari
    SPKG_INSTALL_REQUIRES_cysignals
]
build-backend = "mesonpy"

[project]
name = "sagemath-symbolics"
description = "Sage: Open Source Mathematics Software: Symbolic calculus"
dependencies = [
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cypari
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_numpy
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_flint
    SPKG_INSTALL_REQUIRES_sagemath_modules
    SPKG_INSTALL_REQUIRES_sagemath_mpmath
    SPKG_INSTALL_REQUIRES_sagemath_ntl
    SPKG_INSTALL_REQUIRES_sage_conf
    SPKG_INSTALL_REQUIRES_sympy
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
test            = [SPKG_INSTALL_REQUIRES_sagemath_repl]

# extras by libraries
axiom           = []  # FIXME
giac            = [SPKG_INSTALL_REQUIRES_sagemath_giac]
ginac           = []  # no extra needed, same as pynac
maxima          = []  # no extra needed
ntl             = []  # no extra needed
primecount      = [SPKG_INSTALL_REQUIRES_primecountpy]
pynac           = []  # no extra needed
singular        = []  # no extra needed
sympy           = []  # no extra needed

# extras by other features
plot            = [SPKG_INSTALL_REQUIRES_sagemath_plot]

[tool.setuptools]
include-package-data = false

[tool.setuptools.package-data]
"sage.interfaces" = ["sage-maxima.lisp"]
sage = [
    "ext_data/*",
    "ext_data/kenzo/*",
    "ext_data/singular/*",
    "ext_data/singular/function_field/*",
    "ext_data/magma/*",
    "ext_data/magma/latex/*",
    "ext_data/magma/sage/*",
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
  "pkg:generic/ecl",
  "pkg:generic/gmp",
  "pkg:generic/maxima",
  "pkg:generic/mpc",
  "pkg:generic/mpfr",
  "pkg:generic/singular",
]

dependencies = [
]
