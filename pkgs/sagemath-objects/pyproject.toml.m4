include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_setuptools
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_pkgconfig
]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-objects"
description = "Sage: Open Source Mathematics Software: Sage objects, elements, parents, categories, coercion, metaclasses"
dependencies = [
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.optional-dependencies]
# Currently we do not use the sage doctester to test sagemath-objects,
# so we do not list sagemath-repl here.
test = []

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[tool.setuptools]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}

[tool.setuptools.package-data]
"sage.cpython" = [
    "pycore_long.h",
    "pyx_visit.h",
    "string_impl.h",
    "cython_metaclass.h",
    "python_debug.h",
]
"sage.ext" = [
    "mod_int.h",
]
"sage.rings" = ["integer_fake.h"]

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
