include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_setuptools
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_objects
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_pkgconfig
]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-coxeter3"
description = "Sage: Open Source Mathematics Software: Coxeter groups, Bruhat ordering, Kazhdan-Lusztig polynomials with coxeter3"
dependencies = []
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
test = [
    SPKG_INSTALL_REQUIRES_sagemath_repl
    SPKG_INSTALL_REQUIRES_sagemath_combinat
    SPKG_INSTALL_REQUIRES_sagemath_graphs
    SPKG_INSTALL_REQUIRES_sagemath_groups
    SPKG_INSTALL_REQUIRES_sagemath_modules
]

[tool.setuptools]
packages = ["sage.libs.coxeter3"]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}

[tool.setuptools.package-data]
"sage.libs.coxeter3" = [
    "coxeter.pxd",
    "decl.pxd",
]
