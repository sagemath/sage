include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_setuptools_wheel
    SPKG_INSTALL_REQUIRES_pkgconfig
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_objects
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_cysignals
]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-glpk"
description = "Sage: Open Source Mathematics Software: Linear and mixed integer linear optimization backend using GLPK"
dependencies = [
    SPKG_INSTALL_REQUIRES_sagemath_objects
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_cysignals
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
test = [
    # No test requirements; see comment in tox.ini
    # SPKG_INSTALL_REQUIRES_sagemath_repl
    # SPKG_INSTALL_REQUIRES_sagemath_polyhedra
]

[tool.setuptools]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
