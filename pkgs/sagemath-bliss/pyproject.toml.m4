include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_setuptools
    SPKG_INSTALL_REQUIRES_sagelib
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_pkgconfig
]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-bliss"
description = "Sage: Open Source Mathematics Software: Graph (iso/auto)morphisms with bliss"
dependencies = [
    SPKG_INSTALL_REQUIRES_cysignals
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
test = [SPKG_INSTALL_REQUIRES_sagemath_repl]

[tool.setuptools]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
