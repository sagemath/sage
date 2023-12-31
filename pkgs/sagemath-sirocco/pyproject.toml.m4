include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_setuptools
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_cypari
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_pkgconfig
]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-sirocco"
description = "Sage: Open Source Mathematics Software: Certified root continuation with sirocco"
dependencies = [
    SPKG_INSTALL_REQUIRES_cypari
    SPKG_INSTALL_REQUIRES_cysignals
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[tool.setuptools]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
