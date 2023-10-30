include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_setuptools_wheel
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_pkgconfig
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_categories
]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-groups"
description = "Sage: Open Source Mathematics Software: Groups and Invariant Theory"
dependencies = [
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_gap
    SPKG_INSTALL_REQUIRES_sagemath_modules
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
test            = [SPKG_INSTALL_REQUIRES_sagemath_repl]

# extras by packages
coxeter3        = [SPKG_INSTALL_REQUIRES_sagemath_coxeter3]
gap             = []  # no extra needed

# extras by groups_catalog
additive        = []  # no extra needed
affine          = []  # no extra needed
lie             = []  # FIXME
matrix          = []  # no extra needed
permutation     = []  # no extra needed
presentation    = []  # no extra needed

# extras by other features
representations = [SPKG_INSTALL_REQUIRES_sagemath_combinat]
semigroups      = [SPKG_INSTALL_REQUIRES_sagemath_combinat]

# the whole package
standard        = ["sagemath-groups[additive,matrix,representations,semigroups]"]

[tool.setuptools]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
