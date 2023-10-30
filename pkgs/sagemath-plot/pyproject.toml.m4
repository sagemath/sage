include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_setuptools_wheel
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_pkgconfig
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_modules
    SPKG_INSTALL_REQUIRES_numpy
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-plot"
description = "Sage: Open Source Mathematics Software: Plotting and graphics with Matplotlib, Three.JS, etc."
dependencies = [
    SPKG_INSTALL_REQUIRES_sage_conf
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_modules
    SPKG_INSTALL_REQUIRES_numpy
    SPKG_INSTALL_REQUIRES_scipy
    SPKG_INSTALL_REQUIRES_pillow
    SPKG_INSTALL_REQUIRES_matplotlib
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
test            = [SPKG_INSTALL_REQUIRES_sagemath_repl]

# extras by libraries
jsmol           = [SPKG_INSTALL_REQUIRES_jupyter_jsmol]
matplotlib      = []  # no extra needed

# extras by other features
graphs          = [SPKG_INSTALL_REQUIRES_sagemath_graphs]
polyhedra       = [SPKG_INSTALL_REQUIRES_sagemath_polyhedra]
symbolics       = [SPKG_INSTALL_REQUIRES_sagemath_symbolics]

[tool.setuptools]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
