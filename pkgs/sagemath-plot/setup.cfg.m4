include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-plot
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Plotting and graphics with Matplotlib, Three.JS, etc.
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
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

[options.extras_require]
test            = SPKG_INSTALL_REQUIRES_sagemath_repl

# extras by libraries
jsmol           = SPKG_INSTALL_REQUIRES_jupyter_jsmol
matplotlib      = # no extra needed

# extras by other features
graphs          = SPKG_INSTALL_REQUIRES_sagemath_graphs
polyhedra       = SPKG_INSTALL_REQUIRES_sagemath_polyhedra
symbolics       = SPKG_INSTALL_REQUIRES_sagemath_symbolics
