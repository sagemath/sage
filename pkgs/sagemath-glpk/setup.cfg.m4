include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-glpk
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Linear and mixed integer linear optimization backend using GLPK
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    SPKG_INSTALL_REQUIRES_sagemath_objects
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_cysignals

[options.extras_require]

# No test requirements; see comment in tox.ini
test =
##     SPKG_INSTALL_REQUIRES_sagemath_repl
##     SPKG_INSTALL_REQUIRES_sagemath_polyhedra
