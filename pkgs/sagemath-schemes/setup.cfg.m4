include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-schemes
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Schemes, varieties, elliptic curves, algebraic Riemann surfaces, modular forms, arithmetic dynamics
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_scipy
    SPKG_INSTALL_REQUIRES_sagemath_modules

[options.extras_require]
test = SPKG_INSTALL_REQUIRES_sagemath_repl
