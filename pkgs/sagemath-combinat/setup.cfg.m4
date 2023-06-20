include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-combinat
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Algebraic combinatorics, combinatorial representation theory
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_sagemath_categories

[options.extras_require]
test            = SPKG_INSTALL_REQUIRES_sagemath_repl

# by library
lrcalc          = SPKG_INSTALL_REQUIRES_lrcalc_python
symmetrica      =

# by feature
graphs          = SPKG_INSTALL_REQUIRES_sagemath_graphs
modules         = SPKG_INSTALL_REQUIRES_sagemath_modules

# everything
standard        = sagemath-combinat[lrcalc,symmetrica,graphs,modules]
