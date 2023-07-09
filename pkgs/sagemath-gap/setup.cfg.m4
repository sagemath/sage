include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-gap
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Computational Group Theory with GAP
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    SPKG_INSTALL_REQUIRES_cypari
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_categories

[options.package_data]

sage.libs.gap =
    sage.gaprc

sage.ext_data =
    *
