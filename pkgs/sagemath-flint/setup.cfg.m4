include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-flint
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Fast computations with FLINT and arb
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    SPKG_INSTALL_REQUIRES_sagemath_categories

packages =
    sage.libs.arb
    sage.libs.flint
    sage.matrix
    sage.rings
    sage.rings.number_field
    sage.rings.padics
    sage.rings.polynomial

[options.package_data]

sage.libs.arb =
    *.pxd

sage.libs.flint =
    *.pxd

sage.matrix =
    matrix_complex_ball_dense.pxd

sage.rings =
    *_arb.pxd

sage.rings.number_field =
    number_field_element_quadratic.pxd

sage.rings.padics =
    *_flint_*.pxd

sage.rings.polynomial =
    *_flint.pxd
    *_arb.pxd
