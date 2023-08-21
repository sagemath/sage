include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-meataxe
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Matrices over small finite fields with meataxe
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12

packages =
    sage.libs
    sage.matrix

[options.package_data]
sage.libs =
    meataxe.pxd

sage.matrix =
    matrix_gfpn_dense.pxd
