include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-coxeter3
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Coxeter groups, Bruhat ordering, Kazhdan-Lusztig polynomials with coxeter3
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =

packages =
    sage.libs.coxeter3

[options.package_data]
sage.libs.coxeter3 =
    coxeter.pxd
    decl.pxd
