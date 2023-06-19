# -*- conf-unix -*-
[metadata]
name = sagemath-giac
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Graph (iso/auto)morphisms with giac
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    esyscmd(`sage-get-system-packages install-requires \
        sagemath_categories \
        | sed "2,\$s/^/    /;"')dnl

packages =
    sage.libs.giac
    sage.rings

[options.package_data]

sage.libs.giac =
    *.pxd

sage.rings =
    *.pxd
