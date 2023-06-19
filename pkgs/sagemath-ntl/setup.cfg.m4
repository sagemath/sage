# -*- conf-unix -*-
[metadata]
name = sagemath-ntl
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Computational Number Theory with NTL
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    esyscmd(`sage-get-system-packages install-requires \
        cypari         \
        cysignals      \
        memory_allocator \
        sagemath_environment \
        sagemath_categories \
        sagemath_pari \
        | sed "2,\$s/^/    /;"')dnl

packages =
    sage.libs.ntl
    sage.rings.number_field
