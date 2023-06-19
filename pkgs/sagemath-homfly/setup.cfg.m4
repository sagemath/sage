# -*- conf-unix -*-
[metadata]
name = sagemath-homfly
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Homfly polynomials of knots/links with libhomfly
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    esyscmd(`sage-get-system-packages install-requires \
        cysignals \
        | sed "2,\$s/^/    /;"')dnl
