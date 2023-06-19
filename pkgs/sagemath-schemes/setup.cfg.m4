# -*- conf-unix -*-
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
    esyscmd(`sage-get-system-packages install-requires \
        gmpy2          \
        cysignals      \
        memory_allocator \
        scipy            \
        sagemath_modules \
        | sed "2,\$s/^/    /;"')dnl

[options.extras_require]
test = esyscmd(`sage-get-system-packages install-requires sagemath_repl')
