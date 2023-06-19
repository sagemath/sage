# -*- conf-unix -*-
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
    esyscmd(`sage-get-system-packages install-requires \
        gmpy2          \
        cysignals      \
        memory_allocator \
        sagemath_categories \
        | sed "2,\$s/^/    /;"')dnl

[options.extras_require]
test = esyscmd(`sage-get-system-packages install-requires sagemath_repl')

lrcalc = esyscmd(`sage-get-system-packages install-requires lrcalc_python')
symmetrica =

graphs = esyscmd(`sage-get-system-packages install-requires sagemath_graphs')
modules = esyscmd(`sage-get-system-packages install-requires sagemath_modules')

standard = sagemath-combinat[lrcalc,symmetrica,graphs,modules]
