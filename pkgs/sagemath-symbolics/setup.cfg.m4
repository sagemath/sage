# -*- conf-unix -*-
[metadata]
name = sagemath-symbolics
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Symbolic calculus
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    esyscmd(`sage-get-system-packages install-requires \
        gmpy2          \
        cypari         \
        cysignals      \
        mpmath         \
        numpy          \
        sagemath_categories \
        sagemath_modules \
        sage_conf \
        | sed "2,\$s/^/    /;"')dnl

[options.package_data]

sage.interfaces =
    sage-maxima.lisp

sage =
    ext_data/*
    ext_data/kenzo/*
    ext_data/singular/*
    ext_data/singular/function_field/*
    ext_data/magma/*
    ext_data/magma/latex/*
    ext_data/magma/sage/*

[options.extras_require]
test = esyscmd(`sage-get-system-packages install-requires sagemath_repl')

# extras by libraries
axiom           = # FIXME
giac            = esyscmd(`sage-get-system-packages install-requires sagemath_giac')
ginac           = # no extra needed, same as pynac
maxima          = # no extra needed
primecount      = esyscmd(`sage-get-system-packages install-requires primecountpy')
pynac           = # no extra needed
singular        = # no extra needed
sympy           = esyscmd(`sage-get-system-packages install-requires sympy')

# extras by other features
plot            = esyscmd(`sage-get-system-packages install-requires sagemath_plot')
