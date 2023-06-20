include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
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
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cypari
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_mpmath
    SPKG_INSTALL_REQUIRES_numpy
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_modules
    SPKG_INSTALL_REQUIRES_sage_conf

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
test            = SPKG_INSTALL_REQUIRES_sagemath_repl

# extras by libraries
axiom           = # FIXME
giac            = SPKG_INSTALL_REQUIRES_sagemath_giac
ginac           = # no extra needed, same as pynac
maxima          = # no extra needed
primecount      = SPKG_INSTALL_REQUIRES_primecountpy
pynac           = # no extra needed
singular        = # no extra needed
sympy           = SPKG_INSTALL_REQUIRES_sympy

# extras by other features
plot            = SPKG_INSTALL_REQUIRES_sagemath_plot
