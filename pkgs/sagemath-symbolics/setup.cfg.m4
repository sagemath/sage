# -*- conf-unix -*-
[metadata]
name = sagemath-symbolics
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Symbolic calculus
long_description = file: README.rst
long_description_content_type = text/x-rst
license = GNU General Public License (GPL) v2 or later
author = The Sage Developers
author_email = sage-support@googlegroups.com
url = https://www.sagemath.org

classifiers =
    Development Status :: 6 - Mature
    Intended Audience :: Education
    Intended Audience :: Science/Research
    License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
    Programming Language :: Python :: 3 :: Only
    Programming Language :: Python :: 3.8
    Programming Language :: Python :: 3.9
    Programming Language :: Python :: 3.10
    Programming Language :: Python :: 3.11
    Programming Language :: Python :: Implementation :: CPython
    Topic :: Scientific/Engineering :: Mathematics

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
        sagemath_polyhedra \
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
