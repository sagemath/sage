# -*- conf-unix -*-
[metadata]
name = sagemath-flint
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Fast computations with FLINT and arb
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
        sagemath_categories \
        | sed "2,\$s/^/    /;"')dnl

packages =
    sage.libs.arb
    sage.libs.flint
    sage.matrix
    sage.rings
    sage.rings.number_field
    sage.rings.padics
    sage.rings.polynomial

[options.package_data]

sage.libs.arb =
    *.pxd

sage.libs.flint =
    *.pxd

sage.matrix =
    matrix_complex_ball_dense.pxd

sage.rings =
    *_arb.pxd

sage.rings.number_field =
    number_field_element_quadratic.pxd

sage.rings.padics =
    *_flint_*.pxd

sage.rings.polynomial =
    *_flint.pxd
    *_arb.pxd
