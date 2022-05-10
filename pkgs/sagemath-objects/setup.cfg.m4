# -*- conf-unix -*-
[metadata]
name = sagemath-objects
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Sage objects, elements, parents, categories, coercion, metaclasses
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
    Programming Language :: Python :: 3.7
    Programming Language :: Python :: 3.8
    Programming Language :: Python :: 3.9
    Programming Language :: Python :: 3.10
    Programming Language :: Python :: Implementation :: CPython
    Topic :: Scientific/Engineering :: Mathematics

[options]
python_requires = >=3.7, <3.11
install_requires =
    esyscmd(`sage-get-system-packages install-requires \
        gmpy2          \
        cysignals      \
        | sed "2,\$s/^/    /;"')dnl

[options.extras_require]
# Currently we do not use the sage doctester to test sagemath-objects,
# so we do not list sagemath-repl here.
test =


[options.package_data]
sage.cpython =
    pyx_visit.h
    string_impl.h
    cython_metaclass.h
    python_debug.h

sage.rings =
    integer_fake.h
