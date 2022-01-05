# -*- conf-unix -*-
[metadata]
name = sagemath-standard
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Standard Python Library
long_description = file: README.rst
long_description_content_type = text/x-rst
license = GNU General Public License (GPL) v2 or later
license_file = LICENSE.txt
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
    Programming Language :: Python :: Implementation :: CPython
    Topic :: Scientific/Engineering :: Mathematics

[options]
python_requires = >=3.7, <3.10
install_requires =
    sage_conf
    esyscmd(`sage-get-system-packages install-requires \
        six \
        | sed "2,\$s/^/    /;"')dnl
dnl From build/pkgs/sagelib/dependencies
    esyscmd(`sage-get-system-packages install-requires \
        cypari         \
        cysignals      \
        cython         \
        gmpy2          \
        jinja2         \
        jupyter_core   \
        numpy          \
        pkgconfig      \
        pplpy          \
        | sed "2,\$s/^/    /;"')dnl
dnl From Makefile.in: SAGERUNTIME
    esyscmd(`sage-get-system-packages install-requires \
        ipython        \
        pexpect        \
        psutil         \
        | sed "2,\$s/^/    /;"')dnl
dnl From Makefile.in: DOC_DEPENDENCIES
    esyscmd(`sage-get-system-packages install-requires \
        sphinx         \
        networkx       \
        scipy          \
        sympy          \
        matplotlib     \
        pillow         \
        mpmath         \
        ipykernel      \
        jupyter_client \
        ipywidgets     \
        | sed "2,\$s/^/    /;"')dnl
dnl Other Python packages that are standard spkg, used in doctests
    esyscmd(`sage-get-system-packages install-requires \
        rpy2           \
        fpylll         \
        | sed "2,\$s/^/    /;"')dnl
dnl pycryptosat  # Sage distribution installs it as part of cryptominisat. According to its README on https://pypi.org/project/pycryptosat/: "The pycryptosat python package compiles while compiling CryptoMiniSat. It cannot be compiled on its own, it must be compiled at the same time as CryptoMiniSat."
