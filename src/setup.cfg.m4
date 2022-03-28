# -*- conf-unix -*-
[metadata]
name = sagemath-standard
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Standard Python Library
long_description = file: README.rst
long_description_content_type = text/x-rst
license = GNU General Public License (GPL) v2 or later
license_files = LICENSE.txt
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
python_requires = >=3.7, <3.11
install_requires =
    esyscmd(`sage-get-system-packages install-requires \
        sage_conf \
        six \
        | sed "2,\$s/^/    /;"')dnl'
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
        memory_allocator \
        requests       \
        | sed "2,\$s/^/    /;"')dnl'
dnl From Makefile.in: SAGERUNTIME
    esyscmd(`sage-get-system-packages install-requires \
        ipython        \
        pexpect        \
        | sed "2,\$s/^/    /;"')dnl'
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
        | sed "2,\$s/^/    /;"')dnl'
dnl Other Python packages that are standard spkg, used in doctests
    esyscmd(`sage-get-system-packages install-requires \
        rpy2           \
        fpylll         \
        | sed "2,\$s/^/    /;"')dnl'
dnl pycryptosat  # Sage distribution installs it as part of cryptominisat. According to its README on https://pypi.org/project/pycryptosat/: "The pycryptosat python package compiles while compiling CryptoMiniSat. It cannot be compiled on its own, it must be compiled at the same time as CryptoMiniSat."
dnl Packages with important upper version bounds
    esyscmd(`sage-get-system-packages install-requires \
        ptyprocess     \
        | sed "2,\$s/^/    /;"')dnl'

scripts =
    # The sage script
    bin/sage
    # Other scripts that should be in the path also for OS packaging of sage:
    bin/sage-eval
    # Included because it is useful for doctesting/coverage testing user scripts too:
    bin/sage-runtests
    bin/sage-fixdoctests
    bin/sage-coverage
    # The following is deprecated but might still be used in user package install scripts
    bin/sage-cython
    # Helper scripts invoked by sage script
    # (they would actually belong to something like libexec)
    bin/sage-cachegrind
    bin/sage-callgrind
    bin/sage-massif
    bin/sage-omega
    bin/sage-valgrind
    bin/sage-venv-config
    bin/sage-version.sh
    bin/sage-cleaner
    # Only makes sense in sage-the-distribution. TODO: Move to another installation script.
    bin/sage-list-packages
    # Uncategorized scripts in alphabetical order
    bin/math-readline
    bin/sage-env
    # sage-env-config -- installed by sage_conf
    # sage-env-config.in -- not to be installed
    bin/sage-gdb-commands
    bin/sage-grep
    bin/sage-grepdoc
    bin/sage-inline-fortran
    bin/sage-ipynb2rst
    bin/sage-ipython
    bin/sage-notebook
    bin/sage-num-threads.py
    bin/sage-preparse
    bin/sage-python
    bin/sage-rebase.bat
    bin/sage-rebase.sh
    bin/sage-rebaseall.bat
    bin/sage-rebaseall.sh
    bin/sage-run
    bin/sage-run-cython
    bin/sage-startuptime.py
    bin/sage-update-version

[options.package_data]

sage.libs.gap =
    sage.gaprc

sage.interfaces =
    sage-maxima.lisp

sage.doctest =
    tests/*

sage.repl.rich_output =
    example*

sage =
    ext_data/*
    ext_data/kenzo/*
    ext_data/singular/*
    ext_data/singular/function_field/*
    ext_data/images/*
    ext_data/doctest/*
    ext_data/doctest/invalid/*
    ext_data/gap/*
    ext_data/gap/joyner/*
    ext_data/mwrank/*
    ext_data/notebook-ipython/*
    ext_data/nbconvert/*
    ext_data/graphs/*
    ext_data/pari/*
    ext_data/pari/dokchitser/*
    ext_data/pari/buzzard/*
    ext_data/pari/simon/*
    ext_data/magma/*
    ext_data/magma/latex/*
    ext_data/magma/sage/*
    ext_data/valgrind/*
    ext_data/threejs/*
