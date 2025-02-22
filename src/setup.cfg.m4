include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-standard
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Standard Python Library
long_description = file: README.rst
long_description_content_type = text/x-rst
license_files = LICENSE.txt
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.11, <3.14
install_requires =
    SPKG_INSTALL_REQUIRES_six
dnl From build/pkgs/sagelib/dependencies
    SPKG_INSTALL_REQUIRES_conway_polynomials
    SPKG_INSTALL_REQUIRES_cypari
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_importlib_metadata
    SPKG_INSTALL_REQUIRES_importlib_resources
    SPKG_INSTALL_REQUIRES_jupyter_core
    SPKG_INSTALL_REQUIRES_lrcalc_python
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_numpy
    SPKG_INSTALL_REQUIRES_pkgconfig
    SPKG_INSTALL_REQUIRES_pplpy
    SPKG_INSTALL_REQUIRES_primecountpy
    SPKG_INSTALL_REQUIRES_requests
    SPKG_INSTALL_REQUIRES_typing_extensions
dnl From Makefile.in: SAGERUNTIME
    SPKG_INSTALL_REQUIRES_ipython
    SPKG_INSTALL_REQUIRES_pexpect
dnl From Makefile.in: DOC_DEPENDENCIES
    sphinx >=5.2, <9
    SPKG_INSTALL_REQUIRES_networkx
    SPKG_INSTALL_REQUIRES_scipy
    SPKG_INSTALL_REQUIRES_sympy
    SPKG_INSTALL_REQUIRES_matplotlib
    SPKG_INSTALL_REQUIRES_pillow
    SPKG_INSTALL_REQUIRES_mpmath
    SPKG_INSTALL_REQUIRES_ipykernel
    SPKG_INSTALL_REQUIRES_jupyter_client
    SPKG_INSTALL_REQUIRES_ipywidgets
    SPKG_INSTALL_REQUIRES_fpylll
dnl pycryptosat  # Sage distribution installs it as part of cryptominisat. According to its README on https://pypi.org/project/pycryptosat/: "The pycryptosat python package compiles while compiling CryptoMiniSat. It cannot be compiled on its own, it must be compiled at the same time as CryptoMiniSat."
dnl Packages with important upper version bounds
    SPKG_INSTALL_REQUIRES_ptyprocess

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
    bin/sage-grep
    bin/sage-grepdoc
    bin/sage-inline-fortran
    bin/sage-ipynb2rst
    bin/sage-ipython
    bin/sage-notebook
    bin/sage-num-threads.py
    bin/sage-preparse
    bin/sage-python
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

[options.extras_require]
R        = SPKG_INSTALL_REQUIRES_rpy2
bliss    = SPKG_INSTALL_REQUIRES_sagemath_bliss
coxeter3 = SPKG_INSTALL_REQUIRES_sagemath_coxeter3
mcqd     = SPKG_INSTALL_REQUIRES_sagemath_mcqd
meataxe  = SPKG_INSTALL_REQUIRES_sagemath_meataxe
sirocco  = SPKG_INSTALL_REQUIRES_sagemath_sirocco
tdlib    = SPKG_INSTALL_REQUIRES_sagemath_tdlib
