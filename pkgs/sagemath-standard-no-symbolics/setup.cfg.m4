include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-standard-no-symbolics
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Sage library without the symbolics subsystem
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    SPKG_INSTALL_REQUIRES_sage_conf
    SPKG_INSTALL_REQUIRES_six
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_combinat
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_flint
    SPKG_INSTALL_REQUIRES_sagemath_gap
    SPKG_INSTALL_REQUIRES_sagemath_glpk
    SPKG_INSTALL_REQUIRES_sagemath_graphs
    SPKG_INSTALL_REQUIRES_sagemath_groups
    SPKG_INSTALL_REQUIRES_sagemath_homfly
    SPKG_INSTALL_REQUIRES_sagemath_modules
    SPKG_INSTALL_REQUIRES_sagemath_mpmath
    SPKG_INSTALL_REQUIRES_sagemath_ntl
    SPKG_INSTALL_REQUIRES_sagemath_objects
    SPKG_INSTALL_REQUIRES_sagemath_pari
    SPKG_INSTALL_REQUIRES_sagemath_polyhedra
    SPKG_INSTALL_REQUIRES_sagemath_repl
    SPKG_INSTALL_REQUIRES_sagemath_schemes
    SPKG_INSTALL_REQUIRES_sagemath_singular
dnl From build/pkgs/sagelib/dependencies
    SPKG_INSTALL_REQUIRES_cypari
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_jinja2
    SPKG_INSTALL_REQUIRES_jupyter_core
    SPKG_INSTALL_REQUIRES_lrcalc_python
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_numpy
    SPKG_INSTALL_REQUIRES_pkgconfig
    SPKG_INSTALL_REQUIRES_pplpy
    SPKG_INSTALL_REQUIRES_primecountpy
    SPKG_INSTALL_REQUIRES_requests
dnl From Makefile.in: SAGERUNTIME
    SPKG_INSTALL_REQUIRES_ipython
    SPKG_INSTALL_REQUIRES_pexpect
dnl From Makefile.in: DOC_DEPENDENCIES
    SPKG_INSTALL_REQUIRES_sphinx
    SPKG_INSTALL_REQUIRES_networkx
    SPKG_INSTALL_REQUIRES_scipy
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
    bin/sage-cython
    # Only makes sense in sage-the-distribution. TODO: Move to another installation script.
    bin/sage-list-packages
    # Uncategorized scripts in alphabetical order
    bin/sage-grep
    bin/sage-grepdoc
    bin/sage-rebase.bat
    bin/sage-rebase.sh
    bin/sage-rebaseall.bat
    bin/sage-rebaseall.sh
    bin/sage-update-version

[options.package_data]

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
