include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    # Some version of sage-conf is required.
    # Note that PEP517/518 have no notion of optional sage_spkg dependencies:
    # https://github.com/pypa/pip/issues/6144
    SPKG_INSTALL_REQUIRES_sage_conf
    SPKG_INSTALL_REQUIRES_meson_python
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_cypari
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_jinja2
    SPKG_INSTALL_REQUIRES_jupyter_core
    SPKG_INSTALL_REQUIRES_numpy
    SPKG_INSTALL_REQUIRES_pkgconfig
    SPKG_INSTALL_REQUIRES_pplpy
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_objects
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_modules
    SPKG_INSTALL_REQUIRES_sagemath_polyhedra
    SPKG_INSTALL_REQUIRES_sagemath_graphs
    SPKG_INSTALL_REQUIRES_sagemath_combinat
    SPKG_INSTALL_REQUIRES_sagemath_ntl
    SPKG_INSTALL_REQUIRES_sagemath_pari
    SPKG_INSTALL_REQUIRES_sagemath_repl
    SPKG_INSTALL_REQUIRES_sagemath_schemes
    SPKG_INSTALL_REQUIRES_sagemath_singular
]
build-backend = "mesonpy"

[project]
name = "sagemath-standard-no-symbolics"
description = "Sage: Open Source Mathematics Software: Sage library without the symbolics subsystem"
dependencies = [
    SPKG_INSTALL_REQUIRES_sage_conf
    SPKG_INSTALL_REQUIRES_six
    SPKG_INSTALL_REQUIRES_sagemath_brial
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_combinat
    SPKG_INSTALL_REQUIRES_sagemath_eclib
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_flint
    SPKG_INSTALL_REQUIRES_sagemath_gap
    SPKG_INSTALL_REQUIRES_sagemath_glpk
    SPKG_INSTALL_REQUIRES_sagemath_graphs
    SPKG_INSTALL_REQUIRES_sagemath_groups
    SPKG_INSTALL_REQUIRES_sagemath_homfly
    SPKG_INSTALL_REQUIRES_sagemath_lcalc
    SPKG_INSTALL_REQUIRES_sagemath_libbraiding
    SPKG_INSTALL_REQUIRES_sagemath_libecm
    SPKG_INSTALL_REQUIRES_sagemath_linbox
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
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[tool.setuptools]
script-files = [
    "bin/sage-cython",
    # Only makes sense in sage-the-distribution. TODO: Move to another installation script.
    "bin/sage-list-packages",
    # Uncategorized scripts in alphabetical order
    "bin/sage-grep",
    "bin/sage-grepdoc",
    "bin/sage-update-version",
]
include-package-data = false

[tool.setuptools.package-data]
sage = [
    "ext_data/*",
    "ext_data/kenzo/*",
    "ext_data/images/*",
    "ext_data/doctest/*",
    "ext_data/doctest/invalid/*",
    "ext_data/gap/*",
    "ext_data/gap/joyner/*",
    "ext_data/mwrank/*",
    "ext_data/notebook-ipython/*",
    "ext_data/nbconvert/*",
    "ext_data/graphs/*",
    "ext_data/valgrind/*",
]

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
