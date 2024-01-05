include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    "sage_setup[autogen]",
    # Some version of sage-conf is required.
    # Note that PEP517/518 have no notion of optional sage_spkg dependencies:
    # https://github.com/pypa/pip/issues/6144
    SPKG_INSTALL_REQUIRES_sage_conf
    SPKG_INSTALL_REQUIRES_setuptools
    SPKG_INSTALL_REQUIRES_wheel
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_cypari
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_jupyter_core
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_numpy
    SPKG_INSTALL_REQUIRES_pkgconfig
    SPKG_INSTALL_REQUIRES_pplpy
]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-standard"
description = "Sage: Open Source Mathematics Software: Standard Python Library"
dependencies = [
    SPKG_INSTALL_REQUIRES_sage_conf
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
    SPKG_INSTALL_REQUIRES_sphinx
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
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.optional-dependencies]
R = [
    SPKG_INSTALL_REQUIRES_rpy2
]

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[tool.conda-lock]
platforms = [
    'osx-64', 'linux-64', 'linux-aarch64', 'osx-arm64'
]

[tool.setuptools]
script-files = [
    # The sage script
    "bin/sage",
    # Other scripts that should be in the path also for OS packaging of sage:
    "bin/sage-eval",
    # Included because it is useful for doctesting/coverage testing user scripts too:
    "bin/sage-runtests",
    "bin/sage-fixdoctests",
    "bin/sage-coverage",
    # The following is deprecated but might still be used in user package install scripts
    "bin/sage-cython",
    # Helper scripts invoked by sage script
    # (they would actually belong to something like libexec)
    "bin/sage-cachegrind",
    "bin/sage-callgrind",
    "bin/sage-massif",
    "bin/sage-omega",
    "bin/sage-valgrind",
    "bin/sage-venv-config",
    "bin/sage-version.sh",
    "bin/sage-cleaner",
    # Only makes sense in sage-the-distribution. TODO: Move to another installation script.
    "bin/sage-list-packages",
    # Uncategorized scripts in alphabetical order
    "bin/math-readline",
    "bin/sage-env",
    # sage-env-config -- installed by sage_conf
    # sage-env-config.in -- not to be installed
    "bin/sage-grep",
    "bin/sage-grepdoc",
    "bin/sage-inline-fortran",
    "bin/sage-ipynb2rst",
    "bin/sage-ipython",
    "bin/sage-notebook",
    "bin/sage-num-threads.py",
    "bin/sage-preparse",
    "bin/sage-python",
    "bin/sage-run",
    "bin/sage-run-cython",
    "bin/sage-startuptime.py",
    "bin/sage-update-version",
]
license-files = ["LICENSE.txt"]
include-package-data = false

[tool.setuptools.package-data]
"sage.libs.gap" = ["sage.gaprc"]
"sage.interfaces" = ["sage-maxima.lisp"]
"sage.doctest" = ["tests/*"]
"sage.repl.rich_output" = ["example*"]
sage = [
    "ext_data/*",
    "ext_data/kenzo/*",
    "ext_data/singular/*",
    "ext_data/singular/function_field/*",
    "ext_data/images/*",
    "ext_data/doctest/*",
    "ext_data/doctest/invalid/*",
    "ext_data/gap/*",
    "ext_data/gap/joyner/*",
    "ext_data/mwrank/*",
    "ext_data/notebook-ipython/*",
    "ext_data/nbconvert/*",
    "ext_data/graphs/*",
    "ext_data/pari/*",
    "ext_data/pari/dokchitser/*",
    "ext_data/pari/buzzard/*",
    "ext_data/pari/simon/*",
    "ext_data/magma/*",
    "ext_data/magma/latex/*",
    "ext_data/magma/sage/*",
    "ext_data/valgrind/*",
    "ext_data/threejs/*",
]

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
