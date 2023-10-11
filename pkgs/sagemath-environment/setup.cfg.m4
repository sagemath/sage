include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-environment
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: System and software environment
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.9, <3.12
install_requires =

py_modules =
    sage.all__sagemath_environment
    sage.env
    sage.version
    sage.misc.all__sagemath_environment
    sage.misc.package
    sage.misc.package_dir
    sage.misc.temporary_file
    sage.misc.viewer

packages =
    sage.features

scripts =
    # The sage script
    bin/sage
    # Auxiliary scripts for setting up the environment
    bin/sage-env
    bin/sage-num-threads.py
    bin/sage-venv-config
    bin/sage-version.sh
    # Auxiliary script for invoking Python in the Sage environment
    bin/sage-python
    # Not included:
    # - bin/sage-env-config                  -- installed by sage_conf
    # - bin/sage-env-config.in               -- not to be installed
    # - bin/sage-run, bin/sage-runtests, ... -- installed by sagemath-repl
    # - bin/sage-ipython  -- uses sage.repl, so installed by sagemath-repl

[options.extras_require]
# sage.env can optionally use sage_conf
conf = SPKG_INSTALL_REQUIRES_sage_conf

# For "sage --docbuild"
docbuild = SPKG_INSTALL_REQUIRES_sage_docbuild

# For "sage", "sage -t", ...
sage = SPKG_INSTALL_REQUIRES_sagelib

# For "sage --cython"
cython = SPKG_INSTALL_REQUIRES_cython

# For "sage --pytest"
pytest = SPKG_INSTALL_REQUIRES_pytest

# For "sage --rst2ipynb"
rst2ipynb = SPKG_INSTALL_REQUIRES_rst2ipynb

# For "sage --tox"
tox = SPKG_INSTALL_REQUIRES_tox

# For "sage --sws2rst"
sws2rst = SPKG_INSTALL_REQUIRES_sage_sws2rst
