include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_setuptools
    SPKG_INSTALL_REQUIRES_wheel
]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-environment"
description = "Sage: Open Source Mathematics Software: System and software environment"
dependencies = []
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.optional-dependencies]
conf      = [SPKG_INSTALL_REQUIRES_sage_conf]           # sage.env can optionally use sage_conf
docbuild  = [SPKG_INSTALL_REQUIRES_sage_docbuild]       # For "sage --docbuild"
sage      = [SPKG_INSTALL_REQUIRES_sagelib]             # For "sage", "sage -t", ...
cython    = [SPKG_INSTALL_REQUIRES_cython]              # For "sage --cython"
pytest    = [SPKG_INSTALL_REQUIRES_pytest]              # For "sage --pytest"
rst2ipynb = [SPKG_INSTALL_REQUIRES_rst2ipynb]           # For "sage --rst2ipynb"
tox       = [SPKG_INSTALL_REQUIRES_tox]                 # For "sage --tox"
sws2rst   = [SPKG_INSTALL_REQUIRES_sage_sws2rst]        # For "sage --sws2rst"

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[tool.setuptools]
py-modules = [
    "sage.all__sagemath_environment",
    "sage.env",
    "sage.version",
    "sage.misc.all__sagemath_environment",
    "sage.misc.package",
    "sage.misc.package_dir",
    "sage.misc.temporary_file",
    "sage.misc.viewer",
]
packages = ["sage.features"]
script-files = [
    # The sage script
    "bin/sage",
    # Auxiliary scripts for setting up the environment
    "bin/sage-env",
    "bin/sage-num-threads.py",
    "bin/sage-venv-config",
    "bin/sage-version.sh",
    # Auxiliary script for invoking Python in the Sage environment
    "bin/sage-python",
    # Not included:
    # - bin/sage-env-config                  -- installed by sage_conf
    # - bin/sage-env-config.in               -- not to be installed
    # - bin/sage-run, bin/sage-runtests, ... -- installed by sagemath-repl
    # - bin/sage-ipython  -- uses sage.repl, so installed by sagemath-repl
]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
