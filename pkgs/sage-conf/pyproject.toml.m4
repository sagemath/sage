include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
requires = [
    SPKG_INSTALL_REQUIRES_setuptools
]
build-backend = "setuptools.build_meta"

[project]
name = "sage-conf"
description = "Sage: Open Source Mathematics Software: Configuration module for the SageMath library"
readme = "README.rst"
dnl Not including the standard metadata from pyproject_toml_metadata.m4
dnl because sage-conf is GPL v3+.
license = {text = "GNU General Public License (GPL) v3 or later"}
authors = [{name = "The Sage Developers", email = "sage-support@googlegroups.com"}]
classifiers = [
    "Development Status :: 6 - Mature",
    "Intended Audience :: Education",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
    "Operating System :: POSIX",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3 :: Only",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: Implementation :: CPython",
    "Topic :: Scientific/Engineering :: Mathematics",
]
urls = {Homepage = "https://www.sagemath.org"}
requires-python = ">=3.9, <3.12"
dynamic = ["version"]

[project.scripts]
sage-config = "sage_conf:_main"

[tool.setuptools]
packages = ["_sage_conf"]
py-modules = ["sage_conf"]
script-files = ["bin/sage-env-config"]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
