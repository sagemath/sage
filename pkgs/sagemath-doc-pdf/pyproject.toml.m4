include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_meson_python
    SPKG_INSTALL_REQUIRES_sage_docbuild
    SPKG_INSTALL_REQUIRES_sphinx
    SPKG_INSTALL_REQUIRES_sphinx_copybutton
    SPKG_INSTALL_REQUIRES_sphinx_inline_tabs
    SPKG_INSTALL_REQUIRES_furo
    SPKG_INSTALL_REQUIRES_jupyter_sphinx
    SPKG_INSTALL_REQUIRES_sagelib
]
build-backend = 'mesonpy'

[project]
name = "sagemath-doc-pdf"
description = "Sage: Open Source Mathematics Software: Documentation in PDF format"
dependencies = []
license = {text = "GNU General Public License (GPL) v2 or later; Creative Commons Attribution-ShareAlike 3.0 Unported"}
authors = [{name = "The Sage Developers", email = "sage-support@googlegroups.com"}]
classifiers = [
    "Development Status :: 6 - Mature",
    "Intended Audience :: Education",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved",
    "Programming Language :: Python :: 3 :: Only",
    "Topic :: Scientific/Engineering :: Mathematics",
]
urls = {Homepage = "https://www.sagemath.org"}
dynamic = ["version"]

[project.readme]
file = "README.rst"
content-type = "text/x-rst"
