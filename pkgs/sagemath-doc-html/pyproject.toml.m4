include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_setuptools_wheel
    SPKG_INSTALL_REQUIRES_sage_docbuild
    SPKG_INSTALL_REQUIRES_sphinx
    SPKG_INSTALL_REQUIRES_sphinx_copybutton
    SPKG_INSTALL_REQUIRES_furo
    SPKG_INSTALL_REQUIRES_jupyter_sphinx
    SPKG_INSTALL_REQUIRES_sagelib
]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-doc-html"
description = "Sage: Open Source Mathematics Software: Documentation in HTML format"
dependencies = []
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[tool.setuptools]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
