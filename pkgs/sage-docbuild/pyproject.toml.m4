include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
requires = [
    SPKG_INSTALL_REQUIRES_setuptools
]
build-backend = "setuptools.build_meta"

[project]
name = "sage-docbuild"
description = "Sage: Open Source Mathematics Software: Build system of the Sage documentation"
dependencies = [
    SPKG_INSTALL_REQUIRES_sphinx
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[tool.setuptools]
packages = [
    "sage_docbuild",
    "sage_docbuild.ext",
]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
