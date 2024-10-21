include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
requires = [
    SPKG_INSTALL_REQUIRES_setuptools
]
build-backend = "setuptools.build_meta"

[project]
name = "sage-setup"
description = "Sage: Open Source Mathematics Software: Build system of the Sage library"
dependencies = []
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
autogen = [
    SPKG_INSTALL_REQUIRES_jinja2
]

[project.scripts]
sage-generate-meson = "sage_setup.autogen.meson:generate_meson"

[tool.setuptools]
packages = [
    "sage_setup",
    "sage_setup.autogen",
    "sage_setup.autogen.interpreters",
    "sage_setup.autogen.interpreters.specs",
    "sage_setup.autogen.meson",
    "sage_setup.command",
]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
