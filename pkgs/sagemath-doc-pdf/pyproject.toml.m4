include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_meson_python
    SPKG_INSTALL_REQUIRES_sage_docbuild
    SPKG_INSTALL_REQUIRES_sphinx
    SPKG_INSTALL_REQUIRES_sphinx_copybutton
    SPKG_INSTALL_REQUIRES_furo
    SPKG_INSTALL_REQUIRES_jupyter_sphinx
    SPKG_INSTALL_REQUIRES_sagelib
]
build-backend = 'mesonpy'

[project]
name = "sagemath-doc-pdf"
version = "10.2"
description = "Sage: Open Source Mathematics Software: Documentation in PDF format"
dependencies = []

include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"
