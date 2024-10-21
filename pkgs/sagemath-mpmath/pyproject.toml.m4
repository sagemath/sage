include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_meson_python
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_pkgconfig
]
build-backend = "mesonpy"

[project]
name = "sagemath-mpmath"
description = "Sage: Open Source Mathematics Software: Vendored copy of mpmath using the Sage backend"
dependencies = [
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_sagemath_categories
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[tool.setuptools]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}

[tool.vendoring]
# Following example at https://github.com/pypa/pip/blob/main/pyproject.toml#L30

destination = "sage/libs/mpmath/_vendor/"
requirements = "vendor.txt"
namespace = "sage.libs.mpmath._vendor"

protected-files = ["vendor.txt", "__init__.py", "nodoctest.py"]
patches-dir = "vendoring_patches"

[tool.vendoring.transformations]
substitute = [
  {match='sage[.]all', replace='sage.libs.mpmath.hooks'},
  {match="'MPMATH_NOSAGE' not in os.environ and", replace='True or'},
  {match='from mpmath', replace='from sage.libs.mpmath._vendor.mpmath'},
]
