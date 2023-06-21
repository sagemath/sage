include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_setuptools
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_sagemath_categories
]
build-backend = "setuptools.build_meta"


[tool.vendoring]
# Following example at https://github.com/pypa/pip/blob/main/pyproject.toml#L30

destination = "sage/libs/mpmath/_vendor/"
requirements = "vendor.txt"
namespace = "sage.libs.mpmath._vendor"

protected-files = ["vendor.txt", "all__sagemath_mpmath.py"]
patches-dir = "vendoring_patches"

[tool.vendoring.transformations]
substitute = [
  {match='sage[.]all', replace='sage.all__sagemath_categories'},
]
## try running 'vendoring' as part of building the sdist!  PIP_FIND_LINKS=upstream
