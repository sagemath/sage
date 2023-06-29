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

protected-files = ["vendor.txt", "__init__.py"]
patches-dir = "vendoring_patches"

[tool.vendoring.transformations]
substitute = [
  {match='sage[.]all', replace='sage.libs.mpmath.hooks'},
  {match="'MPMATH_NOSAGE' not in os.environ and", replace='True or'},
  {match='from mpmath', replace='from sage.libs.mpmath._vendor.mpmath'},
]
