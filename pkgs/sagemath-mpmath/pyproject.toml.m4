[build-system]
# Minimum requirements for the build system to execute.
requires = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        sage_setup     \
        sagemath_environment \
        sagemath_categories \
        cython         \
        cysignals      \
                    ')]
build-backend = "setuptools.build_meta"


[tool.vendoring]
# Following example at https://github.com/pypa/pip/blob/main/pyproject.toml#L30

destination = "sage/libs/mpmath/_vendor/"
requirements = "vendor.txt"
namespace = "sage.libs.mpmath._vendor"

protected-files = ["all.py", "vendor.txt", "nodoctest.py"]
patches-dir = "vendoring_patches"

[tool.vendoring.transformations]
substitute = [
  {match='sage[.]all', replace='sage.libs.mpmath.hooks'},
  {match="(?ms)if 'MPMATH_NOGMPY'.*?except:\\s*pass", replace=''},
  {match="(?ms)if [(]*'MPMATH_NOSAGE'.*?:(.*?)except:", replace='if True:\1except ():'},
  {match='from mpmath', replace='from sage.libs.mpmath._vendor.mpmath'},
]
