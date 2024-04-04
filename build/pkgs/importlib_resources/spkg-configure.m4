SAGE_SPKG_CONFIGURE([importlib_resources], [
  SAGE_PYTHON_PACKAGE_CHECK([importlib_resources])
],[
  # Three of our python packages are backport packages providing
  # python-3.11 features (see coding_in_python.rst):
  #
  #   * importlib_metadata
  #   * importlib_resources
  #   * typing_extensions
  #
  # These packages are therefore not needed with >=python-3.11. Here
  # we test for a python minor version component greater than or equal
  # to 11, and mark this package as "not required" if we succeed.
  AC_MSG_CHECKING([for >=python-3.11])

  # Keep in mind that False (~ zero) in python is success in the shell
  AS_IF(["${PYTHON_FOR_VENV}" -c "import sys; sys.exit(sys.version_info.minor < 11)"],[
    AC_MSG_RESULT([yes])
    sage_require_importlib_resources="no"
  ],[
    AC_MSG_RESULT([no])
  ])
])
