SAGE_SPKG_CONFIGURE([importlib_metadata], [dnl
  SAGE_PYTHON_PACKAGE_CHECK([importlib_metadata])
],[
  dnl Three of our python packages are backport packages providing
  dnl python-3.11 features (see coding_in_python.rst):
  dnl
  dnl   * importlib_metadata
  dnl   * importlib_resources
  dnl   * typing_extensions
  dnl
  dnl These packages are therefore not needed with >=python-3.11. Here
  dnl we test for a python minor version component greater than or equal
  dnl to 11, and mark this package as "not required" if we succeed.
  AS_IF([test -z "${PYTHON_FOR_VENV}"], [dnl
    dnl Python from our SPKG is new enough, no need for the backport package.
    sage_require_importlib_metadata="no"
  ], [dnl
    AC_MSG_CHECKING([for >=python-3.11])
    dnl Keep in mind that False (~ zero) in python is success in the shell
    AS_IF(["${PYTHON_FOR_VENV}" -c "import sys; sys.exit(sys.version_info.minor < 11)"], [dnl
      AC_MSG_RESULT([yes])
      sage_require_importlib_metadata="no"
    ],[dnl
      AC_MSG_RESULT([no])
    ])
  ])
])
