#
# SYNOPSIS
#
#   SAGE_PYTHON_PACKAGE_CHECK(package)
#
# DESCRIPTION
#
#   Determine if the system copy of a python package can be used by sage.
#
#   This macro uses setuptools.version's pkg_resources to check that the
#   "install-requires.txt" file for the named package is satisfied, and
#   it can typically fail in four ways:
#
#     1. If --enable-system-site-packages was not passed to ./configure,
#
#     2. If we are not using the system python (no $PYTHON_FOR_VENV),
#
#     3. If setuptools is not available to the system python,
#
#     4. If the contents of install-requires.txt are not met (wrong
#        version, no version, etc.) by the system python.
#
#   In any of those cases, we set sage_spkg_install_$package to "yes"
#   so that the corresponding SPKG is installed. Otherwise, we do
#   nothing, since the default value of sage_spkg_install_$package
#   is "no" (to use the system copy).
#
#   The SAGE_SPKG_CONFIGURE_PYTHON3() macro is AC_REQUIRE'd to ensure
#   that $PYTHON_FOR_VENV is available, if it is going to be available.
#

AC_DEFUN([SAGE_PYTHON_PACKAGE_CHECK], [
  AS_IF([test "${enable_system_site_packages}" = "yes"], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PYTHON3])

    dnl strip all comments from install-requires.txt; this should leave
    dnl only a single line containing the version specification for this
    dnl package.
    SAGE_PKG_VERSPEC=$(sed '/^#/d' "./build/pkgs/$1/install-requires.txt")
    AC_MSG_CHECKING([for python package $1 ("${SAGE_PKG_VERSPEC}")])
    AS_IF(
      ["${PYTHON_FOR_VENV}" -c "from setuptools.version import pkg_resources; pkg_resources.require('${SAGE_PKG_VERSPEC}'.splitlines())"],
      [AC_MSG_RESULT(yes)],
      [AC_MSG_RESULT(no); sage_spkg_install_$1=yes]
    )
  ], [
    sage_spkg_install_$1=yes
  ])
])
