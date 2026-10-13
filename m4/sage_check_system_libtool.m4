# System autotools are needed only for the macOS SPKG builds below.
AC_DEFUN([SAGE_CHECK_SYSTEM_LIBTOOL], [
  AS_CASE([$host_os], [darwin*], [
    AS_IF([{ test "$sage_spkg_install_m4ri" = yes ||
             test "$sage_spkg_install_m4rie" = yes ||
             test "$sage_spkg_install_fflas_ffpack" = yes ||
             test "$sage_spkg_install_linbox" = yes; }], [
      AC_PATH_PROG([SAGE_AUTORECONF], [autoreconf])
      AC_PATH_PROG([SAGE_AUTOCONF], [autoconf])
      AC_PATH_PROG([SAGE_AUTOMAKE], [automake])
      AC_PATH_PROG([SAGE_ACLOCAL], [aclocal])
      AC_PATH_PROGS([SAGE_LIBTOOLIZE], [glibtoolize libtoolize])
      AS_IF([test -z "$SAGE_AUTORECONF" || test -z "$SAGE_AUTOCONF" ||
             test -z "$SAGE_AUTOMAKE" || test -z "$SAGE_ACLOCAL" ||
             test -z "$SAGE_LIBTOOLIZE"], [
        AC_MSG_ERROR([Building M4RI, M4RIE, FFLAS-FFPACK or LinBox on macOS requires system autoconf, automake and GNU libtool >= 2.6.2. With Homebrew, run: brew install autoconf automake libtool])
      ])
      AC_MSG_CHECKING([for GNU libtool >= 2.6.2])
      sage_libtool_version=$("$SAGE_LIBTOOLIZE" --version | sed -n '1s/^.*(GNU libtool) //p')
      AS_IF([test -z "$sage_libtool_version"], [
        AC_MSG_ERROR([SAGE_LIBTOOLIZE must name GNU libtoolize (glibtoolize on macOS)])
      ])
      AC_MSG_RESULT([$sage_libtool_version])
      AS_VERSION_COMPARE([$sage_libtool_version], [2.6.2], [
        AC_MSG_ERROR([GNU libtool >= 2.6.2 is required to handle Apple Clang OpenMP flags, including flags from external dependencies. Upgrade libtool and rerun configure.])
      ])
    ])
  ])
  AC_SUBST([SAGE_AUTORECONF])
  AC_SUBST([SAGE_AUTOCONF])
  AC_SUBST([SAGE_AUTOMAKE])
  AC_SUBST([SAGE_ACLOCAL])
  AC_SUBST([SAGE_LIBTOOLIZE])
])
