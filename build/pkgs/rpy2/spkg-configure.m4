SAGE_SPKG_CONFIGURE([rpy2], [
  SAGE_PYTHON_PACKAGE_CHECK([rpy2])
], [dnl REQUIRED-CHECK
  dnl rpy2 is only needed when there is a usable system R
  dnl Check for the R installation and version
  dnl https://rpy2.github.io/doc/v3.4.x/html/overview.html#requirements
  m4_pushdef([SAGE_R_MINVER], ["3.5"])
  PKG_CHECK_MODULES([R], [libR >= SAGE_R_MINVER], [dnl
    AC_PATH_PROG([R_EXECUTABLE], [R])
    AS_IF([test "x$R_EXECUTABLE" = x], [dnl
      AC_MSG_NOTICE([R is not found])
      dnl No R found, so do not require rpy2 package
      AS_VAR_SET([SPKG_REQUIRE], [no])
    ], [dnl Extract R version
      AC_MSG_CHECKING([for version of R executable])
      R_VERSION=$($R_EXECUTABLE --version | sed -n 's/^R version \([[0-9.]]*\).*/\1/p')
      AC_MSG_RESULT([$R_VERSION])
      dnl Extract libR version
      AC_MSG_CHECKING([for version of libR])
      LIBR_VERSION=$(pkg-config --modversion libR)
      AC_MSG_RESULT([$LIBR_VERSION])
      dnl Compare R and libR versions
      AS_IF([test "x$R_VERSION" = "x$LIBR_VERSION"], [dnl
        AC_MSG_NOTICE([R and libR versions match ($R_VERSION)])
        dnl Good system R is found, require rpy2 package
        AS_VAR_SET([SPKG_REQUIRE], [yes])
      ], [dnl R and libR versions do not match
        AC_MSG_NOTICE([R version ($R_VERSION) does not match libR version ($LIBR_VERSION)])
        AS_VAR_SET([SPKG_REQUIRE], [no])
      ])
    ])
  ], [dnl libR not found or outdated
    AS_VAR_SET([SPKG_REQUIRE], [no])
  ])
  m4_popdef([SAGE_R_MINVER])
])
