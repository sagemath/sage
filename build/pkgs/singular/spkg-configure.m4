SAGE_SPKG_CONFIGURE([singular], [
  SAGE_SPKG_DEPCHECK([gmp ntl flint readline mpfr cddlib], [

    AC_PATH_PROG([SINGULAR_BIN], [Singular])
    AC_SUBST([SINGULAR_BIN])
    AS_IF([test -z "${SINGULAR_BIN}"], [sage_spkg_install_singular=yes], [
      dnl Use pkg-config to ensure that Singular is new enough.
      PKG_CHECK_MODULES([SINGULAR], [Singular >= 4.2.1], [
        AC_MSG_CHECKING([whether Singular is built with FLINT])
        AC_COMPILE_IFELSE([
          AC_LANG_PROGRAM([
            #include <singular/singularconfig.h>
            #if !defined(HAVE_FLINT)
            #  error "Need Singular compiled with FLINT"
            #endif
          ], [])
        ], [
          AC_MSG_RESULT([yes])
          AC_MSG_CHECKING([that Singular's help is working])
          AS_IF([test x`printf "system(\"--browser\", \"builtin\"); \n help;" | Singular 2>&1 | grep "error occurred"` = x], [
            AC_MSG_RESULT(yes)
          ], [
            AC_MSG_RESULT(no)
            sage_spkg_install_singular=yes
          ])
        ], [
          AC_MSG_RESULT([no])
          sage_spkg_install_singular=yes
        ])
      ], [
      dnl pkg-config version check failed
      sage_spkg_install_singular=yes
      ])
    ])
  ])
])
