SAGE_SPKG_CONFIGURE([mpfi], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_MPFR])
    AC_MSG_CHECKING([installing mpfr? ])
    if test x$sage_spkg_install_mpfr = xyes; then
        AC_MSG_RESULT([yes; install mpfi as well])
        sage_spkg_install_mpfi=yes
    else
        AC_MSG_RESULT([no])
    fi

    m4_pushdef(SAGE_MPFI_VERSION_MAJOR, [1])
    m4_pushdef(SAGE_MPFI_VERSION_MINOR, [5])

    if test x$sage_spkg_install_mpfi != xyes; then
        AC_CHECK_HEADER([mpfi.h], [], [sage_spkg_install_mpfi=yes])
        AC_SEARCH_LIBS([mpfi_diam_abs], [mpfi], [
          AC_LANG_PUSH(C)
          AC_MSG_CHECKING([MPFI version >= ]SAGE_MPFI_VERSION_MAJOR[.]SAGE_MPFI_VERSION_MINOR)
          AC_RUN_IFELSE([
            AC_LANG_PROGRAM(
            [[#include <mpfi.h>
              #include <stdio.h>
            ]], [[
              fprintf(stderr, "%s\n", MPFI_VERSION_STRING);
              if (MPFI_VERSION_MAJOR >]] SAGE_MPFI_VERSION_MAJOR[[) return 0;
              else if (MPFI_VERSION_MAJOR ==]] SAGE_MPFI_VERSION_MAJOR[[ &&
                       MPFI_VERSION_MINOR >=]] SAGE_MPFI_VERSION_MINOR[[) return 0;
              else return 1;
            ]])], [AC_MSG_RESULT([yes])], [
                   AC_MSG_RESULT([no])
                   sage_spkg_install_mpfi=yes])
        AC_LANG_POP(C)], [sage_spkg_install_mpfi=yes])
    fi

    m4_popdef([SAGE_MPFI_VERSION_MAJOR])
    m4_popdef([SAGE_MPFI_VERSION_MINOR])
])

