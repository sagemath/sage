SAGE_SPKG_CONFIGURE([ntl], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_GMP])
    AC_MSG_CHECKING([installing gmp/mpir? ])
    if test x$sage_spkg_install_mpir = xyes -o x$sage_spkg_install_gmp = xyes; then
        AC_MSG_RESULT([yes; install ntl as well])
        sage_spkg_install_ntl=yes
    else
        AC_MSG_RESULT([no])
    fi

    m4_pushdef(SAGE_NTL_VERSION_MAJOR, [10])
    m4_pushdef(SAGE_NTL_VERSION_MINOR, [3])

    if test x$sage_spkg_install_ntl != xyes; then
        AC_CHECK_HEADER([NTL/ZZ.h], [], [sage_spkg_install_ntl=yes])
        AC_MSG_CHECKING([whether we can link a program using NTL])
        NTL_SAVED_LIBS=$LIBS
        LIBS="$LIBS -lntl"
        AC_LINK_IFELSE([
            AC_LANG_PROGRAM([[#include <NTL/ZZ.h>]],
                            [[NTL::ZZ a;]]
            )], [AC_MSG_RESULT([yes])], [
            AC_MSG_RESULT([no]); sage_spkg_install_ntl=yes
            LIBS=$NTL_SAVED_LIBS
        ])
        AC_MSG_CHECKING([NTL version >= ]SAGE_NTL_VERSION_MAJOR[.]SAGE_NTL_VERSION_MINOR)
        AC_RUN_IFELSE([
            AC_LANG_PROGRAM(
            [[#include <NTL/version.h>
              #include <stdio.h>
            ]], [[
              printf("%s\n", NTL_VERSION);
              if (NTL_MAJOR_VERSION >]] SAGE_NTL_VERSION_MAJOR[[) return 0;
              else if (NTL_MAJOR_VERSION ==]] SAGE_NTL_VERSION_MAJOR[[ &&
                       NTL_MINOR_VERSION >=]] SAGE_NTL_VERSION_MINOR[[) return 0;
              else return 1;
            ]])], [], [sage_spkg_install_ntl=yes])
    fi

    m4_popdef([SAGE_NTL_VERSION_MAJOR])
    m4_popdef([SAGE_NTL_VERSION_MINOR])
], [], [], [
    if test x$sage_spkg_install_ntl = xyes; then
        AC_SUBST(SAGE_NTL_PREFIX, ['$SAGE_LOCAL'])
        AC_MSG_RESULT([using Sage's ntl SPKG])
    else
        AC_SUBST(SAGE_NTL_PREFIX, [''])
        AC_MSG_RESULT([using ntl library from the system])
    fi
])

