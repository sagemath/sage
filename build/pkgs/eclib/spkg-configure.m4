SAGE_SPKG_CONFIGURE([eclib], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_NTL])
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI])
    AC_MSG_CHECKING([installing ntl or pari? ])
    if test x$sage_spkg_install_ntl = xyes -o x$sage_spkg_install_pari = xyes; then
        AC_MSG_RESULT([yes; install eclib as well])
        sage_spkg_install_eclib=yes
    else
        dnl header types.h appeared in v20180710
        AC_CHECK_HEADER([eclib/types.h], [
          AC_MSG_CHECKING([whether we can link and run a program using eclib])
          ECLIB_SAVED_LIBS=$LIBS
          LIBS="$LIBS -lec"
          AC_RUN_IFELSE([
            AC_LANG_PROGRAM([[#include <eclib/version.h>]
                             [#include <eclib/interface.h>]],
                      [[set_bit_precision(42); /* test for versions >= v20190226 */
                        show_version();
                        return 0;]]
            )], [AC_MSG_RESULT([yes; use eclib from the system])], [
            AC_MSG_RESULT([no; install eclib])
            sage_spkg_install_eclib=yes
            LIBS=$ECLIB_SAVED_LIBS
          ])
        ], [sage_spkg_install_eclib=yes])
    fi
])
