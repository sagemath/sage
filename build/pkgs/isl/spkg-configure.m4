SAGE_SPKG_CONFIGURE([isl], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_GMP])
    AC_MSG_CHECKING([installing gmp/mpir? ])
    if test x$sage_spkg_install_mpir = xyes -o x$sage_spkg_install_gmp = xyes; then
        AC_MSG_RESULT([yes; install isl as well])
        sage_spkg_install_isl=yes
    else
        AC_MSG_RESULT([no])
        PKG_CHECK_MODULES([ISL], [isl >= 0.20], [], [
           sage_spkg_install_isl=yes])
    fi
])

