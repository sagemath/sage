SAGE_SPKG_CONFIGURE([mpfr], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_GMP])
    AC_MSG_CHECKING([installing gmp/mpir? ])
    if test x$sage_spkg_install_mpir = xyes -o x$sage_spkg_install_gmp = xyes; then
        AC_MSG_RESULT([yes; install mpfr as well])
        sage_spkg_install_mpfr=yes
    else
        AC_MSG_RESULT([no])
        AC_CHECK_HEADER(mpfr.h, [], [sage_spkg_install_mpfr=yes])
        dnl mpfr_free_pool appeared in r11922 (Dec 2017) on MPFR svn
        AC_SEARCH_LIBS([mpfr_free_pool], [mpfr], [], [sage_spkg_install_mpfr=yes])
    fi
], [], [], [
    if test x$sage_spkg_install_mpfr = xyes; then
        AC_SUBST(SAGE_MPFR_PREFIX, ['$SAGE_LOCAL'])
        AC_MSG_RESULT([using Sage's mpfr SPKG])
    else
       AC_SUBST(SAGE_MPFR_PREFIX, [''])
       AC_MSG_RESULT([using mpfr library from the system])
    fi
])
