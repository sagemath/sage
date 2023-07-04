SAGE_SPKG_CONFIGURE([mpfr], [
    SAGE_SPKG_DEPCHECK([gmp], [
        AC_CHECK_HEADER(mpfr.h, [], [sage_spkg_install_mpfr=yes])
        dnl mpfr_free_pool appeared in r11922 (Dec 2017) on MPFR svn
        AC_SEARCH_LIBS([mpfr_free_pool], [mpfr], [], [sage_spkg_install_mpfr=yes])
    ])
], [], [], [
    if test x$sage_spkg_install_mpfr = xyes; then
        AC_SUBST(SAGE_MPFR_PREFIX, ['$SAGE_LOCAL'])
        AC_MSG_RESULT([using Sage's mpfr SPKG])
    else
       AC_SUBST(SAGE_MPFR_PREFIX, [''])
       AC_MSG_RESULT([using mpfr library from the system])
    fi
])
