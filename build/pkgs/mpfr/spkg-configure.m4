SAGE_SPKG_CONFIGURE([mpfr], [
    SAGE_SPKG_DEPCHECK([gmp], [
        AC_CHECK_HEADER(mpfr.h, [], [sage_spkg_install_mpfr=yes])
        dnl gmpy2 2.2 needs MPFR >= 4.1.0 according https://github.com/aleaxit/gmpy/blob/master/src/gmpy2.h#L86
        dnl mpfr_cmpabs_ui was added in 4.1.0 according to https://github.com/BrianGladman/mpfr/blob/master/NEWS#L26
        AC_SEARCH_LIBS([mpfr_cmpabs_ui], [mpfr], [], [sage_spkg_install_mpfr=yes])
    ])
], [], [], [
    if test x$sage_spkg_install_mpfr = xyes; then
        AC_SUBST(SAGE_MPFR_PREFIX, ['$SAGE_LOCAL'])
    else
       AC_SUBST(SAGE_MPFR_PREFIX, [''])
    fi
])
