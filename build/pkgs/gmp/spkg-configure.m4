SAGE_SPKG_CONFIGURE([gmp], [
           sage_spkg_install_gmp=no
            AC_CHECK_HEADER(gmp.h, [], [sage_spkg_install_gmp=yes])
            AC_CHECK_HEADER(gmpxx.h, [], [sage_spkg_install_gmp=yes])
            dnl mpn_gcd_11 appeared in GMP 6.2.1
            dnl It is undocumented but is used by Flint when built with default
            dnl flags.
            AC_SEARCH_LIBS([__gmpn_gcd_11], [gmp], [],
                [sage_spkg_install_gmp=yes])
], [], [], [
    if test x$sage_spkg_install_gmp = xyes; then
        AC_SUBST(SAGE_GMP_PREFIX, ['$SAGE_LOCAL'])
        AC_SUBST(SAGE_GMP_INCLUDE, ['$SAGE_LOCAL/include'])
    else
        dnl If found, we want to get the absolute path to where we
        dnl found it for use with some packages (e.g. iml) that need
        dnl this information at configure time
        AX_ABSOLUTE_HEADER([gmp.h])
        if test x$gl_cv_absolute_gmp_h = x; then
            AC_MSG_ERROR(m4_normalize([
                failed to find absolute path to gmp.h despite it being reported
                found
            ]))
        fi
        AC_SUBST(SAGE_GMP_INCLUDE, [`AS_DIRNAME($gl_cv_absolute_gmp_h)`])
        AC_SUBST(SAGE_GMP_PREFIX, [''])
    fi
])
