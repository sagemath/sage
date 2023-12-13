SAGE_SPKG_CONFIGURE([mpc], [
    SAGE_SPKG_DEPCHECK([mpfr], [
        AC_CHECK_HEADER(mpc.h, [], [sage_spkg_install_mpc=yes])
        dnl gmpy2 2.2 needs MPC >= 1.2.1 according https://github.com/aleaxit/gmpy/blob/master/src/gmpy2.h#L86
        dnl mpc_sum was added in MPC 1.2.0 according to https://www.multiprecision.org/mpc/olds.html
        AC_SEARCH_LIBS([mpc_sum], [mpc], [], [sage_spkg_install_mpc=yes])
    ])
], [], [], [
    if test x$sage_spkg_install_mpc = xyes; then
        AC_SUBST(SAGE_MPC_PREFIX, ['$SAGE_LOCAL'])
    else
       AC_SUBST(SAGE_MPC_PREFIX, [''])
    fi
])
