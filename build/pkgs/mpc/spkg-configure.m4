SAGE_SPKG_CONFIGURE([mpc], [
    SAGE_SPKG_DEPCHECK([mpfr], [
        AC_CHECK_HEADER(mpc.h, [], [sage_spkg_install_mpc=yes])
        dnl mpc_cmp_abs appeared in MPC 1.1.0
        AC_SEARCH_LIBS([mpc_cmp_abs], [mpc], [], [sage_spkg_install_mpc=yes])
    ])
], [], [], [
    if test x$sage_spkg_install_mpc = xyes; then
        AC_SUBST(SAGE_MPC_PREFIX, ['$SAGE_LOCAL'])
        AC_MSG_RESULT([using Sage's mpc SPKG])
    else
       AC_SUBST(SAGE_MPC_PREFIX, [''])
       AC_MSG_RESULT([using mpc library from the system])
    fi
])
