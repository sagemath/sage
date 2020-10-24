SAGE_SPKG_CONFIGURE([mpc], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_MPFR])
    AC_MSG_CHECKING([installing mpfr? ])
    if test x$sage_spkg_install_mpfr = xyes; then
        AC_MSG_RESULT([yes; install mpc as well])
        sage_spkg_install_mpc=yes
    else
        AC_MSG_RESULT([no])
        AC_CHECK_HEADER(mpc.h, [], [sage_spkg_install_mpc=yes])
        dnl mpc_cmp_abs appeared in MPC 1.1.0
        AC_SEARCH_LIBS([mpc_cmp_abs], [mpc], [], [sage_spkg_install_mpc=yes])
    fi
], [], [], [
    if test x$sage_spkg_install_mpc = xyes; then
        AC_SUBST(SAGE_MPC_PREFIX, ['$SAGE_LOCAL'])
        AC_MSG_RESULT([using Sage's mpc SPKG])
    else
       AC_SUBST(SAGE_MPC_PREFIX, [''])
       AC_MSG_RESULT([using mpc library from the system])
    fi
])
