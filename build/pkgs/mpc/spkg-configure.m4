SAGE_SPKG_CONFIGURE([mpc], [
    SAGE_SPKG_DEPCHECK([mpfr], [
        dnl gmpy2 2.2 needs MPC >= 1.2.1 according https://github.com/aleaxit/gmpy/blob/master/src/gmpy2.h#L86
        AC_MSG_CHECKING([for MPC >= 1.2.1])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <mpc.h>]], [[
#if (!defined(MPC_VERSION) || (MPC_VERSION < MPC_VERSION_NUM(1,2,1)))
#  error "Sage requires MPC 1.2.1 or later (for gmpy2 2.2)."
#endif
        ]])], [
            AC_MSG_RESULT([yes])
            dnl mpc_sum was added in MPC 1.2.0 according to https://www.multiprecision.org/mpc/olds.html
            AC_SEARCH_LIBS([mpc_sum], [mpc], [], [sage_spkg_install_mpc=yes])
        ], [
            AC_MSG_RESULT([no])
            sage_spkg_install_mpc=yes
        ])
    ])
], [], [], [
    if test x$sage_spkg_install_mpc = xyes; then
        AC_SUBST(SAGE_MPC_PREFIX, ['$SAGE_LOCAL'])
    else
       AC_SUBST(SAGE_MPC_PREFIX, [''])
    fi
])
