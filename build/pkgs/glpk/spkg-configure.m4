SAGE_SPKG_CONFIGURE([glpk], [
    m4_pushdef([SAGE_GLPK_MINVER],["4.63"])
    SAGE_SPKG_DEPCHECK([gmp zlib], [
        AC_PATH_PROG([GLPSOL], [glpsol])
        AS_IF([test x$GLPSOL = x], [
           AC_MSG_NOTICE([glpsol not found. Installing glpk])
           sage_spkg_install_glpk=yes], [
           glpk_ver=`$GLPSOL --version | grep ^GLPSOL | $SED -e 's/GLPSOL.*ver, v//g' 2>> config.log`
           AX_COMPARE_VERSION([$glpk_ver], [ge], [$SAGE_GLPK_MINVER], [
               AC_CHECK_HEADER([glpk.h], [], [sage_spkg_install_glpk=yes])
               AC_SEARCH_LIBS([glp_config], [glpk], [], [sage_spkg_install_glpk=yes])dnl end-AC_SEARCH_LIBS
           ], [sage_spkg_install_glpk=yes])dnl end-AX_COMPARE_VERSION
         ])dnl end-AS_IF
    ])
    m4_popdef([SAGE_GLPK_MINVER])
], [], [], [
    AS_IF([test x$sage_spkg_install_glpk = xyes], [
        AC_SUBST(SAGE_GLPK_PREFIX, ['$SAGE_LOCAL'])
      ], [
        AC_SUBST(SAGE_GLPK_PREFIX, [''])
      ])
])
