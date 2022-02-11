SAGE_SPKG_CONFIGURE([suitesparse], [
  SAGE_SPKG_DEPCHECK([openblas], [
        AC_SEARCH_LIBS([cholmod_speye], [cholmod], [
         AC_SEARCH_LIBS([umfpack_di_solve], [umfpack], [
          AC_SEARCH_LIBS([SuiteSparse_version], [suitesparseconfig], [
            AC_CHECK_HEADERS([suitesparse/SuiteSparse_config.h SuiteSparse_config.h], [
               AC_CHECK_HEADERS([suitesparse/amd.h amd.h], [sage_spkg_install_suitesparse=no; break],
                                [sage_spkg_install_suitesparse=yes])
               break
               ], [sage_spkg_install_suitesparse=yes])
          ], [sage_spkg_install_suitesparse=yes])
         ], [
         sage_spkg_install_suitesparse=yes])
        ], [
        sage_spkg_install_suitesparse=yes])
      ], [
      sage_spkg_install_suitesparse=yes])
], [], [], [
    AS_IF([test x$sage_spkg_install_suitesparse = xyes], [
       AC_SUBST(SAGE_SUITESPARSE_PREFIX, ['$SAGE_LOCAL'])
    ], [
       AC_SUBST(SAGE_SUITESPARSE_PREFIX, [''])
    ])
])
