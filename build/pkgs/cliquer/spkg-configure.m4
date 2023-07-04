SAGE_SPKG_CONFIGURE([cliquer], [
    AC_SEARCH_LIBS([clique_unweighted_max_weight], [cliquer], [
      AC_LANG_PUSH([C])
      AC_CHECK_HEADER([cliquer/cliquer.h], [
         AC_PATH_PROG([CL], [cl])
         AS_IF([test x$CL = x], [
           AC_MSG_NOTICE([cl (cliquer's CLI) is not found.])
           AC_MSG_NOTICE([No action, as Sage does not need it. cl might be named cliquer on this system.])
         ])
      ], [sage_spkg_install_cliquer=yes])
      AC_LANG_POP()
    ], [sage_spkg_install_cliquer=yes])
])
