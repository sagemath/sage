SAGE_SPKG_CONFIGURE([qhull], [
    m4_pushdef([SAGE_QHULL_MINVER],["8.0.2"])
    AC_PATH_PROG([QHULL], [qhull])
    AS_IF([test x$QHULL = x], [
           AC_MSG_NOTICE([qhull not found. Installing qhull])
           sage_spkg_install_qhull=yes
    ], [
           AC_MSG_CHECKING([whether qhull's version is good enough])
           qhull_ver=$($QHULL -V | cut -d' ' -f2 2>& AS_MESSAGE_LOG_FD)
           AX_COMPARE_VERSION([$qhull_ver], [ge], [$SAGE_QHULL_MINVER], [
               AC_MSG_RESULT([yes])
               AC_CHECK_HEADER([libqhull_r/libqhull_r.h], [
                  AC_SEARCH_LIBS([qh_distplane], [qhull_r], [
                  ], [
                    sage_spkg_install_qhull=yes
                  ])dnl AC_SEARCH_LIBS
               ], [
                  sage_spkg_install_qhull=yes
               ])dnl AC_CHECK_HEADER
           ], [
               sage_spkg_install_qhull=yes
           ])dnl AX_COMPARE_VERSION
    ])dnl IF
    m4_popdef([SAGE_QHULL_MINVER])
])
