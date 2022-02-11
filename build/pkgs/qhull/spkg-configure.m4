SAGE_SPKG_CONFIGURE([qhull], [
    m4_pushdef([SAGE_QHULL_MINVER],["8.0.2"])
    AC_PATH_PROG([QHULL], [qhull])
    AS_IF([test x$QHULL = x], [
           AC_MSG_NOTICE([qhull not found. Installing qhull])
           sage_spkg_install_qhull=yes
    ], [
           AC_MSG_CHECKING([is qhull's version good enough? ])
           qhull_ver=`$QHULL -V | cut -d' ' -f2 2>> config.log`
           AX_COMPARE_VERSION([$qhull_ver], [ge], [$SAGE_QHULL_MINVER], [
               AC_MSG_RESULT([yes.])
               AC_MSG_CHECKING([is qhull_r library and headers installed? ])
               AC_CHECK_HEADER([libqhull_r/libqhull_r.h], [
                  AC_SEARCH_LIBS([qh_distplane], [qhull_r], [
                     AC_MSG_RESULT([yes. Use system's qhull])
                  ], [
                  AC_MSG_RESULT([no. Install qhull])
                  sage_spkg_install_qhull=yes]) dnl SEARCH_LIBS
               ], [
                  AC_MSG_RESULT([no. Install qhull])
                  sage_spkg_install_qhull=yes
               ]), dnl CHECK_HEADER
           ], [
               AC_MSG_RESULT([no. Install qhull])
               sage_spkg_install_qhull=yes
           ]) dnl AX_COMPARE_VERSION
    ]) dnl IF
    m4_popdef([SAGE_QHULL_MINVER])
])
