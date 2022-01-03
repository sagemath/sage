SAGE_SPKG_CONFIGURE([r], [
    m4_pushdef([SAGE_R_MINVER],["3.4.4"])
    SAGE_SPKG_DEPCHECK([openblas iconv readline bzip2 liblzma pcre curl], [
      AS_CASE([$host],
        [*-*-cygwin*], [
          dnl #29486: rpy2 2.8.x does not build against system R on cygwin.
          sage_spkg_install_r=yes
        ], [
          PKG_CHECK_MODULES([R], [libR >= SAGE_R_MINVER], [
            AC_PATH_PROG([R], [R])
              AS_IF([test "x$R" = x], [
                AC_MSG_NOTICE([R is not found])
                sage_spkg_install_r=yes
                ], [
                dnl TODO: check that versions of R and libR match
                sage_spkg_install_r=no
                ])
         ], [sage_spkg_install_r=yes])
        ])
    ])
    m4_popdef([SAGE_R_MINVER])
])
