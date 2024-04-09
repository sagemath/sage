SAGE_SPKG_CONFIGURE([ecm], [dnl CHECK - test whether the package is already installed
    m4_pushdef([SAGE_ECM_MINVER],[7.0.4])
    SAGE_SPKG_DEPCHECK([gmp], [
        AC_CHECK_HEADER(ecm.h, [
            AX_ABSOLUTE_HEADER([ecm.h])
            if test x$gl_cv_absolute_ecm_h = x; then
                AC_MSG_ERROR(m4_normalize([
                    failed to find absolute path to ecm.h despite it being reported found
                ]))
                sage_spkg_install_ecm=yes
            else
                dnl check that the version is at least $SAGE_ECM_MINVER
                ecm_version=`grep ECM_VERSION $gl_cv_absolute_ecm_h |
                  $SED -n -e 's/\#define ECM_VERSION "*\([[0-9]]*\.[[0-9]]*\.[[0-9]]*\)"/\1/p'`
                AS_IF([test -n "$ecm_version"], [
                    AX_COMPARE_VERSION([$ecm_version], [ge], [$SAGE_ECM_MINVER], [
                        ac_cv_ECM="$ecm_version"
                        AC_SEARCH_LIBS([ecm_factor], [ecm], [], [sage_spkg_install_ecm=yes])
                    ])
                ])
                AC_PATH_PROGS([ECMBIN], [ecm gmp-ecm])
                AS_IF([test x$ECMBIN != x], [
                    ecmbin_version=`echo 121 | $ECMBIN 4 | grep ^GMP |
                      $SED -n -e 's/GMP\-ECM \([[0-9]]*\.[[0-9]]*\.[[0-9]]*\).*/\1/p'`
                    AS_IF([test -n "$ecmbin_version"], [
                        AX_COMPARE_VERSION([$ecmbin_version], [ge], [$SAGE_ECM_MINVER], [
                          ac_cv_ECMBIN="$ecmbin_version"
                        ])
                    ])
		], [ECMBIN=ecm])
            fi
            AS_IF([test -z "$ac_cv_ECM"], [sage_spkg_install_ecm=yes])
            AS_IF([test -z "$ac_cv_ECMBIN"], [sage_spkg_install_ecm=yes])
        ], [sage_spkg_install_ecm=yes])
    ])
    m4_popdef([SAGE_ECM_MINVER])
], [dnl REQUIRED-check is empty
], [dnl PRE - always perform
    ECMBIN=ecm
], [dnl POST - always perform
    AC_SUBST([SAGE_ECMBIN], ["$ECMBIN"])
])
