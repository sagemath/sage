SAGE_SPKG_CONFIGURE(
    [pkgconf], [
	AC_CACHE_CHECK([for pkg-config >= 0.29], [ac_cv_path_PKGCONF], [
        AC_PATH_PROGS_FEATURE_CHECK([PKGCONF], [pkg-config], [
            pkgconf_version=`$ac_path_PKGCONF --version 2>&1`
            AS_IF([test -n "$pkgconf_version"], [
                AX_COMPARE_VERSION([$pkgconf_version], [ge], [0.29], [
                    ac_cv_path_PKGCONF="$ac_path_PKGCONF"
                ])
            ])
        ])
    ])
    AS_IF([test -z "$ac_cv_path_PKGCONF"], [
       sage_spkg_install_pkgconf=yes
       AC_SUBST(SAGE_PKG_CONFIG_PATH, [''])
       AC_MSG_RESULT([installing pkgconf spkg])], [
dnl the following as needed as long as Sage creates .pc files during build and/or configure
       AC_SUBST(SAGE_PKG_CONFIG_PATH, ['$SAGE_LOCAL/lib/pkgconfig'])
       AC_MSG_RESULT([using pkg-config from the system])
    ])
])
