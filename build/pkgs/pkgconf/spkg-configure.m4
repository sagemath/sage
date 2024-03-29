SAGE_SPKG_CONFIGURE(
    [pkgconf], [
    dnl Check is in configure.ac
    AS_IF([test -z "$PKG_CONFIG"], [
       sage_spkg_install_pkgconf=yes
       AC_SUBST(SAGE_PKG_CONFIG_PATH, [''])
    ], [
dnl the following as needed as long as Sage creates .pc files during build and/or configure
       AC_SUBST(SAGE_PKG_CONFIG_PATH, ['$SAGE_LOCAL/lib/pkgconfig'])
    ])
])
