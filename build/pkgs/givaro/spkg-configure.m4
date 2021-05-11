SAGE_SPKG_CONFIGURE([givaro], [
    m4_pushdef([SAGE_GIVARO_MINVER],["40101"])
    SAGE_SPKG_DEPCHECK([gmp mpir], [
        AC_PATH_PROG([GIVAROCONFIG], [givaro-config])
        AS_IF([test x$GIVAROCONFIG = x], [
           AC_MSG_NOTICE([givaro-config not found. Installing givaro])
           sage_spkg_install_givaro=yes], [
           AC_MSG_CHECKING([is givaro's version good enough? ])
           givaro_ver=`$GIVAROCONFIG --decimal-version 2>> config.log`
           AX_COMPARE_VERSION([$givaro_ver], [ge], [$SAGE_GIVARO_MINVER], [
               AC_MSG_RESULT([yes. Use system's givaro])], [
               AC_MSG_RESULT([no. Install givaro])
               sage_spkg_install_givaro=yes])
        ])
    ])
    m4_popdef([SAGE_GIVARO_MINVER])
])
