SAGE_SPKG_CONFIGURE([givaro], [
    PKG_CHECK_MODULES([GIVARO],
                      [givaro >= 4.1.1],dnl The version test is refined in linbox/spkg-configure.m4
                      [sage_spkg_install_givaro=no],
                      [sage_spkg_install_givaro=yes])
])
