SAGE_SPKG_CONFIGURE([givaro], [
    PKG_CHECK_MODULES([GIVARO],
                      [givaro >= 4.1.1],
                      [sage_spkg_install_givaro=no],
                      [sage_spkg_install_givaro=yes])
])
