SAGE_SPKG_CONFIGURE([isl], [
    SAGE_SPKG_DEPCHECK([gmp], [
        PKG_CHECK_MODULES([ISL], [isl >= 0.20], [], [
           sage_spkg_install_isl=yes])
    ])
])

