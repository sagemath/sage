SAGE_SPKG_CONFIGURE([primesieve], [
    m4_pushdef([SAGE_PRIMESIEVE_MINVER],[8.0])
    dnl Checking for primesieve with pkg-config
    PKG_CHECK_MODULES([PRIMESIEVE], [primesieve >= SAGE_PRIMESIEVE_MINVER], [ ], [
        sage_spkg_install_primesieve=yes])
    m4_popdef([SAGE_PRIMESIEVE_MINVER])
])

