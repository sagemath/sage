SAGE_SPKG_CONFIGURE([cbc], [
    SAGE_SPKG_DEPCHECK([openblas], [
        dnl  checking with pkg-config
        PKG_CHECK_MODULES([CBC], [cbc >= 2.9.4], [], [sage_spkg_install_cbc=yes])
    ])
])
