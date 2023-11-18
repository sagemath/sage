SAGE_SPKG_CONFIGURE([libjpeg], [
    PKG_CHECK_MODULES([libjpeg], [libjpeg], [], [sage_spkg_install_libjpeg=yes])
])
