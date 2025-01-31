SAGE_SPKG_CONFIGURE([libjpeg], [
    PKG_CHECK_MODULES([libjpeg], [libjpeg], [
        AC_SUBST([SAGE_HAVE_LIBJPEG], [1])
    ], [
        sage_spkg_install_libjpeg=yes
    ])
])
