SAGE_SPKG_CONFIGURE([zlib], [
    PKG_CHECK_MODULES([zlib], [zlib >= 1.2.11], [zlib_cv_pc=yes], [zlib_cv_pc=no])
    AC_CHECK_LIB([z], [inflateEnd], [zlib_cv_libz=yes], [zlib_cv_libz=no])
    AC_CHECK_HEADER([zlib.h], [zlib_cv_zlib_h=yes], [zlib_cv_zlib_h=no])
    if test "$zlib_cv_pc" = "yes" && test "$zlib_cv_libz" = "yes" && test "$zlib_cv_zlib_h" = "yes"; then
        PKG_CHECK_MODULES([LIBPNG], [libpng >= 1.2], [], [
            dnl inflateValidate is needed for Sage's libpng, newer than 1.2; this ensures
            dnl we have the minimum required for building zlib version
            AC_CHECK_LIB([z], [inflateValidate], [], [sage_spkg_install_zlib=yes])
        ])
    else
        sage_spkg_install_zlib=yes
    fi
])
