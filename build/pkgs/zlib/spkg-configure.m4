SAGE_SPKG_CONFIGURE([zlib], [
    dnl inflateValidate is needed for libpng at least; checking this ensures
    dnl we have the minimum required zlib version
    AC_CHECK_LIB([z], [inflateValidate], [
        AX_CHECK_ZLIB([], [zlib_cv_libz=no])
    ])
    AS_IF([test "x$zlib_cv_libz" != "xyes"], [sage_spkg_install_zlib=yes])
])
