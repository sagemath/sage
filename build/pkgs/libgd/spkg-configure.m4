SAGE_SPKG_CONFIGURE([libgd], [
    dnl Issue #31624: Avoid C++ ABI issues
    SAGE_SPKG_DEPCHECK([gcc libpng freetype], [
        PKG_CHECK_MODULES([LIBGD], [gdlib >= 2.1], [], [sage_spkg_install_libgd=yes])
    ])
])
