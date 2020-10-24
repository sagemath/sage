SAGE_SPKG_CONFIGURE([libffi], [
    dnl First try checking for libffi with pkg-config
    PKG_CHECK_MODULES([LIBFFI], [libffi], [], [
        dnl Fallback to manually grubbing around for headers and libs
        AC_CHECK_HEADERS([ffi/ffi.h ffi.h], [sage_spkg_install_libffi=no; break], [sage_spkg_install_libffi=yes])
        AC_SEARCH_LIBS([ffi_call], [ffi], [], [sage_spkg_install_libffi=yes])
    ])
])

