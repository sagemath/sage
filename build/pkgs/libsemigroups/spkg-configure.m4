SAGE_SPKG_CONFIGURE([libsemigroups], [
    dnl  checking with pkg-config
    PKG_CHECK_MODULES([LIBSEMIGROUPS], [libsemigroups >= 3.5.5], [], [sage_spkg_install_libsemigroups=yes])
])
