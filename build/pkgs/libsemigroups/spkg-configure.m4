SAGE_SPKG_CONFIGURE([libsemigroups], [
    dnl  checking with pkg-config
    PKG_CHECK_MODULES([LIBSEMIGROUPS], [libsemigroups >= 0.6.7], [], [sage_spkg_install_libsemigroups=yes])
])
