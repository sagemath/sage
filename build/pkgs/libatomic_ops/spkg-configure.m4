SAGE_SPKG_CONFIGURE([libatomic_ops], [
  PKG_CHECK_MODULES([LIBATOMIC_OPS],
                    [atomic_ops >= 7.6.2],
                    [],
                    [sage_spkg_install_libatomic_ops=yes])
])
