SAGE_SPKG_CONFIGURE([openssl], [
  AX_CHECK_OPENSSL([], [sage_spkg_install_openssl=yes])
])
