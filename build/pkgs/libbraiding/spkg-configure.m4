SAGE_SPKG_CONFIGURE([libbraiding], [
  PKG_CHECK_MODULES([libbraiding], [libbraiding >= 1.3.1], [
    AC_MSG_RESULT([found libbraiding >= 1.3.1, will use installed version])
    sage_spkg_install_libbraiding=no
  ],[
    AC_MSG_RESULT([couldn't find libbraiding >= 1.3.1. It will be installed])
    sage_spkg_install_libbraiding=yes
  ])
])
