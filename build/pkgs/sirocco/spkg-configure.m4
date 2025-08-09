SAGE_SPKG_CONFIGURE([sirocco], [
 SAGE_SPKG_DEPCHECK([mpfr], [
  PKG_CHECK_MODULES([libsirocco], [libsirocco >= 2.1.0], [
    AC_MSG_RESULT([found libsirocco >= 2.1.0, will use installed version])
  ],[
    AC_MSG_RESULT([couldn't find libsirocco >= 2.1.0. It will be installed])
    sage_spkg_install_sirocco=yes
  ])
 ])
])
