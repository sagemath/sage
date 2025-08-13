SAGE_SPKG_CONFIGURE([sirocco], [
 SAGE_SPKG_DEPCHECK([mpfr], [
  PKG_CHECK_MODULES([libsirocco], [libsirocco >= 2.1.0], [
    AC_MSG_RESULT([found libsirocco >= 2.1.0, will use installed version])
  ],[
    AC_MSG_RESULT([couldn't find libsirocco >= 2.1.0. Trying to find an instance without pkg-config])
    AC_CHECK_HEADER([sirocco.h], [
         BACKUP_LIBS=${LIBS}
	 LIBS="${LIBS} -lmpfr -lsirocco"
         AC_LINK_IFELSE([
            AC_LANG_PROGRAM(
            [[#include <sirocco.h>
            ]], [[
	      double tcoefs[]={-2.,2.};
	      homotopyPath(2, tcoefs, 1., 0.);
            ]])], [
                AC_MSG_RESULT([yes])
            ], [
                AC_MSG_RESULT([no])
                sage_spkg_install_sirocco=yes
            ], [
                dnl assume that the person running cross-compiling
                dnl knows what they are doing
                AC_MSG_RESULT([yes])
            ]) 
	 LIBS=${BACKUP_LIBS}
  ], [sage_spkg_install_sirocco=yes])
  ])
 ])
])
