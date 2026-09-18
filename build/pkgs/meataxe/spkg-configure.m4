SAGE_SPKG_CONFIGURE([meataxe],
 [PKG_CHECK_MODULES([MEATAXE], [libmtx], [], [sage_spkg_install_meataxe=yes])
 ])
