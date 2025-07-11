SAGE_SPKG_CONFIGURE([planarity], [
     AC_LANG_PUSH([C])
     AC_CHECK_LIB([planarity], [gp_InitGraph], [
       AC_CHECK_HEADERS([planarity/graphLib.h planarity/graph.h], [
       ], [sage_spkg_install_planarity=yes])dnl have not found planarity 3.* or newer headers
     ], [sage_spkg_install_planarity=yes])dnl have not found planarity dylib
     AC_LANG_POP()
])
