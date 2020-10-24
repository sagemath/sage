SAGE_SPKG_CONFIGURE([planarity], [
     AC_LANG_PUSH([C])
     AC_CHECK_HEADER([planarity/planarity.h], [
        AC_CHECK_LIB([planarity], [gp_InitGraph], [
         AC_MSG_CHECKING([for planarity version 3.0 or later])
         AC_COMPILE_IFELSE(
          [AC_LANG_PROGRAM(
            [[#include <planarity/graphStructures.h>]],
            [[vertexRec v;]
             [v.link[0]=1;]])
          ], [
           AC_MSG_RESULT([yes])
          ], [
           AC_MSG_RESULT([no])
           sage_spkg_install_planarity=yes
        ])
       ], [sage_spkg_install_planarity=yes])
        ], [sage_spkg_install_planarity=yes])
     AC_LANG_POP()
])
