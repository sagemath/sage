SAGE_SPKG_CONFIGURE([symmetrica], [
    AC_LANG_PUSH(C)
    AC_CHECK_HEADER([symmetrica/def.h], [
dnl check for one of its many functions 
     AC_SEARCH_LIBS([zykelind_tetraeder_edges_extended], [symmetrica], [
     AC_MSG_CHECKING([whether we have a properly patched Symmetrica version])
     AC_RUN_IFELSE([AC_LANG_PROGRAM([dnl this crashes on unpatched Symmetrica
         [#include "symmetrica/def.h"]
         [#include "symmetrica/macro.h"]],
         [[OP b,n;]
          [int i;]
          [anfang();]
          [for (i=1; i<6; i++) {]
            [n = callocobject();]
            [b = callocobject();]
            [M_I_I(i, n);]
            [kostka_tafel(n, b);]
            [fprintln(stderr, b);]
            [freeall(n);]
            [freeall(b);}]
          [ende();]])],
         [AC_MSG_RESULT([appears to be a well-patched version.])],
         [AC_MSG_RESULT([buggy version. Sage will build its own.])
                         sage_spkg_install_symmetrica=yes],
         [AC_MSG_RESULT([cross compiling. Assume not buggy.])])
     ], [sage_spkg_install_symmetrica=yes])
    ], [sage_spkg_install_symmetrica=yes])
    AC_LANG_POP(C)
])
