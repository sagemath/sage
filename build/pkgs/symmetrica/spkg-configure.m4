SAGE_SPKG_CONFIGURE([symmetrica], [
    AC_LANG_PUSH(C)
    AC_CHECK_HEADER([symmetrica/def.h], [
dnl check for one of its many functions 
     AC_SEARCH_LIBS([zykelind_tetraeder_edges_extended], [symmetrica], [
     AC_MSG_CHECKING([that we have a properly patched Symmetrica version... ])
     AC_RUN_IFELSE([AC_LANG_PROGRAM([dnl this crashes on unpatched Symmetrica
	 [#include "symmetrica/def.h"]
         [#include "symmetrica/macro.h"]],
	 [[OP b,n;]
          [anfang();]
          [n = callocobject();]
	  [b = callocobject();]
          [sscan_integer("4",n);]
          [kostka_tafel(n, b);]
          [println(b);]
          [ende();]])],
         [AC_MSG_RESULT([appears to be a well-patched version.])],
	 [AC_MSG_RESULT([buggy version. Sage will buld its own.])
		         sage_spkg_install_symmetrica=yes])
     ], [sage_spkg_install_symmetrica=yes])
    ], [sage_spkg_install_symmetrica=yes])
    AC_LANG_POP(C)
])
