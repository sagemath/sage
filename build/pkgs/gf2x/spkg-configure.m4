SAGE_SPKG_CONFIGURE([gf2x], [
   PKG_CHECK_MODULES([GF2X], [gf2x >= 1.2], [], [dnl for 1.2 we have to do more work
    AC_CHECK_HEADER(gf2x.h, [], [sage_spkg_install_gf2x=yes])
dnl gf2x_mul_r appeared in version 1.1 of GF2X
    AC_SEARCH_LIBS([gf2x_mul_r], [gf2x], [
        AC_MSG_CHECKING([for GF2X 1.2])
        AC_LANG_PUSH(C)
        old_CFLAGS="$CFLAGS"
        CFLAGS="$CFLAGS -Werror=incompatible-pointer-types"
        AC_COMPILE_IFELSE(
          [AC_LANG_PROGRAM(
            [[#include <gf2x.h>]],
            [[const void (*fptr)(unsigned long *, const unsigned long *, unsigned long,]
             [                   const unsigned long *, unsigned long, gf2x_mul_pool_t);]
             [fptr = gf2x_mul_r;]])
          ],[
           AC_MSG_RESULT([yes])
           sage_spkg_install_gf2x=no
          ],[
           AC_MSG_RESULT([no])
           sage_spkg_install_gf2x=yes
        ])
        CFLAGS="$old_CFLAGS"
        AC_LANG_POP(C)],
      [sage_spkg_install_gf2x=yes])
     ])
])
