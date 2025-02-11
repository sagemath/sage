SAGE_SPKG_CONFIGURE([fplll], [
  SAGE_SPKG_DEPCHECK([gcc mpfr], [
    dnl If we're using the system mpfr, use pkgconfig to determine
    dnl if there's a usable system copy of fplll. Unless there's
    dnl a system that ships fplll without fplll.pc file, falling
    dnl back to a manual header/library search is pointless.
    dnl
    dnl Issue #31025: FPLLL/FPyLLL make no guarantee regarding compatibility
    dnl other than "whatever versions were released at the same time should work together"
    PKG_CHECK_MODULES([FPLLL],
        [fplll >= 5.5.0],
        [
        AC_MSG_CHECKING([whether BKZ default strategy JSON is installed])
        AC_LANG_PUSH([C++])
        FPLLL_SAVED_LIBS=$LIBS
        LIBS="$LIBS -lfplll"
        AC_RUN_IFELSE([
            AC_LANG_PROGRAM(
            [[#include <fstream>
              #include <fplll/fplll.h>
              #include <fplll/bkz_param.h>
            ]], [[
              std::ifstream fs;
              fs.open(fplll::default_strategy());
              if (fs) return 0;
              return 1;
            ]])], [
                AC_MSG_RESULT([yes])
            ], [
                AC_MSG_RESULT([no])
                sage_spkg_install_fplll=yes
            ], [
                dnl assume that the person running cross-compiling
                dnl knows what they are doing
                AC_MSG_RESULT([yes])
            ])
        LIBS=$FPLLL_SAVED_LIBS
        AC_LANG_POP()
        ],
        [sage_spkg_install_fplll=yes])
  ])
])
