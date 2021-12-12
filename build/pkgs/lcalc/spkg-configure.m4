SAGE_SPKG_CONFIGURE([lcalc], [
  SAGE_SPKG_DEPCHECK([pari], [

    dnl Ensure that the "lcalc" program is in the user's executable PATH,
    dnl and that our headers are in his compiler's "include" path.
    AC_PATH_PROG([LCALC_BIN], [lcalc])
    AS_IF([test -z "${LCALC_BIN}"], [sage_spkg_install_lcalc=yes])
    AC_CHECK_HEADER([lcalc/L.h], [], [sage_spkg_install_lcalc=yes])

    dnl Check for the lcalc-2.x API that we now use, ensure that
    dnl that HAVE_LIBPARI is defined in L.h (via lcalc/config.h),
    dnl and check for PRECISION_DOUBLE in the same place.
    AC_MSG_CHECKING([for double-precision libLfunction >= 2.0.0 with pari support])
    AC_LANG_PUSH([C++])
    LCALC_SAVED_LIBS="${LIBS}"
    LIBS="${LIBS} -lLfunction"
    AC_LINK_IFELSE([
      AC_LANG_PROGRAM([[#include <lcalc/L.h>
                        #if !HAVE_LIBPARI
                        #error libLfunction missing PARI support
                        #endif
                        #if !PRECISION_DOUBLE
                        #error libLfunction must use double-precision
                        #endif]],
                      [[initialize_globals();
                        vector<Double> zeros;
                        L_function<int> zeta;
                        zeta.find_zeros(1, 0, 1025, -1, "", &zeros);
                        return 0;]])
    ], [
      AC_MSG_RESULT([found; using lcalc from the system])
    ], [
      AC_MSG_RESULT([not found; installing the lcalc SPKG])
      sage_spkg_install_lcalc=yes
      LIBS="${LCALC_SAVED_LIBS}"
    ])
    AC_LANG_POP([C++])

  ])
])
