SAGE_SPKG_CONFIGURE([bliss], [
    m4_pushdef([SAGE_BLISS_MINVER],[0.77])
    m4_pushdef([SAGE_BLISS_MAJOR],[0])
    m4_pushdef([SAGE_BLISS_MINOR],[77])
    AC_CHECK_HEADER([bliss/bliss_C.h], [
           AC_SEARCH_LIBS([bliss_new], [bliss], [
             AC_MSG_CHECKING([checking bliss version directly])
             AC_RUN_IFELSE([AC_LANG_PROGRAM([
                      [#include <bliss/defs.hh>
                      ]],[[
                       if (BLISS_VERSION_MAJOR > ]] SAGE_BLISS_MAJOR [[ ) return 0;
                       if (BLISS_VERSION_MAJOR == ]] SAGE_BLISS_MAJOR [[  &&
                           BLISS_VERSION_MINOR >= ]] SAGE_BLISS_MINOR [[ ) return 0;
                       else return 1;
                      ]])],
                     [AC_MSG_RESULT([Good.])],
                     [AC_MSG_RESULT([Too old.])
                      sage_spkg_install_bliss=yes],
                     []) dnl cross-compilation - noop
           ],
              [sage_spkg_install_bliss=yes])
    ], [sage_spkg_install_bliss=yes])
    m4_popdef([SAGE_BLISS_MINVER])
    m4_popdef([SAGE_BLISS_MAJOR])
    m4_popdef([SAGE_BLISS_MINOR])
])

