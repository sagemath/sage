SAGE_SPKG_CONFIGURE([sqlite], [
  m4_pushdef([SAGE_SQLITE3_MIN_VERSION_MAJOR], [3])
  m4_pushdef([SAGE_SQLITE3_MIN_VERSION_MINOR], [8])
  m4_pushdef([SAGE_SQLITE3_MIN_VERSION_MICRO], [7])
  m4_pushdef([SAGE_SQLITE3_MIN_VERSION], [SAGE_SQLITE3_MIN_VERSION_MAJOR.SAGE_SQLITE3_MIN_VERSION_MINOR.SAGE_SQLITE3_MIN_VERSION_MICRO])
  AC_MSG_CHECKING([libsqlite3 >= sqlite3_min_version])
                     dnl https://www.sqlite.org/c3ref/libversion.html
                     dnl https://www.sqlite.org/c3ref/c_source_id.html
                     SQLITE_SAVED_LIBS="$LIBS"
                     LIBS="$LIBS -lsqlite3"
                     AC_RUN_IFELSE([
                       AC_LANG_PROGRAM([[
                                         #include <sqlite3.h>
                                         #include <assert.h>
                                         #include <stdlib.h>
                                         #include <string.h>
                                       ]],
                                       [[
                                         assert( strcmp(sqlite3_libversion(),SQLITE_VERSION)==0 );
                                         if (SQLITE_VERSION_NUMBER < ]]SAGE_SQLITE3_MIN_VERSION_MAJOR[[*1000000 + ]]SAGE_SQLITE3_MIN_VERSION_MINOR[[*1000 + ]]SAGE_SQLITE3_MIN_VERSION_MICRO[[)
                                           exit(1);
                                         else
                                           exit(0);
                                       ]])
                       ],
                       [AC_MSG_RESULT([yes])],
                       [AC_MSG_RESULT([no])
                        LIBS="$SQLITE_SAVED_LIBS"
                        sage_spkg_install_sqlite=yes],
                       [AC_MSG_RESULT([cross compiling. assume yes])])
  m4_popdef([SAGE_SQLITE3_MIN_VERSION_MAJOR])
  m4_popdef([SAGE_SQLITE3_MIN_VERSION_MINOR])
  m4_popdef([SAGE_SQLITE3_MIN_VERSION_MICRO])
  m4_popdef([SAGE_SQLITE3_MIN_VERSION])
], [dnl REQUIRED-CHECK
  AS_CASE([$host],
          [*-*-cygwin*], [
            dnl sqlite SetDllDirectory in sage_ostools.pyx
            sage_require_sqlite=yes
          ], [
            AC_REQUIRE([SAGE_SPKG_CONFIGURE_PYTHON3])
            AS_IF([test x$sage_spkg_install_python3 = xno], [
                sage_require_sqlite=no
            ])
          ])
]
)
