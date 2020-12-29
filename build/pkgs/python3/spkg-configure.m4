SAGE_SPKG_CONFIGURE([python3], [
   AC_ARG_WITH([python],
               [AS_HELP_STRING([--with-python=PYTHON3],
                               [Python 3 executable to use for the Sage venv; default: python3])])
   AS_IF([test x"$with_python" = x2], [AC_MSG_ERROR([Sage cannot be built on Python 2. Exiting.])])
   AS_IF([test x"$with_python" = x3], [
       AC_MSG_NOTICE([The meaning of the option --with-python has changed in Sage 9.2. Ignoring.])
       with_python=''
       ])
   AS_IF([test x"$with_python" = x"no"],
         [AC_MSG_ERROR([building Sage --without-python is not supported])])
   ac_path_PYTHON3="$with_python"

   dnl Trac #30559:  Removed the DEPCHECK for sqlite.  We never use libsqlite3 from SPKG for anything
   dnl other than building the python3 SPKG; so our libsqlite3 cannot create shared library conflicts.
   dnl
   dnl However, if we add another package (providing a shared library linked into a Python module)
   dnl that also uses libsqlite3, then we will have to put the DEPCHECK back in.
   SAGE_SPKG_DEPCHECK([libpng bzip2 xz libffi], [
      dnl Check if we can do venv with a system python3
      dnl instead of building our own copy.
      check_modules="sqlite3, ctypes, math, hashlib, crypt, readline, socket, zlib, distutils.core"
      m4_pushdef([MIN_VERSION], [3.6.0])
      m4_pushdef([LT_VERSION],  [3.10.0])
      AC_CACHE_CHECK([for python3 >= ]MIN_VERSION[, < ]LT_VERSION[ with modules $check_modules], [ac_cv_path_PYTHON3], [
        AS_IF([test x"$ac_path_PYTHON3" != x], [dnl checking explicitly specified $with_python
           AC_MSG_RESULT([])
           AC_PATH_PROG([ac_path_PYTHON3], [$ac_path_PYTHON3])
           SAGE_CHECK_PYTHON_FOR_VENV([$ac_path_PYTHON3],
                                           MIN_VERSION, LT_VERSION,
                                           $check_modules, [
                    dnl It is good
                    ac_cv_path_PYTHON3="$ac_path_PYTHON3"
                    ac_path_PYTHON3_found=:
                    AC_MSG_RESULT([yes])
                    dnl introduction for AC_MSG_RESULT printed by AC_CACHE_CHECK
                    AC_MSG_CHECKING([for python3 >= ]MIN_VERSION[, < ]LT_VERSION[ with modules $check_modules])
           ])
           AS_IF([test -z "$ac_cv_path_PYTHON3"], [
               AC_MSG_ERROR([the python3 selected using --with-python=$with_python is not suitable])
           ])
	], [dnl checking the default system python3
           AC_MSG_RESULT([])
           AC_PATH_PROGS_FEATURE_CHECK([PYTHON3], [python3], [
                SAGE_CHECK_PYTHON_FOR_VENV([$ac_path_PYTHON3],
                                           MIN_VERSION, LT_VERSION,
                                           $check_modules, [
                    dnl It is good
                    ac_cv_path_PYTHON3="$ac_path_PYTHON3"
                    ac_path_PYTHON3_found=:
                    AC_MSG_RESULT([yes])
                    dnl introduction for AC_MSG_RESULT printed by AC_CACHE_CHECK
                    AC_MSG_CHECKING([for python3 >= ]MIN_VERSION[, < ]LT_VERSION[ with modules $check_modules])
                ])
            ])
	])
      ])
      AS_IF([test -z "$ac_cv_path_PYTHON3"], [sage_spkg_install_python3=yes])
      m4_popdef([MIN_VERSION])
      m4_popdef([LT_VERSION])
    ])
],, [
    dnl PRE
], [
    dnl POST
    AS_IF([test x$sage_spkg_install_python3 = xno],
          [PYTHON_FOR_VENV="$ac_cv_path_PYTHON3"],
          [SAGE_MACOSX_DEPLOYMENT_TARGET=legacy])
    AC_SUBST([PYTHON_FOR_VENV])
    AC_SUBST([SAGE_MACOSX_DEPLOYMENT_TARGET])

    dnl These temporary directories are created by the check above
    dnl and need to be cleaned up to prevent the "rm -f conftest*"
    dnl (that a bunch of other checks do) from emitting warnings about
    dnl conftest.dir and conftest_venv being directories.
    rm -rf conftest.dir conftest_venv
])
