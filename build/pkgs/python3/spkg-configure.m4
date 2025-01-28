SAGE_SPKG_CONFIGURE([python3], [
   m4_pushdef([MIN_VERSION],               [3.9.0])
   m4_pushdef([MIN_NONDEPRECATED_VERSION], [3.9.0])
   m4_pushdef([LT_STABLE_VERSION],         [3.13.0])
   m4_pushdef([LT_VERSION],                [3.13.0])
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

   dnl Issue #30559:  Removed the DEPCHECK for sqlite.  We never use libsqlite3 from SPKG for anything
   dnl other than building the python3 SPKG; so our libsqlite3 cannot create shared library conflicts.
   dnl
   dnl However, if we add another package (providing a shared library linked into a Python module)
   dnl that also uses libsqlite3, then we will have to put the DEPCHECK back in.
   SAGE_SPKG_DEPCHECK([bzip2 liblzma libffi zlib], [
      dnl Check if we can do venv with a system python3
      dnl instead of building our own copy.
      dnl  Issue #31160: We no longer check for readline here.
      check_modules="sqlite3, ctypes, math, hashlib, socket, zlib, ssl, ensurepip"
      AC_CACHE_CHECK([for python3 >= ]MIN_VERSION[, < ]LT_VERSION[ with modules $check_modules and setuptools/distutils], [ac_cv_path_PYTHON3], [
        AS_IF([test x"$ac_path_PYTHON3" != x], [dnl checking explicitly specified $with_python
           AC_MSG_RESULT([])
           AC_PATH_PROG([ac_path_PYTHON3], [$ac_path_PYTHON3])
           SAGE_CHECK_PYTHON_FOR_VENV([$ac_path_PYTHON3],
                                           MIN_VERSION, LT_VERSION,
                                           $check_modules, [
                    AS_IF([[conftest_venv/bin/python3 -m sysconfig | grep '^\sw*\(C\|LD\)FLAGS *=.*[" ]-[IL] *[^.]' ]] [>& AS_MESSAGE_LOG_FD 2>&1 ], [
                        AC_MSG_WARN([this is a misconfigured Python whose sysconfig compiler/linker flags contain -I or -L options, which may cause wrong versions of libraries to leak into the build of Python packages - see https://github.com/sagemath/sage/issues/31132])
                    ])
                    dnl It is good
                    ac_cv_path_PYTHON3="$ac_path_PYTHON3"
                    ac_path_PYTHON3_found=:
                    AC_MSG_RESULT([yes])
                    dnl introduction for AC_MSG_RESULT printed by AC_CACHE_CHECK
                    AC_MSG_CHECKING([for python3 >= ]MIN_VERSION[, < ]LT_VERSION[ with modules $check_modules and setuptools/distutils])
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
                    AS_IF([[conftest_venv/bin/python3 -m sysconfig | grep '^\sw*\(C\|LD\)FLAGS *=.*[" ]-[IL] *[^.]' ]] [>& AS_MESSAGE_LOG_FD 2>&1 ], [
                        AC_MSG_RESULT([no, this is a misconfigured Python whose sysconfig compiler/linker flags contain -I or -L options, which may cause wrong versions of libraries to leak into the build of Python packages - see https://github.com/sagemath/sage/issues/31132; to use it anyway, use ./configure --with-python=$ac_path_PYTHON3])
                    ], [
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
      ])
      AS_IF([test -z "$ac_cv_path_PYTHON3"], [
          AC_MSG_NOTICE([to try to use a different system python, use ./configure --with-python=/path/to/python])
          sage_spkg_install_python3=yes
      ])
    ])
],, [
    dnl PRE
], [
    dnl POST
    AS_IF([test x$sage_spkg_install_python3 = xno], [
        PYTHON_FOR_VENV="$ac_cv_path_PYTHON3"
        AS_IF([test "$SAGE_ARCHFLAGS" = "unset"], [
           AC_MSG_CHECKING([whether $PYTHON_FOR_VENV is configured to build multiarch extensions])
           AS_IF([[CC="$CC" CXX="$CXX" conftest_venv/bin/python3 -m sysconfig | grep '^\sw*\(C\|LD\)FLAGS *=.*[" ]-arch.* -arch' ]] [>& AS_MESSAGE_LOG_FD 2>&1 ], [
               AC_MSG_RESULT([yes; disabling it by setting ARCHFLAGS])
               SAGE_ARCHFLAGS=""
           ], [
               AC_MSG_RESULT([no])
           ])
        ])
        AS_IF([test "$SAGE_ARCHFLAGS" != "unset"], [
            ARCHFLAGS="$SAGE_ARCHFLAGS"
            export ARCHFLAGS
        ])
        AS_IF([test -n "$CFLAGS_MARCH"], [
            dnl Issue #31228
            AC_MSG_CHECKING([whether "$CFLAGS_MARCH" works with the C/C++ compilers configured for building extensions for $PYTHON_FOR_VENV])
            SAGE_PYTHON_CHECK_DISTUTILS([CC="$CC" CXX="$CXX" CFLAGS="$CFLAGS_MARCH" conftest_venv/bin/python3], [
                AC_MSG_RESULT([yes])
            ], [
                AC_MSG_RESULT([no, with these flags, $reason; disabling use of "$CFLAGS_MARCH"])
                CFLAGS_MARCH=""
            ])
        ])

        AS_IF([test -n "$OPENMP_CFLAGS$OPENMP_CXXFLAGS"], [
            AC_MSG_CHECKING([whether OpenMP works with the C/C++ compilers configured for building extensions for $PYTHON_FOR_VENV])
            SAGE_PYTHON_CHECK_DISTUTILS([CC="$CC" CXX="$CXX" CFLAGS="$CFLAGS $OPENMP_CFLAGS" CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS" conftest_venv/bin/python3], [
                AC_MSG_RESULT([yes])
            ], [
                AC_MSG_RESULT([no, $reason; disabling use OpenMP])
                OPENMP_CFLAGS=""
                OPENMP_CXXFLAGS=""
            ])
        ])

        AX_COMPARE_VERSION([$python3_version], [lt], MIN_NONDEPRECATED_VERSION, [
            AC_MSG_NOTICE([deprecation notice: Support for system python < MIN_NONDEPRECATED_VERSION is deprecated
and will be removed in the next development cycle.  Consider using a newer version of Python
that may be available on your system or can be installed using the system package manager.
To build Sage with a different system python, use ./configure --with-python=/path/to/python])
        ])
        AX_COMPARE_VERSION([$python3_version], [ge], LT_STABLE_VERSION, [
            AC_MSG_WARN([Support for system python >= LT_STABLE_VERSION is experimental.
To build Sage with a different system python, use ./configure --with-python=/path/to/python])
        ])
    ])
    AC_SUBST([PYTHON_FOR_VENV])

    AS_IF([test -n "$PYTHON_FOR_VENV"],
          [PYTHON_VERSION=$("$PYTHON_FOR_VENV" -c "import sysconfig; print(sysconfig.get_python_version())")],
          [PYTHON_VERSION=$(echo $(cat build/pkgs/python3/package-version.txt))])
    AC_SUBST([PYTHON_MINOR], [$(echo $PYTHON_VERSION | cut -d. -f2)])
    export PYTHON_MINOR  # for sage-get-system-packages

    AS_VAR_IF([SAGE_VENV], [auto], [SAGE_VENV=$SAGE_VENV_AUTO])
    AS_CASE([$SAGE_VENV],
        [no],  [SAGE_VENV='${SAGE_LOCAL}'],dnl Quoted so that it is resolved at build time by shell/Makefile
        [yes], [SAGE_VENV='${SAGE_LOCAL}'/var/lib/sage/venv-python$PYTHON_VERSION]
    )
    dnl These temporary directories are created by the check above
    dnl and need to be cleaned up to prevent the "rm -f conftest*"
    dnl (that a bunch of other checks do) from emitting warnings about
    dnl conftest.dir and conftest_venv being directories.
    rm -rf conftest.dir conftest_venv

    m4_popdef([MIN_VERSION])
    m4_popdef([MIN_NONDEPRECATED_VERSION])
    m4_popdef([LT_VERSION])
])
