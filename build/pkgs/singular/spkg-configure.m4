SAGE_SPKG_CONFIGURE([singular], [
  SAGE_SPKG_DEPCHECK([gmp ntl flint readline mpfr cddlib], [

    AC_PATH_PROG([SINGULAR_BIN], [Singular])
    AC_SUBST([SINGULAR_BIN])
    AS_IF([test -z "${SINGULAR_BIN}"], [sage_spkg_install_singular=yes], [
      dnl Use pkg-config to ensure that Singular is new enough.
      PKG_CHECK_MODULES([SINGULAR], [Singular >= 4.2.1], [
       AC_MSG_CHECKING([that Singular's help is working])
       AS_IF([test x`printf "system(\"--browser\", \"builtin\"); \n help;" | Singular 2>&1 | grep "error\ occurred"` = x], [
        AC_MSG_RESULT(yes)
        dnl We have Singular. Now determine the shared library path on
        dnl platforms on which sage.libs.singular needs to reload the library with RTLD_GLOBAL.
        AS_CASE([$host_os],
          [cygwin*], [dnl Nothing to do
                     ],
                     [dnl Use pkg-config to get singular's libdir while we're at it. As a
                      dnl moral compromise for using pkg-config, this ultimately allows us
                      dnl to pass an absolute path to dlopen(), which is the only approach
                      dnl that POSIX guarantees will work.
                      PKG_CHECK_VAR([SINGULAR_LIB_DIR], [Singular], [libdir])
                      dnl The acl_shlibext variable is set in the top-level configure.ac.
                      AC_LANG_PUSH(C)
                      ORIG_LIBS="${LIBS}"
                      LIBS="${LIBS} -ldl"
                      AC_MSG_CHECKING([if we can dlopen($LIBSINGULAR_PATH)])
                      LIBSINGULAR_PATH="${SINGULAR_LIB_DIR}/libSingular.${acl_shlibext}"

                      dnl if we can dlopen() it, substitute the name for sage_conf;
                      dnl otherwise, fall back to using the SPKG.
                      AC_RUN_IFELSE(
                        [AC_LANG_PROGRAM(
                          [[#include <dlfcn.h>]],
                          [[void* h = dlopen("${LIBSINGULAR_PATH}", RTLD_LAZY | RTLD_GLOBAL);
                            if (h == 0) { return 1; } else { return dlclose(h); }]]
                        )], [
                          AC_MSG_RESULT(yes)
                        ], [
                          dnl try Debian-specific name
                          LIBSINGULAR_PATH="${SINGULAR_LIB_DIR}/libsingular-Singular.${acl_shlibext}"
                          AC_RUN_IFELSE(
                           [AC_LANG_PROGRAM(
                             [[#include <dlfcn.h>]],
                             [[void* h = dlopen("${LIBSINGULAR_PATH}", RTLD_LAZY | RTLD_GLOBAL);
                               if (h == 0) { return 1; } else { return dlclose(h); }]]
                           )], [
                             AC_MSG_RESULT(yes)
                           ], [
                            AC_MSG_RESULT(no)
                            sage_spkg_install_singular=yes
                          ], [AC_MSG_RESULT(yes)])
                        ], [AC_MSG_RESULT(yes)])

                      AC_LANG_POP()
                      LIBS="${ORIG_LIBS}"
                     ]
       )], [
       AC_MSG_RESULT(no)
       sage_spkg_install_singular=yes
       ]
      )], [
      dnl pkg-config version check failed
      sage_spkg_install_singular=yes
      ])
    ])
  ])
],[],[],[
  dnl Post-check phase
  dnl We make the sage_conf substitutions here, because the "default"
  dnl substitution needs to be made even if we skipped the system-Singular
  dnl checks themselves.
  AS_IF([test "x${sage_spkg_install_singular}" = "xyes"], [
    AS_CASE([$host_os],
      [cygwin*], [dnl Nothing to do
                 ],
                 [dnl Set shared library path, needed for reloading the library with RTLD_GLOBAL
                  LIBSINGULAR_PATH="\$SAGE_LOCAL/lib/libSingular.${acl_shlibext}"
                 ]
    )
  ])

  AC_SUBST(LIBSINGULAR_PATH, "${LIBSINGULAR_PATH}")
])
