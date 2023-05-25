SAGE_SPKG_CONFIGURE([polymake], [
    dnl Need polymake >= 3.5 for "improve callable library compatibility with threads"
    m4_pushdef([POLYMAKE_VERSION_MIN], [3.5])
    AC_CACHE_CHECK([for polymake-config >= ]POLYMAKE_VERSION_MIN, [ac_cv_path_POLYMAKE_CONFIG], [
        AC_PATH_PROGS_FEATURE_CHECK([POLYMAKE_CONFIG], [polymake-config], [
            polymake_config_version=$($ac_path_POLYMAKE_CONFIG --version 2>&1)
            AS_IF([test -n "$polymake_config_version"], [
                AX_COMPARE_VERSION([$polymake_config_version], [ge], POLYMAKE_VERSION_MIN, [
                    ac_cv_path_POLYMAKE_CONFIG="$ac_path_POLYMAKE_CONFIG"
                ])
            ])
        ])
    ])
    AS_IF([test -z "$ac_cv_path_POLYMAKE_CONFIG"],
          [sage_spkg_install_polymake=yes],
          [AC_MSG_CHECKING([whether libpolymake works])
           AC_LANG_PUSH([C++])
           saved_CXX="$CXX"
           saved_CPPFLAGS="$CPPFLAGS"
           saved_LDFLAGS="$LDFLAGS"
           CXX="$($ac_cv_path_POLYMAKE_CONFIG --cc)"
           CPPFLAGS="$($ac_cv_path_POLYMAKE_CONFIG --includes) $($ac_cv_path_POLYMAKE_CONFIG --cflags) $CPPFLAGS"
           LDFLAGS="-lpolymake $($ac_cv_path_POLYMAKE_CONFIG --ldflags) $LDFLAGS"
           AC_RUN_IFELSE([AC_LANG_PROGRAM([[
               #include <polymake/Main.h>
             ]], [[
               polymake::Main* main_polymake_session = new polymake::Main;
               main_polymake_session->shell_enable();
               main_polymake_session->set_application("polytope");
             ]]
           )], [
             AC_MSG_RESULT([yes])
             sage_require_perl_cpan_polymake_prereq=no
             sage_require_perl_term_readline_gnu=no
           ], [
             AC_MSG_RESULT([no])
             sage_spkg_install_polymake=yes
           ])
           CXX="$saved_CXX"
           CPPFLAGS="$saved_CPPFLAGS"
           LDFLAGS="$saved_LDFLAGS"
           AC_LANG_POP([C++])
          ])
    m4_popdef([POLYMAKE_VERSION_MIN])
])
