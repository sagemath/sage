SAGE_SPKG_CONFIGURE([tox], [
       dnl src/tox.ini has only minimal version requirements. We use 2.5.0 just to set a baseline.
       dnl (SAGE_ROOT/tox.ini needs negated factor conditions introduced in 3.0.0, but it is
       dnl best to run it with system tox anyway.)
       m4_pushdef([TOX_MIN_VERSION], [2.5.0])
       AC_CACHE_CHECK([for tox >= ]TOX_MIN_VERSION, [ac_cv_path_TOX], [
         AC_PATH_PROGS_FEATURE_CHECK([TOX], [tox], [
            tox_version=$($ac_path_TOX --version 2> /dev/null | tail -1)
            AS_IF([test -n "$tox_version"], [
                AX_COMPARE_VERSION([$tox_version], [ge], TOX_MIN_VERSION, [
                    ac_cv_path_TOX="$ac_path_TOX"
                    ac_path_TOX_found=:
                ])
            ])
         ])
       ])
       AS_IF([test -z "$ac_cv_path_TOX"],
             [sage_spkg_install_tox=yes])
])
