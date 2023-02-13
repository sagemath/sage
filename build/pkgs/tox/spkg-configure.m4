SAGE_SPKG_CONFIGURE([tox], [
       dnl Use non-ancient tox with full support for PEP 517.
       m4_pushdef([TOX3_MIN_VERSION], [3.21.4])
       dnl Early 4.0.x versions have bugs regarding complex factor conditions
       m4_pushdef([TOX4_MIN_VERSION], [4.0.15])
       AC_CACHE_CHECK([for tox 3 >= ]TOX3_MIN_VERSION[ or tox 4 >= ]TOX4_MIN_VERSION, [ac_cv_path_TOX], [
         AC_PATH_PROGS_FEATURE_CHECK([TOX], [tox], [
            tox_version=$($ac_path_TOX --version 2> /dev/null | tail -1)
            AS_IF([test -n "$tox_version"], [
                AX_COMPARE_VERSION([$tox_version], [lt], [4], [
                    AX_COMPARE_VERSION([$tox_version], [ge], TOX3_MIN_VERSION, [
                        ac_cv_path_TOX="$ac_path_TOX"
                        ac_path_TOX_found=:
                    ])
                ], [
                    AX_COMPARE_VERSION([$tox_version], [ge], TOX4_MIN_VERSION, [
                        ac_cv_path_TOX="$ac_path_TOX"
                        ac_path_TOX_found=:
                    ])
                ])
            ])
         ])
       ])
       AS_IF([test -z "$ac_cv_path_TOX"],
             [sage_spkg_install_tox=yes])
       m4_popdef([TOX4_MIN_VERSION])
       m4_popdef([TOX3_MIN_VERSION])
])
