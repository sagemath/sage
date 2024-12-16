SAGE_SPKG_CONFIGURE([ninja_build], [dnl
    dnl meson_python needs 1.8.2 or later
    AC_CACHE_CHECK([for ninja >= 1.8.2], [ac_cv_path_NINJA], [dnl
        dnl Do not accept ninja installed from https://pypi.org/project/ninja/
        dnl in the default user scheme; it will not work in our venv because
        dnl we set PYTHONUSERBASE in sage-env.
        WITH_SAGE_PYTHONUSERBASE([dnl
            AC_PATH_PROGS_FEATURE_CHECK([NINJA], [ninja], [dnl
                dnl support both two- and three-component version schemes
                dnl since samurai (a ninja alternative) uses two
                ninja_version=`$ac_path_NINJA --version 2>&1 \
                    | $SED -n -e 's/\([[0-9]]*\(\.[[0-9]]*\)\{1,2\}\).*/\1/p'`
                AS_IF([test -n "$ninja_version"], [dnl
                    AX_COMPARE_VERSION([$ninja_version], [ge], [1.8.2], [
                        ac_cv_path_NINJA="$ac_path_NINJA"
                        ac_path_NINJA_found=:
                    ])
                ])
            ])
        ])
    ])
    AS_IF([test -z "$ac_cv_path_NINJA"], [sage_spkg_install_ninja_build=yes])
])
