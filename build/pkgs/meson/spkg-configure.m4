SAGE_SPKG_CONFIGURE([meson], [dnl
    dnl scipy 1.15.2 needs meson >= 1.5.0
    AC_CACHE_CHECK([for meson >= 1.5.0], [ac_cv_path_MESON], [dnl
        dnl Do not accept meson installed in the default user scheme;
        dnl it will not work in our venv because we set PYTHONUSERBASE
        dnl in sage-env.
        WITH_SAGE_PYTHONUSERBASE([dnl
            AC_PATH_PROGS_FEATURE_CHECK([MESON], [meson], [dnl
                AS_IF([meson_version=$($ac_path_MESON --version 2>&1)], [dnl
                    AS_IF([test -n "$meson_version"], [dnl
                        AX_COMPARE_VERSION([$meson_version], [ge], [1.5.0], [dnl
                            ac_cv_path_MESON="$ac_path_MESON"
                            ac_path_MESON_found=:
                        ])
                    ])
                ])
            ])
        ])
    ])
    AS_IF([test -z "$ac_cv_path_MESON"], [sage_spkg_install_meson=yes])
])
