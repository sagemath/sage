SAGE_SPKG_CONFIGURE([cmake], [dnl
        AC_CACHE_CHECK([for cmake >= 3.22], [ac_cv_path_CMAKE], [dnl
        dnl Do not accept cmake installed via https://pypi.org/project/cmake/
        dnl in the default user scheme; it will not work in our venv because
        dnl we set PYTHONUSERBASE in sage-env.
        WITH_SAGE_PYTHONUSERBASE([dnl
            AC_PATH_PROGS_FEATURE_CHECK([CMAKE], [cmake], [dnl
                cmake_version=`$ac_path_CMAKE --version 2>&1 \
                    | $SED -n -e 's/cmake version *\([[0-9]]*\.[[0-9]]*\.[[0-9]]*\)/\1/p'`
                AS_IF([test -n "$cmake_version"], [dnl
                    AX_COMPARE_VERSION([$cmake_version], [ge], [3.22], [dnl
                        ac_cv_path_CMAKE="$ac_path_CMAKE"
                        ac_path_CMAKE_found=:
                    ])
                ])
            ])
        ])
    ])
    AS_IF([test -z "$ac_cv_path_CMAKE"], [sage_spkg_install_cmake=yes])
])
