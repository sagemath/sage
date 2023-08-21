SAGE_SPKG_CONFIGURE([giac], [
    SAGE_SPKG_DEPCHECK([pari], [
       dnl giac does not seem to reveal its patchlevel
       m4_pushdef([GIAC_MIN_VERSION], [1.9.0])
       m4_pushdef([GIAC_MAX_VERSION], [1.9.999])
       AC_CACHE_CHECK([for giac >= ]GIAC_MIN_VERSION[, <= ]GIAC_MAX_VERSION, [ac_cv_path_GIAC], [
         AC_PATH_PROGS_FEATURE_CHECK([GIAC], [giac], [
            giac_version=$($ac_path_GIAC --version 2> /dev/null | tail -n 1)
            AS_IF([test -n "$giac_version"], [
                AX_COMPARE_VERSION([$giac_version], [ge], GIAC_MIN_VERSION, [
                    AX_COMPARE_VERSION([$giac_version], [le], GIAC_MAX_VERSION, [
                        ac_cv_path_GIAC="$ac_path_GIAC"
                        ac_path_GIAC_found=:
                    ])
                ])
            ])
         ])
       ])
       AS_IF([test -z "$ac_cv_path_GIAC"],
             [sage_spkg_install_giac=yes])
       AC_CHECK_HEADER([giac/giac.h], [
        AC_SEARCH_LIBS([ConvertUTF16toUTF8], [giac], [
        ], [sage_spkg_install_giac=yes])
       ], [sage_spkg_install_giac=yes])
       m4_popdef([GIAC_MIN_VERSION])
       m4_popdef([GIAC_MAX_VERSION])
    ])
])
