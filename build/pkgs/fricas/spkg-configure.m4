SAGE_SPKG_CONFIGURE(
    [fricas], [
        AC_CACHE_CHECK([for FriCAS >= 1.3.8], [ac_cv_path_FRICAS], [
        AC_PATH_PROGS_FEATURE_CHECK([FRICAS], [fricas], [
            fricas_version=`$ac_path_FRICAS -nox -noclef -eval ")quit" | grep Version | tail -1 2>&1 \
                | $SED -n -e 's/.* Version: FriCAS //p'`
            echo "$fricas_version"
            AS_IF([test -n "$fricas_version"], [
                AX_COMPARE_VERSION([$fricas_version], [ge], [1.3.8], [
                    ac_cv_path_FRICAS="$ac_path_FRICAS"
                    ac_path_FRICAS_found=:
                ])
            ])
        ])
    ])
    AS_IF([test -z "$ac_cv_path_FRICAS"], [sage_spkg_install_fricas=yes])
])
