SAGE_SPKG_CONFIGURE([highs], [
    m4_pushdef([SAGE_HIGHS_MINVER],["1.12.0"])
    dnl First check if HiGHS is available via pkg-config
    PKG_CHECK_MODULES([HIGHS], [highs >= $SAGE_HIGHS_MINVER], [
        dnl pkg-config found HiGHS, now verify the headers are actually usable
        dnl Some distro packages (e.g., Debian) have broken pkg-config files
        dnl with incorrect include paths from the build directory
        save_CPPFLAGS="$CPPFLAGS"
        CPPFLAGS="$CPPFLAGS $HIGHS_CFLAGS"
        AC_MSG_CHECKING([whether HiGHS headers are usable])
        AC_COMPILE_IFELSE([
            AC_LANG_PROGRAM([[
                #include <interfaces/highs_c_api.h>
                typedef char sage_highsint_must_match_int[
                    (sizeof(HighsInt) == sizeof(int)) ? 1 : -1
                ];
            ]], [[
                Highs_create();
            ]])
        ], [
            AC_MSG_RESULT([yes])
        ], [
            AC_MSG_RESULT([no])
            AC_MSG_NOTICE([HiGHS pkg-config found but headers not usable.])
            AC_MSG_NOTICE([This may indicate a broken system package with incorrect include paths.])
            AC_MSG_NOTICE([It may also indicate that HiGHS was built with 64-bit HighsInt,])
            AC_MSG_NOTICE([which is not compatible with Sage's current HiGHS Cython wrapper.])
            AC_MSG_NOTICE([Will build HiGHS from source.])
            sage_spkg_install_highs=yes
        ])
        CPPFLAGS="$save_CPPFLAGS"
    ], [
        dnl pkg-config did not find HiGHS
        sage_spkg_install_highs=yes
    ])
    m4_popdef([SAGE_HIGHS_MINVER])
])
