SAGE_SPKG_CONFIGURE([ntl], [
    m4_pushdef(SAGE_NTL_VERSION_MAJOR, [10])
    m4_pushdef(SAGE_NTL_VERSION_MINOR, [3])
    SAGE_SPKG_DEPCHECK([gmp gcc], [
        AC_CHECK_HEADER([NTL/ZZ.h], [], [sage_spkg_install_ntl=yes])
        AC_MSG_CHECKING([whether we can link a program using NTL])
        NTL_SAVED_LIBS=$LIBS
        LIBS="$LIBS -lntl"
        AC_LINK_IFELSE([
            AC_LANG_PROGRAM([[#include <NTL/ZZ.h>]],
                            [[NTL::ZZ a;]]
            )], [AC_MSG_RESULT([yes])], [
            AC_MSG_RESULT([no]); sage_spkg_install_ntl=yes
            LIBS=$NTL_SAVED_LIBS
        ])
        AC_MSG_CHECKING([NTL version >= ]SAGE_NTL_VERSION_MAJOR[.]SAGE_NTL_VERSION_MINOR)
        AC_RUN_IFELSE([
            AC_LANG_PROGRAM(
            [[#include <NTL/version.h>
              #include <stdio.h>
            ]], [[
              fprintf(stderr, "%s\n", NTL_VERSION);
              if (NTL_MAJOR_VERSION >]] SAGE_NTL_VERSION_MAJOR[[) return 0;
              else if (NTL_MAJOR_VERSION ==]] SAGE_NTL_VERSION_MAJOR[[ &&
                       NTL_MINOR_VERSION >=]] SAGE_NTL_VERSION_MINOR[[) return 0;
              else return 1;
            ]])], [
                AC_MSG_RESULT([yes])
            ], [
                AC_MSG_RESULT([no])
                sage_spkg_install_ntl=yes
            ], [
                dnl assume that the person running cross-compiling
                dnl knows what they are doing
                AC_MSG_RESULT([yes])
            ])
    ])

    m4_popdef([SAGE_NTL_VERSION_MAJOR])
    m4_popdef([SAGE_NTL_VERSION_MINOR])
], [], [], [
    if test x$sage_spkg_install_ntl = xyes; then
        AC_SUBST(SAGE_NTL_PREFIX, ['$SAGE_LOCAL'])
    else
        AC_SUBST(SAGE_NTL_PREFIX, [''])
        AX_ABSOLUTE_HEADER([NTL/ZZ.h])
        ntl_inc_ntl_dir=`AS_DIRNAME(["$gl_cv_absolute_NTL_ZZ_h"])`
        ntl_inc_dir=`AS_DIRNAME(["$ntl_inc_ntl_dir"])`
        ntl_prefix=`AS_DIRNAME(["$ntl_inc_dir"])`
        PKG_CHECK_VAR([ntl_includedir], [ntl], [includedir])
        AS_IF([test "x$ntl_includedir" = x], [AC_SUBST(NTL_INCDIR, [$ntl_inc_dir])],
              [AC_SUBST(NTL_INCDIR, [$ntl_includedir])]
             )
        PKG_CHECK_VAR([ntl_libdir], [ntl], [libdir])
        AS_IF([test "x$ntl_libdir" = x], [AC_SUBST(NTL_LIBDIR, [$ntl_prefix/lib])],
              [AC_SUBST(NTL_LIBDIR, [$ntl_libdir])]
             )
    fi
])

