dnl Usage: SAGE_SHOULD_INSTALL_GFORTRAN(reason)
dnl
dnl Use this macro to indicate that we SHOULD install GFORTRAN.
dnl In this case, GFORTRAN will be installed unless SAGE_INSTALL_GFORTRAN=no.
dnl In the latter case, a warning is given.
AC_DEFUN([SAGE_SHOULD_INSTALL_GFORTRAN], [
    if test x$SAGE_INSTALL_GFORTRAN = xexists; then
        # Already installed in Sage, but it should remain selected
        true
    elif test x$SAGE_INSTALL_GFORTRAN = xno; then
        AC_MSG_WARN([$1])
    else
        AC_MSG_NOTICE([Installing gfortran because $1])
        sage_spkg_install_gfortran=yes
    fi
])

dnl Usage: SAGE_MUST_INSTALL_GFORTRAN(reason)
dnl
dnl Use this macro to indicate that we MUST install GFORTRAN.
dnl In this case, it is an error if SAGE_INSTALL_GFORTRAN=no.
AC_DEFUN([SAGE_MUST_INSTALL_GFORTRAN], [
    if test x$SAGE_INSTALL_GFORTRAN = xexists; then
        # Already installed in Sage, but it should remain selected
        true
    elif test x$SAGE_INSTALL_GFORTRAN = xno; then
        AC_MSG_ERROR([SAGE_INSTALL_GFORTRAN is set to 'no', but $1])
    else
        AC_MSG_NOTICE([Installing gfortran because $1])
        sage_spkg_install_gfortran=yes
    fi
])

dnl This macro saves current FCFLAGS for later use.
AC_DEFUN([SAGE_SAVE_FCFLAGS], [
    sage_saved_fcflags="$FCFLAGS"
])

dnl This macro restores saved FCFLAGS.
AC_DEFUN([SAGE_RESTORE_FCFLAGS], [
    FCFLAGS="$sage_saved_fcflags"
])


SAGE_SPKG_CONFIGURE([gfortran], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_GCC])
    AC_REQUIRE([AC_PROG_FC])
    # Check that the Fortran compiler accepts free-format source code (as
    # opposed to the older fixed-format style from Fortran 77).
    # This helps verify the compiler works too, so if some idiot sets FC to
    # /usr/bin/ls, we will at least know it's not a working Fortran
    # compiler.
    AC_REQUIRE([SAGE_SAVE_FCFLAGS])
    AC_FC_FREEFORM([SAGE_HAVE_FC_FREEFORM=yes], [
        AC_MSG_NOTICE([Your Fortran compiler does not accept free-format source code])
        AC_MSG_NOTICE([which means the compiler is either seriously broken, or])
        AC_MSG_NOTICE([is too old to build Sage.])
        SAGE_HAVE_FC_FREEFORM=no])

    # AC_FC_FREEFORM may have added flags.
    # However, it is up to the individual package how they invoke the
    # Fortran compiler.
    # We only check here, whether the compiler is suitable.
    AC_REQUIRE([SAGE_RESTORE_FCFLAGS])

    AS_VAR_IF(SAGE_HAVE_FC_FREEFORM, [no], [
        AS_VAR_SET(sage_spkg_install_gfortran, [yes])
    ])

    # Special case: If we are already installing gcc then don't install
    # gfortran since it's included
    if test "x$sage_spkg_install_gcc" = "xyes" -o x$SAGE_INSTALL_GCC = xexists; then
        sage_spkg_install_gfortran=no

    else
        # AX_COMPILER_VENDOR does not work for Fortran. So we just match the name of the executable
        AS_CASE(["$FC"],
            [*gfortran*], [
                AC_MSG_CHECKING([the version of $FC])
                GFORTRAN_VERSION="`$FC -dumpversion`"
                AC_MSG_RESULT([$GFORTRAN_VERSION])
                # Add the .0 because Debian/Ubuntu gives version numbers like
                # 4.6 instead of 4.6.4 (Issue #18885)
                AS_CASE(["$GFORTRAN_VERSION.0"],
                    [[[0-3]].*|4.[[0-7]].*], [
                        # Install our own gfortran if the system-provided one is older than gcc-4.8.
                        SAGE_SHOULD_INSTALL_GFORTRAN([$FC is version $GFORTRAN_VERSION, which is quite old])
                    ],
                    [1[[6-9]].*], [
                        # Install our own gfortran if the system-provided one is newer than 15.x.
                        # See https://github.com/sagemath/sage/issues/29456, https://github.com/sagemath/sage/issues/31838
                        SAGE_MUST_INSTALL_GFORTRAN([$FC is version $GFORTRAN_VERSION, which is too recent for this version of Sage])
                    ])
            ])
    fi
])
