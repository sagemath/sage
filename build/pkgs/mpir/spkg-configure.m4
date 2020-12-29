SAGE_SPKG_CONFIGURE([mpir], [
dnl Implement cases for what to do on different options here
    _sage_spkg_install_gmp=no
    case "$with_mp" in
        system)
            AC_CHECK_HEADER(gmp.h, [], [_sage_spkg_install_gmp=yes])
            AC_CHECK_HEADER(gmpxx.h, [], [_sage_spkg_install_gmp=yes])
            dnl mpq_cmp_z appeared in GMP 6.1.0 and is used by pynac
            AC_SEARCH_LIBS([__gmpq_cmp_z], [gmp], [],
                [_sage_spkg_install_gmp=yes])
            SAGE_MP_LIBRARY=gmp
            ;;
        mpir)
            sage_spkg_install_mpir=yes
            SAGE_MP_LIBRARY=mpir
            ;;
        gmp)
            _sage_spkg_install_gmp=yes
            SAGE_MP_LIBRARY=gmp
            ;;
    esac
], [], [
    AC_ARG_WITH([mp],
    [AS_HELP_STRING([--with-mp=system],
                    [use the system GMP as multiprecision library, if possible (default)])
dnl Not indented because whitespace ends up in the help text
AS_HELP_STRING([--with-mp=mpir],
               [use the Sage SPKG for MPIR as multiprecision library])
AS_HELP_STRING([--with-mp=gmp],
               [use the Sage SPKG for GMP as multiprecision library])])

dnl Just parse the options here
    case "$with_mp" in
        system) ;;
        MPIR|mpir) with_mp=mpir;;
        GMP|gmp) with_mp=gmp;;
        "") with_mp=system;;
        *)
            AC_MSG_ERROR([allowed values for --with-mp are system, mpir, or gmp]);;
    esac

dnl Set SAGE_MP_LIBRARY depending on the with_mp option
    case "$with_mp" in
    mpir)
        SAGE_MP_LIBRARY=mpir
        ;;
    gmp|system)
        SAGE_MP_LIBRARY=gmp
        ;;
    esac

    AC_SUBST([SAGE_MP_LIBRARY], [$SAGE_MP_LIBRARY])
], [
    if test x$sage_spkg_install_mpir = xyes -o x$_sage_spkg_install_gmp = xyes; then
        AC_SUBST(SAGE_GMP_PREFIX, ['$SAGE_LOCAL'])
        AC_SUBST(SAGE_GMP_INCLUDE, ['$SAGE_LOCAL/include'])
        AC_MSG_RESULT([using $SAGE_MP_LIBRARY SPKG (via --with-mp=$SAGE_MP_LIBRARY)])
    else
        dnl If found, we want to get the absolute path to where we
        dnl found it for use with some packages (e.g. iml) that need
        dnl this information at configure time
        AX_ABSOLUTE_HEADER([gmp.h])
        if test x$gl_cv_absolute_gmp_h = x; then
            AC_MSG_ERROR(m4_normalize([
                failed to find absolute path to gmp.h despite it being reported
                found
            ]))
        fi
        AC_SUBST(SAGE_GMP_INCLUDE, [`AS_DIRNAME($gl_cv_absolute_gmp_h)`])
        AC_SUBST(SAGE_GMP_PREFIX, [''])
        AC_MSG_RESULT([using GMP-compatible library from the system])
    fi
])

