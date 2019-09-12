SAGE_SPKG_CONFIGURE([gp2c], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI])
    AC_MSG_CHECKING([installing pari? ])
    if test x$sage_spkg_install_pari = xyes; then
        AC_MSG_RESULT([yes; install gp2c as well])
        sage_spkg_install_gp2c=yes
        libpari_pari_cfg='$SAGE_LOCAL/lib/pari/pari.cfg'
        AC_MSG_NOTICE([pari.cfg is $libpari_pari_cfg ])
    else
        AC_MSG_RESULT([no])
        AC_PATH_PROG([GP2C], [gp2c])
        if test x$GP2C = x; then
          AC_MSG_NOTICE([using pari/gp from the system, but building gp2c])
          AC_MSG_NOTICE([one might prefer to install a system-wide gp2c, instead])
          AC_MSG_NOTICE([and re-run configure.])
          dnl need to figure out libpari prefix
          gp_prefix=`dirname $GP`
          gp_prefix=`dirname $gp_prefix`
          AC_MSG_NOTICE([gp prefix is $gp_prefix ])
          libpari_pari_cfg=`find $gp_prefix -name pari.cfg`
          AC_MSG_NOTICE([pari.cfg is $libpari_pari_cfg ])
          sage_spkg_install_gp2c=yes
        fi
    fi
], [], [], [
    if test x$sage_spkg_install_gp2c = xyes; then
        AC_SUBST(SAGE_PARI_CFG, [$libpari_pari_cfg])
        AC_MSG_RESULT([using Sage's gp2c SPKG])
    else
        AC_SUBST(SAGE_PARI_CFG, [''])
        AC_MSG_RESULT([using gp2c from the system])
    fi
])
