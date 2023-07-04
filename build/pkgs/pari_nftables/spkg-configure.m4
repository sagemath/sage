SAGE_SPKG_CONFIGURE([pari_nftables], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI])
    AC_MSG_CHECKING([installing pari? ])
    if test x$sage_spkg_install_pari = xyes; then
        AC_MSG_RESULT([yes; install pari_nftables as well])
        sage_spkg_install_pari_nftables=yes
    else
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([Installing nftables should be done by the user,])
        AC_MSG_NOTICE([cf. http://pari.math.u-bordeaux1.fr/pub/pari/packages/nftables/README.txt])
    fi
])
