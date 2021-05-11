AC_DEFUN([SAGE_CHECK_PARI_DEP],[
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI])
    AC_MSG_CHECKING([installing pari? ])
    if test x$sage_spkg_install_pari = xyes; then
        AC_MSG_RESULT([yes; install $1 as well])
        [sage_spkg_install_]$1=yes
    else
        AC_MSG_RESULT([no])
    fi
])
