SAGE_SPKG_CONFIGURE([gmp], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_MPIR])
    if test x"$with_mp" = xgmp -o x"$_sage_spkg_install_gmp" = xyes; then
        if test x"$SAGE_MP_LIBRARY" = xgmp; then
            sage_spkg_install_gmp=yes
	fi
    fi
])
