SAGE_SPKG_CONFIGURE([zn_poly], [
       SAGE_SPKG_DEPCHECK([gmp], [
          AC_CHECK_HEADER([zn_poly/zn_poly.h], [
           AC_SEARCH_LIBS([zn_mod_init], [zn_poly], [
           ], [sage_spkg_install_zn_poly=yes])
          ], [sage_spkg_install_zn_poly=yes])
       ])
])
