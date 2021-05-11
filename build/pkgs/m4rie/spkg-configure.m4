SAGE_SPKG_CONFIGURE([m4rie], [
    SAGE_SPKG_DEPCHECK([m4ri], [
       AC_CHECK_HEADER([m4rie/m4rie.h], [
        AC_SEARCH_LIBS([gf2e_init], [m4rie], [
        ], [sage_spkg_install_m4rie=yes])
       ], [sage_spkg_install_m4rie=yes])
    ])
])
