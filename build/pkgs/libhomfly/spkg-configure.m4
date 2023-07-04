SAGE_SPKG_CONFIGURE([libhomfly], [
    SAGE_SPKG_DEPCHECK([gc], [
       AC_CHECK_HEADER([homfly.h], [
        AC_SEARCH_LIBS([homfly], [homfly], [
        ], [sage_spkg_install_libhomfly=yes])
       ], [sage_spkg_install_libhomfly=yes])
    ])
])
