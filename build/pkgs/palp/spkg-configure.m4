SAGE_SPKG_CONFIGURE([palp], [
   m4_foreach([palpprog], [[poly], [class], [nef], [cws]], [
       AC_PATH_PROG(PALP[]palpprog, [palpprog.x])
       AS_IF([test "x$PALP[]palpprog" = "x"], [sage_spkg_install_palp=yes])
       m4_foreach([suff], [4, 5, 6, 11], [
           AC_PATH_PROG(PALP[]palpprog[]suff, [palpprog[-]suff[d.x]])
           AS_IF([test "x$PALP[]palpprog[]suff" = "x"], [sage_spkg_install_palp=yes])
       ])
   ])
])
