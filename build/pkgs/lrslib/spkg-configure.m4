SAGE_SPKG_CONFIGURE([lrslib], [
    dnl System lrslib may already be 7.x, which may be compiled with FLINT
    SAGE_SPKG_DEPCHECK([gmp mpir flint], [
        AC_CHECK_PROGS([LRSNASH], [lrsnash])
        AS_IF([test -z "$LRSNASH"], [sage_spkg_install_lrslib=yes])
    ])
])
