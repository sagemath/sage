SAGE_SPKG_CONFIGURE([lrslib], [
    dnl Although lrs and lrsnash binaries are the only interfaces to lrslib that
    dnl sagelib uses, the optional polymake package uses lrslib as a library, so
    dnl the following DEPCHECK is needed.
    dnl System lrslib may already be 7.x, which may be compiled with FLINT
    SAGE_SPKG_DEPCHECK([gmp flint], [
        AC_CHECK_PROGS([LRSNASH], [lrsnash])
        AS_IF([test -z "$LRSNASH"], [
            sage_spkg_install_lrslib=yes
        ], [
            AC_MSG_CHECKING([whether $LRSNASH can handle the new input format])
            cat > conftest.lrsnash <<EOF
1 1

0

0
EOF
            AS_IF([$LRSNASH conftest.lrsnash >& AS_MESSAGE_LOG_FD 2>&1], [
                AC_MSG_RESULT([yes])
            ], [
                AC_MSG_RESULT([no])
                sage_spkg_install_lrslib=yes
            ])
            rm -f conftest.lrsnash
        ])
    ])
])
