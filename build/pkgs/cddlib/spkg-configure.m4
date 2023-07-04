SAGE_SPKG_CONFIGURE([cddlib], [
  SAGE_SPKG_DEPCHECK([gmp], [
    dnl The sage library uses BOTH cddexec and cddexec_gmp.
    dnl These two executables were introduced in cddlib-094j.
    AC_CHECK_PROGS([CDDEXEC], [cddexec])
    AS_IF([test x$CDDEXEC = x], [sage_spkg_install_cddlib=yes])

    AC_CHECK_PROGS([CDDEXECGMP], [cddexec_gmp])
    AS_IF([test x$CDDEXECGMP = x], [sage_spkg_install_cddlib=yes])

    dnl LattE needs redcheck_gmp...
    AC_CHECK_PROGS([REDCHECKGMP], [redcheck_gmp])
    AS_IF([test x$REDCHECKGMP = x], [sage_spkg_install_cddlib=yes])

    dnl and EITHER scdd or scdd_gmp.
    AC_CHECK_PROGS(SCDD, [scdd_gmp scdd])
    AS_IF([test x$SCDD = x], [sage_spkg_install_cddlib=yes])

    dnl https://github.com/sagemath/sage/issues/30319
    AS_IF([test -n "$CDDEXEC"], [
        AC_MSG_CHECKING([whether $CDDEXEC --redcheck works correctly for real input])
        cat > conftest.ine <<EOF
H-representation
linearity 2 1 2
begin
 5 4 real
 0.0 0.0 1.0 1.0
 -1.0 0.0 1.0 1.0
 -1.0 2.0 -1.0 0.0
 -1.0 4.0 -2.0 0.0
 1 0 0 0
end
EOF
        rm -f conftest.out
        $CDDEXEC --redcheck <conftest.ine >conftest.out 2>& AS_MESSAGE_LOG_FD
        AS_IF([grep -q "^Redundant rows.*1" conftest.out 2>& AS_MESSAGE_LOG_FD], [
            AC_MSG_RESULT([no])
            sage_spkg_install_cddlib=yes
        ], [
            AC_MSG_RESULT([yes])
        ])
    ])
    dnl Recent versions (>= 0.94k) of cddlib put cddlib's headers in
    dnl a "cddlib" subdirectory.
    AC_CHECK_HEADER([cddlib/cdd.h],[],[sage_spkg_install_cddlib=yes],[
      #include <cddlib/setoper.h>
      #include <cddlib/cddmp.h>
    ])

    dnl Both lattE and gfan try to link against libcddgmp (as
    dnl opposed to libcdd).
    AC_SEARCH_LIBS([dd_abs],[cddgmp],[],[sage_spkg_install_cddlib=yes])
  ])
])
