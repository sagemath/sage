SAGE_SPKG_CONFIGURE([linbox], [
  SAGE_SPKG_DEPCHECK([fflas_ffpack flint fplll givaro gmp iml m4ri m4rie mpfr ntl], [
    PKG_CHECK_MODULES([LINBOX],dnl Check for a set of matching old versions
      [linbox >= 1.7.0 linbox < 1.8.0 fflas-ffpack >= 2.5.0 fflas-ffpack < 2.6.0 givaro >= 4.2.0 givaro < 4.3.0],
      [sage_spkg_install_linbox=no],
      [PKG_CHECK_MODULES([LINBOX],dnl Check for a set of matching new versions
        [linbox >= 1.8.0 linbox <= 1.9.0 fflas-ffpack >= 2.6.0 givaro >= 4.3.0 givaro < 4.4.0],
        [sage_spkg_install_linbox=no],
        [sage_spkg_install_linbox=yes])])
  ])
], [dnl REQUIRED_CHECK
], [dnl PRE
], [dnl POST
   sage_spkg_install_fflas_ffpack=$sage_spkg_install_linbox
   sage_spkg_install_givaro=$sage_spkg_install_linbox
])
