SAGE_SPKG_CONFIGURE([fplll], [
  SAGE_SPKG_DEPCHECK([gcc mpfr], [
    dnl If we're using the system mpfr, use pkgconfig to determine
    dnl if there's a usable system copy of fplll. Unless there's
    dnl a system that ships fplll without fplll.pc file, falling
    dnl back to a manual header/library search is pointless.
    dnl
    dnl Trac #31025: FPLLL/FPyLLL make no guarantee regarding compatibility
    dnl other than "whatever versions were released at the same time should work together"
    PKG_CHECK_MODULES([FPLLL],
                      [fplll = 5.4.0],
                      [sage_spkg_install_fplll=no],
                      [sage_spkg_install_fplll=yes])
  ])
])
