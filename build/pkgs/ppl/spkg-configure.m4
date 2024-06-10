SAGE_SPKG_CONFIGURE([ppl], [
  SAGE_SPKG_DEPCHECK([gcc glpk gmp], [
    # If our dependencies come from the system, then we can use the
    # system ppl, too. This macro works sort-of like the
    # PKG_CHECK_MODULES macro, defining e.g. PPL_CFLAGS when a
    # suitable version of PPL is detected. The upstream version fails
    # to differentiate between LDFLAGS and LIBS (which is in turn the
    # fault of the ppl-config program), leading to argument-order
    # problems on the command line. Our version of the macro
    # defines PPL_LIBS separately so that we can distinguish the two.
    AM_PATH_PPL([1.2], [
                  LIBS="$LIBS $PPL_LIBS"
                  sage_spkg_install_ppl=no
                ],
                [sage_spkg_install_ppl=yes])
  ])
])
