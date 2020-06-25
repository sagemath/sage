SAGE_SPKG_CONFIGURE([flintqs], [
  # The QuadraticSieve program is the only interface to FlintQS that
  # sagelib uses. As a result, we don't need to call SAGE_SPKG_DEPCHECK
  # here because there's no possibility for a library conflict.
  AC_CHECK_PROG(HAVE_QUADRATICSIEVE, QuadraticSieve, yes, no)

  # If we try to just do the obvious thing and swap the return value
  # of AC_CHECK_PROG, then ./configure outputs
  #
  #   checking for QuadraticSieve... no
  #
  # when QuadraticSieve is found... which is not great.
  #
  AS_IF([test "x$HAVE_QUADRATICSIEVE" = "xyes"],
        [sage_spkg_install_flintqs=no],
        [sage_spkg_install_flintqs=yes])
])
