SAGE_SPKG_CONFIGURE([cddlib], [
  SAGE_SPKG_DEPCHECK([gmp mpir], [
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

    dnl Future versions (>= 0.94k) of cddlib will put these headers in
    dnl a "cddlib" subdirectory, and Debian currently relocates them
    dnl under "cdd". But for now they're at the top-level, in e.g.
    dnl /usr/include/cdd.h. The lattE and gfan packages within
    dnl SageMath both look for them there, so that's where we have to
    dnl check, passing up a chance to detect cddlib on Fedora and Debian
    dnl for now. Once all of cddlib's consumers know about the new (or
    dnl both) locations, we can update this check to support them.
    AC_CHECK_HEADER([cdd.h],[],[sage_spkg_install_cddlib=yes],[
      #include <setoper.h>
      #include <cddmp.h>
    ])

    dnl Both lattE and gfan try to link against libcddgmp (as
    dnl opposed to libcdd).
    AC_SEARCH_LIBS([dd_abs],[cddgmp],[],[sage_spkg_install_cddlib=yes])
  ])
])
