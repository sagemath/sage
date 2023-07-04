SAGE_SPKG_CONFIGURE([rw], [
  # Check for "rw.h" in the system's include directory...
  AC_CHECK_HEADER([rw.h],
                  [sage_spkg_install_rw=no],
                  [sage_spkg_install_rw=yes])

  #...and ensure that we have at least one function "calculate_level"
  # that we need from librw. If either check fails, we want to set
  # sage_spkg_install_rw=yes. However if both checks succeed, the
  # first will set sage_spkg_install_rw=no and the second will do
  # nothing.
  AC_SEARCH_LIBS([calculate_level],
                 [rw],
                 [],
                 [sage_spkg_install_rw=yes])
])
