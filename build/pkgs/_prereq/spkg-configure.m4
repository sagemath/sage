dnl We cannot check prerequisites because they are required
dnl already for running the configure script.
SAGE_SPKG_CONFIGURE([_prereq], [
  dnl Just assume that they are present.
  sage_spkg_install__prereq=no
])
