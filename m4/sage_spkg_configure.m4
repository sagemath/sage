# SYNOPSIS
#
#   SAGE_SPKG_CONFIGURE(PACKAGE-NAME,[CHECK],[REQUIRED-CHECK])
#
# DESCRIPTION
#
#   This macro should be used in the build/<spkg>/spkg-configure.m4 templates 
#   for each SPKG (if defined) to specify how to check whether or not it is
#   required to be installed, and whether or not it's already installed.
#
#   The macro takes three arguments.  The first, PACKAGE-NAME, is simply the
#   base name of the SPKG.  The second two arguments, both optional, implement
#   two different kinds of checks (the first of which is more common):
#
#   - CHECK - this should implement a test for whether the package is already
#     available on the system and/or meets any feature tests required for
#     Sage.  If this test succeeds then the shell variable
#     sage_spkg_install_<packagename> is set to "yes".  Otherwise it is set to
#     "no".  In the case of "yes", this implies that Sage may not need to
#     install the package.
#
#   - REQUIRED-CHECK - this checks whether or not the package is a required
#     dependency of Sage at all, depending typically on the platform.  Some
#     packages (e.g. yasm, among others) are only dependencies on certain
#     platforms, and otherwise do not need to be checked for at all.  If
#     a REQUIRED-CHECK determines that the package is not required it sets
#     sage_spkg_install_<packagename>="no".
#
AC_DEFUN([SAGE_SPKG_CONFIGURE], [
AC_DEFUN_ONCE([SAGE_SPKG_CONFIGURE_]m4_toupper($1), [
m4_pushdef([SPKG_NAME], [$1])
m4_pushdef([SPKG_INSTALL_VAR], [sage_spkg_install_]SPKG_NAME)
m4_pushdef([SPKG_REQUIRE_VAR], [sage_require_]SPKG_NAME)
# BEGIN SAGE_SPKG_CONFIGURE_]m4_toupper($1)[
AC_MSG_NOTICE([=== checking whether to install the $1 SPKG ===])
m4_ifval(
[$2],
[AS_VAR_SET_IF(SPKG_INSTALL_VAR, [], SPKG_INSTALL_VAR[=no])],
[AS_VAR_SET_IF(SPKG_INSTALL_VAR, [], SPKG_INSTALL_VAR[=yes])])
m4_ifval([$3], [
AS_VAR_SET_IF(SPKG_REQUIRE_VAR, [], SPKG_REQUIRE_VAR[=no])
$3
],
[AS_VAR_SET_IF(SPKG_REQUIRE_VAR, [], SPKG_REQUIRE_VAR[=yes])])
AS_VAR_IF(SPKG_REQUIRE_VAR, [yes], [$2], SPKG_INSTALL_VAR[=no])
# END SAGE_SPKG_CONFIGURE_]m4_toupper($1)[
m4_popdef([SPKG_REQUIRE_VAR])
m4_popdef([SPKG_INSTALL_VAR])
m4_popdef([SPKG_NAME])
])
])
