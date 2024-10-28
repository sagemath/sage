SAGE_SPKG_CONFIGURE([pip], [
  dnl always run this macro because it changes the default value of
  dnl the --with-system-<package> option.
  SAGE_PYTHON_PACKAGE_CHECK([pip])

  dnl if we might not install the spkg, make sure that "pip" is in
  dnl the user's PATH, too.
  AS_IF([test "x$sage_spkg_install_pip" != "xyes"], [
    AC_CHECK_PROG(HAVE_PIP, pip, yes, no)
    AS_IF([test "x$HAVE_PIP" = "xno"], [sage_spkg_install_pip=yes])
  ])
])
