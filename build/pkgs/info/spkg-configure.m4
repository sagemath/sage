SAGE_SPKG_CONFIGURE([info], [
  AC_PATH_PROG(INFO, info)
  AS_IF([test -z "${INFO}"], [sage_spkg_install_info=yes])
])
