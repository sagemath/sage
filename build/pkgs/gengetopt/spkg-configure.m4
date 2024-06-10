SAGE_SPKG_CONFIGURE([gengetopt], [
  AC_PATH_PROG(GENGETOPT, gengetopt)
  AS_IF([test -z "${GENGETOPT}"], [sage_spkg_install_gengetopt=yes])
])
