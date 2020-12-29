SAGE_SPKG_CONFIGURE([libxml2], [
  AC_PATH_PROG(XML2_CONFIG, xml2-config, no)
  AS_IF([test "$XML2_CONFIG" = no],
        [sage_spkg_install_libxml2=yes])
])
