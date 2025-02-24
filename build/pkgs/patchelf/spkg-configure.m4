SAGE_SPKG_CONFIGURE([patchelf],[
  AC_PATH_PROG(PATCHELF, patchelf)
  AS_IF([test -z "${PATCHELF}"], [sage_spkg_install_patchelf=yes])
])
