SAGE_SPKG_CONFIGURE([igraph], [
  SAGE_SPKG_DEPCHECK([glpk openblas gmp], [
    dnl check for igraph with pkg-config
    PKG_CHECK_MODULES([IGRAPH], [igraph >= 0.8.3], [], [
        sage_spkg_install_igraph=yes])
  ])
])

