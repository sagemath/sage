SAGE_SPKG_CONFIGURE([igraph], [
  SAGE_SPKG_DEPCHECK([glpk openblas gmp], [
    dnl check for igraph with pkg-config
    PKG_CHECK_MODULES([IGRAPH], [igraph >= 0.10 igraph < 0.11], [], [
        sage_spkg_install_igraph=yes])
  ])
])

