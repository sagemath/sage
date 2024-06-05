SAGE_SPKG_CONFIGURE([igraph], [
  SAGE_SPKG_DEPCHECK([glpk openblas gmp], [
    dnl check for igraph with pkg-config
    dnl Per upstream in https://github.com/sagemath/sage/pull/36750#issuecomment-1826998762:
    dnl each python-igraph release is only guaranteed to be compatible with the same C/igraph that it bundles
    PKG_CHECK_MODULES([IGRAPH], [igraph >= 0.10.12 igraph < 0.10.13], [], [
        sage_spkg_install_igraph=yes])
  ])
])

