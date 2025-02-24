SAGE_SPKG_CONFIGURE([cypari], [
  SAGE_SPKG_DEPCHECK([cysignals pari], [
    SAGE_PYTHON_PACKAGE_CHECK([cypari])
  ])
])
