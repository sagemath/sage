SAGE_SPKG_CONFIGURE([numpy], [
  SAGE_SPKG_DEPCHECK([openblas], [
   SAGE_PYTHON_PACKAGE_CHECK([numpy])
  ])
])
