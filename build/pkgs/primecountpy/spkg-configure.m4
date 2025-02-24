SAGE_SPKG_CONFIGURE([primecountpy], [
  SAGE_SPKG_DEPCHECK([cysignals primecount], [
    SAGE_PYTHON_PACKAGE_CHECK([primecountpy])
  ])
])
