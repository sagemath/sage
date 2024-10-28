SAGE_SPKG_CONFIGURE([fpylll], [
  SAGE_SPKG_DEPCHECK([cysignals fplll numpy], [
    SAGE_PYTHON_PACKAGE_CHECK([fpylll])
  ])
])
