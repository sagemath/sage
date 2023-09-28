SAGE_SPKG_CONFIGURE([scipy], [
  SAGE_SPKG_DEPCHECK([openblas], [
    SAGE_PYTHON_PACKAGE_CHECK([scipy])
  ])
])
