SAGE_SPKG_CONFIGURE([pplpy], [
  SAGE_SPKG_DEPCHECK([cysignals gmpy2 ppl], [
    SAGE_PYTHON_PACKAGE_CHECK([pplpy])
  ])
])
