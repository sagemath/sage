SAGE_SPKG_CONFIGURE([cvxopt], [
  SAGE_SPKG_DEPCHECK([gsl glpk suitesparse], [
    SAGE_PYTHON_PACKAGE_CHECK([cvxopt])
  ])
])
