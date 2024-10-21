include(`sage_spkg_versions.m4')dnl
dnl Same as setup.cfg.m4 install_requires (+ their install-requires)
dnl FIXME: should pin to built wheels.
SPKG_INSTALL_REQUIRES_gmpy2
SPKG_INSTALL_REQUIRES_cysignals
SPKG_INSTALL_REQUIRES_memory_allocator
SPKG_INSTALL_REQUIRES_ipython
SPKG_INSTALL_REQUIRES_ipywidgets
# -e ../sagemath-brial
# -e ../sagemath-categories
# -e ../sagemath-combinat
# -e ../sagemath-eclib
# -e ../sagemath-environment
# -e ../sagemath-flint
# -e ../sagemath-gap
# -e ../sagemath-glpk
# -e ../sagemath-graphs
# -e ../sagemath-groups
# -e ../sagemath-homfly
# -e ../sagemath-lcalc
# -e ../sagemath-libbraiding
# -e ../sagemath-libecm
# -e ../sagemath-linbox
# -e ../sagemath-modules
# -e ../sagemath-mpmath
# -e ../sagemath-ntl
# -e ../sagemath-objects
# -e ../sagemath-pari
# -e ../sagemath-polyhedra
# -e ../sagemath-repl
# -e ../sagemath-schemes
# -e ../sagemath-singular
