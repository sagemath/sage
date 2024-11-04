include(`sage_spkg_versions.m4')dnl
dnl Same as setup.cfg.m4 install_requires; FIXME: should pin to built wheels.
SPKG_INSTALL_REQUIRES_gmpy2
SPKG_INSTALL_REQUIRES_cysignals
-e ../sagemath-objects
-e ../sagemath-categories
