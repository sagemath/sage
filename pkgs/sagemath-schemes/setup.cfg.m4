include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-schemes
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Schemes, varieties, elliptic curves, algebraic Riemann surfaces, modular forms, arithmetic dynamics
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_scipy
    SPKG_INSTALL_REQUIRES_sagemath_modules

[options.extras_require]
test = SPKG_INSTALL_REQUIRES_sagemath_repl

# extras by packages (same as sagemath-modules)
flint   = SPKG_INSTALL_REQUIRES_sagemath_flint
linbox  = SPKG_INSTALL_REQUIRES_sagemath_linbox
m4ri    = # FIXME
m4rie   = # FIXME
meataxe = SPKG_INSTALL_REQUIRES_sagemath_meataxe
mpfi    = # FIXME
ntl     = # FIXME
numpy   = SPKG_INSTALL_REQUIRES_numpy

# extras by packages (specific to sagemath-schemes)


# extras by rings; same as in sagemath-modules
RDF     = sagemath-schemes[numpy]
CDF     = sagemath-schemes[numpy]
RR      =                            # no extra needed
CC      =                            # no extra needed
RIF     =
CIF     =
RBF     = sagemath-schemes[flint]
CBF     = sagemath-schemes[flint]
GF      = sagemath-schemes[linbox]
GF2     = sagemath-schemes[m4ri]
GF2e    = sagemath-schemes[m4rie]
GF2n    = sagemath-schemes[m4rie]
GFpn    = sagemath-schemes[meataxe]
QQbar   = sagemath-schemes[NumberField]
AA      = sagemath-schemes[NumberField]
UCF     = sagemath-schemes[NumberField]
Zp      = # FIXME
Qp      = sagemath-schemes[Zp]
Zq      = sagemath-schemes[Zp]
Qq      = sagemath-schemes[Zp]
FiniteField     = sagemath-schemes[GF]
NumberField     = # FIXME
QuadraticField  = sagemath-schemes[NumberField]
CyclotomicField = sagemath-schemes[NumberField]

# extras by features
toric           = SPKG_INSTALL_REQUIRES_sagemath_polyhedra
                  SPKG_INSTALL_REQUIRES_sagemath_graphs
padics          = sagemath-schemes[Zp]

# the whole package
standard        = sagemath-schemes[toric,padics,NumberField,FiniteField,flint,linbox,mpfi,ntl,numpy]
