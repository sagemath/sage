include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-graphs
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Graphs, posets, hypergraphs, designs, abstract complexes, combinatorial polyhedra, abelian sandpiles, quivers
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_sagemath_categories

[options.extras_require]
test        = SPKG_INSTALL_REQUIRES_sagemath_repl

# libraries
bliss       = SPKG_INSTALL_REQUIRES_sagemath_bliss
gap         = SPKG_INSTALL_REQUIRES_sagemath_gap
igraph      = SPKG_INSTALL_REQUIRES_python_igraph
mcqd        = SPKG_INSTALL_REQUIRES_sagemath_mcqd
networkx    = SPKG_INSTALL_REQUIRES_networkx
tdlib       = SPKG_INSTALL_REQUIRES_sagemath_tdlib

# features
combinat    = SPKG_INSTALL_REQUIRES_sagemath_combinat
databases   = # FIXME
editor      = SPKG_INSTALL_REQUIRES_phitigra
mip         = SPKG_INSTALL_REQUIRES_sagemath_polyhedra
modules     = SPKG_INSTALL_REQUIRES_sagemath_modules
plot        = SPKG_INSTALL_REQUIRES_sagemath_plot
polyhedra   = SPKG_INSTALL_REQUIRES_sagemath_polyhedra
repl        = SPKG_INSTALL_REQUIRES_sagemath_repl
sat         = SPKG_INSTALL_REQUIRES_sagemath_combinat

standard    = sagemath-graphs[combinat,databases,mip,modules,plot,polyhedra,repl]
