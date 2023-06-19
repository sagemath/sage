# -*- conf-unix -*-
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
    esyscmd(`sage-get-system-packages install-requires \
        gmpy2          \
        cysignals      \
        memory_allocator \
        sagemath_categories \
        | sed "2,\$s/^/    /;"')dnl

[options.extras_require]
test = esyscmd(`sage-get-system-packages install-requires sagemath_repl')

bliss = esyscmd(`sage-get-system-packages install-requires sagemath_bliss')
mcqd = esyscmd(`sage-get-system-packages install-requires sagemath_mcqd')
tdlib = esyscmd(`sage-get-system-packages install-requires sagemath_tdlib')
igraph = esyscmd(`sage-get-system-packages install-requires python_igraph')

networkx = esyscmd(`sage-get-system-packages install-requires networkx')
gap = esyscmd(`sage-get-system-packages install-requires sagemath_gap')

combinat =
databases =
editor = esyscmd(`sage-get-system-packages install-requires phitigra')
mip =
plot =
polyhedra = esyscmd(`sage-get-system-packages install-requires sagemath_polyhedra')
repl = esyscmd(`sage-get-system-packages install-requires sagemath_repl')
sat =

standard =
    sagemath-graphs[combinat,databases,mip,plot,polyhedra,repl]
