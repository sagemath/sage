# -*- conf-unix -*-
[metadata]
name = sagemath-graphs
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Graphs, posets, hypergraphs, designs, abstract complexes, combinatorial polyhedra, abelian sandpiles, quivers
long_description = file: README.rst
long_description_content_type = text/x-rst
license = GNU General Public License (GPL) v2 or later
author = The Sage Developers
author_email = sage-support@googlegroups.com
url = https://www.sagemath.org

classifiers =
    Development Status :: 6 - Mature
    Intended Audience :: Education
    Intended Audience :: Science/Research
    License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
    Programming Language :: Python :: 3 :: Only
    Programming Language :: Python :: 3.8
    Programming Language :: Python :: 3.9
    Programming Language :: Python :: 3.10
    Programming Language :: Python :: 3.11
    Programming Language :: Python :: Implementation :: CPython
    Topic :: Scientific/Engineering :: Mathematics

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
