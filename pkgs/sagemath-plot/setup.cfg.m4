# -*- conf-unix -*-
[metadata]
name = sagemath-plot
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Plotting and graphics with Matplotlib, Three.JS, etc.
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    esyscmd(`sage-get-system-packages install-requires \
        sage_conf \
        gmpy2          \
        cysignals      \
        memory_allocator \
        sagemath_categories \
        sagemath_modules \
        numpy \
        scipy \
        pillow \
        matplotlib \
        | sed "2,\$s/^/    /;"')dnl

[options.extras_require]
test = esyscmd(`sage-get-system-packages install-requires sagemath_repl')

# extras by libraries
jsmol = esyscmd(`sage-get-system-packages install-requires jupyter_jsmol')
matplotlib =                    # no extra needed

# extras by other features
graphs = esyscmd(`sage-get-system-packages install-requires sagemath_graphs')
polyhedra = esyscmd(`sage-get-system-packages install-requires sagemath_polyhedra')
symbolics = esyscmd(`sage-get-system-packages install-requires sagemath_symbolics')
