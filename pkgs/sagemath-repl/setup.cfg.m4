# -*- conf-unix -*-
[metadata]
name = sagemath-repl
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: System and software environment
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    esyscmd(`sage-get-system-packages install-requires \
        sagemath_objects \
        sagemath_environment \
        ipython \
        ipywidgets \
        | sed "2,\$s/^/    /;"')dnl

py_modules =
    sage.all__sagemath_repl
    sage.misc.all__sagemath_repl
    sage.misc.banner
    sage.misc.sagedoc
    sage.misc.sage_input
    sage.misc.sage_eval

packages =
    sage.doctest
    sage.repl
    sage.repl.display
    sage.repl.ipython_kernel
    sage.repl.rich_output

scripts =
    # Other scripts that should be in the path also for OS packaging of sage:
    bin/sage-eval
    # Included because it is useful for doctesting/coverage testing user scripts too:
    bin/sage-runtests
    bin/sage-fixdoctests
    bin/sage-coverage
    # Helper scripts invoked by sage script
    # (they would actually belong to something like libexec)
    bin/sage-cachegrind
    bin/sage-callgrind
    bin/sage-massif
    bin/sage-omega
    bin/sage-valgrind
    bin/sage-cleaner
    # Uncategorized scripts in alphabetical order
    bin/sage-inline-fortran
    bin/sage-ipynb2rst
    bin/sage-ipython
    bin/sage-notebook
    bin/sage-preparse
    bin/sage-run
    bin/sage-run-cython
    bin/sage-startuptime.py

[options.package_data]

sage.doctest =
    tests/*

sage.repl.rich_output =
    example*
