include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    SPKG_INSTALL_REQUIRES_setuptools
    SPKG_INSTALL_REQUIRES_wheel
]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-repl"
description = "Sage: Open Source Mathematics Software: IPython kernel, Sage preparser, doctester"
dependencies = [
    SPKG_INSTALL_REQUIRES_sagemath_objects
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_ipykernel
    SPKG_INSTALL_REQUIRES_ipython
    SPKG_INSTALL_REQUIRES_ipywidgets
    SPKG_INSTALL_REQUIRES_jupyter_client
]
dynamic = ["version"]
include(`pyproject_toml_metadata.m4')dnl'

[project.optional-dependencies]
# Improved formatting of docstrings in the help system
sphinx = [
    SPKG_INSTALL_REQUIRES_sphinx
]

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[tool.setuptools]
py-modules = [
    "sage.all__sagemath_repl",
    "sage.misc.all__sagemath_repl",
    "sage.misc.banner",
    "sage.misc.sagedoc",
    "sage.misc.sage_input",
    "sage.misc.sage_eval",
]
packages = [
    "sage.doctest",
    "sage.interfaces",
    "sage.repl",
    "sage.repl.display",
    "sage.repl.ipython_kernel",
    "sage.repl.rich_output",
]
script-files = [
    # Other scripts that should be in the path also for OS packaging of sage:
    "bin/sage-eval",
    # Included because it is useful for doctesting/coverage testing user scripts too:
    "bin/sage-runtests",
    "bin/sage-fixdoctests",
    "bin/sage-coverage",
    # Helper scripts invoked by sage script
    # (they would actually belong to something like libexec)
    "bin/sage-cachegrind",
    "bin/sage-callgrind",
    "bin/sage-massif",
    "bin/sage-omega",
    "bin/sage-valgrind",
    "bin/sage-cleaner",
    # Uncategorized scripts in alphabetical order
    "bin/sage-inline-fortran",
    "bin/sage-ipynb2rst",
    "bin/sage-ipython",
    "bin/sage-notebook",
    "bin/sage-preparse",
    "bin/sage-run",
    "bin/sage-run-cython",
    "bin/sage-startuptime.py",
]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}

[tool.setuptools.package-data]
"sage.doctest" = ["tests/*"]
"sage.repl.rich_output" = ["example*"]
