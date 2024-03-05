## requirements.txt for creating venvs with sagelib
##
## Usage:
##
##                   $ ../sage -sh
##         (sage-sh) $ python3 -m venv venv1
##         (sage-sh) $ source venv1/bin/activate
## (venv1) (sage-sh) $ pip install -r requirements.txt
## (venv1) (sage-sh) $ pip install -e .

dnl FIXME: Including the whole package-version.txt does not work for packages that have a patchlevel....
dnl We need a better tool to format this information.

sage-conf==esyscmd(`printf $(sed "s/[.]p.*//;" ../sage_conf/package-version.txt)')
dnl sage_setup     # Will be split out later.

dnl From build/pkgs/sagelib/dependencies
cypari2==esyscmd(`printf $(sed "s/[.]p.*//;" ../cypari/package-version.txt)')
dnl ... but building bdist_wheel of cypari2 fails with recent pip... https://github.com/sagemath/cypari2/issues/93
cysignals==esyscmd(`printf $(sed "s/[.]p.*//;" ../cysignals/package-version.txt)')
Cython==esyscmd(`printf $(sed "s/[.]p.*//;" ../cython/package-version.txt)')
gmpy2==esyscmd(`printf $(sed "s/[.]p.*//;" ../gmpy2/package-version.txt)')
importlib_metadata==esyscmd(`printf $(sed "s/[.]p.*//;" ../importlib_metadata/package-version.txt)')
importlib_resources==esyscmd(`printf $(sed "s/[.]p.*//;" ../importlib_resources/package-version.txt)')
jinja2==esyscmd(`printf $(sed "s/[.]p.*//;" ../jinja2/package-version.txt)')
dnl ... for sage_setup.autogen.interpreters
jupyter_core==esyscmd(`printf $(sed "s/[.]p.*//;" ../jupyter_core/package-version.txt)')
lrcalc==esyscmd(`printf $(sed "s/[.]p.*//;" ../lrcalc_python/package-version.txt)')
memory_allocator==esyscmd(`printf $(sed "s/[.]p.*//;" ../memory_allocator/package-version.txt)')
numpy==esyscmd(`printf $(sed "s/[.]p.*//;" ../numpy/package-version.txt)')
dnl ... already needed by sage.env
pkgconfig==esyscmd(`printf $(sed "s/[.]p.*//;" ../pkgconfig/package-version.txt)')
pplpy==esyscmd(`printf $(sed "s/[.]p.*//;" ../pplpy/package-version.txt)')
primecountpy==esyscmd(`printf $(sed "s/[.]p.*//;" ../primecountpy/package-version.txt)')
requests==esyscmd(`printf $(sed "s/[.]p.*//;" ../requests/package-version.txt)')
typing_extensions==esyscmd(`printf $(sed "s/[.]p.*//;" ../typing_extensions/package-version.txt)')

dnl From Makefile.in: SAGERUNTIME
ipython==esyscmd(`printf $(sed "s/[.]p.*//;" ../ipython/package-version.txt)')
pexpect==esyscmd(`printf $(sed "s/[.]p.*//;" ../pexpect/package-version.txt)')

dnl From Makefile.in: DOC_DEPENDENCIES
sphinx==esyscmd(`printf $(sed "s/[.]p.*//;" ../sphinx/package-version.txt)')
networkx==esyscmd(`printf $(sed "s/[.]p.*//;" ../networkx/package-version.txt)')
scipy==esyscmd(`printf $(sed "s/[.]p.*//;" ../scipy/package-version.txt)')
sympy==esyscmd(`printf $(sed "s/[.]p.*//;" ../sympy/package-version.txt)')
matplotlib==esyscmd(`printf $(sed "s/[.]p.*//;" ../matplotlib/package-version.txt)')
pillow==esyscmd(`printf $(sed "s/[.]p.*//;" ../pillow/package-version.txt)')
mpmath==esyscmd(`printf $(sed "s/[.]p.*//;" ../mpmath/package-version.txt)')
ipykernel==esyscmd(`printf $(sed "s/[.]p.*//;" ../ipykernel/package-version.txt)')
jupyter_client==esyscmd(`printf $(sed "s/[.]p.*//;" ../jupyter_client/package-version.txt)')
ipywidgets==esyscmd(`printf $(sed "s/[.]p.*//;" ../ipywidgets/package-version.txt)')

dnl Other Python packages that are standard spkg, used in doctests
cvxopt==esyscmd(`printf $(sed "s/[.]p.*//;" ../cvxopt/package-version.txt)')
rpy2==esyscmd(`printf $(sed "s/[.]p.*//;" ../rpy2/package-version.txt)')
fpylll==esyscmd(`printf $(sed "s/[.]p.*//;" ../fpylll/package-version.txt)')
dnl pycryptosat  # Sage distribution installs it as part of cryptominisat. According to its README on https://pypi.org/project/pycryptosat/: "The pycryptosat python package compiles while compiling CryptoMiniSat. It cannot be compiled on its own, it must be compiled at the same time as CryptoMiniSat."
