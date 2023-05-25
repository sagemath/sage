## Pipfile with all packages in the Sage distribution and version information locked
## FIXME: Many packages still missing.
[[source]]
name = "pypi"
url = "https://pypi.org/simple"
verify_ssl = true

[dev-packages]
pkgconfig = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../pkgconfig/package-version.txt)')"
cython = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../cython/package-version.txt)')"
pycodestyle = "*"
ipykernel = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../ipykernel/package-version.txt)')"
tox = "*"
jinja2 = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../jinja2/package-version.txt)')"
pytest = "*"
ipywidgets = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../ipywidgets/package-version.txt)')"
sphinx = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../sphinx/package-version.txt)')"
rope = "*"
six = "*"
jupyter-core = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../jupyter_core/package-version.txt)')"

[packages]
numpy = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../numpy/package-version.txt)')"
cysignals = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../cysignals/package-version.txt)')"
cypari2 = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../cypari/package-version.txt)')"
gmpy2 = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../gmpy2/package-version.txt)')"
pexpect = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../pexpect/package-version.txt)')"
ipython = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../ipython/package-version.txt)')"
sympy = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../sympy/package-version.txt)')"
scipy = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../scipy/package-version.txt)')"
pplpy = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../pplpy/package-version.txt)')"
matplotlib = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../matplotlib/package-version.txt)')"
cvxopt = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../cvxopt/package-version.txt)')"
rpy2 = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../rpy2/package-version.txt)')"
networkx = "==esyscmd(`printf $(sed "s/[.]p.*//;" ../networkx/package-version.txt)')"

sagemath-standard = { path = "src" }

[requires]
python_version = "3.9"
