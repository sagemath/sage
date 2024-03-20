dnl Include all from sagemath-categories (via m4 include)
include(`../sagemath_categories/src/requirements.txt.m4')
pplpy==esyscmd(`printf $(sed "s/[.]p.*//;" ../pplpy/package-version.txt)')
