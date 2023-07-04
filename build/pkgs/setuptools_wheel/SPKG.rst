setuptools_wheel: Build the setuptools package as a wheel
=========================================================

After installing setuptools and wheel, we build a wheel of setuptools
to complete the set of wheels stored in our wheelhouse.

This version of setuptools is suitable for PEP 517/518/660 builds,
but it is not suitable for building ``numpy``.
