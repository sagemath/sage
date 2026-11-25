Numerical Tools
===============

Sage has many different components that may be useful for numerical
analysis. In particular three packages deserve mention, they are
numpy, SciPy, and cvxopt.  Numpy is an excellent package that provides
fast array facilities to python. It includes some basic linear algebra
routines, vectorized math routines, random number generators, etc. It
supports a programming style similar to one would use in matlab and
most matlab techniques have an analogue in numpy. SciPy builds on
numpy and provides many different packages for optimization, root
finding, statistics, linear algebra, interpolation, FFT and dsp tools,
etc. Finally cvxopt is an optimization package which can solve linear
and quadratic programming problems and also has a nice linear algebra
interface. Now we will spend a bit more time on each of these
packages.

Before we start let us point out
https://numpy.org/doc/stable/user/numpy-for-matlab-users.html, which has a
comparison between matlab and numpy and gives numpy equivalents of
matlab commands. If you are not familiar with matlab, that's fine, even
better, it means you won't have any pre-conceived notions of how
things should work.  Also this
https://docs.scipy.org/doc/scipy-1.8.1/scipy-ref-1.8.1.pdf
is a very nice tutorial on SciPy and numpy which is more comprehensive
than ours.

.. toctree::
   :maxdepth: 2

   numpy
   scipy
   cvxopt
